/* lhsight2 (Steam Frame): lighthouse-blob detector, rebuild v1.
 *
 * Replaces lhsight's "saturated blobs on short frames" with a temporal
 * classifier tuned for the Frame's SLAM cameras:
 *
 *   - LONG frames are scanned (short frames are the empty ones in a dark
 *     room; stock lhsight scanned only shorts and starved). Short frames are
 *     still reported as "F cam seq t_us mean -1" so the PC's frame-grid
 *     learning keeps its regular cadence.
 *   - Blobs are tracked across frames per camera (nearest within RADIUS px).
 *   - Each track's amplitude (blob peak, 0 when absent) is sampled once per
 *     long frame (~33.3 Hz). A track is CLASSIFIED when its amplitude spectrum
 *     over the last DFT_WIN long frames has a single stable dominant peak:
 *       presence >= PRES_MIN, size <= MAXNPX.
 *     - Rotor phase drifts against the frame grid, so a real base station's
 *       dot AM-modulates (measured: ~9.8 Hz, 73% presence). The rig's own IR
 *       illuminators are strobe-locked to long frames: constant amplitude,
 *       no peak -> rejected. Mains lamps are common-mode; the per-track mean
 *       is subtracted and huge blobs never enter tracks.
 *   - Only CLASSIFIED tracks' blobs go out on the wire, in the stock
 *     4-tuple format ("F cam seq t_us mean nb x y npx peak ..."), so the PC
 *     driver works unchanged. Diagnostics go on "B" lines (ignored by the PC).
 *
 * Lines:
 *   F cam seq t_us mean*10 nb [x y npx peak]...   nb=-1 on short frames
 *   B cam id x y f_hz q pres hits state           once a second per track
 *   T t_us frames scans scan_us torn polls        once a second
 *   I/W/E messages
 *
 * Build: cc -O2 -o lhsight2 lhsight2.c
 * The borrow protocol is identical to lhsight.c (questlhsync_frame driver,
 * SCM_RIGHTS over $XDG_RUNTIME_DIR/questlhsync/xrfds.sock).
 */
#define _GNU_SOURCE
#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <linux/dma-buf.h>
#include <linux/videodev2.h>
#include <math.h>
#include <poll.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ioctl.h>
#include <sys/mman.h>
#include <sys/socket.h>
#include <sys/syscall.h>
#include <sys/un.h>
#include <time.h>
#include <unistd.h>

#define NCAM 4
#define NBUF 16
#define MAXB 48
#define NPITCH 8

/* detector tuning */
#define RADIUS      24.0   /* px: blob <-> track association radius */
#define RADIUS2     (RADIUS*RADIUS)
#define MAXTRACKS   8      /* per camera (4 stations max + margin) */
#define DFT_WIN     768    /* amplitude ring: long frames (~23 s at 30 Hz) */
#define FAST_WIN    384    /* fast-beat window (~11.5 s): 60 Hz-mode beats */
#define DFT_MIN     96     /* minimum samples before classifying (~3.2 s) */
#define LONG_MIN    640    /* minimum samples for the slow-beat window (~21 s) */
#define DFT_STEP    64     /* re-analyse every N long frames */
#define NFREQ       48     /* fast window: 0.5 .. 24.0 Hz in 0.5 Hz steps */
#define NFREQ_SLOW  64     /* slow window: 0.1 .. 6.4 Hz in 0.1 Hz steps */
#define FMIN        0.5f
#define FSTEP       0.5f
#define FMIN_SLOW   0.1f
#define FSTEP_SLOW  0.1f
#define Q_MIN       5.0f   /* peak / median floor of the spectrum */
#define F_STAB      0.75f  /* Hz: first-half vs second-half peak agreement */
#define F_STAB_SLOW 0.25f  /* slow beats need finer stability */
#define PRES_MIN    0.12f  /* fraction of long frames with the dot present */
#define MAXNPX      150    /* bigger blobs never enter tracks */
#define PK_MIN      180    /* blob peak needed to count as "present" */
#define MIN_HITS    12     /* real dot hits before a track may classify */
#define TRACK_TTL   2.5    /* s without any hit: drop the track */
#define EMIT_MIN    0.0    /* (classification gates emission anyway) */

typedef uint8_t u8; typedef uint16_t u16; typedef uint32_t u32; typedef uint64_t u64;

enum { T_NEW = 0, T_ACTIVE = 1 };

struct track {
  int used, id, state;
  float x, y;              /* px, EMA of hit positions */
  double last_hit;         /* monotonic s of last real hit */
  int hits;                /* real hits total */
  int n;                   /* samples in ring (long frames) */
  u8 amp[DFT_WIN];         /* amplitude ring (peak, 0 = absent) */
  u8 npx[DFT_WIN];
  float f0; float q; float pres;   /* last analysis */
  int tick;                        /* analysis pacing once the ring is full */
  double since_classified;
};

struct cam {
  char dev[64];
  int vfd, w, h, nb;
  int stride, bpl;
  int pitch[NPITCH], npitch, vote, votes, tries;
  u32 k;
  int xfd[NBUF], fd[NBUF];
  u8 *map[NBUF]; size_t len[NBUF];
  unsigned seq[NBUF];
  int m1, m2;
  int rmax;                /* decaying running max of mean (exposure class) */
  double last, prev, prev2, polled;
  struct track tr[MAXTRACKS];
  int next_id;
  u32 longs;               /* long frames processed (analysis pacing) */
  double lfts[32];         /* recent long-frame times (s): the real sample rate */
  int nlfts;
  int nlfts_reported;
  float last_fs;           /* the rate check_rate() last saw (mode-switch detect) */
  float cand_fs; int cand_n;
};
static struct cam C[NCAM];
static u8 *scratch;
static u16 *cells; static int *lab; static u8 *hot;
static u32 nscan, ntorn, npoll;
static u64 scan_ns;
static clockid_t ts_clock = -1;

static void out(const char *fmt, ...) {
  char b[2600]; va_list ap; va_start(ap, fmt);
  int n = vsnprintf(b, sizeof b - 1, fmt, ap); va_end(ap);
  if (n < 0) return; if (n > (int)sizeof b - 2) n = sizeof b - 2; b[n++] = '\n';
  for (int o = 0; o < n;) { ssize_t w = write(1, b+o, n-o); if (w <= 0) { if (errno == EINTR) continue; _exit(1); } o += w; }
}
static double now_s(clockid_t c) { struct timespec t; clock_gettime(c, &t); return t.tv_sec + t.tv_nsec/1e9; }
static u64 now_ns(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return (u64)t.tv_sec*1000000000ull + t.tv_nsec; }

/* ---------------------------------------------------------------- discovery */
static pid_t xrservice_pid(void) {
  DIR *d = opendir("/proc"); struct dirent *e; pid_t found = 0;
  if (!d) return 0;
  while (!found && (e = readdir(d))) {
    if (e->d_name[0] < '0' || e->d_name[0] > '9') continue;
    char p[300], cmd[512]; snprintf(p, sizeof p, "/proc/%s/cmdline", e->d_name);
    int fd = open(p, O_RDONLY); if (fd < 0) continue;
    ssize_t n = read(fd, cmd, sizeof cmd - 1); close(fd);
    if (n <= 0) continue; cmd[n] = 0;
    const char *base = strrchr(cmd, '/') ? strrchr(cmd, '/') + 1 : cmd;
    if (!strcmp(base, "XRService")) found = atoi(e->d_name);
  }
  closedir(d); return found;
}
static int camera_nodes(pid_t xr) {
  char dir[64], log[512] = "";
  snprintf(dir, sizeof dir, "/proc/%d/fd", xr);
  DIR *d = opendir(dir); if (!d) return 0;
  struct dirent *e;
  while ((e = readdir(d))) {
    char p[400], l[512];
    snprintf(p, sizeof p, "%s/%s", dir, e->d_name);
    ssize_t n = readlink(p, l, sizeof l - 1); if (n <= 0) continue; l[n] = 0;
    if (strstr(l, "/XRService-") && n > 4 && !strcmp(l+n-4, ".log")) { snprintf(log, sizeof log, "%s", l); break; }
  }
  closedir(d);
  FILE *f = log[0] ? fopen(log, "r") : NULL; if (!f) return 0;
  char line[1024]; int found = 0;
  while (fgets(line, sizeof line, f)) {
    const char *s = strstr(line, "TrackingCameraInit: index: ");
    int i; char dev[64];
    if (!s || sscanf(s, "TrackingCameraInit: index: %d. video device: %63[^. \n]", &i, dev) != 2 || i < 0 || i >= NCAM) continue;
    if (!C[i].dev[0]) found++;
    snprintf(C[i].dev, sizeof C[i].dev, "%s", dev);
  }
  fclose(f); return found;
}

/* ---------------------------------------------------------------- v4l2 */
static int querybuf(struct cam *c, int i, struct v4l2_buffer *b, struct v4l2_plane *pl) {
  memset(b, 0, sizeof *b); memset(pl, 0, VIDEO_MAX_PLANES * sizeof *pl);
  b->type = V4L2_BUF_TYPE_VIDEO_CAPTURE_MPLANE; b->index = i; b->m.planes = pl; b->length = VIDEO_MAX_PLANES;
  return ioctl(c->vfd, VIDIOC_QUERYBUF, b);
}
static void add_pitch(struct cam *c, size_t s, size_t len) {
  if (s < (size_t)c->w || s * (size_t)(c->h-1) + (size_t)c->w > len || c->npitch == NPITCH) return;
  for (int i = 0; i < c->npitch; i++) if (c->pitch[i] == (int)s) return;
  c->pitch[c->npitch++] = (int)s;
}
static double row_step(struct cam *c, const u8 *p, int s, int rows) {
  u32 sum = 0, n = 0;
  for (int y = 0; y+1 < rows; y += 4)
    for (int x = 0; x < c->w; x += 4) {
      int dd = p[(size_t)y*s+x] - p[(size_t)(y+1)*s+x]; sum += (u32)(dd < 0 ? -dd : dd); n++;
    }
  return n ? (double)sum/n : 0;
}
static int pick_pitch(struct cam *c, const u8 *p) {
  int lo = c->pitch[0], hi = c->pitch[0];
  for (int i = 1; i < c->npitch; i++) { if (c->pitch[i] < lo) lo = c->pitch[i]; if (c->pitch[i] > hi) hi = c->pitch[i]; }
  int rows = (int)((long)(c->h-1)*lo/hi), best = -1, first = 1; double d[NPITCH], second = 0;
  for (int i = 0; i < c->npitch; i++) { d[i] = row_step(c, p, c->pitch[i], rows); if (best < 0 || d[i] < d[best]) best = i; }
  for (int i = 0; i < c->npitch; i++) if (i != best && (first || d[i] < second)) { second = d[i]; first = 0; }
  return (second >= 1.3*d[best] && second - d[best] >= 2) ? best : -1;
}
static void count_pitch(struct cam *c, int v) {
  if (++c->tries == 300) out("W %s: row pitch still unknown (dark?)", c->dev);
  if (v < 0) return;
  if (v != c->vote) { c->vote = v; c->votes = 0; }
  if (++c->votes < 3) return;
  c->stride = c->pitch[v];
  out("I %s: stride %d (v4l2 %d)", c->dev, c->stride, c->bpl);
}
static int open_camera(struct cam *c) {
  c->vfd = open(c->dev, O_RDWR | O_NONBLOCK | O_CLOEXEC);
  if (c->vfd < 0) { out("E open %s: %s", c->dev, strerror(errno)); return -1; }
  struct v4l2_format f; memset(&f, 0, sizeof f); f.type = V4L2_BUF_TYPE_VIDEO_CAPTURE_MPLANE;
  if (ioctl(c->vfd, VIDIOC_G_FMT, &f)) { out("E %s G_FMT: %s", c->dev, strerror(errno)); return -1; }
  c->w = f.fmt.pix_mp.width; c->h = f.fmt.pix_mp.height; c->bpl = f.fmt.pix_mp.plane_fmt[0].bytesperline;
  size_t len = 0;
  for (c->nb = 0; c->nb < NBUF; c->nb++) {
    struct v4l2_buffer b; struct v4l2_plane pl[VIDEO_MAX_PLANES];
    if (querybuf(c, c->nb, &b, pl)) break;
    if (b.memory != V4L2_MEMORY_DMABUF) { out("E %s: buffers not dmabuf", c->dev); return -1; }
    c->xfd[c->nb] = pl[0].m.fd; c->len[c->nb] = pl[0].length; c->seq[c->nb] = b.sequence;
    if (!len || pl[0].length < len) len = pl[0].length;
  }
  if (!c->nb) { out("E %s: no buffers (XRService not streaming?)", c->dev); return -1; }
  c->npitch = 0;
  if (c->bpl > 0) add_pitch(c, (size_t)c->bpl, len);
  if (len % ((size_t)c->h*3/2) == 0) add_pitch(c, len/((size_t)c->h*3/2), len);
  if (len % (size_t)c->h == 0) add_pitch(c, len/(size_t)c->h, len);
  for (size_t a = 32; a <= 256; a *= 2) add_pitch(c, ((size_t)c->w + a - 1)/a*a, len);
  c->stride = c->npitch == 1 ? c->pitch[0] : 0;
  c->vote = -1;
  return 0;
}
static int borrow(pid_t xr) {
  int s = socket(AF_UNIX, SOCK_STREAM | SOCK_CLOEXEC, 0);
  struct sockaddr_un a = {.sun_family = AF_UNIX};
  const char *rt = getenv("XDG_RUNTIME_DIR");
  snprintf(a.sun_path, sizeof a.sun_path, "%s/questlhsync/xrfds.sock", rt ? rt : "/run/user/1000");
  if (connect(s, (struct sockaddr *)&a, sizeof a)) { out("E no xrfds.sock: %s (questlhsync_frame driver running?)", strerror(errno)); return -1; }
  char req[1024]; int n = snprintf(req, sizeof req, "G"), want = 0;
  for (int k = 0; k < NCAM; k++) for (int i = 0; i < C[k].nb; i++) n += snprintf(req+n, sizeof req-n, " %d", C[k].xfd[i]), want++;
  req[n++] = '\n';
  if (send(s, req, n, MSG_NOSIGNAL) != n) { out("E driver send"); return -1; }
  char txt[256] = ""; char ctl[CMSG_SPACE(sizeof(int)*NCAM*NBUF)];
  struct iovec io = {txt, sizeof txt-1};
  struct msghdr m = {.msg_iov = &io, .msg_iovlen = 1, .msg_control = ctl, .msg_controllen = sizeof ctl};
  struct pollfd pf = {s, POLLIN, 0};
  ssize_t r = poll(&pf, 1, 3000) == 1 ? recvmsg(s, &m, MSG_CMSG_CLOEXEC) : -1;
  close(s);
  if (r <= 0) { out("E driver answer"); return -1; }
  txt[r] = 0;
  int pid = 0, got = 0;
  if (sscanf(txt, "OK %d %d", &pid, &got) != 2) { out("E driver: %s", txt); return -1; }
  struct cmsghdr *h = CMSG_FIRSTHDR(&m);
  int nfd = h && h->cmsg_type == SCM_RIGHTS ? (int)((h->cmsg_len - CMSG_LEN(0))/sizeof(int)) : 0;
  int fds[NCAM*NBUF];
  if (nfd > 0) memcpy(fds, CMSG_DATA(h), sizeof(int)*(nfd < NCAM*NBUF ? nfd : NCAM*NBUF));
  if (pid != xr || got != want || nfd != want) { out("E driver lent %d of %d", nfd, want); return -1; }
  int j = 0;
  for (int k = 0; k < NCAM; k++) for (int i = 0; i < C[k].nb; i++, j++) {
    C[k].fd[i] = fds[j];
    C[k].map[i] = mmap(NULL, C[k].len[i], PROT_READ, MAP_SHARED, fds[j], 0);
    if (C[k].map[i] == MAP_FAILED) { out("E mmap"); return -1; }
  }
  return 0;
}
static void cpu_access(int fd, int end) {
  struct dma_buf_sync s = {.flags = DMA_BUF_SYNC_READ | (end ? DMA_BUF_SYNC_END : DMA_BUF_SYNC_START)};
  ioctl(fd, DMA_BUF_IOCTL_SYNC, &s);
}
static int smean(const struct cam *c, const u8 *p) {
  u32 n = (u32)(c->w*c->h), step = n/4096, s = 0, k = 0;
  for (u32 o = step/2; o < n; o += step) { s += p[(o/c->w)*c->stride + o%c->w]; k++; }
  return (int)(s*10/k);
}

/* ---------------------------------------------------------------- blobs */
static int find(int x) { while (lab[x] != x) { lab[x] = lab[lab[x]]; x = lab[x]; } return x; }
static int blobs(int w, int h, int T, int mark, int nbl, u64 *sx, u64 *sy, u32 *np, int *pk, u32 *sat) {
  int cw = w>>4, ch = h>>4, nc = cw*ch;
  for (int q = 0; q < nc; q++) cells[q] = 0;
  const u64 add = 0x0101010101010101ull * (u64)(127 - (T-1 > 127 ? 127 : T-1));
  const u64 *q8 = (const u64 *)scratch;
  u32 n = (u32)(w*h), nw = n>>3;
  for (u32 j = 0; j < nw; j++) {
    u64 x = q8[j];
    if (!(((x+add)|x) & 0x8080808080808080ull)) continue;
    for (int bb = 0; bb < 8; bb++) {
      if ((int)((x>>(8*bb)) & 0xFF) < T) continue;
      u32 px = (j<<3)+bb; int y = (int)(px/(u32)w), xx = (int)(px - (u32)y*(u32)w);
      int q = (y>>4)*cw + (xx>>4);
      if (hot[q]) continue;
      if (cells[q] < 65535) cells[q]++;
    }
  }
  for (int q = 0; q < nc; q++) lab[q] = q;
  for (int cy = 0; cy < ch; cy++) for (int cx = 0; cx < cw; cx++) {
    int q = cy*cw+cx; if (!cells[q]) continue;
    if (mark) for (int dy = -1; dy <= 1; dy++) for (int dx = -1; dx <= 1; dx++)
      if (cy+dy >= 0 && cy+dy < ch && cx+dx >= 0 && cx+dx < cw) hot[q+dy*cw+dx] = 1;
    int nb[4] = {cx>0?q-1:-1, cy>0?q-cw:-1, (cy>0&&cx>0)?q-cw-1:-1, (cy>0&&cx<cw-1)?q-cw+1:-1};
    for (int e = 0; e < 4; e++) if (nb[e] >= 0 && cells[nb[e]]) { int r1 = find(q), r2 = find(nb[e]); if (r1 != r2) lab[r1] = r2; }
  }
  int roots[MAXB], b0 = nbl;
  for (int q = 0; q < nc; q++) {
    if (!cells[q]) continue;
    int r = find(q), bi = -1;
    for (int e = b0; e < nbl; e++) if (roots[e] == r) { bi = e; break; }
    if (bi < 0) { if (nbl == MAXB) continue; bi = nbl++; roots[bi] = r; sx[bi]=sy[bi]=0; np[bi]=0; pk[bi]=0; sat[bi]=0; }
    int x0 = (q%cw)*16, y0 = (q/cw)*16;
    for (int y = y0; y < y0+16; y++) {
      const u8 *row = scratch + (size_t)y*w;
      for (int x = x0; x < x0+16; x++) {
        int v = row[x]; if (v > pk[bi]) pk[bi] = v;
        if (v < T) continue;
        sx[bi] += (u64)(x*10+5); sy[bi] += (u64)(y*10+5); np[bi]++;
        if (v >= 250) sat[bi]++;
      }
    }
  }
  return nbl;
}

/* ---------------------------------------------------------------- tracks */
static struct track *track_new(struct cam *c, float x, float y, double t) {
  for (int i = 0; i < MAXTRACKS; i++) {
    struct track *tr = &c->tr[i];
    if (tr->used) continue;
    memset(tr, 0, sizeof *tr);
    tr->used = 1; tr->id = c->next_id++; tr->state = T_NEW;
    tr->x = x; tr->y = y; tr->last_hit = t;
    return tr;
  }
  /* recycle the stalest NEW track when full: actives are never evicted here */
  struct track *best = NULL;
  for (int i = 0; i < MAXTRACKS; i++) {
    struct track *tr = &c->tr[i];
    if (tr->state != T_ACTIVE && (!best || tr->last_hit < best->last_hit)) best = tr;
  }
  if (best) { memset(best, 0, sizeof *best); best->used = 1; best->id = c->next_id++; best->x = x; best->y = y; best->last_hit = t; }
  return best;
}

/* Goertzel magnitude of x[0..n) at normalized frequency f (cycles/sample) */
static float goertzel(const float *x, int n, float f) {
  float w = 2.0f * (float)M_PI * f;
  float cw = cosf(w), sw = sinf(w), coeff = 2.0f * cw;
  float s0 = 0, s1 = 0, s2 = 0;
  for (int i = 0; i < n; i++) { s0 = x[i] + coeff*s1 - s2; s2 = s1; s1 = s0; }
  (void)sw;
  return sqrtf(s1*s1 + s2*s2 - coeff*s1*s2);
}

/* long-frame sample rate from the recent long-frame timestamps (median gap);
 * 0 until enough are seen. *reliable = 0 when the gaps are wildly spread —
 * frames are being dropped, and a spectrum over dropped-frame samples is
 * garbage (bin frequencies alias unpredictably), so classification must wait. */
static float longfps(struct cam *c, int *reliable) {
  double g[31]; int n = 0;
  for (int i = 1; i < c->nlfts && n < 31; i++) {
    double d = c->lfts[i] - c->lfts[i-1];
    if (d > 0.005 && d < 0.2) g[n++] = d;
  }
  if (n < 8) { if (reliable) *reliable = 0; return 0; }
  for (int i = 1; i < n; i++) { double v = g[i]; int j = i-1; while (j >= 0 && g[j] > v) { g[j+1] = g[j]; j--; } g[j+1] = v; }
  double med = g[n/2];
  if (reliable) *reliable = (g[n-1] - g[0]) < 0.35 * med;
  return (float)(1.0/med);
}

/* spectral test of one window of the amplitude ring; fills f0/q/pres and
 * returns 1 when it looks like a lighthouse dot's beat. */
static int spectral_test(const u8 *amp, int n, float fs, float fmin, float fstep, int nf,
                         float q_min, float f_stab, int skip_dc,
                         float *f0_out, float *q_out, float *pres_out) {
  static float win[768];
  float mean = 0; int nz = 0;
  for (int i = 0; i < n; i++) { mean += amp[i]; if (amp[i] >= PK_MIN) nz++; }
  mean /= n;
  for (int i = 0; i < n; i++) win[i] = (amp[i] - mean) * 0.5f * (1 - cosf(2.0f*(float)M_PI*i/(n-1)));
  float floor_sum = 0; int best = -1, nb = nf; float bestm = 0;
  while (nb > 1 && (fmin + fstep*(nb-1)) > 0.45f*fs) nb--;   /* bins above Nyquist alias: drop them */
  for (int k = skip_dc; k < nb; k++) {
    float m = goertzel(win, n, (fmin + fstep*k)/fs);   /* cycles per sample */
    floor_sum += m;
    if (m > bestm) { bestm = m; best = k; }
  }
  if (best < (skip_dc ? 1 : 0)) return 0;
  float floormed = floor_sum/nb;
  if (floormed < 1e-3f) return 0;
  float q = bestm/floormed;
  float f0 = fmin + fstep*best;
  /* stability: peak frequency of the first vs second half of the window */
  int h = n/2;
  float b1 = 0, b2 = 0; int k1 = -1, k2 = -1;
  for (int k = skip_dc; k < nb; k++) {
    float f = (fmin + fstep*k)/fs;
    float m1 = goertzel(win, h, f);
    float m2 = goertzel(win+h, n-h, f);
    if (m1 > b1) { b1 = m1; k1 = k; }
    if (m2 > b2) { b2 = m2; k2 = k; }
  }
  float pres = (float)nz/n;
  int stable = k1 >= 0 && k2 >= 0 && fabsf(fstep*k1 - fstep*k2) <= f_stab && fabsf(fstep*k1 - (f0 - fmin)) <= f_stab;
  *f0_out = f0; *q_out = q; *pres_out = pres;
  int nyquist = fabsf(f0 - fs/2) < 0.3f;   /* perfect alternation = frame artifact */
  return q >= q_min && stable && pres >= PRES_MIN && !nyquist;
}

/* flicker-mode switches step the long-frame rate between ~30.0 and ~33.3 Hz:
 * samples taken across the switch don't share a spectrum, so drop everything
 * and relearn instead of classifying through the seam. Dropped frames can also
 * halve the measured rate (15 Hz) — that is NOT a mode switch, so only rates
 * that look like a real mode, confirmed on consecutive checks, flush. */
static void check_rate(struct cam *c) {
  float fs = longfps(c, 0);
  if (fs <= 0) return;
  int plausible = (fs > 28.5f && fs < 31.5f) || (fs > 32.3f && fs < 34.3f);
  int changed = c->last_fs && fabsf(fs - c->last_fs) > 1.5f;
  if (!changed) { c->last_fs = fs; c->cand_n = 0; return; }
  int confirmed = 0;
  if (plausible && c->cand_fs == fs) confirmed = ++c->cand_n >= 2;
  else { c->cand_fs = fs; c->cand_n = 1; }
  if (!confirmed) return;
  for (int i = 0; i < MAXTRACKS; i++) {
    struct track *tr = &c->tr[i];
    if (tr->used && tr->state == T_ACTIVE)
      out("I cam%d track %d dropped (frame rate changed)", (int)(c-C), tr->id);
    tr->used = 0;
  }
  c->nlfts = 0;
  out("I cam%d frame rate changed %.2f -> %.2f Hz (flicker mode switch?): relearning", (int)(c-C), c->last_fs, fs);
  c->last_fs = fs;
  c->cand_n = 0;
}

static void analyse(struct cam *c, struct track *tr) {
  if (tr->n < DFT_MIN || tr->hits < MIN_HITS) return;
  int reliable = 0;
  float fs = longfps(c, &reliable);
  if (fs < 25.0f || fs > 36.0f || !reliable) return;   /* 30.0 / 33.33 Hz modes only, clean timing */
  if (!c->nlfts_reported) {
    c->nlfts_reported = 1;
    out("I cam%d long frames every %.2f ms (%.2f Hz sample rate)", (int)(c-C), 1000.0f/fs, fs);
  }
  float f0 = 0, q = 0, pres = 0, qfail = 0;
  int good = 0; const char *how = "";
  /* fast window (recent FAST_WIN samples): the 10-12 Hz beats of 60 Hz mode */
  {
    int off = tr->n > FAST_WIN ? tr->n - FAST_WIN : 0;
    float a, b, p;
    if (spectral_test(tr->amp + off, tr->n - off, fs, FMIN, FSTEP, NFREQ, Q_MIN, F_STAB, 1, &a, &b, &p)) {
      good = 1; how = "fast"; f0 = a; q = b; pres = p;
    } else qfail = b;
  }
  /* long window (whole ring): 50 Hz mode aliases the beat to 0.1-2.2 Hz,
   * which needs the full ~23 s to resolve */
  if (!good && tr->n >= LONG_MIN) {
    float a, b, p;
    if (spectral_test(tr->amp, tr->n, fs, FMIN_SLOW, FSTEP_SLOW, NFREQ_SLOW, Q_MIN, F_STAB_SLOW, 0, &a, &b, &p)) {
      good = 1; how = "slow"; f0 = a; q = b; pres = p;
    } else if (b > qfail) qfail = b;
  }
  tr->f0 = f0; tr->q = q ? q : qfail; tr->pres = pres;
  if (good && tr->state != T_ACTIVE) {
    tr->state = T_ACTIVE;
    tr->since_classified = now_s(CLOCK_MONOTONIC);
    out("I cam%d track %d CLASSIFIED lighthouse-dot (%s): f0 %.2f Hz q %.1f pres %.2f at (%.1f,%.1f)",
        (int)(c-C), tr->id, how, f0, q, pres, tr->x, tr->y);
  } else if (!good && tr->state == T_ACTIVE) {
    /* keep emitting through brief spectrum loss; deep loss drops it */
    if (qfail < Q_MIN*0.5f) {
      tr->state = T_NEW;
      out("I cam%d track %d lost the beat (q %.1f)", (int)(c-C), tr->id, qfail);
    }
  }
  out("B %d %d %.1f %.1f %.2f %.1f %.2f %d %d %.2f %s", (int)(c-C), tr->id, tr->x, tr->y, f0, q ? q : qfail, pres, tr->hits, tr->state, fs, how);
}

static void track_update(struct cam *c, double t, const u64 *sx, const u64 *sy, const u32 *np, const int *pk, int nbl) {
  struct track *tr;
  u8 taken[MAXTRACKS] = {0};
  /* age out dead tracks */
  for (int i = 0; i < MAXTRACKS; i++) {
    tr = &c->tr[i];
    if (tr->used && t - tr->last_hit > TRACK_TTL) {
      if (tr->state == T_ACTIVE) out("I cam%d track %d timed out", (int)(c-C), tr->id);
      tr->used = 0;
    }
  }
  for (int b = 0; b < nbl; b++) {
    if (!np[b]) continue;
    float bx = sx[b]/(float)np[b]/10.0f, by = sy[b]/(float)np[b]/10.0f;
    u32 npx = np[b]; int pkv = pk[b];
    if ((int)npx > MAXNPX || pkv < PK_MIN) continue;
    struct track *best = NULL; float bd = RADIUS2;
    for (int i = 0; i < MAXTRACKS; i++) {
      tr = &c->tr[i];
      if (!tr->used || taken[i]) continue;
      float dx = tr->x-bx, dy = tr->y-by, d = dx*dx+dy*dy;
      if (d < bd) { bd = d; best = tr; }
    }
    if (!best) best = track_new(c, bx, by, t);
    if (!best) continue;
    int bi = (int)(best - c->tr);
    taken[bi] = 1;
    best->x += 0.3f*(bx-best->x); best->y += 0.3f*(by-best->y);
    best->last_hit = t; best->hits++;
    /* push amplitude sample (one per long frame) */
    if (best->n < DFT_WIN) { best->amp[best->n] = (u8)(pkv > 255 ? 255 : pkv); best->npx[best->n] = (u8)(npx > 255 ? 255 : npx); best->n++; }
    else { memmove(best->amp, best->amp+1, DFT_WIN-1); memmove(best->npx, best->npx+1, DFT_WIN-1); best->amp[DFT_WIN-1] = (u8)pkv; best->npx[DFT_WIN-1] = (u8)(npx>255?255:npx); }
    if (best->n == DFT_MIN || (best->n % DFT_STEP == 0 && best->n < DFT_WIN) ||
        (best->n == DFT_WIN && (++best->tick & 15) == 0)) analyse(c, best);
  }
  /* absent tracks get a 0 sample so presence stats work */
  for (int i = 0; i < MAXTRACKS; i++) {
    tr = &c->tr[i];
    if (!tr->used || taken[i]) continue;
    if (tr->n < DFT_WIN) { tr->amp[tr->n] = 0; tr->npx[tr->n] = 0; tr->n++; }
    else { memmove(tr->amp, tr->amp+1, DFT_WIN-1); memmove(tr->npx, tr->npx+1, DFT_WIN-1); tr->amp[DFT_WIN-1] = 0; tr->npx[DFT_WIN-1] = 0; }
    if (tr->n == DFT_MIN || (tr->n % DFT_STEP == 0 && tr->n < DFT_WIN) ||
        (tr->n == DFT_WIN && (++tr->tick & 15) == 0)) analyse(c, tr);
  }
}

/* ---------------------------------------------------------------- frames */
static int torn(struct cam *c, int i, const struct v4l2_buffer *b) {
  struct v4l2_buffer b2; struct v4l2_plane pl[VIDEO_MAX_PLANES];
  if (!querybuf(c, i, &b2, pl) && !(b2.flags & V4L2_BUF_FLAG_QUEUED) && b2.sequence == b->sequence) return 0;
  ntorn++;
  return 1;
}
/* motion gate: a coarse still/panning test from cam0's long frames. The dots
 * worth sending come from steady head poses — while panning, every dot's ray
 * is bent by the pose-timing error (12.9 ms +- 3.8 on this link: 30 deg/s of
 * turn bends a 3 m ray ~1.5 cm), and the bend is coherent, so the PC's
 * scatter gates can't see it. Hold dots back while panning. */
#define MDS       2      /* motion downsample: 1 ds px = 2 full-res px */
#define MR_RANGE  8      /* search +-8 ds px (about +-48 deg/s between frames) */
#define MR_STEP   2
#define MOTION_GATE 2    /* ds px of clear shift that counts as panning (~21 deg/s) */
static u8 mcur[2048*2048/(MDS*MDS)], mprev[2048*2048/(MDS*MDS)];
static int mhave, mpan_run, mstill_run, g_panning;

static int motion_est(int w, int h) {
  int dw = w/MDS, dh = h/MDS;
  if (dw > 2048/MDS || dh > 2048/MDS) return g_panning;
  for (int y = 0; y < dh; y++)
    for (int x = 0; x < dw; x++)
      mcur[y*dw+x] = scratch[((size_t)(y*MDS)+(MDS/2))*w + x*MDS];
  if (!mhave) { memcpy(mprev, mcur, (size_t)dw*dh); mhave = 1; return g_panning; }
  int s00 = 0;
  for (int y = 2; y < dh-2; y += 2)
    for (int x = 2; x < dw-2; x += 2)
      s00 += abs(mcur[y*dw+x] - mprev[y*dw+x]);
  int best = 1<<30, bx = 0, by = 0;
  for (int dy = -MR_RANGE; dy <= MR_RANGE; dy += MR_STEP)
    for (int dx = -MR_RANGE; dx <= MR_RANGE; dx += MR_STEP) {
      if (!dx && !dy) continue;
      int sum = 0;
      for (int y = 2; y < dh-2; y += 2)
        for (int x = 2; x < dw-2; x += 2) {
          int yy = y+dy, xx = x+dx;
          if ((unsigned)yy >= (unsigned)(dh-4) || (unsigned)xx >= (unsigned)(dw-4)) continue;
          sum += abs(mcur[y*dw+x] - mprev[yy*dw+xx]);
        }
      if (sum < best) { best = sum; bx = dx; by = dy; }
    }
  memcpy(mprev, mcur, (size_t)dw*dh);
  int n = ((dh-4)/2)*((dw-4)/2);
  if (n <= 0) return g_panning;
  float sad0 = (float)s00/n, sadm = (float)best/n;
  /* no clear shift, or a featureless frame: treat as still (the PC's fast
   * gate backstops); a clear minimum away from 0 is panning */
  int panning;
  if (sadm >= 0.95f*sad0 || sad0 < 2.0f) panning = 0;
  else panning = (abs(bx) >= MOTION_GATE || abs(by) >= MOTION_GATE);
  if (panning) {
    mpan_run++; mstill_run = 0;
    if (mpan_run >= 2 && !g_panning) { g_panning = 1; out("I motion: panning, dots held back"); }
  } else {
    mstill_run++; mpan_run = 0;
    if (mstill_run >= 2 && g_panning) { g_panning = 0; out("I motion: steady, dots resumed"); }
  }
  return g_panning;
}

static void frame(int k, int i, const struct v4l2_buffer *b) {
  struct cam *c = &C[k];
  double ts = b->timestamp.tv_sec + b->timestamp.tv_usec/1e6;
  if (ts_clock < 0) {
    double raw = now_s(CLOCK_MONOTONIC_RAW) - ts, mono = now_s(CLOCK_MONOTONIC) - ts;
    ts_clock = raw >= 0 && raw < 0.2 ? CLOCK_MONOTONIC_RAW : CLOCK_MONOTONIC;
    out("I frame timestamps CLOCK_%s", ts_clock == CLOCK_MONOTONIC_RAW ? "MONOTONIC_RAW" : "MONOTONIC");
  }
  double cnow = now_s(ts_clock), mnow = now_s(CLOCK_MONOTONIC);
  u64 t_us = (u64)((ts + (mnow-cnow))*1e6);
  cpu_access(c->fd[i], 0);
  const u8 *p = c->map[i];
  if (!c->stride) {
    char why[200]; why[0]=0;
    int v = pick_pitch(c, p); (void)why;
    cpu_access(c->fd[i], 1);
    if (!torn(c, i, b)) count_pitch(c, v);
    return;
  }
  int mean = smean(c, p);
  u32 kk = c->k++;
  int hi = c->m1 > c->m2 ? c->m1 : c->m2;
  c->m2 = c->m1; c->m1 = mean;
  /* long/short class: the bright frames of the exposure pair. hi is a slowly
   * decaying running max, so a few missed frames don't break the split. */
  if (c->rmax < hi) c->rmax = hi; else c->rmax = (int)(c->rmax * 0.97);
  int is_short = kk >= 2 && mean*10 < c->rmax*6;
  if (!is_short) {
    for (int y = 0; y < c->h; y++) memcpy(scratch + (size_t)y*c->w, p + (size_t)y*c->stride, c->w);
  }
  cpu_access(c->fd[i], 1);
  if (torn(c, i, b)) return;

  if (is_short) {
    /* short frame: keep the cadence visible to the PC, no blobs */
    out("F %d %u %llu %d -1", k, b->sequence, (unsigned long long)t_us, mean);
    return;
  }
  c->longs++;
  if (k == 0) motion_est(c->w, c->h);
  check_rate(c);
  if (g_panning) {
    /* hold dots back while the head pans: cadence only */
    out("F %d %u %llu %d -1", k, b->sequence, (unsigned long long)t_us, mean);
    return;
  }
  u64 t0 = now_ns();
  u64 sx[MAXB], sy[MAXB]; u32 np[MAXB]; int pk[MAXB]; u32 sat[MAXB];
  int T = mean*4/10; if (T < 64) T = 64;
  int nc = (c->w>>4)*(c->h>>4);
  memset(hot, 0, (size_t)nc);
  int nbl = blobs(c->w, c->h, 250, 1, 0, sx, sy, np, pk, sat);
  if (T < 250) nbl = blobs(c->w, c->h, T, 0, nbl, sx, sy, np, pk, sat);
  scan_ns += now_ns() - t0; nscan++;
  double tn = t_us/1e6;
  if (c->nlfts < 32) c->lfts[c->nlfts++] = tn;
  else { memmove(c->lfts, c->lfts+1, 31*sizeof *c->lfts); c->lfts[31] = tn; }
  check_rate(c);
  track_update(c, tn, sx, sy, np, pk, nbl);
  /* Brightness/intensity gate, then raw emission: the aperture saturates a
   * SOLID core of pixels, while laser spatter on walls is dimmer and sparser
   * and shares the rotor modulation (so beat tests can't reject it). Keep
   * saturated compact cores, brightest first; the PC's geometric gates do the
   * rest. The beat classifier stays up as B diagnostics only. */
  int cand[MAXB], nc2 = 0;
  for (int e = 0; e < nbl; e++) {
    if (pk[e] < 250 || (int)np[e] > 150 || sat[e] < 4) continue;
    cand[nc2++] = e;
  }
  for (int a = 1; a < nc2; a++) {
    int v = cand[a], b2 = a-1;
    while (b2 >= 0 && sat[cand[b2]] < sat[v]) { cand[b2+1] = cand[b2]; b2--; }
    cand[b2+1] = v;
  }
  char bl[1800] = {0}; int L2 = 0, n = 0;
  for (int e = 0; e < nc2 && n < 16; e++) {
    L2 += snprintf(bl+L2, sizeof bl-L2, " %lld %lld %u %d", (long long)sx[cand[e]]/np[cand[e]], (long long)sy[cand[e]]/np[cand[e]], np[cand[e]], pk[cand[e]]);
    n++;
  }
  out("F %d %u %llu %d %d%s", k, b->sequence, (unsigned long long)t_us, mean, n, bl);
}
static void poll_camera(int k) {
  struct cam *c = &C[k];
  int idx[NBUF], n = 0;
  struct v4l2_buffer b[NBUF]; struct v4l2_plane pl[NBUF][VIDEO_MAX_PLANES];
  for (int i = 0; i < c->nb; i++) {
    if (querybuf(c, i, &b[i], pl[i])) { out("W %s QUERYBUF: %s", c->dev, strerror(errno)); exit(4); }
    if (pl[i][0].m.fd != c->xfd[i]) { out("W %s: XRService replaced its buffers", c->dev); exit(4); }
    if (!(b[i].flags & V4L2_BUF_FLAG_QUEUED) && b[i].sequence != c->seq[i]) idx[n++] = i;
  }
  for (int a = 1; a < n; a++) for (int j = a; j > 0 && b[idx[j]].sequence < b[idx[j-1]].sequence; j--)
    { int t = idx[j]; idx[j] = idx[j-1]; idx[j-1] = t; }
  for (int j = 0; j < n; j++) { c->seq[idx[j]] = b[idx[j]].sequence; frame(k, idx[j], &b[idx[j]]); }
  if (n) { c->prev2 = c->prev; c->prev = c->last; c->last = now_s(CLOCK_MONOTONIC); }
}
static double due(const struct cam *c, double t) {
  if (!c->last || t - c->last > 0.25) return c->polled + 0.02;
  double gap = c->prev2 ? c->prev - c->prev2 : 0;
  double d = c->last + (gap > 0.0075 && gap < 0.1 ? gap - 0.0015 : 0.006);
  return d > c->polled + 0.001 ? d : c->polled + 0.001;
}

int lhsight_main(void) {
  const char *runenv = getenv("LHSIGHT2_RUN");
  double run = runenv ? atof(runenv) : 0.0;   /* 0 = forever (the service) */
  double t_start = now_s(CLOCK_MONOTONIC);
  pid_t xr = xrservice_pid();
  if (!xr) { out("E XRService isn't running"); return 2; }
  if (camera_nodes(xr) != NCAM) { out("E XRService's log doesn't name %d cameras", NCAM); return 2; }
  int maxw = 0, maxh = 0;
  for (int k = 0; k < NCAM; k++) {
    if (open_camera(&C[k])) return 2;
    if (C[k].w > maxw) maxw = C[k].w;
    if (C[k].h > maxh) maxh = C[k].h;
  }
  if (borrow(xr)) return 2;
  int pidfd = (int)syscall(SYS_pidfd_open, xr, 0);
  scratch = aligned_alloc(64, (size_t)maxw*maxh);
  cells = malloc((size_t)(maxw/16)*(maxh/16)*sizeof *cells);
  lab = malloc((size_t)(maxw/16)*(maxh/16)*sizeof *lab);
  hot = malloc((size_t)(maxw/16)*(maxh/16));
  for (int k = 0; k < NCAM; k++)
    out("I lhsight2: cam%d %s %dx%d stride %d", k, C[k].dev, C[k].w, C[k].h, C[k].stride);
  double next_beat = now_s(CLOCK_MONOTONIC) + 1;
  for (;;) {
    double t = now_s(CLOCK_MONOTONIC), wake = t + 0.02;
    for (int k = 0; k < NCAM; k++) {
      struct cam *c = &C[k];
      if (t >= due(c, t)) { poll_camera(k); c->polled = t = now_s(CLOCK_MONOTONIC); }
      double d = due(c, t); if (d < wake) wake = d;
    }
    t = now_s(CLOCK_MONOTONIC);
    if (run > 0 && t - t_start >= run) { out("I done"); return 0; }
    if (t >= next_beat) {
      next_beat = t + 1;
      struct pollfd pf = {pidfd, POLLIN, 0};
      if (pidfd >= 0 && poll(&pf, 1, 0) == 1) { out("W XRService exited"); return 3; }
      out("T %llu frames %u scans %u scan_us %.0f torn %u polls %u",
          (unsigned long long)(t*1e6), (unsigned)(C[0].k+C[1].k+C[2].k+C[3].k), nscan,
          nscan ? scan_ns/1e3/nscan : 0.0, ntorn, npoll);
    }
    if (wake > t) usleep((useconds_t)((wake-t)*1e6));
  }
}
