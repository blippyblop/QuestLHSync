// qlhs_presence: ground truth for the direction-cell presence rings, straight from raw pixels.
//
// The driver's presence rings (SpotCell::seqL / seqS) are meant to answer, for one camera's
// ~2.4 deg quantised direction: "of the frames in which this direction was in the image, on
// how many did a blob land in it?" A real base station's answer is the L figure below, and
// the S figure is the S-gate's evidence. This computes both independently of the driver --
// same optics (so the same unprojection and extrinsics), but no cells, no rings, no solver --
// so the driver's log can be checked against it.
//
// The station's place in the ray space comes from the recording's S lines mapped by the
// accepted x: x(S) = Ry(yaw)*S + t. --station overrides it.
//
// usage: qlhs_presence <recording> --calib <calib.json> --yaw <deg> --t <x> <y> <z>
//                       [--station SERIAL=x,y,z]... [--dt <s>] [--no-scan]
//
// A frame's arrival time is not its exposure time, so the pose must be read at exposure:
// without --dt the tool scans the offset and reports at the one the long frames like best.
//
// The "in view" column is the OLD denominator (Sync::OnLine's dot(w, axis) > -0.2, a ~101 deg
// cone about the optical axis); "in image" is the tighter one. The ratio of the two is how
// much the old rule dilutes the ring.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "../driver/sync.h"

namespace {

struct Pose { double t; M3 R; V3 p; };
struct Rec { double pc; int cam; char cls; std::vector<std::pair<double, double>> bl; };

// the driver's cell: camera + the direction quantised to ~2.4 deg (Sync::CellKey)
uint64_t cell_key(int cam, V3 w) {
  uint64_t h = 1469598103934665603ull;
  auto mix = [&](uint64_t v) { h ^= v + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2); };
  mix((uint64_t)(uint32_t)cam);
  mix((uint64_t)(int64_t)std::lround(w.x * 24));
  mix((uint64_t)(int64_t)std::lround(w.y * 24));
  mix((uint64_t)(int64_t)std::lround(w.z * 24));
  return h;
}

struct Tally {
  long considered = 0, inview = 0, inimage = 0, incell = 0, decoy = 0;
  long px3 = 0, px6 = 0, px12 = 0, px25 = 0;
};

// per FIXED cell key: what a ring for that key would read. The ring is per key, not per
// station: the frames it covers are the ones that key was in image on, and the station is
// only IN that key while the head sits in the one pocket of places that puts it there.
struct Key { long Ln = 0, Lhit = 0, Sn = 0, Shit = 0; };

}  // namespace

int main(int argc, char **argv) {
  std::string log, calib, serial;
  double yawdeg = 0, dt = 0, tvec[3] = {0, 0, 0}, from = -1e18, to = 1e18;
  bool have_t = false, scan = true, have_dt = false;
  std::map<std::string, V3> over;
  for (int i = 1; i < argc; i++) {
    std::string a = argv[i];
    if (a == "--calib" && i + 1 < argc) calib = argv[++i];
    else if (a == "--yaw" && i + 1 < argc) yawdeg = atof(argv[++i]);
    else if (a == "--dt" && i + 1 < argc) { dt = atof(argv[++i]); have_dt = true; }
    else if (a == "--no-scan") scan = false;
    else if (a == "--from" && i + 1 < argc) from = atof(argv[++i]);
    else if (a == "--to" && i + 1 < argc) to = atof(argv[++i]);
    else if (a == "--t" && i + 3 < argc) { tvec[0] = atof(argv[++i]); tvec[1] = atof(argv[++i]); tvec[2] = atof(argv[++i]); have_t = true; }
    else if (a == "--station" && i + 1 < argc) {
      std::string v = argv[++i];
      size_t eq = v.find('=');
      if (eq == std::string::npos) { fprintf(stderr, "--station wants SERIAL=x,y,z\n"); return 2; }
      double x, y, z;
      if (sscanf(v.c_str() + eq + 1, "%lf,%lf,%lf", &x, &y, &z) != 3) { fprintf(stderr, "--station wants SERIAL=x,y,z\n"); return 2; }
      over[v.substr(0, eq)] = V3{x, y, z};
    } else if (!a.empty() && a[0] == '-') { fprintf(stderr, "unknown option %s\n", a.c_str()); return 2; }
    else log = a;
  }
  if (log.empty() || calib.empty()) {
    fprintf(stderr, "usage: qlhs_presence <recording> --calib <file> --yaw <deg> --t <x> <y> <z> [--station S=x,y,z] [--dt s]\n");
    return 2;
  }
  if (!have_t && over.empty()) { fprintf(stderr, "need --t (and --yaw) or --station for every station\n"); return 2; }

  std::ifstream cf(calib, std::ios::binary);
  std::stringstream cs;
  cs << cf.rdbuf();
  Optics opt;
  std::string err;
  if (!opt.Load(cs.str(), &err)) { fprintf(stderr, "calibration: %s\n", err.c_str()); return 1; }
  V3 tt{tvec[0], tvec[1], tvec[2]};
  double yaw = yawdeg / kDeg;

  std::map<std::string, V3> Sraw;
  std::vector<Pose> poses;
  std::vector<Rec> frames;
  std::ifstream f(log.c_str());
  if (!f) { fprintf(stderr, "can't open %s\n", log.c_str()); return 1; }
  std::string line;
  while (std::getline(f, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    char *sp = strchr(const_cast<char *>(line.c_str()), ' ');
    if (!sp || sp[1] == 0 || sp[2] != ' ') continue;
    double pc = atoll(line.c_str()) / 1e9;
    char kind = sp[1];
    const char *rest = sp + 3;
    if (kind == 'P') {
      int idx; double m[12];
      if (sscanf(rest, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &idx, &m[0], &m[1], &m[2], &m[3], &m[4],
                 &m[5], &m[6], &m[7], &m[8], &m[9], &m[10], &m[11]) == 13 && idx == 0) {
        M3 R;
        for (int r = 0; r < 3; r++)
          for (int c = 0; c < 3; c++) R.m[r][c] = m[r * 4 + c];
        poses.push_back({pc, R, V3{m[3], m[7], m[11]}});
      }
    } else if (kind == 'S') {
      char sn[64]; double v[7]; int ch = 0;
      if (sscanf(rest, "%63s %lf %lf %lf %lf %lf %lf %lf %d", sn, &v[0], &v[1], &v[2], &v[3], &v[4], &v[5], &v[6], &ch) >= 8)
        Sraw[sn] = V3{v[0], v[1], v[2]};
    } else if (kind == 'F') {
      int cam, nb; unsigned k; long long tus; int mean, off = 0;
      if (sscanf(rest, "%d %u %lld %d %d%n", &cam, &k, &tus, &mean, &nb, &off) < 5) continue;
      size_t L = strlen(rest);
      while (L && rest[L - 1] == ' ') L--;
      char cls = (L >= 2 && rest[L - 2] == ' ' && (rest[L - 1] == 'L' || rest[L - 1] == 'S')) ? rest[L - 1] : '.';
      if (cls == '.') continue;   // untagged logs carry no class: no ground truth to state for them
      Rec r;
      r.pc = pc; r.cam = cam; r.cls = cls;
      const char *s = rest + off;
      for (int i = 0; i < nb; i++) {
        int x10, y10, npx, peak, used = 0;
        if (sscanf(s, " %d %d %d %d%n", &x10, &y10, &npx, &peak, &used) < 4) break;
        s += used;
        r.bl.push_back({x10 / 10.0 - 0.5, y10 / 10.0 - 0.5});
      }
      frames.push_back(std::move(r));
    }
  }
  std::map<std::string, V3> S;
  for (auto &kv : Sraw) S[kv.first] = Ry(yaw) * kv.second + tt;
  for (auto &kv : over) S[kv.first] = kv.second;
  if (S.empty()) { fprintf(stderr, "no station in the log and none given\n"); return 1; }

  // --from/--to are seconds from the first line: the driver's sweeps are every 30 s, so a
  // 20 s tail is what its PresenceAt window holds. Compare like with like.
  double f0 = frames.empty() ? 0 : frames.front().pc;
  std::vector<const Rec *> byL, byS;
  long skipped = 0;
  for (auto &r : frames) {
    double rel = r.pc - f0;
    if (rel < from || rel > to) continue;
    // the SAME cameras the driver runs: with the upper-camera cut on it never even sees the lower
    // pair's frames (OnLine drops them), so pooling all four here would not be the same statistic
    if (!opt.Used(r.cam)) { skipped++; continue; }
    (r.cls == 'L' ? byL : byS).push_back(&r);
  }
  if (skipped) printf("left out %ld frames from the cameras the driver skips\n", skipped);
  if (from > -1e17 || to < 1e17) printf("window %.0f..%.0f s of the log: %zu L, %zu S frames\n", from, to, byL.size(), byS.size());

  // one pass: for every frame, every station, is the station in image, in the old cone, and did a blob
  // land in its cell? -- decoy = the same frame with the direction turned 10 deg about the optical axis.
  // the camera's forward model, the mirror of Optics::InImage (the tool may not touch the driver's
  // private CamCal, so it re-derives the pixel from the intrinsics it parsed itself)
  auto proj = [&](int cam, V3 dc, double &px, double &py) {
    if (!(dc.z > 0.05)) return false;
    return opt.Project(cam, dc, px, py);
  };
  auto pass = [&](double d, std::map<std::string, std::pair<Tally, Tally>> &out) {
    out.clear();
    for (int half = 0; half < 2; half++) {
      for (const Rec *rp : (half ? byS : byL)) {
        const Rec &r = *rp;
        double want = r.pc + d;
        size_t lo = 0, hi = poses.size();
        while (lo < hi) { size_t mid = (lo + hi) / 2; if (poses[mid].t < want) lo = mid + 1; else hi = mid; }
        if (lo == 0 || lo >= poses.size()) continue;
        const Pose &P = (std::fabs(poses[lo].t - want) < std::fabs(poses[lo - 1].t - want)) ? poses[lo] : poses[lo - 1];
        if (std::fabs(P.t - want) > 0.03) continue;
        M3 Rc; V3 tc;
        if (!opt.Pose(r.cam, Rc, tc)) continue;
        V3 C = P.p + P.R * tc;                 // the camera's place in the ray space
        V3 axis = P.R * (Rc * V3{0, 0, 1});    // its optical axis there
        std::vector<std::pair<uint64_t, std::pair<double, double>>> blobs;   // cell key + pixel (for the decoy)
        blobs.reserve(r.bl.size());
        for (auto &b : r.bl) {
          V3 o, dd;
          if (!opt.Ray(r.cam, b.first, b.second, o, dd)) continue;
          blobs.push_back({cell_key(r.cam, P.R * dd), b});
        }
        for (auto &kv : S) {
          V3 W = kv.second - C;
          double dist = norm(W);
          if (dist < 0.2) continue;
          V3 w = W * (1 / dist);
          Tally &t = (half ? out[kv.first].second : out[kv.first].first);
          t.considered++;
          if (dot(w, axis) > -0.2) t.inview++;
          V3 Pcam = T(Rc) * (T(P.R) * W);                       // station, camera frame
          if (!opt.InImage(r.cam, Pcam, 0)) continue;
          t.inimage++;
          uint64_t ck = cell_key(r.cam, w);
          bool hit = false;
          for (auto &bk : blobs)
            if (bk.first == ck) { hit = true; break; }
          if (hit) t.incell++;
          // nearest blob to the PROJECTED pixel (needs the camera's forward model)
          double px, py, d2 = 1e18;
          if (proj(r.cam, Pcam, px, py))
            for (auto &b : r.bl) {
              double dx = b.first - px, dy = b.second - py;
              if (dx * dx + dy * dy < d2) d2 = dx * dx + dy * dy;
            }
          else d2 = 1e18;
          double dpx = std::sqrt(d2);
          if (dpx <= 3) t.px3++;
          if (dpx <= 6) t.px6++;
          if (dpx <= 12) t.px12++;
          if (dpx <= 25) t.px25++;
          // decoy: the station's direction turned 10 deg about the camera's optical axis
          double g = 10.0 / kDeg;
          V3 Q{Pcam.x * std::cos(g) - Pcam.y * std::sin(g), Pcam.x * std::sin(g) + Pcam.y * std::cos(g), Pcam.z};
          V3 Wd = P.R * (Rc * Q);
          uint64_t dk = cell_key(r.cam, Wd * (1.0 / norm(Wd)));
          for (auto &bk : blobs)
            if (bk.first == dk) { t.decoy++; break; }
        }
      }
    }
  };

  // ---- per fixed key, the statistic a ring actually holds
  std::map<std::string, std::map<std::pair<int, uint64_t>, Key>> keys;
  auto keypass = [&](double d) {
    keys.clear();
    for (int half = 0; half < 2; half++) {
      for (const Rec *rp : (half ? byS : byL)) {
        const Rec &r = *rp;
        double want = r.pc + d;
        size_t lo = 0, hi = poses.size();
        while (lo < hi) { size_t mid = (lo + hi) / 2; if (poses[mid].t < want) lo = mid + 1; else hi = mid; }
        if (lo == 0 || lo >= poses.size()) continue;
        const Pose &P = (std::fabs(poses[lo].t - want) < std::fabs(poses[lo - 1].t - want)) ? poses[lo] : poses[lo - 1];
        if (std::fabs(P.t - want) > 0.03) continue;
        M3 Rc; V3 tc;
        if (!opt.Pose(r.cam, Rc, tc)) continue;
        V3 C = P.p + P.R * tc;
        std::vector<uint64_t> blobs;
        for (auto &b : r.bl) {
          V3 o, dd;
          if (opt.Ray(r.cam, b.first, b.second, o, dd)) blobs.push_back(cell_key(r.cam, P.R * dd));
        }
        for (auto &kv : S) {
          V3 W = kv.second - C;
          double dist = norm(W);
          if (dist < 0.2) continue;
          V3 w = W * (1 / dist);
          V3 Pcam = T(Rc) * (T(P.R) * W);
          if (!opt.InImage(r.cam, Pcam, 0)) continue;
          uint64_t ck = cell_key(r.cam, w);
          Key &K = keys[kv.first][{r.cam, ck}];
          bool on = false;
          for (uint64_t b : blobs) if (b == ck) { on = true; break; }
          if (half) { K.Sn++; if (on) K.Shit++; } else { K.Ln++; if (on) K.Lhit++; }
        }
      }
    }
  };

  std::map<std::string, std::pair<Tally, Tally>> best;
  if (!have_dt) {
    double bd = 0, bh = -1;
    std::map<std::string, std::pair<Tally, Tally>> cur;
    for (double d = -0.15; d <= 0.1501; d += 0.005) {
      pass(d, cur);
      double h = 0;
      for (auto &kv : cur) if (kv.second.first.inimage) h += (double)kv.second.first.incell / kv.second.first.inimage;
      if (h > bh) { bh = h; bd = d; best = cur; }
    }
    if (!scan) { /* still scanned to place the pose, as the pose is only readable at exposure */ }
    dt = bd;
    printf("pose offset: %+.0f ms (the long frames' best)\n", dt * 1000);
    pass(dt, best);
  } else {
    pass(dt, best);
  }

  printf("\n%-18s %s   %-9s %-9s %-9s %-9s %-9s\n", "station", "cls", "frames", "in view", "in image", "in cell", "decoy 10d");
  printf("%-18s %s   %-9s %-9s %-9s %-9s %-9s\n", "", "", "(old denom)", "(tighter)", "(= the ring)", "(control)");
  for (auto &kv : best) {
    const Tally *t[2] = {&kv.second.first, &kv.second.second};
    for (int c = 0; c < 2; c++) {
      const Tally &x = *t[c];
      if (!x.considered) continue;
      printf("%-18s %-3s %9ld %8.0f%% %9ld %8.1f%% %8.1f%% %9.1f%%  |  %5.1f%% %5.1f%% %5.1f%% %5.1f%%\n", kv.first.c_str(),
             c ? "S" : "L", x.considered, 100.0 * x.inview / x.considered, x.inimage, 100.0 * x.inimage / x.considered,
             100.0 * x.incell / x.inimage, 100.0 * x.decoy / x.inimage, 100.0 * x.px3 / x.inimage, 100.0 * x.px6 / x.inimage,
             100.0 * x.px12 / x.inimage, 100.0 * x.px25 / x.inimage);
    }
  }
  keypass(dt);
  printf("\nper FIXED cell key (>=32 L and >=32 S frames) -- what a ring for that key actually reads:\n");
  for (auto &kv : keys) {
    std::vector<double> lr, sr;
    for (auto &e : kv.second) {
      const Key &K = e.second;
      if (K.Ln < 32 || K.Sn < 32) continue;
      lr.push_back(100.0 * K.Lhit / K.Ln);
      sr.push_back(100.0 * K.Shit / K.Sn);
    }
    if (lr.empty()) { printf("  %-16s no key with enough frames\n", kv.first.c_str()); continue; }
    std::sort(lr.begin(), lr.end());
    std::sort(sr.begin(), sr.end());
    printf("  %-16s keys %3zu   L: median %4.0f%%  p90 %4.0f%%  max %4.0f%%   |   S: median %4.0f%%  p90 %4.0f%%  max %4.0f%%\n",
           kv.first.c_str(), lr.size(), lr[lr.size() / 2], lr[(size_t)(lr.size() * 0.9)], lr.back(), sr[sr.size() / 2],
           sr[(size_t)(sr.size() * 0.9)], sr.back());
  }
  printf("\n\"in cell\" is the station's own cell frame by frame; a station ROLLS across cells as the head moves,\n");
  printf("\"decoy\" is the same frames with the direction turned 10 deg: the clutter's own rate.\n");
  return 0;
}
