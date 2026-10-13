#include "sync.h"

#include <algorithm>
#include <cctype>
#include <cstdarg>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <unordered_map>

#include "json.h"

static const double kInf = std::numeric_limits<double>::infinity();

static std::string Fmt(const char *fmt, ...) {
  char buf[1024];
  va_list ap;
  va_start(ap, fmt);
  vsnprintf(buf, sizeof buf, fmt, ap);
  va_end(ap);
  return buf;
}

static double Median(std::vector<double> v) {
  if (v.empty()) return 0;
  size_t n = v.size(), h = n / 2;
  std::nth_element(v.begin(), v.begin() + h, v.end());
  double hi = v[h];
  if (n & 1) return hi;
  double lo = *std::max_element(v.begin(), v.begin() + h);
  return (lo + hi) / 2;
}

static bool ReadFile(const std::string &path, std::string &out) {
  std::ifstream f(path, std::ios::binary);
  if (!f) return false;
  std::stringstream ss;
  ss << f.rdbuf();
  out = ss.str();
  return true;
}

static bool WriteFileAtomic(const std::string &path, const std::string &text) {
  std::string tmp = path + ".tmp";
  {
    std::ofstream f(tmp, std::ios::binary | std::ios::trunc);
    if (!f) return false;
    f << text;
    if (!f) return false;
  }
  std::remove(path.c_str());
  return std::rename(tmp.c_str(), path.c_str()) == 0;
}

static std::string Now() {
  time_t t = time(nullptr);
  struct tm tm;
#ifdef _WIN32
  localtime_s(&tm, &t);
#else
  localtime_r(&t, &tm);
#endif
  char b[32];
  strftime(b, sizeof b, "%Y-%m-%d %H:%M:%S", &tm);
  return b;
}

// ================================================================ optics
void CamCal::Dist(double a, double b, double &u, double &v) const {
  double r = std::hypot(a, b), th = std::atan(r), th2 = th * th;
  double thd = th * (1 + th2 * (k[0] + th2 * (k[1] + th2 * (k[2] + th2 * (k[3] + th2 * (k[4] + th2 * k[5]))))));
  double s = r > 1e-12 ? thd / r : 1.0;
  double uu = a * s, vv = b * s, r2 = uu * uu + vv * vv, d = 2 * (uu * p[0] + vv * p[1]);
  u = uu + d * uu + r2 * p[0];
  v = vv + d * vv + r2 * p[1];
}

// Newton on the forward model
bool CamCal::Unproject(double px, double py, V3 &d) const {
  double tu = (px - cx) / f, tv = (py - cy) / fy;
  double rd = std::hypot(tu, tv), th = std::min(std::max(rd, 0.0), 1.5);
  double s = rd > 1e-12 ? std::tan(th) / rd : 1.0, a = tu * s, b = tv * s;
  const double h = 1e-6;
  for (int i = 0; i < 30; i++) {
    double u, v, ua, va, ub, vb;
    Dist(a, b, u, v);
    Dist(a + h, b, ua, va);
    Dist(a, b + h, ub, vb);
    double eu = u - tu, ev = v - tv;
    double J11 = (ua - u) / h, J21 = (va - v) / h, J12 = (ub - u) / h, J22 = (vb - v) / h;
    double det = J11 * J22 - J12 * J21;
    a -= (J22 * eu - J12 * ev) / det;
    b -= (-J21 * eu + J11 * ev) / det;
  }
  double n = std::sqrt(a * a + b * b + 1);
  if (!std::isfinite(n) || n <= 0) return false;
  d = {a / n, b / n, 1 / n};
  return true;
}

bool Optics::Load(const std::string &json, std::string *err) {
  JVal root;
  if (!JParse(json, root)) { if (err) *err = "calibration isn't JSON"; return false; }
  const JVal *kind = root.get("kind");
  if (kind && kind->t == JVal::Str && kind->s == "steam_frame") return LoadFrame(root, err);
  const JVal *cams = root.get("CameraCalibration");
  if (!cams || cams->t != JVal::Arr) { if (err) *err = "no CameraCalibration"; return false; }
  CamCal out[kMaxCam];
  int n = 0;
  for (auto &c : cams->a) {
    const JVal *id = c.get("Id"), *size = c.get("ImageSize"), *dfc = c.get("DeviceFromCamera");
    const JVal *proj = c.get("Projection"), *dist = c.get("Distortion"), *shutter = c.get("Shutter");
    const JVal *sensor = c.get("SensorType");
    if (!id || !size || !dfc || !proj || !dist) continue;
    int i = id->t == JVal::Str ? atoi(id->s.c_str()) : (int)id->num(-1);  // "Id": "2" in Meta's files
    if (i < 0 || i >= kMaxCam) continue;
    const JVal *type = shutter ? shutter->get("Type") : nullptr;
    if (type && type->t == JVal::Str && type->s == "Rolling") continue;  // the colour cameras: not tracking ones
    auto sz = size->nums(), T = dfc->nums();
    auto pc = proj->get("Coefficients") ? proj->get("Coefficients")->nums() : std::vector<double>();
    auto dc = dist->get("Coefficients") ? dist->get("Coefficients")->nums() : std::vector<double>();
    if (sz.size() < 2 || T.size() < 16 || pc.size() < 3 || dc.size() < 8) continue;
    CamCal &k = out[i];
    k.w = (int)sz[0]; k.h = (int)sz[1];
    for (int r = 0; r < 3; r++)
      for (int q = 0; q < 3; q++) k.R.m[r][q] = T[r * 4 + q];
    k.t = {T[3], T[7], T[11]};
    k.f = k.fy = pc[0]; k.cx = pc[1]; k.cy = pc[2];
    for (int j = 0; j < 6; j++) k.k[j] = dc[j];
    k.p[0] = dc[6]; k.p[1] = dc[7];
    k.ir = sensor && sensor->t == JVal::Str && sensor->s == "OV7251" ? kIR : 1.0;
    k.valid = true;
    n++;
  }
  if (!n) { if (err) *err = "no tracking camera calibration"; return false; }
  for (int i = 0; i < kMaxCam; i++) cams_[i] = out[i];
  const JVal *dev = root.get("Device"), *type = dev ? dev->get("DeviceType") : nullptr;
  device_ = type && type->t == JVal::Str ? type->s : "";
  exact_time_ = false;
  return true;
}

// a pose in the Frame's calibration files: the child frame's axes and origin, in its parent
static bool FramePose(const JVal *j, double scale, M3 &R, V3 &t) {
  if (!j) return false;
  const JVal *px = j->get("plus_x"), *pz = j->get("plus_z"), *pos = j->get("position");
  if (!px || !pz || !pos) return false;
  auto x = px->nums(), z = pz->nums(), p = pos->nums();
  if (x.size() != 3 || z.size() != 3 || p.size() != 3) return false;
  V3 X{x[0], x[1], x[2]}, Z{z[0], z[1], z[2]}, Y = cross(Z, X);
  for (int r = 0; r < 3; r++) { R.m[r][0] = X[r]; R.m[r][1] = Y[r]; R.m[r][2] = Z[r]; }
  t = V3{p[0], p[1], p[2]} * scale;
  return true;
}

// The Steam Frame's lhsyncd sends {"kind": "steam_frame", "xrservice": /persist/xrservice.json, "device_config":
// /persist/device_config.json}. xrservice.json: each tracking camera's KB4 intrinsics (Kannala-Brandt: Fisheye62 with
// only the 4 radial terms) and its pose against camera 0 (mm). device_config.json: camera 0's pose in the CAD frame
// (cv.cad_from_cal) and the head's (head), the frame SteamVR's HMD pose is of (m). Camera i's ray:
// head <- cad <- camera 0 <- camera i.
bool Optics::LoadFrame(const JVal &root, std::string *err) {
  const JVal *xr = root.get("xrservice"), *dc = root.get("device_config");
  const JVal *cams = xr ? xr->get("cameras") : nullptr, *cv = dc ? dc->get("cv") : nullptr;
  M3 Rcc, Rch;
  V3 tcc, tch;
  if (!cams || cams->t != JVal::Arr || !cv || !FramePose(cv->get("cad_from_cal"), 1.0, Rcc, tcc) ||
      !FramePose(dc->get("head"), 1.0, Rch, tch)) {
    if (err) *err = "not a Steam Frame calibration (xrservice.json cameras, device_config.json cv.cad_from_cal, head)";
    return false;
  }
  M3 Rhc = T(Rch) * Rcc;  // head <- camera 0
  V3 thc = T(Rch) * (tcc - tch);
  CamCal out[kMaxCam];
  int n = 0;
  for (size_t i = 0; i < cams->a.size() && i < (size_t)kMaxCam; i++) {
    const JVal &c = cams->a[i];
    const JVal *ex = c.get("extrinsics"), *in = c.get("intrinsics");
    const JVal *units = ex ? ex->get("__units") : nullptr;
    double scale = units && units->t == JVal::Str && units->s == "mm" ? 1e-3 : 1.0;
    M3 Rk;
    V3 tk;
    if (!FramePose(ex, scale, Rk, tk) || !in || in->t != JVal::Arr) continue;
    const JVal *kb = nullptr;
    for (auto &m : in->a)
      if (m.get("cameraModel") && m.get("cameraModel")->s == "kb") kb = &m;
    if (!kb || !c.get("width") || !c.get("height")) continue;
    auto num = [&](const char *k) { return kb->get(k) ? kb->get(k)->num(NAN) : NAN; };
    CamCal &k = out[i];
    k.w = (int)c.get("width")->num(); k.h = (int)c.get("height")->num();
    k.f = num("fx"); k.fy = num("fy"); k.cx = num("cx"); k.cy = num("cy");
    k.k[0] = num("k1"); k.k[1] = num("k2"); k.k[2] = num("k3"); k.k[3] = num("k4");
    bool finite = std::isfinite(k.f + k.fy + k.cx + k.cy + k.k[0] + k.k[1] + k.k[2] + k.k[3]) && k.f > 0 && k.fy > 0;
    if (!finite || k.w <= 0 || k.h <= 0) { k = CamCal(); continue; }
    k.R = Rhc * Rk;
    k.t = Rhc * tk + thc;
    k.valid = true;
    n++;
  }
  if (!n) {
    if (err) *err = "no usable camera in xrservice.json (KB4 intrinsics, extrinsics)";
    return false;
  }
  // The lower pair looks down and forward, at the floor and the hands; the base stations hang above. In two Frame
  // recordings only the two looking up ever matched a station, while the others sent most of the bright spots (lamps,
  // glints). They're left out, whatever their index (xrservice.json names them upper_left/upper_right, but the
  // geometry is what counts). If no camera qualifies the file is odd: all are used as before.
  int up = 0;
  for (int i = 0; i < kMaxCam; i++) up += out[i].valid && out[i].R.m[1][2] > kUpward;
  if (up)
    for (int i = 0; i < kMaxCam; i++) out[i].use = out[i].valid && out[i].R.m[1][2] > kUpward;
  for (int i = 0; i < kMaxCam; i++) cams_[i] = out[i];
  const JVal *model = dc->get("model_number");
  device_ = model && model->t == JVal::Str ? model->s : "Steam Frame";
  exact_time_ = true;
  return true;
}

std::string Optics::Unused() const {
  if (!upper_only_) return "";
  std::string s;
  for (int i = 0; i < kMaxCam; i++)
    if (cams_[i].valid && !cams_[i].use) s += (s.empty() ? "" : ", ") + std::to_string(i);
  return s;
}

bool Optics::Ray(int cam, double x, double y, V3 &o, V3 &d) const {
  if (cam < 0 || cam >= kMaxCam || !cams_[cam].valid) return false;
  const CamCal &c = cams_[cam];
  double k = c.ir;
  V3 dc;
  if (!c.Unproject(c.cx + (x - c.cx) * k, c.cy + (y - c.cy) * k, dc)) return false;
  d = c.R * dc;
  double n = norm(d);
  if (!std::isfinite(n) || n <= 0) return false;
  d = d * (1 / n);
  o = c.t;
  return true;
}

bool Optics::Pose(int cam, M3 &R, V3 &t) const {
  if (cam < 0 || cam >= kMaxCam || !cams_[cam].valid) return false;
  R = cams_[cam].R;
  t = cams_[cam].t;
  return true;
}

// a camera-frame direction -> the pixel it lands on, or false if it is off the image (or behind)
bool Optics::Project(int cam, V3 dc, double &px, double &py) const {
  if (cam < 0 || cam >= kMaxCam || !cams_[cam].valid || !(dc.z > 0.05)) return false;
  const CamCal &c = cams_[cam];
  double u, v;
  c.Dist(dc.x / dc.z, dc.y / dc.z, u, v);
  px = c.cx + c.f * u / c.ir;
  py = c.cy + c.fy * v / c.ir;
  return px >= 0 && py >= 0 && px <= c.w && py <= c.h;
}

bool Optics::InImage(int cam, V3 dc, double margin) const {
  double x, y;
  if (!Project(cam, dc, x, y)) return false;
  const CamCal &c = cams_[cam];
  return x >= margin && y >= margin && x <= c.w - margin && y <= c.h - margin;
}

// ================================================================ frame grid, clock
static double PyMod(double a, double m) {
  double r = std::fmod(a, m);
  return r < 0 ? r + m : r;
}

// the widest empty stretch of the detection times folded at P: the poll lags fill the rest, so it ends at the grid
static double FoldGap(const std::deque<double> &h, double P, double *grid) {
  std::vector<double> ph;
  ph.reserve(h.size());
  for (double v : h) ph.push_back(PyMod(v, P));
  std::sort(ph.begin(), ph.end());
  size_t n = ph.size(), best = 0;
  double bg = -1;
  for (size_t i = 0; i < n; i++) {
    double g = (i + 1 < n ? ph[i + 1] : ph[0] + P) - ph[i];
    if (g > bg) { bg = g; best = i; }
  }
  if (grid) *grid = ph[(best + 1) % n];
  return bg;
}

// how tightly the detection times fold at P: the narrowest phase band holding 90% of them (stragglers, a log flush
// tens of ms late, don't count)
static double FoldSpread(const std::deque<double> &h, double P) {
  std::vector<double> ph;
  ph.reserve(h.size() * 2);
  for (double v : h) ph.push_back(PyMod(v, P));
  std::sort(ph.begin(), ph.end());
  size_t n = ph.size(), m = (n * 9 + 9) / 10;
  for (size_t i = 0; i < n; i++) ph.push_back(ph[i] + P);
  double w = P;
  for (size_t i = 0; i < n; i++) w = std::min(w, ph[i + m - 1] - ph[i]);
  return w;
}

// the short frames' period from their detection times (each a few ms after its frame): near the median gap (or the
// period found before), the one they fold at most tightly. 0 if none folds cleanly (frame times too irregular to use)
double FrameGrid::Learn(const std::deque<double> &h, double prev) {
  std::vector<double> d;
  for (size_t i = 1; i < h.size(); i++)
    if (h[i] - h[i - 1] > 0.005 && h[i] - h[i - 1] < 0.5) d.push_back(h[i] - h[i - 1]);
  if (d.size() < 20) return 0;
  double p0 = Median(d), span = h.back() - h.front();
  if (span < 20 * p0) return 0;
  // a step of the period moves the window's far end span/P times as much: steps of 2 ms there, then of 0.05 ms
  double best = p0, bw = kInf, coarse = 0.002 * p0 / span, fine = 0.00005 * p0 / span;
  auto check = [&](double p) {
    double w = FoldSpread(h, p);
    if (w < bw) { bw = w; best = p; }
  };
  double lo = prev > 0 ? 0.995 * prev : 0.9 * p0, hi = prev > 0 ? 1.005 * prev : 1.1 * p0;
  for (double p = lo; p <= hi; p += coarse) check(p);
  double c = best;
  for (double p = c - coarse; p <= c + coarse; p += fine) check(p);
  return bw <= 0.5 * best ? best : 0;
}

double FrameGrid::period(int cam) const {
  auto it = c_.find(cam);
  return it == c_.end() ? 0 : it->second.p;
}

bool FrameGrid::Lag(int cam, double t, double &lag) {
  auto &c = c_[cam];
  auto &h = c.h;
  h.push_back(t);
  while (!h.empty() && h.front() < t - WIN) h.pop_front();
  auto &lh = c.longer;
  lh.push_back(t);
  while (!lh.empty() && lh.front() < t - LEARN_WIN) lh.pop_front();
  if (h.size() < 40) return false;
  // a Quest: whichever of its two periods the frames fold at, if either folds cleanly (otherwise it's learned)
  if (quest_ && !c.p && !c.learn) {
    double a = FoldSpread(h, kQuest50), b = FoldSpread(h, kQuest60), p = b < a ? kQuest60 : kQuest50;
    if (std::min(a, b) <= 0.5 * p) c.p = p;
    else c.learn = true;
  }
  // learned: from the frames so far once the grid's window is full, then refined every WIN over up to LEARN_WIN (a
  // longer stretch pins the period closer). A stretch too messy to fold (frames misread as short ones) keeps the period
  // found before, unless the camera's rate changed: the Frame's flicker mode (automatic by default) can switch between
  // 50 and 60 Hz mid-session. Then the last WIN alone no longer folds at the old period, but does at a new one, found
  // afresh, and the frames from before the switch go
  if (!quest_ || c.learn) {
    if (t >= c.next && h.back() - h.front() >= 0.9 * WIN) {
      double p = Learn(lh, c.p);
      if (p > 0) c.p = p;
      else if (c.p > 0 && FoldSpread(h, c.p) > 0.5 * c.p) {
        double q = Learn(h, 0);
        if (q > 0 && std::fabs(q - c.p) > 0.005 * c.p) {
          c.p = q;
          lh = h;
        }
      }
      c.next = t + WIN;
    }
  }
  double P = c.p;
  if (P <= 0) return false;
  double grid;
  FoldGap(h, P, &grid);
  double l = PyMod(t - grid, P);
  lag = l > P - std::min(0.01, P / 8) ? l - P : l;
  return true;
}

void Clock::AddArrival(double pc, double hs) {
  double d = pc - hs;
  long long k = (long long)std::floor(hs / 10);
  auto it = env_.find(k);
  if (it == env_.end() || d < it->second.d) { env_[k] = {hs, d, 0}; dirty_ = true; }
}

void Clock::AddPing(double pc_send, double pc_recv, double hs) {
  double rtt = pc_recv - pc_send, d = (pc_send + pc_recv) / 2 - hs;
  if (rtt < 0 || rtt > 1.0) return;
  long long k = (long long)std::floor(hs / 10);
  auto it = env_.find(k);
  if (it == env_.end() || rtt < it->second.rtt) { env_[k] = {hs, d, rtt}; dirty_ = true; }
  rtt_ = rtt_ <= 0 ? rtt : 0.9 * rtt_ + 0.1 * rtt;
}

void Clock::Fit() {
  dirty_ = false;
  while (env_.size() > 30) env_.erase(env_.begin());
  if (env_.empty()) return;
  std::vector<double> x, y;
  for (auto &kv : env_) { x.push_back(kv.second.hs); y.push_back(kv.second.d); }
  x0_ = x.back();
  if (x.size() >= 3) {
    double n = (double)x.size(), sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (size_t i = 0; i < x.size(); i++) {
      double u = x[i] - x0_;
      sx += u; sy += y[i]; sxx += u * u; sxy += u * y[i];
    }
    double den = n * sxx - sx * sx;
    a_ = den != 0 ? (n * sxy - sx * sy) / den : 0;
    b_ = (sy - a_ * sx) / n;
  } else {
    a_ = 0;
    b_ = *std::min_element(y.begin(), y.end());
  }
  has_ = true;
}

bool Clock::Map(double hs, double &pc) {
  if (dirty_) Fit();
  if (!has_) return false;
  pc = hs + a_ * (hs - x0_) + b_;
  return true;
}

// ================================================================ poses
// Faster than this from a to b, the head didn't get there: the pose stream jumped
static bool PoseJump(const PoseS &a, const PoseS &b, bool *rot = nullptr, bool *pos = nullptr) {
  double dt = b.t - a.t;
  bool r = QuatDeg(a.q, b.q) > 5.0 + 500.0 * dt, p = norm(b.p - a.p) > 0.05 + 3.0 * dt;
  if (rot) *rot = r;
  if (pos) *pos = p;
  return r || p;
}

// A jump waits kBrkSettle s for the poses after it, since streamed poses come late, early and stale. One that comes
// back near the poses just before it made it a glitch; the head's own motion either side (with kBrkSlack s of
// timestamp slack) or a stall just before it (the jump catches up with the head) explain it too. Only the rest is a
// break: the Quest may have moved its space under the older sightings.
static constexpr double kBrkSettle = 0.2, kBrkSpan = 0.12, kBrkSlack = 0.05, kBrkStall = 0.05, kBrkBefore = 0.15;

void PoseHist::Add(double t, const Quat &q, V3 p) {
  std::lock_guard<std::mutex> g(m_);
  PoseS s{t, q, p};
  if (i_) {
    const PoseS &l = ring_[(i_ - 1) % N];
    if (t - l.t <= 0) return;  // poses must move forward in time
    if (t - l.t > 1.0) {  // a stall this long is a break, whatever the poses
      pend_ = 0;
      Break(l, s);
    } else if (pend_) {
      if (Back(s)) { i_ = pend_; pend_ = 0; }  // a glitch: its poses go
    } else if (PoseJump(l, s)) {
      pend_ = i_;
    }
  }
  ring_[i_ % N] = s;
  i_++;
  if (pend_ && t - ring_[pend_ % N].t >= kBrkSettle) Settle();
}

void PoseHist::Break(const PoseS &a, const PoseS &b) {
  brk_t_ = b.t;
  brk_what_ = Fmt("%.2f s, %.0f cm, %.1f deg", b.t - a.t, norm(b.p - a.p) * 100, QuatDeg(a.q, b.q));
  brk_id_++;
}

// s back near the poses just before the pending jump, and well off where the jump went?
bool PoseHist::Back(const PoseS &s) const {
  const PoseS &a = ring_[(pend_ - 1) % N], &c = ring_[pend_ % N];
  bool jr, jp;
  PoseJump(a, c, &jr, &jp);
  double da = QuatDeg(a.q, c.q), dp = norm(c.p - a.p);
  uint64_t lo = i_ > (uint64_t)N ? i_ - N : 0;
  for (uint64_t k = pend_; k-- > lo;) {
    const PoseS &b = ring_[k % N];
    if (b.t < a.t - kBrkSlack) break;
    if (!PoseJump(b, s) && (!jr || QuatDeg(b.q, s.q) < 0.5 * da) && (!jp || norm(s.p - b.p) < 0.5 * dp)) return true;
  }
  return false;
}

void PoseHist::Settle() {
  uint64_t lo = i_ > (uint64_t)N ? i_ - N : 0, ia = pend_ - 1, ic = pend_;
  pend_ = 0;
  const PoseS &a = ring_[ia % N], &c = ring_[ic % N];
  double dt = c.t - a.t, da = QuatDeg(a.q, c.q), dp = norm(c.p - a.p);
  // the head's speed over kBrkSpan before the jump and after it
  double w = 0, v = 0;
  auto speed = [&](uint64_t i, uint64_t j) {
    const PoseS &x = ring_[i % N], &y = ring_[j % N];
    double d = std::max(1e-3, y.t - x.t);
    w = std::max(w, QuatDeg(x.q, y.q) / d);
    v = std::max(v, norm(y.p - x.p) / d);
  };
  uint64_t ib, id;
  ib = Find(a.t - kBrkSpan, ib) ? ib + 1 : lo;
  if (ib < ia) speed(ib, ia);
  else if (ia > lo) speed(ia - 1, ia);
  id = Find(c.t + kBrkSpan, id) ? std::min(id + 1, i_ - 1) : i_ - 1;
  speed(ic, std::max(id, ic + 1));
  bool moving = da <= 5.0 + w * (dt + kBrkSlack) && dp <= 0.05 + v * (dt + kBrkSlack);
  double span = 0;  // from the pose before the oldest stall within kBrkBefore
  for (uint64_t k = ic; k > lo && ring_[k % N].t >= c.t - kBrkBefore; k--) {
    const PoseS &y = ring_[k % N], &x = ring_[(k - 1) % N];
    if (y.t - x.t > kBrkStall) span = c.t - x.t;
  }
  bool stale = span > 0 && da <= 5.0 + 800.0 * span && dp <= 0.05 + 2.0 * span;
  if (!moving && !stale) Break(a, c);
}

bool PoseHist::Find(double t, uint64_t &idx) const {
  if (!i_) return false;
  uint64_t lo = i_ > (uint64_t)N ? i_ - N : 0, hi = i_;  // valid [lo, hi)
  if (ring_[lo % N].t > t) return false;
  while (hi - lo > 1) {  // newest index with time <= t
    uint64_t mid = lo + (hi - lo) / 2;
    if (ring_[mid % N].t <= t) lo = mid; else hi = mid;
  }
  idx = lo;
  return true;
}

bool PoseHist::AtQ(double t, Quat &q, V3 &p) const {
  std::lock_guard<std::mutex> g(m_);
  uint64_t i;
  if (!Find(t, i) || i == i_ - 1) return false;  // no extrapolation past the newest pose
  const PoseS &a = ring_[i % N], &b = ring_[(i + 1) % N];
  if (b.t - a.t > 0.05) return false;
  double s = (t - a.t) / (b.t - a.t);
  q = Slerp(a.q, b.q, s);
  p = a.p + (b.p - a.p) * s;
  return true;
}

bool PoseHist::At(double t, M3 &R, V3 &p) const {
  Quat q;
  if (!AtQ(t, q, p)) return false;
  R = ToM3(q);
  return true;
}

bool PoseHist::Still(double t, double span, double rot, double pos) const {
  std::lock_guard<std::mutex> g(m_);
  uint64_t hi;
  if (!Find(t, hi)) return true;
  uint64_t lo;
  if (!Find(t - span, lo)) lo = i_ > (uint64_t)N ? i_ - N : 0;
  else if (ring_[lo % N].t < t - span) lo++;
  if (hi < lo || hi - lo + 1 < 10 || ring_[lo % N].t > t - 0.9 * span) return true;  // too little history to tell
  const PoseS &l = ring_[hi % N];
  uint64_t stride = std::max<uint64_t>(1, (hi - lo + 1) / 512);
  for (uint64_t i = lo; i <= hi; i += stride) {
    const PoseS &s = ring_[i % N];
    if (QuatDeg(s.q, l.q) >= rot || norm(s.p - l.p) >= pos) return false;
  }
  return true;
}

bool PoseHist::Speed(double t, double &w, double span) const {
  Quat a, b;
  V3 pa, pb;
  if (!AtQ(t - span, a, pa) || !AtQ(t, b, pb)) return false;
  w = QuatDeg(a, b) / span;
  return true;
}

bool PoseHist::Latest(V3 &p, double *t) const {
  std::lock_guard<std::mutex> g(m_);
  if (!i_) return false;
  p = ring_[(i_ - 1) % N].p;
  if (t) *t = ring_[(i_ - 1) % N].t;
  return true;
}

// ================================================================ reference frame
static Quat QuatFromJ(const JVal *v) {
  auto q = v ? v->nums() : std::vector<double>();
  if (q.size() < 4) return {};
  return {q[0], q[1], q[2], q[3]};  // w x y z
}

bool StationsFile::Load(const std::string &path) {
  *this = StationsFile();
  std::string text;
  JVal d;
  if (!ReadFile(path, text) || !JParse(text, d)) return false;
  const JVal *a = d.get("anchor"), *ap = d.get("anchor_p"), *aq = d.get("anchor_q");
  if (!a || a->t != JVal::Str || !ap || !aq) return false;
  auto p = ap->nums();
  if (p.size() < 3) return false;
  anchor = a->s;
  anchor_p = {p[0], p[1], p[2]};
  anchor_R = ToM3(QuatFromJ(aq));
  if (const JVal *lq = d.get("level")) level = ToM3(QuatFromJ(lq));
  if (const JVal *st = d.get("stations"))
    for (auto &kv : st->o) {
      if (const JVal *pos = kv.second.get("pos")) {
        auto v = pos->nums();
        if (v.size() >= 3) fix[kv.first] = {v[0], v[1], v[2]};
      }
      if (const JVal *rp = kv.second.get("ref")) {
        auto v = rp->nums();
        if (v.size() >= 3) ref[kv.first] = {V3{v[0], v[1], v[2]}, ToM3(QuatFromJ(kv.second.get("ref_q")))};
      }
      const JVal *pl = kv.second.get("place"), *ps = kv.second.get("place_sd");
      if (pl && ps) {
        auto v = pl->nums(), s = ps->nums();
        if (v.size() >= 2 && s.size() >= 2 && s[0] > 0 && s[1] > 0) place[kv.first] = {v[0], v[1], s[0], s[1]};
      }
    }
  const JVal *ag = d.get("auto");
  autogen = ag && ag->t == JVal::Bool && ag->b;
  loaded = true;
  return true;
}

bool StationsFile::SaveAuto(const std::string &path) const {
  Quat q = ToQuat(anchor_R), l = ToQuat(level);
  std::string st;
  std::set<std::string> keys;
  for (auto &kv : ref) keys.insert(kv.first);
  for (auto &kv : place) keys.insert(kv.first);
  for (const std::string &k : keys) {
    std::string e;
    auto r = ref.find(k);
    if (r != ref.end()) {
      Quat rq = ToQuat(r->second.second);
      const V3 &p = r->second.first;
      e = Fmt("\"ref\": [%.6f, %.6f, %.6f], \"ref_q\": [%.6f, %.6f, %.6f, %.6f]", p.x, p.y, p.z, rq.w, rq.x, rq.y, rq.z);
    }
    auto pl = place.find(k);
    if (pl != place.end())
      e += Fmt("%s\"place\": [%.4f, %.4f], \"place_sd\": [%.4f, %.4f]", e.empty() ? "" : ", ", pl->second.d,
               pl->second.h, pl->second.sd_d, pl->second.sd_h);
    st += Fmt("%s\n  \"%s\": {%s}", st.empty() ? "" : ",", k.c_str(), e.c_str());
  }
  if (!st.empty()) st += "\n ";
  std::string s = Fmt(
      "{\n \"auto\": true,\n \"note\": \"QuestLHSync's reference frame: %s's pose when it was first seen. SteamVR re-tilts "
      "its lighthouse universe at every start; this pins it. stations: each base station's pose in it; with three or "
      "more, the frame follows the best fit of them all. place: with two, the other's distance and height from the "
      "anchor (m) as the headset's cameras measured them. level: the turn that levels it with the headset's gravity, "
      "found from the base stations' sightings. Delete to start over.\",\n", anchor.c_str());
  s += Fmt(" \"anchor\": \"%s\",\n \"anchor_p\": [%.6f, %.6f, %.6f],\n \"anchor_q\": [%.6f, %.6f, %.6f, %.6f],\n"
           " \"level\": [%.8f, %.8f, %.8f, %.8f],\n \"stations\": {", anchor.c_str(), anchor_p.x, anchor_p.y,
           anchor_p.z, q.w, q.x, q.y, q.z, l.w, l.x, l.y, l.z);
  s += st + Fmt("},\n \"saved\": \"%s\"\n}\n", Now().c_str());
  return WriteFileAtomic(path, s);
}

// symmetric 4 x 4: the eigenvector of the largest eigenvalue (Jacobi)
static void MaxEigVec4(double A[4][4], double v[4]) {
  double a[4][4], V[4][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
  memcpy(a, A, sizeof a);
  for (int sweep = 0; sweep < 50; sweep++) {
    double off = 0;
    for (int i = 0; i < 4; i++)
      for (int j = i + 1; j < 4; j++) off += a[i][j] * a[i][j];
    if (off < 1e-30) break;
    for (int p = 0; p < 4; p++)
      for (int q = p + 1; q < 4; q++) {
        if (std::fabs(a[p][q]) < 1e-300) continue;
        double th = (a[q][q] - a[p][p]) / (2 * a[p][q]);
        double t = (th >= 0 ? 1 : -1) / (std::fabs(th) + std::sqrt(th * th + 1));
        double c = 1 / std::sqrt(t * t + 1), s = t * c;
        for (int k = 0; k < 4; k++) {
          double akp = a[k][p], akq = a[k][q];
          a[k][p] = c * akp - s * akq;
          a[k][q] = s * akp + c * akq;
          double vkp = V[k][p], vkq = V[k][q];
          V[k][p] = c * vkp - s * vkq;
          V[k][q] = s * vkp + c * vkq;
        }
        for (int k = 0; k < 4; k++) {
          double apk = a[p][k], aqk = a[q][k];
          a[p][k] = c * apk - s * aqk;
          a[q][k] = s * apk + c * aqk;
        }
      }
  }
  int b = 0;
  for (int i = 1; i < 4; i++)
    if (a[i][i] > a[b][b]) b = i;
  for (int k = 0; k < 4; k++) v[k] = V[k][b];
}

// the turn R that puts points a best on points b about their centres: R (a - ca) ~ b - cb. Three or more off one line:
// Horn's quaternion method; fewer: a turn about the vertical alone, keeping `prev`'s tilt (one: `prev` itself).
// Returns the worst miss.
static double FitTurn(const std::vector<V3> &a, const std::vector<V3> &b, const std::vector<char> &use, const M3 &prev,
                      M3 &R, V3 &ca, V3 &cb) {
  std::vector<size_t> u;
  for (size_t i = 0; i < a.size(); i++)
    if (use[i]) u.push_back(i);
  ca = cb = {};
  for (size_t i : u) { ca += a[i]; cb += b[i]; }
  ca = ca * (1.0 / u.size());
  cb = cb * (1.0 / u.size());
  double line = 0;  // how far the stations reach off the line through the two furthest apart
  if (u.size() >= 3) {
    size_t i0 = u[0], i1 = u[1];
    for (size_t i : u)
      for (size_t j : u)
        if (norm(a[i] - a[j]) > norm(a[i0] - a[i1])) { i0 = i; i1 = j; }
    V3 d = a[i1] - a[i0];
    d = d * (1.0 / norm(d));
    for (size_t i : u) line = std::max(line, norm(cross(a[i] - a[i0], d)));
  }
  if (line >= kLine) {
    double S[3][3] = {};
    for (size_t i : u) {
      V3 p = a[i] - ca, q = b[i] - cb;
      for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++) S[r][c] += p[r] * q[c];
    }
    double N[4][4] = {
        {S[0][0] + S[1][1] + S[2][2], S[1][2] - S[2][1], S[2][0] - S[0][2], S[0][1] - S[1][0]},
        {S[1][2] - S[2][1], S[0][0] - S[1][1] - S[2][2], S[0][1] + S[1][0], S[2][0] + S[0][2]},
        {S[2][0] - S[0][2], S[0][1] + S[1][0], -S[0][0] + S[1][1] - S[2][2], S[1][2] + S[2][1]},
        {S[0][1] - S[1][0], S[2][0] + S[0][2], S[1][2] + S[2][1], -S[0][0] - S[1][1] + S[2][2]}};
    double v[4];
    MaxEigVec4(N, v);
    R = ToM3(Quat{v[0], v[1], v[2], v[3]});
  } else {
    const auto &m = prev.m;
    M3 tilt = Ry(-std::atan2(m[0][2] - m[2][0], m[0][0] + m[2][2])) * prev;
    double sn = 0, cs = 0;
    for (size_t i : u) {
      V3 p = tilt * (a[i] - ca), q = b[i] - cb;
      sn += p.z * q.x - p.x * q.z;
      cs += p.x * q.x + p.z * q.z;
    }
    R = u.size() >= 2 ? Ry(std::atan2(sn, cs)) * tilt : prev;
  }
  double worst = 0;
  for (size_t i : u) worst = std::max(worst, norm(R * (a[i] - ca) - (b[i] - cb)));
  return worst;
}

// An automatic frame of three or more stations: SteamVR's frame -> the reference frame is the best fit of all the
// stations' places, not one station's pose. SteamVR re-solves the stations now and then, each a few cm and a degree or
// so its own way; one station's pose then turns the frame, and the level with it, while the fit of them all hardly
// moves. A station kMoved off when the rest fit within it is left out (`left`, `off`: how far); `miss`: the worst miss
// of the stations kept.
bool Frame::Layout(const std::map<std::string, std::pair<V3, M3>> &raw, M3 &Cs, V3 &o, std::string &left,
                   double &off, double &miss) const {
  std::vector<std::string> keys;
  std::vector<V3> a, b;
  M3 prev = Cs_;
  for (auto &kv : raw) {
    auto r = f_.ref.find(kv.first);
    if (r == f_.ref.end()) continue;
    if (keys.empty() && ok_ != 1) prev = r->second.second * T(kv.second.second);  // no fit yet: a station's own turn
    keys.push_back(kv.first);
    a.push_back(kv.second.first);
    b.push_back(r->second.first);
  }
  if (keys.empty()) return false;
  std::vector<char> use(keys.size(), 1);
  V3 ca, cb;
  double worst = FitTurn(a, b, use, prev, Cs, ca, cb);
  left.clear();
  off = 0;
  if (keys.size() >= 3 && worst > kMoved) {  // one station moved: the fit of the rest that misses least
    double best = kMoved;
    M3 R;
    V3 pa, pb;
    for (size_t i = 0; i < keys.size(); i++) {
      use.assign(keys.size(), 1);
      use[i] = 0;
      double w = FitTurn(a, b, use, prev, R, pa, pb);
      if (w < best) {
        best = w;
        left = keys[i];
        off = norm(R * (a[i] - pa) - (b[i] - pb));
        Cs = R; ca = pa; cb = pb;
      }
    }
    if (!left.empty()) worst = best;
  }
  miss = worst;
  o = ca + T(Cs) * (f_.anchor_p - cb);  // about the anchor's place, as the level turns it
  return true;
}

void Frame::Update(const std::map<std::string, std::pair<V3, M3>> &raw) {
  if (!f_.loaded) return;
  M3 Cs;
  V3 o;
  bool lay = layout();
  if (lay) {
    std::string left;
    double off;
    if (!Layout(raw, Cs, o, left, off, miss_)) return;
    if (left != left_ && !left.empty())
      log_(Fmt("lighthouse frame: base station %s is %.0f cm off against the others: left out of the fit", left.c_str(),
               off * 100));
    left_ = left;
  } else {
    auto it = raw.find(f_.anchor);
    if (it == raw.end()) return;
    Cs = f_.anchor_R * T(it->second.second);
    o = it->second.first;
  }
  double ang = RotDeg(Cs);
  int ok = lay || ang < kMaxRot;
  if (ok != ok_ || (ok && RotDeg(Cs * T(Cs_)) > 0.01))
    log_(Fmt("lighthouse frame: SteamVR's is %.2f deg from the reference%s", ang,
             ok ? "" : " (over 3 deg: the anchor moved or a new universe; SteamVR's frame as it is)"));
  ok_ = ok;
  deg_ = ang;
  if (ok) {
    Cs_ = Cs;
    GC_ = grav_ ? G_ : M3();
    C_ = GC_ * f_.level * Cs;
    o_ = o;
    oref_ = f_.anchor_p;
  } else {
    Cs_ = C_ = GC_ = M3();
    o_ = {};
    oref_ = {};
  }
}

void Frame::Reset() {
  Cs_ = C_ = GC_ = M3();
  o_ = oref_ = V3{};
  ok_ = -1;
  deg_ = miss_ = 0;
  left_.clear();
}

// ================================================================ solver
Solver::Solver(LogFn log, uint64_t seed, std::map<std::string, V3> fix)
    : fix_(std::move(fix)), log_(std::move(log)), rng_(seed) {}

void Solver::Reset(const X4 &x) {
  x_ = x; has_x_ = true;
  resets_++;
  anchor_ = x; has_anchor_ = true;
  id_said_.clear();
}

// support for a candidate x: recent rays within gate of any predicted station direction
int Solver::ScoreX(const X4 &x, double gate) {
  std::vector<std::string> keys;
  std::vector<V3> S, Z, P, Zq;
  Stations(keys, S, Z);
  Predict(x, S, Z, P, Zq);
  Rays r = GetRays(std::max(seen_ - kKeep, since_));
  std::vector<int> cnt(S.size(), 0);
  for (size_t n = 0; n < r.size(); n++) {
    if (!r.B[n] && !starved_) continue;
    int k;
    double a;
    Nearest(P, Zq, r.O[n], r.D[n], k, a);
    if (a < gate) cnt[k]++;
  }
  // yaw needs TWO stations: score by the second-best station's count, or a fit that piles every
  // ray on one station (a lamp the alignment happens to sit on) reads as a strong solution.
  return Score(cnt);
}

// the sightings, the fit and what hung on them: gone (Sync::ForceRecalibrate). The stations stay.
void Solver::Clear() {
  std::lock_guard<std::mutex> g(m_);
  rt_.clear(); rg_.clear(); ro_.clear(); rd_.clear(); rc_.clear(); rb_.clear();
  frames_.clear();
  has_x_ = has_anchor_ = has_acq_x_ = false;
  x_ = anchor_ = acq_x_ = X4{};
  since_ = brk_ = last_acq_ = seen_ = last_id_ = last_layout_ = -1e18;
  id_said_.clear();
  stat_ = StepStat();
  has_stat_ = false;
  acq_forced_ = acq_dim_ = false;
  starved_ = false;
  amb_state_ = 0;
  kept_said_ = -1e18;
}

void Solver::Add(double t, V3 o, V3 d, double gt, int cam, bool bright) {
  std::lock_guard<std::mutex> g(m_);
  rt_.push_back(t); ro_.push_back(o); rd_.push_back(d); rg_.push_back(gt); rc_.push_back(cam); rb_.push_back(bright);
}

void Solver::Replace(double from, const std::vector<double> &T, const std::vector<V3> &O, const std::vector<V3> &D,
                     const std::vector<double> &G, const std::vector<int> &C, const std::vector<char> &B) {
  std::lock_guard<std::mutex> g(m_);
  size_t j = std::lower_bound(rt_.begin(), rt_.end(), from) - rt_.begin();
  rt_.resize(j); ro_.resize(j); rd_.resize(j); rg_.resize(j); rc_.resize(j); rb_.resize(j);
  rt_.insert(rt_.end(), T.begin(), T.end());
  ro_.insert(ro_.end(), O.begin(), O.end());
  rd_.insert(rd_.end(), D.begin(), D.end());
  rg_.insert(rg_.end(), G.begin(), G.end());
  rc_.insert(rc_.end(), C.begin(), C.end());
  rb_.insert(rb_.end(), B.begin(), B.end());
}

void Solver::SetChannel(const std::string &serial, int channel) {
  std::lock_guard<std::mutex> g(m_);
  ch_[serial] = channel;
}

void Solver::AddFrame(double t, double gt, int cam, V3 o, const M3 &R) {
  std::lock_guard<std::mutex> g(m_);
  frames_.push_back({t, gt, cam, o, R});
  while (frames_.front().t < t - kKeep) frames_.pop_front();
}

void Solver::SetStation(const std::string &serial, V3 p, const M3 &R) {
  auto f = fix_.find(serial);
  if (f != fix_.end()) {
    double d = norm(p - f->second);
    if (d < kMoved) { p = f->second; stale_.erase(serial); }
    else if (!stale_.count(serial)) {
      stale_.insert(serial);
      log_(Fmt("%s: SteamVR has it %.0f cm from where it was measured: moved? using SteamVR's pose", serial.c_str(), d * 100));
    }
  }
  std::lock_guard<std::mutex> g(m_);
  S_[serial] = {p, R};
}

Solver::Rays Solver::GetRays(double t0) {
  std::lock_guard<std::mutex> g(m_);
  if (!rt_.empty() && rt_.front() < rt_.back() - kKeep) {
    size_t j = std::lower_bound(rt_.begin(), rt_.end(), rt_.back() - kKeep) - rt_.begin();
    rt_.erase(rt_.begin(), rt_.begin() + j);
    ro_.erase(ro_.begin(), ro_.begin() + j);
    rd_.erase(rd_.begin(), rd_.begin() + j);
    rg_.erase(rg_.begin(), rg_.begin() + j);
    rc_.erase(rc_.begin(), rc_.begin() + j);
    rb_.erase(rb_.begin(), rb_.begin() + j);
  }
  Rays r;
  size_t j = std::lower_bound(rt_.begin(), rt_.end(), t0) - rt_.begin();
  r.T.assign(rt_.begin() + j, rt_.end());
  r.O.assign(ro_.begin() + j, ro_.end());
  r.D.assign(rd_.begin() + j, rd_.end());
  r.G.assign(rg_.begin() + j, rg_.end());
  r.C.assign(rc_.begin() + j, rc_.end());
  r.B.assign(rb_.begin() + j, rb_.end());
  return r;
}

// the sightings since t0: their frame-grid times, cameras and world directions (Sync names candidates with them)
void Solver::SightingsSince(double t0, std::vector<double> &G, std::vector<int> &C, std::vector<V3> &O, std::vector<V3> &D) {
  std::lock_guard<std::mutex> g(m_);
  size_t j = std::lower_bound(rt_.begin(), rt_.end(), t0) - rt_.begin();
  G.assign(rg_.begin() + j, rg_.end());
  C.assign(rc_.begin() + j, rc_.end());
  O.assign(ro_.begin() + j, ro_.end());
  D.assign(rd_.begin() + j, rd_.end());
}

// the camera frames since t0: which camera, its grid time, and its origin and rotation in the world (the in-view test)
void Solver::FramesSince(double t0, std::vector<int> &C, std::vector<double> &G, std::vector<V3> &O, std::vector<M3> &R) {
  std::lock_guard<std::mutex> g(m_);
  C.clear();
  G.clear();
  O.clear();
  R.clear();
  for (auto &f : frames_)
    if (f.t >= t0) { C.push_back(f.cam); G.push_back(f.g); O.push_back(f.o); R.push_back(f.R); }
}

void Solver::Stations(std::vector<std::string> &keys, std::vector<V3> &S, std::vector<V3> &Z) const {
  std::lock_guard<std::mutex> g(m_);
  keys.clear(); S.clear(); Z.clear();
  for (auto &kv : S_) { keys.push_back(kv.first); S.push_back(kv.second.first); Z.push_back(kv.second.second.col(2)); }
}

static bool Inv3(const M3 &a, M3 &o) {
  const double(*m)[3] = a.m;
  double c00 = m[1][1] * m[2][2] - m[1][2] * m[2][1], c01 = m[1][2] * m[2][0] - m[1][0] * m[2][2],
         c02 = m[1][0] * m[2][1] - m[1][1] * m[2][0];
  double det = m[0][0] * c00 + m[0][1] * c01 + m[0][2] * c02;
  if (!(std::fabs(det) > 1e-300)) return false;
  double k = 1 / det;
  o.m[0][0] = c00 * k; o.m[1][0] = c01 * k; o.m[2][0] = c02 * k;
  o.m[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * k;
  o.m[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * k;
  o.m[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * k;
  o.m[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * k;
  o.m[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * k;
  o.m[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * k;
  return true;
}

// Each station's place from its sightings alone (headset space): the point nearest their rays by angle (reweighted:
// a ray kTriSig off counts half). The rays are the ones nearest it within GATE of where the alignment puts it, then
// within kPick of the point itself: the alignment can be a degree off a station whose place is what's wrong, and a
// narrow cone about it would pull the point its way (a recording: 4 cm in height). The covariance takes
// kTriSig per ray, scaled up by the rays per 5 cm head cell: rays from one spot add no parallax, and a still head's
// depth is barely known (recordings: windows a head mostly sat in came out 10-20 cm off, with a covariance saying so)
bool Solver::Triangulate(std::map<std::string, Seen> &out) {
  constexpr double kTriSig = 0.25 / kDeg, kPick = 1.0 / kDeg, kIn = 0.5, kCell = 0.05;
  out.clear();
  if (!has_x_) return false;
  std::vector<std::string> keys;
  std::vector<V3> S, Z, P, Zq;
  Stations(keys, S, Z);
  Rays r = GetRays(std::max(seen_ - kKeep, since_));
  Predict(x_, S, Z, P, Zq);
  bool dim = starved_;
  std::vector<std::vector<size_t>> pick(S.size());
  for (size_t n = 0; n < r.size(); n++) {
    if (!r.B[n] && !dim) continue;
    int k;
    double a;
    Nearest(P, Zq, r.O[n], r.D[n], k, a);
    if (a < GATE) pick[k].push_back(n);
  }
  for (size_t k = 0; k < S.size(); k++) {
    if (pick[k].size() < 20) continue;
    V3 q = P[k];
    bool ok = true;
    for (int it = 0; it < 10 && ok; it++) {
      M3 A;
      A.m[0][0] = A.m[1][1] = A.m[2][2] = 0;
      double bv[3] = {0, 0, 0};
      for (size_t n : pick[k]) {
        V3 u = q - r.O[n], d = r.D[n];
        double rr = norm(u), ang = std::acos(std::max(-1.0, std::min(1.0, dot(u, d) / rr)));
        if (it >= 2 && ang > kPick) continue;
        double w = 1 / (1 + (ang / kTriSig) * (ang / kTriSig)) / (rr * rr);
        double dv[3] = {d.x, d.y, d.z}, ov[3] = {r.O[n].x, r.O[n].y, r.O[n].z};
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++) {
            double pij = (i == j ? 1 : 0) - dv[i] * dv[j];
            A.m[i][j] += w * pij;
            bv[i] += w * pij * ov[j];
          }
      }
      M3 Ai;
      if (!(ok = Inv3(A, Ai))) break;
      q = Ai * V3{bv[0], bv[1], bv[2]};
    }
    if (!ok) continue;
    M3 A;
    A.m[0][0] = A.m[1][1] = A.m[2][2] = 0;
    std::set<std::array<long long, 3>> cells;
    int nin = 0;
    double t0 = 1e300;
    for (size_t n : pick[k]) {
      V3 u = q - r.O[n], d = r.D[n];
      double rr = norm(u), ang = std::acos(std::max(-1.0, std::min(1.0, dot(u, d) / rr)));
      if (ang * kDeg >= kIn) continue;
      double w = 1 / ((kTriSig * rr) * (kTriSig * rr));
      double dv[3] = {d.x, d.y, d.z};
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) A.m[i][j] += w * ((i == j ? 1 : 0) - dv[i] * dv[j]);
      cells.insert({(long long)std::floor(r.O[n].x / kCell), (long long)std::floor(r.O[n].y / kCell),
                    (long long)std::floor(r.O[n].z / kCell)});
      nin++;
      t0 = std::min(t0, r.T[n]);
    }
    M3 C;
    if (nin < 20 || !Inv3(A, C)) continue;
    // and by how far the rays scatter about it, when more than kTriSig (lit rooms, an unlearned timing: up to twice)
    std::vector<double> angs;
    for (size_t n : pick[k]) {
      V3 u = q - r.O[n];
      double ang = std::acos(std::max(-1.0, std::min(1.0, dot(u, r.D[n]) / norm(u))));
      if (ang < kPick) angs.push_back(ang);
    }
    double sig = std::max(kTriSig, Median(angs) / 1.177);  // the median of a 2D error's size: 1.177 sigma
    double s = std::max(1.0, (double)nin / cells.size()) * (sig / kTriSig) * (sig / kTriSig);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) C.m[i][j] *= s;
    out[keys[k]] = {q, C, nin, (int)cells.size(), t0};
  }
  return !out.empty();
}

void Solver::Predict(const X4 &x, const std::vector<V3> &S, const std::vector<V3> &Z, std::vector<V3> &P, std::vector<V3> &Zq) {
  double c = std::cos(x[0]), s = std::sin(x[0]);
  V3 t{x[1], x[2], x[3]};
  P.resize(S.size()); Zq.resize(S.size());
  for (size_t k = 0; k < S.size(); k++) { P[k] = RyMul(c, s, S[k]) + t; Zq[k] = RyMul(c, s, Z[k]); }
}

// nearest station to a ray by angle (deg); inf where the head is behind every station
void Solver::Nearest(const std::vector<V3> &P, const std::vector<V3> &Zq, V3 o, V3 d, int &k, double &a) {
  k = 0; a = kInf;
  for (size_t j = 0; j < P.size(); j++) {
    V3 U = P[j] - o;
    double n = norm(U);
    double ang = kInf;
    if (dot(U, Zq[j]) >= 0.2 * n) {
      double c = dot(U, d) / n;
      ang = std::acos(c < -1 ? -1 : c > 1 ? 1 : c) * kDeg;
    }
    if (ang < a) { a = ang; k = (int)j; }
  }
}

void Solver::Support(const X4 &x, const std::vector<V3> &S, const std::vector<V3> &Z, const Rays &r, double gate,
                     std::vector<int> &cnt) const {
  std::vector<V3> P, Zq;
  Predict(x, S, Z, P, Zq);
  cnt.assign(S.size(), 0);
  for (size_t n = 0; n < r.size(); n++) {
    int k; double a;
    Nearest(P, Zq, r.O[n], r.D[n], k, a);
    if (a < gate) cnt[k]++;
  }
}

// yaw needs two stations: every station's support but the best seen one's. Only the second best's, a fit of two
// stations scored the same as one a third station's sightings agree with too (a Steam Frame recording: the sightings
// of a third station told the true fit from one 5 deg and 0.7 m off it, which two stations fit as well)
int Solver::Score(const std::vector<int> &c) {
  if (c.size() < 2) return 0;
  int best = 0, sum = 0;
  for (int n : c) { sum += n; best = std::max(best, n); }
  return sum - best;
}

struct CellKey {
  int64_t k[6];
  bool operator==(const CellKey &o) const { return !memcmp(k, o.k, sizeof k); }
};
struct CellHash {
  size_t operator()(const CellKey &c) const {
    uint64_t h = 1469598103934665603ull;
    for (int i = 0; i < 6; i++) { h ^= (uint64_t)c.k[i]; h *= 1099511628211ull; }
    return (size_t)h;
  }
};
// (head 5 cm, direction ~2 deg) cell of a sighting; np.round is round-half-even, so is nearbyint
static CellKey KeyOf(V3 o, V3 d) {
  return {{(int64_t)std::nearbyint(o.x / 0.05), (int64_t)std::nearbyint(o.y / 0.05), (int64_t)std::nearbyint(o.z / 0.05),
           (int64_t)std::nearbyint(d.x * 30), (int64_t)std::nearbyint(d.y * 30), (int64_t)std::nearbyint(d.z * 30)}};
}

void Solver::Cells(const Rays &r, const std::vector<char> &use, std::vector<int> &cnt) {
  std::unordered_map<CellKey, int, CellHash> m;
  for (size_t n = 0; n < r.size(); n++)
    if (use[n]) m[KeyOf(r.O[n], r.D[n])]++;
  cnt.assign(r.size(), 0);
  for (size_t n = 0; n < r.size(); n++)
    if (use[n]) cnt[n] = m[KeyOf(r.O[n], r.D[n])];
}

constexpr int kMaxDof = 6;  // Fit's unknowns: yaw, translation and, when levelling, a tilt

// symmetric n x n (n <= kMaxDof): smallest eigenvalue (Jacobi)
static double MinEig(double A[kMaxDof][kMaxDof], int n) {
  double a[kMaxDof][kMaxDof];
  memcpy(a, A, sizeof a);
  for (int sweep = 0; sweep < 50; sweep++) {
    double off = 0;
    for (int i = 0; i < n; i++)
      for (int j = i + 1; j < n; j++) off += a[i][j] * a[i][j];
    if (off < 1e-30) break;
    for (int p = 0; p < n; p++)
      for (int q = p + 1; q < n; q++) {
        if (std::fabs(a[p][q]) < 1e-300) continue;
        double th = (a[q][q] - a[p][p]) / (2 * a[p][q]);
        double t = (th >= 0 ? 1 : -1) / (std::fabs(th) + std::sqrt(th * th + 1));
        double c = 1 / std::sqrt(t * t + 1), s = t * c;
        for (int k = 0; k < n; k++) {
          double akp = a[k][p], akq = a[k][q];
          a[k][p] = c * akp - s * akq;
          a[k][q] = s * akp + c * akq;
        }
        for (int k = 0; k < n; k++) {
          double apk = a[p][k], aqk = a[q][k];
          a[p][k] = c * apk - s * aqk;
          a[q][k] = s * apk + c * aqk;
        }
      }
  }
  double m = a[0][0];
  for (int i = 1; i < n; i++) m = std::min(m, a[i][i]);
  return m;
}

static bool Solve(double A[kMaxDof][kMaxDof], const double *b, double *x, int n) {
  double M[kMaxDof][kMaxDof + 1];
  for (int i = 0; i < n; i++) { for (int j = 0; j < n; j++) M[i][j] = A[i][j]; M[i][n] = b[i]; }
  for (int c = 0; c < n; c++) {
    int p = c;
    for (int r = c + 1; r < n; r++) if (std::fabs(M[r][c]) > std::fabs(M[p][c])) p = r;
    if (std::fabs(M[p][c]) < 1e-300) return false;
    if (p != c) for (int j = 0; j <= n; j++) std::swap(M[p][j], M[c][j]);
    for (int r = 0; r < n; r++) {
      if (r == c) continue;
      double f = M[r][c] / M[c][c];
      for (int j = c; j <= n; j++) M[r][j] -= f * M[c][j];
    }
  }
  for (int i = 0; i < n; i++) x[i] = M[i][n] / M[i][i];
  return true;
}

// Robust fit (soft_l1 on the sightings, quadratic prior rows), sightings assigned to the nearest station under x0; yaw
// about the station centroid. A (head, direction) cell weighs one sighting however many it holds, so a still head adds
// nothing by looking longer. xa: the prior's center (default x0). cond: the sightings alone pin every DOF.
// 4 DOF; with a pivot 6: the stations also tilt about it (Tilt(y[4], y[5])), and the prior holds only the tilt.
Solver::FitR Solver::Fit(const X4 &x0, const std::vector<V3> &S, const std::vector<V3> &Z, const Rays &r, double now,
                         double gate, const X4 *xa_in, const V3 *pivot, double dim_gate, const V3 *axis,
                         const std::vector<int> *who) {
  FitR out;
  std::vector<V3> P, Zq;
  Predict(x0, S, Z, P, Zq);
  std::vector<char> use(r.size(), 0);
  std::vector<int> kk(r.size(), 0);
  int m = 0;
  for (size_t n = 0; n < r.size(); n++) {
    double a;
    // who: the station is TOLD (naming by rotor timing), so there is no gate to pass -- a sighting the
    // naming left out (-1) is dropped. Without who, the nearest station within the gate, as before.
    if (who) {
      if ((*who)[n] < 0) continue;
      kk[n] = (*who)[n];
      use[n] = 1; m++;
      continue;
    }
    Nearest(P, Zq, r.O[n], r.D[n], kk[n], a);
    if (a < (r.B.empty() || r.B[n] ? gate : std::min(gate, dim_gate))) { use[n] = 1; m++; }
  }
  out.n = m;
  if (m < 15) return out;
  std::vector<int> cnt;
  Cells(r, use, cnt);
  struct Row { int k; V3 o, d; double sw; };
  std::vector<Row> rows;
  rows.reserve(m);
  for (size_t n = 0; n < r.size(); n++)
    if (use[n]) rows.push_back({kk[n], r.O[n], r.D[n], std::sqrt(std::exp(-(now - r.T[n]) / TAU) / cnt[n])});
  const X4 xa = xa_in ? *xa_in : x0;
  V3 c;
  for (auto &s : S) c += s;
  c = c * (1.0 / S.size());
  V3 cq = Ry(x0[0]) * c + V3{x0[1], x0[2], x0[3]};
  V3 cqa = Ry(xa[0]) * c + V3{xa[1], xa[2], xa[3]};
  std::vector<V3> Sc(S.size());  // about the centroid; with a pivot, about the pivot (pc: the pivot about the centroid)
  for (size_t k = 0; k < S.size(); k++) Sc[k] = S[k] - (pivot ? *pivot : c);
  const V3 pc = pivot ? *pivot - c : V3{};
  if (!pivot) axis = nullptr;
  const int np = pivot ? (axis ? 5 : 6) : 4;
  const double f = FSCALE;
  const size_t nrob = rows.size() * 3, nr = nrob + (pivot ? (axis ? 1 : 2) : 4);

  auto res = [&](const double *y, std::vector<double> &rv) {
    rv.resize(nr);
    double cy = std::cos(y[0]), sy = std::sin(y[0]);
    V3 sh = cq + V3{y[1], y[2], y[3]};
    M3 Rt;
    if (pivot) Rt = axis ? Tilt(y[4] * axis->x, y[4] * axis->z) : Tilt(y[4], y[5]);
    for (size_t i = 0; i < rows.size(); i++) {
      const Row &w = rows[i];
      V3 U = RyMul(cy, sy, pivot ? pc + Rt * Sc[w.k] : Sc[w.k]) + sh - w.o;
      U = U * (1 / norm(U));
      V3 e = cross(U, w.d) * w.sw;
      rv[3 * i] = e.x; rv[3 * i + 1] = e.y; rv[3 * i + 2] = e.z;
    }
    if (pivot) {
      rv[nrob] = f * y[4] / SIG_LEVEL;
      if (!axis) rv[nrob + 1] = f * y[5] / SIG_LEVEL;
      return;
    }
    rv[nrob] = f * Wrap(y[0] - xa[0]) / SIG_YAW;
    rv[nrob + 1] = f * (sh.x - cqa.x) / SIG_T;
    rv[nrob + 2] = f * (sh.y - cqa.y) / SIG_T;
    rv[nrob + 3] = f * (sh.z - cqa.z) / SIG_T;
  };
  auto cost = [&](const std::vector<double> &rv) {
    double s = 0;
    for (size_t i = 0; i < nrob; i++) { double z = rv[i] / f; s += 2 * f * f * (std::sqrt(1 + z * z) - 1); }
    for (size_t i = nrob; i < nr; i++) s += rv[i] * rv[i];
    return s;
  };
  auto jac = [&](const double *y, std::vector<double> &J) {  // central differences, row-major nr x np
    J.assign(nr * np, 0);
    std::vector<double> rp, rm;
    for (int j = 0; j < np; j++) {
      double h = 1e-7, yp[kMaxDof], ym[kMaxDof];
      memcpy(yp, y, sizeof yp); memcpy(ym, y, sizeof ym);
      yp[j] += h; ym[j] -= h;
      res(yp, rp); res(ym, rm);
      for (size_t i = 0; i < nr; i++) J[i * np + j] = (rp[i] - rm[i]) / (2 * h);
    }
  };

  double y[kMaxDof] = {x0[0], 0, 0, 0, 0, 0};
  std::vector<double> rv, J, rn;
  res(y, rv);
  double cst = cost(rv), lam = 1e-3;
  for (int it = 0; it < 100; it++) {
    jac(y, J);
    double A[kMaxDof][kMaxDof] = {}, g[kMaxDof] = {};
    for (size_t i = 0; i < nr; i++) {
      double w = 1;
      if (i < nrob) { double z = rv[i] / f; w = 1 / std::sqrt(1 + z * z); }  // IRLS weight of soft_l1
      const double *Ji = &J[i * np];
      for (int a = 0; a < np; a++) {
        g[a] += w * Ji[a] * rv[i];
        for (int b = 0; b < np; b++) A[a][b] += w * Ji[a] * Ji[b];
      }
    }
    bool moved = false, done = false;
    for (int tries = 0; tries < 20; tries++) {
      double Al[kMaxDof][kMaxDof], ng[kMaxDof], d[kMaxDof];
      for (int a = 0; a < np; a++) {
        for (int b = 0; b < np; b++) Al[a][b] = A[a][b];
        Al[a][a] += lam * std::max(A[a][a], 1e-12);
        ng[a] = -g[a];
      }
      if (!Solve(Al, ng, d, np)) { lam *= 10; continue; }
      double yn[kMaxDof], step = 0;
      memcpy(yn, y, sizeof yn);
      for (int a = 0; a < np; a++) yn[a] += d[a];
      res(yn, rn);
      double cn = cost(rn);
      if (cn <= cst) {
        for (int a = 0; a < np; a++) step += std::fabs(d[a]);
        done = step < 1e-11 || cst - cn <= 1e-14 * std::max(cst, 1e-30);
        memcpy(y, yn, sizeof y);
        rv.swap(rn);
        cst = cn;
        lam = std::max(lam / 3, 1e-9);
        moved = true;
        break;
      }
      lam *= 4;
      if (lam > 1e12) break;
    }
    if (!moved || done) break;
  }
  // conditioning: scipy's loss-scaled Jacobian of the sighting rows ((1+z)^-3/2 for soft_l1), per-DOF scale / SIG_R
  jac(y, J);
  double H[kMaxDof][kMaxDof] = {};
  const double sc[kMaxDof] = {COND_YAW, COND_T, COND_T, COND_T, COND_LEVEL, COND_LEVEL};
  for (size_t i = 0; i < nrob; i++) {
    double z = rv[i] / f, w = std::pow(1 + z * z, -1.5);
    const double *Ji = &J[i * np];
    for (int a = 0; a < np; a++)
      for (int b = 0; b < np; b++) H[a][b] += w * Ji[a] * Ji[b] * sc[a] * sc[b] / (SIG_R * SIG_R);
  }
  out.cond = MinEig(H, np) >= 1.0;
  if (pivot) out.tilt = axis ? Tilt(y[4] * axis->x, y[4] * axis->z) : Tilt(y[4], y[5]);
  if (axis) out.theta = y[4];
  double yaw = Wrap(y[0]);
  V3 t = cq + V3{y[1], y[2], y[3]} - Ry(yaw) * c;
  out.x = {yaw, t.x, t.y, t.z};
  out.ok = true;
  return out;
}

// acquisition: 2-ray minimal solver (one ray per station, both through their stations) + consensus
void Solver::Hypotheses(const Rays &r, const std::vector<V3> &S, std::vector<X4> &H, int M) {
  H.clear();
  size_t N = r.size(), K = S.size();
  if (N < 2 || K < 2) return;
  std::uniform_int_distribution<size_t> U(0, N - 1);
  const double c10 = std::cos(10 / kDeg);
  std::vector<std::pair<size_t, size_t>> pr;
  for (int t = 0; t < M; t++) {
    size_t i = U(rng_), j = U(rng_);
    if (dot(r.D[i], r.D[j]) < c10) pr.push_back({i, j});
  }
  for (size_t ka = 0; ka < K; ka++)
    for (size_t kb = 0; kb < K; kb++) {
      if (ka == kb) continue;
      V3 Dv = S[ka] - S[kb];  // R Dv = (oi + l di) - (oj + m dj)
      for (auto &ij : pr) {
        V3 oi = r.O[ij.first], di = r.D[ij.first], oj = r.O[ij.second], dj = r.D[ij.second];
        V3 e = oi - oj;
        double a = di.y, b = -dj.y, c = Dv.y - e.y, n2 = a * a + b * b;
        if (!(n2 > 1e-6)) continue;
        double l0 = c * a / n2, m0 = c * b / n2, dl = -b, dm = a;
        double Ax = e.x + l0 * di.x - m0 * dj.x, Az = e.z + l0 * di.z - m0 * dj.z;
        double Bx = dl * di.x - dm * dj.x, Bz = dl * di.z - dm * dj.z;
        double qa = Bx * Bx + Bz * Bz, qb = 2 * (Ax * Bx + Az * Bz), qc = Ax * Ax + Az * Az - Dv.x * Dv.x - Dv.z * Dv.z;
        double disc = qb * qb - 4 * qa * qc;
        if (!(disc >= 0) || !(qa > 1e-9)) continue;
        double sq = std::sqrt(disc);
        for (int sg = -1; sg <= 1; sg += 2) {
          double s = (-qb + sg * sq) / (2 * qa);
          double lam = l0 + s * dl, mu = m0 + s * dm;
          if (!(lam > 0.3 && mu > 0.3 && lam < 15 && mu < 15)) continue;
          V3 v = e + di * lam - dj * mu;
          double yaw = std::atan2(v.x, v.z) - std::atan2(Dv.x, Dv.z);
          double cy = std::cos(yaw), sy = std::sin(yaw);
          V3 t = oi + di * lam - RyMul(cy, sy, S[ka]);
          H.push_back({yaw, t.x, t.y, t.z});
        }
      }
    }
}

// (H, K) support counts; each ray counts for the station it's most in front of / closest to
void Solver::BatchSupport(const std::vector<X4> &H, const std::vector<V3> &S, const std::vector<V3> &Z, const Rays &r,
                          const std::vector<int> &idx, double gate, std::vector<int> &out) const {
  size_t K = S.size();
  out.assign(H.size() * K, 0);
  double cg = std::cos(gate / kDeg);
  std::vector<V3> P(K), Zq(K);
  std::vector<double> cs(K);
  std::vector<char> fr(K);
  for (size_t h = 0; h < H.size(); h++) {
    double c = std::cos(H[h][0]), s = std::sin(H[h][0]);
    V3 t{H[h][1], H[h][2], H[h][3]};
    for (size_t k = 0; k < K; k++) { P[k] = RyMul(c, s, S[k]) + t; Zq[k] = RyMul(c, s, Z[k]); }
    int *o = &out[h * K];
    for (int n : idx) {
      int best = 0;
      double bv = -kInf;
      for (size_t k = 0; k < K; k++) {
        V3 U = P[k] - r.O[n];
        double nn = norm(U);
        cs[k] = dot(U, r.D[n]) / nn;
        fr[k] = dot(U, Zq[k]) > 0.2 * nn;
        double v = fr[k] ? cs[k] : -2;
        if (v > bv) { bv = v; best = (int)k; }
      }
      if (cs[best] > cg && fr[best]) o[best]++;
    }
  }
}

bool Solver::Starved(const Rays &r) {
  size_t nb = 0;
  for (char b : r.B) nb += b != 0;
  return !r.B.empty() && r.size() >= 200 && nb * 20 < r.size();
}

Solver::Rays Solver::Bright(const Rays &r) {
  size_t nb = 0;
  for (char b : r.B) nb += b != 0;
  if (r.B.empty() || nb == r.size() || Starved(r)) return r;
  Rays o;
  for (size_t n = 0; n < r.size(); n++) {
    if (!r.B[n]) continue;
    o.T.push_back(r.T[n]); o.G.push_back(r.G[n]); o.O.push_back(r.O[n]); o.D.push_back(r.D[n]);
    o.C.push_back(r.C[n]); o.B.push_back(1);
  }
  return o;
}

// one ray per (5 cm head cell, ~2 deg world direction): a lamp seen from a still head collapses to a few rays,
// so consensus only counts parallax
Solver::Rays Solver::Thin(const Rays &r) {
  std::unordered_map<CellKey, size_t, CellHash> first;
  std::vector<size_t> keep;
  for (size_t n = 0; n < r.size(); n++)
    if (first.emplace(KeyOf(r.O[n], r.D[n]), n).second) keep.push_back(n);
  Rays o;
  for (size_t n : keep) {
    o.T.push_back(r.T[n]); o.O.push_back(r.O[n]); o.D.push_back(r.D[n]);
    if (!r.B.empty()) o.B.push_back(r.B[n]);
  }
  return o;
}

static int SecondBest(const int *c, size_t K) {  // Solver::Score's: all but the best
  if (K < 2) return 0;
  int best = 0, sum = 0;
  for (size_t k = 0; k < K; k++) { sum += c[k]; best = std::max(best, c[k]); }
  return sum - best;
}

// best hypothesis from scratch on the last ACQ_WIN s: score = thinned support, tight = the same within TIGHT_DEG
bool Solver::Acquire(double now, X4 &best, int &bs, int &tight) {
  std::vector<std::string> keys;
  std::vector<V3> S, Z;
  Stations(keys, S, Z);
  Rays all = GetRays(std::max(now - ACQ_WIN, since_)), r = Bright(all);
  acq_dim_ = r.size() == all.size() && std::count(all.B.begin(), all.B.end(), 0) > 0;
  if (S.size() < 2 || r.size() < 30) return false;
  Rays t = Thin(r);
  std::vector<X4> H;
  Hypotheses(t, S, H);
  std::vector<X4> cand;
  if (has_acq_x_) cand.push_back(acq_x_);
  size_t K = S.size();
  if (!H.empty()) {
    std::vector<int> sub(t.size());
    for (size_t i = 0; i < sub.size(); i++) sub[i] = (int)i;
    size_t ns = std::min<size_t>(t.size(), 400);
    for (size_t i = 0; i < ns; i++) {  // partial Fisher-Yates: a random subset without replacement
      std::uniform_int_distribution<size_t> u(i, sub.size() - 1);
      std::swap(sub[i], sub[u(rng_)]);
    }
    sub.resize(ns);
    std::vector<int> c1;
    BatchSupport(H, S, Z, t, sub, 1.5, c1);
    std::vector<size_t> ord(H.size());
    for (size_t i = 0; i < ord.size(); i++) ord[i] = i;
    size_t n1 = std::min<size_t>(40, ord.size());
    std::partial_sort(ord.begin(), ord.begin() + n1, ord.end(),
                      [&](size_t a, size_t b) { return SecondBest(&c1[a * K], K) > SecondBest(&c1[b * K], K); });
    std::vector<X4> top;
    for (size_t i = 0; i < n1; i++) top.push_back(H[ord[i]]);
    std::vector<int> all(t.size()), c2;
    for (size_t i = 0; i < all.size(); i++) all[i] = (int)i;
    BatchSupport(top, S, Z, t, all, 1.5, c2);
    std::vector<size_t> o2(top.size());
    for (size_t i = 0; i < o2.size(); i++) o2[i] = i;
    size_t n2 = std::min<size_t>(5, o2.size());
    std::partial_sort(o2.begin(), o2.begin() + n2, o2.end(),
                      [&](size_t a, size_t b) { return SecondBest(&c2[a * K], K) > SecondBest(&c2[b * K], K); });
    for (size_t i = 0; i < n2; i++) cand.push_back(top[o2[i]]);
  }
  bs = -1;
  bool have = false;
  for (X4 x : cand) {
    for (double gate : {1.5, 0.6}) {
      FitR f = Fit(x, S, Z, r, now, gate, nullptr);
      if (!f.ok) break;
      x = f.x;
    }
    std::vector<int> c;
    Support(x, S, Z, t, INLIER, c);
    int sc = Score(c);
    if (sc > bs) { best = x; bs = sc; have = true; }
  }
  if (!have) return false;
  CheckMirror(best, bs, S, Z, r, t, now);
  acq_x_ = best; has_acq_x_ = true;
  std::vector<int> c;
  Support(best, S, Z, t, TIGHT_DEG, c);
  tight = Score(c);
  return true;
}

int Solver::ThinnedScore(const X4 &x, double now) {
  std::vector<std::string> keys;
  std::vector<V3> S, Z;
  Stations(keys, S, Z);
  Rays t = Thin(Bright(GetRays(std::max(now - ACQ_WIN, since_))));
  if (!t.size()) return 0;
  std::vector<int> c;
  Support(x, S, Z, t, INLIER, c);
  return Score(c);
}

// x turned 180 deg about the vertical through the midpoint of stations a and b (reference frame): their places swap
static X4 Mirror(const X4 &x, V3 a, V3 b) {
  V3 m = (a + b) * 0.5, d = RyMul(std::cos(x[0]), std::sin(x[0]), V3{2 * m.x, 0, 2 * m.z});
  return {Wrap(x[0] + kPi), x[1] + d.x, x[2] + d.y, x[3] + d.z};
}

// how far apart two alignments put the stations (m, the farthest one)
static double Apart(const X4 &a, const X4 &b, const std::vector<V3> &S) {
  double ca = std::cos(a[0]), sa = std::sin(a[0]), cb = std::cos(b[0]), sb = std::sin(b[0]), far = 0;
  for (const V3 &s : S)
    far = std::max(far, norm(RyMul(ca, sa, s) + V3{a[1], a[2], a[3]} - RyMul(cb, sb, s) - V3{b[1], b[2], b[3]}));
  return far;
}

static std::string Meters(double d) { return std::isfinite(d) ? Fmt("%.1f m", d) : std::string("none"); }

// Their motion first: worn or held, they go where the head goes, and a fit and its mirror (turned 180 deg) move them
// opposite ways, wherever one stands and whatever lies about. Then how near the head they stay.
int Solver::BodyFavours(const X4 &a, const X4 &b, double &da, double &db, std::string &how) const {
  how.clear();
  da = body_ ? body_(a) : NAN;
  db = body_ ? body_(b) : NAN;
  double moved = 0, ma = motion_ ? motion_(a, moved) : NAN, mb = motion_ ? motion_(b, moved) : NAN;
  int by_motion = 0, by_place = 0;
  std::string how_motion, how_place;
  if (std::isfinite(ma) && std::isfinite(mb) && moved >= kBodyMoved && std::max(ma, mb) >= kBodyCorr &&
      std::fabs(ma - mb) >= 2 * kBodyCorr) {
    by_motion = ma > mb ? 1 : 2;
    how_motion = Fmt("the lighthouse devices worn or held move with the head one way (%+.2f against %+.2f)",
                     std::max(ma, mb), std::min(ma, mb));
  }
  if (std::isfinite(da) && std::isfinite(db) && std::min(da, db) < kBodyNear && std::fabs(da - db) > kBodyGap) {
    by_place = da < db ? 1 : 2;
    how_place = Fmt("the lighthouse devices worn or held stay near the head one way (%s from it, %s the other way)",
                    Meters(std::min(da, db)).c_str(), Meters(std::max(da, db)).c_str());
  }
  // the two kinds of evidence pointing opposite ways say neither is to be trusted (devices carried about without the
  // head, or lying in the room: their motion correlates by chance, +0.23 against -0.24 flipped a well supported fit to
  // its mirror, seven metres away, for a minute): the current fit stays
  if (by_motion && by_place && by_motion != by_place) {
    how = Fmt("the lighthouse devices worn or held disagree with themselves (their motion %+.2f against %+.2f, their "
              "distance from the head %s against %s)", ma, mb, Meters(da).c_str(), Meters(db).c_str());
    return 0;
  }
  how = by_motion ? how_motion : how_place;
  return by_motion ? by_motion : by_place;
}

// Two stations facing each other at about the same height look the same with their places swapped: the fit turned
// 180 deg about their midpoint explains the sightings almost as well, and puts every other lighthouse device across the
// room. So the mirror of the best fit is fitted too. When the two score about the same, the lighthouse devices worn or
// held (they stay near the head) decide; without them the current fit stays.
void Solver::CheckMirror(X4 &best, int &bs, const std::vector<V3> &S, const std::vector<V3> &Z, const Rays &r,
                         const Rays &t, double now) {
  acq_forced_ = false;
  std::vector<int> c;
  Support(best, S, Z, t, INLIER, c);
  size_t a = 0, b = 1;  // the two best-seen stations
  if (c[b] > c[a]) std::swap(a, b);
  for (size_t k = 2; k < c.size(); k++) {
    if (c[k] > c[a]) { b = a; a = k; }
    else if (c[k] > c[b]) b = k;
  }
  // the stations' own spacing and heights differ a little from where the cameras see their lasers: the exact mirror
  // starts a few degrees off, so its refit starts wide
  X4 m = Mirror(best, S[a], S[b]);
  for (double gate : {2.0, 1.0, 0.6}) {
    FitR f = Fit(m, S, Z, r, now, gate, nullptr);
    if (!f.ok) break;
    m = f.x;
  }
  if (Apart(m, best, S) < 0.3) return;  // the mirror slid back into the best fit
  std::vector<int> cm;
  Support(m, S, Z, t, INLIER, cm);
  int sm = Score(cm);
  // thinned support is a small, noisy count: within a factor 2 the two are a tie
  if (sm * 2 < bs) { amb_state_ = 0; return; }  // clearly worse: not ambiguous
  if (bs * 2 < sm) { amb_state_ = 0; best = m; bs = sm; return; }
  double db, dm;
  std::string how;
  int fav = 0;
  bool pick_m = false;
  int state;
  std::string why;
  // geometry prior, before wearable evidence: a station cannot stand inside
  // the volume the head actually moved through. The mirror slides one onto
  // the walked path (the head positions are fixed data the 180-deg rotation
  // does not follow), and the user's own body would have hit it -- so a
  // supported station within arm's reach of a real head position loses,
  // however the worn devices read.
  auto intrusion = [&](const X4 &x, const std::vector<int> &cs) {
    std::vector<V3> Pq, Zq;
    Predict(x, S, Z, Pq, Zq);
    double worst = 1e9;
    for (size_t k = 0; k < Pq.size() && k < cs.size(); k++) {
      if (cs[k] < 8) continue;   // only stations this candidate actually supports
      for (size_t n = 0; n < r.size(); n++) {
        double dz = Pq[k].z - r.O[n].z;
        if (std::fabs(dz) > 0.8) continue;   // well above or below the head: walking under it is fine
        double hd = std::hypot(Pq[k].x - r.O[n].x, Pq[k].y - r.O[n].y);
        if (hd < worst) worst = hd;
      }
    }
    return worst;   // min horizontal clearance between a supported station and any head position
  };
  double clear_b = intrusion(best, c), clear_m = intrusion(m, cm);
  bool hit_b = clear_b < 0.5, hit_m = clear_m < 0.5;
  if (hit_b != hit_m) {
    pick_m = hit_b;   // the candidate without a station in the walked volume wins
    state = pick_m ? 7 : 8;   // 7/8: the walked-volume prior; 1-6 are the witness's, below
    db = clear_b;
    dm = clear_m;
    why = Fmt("a station %.1f m from walked head positions on the %s, %.1f m on the other: stations do not stand inside the play volume",
              pick_m ? clear_m : clear_b, pick_m ? "mirror" : "fit", pick_m ? clear_b : clear_m);
  } else {
    fav = BodyFavours(best, m, db, dm, how);
    if (fav) {
      pick_m = fav == 2;
      state = pick_m ? 1 : 2;
      why = how + (pick_m ? ": the mirror" : ": the fit");
    } else if (has_x_) {
      pick_m = Apart(m, x_, S) < Apart(best, x_, S);
      state = how.empty() ? 3 : 5;
      why = (how.empty() ? std::string("no lighthouse device near the head tells them apart") : how) +
            ": the current fit stays";
    } else {
      state = how.empty() ? 4 : 6;
      why = (how.empty() ? std::string("no lighthouse device near the head tells them apart") : how) +
            ": wear or hold one";
    }
  }
  if (state != amb_state_)
    log_(Fmt("mirror check: yaw %+.1f (support %d, devices %s from the head) and its mirror yaw %+.1f (support %d, %s) "
             "fit about as well: %s", best[0] * kDeg, bs, Meters(db).c_str(), m[0] * kDeg, sm, Meters(dm).c_str(),
             why.c_str()));
  amb_state_ = state;
  X4 other = pick_m ? best : m;
  if (pick_m) { best = m; bs = sm; }
  // the devices' evidence flips the current fit to its own mirror; a pair without the current fit in it is just
  // another candidate
  acq_forced_ = (state == 1 || state == 2) && has_x_ && Apart(best, x_, S) > 0.3 && Apart(other, x_, S) < 1.0;
}

// the HMD pose stream broke: the Quest's space may have changed under the older sightings (boundary reset)
void Solver::PoseBreak(double t) {
  brk_ = since_ = std::max(since_, t);
  last_acq_ = -1e18;
  has_acq_x_ = false;
}

// ================================================================ station identity
// A 2.0 base station's rotor period is set by its channel (ticks of 48 MHz). A short frame sees a station only in
// part of its rotor turn, so which frames see it tells the stations apart where the geometry can't.
static const int kRotorTicks[16] = {959000, 957000, 953000, 949000, 947000, 943000, 941000, 939000,
                                    937000, 929000, 919000, 911000, 907000, 901000, 893000, 887000};

double RotorPeriod(int channel) { return channel >= 1 && channel <= 16 ? kRotorTicks[channel - 1] / 48e6 : 0; }

struct InView { int cam; double g; bool seen; };  // a short frame a station was in view of

// turns/s the phase at period tau moves from one short frame to the next
static double Drift(double tau, double P) {
  double f = std::fmod(P / tau, 1.0);
  return std::min(f, 1 - f) / P;
}

// How much seeing a station (in the frames it was in view of, sorted by camera, time) follows the phase at each
// period, over fake periods near it. Both on the same per-camera windows of whole turns of the slower drift. False if
// too few frames.
// one score per rotor period in tau: how much better the sightings follow it than a period near it. The two-period
// form was the identity check's; CheckLayout needs one per station (upstream 25486ca)
static bool SeenScores(const std::vector<InView> &s, const std::vector<double> &tau, double P, std::vector<double> &sc) {
  const double kTurns = 2, kEps[] = {0.0015, 0.0025, 0.0035, 0.0045, 0.0055, 0.0065, 0.0075, 0.0085};
  const int kMinN = 60, kMinSeen = 20;
  sc.assign(tau.size(), 0.0);
  double rate = kInf;
  for (double t : tau) rate = std::min(rate, Drift(t, P));
  if (tau.empty() || !(rate > 0)) return false;
  double W = std::max(20.0, (kTurns + 1) / rate);
  struct Win { size_t a, b; double p; };  // frames [a, b), share seen
  std::vector<Win> win;
  int n = 0, ns = 0;
  for (size_t i = 0; i < s.size();) {
    size_t j = i;
    while (j < s.size() && s[j].cam == s[i].cam && s[j].g < s[i].g + W) j++;
    double turns = std::floor((s[j - 1].g - s[i].g) * rate);
    if (turns >= kTurns) {
      size_t k = i;
      int seen = 0;
      for (; k < j && s[k].g < s[i].g + turns / rate; k++) seen += s[k].seen;
      if (seen > 0 && seen < (int)(k - i)) {
        win.push_back({i, k, seen / (double)(k - i)});
        n += (int)(k - i);
        ns += seen;
      }
    }
    i = j;
  }
  if (n < kMinN || ns < kMinSeen) return false;
  auto R = [&](double T) {
    double tot = 0;
    for (auto &w : win) {
      double c = 0, sn = 0, g0 = s[w.a].g;
      for (size_t k = w.a; k < w.b; k++) {
        double a = 2 * kPi * (s[k].g - g0) / T, e = s[k].seen - w.p;
        c += e * std::cos(a);
        sn += e * std::sin(a);
      }
      tot += std::hypot(c, sn);
    }
    return tot / n;
  };
  auto score = [&](double tau) {
    double null = 0;
    for (double e : kEps) null += R(tau * (1 + e)) + R(tau * (1 - e));
    return std::max(0.0, R(tau) - null / (2 * (sizeof kEps / sizeof kEps[0])));
  };
  for (size_t i = 0; i < tau.size(); i++) sc[i] = score(tau[i]);
  return true;
}

// yaw + translation taking points S onto Q, least squares
static X4 YawFit(const std::vector<V3> &S, const std::vector<V3> &Q) {
  V3 ms, mq;
  for (size_t i = 0; i < S.size(); i++) { ms += S[i]; mq += Q[i]; }
  ms = ms * (1.0 / S.size()); mq = mq * (1.0 / S.size());
  double sn = 0, cs = 0;
  for (size_t i = 0; i < S.size(); i++) {
    V3 a = S[i] - ms, b = Q[i] - mq;
    sn += a.z * b.x - a.x * b.z;
    cs += a.z * b.z + a.x * b.x;
  }
  double yaw = std::atan2(sn, cs);
  V3 t = mq - Ry(yaw) * ms;
  return {yaw, t.x, t.y, t.z};
}

// A pair of stations whose sightings follow each other's rotor periods better than their own is swapped: switch to
// the alignment with the two exchanged if it fits the rays about as well.
void Solver::CheckIdentity(double now, const std::vector<std::string> &keys, const std::vector<V3> &S,
                           const std::vector<V3> &Z) {
  double P = frame_p_;
  if (!(P > 0) || keys.size() < 2) return;
  std::vector<double> tau(keys.size(), 0.0);
  {
    std::lock_guard<std::mutex> g(m_);
    for (size_t k = 0; k < keys.size(); k++) {
      auto it = ch_.find(keys[k]);
      if (it != ch_.end()) tau[k] = RotorPeriod(it->second);
    }
  }
  if (!in_image_) return;
  Rays r = Bright(GetRays(since_));
  std::vector<CamFrame> fr;
  {
    std::lock_guard<std::mutex> g(m_);
    for (auto &f : frames_)
      if (f.t >= since_) fr.push_back(f);
  }
  std::vector<V3> Pq, Zq;
  Predict(x_, S, Z, Pq, Zq);
  std::vector<std::set<std::pair<int, double>>> saw(keys.size());  // (camera, frame) with a sighting within 1 deg
  for (size_t n = 0; n < r.size(); n++) {
    if (!std::isfinite(r.G[n])) continue;
    int k;
    double a;
    Nearest(Pq, Zq, r.O[n], r.D[n], k, a);
    if (a < 1.0) saw[k].insert({r.C[n], r.G[n]});
  }
  std::vector<std::vector<InView>> in(keys.size());
  for (auto &f : fr)
    for (size_t k = 0; k < keys.size(); k++) {
      V3 U = Pq[k] - f.o;
      double d = norm(U);
      if (!(dot(U, Zq[k]) >= 0.2 * d) || !in_image_(f.cam, T(f.R) * (U * (1 / d)))) continue;
      in[k].push_back({f.cam, f.g, saw[k].count({f.cam, f.g}) > 0});
    }
  for (auto &v : in)
    std::sort(v.begin(), v.end(), [](const InView &a, const InView &b) { return a.cam != b.cam ? a.cam < b.cam : a.g < b.g; });
  double worst = 0, best = 0;
  size_t wa = 0, wb = 0, ba = 0, bb = 0;
  for (size_t a = 0; a < keys.size(); a++)
    for (size_t b = a + 1; b < keys.size(); b++) {
      if (!tau[a] || !tau[b] || tau[a] == tau[b]) continue;
      double aa, ab, bb_, ba_;
      std::vector<double> sc;
      if (!SeenScores(in[a], {tau[a], tau[b]}, P, sc) || sc.size() < 2) continue;
      aa = sc[0]; ab = sc[1];
      if (!SeenScores(in[b], {tau[b], tau[a]}, P, sc) || sc.size() < 2) continue;
      bb_ = sc[0]; ba_ = sc[1];
      double gain = ab + ba_ - aa - bb_;
      if (gain > worst) { worst = gain; wa = a; wb = b; }
      if (gain < best) { best = gain; ba = a; bb = b; }
    }
  auto say = [&](const std::string &key, const std::string &s) {  // once per alignment
    if (key == id_said_) return;
    id_said_ = key;
    log_(s);
  };
  if (worst > kSwapGain) {
    std::vector<V3> Q = Pq;
    std::swap(Q[wa], Q[wb]);
    X4 xs = YawFit(S, Q);
    for (double gate : {1.5, 0.6}) {
      FitR f = Fit(xs, S, Z, r, now, gate, nullptr);
      if (!f.ok) break;
      xs = f.x;
    }
    int cur = ThinnedScore(x_, now), sw = ThinnedScore(xs, now);
    std::string pair = keys[wa] + " and " + keys[wb];
    // the devices worn or held are the stronger evidence, and the mirror check goes by them: a switch against them
    // would only be switched back
    double dc, ds;
    std::string how;
    if (BodyFavours(x_, xs, dc, ds, how) == 1) {
      say("body " + pair, Fmt("channel check: %s look the wrong way round (by %.2f), but %s: kept as it is",
                              pair.c_str(), worst, how.c_str()));
    } else if (cur > 0 && sw >= kSwapFit * cur) {
      log_(Fmt("channel check: %s were the wrong way round (by %.2f): yaw %+.2f deg t [%.3f %.3f %.3f] support %d (was %d)",
               pair.c_str(), worst, xs[0] * kDeg, xs[1], xs[2], xs[3], sw, cur));
      Reset(xs);
      has_acq_x_ = false;
    } else {
      say("kept " + pair, Fmt("channel check: %s look the wrong way round (by %.2f), but the rays fit that way worse "
                              "(support %d vs %d): kept as it is", pair.c_str(), worst, sw, cur));
    }
  } else if (best < -kSwapGain) {
    say("ok " + keys[ba] + keys[bb], Fmt("channel check: %s are the right way round (by %.2f)",
                                         (keys[ba] + " and " + keys[bb]).c_str(), -best));
  }
}

// Three or more stations: a station SteamVR's map puts beside the others and the cameras see kMoved or more elsewhere
// is logged. SteamVR's base station map is what the alignment maps onto, so a station placed wrong there poisons the
// fit; naming them by rotor timing (CheckIdentity) does not help when the map itself is off. Needs three stations with
// a unique channel each, so it does nothing on a two-station setup. Upstream 25486ca, ported for parity.
void Solver::CheckLayout(const std::vector<std::string> &keys, const std::vector<V3> &S, const std::vector<V3> &Z) {
  constexpr double kVox = 0.10, kNear = 1.0, kPar = 10.0;
  constexpr int kView = 200, kLight = 60;
  double P = frame_p_;
  if (!has_x_ || !(P > 0) || keys.size() < 3 || !in_image_) return;
  std::vector<double> tau(keys.size(), 0.0);
  {
    std::lock_guard<std::mutex> g(m_);
    for (size_t k = 0; k < keys.size(); k++) {
      auto it = ch_.find(keys[k]);
      if (it != ch_.end()) tau[k] = RotorPeriod(it->second);
    }
  }
  std::vector<size_t> u;  // the stations with a channel of their own
  std::vector<double> tu;
  for (size_t k = 0; k < keys.size(); k++)
    if (tau[k] && std::count(tau.begin(), tau.end(), tau[k]) == 1) { u.push_back(k); tu.push_back(tau[k]); }
  if (u.size() < 3) return;
  Rays r = Bright(GetRays(since_));
  std::vector<CamFrame> fr;
  {
    std::lock_guard<std::mutex> g(m_);
    for (auto &f : frames_)
      if (f.t >= since_) fr.push_back(f);
  }
  std::vector<V3> Pq, Zq;
  Predict(x_, S, Z, Pq, Zq);
  std::vector<int> at(r.size(), -1);  // the station whose place a sighting is within kNear deg of; -2: none within GATE
  for (size_t n = 0; n < r.size(); n++) {
    if (!std::isfinite(r.G[n])) continue;
    int k;
    double a;
    Nearest(Pq, Zq, r.O[n], r.D[n], k, a);
    if (a < kNear) at[n] = k;
    else if (!(a < GATE)) at[n] = -2;
  }
  // the short frames a place was in view of (facing z, if given), and which of them saw it
  auto view = [&](V3 q, const V3 *z, const std::set<std::pair<int, double>> &saw) {
    std::vector<InView> v;
    for (auto &f : fr) {
      V3 U = q - f.o;
      double d = norm(U);
      if ((z && !(dot(U, *z) >= 0.2 * d)) || !in_image_(f.cam, T(f.R) * (U * (1 / d)))) continue;
      v.push_back({f.cam, f.g, saw.count({f.cam, f.g}) > 0});
    }
    std::sort(v.begin(), v.end(), [](const InView &a, const InView &b) { return a.cam != b.cam ? a.cam < b.cam : a.g < b.g; });
    return v;
  };
  auto named = [&](const std::vector<double> &sc) {  // the rotor (of u) followed by kNamed over the others, or -1
    int b = (int)(std::max_element(sc.begin(), sc.end()) - sc.begin());
    double other = 0;
    for (size_t i = 0; i < sc.size(); i++)
      if ((int)i != b) other = std::max(other, sc[i]);
    return sc[b] - other >= kNamed ? b : -1;
  };
  std::vector<std::set<std::pair<int, double>>> saw(keys.size());
  for (size_t n = 0; n < r.size(); n++)
    if (at[n] >= 0) saw[at[n]].insert({r.C[n], r.G[n]});
  std::vector<size_t> ok, gone;  // of u: seen at their places, following their own rotors; in view there, unseen
  for (size_t i = 0; i < u.size(); i++) {
    std::vector<InView> v = view(Pq[u[i]], &Zq[u[i]], saw[u[i]]);
    int seen = 0;
    for (auto &e : v) seen += e.seen;
    std::vector<double> sc;
    if ((int)v.size() >= kView && seen * 50 <= (int)v.size()) gone.push_back(i);
    else if (SeenScores(v, tu, P, sc) && named(sc) == (int)i) ok.push_back(i);
  }
  if (ok.size() < 2 || gone.empty()) return;
  std::vector<size_t> fn, vn;  // the sightings no place explains; one of them per (head, direction) cell votes
  {
    std::unordered_map<CellKey, int, CellHash> cells;
    for (size_t n = 0; n < r.size(); n++)
      if (at[n] == -2) {
        fn.push_back(n);
        if (cells.emplace(KeyOf(r.O[n], r.D[n]), 0).second) vn.push_back(n);
      }
  }
  if ((int)fn.size() < kLight) return;
  auto off = [&](V3 q, size_t n) {  // deg
    V3 w = q - r.O[n];
    double d = norm(w);
    return d > 0.4 ? std::acos(std::max(-1.0, std::min(1.0, dot(w, r.D[n]) / d))) * kDeg : kInf;
  };
  const int64_t kO = 1 << 20;  // cells: 21 bits a coordinate
  std::vector<char> left(r.size(), 0);
  for (size_t n : fn) left[n] = 1;
  for (int it = 0; it < 6; it++) {
    std::unordered_map<int64_t, int> vote;
    for (size_t n : vn) {
      if (!left[n]) continue;
      int64_t last = -1;
      for (double d = 0.6; d < 7.0; d += kVox / 2) {
        V3 p = r.O[n] + r.D[n] * d;
        int64_t c = ((int64_t)std::floor(p.x / kVox) + kO) << 42 | ((int64_t)std::floor(p.y / kVox) + kO) << 21 |
                    ((int64_t)std::floor(p.z / kVox) + kO);
        if (c != last) vote[c]++;
        last = c;
      }
    }
    int64_t best = -1;
    int bv = 0;
    for (auto &kv : vote)
      if (kv.second > bv || (kv.second == bv && kv.first < best)) { best = kv.first; bv = kv.second; }
    if (bv < 10) break;
    V3 q{((best >> 42) - kO + 0.5) * kVox, (((best >> 21) & 0x1FFFFF) - kO + 0.5) * kVox, ((best & 0x1FFFFF) - kO + 0.5) * kVox};
    std::vector<size_t> near;
    for (int rep = 0; rep < 6; rep++) {
      near.clear();
      for (size_t n : fn)
        if (left[n] && off(q, n) < (rep < 3 ? 2.0 : 0.7)) near.push_back(n);
      if (near.size() < 10) break;
      M3 A;
      A.m[0][0] = A.m[1][1] = A.m[2][2] = 0;
      double bq[3] = {0, 0, 0};
      for (size_t n : near) {
        double dv[3] = {r.D[n].x, r.D[n].y, r.D[n].z}, ov[3] = {r.O[n].x, r.O[n].y, r.O[n].z};
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++) {
            double pij = (i == j ? 1 : 0) - dv[i] * dv[j];
            A.m[i][j] += pij;
            bq[i] += pij * ov[j];
          }
      }
      M3 Ai;
      if (!Inv3(A, Ai)) { near.clear(); break; }
      q = Ai * V3{bq[0], bq[1], bq[2]};
    }
    for (size_t n : fn)
      if (left[n] && off(q, n) < 2.0) left[n] = 0;
    if ((int)near.size() < kLight) continue;
    // seen from directions kPar apart or more: its depth is known
    V3 m{};
    std::vector<V3> w;
    for (size_t n : near) {
      V3 e = r.O[n] - q;
      w.push_back(e * (1 / norm(e)));
      m += w.back();
    }
    m = m * (1 / norm(m));
    std::vector<double> pa;
    for (auto &e : w) pa.push_back(std::acos(std::max(-1.0, std::min(1.0, dot(e, m)))) * kDeg);
    std::sort(pa.begin(), pa.end());
    double par = pa[pa.size() * 9 / 10];
    std::set<std::pair<int, double>> sq;
    for (size_t n : fn)
      if (off(q, n) < kNear) sq.insert({r.C[n], r.G[n]});
    std::vector<double> sc;
    bool scored = SeenScores(view(q, nullptr, sq), tu, P, sc);
    int i = scored ? named(sc) : -1;
    if (i < 0 || par < kPar || std::find(gone.begin(), gone.end(), (size_t)i) == gone.end()) continue;
    size_t k = u[i];
    double worst = 0, moved = norm(q - Pq[k]);
    std::string js, ds, dc;
    for (size_t j = 0; j < ok.size(); j++) {
      size_t kj = u[ok[j]];
      double a = norm(S[k] - S[kj]), b = norm(q - Pq[kj]);
      worst = std::max(worst, std::fabs(b - a));
      const char *sep = j == 0 ? "" : j + 1 == ok.size() ? " and " : ", ";
      js += sep + keys[kj];
      ds += sep + Fmt("%.2f", a);
      dc += sep + Fmt("%.2f", b);
    }
    if (worst < kMoved) continue;
    auto said = layout_said_.find(keys[k]);
    if (said != layout_said_.end() && std::fabs(said->second - moved) < kMoved) continue;
    layout_said_[keys[k]] = moved;
    log_("base station layout: SteamVR has " + keys[k] + " " + ds + " m from " + js + "; the cameras see it " + dc +
         " m from them" + Fmt(" (%.2f m from where SteamVR has it)", moved) +
         ": SteamVR's map of its base stations is off");
  }
}

// call ~1 Hz. The solver's clock is the newest sighting: with nothing new seen nothing ages out or is refitted
StepStat Solver::Step(double /*now*/) {
  std::vector<std::string> keys;
  std::vector<V3> S, Z;
  Stations(keys, S, Z);
  double newest;
  {
    std::lock_guard<std::mutex> g(m_);
    newest = rt_.empty() ? -kInf : rt_.back();
  }
  if (!(newest > seen_)) {
    StepStat st = has_stat_ ? stat_ : StepStat();
    st.jump = false; st.has_acq = false; st.cond = false;
    stat_ = st; has_stat_ = true;
    return st;
  }
  double now = seen_ = newest;
  // starved_ must be kept up to date even while unacquired: the candidate enumeration's ScoreX
  // counts dim rays only when the room is too light for them to saturate, and on the Frame this is
  // the only place it is set (RANSAC is retired). Left stale-false, a lit room's station dots stay
  // "dim" and the enumeration scores nothing. Set it BEFORE the manual-acq return below.
  starved_ = Starved(GetRays(std::max(now - ACQ_WIN, since_)));
  if (manual_acq_ && !has_x_) return StepStat();   // the RANSAC acquirer is retired on the Frame, cold AND re-acquire (see the last_acq_ block below): x is only ever set by the candidate enumeration (Sync::AcquireFromCandidates); if that cannot earn an unambiguous answer, the system stays unacquired -- held is better than lied to
  StepStat st;
  if (!S.empty() && has_x_) {
    Rays r = GetRays(std::max(now - WIN, since_));
    if (r.size()) {
      std::vector<V3> P, Zq;
      Predict(x_, S, Z, P, Zq);
      std::vector<int> k(r.size());
      std::vector<double> ak(r.size());
      for (size_t n = 0; n < r.size(); n++) Nearest(P, Zq, r.O[n], r.D[n], k[n], ak[n]);
      for (size_t kk = 0; kk < S.size(); kk++) {  // a station's recent sightings moved away from its older ones:
        std::vector<double> rec, old;          // the Quest's tracking jumped, forget the older sightings
        for (size_t n = 0; n < r.size(); n++)
          if (k[n] == (int)kk && ak[n] < (r.B[n] ? GATE : DimGate())) (r.T[n] > now - 4 ? rec : old).push_back(ak[n]);
        if (rec.size() >= 20 && old.size() >= 20 && Median(rec) - Median(old) > JUMP) {
          since_ = now - 4;
          st.jump = true;
          r = GetRays(since_);
          break;
        }
      }
      FitR f = Fit(x_, S, Z, r, now, GATE, has_anchor_ ? &anchor_ : nullptr, nullptr, DimGate());
      if (f.ok) {
        // BOUNDED REFINEMENT -- lighthouses are anchors: they are the fixed points the whole
        // system exists around, so the fit may polish centimeters but may never wander. A
        // correction past the slew (5 cm/update) or drift past the cap (30 cm from the blessed
        // anchor) is held; a run of held corrections is the lock itself confessing -- drop x and
        // re-earn through the candidate pipeline. No bound on re-earning: there is no motion to
        // model, only lies to evict.
        double dstep = Apart(f.x, x_, S), danchor = has_anchor_ ? Apart(f.x, anchor_, S) : 0.0;
        if (now - x_since_ > 10 && (dstep > 0.05 || danchor > 0.30)) {
          if (++sat_count_ == 1)
            log_(Fmt("the fit pushed past its bound (%.0f cm step, %.0f cm from the anchor): holding the old pose", dstep * 100, danchor * 100));
          if (sat_count_ >= 10) {
            log_("the fit kept pushing past its bound for 10 s: the lock is wrong -- re-earning through the candidate pipeline");
            sat_count_ = 0;
            PoseBreak(now);
          }
        } else {
          sat_count_ = 0;
          x_ = f.x;
          if (f.cond || !has_anchor_) { anchor_ = f.x; has_anchor_ = true; }
        }
      }
      st.cond = f.cond;
      std::vector<int> cnt;
      Support(x_, S, Z, r, INLIER, cnt);
      Predict(x_, S, Z, P, Zq);
      std::vector<double> in;
      for (size_t n = 0; n < r.size(); n++) {
        int kk; double a;
        Nearest(P, Zq, r.O[n], r.D[n], kk, a);
        if (a < INLIER) {
          in.push_back(a);
          double &l = st.last[keys[kk]];
          l = std::max(l, r.T[n]);
        }
      }
      st.n = f.n;
      st.has_per = true;
      for (size_t kk = 0; kk < keys.size(); kk++) st.per[keys[kk]] = cnt[kk];
      if (!in.empty()) { st.has_med = true; st.med = Median(in); }
    }
  }
  if (has_x_ && now - last_id_ >= kIdEvery) {
    last_id_ = now;
    CheckIdentity(now, keys, S, Z);
  }
  if (has_x_ && now - last_layout_ >= kLayoutEvery) {
    last_layout_ = now;
    CheckLayout(keys, S, Z);
  }
  bool locked = false;
  if (has_x_ && st.has_per) {
    std::vector<int> c;
    for (auto &kv : st.per) c.push_back(kv.second);
    locked = Score(c) >= 10;
  }
  if (!manual_acq_ && now - last_acq_ > (locked ? 10 : 3)) {
    // manual_acq_: on the Frame the RANSAC path is retired for re-acquisition too, not just cold
    // start. Once a lock exists it can only be re-earned through the candidate pipeline
    // (Sync::AcquireFromCandidates); this block would otherwise be a second, ungated acquirer
    // quietly Reset()ting x behind the enumeration's back. Quest (manual_acq_ false) keeps stock.
    last_acq_ = now;
    X4 xa;
    int sa, tight;
    if (Acquire(now, xa, sa, tight)) {
      int cur = has_x_ ? ThinnedScore(x_, now) : 0;
      st.has_acq = true; st.acq_sa = sa; st.acq_cur = cur; st.acq_tight = tight;
      // the mirror check's flip has the worn devices for evidence as well: half the support will do. A fit far off the
      // current one (four stations in a square fit about as well turned by 90 deg) goes by the worn or held devices
      // too; without them it takes twice the support, or tight sightings: a wrong fit that lines stations up with other
      // stations or lamps gathers a dozen, few of them tight, and after a pose break the current one has few yet
      int k = acq_dim_ ? 2 : 1;  // dim rays: other lights line up with a wrong alignment more easily
      double dn, dc;
      std::string how;
      bool far = has_x_ && !acq_forced_ && Apart(xa, x_, S) > 1.0;
      int fav = far ? BodyFavours(xa, x_, dn, dc, how) : 0;
      int need = far && !fav ? 2 * ACQ : ACQ;
      bool enough = sa >= k * need || tight >= k * TIGHT_N || (acq_forced_ && sa >= k * ACQ / 2);
      if (!enough || !(sa > 2 * cur + 5 || acq_forced_)) {
      } else if (fav == 2) {
        if (now - kept_said_ > 60) {
          kept_said_ = now;
          log_(Fmt("acquisition: yaw %+.1f deg (support %d, was %d) would fit as well, but %s: the current fit stays",
                   xa[0] * kDeg, sa, cur, how.c_str()));
        }
      } else {
        log_(Fmt("acquired: yaw %+.2f deg t [%.3f %.3f %.3f] support %d (%d within %.1f deg, was %d)%s", xa[0] * kDeg, xa[1],
                 xa[2], xa[3], sa, tight, TIGHT_DEG, cur, acq_forced_ ? ", by the mirror check" : ""));
        Reset(xa);
        since_ = brk_;
      }
    }
  }
  stat_ = st; has_stat_ = true;
  return st;
}

// The 4-DOF fit's stations tilted about pivot as well, on the same sightings. Its first gate is wide: with the frame
// 1.5 deg off level the fit of two stations puts the third 4 deg away from where the cameras see it. Kept only if the
// tilt is pinned (cond) by enough stations, each seen from several (head, direction) cells within INLIER.
bool Solver::Level(V3 pivot, LevelR &out) {
  if (!has_x_) return false;
  std::vector<std::string> keys;
  std::vector<V3> S0, Z0;
  Stations(keys, S0, Z0);
  if ((int)S0.size() < LEVEL_STATIONS) return false;
  const double now = seen_;
  Rays r = GetRays(std::max(now - WIN, since_));
  X4 x = x_;
  M3 L;
  std::vector<V3> S(S0.size()), Z(Z0.size());
  auto level = [&] {
    for (size_t k = 0; k < S0.size(); k++) { S[k] = pivot + L * (S0[k] - pivot); Z[k] = L * Z0[k]; }
  };
  FitR f;
  for (double gate : {6.0, GATE, 0.6}) {
    level();
    f = Fit(x, S, Z, r, now, gate, nullptr, &pivot, DimGate());
    if (!f.ok) return false;
    x = f.x;
    L = f.tilt * L;
  }
  if (!f.cond) return false;
  level();
  std::vector<int> cells;
  Support(x, S, Z, Thin(r), INLIER, cells);
  int held = 0;
  for (int c : cells) held += c >= LEVEL_CELLS;
  std::vector<V3> P, Zq;
  Predict(x, S, Z, P, Zq);
  std::vector<double> in;
  for (size_t n = 0; n < r.size(); n++) {
    int k; double a;
    Nearest(P, Zq, r.O[n], r.D[n], k, a);
    if (a < INLIER) in.push_back(a);
  }
  if (held < LEVEL_STATIONS || in.empty()) return false;
  double med = Median(in);
  if (med > LEVEL_MED) return false;
  out.x = x;
  out.tilt = L;
  out.stations = held;
  out.med = med;
  return true;
}

// With two stations the cameras can't see a tilt about the line through them (both stay put), but they do see one
// across it: that raises one station and lowers the other. The 4-DOF fit's stations tilted about axis through
// pivot, on the last WIN s of sightings, and on four time slices of them: their scatter says how well it's known.
bool Solver::LevelAxis(V3 pivot, V3 axis, AxisR &out) {
  if (!has_x_) return false;
  std::vector<std::string> keys;
  std::vector<V3> S0, Z0;
  Stations(keys, S0, Z0);
  if (S0.size() < 2) return false;
  const double now = seen_;
  Rays r = GetRays(std::max(now - WIN, since_));
  if (r.size() < 200) return false;
  auto fit = [&](const Rays &rr, double &th, double *med) {
    X4 x = x_;
    double t = 0;
    std::vector<V3> S(S0.size()), Z(Z0.size());
    for (double gate : {GATE, 0.6}) {
      M3 L = Tilt(t * axis.x, t * axis.z);
      for (size_t k = 0; k < S0.size(); k++) { S[k] = pivot + L * (S0[k] - pivot); Z[k] = L * Z0[k]; }
      FitR f = Fit(x, S, Z, rr, now, gate, nullptr, &pivot, DimGate(), &axis);
      if (!f.ok) return false;
      x = f.x;
      t += f.theta;
    }
    th = t;
    if (med) {  // and each station seen from enough (head, direction) cells
      M3 L = Tilt(t * axis.x, t * axis.z);
      for (size_t k = 0; k < S0.size(); k++) { S[k] = pivot + L * (S0[k] - pivot); Z[k] = L * Z0[k]; }
      std::vector<int> cells;
      Support(x, S, Z, Thin(rr), INLIER, cells);
      for (int c : cells)
        if (c < LEVEL_CELLS) return false;
      std::vector<V3> P, Zq;
      Predict(x, S, Z, P, Zq);
      std::vector<double> in;
      for (size_t n = 0; n < rr.size(); n++) {
        int k;
        double a;
        Nearest(P, Zq, rr.O[n], rr.D[n], k, a);
        if (a < INLIER) in.push_back(a);
      }
      if (in.empty()) return false;
      *med = Median(in);
    }
    return true;
  };
  if (!fit(r, out.theta, &out.med) || out.med > LEVEL_MED) return false;
  std::vector<double> ts;
  double t0 = r.T.front(), t1 = r.T.back();
  for (int q = 0; q < 4; q++) {
    Rays s;
    for (size_t n = 0; n < r.size(); n++)
      if (r.T[n] >= t0 + (t1 - t0) * q / 4 && r.T[n] < t0 + (t1 - t0) * (q + 1) / 4 + (q == 3 ? 1 : 0)) {
        s.T.push_back(r.T[n]); s.G.push_back(r.G[n]); s.O.push_back(r.O[n]); s.D.push_back(r.D[n]);
        s.C.push_back(r.C[n]); s.B.push_back(r.B[n]);
      }
    double th;
    if (s.size() >= 50 && fit(s, th, nullptr)) ts.push_back(th);
  }
  out.slices = (int)ts.size();
  if (ts.size() < 3) return false;
  double m = 0, v = 0;
  for (double t : ts) m += t;
  m /= ts.size();
  for (double t : ts) v += (t - m) * (t - m);
  out.sd = std::sqrt(v / (ts.size() - 1) / ts.size());
  return true;
}

void Solver::Relevel(const X4 &x) {
  x_ = anchor_ = x;
  has_x_ = has_anchor_ = true;
  has_acq_x_ = false;  // found on the stations before the level
}

// ================================================================ timing
void Timing::Record(double tg, V3 od, V3 dd, double hg, int cam, bool bright) {
  rec_.push_back({tg, od, dd, hg, cam, bright});
  while (!rec_.empty() && (rec_.front().tg < tg - 120 || rec_.size() > 20000)) rec_.pop_front();
}

// The alignment to test timings against, refit from slow-head sightings only (a timing error moves those < 0.1 deg):
// the solver's own alignment was shaped by the old timing, or was just seeded from the last session, and holding it
// pulled estimates toward it (3 ms instead of ~18 right after a seed). IRLS Gauss-Newton, Cauchy 0.3 deg, a weak pull
// toward the start for directions the sightings don't pin. False unless both stations have slow sightings.
static bool RefitSlow(X4 &x, const std::vector<V3> &S, const std::vector<V3> &Z, const std::vector<V3> &O,
                      const std::vector<V3> &D) {
  if (S.size() < 2 || O.size() < 30) return false;
  const X4 x0 = x;
  std::vector<V3> P, Zq;
  Solver::Predict(x, S, Z, P, Zq);
  std::vector<int> ks;
  std::vector<size_t> use;
  std::vector<int> per(S.size(), 0);
  for (size_t i = 0; i < O.size(); i++) {
    int k;
    double a;
    Solver::Nearest(P, Zq, O[i], D[i], k, a);
    if (!(a < 1.0)) continue;
    use.push_back(i);
    ks.push_back(k);
    per[k]++;
  }
  for (int c : per)
    if (c < 10) return false;
  const double r0 = 0.3 / kDeg, c2 = r0 * r0;
  const double sig[4] = {0.5 / kDeg, 0.05, 0.05, 0.05};  // the pull toward the start: one sighting's worth at this far
  auto res = [&](const X4 &y, std::vector<V3> &E) {
    std::vector<V3> Py, Zy;
    Solver::Predict(y, S, Z, Py, Zy);
    E.resize(use.size());
    for (size_t j = 0; j < use.size(); j++) {
      V3 U = Py[ks[j]] - O[use[j]];
      E[j] = cross(U * (1 / norm(U)), D[use[j]]);
    }
  };
  std::vector<V3> E, Ep, Em;
  for (int it = 0; it < 10; it++) {
    res(x, E);
    double H[4][4] = {}, g[4] = {};
    std::vector<double> w(E.size());
    for (size_t j = 0; j < E.size(); j++) w[j] = 1 / (1 + dot(E[j], E[j]) / c2);
    std::vector<std::array<V3, 4>> J(E.size());
    for (int q = 0; q < 4; q++) {
      const double h = q == 0 ? 1e-6 : 1e-5;
      X4 xp = x, xm = x;
      xp[q] += h;
      xm[q] -= h;
      res(xp, Ep);
      res(xm, Em);
      for (size_t j = 0; j < E.size(); j++) J[j][q] = (Ep[j] - Em[j]) * (1 / (2 * h));
    }
    for (size_t j = 0; j < E.size(); j++)
      for (int a = 0; a < 4; a++) {
        g[a] += w[j] * dot(J[j][a], E[j]);
        for (int b = 0; b < 4; b++) H[a][b] += w[j] * dot(J[j][a], J[j][b]);
      }
    for (int a = 0; a < 4; a++) {
      double pw = c2 / (sig[a] * sig[a]);
      H[a][a] += pw;
      g[a] += pw * (a == 0 ? Wrap(x[0] - x0[0]) : x[a] - x0[a]);
    }
    double dx[4];  // solve H dx = -g (Gaussian elimination, partial pivoting)
    for (int a = 0; a < 4; a++) {
      int p = a;
      for (int b = a + 1; b < 4; b++)
        if (std::fabs(H[b][a]) > std::fabs(H[p][a])) p = b;
      if (std::fabs(H[p][a]) < 1e-18) return false;
      if (p != a) {
        for (int c = 0; c < 4; c++) std::swap(H[a][c], H[p][c]);
        std::swap(g[a], g[p]);
      }
      for (int b = a + 1; b < 4; b++) {
        double f = H[b][a] / H[a][a];
        for (int c = a; c < 4; c++) H[b][c] -= f * H[a][c];
        g[b] -= f * g[a];
      }
    }
    for (int a = 3; a >= 0; a--) {
      double s = -g[a];
      for (int c = a + 1; c < 4; c++) s -= H[a][c] * dx[c];
      dx[a] = s / H[a][a];
    }
    for (int a = 0; a < 4; a++) x[a] += dx[a];
    if (std::fabs(dx[0]) < 1e-7 && std::fabs(dx[1]) + std::fabs(dx[2]) + std::fabs(dx[3]) < 1e-5) break;
  }
  return std::isfinite(x[0] + x[1] + x[2] + x[3]);
}

// The time a sighting's HMD pose is looked up at is (frame grid point - offset). The offset is the camera's own
// delay plus how the streamer times the poses it hands SteamVR (some predict ahead). A wrong offset shows only while
// the head turns: grid search the robust error of fast-head sightings against an alignment from the slow ones. A
// later estimate searches near the current timing, and the whole range when its minimum lies at that search's edge:
// a timing saved far off (a Frame's at -6.2 ms) could otherwise never move
bool Timing::Estimate(const PoseHist &poses, const X4 &x, const std::vector<V3> &S, const std::vector<V3> &Z, double cur,
                      bool wide, double &best, int &n, Split &sp) {
  bool edge = false;
  if (Search(poses, x, S, Z, cur, wide, best, n, sp, edge)) return true;
  return !wide && edge && Search(poses, x, S, Z, cur, true, best, n, sp, edge);
}

bool Timing::Search(const PoseHist &poses, const X4 &x, const std::vector<V3> &S, const std::vector<V3> &Z, double cur,
                    bool wide, double &best, int &n, Split &sp, bool &edge) {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  sp = {nan, nan, 0, 0, nan};
  edge = false;
  n = 0;
  if (rec_.empty() || S.size() < 2) return false;
  double lo = wide ? -0.04 : cur - 0.015, hi = wide ? 0.10 : cur + 0.015, step = wide ? 0.0025 : 0.001;
  std::vector<const Rec *> fast;
  std::vector<V3> so, sd;
  double newest = rec_.back().tg;
  for (auto &r : rec_) {
    if (r.tg < newest - 120) continue;
    double w;
    M3 R; V3 p;
    if (!poses.Speed(r.tg - cur, w)) continue;
    if (w < 10) {
      if (!poses.Still(r.tg - cur) && poses.At(r.tg - cur, R, p)) { so.push_back(p + R * r.od); sd.push_back(R * r.dd); }
      continue;
    }
    if (w < 15 || w > 150) continue;
    if (!poses.At(r.tg - lo, R, p) || !poses.At(r.tg - hi, R, p)) continue;
    fast.push_back(&r);
  }
  X4 xs = x;
  if (!RefitSlow(xs, S, Z, so, sd)) return false;
  std::vector<V3> P, Zq;
  Solver::Predict(xs, S, Z, P, Zq);
  // only sightings of a base station count: other lights would pad n and flatten the error curve. kNearDeg is wide
  // enough for any timing in the search at 150 deg/s
  std::vector<const Rec *> near;
  for (const Rec *r : fast) {
    M3 R; V3 p;
    if (!poses.At(r->tg - cur, R, p)) continue;
    int k; double a;
    Solver::Nearest(P, Zq, p + R * r->od, R * r->dd, k, a);
    if (std::isfinite(a) && a < kNearDeg) near.push_back(r);
  }
  fast.swap(near);
  n = (int)fast.size();
  // the first (wide) estimate wants plenty of turning
  if (n < (wide ? kTimingFirstN : kTimingNextN)) return false;
  // each sighting's error at every timing tried, and which way the head turned about the vertical then
  std::vector<double> es;
  for (double e = lo; e <= hi + 1e-9; e += step) es.push_back(e);
  const size_t m = es.size();
  const double s2 = 0.3 * 0.3;
  std::vector<double> cost(fast.size() * m);
  std::vector<int> way(fast.size(), 0);
  for (size_t i = 0; i < fast.size(); i++) {
    const Rec *r = fast[i];
    Quat qa, qb;
    V3 pa, pb;
    if (poses.AtQ(r->tg - cur - 0.02, qa, pa) && poses.AtQ(r->tg - cur + 0.02, qb, pb)) {
      Quat d = qb * Quat{qa.w, -qa.x, -qa.y, -qa.z};  // the turn between, in the room
      double yaw = 2 * std::asin(std::max(-1.0, std::min(1.0, d.w < 0 ? -d.y : d.y))) * kDeg / 0.04;  // deg/s
      way[i] = yaw > kTimingTurn ? 1 : yaw < -kTimingTurn ? -1 : 0;
    }
    for (size_t j = 0; j < m; j++) {
      M3 R; V3 p;
      double &c = cost[i * m + j];
      if (!poses.At(r->tg - es[j], R, p)) { c = 1; continue; }
      int k; double a;
      Solver::Nearest(P, Zq, p + R * r->od, R * r->dd, k, a);
      c = std::isfinite(a) ? a * a / (a * a + s2) : 1;
    }
  }
  // the error curve's minimum on the sightings turning one way (1 left, -1 right, 0 all): none at the edge of the
  // search, or when it doesn't stand out of the error floor (a flat curve: no fast sightings of stations)
  auto at = [&](const std::vector<double> &cs, size_t &i) {  // the curve's lowest point, between grid steps
    i = std::min_element(cs.begin(), cs.end()) - cs.begin();
    if (i == 0 || i + 1 == m) return es[i];
    double c0 = cs[i - 1], c1 = cs[i], c2 = cs[i + 1], den = c0 - 2 * c1 + c2;
    double sh = den > 0 ? 0.5 * (c0 - c2) / den : 0;
    return es[i] + std::max(-1.0, std::min(1.0, sh)) * step;
  };
  auto minimum = [&](int w, double &b, int &k) {
    std::vector<double> cs(m, 0.0);
    k = 0;
    for (size_t i = 0; i < fast.size(); i++)
      if (!w || way[i] == w) {
        k++;
        for (size_t j = 0; j < m; j++) cs[j] += cost[i * m + j];
      }
    if (k < (wide ? kTimingFirstN : kTimingNextN) / (w ? 4 : 1)) return false;
    size_t i;
    b = at(cs, i);
    if (i == 0 || i + 1 == m) { edge = true; return false; }
    double mx = *std::max_element(cs.begin(), cs.end());
    return mx - cs[i] > 0.05 * mx;
  };
  // a camera a little off its calibration looks like a timing error one way while the head turns left and the other
  // way turning right: the timing is the two turns' mean (a Frame recording with another Frame's calibration: -1 and
  // 20 ms, 9.6 ms, where all of it put 17.9; Quest recordings' turns lie 0.6-5.5 ms apart)
  int k0;
  bool all = minimum(0, best, k0), l = minimum(1, sp.left, sp.nl), r = minimum(-1, sp.right, sp.nr);
  if (!all || !l || !r) {
    sp.left = sp.right = nan;
    return false;
  }
  if (std::fabs(sp.left - sp.right) > kTimingSplit) return false;
  best = (sp.left + sp.right) / 2;
  // how well the sightings pin it: the mean again on resampled sightings, each way's drawn with replacement
  std::vector<size_t> L, Rt;
  for (size_t i = 0; i < fast.size(); i++)
    if (way[i]) (way[i] > 0 ? L : Rt).push_back(i);
  std::mt19937 rng(1);
  std::vector<double> bs;
  for (int b = 0; b < kTimingBoot; b++) {
    double side[2];
    for (int s = 0; s < 2; s++) {
      const std::vector<size_t> &I = s ? Rt : L;
      std::uniform_int_distribution<size_t> u(0, I.size() - 1);
      std::vector<double> cs(m, 0.0);
      for (size_t k = 0; k < I.size(); k++) {
        const double *c = &cost[I[u(rng)] * m];
        for (size_t j = 0; j < m; j++) cs[j] += c[j];
      }
      size_t i;
      side[s] = at(cs, i);
    }
    bs.push_back((side[0] + side[1]) / 2);
  }
  std::sort(bs.begin(), bs.end());
  sp.ci = bs[bs.size() * 95 / 100] - bs[bs.size() * 5 / 100];
  return sp.ci <= kTimingCI;
}

// ================================================================ the app
// The solver's recent sightings were placed with the old timing: place them again with the new one. Otherwise the
// next fit still leans on the old timing, and the next estimate (made against that fit) is pulled back toward it.
void Sync::Retime(double old_e, double new_e) {
  const auto &rec = timing_.recs();
  if (rec.empty()) return;
  std::vector<double> T, G;
  std::vector<V3> O, D;
  std::vector<int> C;
  std::vector<char> B;
  for (auto &r : rec) {
    double t = r.tg - new_e, w;
    M3 R;
    V3 p;
    if (!poses_.At(t, R, p) || !poses_.Speed(t, w) || poses_.Still(t) || w > kWmax) continue;
    T.push_back(t);
    O.push_back(p + R * r.od);
    D.push_back(R * r.dd);
    G.push_back(r.hg);
    C.push_back(r.cam);
    B.push_back(r.bright);
  }
  solver_.Replace(rec.front().tg - old_e, T, O, D, G, C, B);
}

static std::string Join(const std::string &dir, const char *name) {
  if (dir.empty()) return name;
  char e = dir.back();
  return dir + ((e == '\\' || e == '/') ? "" : "\\") + name;
}

static StationsFile LoadStations(const std::string &dir) {
  StationsFile f;
  f.Load(Join(dir, "stations.json"));
  return f;
}

Sync::Sync(SyncConfig cfg, LogFn log)
    : cfg_(std::move(cfg)), log_(std::move(log)), sfile_(LoadStations(cfg_.dir)), frame_(sfile_, log_),
      solver_(log_, 1, sfile_.autogen ? std::map<std::string, V3>() : sfile_.fix) {
  if (cfg_.arrival_clock) expo_ = kExpoArrival;
  solver_.SetBody([this](const X4 &x) { return BodyDist(x); });
  solver_.SetMotion([this](const X4 &x, double &moved) { return BodyMotion(x, moved); });
  solver_.SetInImage([this](int cam, V3 d) {
    std::lock_guard<std::mutex> g(net_);
    return have_optics_ && optics_.InImage(cam, d, 8.0);
  });
  if (sfile_.loaded)
    log_(Fmt("reference frame: anchor %s%s, %d measured station(s)%s", sfile_.anchor.c_str(),
             sfile_.autogen ? " (automatic)" : "", (int)sfile_.fix.size(),
             frame_.layout() ? Fmt(", the fit of %d stations", (int)sfile_.ref.size()).c_str() : ""));
  else
    log_("no stations.json yet: the reference frame is set the first time base stations show up");
  if (cfg_.per_headset) {
    seeded_ = true;  // this headset's state.json is read once it's known (ApplyHeadset)
  } else {
    state_path_ = Join(cfg_.dir, "state.json");
    LoadState();
  }
  // the camera-measured 'place' prior is retired (see PlaceStep): SteamVR's layout is the prior
}

Sync::~Sync() {}

void Sync::LoadState() {
  quest_.clear();
  {
    std::lock_guard<std::mutex> g(st_m_);
    timings_.clear();
  }
  saved_cs_.clear();
  saved_timings_.clear();
  has_saved_ = false;
  seeded_ = true;
  std::string text;
  JVal d;
  if (state_path_.empty() || !ReadFile(state_path_, text) || !JParse(text, d)) return;
  if (const JVal *q = d.get("quest_stations"))
    for (auto &kv : q->o) {
      auto v = kv.second.nums();
      if (v.size() >= 3) quest_[kv.first] = {v[0], v[1], v[2]};
    }
  // timings checked as Timing::Estimate does since 1.13 ("timing", before, could be a camera's error or noise: one sat
  // at -6.2 ms for a Frame and kept it from ever finding the base stations)
  if (const JVal *t = d.get("pose_timing")) {
    std::lock_guard<std::mutex> g(st_m_);
    for (auto &kv : t->o)
      if (kv.second.t == JVal::Num) timings_[kv.first] = kv.second.n;
  }
  seeded_ = quest_.empty();
}

void Sync::SaveState(const X4 &x, const std::map<std::string, V3> &cs, bool force) {
  if (state_path_.empty()) return;  // the headset isn't known yet
  bool moved = force || !has_saved_ || std::fabs(x[0] - saved_x_[0]) * kDeg > 0.01 || cs.size() != saved_cs_.size() ||
               timings_.size() != saved_timings_.size();
  for (int i = 1; i < 4; i++) moved = moved || std::fabs(x[i] - saved_x_[i]) > 0.001;
  for (auto &kv : cs) {
    auto it = saved_cs_.find(kv.first);
    moved = moved || it == saved_cs_.end() || norm(kv.second - it->second) > 0.001;
  }
  for (auto &kv : timings_) {
    auto it = saved_timings_.find(kv.first);
    moved = moved || it == saved_timings_.end() || std::fabs(kv.second - it->second) > 0.0001;
  }
  if (!moved) return;  // under 0.01 deg / 1 mm / 0.1 ms: the file already has it
  saved_x_ = x;
  saved_cs_ = cs;
  saved_timings_ = timings_;
  has_saved_ = true;
  std::string s = "{\n \"quest_stations\": {";
  bool first = true;
  M3 R = Ry(x[0]);
  for (auto &kv : cs) {
    V3 q = R * kv.second + V3{x[1], x[2], x[3]};
    s += Fmt("%s\n  \"%s\": [%.5f, %.5f, %.5f]", first ? "" : ",", kv.first.c_str(), q.x, q.y, q.z);
    first = false;
  }
  s += Fmt("%s},\n \"x\": [%.9f, %.9f, %.9f, %.9f],\n \"pose_timing\": {", cs.empty() ? "" : "\n ", x[0], x[1], x[2], x[3]);
  first = true;
  for (auto &kv : timings_) { s += Fmt("%s\n  \"%s\": %.5f", first ? "" : ",", kv.first.c_str(), kv.second); first = false; }
  s += Fmt("%s},\n \"saved\": \"%s\"\n}\n", timings_.empty() ? "" : "\n ", Now().c_str());
  WriteFileAtomic(state_path_, s);
}

// 4-DOF fit raw -> Quest from the stations both know
bool Sync::Seed(const std::map<std::string, V3> &raw, X4 &x, double &miss) {
  std::vector<V3> S, Q;
  for (auto &kv : raw) {
    auto it = quest_.find(kv.first);
    if (it != quest_.end()) { S.push_back(kv.second); Q.push_back(it->second); }
  }
  if (S.size() < 2) return false;
  x = YawFit(S, Q);
  V3 t{x[1], x[2], x[3]};
  miss = 0;
  for (size_t i = 0; i < S.size(); i++) miss = std::max(miss, norm(Ry(x[0]) * S[i] + t - Q[i]));
  return true;
}

void Sync::SetUpperOnly(bool on) {
  std::lock_guard<std::mutex> g(net_);
  if (on == optics_.upper_only()) return;
  optics_.SetUpperOnly(on);
  if (!have_optics_) return;
  std::string left = optics_.Unused();
  if (!left.empty())
    log_("cameras " + left + " left out: they look down, the base stations hang above");
  else
    log_("cameras: all of them count again (the Frame's lower pair included)");
}

bool Sync::SetCalibration(const std::string &json, std::string *err) {
  std::lock_guard<std::mutex> g(net_);
  if (!optics_.Load(json, err)) return false;
  have_optics_ = true;
  std::string left_out = optics_.Unused();
  if (!left_out.empty()) log_("cameras " + left_out + " left out: they look down, the base stations hang above");
  else if (optics_.exact_time() && !optics_.upper_only()) log_("cameras: every camera counts (the upper-only cut is off)");
  // a Quest's frame period is one of two, the Frame's is learned from its frames
  grid_quest_ = !optics_.exact_time() && !cfg_.learn_grid;
  grid_ = FrameGrid(grid_quest_);
  grid_logged_.clear();
  return true;
}

void Sync::HeadsetReset() {
  std::lock_guard<std::mutex> g(net_);
  clock_.Clear();
  grid_ = FrameGrid(grid_quest_);
  grid_logged_.clear();
  waiting_.clear();
}

void Sync::SetStreamer(const std::string &system) {
  std::lock_guard<std::mutex> g(st_m_);
  if (system == streamer_) return;
  streamer_ = system;
  auto it = timings_.find(system);
  if (it != timings_.end() && cfg_.learn_timing) {
    expo_ = it->second;
    timing_learned_ = true;
    log_(Fmt("timing for %s: %.1f ms (learned before)", system.c_str(), it->second * 1000));
  }
}

// "Steam-Frame-1a2b3c4d": the model for the eye, a hash of the serial to tell two of a kind apart
void Sync::RecGravity(double t, bool on, const M3 &tilt) {
  if (!recording()) return;
  Rec(t, "G %d %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f", on ? 1 : 0, tilt.m[0][0], tilt.m[0][1], tilt.m[0][2],
      tilt.m[1][0], tilt.m[1][1], tilt.m[1][2], tilt.m[2][0], tilt.m[2][1], tilt.m[2][2]);
}

void Sync::SetHeadset(const std::string &serial, const std::string &model) {
  std::string name;
  for (char c : model.substr(0, 32)) {
    if (std::isalnum((unsigned char)c)) name += c;
    else if (!name.empty() && name.back() != '-') name += '-';
  }
  while (!name.empty() && name.back() == '-') name.pop_back();
  if (name.empty()) name = "headset";
  uint32_t h = 2166136261u;  // FNV-1a
  for (unsigned char c : serial) { h ^= c; h *= 16777619u; }
  std::string key = serial.empty() || serial == "unknown" ? name : Fmt("%s-%08x", name.c_str(), h);
  std::lock_guard<std::mutex> g(hs_m_);
  hs_pending_ = key;
  hs_new_ = true;
}

// A headset other than the one before has its own state.json (headsets/<key>/): what the last one's cameras saw says
// nothing of this one's tracking space, or of how its streamer times its poses. The single state.json of earlier
// versions goes to the first headset that connects.
void Sync::ApplyHeadset() {
  if (!hs_new_.exchange(false)) return;
  std::string key;
  {
    std::lock_guard<std::mutex> g(hs_m_);
    key = hs_pending_;
  }
  if (key == hs_key_) return;  // the same headset again: it reconnected
  namespace fs = std::filesystem;
  const bool first = hs_key_.empty();
  hs_key_ = key;
  std::error_code ec;
  fs::path dir = fs::path(cfg_.dir) / "headsets" / key, legacy = fs::path(cfg_.dir) / "state.json";
  fs::create_directories(dir, ec);
  const std::string path = (dir / "state.json").string();
  bool adopted = false;
  if (first && !fs::exists(path, ec) && fs::exists(legacy, ec)) {
    fs::rename(legacy, path, ec);
    adopted = !ec;
  }
  state_path_ = path;
  LoadState();
  if (!first) ClearAlignment();  // the last headset's sightings and fit mean nothing to this one
  {
    std::lock_guard<std::mutex> g(st_m_);
    ApplyTiming();
  }
  log_(Fmt("state for %s: headsets\\%s\\state.json%s", key.c_str(), key.c_str(),
           adopted ? " (the state.json of earlier versions, moved there)"
                   : quest_.empty() ? " (nothing saved yet)" : Fmt(" (%d base stations saved)", (int)quest_.size()).c_str()));
}

// expo_ for the current streamer from timings_ (st_m_ held); the default when it has none
void Sync::ApplyTiming() {
  auto it = timings_.find(streamer_);
  if (!streamer_.empty() && it != timings_.end() && cfg_.learn_timing) {
    expo_ = it->second;
    timing_learned_ = true;
    log_(Fmt("timing for %s: %.1f ms (learned before)", streamer_.c_str(), it->second * 1000));
  } else {
    expo_ = cfg_.arrival_clock ? kExpoArrival : kExpoDefault;
    timing_learned_ = false;
  }
}

// the sightings, the fit and what hung on them: gone, the lighthouse devices uncorrected until the next acquisition
void Sync::ClearAlignment() {
  solver_.Clear();
  {
    std::lock_guard<std::mutex> g(st_m_);
    last_st_ = StepStat();
    has_last_st_ = false;
  }
  have_applied_ = false;
  applied_ = X4{};
  lock_time_ = -1;
  resets_seen_ = solver_.resets();
  sum_t_ = -1;
  sum_n_ = 0;
}

void Sync::ForceRecalibrate() {
  ClearAlignment();
  quest_.clear();
  seeded_ = true;
  {
    std::lock_guard<std::mutex> g(st_m_);
    timings_.clear();
    ApplyTiming();
  }
  // an automatic reference frame is pinned again from SteamVR's stations as they are now (and its level found again);
  // a measured one is the user's own and stays
  if (sfile_.autogen || !sfile_.loaded) {
    for (auto &kv : sfile_.place) solver_.DropFix(kv.first);
    sfile_ = StationsFile();
    frame_.Reset();
    std::remove(Join(cfg_.dir, "stations.json").c_str());
    place_prior_.clear();
    place_sess_.clear();
    place_off_.clear();
    left_seen_.clear();
    left_since_ = 0;
    bad_since_ = tilt_since_ = -1;
    place_saved_ = -1e18;
  }
  SaveState(X4{}, {}, true);  // the saved alignment goes too, or a restart would seed from it
  Rec(last_step_, "I recalibrate");
  log_("force recalibrate: everything learned about the alignment is dropped, starting over");
}

void Sync::OnHmdPose(double t, const Quat &q, V3 p) {
  poses_.Add(t, q, p);
  if (recording_ && t - last_prec_ >= 0.004) {  // recordings: the HMD at ~250 Hz like the old poller, 3x4 row-major
    last_prec_ = t;
    M3 R = ToM3(q);
    Rec(t, "P 0 %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f", R.m[0][0], R.m[0][1], R.m[0][2], p.x,
        R.m[1][0], R.m[1][1], R.m[1][2], p.y, R.m[2][0], R.m[2][1], R.m[2][2], p.z);
  }
}

void Sync::OnBodyPose(int dev, double t, V3 raw) {
  std::lock_guard<std::mutex> g(body_m_);
  double &last = body_last_[dev];
  if (t - last < 0.1) return;  // 10 Hz is plenty
  last = t;
  body_.push_back({t, dev, raw});
  if (recording_) body_new_.push_back({t, dev, raw});
  while (!body_.empty() && body_.front().t < t - Solver::ACQ_WIN) body_.pop_front();
}

// the median, over 0.1 s steps, of the nearest worn or held lighthouse device's horizontal distance to the head, were x
// the alignment (Tick's thread: frame_ is stable here)
double Sync::BodyDist(const X4 &x) {
  std::vector<BodyS> b;
  {
    std::lock_guard<std::mutex> g(body_m_);
    b.assign(body_.begin(), body_.end());
  }
  M3 C;
  V3 tc;
  frame_.Rotation(C, tc);
  double c = std::cos(x[0]), s = std::sin(x[0]);
  std::map<long long, double> nearest;
  for (const BodyS &e : b) {
    M3 R;
    V3 h;
    if (!poses_.At(e.t, R, h)) continue;
    V3 q = RyMul(c, s, C * e.p + tc) + V3{x[1], x[2], x[3]};
    double d = std::hypot(q.x - h.x, q.z - h.z);
    auto it = nearest.emplace((long long)std::floor(e.t * 10), d).first;
    it->second = std::min(it->second, d);
  }
  if (nearest.size() < 20) return NAN;
  std::vector<double> v;
  for (auto &kv : nearest) v.push_back(kv.second);
  return Median(v);
}

// How the worn or held lighthouse devices' horizontal motion follows the head's, were x the alignment: the
// correlation of their steps (about 0.3 s each) with the head's over the same times, all devices together. +1 every
// step along, -1 every step against; devices left lying add nothing. moved: the head's own path meanwhile (m).
double Sync::BodyMotion(const X4 &x, double &moved) {
  std::vector<BodyS> b;
  {
    std::lock_guard<std::mutex> g(body_m_);
    b.assign(body_.begin(), body_.end());
  }
  moved = 0;
  if (b.size() < 20) return NAN;
  M3 C;
  V3 tc;
  frame_.Rotation(C, tc);
  double c = std::cos(x[0]), s = std::sin(x[0]), num = 0, nh = 0, nd = 0;
  std::map<int, std::vector<const BodyS *>> per;
  for (const BodyS &e : b) per[e.dev].push_back(&e);
  for (auto &kv : per) {
    const std::vector<const BodyS *> &v = kv.second;
    for (size_t i = 0, j = 0; i < v.size(); i = j) {
      for (j = i + 1; j < v.size() && v[j]->t - v[i]->t < 0.25; j++) {}
      if (j >= v.size()) break;
      if (v[j]->t - v[i]->t > 0.5) continue;
      M3 R;
      V3 h0, h1;
      if (!poses_.At(v[i]->t, R, h0) || !poses_.At(v[j]->t, R, h1)) continue;
      V3 dq = RyMul(c, s, C * (v[j]->p - v[i]->p)), dh = h1 - h0;
      num += dh.x * dq.x + dh.z * dq.z;
      nh += dh.x * dh.x + dh.z * dh.z;
      nd += dq.x * dq.x + dq.z * dq.z;
    }
  }
  M3 R;
  V3 h0, h1;
  for (double t = b.front().t; t + 0.25 <= b.back().t; t += 0.25)
    if (poses_.At(t, R, h0) && poses_.At(t + 0.25, R, h1)) moved += std::hypot(h1.x - h0.x, h1.z - h0.z);
  return nh > 0 && nd > 0 ? num / std::sqrt(nh * nd) : NAN;
}

void Sync::SetRecord(FILE *f) {
  std::lock_guard<std::mutex> g(rec_m_);
  if (rec_) fclose(rec_);
  rec_ = f;
  recording_ = f != nullptr;
}

void Sync::Rec(double t, const char *fmt, ...) {
  if (!recording_) return;
  std::lock_guard<std::mutex> g(rec_m_);
  if (!rec_) return;
  fprintf(rec_, "%lld ", (long long)std::llround(t * 1e9));
  va_list ap;
  va_start(ap, fmt);
  vfprintf(rec_, fmt, ap);
  va_end(ap);
  fputc('\n', rec_);
}

FILE *g_ray_dump = nullptr;

// one lhsight line: "F cam k t_us mean nblobs [x y npx peak]..." (x, y in pixels * 10)
void Sync::OnLine(double pc, const char *line) {
  if (line[0] != 'F' || line[1] != ' ') return;
  Rec(pc, "%s", line);
  int cam, nb;
  unsigned k;
  long long tus;
  int mean, off = 0;
  if (sscanf(line, "F %d %u %lld %d %d%n", &cam, &k, &tus, &mean, &nb, &off) < 5) return;
  // a calibrated camera left out (the Frame's lower pair, which looks at the floor): not worth a frame grid, cells or
  // rays. An UNcalibrated camera isn't skipped -- the Quest's frame grid and timing still want its frames
  if (have_optics_ && optics_.Skip(cam)) return;
  // the class tag is LAST on the line ("... [x y npx peak]... L|S"), after the blobs -- not right after nb. Reading it
  // from just past nb only matched blob-less lines, so every line with dots was classed 'S': both classes then went
  // into the frame grid (its period could never fold: 49 frames/s measured), the L dots marked seqS, and seqL stayed
  // empty -- which also left the dot marks off both rings (the marking needs both non-empty). That is why the
  // presence rings read all-zero and the S-gate and the channel check were dead. Tags are one letter: take the last
  // token if it is exactly that and preceded by a space (a Quest's line ends on a blob's peak, a digit).
  char xcls = 0;
  {
    size_t L = strlen(line);
    while (L && line[L - 1] == ' ') L--;
    if (L >= 2 && line[L - 2] == ' ' && (line[L - 1] == 'L' || line[L - 1] == 'S')) xcls = line[L - 1];
  }
  std::lock_guard<std::mutex> g(net_);
  nframes_++;
  double hs = tus / 1e6;
  if (cfg_.arrival_clock) clock_.AddArrival(pc, hs);
  double lag = 0;
  // the Frame stamps its frames itself; the Quest's are when lhsight found them, the grid takes the poll lag off. On
  // the Frame the grid still learns the frame period, for the channel check. Tagged-S lines stay out of the grid too:
  // interleaved with L their 8.4 ms offset halves the median gap and Learn can never fold them (period stuck at 0 ->
  // CheckIdentity dead). L lines alone give the true frame period; S lines are evidence only (below).
  bool grid_lag = nb >= 0 && xcls != 'S' && grid_.Lag(cam, hs, lag);
  if (optics_.exact_time()) lag = 0;
  bool have_lag = nb >= 0 && (optics_.exact_time() || grid_lag);
  if (nb >= 0) {  // the frame period: log it when it's found or moves
    double p = grid_.period(cam), &was = grid_logged_[cam];
    if (p > 0 && std::fabs(p - was) > 0.0001) {
      log_(Fmt("camera %d: short frames every %.3f ms", cam, p * 1000));
      was = p;
    }
  }
  // The presence rings are NOT advanced here any more. They are advanced in Use, with the FRAME'S
  // OWN pose and only for cells whose direction was actually on that camera's image: see
  // Sync::MarkPresence. Here we only classify the line and queue it.
  char cls = xcls ? xcls : (nb > 0 ? 'S' : 'L');   // evidence class: untagged dots sit on short frames (Quest)
  bool solve = !xcls || xcls == 'L';               // untagged lines are stock Quest semantics: dots are rays
  std::vector<Spot> bl;
  {
    // how many numbers per blob: the Frame sends 8 (x10 y10 npx peak sat flux w h), recordings from
    // before that send 4. nb is on the line, so counting the remaining tokens settles it exactly --
    // guessing per-blob would let an 8-wide read swallow two 4-wide records.
    int per = 4, ntok = 0;
    {
      bool in = false;
      for (const char *q = line + off; *q; q++) {
        if (*q == ' ') in = false;
        else if (!in) in = true, ntok++;
      }
      if (xcls) ntok--;   // the class letter is the last token
      if (nb > 0 && ntok == nb * 8) per = 8;
    }
    const char *s = line + off;
    for (int i = 0; i < nb; i++) {
      int x10, y10, npx, peak, sat = 0, fl = 0, bw = 0, bh = 0, used = 0;
      int got = per == 8 ? sscanf(s, " %d %d %d %d %d %d %d %d%n", &x10, &y10, &npx, &peak, &sat, &fl, &bw, &bh, &used)
                         : sscanf(s, " %d %d %d %d%n", &x10, &y10, &npx, &peak, &used);
      if (got < 4) break;
      s += used;
      bl.push_back({x10, y10, npx, peak, sat, fl, bw, bh});
    }
  }
  nb = (int)bl.size();
  // stage 1: collection. What each camera and class actually delivered, before any gate touches it.
  if (cam >= 0 && cam < kMaxCam) {
    CamSpots &cs = camsp_[cam];
    int ci = (cls == 'S' ? 1 : 0);
    cs.fr[ci]++;
    if (nb) {
      cs.dots[ci]++;
      cs.blobs[ci] += nb;
      for (const Spot &b : bl) {
        if (b.peak > cs.pk[ci]) cs.pk[ci] = b.peak;
        int f = b.sum ? b.sum : b.npx * b.peak;
        cs.fl[ci] += f;
        if (f > cs.flmax[ci]) cs.flmax[ci] = f;
      }
    }
  }
  if (!solve) {
    // tagged-S dots are evidence, not rays: no grid, no solver. They are still queued, because the
    // class rings are marked from the frame's own pose like the long ones -- that pose is not known
    // until the frame comes round in Use, so an S line has to wait its turn like an L line.
    // Untagged lines (stock lhsight on Quest) are short-frame dots that DO solve -- only
    // lhsight2's explicit 'S' tag makes a line evidence-only.
    spots_.sframe += nb;
    if (!nb || !have_optics_ || !optics_.exact_time()) return;
    double sgrid;
    if (have_lag && clock_.Map(hs, sgrid)) {
      waiting_.push_back({pc, sgrid - expo_, hs, sgrid, cam, have_lag, 'S', std::move(bl)});
      while (!waiting_.empty()) {
        V3 lp;
        double lt;
        const Shot &f = waiting_.front();
        if (!(pc - f.pc >= kPoseWait || (poses_.Latest(lp, &lt) && lt > f.t))) break;
        Use(f);
        waiting_.pop_front();
      }
    }
    return;
  }
  spots_.frames++;
  spots_.seen += nb;
  if (!nb) return;   // dotless frames feed the grid above only
  if (!have_optics_) { spots_.other += nb; return; }
  double hg = have_lag ? hs - lag : std::numeric_limits<double>::quiet_NaN();  // frame time, headset clock
  if (have_lag && grid_.period(cam) > 0) solver_.SetFramePeriod(grid_.period(cam));
  double grid_pc;
  if (!clock_.Map(hs - (have_lag ? lag : 0.0), grid_pc)) { spots_.other += nb; return; }
  double expo = expo_;
  double t = grid_pc - expo - (have_lag ? 0.0 : kLagFallback);
  waiting_.push_back({pc, t, hg, grid_pc, cam, have_lag, 'L', std::move(bl)});
  // everything that reaches here solves: tagged-S was returned above, untagged
  // (Quest) dots solve regardless of their evidence class -- hence the literal
  // the frames whose head pose is in (one after their pose time), or that waited long enough, oldest first
  while (!waiting_.empty()) {
    V3 lp;
    double lt;
    const Shot &f = waiting_.front();
    if (!(pc - f.pc >= kPoseWait || (poses_.Latest(lp, &lt) && lt > f.t))) break;
    Use(f);
    waiting_.pop_front();
  }
}

// quantise a world direction to ~2.4 deg cells; camera in the key
uint64_t Sync::CellKey(int cam, const V3 &w) {
  uint64_t h = 1469598103934665603ull;
  auto mix = [&](uint64_t v) { h ^= v + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2); };
  mix((uint64_t)(uint32_t)cam);
  mix((uint64_t)(int64_t)std::lround(w.x * 24));
  mix((uint64_t)(int64_t)std::lround(w.y * 24));
  mix((uint64_t)(int64_t)std::lround(w.z * 24));
  return h;
}

void Sync::CellSeen(int cam, const V3 &w, float peak, double t) {
  uint64_t k = CellKey(cam, w);
  SpotCell &cl = cells_[k];
  if (!cl.inact) { cl.inact = true; actk_[cam].push_back(k); }
  if (cl.seen == 0) { cl.cam = cam; cl.w = w; cl.w0 = w; }
  else cl.w = cl.w + (w - cl.w) * 0.1;
  cl.seen++; cl.last = t;
  // ampL is NOT pushed here: it is advanced by MarkPresence, one sample per long frame, 0 when the
  // dot is absent. A run of only the hits is not a time series and has no rotor beat in it.
  cl.amax = std::max(cl.amax, peak);
  cl.amin = std::min(cl.amin, peak);
}

// The class presence rings: one sample per frame per direction cell whose direction was ON THIS
// CAMERA'S IMAGE that frame, set when one of the frame's blobs lands in the same cell.
//
// Two things this fixes, both measured (qlhs_presence on the 2026-10-11 recording):
//  - the pose. OnLine marked with the camera's "last known pose", which is the PREVIOUS frame's:
//    under head motion the mark lands one cell over, and 32% of all marks missed the cell they
//    belonged to (23,562 of 74,772). Here R is the frame's own pose, the same one its rays use.
//  - the denominator. OnLine counted a frame as coverage when the cell was merely within a ~101 deg
//    cone of the optical axis (dot(w, axis) > -0.2, i.e. in front of the camera, mostly off-image).
//    That is 56-67% of frames where the direction is really on the image only 27-40% of the time,
//    so every cell collected ~2x too many zeros. InImage is the honest test.
// The two together read a base station's own cell at ~4%; its real rate is 22-40% (L) and 5-7% (S).
// The amplitude ring is one byte per frame and the flux spans roughly 1e2..1e6 (a bloomed blob is
// thousands of saturated pixels), so compress it logarithmically -- the beat is a RELATIVE
// modulation, and 8 levels per doubling keeps ~0.4 dB of resolution across the whole range.
static uint8_t AmpFlux(int npx, int peak, int sum) {
  uint32_t f = sum ? (uint32_t)sum : (uint32_t)npx * (uint32_t)peak;  // old recordings: approximate
  if (!f) return 0;
  int a = 0;
  while (f > 32 && a < 255) { f >>= 1; a += 8; }
  return (uint8_t)a;
}

void Sync::MarkPresence(int cam, const M3 &R, double t, char cls, const std::vector<Spot> &bl) {
  (void)t;
  M3 Rc; V3 tc;
  if (!optics_.Pose(cam, Rc, tc)) return;
  // the cells this frame's blobs land in (key -> the brightest FLUX, for the amplitude ring)
  std::vector<std::pair<uint64_t, int>> hit;
  for (const Spot &b : bl) {
    V3 o, d;
    if (!optics_.Ray(cam, b.x10 / 10.0 - 0.5, b.y10 / 10.0 - 0.5, o, d)) continue;
    hit.push_back({CellKey(cam, R * d), AmpFlux(b.npx, b.peak, b.sum)});
  }
  const M3 Rt = T(R), Rct = T(Rc);
  for (uint64_t ak : actk_[cam]) {
    auto it = cells_.find(ak);
    if (it == cells_.end()) continue;
    SpotCell &cl = it->second;
    V3 dc = Rct * (Rt * cl.w);   // the cell's direction, in this camera's frame
    if (!optics_.InImage(cam, dc, 0)) continue;
    int pk = 0;
    bool on = false;
    for (auto &e : hit)
      if (e.first == ak) { on = true; if (e.second > pk) pk = e.second; }
    std::vector<uint8_t> &ring = cls == 'L' ? cl.seqL : cl.seqS;
    ring.push_back(on ? 1 : 0);
    if (ring.size() > 64) ring.erase(ring.begin());
    if (cls == 'L') {
      cl.vis++;   // the never-blinks rate's denominator: long frames with the direction on the image
      cl.ampL.push_back(on ? (uint8_t)std::min(255, pk) : 0);   // zero-filled by frame: a hit-run is not a beat
      if (cl.ampL.size() > 64) cl.ampL.erase(cl.ampL.begin());
    }
  }
}

// Stage by stage, every 10 s: what the cameras delivered, what the gates kept, what the cells and the
// candidate stage made of it. The funnel the whole detection rests on -- if a stage is empty the ones
// after it cannot be right, and this is what says which.
void Sync::LogStages(double now) {
  if (now - stages_t_ < 10) return;
  stages_t_ = now;
  static const char *cn[2] = {"L", "S"};
  std::string c;
  for (int k = 0; k < kMaxCam; k++) {
    for (int i = 0; i < 2; i++) {
      long fr = camsp_[k].fr[i] - camsp_last_[k].fr[i], dt = camsp_[k].dots[i] - camsp_last_[k].dots[i];
      long bl = camsp_[k].blobs[i] - camsp_last_[k].blobs[i];
      if (!fr && !dt && !bl) continue;
      c += Fmt("%sc%d%s %ldf %ldd %ldb pk%d fl~%ld max%d", c.empty() ? "" : " | ", k, cn[i], fr, dt, bl, camsp_[k].pk[i],
             dt ? camsp_[k].fl[i] / dt : 0, camsp_[k].flmax[i]);
    }
    camsp_last_[k] = camsp_[k];
  }
  if (!c.empty()) {
    log_("collect (frames/dots/blobs since last, peak): " + c);
    // RECORDED, so a replay's own numbers can be diffed against the driver's. The recording carries
    // every input, but a replay runs it on a different schedule -- the driver's periodic work fires on
    // a wall clock interleaved with arriving data -- and until both sides' accounts of the SAME
    // interval can be put side by side, that difference cannot be attributed. One line per ten seconds
    // against ~240 frame lines a second: nothing.
    Rec(now, "X collect %s", c.c_str());
  }
  Spots sp = spots();
  long seen = sp.seen - stages_sp_.seen, used = sp.used - stages_sp_.used;
  long cells = 0, tri = 0, fuzz = 0, cand = 0;
  {
    std::lock_guard<std::mutex> g(net_);
    cells = (long)cells_.size();
    for (auto &kv : cells_)
      if (kv.second.tri) { tri++; if (kv.second.fuzzy) fuzz++; }
    cand = (long)cands_.size();
  }
  if (seen || used) {
    // "dots" is the LONG frames' blobs -- the only ones that can be a ray. The short frames' dots are
    // counted separately (sp.sframe): they are evidence for the gate, never rays, so they must not sit
    // in the denominator as if they had been dropped.
    std::string sl = Fmt("stages: %ld long frames -> %ld dots -> %ld rays (%.0f%%) | dropped: still %ld, fast %ld, big %ld (kept big %ld), "
             "nopose %ld, dim %ld, S-gated %ld, fuzzy %ld, unmatched %ld | +%ld short-frame evidence dots | cells %ld (3D %ld, fuzzy %ld) -> candidates %ld",
             sp.frames - stages_sp_.frames, seen, used, seen ? 100.0 * used / seen : 0.0, sp.still - stages_sp_.still,
             sp.fast - stages_sp_.fast, sp.big - stages_sp_.big, sp.bigkeep - stages_sp_.bigkeep, sp.other - stages_sp_.other, sp.dim - stages_sp_.dim,
             sp.cellskip - stages_sp_.cellskip, sp.fuzz - stages_sp_.fuzz, sp.matchskip - stages_sp_.matchskip,
             sp.sframe - stages_sp_.sframe, cells, tri, fuzz, cand);
    log_(sl);
    // ... and RECORDED, so a replay's own account of the same interval can be put beside the driver's.
    // The recording carries every input; what it could not carry is the driver's schedule, and until
    // both sides account for the same ten seconds line by line that difference cannot be attributed.
    Rec(now, "X %s", sl.c_str());
  }
  stages_sp_ = sp;
}

// The recent frames, for the projection test. A candidate is a 3D PLACE, and a place is tested by
// projecting it -- not by keying a direction cell, which a real source rolls across as the head
// translates (a fixed cell key reads 16-37% where the place itself is lit on 80% of long frames).
void Sync::Remember(double t, int cam, char cls, const M3 &R, const V3 &p, const std::vector<Spot> &bl) {
  FRec f;
  f.t = t; f.cam = cam; f.cls = cls; f.R = R; f.p = p; f.bl = bl;
  frecs_.push_back(std::move(f));
  for (int guard = 0; guard < 64 && frecs_.size() > 8; guard++) {
    if (t - frecs_.front().t <= kPresWin) break;
    frecs_.pop_front();
  }
}

// Projection-based presence: for each candidate, the frames in which it was on a camera's image and
// whether a blob sat on the pixel it projects to, per exposure class. This is the honest version of
// the rings: no denominator wasted on frames where the direction was only loosely "in front", no
// marker sweeping off the cell as the head moves. MEASURED for a base station (qlhs_presence,
// 2026-10-11): 81-84% of the long frames it is in image for -- the design's 70-85% -- and 12-17% of
// the short ones, so the S verdict has to be the CONTRAST of the two, never the absolute short rate.
void Sync::PresenceAt(std::vector<Cand> &c, double t) {
  if (frecs_.empty()) return;
  for (Cand &cd : c) {
    for (const FRec &f : frecs_) {
      if (t - f.t > kPresWin) continue;
      M3 Rc; V3 tc;
      if (!optics_.Pose(f.cam, Rc, tc)) continue;
      V3 W = cd.q - (f.p + f.R * tc);
      if (norm(W) < 0.2) continue;
      double px, py;
      if (!optics_.Project(f.cam, T(Rc) * (T(f.R) * W), px, py)) continue;
      bool on = false;
      for (const Spot &b : f.bl) {
        double dx = b.x10 / 10.0 - 0.5 - px, dy = b.y10 / 10.0 - 0.5 - py;
        if (dx * dx + dy * dy <= kPresPx * kPresPx) { on = true; break; }
      }
      if (f.cls == 'L') { cd.pln++; cd.pl += on ? 1 : 0; } else { cd.psn++; cd.ps += on ? 1 : 0; }
    }
  }
}

// periodic: promote near-continuous fully-saturated directions to "lit"
// (lamp/window: present in ~every in-view frame, peak pinned at 255). A
// laser dot's dropouts keep its rate at ~0.7-0.8, well below the line -- and
// in a dim room its presence can peg anyway, so this only NAMES steady
// directions for the log; the short-frame gate below does the skipping.
// rotor-beat prominence over a candidate's presence ring (Goertzel, ~0.5-13 Hz): a station's
// in-view long-frame presence oscillates rotor-locked, a lamp's is flat, a glint's stochastic.
// EVIDENCE ACCUMULATION: the score sharpens as the ring grows. A score, not a gate (yet).
static double BeatScore(const std::vector<uint8_t> &ring, double fs) {
  int n = (int)ring.size();
  if (n < 32) return 0;
  double mean = 0;
  for (uint8_t v : ring) mean += v;
  mean /= n;
  static double win[64];
  for (int i = 0; i < n; i++) win[i] = (ring[i] - mean) * (0.5 - 0.5 * std::cos(2 * 3.14159265358979 * i / (n - 1)));
  double best = 0, sum = 0;
  for (int k = 0; k < 24; k++) {
    double f = (0.5 + 0.55 * k) / fs;
    double coef = 2 * std::cos(2 * 3.14159265358979 * f), s0 = 0, s1 = 0, s2 = 0;
    for (int i = 0; i < n; i++) { s0 = win[i] + coef * s1 - s2; s2 = s1; s1 = s0; }
    double m = std::sqrt(std::max(0.0, s0 * s0 + s1 * s1 - coef * s0 * s1));
    sum += m;
    if (m > best) best = m;
  }
  return sum > 0 ? best / (sum / 24) : 0;
}
// unit vector, guarding a degenerate input (a zero direction would make every angle NaN)
static V3 Unit(V3 v) { double n = norm(v); return n > 1e-9 ? v * (1.0 / n) : V3{0, 0, 1}; }
static double AngDeg(V3 a, V3 b) {
  return std::acos(std::max(-1.0, std::min(1.0, dot(Unit(a), Unit(b))))) * kDeg;
}
// Pearson correlation of two presence rings, aligned at the tail: two cells of one camera append
// a sample per in-view frame in lockstep, so their tails are the same frames. NAN when there is
// too little overlap or either ring is flat -- correlation is then UNDEFINED, not zero, and the
// caller falls to the parallax tie-break instead of blocking.
static double RingCorr(const std::vector<uint8_t> &a, const std::vector<uint8_t> &b) {
  int n = (int)std::min(a.size(), b.size());
  if (n < 24) return std::numeric_limits<double>::quiet_NaN();
  const uint8_t *A = a.data() + (a.size() - n), *B = b.data() + (b.size() - n);
  double ma = 0, mb = 0;
  for (int i = 0; i < n; i++) { ma += A[i]; mb += B[i]; }
  ma /= n; mb /= n;
  double sa = 0, sb = 0, sab = 0;
  for (int i = 0; i < n; i++) { double x = A[i] - ma, y = B[i] - mb; sa += x * x; sb += y * y; sab += x * y; }
  if (sa < 1e-9 || sb < 1e-9) return std::numeric_limits<double>::quiet_NaN();
  return sab / std::sqrt(sa * sb);
}
static double PointResidual(const V3 *O, const V3 *D, int n, V3 *qout = nullptr) {
  auto fit = [&](const int *idx, int k, double out[3]) {
    double m[9] = {0}, v[3] = {0};
    for (int i = 0; i < k; i++) {
      double o[3] = {O[idx[i]].x, O[idx[i]].y, O[idx[i]].z};
      double d[3] = {D[idx[i]].x, D[idx[i]].y, D[idx[i]].z};
      double doto = d[0]*o[0] + d[1]*o[1] + d[2]*o[2];
      for (int r = 0; r < 3; r++) {
        for (int c = 0; c < 3; c++) m[r*3+c] += (r == c) - d[r]*d[c];
        v[r] += o[r] - d[r]*doto;
      }
    }
    double cf[9];
    cf[0] =  m[4]*m[8] - m[5]*m[7];  cf[1] = -m[3]*m[8] + m[5]*m[6];  cf[2] =  m[3]*m[7] - m[4]*m[6];
    cf[3] = -m[1]*m[8] + m[2]*m[7];  cf[4] =  m[0]*m[8] - m[2]*m[6];  cf[5] = -m[0]*m[7] + m[1]*m[6];
    cf[6] =  m[1]*m[5] - m[2]*m[4];  cf[7] = -m[0]*m[5] + m[2]*m[3];  cf[8] =  m[0]*m[4] - m[1]*m[3];
    double det = m[0]*cf[0] + m[1]*cf[1] + m[2]*cf[2];
    if (std::fabs(det) < 1e-12) return false;
    for (int r = 0; r < 3; r++) out[r] = (cf[r*3]*v[0] + cf[r*3+1]*v[1] + cf[r*3+2]*v[2]) / det;
    return true;
  };
  int idx[16], keep[16], kn = 0;
  for (int i = 0; i < n; i++) idx[i] = i;
  double p[3];
  if (!fit(idx, n, p)) return 1e9f;
  double res[16], sorted[16];
  for (int i = 0; i < n; i++) {
    double wx = p[0] - O[i].x, wy = p[1] - O[i].y, wz = p[2] - O[i].z;
    double d[3] = {D[i].x, D[i].y, D[i].z};
    double along = wx*d[0] + wy*d[1] + wz*d[2];
    sorted[i] = res[i] = std::sqrt(std::max(0.0, wx*wx + wy*wy + wz*wz - along*along));
  }
  std::sort(sorted, sorted + n);
  double lim = sorted[n - 1 - n / 4];   // drop the worst quarter once: one outlier-resistant pass
  for (int i = 0; i < n; i++) if (res[i] <= lim && kn < 16) keep[kn++] = i;
  if (kn >= 4 && fit(keep, kn, p)) {
    double s = 0;    for (int i = 0; i < kn; i++) {
      double wx = p[0] - O[keep[i]].x, wy = p[1] - O[keep[i]].y, wz = p[2] - O[keep[i]].z;
      double d[3] = {D[keep[i]].x, D[keep[i]].y, D[keep[i]].z};
      double along = wx*d[0] + wy*d[1] + wz*d[2];
      s += std::max(0.0, wx*wx + wy*wy + wz*wz - along*along);
    }
    if (qout) *qout = V3{(float)p[0], (float)p[1], (float)p[2]};
    return std::sqrt(s / kn);
  }
  double s = 0;  for (int i = 0; i < n; i++) s += res[i] * res[i];
  return std::sqrt(s / n);
}

// The three-state verdict on a candidate, from the evidence it already carries. Only POSITIVE
// evidence rejects. See Sync::QKind for why a flat beat is not evidence and why the test is the
// 0.3 ms exposure: a source that appears in the short frames nearly every time they look is
// radiating continuously -- a lamp, a window, a glint -- because a beam sweeping past cannot.
static int QKindOf(int ps, int psn) {
  if (psn < Sync::kQKindN) return Sync::QK_UNKNOWN;      // too few short-frame looks to judge either way
  if (ps <= 0) return Sync::QK_UNKNOWN;                  // absent from the short frames: ambiguous, not evidence
  double r = (double)ps / psn;
  if (r >= Sync::kQMains) return Sync::QK_REJECT;
  if (r <= Sync::kQSpeck) return Sync::QK_SOURCE;
  return Sync::QK_UNKNOWN;                               // in between: refuse to judge
}
static double QWeight(int kind) { return kind == Sync::QK_SOURCE ? 1.0 : Sync::kQUnknownW; }
static const char *QName(int kind) {
  return kind == Sync::QK_SOURCE ? "source-like" : kind == Sync::QK_REJECT ? "CONTINUOUS" : "unknown";
}

void Sync::CellsSweep(double t) {
  long lit = 0, watched = 0; double best_rate = 0; long best_seen = 0, best_vis = 0;
  cands_.clear();   // rebuilt every sweep: the layout matcher reads this snapshot
  std::vector<Cand> raw;
  for (int cam = 0; cam < 4; cam++) actk_[cam].clear();
  for (auto it = cells_.begin(); it != cells_.end();) {
    SpotCell &cl = it->second;
    cl.matched = false;
    cl.mains = false;
    if (cl.last < t - 45) { it = cells_.erase(it); continue; }
    if (cl.last >= t - 10) { cl.inact = true; actk_[cl.cam].push_back(CellKey(cl.cam, cl.w)); }
    else cl.inact = false;
    ++it;
  }
  if (cells_.size() > 6000) {  // shed the coldest, not unordered_map's arbitrary front
    size_t shed = cells_.size() - 6000;
    std::vector<uint64_t> keys;
    keys.reserve(cells_.size());
    for (auto &kv : cells_) keys.push_back(kv.first);
    std::nth_element(keys.begin(), keys.begin() + shed, keys.end(),
                     [&](uint64_t a, uint64_t b) { return cells_.at(a).last < cells_.at(b).last; });
    keys.resize(shed);
    for (uint64_t k : keys) cells_.erase(k);
  }
  for (auto it = cells_.begin(); it != cells_.end();) {
    SpotCell &cl = it->second;
    if (cl.last < t - 120) { it = cells_.erase(it); continue; }
    if (cl.nO >= 6) {   // enough travel-spread samples: does this direction converge on ONE 3D point?
      V3 c{0, 0, 0};
      for (int i = 0; i < cl.nO; i++) c = c + cl.sO[i];
      c = c * (1.0f / cl.nO);
      double spread = 0;
      for (int i = 0; i < cl.nO; i++) spread = std::max(spread, norm(cl.sO[i] - c));
      cl.trav = spread >= kParallaxTravel;   // the head moved enough for the angular-offset (parallax) test to mean something
      if (spread >= 0.25) {   // enough head travel to triangulate at all
        cl.res3d = (float)PointResidual(cl.sO, cl.sD, cl.nO, &cl.q3);
        cl.fuzzy = cl.res3d > 0.08;   // a station aperture stays cm-tight; glow, edges and big lamps smear
        cl.tri = true;
      }
    }
    if (cl.tri && !cl.fuzzy && cl.nO >= 6 && cl.seen >= 20 && t - cl.last < 10) {
      int sh = 0;
      for (uint8_t v : cl.seqS) sh += v;
      Cand cd;
      cd.q = cl.q3; cd.res = cl.res3d; cd.seen = cl.seen; cd.key = it->first;
      cd.beat = (float)std::max(BeatScore(cl.ampL, 33.0), BeatScore(cl.seqL, 33.0));
      cd.beatm = std::max(cl.ampL.size(), cl.seqL.size()) >= 32;   // BeatScore needs 32 samples; below that its 0 is unknown
      cd.keys = {it->first};
      cd.sh = sh; cd.sf = (int)cl.seqS.size();
      cd.dir = Unit(cl.w); cd.dir0 = Unit(cl.w0); cd.cam = cl.cam; cd.trav = cl.trav; cd.ring = cl.seqL;
      cd.snpx = cl.snpx; cd.sn = cl.sn;
      raw.push_back(std::move(cd));
    }
    if (cl.vis >= 60 && cl.seen >= 16) {      watched++;
      double rate = (double)cl.seen / cl.vis;
      // lifetime peaks: a lamp's sightings saturate and stay saturated; a dot
      // dipping through beat nulls swings (and, in a dim room, its long-frame
      // presence can peg anyway -- which is why nothing here GATES, it only
      // names the steady directions; short-frame presence does the gating)
      bool tight = cl.amin >= 240 && cl.amax - cl.amin < 48;
      if (!cl.lit && rate >= 0.85 && tight) {
        cl.lit = true;
        log_(Fmt("cameras: a direction present %.0f%% of %ld in-view frames (peak %ld..%ld) never blinks: a lamp, window or other steady source",
                 rate * 100, cl.vis, (long)cl.amin, (long)cl.amax));
      }
      if (rate > best_rate && cl.vis >= 120) { best_rate = rate; best_seen = cl.seen; best_vis = cl.vis; }
      if (cl.lit) lit++;
    }
    ++it;
  }
  PresenceAt(raw, t);   // the projection test, before the merge: the arms of one place pool their evidence
  // merge starburst fragments IN ANGLE: a starburst is an optical effect, so its arms sit at a
  // roughly fixed ANGULAR offset from the core (a few degrees, whatever the range) and modulate
  // together -- a fixed 3D radius is a different angle at every range (17 deg at 1 m, 2 deg at
  // 8 m) and arm 3D points have poor parallax anyway. Greedy, core-first (the most-seen direction
  // is the core). A pair merges when the directions are within kMergeDeg AND the evidence says one
  // source: same camera + correlated presence rings; never across contrasting signatures (one
  // rotor beat, one flat -- a station beside a lamp); 3D distance is a VETO, not the criterion
  // (both well-conditioned and clearly far apart = two sources); and when the correlation is
  // unavailable (different camera, or a flat ring) the parallax tie-break decides: a starburst's
  // angular offset holds as the head travels, two real sources at different depths drift apart.
  // The merged candidate carries every member key (all arms gate together) and the best beat.
  std::sort(raw.begin(), raw.end(), [](const Cand &a, const Cand &b) { return a.seen > b.seen; });
  int mlog = 0, rlog = 0;
  for (const Cand &c : raw) {
    bool merged = false;
    for (Cand &m : cands_) {
      double sep = AngDeg(c.dir, m.dir);
      if (sep > kMergeDeg) continue;
      double r = (c.cam >= 0 && c.cam == m.cam) ? RingCorr(c.ring, m.ring) : std::numeric_limits<double>::quiet_NaN();
      double bh = std::max(c.beat, m.beat), bl = std::min(c.beat, m.beat);
      double d3 = norm(c.q - m.q);
      double drift = std::fabs(sep - AngDeg(c.dir0, m.dir0));
      const char *why = nullptr;
      // the signature guard needs BOTH beats MEASURED: a short ring scores 0, and reading that as a
      // flat lamp blocked a station from merging with its own fragments (beats 0.00/7.89 at 0.5 deg).
      if (c.beatm && m.beatm && bh >= kBeatHas && bl < kBeatFlat) why = "contrasting signatures (one beats, one is flat)";
      else if (std::isfinite(r) && r < kMergeCorr) why = "presence rings not correlated";
      else if (c.res <= kCondRes && m.res <= kCondRes && d3 > kMergeFar) why = "both well-conditioned and far apart in 3D";
      else if (!std::isfinite(r) && c.trav && m.trav && drift > kDriftDeg) why = "angular offset drifts with head travel (parallax)";
      std::string rs = std::isfinite(r) ? Fmt("%.2f", r) : std::string("n/a");
      std::string dr = (c.trav && m.trav) ? Fmt(", drift %.2f deg", drift) : std::string();
      if (why) {
        if (rlog < 6) {
          rlog++;
          log_(Fmt("merge: held q[%.1f %.1f %.1f] apart from q[%.1f %.1f %.1f]: %.1f deg, corr %s, beats %.2f/%.2f, 3D %.2f m%s -- %s",
                   c.q.x, c.q.y, c.q.z, m.q.x, m.q.y, m.q.z, sep, rs.c_str(), c.beat, m.beat, d3, dr.c_str(), why));
        }
        continue;
      }
      if (mlog < 10) {
        mlog++;
        log_(Fmt("merge: q[%.1f %.1f %.1f] <- q[%.1f %.1f %.1f]: %.1f deg, corr %s, beats %.2f/%.2f, 3D %.2f m%s",
                 m.q.x, m.q.y, m.q.z, c.q.x, c.q.y, c.q.z, sep, rs.c_str(), c.beat, m.beat, d3, dr.c_str()));
      }
      m.keys.push_back(c.keys[0]);
      if (c.beat > m.beat) m.beat = c.beat;
      m.beatm = m.beatm || c.beatm;
      if (c.seen > m.seen) m.seen = c.seen;
      m.snpx = std::max(m.snpx, c.snpx); m.sn += c.sn;
      m.sh += c.sh;
      m.sf += c.sf;
      m.pl += c.pl; m.pln += c.pln; m.ps += c.ps; m.psn += c.psn;
      merged = true;
      break;
    }
    if (!merged) cands_.push_back(c);
  }
  // the S-frame verdict, at CANDIDATE level: a source's S-evidence pools across its merged
  // starburst members, so a lamp dies WHOLE -- before, arms with thin rings could keep a gated
  // core's family alive. Mains candidates are removed from the matcher's list entirely and their
  // member cells carry the verdict to the ray gate in Use.
  //
  // PARKED: the verdict is computed and logged, but not applied. The gate's evidence is the very
  // ring that was broken until now, and the ground truth (qlhs_presence) says a real base station
  // is NOT absent from the short frames -- it reads 5-7% (LHB-904ECE4F 4.9, LHB-EC8BB104 6.8) where
  // this threshold asks for 15%. That is a 2x margin, not the "0.00 against 0.15-0.5" the rule was
  // written for, so it earns its place only once the rings read correctly and the margin has been
  // measured on more rooms. Flip kApplyMainsGate when it has.
  constexpr bool kApplyMainsGate = false;
  int dropped_mains = 0, would_mains = 0;
  for (auto ci = cands_.begin(); ci != cands_.end();) {
    bool mains = ci->sf >= 8 && (double)ci->sh / ci->sf >= 0.15;
    if (mains) would_mains++;
    if (!kApplyMainsGate) { ci->mains = mains; ++ci; continue; }
    for (uint64_t k : ci->keys) cells_[k].mains = mains;
    if (mains) { ci = cands_.erase(ci); dropped_mains++; }
    else ++ci;
  }
  // the S-gate's own activity: it is only live if the class rings are real (the tag parse bug left seqS empty of 1s,
  // so nothing was ever dropped and the gate read as "working" while doing nothing). Parked, it still says what it
  // would have done -- that is the number to watch while the rings are validated.
  if (dropped_mains)
    log_(Fmt("S-gate: %d candidate(s) dropped as steady sources (seen on 8+ short frames, %.0f%%+ short-presence)",
             dropped_mains, 15.0));
  else if (kApplyMainsGate && would_mains)
    log_(Fmt("S-gate: %d candidate(s) carry the steady-source verdict (8+ short frames, %.0f%%+ short-presence)", would_mains,
             15.0));
  // NAMING, diagnostic only (no decision hangs on it yet). A 2.0 base station's rotors sweep at a fixed period, so
  // the frames that see one station's dot modulate at ITS channel; matching a candidate's sighting train to the known
  // channel periods names it. Upstream names stations only after a lock (CheckIdentity); naming the CANDIDATE is what
  // would let acquisition tell one station from the other -- i.e. kill a mirror without the worn or held devices. The
  // train is built the way upstream's CheckLayout builds it: the frames a direction was in view of, and which of them
  // saw it (a sighting whose world direction keys to this cell). Everything is accumulated from what has arrived.
  if (!cands_.empty()) {
    std::map<std::string, int> ch;
    {
      std::lock_guard<std::mutex> g(raw_m_);
      ch = channels_;
    }
    std::vector<double> tau;
    std::vector<std::string> cn;
    for (auto &kv : ch) {
      double p = RotorPeriod(kv.second);
      if (p > 0) { tau.push_back(p); cn.push_back(kv.first); }
    }
    if (!tau.empty()) {
      // the train is built about the candidate's PLACE, not its direction cell (upstream's CheckLayout): a train keyed
      // to a 2.4 deg cell is chopped up by the head sweeping the dot across cell borders, which buries the modulation.
      // A frame sees the place if it lands in that camera's image facing it; a sighting counts for it if the ray passes
      // within kNameNear of it, and its frame is the one that saw the dot.
      std::vector<double> G;
      std::vector<int> C;
      std::vector<V3> O, D;
      solver_.SightingsSince(t - 120.0, G, C, O, D);
      std::vector<int> fc;
      std::vector<double> fg;
      std::vector<V3> fo;
      std::vector<M3> fR;
      solver_.FramesSince(t - 120.0, fc, fg, fo, fR);
      const double kNameNear = 1.0;   // deg: a ray this close to the place saw it
      if (!(grid_.period(cands_[0].cam) > 0))
        log_(Fmt("candidate naming: no frame period learned for cam %d yet: a train has no clock to sit on", cands_[0].cam));
      for (size_t i = 0; i < cands_.size() && i < 8; i++) {
        const Cand &cd = cands_[i];
        double P = grid_.period(cd.cam);
        if (!(P > 0)) continue;
        std::set<std::pair<int, double>> sw;
        for (size_t n = 0; n < G.size(); n++) {
          V3 w = cd.q - O[n];
          double d = norm(w);
          if (d < 0.4) continue;
          if (std::acos(std::max(-1.0, std::min(1.0, dot(w, D[n]) / d))) * kDeg < kNameNear) sw.insert({C[n], G[n]});
        }
        std::vector<InView> v;
        for (size_t n = 0; n < fc.size(); n++) {
          if (fc[n] != cd.cam || !std::isfinite(fg[n])) continue;
          V3 U = cd.q - fo[n];
          double d = norm(U);
          if (!(d > 0.3) || !(dot(U, fR[n] * V3{0, 0, 1}) >= 0.2 * d)) continue;   // facing it
          // net_ is already held: CellsSweep runs from Use, which OnLine calls under it (do NOT re-lock)
          if (!(have_optics_ && optics_.InImage(fc[n], T(fR[n]) * (U * (1 / d)), 8.0))) continue;
          v.push_back({fc[n], fg[n], sw.count({fc[n], fg[n]}) > 0});
        }
        int seen = 0;
        for (auto &e : v) seen += e.seen;
        std::vector<double> sc;
        if (!SeenScores(v, tau, P, sc)) {
          log_(Fmt("candidate naming: q[%.1f %.1f %.1f] seen %ld: not enough of a train (%d frames, %d seen, period %.1f ms)",
                   cd.q.x, cd.q.y, cd.q.z, cd.seen, (int)v.size(), seen, P * 1000));
          continue;
        }
        int b = (int)(std::max_element(sc.begin(), sc.end()) - sc.begin());
        double other = 0;
        for (size_t j = 0; j < sc.size(); j++)
          if ((int)j != b) other = std::max(other, sc[j]);
        bool clear = sc[b] - other >= kNamed;
        int l1 = 0;
        for (uint8_t e : cd.ring) l1 += e;
        log_(Fmt("candidate naming: q[%.1f %.1f %.1f] seen %ld: %s (%.3f, next %.3f) on %d frames, %d seen; L %d/%zu S %d/%d; "
                 "proj L %d/%d S %d/%d; short footprint %.0f px on %d; %s (%d/%d short)",
                 cd.q.x, cd.q.y, cd.q.z, cd.seen, clear ? cn[b].c_str() : "no channel", sc[b], other, (int)v.size(), seen,
                 l1, cd.ring.size(), cd.sh, cd.sf, cd.pl, cd.pln, cd.ps, cd.psn, cd.snpx, cd.sn, QName(cd.kind), cd.ps, cd.psn));
      }
    }
  }
  if (watched > 0 && best_seen > 0)
    log_(Fmt("cameras: %ld presence checks, highest %.0f%% (%ld/%ld), %ld judged steady -- a real station sits near 70-85%%",
             watched, best_rate * 100, best_seen, best_vis, lit));
  else
    log_(Fmt("cameras: %zu direction cells tracked, none old enough to judge yet", cells_.size()));
}

// layout match every sweep: accepted candidates' cells become the only ray sources (Use gates on
// cl.matched); three sweeps without an accepted correspondence means x is fitted onto decoys --
// drop the sightings and re-acquire rather than keep a lock the arrangement says is wrong.
// acquisition enumeration: candidate pairs -> closed-form 4-DOF per assignment (YawFit maps the
// pair onto the SteamVR pair) -> support score over recent rays. Accept only an unambiguous
// winner: real support and twice the runner-up. The rotor-beat score softly tips ties between a
// decoy pair and the real pair (stations beat, lamps and reflections do not). Ambiguity holds:
// x stays unset and the sweep comes back around.
void Sync::AcquireFromCandidates(double t) {
  if ((int)cands_.size() < 2) return;
  std::vector<std::string> keys;
  std::vector<V3> S, Z;
  solver_.Stations(keys, S, Z);
  if (S.size() != 2) return;
  // LOOSE FRONT, STRICT BACK.
  //
  // The baseline is a hard fact: the two stations sit a measured distance apart in the Lighthouse
  // frame (SteamVR's own S lines) and x is rigid, so a valid pair of candidate places must be that
  // same distance apart in SLAM space. But holding the RATCHET at both ends -- only considering pairs
  // within 1 cm, then 2, ... -- meant the enumeration sat out the start of every session waiting for
  // the camera-measured places to become accurate enough to enter the band at all: it took three
  // minutes to lock on the 2026-10-13 walk, and that is the cost being paid here.
  //
  // So generation is generous and verification is not. Any pair within kBaseGenTol is a candidate for
  // SCORING (the nearest-baseline kBaseGenMaxPairs of them, to bound the cost of ScoreX over
  // thousands of rays), and the winner is then VERIFIED against kBaseAccept before anything is
  // installed. Being generous here buys time, not confidence: every discriminator that decides is
  // downstream (the baseline check, the worn/held devices, the runner-up rule), and this way they get
  // something to decide on minutes sooner. A fit that only exists outside kBaseAccept is logged and
  // refused, not installed.
  // 5 cm, not 15: a 15 cm window let the highest-SUPPORT pair be one whose separation was 14 cm off
  // the baseline -- geometrically implausible on its face -- and made the ranking argue about pairs
  // that could never be installed anyway. 5 cm is the window that is physically meaningful here, so
  // the ranking only ever sees pairs that could actually win.
  const double kBaseGenTol = 0.05, kBaseAccept = 0.05;
  const size_t kBaseGenMaxPairs = 400;
  const double base_s = norm(S[1] - S[0]);
  struct Pair { int i, j; double e; };   // e = how far this pair's separation is from the baseline
  std::vector<Pair> pairs;
  int rej_pairs = 0, far_pairs = 0;
  const int n = (int)cands_.size();
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i == j) continue;
      V3 vq = cands_[j].q - cands_[i].q, vs = S[1] - S[0];
      if (std::hypot(vq.x, vq.z) < 0.5 || std::hypot(vs.x, vs.z) < 0.5) continue;
      // EVIDENCE GATE, deliberately soft: a pair with POSITIVE evidence of continuous emission on
      // either side is not scored at all. Everything else is scored, and source-like candidates
      // outrank unknown ones, so a wrong gate degrades the ORDER of the pairs rather than making
      // acquisition impossible, and an unknown candidate is never punished for silence.
      if (cands_[i].kind == QK_REJECT || cands_[j].kind == QK_REJECT) { rej_pairs++; continue; }
      double e = std::fabs(norm(vq) - base_s);
      if (e > kBaseGenTol) { far_pairs++; continue; }
      pairs.push_back({i, j, e});
    }
  }
  // nearest the baseline first, so if the cap bites we keep the most plausible pairs, not an
  // arbitrary slice of them
  std::sort(pairs.begin(), pairs.end(), [](const Pair &a, const Pair &b) { return a.e < b.e; });
  int gen_dropped = 0;
  if (pairs.size() > kBaseGenMaxPairs) { gen_dropped = (int)(pairs.size() - kBaseGenMaxPairs); pairs.resize(kBaseGenMaxPairs); }

  // YawFit(A,B) returns the transform mapping A -> B. x maps the LIGHTHOUSE universe onto the
  // head/SLAM one (x(S) = Ry(yaw)*S + t), so the stations go first and the candidates second.
  // Reversed, x is the inverse map and ScoreX measures nothing (measured: 0-2 support over ~16k
  // rays). ScoreX's gate is in DEGREES (Nearest returns acos*kDeg); 1.0/kDeg made it 0.0175 deg.
  std::vector<X4> px(pairs.size());
  std::vector<int> psc(pairs.size());
  for (size_t k = 0; k < pairs.size(); k++) {
    px[k] = YawFit({S[0], S[1]}, {cands_[pairs[k].i].q, cands_[pairs[k].j].q});
    psc[k] = solver_.ScoreX(px[k], 1.0);
  }
  int bi = -1, bsc = 0, ssc = 0, ri = -1, rj = -1;
  double bw = -1;
  X4 xb{}, xr{};
  for (size_t k = 0; k < pairs.size(); k++) {
    // the beat only BREAKS NEAR-TIES (as a multiplicative weight it once let a high-beat pair win
    // over one with more support, and put the runner-up above the best)
    double w = QWeight(cands_[pairs[k].i].kind) * QWeight(cands_[pairs[k].j].kind) * psc[k] +
               0.05 * std::min(cands_[pairs[k].i].beat, cands_[pairs[k].j].beat);
    if (w > bw) { bw = w; bsc = psc[k]; bi = (int)k; xb = px[k]; }
  }
  const double kSameX = 1.0;   // m: two fits predicting the stations within this are ONE solution
  if (bi >= 0) {
    for (size_t k = 0; k < pairs.size(); k++) {
      if (Apart(px[k], xb, S) <= kSameX) continue;   // the same answer, not a rival
      if (psc[k] > ssc) { ssc = psc[k]; ri = pairs[k].i; rj = pairs[k].j; xr = px[k]; }
    }
  }
  // ---- verification. Nothing is installed on generation evidence alone. ----
  bool accept = false;
  int wi = -1, wj = -1, wsc = bsc;
  X4 x{};
  std::string why;
  double band = bi >= 0 ? pairs[bi].e : 0;
  if (bi >= 0 && band > kBaseAccept) {
    why = Fmt("the best pair is %.1f cm off the baseline, past the %.0f cm a fit must hold to be installed",
              band * 100, kBaseAccept * 100);
  } else if (bi >= 0 && bsc >= 20) {
    // A fit that needs the WITNESS is a close call by definition: the 2x rule could not separate it, so
    // a motion reading decides -- and that reading has the same shape for a right fit and a wrong one.
    // Measured: +0.84/-0.84 on a user-confirmed-correct lock (support 767) and +0.82/-0.82 on a wrong
    // one four minutes later (support 80). What DID separate them is how much light agrees: correct
    // locks scored 628 and 767 rays inside the gate, every wrong one 80-200. So the witness may only
    // break a tie that is worth breaking, while a low-support fit may still lock when it is UNAMBIGUOUS
    // on the 2x rule alone -- the evidence then decides, not the devices.
    const int kStrongSupport = 300;
    if (ssc == 0 || ssc < 20 || bsc >= 2 * ssc) {
      accept = true; wi = pairs[bi].i; wj = pairs[bi].j; x = xb; wsc = bsc;
    } else if (ri >= 0 && bsc < kStrongSupport) {
      why = Fmt("too close to call (%d against %d) and only %d rays agree: not enough light to break the "
                "tie with the worn or held devices (needs %d at this kind of margin)", bsc, ssc, bsc, kStrongSupport);
    } else if (ri >= 0) {
      // Two real solutions -- the mirror case. Two stations facing each other fit as well turned
      // 180 deg about their midpoint, and the ORIENTATIONS cannot break the tie (measured: the mirror
      // predicts axes 7.9 deg from the truth at both stations). The devices the user wears or holds
      // stay near the head, so they decide. A marginal reading must never INSTALL the runner-up.
      double da, db;
      std::string how;
      int fav = solver_.WitnessFavours(xb, xr, da, db, how);
      if (fav == 1) { accept = true; wi = pairs[bi].i; wj = pairs[bi].j; x = xb; wsc = bsc; why = how; }
      else if (fav == 2) why = Fmt("the worn or held devices favour the rival (%s)", how.c_str());
      else why = "two real solutions, no witness";
    }
  }
  if (t > cand_log_) {   // shares the 60 s evidence-watch throttle
    cand_log_ = t + 60;
    std::string s2 = ri >= 0 ? Fmt(" [rival %d at yaw %+.1f, apart %.1f m]", ssc, xr[0] * kDeg, Apart(xr, xb, S)) : std::string();
    log_(Fmt("acquisition: %d pairs generated within %.0f cm of the baseline (%.0f-%.1f cm of separation error, "
             "%d skipped: continuous evidence, %d beyond the band, %d dropped over the pair cap), best support %d "
             "at %.1f cm%s%s", (int)pairs.size(), kBaseGenTol * 100,
             pairs.empty() ? 0.0 : pairs.front().e * 100, pairs.empty() ? 0.0 : pairs.back().e * 100,
             rej_pairs, far_pairs, gen_dropped, wsc, band * 100, s2.c_str(),
             accept ? Fmt(" -- accepted%s", why.empty() ? "" : Fmt(" (%s)", why.c_str()).c_str()).c_str()
                    : Fmt(" -- holding: %s", why.empty() ? "ambiguous or unsupported" : why.c_str()).c_str()));
  }
  if (!accept) return;
  {   // What would the enumeration have preferred WITHOUT the evidence gate, at the same precision?
      // Removing clutter pairs also shrinks the runner-up, which flatters the 2x rule: the fit would
      // then be installed on apparent confidence that the gate created. Computed only on the accept
      // path, so the cost is paid once per lock, and said out loud when the gate was decisive.
      int ub = 0, ubi = -1, ubj = -1, usr = 0, usri = -1, usrj = -1;
      for (const Pair &p : pairs) {
        if (p.e > kBaseAccept) continue;
        X4 xu = YawFit({S[0], S[1]}, {cands_[p.i].q, cands_[p.j].q});
        int sc = solver_.ScoreX(xu, 1.0);
        if (sc > ub) { ub = sc; ubi = p.i; ubj = p.j; }
        if (Apart(xu, xb, S) > kSameX && sc > usr) { usr = sc; usri = p.i; usrj = p.j; }
      }
      bool decisive = ub > wsc || usr > ssc;
      log_(Fmt("the evidence gate: gated best %d, gated rival %d | UNGATED best %d%s, ungated rival %d%s%s",
               wsc, ssc, ub, ubi >= 0 ? Fmt(" (candidates %d,%d %s/%s)", ubi, ubj, QName(cands_[ubi].kind), QName(cands_[ubj].kind)).c_str() : "",
               usr, usri >= 0 ? Fmt(" (candidates %d,%d %s/%s)", usri, usrj, QName(cands_[usri].kind), QName(cands_[usrj].kind)).c_str() : "",
               decisive ? " -- the gate was DECISIVE, so it carries the blame if this fit is wrong" : " -- the gate changed nothing here"));
  }
  solver_.Reset(x);
  solver_.SetLockTime(t);
  log_(Fmt("acquired by enumeration: candidates %d,%d (%s/%s) support %d (rival %d, beats %.2f/%.2f) "
           "yaw %+.2f deg t [%.2f %.2f %.2f], baseline %.1f cm of %.3f m (verified, <= %.0f cm)%s",
           wi, wj, QName(cands_[wi].kind), QName(cands_[wj].kind), wsc, ssc, cands_[wi].beat, cands_[wj].beat,
           x[0] * 180 / kPi, x[1], x[2], x[3], band * 100, base_s, kBaseAccept * 100,
           why.empty() ? "" : Fmt(" (%s)", why.c_str()).c_str()));
}

void Sync::MatchSweep(double t) {
  if (!solver_.has_x() || (int)cands_.size() < 2) {
    if (!solver_.has_x()) { match_fail_ = 0; AcquireFromCandidates(t); }
    return;
  }
  std::vector<std::string> keys;
  std::vector<V3> S, Z, P, Zq;
  solver_.Stations(keys, S, Z);
  Solver::Predict(solver_.x(), S, Z, P, Zq);
  MatchReport mr = MatchLayout(P.data(), (int)P.size());
  if (mr.ok) {
    match_fail_ = 0;
    for (size_t k = 0; k < keys.size() && k < 2; k++) {
      int ci = k == 0 ? mr.a0 : mr.a1;
      if (ci >= 0 && ci < (int)cands_.size()) cells_[cands_[ci].key].matched = true;
    }
  } else if (++match_fail_ == 3) {
    log_("layout match: 3 sweeps without an accepted correspondence -- the fit is not on the stations, re-acquiring");
    solver_.PoseBreak(t);
  }
  if (t > cand_log_) {   // evidence watch: per-candidate beat prominence, the lighthouse-vs-light signal
    cand_log_ = t + 60;
    std::string s;
    for (size_t i = 0; i < cands_.size() && i < 6; i++)
      s += Fmt("  q[%.1f %.1f %.1f] beat %.2f seen %ld%s", cands_[i].q.x, cands_[i].q.y, cands_[i].q.z,
               cands_[i].beat, cands_[i].seen, cands_[i].keys.size() > 1 ? Fmt(" (%d arms)", (int)cands_[i].keys.size() - 1).c_str() : "");
    int nsrc = 0, nunk = 0, nrej = 0;
    for (const Cand &c : cands_) { if (c.kind == QK_SOURCE) nsrc++; else if (c.kind == QK_REJECT) nrej++; else nunk++; }
    log_(Fmt("candidates: %zu (%d source-like, %d unknown, %d continuous/rejected).%s", cands_.size(), nsrc, nunk, nrej, s.c_str()));  }
}

void Sync::Use(const Shot &f) {
  const double t = f.t, hg = f.hg, grid_pc = f.grid_pc;
  const int cam = f.cam, nb = (int)f.bl.size();
  const bool have_lag = f.have_lag;
  const std::vector<Spot> &bl = f.bl;
  M3 R; V3 p;
  double w;
  const bool have_pose = poses_.At(t, R, p);
  if (f.cls == 'S') {   // evidence only: the short frames mark the class rings and carry no rays
    spots_.sframe += nb;
    if (have_pose) {
      if (optics_.exact_time()) {
        MarkPresence(cam, R, t, 'S', bl);
        // ... and their FOOTPRINT per direction. In a ~0.3 ms exposure a source that radiates
        // continuously is still extended (a lamp, a window) while a laser's speck is a dot, so this
        // is the one size test that the long exposure's bloom cannot corrupt. Never creates a cell:
        // a short-frame flash on its own is not a light we have located.
        for (const Spot &b : bl) {
          V3 o, d;
          if (!optics_.Ray(cam, b.x10 / 10.0 - 0.5, b.y10 / 10.0 - 0.5, o, d)) continue;
          auto it = cells_.find(CellKey(cam, R * d));
          if (it == cells_.end()) continue;
          SpotCell &cl = it->second;
          cl.sn++;
          cl.snpx = b.npx > cl.snpx ? (float)b.npx : cl.snpx * 0.95f;  // decaying max: stop being extended, forget it
        }
      }
      Remember(t, cam, 'S', R, p, bl);
    }
    return;
  }
  if (!have_pose || !poses_.Speed(t, w)) { spots_.other += nb; return; }
  bool still = poses_.Still(t);
  if (optics_.exact_time() && t > cells_next_) { cells_next_ = t + 30; CellsSweep(t); MatchSweep(t); }  // Quest: nothing to sweep or match
  // presence, with the frame's own pose, for both classes (see MarkPresence)
  if (optics_.exact_time()) MarkPresence(cam, R, t, 'L', bl);
  Remember(t, cam, 'L', R, p, bl);
  {  // while SteamVR's headset stands still, do the lamps and windows the cameras see stand still too?
    std::vector<std::pair<double, double>> big;
    for (const Spot &b : bl)
      if (b.npx >= kLampPx) big.push_back({b.x10 / 10.0, b.y10 / 10.0});
    auto &prev = prev_big_[cam];
    if (!still) img_moved_ = img_still_ = 0;
    else if (!prev.empty())
      for (auto &q : big) {
        double d = 1e9;
        for (auto &o : prev) d = std::min(d, std::hypot(q.first - o.first, q.second - o.second));
        img_still_ += d <= 3;
        img_moved_ += d > 10;
      }
    prev = big;
  }
  if (still) { spots_.still += nb; return; }
  M3 Rc;
  V3 tc;
  if (have_lag && w <= kWmax && optics_.Pose(cam, Rc, tc)) solver_.AddFrame(t, hg, cam, p + R * tc, R * Rc);
  for (const Spot &b : bl) {
    int x10 = b.x10, y10 = b.y10, npx = b.npx;
    V3 o, d;
    if (!optics_.Ray(cam, x10 / 10.0 - 0.5, y10 / 10.0 - 0.5, o, d)) { spots_.other++; continue; }
    if (npx > kMaxPx) {
      // A long exposure blooms a source to hundreds or thousands of pixels (measured 50-150x against
      // the short frame, at the same pixel), so size on THIS class does not tell a station from a
      // lamp -- and dropping on it deleted the sources. Keep the blob, and veto only when the SHORT
      // frame shows the same direction extended too, i.e. it radiates continuously. The rotor beat
      // (ampL) and the S/L contrast then carry the lamp-versus-station decision, as intended.
      if (!optics_.exact_time()) { spots_.big++; continue; }  // Quest: the size filter still stands
      auto it = cells_.find(CellKey(cam, R * d));
      if (it != cells_.end() && it->second.sn >= kShortN && it->second.snpx >= kShortBig) { spots_.big++; continue; }
      spots_.bigkeep++;
    }
    if (g_ray_dump) {
      V3 O = p + R * o, D = R * d;
      fprintf(g_ray_dump, "R %.6f %d %.6f %.6f %.6f %.7f %.7f %.7f %.6f %.6f %.6f %.1f %d\n", t, cam, O.x, O.y, O.z, D.x,
              D.y, D.z, p.x, p.y, p.z, w, b.peak);
    }
    if (optics_.exact_time()) {  // Steam Frame only; Quest keeps stock behavior
      V3 W = R * d;
      CellSeen(cam, W, (float)b.peak, t);      SpotCell &cl = cells_[CellKey(cam, W)];
      {   // spread ray samples across head travel for the 3D point test (fuzzy gate below, computed in CellsSweep)
        V3 O = p + R * o;
        double dx = O.x - cl.lastO.x, dy = O.y - cl.lastO.y, dz = O.z - cl.lastO.z;
        if (cl.nO == 0 || dx*dx + dy*dy + dz*dz >= 0.0036) {
          cl.lastO = O;
          int s = cl.nO < 12 ? cl.nO++ : (cl.slot_++ % 12);
          cl.sO[s] = O;
          cl.sD[s] = R * d;
        }
      }
      // PRIMARY filter: short-frame presence. The lighthouse's sweeps are
      // microsecond-wide and unsynced with the S exposure, so a real dot is
      // ~absent from S frames (measured 0.00 over thousands of frames). A
      // continuously-radiating source (lamp, window, steady IR) shows up on
      // S too: skip the direction. Measured on the lhscan captures: direct
      // lamp images sit at 0.15-0.5 S-presence (some PWM lamps lower -- those
      // are left for kMaxPx and the solver's geometry), and the two
      // steady-IR sources in the lighthouse-only capture land at 0.17-0.22.
      //
      // The L-presence continuity test is deliberately GONE: in a dim room
      // the long exposure stretches past the sweep period, every long frame
      // catches a flash, the laser's presence pegs at 1.0 and ANY
      // presence-threshold on the L ring false-negatives the real station.
      // Short-frame presence has the clean measured margin, so it alone
      // gates -- at CANDIDATE level (CellsSweep): a source's S-evidence pools
      // across its merged starburst members, so a lamp dies whole instead of
      // surviving through arms with thin rings. The verdict rides the cell.
      if (cl.mains) { spots_.cellskip++; continue; }
      if (cl.fuzzy) { spots_.fuzz++; continue; }   // the direction's rays do not converge on one 3D point: glow, an edge, a lamp face
      if (match_ok_ && !cl.matched) { spots_.matchskip++; continue; }   // the layout matcher accepted other directions: these rays do not reach the solve
    }
    if (have_lag && cfg_.learn_timing) timing_.Record(grid_pc, o, d, hg, cam, b.peak >= kBright);  // gated rays only on the Frame: lamps must not steer exposure learning; fast turns still feed it (Quest stock order)
    if (w > kWmax) { spots_.fast++; continue; }
    bool bright = b.peak >= kBright;
    solver_.Add(t, p + R * o, R * d, hg, cam, bright);  // dim ones too: they tell when the room is too light
    if (bright || solver_.starved()) spots_.used++;
    else spots_.dim++;
  }
}

std::string Sync::Describe(const Spots &a, const Spots &b) {
  long n = b.seen - a.seen;
  if (!n) return Fmt("no bright spots on %ld camera frames", b.frames - a.frames);
  std::string why;
  auto part = [&](long k, const char *what) { if (k > 0) why += Fmt("%s%ld %s", why.empty() ? "" : ", ", k, what); };
  part(b.still - a.still, "while the headset was still");
  part(b.fast - a.fast, "while it turned fast");
  part(b.big - a.big, "lamps or windows (extended on the short frames too)");
  part(b.bigkeep - a.bigkeep, "kept though big: compact on the short frames, so a source");
  part(b.other - a.other, "without a headset pose");
  part(b.dim - a.dim, "dim, kept for a room too light for the dots to saturate");
  part(b.cellskip - a.cellskip, "directions also lit on short frames (lamps, steady IR)");
  part(b.fuzz - a.fuzz, "directions that fail the 3D point test (glow, edges)");
  part(b.matchskip - a.matchskip, "directions the layout matcher did not accept");  part(b.sframe - a.sframe, "short-frame evidence dots");
  return Fmt("%ld bright spots, %ld used", n, b.used - a.used) + (why.empty() ? "" : " (" + why + ")");
}

void Sync::SetStationsRaw(const std::map<std::string, std::pair<V3, M3>> &raw) {
  std::map<std::string, std::pair<V3, M3>> ok;  // a rotation, finite: anything else would turn the frame to NaN
  for (auto &kv : raw) {
    const M3 &R = kv.second.second;
    V3 p = kv.second.first;
    bool good = std::isfinite(p.x + p.y + p.z);
    for (int c = 0; c < 3 && good; c++) {
      double n = R.m[0][c] * R.m[0][c] + R.m[1][c] * R.m[1][c] + R.m[2][c] * R.m[2][c];
      good = std::isfinite(n) && std::fabs(n - 1) < 0.01;
    }
    if (good) ok.insert(kv);
  }
  std::lock_guard<std::mutex> g(raw_m_);
  raw_ = ok;
}

void Sync::SetChannels(const std::map<std::string, int> &ch) {
  std::lock_guard<std::mutex> g(raw_m_);
  channels_ = ch;
}

// One or two base stations leave the level to SteamVR, and the cameras can't fix it (gravity.h); with three, an
// automatic frame is levelled by the cameras (LevelStep)
bool Sync::GravityFrame(M3 &C, M3 &turn, std::string &key) {
  if (!sfile_.loaded || frame_.ok() != 1) return false;
  size_t n;
  {
    std::lock_guard<std::mutex> g(raw_m_);
    n = raw_.size();
  }
  if (sfile_.autogen && n >= (size_t)Solver::LEVEL_STATIONS) return false;
  C = frame_.Anchored();
  V3 t;
  frame_.Rotation(turn, t);
  turn = turn * T(C);
  Quat q = ToQuat(sfile_.anchor_R);
  key = Fmt("%s %.5f %.5f %.5f %.5f", sfile_.anchor.c_str(), q.w, q.x, q.y, q.z);
  return true;
}

// An automatic reference frame is SteamVR's frame as it was when the anchor was first seen, tilt and all. The headset
// knows gravity: level the frame by it whenever the sightings pin its tilt (Solver::Level), and keep the level in
// stations.json. A measured frame is the user's own and stays as it is. A wrong alignment (two stations swapped) can
// be levelled into fitting three stations too, so only a settled one whose worn or held devices sit at the head counts.
void Sync::LevelStep(const std::map<std::string, std::pair<V3, M3>> &raw, double now) {
  if (!sfile_.autogen || frame_.ok() != 1 || !solver_.has_x() || now - lock_time_ < kLevelSettle) return;
  if (!(BodyDist(solver_.x()) < kLevelBody)) return;
  Solver::LevelR lv;
  if (!solver_.Level(sfile_.anchor_p, lv)) return;
  double step = RotDeg(lv.tilt);
  M3 L = lv.tilt * sfile_.level;
  if (step < 0.05 || RotDeg(L) > kMaxRot) return;  // within the noise; or past what a level can be
  sfile_.level = L;
  sfile_.SaveAuto(Join(cfg_.dir, "stations.json"));
  frame_.Update(raw);
  for (auto &kv : raw) {
    V3 p; M3 R;
    frame_.Apply(kv.second.first, kv.second.second, p, R);
    solver_.SetStation(kv.first, p, R);
  }
  solver_.Relevel(lv.x);
  log_(Fmt("lighthouse frame levelled with the headset: tilted %.2f deg (%.2f deg in all), %d base stations within "
           "%.2f deg", step, RotDeg(L), lv.stations, lv.med));
}

// An automatic frame takes its stations' poses in it once three show up (Frame::Layout follows them from then on) and
// adds new ones as they come. A known station's place is never retaken from a fit of the others: SteamVR re-solves
// them tens of cm off at times, and that bent the layout for good. The frame starts over from SteamVR's (its level
// found again) when for kLayoutMoved s the fit leaves a station out, no fit gets within kMoved, or the fit tilts the
// frame over kMaxRot off level.
void Sync::LayoutStep(const std::map<std::string, std::pair<V3, M3>> &raw, double now) {
  if (!sfile_.autogen || frame_.ok() != 1) return;
  bool was = frame_.layout();
  if (!was && raw.size() < 3) return;
  if (!was || frame_.miss() <= kMoved) bad_since_ = -1;
  else if (bad_since_ < 0) bad_since_ = now;
  double tilt = TiltDeg(frame_.Anchored());
  if (!was || tilt <= kMaxRot) tilt_since_ = -1;
  else if (tilt_since_ < 0) tilt_since_ = now;
  std::string left = frame_.left();
  if (left != left_seen_) { left_seen_ = left; left_since_ = now; }
  std::string why;
  if (bad_since_ >= 0 && now - bad_since_ >= kLayoutMoved)
    why = Fmt("base stations moved (%.0f cm off the fit)", frame_.miss() * 100);
  else if (was && !left.empty() && now - left_since_ >= kLayoutMoved)
    why = "base station " + left + " moved";
  else if (tilt_since_ >= 0 && now - tilt_since_ >= kLayoutMoved)
    why = Fmt("the base stations' fit tilts it %.1f deg off SteamVR's level", tilt);
  if (!why.empty() && raw.size() >= 3) {
    sfile_.ref.clear();
    sfile_.ref.insert(raw.begin(), raw.end());
    sfile_.level = M3();  // SteamVR's frame as it is now: its level is found again
    sfile_.SaveAuto(Join(cfg_.dir, "stations.json"));
    bad_since_ = tilt_since_ = -1;
    left_since_ = now;
    log_("lighthouse frame: " + why + ": reference frame reset");
    frame_.Update(raw);
    return;
  }
  int n = 0;
  for (auto &kv : raw) {
    if (sfile_.ref.count(kv.first)) continue;
    V3 p; M3 R;
    frame_.Unlevelled(kv.second.first, kv.second.second, p, R);
    sfile_.ref[kv.first] = {p, R};
    n++;
    if (was) log_(Fmt("lighthouse frame: base station %s added to the fit", kv.first.c_str()));
  }
  if (!n) return;
  if (!was) log_(Fmt("lighthouse frame: follows the fit of %d base stations from now on", n));
  sfile_.SaveAuto(Join(cfg_.dir, "stations.json"));
  frame_.Update(raw);
}

// Two stations: the level the cameras see across them, against gravity's (gravity.h) about the same axis. Logged every
// few minutes (and recorded), so the two can be compared; nothing is applied from it.
// an automatic frame with the anchor and one other station in it (three or more: the layout fit places them)
bool Sync::TwoStations(const std::map<std::string, std::pair<V3, M3>> &raw, std::string &other) const {
  if (!sfile_.autogen || frame_.layout() || frame_.ok() != 1 || raw.size() != 2 || !raw.count(sfile_.anchor))
    return false;
  for (auto &kv : raw)
    if (kv.first != sfile_.anchor) other = kv.first;
  return true;
}

// The other station of an automatic two-station frame where the cameras see it from the anchor (Tick hands it to the
// solver): this session's measurements (of windows that overlap, the better), with what earlier sessions kept unless
// three in a row say the station moved
// The layout matcher: candidates from the 3D cell layer vs the predicted slam-space station
// positions x^-1(S_k). Enumerates assignments (n <= 4), gates each candidate pair at a FIXED
// 0.5 m (floors, not covariance -- the forged +-1.4 cm event), scores by total residual, and
// accepts only an unambiguous winner (single survivor, or margin >= 0.25 m over the runner-up).
Sync::MatchReport Sync::MatchLayout(const V3 *P, int n) {
  MatchReport mr;
  if (n < 2 || n > 4 || (int)cands_.size() < n) return mr;
  const double kTol = 0.5, kMargin = 0.25;
  std::vector<int> perm(n);
  for (int i = 0; i < n; i++) perm[i] = i;
  int best[4] = {-1, -1, -1, -1};
  double bs = 1e9, sb = 1e9;
  do {
    double s = 0;
    bool ok = true;
    for (int k = 0; k < n && ok; k++) {
      double e = norm(cands_[perm[k]].q - P[k]);
      if (e > kTol) { ok = false; break; }
      // ... and it must not land on a candidate the short frames say radiates continuously: a
      // predicted station on a lamp is not a correspondence, it is a coincidence. This is the guard
      // that catches a wrong acquisition quickly, so it is the one place the reject state is hard.
      if (cands_[perm[k]].kind == QK_REJECT) { ok = false; break; }
      s += e;
    }
    if (!ok) continue;
    mr.npass++;
    if (s < bs) {
      sb = bs; bs = s;
      for (int k = 0; k < n; k++) best[k] = perm[k];
    } else if (s < sb) {
      sb = s;
    }
  } while (std::next_permutation(perm.begin(), perm.end()));
  mr.best = mr.npass ? bs : -1;
  mr.second = mr.npass > 1 ? sb : 0;
  mr.a0 = best[0];
  mr.a1 = best[1];
  mr.ok = mr.npass >= 1 && (mr.npass == 1 || sb - bs >= kMargin);
  match_ok_ = mr.ok;
  return mr;
}

void Sync::PlaceStep(const std::map<std::string, std::pair<V3, M3>> &raw, double now) {
  std::string other;
  if (!TwoStations(raw, other)) return;
  // Placement is x and SteamVR's own layout, period: the camera-measured baseline override is
  // RETIRED. It had flip-flopped 4.41 -> 4.85 m (a forged +-1.2 cm) while SteamVR held 4.92 m:
  // the baseline is SteamVR's to state, the cameras only ever produce x. What remains here is the
  // layout matcher's verdict, logged for the evidence watch. Ambiguous or nothing passing leaves
  // x untouched -- there is no second placement path, and no code models a station moving.
  std::vector<std::string> keys;
  std::vector<V3> S, Z, P, Zq;
  solver_.Stations(keys, S, Z);
  Solver::Predict(solver_.x(), S, Z, P, Zq);
  MatchReport mr = MatchLayout(P.data(), (int)P.size());
  std::string sig = Fmt("%d|%d|%.2f", (int)mr.ok, mr.npass, mr.best);
  if (sig != match_sig_ || now > match_next_) {
    match_sig_ = sig;
    match_next_ = now + 60;
    if (mr.npass)
      log_(Fmt("layout match: %d assignment(s) pass, best %.2f m, second %.2f m, margin %.2f m -- %s",
               mr.npass, mr.best, mr.second, mr.margin, mr.ok ? "accepted" : "AMBIGUOUS"));
    else
      log_("layout match: no candidate within 0.5 m of the predicted stations -- holding the current x");
  }
}

void Sync::LevelCheck(double now) {
  std::vector<std::string> keys;
  std::vector<V3> S, Z;
  solver_.Stations(keys, S, Z);
  if (S.size() != 2) return;  // with three, LevelStep levels the frame by the cameras themselves
  if (!match_ok_) return;     // levelling waits for an unambiguous layout match: a decoy cluster must not tilt the frame
  // a measured station stays where the cameras put it, gravity's levelling or not: they'd only check themselves
  if (solver_.measured(keys[0]) || solver_.measured(keys[1])) return;
  V3 b = S[1] - S[0];
  double h = std::hypot(b.x, b.z);
  if (h < 1) return;
  V3 axis{-b.z / h, 0, b.x / h}, pivot = (S[0] + S[1]) * 0.5;
  Solver::AxisR a;
  if (!solver_.LevelAxis(pivot, axis, a)) return;
  // gravity's tilt in the frame now, about the same axis: the cameras' view of the frame before it is the two together
  const M3 &G = frame_.gravity();
  V3 rg{G.m[2][1] - G.m[1][2], G.m[0][2] - G.m[2][0], G.m[1][0] - G.m[0][1]};
  double g = dot(rg, axis) * 0.5, cam = a.theta + g;
  Rec(now, "I level %.5f %.5f %.5f %.4f", cam * kDeg, a.sd * kDeg, g * kDeg, a.med);
  if (now - check_said_ < 180) return;
  check_said_ = now;
  log_(Fmt("level check across %s and %s: the cameras see the frame %+.2f deg off level (within %.2f) before "
           "gravity's levelling, which turns it %+.2f deg about the same axis%s", keys[0].c_str(), keys[1].c_str(),
           cam * kDeg, a.sd * kDeg, g * kDeg, std::fabs(a.theta) > 3 * a.sd ? "" : ": they agree"));
}

// two stations have 10 sightings that fit
static bool Locked(bool has_x, const StepStat &st) {
  if (!has_x || !st.has_per) return false;
  std::vector<int> c;
  for (auto &kv : st.per) c.push_back(kv.second);
  std::sort(c.begin(), c.end());
  return c.size() >= 2 && c[c.size() - 2] >= 10;
}

Transform Sync::Tick(double now) {
  ApplyHeadset();
  if (optics_.exact_time() && !solver_.has_x() && now > acq_next_) {
    acq_next_ = now + 15;   // the enumerator gets its own fair window, not just sweep ticks that need a pose
    AcquireFromCandidates(now);
  }
  solver_.SetManualAcq(optics_.exact_time());   // lighthouse axiom: locks are only ever earned through the candidate pipeline; there is no fallback acquirer
  if (now - last_step_ >= 1.0) {
    last_step_ = now;
    int b = poses_.brk_id();
    if (b != brk_seen_) {
      brk_seen_ = b;
      solver_.PoseBreak(poses_.brk_t());
      log_("HMD pose break (" + poses_.brk_what() + "): older sightings dropped, acquiring");
    }
    std::map<std::string, std::pair<V3, M3>> raw;
    std::map<std::string, int> chans;
    {
      std::lock_guard<std::mutex> g(raw_m_);
      raw = raw_;
      chans = channels_;
    }
    for (auto &kv : chans) solver_.SetChannel(kv.first, kv.second);
    if (!sfile_.loaded && !raw.empty()) {  // first run: pin the reference frame to the station at SteamVR's origin
      auto best = raw.begin();
      for (auto it = raw.begin(); it != raw.end(); ++it)
        if (norm(it->second.first) < norm(best->second.first)) best = it;
      sfile_.loaded = sfile_.autogen = true;
      sfile_.anchor = best->first;
      sfile_.anchor_p = best->second.first;
      sfile_.anchor_R = best->second.second;
      sfile_.SaveAuto(Join(cfg_.dir, "stations.json"));
      log_("reference frame set: anchor " + best->first);
    } else if (sfile_.autogen && !frame_.layout() && raw.count(sfile_.anchor)) {  // follows a moved anchor
      auto &a = raw[sfile_.anchor];
      if (RotDeg(sfile_.anchor_R * T(a.second)) > kMaxRot || norm(a.first - sfile_.anchor_p) > kMoved) {
        sfile_.anchor_p = a.first;
        sfile_.anchor_R = a.second;
        sfile_.level = M3();  // SteamVR's frame as it is now: its level is found again
        sfile_.SaveAuto(Join(cfg_.dir, "stations.json"));
        log_("anchor " + sfile_.anchor + " moved: reference frame reset");
      }
    }
    frame_.Update(raw);
    LayoutStep(raw, now);
    if (sfile_.autogen) {  // two stations: the other where the cameras measured it, in SteamVR's direction from the anchor
      std::string other;
      bool two = TwoStations(raw, other);
      for (auto &kv : sfile_.place) {
        const StationsFile::Place &pl = kv.second;
        auto it = raw.find(kv.first);
        if (!two || kv.first != other || it == raw.end() || pl.sd_d > kPlaceSdD || pl.sd_h > kPlaceSdH) {
          solver_.DropFix(kv.first);
          continue;
        }
        const auto &an = raw.at(sfile_.anchor);
        V3 pa, pb;
        M3 R;
        frame_.Apply(an.first, an.second, pa, R);
        frame_.Apply(it->second.first, it->second.second, pb, R);
        V3 v = pb - pa;
        double h = std::hypot(v.x, v.z);
        if (h < 1) solver_.DropFix(kv.first);
        else solver_.SetFix(kv.first, pa + V3{v.x / h * pl.d, pl.h, v.z / h * pl.d});
      }
    }
    std::map<std::string, V3> cs;
    for (auto &kv : raw) {
      V3 p; M3 R;
      frame_.Apply(kv.second.first, kv.second.second, p, R);
      solver_.SetStation(kv.first, p, R);
      if (now - last_s_ >= 10) {
        frame_.Unlevelled(kv.second.first, kv.second.second, p, R);
        Quat q = ToQuat(R);
        auto c = chans.find(kv.first);
        Rec(now, "S %s %.5f %.5f %.5f %.5f %.5f %.5f %.5f %d", kv.first.c_str(), p.x, p.y, p.z, q.w, q.x, q.y, q.z,
            c == chans.end() ? 0 : c->second);
      }
    }
    if (now - last_s_ >= 10 && !raw.empty()) last_s_ = now;
    if (recording_) {  // worn and held lighthouse devices, in the reference frame like the stations (before the level)
      std::deque<BodyS> nb;
      {
        std::lock_guard<std::mutex> g(body_m_);
        nb.swap(body_new_);
      }
      for (const BodyS &e : nb) {
        V3 p;
        M3 R;
        frame_.Unlevelled(e.p, M3(), p, R);
        Rec(e.t, "B %d %.4f %.4f %.4f", e.dev, p.x, p.y, p.z);
      }
    }
    {
      std::vector<std::string> keys;
      std::vector<V3> S, Z;
      solver_.Stations(keys, S, Z);
      for (size_t i = 0; i < keys.size(); i++)
        if (raw.count(keys[i])) cs[keys[i]] = S[i];
    }
    if (!seeded_ && cs.size() >= 2) {
      seeded_ = true;
      X4 x;
      double miss;
      if (Seed(cs, x, miss) && miss < 0.05) {
        solver_.Reset(x);
        log_(Fmt("seeded from saved stations: yaw %+.3f t [%.4f, %.4f, %.4f] miss %.1f cm", x[0] * kDeg, x[1], x[2], x[3], miss * 100));
      } else {
        log_("saved stations don't fit: acquiring");
      }
    }
    StepStat st = solver_.Step(now);
    if (solver_.resets() != resets_seen_) { resets_seen_ = solver_.resets(); lock_time_ = now; }
    {
      std::lock_guard<std::mutex> g(st_m_);
      last_st_ = st;
      has_last_st_ = true;
    }
    LogStages(now);
    // not acquired: what became of the cameras' bright spots, after a minute and then every 5
    if (Locked(solver_.has_x(), st)) sum_t_ = -1;
    else {
      Spots sp = spots();
      if (sum_t_ < 0) { sum_t_ = now; sum_sp_ = sp; sum_n_ = 0; }
      else if (now - sum_t_ >= (sum_n_ ? 300 : 60)) {
        if (sp.frames > sum_sp_.frames) log_(Fmt("acquiring, last %.0f s: ", now - sum_t_) + Describe(sum_sp_, sp));
        sum_t_ = now; sum_sp_ = sp; sum_n_++;
      }
    }
    X4 x = solver_.x();
    if (recording_) {
      std::string per;
      for (auto &kv : st.per) per += Fmt("%s'%s': %d", per.empty() ? "{" : ", ", kv.first.c_str(), kv.second);
      per = st.has_per ? per + "}" : "None";
      Rec(now, "A %s | %s | %d %s %s%s%s", have_applied_ ? Fmt("%+.4f %.4f %.4f %.4f", applied_[0] * kDeg, applied_[1], applied_[2], applied_[3]).c_str() : "-",
          solver_.has_x() ? Fmt("%+.4f %.4f %.4f %.4f", x[0] * kDeg, x[1], x[2], x[3]).c_str() : "-", st.n,
          st.has_med ? Fmt("%.4f", st.med).c_str() : "None", per.c_str(), st.cond ? " cond" : "", st.jump ? " JUMP" : "");
    }
    if (Locked(solver_.has_x(), st) && st.cond && st.has_med && st.med < 0.3 && now - last_save_ > 10) {
      last_save_ = now;
      SaveState(x, cs);
    }
    // timing: every 15 s on a well-pinned alignment (under 0.4 deg: a wrong timing spreads the sightings itself, and
    // a Frame's at -6.2 ms sat at 0.33, where 0.3 kept it from ever being corrected)
    if (cfg_.learn_timing && solver_.has_x() && st.cond && st.has_med && st.med < 0.4 && now - last_timing_ > 15) {
      last_timing_ = now;
      std::vector<std::string> keys;
      std::vector<V3> S, Z;
      solver_.Stations(keys, S, Z);
      double best;
      int n;
      Timing::Split sp;
      bool ok;
      double cur = expo_;
      {
        std::lock_guard<std::mutex> g(net_);
        ok = timing_.Estimate(poses_, x, S, Z, cur, !timing_learned_, best, n, sp);
      }
      if (!ok && std::isfinite(sp.left) && now - split_said_ >= 300) {
        split_said_ = now;
        log_(Fmt("timing: turning left fits %.1f ms, turning right %.1f ms%s: too unsure to take, it stays %.1f ms",
                 sp.left * 1000, sp.right * 1000,
                 std::isfinite(sp.ci) ? Fmt(", their mean within %.1f ms", sp.ci * 1000).c_str() : "", cur * 1000));
      }
      if (ok) {
        // the step trusts the estimate by how much turning it saw, both ways
        double e = timing_learned_ ? cur + std::min(1.0, 2.0 * std::min(sp.nl, sp.nr) / kTimingFullN) * (best - cur)
                                   : best;
        e = std::max(-0.05, std::min(0.12, e));
        expo_ = e;
        if (std::fabs(e - cur) > 1e-4) {
          std::lock_guard<std::mutex> g(net_);
          Retime(cur, e);
        }
        bool first = !timing_learned_;
        timing_learned_ = true;
        std::string sys;
        {
          std::lock_guard<std::mutex> g(st_m_);
          sys = streamer_;
        }
        if (!sys.empty()) {
          std::lock_guard<std::mutex> g(st_m_);
          timings_[sys] = e;
        }
        if (first || std::fabs(best - cur) > 0.002)
          log_(Fmt("timing%s%s: %.1f ms (best fit %.1f ms within %.1f, on %d fast-head sightings: %.1f turning left, "
                   "%.1f right)", sys.empty() ? "" : " for ", sys.c_str(), e * 1000, best * 1000, sp.ci * 1000, n,
                   sp.left * 1000, sp.right * 1000));
        Rec(now, "I timing %.5f best %.5f ci %.5f n %d left %.5f %d right %.5f %d", e, best, sp.ci, n, sp.left, sp.nl,
            sp.right, sp.nr);
      }
    }
    if (now - last_level_ >= 10) {
      last_level_ = now;
      LevelStep(raw, now);
    }
    if (now - last_place_ >= kPlaceEvery && Locked(solver_.has_x(), st) && st.cond) {
      last_place_ = now;
      PlaceStep(raw, now);
    }
    if (now - last_check_ >= 60 && Locked(solver_.has_x(), st) && st.cond) {
      last_check_ = now;
      LevelCheck(now);
    }
  }
  // slew toward the solution, faster while the head turns or walks (kSlew*), jump if far (acquisition)
  Transform out;
  if (!solver_.has_x()) return out;
  X4 x = solver_.x();
  if (!have_applied_) { applied_ = x; have_applied_ = true; }
  else if (!paused_) {
    V3 h{0, 1.5, 0};
    double th, w = 0, v = 0;
    if (poses_.Latest(h, &th)) {
      M3 R0;
      V3 p0;
      if (!poses_.Speed(th, w)) w = 0;
      if (poses_.At(th - 0.1, R0, p0)) v = norm(h - p0) / 0.1;
    }
    X4 dx = {Wrap(x[0] - applied_[0]), x[1] - applied_[1], x[2] - applied_[2], x[3] - applied_[3]};
    V3 at{applied_[1], applied_[2], applied_[3]};
    double d = norm(Ry(dx[0]) * (h - at) + at + V3{dx[1], dx[2], dx[3]} - h);
    if (d > kJumpApply) applied_ = x;
    else if (d > 1e-6) {
      double s = std::min(1.0, (kSlewStill + kSlewTurn * w + kSlewWalk * v) * 0.05 / d);
      for (int i = 0; i < 4; i++) applied_[i] += dx[i] * s;
    }
  }
  M3 C; V3 tc;
  frame_.Rotation(C, tc);
  M3 R = Ry(applied_[0]);
  out.active = true;
  out.q = ToQuat(R * C);
  out.t = R * tc + V3{applied_[1], applied_[2], applied_[3]};
  return out;
}

Sync::Status Sync::GetStatus(double now) {
  Status s;
  Spots sp = spots();
  if (now - rate_t_ >= 2.0) {
    long f = nframes_;
    if (rate_t_ > 0) {
      double dt = now - rate_t_;
      cam_fps_ = (f - rate_f_) / dt;
      spot_rate_ = (sp.seen - rate_sp_.seen) / dt;
      sight_rate_ = (sp.used - rate_sp_.used) / dt;
    }
    rate_t_ = now; rate_f_ = f; rate_sp_ = sp;
  }
  s.cam_fps = cam_fps_;
  s.spot_rate = spot_rate_;
  s.sight_rate = sight_rate_;
  double th;
  V3 h;
  s.head_still = poses_.Latest(h, &th) && poses_.Still(th);
  // still for kFrozen s while the lamps and windows the cameras see keep moving: SteamVR's pose isn't the head's
  if (!s.head_still || !(cam_fps_ > 5)) still_since_ = -1;
  else if (still_since_ < 0) still_since_ = now;
  long moved, stayed;
  {
    std::lock_guard<std::mutex> g(net_);
    moved = img_moved_;
    stayed = img_still_;
  }
  bool frozen = still_since_ >= 0 && now - still_since_ >= kFrozen && moved >= 20 && moved > 3 * stayed;
  if (frozen) {
    s.head_still = 2;
    if (!frozen_said_)
      log_(Fmt("SteamVR's headset pose hasn't moved for %.0f s while the cameras see the room move (%ld lamp or window "
               "spots moved, %ld stayed): SteamVR isn't getting the head's motion. Is its view in the headset?",
               now - still_since_, moved, stayed));
  }
  frozen_said_ = frozen || (frozen_said_ && still_since_ >= 0);
  s.has_x = solver_.has_x();
  s.x = solver_.x();
  s.expo = expo_;
  s.timing_learned = timing_learned_;
  {
    std::lock_guard<std::mutex> g(net_);
    s.rtt = clock_.rtt();
  }
  StepStat st;
  {
    std::lock_guard<std::mutex> g(st_m_);
    st = last_st_;
  }
  s.cond = st.cond;
  s.n = st.n;
  s.med = st.has_med ? st.med : -1;
  s.locked = Locked(s.has_x, st);
  s.locked_for = lock_time_ >= 0 ? now - lock_time_ : -1;
  if (s.has_x && have_applied_) {
    V3 h{0, 1.5, 0};
    poses_.Latest(h);
    X4 dx = {Wrap(s.x[0] - applied_[0]), s.x[1] - applied_[1], s.x[2] - applied_[2], s.x[3] - applied_[3]};
    V3 at{applied_[1], applied_[2], applied_[3]};
    s.lag_cm = norm(Ry(dx[0]) * (h - at) + at + V3{dx[1], dx[2], dx[3]} - h) * 100;
  }
  std::vector<std::string> keys;
  std::vector<V3> S, Z;
  solver_.Stations(keys, S, Z);
  s.nstations = (int)keys.size();
  V3 head;
  bool hh = poses_.Latest(head);
  for (size_t i = 0; i < keys.size(); i++) {
    Status::St e;
    e.serial = keys[i];
    e.anchor = sfile_.loaded && keys[i] == sfile_.anchor;
    e.measured = solver_.measured(keys[i]);
    auto it = st.per.find(keys[i]);
    e.support = it == st.per.end() ? 0 : it->second;
    auto lt = st.last.find(keys[i]);
    e.last_seen = lt == st.last.end() ? -1 : std::max(0.0, now - lt->second);
    if (s.has_x && hh) e.dist = norm(Ry(s.x[0]) * S[i] + V3{s.x[1], s.x[2], s.x[3]} - head);
    s.st.push_back(e);
  }
  return s;
}
