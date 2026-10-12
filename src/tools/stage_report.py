#!/usr/bin/env python3
# Stage report: how well a recorded session fares at each stage of the pipeline.
#
# Runs the recording through the current driver code offline (qlhs_replay) and through the
# independent ground-truth tool (qlhs_presence), then lays the two side by side:
#
#   1 collection  what each camera and class delivered              (raw, nothing gated)
#   2 filtering   what each gate took, and the yield to a ray       (dots -> rays)
#   3 places      cells -> 3D fits -> candidates                    (directions -> 3D points)
#   4 solve       did it acquire, onto what, on how much support    (against the ground truth)
#
# usage: stage_report.py <recording> [--calib f] [--replay f] [--presence f]
#                        [--frame-dir fresh] [--top N]
import argparse, os, re, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
DEF_REPLAY = "/tmp/qlhs_replay"
DEF_PRESENCE = "/tmp/qlhs_presence"


def run(cmd, **kw):
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=1800, **kw)
    return p.stdout + p.stderr


def parse_replay(text):
    d = {"collect": [], "stages": [], "acquired": [], "acq": [], "naming": [], "sgates": 0,
         "layout": [], "merges": 0, "presence": [], "final": ""}
    for ln in text.splitlines():
        s = ln.split("s ", 1)[-1] if re.match(r"^\s*\d+s ", ln) else ln
        s = s.strip()
        if s.startswith("collect"): d["collect"].append(s)
        elif s.startswith("stages:"): d["stages"].append(s)
        elif s.startswith("acquired"): d["acquired"].append(s)
        elif s.startswith("acquisition:"): d["acq"].append(s)
        elif s.startswith("candidate naming"): d["naming"].append(s)
        elif s.startswith("S-gate"): d["sgates"] += 1
        elif s.startswith("layout match"): d["layout"].append(s)
        elif s.startswith("merge:"): d["merges"] += 1
        elif s.startswith("cameras:"): d["presence"].append(s)
        elif s.startswith("replayed"): d["final"] = s
    return d


def sum_collect(lines):
    """collect lines are per-interval; sum them into per camera/class totals."""
    tot = {}
    for ln in lines:
        for m in re.finditer(r"c(\d)([LS]) (\d+)f (\d+)d (\d+)b pk(\d+)", ln):
            k = (int(m.group(1)), m.group(2))
            a = tot.setdefault(k, [0, 0, 0, 0])
            a[0] += int(m.group(3)); a[1] += int(m.group(4)); a[2] += int(m.group(5))
            a[3] = max(a[3], int(m.group(6)))
    return tot


def sum_stages(lines):
    """The funnel counters, summed over intervals. Parsed generically ("name 1234") so the
    labels can move in the driver's line without the report silently mis-labelling a stage."""
    tot = {}
    for ln in lines:
        m = re.search(r"(\d+) long frames -> (\d+) dots -> (\d+) rays", ln)
        if m:
            for k, v in zip(["frames", "dots", "rays"], m.groups()):
                tot[k] = tot.get(k, 0) + int(v)
        seg = ln.split("dropped:", 1)[-1].split("|", 1)[0]
        for k, v in re.findall(r"(still|fast|big|nopose|dim|S-gated|fuzzy|unmatched) (\d+)", seg):
            k = {"S-gated": "sgated"}.get(k, k)
            tot[k] = tot.get(k, 0) + int(v)
        m = re.search(r"\+(\d+) short-frame evidence dots", ln)
        if m:
            tot["short"] = tot.get("short", 0) + int(m.group(1))
    return tot


def last_cells(lines):
    out = []
    for ln in lines:
        m = re.search(r"cells (\d+) \(3D (\d+), fuzzy (\d+)\) -> candidates (\d+)", ln)
        if m: out = [int(x) for x in m.groups()]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("recording")
    ap.add_argument("--calib", default="/tmp/calib.json")
    ap.add_argument("--replay", default=DEF_REPLAY)
    ap.add_argument("--presence", default=DEF_PRESENCE)
    ap.add_argument("--frame-dir", default="")
    ap.add_argument("--top", type=int, default=6)
    ap.add_argument("--x", default="", help='"yaw tx ty tz" of a KNOWN-GOOD fit, for the ground truth')
    a = ap.parse_args()

    # the replay --dir must be fresh: state.json in it seeds the solver, and a reused dir replays a
    # different session than the recording describes.
    dr = a.frame_dir or tempfile.mkdtemp(prefix="qlhs_stage_")
    if a.frame_dir:
        shutil.rmtree(dr, ignore_errors=True)
        os.makedirs(dr, exist_ok=True)
    rays = os.path.join(dr, "rays.txt")
    print("== replaying", os.path.basename(a.recording), "->", dr)
    txt = run([a.replay, a.recording, "--calib", a.calib, "--dir", dr, "--rays", rays])
    d = parse_replay(txt)

    print("\n" + "=" * 100)
    print("1  COLLECTION  (what the cameras delivered, nothing gated)")
    print("=" * 100)
    c = sum_collect(d["collect"])
    if not c:
        print("  no collect lines -- is this the instrumented build?")
    tot_blobs = 0
    for (cam, cls) in sorted(c):
        fr, dt, bl, pk = c[(cam, cls)]
        tot_blobs += bl
        print("  cam %d %s: %6d frames, %6d carried a blob (%3.0f%%), %7d blobs, peak %d, %.2f blobs/frame"
              % (cam, cls, fr, dt, 100.0 * dt / max(1, fr), bl, pk, bl / max(1, fr)))

    print("\n" + "=" * 100)
    print("2  FILTERING  (every gate a dot meets, in order)")
    print("=" * 100)
    s = sum_stages(d["stages"])
    if s["dots"]:
        order = [("rays", "kept, became a ray"), ("big", "lamps or windows (too many pixels)"),
                 ("dim", "dim (room too light to saturate)"), ("fast", "head turning too fast"),
                 ("nopose", "no headset pose"), ("still", "headset still (lamps move too)"),
                 ("sgated", "S-gate: steady sources (parked, on purpose)"),
                 ("fuzzy", "failed the 3D point test"), ("unmatched", "layout matcher did not accept")]
        for k, what in order:
            print("  %8d  %5.1f%% of the long dots   %s" % (s[k], 100.0 * s[k] / max(1, s["dots"]), what))
        print("  --- %d long dots from %d frames; %.0f%% arrive as a ray"
              % (s["dots"], s["frames"], 100.0 * s["rays"] / max(1, s["dots"])))
        print("  --- %d short-frame dots, kept as evidence only (never rays by design)" % s["short"])
    else:
        print("  no stages lines")

    print("\n" + "=" * 100)
    print("3  PLACES  (cells -> 3D fits -> candidates)")
    print("=" * 100)
    lc = last_cells(d["stages"])
    if lc:
        cells, tri, fuzz, cand = lc
        print("  cells %d, ever got a 3D fit %d (%.1f%%), of those fuzzy %d -> candidates %d"
              % (cells, tri, 100.0 * tri / max(1, cells), fuzz, cand))
        print("  merges logged: %d" % d["merges"])
    else:
        print("  no cell/candidate lines")

    print("\n" + "=" * 100)
    print("4  SOLVE")
    print("=" * 100)
    for t in d["acq"][-3:]:
        print("  " + t)
    if d["acquired"]:
        print("  ACQUIRED: " + d["acquired"][-1])
    else:
        print("  never acquired in this session")
    for t in d["layout"][-2:]:
        print("  " + t)
    if d["sgates"]:
        print("  S-gate verdict lines: %d (parked: computed, not applied)" % d["sgates"])

    # the ground truth, from the accepted x if there is one
    x = None
    from_run = False
    for ln in reversed(open(rays).read().splitlines()):
        if ln.startswith("X "):
            p = ln.split()
            x = (float(p[1]), float(p[2]), float(p[3]), float(p[4]))
            from_run = True
            break
    if a.x:
        x = tuple(float(v) for v in a.x.split())
        from_run = False
    # the station places the ground truth needs come from x. If it came from THIS run and the layout
    # matcher then said the fit is not on the stations, that x is wrong and so is the ground truth.
    bad = any("not on the stations" in t or "no candidate within" in t for t in d["layout"])
    if x and os.path.exists(a.presence):
        print("\n  x used for the station places: yaw %+.2f t [%.3f %.3f %.3f]%s"
              % (x[0], x[1], x[2], x[3], "  (from this run)" if from_run else "  (supplied)"))
        if from_run and bad:
            print("  !! this run's own fit was rejected by the layout matcher, so the places below are"
                  " downstream of a wrong x: pass --x with a known-good fit")
        print("\n" + "=" * 100)
        print("GROUND TRUTH  (qlhs_presence: the raw pixels, independent of the driver)")
        print("=" * 100)
        pt = run([a.presence, a.recording, "--calib", a.calib, "--yaw", "%f" % x[0], "--t", "%f" % x[1], "%f" % x[2],
                  "%f" % x[3]])
        for ln in pt.splitlines():
            if re.search(r"LHB|pose offset|in cell", ln):
                print("  " + ln.strip())
    print("\n" + d["final"])
    if not a.frame_dir:
        shutil.rmtree(dr, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
