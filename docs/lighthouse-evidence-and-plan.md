# Lighthouse detection: field evidence & approved next steps (2026-10-11)

Status: analysis record. NOTHING here is implemented yet — waiting for a go.

## Confirmed by live runs

1. **Core alignment (x) is solid.** The worn/held tracker lands in the right place in
   sync space even in sessions where a station device lands somewhere nonsensical.
   The 4-DOF fit is not the problem.

2. **Rotor-beat separation works.** First in-the-wild evidence (18:46 candidates line):

       q[2.5 0.7 0.5] beat 9.29   q[0.0 0.0 0.0] beat 6.74 (10 arms)
       q[0.3 0.8 2.1] beat 5.25   q[1.1 -0.3 5.7] beat 3.23
       q[0.5 0.7 1.7] beat 0.00   q[3.1 0.7 4.0] beat 0.00

   The amplitude ring (peak dips at rotor-beat nulls) carries the signal; the binary
   presence ring pegs in lit rooms and sees nothing (that was commit c8b5ca8's fix).

3. **The PlaceStep camera override is the placement culprit.** Same session, its own log:

       18:05  cameras: 4.406 m baseline (within 1.4 cm)   SteamVR: 4.918 m
       18:06  cameras: 4.422 m baseline (within 1.2 cm)   SteamVR: 4.918 m
       18:46  cameras: 4.853 m baseline (within 1.2 cm)   SteamVR: 4.915 m

   A 45 cm swing with forged certainty, SteamVR steady at ~4.92 m. Per-serial
   Triangulate is polluted by the decoy family near the lamp-adjacent station, and the
   override then places the device at the corrupted baseline. This is the "one lighthouse
   in a nonsense position" failure.

4. **Acquisition is the remaining root problem.** This session acquired +147, +149,
   +161 deg — three decoy locks, each rejected by the layout matcher and killed by the
   re-acquire escalation, then re-locked elsewhere. The immune system works; the
   infection source is that initial acquisition (no x yet) eats every ray.

## Current pipeline (all shipped)

harvest (3D cell fits) -> starburst merge (0.3 m, family keys) -> S-frame verdict
(pooled per source; mains families erased) -> beat evidence (shadow-logged) ->
layout matcher (fixed 0.5 m floors, shadow for placement, gates tracking rays,
3-fail re-acquire) -> solver.

## Approved, NOT built (waiting for go)

1. **Acquisition enumeration.** Before accepting any x: enumerate candidate-pair
   assignments (merged, S-trimmed, beat-scored candidates), solve closed-form 4-DOF
   (yaw + translation Procrustes) per assignment, score by ray support, accept only an
   unambiguous winner (single survivor or clear margin). Rotor-beat score is the
   tie-breaker when a decoy pair and the real pair both fit the arrangement.
   This replaces "acquire onto whatever RANSAC finds first."

2. **Retire the PlaceStep override.** Axiom: SteamVR's station layout S_k is known
   good and is never modified, estimated, or overridden by camera data; the only
   unknown is x. Flip kMatcherGates, return PlaceStep before the Triangulate override,
   and wipe the stale `place` key from stations.json (it holds the corrupted 4.332 m
   measurement from the 18:0x decoy lock).

3. **Piece 3 lifecycle** (after 1+2): candidates move unknown -> candidate ->
   promoted -> rejected on sustained multi-signal evidence (beat, S-absence,
   saturation stability, 3D residual, persistence) with hysteresis; matcher and
   solver read promoted candidates only.

## Known knobs / facts for that build

- Merge radius 0.3 m (user's room has a station ~0.5 m from a ceiling light).
- Matcher gates: 0.5 m fixed floor per pairing, 0.25 m unambiguity margin.
- Beat scores observed: stations 3.2-9.3, flat decoys 0.00 (small sample; calibrate
  from the candidates: log lines before wiring promotion thresholds).
- ring fs assumption: 33 Hz long-frame cadence.
