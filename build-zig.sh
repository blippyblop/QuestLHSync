#!/bin/bash
# Cross-builds driver/questlhsync/bin/win64/driver_questlhsync.dll without MSVC,
# using the zig toolchain (tools/zig) as x86_64-windows-gnu. Same sources, defines
# and link libs as build.bat's driver step; different toolchain (clang+mingw-w64
# instead of MSVC), so the bytes/size differ from a build.bat build. build.bat on
# Windows remains the canonical build; use this when only this container is at hand.
# Freshness check afterwards: size + string presence (see AGENTS notes / release.py).
set -e
cd "$(dirname "$0")"
ZIG=${ZIG:-tools/zig/zig}
export PATH="$(dirname "$ZIG"):$PATH"
OUT=driver/questlhsync/bin/win64
mkdir -p build/driver "$OUT"
# MinHook (plain C)
for f in buffer hook trampoline; do
  zig cc -c -O2 -DNOMINMAX -DWIN32_LEAN_AND_MEAN -D_CRT_SECURE_NO_WARNINGS \
    -Ithird_party/minhook/include -target x86_64-windows-gnu \
    third_party/minhook/src/$f.c -o build/driver/$f.o
done
zig cc -c -O2 -DNOMINMAX -DWIN32_LEAN_AND_MEAN -Ithird_party/minhook/include \
  -target x86_64-windows-gnu third_party/minhook/src/hde/hde64.c -o build/driver/hde64.o
# Driver (C++)
for f in driver_main sync net gravity relations; do
  zig c++ -c -O2 -fexceptions -std=c++17 -DNOMINMAX -DWIN32_LEAN_AND_MEAN \
    -D_CRT_SECURE_NO_WARNINGS -Ithird_party/openvr/headers -Ithird_party/minhook/include \
    -target x86_64-windows-gnu src/driver/$f.cpp -o build/driver/$f.o
done
# Link
rm -f "$OUT"/driver_questlhsync.dll
zig c++ -shared -target x86_64-windows-gnu -o "$OUT"/driver_questlhsync.dll \
  build/driver/*.o -lws2_32 -liphlpapi -lshell32 -lole32 -lsetupapi -lhid
rm -f "$OUT"/*.lib "$OUT"/*.pdb
echo "built $OUT/driver_questlhsync.dll ($(stat -c%s "$OUT/driver_questlhsync.dll") bytes)"
echo "freshness: grep -ac 'presence checks, highest' '$OUT/driver_questlhsync.dll'  # expect >= 1"
