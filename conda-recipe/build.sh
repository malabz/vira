#!/usr/bin/env bash
set -euo pipefail

make clean
make -j"${CPU_COUNT:-1}" \
  CXX="${CXX}" \
  CPPFLAGS="${CPPFLAGS:-} -Iinclude" \
  CXXFLAGS="${CXXFLAGS:-} -std=c++17 -O3 -pthread -Wall -Wextra" \
  LDFLAGS="${LDFLAGS:-}"
mkdir -p "${PREFIX}/bin"
install -m 0755 vira "${PREFIX}/bin/vira"
