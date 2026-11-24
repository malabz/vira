#!/usr/bin/env bash
set -euo pipefail

make release
mkdir -p "${PREFIX}/bin"
install -m 0755 vira "${PREFIX}/bin/vira"
