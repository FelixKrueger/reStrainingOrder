#!/usr/bin/env bash
# Build the release binary and expose it under the three classic names in
# rust/target/impl/, which is the directory test/run_tests.py --impl takes.
set -euo pipefail
cd "$(dirname "$0")"
cargo build --release --locked
mkdir -p target/impl
for name in reStraining reStrainingOrder reStrainingReport; do
    ln -sf ../release/restrainingorder "target/impl/$name"
done
echo "target/impl/ ready: reStraining, reStrainingOrder, reStrainingReport"
