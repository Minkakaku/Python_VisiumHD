#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$project_dir/logs"
for step in 00_setup 01_load 02_qc; do
  python "$project_dir/workflow/$step.py" 2>&1 | tee "$project_dir/logs/$step.log"
done
printf '%s\n' '请检查 QC 阈值，再单独执行聚类和 marker；人工注释后再运行下游分析。'
