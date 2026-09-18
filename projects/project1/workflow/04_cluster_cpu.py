# 在 CPU 上执行 Visium HD 归一化、降维与聚类。

from pathlib import Path

from visium_hd_workflow.project import load_project
from visium_hd_workflow.steps import step_cluster_cpu

# 由脚本所在位置定位项目目录，避免依赖启动命令时的工作目录。
step_cluster_cpu(load_project(Path(__file__).resolve().parents[1]))
