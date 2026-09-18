# 计算样本内最近邻距离，汇总为样本层面的结果。

from pathlib import Path

from visium_hd_workflow.project import load_project
from visium_hd_workflow.steps import step_nearest_neighbor

# 由脚本所在位置定位项目目录，避免依赖启动命令时的工作目录。
step_nearest_neighbor(load_project(Path(__file__).resolve().parents[1]))
