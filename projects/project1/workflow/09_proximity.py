# 按配置半径和距离环带汇总样本内的邻近关系。

from pathlib import Path

from visium_hd_workflow.project import load_project
from visium_hd_workflow.steps import step_proximity

# 由脚本所在位置定位项目目录，避免依赖启动命令时的工作目录。
step_proximity(load_project(Path(__file__).resolve().parents[1]))
