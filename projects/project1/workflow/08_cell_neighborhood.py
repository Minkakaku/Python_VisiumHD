# 按样本计算局部邻域组成，再识别 CN 并保存标签。

from pathlib import Path

from visium_hd_workflow.project import load_project
from visium_hd_workflow.steps import step_cell_neighborhood

# 由脚本所在位置定位项目目录，避免依赖启动命令时的工作目录。
step_cell_neighborhood(load_project(Path(__file__).resolve().parents[1]))
