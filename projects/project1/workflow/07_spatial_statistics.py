# 读取注释对象，按样本执行空间统计并输出图表。

from pathlib import Path

from visium_hd_workflow.project import load_project
from visium_hd_workflow.steps import step_spatial_statistics

# 由脚本所在位置定位项目目录，避免依赖启动命令时的工作目录。
step_spatial_statistics(load_project(Path(__file__).resolve().parents[1]))
