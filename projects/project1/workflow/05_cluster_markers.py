# 读取已聚类对象，输出聚类展示图和 marker 分析结果。

from pathlib import Path

from visium_hd_workflow.project import load_project
from visium_hd_workflow.steps import step_cluster_plots_and_markers

# 由脚本所在位置定位项目目录，避免依赖启动命令时的工作目录。
step_cluster_plots_and_markers(load_project(Path(__file__).resolve().parents[1]))
