# 使用人工检查后的逐样本阈值过滤，保存后续聚类输入。

from pathlib import Path

from visium_hd_workflow.project import load_project
from visium_hd_workflow.steps import step_filter

# 由脚本所在位置定位项目目录，避免依赖启动命令时的工作目录。
step_filter(load_project(Path(__file__).resolve().parents[1]))
