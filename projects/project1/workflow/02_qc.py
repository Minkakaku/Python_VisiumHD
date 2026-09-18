# 计算过滤前 QC 并生成候选阈值表，供人工检查。

from pathlib import Path

from visium_hd_workflow.project import load_project
from visium_hd_workflow.steps import step_qc

# 由脚本所在位置定位项目目录，避免依赖启动命令时的工作目录。
step_qc(load_project(Path(__file__).resolve().parents[1]))
