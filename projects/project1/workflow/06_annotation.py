# 按项目选择的注释方法或人工映射赋予细胞类型标签并保存对象。

from pathlib import Path

from visium_hd_workflow.project import load_project
from visium_hd_workflow.steps import step_annotation

# 由脚本所在位置定位项目目录，避免依赖启动命令时的工作目录。
step_annotation(load_project(Path(__file__).resolve().parents[1]))
