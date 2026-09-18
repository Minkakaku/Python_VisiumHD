# 加载项目配置并创建数据、结果和日志目录。

from pathlib import Path

from visium_hd_workflow.project import load_project
from visium_hd_workflow.runtime import prepare_project_dirs

# 由脚本所在位置定位项目目录，避免依赖启动命令时的工作目录。
project = load_project(Path(__file__).resolve().parents[1])
prepare_project_dirs(project)
print(f"Project ready: {project.PROJECT_DIR}")
