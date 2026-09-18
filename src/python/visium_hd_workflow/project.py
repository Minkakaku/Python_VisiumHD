"""加载公共配置，并用项目配置覆盖同名参数。"""

from pathlib import Path
from types import SimpleNamespace

import yaml


# 加载公共与项目 YAML，按顶层键覆盖配置并建立该项目的输入输出路径。
def load_project(project_dir):
    project_dir = Path(project_dir).resolve()
    workspace_dir = project_dir.parents[1]
    # 后读取的同名顶层键覆盖前值；嵌套字典在此整体替换，不递归合并。
    settings = {}
    for path in (
        workspace_dir / "config" / "default.yaml",
        workspace_dir / "config" / "reference.yaml",
        project_dir / "config" / "config.yaml",
    ):
        settings.update(yaml.safe_load(path.read_text(encoding="utf-8")) or {})

    # 配置键转为大写属性，以兼容各分析步骤使用的 project.PARAMETER 接口。
    project = SimpleNamespace(**{key.upper(): value for key, value in settings.items()})
    project.WORKSPACE_DIR = workspace_dir
    project.PROJECT_DIR = project_dir
    project.CONFIG_DIR = project_dir / "config"
    project.DATA_DIR = project_dir / "data"
    project.RAW_DIR = project.DATA_DIR / "raw"
    project.PROCESSED_DIR = project.DATA_DIR / "processed"
    project.REFERENCE_DIR = project.DATA_DIR / "reference"
    project.RESULTS_DIR = project_dir / "results"
    project.LOGS_DIR = project_dir / "logs"
    project.SAMPLE_INFO = project.CONFIG_DIR / "samples.tsv"
    project.QC_THRESHOLDS = project.CONFIG_DIR / "qc_thresholds.tsv"
    project.ANNOTATION_FILE = project.CONFIG_DIR / "cluster_annotation.tsv"
    project.SPACERANGER_CONFIG = project.CONFIG_DIR / "spaceranger.yaml"
    project.ADATA_RAW = project.PROCESSED_DIR / "00_adata_raw.h5ad"
    project.ADATA_QC = project.PROCESSED_DIR / "01_adata_qc.h5ad"
    project.ADATA_FILTERED = project.PROCESSED_DIR / "02_adata_filtered.h5ad"
    project.ADATA_CLUSTERED = project.PROCESSED_DIR / "03_adata_clustered.h5ad"
    project.ADATA_ANNOTATED = project.PROCESSED_DIR / "04_adata_annotated.h5ad"
    return project
