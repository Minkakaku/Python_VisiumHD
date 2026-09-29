"""Run with python tests/test_project_layout.py; requires only PyYAML."""
from pathlib import Path
import runpy
import tempfile
import shutil
import yaml

ROOT = Path(__file__).resolve().parents[1]
api = runpy.run_path(str(ROOT / "src/python/visium_hd_workflow/project.py"))
with tempfile.TemporaryDirectory() as tmp:
    workspace = Path(tmp).resolve() / "workspace"
    shutil.copytree(ROOT / "config", workspace / "config")
    for name in ("a", "b"):
        shutil.copytree(ROOT / "projects/project1/config", workspace / "projects" / name / "config")
    a, b = workspace / "projects/a", workspace / "projects/b"
    cfg = yaml.safe_load((a / "config/config.yaml").read_text()) or {}
    cfg["group_colors"] = {"CON": "#000000"}
    (a / "config/config.yaml").write_text(yaml.safe_dump(cfg))
    first, second = api["load_project"](a), api["load_project"](b)
    assert first.PROJECT_DIR == a and second.PROJECT_DIR == b
    assert first.RESULTS_DIR != second.RESULTS_DIR
    assert first.SAMPLE_INFO == a / "config/samples.tsv"
    assert first.GROUP_COLORS["CON"] == "#000000"
    assert first.GROUP_COLORS["EXP"] == second.GROUP_COLORS["EXP"]
    assert api["merge_config"]({"a": {"x": 1, "y": [1, 2]}}, {"a": {"y": [3]}}) == {"a": {"x": 1, "y": [3]}}
print("Project isolation and configuration override checks passed")
