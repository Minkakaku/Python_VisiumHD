from visium_hd_workflow.project import load_project


# 构造临时项目，验证项目参数覆盖公共默认，同时保留其他默认值和路径。
def test_project_config_overrides_public_defaults(tmp_path):
    workspace = tmp_path / "workspace"
    project_dir = workspace / "projects" / "project1"
    (workspace / "config").mkdir(parents=True)
    (project_dir / "config").mkdir(parents=True)
    (workspace / "config" / "default.yaml").write_text("bin_size: 8\nrandom_seed: 10\n")
    (workspace / "config" / "reference.yaml").write_text("organism: Human\n")
    (project_dir / "config" / "config.yaml").write_text("bin_size: 2\n")

    project = load_project(project_dir)

    assert project.BIN_SIZE == 2
    assert project.RANDOM_SEED == 10
    assert project.ORGANISM == "Human"
    assert project.SAMPLE_INFO == project_dir / "config" / "samples.tsv"
