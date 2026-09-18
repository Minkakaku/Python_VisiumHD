# project1

这是一个可复制的新项目入口。样本、阈值和人工注释放在 `config/`，短入口放在 `workflow/`，中间 H5AD 写入 `data/processed/`，最终图表写入 `results/`。

在仓库根目录安装公共环境和包：

```bash
conda env create -f envs/environment.yaml
conda activate visium-hd-workflow
pip install -e .
```

修改 `config/samples.tsv` 和 `config/config.yaml` 后，从仓库根目录逐步运行：

```bash
python projects/project1/workflow/00_setup.py
python projects/project1/workflow/01_load.py
python projects/project1/workflow/02_qc.py
```

QC 后人工确认 `config/qc_thresholds.tsv`，再运行 `03_filter.py`。marker 完成后人工确认 `config/cluster_annotation.tsv`，再运行注释及空间分析。已运行的 workflow 文件不要覆盖；参数改变时新增下一编号脚本。
