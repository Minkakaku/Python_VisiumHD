# Python_VisiumHD

跨项目复用的 Python 公共分析代码库。仓库管理公共方法、默认配置、环境和项目示例；实际项目在 `projects/<项目名>/` 使用这些方法，数据和运行产物不提交 Git。

## 目录

```text
src/python/       公共分析方法与配置加载器
config/            公共默认参数与参考配置
envs/              共享环境定义
projects/project1/ 示例项目（未运行）
  config/          项目覆盖、样本表、QC/注释等人工决策
  workflow/        编号短脚本与分阶段 runner
  data/processed/  中间对象
  data/reference/  项目参考资源
  results/         按分析点保存图表与中文说明
  logs/            运行记录
tests/             路径和项目配置检查
```

## 新项目

在仓库根目录执行：

```bash
conda env create -f envs/environment.yaml
python -m pip install -e .
cp -R projects/project1 projects/my_project
```

修改 `projects/my_project/config/config.yaml` 和样本/注释表，准备输入对象与参考文件。示例中的物种、reference、root、分组和文件名必须按真实项目核对，不能直接作为生物学结论。

每个编号脚本可单独执行；脚本从自身位置定位项目，因此复制后不会继续写入 project1。Python 包项目先完成 editable install。

```bash
python projects/my_project/workflow/00_setup.py
```

可用步骤：
- `00_setup.py`
- `01_load.py`
- `02_qc.py`
- `03_filter.py`
- `04_cluster_cpu.py`
- `04_cluster_gpu.py`
- `05_cluster_markers.py`
- `06_annotation.py`
- `07_spatial_statistics.py`
- `08_cell_neighborhood.py`
- `09_proximity.py`
- `10_nearest_neighbor.py`

`run_all.sh` 对主流程只运行到 QC 检查点；独立分析方法需显式选择步骤，避免把不同候选方法混为同一条流程。已运行脚本不改写，重跑新增编号并使用新输出目录。

## 配置与路径

公共 `config/default.yaml` → `config/reference.yaml` → 项目 `config/config.yaml`，后者覆盖前者。嵌套映射按键合并，列表和标量整体替换。共享方法接收配置/参数，不包含实际项目路径与样本名。项目相对输入输出路径以项目目录为基准；参考文件路径写在配置中。

## 结果与解释

最终图保存同名 PNG/PDF，配套 CSV/TSV 和 `README_results_CN.md` 放在分析点结果目录。中间对象放在 `data/processed/`。细胞或空间 bins 不作为独立生物学重复；条件比较用样本级汇总。方法特有约束见 [AGENTS.md](AGENTS.md) 和公共代码注释。

## 验证

```bash
python tests/test_project_layout.py
```

此检查验证配置、路径或安全导入，不运行真实组学分析。完整分析需要对应生信依赖、真实输入数据和人工 QC/注释。

## 可选 Bin2Cell 与邻近统计

`envs/bin2cell.yaml` 是可选环境定义；现有模板仍默认读取 8 µm bins。实际 2 µm Bin2Cell 项目的专用 workflow、注释和数据保留在项目内。

`proximity_sample_summary_fast` 按样本独立使用 KD 树累计计数，坐标须为 µm；半径内先逐细胞计算邻居类型比例再取来源类型均值，无邻居计零。环带为 `(low, high]`，报告逐细胞计数的中位数，不按面积归一化。可返回逐细胞比例用于描述性图；输出列沿用 Xenium 邻近统计口径，区别于原有 `proximity_sample_summary` 的列名。`tests/test_proximity_fast.py` 用直接欧氏距离检查边界和隔离。
