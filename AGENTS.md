# Python Visium HD 项目代码约定

- 正式平台是 Python/Scanpy Visium HD，默认读取 Space Ranger 的 8 µm square-bin 输出。
- 跨项目复用函数放在 `src/python/`；公共默认放在根目录 `config/`；样本、阈值和人工注释放在 `projects/<项目>/config/`。
- 每个项目的 `workflow/` 只放按执行顺序编号的短入口；已运行脚本不覆盖，重跑或改参时新增下一编号脚本。
- 正式入口只用 `.py`；原始 notebook 和脚本逐字节保存在 `legacy/original/`，不得直接修改。
- `visium_hd_workflow.cluster_gpu` 是独立 RAPIDS/CuPy/RMM 聚类入口，CPU 主流程不强制安装 GPU 依赖。
- 不在 `src/python/` 写项目编号、样本名、绝对路径、固定 cluster 注释或固定 ROI。
- 多样本对象使用 `sample:barcode` 唯一索引；空间邻接、距离和 CN 必须按样本分别计算。
- 组间统计先汇总到样本层面；Visium HD bins 不是独立生物学重复，也不等同于分割细胞。
- 中间 H5AD 写入项目 `data/processed/`；最终图、表和中文说明写入项目对应的 `results/` 分析点目录。
- 正式图同时保存 PNG 和 PDF，绘图输入表保存 TSV；函数与关键算法块使用简洁中文注释。
- 优先保留原始分析意图，只修复明确错误、过时接口、项目硬编码和影响统计结论的问题。
