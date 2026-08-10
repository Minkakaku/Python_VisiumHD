import os
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import gseapy as gp
import squidpy as sq
import seaborn as sns
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import MiniBatchKMeans
from scipy.spatial import KDTree
from scipy.stats import mannwhitneyu, ttest_ind
import multiprocessing
from spatialdata_io import xenium
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context

import warnings
warnings.filterwarnings("ignore")

plt.style.use('default')
plt.rcParams['figure.facecolor'] = 'white'
sc.settings.set_figure_params(dpi=100, dpi_save=200, figsize=(5, 5), facecolor='white')


os.environ['OPENBLAS_NUM_THREADS'] = '0'
# 推荐配置：根据可用CPU核心数设置
num_cores = os.cpu_count()  # 获取系统CPU核心数
optimal_threads = max(1, num_cores - 2)  # 保留2个核心给系统
# 设置OpenBLAS线程数（推荐方式）
# os.environ['OPENBLAS_NUM_THREADS'] = 1
os.environ['OPENBLAS_NUM_THREADS'] = str(optimal_threads)


data_dir = "/home/hf/VH2604022693-01"
sample_info = '/home/hf/VH2604022693-01/sample_info.txt'
output_dir = '/home/hf/VH2604022693-01/OUT/'
bin_size = '008'

if not os.path.exists(output_dir):
    os.makedirs(os.path.join(output_dir, 'data'), exist_ok=True)
    os.makedirs(os.path.join(output_dir, 'report'), exist_ok=True)

for i in ['01-Data_QC', 
          '02-Cell_filter', 
          '03-Dimension_reduction_cluster', 
          '04-Marker_gene_analysis', 
          '05-Cell_type_annotation', 
          '06-Computation_spatial_statistics',
          '07-Cell_Neighborhood_analysis',
          '08-Proximity_analysis'
          ]:
    os.makedirs(os.path.join(output_dir, 'report', i), exist_ok=True)

def cell_neighbors(adata, column, n_neighbors=10, k=10, n_clusters=10, out_dir=None):
    spatial_coords = adata.obsm['spatial']
    onehot_encoding = pd.get_dummies(adata.obs[column])
    
    # 计算最近邻居
    nbrs = NearestNeighbors(n_neighbors=n_neighbors).fit(spatial_coords)
    distances, indices = nbrs.kneighbors(spatial_coords)

    # 对每个cell，计算其邻居的类型分布
    cluster_cols = adata.obs[column].unique()
    values = onehot_encoding[cluster_cols].values

    def compute_window_sums(indices, values, k):
        windows = []
        for idx in range(indices.shape[0]):
            neighbors = indices[idx, :k]
            window_sum = values[neighbors].sum(axis=0)
            windows.append(window_sum)
        return np.array(windows)
    
    windows = compute_window_sums(indices, values, k)
    km = MiniBatchKMeans(n_clusters=n_clusters, random_state=0)
    labels = km.fit_predict(windows)
    adata.obs[f'CNs_{n_clusters}'] = [f'CN{i}' for i in labels]

    adata.obs.to_csv(os.path.join(out_dir, f'CNs_{n_clusters}_all_cell_features.xls'), sep='\t', index=None)

    # 保存聚类中心
    k_centroids = km.cluster_centers_
    
    # 计算Fold Change (fc)
    tissue_avgs = values.mean(axis=0)
    fc = np.log2(((k_centroids + tissue_avgs) / (k_centroids + tissue_avgs).sum(axis=1, keepdims=True)) / tissue_avgs)  
    fc = pd.DataFrame(fc, columns=adata.obs[column].unique())
    fc.index = 'CN' + fc.index.astype(str)
    
    plt.figure(figsize=(20, 8))
    sns.set_style("white")
    g = sns.clustermap(
        fc,
        vmin=-2,
        vmax=2,
        cmap="coolwarm",
        row_cluster=False,
        col_cluster=False,   # ← 关键！！
        linewidths=0,
        figsize=(0.6 * fc.shape[1] + 4, 6)
    )

    # 强制显示全部列标签（现在不会乱）
    g.ax_heatmap.set_xticks(np.arange(fc.shape[1]) + 0.5)
    g.ax_heatmap.set_xticklabels(
        fc.columns,
        rotation=90,
        ha="center",
        fontsize=8
    )

    g.ax_heatmap.set_xlabel("Cell type")
    g.ax_heatmap.set_ylabel("Cell Neighborhood (CN)")

    plt.savefig(
        os.path.join(out_dir, f"CNs_{n_clusters}_enrichment_full_labels.png"),
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()
    
    # 统计柱状图
    # 计算每个样本中每个聚类的细胞数量
    cluster_counts = adata.obs.groupby(['sample', f'CNs_{n_clusters}']).size().unstack(fill_value=0)

    # 计算每个样本中每个聚类的细胞占比
    cluster_proportions = cluster_counts.div(cluster_counts.sum(axis=1), axis=0)

    # 转换为长格式以便绘图
    cluster_proportions_long = cluster_proportions.reset_index().melt(id_vars='sample', var_name='Cluster', value_name='Proportion')

    # 添加组信息
    cluster_proportions_long = cluster_proportions_long.merge(adata.obs[['sample', 'group']].drop_duplicates(), on='sample')

    # 绘制箱线图
    plt.figure(figsize=(12, 6))
    sns.boxplot(x='Cluster', y='Proportion', hue='group', data=cluster_proportions_long)
    sns.stripplot(x='Cluster', y='Proportion', hue='group', data=cluster_proportions_long, dodge=True, jitter=False, color='gray', alpha=0.3, legend=False)
    plt.savefig(os.path.join(out_dir, f"CNs_{n_clusters}_cell_proportion.png"), dpi=300, bbox_inches='tight')

    library_names = adata.obs['sample'].drop_duplicates()
    for library in library_names:
        with rc_context({'figure.figsize': (15, 8)}):
            sq.pl.spatial_scatter(
                adata[adata.obs['sample'] == library,:],
                library_id="spatial",
                color=[f'CNs_{n_clusters}'],
                title=f'{library} CNs',
                shape=None,
                # size=1,
                img=False,
                use_raw=True
            )
            plt.savefig(os.path.join(out_dir, f'{library}_spatial_CNs_{n_clusters}.png'), facecolor='white', bbox_inches="tight", dpi=300)

    return adata



adata = sc.read_h5ad(os.path.join(output_dir, 'data', 'adata_celltype.h5ad'))

adata = cell_neighbors(
    adata,
    'celltype', 
    n_neighbors=10, 
    k=10, 
    n_clusters=10, 
    out_dir=os.path.join(output_dir, 'report', '07-Cell_Neighborhood_analysis')
)

adata = cell_neighbors(
    adata,
    'celltype', 
    n_neighbors=10, 
    k=10, 
    n_clusters=15, 
    out_dir=os.path.join(output_dir, 'report', '07-Cell_Neighborhood_analysis')
)