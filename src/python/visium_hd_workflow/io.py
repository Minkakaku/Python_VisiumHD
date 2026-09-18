import json
from pathlib import Path

import numpy as np
import pandas as pd


# 功能：兼容样本目录和直接指向 outs 的路径写法。
def resolve_outs(path):
    path = Path(path)
    return path / "outs" if (path / "outs").is_dir() else path


# 功能：读取单个 Visium HD square-bin 表达矩阵，并显式对齐 barcode 与空间坐标。
def read_binned_sample(spaceranger_dir, sample, group, bin_size=8):
    import scanpy as sc

    outs = resolve_outs(spaceranger_dir)
    bin_dir = outs / "binned_outputs" / f"square_{int(bin_size):03d}um"
    adata = sc.read_10x_h5(bin_dir / "filtered_feature_bc_matrix.h5")
    adata.var_names_make_unique()

    positions = pd.read_parquet(bin_dir / "spatial" / "tissue_positions.parquet")
    positions = positions.set_index("barcode")
    keep = adata.obs_names.intersection(positions.index, sort=False)
    adata = adata[keep].copy()
    positions = positions.reindex(adata.obs_names)

    for column in positions.columns:
        adata.obs[column] = positions[column].to_numpy()
    adata.obs["barcode"] = adata.obs_names.astype(str)
    adata.obs["sample"] = str(sample)
    adata.obs["group"] = str(group)

    # 图像像素坐标用于组织图叠加；用 Space Ranger 比例换成 µm 后再做距离统计。
    adata.obsm["spatial"] = positions[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(float)
    adata.obsm["spatial_um"] = adata.obsm["spatial"] * read_microns_per_pixel(outs)
    adata.layers["counts"] = adata.X.copy()
    adata.obs_names = pd.Index([f"{sample}:{x}" for x in adata.obs["barcode"]])
    adata.uns["visium_hd"] = {
        "sample": str(sample),
        "group": str(group),
        "bin_size_um": int(bin_size),
        "outs": str(outs),
        "analysis_unit": "square_bin",
    }
    return adata


# 功能：从 Space Ranger scalefactors 中读取像素到 µm 的换算；显式配置可覆盖该值。
def read_microns_per_pixel(outs):
    candidates = sorted((Path(outs) / "binned_outputs").glob("square_*um/spatial/scalefactors_json.json"))
    for path in candidates:
        value = json.loads(path.read_text(encoding="utf-8")).get("microns_per_pixel")
        if value is not None:
            return float(value)
    raise FileNotFoundError("未找到 microns_per_pixel；请在项目 samples.tsv 填写 microns_per_unit。")

# 功能：读取 Space Ranger v4 分割结果，并按 instance_id 对齐 polygon 质心和表达表。
def read_segmented_sample(spaceranger_dir, sample, group, table_key, shape_key, microns_per_unit=""):
    from spatialdata_io import visium_hd

    outs = resolve_outs(spaceranger_dir)
    sdata = visium_hd(
        outs,
        dataset_id=str(sample),
        load_segmentations_only=True,
        load_all_images=False,
    )
    adata = sdata.tables[str(table_key)].copy()
    shapes = sdata.shapes[str(shape_key)]
    instance_key = adata.uns.get("spatialdata_attrs", {}).get("instance_key", "instance_id")
    centroids = shapes.geometry.centroid
    centroid_table = pd.DataFrame(
        {"x": centroids.x.to_numpy(), "y": centroids.y.to_numpy()},
        index=shapes.index.astype(str),
    )
    instance_ids = adata.obs[instance_key].astype(str)
    coords = centroid_table.reindex(instance_ids).to_numpy(float)
    if np.isnan(coords).any():
        raise ValueError("分割 polygon 与表达表的 instance_id 未完全对齐。")

    adata.var_names_make_unique()
    adata.obs["barcode"] = adata.obs_names.astype(str)
    adata.obs["sample"] = str(sample)
    adata.obs["group"] = str(group)
    scale = float(microns_per_unit) if str(microns_per_unit).strip() else read_microns_per_pixel(outs)
    adata.obsm["spatial"] = coords
    adata.obsm["spatial_um"] = coords * scale
    adata.layers["counts"] = adata.X.copy()
    adata.obs_names = pd.Index([f"{sample}:{x}" for x in adata.obs["barcode"]])
    adata.uns["visium_hd"] = {
        "sample": str(sample),
        "group": str(group),
        "outs": str(outs),
        "analysis_unit": "segmented_cell",
        "table_key": str(table_key),
        "shape_key": str(shape_key),
        "microns_per_unit": scale,
    }
    return adata


# 功能：按样本表读取并合并多个切片，保证索引唯一且保留每个样本的空间元数据。
def load_visium_hd_data(sample_info, bin_size=8, base_dir=None):
    import anndata as ad

    sample_table = pd.read_csv(sample_info, sep="\t", dtype=str).fillna("")
    base_dir = Path(base_dir or ".")
    objects = []
    sample_metadata = {}
    for row in sample_table.itertuples(index=False):
        spaceranger_dir = Path(row.spaceranger_dir)
        if not spaceranger_dir.is_absolute():
            spaceranger_dir = base_dir / spaceranger_dir
        mode = (getattr(row, "input_mode", "") or "bins").lower()
        if mode == "segmented":
            obj = read_segmented_sample(
                spaceranger_dir,
                row.sample,
                row.group,
                row.table_key,
                row.shape_key,
                getattr(row, "microns_per_unit", ""),
            )
        else:
            obj = read_binned_sample(spaceranger_dir, row.sample, row.group, bin_size)
        sample_metadata[str(row.sample)] = obj.uns["visium_hd"]
        objects.append(obj)

    combined = ad.concat(objects, join="inner", merge="same", index_unique=None)
    combined.obs_names_make_unique()
    combined.uns["visium_hd_samples"] = sample_metadata
    combined.obs["sample"] = pd.Categorical(combined.obs["sample"], categories=sample_table["sample"])
    combined.obs["group"] = pd.Categorical(combined.obs["group"])
    return combined
