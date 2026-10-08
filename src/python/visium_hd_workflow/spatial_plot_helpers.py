"""项目中提取的空间绘图工具；坐标、细胞标签和颜色由调用者提供。"""
from __future__ import annotations
from collections import defaultdict
from pathlib import Path
import numpy as np
from matplotlib.path import Path as PlotPath
from PIL import Image
from . import runtime  # 导入公有 NPG 默认设置；项目仍可显式覆盖颜色。


# 判断细胞中心是否位于给定 ROI 内，不改变原始坐标。
def inside_roi(centers: np.ndarray, bounds: tuple[float, float, float, float]) -> np.ndarray:
    x_min, x_max, y_min, y_max = bounds
    return (
        (centers[:, 0] >= x_min)
        & (centers[:, 0] < x_max)
        & (centers[:, 1] >= y_min)
        & (centers[:, 1] < y_max)
    )

# 抵消相邻方格的共享边，仅保留真实轮廓边。
def _toggle_edge(edges: dict, start: tuple[int, int], stop: tuple[int, int]) -> None:
    key = (start, stop) if start < stop else (stop, start)
    if key in edges:
        del edges[key]
    else:
        edges[key] = (start, stop)

# 描绘方格并集的外环和内孔，保留分割轮廓。
def polyomino_loops(rows: np.ndarray, cols: np.ndarray) -> list[np.ndarray]:
    """Trace exact loops around a union of square array bins in grid units."""
    rows, cols = np.asarray(rows), np.asarray(cols)
    if rows.ndim != 1 or rows.shape != cols.shape or not np.isfinite(np.r_[rows, cols]).all():
        raise ValueError("Provide matching finite row/column arrays")
    if not (np.equal(rows, np.floor(rows)).all() and np.equal(cols, np.floor(cols)).all()):
        raise ValueError("Array bins must have integer coordinates")
    if len(np.unique(np.column_stack((rows, cols)), axis=0)) != len(rows):
        raise ValueError("Duplicate array bins would cancel their boundaries")
    edges = {}
    for row, col in zip(rows.tolist(), cols.tolist(), strict=True):
        top_left = (2 * col - 1, 2 * row - 1)
        top_right = (2 * col + 1, 2 * row - 1)
        bottom_right = (2 * col + 1, 2 * row + 1)
        bottom_left = (2 * col - 1, 2 * row + 1)
        _toggle_edge(edges, top_left, top_right)
        _toggle_edge(edges, top_right, bottom_right)
        _toggle_edge(edges, bottom_right, bottom_left)
        _toggle_edge(edges, bottom_left, top_left)

    outgoing = defaultdict(list)
    for start, stop in edges.values():
        outgoing[start].append(stop)
    loops = []
    remaining = sum(len(stops) for stops in outgoing.values())
    while remaining:
        start = next(vertex for vertex, stops in outgoing.items() if stops)
        current = start
        loop = [start]
        maximum = remaining + 1
        for _ in range(maximum):
            if not outgoing[current]:
                break
            stop = outgoing[current].pop()
            remaining -= 1
            loop.append(stop)
            current = stop
            if current == start:
                break
        if len(loop) >= 4 and loop[-1] == start:
            loops.append(np.asarray(loop, dtype=np.float64) / 2.0)
    return loops

# 拟合阵列坐标到 H&E 微米坐标，并返回配准误差。
def fit_array_to_he_microns(
    rows: np.ndarray,
    cols: np.ndarray,
    spatial_cropped: np.ndarray,
    scalefactor: float,
    mpp: float,
    max_bins: int,
) -> tuple[np.ndarray, dict]:
    """Fit array column/row coordinates to the rendered H&E micron coordinates."""
    rows, cols, spatial_cropped = np.asarray(rows), np.asarray(cols), np.asarray(spatial_cropped)
    if rows.ndim != 1 or rows.shape != cols.shape or spatial_cropped.shape != (len(rows), 2):
        raise ValueError("Array coordinates and H&E xy coordinates must align")
    if len(rows) < 3 or max_bins < 3 or scalefactor <= 0 or mpp <= 0:
        raise ValueError("Affine registration requires at least three bins and positive scales")
    if not np.isfinite(np.r_[rows, cols, spatial_cropped.ravel(), scalefactor, mpp]).all():
        raise ValueError("Registration coordinates and scales must be finite")
    fit_count = min(rows.size, int(max_bins))
    positions = np.linspace(0, rows.size - 1, fit_count, dtype=np.int64)
    design = np.column_stack(
        [cols[positions].astype(np.float64), rows[positions].astype(np.float64), np.ones(fit_count)]
    )
    if np.linalg.matrix_rank(design) < 3:
        raise ValueError("Registration bins must span both spatial dimensions")
    targets = spatial_cropped[positions].astype(np.float64) * float(scalefactor) * float(mpp)
    coefficients, _, _, _ = np.linalg.lstsq(design, targets, rcond=None)
    residual = design @ coefficients - targets
    distances = np.sqrt(np.square(residual).sum(axis=1))
    diagnostics = {
        "affine_fit_bins": int(fit_count),
        "affine_rmse_um": float(np.sqrt(np.mean(np.square(distances)))),
        "affine_max_error_um": float(distances.max()),
        "array_col_to_he_x": float(coefficients[0, 0]),
        "array_row_to_he_x": float(coefficients[1, 0]),
        "he_x_intercept": float(coefficients[2, 0]),
        "array_col_to_he_y": float(coefficients[0, 1]),
        "array_row_to_he_y": float(coefficients[1, 1]),
        "he_y_intercept": float(coefficients[2, 1]),
    }
    return coefficients, diagnostics

# 每个细胞的外环与孔洞共用一个矢量 Path，填色与轮廓使用相同几何。
def compound_paths(polygons, polygon_cells):
    """同一细胞的外环/内孔组成一个 Path；填色与轮廓严格使用同一几何。"""
    grouped = {}
    for polygon, cell in zip(polygons, polygon_cells, strict=True):
        grouped.setdefault(int(cell), []).append(PlotPath(polygon, closed=True))
    cells = np.asarray(sorted(grouped), dtype=int)
    return [PlotPath.make_compound_path(*grouped[int(cell)]) for cell in cells], cells

# 按位置随机抽取非空等面积窗口，不按基因或细胞类型选窗。
def random_fields(centers, bounds, side, count, seed, margin=0, stratified=False, padding=0):
    """在 ROI 内不重叠的等面积网格中等概率抽样；不读取基因或细胞类型。

    仅剔除没有任何分析保留细胞的空窗口，不增加 tissue_mask，不按密度加权。
    返回坐标单位为 µm；同一样本的两个基因共用该结果。
    新版可限制边缘留白，并按四象限各抽一窗，中心间距至少为窗宽的两倍。
    固定 seed 后，可重复生成相同的随机视野。
    padding 仅用于阳性示例：要求传入的阳性细胞中心远离局部窗边界。
    """
    centers = np.asarray(centers, dtype=float)
    if centers.ndim != 2 or centers.shape[1] != 2 or not np.isfinite(centers).all():
        raise ValueError("Provide finite xy cell centers")
    if not np.isfinite([*bounds, side, count, margin, padding]).all() or side <= 0 or count < 1 or int(count) != count:
        raise ValueError("Field size and integer count must be positive and finite")
    if margin < 0 or padding < 0 or padding >= side / 2:
        raise ValueError("Margins must be nonnegative and padding smaller than half a field")
    count = int(count)
    x0, x1, y0, y1 = bounds
    nx, ny = int((x1 - x0) // side), int((y1 - y0) // side)
    if min(nx, ny) < 1:
        raise ValueError("Field size exceeds ROI")
    left, top = x0 + ((x1 - x0) - nx * side) / 2, y0 + ((y1 - y0) - ny * side) / 2
    candidates = []
    for iy in range(ny):
        for ix in range(nx):
            box = (left + ix * side, left + (ix + 1) * side,
                   top + iy * side, top + (iy + 1) * side)
            if min(box[0] - x0, x1 - box[1], box[2] - y0, y1 - box[3]) < margin:
                continue
            interior = (box[0] + padding, box[1] - padding, box[2] + padding, box[3] - padding)
            if inside_roi(centers, interior).any():
                candidates.append(box)
    if len(candidates) < count:
        raise ValueError("Insufficient nonempty, nonoverlapping random fields")
    rng = np.random.default_rng(seed)
    if stratified:
        if count != 4:
            raise ValueError("Quadrant-stratified fields require exactly four fields")
        middle = np.array([[(b[0] + b[1]) / 2, (b[2] + b[3]) / 2] for b in candidates])
        quadrant = (middle[:, 0] >= (x0 + x1) / 2).astype(int) + 2 * (middle[:, 1] >= (y0 + y1) / 2)
        groups = [np.flatnonzero(quadrant == q) for q in range(4)]
        if any(len(group) == 0 for group in groups):
            raise ValueError("A spatial quadrant has no eligible field")
        for _ in range(1000):
            picks = [int(rng.choice(group)) for group in groups]
            locations = middle[picks]
            distances = np.linalg.norm(locations[:, None] - locations[None, :], axis=2)
            if (distances[np.triu_indices(4, 1)] >= 2 * side).all():
                break
        else:
            raise ValueError("Cannot sample four spatially separated interior fields")
    else:
        picks = rng.choice(len(candidates), size=count, replace=False)
    return [candidates[int(i)] for i in picks], len(candidates)

# 校验用户指定的等面积窗口，不改变分析 ROI。
def manual_fields(centers, bounds, side, count):
    """以用户指定中心生成等大方框；只修改展示视野，不改变分析 ROI。"""
    if count < 1 or int(count) != count or not np.isfinite(bounds).all():
        raise ValueError("Provide finite bounds and a positive integer count")
    centers = np.asarray(centers, dtype=float)
    if centers.shape != (count, 2) or not np.isfinite(centers).all():
        raise ValueError("Manual fields require one finite xy center per field")
    if not np.isfinite(side) or side <= 0 or len(np.unique(centers, axis=0)) != count:
        raise ValueError("Manual field size must be positive and centers distinct")
    boxes = np.column_stack((centers[:, 0] - side / 2, centers[:, 0] + side / 2,
                             centers[:, 1] - side / 2, centers[:, 1] + side / 2))
    x0, x1, y0, y1 = bounds
    if not ((boxes[:, 0] >= x0) & (boxes[:, 1] <= x1)
            & (boxes[:, 2] >= y0) & (boxes[:, 3] <= y1)).all():
        raise ValueError("A manual field extends outside the approved ROI")
    return boxes.tolist()

# 按微米范围裁剪 H&E，返回与空间坐标一致的图像 extent。
def load_he_roi(
    path: Path,
    bounds: tuple[float, float, float, float],
    mpp: float,
    max_pixels: int,
    white_mix: float,
) -> tuple[np.ndarray, tuple[float, float, float, float]]:
    if not np.isfinite([*bounds, mpp, white_mix]).all() or mpp <= 0 or max_pixels < 1 or not 0 <= white_mix <= 1:
        raise ValueError("Provide positive image scales and white_mix between zero and one")
    x_min, x_max, y_min, y_max = bounds
    with Image.open(path) as source:
        left = max(0, int(np.floor(x_min / mpp)))
        right = min(source.width, int(np.ceil(x_max / mpp)))
        top = max(0, int(np.floor(y_min / mpp)))
        bottom = min(source.height, int(np.ceil(y_max / mpp)))
        if left >= right or top >= bottom:
            raise ValueError(f"Empty H&E crop for {path}: {(left, top, right, bottom)}")
        image = source.crop((left, top, right, bottom)).convert("RGB")
        image.thumbnail((int(max_pixels), int(max_pixels)), Image.Resampling.LANCZOS)
        array = np.asarray(image, dtype=np.float32) / 255.0
    mix = float(white_mix)
    extent = (left * mpp, right * mpp, bottom * mpp, top * mpp)
    return np.clip(array * (1.0 - mix) + mix, 0.0, 1.0), extent

# 为给定空间边界增加微米比例尺。
def add_scalebar(
    axis,
    bounds: tuple[float, float, float, float],
    length_um: float,
    font_size: float,
) -> None:
    x_min, x_max, y_min, y_max = bounds
    width = x_max - x_min
    height = y_max - y_min
    if not np.isfinite([*bounds, length_um, font_size]).all() or height <= 0 or not 0 < length_um <= width * 0.9 or font_size <= 0:
        raise ValueError("Scale bar must fit the finite spatial bounds")
    x0 = x_min + 0.055 * width
    y0 = y_max - 0.065 * height
    axis.plot([x0, x0 + length_um], [y0, y0], color="white", linewidth=5.2, solid_capstyle="butt", zorder=6)
    axis.plot([x0, x0 + length_um], [y0, y0], color="#111111", linewidth=2.6, solid_capstyle="butt", zorder=7)
    axis.text(
        x0 + length_um / 2.0,
        y0 - 0.018 * height,
        f"{length_um:g} µm",
        ha="center",
        va="bottom",
        fontsize=font_size,
        color="#111111",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 1.0},
        zorder=8,
    )

# 合并所有面板表达值计算共享上限，全部为零时使用 1。
def gene_vmax(values: list[np.ndarray], quantile: float) -> float:
    if not values or not np.isfinite(quantile) or not 0 <= quantile <= 1:
        raise ValueError("Provide expression arrays and a quantile between zero and one")
    pooled = np.concatenate(values)
    if not pooled.size or not np.isfinite(pooled).all() or (pooled < 0).any():
        raise ValueError("Expression must be finite, nonnegative and nonempty")
    vmax = float(np.quantile(pooled, quantile))
    if vmax <= 0:
        positive = pooled[pooled > 0]
        vmax = float(np.quantile(positive, quantile)) if positive.size else 1.0
    return vmax

# 输出带标签和隐藏标签两版，固定相同画布，并恢复原图的可见性。
def save_label_variants(fig, stem, text_artists, dpi=600):
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    artists = list(text_artists)
    visibility = [artist.get_visible() for artist in artists]
    fig.canvas.draw()
    bounds = fig.get_tightbbox(fig.canvas.get_renderer()).padded(0.1)
    stems = []
    try:
        for labelled in (True, False):
            for artist, visible in zip(artists, visibility):
                artist.set_visible(visible if labelled else False)
            output = stem.with_name(stem.name + ("_with_labels" if labelled else "_without_labels"))
            fig.savefig(str(output) + ".png", dpi=dpi, facecolor="white", bbox_inches=bounds)
            fig.savefig(str(output) + ".pdf", facecolor="white", bbox_inches=bounds)
            stems.append(output)
    finally:
        for artist, visible in zip(artists, visibility):
            artist.set_visible(visible)
    return stems
