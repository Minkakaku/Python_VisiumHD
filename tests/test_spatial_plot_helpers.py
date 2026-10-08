"""用合成几何和图形核对提取的方法，不读取项目数据或运行分割。"""
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import pytest
from visium_hd_workflow.spatial_plot_helpers import (
    polyomino_loops, compound_paths, fit_array_to_he_microns, random_fields,
    manual_fields, load_he_roi, gene_vmax, save_label_variants,
)


# 覆盖孔洞、配准、盲抽样、统一色阶和两版导出的一条可重跑检查。
def test_project_plotting_geometry_and_exports(tmp_path):
    bins = np.array([(r, c) for r in range(3) for c in range(3) if (r, c) != (1, 1)])
    loops = polyomino_loops(bins[:, 0], bins[:, 1])
    areas = [abs(np.dot(p[:-1, 0], p[1:, 1]) - np.dot(p[1:, 0], p[:-1, 1])) / 2 for p in loops]
    assert sorted(areas) == [1, 9]
    paths, cells = compound_paths(loops, [7] * len(loops))
    assert len(paths) == 1 and cells.tolist() == [7]
    with pytest.raises(ValueError, match="Duplicate"):
        polyomino_loops(np.array([0, 0]), np.array([0, 0]))
    rows, cols = np.array([0, 0, 1, 1]), np.array([0, 1, 0, 1])
    xy = np.column_stack((cols * 2 + 3, rows * 4 + 5))
    affine, audit = fit_array_to_he_microns(rows, cols, xy, 1, 1, 4)
    assert np.allclose(np.column_stack((cols, rows, np.ones(4))) @ affine, xy)
    assert audit["affine_max_error_um"] < 1e-12
    centers = np.array([(x, y) for x in range(5, 100, 10) for y in range(5, 100, 10)])
    bounds = (0, 100, 0, 100)
    first, eligible = random_fields(centers, bounds, 20, 4, 17)
    assert first == random_fields(centers, bounds, 20, 4, 17)[0] and eligible == 25
    assert len({tuple(box) for box in first}) == 4
    assert manual_fields([[50, 50]], bounds, 20, 1) == [[40, 60, 40, 60]]
    with pytest.raises(ValueError, match="outside"):
        manual_fields([[0, 0]], bounds, 20, 1)
    image_path = tmp_path / "he.png"
    Image.new("RGB", (10, 10), "black").save(image_path)
    image, extent = load_he_roi(image_path, (2, 8, 2, 8), 1, 3, 0.5)
    assert image.shape == (3, 3, 3) and extent == (2, 8, 8, 2) and np.allclose(image, 0.5)
    assert gene_vmax([np.zeros(3)], 0.99) == 1
    assert gene_vmax([np.array([0, 1]), np.array([10])], 1) == 10
    fig, ax = plt.subplots()
    artist = ax.set_title("sample label")
    try:
        stems = save_label_variants(fig, tmp_path / "panel", [artist], dpi=50)
        assert artist.get_visible()
        assert all(s.with_suffix(".png").exists() and s.with_suffix(".pdf").exists() for s in stems)
        assert Image.open(stems[0].with_suffix(".png")).size == Image.open(stems[1].with_suffix(".png")).size
    finally:
        plt.close(fig)
