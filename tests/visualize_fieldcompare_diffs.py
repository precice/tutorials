#!/usr/bin/env python3
"""Render fieldcompare VTK diff fields as PNG images."""

from __future__ import annotations

import argparse
import re
import sys
from collections.abc import Iterator
from pathlib import Path

import numpy as np
import pyvista as pv


SUPPORTED_SUFFIXES = {".vtk", ".vtp", ".vtu"}
WINDOW_SIZE = (1024, 768)


def _scalar_values(values: np.ndarray) -> np.ndarray | None:
    """Return scalar values, using the magnitude for vectors and tensors."""
    array = np.asarray(values)
    if not np.issubdtype(array.dtype, np.number):
        return None
    if array.ndim == 1:
        return array
    if array.ndim == 2:
        return np.linalg.norm(array, axis=1)
    return None


def _fields(
    dataset: pv.DataSet,
) -> Iterator[tuple[str, np.ndarray, np.ndarray]]:
    """Yield field names, point locations, and scalar values."""
    locations = np.asarray(dataset.points)
    for field_name in dataset.point_data.keys():
        values = _scalar_values(np.asarray(dataset.point_data[field_name]))
        if values is None or len(values) != len(locations):
            continue
        finite = np.isfinite(values)
        if finite.any():
            yield field_name, locations[finite], values[finite]


def _glyph_radius(points: np.ndarray) -> float:
    """Return a radius based on representative nearest-neighbor distances."""
    if len(points) < 2:
        return 1.0

    extent = float(np.max(np.ptp(points, axis=0)))
    if extent <= 0:
        return 1.0

    sample = points[np.linspace(0, len(points) - 1, min(len(points), 64), dtype=int)]
    nearest_distances = []
    for point in sample:
        distances = np.linalg.norm(points - point, axis=1)
        distances = distances[distances > extent * 1e-12]
        if len(distances):
            nearest_distances.append(np.min(distances))
    return 0.2 * float(np.median(nearest_distances)) if nearest_distances else 1.0


def _set_camera(plotter: pv.Plotter, points: np.ndarray) -> None:
    """Use a face-on view for planar data and an isometric view otherwise."""
    extents = np.ptp(points, axis=0)
    max_extent = float(np.max(extents))
    flat_axis = int(np.argmin(extents))
    if max_extent > 0 and extents[flat_axis] <= max_extent * 1e-6:
        (plotter.view_yz, plotter.view_xz, plotter.view_xy)[flat_axis]()
    else:
        plotter.view_isometric()
    plotter.reset_camera()


def render_field(
    source_file: Path,
    output_file: Path,
    field_name: str,
    points: np.ndarray,
    values: np.ndarray,
) -> None:
    """Render one field using sphere glyphs colored by its diff values."""
    point_cloud = pv.PolyData(points)
    scalar_name = "difference"
    point_cloud.point_data[scalar_name] = values
    sphere = pv.Sphere(
        radius=_glyph_radius(points),
        theta_resolution=8,
        phi_resolution=8,
    )
    glyphs = point_cloud.glyph(orient=False, scale=False, geom=sphere)

    output_file.parent.mkdir(parents=True, exist_ok=True)
    vmin = float(np.min(values))
    vmax = float(np.max(values))
    plotter = pv.Plotter(off_screen=True, window_size=WINDOW_SIZE)
    try:
        plotter.set_background("white")
        plotter.add_mesh(
            glyphs,
            scalars=scalar_name,
            cmap="coolwarm",
            clim=(vmin, vmax),
            scalar_bar_args={"title": field_name},
        )
        plotter.add_text(
            (
                f"{source_file.name}\n"
                f"point field: {field_name}\n"
                f"Difference: computed - reference\n"
                f"range: {vmin:.6e} to {vmax:.6e}"
            ),
            font_size=10,
            color="black",
        )
        _set_camera(plotter, points)
        plotter.show(screenshot=str(output_file))
    finally:
        plotter.close()


def visualize_diff_file(
    diff_file: Path,
    diff_results_dir: Path,
    output_dir: Path,
) -> list[Path]:
    """Render every numeric point field in one diff VTK file."""
    dataset = pv.read(diff_file)
    if not isinstance(dataset, pv.DataSet):
        raise TypeError(f"Unsupported VTK dataset in {diff_file}")

    relative = diff_file.relative_to(diff_results_dir)
    safe_stem = (
        re.sub(r"[^\w.-]+", "_", relative.stem, flags=re.UNICODE).strip("_.")
        or "unnamed"
    )
    file_output_dir = output_dir / relative.parent / safe_stem
    generated: list[Path] = []
    for field_name, points, values in _fields(dataset):
        safe_field = (
            re.sub(r"[^\w.-]+", "_", field_name, flags=re.UNICODE).strip("_.")
            or "unnamed"
        )
        output_file = file_output_dir / f"point_{safe_field}.png"
        render_field(diff_file, output_file, field_name, points, values)
        generated.append(output_file)
    if not generated:
        raise ValueError(f"No numeric point fields found in {diff_file}")
    return generated


def visualize_diff_results(
    diff_results_dir: Path,
) -> tuple[list[Path], list[str]]:
    """Render all supported fieldcompare diff files below a directory."""
    diff_results_dir = diff_results_dir.resolve()
    output_dir = diff_results_dir / "visualizations"
    generated: list[Path] = []
    errors: list[str] = []
    diff_files = sorted(
        path
        for path in diff_results_dir.rglob("*")
        if path.is_file()
        and path.suffix.lower() in SUPPORTED_SUFFIXES
        and "diff" in path.name.lower()
    )
    for diff_file in diff_files:
        try:
            generated.extend(
                visualize_diff_file(diff_file, diff_results_dir, output_dir)
            )
        except Exception as error:
            errors.append(f"Could not visualize {diff_file}: {error}")
    return generated, errors


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Render fieldcompare VTK diff fields as PNG images"
    )
    parser.add_argument(
        "diff_results_dir",
        type=Path,
        help="Directory containing archived fieldcompare diff VTK files",
    )
    args = parser.parse_args()

    if not args.diff_results_dir.is_dir():
        parser.error(f"Not a directory: {args.diff_results_dir}")

    generated, errors = visualize_diff_results(args.diff_results_dir)
    for output_file in generated:
        print(f"Wrote {output_file}")
    for error in errors:
        print(f"WARNING: {error}", file=sys.stderr)

    if generated:
        print(f"Wrote {len(generated)} diff visualization(s)")
    elif not errors:
        print("No fieldcompare diff VTK files found")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
