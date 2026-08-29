#!/usr/bin/env python3
"""Render fieldcompare VTK diff fields as PNG images."""

from __future__ import annotations

import argparse
import os
import re
import sys
from collections.abc import Iterator
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import get_context
from pathlib import Path

import numpy as np
import pyvista as pv


SUPPORTED_SUFFIXES = {".vtk", ".vtp", ".vtu"}
WINDOW_SIZE = (1024, 768)


def _log(message: str) -> None:
    """Print immediately so CI and redirected local runs show progress."""
    print(message, flush=True)


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


def _point_size(n_points: int) -> float:
    """Screen-space point size that stays visible without filling the window."""
    return float(max(4.0, min(18.0, 350.0 / max(np.sqrt(n_points), 1.0))))


def _has_nonzero_diff(values: np.ndarray) -> bool:
    """Return whether the field contains any non-zero difference."""
    return bool(np.any(np.abs(values) > 0))


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
) -> bool:
    """Render one field as colored point sprites. Returns True if a PNG was written."""
    if not _has_nonzero_diff(values):
        return False

    point_cloud = pv.PolyData(points)
    scalar_name = "difference"
    point_cloud.point_data[scalar_name] = values

    output_file.parent.mkdir(parents=True, exist_ok=True)
    vmin = float(np.min(values))
    vmax = float(np.max(values))
    plotter = pv.Plotter(off_screen=True, window_size=WINDOW_SIZE)
    try:
        plotter.set_background("white")
        plotter.add_mesh(
            point_cloud,
            scalars=scalar_name,
            cmap="coolwarm",
            clim=(vmin, vmax),
            scalar_bar_args={"title": field_name},
            render_points_as_spheres=True,
            point_size=_point_size(len(points)),
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
    return True


def visualize_diff_file(
    diff_file: Path,
    diff_results_dir: Path,
    output_dir: Path,
) -> list[Path]:
    """Render every numeric point field with a non-zero diff in one VTK file."""
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
    saw_fields = False
    for field_name, points, values in _fields(dataset):
        saw_fields = True
        safe_field = (
            re.sub(r"[^\w.-]+", "_", field_name, flags=re.UNICODE).strip("_.")
            or "unnamed"
        )
        output_file = file_output_dir / f"point_{safe_field}.png"
        if render_field(diff_file, output_file, field_name, points, values):
            generated.append(output_file)
    if not saw_fields:
        raise ValueError(f"No numeric point fields found in {diff_file}")
    return generated


def _worker_count() -> int:
    cpu_count = os.cpu_count() or 1
    return max(1, min(cpu_count, 8))


def _visualize_one_file(
    diff_file: str,
    diff_results_dir: str,
    output_dir: str,
    index: int,
    total: int,
) -> tuple[list[str], str | None]:
    """Render one file in a worker process. Returns (output paths, error)."""
    path = Path(diff_file)
    _log(f"[{index}/{total}] Rendering {path.name}")
    try:
        generated = visualize_diff_file(
            path, Path(diff_results_dir), Path(output_dir)
        )
        if generated:
            names = ", ".join(p.name for p in generated)
            _log(
                f"[{index}/{total}] Wrote {len(generated)} image(s) from "
                f"{path.name}: {names}"
            )
        else:
            _log(
                f"[{index}/{total}] Skipped {path.name} "
                "(only zero-difference fields)"
            )
        return [str(p) for p in generated], None
    except Exception as error:
        message = f"Could not visualize {path}: {error}"
        _log(f"[{index}/{total}] WARNING: {message}")
        return [], message


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
    total = len(diff_files)
    workers = min(_worker_count(), total) if total else 1
    _log(f"Found {total} fieldcompare diff VTK file(s), using {workers} worker(s)")
    if not diff_files:
        return generated, errors

    if workers == 1:
        for index, diff_file in enumerate(diff_files, start=1):
            paths, error = _visualize_one_file(
                str(diff_file),
                str(diff_results_dir),
                str(output_dir),
                index,
                total,
            )
            generated.extend(Path(p) for p in paths)
            if error:
                errors.append(error)
        return generated, errors

    with ProcessPoolExecutor(
        max_workers=workers,
        mp_context=get_context("spawn"),
    ) as executor:
        futures = {
            executor.submit(
                _visualize_one_file,
                str(diff_file),
                str(diff_results_dir),
                str(output_dir),
                index,
                total,
            ): diff_file
            for index, diff_file in enumerate(diff_files, start=1)
        }
        for future in as_completed(futures):
            diff_file = futures[future]
            try:
                paths, error = future.result()
            except Exception as error:
                errors.append(f"Could not visualize {diff_file}: {error}")
                continue
            generated.extend(Path(p) for p in paths)
            if error:
                errors.append(error)
    generated.sort()
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
    for error in errors:
        print(f"WARNING: {error}", file=sys.stderr, flush=True)

    if generated:
        _log(f"Wrote {len(generated)} diff visualization(s)")
    elif not errors:
        _log("No fieldcompare diff VTK files found")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
