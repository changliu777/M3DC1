# AGENTS.md

## Project Scope
- This repository ports IDL `.pro` routines to Python under the `m3dc1` module.
- Keep Python file/function structure close to IDL call structure (1-to-1 where practical).

## User Preferences (Persistent)
- Use Homebrew Python 3.14 for commands:
  - `/opt/homebrew/bin/python3.14`
- Implement module-style code only:
  - No `main` function / no script-only entry requirement.
  - Must support `import m3dc1.<module>` usage.
- When IDL calls a subroutine/function, create a separate Python file for it and call it from other modules.
- Reuse existing Python implementations if already present (do not reimplement duplicates), e.g. `get_normalizations`, `convert_units`.
- For missing field names in HDF5 reads:
  - `read_field.py` should raise `KeyError`.
  - `plot_field.py` should raise `KeyError` when underlying field is missing.

## Plot/Return Behavior
- `plot_mesh.py` should return `(figure, axis)`.
- `plot_field.py` should return `(figure, axis)`.
- `plot_flux_average.py` should return `(figure, axis)`.
- Support `iso` option in plotting functions via `ax.set_aspect('equal')`.
- Default `contourf` levels should be `100`.
- In `contour_and_legend.py`, use:
  - `plasma` colormap if data is all positive or all negative.
  - `turbo` colormap if data crosses zero.

## String/Label Conventions
- Remove IDL text formatting escapes like `!X`, `!6` in Python label strings.
- Convert IDL-style math labels to LaTeX-like strings directly in-place (no extra wrapper function), using `$...$` where needed.
- Preserve intended subscript semantics (e.g. `I_{tot}`), not flattened alternatives.
- Apply similar label conversion for `read_scalar.py`, `read_field.py`, and unit parsing labels.

## Global Migration Rules
- Preserve IDL-style logging/print behavior across all migrated files.
- This rule applies to all future IDL-to-Python migration work in this repository unless explicitly overridden by the user.

## read_field / eval_field Notes
- `eval_field.py` is currently Numba-based (`njit` + `prange`) for element-loop acceleration.
- Worker control is environment-driven:
  - `M3DC1_EVAL_WORKERS` controls requested threads.
  - Clamp to valid Numba thread range internally.

## Performance Testing Conventions
- Prefer direct benchmark comparisons with fixed `points`, reporting:
  - runtime for each worker setting
  - speedup ratio
  - max absolute difference between outputs

## Geometry Mapping Notes (Persistent)
- `read_field.py` now implements `igeometry=1` mapping behavior:
  - For primitive reads, when `igeometry > 0` and `logical=False`, it reads `rst`/`zst`, builds map indices with `create_map`, and remaps data with `map_field`.
  - `logical=True` must bypass geometry remapping.
- IDL subroutine split rule is required here:
  - Keep mapping logic in separate modules: `create_map.py` and `map_field.py`, called from `read_field.py`.
- `create_map.py` implementation requirement:
  - Keep the Numba-decorated implementation.
  - It must still run when Numba is unavailable by using a no-op `njit` fallback (same function body executes in Python mode).
- `map_field.py` implementation requirement:
  - Use `scipy.interpolate.RegularGridInterpolator` (not custom `_interp2`).
  - Use vectorized interpolation calls (e.g., interpolate all `(mx,my)` points at once), not per-point Python-loop interpolation calls.

## Session Decisions (Persistent)
- `eval_field.py`:
  - Keep Numba thread chunksize fixed at `16` via `set_parallel_chunksize(16)`.
  - Keep `debug_spin` support active in `eval_field` and kernel path.
  - In the `debug_spin` loop, call `_eval_poly_numba(...)` (current debug workload behavior).
- `read_field.py` / `plot_field.py` interface updates:
  - Use `timeslices` argument name (not `slices`/`time`).
  - `read_field.py` no longer uses `x,y,t` positional arguments.
  - `plot_field.py` no longer uses `x,y` arguments and uses `timeslices`.
- Multi-timeslice behavior:
  - `read_field.py` supports multiple `timeslices` including `diff` behavior.
  - `plot_field.py` should pass multi-`timeslices` through to `read_field.py`.
- Numba runtime/testing conventions used in this repo:
  - For Numba tests in this environment, use:
    - `source /etc/profile`
    - `module load py-numba`
  - Common benchmark settings used in this session include larger `points` (e.g., `1000`/`10000`) and worker comparisons (`workers=1` vs `workers=4`).
- Plot mesh remap details:
  - In `plot_mesh.py`, when handling `rst_meta`/`zst_meta`, fill masked points using 8-neighbor averaging.
  - Fill target condition is `mask[i, j] == 1` (use `out_mask` from meta).
  - Neighbor inclusion condition is `mask[ii, jj] != 1` (with finite-value check).
- `plot_mesh.py` interface update:
  - Add `phi` argument and pass it to `read_field("rst", ...)` / `read_field("zst", ...)`.
- `plot_field.py` interface update:
  - Add explicit `logical` argument in `plot_field(...)`.
  - Pass `logical` through to `read_field(...)`.
  - Pass `logical` through to all internal `plot_mesh(...)` calls.
- Magnetic probe plotting naming:
  - Python entry is `plot_mag_probes(...)` (plural), corresponding to IDL `plot_mag_probes.pro` in the upper-level IDL directory.
  - It returns `(tdata, data)` and wraps `plot_signals("mag_probes", ...)`.
