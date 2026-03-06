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
