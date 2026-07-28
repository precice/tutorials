# System Tests Coverage Matrix (Baseline)

Concrete tracking artifact for [#448](https://github.com/precice/tutorials/issues/448).
This file complements `tools/tests/tests.yaml` and is intentionally easy to maintain.

- Last sync: 2026-07-28
- Generated from: `tests.yaml`, tutorial `metadata.yaml`, issue `#448`

## Status legend

- `Included`: explicitly covered in `tests.yaml` (typically in `release_test`).
- `Partial`: test exists, but regression/reference coverage is incomplete.
- `Missing`: tutorial cases exist, but are not in release coverage yet.
- `Blocked`: known external blocker (license, unmaintained dependency, etc.).
- `Needs verification`: listed in issue discussion, but no direct suite evidence recorded here yet.

## Component/feature coverage

| Component / feature | Status | Source of truth | Current evidence | Next step / note |
| --- | --- | --- | --- | --- |
| bare (preCICE core only) | Included | `tests.yaml` | `quickstart` (`fluid-openfoam` + `solid-cpp`), `elastic-tube-1d` (`fluid-cpp` + `solid-cpp`) in `release_test` | Keep as baseline |
| Python bindings | Included | `tests.yaml` | `elastic-tube-1d` (`fluid-python`, `solid-python`) in `release_test` | Keep as baseline |
| OpenFOAM adapter | Included | `tests.yaml` | `quickstart`, `flow-over-heated-plate`, `perpendicular-flap`, `multiple-perpendicular-flaps` | Track partial external/monolithic variants separately |
| CalculiX adapter | Included | `tests.yaml` | `perpendicular-flap`, `heat-exchanger-simplified` | Keep as baseline |
| deal.II adapter | Included | `tests.yaml` | `perpendicular-flap`, `multiple-perpendicular-flaps` | Keep as baseline |
| FEniCS adapter | Included | `tests.yaml` + issue `#448` notes | `flow-over-heated-plate`, `perpendicular-flap` in `fenics_test`/`release_test` | `elastic-tube-3d` OpenFOAM/FEniCS combo still noted as problematic |
| DUNE / DUNE-FEM / DuMux adapters | Included | `tests.yaml` | `free-flow-over-porous-media`, `two-scale-heat-conduction` suites | Keep sensitivity/tolerance notes reviewed |
| Micro-Manager | Included | `tests.yaml` | `two-scale-heat-conduction` in `micro_manager_test` | Keep as baseline |
| SU2 adapter | Included | `tests.yaml` | `perpendicular-flap` (`fluid-su2` + `solid-fenics`) | Keep as baseline |
| FEniCSx adapter | Needs verification | issue `#448` | Mentioned as included in issue list | Add explicit suite/case evidence in next update |
| FMI runner | Needs verification | issue `#448` | Mentioned as included in issue list | Add explicit suite/case evidence in next update |
| MercuryDPM | Needs verification | issue `#448` | Mentioned as included in issue list | Add explicit suite/case evidence in next update |
| ASTE | Partial | issue `#448` + [#877](https://github.com/precice/tutorials/issues/877) | `aste-turbine` noted as tested without full regression checks | Add/maintain robust regression references |
| Julia bindings | Missing | issue `#448` | `resonant-circuit/*-julia` cases listed in issue | Land version override support ([#886](https://github.com/precice/tutorials/pull/886)) |
| Rust bindings | Missing | tutorial metadata + issue `#448` | `elastic-tube-1d` includes `fluid-rust`, `solid-rust` in metadata | Land version override support ([#885](https://github.com/precice/tutorials/pull/885)) |
| G+Smo adapter (external) | Missing | issue `#448` | `partitioned-heat-conduction*`, `perpendicular-flap-stress` listed in issue | Add executable-path based integration workflow |
| MATLAB bindings | Blocked | issue `#448` | listed in issue | Needs licensed runner |
| code_aster adapter | Blocked | issue `#448` + [code_aster-adapter#29](https://github.com/precice/code_aster-adapter/issues/29) | listed in issue | Unmaintained adapter |

## Regression coverage gaps explicitly called out in #448

| Case group | Status | Source of truth | Note |
| --- | --- | --- | --- |
| `aste-turbine` | Partial | issue `#448` | Tested, but full regression checks still incomplete |
| external heat-conduction OpenFOAM variants | Partial | issue `#448` | Listed as external cases without reference results |
| monolithic variants (`partitioned-burgers-1d`, `partitioned-pipe-two-phase`) | Partial | issue `#448` | Listed as tested without full regression checks |

## Quick update workflow

Run these before updating this file:

```bash
python tools/tests/print_test_suites.py
python tools/tests/print_metadata.py
```

Then update only:

1. `Status`
2. `Current evidence`
3. `Next step / note`
4. `Last sync` date
