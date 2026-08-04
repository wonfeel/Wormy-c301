# Wormy

![C++](https://img.shields.io/badge/C%2B%2B-17-blue?style=flat)
![Platform](https://img.shields.io/badge/platform-Windows-0078d7?style=flat)

A C. elegans body + nervous system simulation. The real 401-neuron connectome
([Cook et al. 2019](https://doi.org/10.1038/s41586-019-1352-7)) drives a
physically simulated body through resistive-force-theory drag on an
anisotropic substrate — no scripted locomotion, movement is a straight
physical consequence of neural activity acting on body curvature. Built on
[Tessera](https://github.com/wonfeel/Tessera).

![demo](assets/demo.gif)

![screenshot](assets/screenshot.png)

## What it does

- Loads the real C. elegans connectome (401 neurons, chemical + gap-junction
  synapses) into a sparse non-spiking network integrator
- One-directional pipeline: network activity → muscle curvature → segmented
  body → RFT physics → position/heading — nothing shortcuts this chain
- Live ImGui panels for network gains, locomotion parameters, per-neuron
  activity, body-angle profile (press **F10** to hide them for a clean view,
  **F9** to record a clip)
- `tests/worm_*_calibration` — calibration/validation harnesses, findings
  (including negative ones) documented alongside the code

## Build

```bash
cmake -S . -B out/build -G Ninja
cmake --build out/build --target Demo_worm
```

Pulls the [Tessera](https://github.com/wonfeel/Tessera) engine automatically
via CMake `FetchContent` — nothing to clone by hand. Windows + MSVC, same
requirements as Tessera itself.

## Docs

- [WORM.md](demo/worm/WORM.md) — full reference: network model, body physics,
  every `WormSim::Params` field
- [NEURONS.md](demo/worm/NEURONS.md) — neuron naming/nomenclature
- `demo/worm/WORM_V*_RESULTS.md` — calibration history
