# FlutterWing figure style

Every figure uses one style. The tokens (colours, line widths, font and figure sizes) live in
[`util/fw_style.m`](../util/fw_style.m); `help fw_style` lists them. A new figure goes through
[`fw_figure`](../util/fw_figure.m) and [`fw_export`](../util/fw_export.m):

```matlab
S = fw_style();
fig = fw_figure(S.size.standard(1), S.size.standard(2), 'Name', 'my figure');
plot(V_inf, damping);  xlabel('$V_\infty$ (m/s)');  ylabel('Damping (\%)');
fw_export(fig, 'my_figure.png');
```

`Vg_plot`, `pole_plot`, `wing_scene`, `mimo_nyquist`, `plot_wing_layout` and `animate_wing`
follow these rules already.

## Rules

1. A colour appears only with its meaning: ink `#1D1D1F` for text and axes, blue `#0066CC` for
   data, muted `#6E6E73` for secondary text and a second series, grey `#D2D2D7` for control
   surfaces at rest, coral `#E35336` only for instability (a flutter crossing, an unstable pole).
2. All text goes through the LaTeX interpreter, in three sizes: 8 pt small, 9 pt body, 14 pt for
   animation labels. No bold.
3. Three figure sizes: standard 14 x 9.8 cm, wide 14 x 6 cm, square 9.8 x 9.8 cm.
4. Light only: an opaque white background, no transparency, no light/dark pairs. An opaque asset
   stays readable on GitHub's dark theme.

## Generators

| Files | Generator |
|---|---|
| `docs/figures/tutorial/*.png`, `goland_flutter_mode.gif`, `numbers.txt` | `docs/make_tutorial_figures.m` |
| `docs/figures/tutorial/generalized_plant.svg` | `docs/figures/make_generalized_plant.py` |
| `docs/figures/DLM_FEM_Coupling.svg`, `.png` | `docs/figures/make_coupling_diagram.py` |
| `docs/figures/quickstart_vg.png`, `afs_hero.gif` | `docs/figures/make_readme_figures.m` |
| `docs/figures/virtual_flight.gif` | `docs/figures/make_virtual_flight_gif.m` |
| `docs/figures/wordmark.svg`, `mark.svg`, `social_preview.svg`, `.png` | `docs/figures/make_wordmark.py` (PNG command in its header) |
| `docs/figures/wing_layout.svg` | `docs/figures/src/wing_layout.tex` (recipe in its header), then `src/relabel_wing_layout.py` |
| `docs/figures/pipeline.svg` | written by hand |
| `live/QUICKSTART_live.m`, `live/TUTORIAL_live.m` | `docs/make_live_scripts.m` |

Run the MATLAB generators from the repo root after `startup.m`. Run the Python generators with
`uv run --with fonttools --with uharfbuzz --with resvg-py docs/figures/<script>.py`.
