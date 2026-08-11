# Montserrat

Vendored because the figures and Word outputs are set in Montserrat and nothing
resolves the family for us. Naming a font in an SVG only asks for it: the
manuscript builder rasterises every figure through cairosvg, and CI installs
`libcairo2` and no fonts at all, so an unresolved family falls back to DejaVu
Sans without failing anything. The fallback is invisible to the release gates
and visible in the shipped PDF, which is the wrong way round.

Static TrueType builds taken from the upstream project at
<https://github.com/JulietaUla/Montserrat>, `fonts/ttf/`.

`Montserrat-Medium.ttf` and `Montserrat-SemiBold.ttf` register under their own
family names (`Montserrat Medium`, `Montserrat SemiBold`) rather than as weights
of `Montserrat`, which is how static exports work. Ask for them by full family
name; asking `Montserrat` for a mid weight silently gets you Regular or Bold.

`tools/brand.py` points fontconfig at this directory, so nothing needs to be
installed system-wide for figure generation. Rendering the Word files through
LibreOffice is the exception, since it is a separate process that reads the
ambient font configuration:

```bash
mkdir -p ~/.local/share/fonts
cp assets/fonts/montserrat/*.ttf ~/.local/share/fonts/
fc-cache -f
```

Licensed under the SIL Open Font License 1.1; see `OFL.txt`. The licence allows
redistribution alongside this repository provided that file travels with the
fonts.
