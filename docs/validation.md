# Validation

What each scenario gives against published data, or against the everyday
experience a value was calibrated to. The table is printed by
`./build/lqdemo validate`; regenerate it after changing the physics or the
data. Rows marked *calibration point* are not independent checks: they are
where a constant was set (see the provenance in `data/` and
[physics.md](physics.md)).

<!-- validation:begin -->
| Case | Reference | Model | Source |
|---|---|---|---|
| Terminal speed, 0.5 mm water drop | 2.06 m/s | 2.03 m/s | Gunn & Kinzer (1949) |
| Terminal speed, 1.0 mm water drop | 4.03 m/s | 3.86 m/s | Gunn & Kinzer (1949) |
| Terminal speed, 2.0 mm water drop | 6.49 m/s | 6.93 m/s | Gunn & Kinzer (1949) |
| Blood viscosity at 1 / 100 / 1000 1/s | Carreau fit | 27.1 / 4.71 / 3.74 mPa s | Cho & Kensey (1991) |
| Dam-break surge front Z at T = 1.8 / 2.5 | 2.0 / 3.0 | 2.19 / 3.00 | Martin & Moyce (1952) |
| Penetration in 10% gelatin, 9 mm FMJ | 0.60-0.75 m | 0.65 m | FBI-protocol tests (calibration point) |
| Penetration in 10% gelatin, .44 Magnum | - | 0.88 m | model |
| Penetration in 10% gelatin, 7.62 x 51 mm ball | - | 2.12 m | model |
| Bordeaux wine glass: mass, brimful capacity | 140-220 g, - | 152 g, 631 ml | machine-made stemware |
| Bordeaux wine bottle: mass, brimful capacity | 400-550 g, 765-775 ml | 472 g, 770 ml | 750 ml bottle norms |
| Wine glass bowl at 0.3 m/s on Oak table top | never breaks | P(break) 0.001 (88 MPa) | calibration point |
| Wine glass bowl at 1.0 m/s on Oak table top | sometimes | P(break) 0.301 (253 MPa) | calibration point |
| Wine glass bowl at 3.8 m/s on Glazed ceramic tile | breaks | P(break) 1.000 (1356 MPa) | calibration point |
| Wine glass bowl at 3.8 m/s on Cut-pile carpet | rarely | P(break) 0.001 (101 MPa) | calibration point |
| Area of origin from a pane 0.7 m behind the exit wound | 0.70 m | 0.58 m (convergence off by 36 mm) | tangent method |
| Area of origin from the floor below a 1.3 m wound | 1.30 m | 3.61 m | tangent method: overestimates (Carter 2001) |
<!-- validation:end -->

## Notes on the rows

- **Terminal speeds.** Gunn & Kinzer (1949) measured falling water drops in
  still air. The model is a rigid sphere on the standard drag curve; real
  drops above about 1.5 mm flatten and fall more slowly, so the 2 mm row runs
  about 7% fast. Spatter drops are mostly under 1 mm, where the agreement is
  within 5%.
- **Dam break.** Martin & Moyce (1952) released a column of water of width a
  and height 2a; the front position Z = x/a is plotted against
  T = t sqrt(2g/a). The scene is the 0.146 m column of Koshizuka & Oka (1996)
  between glass walls a few particles apart. Inviscid SPH usually runs a little
  ahead of the experiment early on, as here.
- **Penetration.** The Poncelet resistance of the simulant is fitted to one
  point (9 mm FMJ, 0.65 m); the other rounds follow from it. Real rifle bullets
  yaw and tumble after 10-20 cm of tissue and stop far sooner than the 2.1 m
  given for a round that stays nose-first: yaw is not modelled.
- **Glass strikes.** The Weibull scale of the stemware is set so that the four
  everyday outcomes hold. No published drop-test series of stemware was found
  to fit to; the rows show what the model does at those points.
- **Stain analysis.** On a pane 0.7 m behind the exit wound the strings
  converge within a few centimetres of the true source and put it a little
  low; on the floor, 1.3 m below the wounds, the straight strings put the
  origin far too high because the drops arrive on falling arcs. Both are the
  known behaviour of the tangent method (Carter 2001; de Bruin et al. 2011
  for the impact-angle relation on blood); the panel in the lab shows the same
  comparison for every run.

## Tests

`./build/tests/lqd_tests` runs 33 cases. The scene cases run at full length in
optimised builds and shortened in Debug and sanitizer builds (same invariants,
fewer simulated seconds). Among them:

- the liquid's volume is conserved to 1e-9 through every scene: poured volume =
  bulk + drops in the air + stains + absorbed + lost;
- the same seed gives the same stains whether the host advances the world at
  60 Hz or 240 Hz;
- the glTF export parses, its chunk lengths agree, and it holds the liquid
  flipbook, the glass bodies and the stains (the files also pass the Khronos
  glTF validator with no errors or warnings).
