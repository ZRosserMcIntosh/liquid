# Physics

How each part of the core is modelled, how far it goes, and its main gap.
Units are SI; the world frame is glTF's (metres, right-handed, +y up).

## The step

`World::advance(t)` runs whole internal steps until the clock reaches `t`.
The step is chosen from the state alone: the liquid's CFL limit (0.4 of the
particle spacing per step at the fastest particle) and its explicit viscous
limit, 2 ms at most, 50 µs while a round is in flight and 0.2 ms while fresh
spatter is airborne. A step does, in order: rounds (events at their exact
times inside the step), glass bodies and their strikes, the bulk liquid,
the drops (each on its own clock, sub-stepped to 4 mm of travel), the
exchange between bulk and drops, and the wounds' bleeding.

## Bulk liquid

**Solver.** Position-based fluids (Macklin & Müller 2013): the density
constraint C_i = ρ_i/ρ0 − 1 ≥ 0 is solved by 3-5 Jacobi iterations per
step on the cubic B-spline kernel with support 2 × spacing. A correction is
bounded per iteration in length (0.2 spacing) and as a speed (5 m/s), so a
particle squeezed between two moving shards comes out over a few steps
rather than taking a kick that would grow as the step shrinks. The particle
mass carries a relaxed-packing factor of 1.05: liquid that has relaxed under
the constraint packs about 5% more particles into a volume than the lattice it
was poured as, at the same kernel sum (calibrated so a resting column keeps
its poured height, `tests/test_scenes.cpp`).

**Walls.** No wall particles: a particle at signed distance d from a solid
sees the fraction Φ(d) of its kernel inside it (the half-space integral of the
kernel, tabulated once), which adds ρ0 Φ to its density and pushes along the
normal with the kernel's cross-section A(d). Several walls combine by
inclusion-exclusion, Φ = 1 − Π(1 − Φ_k), exact for a separable kernel and close
for the spline, so corners do not count the overlap twice. Φ is shifted so a
lattice's first layer rests half a spacing from the wall, as it would against
a mirrored lattice. A hard projection keeps every particle 0.4 spacing out.
Walls are free-slip (a water boundary layer is far thinner than a particle);
contact-angle hysteresis gives the slip a pinning force σ (cos θ_r − cos θ_a)
per unit length of contact line, per surface, so a spill stops spreading; a
particle nearly alone on a surface is a landed drop and stays.

**Viscosity.** The SPH Laplacian of Morris et al. (1997) with the liquid's
dynamic viscosity. For blood the viscosity follows the Carreau-Yasuda fit of
Cho & Kensey (1991) on the local strain rate (from the SPH velocity gradient),
56 mPa s at rest falling to 3.45 at high shear. A small XSPH smoothing
(0.01) damps particle noise; it is numerical, and at the coarse tiers it is
larger than wine's own viscosity.

**Surface tension.** The continuum surface force σ κ n (Brackbill, Kothe &
Zemach 1992; in SPH after Morris 2000) with the liquid's real σ: the colour
gradient marks the surface, its curvature is the divergence of the unit
normal (Shepard-normalised over the neighbours that carry one). A flat surface
feels nothing. Akinci et al.'s cohesion model was tried first and rejected:
its coefficient has no physical scale and at the particle sizes used it
pulled even a flat surface inward at about 2 g, pumping energy into a
resting tank.

**Coupling.** The constraint's wall pushes and the hard projection's velocity
changes are momentum taken from the wall's body: they go back to it as an
impulse at the contact point, applied at the next step, relaxed by
m / (m + m_liquid in contact) (an implicit added-mass correction for the
explicit coupling). A body only a few particle masses heavy cannot be
resolved by the liquid round it: it takes back the share (n − 20)/180 of the
push for n particle masses, from one-way below 20 to two-way above 200.

**Main gap.** Below the particle size (3-8 mm) the liquid is not resolved:
films thinner than a particle, the crown of a splash and the neck of a pour
are drops, not liquid. Wine's equilibrium puddle (2-3 mm) is thinner than a
particle, so a spill covers about the right area but stands a few millimetres
proud.

## Drops

A drop is a sphere in still air (ρ 1.204 kg/m³, μ 1.813e-5 Pa s, ISO 2533 at
20 °C) with gravity and the drag of Schiller & Naumann below Re 1000 and
Newton's 0.44 above, integrated semi-implicitly (v' = (v + g dt)/(1 + k dt)),
which is stable for the finest mist at any step. Above an air Weber number of
12 a drop breaks up after Pilch & Erdman's (1987) time 6 t* (We − 12)^−¼,
t* = (d/v) sqrt(ρ_l/ρ_a), into drops at a Weber number of 6 (at most 32). Drops
landing on the bulk liquid join it; blobs of the bulk liquid that went
airborne rejoin it where they land, or, if they splash, are thrown off as
spatter with the probability of the splash's satellite fraction (volume
exact). Drops do not collide with each other, and deformation does not add
drag.

## Stains

On the normal component of the impact velocity (Re and We of the drop): the
maximum spreading diameter of Laan et al. (2014), D/d = Re^⅕ P^½ / (1.24 + P^½),
P = We Re^−⅖, fitted on blood and water, at least the drop's own size; for a
shear-thinning liquid at its infinite-shear viscosity. Length W / sin α, the
relation bloodstain analysis inverts; tail along the tangential velocity.
Splashing above K = We^½ Re^¼ > 57.7 (Mundo, Sommerfeld & Tropea 1995),
times the surface's splash factor (rough concrete splashes earlier, smooth
glass later): spines and 2-14 satellite drops carrying 4-25% of the volume
onward along the tail. Absorbent surfaces wick the stain wider and keep the
satellites. Stains do not run down walls (no flow patterns), and a stain's
outline is an ellipse plus a spine count, not a traced shape.

## Glassware and strikes

**Geometry.** A vessel is the mid-surface profile of its wall, (r, y,
thickness) from the rim down the inside to the base, revolved about +y. The
signed distance to the glass is the distance in the meridian plane to the
profile (and its mirror across the axis) minus half the thickness; a solid
stem is a rod written as r = t/2, which Pappus' theorem then gives exactly.
The cavity is the polygon of the inner surface closed along the axis; its
volume below any level is 2π ∫∫ r dA by the polygon formula. Mass and inertia
come from Pappus on the wall's sections.

**Contact.** Rigid bodies (the vessel, later its shards) collide with the
static scene through points on their glass surface: sequential impulses with
speculative contacts, Coulomb friction, restitution from the surface, and
positional correction that moves the body out without adding velocity. Shards
do not collide with each other.

**Breaking.** At every strike faster than 0.25 m/s, at the struck point: the
effective mass m_eff = 1 / (1/m + (r×n)·I⁻¹(r×n)); the peak force from the
energy ½ m_eff v² stored in two springs in series, Hertz contact on the local
radius of the glass (√(R1 R2) of the wall; half the thickness on a rim or
edge) against the surface, and the wall's own point-load stiffness, Reissner's
4 E t² / (R √(3(1−ν²))); a compliant floor (carpet pile on underlay) is a
Winkler foundation of stiffness E/layer that spreads the load over its own
thickness and bottoms out onto concrete. The stress: bending from the part of
the force normal to the wall, 3(1+ν)F/(2πt²)(1 + ln(reach/a)) (a load over a
small circle on a plate, the logarithmic term cut off at the shell's bending
reach 0.76 √(R t)), ×1.3 on a free edge, plus membrane stress from the part
along the wall. The wall breaks with Weibull probability
1 − exp(−(σ/σ0)^m), drawn from the seeded generator. σ0 = 300 MPa for the
stemware (m = 6) is calibrated to everyday outcomes (a 0.3 m/s knock never
breaks it, tipping it over on a table sometimes does, a fall onto tile always
does and onto carpet rarely); it is plausible for fire-polished glass over
the few square millimetres a strike bends, far above large-area design values.

**Shards.** 4-64 Voronoi seeds in the developed wall (arc length along the
profile, azimuth times radius), 75% crowded round the strike within an
exponential spread, the rest anywhere by area; each wall cell of the mesh
belongs to the nearest seed. Each shard is a rigid body with the mass and
inertia of its cells, the vessel's velocity at its centre plus a small burst
away from the strike; the liquid sees each shard as the vessel wall
restricted to its cells, so it escapes through the gaps. Shards do not break
again.

**Main gap.** Strikes are judged one at a time; repeated sub-critical
strikes do not weaken the glass, and the crack pattern is a Voronoi
partition, not grown cracks (Glass Dynamics grows cracks in flat panes; it
does not do curved vessels).

## Rounds, the dummy and spatter

**Penetration.** A round crosses the union of the dummy's capsules along its
path with Poncelet's law, m dv/dx = −(R A + ½ ρ C_d A v²): R = 4.9 MPa is
fitted so a 9 mm 124 gr FMJ at 360 m/s penetrates 0.65 m of 10% gelatin (the
FBI-protocol range is 0.6-0.75 m); ρ = 1060 kg/m³ and C_d 0.25-0.4 per round.
The time and path inside follow exactly from v = √(a/b) tan(θ0 − ω t). A round
slower than 25 m/s inside the dummy lodges. Rounds do not yaw, deform or
fragment.

**Spatter.** At the entry wound, back spatter; at the exit, forward spatter:
a volume proportional to the energy the round deposited, in drops of a
log-normal size distribution, at speeds up to a fraction of the round's,
in a cone round the reversed path (back) or the path (forward), tilted
towards the surface normal. These are orders of magnitude from the published
high-speed films of back and forward spatter (Comiskey, Yarin, Kim &
Attinger 2016-2019; Kabaliuk et al. 2014), not fits: forward spatter carries
more liquid, in larger, faster drops. Each wound then drips (0.8 ml/s,
decaying over 8 s): drops of Tate's volume 0.6 · 2π r σ / (ρ g) from a 1.5 mm
wound edge, about 30 µl.

**Panes.** A round reaching a glass pane is handed to the host (`deferPanes`,
the browser lab hands it to Glass Dynamics), to a callback (`paneResolver`),
or resolved by the built-in estimate: 18 J per mm of glass absorbed, a hole of
1.5 calibres. The answer (penetrated, speed out, hole radius, the pane gone)
comes back through `resolveStrike`; drops pass through holes, and a pane that
has fallen is disabled.

## Stain analysis

For the stains of one surface face (the majority normal): stains under 1 mm
long, rounder than W/L 0.94 or grazing below 10° are left out, as an analyst
would; satellites are left out unless asked for. Widths and lengths can be
read to a ruler's resolution (0.25 mm in the lab). α = asin(W/L); the
direction is the tail's. The area of convergence is the point of the surface
plane closest in least squares to the lines back through the stains; each
stain's height is its distance to it times tan α; the area of origin is the
convergence raised by the median height (with the interquartile range). The
truth is the source of most of the stains used, and the report gives the
convergence's miss in the plane and the height's error. The straight-line
method ignores gravity and drag, so origins far from a horizontal surface come
out too high; this is reported, not corrected.

## Export

`exportGlb` writes one binary glTF: per recorded frame a mesh of the liquid's
surface (marching tetrahedra over a smooth splat of the particles, shared
vertices, area-weighted normals) and the drops as small octahedra, shown for
its frame alone by a STEP scale animation; every glass body and shard as a
mesh in its own frame with LINEAR translation and rotation tracks, hidden
before it exists and after it breaks; the stains as flat ellipses that appear
in the frame they land; and the static scene. Accessor bounds are written
exactly (9 significant digits round-trip a float32). `exportPly` writes one
frame as a binary PLY point cloud (position, pscale, kind: bulk or drop);
`exportStainsJson` lists the stains for placing decals.
