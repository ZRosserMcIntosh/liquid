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

**Splashes and pooling.** Bulk liquid that ends up alone in the air leaves
the bulk. Faster than 2 m/s relative to the liquid's mean it tears into drops,
log-normal about 2.2 mm (estimated: the ligament drops of millimetre-scale
splashes, 1-4 mm); slower, it flies as one blob that rejoins the bulk where it
lands next to it, or stays as a stain anywhere else. Every surface keeps a film
0.3 mm thick in each square centimetre where drops land (estimated); what lands
beyond that pools back into bulk particles there, so a spill makes a puddle, a
heavy spray on a wall makes a run, and light spatter stays as stains. Drops
that land in the bulk liquid, or on a glass body (they run off), join it a
particle at a time. Every transfer is exact in volume.

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

**Shot loads.** A shell fires its pellets (00 buck: nine 3.5 g pellets at
395 m/s) as separate rounds, each deflected from the aim by a Gaussian of rms
angle 6.5 mrad (a cylinder-bore pattern, about 1 in per yard) and launched up
to a shot string of 0.5 m behind the first; each meets panes and the dummy on
its own.

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

**More than life.** `spatter_boost` (1 to 5, the World's `spatterBoost`)
turns it up for films and games: the spatter's volume scales by it, the drop
sizes by its fourth root and the speeds by its 0.15 power, the bleeding by its
square root. The lab's "a little Tarantino" is 3: three times the liquid,
drops a third larger. 1 (the default) is as measured.

**Settling.** A drop that has slowed to its terminal speed (within 3%, its
sideways speed under 1 cm/s plus 3% of it) falls straight down from there in
still air; it is landed at once where it would land: straight down, after
coasting on sideways over its drag relaxation time vt/g, onto the first
static surface below, found by ray. Most fine spatter spends a second or two
at under 1 m/s in the air: this lands it within a fraction of a millimetre of
where flying would (tested), and the air empties three times faster. A drop
over something that may move (the ragdoll, a falling piece of glass) is
flown on. `Spray::settle` turns it off.

**The drop budget.** A wound may use a third of the drops the budget has left
and always at least 48, so every pellet of a shot load and every later shot
throws its spatter. Spatter with more drops than that is sampled: each drop
flown stands for several identical ones (its `weight`, the super-droplet
method of Shima et al. 2009), carrying their volume, so volumes stay exact
and the pattern keeps its size and direction statistics, only sparser. A
full budget is overrun by those few dozen drops a wound, never more.

**Panes.** A round reaching a glass pane is handed to the host (`deferPanes`,
the browser lab hands it to Glass Dynamics), to a callback (`paneResolver`),
or resolved by the built-in estimate: 18 J per mm of glass absorbed, a hole of
1.5 calibres. The answer (penetrated, speed out, hole radius, the pane gone)
comes back through `resolveStrike`; drops pass through holes.

**Glass in pieces.** A host running Glass Dynamics also hands over the
pieces of a broken pane (see docs/integration.md): Glass Dynamics' map of
which piece each cell of the pane holds, each piece's outline and thickness,
and every piece's pose as it moves. A monolithic pane is then open wherever
its piece has gone (a crazed tempered pane that falls out lets the spatter
through behind it); a laminate is held by its interlayer and stays closed but
for the holes. Loose pieces are slab colliders (the outline, with its holes,
extruded to the thickness) that fly on ballistically between the host's
updates, found through a 5 cm hash grid. Drops land on them with their
velocity relative to the glass, and a stain on any piece, held or loose,
rides with it; a piece that breaks again hands its stains to the new piece
under them, and one that crumbles into tempered dice takes them with it. The
stain analysis uses only stains still where they landed. The pieces fall and
break by Glass Dynamics' physics alone: the liquid on them does not weigh
them down.


## Ragdoll

The ragdoll is optional (the scene key `ragdoll`, off by default): without it
the dummy stands where it is and the rounds pass through it.

**Joints.** The thorax (both chests) is one rigid piece; the neck's base is
fixed to it and its top and the head nod and turn within about ±25°; the
pelvis flexes against the thorax at the lumbar spine within about ±15°.
Shoulders and hips swing within elliptical cones in the trunk's frame
(shoulder flexion −50° to 170°, abduction −30° to 85°; hip flexion −20° to
120°, abduction −20° to 45°: the usual adult ranges, AAOS and Norkin &
White, a little narrowed). Knees and elbows are hinges about the left-right
axis of what they hang from: they bend one way, not past straight, fold no
further than about 120°, with some play out of their plane (the segments
cannot twist, so a hip's or a shoulder's rotation shows as play there).
Ankles rock about ±25°. Limbs collide with the trunk and each other. Hard
hinges and folds fight the floor and the other limits in a position-based
solver and pump energy in, so they are soft (a third of the way back each
pass), and each 2 ms substep may not leave the body with more mechanical
energy than it began with plus what a hit gave it: a limp body with friction
only loses energy, and that keeps it so.

**Hits.** A round gives the struck capsule the momentum it loses, at the point
it strikes, shared between the capsule's ends as a lever, so a hit off the
centre turns the body. As in life that barely moves it: a 9 mm round stopping
in 78 kg moves it about 3 cm/s, and people shot do not fly backwards.
`knockback` (0 to 0.3) is the cinema's version: at least that fraction of the
energy the round deposits becomes the body's motion (a load of shot shares it
out over its pellets), J = √(2 M k E). The lab's tarantinoness slider sets it
from 0 to 0.1 (at 0.05, a 9 mm hit puts 12 J into the body, about 0.6 m/s; a
load of buckshot about 1.6 m/s).

The dummy's capsules are two point masses each (the dummy's 78.2 kg, the
Hybrid III 50th-percentile mass, shared by capsule volume), held at their
lengths; capsules that meet are joined (ball joints where the ends coincide,
pinned shoulders and hips where they are apart), the chest, pelvis, neck and
head are one rigid cluster, and knees and elbows cannot fold flat. It is
integrated by position-based dynamics at 480 Hz with 10 constraint
iterations, the ends and middle of every capsule colliding with the static
scene with Coulomb friction, and it settles by static friction below 3 cm/s.
Until a round hits it, it stands on a release stand (a stunt rig); a hit
releases it and gives the struck capsule the momentum the round lost, shared
between its ends as a lever. Wounds and stains on it are attached to their
capsule by position along the axis and angle round it from the trunk's facing,
and move with it; drips leave with the wound's velocity. There is no
self-collision between limbs and no joint motor or muscle tone: it falls as a
limp body.

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
before it exists and after it breaks; the pieces of a broken pane the same
way (each its outline with its holes, ear-clipped and extruded), with the
stains on them as their children; the stains as flat ellipses that appear in
the frame they land; and the static scene. Accessor bounds are written
exactly (9 significant digits round-trip a float32). `exportPly` writes one
frame as a binary PLY point cloud (position, pscale, kind: bulk or drop);
`exportStainsJson` lists the stains for placing decals.
