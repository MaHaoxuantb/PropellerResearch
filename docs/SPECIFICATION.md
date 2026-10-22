# Parametric Toroidal Propeller: Geometry, CFD and Optimization Proposal

**Specification:** Version 1.0  
**Status:** Proposed implementation baseline  
**Primary technologies:** Python, OpenFOAM, NumPy, SciPy, Trimesh, simulated annealing  
**Application:** Automated aerodynamic optimization of a rotationally symmetric, loop-bladed propeller

---

## 1. Project objectives

### 1.1 Primary objective

Develop an automated system to generate, evaluate, and optimize the three-dimensional geometry of a toroidal or loop-bladed propeller.

The system must:

1. Represent complex, three-dimensional blade shapes with relatively few design parameters.
2. Generate physically coherent geometry while minimizing invalid candidate designs.
3. Support adjustable blade count, blade curvature, chord distribution, twist, and airfoil geometry.
4. Export watertight geometry suitable for computational fluid dynamics.
5. Evaluate thrust, torque, shaft power, and aerodynamic performance using OpenFOAM.
6. Optimize shape parameters using simulated annealing.
7. Support increasingly flexible geometry representations without requiring a redesign of the optimization pipeline.
8. Preserve a reproducible record of every evaluated geometry, simulation, and result.

### 1.2 Optimization target

For the initial project, use a fixed operating condition with prescribed rotational speed, fluid properties, and maximum propeller diameter.

Minimize shaft power while producing at least a target thrust:

\[
\boxed{
\begin{aligned}
\min_{\mathbf p}\quad &P(\mathbf p)\\
\text{subject to}\quad
&T(\mathbf p)\ge T_{\text{target}},\\
&2R_{\max}(\mathbf p)\le D_{\max},\\
&\mathbf p\in\mathcal G_{\text{valid}}.
\end{aligned}}
\]

Here:

- \(\mathbf p\) is the vector of geometric design parameters.
- \(P\) is shaft power.
- \(T\) is thrust.
- \(D_{\max}\) is the allowed propeller diameter.
- \(\mathcal G_{\text{valid}}\) is the feasible set of mechanically and geometrically acceptable designs.

Mass, structural stress, and noise may be added as constraints or secondary optimization objectives.

The fixed operating condition is essential. Comparing unrelated thrust, RPM, and diameter combinations can lead to misleading optimization results.

### 1.3 Initial assumptions

The mathematical model assumes:

- The rotor axis is the global \(z\)-axis.
- The propeller consists of \(N_b\) identical, rotationally repeated blade units.
- Each blade unit is a continuous U-shaped or loop-shaped surface, attached near the hub at two ends.
- The hub is initially fixed.
- All repeated blades have identical geometry.
- Blade flexibility and aeroelastic deformation are neglected in the initial CFD model.
- The airfoil section is initially fixed and is optimized only in a later development stage.

If the intended physical design instead comprises separate inner and outer aerodynamic strips connected at their tips, the geometry generator must represent those strips explicitly. The overall optimization architecture remains applicable.

---

## 2. Design philosophy

The central design principle is:

**Introduce only degrees of freedom that can produce meaningful aerodynamic improvements, and constrain them in ways that make geometrically invalid candidates uncommon.**

Unrestricted polynomial models are not preferred.

They can introduce:

- Unwanted oscillations.
- Negative thickness or width.
- Local curvature spikes.
- Self-intersections.
- Poor mesh quality.
- Unnecessarily large optimization search spaces.

The preferred mathematical tools are:

1. Analytic baseline curves for the main loop shape.
2. Bernstein polynomials or B-splines for bounded shape distributions.
3. Local orthonormal frames for three-dimensional blade construction.
4. Coupled surface definitions rather than independently generated, unrelated surfaces.
5. Explicit geometric feasibility checks.

The baseline geometry shall remain compact, with optional higher-order deformation parameters introduced only after an initial performance study.

---

## 3. Coordinate system and geometry definitions

Define:

- \(O=(0,0,0)\): rotor center.
- \(z\): rotational axis.
- \(x,y\): rotor plane.
- \(D_{\max}\): maximum permitted diameter.
- \(R_{\max}\): maximum actual radial extent of the solid.
- \(N_b\): integer number of blade units.
- \(u\in[0,1]\): parameter along a single U-shaped blade.
- \(v\in[0,1]\): parameter along the local chord direction.

The blade is constructed in the following order:

1. Generate the planar reference loop.
2. Convert the reference parameter into normalized arc length.
3. Add three-dimensional axial curvature.
4. Construct a continuously varying local coordinate frame.
5. Generate a chordwise width distribution.
6. Orient the chord using the pitch/twist function.
7. Apply an airfoil profile.
8. Join both blade roots to the hub.
9. Repeat the completed blade around the rotor axis.

The reference curve is defined as the inner planform edge of the blade.

The final blade is a finite-thickness solid, not merely an offset pair of two-dimensional curves.

---

## 4. Reference loop geometry

### 4.1 Baseline equation

The original proposed model was:

\[
r(\theta)=h\cos(a\theta).
\]

It supplies a useful family of smooth curves but couples angular opening and loop shape. When its domain is chosen so that the cosine reaches zero at both ends, both endpoints coincide with the polar origin.

The revised formulation separates these geometric controls.

Let:

\[
\theta(u)=\beta(2u-1)
\]

and

\[
\boxed{
r(u)=r_h+H\sin^p(\pi u)
}
\]

where:

- \(r_h\): radial location of the root attachment points.
- \(H\): radial extension of the loop relative to the root radius.
- \(\beta\): half-angle of the loop opening.
- \(p\): loop-shape exponent.

Then:

\[
x(u)=r(u)\cos\theta(u),
\]

\[
y(u)=r(u)\sin\theta(u).
\]

This gives the planar reference curve:

\[
\boxed{
\mathbf C_{2D}(u)=
\begin{bmatrix}
r(u)\cos\theta(u)\\
r(u)\sin\theta(u)
\end{bmatrix}
}
\]

### 4.2 Endpoint behavior

At \(u=0\):

\[
r(0)=r_h,\qquad \theta(0)=-\beta.
\]

At \(u=1\):

\[
r(1)=r_h,\qquad \theta(1)=+\beta.
\]

Therefore, the two attachment points are separated when \(r_h>0\) and \(\beta>0\).

At \(u=1/2\):

\[
r(1/2)=r_h+H.
\]

This is the radial maximum of the reference curve, provided \(H>0\).

The final solid can extend beyond this radius because of chord width and section thickness. Consequently, the diameter constraint must be applied to the finished solid, not only to the reference curve.

### 4.3 Shape exponent

Changing \(p\) changes the radial profile:

- \(p=1\): the original cosine-derived profile.
- \(p>1\): relatively more concentrated radial extension near the middle of the loop.
- Increasing \(p\): changes the curvature distribution and generally narrows the region near maximum radius.

The appropriate range should be established by geometry tests. An initial experimental range is:

\[
1\le p\le3.
\]

The generated curve is smooth over its parameter interval, but its tangents at the roots must still be connected smoothly to the hub.

### 4.4 Additional shape freedom

After validating the baseline model, use:

\[
r(u)=r_h+
H\sin^p(\pi u)
\left[
1+c_1(2u-1)+c_2(2u-1)^2
\right].
\]

Here:

- \(c_1\) introduces asymmetric deformation along the loop.
- \(c_2\) modifies the radial distribution more generally.

Initially:

\[
c_1=c_2=0.
\]

For later optimization, constrain the bracketed factor to remain strictly positive. For example, the conservative condition

\[
|c_1|+|c_2|\le 0.5
\]

ensures that the factor remains at least \(0.5\).

This condition does not guarantee that the finished blade is free of self-intersections. Those must still be checked independently.

### 4.5 Future spline replacement

A B-spline may ultimately replace the analytic radial profile if the cosine-derived family is insufficient.

The preferred approach is to use a small number of control points with constrained endpoints and smoothness, rather than increasing polynomial degree indefinitely.

The optimized geometry should be tested against the baseline to establish whether the added flexibility delivers statistically meaningful aerodynamic improvement.

---

## 5. Arc-length parameterization

Width, twist, and other shape distributions should ideally be defined using normalized physical arc length rather than the original polar parameter.

For the planar reference curve:

\[
L(u)=\int_0^u
\left\|
\frac{d\mathbf C_{2D}(\xi)}{d\xi}
\right\|d\xi.
\]

Define:

\[
\boxed{
s(u)=\frac{L(u)}{L(1)}
}
\]

where:

\[
s\in[0,1].
\]

This produces a dimensionless arc-length coordinate.

It allows width and shape coefficients to retain more consistent physical meaning when \(H\), \(\beta\), or \(p\) change.

The mapping may be calculated numerically using sampled points and cumulative arc length.

For functions controlling three-dimensional geometry, an alternative is to use arc length of the fully three-dimensional reference curve. The choice must remain consistent across all evaluations.

Version 1 will use normalized planar arc length for width controls and \(u\) for the analytic axial-rise function.

---

## 6. Blade width distribution

### 6.1 Replacing the unrestricted quartic

The original proposed width equation was:

\[
t(x)=-Ax^4+Bx^2+C.
\]

Its symmetry may exclude useful width distributions, and its coefficients do not independently control intuitive dimensions.

For the optimization system, distinguish the following quantities:

- \(w(s)\): local chord or planform width.
- \(h(s,v)\): physical thickness of the airfoil section.

These must not be represented by the same variable.

### 6.2 Bernstein polynomial representation

Use:

\[
\boxed{
w(s)=\sum_{i=0}^{m}
b_i\binom{m}{i}s^i(1-s)^{m-i}
}
\]

with \(m=2\) initially:

\[
\boxed{
w(s)=b_0(1-s)^2+
2b_1s(1-s)+
b_2s^2
}
\]

This requires three parameters.

The values represent the control coefficients for local chord width along the reference curve.

Important mathematical properties:

- The basis functions are nonnegative over \([0,1]\).
- Their sum is one.
- The width is therefore bounded by the minimum and maximum coefficients.
- Different coefficients allow asymmetric distributions.
- The function is smooth and relatively inexpensive to evaluate.

If:

\[
0<w_{\min}\le b_i\le w_{\max},
\]

then:

\[
\boxed{
w_{\min}\le w(s)\le w_{\max}
}
\]

for all \(s\in[0,1]\).

This provides a useful built-in positivity guarantee.

### 6.3 Higher-order extension

A quartic Bernstein polynomial with five coefficients may be used in Version 2.

Increase the degree only if a lower-order model demonstrably prevents useful variation.

To discourage unnecessary high-frequency geometry variations, introduce a smoothness penalty such as:

\[
E_w=\int_0^1
\left(\frac{d^2w}{ds^2}\right)^2ds.
\]

This penalty should not dominate aerodynamic performance.

### 6.4 Root width

The blade roots must have enough material for the expected loads and for attachment to the hub.

Do not impose zero root width unless the geometric construction explicitly transitions into a separate, structurally adequate attachment.

Specify:

\[
w(0)\ge w_{\text{root,min}},
\]

\[
w(1)\ge w_{\text{root,min}}.
\]

Additional root fillets and smooth transitions must be generated when the blade joins the hub.

---

## 7. Three-dimensional loop geometry

A complete toroidal blade requires a three-dimensional reference curve.

### 7.1 Axial-rise function

Use:

\[
\boxed{
z(u)=Z\sin^2(\pi u)
}
\]

where \(Z\) controls maximum axial displacement from the root attachment plane.

This gives:

\[
z(0)=z(1)=0,
\]

\[
z(1/2)=Z.
\]

The three-dimensional reference curve is:

\[
\boxed{
\mathbf C(u)=
\begin{bmatrix}
r(u)\cos\theta(u)\\
r(u)\sin\theta(u)\\
Z\sin^2(\pi u)
\end{bmatrix}
}
\]

This introduces axial curvature without introducing additional endpoint attachment positions.

### 7.2 Future axial deformation

A second coefficient may later be added:

\[
z(u)=
Z_1\sin^2(\pi u)
+
Z_2(2u-1)\sin^2(\pi u).
\]

The second term introduces an asymmetric axial displacement while preserving the original endpoint positions.

Initially:

\[
Z_2=0.
\]

### 7.3 Local moving frame

Calculate the normalized tangent:

\[
\mathbf T(u)=
\frac{\mathbf C'(u)}
{\|\mathbf C'(u)\|}.
\]

Define a planar normal to the projected reference curve:

\[
\mathbf e_1(u)=
\frac{
[-y'(u),x'(u),0]^T
}{
\sqrt{x'(u)^2+y'(u)^2}
}.
\]

Choose its sign consistently so that the chord extends into the intended blade region.

Define:

\[
\mathbf e_2(u)=
\mathbf T(u)\times\mathbf e_1(u).
\]

The vectors \(\mathbf T,\mathbf e_1,\mathbf e_2\) form a local orthonormal frame as long as the projected curve has a nonzero tangent.

This particular frame is suitable because the reference curve remains regular in its planar projection.

If a future geometry permits nearly vertical reference tangents, implement a more general rotation-minimizing frame to avoid numerical instability.

### 7.4 Pitch/twist distribution

Initially define:

\[
\boxed{
\alpha(u)=\alpha_0+
\alpha_1(2u-1)
}
\]

where:

- \(\alpha_0\): baseline chord orientation.
- \(\alpha_1\): linear twist coefficient.

The local chord direction is:

\[
\boxed{
\mathbf d(u)=
\cos\alpha(u)\,\mathbf e_1(u)
+
\sin\alpha(u)\,\mathbf e_2(u)
}
\]

The definition of pitch is geometric: it rotates the chord about the local tangent.

It is not automatically equal to the aerodynamic angle of attack.

That angle depends on the local relative airflow, induced velocity, and section orientation.

The baseline geometry must be inspected to ensure that its leading-edge orientation is appropriate on both sides of the U-shaped loop.

### 7.5 Blade mid-surface

Let \(v\in[0,1]\) be the local chord coordinate.

Then:

\[
\boxed{
\mathbf S(u,v)=
\mathbf C(u)+
v\,w(s(u))\,\mathbf d(u)
}
\]

This generates the three-dimensional blade mid-surface.

The inner reference edge is:

\[
\mathbf S(u,0)=\mathbf C(u).
\]

The opposite planform edge is:

\[
\mathbf S(u,1)=
\mathbf C(u)+w(s(u))\mathbf d(u).
\]

The result is a single continuous, swept blade surface.

The outer blade edge is therefore coupled directly to the reference curve, width distribution, and pitch.

Independent polynomials for the two planform edges are not required.

---

## 8. Airfoil section and physical thickness

The blade mid-surface is not yet a physical solid.

A finite-thickness aerodynamic section must be constructed around it.

### 8.1 Initial airfoil choice

For Version 1, use a fixed conventional, symmetric airfoil profile as a geometric baseline, such as NACA 0012.

This is a convenient initial test profile, not a claim that it is optimal for the proposed rotor.

Its local chord length is:

\[
c(u)=w(s(u)).
\]

A constant nominal thickness-to-chord ratio may be used initially.

The choice must remain appropriate for the expected Reynolds number, operating conditions, and manufacturing process.

### 8.2 Mid-surface normal

Calculate:

\[
\mathbf n_S(u,v)=
\frac{
\partial_u\mathbf S\times
\partial_v\mathbf S
}{
\|\partial_u\mathbf S\times
\partial_v\mathbf S\|
}.
\]

Apply the airfoil's upper and lower normal displacements to the mid-surface.

For a symmetric section with dimensionless thickness function \(f_t(v)\):

\[
\mathbf S_{\pm}(u,v)=
\mathbf S(u,v)
\pm
\frac{c(u)\tau(u)f_t(v)}{2}
\mathbf n_S(u,v).
\]

Here \(\tau(u)\) is the section's thickness-to-chord ratio.

A cambered airfoil additionally requires a camber displacement function.

### 8.3 General camber/thickness formulation

Define dimensionless camber \(f_c(v)\) and thickness \(f_t(v)\).

The upper and lower surfaces are:

\[
\mathbf S_{\pm}(u,v)=
\mathbf S(u,v)+
c(u)
\left[
f_c(v)\pm
\frac{\tau(u)f_t(v)}{2}
\right]
\mathbf n_S(u,v).
\]

This is a local surface-normal construction.

For heavily twisted, curved, or thick geometries, the resulting surfaces must still be checked for offset singularities, crossings, and section distortion.

A production implementation may instead generate explicit two-dimensional airfoil coordinates and map them into each local frame before lofting.

### 8.4 Future CST representation

For airfoil optimization, adopt Class-Shape Transformation (CST), which combines a class function with a Bernstein-polynomial shape function.

A standard form is:

\[
y(x)=
x^{N_1}(1-x)^{N_2}
\sum_{i=0}^{m}
a_iB_{i,m}(x),
\]

where \(x\in[0,1]\) is normalized chord position.

A conventional airfoil class uses:

\[
N_1=\frac12,\qquad N_2=1.
\]

Upper and lower airfoil geometries may be represented separately, with appropriate leading-edge, trailing-edge, and minimum-thickness constraints.

CST coefficients should be optimized only after the principal three-dimensional blade parameters are reasonably established.

### 8.5 Surface completion

The geometry generator must additionally create:

- Physical leading-edge and trailing-edge connections.
- End/root surfaces.
- Smooth blade-to-hub transitions.
- Appropriate fillets.
- A closed, consistently oriented surface mesh.

A blade cannot be accepted for CFD simply because its mathematical mid-surface is smooth.

---

## 9. Rotational replication

The first blade unit is generated at a reference angular position.

For blade index:

\[
k=0,1,\ldots,N_b-1,
\]

define:

\[
\phi_k=\frac{2\pi k}{N_b}.
\]

The three-dimensional rotation matrix is:

\[
\boxed{
\mathbf R_z(\phi_k)=
\begin{bmatrix}
\cos\phi_k&-\sin\phi_k&0\\
\sin\phi_k&\cos\phi_k&0\\
0&0&1
\end{bmatrix}
}
\]

All points belonging to the blade are transformed using:

\[
\mathbf X_k=
\mathbf R_z(\phi_k)\mathbf X_0.
\]

The repeated blade solids are then combined with the hub.

The resulting rotor must be validated for collisions between adjacent blades and unintended geometry overlaps.

The optimizer should not assume that angular sectors are collision-free merely because the centerline opening is smaller than \(2\pi/N_b\). Actual blade width, twist, and airfoil thickness may cause adjacent solids to intersect.

### Blade count

Optimize a small discrete set of blade counts separately, such as:

\[
N_b\in\{2,3,4\}.
\]

This range is an initial study configuration, not a physical limitation.

Each blade count defines a different feasible design space.

Separate optimization runs are easier to interpret than including blade count in the initial continuous SA proposal mechanism.

---

## 10. Initial design parameter vector

For Version 1, fix the hub dimensions, blade count within each run, and airfoil profile.

Optimize:

\[
\boxed{
\mathbf p=
[
H,\beta,p,
b_0,b_1,b_2,
Z,\alpha_0,\alpha_1
]^T
}
\]

This gives **nine continuous design parameters**.

| Parameter | Meaning | Units |
|---|---|---|
| \(H\) | Radial loop extension | m |
| \(\beta\) | Loop half-opening | rad |
| \(p\) | Radial shape exponent | Dimensionless |
| \(b_0\) | Chord-width control coefficient 0 | m |
| \(b_1\) | Chord-width control coefficient 1 | m |
| \(b_2\) | Chord-width control coefficient 2 | m |
| \(Z\) | Maximum axial loop rise | m |
| \(\alpha_0\) | Baseline chord orientation | rad |
| \(\alpha_1\) | Linear twist coefficient | rad |

The fixed hub-root radius \(r_h\) may become a tenth continuous parameter in a later version.

### 10.1 Proposed starting search ranges

Use the following as provisional, dimensionless exploration ranges only:

| Parameter | Exploratory range |
|---|---|
| \(H/D_{\max}\) | 0.20 to 0.35 |
| \(\beta\) | \(0.2\pi/N_b\) to \(0.85\pi/N_b\) |
| \(p\) | 1.0 to 3.0 |
| \(b_i/D_{\max}\) | 0.025 to 0.12 |
| \(Z/D_{\max}\) | 0.00 to 0.20 |
| \(\alpha_0\) | 0° to 45° |
| \(\alpha_1\) | -20° to 20° |

These are initialization guidelines, not validated aerodynamic limits.

The actual interval values must be adjusted to the rotor's size, root radius, manufacturing limitations, and reference blade.

A candidate satisfying these individual parameter ranges is not automatically valid.

### 10.2 Parameter normalization

Normalize each optimization variable:

\[
q_i=
\frac{p_i-p_{i,\min}}
{p_{i,\max}-p_{i,\min}}.
\]

Then:

\[
q_i\in[0,1].
\]

SA will operate on the normalized vector \(\mathbf q\).

This prevents parameters measured in radians, meters, and dimensionless coefficients from being perturbed at incomparable scales.

---

## 11. Geometric and mechanical constraints

Geometry validation must occur before CFD.

### 11.1 Mandatory geometric constraints

Reject candidates when any of the following occurs:

1. Non-finite coordinates or derivatives.
2. Zero-length or numerically degenerate tangent vectors.
3. Invalid airfoil sections.
4. Negative or excessively small width.
5. Self-intersecting surfaces.
6. Intersections between separate blades outside their intended hub connection.
7. Non-manifold solid geometry.
8. Open edges in the final exported solid.
9. Inconsistent surface orientation.
10. Failure to satisfy maximum diameter.
11. Insufficient root attachment dimensions.
12. Violations of required minimum solid thickness.
13. Failed hub unions or invalid fillet construction.
14. Local section offsets that collapse or fold the blade surface.

### 11.2 Curvature constraints

Compute curvature for the reference curve:

\[
\kappa(u)=
\frac{
\|\mathbf C'(u)\times\mathbf C''(u)\|
}{
\|\mathbf C'(u)\|^3
}.
\]

Large curvature relative to chord and thickness is a warning for offset-surface folding, poor mesh quality, and stress concentration.

Introduce curvature constraints appropriate to the local geometry and manufacturing method.

A universal curvature limit is not recommended.

### 11.3 Structural constraints

For mechanically promising candidates, evaluate:

- Blade mass.
- Centrifugal loading.
- Root reaction loads.
- Peak material stress.
- Deflection.
- Resonance risk.
- Fatigue sensitivity.
- Hub attachment integrity.

The centrifugal body-force magnitude per unit mass is locally:

\[
a_c=\Omega^2r.
\]

A first approximation to the total centrifugal-force magnitude is:

\[
F_c\sim\int_{\text{blade}}\Omega^2r\,dm.
\]

The vector directions and internal load path must be evaluated from the actual geometry.

For later structural analysis, include distributed aerodynamic surface pressures and centrifugal body forces.

The rotor should also be checked for imbalance, especially after manufacturing.

### 11.4 Manufacturability

For 3D-printed prototype geometries, introduce:

- Minimum printable wall thickness.
- Minimum printable edge radius.
- Root reinforcement.
- Print orientation considerations.
- Support-removal accessibility.
- Dimensional tolerance.
- Surface finishing requirements.

Manufacturing constraints should be configurable independently of the aerodynamic parameterization.

No numerical stress limit or minimum thickness should be chosen without specifying the material and fabrication method.

---

## 12. Computational fluid dynamics

### 12.1 Governing physics

For low-speed airflow that can reasonably be treated as incompressible:

\[
\nabla\cdot\mathbf U=0.
\]

The momentum equation in an inertial frame is:

\[
\rho
\left(
\frac{\partial\mathbf U}{\partial t}
+
\mathbf U\cdot\nabla\mathbf U
\right)
=
-\nabla p+
\nabla\cdot\boldsymbol\tau.
\]

For Newtonian, constant-viscosity, incompressible flow:

\[
\boldsymbol\tau=
\mu
\left(
\nabla\mathbf U+
\nabla\mathbf U^T
\right).
\]

Turbulent simulations require an additional turbulence closure. In a rotating-reference-frame formulation, the governing momentum equation is transformed and appropriate rotational source terms are included.

### 12.2 Compressible-flow extension

If compressibility becomes important, use a compressible solver with:

Mass conservation:

\[
\frac{\partial\rho}{\partial t}
+
\nabla\cdot(\rho\mathbf U)=0.
\]

Momentum conservation:

\[
\frac{\partial(\rho\mathbf U)}{\partial t}
+
\nabla\cdot(\rho\mathbf U\otimes\mathbf U)
=
-\nabla p+\nabla\cdot\boldsymbol\tau.
\]

Total-energy conservation:

\[
\frac{\partial(\rho E)}{\partial t}
+
\nabla\cdot[(\rho E+p)\mathbf U]
=
\nabla\cdot(\boldsymbol\tau\cdot\mathbf U)
-\nabla\cdot\mathbf q.
\]

For a calorically perfect gas:

\[
p=\rho RT,
\]

\[
E=e+\frac12|\mathbf U|^2.
\]

The turbulence model contributes additional modeled stresses and transport equations.

### 12.3 Flow-regime selection

Evaluate characteristic Reynolds and tip Mach numbers:

\[
Re=\frac{\rho U_{\text{rel}}c}{\mu},
\]

\[
M_{\text{tip}}\approx
\frac{\sqrt{(\Omega R_{\max})^2+V_\infty^2}}
{a_\infty}.
\]

The incompressible approximation may be reasonable for low Mach numbers; \(M\lesssim0.3\) is a common initial screening heuristic rather than a universal guarantee.

The transition between laminar, transitional, and turbulent behavior depends on local Reynolds number, pressure gradients, turbulence intensity, and blade geometry.

A fully turbulent model should not automatically be assumed accurate for low-Reynolds-number rotor sections.

### 12.4 Operating points

Support at least:

**Static-thrust operation**

\[
V_\infty=0.
\]

**Forward-flight operation**

\[
V_\infty>0.
\]

For Version 1, choose one operating point and keep it constant throughout a given optimization experiment.

Later, support multi-point optimization across RPM and incoming flow velocity.

---

## 13. OpenFOAM implementation

### 13.1 Baseline analysis method

Use a steady rotating-reference-frame or Multiple Reference Frame (MRF) approach as the initial screening method.

This avoids physically rotating the computational mesh during every time step.

MRF is useful for inexpensive approximate assessment, but it cannot reproduce every time-dependent blade-wake interaction.

It should not be treated as a substitute for transient validation of the final designs.

### 13.2 Higher-fidelity analysis

For selected candidates, use transient rotating-mesh simulations with a suitable sliding or arbitrary mesh interface.

A PIMPLE-based solver is one practical route, depending on the installed OpenFOAM distribution and version.

The exact solver and case dictionaries must be selected after confirming the OpenFOAM version and the intended compressibility model.

### 13.3 Turbulence model

A practical starting point is \(k\)-\(\omega\) SST.

It is commonly used for aerodynamic flows with adverse pressure gradients and possible separation.

However, a loop-bladed rotor may exhibit significant three-dimensional vortical flow, separated wakes, and transition effects. The model must therefore be validated against mesh refinement and, ideally, experiment.

Consider alternative transition models or higher-fidelity turbulence simulation only when the results justify the added cost.

### 13.4 Computational domain

Use a sufficiently large outer domain to reduce artificial boundary effects.

An initial full-rotor domain-sizing experiment may start with:

- Radial boundary several rotor diameters from the axis.
- Upstream boundary several diameters ahead of the rotor.
- Downstream boundary substantially farther away to capture the wake.

Illustrative starting distances are roughly 3–5 diameters for upstream/radial boundaries and 8–10 diameters downstream.

These are not universal CFD requirements. Domain-size convergence must be established.

Static-thrust conditions require particular attention to entrainment and pressure boundaries. A boundary condition that suppresses inward flow can distort the predicted thrust.

### 13.5 Meshing pipeline

Use:

1. Watertight STL or another supported geometry representation.
2. `blockMesh` for the background mesh.
3. Feature extraction as required.
4. `snappyHexMesh` for geometry-conforming mesh generation.
5. Local blade, tip, hub, and wake refinement.
6. Prism layers where required for boundary-layer resolution.
7. `checkMesh` and other mesh-quality validation.

`snappyHexMesh` supports automated processing of triangulated geometry and parallel meshing, making it suitable for this optimization workflow.

### 13.6 Wall resolution

Choose a consistent wall treatment.

For a wall-resolved SST calculation, aim for first-cell \(y^+\) near or below one, with sufficient boundary-layer resolution.

For wall-function calculations, use a mesh and boundary conditions compatible with the chosen wall-function treatment.

Do not combine incompatible wall-function assumptions and first-cell spacing.

Use the OpenFOAM `yPlus` function object to verify the actual result.

### 13.7 Mesh convergence

At minimum, evaluate the final baseline and selected optimized shapes with:

- Coarse mesh.
- Medium mesh.
- Fine mesh.

Record:

\[
T_{\text{coarse}},
\quad T_{\text{medium}},
\quad T_{\text{fine}},
\]

and corresponding torque and power values.

A provisional acceptance target is a small change, for example 2–3%, in both thrust and torque between the two finest meshes.

This is an engineering screening threshold, not formal proof of grid independence.

Where results are sensitive, estimate discretization uncertainty more rigorously and refine further.

### 13.8 Convergence criteria

Do not declare a CFD run successful merely because its residuals decreased.

Require, at minimum:

- No solver divergence or numerical instability.
- Sufficient reduction or stabilization of residuals.
- Stable thrust and torque histories.
- Acceptable continuity and conservation errors.
- Consistent mesh quality.
- No unresolved large-scale numerical oscillations masquerading as convergence.

For transient simulations, compute averages over a sufficiently developed, repeatable interval and report fluctuations.

---

## 14. Aerodynamic outputs

### 14.1 Forces and torque

Integrate pressure and viscous traction over the solid rotor surfaces.

Extract:

\[
T=\mathbf F\cdot\mathbf e_z.
\]

Calculate shaft-axis torque:

\[
Q=\mathbf M\cdot\mathbf e_z,
\]

using a consistent sign convention for rotation and shaft resistance.

Use the magnitude of resisting torque where appropriate when computing consumed shaft power.

The shaft power is:

\[
\boxed{
P=Q_{\text{resist}}\Omega
}
\]

where:

\[
\Omega=\frac{2\pi\,RPM}{60}.
\]

### 14.2 Nondimensional coefficients

Define:

\[
n=\frac{RPM}{60}.
\]

The thrust coefficient is:

\[
\boxed{
C_T=\frac{T}{\rho n^2D^4}
}
\]

The shaft-power coefficient is:

\[
\boxed{
C_P=\frac{P}{\rho n^3D^5}
}
\]

The advance ratio is:

\[
\boxed{
J=\frac{V_\infty}{nD}.
}
\]

Here \(D\) is the actual rotor diameter for the reported design.

### 14.3 Propulsive efficiency

For forward flight with useful positive thrust:

\[
\boxed{
\eta=\frac{TV_\infty}{P}
=
\frac{JC_T}{C_P}.
}
\]

At zero freestream velocity, this conventional propulsive-efficiency expression is zero and is not a useful hover-performance objective.

### 14.4 Static-thrust performance

For static operation, the primary comparison is shaft power required to achieve a specified thrust.

An additional metric is the conventional hover figure of merit:

\[
FM=
\frac{T^{3/2}}
{P\sqrt{2\rho A}},
\]

where:

\[
A=\pi(D/2)^2.
\]

This compares actual shaft power with ideal momentum-theory induced power for the same thrust and disk area.

Use it as a comparative metric, with the limitations of ideal actuator-disk theory understood.

### 14.5 Additional outputs

Record:

- Thrust.
- Torque.
- Shaft power.
- Thrust coefficient.
- Power coefficient.
- Actual rotor diameter.
- Blade surface area.
- Rotor solid volume.
- Estimated rotor mass.
- Axial and tangential wake velocities.
- Blade surface pressure distribution.
- Blade loading distribution.
- Maximum and representative \(y^+\).
- Mesh statistics.
- Solver convergence status.
- Computational cost.

For later acoustic optimization, store unsteady pressure and loading histories.

---

## 15. Simulated annealing

### 15.1 Optimization representation

SA will modify the normalized parameter vector:

\[
\mathbf q\in[0,1]^9.
\]

Each candidate is converted to physical design parameters, then passed into the geometry-validation and CFD pipeline.

### 15.2 Proposal generation

For a continuous variable:

\[
q_i'=q_i+\sigma_i\epsilon_i,
\]

where:

\[
\epsilon_i\sim\mathcal N(0,1).
\]

Use reflected or otherwise bounded proposals when a candidate leaves its permitted interval.

The proposal distribution should include:

- Single-parameter perturbations.
- Small, correlated multiparameter perturbations.
- Occasional larger exploratory jumps.

Initially, select each variable at random and use a modest proposal scale relative to its normalized range.

Adjust proposal scales based on observed acceptance rates, without modifying the target acceptance rule incorrectly.

### 15.3 Penalized objective

During optimization, an infeasible thrust value should not automatically end the search.

Use a penalty such as:

\[
\boxed{
F(\mathbf p)=
\frac{P(\mathbf p)}{P_{\text{ref}}}
+
\lambda_T
\left[
\max\left(
0,1-\frac{T(\mathbf p)}{T_{\text{target}}}
\right)
\right]^2
}
\]

where \(P_{\text{ref}}\) is a fixed reference shaft power.

Optional terms can penalize excessive mass, noise, or other soft constraints.

The penalty coefficient \(\lambda_T\) must be calibrated so that the optimizer prioritizes satisfying the target thrust while still distinguishing useful improvements.

Strictly invalid geometry, such as intersecting blade surfaces, should be rejected before CFD rather than evaluated with an arbitrary aerodynamic penalty.

Final selection must use the original physical feasibility constraints, regardless of the numerical penalty score.

### 15.4 Acceptance rule

Let:

\[
\Delta F=F(\mathbf p')-F(\mathbf p).
\]

Accept improving candidates automatically:

\[
\Delta F\le0.
\]

For worse candidates, accept with probability:

\[
\boxed{
P_{\text{accept}}=
\exp\left(-\frac{\Delta F}{\Theta}\right)
}
\]

where \(\Theta\) is the artificial SA temperature.

The artificial temperature has no connection with fluid or material temperature.

### 15.5 Cooling schedule

A possible geometric cooling schedule is:

\[
\Theta_k=\Theta_0\gamma^k,
\]

with:

\[
0<\gamma<1.
\]

Choose \(\Theta_0\) from an initial distribution of objective-score differences so that the initial acceptance probability for representative uphill moves is sufficiently high.

A value of \(\gamma\) near one produces slower cooling and more evaluations.

Cooling should be based on completed candidate evaluations rather than solver iterations.

### 15.6 Multiple restarts

Run SA from:

- The original baseline geometry.
- Several feasible, randomly generated geometries.
- Previously discovered promising designs.

For each blade count, retain the best feasible solution across independent runs.

This reduces dependence on a single initial configuration.

### 15.7 Comparison against other optimizers

SA is the initial preferred optimization algorithm, not the only algorithm permitted.

Benchmark it against:

- Differential evolution.
- CMA-ES.
- Bayesian optimization when CFD evaluations are expensive.

Use the same evaluation budget and objective definitions for comparisons.

The final optimizer should be chosen from measured performance rather than familiarity.

---

## 16. Multi-fidelity optimization

Full-resolution CFD for every proposed SA candidate is likely too expensive.

Use three evaluation levels.

### Level 0: Geometry screening

Very low cost.

Evaluate:

- Mathematical parameter validity.
- Surface intersections.
- Watertightness.
- Blade clearance.
- Root geometry.
- Diameter limits.
- Manufacturability limits.
- Approximate geometric curvature.

Reject candidates failing mandatory geometric constraints.

### Level 1: Low-fidelity aerodynamic screening

Use a repeatable, relatively coarse CFD configuration.

Possible simplifications include:

- Coarser meshes.
- Lower computational cost.
- Steady MRF.
- Reduced output frequency.
- Limited iterative convergence targets, provided the resulting uncertainty is characterized.

Traditional blade-element momentum methods may not adequately model the highly coupled three-dimensional flow around loop-shaped blades.

They should not automatically be trusted as the low-fidelity model for toroidal geometry.

The initial low-fidelity model should therefore be a coarse version of a validated CFD setup.

### Level 2: High-fidelity evaluation

For promising candidates:

- Refine the mesh.
- Improve near-wall resolution.
- Increase numerical convergence.
- Evaluate transient flow where necessary.
- Quantify numerical uncertainty.
- Verify thrust and torque stability.

### Surrogate assistance

Once enough evaluated candidates are available, train a model predicting:

\[
\widehat T(\mathbf p),
\qquad
\widehat P(\mathbf p).
\]

Useful methods include Gaussian processes, radial-basis models, and ensemble regression.

Train on confirmed CFD data.

Use the surrogate to rank candidates or estimate uncertainty, but do not substitute unvalidated predictions for final CFD results.

Maintain an exploration mechanism so that the optimizer does not focus only on regions favored by a potentially inaccurate surrogate.

---

## 17. CFD workflow and experiment control

Every candidate must receive a unique identifier.

Each evaluation is performed as a reproducible operation.

### 17.1 Candidate evaluation

1. Receive parameter vector.
2. Convert normalized variables to physical values.
3. Generate reference curve.
4. Generate axial-rise and twist functions.
5. Construct blade surfaces.
6. Construct airfoil geometry and hub transitions.
7. Generate rotationally repeated blade units.
8. Validate solid geometry.
9. Export mesh-ready geometry.
10. Generate or update an OpenFOAM case.
11. Create CFD mesh.
12. Validate CFD mesh.
13. Run solver.
14. Monitor convergence.
15. Extract thrust and torque.
16. Calculate aerodynamic performance.
17. Store logs, inputs, geometry, and results.
18. Return objective score and evaluation status.

### 17.2 Failure handling

Define explicit evaluation states:

- `GEOMETRY_INVALID`
- `EXPORT_FAILED`
- `MESH_FAILED`
- `SOLVER_FAILED`
- `NOT_CONVERGED`
- `SUCCESS`

A failed evaluation must not be interpreted as zero power, zero thrust, or a physically meaningful aerodynamic result.

Record the failure reason and, where practical, retry only numerical failures under a controlled fallback configuration.

### 17.3 Caching

Cache results using a stable identifier based on:

- Exact geometry parameters.
- Geometry-generator version.
- CFD configuration.
- Operating conditions.
- Mesh configuration.
- Solver version.
- Numerical settings.

The cache must not silently treat physically different geometries as identical.

### 17.4 Reproducibility

Record:

- SA random seed.
- Starting vector.
- Proposal history.
- Candidate acceptance decisions.
- Cooling schedule.
- CFD case settings.
- Software versions.
- Solver logs.
- Mesh quality data.

Identical inputs and software settings should reproduce the same geometry and nominal evaluation procedure.

Numerical results may still differ slightly across hardware, parallel decompositions, or solver environments.

---

## 18. Parallel computing

The system shall distinguish two kinds of parallelism.

### 18.1 Candidate-level parallelism

Multiple independent geometry/CFD evaluations may run simultaneously.

This is useful for:

- Different blade counts.
- Independent SA restarts.
- High-fidelity validation.
- Surrogate training data generation.

Standard SA is sequential because each acceptance decision depends on the current state. Parallel candidate evaluation therefore requires a deliberate asynchronous or batched optimization design.

### 18.2 Solver-level parallelism

One CFD case may be decomposed into subdomains and evaluated using MPI.

OpenFOAM supports decomposed parallel cases and parallel meshing.

This permits the simulation memory to be distributed over multiple machines, although communication, mesh interfaces, and data movement impose overhead.

Memory usage and CPU efficiency will not scale perfectly with the number of machines.

### 18.3 Scheduling

Build a resource-aware scheduler that tracks:

- Available cores.
- Available RAM.
- Simulation queue.
- Running jobs.
- Estimated computational requirements.
- Completed evaluations.
- Failed jobs.

Avoid running too many CFD jobs simultaneously and causing memory exhaustion or severe CPU contention.

---

## 19. Software architecture

Use Python as the main orchestration and geometry-definition language.

Suggested project structure:

```text
propeller-optimization/
├── config/
│   ├── design.yaml
│   ├── constraints.yaml
│   ├── operating_point.yaml
│   └── cfd.yaml
├── geometry/
│   ├── loop.py
│   ├── arclength.py
│   ├── width.py
│   ├── frames.py
│   ├── twist.py
│   ├── airfoil.py
│   ├── loft.py
│   ├── hub.py
│   ├── assembly.py
│   └── validation.py
├── simulation/
│   ├── case_builder.py
│   ├── mesh.py
│   ├── runner.py
│   ├── parser.py
│   └── convergence.py
├── optimization/
│   ├── objective.py
│   ├── annealing.py
│   ├── proposals.py
│   ├── surrogate.py
│   └── benchmark.py
├── experiments/
│   ├── database.py
│   ├── scheduler.py
│   └── analysis.py
├── templates/
│   └── openfoam/
├── tests/
├── runs/
└── main.py
```

### 19.1 Configuration

Use structured configuration files to specify:

- Rotor diameter.
- Hub geometry.
- Blade-count candidates.
- Fluid properties.
- Operating RPM.
- Incoming velocity.
- Target thrust.
- Material properties.
- Geometric limits.
- Optimization parameters.
- CFD fidelity.

All geometry dimensions should use SI units internally.

### 19.2 Geometry generator

Core interface:

```python
geometry = generate_rotor(parameters, configuration)
```

The result should include:

- Reference curves.
- Blade mid-surface.
- Complete rotor surface mesh.
- Optional CAD export.
- Surface-region metadata.
- Geometric measurements.

### 19.3 Validation

Core interface:

```python
validation = validate_geometry(
    geometry,
    constraints
)
```

The result should return a detailed status and structured constraint violations.

### 19.4 CFD

Core interface:

```python
result = evaluate_cfd(
    geometry,
    operating_point,
    fidelity
)
```

It should return thrust, torque, power, coefficients, solver status, and diagnostic information.

### 19.5 Optimization

Core interface:

```python
best = optimize(
    initial_parameters,
    parameter_bounds,
    objective,
    optimizer_config
)
```

Geometry generation and CFD evaluation should remain independent of the optimization algorithm.

This permits replacement of SA without rewriting the CAD or simulation system.

---

## 20. Development and verification plan

### Phase A: Geometry-only prototype

Implement:

- Parametric U-shaped curve.
- Width distribution.
- Axial-rise function.
- Local moving frame.
- Chord orientation.
- Rotor replication.
- Three-dimensional surface visualization.

Verify:

- Endpoint positions.
- Symmetry.
- Parameter influence.
- Surface continuity.
- Correct rotation.
- Proper physical dimensions.

No CFD is required during this phase.

**Exit criterion:** The system generates a range of geometrically coherent blades reproducibly.

### Phase B: Solid geometry

Implement:

- Fixed airfoil mapping.
- Solid thickening or section lofting.
- Hub union.
- Blade-root blending.
- Surface triangulation.
- Collision and watertightness tests.

**Exit criterion:** The generator produces valid rotor solids that can be meshed automatically.

### Phase C: CFD baseline

Build one fixed OpenFOAM case.

Verify:

- Stable meshing.
- Correct rotor coordinate system.
- Correct MRF angular velocity.
- Correct thrust and torque sign conventions.
- Mesh and domain-size sensitivity.
- Conservation and solver convergence.

Evaluate at least one conventional propeller reference and one loop-bladed baseline where practical.

**Exit criterion:** Consistent aerodynamic metrics can be extracted repeatedly.

### Phase D: Initial optimization

Run SA with the nine proposed continuous parameters.

Keep blade count fixed within each run.

Start from a validated baseline.

Evaluate several independent runs and compare their best feasible outcomes.

**Exit criterion:** The system finds valid candidates with better evaluated objective values than the original baseline, or establishes that no improvement was found within the explored budget.

### Phase E: Parameterization comparison

Add one or two asymmetric loop-shape coefficients.

Compare:

- Baseline cosine-derived geometry.
- Asymmetric analytic geometry.
- Optional spline-based geometry.

Use the same CFD fidelity and computational budget.

**Exit criterion:** Determine whether extra shape freedom provides an improvement larger than the observed numerical uncertainty.

### Phase F: High-fidelity confirmation

For the best candidates:

- Refine CFD meshes.
- Perform transient simulations where appropriate.
- Recompute performance metrics.
- Evaluate structural integrity.
- Test parameter sensitivity.
- Reject candidates whose claimed improvement disappears under higher fidelity.

**Exit criterion:** Select a physically plausible final rotor design with credible numerical performance evidence.

### Phase G: Experimental validation

Manufacture prototypes and measure performance on an instrumented test rig.

Collect:

- RPM.
- Static thrust.
- Shaft torque if available.
- Mechanical shaft power.
- Electrical input power.
- Vibration.
- Acoustic measurements where possible.

Electrical input power is not identical to shaft power and must not be compared directly with CFD shaft power without accounting for motor and drivetrain losses.

Use suitable containment, remote operation, and overspeed precautions for rotating prototypes.

Compare experiment and CFD at matched operating conditions.

**Exit criterion:** Quantify model error and update the optimization pipeline based on measured performance.

---

## 21. Computational experiment design

### 21.1 Initial optimization comparison

Run separate studies for:

\[
N_b=2,\quad N_b=3,\quad N_b=4.
\]

Use comparable CFD budgets for each blade count.

Store both the best feasible objective and the distribution of obtained results across restarts.

### 21.2 Initial computational budget

A provisional pilot experiment may include:

- Geometry validation of hundreds or thousands of sampled parameter vectors.
- A smaller set of successful low-fidelity CFD evaluations.
- Several independent optimization restarts.
- A small high-fidelity finalist set.

The precise budget depends on measured CFD runtime.

Do not commit to thousands of full CFD solves before establishing runtime, failure rate, and numerical uncertainty.

### 21.3 Sensitivity study

Perturb parameters individually around a validated baseline.

Estimate:

\[
S_i\approx
\frac{
F(\mathbf p+\Delta p_i\mathbf e_i)
-
F(\mathbf p-\Delta p_i\mathbf e_i)
}{
2\Delta p_i
}.
\]

For consistent interpretation, normalize the perturbation size and objective.

This is a local finite-difference sensitivity estimate, not proof of global parameter importance.

Also test selected parameter interactions.

The study should identify which degrees of freedom affect performance substantially and which may be safely fixed.

### 21.4 Robustness study

Test promising designs under small changes in:

- RPM.
- Incoming velocity.
- Geometry dimensions.
- Airfoil thickness.
- Surface roughness where modeling permits.

For uncertain inputs, consider an expected-performance objective or a worst-case penalty.

This reduces the risk of selecting an extremely sensitive design that performs well only at a single exact numerical condition.

---

## 22. Reporting requirements

Each optimization report should contain:

1. Project configuration and operating conditions.
2. Baseline geometry.
3. Optimization algorithm and computational budget.
4. Parameter bounds.
5. Number of valid and invalid candidates.
6. Mesh and solver failure statistics.
7. Best candidate parameter vector.
8. Baseline and optimized aerodynamic metrics.
9. Optimization score history.
10. Mesh sensitivity and convergence results.
11. Thrust and torque histories.
12. Surface pressure and wake visualizations.
13. Structural and manufacturing considerations.
14. High-fidelity validation results.
15. Numerical uncertainty and limitations.

Plot, where useful:

- Objective versus CFD evaluation number.
- Thrust versus shaft power.
- Parameter distributions.
- Parameter sensitivity.
- CFD convergence.
- Rotor geometry comparisons.
- Pressure and wake fields.

Store machine-readable results as well as human-readable reports.

---

## 23. Risks and mitigation

| Risk | Mitigation |
|---|---|
| Invalid geometry dominates proposals | Enforce positivity, bounds, collision checks, and coupled geometry construction |
| Excessively restrictive cosine model | Add asymmetric shape coefficients or splines after baseline testing |
| Too many optimization parameters | Use staged optimization and parameter sensitivity analysis |
| Unrealistic CFD improvements | Use high-fidelity mesh and model checks |
| CFD too expensive for SA | Introduce low-fidelity screening, caching, and surrogates |
| Numerical noise changes rankings | Re-evaluate finalists and quantify uncertainty |
| Low-Re turbulence mismatch | Check flow regime and compare suitable turbulence/transition models |
| Blade-root mechanical failure | Include root reinforcement and structural screening |
| Geometry valid but impossible to manufacture | Enforce process-specific feature constraints |
| Incorrect rotor pitch orientation | Inspect leading edges and local relative airflow before optimization |
| Optimizer exploits a bad numerical model | Validate independently against reference cases and experiments |
| Different designs compared unfairly | Standardize operating points and constraints |

---

## 24. Acceptance criteria

### Geometry

- Every accepted candidate satisfies the geometry validator.
- No unintended self-intersections.
- Closed, consistently oriented rotor surface mesh.
- Correct blade-count symmetry.
- Correct physical units and diameter constraints.

### CFD

- Baseline successfully meshes and runs.
- Reliable thrust and torque extraction.
- Appropriate convergence checks.
- Demonstrated sensitivity to mesh and domain size.
- High-fidelity validation of selected candidates.

### Optimization

- Deterministic geometry generation from fixed parameters.
- Repeatable optimization runs using recorded random seeds and settings.
- Multiple starts supported.
- Invalid candidates rejected before expensive CFD.
- Best designs verified against the target thrust and all hard constraints.
- Improvements compared against numerical uncertainty.

### Engineering

- Structural constraints assessed before manufacture.
- Rotor balance and operational safety considered.
- Experimental measurements used to evaluate the final simulation-based claims.

---

## 25. Final technical decision

The initial implementation will use:

**Main geometry**

\[
\theta(u)=\beta(2u-1),
\]

\[
r(u)=r_h+H\sin^p(\pi u),
\]

\[
z(u)=Z\sin^2(\pi u).
\]

**Width distribution**

\[
w(s)=b_0(1-s)^2+
2b_1s(1-s)+b_2s^2.
\]

**Twist**

\[
\alpha(u)=\alpha_0+
\alpha_1(2u-1).
\]

**Three-dimensional reference curve**

\[
\mathbf C(u)=
\begin{bmatrix}
r(u)\cos\theta(u)\\
r(u)\sin\theta(u)\\
z(u)
\end{bmatrix}.
\]

**Chord direction**

\[
\mathbf d(u)=
\cos\alpha(u)\mathbf e_1(u)+
\sin\alpha(u)\mathbf e_2(u).
\]

**Blade mid-surface**

\[
\boxed{
\mathbf S(u,v)=
\mathbf C(u)+
v\,w(s(u))\mathbf d(u)
}
\]

where:

\[
u,v\in[0,1].
\]

The mid-surface will be combined with a fixed airfoil profile, turned into a valid solid, joined to a hub, and replicated using the three-dimensional rotation matrix.

The nine-parameter initial vector is:

\[
\boxed{
\mathbf p=
(H,\beta,p,b_0,b_1,b_2,Z,\alpha_0,\alpha_1).
}
\]

The blade count is handled through separate optimization experiments.

**Final optimization decision:** Begin with simulated annealing, hard geometric screening, low-fidelity OpenFOAM evaluation, and high-fidelity confirmation of finalists.

Additional spline, CST, camber, and structural variables will only be introduced after the initial system is validated.

---

## 26. Reference material

- OpenFOAM User Guide: https://www.openfoam.com/documentation/user-guide
- OpenFOAM `snappyHexMesh`: https://doc.openfoam.com/2606/tools/pre-processing/mesh/generation/snappyhexmesh/
- OpenFOAM `forces` function object: https://doc.openfoam.com/2606/tools/post-processing/function-objects/forces/forces/
- OpenFOAM `kOmegaSST` turbulence model: https://doc.openfoam.com/2606/tools/processing/models/turbulence/ras/linear-evm/rtm/kOmegaSST/
- OpenFOAM wall functions: https://doc.openfoam.com/2606/tools/processing/models/turbulence/ras/wall-functions/
- OpenFOAM propeller post-processing: https://www.openfoam.com/news/main-news/openfoam-v2112/post-processing
- NASA propeller-performance coefficient reference: https://ntrs.nasa.gov/api/citations/20040070782/downloads/20040070782.pdf
- Kulfan CST geometry representation: https://www.cambridge.org/core/journals/aeronautical-journal/article/abs/recent-extensions-and-applications-of-the-cst-universal-parametric-geometry-representation-method/27BF812E1DE41B5708BA1DB6C65DCD01
