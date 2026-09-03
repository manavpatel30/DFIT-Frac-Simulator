# Model Overview

## Purpose

The simulator exposes the physical sequence of a simplified sleeve-operated DFIT rather than treating the calculation as a black box. It separates the event into three stages and carries the final state of each stage into the next.

All engine calculations use SI units. The dashboard converts selected inputs and outputs to practical units such as GPa, MPa, cP, minutes, hours, and millimeters.

## Stage 0: sealed-wellbore pressurization

Before the sleeve opens, the wellbore is treated as a rigid compressible volume with no communication to the formation. The rate ramps linearly until it reaches the specified pressurization rate:

$$
q(t)=q_{\mathrm{press}}\min\left(\frac{t}{t_{\mathrm{ramp}}},1\right).
$$

Pressure is advanced from the wellbore-storage balance:

$$
C_{wb}\frac{dp}{dt}=q(t),
$$

or, in the time-stepping form used by the code,

$$
p_{k}=p_{k-1}+\frac{q_k\Delta t}{C_{wb}}.
$$

The stage ends at `t_sleeve_open`. The resulting pressure becomes the pressure available when the formation is exposed.

## Breakdown check

The reported breakdown threshold is

$$
p_{bd}=\sigma_{h\min}+T_0,
$$

where $\sigma_{h\min}$ is minimum horizontal stress and $T_0$ is tensile strength. The code reports fracture initiation when the pressure at sleeve opening is at least this threshold.

Important implementation detail: this comparison sets the `initiated` flag and may print a note, but the current source still calculates the Stage 1 PKN arrays when the criterion is not met. A non-initiated result therefore should not be read as a physically valid propagated fracture case.

## Stage 1: fracture propagation during pumping

The fracture is treated as a symmetric two-wing, fixed-height PKN fracture. Plane-strain modulus is

$$
E'=\frac{E}{1-\nu^2},
$$

and the injection rate supplied to each wing is one-half of the total rate.

The average width relation is

$$
\bar{w}=2.24\left(\frac{\mu q_i x_f}{E'}\right)^{1/4},
$$

with $q_i=q/2$. The wellbore net pressure is

$$
p_{net}=\frac{E'\bar{w}}{2h_f\gamma},
\qquad \gamma=\frac{\pi}{5}.
$$

Fracture volume and exposed area are represented by

$$
V_f=2x_fh_f\bar{w},
\qquad
A=4x_fh_f.
$$

The lumped Carter leakoff rate is

$$
q_L=\frac{A C_L}{\sqrt{t}},
$$

with a one-second minimum applied to the time in the numerical implementation. At each step, the code updates target fracture volume using injected volume minus leakoff, then numerically inverts the volume relationship to obtain fracture half-length.

### Treating pressure, friction, and ISIP

The calibrated friction-pressure loss is

$$
\Delta p_{fric}=\Delta p_{total}
\left(\frac{q}{q_{ref}}\right)^{1.5}.
$$

During constant-rate pumping, $q_{ref}$ equals the specified total injection rate, so the loss equals `dP_total_MPa`. The reported treating pressure is

$$
p_{treat}=\sigma_{h\min}+p_{net}+\Delta p_{fric}.
$$

At shut-in, rate and friction fall to zero. The instantaneous shut-in pressure is therefore

$$
ISIP=\sigma_{h\min}+p_{net}.
$$

## Stage 2: post-shut-in decline

Fracture half-length is fixed at its end-of-pumping value and discretized along one wing. Both ends of the one-dimensional domain are sealed for internal slot flow.

The aperture-pressure relationship is

$$
w=\max\left[E_{res}+\frac{\pi h_f}{2E'}
\left(p-\sigma_{h\min}\right),\ E_{res}\right].
$$

This is a linear-elastic branch with a prescribed residual aperture. It does not represent progressive asperity contact or a nonlinear fracture-compliance law.

Leakoff is modified by local exposure time and pressure drive:

$$
v_L=\frac{C_L}{\sqrt{\tau+\Delta t}}
\frac{p-p_{res}}{p_{SI}-p_{res}}.
$$

The pressure is bounded below by reservoir pressure. Internal fracture flow uses a slot-flow transmissibility proportional to

$$
T\propto\frac{h_fw^3}{12\mu}.
$$

The nonlinear implicit finite-difference residual is solved at each time step with SciPy's bounded least-squares solver.

## Primary assumptions and limitations

- Single, planar, symmetric two-wing fracture.
- PKN geometry with constant prescribed height.
- Homogeneous, isotropic elastic properties.
- Constant pumping rate during Stage 1.
- Lumped Carter leakoff during propagation.
- Pressure-scaled Carter-type leakoff after shut-in.
- Empirical rate exponent for combined wellbore, perforation, and near-wellbore friction.
- Fixed fracture length after shut-in.
- Linear-elastic aperture with a residual-width floor.
- No explicit wellbore storage after shut-in.
- No layered stress or rock-property contrasts.
- No height growth, natural-fracture interaction, branching, or fracture-network behavior.
- No proppant transport, thermal coupling, poroelastic stress change, non-Darcy reservoir flow, or multiphase effects.
- No automated calibration, uncertainty quantification, or field-data history matching.
- User-entered dashboard values are converted to numbers, but physical ranges are not comprehensively validated.

These assumptions make the project useful for learning parameter-response relationships and testing a transparent conceptual workflow. Quantitative field interpretation requires calibration, validation, and additional physics appropriate to the formation and test design.
