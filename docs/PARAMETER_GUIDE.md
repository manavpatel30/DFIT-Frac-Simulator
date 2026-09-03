# Parameter Guide

The dashboard starts with the following values from the `PARAMETERS` dictionary in `dfit_app.py`. Effects listed below describe the direct behavior of the current equations; coupled responses can modify the final result.

## Rock and fluid

| Dashboard key | Meaning | Default | Main effect in this model |
| --- | --- | ---: | --- |
| `E_GPa` | Young's modulus | 25.0 GPa | Increases plane-strain stiffness; generally raises net pressure and reduces elastic width for a given state. |
| `nu` | Poisson ratio | 0.25 | Changes plane-strain modulus $E'=E/(1-\nu^2)$ and therefore pressure, width, and compliance. |
| `h_f_m` | Fixed fracture height | 20.0 m | Changes fracture volume, exposed leakoff area, compliance, net pressure, and slot-flow capacity. Height does not grow dynamically. |
| `sigma_hmin_MPa` | Minimum horizontal stress | 35.0 MPa | Shifts breakdown pressure, fracture pressure, ISIP, and the aperture response relative to closure stress. |
| `tensile_strength_MPa` | Tensile strength $T_0$ | 2.0 MPa | Adds to the reported breakdown threshold. It does not otherwise modify the Stage 1 propagation equations. |
| `mu_cP` | Fluid viscosity | 1.0 cP | Affects PKN width and pressure and reduces slot-flow transmissibility as viscosity increases. |

## Stage 0: sleeve closed

| Dashboard key | Meaning | Default | Main effect in this model |
| --- | --- | ---: | --- |
| `p_initial_MPa` | Initial wellbore pressure | 28.0 MPa | Sets the starting level of the pressurization curve. |
| `q_press_m3s` | Pressurization rate | $5\times10^{-5}$ m³/s | A larger value increases the pressure-rise rate. |
| `C_wb_m3Pa` | Wellbore storage coefficient | $1.67\times10^{-9}$ m³/Pa | A larger value reduces pressure rise for the same injected rate and time. |
| `ramp_time_s` | Time to reach full pressurization rate | 30 s | A longer ramp gives less injected volume and lower pressure at a fixed sleeve-opening time. |
| `t_sleeve_open_s` | Sleeve-opening time | 360 s | Extends or shortens sealed pressurization before the breakdown check. |

## Stage 1: pumping and propagation

| Dashboard key | Meaning | Default | Main effect in this model |
| --- | --- | ---: | --- |
| `q_inj_m3s` | Total injection rate for both wings | 0.01 m³/s | Controls injected volume and affects calculated width, length, pressure, and leakoff competition. |
| `t_pump_min` | Pumping duration after sleeve opening | 30 min | A longer duration generally increases injected volume and fracture extent. |
| `C_L_m_sqrts` | Carter leakoff coefficient | $2\times10^{-6}$ m/√s | A larger value removes fluid faster, leaving less volume for fracture growth and accelerating shut-in decline. |
| `dP_total_MPa` | Total friction loss at the specified injection rate | 3.0 MPa | Raises treating pressure during pumping and sets the instantaneous treating-pressure-to-ISIP drop; it does not add to ISIP. |

## Stage 2: shut-in

| Dashboard key | Meaning | Default | Main effect in this model |
| --- | --- | ---: | --- |
| `p_res_MPa` | Reservoir pore pressure | 30.0 MPa | Sets the pressure floor/asymptote and scales the pressure-driven leakoff term. |
| `t_shutin_hr` | Shut-in observation time | 12 hr | Controls the duration of the simulated decline. |
| `Eres_mm` | Residual fracture aperture | 0.10 mm | Prevents complete numerical closure and preserves a minimum slot-flow capacity. |

## Numerical controls

| Dashboard key | Meaning | Default | Main effect in this model |
| --- | --- | ---: | --- |
| `nx` | Spatial nodes along one fracture wing | 61 | More nodes improve spatial resolution but increase nonlinear-solver cost. Values below 2 are not meaningful for the current grid calculation. |
| `dt_shutin_s` | Shut-in time step | 60 s | Smaller steps provide finer time resolution but require more nonlinear solves. |

## Additional engine defaults

`dfit_v4.py` also defines controls that are not exposed in the dashboard, including the initial fracture half-length, the fracture exposure-time exponent, the propagation time step, the residual-solver tolerances, and its maximum function-evaluation count.

## Practical checks before interpreting a run

1. Confirm that pressure at sleeve opening reaches or exceeds $\sigma_{h\min}+T_0$.
2. Keep reservoir pressure below ISIP for the intended pressure-decline interpretation.
3. Check that the selected time step and node count give a stable response when refined.
4. Treat the residual aperture as a modeling assumption, not a directly measured closure width unless it has been calibrated.
5. Compare treating pressure and ISIP separately because the friction input affects only the pumping-pressure offset.

