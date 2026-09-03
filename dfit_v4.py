"""
Simplified Sleeve-based DFIT Simulator
=======================================

Physical sequence
-----------------
  Stage 0  SLEEVE CLOSED  — Wellbore pressurises as a sealed rigid volume.
                            No fluid contacts the formation.
                            dp/dt = q(t) / C_wb
                            Ends when the operator opens the sleeve.

  Stage 1  SLEEVE OPEN / FRACTURE PROPAGATION
                            Fluid enters the formation.  If wellbore pressure
                            exceeds the breakdown pressure (σ_hmin + T₀) the
                            rock fractures and PKN propagation begins.
                            Otherwise, pressure simply equilibrates.

  Stage 2  SHUT-IN / PRESSURE DECAY
                            Injection stops.  Fracture pressure decays as fluid
                            leaks off into the formation via Carter's model.
                            Pressure asymptotes to reservoir pore pressure (p_res).

References: Valkó & Economides, Hydraulic Fracture Mechanics, 1995.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.optimize import brentq, least_squares
from dataclasses import dataclass
from typing import Optional

_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))


# ============================================================
# INPUT DATACLASSES
# ============================================================

@dataclass
class SleeveInputs:
    """
    Stage 0: sealed wellbore pressurisation (sleeve closed).

    Wellbore is a rigid compressible volume; no formation contact.
        C_wb · dp/dt = q(t)
    Rate ramps linearly from 0 to q_press over ramp_time seconds.
    Simulation ends when the sleeve is opened (t = t_sleeve_open).
    """
    sigma_hmin:       float           # Pa — minimum horizontal stress
    p_initial:        float           # Pa — initial wellbore pressure (e.g. hydrostatic)
    q_press:          float = 5e-5    # m³/s — pressurisation rate
    C_wb:             float = 1.67e-9 # m³/Pa — wellbore storage (set for realistic ramp)
    ramp_time:        float = 30.0    # s — linear rate ramp-up duration
    t_sleeve_open:    float = 300.0   # s — time at which sleeve is opened
    dt:               float = 1.0     # s


@dataclass
class FractureInputs:
    """
    Stage 1: fracture initiation and PKN propagation (sleeve open).

    As soon as the sleeve opens the wellbore pressure acts on the formation.
    If p_wb >= p_breakdown, a fracture initiates and propagates per PKN model.
    """
    E:                float           # Pa — Young's modulus
    nu:               float           # — Poisson ratio
    h_f:              float           # m — fracture height (pay zone)
    sigma_hmin:       float           # Pa — minimum horizontal stress
    tensile_strength: float           # Pa — rock tensile strength (T₀)
    mu:               float           # Pa·s — fluid viscosity
    q_inj_total:      float           # m³/s — injection rate (both wings)
    t_end:            float           # s — injection duration after sleeve opens
    dt:               float = 5.0     # s
    C_L:              float = 2e-6    # m/√s — Carter leakoff coefficient
    x0:               float = 0.5     # m — initial fracture half-length at initiation

    # ── Wellbore friction & near-wellbore pressure losses ────────────────
    # During pumping:  p_wellbore = p_fracture + friction_dp(q)
    # At shut-in:      q → 0 ⇒ friction vanishes ⇒ sharp pressure drop to ISIP
    # Set dP_total_MPa to your estimated total friction at q_inj_total.
    # (pipe + perforation + near-wellbore tortuosity combined)
    dP_total_MPa:     float = 3.0     # MPa — total friction at full injection rate

    @property
    def p_breakdown(self):
        return self.sigma_hmin + self.tensile_strength


@dataclass
class ShutInInputs:
    """
    Stage 2: shut-in pressure decay.

    Fracture geometry is fixed at end-of-pumping values.
    Pressure decays via Carter leakoff, asymptoting to p_res.
    """
    E:                float
    nu:               float
    h_f:              float
    sigma_hmin:       float
    mu:               float
    p_res:            float           # Pa — reservoir pore pressure (asymptote)
    t_end_shutin:     float           # s — observation window
    dt:               float = 60.0   # s
    nx:               int   = 61     # spatial nodes along fracture
    C_L:              float = 2e-6   # m/√s — must match Stage 1
    Eres:             float = 1e-4   # m — residual aperture (keeps fracture open)
    growth_exponent_m: float = 0.8   # PKN Newtonian exponent for τ_open profile
    lsq_ftol: float = 1e-8
    lsq_xtol: float = 1e-8
    lsq_gtol: float = 1e-8
    lsq_max_nfev: int = 10000


# ============================================================
# SHARED PHYSICS FUNCTIONS
# ============================================================

PKN_GAMMA = np.pi / 5.0          # PKN shape factor γ = w_avg/w_wellbore

def _Ep(E, nu):    return E / (1.0 - nu**2)
def _qi(q):        return 0.5 * q                      # rate per wing
def _comp(Ep, hf): return np.pi * hf / (2.0 * Ep)     # PKN compliance [m/Pa]

def pkn_avg_width(xf, mu, q, Ep):
    """Valko Eq. 9.12:  w_avg = 2.24·(µ·i·xf/E')^(1/4)"""
    return 2.24 * (mu * _qi(q) * xf / Ep) ** 0.25

def pkn_net_pressure(w_avg, Ep, hf):
    """Wellbore net pressure (Valko Eq. 9.1 + 9.7)."""
    return Ep * w_avg / (2.0 * hf * PKN_GAMMA)

def frac_vol(xf, hf, w):  return 2.0 * xf * hf * w
def exp_area(xf, hf):     return 4.0 * xf * hf

def carter_ql(xf, hf, CL, t, tmin=1.0):
    """Carter lumped leakoff rate [m³/s]."""
    return exp_area(xf, hf) * CL / np.sqrt(max(t, tmin))

def _aperture(p, sig_n, Ep, hf, Eres):
    """Single PKN linear-elastic branch (no closure mechanics)."""
    return np.maximum(Eres + _comp(Ep, hf) * (np.asarray(p, float) - sig_n), Eres)

def _slot_T(w, hf, mu, wmin=1e-7):
    """Slot-flow transmissibility with minimum-width floor."""
    return hf * max(w, wmin)**3 / (12.0 * mu)

def _avg_w(x, w):
    L = x[-1] - x[0]
    return _trapz(w, x) / L if L > 0 else float(np.mean(w))

def _tot_vol(x, w, hf):
    return 2.0 * hf * _trapz(w, x)

def _tau_profile(x, xf, tp, m):
    """Exposure-time profile τ_open(x) at end of pumping."""
    xh = np.clip(x / max(xf, 1e-12), 0.0, 1.0)
    return np.maximum(tp * (1.0 - xh ** (1.0 / max(m, 1e-12))), 0.0)

def _solve_xf(Vt, mu, q, Ep, hf):
    """Invert fracture volume equation for xf."""
    def res(xf):
        return frac_vol(xf, hf, pkn_avg_width(xf, mu, q, Ep)) - Vt
    fl, fh = res(1e-8), res(1e4)
    n = 0
    while fl * fh > 0 and n < 30:
        fh = res(1e4 * 10**n); n += 1
    return brentq(res, 1e-8, 1e4 * 10**(n-1))


# ============================================================
# STAGE 0 — SEALED WELLBORE PRESSURISATION
# ============================================================

def solve_stage0(params: SleeveInputs):
    """
    Pressurise wellbore with sleeve closed.
    No formation contact — pure wellbore storage: C_wb · dp/dt = q(t).
    """
    t_arr = np.arange(0.0, params.t_sleeve_open + params.dt, params.dt)
    p_arr = np.empty(len(t_arr))
    q_arr = np.empty(len(t_arr))

    p = params.p_initial
    for k, t in enumerate(t_arr):
        ramp = min(t / max(params.ramp_time, params.dt), 1.0)
        q = params.q_press * ramp
        if k > 0:
            p = p + q * params.dt / params.C_wb
        p_arr[k] = p
        q_arr[k] = q

    return dict(
        time=t_arr, p_wb=p_arr, q=q_arr,
        p_at_sleeve_open=float(p_arr[-1]),
        t_end=float(t_arr[-1]),
    )


# ============================================================
# FRICTION / WELLBORE PRESSURE LOSSES
# ============================================================

def friction_dp(q: float, params) -> float:
    """
    Total wellbore friction pressure drop (Pa) at injection rate q (m³/s).

    During pumping:  p_wellbore = p_fracture + friction_dp(q)
    At shut-in (q→0): friction_dp → 0, so p_wb drops instantly to p_fracture.
    This creates the sharp ISIP step seen in field DFIT pressure records.

    Scales as q^1.5 (between laminar q^1 and turbulent q^2), calibrated so
    friction_dp(q_inj_total) = dP_total_MPa.
    """
    if q <= 0.0 or params.dP_total_MPa <= 0.0:
        return 0.0
    q_ref = max(params.q_inj_total, 1e-12)
    return params.dP_total_MPa * 1e6 * (q / q_ref) ** 1.5


# ============================================================
# STAGE 1 — FRACTURE INITIATION & PKN PROPAGATION
# ============================================================

def solve_stage1(res0, params: FractureInputs):
    """
    After the sleeve opens, wellbore pressure acts on the formation.
    If p_wb ≥ p_breakdown → fracture initiates and PKN model runs.
    """
    Ep    = _Ep(params.E, params.nu)
    p_bd  = params.p_breakdown
    p_wb0 = res0["p_at_sleeve_open"]

    initiated = p_wb0 >= p_bd
    if not initiated:
        print(f"  NOTE: p_wb at sleeve open ({p_wb0/1e6:.2f} MPa) < "
              f"breakdown ({p_bd/1e6:.2f} MPa). "
              f"Increase q_press, C_wb, or t_sleeve_open.")

    times = np.arange(0.0, params.t_end + params.dt, params.dt)
    xf_h, w_h, V_h, pnet_h, pf_h, qL_h, Vi_h, Vl_h = ([] for _ in range(8))
    cum_i = cum_l = 0.0

    for k, t in enumerate(times):
        if k == 0:
            xf = params.x0
            w  = pkn_avg_width(xf, params.mu, params.q_inj_total, Ep)
            V  = frac_vol(xf, params.h_f, w)
            qL = carter_ql(xf, params.h_f, params.C_L, t)
        else:
            qL = carter_ql(xf_h[-1], params.h_f, params.C_L, t)
            Vt = max(V_h[-1] + (params.q_inj_total - qL) * params.dt, 1e-12)
            xf = _solve_xf(Vt, params.mu, params.q_inj_total, Ep, params.h_f)
            w  = pkn_avg_width(xf, params.mu, params.q_inj_total, Ep)
            V  = frac_vol(xf, params.h_f, w)
            cum_i += params.q_inj_total * params.dt
            cum_l += qL * params.dt

        pn      = pkn_net_pressure(w, Ep, params.h_f)
        p_frac  = pn + params.sigma_hmin           # true fracture pressure
        dP_fric = friction_dp(params.q_inj_total, params)
        pf      = p_frac + dP_fric                 # treating pressure (includes friction)

        xf_h.append(xf); w_h.append(w);   V_h.append(V)
        pnet_h.append(pn); pf_h.append(pf); qL_h.append(qL)
        Vi_h.append(cum_i); Vl_h.append(cum_l)

    # ISIP = fracture pressure at shut-in (friction = 0 at q=0)
    isip = float(pnet_h[-1]) + params.sigma_hmin

    return dict(
        time=np.array(times), xf=np.array(xf_h), w_avg=np.array(w_h),
        Vf=np.array(V_h), p_net=np.array(pnet_h), p_f=np.array(pf_h),
        q_leakoff=np.array(qL_h), Vinj_cum=np.array(Vi_h),
        Vleak_cum=np.array(Vl_h), Eprime=Ep,
        p_breakdown=p_bd, initiated=initiated,
        isip=isip, dP_friction=dP_fric,
    )


# ============================================================
# STAGE 2 — SHUT-IN PRESSURE DECAY
# ============================================================

def _shutin_residual(p_new, p_old, dx, hf, mu, sig_n, Ep, Eres,
                     CL, tau, t_rel, dt, p_res, p_si):
    """
    Implicit FD residual (pressure units [Pa] via dt/comp scaling).

    Leakoff: v_L = [C_L/√(τ+Δt)] · (p−p_res)/(p_si−p_res)
    Sealed ends at both x=0 and x=xf.
    """
    comp = _comp(Ep, hf)
    ker  = CL / np.sqrt(np.maximum(tau + t_rel, 1.0))
    drv  = np.maximum(p_new - p_res, 0.0) / max(p_si - p_res, 1.0)
    vL   = ker * drv

    w    = _aperture(p_new, sig_n, Ep, hf, Eres)
    n    = len(p_new)
    qf   = np.zeros(n + 1)
    for i in range(n - 1):
        wf     = 0.5 * (w[i] + w[i + 1])
        T      = _slot_T(wf, hf, mu)
        qf[i+1] = -T * (p_new[i+1] - p_new[i]) / dx
    qf[0] = qf[n] = 0.0  # sealed at both ends

    R = np.empty(n)
    for i in range(n):
        R[i] = (comp * (p_new[i] - p_old[i]) / dt
                + (qf[i+1] - qf[i]) / (dx * hf)
                + 2.0 * vL[i])
    # Scale to pressure units so Jacobian ≈ O(1)
    return R * (dt / comp)


def solve_stage2(res1, params: ShutInInputs, t_abs_offset: float = 0.0):
    """
    Post-shut-in pressure decay with p → p_res as t → ∞.
    """
    Ep         = _Ep(params.E, params.nu)
    xf_fin     = float(res1["xf"][-1])
    p_si       = float(res1["isip"])  # ISIP = fracture pressure at shut-in (no friction)
    # t_pump_end is the ABSOLUTE time of shut-in.
    # res1["time"] is LOCAL to Stage 1 (starts at 0), so add t_abs_offset
    # (= Stage 0 duration) to get the correct absolute timestamp.
    t_pump_end = float(res1["time"][-1]) + t_abs_offset

    nx   = params.nx
    x    = np.linspace(0.0, xf_fin, nx)
    dx   = (x[1] - x[0]) if nx > 1 else xf_fin
    tau  = _tau_profile(x, xf_fin, t_pump_end, params.growth_exponent_m)

    times = np.arange(0.0, params.t_end_shutin + params.dt, params.dt)
    lower = np.full(nx, params.p_res)   # hard floor = p_res → correct asymptote

    p_old = np.full(nx, p_si)
    w0    = _aperture(p_old, params.sigma_hmin, Ep, params.h_f, params.Eres)

    p_hist = [p_old.copy()]; w_hist = [w0.copy()]
    ph_h   = [float(p_old[0])]; pa_h = [float(np.mean(p_old))]
    wa_h   = [_avg_w(x, w0)];   V_h  = [_tot_vol(x, w0, params.h_f)]
    ql_h   = [0.0]; tabs_h = [t_pump_end]; cl_h = [0.0]

    for step in range(1, len(times)):
        t_rel  = times[step]
        p_snap = p_old.copy()

        def res(pn, _po=p_snap, _tr=t_rel):
            return _shutin_residual(pn, _po, dx, params.h_f, params.mu,
                                    params.sigma_hmin, Ep, params.Eres,
                                    params.C_L, tau, _tr, params.dt,
                                    params.p_res, p_si)

        sol = least_squares(
            res, np.maximum(p_old.copy(), lower + 1.0),
            bounds=(lower, np.full(nx, np.inf)),
            method="trf",
            ftol=params.lsq_ftol, xtol=params.lsq_xtol,
            gtol=params.lsq_gtol, max_nfev=params.lsq_max_nfev,
            x_scale="jac")

        p_new = sol.x.copy()
        w_new = _aperture(p_new, params.sigma_hmin, Ep, params.h_f, params.Eres)
        ker   = params.C_L / np.sqrt(np.maximum(tau + t_rel, 1.0))
        drv   = np.maximum(p_new - params.p_res, 0.0) / max(p_si - params.p_res, 1.0)
        vL    = ker * drv
        qL    = 4.0 * params.h_f * _trapz(vL, x)

        p_hist.append(p_new.copy()); w_hist.append(w_new.copy())
        ph_h.append(float(p_new[0])); pa_h.append(float(np.mean(p_new)))
        wa_h.append(_avg_w(x, w_new)); V_h.append(_tot_vol(x, w_new, params.h_f))
        ql_h.append(qL); tabs_h.append(t_pump_end + t_rel)
        cl_h.append(cl_h[-1] + qL * params.dt)
        p_old = p_new

    return dict(
        x=x, tau_open=tau,
        time_since_shutin=times, time_absolute=np.array(tabs_h),
        p_profile=np.array(p_hist), w_profile=np.array(w_hist),
        p_avg=np.array(pa_h), p_heel=np.array(ph_h), p_wb=np.array(ph_h),
        w_avg=np.array(wa_h), Vf=np.array(V_h),
        q_leakoff=np.array(ql_h), Vleak_cum=np.array(cl_h),
        xf_fixed=xf_fin, p_si=p_si, p_res=params.p_res, Eprime=Ep,
    )


# ============================================================
# PLOTTING  — focused on fracture geometry
# ============================================================

def plot_results(res0, res1, res2,
                 p0: SleeveInputs, p1: FractureInputs, p2: ShutInInputs):

    plt.rcParams.update({"font.size": 11, "lines.linewidth": 2.0,
                         "axes.spines.top": False, "axes.spines.right": False})

    clmpa  = p1.sigma_hmin / 1e6
    p_bd   = p1.p_breakdown / 1e6
    p_res  = p2.p_res / 1e6
    t0_off = 0.0
    t1_off = res0["t_end"]
    t2_off = t1_off + float(res1["time"][-1])

    # Assemble full wellbore-pressure time series
    t_s0   = res0["time"]
    t_s1   = t1_off + res1["time"]
    t_s2   = np.array(res2["time_absolute"])

    p_s0   = res0["p_wb"] / 1e6
    p_s1   = res1["p_f"]  / 1e6
    p_s2   = res2["p_wb"] / 1e6

    q_s0   = res0["q"] * 60.0                                         # m³/min
    q_s1   = np.full(len(t_s1), p1.q_inj_total * 60.0)
    q_s2   = np.zeros(len(t_s2))

    tall   = np.concatenate([t_s0, t_s1, t_s2]) / 60.0               # minutes
    pall   = np.concatenate([p_s0, p_s1, p_s2])
    qall   = np.concatenate([q_s0, q_s1, q_s2])

    col_s0 = "#2196F3"   # blue  — pressurisation
    col_s1 = "#E53935"   # red   — propagation
    col_s2 = "#43A047"   # green — decay

    # ── Figure 1: Full operational history ──────────────────────────────
    fig1, (ax_p, ax_q) = plt.subplots(2, 1, figsize=(12, 7),
                                       gridspec_kw={"height_ratios": [3, 1]},
                                       sharex=True)

    # Colour the background per stage
    t_s1_min = t1_off / 60.0
    t_s2_min = t2_off / 60.0
    t_end_min = tall[-1]
    ax_p.axvspan(0,          t_s1_min, alpha=0.07, color=col_s0, zorder=0)
    ax_p.axvspan(t_s1_min,  t_s2_min, alpha=0.07, color=col_s1, zorder=0)
    ax_p.axvspan(t_s2_min, t_end_min, alpha=0.07, color=col_s2, zorder=0)

    ax_p.plot(tall, pall, color="navy", lw=2.2, zorder=3)
    isip_mpa = res1["isip"] / 1e6
    ax_p.axhline(clmpa,    color="crimson",    ls="--", lw=1.4, label=f"Closure stress  {clmpa:.1f} MPa")
    ax_p.axhline(p_bd,     color="darkorange", ls=":",  lw=1.4, label=f"Breakdown  {p_bd:.1f} MPa")
    ax_p.axhline(isip_mpa, color="gold",       ls="--", lw=1.2, label=f"ISIP  {isip_mpa:.1f} MPa")
    ax_p.axhline(p_res,    color="purple",     ls="-.", lw=1.4, label=f"Reservoir pressure  {p_res:.1f} MPa")

    # Stage boundary markers
    for t_v, lbl in [(t_s1_min, "Sleeve opens\n(frac. initiation)"),
                     (t_s2_min, "Shut-in")]:
        ax_p.axvline(t_v, color="gray", ls="--", lw=1.0, alpha=0.8)
        ax_p.text(t_v + 0.5, ax_p.get_ylim()[0] if ax_p.get_ylim()[0] else clmpa * 0.97,
                  lbl, fontsize=8.5, color="gray", va="bottom")

    # Annotate the friction pressure drop at shut-in
    treat_mpa = res1["p_f"][-1] / 1e6
    isip_mpa2 = res1["isip"] / 1e6
    if t_s2_min >= ax_p.get_xlim()[0] and t_s2_min <= ax_p.get_xlim()[1]:
        ax_p.annotate(
            "", xy=(t_s2_min + 0.5, isip_mpa2),
            xytext=(t_s2_min + 0.5, treat_mpa),
            arrowprops=dict(arrowstyle="<->", color="gold", lw=1.5))
        dp_label = f"dP_fric\n{(treat_mpa-isip_mpa2):.1f} MPa"
        ax_p.text(t_s2_min + 0.6,
                  (treat_mpa + isip_mpa2) / 2,
                  dp_label,
                  fontsize=8, color="gold", va="center")

    # Stage labels
    mid0 = t_s1_min / 2
    mid1 = (t_s1_min + t_s2_min) / 2
    mid2 = (t_s2_min + t_end_min) / 2
    for mid, lbl, col in [(mid0, "Pressurisation\n(sleeve closed)", col_s0),
                           (mid1, "Fracture\nPropagation",            col_s1),
                           (mid2, "Pressure\nDecay",                  col_s2)]:
        ax_p.text(mid, ax_p.get_ylim()[1] if ax_p.get_ylim()[1] else clmpa*1.05,
                  lbl, ha="center", fontsize=9, color=col,
                  fontweight="bold", va="top")

    ax_p.set_ylabel("Wellbore Pressure (MPa)")
    ax_p.set_title("Sleeve-Based DFIT: Operational History", fontweight="bold")
    ax_p.legend(loc="lower right", fontsize=9)
    ax_p.grid(True, alpha=0.25)

    ax_q.fill_between(tall, qall, color=col_s1, alpha=0.5, step="pre")
    ax_q.axvspan(0, t_s1_min, alpha=0.07, color=col_s0)
    ax_q.axvspan(t_s1_min, t_s2_min, alpha=0.07, color=col_s1)
    ax_q.set_ylabel("Rate\n(m³/min)")
    ax_q.set_xlabel("Absolute Time (min)")
    ax_q.grid(True, alpha=0.25)
    fig1.tight_layout()

    # ── Figure 2: Fracture geometry evolution during pumping ────────────
    fig2, axes = plt.subplots(1, 3, figsize=(15, 5))
    t1_min = res1["time"] / 60.0

    axes[0].plot(t1_min, res1["xf"], color=col_s1, lw=2.2)
    axes[0].set_xlabel("Pumping Time (min)")
    axes[0].set_ylabel("Half-Length  xf (m)")
    axes[0].set_title("Fracture Half-Length Growth")
    axes[0].fill_between(t1_min, res1["xf"], alpha=0.15, color=col_s1)
    axes[0].grid(True, alpha=0.25)

    axes[1].plot(t1_min, res1["w_avg"] * 1e3, color="#8E24AA", lw=2.2)
    axes[1].set_xlabel("Pumping Time (min)")
    axes[1].set_ylabel("Average Width  w̄ (mm)")
    axes[1].set_title("Average Fracture Width Growth")
    axes[1].fill_between(t1_min, res1["w_avg"] * 1e3, alpha=0.15, color="#8E24AA")
    axes[1].grid(True, alpha=0.25)

    # Fracture height is constant in PKN model — show alongside width for completeness
    h_arr = np.full(len(t1_min), p1.h_f)
    axes[2].plot(t1_min, h_arr, color="#00897B", lw=2.2, ls="--",
                 label=f"Height  h_f = {p1.h_f:.0f} m  (PKN fixed)")
    axes[2].plot(t1_min, res1["xf"], color=col_s1, lw=2.2, label="Half-length  xf (m)")
    axes[2].plot(t1_min, res1["w_avg"] * 1e3 * 10, color="#8E24AA", lw=2.0, ls="-.",
                 label="10 × Avg width (mm)")   # scaled for visibility
    axes[2].set_xlabel("Pumping Time (min)")
    axes[2].set_ylabel("Dimension")
    axes[2].set_title("All Fracture Dimensions\n(width ×10 for visibility)")
    axes[2].legend(fontsize=9)
    axes[2].grid(True, alpha=0.25)
    fig2.suptitle("Fracture Geometry Evolution During Pumping", fontweight="bold")
    fig2.tight_layout()

    # ── Figure 3: Fracture footprint cross-section snapshots ────────────
    # Show the fracture as a 2-D bird's-eye elliptical footprint at
    # 3 pumping times: early, mid, end-of-pumping.
    fig3, axes3 = plt.subplots(1, 3, figsize=(15, 5))
    n1 = len(res1["time"])
    snap_idx = [max(n1//10, 1), n1//2, n1-1]
    snap_lbls = ["Early", "Mid", "End of pumping"]
    cmap_frac  = plt.cm.plasma

    for col_ax, si, lbl in zip(axes3, snap_idx, snap_lbls):
        xf_s = res1["xf"][si]
        w_s  = res1["w_avg"][si]
        t_s  = res1["time"][si] / 60.0
        hf_s = p1.h_f

        # Bird's eye: ellipse with semi-axes xf (lateral) and hf/2 (vertical)
        theta = np.linspace(0, 2*np.pi, 300)
        x_el  = xf_s * np.cos(theta)
        y_el  = (hf_s / 2) * np.sin(theta)

        # Colour fill from blue (tip) to red (wellbore) by width profile
        # Width is roughly proportional to (1 - (x/xf)^(4/3))^(1/4) for PKN
        xi    = np.linspace(-xf_s, xf_s, 200)
        xi_hat = np.abs(xi) / max(xf_s, 1e-9)
        w_xi  = w_s * (1.0 - xi_hat**(4.0/3.0))**0.25  # approximate PKN profile
        w_xi  = np.maximum(w_xi, 0)

        # Plot width as colour-coded horizontal slices
        for j in range(len(xi) - 1):
            frac = w_xi[j] / max(w_s, 1e-9)
            yi_max = (hf_s/2) * np.sqrt(max(1 - (xi[j]/xf_s)**2, 0))
            col_ax.fill_betweenx([-yi_max, yi_max],
                                  xi[j], xi[j+1],
                                  color=cmap_frac(0.2 + 0.7*frac), alpha=0.8)

        col_ax.plot(x_el, y_el, "k-", lw=1.5)
        col_ax.axvline(0, color="gray", ls="--", lw=1, alpha=0.6)
        col_ax.set_aspect("equal")
        col_ax.set_xlabel("Lateral extent (m)")
        col_ax.set_ylabel("Height (m)")
        col_ax.set_title(f"{lbl}\nt = {t_s:.1f} min\n"
                         f"xf = {xf_s:.0f} m   w̄ = {w_s*1e3:.2f} mm")
        col_ax.grid(True, alpha=0.2)
        col_ax.text(0, hf_s/2 * 1.08, "Wellbore", ha="center", fontsize=8, color="gray")

    sm = plt.cm.ScalarMappable(cmap=cmap_frac, norm=plt.Normalize(0, 1))
    sm.set_array([])
    cbar = fig3.colorbar(sm, ax=axes3, orientation="vertical", shrink=0.7, pad=0.02)
    cbar.set_label("Relative Width\n(0=tip, 1=max)", fontsize=9)
    fig3.suptitle("Fracture Footprint Cross-Section  (bird's-eye view)",
                  fontweight="bold")
    fig3.tight_layout()

    # ── Figure 4: Material balance ───────────────────────────────────────
    fig4, (ax_mb, ax_ql) = plt.subplots(1, 2, figsize=(13, 5))

    ax_mb.stackplot(t1_min,
                    [res1["Vf"], res1["Vleak_cum"]],
                    labels=["Fracture volume", "Cumulative leakoff"],
                    colors=[col_s1, "#90A4AE"], alpha=0.8)
    ax_mb.plot(t1_min, res1["Vinj_cum"], "k--", lw=2, label="Injected volume")
    ax_mb.set_xlabel("Pumping Time (min)")
    ax_mb.set_ylabel("Volume (m³)")
    ax_mb.set_title("Material Balance During Pumping")
    ax_mb.legend(fontsize=9)
    ax_mb.grid(True, alpha=0.25)

    eff = res1["Vf"] / np.maximum(res1["Vinj_cum"], 1e-12) * 100
    ax_ql.plot(t1_min, eff, color=col_s1, lw=2.2)
    ax_ql.set_xlabel("Pumping Time (min)")
    ax_ql.set_ylabel("Fluid Efficiency (%)")
    ax_ql.set_title("Fluid Efficiency  =  Vf / V_injected")
    ax_ql.set_ylim(0, 105)
    ax_ql.grid(True, alpha=0.25)
    ax_ql.fill_between(t1_min, eff, alpha=0.15, color=col_s1)
    fig4.tight_layout()

    # ── Figure 5: Shut-in pressure and width spatial decay ──────────────
    fig5, (ax5l, ax5r) = plt.subplots(1, 2, figsize=(13, 5))
    ns   = len(res2["time_since_shutin"])
    snap = [0, ns//4, ns//2, ns*3//4, ns-1]
    cols5 = ["#1565C0", "#1976D2", "#42A5F5", "#90CAF9", "#E3F2FD"]

    for si, c in zip(snap, cols5):
        tm = res2["time_since_shutin"][si] / 60.0
        ax5l.plot(res2["x"], res2["p_profile"][si]/1e6, color=c, lw=1.8,
                  label=f"t = {tm:.0f} min")
        ax5r.plot(res2["x"], res2["w_profile"][si]*1e3, color=c, lw=1.8,
                  label=f"t = {tm:.0f} min")

    ax5l.axhline(clmpa, color="crimson", ls="--", lw=1.2, label="Closure stress")
    ax5l.axhline(p_res,  color="purple",  ls="-.", lw=1.2, label="Reservoir pressure")
    ax5l.set_xlabel("Distance from Wellbore (m)")
    ax5l.set_ylabel("Pressure (MPa)")
    ax5l.set_title("Pressure Profile Along Fracture\n(Shut-in decay)")
    ax5l.legend(fontsize=8.5, loc="upper right")
    ax5l.grid(True, alpha=0.25)

    ax5r.set_xlabel("Distance from Wellbore (m)")
    ax5r.set_ylabel("Aperture Width (mm)")
    ax5r.set_title("Aperture Profile Along Fracture\n(Shut-in decay)")
    ax5r.legend(fontsize=8.5, loc="upper right")
    ax5r.grid(True, alpha=0.25)
    fig5.tight_layout()

    # ── Figure 6: Wellbore pressure decay + fracture dimensions vs time ─
    fig6, axes6 = plt.subplots(1, 3, figsize=(15, 5))
    t2_min = res2["time_since_shutin"] / 60.0

    axes6[0].plot(t2_min, res2["p_wb"]/1e6, color=col_s2, lw=2.2)
    axes6[0].axhline(clmpa, color="crimson",    ls="--", lw=1.2, label="Closure stress")
    axes6[0].axhline(p_res,  color="purple",     ls="-.", lw=1.2, label="p_res (asymptote)")
    axes6[0].set_xlabel("Time Since Shut-in (min)")
    axes6[0].set_ylabel("Wellbore Pressure (MPa)")
    axes6[0].set_title("Pressure Decay\n(asymptotes to p_res)")
    axes6[0].legend(fontsize=9)
    axes6[0].grid(True, alpha=0.25)

    axes6[1].plot(t2_min, res2["w_avg"]*1e3, color="#8E24AA", lw=2.2)
    axes6[1].set_xlabel("Time Since Shut-in (min)")
    axes6[1].set_ylabel("Spatially-Averaged Width (mm)")
    axes6[1].set_title("Average Fracture Width Decay")
    axes6[1].fill_between(t2_min, res2["w_avg"]*1e3, alpha=0.15, color="#8E24AA")
    axes6[1].grid(True, alpha=0.25)

    axes6[2].plot(t2_min, np.array(res2["q_leakoff"])*1e3, color="#F57C00", lw=2.2)
    axes6[2].set_xlabel("Time Since Shut-in (min)")
    axes6[2].set_ylabel("Total Leakoff Rate (×10⁻³ m³/s)")
    axes6[2].set_title("Leakoff Rate Decay\n(∝ 1/√t · pressure drive)")
    axes6[2].fill_between(t2_min, np.array(res2["q_leakoff"])*1e3,
                          alpha=0.15, color="#F57C00")
    axes6[2].grid(True, alpha=0.25)
    fig6.suptitle("Shut-in Monitoring", fontweight="bold")
    fig6.tight_layout()

    plt.show()


# ============================================================
# EXAMPLE RUN
# ============================================================

if __name__ == "__main__":
    hours = 3600.0

    # ── Stage 0: Sealed wellbore pressurisation ─────────────────────────
    # Sleeves are closed — no formation contact.
    # Inject at low rate; pressure rises until sleeve is opened.
    p0 = SleeveInputs(
        sigma_hmin    = 35e6,       # Pa  (5075 psi)
        p_initial     = 28e6,       # Pa  — hydrostatic wellbore pressure at start
        q_press       = 5e-5,       # m³/s ≈ 0.019 bbl/min (slow pressurisation)
        C_wb          = 1.67e-9,    # m³/Pa — wellbore storage
        ramp_time     = 30.0,       # s  — rate ramp duration
        t_sleeve_open = 360.0,      # s  — open sleeve after 6 min (reaches breakdown)
        dt            = 1.0,        # s
    )
    res0 = solve_stage0(p0)
    print("Stage 0 — Pressurisation:")
    print(f"  Duration              : {res0['t_end']:.0f} s  ({res0['t_end']/60:.1f} min)")
    print(f"  Pressure at sleeve open: {res0['p_at_sleeve_open']/1e6:.2f} MPa")

    # ── Stage 1: Fracture propagation ───────────────────────────────────
    # Sleeve opens → if p_wb ≥ σ_hmin + T₀, fracture initiates and propagates.
    p1 = FractureInputs(
        E                = 25e9,       # Pa
        nu               = 0.25,
        h_f              = 20.0,       # m  — fixed fracture height (PKN)
        sigma_hmin       = 35e6,       # Pa
        tensile_strength = 2.0e6,      # Pa → breakdown at 37 MPa
        mu               = 0.001,      # Pa·s (water)
        q_inj_total      = 0.01,       # m³/s ≈ 3.8 bbl/min (both wings)
        t_end            = 30 * 60.0,  # s  — 30 min of pumping
        dt               = 5.0,        # s
        C_L              = 2e-6,       # m/√s
        dP_total_MPa     = 3.0,        # MPa — total friction at q_inj_total
                                       #   (pipe + perf + near-wellbore tortuosity)
                                       #   Creates sharp ISIP drop at shut-in
    )
    res1 = solve_stage1(res0, p1)
    print("\nStage 1 — Fracture Propagation:")
    print(f"  Breakdown pressure : {res1['p_breakdown']/1e6:.2f} MPa")
    print(f"  Fracture initiated : {res1['initiated']}")
    print(f"  xf (half-length)   : {res1['xf'][-1]:.1f} m")
    print(f"  h_f (height)       : {p1.h_f:.1f} m  (PKN model — fixed)")
    print(f"  w_avg (width)      : {res1['w_avg'][-1]*1e3:.3f} mm")
    print(f"  Treating pressure  : {res1['p_f'][-1]/1e6:.3f} MPa  (fracture + friction)")
    print(f"  Friction ΔP        : {res1['dP_friction']/1e6:.3f} MPa  (drops to 0 at shut-in)")
    print(f"  ISIP               : {res1['isip']/1e6:.3f} MPa  (treating − friction)")
    print(f"  Fluid efficiency   : {res1['Vf'][-1]/max(res1['Vinj_cum'][-1],1e-9)*100:.1f}%")

    # ── Stage 2: Shut-in pressure decay ─────────────────────────────────
    # Pressure asymptotes to p_res (pore pressure).
    p2 = ShutInInputs(
        E            = 25e9,
        nu           = 0.25,
        h_f          = 20.0,
        sigma_hmin   = 35e6,
        mu           = 0.001,
        p_res        = 30.0e6,       # Pa — reservoir pore pressure
        t_end_shutin = 24 * hours,   # s  — 24 h observation
        dt           = 60.0,         # s
        nx           = 61,
        C_L          = 2e-6,         # m/√s — must match Stage 1
        Eres         = 1e-4,         # m   — 0.1 mm residual aperture
    )
    res2 = solve_stage2(res1, p2, t_abs_offset=res0["t_end"])
    pwb = res2["p_wb"]
    print("\nStage 2 — Shut-in Decay:")
    print(f"  ISIP               : {pwb[0]/1e6:.3f} MPa")
    print(f"  p at  6 h          : {pwb[min(360,len(pwb)-1)]/1e6:.3f} MPa")
    print(f"  p at 24 h          : {pwb[-1]/1e6:.3f} MPa")
    print(f"  Reservoir pressure : {p2.p_res/1e6:.1f} MPa  ← asymptote")
    print(f"  Monotone decline   : {bool(np.all(np.diff(pwb) <= 1e3))}")

    plot_results(res0, res1, res2, p0, p1, p2)
