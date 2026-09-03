"""
DFIT Live Dashboard  — dfit_app.py
===================================
Run:  python dfit_app.py
Then open:  http://localhost:5050  in your browser.

Edit any parameter in the PARAMETERS section below and click
"Re-run Simulation" in the browser — the charts update immediately.

No external libraries required (only Python stdlib + numpy + scipy).
"""

# ======================================================================
# ██████╗  █████╗ ██████╗  █████╗ ███╗   ███╗███████╗████████╗███████╗
# ██╔══██╗██╔══██╗██╔══██╗██╔══██╗████╗ ████║██╔════╝╚══██╔══╝██╔════╝
# ██████╔╝███████║██████╔╝███████║██╔████╔██║█████╗     ██║   █████╗
# ██╔═══╝ ██╔══██║██╔══██╗██╔══██║██║╚██╔╝██║██╔══╝     ██║   ██╔══╝
# ██║     ██║  ██║██║  ██║██║  ██║██║ ╚═╝ ██║███████╗   ██║   ███████╗
# ╚═╝     ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝╚══════╝   ╚═╝   ╚══════╝
# ======================================================================
# Edit anything in this section and click Re-run in the browser.

PARAMETERS = {

    # ── Rock & fluid ────────────────────────────────────────────────────
    "E_GPa":               25.0,    # Young's modulus (GPa)
    "nu":                  0.25,    # Poisson ratio
    "h_f_m":               20.0,    # Fracture height — PKN fixed (m)
    "sigma_hmin_MPa":      35.0,    # Minimum horizontal stress / closure (MPa)
    "tensile_strength_MPa": 2.0,    # Rock tensile strength T₀ (MPa)
    "mu_cP":               1.0,     # Fluid viscosity (cP)

    # ── Stage 0: Sealed wellbore pressurisation ──────────────────────────
    "p_initial_MPa":       28.0,    # Initial wellbore pressure (MPa)
    "q_press_m3s":         5e-5,    # Pressurisation rate (m³/s)
    "C_wb_m3Pa":           1.67e-9, # Wellbore storage coefficient (m³/Pa)
    "ramp_time_s":         30.0,    # Rate ramp-up duration (s)
    "t_sleeve_open_s":     360.0,   # When sleeve is opened (s)

    # ── Stage 1: Fracture propagation ────────────────────────────────────
    "q_inj_m3s":           0.01,    # Injection rate — both wings (m³/s)
    "t_pump_min":          30.0,    # Pumping duration (min)
    "C_L_m_sqrts":         2e-6,    # Carter leakoff coefficient (m/√s)
    "dP_total_MPa":        3.0,     # Total friction pressure at q_inj (MPa)

    # ── Stage 2: Shut-in pressure decay ──────────────────────────────────
    "p_res_MPa":           30.0,    # Reservoir pore pressure — asymptote (MPa)
    "t_shutin_hr":         12.0,    # Observation window (hours)
    "Eres_mm":             0.10,    # Residual fracture aperture (mm)

    # ── Solver ───────────────────────────────────────────────────────────
    "nx":                  61,      # Spatial nodes along fracture
    "dt_shutin_s":         60.0,    # Shut-in time step (s)
}

# ======================================================================
# SIMULATION ENGINE  (imported from dfit_v3.py)
# ======================================================================
import sys, io, json, traceback, threading
import numpy as np
from http.server import HTTPServer, BaseHTTPRequestHandler
from urllib.parse import urlparse, parse_qs

# Load simulation functions from dfit_v3.py
_src = open('dfit_v4.py').read().split('if __name__')[0]
exec(_src, globals())


def run_simulation(params):
    """Run all three stages and return JSON-serialisable data dict."""

    E          = params["E_GPa"]       * 1e9
    nu         = params["nu"]
    h_f        = params["h_f_m"]
    sig_hmin   = params["sigma_hmin_MPa"] * 1e6
    T0         = params["tensile_strength_MPa"] * 1e6
    mu         = params["mu_cP"]       * 1e-3
    p_res      = params["p_res_MPa"]   * 1e6
    C_L        = params["C_L_m_sqrts"]
    Eres       = params["Eres_mm"]     * 1e-3
    hours      = 3600.0

    # Stage 0
    p0 = SleeveInputs(
        sigma_hmin       = sig_hmin,
        p_initial        = params["p_initial_MPa"] * 1e6,
        q_press          = params["q_press_m3s"],
        C_wb             = params["C_wb_m3Pa"],
        ramp_time        = params["ramp_time_s"],
        t_sleeve_open    = params["t_sleeve_open_s"],
        dt               = 1.0,
    )
    res0 = solve_stage0(p0)

    # Stage 1
    p1 = FractureInputs(
        E=E, nu=nu, h_f=h_f, sigma_hmin=sig_hmin,
        tensile_strength = T0,
        mu               = mu,
        q_inj_total      = params["q_inj_m3s"],
        t_end            = params["t_pump_min"] * 60.0,
        dt               = 5.0,
        C_L              = C_L,
        dP_total_MPa     = params.get("dP_total_MPa", 3.0),
    )
    res1 = solve_stage1(res0, p1)

    # Stage 2
    p2 = ShutInInputs(
        E=E, nu=nu, h_f=h_f, sigma_hmin=sig_hmin, mu=mu,
        p_res        = p_res,
        t_end_shutin = params["t_shutin_hr"] * hours,
        dt           = params["dt_shutin_s"],
        nx           = int(params["nx"]),
        C_L          = C_L,
        Eres         = Eres,
    )
    res2 = solve_stage2(res1, p2, t_abs_offset=res0["t_end"])

    # ── Build animation frames ──────────────────────────────────────────
    t0_off = 0.0
    t1_off = res0["t_end"]
    t2_off = t1_off + float(res1["time"][-1])  # absolute time of shut-in

    def sub(n_total, n_want):
        return [int(i) for i in np.linspace(0, n_total-1, min(n_want, n_total))]

    frames = []
    for i in sub(len(res1["time"]), 50):
        xf_i = float(res1["xf"][i])
        w_i  = float(res1["w_avg"][i])
        xi   = np.linspace(0, max(xf_i, 1e-3), 20)
        xi_h = xi / max(xf_i, 1e-9)
        w_pr = [round(float(w_i * max(1-v**(4/3),0)**0.25 * 1e3), 5) for v in xi_h]
        frames.append({"t_abs": round((t1_off+res1["time"][i])/60,2),
                        "xf": round(xf_i,2), "w_avg": round(w_i*1e3,4),
                        "hf": h_f, "p_wb": round(float(res1["p_f"][i])/1e6,3),
                        "stage": 1,
                        "x_frac": [round(float(v),1) for v in xi],
                        "w_prof": w_pr})

    for i in sub(len(res2["time_since_shutin"]), 70):
        x_n   = res2["x"]
        w_n   = res2["w_profile"][i]
        step  = max(1, len(x_n)//20)
        frames.append({"t_abs": round((t2_off+res2["time_since_shutin"][i])/60,2),
                        "xf": round(float(res2["xf_fixed"]),2),
                        "w_avg": round(float(res2["w_avg"][i])*1e3,4),
                        "hf": h_f, "p_wb": round(float(res2["p_wb"][i])/1e6,3),
                        "stage": 2,
                        "x_frac": [round(float(v),1) for v in x_n[::step]],
                        "w_prof": [round(float(v)*1e3,5) for v in w_n[::step]]})

    # ── Downsample history ──────────────────────────────────────────────
    t_all = np.concatenate([res0["time"],
                             t1_off + res1["time"],
                             np.array(res2["time_absolute"])])
    p_all = np.concatenate([res0["p_wb"], res1["p_f"], res2["p_wb"]])
    q_all = np.concatenate([res0["q"],
                             np.full(len(res1["time"]), p1.q_inj_total),
                             np.zeros(len(res2["time_absolute"]))])

    def ds(arr, n=500):
        if len(arr) <= n: return arr.tolist()
        idx = np.round(np.linspace(0, len(arr)-1, n)).astype(int)
        return arr[idx].tolist()

    return {
        "ok": True,
        "meta": {
            "sigma_hmin_mpa":    round(sig_hmin/1e6, 2),
            "p_breakdown_mpa":   round((sig_hmin+T0)/1e6, 2),
            "p_res_mpa":         round(p_res/1e6, 2),
            "t_sleeve_open_min": round(t1_off/60, 2),
            "t_shutin_min":      round(t2_off/60, 2),
            "xf_final_m":        round(float(res1["xf"][-1]), 1),
            "hf_m":              h_f,
            "initiated":         bool(res1["initiated"]),
            "eff_pct":           round(float(res1["Vf"][-1])/max(float(res1["Vinj_cum"][-1]),1e-9)*100, 1),
            "treating_mpa":      round(float(res1["p_f"][-1])/1e6, 3),
            "isip_mpa":          round(float(res1.get("isip", res1["p_f"][-1]))/1e6, 3),
            "dp_friction_mpa":   round(float(res1.get("dP_friction", 0.0))/1e6, 3),
            "p_final_mpa":       round(float(res2["p_wb"][-1])/1e6, 3),
        },
        "history": {
            "t_min":   [round(v/60,2) for v in ds(t_all)],
            "p_mpa":   [round(v/1e6,3) for v in ds(p_all)],
            "q_m3pm":  [round(v*60,4)  for v in ds(q_all)],
        },
        "pumping": {
            "t_min":   [round(v/60,2) for v in ds(t1_off+res1["time"])],
            "xf_m":    [round(float(v),2) for v in ds(res1["xf"])],
            "w_mm":    [round(float(v)*1e3,4) for v in ds(res1["w_avg"])],
            "eff_pct": [round(float(res1["Vf"][i])/max(float(res1["Vinj_cum"][i]),1e-9)*100,1)
                        for i in range(len(res1["time"]))],
            "Vf":      [round(float(v),3) for v in ds(res1["Vf"])],
            "Vi":      [round(float(v),3) for v in ds(res1["Vinj_cum"])],
            "Vl":      [round(float(v),3) for v in ds(res1["Vleak_cum"])],
            "ql_pump": [round(float(v)*1e3,5) for v in ds(res1["q_leakoff"])],
        },
        "shutin": {
            "t_min":  [round(float(v)/60,2) for v in ds(res2["time_since_shutin"])],
            "p_mpa":  [round(float(v)/1e6,4) for v in ds(res2["p_wb"])],
            "w_mm":   [round(float(v)*1e3,5) for v in ds(res2["w_avg"])],
            "ql":     [round(float(v)*1e3,5) for v in ds(np.array(res2["q_leakoff"]))],
        },
        "frames": frames,
    }


# ======================================================================
# HTML PAGE  (returned once on GET /)
# ======================================================================
PAGE_HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>DFIT Live Dashboard</title>
<style>
*{box-sizing:border-box;margin:0;padding:0}
body{font-family:'Segoe UI',Arial,sans-serif;background:#0f1117;color:#e0e0e0;display:flex;height:100vh;overflow:hidden}

/* LEFT SIDEBAR */
#sidebar{width:260px;min-width:220px;background:#13161f;border-right:1px solid #2a2f45;display:flex;flex-direction:column;overflow:hidden}
#sidebar h2{padding:14px 14px 8px;font-size:.85rem;color:#90a4ae;text-transform:uppercase;letter-spacing:.07em;border-bottom:1px solid #1e2535}
#params{flex:1;overflow-y:auto;padding:10px}
.pg{margin-bottom:10px}
.pg-title{font-size:.72rem;color:#546e7a;text-transform:uppercase;letter-spacing:.06em;margin-bottom:6px;padding:2px 0;border-bottom:1px solid #1e2535}
.row{display:flex;align-items:center;margin-bottom:5px;gap:6px}
.row label{font-size:.75rem;color:#90a4ae;flex:1;min-width:0}
.row input{width:90px;padding:3px 6px;border-radius:4px;border:1px solid #2a2f45;background:#0d1117;color:#e0e0e0;font-size:.78rem;text-align:right}
.row input:focus{outline:none;border-color:#42a5f5}
#run-btn{margin:10px;padding:9px;border-radius:6px;border:none;background:#1565c0;color:#fff;font-size:.85rem;font-weight:600;cursor:pointer;transition:background .15s}
#run-btn:hover{background:#1976d2}
#run-btn:disabled{background:#37474f;cursor:not-allowed}
#status{margin:0 10px 10px;font-size:.75rem;padding:6px 8px;border-radius:4px;min-height:28px;line-height:1.5}
.st-ok{background:#1b3a1b;color:#81c784}
.st-run{background:#1a2744;color:#64b5f6;animation:pulse 1s infinite}
.st-err{background:#3b1212;color:#ef9a9a}
@keyframes pulse{0%,100%{opacity:1}50%{opacity:.5}}
#summary{margin:0 10px 10px;font-size:.72rem;color:#78909c;line-height:1.7;background:#12151f;border-radius:5px;padding:6px 8px}

/* MAIN AREA */
#main{flex:1;display:flex;flex-direction:column;overflow:hidden}
#charts{flex:1;display:grid;grid-template-columns:1fr 1fr;grid-template-rows:1fr 1fr;gap:6px;padding:6px;min-height:0}
.panel{background:#1a1e2e;border-radius:8px;padding:10px;border:1px solid #2a2f45;display:flex;flex-direction:column;min-height:0}
.panel h3{font-size:.75rem;color:#90a4ae;text-transform:uppercase;letter-spacing:.06em;margin-bottom:4px;flex-shrink:0}
.chart-wrap{flex:1;position:relative;min-height:0}
.chart-wrap canvas{position:absolute;inset:0;width:100%;height:100%}
.xctrl{display:flex;align-items:center;gap:5px;margin-top:4px;flex-shrink:0}
.xctrl label{font-size:.68rem;color:#546e7a;white-space:nowrap}
.xctrl input[type=range]{flex:1;accent-color:#42a5f5;height:4px}
.xctrl span{font-size:.68rem;color:#ccc;min-width:34px;text-align:right}

/* ANIMATION ROW */
#anim-row{height:300px;display:flex;gap:6px;padding:0 6px 6px;flex-shrink:0}
.anim-panel{background:#1a1e2e;border-radius:8px;border:1px solid #2a2f45;display:flex;flex-direction:column;padding:8px}
#anim-left{flex:1}
#anim-right{flex:1}
#anim-controls{width:240px;display:flex;flex-direction:column;justify-content:center;gap:8px;padding:6px 0}
.anim-panel h3{font-size:.72rem;color:#90a4ae;text-transform:uppercase;letter-spacing:.06em;margin-bottom:4px;flex-shrink:0}
.anim-canvas-wrap{flex:1;position:relative;min-height:0}
.anim-canvas-wrap canvas{position:absolute;inset:0;width:100%;height:100%}
.btn-row{display:flex;gap:5px;flex-wrap:wrap}
button.ctrl-btn{padding:4px 10px;border-radius:4px;border:none;cursor:pointer;font-size:.78rem;background:#1565c0;color:#fff}
button.ctrl-btn:hover{background:#1976d2}
button.ctrl-btn.active{background:#c62828}
.sctrl{display:flex;align-items:center;gap:5px}
.sctrl label{font-size:.72rem;color:#546e7a}
.sctrl input[type=range]{flex:1;accent-color:#42a5f5}
#info-box{font-size:.72rem;color:#b0bec5;background:#0d1117;border-radius:5px;padding:5px 7px;line-height:1.6;min-height:36px}
.tag{display:inline-block;padding:1px 6px;border-radius:8px;font-size:.65rem;font-weight:700}
.s1{background:#c62828;color:#fff}.s2{background:#2e7d32;color:#fff}

/* Legend */
.legend{display:flex;gap:8px;flex-wrap:wrap;flex-shrink:0}
.leg{display:flex;align-items:center;gap:4px;font-size:.65rem;color:#78909c}
.sw{width:14px;height:3px;border-radius:1px}
</style>
</head>
<body>

<!-- ── SIDEBAR ─────────────────────────────────────────────────────────── -->
<div id="sidebar">
  <h2>⚙ Parameters</h2>
  <div id="params">

    <div class="pg">
      <div class="pg-title">Rock &amp; Fluid</div>
      __PARAM_ROWS_ROCK__
    </div>

    <div class="pg">
      <div class="pg-title">Stage 0 — Pressurisation</div>
      __PARAM_ROWS_S0__
    </div>

    <div class="pg">
      <div class="pg-title">Stage 1 — Propagation</div>
      __PARAM_ROWS_S1__
    </div>

    <div class="pg">
      <div class="pg-title">Stage 2 — Shut-in Decay</div>
      __PARAM_ROWS_S2__
    </div>

    <div class="pg">
      <div class="pg-title">Solver</div>
      __PARAM_ROWS_SOL__
    </div>

  </div>

  <button id="run-btn">▶ Re-run Simulation</button>
  <div id="status" class="st-ok">Ready</div>
  <div id="summary"></div>
</div>

<!-- ── MAIN ─────────────────────────────────────────────────────────────── -->
<div id="main">
  <div id="charts">

    <!-- Panel 1: Operational History -->
    <div class="panel">
      <h3>Operational History</h3>
      <div class="legend">
        <div class="leg"><div class="sw" style="background:#42a5f5"></div>p_wb</div>
        <div class="leg"><div class="sw" style="background:#66bb6a;opacity:.5;height:8px"></div>Rate</div>
        <div class="leg"><div class="sw" style="background:#ef5350"></div>σ_hmin</div>
        <div class="leg"><div class="sw" style="background:#ffa726"></div>p_bd</div>
        <div class="leg"><div class="sw" style="background:#ba68c8"></div>p_res</div>
      </div>
      <div class="chart-wrap"><canvas id="c_hist"></canvas></div>
      <div class="xctrl">
        <label>X max (zoom)</label>
        <input type="range" id="h_xmax" min="0" max="1" step="0.1" value="1">
        <span id="h_xmax_v">—</span>
      </div>
    </div>

    <!-- Panel 2: Fracture dimensions vs time -->
    <div class="panel">
      <h3>Fracture Dimensions vs Time</h3>
      <div class="legend">
        <div class="leg"><div class="sw" style="background:#ef5350"></div>Half-length xf (m)</div>
        <div class="leg"><div class="sw" style="background:#26a69a;border-top:2px dashed #26a69a"></div>Height h_f (m, fixed)</div>
        <div class="leg"><div class="sw" style="background:#ffd54f"></div>Leakoff rate ×1000 (m³/s)</div>
      </div>
      <div class="chart-wrap"><canvas id="c_grow"></canvas></div>
      <div class="xctrl">
        <label>X max (zoom)</label>
        <input type="range" id="g_xmax" min="0" max="1" step="0.1" value="1">
        <span id="g_xmax_v">—</span>
      </div>
    </div>

    <!-- Panel 3: Shut-in Pressure -->
    <div class="panel">
      <h3>Shut-in Pressure Decay</h3>
      <div class="legend">
        <div class="leg"><div class="sw" style="background:#66bb6a"></div>p_wb</div>
        <div class="leg"><div class="sw" style="background:#ef5350"></div>σ_hmin</div>
        <div class="leg"><div class="sw" style="background:#ba68c8"></div>p_res</div>
      </div>
      <div class="chart-wrap"><canvas id="c_shut"></canvas></div>
      <div class="xctrl">
        <label>X max (zoom)</label>
        <input type="range" id="s_xmax" min="0" max="1" step="1" value="1">
        <span id="s_xmax_v">—</span>
      </div>
    </div>

    <!-- Panel 4: Material Balance -->
    <div class="panel">
      <h3>Material Balance &amp; Efficiency</h3>
      <div class="legend">
        <div class="leg"><div class="sw" style="background:#ef5350"></div>Fracture vol</div>
        <div class="leg"><div class="sw" style="background:#546e7a"></div>Leakoff vol</div>
        <div class="leg"><div class="sw" style="background:#000;border:1px dashed #ccc"></div>Injected</div>
        <div class="leg"><div class="sw" style="background:#ffa726"></div>Efficiency %</div>
      </div>
      <div class="chart-wrap"><canvas id="c_mb"></canvas></div>
      <div class="xctrl">
        <label>X max (zoom)</label>
        <input type="range" id="m_xmax" min="0" max="1" step="0.1" value="1">
        <span id="m_xmax_v">—</span>
      </div>
    </div>

  </div><!-- /charts -->

  <!-- ── ANIMATION ROW ─────────────────────────────────────────────── -->
  <div id="anim-row">
    <div class="anim-panel" id="anim-left">
      <h3>Fracture Cross-Section (Side View)</h3>
      <div class="anim-canvas-wrap"><canvas id="c_anim"></canvas></div>
    </div>
    <div class="anim-panel" id="anim-right">
      <h3>Aperture profile along fracture</h3>
      <div class="anim-canvas-wrap"><canvas id="c_prof"></canvas></div>
    </div>
    <div id="anim-controls">
      <div class="btn-row">
        <button class="ctrl-btn" id="btn_play">▶ Play</button>
        <button class="ctrl-btn" id="btn_step">⏭</button>
        <button class="ctrl-btn" id="btn_reset">⏮</button>
      </div>
      <div class="sctrl">
        <label>Speed</label>
        <input type="range" id="speed" min="1" max="20" value="6" style="width:80px">
        <span id="speed_v" style="font-size:.72rem;color:#ccc">6×</span>
      </div>
      <div class="sctrl">
        <label>Frame</label>
        <input type="range" id="fslider" min="0" max="0" value="0" style="flex:1">
        <span id="f_v" style="font-size:.72rem;color:#ccc">0</span>
      </div>
      <div id="info-box">Press ▶ Play or Re-run to start</div>
    </div>
  </div>

</div><!-- /main -->

<script>
// ── Global state ─────────────────────────────────────────────────────────────
let DATA = null;
const DPR = window.devicePixelRatio || 1;

// ── Utilities ────────────────────────────────────────────────────────────────
function scaleCanvas(c) {
  const r = c.getBoundingClientRect();
  c.width  = r.width  * DPR;
  c.height = r.height * DPR;
  const ctx = c.getContext('2d');
  ctx.scale(DPR, DPR);
  return [ctx, r.width, r.height];
}
function lin(d, r) { return v => r[0] + (v-d[0])/(d[1]-d[0])*(r[1]-r[0]); }
function clamp(v,a,b){return Math.max(a,Math.min(b,v));}

function drawGrid(ctx,W,H,pad,xMn,xMx,yMn,yMx,xLab,yLab,nX,nY){
  const [L,R,T,B]=[pad[0],W-pad[1],pad[2],H-pad[3]];
  const sx=lin([xMn,xMx],[L,R]), sy=lin([yMx,yMn],[T,B]);
  ctx.fillStyle='#0d1117'; ctx.fillRect(0,0,W,H);
  ctx.fillStyle='#12151f'; ctx.fillRect(L,T,R-L,B-T);
  ctx.strokeStyle='#1e2535'; ctx.lineWidth=1;
  for(let i=0;i<=nX;i++){const px=sx(xMn+i*(xMx-xMn)/nX);
    ctx.beginPath();ctx.moveTo(px,T);ctx.lineTo(px,B);ctx.stroke();}
  for(let i=0;i<=nY;i++){const py=sy(yMn+i*(yMx-yMn)/nY);
    ctx.beginPath();ctx.moveTo(L,py);ctx.lineTo(R,py);ctx.stroke();}
  ctx.strokeStyle='#37474f';ctx.lineWidth=1.5;ctx.strokeRect(L,T,R-L,B-T);
  ctx.fillStyle='#78909c';ctx.font='10px Arial';ctx.textAlign='center';
  for(let i=0;i<=nX;i++){const v=xMn+i*(xMx-xMn)/nX;
    ctx.fillText(v<10?v.toFixed(1):v.toFixed(0),sx(v),B+13);}
  ctx.textAlign='right';
  for(let i=0;i<=nY;i++){const v=yMn+i*(yMx-yMn)/nY;
    ctx.fillText(v.toFixed(v<5?2:1),L-4,sy(v)+3);}
  ctx.fillStyle='#90a4ae';ctx.font='bold 10px Arial';
  ctx.textAlign='center';ctx.fillText(xLab,L+(R-L)/2,H-2);
  ctx.save();ctx.translate(10,T+(B-T)/2);ctx.rotate(-Math.PI/2);
  ctx.fillText(yLab,0,0);ctx.restore();
  return {sx,sy,L,R,T,B};
}

function polyline(ctx,xs,ys,sx,sy,col,lw=2,dash=[]){
  if(!xs||!xs.length)return;
  ctx.setLineDash(dash);ctx.beginPath();ctx.strokeStyle=col;ctx.lineWidth=lw;
  ctx.moveTo(sx(xs[0]),sy(ys[0]));
  for(let i=1;i<xs.length;i++)ctx.lineTo(sx(xs[i]),sy(ys[i]));
  ctx.stroke();ctx.setLineDash([]);
}
function fillArea(ctx,xs,ys,sx,sy,col,al=.2,base=0){
  if(!xs||!xs.length)return;
  ctx.beginPath();ctx.moveTo(sx(xs[0]),sy(base));
  xs.forEach((x,i)=>ctx.lineTo(sx(x),sy(ys[i])));
  ctx.lineTo(sx(xs[xs.length-1]),sy(base));ctx.closePath();
  ctx.fillStyle=col;ctx.globalAlpha=al;ctx.fill();ctx.globalAlpha=1;
}
function hline(ctx,y,sy,L,R,col,dash=[]){
  ctx.setLineDash(dash);ctx.strokeStyle=col;ctx.lineWidth=1.5;
  ctx.beginPath();ctx.moveTo(L,sy(y));ctx.lineTo(R,sy(y));ctx.stroke();ctx.setLineDash([]);
}
function vline(ctx,x,sx,T,B,col){
  ctx.setLineDash([4,3]);ctx.strokeStyle=col;ctx.lineWidth=1;
  ctx.beginPath();ctx.moveTo(sx(x),T);ctx.lineTo(sx(x),B);ctx.stroke();ctx.setLineDash([]);
}
function stageBg(ctx,sx,T,B,t1,t2,xMn,xMx){
  const segs=[[Math.max(xMn,0),Math.min(xMx,t1),'#1565c0',.07],
               [Math.max(xMn,t1),Math.min(xMx,t2),'#c62828',.07],
               [Math.max(xMn,t2),xMx,'#2e7d32',.07]];
  segs.forEach(([a,b,c,al])=>{if(b>a){ctx.fillStyle=c;ctx.globalAlpha=al;
    ctx.fillRect(sx(a),T,sx(b)-sx(a),B-T);ctx.globalAlpha=1;}});
}

// ── Slider helper ────────────────────────────────────────────────────────────
function makeSliders(pfx, getMax, draw) {
  const s2=document.getElementById(pfx+'_xmax');
  const v2=document.getElementById(pfx+'_xmax_v');
  function update(){
    const mx=getMax(); s2.max=mx;
    const b=+s2.value;
    v2.textContent=b.toFixed(b<10?1:0);
    draw();
  }
  s2.addEventListener('input',update);
  return {
    reset(mx,step){s2.max=mx;s2.step=step;s2.value=mx;update();},
    get(){return [0,+s2.value];}
  };
}

// ── DRAW FUNCTIONS ───────────────────────────────────────────────────────────
const PAD=[44,12,18,28];

let sl_h, sl_g, sl_s, sl_m;

function drawHist(){
  if(!DATA)return;
  const [xMn,xMx]=sl_h.get(); if(xMn>=xMx)return;
  const c=document.getElementById('c_hist');
  const [ctx,W,H]=scaleCanvas(c);
  const m=DATA.meta, ht=DATA.history;
  const pAll=[...ht.p_mpa,m.sigma_hmin_mpa,m.p_breakdown_mpa,m.p_res_mpa];
  const pMx=Math.max(...pAll)*1.04, pMn=Math.min(...pAll)*0.96;
  const {sx,sy,L,R,T,B}=drawGrid(ctx,W,H,PAD,xMn,xMx,pMn,pMx,'Time (min)','Pressure (MPa)',6,5);
  stageBg(ctx,sx,T,B,m.t_sleeve_open_min,m.t_shutin_min,xMn,xMx);
  const qMx=Math.max(...ht.q_m3pm)||1;
  const toP=v=>pMn+v/qMx*(pMx-pMn)*0.28;
  fillArea(ctx,ht.t_min,ht.q_m3pm.map(toP),sx,sy,'#66bb6a',.35,pMn);
  hline(ctx,m.sigma_hmin_mpa,sy,L,R,'#ef5350',[6,3]);
  hline(ctx,m.p_breakdown_mpa,sy,L,R,'#ffa726',[3,3]);
  if(m.isip_mpa) hline(ctx,m.isip_mpa,sy,L,R,'#ffd54f',[4,4]);
  hline(ctx,m.p_res_mpa,sy,L,R,'#ba68c8',[2,5]);
  if(m.t_sleeve_open_min>xMn&&m.t_sleeve_open_min<xMx)vline(ctx,m.t_sleeve_open_min,sx,T,B,'#ffa726');
  if(m.t_shutin_min>xMn&&m.t_shutin_min<xMx)vline(ctx,m.t_shutin_min,sx,T,B,'#66bb6a');
  polyline(ctx,ht.t_min,ht.p_mpa,sx,sy,'#42a5f5',2.2);
}

function drawGrow(){
  if(!DATA)return;
  const [xMn,xMx]=sl_g.get(); if(xMn>=xMx)return;
  const c=document.getElementById('c_grow');
  const [ctx,W,H]=scaleCanvas(c);
  const pu=DATA.pumping, m=DATA.meta;
  const tRel=pu.t_min.map(v=>v-pu.t_min[0]);
  const xfMx=Math.max(...pu.xf_m,m.hf_m)*1.12;
  const {sx,sy,L,R,T,B}=drawGrid(ctx,W,H,PAD,xMn,xMx,0,xfMx,'Pumping Time (min)','Length (m)',6,5);
  // Height as dashed line
  hline(ctx,m.hf_m,sy,L,R,'#26a69a',[6,3]);
  ctx.fillStyle='#26a69a';ctx.font='9px Arial';ctx.textAlign='left';
  ctx.fillText('h_f = '+m.hf_m+' m',sx(xMn)+3,sy(m.hf_m)-4);
  // Half-length filled area
  fillArea(ctx,tRel,pu.xf_m,sx,sy,'#ef5350',.15,0);
  polyline(ctx,tRel,pu.xf_m,sx,sy,'#ef5350',2.5);
  // Leakoff rate (scaled to fit same axis)
  const qlMax=Math.max(...pu.ql_pump||[1e-9]);
  const qlScale=pu.ql_pump?pu.ql_pump.map(v=>v/qlMax*xfMx*0.35):[];
  if(qlScale.length) polyline(ctx,tRel,qlScale,sx,sy,'#ffd54f',1.8,[4,3]);
  // Leakoff axis label (right side)
  if(pu.ql_pump&&pu.ql_pump.length){
    ctx.fillStyle='#ffd54f';ctx.textAlign='right';ctx.font='9px Arial';
    ctx.fillText('Leakoff →',R-2,sy(qlScale[qlScale.length-1])-4);
  }
}

function drawShut(){
  if(!DATA)return;
  const [xMn,xMx]=sl_s.get(); if(xMn>=xMx)return;
  const c=document.getElementById('c_shut');
  const [ctx,W,H]=scaleCanvas(c);
  const si=DATA.shutin, m=DATA.meta;
  const pMx=Math.max(...si.p_mpa)*1.02;
  const pMn=Math.min(m.p_res_mpa*.98,...si.p_mpa)*.98;
  const {sx,sy,L,R,T,B}=drawGrid(ctx,W,H,PAD,xMn,xMx,pMn,pMx,'Time since shut-in (min)','Pressure (MPa)',6,5);
  hline(ctx,m.sigma_hmin_mpa,sy,L,R,'#ef5350',[6,3]);
  hline(ctx,m.p_res_mpa,sy,L,R,'#ba68c8',[3,3]);
  polyline(ctx,si.t_min,si.p_mpa,sx,sy,'#66bb6a',2.2);
}

function drawMB(){
  if(!DATA)return;
  const [xMn,xMx]=sl_m.get(); if(xMn>=xMx)return;
  const c=document.getElementById('c_mb');
  const [ctx,W,H]=scaleCanvas(c);
  const pu=DATA.pumping;
  const tRel=pu.t_min.map(v=>v-pu.t_min[0]);
  const vMx=Math.max(...pu.Vi)*1.05;
  const {sx,sy,L,R,T,B}=drawGrid(ctx,W,H,PAD,xMn,xMx,0,vMx,'Pumping Time (min)','Volume (m³)',6,5);
  // Stacked: leakoff on top of fracture
  const Vstack=pu.Vf.map((v,i)=>v+pu.Vl[i]);
  fillArea(ctx,tRel,Vstack,sx,sy,'#546e7a',.5,0);
  fillArea(ctx,tRel,pu.Vf,sx,sy,'#ef5350',.5,0);
  polyline(ctx,tRel,pu.Vi,sx,sy,'#cfd8dc',1.5,[6,3]);
  // Efficiency on secondary axis (0-100%)
  const sy2=lin([100,0],[T,B]);
  polyline(ctx,tRel,pu.eff_pct,sx,sy2,'#ffa726',1.5);
  ctx.fillStyle='#ffa726';ctx.font='9px Arial';ctx.textAlign='right';
  ctx.fillText('Eff %',R-2,sy2(pu.eff_pct[pu.eff_pct.length-1])-4);
}

// ── Animation ────────────────────────────────────────────────────────────────
let frameIdx=0, playing=false, animRAF=null, lastTs=0;
const btnPlay=document.getElementById('btn_play');
const fslider=document.getElementById('fslider');
const speedSl=document.getElementById('speed');
speedSl.addEventListener('input',()=>{document.getElementById('speed_v').textContent=speedSl.value+'×';});

function drawFrame(fi){
  if(!DATA||!DATA.frames.length)return;
  fi=clamp(fi,0,DATA.frames.length-1);
  frameIdx=fi; fslider.value=fi;
  document.getElementById('f_v').textContent=fi;
  const fr=DATA.frames[fi], m=DATA.meta;

  // Info
  const sl=fr.stage===1?'<span class="tag s1">Propagation</span>':'<span class="tag s2">Shut-in</span>';
  document.getElementById('info-box').innerHTML=
    `t=<b>${fr.t_abs.toFixed(1)} min</b> ${sl} &nbsp;xf=<b>${fr.xf.toFixed(0)} m</b>`+
    ` &nbsp;w̄=<b>${fr.w_avg.toFixed(3)} mm</b> &nbsp;p=<b>${fr.p_wb.toFixed(2)} MPa</b>`;

  // ── Fracture Cross-Section: side view (depth vs lateral distance) ─────────
  const cA=document.getElementById('c_anim');
  const [ctxA,WA,HA]=scaleCanvas(cA);
  const mg=44, mgT=30, mgB=28;
  const xfFin=m.xf_final_m, hf=m.hf_m;
  const xfC=fr.xf, wC=fr.w_avg;
  const wMaxAll=Math.max(...DATA.frames.map(f=>f.w_avg));

  // Depth: wellbore centre at depth_centre; fracture spans ±hf/2
  const depth_centre = 2740; // m ≈ 9000 ft (display depth)
  const depth_pad    = hf * 3.5;
  const yMin_d = depth_centre - depth_pad;
  const yMax_d = depth_centre + depth_pad;

  // Axes: x = lateral distance from wellbore (0 → xf_final), y = depth (yMin→yMax)
  const sxA = lin([0, xfFin],[mg, WA-10]);
  const syA = lin([yMin_d, yMax_d],[mgT, HA-mgB]);

  ctxA.fillStyle='#0d1117'; ctxA.fillRect(0,0,WA,HA);
  ctxA.fillStyle='#12151f'; ctxA.fillRect(mg,mgT,WA-mg-10,HA-mgT-mgB);

  // Background strata lines (formation layers)
  const strata=[
    [yMin_d + depth_pad*0.3,'#1a2030','Overburden'],
    [depth_centre - hf/2,'#1e2840','Pay Zone'],
    [depth_centre + hf/2,'#1a2030','Underburden'],
  ];
  let prev_d=yMin_d;
  const strata_cols=['#111827','#162032','#111827'];
  for(let si=0;si<strata.length;si++){
    const y0=si===0?mgT:syA(strata[si-1][0]);
    const y1=si===strata.length-1?HA-mgB:syA(strata[si][0]);
    ctxA.fillStyle=strata_cols[si]; ctxA.fillRect(mg,y0,WA-mg-10,y1-y0);
    if(si<strata.length-1){
      ctxA.strokeStyle='#263040';ctxA.lineWidth=1;ctxA.setLineDash([4,3]);
      ctxA.beginPath();ctxA.moveTo(mg,syA(strata[si][0]));ctxA.lineTo(WA-10,syA(strata[si][0]));ctxA.stroke();
      ctxA.setLineDash([]);
    }
  }

  // Pay zone highlight
  const y_top=syA(depth_centre-hf/2), y_bot=syA(depth_centre+hf/2);
  ctxA.fillStyle='rgba(33,96,160,0.12)';
  ctxA.fillRect(mg, y_top, WA-mg-10, y_bot-y_top);

  // Grid lines (x only – lateral distance)
  ctxA.strokeStyle='#1e2535';ctxA.lineWidth=1;
  for(let i=1;i<=5;i++){
    const px=sxA(xfFin*i/5);
    ctxA.beginPath();ctxA.moveTo(px,mgT);ctxA.lineTo(px,HA-mgB);ctxA.stroke();
  }

  // ── PKN Fracture: TRUE RECTANGLE — constant h_f, perfectly flat tip ─────
  // PKN model: height is FIXED by pay zone at every lateral position.
  // Aperture (opening width) varies along x per the PKN profile but HEIGHT
  // is constant → the cross-section is a perfect rectangle with a flat tip.
  const nS   = 160;
  const alpha = Math.max(wC / wMaxAll, 0.05);
  const py_top_rect = syA(depth_centre - hf/2);
  const py_bot_rect = syA(depth_centre + hf/2);
  const rect_h      = py_bot_rect - py_top_rect;
  const px_left     = sxA(0);
  const px_right    = sxA(xfC);

  // 1. Fill with aperture gradient (column by column, full height each slice)
  for(let s=0; s<nS; s++){
    const xi   = xfC * s / nS;
    const xi1  = xfC * (s+1) / nS;
    const xhat = (xi + xi1) * 0.5 / Math.max(xfC, 0.001);
    let r, g, b;
    if(fr.stage === 1){
      // Propagation: warm red-orange at wellbore → cool steel-blue at tip
      r = Math.round(200*(1-xhat) + 30*xhat);
      g = Math.round(50*(1-xhat)  + 80*xhat);
      b = Math.round(30*(1-xhat)  + 170*xhat);
    } else {
      // Shut-in: teal-green at wellbore → near-black at tip
      r = Math.round(15*(1-xhat)  + 5*xhat);
      g = Math.round(150*(1-xhat) + 40*xhat);
      b = Math.round(100*(1-xhat) + 30*xhat);
    }
    ctxA.fillStyle = `rgba(${r},${g},${b},${alpha*0.93})`;
    ctxA.fillRect(sxA(xi), py_top_rect, Math.max(sxA(xi1)-sxA(xi), 1), rect_h);
  }

  // 2. Perfect rectangle outline — strokeRect gives exactly 4 straight sides
  const col_outline = fr.stage===1 ? '#ef5350' : '#66bb6a';
  ctxA.strokeStyle = col_outline;
  ctxA.lineWidth   = 2.5;
  ctxA.setLineDash([]);
  ctxA.strokeRect(px_left, py_top_rect, px_right - px_left, rect_h);

  // 3. Faint vertical hint-lines encoding aperture profile
  for(let i=1; i<=10; i++){
    const xhat2 = i / 11;
    const wf    = Math.max(1 - xhat2**(4/3), 0)**0.25;
    const px    = sxA(xfC * xhat2);
    ctxA.strokeStyle = `rgba(255,255,255,${wf * 0.06 * alpha})`;
    ctxA.lineWidth   = 1;
    ctxA.beginPath(); ctxA.moveTo(px, py_top_rect); ctxA.lineTo(px, py_bot_rect);
    ctxA.stroke();
  }

  // 4. Bright wellbore-face line (fracture mouth, maximum aperture)
  ctxA.strokeStyle = fr.stage===1 ? 'rgba(255,160,80,0.65)' : 'rgba(80,220,120,0.65)';
  ctxA.lineWidth   = 3;
  ctxA.beginPath(); ctxA.moveTo(px_left, py_top_rect); ctxA.lineTo(px_left, py_bot_rect);
  ctxA.stroke();
  // Wellbore line
  ctxA.strokeStyle='#ffa726';ctxA.lineWidth=3;
  ctxA.beginPath();ctxA.moveTo(mg,mgT);ctxA.lineTo(mg,HA-mgB);ctxA.stroke();

  // Axis labels
  ctxA.fillStyle='#546e7a';ctxA.font='10px Arial';ctxA.textAlign='center';
  ctxA.fillText(`Lateral Distance: ${xfC.toFixed(0)} m`,mg+(WA-mg-10)/2,HA-5);
  // X axis ticks
  ctxA.fillStyle='#546e7a';ctxA.font='9px Arial';
  for(let i=0;i<=5;i++){
    const xv=xfFin*i/5;
    ctxA.fillText(Math.round(xv)+'',sxA(xv),HA-mgB+11);
  }
  // Y axis ticks (depth in ft)
  ctxA.textAlign='right';
  for(let d=Math.ceil(yMin_d/50)*50;d<=yMax_d;d+=50){
    const py=syA(d); if(py<mgT||py>HA-mgB) continue;
    const ft=Math.round(d*3.28084/10)*10;
    ctxA.fillText(ft+'',mg-3,py+3);
  }
  ctxA.save();ctxA.translate(10,mgT+(HA-mgT-mgB)/2);ctxA.rotate(-Math.PI/2);
  ctxA.textAlign='center';ctxA.fillStyle='#90a4ae';ctxA.font='bold 10px Arial';
  ctxA.fillText('Depth (ft)',0,0);ctxA.restore();

  // Labels
  ctxA.textAlign='right';ctxA.fillStyle='#ffa726';ctxA.font='bold 10px Arial';
  ctxA.fillText('Wellbore',mg-3,mgT+10);
  ctxA.textAlign='center';
  ctxA.fillStyle='rgba(33,96,160,0.7)';ctxA.font='9px Arial';
  ctxA.fillText('Pay Zone',mg+50,y_top+10);

  ctxA.fillStyle=fr.stage===1?'#ef5350':'#66bb6a';
  ctxA.font='bold 11px Arial';
  ctxA.fillText(fr.stage===1?'▶ Propagating':'◀ Closing / Aperture Narrowing',mg+(WA-mg)/2,mgT-8);

  // Aperture profile
  const cP=document.getElementById('c_prof');
  const [ctxP,WP,HP]=scaleCanvas(cP);
  const wMxP=Math.max(...DATA.frames.map(f=>Math.max(...f.w_prof)));
  const {sx:spx,sy:spy,L:Lp,R:Rp}=drawGrid(ctxP,WP,HP,PAD,0,xfFin,0,wMxP*1.1,
    'Distance from wellbore (m)','Width (mm)',5,5);
  const xp=fr.x_frac, wp=fr.w_prof;
  fillArea(ctxP,xp,wp,spx,spy,fr.stage===1?'rgba(239,83,80,.3)':'rgba(102,187,106,.3)',.9,0);
  polyline(ctxP,xp,wp,spx,spy,fr.stage===1?'#ef5350':'#66bb6a',2.5);
  ctxP.setLineDash([4,3]);ctxP.strokeStyle='#ffa726';ctxP.lineWidth=1.5;
  ctxP.beginPath();ctxP.moveTo(Lp,spy(fr.w_avg));ctxP.lineTo(Rp,spy(fr.w_avg));ctxP.stroke();ctxP.setLineDash([]);
  ctxP.fillStyle='#ffa726';ctxP.font='9px Arial';ctxP.textAlign='left';
  ctxP.fillText(`w̄=${fr.w_avg.toFixed(3)} mm`,Lp+3,spy(fr.w_avg)-3);
}

function animTick(ts){
  if(!playing)return;
  const ms=1000/((+speedSl.value)*15);
  if(ts-lastTs>ms){lastTs=ts;frameIdx=(frameIdx+1)%DATA.frames.length;drawFrame(frameIdx);}
  animRAF=requestAnimationFrame(animTick);
}
btnPlay.addEventListener('click',()=>{
  playing=!playing;btnPlay.textContent=playing?'⏸ Pause':'▶ Play';
  btnPlay.classList.toggle('active',playing);
  if(playing)animRAF=requestAnimationFrame(animTick);
  else cancelAnimationFrame(animRAF);
});
document.getElementById('btn_step').addEventListener('click',()=>{frameIdx=Math.min(frameIdx+1,DATA.frames.length-1);drawFrame(frameIdx);});
document.getElementById('btn_reset').addEventListener('click',()=>{
  playing=false;btnPlay.textContent='▶ Play';btnPlay.classList.remove('active');
  cancelAnimationFrame(animRAF);frameIdx=0;drawFrame(0);
});
fslider.addEventListener('input',()=>{frameIdx=+fslider.value;drawFrame(frameIdx);});

// ── Fetch simulation data from server ────────────────────────────────────────
function loadData(d){
  DATA=d;
  const m=d.meta;
  // Set slider ranges
  const tMax=d.history.t_min[d.history.t_min.length-1];
  const t1Dur=d.pumping.t_min[d.pumping.t_min.length-1]-d.pumping.t_min[0];
  const siMax=d.shutin.t_min[d.shutin.t_min.length-1];
  sl_h.reset(tMax.toFixed(1),0.5);
  sl_g.reset(t1Dur.toFixed(1),0.5);
  sl_s.reset(siMax.toFixed(0),5);
  sl_m.reset(t1Dur.toFixed(1),0.5);
  fslider.max=DATA.frames.length-1; fslider.value=0;
  frameIdx=0;
  drawHist();drawGrow();drawShut();drawMB();drawFrame(0);
  document.getElementById('summary').innerHTML=
    `<b>Results:</b><br>`+
    `xf = ${m.xf_final_m} m &nbsp;|&nbsp; h_f = ${m.hf_m} m<br>`+
    `Treating = ${m.treating_mpa} MPa &nbsp;|&nbsp; <span style="color:#ffd54f">ΔP fric = ${m.dp_friction_mpa} MPa</span><br>`+
    `ISIP = ${m.isip_mpa} MPa &nbsp;|&nbsp; Eff = ${m.eff_pct}%<br>`+
    `Initiated: ${m.initiated?'✓ Yes':'✗ No (check pressurisation)'}<br>`+
    `p at end = ${m.p_final_mpa} MPa (→ ${m.p_res_mpa} MPa)`;
}

// ── Run button ────────────────────────────────────────────────────────────────
document.getElementById('run-btn').addEventListener('click', async ()=>{
  // Collect current input values
  const inputs={};
  document.querySelectorAll('#params input[data-key]').forEach(el=>{
    inputs[el.dataset.key]=el.value;
  });
  const btn=document.getElementById('run-btn');
  const st=document.getElementById('status');
  btn.disabled=true; btn.textContent='⏳ Running...';
  st.className='st-run'; st.textContent='Simulation running…';
  try {
    const resp=await fetch('/run', {method:'POST',
      headers:{'Content-Type':'application/json'},
      body:JSON.stringify(inputs)});
    const d=await resp.json();
    if(d.ok){
      st.className='st-ok'; st.textContent='✓ Done';
      loadData(d);
    } else {
      st.className='st-err'; st.textContent='Error: '+d.error;
    }
  } catch(e){
    st.className='st-err'; st.textContent='Connection error: '+e;
  }
  btn.disabled=false; btn.textContent='▶ Re-run Simulation';
});

// ── Init sliders after DOM ready ─────────────────────────────────────────────
window.addEventListener('load', ()=>{
  sl_h = makeSliders('h', ()=>DATA?DATA.history.t_min[DATA.history.t_min.length-1]:1, drawHist);
  sl_g = makeSliders('g', ()=>DATA?DATA.pumping.t_min[DATA.pumping.t_min.length-1]-DATA.pumping.t_min[0]:1, drawGrow);
  sl_s = makeSliders('s', ()=>DATA?DATA.shutin.t_min[DATA.shutin.t_min.length-1]:1, drawShut);
  sl_m = makeSliders('m', ()=>DATA?DATA.pumping.t_min[DATA.pumping.t_min.length-1]-DATA.pumping.t_min[0]:1, drawMB);

  // Auto-run on load
  document.getElementById('run-btn').click();
});
window.addEventListener('resize', ()=>{drawHist();drawGrow();drawShut();drawMB();if(DATA)drawFrame(frameIdx);});
</script>
</body>
</html>
"""

# ======================================================================
# HTTP SERVER
# ======================================================================
# Parameter metadata: key -> (label, group, step)
PARAM_META = {
    "E_GPa":               ("E (GPa)",                   "rock", "0.5"),
    "nu":                  ("ν (Poisson ratio)",          "rock", "0.01"),
    "h_f_m":               ("h_f (m)",                   "rock", "1"),
    "sigma_hmin_MPa":      ("σ_hmin (MPa)",               "rock", "0.5"),
    "tensile_strength_MPa":("Tensile strength T₀ (MPa)",  "rock", "0.5"),
    "mu_cP":               ("Viscosity μ (cP)",            "rock", "0.1"),
    "p_initial_MPa":       ("p_initial (MPa)",             "s0",   "0.5"),
    "q_press_m3s":         ("q_press (m³/s)",              "s0",   "1e-5"),
    "C_wb_m3Pa":           ("C_wb (m³/Pa)",                "s0",   "1e-10"),
    "ramp_time_s":         ("Ramp time (s)",               "s0",   "5"),
    "t_sleeve_open_s":     ("t_sleeve_open (s)",           "s0",   "10"),
    "q_inj_m3s":           ("q_inject (m³/s)",             "s1",   "0.001"),
    "t_pump_min":          ("t_pump (min)",                "s1",   "1"),
    "C_L_m_sqrts":         ("C_L (m/√s)",                 "s1",   "1e-7"),
    "dP_total_MPa":        ("dP_friction (MPa)",           "s1",   "0.5"),
    "p_res_MPa":           ("p_res (MPa)",                 "s2",   "0.5"),
    "t_shutin_hr":         ("t_shutin (hr)",               "s2",   "1"),
    "Eres_mm":             ("E_res (mm)",                  "s2",   "0.01"),
    "nx":                  ("nx (nodes)",                  "sol",  "1"),
    "dt_shutin_s":         ("dt_shutin (s)",               "sol",  "10"),
}
GROUP_PLACEHOLDER = {"rock":"__PARAM_ROWS_ROCK__","s0":"__PARAM_ROWS_S0__",
                     "s1":"__PARAM_ROWS_S1__","s2":"__PARAM_ROWS_S2__","sol":"__PARAM_ROWS_SOL__"}

def build_param_rows():
    groups = {g: [] for g in ["rock","s0","s1","s2","sol"]}
    for key, (label, grp, step) in PARAM_META.items():
        val = PARAMETERS[key]
        row = (f'<div class="row"><label>{label}</label>'
               f'<input type="text" data-key="{key}" value="{val}" step="{step}"></div>')
        groups[grp].append(row)
    html = PAGE_HTML
    for grp, ph in GROUP_PLACEHOLDER.items():
        html = html.replace(ph, '\n'.join(groups[grp]))
    return html

PAGE = build_param_rows()


class Handler(BaseHTTPRequestHandler):
    def log_message(self, fmt, *args):
        pass  # silence access log

    def do_GET(self):
        self.send_response(200)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.end_headers()
        self.wfile.write(PAGE.encode())

    def do_POST(self):
        length = int(self.headers.get("Content-Length", 0))
        body   = self.rfile.read(length)
        try:
            raw = json.loads(body)
            # Convert string values to float/int
            params = {}
            for k, v in raw.items():
                try:   params[k] = int(v) if k in ("nx",) else float(v)
                except: params[k] = v
            result = run_simulation(params)
        except Exception as e:
            result = {"ok": False, "error": traceback.format_exc()[-400:]}

        out = json.dumps(result, separators=(",", ":")).encode()
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", len(out))
        self.end_headers()
        self.wfile.write(out)


if __name__ == "__main__":
    PORT = 5050
    server = HTTPServer(("127.0.0.1", PORT), Handler)
    print(f"""
╔══════════════════════════════════════════════════════╗
║        DFIT Live Dashboard  —  ready                 ║
╠══════════════════════════════════════════════════════╣
║  Open in browser:  http://localhost:{PORT}              ║
║                                                      ║
║  Edit parameters in the PARAMETERS dict at the top   ║
║  of this file, then click "Re-run Simulation" in     ║
║  the browser — charts update instantly.              ║
║                                                      ║
║  Press  Ctrl+C  to stop the server.                  ║
╚══════════════════════════════════════════════════════╝
""")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nServer stopped.")
