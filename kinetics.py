"""CoreForge transient module — point reactor kinetics & accidents.

Time-dependent core power for reactivity events (rod ejection, rod drop,
ramps) and full design-basis accident sequences:

    dn/dt   = (rho(t) - beta)/Lambda * n + sum_i lambda_i c_i
    dc_i/dt = beta_i/Lambda * n - lambda_i c_i
    rho(t)  = rho_ext(t) + rho_scram(t) + rho_feedback(T_fuel, T_mod)

Reactivity feedback (optional, see simulate):
  * Doppler, linear  alpha_D (T_f - T_f0)  or physical resonance-broadening
    sqrt form  A (sqrt(T_f) - sqrt(T_f0)) with A matched to alpha_D at T_f0;
  * moderator (MTC)  alpha_M (T_m - T_m0) on a TWO-NODE fuel+moderator
    thermal model whose nominal temperatures follow from the steady heat
    balance at P0 (so both feedbacks are zero at nominal).
Reactor protection (optional): a high-flux / short-period trip latches the
first setpoint crossing and inserts a scram bank over the rod-drop time,
turning an initiating event into a protected SEQUENCE (REA) — or, with the
trip removed, an anticipated-transient-without-scram (ATWS).

Numerics: over each output step the kinetics matrix is FROZEN and the
exact solution y(t+dt) = expm(A dt) y(t) is applied via eigen-
decomposition (numpy only) — unconditionally stable, no stiffness
issues from the prompt mode.  The thermal ODEs are advanced with the
same step (their timescale is seconds, the step is milliseconds).

Validation (see verify.py): for constant rho the simulated asymptotic
period must match the analytic INHOUR root; a dedicated accident suite
checks trip→scram shutdown, ATWS vs REA, feedback signs and the
prompt-critical / energy metrics.

Delayed-neutron data: 6-group U-235 thermal fission set (Keepin),
beta_eff = 0.0065; prompt generation time Lambda user-settable
(PWR ~ 2e-5 s).
"""
import numpy as np

# 6-group U-235 delayed data (Keepin, thermal fission)
BETA_I = np.array([0.000215, 0.001424, 0.001274,
                   0.002568, 0.000748, 0.000273])
LAM_I = np.array([0.0124, 0.0305, 0.111, 0.301, 1.14, 3.01])
BETA = float(BETA_I.sum())          # 0.006502


def _pk_matrix(rho, Lambda, betas=BETA_I, lams=LAM_I):
    """Kinetics matrix A for y = [n, c1..c6] (c_i scaled to n units)."""
    m = len(betas)
    beta = betas.sum()
    A = np.zeros((m + 1, m + 1))
    A[0, 0] = (rho - beta) / Lambda
    A[0, 1:] = lams
    A[1:, 0] = betas / Lambda
    A[1:, 1:] = -np.diag(lams)
    return A


def _expm(A, dt):
    """Matrix exponential via eigen-decomposition (A is diagonalizable
    for physical kinetics parameters)."""
    w, V = np.linalg.eig(A)
    return (V @ np.diag(np.exp(w * dt)) @ np.linalg.inv(V)).real


def equilibrium_y(n0=1.0, Lambda=2e-5, betas=BETA_I, lams=LAM_I):
    """Steady-state precursor populations for power n0 at rho = 0."""
    y = np.empty(len(betas) + 1)
    y[0] = n0
    y[1:] = betas * n0 / (Lambda * lams)
    return y


def inhour_omega(rho, Lambda=2e-5, betas=BETA_I, lams=LAM_I):
    """Largest root omega of the inhour equation
        rho = omega*Lambda + omega * sum_i beta_i/(omega + lambda_i)
    (asymptotic period T = 1/omega).  Bisection on the physical branch."""
    def f(w):
        return w * Lambda + np.sum(betas * w / (w + lams)) - rho

    if abs(rho) < 1e-12:
        return 0.0
    if rho > 0:
        lo, hi = 1e-12, 1.0
        while f(hi) < 0:
            hi *= 2.0
            if hi > 1e8:
                raise RuntimeError("inhour: no positive root found")
    else:
        # largest negative root lies in (-lambda_min, 0)
        lo, hi = -min(lams) + 1e-9, -1e-12
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ---------------------------------------------------------------------
# Xenon-135 transient (iodine pit / load-follow), 0-D
# ---------------------------------------------------------------------
LAM_I135 = 2.87e-5           # I-135 decay  [1/s]
LAM_XE135 = 2.093e-5         # Xe-135 decay [1/s]
Y_I135 = 0.0639              # direct I-135 yield (dominant path to Xe)
Y_XE_DIR = 0.0025            # small direct Xe-135 yield


def xenon_transient(power_frac_final=0.0, t_end_h=72.0,
                    sigphi0=5.3e-5, rho_eq_pcm=-2800.0, npts=800):
    """Xe-135 reactivity after a power step from 100% to
    `power_frac_final` (0 = shutdown -> the classic iodine pit).

    sigphi0    : sigma_Xe * phi at 100% power [1/s]:
                 2.65e6 b = 2.65e-18 cm2, phi_th ~ 2e13  ->  5.3e-5;
                 PWR band ~ 3e-5 .. 1.2e-4
    rho_eq_pcm : equilibrium-xenon worth at 100% power (calibrates the
                 reactivity scale; CoreForge's 0-D chain gives -3178).
    Returns t_h, rho_pcm, Xe/Xe_eq plus pit depth/time and, for
    shutdown, the time xenon decays back to its pre-shutdown worth.
    Validation: the analytic pit time  t_peak = ln(lI/lXe_eff)/(lI-lXe)
    structure is checked in verify.py (expected ~9-12 h)."""
    f = max(0.0, power_frac_final)
    lam_xe_eff0 = LAM_XE135 + sigphi0                 # burnout at 100%
    # equilibrium (relative units, production normalised to 1 at 100%)
    I0 = Y_I135 / LAM_I135
    Xe0 = (Y_XE_DIR + Y_I135) / lam_xe_eff0
    t = np.linspace(0.0, t_end_h * 3600.0, npts)
    dt = t[1] - t[0]
    I, Xe = I0, Xe0
    lam_xe_eff = LAM_XE135 + sigphi0 * f
    Ih, Xeh = [I0], [Xe0]
    for _ in range(npts - 1):
        # exact linear update over dt (I decoupled; Xe linear given I)
        Inew = (I - Y_I135 * f / LAM_I135) * np.exp(-LAM_I135 * dt) \
            + Y_I135 * f / LAM_I135
        src = Y_XE_DIR * f
        # integrate dXe/dt = src + lI*I(t) - lam_eff*Xe with I linearised
        Imid = 0.5 * (I + Inew)
        Xe = (Xe - (src + LAM_I135 * Imid) / lam_xe_eff) \
            * np.exp(-lam_xe_eff * dt) \
            + (src + LAM_I135 * Imid) / lam_xe_eff
        I = Inew
        Ih.append(I)
        Xeh.append(Xe)
    Xe_rel = np.array(Xeh) / Xe0
    rho = rho_eq_pcm * Xe_rel
    t_h = t / 3600.0
    ipk = int(np.argmax(Xe_rel)) if f < 1.0 else 0
    out = dict(t_h=t_h, rho_pcm=rho, xe_rel=Xe_rel,
               pit_time_h=float(t_h[ipk]),
               pit_rho_pcm=float(rho[ipk]),
               pit_ratio=float(Xe_rel[ipk]))
    if f == 0.0:
        back = np.where((t_h > out["pit_time_h"]) & (Xe_rel <= 1.0))[0]
        out["recover_h"] = float(t_h[back[0]]) if len(back) else None
    return out


def simulate(scenario, t_end=60.0, dt=1e-3, P0_MW=160.0, Lambda=2e-5,
             betas=BETA_I, lams=LAM_I, feedback=None, out_every=None,
             trip=None, fuel_mass_kg=None):
    """Run a point-kinetics transient with optional reactivity feedback
    and a reactor-protection trip (automatic scram).

    scenario : dict describing the EXTERNAL reactivity in pcm:
        {"type": "step", "rho_pcm": 50.0}
        {"type": "ramp", "rho_pcm": 100.0, "t_ramp": 10.0}
        {"type": "scram", "rho_pcm": -4000.0, "t_delay": 1.0}
    feedback : None, or dict with a thermal-reactivity model.
        Single-node linear Doppler (back-compatible default):
          {"alpha_pcm_K": -2.5, "mcp_MJ_K": 6.0, "tau_c": 5.0, "T0": 580.0}
        Add "doppler_mode": "sqrt" for physical resonance-broadening
          Doppler  Δρ = A(√T_f − √T_f0),  A matched to alpha_pcm_K at T_f0.
        Add "alpha_mod_pcm_K" (moderator temperature coefficient, MTC) to
          switch on a TWO-NODE fuel+moderator thermal model:
          {"alpha_pcm_K":-3, "doppler_mode":"sqrt", "mcp_MJ_K":6.0,
           "alpha_mod_pcm_K":-15, "mcp_mod_MJ_K":12.0,
           "tau_fm":6.0, "tau_ms":3.0, "T_sink":565.0}
          Nominal fuel/moderator temperatures are derived from the steady
          heat balance at P0, so Doppler and MTC are both zero at nominal.
    trip : None, or a reactor-protection description that inserts a scram
        bank automatically when a setpoint is crossed:
          {"power_frac": 1.18,      # high-flux trip (× nominal)
           "period_s": 10.0,        # optional short-period trip [s]
           "delay_s": 0.5,          # signal → rod-motion delay
           "scram_rho_pcm": -6000., # total scram-bank worth
           "scram_ramp_s": 2.0}     # rod insertion (drop) time
    fuel_mass_kg : optional fuel mass; if given the peak specific enthalpy
        rise is reported in cal/g (the RIA safety metric).

    Returns dict(t, P_MW, rho_pcm, T_fuel, [T_mod], peak_MW, energy_MJ,
    prompt_critical, t_trip, period_inhour_s ...).
    """
    beta = betas.sum()
    y = equilibrium_y(1.0, Lambda, betas, lams)
    fb = feedback or {}

    # --- thermal / feedback model -------------------------------------
    two_node = ("alpha_mod_pcm_K" in fb) or ("mcp_mod_MJ_K" in fb)
    alpha_dop = fb.get("alpha_pcm_K", 0.0) * 1e-5      # linear Doppler [Δk/K]
    dop_sqrt = fb.get("doppler_mode") == "sqrt"
    Cf = fb.get("mcp_MJ_K", 6.0) * 1e6                 # J/K, fuel
    P0_W = P0_MW * 1e6

    if two_node:
        alpha_mod = fb.get("alpha_mod_pcm_K", 0.0) * 1e-5   # MTC [Δk/K]
        Cm = fb.get("mcp_mod_MJ_K", 12.0) * 1e6             # J/K, moderator
        tau_fm = max(fb.get("tau_fm", 6.0), 1e-6)           # fuel→mod [s]
        tau_ms = max(fb.get("tau_ms", 3.0), 1e-6)           # mod→sink [s]
        R_fm = tau_fm / Cf                                  # K/W
        R_ms = tau_ms / Cm                                  # K/W
        T_sink = fb.get("T_sink", 565.0)
        Tm0 = T_sink + P0_W * R_ms                          # nominal mod T
        Tf0 = Tm0 + P0_W * R_fm                             # nominal fuel T
        Tf, Tm = Tf0, Tm0
    else:
        alpha_mod = 0.0
        tau = max(fb.get("tau_c", 5.0), 1e-6)
        Tf0 = fb.get("T0", 580.0)
        Tf = Tf0
        Tm0 = Tm = Tf0
    # sqrt-Doppler coefficient matched to the linear one at T_f0
    A_dop = alpha_dop * 2.0 * np.sqrt(Tf0) if dop_sqrt else 0.0

    def rho_feedback(Tf, Tm):
        if dop_sqrt:
            rD = A_dop * (np.sqrt(Tf) - np.sqrt(Tf0))
        else:
            rD = alpha_dop * (Tf - Tf0)
        rM = alpha_mod * (Tm - Tm0) if two_node else 0.0
        return rD + rM

    def rho_ext(t):
        s = scenario
        r = s.get("rho_pcm", 0.0) * 1e-5
        if s["type"] == "step":
            return r
        if s["type"] == "ramp":
            tr = max(s.get("t_ramp", 1.0), 1e-9)
            return r * min(t / tr, 1.0)
        if s["type"] == "scram":
            return r if t >= s.get("t_delay", 0.0) else 0.0
        raise ValueError(f"unknown scenario type {s['type']}")

    # --- reactor-protection trip (automatic scram) --------------------
    tp = trip or {}
    trip_power = tp.get("power_frac", None)
    trip_period = tp.get("period_s", None)
    trip_delay = tp.get("delay_s", 0.5)
    trip_worth = tp.get("scram_rho_pcm", -6000.0) * 1e-5
    trip_ramp = max(tp.get("scram_ramp_s", 2.0), 1e-9)
    t_trip = [None]

    def rho_scram(t):
        if t_trip[0] is None:
            return 0.0
        t0 = t_trip[0] + trip_delay
        if t < t0:
            return 0.0
        return trip_worth * min((t - t0) / trip_ramp, 1.0)

    if out_every is None:
        out_every = max(1, int((t_end / dt) / 2000))
    ts, Ps, rhos, Tfs, Tms = [], [], [], [], []
    nstep = int(round(t_end / dt))
    cacheA = {}
    t = 0.0
    n_prev = y[0]
    energy_MJ = 0.0
    prompt_critical = False
    for i in range(nstep + 1):
        rho = rho_ext(t) + rho_scram(t) + rho_feedback(Tf, Tm)
        if rho_ext(t) > beta:
            prompt_critical = True
        if i % out_every == 0:
            ts.append(t)
            Ps.append(y[0] * P0_MW)
            rhos.append(rho * 1e5)
            Tfs.append(Tf)
            Tms.append(Tm)
        if i == nstep:
            break
        # trip logic: latch the first setpoint crossing
        if trip and t_trip[0] is None:
            if trip_power is not None and y[0] >= trip_power:
                t_trip[0] = t
            elif trip_period is not None and t > 0:
                dndt = (y[0] - n_prev) / dt
                if dndt > 0 and y[0] / dndt < trip_period:
                    t_trip[0] = t
        n_prev = y[0]
        # frozen-rho exact step (cache expm on a rounded-rho key)
        key = round(rho * 1e5, 1)
        if key not in cacheA:
            if len(cacheA) > 4000:
                cacheA.clear()
            cacheA[key] = _expm(_pk_matrix(key * 1e-5, Lambda, betas, lams),
                                dt)
        y = cacheA[key] @ y
        y = np.maximum(y, 0.0)
        energy_MJ += y[0] * P0_MW * dt
        # thermal update (only if a feedback model is supplied)
        if feedback:
            P = y[0] * P0_W                        # W
            if two_node:
                dTf = ((P - (Tf - Tm) / R_fm) / Cf) * dt
                dTm = (((Tf - Tm) / R_fm - (Tm - T_sink) / R_ms) / Cm) * dt
                Tf += dTf
                Tm += dTm
            else:
                Tf = Tf + dt * ((P - P0_W) / Cf - (Tf - Tf0) / tau)
        t += dt

    Ps = np.array(Ps)
    Tfs = np.array(Tfs)
    out = dict(t=np.array(ts), P_MW=Ps, rho_pcm=np.array(rhos),
               T_fuel=Tfs, peak_MW=float(Ps.max()),
               final_MW=float(Ps[-1]), final_T=float(Tfs[-1]),
               peak_time_s=float(ts[int(np.argmax(Ps))]),
               energy_MJ=float(energy_MJ),
               rho_dollars_ext=float(scenario.get("rho_pcm", 0.0)
                                     * 1e-5 / beta),
               prompt_critical=bool(prompt_critical),
               t_trip=(None if t_trip[0] is None else float(t_trip[0])),
               beta=beta, Lambda=Lambda)
    if two_node:
        out["T_mod"] = np.array(Tms)
    dTf_peak = float(Tfs.max() - Tfs[0]) if feedback else 0.0
    out["dT_fuel_peak"] = dTf_peak
    if feedback and fuel_mass_kg:
        out["enthalpy_cal_g"] = (Cf * dTf_peak
                                 / (fuel_mass_kg * 1000.0) / 4.184)
    # asymptotic period (meaningful only for a constant final rho: a pure
    # step/scram with NO feedback and NO trip)
    if not feedback and not trip and scenario["type"] in ("step", "scram"):
        w = inhour_omega(scenario.get("rho_pcm", 0.0) * 1e-5,
                         Lambda, betas, lams)
        out["period_inhour_s"] = (1.0 / w) if w != 0 else np.inf
        # measured from the tail of the simulation
        tail = Ps[-max(4, len(Ps) // 10):]
        ttail = out["t"][-len(tail):]
        if np.all(tail > 0) and tail[-1] != tail[0]:
            slope = np.polyfit(ttail, np.log(tail), 1)[0]
            out["period_sim_s"] = (1.0 / slope) if slope != 0 else np.inf
    return out


# ---------------------------------------------------------------------
# Canonical reactivity-accident presets (scenario + feedback + trip)
# ---------------------------------------------------------------------
def accident_preset(name, P0_MW=160.0):
    """Ready-to-run design-basis reactivity accidents.  Returns a kwargs
    dict for simulate() (scenario, feedback, trip, t_end, fuel_mass_kg).

    'rea'  — rod-ejection accident: a ~1.2$ ejected-rod step, √T Doppler +
             moderator feedback, protected by a high-flux trip + scram.
    'atws' — the SAME ejection but the scram FAILS (trip removed): the
             excursion is held only by inherent (Doppler+MTC) feedback,
             settling at an elevated quasi-steady power.
    'rod_withdrawal' — uncontrolled bank withdrawal: a slow +ramp caught
             by the high-flux trip.
    """
    name = name.lower()
    fb = dict(alpha_pcm_K=-3.0, doppler_mode="sqrt", mcp_MJ_K=6.0,
              alpha_mod_pcm_K=-15.0, mcp_mod_MJ_K=12.0,
              tau_fm=6.0, tau_ms=3.0, T_sink=565.0)
    if name in ("rea", "rod_ejection"):
        return dict(scenario={"type": "step", "rho_pcm": 800.0},
                    feedback=fb,
                    trip=dict(power_frac=1.18, delay_s=0.5,
                              scram_rho_pcm=-6000.0, scram_ramp_s=2.0),
                    t_end=15.0, fuel_mass_kg=None)
    if name == "atws":
        return dict(scenario={"type": "step", "rho_pcm": 800.0},
                    feedback=fb, trip=None, t_end=30.0, fuel_mass_kg=None)
    if name in ("rod_withdrawal", "withdrawal"):
        return dict(scenario={"type": "ramp", "rho_pcm": 300.0,
                              "t_ramp": 30.0},
                    feedback=fb,
                    trip=dict(power_frac=1.15, delay_s=0.7,
                              scram_rho_pcm=-6000.0, scram_ramp_s=2.5),
                    t_end=45.0, fuel_mass_kg=None)
    raise ValueError(f"unknown accident preset {name!r}")
