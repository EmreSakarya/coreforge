"""CoreForge transient module — point reactor kinetics.

Time-dependent core power for reactivity events (rod ejection, rod drop,
ramps), with optional Doppler fuel-temperature feedback:

    dn/dt  = (rho(t) - beta)/Lambda * n + sum_i lambda_i c_i
    dc_i/dt = beta_i/Lambda * n - lambda_i c_i
    M_cp dT/dt = (P - P0) - M_cp (T - T0)/tau_c        [feedback model]
    rho(t) = rho_ext(t) + alpha_D (T - T0)

Numerics: over each output step the kinetics matrix is FROZEN and the
exact solution y(t+dt) = expm(A dt) y(t) is applied via eigen-
decomposition (numpy only) — unconditionally stable, no stiffness
issues from the prompt mode.  The temperature ODE is advanced with the
same step (its timescale is seconds, the step is milliseconds).

Validation (see verify.py): for constant rho the simulated asymptotic
period must match the analytic INHOUR equation root.

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
             betas=BETA_I, lams=LAM_I, feedback=None, out_every=None):
    """Run a point-kinetics transient.

    scenario : dict describing the EXTERNAL reactivity in pcm:
        {"type": "step", "rho_pcm": 50.0}
        {"type": "ramp", "rho_pcm": 100.0, "t_ramp": 10.0}
        {"type": "scram", "rho_pcm": -4000.0, "t_delay": 1.0}
    feedback : None, or dict with Doppler fuel model:
        {"alpha_pcm_K": -2.5, "mcp_MJ_K": 6.0, "tau_c": 5.0, "T0": 580.0}
    Returns dict(t, P_MW, rho_pcm, T_fuel, peak_MW, period_inhour_s ...)
    """
    beta = betas.sum()
    n = 1.0
    y = equilibrium_y(1.0, Lambda, betas, lams)
    fb = feedback or {}
    T = fb.get("T0", 580.0)
    T0 = T
    alpha = fb.get("alpha_pcm_K", 0.0) * 1e-5
    mcp = fb.get("mcp_MJ_K", 6.0) * 1e6            # J/K
    tau = fb.get("tau_c", 5.0)

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

    if out_every is None:
        out_every = max(1, int((t_end / dt) / 2000))
    ts, Ps, rhos, Ts = [], [], [], []
    nstep = int(round(t_end / dt))
    cacheA = {}
    t = 0.0
    for i in range(nstep + 1):
        rho = rho_ext(t) + alpha * (T - T0)
        if i % out_every == 0:
            ts.append(t)
            Ps.append(y[0] * P0_MW)
            rhos.append(rho * 1e5)
            Ts.append(T)
        if i == nstep:
            break
        # frozen-rho exact step (cache expm on a rounded-rho key)
        key = round(rho * 1e5, 1)
        if key not in cacheA:
            if len(cacheA) > 4000:
                cacheA.clear()
            cacheA[key] = _expm(_pk_matrix(key * 1e-5, Lambda, betas, lams),
                                dt)
        y = cacheA[key] @ y
        y = np.maximum(y, 0.0)
        # fuel temperature (only if feedback enabled)
        if feedback:
            P = y[0] * P0_MW * 1e6                 # W
            T = T + dt * ((P - P0_MW * 1e6) / mcp - (T - T0) / tau)
        t += dt

    Ps = np.array(Ps)
    out = dict(t=np.array(ts), P_MW=Ps, rho_pcm=np.array(rhos),
               T_fuel=np.array(Ts), peak_MW=float(Ps.max()),
               final_MW=float(Ps[-1]), final_T=float(Ts[-1]),
               beta=beta, Lambda=Lambda)
    # asymptotic period (only meaningful for constant final rho, no fb)
    if not feedback and scenario["type"] in ("step", "scram"):
        w = inhour_omega(scenario.get("rho_pcm", 0.0) * 1e-5
                         if scenario["type"] == "step"
                         else scenario.get("rho_pcm", 0.0) * 1e-5,
                         Lambda, betas, lams)
        out["period_inhour_s"] = (1.0 / w) if w != 0 else np.inf
        # measured from the tail of the simulation
        tail = Ps[-max(4, len(Ps) // 10):]
        ttail = out["t"][-len(tail):]
        if np.all(tail > 0) and tail[-1] != tail[0]:
            slope = np.polyfit(ttail, np.log(tail), 1)[0]
            out["period_sim_s"] = (1.0 / slope) if slope != 0 else np.inf
    return out
