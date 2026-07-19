"""CoreForge burnup (depletion) module.

Quasi-static fuel-cycle simulation on top of the Fortran flux engine:

    solve flux  ->  normalise to the specified power  ->  integrate the
    nuclide chain over the step  ->  rebuild cross sections  ->  repeat

Chain (per burnable block, cell-homogenised densities):

    U-235  --(a)-->  removed                    (thermal+fast absorption)
    U-238  --(capture)-->  Pu-239               (2 fast betas folded)
    Pu-239 --(capture)-->  Pu-240 --(c)--> Pu-241 --(a)--> removed
    fission --> Xe-135 (equilibrium), Sm-149 (equilibrium), lumped FP

Power normalisation converts the engine's relative flux to absolute
n/cm^2/s from the specified specific power [W/gU], so all reaction
rates and burnup [MWd/kgU] are absolute.

Limits (stated honestly): 2-group frozen spectrum; no Doppler/T-H
feedback; no U-236/Np chain; no discrete burnable absorbers;
equilibrium Xe only; RADIAL (2-D) cores only in this version.
"""
import copy

import numpy as np

import runner
import xslib

KAPPA_J = 3.204e-11          # J per fission (200 MeV)
Y_XE = 0.066                 # cumulative I-135 -> Xe-135 fission yield
Y_SM = 0.0113                # cumulative Pm-149 -> Sm-149 yield
Y_FP = 1.0                   # lumped FP pairs per fission
LAM_XE = 2.093e-5            # Xe-135 decay [1/s]
HEAVY = ["U235", "U238", "Pu239", "Pu240", "Pu241"]


def _u_mass_density(N):
    """Initial-heavy-metal mass density [g/cm^3], cell-homogenised."""
    return sum(N[k] * xslib.M_NUC[k] for k in HEAVY) / xslib.NA


def _micro(kind, nuc):
    s1, s2200 = (xslib.SIG_A if kind == "a" else xslib.SIG_F)[nuc]
    return s1 * 1e-24, xslib.S_TH * s2200 * 1e-24     # [cm^2] eff. pair


def _rates(N, phi1, phi2):
    """Per-nuclide absorption/fission rates [1/s per atom] and the cell
    fission-rate density [fissions/cm^3/s]."""
    ra, rf, F = {}, {}, 0.0
    for k in xslib.NUCS:
        a1, a2 = _micro("a", k)
        f1, f2 = _micro("f", k)
        ra[k] = a1 * phi1 + a2 * phi2
        rf[k] = f1 * phi1 + f2 * phi2
        F += rf[k] * N[k] * 1e24                      # N is atoms/b-cm
    return ra, rf, F


def _capture(ra, rf, k):
    return max(ra[k] - rf[k], 0.0)


def equilibrium_poisons(N, phi1, phi2):
    """Set Xe-135 and Sm-149 to their at-power equilibrium densities."""
    ra, rf, F = _rates(N, phi1, phi2)
    F_bcm = F * 1e-24                                  # fissions/(b-cm)/s
    axe = ra["Xe135"]
    N["Xe135"] = Y_XE * F_bcm / (LAM_XE + axe) if (LAM_XE + axe) > 0 else 0.0
    asm = ra["Sm149"]
    N["Sm149"] = Y_SM * F_bcm / asm if asm > 0 else 0.0
    return N


def deplete_step(N, phi1, phi2, dt_s, nsub=40):
    """Integrate the heavy chain + lumped FP over dt_s [s] (RK4 on the
    small ODE system, poisons re-equilibrated at the end)."""
    y = np.array([N[k] for k in HEAVY] + [N["FP"]])

    ra, rf, _ = _rates(N, phi1, phi2)      # rates fixed over the step

    def f(v):
        n5, n8, n9, n0, n1, nfp = v
        c8 = _capture(ra, rf, "U238")
        c9 = _capture(ra, rf, "Pu239")
        c0 = _capture(ra, rf, "Pu240")
        d5 = -ra["U235"] * n5
        d8 = -ra["U238"] * n8
        d9 = c8 * n8 - ra["Pu239"] * n9
        d0 = c9 * n9 - ra["Pu240"] * n0
        d1 = c0 * n0 - ra["Pu241"] * n1
        Fb = (rf["U235"] * n5 + rf["U238"] * n8 + rf["Pu239"] * n9 +
              rf["Pu240"] * n0 + rf["Pu241"] * n1)
        dfp = Y_FP * Fb - ra["FP"] * nfp
        return np.array([d5, d8, d9, d0, d1, dfp])

    h = dt_s / nsub
    for _ in range(nsub):
        k1 = f(y)
        k2 = f(y + 0.5 * h * k1)
        k3 = f(y + 0.5 * h * k2)
        k4 = f(y + h * k3)
        y = y + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        y = np.maximum(y, 0.0)

    for i, k in enumerate(HEAVY):
        N[k] = float(y[i])
    N["FP"] = float(y[5])
    return equilibrium_poisons(N, phi1, phi2)


# ----------------------------------------------------------------------
def _block_flux(res, cfg):
    """Block-averaged (phi1, phi2) from an engine result, top-first."""
    div = cfg["div"]
    nby, nbx = np.asarray(cfg["blocks"]).shape
    out = []
    for g in range(cfg["ng"]):
        p = np.flipud(res["phi"][g])                   # top-first
        out.append(p.reshape(nby, div, nbx, div).mean(axis=(1, 3)))
    return out[0], out[1]


def reload_core(prev_state, fresh_fraction=0.34):
    """Multi-cycle reload: return an init_state for the NEXT cycle in
    which the most-burned `fresh_fraction` of the fuel blocks are
    replaced by FRESH fuel (their original composition), while the rest
    carry their end-of-cycle nuclide vectors and cumulative burnup —
    the standard batch-reload strategy of a core designer."""
    items = sorted(prev_state.items(), key=lambda kv: -kv[1]["bu"])
    n_fresh = max(1, int(round(fresh_fraction * len(items))))
    init = {}
    for rank, (key, s) in enumerate(items):
        if rank < n_fresh:
            init[key] = dict(N=dict(s["N0"]), bu=0.0)      # fresh batch
        else:
            N = dict(s["N"])                               # carry over
            # the outage decays Xe-135/Sm-149 away; the next cycle starts
            # xenon-free and rebuilds equilibrium, so zero them here (keeps
            # the BOC 'no-xenon' reference and its xenon worth physical)
            N["Xe135"] = 0.0
            N["Sm149"] = 0.0
            init[key] = dict(N=N, bu=s["bu"])
    return init


def deplete_core(cfg, specific_power=38.0, bu_target=36.0, dbu=1.5,
                 letdown=False, ppm_fixed=None, threads=None,
                 callback=None, tol_pcm=10.0, k_eoc=1.0, bu_cap=60.0,
                 init_state=None):
    """Block-wise core depletion.

    cfg            : a core whose burnable materials are designer-made
    specific_power : core-average specific power [W/gU]
    bu_target      : end-of-cycle core-average burnup [MWd/kgU], or None
                     for AUTO end-of-cycle: march until the core can no
                     longer hold k_eff >= k_eoc (with letdown: until the
                     critical boron reaches 0 ppm and k falls below k_eoc)
                     — the physics then determines the cycle length.
    dbu            : step size [MWd/kgU]
    letdown        : if True, find critical boron each step (letdown curve)
    ppm_fixed      : operating boron when letdown=False (default: as built)
    k_eoc          : end-of-cycle criterion for the auto mode
    bu_cap         : hard stop for the auto mode [MWd/kgU]

    Returns dict(history=..., bumap=..., maps=..., state=..., efpd=...).
    """
    if cfg.get("axial"):
        raise ValueError("burnup is 2-D (radial) only in this version — "
                         "run the radial model of the core")
    cfg = copy.deepcopy(cfg)
    blocks = np.asarray(cfg["blocks"])
    nby, nbx = blocks.shape

    mats = cfg["materials"]
    burnable_ids = [i + 1 for i, m in enumerate(mats)
                    if "designer" in m and
                    m["designer"].get("N", {}).get("U238", 0.0) > 0.0]
    if not burnable_ids:
        raise ValueError("no designer fuel materials to deplete")
    water_ids = [i + 1 for i, m in enumerate(mats) if i + 1 not in burnable_ids]

    state = {}
    for j in range(nby):
        for i in range(nbx):
            mid = int(blocks[j, i])
            if mid in burnable_ids:
                d = mats[mid - 1]["designer"]
                st0 = (init_state or {}).get((j, i), {})
                state[(j, i)] = dict(N=dict(st0.get("N", d["N"])),
                                     bu=float(st0.get("bu", 0.0)),
                                     N0=dict(d["N"]),   # fresh, for reload
                                     rho_mod=d["rho_mod"], r_fuel=d["r_fuel"],
                                     pitch=d["pitch"])
    rhoU = {k: _u_mass_density(s["N0"]) for k, s in state.items()}

    ppm0 = ppm_fixed
    if ppm0 is None:
        ppm0 = next(m["designer"]["ppm"] for m in mats if "designer" in m)

    def build_cfg(ppm):
        c = copy.deepcopy(cfg)
        newmats = [copy.deepcopy(mats[i - 1]) for i in water_ids]
        for m in newmats:
            if "designer" in m:
                d = m["designer"]
                mm = xslib.cell_from_N(d["N"], ppm, d["rho_mod"],
                                       d["r_fuel"], d["pitch"])
                m.update({k: mm[k] for k in
                          ("D", "Sa", "nuSf", "chi", "scat")})
        idmap = {mid: k + 1 for k, mid in enumerate(water_ids)}
        newblocks = [[0] * nbx for _ in range(nby)]
        for (j, i), s in state.items():
            mm = xslib.cell_from_N(s["N"], ppm, s["rho_mod"],
                                   s["r_fuel"], s["pitch"],
                                   name=f"b{j}-{i}")
            newmats.append(mm)
            newblocks[j][i] = len(newmats)
        for j in range(nby):
            for i in range(nbx):
                if (j, i) not in state:
                    newblocks[j][i] = idmap[int(blocks[j, i])]
        c["materials"] = newmats
        c["blocks"] = newblocks
        return c

    def solve_at(ppm):
        return runner.run_case(build_cfg(ppm), threads=threads)

    def find_critical_ppm(start):
        r0 = solve_at(0.0)
        if r0["keff"] <= 1.0 + tol_pcm * 1e-5:
            return 0.0
        lo, hi = 0.0, max(1000.0, start)
        khi = solve_at(hi)["keff"]
        while khi > 1.0 and hi < 40000.0:
            lo, hi = hi, hi * 2.0
            khi = solve_at(hi)["keff"]
        mid = hi
        for _ in range(14):
            mid = 0.5 * (lo + hi)
            km = solve_at(mid)["keff"]
            if abs(runner.rho_pcm(km)) < tol_pcm:
                return mid
            if km > 1.0:
                lo = mid
            else:
                hi = mid
        return mid

    def block_scales(r, ppm):
        """Absolute per-block flux scale so that block specific power
        matches SP * (relative block power): F_abs = SP_b*rhoU/kappa."""
        ph1, ph2 = _block_flux(r, cfg)
        div = cfg["div"]
        pfine = np.flipud(r["power"])                       # top-first fine
        pmap = pfine.reshape(nby, div, nbx, div).mean(axis=(1, 3))
        pavg = np.mean([pmap[j, i] for (j, i) in state])
        sc = {}
        for (j, i), s in state.items():
            sp_b = specific_power * pmap[j, i] / pavg
            F_abs = sp_b * rhoU[(j, i)] / KAPPA_J          # fis/cm^3/s
            _, rfv, _ = _rates(s["N"], ph1[j, i], ph2[j, i])
            Frel = sum(rfv[k] * s["N"][k] for k in xslib.NUCS) * 1e24
            sc[(j, i)] = F_abs / Frel if Frel > 0 else 0.0
        return sc, ph1, ph2, pmap, pavg

    auto = bu_target is None
    bu_end = bu_cap if auto else bu_target
    history = []
    bu = 0.0
    efpd = 0.0
    ppm = ppm0
    maps = {}
    k_noxe = None
    nstep = 0
    while True:
        # ---- flux at this burnup point ---------------------------------
        r = solve_at(ppm)
        if nstep == 0:
            k_noxe = r["keff"]                     # BOC, no xenon yet

        # ---- equilibrate Xe/Sm at the current flux ----------------------
        sc, ph1, ph2, _, _ = block_scales(r, ppm)
        for (j, i), s in state.items():
            equilibrium_poisons(s["N"], sc[(j, i)] * ph1[j, i],
                                sc[(j, i)] * ph2[j, i])
        if nstep == 0:
            # xenon worth at unchanged boron (before any letdown)
            k_xe = solve_at(ppm)["keff"]
            maps["xenon_worth_pcm"] = (runner.rho_pcm(k_xe) -
                                       runner.rho_pcm(k_noxe))
        if letdown:
            ppm = find_critical_ppm(ppm if ppm > 0 else 1000.0)

        # re-solve with equilibrium poisons (and letdown boron) in place
        r = solve_at(ppm)
        sc, ph1, ph2, pmap, pavg = block_scales(r, ppm)

        hrow = dict(bu=bu, efpd=efpd, keff=r["keff"], ppm=ppm, fxy=r["fxy"],
                    rho_pcm=runner.rho_pcm(r["keff"]))
        history.append(hrow)
        if callback:
            callback(nstep, hrow)
        if nstep == 0:
            maps["power_boc"] = r["power"]
            maps["res_boc"] = r

        # ---- end-of-cycle tests ----------------------------------------
        eoc = False
        if auto:
            depleted_out = (r["keff"] < k_eoc - tol_pcm * 1e-5 and
                            (not letdown or ppm <= 1e-9))
            eoc = depleted_out or bu >= bu_end - 1e-9
            if eoc:
                maps["eoc_reason"] = ("reactivity-limited"
                                      if depleted_out else
                                      f"bu_cap {bu_cap} reached")
        else:
            eoc = bu >= bu_end - 1e-9
        if eoc:
            maps["power_eol"] = r["power"]
            maps["res_eol"] = r
            break

        # ---- deplete one step ------------------------------------------
        step = min(dbu, bu_end - bu)
        dt_days = 1000.0 * step / specific_power
        dt_s = dt_days * 86400.0
        for (j, i), s in state.items():
            sp_b = specific_power * pmap[j, i] / pavg
            deplete_step(s["N"], sc[(j, i)] * ph1[j, i],
                         sc[(j, i)] * ph2[j, i], dt_s)
            s["bu"] += sp_b * dt_days / 1000.0
        bu += step
        efpd += dt_days
        nstep += 1

    bumap = np.full((nby, nbx), np.nan)
    for (j, i), s in state.items():
        bumap[j, i] = s["bu"]
    return dict(history=history, bumap=bumap, maps=maps, state=state,
                ppm_final=ppm, efpd=efpd,
                eoc_reason=maps.get("eoc_reason"))


# ----------------------------------------------------------------------
def deplete_cell_kinf(enrich=3.1, specific_power=38.0, bu_target=45.0,
                      dbu=1.5, ppm=0.0):
    """0-D (infinite-medium) depletion of one pin cell — engine-free.
    Returns list of dicts (bu, kinf, N...) — used by verify.py."""
    N = xslib.fresh_N(enrich)
    rhoU = _u_mass_density(N)
    out = []
    bu = 0.0
    while True:
        m = xslib.cell_from_N(N, ppm)
        sa2 = m["Sa"][1]
        s12 = m["scat"][0][1]
        ratio = s12 / sa2                          # phi2/phi1 (inf. medium)
        F_abs = specific_power * rhoU / KAPPA_J
        _, rfv, _ = _rates(N, 1.0, ratio)
        Frel = sum(rfv[k] * N[k] for k in xslib.NUCS) * 1e24
        sc = F_abs / Frel if Frel > 0 else 0.0
        equilibrium_poisons(N, sc, sc * ratio)
        m = xslib.cell_from_N(N, ppm)
        out.append(dict(bu=bu, kinf=m["designer"]["kinf"],
                        pu239=N["Pu239"], u235=N["U235"]))
        if bu >= bu_target - 1e-9:
            break
        step = min(dbu, bu_target - bu)
        dt_s = 1000.0 * step / specific_power * 86400.0
        deplete_step(N, sc, sc * ratio, dt_s)
        bu += step
    return out
