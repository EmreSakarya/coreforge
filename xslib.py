"""CoreForge fuel designer — semi-empirical two-group PWR pin-cell
cross-section generator.

Turns PHYSICAL inputs (U-235 enrichment, soluble boron, moderator
density, pin geometry) into homogenised two-group constants for the
diffusion engine, the way a lattice code feeds a core simulator.

Model (documented honestly):
  * Number densities from real densities & atomic masses.
  * THERMAL group: true 2200 m/s microscopic data per nuclide
    (U-235: sigma_f 582.6 b, sigma_a 680.9 b, nu 2.432; U-238 2.68 b;
    H 0.3326 b; nat-B 759 b), scaled by ONE spectrum factor S_TH that
    lumps Maxwell averaging + in-fuel flux depression.  Relative
    nuclide physics (e.g. boron-to-fuel worth) is therefore preserved.
  * FAST group & slowing-down: four semi-empirical per-atom constants
    (U-235 epithermal+fast fission, U-238 self-shielded resonance
    capture, U-238 fast fission, hydrogen slowing-down power).
  * Transport/diffusion: per-nuclide transport constants -> D = 1/(3*Str).
  * All tuned constants are CALIBRATED ONCE to a nominal fresh PWR cell
    (e = 3.1 w/o, 0 ppm, rho_mod = 0.71 g/cc -> the well-known
    two-group ballpark D=[1.43,0.38], Sa=[0.0105,~0.10], nuSf=[~0.007,
    ~0.145], S12=0.0165) and then NEVER touched: enrichment/boron/
    density response comes from the physics, not from refitting.

This is an educational/parametric tool ("textbook lattice physics"),
not a lattice code: absolute k is representative, trends are physical.
"""
import math

NA = 0.6022140857          # Avogadro x 1e-24  (atoms/b-cm per g/cc/M)

# ---- fixed nuclear data (2200 m/s thermal microscopics, barns) --------
SF5_TH, SA5_TH, NU5 = 582.6, 680.9, 2.432
SA8_TH = 2.68
SAH_TH = 0.3326
SAO_TH = 1.9e-4
SAB_TH = 759.0             # natural boron
NU8 = 2.80                 # fast fission of U-238

# plutonium / poison chain data (2200 m/s table values; nu thermal)
SF9_TH, SA9_TH, NU9 = 747.4, 1017.9, 2.871          # Pu-239
SF0_TH, SA0_TH      = 0.06, 289.5                    # Pu-240 (thermal absorber)
SF1_TH, SA1_TH, NU1 = 1012.3, 1374.9, 2.917          # Pu-241
SAXE_TH = 2.65e6                                     # Xe-135
SASM_TH = 4.014e4                                    # Sm-149
# lumped fission-product PAIR: one absorber pair per fission, covering
# everything except Xe-135/Sm-149 (tracked separately).  The effective
# (thermal, fast/epithermal) pair below is CALIBRATED once so that the
# 0-D reactivity-limited burnup of 3.1 w/o fuel lands at ~40 MWd/kgU
# (literature 28-45) — larger than the classic 45-50 b/pair value
# because strong individual absorbers (Rh-103, Xe-131, Nd-143, ...) and
# their epithermal resonance integrals are folded into the single pair.
SAFP_TH = 90.0
# fast-group effective microscopics for the chain (documented estimates;
# thermal reactions dominate the depletion reactivity in a PWR cell)
SIGF9_FAST, SIGA9_FAST = 10.0, 11.6                  # Pu-239
SIGF0_FAST, SIGA0_FAST = 0.6, 1.6                    # Pu-240
SIGF1_FAST, SIGA1_FAST = 11.0, 12.6                  # Pu-241
SAXE_FAST = 0.0
SASM_FAST = 0.0
SAFP_FAST = 4.0            # calibrated with SAFP_TH (see above)

# nuclide ordering of the depletion vector (homogenised cell densities)
NUCS = ["U235", "U238", "Pu239", "Pu240", "Pu241", "Xe135", "Sm149", "FP"]
M_NUC = dict(U235=235.044, U238=238.051, Pu239=239.052, Pu240=240.054,
             Pu241=241.057, Xe135=134.907, Sm149=148.917, FP=117.0)
NU_TH = dict(U235=NU5, Pu239=NU9, Pu241=NU1)

# ---- calibrated constants (fit once to the nominal cell: e = 3.1 w/o,
#      0 ppm, rho_mod = 0.71  ->  D = [1.43, 0.38], S12 = 0.0165,
#      nuSf = [0.0070, 0.145], k_inf = 1.320; never refit afterwards) ---
S_TH   = 0.42356           # thermal spectrum/self-shielding factor
SIGF5_FAST = 8.7148        # b, U-235 epithermal+fast fission (effective)
SIGA5_FAST = 10.2148       # b, U-235 fast absorption (fission + capture)
SIGA8_FAST = 1.10952       # b, U-238 self-shielded resonance capture
SIGF8_FAST = 0.09000       # b, U-238 fast fission (effective)
C_SD   = 0.52023           # b, per H atom slowing-down (Sigma_1->2)
TR1_U, TR1_H2O = 11.614, 9.061   # b/atom fast transport (U, per-H2O)
TR2_U  = 14.0                    # b/atom thermal transport of U (fixed)
TR2_H2O = 48.5175                # solved from the nominal-cell D2 target
# pure-water cells use measured water diffusion constants instead of the
# cell-calibrated transport fit (the fit is only valid inside fuel cells):
D_WATER = [1.55, 0.16]

# per-nuclide (fast, thermal-2200) absorption & fission microscopics
# (assembled here, after the calibrated fast-group constants exist)
SIG_A = dict(U235=(SIGA5_FAST, SA5_TH),  U238=(SIGA8_FAST, SA8_TH),
             Pu239=(SIGA9_FAST, SA9_TH), Pu240=(SIGA0_FAST, SA0_TH),
             Pu241=(SIGA1_FAST, SA1_TH), Xe135=(SAXE_FAST, SAXE_TH),
             Sm149=(SASM_FAST, SASM_TH), FP=(SAFP_FAST, SAFP_TH))
SIG_F = dict(U235=(SIGF5_FAST, SF5_TH),  U238=(SIGF8_FAST, 0.0),
             Pu239=(SIGF9_FAST, SF9_TH), Pu240=(SIGF0_FAST, SF0_TH),
             Pu241=(SIGF1_FAST, SF1_TH), Xe135=(0.0, 0.0),
             Sm149=(0.0, 0.0), FP=(0.0, 0.0))
NU_FAST = dict(U235=NU5, U238=NU8, Pu239=2.90, Pu240=2.80, Pu241=2.95)

RHO_UO2 = 10.4             # g/cc
M_U5, M_U8, M_O, M_H2O, M_B = 235.044, 238.051, 15.999, 18.015, 10.811


def _densities(enrich_wo, boron_ppm, rho_mod):
    """Number densities [atoms/b-cm] for the fuel pellet and the water."""
    w5 = enrich_wo / 100.0
    M_U = 1.0 / (w5 / M_U5 + (1 - w5) / M_U8)
    M_UO2 = M_U + 2 * M_O
    N_UO2 = RHO_UO2 * NA / M_UO2
    N5 = N_UO2 * (w5 / M_U5) * M_U
    N8 = N_UO2 * ((1 - w5) / M_U8) * M_U
    NO_f = 2 * N_UO2
    N_H2O = rho_mod * NA / M_H2O
    NH = 2 * N_H2O
    NO_w = N_H2O
    NB = rho_mod * (boron_ppm * 1e-6) * NA / M_B
    return dict(N5=N5, N8=N8, NO_f=NO_f, N_H2O=N_H2O, NH=NH, NO_w=NO_w, NB=NB)


def _water_densities(boron_ppm, rho_mod):
    """Water-side number densities only (enrichment-independent)."""
    return _densities(3.1, boron_ppm, rho_mod)    # 3.1 is a dummy: only the
                                                  # NH/N_H2O/NO_w/NB keys of
                                                  # the result are used


def fresh_N(enrich_wo, r_fuel=0.4095, pitch=1.26):
    """Homogenised (cell-averaged) heavy-nuclide vector of fresh UO2."""
    n = _densities(enrich_wo, 0.0, 0.71)
    f = math.pi * r_fuel**2 / pitch**2 if (pitch > 0 and r_fuel > 0) else 0.0
    N = {k: 0.0 for k in NUCS}
    N["U235"], N["U238"] = f * n["N5"], f * n["N8"]
    return N


def cell_from_N(N, boron_ppm=0.0, rho_mod=0.71, r_fuel=0.4095, pitch=1.26,
                name=None, enrich_label=None):
    """Homogenised two-group constants for a pin cell with an ARBITRARY
    heavy-nuclide/poison vector N (cell-averaged densities, atoms/b-cm).
    This is the depletion-time path; fresh fuel reproduces pincell_xs
    exactly.  Set r_fuel=0 (with all-zero N) for a pure water cell."""
    nn = _water_densities(boron_ppm, rho_mod)
    f = math.pi * r_fuel**2 / pitch**2 if (pitch > 0 and r_fuel > 0) else 0.0
    w = 1.0 - f
    heavy = ["U235", "U238", "Pu239", "Pu240", "Pu241"]

    # ---- fast group -----------------------------------------------------
    sa1 = sum(SIG_A[k][0] * N[k] for k in NUCS)
    nf1 = sum(NU_FAST.get(k, 0.0) * SIG_F[k][0] * N[k] for k in heavy)
    s12 = w * C_SD * nn["NH"]
    NHV = sum(N[k] for k in heavy)                 # heavy atoms (cell-avg)
    tr1 = TR1_U * NHV + w * TR1_H2O * nn["N_H2O"]

    # ---- thermal group (real microscopics x spectrum factor) ------------
    NO_fuel = 2.0 * NHV                            # oxide oxygen ~ 2/heavy
    sa2 = S_TH * (sum(SIG_A[k][1] * N[k] for k in NUCS)
                  + SAO_TH * NO_fuel
                  + w * (SAH_TH * nn["NH"] + SAO_TH * nn["NO_w"]
                         + SAB_TH * nn["NB"]))
    nf2 = S_TH * sum(NU_TH.get(k, 0.0) * SIG_F[k][1] * N[k] for k in heavy)
    tr2 = TR2_U * NHV + w * TR2_H2O * nn["N_H2O"]

    if f > 0.0:
        D1, D2 = 1.0 / (3.0 * tr1), 1.0 / (3.0 * tr2)
    else:
        D1, D2 = D_WATER
    dsa2_dppm = S_TH * w * SAB_TH * (rho_mod * 1e-6 * NA / M_B)
    kinf = (nf1 + nf2 * s12 / sa2) / (sa1 + s12) if sa2 > 0 else 0.0
    label = name or "cell"
    return dict(name=label,
                D=[D1, D2], Sa=[sa1, sa2], nuSf=[nf1, nf2],
                chi=[1.0, 0.0], scat=[[0.0, s12], [0.0, 0.0]],
                designer=dict(enrich=enrich_label, ppm=boron_ppm,
                              rho_mod=rho_mod, r_fuel=r_fuel, pitch=pitch,
                              dsa2_dppm=dsa2_dppm, kinf=kinf,
                              N=dict(N)))


def pincell_xs(enrich_wo=3.1, boron_ppm=0.0, rho_mod=0.71,
               r_fuel=0.4095, pitch=1.26, name=None):
    """Homogenised two-group constants for one FRESH PWR pin cell.
    Set r_fuel=0 for a pure (borated) water cell."""
    N = fresh_N(enrich_wo, r_fuel, pitch)
    label = name or (f"UO₂ {enrich_wo:.2f}% · {boron_ppm:.0f} ppm"
                     if r_fuel > 0 else f"Water · {boron_ppm:.0f} ppm")
    m = cell_from_N(N, boron_ppm, rho_mod, r_fuel, pitch,
                    name=label, enrich_label=enrich_wo)
    return m


def kinf_of(mat):
    """Infinite-medium k of a two-group material dict."""
    sa1, sa2 = mat["Sa"]
    nf1, nf2 = mat["nuSf"]
    s12 = mat["scat"][0][1]
    return (nf1 + nf2 * s12 / sa2) / (sa1 + s12) if sa2 > 0 else 0.0


# ----------------------------------------------------------------------
def match_fuel_to_xs(target, rho_mod=0.71, r_fuel=0.4095, pitch=1.26):
    """Inverse designer: find the FRESH fuel (enrichment, soluble boron)
    whose generated two-group set best matches a user-supplied macroscopic
    set `target` (a material dict with 2-group D/Sa/nuSf/scat).

    Sequential exact matching on the two constants the model controls
    directly and monotonically:
      1. enrichment  <- nuSf_2  (thermal fission production)
      2. boron [ppm] <- Sa_2    (thermal absorption)
    The remaining constants (D1, D2, Sa1, nuSf1, S12) are REPORTED as
    residuals — they are properties of the generator's calibration, not
    fit for.  Returns dict(enrich, ppm, material, residuals, note).
    """
    tg_nf2 = target["nuSf"][1]
    tg_sa2 = target["Sa"][1]
    if tg_nf2 <= 0:
        raise ValueError("target has no thermal fission (nuSf2 = 0) — "
                         "not a fuel")

    def nf2_of(e):
        return pincell_xs(e, 0.0, rho_mod, r_fuel, pitch)["nuSf"][1]

    lo, hi = 0.7, 5.0
    if not (nf2_of(lo) <= tg_nf2 <= nf2_of(hi)):
        raise ValueError(f"target nuSf2={tg_nf2:.4g} outside the "
                         f"0.7-5.0 w/o range of the generator")
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if nf2_of(mid) < tg_nf2:
            lo = mid
        else:
            hi = mid
    e_eq = 0.5 * (lo + hi)

    def sa2_of(ppm):
        return pincell_xs(e_eq, ppm, rho_mod, r_fuel, pitch)["Sa"][1]

    note = ""
    if sa2_of(0.0) >= tg_sa2:
        ppm_eq = 0.0
        if sa2_of(0.0) > tg_sa2 * (1.0 + 1e-12):
            note = ("target thermal absorption is below the unborated "
                    "cell — matched at 0 ppm, Sa2 residual remains")
    else:
        plo, phi = 0.0, 5000.0
        if sa2_of(phi) < tg_sa2:
            ppm_eq = phi
            note = "Sa2 match capped at 5000 ppm"
        else:
            for _ in range(60):
                pm = 0.5 * (plo + phi)
                if sa2_of(pm) < tg_sa2:
                    plo = pm
                else:
                    phi = pm
            ppm_eq = 0.5 * (plo + phi)

    m = pincell_xs(e_eq, ppm_eq, rho_mod, r_fuel, pitch,
                   name=(target.get("name", "fuel") + " (eq-fuel)"))
    res = {}
    pairs = [("D1", m["D"][0], target["D"][0]),
             ("D2", m["D"][1], target["D"][1]),
             ("Sa1", m["Sa"][0], target["Sa"][0]),
             ("Sa2", m["Sa"][1], target["Sa"][1]),
             ("nuSf1", m["nuSf"][0], target["nuSf"][0]),
             ("nuSf2", m["nuSf"][1], target["nuSf"][1]),
             ("S12", m["scat"][0][1], target["scat"][0][1])]
    for key, got, want in pairs:
        res[key] = (100.0 * (got - want) / want) if want != 0 else None
    res["kinf_target"] = kinf_of(target)
    res["kinf_eq"] = m["designer"]["kinf"]
    return dict(enrich=e_eq, ppm=ppm_eq, material=m, residuals=res,
                note=note)


def reliability_flag(residuals, threshold=15.0):
    """Reactivity-relevant reliability check on an inverse-designer fit.

    Sa1/Sa2 residuals drive removal/absorption and therefore k directly;
    a large one (e.g. an inserted control rod, whose absorption spectrum
    a boron-water model cannot reproduce) means results computed with
    the equivalent — especially k_eff and reactivity worths — should NOT
    be trusted for THAT material, even though nuSf2/Sa2 look "matched".
    Returns (reliable: bool, worst_key: str, worst_pct: float)."""
    watch = {k: residuals[k] for k in ("Sa1", "Sa2")
            if isinstance(residuals.get(k), float)}
    if not watch:
        return True, None, 0.0
    worst_key = max(watch, key=lambda k: abs(watch[k]))
    worst = watch[worst_key]
    return abs(worst) <= threshold, worst_key, worst
