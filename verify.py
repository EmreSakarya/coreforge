#!/usr/bin/env python3
"""CoreForge — headless verification suite.

Runs the bundled benchmark presets through the compiled engine and checks
the eigenvalues against their references.  Run it right after building:

    python3 verify.py              # quick
    python3 verify.py --fine       # adds the fine-mesh cases
    python3 verify.py --no-engine  # pure-Python physics subset (no Fortran)
"""
import math
import sys

import numpy as np

import burnup
import kinetics
import presets
import runner
import xslib


def bare_square_1g():
    """1-group bare square, zero-flux BC: k = nuSf/(Sa + D*B^2) analytic."""
    return dict(ng=1, pitch=10.0, div=10,
                blocks=[[1] * 12 for _ in range(12)],
                materials=[dict(name="bare", D=[1.0], Sa=[0.07],
                                nuSf=[0.09], chi=[1.0], scat=[[0.0]])],
                bc=["vacuum"] * 4, gamma=1e9,
                ref_keff=0.09 / (0.07 + 2.0 * (math.pi / 120.0) ** 2),
                ref_source="analytic buckling formula")


def bare_cube_1g():
    """3-D: 1-group bare cube, zero-flux BC, B^2 = 3*(pi/L)^2 analytic."""
    return dict(ng=1, pitch=10.0, div=2,
                blocks=[[1] * 12 for _ in range(12)],
                materials=[dict(name="bare", D=[1.0], Sa=[0.07],
                                nuSf=[0.09], chi=[1.0], scat=[[0.0]])],
                bc=["vacuum"] * 6, gamma=1e9,
                axial=dict(dz=10.0, divz=2,
                           zones=[dict(label="cube", layers=12,
                                       blocks=[[1] * 12 for _ in range(12)])]),
                ref_keff=0.09 / (0.07 + 3.0 * (math.pi / 120.0) ** 2),
                ref_source="analytic buckling formula (3-D)")


def equiv_2d_3d():
    """Invariance: a 2-D core and the same core extruded into ONE axial
    zone with reflective bottom/top must give the same k_eff — this
    exercises the whole 3-D input/solve/parse path against the 2-D one."""
    cfg2 = presets.preset_iaea2d()
    cfg2["div"] = 2
    r2 = runner.run_case(cfg2)
    cfg3 = presets.preset_iaea2d()
    cfg3["div"] = 2
    cfg3["bc"] = cfg3["bc"] + ["reflective", "reflective"]
    cfg3["axial"] = dict(dz=20.0, divz=1,
                         zones=[dict(label="extruded", layers=4,
                                     blocks=cfg3["blocks"])])
    r3 = runner.run_case(cfg3)
    d = abs(runner.rho_pcm(r3["keff"]) - runner.rho_pcm(r2["keff"]))
    ok = d < 0.5
    print(f"{'2D == extruded-3D (reflective B/T)':<40}"
          f"{r2['keff']:>12.7f} vs {r3['keff']:.7f}  "
          f"{'PASS' if ok else 'FAIL'}")
    return ok


def designer_trends():
    """Physics trends of the fuel-designer XS generator (no engine)."""
    ks = [xslib.pincell_xs(e, 0)["designer"]["kinf"] for e in (2.1, 3.1, 4.5)]
    mono = ks[0] < ks[1] < ks[2]
    k0 = xslib.pincell_xs(3.1, 590)["designer"]["kinf"]
    k1 = xslib.pincell_xs(3.1, 610)["designer"]["kinf"]
    bw = ((k1 - 1) / k1 - (k0 - 1) / k0) * 1e5 / 20.0
    okbw = -14.0 <= bw <= -4.0
    print(f"{'designer: kinf(e) rising':<40}"
          f"{' / '.join(f'{k:.4f}' for k in ks):>28}  "
          f"{'PASS' if mono else 'FAIL'}")
    print(f"{'designer: boron worth in [-14,-4]':<40}"
          f"{bw:>22.2f} pcm/ppm  {'PASS' if okbw else 'FAIL'}")
    return mono and okbw


def designer_vs_iaea():
    """Independent-source cross-check: the xslib nominal cell against the
    IAEA-2D benchmark fuel constants (1976, completely unrelated origin)."""
    m = xslib.pincell_xs(3.1, 0)
    ref = dict(D1=1.5, D2=0.4, Sa1=0.010, Sa2=0.080, S12=0.020)
    got = dict(D1=m["D"][0], D2=m["D"][1], Sa1=m["Sa"][0], Sa2=m["Sa"][1],
               S12=m["scat"][0][1])
    devs = {k: abs(got[k] - ref[k]) / ref[k] * 100.0 for k in ref}
    worst = max(devs, key=devs.get)
    ok = devs[worst] <= 25.0
    print(f"{'designer vs IAEA-2D fuel (indep. src)':<40}"
          f"{'worst ' + worst + f' {devs[worst]:.1f}%':>28}  "
          f"{'PASS' if ok else 'FAIL'}")
    return ok


def inverse_designer_check():
    """Reconstruct an equivalent fuel from the IAEA-2D fuel-1 macroscopic
    set: the matched constants must be exact, the enrichment plausible."""
    tgt = dict(name="IAEA fuel 1", D=[1.5, 0.4], Sa=[0.01012, 0.080032],
               nuSf=[0.0, 0.135], chi=[1.0, 0.0],
               scat=[[0.0, 0.02], [0.0, 0.0]])
    fit = xslib.match_fuel_to_xs(tgt)
    ok = (2.0 <= fit["enrich"] <= 3.5 and
          abs(fit["residuals"]["nuSf2"]) < 0.1 and
          abs(fit["residuals"]["Sa2"]) < 0.1)
    print(f"{'inverse designer: IAEA fuel-1 -> e_eq':<40}"
          f"{fit['enrich']:>14.3f} w/o, {fit['ppm']:>4.0f} ppm  "
          f"{'PASS' if ok else 'FAIL'}")
    return ok


def burnup_trends():
    """0-D depletion physics of the chain (engine-free, ~2 s)."""
    h = burnup.deplete_cell_kinf(3.1, 38.0, 55.0, 2.5)
    bus = [r["bu"] for r in h]
    ks = [r["kinf"] for r in h]
    bu1 = float(np.interp(1.0, ks[::-1], bus[::-1]))
    ok1 = 28.0 <= bu1 <= 46.0
    print(f"{'burnup: reactivity-limited BU (3.1%)':<40}"
          f"{bu1:>18.1f} MWd/kgU  {'PASS' if ok1 else 'FAIL'}")
    kfresh = xslib.pincell_xs(3.1, 0)["designer"]["kinf"]
    xw = (1.0 / kfresh - 1.0 / ks[0]) * 1e5      # rho(Xe) - rho(fresh) < 0
    ok2 = -4500.0 <= xw <= -2000.0
    print(f"{'burnup: Xe+Sm equilibrium worth':<40}"
          f"{xw:>22.0f} pcm  {'PASS' if ok2 else 'FAIL'}")
    ok3 = h[-1]["pu239"] > 5e-5 and h[-1]["u235"] < 0.4 * h[0]["u235"]
    print(f"{'burnup: U235 burns, Pu239 builds':<40}"
          f"{'':>26}  {'PASS' if ok3 else 'FAIL'}")
    return ok1 and ok2 and ok3


def eqfuel_reliability_check():
    """The reliability flag must fire on a strongly-absorbing material
    (SMR control-rod fuel) and stay silent on plain fuel; the flagged
    case's k_eff sensitivity must be large (reproduces the +2446 pcm
    finding), confirming the warning is not just noise."""
    import copy

    import xslib
    smr = presets.preset_smr()
    plain = xslib.match_fuel_to_xs(smr["materials"][1])   # 3.1%, no rod
    cra = xslib.match_fuel_to_xs(smr["materials"][3])      # 3.1% + CRA
    ok_plain, *_ = xslib.reliability_flag(plain["residuals"])
    ok_cra, wkey, wval = xslib.reliability_flag(cra["residuals"])

    c1 = copy.deepcopy(smr)
    c1["materials"][3] = dict(cra["material"], name="cra-eq")
    k_shift = runner.rho_pcm(runner.run_case(c1)["keff"]) - \
        runner.rho_pcm(runner.run_case(smr)["keff"])

    ok = ok_plain and not ok_cra and abs(k_shift) > 1000.0
    print(f"{'eq-fuel reliability: flags CRA, not plain fuel':<52}"
          f"{wkey}={wval:+.1f}%, Δk={k_shift:+.0f}pcm  "
          f"{'PASS' if ok else 'FAIL'}")
    return ok


def integrity_guard_check():
    """QA fingerprint: mesh refinement must NOT invalidate a benchmark
    reference; any physics change (pitch, XS, BC, axial zones) MUST."""
    import copy
    base = presets.preset_iaea3d()
    fp0 = presets.fingerprint_of(base)
    same = presets.fingerprint_of(copy.deepcopy(base))
    c_div = copy.deepcopy(base)
    c_div["div"] = 7
    c_div["axial"]["divz"] = 3
    c_pitch = copy.deepcopy(base)
    c_pitch["pitch"] = 12.0
    c_xs = copy.deepcopy(base)
    c_xs["materials"][0]["Sa"][1] += 1e-3
    c_bc = copy.deepcopy(base)
    c_bc["bc"][1] = "reflective"
    ok = (same == fp0 and
          presets.fingerprint_of(c_div) == fp0 and
          presets.fingerprint_of(c_pitch) != fp0 and
          presets.fingerprint_of(c_xs) != fp0 and
          presets.fingerprint_of(c_bc) != fp0)
    print(f"{'QA fingerprint: mesh-invariant, physics-sensitive':<52}"
          f"{fp0:>16}  {'PASS' if ok else 'FAIL'}")
    return ok


def flux_calibration_check():
    """Absolute-flux scale S = q_avg/kappa (Operating point panel & HTML
    report): reproducing a stated core power from the calibrated flux
    must close the energy balance exactly (1-group test, nu=1 so nuSf=Sf
    directly, isolating the calibration formula from nu-bar ambiguity)."""
    cfg = dict(ng=1, pitch=10.0, div=6, blocks=[[1] * 10 for _ in range(10)],
               materials=[dict(name="test", D=[1.2], Sa=[0.06],
                               nuSf=[0.08], chi=[1.0], scat=[[0.0]])],
               bc=["reflective"] * 4, gamma=0.4692)
    r = runner.run_case(cfg)
    KAPPA = 3.204e-11
    P_MW = 50.0
    fuelmask = r["power"] > 1e-12
    Vf = int(fuelmask.sum()) * r["dx"] ** 2 * 200.0
    q_avg = P_MW * 1e6 / Vf
    S = q_avg / KAPPA                          # the formula under test
    Sf = 0.08                                  # = nuSf since nu=1 here
    implied_MW = (Sf * np.sum(r["phi"][0][fuelmask] * S)
                 * r["dx"] ** 2 * 200.0) * KAPPA / 1e6
    ok = abs(implied_MW / P_MW - 1.0) < 0.01
    print(f"{'flux calibration: energy closure':<40}"
          f"{implied_MW:>10.3f} MW vs {P_MW:.1f} MW target  "
          f"{'PASS' if ok else 'FAIL'}")
    return ok


def thermal_checks():
    """Closed-channel T-H: exact energy closure + hand-verifiable
    conduction/convection formulas (engine-free)."""
    import math

    import thermal
    one = thermal.channel(4.32, 2.0, mflow_kg_s=70.0, t_in_C=258.0)
    Q = 70.0 * thermal.CP_W * (one["T_out"] - 258.0) / 1e6
    ok1 = abs(Q - 4.32) < 1e-9
    dTf = 20000.0 / (4 * math.pi * thermal.K_FUEL)
    ok2 = abs(dTf - 530.5) < 1.0
    ok3 = one["Re"] > 1e4 and one["htc"] > 1e4
    # T_sat(P) correlation + fz-clamp flag (audit fixes)
    ok4 = (abs(thermal.t_sat(15.5) - 344.8) < 0.5 and
           abs(thermal.t_sat(7.0) - 285.8) < 1.5)
    hi = thermal.channel(4.32, 2.0, mflow_kg_s=70.0, fz=2.2)
    lo = thermal.channel(4.32, 2.0, mflow_kg_s=70.0, fz=1.55)
    ok5 = hi["fz_clamped"] and not lo["fz_clamped"]
    ok = ok1 and ok2 and ok3 and ok4 and ok5
    print(f"{'thermal: energy + formulas + Tsat(P) + fz-flag':<52}"
          f"Q={Q:.4f}, Tsat={thermal.t_sat(15.5):.1f}  "
          f"{'PASS' if ok else 'FAIL'}")
    return ok


def xenon_checks():
    """Iodine-pit physics window (0-D, engine-free)."""
    s = kinetics.xenon_transient(0.0, 72.0)
    ok = (6.0 <= s["pit_time_h"] <= 13.0 and
          1.2 <= s["pit_ratio"] <= 4.0 and
          s.get("recover_h") is not None)
    print(f"{'xenon: iodine pit time/depth in PWR band':<52}"
          f"t={s['pit_time_h']:.1f} h, ×{s['pit_ratio']:.2f}  "
          f"{'PASS' if ok else 'FAIL'}")
    return ok


def multicycle_check():
    """Batch reload: a partially-burned reloaded core must be less
    reactive at BOC than the fresh first core (engine, ~10 solves)."""
    cfg = presets.preset_designer(1200.0)
    c1 = burnup.deplete_core(cfg, bu_target=4.0, dbu=2.0, letdown=False)
    init = burnup.reload_core(c1["state"], 0.34)
    c2 = burnup.deplete_core(cfg, bu_target=4.0, dbu=2.0, letdown=False,
                             init_state=init)
    k1, k2 = c1["history"][0]["keff"], c2["history"][0]["keff"]
    # reload zeros carried Xe/Sm, so cycle-2 xenon worth is now physical
    # (was understated ~22% before the audit fix)
    xw2 = c2["maps"]["xenon_worth_pcm"]
    ok = (k2 < k1 and c2["history"][-1]["bu"] == 4.0 and xw2 < -3000.0)
    print(f"{'multi-cycle: reloaded BOC + physical Xe worth':<52}"
          f"k1={k1:.5f}>k2={k2:.5f}, Xe2={xw2:.0f}  "
          f"{'PASS' if ok else 'FAIL'}")
    return ok


def rod_position_check():
    """3-D critical-rod-position search hits a mid-range target k
    (uses the IAEA-3D rod-5 depth parametrisation, coarse mesh)."""
    def k_at(depth):
        c = presets.preset_iaea3d()
        c["axial"] = presets.iaea3d_axial(float(depth), dz=10.0)
        return runner.run_case(c)["keff"]

    k0, k1 = k_at(0.0), k_at(340.0)
    tgt = 0.5 * (k0 + k1)
    lo, hi = 0.0, 340.0
    d = 0.0
    for _ in range(10):
        d = 0.5 * (lo + hi)
        km = k_at(d)
        if abs(runner.rho_pcm(km) - runner.rho_pcm(tgt)) < 3.0:
            break
        if km > tgt:
            lo = d
        else:
            hi = d
    # depth is quantised to whole dz layers (~29 pcm per 10 cm step for
    # this rod) — the search must land within one layer of the target
    ok = k1 < tgt < k0 and abs(runner.rho_pcm(km) -
                               runner.rho_pcm(tgt)) < 25.0
    print(f"{'rod-position search: hits mid-range target':<52}"
          f"depth={d:.0f} cm, k={km:.5f}  {'PASS' if ok else 'FAIL'}")
    return ok


def kinetics_trends():
    """Point-kinetics vs the analytic inhour equation (engine-free)."""
    ok = True
    for rho, t_end, tol in ((+100.0, 120.0, 1.0), (-100.0, 900.0, 1.0)):
        s = kinetics.simulate({"type": "step", "rho_pcm": rho},
                              t_end=t_end, dt=2e-3)
        Ti, Ts = s["period_inhour_s"], s["period_sim_s"]
        err = 100.0 * abs(Ts - Ti) / abs(Ti)
        good = err <= tol
        ok &= good
        print(f"{'kinetics: inhour period, ' + f'{rho:+.0f} pcm':<40}"
              f"{Ti:>10.1f} vs {Ts:>8.1f} s ({err:.2f}%)  "
              f"{'PASS' if good else 'FAIL'}")
    return ok


def version_consistency_check():
    """Version-string lockstep guard (recurring release-drift bug).
    The engine header + --version banner + run banner, report.py VERSION,
    and app.py APP_VERSION + footer must ALL declare the same version.
    Historically these drifted silently (the footer showed one version
    while the engine reported another); this makes any mismatch a hard
    failure.  Engine-free — just parses the source files."""
    import re

    def grab(path, pat, label):
        with open(path, encoding="utf-8") as f:
            hits = re.findall(pat, f.read())
        return [(f"{label}[{i}]", v) for i, v in enumerate(hits)]

    sites = []
    sites += grab("app.py", r'APP_VERSION\s*=\s*"([\d.]+)"', "app.APP_VERSION")
    sites += grab("app.py", r'CoreForge v([\d.]+)', "app.footer")
    sites += grab("report.py", r'VERSION\s*=\s*"([\d.]+)"', "report.VERSION")
    sites += grab("solver/coreforge.f90", r'COREFORGE v?([\d.]+)', "engine")
    versions = sorted({v for _, v in sites})
    ok = len(versions) == 1 and len(sites) >= 5
    tag = versions[0] if len(versions) == 1 else "/".join(versions)
    print(f"{'version lockstep: engine == report == app':<52}"
          f"v{tag} ({len(sites)} sites)  {'PASS' if ok else 'FAIL'}")
    return ok


def main():
    fine = "--fine" in sys.argv
    no_engine = "--no-engine" in sys.argv
    if not no_engine and runner.engine_path() is None:
        print("engine not built and auto-build failed — build manually:")
        print("  ifx -O3 -qopenmp solver/coreforge.f90 -o solver/coreforge")
        print("  gfortran -O3 -fopenmp solver/coreforge.f90 -o solver/coreforge")
        print("or run the engine-free physics subset:  python3 verify.py --no-engine")
        return 2

    if no_engine:
        # engine-free subset: pure-Python physics checks (no Fortran needed).
        # Useful for independent review in sandboxes without a compiler.
        ok = True
        ok &= version_consistency_check()
        ok &= designer_trends()
        ok &= designer_vs_iaea()
        ok &= inverse_designer_check()
        ok &= burnup_trends()
        ok &= integrity_guard_check()
        ok &= thermal_checks()
        ok &= xenon_checks()
        ok &= kinetics_trends()
        print("-" * 78)
        print("all engine-free checks passed" if ok
              else "SOME CHECKS FAILED")
        return 0 if ok else 1

    cases = []
    cases.append(("Homogeneous k-inf", presets.preset_homogeneous(), 2.0))
    cases.append(("Bare square 1-group", bare_square_1g(), 2.0))
    cases.append(("Bare cube 1-group 3D", bare_cube_1g(), 6.0))
    cases.append(("IAEA-2D  (h=2 cm)", presets.preset_iaea2d(), 20.0))
    cases.append(("IAEA-3D  (h=5,dz=10)",
                  dict(presets.preset_iaea3d(),
                       axial=dict(presets.iaea3d_axial(80.0), divz=2)),
                  60.0))
    cases.append(("C5G7 demo (h=1.26)", presets.preset_c5g7(), 300.0))
    cases.append(("SMR-class demo core", presets.preset_smr(), None))

    if fine:
        cfg = presets.preset_iaea2d(); cfg["div"] = 10
        cases.append(("IAEA-2D  (h=1 cm)", cfg, 10.0))
        cfg = presets.preset_iaea3d(); cfg["div"] = 4
        cfg["axial"] = dict(presets.iaea3d_axial(80.0), divz=4)
        cases.append(("IAEA-3D  (h=2.5,dz=5)", cfg, 30.0))
        cfg = presets.preset_c5g7(); cfg["div"] = 2
        cases.append(("C5G7 demo (h=0.63)", cfg, 300.0))

    print(f"{'case':<22}{'k_eff':>12}{'reference':>12}{'diff pcm':>10}"
          f"{'outers':>8}{'time s':>8}  status")
    print("-" * 78)
    ok = True
    for name, cfg, tol_pcm in cases:
        r = runner.run_case(cfg)
        if cfg.get("ref_keff") is None:
            ok &= r.get("converged", False)
            print(f"{name:<22}{r['keff']:>12.7f}{'—':>12}{'—':>10}"
                  f"{r['outers']:>8}{r['time_s']:>8.2f}  "
                  f"{'INFO' if r.get('converged') else 'FAIL'}")
            continue
        d = runner.rho_pcm(r["keff"]) - runner.rho_pcm(cfg["ref_keff"])
        passed = abs(d) <= tol_pcm and r.get("converged", False)
        ok &= passed
        print(f"{name:<22}{r['keff']:>12.7f}{cfg['ref_keff']:>12.7f}"
              f"{d:>+10.1f}{r['outers']:>8}{r['time_s']:>8.2f}  "
              f"{'PASS' if passed else 'FAIL'}")
    ok &= version_consistency_check()
    ok &= equiv_2d_3d()
    ok &= designer_trends()
    ok &= designer_vs_iaea()
    ok &= inverse_designer_check()
    ok &= burnup_trends()
    ok &= flux_calibration_check()
    ok &= integrity_guard_check()
    ok &= eqfuel_reliability_check()
    ok &= thermal_checks()
    ok &= xenon_checks()
    ok &= multicycle_check()
    ok &= rod_position_check()
    ok &= kinetics_trends()
    print("-" * 78)
    print("all checks passed" if ok else "SOME CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
