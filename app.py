"""CoreForge — interactive 2-D/3-D multigroup reactor core simulator.

Streamlit front-end over the Fortran diffusion engine (solver/coreforge.f90).
Workflow, in the spirit of production core simulators:
  materials -> core loading pattern -> solve (2-D or 3-D) -> analyse
  (+ fuel designer, burnup, physics tools, project save/load).

Run:  streamlit run app.py
"""
import copy
import json
import os

import numpy as np
import pandas as pd
import streamlit as st

import burnup
import kinetics
import plots
import presets
import report
import runner
import thermal
import ui
import xslib

APP_VERSION = "8.4"

st.set_page_config(page_title="CoreForge · reactor core simulator",
                   page_icon="⚛️", layout="wide",
                   initial_sidebar_state="expanded")
ui.inject()

# ======================================================================
# session state
# ======================================================================
def load_preset(p):
    ss = st.session_state
    ss.ng = p["ng"]
    ss.pitch = float(p["pitch"])
    ss.div = int(p["div"])
    ss.blocks = [row[:] for row in p["blocks"]]
    ss.materials = copy.deepcopy(p["materials"])
    ss.bc = list(p["bc"])
    ss.gamma = float(p["gamma"])
    ss.axial = copy.deepcopy(p.get("axial"))
    ss.rod_meta = copy.deepcopy(p.get("rod_meta"))
    ss.divz = int(p.get("axial", {}).get("divz", 1)) if p.get("axial") else 1
    ss.ref_keff = p.get("ref_keff")
    ss.ref_source = p.get("ref_source")
    ss.ref_fp = (presets.fingerprint_of(p) if p.get("ref_keff") else None)
    ss.preset_title = p["title"]
    ss.rev = ss.get("rev", 0) + 1
    ss.pop("last", None)
    ss.pop("burn", None)
    ss.pop("burn_params", None)
    ss.pop("eqfit", None)
    ss.pop("eqfit_log", None)
    ss.pop("th", None)          # stale results from ANOTHER core must
    ss.pop("transient", None)   # never survive a preset change
    ss.pop("xet", None)
    ss.pop("cycles", None)


def current_fp():
    """Physics fingerprint of the CURRENT session configuration."""
    ss = st.session_state
    return presets.physics_fingerprint(ss.ng, ss.pitch, ss.blocks,
                                       ss.materials, ss.bc, ss.gamma,
                                       ss.get("axial"))


def ref_snapshot():
    """Reference info frozen at solve time, with QA validity."""
    ss = st.session_state
    ok = bool(ss.ref_keff) and ss.get("ref_fp") == current_fp()
    return dict(keff=ss.ref_keff, source=ss.ref_source, ok=ok)


def project_json():
    """Serialise the full session configuration as a JSON project file."""
    ss = st.session_state
    return json.dumps(dict(
        coreforge_project="5.1",
        title=ss.preset_title, ng=ss.ng, pitch=ss.pitch, div=ss.div,
        divz=int(ss.get("divz", 1)), blocks=ss.blocks,
        materials=ss.materials, bc=ss.bc, gamma=ss.gamma,
        axial=ss.get("axial"), rod_meta=ss.get("rod_meta"),
        ref_keff=ss.get("ref_keff"), ref_source=ss.get("ref_source"),
        ref_fp=ss.get("ref_fp"),
        omega=ss.omega, ninner=int(ss.ninner),
        tolk=ss.tolk, tols=ss.tols), indent=1)


def load_project(d):
    """Load a project dict saved by project_json()."""
    ss = st.session_state
    if "coreforge_project" not in d:
        raise ValueError("not a CoreForge project file")
    p = dict(key="project", title=d.get("title", "Loaded project"),
             description="", ng=int(d["ng"]), pitch=float(d["pitch"]),
             div=int(d["div"]), blocks=d["blocks"],
             materials=d["materials"], bc=list(d["bc"]),
             gamma=float(d["gamma"]), axial=d.get("axial"),
             rod_meta=d.get("rod_meta"),
             ref_keff=d.get("ref_keff"), ref_source=d.get("ref_source"))
    load_preset(p)
    # trust the SAVED fingerprint, not one recomputed from the file: a
    # project saved after modifying a benchmark must stay 'invalidated'
    ss.ref_fp = d.get("ref_fp") if d.get("ref_keff") else None
    if d.get("axial"):
        ss.divz = int(d.get("divz", 1))
    for k in ("omega", "tolk", "tols"):
        if k in d:
            setattr(ss, k, float(d[k]))
    if "ninner" in d:
        ss.ninner = int(d["ninner"])


def blank_project():
    """A clean 9x9 two-material 2-D core to start a design from scratch."""
    return dict(
        key="blank", title="New core (untitled)", description="",
        ng=2, pitch=20.0, div=4,
        blocks=[[2]*9] + [[2] + [1]*7 + [2] for _ in range(7)] + [[2]*9],
        materials=[
            dict(name="Fuel", D=[1.4, 0.4], Sa=[0.01, 0.10],
                 nuSf=[0.007, 0.13], chi=[1.0, 0.0],
                 scat=[[0.0, 0.018], [0.0, 0.0]]),
            dict(name="Water reflector", D=[1.5, 0.3], Sa=[0.0004, 0.011],
                 nuSf=[0.0, 0.0], chi=[0.0, 0.0],
                 scat=[[0.0, 0.045], [0.0, 0.0]]),
        ],
        bc=["vacuum"] * 4, gamma=0.4692,
        ref_keff=None, ref_source=None)


if "rev" not in st.session_state:
    st.session_state.omega = 1.6
    st.session_state.ninner = 4
    st.session_state.tolk = 1e-7
    st.session_state.tols = 1e-5
    load_preset(presets.preset_smr())     # SMR-first (Teknofest focus)

ss = st.session_state
REV = ss.rev
IS3D = bool(ss.get("axial"))


def conform_materials(mats, ng):
    """Pad/truncate every material's group data to NG entries."""
    for m in mats:
        for key, fill in (("D", 1.0), ("Sa", 0.0), ("nuSf", 0.0), ("chi", 0.0)):
            v = list(m[key])[:ng]
            v += [fill] * (ng - len(v))
            m[key] = v
        sc = [row[:ng] + [0.0] * (ng - len(row)) for row in m["scat"][:ng]]
        while len(sc) < ng:
            sc.append([0.0] * ng)
        m["scat"] = sc
    return mats


def current_cfg():
    cfg = dict(ng=ss.ng, pitch=ss.pitch, div=ss.div,
               blocks=[row[:] for row in ss.blocks],
               materials=copy.deepcopy(ss.materials),
               bc=list(ss.bc), gamma=ss.gamma,
               omega=ss.omega, ninner=int(ss.ninner),
               tolk=ss.tolk, tols=ss.tols)
    if ss.get("axial"):
        ax = copy.deepcopy(ss.axial)
        ax["divz"] = int(ss.divz)
        cfg["axial"] = ax
    return cfg


# ======================================================================
# sidebar — engine & solver settings
# ======================================================================
with st.sidebar:
    st.markdown("### ⚛️ CoreForge")
    st.caption("Fortran 2-D/3-D diffusion engine · web console")

    eng = runner.engine_path(autobuild=False)
    if not eng:
        # first run on a fresh host (e.g. Streamlit Cloud): build from source
        try:
            with st.spinner("Fortran engine building from source (first run)…"):
                eng = runner.engine_path()
        except Exception:
            eng = None
    if not eng:
        st.error("Engine not built and auto-build failed "
                 "(see `solver/build.log`). In the project folder run:\n\n"
                 "`ifx -O3 -qopenmp solver/coreforge.f90 -o solver/coreforge`\n\n"
                 "or\n\n"
                 "`gfortran -O3 -fopenmp solver/coreforge.f90 -o solver/coreforge`",
                 icon="🛑")

    ncpu = os.cpu_count() or 1
    if st.checkbox("auto threads (recommended)", True,
                   help="Matches thread count to mesh size — small meshes "
                        "run FASTER on fewer threads."):
        threads = None
    else:
        threads = st.slider("OpenMP threads", 1, ncpu, ncpu)

    dim = st.radio("Solver model", ["2-D (x-y)", "3-D (x-y-z)"],
                   index=1 if IS3D else 0, key=f"dim_{REV}", horizontal=True,
                   help="Switching to 3-D extrudes the current radial map "
                        "into a single axial zone — then add reflector/rod "
                        "zones in the Core builder. Switching to 2-D keeps "
                        "the first zone's radial map.")
    want3d = dim.startswith("3")
    if want3d != IS3D:
        if want3d:
            ss.axial = dict(dz=20.0, divz=1,
                            zones=[dict(label="core",
                                        layers=10,
                                        blocks=[row[:] for row in ss.blocks])])
            ss.divz = 1
            while len(ss.bc) < 6:
                ss.bc.append("vacuum")
        else:
            # keep the most-fissile zone's radial map (zone 0 may be a
            # bottom reflector — taking it would leave a fuel-free core)
            fids = set(runner.fissile_ids(ss.materials))
            zbest = max(ss.axial["zones"],
                        key=lambda z: sum(1 for row in z["blocks"]
                                          for v in row if v in fids))
            ss.blocks = [row[:] for row in zbest["blocks"]]
            ss.axial = None
            ss.rod_meta = None
            ss.bc = ss.bc[:4]
        ss.rev += 1
        st.rerun()

    ss.div = st.slider("Cells per block (radial refinement)", 1, 10,
                       int(ss.div), key=f"div_{REV}")
    nby = len(ss.blocks)
    nbx = len(ss.blocks[0])
    if IS3D:
        ss.divz = st.slider("Axial refinement (divz)", 1, 4,
                            int(ss.divz), key=f"divz_{REV}")
        nzl = sum(int(z["layers"]) for z in ss.axial["zones"]) * int(ss.divz)
        ncells = nbx * ss.div * nby * ss.div * nzl
        st.caption(f"mesh **{nbx*ss.div} × {nby*ss.div} × {nzl}** cells · "
                   f"Δxy = {ss.pitch/ss.div:.2f} cm · "
                   f"Δz = {ss.axial['dz']/ss.divz:.2f} cm")
    else:
        ncells = nbx * ss.div * nby * ss.div
        st.caption(f"mesh **{nbx*ss.div} × {nby*ss.div}** cells · "
                   f"Δ = {ss.pitch/ss.div:.3f} cm")
    if ncells * ss.ng > 2_000_000:
        st.warning("large problem — the solve may take a while", icon="⏳")

    st.subheader("Boundary conditions")
    bc_names = ["vacuum", "reflective"]
    labels = ["West", "East", "South", "North"] + \
             (["Bottom", "Top"] if IS3D else [])
    while len(ss.bc) < len(labels):
        ss.bc.append("vacuum")
    newbc = []
    cols = st.columns(2)
    for kk, lab in enumerate(labels):
        with cols[kk % 2]:
            newbc.append(st.selectbox(lab, bc_names,
                                      index=bc_names.index(ss.bc[kk]),
                                      key=f"bc{kk}_{REV}"))
    ss.bc = newbc

    with st.expander("Advanced solver settings"):
        ss.gamma = st.number_input("vacuum γ (J/φ at face)", 0.01, 1e9,
                                   float(ss.gamma), format="%.4f",
                                   key=f"gam_{REV}",
                                   help="0.4692 transport-corrected, 0.5 Marshak")
        ss.omega = st.number_input("SOR ω", 0.5, 1.95, float(ss.omega), 0.05)
        ss.ninner = st.number_input("inner sweeps / group", 1, 20, int(ss.ninner))
        ss.tolk = st.number_input("k tolerance", 1e-9, 1e-4, float(ss.tolk),
                                  format="%.1e")
        ss.tols = st.number_input("source tolerance", 1e-8, 1e-3, float(ss.tols),
                                  format="%.1e")

    st.divider()
    st.subheader("Project")
    if st.button("🆕 New blank core", width="stretch"):
        load_preset(blank_project())
        st.rerun()
    st.download_button("💾 Save project (.json)", project_json(),
                       "coreforge_project.json", mime="application/json",
                       width="stretch")
    up = st.file_uploader("📂 Load project", type=["json"],
                          key="proj_upload")
    if up is not None:
        sig = (up.name, up.size)
        if ss.get("_proj_sig") != sig:
            try:
                load_project(json.loads(up.getvalue().decode("utf-8")))
                ss._proj_sig = sig
                st.rerun()
            except Exception as e:
                ss._proj_sig = sig
                st.error(f"could not load project: {e}")

    st.divider()
    st.caption(f"loaded: **{ss.preset_title}**")
    if ss.ref_keff:
        st.caption(f"reference k_eff = {ss.ref_keff}")

# ======================================================================
ui.hero(eng is not None, APP_VERSION,
        f"SMR/MMR neutronics · loaded: {ss.preset_title}"
        + ("  ·  3-D (x-y-z)" if IS3D else "  ·  2-D (x-y)"))

tabs = st.tabs(["📚 Benchmarks", "🧬 Fuel designer", "🧪 Materials",
                "🧱 Core builder", "⚡ Solve & results",
                "🌡 Thermal-hydraulics", "🔥 Burnup",
                "⏱ Transient", "🔧 Physics tools"])

# ----------------------------------------------------------------------
# 📚 Benchmarks
# ----------------------------------------------------------------------
with tabs[0]:
    ui.section("Bundled cores & benchmarks",
               "Load a validated configuration with one click, then edit "
               "anything. Cards tagged 🧊 3-D run the full x-y-z engine "
               "automatically.")
    with st.expander("▶  New here? 60-second quick start", expanded=False):
        qa, qb, qc = st.columns(3)
        qa.markdown("**1 · Load** a core below (start with **SMR-class "
                    "PWR**) → jump to **⚡ Solve** and press *Solve k_eff*.")
        qb.markdown("**2 · Analyse**: **🌡 Thermal-hydraulics** for "
                    "temperatures, **🔥 Burnup** for the fuel cycle, "
                    "**⏱ Transient** for rod ejection / xenon.")
        qc.markdown("**3 · Export**: the **📄 HTML report** button on the "
                    "Solve tab bundles everything into one shareable file. "
                    "**💾 Save project** keeps your design.")
    _plist = [pf() for pf in presets.ALL_PRESETS]
    _rows = [_plist[i:i + 3] for i in range(0, len(_plist), 3)]
    for _row in _rows:
        cols = st.columns(3)
        for col, p in zip(cols, _row):
            with col, st.container(border=True):
                dim_tag = "🧊 3-D" if p.get("axial") else "▦ 2-D"
                st.subheader(f"{p['title']}")
                st.caption(dim_tag)
                blk = np.flipud(np.asarray(p["blocks"]))
                xc = (np.arange(blk.shape[1]) + 0.5) * p["pitch"]
                yc = (np.arange(blk.shape[0]) + 0.5) * p["pitch"]
                st.plotly_chart(
                    plots.material_map_fig(blk, p["materials"], xc, yc,
                                           mini=True),
                    width="stretch", key=f"thumb_{p['key']}",
                    config={"displayModeBar": False})
                st.write(p["description"])
                if p["ref_keff"]:
                    st.caption(f"reference k_eff = **{p['ref_keff']}** "
                               f"({p['ref_source']})")
                if st.button("Load", key=f"load_{p['key']}", type="primary"):
                    load_preset(p)
                    st.rerun()

    st.subheader("Validation (this engine)")
    st.dataframe(pd.DataFrame([
        dict(case="Homogeneous k∞ (all reflective)", model="2-D", mesh="8×8",
             keff="1.0857142", reference="1.0857143 (analytic)", diff="−0.0 pcm"),
        dict(case="Bare square, 1 group (analytic B²)", model="2-D",
             mesh="120×120", keff="1.2610188", reference="1.2610203",
             diff="−0.1 pcm"),
        dict(case="Bare cube, 1 group (analytic B²)", model="3-D",
             mesh="24×24×24", keff="≈1.24897", reference="1.2489663",
             diff="~±5 pcm"),
        dict(case="IAEA-2D PWR (11-A2)", model="2-D", mesh="170×170 (h=1 cm)",
             keff="1.0295785", reference="1.02959", diff="−1.1 pcm"),
        dict(case="IAEA-3D PWR (problem 11)", model="3-D",
             mesh="85×85×95 (h=2, dz=4)", keff="1.0289590",
             reference="1.02903", diff="−6.7 pcm"),
        dict(case="C5G7 MOX, pin-cell diffusion demo", model="2-D",
             mesh="102×102 (h=0.63 cm)", keff="1.1863920",
             reference="1.18655 (MCNP transport)", diff="−11 pcm *"),
    ]), hide_index=True, width="stretch")
    st.caption("Monotone mesh convergence toward every published eigenvalue; "
               "vacuum faces use the transport-corrected Robin condition "
               "(γ = 0.4692). IAEA-3D converges −27 → −12.6 → −6.7 pcm as "
               "the mesh refines. *C5G7 is a transport benchmark solved with "
               "volume-homogenised pin cells: the eigenvalue agrees thanks "
               "to favourable error cancellation while local pin powers "
               "deviate — exactly the transport effect the demo illustrates.")

# ----------------------------------------------------------------------
# 🧬 Fuel designer
# ----------------------------------------------------------------------
with tabs[1]:
    st.header("Fuel designer — physical inputs → two-group constants")
    st.write("A built-in pin-cell model turns **enrichment, soluble boron "
             "and moderator density** into homogenised two-group cross "
             "sections (real 2200 m/s microscopic data; fast-group and "
             "spectrum constants calibrated once to a nominal PWR cell — "
             "see README). The way a lattice code feeds a core simulator.")
    dc1, dc2, dc3, dc4 = st.columns(4)
    with dc1:
        d_wat = st.checkbox("pure water cell (no fuel)")
        d_e = st.number_input("U-235 enrichment [w/o]", 0.7, 5.0, 3.1, 0.05,
                              disabled=d_wat)
    with dc2:
        d_ppm = st.number_input("soluble boron [ppm]", 0.0, 5000.0, 0.0, 50.0)
    with dc3:
        d_rho = st.number_input("moderator density [g/cm³]", 0.50, 1.00,
                                0.71, 0.01)
    with dc4:
        d_rf = st.number_input("fuel radius [cm]", 0.30, 0.60, 0.4095, 0.005,
                               disabled=d_wat)

    dm = xslib.pincell_xs(0.0 if d_wat else d_e, d_ppm, d_rho,
                          0.0 if d_wat else d_rf)
    k1c, k2c, k3c = st.columns(3)
    k1c.metric("k∞ of this cell", f"{dm['designer']['kinf']:.4f}"
               if not d_wat else "—")
    k2c.metric("Σa2 boron sensitivity",
               f"{dm['designer']['dsa2_dppm']*1e6:.2f} ×10⁻⁶ cm⁻¹/ppm")
    k3c.metric("Σs 1→2", f"{dm['scat'][0][1]:.4f} cm⁻¹")
    st.dataframe(pd.DataFrame(
        dict(group=[1, 2], D=dm["D"], Sa=dm["Sa"], nuSf=dm["nuSf"],
             chi=dm["chi"])), hide_index=True, width="stretch")
    d_name = st.text_input("material name", dm["name"])
    if st.button("➕  Add to materials", type="primary"):
        if ss.ng != 2:
            st.error("The fuel designer generates 2-group materials — "
                     "set NG = 2 first (Materials tab).")
        elif len(ss.materials) >= 8:
            st.error("Material limit (8) reached.")
        else:
            dm["name"] = d_name
            ss.materials.append(dm)
            ss.rev += 1
            st.toast(f"added: {d_name}")
            st.rerun()
    st.caption("Educational lattice physics: trends (enrichment, boron, "
               "density) are physical; absolute k is representative, not "
               "lattice-code accuracy.")

    st.divider()
    st.subheader("🔁 Equivalent fuel — match a composition to existing XS")
    st.write("Have a material defined only by macroscopic cross sections "
             "(e.g. a benchmark fuel)? The inverse designer finds the "
             "fresh fuel (enrichment + boron) whose generated set matches "
             "its thermal fission and absorption **exactly**, and reports "
             "the residuals on everything else. Attaching the equivalent "
             "makes the material **depletable** in the 🔥 Burnup tab.")
    cand = [(i, m) for i, m in enumerate(ss.materials)
            if m["nuSf"][min(1, ss.ng - 1)] > 0 and "designer" not in m]
    if ss.ng != 2:
        st.info("The inverse designer works on 2-group materials (NG = 2).")
    elif not cand:
        st.info("No XS-only fuels in the current core — every fissile "
                "material already carries a composition.")
    else:
        names_c = [f"{i+1} · {m['name']}" for i, m in cand]
        sel = st.selectbox("material to match", names_c)
        if st.button("Fit equivalent fuel"):
            idx = cand[names_c.index(sel)][0]
            try:
                fit = xslib.match_fuel_to_xs(ss.materials[idx])
                ss.eqfit = (idx, fit)
            except Exception as e:
                st.error(str(e))
        if "eqfit" in ss:
            idx, fit = ss.eqfit
            f1, f2, f3 = st.columns(3)
            f1.metric("equivalent enrichment", f"{fit['enrich']:.3f} w/o")
            f2.metric("equivalent boron", f"{fit['ppm']:.0f} ppm")
            f3.metric("k∞ target → equivalent",
                      f"{fit['residuals']['kinf_target']:.4f} → "
                      f"{fit['residuals']['kinf_eq']:.4f}")
            rdf = pd.DataFrame([{k: v for k, v in fit["residuals"].items()
                                 if isinstance(v, float) and
                                 not k.startswith("kinf")}])
            st.dataframe(rdf.round(2), hide_index=True, width="stretch")
            st.caption("residuals in % (matched: nuSf2, Sa2 ≈ 0). The k∞ "
                       "difference reflects the generator's own fast-group "
                       "calibration — the equivalent is an approximate, "
                       "clearly-labelled stand-in that enables depletion.")
            if fit.get("note"):
                st.warning(fit["note"])
            reliable, wkey, wval = xslib.reliability_flag(fit["residuals"])
            if not reliable:
                st.error(
                    f"⚠️ **Low reliability for reactivity-sensitive use** "
                    f"— {wkey} residual is {wval:+.1f}% (threshold ±15%). "
                    f"This usually means the target material has an "
                    f"absorption signature a boron-water model cannot "
                    f"reproduce (e.g. an **inserted control rod / "
                    f"burnable poison**, not a plain fuel cell). k_eff "
                    f"and reactivity worths computed with this equivalent "
                    f"will be measurably off for cores that depend on "
                    f"this material's worth — verified case: substituting "
                    f"an SMR's CRA (control-rod) material this way shifted "
                    f"core k_eff by **+2,446 pcm** on its own. Safe to "
                    f"attach for depletion bookkeeping; do NOT trust "
                    f"rods-in/rods-out k_eff comparisons that depend on "
                    f"this material.")
            if st.button("✅ Attach (replace XS with the equivalent)",
                         type="primary"):
                nm = fit["material"]
                nm["name"] = ss.materials[idx]["name"] + " (eq-fuel)"
                ss.setdefault("eqfit_log", []).append(dict(
                    name=nm["name"],
                    target=copy.deepcopy(ss.materials[idx]),
                    material=copy.deepcopy(nm),
                    enrich=fit["enrich"], ppm=fit["ppm"],
                    residuals=dict(fit["residuals"]),
                    note=fit.get("note", ""),
                    reliable=reliable, worst_key=wkey, worst_pct=wval))
                ss.materials[idx] = nm
                ss.pop("eqfit")
                ss.rev += 1
                st.toast("equivalent fuel attached — now depletable")
                st.rerun()

# ----------------------------------------------------------------------
# 🧪 Materials
# ----------------------------------------------------------------------
with tabs[2]:
    st.header("Materials & cross sections")

    c1, c2 = st.columns(2)
    with c1:
        new_ng = st.number_input("Energy groups (NG)", 1, 8, int(ss.ng),
                                 key=f"ng_{REV}")
    with c2:
        new_nm = st.number_input("Number of materials", 1, 8,
                                 len(ss.materials), key=f"nm_{REV}")

    if new_ng != ss.ng:
        ss.ng = int(new_ng)
        conform_materials(ss.materials, ss.ng)
        ss.rev += 1
        st.rerun()
    while len(ss.materials) < new_nm:
        ss.materials.append(dict(
            name=f"material {len(ss.materials)+1}",
            D=[1.0] * ss.ng, Sa=[0.01] * ss.ng, nuSf=[0.0] * ss.ng,
            chi=[1.0] + [0.0] * (ss.ng - 1),
            scat=[[0.0] * ss.ng for _ in range(ss.ng)]))
    if len(ss.materials) > new_nm:
        ss.materials = ss.materials[:int(new_nm)]
        ss.blocks = [[min(v, int(new_nm)) for v in row] for row in ss.blocks]

    conform_materials(ss.materials, ss.ng)

    rows = []
    for m in ss.materials:
        r = {"name": m["name"]}
        for g in range(ss.ng):
            r[f"D{g+1}"] = m["D"][g]
        for g in range(ss.ng):
            r[f"Sa{g+1}"] = m["Sa"][g]
        for g in range(ss.ng):
            r[f"nuSf{g+1}"] = m["nuSf"][g]
        for g in range(ss.ng):
            r[f"chi{g+1}"] = m["chi"][g]
        rows.append(r)
    df = pd.DataFrame(rows, index=[f"mat {i+1}" for i in range(len(rows))])
    st.markdown(" ".join(
        f'<span style="display:inline-block;width:0.85em;height:0.85em;'
        f'background:{plots.CAT[i % len(plots.CAT)]};border-radius:2px;'
        f'margin-right:0.25em"></span><small>{i+1}·{m["name"]}</small>&nbsp;&nbsp;'
        for i, m in enumerate(ss.materials)), unsafe_allow_html=True)
    st.caption("D = diffusion coefficient [cm] · Sa = Σ_absorption [cm⁻¹] · "
               "nuSf = ν·Σ_fission [cm⁻¹] · chi = fission spectrum")
    edited = st.data_editor(
        df, width="stretch", key=f"mats_{REV}",
        column_config={"name": st.column_config.TextColumn("name")})
    for i, m in enumerate(ss.materials):
        m["name"] = str(edited.iloc[i]["name"])
        m["D"] = [float(edited.iloc[i][f"D{g+1}"]) for g in range(ss.ng)]
        m["Sa"] = [float(edited.iloc[i][f"Sa{g+1}"]) for g in range(ss.ng)]
        m["nuSf"] = [float(edited.iloc[i][f"nuSf{g+1}"]) for g in range(ss.ng)]
        m["chi"] = [float(edited.iloc[i][f"chi{g+1}"]) for g in range(ss.ng)]

    st.subheader("Scattering matrices Σs(g → g′)")
    st.caption("Row = source group, column = destination group. "
               "Up-scatter is allowed; the g→g diagonal is ignored "
               "(within-group scatter cancels in removal).")
    for i, m in enumerate(ss.materials):
        with st.expander(f"mat {i+1} · {m['name']}"):
            sdf = pd.DataFrame(
                m["scat"],
                index=[f"from g{g+1}" for g in range(ss.ng)],
                columns=[f"to g{g+1}" for g in range(ss.ng)])
            sed = st.data_editor(sdf, key=f"scat_{i}_{REV}",
                                 width="stretch")
            m["scat"] = [[float(sed.iloc[a, b]) for b in range(ss.ng)]
                         for a in range(ss.ng)]

# ----------------------------------------------------------------------
# 🧱 Core builder
# ----------------------------------------------------------------------
def _resize_map(mp, nrows, ncols, fill=1):
    b = [row[:ncols] + [fill] * max(0, ncols - len(row)) for row in mp]
    return b[:nrows] + [[fill] * ncols for _ in range(max(0, nrows - len(b)))]


with tabs[3]:
    st.header("Core loading pattern")
    c1, c2, c3 = st.columns(3)
    with c1:
        ss.pitch = st.number_input("Block pitch [cm]", 0.5, 100.0,
                                   float(ss.pitch), key=f"pitch_{REV}")
    with c2:
        nrows = st.number_input("Rows", 3, 64, len(ss.blocks),
                                key=f"nr_{REV}")
    with c3:
        ncols = st.number_input("Columns", 3, 64, len(ss.blocks[0]),
                                key=f"nc_{REV}")
    if (int(nrows), int(ncols)) != (len(ss.blocks), len(ss.blocks[0])):
        ss.blocks = _resize_map(ss.blocks, int(nrows), int(ncols))
        if IS3D:
            for z in ss.axial["zones"]:
                z["blocks"] = _resize_map(z["blocks"], int(nrows), int(ncols))
        ss.rev += 1
        st.rerun()

    if IS3D:
        st.write("**3-D core** — stack axial zones (bottom → top), each "
                 "with its own radial block map. Add reflector zones, "
                 "rodded zones, anything; the engine represents every "
                 "layer exactly.")
        ax = ss.axial
        zc1, zc2 = st.columns([1, 3])
        with zc1:
            newdz = st.number_input("layer thickness dz [cm]", 1.0, 50.0,
                                    float(ax["dz"]), key=f"dz_{REV}")
            if newdz != ax["dz"]:
                ax["dz"] = float(newdz)
        nzl_tot = sum(int(z["layers"]) for z in ax["zones"])
        with zc2:
            st.metric("total height",
                      f"{nzl_tot * ax['dz']:.0f} cm  ·  {nzl_tot} layers × "
                      f"{ax['dz']:.0f} cm  ·  {len(ax['zones'])} zones")

        # zone selector with physical z-ranges (list is bottom -> top)
        zlabels, zcur = [], 0.0
        for zi_, z in enumerate(ax["zones"]):
            h = z["layers"] * ax["dz"]
            zlabels.append(f"{zi_+1} · {z['label']}  "
                           f"[z = {zcur:.0f}–{zcur+h:.0f} cm]")
            zcur += h
        zsel = st.selectbox("Active axial zone (1 = bottom)", zlabels,
                            index=min(1, len(zlabels) - 1),
                            key=f"zone_{REV}")
        zi = zlabels.index(zsel)
        zone = ax["zones"][zi]

        e1, e2, e3 = st.columns([2, 1, 3])
        with e1:
            zone["label"] = st.text_input("zone label", zone["label"],
                                          key=f"zl_{zi}_{REV}")
        with e2:
            zone["layers"] = int(st.number_input(
                "layers", 1, 200, int(zone["layers"]), key=f"zn_{zi}_{REV}"))
        with e3:
            st.write("")
            b1, b2, b3, b4, b5 = st.columns(5)
            if b1.button("➕ above", key=f"za_{zi}_{REV}",
                         help="insert a copy of this zone above (toward top)"):
                ax["zones"].insert(zi + 1, copy.deepcopy(zone))
                ss.rev += 1
                st.rerun()
            if b2.button("➕ below", key=f"zb_{zi}_{REV}",
                         help="insert a copy of this zone below (toward bottom)"):
                ax["zones"].insert(zi, copy.deepcopy(zone))
                ss.rev += 1
                st.rerun()
            if b3.button("⬆", key=f"zu_{zi}_{REV}", help="move toward top",
                         disabled=zi == len(ax["zones"]) - 1):
                ax["zones"][zi], ax["zones"][zi + 1] = \
                    ax["zones"][zi + 1], ax["zones"][zi]
                ss.rev += 1
                st.rerun()
            if b4.button("⬇", key=f"zd_{zi}_{REV}", help="move toward bottom",
                         disabled=zi == 0):
                ax["zones"][zi], ax["zones"][zi - 1] = \
                    ax["zones"][zi - 1], ax["zones"][zi]
                ss.rev += 1
                st.rerun()
            if b5.button("🗑", key=f"zx_{zi}_{REV}", help="delete this zone",
                         disabled=len(ax["zones"]) == 1):
                ax["zones"].pop(zi)
                ss.rev += 1
                st.rerun()

        f1, f2, f3 = st.columns([2, 1, 2])
        with f1:
            src = st.selectbox("copy radial map from",
                               [l for j, l in enumerate(zlabels) if j != zi],
                               key=f"cpsrc_{zi}_{REV}")
        with f2:
            st.write("")
            if st.button("↷ copy", key=f"cp_{zi}_{REV}"):
                zone["blocks"] = copy.deepcopy(
                    ax["zones"][zlabels.index(src)]["blocks"])
                ss.rev += 1
                st.rerun()
        with f3:
            fm1, fm2 = st.columns([2, 1])
            fillid = fm1.number_input("fill material id", 1,
                                      len(ss.materials), 1,
                                      key=f"fid_{zi}_{REV}")
            fm2.write("")
            if fm2.button("▩ fill", key=f"fill_{zi}_{REV}",
                          help="fill the whole zone with this material"):
                zone["blocks"] = [[int(fillid)] * len(ss.blocks[0])
                                  for _ in range(len(ss.blocks))]
                ss.rev += 1
                st.rerun()
        edit_blocks = zone["blocks"]
    else:
        st.write("The core is a lattice of **blocks** (assemblies / coarse "
                 "zones). Type a material id into each block — the "
                 "refinement slider in the sidebar subdivides every block "
                 "into cells for the solver. Flip the sidebar switch to "
                 "extrude this core into 3-D.")
        edit_blocks = ss.blocks

    left, right = st.columns([1.15, 1.0])
    with left:
        bdf = pd.DataFrame(edit_blocks,
                           index=[f"r{r+1}" for r in range(len(edit_blocks))],
                           columns=[f"c{c+1}"
                                    for c in range(len(edit_blocks[0]))])
        key_sfx = f"z{zi}_" if IS3D else ""
        bed = st.data_editor(bdf, width="stretch",
                             key=f"blocks_{key_sfx}{REV}")
        nmat = len(ss.materials)
        newb = [[int(np.clip(int(bed.iloc[r, c]), 1, nmat))
                 for c in range(bed.shape[1])]
                for r in range(bed.shape[0])]
        if IS3D:
            ss.axial["zones"][zi]["blocks"] = newb
        else:
            ss.blocks = newb
        st.caption("row 1 = top of the core (max y)")
    with right:
        blk = np.flipud(np.asarray(newb))
        xc = (np.arange(blk.shape[1]) + 0.5) * ss.pitch
        yc = (np.arange(blk.shape[0]) + 0.5) * ss.pitch
        st.plotly_chart(
            plots.material_map_fig(blk, ss.materials, xc, yc,
                                   title="Loading pattern"
                                         + (f" — zone {zi+1}" if IS3D else "")),
            width="stretch")

# ----------------------------------------------------------------------
# ⚡ Solve & results
# ----------------------------------------------------------------------
with tabs[4]:
    st.header("Solve")
    nby, nbx = len(ss.blocks), len(ss.blocks[0])
    dim_note = ""
    if IS3D:
        nzl = sum(int(z["layers"]) for z in ss.axial["zones"]) * int(ss.divz)
        dim_note = f" · **3-D**: {nzl} layers × {ss.axial['dz']/ss.divz:.1f} cm"
    st.write(f"**{ss.preset_title}** — {nbx}×{nby} blocks · "
             f"{nbx*ss.div}×{nby*ss.div} radial mesh{dim_note} · "
             f"{ss.ng} groups · {len(ss.materials)} materials · "
             f"BC = {', '.join(ss.bc)}")

    if st.button("▶  Solve k_eff", type="primary", disabled=eng is None):
        try:
            with st.spinner("Fortran engine running…"):
                cfg = current_cfg()
                res = runner.run_case(cfg, threads=threads)
            refinfo = ref_snapshot()
            ss.last = (res, cfg, refinfo)
            ss.setdefault("history", []).append(dict(
                core=ss.preset_title,
                model="3D" if cfg.get("axial") else "2D",
                mesh=f"{len(cfg['blocks'][0])*cfg['div']}×"
                     f"{len(cfg['blocks'])*cfg['div']}"
                     + (f"×{res['nz']}" if res.get("nz", 1) > 1 else ""),
                groups=cfg["ng"], keff=round(res["keff"], 6),
                rho_pcm=round(runner.rho_pcm(res["keff"]), 1),
                Fq=round(res["fxy"], 3), time_s=res["time_s"],
                ref=("ok" if refinfo["ok"] else
                     ("MODIFIED" if refinfo["keff"] else "—"))))
        except Exception as e:
            st.error(f"Solver failed:\n\n```\n{e}\n```")

    if "last" in ss:
        res, cfg, refinfo = ss.last
        k = res["keff"]
        is3d_res = res.get("nz", 1) > 1
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("k_eff", f"{k:.6f}")
        m2.metric("reactivity ρ", f"{runner.rho_pcm(k):+,.0f} pcm")
        if refinfo["keff"] and refinfo["ok"]:
            m3.metric("Δ vs reference",
                      f"{runner.rho_pcm(k)-runner.rho_pcm(refinfo['keff']):+.1f} pcm",
                      help=str(refinfo["source"]))
        elif refinfo["keff"]:
            m3.metric("Δ vs reference", "invalidated",
                      help="Geometry, materials or BCs were modified "
                           "after loading the benchmark — the published "
                           "reference no longer applies to this case. "
                           "(Mesh refinement alone does NOT invalidate.)")
        else:
            m3.metric("state", "supercritical" if k > 1.001 else
                      ("critical" if abs(k - 1) <= 1e-3 else "subcritical"))
        m4.metric("F_q (cell peaking)", f"{res['fxy']:.3f}")
        if refinfo["keff"] and not refinfo["ok"]:
            st.warning("⚠️ QA: benchmark configuration was **modified** — "
                       "the reference comparison is suppressed. Reload "
                       "the preset to restore a valid comparison.")
        m5, m6, m7, m8 = st.columns(4)
        m5.metric("outer iterations", res["outers"])
        m6.metric("F_z (axial peaking)" if is3d_res else "dominance ratio",
                  f"{res['fz']:.3f}" if is3d_res
                  else f"{res.get('dr', 0):.3f}")
        m7.metric("solve time", f"{res['time_s']:.2f} s")
        m8.metric("converged", "yes" if res.get("converged") else "NO ⚠️")

        fuel_ids = runner.fissile_ids(cfg["materials"])
        fuelmask = np.isin(res["mat_fine"], fuel_ids)

        r1a, r1b = st.columns(2)
        with r1a:
            st.plotly_chart(plots.material_map_fig(
                res["mat_fine"], cfg["materials"], res["x"], res["y"],
                title="Core layout (radial)"), width="stretch")
        with r1b:
            st.plotly_chart(plots.power_fig(
                res["power"], fuelmask, res["x"], res["y"],
                title="Relative power"
                      + (" (axial average)" if is3d_res else "")),
                width="stretch")

        gcols = st.columns(min(cfg["ng"], 4))
        for g in range(cfg["ng"]):
            with gcols[g % len(gcols)]:
                st.plotly_chart(plots.field_fig(
                    res["phi"][g], res["x"], res["y"],
                    title=f"Flux — group {g+1}"
                          + (" (axial avg.)" if is3d_res else ""),
                    ramp=plots.GROUP_RAMPS[g % len(plots.GROUP_RAMPS)]),
                    width="stretch")

        if is3d_res:
            st.subheader("Axial view (3-D)")
            a1, a2 = st.columns(2)
            with a1:
                st.plotly_chart(plots.traverse_fig(
                    res["z"], [("axial power", res["axial_power"])],
                    title="Axial power profile", xlabel="z [cm]",
                    ylabel="P/P̄ (layer)"), width="stretch")
            with a2:
                lay = st.slider("Layer (z)", 1, res["nz"],
                                max(1, res["nz"] // 2))
                zc = res["z"][lay - 1]
                st.plotly_chart(plots.field_fig(
                    res["phi3d"][cfg["ng"] - 1][lay - 1],
                    res["x"], res["y"],
                    title=f"Thermal flux at z = {zc:.0f} cm",
                    ramp=plots.GROUP_RAMPS[(cfg["ng"] - 1)
                                           % len(plots.GROUP_RAMPS)]),
                    width="stretch")

        r3a, r3b = st.columns(2)
        with r3a:
            bp = runner.block_powers(res, cfg)
            st.plotly_chart(plots.block_power_fig(bp, cfg["pitch"]),
                            width="stretch")
        with r3b:
            ysel = st.slider("Traverse at y [cm]", float(res["y"][0]),
                             float(res["y"][-1]),
                             float(res["y"][len(res["y"]) // 2]))
            jrow = int(np.argmin(np.abs(res["y"] - ysel)))
            series = [(f"group {g+1}", res["phi"][g][jrow, :])
                      for g in range(cfg["ng"])]
            st.plotly_chart(plots.traverse_fig(
                res["x"], series,
                title=f"Flux traverse at y = {res['y'][jrow]:.1f} cm"),
                width="stretch")

        with st.expander("🔌 Operating point — absolute units "
                         "(flux, power density, energy)"):
            op1, op2 = st.columns(2)
            with op1:
                P_MW = st.number_input("core thermal power [MW]", 0.1,
                                       5000.0, 160.0, key="op_pmw")
            with op2:
                if is3d_res:
                    Hcm = None
                    st.caption("fuel height taken from the 3-D geometry")
                else:
                    Hcm = st.number_input("active fuel height [cm]", 10.0,
                                          500.0, 200.0, key="op_h")
            # fuel volume and scaling: engine normalises fuel-cell-average
            # relative fission power to 1
            if is3d_res:
                dzf = res["z"][1] - res["z"][0] if res["nz"] > 1 else 1.0
                nfuel3 = int((res["power3d"] > 1e-12).sum())
                Vf = nfuel3 * res["dx"] * res["dx"] * dzf          # cm^3
            else:
                nfuel2 = int(fuelmask.sum())
                Vf = nfuel2 * res["dx"] * res["dx"] * Hcm
            q_avg = P_MW * 1e6 / Vf                                # W/cm^3
            q_peak = q_avg * res["fxy"]
            # Absolute flux scale S: derived from total-power conservation
            # (Sum_cells [Sum_g nuSf_g * phi_g_rel] * V_cell) * S = P/kappa.
            # Since the engine normalises phi such that the fuel-cell
            # AVERAGE relative fission source is 1, the cell sum reduces to
            # exactly V_fuel, giving S = q_avg/kappa -- NO nu-bar factor
            # (a second nu would double-count; nu is already inside nuSf).
            # Verified against a 1-group test case: this reproduces the
            # target power to 4 significant digits.
            KAPPA = 3.204e-11
            s_flux = q_avg / KAPPA
            phi_abs = [float(res["phi"][g].max()) * s_flux
                       for g in range(cfg["ng"])]
            o1, o2, o3, o4 = st.columns(4)
            o1.metric("fuel volume", f"{Vf/1e6:.2f} m³")
            o2.metric("avg power density", f"{q_avg:.1f} W/cm³")
            o3.metric("peak power density", f"{q_peak:.1f} W/cm³")
            o4.metric("peak thermal flux",
                      f"{phi_abs[-1]:.2e} n/cm²·s")
            bp_abs = runner.block_powers(res, cfg)
            nfb = int(np.isfinite(bp_abs).sum())
            o5, o6, o7, o8 = st.columns(4)
            o5.metric("fuel blocks", nfb)
            o6.metric("avg assembly power", f"{P_MW/nfb:.2f} MW")
            o7.metric("hottest assembly",
                      f"{np.nanmax(bp_abs)*P_MW/nfb:.2f} MW")
            o8.metric("peak fast flux", f"{phi_abs[0]:.2e} n/cm²·s")
            st.caption("Absolute scaling from the stated core power: "
                       "q‴ = P/V_fuel; flux scale S = q‴/κ (κ = 200 "
                       "MeV/fission) is fixed by total fission-rate "
                       "conservation, applied uniformly to every group "
                       "(no separate ν̄ factor — ν is already inside νΣf). "
                       "Typical PWR checks: q‴ ~ 50–110 W/cm³, thermal "
                       "flux ~ 3–6 ×10¹³ n/cm²·s.")

        if "balance" in res:
            st.subheader("Neutron balance")
            bdf2 = pd.DataFrame(res["balance"]["rows"],
                                columns=res["balance"]["columns"])
            st.dataframe(bdf2, hide_index=True, width="stretch")
            st.caption("production = χ·F/k emitted in the group; net leakage "
                       "closes the balance. All values in relative units.")

        d1, d2, d3, d4 = st.columns(4)
        d1.download_button("⬇ input.txt", res["input_text"], "input.txt")
        with open(os.path.join(res["workdir"], "flux.csv"), "rb") as f:
            d2.download_button("⬇ flux.csv", f.read(), "flux.csv")
        with open(os.path.join(res["workdir"], "power.csv"), "rb") as f:
            d3.download_button("⬇ power.csv", f.read(), "power.csv")
        _ref_note = None
        if refinfo["keff"] and not refinfo["ok"]:
            _ref_note = ("Benchmark configuration was modified after "
                         "loading the preset — the published reference "
                         f"({refinfo['keff']}) no longer applies and the "
                         "comparison is suppressed.")
        d4.download_button(
            "📄 HTML report",
            report.build_html(cfg, ss.preset_title,
                              refinfo["keff"] if refinfo["ok"] else None,
                              res,
                              ss.burn[0] if "burn" in ss else None,
                              power_MW=float(ss.get("op_pmw", 160.0)),
                              fuel_height_cm=float(ss.get("op_h", 200.0)),
                              ref_note=_ref_note,
                              eqfits=ss.get("eqfit_log"),
                              th=ss.get("th")),
            "coreforge_report.html", mime="text/html")
        with st.expander("Engine log"):
            st.code(res["log"])

    if ss.get("history"):
        st.subheader("Run history (this session)")
        hc1, hc2 = st.columns([5, 1])
        with hc1:
            st.dataframe(pd.DataFrame(ss.history), hide_index=True,
                         width="stretch")
        with hc2:
            if st.button("clear", key="clear_hist"):
                ss.history = []
                st.rerun()

# ----------------------------------------------------------------------
# 🌡 Thermal-hydraulics (closed channel, steady state)
# ----------------------------------------------------------------------
with tabs[5]:
    st.header("Thermal-hydraulics — closed-channel steady state")
    st.write("Single-phase PWR channel analysis driven by the neutronics "
             "solution: coolant heat-up along the channel, clad surface "
             "and **fuel centreline temperatures** (Dittus-Boelter + "
             "conduction resistances), friction **pressure drop** and the "
             "saturation margin — for the **average** and the **hot** "
             "assembly (radial peaking from the block-power map, axial "
             "shape from the 3-D profile or a chopped cosine).")
    if "last" not in ss:
        st.info("Solve the core first (⚡ tab) — the channel analysis "
                "uses its power distribution.")
    else:
        res_t, cfg_t, _ri = ss.last
        bp_t = runner.block_powers(res_t, cfg_t)
        fdh = float(np.nanmax(bp_t))
        t1, t2, t3, t4 = st.columns(4)
        with t1:
            th_P = st.number_input("core power [MW]", 1.0, 5000.0,
                                   float(ss.get("op_pmw", 160.0)))
            th_H = st.number_input(
                "active height [m]", 0.5, 6.0,
                (float(np.count_nonzero(res_t["axial_power"] > 1e-12)
                       * (res_t["z"][1] - res_t["z"][0]) / 100.0)
                 if res_t.get("nz", 1) > 1 else 2.0))
        with t2:
            th_m = st.number_input("total core flow [kg/s]", 10.0,
                                   30000.0, 37 * 70.0)
            th_tin = st.number_input("inlet temperature [°C]", 200.0,
                                     330.0, 258.0)
        with t3:
            th_np = st.number_input("pins / assembly", 10, 400, 264)
            th_d = st.number_input("pin OD [mm]", 6.0, 15.0, 9.5)
            th_pr = st.number_input("system pressure [MPa]", 5.0, 18.0,
                                    15.5, 0.5,
                                    help="sets the saturation temperature "
                                         "used for the boiling margin")
        with t4:
            th_pp = st.number_input("pin pitch [mm]", 7.0, 20.0, 12.6)
            # a symmetric chopped cosine cannot exceed peak/mean = pi/2;
            # cap the input honestly instead of silently clamping
            th_fz = st.number_input(
                "axial peaking F_z (2-D cores)", 1.0, 1.57,
                min(float(res_t.get("fz", 1.4)), 1.57),
                disabled=res_t.get("nz", 1) > 1,
                help="chopped-cosine ceiling ≈ 1.571; for higher axial "
                     "peaking use a 3-D core (real axial profile)")
        st.caption(f"hot-assembly radial peaking from neutronics: "
                   f"F_ΔH = {fdh:.3f}  ·  T_sat @ {th_pr:.1f} MPa = "
                   f"{thermal.t_sat(float(th_pr)):.1f} °C")
        if st.button("🌡  Run channel analysis", type="primary"):
            nfa = int(np.isfinite(bp_t).sum())
            th_res = thermal.core_th(
                float(th_P), nfa, fdh, float(th_H),
                axial_profile=(res_t["axial_power"]
                               if res_t.get("nz", 1) > 1 else None),
                mflow_total_kg_s=float(th_m), t_in_C=float(th_tin),
                n_pins=int(th_np), pin_d_mm=float(th_d),
                pin_pitch_mm=float(th_pp), fz=float(th_fz),
                pressure_MPa=float(th_pr))
            ss.th = th_res

    if "th" in ss:
        th = ss.th
        a, h = th["avg"], th["hot"]
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("outlet T (avg / hot)",
                  f"{a['T_out']:.1f} / {h['T_out']:.1f} °C")
        c2.metric("max clad surface (hot)", f"{h['T_clad_max']:.1f} °C")
        c3.metric("max fuel centreline (hot)", f"{h['T_fuel_max']:.0f} °C")
        c4.metric("saturation margin (hot)",
                  f"{h['sat_margin']:+.1f} K",
                  delta_color="normal" if h["sat_margin"] > 0 else "inverse")
        c5, c6, c7, c8 = st.columns(4)
        c5.metric("core ΔP (friction)", f"{a['dP_kPa']:.1f} kPa")
        c6.metric("coolant velocity", f"{a['v_cool']:.2f} m/s")
        c7.metric("Re", f"{a['Re']:.2e}")
        c8.metric("energy closure",
                  f"{th['Q_balance_MW']:.1f} MW")
        if h["sat_margin"] <= 0:
            st.error("⚠️ Clad surface exceeds saturation — the "
                     "single-phase model is out of its validity range "
                     "(boiling would start). Increase flow or reduce "
                     "power.")
        g1, g2 = st.columns(2)
        with g1:
            st.plotly_chart(plots.traverse_fig(
                a["z"], [("avg channel", a["T_cool"]),
                         ("hot channel", h["T_cool"]),
                         ("hot clad surface", h["T_clad"])],
                title="Coolant & clad temperatures", xlabel="z [m]",
                ylabel="T [°C]"), width="stretch")
        with g2:
            st.plotly_chart(plots.traverse_fig(
                a["z"], [("avg channel", a["T_fuel"]),
                         ("hot channel", h["T_fuel"])],
                title="Fuel centreline temperature", xlabel="z [m]",
                ylabel="T [°C]"), width="stretch")
        if h.get("fz_clamped"):
            st.warning(f"⚠️ Requested F_z exceeded the symmetric "
                       f"chopped-cosine ceiling ({h['fz_max']:.3f}) and "
                       f"was clamped — peak temperatures here are "
                       f"UNDER-predicted. Use a 3-D core for higher "
                       f"axial peaking.")
        st.caption(f"Single-phase model; constant water properties; "
                   f"saturation margin uses T_sat = {h.get('t_sat',344.8):.1f} "
                   f"°C at {h.get('p_MPa',15.5):.1f} MPa; uniform flow split "
                   f"between assemblies; UO₂ k=3 W/m·K, gap 6 kW/m²·K, "
                   f"Zry k=17 W/m·K. Energy conservation exact; formulas "
                   f"hand-verified in verify.py. Boiling/DNB not modelled "
                   f"— the saturation margin flags the limit.")

# ----------------------------------------------------------------------
# 🔥 Burnup
# ----------------------------------------------------------------------
with tabs[6]:
    st.header("Fuel-cycle simulation (block-wise depletion)")
    st.write("Quasi-static burnup: the Fortran engine solves the flux, the "
             "nuclide chain (U-235/238 → Pu-239/240/241, equilibrium "
             "Xe-135/Sm-149, lumped FPs) is integrated per block, cross "
             "sections are rebuilt, and the cycle marches on. Optionally "
             "with a **critical-boron letdown** each step.")
    if IS3D:
        st.info("Burnup runs on the 2-D (radial) model in this version — "
                "load a 2-D core (e.g. Designer PWR). 3-D depletion is on "
                "the roadmap.")
    n_des_b = sum(1 for m in ss.materials
                  if "designer" in m and
                  m["designer"].get("N", {}).get("U238", 0) > 0)
    if n_des_b == 0:
        st.info("Depletion needs designer-made fuels — load the "
                "**Designer PWR** preset, build fuels in the 🧬 tab, or "
                "attach an equivalent fuel to a benchmark material.")
    elif not IS3D:
        b_auto = st.checkbox(
            "🔚 auto end-of-cycle — deplete until the core can no longer "
            "stay critical (the physics sets the cycle length)", True)
        bc1, bc2, bc3, bc4 = st.columns(4)
        with bc1:
            b_sp = st.number_input("specific power [W/gU]", 5.0, 60.0, 38.0)
        with bc2:
            b_target = st.number_input("cycle burnup [MWd/kgU]", 1.0, 60.0,
                                       15.0, 1.0, disabled=b_auto)
        with bc3:
            b_step = st.number_input("step [MWd/kgU]", 0.5, 5.0, 3.0, 0.5)
        with bc4:
            b_mode = st.radio("boron", ["fixed (as built)",
                                        "critical letdown"], index=1)
        if st.button("🔥  Run fuel cycle", type="primary",
                     disabled=eng is None):
            horizon = 60.0 if b_auto else b_target
            nsteps = int(horizon / b_step) + 1
            bprog = st.progress(0.0, text="starting…")

            def bucb(n, row):
                bprog.progress(min((n + 1) / (nsteps + 1), 1.0),
                               text=f"BU={row['bu']:.1f} MWd/kgU  "
                                    f"EFPD={row['efpd']:.0f}  "
                                    f"k={row['keff']:.5f}  "
                                    f"boron={row['ppm']:.0f} ppm")

            try:
                bres = burnup.deplete_core(
                    current_cfg(), specific_power=b_sp,
                    bu_target=(None if b_auto else b_target), dbu=b_step,
                    letdown=(b_mode == "critical letdown"),
                    threads=threads, callback=bucb)
                ss.burn = (bres, current_cfg())
                ss.burn_params = dict(
                    specific_power=float(b_sp),
                    bu_target=(None if b_auto else float(b_target)),
                    dbu=float(b_step),
                    letdown=(b_mode == "critical letdown"))
                bprog.progress(1.0, text="cycle complete")
            except Exception as e:
                st.error(f"Burnup failed:\n\n```\n{e}\n```")

    if "burn" in ss:
        bres, bcfg = ss.burn
        hist = bres["history"]
        if bres.get("eoc_reason"):
            st.success(f"**End of cycle ({bres['eoc_reason']})** — "
                       f"achievable cycle burnup "
                       f"**{hist[-1]['bu']:.1f} MWd/kgU** = "
                       f"**{bres['efpd']:.0f} EFPD** "
                       f"({bres['efpd']/365.25:.2f} full-power years)")
        hm1, hm2, hm3, hm4 = st.columns(4)
        hm1.metric("xenon worth (BOC)",
                   f"{bres['maps']['xenon_worth_pcm']:+,.0f} pcm")
        hm2.metric("EOC k_eff", f"{hist[-1]['keff']:.5f}")
        hm3.metric("EOC boron", f"{hist[-1]['ppm']:.0f} ppm")
        hm4.metric("F_xy BOC → EOC",
                   f"{hist[0]['fxy']:.2f} → {hist[-1]['fxy']:.2f}")
        cA, cB = st.columns(2)
        bus = [h["bu"] for h in hist]
        with cA:
            st.plotly_chart(plots.search_fig(
                list(zip(bus, [h["keff"] for h in hist])), 1.0,
                title="k_eff over the cycle", xlabel="burnup [MWd/kgU]"),
                width="stretch")
        with cB:
            st.plotly_chart(plots.traverse_fig(
                bus, [("boron [ppm]", [h["ppm"] for h in hist])],
                title="Boron letdown curve", xlabel="burnup [MWd/kgU]",
                ylabel="boron [ppm]"), width="stretch")
        cC, cD = st.columns(2)
        with cC:
            st.plotly_chart(plots.block_power_fig(
                bres["bumap"], bcfg["pitch"],
                title="EOC burnup map [MWd/kgU]", zmid=None,
                colorscale="Oranges", fmt=1, unit="BU"), width="stretch")
        with cD:
            bp_boc = runner.block_powers(bres["maps"]["res_boc"], bcfg)
            bp_eol = runner.block_powers(bres["maps"]["res_eol"], bcfg)
            st.plotly_chart(plots.block_power_fig(
                bp_eol - bp_boc, bcfg["pitch"],
                title="Power shift EOC − BOC (flattening)", zmid=0.0,
                fmt=2, unit="ΔP/P̄"), width="stretch")
        hdf = pd.DataFrame(hist)
        st.dataframe(hdf, hide_index=True, width="stretch")
        st.download_button("⬇ cycle_history.csv",
                           hdf.to_csv(index=False), "cycle_history.csv")

        with st.expander("🔁 Multi-cycle: reload & run the NEXT cycle"):
            st.write("Batch reload, the core-designer way: the most-burned "
                     "fraction of the fuel blocks is replaced by fresh "
                     "assemblies (their original composition); the rest "
                     "carry their end-of-cycle isotopics and cumulative "
                     "burnup into the next cycle.")
            rl1, rl2 = st.columns([1, 1])
            with rl1:
                rl_frac = st.slider("fresh batch fraction", 0.1, 0.9,
                                    0.34, 0.02)
            with rl2:
                st.write("")
                go_nc = st.button("🔁 Run next cycle", type="primary",
                                  disabled=eng is None)
            if go_nc:
                bp2 = st.progress(0.0, text="next cycle…")

                def ncb(n, row):
                    bp2.progress(min((n + 1) / 22, 1.0),
                                 text=f"BU={row['bu']:.1f}  "
                                      f"k={row['keff']:.5f}  "
                                      f"boron={row['ppm']:.0f} ppm")
                try:
                    init = burnup.reload_core(bres["state"],
                                              float(rl_frac))
                    bp_prm = ss.get("burn_params") or dict(
                        specific_power=38.0, bu_target=None, dbu=3.0,
                        letdown=True)
                    bres2 = burnup.deplete_core(
                        bcfg, threads=threads, callback=ncb,
                        init_state=init, **bp_prm)
                    ss.setdefault("cycles", []).append(dict(
                        cycle=len(ss.get("cycles", [])) + 2,
                        fresh_frac=float(rl_frac),
                        cycle_BU=round(bres2["history"][-1]["bu"], 1),
                        EFPD=round(bres2["efpd"], 0),
                        BOC_ppm=round(bres2["history"][0]["ppm"], 0),
                        EOC_reason=bres2.get("eoc_reason") or "target"))
                    ss.burn = (bres2, bcfg)
                    bp2.progress(1.0, text="cycle complete")
                    st.rerun()
                except Exception as e:
                    st.error(str(e))
            if ss.get("cycles"):
                st.dataframe(pd.DataFrame(ss.cycles), hide_index=True,
                             width="stretch")
                st.caption("Cycle 1 is the run above; each row here is a "
                           "subsequent reload cycle. Watch the cycle "
                           "length approach the equilibrium-cycle value "
                           "as the batch pattern converges.")

# ----------------------------------------------------------------------
# ⏱ Transient (point kinetics)
# ----------------------------------------------------------------------
with tabs[7]:
    st.header("Time-dependent power — point kinetics & accidents")
    st.write("Reactivity events and full **design-basis accident "
             "sequences** solved with **6-delayed-group point kinetics** "
             "and an **exact matrix-exponential** step, with optional "
             "**√T Doppler + moderator (MTC)** feedback and an automatic "
             "**reactor-protection trip → scram**. Validated against the "
             "analytic **inhour equation** and an accident-physics suite "
             "in verify.py. Reactivities come straight from the worth "
             "tools (e.g. CRA bank worth).")

    _PRESETS = {"custom (manual)": None,
                "REA — rod ejection (protected)": "rea",
                "ATWS — rod ejection, scram fails": "atws",
                "rod withdrawal (protected ramp)": "rod_withdrawal"}
    t_preset = st.selectbox(
        "scenario", list(_PRESETS), index=0,
        help="Named presets load canonical design-basis accidents; "
             "'custom' exposes every control below.")
    is_custom = _PRESETS[t_preset] is None

    tc1, tc2, tc3, tc4 = st.columns(4)
    with tc1:
        t_P0 = st.number_input("initial power P₀ [MW]", 0.1, 5000.0, 160.0)
        t_Lam = st.number_input("Λ prompt gen. time [s]", 1e-6, 1e-3,
                                2e-5, format="%.1e")
    with tc2:
        t_kind = st.radio("event", ["step", "ramp", "scram"], index=0,
                          disabled=not is_custom)
        t_rho = st.number_input(
            "external ρ [pcm]", -10000.0, 2000.0, 150.0, 10.0,
            disabled=not is_custom,
            help="β = 650 pcm; ρ > β is prompt-critical — allowed here "
                 "for rod-ejection studies")
    with tc3:
        t_aux = st.number_input(
            "ramp time / scram delay [s]", 0.0, 120.0, 5.0,
            disabled=(not is_custom or t_kind == "step"))
        t_end = st.number_input("simulate for [s]", 0.1, 900.0, 15.0)
    with tc4:
        t_fb = st.checkbox("feedback", True, disabled=not is_custom)
        t_alpha = st.number_input("Doppler α_D [pcm/K]", -10.0, 0.0, -3.0,
                                  0.1, disabled=not (is_custom and t_fb))
        t_sqrt = st.checkbox("√T Doppler", True,
                             disabled=not (is_custom and t_fb))
    with st.expander("feedback & reactor protection (custom)"):
        g1, g2, g3 = st.columns(3)
        t_mtc = g1.number_input("moderator α_M / MTC [pcm/K]", -80.0, 5.0,
                                -15.0, disabled=not is_custom)
        t_mcpf = g2.number_input("fuel heat cap. [MJ/K]", 0.5, 100.0, 6.0,
                                 disabled=not is_custom)
        t_mcpm = g3.number_input("moderator heat cap. [MJ/K]", 0.5, 200.0,
                                 12.0, disabled=not is_custom)
        h1, h2, h3 = st.columns(3)
        t_trip_on = h1.checkbox("reactor trip → scram", True,
                                disabled=not is_custom)
        t_setpt = h2.number_input("high-flux trip [× nominal]", 1.01, 5.0,
                                  1.18,
                                  disabled=not (is_custom and t_trip_on))
        t_scrw = h3.number_input("scram-bank worth [pcm]", -12000.0,
                                 -500.0, -6000.0,
                                 disabled=not (is_custom and t_trip_on))
        t_fmass = st.number_input(
            "fuel mass [kg]  (0 = skip cal/g)", 0.0, 1e5, 0.0,
            disabled=not is_custom,
            help="If set, the peak specific fuel-enthalpy rise is "
                 "reported in cal/g — the RIA safety metric.")

    if st.button("⏱  Run transient / accident", type="primary"):
        if not is_custom:
            cfg = kinetics.accident_preset(_PRESETS[t_preset],
                                           P0_MW=float(t_P0))
            with st.spinner("integrating accident sequence…"):
                sim = kinetics.simulate(P0_MW=float(t_P0),
                                        Lambda=float(t_Lam), **cfg)
            ss.transient = (sim, cfg["scenario"], True)
        else:
            scen = {"type": t_kind, "rho_pcm": float(t_rho)}
            if t_kind == "ramp":
                scen["t_ramp"] = float(t_aux)
            if t_kind == "scram":
                scen["t_delay"] = float(t_aux)
            fb = None
            if t_fb:
                fb = dict(alpha_pcm_K=float(t_alpha), mcp_MJ_K=float(t_mcpf),
                          alpha_mod_pcm_K=float(t_mtc),
                          mcp_mod_MJ_K=float(t_mcpm),
                          tau_fm=6.0, tau_ms=3.0, T_sink=565.0)
                if t_sqrt:
                    fb["doppler_mode"] = "sqrt"
            trip = (dict(power_frac=float(t_setpt), delay_s=0.5,
                         scram_rho_pcm=float(t_scrw), scram_ramp_s=2.0)
                    if t_trip_on else None)
            with st.spinner("integrating kinetics…"):
                sim = kinetics.simulate(
                    scen, t_end=float(t_end), dt=1e-3, P0_MW=float(t_P0),
                    Lambda=float(t_Lam), feedback=fb, trip=trip,
                    fuel_mass_kg=(float(t_fmass) or None))
            ss.transient = (sim, scen, t_fb)

    if "transient" in ss:
        sim, scen, had_fb = ss.transient
        m = st.columns(4)
        m[0].metric("peak power", f"{sim['peak_MW']:.1f} MW",
                    f"{sim['peak_MW']/sim['P_MW'][0]:.0f}× nom"
                    if sim["P_MW"][0] > 0 else None)
        m[1].metric("power at end", f"{sim['final_MW']:.2f} MW")
        if sim.get("t_trip") is not None:
            m[2].metric("reactor trip", f"{sim['t_trip']:.2f} s")
        elif "period_inhour_s" in sim:
            m[2].metric("asymptotic period", f"{sim['period_inhour_s']:.1f} s")
        else:
            m[2].metric("net ρ at end", f"{sim['rho_pcm'][-1]:+.0f} pcm")
        m[3].metric("energy released", f"{sim['energy_MJ']:.0f} MJ")
        s = st.columns(4)
        s[0].metric("external reactivity", f"{sim['rho_dollars_ext']:+.2f} $")
        s[1].metric("prompt-critical",
                    "YES ⚠️" if sim["prompt_critical"] else "no")
        if had_fb:
            s[2].metric("peak fuel ΔT", f"{sim['dT_fuel_peak']:.0f} K")
        if "enthalpy_cal_g" in sim:
            s[3].metric("peak fuel enthalpy",
                        f"{sim['enthalpy_cal_g']:.0f} cal/g")
        pa, pb = st.columns(2)
        with pa:
            st.plotly_chart(plots.traverse_fig(
                sim["t"], [("P [MW]", sim["P_MW"])],
                title="Core power vs time", xlabel="t [s]",
                ylabel="P [MW]"), width="stretch")
        with pb:
            series = [("net ρ [pcm]", sim["rho_pcm"])]
            if had_fb:
                series.append(("T_fuel − T₀ [K]",
                               sim["T_fuel"] - sim["T_fuel"][0]))
                if "T_mod" in sim:
                    series.append(("T_mod − T₀ [K]",
                                   sim["T_mod"] - sim["T_mod"][0]))
            st.plotly_chart(plots.traverse_fig(
                sim["t"], series, title="Reactivity & temperatures",
                xlabel="t [s]", ylabel=""), width="stretch")
        st.caption("Point-kinetics approximation (spatial flux shape "
                   "frozen — space-time kinetics is on the roadmap). "
                   "6-group U-235 delayed data, β = 650 pcm. Automatic "
                   "scram = statically-computed bank worth inserted over "
                   "the drop time; decay heat is NOT modelled. RIA cal/g "
                   "uses the lumped fuel node (integral estimate, not a "
                   "radial enthalpy profile).")

    st.divider()
    st.subheader("☁️ Xenon-135 transient (iodine pit / load-follow)")
    st.write("Xe-135 reactivity after a power-level change — the classic "
             "**iodine pit** after shutdown and the transient peak on "
             "load reduction. I-135/Xe-135 balance solved with exact "
             "per-step linear updates; the reactivity scale is anchored "
             "to the equilibrium-xenon worth.")
    x1, x2, x3, x4 = st.columns(4)
    with x1:
        xe_pf = st.number_input("final power [% of nominal]", 0.0, 100.0,
                                0.0, 5.0) / 100.0
    with x2:
        xe_sp = st.number_input("σφ at 100% [10⁻⁵/s]", 1.0, 20.0, 5.3,
                                help="2.65 Mb × φ_thermal; PWR ~3–12")
    with x3:
        xe_rho = st.number_input("equilibrium Xe worth [pcm]", -6000.0,
                                 -500.0, -2800.0, 100.0)
    with x4:
        xe_t = st.number_input("simulate [h]", 12.0, 200.0, 72.0)
    if st.button("☁️  Run xenon transient"):
        ss.xet = kinetics.xenon_transient(xe_pf, float(xe_t),
                                          sigphi0=float(xe_sp) * 1e-5,
                                          rho_eq_pcm=float(xe_rho))
    if "xet" in ss:
        xe = ss.xet
        xc1, xc2, xc3 = st.columns(3)
        xc1.metric("pit depth", f"{xe['pit_rho_pcm']:.0f} pcm "
                                f"(×{xe['pit_ratio']:.2f})")
        xc2.metric("pit time", f"{xe['pit_time_h']:.1f} h")
        xc3.metric("back to pre-event worth",
                   f"{xe['recover_h']:.1f} h" if xe.get("recover_h")
                   else "—")
        st.plotly_chart(plots.traverse_fig(
            xe["t_h"], [("ρ_Xe [pcm]", xe["rho_pcm"])],
            title="Xenon reactivity vs time", xlabel="t [h]",
            ylabel="ρ_Xe [pcm]"), width="stretch")
        st.caption("Restart inside the pit needs the pit-depth worth of "
                   "excess reactivity — the practical meaning of the "
                   "'restart window'. 0-D model; scale anchored to the "
                   "stated equilibrium worth.")

# ----------------------------------------------------------------------
# 🔧 Physics tools
# ----------------------------------------------------------------------
with tabs[8]:
    st.header("Physics tools")
    st.caption("Each tool runs the Fortran engine on the *current* core "
               "configuration (possibly several times).")

    st.subheader("1 · Material-swap reactivity worth")
    st.write("Replace one material by another everywhere in the map "
             "(including all axial zones in 3-D) and measure the "
             "reactivity change — e.g. control-rod worth.")
    names = [f"{i+1} · {m['name']}" for i, m in enumerate(ss.materials)]
    w1, w2, w3 = st.columns([1, 1, 1])
    with w1:
        ma = st.selectbox("replace material", names, index=min(2, len(names)-1))
    with w2:
        mb = st.selectbox("with material", names, index=min(1, len(names)-1))
    with w3:
        st.write("")
        go_w = st.button("Compute worth", disabled=eng is None)
    if go_w:
        ia, ib = names.index(ma) + 1, names.index(mb) + 1
        if ia == ib:
            st.warning("Pick two different materials.")
        else:
            with st.spinner("two eigenvalue runs…"):
                w = runner.material_swap_worth(current_cfg(), ia, ib,
                                               threads=threads)
            c1, c2, c3 = st.columns(3)
            c1.metric("k (base)", f"{w['k1']:.6f}")
            c2.metric("k (swapped)", f"{w['k2']:.6f}")
            c3.metric("worth Δρ", f"{w['worth_pcm']:+,.0f} pcm")
            st.caption("Δρ = ρ(swapped) − ρ(base),  ρ = (k−1)/k")

    if ss.get("rod_meta") and IS3D:
        st.divider()
        st.subheader("2 · Rod insertion sweep (3-D only) — the S-curve")
        st.write("Vary the insertion depth of **rod 5** from fully "
                 "withdrawn to fully inserted and watch the classic "
                 "S-shaped worth curve — differential worth peaks where "
                 "the axial flux peaks. Physics a 2-D model cannot show.")
        nstep_rod = st.slider("depths to compute", 3, 11, 6)
        if st.button("Sweep rod insertion", type="primary",
                     disabled=eng is None):
            core_h = float(ss.rod_meta.get("core_h", 340.0))
            depths = np.linspace(0.0, core_h, nstep_rod)
            prog = st.progress(0.0, text="sweeping…")
            hist_rod = []
            try:
                for irun, dpt in enumerate(depths):
                    cfg = current_cfg()
                    cfg["axial"] = dict(
                        presets.iaea3d_axial(float(dpt),
                                             dz=ss.axial["dz"]),
                        divz=int(ss.divz))
                    r = runner.run_case(cfg, threads=threads)
                    hist_rod.append((float(dpt), r["keff"]))
                    prog.progress((irun + 1) / len(depths),
                                  text=f"depth {dpt:.0f} cm → "
                                       f"k = {r['keff']:.6f}")
                worth = (runner.rho_pcm(hist_rod[-1][1]) -
                         runner.rho_pcm(hist_rod[0][1]))
                c1, c2 = st.columns(2)
                c1.metric("total rod-5 worth", f"{worth:+,.0f} pcm")
                c2.metric("k (out → in)",
                          f"{hist_rod[0][1]:.5f} → {hist_rod[-1][1]:.5f}")
                st.plotly_chart(plots.search_fig(
                    hist_rod, 1.0, title="Rod-5 insertion S-curve",
                    xlabel="insertion depth from top of core [cm]"),
                    width="stretch")
            except Exception as e:
                st.error(str(e))

        st.markdown("**Critical rod position** — bisect the insertion "
                    "depth until k_eff hits a target (the 3-D analogue "
                    "of a boron search).")
        rc1, rc2 = st.columns([1, 1])
        with rc1:
            rk_tgt = st.number_input("target k_eff", 0.95, 1.10, 1.025,
                                     format="%.4f", key="rodtgt")
        with rc2:
            st.write("")
            go_rp = st.button("Find rod position", disabled=eng is None)
        if go_rp:
            core_h = float(ss.rod_meta.get("core_h", 340.0))

            def k_at(depth):
                c = current_cfg()
                c["axial"] = dict(
                    presets.iaea3d_axial(float(depth),
                                         dz=ss.axial["dz"]),
                    divz=int(ss.divz))
                return runner.run_case(c, threads=threads)["keff"]

            prog2 = st.progress(0.0, text="bisecting…")
            lo_d, hi_d = 0.0, core_h
            k_lo, k_hi = k_at(lo_d), k_at(hi_d)     # k decreases w/ depth
            if not (min(k_lo, k_hi) <= rk_tgt <= max(k_lo, k_hi)):
                st.error(f"target outside the achievable range "
                         f"[{k_hi:.5f}, {k_lo:.5f}]")
            else:
                d_mid, k_mid = lo_d, k_lo
                for it2 in range(10):
                    d_mid = 0.5 * (lo_d + hi_d)
                    k_mid = k_at(d_mid)
                    prog2.progress((it2 + 1) / 10,
                                   text=f"depth {d_mid:.0f} cm → "
                                        f"k = {k_mid:.6f}")
                    if abs(runner.rho_pcm(k_mid) -
                           runner.rho_pcm(rk_tgt)) < 3.0:
                        break
                    if k_mid > rk_tgt:
                        lo_d = d_mid
                    else:
                        hi_d = d_mid
                cc1, cc2 = st.columns(2)
                cc1.metric("critical rod position",
                           f"{d_mid:.0f} cm inserted")
                cc2.metric("k at that depth", f"{k_mid:.6f}")
                st.caption(f"depth resolution = one axial layer "
                           f"(dz/divz = {ss.axial['dz']/ss.divz:.0f} cm) "
                           f"— refine the axial mesh for finer rod "
                           f"positioning.")

    st.divider()
    st.subheader("3 · Critical boron search [ppm]")
    n_des = sum(1 for m in ss.materials if "designer" in m)
    st.write("Bisection on the **soluble-boron concentration**: every "
             "designer-made material (🧬 tab) is regenerated at the trial "
             "ppm until k_eff hits the target — a real critical-boron "
             "search, like a core simulator.")
    if n_des == 0:
        st.info("No designer materials in this core. Load the "
                "**Designer PWR** preset or add fuels from the 🧬 tab.")
    else:
        b1, b2 = st.columns([1, 1])
        with b1:
            b_target2 = st.number_input("target k", 0.9, 1.1, 1.0,
                                        format="%.4f", key="btgt")
        with b2:
            st.write("")
            go_b = st.button("Search critical boron",
                             disabled=eng is None, type="primary")
        if go_b:
            bprog2 = st.progress(0.0, text="searching…")

            def bcb(nn, ppm, kk):
                bprog2.progress(min(nn / 20, 1.0),
                                text=f"run {nn}:  {ppm:.0f} ppm → "
                                     f"k = {kk:.6f}")

            try:
                sb = runner.critical_boron_search(
                    current_cfg(), target=float(b_target2),
                    threads=threads, callback=bcb)
                bprog2.progress(1.0, text="done")
                c1, c2, c3 = st.columns(3)
                c1.metric("critical boron", f"{sb['ppm']:.0f} ppm")
                c2.metric("final k_eff", f"{sb['k']:.6f}")
                c3.metric("engine runs", sb["runs"])
                if sb.get("note"):
                    st.warning(sb["note"])
                st.plotly_chart(plots.search_fig(
                    sb["history"], float(b_target2),
                    title="Critical boron search",
                    xlabel="soluble boron [ppm]"), width="stretch")
            except Exception as e:
                st.error(str(e))

    st.divider()
    st.subheader("4 · Criticality search (generic ΔΣa)")
    st.write("Bisection on a uniform ΔΣa added to the chosen group of the "
             "chosen materials until k_eff hits the target.")
    fis_ids = runner.fissile_ids(ss.materials)
    s1, s2, s3, s4 = st.columns([2, 1, 1, 1])
    with s1:
        tgt_mats = st.multiselect("materials to poison", names,
                                  default=[names[i - 1] for i in fis_ids])
    with s2:
        grp = st.number_input("group", 1, ss.ng, ss.ng)
    with s3:
        ktarget = st.number_input("target k", 0.5, 1.5, 1.0, format="%.4f")
    with s4:
        st.write("")
        go_s = st.button("Search", disabled=eng is None or not tgt_mats)
    if go_s:
        ids = [names.index(n) + 1 for n in tgt_mats]
        prog = st.progress(0.0, text="searching…")

        def cb(n, dsa, kk):
            prog.progress(min(n / 20, 1.0),
                          text=f"run {n}:  ΔΣa = {dsa:+.4e}  →  k = {kk:.6f}")

        try:
            s = runner.criticality_search(current_cfg(), ids, group=int(grp),
                                          target=float(ktarget),
                                          threads=threads, callback=cb)
            prog.progress(1.0, text="done")
            c1, c2, c3 = st.columns(3)
            c1.metric("critical ΔΣa", f"{s['dsa']:+.4e} cm⁻¹")
            c2.metric("final k_eff", f"{s['k']:.6f}")
            c3.metric("engine runs", s["runs"])
            st.plotly_chart(plots.search_fig(s["history"], float(ktarget)),
                            width="stretch")
        except Exception as e:
            st.error(str(e))

st.divider()
st.caption("CoreForge v8.4 — SMR/MMR neutronics design & analysis code "
           "system · Fortran 2-D/3-D multigroup diffusion (red-black SOR, "
           "OpenMP) · 3-D core builder + project save/load · IAEA-2D −1 "
           "pcm · IAEA-3D −7 pcm · C5G7 demo · fuel designer + inverse "
           "matching · burnup with letdown & auto-EOC · point-kinetics "
           "transients (inhour-validated) · operating-point absolute units "
           "· QA integrity guard (benchmark fingerprinting) · versioned "
           "HTML reports with traceability · Türkçe dokümantasyon: docs/ "
           "· MIT license")
