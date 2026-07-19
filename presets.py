"""CoreForge — bundled benchmark presets and the physics fingerprint.

Each preset is a complete core description on a *block* lattice:
`blocks` is a list of rows (row 0 = TOP of the core), each entry a
material id.  The UI/runner expands every block into `div x div`
fine-mesh cells (cell size = pitch/div), so the geometry stays exact
under refinement.  3-D presets add `axial` (bottom-to-top zones, each
with its own radial block map) and 6 boundary conditions.
"""

import hashlib
import json


def physics_fingerprint(ng, pitch, blocks, materials, bc, gamma, axial=None):
    """Short hash of everything that defines the PHYSICS of a case:
    groups, lattice pitch, block maps, cross sections, boundary
    conditions, vacuum gamma and the axial zone structure.

    Deliberately EXCLUDED: mesh refinement (div/divz) — refining the
    mesh is a legitimate convergence study and keeps a benchmark
    reference meaningful — and cosmetic fields (material names).
    Used by the QA layer: a benchmark's 'diff vs reference' is only
    reported while the loaded preset's fingerprint still matches."""
    mats = [[[round(float(v), 12) for v in m["D"]],
             [round(float(v), 12) for v in m["Sa"]],
             [round(float(v), 12) for v in m["nuSf"]],
             [round(float(v), 12) for v in m["chi"]],
             [[round(float(v), 12) for v in row] for row in m["scat"]]]
            for m in materials]
    ax = None
    if axial:
        ax = [round(float(axial["dz"]), 9),
              [[z.get("label", ""), int(z["layers"]),
                [list(map(int, r)) for r in z["blocks"]]]
               for z in axial["zones"]]]
    payload = json.dumps([int(ng), round(float(pitch), 9),
                          [list(map(int, r)) for r in blocks],
                          mats, list(bc), round(float(gamma), 9), ax],
                         sort_keys=True)
    return hashlib.sha1(payload.encode()).hexdigest()[:16]


def fingerprint_of(p):
    """Fingerprint of a preset/config-like dict."""
    return physics_fingerprint(p["ng"], p["pitch"], p["blocks"],
                               p["materials"], p["bc"], p["gamma"],
                               p.get("axial"))


# ----------------------------------------------------------------------
# IAEA-2D PWR benchmark (ANL benchmark book 11-A2, "2D IAEA")
# Quarter core, 10 cm blocks (17x17), two groups, axial buckling
# Bg2 = 0.8e-4 folded into Sigma_a (Sa += D*Bg2) as the spec requires.
# Geometry transcribed from the official polygon description and
# rasterised exactly (every polygon edge is a multiple of 10 cm).
# Reference: k_eff = 1.02959.
# ----------------------------------------------------------------------
_BG2 = 0.8e-4

IAEA2D_BLOCKS = [
    [4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4],
    [4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4],
    [1,1,1,1,1,4,4,4,4,4,4,4,4,4,4,4,4],
    [1,1,1,1,1,4,4,4,4,4,4,4,4,4,4,4,4],
    [2,2,2,1,1,1,1,1,1,4,4,4,4,4,4,4,4],
    [2,2,2,1,1,1,1,1,1,4,4,4,4,4,4,4,4],
    [2,2,2,2,2,2,2,1,1,1,1,4,4,4,4,4,4],
    [2,2,2,2,2,2,2,1,1,1,1,4,4,4,4,4,4],
    [3,2,2,2,2,2,2,3,3,1,1,1,1,4,4,4,4],
    [3,2,2,2,2,2,2,3,3,1,1,1,1,4,4,4,4],
    [2,2,2,2,2,2,2,2,2,2,2,1,1,4,4,4,4],
    [2,2,2,2,2,2,2,2,2,2,2,1,1,4,4,4,4],
    [2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,4,4],
    [2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,4,4],
    [2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,4,4],
    [2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,4,4],
    [3,2,2,2,2,2,2,3,3,2,2,2,2,1,1,4,4],
]

IAEA2D_MATERIALS = [
    dict(name="Fuel 1 (outer)",
         D=[1.5, 0.4], Sa=[0.010 + 1.5*_BG2, 0.080 + 0.4*_BG2],
         nuSf=[0.0, 0.135], chi=[1.0, 0.0], scat=[[0.0, 0.02], [0.0, 0.0]]),
    dict(name="Fuel 2 (inner)",
         D=[1.5, 0.4], Sa=[0.010 + 1.5*_BG2, 0.085 + 0.4*_BG2],
         nuSf=[0.0, 0.135], chi=[1.0, 0.0], scat=[[0.0, 0.02], [0.0, 0.0]]),
    dict(name="Fuel 2 + control rod",
         D=[1.5, 0.4], Sa=[0.010 + 1.5*_BG2, 0.130 + 0.4*_BG2],
         nuSf=[0.0, 0.135], chi=[1.0, 0.0], scat=[[0.0, 0.02], [0.0, 0.0]]),
    dict(name="Reflector (water)",
         D=[2.0, 0.3], Sa=[0.000 + 2.0*_BG2, 0.010 + 0.3*_BG2],
         nuSf=[0.0, 0.0], chi=[0.0, 0.0], scat=[[0.0, 0.04], [0.0, 0.0]]),
]

# ----------------------------------------------------------------------
# Homogeneous 2-group infinite medium: analytic k_inf = 1.085714...
# ----------------------------------------------------------------------
HOMOG_MATERIALS = [
    dict(name="Homogeneous fuel",
         D=[1.4, 0.4], Sa=[0.01, 0.10],
         nuSf=[0.007, 0.13], chi=[1.0, 0.0], scat=[[0.0, 0.018], [0.0, 0.0]]),
]

# ----------------------------------------------------------------------
# Demo PWR quarter mini-core (uniform 20 cm blocks)
# ----------------------------------------------------------------------
MINICORE_BLOCKS = [
    [4,4,4,4,4,4,4,4,4],
    [1,1,4,4,4,4,4,4,4],
    [1,1,1,1,4,4,4,4,4],
    [2,2,1,1,1,1,4,4,4],
    [2,2,2,2,1,1,1,4,4],
    [3,2,2,2,2,1,1,4,4],
    [2,2,2,2,2,2,1,1,4],
    [2,2,3,2,2,2,1,1,4],
    [3,2,2,2,3,2,2,1,4],
]


def preset_iaea2d():
    return dict(
        key="iaea2d",
        title="IAEA-2D PWR benchmark (11-A2)",
        description=(
            "The classic 2-group quarter-core LWR benchmark (ANL benchmark "
            "book 11-A2). Two fuel zones, four control-rod regions, 20 cm "
            "water reflector, axial buckling 0.8e-4 folded into absorption. "
            "Reference k_eff = 1.02959."
        ),
        ng=2, pitch=10.0, div=5,
        blocks=[row[:] for row in IAEA2D_BLOCKS],
        materials=[dict(m, Sa=list(m["Sa"]), D=list(m["D"]), nuSf=list(m["nuSf"]),
                        chi=list(m["chi"]), scat=[r[:] for r in m["scat"]])
                   for m in IAEA2D_MATERIALS],
        bc=["reflective", "vacuum", "reflective", "vacuum"],   # W E S N
        gamma=0.4692,
        ref_keff=1.02959,
        ref_source="ANL-7416 benchmark book, problem 11-A2 (2D IAEA)",
    )


def preset_homogeneous():
    return dict(
        key="homog",
        title="Homogeneous k∞ verification",
        description=(
            "Single homogeneous 2-group material, all boundaries reflective. "
            "The solver must reproduce the analytic infinite-medium value "
            "k∞ = (νΣf1 + νΣf2·Σ12/Σa2)/(Σa1+Σ12) = 1.085714 exactly."
        ),
        ng=2, pitch=20.0, div=1,
        blocks=[[1]*8 for _ in range(8)],
        materials=[dict(HOMOG_MATERIALS[0], Sa=list(HOMOG_MATERIALS[0]["Sa"]))],
        bc=["reflective", "reflective", "reflective", "reflective"],
        gamma=0.4692,
        ref_keff=1.0857142857,
        ref_source="analytic infinite-medium two-group formula",
    )


def preset_minicore():
    return dict(
        key="minicore",
        title="Demo PWR quarter mini-core",
        description=(
            "A 9×9-block quarter mini-core (20 cm blocks) built from the "
            "IAEA-2D material set: two fuel zones, four control rods, water "
            "reflector. A playground for the core builder and the tools tab "
            "— no published reference."
        ),
        ng=2, pitch=20.0, div=5,
        blocks=[row[:] for row in MINICORE_BLOCKS],
        materials=[dict(m, Sa=list(m["Sa"]), D=list(m["D"]), nuSf=list(m["nuSf"]),
                        chi=list(m["chi"]), scat=[r[:] for r in m["scat"]])
                   for m in IAEA2D_MATERIALS],
        bc=["reflective", "vacuum", "reflective", "vacuum"],
        gamma=0.4692,
        ref_keff=None,
        ref_source=None,
    )


# ----------------------------------------------------------------------
# C5G7 (OECD/NEA) — 7-group MOX quarter core, *pin-cell diffusion demo*.
# ----------------------------------------------------------------------
def _c5g7_blocks():
    import c5g7_data as d
    W = 7
    rows = []
    for r in range(17):
        rows.append(list(d.UO2_MAP[r]) + list(d.MOX_MAP[r]) + [W] * 17)
    for r in range(17):
        rows.append(list(d.MOX_MAP[r]) + list(d.UO2_MAP[r]) + [W] * 17)
    for _ in range(17):
        rows.append([W] * 51)
    return rows


def preset_c5g7():
    import math

    import c5g7_data as d

    # Each 1.26 cm pin cell holds an r = 0.54 cm pin surrounded by
    # moderator -> homogenise every pin cell volume-weighted.
    f = math.pi * 0.54**2 / 1.26**2          # pin volume fraction = 0.5770
    names = ["UO₂ pin cell", "MOX 4.3% pin cell", "MOX 7.0% pin cell",
             "MOX 8.7% pin cell", "Fission-chamber cell", "Guide-tube cell",
             "Moderator (water)"]
    MOD = 6
    mats = []
    for m in range(7):
        w = 1.0 if m == MOD else f
        st = [w * d.ST[m][g] + (1 - w) * d.ST[MOD][g] for g in range(7)]
        nsf = [w * d.NSF[m][g] for g in range(7)]
        sc = [[w * d.SS[m][a][b] + (1 - w) * d.SS[MOD][a][b]
               for b in range(7)] for a in range(7)]
        D = [1.0 / (3.0 * st[g]) for g in range(7)]
        Sa = [st[g] - sum(sc[g]) for g in range(7)]
        mats.append(dict(name=names[m], D=D, Sa=Sa, nuSf=nsf,
                         chi=list(d.CHI[m]), scat=sc))
    return dict(
        key="c5g7",
        title="C5G7 MOX — pin-cell diffusion demo",
        description=(
            "The OECD/NEA C5G7 7-group MOX benchmark solved with **pin-cell "
            "diffusion** (D = 1/3Σt, volume-homogenised pin cells). The real "
            "problem is a transport benchmark (MCNP k_eff = 1.18655; my S_N "
            "solver: −182 pcm) — this demo shows how far diffusion gets on a "
            "strongly heterogeneous pin lattice, and why transport is needed."
        ),
        ng=7, pitch=1.26, div=1,
        blocks=_c5g7_blocks(),
        materials=mats,
        bc=["reflective", "vacuum", "vacuum", "reflective"],   # W E S N
        gamma=0.4692,
        ref_keff=1.18655,
        ref_source="MCNP transport reference (NEA/NSC/DOC(2003)16) — "
                   "diffusion is EXPECTED to deviate here",
    )


# ----------------------------------------------------------------------
# SMR-class PWR (NuScale-style 37-assembly core, full core).
# ----------------------------------------------------------------------
_SMR_BZ2 = 2.25e-4          # axial buckling for a short (H ~ 2 m) core

SMR_BLOCKS = [
    [5,5,5,5,5,5,5,5,5],
    [5,5,5,3,3,3,5,5,5],
    [5,5,3,2,2,2,3,5,5],
    [5,3,2,4,1,4,2,3,5],
    [5,3,2,1,4,1,2,3,5],
    [5,3,2,4,1,4,2,3,5],
    [5,5,3,2,2,2,3,5,5],
    [5,5,5,3,3,3,5,5,5],
    [5,5,5,5,5,5,5,5,5],
]

SMR_MATERIALS = [
    dict(name="Fuel 2.6% (inner)",
         D=[1.43, 0.37], Sa=[0.0100 + 1.43*_SMR_BZ2, 0.0820 + 0.37*_SMR_BZ2],
         nuSf=[0.0058, 0.1240], chi=[1.0, 0.0],
         scat=[[0.0, 0.0165], [0.0, 0.0]]),
    dict(name="Fuel 3.1% (middle)",
         D=[1.43, 0.37], Sa=[0.0102 + 1.43*_SMR_BZ2, 0.0860 + 0.37*_SMR_BZ2],
         nuSf=[0.0062, 0.1340], chi=[1.0, 0.0],
         scat=[[0.0, 0.0165], [0.0, 0.0]]),
    dict(name="Fuel 3.8% (periphery)",
         D=[1.43, 0.37], Sa=[0.0104 + 1.43*_SMR_BZ2, 0.0910 + 0.37*_SMR_BZ2],
         nuSf=[0.0068, 0.1460], chi=[1.0, 0.0],
         scat=[[0.0, 0.0165], [0.0, 0.0]]),
    dict(name="Fuel 3.1% + CRA (rods in)",
         D=[1.43, 0.37], Sa=[0.0162 + 1.43*_SMR_BZ2, 0.1260 + 0.37*_SMR_BZ2],
         nuSf=[0.0062, 0.1340], chi=[1.0, 0.0],
         scat=[[0.0, 0.0165], [0.0, 0.0]]),
    dict(name="Heavy reflector (SS/H₂O)",
         D=[1.30, 0.35], Sa=[0.0025 + 1.30*_SMR_BZ2, 0.0280 + 0.35*_SMR_BZ2],
         nuSf=[0.0, 0.0], chi=[0.0, 0.0],
         scat=[[0.0, 0.0250], [0.0, 0.0]]),
    dict(name="Water reflector (for swap studies)",
         D=[1.50, 0.30], Sa=[0.0004 + 1.50*_SMR_BZ2, 0.0110 + 0.30*_SMR_BZ2],
         nuSf=[0.0, 0.0], chi=[0.0, 0.0],
         scat=[[0.0, 0.0450], [0.0, 0.0]]),
]


def preset_smr():
    return dict(
        key="smr",
        title="SMR-class PWR (37-assembly, full core)",
        description=(
            "A NuScale-style small modular PWR: 37 assemblies "
            "(3-5-7-7-7-5-3), out-in three-zone loading, five CRAs "
            "(rods **in** as loaded), heavy SS/H₂O reflector, axial "
            "buckling for a 2 m core folded into Sa. Representative "
            "two-group constants — no public SMR diffusion benchmark "
            "exists, so no reference k is claimed. Try the tools: CRA "
            "bank worth (swap 4→2), heavy-vs-water reflector worth "
            "(swap 5→6), poison-to-critical search."
        ),
        ng=2, pitch=21.42, div=6,
        blocks=[row[:] for row in SMR_BLOCKS],
        materials=[dict(m, D=list(m["D"]), Sa=list(m["Sa"]), nuSf=list(m["nuSf"]),
                        chi=list(m["chi"]), scat=[r[:] for r in m["scat"]])
                   for m in SMR_MATERIALS],
        bc=["vacuum", "vacuum", "vacuum", "vacuum"],
        gamma=0.4692,
        ref_keff=None,
        ref_source=None,
    )


# ----------------------------------------------------------------------
# Fuel-designer demo core (materials generated by xslib.pincell_xs)
# ----------------------------------------------------------------------
DESIGNER_BLOCKS = [
    [4,4,4,4,4,4,4,4,4,4],
    [3,3,3,4,4,4,4,4,4,4],
    [3,2,3,3,3,4,4,4,4,4],
    [2,2,2,2,3,3,4,4,4,4],
    [2,1,1,2,2,3,3,4,4,4],
    [1,1,1,1,2,2,3,3,4,4],
    [1,1,1,1,1,2,2,3,3,4],
    [1,1,1,1,1,2,2,3,3,4],
    [1,1,1,1,1,1,2,2,3,4],
    [1,1,1,1,1,1,2,2,3,4],
]


def preset_designer(ppm=1000.0):
    import xslib
    mats = [
        xslib.pincell_xs(1.8, ppm, name="UO₂ 1.8% (inner)"),
        xslib.pincell_xs(2.4, ppm, name="UO₂ 2.4% (middle)"),
        xslib.pincell_xs(3.0, ppm, name="UO₂ 3.0% (periphery)"),
        xslib.pincell_xs(0.0, ppm, r_fuel=0.0, name="Borated water refl."),
    ]
    return dict(
        key="designer",
        title="Designer PWR (physics-generated XS)",
        description=(
            "A quarter PWR whose cross sections are **generated from "
            "physical inputs** — enrichment (1.8/2.4/3.0 %), soluble boron "
            f"({ppm:.0f} ppm), moderator density — by the built-in "
            "two-group pin-cell model. Change the fuel in the 🧬 tab, run "
            "the **critical boron search** or a full **burnup cycle**."
        ),
        ng=2, pitch=21.42, div=6,
        blocks=[row[:] for row in DESIGNER_BLOCKS],
        materials=mats,
        bc=["reflective", "vacuum", "reflective", "vacuum"],
        gamma=0.4692,
        ref_keff=None, ref_source=None,
    )


# ----------------------------------------------------------------------
# IAEA-3D PWR benchmark (ANL benchmark book, problem 11 / "3D IAEA").
#
# Full x-y-z quarter core: 170x170 cm radial (the SAME validated 17x17
# block lattice as IAEA-2D), 380 cm high: 20 cm lower reflector, 340 cm
# core, 20 cm upper reflector.  Four control rods are fully inserted
# (they exist at every core elevation — exactly the 2-D map); a FIFTH
# rod at (30-50, 30-50) cm is inserted only 80 cm into the core from
# the top; in the upper reflector every rod position is 'reflrod'.
# Geometry transcribed from the official description (box coordinates
# cross-checked against the FeenoX example).  NO axial buckling here —
# the axial dimension is explicit.  Reference: k_eff = 1.02903.
# ----------------------------------------------------------------------
_ROD_BOXES_FULL = [(0, 10, 0, 10), (70, 90, 0, 10),
                   (0, 10, 70, 90), (70, 90, 70, 90)]
_ROD5_BOX = (30, 50, 30, 50)

IAEA3D_MATERIALS = [
    dict(name="Fuel 1 (outer)",
         D=[1.5, 0.4], Sa=[0.010, 0.080],
         nuSf=[0.0, 0.135], chi=[1.0, 0.0], scat=[[0.0, 0.02], [0.0, 0.0]]),
    dict(name="Fuel 2 (inner)",
         D=[1.5, 0.4], Sa=[0.010, 0.085],
         nuSf=[0.0, 0.135], chi=[1.0, 0.0], scat=[[0.0, 0.02], [0.0, 0.0]]),
    dict(name="Fuel 2 + control rod",
         D=[1.5, 0.4], Sa=[0.010, 0.130],
         nuSf=[0.0, 0.135], chi=[1.0, 0.0], scat=[[0.0, 0.02], [0.0, 0.0]]),
    dict(name="Reflector (water)",
         D=[2.0, 0.3], Sa=[0.000, 0.010],
         nuSf=[0.0, 0.0], chi=[0.0, 0.0], scat=[[0.0, 0.04], [0.0, 0.0]]),
    dict(name="Reflector + rod guide",
         D=[2.0, 0.3], Sa=[0.000, 0.055],
         nuSf=[0.0, 0.0], chi=[0.0, 0.0], scat=[[0.0, 0.04], [0.0, 0.0]]),
]


def _stamp(blocks, box, mid):
    """Set material mid on every 10-cm block whose centre is in box."""
    x0, x1, y0, y1 = box
    out = [row[:] for row in blocks]
    nby = len(out)
    for c in range(x0 // 10, x1 // 10):
        for jb in range(y0 // 10, y1 // 10):        # bottom-based row
            out[nby - 1 - jb][c] = mid
    return out


def iaea3d_axial(rod5_depth_cm=80.0, dz=20.0):
    """Axial zones (bottom->top) of IAEA-3D with rod 5 inserted
    `rod5_depth_cm` into the core from the top.  Depth is quantised to
    the layer size dz.  Core span 20..360 cm, total height 380 cm."""
    core_h = 340.0
    depth = max(0.0, min(rod5_depth_cm, core_h))
    n_dep = int(round(depth / dz))                 # rodded layers (top)
    n_core = int(round(core_h / dz))
    base = [row[:] for row in IAEA2D_BLOCKS]       # 4 full rods present
    unrod5 = base
    rod5 = _stamp(base, _ROD5_BOX, 3)
    refl = [[4] * 17 for _ in range(17)]
    reflrod = refl
    for b in _ROD_BOXES_FULL + [_ROD5_BOX]:
        reflrod = _stamp(reflrod, b, 5)
    zones = [dict(label="lower reflector", layers=int(round(20.0 / dz)),
                  blocks=refl)]
    if n_core - n_dep > 0:
        zones.append(dict(label="core (4 full rods)", layers=n_core - n_dep,
                          blocks=unrod5))
    if n_dep > 0:
        zones.append(dict(label="core + rod 5", layers=n_dep, blocks=rod5))
    zones.append(dict(label="upper reflector (+rod guides)",
                      layers=int(round(20.0 / dz)), blocks=reflrod))
    return dict(dz=dz, divz=1, zones=zones)


def preset_iaea3d():
    return dict(
        key="iaea3d",
        title="IAEA-3D PWR benchmark (full x-y-z)",
        description=(
            "THE classic 3-D coarse-mesh benchmark (ANL benchmark book, "
            "problem 11): the IAEA-2D radial lattice extruded to 380 cm "
            "with axial reflectors, four fully-inserted control rods and a "
            "**fifth rod inserted only 80 cm from the top** — physics no "
            "2-D model can represent. Reference k_eff = 1.02903. Try the "
            "**rod insertion sweep** in 🔧 Physics tools for the classic "
            "S-curve."
        ),
        ng=2, pitch=10.0, div=2,
        blocks=[row[:] for row in IAEA2D_BLOCKS],
        materials=[dict(m, D=list(m["D"]), Sa=list(m["Sa"]),
                        nuSf=list(m["nuSf"]), chi=list(m["chi"]),
                        scat=[r[:] for r in m["scat"]])
                   for m in IAEA3D_MATERIALS],
        bc=["reflective", "vacuum", "reflective", "vacuum",
            "vacuum", "vacuum"],                    # W E S N Bottom Top
        gamma=0.4692,
        axial=iaea3d_axial(80.0),
        rod_meta=dict(kind="rod5_insertion", core_h=340.0,
                      description="rod 5 at (30-50, 30-50) cm"),
        ref_keff=1.02903,
        ref_source="ANL-7416 benchmark book, problem 11 (3D IAEA)",
    )


ALL_PRESETS = [preset_smr, preset_iaea3d, preset_iaea2d, preset_c5g7,
               preset_designer, preset_minicore, preset_homogeneous]
