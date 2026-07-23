"""CoreForge — live interactive-core mechanisms.

Pure config-in / config-out transformations behind the 🕹 Live core tab,
kept OUT of app.py so verify.py can lock their physics in as permanent
checks (bank-out worth positive, deeper rod -> lower k, more boron ->
lower k, higher enrichment -> higher k, poison monotone):

  * bank_out_cfg    — swap an as-loaded rodded material for its unrodded
                      partner everywhere (2-D and every 3-D zone);
  * rod5_cfg        — rebuild the IAEA-3D axial stack with rod 5 at an
                      arbitrary insertion depth (the S-curve mechanism);
  * boron_cfg       — rebuild every designer material at a new soluble-
                      boron concentration (exact xslib path, the same
                      mechanism as the critical-boron search);
  * enrichment_cfg  — rebuild one FRESH designer fuel at a new enrichment
                      (keeps its boron/geometry);
  * poison_cfg      — uniform delta-Sigma_a on the fuel materials' last
                      group (generic absorber for benchmark-constant
                      cores, where ppm/enrichment are not defined).

`describe_controls` inspects a config and reports which of these controls
are physically meaningful for it, so the UI adapts to ANY core — bundled
benchmark or user-built.
"""
import copy

import runner
import xslib
import presets as _presets


# ----------------------------------------------------------------------
def rod_swap_candidates(materials):
    """Detect as-loaded control-bank materials and their unrodded partners.

    A candidate is a material whose name mentions a rod/CRA and for which
    another material exists with IDENTICAL D and nuSf but a lower
    last-group Sigma_a (the same fuel without the absorber).  Returns a
    list of (rodded_id, unrodded_id, label), ids 1-based; the partner
    with the CLOSEST (largest) Sa below the rodded one is chosen."""
    out = []
    for i, m in enumerate(materials):
        nm = str(m.get("name", "")).lower()
        if not any(t in nm for t in ("rod", "cra")):
            continue
        best = None
        for j, u in enumerate(materials):
            if j == i:
                continue
            if (list(u["nuSf"]) == list(m["nuSf"])
                    and list(u["D"]) == list(m["D"])
                    and u["Sa"][-1] < m["Sa"][-1]):
                if best is None or u["Sa"][-1] > materials[best]["Sa"][-1]:
                    best = j
        if best is not None:
            out.append((i + 1, best + 1, m.get("name", f"mat {i+1}")))
    return out


def bank_out_cfg(cfg, rodded_id, unrodded_id):
    """Withdraw a bank: replace `rodded_id` by `unrodded_id` in the radial
    map and in every 3-D axial zone (the material_swap_worth mechanism)."""
    c = copy.deepcopy(cfg)
    c["blocks"] = [[unrodded_id if v == rodded_id else v for v in row]
                   for row in c["blocks"]]
    if c.get("axial"):
        for z in c["axial"]["zones"]:
            z["blocks"] = [[unrodded_id if v == rodded_id else v
                            for v in row] for row in z["blocks"]]
    return c


# ----------------------------------------------------------------------
def rod5_cfg(cfg, rod_meta, depth_cm):
    """IAEA-3D rod-5 at `depth_cm` from the top of the core: rebuild the
    axial stack with the preset's geometry builder, preserving the
    current dz and divz.  Only valid for rod_meta kind 'rod5_insertion'."""
    if not rod_meta or rod_meta.get("kind") != "rod5_insertion":
        raise ValueError("this core has no parameterised rod-5 "
                         "(rod_meta kind 'rod5_insertion' required)")
    c = copy.deepcopy(cfg)
    dz = float(c["axial"]["dz"])
    divz = int(c["axial"].get("divz", 1))
    c["axial"] = dict(_presets.iaea3d_axial(float(depth_cm), dz), divz=divz)
    return c


def rod_geometry(rod_meta, depth_cm):
    """Absolute rod boxes for the 3-D view: [(x0,x1,y0,y1,depth_cm), ...].
    For IAEA-3D: the four full-length rods plus rod 5 at `depth_cm`."""
    if not rod_meta or rod_meta.get("kind") != "rod5_insertion":
        return []
    core_h = float(rod_meta.get("core_h", 340.0))
    rods = [(x0, x1, y0, y1, core_h)
            for (x0, x1, y0, y1) in _presets._ROD_BOXES_FULL]
    x0, x1, y0, y1 = _presets._ROD5_BOX
    rods.append((x0, x1, y0, y1, max(0.0, min(float(depth_cm), core_h))))
    return rods


# ----------------------------------------------------------------------
def designer_ids(materials):
    """1-based ids of materials carrying designer (xslib) metadata."""
    return [i + 1 for i, m in enumerate(materials) if "designer" in m]


def boron_cfg(cfg, ppm):
    """Rebuild every designer material at soluble boron `ppm` (exact
    pin-cell physics — same mechanism as the critical-boron search).
    Returns (new_cfg, n_rebuilt); non-designer materials are untouched."""
    c = copy.deepcopy(cfg)
    n = 0
    for i, m in enumerate(c["materials"]):
        d = m.get("designer")
        if not d:
            continue
        nm = xslib.cell_from_N(d["N"], float(ppm), d["rho_mod"],
                               d["r_fuel"], d["pitch"],
                               name=m["name"],
                               enrich_label=d.get("enrich"))
        c["materials"][i] = nm
        n += 1
    return c, n


def enrichment_cfg(cfg, mat_id, enrich_wo, ppm=None):
    """Rebuild ONE fresh designer fuel at a new enrichment [w/o U-235],
    keeping its boron (unless `ppm` overrides), moderator density and
    cell geometry.  Only fresh fuels (enrich label present, r_fuel>0)
    are re-designable — depleted/arbitrary-N cells are not."""
    m = cfg["materials"][mat_id - 1]
    d = m.get("designer")
    if not d or d.get("enrich") is None or d.get("r_fuel", 0.0) <= 0.0:
        raise ValueError(f"material {mat_id} is not a fresh designer fuel")
    c = copy.deepcopy(cfg)
    use_ppm = float(d["ppm"] if ppm is None else ppm)
    c["materials"][mat_id - 1] = xslib.pincell_xs(
        float(enrich_wo), use_ppm, d["rho_mod"], d["r_fuel"], d["pitch"])
    return c


# ----------------------------------------------------------------------
def poison_cfg(cfg, dsa):
    """Uniform generic absorber: add `dsa` [1/cm] to the LAST group's
    Sigma_a of every fissile material (boron-like control for cores whose
    constants are plain benchmark numbers).  Guards Sa > 0."""
    c = copy.deepcopy(cfg)
    for i in runner.fissile_ids(c["materials"]):
        sa = c["materials"][i - 1]["Sa"]
        if sa[-1] + dsa <= 0.0:
            raise ValueError("poison would make Sigma_a non-positive")
        sa[-1] = sa[-1] + dsa
    return c


# ----------------------------------------------------------------------
# 2-D -> 3-D lift: give ANY 2-D core a real axial dimension so control
# rods can be moved CONTINUOUSLY on it.
# ----------------------------------------------------------------------
def lift_to_3d(cfg, axial3d=None, mode="physical"):
    """Turn a 2-D core into an explicit 3-D core.

    mode='physical' (requires `axial3d` metadata: b2, core_h, dz,
    refl_mat, refl_cm): UNFOLD the axial buckling from every material
    (Sa -= D*b2 — the fold was the fundamental-mode axial-leakage model),
    stack refl/core/refl zones and open vacuum bottom/top.  This
    reconstructs the honest 3-D counterpart of the folded 2-D model
    (IAEA-2D lifted this way is essentially IAEA-3D's geometry).

    mode='buckled': keep the folded Sa, reflective bottom/top — exactly
    equivalent to the 2-D solve (verify.py's 2D==extruded-3D identity),
    endpoint-exact but with a flat axial flux.
    """
    if cfg.get("axial"):
        raise ValueError("core is already 3-D")
    c = copy.deepcopy(cfg)
    if mode == "physical":
        if not axial3d:
            raise ValueError("physical lift needs axial3d metadata "
                             "(b2, core_h, dz, refl_mat, refl_cm)")
        b2 = float(axial3d.get("b2", 0.0))
        core_h = float(axial3d["core_h"])
        dz = float(axial3d.get("dz", 20.0))
        rmat = int(axial3d["refl_mat"])
        rcm = float(axial3d.get("refl_cm", 20.0))
        if b2 > 0.0:
            for m in c["materials"]:
                for g in range(c["ng"]):
                    m["Sa"][g] = m["Sa"][g] - m["D"][g] * b2
                    if m["Sa"][g] <= 0.0:
                        m["Sa"][g] = 1e-6
        nby = len(c["blocks"]); nbx = len(c["blocks"][0])
        refl = [[rmat] * nbx for _ in range(nby)]
        zones = []
        n_r = int(round(rcm / dz))
        if n_r > 0:
            zones.append(dict(label="lower reflector", layers=n_r,
                              blocks=refl))
        zones.append(dict(label="core", layers=int(round(core_h / dz)),
                          blocks=[row[:] for row in c["blocks"]]))
        if n_r > 0:
            zones.append(dict(label="upper reflector", layers=n_r,
                              blocks=[row[:] for row in refl]))
        c["axial"] = dict(dz=dz, divz=1, zones=zones)
        while len(c["bc"]) < 6:
            c["bc"].append("vacuum")
        c["bc"][4] = "vacuum"; c["bc"][5] = "vacuum"
    else:                                   # buckled, endpoint-exact
        dz = float((axial3d or {}).get("dz", 20.0))
        core_h = float((axial3d or {}).get("core_h", 200.0))
        c["axial"] = dict(dz=dz, divz=1, zones=[
            dict(label="core", layers=int(round(core_h / dz)),
                 blocks=[row[:] for row in c["blocks"]])])
        while len(c["bc"]) < 6:
            c["bc"].append("reflective")
        c["bc"][4] = "reflective"; c["bc"][5] = "reflective"
    return c


def _flatten_layers(axial):
    """Zones -> flat per-layer list of (label, blocks_ref). Bottom-first."""
    out = []
    for z in axial["zones"]:
        for _ in range(int(z["layers"])):
            out.append((z.get("label", "zone"), z["blocks"]))
    return out


def _regroup_layers(layers, dz, divz):
    """Flat per-layer list -> zones (merging identical consecutive maps)."""
    zones = []
    for label, blocks in layers:
        if zones and zones[-1]["blocks"] is blocks:
            zones[-1]["layers"] += 1
        else:
            zones.append(dict(label=label, layers=1, blocks=blocks))
    return dict(dz=dz, divz=divz, zones=zones)


def fissile_span(cfg):
    """(z_bot, z_top) of the fissile region [cm] of a 3-D core."""
    fids = set(runner.fissile_ids(cfg["materials"]))
    dz = float(cfg["axial"]["dz"])
    layers = _flatten_layers(cfg["axial"])
    idx = [i for i, (_, b) in enumerate(layers)
           if any(v in fids for row in b for v in row)]
    if not idx:
        raise ValueError("no fissile layers in this core")
    return min(idx) * dz, (max(idx) + 1) * dz


def bank_depth_cfg(cfg, rodded_id, unrodded_id, depth_cm):
    """CONTINUOUS bank insertion for ANY 3-D core: within the fissile
    span, layers in the top `depth_cm` keep the rodded material; below
    it every `rodded_id` becomes `unrodded_id`.  depth is quantised to
    the layer thickness dz.  Generalises the IAEA-3D rod-5 mechanism to
    every core that has a rodded/unrodded material pair."""
    if not cfg.get("axial"):
        raise ValueError("bank_depth_cfg needs a 3-D core "
                         "(lift a 2-D core with lift_to_3d first)")
    c = copy.deepcopy(cfg)
    dz = float(c["axial"]["dz"])
    divz = int(c["axial"].get("divz", 1))
    z_bot, z_top = fissile_span(c)
    depth = max(0.0, min(float(depth_cm), z_top - z_bot))
    swapped = {}                            # id(blocks) -> swapped copy
    layers = []
    for i, (label, blocks) in enumerate(_flatten_layers(c["axial"])):
        z_mid = (i + 0.5) * dz
        in_core = z_bot < z_mid < z_top
        rodded = in_core and z_mid > z_top - depth
        if in_core and not rodded:
            key = id(blocks)
            if key not in swapped:
                swapped[key] = [[unrodded_id if v == rodded_id else v
                                 for v in row] for row in blocks]
            layers.append((label + " (bank out)", swapped[key]))
        else:
            layers.append((label, blocks))
    c["axial"] = _regroup_layers(layers, dz, divz)
    return c


def rod_boxes_for(cfg, rodded_id, depth_cm, core_h):
    """Rod columns for the 3-D view of a generic core: one box per block
    holding `rodded_id` in the (unswapped) radial map, descending
    `depth_cm` from the top of the fissile span."""
    pitch = float(cfg["pitch"])
    if cfg.get("axial"):
        fids = set(runner.fissile_ids(cfg["materials"]))
        base = None                          # densest fissile map
        bestn = -1
        for z in cfg["axial"]["zones"]:
            n = sum(1 for row in z["blocks"] for v in row if v in fids)
            if n > bestn:
                base, bestn = z["blocks"], n
    else:
        base = cfg["blocks"]
    nby = len(base)
    out = []
    for r, row in enumerate(base):           # rows are top-first
        for cix, v in enumerate(row):
            if v == rodded_id:
                yb = (nby - 1 - r) * pitch   # bottom-based y
                out.append((cix * pitch, (cix + 1) * pitch,
                            yb, yb + pitch,
                            max(0.0, min(float(depth_cm), core_h))))
    return out


# ----------------------------------------------------------------------
def describe_controls(cfg, rod_meta=None, axial3d=None):
    """Which live controls make sense for this core.

    Bank pairs split by physics: a FISSILE rodded material (fuel+rod,
    CRA) supports CONTINUOUS insertion depth (`depth_banks`); a
    non-fissile one (reflector rod guide) is a binary in/out swap
    (`binary_banks`) — partial depth has no meaning outside the fuel."""
    mats = cfg["materials"]
    des = designer_ids(mats)
    fresh = [i for i in des
             if mats[i - 1]["designer"].get("enrich") is not None
             and mats[i - 1]["designer"].get("r_fuel", 0.0) > 0.0]
    swaps = rod_swap_candidates(mats)
    depth_banks = [s for s in swaps
                   if any(v > 0 for v in mats[s[0] - 1]["nuSf"])]
    binary_banks = [s for s in swaps if s not in depth_banks]
    liftable = (not cfg.get("axial")) and bool(depth_banks)
    return dict(
        swaps=swaps,
        depth_banks=depth_banks,
        binary_banks=binary_banks,
        designer=des,
        fresh_fuels=fresh,
        has_rod5=bool(rod_meta and rod_meta.get("kind") == "rod5_insertion"
                      and cfg.get("axial")),
        core_h=float(rod_meta.get("core_h", 340.0)) if rod_meta else None,
        liftable=liftable,
        lift_mode=("physical" if (liftable and axial3d) else
                   "buckled" if liftable else None),
        ppm_now=(mats[des[0] - 1]["designer"].get("ppm", 0.0)
                 if des else None),
    )
