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
def describe_controls(cfg, rod_meta=None):
    """Which live controls make sense for this core."""
    mats = cfg["materials"]
    des = designer_ids(mats)
    fresh = [i for i in des
             if mats[i - 1]["designer"].get("enrich") is not None
             and mats[i - 1]["designer"].get("r_fuel", 0.0) > 0.0]
    return dict(
        swaps=rod_swap_candidates(mats),
        designer=des,
        fresh_fuels=fresh,
        has_rod5=bool(rod_meta and rod_meta.get("kind") == "rod5_insertion"
                      and cfg.get("axial")),
        core_h=float(rod_meta.get("core_h", 340.0)) if rod_meta else None,
        ppm_now=(mats[des[0] - 1]["designer"].get("ppm", 0.0)
                 if des else None),
    )
