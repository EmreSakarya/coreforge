"""CoreForge — plotly figures.

Color roles follow one system:
  * materials (identity)  -> fixed-order categorical palette (CVD-validated)
  * flux (magnitude)      -> one single-hue sequential ramp per group
  * power (polarity vs 1) -> diverging blue-red around the core average
"""
import numpy as np
import plotly.graph_objects as go

# fixed-order categorical palette (slot order is the CVD-safety mechanism)
CAT = ["#2a78d6", "#1baf7a", "#eda100", "#008300",
       "#4a3aa7", "#e34948", "#e87ba4", "#eb6834"]
# one single-hue ramp per energy group (fast -> thermal -> ...)
GROUP_RAMPS = ["Blues", "Teal", "Purples", "Greens", "Oranges", "Reds",
               "Burg", "Mint"]

_INK = "#0b0b0b"
_MUTED = "#898781"
_GRID = "#e1e0d9"


def _layout(fig, title, equal=True, h=430):
    fig.update_layout(
        title=dict(text=title, font=dict(size=15, color=_INK)),
        template="plotly_white",
        font=dict(family='system-ui, -apple-system, "Segoe UI", sans-serif',
                  color=_INK, size=12),
        margin=dict(l=10, r=10, t=45, b=10),
        height=h,
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="#fcfcfb",
    )
    fig.update_xaxes(title="x [cm]", gridcolor=_GRID, zeroline=False,
                     color=_MUTED, title_font_color=_MUTED)
    fig.update_yaxes(title="y [cm]", gridcolor=_GRID, zeroline=False,
                     color=_MUTED, title_font_color=_MUTED)
    if equal:
        fig.update_yaxes(scaleanchor="x", scaleratio=1)
    return fig


def material_map_fig(mat2d_bottomfirst, materials, x, y, title="Core layout",
                     mini=False):
    """Discrete material-id map with a legend (names, fixed colors).
    mini=True gives a compact thumbnail (no legend / axis titles)."""
    m = np.asarray(mat2d_bottomfirst, dtype=float)
    nmat = len(materials)
    steps = []
    for i in range(nmat):
        c = CAT[i % len(CAT)]
        steps.append([i / nmat, c])
        steps.append([(i + 1) / nmat, c])
    names = np.empty(m.shape, dtype=object)
    for i in range(nmat):
        names[m == i + 1] = materials[i]["name"]
    fig = go.Figure(go.Heatmap(
        z=m, x=x, y=y, colorscale=steps, zmin=0.5, zmax=nmat + 0.5,
        showscale=False, customdata=names,
        hovertemplate="x=%{x:.1f} cm<br>y=%{y:.1f} cm<br>%{customdata}<extra></extra>"))
    if not mini:
        for i, mt in enumerate(materials):
            fig.add_trace(go.Scatter(
                x=[None], y=[None], mode="markers", name=f"{i+1} · {mt['name']}",
                marker=dict(size=12, color=CAT[i % len(CAT)], symbol="square")))
        # legend BELOW the plot so it never collides with the title when
        # many materials wrap the row
        fig.update_layout(legend=dict(orientation="h", yanchor="top",
                                      y=-0.16, x=0, font=dict(size=11)))
    _layout(fig, title if not mini else "", h=210 if mini else 470)
    if not mini:
        fig.update_layout(margin=dict(l=10, r=10, t=45, b=70))
    if mini:
        fig.update_xaxes(title=None, showticklabels=False)
        fig.update_yaxes(title=None, showticklabels=False)
        fig.update_layout(margin=dict(l=4, r=4, t=8, b=4))
    return fig


def field_fig(grid_bottomfirst, x, y, title, ramp="Blues", unit=""):
    fig = go.Figure(go.Heatmap(
        z=grid_bottomfirst, x=x, y=y, colorscale=ramp,
        colorbar=dict(thickness=12, outlinewidth=0, title=unit,
                      tickfont=dict(size=10, color=_MUTED)),
        hovertemplate="x=%{x:.1f} cm<br>y=%{y:.1f} cm<br>%{z:.4e}<extra></extra>"))
    return _layout(fig, title)


def power_fig(power_bottomfirst, fuelmask_bottomfirst, x, y,
              title="Relative power (fuel cells)"):
    p = np.where(fuelmask_bottomfirst, power_bottomfirst, np.nan)
    zmax = float(np.nanmax(p)) if np.isfinite(p).any() else 2.0
    span = max(zmax - 1.0, 1.0 - float(np.nanmin(p)), 0.05)
    fig = go.Figure(go.Heatmap(
        z=p, x=x, y=y, colorscale="RdBu_r", zmid=1.0,
        zmin=1.0 - span, zmax=1.0 + span,
        colorbar=dict(thickness=12, outlinewidth=0, title="P/P̄",
                      tickfont=dict(size=10, color=_MUTED)),
        hovertemplate="x=%{x:.1f} cm<br>y=%{y:.1f} cm<br>P/P̄=%{z:.3f}<extra></extra>"))
    return _layout(fig, title)


def block_power_fig(bp_topfirst, pitch, title="Block powers (P/P̄)",
                    zmid=1.0, colorscale="RdBu_r", fmt=3, unit="P/P̄"):
    """Annotated block map. zmid=1 -> diverging around the average;
    zmid=None -> plain sequential magnitude (e.g. burnup)."""
    bp = np.flipud(np.asarray(bp_topfirst, dtype=float))     # bottom-first
    nby, nbx = bp.shape
    xc = (np.arange(nbx) + 0.5) * pitch
    yc = (np.arange(nby) + 0.5) * pitch
    txt = np.where(np.isfinite(bp), np.round(bp, fmt).astype(str), "")
    kw = {}
    if zmid is not None:
        span = max(np.nanmax(bp) - zmid if np.isfinite(bp).any() else 0.5,
                   zmid - (np.nanmin(bp) if np.isfinite(bp).any() else 0.5),
                   0.05)
        kw = dict(zmid=zmid, zmin=zmid - span, zmax=zmid + span)
    fig = go.Figure(go.Heatmap(
        z=bp, x=xc, y=yc, colorscale=colorscale, **kw,
        text=txt, texttemplate="%{text}", textfont=dict(size=10),
        xgap=2, ygap=2,
        colorbar=dict(thickness=12, outlinewidth=0,
                      tickfont=dict(size=10, color=_MUTED)),
        hovertemplate="x=%{x:.0f} cm<br>y=%{y:.0f} cm<br>" + unit +
                      "=%{z:.3f}<extra></extra>"))
    return _layout(fig, title)


def traverse_fig(x, series, title="Flux traverse", xlabel="x [cm]",
                 ylabel="flux [rel.]"):
    """series: list of (label, values). One categorical slot per series."""
    fig = go.Figure()
    for i, (label, v) in enumerate(series):
        fig.add_trace(go.Scatter(
            x=x, y=v, mode="lines", name=label,
            line=dict(width=2, color=CAT[i % len(CAT)])))
    fig.update_layout(hovermode="x unified",
                      legend=dict(orientation="h", yanchor="bottom", y=1.02))
    _layout(fig, title, equal=False, h=360)
    fig.update_xaxes(title=xlabel)
    fig.update_yaxes(title=ylabel, rangemode="tozero")
    return fig


def residual_fig(residuals, title="Equivalent-fuel residuals"):
    """Bar chart of the inverse-designer residuals [%] per group
    constant (matched constants sit at ~0; sign shown around a zero
    baseline)."""
    keys = [k for k, v in residuals.items()
            if isinstance(v, float) and not k.startswith("kinf")]
    vals = [residuals[k] for k in keys]
    fig = go.Figure(go.Bar(
        x=keys, y=vals, marker_color=CAT[0],
        hovertemplate="%{x}: %{y:+.2f}%<extra></extra>"))
    fig.add_hline(y=0, line_color=_MUTED, line_width=1)
    _layout(fig, title, equal=False, h=320)
    fig.update_xaxes(title="group constant")
    fig.update_yaxes(title="residual (equivalent − target) [%]")
    return fig


def search_fig(history, target, title="Criticality search",
               xlabel="ΔΣa [cm⁻¹]"):
    h = sorted(history)
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=[d for d, _ in h], y=[k for _, k in h],
        mode="lines+markers", name="k_eff",
        line=dict(width=2, color=CAT[0]), marker=dict(size=8)))
    fig.add_hline(y=target, line_dash="dot", line_color=_MUTED,
                  annotation_text=f"target k = {target}")
    _layout(fig, title, equal=False, h=360)
    fig.update_xaxes(title=xlabel)
    fig.update_yaxes(title="k_eff")
    return fig
