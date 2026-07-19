"""CoreForge — UI polish layer (theme CSS + header components).

Kept separate from app.py so the styling can evolve without touching the
logic.  The CSS is deliberately scoped to Streamlit's STABLE data-testid
hooks (stTabs, stMetric, stSidebar, stButton, ...) plus generic element
styling — no dependence on hashed/internal class names — so it degrades
gracefully across Streamlit versions.
"""
import streamlit as st

# palette shared with plots.py (single source of visual truth)
BLUE = "#2a78d6"
BLUE_DK = "#1c5cab"
AQUA = "#1baf7a"
INK = "#0b1220"
MUTED = "#5b6675"
SURF = "#ffffff"
SURF2 = "#eef2f8"
LINE = "#dbe3ef"

_CSS = f"""
<style>
/* ---- typography & base ------------------------------------------- */
html, body, [class*="css"] {{
  font-family: system-ui, -apple-system, "Segoe UI", Roboto, sans-serif;
}}
.block-container {{ padding-top: 1.1rem; max-width: 1500px; }}
h1, h2, h3 {{ color: {INK}; letter-spacing: -0.01em; }}
a {{ color: {BLUE_DK}; }}

/* ---- hero header -------------------------------------------------- */
.cf-hero {{
  background: linear-gradient(120deg, {BLUE} 0%, {BLUE_DK} 55%, #123, 100%);
  background: linear-gradient(120deg, {BLUE} 0%, {BLUE_DK} 60%, #17305c 100%);
  border-radius: 16px; padding: 1.15rem 1.5rem; margin-bottom: 1.1rem;
  color: #fff; box-shadow: 0 6px 24px rgba(28,92,171,0.22);
  display: flex; align-items: center; gap: 1.1rem; flex-wrap: wrap;
}}
.cf-hero .logo {{ font-size: 2.6rem; line-height: 1; filter: drop-shadow(0 2px 4px rgba(0,0,0,.25)); }}
.cf-hero .title {{ font-size: 1.65rem; font-weight: 800; margin: 0; }}
.cf-hero .tag {{ font-size: .93rem; opacity: .92; margin-top: .15rem; }}
.cf-hero .spacer {{ flex: 1; }}
.cf-chip {{
  display: inline-flex; align-items: center; gap: .4rem;
  background: rgba(255,255,255,.16); border: 1px solid rgba(255,255,255,.28);
  padding: .3rem .7rem; border-radius: 999px; font-size: .82rem;
  font-weight: 600; white-space: nowrap;
}}
.cf-chip.ok {{ background: rgba(27,175,122,.28); border-color: rgba(27,175,122,.5); }}
.cf-chip.bad {{ background: rgba(224,59,59,.30); border-color: rgba(224,59,59,.55); }}

/* ---- tabs: pill-style nav ---------------------------------------- */
.stTabs [data-baseweb="tab-list"] {{
  gap: .35rem; background: {SURF2}; padding: .35rem; border-radius: 12px;
  border: 1px solid {LINE}; flex-wrap: wrap;
}}
.stTabs [data-baseweb="tab"] {{
  height: auto; padding: .45rem .85rem; border-radius: 9px;
  background: transparent; color: {MUTED}; font-weight: 600; font-size: .9rem;
  border: 1px solid transparent;
}}
.stTabs [data-baseweb="tab"]:hover {{ background: rgba(42,120,214,.08); color: {BLUE_DK}; }}
.stTabs [aria-selected="true"] {{
  background: {SURF} !important; color: {BLUE_DK} !important;
  border: 1px solid {LINE}; box-shadow: 0 2px 6px rgba(11,18,32,.06);
}}

/* ---- metrics: card look ------------------------------------------ */
[data-testid="stMetric"] {{
  background: {SURF}; border: 1px solid {LINE}; border-radius: 12px;
  padding: .7rem .9rem; box-shadow: 0 1px 3px rgba(11,18,32,.04);
}}
[data-testid="stMetricLabel"] {{ color: {MUTED}; font-weight: 600; }}
[data-testid="stMetricValue"] {{
  color: {INK}; font-weight: 750; font-size: 1.55rem; line-height: 1.15;
}}
/* let long metric values (e.g. '+4,210 pcm', 'supercritical') wrap
   instead of truncating with an ellipsis */
[data-testid="stMetricValue"] > div {{
  white-space: normal; overflow: visible; text-overflow: clip;
}}

/* ---- buttons ------------------------------------------------------ */
.stButton > button, .stDownloadButton > button {{
  border-radius: 10px; font-weight: 650; border: 1px solid {LINE};
  transition: transform .04s ease;
}}
.stButton > button:active {{ transform: translateY(1px); }}
.stButton > button[kind="primary"] {{
  background: {BLUE}; border-color: {BLUE}; box-shadow: 0 2px 8px rgba(42,120,214,.28);
}}

/* ---- containers / expanders / dataframes ------------------------- */
[data-testid="stExpander"] {{ border-radius: 12px; border: 1px solid {LINE}; }}
[data-testid="stVerticalBlockBorderWrapper"] > div {{ border-radius: 12px; }}
hr {{ margin: .9rem 0; border-color: {LINE}; }}

/* ---- sidebar ------------------------------------------------------ */
[data-testid="stSidebar"] {{ background: {SURF2}; border-right: 1px solid {LINE}; }}
[data-testid="stSidebar"] .block-container {{ padding-top: 1.2rem; }}

/* ---- section header helper --------------------------------------- */
.cf-sec {{ margin: .3rem 0 .2rem 0; }}
.cf-sec .h {{ font-size: 1.15rem; font-weight: 750; color: {INK}; }}
.cf-sec .s {{ font-size: .9rem; color: {MUTED}; margin-top: .1rem; }}
.cf-step {{
  display:inline-block; width:1.4rem; height:1.4rem; line-height:1.4rem;
  text-align:center; border-radius:50%; background:{BLUE}; color:#fff;
  font-weight:700; font-size:.82rem; margin-right:.5rem;
}}
</style>
"""


def inject():
    st.markdown(_CSS, unsafe_allow_html=True)


def hero(engine_ok, version, subtitle):
    chip = ('<span class="cf-chip ok">● engine ready</span>' if engine_ok
            else '<span class="cf-chip bad">● engine not built</span>')
    st.markdown(f"""
<div class="cf-hero">
  <div class="logo">⚛️</div>
  <div>
    <div class="title">CoreForge</div>
    <div class="tag">{subtitle}</div>
  </div>
  <div class="spacer"></div>
  {chip}
  <span class="cf-chip">v{version}</span>
  <span class="cf-chip">2-D / 3-D diffusion</span>
</div>
""", unsafe_allow_html=True)


def section(title, subtitle="", step=None):
    s = f'<span class="cf-step">{step}</span>' if step else ""
    sub = f'<div class="s">{subtitle}</div>' if subtitle else ""
    st.markdown(f'<div class="cf-sec"><div class="h">{s}{title}</div>{sub}</div>',
                unsafe_allow_html=True)
