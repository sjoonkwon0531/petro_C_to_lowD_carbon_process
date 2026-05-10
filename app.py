"""
Petro-Carbon to Low-Dimensional Carbon Semiconductor Digital Twin
==================================================================
Streamlit application for process simulation, LCA, LCOE, and TEA.
Covers: Graphene (CVD), CNT, GNR, SiC, SEG (Epitaxial Graphene on SiC)

To run locally:
    pip install streamlit numpy pandas plotly scipy
    streamlit run app.py

To deploy:
    Push to GitHub, connect to Streamlit Community Cloud.
"""

import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from dataclasses import dataclass, field
from typing import Dict, List, Tuple
import math

# ============================================================
# CONSTANTS & DATA
# ============================================================

kB_eV = 8.617e-5  # Boltzmann constant eV/K

# --- Feedstock Database ---
FEEDSTOCKS = {
    "Methane (CH4)": {"C_yield": 0.75, "purity": 0.999, "cost_aramco": 0.18, "cost_market": 0.50, "co2_kg_per_kg": 2.75, "Ea_decomp_eV": 2.74, "MW": 16.04},
    "Ethylene (C2H4)": {"C_yield": 0.857, "purity": 0.995, "cost_aramco": 0.45, "cost_market": 1.00, "co2_kg_per_kg": 1.73, "Ea_decomp_eV": 2.50, "MW": 28.05},
    "Propane (C3H8)": {"C_yield": 0.818, "purity": 0.995, "cost_aramco": 0.22, "cost_market": 0.55, "co2_kg_per_kg": 2.99, "Ea_decomp_eV": 2.40, "MW": 44.10},
    "Toluene (C7H8)": {"C_yield": 0.913, "purity": 0.99, "cost_aramco": 0.40, "cost_market": 0.90, "co2_kg_per_kg": 3.35, "Ea_decomp_eV": 2.10, "MW": 92.14},
    "Pet-Coke": {"C_yield": 0.90, "purity": 0.95, "cost_aramco": 0.05, "cost_market": 0.15, "co2_kg_per_kg": 3.50, "Ea_decomp_eV": None, "MW": None},
}

# --- CVD Methods ---
CVD_METHODS = {
    "Thermal CVD (Hot-Wall)": {"T_range": (900, 1100), "P_range": (0.1, 760), "E_fixed_kWh": 8.0, "E_marginal_kWh_m2": 5.0, "quality_base": 0.95, "capex_M": 8, "trl": "7-8"},
    "PECVD": {"T_range": (300, 700), "P_range": (0.01, 20), "E_fixed_kWh": 5.0, "E_marginal_kWh_m2": 15.0, "quality_base": 0.80, "capex_M": 12, "trl": "5-6"},
    "Solar-Thermal CVD": {"T_range": (800, 1200), "P_range": (0.1, 760), "E_fixed_kWh": 0.5, "E_marginal_kWh_m2": 1.0, "quality_base": 0.90, "capex_M": 6, "trl": "3-4"},
    "Joule Heating R2R": {"T_range": (900, 1050), "P_range": (0.1, 5), "E_fixed_kWh": 2.0, "E_marginal_kWh_m2": 8.0, "quality_base": 0.85, "capex_M": 20, "trl": "5-7"},
    "EM Induction Heating": {"T_range": (1000, 1400), "P_range": (0.1, 760), "E_fixed_kWh": 4.0, "E_marginal_kWh_m2": 10.0, "quality_base": 0.97, "capex_M": 10, "trl": "4-5"},
}

SUBSTRATES = {
    "Cu Polycrystalline": {"cost_m2": 25, "reusable": 3, "mu_sub": 25000, "transfer": True},
    "Cu(111) Single Crystal": {"cost_m2": 500, "reusable": 10, "mu_sub": 45000, "transfer": True},
    "c-Sapphire": {"cost_m2": 150, "reusable": 50, "mu_sub": 15000, "transfer": False},
    "SiC (Epitaxial)": {"cost_m2": 800, "reusable": 100, "mu_sub": 30000, "transfer": False},
}

# --- SiC/SEG Process from PPT ---
SIC_STAGES = {
    "Stage 1: Acheson (Pet-Coke to SiC Powder)": {
        "T_C": 2200, "input_petcoke_ton": 12, "input_silica_ton": 18,
        "output_SiC_ton": 2.4, "yield_pct": 20,
        "capex_M": 16, "opex_month_K": 292, "energy_kWh_month": 300000,
        "trl": "9", "co2_kg_per_ton_SiC": 8500
    },
    "Stage 2: PVT (SiC Powder to Boule)": {
        "T_C": 2050, "input_SiC_ton": 2.4, "output_boules": 304,
        "boule_mass_kg": 7.0, "growth_rate_mm_hr": 0.4,
        "capex_M": 176, "opex_month_K": 830, "energy_kWh_month": 400000,
        "trl": "8", "reactors": 47, "reactor_cost_M": 3.0
    },
    "Stage 3: Wafering (Boule to Wafer)": {
        "input_boules": 304, "output_wafers": 10000,
        "wafers_per_boule": 33, "kerf_loss_pct": 40,
        "capex_M": 33, "opex_month_K": 826,
        "trl": "9"
    },
    "Stage 4: SEG (SiC Wafer to Graphene)": {
        "T_C": 1400, "atmosphere": "Ar", "time_hr": 3,
        "output_wafers": 10000, "wafers_per_batch": 10,
        "batches_per_month": 1000, "furnaces": 5,
        "capex_M": 26, "opex_month_K": 467,
        "trl": "4-5", "bandgap_eV": 0.6, "mobility_cm2Vs": 5000
    }
}

# --- Electricity Costs by Region ---
ELECTRICITY = {
    "Saudi Arabia": 0.048,
    "USA": 0.10,
    "EU": 0.12,
    "China": 0.08,
    "Japan/Korea": 0.14,
}

# --- Grid Emission Factors (kg CO2/kWh) ---
GRID_EF = {
    "Saudi Arabia": 0.40,
    "USA": 0.39,
    "EU": 0.23,
    "China": 0.55,
    "Japan/Korea": 0.45,
    "Solar-Thermal": 0.02,
}


# ============================================================
# PHYSICS ENGINE (v2.1 from previous work)
# ============================================================

def simulate_graphene_cvd(T_C, P_torr, ch4_h2, t_min, o2_ppm, method_key, feed_key, sub_key, wafer_diam_cm):
    """Simulate graphene CVD growth with v2.1 corrected physics."""
    cvd = CVD_METHODS[method_key]
    feed = FEEDSTOCKS[feed_key]
    sub = SUBSTRATES[sub_key]
    T_K = T_C + 273.15

    # Nucleation (dual-regime)
    logP = np.log10(max(0.1, P_torr))
    E_nuc_desorp = 4.0 + 5.0 * min(1, max(0, (logP + 1) / 3.88))
    T_cross_K = 1050 + 273.15
    A_nuc = 1e18
    if T_K < T_cross_K:
        nuc_density = A_nuc * np.exp(-2.7 / (kB_eV * T_K)) * ch4_h2**1.0 * (1 + 0.05 * o2_ppm)**2 * np.exp(-(T_cross_K - T_K) * 0.002)
    else:
        nuc_density = A_nuc * np.exp(-E_nuc_desorp / (kB_eV * T_K)) * ch4_h2**1.2 * (1 + 0.05 * o2_ppm)**2

    # Surface diffusion
    D_surf = 1e-4 * np.exp(-0.7 * 96.485 / (8.314e-3 * T_K))

    # Growth rate: calibrated empirical model WITH feedstock dependence
    # Ea_growth depends on feedstock C-H bond dissociation energy
    # Methane: Ea=2.74 eV (hardest to crack), Toluene: Ea=2.10 eV (easier)
    Ea_growth = feed.get("Ea_decomp_eV", 2.5)  # feedstock-specific activation energy
    # Plasma enhancement: plasma methods lower the effective Ea by providing
    # radical species that bypass thermal decomposition
    if "PECVD" in method_key or "Plasma" in method_key:
        Ea_growth = Ea_growth * 0.45  # plasma reduces Ea by ~55% (radical-driven)
    C_yield_factor = feed["C_yield"]  # higher C yield = more carbon per molecule
    purity_factor = feed["purity"]  # impurities reduce effective flux
    G0 = 5e8  # um/s, pre-exponential (calibrated for methane at 1050C)
    G_attach_um_s = G0 * np.exp(-Ea_growth / (kB_eV * T_K)) * np.sqrt(ch4_h2) * P_torr**0.3 * C_yield_factor * purity_factor
    P_H2 = P_torr * (1 - ch4_h2) / (1 + ch4_h2)
    # H2 etching: independent of feedstock (only depends on H2 and temperature)
    G_etch_um_s = 1e7 * np.exp(-4.0 / (kB_eV * T_K)) * P_H2**0.5
    G_net_um_s = max(0, G_attach_um_s - G_etch_um_s)
    G_radial = G_net_um_s * 1e-4  # convert um/s to cm/s

    # CMC check
    CMC = 1e-6 * max(0.01, P_H2)**1.5 * np.exp(-3.4 * 96.485 / (8.314e-3 * T_K))
    above_CMC = (ch4_h2 * P_torr) > CMC

    # JMAK coverage
    t_sec = t_min * 60
    coverage_raw = 1 - np.exp(-np.pi * nuc_density * G_radial**2 * t_sec**2)
    coverage = min(1.0, max(0, coverage_raw)) if above_CMC else max(0, coverage_raw * 0.3)

    # Domain size
    domain_um = min(1000, G_radial * t_sec * 1e4)
    domain_nuc = 1e4 / np.sqrt(max(1, nuc_density))
    eff_domain = min(domain_um, domain_nuc)

    # Layer count (Lewis phase diagram)
    if P_torr > 100:
        layers = 1.0 if ch4_h2 < 0.005 else 1 + (ch4_h2 - 0.005) * 30
    else:
        layers = 1.0 if ch4_h2 < 0.02 else 1 + (ch4_h2 - 0.02) * 12.5 if ch4_h2 < 0.1 else 2 + (ch4_h2 - 0.1) * 10

    # I_D/I_G (structured model with feedstock purity effect)
    gb = 0.15 / max(0.01, eff_domain)
    o2_t = o2_ppm * 0.002
    anneal = 0.5 * np.exp(-3.5 * 96.485 / (8.314e-3 * T_K))
    G_opt = G_attach_um_s * 0.7
    rate_t = 0.1 * (G_net_um_s / max(1e-10, G_opt) - 1)**2 if G_net_um_s > 0 else 0.3
    method_t = (1 - cvd["quality_base"]) * 0.8
    # Feedstock purity: impurities (1 - purity) introduce point defects
    purity_defect = (1 - feed["purity"]) * 5.0  # 99.9% purity -> 0.005; 95% -> 0.25
    ID_IG = max(0.005, min(3.0, (gb + o2_t + anneal + rate_t + method_t + purity_defect) / sub.get("quality_mult", 1.0)))

    # Cancado defect density
    lam = 514.5
    LD_sq = 1.8e-9 * lam**4 / max(0.001, ID_IG)
    LD_nm = np.sqrt(LD_sq)
    nD = 1e14 / (np.pi * LD_sq) if LD_sq > 0 else 0

    # Mobility (Matthiessen)
    mu = sub["mu_sub"] / (1 + nD / 1e12)

    # Sheet resistance (percolation)
    theta_c = 0.67
    if coverage >= theta_c:
        sigma = 1.6e-19 * 3e12 * mu * ((coverage - theta_c) / (1 - theta_c))**1.3
        Rs = 1 / sigma if sigma > 0 else 1e8
    else:
        Rs = 1e8

    # Grade
    if ID_IG < 0.10 and coverage > 0.99:
        grade = "Pristine"
    elif ID_IG < 0.30 and coverage > 0.90:
        grade = "Controlled"
    else:
        grade = "Low-Grade"

    wafer_area_cm2 = np.pi * (wafer_diam_cm / 2)**2

    return {
        "coverage": coverage, "ID_IG": ID_IG, "LD_nm": LD_nm, "nD": nD,
        "mu": mu, "Rs": Rs, "domain_um": eff_domain, "layers": layers,
        "grade": grade, "above_CMC": above_CMC, "nuc_density": nuc_density,
        "G_radial": G_radial, "wafer_area_cm2": wafer_area_cm2,
    }


def compute_economics(params, growth, method_key, feed_key, sub_key, region, use_aramco):
    """Compute per-wafer economics."""
    cvd = CVD_METHODS[method_key]
    feed = FEEDSTOCKS[feed_key]
    sub = SUBSTRATES[sub_key]
    elec = ELECTRICITY[region]
    A_m2 = growth["wafer_area_cm2"] * 1e-4
    cycle_hr = params["t_min"] / 60 + 0.5

    E_kWh = cvd["E_fixed_kWh"] + cvd["E_marginal_kWh_m2"] * A_m2
    energy_cost = E_kWh * elec
    C_mass = A_m2 * 7.6e-7 * growth["layers"]
    fc = feed["cost_aramco"] if use_aramco else feed["cost_market"]
    feed_cost = (C_mass / feed["C_yield"]) * fc * 10
    sub_cost = sub["cost_m2"] * A_m2 / sub["reusable"]
    gas_cost = 0.5 * A_m2 + 0.2
    transfer_cost = (15 * A_m2 + 2) if sub["transfer"] else 0
    labor_hr = 35 if use_aramco else 65
    labor_cost = labor_hr * cycle_hr / 8
    total = energy_cost + feed_cost + sub_cost + gas_cost + transfer_cost + labor_cost

    # CO2
    ef = GRID_EF.get(region, 0.40)
    if "Solar" in method_key:
        ef = GRID_EF["Solar-Thermal"]
    co2_energy = E_kWh * ef
    co2_feed = (C_mass / feed["C_yield"]) * feed["co2_kg_per_kg"]
    co2_total = co2_energy + co2_feed

    # Revenue
    P0 = {"Pristine": 80, "Controlled": 25, "Low-Grade": 5}[growth["grade"]]
    floor = {"Pristine": 15, "Controlled": 3, "Low-Grade": 0.5}[growth["grade"]]
    annual_wafers = (8760 * 0.80) / cycle_hr
    ann_area = annual_wafers * growth["wafer_area_cm2"]
    price_cm2 = max(floor, P0 * max(1, ann_area / 100)**(-0.15))
    revenue = price_cm2 * growth["wafer_area_cm2"]
    margin = (revenue - total) / revenue * 100 if revenue > 0 else -100

    return {
        "energy_cost": energy_cost, "feed_cost": feed_cost, "sub_cost": sub_cost,
        "gas_cost": gas_cost, "transfer_cost": transfer_cost, "labor_cost": labor_cost,
        "total_cost": total, "co2_energy": co2_energy, "co2_feed": co2_feed,
        "co2_total": co2_total, "E_kWh": E_kWh, "price_cm2": price_cm2,
        "revenue": revenue, "margin": margin, "annual_wafers": annual_wafers,
    }


def compute_sic_seg_economics(wafers_month, region, use_aramco):
    """Compute SiC/SEG process chain economics from PPT data."""
    elec = ELECTRICITY[region]
    ef = GRID_EF.get(region, 0.40)
    scale = wafers_month / 10000  # PPT baseline is 10,000 wafers/month

    s1 = SIC_STAGES["Stage 1: Acheson (Pet-Coke to SiC Powder)"]
    s2 = SIC_STAGES["Stage 2: PVT (SiC Powder to Boule)"]
    s3 = SIC_STAGES["Stage 3: Wafering (Boule to Wafer)"]
    s4 = SIC_STAGES["Stage 4: SEG (SiC Wafer to Graphene)"]

    petcoke_cost = s1["input_petcoke_ton"] * scale * (50 if use_aramco else 150)
    silica_cost = s1["input_silica_ton"] * scale * 30
    energy_s1 = s1["energy_kWh_month"] * scale * elec
    energy_s2 = s2["energy_kWh_month"] * scale * elec

    capex_total = s1["capex_M"] + s2["capex_M"] + s3["capex_M"] + s4["capex_M"] + 5 + 75  # AI + Infra
    opex_month = (s1["opex_month_K"] + s2["opex_month_K"] + s3["opex_month_K"] + s4["opex_month_K"]) * scale / 1000  # $M

    depreciation_month = capex_total / (10 * 12)  # 10yr straight-line

    # Revenue
    bare_sic_price = 800
    seg_price = 1200
    seg_fraction = 0.7  # 70% SEG, 30% bare SiC
    revenue_month = wafers_month * (seg_fraction * seg_price + (1 - seg_fraction) * bare_sic_price) / 1e6

    ebitda_month = revenue_month - opex_month
    net_income_month = ebitda_month - depreciation_month

    # CO2
    co2_acheson = s1["co2_kg_per_ton_SiC"] * s1["output_SiC_ton"] * scale  # kg/month
    co2_energy = (s1["energy_kWh_month"] + s2["energy_kWh_month"]) * scale * ef
    co2_per_wafer = (co2_acheson + co2_energy) / wafers_month

    cost_per_wafer = opex_month * 1e6 / wafers_month + depreciation_month * 1e6 / wafers_month

    return {
        "capex_M": capex_total, "opex_month_M": opex_month,
        "revenue_month_M": revenue_month, "ebitda_month_M": ebitda_month,
        "net_income_month_M": net_income_month,
        "cost_per_wafer": cost_per_wafer,
        "co2_per_wafer_kg": co2_per_wafer,
        "payback_years": capex_total / max(0.01, ebitda_month * 12),
        "petcoke_cost_month": petcoke_cost,
        "annual_revenue_M": revenue_month * 12,
        "annual_ebitda_M": ebitda_month * 12,
    }


def compute_lca(method_key, feed_key, region, E_kWh, co2_feed, co2_energy, wafer_area_cm2):
    """Compute LCA indicators per cm2."""
    ef = GRID_EF.get(region, 0.40)
    if "Solar" in method_key:
        ef = GRID_EF["Solar-Thermal"]

    gwp = (co2_energy + co2_feed) / max(1, wafer_area_cm2) * 1000  # g CO2-eq / cm2
    ced = E_kWh / max(1, wafer_area_cm2) * 1000  # Wh/cm2
    water = E_kWh * 2 / max(1, wafer_area_cm2)  # mL/cm2 (rough estimate)
    acid = gwp * 0.003  # rough SO2-eq proxy

    return {"gwp_g_cm2": gwp, "ced_Wh_cm2": ced, "water_mL_cm2": water, "acid_ug_cm2": acid * 1000}


def npv_irr_analysis(capex_M, annual_revenue_M, annual_opex_M, years=10, wacc=0.10):
    """Simple NPV and IRR calculation."""
    annual_cf = annual_revenue_M - annual_opex_M
    npv = -capex_M + sum(annual_cf / (1 + wacc)**t for t in range(1, years + 1))
    # IRR by bisection
    irr = None
    for r in np.linspace(-0.5, 3.0, 1000):
        test_npv = -capex_M + sum(annual_cf / (1 + r)**t for t in range(1, years + 1))
        if test_npv <= 0:
            irr = r
            break
    payback = capex_M / max(0.01, annual_cf) if annual_cf > 0 else float('inf')
    return {"npv_M": npv, "irr_pct": (irr or 0) * 100, "payback_yr": payback}


# ============================================================
# STREAMLIT APP
# ============================================================

st.set_page_config(page_title="Petro-C Digital Twin", page_icon="C", layout="wide")

st.markdown("""
<style>
    .block-container {padding: 1rem 2rem;}
    .metric-card {background: #f0f8ff; border-radius: 8px; padding: 12px; text-align: center; border: 1px solid #dde;}
    .stTabs [data-baseweb="tab-list"] {gap: 4px;}
</style>
""", unsafe_allow_html=True)

st.title("Petro-Carbon to Carbon Semiconductor: Digital Twin")
st.caption("Process Simulation | LCA | LCOE | Techno-Economic Analysis | Aramco Advanced Materials R&D")

tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "Graphene CVD Simulator",
    "SiC / SEG Process Chain",
    "LCA & Environmental",
    "Techno-Economic Analysis",
    "Multi-Product Comparison"
])

# ============================================
# TAB 1: GRAPHENE CVD
# ============================================
with tab1:
    col_ctrl, col_res = st.columns([1, 2])

    with col_ctrl:
        st.subheader("Process Configuration")
        method = st.selectbox("Energy Source", list(CVD_METHODS.keys()))
        feedstock = st.selectbox("Carbon Feedstock", [k for k in FEEDSTOCKS if k != "Pet-Coke"])
        substrate = st.selectbox("Substrate", list(SUBSTRATES.keys()))
        region = st.selectbox("Production Region", list(ELECTRICITY.keys()))
        use_aramco = st.checkbox("Use Aramco internal pricing", value=True)

        st.subheader("Growth Parameters")
        cvd = CVD_METHODS[method]
        T_lo, T_hi = cvd["T_range"]
        T = st.slider("Temperature (C)", int(T_lo), int(T_hi), int((T_lo + T_hi) / 2), 10)
        P_lo, P_hi = float(cvd["P_range"][0]), float(cvd["P_range"][1])
        P_lo = max(0.1, P_lo)  # ensure min >= step size (0.1) to avoid Streamlit errors
        P_max_display = min(100.0, P_hi)
        P_default = min(0.5, P_max_display)  # ensure default is within range
        P_default = max(P_lo, P_default)     # ensure default >= min
        P = st.slider("Pressure (Torr)", P_lo, P_max_display, P_default, 0.1)
        ch4_h2 = st.slider("CH4/H2 Ratio", 0.001, 0.200, 0.015, 0.001, format="%.3f")
        t_min = st.slider("Growth Time (min)", 1, 180, 30)
        o2_ppm = st.slider("Residual O2 (ppm)", 0.0, 50.0, 2.0, 0.5)
        wafer_d = st.slider("Wafer Diameter (cm)", 2.5, 30.0, 15.0, 0.5)

    params = {"T_C": T, "P_torr": P, "ch4_h2": ch4_h2, "t_min": t_min, "o2_ppm": o2_ppm, "wafer_diam_cm": wafer_d}
    g = simulate_graphene_cvd(T, P, ch4_h2, t_min, o2_ppm, method, feedstock, substrate, wafer_d)
    e = compute_economics(params, g, method, feedstock, substrate, region, use_aramco)
    lca = compute_lca(method, feedstock, region, e["E_kWh"], e["co2_feed"], e["co2_energy"], g["wafer_area_cm2"])

    with col_res:
        # Grade badge
        grade_colors = {"Pristine": "green", "Controlled": "blue", "Low-Grade": "orange"}
        cmc_warn = " | BELOW CMC" if not g["above_CMC"] else ""
        st.markdown(f"### :{grade_colors[g['grade']]}[{g['grade']} Grade]{cmc_warn}")

        # Key metrics
        m1, m2, m3, m4, m5, m6 = st.columns(6)
        m1.metric("I_D/I_G", f"{g['ID_IG']:.3f}")
        m2.metric("Coverage", f"{g['coverage']*100:.1f}%")
        m3.metric("Mobility", f"{g['mu']:.0f} cm2/Vs")
        m4.metric("Domain", f"{g['domain_um']:.1f} um")
        m5.metric("Cost/Wafer", f"${e['total_cost']:.2f}")
        m6.metric("Margin", f"{e['margin']:.1f}%")

        # Temp sweep chart
        st.subheader("Temperature Sweep")
        temps = np.linspace(cvd["T_range"][0], cvd["T_range"][1], 30)
        sweep = [simulate_graphene_cvd(t, P, ch4_h2, t_min, o2_ppm, method, feedstock, substrate, wafer_d) for t in temps]
        fig_sweep = make_subplots(specs=[[{"secondary_y": True}]])
        fig_sweep.add_trace(go.Scatter(x=temps, y=[s["ID_IG"] for s in sweep], name="I_D/I_G", line=dict(color="#e74c3c")), secondary_y=False)
        fig_sweep.add_trace(go.Scatter(x=temps, y=[s["coverage"]*100 for s in sweep], name="Coverage %", line=dict(color="#2ecc71"), fill="tozeroy"), secondary_y=True)
        fig_sweep.update_layout(height=300, margin=dict(t=30, b=30), xaxis_title="Temperature (C)")
        fig_sweep.update_yaxes(title_text="I_D/I_G", secondary_y=False)
        fig_sweep.update_yaxes(title_text="Coverage %", secondary_y=True)
        st.plotly_chart(fig_sweep, use_container_width=True)

        # Cost breakdown
        cc1, cc2 = st.columns(2)
        with cc1:
            st.subheader("Cost Breakdown")
            cost_df = pd.DataFrame({
                "Component": ["Energy", "Feedstock", "Substrate", "Gas", "Transfer", "Labor"],
                "Cost ($)": [e["energy_cost"], e["feed_cost"], e["sub_cost"], e["gas_cost"], e["transfer_cost"], e["labor_cost"]]
            })
            fig_cost = px.bar(cost_df, x="Component", y="Cost ($)", color="Component", height=280)
            fig_cost.update_layout(showlegend=False, margin=dict(t=10, b=30))
            st.plotly_chart(fig_cost, use_container_width=True)
        with cc2:
            st.subheader("LCA Indicators (per cm2)")
            lca_df = pd.DataFrame({
                "Indicator": ["GWP (g CO2-eq)", "Energy (Wh)", "Water (mL)", "Acidification (ug SO2-eq)"],
                "Value": [lca["gwp_g_cm2"], lca["ced_Wh_cm2"], lca["water_mL_cm2"], lca["acid_ug_cm2"]]
            })
            st.dataframe(lca_df, use_container_width=True, hide_index=True)
            st.metric("CO2 per wafer", f"{e['co2_total']*1000:.1f} g")

# ============================================
# TAB 2: SiC / SEG
# ============================================
with tab2:
    st.subheader("Pet-Coke to SiC to SEG Process Chain")
    st.caption("Based on PPT data: 4-stage process, baseline 10,000 wafers/month")

    sc1, sc2 = st.columns([1, 2])
    with sc1:
        wafers_m = st.slider("Monthly Wafer Output", 1000, 30000, 10000, 1000)
        seg_region = st.selectbox("Region", list(ELECTRICITY.keys()), key="seg_region")
        seg_aramco = st.checkbox("Aramco pricing (pet-coke $50/ton)", value=True, key="seg_aramco")
        seg_price = st.slider("SEG Wafer Price ($/wafer)", 500, 3000, 1200, 50)
        sic_price = st.slider("Bare SiC Price ($/wafer)", 300, 1500, 800, 50)
        seg_frac = st.slider("SEG fraction of output", 0.0, 1.0, 0.7, 0.05)

    sic = compute_sic_seg_economics(wafers_m, seg_region, seg_aramco)
    # Override prices
    rev_override = wafers_m * (seg_frac * seg_price + (1 - seg_frac) * sic_price) / 1e6
    ebitda_override = rev_override - sic["opex_month_M"]

    with sc2:
        st.markdown("### Financial Summary (Monthly)")
        f1, f2, f3, f4 = st.columns(4)
        f1.metric("CAPEX Total", f"${sic['capex_M']:.0f}M")
        f2.metric("Monthly Revenue", f"${rev_override:.1f}M")
        f3.metric("Monthly EBITDA", f"${ebitda_override:.1f}M")
        f4.metric("Payback", f"{sic['capex_M'] / max(0.01, ebitda_override * 12):.1f} yr")

        st.markdown("### Stage-by-Stage CAPEX")
        capex_data = pd.DataFrame({
            "Stage": ["Acheson", "PVT Growth", "Wafering", "SEG", "AI Platform", "Infrastructure"],
            "CAPEX ($M)": [16, 176, 33, 26, 5, 75],
            "TRL": ["9", "8", "9", "4-5", "N/A", "N/A"]
        })
        fig_capex = px.bar(capex_data, x="Stage", y="CAPEX ($M)", color="Stage", text="CAPEX ($M)", height=300)
        fig_capex.update_layout(showlegend=False, margin=dict(t=10, b=30))
        st.plotly_chart(fig_capex, use_container_width=True)

        st.markdown("### Unit Economics")
        ue1, ue2, ue3 = st.columns(3)
        ue1.metric("Cost/Wafer", f"${sic['cost_per_wafer']:.0f}")
        ue2.metric("CO2/Wafer", f"{sic['co2_per_wafer_kg']:.1f} kg")
        ue3.metric("Annual Revenue", f"${rev_override*12:.0f}M")

# ============================================
# TAB 3: LCA & ENVIRONMENTAL
# ============================================
with tab3:
    st.subheader("Life-Cycle Assessment Comparison")

    # Compare all methods
    lca_compare = []
    for mk, mv in CVD_METHODS.items():
        T_mid = int(np.mean(mv["T_range"]))
        g_tmp = simulate_graphene_cvd(T_mid, 0.5, 0.015, 30, 2, mk, "Methane (CH4)", "Cu Polycrystalline", 15)
        e_tmp = compute_economics({"t_min": 30}, g_tmp, mk, "Methane (CH4)", "Cu Polycrystalline", "Saudi Arabia", True)
        l_tmp = compute_lca(mk, "Methane (CH4)", "Saudi Arabia", e_tmp["E_kWh"], e_tmp["co2_feed"], e_tmp["co2_energy"], g_tmp["wafer_area_cm2"])
        lca_compare.append({"Method": mk.split("(")[0].strip(), "GWP (g CO2/cm2)": l_tmp["gwp_g_cm2"],
                           "Energy (Wh/cm2)": l_tmp["ced_Wh_cm2"], "Water (mL/cm2)": l_tmp["water_mL_cm2"]})

    # Add Si reference
    lca_compare.append({"Method": "Si CZ (reference)", "GWP (g CO2/cm2)": 3.2, "Energy (Wh/cm2)": 85, "Water (mL/cm2)": 18})

    lca_df = pd.DataFrame(lca_compare)

    lc1, lc2 = st.columns(2)
    with lc1:
        fig_gwp = px.bar(lca_df, x="Method", y="GWP (g CO2/cm2)", color="Method", title="Global Warming Potential", height=350)
        fig_gwp.update_layout(showlegend=False, margin=dict(t=40, b=60), xaxis_tickangle=-30)
        st.plotly_chart(fig_gwp, use_container_width=True)
    with lc2:
        fig_ced = px.bar(lca_df, x="Method", y="Energy (Wh/cm2)", color="Method", title="Cumulative Energy Demand", height=350)
        fig_ced.update_layout(showlegend=False, margin=dict(t=40, b=60), xaxis_tickangle=-30)
        st.plotly_chart(fig_ced, use_container_width=True)

    st.subheader("Regional Comparison: Graphene CVD CO2 Footprint")
    regions_data = []
    for rg in ELECTRICITY:
        e_rg = compute_economics({"t_min": 30}, g, method, feedstock, substrate, rg, rg == "Saudi Arabia")
        l_rg = compute_lca(method, feedstock, rg, e_rg["E_kWh"], e_rg["co2_feed"], e_rg["co2_energy"], g["wafer_area_cm2"])
        regions_data.append({"Region": rg, "GWP (g CO2/cm2)": l_rg["gwp_g_cm2"], "Cost/Wafer ($)": e_rg["total_cost"]})
    reg_df = pd.DataFrame(regions_data)
    fig_reg = px.bar(reg_df, x="Region", y=["GWP (g CO2/cm2)", "Cost/Wafer ($)"], barmode="group", height=300)
    st.plotly_chart(fig_reg, use_container_width=True)

# ============================================
# TAB 4: TEA
# ============================================
with tab4:
    st.subheader("Techno-Economic Analysis: Investment Scenarios")

    tea1, tea2 = st.columns([1, 2])
    with tea1:
        st.markdown("#### Graphene Business")
        gr_capex = st.slider("Graphene CAPEX ($M)", 5, 100, 35, 5)
        gr_rev = st.slider("Annual Revenue ($M)", 10, 200, 80, 5)
        gr_opex = st.slider("Annual OPEX ($M)", 5, 50, 15, 2)
        gr_wacc = st.slider("Discount Rate (WACC %)", 5, 20, 10) / 100
        gr_years = st.slider("Project Horizon (years)", 5, 20, 10)

    gr_npv = npv_irr_analysis(gr_capex, gr_rev, gr_opex, gr_years, gr_wacc)

    with tea2:
        n1, n2, n3 = st.columns(3)
        n1.metric("NPV", f"${gr_npv['npv_M']:.0f}M", delta=f"{'Positive' if gr_npv['npv_M'] > 0 else 'Negative'}")
        n2.metric("IRR", f"{gr_npv['irr_pct']:.1f}%")
        n3.metric("Payback", f"{gr_npv['payback_yr']:.1f} yr")

        # Sensitivity tornado
        st.subheader("Sensitivity Analysis (NPV Impact)")
        base_npv = gr_npv["npv_M"]
        sensitivities = []
        for var, delta, fn in [
            ("Revenue +20%", 0.2, lambda: npv_irr_analysis(gr_capex, gr_rev*1.2, gr_opex, gr_years, gr_wacc)["npv_M"]),
            ("Revenue -20%", -0.2, lambda: npv_irr_analysis(gr_capex, gr_rev*0.8, gr_opex, gr_years, gr_wacc)["npv_M"]),
            ("OPEX +30%", 0.3, lambda: npv_irr_analysis(gr_capex, gr_rev, gr_opex*1.3, gr_years, gr_wacc)["npv_M"]),
            ("OPEX -30%", -0.3, lambda: npv_irr_analysis(gr_capex, gr_rev, gr_opex*0.7, gr_years, gr_wacc)["npv_M"]),
            ("CAPEX +50%", 0.5, lambda: npv_irr_analysis(gr_capex*1.5, gr_rev, gr_opex, gr_years, gr_wacc)["npv_M"]),
            ("WACC +5%", 0.05, lambda: npv_irr_analysis(gr_capex, gr_rev, gr_opex, gr_years, gr_wacc+0.05)["npv_M"]),
        ]:
            sensitivities.append({"Variable": var, "NPV ($M)": fn(), "Delta": fn() - base_npv})

        sens_df = pd.DataFrame(sensitivities)
        fig_sens = px.bar(sens_df, x="Delta", y="Variable", orientation="h", color="Delta",
                         color_continuous_scale=["#e74c3c", "#f0f0f0", "#2ecc71"], height=280,
                         title="NPV Change from Baseline")
        fig_sens.update_layout(margin=dict(t=40, b=30))
        st.plotly_chart(fig_sens, use_container_width=True)

# ============================================
# TAB 5: MULTI-PRODUCT COMPARISON
# ============================================
with tab5:
    st.subheader("Multi-Product Portfolio Comparison")
    st.caption("Side-by-side comparison of carbon semiconductor product lines")

    products = [
        {"Product": "Graphene (Pristine)", "Market 2030 ($M)": "200-500", "Aramco CAPEX ($M)": "8-15",
         "OPEX ($/unit)": "$0.02/cm2", "ASP": "$15-80/cm2", "TRL": "5-7", "Timeline": "2027-2029",
         "Key Risk": "Yield consistency"},
        {"Product": "Graphene (Low-Grade)", "Market 2030 ($M)": "500-1500", "Aramco CAPEX ($M)": "20-40",
         "OPEX ($/unit)": "$0.007/cm2", "ASP": "$0.5-5/cm2", "TRL": "7-8", "Timeline": "2026-2028",
         "Key Risk": "Price competition"},
        {"Product": "CNT (Multi-Wall)", "Market 2030 ($M)": "800-1200", "Aramco CAPEX ($M)": "15-30",
         "OPEX ($/unit)": "$20-60/kg", "ASP": "$50-200/kg", "TRL": "7-8", "Timeline": "2026-2027",
         "Key Risk": "Commoditization"},
        {"Product": "GNR (Nanoribbons)", "Market 2030 ($M)": "25-61", "Aramco CAPEX ($M)": "5-10",
         "OPEX ($/unit)": "$500-2000/g", "ASP": "$5k-50k/g", "TRL": "3-4", "Timeline": "2030-2033",
         "Key Risk": "Width control"},
        {"Product": "SiC Wafer (6-in)", "Market 2030 ($M)": "5000-8000", "Aramco CAPEX ($M)": "300-330",
         "OPEX ($/unit)": "$520/wafer", "ASP": "$800/wafer", "TRL": "8-9", "Timeline": "2028-2030",
         "Key Risk": "Crystal yield"},
        {"Product": "SEG on SiC", "Market 2030 ($M)": "50-200 (est.)", "Aramco CAPEX ($M)": "330",
         "OPEX ($/unit)": "$520/wafer", "ASP": "$1200/wafer", "TRL": "4-5", "Timeline": "2032+",
         "Key Risk": "Market existence"},
    ]

    st.dataframe(pd.DataFrame(products), use_container_width=True, hide_index=True)

    st.subheader("Aramco's Phased Entry Strategy (from PPT)")
    phase_data = pd.DataFrame({
        "Phase": ["Phase 1: R&D", "Phase 2: Pilot", "Phase 3: Scale"],
        "Timeline": ["2026-2028", "2029-2031", "2032-2035"],
        "CAPEX ($M)": [15, 63, 160],
        "Target Revenue ($M/yr)": [0, 80, 300],
        "Kill Point": ["Yield <5% at 2027 Q4", "Revenue <$30M at 2030 Q4", "N/A"]
    })
    st.dataframe(phase_data, use_container_width=True, hide_index=True)

    fig_phase = px.bar(phase_data, x="Phase", y="CAPEX ($M)", color="Phase",
                      text="CAPEX ($M)", title="Investment by Phase", height=300)
    fig_phase.update_layout(showlegend=False)
    st.plotly_chart(fig_phase, use_container_width=True)

# Footer
st.divider()
st.caption("Petro-Carbon to Carbon Semiconductor Digital Twin v1.0 | Aramco Advanced Materials R&D | Physics: Vlassiouk nucleation, Cancado defects, JMAK kinetics, Celebi CMC, Percolation transport | Process data: SKKU-SPMDL PPT (April 2026)")
