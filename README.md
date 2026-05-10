# Petro-Carbon to Low-Dimensional Carbon Semiconductor: Digital Twin

A Streamlit-based process simulator and techno-economic analysis tool for evaluating the production of graphene, CNTs, GNRs, SiC wafers, and epitaxial graphene (SEG) from petroleum-derived carbon feedstocks.

## Overview

This Digital Twin provides interactive simulation of:

- **Graphene CVD Process**: Physics-based growth model (v2.1) with 5 validated kinetic sub-models
- **SiC/SEG Process Chain**: 4-stage industrial process from pet-coke to epitaxial graphene on SiC
- **Life-Cycle Assessment (LCA)**: ISO 14040-aligned environmental impact comparison
- **Techno-Economic Analysis (TEA)**: NPV, IRR, payback, and sensitivity analysis
- **Multi-Product Portfolio**: Side-by-side comparison of 6 carbon semiconductor product lines

## Physics Engine (v2.1)

The graphene CVD simulator integrates five peer-reviewed kinetic models:

1. **Vlassiouk nucleation** (J. Phys. Chem. C, 2013): Pressure-dependent activation energy, dual-regime temperature model
2. **JMAK coverage** (Johnson-Mehl-Avrami-Kolmogorov): 2D island coalescence kinetics
3. **Celebi CMC** (J. Mater. Chem. C, 2014): Critical methane concentration for growth/etch boundary
4. **Cancado defect density** (Nano Lett., 2011): Raman I\_D/I\_G to inter-defect distance conversion
5. **Percolation transport**: Coverage-dependent sheet resistance with 2D percolation threshold

## Quick Start

### Local

```bash
git clone https://github.com/YOUR_USERNAME/petro-carbon-digital-twin.git
cd petro-carbon-digital-twin
pip install -r requirements.txt
streamlit run app.py
```

### Streamlit Community Cloud

1. Push this repository to GitHub
2. Go to [share.streamlit.io](https://share.streamlit.io)
3. Connect your GitHub account
4. Select this repository, branch `main`, and file `app.py`
5. Click **Deploy**

## Project Structure

```
petro-carbon-digital-twin/
|-- app.py                  # Main Streamlit application
|-- requirements.txt        # Python dependencies
|-- .streamlit/
|   |-- config.toml         # Streamlit theme and server config
|-- README.md               # This file
|-- LICENSE                  # MIT License
|-- .gitignore              # Git ignore rules
```

## Configuration Parameters

### CVD Energy Sources
| Method | Temperature Range | TRL | Key Advantage |
|--------|------------------|-----|---------------|
| Thermal CVD (Hot-Wall) | 900-1100 C | 7-8 | Highest quality, most mature |
| PECVD | 300-700 C | 5-6 | Low temperature, flexible substrates |
| Solar-Thermal CVD | 800-1200 C | 3-4 | Near-zero electricity cost |
| Joule Heating R2R | 900-1050 C | 5-7 | Highest throughput |
| EM Induction Heating | 1000-1400 C | 4-5 | Best crystal quality |

### Carbon Feedstocks
| Feedstock | C Yield | Aramco Cost | Market Cost |
|-----------|---------|-------------|-------------|
| Methane (CH4) | 75% | $0.18/kg | $0.50/kg |
| Ethylene (C2H4) | 85.7% | $0.45/kg | $1.00/kg |
| Propane (C3H8) | 81.8% | $0.22/kg | $0.55/kg |
| Toluene (C7H8) | 91.3% | $0.40/kg | $0.90/kg |

### SiC/SEG Process (from SKKU-SPMDL analysis)
| Stage | Process | Temperature | CAPEX |
|-------|---------|------------|-------|
| 1 | Acheson (Pet-Coke to SiC) | 2,200 C | $16M |
| 2 | PVT (SiC to Boule) | 2,050 C | $176M |
| 3 | Wafering (Boule to Wafer) | Room T | $33M |
| 4 | SEG (SiC to Graphene) | 1,400 C | $26M |

## Known Limitations

1. **Growth rate calibration (RF-1)**: The radial growth rate is empirically calibrated but not fitted to a comprehensive experimental dataset. Above T\_crossover (~1050 C), nucleation density drops sharply, which is physically correct but may overestimate the growth difficulty at high temperatures.

2. **I\_D/I\_G coefficients (RF-2)**: The defect model uses semi-empirical coefficients. While physically structured (4 independent terms), the individual coefficients have not been systematically fitted to a large Raman dataset.

3. **Feedstock-specific kinetics**: The simulator uses methane decomposition parameters for all feedstocks. Ethylene, propane, and toluene have different decomposition pathways that require experimental validation.

4. **SiC/SEG economics**: Based on the SKKU-SPMDL presentation (April 2026). Revenue projections assume SEG wafer pricing of $1,200/wafer, which is an estimate for a market that does not yet exist.

## References

- Vlassiouk, I. et al. *J. Phys. Chem. C* **117**, 18919 (2013)
- Cancado, L.G. et al. *Nano Lett.* **11**, 3190 (2011)
- Celebi, K. et al. *J. Mater. Chem. C* **2**, 2277 (2014)
- Kim, H. et al. *ACS Nano* **6**, 3614 (2012)
- Amontree, J. et al. *Nature* **630**, 636 (2024)
- Alghfeli, A. & Fisher, T.S. *Sci. Rep.* **14**, 3660 (2024)
- Zhao, J. et al. *Nature* **625**, 60 (2024) [SEG on SiC]

## License

MIT License. See [LICENSE](LICENSE) for details.

## Acknowledgments

Developed as part of the Aramco Advanced Materials R&D program in collaboration with SKKU Smart Process & Materials Design Lab (SPMDL).
