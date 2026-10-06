# Replication Code: Nonparametric Regression via a Data-Driven Gram-Schmidt Method

This repository contains the R replication code and data snapshot for the thesis:
> *"Nonparametric Regression via a Data-Driven Gram-Schmidt Method"*

## Requirements

The analysis was executed in R (version 4.5.3). The required packages can be installed via:

```r
install.packages(c("quantmod", "sandwich", "ggplot2", "tikzDevice"))
```

## How to Run

1. Clone or download this repository.
2. Set your R working directory to the project root.
3. Run the scripts:
   - `empirical_analysis.R`: Runs the empirical asset pricing models (AAPL & NVDA), degree selection (CV/GCV), Newey–West HAC inference, and Wald linearity tests.
   - `mc_simulation.R`: Runs the Monte Carlo simulation ($M = 1000$ replications) comparing the Raw, Theoretical, and Empirical polynomial representations.

All tables (`.csv`) and figures (`.tex`) are automatically saved to the `output/` directory.

## Data

The empirical application uses daily return data for AAPL, NVDA, and the S&P 500 from July 5, 2016, to June 30, 2026 ($n = 2511$).  
The cleaned data is provided in `output/market_returns_2016_2026.csv` so the entire pipeline runs offline without depending on live API connectivity.
