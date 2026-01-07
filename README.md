# Exploring new statistical metrics to evaluate the magnitude distribution
Code to reproduce the analysis presented in the article *Exploring new statistical metrics to evaluate the magnitude distribution of earthquake forecasting models*

The analysis is performed in R, while the plotting is done using Python. Specifically, the file 
* `magtest_fun.R` contains the functions to calculate the alternative metrics presented in Section 4.
* `test_magtest.R` contains the functions to analyse the ability of the metrics to identify inconsistencies in the body of the distribution.
* `test_magtest_tapered.R` contains the functions to analyse the ability of the metrics to identify inconsistencies in the tail of the distribution (under/overestimation).
* `plot_results.ipynb` contains the code to produce Figures 1, 5, 6, 7, 8, 9, 10, 11, 12, 13, A1.
* `M_test_new_fun.ipynb`contains the functions used to calculate the resampled M-test for the Swiss and European forecasts.

Figures 4, 15 on the European model have been produced by Marta Han by modifying the code relative to the article [*Towards a harmonized operational earthquake forecasting model for Europe*](https://nhess.copernicus.org/articles/25/991/2025/), while Figures 2, 3, 14 by Leila Mizrahi modifying the code relative to the article [*suiETAS: Developing and Testing ETAS‐Based Earthquake Forecasting Models for Switzerland*](https://pubs.geoscienceworld.org/ssa/bssa/article/114/5/2591/644286/suiETAS-Developing-and-Testing-ETAS-Based)
