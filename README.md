# circular regression
Briefly summarize the circular regression and the circular correlation.

# References
L. Dryden, T.Kent, Geometry Driven Statistics, Parameteric circular-circular regression and diagnostic analysis, pp 115-128

R. Jammalamadaka, TOPICS IN CIRCULAR STATISTICS, Circular Correlation and regression, pp 175-203 

# purpose 
We introduce a new approach to a circular time series modeling that based on the projected normal distributions with circular-circular regression. 

```
Circular_regression/  
├── README.md
├── circular\ regression.Rproj
├── circular_reg_analysis.Rmd 
├── circular_reg_optim.Rmd (from 塩濱先生, optim関数による最適化)
├── circular_reg_sim.Rmd (simulationによる, モデルの妥当性検討)
├── circular_test.Rmd 
├── csv
│   └── nerima_wind.csv
├── script
│   ├── PN_MLE.R
│   ├── PN_stan.R
│   └── function.R
└── stan
    ├── circularVAR.stan
    ├── test.stan
    └── test1.stan

```
