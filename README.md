# circular regression

Briefly summarize the circular regression and the circular correlation.

# References
L. Dryden, T.Kent, Geometry Driven Statistics, Parameteric circular-circular regression and diagnostic analysis, pp 115-128

R. Jammalamadaka, TOPICS IN CIRCULAR STATISTICS, Circular Correlation and regression, pp 175-203 

# purpose 

We introduce a new approach to a circular time series modeling that based on the projected normal distributions with circular-circular regression. 
 
# Tree

```
├── README.md
├── circular regression.Rproj
├── circular_reg_analysis.Rmd (実データ分析用)
├── circular_reg_analysis.html 
├── circular_reg_optim.Rmd (実データ分析用, optim関数, 先生作)
├── circular_reg_sim.Rmd (simulation data 分析用, 作成途中)
├── circular_reg_stan.Rmd (MCMCによる推定, VAR(p) model まで一通り作成)
├── circular_test.Rmd 
├── circular_test.html
├── csv
│   └── nerima_wind.csv
├── fit
├── model
├── script
│   ├── PN_MLE.R (Rsolnp package による実装, 実行不可)
│   ├── VAR.R 
│   └── function.R (確率密度関数等, 種々の関数保存ファイル)
└── stan
    ├── circularVAR_1.stan (VAR(1) model)
    ├── circularVAR_p.stan (VAR(p) model)
    ├── test.stan
    └── test1.stan
```
