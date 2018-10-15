# circular regression

- Briefly summarize the circular regression and the circular correlation.

## References

- L. Dryden, T.Kent, Geometry Driven Statistics, Parameteric circular-circular regression and diagnostic analysis, pp 115-128
- R. Jammalamadaka, TOPICS IN CIRCULAR STATISTICS, Circular Correlation and regression, pp 175-203 

## purpose 

We introduce a new approach to a circular time series modeling that based on the projected normal distributions with circular-circular regression. 
 
## memo 

- 時系列モデルの構造をPN2(mu,Sigma)と仮定した時, どちらかの変数になんらかの制約を入れる入れる必要がある. パラメータが爆発的に増加する
- 時系列モデルの平均のような, 中心線はalpha_0の角度 (得られた値を消去し, c(0,0)として計算してみた)

# About files 

## Rmd file

- circular_reg_analysis.Rmd: 気象庁データを活用. 風向データは16方位のため, 角度データに変更を行っている.
- circular_reg_optim.Rmd : 実データ分析用, optim関数, Rsolnp関数, 制約付き実行可能
- circular_reg_stan.Rmd: MCMCによる推定, VAR(p) model まで一通り作成
- circular_reg_sim.Rmd: simulation data 分析用, 作成途中, 要検討

## stan file

`circularVAR_p_penalty.stan`におけるpenalty function について記載する.

罰則項一覧: 
- fabs(mu' * mu - 1) : ノルムの絶対値を1に近づける
- square(mu' * mu - 1) : ノルムの二乗誤差を0に近づける
- (mu' * mu) : パラメーターを0付近にする, Lasso 

入力の`lambda`を0とすれば, 罰則なしとして計算可能

### each codes

- circularVAR_p_.stan : optim関数などの最適化手法と同等の結果を得られるようにすれば, 各パラメータに対して必要なpriorがわかるかも
- circularVAR_noconst.stan : 定数項を除いたもの, まあまあ綺麗に書ける
- circularVAR_penalty.stan : うまく制約を入れれば機能する??
- circularVAR_sigma_I.stan : det(sigma) = 1 の制約のみ,パラメータの増大防げる 

# Directrory

```
├── README.md
├── circular_regression.Rproj
├── Rmd
│   ├── circular_reg_analysis.Rmd
│   ├── circular_reg_optim.Rmd 
│   ├── circular_reg_sim.Rmd 
│   └── circular_reg_stan.Rmd
├── csv
│   └── nerima_wind.csv
├── model
│   ├── circularVAR_p_.rds
│   ├── circularVAR_p_noconst.rds
│   ├── circularVAR_p_penalty.rds
│   └── circularVAR_p_sigma_I.rds
├── script
│   ├── PN_MLE.R
│   ├── SubmitPlot.R
│   ├── VAR.R
│   ├── function.R
│   ├── misc.R
│   └── stan_plot_label_helper.R
└── stan
    ├── circularVAR_p_.stan
    ├── circularVAR_p_noconst.stan
    ├── circularVAR_p_penalty.stan
    ├── circularVAR_p_sigma_I.stan
    ├── circular_AR.stan
    ├── circular_AR_.stan
    ├── observation_circular_VAR.stan
    ├── state_circular_VAR.stan
    ├── state_circular_atanh_VAR.stan
    ├── state_obs_circular_VAR.stan
    ├── state_obs_circular_atanh_VAR.stan
    └── test.stan
```

