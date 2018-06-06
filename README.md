# circular regression

Briefly summarize the circular regression and the circular correlation.

# References
L. Dryden, T.Kent, Geometry Driven Statistics, Parameteric circular-circular regression and diagnostic analysis, pp 115-128

R. Jammalamadaka, TOPICS IN CIRCULAR STATISTICS, Circular Correlation and regression, pp 175-203 

# purpose 

We introduce a new approach to a circular time series modeling that based on the projected normal distributions with circular-circular regression. 
 
# memo 
 
- 時系列モデルの構造をPN2(mu,Sigma)と仮定した時, どちらかの変数になんらかの制約を入れる入れる必要がある. パラメータが爆発的に増加する
- 時系列モデルの平均のような, 中心線はalpha_0の角度 (得られた値を消去し, c(0,0)として計算してみた)

# Penalty function method 

`circularVAR_p_penalty.stan`におけるpenalty function について記載する.

罰則項一覧: 
- fabs(mu' * mu - 1) : ノルムの絶対値を1に近づける
- square(mu' * mu - 1) : ノルムの二乗誤差を0に近づける
- (mu' * mu) : パラメーターを0付近にする, Lasso 

入力の`lambda`を0とすれば, 罰則なしとして計算可能

# stan code

- circularVAR_p.stan : optim関数などの最適化手法と同等の結果を得られるようにすれば, 各パラメータに対して必要なpriorがわかるかも

-

# Tree

```
├── README.md
├── circular regression.Rproj
├── circular_reg_analysis.Rmd (実データ分析用)
├── circular_reg_analysis.html 
├── circular_reg_optim.Rmd (実データ分析用, optim関数, Rsolnp関数, 制約付き実行可能)
├── circular_reg_sim.Rmd (simulation data 分析用, 作成途中)
├── circular_reg_stan.Rmd (MCMCによる推定, VAR(p) model まで一通り作成)
├── circular_test.Rmd 
├── circular_test.html
├── csv
│   └── nerima_wind.csv
├── fit
├── model
│   ├── circularVAR_p_1.rds (circularVAR_p.stanによるstan model) 
│   └── circularVAR_p_2.rds (circularVAR_p_.stanによるstan model)
├── script
│   ├── PN_MLE.R (Rsolnp package による実装, 実行不可)
│   ├── VAR.R (Vector auto regressive関数, 自作関数を作成する参考になるかも)
│   └── function.R (確率密度関数等, 種々の関数保存ファイル)
└── stan
    ├── circularVAR_1.stan (VAR(1) model, errorあり, p次元に拡張したモデルあるので不要)
    ├── circularVAR_p.stan (VAR(p) model, const vectorにnorm1となる制約を付与) 
    ├── circularVAR_p_.stan (VAR(p) model, 全てのパラメータの事前分布を正規分布とする. 平均方向に関するパラメータに-1~1の制約を付与)
    ├── circularVAR_p_norm.stan (VAR(p) model, 平均方向のベクトルのノルムを1にする, 実現不可能なので使わない) 
    ├── circularVAR_p_I.stan ( VAR(p) model, 分散構造を単位行列とする)
    ├── circularVAR_p_penalty.stan (平均ベクトルが, 単位円周上付近に存在するとして罰則項を与える)
    ├── circularVAR_p_noconst.stan (const vectorを c(0,0)とする) 
    ├── test.stan (stan codeをR内で使用する際に, 関数式として使用)
    └── test1.stan (test code 2, 不要)
```

