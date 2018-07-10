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

- circularVAR_p_.stan : optim関数などの最適化手法と同等の結果を得られるようにすれば, 各パラメータに対して必要なpriorがわかるかも
- circularVAR_noconst.stan : 定数項を除いたもの, まあまあ綺麗に書ける
- circularVAR_penalty.stan : うまく制約を入れれば機能する??
- circularVAR_sigma_I.stan : det(sigma) = 1 の制約のみ,パラメータの増大防げる 
# Tree

```
├── README.md
├── circular regression.Rproj
├── circular_reg_analysis.Rmd (実データ分析用)
├── circular_reg_analysis.html 
├── circular_reg_optim.Rmd (実データ分析用, optim関数, Rsolnp関数, 制約付き実行可能)
├── circular_reg_sim.Rmd (simulation data 分析用, 作成途中)
├── circular_reg_stan.Rmd (MCMCによる推定, VAR(p) model まで一通り作成)
├── csv
│   └── nerima_wind.csv
├── fit (Not push)
├── model (Not push)
│   ├── circularVAR_p_.rds
│   ├── circularVAR_p_noconst.rds
│   ├── circularVAR_p_penalty.rds
│   ├── circularVAR_p_sigma_I.rds
│   ├── circular_AR.rds
│   ├── obs_circular_VAR.rds
│   ├── state_circular_VAR.rds
│   ├── state_circular_atanh_VAR.rds
│   ├── state_obs_circular_VAR.rds
│   └── state_obs_circular_atanh_VAR.rds 
├── script
│   ├── PN_MLE.R (Rsolnp package による実装, 実行不可)
│   ├── VAR.R (Vector auto regressive関数, 自作関数を作成する参考になるかも)
│   └── function.R (確率密度関数等, 種々の関数保存ファイル)
└── stan
    ├── circularVAR_p_.stan (VAR(p) model, 全てのパラメータの事前分布を正規分布とする. constに関するパラメータの事前分布をN(0,10), Simgaに関するパラメータの事前分布を小さく設定する)
    ├── circularVAR_p_penalty.stan (平均ベクトルが, 単位円周上付近に存在するとして罰則項を与える)
    ├── circularVAR_p_sigma_I.stan (det(Sigma)=1として, 制約を置く, その他のパラメータは大きめの正規分布にしてみる) 
    ├── circular_AR.stan (AR model 角度データだけでなく, 実空間上の値にも適用可能) 
    ├── circular_AR_.stan (アヒル本を参考に, 予測曲線を描く式も導入)
    ├── observation_circular_VAR.stan
    ├── state_circular_VAR.stan
    ├── state_circular_atanh_VAR.stan
    ├── state_obs_circular_VAR.stan
    ├── state_obs_circular_atanh_VAR.stan
    └── test.stan  (stan codeをR内で使用する際に, 関数式として使用)    
```

