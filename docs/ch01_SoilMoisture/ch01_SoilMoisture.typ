#import "@preview/modern-cug-report:0.1.3": *
#show: doc => template(doc, footer: "ModernBEPS.jl", header: "")

#show "{": "{"
#show "}": "}"

= 1 土壤水运动

== 1.1 土壤水运动模型结构

#h(2em)
ModernBEPS采用一维土壤水模型，基于达西定律与Richards方程进行求解。土壤默认划分为5层，每层厚度为$Delta z = [0.05, 0.10, 0.20, 0.40, 1.25] "m"$，土壤水运动模块的主要变量见图#[@fig_SM]与表#[@table_module_SM]。

#figure(
  image("../Figures/Figure_Soil_BEPS.pdf", width: 70%),
  caption: [
    土壤水分层示意图。
    $theta$为体积含水量，$psi$为吸力水头，$K$为水力传导系数，$Q$为层间水分通量。$E$为土壤蒸发，主要发生在土壤的第一层；$T$为植被蒸腾，可能发生在土壤的每一层。
    到达地表的水量一部分下渗（$I$）、一部分形成地表径流（$R_s$）。其中，到达地表的水量指扣除冠层截留等过程后的降水输入。
  ],
) <fig_SM>

#h(2em)
为便于描述模型计算流程，表中沿用 State、Param 和 Output 三类。State 表示随时间更新的状态量，Param 表示模拟期间保持不变或预设的模型量，Output 表示每个时间步计算得到的过程量。

#figure(
  caption: [土壤水运动模块主要变量。],
  table(
    columns: (2cm, 2cm, 3cm, 8cm, 2cm),
    align: (horizon, horizon, horizon, horizon, horizon),
    inset: 4pt,
    rows: 0.75cm,
    // stroke: frame(rgb("21222C")),
    [*类别*], [*符号*], [*模型变量*], [*含义*], [*单位*],
    [State], [$P_g$], [`r_rain_g`], [到达地表的降水或穿透雨速率], [$"m s"^(-1)$],
    [State], [$z_w$], [`z_water`], [地表积水深度], [$"m"$],
    [State], [$theta_i$], [`θ[i]`], [第 $i$ 层体积含水量], [$"m"^3 " m"^(-3)$],
    [Param], [$Delta t$], [`kstep`], [土壤水分更新步长，默认 360 s], [$"s"$],
    [Param], [$r_d$], [`r_drainage`], [未下渗水量的地表滞留比例], [-],
    [Param], [$theta_("sat",i)$], [`θ_sat[i]`], [饱和体积含水量], [$"m"^3 " m"^(-3)$],
    [Param], [$theta_("vwp",i)$], [`θ_vwp[i]`], [凋萎点含水量，下限截断值], [$"m"^3 " m"^(-3)$],
    [Param], [$K_("sat",i)$], [`K_sat[i]`], [饱和水力传导系数], [$"cm h"^(-1)$],
    [Param], [$psi_("sat",i)$], [`ψ_sat[i]`], [饱和吸力水头参数（负值）], [$"cm"$],
    [Param], [$b_i$], [`b[i]`], [Campbell 参数], [-],
    [Output], [$S_i$], [`ETi[i]`], [第 $i$ 层蒸散耗水汇项], [$"m s"^(-1)$],
    [Output], [$Q_i$], [`r_waterflow[i]`], [第 $i$ 与 $i+1$ 层之间的水分通量], [$"m s"^(-1)$],
  ),
) <table_module_SM>

#table-note[
  饱和水力传导系数进入方程前单位从$"cm h"^(-1)$转换为 $"m s"^(-1)$，记为 $K^"SI"$，$K_("sat",i)^"SI" = K_("sat",i) / (3600 times 100)$。
  代码内部$psi_i$与$psi_("sat",i)$按负的压力水头存储，单位为$"cm"$；涉及梯度计算时，层厚$Delta z$需换算到相同单位。
]


== 1.2 地表入渗与径流

=== 1.2.1 最大下渗速率

#h(2em)
根据达西定律，土壤最大入渗能力$I_max$为：
$ I_("max") approx K_("sat",1)^"SI" [1 + (psi_s - psi_1) / (Delta z_1^
"cm")] $ <eq_Imax>

其中，$psi_1$ 是表层吸力水头，$psi_s$ 是地表近饱和吸力水头，二者在模型中均按负压头处理；$Delta z_1$为地表到第一层土壤中心的厚度，且$Delta z_1^"cm" = 100 Delta z_1$。$psi_s approx psi_("sat",1) < 0$ 为地表近饱和吸力水头，$K$ 在入渗能力估计中取表层饱和导水率。

基于Campbell土壤转换函数：

$ psi(theta) = psi_("sat",1) (theta / theta_("sat",1))^(-b_1) $

$ dv(psi, theta)|_(theta = theta_("sat",1)) = - b_1 psi_("sat",1) / theta_("sat",1) $


可得$psi_s - psi_1 approx - b_1 psi_("sat",1) (theta_("sat",1) - theta_1) / theta_("sat",1)$。
代回式：

$
  I_("max") = K_("sat",1)^"SI" [1 - (theta_("sat",1) - theta_1) / Delta z_1^"cm" dot (psi_("sat",1) b_1) / theta_("sat",1)]
$

由于$psi_("sat",1) < 0$，括号中的第二项实际为正值，表示由表层未饱和亏缺产生的毛管吸力增强；表层越干，最大入渗能力越大。

=== 1.2.2 实际下渗速率
#h(2em)
实际入渗速率同时受限于当前水资源量$A$与土壤质地控制下的最大下渗能力$I_max$：

$ I = min(A, I_("max")) $

$ A = max(z_w^t / (Delta t) + P_g, 0) $

其中，$z_w^t$ 是当前地表积水深度，$P_g$ 是到达地表的降水输入，$Delta t$是时间步长。$A$是当前地表可供水量折算得到的供水速率上限，代表一个时间步长内$Delta t$地表积水和降水全部转为下渗对应的最大供水速率。

=== 1.2.3 地表积水

#h(2em)
地表积水$z_w^t$一部分下渗$I$，扣除下渗之后的剩余水量为：

$ D_("ex") = (A - I) Delta t $

剩余水量$D_"ex"$一部分转为地表径流$R_s$，剩余部分为下一时刻的积水深度$z_w^(t+1)$：

$ R_s = D_("ex") (1 - r_d) $

$ z_w^(t+1) = D_("ex") r_d $

其中，`r_drainage`为未下渗水量的地表滞留系数，其值越大，水在地表停留越久。

== 1.3 土壤水层间导水系数与垂向通量

#h(2em)
BEPS 使用 Campbell 形式的保水曲线和非饱和水力传导率：

$
  psi_i = psi_("sat",i) (theta_i / theta_("sat",i))^(-b_i), K_i = K_("sat",i)^"SI" (theta_i / theta_("sat",i))^(2 b_i + 3)
$

第 $i$ 层与第 $i+1$ 层之间的界面导水率为：

$
  K_(i+1/2) = underbrace((K_i psi_i + K_(i+1) psi_(i+1)) / (psi_i + psi_(i+1)), "吸力水头加权平均") dot underbrace((b_i + b_(i+1)) / (b_i + b_(i+1) + 6), "Campbell 非线性修正")
$

两层土壤之间垂向上的水分通量$Q_i$为：

$
  Q_i = K_(i+1/2) ( 1 - (psi_(i+1) - psi_i) / (Delta z_(i+1/2)^"cm")),
  Delta z_(i+1/2)^"cm" = 100 (Delta z_i + Delta z_(i+1))/2
$

其中重力项为 $+1$，而压头梯度项在当前负压头约定下写为 $- (psi_(i+1) - psi_i) / Delta z_(i+1/2)^"cm"$。$Q_i > 0$ 表示向下流动，$Q_i < 0$ 表示向上流动。为避免下层在一个外部步长内超过饱和，向下通量受限于：

$ Q_(i,"max") = ((theta_("sat",i+1) - theta_(i+1)) Delta z_(i+1)) / (Delta t) + S_(i+1) $

$ Q_i arrow.r.long min(Q_i, Q_(i,"max")) $

其中$S_i$为第$i$层的蒸散汇项，详见后续章节_*蒸散汇项*_。

== 1.4 含水量更新

#h(2em)
模型中土壤的上边界通量$Q_0 = I$，下边界通量$Q_N = 0$。再叠加上一个时刻的状态，就可以通过控制方程更新每层的含水量：

$ (dif theta_i) / (dif t) = (Q_(i-1) - Q_i - S_i) / Delta z_i $

显式更新为：

$ theta_1^(t+delta t) = theta_1^t + (I - Q_1 - S_1) (delta t) / (Delta z_1) $

$ theta_i^(t+delta t) = theta_i^t + (Q_(i-1) - Q_i - S_i) (delta t) / (Delta z_i), quad i = 2, dots, N $

每个子步后截断：

$ theta_i arrow.r.long min(max(theta_i, theta_("vwp",i)), theta_("sat",i)) $

该方案以凋萎点作为水分下限，而不是残余含水量。

=== 1.4.1 自适应子步长

外部步长通常为 $Delta t = 360 "s"$。每个子步重新计算 $psi_i$、$K_i$、$K_(i+1/2)$ 和 $Q_i$，并按最大通量选择子步长：

$
  delta t = cases(
    1 "s" & "if" max_i abs(Q_i) > 10^(-5) "m s"^(-1),
    30 "s" & "if" max_i abs(Q_i) > 10^(-6) "m s"^(-1),
    360 "s" & "otherwise"
  )
$

最后一个子步若超过外部步长，则截断到剩余时间。

== 1.5 蒸散汇项

#h(2em)
蒸散汇项$S_i$表示第$i$层的蒸散耗水速率，单位为$"m s"^(-1)$。$S_i$包括植被蒸腾$T_i$与土壤蒸发$E_s$。土壤蒸发主要发生在表层土壤，直接从第一层土壤扣水。植被蒸腾可能发生在每一层土壤。模型先根据气象驱动和冠层状态计算总蒸腾量$T$，再根据每层土壤的水分状态、温度状态和根系分布，将$T$分配到每层得到$T_i$。

$ S_1 = T / rho_w dot w_1 + E_s / rho_w $ <eq_S1>

$ S_i = T / rho_w dot w_i, quad i = 2, dots, N $ <eq_Si>

其中，$T_o$ 和 $T_u$ 为两层植被蒸腾，$E_s$ 为土壤蒸发。$w_i$为第$i$层归一化根系吸水权重，$w_i$是土壤的水文状态、温度状态和根系分布的综合加权。

同时考虑土壤水分和温度限制，胁迫因子可写为：

$ f_("stress",i) = f_(psi,i) f_(T,i) $

令$chi_i = -psi_i > 0$表示吸力水头大小，则水分胁迫函数写为：

$
  f_(psi,i) = cases(
    1 & "if" chi_i <= psi_("min"),
    [1 + ((chi_i - psi_("min")) / psi_("min"))^alpha]^(-1) & "if" chi_i > psi_("min")
  )
$

$
  f_(T,i) = cases(
    1 - exp(-0.02 thin T_("soil",i)^2) & "if" T_("soil",i) > 0,
    0 & "if" T_("soil",i) <= 0
  )
$

其中，$psi_("min")$为植被开始感受到水分胁迫的吸力水头阈值，按正值大小使用；$alpha$控制胁迫随吸力水头增加而增强的速率。落叶阔叶林（DBF）和常绿阔叶林（EBF）取$psi_("min") = 1000 "cm"$、$alpha = 1.5$，其它植被类型取$psi_("min") = 3300 "cm"$、$alpha = 0.4$。对应到更直观的长度单位分别约为$10 "m"$与$33 "m"$，对应压力约为$0.10 "MPa"$与$0.33 "MPa"$。

#beamer-block[
  *$psi_("min")$取值说明*：田间持水量对应的基质势约为$-0.033 "MPa"$，折合吸力水头约$3.4 "m"$；永久凋萎点对应的基质势约为$-1.5 "MPa"$，折合吸力水头约$153 "m"$。多数树木气孔开始关闭时的叶片或根区水势约为$0.66 "MPa"$到$0.74 "MPa"$，完全关闭约为$2.55 "MPa"$到$2.75 "MPa"$，不同植被功能型会有所差异。

  $psi_("min") = 1000 "cm"$或$3300 "cm"$并不等同于气孔关闭水势，而是经验性的土壤水分胁迫启动阈值。该阈值相对偏小，表示土壤尚未接近永久凋萎点之前，就开始降低根系吸水权重和冠层水分可利用性。
]
// 在土壤水分参数中，田间持水量由$theta_("vfc")$表示，永久凋萎点含水量由$theta_("vwp")$表示。BEPS 土壤参数库中$theta_("vfc")$随土壤质地约为$0.09$到$0.40 "m"^3 "m"^(-3)$，$theta_("vwp")$约为$0.03$到$0.32 "m"^3 "m"^(-3)$。前文含水量更新中的下限截断采用$theta_("vwp")$，而水分胁迫函数中的$psi_("min")$用于控制蒸腾吸水权重，两者分别作用于水量边界和植被胁迫过程。

#h(2em)
最后，再根据根系分布进行加权，得到未归一化吸水权重。

$
  a_i = f_("root",i) f_("stress",i), quad
  w_i = a_i / sum_(j=1)^N a_j
$ <eq_wi>

#h(2em)
其中，根系分布权重$f_("root",i)$，通过指数公式计算得到：

$
  f_("root",i) = cases(
    1 - beta^(Z_1) & "if" i = 1,
    beta^(Z_(i-1)) - beta^(Z_i) & "if" 2 <= i < N,
    beta^(Z_(N-1)) & "if" i = N
  ),
  quad Z_i = 100 sum_(j=1)^i Delta z_j
$

// $beta$为根系分布衰减系数，$Z_i$为从地表到第$i$层底部的累计深度（cm）。$beta$不同取值情况下的根系分布见图#[@fig_RootFraction]。将最终计算的$w_i$（式#[@eq_wi]）代回式#[@eq_S1]和#[@eq_Si]，即可得到每层土壤的蒸散汇项$S_i$。

#v(-1em)
#figure(
  image("../Figures/Figure_root_fraction.pdf", width: 100%),
  caption: [
    根系分布示意图。
  ],
) <fig_RootFraction>


// 当$sum_(j=1)^N a_j$过小时，模型会采用下限保护，避免归一化权重不稳定。
// 若不考虑逐层胁迫，则可令$f_("stress",i)=1$，此时$w_i$仅由根系分布决定。
