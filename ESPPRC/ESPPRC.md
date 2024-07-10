# ESPPRC

找到 $p\to d$ 的最小 cost 路径，且每个节点只能访问一遍 

- $G=(V,A)$  包括一个起点 $p$ 和终点 $d$
- $c_{ij}$
- $L$ The number of resources  
- $d_{ij}^{l}\geq0$ The consumption of resource $l$ along $(i,j)$, 满足三角不等式 
- 每个 node $i$ 和 resource $l$ 对应区间 $[a_{i}^{l},b_{i}^{l}]$, 代表从 $p$ 到 $i$ $l$ 消耗的区间, 可以表示 Time Window约束，或者用 $[0,Q]$ 表示 capacity constraints, 消耗量用 $t_{i}^{l}$ 表示

$$
\begin{aligned}
	\min & \quad \sum_{(i,j)\in A}c_{ij}x_{ij}\\
	\text{s.t.} &\quad \sum_{(i,j)\in A}x_{ij} - \sum_{(j,i) \in A}x_{ji} = 0 & \forall i \in V \setminus \{ p,d \}\\
	&\quad \sum_{(p,j) \in A}x_{pj}=1\\ 
	& \quad \sum_{(i,d)\in A}x_{id} = 1 \\ 
	& \quad t_{i}^{l} + d_{ij}^{l} - t_{j}^{l} \leq M(1-x_{ij})  & \forall l \in \{ 1,\dots,L \}, (i,j) \in A\\  
	& \quad a_{i}^{l} \leq t_{i}^{l} \leq b_{i}^{l} & \forall l \in \{ 1,\dots,L \}, i \in V\\ 
	& \quad x_{ij} \in \{ 0,1 \} & \forall (i,j) \in A
\end{aligned}
$$


# 1. Label Correcting Algorithm  for SPPRC 

- $X_{pi}$ $p\to i$ 的路径，associated a state $(T_{i}^{1},\dots,T_{i}^{L})$ (每种资源的消耗量) 和 cost $C(T_{i}^{1},\dots,T_{i}^{L})$
- $R_{i}$ $=(T_{i}^{1},\dots,T_{i}^{L})$, $C_{i}$ $=C(T_{i}^{1},\dots,T_{i}^{L})$


**Dominance Rule**:

$X_{pi}'$ 和 $X_{p i }^{*}$ 是两个不同的 path, 满足以下两个条件 $X_{p i}'$ dominates $X_{pi}^{*}$
- $C_{i}'\leq C_{i}^{*}$ 且 $T_{i}'^{l} \leq T_{i}^{*l}$ 
- $ (R'_{i}, C'_{i}) \neq (R^{*}_{i},C^{*}_{i}) $

Dominance Rule 的意义就是，$X_{p i }'$ extend 出来的 $X_{p j}'$ 必定 dominate 或等于 $X_{p j }^{*}$, 是显然的




# 2. Adaption to ESPPRC  




**Unreachable Node**


A node $v_{k}$ is said to be unreachable if it is included in $X_{p i}$ or if there exists a resource $l \in \{ 1,\dots,L \}$ satisfying $T_{i}^{l}+d_{ik}^{l} > b_{k}^{l}$ 

**New Definition of Label** 


$R_{i}=(T_{i}^{l},\dots,T_{i}^{L},s_{i},V_{i}^{1},\dots,V_{i}^{n})$ ,对应资源的消耗量，unreachable nodes 数量，unreachable node 向量

**Dominance Rule**
$X_{pi}'$ 和 $X_{p i }^{*}$ 是两个不同的 path, 满足以下两个条件 $X_{p i}'$ dominates $X_{pi}^{*}$
- $C_{i}'\leq C_{i}^{*}$ 且 $T_{i}'^{l} \leq T_{i}^{*l}$ 且 $s_{i}'\leq s_{i}^{*}$
- $(R'_{i}, C'_{i})\neq (R^{*}_{i},C^{*}_{i})$
- $V_{i}'^{k}\leq V_{i}^{*k}$

>**All nondominated paths are the extension of a nondominated path**

**Claim** 

During the execution of the modified algorithm, we need only to consider **nondominated paths** 

**Proof**


证明很容易：

设 $X_{p i}+(i,j)$ 是一个 $p\to j$ 的feasible elementary path $X_{pj}$ ，只需证明 $X_{p i}'+(i,j)$ 也是一个 $p\to j$ 的feasible elementary path $X_{pj}'$，且 dominates or is equal to $X_{pj}$
首先由
- $C_{i}'\leq C_{i}^{*}$ 且 $ T_{i}'^{l} \leq T_{i}^{*l} $ 且 $ s_{i}'\leq s_{i}^{*} $
- $ (R'_{i}, C'_{i})\neq (R^{*}_{i},C^{*}_{i}) $
知 $X_{pj}'$ 必定可行

易知 $T_{j}'^{l} \leq T_{j}^{l}$, 且 $V_{j}^{'l}=V_{j}^{j}=1$, $C_{j}'\leq C_{j}$, 只需证明 $V_{j}'^{k}\leq V_{j}^k\implies s_{j}^{'}\leq s_{j}$, 只需证明 $X_{p j}'$ reach 不到的 $X_{pj}$ 也 reach 不到，假设 $k$ 是 $X_{pj}'$ 的一个 unreachable node ,有两种情况
1. $k$ is included in $X_{pj}'$  因为 $V_{i}'^{k}=1\implies V_{j}^{k}\geq V_{i}^{k} \geq V_{i}'^{k}=1$
2. 某个资源 $l$ $T_{j}'^{l} + d_{jk}^{l}  >b_{k}^{l}$  $T_{j}^{l} + d_{jk}^{l} >T_{j}'^{l} +d_{jk}^{l}>b_{k}$ 所以 unreachable 

