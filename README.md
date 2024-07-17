<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://user-images.githubusercontent.com/25423296/163456776-7f95b81a-f1ed-45f7-b7ab-8fa810d529fa.png">
  <source media="(prefers-color-scheme: light)" srcset="https://user-images.githubusercontent.com/25423296/163456779-a8556205-d0a5-45e2-ac17-42d089e3c3f8.png">
  <img alt="Shows an illustrated sun in light mode and a moon with stars in dark mode." src="https://user-images.githubusercontent.com/25423296/163456779-a8556205-d0a5-45e2-ac17-42d089e3c3f8.png">
</picture>

来和我一起恶补VRP吧!!

# CVRP

## Three Index Formulation

[CVPR-3-index代码实现](CVRP/CVPR.py)

有向图 $G(V,A)$, depot $0$ 可以用两个点 $o$ 和 $d$ 表示 


$$
V:=N \cup \{ o,d \} \ \text{and}\  A:= (V \setminus \{ d \}) \times(V \setminus \{ o \})
$$

决策变量：
- $x_{ijk}$ 车 $k$ 是否经过 $(i, j)\in A$ 
- $y_{ik}$ 车 $k$ 是否经过点 $i \in V$
- $u_{ik}$ $k$ 车在到达 $i$ 时的 accumulative demand 

同样定义 $q_{o}=q_{d}=0$


$$
\begin{aligned}
	\min & \quad \sum_{k \in K} c^{\top}x_{k}\\
	\text{s.t.} &\quad \sum_{k \in K} y_{ik}=1 & \forall i \in N\\
	&\quad x_{k}(\delta^{+}(i))-x_{k}(\delta^{-}(i))=\left\{ \begin{array}{l}
 1 & i=o  \\
0 & i \in N  
 \end{array} \right. & \forall i \in V \setminus \{ d \}  , k \in K\\ 
&  y_{ik} = x_{k}(\delta^{+}(i))  & \forall i \in V \setminus \{ d \}, k \in K  \\ 
& y_{dk} = x_{k}(\delta^{-}(d)) & \forall k \in K \\ 
& u_{ik} - u_{jk} + Q x_{ijk} \leq Q - q_{j}  & \forall(i,j) \in A, k \in K\\ 
& q_{i} \leq u_{ik}\leq Q & \forall i \in V, k \in K \\ 
& x = (x_{k})\in \{ 0,1 \}^{K \times A}  \\ 
& y = (y_{k})\in \{ 0,1 \}^{K\times V}
\end{aligned}
$$

代码中所用的数据集为Solomon数据集

## Extensive Formulation 

直接对 feasible routes 进行操作

-  $\Omega$: The set of feasible CVPR routes. 每个 $r \in \Omega$ 都是 $(i_{0},\dots,i_{s+1})$ 的形式，其中 $i_{0}=o,i_{s+1}=d$

决策变量
- $a_{ir}\in \{ 0,1 \}$ route $r$ 是否经过 $i \in N$
- $\lambda_{r}$： route $r$ 是否被选中 

$$
\begin{aligned}
	\min & \quad c^{\top}\lambda\\
	\text{s.t.} &\quad \sum_{r \in \Omega}a_{ir}\lambda_{r}=1\\
	&\quad \mathbb{1}^{\top}\lambda=\lvert K \rvert \\ 
	& \lambda \in \{ 0,1 \}^{\Omega}
\end{aligned}
$$

如果 cost 满足三角不等式 $c_{hi}+c_{ij}\geq c_{hi}$，第一个约束可以改成 $\geq 1$
同理如果不需要全部 $\lvert K \rvert$ 个，就用 $\leq \lvert K \rvert$


# VRPTW

Time widow can be **soft or hard**, **one-sided or two-sided**

- 用一个有向图 $G=(V,A)$，建模， $N = V \setminus \{ 0,n+1 \}$ 是 customer 节点
-  $q_{i}$  是 **demand**,  $s_{i}$  是 **service time**
-  

$$
q_{0}=q_{n+1}=s_{0}=s_{n+1}=0
$$


$$
[a_{0},b_{0}]=[a_{n+1},b_{n+1}]
$$


其中 $a_{0},b_{0}$ 是 earliest possible departure time from the depot and the latest possible arrival time at the depot.

假设 travel time 满足三角不等式，feasible solutions exists only if 


$$
\begin{align}
a_{0}\leq& \min_{i \in V\setminus \{ 0 \}}\{ b_{i}-t_{0i} \} \\
b_{0}\geq& \max_{i \in V\setminus \{ 0 \}} \{  \max\{ a_{0}+t_{0i},a_{i} \}+s_{i}+t_{i,n+1}\}
\end{align}
$$


如果 $a_{i}+s_{i}+t_{ij}>b_{i}$ 或者 $q_{i}+q_{j}>Q$, arc $(i,j)$ 就可以直接被省略 

- 如果允许待在 depot

$$
c_{0,n+1}=t_{0,n+1}=0
$$

##  MIP Formulation

[VRPTW代码实现](VRPTW/VRPTW.py)


决策变量
- $x_{ijk}$
- $T_{ik}$


$$
\begin{align}
	\min & \quad \sum_{k \in K} \sum_{(i,j)\in A}c_{ij}x_{ijk} \\
	\text{s.t.} &\quad \sum_{k\in K}\sum_{j \in \delta^{+}(i)}x_{ijk}=1 & \forall i \in N\\ 
	& \quad \sum_{j \in \delta^{+}(0)}x_{0jk}=1 & \forall k \in K \\
	&\quad \sum_{i \in \delta^{-}(j)}x_{ijk} - \sum_{i \in \delta^{+}(j)}x_{jik} = 0 & \forall k \in K,j \in N  \\  \\
& \quad \sum_{i \in \delta^{-}(n+1)}x_{i,n+1,k}=1 & \forall k \in K  \\
& \quad x_{ijk}(T_{ik}+t_{ij}+s_{i}-T_{jk})\leq 0  & \forall k \in K,(i,j) \in A   \\
&\quad a_{i} \leq T_{ik} \leq b_{i} & \forall k \in K, i \in V  \\
& \quad \sum_{i \in N}q_{i} \sum_{j \in \delta^{+}(i)}x_{ijk} \leq Q & \forall k \in K  \\
& \quad x_{ijk}\in \{ 0,1 \} & \forall k \in K,(i,j) \in A
\end{align}
$$

第四条约束可以用 MTZ 线性化

$$
T_{ik}+t_{ij}+s_{i}-T_{jk} \leq (1-x_{ijk})M_{ij} \quad \forall k \in K, (i,j) \in A
$$


其中 $M_{ij}$ 是一个足够大的常数，可以设成  $\max\{ b_{i}+s_{i}+t_{ij}-a_{j},0 \}$




## Set partitioning model

$$
\Omega: \text{The set of all feasible routes}
$$

- $c_{r}$ The cost of route $r$ 
- $a_{ir}\in \{ 0,1 \}$ The number of visits to customer $i$ in route $r \in \Omega$ 
- $y_{r}$ route $r$ 是否被选中 

$$
\begin{aligned}
	\min & \quad \sum_{r\in \Omega}c_{r}y_{r}\\
	\text{s.t.} &\quad \sum_{r\in \Omega}a_{ir}y_{r}=1\\
	&\quad y_{r}\in \{  0,1 \}\\
\end{aligned}
$$

子问题是ESPPRC问题，需要用到[动态规划求解ESPPRC](ESPPRC/ESPPRC.md)的知识


## Branch and Price

[[Branch and Price|这里有简单的介绍]]

### 主问题和列生成

使用 **Set covering Model**, 它的松弛可以得到**更好的下界** 
$$
\begin{aligned}
	\min & \quad \sum_{r \in \Omega}c_{r}y_{r}\\
	\text{s.t.} &\quad \sum_{r \in \Omega}a_{ir}y_{r} \geq 1& \forall i \in N\\
	&\quad \sum_{r \in \Omega}y_{r} \leq U\\
	 & \quad y_{r}  \in \mathbb{N}  & \forall   r \in \Omega
\end{aligned}
$$
注意到这里 $y_{r}$ 不是 $0-1$ 的，这样进行线性松弛的时候就可以避免 $y_{r} \leq 1$ 这个约束
> 显然 $y_{r}\geq 2$ 的时候不会是最优的


这个模型的线性松弛的解**可以转化为MIP model 线性松弛的解**，反过来就不一定了 

由于 $\Omega$ 的大小随着 $n$ 的大小是指数增长的，必须用列生成求解 


**MP**是原问题的线性松弛, 再将 $\Omega$ 限制为 $\Omega_{1} \subseteq \Omega$ 就得到 **RMP** $MP(\Omega_{1})$

$$
\begin{aligned}
 & MP(\Omega_{1})\\
	\min & \quad \sum_{r \in \Omega}c_{r}y_{r}\\
	\text{s.t.} &\quad  \sum_{r \in \Omega_{1}}a_{ir}y_{r} \geq 1 & \forall i \in N\\
	&\quad \sum_{r \in \Omega_{1}}y_{r}\leq \lvert K \rvert  \\
	& \quad y_{r}\geq 0 & \forall r \in \Omega_{1}
\end{aligned}
$$



$MP(\Omega_{1})$ 的对偶 $D(\Omega_{1})$, 第一个约束的对偶变量为 $\alpha_{i}, i \in N$, 第二个约束的对偶变量为 $\lambda_{0}$

$$
\begin{aligned}
& D(\Omega_{1})\\
	\max & \quad \sum_{i \in N} \alpha_{i} + \lvert K \rvert \lambda_{0} \\
	\text{s.t.} &\quad \sum_{i \in N}a_{ir}\alpha_{i} +\lambda_{0} \leq c_{r} & (r \in \Omega_{1}) \\
	&\quad \alpha_{i}  \geq0 & (i \in N)\\ 
	&\quad \lambda_{0}\leq 0
\end{aligned}
$$

如果使用 Set partitioning model, $\alpha$ 就会变成自由变量，会使得收敛变慢


>$\Omega_{1} \subseteq \Omega$, 设 $\lambda ^{*}$ 是 $D(\Omega_{1})$ 的最优解，如果 $\lambda ^{*}$ 对 $D(\Omega)$ 是可行的，那么 $\lambda ^{*}$ 也是 $D(\Omega)$ 的最优解 


考虑 $D(\Omega_{1})$ 的最优解 $\lambda ^{*}$, 如果对 $D(\Omega)$ 可行，那么 $\lambda ^{*}$ 同时解决了 $D(\Omega)$ 和 $MP(\Omega)$, 否则说明若干个从 $\Omega \setminus \Omega_{1}$ 推导出的约束 are violated. 对偶问题有不满足的约束，也就是 $MP$ 问题有 $<0$ 的 reduced cost. 


### 2.1.2. 子问题

**Subproblem** 就是要找到这样的路径 $r \in \Omega \setminus \Omega_{1}$，再将这条路径加入 $\Omega_{1}$ 中重新求解 

设 $\alpha_{i}$ 是约束 $\sum_{r}a_{ir}y_{r}=1, i \in N$ 的对偶变量, 且令 $\alpha_{0}=0$ 那么 $y_{r}$ 的 reduced cost 就是 $c_{r}-\sum_{i \in N}a_{ir}\alpha_{i}-\lambda_{0}^{\ast}$，等价形式为 $\sum_{(i,j) \in A}(c_{ij}-\alpha_{i})x_{ij}$ 

$$
\begin{aligned}
&\quad (SP)\\
	\min & \quad \sum_{(i,j) \in A}(c_{ij}-\alpha_{i})x_{ij}\\
	\text{s.t.} &\quad \sum_{j \in \delta^{+}(0)}x_{0j}=1\\
	&\quad \sum_{i \in \delta^{-}(j)}x_{ij} - \sum_{i \in \delta^{+}(j)}x_{ji} = 0   & \forall j \in N\\ 
	& \sum_{i \in \delta^{-}(n+1)}x_{i,n+1}=1 \\ 
	& x_{ij}(T_{i}+s_{i}+t_{ij}-T_{j})\leq 0  & \forall (i,j) \in A \\ 
	& \quad a_{i} \leq T_{i} \leq b_{i} & \forall i \in V \\ 
	& \quad \sum_{i \in N}q_{i}\sum_{j \in \delta^{+}(i)} x_{ij} \leq Q, \\ 
	& \quad x_{ij}\in \{ 0,1 \} & \forall(i,j)\in A
\end{aligned}
$$


这是一个 ESPPRC 问题

### 2.1.3. Branching scheme 

需要有一个精妙的分支结构，, 步骤如下：
- 求解了 $MP(\Omega_{1})$ 后，选择一个 $(i, j)\in A$, 使得 $0<x_{ij}< 1$, 这里 $x_{ij}$ 可以通过 $y_{r}, r\in \Omega_{1}$ 求得
- 生成两个 branches:
	- 删除这条边 $x_{ij}=0$
	- 保留这条边 $x_{ij}=1$

删除: 把 MP 中包含 $(i,j)$ 的 route 删去, 子问题中 $(i,j)$ 也要被删除，确保不会有包含 $(i,j)$ 的 route 生成
保留：禁止 $(i,l),l\neq j$ 的边

[VRPTW的代码实现](VRPTW/VRPTW.py)
