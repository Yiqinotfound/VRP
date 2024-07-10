为了科研，需要恶补一些关于VRP的知识!

# CVRP

## Three Index Formulation

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

- 用一个有向图 $G=(V,A)$，建模，$N=V\setminus \{ 0,n+1 \}$ 是 customer 节点
- $q_{i}$ 是 **demand**, $s_{i}$ 是 **service time**
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


其中 $M_{ij}$ 是一个足够大的常数，可以设成$\max\{ b_{i}+s_{i}+t_{ij}-a_{j},0 \}$




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
