为了科研，需要学习一些关于VRP的知识

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