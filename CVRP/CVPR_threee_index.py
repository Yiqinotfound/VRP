# 用 three_index 模型求解CVPR问题,不用考虑时间窗
import pandas as pd
import numpy as np
from typing import List
from gurobipy import Model, GRB, quicksum
from matplotlib import pyplot as plt
import scienceplots


class Node:
    _index_counter = 0

    def __init__(
        self,
        x: float,
        y: float,
        demand: float,
        ready_time: float,
        due: float,
        service_time: float = 0.0,
    ):
        self.index = Node._index_counter
        Node._index_counter += 1

        self.x: float = x
        self.y: float = y
        self.demand: float = demand
        self.ready_time: float = ready_time
        self.due: float = due
        self.service_time: float = service_time


class DiGraph:
    def __init__(self, customer_nodes: List[Node], depot: Node):
        self.nodes: List[Node] = [depot] + customer_nodes + [depot]
        self.depot = depot
        self.customer_num = len(customer_nodes)
        self.node_num = len(self.nodes)

        # 计算邻接矩阵
        self.cost_matrix = np.zeros((self.node_num, self.node_num))
        self.cal_adj_matrix()

    def cal_adj_matrix(self):
        for i in range(len(self.nodes)):
            for j in range(len(self.nodes)):
                self.cost_matrix[i][j] = np.sqrt(
                    (self.nodes[i].x - self.nodes[j].x) ** 2
                    + (self.nodes[i].y - self.nodes[j].y) ** 2
                )

    def cost(self, i: int, j: int):
        return self.cost_matrix[i][j]


def initialize_graph(data_path: str, customer_num: int = 15):
    data = pd.read_csv(data_path, sep="\s+")

    customer_nodes = []
    
    # 将第一个customer node作为depot
    depot = Node(
        data["XCOORD."][0],
        data["YCOORD."][0],
        data["DEMAND"][0],
        data["READY_TIME"][0],
        data["DUE_DATE"][0],
        data["SERVICE_TIME"][0],
    )
    for i in range(1, customer_num + 1):
        customer_nodes.append(
            Node(
                data["XCOORD."][i],
                data["YCOORD."][i],
                data["DEMAND"][i],
                data["READY_TIME"][i],
                data["DUE_DATE"][i],
                data["SERVICE_TIME"][i],
            )
        )

    return DiGraph(customer_nodes, depot)


class CVPR:
    
    class Route:
        def __init__(self,nodes:List[Node],cost):
            self.nodes = nodes
            self.cost = cost
    
    def __init__(self, digraph: DiGraph, vehicle_num: int, vehicle_capacity: int):
        self.digraph = digraph
        self.vehicle_num = vehicle_num
        self.vehicle_capacity = vehicle_capacity
        self.print_problem_information()

    def solve(self):

        o = 0  # depot
        d = self.digraph.node_num - 1  # end

        # 总节点列表
        V = list(range(self.digraph.node_num))

        # 客户节点列表,不包括depot和end
        N = V[1:-1]
        
        # i,j \in A时的集合
        I_set = V[:-1]
        J_set = V[1:]

        # 车辆列表
        K = list(range(self.vehicle_num))

        # 创建模型
        m = Model("CVRP_three_index")

        # 设置GAP
        m.setParam("MIPGap", 0.05)

        # 创建变量

        ## x_ijk,车k是否从i点到j点 i \in V[:-1], j \in V[1:], k \in K
        x = m.addVars(I_set, J_set, K, vtype=GRB.BINARY, name="x")

        ## y_ik,车k是否经过节点 i \in V, k \in K
        y = m.addVars(V, K, vtype=GRB.BINARY, name="y")

        ## u_ik,车k在i点的accumaleted demand i \in V, k \in K
        u = m.addVars(V, K, vtype=GRB.CONTINUOUS, name="u")

        # 目标函数
        m.setObjective(
            quicksum(
                self.digraph.cost(i, j) * x[i, j, k]
                for i in I_set
                for j in J_set
                for k in K
            ),
            GRB.MINIMIZE,
        )

        # 设置约束

        ## 每个客户节点只被访问一次
        m.addConstrs(quicksum(y[i, k] for k in K) == 1 for i in N)

        ## path flow constraint
        m.addConstrs(quicksum(x[o, j, k] for j in J_set if j != o) == 1 for k in K)  # depot 节点,出度-入度=1
        m.addConstrs(
            quicksum(x[i, j, k] for j in J_set if j != i)
            - quicksum(x[j, i, k] for j in I_set if j != i)
            == 0
            for i in N
            for k in K
        )  # customer 节点，出度-入度=0

        ## y和x的关系
        m.addConstrs(
            y[i, k] == quicksum(x[i, j, k] for j in J_set if j != i)
            for i in I_set
            for k in K
        )  # 非end节点
        m.addConstrs(
            y[d, k] == quicksum(x[i, d, k] for i in I_set if i != d) for k in K
        )  # end节点

        ## 容量+MTZ约束
        Q = self.vehicle_capacity
        m.addConstrs(
            u[i, k] - u[j, k] + Q * x[i, j, k] <= Q - self.digraph.nodes[j].demand
            for i in I_set
            for j in J_set
            if j != i
            for k in K
        )
        m.addConstrs(u[i, k] >= self.digraph.nodes[i].demand for i in V for k in K)
        m.addConstrs(u[i, k] <= Q for i in V for k in K)

        # 求解
        m.optimize()

        self.objective = m.objVal
        self.routes = self.get_routes(m)

    def get_routes(self, m):
        routes = []
        for k in range(self.vehicle_num):
            route = []
            i = 0
            route.append(i)
            while True:
                if i == self.digraph.node_num - 1:
                    break 
                for j in range(1,self.digraph.node_num):
                    if j != i and m.getVarByName(f"x[{i},{j},{k}]").x > 0.5:
                        route.append(j)
                        i = j
                        break
            routes.append(route)
        return routes
                

    def cal_route_distance(self, route):
        distance = 0
        for i in range(len(route) - 1):
            distance += self.digraph.cost(route[i], route[i + 1])
        return distance

    def plot_solution(self):
        plt.style.use(["science"])
        plt.figure(figsize=(10, 10))
        plt.scatter(
            [node.x for node in self.digraph.nodes],
            [node.y for node in self.digraph.nodes],
            c="r",
            label="customer",
        )
        plt.annotate("depot", (self.digraph.nodes[0].x, self.digraph.nodes[0].y))
        for i in range(1, self.digraph.node_num - 1):
            plt.annotate(
                self.digraph.nodes[i].index,
                (self.digraph.nodes[i].x, self.digraph.nodes[i].y),
            )
        labeld = 0
        for route in self.routes:
            if labeld == 0:
                plt.plot(
                    [self.digraph.nodes[i].x for i in route],
                    [self.digraph.nodes[i].y for i in route],
                    marker="o",
                    label="route",
                )
                labeld = 1
            x = [self.digraph.nodes[i].x for i in route]
            y = [self.digraph.nodes[i].y for i in route]
            plt.plot(x, y, marker="o")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title(f"CVPR With {self.vehicle_num} Vehicles, {self.vehicle_capacity} Capacity, {self.digraph.customer_num} Customers")
        plt.legend()
        plt.savefig(f"CVRP\Figs\CVPR_3_index With {self.vehicle_num} Vehicles, {self.vehicle_capacity} Capacity, {self.digraph.customer_num} Customers.png")
        plt.show()

    def print_problem_information(self):
        print("-" * 20, "Problem Information", "-" * 20)
        print(f"节点总数: {self.digraph.node_num}")
        print(f"客户点总数: {self.digraph.customer_num}")
        print(f"车辆总数: {self.vehicle_num}")
        print(f"车容量: {self.vehicle_capacity}")


def main():
    CUSTOMER_NUM = 5
    VEHICLE_NUM = 4
    VEHICLE_CAPACITY = 100
    digraph = initialize_graph(data_path="dataset/Solomon/R101.txt", customer_num=CUSTOMER_NUM)
    model = CVPR(digraph, vehicle_num=VEHICLE_NUM, vehicle_capacity=VEHICLE_CAPACITY)
    model.solve()
    model.plot_solution()

if __name__ == "__main__":
    main()


