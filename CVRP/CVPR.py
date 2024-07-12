# 用 three_index 模型求解CVPR问题,不用考虑时间窗
from gurobipy import Model, GRB, quicksum
from matplotlib import pyplot as plt
import scienceplots
import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)
from env import *


class CVPR:
    def __init__(self, digraph: DiGraph, vehicle_num: int, vehicle_capacity: int):
        self.digraph = digraph
        self.vehicle_num = vehicle_num
        self.vehicle_capacity = vehicle_capacity
        self.print_problem_information()
        self.o = 0  # depot
        self.d = self.digraph.node_num - 1  # end

        # 总节点列表
        self.V = list(range(self.digraph.node_num))

        # 客户节点列表,不包括depot和end
        self.N = self.V[1:-1]

        # i,j \in A时的集合
        self.I_set = self.V[:-1]
        self.J_set = self.V[1:]

        # 车辆列表
        self.K = list(range(self.vehicle_num))

    def violently_solve(self):

        # 创建模型
        m = Model("CVRP_three_index")

        # 设置GAP
        m.setParam("MIPGap", 0.05)

        # 创建变量

        ## x_ijk,车k是否从i点到j点 i \in V[:-1], j \in V[1:], k \in K
        x = m.addVars(self.I_set, self.J_set, self.K, vtype=GRB.BINARY, name="x")

        ## y_ik,车k是否经过节点 i \in V, k \in K
        y = m.addVars(self.V, self.K, vtype=GRB.BINARY, name="y")

        ## u_ik,车k在i点的accumaleted demand i \in V, k \in K
        u = m.addVars(self.V, self.K, vtype=GRB.CONTINUOUS, name="u")

        # 目标函数
        m.setObjective(
            quicksum(
                self.digraph.cost(i, j) * x[i, j, k]
                for i in self.I_set
                for j in self.J_set
                for k in self.K
            ),
            GRB.MINIMIZE,
        )

        # 设置约束

        ## 每个客户节点只被访问一次
        m.addConstrs(quicksum(y[i, k] for k in self.K) == 1 for i in self.N)

        ## path flow constraint
        m.addConstrs(
            quicksum(x[self.o, j, k] for j in self.J_set if j != self.o) == 1 for k in self.K
        )  # depot 节点,出度-入度=1
        m.addConstrs(
            quicksum(x[i, j, k] for j in self.J_set if j != i)
            - quicksum(x[j, i, k] for j in self.I_set if j != i)
            == 0
            for i in self.N
            for k in self.K
        )  # customer 节点，出度-入度=0

        ## y和x的关系
        m.addConstrs(
            y[i, k] == quicksum(x[i, j, k] for j in self.J_set if j != i)
            for i in self.I_set
            for k in self.K
        )  # 非end节点
        m.addConstrs(
            y[self.d, k] == quicksum(x[i, self.d, k] for i in self.I_set if i != self.d) for k in self.K
        )  # end节点

        ## 容量+MTZ约束
        Q = self.vehicle_capacity
        m.addConstrs(
            u[i, k] - u[j, k] + Q * x[i, j, k] <= Q - self.digraph.nodes[j].demand
            for i in self.I_set
            for j in self.J_set
            if j != i
            for k in self.K
        )
        m.addConstrs(u[i, k] >= self.digraph.nodes[i].demand for i in self.V for k in self.K)
        m.addConstrs(u[i, k] <= Q for i in self.V for k in self.K)

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
                for j in range(1, self.digraph.node_num):
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
        # 用箭头画出路径
        cmap = plt.get_cmap("Dark2")  # 使用tab10颜色映射，可选择其他颜色映射
        for idx, route in enumerate(self.routes):
            route_color = cmap(idx % cmap.N)  # 使用余数操作确保颜色循环重复
            for i in range(len(route) - 1):
                x1, y1 = self.digraph.nodes[route[i]].x, self.digraph.nodes[route[i]].y
                x2, y2 = (
                    self.digraph.nodes[route[i + 1]].x,
                    self.digraph.nodes[route[i + 1]].y,
                )
                plt.arrow(
                    x1,
                    y1,
                    x2 - x1,
                    y2 - y1,
                    head_width=0.5,
                    head_length=0.5,
                    ec=route_color,
                    fc="black",
                )
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title(
            f"CVPR With {self.vehicle_num} Vehicles, {self.vehicle_capacity} Capacity, {self.digraph.customer_num} Customers"
        )
        plt.legend()
        plt.savefig(
            f"CVRP\Figs\CVPR With {self.vehicle_num} Vehicles, {self.vehicle_capacity} Capacity, {self.digraph.customer_num} Customers.png"
        )
        plt.show()

    def print_problem_information(self):
        print("-" * 20, "Problem Information", "-" * 20)
        print(f"节点总数: {self.digraph.node_num}")
        print(f"客户点总数: {self.digraph.customer_num}")
        print(f"车辆总数: {self.vehicle_num}")
        print(f"车容量: {self.vehicle_capacity}")


def main():
    CUSTOMER_NUM = 10
    VEHICLE_NUM = 4
    VEHICLE_CAPACITY = 50
    digraph = initialize_graph(
        data_path="dataset/Solomon/R101.txt", customer_num=CUSTOMER_NUM
    )
    model = CVPR(digraph, vehicle_num=VEHICLE_NUM, vehicle_capacity=VEHICLE_CAPACITY)
    model.violently_solve()
    model.plot_solution()


if __name__ == "__main__":
    main()
