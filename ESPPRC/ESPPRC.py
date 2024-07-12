import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import deque
from gurobipy import GRB, Model, quicksum
import gurobipy
import sys
import os
import time
import scienceplots

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)
from env import *


class ESPPRC:
    def __init__(self, digraph: DiGraph, vehicle_capacity: float, duals: List[float]):
        self.digraph = digraph
        self.nodes: Set[Node] = set(self.digraph.nodes)
        self.vehicle_capacity = vehicle_capacity
        self.duals = duals
        self.depot = self.digraph.depot
        self.end = self.digraph.end
        self.customer_num = self.digraph.customer_num
        self.node_num = self.digraph.node_num

        self.o = 0  # depot
        self.d = self.node_num - 1  # end

        # 总节点列表
        self.V = list(range(self.digraph.node_num))

        # 客户节点列表
        self.N = self.V[1:-1]

        # (i,j) \in A时的集合
        self.I_set = self.V[:-1]
        self.J_set = self.V[1:]

        self.M = np.zeros((self.digraph.node_num, self.digraph.node_num))
        self.cal_M()

    def feasible_labels_from(self, from_label: Label):
        to_labels = []
        for to_node in self.nodes - from_label.unreachable_nodes - {from_label.node}:
            to_label = self.extended_label(from_label, to_node)
            if not to_label:
                from_label.unreachable_nodes.add(to_node)
            else:
                to_labels.append(to_label)
        return to_labels

    def extended_label(self, from_label: Label, to_node: Node):
        load = from_label.load + to_node.demand
        if load > self.vehicle_capacity:
            return

        from_node = from_label.node
        time = max(
            from_label.time
            + from_node.service_time
            + self.digraph.cost(from_node.index, to_node.index),
            to_node.ready_time,
        )
        if time > to_node.due:
            return

        cost = (
            from_label.cost
            + self.digraph.cost(from_node.index, to_node.index)
            - self.duals[from_node.index]
        )
        path = from_label.path + [to_node.index]
        return Label(to_node, path, load, time, cost)

    def violently_solve(self):

        # 创建模型
        m = Model("ESPPRC")

        # 设置GAP
        m.setParam("MIPGap", 0.01)
        
        # 不要输出log
        m.setParam("OutputFlag", 0)

        # 创建变量

        ## x_ij, (i,j) \in A
        x = m.addVars(self.I_set, self.J_set, vtype=GRB.BINARY, name="x")

        ## T_i, 节点i的到达时间
        T = m.addVars(self.V, vtype=GRB.CONTINUOUS, name="T")

        # 目标函数
        m.setObjective(
            quicksum(
                (self.digraph.cost(i, j) - self.duals[i]) * x[i, j]
                for i in self.I_set
                for j in self.J_set
                if i != j
            ),
            GRB.MINIMIZE,
        )

        # 约束条件

        ## o 的出度为1
        m.addConstr(quicksum(x[self.o, j] for j in self.J_set) == 1)

        ## d 的入度为1
        m.addConstr(quicksum(x[i, self.d] for i in self.I_set) == 1)

        ## flow conservation
        m.addConstrs(
            quicksum(x[i, j] for i in self.I_set if i != j)
            - quicksum(x[j, i] for i in self.J_set if i != j)
            == 0
            for j in self.N
        )

        ## time window + MTZ约束
        m.addConstrs(
            T[i] + self.digraph.cost(i, j) + self.digraph.nodes[i].service_time - T[j]
            <= (1 - x[i, j]) * self.M[i, j]
            for i in self.I_set
            for j in self.J_set
        )

        m

        m.addConstrs(T[i] >= self.digraph.nodes[i].ready_time for i in self.V)
        m.addConstrs(T[i] <= self.digraph.nodes[i].due for i in self.V)

        ## 容量约束
        m.addConstr(
            quicksum(
                self.digraph.nodes[i].demand
                * quicksum(x[i, j] for j in self.J_set if j != i)
                for i in self.I_set
            )
            <= self.vehicle_capacity
        )
        m.optimize()
        self.MIP_objective = m.objVal
        self.MIP_routes = self.get_routes_MIP_model(m)
        self.m = m

    # 计算MTZ约束中的M
    def cal_M(self):
        self.M = np.zeros((self.digraph.node_num, self.digraph.node_num))
        for i in self.I_set:
            for j in self.J_set:
                self.M[i, j] = max(
                    self.digraph.nodes[i].due
                    + self.digraph.nodes[i].service_time
                    + self.digraph.cost(i, j)
                    - self.digraph.nodes[j].ready_time,
                    0,
                )
    
    def get_routes_MIP_model(self, m):
        route = []
        i = 0
        route.append(i)
        while True:
            if i == self.digraph.node_num - 1:
                break
            for j in range(1, self.digraph.node_num):
                if j != i and m.getVarByName(f"x[{i},{j}]").x > 0.5:
                    route.append(j)
                    i = j
                    break
        return route

    def dp_solve(self):
        for node in self.digraph.nodes:
            node.labels = []
        depot_label = Label(self.depot, path=[self.depot.index], load=0, time=0, cost=0)
        to_be_extended = deque([depot_label])
        while to_be_extended:
            from_label = to_be_extended.popleft()
            if from_label.dominated:
                continue
            to_labels: List[Label] = self.feasible_labels_from(from_label)
            for to_label in to_labels:
                to_node = to_label.node

                # 到达end节点时，不再扩展
                if to_node is not self.end:
                    to_label.unreachable_nodes.update(from_label.unreachable_nodes)
                    if to_label.is_dominated():
                        continue
                    to_label.filter_dominated()
                    to_be_extended.append(to_label)
                to_node.labels.append(to_label)
        self.end.labels.sort(key=lambda x: x.cost)
        self.dp_path = self.end.labels[0].path
        self.DP_objective = self.end.labels[0].cost

    def cal_route_cost(self, route: List[int]):
        cost = 0
        for i in range(len(route) - 1):
            cost += self.digraph.cost(route[i], route[i + 1]) - self.duals[route[i]]
        return cost

    def plot_solution(self,routes):
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
        for idx, route in enumerate(routes):
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
            f"ESPPRC With {self.vehicle_capacity} Capacity, {self.digraph.customer_num} Customers"
        )
        plt.legend()
        plt.show()


def main():
    VEHICLE_CAPACITY = 100
    CUSTOMER_NUM = 10
    digraph = initialize_graph(
        data_path="dataset/Solomon/R101.txt", customer_num=CUSTOMER_NUM
    )

    duals = [50 for _ in range(digraph.node_num)]
    model = ESPPRC(digraph=digraph, vehicle_capacity=VEHICLE_CAPACITY, duals=duals)

    start = time.time()
    model.violently_solve()
    end = time.time()
    print("MIP Path:",model.MIP_routes)
    print("MIP Cost:",model.MIP_objective)
    print("MIP Time:",end-start)

    start = time.time()
    model.dp_solve()
    end = time.time()
    print("DP Path:",model.dp_path)
    print("DP Cost",model.DP_objective)
    print("DP Time:",end-start)
    model.plot_solution([model.dp_path])

if __name__ == "__main__":
    main()
