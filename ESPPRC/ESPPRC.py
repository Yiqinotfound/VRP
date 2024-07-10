import sys
import os
from gurobipy import Model, GRB, quicksum
import matplotlib.pyplot as plt
import scienceplots
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_dir)
from env import *


class ESPPRC:
    def __init__(self, digraph: DiGraph, vehicle_capacity: float):
        self.digraph = digraph
        self.vehicle_capacity = vehicle_capacity
        self.node_num = digraph.node_num
        self.customer_num = digraph.customer_num
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

        self.routes = None

    def violently_solve(self):

        # 创建模型
        m = Model("ESPPRC")

        # 创建变量

        ## x_ij, (i,j) \in A
        x = m.addVars(self.I_set, self.J_set, vtype=GRB.BINARY, name="x")

        ## T_i, 节点i的到达时间
        T = m.addVars(self.V, vtype=GRB.CONTINUOUS, name="T")

        # 目标函数
        m.setObjective(
            quicksum(
                self.digraph.cost(i, j) * x[i, j]
                for i in self.I_set
                for j in self.J_set
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
            quicksum(x[i,j] for i in self.I_set if i != j)
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

        m.addConstrs(T[i] >= self.digraph.nodes[i].ready_time for i in self.V)
        m.addConstrs(T[i] <= self.digraph.nodes[i].due for i in self.V)

        ## 容量约束
        m.addConstr(
            quicksum(
                self.digraph.nodes[i].demand
                * quicksum(x[i, j] for j in self.J_set if j != i)
                for i in self.N
            )
            <= self.vehicle_capacity
        )
        m.optimize()
        self.objective = m.objVal
        self.routes = self.get_routes(m)
    
    def dp_solve(self):
        pass

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

    def get_routes(self, m):
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
        return [route]

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
        plt.title(
            f"ESPPRC With {self.vehicle_capacity} Capacity, {self.digraph.customer_num} Customers"
        )
        plt.legend()
        plt.show()


def main():
    CUSTOMER_NUM = 10
    VEHICLE_CAPACITY = 100
    digraph = initialize_graph("dataset/Solomon/R101.txt", CUSTOMER_NUM)
    model = ESPPRC(digraph=digraph, vehicle_capacity=VEHICLE_CAPACITY)
    model.violently_solve()
    model.plot_solution()
    print(model.routes)

if __name__ == "__main__":
    main()
