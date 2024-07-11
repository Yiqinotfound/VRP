import sys
import time
import os
from gurobipy import Model, GRB, quicksum
import matplotlib.pyplot as plt
import scienceplots

parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_dir)
from env import *


class ESPPRC:
    def __init__(self, digraph: DiGraph, alpha: np.ndarray, vehicle_capacity: float):
        self.digraph = digraph
        self.alpha = alpha
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

        self.MTZ_routes = None
        
    
        

    def violently_solve(self):

        # 创建模型
        m = Model("ESPPRC")

        # 设置GAP 
        m.setParam("MIPGap", 0.01)

        # 创建变量

        ## x_ij, (i,j) \in A
        x = m.addVars(self.I_set, self.J_set, vtype=GRB.BINARY, name="x")

        ## T_i, 节点i的到达时间
        T = m.addVars(self.V, vtype=GRB.CONTINUOUS, name="T")

        # 目标函数
        m.setObjective(
            quicksum(
                (self.digraph.cost(i, j) - self.alpha[i]) * x[i, j]
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
        self.objective = m.objVal
        self.MTZ_routes = self.get_routes_MTZ_model(m)
        self.m = m

    def dp_solve(self):
        def expand(label: Label, from_node: Node, to_node: Node):
            """
            Function that returns the label resulting from the extension of 'label' to 'node' when the extension is possible, nothing otherwise.
            """
            # 更新resource的消耗，检查是否满足
            # [0]是时间，[1]是容量
            resource_consumed = [
                label.resources_consumed[0]
                + self.digraph.cost(from_node.index, to_node.index),
                label.resources_consumed[1] + to_node.demand,
            ]
            if (
                resource_consumed[0] <= to_node.due
                and resource_consumed[1] <= self.vehicle_capacity
            ):
                if resource_consumed[0] < to_node.ready_time:
                    resource_consumed[0] = to_node.ready_time
                    
                # 更新unreachable_nodes_vector和unreachable_nodes_num
                unreachable_nodes_vector = copy.deepcopy(label.unreachable_nodes_vector)
                unreachable_nodes_vector[to_node.index] = True
                unreachable_nodes_num = label.unreachable_nodes_num + 1
                if to_node.index == self.digraph.end.index:
                    unreachable_nodes_num = self.digraph.node_num
                    unreachable_nodes_vector = [
                        True for _ in range(self.digraph.node_num)
                    ]
                cost = (
                    label.cost
                    + self.digraph.cost(from_node.index, to_node.index)
                    - self.alpha[from_node.index]
                )
                new_label = Label(
                    path=label.path + [to_node.index],
                    resources_consumed=resource_consumed,
                    unreachable_nodes_num=unreachable_nodes_num,
                    unreachable_nodes_vector=unreachable_nodes_vector,
                    cost=cost,
                )
                return new_label
            return None

        def EEF(extended_labels, original_labels):
            total_labels: List[Label] = original_labels + extended_labels
            # print('original_labels:', [label.cost for label in original_labels])
            # print('extended_labels:', [label.cost for label in extended_labels])
            if len(total_labels) == 0:
                return total_labels
            label_to_remove = []

            # 遍历所有label,如果有被支配的，就删除
            for i in range(len(total_labels)):
                for j in range(len(total_labels)):
                    if i != j:
                        if total_labels[i].dominated_by(total_labels[j]) and total_labels[i] not in label_to_remove:
                            label_to_remove.append(total_labels[i])
            for label in label_to_remove:
                total_labels.remove(label)

            return total_labels

        E: List[Node] = []
        F = [
            [[] for _ in range(self.digraph.node_num)]
            for _ in range(self.digraph.node_num)
        ]

        self.digraph.depot.labels.append(
            Label(
                path=[0],
                resources_consumed=[0, 0],
                unreachable_nodes_num=1,
                unreachable_nodes_vector=[True]
                + [False for _ in range(self.digraph.node_num - 1)],
                cost=0,
            )
        )
        for node in self.digraph.nodes:
            if node.index != self.digraph.depot.index:
                node.labels = []
        E.append(self.digraph.depot)

        while True > 0:
            # print('要处理的节点有:', [node.index for node in E])
            # print('1,2,3节点的labels数量:', [len(node.labels) for node in self.digraph.nodes[1:]])
            if len(E) == 0:
                break
            node = E.pop(0)
            for successor_index in node.successor_indexes:
                F:List[Label] = []
                for label in node.labels:
                    if label.unreachable_nodes_vector[successor_index] == False:
                        new_label = expand(
                            label, node, self.digraph.nodes[successor_index]
                        )
                        if new_label is not None:
                            F.append(new_label)
                original_labels = copy.deepcopy(
                    self.digraph.nodes[successor_index].labels
                )
                self.digraph.nodes[successor_index].labels = EEF(
                    F,
                    self.digraph.nodes[successor_index].labels,
                )
                

                # 判断是否有变化
                if (
                    original_labels != self.digraph.nodes[successor_index].labels
                    and successor_index != self.digraph.end.index
                ):
                    E.append(self.digraph.nodes[successor_index])
            # print('处理之后1,2,3节点的labels数量:', [len(node.labels) for node in self.digraph.nodes[1:]])
            # print('剩下处理节点为:', [node.index for node in E])
            # print('-----------------------------------')

        print(min([label.cost for label in self.digraph.end.labels]))
        # 打印出最优路径
        routes = []
        for label in self.digraph.end.labels:
            if label.cost == min([label.cost for label in self.digraph.end.labels]):
                routes.append(label.path)
        self.DP_routes = routes 
            

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

    def get_routes_MTZ_model(self, m):
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

    def cal_route_cost(self, route):
        distance = 0
        for i in range(len(route) - 1):
            distance += self.digraph.cost(route[i], route[i + 1]) - self.alpha[route[i]]
        return distance
    
    # 判断route是否可行
    def route_is_feasible(self,route):
        current_time = 0 # 从 depot出发时间为0
        current_load = 0 # 从depot出发load 为 0
        for i in range(len(route) - 1):

            # 更新current_time
            current_time += self.digraph.cost(route[i], route[i + 1])
            
            # 更新load
            current_load += self.digraph.nodes[route[i]].demand
            if current_load > self.vehicle_capacity:
                return False
            if current_time > self.digraph.nodes[route[i+1]].due:
                return False
            elif current_time < self.digraph.nodes[route[i+1]].ready_time:
                current_time = self.digraph.nodes[route[i+1]].ready_time
            
            # 更新current_time
            current_time += self.digraph.nodes[route[i]].service_time
        return True


def main():
    CUSTOMER_NUM = 15
    VEHICLE_CAPACITY = 200
    digraph = initialize_graph("dataset/Solomon/R101.txt", CUSTOMER_NUM)

    alpha = np.zeros(digraph.node_num)
    alpha[:] = 20
    alpha[0] = 0

    model = ESPPRC(digraph=digraph, alpha=alpha, vehicle_capacity=VEHICLE_CAPACITY)

    start = time.time()
    model.dp_solve()
    end = time.time()
    print("DP solve time:", end - start)
    model.plot_solution(model.DP_routes)
    


if __name__ == "__main__":
    main()
