# 用three_index模型求解VRPTW问题
# t_ij = c_ij

from gurobipy import Model, GRB, quicksum
import matplotlib.pyplot as plt
import scienceplots
import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)
from env import *


class VRPTW:
    def __init__(
        self, digraph: DiGraph, vehicle_num: int, vehicle_capacity: float = 100
    ):
        self.digraph = digraph
        self.vehicle_num = vehicle_num
        self.vehicle_capacity = vehicle_capacity
        self.print_problem_information()
        self.depot = self.digraph.depot 
        self.end = self.digraph.end 
        self.node_num = self.digraph.node_num
        self.nodes = self.digraph.nodes
        
        self.o = 0  # depot
        self.d = self.digraph.node_num - 1  # end

        # 总节点列表
        self.V = list(range(self.digraph.node_num))

        # 客户节点列表
        self.N = self.V[1:-1]

        # (i,j) \in A时的集合
        self.I_set = self.V[:-1]
        self.J_set = self.V[1:]

        # 车辆列表
        self.K = list(range(self.vehicle_num))

        self.M = np.zeros((self.digraph.node_num, self.digraph.node_num))
        self.cal_M()

        self.MIP_routes = None

    def violently_solve(self):

        # 创建模型
        MIP_model = Model("VRPTW_three_index")

        # 设置Gap
        MIP_model.setParam("MIPGAP", 0.05)

        # 创建变量

        ## x_ijk,车 k 是否经过 (i,j) \in A
        x = MIP_model.addVars(self.I_set, self.J_set, self.K, vtype=GRB.BINARY, name="x")

        ## T_ik, 车 k 到达 i \in V 的时间
        T = MIP_model.addVars(self.V, self.K, vtype=GRB.CONTINUOUS, name="T")

        # 目标函数
        MIP_model.setObjective(
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
        MIP_model.addConstrs(
            quicksum(x[i, j, k] for k in self.K for j in self.J_set if j != i) == 1
            for i in self.N
        )

        ## o 的出度为1
        MIP_model.addConstrs(
            quicksum(x[self.o, j, k] for j in self.J_set if j != self.o) == 1
            for k in self.K
        )

        ## d 的入度为1
        MIP_model.addConstrs(
            quicksum(x[i, self.d, k] for i in self.I_set if i != self.d) == 1
            for k in self.K
        )

        ## flow constraint
        MIP_model.addConstrs(
            quicksum(x[i, j, k] for i in self.I_set if i != j)
            - quicksum(x[j, i, k] for i in self.J_set if i != j)
            == 0
            for k in self.K
            for j in self.N
        )

        ## time window + MTZ 约束
        MIP_model.addConstrs(
            T[i, k]
            + self.digraph.time(i, j)
            + self.digraph.nodes[i].service_time
            - T[j, k]
            <= (1 - x[i, j, k]) * self.M[i, j]
            for k in self.K
            for i in self.I_set
            for j in self.J_set
        )

        MIP_model.addConstrs(
            T[i, k] >= self.digraph.nodes[i].ready_time for k in self.K for i in self.V
        )
        MIP_model.addConstrs(
            T[i, k] <= self.digraph.nodes[i].due for k in self.K for i in self.V
        )

        ## 容量约束
        MIP_model.addConstrs(
            quicksum(
                self.digraph.nodes[i].demand
                * quicksum(x[i, j, k] for j in self.J_set if j != i)
                for i in self.N
            )
            <= self.vehicle_capacity
            for k in self.K
        )

        # 求解
        MIP_model.optimize()
        self.objective = MIP_model.ObjVal
        self.MIP_routes = self.get_routes(MIP_model)
    
    def column_generation_solve(self):
        
        # 模型里的a_ir
        initial_routes = self.get_initial_routes()
        self.routes_vector = self.convert_routes_to_vectors(initial_routes)
        self.route_costs = np.array([self.cal_route_cost(route) for route in initial_routes])
        self.routes_num = len(initial_routes)
        
        
        
        # 创建模型
        m = Model("VRPTW_column_generation")
        
        # 设置Gap
        m.setParam("MIPGAP", 0.05)
        
        # 创建变量
        
        ## y_r = 1, if route r is selected
        y = m.addVars(self.routes_num, vtype=GRB.CONTINUOUS, name="y")
        
        # 目标函数
        m.setObjective(quicksum(self.route_costs[r] * y[r] for r in range(self.routes_num)),GRB.MINIMIZE)
        
        # 设置约束
        
        ## 每个节点只被访问一次
        visit_constrs = m.addConstrs(quicksum(self.routes_vector[i,r] * y[r] for r in range(self.routes_num)) == 1 for i in self.N)

        ## 选择的routes数小于等于车辆数
        m.addConstr(quicksum(y[r] for r in range(self.routes_num)) <= self.vehicle_num)
        
        # 求解
        m.optimize()

        

        # 第一个约束的对偶变量，就是每个节点的价格，令depot节点的价格为0
        duals = np.array([visit_constrs[i].Pi for i in self.N])
        duals = np.insert(duals, 0, 0)
        
        for r in range(self.routes_num):
            print(y[r].x)
        
    
    # 用贪心法获取初始routes    
    def get_initial_routes(self):
        routes = []
        unserved_customer = self.digraph.nodes[1:-1]
        while unserved_customer:
            current_time = 0
            current_load = 0
            current_node = self.depot 
            route = [current_node.index]
            while True:
                unserved_customer.sort(key=lambda x: x.ready_time)
                to_customer = None
                for customer in unserved_customer:
                    time = max(current_time + self.digraph.time(current_node.index, customer.index)+current_node.service_time, customer.ready_time)
                    load = current_load + customer.demand
                    if time <= customer.due and load <= self.vehicle_capacity:
                        current_time = time 
                        current_load = load
                        to_customer = customer
                        break
                if to_customer is None:
                    route.append(self.end.index)
                    break 
                current_node = to_customer
                route.append(to_customer.index)
                unserved_customer.remove(to_customer)
            routes.append(route)
        return routes 
    
    
     
     # 判断路径是否可行
    def route_is_feasible(self,route:List[int]):
        current_time = 0 
        current_load = 0 
        for i in range(len(route) - 1):
            current_time = max(current_time + self.digraph.time(route[i],route[i+1]) + self.digraph.nodes[route[i]].service_time, self.digraph.nodes[route[i+1]].ready_time)
            current_load += self.digraph.nodes[route[i]].demand
            if current_time > self.digraph.nodes[route[i+1]].due or current_load > self.vehicle_capacity:
                if current_time > self.digraph.nodes[route[i+1]].due:
                    print(f"节点{route[i+1]}超时")
                if current_load > self.vehicle_capacity:
                    print(f"节点{route[i+1]}超载")
                return False
        return True
    
    # 计算MTZ约束中的M
    def cal_M(self):
        self.M = np.zeros((self.digraph.node_num, self.digraph.node_num))
        for i in self.I_set:
            for j in self.J_set:
                self.M[i, j] = max(
                    self.digraph.nodes[i].due
                    + self.digraph.nodes[i].service_time
                    + self.digraph.time(i, j)
                    - self.digraph.nodes[j].ready_time,
                    0,
                )
    
    
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
    
    def cal_route_cost(self,route):
        cost = 0
        for i in range(len(route) - 1):
            cost += self.digraph.cost(route[i],route[i+1])
        return cost
    
    def convert_routes_to_vectors(self,routes):
        vectors = np.zeros((self.node_num,len(routes)))
        for i,route in enumerate(routes):
            for j in route:
                vectors[j,i] = 1
        return vectors

    def plot_solution(self):
        self.plot_routes(self.MIP_routes)
        
    def plot_routes(self,routes):
        plt.style.use(["science"])
        plt.figure(figsize=(10, 10))
        plt.scatter(
            [node.x for node in self.digraph.nodes],
            [node.y for node in self.digraph.nodes],
            c="r",
            label="customer",
            s=100,
        )
        plt.annotate("depot", (self.digraph.nodes[0].x, self.digraph.nodes[0].y))
        for i in range(1, self.digraph.node_num - 1):
            plt.annotate(
                self.digraph.nodes[i].index,
                (self.digraph.nodes[i].x, self.digraph.nodes[i].y),
            )
        # 用箭头画出路径
        cmap = plt.get_cmap("Dark2")  # 使用tab10颜色映射，可选择其他颜色映射
        for idx,route in enumerate(routes):
            route_color = cmap(idx % cmap.N)  # 使用余数操作确保颜色循环重复
            for i in range(len(route) - 1):
                x1, y1 = self.digraph.nodes[route[i]].x, self.digraph.nodes[route[i]].y
                x2, y2 = (
                    self.digraph.nodes[route[i + 1]].x,
                    self.digraph.nodes[route[i + 1]].y,
                )
                plt.arrow(x1, y1, x2 - x1, y2 - y1, head_width=0.5, head_length=0.5,ec=route_color,fc='black')

        plt.xlabel("x")
        plt.ylabel("y")
        plt.title(
            f"VPRTW With {self.vehicle_num} Vehicles, {self.vehicle_capacity} Capacity, {self.digraph.customer_num} Customers"
        )
        plt.legend()
        plt.show()
        

    def print_problem_information(self):
        print("-" * 20, "Problem Information", "-" * 20)
        print(f"节点总数: {self.digraph.node_num}")
        print(f"客户点总数: {self.digraph.customer_num}")
        print(f"车辆总数: {self.vehicle_num}")
        print(f"车容量: {self.vehicle_capacity}")


def main():
    CUSTOMER_NUM = 10
    VEHICLE_NUM = 6
    VEHICLE_CAPACITY = 200
    digraph = initialize_graph(
        data_path="dataset/Solomon/R101.txt", customer_num=CUSTOMER_NUM
    )
    model = VRPTW(digraph, vehicle_num=VEHICLE_NUM, vehicle_capacity=VEHICLE_CAPACITY)
    model.violently_solve()
    model.plot_routes(model.get_initial_routes())
    print(model.get_initial_routes())
    model.column_generation_solve()
    

if __name__ == "__main__":
    main()
