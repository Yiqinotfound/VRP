# 用three_index模型求解VRPTW问题
# t_ij = c_ij
import sys
import os
import time

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)
from gurobipy import Model, GRB, quicksum
import matplotlib.pyplot as plt
import scienceplots
from ESPPRC.ESPPRC import ESPPRC
from env import *
import copy


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

        self.mip_paths = None

    def violently_solve(self):

        # 创建模型
        MIP_model = Model("VRPTW_three_index")

        # 设置Gap
        MIP_model.setParam("MIPGAP", 0)

        # 创建变量

        ## x_ijk,车 k 是否经过 (i,j) \in A
        x = MIP_model.addVars(
            self.I_set, self.J_set, self.K, vtype=GRB.BINARY, name="x"
        )

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
        self.mip_objective = MIP_model.ObjVal
        self.mip_paths = self.get_paths(MIP_model)

    def column_generation_solve(self):

        # 获取初始routes
        initial_paths = self.get_naive_paths() 

        route_costs = np.array([self.cal_path_cost(path) for path in initial_paths])
        arc_matrixes = [self.get_arc_matrix_from_path(path) for path in initial_paths]
        routes = [
            Route(initial_paths[i], arc_matrixes[i], route_costs[i])
            for i in range(len(initial_paths))
        ]
        self.initial_solution_node = SolutionNode(
            digraph=self.digraph,
            vehicle_capcacity=self.vehicle_capacity,
            vehicle_num=self.vehicle_num,
            routes=routes,
        )
        self.initial_solution_node.solve()

    # 用贪心法获取初始paths
    def get_greedy_paths(self):
        paths = []
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
                    time = max(
                        current_time
                        + self.digraph.time(current_node.index, customer.index)
                        + current_node.service_time,
                        customer.ready_time,
                    )
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
            paths.append(route)
        return paths

    def get_naive_paths(self):
        paths = []
        for i in range(1, self.digraph.node_num - 1):
            paths.append([0, i, self.digraph.node_num - 1])
        return paths

    # 判断path是否可行
    def path_is_feasible(self, route: List[int]):
        current_time = 0
        current_load = 0
        for i in range(len(route) - 1):
            current_time = max(
                current_time
                + self.digraph.time(route[i], route[i + 1])
                + self.digraph.nodes[route[i]].service_time,
                self.digraph.nodes[route[i + 1]].ready_time,
            )
            current_load += self.digraph.nodes[route[i]].demand
            if (
                current_time > self.digraph.nodes[route[i + 1]].due
                or current_load > self.vehicle_capacity
            ):
                if current_time > self.digraph.nodes[route[i + 1]].due:
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

    def get_paths(self, m):
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

    def cal_path_cost(self, path):
        cost = 0
        for i in range(len(path) - 1):
            cost += self.digraph.cost(path[i], path[i + 1])
        return cost

    def get_arc_matrix_from_path(self, path):
        arc_matrix = np.zeros((self.node_num, self.node_num))
        for i in range(len(path) - 1):
            arc_matrix[path[i]][path[i + 1]] = 1
        return arc_matrix

    def plot_solution(self):
        self.plot_paths(self.mip_paths)

    def plot_paths(self, routes):
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
        cmap = plt.get_cmap("Dark2")
        for idx, route in enumerate(routes):
            route_color = cmap(idx % cmap.N)
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


class SolutionNode:
    _id_counter = 0

    def __init__(
        self,
        digraph: DiGraph,
        vehicle_num: int,
        vehicle_capcacity: float,
        routes: List[Route],
        forbidden_arcs: List[tuple] = None,
        retained_arcs: List[tuple] = None,
    ):

        self.digraph = digraph
        self.node_num = self.digraph.node_num

        # 禁止的arc [(i,j),...]
        self.forbidden_arcs: List[Tuple] = forbidden_arcs
        if self.forbidden_arcs is None:
            self.forbidden_arcs = []

        # 保留的arc [(i,j),...]
        self.retained_arcs: List[tuple] = retained_arcs
        if self.retained_arcs is None:
            self.retained_arcs = []

        # 路径列表，已经除去了分支要求的路径
        self.routes: List[Route] = routes

        # 车辆数和容量
        self.vehicle_num = vehicle_num
        self.vehicle_capacity = vehicle_capcacity

        self.V = list(range(self.node_num))
        self.N = self.V[1:-1]

        self.generated = 0

        # 是否求解RMP达到了最优解
        self.optimal = False

        # 最优解
        self.objective = None

        # RMP是否可行
        self.feasible = True

        # 流量矩阵flow_matrix[i,j]表示(i,j)的流
        self.flow_matrix = np.zeros((self.node_num, self.node_num))

        self.solution_routes = None

        self.forbid_arc_node = None
        self.retain_arc_node = None

        self.id = SolutionNode._id_counter
        SolutionNode._id_counter += 1
        print("创建节点", self.id)
        self.print_routes()

    def solve(self):

        # a_ir = 1 if route r contains node i
        self.visits_vector = cal_routes_visits_vector(self.node_num, self.routes)

        # 构建Master Problem
        self.generate_MP()
        self.master_problem.optimize()

        # 判断是否可行
        if self.master_problem.Status != GRB.OPTIMAL:
            self.feasible = False
            return
        print('-'*40)
        print(self.id, "节点RMP求解成功", self.master_problem.objVal)
        self.y = self.master_problem.getAttr("x", self.master_problem.getVars())
        print('y=',self.y)
        print('-'*40)

        A = self.master_problem.getA().todense()
        B = np.empty((A.shape[0],0), dtype=A.dtype)
        column_num = self.node_num - 2
        for i in range(len(self.y)):
            if self.y[i] > 0:
                column = A[:,i]
                B = np.column_stack((B,column))
                column_num -= 1
                if column_num == 0:
                    break 
        if column_num:
            for i in range(len(self.y)):
                if self.y[i] == 0:
                    column = A[:,i]
                    B = np.column_stack((B,column))
                    column_num -= 1
                    if column_num == 0:
                        break 

        # visit_constraints的对偶变量
        self.duals = self.master_problem.getAttr("pi", self.visit_constraints)
        self.duals = [self.duals[i] for i in self.N]
        self.duals.insert(0, 0)

        # 构建SubProblem
        self.generate_SP()
        self.sub_problem.dp_solve()

        # 最小的reduced cost
        min_rc = self.sub_problem.dp_objective
        if min_rc >= 0:  # reduced_cost >= 0

            # y如果都是整数，说明已经是最优解
            if all(y == 0 or y == 1 for y in self.y):
                self.optimal = True
                self.objective = self.master_problem.objVal
                self.solution_routes = [
                    self.routes[i] for i in range(self.route_num()) if self.y[i] == 1
                ]
                return

            # 计算流量矩阵
            self.flow_matrix = sum(
                self.y[i] * self.routes[i].arc_matrix for i in range(self.route_num())
            )

            # 找到一个非0-1流量的arc
            arc = self.get_fractional_arc(self.flow_matrix)

            new_forbidden_arcs = self.forbidden_arcs.copy()
            new_forbidden_arcs.append(arc)
            new_retained_arcs = self.retained_arcs.copy()
            new_retained_arcs.append(arc)

            self.forbid_arc_node = SolutionNode(
                digraph=self.digraph,
                vehicle_num=self.vehicle_num,
                vehicle_capcacity=self.vehicle_capacity,
                routes=self.filter_routes(routes=self.routes, forbidden_arc=arc),
                forbidden_arcs=new_forbidden_arcs,
                retained_arcs=self.retained_arcs,
            )
            self.retain_arc_node = SolutionNode(
                digraph=self.digraph,
                vehicle_num=self.vehicle_num,
                vehicle_capcacity=self.vehicle_capacity,
                routes=self.filter_routes(routes=self.routes, retained_arc=arc),
                forbidden_arcs=self.forbidden_arcs,
                retained_arcs=new_retained_arcs,
            )
            self.forbid_arc_node.solve()
            self.retain_arc_node.solve()

            self.retain_objective = (
                self.retain_arc_node.objective
                if self.retain_arc_node.objective is not None
                else float("inf")
            )
            self.forbid_objective = (
                self.forbid_arc_node.objective
                if self.forbid_arc_node.objective is not None
                else float("inf")
            )

            self.objective = min(self.retain_objective, self.forbid_objective)
        else:
            print('min_rc=',min_rc)
            new_path = self.sub_problem.dp_path
            new_path_cost = self.cal_path_cost(new_path)
            new_path_arc_matrix = self.get_arc_matrix_from_path(new_path)
            new_route = Route(new_path, new_path_arc_matrix, new_path_cost)
            self.routes.append(new_route)
            print('节点',self.id,'增加了',new_path,'条路径')
            self.print_routes()
            self.generated += 1
            self.solve()

    # 根据routes、forbidden_arcs、retained_arcs生成MP
    def generate_MP(self):

        # 创建模型
        self.master_problem = Model("MP")
        self.master_problem.setParam("OutputFlag", 0)

        # 创建变量
        y = self.master_problem.addVars(
            self.route_num(), vtype=GRB.CONTINUOUS, name="y"
        )

        # 目标函数
        self.master_problem.setObjective(
            quicksum(self.routes[r].cost * y[r] for r in range(self.route_num())),
            GRB.MINIMIZE,
        )

        # 访问约束
        self.visit_constraints = self.master_problem.addConstrs(
            quicksum(self.visits_vector[i, r] * y[r] for r in range(len(self.routes)))
            == 1
            for i in self.N
        )

        if len(self.routes) == 7:
            self.master_problem.write('problem7.lp')
        if len(self.routes) == 6:
            self.master_problem.write('problem6.lp')

    def generate_SP(self):
        self.sub_problem = ESPPRC(
            self.digraph,
            self.vehicle_capacity,
            self.duals,
            self.forbidden_arcs,
            self.retained_arcs,
            self.routes,
        )

    def cal_path_cost(self, path):
        cost = 0
        for i in range(len(path) - 1):
            cost += self.digraph.cost(path[i], path[i + 1])
        return cost

    def get_arc_matrix_from_path(self, path):
        arc_matrix = np.zeros((self.node_num, self.node_num))
        for i in range(len(path) - 1):
            arc_matrix[path[i]][path[i + 1]] = 1
        return arc_matrix

    def route_num(self):
        return len(self.routes)

    def get_fractional_arc(self, flow_matrix: np.ndarray):
        self.fractional_arcs: List[Tuple] = []
        arc = None
        for i in range(1, len(flow_matrix)):
            for j in range(len(flow_matrix) - 1):
                if flow_matrix[i, j] != 0 and self.flow_matrix[i, j] != 1:
                    arc = (i, j)
                    break
            if arc is not None:
                break
        return arc

    def filter_routes(
        self,
        routes: List[Route],
        forbidden_arc: Tuple = None,
        retained_arc: Tuple = None,
    ):
        new_routes: List[Route] = []
        for route in routes:
            flag = True
            for i in range(len(route.path) - 1):
                arc = (route.path[i], route.path[i + 1])
                if forbidden_arc is not None and arc == forbidden_arc:
                    flag = False
                    break
                if retained_arc is not None:
                    if retained_arc[0] != 0:
                        if (
                            arc[0] == retained_arc[0] and arc[1] != retained_arc[1]
                        ) or (arc[1] == retained_arc[1] and arc[0] != retained_arc[0]):
                            flag = False
                            break
                    else:
                        if arc[1] == retained_arc[1] and arc[0] != 0:
                            flag = False
                            break
            if flag:
                new_routes.append(route)
        return new_routes

    def print_routes(self):
        print("共有", self.route_num(), "条路径")
        for route in self.routes:
            print(route)


def main():
    CUSTOMER_NUM = 5
    VEHICLE_NUM = 50
    VEHICLE_CAPACITY = 200
    digraph = initialize_graph(
        data_path="dataset/Solomon/R101.txt", customer_num=CUSTOMER_NUM
    )
    model = VRPTW(
        digraph=digraph, vehicle_num=VEHICLE_NUM, vehicle_capacity=VEHICLE_CAPACITY
    )
    start = time.time()
    model.violently_solve()
    end = time.time()
    mip_time = end - start
    start = time.time()
    model.column_generation_solve()
    end = time.time()
    print("MIP求解时间:", mip_time)
    print("CG求解时间:", end - start)

    print(model.initial_solution_node.objective,model.mip_objective)


if __name__ == "__main__":
    main()
