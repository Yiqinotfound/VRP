# 定义一些VRP问题基本的数据类型

from typing import List
import numpy as np
import pandas as pd
import copy


class Node:
    _index_counter = 0

    def __init__(
        self,
        x: float,
        y: float,
        demand: float,
        ready_time: float,
        due: float,
        service_time: float,
    ):
        self.index = Node._index_counter
        Node._index_counter += 1

        self.x: float = x
        self.y: float = y
        self.demand: float = demand
        self.ready_time: float = ready_time
        self.due: float = due
        self.service_time: float = service_time
        self.labels: List[Label] = []
        self.successor_indexes: List[int] = []
        
    def print(self):
        print("Index:", self.index)
        print("X:", self.x)
        print("Y:", self.y)
        print("Demand:", self.demand)
        print("Ready time:", self.ready_time)
        print("Due:", self.due)
        print("Service time:", self.service_time)
        print('-----------------------')
       


class Label:
    def __init__(
        self,
        path: List[int],
        resources_consumed: List[float],
        unreachable_nodes_num: int,
        unreachable_nodes_vector: List[bool],
        cost: float,
    ):
        self.path: List[int] = path
        self.resources_consumed: List[float] = resources_consumed
        self.unreachable_nodes_num: int = unreachable_nodes_num
        self.unreachable_nodes_vector: List[bool] = unreachable_nodes_vector
        self.cost = cost

    def __eq__(self, other):
        if (
            self.path == other.path
            and self.resources_consumed == other.resources_consumed
            and self.unreachable_nodes_num == other.unreachable_nodes_num
            and self.unreachable_nodes_vector == other.unreachable_nodes_vector
            and self.cost == other.cost
        ):
            return True
        return False

    def dominated_by(self, other):
        if self.__eq__(other):
            return True
        if other.cost > self.cost:
            return False
        for i in range(len(self.resources_consumed)):
            if other.resources_consumed[i] >= self.resources_consumed[i]:
                return False
        if other.unreachable_nodes_num >= self.unreachable_nodes_num:
            return False
        for i in range(len(self.unreachable_nodes_vector)):
            if (
                other.unreachable_nodes_vector[i] == True
                and self.unreachable_nodes_vector[i] == False
            ):
                return False
        return True

    def print(self):
        print("Path:", self.path)
        print("Resources consumed:", self.resources_consumed)
        print("Unreachable nodes num:", self.unreachable_nodes_num)
        print("Unreachable nodes vector:", self.unreachable_nodes_vector)
        print("Cost:", self.cost)

    def __eq__(self, other):
        if (
            self.resources_consumed == other.resources_consumed
            and self.unreachable_nodes_num == other.unreachable_nodes_num
            and self.unreachable_nodes_vector == other.unreachable_nodes_vector
            and self.cost == other.cost
        ):
            return True
        return False


class Route:
    def __init__(self, nodes: List[Node], cost):
        self.nodes = nodes
        self.cost = cost


class DiGraph:
    def __init__(self, customer_nodes: List[Node], depot: Node):
        self.customer_num = len(customer_nodes)
        self.end = copy.deepcopy(depot)
        self.end.index = len(customer_nodes) + 1
        self.nodes: List[Node] = [depot] + customer_nodes + [self.end]
        self.depot = depot
        self.node_num = len(self.nodes)

        # 计算邻接矩阵
        self.cost_matrix = np.zeros((self.node_num, self.node_num))
        self.cal_adj_matrix()

        # 给每个节点添加successor_indexes
        for node in self.nodes:
            successor_indexes = []
            if node.index == self.end.index:
                node.successor_indexes = successor_indexes
                continue
            else:
                for i in range(self.node_num):
                    if i != node.index and i != self.depot.index:
                        successor_indexes.append(i)
            node.successor_indexes = successor_indexes

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
