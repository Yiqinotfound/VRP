# 定义一些VRP问题基本的数据类型

from typing import List, Set
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

    def __eq__(self, other):
        return self.index == other.index

    def __hash__(self):
        return hash(self.index)
    
    def __repr__(self):
        return f"Node(index:{self.index}, x:{self.x}, y:{self.y}, demand:{self.demand}, ready time:{self.ready_time}, due time:{self.due}, service time:{self.service_time})"


class Label:
    def __init__(
        self,
        node: Node,
        path: List[int],
        load: float,
        time: float,
        cost: float,
    ):
        self.node: Node = node
        self.path: List[int] = path
        self.load: float = load
        self.time: float = time
        self.unreachable_nodes: Set[Node] = {self.node}
        self.cost: float = cost
        self.dominated = False

    def __eq__(self, other):
        if (
            self.node == other.node
            and self.path == other.path
            and self.load == other.load
            and self.time == other.time
            and self.unreachable_nodes == other.unreachable_nodes
            and self.cost == other.cost
        ):
            return True
        else:
            return False

    def dominates(self, other):
        return (
            self.cost <= other.cost
            and self.load <= other.load
            and self.time <= other.time
            and self.unreachable_nodes.issubset(other.unreachable_nodes)
        )

    def is_dominated(self):
        """Returns True if this label is dominated by any of the labels 
           associated with the same customer, False otherwise"""

        for label in self.node.labels:
            if label.dominates(self):
                return True 
        return False
    
    def filter_dominated(self):
        """Removes labels dominated by this label on its customer.
        """
        labels = []
        for label in self.node.labels:
            if self.dominates(label):
                label.dominated = True
            else:
                labels.append(label)
        self.node.labels = labels

    def unreachable_nodes_num(self):
        return len(self.unreachable_nodes)

    def __repr__(self):
        return f"Label(customer:{self.node.index}, cost:{self.cost}, load:{self.load}, time:{self.time})"


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
        self.adj_matrix = np.zeros((self.node_num, self.node_num))
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
                self.adj_matrix[i][j] = np.sqrt(
                    (self.nodes[i].x - self.nodes[j].x) ** 2
                    + (self.nodes[i].y - self.nodes[j].y) ** 2
                )

    def distance(self, i: int, j: int):
        return self.adj_matrix[i][j]
    
    def time(self,i:int,j:int):
        # 从i到j的时间暂时等于距离，后续可以加入速度等信息
        return self.adj_matrix[i][j] 
    
    def cost(self,i:int,j:int):
        # 从i到j的成本暂时等于距离，后续可以加入速度等信息
        return self.adj_matrix[i][j]


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
