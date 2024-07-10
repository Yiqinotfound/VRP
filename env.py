# 定义一些VRP问题基本的数据类型

from typing import List 
import numpy as np 
import pandas as pd 

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

class Route:
    def __init__(self,nodes:List[Node],cost):
        self.nodes = nodes
        self.cost = cost

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
