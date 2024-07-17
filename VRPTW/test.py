import gurobipy as gp
from gurobipy import GRB

# 创建一个模型
model = gp.Model()

# 添加变量
x = model.addVar(name="x")
y = model.addVar(name="y")
z = model.addVar(name="z")

# 设置目标函数
model.setObjective(2 * x + 3 * y + z, GRB.MAXIMIZE)

# 添加约束
c0 = model.addConstr(x + y + z <= 4, "c0")
c1 = model.addConstr(x - y + z >= 1, "c1")
c2 = model.addConstr(x + 2 * y + 3 * z == 3, "c2")

# 优化模型
model.optimize()

# 检查模型是否最优求解
if model.status == GRB.OPTIMAL:
    # 输出标准形式的约束矩阵 A
    for constr in model.getConstrs():
        row = [
            constr.getCoefficient(x),
            constr.getCoefficient(y),
            constr.getCoefficient(z),
        ]
        print(row)
else:
    print("Model is not optimal.")
