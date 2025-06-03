import numpy as np
from math import sqrt

def parse_fe_output(file_path):
    """
    解析 fe.out 文件，提取节点坐标、节点位移和单元连接信息
    """
    node_coords = {}
    node_displacements = {}
    elements = []
    current_section = None
    material_params = {}
    
    with open(file_path, 'r') as file:
        for line in file:
            # 识别当前解析的部分
            if "N O D A L   P O I N T   D A T A" in line:
                current_section = "node_data"
                # 跳过后续两行标题
                next(file)
                next(file)
                continue
            elif "E L E M E N T   I N F O R M A T I O N" in line:
                current_section = "element_info"
                # 找到单元信息表格开始位置
                while "ELEMENT     NODE     NODE     NODE" not in next(file):
                    pass
                next(file)  # 跳过标题行
                continue
            elif "D I S P L A C E M E N T S" in line:
                current_section = "displacements"
                # 跳过后续两行标题
                next(file)
                next(file)
                continue
            elif "YOUNG'S     POISSON'S     THICKNESS" in line:
                current_section = "material"
                next(file)  # 跳过标题行
                continue
            elif "S T R E S S" in line or "S O L U T I O N" in line:
                current_section = None
            
            # 空行或注释行跳过
            if not line.strip():
                current_section = None
                continue
            
            # 处理不同部分的数据
            if current_section == "node_data":
                # 处理节点坐标数据
                parts = line.split()
                if len(parts) >= 6:
                    node_id = int(parts[0])
                    # 读取坐标 (x, y, z)
                    x = float(parts[5].replace('D', 'e'))
                    y = float(parts[6].replace('D', 'e'))
                    z = float(parts[7].replace('D', 'e'))
                    node_coords[node_id] = (x, y)
            
            elif current_section == "element_info":
                # 处理单元连接信息
                parts = line.split()
                if len(parts) >= 5:
                    # 忽略第一个元素（单元ID），读取三个节点ID
                    nodes = [int(parts[1]), int(parts[2]), int(parts[3])]
                    elements.append(nodes)
            
            elif current_section == "displacements":
                # 处理节点位移数据
                parts = line.split()
                if len(parts) >= 4:
                    node_id = int(parts[0])
                    ux = float(parts[1].replace('D', 'e'))
                    uy = float(parts[2].replace('D', 'e'))
                    uz = float(parts[3].replace('D', 'e'))
                    node_displacements[node_id] = (ux, uy)
            
            elif current_section == "material":
                # 读取材料参数
                parts = line.split()
                if len(parts) >= 3:
                    material_params["E"] = float(parts[1].replace('D', 'e'))
                    material_params["nu"] = float(parts[2].replace('D', 'e'))
    
    return node_coords, node_displacements, elements, material_params

# 解析解函数 - 需要材料参数
def u_exact(x, y, E, I):
    F, h = 1.0, 1.0
    return (F * h) / (E * I) * x * y

def v_exact(x, E, I):
    F, h = 1.0, 1.0
    return -(F * h) / (2 * E * I) * x**2

# 三角形形函数
def shape_functions(xi, eta):
    return np.array([xi, eta, 1 - xi - eta])

# 雅可比行列式计算（返回行列式值）
def jacobian_det(nodes):
    # 节点坐标矩阵
    nodes = np.array(nodes)
    # 参考单元到物理单元的映射雅可比矩阵
    J = np.zeros((2, 2))
    for i in range(2):
        J[0, i] = nodes[0, 0] - nodes[2, 0] if i == 0 else nodes[1, 0] - nodes[2, 0]
        J[1, i] = nodes[0, 1] - nodes[2, 1] if i == 0 else nodes[1, 1] - nodes[2, 1]
    return abs(np.linalg.det(J))

def calculate_l2_error(file_path):
    # 解析输出文件
    node_coords, node_displacements, elements, material_params = parse_fe_output(file_path)
    
    # 获取材料参数
    E = material_params["E"] if "E" in material_params else 1e4
    nu = material_params["nu"] if "nu" in material_params else 0.0
    
    # 计算截面惯性矩
    b, h = 1.0, 1.0  # 根据问题描述
    I = b * h**3 / 12
    
    # 6点高斯积分点（精确积分四次多项式）
    gauss_points = [
        # (xi, eta, weight)
        (0.091576213509770, 0.091576213509770, 0.054975871827661),
        (0.091576213509770, 0.816847572980459, 0.054975871827661),
        (0.816847572980459, 0.091576213509770, 0.054975871827661),
        (0.108103018168070, 0.445948490915965, 0.111690794839005),
        (0.445948490915965, 0.108103018168070, 0.111690794839005),
        (0.445948490915965, 0.445948490915965, 0.111690794839005)
    ]
    
    # 初始化误差平方和
    error_squared = 0.0
    
    # 遍历每个单元
    for elem in elements:
        # 获取单元节点坐标和位移
        node_ids = elem
        nodes = [node_coords[n] for n in node_ids]
        disps = [node_displacements[n] for n in node_ids]
        
        # 计算单元雅可比行列式（用于面积计算）
        jac_det = jacobian_det(nodes)
        area = jac_det / 2.0
        
        # 遍历高斯积分点
        for xi, eta, weight in gauss_points:
            # 计算形函数
            N = shape_functions(xi, eta)
            
            # 计算全局坐标 (x, y)
            x = sum(N[i] * nodes[i][0] for i in range(3))
            y = sum(N[i] * nodes[i][1] for i in range(3))
            
            # 计算有限元解
            u_fem = sum(N[i] * disps[i][0] for i in range(3))
            v_fem = sum(N[i] * disps[i][1] for i in range(3))
            
            # 计算解析解
            u_true = u_exact(x, y, E, I)
            v_true = v_exact(x, E, I)
            
            # 计算误差
            error_u = u_fem - u_true
            error_v = v_fem - v_true
            
            # 累积误差平方和
            error_squared += weight * jac_det * (error_u**2 + error_v**2)
    
    # 计算L₂范数误差
    l2_error = sqrt(error_squared)
    
    # 输出结果
    print(f"材料参数: E = {E:.1e}, nu = {nu:.1f}")
    print(f"网格尺寸: h = {h:.1f}")
    print(f"加载力: F = {1.0:.1f}")
    print(f"截面惯性矩: I = {I:.5e}")
    print(f"\n位移L₂范数误差: {l2_error:.5e}")
    return l2_error

# 主程序
if __name__ == "__main__":
    fe_file = "/root/homework/finite_element/STAPpp/data/convergence/dats/1.out"  # 输入文件路径
    try:
        error = calculate_l2_error(fe_file)
        print(f"误差值: {error:.5e} (科学计数法)")
    except Exception as e:
        print(f"错误: {str(e)}")