import numpy as np
from math import sqrt
import os

def debug_file_structure(file_path):
    """调试函数：查看文件结构"""
    print(f"调试文件: {file_path}")
    with open(file_path, 'r') as file:
        for i, line in enumerate(file):
            if i < 50:  # 只看前50行
                print(f"Line {i+1}: {line.strip()}")
            if "N O D A L" in line or "E L E M E N T" in line or "D I S P L A C E M E N T S" in line:
                print(f"Line {i+1}: {line.strip()}")
                # 显示接下来的几行
                for j in range(3):
                    try:
                        next_line = next(file)
                        print(f"Line {i+j+2}: {next_line.strip()}")
                    except StopIteration:
                        break

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
        lines = file.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        # 识别当前解析的部分
        if "N O D A L   P O I N T   D A T A" in line:
            current_section = "node_data"
            i += 1
            # 跳过标题行
            while i < len(lines) and ("NODE" not in lines[i] or "BOUNDARY" not in lines[i]):
                i += 1
            i += 1  # 跳过标题行
            
            # 读取节点数据
            while i < len(lines):
                line = lines[i].strip()
                if not line or "EQUATION" in line or "E L E M E N T" in line:
                    break
                
                parts = line.split()
                if len(parts) >= 7:  # 至少需要7个部分：node_id, 3个边界条件码, 3个坐标
                    try:
                        node_id = int(parts[0])
                        # 坐标在最后三个位置
                        x = float(parts[-3].replace('D', 'e'))
                        y = float(parts[-2].replace('D', 'e'))
                        node_coords[node_id] = (x, y)
                    except (ValueError, IndexError) as e:
                        # 只在调试模式下输出错误
                        pass
                i += 1
            continue
            
        elif "E L E M E N T   I N F O R M A T I O N" in line:
            current_section = "element_info"
            i += 1
            # 寻找单元表格开始 - 跳过标题行
            while i < len(lines) and not ("ELEMENT" in lines[i] and "NODE" in lines[i] and "MATERIAL" in lines[i]):
                i += 1
            i += 1  # 跳过标题行 "NUMBER-N      I        J        K       SET NUMBER"
            
            # 读取单元数据
            while i < len(lines):
                line = lines[i].strip()
                if not line or "YOUNG'S" in line or "D I S P L A C E M E N T S" in line or "M A T E R I A L" in line:
                    break
                
                parts = line.split()
                # 跳过标题行
                if parts and parts[0] == "NUMBER-N":
                    i += 1
                    continue
                    
                if len(parts) >= 5:  # 单元号 + 3个节点 + 材料集
                    try:
                        # 检查是否为数字
                        int(parts[0])  # 测试第一个是否为数字
                        # 读取三个节点ID (跳过单元号)
                        nodes = [int(parts[1]), int(parts[2]), int(parts[3])]
                        elements.append(nodes)
                    except (ValueError, IndexError):
                        # 跳过非数据行
                        pass
                i += 1
            continue
            
        elif "D I S P L A C E M E N T S" in line:
            current_section = "displacements"
            i += 1
            # 跳过标题行
            while i < len(lines) and "NODE" not in lines[i]:
                i += 1
            i += 1  # 跳过标题行
            
            # 读取位移数据
            while i < len(lines):
                line = lines[i].strip()
                if not line or "S T R E S S" in line:
                    break
                
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        node_id = int(parts[0])
                        ux = float(parts[1].replace('D', 'e'))
                        uy = float(parts[2].replace('D', 'e'))
                        node_displacements[node_id] = (ux, uy)
                    except (ValueError, IndexError):
                        pass
                i += 1
            continue
            
        elif "YOUNG'S     POISSON'S     THICKNESS" in line:
            current_section = "material"
            i += 1
            # 跳过可能的额外标题行
            while i < len(lines) and ("NUMBER" in lines[i] or "MODULUS" in lines[i] or "E" == lines[i].strip()):
                i += 1
            
            # 读取材料参数
            if i < len(lines):
                line = lines[i].strip()
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        # 检查是否包含数字
                        if 'e' in parts[1] or parts[1].replace('.','').replace('-','').isdigit():
                            material_params["E"] = float(parts[1].replace('D', 'e'))
                        if len(parts) > 2 and ('e' in parts[2] or parts[2].replace('.','').replace('-','').isdigit()):
                            material_params["nu"] = float(parts[2].replace('D', 'e'))
                    except (ValueError, IndexError):
                        pass
        
        i += 1
    
    return node_coords, node_displacements, elements, material_params

# 解析解函数 - 根据您提供的公式
def u_exact(x, y):
    L = 10
    h = 2
    E = 1e4
    nu = 0
    F = 1
    b = 1
    I = b * h**3 / 12
    return F * h / (E * I) * x * y

def v_exact(x):
    L = 10
    h = 2
    E = 1e4
    nu = 0
    F = 1
    b = 1
    I = b * h**3 / 12
    return -F * h / (2 * E * I) * x**2

# 三角形形函数
def shape_functions(xi, eta):
    return np.array([xi, eta, 1 - xi - eta])

# 雅可比行列式计算
def jacobian_det(nodes):
    nodes = np.array(nodes)
    J = np.zeros((2, 2))
    J[0, 0] = nodes[0, 0] - nodes[2, 0]
    J[0, 1] = nodes[1, 0] - nodes[2, 0]
    J[1, 0] = nodes[0, 1] - nodes[2, 1]
    J[1, 1] = nodes[1, 1] - nodes[2, 1]
    return abs(np.linalg.det(J))

def calculate_l2_error(file_path):
    # 解析输出文件
    node_coords, node_displacements, elements, material_params = parse_fe_output(file_path)
    
    print(f"解析文件: {file_path}")
    print(f"节点数量: {len(node_coords)}")
    print(f"单元数量: {len(elements)}")
    
    if len(node_coords) == 0 or len(elements) == 0:
        print(f"警告: 未能正确解析文件数据")
        return None
    
    # 材料参数（使用您提供的值）
    L = 2
    h = 1
    E = 1e4
    nu = 0
    F = 1
    b = 1
    I = b * h**3 / 12
    
    # 6点高斯积分点（精确积分四次多项式）
    gauss_points = [
        # (xi, eta, weight)
        (0.091576213509770, 0.091576213509770, 0.054975871827661),
        (0.816847572980459, 0.091576213509770, 0.054975871827661),
        (0.091576213509770, 0.816847572980459, 0.054975871827661),
        (0.108103018168070, 0.445948490915965, 0.111690794839005),
        (0.445948490915965, 0.108103018168070, 0.111690794839005),
        (0.445948490915965, 0.445948490915965, 0.111690794839005)
    ]
    
    # 初始化误差平方和
    error_squared = 0.0
    
    # 遍历每个单元
    for elem_id, elem in enumerate(elements):
        # 获取单元节点坐标和位移
        node_ids = elem
        
        # 检查节点是否存在
        if not all(nid in node_coords and nid in node_displacements for nid in node_ids):
            print(f"警告: 单元 {elem_id} 的节点数据不完整")
            continue
            
        nodes = [node_coords[n] for n in node_ids]
        disps = [node_displacements[n] for n in node_ids]
        
        # 计算单元雅可比行列式（用于面积计算）
        jac_det = jacobian_det(nodes)
        
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
            u_true = u_exact(x, y)
            v_true = v_exact(x)
            
            # 计算误差
            error_u = u_fem - u_true
            error_v = v_fem - v_true
            
            # 累积误差平方和
            error_squared += weight * jac_det * (error_u**2 + error_v**2)
    
    # 计算L₂范数误差
    l2_error = sqrt(error_squared)
    
    # 输出结果
    print(f"材料参数: E = {E:.1e}, nu = {nu:.1f}")
    print(f"加载力: F = {F:.1f}")
    print(f"截面惯性矩: I = {I:.5e}")
    print(f"\n位移L₂范数误差: {l2_error:.5e}")
    return l2_error

def process_all_files():
    """处理dats文件夹中的所有.out文件"""
    dats_folder = "/root/homework/finite_element/STAPpp/data/convergence/dats"
    results = {}
    
    # 获取所有.out文件
    out_files = [f for f in os.listdir(dats_folder) if f.endswith('.out')]
    out_files.sort()  # 按文件名排序
    
    print("="*60)
    print("有限元L₂范数误差计算")
    print("="*60)
    
    for out_file in out_files:
        file_path = os.path.join(dats_folder, out_file)
        print(f"\n处理文件: {out_file}")
        print("-"*40)
        
        try:
            error = calculate_l2_error(file_path)
            results[out_file] = error
        except Exception as e:
            print(f"处理文件 {out_file} 时出错: {str(e)}")
            # 调试文件结构
            try:
                debug_file_structure(file_path)
            except:
                pass
            results[out_file] = None
    
    # 输出汇总结果
    print("\n" + "="*60)
    print("汇总结果")
    print("="*60)
    for file_name, error in results.items():
        if error is not None:
            print(f"{file_name:<15}: L₂误差 = {error:.5e}")
        else:
            print(f"{file_name:<15}: 计算失败")
    
    return results

def calculate_convergence_rate(results):
    """计算收敛率"""
    print("\n" + "="*60)
    print("收敛率分析")
    print("="*60)
    
    # 提取有效的误差值
    valid_results = [(name, error) for name, error in results.items() if error is not None]
    valid_results.sort()  # 按文件名排序
    
    if len(valid_results) < 2:
        print("需要至少两个有效结果才能计算收敛率")
        return
    
    # 假设网格尺寸比例为 h1:h2:h3 = 1:0.5:0.25
    h_values = [1.0, 0.5, 0.25, 0.125]  # 对应1.out, 2.out, 3.out的网格尺寸
    errors = [result[1] for result in valid_results]
    
    print("网格尺寸h\t误差")
    print("-" * 30)
    for i, (name, error) in enumerate(valid_results):
        if i < len(h_values):
            print(f"{h_values[i]:.3f}\t\t{error:.5e}")
    
    # 计算收敛率
    if len(valid_results) >= 2:
        for i in range(len(valid_results) - 1):
            if i + 1 < len(h_values):
                h1, h2 = h_values[i], h_values[i + 1]
                e1, e2 = errors[i], errors[i + 1]
                
                # 收敛率 p = log(e1/e2) / log(h1/h2)
                convergence_rate = np.log(e1/e2) / np.log(h1/h2)
                print(f"\n从{valid_results[i][0]}到{valid_results[i+1][0]}的收敛率: {convergence_rate:.2f}")
    
    # 总体收敛率（使用最粗和最细网格）
    if len(valid_results) >= 3:
        h_coarse, h_fine = h_values[0], h_values[-1]
        e_coarse, e_fine = errors[0], errors[-1]
        overall_rate = np.log(e_coarse/e_fine) / np.log(h_coarse/h_fine)
        print(f"\n总体收敛率: {overall_rate:.2f}")
        print(f"理论收敛率（线性单元）: 2.0")
        
        if overall_rate >= 1.8:
            print("✓ 收敛率良好，接近理论值")
        else:
            print("⚠ 收敛率偏低，可能需要检查网格质量或边界条件")

# 主程序
if __name__ == "__main__":
    # 可以选择处理单个文件或所有文件
    
    # 选项1: 处理单个文件并调试
    # single_file = "/root/homework/finite_element/STAPpp/data/convergence/dats/1.out"
    # try:
    #     print("调试单个文件:")
    #     debug_file_structure(single_file)
    #     print("\n" + "="*60)
    #     error = calculate_l2_error(single_file)
    #     print(f"单文件误差值: {error:.5e}" if error else "计算失败")
    # except Exception as e:
    #     print(f"单文件计算出错: {str(e)}")
    
    # print("\n" + "="*80)
    
    # 选项2: 处理所有文件
    try:
        results = process_all_files()
        # 添加收敛率分析
        calculate_convergence_rate(results)
    except Exception as e:
        print(f"批量处理出错: {str(e)}")