import math
import os

# ======================= 参数配置 =======================

# 粒子数量 (256*256)
NUM_PARTICLES_X = 64
NUM_PARTICLES_Y = 64
TOTAL_PARTICLES = NUM_PARTICLES_X * NUM_PARTICLES_Y

# 粒子半径
RADIUS = 1.0

# 目标体积分数 (Packing Fraction)
PACKING_FRACTION = 0.72

# 输出文件名（默认与 PACKING_FRACTION 耦合，可在函数调用时覆盖）
OUTPUT_FILENAME = f"config_phi{PACKING_FRACTION:.3f}_centered.dat"

# ==========================================================

def generate_centered_lattice_configuration(packing_fraction=PACKING_FRACTION, output_filename=None):
    """
    基于给定的物理参数，生成一个以原点为中心的 [-L/2, L/2) 
    正方形晶格初始位型文件。
    """
    print("开始生成初始位型...")
    print(f"  粒子总数: {TOTAL_PARTICLES}")
    print(f"  粒子半径: {RADIUS}")
    print(f"  体积分数: {packing_fraction}")
    print(f"  目标坐标系: [-L/2, L/2)")
    
    # 1. 计算单个粒子的面积
    area_per_particle = math.pi * RADIUS**2
    
    # 2. 计算所有粒子的总面积
    total_particle_area = TOTAL_PARTICLES * area_per_particle
    
    # 3. 根据体积分数，计算所需的总盒子面积
    total_box_area = total_particle_area / packing_fraction
    
    # 4. 计算正方形盒子的边长 (Lx = Ly = L)
    box_side_length = math.sqrt(total_box_area)
    half_box_side = 0.5 * box_side_length  # 预先乘法，循环内少一次除法
    
    print(f"计算出的盒子边长 L: {box_side_length:.8f}")
    
    # 5. 计算晶格间距
    lattice_spacing = box_side_length / NUM_PARTICLES_X
    
    # 安全性检查
    if lattice_spacing < 2.0 * RADIUS:
        print("\n!!! 警告: 晶格间距小于粒子直径，请检查参数。")
    else:
        print(f"计算出的晶格间距: {lattice_spacing:.8f} (安全)")

    # 6. 生成文件
    if output_filename is None:
        output_filename = f"config_phi{packing_fraction:.3f}_centered.dat"

    try:
        with open(output_filename, 'w') as f:
            # 写入第一行：盒子尺寸
            f.write(f"{box_side_length} {box_side_length}\n")
            
            # 遍历晶格并放置粒子
            particle_count = 0
            for i in range(NUM_PARTICLES_X):
                for j in range(NUM_PARTICLES_Y):
                    # a. 先在 [0, L) 坐标系中计算位置
                    x_positive = (i + 0.5) * lattice_spacing
                    y_positive = (j + 0.5) * lattice_spacing
                    
                    # b. === 关键修正: 将坐标平移到 [-L/2, L/2) ===
                    x_centered = x_positive - half_box_side
                    y_centered = y_positive - half_box_side
                    
                    # 写入修正后的粒子信息
                    f.write(f"{x_centered} {y_centered} {RADIUS}\n")
                    particle_count += 1
            
        print("\n生成成功！")
        print(f"总共写入 {particle_count} 个粒子。")
        print(f"配置文件已保存到: '{os.path.abspath(output_filename)}'")
        
    except Exception as e:
        print(f"\n生成文件时发生错误: {e}")

# --- 运行脚本 ---
if __name__ == "__main__":
    # 生成范围：0.69 到 0.72，步长 0.002（包含两端）
    start_phi = 0.69
    end_phi = 0.72
    step = 0.002
    n_steps = int(round((end_phi - start_phi) / step)) + 1

    for i in range(n_steps):
        phi = round(start_phi + i * step, 3)
        out_name = f"config_phi{phi:.3f}_centered.dat"
        print(f"\nGenerating configuration for phi={phi:.3f} -> {out_name}")
        generate_centered_lattice_configuration(packing_fraction=phi, output_filename=out_name)
