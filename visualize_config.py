"""
可视化 ECMC 测试的粒子位型
读取 .dat 文件并绘制粒子圆盘
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl
from matplotlib.collections import EllipseCollection
import numpy as np
import sys

def read_configuration(filename):
    """读取位型文件
    
    格式：
        第1行: Lx Ly
        后续行: x y radius
    """
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # 第一行：盒子大小
    Lx, Ly = map(float, lines[0].split())
    
    # 读取粒子数据
    particles = []
    for line in lines[1:]:
        line = line.strip()
        if line == '' or line.startswith('#'):
            continue
        x, y, r = map(float, line.split())
        particles.append((x, y, r))
    
    return Lx, Ly, particles

def plot_configuration(filename, title="Particle Configuration", save_as=None):
    """绘制位型"""
    Lx, Ly, particles = read_configuration(filename)
    
    # Use GridSpec to give the colorbar its own narrow axis so it won't cover data
    fig = plt.figure(figsize=(10, 10))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 0.05], wspace=0.05)
    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[0, 1])
    
    # 设置坐标轴范围（中心对称盒子：[-L/2, L/2]）
    ax.set_xlim(-Lx/2, Lx/2)
    ax.set_ylim(-Ly/2, Ly/2)
    ax.set_aspect('equal')
    
    # 绘制粒子
    if particles:
        particles_arr = np.array(particles)
        # EllipseCollection 需要直径
        diameters = 2 * particles_arr[:, 2]
        radii = particles_arr[:, 2]
        offsets = particles_arr[:, :2]
        
        ec = EllipseCollection(
            widths=diameters, 
            heights=diameters, 
            angles=0, 
            units='x',
            offsets=offsets, 
            transOffset=ax.transData,
            edgecolor='none', 
            cmap='coolwarm',
            array=radii,
            linewidth=0,
            alpha=0.7
        )
        ax.add_collection(ec)
        cbar = fig.colorbar(ec, cax=cax)
        cbar.set_label('Radius')
    
    # 绘制盒子边界（中心对称）
    box = patches.Rectangle((-Lx/2, -Ly/2), Lx, Ly, 
                            linewidth=2, 
                            edgecolor='red', 
                            facecolor='none')
    ax.add_patch(box)
    
    ax.set_xlabel('X', fontsize=14)
    ax.set_ylabel('Y', fontsize=14)
    ax.set_title(f'{title}\nN={len(particles)}, Box=[{Lx:.1f}, {Ly:.1f}] (centered at origin)', 
                 fontsize=16)
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax.axvline(x=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    
    fig.tight_layout()
    
    if save_as:
        plt.savefig(save_as, dpi=150, bbox_inches='tight')
        print(f"Figure saved to {save_as}")
    
    plt.show()

def compare_configurations(file1, file2, title1="Initial", title2="Final", save_as=None):
    """对比两个位型"""
    Lx1, Ly1, particles1 = read_configuration(file1)
    Lx2, Ly2, particles2 = read_configuration(file2)
    
    # 计算全局半径范围以统一颜色映射
    all_radii = []
    if particles1: all_radii.extend([p[2] for p in particles1])
    if particles2: all_radii.extend([p[2] for p in particles2])
    
    norm = None
    cmap = 'coolwarm'
    if all_radii:
        vmin, vmax = min(all_radii), max(all_radii)
        if vmin == vmax:
            pad = abs(vmin) * 0.01 or 1e-6  # 避免色标上下限相等导致 colorbar 出错
            vmin -= pad
            vmax += pad
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    
    fig = plt.figure(figsize=(16, 8))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.05], wspace=0.12)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    cax = fig.add_subplot(gs[0, 2])
    ec_list = []  # store collections for shared colorbar
    
    for ax, Lx, Ly, particles, title in [
        (ax1, Lx1, Ly1, particles1, title1),
        (ax2, Lx2, Ly2, particles2, title2)
    ]:
        # 设置坐标轴范围（中心对称盒子：[-L/2, L/2]）
        ax.set_xlim(-Lx/2, Lx/2)
        ax.set_ylim(-Ly/2, Ly/2)
        ax.set_aspect('equal')
        
        # 绘制粒子
        if particles:
            particles_arr = np.array(particles)
            diameters = 2 * particles_arr[:, 2]
            radii = particles_arr[:, 2]
            offsets = particles_arr[:, :2]
            
            ec = EllipseCollection(
                widths=diameters, 
                heights=diameters, 
                angles=0, 
                units='x',
                offsets=offsets, 
                transOffset=ax.transData,
                edgecolor='none', 
                cmap=cmap,
                norm=norm,
                array=radii,
                linewidth=0,
                alpha=0.7
            )
            ax.add_collection(ec)
            ec_list.append(ec)
        
        # 绘制盒子边界（中心对称）
        box = patches.Rectangle((-Lx/2, -Ly/2), Lx, Ly, 
                                linewidth=2, 
                                edgecolor='red', 
                                facecolor='none')
        ax.add_patch(box)
        
        ax.set_xlabel('X', fontsize=12)
        ax.set_ylabel('Y', fontsize=12)
        ax.set_title(f'{title}\nN={len(particles)}', fontsize=14)
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
        ax.axvline(x=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)

    if norm and ec_list:
        sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])  # required for older matplotlib versions
        cbar = fig.colorbar(sm, cax=cax)
        cbar.set_label('Radius')
    
    fig.tight_layout()
    
    if save_as:
        plt.savefig(save_as, dpi=150, bbox_inches='tight')
        print(f"Comparison saved to {save_as}")
    
    plt.show()

if __name__ == "__main__":
    import os

    def auto_search():
        """Find default initial/final configuration files (no CLI args needed)."""
        search_dirs = [
            ".",                             # 当前目录
            "./snapshot",                    # 可选 snapshot 目录
            os.path.join(os.path.dirname(__file__), "build"),
        ]
        init_file = None
        final_file = None
        for dir_path in search_dirs:
            cand_init = os.path.join(dir_path, "config_phi0.690_centered.dat")
            cand_final = os.path.join(dir_path, "finalconf.dat")
            if os.path.exists(cand_init):
                init_file = cand_init
            if os.path.exists(cand_final):
                final_file = cand_final
            if init_file and final_file:
                break
        return init_file, final_file, search_dirs

    initial_file, final_file, search_dirs = auto_search()

    if initial_file and final_file:
        print(f"Found both configurations in: {os.path.dirname(initial_file)}")
        print("Plotting comparison...")
        compare_configurations(initial_file, final_file,
                              title1="Initial Configuration",
                              title2="Final Configuration (after ECMC)",
                              save_as="comparison.png")
    elif final_file:
        print(f"Found final configuration in: {os.path.dirname(final_file)}")
        plot_configuration(final_file, 
                          title="Final Configuration",
                          save_as="final.png")
    elif initial_file:
        print(f"Found initial configuration in: {os.path.dirname(initial_file)}")
        plot_configuration(initial_file, 
                          title="Initial Configuration",
                          save_as="initial.png")
    else:
        print("No configuration files found!")
        print(f"Searched in: {search_dirs}")
        print("Please run the test_high_density executable or place initialconf.dat/finalconf.dat in the search paths.")
        sys.exit(1)
