import numpy as np

def read_ecmc_config(filename):
    """读取ECMC配置文件"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # 第一行是盒子尺寸
    box_size = list(map(float, lines[0].strip().split()))
    Lx, Ly = box_size[0], box_size[1]
    
    # 后续行是粒子数据
    particles = []
    for line in lines[1:]:
        if line.strip():
            x, y, radius = map(float, line.strip().split())
            particles.append((x, y, radius))
    
    return Lx, Ly, particles

def calculate_volume_fraction(filename):
    """计算配置文件的体积分数"""
    Lx, Ly, particles = read_ecmc_config(filename)
    
    # 计算体积分数
    total_area = sum(np.pi * r**2 for _, _, r in particles)
    volume_fraction = total_area / (Lx * Ly)
    
    return volume_fraction, len(particles), Lx, Ly

def main():
    """主函数，检查所有测试配置的体积分数"""
    config_files = [
        'test_config1.dat',
        'test_config2.dat', 
        'test_config3.dat',
        'test_config4.dat'
    ]
    
    print("Checking ECMC configuration volume fractions...")
    print("=" * 60)
    
    for filename in config_files:
        try:
            volume_fraction, num_particles, Lx, Ly = calculate_volume_fraction(filename)
            print(f"\n{filename}:")
            print(f"  Box size: {Lx:.4f} x {Ly:.4f}")
            print(f"  Particles: {num_particles}")
            print(f"  Volume fraction: {volume_fraction:.4f}")
        except Exception as e:
            print(f"  Error checking {filename}: {e}")

if __name__ == "__main__":
    main()
