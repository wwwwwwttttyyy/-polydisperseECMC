# ECMC - Event-Chain Monte Carlo for Hard Disks

高效的二维硬盘系统 Event-Chain Monte Carlo (ECMC) 模拟程序,支持多分散系统和压强测量。

## 特性

- ✅ **高效的 Cell Grid 算法** - O(1) 碰撞检测
- ✅ **周期性边界条件** - 正确处理粒子跨边界运动
- ✅ **压强测量** - 基于碰撞累积的精确压强计算
- ✅ **多分散系统** - 支持粒子半径分布
- ✅ **快照输出** - 定期保存系统构型
- ✅ **统计分析** - 自动区分平衡期和产出期
- ✅ **二进制I/O** - 节省9%存储空间的二进制格式

## 构建

### 要求

- CMake 3.15+
- C++11 编译器 (GCC, Clang, MSVC)
- Ninja (可选,推荐)

### 编译步骤

```bash
# Windows (PowerShell)
mkdir build
cd build
cmake .. -G Ninja -DCMAKE_BUILD_TYPE=Release
ninja

# Linux/Mac
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j
```

## 使用方法

### 准备初始位型文件

程序需要一个初始构型文件作为输入。项目提供了示例文件 `config_initial.dat` (256个粒子, φ=0.70)。

**文本格式**:
```text
33.8958 33.8958              # 第1行: 盒子尺寸 Lx Ly
-15.8887 -15.8887 1.01186    # 第2行开始: x y radius (每个粒子一行)
-15.8887 -13.7702 0.987337
-15.8887 -11.6517 1.01119
...
```

**注意**: 坐标系统以盒子中心为原点,范围 [-Lx/2, Lx/2) × [-Ly/2, Ly/2)

### 1. 使用参数文件 (推荐)

创建参数文件 `params.txt`:

```text
input_file = config_initial.dat
n_chains = 1000000
sample_interval_pressure = 50000
sample_interval_snapshot = 100000
n_equilibration_pressure = 100000
rng_seed = 42
```

运行模拟:

```bash
./ecmc_main -p params.txt
```

### 2. 纯命令行模式

```bash
./ecmc_main -i initial.dat -n 1000000 -s 50000 -e 100000
```

### 3. 混合模式 (命令行覆盖参数文件)

```bash
./ecmc_main -p params.txt -n 2000000 -s 100000
```

## 参数说明

| 参数 | 短选项 | 默认值 | 说明 |
|------|--------|--------|------|
| `input_file` | `-i` | *必需* | 初始构型文件 |
| `output_file` | `-o` | `config_final.dat` | 输出构型文件 |
| `n_chains` | `-n` | 100000 | 事件链总数 |
| `chain_length` | `-l` | 1.0 | 链长度 (盒子尺寸倍数) |
| `sample_interval_pressure` | `-s` | 10000 | 压强采样间隔 |
| `sample_interval_snapshot` | `-S` | 0 | 快照输出间隔 (0=禁用) |
| `n_equilibration_pressure` | `-e` | 0 | 压强平衡期链数 |
| `rng_seed` | `-r` | *随机* | 随机数种子 |
| `binary_output` | `-b` | false | 使用二进制输出 |

## 输出文件

### 1. 最终构型 (`config_final.dat`)

包含模拟结束时的粒子位置和半径。

**文本格式**:
```text
256         # 粒子数
33.8958 33.8958    # 盒子尺寸
-3.142 5.678 1.000  # x y radius
...
```

**二进制格式**: 使用 `-b` 选项,节省约9%空间。

### 2. 压强数据 (`pressure_data.txt`)

包含完整的压强统计信息:

```text
# Production Statistics (after equilibration):
#   betaP* = 9.28974 +/- 0.0483076
#   Px     = 9.31525 +/- 0.0632583
#   Py     = 9.26423 +/- 0.0548891
#
# Sample  Chains  betaP*   Px       Py       Phase
0  50000  9.5959  9.6123  9.5795  equilibration
1  50000  9.3621  9.4012  9.3230  production
...
```

Notes:
- `Chains` is the cumulative chain count within the current phase.
- At the start of production sampling (after equilibration), `Chains` is reset and counted again from 0.
- `Phase` indicates whether the sample belongs to `equilibration` or `production`.

### 3. 快照文件 (`snapshot_*.dat`)

定期保存的构型文件,文件名格式: `snapshot_00000100000.dat`

使用 `-S` 选项启用,例如 `-S 100000` 每10万链保存一次。

## 压强计算公式

本程序使用基于碰撞的压强计算方法:

$$
\beta P^* = \left[\rho + \rho \cdot \frac{\sum_i \delta_i}{L}\right] \times (2\bar{r})^2
$$

其中:
- $\rho = N/V$ - 数密度
- $\delta_i$ - 第 $i$ 次碰撞的自由程减少量
- $L = n_{\text{chains}} \times l_{\text{chain}}$ - 总链长
- $\bar{r}$ - 平均粒子半径

## 测试

```bash
cd tests/build
cmake .. -G Ninja
ninja
ctest
```

测试包括:
- 单粒子移动
- 碰撞检测
- 周期性边界条件
- 高密度系统稳定性
- 多分散系统
- I/O 格式验证

## 性能建议

1. **平衡期设置**: `n_equilibration_pressure ≈ 0.1 × n_chains`
2. **采样频率**: `sample_interval_pressure ≈ 0.05 × n_chains`
3. **快照频率**: `sample_interval_snapshot ≈ 0.1 × n_chains` (或禁用以节省I/O)
4. **填充分数**: φ < 0.75 推荐 (避免玻璃化)

## 验证结果

对于 φ = 0.70 的硬盘系统:
- **测试结果**: βP* = 9.29 ± 0.05
- **Carnahan-Starling**: βP* ≈ 9.28
- **相对误差**: < 0.2%

## 参考文献

1. Bernard et al., *Phys. Rev. E* **80**, 056704 (2009)
2. Michel et al., *J. Chem. Phys.* **140**, 054116 (2014)
3. Kampmann et al., *J. Chem. Phys.* **147**, 214110 (2017)

## 许可证

MIT License

## 作者

卫天宇
