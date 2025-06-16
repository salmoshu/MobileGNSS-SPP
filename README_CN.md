# MobileGNSS-SPP

## An EKF-based SPP system optimized for smartphone

[[中]](./README_CN.md) &ensp; [[EN]](./README.md)

![MobileGNSS-SPP框架](https://raw.githubusercontent.com/salmoshu/Winchell-ImgBed/main/img/MobileGNSS-SPP%20Framework-20250616-150022.png)

MobileGNSS-SPP 是一个基于 [RTKLIB](https://www.rtklib.com/) 的开源单点定位（SPP）项目，专门针对智能手机 GNSS 芯片进行了优化。虽然项目最初为特定 GNSS 芯片设计，但其优化思路具有通用性和启发性，可广泛应用于其他 GNSS 硬件。该项目不以算法或框架的先进性为核心，而是聚焦于工程化实现以及算法在多样化场景下的鲁棒性，因此可作为从开源代码到工程化方案的参考路线图。

## 核心特性
- **智能手机优化 GNSS 处理**：针对低质量数据的移动设备，优化算法以提升定位精度。
- **多样场景鲁棒性**：在开阔高速公路到复杂城市峡谷等多种环境中经过广泛测试。
- **全面测试框架**：提供 Python 工具集，用于批量处理、精度评估和深入数据分析。

## 算法优化
项目对 RTKLIB 核心代码（主要在 `rtklib_src/pntpos.c` 文件）进行了多项改进，以提升定位精度和鲁棒性：

- **M估计（抗差估计）**：采用迭代最小二乘法构建权阵 W，使用 Huber 核函数，并对大残差项进行截断处理以提高稳定性。
- **零速修正**：优化静态或低速场景下的定位精度，适用于智能手机典型用例。
- **SNR加权模型**：对伪距和多普勒测量值均应用信噪比（SNR）加权，提升信号质量评估。
- **多路径误差补偿**：基于芯片测试经验，补偿伪距残差中的多路径效应，增强复杂环境下的性能。
- **自适应 Q 矩阵**：通过扩展卡尔曼滤波（EKF）预测速度与鲁棒加权最小二乘（RWLS）速度差值，动态调整 EKF 中速度协方差。
- **基于二次规划的代价最小化**：采用二次规划优化方法，减少定位误差（该部分使用Python实现，仅针对后处理）。

对于更多的技术细节，请查看[在线文档](https://salmoshu.github.io/algorithm/MobileGNSS-SPP/)。

## 1. 编译与运行

### 1.1 编译环境
MobileGNSS-SPP 使用 [CMake](https://cmake.org/) 进行跨平台构建管理，理论上支持 **Linux**、**macOS** 和 **Windows** 环境。为确保最佳兼容性，推荐在 Windows 环境下使用 Microsoft Visual Studio 进行编译。

### 1.2 在 Windows 下编译
`rnx2rtkp` 应用包含一个预配置的 Visual Studio 项目文件（`msc`），可直接打开使用。编译步骤如下：

1. 在 Visual Studio 中打开 `rnx2rtkp` 项目。
2. 将构建配置切换为 **Release** 模式，并设置平台为 **Win32**，以避免兼容性问题。
3. 配置命令行参数：
   - 进入 **配置属性 > 调试 > 命令参数**。
   - 确保配置为 **Release** 模式。
   - 输入以下命令行参数：

```shell
.\rnx2rtkp -x 0 -k ..\conf\rover.conf -o ..\..\..\data\01-opensky\data01\rover.pos ..\..\..\data\01-opensky\data01\rover.obs ..\..\..\data\01-opensky\data01\rover.nav
```

4. 构建解决方案以生成可执行文件。

## 2. 使用 MobileGNSS-SPP

### 2.1 测试场景
MobileGNSS-SPP 在多种环境中进行了严格测试，以确保鲁棒性和可靠性。支持的测试场景包括：

| 测试场景 | 描述 |
|----------|------|
| **开阔环境（高速公路）** | 在高架桥面、城市外环路等开阔直线路段测试，GNSS 信号接收条件优越，适合评估高速移动和无障碍环境下的性能。 |
| **城市街道（树荫遮挡）** | 在有路边树木的城市街道测试，存在一定程度的 GNSS 信号遮挡，评估部分干扰下的性能。 |
| **复杂城市环境（市中心）** | 在高楼林立、建筑密集的市中心测试，GNSS 信号受多路径和遮挡干扰，挑战系统鲁棒性。 |
| **遮挡环境（高架桥下）** | 在高架桥下或附近测试，桥墩、桥面及周边建筑导致信号遮挡和反射，形成复杂信号环境。 |

### 2.2 测试与分析脚本
项目提供一系列 Python 脚本，便于测试、评估和数据分析：

```plaintext
\python
├── rnx2rtkp_batch.py     : 批量执行 rnx2rtkp，处理多个测试场景。
├── scores.py             : 计算算法定位精度的评估指标。
├── scores_batch.py       : 批量处理多个测试场景的精度评估（在 rnx2rtkp_batch.py 后运行）。
├── data_analysis         : 数据分析工具。
│   ├── 2.4#pr_doppler_corr.py : 伪距与多普勒钟漂分析（参考文档 2.4 节）。
│   └── 2.5#prr_tdcp.py       : 多普勒与时间差分载波相位（TDCP）相关性分析（参考文档 2.5 节）。
├── mincost               : 基于二次规划的全局优化器，用于后处理。
└── rtklipy              : RTKLIB 的 Python 实现，增加灵活性。
```

### 2.3 运行测试
1. 在 `data` 目录中准备测试数据，按场景组织（例如 `data/01-opensky/`）。
2. 运行 `rnx2rtkp_batch.py` 处理多个场景的测试数据。
3. 运行 `scores_batch.py` 评估定位精度并生成性能报告。
4. 使用 `data_analysis` 脚本深入分析伪距、多普勒和多路径效应。

## 3. 许可证
MobileGNSS-SPP 采用 [MIT 许可证](LICENSE)。详情见 `LICENSE` 文件。

## 4. 致谢
- 基于 [RTKLIB](https://www.rtklib.com/) 的强大基础开发。
- 感谢谷歌分米挑战大赛参赛选手提供的优秀代码：[@taroz1461](https://www.kaggle.com/taroz1461), [@saitodevel01](https://www.kaggle.com/saitodevel01), [@timeverett](https://www.kaggle.com/timeverett)。
- 感谢 GNSS 研究社区提供的宝贵见解和测试方法。

## 5. 联系方式
如有问题、错误报告或功能请求，请在 [GitHub 仓库](https://github.com/salmoshu/MobileGNSS-SPP) 上提交 Issue。通用问题可联系 [winchell.hu@outlook.com](mailto:winchell.hu@outlook.com)。
