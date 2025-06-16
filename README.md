# MobileGNSS-SPP

## An EKF-based SPP system optimized for smartphone

[[中]](./README_CN.md) &ensp; [[EN]](./README.md)

![MobileGNSS-SPP框架](https://raw.githubusercontent.com/salmoshu/Winchell-ImgBed/main/img/20250616-204621.png)

MobileGNSS-SPP is an advanced, open-source Single Point Positioning (SPP) system built upon the foundation of [RTKLIB](https://www.rtklib.com/), tailored specifically for optimizing GNSS performance on smartphone chipsets. While initially developed for a specific GNSS chipset, the optimization techniques and engineering approaches are designed with generality and adaptability in mind, making them applicable to a wide range of GNSS hardware. The project prioritizes practical engineering solutions and algorithmic robustness across diverse real-world scenarios, serving as a comprehensive roadmap for transitioning from open-source code to production-ready GNSS solutions.

## Key Features
- **Smartphone-Optimized GNSS Processing**: Optimize the algorithm for mobile devices with low-quality data to improve positioning accuracy.
- **Robustness Across Scenarios**: Extensively tested in diverse environments, from open-sky highways to urban canyons with heavy multipath interference.
- **Comprehensive Testing Framework**: Includes Python-based tools for batch processing, accuracy evaluation, and in-depth data analysis.

## Algorithmic Optimizations
The project introduces several key enhancements to the RTKLIB codebase, primarily in the `rtklib_src/pntpos.c` file, to improve positioning accuracy and robustness:

- **M-Estimation (Robust Estimation)**: Employs iterative least squares with a Huber kernel-based weight matrix (W) to mitigate outliers, incorporating truncation for large residuals to enhance stability.
- **Zero-Velocity Correction**: Improves positioning accuracy during static or low-motion scenarios, critical for smartphone use cases.
- **SNR-Weighted Model**: Applies Signal-to-Noise Ratio (SNR) weighting to both pseudorange and Doppler measurements, improving signal quality assessment.
- **Multipath Error Compensation**: Leverages empirical chipset testing to compensate for pseudorange residuals caused by multipath effects, enhancing performance in complex environments.
- **Adaptive Q Matrix**: Dynamically adjusts the velocity covariance in the Extended Kalman Filter (EKF) using the velocity difference between EKF predictions and Robust Weighted Least Squares (RWLS) estimates.
- **Cost Minimization via Quadratic Programming**: Implements a post-processing optimization technique based on quadratic programming to minimize positioning errors. This part is implemented using Python and is only for post-processing.

For more technical details, please refer to [the online documentation](https://salmoshu.github.io/algorithm/MobileGNSS-SPP/).

## 1. Building and Running MobileGNSS-SPP

### 1.1 Build Environment
MobileGNSS-SPP uses [CMake](https://cmake.org/) for cross-platform build management, supporting compilation on **Linux**, **macOS**, and **Windows**. For optimal compatibility and ease of use, we recommend compiling on Windows using Microsoft Visual Studio (VS).

### 1.2 Compiling on Windows
The `rnx2rtkp` application includes a pre-configured Visual Studio project file (`msc`) for seamless integration. Follow these steps to build the project:

1. Open the `rnx2rtkp` project in Visual Studio.
2. Switch the build configuration to **Release** mode and set the platform to **Win32** to avoid compatibility issues.
3. Configure the command-line arguments in Visual Studio:
   - Navigate to **Configuration Properties > Debugging > Command Arguments**.
   - Ensure the configuration is set to **Release**.
   - Add the following command-line arguments:

```shell
.\rnx2rtkp -x 0 -k ..\conf\rover.conf -o ..\..\..\data\01-opensky\data01\rover.pos ..\..\..\data\01-opensky\data01\rover.obs ..\..\..\data\01-opensky\data01\rover.nav
```

4. Build the solution to generate the executable.

## 2. Using MobileGNSS-SPP

### 2.1 Test Scenarios
MobileGNSS-SPP has been rigorously tested across diverse environments to ensure robustness and reliability. The following scenarios are supported for evaluation:

| Test Scenario | Description |
|---------------|-------------|
| **Open-Sky (Highway)** | High-speed testing on elevated highways or outer ring roads with unobstructed GNSS signals, ideal for evaluating performance in optimal conditions. |
| **Urban Streets (Tree-Lined)** | Testing in city streets with tree cover, introducing moderate GNSS signal occlusion for assessing performance under partial interference. |
| **Complex Urban (Downtown)** | Testing in dense urban environments with high-rise buildings, heavy multipath, and signal obstructions, challenging the system's robustness. |
| **Occluded Environment (Underpass)** | Testing near or under elevated structures (e.g., bridges), where GNSS signals face significant blockage and multipath effects. |

### 2.2 Test and Analysis Scripts
The project includes a suite of Python scripts to facilitate testing, evaluation, and data analysis:

```plaintext
\python
├── rnx2rtkp_batch.py     : Batch processes rnx2rtkp for multiple test scenarios.
├── scores.py             : Computes positioning accuracy metrics for algorithm evaluation.
├── scores_batch.py       : Aggregates accuracy metrics across multiple test runs (run after rnx2rtkp_batch.py).
├── data_analysis         : Tools for in-depth GNSS data analysis.
│   ├── 2.4#pr_doppler_corr.py : Analyzes pseudorange and Doppler clock drift (see Section 2.4 of the documentation).
│   └── 2.5#prr_tdcp.py       : Examines Doppler and Time-Differenced Carrier Phase (TDCP) correlations (see Section 2.5).
├── mincost               : Quadratic programming-based global optimizer for post-processing.
└── rtklipy              : Python implementation of RTKLIB for additional flexibility.
```

### 2.3 Running Tests
1. Prepare test data in the `data` directory, organized by scenario (e.g., `data/01-opensky/`).
2. Execute `rnx2rtkp_batch.py` to process test data across multiple scenarios.
3. Run `scores_batch.py` to evaluate positioning accuracy and generate performance reports.
4. Use the `data_analysis` scripts for detailed insights into pseudorange, Doppler, and multipath effects.

## 3. License
MobileGNSS-SPP is licensed under the [MIT License](LICENSE). See the `LICENSE` file for details.

## 4. Acknowledgments
- Built upon the robust foundation of [RTKLIB](https://www.rtklib.com/).
- Gratitude to the participants of the Google Decimeter Challenge for their excellent code contributions: [@taroz1461](https://www.kaggle.com/taroz1461), [@saitodevel01](https://www.kaggle.com/saitodevel01), [@timeverett](https://www.kaggle.com/timeverett).
- Special thanks to the GNSS research community for providing valuable insights and test methodologies.

## 5. Contact
For questions, bug reports, or feature requests, please open an issue on the [GitHub repository](https://github.com/salmoshu/MobileGNSS-SPP). For general inquiries, contact [winchell.hu@outlook.com](mailto:winchell.hu@outlook.com).
