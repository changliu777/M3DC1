# M3DC1 中文说明

M3DC1 是一个面向磁约束聚变等离子体数值模拟的科研代码库。项目包含 M3D-C1 相关的核心求解、非结构网格支持、SCOREC/PUMI 接口以及用户文档，用于构建和运行高阶有限元、三维磁流体力学及相关物理模型的计算程序。

> 本说明主要帮助中文用户快速了解仓库结构和构建入口。更完整的物理模型、输入文件、运行作业和后处理说明，请参考 `doc/` 目录中的英文用户手册源码。

## 仓库结构

```text
.
├── CMakeLists.txt       # 顶层 CMake 配置
├── LICENSE              # 项目许可证
├── doc/                 # M3D-C1 用户手册的 LaTeX 源文件
├── m3dc1_scorec/        # M3DC1 与 SCOREC/PUMI 相关接口和库
├── skeleton/            # 示例/骨架程序
└── unstructured/        # 非结构网格相关代码、脚本和平台说明
```

## 主要依赖

该项目通常在 Linux/HPC 环境中构建，需要根据目标平台准备编译器、MPI 和科学计算库。常见依赖包括：

- CMake 3.8 或更高版本
- Fortran、C、C++ 编译器
- MPI
- PETSc
- SCOREC/PUMI
- METIS/ParMETIS
- Zoltan
- HDF5、LAPACK 等数值计算库

不同集群和编译器组合的配置差异较大，`unstructured/README/` 下提供了多个平台的构建说明，可作为配置环境变量和 CMake 参数的参考。

## 获取代码

```bash
git clone https://github.com/wjh301/M3DC1.git
cd M3DC1
```

如果使用 fork 参与开发，建议先添加上游仓库：

```bash
git remote add upstream https://github.com/wjh301/M3DC1.git
git fetch upstream
```

## 构建入口

顶层工程使用 CMake，并包含两个主要子目录：

- `m3dc1_scorec`
- `unstructured`

可从独立构建目录开始配置：

```bash
mkdir build
cd build
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_OPENMP=OFF \
  -DENABLE_COMPLEX=OFF \
  -DENABLE_3D=OFF \
  -DENABLE_ST=OFF \
  -DENABLE_PARTICLE=OFF
make -j
```

实际构建通常还需要显式指定 MPI、PETSc、SCOREC/PUMI、Zoltan、ParMETIS 等库的位置。请结合目标平台阅读：

- `m3dc1_scorec/README`
- `doc/installation.tex`
- `unstructured/README/` 下对应平台的说明文件

## 常用 CMake 选项

顶层 `CMakeLists.txt` 提供了一些功能开关：

- `ENABLE_OPENMP`：启用 OpenMP 支持
- `ENABLE_COMPLEX`：构建复数版本
- `ENABLE_3D`：构建三维版本
- `ENABLE_ST`：构建 stellarator 版本
- `ENABLE_PARTICLE`：启用动理学粒子模块

示例：

```bash
cmake .. -DENABLE_OPENMP=ON -DENABLE_3D=ON
```

## 生成用户手册

用户手册的 LaTeX 源文件位于 `doc/`。进入该目录后可使用：

```bash
cd doc
pdflatex M3DC1.tex
```

根据本地 LaTeX 环境，可能需要多次运行 `pdflatex` 或安装参考文献相关工具。

## 许可证

请参考根目录 `LICENSE` 以及 `m3dc1_scorec/LICENSE` 中的许可证文本。

## 贡献

建议通过 fork + pull request 的方式贡献代码或文档：

1. 从上游仓库同步最新代码。
2. 基于 `master` 创建功能分支。
3. 提交清晰、范围集中的修改。
4. 向上游仓库发起 pull request，并说明修改目的和验证方式。

