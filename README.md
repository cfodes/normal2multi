 (cd "$(git rev-parse --show-toplevel)" && git apply --3way <<'EOF' 
diff --git a/README.md b/README.md
new file mode 100644
index 0000000000000000000000000000000000000000..7d005c9796751aa2145a3aa58222e6ae0a875ab6
--- /dev/null
+++ b/README.md
@@ -0,0 +1,42 @@
+# normal2multi
+
+## 项目简介
+`normal2multi` 是一个基于 C++17 的网格多级分区与径向基函数（RBF）变形示例项目，结合 METIS 分区与 Eigen 线性代数库，对包含翼型等壁面边界的网格执行分层分区、静/动边界识别以及分级 RBF 变形。
+
+项目入口位于 `src/main.cpp`，通过构造 `TestDriver` 对象串联起网格读取、预处理、多级分区 RBF 计算与结果输出的完整流程。默认示例使用 `data/` 下的 NACA0012 C 型网格，并将变形后的结果写入 `output/` 目录。【F:src/main.cpp†L1-L22】【F:include/test_driver.hpp†L12-L36】
+
+## 功能概览
+- **多级 METIS 分区：** `multi_partition` 类根据用户配置的分区数量，对壁面单元进行层级化划分，为后续分块 RBF 提供数据结构支持。【F:include/multi_partition.hpp†L16-L70】
+- **分级 RBF 变形：** 项目实现了从全局 RBF 到分块 DRRBF 的多级变形流程，结合节点等级与动静边界距离，提高变形效率与稳定性。【F:include/multi_partition.hpp†L71-L162】
+- **网格预处理与驱动：** `TestDriver` 负责网格读取、物面树构建、节点等级映射、调用多级分区算法以及结果写出，方便集成不同的算例与参数配置。【F:include/test_driver.hpp†L18-L36】
+
+## 目录结构
+- `src/`：主程序及实现文件。
+- `include/`：核心算法类与数据结构的头文件。
+- `data/`：示例输入网格（SU2 格式）。
+- `output/`：默认输出目录，用于保存变形后的网格文件。
+- `CMakeLists.txt`：构建脚本，配置 C++17、METIS 与 Eigen 依赖。【F:CMakeLists.txt†L1-L27】
+
+## 环境依赖
+- C++17 编译器（如 GCC 9+/Clang 10+）。
+- CMake ≥ 3.10。
+- [Eigen 3.3](https://eigen.tuxfamily.org/)。
+- METIS（Linux 下可通过 `libmetis-dev` 安装）。构建脚本会在找不到 METIS 时给出错误提示。【F:CMakeLists.txt†L8-L27】
+
+## 构建与运行
+```bash
+mkdir -p build && cd build
+cmake ..
+make -j
+./MyProject
+```
+可根据需要修改 `src/main.cpp` 中的输入/输出路径与分区、误差容限配置。【F:src/main.cpp†L10-L21】
+
+## 数据与输出
+- 默认输入：`data/naca0012_Cmesh.su2`。
+- 默认输出：`output/naca0012_Cmesh_result.su2`。
+
+运行后可在 `output/` 目录查看生成的 SU2 网格文件，并进一步用于后处理或可视化。【F:src/main.cpp†L10-L21】
+
+## 许可证
+项目当前未包含许可证文件，如需在其他项目中使用请与作者确认授权方式。 
EOF
)
