cmake_minimum_required(VERSION 3.27)
project(AuxDM)

set(CMAKE_CXX_STANDARD 17)

# Append the path for LibTorch
list(APPEND CMAKE_PREFIX_PATH "/Users/argha/ExternalLibraries/libtorch")
# Find LibTorch package
find_package(Torch REQUIRED)

find_library(ACCELERATE_FRAMEWORK Accelerate)

# Set OpenMP paths
set(LIBOMP_PREFIX "/opt/homebrew/opt/libomp")

# Set the OpenMP include and library paths
include_directories(/opt/homebrew/opt/libomp/include)
link_directories(/opt/homebrew/opt/libomp/lib)

# Include directory
include_directories(include)

add_executable(AuxDM main.cpp
        src/SlaterBasisData.cpp
        include/SlaterBasisData.h
        utils/readQuadWeightsFile.cpp
        include/readQuadWeightsFile.h
        src/SlaterPrimitive.cpp
        include/SlaterPrimitive.h
        src/AtomInfo.cpp
        include/AtomInfo.h
        src/SlaterBasisSet.cpp
        include/SlaterBasisSet.h
        src/SphericalHarmonicFunc.cpp
        include/SphericalHarmonicFunc.h
        src/AuxDensityMatrix.cpp
        include/AuxDensityMatrix.h
        src/AuxDensityMatrixSlater.cpp
        include/AuxDensityMatrixSlater.h
)

# Directly link against libomp.dylib and LibTorch
target_link_libraries(AuxDM "${TORCH_LIBRARIES}" "${LIBOMP_PREFIX}/lib/libomp.dylib" "${ACCELERATE_FRAMEWORK}")

# Add OpenMP flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lomp")

# Ensure the C++17 standard is used
set_property(TARGET AuxDM PROPERTY CXX_STANDARD 17)
