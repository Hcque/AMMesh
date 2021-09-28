cmake_minimum_required(VERSION 3.5)

project(ACAM_mesh)

#设置不生成MinSizeRel;RelWithDebInfo
set(CMAKE_CONFIGURATION_TYPES "Debug;Release")

#path
SET(EIGEN_PATH "~/Desktop/eigen/Eigen")
set(CMAKE_BUILD_TYPE Debug)
############################################################
# Create a library
############################################################

#Generate the static library from the library sources
add_library(ACAM_mesh STATIC 
    src/IOManger.cpp
	src/PolyMesh.cpp
	src/PolyMesh_Base.cpp
	src/PolyMeshIterators.cpp
)

#SET(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib) #设置库文件的输出目录

target_include_directories(ACAM_mesh
    PUBLIC 
        ${PROJECT_SOURCE_DIR}/include
)


############################################################
# Create an executable
############################################################

# Add an executable with the above sources
add_executable(ACAM_example 
    src/main.cpp
)

# link the new ACAM_example target with the ACAM_mesh target
target_link_libraries( ACAM_example
    PRIVATE 
        ACAM_mesh
)

############################################################
# Project Homework 1
############################################################

file (GLOB HW1_FILE
    src/hw1/*.cpp
)

# Add an executable with the above sources
add_executable(${PROJECT_NAME}_HW1 ${HW1_FILE})
target_include_directories(${PROJECT_NAME}_HW1
    PRIVATE 
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(${PROJECT_NAME}_HW1
    PRIVATE 
        ACAM_mesh
)


############################################################
# Project Homework 2
############################################################

file (GLOB HW2_FILE
    src/hw2/*.cpp
)

# Add an executable with the above sources
add_executable(${PROJECT_NAME}_HW2 ${HW2_FILE})
target_include_directories(${PROJECT_NAME}_HW2
    PRIVATE 
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(${PROJECT_NAME}_HW2
    PRIVATE 
        ACAM_mesh
)


############################################################
# Project Homework 3
############################################################

file (GLOB HW3_FILE
    src/hw3/*.cpp
)

# Add an executable with the above sources
add_executable(${PROJECT_NAME}_HW3 ${HW3_FILE})
target_include_directories(${PROJECT_NAME}_HW3
    PRIVATE 
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(${PROJECT_NAME}_HW3
    PRIVATE 
        ACAM_mesh
)

############################################################
# Project Homework 9
############################################################

file (GLOB HW9_FILE
    src/hw9/*.cpp
)

# Add an executable with the above sources
add_executable(${PROJECT_NAME}_HW9 ${HW9_FILE})
target_include_directories(${PROJECT_NAME}_HW9
    PRIVATE 
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(${PROJECT_NAME}_HW9
    PRIVATE 
        ACAM_mesh
)

