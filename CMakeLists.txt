cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_STANDARD 20)
project("MeshHoleFilling" C CXX)

add_compile_definitions(NOMINMAX)

find_package(Eigen3 REQUIRED)
find_package(CGAL REQUIRED)
find_package(OpenMP REQUIRED)
include(CGAL_Eigen3_support)

add_library(MeshHoleFilling INTERFACE "MeshHoleFilling.h")
target_link_libraries(MeshHoleFilling INTERFACE CGAL::CGAL CGAL::Eigen3_support OpenMP::OpenMP_CXX)
target_include_directories(MeshHoleFilling INTERFACE .)

add_executable(MeshHoleFillingTest "test.cpp")
target_link_libraries(MeshHoleFillingTest PUBLIC MeshHoleFilling)

