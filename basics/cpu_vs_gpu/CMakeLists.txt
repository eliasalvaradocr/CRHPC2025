cmake_minimum_required(VERSION 3.12)

project(cpu_vs_gpu LANGUAGES CUDA C)

set(CMAKE_C_STANDARD 99)

add_executable(cpu_version_block main_block.c matrix.c)
add_executable(cpu_version_naive main_naive.c matrix.c)
add_executable(gpu_version_naive main.cu matrix.c)
add_executable(gpu_version_coalesced main_coalesced.cu matrix.c)
