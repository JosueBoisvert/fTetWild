cmake_minimum_required(VERSION 3.8)
project(FloatTetwild)

# FloatTetWild library
add_library(ftetwild src/Logger.cpp src/external/WindingNumber.h)

# FloatTetWild CLI executable
set(ftetwild_cli_INCLUDE_DIRS
        3rdparty/cli11/include)

add_executable(ftetwild_cli src/main.cpp)
target_include_directories(ftetwild_cli PUBLIC ${ftetwild_cli_INCLUDE_DIRS})