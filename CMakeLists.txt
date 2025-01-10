cmake_minimum_required(VERSION 3.10)

project(HestonModel)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(HestonModelApp main.cpp heston.cpp)