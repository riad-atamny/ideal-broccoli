cmake_minimum_required(VERSION 3.15)
project(ideal-broccoli)

set(CMAKE_CXX_STANDARD 11)

add_executable(ideal-broccoli main.cpp switch.h rule.h packet.h Controller.h)
#include $DIr/vendor/Makefile