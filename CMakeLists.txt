cmake_minimum_required(VERSION 2.6)
project(raptordebugtest)

add_executable(raptordebugtest main.cpp)

install(TARGETS raptordebugtest RUNTIME DESTINATION bin)
