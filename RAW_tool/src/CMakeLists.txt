# Recent version of cmake
cmake_minimum_required(VERSION 3.0)

# Name of the project
project(libecho)

# compiler flags
add_definitions (-g -lm -lpthread -std=c++11)

# Variable creation for file to compile
#set(SRCS
	#pointer_management.cpp
	#scan_conversion_float.cpp
	#sig_proc_float.cpp
	#TCP_API.cpp
	#methods.cpp
	#)
	
#set(HEADERS
	#pointer_management.hpp
	#scan_conversion_float.hpp
	#sig_proc_float.hpp
	#TCP_API.hpp
	#methods.hpp
	#)

# Library creation
#add_library(echo ${SRCS} ${HEADERS})
add_library(pointer_management pointer_management.cpp pointer_management.hpp)
add_library(scan_conv scan_conversion_float.cpp scan_conversion_float.hpp)
#target_link_libraries(scan_conv pointer_management)
add_library(sig_proc sig_proc_float.cpp sig_proc_float.hpp)
add_library(TCP TCP_API.cpp TCP_API.hpp)
target_link_libraries(TCP pthread)
add_library(methods methods.cpp methods.hpp)
target_link_libraries(methods pointer_management sig_proc scan_conv)

