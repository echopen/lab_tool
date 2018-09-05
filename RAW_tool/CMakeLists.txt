# Minimal version of cmake used to compile
cmake_minimum_required(VERSION 3.0)

# Declaration of the project
project (RAW_tool)

# Compiler option
add_definitions (-g -lpthread -lm )

# Include directory for sources
include_directories(src)

# Define used directories first source, second destination
add_subdirectory(src)

# Notre exécutable
add_executable(RAW_tool RAW_tool.cpp)


# Et pour que l'exécutable fonctionne, il faut lui indiquer la bibliothèque dont il dépend
target_link_libraries(RAW_tool pointer_management scan_conv sig_proc TCP methods)
