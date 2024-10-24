###
cmake_minimum_required (VERSION 3.5)

project(sBIF)

###Adding sources
aux_source_directory(${PROJECT_SOURCE_DIR}/src SRC_LIST)

include_directories(${PROJECT_SOURCE_DIR}/src)
add_executable(sBIF ${SRC_LIST})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

if(UNIX)
    # set(CMAKE_C_COMPILER "/opt/homebrew/opt/llvm/bin/clang")
    # set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/llvm/bin/clang++")
    add_definitions(-std=c++11)
    if(CMAKE_C_COMPILER_ID MATCHES "Clang")
        set(OpenMP_C "${CMAKE_C_COMPILER}")
        set(OpenMP_C_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
        set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
        set(OpenMP_CXX_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
        set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
    endif()
else()
    add_definitions(-D_SCL_SECURE_NO_WARNINGS
    -D_CRT_SECURE_NO_WARNINGS -D_WINSOCK_DEPRECATED_NO_WARNINGS)
endif()

###Running openmp
FIND_PACKAGE(OpenMP REQUIRED)

if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS})
    set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS})
    set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS})
endif()

find_package(PostgreSQL REQUIRED)
if(PostgreSQL_FOUND)
    include_directories(${PostgreSQL_INCLUDE_DIRS})
    target_link_libraries(sBIF ${PostgreSQL_LIBRARIES})
endif()