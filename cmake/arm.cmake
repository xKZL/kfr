# For internal use only

set (CMAKE_SYSTEM_NAME Linux)
set (CMAKE_SYSTEM_VERSION 1)
set (UNIX True)
set (ARM True)
set (CMAKE_SYSTEM_PROCESSOR arm)
set (EMULATOR qemu-arm)

include (CMakeForceCompiler)
CMAKE_FORCE_CXX_COMPILER (/usr/bin/clang++${CLANG_SUFFIX} Clang)
CMAKE_FORCE_C_COMPILER (/usr/bin/clang${CLANG_SUFFIX} Clang)
set (CMAKE_CXX_COMPILER_WORKS TRUE)
set (CMAKE_C_COMPILER_WORKS TRUE)

set(TGT_TRIPLET arm-linux-gnueabihf)

set (ARM_ROOT "/usr/${TGT_TRIPLET}/include")
if (NOT GCC_VER)
    set (GCC_VER 7.5.0)
endif ()
set (SYS_PATHS "-isystem ${ARM_ROOT}/c++/${GCC_VER} -isystem ${ARM_ROOT}/c++/${GCC_VER}/backward -isystem ${ARM_ROOT}/c++/${GCC_VER}/${TGT_TRIPLET} -isystem ${ARM_ROOT}")

set (ARM_COMMON_FLAGS "-target ${TGT_TRIPLET} -mcpu=cortex-a15 -mfpu=neon-vfpv4 -mfloat-abi=hard -static")

set (CMAKE_CXX_FLAGS "-std=gnu++17 ${SYS_PATHS} ${ARM_COMMON_FLAGS}" CACHE STRING "")
set (CMAKE_C_FLAGS " ${SYS_PATHS} ${ARM_COMMON_FLAGS}" CACHE STRING "")

set (CMAKE_CXX_LINK_FLAGS " ${ARM_COMMON_FLAGS} ${CMAKE_CXX_LINK_FLAGS}")
set (CMAKE_C_LINK_FLAGS " ${ARM_COMMON_FLAGS} ${CMAKE_C_LINK_FLAGS}")

message(STATUS "${ARM_COMMON_FLAGS}")

set (CMAKE_FIND_ROOT_PATH_MODE_PROGRAM ONLY)
set (CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set (CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
