cmake_minimum_required(VERSION 2.6)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

#set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/.. )

project (s-rlbwt)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/internal)
include_directories(~/include) #SDSL headers are here

LINK_DIRECTORIES(~/lib) #SDSL lib are here

message("Building in ${CMAKE_BUILD_TYPE} mode")


if(CPP_STATIC)
set(CMAKE_CXX_FLAGS "-static")
endif(CPP_STATIC)


set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g")
set(CMAKE_CXX_FLAGS_RELEASE "-g -ggdb -Ofast -fstrict-aliasing -DNDEBUG -march=native")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -ggdb -Ofast -fstrict-aliasing -march=native")

add_executable(mod_ri-build mod_ri-build.cpp)
TARGET_LINK_LIBRARIES(mod_ri-build sdsl)
TARGET_LINK_LIBRARIES(mod_ri-build divsufsort)
TARGET_LINK_LIBRARIES(mod_ri-build divsufsort64)


add_executable(ri-build ri-build.cpp)
TARGET_LINK_LIBRARIES(ri-build sdsl)
TARGET_LINK_LIBRARIES(ri-build divsufsort)
TARGET_LINK_LIBRARIES(ri-build divsufsort64)

add_executable(ri-locate ri-locate.cpp)
TARGET_LINK_LIBRARIES(ri-locate sdsl)
TARGET_LINK_LIBRARIES(ri-locate divsufsort)
TARGET_LINK_LIBRARIES(ri-locate divsufsort64)

add_executable(ri-count ri-count.cpp)
TARGET_LINK_LIBRARIES(ri-count sdsl)
TARGET_LINK_LIBRARIES(ri-count divsufsort)
TARGET_LINK_LIBRARIES(ri-count divsufsort64)

#add_executable(ri-space ri-space.cpp)
#TARGET_LINK_LIBRARIES(ri-space sdsl)
#TARGET_LINK_LIBRARIES(ri-space divsufsort)
#TARGET_LINK_LIBRARIES(ri-space divsufsort64)
