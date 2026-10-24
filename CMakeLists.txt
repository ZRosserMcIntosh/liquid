cmake_minimum_required(VERSION 3.20)
project(LiquidDynamics VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(LQD_BUILD_TESTS "Build Liquid Dynamics tests" ON)
option(LQD_BUILD_TOOLS "Build the lqdemo tool" ON)
option(LQD_WARNINGS_AS_ERRORS "Treat warnings as errors" OFF)

if(MSVC)
  set(LQD_WARNINGS /W4 /permissive-)
  if(LQD_WARNINGS_AS_ERRORS)
    list(APPEND LQD_WARNINGS /WX)
  endif()
else()
  set(LQD_WARNINGS -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wno-sign-conversion)
  if(LQD_WARNINGS_AS_ERRORS)
    list(APPEND LQD_WARNINGS -Werror)
  endif()
endif()

# ---- Engine-independent fluid core (no engine/renderer headers) ----
add_library(lqd_core
  src/core/src/json.cpp
  src/core/src/library.cpp
  src/core/src/grid.cpp
  src/core/src/vessel.cpp
  src/core/src/collider.cpp
  src/core/src/rigid.cpp
  src/core/src/ragdoll.cpp
  src/core/src/fluid.cpp
  src/core/src/spray.cpp
  src/core/src/bpa.cpp
  src/core/src/world.cpp
  src/core/src/scenarios.cpp
  src/core/src/mesher.cpp
  src/core/src/export.cpp
)
add_library(lqd::core ALIAS lqd_core)
target_include_directories(lqd_core PUBLIC src/core/include)
target_compile_options(lqd_core PRIVATE ${LQD_WARNINGS})
target_compile_definitions(lqd_core PUBLIC LQD_DATA_DIR="${PROJECT_SOURCE_DIR}/data")

if(LQD_BUILD_TOOLS)
  add_executable(lqdemo tools/lqdemo/main.cpp)
  target_link_libraries(lqdemo PRIVATE lqd_core)
  target_compile_options(lqdemo PRIVATE ${LQD_WARNINGS})
endif()

if(LQD_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
