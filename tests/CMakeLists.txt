set(LQD_TEST_SOURCES
  test_main.cpp
  test_foundations.cpp
  test_scenes.cpp
)
add_executable(lqd_tests ${LQD_TEST_SOURCES})
target_link_libraries(lqd_tests PRIVATE lqd_core)
target_compile_options(lqd_tests PRIVATE ${LQD_WARNINGS})

# The runner prints per-case results; pass a substring to run a subset: lqd_tests "spray:"
add_test(NAME lqd_tests COMMAND lqd_tests)
