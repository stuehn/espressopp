add_test(rattle ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_rattle.py)
set_tests_properties(rattle PROPERTIES ENVIRONMENT "${TEST_ENV}")
