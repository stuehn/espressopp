add_test(adress_water ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/adress.py)
add_test(adress_water_compare ${CMAKE_COMMAND} 
    -E compare_files ${CMAKE_CURRENT_SOURCE_DIR}/expected.txt ${CMAKE_CURRENT_SOURCE_DIR}/output.txt)
set_tests_properties(adress_water PROPERTIES ENVIRONMENT "${TEST_ENV}")
