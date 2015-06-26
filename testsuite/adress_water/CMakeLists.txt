if (NOT ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
  add_custom_target(adress_water_testdata ALL)
  foreach(_file adress.py expected.txt adress.espressopp) 
    add_custom_command(TARGET adress_water_testdata COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file})
  endforeach()
endif()
add_test(adress_water ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/adress.py)
add_test(adress_water_compare ${CMAKE_COMMAND} 
    -E compare_files ${CMAKE_CURRENT_SOURCE_DIR}/expected.txt ${CMAKE_CURRENT_SOURCE_DIR}/output.txt)
set_tests_properties(adress_water PROPERTIES ENVIRONMENT "${TEST_ENV}")
