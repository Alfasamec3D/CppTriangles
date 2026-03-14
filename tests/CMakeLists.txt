enable_testing()

find_package(GTest REQUIRED)

include(GoogleTest)

add_executable(unit_test unit/unit_test.cpp)

target_link_libraries(unit_test PRIVATE geom GTest::gtest_main)

gtest_discover_tests(unit_test)


find_program(BASH_EXECUTABLE bash REQUIRED)
file (GLOB testfiles "e2e/io/*.dat")
foreach (file ${testfiles})
    add_test(NAME ${file} 
        #COMMAND ${BASH_EXECUTABLE} 
        #       e2e/run.sh 
        #      "${file}" "$<TARGET_FILE:main>" 
        COMMAND bash -c "${CMAKE_SOURCE_DIR}/tests/e2e/run.sh ${file} ./main"
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        
)
set_tests_properties(${file} PROPERTIES DEPENDS main)
endforeach()