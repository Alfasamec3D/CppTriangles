CMake_minimum_required(VERSION 3.13)

project(triangle)

enable_testing()

add_executable(main main.cpp)

add_library(geom libraries/geom.cpp)

target_include_directories(geom PUBLIC "headers")

target_link_libraries(main geom)



add_executable(test_comp_ps Unit_Test/test_comp_ps.cpp)

add_executable(test_comp_pt Unit_Test/test_comp_pt.cpp)

add_executable(test_comp_ss Unit_Test/test_comp_ss.cpp)

add_executable(test_comp_st Unit_Test/test_comp_st.cpp)

add_executable(test_comp_tt Unit_Test/test_comp_tt.cpp)

add_executable(test_intsec Unit_Test/test_intsec.cpp)

add_executable(test_scalprod Unit_Test/test_scalprod.cpp)

add_executable(test_sqrdist Unit_Test/test_sqrdist.cpp)

add_executable(test_tr_segm Unit_Test/test_tr_segm.cpp)

add_executable(test_vectprod Unit_Test/test_vectprod.cpp)



target_link_libraries(test_comp_ps geom)

target_link_libraries(test_comp_pt geom)

target_link_libraries(test_comp_ss geom)

target_link_libraries(test_comp_st geom)

target_link_libraries(test_comp_tt geom)

target_link_libraries(test_intsec geom)

target_link_libraries(test_scalprod geom)

target_link_libraries(test_sqrdist geom)

target_link_libraries(test_tr_segm geom)

target_link_libraries(test_vectprod geom)


add_test(NAME Test_comp_ps COMMAND test_comp_ps)

add_test(NAME Test_comp_pt COMMAND test_comp_pt)

add_test(NAME Test_comp_ss COMMAND test_comp_ss)

add_test(NAME Test_comp_st COMMAND test_comp_st)

add_test(NAME Test_comp_tt COMMAND test_comp_tt)

add_test(NAME Test_intsec COMMAND test_intsec)

add_test(NAME Test_scalprod COMMAND test_scalprod)

add_test(NAME Test_sqrdist COMMAND test_sqrdist)

add_test(NAME Test_tr_segm COMMAND test_tr_segm)

add_test(NAME Test_vectprod COMMAND test_vectprod)



file(GLOB testfiles "e2e_tests/*.dat")
foreach(file ${testfiles})
  add_test(
    NAME ${file}
    COMMAND bash -c "${CMAKE_SOURCE_DIR}/runtest.sh ${file} ./main"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(${file} PROPERTIES DEPENDS main)
endforeach()
