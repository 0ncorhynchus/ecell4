set(TEST_NAMES
    ODESimulator_test)

set(test_library_dependencies)
if (Boost_UNIT_TEST_FRAMEWORK_FOUND)
    add_definitions(-DBOOST_TEST_DYN_LINK)
    add_definitions(-DUNITTEST_FRAMEWORK_LIBRARY_EXIST)
    set(test_library_dependencies ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
endif()

foreach(TEST_NAME ${TEST_NAMES})
    add_executable(${TEST_NAME} ${TEST_NAME}.cpp)
    target_link_libraries(${TEST_NAME} ecell4-ode ${test_library_dependencies})
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach(TEST_NAME)

