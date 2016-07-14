
# for all files test_*cc do the same
file (GLOB APP_SOURCES test_*.cc)
foreach (sourcefilepath ${APP_SOURCES})
    # make name from path by taking name and cutting the extention
    get_filename_component(sourcename ${sourcefilepath} NAME)
    string( REPLACE ".cc" "" appname ${sourcename} )
    # add executable and dependencies
    add_executable(${appname} ${sourcefilepath} )
    target_link_libraries( ${appname} polesearch )
    target_link_libraries( ${appname} ${ROOT_LIBRARIES})
endforeach (sourcefilepath ${APP_SOURCES})

# tests which need ublas
target_link_libraries(test_ublas_prod ${Boost_LIBRARIES})
target_link_libraries(test_kmatrix ${Boost_LIBRARIES})


# main program
# include_directories(${Libconfig_INCLUDE_DIR})
add_executable(read_model_settings read_model_settings.cc)
target_link_libraries(read_model_settings polesearch )
target_link_libraries(read_model_settings ${ROOT_LIBRARIES})
target_link_libraries(read_model_settings ${Boost_LIBRARIES})
target_link_libraries(read_model_settings config++)

# target_link_libraries(read_model_settings ${Boost_LIBRARIES})
