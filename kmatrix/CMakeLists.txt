
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
add_executable(wavefitter wavefitter.cc)
target_link_libraries(wavefitter polesearch )
target_link_libraries(wavefitter ${ROOT_LIBRARIES})
target_link_libraries(wavefitter ${Boost_LIBRARIES})
target_link_libraries(wavefitter ${Libconfig_LIBS})

include_directories(${Libconfig_INCLUDE_DIR})


# for all files test_*cc do the same
add_executable(testi_MDeck testi_MDeck.cc)
target_link_libraries(testi_MDeck isobars_and_deck)
target_link_libraries(testi_MDeck ${ROOT_LIBRARIES})


# additional test with libconfig
target_link_libraries(test_libconfig ${Libconfig_LIBS})
