# all files
file( GLOB POLESEARCH_ALL_FILES *.cc )
add_library(polesearch ${POLESEARCH_ALL_FILES} )
target_link_libraries(polesearch MathMore)

# isobars and deck
add_library(isobars_and_deck MIsobar.cc MDeck.cc waves.cc MAscoli.cc deflib.cc mintegrate.cc TWigner.cc)
target_link_libraries(isobars_and_deck MathMore)
