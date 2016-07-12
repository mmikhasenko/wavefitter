file( GLOB POLESEARCH_ALL_FILES *.cc )
add_library(polesearch ${POLESEARCH_ALL_FILES} )
target_link_libraries(polesearch MathMore)
