
# only one program
add_executable(add_weight_to_tree add_weight_to_tree.cc)
target_link_libraries(add_weight_to_tree polesearch )
target_link_libraries(add_weight_to_tree ${ROOT_LIBRARIES})
