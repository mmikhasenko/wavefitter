
# only one program
add_executable(add_weight_to_complete_tree add_weight_to_complete_tree.cc)
target_link_libraries(add_weight_to_complete_tree polesearch )
target_link_libraries(add_weight_to_complete_tree ${ROOT_LIBRARIES})

# only one program
add_executable(HitAndMiss HitAndMiss.cc)
target_link_libraries(HitAndMiss polesearch )
target_link_libraries(HitAndMiss ${ROOT_LIBRARIES})

# only one program
add_executable(add_weight_MDeck add_weight_MDeck.cc)
target_link_libraries(add_weight_MDeck isobars_and_deck )
target_link_libraries(add_weight_MDeck ${ROOT_LIBRARIES})
