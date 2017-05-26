
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

# only one program
add_executable(AngularIntegral AngularIntegral.cc)
target_link_libraries(AngularIntegral polesearch )
target_link_libraries(AngularIntegral isobars_and_deck )
target_link_libraries(AngularIntegral ${ROOT_LIBRARIES})

# only one program
add_executable(SymmetrizedPhaseSpace SymmetrizedPhaseSpace.cc)
target_link_libraries(SymmetrizedPhaseSpace polesearch )
target_link_libraries(SymmetrizedPhaseSpace isobars_and_deck )
target_link_libraries(SymmetrizedPhaseSpace ${ROOT_LIBRARIES})

# only one program
add_executable(calculate_phase_space calculate_phase_space.cc)
target_link_libraries(calculate_phase_space polesearch )
target_link_libraries(calculate_phase_space isobars_and_deck )
target_link_libraries(calculate_phase_space ${ROOT_LIBRARIES})

# only one program
add_executable(generate_3pi_angles generate_3pi_angles.cc)
target_link_libraries(generate_3pi_angles polesearch )
target_link_libraries(generate_3pi_angles ${ROOT_LIBRARIES})

# only one program
add_executable(calculate_Deck_projections calculate_Deck_projections.cc)
target_link_libraries(calculate_Deck_projections polesearch )
target_link_libraries(calculate_Deck_projections ${ROOT_LIBRARIES})

# only one program
add_executable(invert_matrix invert_matrix.cc)
target_link_libraries(invert_matrix polesearch )
target_link_libraries(invert_matrix ${ROOT_LIBRARIES})

# only one program
add_executable(add_weight_PWs add_weight_PWs.cc)
target_link_libraries(add_weight_PWs polesearch )
target_link_libraries(add_weight_PWs ${ROOT_LIBRARIES})

# only one program
add_executable(cross_check_deck_shape cross_check_deck_shape.cc)
target_link_libraries(cross_check_deck_shape polesearch )
target_link_libraries(cross_check_deck_shape ${ROOT_LIBRARIES})

# only one program
add_executable(calculate_analytical_projections calculate_analytical_projections.cc)
target_link_libraries(calculate_analytical_projections polesearch )
target_link_libraries(calculate_analytical_projections ${ROOT_LIBRARIES})

# only one program
add_executable(calculate_analytical_exact_2d_projections calculate_analytical_exact_2d_projections.cc)
target_link_libraries(calculate_analytical_exact_2d_projections polesearch )
target_link_libraries(calculate_analytical_exact_2d_projections ${ROOT_LIBRARIES})

# only one program
add_executable(calculate_analytical_exact_2d_projections_summed_over_helicity calculate_analytical_exact_2d_projections_summed_over_helicity.cc)
target_link_libraries(calculate_analytical_exact_2d_projections_summed_over_helicity polesearch )
target_link_libraries(calculate_analytical_exact_2d_projections_summed_over_helicity ${ROOT_LIBRARIES})
