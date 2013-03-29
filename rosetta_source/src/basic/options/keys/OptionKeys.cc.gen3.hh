namespace edensity { BooleanOptionKey const score_sliding_window_context( "edensity:score_sliding_window_context" );  }
namespace edensity { RealOptionKey const whole_structure_ca_wt( "edensity:whole_structure_ca_wt" );  }
namespace edensity { RealOptionKey const whole_structure_allatom_wt( "edensity:whole_structure_allatom_wt" );  }
namespace edensity { BooleanOptionKey const no_edens_in_minimizer( "edensity:no_edens_in_minimizer" );  }
namespace edensity { BooleanOptionKey const debug_derivatives( "edensity:debug_derivatives" );  }
namespace edensity { StringOptionKey const realign( "edensity:realign" );  }
namespace edensity { StringOptionKey const membrane_axis( "edensity:membrane_axis" );  }
namespace edensity { RealOptionKey const atom_mask( "edensity:atom_mask" );  }
namespace edensity { RealOptionKey const ca_mask( "edensity:ca_mask" );  }
namespace edensity { BooleanOptionKey const score_symm_complex( "edensity:score_symm_complex" );  }
namespace edensity { RealOptionKey const sc_scaling( "edensity:sc_scaling" );  }
namespace edensity { RealOptionKey const render_sigma( "edensity:render_sigma" );  }
namespace patterson { BooleanOptionKey const patterson( "patterson" );  }
namespace patterson { BooleanOptionKey const debug( "patterson:debug" );  }
namespace patterson { RealOptionKey const weight( "patterson:weight" );  }
namespace patterson { RealOptionKey const sc_scaling( "patterson:sc_scaling" );  }
namespace patterson { RealVectorOptionKey const radius_cutoffs( "patterson:radius_cutoffs" );  }
namespace patterson { RealVectorOptionKey const resolution_cutoffs( "patterson:resolution_cutoffs" );  }
namespace patterson { RealOptionKey const Bsol( "patterson:Bsol" );  }
namespace patterson { RealOptionKey const Fsol( "patterson:Fsol" );  }
namespace patterson { RealOptionKey const model_B( "patterson:model_B" );  }
namespace patterson { RealOptionKey const rmsd( "patterson:rmsd" );  }
namespace patterson { BooleanOptionKey const no_ecalc( "patterson:no_ecalc" );  }
namespace patterson { IntegerOptionKey const nshells( "patterson:nshells" );  }
namespace patterson { BooleanOptionKey const use_spline_interpolation( "patterson:use_spline_interpolation" );  }
namespace patterson { BooleanOptionKey const use_on_repack( "patterson:use_on_repack" );  }
namespace patterson { BooleanOptionKey const dont_use_symm_in_pcalc( "patterson:dont_use_symm_in_pcalc" );  }
namespace cryst { BooleanOptionKey const cryst( "cryst" );  }
namespace cryst { StringOptionKey const mtzfile( "cryst:mtzfile" );  }
namespace cryst { BooleanOptionKey const crystal_refine( "cryst:crystal_refine" );  }
namespace optE { BooleanOptionKey const optE( "optE" );  }
namespace optE { StringOptionKey const load_from_silent( "optE:load_from_silent" );  }
namespace optE { StringOptionKey const data_in( "optE:data_in" );  }
namespace optE { StringOptionKey const data_out( "optE:data_out" );  }
namespace optE { StringOptionKey const weights( "optE:weights" );  }
namespace optE { StringVectorOptionKey const fix( "optE:fix" );  }
namespace optE { FileOptionKey const free( "optE:free" );  }
namespace optE { FileOptionKey const fixed( "optE:fixed" );  }
namespace optE { FileOptionKey const parse_tagfile( "optE:parse_tagfile" );  }
namespace optE { FileOptionKey const constant_logic_taskops_file( "optE:constant_logic_taskops_file" );  }
namespace optE { BooleanOptionKey const optE_soft_rep( "optE:optE_soft_rep" );  }
namespace optE { BooleanOptionKey const no_hb_env_dependence( "optE:no_hb_env_dependence" );  }
namespace optE { BooleanOptionKey const no_hb_env_dependence_DNA( "optE:no_hb_env_dependence_DNA" );  }
namespace optE { BooleanOptionKey const optE_no_protein_hack_elec( "optE:optE_no_protein_hack_elec" );  }
namespace optE { BooleanOptionKey const design_first( "optE:design_first" );  }
namespace optE { IntegerOptionKey const n_design_cycles( "optE:n_design_cycles" );  }
namespace optE { BooleanOptionKey const recover_nat_rot( "optE:recover_nat_rot" );  }
namespace optE { FileOptionKey const component_weights( "optE:component_weights" );  }
namespace optE { BooleanOptionKey const optimize_nat_aa( "optE:optimize_nat_aa" );  }
namespace optE { BooleanOptionKey const optimize_nat_rot( "optE:optimize_nat_rot" );  }
namespace optE { FileOptionKey const optimize_ligand_rot( "optE:optimize_ligand_rot" );  }
namespace optE { BooleanOptionKey const optimize_pssm( "optE:optimize_pssm" );  }
namespace optE { FileOptionKey const optimize_dGbinding( "optE:optimize_dGbinding" );  }
namespace optE { FileOptionKey const optimize_ddG_bind_correlation( "optE:optimize_ddG_bind_correlation" );  }
namespace optE { FileOptionKey const optimize_ddGmutation( "optE:optimize_ddGmutation" );  }
namespace optE { BooleanOptionKey const optimize_ddGmutation_straight_mean( "optE:optimize_ddGmutation_straight_mean" );  }
namespace optE { BooleanOptionKey const optimize_ddGmutation_boltzman_average( "optE:optimize_ddGmutation_boltzman_average" );  }
namespace optE { RealOptionKey const exclude_badrep_ddGs( "optE:exclude_badrep_ddGs" );  }
namespace optE { BooleanOptionKey const pretend_no_ddG_repulsion( "optE:pretend_no_ddG_repulsion" );  }
namespace optE { FileOptionKey const optimize_decoy_discrimination( "optE:optimize_decoy_discrimination" );  }
namespace optE { StringOptionKey const normalize_decoy_score_spread( "optE:normalize_decoy_score_spread" );  }
namespace optE { BooleanOptionKey const ramp_nativeness( "optE:ramp_nativeness" );  }
namespace optE { IntegerOptionKey const n_top_natives_to_optimize( "optE:n_top_natives_to_optimize" );  }
namespace optE { RealOptionKey const approximate_decoy_entropy( "optE:approximate_decoy_entropy" );  }
namespace optE { BooleanOptionKey const repack_and_minimize_decoys( "optE:repack_and_minimize_decoys" );  }
namespace optE { BooleanOptionKey const repack_and_minimize_input_structures( "optE:repack_and_minimize_input_structures" );  }
namespace optE { IntegerOptionKey const output_top_n_new_decoys( "optE:output_top_n_new_decoys" );  }
namespace optE { FileOptionKey const optimize_ligand_discrimination( "optE:optimize_ligand_discrimination" );  }
namespace optE { BooleanOptionKey const no_design( "optE:no_design" );  }
namespace optE { BooleanOptionKey const sqrt_pssm( "optE:sqrt_pssm" );  }
namespace optE { RealOptionKey const min_decoy_rms_to_native( "optE:min_decoy_rms_to_native" );  }
namespace optE { RealOptionKey const max_rms_from_native( "optE:max_rms_from_native" );  }
namespace optE { BooleanOptionKey const optimize_starting_free_weights( "optE:optimize_starting_free_weights" );  }
namespace optE { FileOptionKey const wrap_dof_optimization( "optE:wrap_dof_optimization" );  }
namespace optE { RealOptionKey const randomly_perturb_starting_free_weights( "optE:randomly_perturb_starting_free_weights" );  }
namespace optE { RealOptionKey const inv_kT_natrot( "optE:inv_kT_natrot" );  }
namespace optE { RealOptionKey const inv_kT_nataa( "optE:inv_kT_nataa" );  }
namespace optE { RealOptionKey const inv_kT_natstruct( "optE:inv_kT_natstruct" );  }
namespace optE { BooleanOptionKey const mpi_weight_minimization( "optE:mpi_weight_minimization" );  }
namespace optE { BooleanOptionKey const dont_use_reference_energies( "optE:dont_use_reference_energies" );  }
namespace optE { IntegerOptionKey const number_of_swarm_particles( "optE:number_of_swarm_particles" );  }
namespace optE { IntegerOptionKey const number_of_swarm_cycles( "optE:number_of_swarm_cycles" );  }
namespace optE { FileOptionKey const constrain_weights( "optE:constrain_weights" );  }
namespace optE { BooleanOptionKey const fit_reference_energies_to_aa_profile_recovery( "optE:fit_reference_energies_to_aa_profile_recovery" );  }
namespace optE { FileOptionKey const starting_refEs( "optE:starting_refEs" );  }
namespace optE { BooleanOptionKey const repeat_swarm_optimization_until_fitness_improves( "optE:repeat_swarm_optimization_until_fitness_improves" );  }
namespace optE { BooleanOptionKey const design_with_minpack( "optE:design_with_minpack" );  }
namespace optE { BooleanOptionKey const limit_bad_scores( "optE:limit_bad_scores" );  }
namespace optE { namespace rescore { BooleanOptionKey const rescore( "optE:rescore" );  } }
namespace optE { namespace rescore { FileOptionKey const weights( "optE:rescore:weights" );  } }
namespace optE { namespace rescore { IntegerOptionKey const context_round( "optE:rescore:context_round" );  } }
namespace optE { namespace rescore { FileOptionKey const outlog( "optE:rescore:outlog" );  } }
namespace optE { namespace rescore { BooleanOptionKey const measure_sequence_recovery( "optE:rescore:measure_sequence_recovery" );  } }
namespace optE { BooleanOptionKey const no_design_pdb_output( "optE:no_design_pdb_output" );  }
namespace backrub { BooleanOptionKey const backrub( "backrub" );  }
namespace backrub { IntegerVectorOptionKey const pivot_residues( "backrub:pivot_residues" );  }
namespace backrub { StringVectorOptionKey const pivot_atoms( "backrub:pivot_atoms" );  }
namespace backrub { IntegerOptionKey const min_atoms( "backrub:min_atoms" );  }
namespace backrub { IntegerOptionKey const max_atoms( "backrub:max_atoms" );  }
namespace bbg { BooleanOptionKey const bbg( "bbg" );  }
namespace bbg { RealOptionKey const factorA( "bbg:factorA" );  }
namespace bbg { RealOptionKey const factorB( "bbg:factorB" );  }
namespace bbg { BooleanOptionKey const ignore_improper_res( "bbg:ignore_improper_res" );  }
namespace bbg { BooleanOptionKey const fix_short_segment( "bbg:fix_short_segment" );  }
namespace flexpack { BooleanOptionKey const flexpack( "flexpack" );  }
namespace flexpack { namespace annealer { BooleanOptionKey const annealer( "flexpack:annealer" );  } }
namespace flexpack { namespace annealer { RealOptionKey const inner_iteration_scale( "flexpack:annealer:inner_iteration_scale" );  } }
namespace flexpack { namespace annealer { RealOptionKey const outer_iteration_scale( "flexpack:annealer:outer_iteration_scale" );  } }
namespace flexpack { namespace annealer { RealOptionKey const fixbb_substitutions_scale( "flexpack:annealer:fixbb_substitutions_scale" );  } }
namespace flexpack { namespace annealer { RealOptionKey const pure_movebb_substitutions_scale( "flexpack:annealer:pure_movebb_substitutions_scale" );  } }
namespace flexpack { namespace annealer { RealOptionKey const rotsub_movebb_substitutions_scale( "flexpack:annealer:rotsub_movebb_substitutions_scale" );  } }
namespace hotspot { BooleanOptionKey const hotspot( "hotspot" );  }
namespace hotspot { BooleanOptionKey const allow_gly( "hotspot:allow_gly" );  }
namespace hotspot { BooleanOptionKey const allow_proline( "hotspot:allow_proline" );  }
namespace hotspot { BooleanOptionKey const benchmark( "hotspot:benchmark" );  }
namespace hotspot { StringVectorOptionKey const residue( "hotspot:residue" );  }
namespace hotspot { FileOptionKey const hashfile( "hotspot:hashfile" );  }
namespace hotspot { FileOptionKey const target( "hotspot:target" );  }
namespace hotspot { IntegerOptionKey const target_res( "hotspot:target_res" );  }
namespace hotspot { RealOptionKey const target_dist( "hotspot:target_dist" );  }
namespace hotspot { FileOptionKey const density( "hotspot:density" );  }
namespace hotspot { FileOptionKey const weighted_density( "hotspot:weighted_density" );  }
namespace hotspot { FileOptionKey const rms_target( "hotspot:rms_target" );  }
namespace hotspot { FileOptionKey const rms_hotspot( "hotspot:rms_hotspot" );  }
namespace hotspot { IntegerOptionKey const rms_hotspot_res( "hotspot:rms_hotspot_res" );  }
namespace hotspot { BooleanOptionKey const rescore( "hotspot:rescore" );  }
namespace hotspot { RealOptionKey const threshold( "hotspot:threshold" );  }
namespace hotspot { BooleanOptionKey const sc_only( "hotspot:sc_only" );  }
namespace hotspot { BooleanOptionKey const fxnal_group( "hotspot:fxnal_group" );  }
namespace hotspot { BooleanOptionKey const cluster( "hotspot:cluster" );  }
namespace hotspot { BooleanOptionKey const colonyE( "hotspot:colonyE" );  }
namespace hotspot { IntegerOptionKey const length( "hotspot:length" );  }
namespace hotspot { BooleanOptionKey const envhb( "hotspot:envhb" );  }
namespace hotspot { RealOptionKey const angle( "hotspot:angle" );  }
namespace hotspot { IntegerOptionKey const angle_res( "hotspot:angle_res" );  }
namespace parser { BooleanOptionKey const parser( "parser" );  }
namespace parser { StringOptionKey const protocol( "parser:protocol" );  }
namespace parser { StringVectorOptionKey const script_vars( "parser:script_vars" );  }
namespace parser { BooleanOptionKey const view( "parser:view" );  }
namespace parser { StringOptionKey const patchdock( "parser:patchdock" );  }
namespace parser { IntegerVectorOptionKey const patchdock_random_entry( "parser:patchdock_random_entry" );  }
namespace DomainAssembly { BooleanOptionKey const DomainAssembly( "DomainAssembly" );  }
namespace DomainAssembly { BooleanOptionKey const da_setup( "DomainAssembly:da_setup" );  }
namespace DomainAssembly { FileOptionKey const da_setup_option_file( "DomainAssembly:da_setup_option_file" );  }
namespace DomainAssembly { FileOptionKey const da_setup_output_pdb( "DomainAssembly:da_setup_output_pdb" );  }
namespace DomainAssembly { FileOptionKey const da_linker_file( "DomainAssembly:da_linker_file" );  }
namespace DomainAssembly { FileOptionKey const da_start_pdb( "DomainAssembly:da_start_pdb" );  }
namespace DomainAssembly { IntegerOptionKey const da_nruns( "DomainAssembly:da_nruns" );  }
namespace DomainAssembly { IntegerOptionKey const da_start_pdb_num( "DomainAssembly:da_start_pdb_num" );  }
namespace DomainAssembly { FileOptionKey const da_linker_file_rna( "DomainAssembly:da_linker_file_rna" );  }
namespace remodel { BooleanOptionKey const remodel( "remodel" );  }
namespace remodel { BooleanOptionKey const help( "remodel:help" );  }
namespace remodel { BooleanOptionKey const autopilot( "remodel:autopilot" );  }
namespace remodel { FileOptionKey const blueprint( "remodel:blueprint" );  }
namespace remodel { FileOptionKey const cstfile( "remodel:cstfile" );  }
namespace remodel { IntegerOptionKey const cstfilter( "remodel:cstfilter" );  }
namespace remodel { StringOptionKey const cen_sfxn( "remodel:cen_sfxn" );  }
namespace remodel { BooleanOptionKey const check_scored_centroid( "remodel:check_scored_centroid" );  }
namespace remodel { IntegerOptionKey const num_trajectory( "remodel:num_trajectory" );  }
namespace remodel { IntegerOptionKey const save_top( "remodel:save_top" );  }
namespace remodel { BooleanOptionKey const swap_refine_confirm_protocols( "remodel:swap_refine_confirm_protocols" );  }
namespace remodel { IntegerOptionKey const num_frag_moves( "remodel:num_frag_moves" );  }
namespace remodel { BooleanOptionKey const bypass_fragments( "remodel:bypass_fragments" );  }
namespace remodel { BooleanOptionKey const use_same_length_fragments( "remodel:use_same_length_fragments" );  }
namespace remodel { BooleanOptionKey const enable_ligand_aa( "remodel:enable_ligand_aa" );  }
namespace remodel { BooleanOptionKey const no_jumps( "remodel:no_jumps" );  }
namespace remodel { BooleanOptionKey const backrub( "remodel:backrub" );  }
namespace remodel { BooleanOptionKey const use_blueprint_sequence ( "remodel:use_blueprint_sequence " );  }
namespace remodel { BooleanOptionKey const randomize_equivalent_fragments ( "remodel:randomize_equivalent_fragments " );  }
namespace remodel { BooleanOptionKey const quick_and_dirty ( "remodel:quick_and_dirty " );  }
namespace remodel { BooleanOptionKey const checkpoint ( "remodel:checkpoint " );  }
namespace remodel { BooleanOptionKey const use_ccd_refine ( "remodel:use_ccd_refine " );  }
namespace remodel { BooleanOptionKey const use_pose_relax ( "remodel:use_pose_relax " );  }
namespace remodel { BooleanOptionKey const use_cart_relax ( "remodel:use_cart_relax " );  }
namespace remodel { BooleanOptionKey const free_relax ( "remodel:free_relax " );  }
namespace remodel { BooleanOptionKey const use_dssp_assignment( "remodel:use_dssp_assignment" );  }
namespace remodel { BooleanOptionKey const keep_jumps_in_minimizer ( "remodel:keep_jumps_in_minimizer " );  }
namespace remodel { FileOptionKey const output_fragfiles( "remodel:output_fragfiles" );  }
namespace remodel { FileOptionKey const read_fragfile( "remodel:read_fragfile" );  }
namespace remodel { StringOptionKey const generic_aa( "remodel:generic_aa" );  }
namespace remodel { RealOptionKey const cluster_radius( "remodel:cluster_radius" );  }
namespace remodel { BooleanOptionKey const use_clusters( "remodel:use_clusters" );  }
namespace remodel { BooleanOptionKey const run_confirmation( "remodel:run_confirmation" );  }
namespace remodel { BooleanOptionKey const cluster_on_entire_pose( "remodel:cluster_on_entire_pose" );  }
namespace remodel { IntegerOptionKey const collect_clustered_top( "remodel:collect_clustered_top" );  }
namespace remodel { IntegerOptionKey const dr_cycles( "remodel:dr_cycles" );  }
namespace remodel { IntegerOptionKey const two_chain_tree( "remodel:two_chain_tree" );  }
namespace remodel { IntegerOptionKey const repeat_structure( "remodel:repeat_structure" );  }
namespace remodel { IntegerOptionKey const lh_ex_limit( "remodel:lh_ex_limit" );  }
namespace remodel { StringVectorOptionKey const lh_filter_string( "remodel:lh_filter_string" );  }
namespace remodel { IntegerOptionKey const lh_cbreak_selection( "remodel:lh_cbreak_selection" );  }
namespace remodel { BooleanOptionKey const lh_closure_filter( "remodel:lh_closure_filter" );  }
namespace remodel { BooleanOptionKey const cen_minimize( "remodel:cen_minimize" );  }
namespace remodel { IntegerOptionKey const core_cutoff( "remodel:core_cutoff" );  }
namespace remodel { IntegerOptionKey const boundary_cutoff( "remodel:boundary_cutoff" );  }
namespace remodel { BooleanOptionKey const resclass_by_sasa( "remodel:resclass_by_sasa" );  }
namespace remodel { RealOptionKey const helical_rise( "remodel:helical_rise" );  }
namespace remodel { RealOptionKey const helical_radius( "remodel:helical_radius" );  }
namespace remodel { RealOptionKey const helical_omega( "remodel:helical_omega" );  }
namespace remodel { BooleanOptionKey const abinitio_like_sampling( "remodel:abinitio_like_sampling" );  }
namespace remodel { FileOptionKey const disallow_sampling_at_pos( "remodel:disallow_sampling_at_pos" );  }
namespace remodel { BooleanOptionKey const use_fragment_sequence( "remodel:use_fragment_sequence" );  }
namespace remodel { StringOptionKey const use_sequence( "remodel:use_sequence" );  }
namespace remodel { namespace domainFusion { BooleanOptionKey const domainFusion( "remodel:domainFusion" );  } }
namespace remodel { namespace domainFusion { FileOptionKey const insert_segment_from_pdb( "remodel:domainFusion:insert_segment_from_pdb" );  } }
namespace remodel { RealOptionKey const vdw( "remodel:vdw" );  }
namespace remodel { RealOptionKey const rama( "remodel:rama" );  }
namespace remodel { RealOptionKey const cbeta( "remodel:cbeta" );  }
namespace remodel { RealOptionKey const cenpack( "remodel:cenpack" );  }
namespace remodel { RealOptionKey const hb_lrbb( "remodel:hb_lrbb" );  }
namespace remodel { RealOptionKey const hb_srbb( "remodel:hb_srbb" );  }
namespace remodel { RealOptionKey const rg( "remodel:rg" );  }
namespace remodel { RealOptionKey const rsigma( "remodel:rsigma" );  }
namespace remodel { RealOptionKey const ss_pair( "remodel:ss_pair" );  }
namespace remodel { BooleanOptionKey const build_disulf( "remodel:build_disulf" );  }
namespace remodel { IntegerOptionKey const max_disulf_allowed( "remodel:max_disulf_allowed" );  }
namespace remodel { RealOptionKey const match_rt_limit( "remodel:match_rt_limit" );  }
namespace remodel { IntegerVectorOptionKey const disulf_landing_range( "remodel:disulf_landing_range" );  }
namespace remodel { namespace design { BooleanOptionKey const design( "remodel:design" );  } }
namespace remodel { namespace design { BooleanOptionKey const no_design ( "remodel:design:no_design " );  } }
namespace remodel { namespace design { BooleanOptionKey const design_all( "remodel:design:design_all" );  } }
namespace remodel { namespace design { BooleanOptionKey const allow_rare_aro_chi( "remodel:design:allow_rare_aro_chi" );  } }
namespace remodel { namespace design { BooleanOptionKey const silent( "remodel:design:silent" );  } }
namespace remodel { namespace design { BooleanOptionKey const skip_partial( "remodel:design:skip_partial" );  } }
namespace remodel { namespace design { BooleanOptionKey const design_neighbors( "remodel:design:design_neighbors" );  } }
namespace remodel { namespace design { BooleanOptionKey const find_neighbors( "remodel:design:find_neighbors" );  } }
namespace remodel { BooleanOptionKey const rank_by_bsasa( "remodel:rank_by_bsasa" );  }
namespace remodel { namespace RemodelLoopMover { BooleanOptionKey const RemodelLoopMover( "remodel:RemodelLoopMover" );  } }
namespace remodel { namespace RemodelLoopMover { RealOptionKey const max_linear_chainbreak( "remodel:RemodelLoopMover:max_linear_chainbreak" );  } }
namespace remodel { namespace RemodelLoopMover { BooleanOptionKey const randomize_loops( "remodel:RemodelLoopMover:randomize_loops" );  } }
namespace remodel { namespace RemodelLoopMover { BooleanOptionKey const use_loop_hash( "remodel:RemodelLoopMover:use_loop_hash" );  } }
namespace remodel { namespace RemodelLoopMover { IntegerOptionKey const allowed_closure_attempts( "remodel:RemodelLoopMover:allowed_closure_attempts" );  } }
namespace remodel { namespace RemodelLoopMover { IntegerOptionKey const loophash_cycles( "remodel:RemodelLoopMover:loophash_cycles" );  } }
namespace remodel { namespace RemodelLoopMover { IntegerOptionKey const simultaneous_cycles( "remodel:RemodelLoopMover:simultaneous_cycles" );  } }
namespace remodel { namespace RemodelLoopMover { IntegerOptionKey const independent_cycles( "remodel:RemodelLoopMover:independent_cycles" );  } }
namespace remodel { namespace RemodelLoopMover { IntegerOptionKey const boost_closure_cycles( "remodel:RemodelLoopMover:boost_closure_cycles" );  } }
namespace remodel { namespace RemodelLoopMover { BooleanOptionKey const force_cutting_N( "remodel:RemodelLoopMover:force_cutting_N" );  } }
namespace remodel { namespace RemodelLoopMover { BooleanOptionKey const bypass_closure( "remodel:RemodelLoopMover:bypass_closure" );  } }
namespace remodel { namespace RemodelLoopMover { BooleanOptionKey const cyclic_peptide( "remodel:RemodelLoopMover:cyclic_peptide" );  } }
namespace remodel { namespace RemodelLoopMover { RealOptionKey const temperature( "remodel:RemodelLoopMover:temperature" );  } }
namespace remodel { namespace RemodelLoopMover { IntegerOptionKey const max_chews( "remodel:RemodelLoopMover:max_chews" );  } }
namespace fold_from_loops { BooleanOptionKey const fold_from_loops( "fold_from_loops" );  }
namespace fold_from_loops { BooleanOptionKey const native_ca_cst( "fold_from_loops:native_ca_cst" );  }
namespace fold_from_loops { FileOptionKey const swap_loops( "fold_from_loops:swap_loops" );  }
namespace fold_from_loops { StringOptionKey const checkpoint( "fold_from_loops:checkpoint" );  }
namespace fold_from_loops { RealOptionKey const ca_csts_dev( "fold_from_loops:ca_csts_dev" );  }
namespace fold_from_loops { IntegerOptionKey const add_relax_cycles( "fold_from_loops:add_relax_cycles" );  }
namespace fold_from_loops { IntegerOptionKey const loop_mov_nterm( "fold_from_loops:loop_mov_nterm" );  }
namespace fold_from_loops { IntegerOptionKey const loop_mov_cterm( "fold_from_loops:loop_mov_cterm" );  }
namespace fold_from_loops { RealOptionKey const ca_rmsd_cutoff( "fold_from_loops:ca_rmsd_cutoff" );  }
namespace fold_from_loops { IntegerVectorOptionKey const res_design_bs( "fold_from_loops:res_design_bs" );  }
namespace fold_from_loops { FileOptionKey const clear_csts( "fold_from_loops:clear_csts" );  }
namespace fold_from_loops { BooleanOptionKey const output_centroid( "fold_from_loops:output_centroid" );  }
namespace fold_from_loops { BooleanOptionKey const add_cst_loop( "fold_from_loops:add_cst_loop" );  }
namespace symmetry { BooleanOptionKey const symmetry( "symmetry" );  }
namespace symmetry { StringOptionKey const symmetry_definition( "symmetry:symmetry_definition" );  }
namespace symmetry { RealOptionKey const reweight_symm_interactions( "symmetry:reweight_symm_interactions" );  }
namespace symmetry { BooleanOptionKey const initialize_rigid_body_dofs( "symmetry:initialize_rigid_body_dofs" );  }
namespace symmetry { RealVectorOptionKey const perturb_rigid_body_dofs( "symmetry:perturb_rigid_body_dofs" );  }
namespace symmetry { BooleanOptionKey const symmetric_rmsd( "symmetry:symmetric_rmsd" );  }
namespace fold_and_dock { BooleanOptionKey const fold_and_dock( "fold_and_dock" );  }
namespace fold_and_dock { BooleanOptionKey const move_anchor_points( "fold_and_dock:move_anchor_points" );  }
namespace fold_and_dock { BooleanOptionKey const set_anchor_at_closest_point( "fold_and_dock:set_anchor_at_closest_point" );  }
namespace fold_and_dock { BooleanOptionKey const rotate_anchor_to_x( "fold_and_dock:rotate_anchor_to_x" );  }
namespace fold_and_dock { RealOptionKey const trans_mag_smooth( "fold_and_dock:trans_mag_smooth" );  }
namespace fold_and_dock { RealOptionKey const rot_mag_smooth( "fold_and_dock:rot_mag_smooth" );  }
namespace fold_and_dock { RealOptionKey const rb_rot_magnitude( "fold_and_dock:rb_rot_magnitude" );  }
namespace fold_and_dock { RealOptionKey const rb_trans_magnitude( "fold_and_dock:rb_trans_magnitude" );  }
namespace fold_and_dock { IntegerOptionKey const rigid_body_cycles( "fold_and_dock:rigid_body_cycles" );  }
namespace fold_and_dock { RealOptionKey const move_anchor_frequency( "fold_and_dock:move_anchor_frequency" );  }
namespace fold_and_dock { RealOptionKey const rigid_body_frequency( "fold_and_dock:rigid_body_frequency" );  }
namespace fold_and_dock { BooleanOptionKey const rigid_body_disable_mc( "fold_and_dock:rigid_body_disable_mc" );  }
namespace fold_and_dock { RealOptionKey const slide_contact_frequency( "fold_and_dock:slide_contact_frequency" );  }
namespace match { BooleanOptionKey const match( "match" );  }
namespace match { StringOptionKey const lig_name( "match:lig_name" );  }
namespace match { RealOptionKey const bump_tolerance( "match:bump_tolerance" );  }
namespace match { FileOptionKey const active_site_definition_by_residue( "match:active_site_definition_by_residue" );  }
namespace match { FileOptionKey const active_site_definition_by_gridlig( "match:active_site_definition_by_gridlig" );  }
namespace match { FileOptionKey const required_active_site_atom_names( "match:required_active_site_atom_names" );  }
namespace match { FileOptionKey const grid_boundary( "match:grid_boundary" );  }
namespace match { FileOptionKey const geometric_constraint_file( "match:geometric_constraint_file" );  }
namespace match { FileOptionKey const scaffold_active_site_residues( "match:scaffold_active_site_residues" );  }
namespace match { FileOptionKey const scaffold_active_site_residues_for_geomcsts( "match:scaffold_active_site_residues_for_geomcsts" );  }
namespace match { RealOptionKey const euclid_bin_size( "match:euclid_bin_size" );  }
namespace match { RealOptionKey const euler_bin_size( "match:euler_bin_size" );  }
namespace match { BooleanOptionKey const consolidate_matches( "match:consolidate_matches" );  }
namespace match { IntegerOptionKey const output_matches_per_group( "match:output_matches_per_group" );  }
namespace match { StringVectorOptionKey const orientation_atoms( "match:orientation_atoms" );  }
namespace match { StringOptionKey const output_format( "match:output_format" );  }
namespace match { StringOptionKey const match_grouper( "match:match_grouper" );  }
namespace match { RealOptionKey const grouper_downstream_rmsd( "match:grouper_downstream_rmsd" );  }
namespace match { BooleanOptionKey const output_matchres_only( "match:output_matchres_only" );  }
namespace match { IntegerVectorOptionKey const geom_csts_downstream_output( "match:geom_csts_downstream_output" );  }
namespace match { BooleanOptionKey const filter_colliding_upstream_residues( "match:filter_colliding_upstream_residues" );  }
namespace match { RealOptionKey const upstream_residue_collision_tolerance( "match:upstream_residue_collision_tolerance" );  }
namespace match { RealOptionKey const upstream_residue_collision_score_cutoff( "match:upstream_residue_collision_score_cutoff" );  }
namespace match { RealOptionKey const upstream_residue_collision_Wfa_atr( "match:upstream_residue_collision_Wfa_atr" );  }
namespace match { RealOptionKey const upstream_residue_collision_Wfa_rep( "match:upstream_residue_collision_Wfa_rep" );  }
namespace match { RealOptionKey const upstream_residue_collision_Wfa_sol( "match:upstream_residue_collision_Wfa_sol" );  }
namespace match { BooleanOptionKey const filter_upstream_downstream_collisions( "match:filter_upstream_downstream_collisions" );  }
namespace match { RealOptionKey const updown_collision_tolerance( "match:updown_collision_tolerance" );  }
namespace match { RealOptionKey const updown_residue_collision_score_cutoff( "match:updown_residue_collision_score_cutoff" );  }
namespace match { RealOptionKey const updown_residue_collision_Wfa_atr( "match:updown_residue_collision_Wfa_atr" );  }
namespace match { RealOptionKey const updown_residue_collision_Wfa_rep( "match:updown_residue_collision_Wfa_rep" );  }
namespace match { RealOptionKey const updown_residue_collision_Wfa_sol( "match:updown_residue_collision_Wfa_sol" );  }
namespace match { BooleanOptionKey const define_match_by_single_downstream_positioning( "match:define_match_by_single_downstream_positioning" );  }
namespace match { IntegerOptionKey const ligand_rotamer_index( "match:ligand_rotamer_index" );  }
namespace match { BooleanOptionKey const enumerate_ligand_rotamers( "match:enumerate_ligand_rotamers" );  }
namespace match { BooleanOptionKey const only_enumerate_non_match_redundant_ligand_rotamers( "match:only_enumerate_non_match_redundant_ligand_rotamers" );  }
namespace match { BooleanOptionKey const dynamic_grid_refinement( "match:dynamic_grid_refinement" );  }
namespace canonical_sampling { BooleanOptionKey const canonical_sampling( "canonical_sampling" );  }
namespace canonical_sampling { namespace probabilities { BooleanOptionKey const probabilities( "canonical_sampling:probabilities" );  } }
namespace canonical_sampling { namespace probabilities { RealOptionKey const sc( "canonical_sampling:probabilities:sc" );  } }
namespace canonical_sampling { namespace probabilities { RealOptionKey const localbb( "canonical_sampling:probabilities:localbb" );  } }
namespace canonical_sampling { namespace probabilities { RealOptionKey const sc_prob_uniform( "canonical_sampling:probabilities:sc_prob_uniform" );  } }
namespace canonical_sampling { namespace probabilities { RealOptionKey const sc_prob_withinrot( "canonical_sampling:probabilities:sc_prob_withinrot" );  } }
namespace canonical_sampling { namespace probabilities { RealOptionKey const sc_prob_perturbcurrent( "canonical_sampling:probabilities:sc_prob_perturbcurrent" );  } }
namespace canonical_sampling { namespace probabilities { BooleanOptionKey const MPI_sync_pools( "canonical_sampling:probabilities:MPI_sync_pools" );  } }
namespace canonical_sampling { namespace probabilities { BooleanOptionKey const MPI_bcast( "canonical_sampling:probabilities:MPI_bcast" );  } }
namespace canonical_sampling { namespace probabilities { BooleanOptionKey const fast_sc_moves( "canonical_sampling:probabilities:fast_sc_moves" );  } }
namespace canonical_sampling { namespace probabilities { RealOptionKey const fast_sc_moves_ntrials( "canonical_sampling:probabilities:fast_sc_moves_ntrials" );  } }
namespace canonical_sampling { namespace probabilities { BooleanOptionKey const no_jd2_output( "canonical_sampling:probabilities:no_jd2_output" );  } }
namespace canonical_sampling { namespace probabilities { BooleanOptionKey const use_hierarchical_clustering( "canonical_sampling:probabilities:use_hierarchical_clustering" );  } }
namespace canonical_sampling { namespace probabilities { IntegerOptionKey const hierarchical_max_cache_size( "canonical_sampling:probabilities:hierarchical_max_cache_size" );  } }
namespace canonical_sampling { namespace probabilities { RealOptionKey const backrub( "canonical_sampling:probabilities:backrub" );  } }
namespace canonical_sampling { namespace probabilities { RealOptionKey const conrot( "canonical_sampling:probabilities:conrot" );  } }
namespace canonical_sampling { namespace sampling { BooleanOptionKey const sampling( "canonical_sampling:sampling" );  } }
namespace canonical_sampling { namespace sampling { BooleanOptionKey const no_detailed_balance( "canonical_sampling:sampling:no_detailed_balance" );  } }
namespace canonical_sampling { namespace sampling { IntegerOptionKey const ntrials( "canonical_sampling:sampling:ntrials" );  } }
namespace canonical_sampling { namespace sampling { RealOptionKey const mc_kt( "canonical_sampling:sampling:mc_kt" );  } }
namespace canonical_sampling { namespace sampling { IntegerOptionKey const interval_pose_dump( "canonical_sampling:sampling:interval_pose_dump" );  } }
namespace canonical_sampling { namespace sampling { IntegerOptionKey const interval_data_dump( "canonical_sampling:sampling:interval_data_dump" );  } }
namespace canonical_sampling { namespace sampling { BooleanOptionKey const output_only_cluster_transitions( "canonical_sampling:sampling:output_only_cluster_transitions" );  } }
namespace canonical_sampling { namespace sampling { RealOptionKey const transition_threshold( "canonical_sampling:sampling:transition_threshold" );  } }
namespace canonical_sampling { namespace sampling { IntegerOptionKey const max_files_per_dir( "canonical_sampling:sampling:max_files_per_dir" );  } }
namespace canonical_sampling { namespace sampling { BooleanOptionKey const save_loops_only( "canonical_sampling:sampling:save_loops_only" );  } }
namespace canonical_sampling { namespace sampling { BooleanOptionKey const dump_loops_only( "canonical_sampling:sampling:dump_loops_only" );  } }
namespace canonical_sampling { namespace out { BooleanOptionKey const out( "canonical_sampling:out" );  } }
namespace canonical_sampling { namespace out { FileOptionKey const new_structures( "canonical_sampling:out:new_structures" );  } }
namespace rdc { BooleanOptionKey const rdc( "rdc" );  }
namespace rdc { BooleanOptionKey const correct_NH_length( "rdc:correct_NH_length" );  }
namespace rdc { BooleanOptionKey const reduced_couplings( "rdc:reduced_couplings" );  }
namespace rdc { FileOptionKey const weights( "rdc:weights" );  }
namespace rdc { RealOptionKey const iterate_weights( "rdc:iterate_weights" );  }
namespace rdc { FileOptionKey const segment_file( "rdc:segment_file" );  }
namespace rdc { StringOptionKey const segment_scoring_mode( "rdc:segment_scoring_mode" );  }
namespace rdc { RealOptionKey const total_weight( "rdc:total_weight" );  }
namespace rdc { RealOptionKey const tensor_weight( "rdc:tensor_weight" );  }
namespace rdc { FileOptionKey const print_rdc_values( "rdc:print_rdc_values" );  }
namespace rdc { RealOptionKey const iterate_tol( "rdc:iterate_tol" );  }
namespace rdc { BooleanOptionKey const iterate_reset( "rdc:iterate_reset" );  }
namespace rdc { FileOptionKey const dump_weight_trajectory( "rdc:dump_weight_trajectory" );  }
namespace rdc { RealVectorOptionKey const fix_normAzz( "rdc:fix_normAzz" );  }
namespace rdc { FileOptionKey const select_residues_file( "rdc:select_residues_file" );  }
namespace rdc { StringOptionKey const fit_method( "rdc:fit_method" );  }
namespace rdc { RealVectorOptionKey const fixDa( "rdc:fixDa" );  }
namespace rdc { RealVectorOptionKey const fixR( "rdc:fixR" );  }
namespace csa { BooleanOptionKey const csa( "csa" );  }
namespace csa { BooleanOptionKey const useZ( "csa:useZ" );  }
namespace dc { BooleanOptionKey const dc( "dc" );  }
namespace dc { BooleanOptionKey const useZ( "dc:useZ" );  }
namespace antibody { BooleanOptionKey const antibody( "antibody" );  }
namespace antibody { BooleanOptionKey const graft_l1( "antibody:graft_l1" );  }
namespace antibody { StringOptionKey const l1_template( "antibody:l1_template" );  }
namespace antibody { BooleanOptionKey const graft_l2( "antibody:graft_l2" );  }
namespace antibody { StringOptionKey const l2_template( "antibody:l2_template" );  }
namespace antibody { BooleanOptionKey const graft_l3( "antibody:graft_l3" );  }
namespace antibody { StringOptionKey const l3_template( "antibody:l3_template" );  }
namespace antibody { BooleanOptionKey const graft_h1( "antibody:graft_h1" );  }
namespace antibody { StringOptionKey const h1_template( "antibody:h1_template" );  }
namespace antibody { BooleanOptionKey const graft_h2( "antibody:graft_h2" );  }
namespace antibody { StringOptionKey const h2_template( "antibody:h2_template" );  }
namespace antibody { BooleanOptionKey const graft_h3( "antibody:graft_h3" );  }
namespace antibody { StringOptionKey const h3_template( "antibody:h3_template" );  }
namespace antibody { BooleanOptionKey const h3_no_stem_graft( "antibody:h3_no_stem_graft" );  }
namespace antibody { BooleanOptionKey const packonly_after_graft( "antibody:packonly_after_graft" );  }
namespace antibody { BooleanOptionKey const model_h3( "antibody:model_h3" );  }
namespace antibody { BooleanOptionKey const snugfit( "antibody:snugfit" );  }
namespace antibody { BooleanOptionKey const refine_h3( "antibody:refine_h3" );  }
namespace antibody { BooleanOptionKey const h3_filter( "antibody:h3_filter" );  }
namespace antibody { RealOptionKey const h3_filter_tolerance( "antibody:h3_filter_tolerance" );  }
namespace antibody { BooleanOptionKey const cter_insert( "antibody:cter_insert" );  }
namespace antibody { BooleanOptionKey const flank_residue_min( "antibody:flank_residue_min" );  }
namespace antibody { BooleanOptionKey const sc_min( "antibody:sc_min" );  }
namespace antibody { BooleanOptionKey const rt_min( "antibody:rt_min" );  }
namespace antibody { BooleanOptionKey const bad_nter( "antibody:bad_nter" );  }
namespace antibody { StringOptionKey const remodel( "antibody:remodel" );  }
namespace antibody { StringOptionKey const refine( "antibody:refine" );  }
namespace antibody { StringOptionKey const centroid_refine( "antibody:centroid_refine" );  }
namespace antibody { BooleanOptionKey const constrain_cter( "antibody:constrain_cter" );  }
namespace antibody { BooleanOptionKey const constrain_vlvh_qq( "antibody:constrain_vlvh_qq" );  }
namespace antibody { BooleanOptionKey const snug_loops( "antibody:snug_loops" );  }
namespace antibody { FileOptionKey const input_fv( "antibody:input_fv" );  }
namespace antibody { BooleanOptionKey const camelid( "antibody:camelid" );  }
namespace antibody { BooleanOptionKey const camelid_constraints( "antibody:camelid_constraints" );  }
namespace antibody { StringOptionKey const numbering_scheme( "antibody:numbering_scheme" );  }
namespace flexPepDocking { BooleanOptionKey const flexPepDocking( "flexPepDocking" );  }
namespace flexPepDocking { StringOptionKey const params_file( "flexPepDocking:params_file" );  }
namespace flexPepDocking { IntegerOptionKey const peptide_anchor( "flexPepDocking:peptide_anchor" );  }
namespace flexPepDocking { StringOptionKey const receptor_chain( "flexPepDocking:receptor_chain" );  }
namespace flexPepDocking { StringOptionKey const peptide_chain( "flexPepDocking:peptide_chain" );  }
namespace flexPepDocking { BooleanOptionKey const pep_fold_only( "flexPepDocking:pep_fold_only" );  }
namespace flexPepDocking { BooleanOptionKey const lowres_abinitio( "flexPepDocking:lowres_abinitio" );  }
namespace flexPepDocking { BooleanOptionKey const lowres_preoptimize( "flexPepDocking:lowres_preoptimize" );  }
namespace flexPepDocking { BooleanOptionKey const flexPepDockingMinimizeOnly( "flexPepDocking:flexPepDockingMinimizeOnly" );  }
namespace flexPepDocking { BooleanOptionKey const extend_peptide( "flexPepDocking:extend_peptide" );  }
namespace flexPepDocking { BooleanOptionKey const pep_refine( "flexPepDocking:pep_refine" );  }
namespace flexPepDocking { BooleanOptionKey const rbMCM( "flexPepDocking:rbMCM" );  }
namespace flexPepDocking { BooleanOptionKey const torsionsMCM( "flexPepDocking:torsionsMCM" );  }
namespace flexPepDocking { BooleanOptionKey const peptide_loop_model( "flexPepDocking:peptide_loop_model" );  }
namespace flexPepDocking { BooleanOptionKey const backrub_peptide( "flexPepDocking:backrub_peptide" );  }
namespace flexPepDocking { BooleanOptionKey const boost_fa_atr( "flexPepDocking:boost_fa_atr" );  }
namespace flexPepDocking { BooleanOptionKey const ramp_fa_rep( "flexPepDocking:ramp_fa_rep" );  }
namespace flexPepDocking { BooleanOptionKey const ramp_rama( "flexPepDocking:ramp_rama" );  }
namespace flexPepDocking { BooleanOptionKey const flexpep_score_only( "flexPepDocking:flexpep_score_only" );  }
namespace flexPepDocking { FileOptionKey const ref_startstruct( "flexPepDocking:ref_startstruct" );  }
namespace flexPepDocking { BooleanOptionKey const use_cen_score( "flexPepDocking:use_cen_score" );  }
namespace flexPepDocking { BooleanOptionKey const design_peptide( "flexPepDocking:design_peptide" );  }
namespace flexPepDocking { IntegerOptionKey const rep_ramp_cycles( "flexPepDocking:rep_ramp_cycles" );  }
namespace flexPepDocking { IntegerOptionKey const mcm_cycles( "flexPepDocking:mcm_cycles" );  }
namespace flexPepDocking { RealOptionKey const random_phi_psi_preturbation( "flexPepDocking:random_phi_psi_preturbation" );  }
namespace flexPepDocking { RealOptionKey const smove_angle_range( "flexPepDocking:smove_angle_range" );  }
namespace flexPepDocking { BooleanOptionKey const min_receptor_bb( "flexPepDocking:min_receptor_bb" );  }
namespace flexPepDocking { RealOptionKey const random_trans_start( "flexPepDocking:random_trans_start" );  }
namespace flexPepDocking { RealOptionKey const random_rot_start( "flexPepDocking:random_rot_start" );  }
namespace flexPepDocking { BooleanOptionKey const flexpep_prepack( "flexPepDocking:flexpep_prepack" );  }
namespace flexPepDocking { BooleanOptionKey const flexpep_noprepack1( "flexPepDocking:flexpep_noprepack1" );  }
namespace flexPepDocking { BooleanOptionKey const flexpep_noprepack2( "flexPepDocking:flexpep_noprepack2" );  }
namespace flexPepDocking { RealOptionKey const score_filter( "flexPepDocking:score_filter" );  }
namespace flexPepDocking { IntegerOptionKey const hb_filter( "flexPepDocking:hb_filter" );  }
namespace flexPepDocking { IntegerOptionKey const hotspot_filter( "flexPepDocking:hotspot_filter" );  }
namespace flexPepDocking { StringOptionKey const frag5( "flexPepDocking:frag5" );  }
namespace flexPepDocking { RealOptionKey const frag9_weight( "flexPepDocking:frag9_weight" );  }
namespace flexPepDocking { RealOptionKey const frag5_weight( "flexPepDocking:frag5_weight" );  }
namespace flexPepDocking { RealOptionKey const frag3_weight( "flexPepDocking:frag3_weight" );  }
namespace flexPepDocking { BooleanOptionKey const pSer2Asp_centroid( "flexPepDocking:pSer2Asp_centroid" );  }
namespace flexPepDocking { BooleanOptionKey const pSer2Glu_centroid( "flexPepDocking:pSer2Glu_centroid" );  }
namespace flexPepDocking { BooleanOptionKey const dumpPDB_abinitio( "flexPepDocking:dumpPDB_abinitio" );  }
namespace flexPepDocking { BooleanOptionKey const dumpPDB_lowres( "flexPepDocking:dumpPDB_lowres" );  }
namespace flexPepDocking { BooleanOptionKey const dumpPDB_hires( "flexPepDocking:dumpPDB_hires" );  }
namespace threadsc { BooleanOptionKey const threadsc( "threadsc" );  }
namespace threadsc { StringOptionKey const src_chain( "threadsc:src_chain" );  }
namespace threadsc { StringOptionKey const trg_chain( "threadsc:trg_chain" );  }
namespace threadsc { IntegerOptionKey const src_first_resid( "threadsc:src_first_resid" );  }
namespace threadsc { IntegerOptionKey const trg_first_resid( "threadsc:trg_first_resid" );  }
namespace threadsc { IntegerOptionKey const nres( "threadsc:nres" );  }
namespace threadsc { IntegerOptionKey const trg_anchor( "threadsc:trg_anchor" );  }
namespace cp { BooleanOptionKey const cp( "cp" );  }
namespace cp { RealOptionKey const cutoff( "cp:cutoff" );  }
namespace cp { StringOptionKey const minimizer( "cp:minimizer" );  }
namespace cp { StringOptionKey const relax_sfxn( "cp:relax_sfxn" );  }
namespace cp { StringOptionKey const pack_sfxn( "cp:pack_sfxn" );  }
namespace cp { RealOptionKey const minimizer_tol( "cp:minimizer_tol" );  }
namespace cp { StringOptionKey const minimizer_score_fxn( "cp:minimizer_score_fxn" );  }
namespace cp { StringOptionKey const output( "cp:output" );  }
namespace cp { IntegerOptionKey const ncycles( "cp:ncycles" );  }
namespace cp { IntegerOptionKey const max_failures( "cp:max_failures" );  }
namespace cp { BooleanOptionKey const print_reports( "cp:print_reports" );  }
namespace cp { StringOptionKey const vipReportFile( "cp:vipReportFile" );  }
namespace cp { StringOptionKey const exclude_file( "cp:exclude_file" );  }
namespace cp { StringOptionKey const relax_mover( "cp:relax_mover" );  }
namespace cp { BooleanOptionKey const skip_relax( "cp:skip_relax" );  }
namespace cp { BooleanOptionKey const local_relax( "cp:local_relax" );  }
namespace cp { BooleanOptionKey const print_intermediate_pdbs( "cp:print_intermediate_pdbs" );  }
namespace archive { BooleanOptionKey const archive( "archive" );  }
namespace archive { BooleanOptionKey const reread_all_structures( "archive:reread_all_structures" );  }
namespace archive { IntegerOptionKey const completion_notify_frequency( "archive:completion_notify_frequency" );  }
namespace optimization { BooleanOptionKey const optimization( "optimization" );  }
namespace optimization { IntegerOptionKey const default_max_cycles( "optimization:default_max_cycles" );  }
namespace optimization { RealOptionKey const armijo_min_stepsize( "optimization:armijo_min_stepsize" );  }
namespace optimization { IntegerOptionKey const lbfgs_M( "optimization:lbfgs_M" );  }
namespace optimization { RealOptionKey const scale_d( "optimization:scale_d" );  }
namespace optimization { RealOptionKey const scale_theta( "optimization:scale_theta" );  }
namespace optimization { RealOptionKey const scale_rb( "optimization:scale_rb" );  }
namespace optimization { RealOptionKey const scale_rbangle( "optimization:scale_rbangle" );  }
namespace optimization { BooleanOptionKey const scmin_nonideal( "optimization:scmin_nonideal" );  }
namespace optimization { BooleanOptionKey const scmin_cartesian( "optimization:scmin_cartesian" );  }
namespace swa { BooleanOptionKey const swa( "swa" );  }
namespace swa { StringVectorOptionKey const s1( "swa:s1" );  }
namespace swa { StringVectorOptionKey const s2( "swa:s2" );  }
namespace swa { StringVectorOptionKey const silent1( "swa:silent1" );  }
namespace swa { StringVectorOptionKey const silent2( "swa:silent2" );  }
namespace swa { StringVectorOptionKey const tags1( "swa:tags1" );  }
namespace swa { StringVectorOptionKey const tags2( "swa:tags2" );  }
namespace swa { IntegerVectorOptionKey const slice_res1( "swa:slice_res1" );  }
namespace swa { IntegerVectorOptionKey const slice_res2( "swa:slice_res2" );  }
namespace swa { IntegerVectorOptionKey const input_res1( "swa:input_res1" );  }
namespace swa { IntegerVectorOptionKey const input_res2( "swa:input_res2" );  }
namespace swa { BooleanOptionKey const backbone_only1( "swa:backbone_only1" );  }
namespace swa { BooleanOptionKey const backbone_only2( "swa:backbone_only2" );  }
namespace ufv { BooleanOptionKey const ufv( "ufv" );  }
namespace ufv { IntegerOptionKey const left( "ufv:left" );  }
namespace ufv { IntegerOptionKey const right( "ufv:right" );  }
namespace ufv { StringOptionKey const ss( "ufv:ss" );  }
namespace ufv { StringOptionKey const aa_during_build( "ufv:aa_during_build" );  }
namespace ufv { StringOptionKey const aa_during_design_refine( "ufv:aa_during_design_refine" );  }
namespace ufv { BooleanOptionKey const keep_junction_torsions( "ufv:keep_junction_torsions" );  }
namespace ufv { FileOptionKey const ufv_loops( "ufv:ufv_loops" );  }
namespace ufv { BooleanOptionKey const use_fullmer( "ufv:use_fullmer" );  }
namespace ufv { StringOptionKey const centroid_loop_mover( "ufv:centroid_loop_mover" );  }
namespace ufv { BooleanOptionKey const no_neighborhood_design( "ufv:no_neighborhood_design" );  }
namespace ufv { IntegerOptionKey const dr_cycles( "ufv:dr_cycles" );  }
namespace ufv { StringOptionKey const centroid_sfx( "ufv:centroid_sfx" );  }
namespace ufv { StringOptionKey const centroid_sfx_patch( "ufv:centroid_sfx_patch" );  }
namespace ufv { StringOptionKey const fullatom_sfx( "ufv:fullatom_sfx" );  }
namespace ufv { StringOptionKey const fullatom_sfx_patch( "ufv:fullatom_sfx_patch" );  }
namespace ufv { namespace insert { BooleanOptionKey const insert( "ufv:insert" );  } }
namespace ufv { namespace insert { FileOptionKey const insert_pdb( "ufv:insert:insert_pdb" );  } }
namespace ufv { namespace insert { FileOptionKey const attached_pdb( "ufv:insert:attached_pdb" );  } }
namespace ufv { namespace insert { StringOptionKey const connection_scheme( "ufv:insert:connection_scheme" );  } }
namespace chrisk { BooleanOptionKey const chrisk( "chrisk" );  }
namespace chrisk { BooleanOptionKey const hb_elec( "chrisk:hb_elec" );  }
namespace rot_anl { BooleanOptionKey const rot_anl( "rot_anl" );  }
namespace rot_anl { StringOptionKey const tag( "rot_anl:tag" );  }
namespace rot_anl { BooleanOptionKey const premin( "rot_anl:premin" );  }
namespace rot_anl { BooleanOptionKey const min( "rot_anl:min" );  }
namespace rot_anl { BooleanOptionKey const diff_to_min( "rot_anl:diff_to_min" );  }
namespace rot_anl { BooleanOptionKey const repack( "rot_anl:repack" );  }
namespace rot_anl { BooleanOptionKey const rtmin( "rot_anl:rtmin" );  }
namespace rot_anl { BooleanOptionKey const scmove( "rot_anl:scmove" );  }
namespace rot_anl { BooleanOptionKey const design( "rot_anl:design" );  }
namespace rot_anl { RealOptionKey const score_tol( "rot_anl:score_tol" );  }
namespace rot_anl { RealOptionKey const rmsd_tol( "rot_anl:rmsd_tol" );  }
namespace rot_anl { BooleanOptionKey const dump_pdb( "rot_anl:dump_pdb" );  }
namespace rot_anl { IntegerOptionKey const nloop_scmove( "rot_anl:nloop_scmove" );  }
namespace helixAssembly { BooleanOptionKey const helixAssembly( "helixAssembly" );  }
namespace helixAssembly { FileOptionKey const query_structure_path( "helixAssembly:query_structure_path" );  }
namespace helixAssembly { IntegerOptionKey const frag1_start( "helixAssembly:frag1_start" );  }
namespace helixAssembly { IntegerOptionKey const frag1_end( "helixAssembly:frag1_end" );  }
namespace helixAssembly { IntegerOptionKey const frag2_start( "helixAssembly:frag2_start" );  }
namespace helixAssembly { IntegerOptionKey const frag2_end( "helixAssembly:frag2_end" );  }
namespace helixAssembly { IntegerOptionKey const minimum_helix_contacts( "helixAssembly:minimum_helix_contacts" );  }
namespace helixAssembly { IntegerOptionKey const helices_to_add( "helixAssembly:helices_to_add" );  }
namespace helixAssembly { RealOptionKey const single_helix_rmsd_cutoff( "helixAssembly:single_helix_rmsd_cutoff" );  }
namespace helixAssembly { RealOptionKey const helix_pair_rmsd_cutoff( "helixAssembly:helix_pair_rmsd_cutoff" );  }
namespace helixAssembly { RealOptionKey const helix_cap_dist_cutoff( "helixAssembly:helix_cap_dist_cutoff" );  }
namespace helixAssembly { RealOptionKey const helix_contact_dist_cutoff( "helixAssembly:helix_contact_dist_cutoff" );  }
namespace helixAssembly { IntegerOptionKey const min_helix_size( "helixAssembly:min_helix_size" );  }
namespace strand_assembly { BooleanOptionKey const strand_assembly( "strand_assembly" );  }
namespace strand_assembly { IntegerOptionKey const min_num_strands_to_deal( "strand_assembly:min_num_strands_to_deal" );  }
namespace strand_assembly { IntegerOptionKey const max_num_strands_to_deal( "strand_assembly:max_num_strands_to_deal" );  }
namespace strand_assembly { BooleanOptionKey const extract_native_only( "strand_assembly:extract_native_only" );  }
namespace strand_assembly { BooleanOptionKey const extract_sandwich( "strand_assembly:extract_sandwich" );  }
namespace strand_assembly { IntegerOptionKey const min_strand_size( "strand_assembly:min_strand_size" );  }
namespace strand_assembly { IntegerOptionKey const max_strand_size( "strand_assembly:max_strand_size" );  }
namespace strand_assembly { RealOptionKey const min_CA_CA_dis( "strand_assembly:min_CA_CA_dis" );  }
namespace strand_assembly { RealOptionKey const max_CA_CA_dis( "strand_assembly:max_CA_CA_dis" );  }
namespace strand_assembly { RealOptionKey const min_O_N_dis( "strand_assembly:min_O_N_dis" );  }
namespace strand_assembly { RealOptionKey const max_O_N_dis( "strand_assembly:max_O_N_dis" );  }
namespace strand_assembly { RealOptionKey const min_C_O_N_angle( "strand_assembly:min_C_O_N_angle" );  }
namespace strand_assembly { RealOptionKey const min_sheet_dis( "strand_assembly:min_sheet_dis" );  }
namespace strand_assembly { RealOptionKey const max_sheet_dis( "strand_assembly:max_sheet_dis" );  }
namespace strand_assembly { RealOptionKey const min_sheet_torsion( "strand_assembly:min_sheet_torsion" );  }
namespace strand_assembly { RealOptionKey const max_sheet_torsion( "strand_assembly:max_sheet_torsion" );  }
namespace strand_assembly { RealOptionKey const min_sheet_torsion_cen_res( "strand_assembly:min_sheet_torsion_cen_res" );  }
namespace strand_assembly { RealOptionKey const max_sheet_torsion_cen_res( "strand_assembly:max_sheet_torsion_cen_res" );  }
namespace strand_assembly { RealOptionKey const min_sheet_angle( "strand_assembly:min_sheet_angle" );  }
namespace strand_assembly { RealOptionKey const max_sheet_angle( "strand_assembly:max_sheet_angle" );  }
namespace strand_assembly { IntegerOptionKey const min_num_strands_in_sheet( "strand_assembly:min_num_strands_in_sheet" );  }
namespace strand_assembly { RealOptionKey const min_inter_sheet_dis_CA_CA( "strand_assembly:min_inter_sheet_dis_CA_CA" );  }
namespace strand_assembly { RealOptionKey const max_inter_sheet_dis_CA_CA( "strand_assembly:max_inter_sheet_dis_CA_CA" );  }
namespace strand_assembly { RealOptionKey const min_shortest_dis_sidechain_inter_sheet( "strand_assembly:min_shortest_dis_sidechain_inter_sheet" );  }
namespace strand_assembly { BooleanOptionKey const write_chain_B_resnum( "strand_assembly:write_chain_B_resnum" );  }
namespace strand_assembly { BooleanOptionKey const write_phi_psi( "strand_assembly:write_phi_psi" );  }
namespace strand_assembly { BooleanOptionKey const no_helix_in_pdb( "strand_assembly:no_helix_in_pdb" );  }
namespace strand_assembly { IntegerOptionKey const max_helix_in_extracted_sw_loop( "strand_assembly:max_helix_in_extracted_sw_loop" );  }
namespace strand_assembly { BooleanOptionKey const no_strand_in_loop_in_extracted_sw( "strand_assembly:no_strand_in_loop_in_extracted_sw" );  }
namespace strand_assembly { BooleanOptionKey const exclude_sandwich_that_is_linked_w_same_direction_strand( "strand_assembly:exclude_sandwich_that_is_linked_w_same_direction_strand" );  }
namespace strand_assembly { RealOptionKey const max_inter_strand_angle_to_not_be_same_direction_strands( "strand_assembly:max_inter_strand_angle_to_not_be_same_direction_strands" );  }
namespace strand_assembly { RealOptionKey const max_abs_inter_strand_dihedral_to_not_be_same_direction_strands( "strand_assembly:max_abs_inter_strand_dihedral_to_not_be_same_direction_strands" );  }
namespace strand_assembly { IntegerOptionKey const max_num_sw_per_pdb( "strand_assembly:max_num_sw_per_pdb" );  }
namespace strand_assembly { StringOptionKey const check_N_to_C_direction_by( "strand_assembly:check_N_to_C_direction_by" );  }
namespace strand_assembly { RealOptionKey const check_canonicalness_cutoff( "strand_assembly:check_canonicalness_cutoff" );  }
namespace pepspec { BooleanOptionKey const pepspec( "pepspec" );  }
namespace pepspec { StringOptionKey const soft_wts( "pepspec:soft_wts" );  }
namespace pepspec { StringOptionKey const cen_wts( "pepspec:cen_wts" );  }
namespace pepspec { BooleanOptionKey const binding_score( "pepspec:binding_score" );  }
namespace pepspec { BooleanOptionKey const no_cen( "pepspec:no_cen" );  }
namespace pepspec { BooleanOptionKey const no_cen_rottrials( "pepspec:no_cen_rottrials" );  }
namespace pepspec { BooleanOptionKey const run_sequential( "pepspec:run_sequential" );  }
namespace pepspec { IntegerOptionKey const pep_anchor( "pepspec:pep_anchor" );  }
namespace pepspec { StringOptionKey const pep_chain( "pepspec:pep_chain" );  }
namespace pepspec { IntegerOptionKey const n_peptides( "pepspec:n_peptides" );  }
namespace pepspec { IntegerOptionKey const n_build_loop( "pepspec:n_build_loop" );  }
namespace pepspec { IntegerOptionKey const n_cgrelax_loop( "pepspec:n_cgrelax_loop" );  }
namespace pepspec { IntegerOptionKey const n_dock_loop( "pepspec:n_dock_loop" );  }
namespace pepspec { RealOptionKey const interface_cutoff( "pepspec:interface_cutoff" );  }
namespace pepspec { BooleanOptionKey const use_input_bb( "pepspec:use_input_bb" );  }
namespace pepspec { BooleanOptionKey const remove_input_bb( "pepspec:remove_input_bb" );  }
namespace pepspec { StringOptionKey const homol_csts( "pepspec:homol_csts" );  }
namespace pepspec { RealOptionKey const p_homol_csts( "pepspec:p_homol_csts" );  }
namespace pepspec { StringOptionKey const frag_file( "pepspec:frag_file" );  }
namespace pepspec { BooleanOptionKey const gen_pep_bb_sequential( "pepspec:gen_pep_bb_sequential" );  }
namespace pepspec { StringOptionKey const input_seq( "pepspec:input_seq" );  }
namespace pepspec { StringOptionKey const ss_type( "pepspec:ss_type" );  }
namespace pepspec { BooleanOptionKey const upweight_interface( "pepspec:upweight_interface" );  }
namespace pepspec { BooleanOptionKey const calc_sasa( "pepspec:calc_sasa" );  }
namespace pepspec { BooleanOptionKey const diversify_pep_seqs( "pepspec:diversify_pep_seqs" );  }
namespace pepspec { IntegerOptionKey const diversify_lvl( "pepspec:diversify_lvl" );  }
namespace pepspec { BooleanOptionKey const dump_cg_bb( "pepspec:dump_cg_bb" );  }
namespace pepspec { BooleanOptionKey const save_low_pdbs( "pepspec:save_low_pdbs" );  }
namespace pepspec { BooleanOptionKey const save_all_pdbs( "pepspec:save_all_pdbs" );  }
namespace pepspec { BooleanOptionKey const no_design( "pepspec:no_design" );  }
namespace pepspec { StringOptionKey const pdb_list( "pepspec:pdb_list" );  }
namespace pepspec { StringOptionKey const ref_pdb_list( "pepspec:ref_pdb_list" );  }
namespace pepspec { BooleanOptionKey const add_buffer_res( "pepspec:add_buffer_res" );  }
namespace pepspec { StringOptionKey const cg_res_type( "pepspec:cg_res_type" );  }
namespace pepspec { IntegerOptionKey const native_pep_anchor( "pepspec:native_pep_anchor" );  }
namespace pepspec { StringOptionKey const native_pep_chain( "pepspec:native_pep_chain" );  }
namespace pepspec { BooleanOptionKey const native_align( "pepspec:native_align" );  }
namespace pepspec { BooleanOptionKey const rmsd_analysis( "pepspec:rmsd_analysis" );  }
namespace pepspec { BooleanOptionKey const phipsi_analysis( "pepspec:phipsi_analysis" );  }
namespace pepspec { StringOptionKey const anchor_type( "pepspec:anchor_type" );  }
namespace pepspec { BooleanOptionKey const no_prepack_prot( "pepspec:no_prepack_prot" );  }
namespace pepspec { BooleanOptionKey const prep_use_ref_rotamers( "pepspec:prep_use_ref_rotamers" );  }
namespace pepspec { IntegerOptionKey const n_prepend( "pepspec:n_prepend" );  }
namespace pepspec { IntegerOptionKey const n_append( "pepspec:n_append" );  }
namespace pepspec { RealOptionKey const clash_cutoff( "pepspec:clash_cutoff" );  }
namespace pepspec { RealOptionKey const n_anchor_dock_std_devs( "pepspec:n_anchor_dock_std_devs" );  }
namespace pepspec { RealOptionKey const prep_trans_std_dev( "pepspec:prep_trans_std_dev" );  }
namespace pepspec { RealOptionKey const prep_rot_std_dev( "pepspec:prep_rot_std_dev" );  }
namespace pepspec { BooleanOptionKey const seq_align( "pepspec:seq_align" );  }
namespace pepspec { StringOptionKey const prep_align_prot_to( "pepspec:prep_align_prot_to" );  }
namespace sicdock { BooleanOptionKey const sicdock( "sicdock" );  }
namespace sicdock { RealOptionKey const clash_dis( "sicdock:clash_dis" );  }
namespace sicdock { RealOptionKey const contact_dis( "sicdock:contact_dis" );  }
namespace sicdock { RealOptionKey const hash_2D_vs_3D( "sicdock:hash_2D_vs_3D" );  }
namespace sicdock { RealOptionKey const term_min_expose( "sicdock:term_min_expose" );  }
namespace sicdock { RealOptionKey const term_max_angle( "sicdock:term_max_angle" );  }
namespace orbitals { BooleanOptionKey const orbitals( "orbitals" );  }
namespace orbitals { BooleanOptionKey const Hpol( "orbitals:Hpol" );  }
namespace orbitals { BooleanOptionKey const Haro( "orbitals:Haro" );  }
namespace orbitals { BooleanOptionKey const bb_stats( "orbitals:bb_stats" );  }
namespace orbitals { BooleanOptionKey const sc_stats( "orbitals:sc_stats" );  }
namespace orbitals { BooleanOptionKey const orb_orb_stats( "orbitals:orb_orb_stats" );  }
namespace orbitals { BooleanOptionKey const sc_bb( "orbitals:sc_bb" );  }
namespace carbohydrates { BooleanOptionKey const carbohydrates( "carbohydrates" );  }
namespace carbohydrates { BooleanOptionKey const lock_rings( "carbohydrates:lock_rings" );  }
namespace dwkulp { BooleanOptionKey const dwkulp( "dwkulp" );  }
namespace dwkulp { StringOptionKey const forcePolyAAfragments( "dwkulp:forcePolyAAfragments" );  }
namespace matdes { BooleanOptionKey const matdes( "matdes" );  }
namespace matdes { IntegerOptionKey const num_subs_building_block( "matdes:num_subs_building_block" );  }
namespace matdes { IntegerOptionKey const num_subs_total( "matdes:num_subs_total" );  }
namespace matdes { StringOptionKey const pdbID( "matdes:pdbID" );  }
namespace matdes { StringOptionKey const prefix( "matdes:prefix" );  }
namespace matdes { RealVectorOptionKey const radial_disp( "matdes:radial_disp" );  }
namespace matdes { RealVectorOptionKey const angle( "matdes:angle" );  }
namespace matdes { StringOptionKey const tag( "matdes:tag" );  }
namespace matdes { namespace dock { BooleanOptionKey const dock( "matdes:dock" );  } }
namespace matdes { namespace dock { RealOptionKey const neg_r( "matdes:dock:neg_r" );  } }
namespace matdes { namespace dock { BooleanOptionKey const dump_pdb( "matdes:dock:dump_pdb" );  } }
namespace matdes { namespace dock { BooleanOptionKey const dump_chainA_only( "matdes:dock:dump_chainA_only" );  } }
namespace matdes { namespace design { BooleanOptionKey const design( "matdes:design" );  } }
namespace matdes { namespace design { RealOptionKey const contact_dist( "matdes:design:contact_dist" );  } }
namespace matdes { namespace design { RealOptionKey const grid_size_angle( "matdes:design:grid_size_angle" );  } }
namespace matdes { namespace design { RealOptionKey const grid_size_radius( "matdes:design:grid_size_radius" );  } }
namespace matdes { namespace design { IntegerOptionKey const grid_nsamp_angle( "matdes:design:grid_nsamp_angle" );  } }
namespace matdes { namespace design { IntegerOptionKey const grid_nsamp_radius( "matdes:design:grid_nsamp_radius" );  } }
namespace matdes { namespace design { RealOptionKey const fav_nat_bonus( "matdes:design:fav_nat_bonus" );  } }
namespace matdes { namespace mutalyze { BooleanOptionKey const mutalyze( "matdes:mutalyze" );  } }
namespace matdes { namespace mutalyze { BooleanOptionKey const calc_rot_boltz( "matdes:mutalyze:calc_rot_boltz" );  } }
namespace matdes { namespace mutalyze { BooleanOptionKey const ala_scan( "matdes:mutalyze:ala_scan" );  } }
namespace matdes { namespace mutalyze { BooleanOptionKey const revert_scan( "matdes:mutalyze:revert_scan" );  } }
namespace matdes { namespace mutalyze { BooleanOptionKey const min_rb( "matdes:mutalyze:min_rb" );  } }
namespace gpu { BooleanOptionKey const gpu( "gpu" );  }
namespace gpu { IntegerOptionKey const device( "gpu:device" );  }
namespace gpu { IntegerOptionKey const threads( "gpu:threads" );  }
namespace pb_potential { BooleanOptionKey const pb_potential( "pb_potential" );  }
namespace pb_potential { IntegerVectorOptionKey const charged_chains( "pb_potential:charged_chains" );  }
namespace pb_potential { BooleanOptionKey const sidechain_only( "pb_potential:sidechain_only" );  }
namespace pb_potential { IntegerVectorOptionKey const revamp_near_chain( "pb_potential:revamp_near_chain" );  }
namespace pb_potential { StringOptionKey const apbs_path( "pb_potential:apbs_path" );  }
namespace pb_potential { RealOptionKey const potential_cap( "pb_potential:potential_cap" );  }
namespace pb_potential { RealOptionKey const epsilon( "pb_potential:epsilon" );  }
namespace pb_potential { IntegerOptionKey const apbs_debug( "pb_potential:apbs_debug" );  }
namespace pb_potential { BooleanOptionKey const calcenergy( "pb_potential:calcenergy" );  }
