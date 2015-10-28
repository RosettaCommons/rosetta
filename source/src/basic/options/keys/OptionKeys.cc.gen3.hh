namespace constraints { BooleanOptionKey const epr_distance( "constraints:epr_distance" );  }
namespace constraints { IntegerOptionKey const combine( "constraints:combine" );  }
namespace constraints { FileOptionKey const combine_exclude_region( "constraints:combine_exclude_region" );  }
namespace constraints { BooleanOptionKey const skip_redundant( "constraints:skip_redundant" );  }
namespace constraints { IntegerOptionKey const skip_redundant_width( "constraints:skip_redundant_width" );  }
namespace corrections { BooleanOptionKey const corrections( "corrections" );  }
namespace corrections { BooleanOptionKey const beta_july15( "corrections:beta_july15" );  }
namespace corrections { BooleanOptionKey const beta_july15_cart( "corrections:beta_july15_cart" );  }
namespace corrections { BooleanOptionKey const newdna( "corrections:newdna" );  }
namespace corrections { BooleanOptionKey const correct( "corrections:correct" );  }
namespace corrections { BooleanOptionKey const hbond_sp2_correction( "corrections:hbond_sp2_correction" );  }
namespace corrections { BooleanOptionKey const facts_default( "corrections:facts_default" );  }
namespace corrections { namespace score { BooleanOptionKey const score( "corrections:score" );  } }
namespace corrections { namespace score { BooleanOptionKey const bbdep_omega( "corrections:score:bbdep_omega" );  } }
namespace corrections { namespace score { BooleanOptionKey const bbdep_bond_params( "corrections:score:bbdep_bond_params" );  } }
namespace corrections { namespace score { BooleanOptionKey const bbdep_bond_devs( "corrections:score:bbdep_bond_devs" );  } }
namespace corrections { namespace score { BooleanOptionKey const no_his_his_pairE( "corrections:score:no_his_his_pairE" );  } }
namespace corrections { namespace score { BooleanOptionKey const no_his_DE_pairE( "corrections:score:no_his_DE_pairE" );  } }
namespace corrections { namespace score { StringOptionKey const p_aa_pp( "corrections:score:p_aa_pp" );  } }
namespace corrections { namespace score { BooleanOptionKey const p_aa_pp_nogridshift( "corrections:score:p_aa_pp_nogridshift" );  } }
namespace corrections { namespace score { BooleanOptionKey const rama_not_squared( "corrections:score:rama_not_squared" );  } }
namespace corrections { namespace score { FileOptionKey const rama_map( "corrections:score:rama_map" );  } }
namespace corrections { namespace score { BooleanOptionKey const cenrot( "corrections:score:cenrot" );  } }
namespace corrections { namespace score { BooleanOptionKey const dun10( "corrections:score:dun10" );  } }
namespace corrections { namespace score { StringOptionKey const dun10_dir( "corrections:score:dun10_dir" );  } }
namespace corrections { namespace score { StringOptionKey const dun02_file( "corrections:score:dun02_file" );  } }
namespace corrections { namespace score { StringOptionKey const ch_o_bond_potential( "corrections:score:ch_o_bond_potential" );  } }
namespace corrections { namespace score { RealOptionKey const lj_hbond_hdis( "corrections:score:lj_hbond_hdis" );  } }
namespace corrections { namespace score { RealOptionKey const lj_hbond_OH_donor_dis( "corrections:score:lj_hbond_OH_donor_dis" );  } }
namespace corrections { namespace score { BooleanOptionKey const score12prime( "corrections:score:score12prime" );  } }
namespace corrections { namespace score { BooleanOptionKey const talaris2014( "corrections:score:talaris2014" );  } }
namespace corrections { namespace score { RealOptionKey const hbond_energy_shift( "corrections:score:hbond_energy_shift" );  } }
namespace corrections { namespace score { RealOptionKey const hb_sp2_BAH180_rise( "corrections:score:hb_sp2_BAH180_rise" );  } }
namespace corrections { namespace score { RealOptionKey const hb_sp2_outer_width( "corrections:score:hb_sp2_outer_width" );  } }
namespace corrections { namespace score { BooleanOptionKey const hb_sp2_chipen( "corrections:score:hb_sp2_chipen" );  } }
namespace corrections { namespace score { BooleanOptionKey const hbond_measure_sp3acc_BAH_from_hvy( "corrections:score:hbond_measure_sp3acc_BAH_from_hvy" );  } }
namespace corrections { namespace score { BooleanOptionKey const hb_fade_energy( "corrections:score:hb_fade_energy" );  } }
namespace corrections { namespace score { BooleanOptionKey const use_bicubic_interpolation( "corrections:score:use_bicubic_interpolation" );  } }
namespace corrections { namespace score { BooleanOptionKey const dun_normsd( "corrections:score:dun_normsd" );  } }
namespace corrections { namespace score { BooleanOptionKey const dun_entropy_correction( "corrections:score:dun_entropy_correction" );  } }
namespace corrections { namespace chemical { BooleanOptionKey const chemical( "corrections:chemical" );  } }
namespace corrections { namespace chemical { BooleanOptionKey const icoor_05_2009( "corrections:chemical:icoor_05_2009" );  } }
namespace corrections { namespace chemical { BooleanOptionKey const parse_charge( "corrections:chemical:parse_charge" );  } }
namespace corrections { namespace chemical { BooleanOptionKey const expand_st_chi2sampling( "corrections:chemical:expand_st_chi2sampling" );  } }
namespace corrections { BooleanOptionKey const shapovalov_lib_fixes_enable( "corrections:shapovalov_lib_fixes_enable" );  }
namespace corrections { namespace shapovalov_lib { BooleanOptionKey const shapovalov_lib( "corrections:shapovalov_lib" );  } }
namespace corrections { namespace shapovalov_lib { BooleanOptionKey const shap_dun10_enable( "corrections:shapovalov_lib:shap_dun10_enable" );  } }
namespace corrections { namespace shapovalov_lib { StringOptionKey const shap_dun10_smooth_level( "corrections:shapovalov_lib:shap_dun10_smooth_level" );  } }
namespace corrections { namespace shapovalov_lib { StringOptionKey const shap_dun10_dir( "corrections:shapovalov_lib:shap_dun10_dir" );  } }
namespace corrections { namespace shapovalov_lib { BooleanOptionKey const shap_dun10_use_minus_log_P_ignore_P( "corrections:shapovalov_lib:shap_dun10_use_minus_log_P_ignore_P" );  } }
namespace corrections { namespace shapovalov_lib { BooleanOptionKey const shap_rama_enable( "corrections:shapovalov_lib:shap_rama_enable" );  } }
namespace corrections { namespace shapovalov_lib { StringOptionKey const shap_rama_smooth_level( "corrections:shapovalov_lib:shap_rama_smooth_level" );  } }
namespace corrections { namespace shapovalov_lib { FileOptionKey const shap_rama_map( "corrections:shapovalov_lib:shap_rama_map" );  } }
namespace corrections { namespace shapovalov_lib { BooleanOptionKey const shap_rama_nogridshift( "corrections:shapovalov_lib:shap_rama_nogridshift" );  } }
namespace corrections { namespace shapovalov_lib { BooleanOptionKey const shap_p_aa_pp_enable( "corrections:shapovalov_lib:shap_p_aa_pp_enable" );  } }
namespace corrections { namespace shapovalov_lib { StringOptionKey const shap_p_aa_pp_smooth_level( "corrections:shapovalov_lib:shap_p_aa_pp_smooth_level" );  } }
namespace corrections { namespace shapovalov_lib { StringOptionKey const shap_p_aa_pp( "corrections:shapovalov_lib:shap_p_aa_pp" );  } }
namespace evaluation { BooleanOptionKey const evaluation( "evaluation" );  }
namespace evaluation { FileVectorOptionKey const rmsd_target( "evaluation:rmsd_target" );  }
namespace evaluation { StringVectorOptionKey const rmsd_column( "evaluation:rmsd_column" );  }
namespace evaluation { FileVectorOptionKey const rmsd_select( "evaluation:rmsd_select" );  }
namespace evaluation { FileVectorOptionKey const align_rmsd_target( "evaluation:align_rmsd_target" );  }
namespace evaluation { FileVectorOptionKey const structural_similarity( "evaluation:structural_similarity" );  }
namespace evaluation { BooleanOptionKey const contact_map( "evaluation:contact_map" );  }
namespace evaluation { StringVectorOptionKey const jscore_evaluator( "evaluation:jscore_evaluator" );  }
namespace evaluation { StringVectorOptionKey const align_rmsd_column( "evaluation:align_rmsd_column" );  }
namespace evaluation { FileVectorOptionKey const align_rmsd_fns( "evaluation:align_rmsd_fns" );  }
namespace evaluation { StringOptionKey const align_rmsd_format( "evaluation:align_rmsd_format" );  }
namespace evaluation { StringOptionKey const predicted_burial_fn( "evaluation:predicted_burial_fn" );  }
namespace evaluation { FileOptionKey const pool( "evaluation:pool" );  }
namespace evaluation { FileVectorOptionKey const rmsd( "evaluation:rmsd" );  }
namespace evaluation { FileVectorOptionKey const chirmsd( "evaluation:chirmsd" );  }
namespace evaluation { BooleanOptionKey const gdtmm( "evaluation:gdtmm" );  }
namespace evaluation { BooleanOptionKey const gdttm( "evaluation:gdttm" );  }
namespace evaluation { BooleanOptionKey const score_with_rmsd( "evaluation:score_with_rmsd" );  }
namespace evaluation { FileVectorOptionKey const constraints( "evaluation:constraints" );  }
namespace evaluation { FileVectorOptionKey const constraints_column( "evaluation:constraints_column" );  }
namespace evaluation { FileVectorOptionKey const combined_constraints( "evaluation:combined_constraints" );  }
namespace evaluation { FileVectorOptionKey const combined_constraints_column( "evaluation:combined_constraints_column" );  }
namespace evaluation { IntegerOptionKey const combine_statistics( "evaluation:combine_statistics" );  }
namespace evaluation { StringVectorOptionKey const chemical_shifts( "evaluation:chemical_shifts" );  }
namespace evaluation { StringOptionKey const sparta_dir( "evaluation:sparta_dir" );  }
namespace evaluation { StringVectorOptionKey const cam_shifts( "evaluation:cam_shifts" );  }
namespace evaluation { StringVectorOptionKey const pales( "evaluation:pales" );  }
namespace evaluation { FileVectorOptionKey const extra_score( "evaluation:extra_score" );  }
namespace evaluation { FileVectorOptionKey const extra_score_patch( "evaluation:extra_score_patch" );  }
namespace evaluation { StringVectorOptionKey const extra_score_column( "evaluation:extra_score_column" );  }
namespace evaluation { FileVectorOptionKey const extra_score_select( "evaluation:extra_score_select" );  }
namespace evaluation { FileVectorOptionKey const rdc_target( "evaluation:rdc_target" );  }
namespace evaluation { BooleanOptionKey const symmetric_rmsd( "evaluation:symmetric_rmsd" );  }
namespace evaluation { StringVectorOptionKey const rdc_column( "evaluation:rdc_column" );  }
namespace evaluation { StringVectorOptionKey const rdc( "evaluation:rdc" );  }
namespace evaluation { StringOptionKey const built_in_rdc( "evaluation:built_in_rdc" );  }
namespace evaluation { BooleanOptionKey const jump_nr( "evaluation:jump_nr" );  }
namespace evaluation { IntegerVectorOptionKey const score_exclude_res( "evaluation:score_exclude_res" );  }
namespace evaluation { IntegerOptionKey const score_sscore_short_helix( "evaluation:score_sscore_short_helix" );  }
namespace evaluation { IntegerOptionKey const score_sscore_maxloop( "evaluation:score_sscore_maxloop" );  }
namespace evaluation { IntegerOptionKey const window_size( "evaluation:window_size" );  }
namespace evaluation { StringOptionKey const I_sc( "evaluation:I_sc" );  }
namespace evaluation { BooleanOptionKey const Irms( "evaluation:Irms" );  }
namespace evaluation { BooleanOptionKey const Ca_Irms( "evaluation:Ca_Irms" );  }
namespace evaluation { BooleanOptionKey const Fnat( "evaluation:Fnat" );  }
namespace evaluation { BooleanOptionKey const Lrmsd( "evaluation:Lrmsd" );  }
namespace evaluation { BooleanOptionKey const Fnonnat( "evaluation:Fnonnat" );  }
namespace evaluation { BooleanOptionKey const DockMetrics( "evaluation:DockMetrics" );  }
namespace filters { BooleanOptionKey const filters( "filters" );  }
namespace filters { BooleanOptionKey const disable_all_filters( "filters:disable_all_filters" );  }
namespace filters { BooleanOptionKey const disable_rg_filter( "filters:disable_rg_filter" );  }
namespace filters { BooleanOptionKey const disable_co_filter( "filters:disable_co_filter" );  }
namespace filters { BooleanOptionKey const disable_sheet_filter( "filters:disable_sheet_filter" );  }
namespace filters { RealOptionKey const set_pddf_filter( "filters:set_pddf_filter" );  }
namespace filters { RealOptionKey const set_saxs_filter( "filters:set_saxs_filter" );  }
namespace gpu { BooleanOptionKey const gpu( "gpu" );  }
namespace gpu { IntegerOptionKey const device( "gpu:device" );  }
namespace gpu { IntegerOptionKey const threads( "gpu:threads" );  }
namespace jumps { BooleanOptionKey const jumps( "jumps" );  }
namespace jumps { BooleanOptionKey const evaluate( "jumps:evaluate" );  }
namespace jumps { FileOptionKey const extra_frags_for_ss( "jumps:extra_frags_for_ss" );  }
namespace jumps { BooleanOptionKey const fix_chainbreak( "jumps:fix_chainbreak" );  }
namespace jumps { FileOptionKey const fix_jumps( "jumps:fix_jumps" );  }
namespace jumps { FileOptionKey const jump_lib( "jumps:jump_lib" );  }
namespace jumps { FileOptionKey const loop_definition_from_file( "jumps:loop_definition_from_file" );  }
namespace jumps { BooleanOptionKey const no_chainbreak_in_relax( "jumps:no_chainbreak_in_relax" );  }
namespace jumps { FileOptionKey const pairing_file( "jumps:pairing_file" );  }
namespace jumps { IntegerVectorOptionKey const random_sheets( "jumps:random_sheets" );  }
namespace jumps { FileOptionKey const residue_pair_jump_file( "jumps:residue_pair_jump_file" );  }
namespace jumps { IntegerVectorOptionKey const sheets( "jumps:sheets" );  }
namespace jumps { FileOptionKey const topology_file( "jumps:topology_file" );  }
namespace jumps { BooleanOptionKey const bb_moves( "jumps:bb_moves" );  }
namespace jumps { BooleanOptionKey const no_wobble( "jumps:no_wobble" );  }
namespace jumps { BooleanOptionKey const no_shear( "jumps:no_shear" );  }
namespace jumps { BooleanOptionKey const no_sample_ss_jumps( "jumps:no_sample_ss_jumps" );  }
namespace jumps { IntegerOptionKey const invrate_jump_move( "jumps:invrate_jump_move" );  }
namespace jumps { RealOptionKey const chainbreak_weight_stage1( "jumps:chainbreak_weight_stage1" );  }
namespace jumps { RealOptionKey const chainbreak_weight_stage2( "jumps:chainbreak_weight_stage2" );  }
namespace jumps { RealOptionKey const chainbreak_weight_stage3( "jumps:chainbreak_weight_stage3" );  }
namespace jumps { RealOptionKey const chainbreak_weight_stage4( "jumps:chainbreak_weight_stage4" );  }
namespace jumps { BooleanOptionKey const ramp_chainbreaks( "jumps:ramp_chainbreaks" );  }
namespace jumps { RealOptionKey const increase_chainbreak( "jumps:increase_chainbreak" );  }
namespace jumps { BooleanOptionKey const overlap_chainbreak( "jumps:overlap_chainbreak" );  }
namespace jumps { RealOptionKey const sep_switch_accelerate( "jumps:sep_switch_accelerate" );  }
namespace jumps { BooleanOptionKey const dump_frags( "jumps:dump_frags" );  }
namespace jumps { IntegerOptionKey const njumps( "jumps:njumps" );  }
namespace jumps { IntegerOptionKey const max_strand_gap_allowed( "jumps:max_strand_gap_allowed" );  }
namespace jumps { RealOptionKey const contact_score( "jumps:contact_score" );  }
namespace jumps { BooleanOptionKey const filter_templates( "jumps:filter_templates" );  }
namespace loops { BooleanOptionKey const loops( "loops" );  }
namespace loops { StringOptionKey const cen_weights( "loops:cen_weights" );  }
namespace loops { StringOptionKey const cen_patch( "loops:cen_patch" );  }
namespace loops { FileOptionKey const input_pdb( "loops:input_pdb" );  }
namespace loops { StringVectorOptionKey const loop_file( "loops:loop_file" );  }
namespace loops { FileOptionKey const extended_loop_file( "loops:extended_loop_file" );  }
namespace loops { FileOptionKey const mm_loop_file( "loops:mm_loop_file" );  }
namespace loops { BooleanOptionKey const fix_natsc( "loops:fix_natsc" );  }
namespace loops { BooleanOptionKey const refine_only( "loops:refine_only" );  }
namespace loops { BooleanOptionKey const fast( "loops:fast" );  }
namespace loops { FileOptionKey const vall_file( "loops:vall_file" );  }
namespace loops { IntegerVectorOptionKey const frag_sizes( "loops:frag_sizes" );  }
namespace loops { FileVectorOptionKey const frag_files( "loops:frag_files" );  }
namespace loops { FileOptionKey const output_pdb( "loops:output_pdb" );  }
namespace loops { BooleanOptionKey const debug( "loops:debug" );  }
namespace loops { BooleanOptionKey const build_initial( "loops:build_initial" );  }
namespace loops { BooleanOptionKey const extended( "loops:extended" );  }
namespace loops { BooleanOptionKey const remove_extended_loops( "loops:remove_extended_loops" );  }
namespace loops { BooleanOptionKey const idealize_after_loop_close( "loops:idealize_after_loop_close" );  }
namespace loops { BooleanOptionKey const idealize_before_loop_close( "loops:idealize_before_loop_close" );  }
namespace loops { IntegerOptionKey const select_best_loop_from( "loops:select_best_loop_from" );  }
namespace loops { IntegerOptionKey const build_attempts( "loops:build_attempts" );  }
namespace loops { IntegerOptionKey const grow_attempts( "loops:grow_attempts" );  }
namespace loops { RealOptionKey const random_grow_loops_by( "loops:random_grow_loops_by" );  }
namespace loops { BooleanOptionKey const accept_aborted_loops( "loops:accept_aborted_loops" );  }
namespace loops { BooleanOptionKey const strict_loops( "loops:strict_loops" );  }
namespace loops { BooleanOptionKey const superimpose_native( "loops:superimpose_native" );  }
namespace loops { IntegerVectorOptionKey const build_specific_loops( "loops:build_specific_loops" );  }
namespace loops { BooleanOptionKey const random_order( "loops:random_order" );  }
namespace loops { BooleanOptionKey const build_all_loops( "loops:build_all_loops" );  }
namespace loops { BooleanOptionKey const fa_closure_protocol( "loops:fa_closure_protocol" );  }
namespace loops { RealOptionKey const combine_rate( "loops:combine_rate" );  }
namespace loops { BooleanOptionKey const skip_initial_loop_build( "loops:skip_initial_loop_build" );  }
namespace loops { StringOptionKey const remodel( "loops:remodel" );  }
namespace loops { StringOptionKey const intermedrelax( "loops:intermedrelax" );  }
namespace loops { StringOptionKey const refine( "loops:refine" );  }
namespace loops { StringOptionKey const relax( "loops:relax" );  }
namespace loops { BooleanOptionKey const final_clean_fastrelax( "loops:final_clean_fastrelax" );  }
namespace loops { BooleanOptionKey const relax_with_foldtree( "loops:relax_with_foldtree" );  }
namespace loops { RealOptionKey const constrain_rigid_segments( "loops:constrain_rigid_segments" );  }
namespace loops { StringOptionKey const loopscores( "loops:loopscores" );  }
namespace loops { BooleanOptionKey const timer( "loops:timer" );  }
namespace loops { BooleanOptionKey const vicinity_sampling( "loops:vicinity_sampling" );  }
namespace loops { RealOptionKey const vicinity_degree( "loops:vicinity_degree" );  }
namespace loops { RealOptionKey const neighbor_dist( "loops:neighbor_dist" );  }
namespace loops { IntegerOptionKey const kic_max_seglen( "loops:kic_max_seglen" );  }
namespace loops { BooleanOptionKey const kic_recover_last( "loops:kic_recover_last" );  }
namespace loops { BooleanOptionKey const kic_min_after_repack( "loops:kic_min_after_repack" );  }
namespace loops { BooleanOptionKey const optimize_only_kic_region_sidechains_after_move( "loops:optimize_only_kic_region_sidechains_after_move" );  }
namespace loops { IntegerOptionKey const max_kic_build_attempts( "loops:max_kic_build_attempts" );  }
namespace loops { IntegerOptionKey const remodel_kic_attempts( "loops:remodel_kic_attempts" );  }
namespace loops { IntegerOptionKey const max_kic_perturber_samples( "loops:max_kic_perturber_samples" );  }
namespace loops { BooleanOptionKey const nonpivot_torsion_sampling( "loops:nonpivot_torsion_sampling" );  }
namespace loops { BooleanOptionKey const fix_ca_bond_angles( "loops:fix_ca_bond_angles" );  }
namespace loops { BooleanOptionKey const kic_use_linear_chainbreak( "loops:kic_use_linear_chainbreak" );  }
namespace loops { BooleanOptionKey const sample_omega_at_pre_prolines( "loops:sample_omega_at_pre_prolines" );  }
namespace loops { BooleanOptionKey const allow_omega_move( "loops:allow_omega_move" );  }
namespace loops { BooleanOptionKey const kic_with_cartmin( "loops:kic_with_cartmin" );  }
namespace loops { BooleanOptionKey const allow_takeoff_torsion_move( "loops:allow_takeoff_torsion_move" );  }
namespace loops { IntegerOptionKey const extend_length( "loops:extend_length" );  }
namespace loops { IntegerOptionKey const perturb_outer_cycles( "loops:perturb_outer_cycles" );  }
namespace loops { IntegerOptionKey const refine_outer_cycles( "loops:refine_outer_cycles" );  }
namespace loops { IntegerOptionKey const max_inner_cycles( "loops:max_inner_cycles" );  }
namespace loops { IntegerOptionKey const repack_period( "loops:repack_period" );  }
namespace loops { RealOptionKey const remodel_init_temp( "loops:remodel_init_temp" );  }
namespace loops { RealOptionKey const remodel_final_temp( "loops:remodel_final_temp" );  }
namespace loops { RealOptionKey const refine_init_temp( "loops:refine_init_temp" );  }
namespace loops { RealOptionKey const refine_final_temp( "loops:refine_final_temp" );  }
namespace loops { IntegerOptionKey const gapspan( "loops:gapspan" );  }
namespace loops { IntegerOptionKey const spread( "loops:spread" );  }
namespace loops { IntegerOptionKey const kinematic_wrapper_cycles( "loops:kinematic_wrapper_cycles" );  }
namespace loops { IntegerOptionKey const kic_num_rotamer_trials( "loops:kic_num_rotamer_trials" );  }
namespace loops { BooleanOptionKey const kic_omega_sampling( "loops:kic_omega_sampling" );  }
