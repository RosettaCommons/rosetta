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
namespace templates { BooleanOptionKey const templates( "templates" );  }
namespace templates { FileOptionKey const config( "templates:config" );  }
namespace templates { BooleanOptionKey const fix_aligned_residues( "templates:fix_aligned_residues" );  }
namespace templates { FileOptionKey const fix_frag_file( "templates:fix_frag_file" );  }
namespace templates { IntegerOptionKey const fix_margin( "templates:fix_margin" );  }
namespace templates { IntegerOptionKey const min_nr_large_frags( "templates:min_nr_large_frags" );  }
namespace templates { IntegerOptionKey const min_nr_small_frags( "templates:min_nr_small_frags" );  }
namespace templates { BooleanOptionKey const no_pick_fragments( "templates:no_pick_fragments" );  }
namespace templates { IntegerOptionKey const nr_large_copies( "templates:nr_large_copies" );  }
namespace templates { IntegerOptionKey const nr_small_copies( "templates:nr_small_copies" );  }
namespace templates { BooleanOptionKey const pairings( "templates:pairings" );  }
namespace templates { BooleanOptionKey const pick_multiple_sizes( "templates:pick_multiple_sizes" );  }
namespace templates { BooleanOptionKey const strand_constraint( "templates:strand_constraint" );  }
namespace templates { BooleanOptionKey const vary_frag_size( "templates:vary_frag_size" );  }
namespace templates { BooleanOptionKey const no_culling( "templates:no_culling" );  }
namespace templates { FileOptionKey const helix_pairings( "templates:helix_pairings" );  }
namespace templates { FileOptionKey const prefix( "templates:prefix" );  }
namespace templates { IntegerOptionKey const change_movemap( "templates:change_movemap" );  }
namespace templates { BooleanOptionKey const force_native_topology( "templates:force_native_topology" );  }
namespace templates { RealOptionKey const topology_rank_cutoff( "templates:topology_rank_cutoff" );  }
namespace templates { IntegerOptionKey const min_frag_size( "templates:min_frag_size" );  }
namespace templates { IntegerOptionKey const max_shrink( "templates:max_shrink" );  }
namespace templates { IntegerOptionKey const shrink_step( "templates:shrink_step" );  }
namespace templates { IntegerOptionKey const shrink_pos_step( "templates:shrink_pos_step" );  }
namespace templates { IntegerOptionKey const min_padding( "templates:min_padding" );  }
namespace templates { IntegerOptionKey const min_align_pos( "templates:min_align_pos" );  }
namespace templates { IntegerOptionKey const max_align_pos( "templates:max_align_pos" );  }
namespace templates { namespace cst { BooleanOptionKey const cst( "templates:cst" );  } }
namespace templates { namespace cst { IntegerOptionKey const topN( "templates:cst:topN" );  } }
namespace templates { namespace cst { RealOptionKey const wTopol( "templates:cst:wTopol" );  } }
namespace templates { namespace cst { RealOptionKey const wExtern( "templates:cst:wExtern" );  } }
namespace templates { namespace fragsteal { BooleanOptionKey const fragsteal( "templates:fragsteal" );  } }
namespace templates { namespace fragsteal { IntegerOptionKey const topN( "templates:fragsteal:topN" );  } }
namespace templates { namespace fragsteal { RealOptionKey const wTopol( "templates:fragsteal:wTopol" );  } }
namespace templates { namespace fragsteal { RealOptionKey const wExtern( "templates:fragsteal:wExtern" );  } }
namespace abrelax { BooleanOptionKey const abrelax( "abrelax" );  }
namespace abrelax { BooleanOptionKey const filters( "abrelax:filters" );  }
namespace abrelax { BooleanOptionKey const fail_unclosed( "abrelax:fail_unclosed" );  }
namespace chemical { BooleanOptionKey const chemical( "chemical" );  }
namespace chemical { StringVectorOptionKey const exclude_patches( "chemical:exclude_patches" );  }
namespace chemical { StringVectorOptionKey const include_patches( "chemical:include_patches" );  }
namespace chemical { BooleanOptionKey const enlarge_H_lj( "chemical:enlarge_H_lj" );  }
namespace chemical { StringVectorOptionKey const add_atom_type_set_parameters( "chemical:add_atom_type_set_parameters" );  }
namespace chemical { StringVectorOptionKey const set_atom_properties( "chemical:set_atom_properties" );  }
namespace score { BooleanOptionKey const score_pose_cutpoint_variants( "score:score_pose_cutpoint_variants" );  }
namespace score { BooleanOptionKey const score( "score" );  }
namespace score { StringOptionKey const weights( "score:weights" );  }
namespace score { StringVectorOptionKey const set_weights( "score:set_weights" );  }
namespace score { StringOptionKey const pack_weights( "score:pack_weights" );  }
namespace score { StringOptionKey const soft_wts( "score:soft_wts" );  }
namespace score { BooleanOptionKey const docking_interface_score( "score:docking_interface_score" );  }
namespace score { RealOptionKey const min_score_score( "score:min_score_score" );  }
namespace score { StringOptionKey const custom_atom_pair( "score:custom_atom_pair" );  }
namespace score { FileVectorOptionKey const patch( "score:patch" );  }
namespace score { BooleanOptionKey const empty( "score:empty" );  }
namespace score { RealOptionKey const fa_max_dis( "score:fa_max_dis" );  }
namespace score { BooleanOptionKey const fa_Hatr( "score:fa_Hatr" );  }
namespace score { BooleanOptionKey const no_smooth_etables( "score:no_smooth_etables" );  }
namespace score { RealOptionKey const etable_lr( "score:etable_lr" );  }
namespace score { BooleanOptionKey const no_lk_polar_desolvation( "score:no_lk_polar_desolvation" );  }
namespace score { StringOptionKey const input_etables( "score:input_etables" );  }
namespace score { StringOptionKey const output_etables( "score:output_etables" );  }
namespace score { BooleanOptionKey const analytic_etable_evaluation( "score:analytic_etable_evaluation" );  }
namespace score { RealOptionKey const rms_target( "score:rms_target" );  }
namespace score { BooleanOptionKey const ramaneighbors( "score:ramaneighbors" );  }
namespace score { StringOptionKey const optH_weights( "score:optH_weights" );  }
namespace score { StringOptionKey const optH_patch( "score:optH_patch" );  }
namespace score { StringOptionKey const hbond_params( "score:hbond_params" );  }
namespace score { BooleanOptionKey const hbond_disable_bbsc_exclusion_rule( "score:hbond_disable_bbsc_exclusion_rule" );  }
namespace score { IntegerOptionKey const symE_units( "score:symE_units" );  }
namespace score { RealOptionKey const symE_bonus( "score:symE_bonus" );  }
namespace score { RealOptionKey const NV_lbound( "score:NV_lbound" );  }
namespace score { RealOptionKey const NV_ubound( "score:NV_ubound" );  }
namespace score { StringOptionKey const NV_table( "score:NV_table" );  }
namespace score { BooleanOptionKey const disable_orientation_dependent_rna_ch_o_bonds( "score:disable_orientation_dependent_rna_ch_o_bonds" );  }
namespace score { StringOptionKey const rna_torsion_potential( "score:rna_torsion_potential" );  }
namespace score { BooleanOptionKey const rna_torsion_skip_chainbreak( "score:rna_torsion_skip_chainbreak" );  }
namespace score { StringOptionKey const rna_chemical_shift_exp_data( "score:rna_chemical_shift_exp_data" );  }
namespace score { StringOptionKey const rna_chemical_shift_H5_prime_mode( "score:rna_chemical_shift_H5_prime_mode" );  }
namespace score { IntegerVectorOptionKey const rna_chemical_shift_include_res( "score:rna_chemical_shift_include_res" );  }
namespace score { BooleanOptionKey const use_2prime_OH_potential( "score:use_2prime_OH_potential" );  }
namespace score { BooleanOptionKey const include_neighbor_base_stacks( "score:include_neighbor_base_stacks" );  }
namespace score { BooleanOptionKey const find_neighbors_3dgrid( "score:find_neighbors_3dgrid" );  }
namespace score { BooleanOptionKey const find_neighbors_stripehash( "score:find_neighbors_stripehash" );  }
namespace score { StringOptionKey const seqdep_refene_fname( "score:seqdep_refene_fname" );  }
namespace score { StringOptionKey const secondary_seqdep_refene_fname( "score:secondary_seqdep_refene_fname" );  }
namespace score { BooleanOptionKey const exact_occ_pairwise( "score:exact_occ_pairwise" );  }
namespace score { BooleanOptionKey const exact_occ_skip_Hbonders( "score:exact_occ_skip_Hbonders" );  }
namespace score { BooleanOptionKey const exact_occ_include_Hbond_contribution( "score:exact_occ_include_Hbond_contribution" );  }
namespace score { BooleanOptionKey const exact_occ_pairwise_by_res( "score:exact_occ_pairwise_by_res" );  }
namespace score { BooleanOptionKey const exact_occ_split_between_res( "score:exact_occ_split_between_res" );  }
namespace score { BooleanOptionKey const exact_occ_self_res_no_occ( "score:exact_occ_self_res_no_occ" );  }
namespace score { RealOptionKey const exact_occ_radius_scaling( "score:exact_occ_radius_scaling" );  }
namespace score { StringVectorOptionKey const ref_offsets( "score:ref_offsets" );  }
namespace score { BooleanOptionKey const output_residue_energies( "score:output_residue_energies" );  }
namespace score { StringOptionKey const fa_custom_pair_distance_file( "score:fa_custom_pair_distance_file" );  }
namespace score { RealOptionKey const disulf_matching_probe( "score:disulf_matching_probe" );  }
namespace score { RealVectorOptionKey const bonded_params( "score:bonded_params" );  }
namespace score { StringOptionKey const bonded_params_dir( "score:bonded_params_dir" );  }
namespace score { StringOptionKey const extra_improper_file( "score:extra_improper_file" );  }
namespace score { RealOptionKey const pro_close_planar_constraint( "score:pro_close_planar_constraint" );  }
namespace score { BooleanOptionKey const linear_bonded_potential( "score:linear_bonded_potential" );  }
namespace score { BooleanOptionKey const geom_sol_correct_acceptor_base( "score:geom_sol_correct_acceptor_base" );  }
namespace score { IntegerVectorOptionKey const rg_local_span( "score:rg_local_span" );  }
namespace score { BooleanOptionKey const unmodifypot( "score:unmodifypot" );  }
namespace score { namespace saxs { BooleanOptionKey const saxs( "score:saxs" );  } }
namespace score { namespace saxs { RealOptionKey const min_score( "score:saxs:min_score" );  } }
namespace score { namespace saxs { StringOptionKey const custom_ff( "score:saxs:custom_ff" );  } }
namespace score { namespace saxs { StringOptionKey const print_i_calc( "score:saxs:print_i_calc" );  } }
namespace score { namespace saxs { FileOptionKey const ref_fa_spectrum( "score:saxs:ref_fa_spectrum" );  } }
namespace score { namespace saxs { FileOptionKey const ref_cen_spectrum( "score:saxs:ref_cen_spectrum" );  } }
namespace score { namespace saxs { FileOptionKey const ref_spectrum( "score:saxs:ref_spectrum" );  } }
namespace score { namespace saxs { FileOptionKey const ref_pddf( "score:saxs:ref_pddf" );  } }
namespace score { namespace saxs { BooleanOptionKey const skip_hydrogens( "score:saxs:skip_hydrogens" );  } }
namespace score { namespace saxs { RealOptionKey const d_min( "score:saxs:d_min" );  } }
namespace score { namespace saxs { RealOptionKey const d_max( "score:saxs:d_max" );  } }
namespace score { namespace saxs { RealOptionKey const d_step( "score:saxs:d_step" );  } }
namespace score { namespace saxs { RealOptionKey const q_min( "score:saxs:q_min" );  } }
namespace score { namespace saxs { RealOptionKey const q_max( "score:saxs:q_max" );  } }
namespace score { namespace saxs { RealOptionKey const q_step( "score:saxs:q_step" );  } }
namespace score { namespace saxs { BooleanOptionKey const fit_pddf_area( "score:saxs:fit_pddf_area" );  } }
namespace score { IntegerVectorOptionKey const sidechain_buried( "score:sidechain_buried" );  }
namespace score { IntegerVectorOptionKey const sidechain_exposed( "score:sidechain_exposed" );  }
namespace score { RealOptionKey const elec_min_dis( "score:elec_min_dis" );  }
namespace score { RealOptionKey const elec_max_dis( "score:elec_max_dis" );  }
namespace score { RealOptionKey const elec_die( "score:elec_die" );  }
namespace score { BooleanOptionKey const elec_r_option( "score:elec_r_option" );  }
namespace score { BooleanOptionKey const smooth_fa_elec( "score:smooth_fa_elec" );  }
namespace score { RealOptionKey const facts_GBpair_cut( "score:facts_GBpair_cut" );  }
namespace score { RealOptionKey const facts_kappa( "score:facts_kappa" );  }
namespace score { IntegerOptionKey const facts_asp_patch( "score:facts_asp_patch" );  }
namespace score { BooleanOptionKey const facts_plane_to_self( "score:facts_plane_to_self" );  }
namespace score { RealOptionKey const facts_saltbridge_correction( "score:facts_saltbridge_correction" );  }
namespace score { RealVectorOptionKey const facts_dshift( "score:facts_dshift" );  }
namespace score { RealOptionKey const facts_die( "score:facts_die" );  }
namespace score { BooleanOptionKey const facts_binding_affinity( "score:facts_binding_affinity" );  }
namespace score { BooleanOptionKey const facts_intrascale_by_level( "score:facts_intrascale_by_level" );  }
namespace score { RealVectorOptionKey const facts_intbb_elec_scale( "score:facts_intbb_elec_scale" );  }
namespace score { RealVectorOptionKey const facts_intbb_solv_scale( "score:facts_intbb_solv_scale" );  }
namespace score { RealVectorOptionKey const facts_adjbb_elec_scale( "score:facts_adjbb_elec_scale" );  }
namespace score { RealVectorOptionKey const facts_adjbb_solv_scale( "score:facts_adjbb_solv_scale" );  }
namespace score { RealVectorOptionKey const facts_intbs_elec_scale( "score:facts_intbs_elec_scale" );  }
namespace score { RealVectorOptionKey const facts_intbs_solv_scale( "score:facts_intbs_solv_scale" );  }
namespace score { RealVectorOptionKey const facts_adjbs_elec_scale( "score:facts_adjbs_elec_scale" );  }
namespace score { RealVectorOptionKey const facts_adjbs_solv_scale( "score:facts_adjbs_solv_scale" );  }
namespace score { RealVectorOptionKey const facts_intsc_elec_scale( "score:facts_intsc_elec_scale" );  }
namespace score { RealVectorOptionKey const facts_intsc_solv_scale( "score:facts_intsc_solv_scale" );  }
namespace score { StringOptionKey const facts_charge_dir( "score:facts_charge_dir" );  }
namespace score { StringOptionKey const facts_eff_charge_dir( "score:facts_eff_charge_dir" );  }
namespace score { StringVectorOptionKey const facts_plane_aa( "score:facts_plane_aa" );  }
namespace score { StringOptionKey const facts_eq_type( "score:facts_eq_type" );  }
namespace score { StringOptionKey const nmer_ref_energies( "score:nmer_ref_energies" );  }
namespace score { StringOptionKey const nmer_ref_energies_list( "score:nmer_ref_energies_list" );  }
namespace score { StringOptionKey const nmer_pssm( "score:nmer_pssm" );  }
namespace score { StringOptionKey const nmer_pssm_list( "score:nmer_pssm_list" );  }
namespace score { RealOptionKey const nmer_pssm_scorecut( "score:nmer_pssm_scorecut" );  }
namespace score { StringOptionKey const nmer_svm( "score:nmer_svm" );  }
namespace score { StringOptionKey const nmer_svm_list( "score:nmer_svm_list" );  }
namespace score { RealOptionKey const nmer_svm_scorecut( "score:nmer_svm_scorecut" );  }
namespace score { StringOptionKey const nmer_svm_aa_matrix( "score:nmer_svm_aa_matrix" );  }
namespace score { IntegerOptionKey const nmer_svm_term_length( "score:nmer_svm_term_length" );  }
namespace score { BooleanOptionKey const nmer_svm_pssm_feat( "score:nmer_svm_pssm_feat" );  }
namespace score { IntegerOptionKey const nmer_ref_seq_length( "score:nmer_ref_seq_length" );  }
namespace score { BooleanOptionKey const just_calc_rmsd( "score:just_calc_rmsd" );  }
namespace ProQ { BooleanOptionKey const ProQ( "ProQ" );  }
namespace ProQ { IntegerOptionKey const svmmodel( "ProQ:svmmodel" );  }
namespace ProQ { StringOptionKey const basename( "ProQ:basename" );  }
namespace ProQ { BooleanOptionKey const membrane( "ProQ:membrane" );  }
namespace ProQ { BooleanOptionKey const prof_bug( "ProQ:prof_bug" );  }
namespace ProQ { BooleanOptionKey const output_feature_vector( "ProQ:output_feature_vector" );  }
namespace ProQ { BooleanOptionKey const output_local_prediction( "ProQ:output_local_prediction" );  }
namespace ProQ { StringOptionKey const prefix( "ProQ:prefix" );  }
namespace ProQ { BooleanOptionKey const use_gzip( "ProQ:use_gzip" );  }
namespace ProQ { RealOptionKey const normalize( "ProQ:normalize" );  }
namespace corrections { BooleanOptionKey const corrections( "corrections" );  }
namespace corrections { BooleanOptionKey const correct( "corrections:correct" );  }
namespace corrections { BooleanOptionKey const hbond_sp2_correction( "corrections:hbond_sp2_correction" );  }
namespace corrections { BooleanOptionKey const facts_default( "corrections:facts_default" );  }
namespace corrections { namespace score { BooleanOptionKey const score( "corrections:score" );  } }
namespace corrections { namespace score { BooleanOptionKey const bbdep_omega( "corrections:score:bbdep_omega" );  } }
namespace corrections { namespace score { BooleanOptionKey const bbdep_bond_params( "corrections:score:bbdep_bond_params" );  } }
namespace corrections { namespace score { BooleanOptionKey const bbdep_bond_devs( "corrections:score:bbdep_bond_devs" );  } }
namespace corrections { namespace score { BooleanOptionKey const no_his_his_pairE( "corrections:score:no_his_his_pairE" );  } }
namespace corrections { namespace score { BooleanOptionKey const no_his_DE_pairE( "corrections:score:no_his_DE_pairE" );  } }
namespace corrections { namespace score { BooleanOptionKey const hbond_His_Phil_fix( "corrections:score:hbond_His_Phil_fix" );  } }
namespace corrections { namespace score { BooleanOptionKey const helix_hb_06_2009( "corrections:score:helix_hb_06_2009" );  } }
namespace corrections { namespace score { BooleanOptionKey const use_incorrect_hbond_deriv( "corrections:score:use_incorrect_hbond_deriv" );  } }
namespace corrections { namespace score { StringOptionKey const p_aa_pp( "corrections:score:p_aa_pp" );  } }
namespace corrections { namespace score { BooleanOptionKey const p_aa_pp_nogridshift( "corrections:score:p_aa_pp_nogridshift" );  } }
namespace corrections { namespace score { BooleanOptionKey const rama_not_squared( "corrections:score:rama_not_squared" );  } }
namespace corrections { namespace score { FileOptionKey const rama_map( "corrections:score:rama_map" );  } }
namespace corrections { namespace score { BooleanOptionKey const cenrot( "corrections:score:cenrot" );  } }
namespace corrections { namespace score { BooleanOptionKey const dun10( "corrections:score:dun10" );  } }
namespace corrections { namespace score { StringOptionKey const dun10_dir( "corrections:score:dun10_dir" );  } }
namespace corrections { namespace score { StringOptionKey const dun02_file( "corrections:score:dun02_file" );  } }
namespace corrections { namespace score { StringOptionKey const ch_o_bond_potential( "corrections:score:ch_o_bond_potential" );  } }
namespace corrections { namespace score { BooleanOptionKey const fa_elec_co_only( "corrections:score:fa_elec_co_only" );  } }
namespace corrections { namespace score { RealOptionKey const lj_hbond_hdis( "corrections:score:lj_hbond_hdis" );  } }
namespace corrections { namespace score { RealOptionKey const lj_hbond_OH_donor_dis( "corrections:score:lj_hbond_OH_donor_dis" );  } }
namespace corrections { namespace score { BooleanOptionKey const score12prime( "corrections:score:score12prime" );  } }
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
namespace mistakes { BooleanOptionKey const mistakes( "mistakes" );  }
namespace mistakes { BooleanOptionKey const restore_pre_talaris_2013_behavior( "mistakes:restore_pre_talaris_2013_behavior" );  }
namespace mistakes { namespace chemical { BooleanOptionKey const chemical( "mistakes:chemical" );  } }
namespace mistakes { namespace chemical { BooleanOptionKey const pre_talaris2013_geometries( "mistakes:chemical:pre_talaris2013_geometries" );  } }
namespace willmatch { BooleanOptionKey const willmatch( "willmatch" );  }
namespace willmatch { RealOptionKey const arg_dun_th( "willmatch:arg_dun_th" );  }
namespace willmatch { RealOptionKey const asp_dun_th( "willmatch:asp_dun_th" );  }
namespace willmatch { RealOptionKey const glu_dun_th( "willmatch:glu_dun_th" );  }
namespace willmatch { RealOptionKey const lys_dun_th( "willmatch:lys_dun_th" );  }
namespace willmatch { BooleanOptionKey const usecache( "willmatch:usecache" );  }
namespace willmatch { StringVectorOptionKey const write_reduced_matchset( "willmatch:write_reduced_matchset" );  }
namespace willmatch { RealOptionKey const interface_size( "willmatch:interface_size" );  }
namespace willmatch { RealOptionKey const max_dis_any( "willmatch:max_dis_any" );  }
namespace willmatch { RealOptionKey const max_dis_all( "willmatch:max_dis_all" );  }
namespace willmatch { RealOptionKey const max_dis_hb( "willmatch:max_dis_hb" );  }
namespace willmatch { RealOptionKey const min_dis_hb( "willmatch:min_dis_hb" );  }
namespace willmatch { RealOptionKey const max_dis_hb_colinear( "willmatch:max_dis_hb_colinear" );  }
namespace willmatch { RealOptionKey const max_dis_metal( "willmatch:max_dis_metal" );  }
namespace willmatch { RealOptionKey const max_ang_metal( "willmatch:max_ang_metal" );  }
namespace willmatch { RealOptionKey const clash_dis( "willmatch:clash_dis" );  }
namespace willmatch { RealOptionKey const c2_linker_dist( "willmatch:c2_linker_dist" );  }
namespace willmatch { RealOptionKey const identical_match_dis( "willmatch:identical_match_dis" );  }
namespace willmatch { RealOptionKey const chi1_increment( "willmatch:chi1_increment" );  }
namespace willmatch { RealOptionKey const chi2_increment( "willmatch:chi2_increment" );  }
namespace willmatch { RealOptionKey const c2_symm_increment( "willmatch:c2_symm_increment" );  }
namespace willmatch { RealOptionKey const cb_sasa_thresh( "willmatch:cb_sasa_thresh" );  }
namespace willmatch { BooleanOptionKey const design_interface( "willmatch:design_interface" );  }
namespace willmatch { FileOptionKey const chilist( "willmatch:chilist" );  }
namespace willmatch { FileOptionKey const fixed_res( "willmatch:fixed_res" );  }
namespace willmatch { FileOptionKey const native1( "willmatch:native1" );  }
namespace willmatch { FileOptionKey const native2( "willmatch:native2" );  }
namespace willmatch { FileOptionKey const exclude_res1( "willmatch:exclude_res1" );  }
namespace willmatch { FileOptionKey const exclude_res2( "willmatch:exclude_res2" );  }
namespace willmatch { FileOptionKey const taglist( "willmatch:taglist" );  }
namespace willmatch { IntegerVectorOptionKey const residues( "willmatch:residues" );  }
namespace willmatch { BooleanOptionKey const symmetry_d2( "willmatch:symmetry_d2" );  }
namespace willmatch { BooleanOptionKey const symmetry_c2_dock( "willmatch:symmetry_c2_dock" );  }
namespace willmatch { IntegerVectorOptionKey const splitwork( "willmatch:splitwork" );  }
namespace willmatch { BooleanOptionKey const exclude_ala( "willmatch:exclude_ala" );  }
namespace willmatch { RealOptionKey const match_overlap_dis( "willmatch:match_overlap_dis" );  }
namespace willmatch { RealOptionKey const match_overlap_ang( "willmatch:match_overlap_ang" );  }
namespace willmatch { IntegerVectorOptionKey const forbid_residues( "willmatch:forbid_residues" );  }
namespace willmatch { RealVectorOptionKey const poi( "willmatch:poi" );  }
namespace willmatch { RealOptionKey const poidis( "willmatch:poidis" );  }
namespace willmatch { BooleanOptionKey const homodimer( "willmatch:homodimer" );  }
namespace willmatch { RealOptionKey const fa_dun_thresh( "willmatch:fa_dun_thresh" );  }
namespace holes { BooleanOptionKey const holes( "holes" );  }
namespace holes { FileOptionKey const dalphaball( "holes:dalphaball" );  }
namespace holes { FileOptionKey const params( "holes:params" );  }
namespace holes { IntegerOptionKey const h_mode( "holes:h_mode" );  }
namespace holes { BooleanOptionKey const water( "holes:water" );  }
namespace holes { BooleanOptionKey const make_pdb( "holes:make_pdb" );  }
namespace holes { BooleanOptionKey const make_voids( "holes:make_voids" );  }
namespace holes { BooleanOptionKey const atom_scores( "holes:atom_scores" );  }
namespace holes { BooleanOptionKey const residue_scores( "holes:residue_scores" );  }
namespace holes { RealOptionKey const cav_shrink( "holes:cav_shrink" );  }
namespace holes { StringOptionKey const minimize( "holes:minimize" );  }
namespace holes { BooleanOptionKey const debug( "holes:debug" );  }
namespace packstat { BooleanOptionKey const packstat( "packstat" );  }
namespace packstat { BooleanOptionKey const include_water( "packstat:include_water" );  }
namespace packstat { IntegerOptionKey const oversample( "packstat:oversample" );  }
namespace packstat { BooleanOptionKey const packstat_pdb( "packstat:packstat_pdb" );  }
namespace packstat { BooleanOptionKey const surface_accessibility( "packstat:surface_accessibility" );  }
namespace packstat { BooleanOptionKey const residue_scores( "packstat:residue_scores" );  }
namespace packstat { RealOptionKey const cavity_burial_probe_radius( "packstat:cavity_burial_probe_radius" );  }
namespace packstat { BooleanOptionKey const raw_stats( "packstat:raw_stats" );  }
namespace packstat { IntegerOptionKey const threads( "packstat:threads" );  }
namespace packstat { RealOptionKey const cluster_min_volume( "packstat:cluster_min_volume" );  }
namespace packstat { RealOptionKey const min_surface_accessibility( "packstat:min_surface_accessibility" );  }
namespace packstat { RealOptionKey const min_cluster_overlap( "packstat:min_cluster_overlap" );  }
namespace packstat { RealOptionKey const min_cav_ball_radius( "packstat:min_cav_ball_radius" );  }
namespace packstat { RealOptionKey const max_cav_ball_radius( "packstat:max_cav_ball_radius" );  }
namespace crossmatch { BooleanOptionKey const crossmatch( "crossmatch" );  }
namespace crossmatch { StringVectorOptionKey const write_reduced_matchset( "crossmatch:write_reduced_matchset" );  }
namespace crossmatch { IntegerOptionKey const interface_size( "crossmatch:interface_size" );  }
namespace crossmatch { RealOptionKey const max_dis_any( "crossmatch:max_dis_any" );  }
namespace crossmatch { RealOptionKey const max_dis_all( "crossmatch:max_dis_all" );  }
namespace crossmatch { RealOptionKey const max_dis_metal( "crossmatch:max_dis_metal" );  }
namespace crossmatch { RealOptionKey const clash_dis( "crossmatch:clash_dis" );  }
namespace crossmatch { RealOptionKey const identical_match_dis( "crossmatch:identical_match_dis" );  }
namespace smhybrid { BooleanOptionKey const smhybrid( "smhybrid" );  }
namespace smhybrid { BooleanOptionKey const add_cavities( "smhybrid:add_cavities" );  }
namespace smhybrid { BooleanOptionKey const abinitio_design( "smhybrid:abinitio_design" );  }
namespace smhybrid { BooleanOptionKey const fa_refine( "smhybrid:fa_refine" );  }
namespace smhybrid { BooleanOptionKey const virtual_nterm( "smhybrid:virtual_nterm" );  }
namespace smhybrid { BooleanOptionKey const debug( "smhybrid:debug" );  }
namespace smhybrid { BooleanOptionKey const refine( "smhybrid:refine" );  }
namespace smhybrid { BooleanOptionKey const filter( "smhybrid:filter" );  }
namespace smhybrid { BooleanOptionKey const floating_scs_rep( "smhybrid:floating_scs_rep" );  }
namespace smhybrid { BooleanOptionKey const flxbb( "smhybrid:flxbb" );  }
namespace smhybrid { BooleanOptionKey const centroid_all_val( "smhybrid:centroid_all_val" );  }
namespace smhybrid { BooleanOptionKey const subsubs_attract( "smhybrid:subsubs_attract" );  }
namespace smhybrid { BooleanOptionKey const linker_cst( "smhybrid:linker_cst" );  }
namespace smhybrid { BooleanOptionKey const pseudosym( "smhybrid:pseudosym" );  }
namespace smhybrid { BooleanOptionKey const design_linker( "smhybrid:design_linker" );  }
namespace smhybrid { BooleanOptionKey const design( "smhybrid:design" );  }
namespace smhybrid { BooleanOptionKey const restrict_design_to_interface( "smhybrid:restrict_design_to_interface" );  }
namespace smhybrid { BooleanOptionKey const restrict_design_to_subsub_interface( "smhybrid:restrict_design_to_subsub_interface" );  }
namespace smhybrid { BooleanOptionKey const design_hydrophobic( "smhybrid:design_hydrophobic" );  }
namespace smhybrid { BooleanOptionKey const add_metal_at_0( "smhybrid:add_metal_at_0" );  }
namespace smhybrid { IntegerOptionKey const nres_mono( "smhybrid:nres_mono" );  }
namespace smhybrid { IntegerOptionKey const abinitio_cycles( "smhybrid:abinitio_cycles" );  }
namespace smhybrid { IntegerOptionKey const primary_subsubunit( "smhybrid:primary_subsubunit" );  }
namespace smhybrid { IntegerOptionKey const minbb( "smhybrid:minbb" );  }
namespace smhybrid { IntegerOptionKey const switch_concert_sub( "smhybrid:switch_concert_sub" );  }
namespace smhybrid { RealOptionKey const temperature( "smhybrid:temperature" );  }
namespace smhybrid { BooleanOptionKey const inter_subsub_cst( "smhybrid:inter_subsub_cst" );  }
namespace smhybrid { RealOptionKey const rb_mag( "smhybrid:rb_mag" );  }
namespace smhybrid { StringOptionKey const ss( "smhybrid:ss" );  }
namespace smhybrid { FileOptionKey const symm_def_template( "smhybrid:symm_def_template" );  }
namespace smhybrid { FileOptionKey const symm_def_template_reduced( "smhybrid:symm_def_template_reduced" );  }
namespace smhybrid { IntegerVectorOptionKey const attach_as_sc( "smhybrid:attach_as_sc" );  }
namespace smhybrid { IntegerVectorOptionKey const attach_as_sc_sub( "smhybrid:attach_as_sc_sub" );  }
namespace smhybrid { IntegerVectorOptionKey const inversion_subs( "smhybrid:inversion_subs" );  }
namespace smhybrid { BooleanVectorOptionKey const chainbreaks( "smhybrid:chainbreaks" );  }
namespace smhybrid { StringVectorOptionKey const design_res_files( "smhybrid:design_res_files" );  }
namespace smhybrid { StringVectorOptionKey const fixed_res_files( "smhybrid:fixed_res_files" );  }
namespace smhybrid { StringVectorOptionKey const frag_res_files( "smhybrid:frag_res_files" );  }
namespace smhybrid { StringVectorOptionKey const scattach_res_files( "smhybrid:scattach_res_files" );  }
namespace smhybrid { StringVectorOptionKey const rep_edge_files( "smhybrid:rep_edge_files" );  }
namespace smhybrid { StringVectorOptionKey const virtual_res_files( "smhybrid:virtual_res_files" );  }
namespace smhybrid { StringVectorOptionKey const jumpcut_files( "smhybrid:jumpcut_files" );  }
namespace smhybrid { StringVectorOptionKey const cst_sub_files( "smhybrid:cst_sub_files" );  }
namespace smhybrid { StringVectorOptionKey const symm_file_tag( "smhybrid:symm_file_tag" );  }
namespace smhybrid { StringVectorOptionKey const attach_atom( "smhybrid:attach_atom" );  }
namespace smhybrid { StringVectorOptionKey const add_res_before( "smhybrid:add_res_before" );  }
namespace smhybrid { StringVectorOptionKey const add_res_after( "smhybrid:add_res_after" );  }
namespace smhybrid { StringVectorOptionKey const add_ss_before( "smhybrid:add_ss_before" );  }
namespace smhybrid { StringVectorOptionKey const add_ss_after( "smhybrid:add_ss_after" );  }
namespace smhybrid { StringVectorOptionKey const add_atom_at_cen( "smhybrid:add_atom_at_cen" );  }
namespace smhybrid { StringVectorOptionKey const attach_rsd( "smhybrid:attach_rsd" );  }
namespace evolution { BooleanOptionKey const evolution( "evolution" );  }
namespace evolution { FileVectorOptionKey const parentlist( "evolution:parentlist" );  }
namespace evolution { FileVectorOptionKey const childlist( "evolution:childlist" );  }
namespace evolution { StringOptionKey const action( "evolution:action" );  }
namespace evolution { RealOptionKey const rms_threshold( "evolution:rms_threshold" );  }
namespace evolution { RealOptionKey const rms_topmargin( "evolution:rms_topmargin" );  }
namespace evolution { StringOptionKey const targetdir( "evolution:targetdir" );  }
namespace evolution { RealOptionKey const padding_score_filter( "evolution:padding_score_filter" );  }
namespace evolution { RealOptionKey const padding_stage2_filter( "evolution:padding_stage2_filter" );  }
namespace cluster { BooleanOptionKey const cluster( "cluster" );  }
namespace cluster { BooleanOptionKey const lite( "cluster:lite" );  }
namespace cluster { RealOptionKey const input_score_filter( "cluster:input_score_filter" );  }
namespace cluster { RealOptionKey const output_score_filter( "cluster:output_score_filter" );  }
namespace cluster { IntegerVectorOptionKey const exclude_res( "cluster:exclude_res" );  }
namespace cluster { RealOptionKey const thinout_factor( "cluster:thinout_factor" );  }
namespace cluster { IntegerOptionKey const max_cluster_seeds( "cluster:max_cluster_seeds" );  }
namespace cluster { RealOptionKey const radius( "cluster:radius" );  }
namespace cluster { IntegerOptionKey const limit_cluster_size( "cluster:limit_cluster_size" );  }
namespace cluster { RealOptionKey const limit_cluster_size_percent( "cluster:limit_cluster_size_percent" );  }
namespace cluster { RealOptionKey const random_limit_cluster_size_percent( "cluster:random_limit_cluster_size_percent" );  }
namespace cluster { IntegerOptionKey const limit_clusters( "cluster:limit_clusters" );  }
namespace cluster { IntegerOptionKey const limit_total_structures( "cluster:limit_total_structures" );  }
namespace cluster { IntegerOptionKey const max_total_cluster( "cluster:max_total_cluster" );  }
namespace cluster { BooleanOptionKey const gdtmm( "cluster:gdtmm" );  }
namespace cluster { BooleanOptionKey const sort_groups_by_energy( "cluster:sort_groups_by_energy" );  }
namespace cluster { BooleanOptionKey const sort_groups_by_size( "cluster:sort_groups_by_size" );  }
namespace cluster { BooleanOptionKey const remove_singletons( "cluster:remove_singletons" );  }
namespace cluster { BooleanOptionKey const export_only_low( "cluster:export_only_low" );  }
namespace cluster { BooleanOptionKey const remove_highest_energy_member( "cluster:remove_highest_energy_member" );  }
namespace cluster { BooleanOptionKey const idealize_final_structures( "cluster:idealize_final_structures" );  }
namespace cluster { IntegerOptionKey const limit_dist_matrix( "cluster:limit_dist_matrix" );  }
namespace cluster { BooleanOptionKey const make_ensemble_cst( "cluster:make_ensemble_cst" );  }
namespace cluster { BooleanOptionKey const hotspot_hash( "cluster:hotspot_hash" );  }
namespace cluster { BooleanOptionKey const loops( "cluster:loops" );  }
namespace cluster { RealOptionKey const population_weight( "cluster:population_weight" );  }
namespace cluster { StringOptionKey const template_scores( "cluster:template_scores" );  }
namespace cluster { IntegerOptionKey const K_level( "cluster:K_level" );  }
namespace cluster { RealVectorOptionKey const K_radius( "cluster:K_radius" );  }
namespace cluster { IntegerVectorOptionKey const K_n_cluster( "cluster:K_n_cluster" );  }
namespace cluster { StringVectorOptionKey const K_style( "cluster:K_style" );  }
namespace cluster { RealOptionKey const K_threshold( "cluster:K_threshold" );  }
namespace cluster { IntegerOptionKey const K_n_sub( "cluster:K_n_sub" );  }
namespace cluster { IntegerOptionKey const K_deque_size( "cluster:K_deque_size" );  }
namespace cluster { IntegerOptionKey const K_deque_level( "cluster:K_deque_level" );  }
namespace cluster { BooleanOptionKey const K_redundant( "cluster:K_redundant" );  }
namespace cluster { BooleanOptionKey const K_not_fit_xyz( "cluster:K_not_fit_xyz" );  }
namespace cluster { BooleanOptionKey const K_save_headers( "cluster:K_save_headers" );  }
namespace cluster { RealOptionKey const score_diff_cut( "cluster:score_diff_cut" );  }
namespace cluster { BooleanOptionKey const auto_tune( "cluster:auto_tune" );  }
namespace rescore { BooleanOptionKey const rescore( "rescore" );  }
namespace rescore { BooleanOptionKey const pose_metrics( "rescore:pose_metrics" );  }
namespace rescore { BooleanOptionKey const assign_ss( "rescore:assign_ss" );  }
namespace rescore { BooleanOptionKey const skip( "rescore:skip" );  }
namespace rescore { BooleanOptionKey const verbose( "rescore:verbose" );  }
namespace rescore { StringOptionKey const msms_analysis( "rescore:msms_analysis" );  }
namespace mc { BooleanOptionKey const mc( "mc" );  }
namespace mc { StringOptionKey const hierarchical_pool( "mc:hierarchical_pool" );  }
namespace mc { FileOptionKey const read_structures_into_pool( "mc:read_structures_into_pool" );  }
namespace mc { IntegerOptionKey const convergence_check_frequency( "mc:convergence_check_frequency" );  }
namespace mc { FileOptionKey const known_structures( "mc:known_structures" );  }
namespace mc { RealOptionKey const max_rmsd_against_known_structures( "mc:max_rmsd_against_known_structures" );  }
namespace mc { IntegerVectorOptionKey const excluded_residues_from_rmsd( "mc:excluded_residues_from_rmsd" );  }
namespace mc { IntegerOptionKey const heat_convergence_check( "mc:heat_convergence_check" );  }
namespace batch_relax { BooleanOptionKey const batch_relax( "batch_relax" );  }
namespace batch_relax { IntegerOptionKey const batch_size( "batch_relax:batch_size" );  }
namespace relax { BooleanOptionKey const relax( "relax" );  }
namespace relax { BooleanOptionKey const fast( "relax:fast" );  }
namespace relax { BooleanOptionKey const thorough( "relax:thorough" );  }
namespace relax { BooleanOptionKey const membrane( "relax:membrane" );  }
namespace relax { BooleanOptionKey const centroid_mode( "relax:centroid_mode" );  }
namespace relax { IntegerOptionKey const default_repeats( "relax:default_repeats" );  }
namespace relax { BooleanOptionKey const dualspace( "relax:dualspace" );  }
namespace relax { BooleanOptionKey const ramady( "relax:ramady" );  }
namespace relax { RealOptionKey const ramady_rms_limit( "relax:ramady_rms_limit" );  }
namespace relax { RealOptionKey const ramady_cutoff( "relax:ramady_cutoff" );  }
namespace relax { IntegerOptionKey const ramady_max_rebuild( "relax:ramady_max_rebuild" );  }
namespace relax { BooleanOptionKey const ramady_force( "relax:ramady_force" );  }
namespace relax { FileOptionKey const script( "relax:script" );  }
namespace relax { IntegerOptionKey const script_max_accept( "relax:script_max_accept" );  }
namespace relax { BooleanOptionKey const superimpose_to_native( "relax:superimpose_to_native" );  }
namespace relax { FileOptionKey const superimpose_to_file( "relax:superimpose_to_file" );  }
namespace relax { BooleanOptionKey const constrain_relax_to_native_coords( "relax:constrain_relax_to_native_coords" );  }
namespace relax { BooleanOptionKey const constrain_relax_to_start_coords( "relax:constrain_relax_to_start_coords" );  }
namespace relax { BooleanOptionKey const coord_constrain_sidechains( "relax:coord_constrain_sidechains" );  }
namespace relax { RealOptionKey const sc_cst_maxdist( "relax:sc_cst_maxdist" );  }
namespace relax { BooleanOptionKey const limit_aroma_chi2( "relax:limit_aroma_chi2" );  }
namespace relax { BooleanOptionKey const respect_resfile( "relax:respect_resfile" );  }
namespace relax { BooleanOptionKey const bb_move( "relax:bb_move" );  }
namespace relax { BooleanOptionKey const chi_move( "relax:chi_move" );  }
namespace relax { BooleanOptionKey const jump_move( "relax:jump_move" );  }
namespace relax { BooleanOptionKey const dna_move( "relax:dna_move" );  }
namespace relax { BooleanOptionKey const fix_omega( "relax:fix_omega" );  }
namespace relax { BooleanOptionKey const minimize_bond_lengths( "relax:minimize_bond_lengths" );  }
namespace relax { BooleanOptionKey const minimize_bond_angles( "relax:minimize_bond_angles" );  }
namespace relax { IntegerOptionKey const minimize_bondlength_subset( "relax:minimize_bondlength_subset" );  }
namespace relax { IntegerOptionKey const minimize_bondangle_subset( "relax:minimize_bondangle_subset" );  }
namespace relax { StringOptionKey const min_type( "relax:min_type" );  }
namespace relax { BooleanOptionKey const cartesian( "relax:cartesian" );  }
namespace relax { RealOptionKey const chainbreak_weight( "relax:chainbreak_weight" );  }
namespace relax { RealOptionKey const linear_chainbreak_weight( "relax:linear_chainbreak_weight" );  }
namespace relax { RealOptionKey const overlap_chainbreak_weight( "relax:overlap_chainbreak_weight" );  }
namespace relax { BooleanOptionKey const classic( "relax:classic" );  }
namespace relax { FileOptionKey const sequence_file( "relax:sequence_file" );  }
namespace relax { BooleanOptionKey const quick( "relax:quick" );  }
namespace relax { BooleanOptionKey const sequence( "relax:sequence" );  }
namespace relax { IntegerOptionKey const minirelax_repeats( "relax:minirelax_repeats" );  }
namespace relax { RealOptionKey const minirelax_sdev( "relax:minirelax_sdev" );  }
namespace relax { BooleanOptionKey const wobblemoves( "relax:wobblemoves" );  }
namespace relax { FileOptionKey const constrain_relax_segments( "relax:constrain_relax_segments" );  }
namespace relax { RealOptionKey const coord_cst_width( "relax:coord_cst_width" );  }
namespace relax { RealOptionKey const coord_cst_stdev( "relax:coord_cst_stdev" );  }
namespace relax { BooleanOptionKey const ramp_constraints( "relax:ramp_constraints" );  }
namespace relax { RealOptionKey const energycut( "relax:energycut" );  }
namespace relax { BooleanOptionKey const mini( "relax:mini" );  }
namespace relax { IntegerOptionKey const stage1_ramp_cycles( "relax:stage1_ramp_cycles" );  }
namespace relax { IntegerOptionKey const stage1_ramp_inner_cycles( "relax:stage1_ramp_inner_cycles" );  }
namespace relax { IntegerOptionKey const stage2_repack_period( "relax:stage2_repack_period" );  }
namespace relax { IntegerOptionKey const stage2_cycles( "relax:stage2_cycles" );  }
namespace relax { RealOptionKey const min_tolerance( "relax:min_tolerance" );  }
namespace relax { IntegerOptionKey const stage3_cycles( "relax:stage3_cycles" );  }
namespace relax { RealOptionKey const cycle_ratio( "relax:cycle_ratio" );  }
namespace relax { RealOptionKey const filter_stage2_beginning( "relax:filter_stage2_beginning" );  }
namespace relax { RealOptionKey const filter_stage2_quarter( "relax:filter_stage2_quarter" );  }
namespace relax { RealOptionKey const filter_stage2_half( "relax:filter_stage2_half" );  }
namespace relax { RealOptionKey const filter_stage2_end( "relax:filter_stage2_end" );  }
namespace relax { namespace centroid { BooleanOptionKey const centroid( "relax:centroid" );  } }
namespace relax { namespace centroid { StringOptionKey const weights( "relax:centroid:weights" );  } }
namespace relax { namespace centroid { BooleanOptionKey const ramp_vdw( "relax:centroid:ramp_vdw" );  } }
namespace relax { namespace centroid { BooleanOptionKey const ramp_rama( "relax:centroid:ramp_rama" );  } }
namespace relax { namespace centroid { StringOptionKey const parameters( "relax:centroid:parameters" );  } }
namespace relax { namespace centroid { BooleanOptionKey const do_final_repack( "relax:centroid:do_final_repack" );  } }
namespace relax { namespace centroid { BooleanOptionKey const increase_vdw_radii( "relax:centroid:increase_vdw_radii" );  } }
namespace enzdes { BooleanOptionKey const enzdes( "enzdes" );  }
namespace enzdes { StringOptionKey const checkpoint( "enzdes:checkpoint" );  }
namespace enzdes { BooleanOptionKey const enz_score( "enzdes:enz_score" );  }
namespace enzdes { BooleanOptionKey const enz_repack( "enzdes:enz_repack" );  }
namespace enzdes { BooleanOptionKey const cst_opt( "enzdes:cst_opt" );  }
namespace enzdes { BooleanOptionKey const cst_predock( "enzdes:cst_predock" );  }
namespace enzdes { RealOptionKey const trans_magnitude( "enzdes:trans_magnitude" );  }
namespace enzdes { RealOptionKey const rot_magnitude( "enzdes:rot_magnitude" );  }
namespace enzdes { RealOptionKey const dock_trials( "enzdes:dock_trials" );  }
namespace enzdes { BooleanOptionKey const cst_min( "enzdes:cst_min" );  }
namespace enzdes { BooleanOptionKey const cst_design( "enzdes:cst_design" );  }
namespace enzdes { IntegerOptionKey const design_min_cycles( "enzdes:design_min_cycles" );  }
namespace enzdes { BooleanOptionKey const make_consensus_mutations( "enzdes:make_consensus_mutations" );  }
namespace enzdes { BooleanOptionKey const bb_min( "enzdes:bb_min" );  }
namespace enzdes { RealOptionKey const bb_min_allowed_dev( "enzdes:bb_min_allowed_dev" );  }
namespace enzdes { RealOptionKey const loop_bb_min_allowed_dev( "enzdes:loop_bb_min_allowed_dev" );  }
namespace enzdes { RealOptionKey const minimize_ligand_torsions( "enzdes:minimize_ligand_torsions" );  }
namespace enzdes { RealOptionKey const minimize_all_ligand_torsions( "enzdes:minimize_all_ligand_torsions" );  }
namespace enzdes { BooleanOptionKey const chi_min( "enzdes:chi_min" );  }
namespace enzdes { BooleanOptionKey const min_all_jumps( "enzdes:min_all_jumps" );  }
namespace enzdes { BooleanOptionKey const cst_dock( "enzdes:cst_dock" );  }
namespace enzdes { BooleanOptionKey const run_ligand_motifs( "enzdes:run_ligand_motifs" );  }
namespace enzdes { BooleanOptionKey const enz_debug( "enzdes:enz_debug" );  }
namespace enzdes { FileOptionKey const cstfile( "enzdes:cstfile" );  }
namespace enzdes { FileOptionKey const enz_loops_file( "enzdes:enz_loops_file" );  }
namespace enzdes { BooleanOptionKey const flexbb_protocol( "enzdes:flexbb_protocol" );  }
namespace enzdes { BooleanOptionKey const remodel_protocol( "enzdes:remodel_protocol" );  }
namespace enzdes { BooleanOptionKey const kic_loop_sampling( "enzdes:kic_loop_sampling" );  }
namespace enzdes { StringOptionKey const dump_loop_samples( "enzdes:dump_loop_samples" );  }
namespace enzdes { BooleanOptionKey const fix_catalytic_aa( "enzdes:fix_catalytic_aa" );  }
namespace enzdes { IntegerOptionKey const additional_packing_ligand_rb_confs( "enzdes:additional_packing_ligand_rb_confs" );  }
namespace enzdes { IntegerOptionKey const ex_catalytic_rot( "enzdes:ex_catalytic_rot" );  }
namespace enzdes { IntegerOptionKey const single_loop_ensemble_size( "enzdes:single_loop_ensemble_size" );  }
namespace enzdes { IntegerOptionKey const loop_generator_trials( "enzdes:loop_generator_trials" );  }
namespace enzdes { BooleanOptionKey const no_catres_min_in_loopgen( "enzdes:no_catres_min_in_loopgen" );  }
namespace enzdes { RealOptionKey const mc_kt_low( "enzdes:mc_kt_low" );  }
namespace enzdes { RealOptionKey const mc_kt_high( "enzdes:mc_kt_high" );  }
namespace enzdes { RealOptionKey const min_cacb_deviation( "enzdes:min_cacb_deviation" );  }
namespace enzdes { RealOptionKey const max_bb_deviation( "enzdes:max_bb_deviation" );  }
namespace enzdes { RealOptionKey const max_bb_deviation_from_startstruct( "enzdes:max_bb_deviation_from_startstruct" );  }
namespace enzdes { IntegerOptionKey const flexbb_outstructs( "enzdes:flexbb_outstructs" );  }
namespace enzdes { IntegerOptionKey const remodel_trials( "enzdes:remodel_trials" );  }
namespace enzdes { BooleanOptionKey const remodel_secmatch( "enzdes:remodel_secmatch" );  }
namespace enzdes { BooleanOptionKey const dump_inverse_rotamers( "enzdes:dump_inverse_rotamers" );  }
namespace enzdes { RealOptionKey const remodel_aggressiveness( "enzdes:remodel_aggressiveness" );  }
namespace enzdes { RealOptionKey const favor_native_res( "enzdes:favor_native_res" );  }
namespace enzdes { BooleanOptionKey const detect_design_interface( "enzdes:detect_design_interface" );  }
namespace enzdes { BooleanOptionKey const include_catres_in_interface_detection( "enzdes:include_catres_in_interface_detection" );  }
namespace enzdes { BooleanOptionKey const arg_sweep_interface( "enzdes:arg_sweep_interface" );  }
namespace enzdes { RealOptionKey const arg_sweep_cutoff( "enzdes:arg_sweep_cutoff" );  }
namespace enzdes { RealOptionKey const cut1( "enzdes:cut1" );  }
namespace enzdes { RealOptionKey const cut2( "enzdes:cut2" );  }
namespace enzdes { RealOptionKey const cut3( "enzdes:cut3" );  }
namespace enzdes { RealOptionKey const cut4( "enzdes:cut4" );  }
namespace enzdes { RealOptionKey const lig_packer_weight( "enzdes:lig_packer_weight" );  }
namespace enzdes { BooleanOptionKey const no_unconstrained_repack( "enzdes:no_unconstrained_repack" );  }
namespace enzdes { RealOptionKey const secmatch_Ecutoff( "enzdes:secmatch_Ecutoff" );  }
namespace enzdes { FileOptionKey const change_lig( "enzdes:change_lig" );  }
namespace enzdes { StringOptionKey const process_ligrot_separately( "enzdes:process_ligrot_separately" );  }
namespace enzdes { BooleanOptionKey const start_from_random_rb_conf( "enzdes:start_from_random_rb_conf" );  }
namespace enzdes { RealOptionKey const bb_bump_cutoff( "enzdes:bb_bump_cutoff" );  }
namespace enzdes { RealOptionKey const sc_sc_bump_cutoff( "enzdes:sc_sc_bump_cutoff" );  }
namespace enzdes { BooleanOptionKey const no_packstat_calculation( "enzdes:no_packstat_calculation" );  }
namespace enzdes { StringOptionKey const compare_native( "enzdes:compare_native" );  }
namespace enzdes { BooleanOptionKey const final_repack_without_ligand( "enzdes:final_repack_without_ligand" );  }
namespace enzdes { BooleanOptionKey const dump_final_repack_without_ligand_pdb( "enzdes:dump_final_repack_without_ligand_pdb" );  }
namespace enzdes { BooleanOptionKey const parser_read_cloud_pdb( "enzdes:parser_read_cloud_pdb" );  }
namespace packing { BooleanOptionKey const packing( "packing" );  }
namespace packing { BooleanOptionKey const repack_only( "packing:repack_only" );  }
namespace packing { BooleanOptionKey const prevent_repacking( "packing:prevent_repacking" );  }
namespace packing { RealOptionKey const cenrot_cutoff( "packing:cenrot_cutoff" );  }
namespace packing { IntegerOptionKey const ndruns( "packing:ndruns" );  }
namespace packing { BooleanOptionKey const soft_rep_design( "packing:soft_rep_design" );  }
namespace packing { BooleanOptionKey const use_electrostatic_repulsion( "packing:use_electrostatic_repulsion" );  }
namespace packing { BooleanOptionKey const dump_rotamer_sets( "packing:dump_rotamer_sets" );  }
namespace packing { RealOptionKey const dunbrack_prob_buried( "packing:dunbrack_prob_buried" );  }
namespace packing { RealOptionKey const dunbrack_prob_nonburied( "packing:dunbrack_prob_nonburied" );  }
namespace packing { RealOptionKey const dunbrack_prob_nonburied_semirotameric( "packing:dunbrack_prob_nonburied_semirotameric" );  }
namespace packing { BooleanOptionKey const no_optH( "packing:no_optH" );  }
namespace packing { BooleanOptionKey const optH_MCA( "packing:optH_MCA" );  }
namespace packing { BooleanOptionKey const pack_missing_sidechains( "packing:pack_missing_sidechains" );  }
namespace packing { BooleanOptionKey const preserve_c_beta( "packing:preserve_c_beta" );  }
namespace packing { BooleanOptionKey const flip_HNQ( "packing:flip_HNQ" );  }
namespace packing { IntegerVectorOptionKey const fix_his_tautomer( "packing:fix_his_tautomer" );  }
namespace packing { BooleanOptionKey const print_pymol_selection( "packing:print_pymol_selection" );  }
namespace packing { namespace ex1 { BooleanOptionKey const ex1( "packing:ex1" );  } }
namespace packing { namespace ex1 { IntegerOptionKey const level( "packing:ex1:level" );  } }
namespace packing { namespace ex1 { BooleanOptionKey const operate( "packing:ex1:operate" );  } }
namespace packing { namespace ex2 { BooleanOptionKey const ex2( "packing:ex2" );  } }
namespace packing { namespace ex2 { IntegerOptionKey const level( "packing:ex2:level" );  } }
namespace packing { namespace ex2 { BooleanOptionKey const operate( "packing:ex2:operate" );  } }
namespace packing { namespace ex3 { BooleanOptionKey const ex3( "packing:ex3" );  } }
namespace packing { namespace ex3 { IntegerOptionKey const level( "packing:ex3:level" );  } }
namespace packing { namespace ex3 { BooleanOptionKey const operate( "packing:ex3:operate" );  } }
namespace packing { namespace ex4 { BooleanOptionKey const ex4( "packing:ex4" );  } }
namespace packing { namespace ex4 { IntegerOptionKey const level( "packing:ex4:level" );  } }
namespace packing { namespace ex4 { BooleanOptionKey const operate( "packing:ex4:operate" );  } }
namespace packing { namespace ex1aro { BooleanOptionKey const ex1aro( "packing:ex1aro" );  } }
namespace packing { namespace ex1aro { IntegerOptionKey const level( "packing:ex1aro:level" );  } }
namespace packing { namespace ex2aro { BooleanOptionKey const ex2aro( "packing:ex2aro" );  } }
namespace packing { namespace ex2aro { IntegerOptionKey const level( "packing:ex2aro:level" );  } }
namespace packing { namespace ex1aro_exposed { BooleanOptionKey const ex1aro_exposed( "packing:ex1aro_exposed" );  } }
namespace packing { namespace ex1aro_exposed { IntegerOptionKey const level( "packing:ex1aro_exposed:level" );  } }
namespace packing { namespace ex2aro_exposed { BooleanOptionKey const ex2aro_exposed( "packing:ex2aro_exposed" );  } }
namespace packing { namespace ex2aro_exposed { IntegerOptionKey const level( "packing:ex2aro_exposed:level" );  } }
namespace packing { namespace exdna { BooleanOptionKey const exdna( "packing:exdna" );  } }
namespace packing { namespace exdna { IntegerOptionKey const level( "packing:exdna:level" );  } }
namespace packing { IntegerOptionKey const extrachi_cutoff( "packing:extrachi_cutoff" );  }
namespace packing { FileVectorOptionKey const resfile( "packing:resfile" );  }
namespace packing { RealOptionKey const outeriterations_scaling( "packing:outeriterations_scaling" );  }
namespace packing { RealOptionKey const inneriterations_scaling( "packing:inneriterations_scaling" );  }
namespace packing { BooleanOptionKey const explicit_h2o( "packing:explicit_h2o" );  }
namespace packing { StringVectorOptionKey const adducts( "packing:adducts" );  }
namespace packing { BooleanOptionKey const solvate( "packing:solvate" );  }
namespace packing { BooleanOptionKey const use_input_sc( "packing:use_input_sc" );  }
namespace packing { FileVectorOptionKey const unboundrot( "packing:unboundrot" );  }
namespace packing { RealOptionKey const max_rotbump_energy( "packing:max_rotbump_energy" );  }
namespace packing { BooleanOptionKey const lazy_ig( "packing:lazy_ig" );  }
namespace packing { BooleanOptionKey const double_lazy_ig( "packing:double_lazy_ig" );  }
namespace packing { IntegerOptionKey const double_lazy_ig_mem_limit( "packing:double_lazy_ig_mem_limit" );  }
namespace packing { IntegerOptionKey const linmem_ig( "packing:linmem_ig" );  }
namespace packing { IntegerOptionKey const multi_cool_annealer( "packing:multi_cool_annealer" );  }
namespace packing { RealVectorOptionKey const minpack_temp_schedule( "packing:minpack_temp_schedule" );  }
namespace packing { IntegerOptionKey const minpack_inner_iteration_scale( "packing:minpack_inner_iteration_scale" );  }
namespace packing { BooleanOptionKey const minpack_disable_bumpcheck( "packing:minpack_disable_bumpcheck" );  }
namespace phil { BooleanOptionKey const phil( "phil" );  }
namespace phil { IntegerOptionKey const nloop( "phil:nloop" );  }
namespace phil { StringOptionKey const vall_file( "phil:vall_file" );  }
namespace phil { StringOptionKey const align_file( "phil:align_file" );  }
namespace wum { BooleanOptionKey const wum( "wum" );  }
namespace wum { IntegerOptionKey const n_slaves_per_master( "wum:n_slaves_per_master" );  }
namespace wum { IntegerOptionKey const n_masters( "wum:n_masters" );  }
namespace wum { IntegerOptionKey const memory_limit( "wum:memory_limit" );  }
namespace wum { StringOptionKey const extra_scorefxn( "wum:extra_scorefxn" );  }
namespace wum { FileOptionKey const extra_scorefxn_ref_structure( "wum:extra_scorefxn_ref_structure" );  }
namespace wum { IntegerOptionKey const extra_scorefxn_relax( "wum:extra_scorefxn_relax" );  }
namespace wum { RealOptionKey const trim_proportion( "wum:trim_proportion" );  }
namespace els { BooleanOptionKey const els( "els" );  }
namespace els { IntegerOptionKey const master_wu_per_send( "els:master_wu_per_send" );  }
namespace els { StringOptionKey const vars( "els:vars" );  }
namespace els { FileOptionKey const script( "els:script" );  }
namespace els { IntegerOptionKey const num_traj( "els:num_traj" );  }
namespace els { IntegerOptionKey const traj_per_master( "els:traj_per_master" );  }
namespace els { IntegerOptionKey const shortest_wu( "els:shortest_wu" );  }
namespace els { BooleanOptionKey const pool( "els:pool" );  }
namespace els { BooleanOptionKey const singlenode( "els:singlenode" );  }
namespace lh { BooleanOptionKey const lh( "lh" );  }
namespace lh { StringOptionKey const db_prefix( "lh:db_prefix" );  }
namespace lh { IntegerVectorOptionKey const loopsizes( "lh:loopsizes" );  }
namespace lh { IntegerOptionKey const num_partitions( "lh:num_partitions" );  }
namespace lh { PathOptionKey const db_path( "lh:db_path" );  }
namespace lh { BooleanOptionKey const exclude_homo( "lh:exclude_homo" );  }
namespace lh { BooleanOptionKey const bss( "lh:bss" );  }
namespace lh { StringOptionKey const refstruct( "lh:refstruct" );  }
namespace lh { StringOptionKey const homo_file( "lh:homo_file" );  }
namespace lh { RealVectorOptionKey const createdb_rms_cutoff( "lh:createdb_rms_cutoff" );  }
namespace lh { RealOptionKey const min_bbrms( "lh:min_bbrms" );  }
namespace lh { RealOptionKey const max_bbrms( "lh:max_bbrms" );  }
namespace lh { RealOptionKey const min_rms( "lh:min_rms" );  }
namespace lh { RealOptionKey const max_rms( "lh:max_rms" );  }
namespace lh { BooleanOptionKey const filter_by_phipsi( "lh:filter_by_phipsi" );  }
namespace lh { IntegerOptionKey const max_radius( "lh:max_radius" );  }
namespace lh { IntegerOptionKey const max_struct( "lh:max_struct" );  }
namespace lh { IntegerOptionKey const max_struct_per_radius( "lh:max_struct_per_radius" );  }
namespace lh { RealOptionKey const grid_space_multiplier( "lh:grid_space_multiplier" );  }
namespace lh { RealOptionKey const grid_angle_multiplier( "lh:grid_angle_multiplier" );  }
namespace lh { IntegerOptionKey const skim_size( "lh:skim_size" );  }
namespace lh { IntegerOptionKey const rounds( "lh:rounds" );  }
namespace lh { StringOptionKey const jobname( "lh:jobname" );  }
namespace lh { IntegerOptionKey const max_lib_size( "lh:max_lib_size" );  }
namespace lh { IntegerOptionKey const max_emperor_lib_size( "lh:max_emperor_lib_size" );  }
namespace lh { IntegerOptionKey const max_emperor_lib_round( "lh:max_emperor_lib_round" );  }
namespace lh { IntegerOptionKey const library_expiry_time( "lh:library_expiry_time" );  }
namespace lh { StringOptionKey const objective_function( "lh:objective_function" );  }
namespace lh { IntegerOptionKey const expire_after_rounds( "lh:expire_after_rounds" );  }
namespace lh { StringOptionKey const mpi_resume( "lh:mpi_resume" );  }
namespace lh { StringOptionKey const mpi_feedback( "lh:mpi_feedback" );  }
namespace lh { IntegerOptionKey const mpi_batch_relax_chunks( "lh:mpi_batch_relax_chunks" );  }
namespace lh { IntegerOptionKey const mpi_batch_relax_absolute_max( "lh:mpi_batch_relax_absolute_max" );  }
namespace lh { IntegerOptionKey const mpi_outbound_wu_buffer_size( "lh:mpi_outbound_wu_buffer_size" );  }
namespace lh { IntegerOptionKey const mpi_loophash_split_size    ( "lh:mpi_loophash_split_size    " );  }
namespace lh { RealOptionKey const mpi_metropolis_temp( "lh:mpi_metropolis_temp" );  }
namespace lh { IntegerOptionKey const mpi_save_state_interval( "lh:mpi_save_state_interval" );  }
namespace lh { BooleanOptionKey const mpi_master_save_score_only( "lh:mpi_master_save_score_only" );  }
namespace lh { IntegerOptionKey const max_loophash_per_structure( "lh:max_loophash_per_structure" );  }
namespace lh { RealOptionKey const rms_limit( "lh:rms_limit" );  }
namespace lh { BooleanOptionKey const centroid_only( "lh:centroid_only" );  }
namespace lh { BooleanOptionKey const write_centroid_structs( "lh:write_centroid_structs" );  }
namespace lh { BooleanOptionKey const write_all_fa_structs( "lh:write_all_fa_structs" );  }
namespace lh { BooleanOptionKey const sandbox( "lh:sandbox" );  }
namespace lh { BooleanOptionKey const create_db( "lh:create_db" );  }
namespace lh { FileOptionKey const sample_weight_file( "lh:sample_weight_file" );  }
namespace lh { namespace fragpdb { BooleanOptionKey const fragpdb( "lh:fragpdb" );  } }
namespace lh { namespace fragpdb { StringOptionKey const out_path( "lh:fragpdb:out_path" );  } }
namespace lh { namespace fragpdb { IntegerVectorOptionKey const indexoffset( "lh:fragpdb:indexoffset" );  } }
namespace lh { namespace fragpdb { StringVectorOptionKey const bin( "lh:fragpdb:bin" );  } }
namespace lh { namespace symfragrm { BooleanOptionKey const symfragrm( "lh:symfragrm" );  } }
namespace lh { namespace symfragrm { FileVectorOptionKey const pdblist( "lh:symfragrm:pdblist" );  } }
namespace rbe { BooleanOptionKey const rbe( "rbe" );  }
namespace rbe { StringOptionKey const server_url( "rbe:server_url" );  }
namespace rbe { StringOptionKey const server_port( "rbe:server_port" );  }
namespace rbe { RealOptionKey const poll_frequency( "rbe:poll_frequency" );  }
namespace blivens { BooleanOptionKey const blivens( "blivens" );  }
namespace blivens { namespace disulfide_scorer { BooleanOptionKey const disulfide_scorer( "blivens:disulfide_scorer" );  } }
namespace blivens { namespace disulfide_scorer { RealOptionKey const nds_prob( "blivens:disulfide_scorer:nds_prob" );  } }
namespace blivens { namespace disulfide_scorer { RealOptionKey const cys_prob( "blivens:disulfide_scorer:cys_prob" );  } }
namespace blivens { StringOptionKey const score_type( "blivens:score_type" );  }
namespace krassk { BooleanOptionKey const krassk( "krassk" );  }
namespace krassk { IntegerOptionKey const left_tail( "krassk:left_tail" );  }
namespace krassk { IntegerOptionKey const right_tail( "krassk:right_tail" );  }
namespace krassk { BooleanOptionKey const tail_mode( "krassk:tail_mode" );  }
namespace krassk { IntegerOptionKey const tail_mode_name( "krassk:tail_mode_name" );  }
namespace krassk { StringOptionKey const tail_output_file_name( "krassk:tail_output_file_name" );  }
namespace rotamerdump { BooleanOptionKey const rotamerdump( "rotamerdump" );  }
namespace rotamerdump { BooleanOptionKey const xyz( "rotamerdump:xyz" );  }
namespace rotamerdump { BooleanOptionKey const one_body( "rotamerdump:one_body" );  }
namespace rotamerdump { BooleanOptionKey const two_body( "rotamerdump:two_body" );  }
namespace rotamerdump { BooleanOptionKey const annealer( "rotamerdump:annealer" );  }
namespace robert { BooleanOptionKey const robert( "robert" );  }
namespace robert { StringOptionKey const pairdata_input_pdb_list( "robert:pairdata_input_pdb_list" );  }
namespace robert { RealOptionKey const pcs_maxsub_filter( "robert:pcs_maxsub_filter" );  }
namespace robert { RealOptionKey const pcs_maxsub_rmsd( "robert:pcs_maxsub_rmsd" );  }
namespace robert { BooleanOptionKey const pcs_dump_cluster( "robert:pcs_dump_cluster" );  }
namespace robert { RealOptionKey const pcs_cluster_coverage( "robert:pcs_cluster_coverage" );  }
namespace robert { BooleanOptionKey const pcs_cluster_lowscoring( "robert:pcs_cluster_lowscoring" );  }
namespace cmiles { BooleanOptionKey const cmiles( "cmiles" );  }
namespace cmiles { namespace kcluster { BooleanOptionKey const kcluster( "cmiles:kcluster" );  } }
namespace cmiles { namespace kcluster { IntegerOptionKey const num_clusters( "cmiles:kcluster:num_clusters" );  } }
namespace cmiles { namespace jumping { BooleanOptionKey const jumping( "cmiles:jumping" );  } }
namespace cmiles { namespace jumping { IntegerOptionKey const resi( "cmiles:jumping:resi" );  } }
