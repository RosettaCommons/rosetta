namespace casp { StringOptionKey const weight_file( "casp:weight_file" );  }
namespace casp { StringOptionKey const refine_res( "casp:refine_res" );  }
namespace pose_metrics { BooleanOptionKey const pose_metrics( "pose_metrics" );  }
namespace pose_metrics { RealOptionKey const atomic_burial_cutoff( "pose_metrics:atomic_burial_cutoff" );  }
namespace pose_metrics { RealOptionKey const sasa_calculator_probe_radius( "pose_metrics:sasa_calculator_probe_radius" );  }
namespace pose_metrics { RealOptionKey const interface_cutoff( "pose_metrics:interface_cutoff" );  }
namespace pose_metrics { IntegerOptionKey const min_sequence_separation( "pose_metrics:min_sequence_separation" );  }
namespace pose_metrics { RealOptionKey const contact_cutoffE( "pose_metrics:contact_cutoffE" );  }
namespace pose_metrics { RealOptionKey const neighbor_by_distance_cutoff( "pose_metrics:neighbor_by_distance_cutoff" );  }
namespace pose_metrics { RealOptionKey const inter_group_neighbors_cutoff( "pose_metrics:inter_group_neighbors_cutoff" );  }
namespace pose_metrics { RealOptionKey const semiex_water_burial_cutoff( "pose_metrics:semiex_water_burial_cutoff" );  }
namespace ddg { BooleanOptionKey const ddg( "ddg" );  }
namespace ddg { BooleanOptionKey const avg_rot_cst_enrg( "ddg:avg_rot_cst_enrg" );  }
namespace ddg { BooleanOptionKey const use_bound_cst( "ddg:use_bound_cst" );  }
namespace ddg { RealOptionKey const cap_rot_cst_enrg( "ddg:cap_rot_cst_enrg" );  }
namespace ddg { BooleanOptionKey const opt_input_structure( "ddg:opt_input_structure" );  }
namespace ddg { BooleanOptionKey const pack_until_converge( "ddg:pack_until_converge" );  }
namespace ddg { BooleanOptionKey const no_constraints( "ddg:no_constraints" );  }
namespace ddg { RealOptionKey const apply_rot_cst_to_mutation_region_only( "ddg:apply_rot_cst_to_mutation_region_only" );  }
namespace ddg { RealOptionKey const rot_cst_enrg_cutoff( "ddg:rot_cst_enrg_cutoff" );  }
namespace ddg { BooleanOptionKey const use_rotamer_constraints_to_native( "ddg:use_rotamer_constraints_to_native" );  }
namespace ddg { BooleanOptionKey const single_res_scoring( "ddg:single_res_scoring" );  }
namespace ddg { BooleanOptionKey const downweight_by_sasa( "ddg:downweight_by_sasa" );  }
namespace ddg { BooleanOptionKey const global( "ddg:global" );  }
namespace ddg { BooleanOptionKey const exclude_solvent_exposed_res( "ddg:exclude_solvent_exposed_res" );  }
namespace ddg { RealOptionKey const radius( "ddg:radius" );  }
namespace ddg { StringOptionKey const wt( "ddg:wt" );  }
namespace ddg { StringOptionKey const mut( "ddg:mut" );  }
namespace ddg { BooleanOptionKey const suppress_checkpointing( "ddg:suppress_checkpointing" );  }
namespace ddg { BooleanOptionKey const wt_only( "ddg:wt_only" );  }
namespace ddg { BooleanOptionKey const mut_only( "ddg:mut_only" );  }
namespace ddg { BooleanOptionKey const output_silent( "ddg:output_silent" );  }
namespace ddg { StringOptionKey const minimization_scorefunction( "ddg:minimization_scorefunction" );  }
namespace ddg { StringOptionKey const minimization_patch( "ddg:minimization_patch" );  }
namespace ddg { BooleanOptionKey const min_cst( "ddg:min_cst" );  }
namespace ddg { IntegerOptionKey const lowest_x_decoys( "ddg:lowest_x_decoys" );  }
namespace ddg { BooleanOptionKey const local_opt_only( "ddg:local_opt_only" );  }
namespace ddg { BooleanOptionKey const print_per_res_diff( "ddg:print_per_res_diff" );  }
namespace ddg { BooleanOptionKey const mean( "ddg:mean" );  }
namespace ddg { BooleanOptionKey const min( "ddg:min" );  }
namespace ddg { BooleanOptionKey const rb_restrict_to_mutsite_nbrs( "ddg:rb_restrict_to_mutsite_nbrs" );  }
namespace ddg { BooleanOptionKey const no_bb_movement( "ddg:no_bb_movement" );  }
namespace ddg { BooleanOptionKey const initial_repack( "ddg:initial_repack" );  }
namespace ddg { StringOptionKey const rb_file( "ddg:rb_file" );  }
namespace ddg { IntegerOptionKey const interface_ddg( "ddg:interface_ddg" );  }
namespace ddg { RealOptionKey const ens_variation( "ddg:ens_variation" );  }
namespace ddg { BooleanOptionKey const sc_min_only( "ddg:sc_min_only" );  }
namespace ddg { StringOptionKey const min_cst_weights( "ddg:min_cst_weights" );  }
namespace ddg { RealOptionKey const opt_radius( "ddg:opt_radius" );  }
namespace ddg { StringOptionKey const output_dir( "ddg:output_dir" );  }
namespace ddg { StringOptionKey const last_accepted_pose_dir( "ddg:last_accepted_pose_dir" );  }
namespace ddg { BooleanOptionKey const min_with_cst( "ddg:min_with_cst" );  }
namespace ddg { RealOptionKey const temperature( "ddg:temperature" );  }
namespace ddg { BooleanOptionKey const ramp_repulsive( "ddg:ramp_repulsive" );  }
namespace ddg { StringOptionKey const mut_file( "ddg:mut_file" );  }
namespace ddg { StringOptionKey const out_pdb_prefix( "ddg:out_pdb_prefix" );  }
namespace ddg { RealOptionKey const constraint_weight( "ddg:constraint_weight" );  }
namespace ddg { RealOptionKey const harmonic_ca_tether( "ddg:harmonic_ca_tether" );  }
namespace ddg { IntegerOptionKey const iterations( "ddg:iterations" );  }
namespace ddg { StringOptionKey const out( "ddg:out" );  }
namespace ddg { BooleanOptionKey const debug_output( "ddg:debug_output" );  }
namespace ddg { BooleanOptionKey const dump_pdbs( "ddg:dump_pdbs" );  }
namespace ddg { StringOptionKey const weight_file( "ddg:weight_file" );  }
namespace ddg { IntegerOptionKey const translate_by( "ddg:translate_by" );  }
namespace murphp { BooleanOptionKey const murphp( "murphp" );  }
namespace murphp { StringOptionKey const inv_kin_lig_loop_design_filename( "murphp:inv_kin_lig_loop_design_filename" );  }
namespace motifs { BooleanOptionKey const motifs( "motifs" );  }
namespace motifs { RealOptionKey const close_enough( "motifs:close_enough" );  }
namespace motifs { IntegerOptionKey const max_depth( "motifs:max_depth" );  }
namespace motifs { BooleanOptionKey const keep_motif_xtal_location( "motifs:keep_motif_xtal_location" );  }
namespace motifs { RealOptionKey const pack_score_cutoff( "motifs:pack_score_cutoff" );  }
namespace motifs { RealOptionKey const hb_score_cutoff( "motifs:hb_score_cutoff" );  }
namespace motifs { RealOptionKey const water_score_cutoff( "motifs:water_score_cutoff" );  }
namespace motifs { StringOptionKey const motif_output_directory( "motifs:motif_output_directory" );  }
namespace motifs { BooleanOptionKey const eliminate_weak_motifs( "motifs:eliminate_weak_motifs" );  }
namespace motifs { RealOptionKey const duplicate_motif_cutoff( "motifs:duplicate_motif_cutoff" );  }
namespace motifs { BooleanOptionKey const preminimize_motif_pdbs( "motifs:preminimize_motif_pdbs" );  }
namespace motifs { BooleanOptionKey const preminimize_motif_pdbs_sconly( "motifs:preminimize_motif_pdbs_sconly" );  }
namespace motifs { BooleanOptionKey const place_adduct_waters( "motifs:place_adduct_waters" );  }
namespace motifs { FileVectorOptionKey const list_motifs( "motifs:list_motifs" );  }
namespace motifs { StringOptionKey const motif_filename( "motifs:motif_filename" );  }
namespace motifs { StringOptionKey const BPData( "motifs:BPData" );  }
namespace motifs { FileVectorOptionKey const list_dnaconformers( "motifs:list_dnaconformers" );  }
namespace motifs { StringVectorOptionKey const target_dna_defs( "motifs:target_dna_defs" );  }
namespace motifs { StringVectorOptionKey const motif_build_defs( "motifs:motif_build_defs" );  }
namespace motifs { StringOptionKey const motif_build_position_chain( "motifs:motif_build_position_chain" );  }
namespace motifs { IntegerVectorOptionKey const motif_build_positions( "motifs:motif_build_positions" );  }
namespace motifs { RealOptionKey const r1( "motifs:r1" );  }
namespace motifs { RealOptionKey const r2( "motifs:r2" );  }
namespace motifs { RealOptionKey const z1( "motifs:z1" );  }
namespace motifs { RealOptionKey const z2( "motifs:z2" );  }
namespace motifs { RealOptionKey const dtest( "motifs:dtest" );  }
namespace motifs { IntegerOptionKey const rotlevel( "motifs:rotlevel" );  }
namespace motifs { IntegerOptionKey const num_repacks( "motifs:num_repacks" );  }
namespace motifs { BooleanOptionKey const minimize( "motifs:minimize" );  }
namespace motifs { BooleanOptionKey const minimize_dna( "motifs:minimize_dna" );  }
namespace motifs { BooleanOptionKey const run_motifs( "motifs:run_motifs" );  }
namespace motifs { BooleanOptionKey const expand_motifs( "motifs:expand_motifs" );  }
namespace motifs { BooleanOptionKey const aromatic_motifs( "motifs:aromatic_motifs" );  }
namespace motifs { BooleanOptionKey const dump_motifs( "motifs:dump_motifs" );  }
namespace motifs { BooleanOptionKey const quick_and_dirty( "motifs:quick_and_dirty" );  }
namespace motifs { RealOptionKey const special_rotweight( "motifs:special_rotweight" );  }
namespace motifs { StringOptionKey const output_file( "motifs:output_file" );  }
namespace motifs { StringOptionKey const data_file( "motifs:data_file" );  }
namespace motifs { RealOptionKey const constraint_max( "motifs:constraint_max" );  }
namespace motifs { BooleanOptionKey const flex_sugar( "motifs:flex_sugar" );  }
namespace motifs { BooleanOptionKey const clear_bprots( "motifs:clear_bprots" );  }
namespace motifs { IntegerOptionKey const rots2add( "motifs:rots2add" );  }
namespace motifs { BooleanOptionKey const restrict_to_wt( "motifs:restrict_to_wt" );  }
namespace motifs { BooleanOptionKey const rerun_motifsearch( "motifs:rerun_motifsearch" );  }
namespace motifs { BooleanOptionKey const no_rotamer_bump( "motifs:no_rotamer_bump" );  }
namespace motifs { RealOptionKey const ligand_motif_sphere( "motifs:ligand_motif_sphere" );  }
namespace constraints { BooleanOptionKey const constraints( "constraints" );  }
namespace constraints { RealOptionKey const CA_tether( "constraints:CA_tether" );  }
namespace constraints { BooleanOptionKey const exit_on_bad_read( "constraints:exit_on_bad_read" );  }
namespace constraints { StringVectorOptionKey const cst_file( "constraints:cst_file" );  }
namespace constraints { RealOptionKey const cst_weight( "constraints:cst_weight" );  }
namespace constraints { RealOptionKey const max_cst_dist( "constraints:max_cst_dist" );  }
namespace constraints { StringVectorOptionKey const cst_fa_file( "constraints:cst_fa_file" );  }
namespace constraints { RealOptionKey const cst_fa_weight( "constraints:cst_fa_weight" );  }
namespace constraints { BooleanOptionKey const normalize_mixture_func( "constraints:normalize_mixture_func" );  }
namespace constraints { BooleanOptionKey const penalize_mixture_func( "constraints:penalize_mixture_func" );  }
namespace constraints { FileOptionKey const forest_file( "constraints:forest_file" );  }
namespace constraints { BooleanOptionKey const compute_total_dist_cst( "constraints:compute_total_dist_cst" );  }
namespace constraints { IntegerOptionKey const cull_with_native( "constraints:cull_with_native" );  }
namespace constraints { FileOptionKey const dump_cst_set( "constraints:dump_cst_set" );  }
namespace constraints { IntegerVectorOptionKey const evaluate_max_seq_sep( "constraints:evaluate_max_seq_sep" );  }
namespace constraints { BooleanOptionKey const named( "constraints:named" );  }
namespace constraints { BooleanOptionKey const no_cst_in_relax( "constraints:no_cst_in_relax" );  }
namespace constraints { BooleanOptionKey const no_linearize_bounded( "constraints:no_linearize_bounded" );  }
namespace constraints { RealOptionKey const pocket_constraint_weight( "constraints:pocket_constraint_weight" );  }
namespace constraints { BooleanOptionKey const pocket_zero_derivatives( "constraints:pocket_zero_derivatives" );  }
namespace constraints { BooleanOptionKey const viol( "constraints:viol" );  }
namespace constraints { IntegerOptionKey const viol_level( "constraints:viol_level" );  }
namespace constraints { StringOptionKey const viol_type( "constraints:viol_type" );  }
namespace constraints { RealOptionKey const sog_cst_param( "constraints:sog_cst_param" );  }
namespace constraints { RealOptionKey const sog_upper_bound( "constraints:sog_upper_bound" );  }
namespace constraints { BooleanOptionKey const epr_distance( "constraints:epr_distance" );  }
namespace constraints { IntegerOptionKey const combine( "constraints:combine" );  }
namespace constraints { FileOptionKey const combine_exclude_region( "constraints:combine_exclude_region" );  }
namespace constraints { BooleanOptionKey const skip_redundant( "constraints:skip_redundant" );  }
namespace constraints { IntegerOptionKey const skip_redundant_width( "constraints:skip_redundant_width" );  }
namespace constraints { RealOptionKey const increase_constraints( "constraints:increase_constraints" );  }
namespace dna { BooleanOptionKey const dna( "dna" );  }
namespace dna { namespace specificity { BooleanOptionKey const specificity( "dna:specificity" );  } }
namespace dna { namespace specificity { BooleanOptionKey const exclude_dna_dna( "dna:specificity:exclude_dna_dna" );  } }
namespace dna { namespace specificity { RealVectorOptionKey const params( "dna:specificity:params" );  } }
namespace dna { namespace specificity { FileVectorOptionKey const frag_files( "dna:specificity:frag_files" );  } }
namespace dna { namespace specificity { BooleanOptionKey const exclude_bb_sc_hbonds( "dna:specificity:exclude_bb_sc_hbonds" );  } }
namespace dna { namespace specificity { BooleanOptionKey const only_repack( "dna:specificity:only_repack" );  } }
namespace dna { namespace specificity { BooleanOptionKey const design_DNA( "dna:specificity:design_DNA" );  } }
namespace dna { namespace specificity { BooleanOptionKey const run_test( "dna:specificity:run_test" );  } }
namespace dna { namespace specificity { BooleanOptionKey const soft_rep( "dna:specificity:soft_rep" );  } }
namespace dna { namespace specificity { BooleanOptionKey const dump_pdbs( "dna:specificity:dump_pdbs" );  } }
namespace dna { namespace specificity { BooleanOptionKey const fast( "dna:specificity:fast" );  } }
namespace dna { namespace specificity { BooleanOptionKey const randomize_motif( "dna:specificity:randomize_motif" );  } }
namespace dna { namespace specificity { RealOptionKey const Wfa_elec( "dna:specificity:Wfa_elec" );  } }
namespace dna { namespace specificity { RealOptionKey const Wdna_bs( "dna:specificity:Wdna_bs" );  } }
namespace dna { namespace specificity { RealOptionKey const Wdna_bp( "dna:specificity:Wdna_bp" );  } }
namespace dna { namespace specificity { RealOptionKey const minimize_tolerance( "dna:specificity:minimize_tolerance" );  } }
namespace dna { namespace specificity { StringOptionKey const weights_tag( "dna:specificity:weights_tag" );  } }
namespace dna { namespace specificity { StringOptionKey const weights_tag_list( "dna:specificity:weights_tag_list" );  } }
namespace dna { namespace specificity { StringOptionKey const min_type( "dna:specificity:min_type" );  } }
namespace dna { namespace specificity { StringOptionKey const tf( "dna:specificity:tf" );  } }
namespace dna { namespace specificity { StringOptionKey const mode( "dna:specificity:mode" );  } }
namespace dna { namespace specificity { StringOptionKey const score_function( "dna:specificity:score_function" );  } }
namespace dna { namespace specificity { BooleanOptionKey const pre_minimize( "dna:specificity:pre_minimize" );  } }
namespace dna { namespace specificity { BooleanOptionKey const post_minimize( "dna:specificity:post_minimize" );  } }
namespace dna { namespace specificity { BooleanOptionKey const pre_pack( "dna:specificity:pre_pack" );  } }
namespace dna { namespace specificity { IntegerOptionKey const nloop( "dna:specificity:nloop" );  } }
namespace dna { namespace specificity { IntegerOptionKey const n_inner( "dna:specificity:n_inner" );  } }
namespace dna { namespace specificity { IntegerOptionKey const n_outer( "dna:specificity:n_outer" );  } }
namespace dna { namespace specificity { IntegerOptionKey const nstep_water( "dna:specificity:nstep_water" );  } }
namespace dna { namespace specificity { IntegerOptionKey const moving_jump( "dna:specificity:moving_jump" );  } }
namespace dna { namespace specificity { IntegerOptionKey const motif_begin( "dna:specificity:motif_begin" );  } }
namespace dna { namespace specificity { IntegerOptionKey const motif_size( "dna:specificity:motif_size" );  } }
namespace dna { namespace specificity { StringVectorOptionKey const pdb_pos( "dna:specificity:pdb_pos" );  } }
namespace dna { namespace specificity { StringVectorOptionKey const methylate( "dna:specificity:methylate" );  } }
namespace dna { namespace design { BooleanOptionKey const design( "dna:design" );  } }
namespace dna { namespace design { BooleanOptionKey const output_initial_pdb( "dna:design:output_initial_pdb" );  } }
namespace dna { namespace design { BooleanOptionKey const output_unbound_pdb( "dna:design:output_unbound_pdb" );  } }
namespace dna { namespace design { RealOptionKey const z_cutoff( "dna:design:z_cutoff" );  } }
namespace dna { namespace design { StringOptionKey const protein_scan( "dna:design:protein_scan" );  } }
namespace dna { namespace design { StringOptionKey const checkpoint( "dna:design:checkpoint" );  } }
namespace dna { namespace design { BooleanOptionKey const minimize( "dna:design:minimize" );  } }
namespace dna { namespace design { IntegerOptionKey const iterations( "dna:design:iterations" );  } }
namespace dna { namespace design { StringOptionKey const bb_moves( "dna:design:bb_moves" );  } }
namespace dna { namespace design { StringVectorOptionKey const dna_defs( "dna:design:dna_defs" );  } }
namespace dna { namespace design { StringOptionKey const dna_defs_file( "dna:design:dna_defs_file" );  } }
namespace dna { namespace design { BooleanOptionKey const preminimize_interface( "dna:design:preminimize_interface" );  } }
namespace dna { namespace design { BooleanOptionKey const prepack_interface( "dna:design:prepack_interface" );  } }
namespace dna { namespace design { BooleanOptionKey const flush( "dna:design:flush" );  } }
namespace dna { namespace design { BooleanOptionKey const nopdb( "dna:design:nopdb" );  } }
namespace dna { namespace design { BooleanOptionKey const nopack( "dna:design:nopack" );  } }
namespace dna { namespace design { BooleanOptionKey const more_stats( "dna:design:more_stats" );  } }
namespace dna { namespace design { BooleanOptionKey const pdb_each_iteration( "dna:design:pdb_each_iteration" );  } }
namespace dna { namespace design { BooleanOptionKey const designable_second_shell( "dna:design:designable_second_shell" );  } }
namespace dna { namespace design { BooleanOptionKey const base_contacts_only( "dna:design:base_contacts_only" );  } }
namespace dna { namespace design { IntegerOptionKey const probe_specificity( "dna:design:probe_specificity" );  } }
namespace dna { namespace design { namespace specificity { BooleanOptionKey const specificity( "dna:design:specificity" );  } } }
namespace dna { namespace design { namespace specificity { BooleanOptionKey const output_structures( "dna:design:specificity:output_structures" );  } } }
namespace dna { namespace design { namespace specificity { BooleanOptionKey const include_dna_potentials( "dna:design:specificity:include_dna_potentials" );  } } }
namespace dna { namespace design { BooleanOptionKey const reversion_scan( "dna:design:reversion_scan" );  } }
namespace dna { namespace design { namespace reversion { BooleanOptionKey const reversion( "dna:design:reversion" );  } } }
namespace dna { namespace design { namespace reversion { RealOptionKey const dscore_cutoff( "dna:design:reversion:dscore_cutoff" );  } } }
namespace dna { namespace design { namespace reversion { RealOptionKey const dspec_cutoff( "dna:design:reversion:dspec_cutoff" );  } } }
namespace dna { namespace design { BooleanOptionKey const binding( "dna:design:binding" );  } }
namespace dna { namespace design { RealOptionKey const Boltz_temp( "dna:design:Boltz_temp" );  } }
namespace dna { namespace design { BooleanOptionKey const repack_only( "dna:design:repack_only" );  } }
namespace dna { namespace design { BooleanOptionKey const sparse_pdb_output( "dna:design:sparse_pdb_output" );  } }
namespace flxbb { BooleanOptionKey const flxbb( "flxbb" );  }
namespace flxbb { BooleanOptionKey const view( "flxbb:view" );  }
namespace flxbb { IntegerOptionKey const ncycle( "flxbb:ncycle" );  }
namespace flxbb { RealOptionKey const constraints_sheet( "flxbb:constraints_sheet" );  }
namespace flxbb { BooleanOptionKey const constraints_sheet_include_cacb_pseudotorsion( "flxbb:constraints_sheet_include_cacb_pseudotorsion" );  }
namespace flxbb { RealOptionKey const constraints_NtoC( "flxbb:constraints_NtoC" );  }
namespace flxbb { IntegerOptionKey const filter_trial( "flxbb:filter_trial" );  }
namespace flxbb { StringOptionKey const filter_type( "flxbb:filter_type" );  }
namespace flxbb { BooleanOptionKey const exclude_Met( "flxbb:exclude_Met" );  }
namespace flxbb { BooleanOptionKey const exclude_Ala( "flxbb:exclude_Ala" );  }
namespace flxbb { FileOptionKey const blueprint( "flxbb:blueprint" );  }
namespace flxbb { BooleanOptionKey const movemap_from_blueprint( "flxbb:movemap_from_blueprint" );  }
namespace flxbb { namespace layer { StringOptionKey const layer( "flxbb:layer" );  } }
namespace flxbb { namespace layer { RealOptionKey const pore_radius( "flxbb:layer:pore_radius" );  } }
namespace flxbb { namespace layer { RealOptionKey const burial( "flxbb:layer:burial" );  } }
namespace flxbb { namespace layer { RealOptionKey const surface( "flxbb:layer:surface" );  } }
namespace fldsgn { BooleanOptionKey const fldsgn( "fldsgn" );  }
namespace fldsgn { BooleanOptionKey const view( "fldsgn:view" );  }
namespace fldsgn { FileVectorOptionKey const blueprint( "fldsgn:blueprint" );  }
namespace fldsgn { IntegerOptionKey const dr_cycles( "fldsgn:dr_cycles" );  }
namespace fldsgn { StringOptionKey const centroid_sfx( "fldsgn:centroid_sfx" );  }
namespace fldsgn { StringOptionKey const centroid_sfx_patch( "fldsgn:centroid_sfx_patch" );  }
namespace fldsgn { StringOptionKey const fullatom_sfx( "fldsgn:fullatom_sfx" );  }
namespace fldsgn { StringOptionKey const fullatom_sfx_patch( "fldsgn:fullatom_sfx_patch" );  }
namespace fldsgn { IntegerOptionKey const run_flxbb( "fldsgn:run_flxbb" );  }
namespace rna { BooleanOptionKey const rna( "rna" );  }
namespace rna { IntegerOptionKey const minimize_rounds( "rna:minimize_rounds" );  }
namespace rna { BooleanOptionKey const corrected_geo( "rna:corrected_geo" );  }
namespace rna { BooleanOptionKey const vary_geometry( "rna:vary_geometry" );  }
namespace rna { BooleanOptionKey const skip_coord_constraints( "rna:skip_coord_constraints" );  }
namespace rna { BooleanOptionKey const skip_o2prime_trials( "rna:skip_o2prime_trials" );  }
namespace rna { StringOptionKey const vall_torsions( "rna:vall_torsions" );  }
namespace rna { BooleanOptionKey const jump_database( "rna:jump_database" );  }
namespace rna { BooleanOptionKey const bps_database( "rna:bps_database" );  }
namespace rna { BooleanOptionKey const rna_prot_erraser( "rna:rna_prot_erraser" );  }
namespace rna { BooleanOptionKey const deriv_check( "rna:deriv_check" );  }
namespace cm { BooleanOptionKey const cm( "cm" );  }
namespace cm { namespace sanitize { BooleanOptionKey const sanitize( "cm:sanitize" );  } }
namespace cm { namespace sanitize { RealOptionKey const bound_delta( "cm:sanitize:bound_delta" );  } }
namespace cm { namespace sanitize { RealOptionKey const bound_sd( "cm:sanitize:bound_sd" );  } }
namespace cm { namespace sanitize { IntegerOptionKey const num_fragments( "cm:sanitize:num_fragments" );  } }
namespace cm { namespace sanitize { RealOptionKey const cst_weight_pair( "cm:sanitize:cst_weight_pair" );  } }
namespace cm { namespace sanitize { RealOptionKey const cst_weight_coord( "cm:sanitize:cst_weight_coord" );  } }
namespace cm { BooleanOptionKey const start_models_only( "cm:start_models_only" );  }
namespace cm { StringOptionKey const aln_format( "cm:aln_format" );  }
namespace cm { BooleanOptionKey const recover_side_chains( "cm:recover_side_chains" );  }
namespace cm { FileVectorOptionKey const steal_extra_residues( "cm:steal_extra_residues" );  }
namespace cm { StringOptionKey const loop_mover( "cm:loop_mover" );  }
namespace cm { IntegerOptionKey const loop_close_level( "cm:loop_close_level" );  }
namespace cm { IntegerOptionKey const min_loop_size( "cm:min_loop_size" );  }
namespace cm { IntegerOptionKey const max_loop_rebuild( "cm:max_loop_rebuild" );  }
namespace cm { RealOptionKey const loop_rebuild_filter( "cm:loop_rebuild_filter" );  }
namespace cm { RealOptionKey const aln_length_filter_quantile( "cm:aln_length_filter_quantile" );  }
namespace cm { IntegerOptionKey const aln_length_filter( "cm:aln_length_filter" );  }
namespace cm { StringVectorOptionKey const template_ids( "cm:template_ids" );  }
namespace cm { FileOptionKey const ligand_pdb( "cm:ligand_pdb" );  }
namespace cm { StringVectorOptionKey const seq_score( "cm:seq_score" );  }
namespace cm { StringOptionKey const aligner( "cm:aligner" );  }
namespace cm { RealOptionKey const min_gap_open( "cm:min_gap_open" );  }
namespace cm { RealOptionKey const max_gap_open( "cm:max_gap_open" );  }
namespace cm { RealOptionKey const min_gap_extend( "cm:min_gap_extend" );  }
namespace cm { RealOptionKey const max_gap_extend( "cm:max_gap_extend" );  }
namespace cm { IntegerOptionKey const nn( "cm:nn" );  }
namespace cm { RealOptionKey const fr_temperature( "cm:fr_temperature" );  }
namespace cm { FileVectorOptionKey const ev_map( "cm:ev_map" );  }
namespace cm { FileVectorOptionKey const hh_map( "cm:hh_map" );  }
namespace cm { namespace hybridize { BooleanOptionKey const hybridize( "cm:hybridize" );  } }
namespace cm { namespace hybridize { FileVectorOptionKey const templates( "cm:hybridize:templates" );  } }
namespace cm { namespace hybridize { FileOptionKey const template_list( "cm:hybridize:template_list" );  } }
namespace cm { namespace hybridize { IntegerVectorOptionKey const starting_template( "cm:hybridize:starting_template" );  } }
namespace cm { namespace hybridize { BooleanOptionKey const realign_domains( "cm:hybridize:realign_domains" );  } }
namespace cm { namespace hybridize { BooleanOptionKey const add_non_init_chunks( "cm:hybridize:add_non_init_chunks" );  } }
namespace cm { namespace hybridize { StringOptionKey const ss( "cm:hybridize:ss" );  } }
namespace cm { namespace hybridize { RealOptionKey const stage1_increase_cycles( "cm:hybridize:stage1_increase_cycles" );  } }
namespace cm { namespace hybridize { RealOptionKey const stage2_increase_cycles( "cm:hybridize:stage2_increase_cycles" );  } }
namespace cm { namespace hybridize { RealOptionKey const stage2min_increase_cycles( "cm:hybridize:stage2min_increase_cycles" );  } }
namespace cm { namespace hybridize { RealOptionKey const stage1_probability( "cm:hybridize:stage1_probability" );  } }
namespace cm { namespace hybridize { StringOptionKey const stage1_weights( "cm:hybridize:stage1_weights" );  } }
namespace cm { namespace hybridize { StringOptionKey const stage1_patch( "cm:hybridize:stage1_patch" );  } }
namespace cm { namespace hybridize { BooleanOptionKey const skip_stage2( "cm:hybridize:skip_stage2" );  } }
namespace cm { namespace hybridize { BooleanOptionKey const no_global_frame( "cm:hybridize:no_global_frame" );  } }
namespace cm { namespace hybridize { BooleanOptionKey const linmin_only( "cm:hybridize:linmin_only" );  } }
namespace cm { namespace hybridize { StringOptionKey const stage2_weights( "cm:hybridize:stage2_weights" );  } }
namespace cm { namespace hybridize { StringOptionKey const stage2_patch( "cm:hybridize:stage2_patch" );  } }
namespace cm { namespace hybridize { IntegerOptionKey const relax( "cm:hybridize:relax" );  } }
namespace cm { namespace hybridize { RealOptionKey const frag_weight_aligned( "cm:hybridize:frag_weight_aligned" );  } }
namespace cm { namespace hybridize { BooleanOptionKey const move_anchor( "cm:hybridize:move_anchor" );  } }
namespace cm { namespace hybridize { IntegerOptionKey const max_registry_shift( "cm:hybridize:max_registry_shift" );  } }
namespace cm { namespace hybridize { BooleanOptionKey const alignment_from_template_seqpos( "cm:hybridize:alignment_from_template_seqpos" );  } }
namespace cm { namespace hybridize { IntegerVectorOptionKey const alignment_from_chunk_mapping( "cm:hybridize:alignment_from_chunk_mapping" );  } }
namespace cm { namespace hybridize { BooleanOptionKey const virtual_loops( "cm:hybridize:virtual_loops" );  } }
namespace cm { namespace hybridize { BooleanOptionKey const revert_real_loops( "cm:hybridize:revert_real_loops" );  } }
namespace cm { namespace hybridize { BooleanOptionKey const realign_domains_stage2( "cm:hybridize:realign_domains_stage2" );  } }
namespace cm { namespace hybridize { RealOptionKey const frag_1mer_insertion_weight( "cm:hybridize:frag_1mer_insertion_weight" );  } }
namespace cm { namespace hybridize { RealOptionKey const small_frag_insertion_weight( "cm:hybridize:small_frag_insertion_weight" );  } }
namespace cm { namespace hybridize { RealOptionKey const big_frag_insertion_weight( "cm:hybridize:big_frag_insertion_weight" );  } }
namespace cm { namespace hybridize { BooleanOptionKey const auto_frag_insertion_weight( "cm:hybridize:auto_frag_insertion_weight" );  } }
namespace cm { namespace hybridize { IntegerOptionKey const stage1_1_cycles( "cm:hybridize:stage1_1_cycles" );  } }
namespace cm { namespace hybridize { IntegerOptionKey const stage1_2_cycles( "cm:hybridize:stage1_2_cycles" );  } }
namespace cm { namespace hybridize { IntegerOptionKey const stage1_3_cycles( "cm:hybridize:stage1_3_cycles" );  } }
namespace cm { namespace hybridize { IntegerOptionKey const stage1_4_cycles( "cm:hybridize:stage1_4_cycles" );  } }
namespace cm { namespace hybridize { RealOptionKey const stage2_temperature( "cm:hybridize:stage2_temperature" );  } }
namespace cm { namespace hybridize { StringOptionKey const stage1_4_cenrot_score( "cm:hybridize:stage1_4_cenrot_score" );  } }
namespace ms { BooleanOptionKey const ms( "ms" );  }
namespace ms { BooleanOptionKey const share_data( "ms:share_data" );  }
namespace ms { BooleanOptionKey const verbose( "ms:verbose" );  }
namespace ms { BooleanOptionKey const restrict_to_canonical( "ms:restrict_to_canonical" );  }
namespace ms { IntegerOptionKey const pop_from_ss( "ms:pop_from_ss" );  }
namespace ms { IntegerOptionKey const pop_size( "ms:pop_size" );  }
namespace ms { IntegerOptionKey const generations( "ms:generations" );  }
namespace ms { IntegerOptionKey const num_packs( "ms:num_packs" );  }
namespace ms { IntegerOptionKey const numresults( "ms:numresults" );  }
namespace ms { RealOptionKey const anchor_offset( "ms:anchor_offset" );  }
namespace ms { RealOptionKey const Boltz_temp( "ms:Boltz_temp" );  }
namespace ms { RealOptionKey const mutate_rate( "ms:mutate_rate" );  }
namespace ms { RealOptionKey const fraction_by_recombination( "ms:fraction_by_recombination" );  }
namespace ms { namespace checkpoint { BooleanOptionKey const checkpoint( "ms:checkpoint" );  } }
namespace ms { namespace checkpoint { StringOptionKey const prefix( "ms:checkpoint:prefix" );  } }
namespace ms { namespace checkpoint { IntegerOptionKey const interval( "ms:checkpoint:interval" );  } }
namespace ms { namespace checkpoint { BooleanOptionKey const gz( "ms:checkpoint:gz" );  } }
namespace ms { namespace checkpoint { BooleanOptionKey const rename( "ms:checkpoint:rename" );  } }
namespace loops { BooleanOptionKey const loops( "loops" );  }
namespace loops { StringOptionKey const cen_weights( "loops:cen_weights" );  }
namespace loops { StringOptionKey const cen_patch( "loops:cen_patch" );  }
namespace loops { FileOptionKey const input_pdb( "loops:input_pdb" );  }
namespace loops { StringVectorOptionKey const loop_file( "loops:loop_file" );  }
namespace loops { FileOptionKey const extended_loop_file( "loops:extended_loop_file" );  }
namespace loops { FileOptionKey const mm_loop_file( "loops:mm_loop_file" );  }
namespace loops { BooleanOptionKey const fix_natsc( "loops:fix_natsc" );  }
namespace loops { BooleanOptionKey const refine_only( "loops:refine_only" );  }
namespace loops { BooleanOptionKey const fa_input( "loops:fa_input" );  }
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
namespace loops { StringOptionKey const remodel( "loops:remodel" );  }
namespace loops { StringOptionKey const intermedrelax( "loops:intermedrelax" );  }
namespace loops { StringOptionKey const refine( "loops:refine" );  }
namespace loops { StringOptionKey const relax( "loops:relax" );  }
namespace loops { IntegerOptionKey const n_rebuild_tries( "loops:n_rebuild_tries" );  }
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
namespace loops { IntegerOptionKey const outer_cycles( "loops:outer_cycles" );  }
namespace loops { IntegerOptionKey const max_inner_cycles( "loops:max_inner_cycles" );  }
namespace loops { IntegerOptionKey const repack_period( "loops:repack_period" );  }
namespace loops { BooleanOptionKey const repack_never( "loops:repack_never" );  }
namespace loops { RealOptionKey const remodel_init_temp( "loops:remodel_init_temp" );  }
namespace loops { RealOptionKey const remodel_final_temp( "loops:remodel_final_temp" );  }
namespace loops { RealOptionKey const refine_init_temp( "loops:refine_init_temp" );  }
namespace loops { RealOptionKey const refine_final_temp( "loops:refine_final_temp" );  }
namespace loops { IntegerOptionKey const gapspan( "loops:gapspan" );  }
namespace loops { IntegerOptionKey const spread( "loops:spread" );  }
namespace loops { IntegerOptionKey const kinematic_wrapper_cycles( "loops:kinematic_wrapper_cycles" );  }
namespace loops { IntegerOptionKey const kic_num_rotamer_trials( "loops:kic_num_rotamer_trials" );  }
namespace loops { BooleanOptionKey const kic_omega_sampling( "loops:kic_omega_sampling" );  }
namespace loops { RealOptionKey const kic_bump_overlap_factor( "loops:kic_bump_overlap_factor" );  }
namespace loops { StringOptionKey const kic_cen_weights( "loops:kic_cen_weights" );  }
namespace loops { StringOptionKey const kic_cen_patch( "loops:kic_cen_patch" );  }
namespace loops { StringOptionKey const restrict_kic_sampling_to_torsion_string( "loops:restrict_kic_sampling_to_torsion_string" );  }
namespace loops { BooleanOptionKey const derive_torsion_string_from_native_pose( "loops:derive_torsion_string_from_native_pose" );  }
namespace loops { BooleanOptionKey const always_remodel_full_loop( "loops:always_remodel_full_loop" );  }
namespace loops { BooleanOptionKey const taboo_sampling( "loops:taboo_sampling" );  }
namespace loops { BooleanOptionKey const taboo_in_fa( "loops:taboo_in_fa" );  }
namespace loops { BooleanOptionKey const ramp_fa_rep( "loops:ramp_fa_rep" );  }
namespace loops { BooleanOptionKey const ramp_rama( "loops:ramp_rama" );  }
namespace loops { BooleanOptionKey const kic_rama2b( "loops:kic_rama2b" );  }
namespace loops { BooleanOptionKey const kic_small_moves( "loops:kic_small_moves" );  }
namespace loops { RealOptionKey const kic_small_move_magnitude( "loops:kic_small_move_magnitude" );  }
namespace loops { BooleanOptionKey const kic_pivot_based( "loops:kic_pivot_based" );  }
namespace loops { BooleanOptionKey const kic_no_centroid_min( "loops:kic_no_centroid_min" );  }
namespace loops { BooleanOptionKey const kic_leave_centroid_after_initial_closure( "loops:kic_leave_centroid_after_initial_closure" );  }
namespace loops { BooleanOptionKey const kic_repack_neighbors_only( "loops:kic_repack_neighbors_only" );  }
namespace loops { BooleanOptionKey const legacy_kic( "loops:legacy_kic" );  }
namespace loops { BooleanOptionKey const alternative_closure_protocol( "loops:alternative_closure_protocol" );  }
namespace loops { RealOptionKey const chainbreak_max_accept( "loops:chainbreak_max_accept" );  }
namespace loops { BooleanOptionKey const debug_loop_closure( "loops:debug_loop_closure" );  }
namespace loops { BooleanOptionKey const non_ideal_loop_closing( "loops:non_ideal_loop_closing" );  }
namespace loops { RealOptionKey const scored_frag_cycles( "loops:scored_frag_cycles" );  }
namespace loops { RealOptionKey const short_frag_cycles( "loops:short_frag_cycles" );  }
namespace loops { RealOptionKey const rmsd_tol( "loops:rmsd_tol" );  }
namespace loops { RealOptionKey const chain_break_tol( "loops:chain_break_tol" );  }
namespace loops { BooleanOptionKey const random_loop( "loops:random_loop" );  }
namespace loops { FileVectorOptionKey const stealfrags( "loops:stealfrags" );  }
namespace loops { IntegerOptionKey const stealfrags_times( "loops:stealfrags_times" );  }
namespace loops { RealOptionKey const coord_cst( "loops:coord_cst" );  }
namespace loops { RealOptionKey const skip_1mers( "loops:skip_1mers" );  }
namespace loops { RealOptionKey const skip_3mers( "loops:skip_3mers" );  }
namespace loops { RealOptionKey const skip_9mers( "loops:skip_9mers" );  }
namespace loops { BooleanOptionKey const loop_model( "loops:loop_model" );  }
namespace loops { RealOptionKey const score_filter_cutoff( "loops:score_filter_cutoff" );  }
namespace loops { BooleanOptionKey const loop_farlx( "loops:loop_farlx" );  }
namespace loops { BooleanOptionKey const ccd_closure( "loops:ccd_closure" );  }
namespace loops { BooleanOptionKey const skip_ccd_moves( "loops:skip_ccd_moves" );  }
namespace loops { BooleanOptionKey const no_randomize_loop( "loops:no_randomize_loop" );  }
namespace loops { BooleanOptionKey const loops_subset( "loops:loops_subset" );  }
namespace loops { IntegerOptionKey const num_desired_loops( "loops:num_desired_loops" );  }
namespace loops { RealOptionKey const loop_combine_rate( "loops:loop_combine_rate" );  }
namespace loops { RealOptionKey const final_score_filter( "loops:final_score_filter" );  }
namespace loops { BooleanOptionKey const no_combine_if_fail( "loops:no_combine_if_fail" );  }
namespace loops { BooleanOptionKey const shorten_long_terminal_loop( "loops:shorten_long_terminal_loop" );  }
namespace loops { IntegerOptionKey const backrub_trials( "loops:backrub_trials" );  }
namespace loops { RealOptionKey const looprlx_cycle_ratio( "loops:looprlx_cycle_ratio" );  }
namespace loops { RealOptionKey const extended_beta( "loops:extended_beta" );  }
namespace loops { BooleanOptionKey const shortrelax( "loops:shortrelax" );  }
namespace loops { BooleanOptionKey const fastrelax( "loops:fastrelax" );  }
namespace loops { BooleanOptionKey const no_looprebuild( "loops:no_looprebuild" );  }
namespace loops { BooleanOptionKey const allow_lig_move( "loops:allow_lig_move" );  }
namespace loops { FileOptionKey const keep_natro( "loops:keep_natro" );  }
namespace loops { IntegerOptionKey const refine_design_iterations( "loops:refine_design_iterations" );  }
namespace loops { namespace loop_closure { BooleanOptionKey const loop_closure( "loops:loop_closure" );  } }
namespace loops { namespace loop_closure { StringOptionKey const loop_insert( "loops:loop_closure:loop_insert" );  } }
namespace loops { namespace loop_closure { StringOptionKey const loop_insert_rclrc( "loops:loop_closure:loop_insert_rclrc" );  } }
namespace loops { namespace loop_closure { StringOptionKey const blueprint( "loops:loop_closure:blueprint" );  } }
namespace assembly { BooleanOptionKey const assembly( "assembly" );  }
namespace assembly { FileOptionKey const pdb1( "assembly:pdb1" );  }
namespace assembly { FileOptionKey const pdb2( "assembly:pdb2" );  }
namespace assembly { StringOptionKey const nterm_seq( "assembly:nterm_seq" );  }
namespace assembly { StringOptionKey const cterm_seq( "assembly:cterm_seq" );  }
namespace assembly { IntegerVectorOptionKey const linkers_pdb1( "assembly:linkers_pdb1" );  }
namespace assembly { IntegerVectorOptionKey const linkers_pdb2( "assembly:linkers_pdb2" );  }
namespace assembly { IntegerVectorOptionKey const preserve_sidechains_pdb1( "assembly:preserve_sidechains_pdb1" );  }
namespace assembly { IntegerVectorOptionKey const preserve_sidechains_pdb2( "assembly:preserve_sidechains_pdb2" );  }
namespace fast_loops { BooleanOptionKey const fast_loops( "fast_loops" );  }
namespace fast_loops { RealOptionKey const window_accept_ratio( "fast_loops:window_accept_ratio" );  }
namespace fast_loops { IntegerOptionKey const nr_scored_sampling_passes( "fast_loops:nr_scored_sampling_passes" );  }
namespace fast_loops { IntegerOptionKey const nr_scored_fragments( "fast_loops:nr_scored_fragments" );  }
namespace fast_loops { IntegerOptionKey const min_breakout_good_loops( "fast_loops:min_breakout_good_loops" );  }
namespace fast_loops { IntegerOptionKey const min_breakout_fast_loops( "fast_loops:min_breakout_fast_loops" );  }
namespace fast_loops { IntegerOptionKey const min_good_loops( "fast_loops:min_good_loops" );  }
namespace fast_loops { IntegerOptionKey const min_fast_loops( "fast_loops:min_fast_loops" );  }
namespace fast_loops { RealOptionKey const vdw_delta( "fast_loops:vdw_delta" );  }
namespace fast_loops { IntegerOptionKey const give_up( "fast_loops:give_up" );  }
namespace fast_loops { RealOptionKey const chainbreak_max( "fast_loops:chainbreak_max" );  }
namespace fast_loops { FileOptionKey const fragsample_score( "fast_loops:fragsample_score" );  }
namespace fast_loops { FileOptionKey const fragsample_patch( "fast_loops:fragsample_patch" );  }
namespace fast_loops { FileOptionKey const overwrite_filter_scorefxn( "fast_loops:overwrite_filter_scorefxn" );  }
namespace fast_loops { FileOptionKey const patch_filter_scorefxn( "fast_loops:patch_filter_scorefxn" );  }
namespace fast_loops { FileOptionKey const filter_cst_file( "fast_loops:filter_cst_file" );  }
namespace fast_loops { RealOptionKey const filter_cst_weight( "fast_loops:filter_cst_weight" );  }
namespace fast_loops { FileOptionKey const fast_relax_sequence_file( "fast_loops:fast_relax_sequence_file" );  }
namespace SSrbrelax { BooleanOptionKey const SSrbrelax( "SSrbrelax" );  }
namespace SSrbrelax { FileOptionKey const input_pdb( "SSrbrelax:input_pdb" );  }
namespace SSrbrelax { FileOptionKey const rb_file( "SSrbrelax:rb_file" );  }
namespace SSrbrelax { FileOptionKey const rb_param_file( "SSrbrelax:rb_param_file" );  }
namespace SSrbrelax { IntegerVectorOptionKey const frag_sizes( "SSrbrelax:frag_sizes" );  }
namespace SSrbrelax { FileVectorOptionKey const frag_files( "SSrbrelax:frag_files" );  }
namespace boinc { BooleanOptionKey const boinc( "boinc" );  }
namespace boinc { BooleanOptionKey const graphics( "boinc:graphics" );  }
namespace boinc { BooleanOptionKey const fullscreen( "boinc:fullscreen" );  }
namespace boinc { IntegerOptionKey const max_fps( "boinc:max_fps" );  }
namespace boinc { IntegerOptionKey const max_cpu( "boinc:max_cpu" );  }
namespace boinc { BooleanOptionKey const noshmem( "boinc:noshmem" );  }
namespace boinc { IntegerOptionKey const cpu_run_time( "boinc:cpu_run_time" );  }
namespace boinc { IntegerOptionKey const max_nstruct( "boinc:max_nstruct" );  }
namespace boinc { RealOptionKey const cpu_frac( "boinc:cpu_frac" );  }
namespace boinc { RealOptionKey const frame_rate( "boinc:frame_rate" );  }
namespace boinc { BooleanOptionKey const watchdog( "boinc:watchdog" );  }
namespace boinc { IntegerOptionKey const watchdog_time( "boinc:watchdog_time" );  }
namespace boinc { IntegerOptionKey const cpu_run_timeout( "boinc:cpu_run_timeout" );  }
namespace boinc { FileOptionKey const description_file( "boinc:description_file" );  }
namespace LoopModel { BooleanOptionKey const LoopModel( "LoopModel" );  }
namespace LoopModel { FileOptionKey const input_pdb( "LoopModel:input_pdb" );  }
namespace LoopModel { FileOptionKey const loop_file( "LoopModel:loop_file" );  }
namespace AnchoredDesign { BooleanOptionKey const AnchoredDesign( "AnchoredDesign" );  }
namespace AnchoredDesign { FileOptionKey const anchor( "AnchoredDesign:anchor" );  }
namespace AnchoredDesign { BooleanOptionKey const allow_anchor_repack( "AnchoredDesign:allow_anchor_repack" );  }
namespace AnchoredDesign { BooleanOptionKey const vary_cutpoints( "AnchoredDesign:vary_cutpoints" );  }
namespace AnchoredDesign { BooleanOptionKey const no_frags( "AnchoredDesign:no_frags" );  }
namespace AnchoredDesign { BooleanOptionKey const debug( "AnchoredDesign:debug" );  }
namespace AnchoredDesign { BooleanOptionKey const show_extended( "AnchoredDesign:show_extended" );  }
namespace AnchoredDesign { BooleanOptionKey const refine_only( "AnchoredDesign:refine_only" );  }
namespace AnchoredDesign { BooleanOptionKey const perturb_show( "AnchoredDesign:perturb_show" );  }
namespace AnchoredDesign { IntegerOptionKey const perturb_cycles( "AnchoredDesign:perturb_cycles" );  }
namespace AnchoredDesign { RealOptionKey const perturb_temp( "AnchoredDesign:perturb_temp" );  }
namespace AnchoredDesign { BooleanOptionKey const perturb_CCD_off( "AnchoredDesign:perturb_CCD_off" );  }
namespace AnchoredDesign { BooleanOptionKey const perturb_KIC_off( "AnchoredDesign:perturb_KIC_off" );  }
namespace AnchoredDesign { BooleanOptionKey const refine_CCD_off( "AnchoredDesign:refine_CCD_off" );  }
namespace AnchoredDesign { BooleanOptionKey const refine_KIC_off( "AnchoredDesign:refine_KIC_off" );  }
namespace AnchoredDesign { IntegerOptionKey const refine_cycles( "AnchoredDesign:refine_cycles" );  }
namespace AnchoredDesign { RealOptionKey const refine_temp( "AnchoredDesign:refine_temp" );  }
namespace AnchoredDesign { IntegerOptionKey const refine_repack_cycles( "AnchoredDesign:refine_repack_cycles" );  }
namespace AnchoredDesign { BooleanOptionKey const rmsd( "AnchoredDesign:rmsd" );  }
namespace AnchoredDesign { BooleanOptionKey const unbound_mode( "AnchoredDesign:unbound_mode" );  }
namespace AnchoredDesign { RealOptionKey const chainbreak_weight( "AnchoredDesign:chainbreak_weight" );  }
namespace AnchoredDesign { namespace filters { BooleanOptionKey const filters( "AnchoredDesign:filters" );  } }
namespace AnchoredDesign { namespace filters { RealOptionKey const score( "AnchoredDesign:filters:score" );  } }
namespace AnchoredDesign { namespace filters { RealOptionKey const sasa( "AnchoredDesign:filters:sasa" );  } }
namespace AnchoredDesign { namespace filters { BooleanOptionKey const omega( "AnchoredDesign:filters:omega" );  } }
namespace AnchoredDesign { namespace akash { BooleanOptionKey const akash( "AnchoredDesign:akash" );  } }
namespace AnchoredDesign { namespace akash { IntegerOptionKey const dyepos( "AnchoredDesign:akash:dyepos" );  } }
namespace AnchoredDesign { namespace testing { BooleanOptionKey const testing( "AnchoredDesign:testing" );  } }
namespace AnchoredDesign { namespace testing { RealOptionKey const VDW_weight( "AnchoredDesign:testing:VDW_weight" );  } }
namespace AnchoredDesign { namespace testing { BooleanOptionKey const anchor_via_constraints( "AnchoredDesign:testing:anchor_via_constraints" );  } }
namespace AnchoredDesign { namespace testing { BooleanOptionKey const delete_interface_native_sidechains( "AnchoredDesign:testing:delete_interface_native_sidechains" );  } }
namespace AnchoredDesign { namespace testing { FileOptionKey const RMSD_only_this( "AnchoredDesign:testing:RMSD_only_this" );  } }
namespace AnchoredDesign { namespace testing { BooleanOptionKey const anchor_noise_constraints_mode( "AnchoredDesign:testing:anchor_noise_constraints_mode" );  } }
namespace AnchoredDesign { namespace testing { BooleanOptionKey const super_secret_fixed_interface_mode( "AnchoredDesign:testing:super_secret_fixed_interface_mode" );  } }
namespace AnchoredDesign { namespace testing { BooleanOptionKey const randomize_input_sequence( "AnchoredDesign:testing:randomize_input_sequence" );  } }
namespace chemically_conjugated_docking { BooleanOptionKey const chemically_conjugated_docking( "chemically_conjugated_docking" );  }
namespace chemically_conjugated_docking { FileOptionKey const UBQpdb( "chemically_conjugated_docking:UBQpdb" );  }
namespace chemically_conjugated_docking { FileOptionKey const E2pdb( "chemically_conjugated_docking:E2pdb" );  }
namespace chemically_conjugated_docking { IntegerOptionKey const E2_residue( "chemically_conjugated_docking:E2_residue" );  }
namespace chemically_conjugated_docking { RealOptionKey const SASAfilter( "chemically_conjugated_docking:SASAfilter" );  }
namespace chemically_conjugated_docking { RealOptionKey const scorefilter( "chemically_conjugated_docking:scorefilter" );  }
namespace chemically_conjugated_docking { BooleanOptionKey const publication( "chemically_conjugated_docking:publication" );  }
namespace chemically_conjugated_docking { IntegerOptionKey const n_tail_res( "chemically_conjugated_docking:n_tail_res" );  }
namespace chemically_conjugated_docking { BooleanOptionKey const two_ubiquitins( "chemically_conjugated_docking:two_ubiquitins" );  }
namespace chemically_conjugated_docking { FileVectorOptionKey const extra_bodies( "chemically_conjugated_docking:extra_bodies" );  }
namespace chemically_conjugated_docking { IntegerOptionKey const UBQ2_lys( "chemically_conjugated_docking:UBQ2_lys" );  }
namespace chemically_conjugated_docking { FileOptionKey const UBQ2_pdb( "chemically_conjugated_docking:UBQ2_pdb" );  }
namespace chemically_conjugated_docking { BooleanOptionKey const dont_minimize_omega( "chemically_conjugated_docking:dont_minimize_omega" );  }
namespace chemically_conjugated_docking { BooleanOptionKey const pdz( "chemically_conjugated_docking:pdz" );  }
namespace chemically_conjugated_docking { FileOptionKey const GTPasepdb( "chemically_conjugated_docking:GTPasepdb" );  }
namespace chemically_conjugated_docking { IntegerOptionKey const GTPase_residue( "chemically_conjugated_docking:GTPase_residue" );  }
namespace FloppyTail { BooleanOptionKey const FloppyTail( "FloppyTail" );  }
namespace FloppyTail { IntegerOptionKey const flexible_start_resnum( "FloppyTail:flexible_start_resnum" );  }
namespace FloppyTail { IntegerOptionKey const flexible_stop_resnum( "FloppyTail:flexible_stop_resnum" );  }
namespace FloppyTail { StringOptionKey const flexible_chain( "FloppyTail:flexible_chain" );  }
namespace FloppyTail { RealOptionKey const shear_on( "FloppyTail:shear_on" );  }
namespace FloppyTail { namespace short_tail { BooleanOptionKey const short_tail( "FloppyTail:short_tail" );  } }
namespace FloppyTail { namespace short_tail { RealOptionKey const short_tail_fraction( "FloppyTail:short_tail:short_tail_fraction" );  } }
namespace FloppyTail { namespace short_tail { RealOptionKey const short_tail_off( "FloppyTail:short_tail:short_tail_off" );  } }
namespace FloppyTail { BooleanOptionKey const pair_off( "FloppyTail:pair_off" );  }
namespace FloppyTail { BooleanOptionKey const publication( "FloppyTail:publication" );  }
namespace FloppyTail { BooleanOptionKey const C_root( "FloppyTail:C_root" );  }
namespace FloppyTail { BooleanOptionKey const force_linear_fold_tree( "FloppyTail:force_linear_fold_tree" );  }
namespace FloppyTail { BooleanOptionKey const debug( "FloppyTail:debug" );  }
namespace FloppyTail { StringOptionKey const cen_weights( "FloppyTail:cen_weights" );  }
namespace FloppyTail { BooleanOptionKey const perturb_show( "FloppyTail:perturb_show" );  }
namespace FloppyTail { IntegerOptionKey const perturb_cycles( "FloppyTail:perturb_cycles" );  }
namespace FloppyTail { RealOptionKey const perturb_temp( "FloppyTail:perturb_temp" );  }
namespace FloppyTail { IntegerOptionKey const refine_cycles( "FloppyTail:refine_cycles" );  }
namespace FloppyTail { RealOptionKey const refine_temp( "FloppyTail:refine_temp" );  }
namespace FloppyTail { IntegerOptionKey const refine_repack_cycles( "FloppyTail:refine_repack_cycles" );  }
namespace DenovoProteinDesign { BooleanOptionKey const DenovoProteinDesign( "DenovoProteinDesign" );  }
namespace DenovoProteinDesign { BooleanOptionKey const redesign_core( "DenovoProteinDesign:redesign_core" );  }
namespace DenovoProteinDesign { BooleanOptionKey const redesign_loops( "DenovoProteinDesign:redesign_loops" );  }
namespace DenovoProteinDesign { BooleanOptionKey const redesign_surface( "DenovoProteinDesign:redesign_surface" );  }
namespace DenovoProteinDesign { BooleanOptionKey const redesign_complete( "DenovoProteinDesign:redesign_complete" );  }
namespace DenovoProteinDesign { BooleanOptionKey const disallow_native_aa( "DenovoProteinDesign:disallow_native_aa" );  }
namespace DenovoProteinDesign { BooleanOptionKey const optimize_loops( "DenovoProteinDesign:optimize_loops" );  }
namespace DenovoProteinDesign { FileOptionKey const secondary_structure_file( "DenovoProteinDesign:secondary_structure_file" );  }
namespace DenovoProteinDesign { FileOptionKey const hydrophobic_polar_pattern( "DenovoProteinDesign:hydrophobic_polar_pattern" );  }
namespace DenovoProteinDesign { BooleanOptionKey const use_template_sequence( "DenovoProteinDesign:use_template_sequence" );  }
namespace DenovoProteinDesign { BooleanOptionKey const use_template_topology( "DenovoProteinDesign:use_template_topology" );  }
namespace DenovoProteinDesign { FileOptionKey const create_from_template_pdb( "DenovoProteinDesign:create_from_template_pdb" );  }
namespace DenovoProteinDesign { BooleanOptionKey const create_from_secondary_structure( "DenovoProteinDesign:create_from_secondary_structure" );  }
namespace RBSegmentRelax { BooleanOptionKey const RBSegmentRelax( "RBSegmentRelax" );  }
namespace RBSegmentRelax { FileOptionKey const input_pdb( "RBSegmentRelax:input_pdb" );  }
namespace RBSegmentRelax { FileOptionKey const rb_file( "RBSegmentRelax:rb_file" );  }
namespace RBSegmentRelax { RealOptionKey const cst_wt( "RBSegmentRelax:cst_wt" );  }
namespace RBSegmentRelax { RealOptionKey const cst_width( "RBSegmentRelax:cst_width" );  }
namespace RBSegmentRelax { StringOptionKey const cst_pdb( "RBSegmentRelax:cst_pdb" );  }
namespace RBSegmentRelax { IntegerOptionKey const nrbmoves( "RBSegmentRelax:nrbmoves" );  }
namespace RBSegmentRelax { IntegerOptionKey const nrboutercycles( "RBSegmentRelax:nrboutercycles" );  }
namespace RBSegmentRelax { StringOptionKey const rb_scorefxn( "RBSegmentRelax:rb_scorefxn" );  }
namespace RBSegmentRelax { BooleanOptionKey const skip_fragment_moves( "RBSegmentRelax:skip_fragment_moves" );  }
namespace RBSegmentRelax { BooleanOptionKey const skip_seqshift_moves( "RBSegmentRelax:skip_seqshift_moves" );  }
namespace RBSegmentRelax { BooleanOptionKey const skip_rb_moves( "RBSegmentRelax:skip_rb_moves" );  }
namespace RBSegmentRelax { RealVectorOptionKey const helical_movement_params( "RBSegmentRelax:helical_movement_params" );  }
namespace RBSegmentRelax { RealVectorOptionKey const strand_movement_params( "RBSegmentRelax:strand_movement_params" );  }
namespace RBSegmentRelax { RealVectorOptionKey const default_movement_params( "RBSegmentRelax:default_movement_params" );  }
namespace RBSegmentRelax { IntegerOptionKey const cst_seqwidth( "RBSegmentRelax:cst_seqwidth" );  }
namespace edensity { BooleanOptionKey const edensity( "edensity" );  }
namespace edensity { BooleanOptionKey const debug( "edensity:debug" );  }
namespace edensity { StringOptionKey const mapfile( "edensity:mapfile" );  }
namespace edensity { RealOptionKey const mapreso( "edensity:mapreso" );  }
namespace edensity { RealOptionKey const grid_spacing( "edensity:grid_spacing" );  }
namespace edensity { RealOptionKey const centroid_density_mass( "edensity:centroid_density_mass" );  }
namespace edensity { IntegerOptionKey const sliding_window( "edensity:sliding_window" );  }
namespace edensity { BooleanOptionKey const cryoem_scatterers( "edensity:cryoem_scatterers" );  }
namespace edensity { RealOptionKey const force_apix( "edensity:force_apix" );  }
namespace edensity { RealOptionKey const fastdens_wt( "edensity:fastdens_wt" );  }
namespace edensity { RealVectorOptionKey const fastdens_params( "edensity:fastdens_params" );  }
namespace edensity { BooleanOptionKey const legacy_fastdens_score( "edensity:legacy_fastdens_score" );  }
namespace edensity { RealOptionKey const sliding_window_wt( "edensity:sliding_window_wt" );  }
namespace edensity { BooleanOptionKey const score_sliding_window_context( "edensity:score_sliding_window_context" );  }
namespace edensity { RealOptionKey const whole_structure_ca_wt( "edensity:whole_structure_ca_wt" );  }
namespace edensity { RealOptionKey const whole_structure_allatom_wt( "edensity:whole_structure_allatom_wt" );  }
namespace edensity { BooleanOptionKey const no_edens_in_minimizer( "edensity:no_edens_in_minimizer" );  }
namespace edensity { BooleanOptionKey const debug_derivatives( "edensity:debug_derivatives" );  }
namespace edensity { StringOptionKey const realign( "edensity:realign" );  }
namespace edensity { StringOptionKey const membrane_axis( "edensity:membrane_axis" );  }
namespace edensity { RealOptionKey const atom_mask( "edensity:atom_mask" );  }
namespace edensity { RealOptionKey const atom_mask_min( "edensity:atom_mask_min" );  }
namespace edensity { RealOptionKey const ca_mask( "edensity:ca_mask" );  }
namespace edensity { BooleanOptionKey const score_symm_complex( "edensity:score_symm_complex" );  }
namespace edensity { RealOptionKey const sc_scaling( "edensity:sc_scaling" );  }
namespace edensity { IntegerOptionKey const n_kbins( "edensity:n_kbins" );  }
namespace edensity { RealOptionKey const render_sigma( "edensity:render_sigma" );  }
namespace edensity { BooleanOptionKey const unmask_bb( "edensity:unmask_bb" );  }
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
namespace optE { BooleanOptionKey const optE_no_protein_fa_elec( "optE:optE_no_protein_fa_elec" );  }
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
