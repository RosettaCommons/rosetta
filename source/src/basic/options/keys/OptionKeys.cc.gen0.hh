namespace in { BooleanOptionKey const in( "in" );  }
namespace in { StringOptionKey const Ntermini( "in:Ntermini" );  }
namespace in { StringOptionKey const Ctermini( "in:Ctermini" );  }
namespace in { BooleanOptionKey const use_truncated_termini( "in:use_truncated_termini" );  }
namespace in { BooleanOptionKey const ignore_unrecognized_res( "in:ignore_unrecognized_res" );  }
namespace in { BooleanOptionKey const ignore_waters( "in:ignore_waters" );  }
namespace in { BooleanOptionKey const add_orbitals( "in:add_orbitals" );  }
namespace in { BooleanOptionKey const show_all_fixes( "in:show_all_fixes" );  }
namespace in { BooleanOptionKey const include_sugars( "in:include_sugars" );  }
namespace in { BooleanOptionKey const include_lipids( "in:include_lipids" );  }
namespace in { BooleanOptionKey const include_surfaces( "in:include_surfaces" );  }
namespace in { BooleanOptionKey const membrane( "in:membrane" );  }
namespace in { BooleanOptionKey const remember_unrecognized_res( "in:remember_unrecognized_res" );  }
namespace in { BooleanOptionKey const remember_unrecognized_water( "in:remember_unrecognized_water" );  }
namespace in { BooleanOptionKey const preserve_crystinfo( "in:preserve_crystinfo" );  }
namespace in { BooleanOptionKey const detect_oops( "in:detect_oops" );  }
namespace in { BooleanOptionKey const detect_disulf( "in:detect_disulf" );  }
namespace in { RealOptionKey const detect_disulf_tolerance( "in:detect_disulf_tolerance" );  }
namespace in { BooleanOptionKey const constraints_from_link_records( "in:constraints_from_link_records" );  }
namespace in { BooleanOptionKey const auto_setup_metals( "in:auto_setup_metals" );  }
namespace in { RealOptionKey const metals_detection_LJ_multiplier( "in:metals_detection_LJ_multiplier" );  }
namespace in { RealOptionKey const metals_distance_constraint_multiplier( "in:metals_distance_constraint_multiplier" );  }
namespace in { RealOptionKey const metals_angle_constraint_multiplier( "in:metals_angle_constraint_multiplier" );  }
namespace in { StringVectorOptionKey const alternate_3_letter_codes( "in:alternate_3_letter_codes" );  }
namespace in { FileOptionKey const fix_disulf( "in:fix_disulf" );  }
namespace in { BooleanOptionKey const missing_density_to_jump( "in:missing_density_to_jump" );  }
namespace in { IntegerVectorOptionKey const target_residues( "in:target_residues" );  }
namespace in { IntegerVectorOptionKey const replonly_residues( "in:replonly_residues" );  }
namespace in { BooleanOptionKey const replonly_loops( "in:replonly_loops" );  }
namespace in { BooleanOptionKey const use_database( "in:use_database" );  }
namespace in { StringVectorOptionKey const select_structures_from_database( "in:select_structures_from_database" );  }
namespace in { namespace dbms { BooleanOptionKey const dbms( "in:dbms" );  } }
namespace in { namespace dbms { StringVectorOptionKey const struct_ids( "in:dbms:struct_ids" );  } }
namespace in { namespace path { PathVectorOptionKey const path( "in:path" );  } }
namespace in { namespace path { PathVectorOptionKey const fragments( "in:path:fragments" );  } }
namespace in { namespace path { PathVectorOptionKey const pdb( "in:path:pdb" );  } }
namespace in { namespace path { PathVectorOptionKey const database( "in:path:database" );  } }
namespace in { namespace file { BooleanOptionKey const file( "in:file" );  } }
namespace in { namespace file { FileVectorOptionKey const s( "in:file:s" );  } }
namespace in { namespace file { FileVectorOptionKey const t( "in:file:t" );  } }
namespace in { namespace file { FileVectorOptionKey const l( "in:file:l" );  } }
namespace in { namespace file { FileVectorOptionKey const list( "in:file:list" );  } }
namespace in { namespace file { FileVectorOptionKey const screening_list( "in:file:screening_list" );  } }
namespace in { namespace file { FileOptionKey const screening_job_file( "in:file:screening_job_file" );  } }
namespace in { namespace file { BooleanOptionKey const shuffle_screening_jobs( "in:file:shuffle_screening_jobs" );  } }
namespace in { namespace file { FileOptionKey const native( "in:file:native" );  } }
namespace in { namespace file { FileOptionKey const torsion_bin_probs( "in:file:torsion_bin_probs" );  } }
namespace in { namespace file { FileOptionKey const PCS_frag_cst( "in:file:PCS_frag_cst" );  } }
namespace in { namespace file { FileOptionKey const talos_phi_psi( "in:file:talos_phi_psi" );  } }
namespace in { namespace file { FileOptionKey const talos_cs( "in:file:talos_cs" );  } }
namespace in { namespace file { FileOptionKey const ambig_talos_cs_A( "in:file:ambig_talos_cs_A" );  } }
namespace in { namespace file { FileOptionKey const ambig_talos_cs_B( "in:file:ambig_talos_cs_B" );  } }
namespace in { namespace file { IntegerVectorOptionKey const native_exclude_res( "in:file:native_exclude_res" );  } }
namespace in { namespace file { StringVectorOptionKey const tags( "in:file:tags" );  } }
namespace in { namespace file { StringVectorOptionKey const user_tags( "in:file:user_tags" );  } }
namespace in { namespace file { FileOptionKey const tagfile( "in:file:tagfile" );  } }
namespace in { namespace file { FileVectorOptionKey const frag_files( "in:file:frag_files" );  } }
namespace in { namespace file { StringVectorOptionKey const remap_pdb_atom_names_for( "in:file:remap_pdb_atom_names_for" );  } }
namespace in { namespace file { FileVectorOptionKey const extra_res( "in:file:extra_res" );  } }
namespace in { namespace file { FileVectorOptionKey const extra_res_fa( "in:file:extra_res_fa" );  } }
namespace in { namespace file { FileVectorOptionKey const extra_res_mol( "in:file:extra_res_mol" );  } }
namespace in { namespace file { StringOptionKey const extra_res_database( "in:file:extra_res_database" );  } }
namespace in { namespace file { StringOptionKey const extra_res_pq_schema( "in:file:extra_res_pq_schema" );  } }
namespace in { namespace file { StringOptionKey const extra_res_database_mode( "in:file:extra_res_database_mode" );  } }
namespace in { namespace file { FileOptionKey const extra_res_database_resname_list( "in:file:extra_res_database_resname_list" );  } }
namespace in { namespace file { FileVectorOptionKey const extra_res_cen( "in:file:extra_res_cen" );  } }
namespace in { namespace file { PathVectorOptionKey const extra_res_path( "in:file:extra_res_path" );  } }
namespace in { namespace file { PathVectorOptionKey const extra_rot_lib_path( "in:file:extra_rot_lib_path" );  } }
namespace in { namespace file { PathVectorOptionKey const extra_res_batch_path( "in:file:extra_res_batch_path" );  } }
namespace in { namespace file { FileVectorOptionKey const extra_patch_fa( "in:file:extra_patch_fa" );  } }
namespace in { namespace file { FileVectorOptionKey const extra_patch_cen( "in:file:extra_patch_cen" );  } }
namespace in { namespace file { StringOptionKey const frag3( "in:file:frag3" );  } }
namespace in { namespace file { StringOptionKey const frag9( "in:file:frag9" );  } }
namespace in { namespace file { StringOptionKey const fragA( "in:file:fragA" );  } }
namespace in { namespace file { StringOptionKey const fragB( "in:file:fragB" );  } }
namespace in { namespace file { StringOptionKey const surface_vectors( "in:file:surface_vectors" );  } }
namespace in { namespace file { StringOptionKey const xyz( "in:file:xyz" );  } }
namespace in { namespace file { BooleanOptionKey const keep_input_scores( "in:file:keep_input_scores" );  } }
namespace in { namespace file { BooleanOptionKey const lazy_silent( "in:file:lazy_silent" );  } }
namespace in { namespace file { FileVectorOptionKey const silent( "in:file:silent" );  } }
namespace in { namespace file { BooleanOptionKey const force_silent_bitflip_on_read( "in:file:force_silent_bitflip_on_read" );  } }
namespace in { namespace file { FileVectorOptionKey const atom_tree_diff( "in:file:atom_tree_diff" );  } }
namespace in { namespace file { StringOptionKey const zip( "in:file:zip" );  } }
namespace in { namespace file { FileVectorOptionKey const boinc_wu_zip( "in:file:boinc_wu_zip" );  } }
namespace in { namespace file { BooleanOptionKey const fullatom( "in:file:fullatom" );  } }
namespace in { namespace file { BooleanOptionKey const centroid_input( "in:file:centroid_input" );  } }
namespace in { namespace file { BooleanOptionKey const centroid( "in:file:centroid" );  } }
namespace in { namespace file { BooleanOptionKey const assign_gasteiger_atom_types( "in:file:assign_gasteiger_atom_types" );  } }
namespace in { namespace file { StringOptionKey const treat_residues_in_these_chains_as_separate_chemical_entities( "in:file:treat_residues_in_these_chains_as_separate_chemical_entities" );  } }
namespace in { namespace file { StringOptionKey const residue_type_set( "in:file:residue_type_set" );  } }
namespace in { namespace file { FileOptionKey const pca( "in:file:pca" );  } }
namespace in { namespace file { RealOptionKey const silent_energy_cut( "in:file:silent_energy_cut" );  } }
namespace in { namespace file { FileVectorOptionKey const silent_list( "in:file:silent_list" );  } }
namespace in { namespace file { BooleanOptionKey const silent_renumber( "in:file:silent_renumber" );  } }
namespace in { namespace file { StringOptionKey const silent_struct_type( "in:file:silent_struct_type" );  } }
namespace in { namespace file { BooleanOptionKey const silent_read_through_errors( "in:file:silent_read_through_errors" );  } }
namespace in { namespace file { StringOptionKey const silent_score_prefix( "in:file:silent_score_prefix" );  } }
namespace in { namespace file { IntegerOptionKey const silent_select_random( "in:file:silent_select_random" );  } }
namespace in { namespace file { IntegerOptionKey const silent_select_range_start( "in:file:silent_select_range_start" );  } }
namespace in { namespace file { IntegerOptionKey const silent_select_range_mul( "in:file:silent_select_range_mul" );  } }
namespace in { namespace file { IntegerOptionKey const silent_select_range_len( "in:file:silent_select_range_len" );  } }
namespace in { namespace file { BooleanOptionKey const skip_failed_simulations( "in:file:skip_failed_simulations" );  } }
namespace in { namespace file { StringVectorOptionKey const silent_scores_wanted( "in:file:silent_scores_wanted" );  } }
namespace in { namespace file { FileVectorOptionKey const fasta( "in:file:fasta" );  } }
namespace in { namespace file { FileVectorOptionKey const pssm( "in:file:pssm" );  } }
namespace in { namespace file { StringVectorOptionKey const seq( "in:file:seq" );  } }
namespace in { namespace file { FileOptionKey const checkpoint( "in:file:checkpoint" );  } }
namespace in { namespace file { FileOptionKey const binary_chk( "in:file:binary_chk" );  } }
namespace in { namespace file { FileVectorOptionKey const alignment( "in:file:alignment" );  } }
namespace in { namespace file { FileVectorOptionKey const alignment2( "in:file:alignment2" );  } }
namespace in { namespace file { FileOptionKey const rama2b_map( "in:file:rama2b_map" );  } }
namespace in { namespace file { FileOptionKey const psipred_ss2( "in:file:psipred_ss2" );  } }
namespace in { namespace file { FileOptionKey const dssp( "in:file:dssp" );  } }
namespace in { namespace file { BooleanOptionKey const fail_on_bad_hbond( "in:file:fail_on_bad_hbond" );  } }
namespace in { namespace file { FileOptionKey const movemap( "in:file:movemap" );  } }
namespace in { namespace file { BooleanOptionKey const repair_sidechains( "in:file:repair_sidechains" );  } }
namespace in { namespace file { BooleanOptionKey const no_binary_dunlib( "in:file:no_binary_dunlib" );  } }
namespace in { namespace file { IntegerOptionKey const extended_pose( "in:file:extended_pose" );  } }
namespace in { namespace file { FileVectorOptionKey const template_pdb( "in:file:template_pdb" );  } }
namespace in { namespace file { FileOptionKey const template_silent( "in:file:template_silent" );  } }
namespace in { namespace file { FileVectorOptionKey const rdc( "in:file:rdc" );  } }
namespace in { namespace file { FileVectorOptionKey const csa( "in:file:csa" );  } }
namespace in { namespace file { FileVectorOptionKey const dc( "in:file:dc" );  } }
namespace in { namespace file { FileVectorOptionKey const burial( "in:file:burial" );  } }
namespace in { namespace file { FileVectorOptionKey const vall( "in:file:vall" );  } }
namespace in { namespace file { BooleanOptionKey const rescore( "in:file:rescore" );  } }
namespace in { namespace file { StringOptionKey const spanfile( "in:file:spanfile" );  } }
namespace in { namespace file { StringOptionKey const lipofile( "in:file:lipofile" );  } }
namespace in { namespace file { StringOptionKey const HDX( "in:file:HDX" );  } }
namespace in { namespace file { RealOptionKey const d2h_sa_reweight( "in:file:d2h_sa_reweight" );  } }
namespace in { namespace file { FileOptionKey const sucker_params( "in:file:sucker_params" );  } }
namespace in { namespace file { FileOptionKey const fold_tree( "in:file:fold_tree" );  } }
namespace in { namespace file { BooleanOptionKey const obey_ENDMDL( "in:file:obey_ENDMDL" );  } }
namespace in { namespace file { BooleanOptionKey const new_chain_order( "in:file:new_chain_order" );  } }
namespace in { namespace file { FileOptionKey const ddg_predictions_file( "in:file:ddg_predictions_file" );  } }
namespace in { namespace file { ResidueChainVectorOptionKey const input_res( "in:file:input_res" );  } }
namespace in { namespace file { IntegerVectorOptionKey const minimize_res( "in:file:minimize_res" );  } }
namespace in { namespace file { StringOptionKey const md_schfile( "in:file:md_schfile" );  } }
namespace in { namespace file { BooleanOptionKey const read_pdb_link_records( "in:file:read_pdb_link_records" );  } }
namespace in { namespace file { FileOptionKey const native_contacts( "in:file:native_contacts" );  } }
namespace in { namespace rdf { BooleanOptionKey const rdf( "in:rdf" );  } }
namespace in { namespace rdf { BooleanOptionKey const sep_bb_ss( "in:rdf:sep_bb_ss" );  } }
namespace inout { BooleanOptionKey const inout( "inout" );  }
namespace inout { BooleanOptionKey const fold_tree_io( "inout:fold_tree_io" );  }
namespace inout { BooleanOptionKey const dump_connect_info( "inout:dump_connect_info" );  }
namespace inout { RealOptionKey const connect_info_cutoff( "inout:connect_info_cutoff" );  }
namespace inout { namespace dbms { BooleanOptionKey const dbms( "inout:dbms" );  } }
namespace inout { namespace dbms { StringOptionKey const mode( "inout:dbms:mode" );  } }
namespace inout { namespace dbms { StringOptionKey const database_name( "inout:dbms:database_name" );  } }
namespace inout { namespace dbms { StringOptionKey const pq_schema( "inout:dbms:pq_schema" );  } }
namespace inout { namespace dbms { StringOptionKey const host( "inout:dbms:host" );  } }
namespace inout { namespace dbms { StringOptionKey const user( "inout:dbms:user" );  } }
namespace inout { namespace dbms { StringOptionKey const password( "inout:dbms:password" );  } }
namespace inout { namespace dbms { IntegerOptionKey const port( "inout:dbms:port" );  } }
namespace inout { namespace dbms { BooleanOptionKey const readonly( "inout:dbms:readonly" );  } }
namespace inout { namespace dbms { BooleanOptionKey const separate_db_per_mpi_process( "inout:dbms:separate_db_per_mpi_process" );  } }
namespace inout { namespace dbms { IntegerOptionKey const database_partition( "inout:dbms:database_partition" );  } }
namespace inout { namespace dbms { BooleanOptionKey const use_compact_residue_schema( "inout:dbms:use_compact_residue_schema" );  } }
namespace inout { namespace dbms { BooleanOptionKey const retry_failed_reads( "inout:dbms:retry_failed_reads" );  } }
namespace inout { namespace dbms { PathOptionKey const path( "inout:dbms:path" );  } }
namespace out { BooleanOptionKey const out( "out" );  }
namespace out { BooleanOptionKey const overwrite( "out:overwrite" );  }
namespace out { IntegerOptionKey const nstruct( "out:nstruct" );  }
namespace out { IntegerOptionKey const shuffle_nstruct( "out:shuffle_nstruct" );  }
namespace out { StringOptionKey const prefix( "out:prefix" );  }
namespace out { StringOptionKey const suffix( "out:suffix" );  }
namespace out { BooleanOptionKey const no_nstruct_label( "out:no_nstruct_label" );  }
namespace out { BooleanOptionKey const pdb_gz( "out:pdb_gz" );  }
namespace out { BooleanOptionKey const pdb( "out:pdb" );  }
namespace out { BooleanOptionKey const silent_gz( "out:silent_gz" );  }
namespace out { BooleanOptionKey const membrane_pdb( "out:membrane_pdb" );  }
namespace out { RealOptionKey const membrane_pdb_thickness( "out:membrane_pdb_thickness" );  }
namespace out { BooleanOptionKey const use_database( "out:use_database" );  }
namespace out { IntegerOptionKey const database_protocol_id( "out:database_protocol_id" );  }
namespace out { StringVectorOptionKey const database_filter( "out:database_filter" );  }
namespace out { IntegerVectorOptionKey const resume_batch( "out:resume_batch" );  }
namespace out { BooleanOptionKey const nooutput( "out:nooutput" );  }
namespace out { BooleanOptionKey const output( "out:output" );  }
namespace out { RealOptionKey const scorecut( "out:scorecut" );  }
namespace out { BooleanOptionKey const show_accessed_options( "out:show_accessed_options" );  }
namespace out { BooleanOptionKey const show_unused_options( "out:show_unused_options" );  }
namespace out { FileOptionKey const sf( "out:sf" );  }
namespace out { StringVectorOptionKey const mute( "out:mute" );  }
namespace out { StringVectorOptionKey const unmute( "out:unmute" );  }
namespace out { IntegerOptionKey const level( "out:level" );  }
namespace out { StringVectorOptionKey const levels( "out:levels" );  }
namespace out { IntegerOptionKey const std_IO_exit_error_code( "out:std_IO_exit_error_code" );  }
namespace out { BooleanOptionKey const chname( "out:chname" );  }
namespace out { BooleanOptionKey const chtimestamp( "out:chtimestamp" );  }
namespace out { BooleanOptionKey const dry_run( "out:dry_run" );  }
namespace out { StringOptionKey const mpi_tracer_to_file( "out:mpi_tracer_to_file" );  }
namespace out { StringOptionKey const user_tag( "out:user_tag" );  }
namespace out { StringOptionKey const output_tag( "out:output_tag" );  }
namespace out { namespace file { BooleanOptionKey const file( "out:file" );  } }
namespace out { namespace file { StringOptionKey const o( "out:file:o" );  } }
namespace out { namespace file { FileOptionKey const design_contrast( "out:file:design_contrast" );  } }
namespace out { namespace file { StringOptionKey const residue_type_set( "out:file:residue_type_set" );  } }
namespace out { namespace file { StringOptionKey const atom_tree_diff( "out:file:atom_tree_diff" );  } }
namespace out { namespace file { IntegerOptionKey const atom_tree_diff_bb( "out:file:atom_tree_diff_bb" );  } }
namespace out { namespace file { IntegerOptionKey const atom_tree_diff_sc( "out:file:atom_tree_diff_sc" );  } }
namespace out { namespace file { IntegerOptionKey const atom_tree_diff_bl( "out:file:atom_tree_diff_bl" );  } }
namespace out { namespace file { StringOptionKey const alignment( "out:file:alignment" );  } }
namespace out { namespace file { StringOptionKey const score_only( "out:file:score_only" );  } }
namespace out { namespace file { StringOptionKey const scorefile( "out:file:scorefile" );  } }
namespace out { namespace file { StringOptionKey const scorefile_format( "out:file:scorefile_format" );  } }
namespace out { namespace file { StringOptionKey const silent( "out:file:silent" );  } }
namespace out { namespace file { StringOptionKey const silent_struct_type( "out:file:silent_struct_type" );  } }
namespace out { namespace file { BooleanOptionKey const silent_print_all_score_headers( "out:file:silent_print_all_score_headers" );  } }
namespace out { namespace file { BooleanOptionKey const raw( "out:file:raw" );  } }
namespace out { namespace file { BooleanOptionKey const weight_silent_scores( "out:file:weight_silent_scores" );  } }
namespace out { namespace file { BooleanOptionKey const silent_preserve_H( "out:file:silent_preserve_H" );  } }
namespace out { namespace file { BooleanOptionKey const fullatom( "out:file:fullatom" );  } }
namespace out { namespace file { BooleanOptionKey const suppress_zero_occ_pdb_output( "out:file:suppress_zero_occ_pdb_output" );  } }
namespace out { namespace file { BooleanOptionKey const output_virtual( "out:file:output_virtual" );  } }
namespace out { namespace file { BooleanOptionKey const output_virtual_zero_occ( "out:file:output_virtual_zero_occ" );  } }
namespace out { namespace file { BooleanOptionKey const no_output_cen( "out:file:no_output_cen" );  } }
namespace out { namespace file { BooleanOptionKey const output_orbitals( "out:file:output_orbitals" );  } }
namespace out { namespace file { BooleanOptionKey const no_scores_in_pdb( "out:file:no_scores_in_pdb" );  } }
