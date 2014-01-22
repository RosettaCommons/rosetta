namespace in { BooleanOptionKey const in( "in" );  }
namespace in { StringOptionKey const Ntermini( "in:Ntermini" );  }
namespace in { StringOptionKey const Ctermini( "in:Ctermini" );  }
namespace in { BooleanOptionKey const use_truncated_termini( "in:use_truncated_termini" );  }
namespace in { BooleanOptionKey const ignore_unrecognized_res( "in:ignore_unrecognized_res" );  }
namespace in { BooleanOptionKey const ignore_waters( "in:ignore_waters" );  }
namespace in { BooleanOptionKey const add_orbitals( "in:add_orbitals" );  }
namespace in { BooleanOptionKey const show_all_fixes( "in:show_all_fixes" );  }
namespace in { BooleanOptionKey const include_sugars( "in:include_sugars" );  }
namespace in { BooleanOptionKey const include_surfaces( "in:include_surfaces" );  }
namespace in { BooleanOptionKey const enable_branching( "in:enable_branching" );  }
namespace in { BooleanOptionKey const remember_unrecognized_res( "in:remember_unrecognized_res" );  }
namespace in { BooleanOptionKey const remember_unrecognized_water( "in:remember_unrecognized_water" );  }
namespace in { BooleanOptionKey const preserve_crystinfo( "in:preserve_crystinfo" );  }
namespace in { BooleanOptionKey const detect_oops( "in:detect_oops" );  }
namespace in { BooleanOptionKey const detect_disulf( "in:detect_disulf" );  }
namespace in { RealOptionKey const detect_disulf_tolerance( "in:detect_disulf_tolerance" );  }
namespace in { FileOptionKey const fix_disulf( "in:fix_disulf" );  }
namespace in { BooleanOptionKey const missing_density_to_jump( "in:missing_density_to_jump" );  }
namespace in { IntegerVectorOptionKey const target_residues( "in:target_residues" );  }
namespace in { IntegerVectorOptionKey const replonly_residues( "in:replonly_residues" );  }
namespace in { BooleanOptionKey const replonly_loops( "in:replonly_loops" );  }
namespace in { BooleanOptionKey const use_database( "in:use_database" );  }
namespace in { namespace dbms { BooleanOptionKey const dbms( "in:dbms" );  } }
namespace in { namespace dbms { StringVectorOptionKey const struct_ids( "in:dbms:struct_ids" );  } }
namespace in { IntegerOptionKey const database_protocol( "in:database_protocol" );  }
namespace in { StringVectorOptionKey const select_structures_from_database( "in:select_structures_from_database" );  }
namespace in { namespace path { PathVectorOptionKey const path( "in:path" );  } }
namespace in { namespace path { PathVectorOptionKey const fragments( "in:path:fragments" );  } }
namespace in { namespace path { PathVectorOptionKey const pdb( "in:path:pdb" );  } }
namespace in { namespace path { PathVectorOptionKey const database( "in:path:database" );  } }
namespace in { namespace file { BooleanOptionKey const file( "in:file" );  } }
namespace in { namespace file { FileVectorOptionKey const s( "in:file:s" );  } }
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
namespace in { namespace file { IntegerVectorOptionKey const frag_sizes( "in:file:frag_sizes" );  } }
namespace in { namespace file { FileVectorOptionKey const extra_res( "in:file:extra_res" );  } }
namespace in { namespace file { FileVectorOptionKey const extra_res_fa( "in:file:extra_res_fa" );  } }
namespace in { namespace file { FileVectorOptionKey const extra_res_mol( "in:file:extra_res_mol" );  } }
namespace in { namespace file { StringOptionKey const extra_res_database( "in:file:extra_res_database" );  } }
namespace in { namespace file { StringOptionKey const extra_res_pq_schema( "in:file:extra_res_pq_schema" );  } }
namespace in { namespace file { StringOptionKey const extra_res_database_mode( "in:file:extra_res_database_mode" );  } }
namespace in { namespace file { FileOptionKey const extra_res_database_resname_list( "in:file:extra_res_database_resname_list" );  } }
namespace in { namespace file { FileVectorOptionKey const extra_res_cen( "in:file:extra_res_cen" );  } }
namespace in { namespace file { PathVectorOptionKey const extra_res_path( "in:file:extra_res_path" );  } }
namespace in { namespace file { PathVectorOptionKey const extra_res_batch_path( "in:file:extra_res_batch_path" );  } }
namespace in { namespace file { FileVectorOptionKey const extra_patch_fa( "in:file:extra_patch_fa" );  } }
namespace in { namespace file { FileVectorOptionKey const extra_patch_cen( "in:file:extra_patch_cen" );  } }
namespace in { namespace file { StringOptionKey const frag3( "in:file:frag3" );  } }
namespace in { namespace file { StringOptionKey const frag9( "in:file:frag9" );  } }
namespace in { namespace file { StringOptionKey const fragA( "in:file:fragA" );  } }
namespace in { namespace file { StringOptionKey const fragB( "in:file:fragB" );  } }
namespace in { namespace file { StringOptionKey const surface_vectors( "in:file:surface_vectors" );  } }
namespace in { namespace file { StringOptionKey const xyz( "in:file:xyz" );  } }
namespace in { namespace file { IntegerOptionKey const fragA_size( "in:file:fragA_size" );  } }
namespace in { namespace file { BooleanOptionKey const keep_input_scores( "in:file:keep_input_scores" );  } }
namespace in { namespace file { BooleanOptionKey const lazy_silent( "in:file:lazy_silent" );  } }
namespace in { namespace file { FileVectorOptionKey const silent( "in:file:silent" );  } }
namespace in { namespace file { FileVectorOptionKey const atom_tree_diff( "in:file:atom_tree_diff" );  } }
namespace in { namespace file { StringOptionKey const zip( "in:file:zip" );  } }
namespace in { namespace file { FileVectorOptionKey const boinc_wu_zip( "in:file:boinc_wu_zip" );  } }
namespace in { namespace file { BooleanOptionKey const fullatom( "in:file:fullatom" );  } }
namespace in { namespace file { BooleanOptionKey const centroid_input( "in:file:centroid_input" );  } }
namespace in { namespace file { BooleanOptionKey const centroid( "in:file:centroid" );  } }
namespace in { namespace file { StringOptionKey const treat_residues_in_these_chains_as_separate_chemical_entities( "in:file:treat_residues_in_these_chains_as_separate_chemical_entities" );  } }
namespace in { namespace file { StringOptionKey const residue_type_set( "in:file:residue_type_set" );  } }
namespace in { namespace file { FileOptionKey const pca( "in:file:pca" );  } }
namespace in { namespace file { RealOptionKey const silent_energy_cut( "in:file:silent_energy_cut" );  } }
namespace in { namespace file { FileVectorOptionKey const silent_list( "in:file:silent_list" );  } }
namespace in { namespace file { BooleanOptionKey const silent_renumber( "in:file:silent_renumber" );  } }
namespace in { namespace file { BooleanOptionKey const silent_optH( "in:file:silent_optH" );  } }
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
namespace in { namespace file { IntegerVectorOptionKey const input_res( "in:file:input_res" );  } }
namespace in { namespace file { IntegerVectorOptionKey const minimize_res( "in:file:minimize_res" );  } }
namespace in { namespace file { StringOptionKey const md_schfile( "in:file:md_schfile" );  } }
namespace in { namespace file { BooleanOptionKey const read_pdb_link_records( "in:file:read_pdb_link_records" );  } }
namespace in { namespace file { FileOptionKey const native_contacts( "in:file:native_contacts" );  } }
namespace in { namespace rdf { BooleanOptionKey const rdf( "in:rdf" );  } }
namespace in { namespace rdf { BooleanOptionKey const sep_bb_ss( "in:rdf:sep_bb_ss" );  } }
namespace inout { BooleanOptionKey const inout( "inout" );  }
namespace inout { BooleanOptionKey const fold_tree_io( "inout:fold_tree_io" );  }
namespace inout { BooleanOptionKey const dump_connect_info( "inout:dump_connect_info" );  }
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
namespace out { BooleanOptionKey const out( "out" );  }
namespace out { BooleanOptionKey const overwrite( "out:overwrite" );  }
namespace out { IntegerOptionKey const nstruct( "out:nstruct" );  }
namespace out { IntegerOptionKey const shuffle_nstruct( "out:shuffle_nstruct" );  }
namespace out { StringOptionKey const prefix( "out:prefix" );  }
namespace out { StringOptionKey const suffix( "out:suffix" );  }
namespace out { StringOptionKey const force_output_name( "out:force_output_name" );  }
namespace out { BooleanOptionKey const no_nstruct_label( "out:no_nstruct_label" );  }
namespace out { BooleanOptionKey const pdb_gz( "out:pdb_gz" );  }
namespace out { BooleanOptionKey const pdb( "out:pdb" );  }
namespace out { BooleanOptionKey const silent_gz( "out:silent_gz" );  }
namespace out { BooleanOptionKey const use_database( "out:use_database" );  }
namespace out { IntegerOptionKey const database_protocol_id( "out:database_protocol_id" );  }
namespace out { StringVectorOptionKey const database_filter( "out:database_filter" );  }
namespace out { IntegerVectorOptionKey const resume_batch( "out:resume_batch" );  }
namespace out { BooleanOptionKey const nooutput( "out:nooutput" );  }
namespace out { BooleanOptionKey const output( "out:output" );  }
namespace out { RealOptionKey const scorecut( "out:scorecut" );  }
namespace out { BooleanOptionKey const show_accessed_options( "out:show_accessed_options" );  }
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
namespace out { namespace file { StringOptionKey const silent( "out:file:silent" );  } }
namespace out { namespace file { StringOptionKey const silent_struct_type( "out:file:silent_struct_type" );  } }
namespace out { namespace file { BooleanOptionKey const silent_print_all_score_headers( "out:file:silent_print_all_score_headers" );  } }
namespace out { namespace file { BooleanOptionKey const silent_decoytime( "out:file:silent_decoytime" );  } }
namespace out { namespace file { IntegerOptionKey const silent_comment_bound( "out:file:silent_comment_bound" );  } }
namespace out { namespace file { BooleanOptionKey const raw( "out:file:raw" );  } }
namespace out { namespace file { BooleanOptionKey const weight_silent_scores( "out:file:weight_silent_scores" );  } }
namespace out { namespace file { BooleanOptionKey const silent_preserve_H( "out:file:silent_preserve_H" );  } }
namespace out { namespace file { BooleanOptionKey const fullatom( "out:file:fullatom" );  } }
namespace out { namespace file { BooleanOptionKey const suppress_zero_occ_pdb_output( "out:file:suppress_zero_occ_pdb_output" );  } }
namespace out { namespace file { BooleanOptionKey const output_virtual( "out:file:output_virtual" );  } }
namespace out { namespace file { BooleanOptionKey const no_output_cen( "out:file:no_output_cen" );  } }
namespace out { namespace file { BooleanOptionKey const output_orbitals( "out:file:output_orbitals" );  } }
namespace out { namespace file { BooleanOptionKey const renumber_pdb( "out:file:renumber_pdb" );  } }
namespace out { namespace file { BooleanOptionKey const pdb_parents( "out:file:pdb_parents" );  } }
namespace out { namespace file { BooleanOptionKey const per_chain_renumbering( "out:file:per_chain_renumbering" );  } }
namespace out { namespace file { BooleanOptionKey const output_torsions( "out:file:output_torsions" );  } }
namespace out { namespace file { BooleanOptionKey const pdb_comments( "out:file:pdb_comments" );  } }
namespace out { namespace file { BooleanOptionKey const force_nonideal_structure( "out:file:force_nonideal_structure" );  } }
namespace out { namespace file { BooleanOptionKey const write_pdb_link_records( "out:file:write_pdb_link_records" );  } }
namespace out { namespace file { BooleanOptionKey const dont_rewrite_dunbrack_database( "out:file:dont_rewrite_dunbrack_database" );  } }
namespace out { namespace file { StringOptionKey const frag_prefix( "out:file:frag_prefix" );  } }
namespace out { namespace path { PathOptionKey const all( "out:path:all" );  } }
namespace out { namespace path { PathOptionKey const path( "out:path" );  } }
namespace out { namespace path { PathOptionKey const pdb( "out:path:pdb" );  } }
namespace out { namespace path { PathOptionKey const score( "out:path:score" );  } }
namespace out { namespace path { PathOptionKey const movie( "out:path:movie" );  } }
namespace out { namespace path { PathOptionKey const scratch( "out:path:scratch" );  } }
namespace out { namespace path { BooleanOptionKey const mpi_rank_dir( "out:path:mpi_rank_dir" );  } }
namespace rigid { BooleanOptionKey const rigid( "rigid" );  }
namespace rigid { RealOptionKey const chainbreak_bias( "rigid:chainbreak_bias" );  }
namespace rigid { BooleanOptionKey const close_loops( "rigid:close_loops" );  }
namespace rigid { IntegerOptionKey const fragment_cycles( "rigid:fragment_cycles" );  }
namespace rigid { BooleanOptionKey const log_accepted_moves( "rigid:log_accepted_moves" );  }
namespace rigid { RealOptionKey const max_ca_ca_dist( "rigid:max_ca_ca_dist" );  }
namespace rigid { IntegerOptionKey const medium_range_seqsep( "rigid:medium_range_seqsep" );  }
namespace rigid { FileOptionKey const patch( "rigid:patch" );  }
namespace rigid { IntegerOptionKey const residues_backbone_move( "rigid:residues_backbone_move" );  }
namespace rigid { RealOptionKey const rotation( "rigid:rotation" );  }
namespace rigid { FileOptionKey const sampling_prob( "rigid:sampling_prob" );  }
namespace rigid { StringOptionKey const score( "rigid:score" );  }
namespace rigid { IntegerOptionKey const sequence_separation( "rigid:sequence_separation" );  }
namespace rigid { IntegerOptionKey const short_range_seqsep( "rigid:short_range_seqsep" );  }
namespace rigid { IntegerOptionKey const small_cycles( "rigid:small_cycles" );  }
namespace rigid { IntegerOptionKey const stages( "rigid:stages" );  }
namespace rigid { RealOptionKey const temperature( "rigid:temperature" );  }
namespace rigid { RealOptionKey const translation( "rigid:translation" );  }
namespace MM { BooleanOptionKey const MM( "MM" );  }
namespace MM { BooleanOptionKey const ignore_missing_bondangle_params( "MM:ignore_missing_bondangle_params" );  }
namespace qsar { BooleanOptionKey const qsar( "qsar" );  }
namespace qsar { StringOptionKey const weights( "qsar:weights" );  }
namespace qsar { StringOptionKey const grid_dir( "qsar:grid_dir" );  }
namespace qsar { IntegerOptionKey const max_grid_cache_size( "qsar:max_grid_cache_size" );  }
namespace residues { BooleanOptionKey const residues( "residues" );  }
namespace residues { StringVectorOptionKey const patch_selectors( "residues:patch_selectors" );  }
namespace PCS { BooleanOptionKey const PCS( "PCS" );  }
namespace PCS { FileOptionKey const write_extra( "PCS:write_extra" );  }
namespace PCS { IntegerOptionKey const normalization_id( "PCS:normalization_id" );  }
namespace pocket_grid { BooleanOptionKey const pocket_grid( "pocket_grid" );  }
namespace pocket_grid { RealOptionKey const pocket_grid_size( "pocket_grid:pocket_grid_size" );  }
namespace pocket_grid { RealOptionKey const pocket_grid_size_x( "pocket_grid:pocket_grid_size_x" );  }
namespace pocket_grid { RealOptionKey const pocket_grid_size_y( "pocket_grid:pocket_grid_size_y" );  }
namespace pocket_grid { RealOptionKey const pocket_grid_size_z( "pocket_grid:pocket_grid_size_z" );  }
namespace pocket_grid { RealOptionKey const pocket_grid_spacing( "pocket_grid:pocket_grid_spacing" );  }
namespace pocket_grid { RealOptionKey const pocket_max_spacing( "pocket_grid:pocket_max_spacing" );  }
namespace pocket_grid { RealOptionKey const pocket_min_size( "pocket_grid:pocket_min_size" );  }
namespace pocket_grid { RealOptionKey const pocket_max_size( "pocket_grid:pocket_max_size" );  }
namespace pocket_grid { RealOptionKey const pocket_probe_radius( "pocket_grid:pocket_probe_radius" );  }
namespace pocket_grid { StringOptionKey const central_relax_pdb_num( "pocket_grid:central_relax_pdb_num" );  }
namespace pocket_grid { IntegerOptionKey const pocket_ntrials( "pocket_grid:pocket_ntrials" );  }
namespace pocket_grid { IntegerOptionKey const pocket_num_angles( "pocket_grid:pocket_num_angles" );  }
namespace pocket_grid { BooleanOptionKey const pocket_side( "pocket_grid:pocket_side" );  }
namespace pocket_grid { BooleanOptionKey const pocket_dump_pdbs( "pocket_grid:pocket_dump_pdbs" );  }
namespace pocket_grid { BooleanOptionKey const pocket_dump_exemplars( "pocket_grid:pocket_dump_exemplars" );  }
namespace pocket_grid { BooleanOptionKey const pocket_filter_by_exemplar( "pocket_grid:pocket_filter_by_exemplar" );  }
namespace pocket_grid { BooleanOptionKey const pocket_dump_rama( "pocket_grid:pocket_dump_rama" );  }
namespace pocket_grid { BooleanOptionKey const pocket_restrict_size( "pocket_grid:pocket_restrict_size" );  }
namespace pocket_grid { BooleanOptionKey const pocket_ignore_buried( "pocket_grid:pocket_ignore_buried" );  }
namespace pocket_grid { BooleanOptionKey const pocket_only_buried( "pocket_grid:pocket_only_buried" );  }
namespace pocket_grid { BooleanOptionKey const pocket_psp( "pocket_grid:pocket_psp" );  }
namespace pocket_grid { BooleanOptionKey const pocket_sps( "pocket_grid:pocket_sps" );  }
namespace pocket_grid { BooleanOptionKey const pocket_search13( "pocket_grid:pocket_search13" );  }
namespace pocket_grid { RealOptionKey const pocket_surface_score( "pocket_grid:pocket_surface_score" );  }
namespace pocket_grid { RealOptionKey const pocket_surface_dist( "pocket_grid:pocket_surface_dist" );  }
namespace pocket_grid { RealOptionKey const pocket_buried_score( "pocket_grid:pocket_buried_score" );  }
namespace pocket_grid { RealOptionKey const pocket_buried_dist( "pocket_grid:pocket_buried_dist" );  }
namespace pocket_grid { RealOptionKey const pocket_exemplar_vdw_pen( "pocket_grid:pocket_exemplar_vdw_pen" );  }
namespace pocket_grid { BooleanOptionKey const pocket_debug_output( "pocket_grid:pocket_debug_output" );  }
namespace pocket_grid { BooleanOptionKey const print_grid( "pocket_grid:print_grid" );  }
namespace pocket_grid { BooleanOptionKey const extend_eggshell( "pocket_grid:extend_eggshell" );  }
namespace pocket_grid { RealOptionKey const extend_eggshell_dist( "pocket_grid:extend_eggshell_dist" );  }
namespace pocket_grid { RealOptionKey const extra_eggshell_dist( "pocket_grid:extra_eggshell_dist" );  }
namespace pocket_grid { RealOptionKey const eggshell_dist( "pocket_grid:eggshell_dist" );  }
namespace pocket_grid { BooleanOptionKey const reduce_rays( "pocket_grid:reduce_rays" );  }
namespace pocket_grid { BooleanOptionKey const pocket_static_grid( "pocket_grid:pocket_static_grid" );  }
namespace fingerprint { BooleanOptionKey const fingerprint( "fingerprint" );  }
namespace fingerprint { BooleanOptionKey const print_eggshell( "fingerprint:print_eggshell" );  }
namespace fingerprint { RealOptionKey const atom_radius_scale( "fingerprint:atom_radius_scale" );  }
namespace fingerprint { RealOptionKey const atom_radius_buffer( "fingerprint:atom_radius_buffer" );  }
namespace fingerprint { RealOptionKey const packing_weight( "fingerprint:packing_weight" );  }
namespace fingerprint { RealOptionKey const dist_cut_off( "fingerprint:dist_cut_off" );  }
namespace fingerprint { BooleanOptionKey const include_hydrogens( "fingerprint:include_hydrogens" );  }
namespace fingerprint { BooleanOptionKey const use_DARC_gpu( "fingerprint:use_DARC_gpu" );  }
namespace fingerprint { BooleanOptionKey const square_score( "fingerprint:square_score" );  }
namespace fingerprint { IntegerOptionKey const set_origin( "fingerprint:set_origin" );  }
namespace fingerprint { IntegerOptionKey const origin_res_num( "fingerprint:origin_res_num" );  }
namespace contactMap { BooleanOptionKey const contactMap( "contactMap" );  }
namespace contactMap { StringOptionKey const prefix( "contactMap:prefix" );  }
namespace contactMap { RealOptionKey const distance_cutoff( "contactMap:distance_cutoff" );  }
namespace contactMap { RealOptionKey const energy_cutoff( "contactMap:energy_cutoff" );  }
namespace contactMap { StringOptionKey const region_def( "contactMap:region_def" );  }
namespace contactMap { BooleanOptionKey const row_format( "contactMap:row_format" );  }
namespace contactMap { BooleanOptionKey const distance_matrix( "contactMap:distance_matrix" );  }
namespace docking { BooleanOptionKey const kick_relax( "docking:kick_relax" );  }
namespace docking { BooleanOptionKey const docking( "docking" );  }
namespace docking { BooleanOptionKey const view( "docking:view" );  }
namespace docking { BooleanOptionKey const no_filters( "docking:no_filters" );  }
namespace docking { StringVectorOptionKey const design_chains( "docking:design_chains" );  }
namespace docking { FileOptionKey const recover_sidechains( "docking:recover_sidechains" );  }
namespace docking { StringOptionKey const partners( "docking:partners" );  }
namespace docking { BooleanOptionKey const docking_local_refine( "docking:docking_local_refine" );  }
namespace docking { BooleanOptionKey const low_res_protocol_only( "docking:low_res_protocol_only" );  }
namespace docking { BooleanOptionKey const randomize1( "docking:randomize1" );  }
namespace docking { BooleanOptionKey const randomize2( "docking:randomize2" );  }
namespace docking { BooleanOptionKey const use_ellipsoidal_randomization( "docking:use_ellipsoidal_randomization" );  }
namespace docking { BooleanOptionKey const spin( "docking:spin" );  }
namespace docking { RealVectorOptionKey const dock_pert( "docking:dock_pert" );  }
namespace docking { RealOptionKey const uniform_trans( "docking:uniform_trans" );  }
namespace docking { BooleanOptionKey const center_at_interface( "docking:center_at_interface" );  }
namespace docking { IntegerOptionKey const dock_mcm_first_cycles( "docking:dock_mcm_first_cycles" );  }
namespace docking { IntegerOptionKey const dock_mcm_second_cycles( "docking:dock_mcm_second_cycles" );  }
namespace docking { IntegerOptionKey const docking_centroid_outer_cycles( "docking:docking_centroid_outer_cycles" );  }
namespace docking { IntegerOptionKey const docking_centroid_inner_cycles( "docking:docking_centroid_inner_cycles" );  }
namespace docking { BooleanOptionKey const dock_min( "docking:dock_min" );  }
namespace docking { StringOptionKey const flexible_bb_docking( "docking:flexible_bb_docking" );  }
namespace docking { RealOptionKey const flexible_bb_docking_interface_dist( "docking:flexible_bb_docking_interface_dist" );  }
namespace docking { StringOptionKey const ensemble1( "docking:ensemble1" );  }
namespace docking { StringOptionKey const ensemble2( "docking:ensemble2" );  }
namespace docking { RealOptionKey const dock_mcm_trans_magnitude( "docking:dock_mcm_trans_magnitude" );  }
namespace docking { RealOptionKey const dock_mcm_rot_magnitude( "docking:dock_mcm_rot_magnitude" );  }
namespace docking { RealOptionKey const minimization_threshold( "docking:minimization_threshold" );  }
namespace docking { RealOptionKey const temperature( "docking:temperature" );  }
namespace docking { IntegerOptionKey const repack_period( "docking:repack_period" );  }
namespace docking { BooleanOptionKey const extra_rottrial( "docking:extra_rottrial" );  }
namespace docking { BooleanOptionKey const dock_rtmin( "docking:dock_rtmin" );  }
namespace docking { BooleanOptionKey const sc_min( "docking:sc_min" );  }
namespace docking { BooleanOptionKey const norepack1( "docking:norepack1" );  }
namespace docking { BooleanOptionKey const norepack2( "docking:norepack2" );  }
namespace docking { IntegerVectorOptionKey const bb_min_res( "docking:bb_min_res" );  }
namespace docking { IntegerVectorOptionKey const sc_min_res( "docking:sc_min_res" );  }
namespace docking { BooleanOptionKey const dock_ppk( "docking:dock_ppk" );  }
namespace docking { IntegerOptionKey const max_repeats( "docking:max_repeats" );  }
namespace docking { RealVectorOptionKey const dock_lowres_filter( "docking:dock_lowres_filter" );  }
namespace docking { IntegerVectorOptionKey const multibody( "docking:multibody" );  }
namespace docking { BooleanOptionKey const ignore_default_docking_task( "docking:ignore_default_docking_task" );  }
namespace docking { StringOptionKey const low_patch( "docking:low_patch" );  }
namespace docking { StringOptionKey const high_patch( "docking:high_patch" );  }
namespace docking { StringOptionKey const high_min_patch( "docking:high_min_patch" );  }
namespace docking { StringOptionKey const pack_patch( "docking:pack_patch" );  }
namespace docking { BooleanOptionKey const use_legacy_protocol( "docking:use_legacy_protocol" );  }
namespace docking { RealOptionKey const docklowres_trans_magnitude( "docking:docklowres_trans_magnitude" );  }
namespace docking { RealOptionKey const docklowres_rot_magnitude( "docking:docklowres_rot_magnitude" );  }
namespace docking { namespace ligand { BooleanOptionKey const ligand( "docking:ligand" );  } }
namespace docking { namespace ligand { StringOptionKey const protocol( "docking:ligand:protocol" );  } }
namespace docking { namespace ligand { BooleanOptionKey const soft_rep( "docking:ligand:soft_rep" );  } }
namespace docking { namespace ligand { BooleanOptionKey const tweak_sxfn( "docking:ligand:tweak_sxfn" );  } }
namespace docking { namespace ligand { BooleanOptionKey const old_estat( "docking:ligand:old_estat" );  } }
namespace docking { namespace ligand { BooleanOptionKey const random_conformer( "docking:ligand:random_conformer" );  } }
namespace docking { namespace ligand { IntegerOptionKey const improve_orientation( "docking:ligand:improve_orientation" );  } }
namespace docking { namespace ligand { BooleanOptionKey const mutate_same_name3( "docking:ligand:mutate_same_name3" );  } }
namespace docking { namespace ligand { RealOptionKey const subset_to_keep( "docking:ligand:subset_to_keep" );  } }
namespace docking { namespace ligand { RealOptionKey const min_rms( "docking:ligand:min_rms" );  } }
namespace docking { namespace ligand { IntegerOptionKey const max_poses( "docking:ligand:max_poses" );  } }
namespace docking { namespace ligand { BooleanOptionKey const minimize_ligand( "docking:ligand:minimize_ligand" );  } }
namespace docking { namespace ligand { RealOptionKey const harmonic_torsions( "docking:ligand:harmonic_torsions" );  } }
namespace docking { namespace ligand { BooleanOptionKey const use_ambig_constraints( "docking:ligand:use_ambig_constraints" );  } }
namespace docking { namespace ligand { IntegerOptionKey const shear_moves( "docking:ligand:shear_moves" );  } }
namespace docking { namespace ligand { BooleanOptionKey const minimize_backbone( "docking:ligand:minimize_backbone" );  } }
namespace docking { namespace ligand { RealOptionKey const harmonic_Calphas( "docking:ligand:harmonic_Calphas" );  } }
namespace docking { namespace ligand { RealOptionKey const tether_ligand( "docking:ligand:tether_ligand" );  } }
namespace docking { namespace ligand { RealVectorOptionKey const start_from( "docking:ligand:start_from" );  } }
namespace docking { namespace ligand { StringOptionKey const option_file( "docking:ligand:option_file" );  } }
namespace docking { namespace ligand { BooleanOptionKey const rescore( "docking:ligand:rescore" );  } }
namespace docking { namespace ligand { namespace grid { BooleanOptionKey const grid( "docking:ligand:grid" );  } } }
namespace docking { namespace ligand { namespace grid { FileOptionKey const grid_kin( "docking:ligand:grid:grid_kin" );  } } }
namespace docking { namespace ligand { namespace grid { FileOptionKey const grid_map( "docking:ligand:grid:grid_map" );  } } }
namespace docking { namespace symmetry { BooleanOptionKey const symmetry( "docking:symmetry" );  } }
namespace docking { namespace symmetry { BooleanOptionKey const minimize_backbone( "docking:symmetry:minimize_backbone" );  } }
namespace docking { namespace symmetry { BooleanOptionKey const minimize_sidechains( "docking:symmetry:minimize_sidechains" );  } }
namespace pH { BooleanOptionKey const pH( "pH" );  }
namespace pH { BooleanOptionKey const pH_mode( "pH:pH_mode" );  }
namespace pH { BooleanOptionKey const keep_input_protonation_state( "pH:keep_input_protonation_state" );  }
namespace pH { RealOptionKey const value_pH( "pH:value_pH" );  }
namespace pH { namespace calc_pka { BooleanOptionKey const calc_pka( "pH:calc_pka" );  } }
namespace pH { namespace calc_pka { BooleanOptionKey const pka_all( "pH:calc_pka:pka_all" );  } }
namespace pH { namespace calc_pka { RealVectorOptionKey const pka_for_resnos( "pH:calc_pka:pka_for_resnos" );  } }
namespace pH { namespace calc_pka { StringOptionKey const pka_for_chainno( "pH:calc_pka:pka_for_chainno" );  } }
namespace pH { namespace calc_pka { BooleanOptionKey const pH_neighbor_pack( "pH:calc_pka:pH_neighbor_pack" );  } }
namespace pH { namespace calc_pka { RealOptionKey const pka_rad( "pH:calc_pka:pka_rad" );  } }
namespace pH { namespace calc_pka { BooleanOptionKey const pH_prepack( "pH:calc_pka:pH_prepack" );  } }
namespace pH { namespace calc_pka { BooleanOptionKey const pH_relax( "pH:calc_pka:pH_relax" );  } }
namespace pH { namespace calc_pka { BooleanOptionKey const rotamer_prot_stats( "pH:calc_pka:rotamer_prot_stats" );  } }
namespace pH { FileVectorOptionKey const pH_unbound( "pH:pH_unbound" );  }
namespace pH { BooleanOptionKey const output_raw_scores( "pH:output_raw_scores" );  }
namespace pH { BooleanOptionKey const pre_process( "pH:pre_process" );  }
namespace pH { StringOptionKey const cognate_partners( "pH:cognate_partners" );  }
namespace pH { FileOptionKey const cognate_pdb( "pH:cognate_pdb" );  }
namespace run { BooleanOptionKey const run( "run" );  }
namespace run { FileVectorOptionKey const batches( "run:batches" );  }
namespace run { BooleanOptionKey const no_prof_info_in_silentout( "run:no_prof_info_in_silentout" );  }
namespace run { BooleanOptionKey const archive( "run:archive" );  }
namespace run { IntegerOptionKey const n_replica( "run:n_replica" );  }
namespace run { BooleanOptionKey const shuffle( "run:shuffle" );  }
namespace run { IntegerOptionKey const n_cycles( "run:n_cycles" );  }
namespace run { IntegerOptionKey const repeat( "run:repeat" );  }
namespace run { IntegerOptionKey const max_min_iter( "run:max_min_iter" );  }
namespace run { IntegerOptionKey const maxruntime( "run:maxruntime" );  }
namespace run { BooleanOptionKey const write_failures( "run:write_failures" );  }
namespace run { BooleanOptionKey const clean( "run:clean" );  }
namespace run { BooleanOptionKey const benchmark( "run:benchmark" );  }
namespace run { BooleanOptionKey const test_cycles( "run:test_cycles" );  }
namespace run { BooleanOptionKey const memory_test_cycles( "run:memory_test_cycles" );  }
namespace run { BooleanOptionKey const dry_run( "run:dry_run" );  }
namespace run { BooleanOptionKey const debug( "run:debug" );  }
namespace run { BooleanOptionKey const profile( "run:profile" );  }
namespace run { IntegerOptionKey const max_retry_job( "run:max_retry_job" );  }
namespace run { IntegerOptionKey const verbosity( "run:verbosity" );  }
namespace run { BooleanOptionKey const version( "run:version" );  }
namespace run { BooleanOptionKey const nodelay( "run:nodelay" );  }
namespace run { IntegerOptionKey const delay( "run:delay" );  }
namespace run { IntegerOptionKey const random_delay( "run:random_delay" );  }
namespace run { BooleanOptionKey const timer( "run:timer" );  }
namespace run { StringOptionKey const series( "run:series" );  }
namespace run { StringOptionKey const protein( "run:protein" );  }
namespace run { StringOptionKey const chain( "run:chain" );  }
namespace run { BooleanOptionKey const score_only( "run:score_only" );  }
namespace run { BooleanOptionKey const silent_input( "run:silent_input" );  }
namespace run { BooleanOptionKey const decoystats( "run:decoystats" );  }
namespace run { BooleanOptionKey const output_hbond_info( "run:output_hbond_info" );  }
namespace run { RealOptionKey const wide_nblist_extension( "run:wide_nblist_extension" );  }
namespace run { BooleanOptionKey const status( "run:status" );  }
namespace run { BooleanOptionKey const constant_seed( "run:constant_seed" );  }
namespace run { IntegerOptionKey const jran( "run:jran" );  }
namespace run { BooleanOptionKey const use_time_as_seed( "run:use_time_as_seed" );  }
namespace run { StringOptionKey const rng_seed_device( "run:rng_seed_device" );  }
namespace run { IntegerOptionKey const seed_offset( "run:seed_offset" );  }
namespace run { StringOptionKey const rng( "run:rng" );  }
namespace run { IntegerOptionKey const run_level( "run:run_level" );  }
namespace run { StringOptionKey const verbose( "run:verbose" );  }
namespace run { BooleanOptionKey const silent( "run:silent" );  }
namespace run { BooleanOptionKey const regions( "run:regions" );  }
namespace run { BooleanOptionKey const find_disulf( "run:find_disulf" );  }
namespace run { BooleanOptionKey const rebuild_disulf( "run:rebuild_disulf" );  }
namespace run { BooleanOptionKey const movie( "run:movie" );  }
namespace run { BooleanOptionKey const trajectory( "run:trajectory" );  }
namespace run { BooleanOptionKey const IUPAC( "run:IUPAC" );  }
namespace run { BooleanOptionKey const preserve_header( "run:preserve_header" );  }
namespace run { BooleanOptionKey const evolution( "run:evolution" );  }
namespace run { BooleanOptionKey const suppress_checkpoints( "run:suppress_checkpoints" );  }
namespace run { BooleanOptionKey const checkpoint( "run:checkpoint" );  }
namespace run { BooleanOptionKey const delete_checkpoints( "run:delete_checkpoints" );  }
namespace run { IntegerOptionKey const checkpoint_interval( "run:checkpoint_interval" );  }
namespace run { StringOptionKey const protocol( "run:protocol" );  }
namespace run { BooleanOptionKey const remove_ss_length_screen( "run:remove_ss_length_screen" );  }
namespace run { StringOptionKey const min_type( "run:min_type" );  }
namespace run { RealOptionKey const min_tolerance( "run:min_tolerance" );  }
namespace run { BooleanOptionKey const nblist_autoupdate( "run:nblist_autoupdate" );  }
namespace run { RealOptionKey const nblist_autoupdate_narrow( "run:nblist_autoupdate_narrow" );  }
namespace run { RealOptionKey const nblist_autoupdate_wide( "run:nblist_autoupdate_wide" );  }
namespace run { BooleanOptionKey const skip_set_reasonable_fold_tree( "run:skip_set_reasonable_fold_tree" );  }
namespace run { BooleanOptionKey const randomize_missing_coords( "run:randomize_missing_coords" );  }
namespace run { BooleanOptionKey const ignore_zero_occupancy( "run:ignore_zero_occupancy" );  }
namespace run { IntegerOptionKey const cycles_outer( "run:cycles_outer" );  }
namespace run { IntegerOptionKey const cycles_inner( "run:cycles_inner" );  }
namespace run { IntegerOptionKey const repack_rate( "run:repack_rate" );  }
namespace run { BooleanOptionKey const reinitialize_mover_for_each_job( "run:reinitialize_mover_for_each_job" );  }
namespace run { BooleanOptionKey const reinitialize_mover_for_new_input( "run:reinitialize_mover_for_new_input" );  }
namespace run { BooleanOptionKey const multiple_processes_writing_to_one_directory( "run:multiple_processes_writing_to_one_directory" );  }
namespace run { StringOptionKey const jobdist_miscfile_ext( "run:jobdist_miscfile_ext" );  }
namespace run { BooleanOptionKey const no_scorefile( "run:no_scorefile" );  }
namespace run { BooleanOptionKey const other_pose_to_scorefile( "run:other_pose_to_scorefile" );  }
namespace run { FileOptionKey const other_pose_scorefile( "run:other_pose_scorefile" );  }
namespace run { BooleanOptionKey const intermediate_scorefiles( "run:intermediate_scorefiles" );  }
namespace run { BooleanOptionKey const intermediate_structures( "run:intermediate_structures" );  }
namespace run { BooleanOptionKey const idealize_before_protocol( "run:idealize_before_protocol" );  }
namespace run { BooleanOptionKey const interactive( "run:interactive" );  }
namespace run { BooleanOptionKey const condor( "run:condor" );  }
namespace run { IntegerOptionKey const nproc( "run:nproc" );  }
namespace run { IntegerOptionKey const proc_id( "run:proc_id" );  }
namespace run { BooleanOptionKey const exit_if_missing_heavy_atoms( "run:exit_if_missing_heavy_atoms" );  }
namespace run { RealOptionKey const show_simulation_in_pymol( "run:show_simulation_in_pymol" );  }
namespace run { BooleanOptionKey const keep_pymol_simulation_history( "run:keep_pymol_simulation_history" );  }
namespace jd2 { BooleanOptionKey const jd2( "jd2" );  }
namespace jd2 { BooleanOptionKey const pose_input_stream( "jd2:pose_input_stream" );  }
namespace jd2 { BooleanOptionKey const lazy_silent_file_reader( "jd2:lazy_silent_file_reader" );  }
namespace jd2 { BooleanOptionKey const mpi_nowait_for_remaining_jobs( "jd2:mpi_nowait_for_remaining_jobs" );  }
namespace jd2 { RealOptionKey const mpi_timeout_factor( "jd2:mpi_timeout_factor" );  }
namespace jd2 { BooleanOptionKey const mpi_work_partition_job_distributor( "jd2:mpi_work_partition_job_distributor" );  }
namespace jd2 { BooleanOptionKey const mpi_file_buf_job_distributor( "jd2:mpi_file_buf_job_distributor" );  }
namespace jd2 { BooleanOptionKey const mpi_filebuf_jobdistributor( "jd2:mpi_filebuf_jobdistributor" );  }
namespace jd2 { BooleanOptionKey const mpi_fast_nonblocking_output( "jd2:mpi_fast_nonblocking_output" );  }
namespace jd2 { BooleanOptionKey const dd_parser( "jd2:dd_parser" );  }
namespace jd2 { IntegerOptionKey const ntrials( "jd2:ntrials" );  }
namespace jd2 { StringOptionKey const generic_job_name( "jd2:generic_job_name" );  }
namespace jd2 { BooleanOptionKey const no_output( "jd2:no_output" );  }
namespace jd2 { BooleanOptionKey const enzdes_out( "jd2:enzdes_out" );  }
namespace jd2 { IntegerOptionKey const buffer_silent_output( "jd2:buffer_silent_output" );  }
namespace jd2 { RealOptionKey const buffer_flush_frequency( "jd2:buffer_flush_frequency" );  }
namespace jd2 { BooleanOptionKey const delete_old_poses( "jd2:delete_old_poses" );  }
namespace jd2 { FileVectorOptionKey const resource_definition_files( "jd2:resource_definition_files" );  }
namespace jd2 { FileOptionKey const checkpoint_file( "jd2:checkpoint_file" );  }
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
namespace evaluation { FileVectorOptionKey const rdc_select( "evaluation:rdc_select" );  }
namespace evaluation { FileVectorOptionKey const rdc_target( "evaluation:rdc_target" );  }
namespace evaluation { BooleanOptionKey const symmetric_rmsd( "evaluation:symmetric_rmsd" );  }
namespace evaluation { StringVectorOptionKey const rdc_column( "evaluation:rdc_column" );  }
namespace evaluation { StringVectorOptionKey const rdc( "evaluation:rdc" );  }
namespace evaluation { StringOptionKey const built_in_rdc( "evaluation:built_in_rdc" );  }
namespace evaluation { BooleanOptionKey const jump_nr( "evaluation:jump_nr" );  }
namespace evaluation { IntegerVectorOptionKey const score_exclude_res( "evaluation:score_exclude_res" );  }
namespace evaluation { IntegerOptionKey const score_sscore_short_helix( "evaluation:score_sscore_short_helix" );  }
namespace evaluation { IntegerOptionKey const score_sscore_maxloop( "evaluation:score_sscore_maxloop" );  }
namespace evaluation { BooleanOptionKey const rpf( "evaluation:rpf" );  }
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
namespace MonteCarlo { BooleanOptionKey const MonteCarlo( "MonteCarlo" );  }
namespace MonteCarlo { RealOptionKey const temp_initial( "MonteCarlo:temp_initial" );  }
namespace MonteCarlo { RealOptionKey const temp_final( "MonteCarlo:temp_final" );  }
namespace frags { BooleanOptionKey const frags( "frags" );  }
namespace frags { IntegerOptionKey const j( "frags:j" );  }
namespace frags { BooleanOptionKey const filter_JC( "frags:filter_JC" );  }
namespace frags { BooleanOptionKey const bounded_protocol( "frags:bounded_protocol" );  }
namespace frags { BooleanOptionKey const keep_all_protocol( "frags:keep_all_protocol" );  }
namespace frags { BooleanOptionKey const quota_protocol( "frags:quota_protocol" );  }
namespace frags { BooleanOptionKey const nonlocal_pairs( "frags:nonlocal_pairs" );  }
namespace frags { BooleanOptionKey const fragment_contacts( "frags:fragment_contacts" );  }
namespace frags { BooleanOptionKey const p_value_selection( "frags:p_value_selection" );  }
namespace frags { IntegerOptionKey const n_frags( "frags:n_frags" );  }
namespace frags { FileOptionKey const allowed_pdb( "frags:allowed_pdb" );  }
namespace frags { StringVectorOptionKey const ss_pred( "frags:ss_pred" );  }
namespace frags { FileOptionKey const spine_x( "frags:spine_x" );  }
namespace frags { FileOptionKey const depth( "frags:depth" );  }
namespace frags { FileOptionKey const denied_pdb( "frags:denied_pdb" );  }
namespace frags { IntegerVectorOptionKey const frag_sizes( "frags:frag_sizes" );  }
namespace frags { BooleanOptionKey const write_ca_coordinates( "frags:write_ca_coordinates" );  }
namespace frags { BooleanOptionKey const write_scores( "frags:write_scores" );  }
namespace frags { BooleanOptionKey const annotate( "frags:annotate" );  }
namespace frags { IntegerOptionKey const nr_large_copies( "frags:nr_large_copies" );  }
namespace frags { IntegerOptionKey const n_candidates( "frags:n_candidates" );  }
namespace frags { BooleanOptionKey const write_rama_tables( "frags:write_rama_tables" );  }
namespace frags { RealOptionKey const rama_C( "frags:rama_C" );  }
namespace frags { RealOptionKey const rama_B( "frags:rama_B" );  }
namespace frags { RealOptionKey const sigmoid_cs_A( "frags:sigmoid_cs_A" );  }
namespace frags { RealOptionKey const sigmoid_cs_B( "frags:sigmoid_cs_B" );  }
namespace frags { RealOptionKey const seqsim_H( "frags:seqsim_H" );  }
namespace frags { RealOptionKey const seqsim_E( "frags:seqsim_E" );  }
namespace frags { RealOptionKey const seqsim_L( "frags:seqsim_L" );  }
namespace frags { RealOptionKey const rama_norm( "frags:rama_norm" );  }
namespace frags { StringOptionKey const describe_fragments( "frags:describe_fragments" );  }
namespace frags { RealOptionKey const picking_old_max_score( "frags:picking_old_max_score" );  }
namespace frags { BooleanOptionKey const write_sequence_only( "frags:write_sequence_only" );  }
namespace frags { BooleanOptionKey const output_silent( "frags:output_silent" );  }
namespace frags { BooleanOptionKey const score_output_silent( "frags:score_output_silent" );  }
namespace frags { namespace scoring { BooleanOptionKey const scoring( "frags:scoring" );  } }
namespace frags { namespace scoring { FileOptionKey const config( "frags:scoring:config" );  } }
namespace frags { namespace scoring { StringOptionKey const profile_score( "frags:scoring:profile_score" );  } }
namespace frags { namespace scoring { FileVectorOptionKey const predicted_secondary( "frags:scoring:predicted_secondary" );  } }
namespace frags { namespace picking { BooleanOptionKey const picking( "frags:picking" );  } }
namespace frags { namespace picking { StringOptionKey const selecting_rule( "frags:picking:selecting_rule" );  } }
namespace frags { namespace picking { StringOptionKey const selecting_scorefxn( "frags:picking:selecting_scorefxn" );  } }
namespace frags { namespace picking { FileOptionKey const quota_config_file( "frags:picking:quota_config_file" );  } }
namespace frags { namespace picking { IntegerVectorOptionKey const query_pos( "frags:picking:query_pos" );  } }
namespace frags { namespace nonlocal { BooleanOptionKey const nonlocal( "frags:nonlocal" );  } }
namespace frags { namespace nonlocal { BooleanOptionKey const relax_input( "frags:nonlocal:relax_input" );  } }
namespace frags { namespace nonlocal { BooleanOptionKey const relax_input_with_coordinate_constraints( "frags:nonlocal:relax_input_with_coordinate_constraints" );  } }
namespace frags { namespace nonlocal { IntegerOptionKey const relax_frags_repeats( "frags:nonlocal:relax_frags_repeats" );  } }
namespace frags { namespace nonlocal { BooleanOptionKey const single_chain( "frags:nonlocal:single_chain" );  } }
namespace frags { namespace nonlocal { RealOptionKey const min_contacts_per_res( "frags:nonlocal:min_contacts_per_res" );  } }
namespace frags { namespace nonlocal { RealOptionKey const max_ddg_score( "frags:nonlocal:max_ddg_score" );  } }
namespace frags { namespace nonlocal { RealOptionKey const max_rmsd_after_relax( "frags:nonlocal:max_rmsd_after_relax" );  } }
namespace frags { namespace nonlocal { BooleanOptionKey const output_frags_pdbs( "frags:nonlocal:output_frags_pdbs" );  } }
namespace frags { namespace nonlocal { BooleanOptionKey const output_idealized( "frags:nonlocal:output_idealized" );  } }
namespace frags { namespace nonlocal { BooleanOptionKey const output_silent( "frags:nonlocal:output_silent" );  } }
namespace frags { namespace contacts { BooleanOptionKey const contacts( "frags:contacts" );  } }
namespace frags { namespace contacts { IntegerOptionKey const min_seq_sep( "frags:contacts:min_seq_sep" );  } }
namespace frags { namespace contacts { RealVectorOptionKey const dist_cutoffs( "frags:contacts:dist_cutoffs" );  } }
namespace frags { namespace contacts { RealOptionKey const centroid_distance_scale_factor( "frags:contacts:centroid_distance_scale_factor" );  } }
namespace frags { namespace contacts { StringVectorOptionKey const type( "frags:contacts:type" );  } }
namespace frags { namespace contacts { IntegerOptionKey const neighbors( "frags:contacts:neighbors" );  } }
namespace frags { namespace contacts { BooleanOptionKey const output_all( "frags:contacts:output_all" );  } }
namespace frags { namespace ABEGO { BooleanOptionKey const ABEGO( "frags:ABEGO" );  } }
namespace frags { namespace ABEGO { RealOptionKey const phi_psi_range_A( "frags:ABEGO:phi_psi_range_A" );  } }
namespace broker { BooleanOptionKey const broker( "broker" );  }
namespace broker { FileVectorOptionKey const setup( "broker:setup" );  }
namespace chunk { BooleanOptionKey const chunk( "chunk" );  }
namespace chunk { FileOptionKey const pdb2( "chunk:pdb2" );  }
namespace chunk { FileOptionKey const loop2( "chunk:loop2" );  }
namespace nonlocal { BooleanOptionKey const nonlocal( "nonlocal" );  }
namespace nonlocal { StringOptionKey const builder( "nonlocal:builder" );  }
namespace nonlocal { FileOptionKey const chunks( "nonlocal:chunks" );  }
namespace nonlocal { IntegerOptionKey const max_chunk_size( "nonlocal:max_chunk_size" );  }
namespace nonlocal { BooleanOptionKey const randomize_missing( "nonlocal:randomize_missing" );  }
namespace nonlocal { RealOptionKey const rdc_weight( "nonlocal:rdc_weight" );  }
namespace abinitio { namespace star { BooleanOptionKey const star( "abinitio:star" );  } }
namespace abinitio { namespace star { RealOptionKey const initial_dist_cutoff( "abinitio:star:initial_dist_cutoff" );  } }
namespace abinitio { namespace star { IntegerOptionKey const min_unaligned_len( "abinitio:star:min_unaligned_len" );  } }
namespace abinitio { namespace star { IntegerOptionKey const short_loop_len( "abinitio:star:short_loop_len" );  } }
namespace abinitio { RealOptionKey const prob_perturb_weights( "abinitio:prob_perturb_weights" );  }
namespace abinitio { BooleanOptionKey const abinitio( "abinitio" );  }
namespace abinitio { BooleanOptionKey const membrane( "abinitio:membrane" );  }
namespace abinitio { FileOptionKey const kill_hairpins( "abinitio:kill_hairpins" );  }
namespace abinitio { RealOptionKey const kill_hairpins_frequency( "abinitio:kill_hairpins_frequency" );  }
namespace abinitio { BooleanOptionKey const smooth_cycles_only( "abinitio:smooth_cycles_only" );  }
namespace abinitio { BooleanOptionKey const relax( "abinitio:relax" );  }
namespace abinitio { BooleanOptionKey const final_clean_relax( "abinitio:final_clean_relax" );  }
namespace abinitio { BooleanOptionKey const fastrelax( "abinitio:fastrelax" );  }
namespace abinitio { BooleanOptionKey const multifastrelax( "abinitio:multifastrelax" );  }
namespace abinitio { BooleanOptionKey const debug( "abinitio:debug" );  }
namespace abinitio { BooleanOptionKey const clear_pose_cache( "abinitio:clear_pose_cache" );  }
namespace abinitio { BooleanOptionKey const explicit_pdb_debug( "abinitio:explicit_pdb_debug" );  }
namespace abinitio { BooleanOptionKey const use_filters( "abinitio:use_filters" );  }
namespace abinitio { RealOptionKey const increase_cycles( "abinitio:increase_cycles" );  }
namespace abinitio { IntegerOptionKey const number_3mer_frags( "abinitio:number_3mer_frags" );  }
namespace abinitio { IntegerOptionKey const number_9mer_frags( "abinitio:number_9mer_frags" );  }
namespace abinitio { RealOptionKey const temperature( "abinitio:temperature" );  }
namespace abinitio { RealOptionKey const rg_reweight( "abinitio:rg_reweight" );  }
namespace abinitio { RealOptionKey const strand_dist_cutoff( "abinitio:strand_dist_cutoff" );  }
namespace abinitio { BooleanOptionKey const stretch_strand_dist_cutoff( "abinitio:stretch_strand_dist_cutoff" );  }
namespace abinitio { RealOptionKey const rsd_wt_helix( "abinitio:rsd_wt_helix" );  }
namespace abinitio { RealOptionKey const rsd_wt_strand( "abinitio:rsd_wt_strand" );  }
namespace abinitio { RealOptionKey const rsd_wt_loop( "abinitio:rsd_wt_loop" );  }
namespace abinitio { BooleanOptionKey const fast( "abinitio:fast" );  }
namespace abinitio { BooleanOptionKey const skip_convergence_check( "abinitio:skip_convergence_check" );  }
namespace abinitio { FileVectorOptionKey const stage1_patch( "abinitio:stage1_patch" );  }
namespace abinitio { FileVectorOptionKey const stage2_patch( "abinitio:stage2_patch" );  }
namespace abinitio { FileVectorOptionKey const stage3a_patch( "abinitio:stage3a_patch" );  }
namespace abinitio { FileVectorOptionKey const stage3b_patch( "abinitio:stage3b_patch" );  }
namespace abinitio { FileVectorOptionKey const stage4_patch( "abinitio:stage4_patch" );  }
namespace abinitio { FileVectorOptionKey const stage5_patch( "abinitio:stage5_patch" );  }
namespace abinitio { BooleanOptionKey const exit_when_converged( "abinitio:exit_when_converged" );  }
namespace abinitio { BooleanOptionKey const steal_3mers( "abinitio:steal_3mers" );  }
namespace abinitio { BooleanOptionKey const steal_9mers( "abinitio:steal_9mers" );  }
namespace abinitio { BooleanOptionKey const no_write_failures( "abinitio:no_write_failures" );  }
namespace abinitio { BooleanOptionKey const relax_failures( "abinitio:relax_failures" );  }
namespace abinitio { BooleanOptionKey const relax_with_jumps( "abinitio:relax_with_jumps" );  }
namespace abinitio { BooleanOptionKey const process_store( "abinitio:process_store" );  }
namespace abinitio { IntegerVectorOptionKey const fix_residues_to_native( "abinitio:fix_residues_to_native" );  }
namespace abinitio { BooleanOptionKey const return_full_atom( "abinitio:return_full_atom" );  }
namespace abinitio { BooleanOptionKey const detect_disulfide_before_relax( "abinitio:detect_disulfide_before_relax" );  }
namespace abinitio { BooleanOptionKey const close_loops( "abinitio:close_loops" );  }
namespace abinitio { BooleanOptionKey const bGDT( "abinitio:bGDT" );  }
namespace abinitio { BooleanOptionKey const dump_frags( "abinitio:dump_frags" );  }
namespace abinitio { BooleanOptionKey const jdist_rerun( "abinitio:jdist_rerun" );  }
namespace abinitio { RealOptionKey const perturb( "abinitio:perturb" );  }
namespace abinitio { BooleanOptionKey const rerun( "abinitio:rerun" );  }
namespace abinitio { IntegerVectorOptionKey const rmsd_residues( "abinitio:rmsd_residues" );  }
namespace abinitio { BooleanOptionKey const start_native( "abinitio:start_native" );  }
namespace abinitio { BooleanOptionKey const debug_structures( "abinitio:debug_structures" );  }
namespace abinitio { FileOptionKey const log_frags( "abinitio:log_frags" );  }
namespace abinitio { BooleanOptionKey const only_stage1( "abinitio:only_stage1" );  }
namespace abinitio { RealOptionKey const end_bias( "abinitio:end_bias" );  }
namespace abinitio { IntegerOptionKey const symmetry_residue( "abinitio:symmetry_residue" );  }
namespace abinitio { RealOptionKey const vdw_weight_stage1( "abinitio:vdw_weight_stage1" );  }
namespace abinitio { BooleanOptionKey const override_vdw_all_stages( "abinitio:override_vdw_all_stages" );  }
namespace abinitio { IntegerVectorOptionKey const recover_low_in_stages( "abinitio:recover_low_in_stages" );  }
namespace abinitio { IntegerVectorOptionKey const skip_stages( "abinitio:skip_stages" );  }
namespace abinitio { BooleanOptionKey const close_chbrk( "abinitio:close_chbrk" );  }
namespace abinitio { BooleanOptionKey const include_stage5( "abinitio:include_stage5" );  }
namespace abinitio { BooleanOptionKey const close_loops_by_idealizing( "abinitio:close_loops_by_idealizing" );  }
namespace abinitio { BooleanOptionKey const optimize_cutpoints_using_kic( "abinitio:optimize_cutpoints_using_kic" );  }
namespace abinitio { IntegerOptionKey const optimize_cutpoints_margin( "abinitio:optimize_cutpoints_margin" );  }
namespace abinitio { FileOptionKey const HD_EX_Info( "abinitio:HD_EX_Info" );  }
namespace abinitio { RealOptionKey const HD_penalty( "abinitio:HD_penalty" );  }
namespace abinitio { RealOptionKey const HD_fa_penalty( "abinitio:HD_fa_penalty" );  }
namespace abinitio { FileOptionKey const sheet_edge_pred( "abinitio:sheet_edge_pred" );  }
namespace abinitio { RealOptionKey const SEP_score_scalling( "abinitio:SEP_score_scalling" );  }
namespace fold_cst { BooleanOptionKey const fold_cst( "fold_cst" );  }
namespace fold_cst { RealOptionKey const constraint_skip_rate( "fold_cst:constraint_skip_rate" );  }
namespace fold_cst { IntegerOptionKey const violation_skip_basis( "fold_cst:violation_skip_basis" );  }
namespace fold_cst { IntegerOptionKey const violation_skip_ignore( "fold_cst:violation_skip_ignore" );  }
namespace fold_cst { BooleanOptionKey const keep_skipped_csts( "fold_cst:keep_skipped_csts" );  }
namespace fold_cst { BooleanOptionKey const no_minimize( "fold_cst:no_minimize" );  }
namespace fold_cst { BooleanOptionKey const force_minimize( "fold_cst:force_minimize" );  }
namespace fold_cst { RealVectorOptionKey const seq_sep_stages( "fold_cst:seq_sep_stages" );  }
namespace fold_cst { IntegerOptionKey const reramp_cst_cycles( "fold_cst:reramp_cst_cycles" );  }
namespace fold_cst { RealOptionKey const reramp_start_cstweight( "fold_cst:reramp_start_cstweight" );  }
namespace fold_cst { IntegerOptionKey const reramp_iterations( "fold_cst:reramp_iterations" );  }
namespace fold_cst { BooleanOptionKey const skip_on_noviolation_in_stage1( "fold_cst:skip_on_noviolation_in_stage1" );  }
namespace fold_cst { RealOptionKey const stage1_ramp_cst_cycle_factor( "fold_cst:stage1_ramp_cst_cycle_factor" );  }
namespace fold_cst { RealOptionKey const stage2_constraint_threshold( "fold_cst:stage2_constraint_threshold" );  }
namespace fold_cst { BooleanOptionKey const ignore_sequence_seperation( "fold_cst:ignore_sequence_seperation" );  }
namespace fold_cst { BooleanOptionKey const no_recover_low_at_constraint_switch( "fold_cst:no_recover_low_at_constraint_switch" );  }
namespace fold_cst { BooleanOptionKey const ramp_coord_cst( "fold_cst:ramp_coord_cst" );  }
namespace resample { BooleanOptionKey const resample( "resample" );  }
namespace resample { FileOptionKey const silent( "resample:silent" );  }
namespace resample { StringOptionKey const tag( "resample:tag" );  }
namespace resample { BooleanOptionKey const stage1( "resample:stage1" );  }
namespace resample { BooleanOptionKey const stage2( "resample:stage2" );  }
namespace resample { BooleanOptionKey const jumps( "resample:jumps" );  }
namespace resample { RealVectorOptionKey const min_max_start_seq_sep( "resample:min_max_start_seq_sep" );  }
namespace loopfcst { BooleanOptionKey const loopfcst( "loopfcst" );  }
namespace loopfcst { RealOptionKey const coord_cst_weight( "loopfcst:coord_cst_weight" );  }
namespace loopfcst { BooleanOptionKey const coord_cst_all_atom( "loopfcst:coord_cst_all_atom" );  }
namespace loopfcst { BooleanOptionKey const use_general_protocol( "loopfcst:use_general_protocol" );  }
namespace loopfcst { FileOptionKey const coord_cst_weight_array( "loopfcst:coord_cst_weight_array" );  }
namespace loopfcst { FileOptionKey const dump_coord_cst_weight_array( "loopfcst:dump_coord_cst_weight_array" );  }
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
