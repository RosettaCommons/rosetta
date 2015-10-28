namespace score { BooleanOptionKey const allow_complex_loop_graph( "score:allow_complex_loop_graph" );  }
namespace score { BooleanOptionKey const compute_mg_sol_for_hydrogens( "score:compute_mg_sol_for_hydrogens" );  }
namespace score { IntegerVectorOptionKey const rg_local_span( "score:rg_local_span" );  }
namespace score { BooleanOptionKey const unmodifypot( "score:unmodifypot" );  }
namespace score { RealOptionKey const conc( "score:conc" );  }
namespace score { RealOptionKey const loop_fixed_cost( "score:loop_fixed_cost" );  }
namespace score { IntegerVectorOptionKey const sidechain_buried( "score:sidechain_buried" );  }
namespace score { IntegerVectorOptionKey const sidechain_exposed( "score:sidechain_exposed" );  }
namespace score { StringOptionKey const aa_composition_setup_file( "score:aa_composition_setup_file" );  }
namespace score { StringOptionKey const aa_repeat_energy_penalty_file( "score:aa_repeat_energy_penalty_file" );  }
namespace score { RealOptionKey const elec_min_dis( "score:elec_min_dis" );  }
namespace score { RealOptionKey const elec_max_dis( "score:elec_max_dis" );  }
namespace score { RealOptionKey const elec_die( "score:elec_die" );  }
namespace score { BooleanOptionKey const elec_r_option( "score:elec_r_option" );  }
namespace score { BooleanOptionKey const elec_sigmoidal_die( "score:elec_sigmoidal_die" );  }
namespace score { RealOptionKey const elec_sigmoidal_die_D( "score:elec_sigmoidal_die_D" );  }
namespace score { RealOptionKey const elec_sigmoidal_die_D0( "score:elec_sigmoidal_die_D0" );  }
namespace score { RealOptionKey const elec_sigmoidal_die_S( "score:elec_sigmoidal_die_S" );  }
namespace score { RealOptionKey const intrares_elec_correction_scale( "score:intrares_elec_correction_scale" );  }
namespace score { BooleanOptionKey const smooth_fa_elec( "score:smooth_fa_elec" );  }
namespace score { StringOptionKey const grpelec_fade_type( "score:grpelec_fade_type" );  }
namespace score { RealOptionKey const grpelec_fade_param1( "score:grpelec_fade_param1" );  }
namespace score { RealOptionKey const grpelec_fade_param2( "score:grpelec_fade_param2" );  }
namespace score { StringOptionKey const elec_group_file( "score:elec_group_file" );  }
namespace score { StringOptionKey const elec_group_extrafile( "score:elec_group_extrafile" );  }
namespace score { BooleanOptionKey const grpelec_fade_hbond( "score:grpelec_fade_hbond" );  }
namespace score { RealVectorOptionKey const grpelec_max_qeps( "score:grpelec_max_qeps" );  }
namespace score { BooleanOptionKey const grpelec_context_dependent( "score:grpelec_context_dependent" );  }
namespace score { BooleanOptionKey const grp_cpfxn( "score:grp_cpfxn" );  }
namespace score { RealVectorOptionKey const grpelec_cpfxn_weight( "score:grpelec_cpfxn_weight" );  }
namespace score { RealOptionKey const elec_context_minstrength( "score:elec_context_minstrength" );  }
namespace score { RealOptionKey const elec_context_minburial( "score:elec_context_minburial" );  }
namespace score { RealOptionKey const elec_context_maxburial( "score:elec_context_maxburial" );  }
namespace score { BooleanOptionKey const use_polarization( "score:use_polarization" );  }
namespace score { BooleanOptionKey const use_gen_kirkwood( "score:use_gen_kirkwood" );  }
namespace score { RealOptionKey const protein_dielectric( "score:protein_dielectric" );  }
namespace score { RealOptionKey const water_dielectric( "score:water_dielectric" );  }
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
namespace score { BooleanOptionKey const length_dep_srbb( "score:length_dep_srbb" );  }
namespace score { RealOptionKey const ldsrbb_low_scale( "score:ldsrbb_low_scale" );  }
namespace score { RealOptionKey const ldsrbb_high_scale( "score:ldsrbb_high_scale" );  }
namespace score { IntegerOptionKey const ldsrbb_minlength( "score:ldsrbb_minlength" );  }
namespace score { IntegerOptionKey const ldsrbb_maxlength( "score:ldsrbb_maxlength" );  }
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
namespace score { BooleanOptionKey const envsmooth_zero_negatives( "score:envsmooth_zero_negatives" );  }
namespace score { namespace saxs { BooleanOptionKey const saxs( "score:saxs" );  } }
namespace score { namespace saxs { RealOptionKey const min_score( "score:saxs:min_score" );  } }
namespace score { namespace saxs { StringOptionKey const custom_ff( "score:saxs:custom_ff" );  } }
namespace score { namespace saxs { StringOptionKey const print_i_calc( "score:saxs:print_i_calc" );  } }
namespace score { namespace saxs { FileOptionKey const ref_fa_spectrum( "score:saxs:ref_fa_spectrum" );  } }
namespace score { namespace saxs { FileOptionKey const ref_cen_spectrum( "score:saxs:ref_cen_spectrum" );  } }
namespace score { namespace saxs { FileOptionKey const ref_spectrum( "score:saxs:ref_spectrum" );  } }
namespace score { namespace saxs { FileOptionKey const ref_pddf( "score:saxs:ref_pddf" );  } }
namespace score { namespace saxs { RealOptionKey const d_min( "score:saxs:d_min" );  } }
namespace score { namespace saxs { RealOptionKey const d_max( "score:saxs:d_max" );  } }
namespace score { namespace saxs { RealOptionKey const d_step( "score:saxs:d_step" );  } }
namespace score { namespace saxs { RealOptionKey const q_min( "score:saxs:q_min" );  } }
namespace score { namespace saxs { RealOptionKey const q_max( "score:saxs:q_max" );  } }
namespace score { namespace saxs { RealOptionKey const q_step( "score:saxs:q_step" );  } }
namespace score { namespace saxs { BooleanOptionKey const fit_pddf_area( "score:saxs:fit_pddf_area" );  } }
namespace score { namespace fiber_diffraction { BooleanOptionKey const fiber_diffraction( "score:fiber_diffraction" );  } }
namespace score { namespace fiber_diffraction { FileOptionKey const layer_lines( "score:fiber_diffraction:layer_lines" );  } }
namespace score { namespace fiber_diffraction { RealOptionKey const a( "score:fiber_diffraction:a" );  } }
namespace score { namespace fiber_diffraction { RealOptionKey const b( "score:fiber_diffraction:b" );  } }
namespace score { namespace fiber_diffraction { RealOptionKey const p( "score:fiber_diffraction:p" );  } }
namespace score { namespace fiber_diffraction { RealOptionKey const radius( "score:fiber_diffraction:radius" );  } }
namespace score { namespace fiber_diffraction { RealOptionKey const resolution_cutoff_low( "score:fiber_diffraction:resolution_cutoff_low" );  } }
namespace score { namespace fiber_diffraction { RealOptionKey const resolution_cutoff_high( "score:fiber_diffraction:resolution_cutoff_high" );  } }
namespace score { namespace fiber_diffraction { IntegerOptionKey const max_bessel_order( "score:fiber_diffraction:max_bessel_order" );  } }
namespace score { namespace fiber_diffraction { IntegerOptionKey const cn_symmetry( "score:fiber_diffraction:cn_symmetry" );  } }
namespace score { namespace fiber_diffraction { RealOptionKey const b_factor( "score:fiber_diffraction:b_factor" );  } }
namespace score { namespace fiber_diffraction { RealOptionKey const b_factor_solv( "score:fiber_diffraction:b_factor_solv" );  } }
namespace score { namespace fiber_diffraction { RealOptionKey const b_factor_solv_K( "score:fiber_diffraction:b_factor_solv_K" );  } }
namespace score { namespace fiber_diffraction { RealOptionKey const grid_reso( "score:fiber_diffraction:grid_reso" );  } }
namespace score { namespace fiber_diffraction { IntegerOptionKey const grid_r( "score:fiber_diffraction:grid_r" );  } }
namespace score { namespace fiber_diffraction { IntegerOptionKey const grid_phi( "score:fiber_diffraction:grid_phi" );  } }
namespace score { namespace fiber_diffraction { IntegerOptionKey const grid_z( "score:fiber_diffraction:grid_z" );  } }
namespace score { namespace fiber_diffraction { RealOptionKey const qfht_K1( "score:fiber_diffraction:qfht_K1" );  } }
namespace score { namespace fiber_diffraction { RealOptionKey const qfht_K2( "score:fiber_diffraction:qfht_K2" );  } }
namespace score { namespace fiber_diffraction { IntegerOptionKey const chi_iterations( "score:fiber_diffraction:chi_iterations" );  } }
namespace score { namespace fiber_diffraction { BooleanOptionKey const rfactor_refinement( "score:fiber_diffraction:rfactor_refinement" );  } }
namespace score { namespace fiber_diffraction { BooleanOptionKey const output_fiber_spectra( "score:fiber_diffraction:output_fiber_spectra" );  } }
namespace score { namespace fiber_diffraction { IntegerOptionKey const gpu_processor( "score:fiber_diffraction:gpu_processor" );  } }
namespace score { namespace fiber_diffraction { RealOptionKey const centroid_density_mass( "score:fiber_diffraction:centroid_density_mass" );  } }
namespace packing { BooleanOptionKey const packing( "packing" );  }
namespace packing { BooleanOptionKey const repack_only( "packing:repack_only" );  }
namespace packing { BooleanOptionKey const prevent_repacking( "packing:prevent_repacking" );  }
namespace packing { RealOptionKey const cenrot_cutoff( "packing:cenrot_cutoff" );  }
namespace packing { BooleanOptionKey const ignore_ligand_chi( "packing:ignore_ligand_chi" );  }
namespace packing { IntegerOptionKey const ndruns( "packing:ndruns" );  }
namespace packing { BooleanOptionKey const soft_rep_design( "packing:soft_rep_design" );  }
namespace packing { RealOptionKey const mainchain_h_rebuild_threshold( "packing:mainchain_h_rebuild_threshold" );  }
namespace packing { BooleanOptionKey const use_electrostatic_repulsion( "packing:use_electrostatic_repulsion" );  }
namespace packing { BooleanOptionKey const dump_rotamer_sets( "packing:dump_rotamer_sets" );  }
namespace packing { RealOptionKey const dunbrack_prob_buried( "packing:dunbrack_prob_buried" );  }
namespace packing { RealOptionKey const dunbrack_prob_nonburied( "packing:dunbrack_prob_nonburied" );  }
namespace packing { BooleanOptionKey const no_optH( "packing:no_optH" );  }
namespace packing { BooleanOptionKey const optH_MCA( "packing:optH_MCA" );  }
namespace packing { BooleanOptionKey const pack_missing_sidechains( "packing:pack_missing_sidechains" );  }
namespace packing { BooleanOptionKey const preserve_c_beta( "packing:preserve_c_beta" );  }
namespace packing { BooleanOptionKey const flip_HNQ( "packing:flip_HNQ" );  }
namespace packing { IntegerVectorOptionKey const fix_his_tautomer( "packing:fix_his_tautomer" );  }
namespace packing { BooleanOptionKey const print_pymol_selection( "packing:print_pymol_selection" );  }
namespace packing { IntegerOptionKey const extrachi_cutoff( "packing:extrachi_cutoff" );  }
namespace packing { FileVectorOptionKey const resfile( "packing:resfile" );  }
namespace packing { RealOptionKey const outeriterations_scaling( "packing:outeriterations_scaling" );  }
namespace packing { RealOptionKey const inneriterations_scaling( "packing:inneriterations_scaling" );  }
namespace packing { StringVectorOptionKey const adducts( "packing:adducts" );  }
namespace packing { BooleanOptionKey const use_input_sc( "packing:use_input_sc" );  }
namespace packing { FileVectorOptionKey const unboundrot( "packing:unboundrot" );  }
namespace packing { RealOptionKey const max_rotbump_energy( "packing:max_rotbump_energy" );  }
namespace packing { BooleanOptionKey const lazy_ig( "packing:lazy_ig" );  }
namespace packing { BooleanOptionKey const double_lazy_ig( "packing:double_lazy_ig" );  }
namespace packing { IntegerOptionKey const linmem_ig( "packing:linmem_ig" );  }
namespace packing { IntegerOptionKey const multi_cool_annealer( "packing:multi_cool_annealer" );  }
namespace packing { RealVectorOptionKey const minpack_temp_schedule( "packing:minpack_temp_schedule" );  }
namespace packing { IntegerOptionKey const minpack_inner_iteration_scale( "packing:minpack_inner_iteration_scale" );  }
namespace packing { BooleanOptionKey const minpack_disable_bumpcheck( "packing:minpack_disable_bumpcheck" );  }
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
namespace archive { BooleanOptionKey const archive( "archive" );  }
namespace archive { BooleanOptionKey const reread_all_structures( "archive:reread_all_structures" );  }
namespace archive { IntegerOptionKey const completion_notify_frequency( "archive:completion_notify_frequency" );  }
namespace carbohydrates { BooleanOptionKey const carbohydrates( "carbohydrates" );  }
namespace carbohydrates { BooleanOptionKey const glycam_pdb_format( "carbohydrates:glycam_pdb_format" );  }
namespace rings { BooleanOptionKey const rings( "rings" );  }
namespace rings { BooleanOptionKey const lock_rings( "rings:lock_rings" );  }
namespace rings { BooleanOptionKey const idealize_rings( "rings:idealize_rings" );  }
namespace rings { BooleanOptionKey const sample_high_energy_conformers( "rings:sample_high_energy_conformers" );  }
namespace chemical { BooleanOptionKey const chemical( "chemical" );  }
namespace chemical { StringVectorOptionKey const exclude_patches( "chemical:exclude_patches" );  }
namespace chemical { StringVectorOptionKey const include_patches( "chemical:include_patches" );  }
namespace chemical { StringVectorOptionKey const add_atom_type_set_parameters( "chemical:add_atom_type_set_parameters" );  }
namespace chemical { StringVectorOptionKey const set_atom_properties( "chemical:set_atom_properties" );  }
namespace chemical { StringVectorOptionKey const patch_selectors( "chemical:patch_selectors" );  }
namespace chemical { BooleanOptionKey const override_rsd_type_limit( "chemical:override_rsd_type_limit" );  }
namespace chemical { StringVectorOptionKey const clone_atom_types( "chemical:clone_atom_types" );  }
namespace chemical { StringVectorOptionKey const reassign_atom_types( "chemical:reassign_atom_types" );  }
namespace chemical { StringVectorOptionKey const reassign_icoor( "chemical:reassign_icoor" );  }
namespace chemical { StringVectorOptionKey const set_atomic_charge( "chemical:set_atomic_charge" );  }
namespace chemical { BooleanOptionKey const enlarge_H_lj( "chemical:enlarge_H_lj" );  }
namespace chemical { BooleanOptionKey const no_hbonds_to_ether_oxygens( "chemical:no_hbonds_to_ether_oxygens" );  }
namespace chemical { BooleanOptionKey const check_rsd_type_finder( "chemical:check_rsd_type_finder" );  }
namespace constraints { BooleanOptionKey const constraints( "constraints" );  }
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
