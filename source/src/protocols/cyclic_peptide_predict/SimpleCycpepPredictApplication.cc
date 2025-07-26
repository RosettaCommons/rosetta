// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication.cc
/// @brief Application-level code for the simple_cycpep_predict app.
/// @details  This application predicts structures of simple backbone-cyclized peptides made of alpha-, beta-, or gamma-amino acids (of any chirality)
/// using generalized kinematic closure (GenKIC) for cyclization, and enforcing user-defined requiresments for numbers of mainchain hydrogen bonds.
/// In April 2017, I also added support for cyclization through disulfide bonds.
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu)

#ifdef BOINC
#include <utility/io/izstream.hh>
#include <utility/boinc/boinc_util.hh>
#include <protocols/boinc/boinc.hh>
#endif // BOINC

// Unit Headers
#include <protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication.hh>
#include <protocols/cyclic_peptide_predict/HierarchicalHybridJD_JobResultsSummary.hh>
#include <protocols/cyclic_peptide_predict/util.hh>

// Package Headers
#include <basic/options/option.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/palette/CustomBaseTypePackerPalette.hh>
#include <core/select/residue_selector/PhiSelector.hh>
#include <core/select/residue_selector/BinSelector.hh>
#include <core/select/residue_selector/NotResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <utility/exit.hh>
#include <basic/Tracer.hh>
#include <core/pose/init_id_map.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/RamaPrePro.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.hh>
#include <protocols/aa_composition/AddCompositionConstraintMover.hh>
#include <protocols/aa_composition/ClearCompositionConstraintsMover.hh>
#include <protocols/cyclic_peptide/OversaturatedHbondAcceptorFilter.hh>
#include <protocols/cyclic_peptide/CycpepSymmetryFilter.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/denovo_design/movers/FastDesign.hh>
#include <protocols/cyclic_peptide/PeptideInternalHbondsFilter.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/scoring/Energies.hh>
#include <utility/file/file_sys_util.hh>

//Constraints
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/func/CharmmPeriodicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>

//Disulfides
#include <protocols/cyclic_peptide/TryDisulfPermutations.hh>
#include <protocols/score_filters/ScoreTypeFilter.hh>

//N-methylation
#include <core/select/residue_selector/ResidueIndexSelector.fwd.hh>
#include <protocols/simple_moves/ModifyVariantTypeMover.hh>

//TBMB
#include <protocols/cyclic_peptide/CrosslinkerMover.hh>
#include <protocols/cyclic_peptide/crosslinker/TBMB_Helper.hh>
#include <protocols/cyclic_peptide/crosslinker/TMA_Helper.hh>
#include <protocols/cyclic_peptide/crosslinker/1_4_BBMB_Helper.hh>
#include <protocols/cyclic_peptide/crosslinker/TrigonalPyramidalMetal_Helper.hh>
#include <protocols/cyclic_peptide/crosslinker/TrigonalPlanarMetal_Helper.hh>
#include <protocols/cyclic_peptide/crosslinker/SquarePyramidalMetal_Helper.hh>
#include <protocols/cyclic_peptide/crosslinker/SquarePlanarMetal_Helper.hh>
#include <protocols/cyclic_peptide/crosslinker/TetrahedralMetal_Helper.hh>
#include <protocols/cyclic_peptide/crosslinker/OctahedralMetal_Helper.hh>

//Thioether cyclization
#include <protocols/cyclic_peptide/crosslinker/thioether_util.hh>

//Lanthipeptide cyclization
#include <protocols/cyclic_peptide/crosslinker/lanthionine_util.hh>

// Project Headers
#include <protocols/cyclic_peptide/PeptideStubMover.hh>
#include <protocols/simple_moves/DeclareBond.hh>
#include <protocols/generalized_kinematic_closure/GeneralizedKIC.hh>

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/cyclic_peptide.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

//numeric headers
#include <numeric/random/random.hh>
#include <numeric/constants.hh>
#include <numeric/angle.functions.hh>

// Utility headers
#include <basic/Tracer.hh>

// C++ headers
#include <cstdio>

#ifdef GL_GRAPHICS
// for graphics
#include <protocols/viewer/viewers.hh>
#endif

#include <utility/fixedsizearray1.hh> // AUTO IWYU For fixedsizearray1
#include <utility/file/file_sys_util.hh>

static basic::Tracer TR( "protocols.cyclic_peptide_predict.SimpleCycpepPredictApplication" );

/// @brief Register the set of options that this application uses (for the help menu).
///
void
protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::sequence_file                        );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::cyclization_type                     );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::use_chainbreak_energy                  );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::genkic_closure_attempts              );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::genkic_min_solution_count            );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::cyclic_permutations                  );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::exclude_residues_from_rms            );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::use_rama_filter                      );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::rama_cutoff                          );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::high_hbond_weight_multiplier         );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::min_genkic_hbonds                    );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::min_final_hbonds                     );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::total_energy_cutoff                  );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::hbond_energy_cutoff                  );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::fast_relax_rounds                    );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::count_sc_hbonds                      );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::checkpoint_job_identifier            );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::default_rama_sampling_table          );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::rama_sampling_table_by_res           );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::checkpoint_file                      );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::rand_checkpoint_file                 );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::require_disulfides                   );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::disulf_cutoff_prerelax               );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::disulf_cutoff_postrelax              );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::user_set_alpha_dihedrals             );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::user_set_alpha_dihedral_perturbation );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::bondlength_perturbation_magnitude    );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::bondangle_perturbation_magnitude     );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::filter_oversaturated_hbond_acceptors );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::hbond_acceptor_energy_cutoff         );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::design_peptide                       );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::allowed_residues_by_position         );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::prohibit_D_at_negative_phi           );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::prohibit_L_at_positive_phi           );
	option.add_relevant( basic::options::OptionKeys::score::aa_composition_setup_file                     );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::L_alpha_comp_file                    );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::D_alpha_comp_file                    );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::L_beta_comp_file                     );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::D_beta_comp_file                     );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::do_not_count_adjacent_res_hbonds     );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::sample_cis_pro_frequency             );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::angle_relax_rounds                   );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::angle_length_relax_rounds            );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::cartesian_relax_rounds               );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::use_classic_rama_for_sampling        );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::n_methyl_positions                   );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::lariat_sidechain_index               );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::lariat_sample_cis               );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::sidechain_isopeptide_indices         );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::lanthionine_positions                       );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::paraBBMB_positions                       );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::use_paraBBMB_filters                     );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::paraBBMB_sidechain_distance_filter_multiplier  );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::paraBBMB_constraints_energy_filter_multiplier  );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::link_all_cys_with_paraBBMB               );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::TBMB_positions                       );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::use_TBMB_filters                     );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::TBMB_sidechain_distance_filter_multiplier  );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::TBMB_constraints_energy_filter_multiplier  );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::link_all_cys_with_TBMB               );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::TMA_positions                        );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::use_TMA_filters                      );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::TMA_sidechain_distance_filter_multiplier   );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::TMA_constraints_energy_filter_multiplier   );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::trigonal_pyramidal_metal_positions          );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::use_trigonal_pyramidal_metal_filters        );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::trigonal_pyramidal_metal_sidechain_distance_filter_multiplier   );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::trigonal_pyramidal_metal_constraints_energy_filter_multiplier   );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::trigonal_planar_metal_positions          );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::use_trigonal_planar_metal_filters        );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::trigonal_planar_metal_sidechain_distance_filter_multiplier   );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::trigonal_planar_metal_constraints_energy_filter_multiplier   );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::square_pyramidal_metal_positions          );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::use_square_pyramidal_metal_filters        );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::square_pyramidal_metal_sidechain_distance_filter_multiplier   );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::square_pyramidal_metal_constraints_energy_filter_multiplier   );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::square_planar_metal_positions          );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::use_square_planar_metal_filters        );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::square_planar_metal_sidechain_distance_filter_multiplier   );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::square_planar_metal_constraints_energy_filter_multiplier   );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::tetrahedral_metal_positions          );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::use_tetrahedral_metal_filters        );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::tetrahedral_metal_sidechain_distance_filter_multiplier   );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::tetrahedral_metal_constraints_energy_filter_multiplier   );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::octahedral_metal_positions          );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::use_octahedral_metal_filters        );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::octahedral_metal_sidechain_distance_filter_multiplier   );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::octahedral_metal_constraints_energy_filter_multiplier   );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::require_symmetry_repeats             );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::require_symmetry_mirroring           );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::require_symmetry_angle_threshold     );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::require_symmetry_perturbation        );
	option.add_relevant( basic::options::OptionKeys::out::path::pdb);
	option.add_relevant( basic::options::OptionKeys::out::path::all);
	option.add_relevant( basic::options::OptionKeys::out::path::path);
	option.add_relevant( basic::options::OptionKeys::out::prefix);
	option.add_relevant( basic::options::OptionKeys::out::suffix);

#ifdef USEMPI //Options that are only needed in the MPI version:
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::MPI_auto_2level_distribution         );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::MPI_processes_by_level               );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::MPI_batchsize_by_level               );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::MPI_sort_by                          );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::MPI_choose_highest                   );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::MPI_output_fraction                  );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::MPI_stop_after_time                  );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::MPI_pnear_lambda                     );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::MPI_pnear_kbt                        );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::compute_rmsd_to_lowest               );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::compute_pnear_to_this_fract          );
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::compute_ensemble_sasa_metrics        );
#ifdef MULTI_THREADED //Options that are only needed in the MPI+threads version:
	option.add_relevant( basic::options::OptionKeys::cyclic_peptide::threads_per_worker                    );
#endif //ifdef MULTI_THREADED
#endif //ifdef USEMPI
	return;
}

/// @brief Given a cyclization type enum, return its name string.
///
protocols::cyclic_peptide_predict::SCPA_cyclization_type
protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::get_cyclization_type_from_name(
	std::string const &name
) {
	using namespace protocols::cyclic_peptide_predict;
	for ( core::Size i(1); i<SCPA_number_of_types; ++i ) {
		if ( name == get_cyclization_name_from_type( static_cast<SCPA_cyclization_type>(i) ) ) {
			return static_cast<SCPA_cyclization_type>(i);
		}
	}
	return SCPA_invalid_type;
}

/// @brief Given a cyclization name string, return its type enum.
///
std::string
protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::get_cyclization_name_from_type(
	protocols::cyclic_peptide_predict::SCPA_cyclization_type const type
) {
	using namespace protocols::cyclic_peptide_predict;
	switch( type ) {
	case SCPA_n_to_c_amide_bond :
		return "n_to_c_amide_bond";
	case SCPA_terminal_disulfide :
		return "terminal_disulfide";
	case SCPA_thioether_lariat :
		return "thioether_lariat";
	case SCPA_lanthipeptide :
		return "lanthipeptide";
	case SCPA_nterm_isopeptide_lariat :
		return "nterm_isopeptide_lariat";
	case SCPA_cterm_isopeptide_lariat :
		return "cterm_isopeptide_lariat";
	case SCPA_sidechain_isopeptide :
		return "sidechain_isopeptide";
	case SCPA_invalid_type :
		return "INVALID";
	}
	return "INVALID";
}


namespace protocols {
namespace cyclic_peptide_predict {


/// @brief Constructor
/// @details If allow_file_read is true, initialization triggers reads
/// from the filesystem.
SimpleCycpepPredictApplication::SimpleCycpepPredictApplication(
	bool const allow_file_read
) :
	my_rank_(0),
	already_completed_job_count_(0),
	cyclization_type_( SCPA_n_to_c_amide_bond ),
	use_chainbreak_energy_(true),
	scorefxn_(),
	suppress_checkpoints_(false),
	silent_out_(false),
	silentlist_out_(false),
	silentlist_(nullptr),
	summarylist_(nullptr),
	native_pose_(),
	out_filename_("S_"),
	out_scorefilename_("default.sc"),
	sequence_file_(""),
	sequence_string_(""),
	sequence_length_(0),
	genkic_closure_attempts_(1),
	genkic_min_solution_count_(1),
	cyclic_permutations_(true),
	use_rama_filter_(true),
	rama_cutoff_(0.3),
	high_hbond_weight_multiplier_(10),
	min_genkic_hbonds_(3.0),
	min_final_hbonds_(0.0),
	total_energy_cutoff_(0.0),
	use_total_energy_cutoff_(false),
	hbond_energy_cutoff_(-0.25),
	fast_relax_rounds_(3),
	count_sc_hbonds_(false),
	native_exists_(false),
	native_filename_(""),
	nstruct_(1),
	checkpoint_job_identifier_(""),
	checkpoint_filename_("checkpoint.txt"),
	default_rama_table_type_( core::scoring::unknown_ramatable_type ),
	rama_table_type_by_res_(),
	rand_checkpoint_file_("rng.state.gz"),
	try_all_disulfides_(false),
	disulf_energy_cutoff_prerelax_(15.0),
	disulf_energy_cutoff_postrelax_(0.5),
	user_set_alpha_dihedrals_(),
	user_set_dihedral_perturbation_(0.0),
	bondlength_perturbation_magnitude_(0.0),
	bondangle_perturbation_magnitude_(0.0),
	filter_oversaturated_hbond_acceptors_(true),
	oversaturated_hbond_cutoff_energy_(-0.1),
	sample_cis_pro_(true),
	sample_cis_pro_frequency_(0.3),
	design_peptide_(false),
	design_filename_(""),
	prevent_design_file_read_( !allow_file_read ),
	allowed_canonicals_by_position_(),
	allowed_noncanonicals_by_position_(),
	prohibit_D_at_negative_phi_(true),
	prohibit_L_at_positive_phi_(true),
	use_aa_comp_(false),
	L_alpha_comp_file_exists_(false),
	D_alpha_comp_file_exists_(false),
	L_beta_comp_file_exists_(false),
	D_beta_comp_file_exists_(false),
	comp_file_contents_L_alpha_(""),
	comp_file_contents_D_alpha_(""),
	comp_file_contents_L_beta_(""),
	comp_file_contents_D_beta_(""),
	abba_bins_binfile_("ABBA.bin_params"),
	do_not_count_adjacent_res_hbonds_(true),
	angle_relax_rounds_(0),
	angle_length_relax_rounds_(0),
	cartesian_relax_rounds_(0),
	use_rama_prepro_for_sampling_(true),
	n_methyl_positions_(),
	lanthionine_positions_(),
	parabbmb_positions_(),
	link_all_cys_with_parabbmb_(false),
	use_parabbmb_filters_(true),
	parabbmb_sidechain_distance_filter_multiplier_(1.0),
	parabbmb_constraints_energy_filter_multiplier_(1.0),
	tbmb_positions_(),
	link_all_cys_with_tbmb_(false),
	use_tbmb_filters_(true),
	tbmb_sidechain_distance_filter_multiplier_(1.0),
	tbmb_constraints_energy_filter_multiplier_(1.0),
	tma_positions_(),
	use_tma_filters_(true),
	tma_sidechain_distance_filter_multiplier_(1.0),
	tma_constraints_energy_filter_multiplier_(1.0),
	trigonal_pyramidal_metal_positions_(),
	use_trigonal_pyramidal_metal_filters_(true),
	trigonal_pyramidal_metal_sidechain_distance_filter_multiplier_(1.0),
	trigonal_pyramidal_metal_constraints_energy_filter_multiplier_(1.0),
	trigonal_planar_metal_positions_(),
	use_trigonal_planar_metal_filters_(true),
	trigonal_planar_metal_sidechain_distance_filter_multiplier_(1.0),
	trigonal_planar_metal_constraints_energy_filter_multiplier_(1.0),
	square_pyramidal_metal_positions_(),
	use_square_pyramidal_metal_filters_(true),
	square_pyramidal_metal_sidechain_distance_filter_multiplier_(1.0),
	square_pyramidal_metal_constraints_energy_filter_multiplier_(1.0),
	square_planar_metal_positions_(),
	use_square_planar_metal_filters_(true),
	square_planar_metal_sidechain_distance_filter_multiplier_(1.0),
	square_planar_metal_constraints_energy_filter_multiplier_(1.0),
	tetrahedral_metal_positions_(),
	use_tetrahedral_metal_filters_(true),
	tetrahedral_metal_sidechain_distance_filter_multiplier_(1.0),
	tetrahedral_metal_constraints_energy_filter_multiplier_(1.0),
	octahedral_metal_positions_(),
	use_octahedral_metal_filters_(true),
	octahedral_metal_sidechain_distance_filter_multiplier_(1.0),
	octahedral_metal_constraints_energy_filter_multiplier_(1.0),
	required_symmetry_repeats_(1),
	required_symmetry_mirroring_(false),
	required_symmetry_angle_threshold_(10.0),
	required_symmetry_perturbation_(0.0),
	exclude_residues_from_rms_(),
	lariat_sidechain_index_(0),
	lariat_sample_cis_(true),
	sidechain_isopeptide_indices_(std::pair<core::Size, core::Size>(0,0))
	//TODO -- initialize variables here.
{
	register_with_citation_manager(); //Defined in util.hh/util.cc.
	initialize_from_options();
}


/// @brief Explicit virtual destructor.
///
SimpleCycpepPredictApplication::~SimpleCycpepPredictApplication() = default;


/// @brief Explicit copy constructor.
///
SimpleCycpepPredictApplication::SimpleCycpepPredictApplication( SimpleCycpepPredictApplication const &src ) :
	VirtualBase( src ),
	my_rank_(src.my_rank_),
	already_completed_job_count_( src.already_completed_job_count_ ),
	cyclization_type_( src.cyclization_type_ ),
	use_chainbreak_energy_( src.use_chainbreak_energy_ ),
	scorefxn_(), //Cloned below
	suppress_checkpoints_(src.suppress_checkpoints_),
	silent_out_(src.silent_out_),
	silentlist_out_(src.silentlist_out_),
	silentlist_(src.silentlist_),
	summarylist_(src.summarylist_),
	native_pose_(src.native_pose_),
	out_filename_(src.out_filename_),
	out_scorefilename_(src.out_scorefilename_),
	sequence_file_(src.sequence_file_),
	sequence_string_(src.sequence_string_),
	sequence_length_(src.sequence_length_),
	genkic_closure_attempts_(src.genkic_closure_attempts_),
	genkic_min_solution_count_(src.genkic_min_solution_count_),
	cyclic_permutations_(src.cyclic_permutations_),
	use_rama_filter_(src.use_rama_filter_),
	rama_cutoff_(src.rama_cutoff_),
	high_hbond_weight_multiplier_(src.high_hbond_weight_multiplier_),
	min_genkic_hbonds_(src.min_genkic_hbonds_),
	min_final_hbonds_(src.min_final_hbonds_),
	total_energy_cutoff_(src.total_energy_cutoff_),
	use_total_energy_cutoff_(src.use_total_energy_cutoff_),
	hbond_energy_cutoff_(src.hbond_energy_cutoff_),
	fast_relax_rounds_(src.fast_relax_rounds_),
	count_sc_hbonds_(src.count_sc_hbonds_),
	native_exists_(src.native_exists_),
	native_filename_(src.native_filename_),
	nstruct_(src.nstruct_),
	checkpoint_job_identifier_(src.checkpoint_job_identifier_),
	checkpoint_filename_(src.checkpoint_filename_),
	default_rama_table_type_( src.default_rama_table_type_ ),
	rama_table_type_by_res_( src.rama_table_type_by_res_ ),
	rand_checkpoint_file_(src.rand_checkpoint_file_),
	try_all_disulfides_(src.try_all_disulfides_),
	disulf_energy_cutoff_prerelax_(src.disulf_energy_cutoff_prerelax_),
	disulf_energy_cutoff_postrelax_(src.disulf_energy_cutoff_postrelax_),
	user_set_alpha_dihedrals_(src.user_set_alpha_dihedrals_),
	user_set_dihedral_perturbation_(src.user_set_dihedral_perturbation_),
	bondlength_perturbation_magnitude_(src.bondlength_perturbation_magnitude_),
	bondangle_perturbation_magnitude_(src.bondangle_perturbation_magnitude_),
	filter_oversaturated_hbond_acceptors_(src.filter_oversaturated_hbond_acceptors_),
	oversaturated_hbond_cutoff_energy_(src.oversaturated_hbond_cutoff_energy_),
	sample_cis_pro_(src.sample_cis_pro_),
	sample_cis_pro_frequency_(src.sample_cis_pro_frequency_),
	design_peptide_(src.design_peptide_),
	design_filename_(src.design_filename_),
	prevent_design_file_read_(src.prevent_design_file_read_),
	allowed_canonicals_by_position_(src.allowed_canonicals_by_position_),
	allowed_noncanonicals_by_position_(src.allowed_noncanonicals_by_position_),
	prohibit_D_at_negative_phi_(src.prohibit_D_at_negative_phi_),
	prohibit_L_at_positive_phi_(src.prohibit_L_at_positive_phi_),
	use_aa_comp_(src.use_aa_comp_),
	L_alpha_comp_file_exists_(src.L_alpha_comp_file_exists_),
	D_alpha_comp_file_exists_(src.D_alpha_comp_file_exists_),
	L_beta_comp_file_exists_(src.L_beta_comp_file_exists_),
	D_beta_comp_file_exists_(src.D_beta_comp_file_exists_),
	comp_file_contents_L_alpha_(src.comp_file_contents_L_alpha_),
	comp_file_contents_D_alpha_(src.comp_file_contents_D_alpha_),
	comp_file_contents_L_beta_(src.comp_file_contents_L_beta_),
	comp_file_contents_D_beta_(src.comp_file_contents_D_beta_),
	abba_bins_binfile_(src.abba_bins_binfile_),
	do_not_count_adjacent_res_hbonds_(src.do_not_count_adjacent_res_hbonds_),
	angle_relax_rounds_(src.angle_relax_rounds_),
	angle_length_relax_rounds_(src.angle_length_relax_rounds_),
	cartesian_relax_rounds_(src.cartesian_relax_rounds_),
	use_rama_prepro_for_sampling_(src.use_rama_prepro_for_sampling_),
	n_methyl_positions_(src.n_methyl_positions_),
	lanthionine_positions_(src.lanthionine_positions_),
	parabbmb_positions_(src.parabbmb_positions_),
	link_all_cys_with_parabbmb_(src.link_all_cys_with_parabbmb_),
	use_parabbmb_filters_(src.use_parabbmb_filters_),
	parabbmb_sidechain_distance_filter_multiplier_(src.parabbmb_sidechain_distance_filter_multiplier_),
	parabbmb_constraints_energy_filter_multiplier_(src.parabbmb_constraints_energy_filter_multiplier_),
	tbmb_positions_(src.tbmb_positions_),
	link_all_cys_with_tbmb_(src.link_all_cys_with_tbmb_),
	use_tbmb_filters_(src.use_tbmb_filters_),
	tbmb_sidechain_distance_filter_multiplier_(src.tbmb_sidechain_distance_filter_multiplier_),
	tbmb_constraints_energy_filter_multiplier_(src.tbmb_constraints_energy_filter_multiplier_),
	tma_positions_(src.tma_positions_),
	use_tma_filters_(src.use_tma_filters_),
	tma_sidechain_distance_filter_multiplier_(src.tma_sidechain_distance_filter_multiplier_),
	tma_constraints_energy_filter_multiplier_(src.tma_constraints_energy_filter_multiplier_),
	trigonal_pyramidal_metal_positions_(src.trigonal_pyramidal_metal_positions_),
	use_trigonal_pyramidal_metal_filters_(src.use_trigonal_pyramidal_metal_filters_),
	trigonal_pyramidal_metal_sidechain_distance_filter_multiplier_(src.trigonal_pyramidal_metal_sidechain_distance_filter_multiplier_),
	trigonal_pyramidal_metal_constraints_energy_filter_multiplier_(src.trigonal_pyramidal_metal_constraints_energy_filter_multiplier_),
	trigonal_planar_metal_positions_(src.trigonal_planar_metal_positions_),
	use_trigonal_planar_metal_filters_(src.use_trigonal_planar_metal_filters_),
	trigonal_planar_metal_sidechain_distance_filter_multiplier_(src.trigonal_planar_metal_sidechain_distance_filter_multiplier_),
	trigonal_planar_metal_constraints_energy_filter_multiplier_(src.trigonal_planar_metal_constraints_energy_filter_multiplier_),
	square_pyramidal_metal_positions_(src.square_pyramidal_metal_positions_),
	use_square_pyramidal_metal_filters_(src.use_square_pyramidal_metal_filters_),
	square_pyramidal_metal_sidechain_distance_filter_multiplier_(src.square_pyramidal_metal_sidechain_distance_filter_multiplier_),
	square_pyramidal_metal_constraints_energy_filter_multiplier_(src.square_pyramidal_metal_constraints_energy_filter_multiplier_),
	square_planar_metal_positions_(src.square_planar_metal_positions_),
	use_square_planar_metal_filters_(src.use_square_planar_metal_filters_),
	square_planar_metal_sidechain_distance_filter_multiplier_(src.square_planar_metal_sidechain_distance_filter_multiplier_),
	square_planar_metal_constraints_energy_filter_multiplier_(src.square_planar_metal_constraints_energy_filter_multiplier_),
	tetrahedral_metal_positions_(src.tetrahedral_metal_positions_),
	use_tetrahedral_metal_filters_(src.use_tetrahedral_metal_filters_),
	tetrahedral_metal_sidechain_distance_filter_multiplier_(src.tetrahedral_metal_sidechain_distance_filter_multiplier_),
	tetrahedral_metal_constraints_energy_filter_multiplier_(src.tetrahedral_metal_constraints_energy_filter_multiplier_),
	octahedral_metal_positions_(src.octahedral_metal_positions_),
	use_octahedral_metal_filters_(src.use_octahedral_metal_filters_),
	octahedral_metal_sidechain_distance_filter_multiplier_(src.octahedral_metal_sidechain_distance_filter_multiplier_),
	octahedral_metal_constraints_energy_filter_multiplier_(src.octahedral_metal_constraints_energy_filter_multiplier_),
	required_symmetry_repeats_(src.required_symmetry_repeats_),
	required_symmetry_mirroring_(src.required_symmetry_mirroring_),
	required_symmetry_angle_threshold_(src.required_symmetry_angle_threshold_),
	required_symmetry_perturbation_(src.required_symmetry_perturbation_),
	exclude_residues_from_rms_(src.exclude_residues_from_rms_),
	lariat_sidechain_index_(src.lariat_sidechain_index_),
	lariat_sample_cis_(src.lariat_sample_cis_),
	sidechain_isopeptide_indices_(src.sidechain_isopeptide_indices_),
	out_prefix_(src.out_prefix_),
	out_suffix_(src.out_suffix_),
	out_path_(src.out_path_)

	//TODO -- copy variables here.
{
	if ( src.scorefxn_ ) scorefxn_ = (src.scorefxn_)->clone();
}

/// @brief Initialize the application.
/// @details Initializes using the option system.
void
SimpleCycpepPredictApplication::initialize_from_options(
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace utility::io;

	//Initial checks:
	runtime_assert_string_msg( option[basic::options::OptionKeys::cyclic_peptide::sequence_file].user(), "Error in simple_cycpep_predict app: the user MUST provide a sequence file using the \"-cyclic_peptide:sequence_file\" flag." );
	runtime_assert_string_msg( option[basic::options::OptionKeys::cyclic_peptide::genkic_closure_attempts]() >= 0, "Error in simple_cycpep_predict app: the number of GeneralizedKIC closure attempts (\"-cyclic_peptide:genkic_closure_attempts\" flag) cannot be negative.  (Note also that setting this to zero is risky, since GenKIC will continue to seek solutions until the minimum number of solutions is reached.)" );
	runtime_assert_string_msg( option[basic::options::OptionKeys::cyclic_peptide::genkic_min_solution_count]() >= 0, "Error in simple_cycpep_predict app: the minimum number of GenKIC solutions (\"-cyclic_peptide:genkic_min_solution_count\" flag) cannot be negative.  (Note also that setting this to zero means no minimum.)" );
	runtime_assert_string_msg( !(option[basic::options::OptionKeys::cyclic_peptide::genkic_closure_attempts]() == 0 && option[basic::options::OptionKeys::cyclic_peptide::genkic_min_solution_count]() == 0), "Error in simple_cycpep_predict app: both the \"-cyclic_peptide:genkic_closure_attempts\" and \"-cyclic_peptide:genkic_min_solution_count\" flags were set to zero.  This would result in GenKIC looping infinitely." );
	runtime_assert_string_msg( option[basic::options::OptionKeys::cyclic_peptide::min_genkic_hbonds]() >= 0.0, "Error in simple_cycpep_predict app: the minimum number of hbonds during GenKIC steps (\"-cyclic_peptide:min_genkic_hbonds\" flag) can be zero, but cannot be negative." );
	runtime_assert_string_msg( option[basic::options::OptionKeys::cyclic_peptide::min_final_hbonds]() >= 0.0, "Error in simple_cycpep_predict app: the minimum number of hbonds after relaxation steps (\"-cyclic_peptide:min_final_hbonds\" flag) can be zero, but cannot be negative." );
	runtime_assert_string_msg( !( option[out::file::silent].user() && option[out::file::o].user() ), "Error in simple_cycpep_predict app: either silent file output (\"-out:file:silent\" flag) or PDB output (\"-out:file:o\") output may be used, but not both." );

	//Set cyclization type:
	SCPA_cyclization_type const type_from_options( get_cyclization_type_from_name( option[basic::options::OptionKeys::cyclic_peptide::cyclization_type]() ) );
	runtime_assert_string_msg( type_from_options != SCPA_invalid_type, "Error in simple_cycpep_predict app: the provided cyclization type is not an allowed type!" );
	set_cyclization_type( type_from_options );
	runtime_assert_string_msg(
		!( (option[basic::options::OptionKeys::cyclic_peptide::link_all_cys_with_TBMB]() || option[basic::options::OptionKeys::cyclic_peptide::link_all_cys_with_paraBBMB]() )
		&& ( cyclization_type() == SCPA_terminal_disulfide || cyclization_type() == SCPA_thioether_lariat || cyclization_type() == SCPA_lanthipeptide ) ),
		"Error in simple_cycpep_predict app: linking all cysteine residues with a crosslinker (\"-cyclic_peptide:link_all_cys_with_TBMB\" or \"-cyclic_peptide:link_all_cys_with_paraBBMB\" flags) is incompatible with teriminal disulfide cyclization, thioether, and lanthipeptide (\"-cyclic_peptide:cyclization_type\" flag)."
	);
	runtime_assert_string_msg(
		!( option[basic::options::OptionKeys::cyclic_peptide::require_symmetry_repeats]() > 2 && (cyclization_type() == SCPA_terminal_disulfide || cyclization_type() == SCPA_thioether_lariat
		|| cyclization_type() == SCPA_nterm_isopeptide_lariat || cyclization_type() == SCPA_cterm_isopeptide_lariat
		|| cyclization_type() == SCPA_sidechain_isopeptide || cyclization_type() == SCPA_lanthipeptide ) ),
		"Error in simple_cycpep_predict app: quasi-symmetric sampling (\"-cyclic_peptide:require_symmetry_repeats\" flag) is incompatible with teriminal disulfide cyclization, thioether lariat cyclization, lanthipeptides, or with isopeptide bond cyclization (\"-cyclic_peptide:cyclization_type\" flag)."
	);
	runtime_assert_string_msg(
		!(option[basic::options::OptionKeys::cyclic_peptide::link_all_cys_with_TBMB]() && option[basic::options::OptionKeys::cyclic_peptide::link_all_cys_with_paraBBMB]() ),
		"Error in simple_cycpep_predict app: The \"-link_all_cys_with_TBMB\" and \"-link_all_cys_with_paraBBMB\" options are mutually incompatible."
	);

	//Set whether to use chainbreak energy:
	set_use_chainbreak_energy( option[basic::options::OptionKeys::cyclic_peptide::use_chainbreak_energy]() );

	//Copy options to private member variables:
	sequence_file_ = option[basic::options::OptionKeys::cyclic_peptide::sequence_file]();
	genkic_closure_attempts_ = static_cast<core::Size>( option[basic::options::OptionKeys::cyclic_peptide::genkic_closure_attempts]() );
	genkic_min_solution_count_ = static_cast<core::Size>( option[basic::options::OptionKeys::cyclic_peptide::genkic_min_solution_count]() );

	if ( cyclization_type() == SCPA_terminal_disulfide || cyclization_type() == SCPA_thioether_lariat || cyclization_type() == SCPA_lanthipeptide || cyclization_type() == SCPA_cterm_isopeptide_lariat || cyclization_type() == SCPA_nterm_isopeptide_lariat || cyclization_type() == SCPA_sidechain_isopeptide ) {
		cyclic_permutations_ = false;
		if ( option[basic::options::OptionKeys::cyclic_peptide::cyclic_permutations]() && TR.Warning.visible() ) {
			TR.Warning << "Warning: the \"-cyclic_peptide:cyclic_permutations\" option was set to \"true\", but this is incompatible with the chosen cyclization mode (\"-cyclic_peptide:cyclization_type\").  Disabling cyclic permutations." << std::endl;
		}
	} else {
		cyclic_permutations_ = option[basic::options::OptionKeys::cyclic_peptide::cyclic_permutations]();
	}

	use_rama_filter_ = option[basic::options::OptionKeys::cyclic_peptide::use_rama_filter]();
	rama_cutoff_ = static_cast<core::Real>( option[basic::options::OptionKeys::cyclic_peptide::rama_cutoff]() );
	high_hbond_weight_multiplier_ = static_cast<core::Real>( option[basic::options::OptionKeys::cyclic_peptide::high_hbond_weight_multiplier]() );
	min_genkic_hbonds_ = static_cast<core::Real>( option[basic::options::OptionKeys::cyclic_peptide::min_genkic_hbonds]() );
	min_final_hbonds_ = static_cast<core::Real>( option[basic::options::OptionKeys::cyclic_peptide::min_final_hbonds]() );
	if ( option[basic::options::OptionKeys::cyclic_peptide::total_energy_cutoff].user() ) {
		set_total_energy_cutoff( option[basic::options::OptionKeys::cyclic_peptide::total_energy_cutoff]() );
	}
	hbond_energy_cutoff_ = static_cast<core::Real>( option[basic::options::OptionKeys::cyclic_peptide::hbond_energy_cutoff]() );
	fast_relax_rounds_ = static_cast<core::Size>( option[basic::options::OptionKeys::cyclic_peptide::fast_relax_rounds]() );
	count_sc_hbonds_ = option[basic::options::OptionKeys::cyclic_peptide::count_sc_hbonds]();
	try_all_disulfides_ = option[basic::options::OptionKeys::cyclic_peptide::require_disulfides].value();
	disulf_energy_cutoff_prerelax_ = static_cast<core::Real>( option[basic::options::OptionKeys::cyclic_peptide::disulf_cutoff_prerelax]() );
	disulf_energy_cutoff_postrelax_ = static_cast<core::Real>( option[basic::options::OptionKeys::cyclic_peptide::disulf_cutoff_postrelax]() );

	//Get the scorefunction:
	if ( !prevent_design_file_read_ ) scorefxn_ = core::scoring::get_score_function(); //Reads from file.  Don't use in MPI mode.

	//Read in the comp files, if any (in non-MPI mode):
	if ( !prevent_design_file_read_ ) {
		if ( option[basic::options::OptionKeys::cyclic_peptide::L_alpha_comp_file].user() ) {
			read_file_into_string( comp_file_contents_L_alpha_, option[basic::options::OptionKeys::cyclic_peptide::L_alpha_comp_file](), false );
			L_alpha_comp_file_exists_=true;
			TR << "Loaded " << option[basic::options::OptionKeys::cyclic_peptide::L_alpha_comp_file]() << "." << std::endl;
		}
		if ( option[basic::options::OptionKeys::cyclic_peptide::D_alpha_comp_file].user() ) {
			read_file_into_string( comp_file_contents_D_alpha_, option[basic::options::OptionKeys::cyclic_peptide::D_alpha_comp_file](), false );
			D_alpha_comp_file_exists_=true;
			TR << "Loaded " << option[basic::options::OptionKeys::cyclic_peptide::D_alpha_comp_file]() << "." << std::endl;
		}
		if ( option[basic::options::OptionKeys::cyclic_peptide::L_beta_comp_file].user() ) {
			read_file_into_string( comp_file_contents_L_beta_, option[basic::options::OptionKeys::cyclic_peptide::L_beta_comp_file](), false );
			L_beta_comp_file_exists_=true;
			TR << "Loaded " << option[basic::options::OptionKeys::cyclic_peptide::L_beta_comp_file]() << "." << std::endl;
		}
		if ( option[basic::options::OptionKeys::cyclic_peptide::D_beta_comp_file].user() ) {
			read_file_into_string( comp_file_contents_D_beta_, option[basic::options::OptionKeys::cyclic_peptide::D_beta_comp_file](), false );
			D_beta_comp_file_exists_=true;
			TR << "Loaded " << option[basic::options::OptionKeys::cyclic_peptide::D_beta_comp_file]() << "." << std::endl;
		}
	}

	//Get the native, if it exists:
	if ( option[in::file::native].user() ) {
		native_exists_ = true;
		native_filename_ = option[in::file::native]();
	} else {
		native_exists_ = false;
		native_filename_ = "";
	}

	//Set up output file names:
	out_filename_ = "S_";
	out_scorefilename_ = "default.sc";
	if ( option[out::file::silent].user() ) {
		out_filename_=option[out::file::silent]();
		silent_out_=true;
		option[ basic::options::OptionKeys::out::file::silent_struct_type ].def( "binary"); //Force binary file output.
	} else if ( option[out::file::o].user() ) {
		out_filename_=option[out::file::o]();
		silent_out_=false;
	}
	if ( option[out::file::scorefile].user() ) { out_scorefilename_=option[out::file::scorefile](); }

	//Setup PDB suffix/directory/etc.
	if ( option[out::prefix].user() ) {
		out_prefix_ = option[out::prefix]();
	}

	if ( option[ out::suffix].user() ) {
		out_suffix_ = option[out::suffix]();
	}
	if ( option[ out::path::pdb ].user() ) {
		out_path_ = option[ out::path::pdb ]().path();
	} else if ( option[ out::path::all ].user() ) {
		out_path_ =  option[ out::path::all ]().path();
	} else if ( option[ out::path::path ].user() ) {
		out_path_ = option[ out::path::path ]().path();
	}

	//Figure out number of structures to try to generate:
	if ( option[out::nstruct].user() ) {
		runtime_assert_string_msg( option[out::nstruct]() > 0, "Error in simple_cycpep_predict app: the \"-out:nstruct\" flag's value cannot be less than 1." );
		nstruct_ = static_cast<core::Size>(option[out::nstruct]());
	} else { nstruct_ = 1; }

	checkpoint_job_identifier_ = option[basic::options::OptionKeys::cyclic_peptide::checkpoint_job_identifier]();
	checkpoint_filename_ = option[basic::options::OptionKeys::cyclic_peptide::checkpoint_file]();
	rand_checkpoint_file_ = option[basic::options::OptionKeys::cyclic_peptide::rand_checkpoint_file]();

	//Figure out what custom Ramachandran tables we're using, if any:
	set_default_rama_table_type( option[basic::options::OptionKeys::cyclic_peptide::default_rama_sampling_table]() );
	set_rama_table_type_by_res( option[basic::options::OptionKeys::cyclic_peptide::rama_sampling_table_by_res]() );

	//If the user has set particular dihedrals, read them in
	if ( option[basic::options::OptionKeys::cyclic_peptide::user_set_alpha_dihedrals].user() ) {
		utility::vector1 < core::Real > const user_set_dihedrals( option[basic::options::OptionKeys::cyclic_peptide::user_set_alpha_dihedrals]() );
		runtime_assert_string_msg( user_set_dihedrals.size() % 4 == 0, "Error in simple_cycpep_predict app: the \"-user_set_alpha_dihedrals\" option must be followed by one or more groups of four numbers.  Each group must consist of a sequence position, then phi/psi/omega values." );
		for ( core::Size i=1, imax=user_set_dihedrals.size(); i<=imax; i+=4 ) {
			auto const seqpos( static_cast< core::Size >( user_set_dihedrals[i] ) );
			runtime_assert_string_msg( user_set_alpha_dihedrals_.count(seqpos) == 0, "Error in simple_cycpep_predict app: a residue index was specified more than once with the \"-user_set_alpha_dihedrals\" option." );
			utility::vector1< core::Real > phipsiomegavect;
			phipsiomegavect.reserve(3);
			phipsiomegavect.push_back( user_set_dihedrals[i+1] );
			phipsiomegavect.push_back( user_set_dihedrals[i+2] );
			phipsiomegavect.push_back( user_set_dihedrals[i+3] );
			user_set_alpha_dihedrals_[seqpos] = phipsiomegavect;
		}
	}
	user_set_dihedral_perturbation_ = option[basic::options::OptionKeys::cyclic_peptide::user_set_alpha_dihedral_perturbation]();
	bondlength_perturbation_magnitude_ = option[basic::options::OptionKeys::cyclic_peptide::bondlength_perturbation_magnitude]();
	bondangle_perturbation_magnitude_ = option[basic::options::OptionKeys::cyclic_peptide::bondangle_perturbation_magnitude]();

	filter_oversaturated_hbond_acceptors_ = option[basic::options::OptionKeys::cyclic_peptide::filter_oversaturated_hbond_acceptors]();
	oversaturated_hbond_cutoff_energy_ = option[basic::options::OptionKeys::cyclic_peptide::hbond_acceptor_energy_cutoff]();

	//Options related to sampling cis prolines:
	if ( option[basic::options::OptionKeys::cyclic_peptide::sample_cis_pro_frequency].user() ) { //Turn on cis proline sampling iff the user specifies it.
		set_sample_cis_pro_frequency( option[basic::options::OptionKeys::cyclic_peptide::sample_cis_pro_frequency]() );
	}

	//Options related to design:
	design_peptide_ = option[basic::options::OptionKeys::cyclic_peptide::design_peptide]();
	if ( option[basic::options::OptionKeys::cyclic_peptide::allowed_residues_by_position].user() ) {
		design_filename_ = option[basic::options::OptionKeys::cyclic_peptide::allowed_residues_by_position]();
	} else {
		prevent_design_file_read_ = true;
	}
	prohibit_D_at_negative_phi_ = option[basic::options::OptionKeys::cyclic_peptide::prohibit_D_at_negative_phi]();
	prohibit_L_at_positive_phi_ = option[basic::options::OptionKeys::cyclic_peptide::prohibit_L_at_positive_phi]();

	//Read design options from file:
	if ( !prevent_design_file_read_ ) {
		read_peptide_design_file( design_filename_, allowed_canonicals_by_position_, allowed_noncanonicals_by_position_ );
		prevent_design_file_read_ = true; //Prevents re-read if reinitialized.
	}

	// If we're designing and the user has specified a comp file, turn on aa_composition for design steps.
	if ( option[basic::options::OptionKeys::cyclic_peptide::design_peptide]() &&
			( option[basic::options::OptionKeys::cyclic_peptide::L_alpha_comp_file].user() ||
			option[basic::options::OptionKeys::cyclic_peptide::D_alpha_comp_file].user() ||
			option[basic::options::OptionKeys::cyclic_peptide::L_beta_comp_file].user() ||
			option[basic::options::OptionKeys::cyclic_peptide::D_beta_comp_file].user() ||
			option[basic::options::OptionKeys::score::aa_composition_setup_file].user()
			)
			) {
		TR << "One or more .comp files have been provided.  The app will turn on the aa_composition score term during design steps." << std::endl;
		use_aa_comp_ = true;
	}

	do_not_count_adjacent_res_hbonds_ = option[basic::options::OptionKeys::cyclic_peptide::do_not_count_adjacent_res_hbonds]();

	if ( option[basic::options::OptionKeys::cyclic_peptide::angle_relax_rounds].user() ) {
		set_angle_relax_rounds( static_cast<core::Size>( option[basic::options::OptionKeys::cyclic_peptide::angle_relax_rounds]() ) );
	}
	if ( option[basic::options::OptionKeys::cyclic_peptide::angle_length_relax_rounds].user() ) {
		set_angle_length_relax_rounds( static_cast<core::Size>( option[basic::options::OptionKeys::cyclic_peptide::angle_length_relax_rounds]() ) );
	}
	if ( option[basic::options::OptionKeys::cyclic_peptide::cartesian_relax_rounds].user() ) {
		set_cartesian_relax_rounds( static_cast<core::Size>( option[basic::options::OptionKeys::cyclic_peptide::cartesian_relax_rounds]() ) );
	}

	if ( option[basic::options::OptionKeys::cyclic_peptide::use_classic_rama_for_sampling].user() ) {
		set_use_rama_prepro_for_sampling( !option[basic::options::OptionKeys::cyclic_peptide::use_classic_rama_for_sampling]() );
	}

	//Store the N-methylated positions.
	if ( option[basic::options::OptionKeys::cyclic_peptide::n_methyl_positions].user() ) {
		core::Size const npos( option[basic::options::OptionKeys::cyclic_peptide::n_methyl_positions]().size() );
		n_methyl_positions_.resize( npos );
		for ( core::Size i=1; i<=npos; ++i ) {
			runtime_assert_string_msg( option[basic::options::OptionKeys::cyclic_peptide::n_methyl_positions]()[i] > 0, "Error in simple_cycpep_predict app: The N-methylated positions must all have indices greater than zero." );
			n_methyl_positions_[i] = static_cast< core::Size >( option[basic::options::OptionKeys::cyclic_peptide::n_methyl_positions]()[i] );
		}
	}

	//CWT store lan and melan positions
	if ( option[basic::options::OptionKeys::cyclic_peptide::lanthionine_positions].user() ) {
		core::Size const nlanres(option[basic::options::OptionKeys::cyclic_peptide::lanthionine_positions]().size());
		runtime_assert_string_msg( nlanres > 0, "Error in simple_cycpep_predict application: The \"-cyclic_peptide:lanthionine_positions\" commandline option must be followed by a list of residues to form lanthionine rings." );
		//runtime_assert_string_msg( nlanres % 2 == 0, "Error in simple_cycpep_predict application: The \"-cyclic_peptide:lanthionine_positions\" commandline option must be followed by a list of residues, where the number of residues in the list is a multiple of two.  Pairs of residues must form lanthionine rings." );
		//Only allowing two specified residues at the moment
		//This is due to limitations of the application only being able to handle rings within rings
		//i.e. cannot handle overlapping or adjacent rings.
		runtime_assert_string_msg( nlanres == 2, "Error in simple_cycpep_predict application: The \"-cyclic_peptide:lanthionine_positions\" commandline option must be followed by a list of two residues. Pairs of residues must form lanthionine rings." );
		core::Size count(0);
		lanthionine_positions_.resize(nlanres / 2);
		for ( core::Size i(1), imax(nlanres / 2); i<=imax; ++i ) {
			utility::vector1 <core::Size> innervect(2);
			for ( core::Size j=1; j<=2; ++j ) {
				++count;
				innervect[j] = option[basic::options::OptionKeys::cyclic_peptide::lanthionine_positions]()[count];
			}
			lanthionine_positions_[i] = innervect;
		}
	}

	//Store the paraBBMB positions.
	if ( option[basic::options::OptionKeys::cyclic_peptide::paraBBMB_positions].user() ) {
		runtime_assert_string_msg( !option[basic::options::OptionKeys::cyclic_peptide::link_all_cys_with_TBMB].user(), "Error in simple_cycpep_predict application: The \"-paraBBMB_positions\" flag and the \"-link_all_cys_with_TBMB\" flag cannot be used together." );
		runtime_assert_string_msg( !option[basic::options::OptionKeys::cyclic_peptide::link_all_cys_with_paraBBMB].user(), "Error in simple_cycpep_predict application: The \"-paraBBMB_positions\" flag and the \"-link_all_cys_with_paraBBMB\" flag cannot be used together." );
		core::Size const nparabbmbres(option[basic::options::OptionKeys::cyclic_peptide::paraBBMB_positions]().size());
		runtime_assert_string_msg( nparabbmbres > 0, "Error in simple_cycpep_predict application: The \"-cyclic_peptide:paraBBMB_positions\" commandline option must be followed by a list of residues to link with 1,4-bis(bromomethyl)benzene." );
		runtime_assert_string_msg( nparabbmbres % 2 == 0, "Error in simple_cycpep_predict application: The \"-cyclic_peptide:paraBBMB_positions\" commandline option must be followed by a list of residues, where the number of residues in the list is a multiple of two.  Paris residues will be linked with 1,4-bis(bromomethyl)benzene." );
		core::Size count(0);
		parabbmb_positions_.resize(nparabbmbres / 2);
		for ( core::Size i(1), imax(nparabbmbres / 2); i<=imax; ++i ) {
			utility::vector1 <core::Size> innervect(2);
			for ( core::Size j=1; j<=2; ++j ) {
				++count;
				innervect[j] = option[basic::options::OptionKeys::cyclic_peptide::paraBBMB_positions]()[count];
			}
			parabbmb_positions_[i] = innervect;
		}
	}
	use_parabbmb_filters_ = option[basic::options::OptionKeys::cyclic_peptide::use_paraBBMB_filters]();
	parabbmb_sidechain_distance_filter_multiplier_ = option[basic::options::OptionKeys::cyclic_peptide::paraBBMB_sidechain_distance_filter_multiplier]();
	parabbmb_constraints_energy_filter_multiplier_ = option[basic::options::OptionKeys::cyclic_peptide::paraBBMB_constraints_energy_filter_multiplier]();
	link_all_cys_with_parabbmb_ = option[basic::options::OptionKeys::cyclic_peptide::link_all_cys_with_paraBBMB]();

	//Store the TBMB positions.
	if ( option[basic::options::OptionKeys::cyclic_peptide::TBMB_positions].user() ) {
		runtime_assert_string_msg( !option[basic::options::OptionKeys::cyclic_peptide::link_all_cys_with_paraBBMB].user(), "Error in simple_cycpep_predict application: The \"-TBMB_positions\" flag and the \"-link_all_cys_with_paraBBMB\" flag cannot be used together." );
		runtime_assert_string_msg( !option[basic::options::OptionKeys::cyclic_peptide::link_all_cys_with_TBMB].user(), "Error in simple_cycpep_predict application: The \"-TBMB_positions\" flag and the \"-link_all_cys_with_TBMB\" flag cannot be used together." );
		core::Size const ntbmbres(option[basic::options::OptionKeys::cyclic_peptide::TBMB_positions]().size());
		runtime_assert_string_msg( ntbmbres > 0, "Error in simple_cycpep_predict application: The \"-cyclic_peptide:TBMB_positions\" commandline option must be followed by a list of residues to link with 1,3,5-tris(bromomethyl)benzene." );
		runtime_assert_string_msg( ntbmbres % 3 == 0, "Error in simple_cycpep_predict application: The \"-cyclic_peptide:TBMB_positions\" commandline option must be followed by a list of residues, where the number of residues in the list is a multiple of three.  Groups of three residues will be linked with 1,3,5-tris(bromomethyl)benzene." );
		core::Size count(0);
		tbmb_positions_.resize(ntbmbres / 3);
		for ( core::Size i(1), imax(ntbmbres / 3); i<=imax; ++i ) {
			utility::vector1 <core::Size> innervect(3);
			for ( core::Size j=1; j<=3; ++j ) {
				++count;
				innervect[j] = option[basic::options::OptionKeys::cyclic_peptide::TBMB_positions]()[count];
			}
			tbmb_positions_[i] = innervect;
		}
	}
	use_tbmb_filters_ = option[basic::options::OptionKeys::cyclic_peptide::use_TBMB_filters]();
	tbmb_sidechain_distance_filter_multiplier_ = option[basic::options::OptionKeys::cyclic_peptide::TBMB_sidechain_distance_filter_multiplier]();
	tbmb_constraints_energy_filter_multiplier_ = option[basic::options::OptionKeys::cyclic_peptide::TBMB_constraints_energy_filter_multiplier]();
	link_all_cys_with_tbmb_ = option[basic::options::OptionKeys::cyclic_peptide::link_all_cys_with_TBMB]();

	//Store the TMA positions.
	if ( option[basic::options::OptionKeys::cyclic_peptide::TMA_positions].user() ) {
		core::Size const n_tma_res(option[basic::options::OptionKeys::cyclic_peptide::TMA_positions]().size());
		runtime_assert_string_msg( n_tma_res > 0, "Error in simple_cycpep_predict application: The \"-cyclic_peptide:TMA_positions\" commandline option must be followed by a list of residues to link with trimesic acid." );
		runtime_assert_string_msg( n_tma_res % 3 == 0, "Error in simple_cycpep_predict application: The \"-cyclic_peptide:TMA_positions\" commandline option must be followed by a list of residues, where the number of residues in the list is a multiple of three.  Groups of three residues will be linked with trimesic acid." );
		core::Size count(0);
		tma_positions_.resize(n_tma_res / 3);
		for ( core::Size i(1), imax(n_tma_res / 3); i<=imax; ++i ) {
			utility::vector1 <core::Size> innervect(3);
			for ( core::Size j=1; j<=3; ++j ) {
				++count;
				innervect[j] = option[basic::options::OptionKeys::cyclic_peptide::TMA_positions]()[count];
			}
			tma_positions_[i] = innervect;
		}
	}
	use_tma_filters_ = option[basic::options::OptionKeys::cyclic_peptide::use_TMA_filters]();
	tma_sidechain_distance_filter_multiplier_ = option[basic::options::OptionKeys::cyclic_peptide::TMA_sidechain_distance_filter_multiplier]();
	tma_constraints_energy_filter_multiplier_ = option[basic::options::OptionKeys::cyclic_peptide::TMA_constraints_energy_filter_multiplier]();

	//Store the trigonal pyramidal metal positions.
	if ( option[basic::options::OptionKeys::cyclic_peptide::trigonal_pyramidal_metal_positions].user() ) {
		set_trigonal_pyramidal_metal_positions_from_string_vector( option[basic::options::OptionKeys::cyclic_peptide::trigonal_pyramidal_metal_positions]() );
	}
	use_trigonal_pyramidal_metal_filters_ = option[basic::options::OptionKeys::cyclic_peptide::use_trigonal_pyramidal_metal_filters]();
	trigonal_pyramidal_metal_sidechain_distance_filter_multiplier_ = option[basic::options::OptionKeys::cyclic_peptide::trigonal_pyramidal_metal_sidechain_distance_filter_multiplier]();
	trigonal_pyramidal_metal_constraints_energy_filter_multiplier_ = option[basic::options::OptionKeys::cyclic_peptide::trigonal_pyramidal_metal_constraints_energy_filter_multiplier]();

	//Store the trigonal planar metal positions.
	if ( option[basic::options::OptionKeys::cyclic_peptide::trigonal_planar_metal_positions].user() ) {
		set_trigonal_planar_metal_positions_from_string_vector( option[basic::options::OptionKeys::cyclic_peptide::trigonal_planar_metal_positions]() );
	}
	use_trigonal_planar_metal_filters_ = option[basic::options::OptionKeys::cyclic_peptide::use_trigonal_planar_metal_filters]();
	trigonal_planar_metal_sidechain_distance_filter_multiplier_ = option[basic::options::OptionKeys::cyclic_peptide::trigonal_planar_metal_sidechain_distance_filter_multiplier]();
	trigonal_planar_metal_constraints_energy_filter_multiplier_ = option[basic::options::OptionKeys::cyclic_peptide::trigonal_planar_metal_constraints_energy_filter_multiplier]();

	//Store the square pyramidal metal positions.
	if ( option[basic::options::OptionKeys::cyclic_peptide::square_pyramidal_metal_positions].user() ) {
		set_square_pyramidal_metal_positions_from_string_vector( option[basic::options::OptionKeys::cyclic_peptide::square_pyramidal_metal_positions]() );
	}
	use_square_pyramidal_metal_filters_ = option[basic::options::OptionKeys::cyclic_peptide::use_square_pyramidal_metal_filters]();
	square_pyramidal_metal_sidechain_distance_filter_multiplier_ = option[basic::options::OptionKeys::cyclic_peptide::square_pyramidal_metal_sidechain_distance_filter_multiplier]();
	square_pyramidal_metal_constraints_energy_filter_multiplier_ = option[basic::options::OptionKeys::cyclic_peptide::square_pyramidal_metal_constraints_energy_filter_multiplier]();

	//Store the square planar metal positions.
	if ( option[basic::options::OptionKeys::cyclic_peptide::square_planar_metal_positions].user() ) {
		set_square_planar_metal_positions_from_string_vector( option[basic::options::OptionKeys::cyclic_peptide::square_planar_metal_positions]() );
	}
	use_square_planar_metal_filters_ = option[basic::options::OptionKeys::cyclic_peptide::use_square_planar_metal_filters]();
	square_planar_metal_sidechain_distance_filter_multiplier_ = option[basic::options::OptionKeys::cyclic_peptide::square_planar_metal_sidechain_distance_filter_multiplier]();
	square_planar_metal_constraints_energy_filter_multiplier_ = option[basic::options::OptionKeys::cyclic_peptide::square_planar_metal_constraints_energy_filter_multiplier]();

	//Store the tetrahedral metal positions.
	if ( option[basic::options::OptionKeys::cyclic_peptide::tetrahedral_metal_positions].user() ) {
		set_tetrahedral_metal_positions_from_string_vector( option[basic::options::OptionKeys::cyclic_peptide::tetrahedral_metal_positions]() );
	}
	use_tetrahedral_metal_filters_ = option[basic::options::OptionKeys::cyclic_peptide::use_tetrahedral_metal_filters]();
	tetrahedral_metal_sidechain_distance_filter_multiplier_ = option[basic::options::OptionKeys::cyclic_peptide::tetrahedral_metal_sidechain_distance_filter_multiplier]();
	tetrahedral_metal_constraints_energy_filter_multiplier_ = option[basic::options::OptionKeys::cyclic_peptide::tetrahedral_metal_constraints_energy_filter_multiplier]();

	//Store the octahedral metal positions.
	if ( option[basic::options::OptionKeys::cyclic_peptide::octahedral_metal_positions].user() ) {
		set_octahedral_metal_positions_from_string_vector( option[basic::options::OptionKeys::cyclic_peptide::octahedral_metal_positions]() );
	}
	use_octahedral_metal_filters_ = option[basic::options::OptionKeys::cyclic_peptide::use_octahedral_metal_filters]();
	octahedral_metal_sidechain_distance_filter_multiplier_ = option[basic::options::OptionKeys::cyclic_peptide::octahedral_metal_sidechain_distance_filter_multiplier]();
	octahedral_metal_constraints_energy_filter_multiplier_ = option[basic::options::OptionKeys::cyclic_peptide::octahedral_metal_constraints_energy_filter_multiplier]();

	//Options for symmetric sampling:
	runtime_assert_string_msg( option[basic::options::OptionKeys::cyclic_peptide::require_symmetry_repeats]() > 0, "Error in simple_cycpep_predict application: The \"-cyclic_peptide:require_symmetry_repeats\" flag must be provided with a positive value." );
	runtime_assert_string_msg( option[basic::options::OptionKeys::cyclic_peptide::require_symmetry_angle_threshold]() > 0.0, "Error in simple_cycpep_predict application: The \"-cyclic_peptide:require_symmetry_angle_threshold\" flag must be provided with a positive value." );
	if ( option[basic::options::OptionKeys::cyclic_peptide::require_symmetry_mirroring]() ) {
		runtime_assert_string_msg( option[basic::options::OptionKeys::cyclic_peptide::require_symmetry_repeats].user() && option[basic::options::OptionKeys::cyclic_peptide::require_symmetry_repeats]() > 1 && option[basic::options::OptionKeys::cyclic_peptide::require_symmetry_repeats]() % 2 == 0,
			"Error in simple_cycpep_predict application: If the \"-cyclic_peptide:require_symmetry_mirroring\" option is used, then the \"-cyclic_peptide:require_symmetry_repeats\" option must be provided, must be set greater than 1, and must be set to a value divisible by 2." );
	}
	required_symmetry_repeats_ = option[basic::options::OptionKeys::cyclic_peptide::require_symmetry_repeats]();
	required_symmetry_mirroring_ = option[basic::options::OptionKeys::cyclic_peptide::require_symmetry_mirroring]();
	required_symmetry_angle_threshold_ = option[basic::options::OptionKeys::cyclic_peptide::require_symmetry_angle_threshold]();
	required_symmetry_perturbation_ = option[basic::options::OptionKeys::cyclic_peptide::require_symmetry_perturbation]();

	if ( option[basic::options::OptionKeys::cyclic_peptide::exclude_residues_from_rms].user() ) {
		core::Size const noptions( option[basic::options::OptionKeys::cyclic_peptide::exclude_residues_from_rms]().size() );
		runtime_assert_string_msg( noptions > 0, "Error in simple_cycpep_predict application: The \"-cyclic_peptide:exclude_residues_from_rms\" option must be followed by at least one residue index.");
		exclude_residues_from_rms_.reserve( noptions );
		for ( core::Size i(1); i<=noptions; ++i ) {
			runtime_assert_string_msg( option[basic::options::OptionKeys::cyclic_peptide::exclude_residues_from_rms]()[i] > 0, "Error in simple_cycpep_predict application: The \"-cyclic_peptide:exclude_residues_from_rms\" option must be followed by a vector of positive integers representing residue indices." );
			exclude_residues_from_rms_.push_back( static_cast< core::Size >(option[basic::options::OptionKeys::cyclic_peptide::exclude_residues_from_rms]()[i]) );
		}
	} else {
		exclude_residues_from_rms_.clear();
	}

	if ( option[basic::options::OptionKeys::cyclic_peptide::lariat_sidechain_index].user() ) {
		runtime_assert_string_msg( is_lariat_type( cyclization_type_ ), "Error in simple_cycpep_predict application: The \"-cyclic_peptide:lariat_sidechain_index\" option may only be used with a lariat cyclization type." );
	}
	lariat_sidechain_index_ = option[basic::options::OptionKeys::cyclic_peptide::lariat_sidechain_index]();
	lariat_sample_cis_ = option[basic::options::OptionKeys::cyclic_peptide::lariat_sample_cis]();

	if ( option[basic::options::OptionKeys::cyclic_peptide::sidechain_isopeptide_indices].user() ) {
		runtime_assert_string_msg( cyclization_type_ == SCPA_sidechain_isopeptide, "Error in simple_cycpep_predict application: The \"-cyclic_peptide:sidechain_isopeptide_indices\" option may only be used with the \"sidechain_isopeptide\" cyclization type." );
		runtime_assert_string_msg( option[basic::options::OptionKeys::cyclic_peptide::sidechain_isopeptide_indices]().size() == 2, "Error in simple_cycpep_predict application: The \"-cyclic_peptide:sidechain_isopeptide_indices\" option must be followed by exactly two residue indices." );
		runtime_assert_string_msg( option[basic::options::OptionKeys::cyclic_peptide::sidechain_isopeptide_indices]()[1] > 0 && option[basic::options::OptionKeys::cyclic_peptide::sidechain_isopeptide_indices]()[2] > 0, "Error in simple_cycpep_predict application: The \"-cyclic_peptide:sidechain_isopeptide_indices\" option must be used to provide residue indices.  Residue indices must be positive." );
		sidechain_isopeptide_indices_.first = static_cast< core::Size >(option[basic::options::OptionKeys::cyclic_peptide::sidechain_isopeptide_indices]()[1]);
		sidechain_isopeptide_indices_.second = static_cast< core::Size >(option[basic::options::OptionKeys::cyclic_peptide::sidechain_isopeptide_indices]()[2]);
	} else { \
		//If not specified:
		sidechain_isopeptide_indices_.first = 0;
		sidechain_isopeptide_indices_.second = 0;
	}

	if ( cyclization_type_ != SCPA_n_to_c_amide_bond ) {
		runtime_assert_string_msg( !bondlength_perturbation_magnitude_, "Error in simple_cycpep_predict application: The \"-bondlength_perturbation_magnitude\" option may only be used with N-to-C cyclization." );
		runtime_assert_string_msg( !bondangle_perturbation_magnitude_, "Error in simple_cycpep_predict application: The \"-bondangle_perturbation_magnitude\" option may only be used with N-to-C cyclization." );
	}

	return;
} //initialize_from_options()


/// @brief Set the cyclization type.
///
void
SimpleCycpepPredictApplication::set_cyclization_type(
	SCPA_cyclization_type const type_in
) {
	runtime_assert_string_msg( type_in != SCPA_invalid_type, "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::set_cyclization_type(): The type provided is invalid!" );
	cyclization_type_ = type_in;
}

/// @brief Set whether we should use the chainbreak energy (true) or constraints (false) to enforce
/// terminal amide bond geometry.
void
SimpleCycpepPredictApplication::set_use_chainbreak_energy(
	bool const setting
) {
	use_chainbreak_energy_ = setting;
}


/// @brief Sets the default scorefunction to use.
/// @details The scorefunction is cloned.  The high-hbond version is constructed
/// from this one.  If necessary, the aa_composition score term will be turned on
/// in that one; it needn't be turned on in this one.
void
SimpleCycpepPredictApplication::set_scorefxn(
	core::scoring::ScoreFunctionCOP sfxn_in
) {
	debug_assert( sfxn_in );
	scorefxn_ = sfxn_in->clone();
}

/// @brief Allows external code to provide a native, so that the SimpleCycpepPredictApplication doesn't have to read
/// directly from disk.
void
SimpleCycpepPredictApplication::set_native(
	core::pose::PoseCOP native
) {
	runtime_assert(native); //Can't be NULL.
	native_exists_=true;
	core::pose::PoseOP native_pose_copy( native->clone() );
	set_up_native( native_pose_copy, 0);
	native_pose_ = native_pose_copy;
}

/// @brief Allows external code to provide a sequence, so that the SimpleCycpepPredictApplication doesn't have to read
/// directly from disk.
void
SimpleCycpepPredictApplication::set_sequence(
	std::string const &seq
) {
	runtime_assert( seq != "" ); //Can't be empty string.
	sequence_string_ = seq;
}

/// @brief Allows external code to set the allowed residues by position, so that this needn't be read directly
/// from disk.
void
SimpleCycpepPredictApplication::set_allowed_residues_by_position (
	std::map< core::Size, utility::vector1< std::string > > const &allowed_canonicals,
	std::map< core::Size, utility::vector1< std::string > > const &allowed_noncanonicals
) {
	allowed_canonicals_by_position_ = allowed_canonicals;
	allowed_noncanonicals_by_position_ = allowed_noncanonicals;
}

/// @brief Allows external code to specify that output should be appended to a list of SilentStructureOPs, so that the
/// SimpleCycpepPredictApplication doesn't have to write directly to disk.
void
SimpleCycpepPredictApplication::set_silentstructure_outputlist(
	utility::vector1 < core::io::silent::SilentStructOP > * silentlist,
	utility::vector1 < HierarchicalHybridJD_JobResultsSummaryOP > * summarylist
) {
	runtime_assert( silentlist && summarylist );
	silentlist_ = silentlist;
	summarylist_ = summarylist;
	silentlist_out_ = true;
	silent_out_ = false;
}

/// @brief Allows external code to suppress checkpointing, to prevent direct file I/O from disk.
/// @details Useful on Blue Gene.
void
SimpleCycpepPredictApplication::set_suppress_checkpoints(
	bool const suppress_checkpoints
) {
	suppress_checkpoints_ = suppress_checkpoints;
}

/// @brief If called by MPI code, the rank of the current process can be stored here.
/// @details Used for output of job summaries.
void
SimpleCycpepPredictApplication::set_my_rank(
	int const rank_in
) {
	my_rank_ = rank_in;
}

/// @brief Set the number of jobs that this process has already completed.
///
void
SimpleCycpepPredictApplication::set_already_completed_job_count(
	core::Size const count_in
) {
	already_completed_job_count_ = count_in;
}

/// @brief Allows external code to override the number of structures that this should generate (otherwise
/// set by options system.
void
SimpleCycpepPredictApplication::set_nstruct(
	core::Size const nstruct_in
) {
	runtime_assert( nstruct_in > 0 );
	nstruct_ = nstruct_in;
}

/// @brief Allows external code to set the file contents for the L-alpha aa_composition file.
///
void
SimpleCycpepPredictApplication::set_L_alpha_compfile_contents(
	std::string const &contents_in
) {
	comp_file_contents_L_alpha_ = contents_in;
	L_alpha_comp_file_exists_=true;
}

/// @brief Allows external code to set the file contents for the D-alpha aa_composition file.
///
void
SimpleCycpepPredictApplication::set_D_alpha_compfile_contents(
	std::string const &contents_in
) {
	comp_file_contents_D_alpha_ = contents_in;
	D_alpha_comp_file_exists_=true;
}

/// @brief Allows external code to set the file contents for the L-beta aa_composition file.
///
void
SimpleCycpepPredictApplication::set_L_beta_compfile_contents(
	std::string const &contents_in
) {
	comp_file_contents_L_beta_ = contents_in;
	L_beta_comp_file_exists_=true;
}

/// @brief Allows external code to set the file contents for the D-beta aa_composition file.
///
void
SimpleCycpepPredictApplication::set_D_beta_compfile_contents(
	std::string const &contents_in
) {
	comp_file_contents_D_beta_ = contents_in;
	D_beta_comp_file_exists_=true;
}

/// @brief Set the bin transitions file.
void
SimpleCycpepPredictApplication::set_abba_bins_binfile(
	std::string const &binfile_in
) {
	abba_bins_binfile_ = binfile_in;
}

/// @brief Set the frequency with which we sample cis proline.
/// @details Implicitly sets sample_cis_pro_ to "true" if freq_in is not 0.0, "false" if it is.
void
SimpleCycpepPredictApplication::set_sample_cis_pro_frequency(
	core::Real const &freq_in
) {
	runtime_assert_string_msg( 0.0 <= freq_in && freq_in <= 1.0, "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::set_sample_cis_pro_frequency(): The frequency must be between 0 and 1." );
	sample_cis_pro_frequency_ = freq_in;
	if ( freq_in > 0.0 ) {
		sample_cis_pro_ = true;
	} else {
		sample_cis_pro_ = false;
		TR << "Disabling cis-proline sampling." << std::endl;
	}
}

/// @brief Set cis proline sampling OFF.
///
void
SimpleCycpepPredictApplication::disable_cis_pro_sampling() {
	sample_cis_pro_ = false;
	sample_cis_pro_frequency_ = 0.0;
}

/// @brief Set the total energy cutoff.
/// @details Also sets use_total_energy_cutoff_ to 'true'.
void
SimpleCycpepPredictApplication::set_total_energy_cutoff(
	core::Real const &value_in
) {
	total_energy_cutoff_ = value_in;
	use_total_energy_cutoff_ = true;
}

/// @brief Sets use_total_energy_cutoff_ to 'false'.
///
void
SimpleCycpepPredictApplication::disable_total_energy_cutoff() {
	use_total_energy_cutoff_ = false;
}

/// @brief Set the number of rounds of relaxation with flexible
/// bond angles.
void
SimpleCycpepPredictApplication::set_angle_relax_rounds(
	core::Size const rounds_in
) {
	angle_relax_rounds_ = rounds_in;
}

/// @brief Set the number of rounds of relaxation with flexible
/// bond angles and bond lengths.
void
SimpleCycpepPredictApplication::set_angle_length_relax_rounds(
	core::Size const rounds_in
) {
	angle_length_relax_rounds_ = rounds_in;
}

/// @brief Set the number of rounds of Cartesian relaxation.
///
void
SimpleCycpepPredictApplication::set_cartesian_relax_rounds(
	core::Size const rounds_in
) {
	cartesian_relax_rounds_ = rounds_in;
}

/// @brief Given an input vector of strings of the form "res1,res2,res3,metal_name", parse this and populate
/// the trigonal_pyramidal_metal_positions_ vector.
/// @details Resets the trigonal_pyramidal_metal_positions_ vector.
void
SimpleCycpepPredictApplication::set_trigonal_pyramidal_metal_positions_from_string_vector(
	utility::vector1< std::string > const &vect
) {
	core::Size const nstrings( vect.size() );
	reset_trigonal_pyramidal_metal_positions();
	if ( !nstrings ) return; //Do nothing more for size 0 vector.

	utility::fixedsizearray1< int, 3 > resnums; //Deliberately ints.
	std::string element;

	for ( core::Size i(1); i<=nstrings; ++i ) {
		std::istringstream ss( vect[i] );
		core::Size counter(0);
		while ( std::getline(ss, element, ',') ) {
			++counter;
			if ( counter < 4 ) {
				resnums[counter] = std::atoi(element.c_str());
				runtime_assert_string_msg( resnums[counter] > 0, "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::set_trigonal_pyramidal_metal_positions_from_string_vector(): All residue numbers must be greater than zero when setting residue indices for metals coordinated with trigonal pyramidal geometry." );
			} else if ( counter > 4 ) {
				utility_exit_with_message( "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::set_trigonal_pyramidal_metal_positions_from_string_vector(): A metal setup input flag for trigonal pyramidal metals was provided that was not of the form \"res1,res2,res3,res4,metal_name\".  Comma-separated lists of five entries (four numbers and a metal name) are required." );
			}
		}
		add_entry_to_trigonal_pyramidal_metal_positions(static_cast<core::Size>(resnums[1]), static_cast<core::Size>(resnums[2]), static_cast<core::Size>(resnums[3]), element);
		TR << "Parsed trigonal pyramidal metal setup where res1=" << resnums[1] << ", res2=" << resnums[2] << ", res3=" << resnums[3] << ", metal=" << element << std::endl;
	}

}

/// @brief Given an input vector of strings of the form "res1,res2,res3,metal_name", parse this and populate
/// the trigonal_planar_metal_positions_ vector.
/// @details Resets the trigonal_planar_metal_positions_ vector.
void
SimpleCycpepPredictApplication::set_trigonal_planar_metal_positions_from_string_vector(
	utility::vector1< std::string > const &vect
) {
	core::Size const nstrings( vect.size() );
	reset_trigonal_planar_metal_positions();
	if ( !nstrings ) return; //Do nothing more for size 0 vector.

	utility::fixedsizearray1< int, 3 > resnums; //Deliberately ints.
	std::string element;

	for ( core::Size i(1); i<=nstrings; ++i ) {
		std::istringstream ss( vect[i] );
		core::Size counter(0);
		while ( std::getline(ss, element, ',') ) {
			++counter;
			if ( counter < 4 ) {
				resnums[counter] = std::atoi(element.c_str());
				runtime_assert_string_msg( resnums[counter] > 0, "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::set_trigonal_planar_metal_positions_from_string_vector(): All residue numbers must be greater than zero when setting residue indices for metals coordinated with trigonal planar geometry." );
			} else if ( counter > 4 ) {
				utility_exit_with_message( "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::set_trigonal_planar_metal_positions_from_string_vector(): A metal setup input flag for trigonal planar metals was provided that was not of the form \"res1,res2,res3,res4,metal_name\".  Comma-separated lists of five entries (four numbers and a metal name) are required." );
			}
		}
		add_entry_to_trigonal_planar_metal_positions(static_cast<core::Size>(resnums[1]), static_cast<core::Size>(resnums[2]), static_cast<core::Size>(resnums[3]), element);
		TR << "Parsed trigonal planar metal setup where res1=" << resnums[1] << ", res2=" << resnums[2] << ", res3=" << resnums[3] << ", metal=" << element << std::endl;
	}
}

/// @brief Given an input vector of strings of the form "res1,res2,res3,res4,res5,metal_name", parse this and populate
/// the square_pyramidal_metal_positions_ vector.
/// @details Resets the square_pyramidal_metal_positions_ vector.
void
SimpleCycpepPredictApplication::set_square_pyramidal_metal_positions_from_string_vector(
	utility::vector1< std::string > const &vect
) {
	core::Size const nstrings( vect.size() );
	reset_square_pyramidal_metal_positions();
	if ( !nstrings ) return; //Do nothing more for size 0 vector.

	utility::fixedsizearray1< int, 5 > resnums; //Deliberately ints.
	std::string element;

	for ( core::Size i(1); i<=nstrings; ++i ) {
		std::istringstream ss( vect[i] );
		core::Size counter(0);
		while ( std::getline(ss, element, ',') ) {
			++counter;
			if ( counter < 6 ) {
				resnums[counter] = std::atoi(element.c_str());
				runtime_assert_string_msg( resnums[counter] > 0, "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::set_square_pyramidal_metal_positions_from_string_vector(): All residue numbers must be greater than zero when setting residue indices for metals coordinated with square pyramidal geometry." );
			} else if ( counter > 6 ) {
				utility_exit_with_message( "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::set_square_pyramidal_metal_positions_from_string_vector(): A metal setup input flag was provided for square pyramidal geometry that was not of the form \"res1,res2,res3,res4,res5,metal_name\".  Comma-separated lists of five entries (four numbers and a metal name) are required." );
			}
		}
		add_entry_to_square_pyramidal_metal_positions(static_cast<core::Size>(resnums[1]), static_cast<core::Size>(resnums[2]), static_cast<core::Size>(resnums[3]), static_cast<core::Size>(resnums[4]), static_cast<core::Size>(resnums[5]), element);
		TR << "Parsed square pyramidal metal setup where res1=" << resnums[1] << ", res2=" << resnums[2] << ", res3=" << resnums[3] << ", res4=" << resnums[4] << ", res5=" << resnums[5] << ", metal=" << element << std::endl;
	}
}

/// @brief Given an input vector of strings of the form "res1,res2,res3,res4,metal_name", parse this and populate
/// the square_planar_metal_positions_ vector.
/// @details Resets the square_planar_metal_positions_ vector.
void
SimpleCycpepPredictApplication::set_square_planar_metal_positions_from_string_vector(
	utility::vector1< std::string > const &vect
) {
	core::Size const nstrings( vect.size() );
	reset_square_planar_metal_positions();
	if ( !nstrings ) return; //Do nothing more for size 0 vector.

	utility::fixedsizearray1< int, 4 > resnums; //Deliberately ints.
	std::string element;

	for ( core::Size i(1); i<=nstrings; ++i ) {
		std::istringstream ss( vect[i] );
		core::Size counter(0);
		while ( std::getline(ss, element, ',') ) {
			++counter;
			if ( counter < 5 ) {
				resnums[counter] = std::atoi(element.c_str());
				runtime_assert_string_msg( resnums[counter] > 0, "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::set_square_planar_metal_positions_from_string_vector(): All residue numbers must be greater than zero when setting residue indices for metals coordinated with square planar geometry." );
			} else if ( counter > 5 ) {
				utility_exit_with_message( "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::set_square_planar_metal_positions_from_string_vector(): A metal setup input flag was provided for square planar geometry that was not of the form \"res1,res2,res3,res4,metal_name\".  Comma-separated lists of five entries (four numbers and a metal name) are required." );
			}
		}
		add_entry_to_square_planar_metal_positions(static_cast<core::Size>(resnums[1]), static_cast<core::Size>(resnums[2]), static_cast<core::Size>(resnums[3]), static_cast<core::Size>(resnums[4]), element);
		TR << "Parsed square planar metal setup where res1=" << resnums[1] << ", res2=" << resnums[2] << ", res3=" << resnums[3] << ", res4=" << resnums[4] << ", metal=" << element << std::endl;
	}
}


/// @brief Given an input vector of strings of the form "res1,res2,res3,res4,metal_name", parse this and populate
/// the tetrahedral_metal_positions_ vector.
/// @details Resets the tetrahedral_metal_positions_ vector.
void
SimpleCycpepPredictApplication::set_tetrahedral_metal_positions_from_string_vector(
	utility::vector1< std::string > const &vect
) {
	core::Size const nstrings( vect.size() );
	reset_tetrahedral_metal_positions();
	if ( !nstrings ) return; //Do nothing more for size 0 vector.

	utility::fixedsizearray1< int, 4 > resnums; //Deliberately ints.
	std::string element;

	for ( core::Size i(1); i<=nstrings; ++i ) {
		std::istringstream ss( vect[i] );
		core::Size counter(0);
		while ( std::getline(ss, element, ',') ) {
			++counter;
			if ( counter < 5 ) {
				resnums[counter] = std::atoi(element.c_str());
				runtime_assert_string_msg( resnums[counter] > 0, "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::set_tetrahedral_metal_positions_from_string_vector(): All residue numbers must be greater than zero when setting residue indices for tetrahedrally-coordinated metals." );
			} else if ( counter > 5 ) {
				utility_exit_with_message( "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::set_tetrahedral_metal_positions_from_string_vector(): A tetrahedrally-coordinated metal setup input flag was provided that was not of the form \"res1,res2,res3,res4,metal_name\".  Comma-separated lists of five entries (four numbers and a metal name) are required." );
			}
		}
		add_entry_to_tetrahedral_metal_positions(static_cast<core::Size>(resnums[1]), static_cast<core::Size>(resnums[2]), static_cast<core::Size>(resnums[3]), static_cast<core::Size>(resnums[4]), element);
		TR << "Parsed tetrahedral metal setup where res1=" << resnums[1] << ", res2=" << resnums[2] << ", res3=" << resnums[3] << ", res4=" << resnums[4] << ", metal=" << element << std::endl;
	}
}

/// @brief Given an input vector of strings of the form "res1,res2,res3,res4,metal_name", parse this and populate
/// the octahedral_metal_positions_ vector.
/// @details Resets the octahedral_metal_positions_ vector.
void
SimpleCycpepPredictApplication::set_octahedral_metal_positions_from_string_vector(
	utility::vector1< std::string > const &vect
) {
	core::Size const nstrings( vect.size() );
	reset_octahedral_metal_positions();
	if ( !nstrings ) return; //Do nothing more for size 0 vector.

	utility::fixedsizearray1< int, 6 > resnums; //Deliberately ints.
	std::string element;

	for ( core::Size i(1); i<=nstrings; ++i ) {
		std::istringstream ss( vect[i] );
		core::Size counter(0);
		while ( std::getline(ss, element, ',') ) {
			++counter;
			if ( counter < 7 ) {
				resnums[counter] = std::atoi(element.c_str());
				runtime_assert_string_msg( resnums[counter] > 0, "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::set_octahedral_metal_positions_from_string_vector(): All residue numbers must be greater than zero when setting residue indices for octahedrally-coordinated metals." );
			} else if ( counter > 7 ) {
				utility_exit_with_message( "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::set_octahedral_metal_positions_from_string_vector(): A octahedrally-coordinated metal setup input flag was provided that was not of the form \"res1,res2,res3,res4,metal_name\".  Comma-separated lists of seven entries (six numbers and a metal name) are required." );
			}
		}
		add_entry_to_octahedral_metal_positions(static_cast<core::Size>(resnums[1]), static_cast<core::Size>(resnums[2]), static_cast<core::Size>(resnums[3]), static_cast<core::Size>(resnums[4]), static_cast<core::Size>(resnums[5]), static_cast<core::Size>(resnums[6]), element);
		TR << "Parsed octahedral metal setup where res1=" << resnums[1] << ", res2=" << resnums[2] << ", res3=" << resnums[3] << ", res4=" << resnums[4] << ", res5=" << resnums[5] << ", res6=" << resnums[6] << ", metal=" << element << std::endl;
	}
}

/// @brief Resets the trigonal_pyramidal_metal_positions_ vector.
void
SimpleCycpepPredictApplication::reset_trigonal_pyramidal_metal_positions() {
	trigonal_pyramidal_metal_positions_.clear();
}

/// @brief Resets the trigonal_planar_metal_positions_ vector.
void
SimpleCycpepPredictApplication::reset_trigonal_planar_metal_positions() {
	trigonal_planar_metal_positions_.clear();
}

/// @brief Resets the square_pyramidal_metal_positions_ vector.
void
SimpleCycpepPredictApplication::reset_square_pyramidal_metal_positions() {
	square_pyramidal_metal_positions_.clear();
}

/// @brief Resets the square_planar_metal_positions_ vector.
void
SimpleCycpepPredictApplication::reset_square_planar_metal_positions() {
	square_planar_metal_positions_.clear();
}

/// @brief Resets the tetrahedral_metal_positions_ vector.
void
SimpleCycpepPredictApplication::reset_tetrahedral_metal_positions() {
	tetrahedral_metal_positions_.clear();
}

/// @brief Resets the octahedral_metal_positions_ vector.
void
SimpleCycpepPredictApplication::reset_octahedral_metal_positions() {
	octahedral_metal_positions_.clear();
}

/// @brief Adds an entry to the trigonal_pyramidal_metal_positions_ vector.
void
SimpleCycpepPredictApplication::add_entry_to_trigonal_pyramidal_metal_positions(
	core::Size const res1,
	core::Size const res2,
	core::Size const res3,
	std::string const & metal_type
) {
	utility::fixedsizearray1< core::Size, 3 > residues;
	runtime_assert(res1 > 0 && res2 > 0 && res3 > 0 );
	runtime_assert(!metal_type.empty());
	residues[1] = res1;
	residues[2] = res2;
	residues[3] = res3;
	trigonal_pyramidal_metal_positions_.push_back( std::make_pair( residues, metal_type) );
}

/// @brief Adds an entry to the trigonal_planar_metal_positions_ vector.
void
SimpleCycpepPredictApplication::add_entry_to_trigonal_planar_metal_positions(
	core::Size const res1,
	core::Size const res2,
	core::Size const res3,
	std::string const & metal_type
) {
	utility::fixedsizearray1< core::Size, 3 > residues;
	runtime_assert(res1 > 0 && res2 > 0 && res3 > 0 );
	runtime_assert(!metal_type.empty());
	residues[1] = res1;
	residues[2] = res2;
	residues[3] = res3;
	trigonal_planar_metal_positions_.push_back( std::make_pair( residues, metal_type) );
}

/// @brief Adds an entry to the square_pyramidal_metal_positions_ vector.
void
SimpleCycpepPredictApplication::add_entry_to_square_pyramidal_metal_positions(
	core::Size const res1,
	core::Size const res2,
	core::Size const res3,
	core::Size const res4,
	core::Size const res5,
	std::string const & metal_type
) {
	utility::fixedsizearray1< core::Size, 5 > residues;
	runtime_assert(res1 > 0 && res2 > 0 && res3 > 0 && res4 > 0 && res5 > 0);
	runtime_assert(!metal_type.empty());
	residues[1] = res1;
	residues[2] = res2;
	residues[3] = res3;
	residues[4] = res4;
	residues[5] = res5;
	square_pyramidal_metal_positions_.push_back( std::make_pair( residues, metal_type) );
}

/// @brief Adds an entry to the square_planar_metal_positions_ vector.
void
SimpleCycpepPredictApplication::add_entry_to_square_planar_metal_positions(
	core::Size const res1,
	core::Size const res2,
	core::Size const res3,
	core::Size const res4,
	std::string const & metal_type
) {
	utility::fixedsizearray1< core::Size, 4 > residues;
	runtime_assert(res1 > 0 && res2 > 0 && res3 > 0 && res4 > 0);
	runtime_assert(!metal_type.empty());
	residues[1] = res1;
	residues[2] = res2;
	residues[3] = res3;
	residues[4] = res4;
	square_planar_metal_positions_.push_back( std::make_pair( residues, metal_type) );
}

/// @brief Adds an entry to the tetrahedral_metal_positions_ vector.
void
SimpleCycpepPredictApplication::add_entry_to_tetrahedral_metal_positions(
	core::Size const res1,
	core::Size const res2,
	core::Size const res3,
	core::Size const res4,
	std::string const & metal_type
) {
	utility::fixedsizearray1< core::Size, 4 > residues;
	runtime_assert(res1 > 0 && res2 > 0 && res3 > 0 && res4 > 0);
	runtime_assert(!metal_type.empty());
	residues[1] = res1;
	residues[2] = res2;
	residues[3] = res3;
	residues[4] = res4;
	tetrahedral_metal_positions_.push_back( std::make_pair( residues, metal_type) );
}

/// @brief Adds an entry to the octahedral_metal_positions_ vector.
void
SimpleCycpepPredictApplication::add_entry_to_octahedral_metal_positions(
	core::Size const res1,
	core::Size const res2,
	core::Size const res3,
	core::Size const res4,
	core::Size const res5,
	core::Size const res6,
	std::string const & metal_type
) {
	utility::fixedsizearray1< core::Size, 6 > residues;
	runtime_assert(res1 > 0 && res2 > 0 && res3 > 0 && res4 > 0 && res5 > 0 && res6 > 0);
	runtime_assert(!metal_type.empty());
	residues[1] = res1;
	residues[2] = res2;
	residues[3] = res3;
	residues[4] = res4;
	residues[5] = res5;
	residues[6] = res6;
	octahedral_metal_positions_.push_back( std::make_pair( residues, metal_type) );
}


/// @brief Set whether we're using RamaPrePro tables for sampling.
/// @details Setting this to "false" lets us use classic rama tables.  True by default.
void
SimpleCycpepPredictApplication::set_use_rama_prepro_for_sampling(
	bool const setting
) {
	use_rama_prepro_for_sampling_ = setting;
}

/// @brief Align pose to native_pose, and return the RMSD between the two poses.
/// @details Assumes that the pose has already been de-permuted (i.e. the native and the pose line up).
/// Only uses alpha-amino acids for the alignment, currently.
core::Real
SimpleCycpepPredictApplication::align_and_calculate_rmsd(
	core::pose::Pose & pose,
	core::pose::Pose const & native_pose,
	bool const skip_seq_comparison/*=false*/
) const {
	core::Size const nres( sequence_length() );
	core::Size res_counter(0); //Residue indices might not match between native pose and pose, due to linkers.

	core::id::AtomID_Map< core::id::AtomID > amap;
	core::pose::initialize_atomid_map(amap, pose, core::id::AtomID::BOGUS_ATOM_ID());
	for ( core::Size ir=1, irmax=native_pose.total_residue(); ir<=irmax; ++ir ) {
		if ( !is_supported_restype( native_pose.residue_type(ir) ) ) continue;
		if ( is_residue_ignored_in_rms(ir) ) continue;
		++res_counter;
		if ( !skip_seq_comparison ) {
			runtime_assert_string_msg( res_counter <= nres, "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::align_and_calculate_rmsd(): The native pose has more residues than the input sequence.  (Native pose: " + std::to_string( irmax ) + " residues; sequence: " + std::to_string( nres ) + " residues.)" );
		}
		for ( core::Size ia=1, iamax=native_pose.residue_type(ir).first_sidechain_atom(); ia<iamax; ++ia ) { //Loop through all mainchain heavyatoms (including atoms coming off mainchain that are not sidechain atoms, like peptide "O").
			if ( native_pose.residue_type(ir).atom_is_hydrogen(ia) ) continue;
			if ( native_pose.residue_type(ir).is_virtual(ia) ) continue; //Ignore virtual atoms
			core::Size const ia_pose( pose.residue_type(res_counter).atom_index( native_pose.residue_type(ir).atom_name(ia) ) );
			//TR << "ir=" << ir << " ia=" << ia << " res_counter=" << res_counter << " native=" << native_pose.residue_type(ir).atom_name(ia) << " pred=" << pose.residue_type(res_counter).atom_name(ia_pose) << std::endl; TR.flush(); //DELETE ME.
			amap[ core::id::AtomID(ia_pose,res_counter) ] = core::id::AtomID(ia,ir);
			//TR << "Adding ia=" << ia << " ir=" << ir << " to map." << std::endl; //DELETE ME
		}
	}

	if ( !skip_seq_comparison ) {
		runtime_assert_string_msg( res_counter == nres - exclude_residues_from_rms_.size(), "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::align_and_calculate_rmsd(): The native pose has fewer residues than the input sequence." );
	}
	core::Real rms_val;
	core::Real offset_val( 1.0e-7 );
	bool failed = true;
	for ( core::Size i(1); i<=3; ++i ) {
		try {
			rms_val = core::scoring::superimpose_pose( pose, native_pose, amap, offset_val, false, true ); //Superimpose the pose and return the RMSD.
		} catch ( ... ) {
			TR.Warning << "Warning in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::align_and_calculate_rmsd(): RMSD calculation failed with offset = " << offset_val << std::endl;
			offset_val *= 10.0;
			continue;
		}
		failed = false;
		break;
	}
	if ( failed ) {
		TR.Warning << "Warning in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::align_and_calculate_rmsd(): RMSD calculation failed!  Returning -1.0 to indicate failure." << std::endl;
		rms_val = -1.0;
	}
	return rms_val;
}

/// @brief Actually run the application.
/// @details The initialize_from_options() function must be called before calling this.  (Called by default constructor.)
void
SimpleCycpepPredictApplication::run() const {

	//Get the scorefunction:
	debug_assert( scorefxn_ );
	core::scoring::ScoreFunctionOP sfxn_default( scorefxn_ );
	//Create a scorefunction variant with upweighted backbone hbond terms:
	core::scoring::ScoreFunctionOP sfxn_highhbond( sfxn_default->clone() );
	sfxn_highhbond->set_weight( core::scoring::hbond_lr_bb, high_hbond_weight_multiplier_ * sfxn_default->get_weight(core::scoring::hbond_lr_bb) ); //Upweight the long-range backbone hbonds
	sfxn_highhbond->set_weight( core::scoring::hbond_sr_bb, high_hbond_weight_multiplier_ * sfxn_default->get_weight(core::scoring::hbond_sr_bb) ); //Upweight the short-range backbone hbonds
	//Turn on aa_composition in this scorefunction if we're doing design and .comp files have been provided.
	if ( sfxn_highhbond->get_weight( core::scoring::aa_composition ) == 0.0 && use_aa_comp_ ) { sfxn_highhbond->set_weight( core::scoring::aa_composition, 1.0 ); }
	//Create variants of the above two scorefunctions with constraint weights turned on:
	core::scoring::ScoreFunctionOP sfxn_default_cst( sfxn_default->clone() ); //Will NOT have aa_compostion turned on.
	core::scoring::ScoreFunctionOP sfxn_highhbond_cst( sfxn_highhbond->clone() ); //Will have aa_composition turned on if we're doing design and .comp files are provided.
	if ( sfxn_default->get_weight( core::scoring::atom_pair_constraint ) == 0.0 ) { sfxn_default_cst->set_weight( core::scoring::atom_pair_constraint, 1.0); }
	if ( sfxn_default->get_weight( core::scoring::angle_constraint ) == 0.0 ) { sfxn_default_cst->set_weight( core::scoring::angle_constraint, 1.0); }
	if ( sfxn_default->get_weight( core::scoring::dihedral_constraint ) == 0.0 ) { sfxn_default_cst->set_weight( core::scoring::dihedral_constraint, 1.0); }
	if ( use_chainbreak_energy() && sfxn_default->get_weight( core::scoring::chainbreak ) == 0.0 ) { sfxn_default_cst->set_weight( core::scoring::chainbreak, 1.0 ); }
	if ( sfxn_highhbond->get_weight( core::scoring::atom_pair_constraint ) == 0.0 ) { sfxn_highhbond_cst->set_weight( core::scoring::atom_pair_constraint, 1.0); }
	if ( sfxn_highhbond->get_weight( core::scoring::angle_constraint ) == 0.0 ) { sfxn_highhbond_cst->set_weight( core::scoring::angle_constraint, 1.0); }
	if ( sfxn_highhbond->get_weight( core::scoring::dihedral_constraint ) == 0.0 ) { sfxn_highhbond_cst->set_weight( core::scoring::dihedral_constraint, 1.0); }
	if ( use_chainbreak_energy() && sfxn_highhbond->get_weight( core::scoring::chainbreak ) == 0.0 ) { sfxn_highhbond_cst->set_weight( core::scoring::chainbreak, 1.0 ); }
	core::scoring::ScoreFunctionOP sfxn_highhbond_cst_cart, sfxn_default_cst_cart;
	if ( angle_relax_rounds() > 0 || angle_length_relax_rounds() > 0 || cartesian_relax_rounds() > 0 ) {
		sfxn_highhbond_cst_cart = sfxn_highhbond_cst->clone();
		debug_assert( sfxn_highhbond_cst_cart );
		if ( sfxn_highhbond_cst_cart->get_weight( core::scoring::cart_bonded ) == 0.0 ) {
			TR << "Activating cart_bonded (with weight=0.5) for high-hbonds scorefunction Cartesian variant, and de-activating pro_close term." << std::endl;
			sfxn_highhbond_cst_cart->set_weight( core::scoring::cart_bonded, 0.5 );
			sfxn_highhbond_cst_cart->set_weight( core::scoring::pro_close, 0.0 );
		}
		sfxn_default_cst_cart = sfxn_default_cst->clone();
		debug_assert( sfxn_default_cst_cart );
		if ( sfxn_default_cst_cart->get_weight( core::scoring::cart_bonded ) == 0.0 ) {
			TR << "Activating cart_bonded (with weight=0.5) for default scorefunction Cartesian variant, and de-activating pro_close term." << std::endl;
			sfxn_default_cst_cart->set_weight( core::scoring::cart_bonded, 0.5 );
			sfxn_default_cst_cart->set_weight( core::scoring::pro_close, 0.0 );
		}
	}

	//Get the sequence that we're considering:
	utility::vector1 < std::string > resnames;
	read_sequence( sequence_file_, resnames );
	sequence_length_ = resnames.size(); //Store the number of residues in the sequence, excluding crosslinkers.

	//Check that, if we're doing cyclization through a sidechain, the loop formed is long enough:
	check_loop_length( resnames );

	//Check that, if we're enforcing symmetry, the sequence length is an integer multiple of the number of symmetry repeats:
	if ( required_symmetry_repeats_ > 1 ) {
		runtime_assert_string_msg(
			sequence_length_ % required_symmetry_repeats_ == 0,
			"Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::run(): Symmetry has been specified by the user, but the number of residues in the peptide is not an integral multiple of the number of symmetry repeats."
		);
	}

	//Check that, if the link_all_cys_with_paraBBMB flag is used, the sequence has exactly two cysteine residues:
	if ( link_all_cys_with_parabbmb_ ) {
		debug_assert( parabbmb_positions_.size() == 0 ); //Should be true
		utility::vector1 < core::Size > cys_positions;
		for ( core::Size i=1; i<=sequence_length(); ++i ) {
			if ( resnames[i] == "CYS" || resnames[i] == "DCYS" ) cys_positions.push_back(i);
		}
		runtime_assert_string_msg( cys_positions.size() == 2, "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::run(): The \"-cyclic_peptide:link_all_cys_with_paraBBMB\" flag was used, but the sequence does not contain exactly two CYS/DCYS residues." );
		parabbmb_positions_.push_back( cys_positions );
	}

	//Check that, if the link_all_cys_with_TBMB flag is used, the sequence has exactly three cysteine residues:
	if ( link_all_cys_with_tbmb_ ) {
		debug_assert( tbmb_positions_.size() == 0 ); //Should be true
		utility::vector1 < core::Size > cys_positions;
		for ( core::Size i=1; i<=sequence_length(); ++i ) {
			if ( resnames[i] == "CYS" || resnames[i] == "DCYS" ) cys_positions.push_back(i);
		}
		runtime_assert_string_msg( cys_positions.size() == 3, "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::run(): The \"-cyclic_peptide:link_all_cys_with_TBMB\" flag was used, but the sequence does not contain exactly three CYS/DCYS residues." );
		tbmb_positions_.push_back( cys_positions );
	}

	//Get the native sequence that we will compare to.
	core::pose::PoseOP native_pose;
	if ( native_exists_ ) {
		if ( native_pose_ ) {
			native_pose=native_pose_->clone();
		} else {
			native_pose=utility::pointer::make_shared< core::pose::Pose >();
			TR << "Importing native structure from " << native_filename_ << "." << std::endl;
			import_and_set_up_native ( native_filename_, native_pose, resnames.size() );
		}
#ifdef BOINC_GRAPHICS
		// set native for graphics
		boinc::Boinc::set_graphics_native_pose( *native_pose );
#endif
	} else {
		TR << "No native structure specified by the user.  No RMSD values will be calculated." << std::endl;
	}

	//Set up a filter for total number of hbonds:
	protocols::filters::FilterOP total_hbond( set_up_hbond_filter( min_genkic_hbonds_) );

	//Get the checkpoint information:
	core::Size success_count(0);
	core::Size curstruct(0);
	initialize_checkpointing( curstruct, success_count );

#ifdef BOINC
	//Max attempts before exiting without success
	// This is used to hopefully prevent finishing without a model due to filtering logic.
	// Avoid setting this too high which may cause jobs to run too long for BOINC users depending on the filtering
	// rate of success.
	core::Size const max_attempts(20);
#endif

	//EVERYTHING ABOVE THIS POINT IS DONE ONCE PER PROGRAM EXECUTION.
	++curstruct;
	for ( core::Size irepeat=curstruct, irepeat_max=nstruct_; irepeat<=irepeat_max; ++irepeat ) { //Loop nstruct times
#ifdef BOINC_GRAPHICS
		{ //Increment the model count for BOINC.
			protocols::boinc::BoincSharedMemory* shmem = protocols::boinc::Boinc::get_shmem();
			shmem->model_count = shmem->model_count + 1;
		}
#endif

		//JAB - Basic checkpointing by file ala JD2 for PDBs.
		//Not sure how to reconcile this with the checkpointing system, which is very much not this.
		if ( ! silent_out_ && ! silentlist_out_ && ! basic::options::option[basic::options::OptionKeys::out::overwrite]() ) {
			char outstring[512];
			snprintf(outstring, sizeof(outstring), "%s%s%s%04lu%s.pdb", out_path_.c_str(), out_prefix_.c_str(), out_filename_.c_str(), static_cast<unsigned long>(irepeat), out_suffix_.c_str() );

			if ( utility::file::file_exists(std::string(outstring)) ) {
				TR << "nstruct " <<irepeat << " exists.  Continueing." << std::endl;
				continue;
			}
		}

		//Cyclic permutation of sequence.
		core::Size cyclic_offset(0);
		utility::vector1 < std::string > resnames_copy;
		if ( cyclic_permutations_ ) {
			cyclic_offset = do_cyclic_permutation( resnames, resnames_copy );
		} else {
			resnames_copy = resnames;
		}
		runtime_assert(cyclic_offset < resnames_copy.size() ); //Should be true.

		//Create the pose:
		core::pose::PoseOP pose( utility::pointer::make_shared< core::pose::Pose >() );
		build_polymer(pose, resnames_copy);

		core::Size thioether_sidechain_index(0);

		if ( cyclization_type() == SCPA_terminal_disulfide ) {
			set_up_terminal_disulfide_variants( pose ); //Add disulfide variants, if we're doing disulfide cyclization.
		} else if ( cyclization_type() == SCPA_thioether_lariat ) {
			thioether_sidechain_index = set_up_terminal_thioether_lariat_variants( pose ); //Add thioether lariat variant to the N-terminus, if we're doing that type of cyclization.
		} else if ( cyclization_type() == SCPA_lanthipeptide ) {
			//need to call for all lan and melan resi pairs
			if ( lanthionine_positions_.size() > 0 ) {
				for ( core::Size i=1, imax=lanthionine_positions_.size(); i<=imax; ++i ) { //Loop through all sets of doubles of residues.
					debug_assert(lanthionine_positions_[i].size() == 2); //Should always be true.
					protocols::cyclic_peptide::crosslinker::set_up_lanthionine_variants(*pose, lanthionine_positions_[i][1], lanthionine_positions_[i][2]);
				}
				core::pose::add_lower_terminus_type_to_pose_residue( *pose, 1 );
				core::pose::add_upper_terminus_type_to_pose_residue( *pose, pose->size() );
			} else {
				runtime_assert_string_msg( lanthionine_positions_.size() > 0 , "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication lanthionine positions are empty");
			}
		} else if ( cyclization_type() == SCPA_cterm_isopeptide_lariat || cyclization_type() == SCPA_nterm_isopeptide_lariat || cyclization_type() == SCPA_sidechain_isopeptide ) {
			set_up_isopeptide_variants( pose );
		}

		//Mover to cyclize the polymer and to update terminal peptide bond O and H atoms:
		protocols::simple_moves::DeclareBondOP termini( utility::pointer::make_shared< protocols::simple_moves::DeclareBond >() );
		set_up_cyclization_mover( termini, pose ); //Handles the cyclization appropriately, contingent on the cyclization type.
		termini->apply(*pose);
		if ( cyclization_type() == SCPA_thioether_lariat ) {
			runtime_assert( thioether_sidechain_index != 0  );
			protocols::cyclic_peptide::crosslinker::correct_thioether_virtuals( *pose, 1, thioether_sidechain_index );
		} else if ( cyclization_type() == SCPA_lanthipeptide ) {
			if ( lanthionine_positions_.size() > 0 ) {
				for ( core::Size i=1, imax=lanthionine_positions_.size(); i<=imax; ++i ) { //Loop through all sets of doubles of residues.
					debug_assert(lanthionine_positions_[i].size() == 2); //Should always be true.
					protocols::cyclic_peptide::crosslinker::correct_lanthionine_virtuals( *pose, lanthionine_positions_[i][1], lanthionine_positions_[i][2]);
				}
			} else {
				runtime_assert_string_msg( lanthionine_positions_.size() > 0 , "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication lanthionine positions are empty");
			}
		}

		//Add cyclic constraints:
		if ( cyclization_type() == SCPA_n_to_c_amide_bond  && use_chainbreak_energy() ) {
			add_cutpoint_variants_at_termini(pose);
		}

		//Add cyclic constraints for appropriate cyclization type:
		add_cyclic_constraints(pose);

		//Set all omega values to 180 and randomize mainchain torsions:
		set_mainchain_torsions(pose, cyclic_offset);

		//Add N-methylation:
		add_n_methylation( pose, cyclic_offset );

#ifdef BOINC_GRAPHICS
		// attach boinc graphics pose observer
		protocols::boinc::Boinc::attach_graphics_current_pose_observer( *pose );
		protocols::boinc::Boinc::attach_graphics_current_pose_ghost_observer( *pose );
#endif

#ifdef GL_GRAPHICS
		core::Vector center_vector = core::pose::get_center_of_mass( *pose );
		protocols::viewer::add_conformation_viewer( pose->conformation(), "current", 500, 500, false, true, center_vector );
#endif
		//Do the kinematic closure:
		bool const success( genkic_close(pose, sfxn_highhbond_cst, sfxn_highhbond_cst_cart, sfxn_default, total_hbond, cyclic_offset) );

#ifdef BOINC_GRAPHICS
		// attach boinc graphics pose observer
		protocols::boinc::Boinc::update_graphics_current( *pose );
		protocols::boinc::Boinc::update_graphics_current_ghost( *pose );
#endif

		if ( !success ) {
			TR << "Closure failed.";
			if ( irepeat < irepeat_max ) {
				TR << "  Continuing to next job." << std::endl;
				checkpoint( irepeat, success_count ); //This job has been attempted and has failed; don't repeat it.
#ifdef BOINC
				//Increment total jobs and check whether it's time to quit.
				if ((success_count || irepeat > max_attempts) && protocols::boinc::Boinc::worker_is_finished( success_count )) break;
#endif
			} else {
				TR << std::endl;
			}
			TR.flush();
			continue;
		}

		//If we reach here, then closure was successful.  Time to relax the pose.

		TR << "Closure successful." << std::endl;

		if ( fast_relax_rounds_ > 0 ) do_final_fastrelax( pose, sfxn_default_cst, fast_relax_rounds_, false, false, false );
		if ( angle_relax_rounds() > 0 ) do_final_fastrelax( pose, sfxn_default_cst_cart, angle_relax_rounds(), true, false, false );
		if ( angle_length_relax_rounds() > 0 ) do_final_fastrelax( pose, sfxn_default_cst_cart, angle_length_relax_rounds(), true, true, false );
		if ( cartesian_relax_rounds() > 0 ) do_final_fastrelax( pose, sfxn_default_cst_cart, cartesian_relax_rounds(), false, false, true );
		if ( angle_relax_rounds() > 0 || angle_length_relax_rounds() > 0 || cartesian_relax_rounds() > 0 ) {
			do_final_fastrelax( pose, sfxn_default_cst, 1, false, false, false ); //Do one more round of regular FastRelax if we've done any Cartesian, just to make sure we're in a pro_close minimum.
		}

		//If we're filtering by symmetry, do so here a final time:
		if ( required_symmetry_repeats_ > 1 ) {
			protocols::cyclic_peptide::CycpepSymmetryFilterOP symmfilter3( utility::pointer::make_shared< protocols::cyclic_peptide::CycpepSymmetryFilter >() );
			symmfilter3->set_symm_repeats( required_symmetry_repeats_ );
			symmfilter3->set_mirror_symm( required_symmetry_mirroring_ );
			symmfilter3->set_angle_threshold( required_symmetry_angle_threshold_ );
			core::select::residue_selector::ResidueIndexSelectorOP iselector( utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >() );
			std::stringstream pep_indices("");
			pep_indices << "1-" << sequence_length();
			iselector->set_index( pep_indices.str() );
			symmfilter3->set_selector( iselector );
			if ( symmfilter3->apply( *pose ) ) {
				TR << "Final symmetry filter passes.  This peptide has " << ( required_symmetry_mirroring_ ? "s" : "c" ) << required_symmetry_repeats_ << "symmetry." << std::endl;
			} else {
				TR << "Final symmetry filter failed.  This peptide lost " << ( required_symmetry_mirroring_ ? "s" : "c") << required_symmetry_repeats_ << "symmetry during the final relaxation." << std::endl;
				if ( irepeat < irepeat_max ) {
					TR << "Continuing to next job." << std::endl;
					checkpoint( irepeat, success_count ) ;
#ifdef BOINC
					//Increment total jobs and check whether it's time to quit.
					if ((success_count || irepeat > max_attempts) && protocols::boinc::Boinc::worker_is_finished( success_count )) break;
#endif
				}
				continue;
			}
		}

		//Undo the cyclic permutation in anticipation of re-aligning to the native:
		if ( cyclic_permutations_ ) {
			depermute( pose, cyclic_offset );
		}

#ifdef BOINC_GRAPHICS
		// attach boinc graphics pose observer
		protocols::boinc::Boinc::update_graphics_current( *pose );
#endif

		core::Real native_rmsd(0.0);

		//Score the pose before output:
		(*sfxn_default)(*pose);

		//Filter based on total energy:
		if ( use_total_energy_cutoff_ && pose->energies().total_energy() > total_energy_cutoff_ ) {
			TR << "Total final pose energy is " << pose->energies().total_energy() << ", which is greater than the cutoff of " << total_energy_cutoff_ << ".  Failing job." << std::endl;
			TR.flush();
			checkpoint( irepeat, success_count ); //This job has been attempted and has failed; don't repeat it.
#ifdef BOINC
			//Increment total jobs and check whether it's time to quit.
			if ((success_count || irepeat > max_attempts) && protocols::boinc::Boinc::worker_is_finished( success_count )) break;
#endif
			continue;
		}

		//Re-filter based on number of Hbonds (using option[min_final_hbonds]()):
		core::Real const final_hbonds( total_hbond->report_sm( *pose ) );
		if ( final_hbonds < static_cast<core::Real>(min_final_hbonds_) ) {
			TR << "Final hbond count is " << final_hbonds << ", which is less than the minimum.  Failing job." << std::endl;
			TR.flush();
			checkpoint( irepeat, success_count ); //This job has been attempted and has failed; don't repeat it.
#ifdef BOINC
			//Increment total jobs and check whether it's time to quit.
			if ((success_count || irepeat > max_attempts) && protocols::boinc::Boinc::worker_is_finished( success_count )) break;
#endif
			continue;
		}

		++success_count; //Increment the count of number of successes.

		if ( native_pose ) {
			native_rmsd = align_and_calculate_rmsd(*pose, *native_pose);
		}

		core::Size const cis_peptide_bonds( count_cis_peptide_bonds( pose ) );

		TR << "Result\tRMSD\tEnergy\tHbonds\tCisPepBonds" << std::endl;
		TR << irepeat << "\t";
		if ( native_pose ) { TR << native_rmsd; }
		else { TR << "--"; }
		TR << "\t" << pose->energies().total_energy() << "\t" << static_cast<core::Size>( std::round( final_hbonds) ) << "\t" << cis_peptide_bonds << std::endl;

		if ( silent_out_ || silentlist_out_ ) { //Writing directly to silent file or to a list of silent file data OPs
			core::io::silent::SilentFileOptions opts;
			opts.in_fullatom(true);
			core::io::silent::SilentStructOP ss( core::io::silent::SilentStructFactory::get_instance()->get_silent_struct("binary", opts ) );
			char tag[512];
			if ( my_rank_ > 0 ) {
				snprintf(tag, sizeof(tag), "result_proc%04lu_%04lu", static_cast<unsigned long>(my_rank_), static_cast<unsigned long>(irepeat+already_completed_job_count_) );
			} else {
				snprintf(tag, sizeof(tag), "result_%04lu", static_cast<unsigned long>(irepeat) );
			}
			ss->fill_struct( *pose, std::string(tag) );
			if ( native_pose ) ss->add_energy( "RMSD", native_rmsd ); //Add the RMSD to the energy to be written out in the silent file.
			ss->add_energy( "HBOND_COUNT", final_hbonds ); //Add the hbond count to be written out in the silent file.
			ss->add_energy( "CIS_PEPTIDE_BOND_COUNT", cis_peptide_bonds ); //Add the cis-peptide bond count to be written out in the silent file.
#ifdef BOINC_GRAPHICS
			protocols::boinc::Boinc::update_graphics_current( *pose );
			protocols::boinc::Boinc::update_graphics_current_ghost( *pose );
			protocols::boinc::Boinc::update_graphics_last_accepted( *pose, pose->energies().total_energy() );
			protocols::boinc::Boinc::update_graphics_low_energy( *pose, pose->energies().total_energy() );
#endif
			if ( silent_out_ ) {
				core::io::silent::SilentFileDataOP silent_file (new core::io::silent::SilentFileData( opts ) );
				silent_file->set_filename( out_filename_ );
				silent_file->write_silent_struct( *ss, out_filename_ );
			}
			if ( silentlist_out_ ) {
				silentlist_->push_back(ss);
				core::Size curjob( summarylist_->size() + 1 );
				summarylist_->push_back( utility::pointer::make_shared< HierarchicalHybridJD_JobResultsSummary >( my_rank_, curjob, pose->energies().total_energy(), (native_pose ? native_rmsd : 0), static_cast< core::Size >( std::round(final_hbonds) ), cis_peptide_bonds ) );
			}
		} else { //if pdb output

			char outstring[512];
			snprintf(outstring, sizeof(outstring), "%s%s%s%04lu%s.pdb", out_path_.c_str(), out_prefix_.c_str(), out_filename_.c_str(), static_cast<unsigned long>(irepeat), out_suffix_.c_str() );

			pose->dump_scored_pdb( std::string(outstring), *sfxn_default );
		}

		TR.flush();
		checkpoint( irepeat, success_count ); //This job has been attempted and has succeeded; don't repeat it.

#ifdef BOINC
		//Increment total jobs and check whether it's time to quit.
		if (protocols::boinc::Boinc::worker_is_finished( success_count )) break;
#endif

#ifdef GL_GRAPHICS
		protocols::viewer::clear_conformation_viewers();
#endif
	} //Looping through nstruct

	TR << nstruct_ << " jobs attempted.  " << success_count << " jobs returned solutions." << std::endl;
	TR.flush();

	end_checkpointing(); //Delete the checkpoint file at this point, since all jobs have completed.
	return;
}

/// @brief Is a residue type supported for macrocycle structure prediction?
/// @details Currently returns true for alpha-, beta-, or gamma-amino acids and for peptoids,
/// false otherwise.  This will be expanded in the future.
bool
SimpleCycpepPredictApplication::is_supported_restype(
	core::chemical::ResidueType const & restype
) const {
	return ( restype.is_alpha_aa() || restype.is_beta_aa() || restype.is_gamma_aa() || restype.is_peptoid() || restype.is_oligourea() );
}

/// @brief Is this residue to be ignored in calculating RMSDs?
bool
SimpleCycpepPredictApplication::is_residue_ignored_in_rms(
	core::Size const res_index
) const {
	if ( exclude_residues_from_rms_.empty() ) return false;
	return exclude_residues_from_rms_.has_value( res_index );
}

///  @brief Is a given cyclization type a lariat type (i.e. one where a side-chain connects to backbone)?
bool
SimpleCycpepPredictApplication::is_lariat_type(
	SCPA_cyclization_type const type_in
) const {
	if ( type_in == SCPA_nterm_isopeptide_lariat || type_in == SCPA_cterm_isopeptide_lariat || type_in == SCPA_thioether_lariat ) return true;
	return false;
}

/// @brief Can a position's backbone be randomized?
/// @details Returns false for disulfide or isopeptide positions, true otherwise.
bool
SimpleCycpepPredictApplication::position_backbone_is_randomizable(
	core::Size const res_index
) const {
	if ( cyclization_type() == SCPA_terminal_disulfide || cyclization_type() == SCPA_sidechain_isopeptide ) {
		if ( res_index == 1 || res_index == sequence_length() ) return false;
	} else if ( cyclization_type() == SCPA_nterm_isopeptide_lariat && res_index == sequence_length() ) {
		return false;
	} else if ( cyclization_type() == SCPA_cterm_isopeptide_lariat && res_index == 1 ) {
		return false;
	} else if ( cyclization_type() == SCPA_lanthipeptide ) {
		return false;
	}
	return true;
}

/// @brief Check that the loop formed is long enough.
void
SimpleCycpepPredictApplication::check_loop_length(
	utility::vector1< std::string > const &resnames
) const {
	if ( cyclization_type() == SCPA_terminal_disulfide ) {
		//Create a temporary pose for analysis.  We need this because the check is based
		//on disulfide-forming residue types.  This is slightly inefficient, but it only
		//happens once:
		core::pose::PoseOP temp_pose( utility::pointer::make_shared< core::pose::Pose >() );
		build_polymer(temp_pose, resnames);
		core::Size const firstdisulf( find_first_disulf_res( temp_pose  ) );
		core::Size const lastdisulf( find_last_disulf_res( temp_pose  ) );
		runtime_assert_string_msg(
			lastdisulf > firstdisulf + 3,
			"Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::check_loop_length(): The sequence provided has terminal disulfide-forming residues that are too close together.  There must be at least three residues between the terminal disulfide-forming residues to allow cyclization."
		);
	} else if ( cyclization_type() == SCPA_nterm_isopeptide_lariat || cyclization_type() == SCPA_cterm_isopeptide_lariat || cyclization_type() == SCPA_sidechain_isopeptide ) {
		//Create a temporary pose for analysis.
		core::pose::PoseOP temp_pose( utility::pointer::make_shared< core::pose::Pose >() );
		build_polymer(temp_pose, resnames);
		core::Size firstres, lastres;
		find_first_and_last_isopeptide_residues( temp_pose, firstres, lastres );
		runtime_assert_string_msg( lastres - firstres > 3, "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::check_loop_length(): The sequence provided has isopeptide bond-forming residues that are too close together.  There must be at least three residues between the isopeptide bond-forming residues to allow cyclization." );
	} else if ( cyclization_type() == SCPA_thioether_lariat ) { //Check that the thioether lariat loop is long enough
		core::pose::PoseOP temp_pose( utility::pointer::make_shared< core::pose::Pose >() );
		build_polymer(temp_pose, resnames);
		core::Size firstres, lastres;
		find_first_and_last_thioether_lariat_residues( temp_pose, firstres, lastres );
		runtime_assert_string_msg( lastres - firstres > 3, "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::check_loop_length(): The sequence provided has thioether lariat bond-forming residues that are too close together.  There must be at least three residues between the isopeptide bond-forming residues to allow cyclization." );
	}
}

/// @brief Given a pose, find the first and last isopeptide bond-forming residues.
/// @details Bases this on lariat_sidechain_index_ for lariat types, unless set to 0, in which case it finds the
/// suitable type closest to the opposite terminus.  Bases this on sidechain_isopeptide_indices_ for sidechain isopeptide
/// cyclization, unless set to 0, in which case it finds the suitable types that are furthest apart.
void
SimpleCycpepPredictApplication::find_first_and_last_isopeptide_residues(
	core::pose::PoseCOP pose,
	core::Size &firstres,
	core::Size &lastres
) const {
	using namespace core::chemical;

	runtime_assert( cyclization_type() == SCPA_cterm_isopeptide_lariat || cyclization_type() == SCPA_nterm_isopeptide_lariat || cyclization_type() == SCPA_sidechain_isopeptide );
	if ( cyclization_type() == SCPA_cterm_isopeptide_lariat ) {
		lastres = sequence_length();
		if ( lariat_sidechain_index_ != 0 ) {
			firstres = lariat_sidechain_index_;
			return;
		}
		//If we reach here, we need to search for the first amine residue:
		for ( core::Size ir(1), irmax(pose->total_residue()); ir < irmax; ++ir ) {
			std::string const & curbasename( pose->residue_type(ir).base_name() );
			if ( is_isopeptide_forming_amide_type(curbasename) ) {
				firstres = ir;
				return;
			}
		}
		//If we reach here, we've failed to find a carboxyl-containing residue:
		utility_exit_with_message("Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::find_first_and_last_isopeptide_residues(): No amine-containing sidechain was found that could form an isopeptide bond with the peptide C-terminus!");
	} else if ( cyclization_type() == SCPA_nterm_isopeptide_lariat ) {
		firstres = 1;
		if ( lariat_sidechain_index_ != 0 ) {
			lastres = lariat_sidechain_index_;
			return;
		}
		//If we reach here, we need to search for the last carboxyl residue:
		for ( core::Size ir(pose->total_residue()); ir > 1; --ir ) {
			AA const curaa( pose->residue(ir).aa() );
			if ( is_isopeptide_forming_carbonyl_type(curaa) ) {
				lastres = ir;
				return;
			}
		}
		//If we reach here, we've failed to find a carboxyl-containing residue:
		utility_exit_with_message("Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::find_first_and_last_isopeptide_residues(): No carboxyl-containing sidechain was found that could form an isopeptide bond with the peptide N-terminus!");
	} else if ( cyclization_type() == SCPA_sidechain_isopeptide ) {
		if ( sidechain_isopeptide_indices_.first != 0 && sidechain_isopeptide_indices_.second != 0 ) {
			firstres = sidechain_isopeptide_indices_.first;
			lastres = sidechain_isopeptide_indices_.second;
			if ( firstres > lastres ) {
				core::Size const temp( firstres );
				firstres = lastres;
				lastres = temp;
			}
			return;
		}
		//If we reach here, we need to search.  Two possibilities: nitrogenous aa first, carboxyl aa second, or the converse.
		core::Size first_nitrog(0), last_nitrog(0), first_carbox(0), last_carbox(0);
		for ( core::Size ir(1), irmax(pose->total_residue()); ir<=irmax; ++ir ) {
			std::string const & curbasename( pose->residue_type(ir).base_name() );
			AA const curaa( pose->residue_type(ir).aa() );
			if ( is_isopeptide_forming_amide_type(curbasename) ) {
				if ( !first_nitrog ) {
					first_nitrog = ir;
				}
				last_nitrog = ir;
			} else if ( is_isopeptide_forming_carbonyl_type(curaa) ) {
				if ( !first_carbox ) {
					first_carbox = ir;
				}
				last_carbox = ir;
			}
		}
		runtime_assert_string_msg( first_carbox && first_nitrog, "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::find_first_and_last_isopeptide_residues(): Two residues able to form an isopeptide bond through their side-chains were not found!" );
		if ( (static_cast<signed int>(last_nitrog) - static_cast<signed int>(first_carbox)) > (static_cast<signed int>(last_carbox) - static_cast<signed int>(first_nitrog)) ) {
			lastres = last_nitrog;
			firstres = first_carbox;
			return;
		}
		lastres = last_carbox;
		firstres = first_nitrog;
		return;
	}
	utility_exit_with_message( "PROGRAM ERROR.  IT SHOULD NOT BE POSSIBLE TO REACH THIS POINT." );
}

/// @brief Given a pose, find the first and last thioether lariat bond-forming residues.
/// @details First residue is 1 by definition (chloroacetyl goes where??); last is
/// TYPICALLY C-term in the classic peptidream approach but does not have to be.
/// n.b. Suga has methods for incorporating additional cysteines into bicyclic
/// peptides, but for the moment "closest to the opposite terminus" is sufficient.
void
SimpleCycpepPredictApplication::find_first_and_last_thioether_lariat_residues(
	core::pose::PoseCOP pose,
	core::Size &firstres,
	core::Size &lastres
) const {
	using namespace core::chemical;

	runtime_assert( cyclization_type() == SCPA_thioether_lariat );

	firstres = 1;
	if ( lariat_sidechain_index_ != 0 ) {
		lastres = lariat_sidechain_index_;
	} else {
		lastres = find_last_disulf_res( pose );
	}
}

/// @brief Given a pose, find the first and last lanthipeptide bond-forming residues.
/// @details Randomly choose based on input residue numbers
void
SimpleCycpepPredictApplication::find_first_and_last_lanthipeptide_residues(
	/*core::pose::PoseCOP pose,*/
	core::Size &firstres,
	core::Size &lastres
) const {
	using namespace core::chemical;

	runtime_assert( cyclization_type() == SCPA_lanthipeptide );

	//only lanthionine rings rn
	if ( lanthionine_positions_.size() > 0 ) {
		core::Size  lan_ring( numeric::random::rg().random_range(1, lanthionine_positions_.size()) );
		if ( lanthionine_positions_[lan_ring][1] < lanthionine_positions_[lan_ring][2] ) {
			firstres = lanthionine_positions_[lan_ring][1];
			lastres = lanthionine_positions_[lan_ring][2];
		} else {
			firstres = lanthionine_positions_[lan_ring][2];
			lastres = lanthionine_positions_[lan_ring][1];
		}
	} else {
		runtime_assert_string_msg( lanthionine_positions_.size() > 0, "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication only lanthionine rn");
	}
	//firstres = 1;
	//if ( lariat_sidechain_index_ != 0 ) {
	// lastres = lariat_sidechain_index_;
	//} else {
	// lastres = find_last_disulf_res( pose );
	//}
}

/// @brief Count the number of cis-peptide bonds in the pose.
/// @details Counts as cis if in the range (-90,90].
core::Size
SimpleCycpepPredictApplication::count_cis_peptide_bonds(
	core::pose::PoseCOP pose
) const {
	core::Size count(0);
	for ( core::Size i=1, imax=sequence_length(); i<=imax; ++i ) {
		core::Real const omegaval( numeric::principal_angle_degrees( pose->omega(i) /*Should handle terminal peptide bonds.*/ ) );
		TR.Debug << "omega" << i << "=" << omegaval << std::endl;
		if ( omegaval <= 90 && omegaval > -90 ) ++count; //Count this as cis if in the interval (-90,90]
	}
	return count;
}

/// @brief Carry out the final FastRelax.
///
void
SimpleCycpepPredictApplication::do_final_fastrelax(
	core::pose::PoseOP pose,
	core::scoring::ScoreFunctionOP sfxn,
	core::Size const relax_rounds,
	bool const angle_min,
	bool const length_min,
	bool const cartesian_min
) const {
	protocols::relax::FastRelaxOP frlx( utility::pointer::make_shared< protocols::relax::FastRelax >(sfxn, 1) );
	if ( cartesian_min ) {
		frlx->cartesian(true);
	} else {
		if ( angle_min ) {
			frlx->minimize_bond_angles(true);
		}
		if ( length_min ) {
			frlx->minimize_bond_lengths(true);
		}
	}

	//Mover to update terminal peptide bond O and H atoms:
	protocols::simple_moves::DeclareBondOP final_termini( utility::pointer::make_shared< protocols::simple_moves::DeclareBond >() );
	set_up_cyclization_mover( final_termini, pose );

	(*sfxn)(*pose);
	core::Real cur_energy( pose->energies().total_energy() );
	for ( core::Size i=1; i<=relax_rounds; ++i ) {
		core::pose::PoseOP pose_copy( pose->clone() );
		if ( TR.visible() ) {
			TR << "Applying final FastRelax, round " << i;
			if ( cartesian_min ) TR << ", with Cartesian minimization.";
			else if ( angle_min && length_min ) TR << ", with bond angle and bond length minimization.";
			else if ( angle_min && !length_min ) TR << ", with bond angle minimization.";
			TR << "." << std::endl;
		}
		frlx->apply( *pose_copy );
		final_termini->apply( *pose_copy );
		(*sfxn)(*pose_copy);
		if ( pose_copy->energies().total_energy() < cur_energy ) {
			cur_energy = pose_copy->energies().total_energy();
			(*pose) = (*pose_copy);
		}
#ifdef BOINC_GRAPHICS
		// attach boinc graphics pose observer
		protocols::boinc::Boinc::update_graphics_current( *pose );
#endif
	}
}

/// @brief Actually build the geometry that we'll be working with.
///
void
SimpleCycpepPredictApplication::build_polymer(
	core::pose::PoseOP pose,
	utility::vector1<std::string> const &restypes
) const {
	using namespace protocols::cyclic_peptide;
	core::Size const nres( restypes.size() );
	runtime_assert(restypes.size() >=4 );

	TR << "Building sequence ";
	for ( core::Size i=1; i<=nres; ++i ) {
		TR << restypes[i];
		if ( i<nres ) TR << " ";
	}
	TR << "." << std::endl;

	PeptideStubMover stubmover;

	stubmover.set_reset_mode(true);
	stubmover.reset_mover_data();
	for ( core::Size i=1; i<=nres; ++i ) {
		stubmover.add_residue( "Append", restypes[i], 0, false, "", 1, 0, nullptr, "" );
	}

	stubmover.apply(*pose);

	TR << "Build successful." << std::endl;

	return;
} //build_polymer()

/// @brief Add N-methylation.
/// @details Must be called after pose is cyclized.
void
SimpleCycpepPredictApplication::add_n_methylation(
	core::pose::PoseOP pose,
	core::Size const cyclic_offset
) const {
	if ( n_methyl_positions_.size() == 0 ) { return; } //Do nothing if there's no N-methylation.

	core::Size const nres(sequence_length());

	TR << "Adding N-methylation" << std::endl;
	protocols::simple_moves::ModifyVariantTypeMover add_nmethyl;
	add_nmethyl.set_additional_type_to_add("N_METHYLATION");
	core::select::residue_selector::ResidueIndexSelectorOP selector( utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >() );
	for ( core::Size i=1, imax=n_methyl_positions_.size(); i<=imax; ++i ) {
		runtime_assert_string_msg( n_methyl_positions_[i] > 0  && n_methyl_positions_[i] <= nres, "Error in simple_cycpep_predict app: The N-methylation position indices must be within the pose!" );
		int permuted_position( static_cast<int>(n_methyl_positions_[i]) - static_cast<int>(cyclic_offset) );
		if ( permuted_position < 1 ) permuted_position += static_cast<int>(nres);
		selector->append_index( static_cast<core::Size>(permuted_position) );
	}
	add_nmethyl.set_residue_selector(selector);
	add_nmethyl.apply(*pose);

}

/// @brief Given the name of a Rama_Table_Type, set the default Rama_Table_Type.
/// @details Error if unknown type.
void
SimpleCycpepPredictApplication::set_default_rama_table_type(
	std::string const &type_name
) {
	if ( type_name=="" ) return; //Default case -- no default rama table provided.

	default_rama_table_type_ = get_rama_table_type_from_name( type_name );
	if ( TR.visible() ) {
		TR << "Set default Rama table type to " << type_name << "." << std::endl;
		TR.flush();
	}
	return;
}

/// @brief Given a string vector that we need to parse, populate the rama_table_type_by_res_ map.
/// @details The string vector must be of the format: [integer] [rama type name] [integer] [rama type name] etc.
/// Throws error if could not parse.
void
SimpleCycpepPredictApplication::set_rama_table_type_by_res(
	utility::vector1 <std::string> const &type_name_vector
) {
	core::Size const vectsize( type_name_vector.size() );
	if ( vectsize == 0 ) return;

	core::Size i(0);
	while ( i<vectsize ) {
		++i;
		std::stringstream ss(type_name_vector[i]);
		core::Size tempval(0);
		ss >> tempval;
		if ( ss.fail() ) {
			std::stringstream msg("");
			msg << "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::set_rama_table_type_by_res(): Could not interpret " << type_name_vector[i] << " as an integer.  The \"-rama_sampling_table_by_res\" flag must be given a series of values of the pattern [res_index] [rama_table_type] [res_index] [rama_table_type] etc." << std::endl;
			utility_exit_with_message(msg.str());
		}
		++i;
		if ( i > vectsize ) {
			std::stringstream msg("");
			msg << "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::set_rama_table_type_by_res(): A residue index was found with no corresponding Rama table type." << std::endl;
			utility_exit_with_message(msg.str());
		}
		std::string tempstring("");
		std::stringstream ss2(type_name_vector[i]);
		ss2 >> tempstring;
		rama_table_type_by_res_[tempval] = get_rama_table_type_from_name(tempstring);
		if ( TR.visible() ) TR << "Set custom Ramachandran table for residue " << tempval << " to " << tempstring << "." << std::endl;
	}

	if ( TR.visible() ) TR.flush();

	return;
}

/// @brief Given a Rama_Table_Type name, return the Rama_Table_Type, or an informative error message on failure.
///
core::scoring::Rama_Table_Type
SimpleCycpepPredictApplication::get_rama_table_type_from_name(
	std::string const &type_name
) const {
	//Get an instance of the Ramachandran object:
	core::scoring::Ramachandran const & rama = core::scoring::ScoringManager::get_instance()->get_Ramachandran();

	//Get the type from the name:
	core::scoring::Rama_Table_Type const type ( rama.get_ramatable_type_by_name(type_name) );

	//Check that the type was parsed sensibly:
	if ( type == core::scoring::unknown_ramatable_type ) {
		std::stringstream err_msg("Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::get_rama_table_type_from_name():");
		err_msg << "  The provided custom Ramachandran table type (\"" << type_name << "\") is unknown.  Allowed options are: ";
		for ( core::Size i=1; i<static_cast<core::Size>(core::scoring::unknown_ramatable_type); ++i ) {
			err_msg << rama.get_ramatable_name_by_type( static_cast<core::scoring::Rama_Table_Type>(i) );
			if ( i < (static_cast<core::Size>(core::scoring::unknown_ramatable_type) - 1) ) {
				err_msg << ", ";
			}
		}
		err_msg << "." << std::endl;
		utility_exit_with_message( err_msg.str() );
	}

	return(type);
}


/// @brief Read a sequence (as a series of full names, separated by whitespace) and store
/// it in a string vector.
void
SimpleCycpepPredictApplication::read_sequence (
	std::string const &seqfile,
	utility::vector1 < std::string > &resnames
) const {
	using namespace utility::io;
	resnames.clear();

	utility::vector1< std::string > lines; //Storing all lines
	if ( sequence_string_ == "" ) {
		izstream infile;
		infile.open( seqfile );
		runtime_assert_string_msg( infile.good(), "Error in read_sequence() in app simple_cycpep_predict:  Unable to open sequence file for read!" );

		TR << "Opened " << seqfile << " for read." << std::endl;

		std::string curline(""); //Buffer for current line.

		//Read the file:
		while ( getline(infile, curline) ) {
			if ( curline.size() < 1 ) continue; //Ignore blank lines.
			lines.push_back( curline );
		}
		infile.close();
	} else {
		lines.push_back( sequence_string_ );
	}

	//Parse the lines:
	for ( core::Size i=1, imax=lines.size(); i<=imax; ++i ) { //Loop through all lines
		if ( TR.Debug.visible() ) TR.Debug << "Parsing \"" << lines[i] << "\"." << std::endl;
		std::istringstream curline(lines[i]);
		std::string oneword("");
		while ( !curline.eof() ) {
			curline >> oneword;
			resnames.push_back( oneword );
		}
	}

	if ( TR.visible() ) {
		TR << "Parsed the following sequence:" << std::endl;
		for ( core::Size i=1, imax=resnames.size(); i<=imax; ++i ) {
			TR << resnames[i];
			if ( i<imax ) TR << ", ";
		}
		TR << "." << std::endl;
	}

	runtime_assert_string_msg( resnames.size() >= 4, "Error in simple_cycpcp_predict app read_sequence() function!  The minimum number of residues for a cyclic peptide is 4.  (GenKIC requires three residues, plus a fourth to serve as an anchor)." );

	return;
}

/// @brief Set up the mover that creates N-to-C amide bonds, and which updates the
/// atoms dependent on the amide bond.
void
SimpleCycpepPredictApplication::set_up_n_to_c_cyclization_mover (
	protocols::simple_moves::DeclareBondOP termini,
	core::pose::PoseCOP pose,
	bool const native,
	core::Size const last_res
) const {
	core::Size const nres(sequence_length());

	core::Size const cterm( last_res == 0 ? nres : last_res );

	runtime_assert_string_msg(pose->residue(1).has_lower_connect(), "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictpplication::set_up_n_to_c_cyclization_mover() function: residue 1 does not have a LOWER_CONNECT.");
	runtime_assert_string_msg(pose->residue(cterm).has_upper_connect(), "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictpplication::set_up_n_to_c_cyclization_mover() function: the final residue does not have an UPPER_CONNECT.");
	std::string firstatom( pose->residue(1).atom_name( pose->residue(1).lower_connect_atom() ) );
	std::string lastatom( pose->residue(cterm).atom_name( pose->residue(cterm).upper_connect_atom() ) );

	if ( native ) {
		TR << "Setting up terminal bond for the native pose between residue 1, atom " << firstatom << " and residue " << cterm << ", atom " << lastatom << "." << std::endl;
	} else {
		TR << "Setting up terminal bond between residue 1, atom " << firstatom << " and residue " << cterm << ", atom " << lastatom << "." << std::endl;
	}

	termini->set( cterm, lastatom, 1, firstatom, false, false, 0, 0, false  );
}

/// @brief Set up the mover that creates terminal disulfide bonds.
///
void
SimpleCycpepPredictApplication::set_up_terminal_disulfide_cyclization_mover (
	protocols::simple_moves::DeclareBondOP termini,
	core::pose::PoseCOP pose,
	bool const native,
	core::Size const last_disulf_res,
	core::Size const first_disulf_res
) const {
	core::Size first_disulf( first_disulf_res == 0 ? find_first_disulf_res( pose ) : first_disulf_res );
	runtime_assert_string_msg( pose->residue_type(first_disulf).is_sidechain_thiol() || pose->residue_type(first_disulf).is_disulfide_bonded(), "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictpplication::set_up_terminal_disulfide_cyclization_mover():  The first residue that is supposed to form the disulfide doesn't contain a sidechain thiol." );

	core::Size last_disulf( last_disulf_res == 0 ? find_last_disulf_res( pose ) : last_disulf_res );
	runtime_assert_string_msg( pose->residue_type(last_disulf).is_sidechain_thiol() || pose->residue_type(last_disulf).is_disulfide_bonded(), "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictpplication::set_up_terminal_disulfide_cyclization_mover():  The second residue that is supposed to form the disulfide doesn't contain a sidechain thiol." );
	runtime_assert_string_msg( first_disulf != last_disulf, "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictpplication::set_up_terminal_disulfide_cyclization_mover(): The first and last disulfide residue were the same.  Are there at least two disulfide-forming types in this peptide?" );

	std::string const firstatom( pose->residue_type(first_disulf).get_disulfide_atom_name() );
	std::string const lastatom( pose->residue_type(last_disulf).get_disulfide_atom_name() );

	TR << "Setting up terminal disulfide bond" << (native ? " for the native pose " : " ") << "between residue " << first_disulf << " atom " << firstatom << " and residue " << last_disulf << " atom " << lastatom << "." << std::endl;

	termini->set( first_disulf, firstatom, last_disulf, lastatom, true, false, 0, 0, false );
}


/// @brief Set up the mover that creates N-terminal isopeptide bonds.
void
SimpleCycpepPredictApplication::set_up_nterm_isopeptide_cyclization_mover(
	protocols::simple_moves::DeclareBondOP termini,
	core::pose::PoseCOP pose
) const {
	debug_assert(termini); //Should already exist.
	debug_assert(cyclization_type() == SCPA_nterm_isopeptide_lariat); //Should be true.

	core::Size firstres, lastres;
	find_first_and_last_isopeptide_residues( pose, firstres, lastres );

	//Get the name of the sidechain connection atom:
	core::chemical::ResidueType const &restype( pose->residue_type(lastres) );
	core::Size const lastres_sc_connection_id( restype.n_possible_residue_connections() ); //We will assume that the highest-numbered residue connection is the sidechain connection ID.
	core::Size const lastres_sc_connection_atom_index( restype.residue_connect_atom_index( lastres_sc_connection_id ) );
	std::string const lastres_sc_connection_atom( restype.atom_name( lastres_sc_connection_atom_index ) );

	TR << "Setting up isopeptide bond from N-terminus to residue " << restype.name3() << lastres << "." << std::endl;

	termini->set( firstres, "N", lastres, lastres_sc_connection_atom, false );
}

/// @brief Set up the mover that creates C-terminal isopeptide bonds.
void
SimpleCycpepPredictApplication::set_up_cterm_isopeptide_cyclization_mover(
	protocols::simple_moves::DeclareBondOP termini,
	core::pose::PoseCOP pose
) const {
	debug_assert(termini); //Should already exist.
	debug_assert(cyclization_type() == SCPA_cterm_isopeptide_lariat); //Should be true.

	core::Size firstres, lastres;
	find_first_and_last_isopeptide_residues( pose, firstres, lastres );

	//Get the name of the sidechain connection atom:
	core::chemical::ResidueType const &restype( pose->residue_type(firstres) );
	core::Size const firstres_sc_connection_id( restype.n_possible_residue_connections() ); //We will assume that the highest-numbered residue connection is the sidechain connection ID.
	core::Size const firstres_sc_connection_atom_index( restype.residue_connect_atom_index( firstres_sc_connection_id ) );
	std::string const firstres_sc_connection_atom( restype.atom_name( firstres_sc_connection_atom_index ) );

	TR << "Setting up isopeptide bond from C-terminus to residue " << restype.name3() << firstres << "." << std::endl;

	termini->set( lastres, "C", firstres, firstres_sc_connection_atom, false );
}

/// @brief Set up the mover that creates sidechain isopeptide bonds.
void
SimpleCycpepPredictApplication::set_up_sidechain_isopeptide_cyclization_mover(
	protocols::simple_moves::DeclareBondOP termini,
	core::pose::PoseCOP pose
) const {
	debug_assert(termini); //Should already exist.
	debug_assert(cyclization_type() == SCPA_sidechain_isopeptide); //Should be true.

	core::Size firstres, lastres;
	find_first_and_last_isopeptide_residues( pose, firstres, lastres );

	//Get the name of the first sidechain connection atom:
	core::chemical::ResidueType const &restype( pose->residue_type(firstres) );
	core::Size const firstres_sc_connection_id( restype.n_possible_residue_connections() ); //We will assume that the highest-numbered residue connection is the sidechain connection ID.
	core::Size const firstres_sc_connection_atom_index( restype.residue_connect_atom_index( firstres_sc_connection_id ) );
	std::string const firstres_sc_connection_atom( restype.atom_name( firstres_sc_connection_atom_index ) );

	//Get the name of the second sidechain connection atom:
	core::chemical::ResidueType const &restype2( pose->residue_type(lastres) );
	core::Size const lastres_sc_connection_id( restype2.n_possible_residue_connections() ); //We will assume that the highest-numbered residue connection is the sidechain connection ID.
	core::Size const lastres_sc_connection_atom_index( restype2.residue_connect_atom_index( lastres_sc_connection_id ) );
	std::string const lastres_sc_connection_atom( restype2.atom_name( lastres_sc_connection_atom_index ) );

	TR << "Setting up isopeptide bond from " << restype.name3() << firstres << " to residue " << restype2.name3() << lastres << "." << std::endl;

	termini->set( lastres, lastres_sc_connection_atom, firstres, firstres_sc_connection_atom, false );
}

/// @brief Set up the mover that creates thioether lariat bonds.
void
SimpleCycpepPredictApplication::set_up_thioether_lariat_cyclization_mover(
	protocols::simple_moves::DeclareBondOP termini,
	core::pose::PoseCOP pose
) const {
	runtime_assert( termini != nullptr ); //Should already exist.
	debug_assert( cyclization_type() == SCPA_thioether_lariat ); //Should be true.

	core::Size firstres, lastres;
	find_first_and_last_thioether_lariat_residues( pose, firstres, lastres );

	protocols::cyclic_peptide::crosslinker::set_up_thioether_bond_mover(
		*termini, *pose, firstres, lastres
	);
}

void
SimpleCycpepPredictApplication::set_up_lanthipeptide_cyclization_mover(
	protocols::simple_moves::DeclareBondOP termini,
	core::pose::PoseCOP pose
) const {
	runtime_assert( termini != nullptr ); //Should already exist.
	debug_assert( cyclization_type() == SCPA_lanthipeptide ); //Should be true.

	if ( lanthionine_positions_.size() > 0 ) {
		for ( core::Size i=1, imax=lanthionine_positions_.size(); i<=imax; ++i ) { //Loop through all sets of doubles of residues.
			debug_assert(lanthionine_positions_[i].size() == 2); //Should always be true.
			protocols::cyclic_peptide::crosslinker::set_up_lanthionine_bond_mover(*termini, *pose, lanthionine_positions_[i][1], lanthionine_positions_[i][2]);
		}
		//termini->set( 1, "C", 2, "N", true );
	} else {
		runtime_assert_string_msg( lanthionine_positions_.size() > 0 , "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication lanthionine positions are empty");
	}
}

/// @brief Set up the DeclareBond mover used to connect the termini, or whatever
/// atoms are involved in the cyclization.  (Handles different cyclization modes).
void
SimpleCycpepPredictApplication::set_up_cyclization_mover (
	protocols::simple_moves::DeclareBondOP termini,
	core::pose::PoseCOP pose,
	bool const native, /*=false*/
	core::Size const last_res, /*=0*/
	core::Size const /*first_res*/ /*=0*/
) const {
	switch( cyclization_type() ) {
	case SCPA_n_to_c_amide_bond :
		set_up_n_to_c_cyclization_mover( termini, pose, native, last_res );
		return;
	case SCPA_terminal_disulfide :
		set_up_terminal_disulfide_cyclization_mover( termini, pose, native, 0, 0 );
		return;
	case SCPA_nterm_isopeptide_lariat :
		set_up_nterm_isopeptide_cyclization_mover( termini, pose );
		return;
	case SCPA_cterm_isopeptide_lariat :
		set_up_cterm_isopeptide_cyclization_mover( termini, pose );
		return;
	case SCPA_sidechain_isopeptide :
		set_up_sidechain_isopeptide_cyclization_mover( termini, pose );
		return;
	case SCPA_thioether_lariat :
		set_up_thioether_lariat_cyclization_mover( termini, pose );
		return;
	case SCPA_lanthipeptide :
		set_up_lanthipeptide_cyclization_mover( termini, pose );
		return;
	case SCPA_invalid_type :
		break;
	}

	utility_exit_with_message( "Error in protocols::cyclic_peptide::SimpleCycpepPredictApplication::set_up_cyclization_mover(): The specified cyclization type is invalid." );

	return;
}


/// @brief Takes a vector of residue names, chooses a random number for cyclic offset, and
/// does a cyclic permutation.
/// @details Returns the offset and stores the new string vector in resnames_copy.
core::Size
SimpleCycpepPredictApplication::do_cyclic_permutation (
	utility::vector1 <std::string> const &resnames,
	utility::vector1 <std::string> &resnames_copy
) const {
	core::Size const nname( resnames.size() );//Number of residue names
	auto const offset( static_cast<core::Size>(numeric::random::rg().random_range(0,resnames.size()-1)) );

	resnames_copy.clear();
	resnames_copy.resize(nname, "");
	core::Size counter(0);
	for ( core::Size i=offset+1; i<=nname; ++i ) {
		++counter;
		resnames_copy[counter] = resnames[i];
	}
	for ( core::Size i=1; i<=offset; ++i ) {
		++counter;
		resnames_copy[counter] = resnames[i];
	}

	TR << "Circularly shifted residue list by " << offset << ".  New list is: ";
	for ( core::Size i=1; i<=nname; ++i ) {
		TR << resnames_copy[i];
		if ( i<nname ) TR << ", ";
	}
	TR << std::endl;
	TR.flush();

	return offset;
}


/// @brief Imports the native pose and sets up a terminial peptide bond.
///
void
SimpleCycpepPredictApplication::import_and_set_up_native (
	std::string const &native_file,
	core::pose::PoseOP native_pose,
	core::Size const expected_residue_count
) const {
	core::import_pose::pose_from_file(*native_pose, native_file, core::import_pose::PDB_file);
	TR << "Improrting native structure from " << native_file << "." << std::endl;

	set_up_native( native_pose, expected_residue_count );

	return;
}

/// @brief Sets up a terminial peptide bond and does some checks.
///
void
SimpleCycpepPredictApplication::set_up_native (
	core::pose::PoseOP native_pose,
	core::Size const expected_residue_count
) const {
	// Count residues and find the last residue.
	core::Size last_res(0), res_count(0);
	for ( core::Size ir=1, irmax=native_pose->total_residue(); ir<=irmax; ++ir ) {
		if ( is_supported_restype( native_pose->residue_type(ir) ) ) {
			++res_count;
			last_res = ir;
		}
	}

	if ( expected_residue_count != 0 ) {
		runtime_assert_string_msg(
			res_count == expected_residue_count,
			"Error in simple_cycpep_predict app!  The imported native pose has a different number of residues than the sequence provided."
		);
	}

	core::Size thioether_sidechain_index(0); //Only used for thioethers.

	if ( cyclization_type() == SCPA_n_to_c_amide_bond || cyclization_type() == SCPA_nterm_isopeptide_lariat || cyclization_type() == SCPA_cterm_isopeptide_lariat ) {
		TR << "Stripping termini from native structure." << std::endl;
		if ( cyclization_type() == SCPA_n_to_c_amide_bond || cyclization_type() == SCPA_nterm_isopeptide_lariat ) {
			if ( native_pose->residue_type(1).is_lower_terminus() ) {
				core::pose::remove_lower_terminus_type_from_pose_residue(*native_pose, 1);
				TR << "Removed lower terminus type from first residue." << std::endl;
			} else {
				TR << "No lower terminus type found on first residue.  (The PDB file may contain a LINK record specifying cyclic geometry.)" << std::endl;
			}
		}
		if ( cyclization_type() == SCPA_n_to_c_amide_bond || cyclization_type() == SCPA_cterm_isopeptide_lariat ) {
			if ( native_pose->residue_type(last_res).is_upper_terminus() ) {
				core::pose::remove_upper_terminus_type_from_pose_residue(*native_pose, last_res);
				TR << "Removed upper terminus type from last residue." << std::endl;
			} else {
				TR << "No upper terminus type found on last residue.  (The PDB file may contain a LINK record specifying cyclic geometry.)" << std::endl;
			}
		}
	} else if ( cyclization_type() == SCPA_thioether_lariat ) {
		thioether_sidechain_index = set_up_terminal_thioether_lariat_variants( native_pose ); // Returns thioether cysteine index.
	} else if ( cyclization_type() == SCPA_terminal_disulfide ) {
		//Set up disulfide variants, if we're doing disulfide cyclization.
		set_up_terminal_disulfide_variants( native_pose );
	}

	//Mover to cyclize the polymer and to update terminal peptide bond O and H atoms:
	protocols::simple_moves::DeclareBondOP termini( utility::pointer::make_shared< protocols::simple_moves::DeclareBond >() );
	set_up_cyclization_mover( termini, native_pose, true, last_res );
	termini->apply(*native_pose);
	if ( cyclization_type() == SCPA_thioether_lariat ) {
		runtime_assert( thioether_sidechain_index != 0 ); //Should be true.
		protocols::cyclic_peptide::crosslinker::correct_thioether_virtuals( *native_pose, 1, thioether_sidechain_index );
	} else if ( cyclization_type() == SCPA_lanthipeptide ) {
		if ( lanthionine_positions_.size() > 0 ) {
			for ( core::Size i=1, imax=lanthionine_positions_.size(); i<=imax; ++i ) { //Loop through all sets of doubles of residues.
				debug_assert(lanthionine_positions_[i].size() == 2); //Should always be true.
				protocols::cyclic_peptide::crosslinker::correct_lanthionine_virtuals(*native_pose, lanthionine_positions_[i][1], lanthionine_positions_[i][2]);
			}
		} else {
			runtime_assert_string_msg( lanthionine_positions_.size() > 0, "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication lanthionine positions are empty");
		}
	}
}

/// @brief Add cutpoint variants to the terminal residues of an N-to-C cyclic peptide.
///
void
SimpleCycpepPredictApplication::add_cutpoint_variants_at_termini(
	core::pose::PoseOP pose
) const {
	TR << "Setting up cutpoint variants." << std::endl;

	core::pose::correctly_remove_variants_incompatible_with_lower_cutpoint_variant( *pose, pose->total_residue() );
	core::pose::correctly_remove_variants_incompatible_with_upper_cutpoint_variant( *pose, 1 );
	core::pose::correctly_add_cutpoint_variants( *pose, pose->total_residue(), false, 1 );
}

/// @brief Function to add cyclic constraints to a pose.
///
void
SimpleCycpepPredictApplication::add_amide_bond_cyclic_constraints (
	core::pose::PoseOP pose,
	core::Size n_index,
	core::Size c_index
) const {
	using namespace core::pose;
	using namespace core::scoring::constraints;
	using namespace core::scoring::func;
	using namespace core::id;

	TR << "Setting up cyclic constraints." << std::endl;

	bool n_term(false), c_term(false);

	if ( !n_index ) {
		n_index = 1;
		n_term = true;
	}
	if ( !c_index ) {
		c_index = sequence_length();
		c_term = true;
	}

	//The residue types of the N and C residues:
	core::chemical::ResidueType const & ntype( pose->residue_type(n_index) );
	core::chemical::ResidueType const & ctype( pose->residue_type(c_index) );

	//The four atoms defining the peptide bond:
	AtomID const atom_a(
		c_term ?
		ctype.icoor(ctype.upper_connect_atom()).stub_atom1().atomno() :
		ctype.icoor(ctype.residue_connect_atom_index(ctype.n_possible_residue_connections())).stub_atom1().atomno()
		, c_index );
	AtomID const atom_b(
		c_term ?
		ctype.upper_connect_atom() :
		ctype.residue_connect_atom_index( ctype.n_possible_residue_connections() )
		, c_index );
	AtomID const atom_c(
		n_term ?
		ntype.lower_connect_atom() :
		ntype.residue_connect_atom_index( ntype.n_possible_residue_connections() )
		, n_index );
	core::Size atom_d_index(0);
	if ( n_term ) {
		for ( core::Size i=1, imax=ntype.mainchain_atoms().size(); i<=imax; ++i ) { //Find the atom index of the first mainchain atom with the lower_connect atom as a parent.
			if ( i == atom_c.atomno() ) continue;
			if ( ntype.icoor(i).stub_atom1().atomno() == atom_c.atomno() ) {
				atom_d_index=i;
				break;
			}
		}
	} else {
		atom_d_index = ntype.icoor( ntype.residue_connect_atom_index( ntype.n_possible_residue_connections() ) ).stub_atom1().atomno();
	}
	AtomID const atom_d( atom_d_index, n_index );

	TR << "The following four atoms define the " << ( c_term && n_term ? "terminal" : "isopeptide" ) << " bond:" << std::endl;
	TR << "1.\tRes=" << atom_a.rsd() << "\tAtom=" << pose->residue(atom_a.rsd()).atom_name(atom_a.atomno()) << std::endl;
	TR << "2.\tRes=" << atom_b.rsd() << "\tAtom=" << pose->residue(atom_b.rsd()).atom_name(atom_b.atomno()) << std::endl;
	TR << "3.\tRes=" << atom_c.rsd() << "\tAtom=" << pose->residue(atom_c.rsd()).atom_name(atom_c.atomno()) << std::endl;
	TR << "4.\tRes=" << atom_d.rsd() << "\tAtom=" << pose->residue(atom_d.rsd()).atom_name(atom_d.atomno()) << std::endl;

	{//Peptide bond length constraint:
		FuncOP harmfunc1( utility::pointer::make_shared< HarmonicFunc >( SimpleCycpepPredictApplication_PEPBOND_LENGTH, 0.01) );
		ConstraintCOP distconst1( utility::pointer::make_shared< AtomPairConstraint >( atom_b, atom_c, harmfunc1 ) );
		pose->add_constraint (distconst1);
	}

	{ //Peptide dihedral angle constraints:
		// (TODO -- change these if we sample a trans-proline.)
		FuncOP circharmfunc1( utility::pointer::make_shared< CircularHarmonicFunc >( numeric::constants::d::pi, 0.02) );
		ConstraintCOP dihedconst1( utility::pointer::make_shared< DihedralConstraint >( atom_a, atom_b, atom_c, atom_d, circharmfunc1) );
		pose->add_constraint (dihedconst1);
	}

	{ //Peptide bond angle constraints:
		FuncOP circharmfunc2a( utility::pointer::make_shared< CircularHarmonicFunc >( SimpleCycpepPredictApplication_PEPBOND_C_ANGLE, 0.02) );
		FuncOP circharmfunc2b( utility::pointer::make_shared< CircularHarmonicFunc >( SimpleCycpepPredictApplication_PEPBOND_N_ANGLE, 0.02) );
		ConstraintCOP angleconst1( utility::pointer::make_shared< AngleConstraint >( atom_a, atom_b, atom_c, circharmfunc2a) );
		ConstraintCOP angleconst2( utility::pointer::make_shared< AngleConstraint >( atom_b, atom_c, atom_d, circharmfunc2b) );
		pose->add_constraint (angleconst1);
		pose->add_constraint (angleconst2);
	}

	TR << "Finished setting up constraints." << std::endl;

	return;
}

/// @brief Function to add thioether lariat cyclic constraints to a pose.
///
void
SimpleCycpepPredictApplication::add_thioether_lariat_cyclic_constraints (
	core::pose::PoseOP pose,
	core::Size n_index,
	core::Size c_index
) const {
	protocols::cyclic_peptide::crosslinker::set_up_thioether_constraints(
		*pose,
		n_index,
		c_index
	);
}

/// @brief Function to add lanthipeptide cyclic constraints to a pose.
///
void
SimpleCycpepPredictApplication::add_lanthipeptide_cyclic_constraints (
	core::pose::PoseOP pose
	//core::Size n_index,
	//core::Size c_index
	//first res should be ala equiv, second cys
) const {
	//protocols::cyclic_peptide::crosslinker::set_up_lanthipeptide_constraints(
	// *pose,
	// n_index,
	// c_index
	//);
	if ( lanthionine_positions_.size() > 0 ) {
		for ( core::Size i=1, imax=lanthionine_positions_.size(); i<=imax; ++i ) { //Loop through all sets of doubles of residues.
			debug_assert(lanthionine_positions_[i].size() == 2); //Should always be true.
			protocols::cyclic_peptide::crosslinker::set_up_lanthionine_constraints(*pose, lanthionine_positions_[i][1], lanthionine_positions_[i][2]);
		}
	} else {
		runtime_assert_string_msg( lanthionine_positions_.size() > 0, "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication lanthionine positions are empty");
	}
}

/// @brief Function to add cyclic constraints to a pose.
/// @details Calls functions that do this for particular cyclization types.
void
SimpleCycpepPredictApplication::add_cyclic_constraints(
	core::pose::PoseOP pose
) const {
	switch( cyclization_type() ) {
	case SCPA_n_to_c_amide_bond :
		add_amide_bond_cyclic_constraints( pose, 0, 0 );
		return;
	case SCPA_terminal_disulfide :
		//Do nothing.  The constraints are handled by the dslf_fa13 score term.
		return;
	case SCPA_thioether_lariat :
		core::Size firstres, lastres;
		find_first_and_last_thioether_lariat_residues( pose, firstres, lastres );
		add_thioether_lariat_cyclic_constraints( pose, firstres, lastres );
		return;
	case SCPA_lanthipeptide :
		add_lanthipeptide_cyclic_constraints(pose);
		return;
	case SCPA_nterm_isopeptide_lariat :
		{
		core::Size firstres, lastres;
		find_first_and_last_isopeptide_residues( pose, firstres, lastres );
		add_amide_bond_cyclic_constraints( pose, 0, lastres );
		return;
	}
	case SCPA_cterm_isopeptide_lariat :
		{
		core::Size firstres, lastres;
		find_first_and_last_isopeptide_residues( pose, firstres, lastres );
		add_amide_bond_cyclic_constraints( pose, firstres, 0 );
		return;
	}
	case SCPA_sidechain_isopeptide :
		{
		core::Size firstres, lastres;
		find_first_and_last_isopeptide_residues( pose, firstres, lastres );
		add_amide_bond_cyclic_constraints( pose, firstres, lastres );
		return;
	}
	case SCPA_invalid_type :
		utility_exit_with_message("Error in SimpleCycpepPredictApplication::add_cyclic_constraints(): Invalid cyclization type!");
	}
}

/// @brief Sets all omega values to 180, and randomizes mainchain torsions.
/// @details For alpha-amino acids, mainchain torsions are randomized by the Ramachandran plot.
/// For other residue types, just randomizes mainchain torsions other than peptide bonds.
void
SimpleCycpepPredictApplication::set_mainchain_torsions (
	core::pose::PoseOP pose,
	core::Size const cyclic_offset
) const {
	TR << "Randomizing mainchain torsions." << std::endl;
	core::Size const nres(sequence_length());
	core::scoring::Ramachandran const & rama(core::scoring::ScoringManager::get_instance()->get_Ramachandran() ); //Get the Rama scoring function
	core::scoring::RamaPrePro const & ramaprepro( core::scoring::ScoringManager::get_instance()->get_RamaPrePro() );

	for ( core::Size i=1; i<=nres; ++i ) { //Loop through all residues
		if ( is_supported_restype( pose->residue_type(i) ) ) {
			utility::vector1< core::Real > rand_torsions( pose->residue(i).mainchain_torsions().size() - 1, 0.0 );
			core::Size const cur_abs_pos( original_position( i, cyclic_offset, pose->size() ) );
			if ( custom_rama_table_defined( cur_abs_pos ) ) {
				if ( pose->residue_type(i).is_alpha_aa() ) {
					rama.draw_random_phi_psi_from_extra_cdf( rama_table_type_by_res(cur_abs_pos), rand_torsions[1], rand_torsions[2]);
				}
			} else if ( default_rama_table_type() != core::scoring::unknown_ramatable_type ) {
				if ( pose->residue_type(i).is_alpha_aa() ) {
					rama.draw_random_phi_psi_from_extra_cdf( default_rama_table_type(), rand_torsions[1], rand_torsions[2]);
				}
			} else {
				if ( use_rama_prepro_for_sampling() ) { //Using rama_prepro tables for sampling:
					// Use a fake alanine for following_rsd for the cterm residue unless you're doing
					// cterm cyclization and thus you can use its upper-connected residue.
					core::chemical::ResidueTypeCOP following_rsd(
						( cyclization_type() == SCPA_terminal_disulfide || cyclization_type() == SCPA_nterm_isopeptide_lariat || cyclization_type() == SCPA_sidechain_isopeptide || cyclization_type() == SCPA_thioether_lariat || cyclization_type() == SCPA_lanthipeptide ) && i == nres
						?
						core::chemical::ResidueTypeFinder( *( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) ) ).residue_base_name("ALA").get_representative_type()
						:
						pose->residue_type_ptr(
						pose->residue(i).residue_connection_partner( pose->residue(i).upper_connect().index() )
						)
					);
					ramaprepro.random_mainchain_torsions( pose->conformation(), pose->residue_type_ptr(i), following_rsd, rand_torsions );
					utility::vector1< core::Size > const & covered_torsions( ramaprepro.get_mainchain_torsions_covered( pose->conformation(), pose->residue_type_ptr(i), following_rsd ) );
					debug_assert( pose->residue(i).mainchain_torsions().size() - 1 == rand_torsions.size() ); //Should be true.
					for ( core::Size j(1), jmax(rand_torsions.size()); j<=jmax; ++j ) {
						if ( !covered_torsions.has_value(j) ) {
							rand_torsions[j] = pose->residue(i).mainchain_torsions()[j]; //Copy torsion values for torsions that are not covered by the mainchain potential -- i.e. don't set these to 0.
						}
					}
				} else { //Using classic rama tables for sampling:
					if ( pose->residue(i).backbone_aa() != core::chemical::aa_unk ) {
						rama.random_phipsi_from_rama(pose->residue(i).backbone_aa(), rand_torsions[1], rand_torsions[2]);
					} else {
						runtime_assert_string_msg( pose->residue_type(i).aa() != core::chemical::aa_unk, "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::set_mainchain_torsions(): Unable to get a suitable classic Ramachandran table for " + pose->residue_type(i).name() + "." );
						rama.random_phipsi_from_rama( pose->residue_type(i).aa(), rand_torsions[1], rand_torsions[2]);
					}
				}

			}
			debug_assert(rand_torsions.size() == pose->residue(i).mainchain_torsions().size() - 1);
			for ( core::Size j(1), jmax(rand_torsions.size()); j<=jmax; ++j ) {
				pose->set_torsion( core::id::TorsionID( i, core::id::BB, j ), rand_torsions[j] );
			}
			if ( i!=nres ) pose->set_omega(i, 180.0);
			if ( pose->residue_type(i).is_oligourea() ) pose->set_mu(i, 180.0);
		} else { //If this is not a recognized type:
			for ( core::Size j=1, jmax=pose->residue(i).mainchain_torsions().size(); j<=jmax; ++j ) { //Loop through all mainchain torsions.
				if ( i==nres && j==jmax ) continue; //Skip the last mainchain torsion (not a DOF).
				core::Real setting(180.0);
				if ( j!=jmax ) {
					setting = numeric::random::rg().uniform()*360.0 - 180.0;
				}
				pose->set_torsion( core::id::TorsionID(i, core::id::BB, j), setting );
			}
		}
	}
	return;
}

/// @brief Set up the filters for the mainchain hydrogen bonds that will
/// be used to discard solutions with too few mainchain hydrogen bonds.
protocols::filters::FilterOP
SimpleCycpepPredictApplication::set_up_hbond_filter(
	core::Size const min_hbonds
) const {
	protocols::cyclic_peptide::PeptideInternalHbondsFilterOP total_hbond(
		utility::pointer::make_shared<protocols::cyclic_peptide::PeptideInternalHbondsFilter>()
	);
	total_hbond->set_exclusion_distance( do_not_count_adjacent_res_hbonds_ ? 1 : 0 );
	total_hbond->set_hbond_cutoff( min_hbonds );
	total_hbond->set_hbond_types(true, count_sc_hbonds_, false);
	total_hbond->set_hbond_energy_cutoff(hbond_energy_cutoff_);

	return total_hbond;
}

/// @brief Set up the logic to close the bond at the cyclization point.
/// @details This version is for N-to-C amide bond cyclization.
void
SimpleCycpepPredictApplication::add_closebond_logic_n_to_c_amide_bond(
	core::pose::PoseCOP pose,
	core::Size const /*cyclization_point_start*/,
	core::Size const /*cyclization_point_end*/,
	protocols::generalized_kinematic_closure::GeneralizedKICOP genkic
) const {
	debug_assert ( cyclization_type() == SCPA_n_to_c_amide_bond );
	std::string const firstatom( pose->residue(1).atom_name( pose->residue(1).lower_connect_atom() ) );
	std::string const lastatom( pose->residue(sequence_length()).atom_name( pose->residue(sequence_length()).upper_connect_atom() ) );
	genkic->close_bond( sequence_length(), lastatom, 1, firstatom, 0, "", 0, "", SimpleCycpepPredictApplication_PEPBOND_LENGTH, SimpleCycpepPredictApplication_PEPBOND_C_ANGLE/numeric::constants::d::pi*180.0, SimpleCycpepPredictApplication_PEPBOND_N_ANGLE/numeric::constants::d::pi*180.0, 180.0, false, false );
}

/// @brief Set up the logic to close the bond at the cyclization point.
/// @details This version is for terminal disulfide cyclization.
void
SimpleCycpepPredictApplication::add_closebond_logic_terminal_disulfide(
	core::pose::PoseCOP pose,
	core::Size const cyclization_point_start,
	core::Size const cyclization_point_end,
	protocols::generalized_kinematic_closure::GeneralizedKICOP genkic
) const {
	debug_assert( cyclization_type() == SCPA_terminal_disulfide );
	//Close S-S bond:
	std::string const startatom( pose->residue_type(cyclization_point_start).get_disulfide_atom_name() );
	std::string const endatom( pose->residue_type(cyclization_point_end).get_disulfide_atom_name() );
	genkic->close_bond( cyclization_point_start, startatom, cyclization_point_end, endatom, 0, "", 0, "",
		SimpleCycpepPredictApplication_DISULFBOND_LENGTH,
		SimpleCycpepPredictApplication_DISULFBOND_ANGLE/numeric::constants::d::pi*180.0,
		SimpleCycpepPredictApplication_DISULFBOND_ANGLE/numeric::constants::d::pi*180.0,
		0 /*Will be sampled later*/,
		false,
		false
	);

	//Randomize S-S bond:
	{
		genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::randomize_dihedral );
		utility::vector1< core::id::NamedAtomID > SSvect;
		SSvect.push_back( core::id::NamedAtomID( startatom, cyclization_point_start ) );
		SSvect.push_back( core::id::NamedAtomID( endatom, cyclization_point_end ) );
		genkic->add_atomset_to_perturber_atomset_list( SSvect );
	}

	//Randomize sidechains:
	for ( core::Size l1(1); l1<=2; ++l1 ) {
		core::Size const curres( l1 == 1 ? cyclization_point_start : cyclization_point_end );
		genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::randomize_dihedral );
		core::Size const first_sc_at( pose->residue_type(curres).first_sidechain_atom() );
		core::Size const first_sc_at_parent( pose->residue_type(curres).icoor(first_sc_at).stub_atom1().atomno() );
		for ( core::Size curat( pose->residue_type(curres).atom_index( pose->residue_type(curres).get_disulfide_atom_name() ) ); curat!=first_sc_at_parent; curat = pose->residue_type(curres).icoor(curat).stub_atom1().atomno() ) { //March up from the thiol atom
			core::Size const otherat( pose->residue_type(curres).icoor(curat).stub_atom1().atomno() ); //Get the parent of the current atom
			utility::vector1< core::id::NamedAtomID > sc_vect;
			sc_vect.push_back( core::id::NamedAtomID( pose->residue_type(curres).atom_name(curat), curres ) );
			sc_vect.push_back( core::id::NamedAtomID( pose->residue_type(curres).atom_name(otherat), curres ) );
			genkic->add_atomset_to_perturber_atomset_list( sc_vect );
		}
	}

	//Randomize upper cysteine backbone:
	{
		genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::randomize_dihedral );
		utility::vector1< core::id::NamedAtomID > bb_vect;
		if ( pose->residue_type( cyclization_point_start ).is_alpha_aa() || pose->residue_type(cyclization_point_start).is_oligourea() || pose->residue_type( cyclization_point_start ).is_peptoid() || pose->residue_type(cyclization_point_start).is_beta_aa() ) {
			bb_vect.push_back( core::id::NamedAtomID( "N", cyclization_point_start ) );
			bb_vect.push_back( core::id::NamedAtomID( "CA", cyclization_point_start ) );
		} else if ( pose->residue_type(cyclization_point_start).is_gamma_aa() ) {
			bb_vect.push_back( core::id::NamedAtomID( "N", cyclization_point_start ) );
			bb_vect.push_back( core::id::NamedAtomID( "C4", cyclization_point_start ) );
		}
		genkic->add_atomset_to_perturber_atomset_list( bb_vect );
	}

	//Randomize lower cysteine backbone:
	{
		genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::randomize_dihedral );
		if ( pose->residue_type( cyclization_point_end ).is_alpha_aa() || pose->residue_type( cyclization_point_end ).is_peptoid() ) {
			utility::vector1< core::id::NamedAtomID > bb_vect;
			bb_vect.push_back( core::id::NamedAtomID( "CA", cyclization_point_end ) );
			bb_vect.push_back( core::id::NamedAtomID( "C", cyclization_point_end ) );
			genkic->add_atomset_to_perturber_atomset_list( bb_vect );
		} else if ( pose->residue_type(cyclization_point_end).is_beta_aa() || pose->residue_type(cyclization_point_end).is_oligourea() ) {
			utility::vector1< core::id::NamedAtomID > bb_vect1, bb_vect2;
			bb_vect1.push_back( core::id::NamedAtomID( "CA", cyclization_point_end ) );
			bb_vect1.push_back( core::id::NamedAtomID( "CM", cyclization_point_end ) );
			genkic->add_atomset_to_perturber_atomset_list( bb_vect1 );
			bb_vect2.push_back( core::id::NamedAtomID( "CM", cyclization_point_end ) );
			bb_vect2.push_back( core::id::NamedAtomID( "C", cyclization_point_end ) );
			genkic->add_atomset_to_perturber_atomset_list( bb_vect2 );
		} else if ( pose->residue_type(cyclization_point_end).is_gamma_aa() ) {
			utility::vector1< core::id::NamedAtomID > bb_vect1, bb_vect2, bb_vect3;
			bb_vect1.push_back( core::id::NamedAtomID( "C4", cyclization_point_end ) );
			bb_vect1.push_back( core::id::NamedAtomID( "C3", cyclization_point_end ) );
			genkic->add_atomset_to_perturber_atomset_list( bb_vect1 );
			bb_vect2.push_back( core::id::NamedAtomID( "C3", cyclization_point_end ) );
			bb_vect2.push_back( core::id::NamedAtomID( "C2", cyclization_point_end ) );
			genkic->add_atomset_to_perturber_atomset_list( bb_vect2 );
			bb_vect3.push_back( core::id::NamedAtomID( "C2", cyclization_point_end ) );
			bb_vect3.push_back( core::id::NamedAtomID( "C", cyclization_point_end ) );
			genkic->add_atomset_to_perturber_atomset_list( bb_vect3 );
		}
	}
}

/// @brief Set up the logic to close the bond at the cyclization point.
/// @details This version is for all types of isopeptide bond cyclization.
/// @note cyclization_point_start > cyclization_point_end
void
SimpleCycpepPredictApplication::add_closebond_logic_isopeptide(
	core::pose::PoseCOP pose,
	core::Size const cyclization_point_start,
	core::Size const cyclization_point_end,
	protocols::generalized_kinematic_closure::GeneralizedKICOP genkic
) const {
	debug_assert( cyclization_type() == SCPA_nterm_isopeptide_lariat || cyclization_type() == SCPA_cterm_isopeptide_lariat || cyclization_type() == SCPA_sidechain_isopeptide );
	debug_assert( cyclization_point_start > cyclization_point_end );

	//Some stuff we'll use:
	core::chemical::ResidueType const & rt_start( pose->residue_type( cyclization_point_start ) );
	core::chemical::ResidueType const & rt_end( pose->residue_type( cyclization_point_end ) );
	bool const backward_bond( cyclization_type() == SCPA_sidechain_isopeptide && is_isopeptide_forming_amide_type( rt_end.name3() ) && is_isopeptide_forming_carbonyl_type( rt_start.aa() ));

	//Close isopeptide bond:
	std::string const startatom (
		cyclization_type() == SCPA_cterm_isopeptide_lariat ?
		rt_start.atom_name( rt_start.upper_connect_atom() ) :
		rt_start.atom_name( rt_start.residue_connect_atom_index( rt_start.n_possible_residue_connections() ) )
	);
	std::string const endatom (
		cyclization_type() == SCPA_nterm_isopeptide_lariat ?
		rt_end.atom_name( rt_end.lower_connect_atom() ) :
		rt_end.atom_name( rt_end.residue_connect_atom_index( rt_end.n_possible_residue_connections() ) )
	);
	if ( backward_bond ) {
		genkic->close_bond( cyclization_point_start, startatom, cyclization_point_end, endatom, 0, "", 0, "", SimpleCycpepPredictApplication_PEPBOND_LENGTH, SimpleCycpepPredictApplication_PEPBOND_C_ANGLE/numeric::constants::d::pi*180.0, SimpleCycpepPredictApplication_PEPBOND_N_ANGLE/numeric::constants::d::pi*180.0, 180.0, false, false );
	} else {
		genkic->close_bond( cyclization_point_start, startatom, cyclization_point_end, endatom, 0, "", 0, "", SimpleCycpepPredictApplication_PEPBOND_LENGTH, SimpleCycpepPredictApplication_PEPBOND_N_ANGLE/numeric::constants::d::pi*180.0, SimpleCycpepPredictApplication_PEPBOND_C_ANGLE/numeric::constants::d::pi*180.0, 180.0, false, false );
	}

	//Randomize sidechains:
	for ( core::Size l1(1); l1<=2; ++l1 ) {
		core::Size const curres( l1 == 1 ? cyclization_point_start : cyclization_point_end );
		if ( curres == cyclization_point_start && cyclization_type() == SCPA_cterm_isopeptide_lariat ) continue;
		if ( curres == cyclization_point_end && cyclization_type() == SCPA_nterm_isopeptide_lariat ) continue;

		genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::randomize_dihedral );
		core::Size const first_sc_at( pose->residue_type(curres).first_sidechain_atom() );
		core::chemical::ResidueType const & curtype( pose->residue_type(curres) );
		for ( core::Size curat( curtype.residue_connect_atom_index( curtype.n_possible_residue_connections() ) ); curat!=first_sc_at; curat = curtype.icoor(curat).stub_atom1().atomno() ) { //March up from the thiol atom
			core::Size const otherat( curtype.icoor(curat).stub_atom1().atomno() ); //Get the parent of the current atom
			utility::vector1< core::id::NamedAtomID > sc_vect;
			sc_vect.push_back( core::id::NamedAtomID( curtype.atom_name(curat), curres ) );
			sc_vect.push_back( core::id::NamedAtomID( curtype.atom_name(otherat), curres ) );
			genkic->add_atomset_to_perturber_atomset_list( sc_vect );
		}
	}

	//Randomize backbone of upper res:
	if ( cyclization_type() != SCPA_cterm_isopeptide_lariat ) {
		genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::randomize_dihedral );
		utility::vector1< core::id::NamedAtomID > bb_vect;
		if ( pose->residue_type( cyclization_point_start ).is_alpha_aa() || pose->residue_type(cyclization_point_start).is_oligourea() || pose->residue_type( cyclization_point_start ).is_peptoid() || pose->residue_type(cyclization_point_start).is_beta_aa() ) {
			bb_vect.push_back( core::id::NamedAtomID( "N", cyclization_point_start ) );
			bb_vect.push_back( core::id::NamedAtomID( "CA", cyclization_point_start ) );
		} else if ( pose->residue_type(cyclization_point_start).is_gamma_aa() ) {
			bb_vect.push_back( core::id::NamedAtomID( "N", cyclization_point_start ) );
			bb_vect.push_back( core::id::NamedAtomID( "C4", cyclization_point_start ) );
		}
		genkic->add_atomset_to_perturber_atomset_list( bb_vect );
	}

	//Randomize lower res backbone:
	if ( cyclization_type() != SCPA_nterm_isopeptide_lariat ) {
		genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::randomize_dihedral );
		if ( pose->residue_type( cyclization_point_end ).is_alpha_aa() || pose->residue_type( cyclization_point_end ).is_peptoid() ) {
			utility::vector1< core::id::NamedAtomID > bb_vect;
			bb_vect.push_back( core::id::NamedAtomID( "CA", cyclization_point_end ) );
			bb_vect.push_back( core::id::NamedAtomID( "C", cyclization_point_end ) );
			genkic->add_atomset_to_perturber_atomset_list( bb_vect );
		} else if ( pose->residue_type(cyclization_point_end).is_beta_aa() || pose->residue_type(cyclization_point_end).is_oligourea() ) {
			utility::vector1< core::id::NamedAtomID > bb_vect1, bb_vect2;
			bb_vect1.push_back( core::id::NamedAtomID( "CA", cyclization_point_end ) );
			bb_vect1.push_back( core::id::NamedAtomID( "CM", cyclization_point_end ) );
			genkic->add_atomset_to_perturber_atomset_list( bb_vect1 );
			bb_vect2.push_back( core::id::NamedAtomID( "CM", cyclization_point_end ) );
			bb_vect2.push_back( core::id::NamedAtomID( "C", cyclization_point_end ) );
			genkic->add_atomset_to_perturber_atomset_list( bb_vect2 );
		} else if ( pose->residue_type(cyclization_point_end).is_gamma_aa() ) {
			utility::vector1< core::id::NamedAtomID > bb_vect1, bb_vect2, bb_vect3;
			bb_vect1.push_back( core::id::NamedAtomID( "C4", cyclization_point_end ) );
			bb_vect1.push_back( core::id::NamedAtomID( "C3", cyclization_point_end ) );
			genkic->add_atomset_to_perturber_atomset_list( bb_vect1 );
			bb_vect2.push_back( core::id::NamedAtomID( "C3", cyclization_point_end ) );
			bb_vect2.push_back( core::id::NamedAtomID( "C2", cyclization_point_end ) );
			genkic->add_atomset_to_perturber_atomset_list( bb_vect2 );
			bb_vect3.push_back( core::id::NamedAtomID( "C2", cyclization_point_end ) );
			bb_vect3.push_back( core::id::NamedAtomID( "C", cyclization_point_end ) );
			genkic->add_atomset_to_perturber_atomset_list( bb_vect3 );
		}
	}

}

/// @brief Set up the logic to close the bond at the cyclization point.
/// @details This version is for thioether lariat cyclization.
/// @note cyclization_point_start > cyclization_point_end
void
SimpleCycpepPredictApplication::add_closebond_logic_thioether_lariat(
	core::pose::PoseCOP pose,
	core::Size const cyclization_point_start,
	core::Size const cyclization_point_end,
	protocols::generalized_kinematic_closure::GeneralizedKICOP genkic
) const {
	debug_assert( cyclization_type() == SCPA_thioether_lariat );
	debug_assert( cyclization_point_start > cyclization_point_end );

	// This means start is cys and end is nterm?

	//Some stuff we'll use:
	core::chemical::ResidueType const & rt_start( pose->residue_type( cyclization_point_start ) );
	// core::chemical::ResidueType const & rt_end( pose->residue_type( cyclization_point_end ) );

	//Close thioether bond:
	std::string const startatom ( rt_start.get_disulfide_atom_name() );
	std::string const endatom ( "CP2" );
	genkic->close_bond(
		cyclization_point_start,
		startatom,
		cyclization_point_end,
		endatom,
		0, "", 0, "",
		protocols::cyclic_peptide::crosslinker::THIOETHER_UTIL_THIOETHER_BOND_LENGTH,
		protocols::cyclic_peptide::crosslinker::THIOETHER_UTIL_THIOETHER_BOND_N_ANGLE/numeric::constants::d::pi*180.0,
		protocols::cyclic_peptide::crosslinker::THIOETHER_UTIL_THIOETHER_BOND_C_ANGLE/numeric::constants::d::pi*180.0,
		180.0, false, false );

	//Randomize sidechains:
	// core::Size const first_sc_at( pose->residue_type(cyclization_point_end).first_sidechain_atom() );
	// core::chemical::ResidueType const & curtype( pose->residue_type(cyclization_point_end) );
	// for ( core::Size curat( curtype.residue_connect_atom_index( curtype.n_possible_residue_connections() ) ); curat!=first_sc_at; curat = curtype.icoor(curat).stub_atom1().atomno() ) { //March up from the thiol atom
	//  core::Size const otherat( curtype.icoor(curat).stub_atom1().atomno() ); //Get the parent of the current atom
	//  utility::vector1< core::id::NamedAtomID > sc_vect;
	//  sc_vect.push_back( core::id::NamedAtomID( curtype.atom_name(curat), cyclization_point_end ) );
	//  sc_vect.push_back( core::id::NamedAtomID( curtype.atom_name(otherat), cyclization_point_end ) );
	//  genkic->add_atomset_to_perturber_atomset_list( sc_vect );
	// }
	// I think the above is strictly unnecessary; no one needs to randomize the un-involved sidechain
	// and I'm going to say that we basically cannot randomize any dihedrals in 2-chloro-acetyl.
	// In theory you could randomize CONN2-CP2-CO-N and CO-N-CA-C but that's tomorrow's problem. AMW TODO

	//This application appears to only produce cis CA-N-CO-CP2 bonds without randomization,
	// trans should be energetically favored, so setting bond to trans here if option specified
	//lariat_sample_cis_ should be false for this option to sample trans peptide bonds
	if ( !lariat_sample_cis_ ) {
		genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::set_dihedral );
		utility::vector1< core::id::NamedAtomID > omega_vect;
		omega_vect.push_back( core::id::NamedAtomID( "N", cyclization_point_end));
		omega_vect.push_back( core::id::NamedAtomID( "CO", cyclization_point_end));
		genkic->add_atomset_to_perturber_atomset_list( omega_vect );
		genkic->add_value_to_perturber_value_list( 180.0 );
	}

	genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::randomize_dihedral );
	core::Size const first_sc_at( pose->residue_type(cyclization_point_start).first_sidechain_atom() );
	core::chemical::ResidueType const & curtype( pose->residue_type(cyclization_point_start) );
	for ( core::Size curat( curtype.residue_connect_atom_index( curtype.n_possible_residue_connections() ) ); curat!=first_sc_at; curat = curtype.icoor(curat).stub_atom1().atomno() ) { //March up from the thiol atom
		core::Size const otherat( curtype.icoor(curat).stub_atom1().atomno() ); //Get the parent of the current atom
		utility::vector1< core::id::NamedAtomID > sc_vect;
		sc_vect.push_back( core::id::NamedAtomID( curtype.atom_name(curat), cyclization_point_start ) );
		sc_vect.push_back( core::id::NamedAtomID( curtype.atom_name(otherat), cyclization_point_start ) );
		genkic->add_atomset_to_perturber_atomset_list( sc_vect );
	}

	//Randomize backbone of upper res:
	genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::randomize_dihedral );
	utility::vector1< core::id::NamedAtomID > bb_vect;
	if ( pose->residue_type( cyclization_point_start ).is_alpha_aa()
			|| pose->residue_type(cyclization_point_start).is_oligourea()
			|| pose->residue_type( cyclization_point_start ).is_peptoid()
			|| pose->residue_type(cyclization_point_start).is_beta_aa() ) {
		bb_vect.push_back( core::id::NamedAtomID( "N", cyclization_point_start ) );
		bb_vect.push_back( core::id::NamedAtomID( "CA", cyclization_point_start ) );
	} else if ( pose->residue_type(cyclization_point_start).is_gamma_aa() ) {
		bb_vect.push_back( core::id::NamedAtomID( "N", cyclization_point_start ) );
		bb_vect.push_back( core::id::NamedAtomID( "C4", cyclization_point_start ) );
	}
	genkic->add_atomset_to_perturber_atomset_list( bb_vect );

	//Randomize lower res backbone:
	genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::randomize_dihedral );
	if ( pose->residue_type( cyclization_point_end ).is_alpha_aa() || pose->residue_type( cyclization_point_end ).is_peptoid() ) {
		utility::vector1< core::id::NamedAtomID > bb_vect;
		bb_vect.push_back( core::id::NamedAtomID( "CA", cyclization_point_end ) );
		bb_vect.push_back( core::id::NamedAtomID( "C", cyclization_point_end ) );
		genkic->add_atomset_to_perturber_atomset_list( bb_vect );
	} else if ( pose->residue_type(cyclization_point_end).is_beta_aa() || pose->residue_type(cyclization_point_end).is_oligourea() ) {
		utility::vector1< core::id::NamedAtomID > bb_vect1, bb_vect2;
		bb_vect1.push_back( core::id::NamedAtomID( "CA", cyclization_point_end ) );
		bb_vect1.push_back( core::id::NamedAtomID( "CM", cyclization_point_end ) );
		genkic->add_atomset_to_perturber_atomset_list( bb_vect1 );
		bb_vect2.push_back( core::id::NamedAtomID( "CM", cyclization_point_end ) );
		bb_vect2.push_back( core::id::NamedAtomID( "C", cyclization_point_end ) );
		genkic->add_atomset_to_perturber_atomset_list( bb_vect2 );
	} else if ( pose->residue_type(cyclization_point_end).is_gamma_aa() ) {
		utility::vector1< core::id::NamedAtomID > bb_vect1, bb_vect2, bb_vect3;
		bb_vect1.push_back( core::id::NamedAtomID( "C4", cyclization_point_end ) );
		bb_vect1.push_back( core::id::NamedAtomID( "C3", cyclization_point_end ) );
		genkic->add_atomset_to_perturber_atomset_list( bb_vect1 );
		bb_vect2.push_back( core::id::NamedAtomID( "C3", cyclization_point_end ) );
		bb_vect2.push_back( core::id::NamedAtomID( "C2", cyclization_point_end ) );
		genkic->add_atomset_to_perturber_atomset_list( bb_vect2 );
		bb_vect3.push_back( core::id::NamedAtomID( "C2", cyclization_point_end ) );
		bb_vect3.push_back( core::id::NamedAtomID( "C", cyclization_point_end ) );
		genkic->add_atomset_to_perturber_atomset_list( bb_vect3 );
	}
}

//Lanthionine for now
void
SimpleCycpepPredictApplication::add_closebond_logic_lanthipeptide(
	core::pose::PoseCOP pose,
	core::Size const cyclization_point_start,
	core::Size const cyclization_point_end,
	protocols::generalized_kinematic_closure::GeneralizedKICOP genkic
) const {
	debug_assert( cyclization_type() == SCPA_lanthipeptide );
	core::Size ala_res = 0;
	core::Size cys_res = 0;
	if ( pose->residue_type(cyclization_point_start).base_name() == "CYS" ||  pose->residue_type(cyclization_point_start).base_name() == "DCYS" ) {
		ala_res = cyclization_point_end;
		cys_res = cyclization_point_start;
	} else {
		ala_res = cyclization_point_start;
		cys_res = cyclization_point_end;
	}
	//} else if( pose->residue_type(cyclization_point_end).base_name() == "DALA" || pose->residue_type(cyclization_point_end).base_name() == "ALA") {
	// ala_res = cyclization_point_end;
	// cys_res = cyclization_point_start;
	//} else {
	// runtime_assert_string_msg(
	//  pose->residue_type(cyclization_point_start).base_name() == "DALA" || pose->residue_type(cyclization_point_start).base_name() == "ALA" ||
	//  pose->residue_type(cyclization_point_end).base_name() == "DALA" || pose->residue_type(cyclization_point_end).base_name() == "ALA",
	//  "Neither " + std::to_string( cyclization_point_end ) + " nor " + std::to_string( cyclization_point_start ) +
	//  "is a D-ALA or ALA."
	// );
	//}
	//Close lanthionine bond
	core::chemical::ResidueType const & rt_start( pose->residue_type( cys_res ) );
	std::string const startatom  ( rt_start.get_disulfide_atom_name() );
	std::string const endatom ( "CB" );
	genkic->close_bond(
		cys_res,
		startatom,
		ala_res,
		endatom,
		0, "", 0, "",
		protocols::cyclic_peptide::crosslinker::LANTHIONINE_UTIL_LANTHIONINE_BOND_LENGTH,
		protocols::cyclic_peptide::crosslinker::LANTHIONINE_UTIL_LANTHIONINE_BOND_C_ANGLE/numeric::constants::d::pi*180.0,
		protocols::cyclic_peptide::crosslinker::LANTHIONINE_UTIL_LANTHIONINE_BOND_S_ANGLE/numeric::constants::d::pi*180.0,
		180.0, false, false );

	//sidechain randomization
	genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::randomize_dihedral );

	core::Size const first_sc_at( pose->residue_type(cys_res).first_sidechain_atom() ); //Value is 5
	core::chemical::ResidueType const & curtype( pose->residue_type(cys_res) );
	for ( core::Size curat( curtype.residue_connect_atom_index( curtype.n_possible_residue_connections() ) );
			curat!=first_sc_at; curat = curtype.icoor(curat).stub_atom1().atomno() ) { //March up from the thiol (n.possible_residue_connections()==SG) atom
		core::Size const otherat( curtype.icoor(curat).stub_atom1().atomno() ); //Get the parent of the current atom
		utility::vector1< core::id::NamedAtomID > sc_vect;
		sc_vect.push_back( core::id::NamedAtomID( curtype.atom_name(curat), cys_res ) );
		sc_vect.push_back( core::id::NamedAtomID( curtype.atom_name(otherat), cys_res ) );
		genkic->add_atomset_to_perturber_atomset_list( sc_vect );
	}

	//Randomize upper res backbone:
	genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::randomize_dihedral );
	utility::vector1< core::id::NamedAtomID > upper_bb_vect;
	upper_bb_vect.push_back( core::id::NamedAtomID( "N", cyclization_point_start ) );
	upper_bb_vect.push_back( core::id::NamedAtomID( "CA", cyclization_point_start ) );
	genkic->add_atomset_to_perturber_atomset_list( upper_bb_vect );

	//Randomize lower res backbone:
	genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::randomize_dihedral );
	utility::vector1< core::id::NamedAtomID > lower_bb_vect;
	lower_bb_vect.push_back( core::id::NamedAtomID( "CA", cyclization_point_end ) );
	lower_bb_vect.push_back( core::id::NamedAtomID( "C", cyclization_point_end ) );
	genkic->add_atomset_to_perturber_atomset_list( lower_bb_vect );
}

/// @brief Set up the logic to close the bond at the cyclization point.
/// @details Calls different functions for different cyclization types.
void
SimpleCycpepPredictApplication::add_closebond_logic(
	core::pose::PoseCOP pose,
	core::Size const cyclization_point_start,
	core::Size const cyclization_point_end,
	protocols::generalized_kinematic_closure::GeneralizedKICOP genkic
) const {
	switch(cyclization_type()) {
	case SCPA_n_to_c_amide_bond :
		add_closebond_logic_n_to_c_amide_bond( pose, cyclization_point_start, cyclization_point_end, genkic );
		return;
	case SCPA_terminal_disulfide :
		add_closebond_logic_terminal_disulfide( pose, cyclization_point_start, cyclization_point_end, genkic );
		return;
		//The following three cases all call the same function:
	case SCPA_nterm_isopeptide_lariat:
	case SCPA_cterm_isopeptide_lariat:
	case SCPA_sidechain_isopeptide :
		add_closebond_logic_isopeptide( pose, cyclization_point_start, cyclization_point_end, genkic );
		return;
	case SCPA_thioether_lariat :
		add_closebond_logic_thioether_lariat( pose, cyclization_point_start, cyclization_point_end, genkic );
		return;
	case SCPA_lanthipeptide :
		add_closebond_logic_lanthipeptide( pose, cyclization_point_start, cyclization_point_end, genkic);
		return;
	case SCPA_invalid_type :
		utility_exit_with_message("Error in SimpleCycpepPredictApplication::add_closebond_logic(): Invalid cyclization type!");
	}
}


/// @brief Use GeneralizedKIC to close the pose.
///
bool
SimpleCycpepPredictApplication::genkic_close(
	core::pose::PoseOP pose,
	core::scoring::ScoreFunctionOP sfxn_highhbond,
	core::scoring::ScoreFunctionOP sfxn_highhbond_cart,
	core::scoring::ScoreFunctionCOP sfxn_default,
	protocols::filters::FilterOP total_hbond,
	core::Size const cyclic_offset
) const {
	using namespace protocols::generalized_kinematic_closure;
	using namespace protocols::generalized_kinematic_closure::selector;

	TR << "Performing GeneralizedKIC closure of loop." << std::endl;

	//Number of residues in the pose:
	core::Size const nres( sequence_length() );
	runtime_assert( nres >= 4 ); //Already checked at sequence load time, so should be true, but let's make sure.
	core::Size const res_per_symm_repeat(
		required_symmetry_repeats_ > 1 ? nres / required_symmetry_repeats_ : nres
	);

	//Randomly pick one of the middle residues to be the anchor residue (or one of the middle residues of the asymmetric unit if this peptide has symmetry):
	core::Size cyclization_point_start( nres );
	core::Size cyclization_point_end(1);

	if ( cyclization_type() == SCPA_terminal_disulfide ) {
		cyclization_point_start = find_last_disulf_res(pose);
		cyclization_point_end = find_first_disulf_res(pose);
	} else if ( cyclization_type() == SCPA_nterm_isopeptide_lariat || cyclization_type() == SCPA_cterm_isopeptide_lariat || cyclization_type() == SCPA_sidechain_isopeptide ) {
		find_first_and_last_isopeptide_residues( pose, cyclization_point_end, cyclization_point_start );
	} else if ( cyclization_type() == SCPA_thioether_lariat ) {
		find_first_and_last_thioether_lariat_residues( pose, cyclization_point_end, cyclization_point_start );
	} else if ( cyclization_type() == SCPA_lanthipeptide ) {
		//Ramdomly choose order of loops to close, close all loops
		//if possible randomly choose if cyc_point is upstream or downstream
		find_first_and_last_lanthipeptide_residues( cyclization_point_end, cyclization_point_start );
	}
	core::Size const anchor_res_min(
		cyclization_type() == SCPA_terminal_disulfide || cyclization_type() == SCPA_nterm_isopeptide_lariat || cyclization_type() == SCPA_cterm_isopeptide_lariat || cyclization_type() == SCPA_sidechain_isopeptide || cyclization_type() == SCPA_thioether_lariat || cyclization_type() == SCPA_lanthipeptide ?
		cyclization_point_end + 2 :
		2
	);
	core::Size const anchor_res_max(
		( cyclization_type() == SCPA_terminal_disulfide || cyclization_type() == SCPA_nterm_isopeptide_lariat || cyclization_type() == SCPA_cterm_isopeptide_lariat || cyclization_type() == SCPA_sidechain_isopeptide || cyclization_type() == SCPA_thioether_lariat || cyclization_type() == SCPA_lanthipeptide ?
		cyclization_point_start - 2 :
		( required_symmetry_repeats_ > 1 ?
		res_per_symm_repeat :
		nres - 1
		)
		)
	);

	core::Size const anchor_res( numeric::random::rg().random_range(anchor_res_min, anchor_res_max) );
	core::Size const first_loop_res( anchor_res + 1 );
	core::Size const last_loop_res( anchor_res - 1 );

	//The following should be guaranteed true:
	if ( cyclization_type() == SCPA_terminal_disulfide || cyclization_type() == SCPA_nterm_isopeptide_lariat || cyclization_type() == SCPA_cterm_isopeptide_lariat || cyclization_type() == SCPA_sidechain_isopeptide || cyclization_type() == SCPA_thioether_lariat || cyclization_type() == SCPA_lanthipeptide ) {
		debug_assert( cyclization_point_end < last_loop_res );
		debug_assert( first_loop_res < cyclization_point_start );
	} else {
		debug_assert( cyclization_point_end <= last_loop_res );
		debug_assert( first_loop_res <= cyclization_point_start );
	}
	debug_assert( last_loop_res < anchor_res );
	debug_assert( anchor_res < first_loop_res );

	//If the anchor res torsions are to be set, set 'em here:
	if ( pose->residue_type(anchor_res).is_alpha_aa() ) {
		core::Size const anchor_res_in_original( original_position(anchor_res, cyclic_offset, pose->total_residue() ) ); //Get the index of the anchor residue in the original pose (prior to any circular permutation).
		if ( user_set_alpha_dihedrals_.count(anchor_res_in_original) ) { //If this position is being set to a particular value...
			utility::vector1< core::Size > const &diheds( user_set_alpha_dihedrals_.at(anchor_res_in_original) );
			debug_assert(diheds.size() == 3); //Should be true
			pose->set_phi( anchor_res_in_original, diheds[1] );
			pose->set_psi( anchor_res_in_original, diheds[2] );
			pose->set_omega( anchor_res_in_original, diheds[3] );
			pose->update_residue_neighbors();
		}
	}

	//Randomly pick a residue to be the middle pivot residue.  Can't be first in loop, last in loop, or anchor res.
	core::Size middle_loop_res( numeric::random::rg().random_range(cyclization_point_end, cyclization_point_start-3 ) );
	if ( middle_loop_res == last_loop_res ) { middle_loop_res += 3; }
	else if ( middle_loop_res == anchor_res ) { middle_loop_res +=2; }
	else if ( middle_loop_res == first_loop_res ) { middle_loop_res +=1; }
	if ( middle_loop_res > nres ) { middle_loop_res -= nres; }

	//Create the pre-selection mover and set options.
	protocols::rosetta_scripts::ParsedProtocolOP pp( utility::pointer::make_shared< protocols::rosetta_scripts::ParsedProtocol >() );

	//Update O and H atoms at the cyclization point:
	protocols::simple_moves::DeclareBondOP update_OH( utility::pointer::make_shared< protocols::simple_moves::DeclareBond >() );
	set_up_cyclization_mover( update_OH, pose );
	if ( cyclization_type() == SCPA_n_to_c_amide_bond || cyclization_type() == SCPA_nterm_isopeptide_lariat || cyclization_type() == SCPA_cterm_isopeptide_lariat || cyclization_type() == SCPA_sidechain_isopeptide || cyclization_type() == SCPA_thioether_lariat ) {
		pp->add_step( update_OH, "Update_cyclization_point_polymer_dependent_atoms_1", nullptr );
	}

	//Filter for total hydrogen bonds:
	if ( min_genkic_hbonds_ > 0.0 ) pp->add_step( nullptr, "Total_Hbonds", total_hbond );

	//Filter out poses with oversaturated hydrogen bond acceptors.
	if ( filter_oversaturated_hbond_acceptors_ ) {
		protocols::cyclic_peptide::OversaturatedHbondAcceptorFilterOP oversat1( utility::pointer::make_shared< protocols::cyclic_peptide::OversaturatedHbondAcceptorFilter >() );
		oversat1->set_scorefxn( sfxn_default );
		oversat1->set_hbond_energy_cutoff( oversaturated_hbond_cutoff_energy_ );
		pp->add_step( nullptr, "Oversaturated_Hbond_Acceptors", oversat1 );
	}

	//If we're filtering by symmetry, do so here:
	if ( required_symmetry_repeats_ > 1 ) {
		protocols::cyclic_peptide::CycpepSymmetryFilterOP symmfilter1( utility::pointer::make_shared< protocols::cyclic_peptide::CycpepSymmetryFilter >() );
		symmfilter1->set_symm_repeats( required_symmetry_repeats_ );
		symmfilter1->set_mirror_symm( required_symmetry_mirroring_ );
		symmfilter1->set_angle_threshold( required_symmetry_angle_threshold_ );
		core::select::residue_selector::ResidueIndexSelectorOP iselector( utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >() );
		std::stringstream pep_indices("");
		pep_indices << "1-" << sequence_length();
		iselector->set_index( pep_indices.str() );
		symmfilter1->set_selector( iselector );
		pp->add_step(nullptr, "Cycpep_Symmetry_Filter_1", symmfilter1);
	}

	//If we're considering paraBBMB, add it here.
	if ( parabbmb_positions_.size() > 0 ) {
		for ( core::Size i=1, imax=parabbmb_positions_.size(); i<=imax; ++i ) { //Loop through all sets of triples of residues.
			debug_assert(parabbmb_positions_[i].size() == 2); //Should always be true.
			std::stringstream cys_indices;
			for ( core::Size j=1; j<=2; ++j ) {
				cys_indices << current_position( parabbmb_positions_[i][j], cyclic_offset, nres );
				if ( j<3 ) cys_indices << ",";
			}
			core::select::residue_selector::ResidueIndexSelectorOP index_selector( utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >() );
			index_selector->set_index( cys_indices.str() );
			protocols::cyclic_peptide::CrosslinkerMoverOP twolinker( utility::pointer::make_shared< protocols::cyclic_peptide::CrosslinkerMover >() );
			twolinker->set_residue_selector(index_selector);
			twolinker->set_linker_name("1_4_BBMB");
			twolinker->set_behaviour( true, true, true, false );
			twolinker->set_filter_behaviour( use_parabbmb_filters_, use_parabbmb_filters_, false, 0.0, parabbmb_sidechain_distance_filter_multiplier_, parabbmb_constraints_energy_filter_multiplier_ );
			twolinker->set_scorefxn( sfxn_highhbond );
			twolinker->set_sidechain_frlx_rounds(3);
			std::stringstream movername;
			movername << "paraBBMB_link_" << i;
			pp->add_step( twolinker, movername.str(), nullptr );
		}
	}

	//If we're considering TBMB, add it here.
	if ( tbmb_positions_.size() > 0 ) {
		for ( core::Size i=1, imax=tbmb_positions_.size(); i<=imax; ++i ) { //Loop through all sets of triples of residues.
			debug_assert(tbmb_positions_[i].size() == 3); //Should always be true.
			std::stringstream cys_indices;
			for ( core::Size j=1; j<=3; ++j ) {
				cys_indices << current_position( tbmb_positions_[i][j], cyclic_offset, nres );
				if ( j<3 ) cys_indices << ",";
			}
			core::select::residue_selector::ResidueIndexSelectorOP index_selector( utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >() );
			index_selector->set_index( cys_indices.str() );
			protocols::cyclic_peptide::CrosslinkerMoverOP threelinker( utility::pointer::make_shared< protocols::cyclic_peptide::CrosslinkerMover >() );
			threelinker->set_residue_selector(index_selector);
			threelinker->set_linker_name("TBMB");
			threelinker->set_behaviour( true, true, true, false );
			threelinker->set_filter_behaviour( use_tbmb_filters_, use_tbmb_filters_, false, 0.0, tbmb_sidechain_distance_filter_multiplier_, tbmb_constraints_energy_filter_multiplier_ );
			threelinker->set_scorefxn( sfxn_highhbond );
			threelinker->set_sidechain_frlx_rounds(3);
			std::stringstream movername;
			movername << "TBMB_link_" << i;
			pp->add_step( threelinker, movername.str(), nullptr );
		}
	}

	//If we're considering TMA, add it here.
	if ( tma_positions_.size() > 0 ) {
		for ( core::Size i(1), imax(tma_positions_.size()); i<=imax; ++i ) {
			debug_assert(tma_positions_[i].size() == 3); //Should be true always.
			std::stringstream tma_indices;
			for ( core::Size j=1; j<=3; ++j ) {
				tma_indices << current_position( tma_positions_[i][j], cyclic_offset, nres );
				if ( j<3 ) tma_indices << ",";
			}
			core::select::residue_selector::ResidueIndexSelectorOP index_selector( utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >() );
			index_selector->set_index( tma_indices.str() );
			protocols::cyclic_peptide::CrosslinkerMoverOP threelinker( utility::pointer::make_shared< protocols::cyclic_peptide::CrosslinkerMover >() );
			threelinker->set_residue_selector(index_selector);
			threelinker->set_linker_name("TMA");
			threelinker->set_behaviour( true, true, true, false );
			threelinker->set_filter_behaviour( use_tma_filters_, use_tma_filters_, false, 0.0, tma_sidechain_distance_filter_multiplier_, tma_constraints_energy_filter_multiplier_ );
			threelinker->set_scorefxn( sfxn_highhbond );
			threelinker->set_sidechain_frlx_rounds(3);
			std::stringstream movername;
			movername << "TMA_link_" << i;
			pp->add_step( threelinker, movername.str(), nullptr );
		}
	}

	//If we're considering trigonal pyramidal metal crosslinks, add 'em here:
	if ( trigonal_pyramidal_metal_positions_.size() > 0 ) {
		for ( core::Size i(1), imax(trigonal_pyramidal_metal_positions_.size()); i<=imax; ++i ) {
			utility::fixedsizearray1< core::Size, 3 > const & resnums( trigonal_pyramidal_metal_positions_[i].first );
			core::select::residue_selector::ResidueIndexSelectorOP index_selector( utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >() );
			for ( core::Size j(1); j<=3; ++j ) index_selector->append_index( current_position( resnums[j], cyclic_offset, nres ) );
			protocols::cyclic_peptide::CrosslinkerMoverOP metallinker( utility::pointer::make_shared< protocols::cyclic_peptide::CrosslinkerMover >() );
			metallinker->set_residue_selector(index_selector);
			metallinker->set_linker_name("trigonal_pyramidal_metal");
			metallinker->set_metal_type(trigonal_pyramidal_metal_positions_[i].second);
			metallinker->set_behaviour(true, true, true, false);
			metallinker->set_filter_behaviour(use_trigonal_pyramidal_metal_filters_, use_trigonal_pyramidal_metal_filters_, false, 0.0, trigonal_pyramidal_metal_sidechain_distance_filter_multiplier_, trigonal_pyramidal_metal_constraints_energy_filter_multiplier_);
			metallinker->set_scorefxn( sfxn_highhbond );
			metallinker->set_sidechain_frlx_rounds(3);
			std::stringstream movername;
			movername << "trigonal_pyramidal_metal_link_" << i;
			pp->add_step( metallinker, movername.str(), nullptr );
		}
	}

	//If we're considering trigonal planar metal crosslinks, add 'em here:
	if ( trigonal_planar_metal_positions_.size() > 0 ) {
		for ( core::Size i(1), imax(trigonal_planar_metal_positions_.size()); i<=imax; ++i ) {
			utility::fixedsizearray1< core::Size, 3 > const & resnums( trigonal_planar_metal_positions_[i].first );
			core::select::residue_selector::ResidueIndexSelectorOP index_selector( utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >() );
			for ( core::Size j(1); j<=3; ++j ) index_selector->append_index( current_position( resnums[j], cyclic_offset, nres ) );
			protocols::cyclic_peptide::CrosslinkerMoverOP metallinker( utility::pointer::make_shared< protocols::cyclic_peptide::CrosslinkerMover >() );
			metallinker->set_residue_selector(index_selector);
			metallinker->set_linker_name("trigonal_planar_metal");
			metallinker->set_metal_type(trigonal_planar_metal_positions_[i].second);
			metallinker->set_behaviour(true, true, true, false);
			metallinker->set_filter_behaviour(use_trigonal_planar_metal_filters_, use_trigonal_planar_metal_filters_, false, 0.0, trigonal_planar_metal_sidechain_distance_filter_multiplier_, trigonal_planar_metal_constraints_energy_filter_multiplier_);
			metallinker->set_scorefxn( sfxn_highhbond );
			metallinker->set_sidechain_frlx_rounds(3);
			std::stringstream movername;
			movername << "trigonal_planar_metal_link_" << i;
			pp->add_step( metallinker, movername.str(), nullptr );
		}
	}

	//If we're considering square pyramidal metal crosslinks, add 'em here:
	if ( square_pyramidal_metal_positions_.size() > 0 ) {
		for ( core::Size i(1), imax(square_pyramidal_metal_positions_.size()); i<=imax; ++i ) {
			utility::fixedsizearray1< core::Size, 5 > const & resnums( square_pyramidal_metal_positions_[i].first );
			core::select::residue_selector::ResidueIndexSelectorOP index_selector( utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >() );
			for ( core::Size j(1); j<=5; ++j ) index_selector->append_index( current_position( resnums[j], cyclic_offset, nres ) );
			protocols::cyclic_peptide::CrosslinkerMoverOP metallinker( utility::pointer::make_shared< protocols::cyclic_peptide::CrosslinkerMover >() );
			metallinker->set_residue_selector(index_selector);
			metallinker->set_linker_name("square_pyramidal_metal");
			metallinker->set_metal_type(square_pyramidal_metal_positions_[i].second);
			metallinker->set_behaviour(true, true, true, false);
			metallinker->set_filter_behaviour(use_square_pyramidal_metal_filters_, use_square_pyramidal_metal_filters_, false, 0.0, square_pyramidal_metal_sidechain_distance_filter_multiplier_, square_pyramidal_metal_constraints_energy_filter_multiplier_);
			metallinker->set_scorefxn( sfxn_highhbond );
			metallinker->set_sidechain_frlx_rounds(3);
			std::stringstream movername;
			movername << "square_pyramidal_metal_link_" << i;
			pp->add_step( metallinker, movername.str(), nullptr );
		}
	}

	//If we're considering square planar metal crosslinks, add 'em here:
	if ( square_planar_metal_positions_.size() > 0 ) {
		for ( core::Size i(1), imax(square_planar_metal_positions_.size()); i<=imax; ++i ) {
			utility::fixedsizearray1< core::Size, 4 > const & resnums( square_planar_metal_positions_[i].first );
			core::select::residue_selector::ResidueIndexSelectorOP index_selector( utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >() );
			for ( core::Size j(1); j<=4; ++j ) index_selector->append_index( current_position( resnums[j], cyclic_offset, nres ) );
			protocols::cyclic_peptide::CrosslinkerMoverOP metallinker( utility::pointer::make_shared< protocols::cyclic_peptide::CrosslinkerMover >() );
			metallinker->set_residue_selector(index_selector);
			metallinker->set_linker_name("square_planar_metal");
			metallinker->set_metal_type(square_planar_metal_positions_[i].second);
			metallinker->set_behaviour(true, true, true, false);
			metallinker->set_filter_behaviour(use_square_planar_metal_filters_, use_square_planar_metal_filters_, false, 0.0, square_planar_metal_sidechain_distance_filter_multiplier_, square_planar_metal_constraints_energy_filter_multiplier_);
			metallinker->set_scorefxn( sfxn_highhbond );
			metallinker->set_sidechain_frlx_rounds(3);
			std::stringstream movername;
			movername << "square_planar_metal_link_" << i;
			pp->add_step( metallinker, movername.str(), nullptr );
		}
	}

	//If we're considering tetrahedral metal crosslinks, add 'em here:
	if ( tetrahedral_metal_positions_.size() > 0 ) {
		for ( core::Size i(1), imax(tetrahedral_metal_positions_.size()); i<=imax; ++i ) {
			utility::fixedsizearray1< core::Size, 4 > const & resnums( tetrahedral_metal_positions_[i].first );
			core::select::residue_selector::ResidueIndexSelectorOP index_selector( utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >() );
			for ( core::Size j(1); j<=4; ++j ) index_selector->append_index( current_position( resnums[j], cyclic_offset, nres ) );
			protocols::cyclic_peptide::CrosslinkerMoverOP metallinker( utility::pointer::make_shared< protocols::cyclic_peptide::CrosslinkerMover >() );
			metallinker->set_residue_selector(index_selector);
			metallinker->set_linker_name("tetrahedral_metal");
			metallinker->set_metal_type(tetrahedral_metal_positions_[i].second);
			metallinker->set_behaviour(true, true, true, false);
			metallinker->set_filter_behaviour(use_tetrahedral_metal_filters_, use_tetrahedral_metal_filters_, false, 0.0, tetrahedral_metal_sidechain_distance_filter_multiplier_, tetrahedral_metal_constraints_energy_filter_multiplier_);
			metallinker->set_scorefxn( sfxn_highhbond );
			metallinker->set_sidechain_frlx_rounds(3);
			std::stringstream movername;
			movername << "tetrahedral_metal_link_" << i;
			pp->add_step( metallinker, movername.str(), nullptr );
		}
	}

	//If we're considering octahedral metal crosslinks, add 'em here:
	if ( octahedral_metal_positions_.size() > 0 ) {
		for ( core::Size i(1), imax(octahedral_metal_positions_.size()); i<=imax; ++i ) {
			utility::fixedsizearray1< core::Size, 6 > const & resnums( octahedral_metal_positions_[i].first );
			core::select::residue_selector::ResidueIndexSelectorOP index_selector( utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >() );
			for ( core::Size j(1); j<=6; ++j ) index_selector->append_index( current_position( resnums[j], cyclic_offset, nres ) );
			protocols::cyclic_peptide::CrosslinkerMoverOP metallinker( utility::pointer::make_shared< protocols::cyclic_peptide::CrosslinkerMover >() );
			metallinker->set_residue_selector(index_selector);
			metallinker->set_linker_name("octahedral_metal");
			metallinker->set_metal_type(octahedral_metal_positions_[i].second);
			metallinker->set_behaviour(true, true, true, false);
			metallinker->set_filter_behaviour(use_octahedral_metal_filters_, use_octahedral_metal_filters_, false, 0.0, octahedral_metal_sidechain_distance_filter_multiplier_, octahedral_metal_constraints_energy_filter_multiplier_);
			metallinker->set_scorefxn( sfxn_highhbond );
			metallinker->set_sidechain_frlx_rounds(3);
			std::stringstream movername;
			movername << "octahedral_metal_link_" << i;
			pp->add_step( metallinker, movername.str(), nullptr );
		}
	}

	core::Size disulf_count(0);
	//If we're considering disulfides, add the TryDisulfPermutations mover and a filter to the ParsedProtocol:
	if ( try_all_disulfides_ ) {
		protocols::cyclic_peptide::TryDisulfPermutationsOP trydisulf( utility::pointer::make_shared< protocols::cyclic_peptide::TryDisulfPermutations >() ); //Default settings should be fine.
		core::Size disulf_res_count(0);
		core::Size const first_disulf_res( cyclization_type() == SCPA_terminal_disulfide ? find_first_disulf_res( pose ) : 0 );
		core::Size const last_disulf_res( cyclization_type() == SCPA_terminal_disulfide ? find_last_disulf_res( pose ) : 0 );
		if ( cyclization_type() == SCPA_terminal_disulfide ) {
			core::select::residue_selector::ResidueIndexSelectorOP sel( utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >() );
			std::stringstream indices;
			indices << first_disulf_res << "," << last_disulf_res ;
			sel->set_index( indices.str() );
			core::select::residue_selector::NotResidueSelectorOP notsel( utility::pointer::make_shared< core::select::residue_selector::NotResidueSelector >() );
			notsel->set_residue_selector(sel);
			trydisulf->set_selector( notsel );
		}
		for ( core::Size ir=1, irmax=pose->size(); ir<=irmax; ++ir ) {
			if ( pose->residue(ir).type().get_disulfide_atom_name() != "NONE" ) ++disulf_res_count; //Count disulfide-forming residues in the pose.
		}
		disulf_count = disulf_res_count / 2; //Div operator -- gives correct number of disulfides even in odd disulfide-forming residue case.
		protocols::score_filters::ScoreTypeFilterOP disulf_filter1( utility::pointer::make_shared< protocols::score_filters::ScoreTypeFilter >( sfxn_highhbond, core::scoring::dslf_fa13, disulf_energy_cutoff_prerelax_ * static_cast<core::Real>(disulf_count) ) );
		pp->add_step( trydisulf, "Try_Disulfide_Permutations", disulf_filter1 );
	}

	//Add the FastRelax with high hbond weight to the pre-selection parsed protocol.
	if ( design_peptide_ ) {
		if ( L_alpha_comp_file_exists_ ) {
			protocols::aa_composition::AddCompositionConstraintMoverOP L_alpha_cst( utility::pointer::make_shared< protocols::aa_composition::AddCompositionConstraintMover >() );
			core::select::residue_selector::BinSelectorOP select_L_alpha( utility::pointer::make_shared< core::select::residue_selector::BinSelector >() );
			select_L_alpha->set_bin_name("A");
			select_L_alpha->set_bin_params_file_name(abba_bins_binfile_);
			select_L_alpha->initialize_and_check();
			L_alpha_cst->create_constraint_from_file_contents( comp_file_contents_L_alpha_ );
			L_alpha_cst->add_residue_selector( select_L_alpha );
			pp->add_step( L_alpha_cst, "Add_L_Alpha_AACompositionConstraints", nullptr );
		}
		if ( D_alpha_comp_file_exists_ ) {
			protocols::aa_composition::AddCompositionConstraintMoverOP D_alpha_cst( utility::pointer::make_shared< protocols::aa_composition::AddCompositionConstraintMover >() );
			core::select::residue_selector::BinSelectorOP select_D_alpha( utility::pointer::make_shared< core::select::residue_selector::BinSelector >() );
			select_D_alpha->set_bin_name("Aprime");
			select_D_alpha->set_bin_params_file_name(abba_bins_binfile_);
			select_D_alpha->initialize_and_check();
			D_alpha_cst->create_constraint_from_file_contents( comp_file_contents_D_alpha_ );
			D_alpha_cst->add_residue_selector( select_D_alpha );
			pp->add_step( D_alpha_cst, "Add_D_Alpha_AACompositionConstraints", nullptr );
		}
		if ( L_beta_comp_file_exists_ ) {
			protocols::aa_composition::AddCompositionConstraintMoverOP L_beta_cst( utility::pointer::make_shared< protocols::aa_composition::AddCompositionConstraintMover >() );
			core::select::residue_selector::BinSelectorOP select_L_beta( utility::pointer::make_shared< core::select::residue_selector::BinSelector >() );
			select_L_beta->set_bin_name("B");
			select_L_beta->set_bin_params_file_name(abba_bins_binfile_);
			select_L_beta->initialize_and_check();
			L_beta_cst->create_constraint_from_file_contents( comp_file_contents_L_beta_ );
			L_beta_cst->add_residue_selector( select_L_beta );
			pp->add_step( L_beta_cst, "Add_L_Beta_AACompositionConstraints", nullptr );
		}
		if ( D_beta_comp_file_exists_ ) {
			protocols::aa_composition::AddCompositionConstraintMoverOP D_beta_cst( utility::pointer::make_shared< protocols::aa_composition::AddCompositionConstraintMover >() );
			core::select::residue_selector::BinSelectorOP select_D_beta( utility::pointer::make_shared< core::select::residue_selector::BinSelector >() );
			select_D_beta->set_bin_name("Bprime");
			select_D_beta->set_bin_params_file_name(abba_bins_binfile_);
			select_D_beta->initialize_and_check();
			D_beta_cst->create_constraint_from_file_contents( comp_file_contents_D_beta_ );
			D_beta_cst->add_residue_selector( select_D_beta );
			pp->add_step( D_beta_cst, "Add_D_Beta_AACompositionConstraints", nullptr );
		}

		if ( fast_relax_rounds_ > 0 ) {
			protocols::denovo_design::movers::FastDesignOP fdes( utility::pointer::make_shared< protocols::denovo_design::movers::FastDesign >(sfxn_highhbond, fast_relax_rounds_) );
			set_up_design_taskoperations( fdes, cyclic_offset, pose->size(), pose );
			pp->add_step( fdes, "High_Hbond_FastDesign", nullptr );
		}
		if ( angle_relax_rounds() > 0 ) {
			protocols::denovo_design::movers::FastDesignOP fdes2( utility::pointer::make_shared< protocols::denovo_design::movers::FastDesign >(sfxn_highhbond_cart, angle_relax_rounds()) );
			fdes2->minimize_bond_angles(true);
			set_up_design_taskoperations( fdes2, cyclic_offset, pose->size(), pose );
			pp->add_step( fdes2, "High_Hbond_FastDesign_angle_relax", nullptr );
		}
		if ( angle_length_relax_rounds() > 0 ) {
			protocols::denovo_design::movers::FastDesignOP fdes3( utility::pointer::make_shared< protocols::denovo_design::movers::FastDesign >(sfxn_highhbond_cart, angle_length_relax_rounds()) );
			fdes3->minimize_bond_angles(true);
			fdes3->minimize_bond_lengths(true);
			set_up_design_taskoperations( fdes3, cyclic_offset, pose->size(), pose );
			pp->add_step( fdes3, "High_Hbond_FastDesign_angle_length_relax", nullptr );
		}
		if ( cartesian_relax_rounds() > 0 ) {
			protocols::denovo_design::movers::FastDesignOP fdes4( utility::pointer::make_shared< protocols::denovo_design::movers::FastDesign >(sfxn_highhbond_cart, cartesian_relax_rounds()) );
			fdes4->cartesian(true);
			set_up_design_taskoperations( fdes4, cyclic_offset, pose->size(), pose );
			pp->add_step( fdes4, "High_Hbond_FastDesign_Cartesian_relax", nullptr );
		}

		protocols::aa_composition::ClearCompositionConstraintsMoverOP clear_aacomp_cst( utility::pointer::make_shared< protocols::aa_composition::ClearCompositionConstraintsMover >() );
		pp->add_step( clear_aacomp_cst, "Clear_AACompositionConstraints", nullptr );
	} else {
		if ( fast_relax_rounds_ > 0 ) {
			protocols::relax::FastRelaxOP frlx( utility::pointer::make_shared< protocols::relax::FastRelax >(sfxn_highhbond, fast_relax_rounds_) );
			pp->add_step( frlx, "High_Hbond_FastRelax", nullptr );
		}
		if ( angle_relax_rounds() > 0 ) {
			protocols::relax::FastRelaxOP frlx2( utility::pointer::make_shared< protocols::relax::FastRelax >(sfxn_highhbond_cart, angle_relax_rounds() ) );
			frlx2->minimize_bond_angles(true);
			pp->add_step( frlx2, "High_Hbond_FastRelax_angles", nullptr );
		}
		if ( angle_length_relax_rounds() > 0 ) {
			protocols::relax::FastRelaxOP frlx3( utility::pointer::make_shared< protocols::relax::FastRelax >(sfxn_highhbond_cart, angle_length_relax_rounds() ) );
			frlx3->minimize_bond_angles(true);
			pp->add_step( frlx3, "High_Hbond_FastRelax_angles_bondlengths", nullptr );
		}
		if ( cartesian_relax_rounds() > 0 ) {
			protocols::relax::FastRelaxOP frlx4( utility::pointer::make_shared< protocols::relax::FastRelax >(sfxn_highhbond_cart, cartesian_relax_rounds() ) );
			frlx4->minimize_bond_angles(true);
			pp->add_step( frlx4, "High_Hbond_FastRelax_Cartesian", nullptr );
		}
	}

	//Update O and H atoms at the cyclization point:
	if ( cyclization_type() == SCPA_n_to_c_amide_bond || cyclization_type() == SCPA_nterm_isopeptide_lariat || cyclization_type() == SCPA_cterm_isopeptide_lariat || cyclization_type() == SCPA_sidechain_isopeptide || cyclization_type() == SCPA_thioether_lariat ) {
		pp->add_step( update_OH, "Update_cyclization_point_polymer_dependent_atoms_2", nullptr );
	}

	//Add more stringent disulfide filtering post-relax:
	if ( try_all_disulfides_ ) {
		protocols::score_filters::ScoreTypeFilterOP disulf_filter2( utility::pointer::make_shared< protocols::score_filters::ScoreTypeFilter >( sfxn_highhbond, core::scoring::dslf_fa13, disulf_energy_cutoff_postrelax_ * static_cast<core::Real>(disulf_count) ) );
		pp->add_step( nullptr, "Postrelax_disulfide_filter", disulf_filter2 );
	}

	if ( filter_oversaturated_hbond_acceptors_ ) {
		protocols::cyclic_peptide::OversaturatedHbondAcceptorFilterOP oversat2( utility::pointer::make_shared< protocols::cyclic_peptide::OversaturatedHbondAcceptorFilter >() );
		oversat2->set_scorefxn( sfxn_default );
		oversat2->set_hbond_energy_cutoff( oversaturated_hbond_cutoff_energy_ );
		pp->add_step( nullptr, "Postrelax_Oversaturated_Hbond_Acceptors", oversat2 );
	}

	//If we're filtering by symmetry, do so again here:
	if ( required_symmetry_repeats_ > 1 ) {
		protocols::cyclic_peptide::CycpepSymmetryFilterOP symmfilter2( utility::pointer::make_shared< protocols::cyclic_peptide::CycpepSymmetryFilter >() );
		symmfilter2->set_symm_repeats( required_symmetry_repeats_ );
		symmfilter2->set_mirror_symm( required_symmetry_mirroring_ );
		symmfilter2->set_angle_threshold( required_symmetry_angle_threshold_ );
		core::select::residue_selector::ResidueIndexSelectorOP iselector( utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >() );
		std::stringstream pep_indices("");
		pep_indices << "1-" << sequence_length();
		iselector->set_index( pep_indices.str() );
		symmfilter2->set_selector( iselector );
		pp->add_step(nullptr, "Cycpep_Symmetry_Filter_2", symmfilter2);
	}

	//Create the mover and set options:
	GeneralizedKICOP genkic( utility::pointer::make_shared< GeneralizedKIC >() );
	genkic->set_selector_type( lowest_energy_selector );
	genkic->set_closure_attempts( genkic_closure_attempts_ );
	genkic->set_min_solution_count( genkic_min_solution_count_ );
	genkic->set_selector_scorefunction( sfxn_highhbond );
	genkic->set_preselection_mover(pp);
	genkic->set_correct_polymer_dependent_atoms(true);

	//If we're using BOINC graphics, let the GenKIC mover update the graphics with a "ghost" of the current
	//conformation being sampled:
#ifdef BOINC_GRAPHICS
	genkic->set_attach_boinc_ghost_observer(true);
#endif


	//Define the loop residues:
	for ( core::Size i=first_loop_res; i<=cyclization_point_start; ++i ) { genkic->add_loop_residue(i); }
	for ( core::Size i=cyclization_point_end; i<=last_loop_res; ++i ) { genkic->add_loop_residue(i); }

	//Add tail residues, if relevant:
	if ( cyclization_type() == SCPA_terminal_disulfide || cyclization_type() == SCPA_nterm_isopeptide_lariat || cyclization_type() == SCPA_cterm_isopeptide_lariat || cyclization_type() == SCPA_sidechain_isopeptide || cyclization_type() == SCPA_thioether_lariat || cyclization_type() == SCPA_lanthipeptide ) {
		if ( cyclization_point_end > 1 ) {
			for ( core::Size i=cyclization_point_end-1; i>=1; --i ) { genkic->add_tail_residue(i); }
		}
		if ( cyclization_point_start < nres ) {
			for ( core::Size i=cyclization_point_start+1; i<=nres; ++i ) { genkic->add_tail_residue(i); }
		}
	}

	//Set pivots:
	std::string at1(""), at2(""), at3("");
	if ( pose->residue(first_loop_res).type().is_alpha_aa() || pose->residue(first_loop_res).is_peptoid() ) { at1="CA"; }
	else if ( pose->residue_type(first_loop_res).is_beta_aa() || pose->residue_type(first_loop_res).is_oligourea() ) { at1="CM"; }
	else if ( pose->residue(first_loop_res).type().is_gamma_aa() ) { at1="C3"; }
	else if ( pose->residue(first_loop_res).type().is_aramid() ) {
		if ( pose->residue(first_loop_res).type().is_pre_methylene_ortho_aramid() ) { at1="CA"; }
		else if ( pose->residue(first_loop_res).type().is_pre_methylene_meta_aramid() ) { at1="CA"; }
		else if ( pose->residue(first_loop_res).type().is_pre_methylene_para_aramid() ) { at1="CA"; }
		else if ( pose->residue(first_loop_res).type().is_post_methylene_ortho_aramid() ) { at1="Ca"; }
		else if ( pose->residue(first_loop_res).type().is_post_methylene_meta_aramid() ) { at1="CA"; }
		else if ( pose->residue(first_loop_res).type().is_post_methylene_para_aramid() ) { at1="CA"; }
		else if ( pose->residue(first_loop_res).type().is_pre_methylene_post_methylene_ortho_aramid() ) { at1="CA"; }
		else if ( pose->residue(first_loop_res).type().is_pre_methylene_post_methylene_meta_aramid() ) { at1="CA"; }
		else if ( pose->residue(first_loop_res).type().is_pre_methylene_post_methylene_para_aramid() ) { at1="CA"; }
	} else { utility_exit_with_message( "Unrecognized residue type at loop start.  Currently, this app only works with alpha, beta, and gamma amino acids, peptoids, and oligoureas. Polyaramids with a free methylene work, too." ); }
	if ( cyclization_type() == SCPA_terminal_disulfide && (middle_loop_res == cyclization_point_start || middle_loop_res == cyclization_point_end )  ) { at2 = pose->residue_type(middle_loop_res).get_disulfide_atom_name(); }
	else if ( pose->residue(middle_loop_res).type().is_alpha_aa() || pose->residue(middle_loop_res).is_peptoid() ) { at2="CA"; }
	else if ( pose->residue_type(middle_loop_res).is_beta_aa() || pose->residue_type(middle_loop_res).is_oligourea() ) { at2="CM"; }
	else if ( pose->residue(middle_loop_res).type().is_gamma_aa() ) { at2="C3"; }
	else if ( pose->residue(middle_loop_res).type().is_aramid() ) {
		if ( pose->residue(middle_loop_res).type().is_pre_methylene_ortho_aramid() ) { at2="CA"; }
		else if ( pose->residue(middle_loop_res).type().is_pre_methylene_meta_aramid() ) { at2="CA"; }
		else if ( pose->residue(middle_loop_res).type().is_pre_methylene_para_aramid() ) { at2="CA"; }
		else if ( pose->residue(middle_loop_res).type().is_post_methylene_ortho_aramid() ) { at2="Ca"; }
		else if ( pose->residue(middle_loop_res).type().is_post_methylene_meta_aramid() ) { at2="CA"; }
		else if ( pose->residue(middle_loop_res).type().is_post_methylene_para_aramid() ) { at2="CA"; }
		else if ( pose->residue(middle_loop_res).type().is_pre_methylene_post_methylene_ortho_aramid() ) { at2="CA"; }
		else if ( pose->residue(middle_loop_res).type().is_pre_methylene_post_methylene_meta_aramid() ) { at2="CA"; }
		else if ( pose->residue(middle_loop_res).type().is_pre_methylene_post_methylene_para_aramid() ) { at2="CA"; }
	} else { utility_exit_with_message( "Unrecognized residue type at loop midpoint.  Currently, this app only works with alpha, beta, and gamma amino acids, peptoids, and oligoureas. Polyaramids with a free methylene work, too." ); }
	if ( pose->residue(last_loop_res).type().is_alpha_aa() || pose->residue(last_loop_res).is_peptoid() ) { at3="CA"; }
	else if ( pose->residue_type(last_loop_res).is_beta_aa() || pose->residue_type(last_loop_res).is_oligourea() ) { at3="CM"; }
	else if ( pose->residue(last_loop_res).type().is_gamma_aa() ) { at3="C3"; }
	else if ( pose->residue(middle_loop_res).type().is_aramid() ) {
		if ( pose->residue(last_loop_res).type().is_pre_methylene_ortho_aramid() ) { at3="CA"; }
		else if ( pose->residue(last_loop_res).type().is_pre_methylene_meta_aramid() ) { at3="CA"; }
		else if ( pose->residue(last_loop_res).type().is_pre_methylene_para_aramid() ) { at3="CA"; }
		else if ( pose->residue(last_loop_res).type().is_post_methylene_ortho_aramid() ) { at3="Ca"; }
		else if ( pose->residue(last_loop_res).type().is_post_methylene_meta_aramid() ) { at3="CA"; }
		else if ( pose->residue(last_loop_res).type().is_post_methylene_para_aramid() ) { at3="CA"; }
		else if ( pose->residue(last_loop_res).type().is_pre_methylene_post_methylene_ortho_aramid() ) { at3="CA"; }
		else if ( pose->residue(last_loop_res).type().is_pre_methylene_post_methylene_meta_aramid() ) { at3="CA"; }
		else if ( pose->residue(last_loop_res).type().is_pre_methylene_post_methylene_para_aramid() ) { at3="CA"; }
	} else { utility_exit_with_message( "Unrecognized residue type at loop midpoint.  Currently, this app only works with alpha, beta, and gamma amino acids, peptoids, and oligoureas. Polyaramids with a free methylene work, too." ); }
	genkic->set_pivot_atoms( first_loop_res, at1, middle_loop_res, at2, last_loop_res, at3 );

	//Close the bond:
	add_closebond_logic( pose, cyclization_point_start, cyclization_point_end, genkic );

	//Add backbone perturbers:
	for ( core::Size i=anchor_res+1; i!=anchor_res; ++i ) {
		if ( i > cyclization_point_start ) i = cyclization_point_end;

		if ( i==anchor_res ) continue; //Can't perturb the anchor residue.
		if ( cyclization_type() == SCPA_terminal_disulfide || cyclization_type() == SCPA_sidechain_isopeptide || cyclization_type() == SCPA_thioether_lariat || cyclization_type() == SCPA_lanthipeptide ) { //If we're cyclizing through a disulfide, don't randomize by rama
			if ( i==cyclization_point_start || i==cyclization_point_end ) {
				continue;
			}
		} else if ( cyclization_type() == SCPA_nterm_isopeptide_lariat && i == cyclization_point_end ) { continue; }
		else if ( cyclization_type() == SCPA_cterm_isopeptide_lariat && i == cyclization_point_start ) { continue; }
		else if ( cyclization_type() == SCPA_thioether_lariat && i == cyclization_point_end ) { continue; }

		if ( pose->residue_type(i).is_alpha_aa() || pose->residue_type(i).is_oligourea() || pose->residue_type(i).is_peptoid() || pose->residue_type(i).is_aramid() || pose->residue_type(i).is_beta_aa() ) {
			core::Size const res_in_original( original_position(i, cyclic_offset, pose->size() ) ); //Get the index of this position in the original pose (prior to any circular permutation).
			if ( user_set_alpha_dihedrals_.count(res_in_original) ) { //If this position is being set to a particular value...
				runtime_assert_string_msg( pose->residue_type(i).is_alpha_aa() || pose->residue_type(i).is_peptoid(), "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication: dihedral setting currently only works for alpha-amino acids and peptoids." );
				genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::set_dihedral );
				core::id::NamedAtomID Natom( "N", i );
				core::id::NamedAtomID CAatom( "CA", i );
				core::id::NamedAtomID Catom( "C", i );
				core::Size const nextres( i==nres ? 1 : i+1 );
				core::id::NamedAtomID Nnextatom( "N", nextres );
				utility::vector1< core::id::NamedAtomID > const phivect = { Natom, CAatom };
				utility::vector1< core::id::NamedAtomID > const psivect = { CAatom, Catom };
				utility::vector1< core::id::NamedAtomID > const omegavect = { Catom, Nnextatom };
				genkic->add_atomset_to_perturber_atomset_list(phivect);
				genkic->add_atomset_to_perturber_atomset_list(psivect);
				if ( nextres != anchor_res ) genkic->add_atomset_to_perturber_atomset_list(omegavect);
				genkic->add_value_to_perturber_value_list( user_set_alpha_dihedrals_.at(res_in_original)[1] );
				genkic->add_value_to_perturber_value_list( user_set_alpha_dihedrals_.at(res_in_original)[2] );
				if ( nextres != anchor_res ) genkic->add_value_to_perturber_value_list( user_set_alpha_dihedrals_.at(res_in_original)[3] );
				if ( user_set_dihedral_perturbation_ ) {
					genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::perturb_dihedral );
					genkic->add_atomset_to_perturber_atomset_list(phivect);
					genkic->add_atomset_to_perturber_atomset_list(psivect);
					//if ( nextres != anchor_res ) genkic->add_atomset_to_perturber_atomset_list(omegavect);
					genkic->add_value_to_perturber_value_list( user_set_dihedral_perturbation_ );
				}
			} else { //If this position is not set, randomize it.
				if ( position_backbone_is_randomizable( i ) ) {
					if ( custom_rama_table_defined( res_in_original ) ) { //If there is a custom rama table defined for sampling at this position, use it.
						runtime_assert_string_msg( pose->residue_type(i).is_alpha_aa(), "Error in SimpleCycpepPredictApplication: the use of custom old-style Ramachandran maps is limited to alpha-amino acids." );
						genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::randomize_alpha_backbone_by_rama );
						genkic->add_residue_to_perturber_residue_list(i);
						genkic->set_perturber_custom_rama_table( rama_table_type_by_res( res_in_original ) );
					} else {
						if ( use_rama_prepro_for_sampling() && default_rama_table_type() == core::scoring::unknown_ramatable_type ) {
							genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::randomize_backbone_by_rama_prepro );
							genkic->add_residue_to_perturber_residue_list(i);
						} else {
							runtime_assert_string_msg( pose->residue_type(i).is_alpha_aa(), "Error in SimpleCycpepPredictApplication: the use of old-style Ramachandran maps is limited to alpha-amino acids." );
							genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::randomize_alpha_backbone_by_rama );
							genkic->add_residue_to_perturber_residue_list(i);
							if ( default_rama_table_type() != core::scoring::unknown_ramatable_type ) {
								genkic->set_perturber_custom_rama_table( default_rama_table_type() );
							}
						}
					}
				}
			}
		} else {
			//TODO Randomize mainchain torsions here for beta- and gamma-amino acids.
			utility_exit_with_message( "Handling of gamma-amino acids, peptoids, and other exotic backbones in setup of the genKIC perturber in the simple_cycpep_predict app has not yet been written.  TODO." );
		}
	}

	//Additional perturber: sampling cis proline.  Must be after the other perturbers.
	if ( sample_cis_pro() ) {
		for ( core::Size i=anchor_res+1; i!=anchor_res; ++i ) {
			if ( i > cyclization_point_start ) i = cyclization_point_end;

			if ( cyclization_type() == SCPA_n_to_c_amide_bond && i==1 && nres==anchor_res ) continue; //Can't perturb the anchor residue.
			if ( i-1 == anchor_res ) continue; //Can't perturb the anchor residue.
			if ( i == cyclization_point_end && ( cyclization_type() == SCPA_nterm_isopeptide_lariat || cyclization_type() == SCPA_cterm_isopeptide_lariat || cyclization_type() == SCPA_terminal_disulfide || cyclization_type() == SCPA_sidechain_isopeptide || cyclization_type() == SCPA_thioether_lariat || cyclization_type() == SCPA_lanthipeptide ) ) {
				continue; //Can't perturb preceding bond if it's not in the loop.
			}
			if ( pose->residue_type(i).aa() == core::chemical::aa_pro || pose->residue_type(i).aa() == core::chemical::aa_dpr || pose->residue_type(i).is_peptoid() || pose->residue_type(i).is_n_methylated() ) {
				genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::sample_cis_peptide_bond );
				genkic->add_value_to_perturber_value_list( sample_cis_pro_frequency() );
				genkic->add_residue_to_perturber_residue_list( i==1 ? nres : i-1 ); //The residue PRECEDING the proline is perturbed
			}
		}
	}

	//Additional perturbers: randomize backbone bond lengths and bond angles.
	if ( bondlength_perturbation_magnitude_ /*If the perturbation is nonzero...*/ ) {
		add_bondlength_perturbation( *genkic, bondlength_perturbation_magnitude_, *pose, anchor_res );
	}
	if ( bondangle_perturbation_magnitude_ /*If the perturbation is nonzero...*/ ) {
		add_bondangle_perturbation( *genkic, bondangle_perturbation_magnitude_, *pose, anchor_res );
	}

	//Very last sampling: copy symmetric positions.  Must be after ALL other perturbers.
	if ( required_symmetry_repeats_ > 1 ) {
		for ( core::Size i(1), imax(pose->total_residue()); i<=imax; ++i ) { //This time, we loop through all residues in order.
			//Since symmetric sampling is only compatible with N-to-C cyclization, we don't need to worry
			//about tail residues.  We only need to exclude the anchor residue:
			if ( i == anchor_res ) continue;
			if ( i > res_per_symm_repeat ) { //This is a symmetry repeat
				core::Size res_to_copy( i % res_per_symm_repeat );
				if ( res_to_copy == 0 ) { res_to_copy = res_per_symm_repeat; }
				if ( required_symmetry_mirroring_ && ( (i-1) / res_per_symm_repeat ) % 2 == 1 ) {
					genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::mirror_backbone_dihedrals );
				} else {
					genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::copy_backbone_dihedrals );
				}
				genkic->add_residue_to_perturber_residue_list( res_to_copy );
				genkic->add_residue_to_perturber_residue_list( i );
				if ( required_symmetry_perturbation_ != 0.0 ) {
					genkic->add_perturber( protocols::generalized_kinematic_closure::perturber::perturb_dihedral );
					if ( pose->residue_type(i).is_alpha_aa() ) {
						core::id::NamedAtomID const Natom( "N", i );
						core::id::NamedAtomID const CAatom( "CA", i );
						core::id::NamedAtomID const Catom( "C", i );
						utility::vector1< core::id::NamedAtomID > phivect; phivect.push_back( Natom ); phivect.push_back( CAatom );
						utility::vector1< core::id::NamedAtomID > psivect; psivect.push_back( CAatom ); psivect.push_back( Catom );
						genkic->add_atomset_to_perturber_atomset_list(phivect);
						genkic->add_atomset_to_perturber_atomset_list(psivect);
					} else if ( pose->residue_type(i).is_oligourea() ) {
						core::id::NamedAtomID const Natom( "N", i );
						core::id::NamedAtomID const CAatom( "CA", i );
						core::id::NamedAtomID const CMatom( "CM", i );
						core::id::NamedAtomID const NUatom( "NU", i );
						utility::vector1< core::id::NamedAtomID > phivect; phivect.push_back( Natom ); phivect.push_back( CAatom );
						utility::vector1< core::id::NamedAtomID > thetavect; thetavect.push_back( CAatom ); thetavect.push_back( CMatom );
						utility::vector1< core::id::NamedAtomID > psivect; psivect.push_back( CMatom ); psivect.push_back( NUatom );
						genkic->add_atomset_to_perturber_atomset_list(phivect);
						genkic->add_atomset_to_perturber_atomset_list(thetavect);
						genkic->add_atomset_to_perturber_atomset_list(psivect);
					}
					genkic->add_value_to_perturber_value_list( required_symmetry_perturbation_ );
				}
			}
		}
	}

	//Add bump check filter:
	genkic->add_filter( protocols::generalized_kinematic_closure::filter::loop_bump_check );

	//Add rama check filters:
	if ( use_rama_filter() ) {
		for ( core::Size i=1; i<=nres; ++i ) {
			if ( i!=first_loop_res && i!=middle_loop_res && i!=last_loop_res ) continue; //Just filter the pivots.
			if ( (cyclization_type() == SCPA_terminal_disulfide || cyclization_type() == SCPA_sidechain_isopeptide || cyclization_type() == SCPA_lanthipeptide ) &&
					( i == cyclization_point_end || i == cyclization_point_start ) ) continue;
			else if ( cyclization_type() == SCPA_nterm_isopeptide_lariat && i == cyclization_point_start ) continue;
			else if ( cyclization_type() == SCPA_cterm_isopeptide_lariat && i == cyclization_point_end ) continue;

			if ( use_rama_prepro_for_sampling() ) {
				genkic->add_filter( protocols::generalized_kinematic_closure::filter::rama_prepro_check );
				genkic->set_filter_resnum(i);
				genkic->set_filter_rama_cutoff_energy( rama_cutoff_ );
				if ( i==first_loop_res ) genkic->set_filter_attach_boinc_ghost_observer(true);
			} else {
				if ( pose->residue(i).type().is_alpha_aa() ) {
					genkic->add_filter( protocols::generalized_kinematic_closure::filter::alpha_aa_rama_check );
					genkic->set_filter_resnum(i);
					genkic->set_filter_rama_cutoff_energy( rama_cutoff_ );
					if ( i==first_loop_res ) genkic->set_filter_attach_boinc_ghost_observer(true);
				}
			}
		}
	}

	//Apply the mover:
	genkic->apply( *pose );

	return genkic->last_run_successful();
}

/// @brief Set up the TaskOperations that conrol the design process, given user inputs.
/// @details Default behaviour is designing all positions with L-canonicals and their
/// D-equivalents EXCEPT cys and met (and gly), unless the user overrides this.
void
SimpleCycpepPredictApplication::set_up_design_taskoperations(
	protocols::denovo_design::movers::FastDesignOP fdes,
	core::Size const cyclic_offset,
	core::Size const nres,
	core::pose::PoseCOP pose
) const {
	using namespace core::pack::palette;
	using namespace core::pack::task;

	CustomBaseTypePackerPaletteOP palette( utility::pointer::make_shared< CustomBaseTypePackerPalette >() );
	static const utility::fixedsizearray1< std::string, 18 > d_names( {"DALA","DASP","DGLU","DPHE","DHIS","DILE","DLYS","DLEU","DMET","DASN","DPRO","DGLN","DARG","DSER","DTHR","DVAL","DTRP","DTYR"} );
	for ( std::string const & name : d_names ) {
		palette->add_type( name );
	}

	TaskFactoryOP tf( utility::pointer::make_shared< TaskFactory >() );
	tf->set_packer_palette( palette );

	std::stringstream L_resfile("");
	std::stringstream D_resfile("");

	L_resfile << "start" << std::endl;
	D_resfile << "start" << std::endl;

	static const std::string L_default( "ADEFHIKLNPQRSTVWY" );
	static const std::string D_default( "X[DALA]X[DASP]X[DGLU]X[DPHE]X[DHIS]X[DILE]X[DLYS]X[DLEU]X[DASN]X[DPRO]X[DGLN]X[DARG]X[DSER]X[DTHR]X[DVAL]X[DTRP]X[DTYR]" );
	for ( core::Size i=1; i<=nres; ++i ) {
		core::Size orig_res( i + cyclic_offset ) ;
		if ( orig_res > nres ) orig_res -= nres;
		debug_assert( orig_res >= 1 && orig_res <= nres ); //Should be true.

		if ( allowed_canonicals_by_position_.count( static_cast<core::Size>(orig_res) ) == 1 ) {
			if ( allowed_canonicals_by_position_.at( static_cast<core::Size>(orig_res) ).size() > 0 ) {
				L_resfile << i << " A PIKAA " << get_oneletter_codes( allowed_canonicals_by_position_.at( static_cast<core::Size>(orig_res) ) ) << std::endl;
			} else {
				L_resfile << i << " A NATAA" << std::endl;
			}
		} else if ( allowed_canonicals_by_position_.count( 0 ) == 1 ) {
			if ( allowed_canonicals_by_position_ .at( 0 ).size() > 0 ) {
				L_resfile << i << " A PIKAA " << get_oneletter_codes( allowed_canonicals_by_position_.at( 0 ) ) << std::endl;
			} else {
				L_resfile << i << " A NATAA" << std::endl;
			}
		} else {
			L_resfile << i << " A PIKAA " << L_default;
			if ( !prohibit_D_at_negative_phi_ ) {
				L_resfile << D_default;
			}
			L_resfile << std::endl;
		}
		if ( allowed_noncanonicals_by_position_.count( static_cast<core::Size>(orig_res) ) == 1 ) {
			if ( allowed_noncanonicals_by_position_.at( static_cast<core::Size>(orig_res) ).size() > 0 ) {
				D_resfile << i << " A PIKAA " << get_nc_name_codes( allowed_noncanonicals_by_position_.at( static_cast<core::Size>(orig_res) ) ) << std::endl;
			}
		} else if ( allowed_noncanonicals_by_position_.count( 0 ) == 1 ) {
			if ( allowed_noncanonicals_by_position_.at( 0 ).size() > 0 ) {
				D_resfile << i << " A PIKAA " << get_nc_name_codes( allowed_noncanonicals_by_position_.at( 0 ) ) << std::endl;
			}
		} else {
			D_resfile << i << " A PIKAA " << D_default;
			if ( !prohibit_L_at_positive_phi_ ) {
				D_resfile << L_default;
			}
			D_resfile << std::endl;
		}
	}

	if ( TR.Debug.visible() ) {
		TR.Debug << "L-resfile:\n" << L_resfile.str() << std::endl;
		TR.Debug << "D-resfile:\n" << D_resfile.str() << std::endl;
	}

	core::pack::task::operation::ReadResfileOP L_resfile_taskop( utility::pointer::make_shared< core::pack::task::operation::ReadResfile >() );
	core::pack::task::operation::ReadResfileOP D_resfile_taskop( utility::pointer::make_shared< core::pack::task::operation::ReadResfile >() );

	L_resfile_taskop->set_cached_resfile( L_resfile.str() );
	D_resfile_taskop->set_cached_resfile( D_resfile.str() );


	core::select::residue_selector::PhiSelectorOP neg_phi_selector( utility::pointer::make_shared< core::select::residue_selector::PhiSelector >() );
	neg_phi_selector->set_select_positive_phi(false);
	L_resfile_taskop->set_residue_selector( neg_phi_selector );
	core::select::residue_selector::PhiSelectorOP pos_phi_selector( utility::pointer::make_shared< core::select::residue_selector::PhiSelector >() );
	pos_phi_selector->set_select_positive_phi(true);
	D_resfile_taskop->set_residue_selector( pos_phi_selector );

	tf->push_back( L_resfile_taskop );
	tf->push_back( D_resfile_taskop );

	//Prohibit design at isopeptide cyclization positions:
	if ( cyclization_type() == SCPA_nterm_isopeptide_lariat || cyclization_type() == SCPA_cterm_isopeptide_lariat || cyclization_type() == SCPA_sidechain_isopeptide ) {
		core::select::residue_selector::ResidueIndexSelectorOP linker_selector( utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >() );
		core::Size firstres, lastres;
		find_first_and_last_isopeptide_residues(pose, firstres, lastres);
		if ( cyclization_type() == SCPA_nterm_isopeptide_lariat || cyclization_type() == SCPA_sidechain_isopeptide ) {
			linker_selector->append_index( lastres );
		}
		if ( cyclization_type() == SCPA_cterm_isopeptide_lariat || cyclization_type() == SCPA_sidechain_isopeptide ) {
			linker_selector->append_index( firstres );
		}

		core::pack::task::operation::RestrictToRepackingRLTOP restrict_repacking( utility::pointer::make_shared< core::pack::task::operation::RestrictToRepackingRLT >() );
		core::pack::task::operation::OperateOnResidueSubsetOP op_on_subset( utility::pointer::make_shared< core::pack::task::operation::OperateOnResidueSubset >( restrict_repacking, linker_selector ) );
		tf->push_back(op_on_subset);
	}

	fdes->set_task_factory( tf );
}

/// @brief Given a vector of full residue names of canonical residues, give me a concatenated list of one-letter codes.
/// @details Does no checking for duplicates.
std::string
SimpleCycpepPredictApplication::get_oneletter_codes(
	utility::vector1< std::string > const &fullnames
) const {
	std::stringstream outstr;
	for ( core::Size i=1, imax=fullnames.size(); i<=imax; ++i ) {
		if ( fullnames[i] == "HIS_D" ) {
			outstr << "H";
		} else { //We can use the fact that the full name and the three-letter code are the same for canonicals:
			outstr << core::chemical::oneletter_code_from_aa( core::chemical::aa_from_name( fullnames[i] ) );
		}
	}

	return outstr.str();
}

/// @brief Given a vector of full residue names, give me a string of the form "X[<fullname1>]X[<fullname2>]X[<fullname3>] ..."
/// @details Does no checking for duplicates.  Will fail gracelessly with invalid names.
std::string
SimpleCycpepPredictApplication::get_nc_name_codes(
	utility::vector1< std::string> const &fullnames
) const {
	using namespace core::chemical;
	std::stringstream outstr;
	ResidueTypeSetCOP restype_set( ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );
	for ( core::Size i=1, imax=fullnames.size(); i<=imax; ++i ) {
		ResidueTypeCOP curtype( ResidueTypeFinder( *restype_set ).residue_base_name( fullnames[i] ).get_representative_type() );
		runtime_assert_string_msg( curtype, "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::get_nc_name_codes(): The name " + fullnames[i] + " corresponds to no known residue type." );
		outstr << "X[" << curtype->base_name() << "]";
	}
	return outstr.str();
}


/// @brief Given a pose, store a list of the disulfides in the pose.
/// @details Clears the old_disulfides list and repopulates it.
void
SimpleCycpepPredictApplication::store_disulfides (
	core::pose::PoseCOP pose,
	utility::vector1 < std::pair < core::Size, core::Size > > &old_disulfides
) const {
	old_disulfides.clear();
	core::conformation::disulfide_bonds( pose->conformation(), old_disulfides );
	return;
}

/// @brief Given a pose and a list of the disulfides in the pose, break the disulfides.
///
void
SimpleCycpepPredictApplication::break_disulfides (
	core::pose::PoseOP pose,
	utility::vector1 < std::pair < core::Size, core::Size > > const &disulfides
) const {
	for ( core::Size i=1, imax=disulfides.size(); i<=imax; ++i ) {
		core::conformation::break_disulfide( pose->conformation(), disulfides[i].first, disulfides[i].second );
	}
	return;
}

/// @brief Given a pose and a list of the disulfides that should be in the pose, form the disulfides.
///
void
SimpleCycpepPredictApplication::rebuild_disulfides (
	core::pose::PoseOP pose,
	utility::vector1 < std::pair < core::Size, core::Size > > const &disulfides
) const {
	for ( core::Size i=1, imax=disulfides.size(); i<=imax; ++i ) {
		core::conformation::form_disulfide( pose->conformation(), disulfides[i].first, disulfides[i].second, true, false );
	}
	return;
}

/// @brief Given a list of old disulfide positions, generate a list of new disulfide positions based on the offset.
/// @details Replaces the new_disulfides list.
void
SimpleCycpepPredictApplication::depermute_disulfide_list(
	utility::vector1 < std::pair < core::Size, core::Size > > const &old_disulfides,
	utility::vector1 < std::pair < core::Size, core::Size > > &new_disulfides,
	core::Size const offset,
	core::Size const nres
) const {
	//TR << "Depermuting disulfide list." << std::endl; //DELETE ME -- FOR DEBUGGING ONLY.

	new_disulfides.clear();
	for ( core::Size i=1, imax=old_disulfides.size(); i<=imax; ++i ) {
		//TR << "Old:\t" << old_disulfides[i].first << "\t" << old_disulfides[i].second << std::endl; //DELETE ME -- FOR DEBUGGING ONLY.
		core::Size res1( old_disulfides[i].first + offset );
		if ( res1 > nres ) res1 -= nres;
		core::Size res2( old_disulfides[i].second + offset );
		if ( res2 > nres ) res2 -= nres;
		new_disulfides.push_back( std::pair< core::Size, core::Size >( res1, res2 ) );
		//TR << "New:\t" << new_disulfides[i].first << "\t" << new_disulfides[i].second << std::endl; //DELETE ME -- FOR DEBUGGING ONLY.
	}

	return;
}

/// @brief Given a pose that has undergone an N-residue cyclic permutation, restore
/// the original pose, without the permutation.
void
SimpleCycpepPredictApplication::depermute (
	core::pose::PoseOP pose,
	core::Size const offset
) const {

	/// 1 2 3 4 5 6 7 8
	/// 2 3 4 5 6 7 8 1
	/// 3 4 5 6 7 8 1 2
	/// 4 5 6 7 8 1 2 3

	if ( offset==0 ) return; //Do nothing if the pose was not offset.

	core::Size const nres( sequence_length() );
	debug_assert(nres > offset);
	core::Size const old_first_res_index( nres-offset+1 );

	//TR << "nres=" << nres << " offset=" << offset << " old_first_res_index=" << old_first_res_index << std::endl; //DELETE ME

	//Store the old disulfides:
	utility::vector1 < std::pair < core::Size, core::Size > > old_disulfides;
	store_disulfides( pose, old_disulfides );
	//Break the old disulfides:
	break_disulfides( pose, old_disulfides );

	core::pose::PoseOP newpose( utility::pointer::make_shared< core::pose::Pose >() );

	for ( core::Size ir=old_first_res_index; ir<=nres; ++ir ) {
		if ( ir == old_first_res_index ) {
			newpose->append_residue_by_jump( *(pose->residue(ir).clone()), 0, "", "", true );
		} else {
			newpose->append_residue_by_bond( *(pose->residue(ir).clone()), false, 0, 0, 0, false, false );
		}
	}

	for ( core::Size ir=1; ir<old_first_res_index; ++ir ) {
		newpose->append_residue_by_bond( *(pose->residue(ir).clone()), false, 0, 0, 0, false, false );
	}

	//Depermute the old disulfide list:
	utility::vector1 < std::pair < core::Size, core::Size > > new_disulfides;
	depermute_disulfide_list( old_disulfides, new_disulfides, offset, nres );

	//Re-form the disulfides:
	rebuild_disulfides( newpose, new_disulfides );

	//Re-append paraBBMB linker residues:
	if ( parabbmb_positions_.size() > 0 ) {
		re_append_linker_residues( pose, newpose, offset, parabbmb_positions_, "BB4" );
	}

	//Re-append TBMB linker residues:
	if ( tbmb_positions_.size() > 0 ) {
		re_append_linker_residues( pose, newpose, offset, tbmb_positions_, "TBM" );
	}

	//Re-append TMA linker residues:
	if ( tma_positions_.size() > 0 ) {
		re_append_linker_residues( pose, newpose, offset, tma_positions_, "TMA" );
	}

	//Re-append trigonal planar linker residues:
	if ( trigonal_planar_metal_positions_.size() > 0 ) {
		re_append_trigonal_metal_residues( pose, newpose, true );
	}

	//Re-append trigonal pyramidal linker residues:
	if ( trigonal_pyramidal_metal_positions_.size() > 0 ) {
		re_append_trigonal_metal_residues( pose, newpose, false );
	}

	//Re-append square_pyramidal metal linker residues:
	if ( square_pyramidal_metal_positions_.size() > 0 ) {
		re_append_square_pyramidal_metal_residues( pose, newpose );
	}

	//Re-append square_planar metal linker residues:
	if ( square_planar_metal_positions_.size() > 0 ) {
		re_append_square_planar_metal_residues( pose, newpose );
	}

	//Re-append tetrahedral metal linker residues:
	if ( tetrahedral_metal_positions_.size() > 0 ) {
		re_append_tetrahedral_metal_residues( pose, newpose );
	}

	//Re-append octahedral metal linker residues:
	if ( octahedral_metal_positions_.size() > 0 ) {
		re_append_octahedral_metal_residues( pose, newpose );
	}

	//I don't bother to set up cyclic constraints, since we won't be doing any more minimization after calling this function.

	//Set up disulfide variants, if we're doing disulfide cyclization.
	if ( cyclization_type() == SCPA_terminal_disulfide ) {
		set_up_terminal_disulfide_variants( newpose );
	}

	//Mover to cyclize the polymer and to update terminal peptide bond O and H atoms:
	protocols::simple_moves::DeclareBondOP termini( utility::pointer::make_shared< protocols::simple_moves::DeclareBond >() );
	set_up_cyclization_mover( termini, newpose );
	termini->apply(*newpose);

	(*pose) = (*newpose);

	return;
}

/// @brief Create a new checkpoint file.
///
void
SimpleCycpepPredictApplication::new_checkpoint_file() const {
	using namespace utility::io;

	if ( !suppress_checkpoints_ ) {
		ozstream outfile;
		outfile.open( checkpoint_filename_ );
		runtime_assert(outfile.good());
		outfile << checkpoint_job_identifier_ << std::endl;
		outfile << "LAST\t0\tSUCCESS\t0" << std::endl;
		outfile.flush();
		outfile.close();
	}

	erase_random_seed_info();
	store_random_seed_info();

	return;
}


/// @brief Initialize checkpointing for this run.
/// @details  This function does several things.  First, it checks for an existing checkpoint
/// file.  If one exists, it checks whether the unique job name in the file matches the current
/// job.  If it does, then this job has already been attempted, and we're somewhere in the middle
/// of it.  The function reads the last attempt number and success count from the checkpoint
/// file, and returns these values.  Otherwise, it creates a new checkpoint file with the current
/// job name and returns (0,0).  If checkpointing is disabled, this function does nothing, and
/// returns (0,0).
/// @param[out] lastjob The index of the last job run.  Set to zero if checkpointing is disabled
/// or if we're creating a new checkpoint file (first job run).
/// @param[out] successes The number of successes so far.  Set to zero if checkpointing is
/// disabled or if we're creating a new checkpoint file (first job run).
void
SimpleCycpepPredictApplication::initialize_checkpointing(
	core::Size &lastjob,
	core::Size &successes
) const {
	using namespace utility::io;

	//If we're not using checkpointing, return 0,0:
	if ( checkpoint_job_identifier_ == "" || suppress_checkpoints_ ) {
		lastjob=0;
		successes=0;
		return;
	}

	//Check for a checkpoint file:
	izstream infile;
	infile.open( checkpoint_filename_ );
	if ( !infile.good() ) {
		//If the checkpoint file doesn't exist/isn't readable, then we're starting a new job and need a new checkpoint file.
		infile.close();
		new_checkpoint_file();
		lastjob=0;
		successes=0;
		return;
	}

	//If we've reached this point, then the checkpoint file IS good, and we need to read the job name and the last job ID/success count:
	std::string curline;
	infile.getline(curline);
	if ( infile.eof() || curline=="" || curline!=checkpoint_job_identifier_ ) {
		//If the checkpoint file isn't readable or is for a different job, then we're starting a new job and need a new checkpoint file.
		infile.close();
		new_checkpoint_file();
		lastjob=0;
		successes=0;
		return;
	}

	//Loop through the checkpoint file and read the job lines:
	while ( !infile.eof() ) {
		infile.getline(curline);
		std::istringstream ss(curline);
		ss >> curline;
		ss >> lastjob;
		ss >> curline;
		ss >> successes;
	}
	infile.close();

	get_random_seed_info();

	if ( TR.Debug.visible() ) {
		TR.Debug << "Initialized job to " << lastjob << ", successes to " << successes << "." << std::endl;
		TR.Debug.flush();
	}

	return;
}

/// @brief Add a checkpoint to the checkpoint file.
/// @details  The checkpoint file must already exist.  Does nothing if checkpointing is disabled.
/// @param[in] curjob The index of the current job just run, for writing to the checkpoint file.
/// @param[in] successes The number of successes so far, for writing to the checkpoint file.
void
SimpleCycpepPredictApplication::checkpoint(
	core::Size const curjob,
	core::Size const successes
) const {
	using namespace utility::io;

	//Do nothing if we're not using checkpointing:
	if ( checkpoint_job_identifier_ == "" || suppress_checkpoints_ ) {
		return;
	}

	ozstream outfile;
	outfile.open_append( checkpoint_filename_ );
	runtime_assert(outfile.good());
	outfile << "LAST\t" << curjob << "\tSUCCESS\t" << successes << std::endl;
	outfile.flush();
	outfile.close();

	store_random_seed_info();

#ifdef BOINC
	protocols::boinc::Boinc::update_pct_complete();
#endif

	return;

}

/// @brief End checkpointing and delete the checkpoint file.
/// @details Does nothing if checkpointing is disabled.
void
SimpleCycpepPredictApplication::end_checkpointing() const {
	using namespace utility::io;

	//Do nothing if we're not using checkpointing:
	if ( checkpoint_job_identifier_ == "" || suppress_checkpoints_ ) {
		return;
	}
	runtime_assert( remove( checkpoint_filename_.c_str() ) == 0 );
	erase_random_seed_info();
	return;
}

/// @brief Restore the state of the random generator from a previous run.
///
void SimpleCycpepPredictApplication::get_random_seed_info() const {
#ifdef BOINC
	boinc_begin_critical_section();
#endif
	if ( utility::file::file_exists(rand_checkpoint_file_) ) {
		utility::io::izstream izs(rand_checkpoint_file_);
		numeric::random::rg().restoreState(izs);
		izs.close();
	}
#ifdef BOINC
	boinc_end_critical_section();
#endif
	return;
}

/// @brief Store the state of the random generator from a previous run.
///
void SimpleCycpepPredictApplication::store_random_seed_info() const {
	if ( suppress_checkpoints_ ) return; //Do nothing if we're not checkpointing
#ifdef BOINC
	boinc_begin_critical_section();
#endif
	utility::io::ozstream ozs(rand_checkpoint_file_);
	numeric::random::rg().saveState(ozs);
	ozs.close();
#ifdef BOINC
	boinc_end_critical_section();
#endif
	return;
}

/// @brief Erase the stored state of the random generator from a previous run.
///
void SimpleCycpepPredictApplication::erase_random_seed_info() const {
	if ( suppress_checkpoints_ ) return; //Do nothing if we're not checkpointing
	if ( utility::file::file_exists(rand_checkpoint_file_) ) {
		utility::file::file_delete(rand_checkpoint_file_);
	}
	return;
}

/// @brief Given a pose with a linker (e.g. TBMB, paraBBMB, TMA) in it and another pose without the linker, copy the linker residues from the first to the second,
/// and add back covalent bonds.
/// @details This function is called at the end of the protocol, and therefore doesn't bother to add back constraints.
void
SimpleCycpepPredictApplication::re_append_linker_residues(
	core::pose::PoseCOP pose,
	core::pose::PoseOP newpose,
	core::Size const offset,
	utility::vector1< utility::vector1< core::Size > > const &linker_positions,
	std::string const & linker_name
) const {
	debug_assert( pose->total_residue() > newpose->total_residue() );

	core::Size lastres(0);
	for ( core::Size i=1, imax=linker_positions.size(); i<=imax; ++i ) { //Loop through all linkers
		//For each one, loop through the sequence and find the next linker of the given type.
		for ( core::Size j=lastres+1, jmax=pose->total_residue(); j<=jmax; ++j ) {
			if ( pose->residue_type(j).name3() == linker_name ) {
				lastres = j;
				break;
			}
		}
		newpose->append_residue_by_jump( pose->residue(lastres), linker_positions[i][1] ); //Jump from the first residue that links this linker, to the linker itself.
		protocols::cyclic_peptide::crosslinker::CrosslinkerMoverHelperOP helper;
		bool twofold_linker(false);
		if ( linker_name == "TBM" ) {
			helper = utility::pointer::make_shared< protocols::cyclic_peptide::crosslinker::TBMB_Helper >();
		} else if ( linker_name == "TMA" ) {
			helper = utility::pointer::make_shared< protocols::cyclic_peptide::crosslinker::TMA_Helper >();
		} else if ( linker_name == "BB4" ) {
			helper = utility::pointer::make_shared< protocols::cyclic_peptide::crosslinker::One_Four_BBMB_Helper >();
			twofold_linker=true;
		}

		utility::vector1< core::Size > res_indices( twofold_linker ? 2 : 3);
		if ( twofold_linker ) {
			if ( offset >= linker_positions[i][2] || offset < linker_positions[i][1] ) {
				res_indices[1] = linker_positions[i][1]; res_indices[2] = linker_positions[i][2];
			} else {
				res_indices[1] = linker_positions[i][2]; res_indices[2] = linker_positions[i][1];
			}
		} else {
			if ( offset >= linker_positions[i][3] || offset < linker_positions[i][1] ) {
				res_indices[1] = linker_positions[i][1]; res_indices[2] = linker_positions[i][2]; res_indices[3] = linker_positions[i][3];
			} else if ( offset >= linker_positions[i][2] ) {
				res_indices[1] = linker_positions[i][3]; res_indices[2] = linker_positions[i][1]; res_indices[3] = linker_positions[i][2];
			} else {
				res_indices[1] = linker_positions[i][2]; res_indices[2] = linker_positions[i][3]; res_indices[3] = linker_positions[i][1];
			}
		}
		helper->add_linker_bonds_asymmetric( *newpose, res_indices, newpose->total_residue() );
	}
}

/// @brief Given a pose with trigonal pyramidal or trigonal planar metal variants and another pose without the variants, add back the variants to the latter.
/// @details This function is called at the end of the protocol, and therefore doesn't bother to add back contraints.  If planar is true, it calls the trigonal planar code; otherwise it calls the trigonal pyramidal.
void
SimpleCycpepPredictApplication::re_append_trigonal_metal_residues(
	core::pose::PoseCOP pose,
	core::pose::PoseOP newpose,
	bool const planar
) const {
	utility::vector1< std::pair< utility::fixedsizearray1< core::Size, 3 >, std::string > > const & positions( planar ? trigonal_planar_metal_positions_ : trigonal_pyramidal_metal_positions_ );
	debug_assert( pose->total_residue() == newpose->total_residue() );
	for ( core::Size i(1), imax(positions.size()); i<=imax; ++i ) {
		utility::fixedsizearray1< core::Size, 3 > const &residues( positions[i].first );
		debug_assert(residues.size() == 3);
		core::select::residue_selector::ResidueSubset selection( newpose->total_residue() );
		for ( core::Size i(1); i<=3; ++i ) {
			selection[residues[i]] = true;
		}
		if ( planar ) {
			protocols::cyclic_peptide::crosslinker::TrigonalPlanarMetal_Helper helper;
			helper.add_linker_asymmetric( *newpose, selection );
		} else {
			protocols::cyclic_peptide::crosslinker::TrigonalPyramidalMetal_Helper helper;
			helper.add_linker_asymmetric( *newpose, selection );
		}
	}
}

/// @brief Given a pose with square pyramidal metal variants and another pose without the variants, add back the variants to the latter.
/// @details This function is called at the end of the protocol, and therefore doesn't bother to add back contraints.
void
SimpleCycpepPredictApplication::re_append_square_pyramidal_metal_residues(
	core::pose::PoseCOP pose,
	core::pose::PoseOP newpose
) const {
	debug_assert( pose->total_residue() == newpose->total_residue() );
	for ( core::Size i(1), imax(square_pyramidal_metal_positions_.size()); i<=imax; ++i ) {
		utility::fixedsizearray1< core::Size, 5 > const &residues( square_pyramidal_metal_positions_[i].first );
		debug_assert(residues.size() == 5);
		protocols::cyclic_peptide::crosslinker::SquarePyramidalMetal_Helper helper;
		core::select::residue_selector::ResidueSubset selection( newpose->total_residue() );
		for ( core::Size i(1); i<=5; ++i ) {
			selection[residues[i]] = true;
		}
		helper.add_linker_asymmetric( *newpose, selection );
	}
}

/// @brief Given a pose with square planar metal variants and another pose without the variants, add back the variants to the latter.
/// @details This function is called at the end of the protocol, and therefore doesn't bother to add back contraints.
void
SimpleCycpepPredictApplication::re_append_square_planar_metal_residues(
	core::pose::PoseCOP pose,
	core::pose::PoseOP newpose
) const {
	debug_assert( pose->total_residue() == newpose->total_residue() );
	for ( core::Size i(1), imax(square_planar_metal_positions_.size()); i<=imax; ++i ) {
		utility::fixedsizearray1< core::Size, 4 > const &residues( square_planar_metal_positions_[i].first );
		debug_assert(residues.size() == 4);
		protocols::cyclic_peptide::crosslinker::SquarePlanarMetal_Helper helper;
		core::select::residue_selector::ResidueSubset selection( newpose->total_residue() );
		for ( core::Size i(1); i<=4; ++i ) {
			selection[residues[i]] = true;
		}
		helper.add_linker_asymmetric( *newpose, selection );
	}
}

/// @brief Given a pose with tetrahedral metal variants and another pose without the variants, add back the variants to the latter.
/// @details This function is called at the end of the protocol, and therefore doesn't bother to add back contraints.
void
SimpleCycpepPredictApplication::re_append_tetrahedral_metal_residues(
	core::pose::PoseCOP pose,
	core::pose::PoseOP newpose
) const {
	debug_assert( pose->total_residue() == newpose->total_residue() );
	for ( core::Size i(1), imax(tetrahedral_metal_positions_.size()); i<=imax; ++i ) {
		utility::fixedsizearray1< core::Size, 4 > const &residues( tetrahedral_metal_positions_[i].first );
		debug_assert(residues.size() == 4);
		protocols::cyclic_peptide::crosslinker::TetrahedralMetal_Helper helper;
		core::select::residue_selector::ResidueSubset selection( newpose->total_residue() );
		for ( core::Size i(1); i<=4; ++i ) {
			selection[residues[i]] = true;
		}
		helper.add_linker_asymmetric( *newpose, selection );
	}
}

/// @brief Given a pose with octahedral metal variants and another pose without the variants, add back the variants to the latter.
/// @details This function is called at the end of the protocol, and therefore doesn't bother to add back contraints.
void
SimpleCycpepPredictApplication::re_append_octahedral_metal_residues(
	core::pose::PoseCOP pose,
	core::pose::PoseOP newpose
) const {
	debug_assert( pose->total_residue() == newpose->total_residue() );
	for ( core::Size i(1), imax(octahedral_metal_positions_.size()); i<=imax; ++i ) {
		utility::fixedsizearray1< core::Size, 6 > const &residues( octahedral_metal_positions_[i].first );
		protocols::cyclic_peptide::crosslinker::OctahedralMetal_Helper helper;
		core::select::residue_selector::ResidueSubset selection( newpose->total_residue() );
		for ( core::Size i(1); i<=6; ++i ) {
			selection[residues[i]] = true;
		}
		helper.add_linker_asymmetric( *newpose, selection );
	}
}

/// @brief Given a pose, return the index of the first residue that can form a disulfide.
/// @details Throws an error if no residue is found.
core::Size
SimpleCycpepPredictApplication::find_first_disulf_res(
	core::pose::PoseCOP pose
) const {
	for ( core::Size i(1), imax(pose->total_residue()); i<=imax; ++i ) {
		if ( pose->residue_type(i).is_sidechain_thiol() || pose->residue_type(i).is_disulfide_bonded() ) return i;
	}

	utility_exit_with_message( "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::find_first_disulf_res(): No disulfide-forming residue could be found!" );
	return 0;
}

/// @brief Given a pose, return the index of the last residue that can form a disulfide.
/// @details Throws an error if no residue is found.
core::Size
SimpleCycpepPredictApplication::find_last_disulf_res(
	core::pose::PoseCOP pose
) const {
	for ( core::Size i(pose->total_residue()); i>0; --i ) {
		if ( pose->residue_type(i).is_sidechain_thiol() || pose->residue_type(i).is_disulfide_bonded() ) return i;
	}

	utility_exit_with_message( "Error in protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::find_last_disulf_res(): No disulfide-forming residue could be found!" );
	return 0;
}

/// @brief Find the first and last polymer residues in a pose.
void
SimpleCycpepPredictApplication::find_first_and_last_polymer_residues(
	core::pose::Pose const & pose,
	core::Size & first_polymer_res,
	core::Size & last_polymer_res
) const {
	std::string const errmsg( "Error in SimpleCycpepPredictApplication::find_first_and_last_polymer_residues(): " );

	for ( core::Size ir(1), irmax(pose.total_residue()); ir<=irmax; ++ir ) {
		if ( pose.residue_type(ir).is_polymer() ) {
			last_polymer_res = ir;
			if ( first_polymer_res == 0 ) first_polymer_res = ir;
		}
	}
	runtime_assert_string_msg( last_polymer_res != 0, errmsg + "No polymeric residue types were found in this pose!  This isn't a macrocycle!" );
	runtime_assert_string_msg( last_polymer_res > first_polymer_res, errmsg + "Only one polymeric residue type was found in this pose!  This isn't a macrocycle!" );
}

/// @brief Given a pose, add disulfide variant types to the first and last cysteine residues in the pose.
/// @details This should ONLY be called on a pose just before a bond is declared between these residues.
void
SimpleCycpepPredictApplication::set_up_terminal_disulfide_variants(
	core::pose::PoseOP pose
) const {

	std::string const errmsg( "Error in SimpleCycpepPredictApplication::set_up_terminal_disulfide_variants():  ");

	core::Size first_polymer_res(0), last_polymer_res(0);
	find_first_and_last_polymer_residues( *pose, first_polymer_res, last_polymer_res );

	protocols::simple_moves::ModifyVariantTypeMover add_disulf_var;
	add_disulf_var.set_additional_type_to_add("DISULFIDE");
	core::select::residue_selector::ResidueIndexSelectorOP selector( utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >() );
	core::Size const first_disulf( find_first_disulf_res(pose) ), last_disulf( find_last_disulf_res(pose) );
	runtime_assert_string_msg( first_disulf != 0 && last_disulf > first_disulf, errmsg + "Indices " + std::to_string(first_disulf) + " and " + std::to_string(last_disulf) + " are invalid for a disulfide bonded pair." );
	selector->append_index( first_disulf );
	selector->append_index( last_disulf );
	add_disulf_var.set_residue_selector(selector);
	add_disulf_var.apply(*pose);
	core::pose::add_lower_terminus_type_to_pose_residue( *pose, first_polymer_res );
	core::pose::add_upper_terminus_type_to_pose_residue( *pose, last_polymer_res );
	TR << "Added DISULFIDE variant to residues " << first_disulf << " and " << last_disulf << "." << std::endl;
	TR << "Added lower and upper terminal variant types to residues " << first_polymer_res << " and " << last_polymer_res << ", respectively." << std::endl;
}

/// @brief Given a pose, add sidechain conjugation variant types to sidechains involved in making an isopeptide
/// bond, and strip termini from termini involved in the isopeptide bond.
void
SimpleCycpepPredictApplication::set_up_isopeptide_variants(
	core::pose::PoseOP pose
) const {
	runtime_assert( cyclization_type() == SCPA_cterm_isopeptide_lariat || cyclization_type() == SCPA_nterm_isopeptide_lariat || cyclization_type() == SCPA_sidechain_isopeptide );

	core::Size firstres, lastres;
	find_first_and_last_isopeptide_residues(pose, firstres, lastres);

	core::Size first_polymer_res(0), last_polymer_res(0);
	find_first_and_last_polymer_residues( *pose, first_polymer_res, last_polymer_res );

	if ( cyclization_type() == SCPA_cterm_isopeptide_lariat ) {
		core::pose::add_variant_type_to_pose_residue(*pose, core::chemical::SIDECHAIN_CONJUGATION, firstres);
		core::pose::remove_upper_terminus_type_from_pose_residue( *pose, lastres );
		core::pose::add_lower_terminus_type_to_pose_residue( *pose, 1 );
	} else if ( cyclization_type() == SCPA_nterm_isopeptide_lariat ) {
		core::pose::add_variant_type_to_pose_residue(*pose, core::chemical::SIDECHAIN_CONJUGATION, lastres);
		core::pose::remove_lower_terminus_type_from_pose_residue( *pose, firstres );
		core::pose::add_upper_terminus_type_to_pose_residue( *pose, last_polymer_res );
	} else if ( cyclization_type() == SCPA_sidechain_isopeptide ) {
		core::pose::add_variant_type_to_pose_residue(*pose, core::chemical::SIDECHAIN_CONJUGATION, firstres);
		core::pose::add_variant_type_to_pose_residue(*pose, core::chemical::SIDECHAIN_CONJUGATION, lastres);
		core::pose::add_lower_terminus_type_to_pose_residue( *pose, 1 );
		core::pose::add_upper_terminus_type_to_pose_residue( *pose, last_polymer_res );
	}
	pose->update_residue_neighbors();
}

/// @brief Given a pose, add sidechain conjugation variant types to the C-terminal
/// cysteine and add a special acetyl terminus to the N-terminal residue.
/// @returns The thioether cysteine index.
core::Size
SimpleCycpepPredictApplication::set_up_terminal_thioether_lariat_variants(
	core::pose::PoseOP pose
) const {
	runtime_assert( cyclization_type() == SCPA_thioether_lariat );

	core::Size firstres, lastres;
	find_first_and_last_thioether_lariat_residues(pose, firstres, lastres);

	core::Size first_polymer_res(0), last_polymer_res(0);
	find_first_and_last_polymer_residues( *pose, first_polymer_res, last_polymer_res );

	core::pose::add_upper_terminus_type_to_pose_residue( *pose, last_polymer_res );
	protocols::cyclic_peptide::crosslinker::set_up_thioether_variants( *pose, firstres, lastres );

	return lastres;
}

/// @brief Given the basename of a residue type, return true if this is a type that can donate the nitrogen to an
/// isopeptide bond, false otherwise.
bool SimpleCycpepPredictApplication::is_isopeptide_forming_amide_type( std::string const &basename ) const {
	return (basename == "LYS" || basename == "ORN" || basename == "DAB" || basename == "DPP" ||
		basename == "DLYS" || basename == "DORN" || basename == "DDAB" || basename == "DDPP");
}

/// @brief Given the AA of a residue type, return true if this is a type that can donate the carbonyl to an
/// isopeptide bond, false otherwise.
bool SimpleCycpepPredictApplication::is_isopeptide_forming_carbonyl_type( core::chemical::AA const aa ) const {
	using namespace core::chemical;
	return (aa == aa_asp || aa == aa_glu || aa == aa_das || aa == aa_dgu);
}

/// @brief Given a GenKIC object, a pose, and a bond length perturbation magnitude, add
/// bond length perturbation to all bond lengths in the pose.
void
SimpleCycpepPredictApplication::add_bondlength_perturbation(
	protocols::generalized_kinematic_closure::GeneralizedKIC & genkic,
	core::Real const bondlength_perturbation_magnitude,
	core::pose::Pose const & pose,
	core::Size const anchor_res
) const {
	if ( !bondlength_perturbation_magnitude ) return;
	debug_assert( cyclization_type_ == SCPA_n_to_c_amide_bond );

	using namespace protocols::generalized_kinematic_closure;
	using namespace protocols::generalized_kinematic_closure::perturber;
	using namespace core::id;

	genkic.add_perturber( perturb_bondlength );
	for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
		if ( i == anchor_res ) continue;
		core::chemical::ResidueType const & restype( pose.residue_type(i) );
		core::Size next_i( i + 1 );
		if ( next_i == anchor_res ) next_i = 0;
		if ( next_i > imax ) next_i = 1;
		if ( next_i == anchor_res ) next_i = 0; //Have to do this twice.
		core::chemical::AtomIndices const & mainchain_atoms( restype.mainchain_atoms() );
		core::Size next_atom_index( 0 );
		if ( next_i != 0 ) { //If we're not adjacent to the anchor residue.
			core::chemical::AtomIndices const & mainchain_atoms2( pose.residue_type(next_i).mainchain_atoms() );
			runtime_assert( mainchain_atoms2.size() >= 1);
			next_atom_index = mainchain_atoms2[1];
		}
		for ( core::Size j(1), jmax(mainchain_atoms.size()); j<=jmax; ++j ) {
			if ( j<jmax ) {
				genkic.add_atomset_to_perturber_atomset_list(
					utility::vector1< NamedAtomID > {
					NamedAtomID( restype.atom_name(mainchain_atoms[j]), i),
					NamedAtomID( restype.atom_name(mainchain_atoms[j+1]), i)
					}
				);
			} else { //j == jmax
				if ( next_i != 0 ) { //If we're not adjacent to the anchor residue.
					genkic.add_atomset_to_perturber_atomset_list(
						utility::vector1< NamedAtomID > {
						NamedAtomID( restype.atom_name(mainchain_atoms[j]), i),
						NamedAtomID( restype.atom_name(next_atom_index), next_i)
						}
					);
				}
			}
		}
	}
	genkic.add_value_to_perturber_value_list( bondlength_perturbation_magnitude );
}

/// @brief Given a GenKIC object, a pose, and a bond angle perturbation magnitude, add
/// bond angle perturbation to all bond angles in the pose.
void
SimpleCycpepPredictApplication::add_bondangle_perturbation(
	protocols::generalized_kinematic_closure::GeneralizedKIC & genkic,
	core::Real const bondangle_perturbation_magnitude,
	core::pose::Pose const & pose,
	core::Size const anchor_res
) const {
	if ( !bondangle_perturbation_magnitude ) return;
	debug_assert( cyclization_type_ == SCPA_n_to_c_amide_bond );

	using namespace protocols::generalized_kinematic_closure;
	using namespace protocols::generalized_kinematic_closure::perturber;
	using namespace core::id;

	genkic.add_perturber( perturb_bondangle );
	for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
		if ( i == anchor_res ) continue;
		core::chemical::ResidueType const & restype( pose.residue_type(i) );
		core::Size next_i( i + 1 ), prev_i( i - 1 );
		if ( next_i == anchor_res ) next_i = 0;
		if ( next_i > imax ) next_i = 1;
		if ( next_i == anchor_res ) next_i = 0; //Have to do this twice.
		if ( i == 1 ) prev_i = imax;
		if ( prev_i == anchor_res ) prev_i = 0;
		core::chemical::AtomIndices const & mainchain_atoms( restype.mainchain_atoms() );
		core::Size next_atom_index( 0 ), prev_atom_index( 0 );
		if ( next_i != 0 ) { //If we're not before the anchor residue.
			core::chemical::AtomIndices const & mainchain_atoms2( pose.residue_type(next_i).mainchain_atoms() );
			runtime_assert( mainchain_atoms2.size() >= 1);
			next_atom_index = mainchain_atoms2[1];
		}
		if ( prev_i != 0 ) { //If we're not after the anchor residue.
			core::chemical::AtomIndices const & mainchain_atoms3( pose.residue_type(prev_i).mainchain_atoms() );
			runtime_assert( mainchain_atoms3.size() >= 1);
			prev_atom_index = mainchain_atoms3[mainchain_atoms3.size()];
		}
		runtime_assert( mainchain_atoms.size() > 1 ); //Doesn't work if we have only 1 mainchain atom.
		for ( core::Size j(1), jmax(mainchain_atoms.size()); j<=jmax; ++j ) {
			if ( j == 1 ) {
				if ( prev_i != 0 ) {
					genkic.add_atomset_to_perturber_atomset_list(
						utility::vector1< NamedAtomID > {
						NamedAtomID( restype.atom_name(prev_atom_index), prev_i),
						NamedAtomID( restype.atom_name(mainchain_atoms[j]), i),
						NamedAtomID( restype.atom_name(mainchain_atoms[j+1]), i)
						}
					);
				}
			} else if ( j<jmax ) { //j > 1 but j < jmax
				genkic.add_atomset_to_perturber_atomset_list(
					utility::vector1< NamedAtomID > {
					NamedAtomID( restype.atom_name(mainchain_atoms[j-1]), i),
					NamedAtomID( restype.atom_name(mainchain_atoms[j]), i),
					NamedAtomID( restype.atom_name(mainchain_atoms[j+1]), i)
					}
				);
			} else { //j == jmax
				if ( next_i != 0 ) { //If we're not adjacent to the anchor residue.
					genkic.add_atomset_to_perturber_atomset_list(
						utility::vector1< NamedAtomID > {
						NamedAtomID( restype.atom_name(mainchain_atoms[j-1]), i),
						NamedAtomID( restype.atom_name(mainchain_atoms[j]), i),
						NamedAtomID( restype.atom_name(next_atom_index), next_i)
						}
					);
				}
			}
		}
	}
	genkic.add_value_to_perturber_value_list( bondangle_perturbation_magnitude );
}


} //cyclic_peptide_predict
} //protocols
