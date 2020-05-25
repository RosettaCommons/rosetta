// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/rna/denovo/movers/RNA_DeNovoProtocolMover.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu
/// @author Andy Watkins, amw579@stanford.edu


#include <protocols/rna/denovo/movers/RNA_DeNovoProtocolMover.hh>
#include <protocols/rna/denovo/movers/RNA_DeNovoProtocolMoverCreator.hh>
#include <core/import_pose/RNA_DeNovoParameters.hh>
#include <core/import_pose/options/RNA_DeNovoProtocolOptions.hh>
#include <core/import_pose/RNA_HelixAssembler.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/util.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/rna/chemical_shift/RNA_ChemicalShiftPotential.hh>
#include <core/scoring/func/FadeFunc.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/io/rna/RNA_DataReader.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/full_model_info/FullModelParameters.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/util.hh>
#include <core/pose/subpose_manipulation_util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/pose/rna/leontis_westhof_util.hh>
#include <core/pose/rna/RNA_SecStruct.hh>

#include <protocols/rna/denovo/RNA_FragmentMonteCarlo.hh>
// #include <core/fragment/rna/FullAtomRNA_Fragments.hh>
#include <protocols/rna/movers/RNA_LoopCloser.hh>
#include <protocols/rna/denovo/movers/RNA_Minimizer.hh>
#include <core/import_pose/RNA_BasePairHandler.hh>
#include <protocols/rna/denovo/RNA_DeNovoPoseInitializer.hh>
// #include <protocols/rna/denovo/movers/RNA_Relaxer.hh>
#include <core/import_pose/libraries/RNA_ChunkLibrary.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/pose/rna/RNA_BasePairClassifier.hh>
#include <protocols/rna/denovo/util.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_VDW_BinChecker.hh>
#include <protocols/scoring/VDW_CachedRepScreenInfo.hh>
#include <protocols/jd2/util.hh>


#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>


#include <utility/options/OptionCollection.hh>
#include <utility/options/keys/OptionKeyList.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

#include <ObjexxFCL/format.hh>

#include <basic/Tracer.hh>

#include <utility/io/izstream.hh>

static basic::Tracer TR( "protocols.rna.denovo.movers.RNA_DeNovoProtocolMover" );
using utility::vector1;
using utility::tools::make_vector1;
using namespace core;
using namespace core::chemical::rna;

/////////////////////////////////////////////////////////////////////////////////////////////////////
// @details
//
// Sets up a pose for FARFAR2 given inputs on command line or (ideally??) RosettaScripts tag
//
// Maybe we're going to have to cache the rna_de_novo_parameters in the pose.
//
//       -- Andy Watkins, May 2020
//
/////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace rna {
namespace denovo {
namespace movers {

using namespace core::pose::rna;
using namespace core::import_pose;

//Constructor
RNA_DeNovoProtocolMover::RNA_DeNovoProtocolMover()
{}

//Destructor
RNA_DeNovoProtocolMover::~RNA_DeNovoProtocolMover() = default;


std::string
RNA_DeNovoProtocolMover::get_name() const {
	return mover_name();
}

std::string
RNA_DeNovoProtocolMover::mover_name() {
	return "RNA_DeNovoProtocol";
}

moves::MoverOP RNA_DeNovoProtocolMover::clone() const {
	return utility::pointer::make_shared< RNA_DeNovoProtocolMover >( *this );
}

moves::MoverOP RNA_DeNovoProtocolMover::fresh_instance() const {
	return utility::pointer::make_shared< RNA_DeNovoProtocolMover >();
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocolMover::initialize_scorefxn( core::pose::Pose & pose ) {

	using namespace core::scoring;
	using namespace basic::options;

	// RNA low-resolution score function.
	denovo_scorefxn_ = ScoreFunctionFactory::create_score_function( options_->lores_scorefxn() );
	apply_set_weights(  denovo_scorefxn_, option[ OptionKeys::rna::denovo::set_lores_weights ]() );

	if ( core::scoring::rna::nonconst_rna_scoring_info_from_pose( pose ).rna_data_info().rna_reactivities().size() > 0 ) {
		denovo_scorefxn_->set_weight( core::scoring::rna_chem_map_lores, option[ OptionKeys::score::rna::rna_chem_map_lores_weight ]() );
	}

	if ( options_->filter_vdw() ) {
		denovo_scorefxn_->set_weight( core::scoring::grid_vdw, options_->grid_vdw_weight() );
	}


	// initial_denovo_scorefxn_ = denovo_scorefxn_->clone();
	if ( options_->chainbreak_weight() > -1.0 ) denovo_scorefxn_->set_weight( chainbreak, options_->chainbreak_weight() );
	if ( options_->linear_chainbreak_weight() > -1.0 ) denovo_scorefxn_->set_weight( linear_chainbreak, options_->linear_chainbreak_weight() );


	// inlined initialize_constraints from RNA_DeNovoProtocol.cc
	using namespace core::scoring;
	if ( pose.constraint_set()->has_constraints() ) {
		denovo_scorefxn_->set_weight( atom_pair_constraint, 1.0 );
		denovo_scorefxn_->set_weight( coordinate_constraint, 1.0 ); // now useable in RNA denovo!
		denovo_scorefxn_->set_weight( base_pair_constraint, 1.0 );
	}
	if ( options_->rmsd_screen() > 0.0 && !denovo_scorefxn_->has_nonzero_weight( coordinate_constraint) ) {
		denovo_scorefxn_->set_weight( coordinate_constraint, 1.0 ); // now useable in RNA denovo!
	}



	// RNA high-resolution score function.
	hires_scorefxn_ = get_rna_hires_scorefxn( rna_params_->is_rna_and_protein() );//->clone();

	if ( options_->model_with_density() ) {
		// Then we should turn on the density score terms in the low-res and the high-res score functions
		if ( denovo_scorefxn_->get_weight( core::scoring::elec_dens_fast ) <= 0 ) {
			denovo_scorefxn_->set_weight( core::scoring::elec_dens_fast, 10.0 );
		}
		if ( hires_scorefxn_->get_weight( core::scoring::elec_dens_fast ) <= 0 ) {
			hires_scorefxn_->set_weight( core::scoring::elec_dens_fast, 10.0 );
		}
	}
}


void RNA_DeNovoProtocolMover::apply( core::pose::Pose & pose ) {
	//

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::import_pose::options;

	// This is actually all of apply, basically.
	pose = *pose_;






	using namespace core::pose;
	using namespace core::scoring;
	using namespace core::io::pdb;
	using namespace core::io::silent;
	using namespace protocols::rna::denovo;

	// bps_moves and new_fold_tree_initializer are mutually exclusive
	if ( options_->bps_moves() && options_->new_fold_tree_initializer() ) {
		utility_exit_with_message( "If you want to supply -new_fold_tree_initializer true, as for an electron density case, you MUST supply -bps_moves false." );
	}

	///////////////////////////////////////////////////////////////////////////
	// A bunch of initialization
	///////////////////////////////////////////////////////////////////////////
	if ( options_->dump_pdb() ) pose.dump_pdb( "init.pdb" );

	// Set up the cached vdw rep screen info in the pose
	// requires further setup later (once user input fragments have been inserted in the pose)
	if ( options_->filter_vdw() ) {
		protocols::scoring::fill_vdw_cached_rep_screen_info_from_command_line( pose );
		vdw_grid_ = utility::pointer::make_shared< protocols::stepwise::modeler::rna::checker::RNA_VDW_BinChecker >( pose );
	}

	// RNA score function (both low-res and high-res).
	initialize_scorefxn( pose );

	auto rna_de_novo_pose_initializer = utility::pointer::make_shared< protocols::rna::denovo::RNA_DeNovoPoseInitializer >( *rna_params_ );
	rna_de_novo_pose_initializer->set_bps_moves( options_->bps_moves() );
	rna_de_novo_pose_initializer->set_root_at_first_rigid_body( options_->root_at_first_rigid_body() );
	rna_de_novo_pose_initializer->set_dock_each_chunk( options_->dock_each_chunk() );
	rna_de_novo_pose_initializer->set_dock_each_chunk_per_chain( options_->dock_each_chunk_per_chain() );
	rna_de_novo_pose_initializer->set_dock_chunks_res( options_->dock_chunks_res() );
	rna_de_novo_pose_initializer->set_center_jumps_in_single_stranded( options_->center_jumps_in_single_stranded() );
	rna_de_novo_pose_initializer->set_new_fold_tree_initializer( options_->new_fold_tree_initializer() );
	rna_de_novo_pose_initializer->set_model_with_density( options_->model_with_density() );
	rna_de_novo_pose_initializer->initialize_for_de_novo_protocol( pose, options_->ignore_secstruct() ); // virtualize phosphates, but no chainbreaks -- PUT HIGHER?

	//Keep a copy for resetting after each decoy.
	// Pose start_pose = pose;
	// This copy of the pose does not have its grid_vdw stuff set up
	// But it will be set up in RNA_FragmentMonteCarlo before any sampling happens

	///////////////////////////////////////////////////////////////////////////
	// Main Loop.
	///////////////////////////////////////////////////////////////////////////
	// core::Size refine_pose_id( 1 );

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace basic::options::OptionKeys::stepwise::monte_carlo;

	using namespace core::import_pose::libraries;

	RNA_ChunkLibraryOP user_input_chunk_library( new RNA_ChunkLibrary( options_->chunk_pdb_files(), options_->chunk_silent_files(), pose,
		options_->input_res(), rna_params_->allow_insert_res() ) );
	RNA_BasePairHandlerOP rna_base_pair_handler( /*refine_pose ? new RNA_BasePairHandler( pose ) :*/ new RNA_BasePairHandler( *rna_params_ ) );

	// main loop
	rna_fragment_monte_carlo_ = utility::pointer::make_shared< RNA_FragmentMonteCarlo >( options_ );

	rna_fragment_monte_carlo_->set_user_input_chunk_library( user_input_chunk_library );

	if ( options_->initial_structures_provided() ) {
		utility::vector1< std::string > silent_files_empty;
		RNA_ChunkLibraryOP user_input_chunk_initialization_library( new RNA_ChunkLibrary( options_->chunk_initialization_pdb_files(),
			silent_files_empty, pose, options_->input_res_initialize(), rna_params_->allow_insert_res() ) );
		rna_fragment_monte_carlo_->set_user_input_chunk_initialization_library( user_input_chunk_initialization_library );
	}
	rna_fragment_monte_carlo_->set_rna_base_pair_handler( rna_base_pair_handler ); // could later have this look inside pose's sec_struct_info

	//the jd will take care of outfile tagging
	//rna_fragment_monte_carlo_->set_out_file_tag( out_file_tag );
	rna_fragment_monte_carlo_->set_native_pose( get_native_pose() );
	rna_fragment_monte_carlo_->set_denovo_scorefxn( denovo_scorefxn_ );
	rna_fragment_monte_carlo_->set_hires_scorefxn( hires_scorefxn_ );

	std::clock_t start_time = clock();

	rna_fragment_monte_carlo_->set_refine_pose( false );
	rna_fragment_monte_carlo_->set_is_rna_and_protein( rna_params_->is_rna_and_protein() ); // need to know this for the high resolution stuff
	if ( options_->filter_vdw() ) rna_fragment_monte_carlo_->set_vdw_grid( vdw_grid_ );
	/*if ( !refine_pose || options_->refine_native_get_good_FT() ) */ rna_fragment_monte_carlo_->set_rna_de_novo_pose_initializer( rna_de_novo_pose_initializer ); // only used for resetting fold-tree & cutpoints on each try.
	rna_fragment_monte_carlo_->set_all_lores_score_final( all_lores_score_final_ );  // used for filtering.
	// if ( outputter != nullptr ) rna_fragment_monte_carlo_->set_outputter( outputter); // accumulate stats in histogram.

	// This calls apply on the rna_fragment_monte_carlo_ object, to which it
	// has a shared_ptr
	// denovo_job_distributor->apply( pose );
	rna_fragment_monte_carlo_->apply( pose );

	// extend
	all_lores_score_final_ = rna_fragment_monte_carlo_->all_lores_score_final(); // might have been updated, used for filtering.
	if ( options_->output_lores_silent_file() ) {
		add_number_base_pairs( *(rna_fragment_monte_carlo_->lores_pose()) );
		if ( get_native_pose() ) {
			add_number_native_base_pairs( *(rna_fragment_monte_carlo_->lores_pose()), *get_native_pose() );
		}

		protocols::jd2::output_intermediate_pose(
			*(rna_fragment_monte_carlo_->lores_pose()), "LORES");
	}
	// outputter = rna_fragment_monte_carlo_->outputter();

	// std::string const out_file_name = out_file_tag + ".pdb";
	// if ( options_->dump_pdb() ) dump_pdb( pose,  out_file_name );

	if ( options_->save_times() ) setPoseExtraScore( pose, "time", static_cast< Real >( clock() - start_time ) / CLOCKS_PER_SEC );

	// jd2 does this.
	// align_and_output_to_silent_file( pose, options_->silent_file(), out_file_tag );
	if ( options_->use_chem_shift_data() ) add_chem_shift_info( pose);

	if ( get_native_pose() ) {
		Real const rmsd = rna_fragment_monte_carlo_->get_rmsd_no_superimpose( pose );
		TR << "All atom rmsd: " << rmsd << " for this model." << std::endl;
		setPoseExtraScore( pose, "rms", rmsd );

		Real const rmsd_stems = rna_fragment_monte_carlo_->get_rmsd_stems_no_superimpose( pose );
		if ( rmsd_stems > 0.0 ) {
			TR << "All atom rmsd over stems: " << rmsd_stems << " for this model." << std::endl;
			setPoseExtraScore( pose, "rms_stem", rmsd_stems );
		}
	}

	add_number_base_pairs( pose );
	if ( get_native_pose() ) {
		add_number_native_base_pairs( pose, *get_native_pose() );
	}

	// hopefully these will end up in silent file...
	if ( options_->output_filters() && ( rna_fragment_monte_carlo_ ) ) {
		setPoseExtraScore( pose, "lores_early", rna_fragment_monte_carlo_->lores_score_early() );
		if ( options_->minimize_structure() ) setPoseExtraScore( pose, "lores_final", rna_fragment_monte_carlo_->lores_score_final() );
	}
}

void
RNA_DeNovoProtocolMover::add_chem_shift_info( core::pose::Pose & not_very_const_pose ) const {

	using namespace core::scoring;
	using namespace core::pose;

	runtime_assert( options_->use_chem_shift_data() );

	pose::Pose chem_shift_pose = not_very_const_pose; //HARD COPY SLOW!
	core::scoring::ScoreFunctionOP temp_scorefxn( new ScoreFunction );
	temp_scorefxn->set_weight( rna_chem_shift  , 1.00 );
	(*temp_scorefxn)(chem_shift_pose);

	EnergyMap const & energy_map = chem_shift_pose.energies().total_energies();

	Real const rosetta_chem_shift_score= energy_map[ rna_chem_shift ];

	//This statement should be very fast except possibly the 1st call.
	core::scoring::rna::chemical_shift::RNA_ChemicalShiftPotential const &
		rna_chemical_shift_potential( core::scoring::ScoringManager::
		get_instance()->get_RNA_ChemicalShiftPotential() );

	core::Size const num_chem_shift_data_points=rna_chemical_shift_potential.get_total_exp_chemical_shift_data_points();

	//rosetta_chem_shift_score --> Sum_square chemical_shift deviation.

	Real const chem_shift_RMSD=sqrt( rosetta_chem_shift_score /
		float(num_chem_shift_data_points) );

	setPoseExtraScore( not_very_const_pose, "chem_shift_RMSD", chem_shift_RMSD );
	//silent_struct.add_energy( "chem_shift_RMSD", chem_shift_RMSD);

	// silent_struct.add_energy( "num_chem_shift_data",
	//  float(num_chem_shift_data_points) );

	setPoseExtraScore( not_very_const_pose, "num_chem_shift_data", num_chem_shift_data_points );

	// if ( silent_struct.has_energy("rna_chem_shift")==false ) {
	//  //If missing this term, then the rna_chem_shift weight is probably
	//  //zero in the weight_file.
	//  silent_struct.add_energy( "rna_chem_shift", 0.0);
	// }
}

void
RNA_DeNovoProtocolMover::set_input_initialization_pdbs(
	utility::vector1< std::string > const & input_initialization_pdbs
) {
	input_initialization_pdbs_ = input_initialization_pdbs;
	options_->set_chunk_initialization_pdb_files( input_initialization_pdbs_  );
}

void RNA_DeNovoProtocolMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & //data
) {
	// We use this pattern because we want users (say, setting this up through pyrosetta)
	// to be able to echo this setup pattern with public setters.
	auto protocol_options = utility::pointer::make_shared< options::RNA_DeNovoProtocolOptions >();
	set_protocol_options( protocol_options );

	if ( tag->hasOption( "fasta_files" ) ) {
		set_fasta_files( utility::string_split_simple( tag->getOption< std::string >( "fasta_files" ), ',' ) );
	}
	if ( tag->hasOption( "input_pdbs" ) ) {
		set_input_pdbs( utility::string_split_simple( tag->getOption< std::string >( "input_pdbs" ), ',' ) );
	}
	if ( tag->hasOption( "input_silent_files" ) ) {
		set_input_silent_files( utility::string_split_simple( tag->getOption< std::string >( "input_silent_files" ), ',' ) );
	}
	if ( tag->hasOption( "input_sequence_strings" ) ) {
		set_sequence_strings( utility::string_split_simple( tag->getOption< std::string >( "input_sequence_strings" ), ',' ) );
	}
	if ( tag->hasOption( "minimize_rna" ) ) {
		set_minimize_rna( tag->getOption< bool >( "minimize_rna" ) );
	}
	if ( tag->hasOption( "helical_substructs" ) ) {
		set_helical_substructs( utility::string_split_simple( tag->getOption< std::string >( "helical_substructs" ), ',' ) );
	}
	if ( tag->hasOption( "dock_chunks" ) ) {
		set_dock_chunks( utility::string_split_simple( tag->getOption< std::string >( "dock_chunks" ), ',' ) );
	}
	if ( tag->hasOption( "initial_structures" ) ) {
		set_input_initialization_pdbs( utility::string_split_simple( tag->getOption< std::string >( "initial_structures" ), ',' ) );
	}

	// We are editing options_, but we installed this handle using a public setter,
	// so this is valid. PyRosetta users could do this.
	protocol_options->initialize_from_tag( tag );
	// protocol_options->initialize_from_options( basic::options::option );

	set_residue_type_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );

	// handled above, using the tag
	// initialize_inputs_from_options( option );

	// in pyrosetta, you could also use the de_novo_setup_from_command_line function
	// or build yourself an OptionCollection. This is a work in progress; more should
	// be exposed soon.
	// de_novo_setup_from_options( basic::options::option );
	de_novo_setup_from_tag( tag );

	// if output_res_num supplied, this will change PDBInfo numbering & chain.
	// set_output_res_and_chain( *pose_, basic::options::option[ basic::options::OptionKeys::rna::denovo::output_res_num ].resnum_and_chain() );
}

void
RNA_DeNovoProtocolMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	// TO DO!
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "fasta_files", xs_string, "Input FASTA files", "" )
		+ XMLSchemaAttribute::attribute_w_default( "input_pdbs", xs_string, "Input PDB files", "" )
		+ XMLSchemaAttribute::attribute_w_default( "input_silent_files", xs_string, "Input silent files", "" )
		+ XMLSchemaAttribute::attribute_w_default( "input_sequence_strings", xs_string, "Input sequence strings", "" )
		// this one is taken care of elsewhere!
		// + XMLSchemaAttribute::required_attribute( "minimize_rna", xsct_rosetta_bool, "Minimize RNA" )
		+ XMLSchemaAttribute::attribute_w_default( "helical_substructs", xs_string, "PDBs of helical substructures", "" )
		// this one is taken care of elsewhere!
		// + XMLSchemaAttribute::attribute_w_default( "dock_chunks", xs_string, "PDBs of RNA to be docked as rigid bodies", "" )
		+ XMLSchemaAttribute::attribute_w_default( "initial_structures", xs_string, "Initial locations for DRRAFTER", "" )
		+ XMLSchemaAttribute::attribute_w_default( "offset", xsct_non_negative_integer, "Sequence numbering offset", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "edensity_mapfile", xs_string, "Electron density mapfile", "" )
		// handled elsewhere?
		// + XMLSchemaAttribute::attribute_w_default( "rna_protein_docking", xsct_rosetta_bool, "RNA protein docking", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "virtual_anchor", xs_string, "Virtual anchor", "" )
		+ XMLSchemaAttribute::attribute_w_default( "cutpoint_open", xs_string, "Cutpoint open", "" ) // a tag-list
		+ XMLSchemaAttribute::attribute_w_default( "secstruct", xs_string, "Secstruct in dot-bracket notation", "" )
		+ XMLSchemaAttribute::attribute_w_default( "secstruct_file", xs_string, "File containing secstruct in dot-bracket notation", "" )
		+ XMLSchemaAttribute::attribute_w_default( "secstruct_general", xs_string, "General secstruct in dot-bracket notation", "" )
		+ XMLSchemaAttribute::attribute_w_default( "secstruct_general_file", xs_string, "File containing general secstruct in dot-bracket notation", "" )
		+ XMLSchemaAttribute::attribute_w_default( "secstruct_legacy", xs_string, "Secstruct legacy", "" )
		+ XMLSchemaAttribute::attribute_w_default( "input_res", xs_string, "Input residues", "" )
		+ XMLSchemaAttribute::attribute_w_default( "input_silent_res", xs_string, "Input silent residues", "" )
		+ XMLSchemaAttribute::attribute_w_default( "obligate_pair", xs_string, "Obligate pair", "" )
		+ XMLSchemaAttribute::attribute_w_default( "obligate_pair_explicit", xs_string, "Obligate pair (expliict)", "" )
		// cared for elsewhere?
		// + XMLSchemaAttribute::attribute_w_default( "extra_minimize_res", xs_string, "", "" )
		+ XMLSchemaAttribute::attribute_w_default( "working_res", xs_string, "Working residues", "" )

		+ XMLSchemaAttribute::attribute_w_default( "cutpoint_closed", xs_string, "Cutpoints closed", "" )
		+ XMLSchemaAttribute::attribute_w_default( "cyclize_res", xs_string, "Cyclize res", "" )
		+ XMLSchemaAttribute::attribute_w_default( "twoprime_res", xs_string, "2' connection residues", "" )
		+ XMLSchemaAttribute::attribute_w_default( "fiveprime_cap_res", xs_string, "5' capped residues", "" )
		+ XMLSchemaAttribute::attribute_w_default( "block_stack_above_res", xs_string, "Block stacking above these residues", "" )
		+ XMLSchemaAttribute::attribute_w_default( "block_stack_below_res", xs_string, "Block stacking below these residues", "" )
		// addressed earlier.
		//+ XMLSchemaAttribute::attribute_w_default( "virtual_anchor", xs_string, "", "" )
		+ XMLSchemaAttribute::attribute_w_default( "in_path", xs_string, "Input path", "" )
		// handled elsewhere
		// + XMLSchemaAttribute::attribute_w_default( "lores_scorefxn", xs_string, "", "" )
		+ XMLSchemaAttribute::attribute_w_default( "native", xs_string, "Native PDB file", "" )
		// handled elsewhere
		//+ XMLSchemaAttribute::attribute_w_default( "edensity_mapfile", xs_string, "", "" )
		+ XMLSchemaAttribute::attribute_w_default( "working_native", xs_string, "Working native PDB file", "" )
		+ XMLSchemaAttribute::attribute_w_default( "cst_file", xs_string, "Constraint file", "" )
		+ XMLSchemaAttribute::attribute_w_default( "rna_data_file", xs_string, "RNA data file", "" )
		+ XMLSchemaAttribute::attribute_w_default( "remove_pair", xs_string, "Remove pairs", "" )
		+ XMLSchemaAttribute::attribute_w_default( "remove_obligate_pair", xs_string, "Remove obligate pairs", "" )
		+ XMLSchemaAttribute::attribute_w_default( "chain_connection", xs_string, "Chain connection", "" )
		+ XMLSchemaAttribute::attribute_w_default( "get_fold_tree_from_silent_file", xs_string, "Silent file for FT setup", "" )
		+ XMLSchemaAttribute::attribute_w_default( "fold_tree_from_silent_file_tag", xs_string, "Tag to base FT on", "" )

		// handled earlier
		// + XMLSchemaAttribute::attribute_w_default( "extra_minimize_res", xs_string, "", "" )
		// + XMLSchemaAttribute::attribute_w_default( "extra_minimize_chi_res", xs_string, "", "" )
		// + XMLSchemaAttribute::attribute_w_default( "output_jump_res", xs_string, "", "" )

		+ XMLSchemaAttribute::attribute_w_default( "force_syn_chi_res_list", xs_string, "Res list for forcing syn chis", "" )
		+ XMLSchemaAttribute::attribute_w_default( "force_anti_chi_res_list", xs_string, "Res list for forcing anti chis", "" );


	core::import_pose::options::RNA_DeNovoProtocolOptions::list_attributes( attlist );

	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Replaces current PDB with a model ready for RNA fragment assembly, preserving "
		"the style of rna_denovo -- at least, its most common features -- as exactly as "
		"possible.",
		attlist );
}

void RNA_DeNovoProtocolMover::set_input_pdbs( utility::vector1< std::string > const & input_pdbs ) {
	input_pdbs_ = input_pdbs;
	runtime_assert( options_ );
	// This must propagate to options.
	options_->set_chunk_pdb_files( input_pdbs_ );
}

void RNA_DeNovoProtocolMover::set_input_silent_files( utility::vector1< std::string > const & input_silent_files ) {
	input_silent_files_ = input_silent_files;
	runtime_assert( options_ );
	// This must propagate to options.
	options_->set_chunk_silent_files( input_silent_files_ );
}

void RNA_DeNovoProtocolMover::set_align_pdb( std::string const & align_pdb ) {
	runtime_assert( options_ );
	// This must propagate to options.
	options_->set_align_pdb( align_pdb );
}

void RNA_DeNovoProtocolMover::set_nstruct( Size const nstruct ) {
	runtime_assert( options_ );
	// This must propagate to options.
	options_->set_nstruct( nstruct );
}

void RNA_DeNovoProtocolMover::set_silent_file( std::string const & silent_file ) {
	runtime_assert( options_ );
	// This must propagate to options.
	options_->set_silent_file( silent_file );
}


///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
// Following is directly adapted from rna_denovo_setup.py
//
// The main update is that for most of the script, we previously
//  kept track of 'conventional' numbers & chains for most
//  residues (working_res, cutpoint_open, extra_minimize_res, etc.)
//
// Now most of those values are held in 'full model' numbering.
//
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocolMover::de_novo_setup_from_command_line()
{
	de_novo_setup_from_options( basic::options::option );
}

void
RNA_DeNovoProtocolMover::list_options_read( utility::options::OptionKeyList & opts )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	opts + OptionKeys::rna::denovo::offset
		+ OptionKeys::full_model::cutpoint_open
		+ full_model::cutpoint_open
		+ full_model::working_res
		+ full_model::cutpoint_closed
		+ full_model::cyclize
		+ full_model::twoprime
		+ full_model::fiveprime_cap
		+ full_model::rna::block_stack_above_res
		+ full_model::rna::block_stack_below_res
		+ OptionKeys::rna::denovo::minimize::extra_minimize_res
		+ OptionKeys::rna::denovo::minimize::extra_minimize_chi_res
		+ in::file::input_res
		+ OptionKeys::rna::denovo::input_silent_res
		+ OptionKeys::rna::denovo::virtual_anchor
		+ OptionKeys::rna::denovo::obligate_pair
		+ OptionKeys::rna::denovo::remove_pair
		+ OptionKeys::rna::denovo::remove_obligate_pair
		+ OptionKeys::rna::denovo::out::output_jump_res
		+ OptionKeys::rna::denovo::secstruct
		+ OptionKeys::rna::denovo::secstruct_file
		+ OptionKeys::rna::denovo::secstruct_general
		+ OptionKeys::rna::denovo::secstruct_general_file
		+ in::file::native
		+ OptionKeys::rna::denovo::working_native
		+ OptionKeys::rna::denovo::obligate_pair_explicit
		+ OptionKeys::rna::denovo::chain_connection
		+ OptionKeys::rna::denovo::secstruct_legacy
		+ in::path::path
		+ OptionKeys::rna::denovo::lores_scorefxn
		+ OptionKeys::rna::denovo::refine_native
		+ OptionKeys::rna::denovo::working_native
		+ OptionKeys::rna::denovo::refine_native
		+ OptionKeys::rna::denovo::working_native
		+ OptionKeys::constraints::cst_file
		+ OptionKeys::rna::data_file
		+ OptionKeys::rna::denovo::minimize_rna
		+ full_model::rna::force_syn_chi_res_list
		+ full_model::rna::force_anti_chi_res_list;

	opts + OptionKeys::out::nstruct
		+ OptionKeys::out::file::silent
		+ OptionKeys::rna::denovo::tag
		+ OptionKeys::rna::denovo::out::output_lores_silent_file
		+ OptionKeys::rna::denovo::output_filters
		+ OptionKeys::out::overwrite
		+ OptionKeys::in::file::s
		+ OptionKeys::in::file::silent
		+ OptionKeys::in::file::input_res
		+ OptionKeys::rna::denovo::refine_silent_file
		+ OptionKeys::rna::denovo::lores_scorefxn
		+ OptionKeys::in::file::silent_struct_type
		+ OptionKeys::rna::denovo::out::binary_output
		+ OptionKeys::out::save_times
		+ OptionKeys::rna::denovo::use_legacy_setup
		+ OptionKeys::rna::denovo::cst_gap
		+ OptionKeys::rna::denovo::dump_stems;
}

void
RNA_DeNovoProtocolMover::initialize_sequence_information(
	int const offset,
	bool const edensity_mapfile_provided,
	bool const rna_protein_docking,
	bool const virtual_anchor_provided,
	bool const cutpoint_open_provided,
	std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > > const & possible_cutpoint_open_numbering,
	core::pose::full_model_info::FullModelParametersOP & full_model_parameters,
	utility::vector1< core::Size > & cutpoint_open_in_full_model )
{
	using namespace basic::options;
	using namespace core::chemical;
	using namespace core::pose::full_model_info;

	if ( !fasta_files_.empty() ) {
		// use fasta readin developed for stepwise application -- also reads in
		// numbers & chains based on fasta header lines.
		if ( sequence_strings_.size() > 0 ) {
			TR << fasta_files_ << std::endl;
			TR << sequence_strings_ << std::endl;
			utility_exit_with_message( "Cannot specify both -sequence and -fasta" );
		}
		if ( fasta_files_.size() != 1 ) utility_exit_with_message( "Please specify exactly one fasta file." );
		bool append_virtual = false;
		// AMW: hey Kalli, is this where we would want to also append a virtual... say, if -please_give_new_FT is true?
		if ( edensity_mapfile_provided || rna_protein_docking || virtual_anchor_provided ) {
			// if a density map was supplied by the user, the virtual residue would have been appended by default anyway
			append_virtual = true;
		}
		full_model_parameters = core::import_pose::get_sequence_information( fasta_files_[ 1 ], cutpoint_open_in_full_model, append_virtual );
		if ( offset != 0 ) {
			vector1< int > new_numbering = full_model_parameters->conventional_numbering();
			for ( Size n = 1; n <= new_numbering.size(); n++ ) { new_numbering[ n ] += offset; }
			full_model_parameters->set_conventional_numbering( new_numbering );
		}

	} else {
		if ( edensity_mapfile_provided || rna_protein_docking ) {

			utility_exit_with_message( "Please use fasta files to specify sequence if providing a density map or doing RNA/protein docking");
		}
		// basic read-in of sequence from command line
		if ( sequence_strings_.size() == 0 ) utility_exit_with_message( "Must specify -sequence or -fasta" );
		std::string sequence( sequence_strings_[1] );
		for ( Size n = 2; n <= sequence_strings_.size(); n++ ) sequence += std::string( " " + sequence_strings_[ n ] );
		cutpoint_open_in_full_model = core::sequence::strip_spacers( sequence );
		std::map< Size, std::string > non_standard_residue_map = core::sequence::parse_out_non_standard_residues( sequence );
		vector1< int > res_numbers_in_pose;
		for ( Size n = 1; n <= sequence.size(); n++ ) res_numbers_in_pose.push_back( n + offset );
		core::import_pose::get_extra_cutpoints_from_names( sequence.size(), cutpoint_open_in_full_model, non_standard_residue_map );
		full_model_parameters = utility::pointer::make_shared< FullModelParameters >( sequence, cutpoint_open_in_full_model, res_numbers_in_pose );
		full_model_parameters->set_non_standard_residue_map( non_standard_residue_map );
		Size chain_num( 1 ), res( 0 );
		utility::vector1< char > chains; utility::vector1< Size > resnum;
		for ( Size n = 1; n <= sequence.size(); n++ ) {
			chains.push_back( chr_chains[ (chain_num - 1)  % chr_chains.size() ] );
			resnum.push_back( ++res );
			if ( cutpoint_open_in_full_model.has_value( n ) ) {
				chain_num++; res = 0;
			}
		}
		full_model_parameters->set_conventional_chains( chains );
		full_model_parameters->set_conventional_numbering( resnum );
	}

	// Finish initializing cutpoints
	if ( cutpoint_open_provided ) {
		utility::vector1< core::Size > cutpoint_open_in_full_model_from_command_line  =
			full_model_parameters->conventional_to_full( possible_cutpoint_open_numbering );
		// check whether these residues are already in cutpoint_open_in_full_model
		// if not, add them
		for ( Size n = 1; n<=cutpoint_open_in_full_model_from_command_line.size(); ++n ) {
			if ( !cutpoint_open_in_full_model.has_value( cutpoint_open_in_full_model_from_command_line[ n ] ) ) {
				cutpoint_open_in_full_model.push_back( cutpoint_open_in_full_model_from_command_line[ n ] );
			}
		}
	}
}

void
RNA_DeNovoProtocolMover::check_secstructs(
	bool const model_with_density,
	bool const rna_protein_docking,
	std::string const & sequence,
	core::pose::rna::RNA_SecStruct const & secstruct,
	core::pose::rna::RNA_SecStruct const & secstruct_general )
{

	//secstruct.check_compatible_with_sequence( sequence, true  /*check_complementarity*/ );

	// At this point the virtual residue for density scoring has already been added
	// so there's a mismatch is sequence and secondary structure length
	// so we need to check against the sequence, minus the last residue
	if ( model_with_density || rna_protein_docking ) {
		std::string sequence_without_last;
		for ( core::Size i=0; i<sequence.size()-1; ++i ) {
			sequence_without_last+=sequence[i];
		}
		// And check that the last residue really is virtual?
		// AMW: we can't do this: the secstruct is too long to be compatible with a shorter sequence!
		/*
		std::cout << "SEQUENCE WITHOUT LAST" << std::endl;
		std::cout << sequence_without_last << std::endl;
		secstruct.check_compatible_with_sequence( sequence_without_last, true );
		*/
	} else {
		secstruct.check_compatible_with_sequence( sequence, true  /*check_complementarity*/ );
	}
	secstruct_general.check_compatible_with_sequence( sequence, false /*check_complementarity*/ );

	if ( !secstruct_general.blank() && !options_->bps_moves() ) utility_exit_with_message("cannot supply secstruct_general without bps_moves");
}

void
RNA_DeNovoProtocolMover::input_pdb_numbering_setup(
	std::string const & sequence,
	core::pose::full_model_info::FullModelParametersOP const & full_model_parameters,
	utility::vector1< core::Size > & input_res,
	utility::vector1< utility::vector1< core::Size > > & helical_substruct_res,
	utility::vector1< core::Size > & dock_chunks_res,
	utility::vector1< utility::vector1< int > > & resnum_list,
	utility::vector1< core::Size > & input_res_user_defined,
	Size & input_res_user_defined_count
) {
	for ( std::string const & pdb : input_pdbs_ ) {
		std::string pdb_seq;
		vector1< int > resnum;
		vector1< char >  chain;
		vector1< std::string >  segid;
		get_seq_and_resnum( pdb, pdb_seq, resnum, chain, segid );
		vector1< Size > resnum_in_full_model;

		if ( input_res_user_defined_count + resnum.size() <= input_res_user_defined.size() ) {
			// input res could have come from user after flag -input_res
			for ( Size q = 1; q <= resnum.size(); q++ ) {
				input_res_user_defined_count++;
				resnum_in_full_model.push_back( input_res_user_defined[ input_res_user_defined_count ] );
			}
		} else {
			// figure out residue numbers from PDB resnum & chain.
			resnum_in_full_model = full_model_parameters->conventional_to_full( std::make_tuple( resnum, chain, segid ) );
		}

		std::string actual_seq = "";
		for ( Size q = 1; q <= resnum_in_full_model.size(); q++ ) {
			Size const i = resnum_in_full_model[ q ];
			if ( input_res.has_value( i ) )  TR << TR.Red << "WARNING! Input residue " << resnum[q] << " " << chain[q] << " exists in two pdb files!!" << std::endl;
			actual_seq += sequence[ i - 1 ];
			input_res.push_back( i );
		}
		if ( pdb_seq != actual_seq ) {
			TR << TR.Red << "   pdb_seq: " << pdb_seq << std::endl;
			TR << TR.Red << "target_seq: " << actual_seq << std::endl;
			for ( Size q = 1; q <= resnum_in_full_model.size(); q++ ) {
				if ( q > actual_seq.size() ) {
					TR << "mismatch in length beyond " << q << std::endl;
					break;
				}
				if ( pdb_seq[q-1] != actual_seq[q-1] ) {
					Size n( resnum_in_full_model[ q ] );
					TR << "mismatch in sequence: pdb  " << pdb_seq[q-1] << " vs target " << actual_seq[q-1] << " at " << full_model_parameters->full_to_conventional( n ) << std::endl;
				}
			}
			utility_exit_with_message("The sequence in "+pdb+" does not match target sequence!!");
		}
		resnum_list.push_back( resnum_in_full_model );

		// check whether this is one of the helical_substructs
		if ( helical_substructs_.contains( pdb ) ) {
			helical_substruct_res.push_back( resnum_in_full_model );
		}
		// check whether this is one of the dock chunks
		if ( dock_chunks_.contains( pdb ) ) {
			for ( auto const r : resnum_in_full_model ) {
				dock_chunks_res.push_back( r );
			}
		}
	}
}

void
RNA_DeNovoProtocolMover::input_initialization_pdb_numbering_setup(
	std::string const & sequence,
	core::pose::full_model_info::FullModelParametersOP const & full_model_parameters,
	utility::vector1< core::Size > & input_res_initialize
) {
	///////
	// Also set up the input res for the initialization pdbs
	///////

	vector1< vector1< int > > resnum_list_initialize;
	//Size input_res_initialize_user_defined_count( 0 );
	for ( auto const & pdb : input_initialization_pdbs_ ) {
		std::string pdb_seq;
		vector1< int > resnum;
		vector1< char >  chain;
		vector1< std::string >  segid;
		get_seq_and_resnum( pdb, pdb_seq, resnum, chain, segid );
		vector1< Size > resnum_in_full_model;
		// don't allow input res from -input_res flag for initialization structures

		// figure out residue numbers from PDB resnum & chain.
		resnum_in_full_model = full_model_parameters->conventional_to_full( std::make_tuple( resnum, chain, segid ) );

		std::string actual_seq = "";
		for ( Size q = 1; q <= resnum_in_full_model.size(); q++ ) {
			Size const i = resnum_in_full_model[ q ];
			if ( input_res_initialize.has_value( i ) )  TR << TR.Red << "WARNING! Input residue " << resnum[q] << " " << chain[q] << " exists in two pdb files!!" << std::endl;
			actual_seq += sequence[ i - 1 ];
			input_res_initialize.push_back( i );
		}
		if ( pdb_seq != actual_seq ) {
			TR << TR.Red << "   pdb_seq: " << pdb_seq << std::endl;
			TR << TR.Red << "target_seq: " << actual_seq << std::endl;
			for ( Size q = 1; q <= resnum_in_full_model.size(); q++ ) {
				if ( q > actual_seq.size() ) {
					TR << "mismatch in length beyond " << q << std::endl;
					break;
				}
				if ( pdb_seq[q-1] != actual_seq[q-1] ) {
					Size n( resnum_in_full_model[ q ] );
					TR << "mismatch in sequence: pdb  " << pdb_seq[q-1] << " vs target " << actual_seq[q-1] << " at " << full_model_parameters->full_to_conventional( n ) << std::endl;
				}
			}
			utility_exit_with_message("The sequence in "+pdb+" does not match target sequence!!");
		}
		resnum_list_initialize.push_back( resnum_in_full_model );
	}
}

void
RNA_DeNovoProtocolMover::input_silent_numbering_setup(
	std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > > const & input_silent_res_resnum_and_chain,
	std::string const & sequence,
	core::pose::full_model_info::FullModelParametersOP const & full_model_parameters,
	utility::vector1< core::Size > const & input_res_user_defined,
	Size const & input_res_user_defined_count,
	utility::vector1< core::Size > & input_res,
	utility::vector1< utility::vector1< core::Size > > & helical_substruct_res,
	utility::vector1< core::Size > & dock_chunks_res,
	utility::vector1< utility::vector1< int > > & resnum_list
) {
	using namespace basic::options;
	////////////////////
	// Step 6
	////////////////////
	// go through silent files for each residue, and sequences match up.
	vector1< Size > input_silent_res_user_defined =
		full_model_parameters->conventional_to_full( input_silent_res_resnum_and_chain );
	Size input_silent_res_user_defined_count( 0 );
	vector1< Size > input_silent_res;
	if ( input_silent_res_user_defined.size() > 0 ) {
		// must have run through input_res.
		runtime_assert( input_res_user_defined_count == input_res_user_defined.size() );
		input_silent_res = input_silent_res_user_defined;
	} else {
		input_silent_res = input_res_user_defined;
		input_silent_res_user_defined_count = input_res_user_defined_count;
	}

	for ( std::string const & silent : input_silent_files_ ) {
		std::string seq = get_silent_seq( silent );
		Size len_seq = seq.size();
		vector1< Size > resnum_in_full_model;
		if ( ( input_silent_res_user_defined_count + len_seq ) <= input_silent_res.size() ) {
			// input res could have come from user after flag -input_res or -input_silent_res
			for ( Size q = 1; q <= len_seq; q++ ) {
				input_silent_res_user_defined_count++;
				resnum_in_full_model.push_back( input_silent_res[ input_silent_res_user_defined_count ] );
			}
		} else {
			vector1<Size> input_silent_res_from_file = full_model_parameters->conventional_to_full( get_silent_resnum( silent ) );
			for ( Size const res : input_silent_res_from_file ) {
				resnum_in_full_model.push_back( res );
			}
		}
		runtime_assert( resnum_in_full_model.size() == len_seq);

		std::string actual_seq;
		for ( Size q = 1; q <= len_seq; q++ ) {
			Size const i = resnum_in_full_model[ q ];
			if ( input_res.has_value( i ) )  TR << TR.Red << "WARNING! Input residue " << i << " exists in two pdb/silent files!!" << std::endl;
			actual_seq += sequence[ i - 1 ];
			input_res.push_back( i );
		}

		if ( seq != actual_seq ) {
			TR << TR.Red << seq << std::endl;
			TR << TR.Red << actual_seq << std::endl;
			utility_exit_with_message("The sequence in "+silent+" does not match input sequence!!");
		}

		resnum_list.push_back( resnum_in_full_model );

		// check whether this is one of the helical_substructs
		if ( helical_substructs_.contains( silent ) ) {
			helical_substruct_res.push_back( resnum_in_full_model );
		}

		// check whether this is one of the dock chunks
		if ( dock_chunks_.contains( silent ) ) {
			for ( auto const r : resnum_in_full_model ) {
				dock_chunks_res.push_back( r );
			}
		}

	}
	runtime_assert( input_silent_res_user_defined_count == input_silent_res.size() );
}

void
RNA_DeNovoProtocolMover::input_numbering_setup(
	std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > > const & input_res_resnum_and_chain,
	std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > > const & input_silent_res_resnum_and_chain,
	std::string const & sequence,
	core::pose::full_model_info::FullModelParametersOP const & full_model_parameters,
	utility::vector1< core::Size > & input_res,
	utility::vector1< utility::vector1< int > > & resnum_list,
	utility::vector1< core::Size > & input_res_initialize,
	utility::vector1< utility::vector1< core::Size > > & helical_substruct_res,
	utility::vector1< core::Size > & dock_chunks_res
) {
	using namespace basic::options;

	vector1< Size > input_res_user_defined =
		full_model_parameters->conventional_to_full( input_res_resnum_and_chain );
	Size input_res_user_defined_count( 0 );
	input_pdb_numbering_setup(
		sequence, full_model_parameters, input_res, helical_substruct_res, dock_chunks_res, resnum_list,
		input_res_user_defined, input_res_user_defined_count
	);
	input_initialization_pdb_numbering_setup(
		sequence, full_model_parameters, input_res_initialize
	);

	input_silent_numbering_setup(
		input_silent_res_resnum_and_chain, sequence, full_model_parameters, input_res_user_defined, input_res_user_defined_count,
		input_res, helical_substruct_res, dock_chunks_res, resnum_list
	);
}

void
RNA_DeNovoProtocolMover::setup_obligate_pair(
	std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > > const & extra_minimize_res_rc,
	core::pose::full_model_info::FullModelParametersOP const & full_model_parameters,
	utility::vector1< utility::vector1< int > > const & resnum_list,
	core::pose::rna::RNA_SecStruct const & secstruct,
	core::pose::rna::RNA_SecStruct const & secstruct_general,
	utility::vector1< core::Size > const & cutpoint_open_in_full_model,
	utility::vector1< core::Size > & obligate_pair,
	utility::vector1< std::string > & obligate_pair_explicit,
	utility::vector1< core::Size > & domain_map )
{
	using namespace basic::options;

	vector1< Size > extra_minimize_res =
		full_model_parameters->conventional_to_full( extra_minimize_res_rc );
	runtime_assert( obligate_pair_explicit.size() % 5 == 0 );
	vector1< Size > obligate_pair_explicit_full_model;
	for ( Size m = 0; m < obligate_pair_explicit.size()/5; m++ ) {
		vector1< int > resnum;
		vector1< char > chains;
		vector1< std::string > segids;
		vector1< Size > resnum_full;

		utility::get_resnum_and_chain_from_one_tag( obligate_pair_explicit[ 5*m + 1 ], resnum, chains, segids );
		resnum_full = full_model_parameters->conventional_to_full( std::make_tuple( vector1<int>( resnum ),
			vector1<char>( chains ), vector1< std::string >( segids ) ) );
		runtime_assert( resnum_full.size() == 1 );
		Size const pos1 = resnum_full[ 1 ];

		resnum.clear(); chains.clear();
		utility::get_resnum_and_chain_from_one_tag( obligate_pair_explicit[ 5*m + 2 ], resnum, chains, segids );
		resnum_full = full_model_parameters->conventional_to_full( std::make_tuple( vector1<int>( resnum ),
			vector1<char>( chains ), vector1< std::string >( segids ) ) );
		runtime_assert( resnum_full.size() == 1 );
		Size const pos2 = resnum_full[ 1 ];

		obligate_pair_explicit_full_model.push_back( pos1 );
		obligate_pair_explicit_full_model.push_back( pos2 );
	}


	vector1< char > const & conventional_chains = full_model_parameters->conventional_chains();
	// Go through each of the inputs, and look for broken chains -- will define obligate pairs.
	for ( Size n = 1; n <= resnum_list.size(); n++ ) {
		vector1< Size > resnum = resnum_list[ n ];
		//Find obligate pairs
		vector1< vector1< Size > > chunks;
		vector1< Size > curr_chunk;
		//Size i( 0 ), j( 0 );
		char c(' '), d( ' ' );
		for ( Size q = 1, i = 0, j = 0; q <= resnum.size(); q++ ) {
			i = resnum[ q ];
			c = conventional_chains[ resnum[ q ] ];
			domain_map[ i ] = n;
			if ( j > 0 &&
					( ( ( i - 1 ) != j ) || d != c || cutpoint_open_in_full_model.has_value( j )  ) ) {
				chunks.push_back(curr_chunk);
				curr_chunk.clear();
			}
			j = i;
			d = c;
			curr_chunk.push_back( i );
		}
		chunks.push_back(curr_chunk);

		Size n_jumps( 0 );
		// look to see if each chunk might already be connected to the next.
		vector1< bool > connected_to_next( chunks.size(), false );
		for ( Size i = 1; i <= chunks.size(); i++ ) {
			Size const i_next = ( i < chunks.size() ) ? ( i + 1 ) : 1;

			// first check if there is already an obligate_pair between these chunks:
			bool found_pair( false );
			for ( Size j = 1; j <= chunks[ i ].size(); j++ ) {
				for ( Size k = 1; k <= chunks[ i_next  ].size(); k++ ) {
					Size const segment1_end   = chunks[ i     ][ j ];
					Size const segment2_start = chunks[ i_next ][ k ];
					Size const new_pos1 = std::min( segment1_end, segment2_start );
					Size const new_pos2 = std::max( segment1_end, segment2_start );
					vector1< Size > new_pair = make_vector1( new_pos1, new_pos2 );
					if ( already_listed_in_obligate_pair( new_pair, obligate_pair, obligate_pair_explicit_full_model ) ) {
						found_pair = true; break;
					}
					if ( found_pair ) break;
				}
			}
			if ( found_pair ) {
				connected_to_next[ i ] = true;
				n_jumps += 1;
			}
		}

		for ( Size i = 1; i <= chunks.size(); i++ ) {
			//obligate pairs -- go from end of one chain to beginning of next chain.
			if ( n_jumps == chunks.size() - 1 ) break;
			if ( connected_to_next[ i ] ) continue;
			Size const i_next = ( i < chunks.size() ) ? ( i + 1 ) : 1;

			Size segment1_end   = chunks[ i     ][ chunks[ i ].size() ];
			Size segment2_start = chunks[ i_next ][ 1 ];

			// avoid placing jump positions at extra_minimize_res
			while ( extra_minimize_res.has_value( segment1_end )   && segment1_end > chunks[ i ][ 1 ] ) segment1_end--;
			while ( extra_minimize_res.has_value( segment2_start ) && segment2_start < chunks[ i_next ][ chunks[ i_next ].size() ] ) segment2_start++;

			Size const new_pos1 = std::min( segment1_end, segment2_start );
			Size const new_pos2 = std::max( segment1_end, segment2_start );
			obligate_pair.push_back( new_pos1 );
			obligate_pair.push_back( new_pos2 );
			//   TR << TR.Cyan << "Creating new obligate pair: " << obligate_pair << " for chunk with residues " << resnum << std::endl;
			n_jumps++;
		}
	}
	vector1< std::pair< Size, Size > > canonical_pairs = secstruct.base_pairs();
	vector1< std::pair< Size, Size > > general_pairs   = secstruct_general.base_pairs();
	for ( Size n = 1; n <= general_pairs.size(); n++ ) {
		std::pair< Size, Size > const & p = general_pairs[ n ];
		if ( !canonical_pairs.has_value( p ) &&
				( domain_map[ p.first ] == 0 || domain_map[ p.second ] == 0 ) ) {
			vector1< Size > const new_pair = make_vector1( p.first, p.second );
			if ( !already_listed_in_obligate_pair( new_pair, obligate_pair, obligate_pair_explicit_full_model ) ) {
				obligate_pair.push_back( p.first );
				obligate_pair.push_back( p.second );
			}
		}
	}
}

void
RNA_DeNovoProtocolMover::initial_pose_setup(
	std::string const & in_path,
	bool const lores_sfxn_provided,
	bool const native_provided,
	bool const edensity_mapfile_provided,
	bool const working_native_provided,
	std::string const & working_native,
	core::pose::full_model_info::FullModelParametersOP const & full_model_parameters,
	utility::vector1< core::Size > const & working_res,
	core::pose::Pose & full_pose,
	bool & is_rna_and_protein )
{
	using namespace basic::options;

	////////////////////////
	// working pose
	////////////////////////
	pose_ = utility::pointer::make_shared< Pose >();
	std::string const full_annotated_sequence = full_model_parameters->full_annotated_sequence();
	make_pose_from_sequence( full_pose, full_annotated_sequence, *rsd_set_ );
	set_output_res_and_chain( full_pose, std::make_tuple( full_model_parameters->conventional_numbering(),
		full_model_parameters->conventional_chains(), full_model_parameters->conventional_segids()  ) );

	// Check whether the sequence contains protein residues
	for ( core::Size r=1; r<=full_pose.total_residue(); ++r ) {
		if ( full_pose.residue( r ).is_protein() ) {
			is_rna_and_protein = true;
			break;
		}
	}

	//
	if ( is_rna_and_protein ) {
		if ( !lores_sfxn_provided ) {
			// set default low-res RNA/protein score function
			options_->set_lores_scorefxn( "rna/denovo/rna_lores_with_rnp" );
		}
	}

	// // Try to set stuff up for density scoring ?
	// if ( edensity_mapfile_provided ) {
	//  //  pose::addVirtualResAsRoot( full_pose );
	// }
	TR.Debug << "THE SEQUENCE:" << std::endl;
	TR.Debug << full_pose.sequence() << std::endl;

	////////////////////////
	// Working native pose
	////////////////////////
	//Read in native if it exists.
	if ( native_provided ) {
		// AMW TODO: later, let this pass align_pose on to the RNA_DeNovoProtocol
		// and the RNA_FragmentMonteCarlo. At the moment, this isn't necessary at
		// all, though.
		core::pose::PoseOP align_pose;
		core::import_pose::initialize_native_and_align_pose( native_pose_, align_pose, rsd_set_, pose_ );

		// if we're using density, append the virtual residue
		if ( edensity_mapfile_provided ) {
			core::pose::addVirtualResAsRoot( *native_pose_ );
		}

		// AMW: this function will ensure that it just copies the native pose if in fact working_res is 'everything'
		pdbslice( *native_pose_, working_res );
	} else if ( working_native_provided ) {
		native_pose_ = utility::pointer::make_shared< Pose >();
		core::import_pose::pose_from_file( *native_pose_, *rsd_set_, in_path + working_native , core::import_pose::PDB_file);
		if ( edensity_mapfile_provided ) {
			core::pose::addVirtualResAsRoot( *native_pose_ );
		}
	}

	if ( ! working_native_provided ) { // usually not defined by user
		pdbslice( *pose_, full_pose, working_res );
	} else {
		// there might still be issues with how the csts are set up here...?
		pose_ = native_pose_->clone();
	}

	TR.Debug << "THE SEQUENCE OF THE POSE:" << std::endl;
	TR.Debug << pose_->sequence() << std::endl;
}

void
RNA_DeNovoProtocolMover::constraint_setup(
	bool const constraint_file_provided,
	std::string const & constraint_file,
	utility::vector1< core::Size > const & working_res,
	core::pose::Pose & full_pose
) {
	using namespace basic::options;
	using namespace core::scoring::constraints;
	using namespace core::scoring::func;
	using namespace core::id;

	////////////////////////
	// working constraints.
	////////////////////////
	if ( constraint_file_provided ) {
		ConstraintSetOP cst_set( new ConstraintSet );
		//opts[ OptionKeys::constraints::force_pdb_info_mapping ].def( true ); // using option as global variable due to difficulty in dealing with static functions.
		cst_set = ConstraintIO::get_instance()->read_constraints( constraint_file, utility::pointer::make_shared< ConstraintSet >(), full_pose, true /*force_pdb_info_mapping*/ );
		full_pose.constraint_set( cst_set );
		id::SequenceMappingOP sequence_map( new id::SequenceMapping( working_res ) );
		sequence_map->reverse();
		pose_->constraint_set( full_pose.constraint_set()->remapped_clone( full_pose, *pose_, sequence_map ) );
	}

	if ( !options_->cst_gap() ) return;
	// also create a constraints file that will help close chains across 2-5 residue gaps...
	// this is a slight hack, but if it works, might be worth putting a term into Rosetta, as
	// well as automated handling of "working_res"
	for ( Size m = 1; m <= working_res.size(); m++ ) {
		Size const gap = working_res[m+1] - working_res[m];
		if ( gap < 1 || gap > 6 ) continue;

		Distance const stdev( 10.0 );
		Distance max_dist = gap * 5.0 + 4;
		Real const bonus = 200.0;
		//        cst_file_outstring +=  " O3* %d  C5* %d   FADE %6.3f  %6.3f  %6.3f %6.3f %6.3f \n" %
		//            ( cst_gap[0], cst_gap[1],  -stdev, max_dist, stdev, -1*bonus, bonus)
		FuncOP distance_func( new FadeFunc( -stdev, max_dist, stdev, -1.0*bonus, bonus ) );
		pose_->add_constraint( ConstraintCOP( utility::pointer::make_shared< AtomPairConstraint >(
			named_atom_id_to_atom_id( NamedAtomID( " O3'",m), *pose_),
			named_atom_id_to_atom_id( NamedAtomID( " C5'",m+1), *pose_),
			distance_func ) ) );
	}
}

void
RNA_DeNovoProtocolMover::refine_working_obligate_pairs(
	std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > > const & remove_obligate_pair_rc,
	core::pose::full_model_info::FullModelParametersOP const & full_model_parameters,
	utility::vector1< core::Size > const & obligate_pair,
	utility::vector1< std::string > & obligate_pair_explicit,
	utility::vector1< core::Size > const & working_res,
	utility::vector1< core::pose::rna::BasePair > & working_obligate_pairs
) {
	using namespace basic::options;
	vector1< Size > remove_obligate_pair =
		full_model_parameters->conventional_to_full( remove_obligate_pair_rc );
	runtime_assert( obligate_pair.size() % 2 == 0 );
	for ( Size m = 0; m < obligate_pair.size()/2; m++ ) {
		Size const pos1 = obligate_pair[ 2*m + 1 ];
		Size const pos2 = obligate_pair[ 2*m + 2 ];
		if ( !working_res.has_value( pos1 ) ) continue;
		if ( !working_res.has_value( pos2 ) ) continue;

		bool rm_pair( false );
		for ( Size q = 1; q <= remove_obligate_pair.size(); q++ ) {
			Size const rm_pos1 = remove_obligate_pair[ 2*q + 1 ];
			Size const rm_pos2 = remove_obligate_pair[ 2*q + 2 ];
			if ( rm_pos1 == pos1 && rm_pos2 == pos2 ) {
				rm_pair = true; break;
			}
			if ( rm_pos1 == pos2 && rm_pos2 == pos1 ) {
				rm_pair = true; break;
			}
		}
		if ( rm_pair ) continue;

		working_obligate_pairs.push_back( BasePair( working_res.index( pos1 ), working_res.index( pos2 ), ANY_BASE_EDGE, ANY_BASE_EDGE, ANY_BASE_DOUBLET_ORIENTATION ) );
	}

	runtime_assert( obligate_pair_explicit.size() % 5 == 0 );
	for ( Size m = 0; m < obligate_pair_explicit.size()/5; m++ ) {

		vector1< int > resnum;
		vector1< char > chains;
		vector1< std::string > segids;
		vector1< Size > resnum_full;

		utility::get_resnum_and_chain_from_one_tag( obligate_pair_explicit[ 5*m + 1 ], resnum, chains, segids );
		resnum_full = full_model_parameters->conventional_to_full( std::make_tuple( vector1<int>( resnum ),
			vector1<char>( chains ), vector1<std::string>( segids ) ) );
		runtime_assert( resnum_full.size() == 1 );
		Size const pos1 = resnum_full[ 1 ];
		if ( !working_res.has_value( pos1 ) ) continue;

		resnum.clear(); chains.clear();
		utility::get_resnum_and_chain_from_one_tag( obligate_pair_explicit[ 5*m + 2 ], resnum, chains, segids );
		resnum_full = full_model_parameters->conventional_to_full( std::make_tuple( vector1<int>( resnum ),
			vector1<char>( chains ), vector1<std::string>( segids ) ) );
		runtime_assert( resnum_full.size() == 1 );
		Size const pos2 = resnum_full[ 1 ];
		if ( !working_res.has_value( pos2 ) ) continue;

		runtime_assert( obligate_pair_explicit[ 5*m + 3 ].size() == 1 );
		BaseEdge const edge1 = get_edge_from_char( obligate_pair_explicit[ 5*m + 3 ][ 0 ] );

		runtime_assert( obligate_pair_explicit[ 5*m + 4 ].size() == 1 );
		BaseEdge const edge2 = get_edge_from_char( obligate_pair_explicit[ 5*m + 4 ][ 0 ] );

		runtime_assert( obligate_pair_explicit[ 5*m + 5 ].size() == 1 );
		char const o_char = obligate_pair_explicit[ 5*m + 5 ][ 0 ];

		BaseDoubletOrientation orientation;
		if ( o_char == 'C' || o_char == 'T' ) { // convert from leontis-westhof cis/trans to antiparallel/antiparallel
			orientation  = get_base_doublet_orientation_from_LW( edge1, edge2, get_LW_orientation_from_char( o_char ) );
		} else {
			orientation = get_orientation_from_char( o_char );
		}

		working_obligate_pairs.push_back( BasePair( working_res.index( pos1 ), working_res.index( pos2 ),
			edge1, edge2, orientation ) );
	}
}

void
RNA_DeNovoProtocolMover::setup_working_chain_connections(
	utility::vector1< std::string > const & chain_connections,
	core::pose::full_model_info::FullModelParametersOP const & full_model_parameters,
	utility::vector1< core::Size > const & working_res,
	utility::vector1< std::pair< utility::vector1< core::Size >, utility::vector1< core::Size > > > & working_chain_connections
) {
	using namespace basic::options;
	//////////////////////////////
	// working chain connections
	//////////////////////////////
	if ( chain_connections.empty() ) return;

	runtime_assert( chain_connections.size() >= 4 );
	if ( !chain_connections.has_value( "SET1" ) ) return;

	// new format is more flexible enough to allow for noncontiguous numbering in set1 vs. set2)
	Size which_set = 0;
	utility::vector1< core::Size > resnum1;
	utility::vector1< core::Size > resnum2;
	for ( Size k = 1; k <= chain_connections.size() + 1; k++ ) {
		if ( k == ( chain_connections.size() + 1 ) || chain_connections[ k ] == "SET1" ) { // at of previous block
			utility::vector1< core::Size > working_resnum1 = working_res_map( resnum1, working_res );
			utility::vector1< core::Size > working_resnum2 = working_res_map( resnum2, working_res );
			if ( working_resnum1.size() > 0 &&
					working_resnum2.size() > 0 ) {
				working_chain_connections.push_back( std::make_pair( working_resnum1, working_resnum2) );
			}
			resnum1.clear();
			resnum2.clear();
			which_set = 1;
			continue;
		}
		if ( chain_connections[ k ] == "SET2" ) {
			which_set = 2;
			continue;
		}
		runtime_assert( which_set > 0 );

		vector1< int > resnum;
		vector1< char > chains;
		vector1< std::string > segids;
		bool ok = utility::get_resnum_and_chain_from_one_tag( chain_connections[ k ], resnum, chains, segids );
		runtime_assert( ok );
		utility::vector1< core::Size > resnum_full = full_model_parameters->conventional_to_full(
			std::make_tuple( vector1< core::Size >( resnum ),
			utility::vector1< char >( chains ),
			utility::vector1< std::string >( segids ) ) );
		utility::vector1< core::Size > & resnum_for_set = ( which_set == 1 ) ? resnum1 : resnum2;
		for ( core::Size n = 1; n <= resnum_full.size(); n++ ) resnum_for_set.push_back( resnum_full[ n ] );
	}
}

void
RNA_DeNovoProtocolMover::setup_pairing_sets(
	utility::vector1< utility::vector1< std::pair< core::Size, core::Size > > > const & working_stems,
	utility::vector1< core::pose::rna::BasePair > & working_obligate_pairs,
	RNA_BasePairList & rna_pairing_list,
	utility::vector1< utility::vector1< core::Size > > & working_obligate_pair_sets,
	utility::vector1< utility::vector1< core::Size > > & working_stem_pairing_sets
) {
	Size count( 0 );
	for ( auto const & working_stem : working_stems ) {
		vector1< Size > stem_pairing_set;
		for ( auto const & elem : working_stem ) {
			count++;
			BasePair base_pair( elem.first, elem.second, WATSON_CRICK, WATSON_CRICK, ANTIPARALLEL );
			runtime_assert( !rna_pairing_list.has_value( base_pair ) );
			rna_pairing_list.push_back( base_pair );
			stem_pairing_set.push_back( count );
		}
		working_stem_pairing_sets.push_back( stem_pairing_set );
	}

	for ( auto const & base_pair : working_obligate_pairs ) {
		count++;
		rna_pairing_list.push_back( base_pair );
		working_obligate_pair_sets.push_back( make_vector1( count ) );
	}
}

void
RNA_DeNovoProtocolMover::setup_ft_from_silent(
	bool const ft_from_silent_provided,
	std::string const & ft_silent_file,
	bool const tag_provided,
	utility::vector1< std::string > const & tag,
	bool & use_fold_tree_from_silent_file,
	core::kinematics::FoldTree & fold_tree_from_silent_file
) {
	if ( ft_from_silent_provided ) {
		// if this works, unify with the refine silent pose list stuff
		vector1<pose::PoseOP> refine_pose_list;
		if ( !ft_silent_file.empty() ) {

			core::import_pose::pose_stream::SilentFilePoseInputStreamOP input;

			if ( tag_provided ) {
				utility::vector1< std::string > filename_vector;
				filename_vector.push_back( ft_silent_file );
				input = utility::pointer::make_shared< core::import_pose::pose_stream::SilentFilePoseInputStream >( filename_vector, tag );
			} else {
				input = utility::pointer::make_shared< core::import_pose::pose_stream::SilentFilePoseInputStream >( ft_silent_file );
			}

			int pose_count = 1;
			input->set_order_by_energy( true );
			while ( input->has_another_pose() && pose_count < 2 ) {

				pose::PoseOP new_pose( new pose::Pose );
				input->fill_pose( *new_pose, *rsd_set_ );
				fold_tree_from_silent_file = new_pose->fold_tree();
				use_fold_tree_from_silent_file = true;
				pose_count += 1;

				//new_pose->dump_pdb( "test_pose_dump_" + ObjexxFCL::string_of( pose_count ) + ".pdb" );

			}
		}

	}
}

void
RNA_DeNovoProtocolMover::setup_final_res_lists(
	utility::vector1< core::Size > const & working_res,
	utility::vector1< core::Size > const & cutpoint_open_in_full_model,
	std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > > const & extra_minimize_res_rc,
	std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > > const & extra_minimize_chi_res_rc,
	std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > > const & output_jump_res_rc,
	bool const minimize_rna_specified,
	std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > > const & syn_chi_rc,
	std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > > const & anti_chi_rc,
	std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > > const & block_stack_above_rc,
	std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > > const & block_stack_below_rc,
	core::pose::full_model_info::FullModelParametersOP const & full_model_parameters // can use a const OP -- not a COP!
) {
	// Moved this up to restore refine_native functionality
	vector1< Size > extra_minimize_res =
		full_model_parameters->conventional_to_full( extra_minimize_res_rc );
	vector1< Size > extra_minimize_chi_res =
		full_model_parameters->conventional_to_full( extra_minimize_chi_res_rc );
	vector1< Size > output_jump_res =
		full_model_parameters->conventional_to_full( output_jump_res_rc );

	if ( !minimize_rna_specified && !minimize_rna_has_been_specified_ ) utility_exit_with_message( "Please specify either '-minimize_rna true' or '-minimize_rna false'." );

	if ( minimize_rna_has_been_specified_ ) {
		options_->set_minimize_structure( minimize_rna_ );
	}

	// runtime_assert( opts[ OptionKeys::score::include_neighbor_base_stacks ].user() ); // user should specify -include_neighbor_base_stacks true or -include_neighbor_base_stacks false.

	// some stuff to update in *options*
	options_->set_extra_minimize_res( working_res_map( extra_minimize_res, working_res ) );
	options_->set_extra_minimize_chi_res( working_res_map( extra_minimize_chi_res, working_res ) );
	options_->set_output_jump_res( working_res_map( output_jump_res, working_res ) );
	////////////////////
	// Step 19
	////////////////////
	using namespace core::pose::full_model_info;
	// could also set up other stuff inside full_model_parameters -- see protocols/stepise/FullModelInfoSetupFromCommandLine.cc
	// a better route, however, would be to *deprecate* this setup code, and instead use that stepwise code + build_full_model to
	// handle FARFAR setup. -- rhiju & amwatkins, dec. 2016.
	TR.Debug << "CUTPOINT OPEN IN FULL MODEL" << std::endl;
	TR.Debug << cutpoint_open_in_full_model << std::endl;
	full_model_parameters->set_parameter_as_res_list( CUTPOINT_OPEN, cutpoint_open_in_full_model );
	full_model_parameters->set_parameter_as_res_list( RNA_SYN_CHI,
		full_model_parameters->conventional_to_full( syn_chi_rc ) );
	full_model_parameters->set_parameter_as_res_list( RNA_ANTI_CHI,
		full_model_parameters->conventional_to_full( anti_chi_rc ) );
	full_model_parameters->set_parameter_as_res_list( RNA_BLOCK_STACK_ABOVE,
		full_model_parameters->conventional_to_full( block_stack_above_rc ) );
	full_model_parameters->set_parameter_as_res_list( RNA_BLOCK_STACK_BELOW,
		full_model_parameters->conventional_to_full( block_stack_below_rc ) );
	full_model_parameters->set_parameter_as_res_list( EXTRA_MINIMIZE, extra_minimize_res );
}

void
RNA_DeNovoProtocolMover::de_novo_setup_from_options( utility::options::OptionCollection const & opts ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::import_pose::options;
	using namespace core::chemical;
	using namespace core::id;
	using namespace core::pose;
	using namespace core::pose::rna;
	using namespace core::pose::full_model_info;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::scoring::func;

	////////////////////
	// Step 1
	////////////////////

	// Sequence setup:
	// FullModelParameters is a nice object that holds sequence, non-standard residues ("Z[Mg]"),
	// and what chains and residue numbers to use.
	FullModelParametersOP full_model_parameters;
	vector1< Size > cutpoint_open_in_full_model;
	initialize_sequence_information(
		opts[ OptionKeys::rna::denovo::offset ](),
		opts[ basic::options::OptionKeys::edensity::mapfile ].user(),
		opts[ basic::options::OptionKeys::rna::denovo::rna_protein_docking ](),
		opts[ OptionKeys::rna::denovo::virtual_anchor ].user(),
		opts[ OptionKeys::full_model::cutpoint_open ].user(),
		opts[ OptionKeys::full_model::cutpoint_open ].resnum_and_chain(),
		full_model_parameters, cutpoint_open_in_full_model );

	std::string const sequence = full_model_parameters->full_sequence();

	////////////////////
	// Step 3
	////////////////////
	// secondary structure setup.
	RNA_SecStruct secstruct( opts[ OptionKeys::rna::denovo::secstruct ](), opts[ OptionKeys::rna::denovo::secstruct_file ](), sequence );
	// "general" secondary structure includes non-canonical pairs that should be connected by jumps during run; used with -bps_moves.
	RNA_SecStruct secstruct_general( opts[ OptionKeys::rna::denovo::secstruct_general ](), opts[ OptionKeys::rna::denovo::secstruct_general_file ](), sequence );
	check_secstructs(
		opts[ basic::options::OptionKeys::edensity::mapfile ].user(),
		opts[ basic::options::OptionKeys::rna::denovo::rna_protein_docking ](),
		sequence, secstruct, secstruct_general );

	vector1< Size > remove_pair =
		full_model_parameters->conventional_to_full( opts[ OptionKeys::rna::denovo::remove_pair ].resnum_and_chain() );
	runtime_assert( remove_pair.size() % 2 == 0 );
	for ( Size n = 1; n <= remove_pair.size(); n += 2 ) {
		secstruct.remove_pair( std::make_pair( remove_pair[ n ], remove_pair[ n + 1 ] ) );
	}

	////////////////////
	// Step 5
	////////////////////
	// in full_model numbering.
	vector1< Size > input_res;
	vector1< vector1< int > > resnum_list;
	vector1< Size > input_res_initialize;
	vector1< vector1< Size > > helical_substruct_res;
	vector1< Size > dock_chunks_res;
	input_numbering_setup(
		opts[ OptionKeys::in::file::input_res ].resnum_and_chain(),
		opts[ OptionKeys::rna::denovo::input_silent_res ].resnum_and_chain(),
		sequence, full_model_parameters, input_res, resnum_list, input_res_initialize, helical_substruct_res, dock_chunks_res );

	////////////////////
	// Step 7
	////////////////////
	vector1< Size > obligate_pair =
		full_model_parameters->conventional_to_full( opts[ OptionKeys::rna::denovo::obligate_pair ].resnum_and_chain() );
	vector1< std::string > obligate_pair_explicit = opts[ OptionKeys::rna::denovo::obligate_pair_explicit ]();
	vector1< Size > domain_map( sequence.size(), 0 );
	setup_obligate_pair( opts[ OptionKeys::rna::denovo::minimize::extra_minimize_res ].resnum_and_chain(),
		full_model_parameters, resnum_list, secstruct, secstruct_general, cutpoint_open_in_full_model,
		obligate_pair, obligate_pair_explicit, domain_map );

	////////////////////
	// Step 8
	////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Figure out working sequence, secstruct, secstruct_general
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	vector1< Size > working_res        =
		full_model_parameters->conventional_to_full( opts[ full_model::working_res ].resnum_and_chain() ); //all working stuff
	std::sort( working_res.begin(), working_res.end() ); // some the following depends on correct order.
	if ( working_res.size() == 0 ) {
		for ( Size n = 1; n <= sequence.size(); n++ ) working_res.push_back( n );
	}

	// previously, kept spaces in sequences to mimic user-input with spaces.
	// don't worry about that anymore -- try to handle above in get_sequence_info.
	std::string working_sequence = working_res_map( sequence, working_res );
	RNA_SecStruct working_secstruct = working_res_map( secstruct, working_res );
	RNA_SecStruct working_secstruct_general = working_res_map( secstruct_general, working_res );

	TR << "Sequence:            " << working_sequence << std::endl;
	TR << "Secstruct:           " << working_secstruct.secstruct() << std::endl;
	if ( !secstruct_general.blank() ) TR << "Secstruct [general]: " << working_secstruct_general.secstruct() << std::endl;

	// Step 8B [good ol' secstruct_legacy]
	std::string secstruct_legacy = opts[ OptionKeys::rna::denovo::secstruct_legacy ]();
	std::string working_secstruct_legacy = working_res_map( secstruct_legacy, working_res );
	if ( secstruct_legacy.size() > 0  ) TR << "Secstruct [legacy]: " << working_secstruct_legacy << std::endl;

	////////////////////
	// Step 9
	////////////////////
	utility::vector1< core::pose::rna::BasePair > working_obligate_pairs;
	vector1< Size > const working_input_res = working_res_map( input_res, working_res );
	vector1< vector1< std::pair< Size, Size > > > working_stems = working_secstruct.stems();
	if ( options_->fixed_stems() ) {
		update_working_obligate_pairs_with_stems( working_obligate_pairs, working_stems, working_input_res );
	}
	///////////////////////////////////////////////////////////////////
	// initialize variables needed for RNA_DeNovoParams (rna_params_)
	///////////////////////////////////////////////////////////////////
	vector1< Size > working_cutpoint_open   = working_res_map( cutpoint_open_in_full_model, working_res, true /*leave out last residue*/ );
	for ( Size n = 1; n < working_res.size(); n++ ) {
		if ( working_res[ n+1 ] > working_res[n] + 1  && !working_cutpoint_open.has_value( n ) )  working_cutpoint_open.push_back( n );
	}

	vector1< Size > const working_cutpoint_closed = working_res_map( full_model_parameters->conventional_to_full( opts[ full_model::cutpoint_closed ].resnum_and_chain() ), working_res );
	vector1< Size > const working_cutpoint_cyclize = working_res_map( full_model_parameters->conventional_to_full( opts[ full_model::cyclize ].resnum_and_chain() ), working_res );
	vector1< Size > const working_twoprime = working_res_map( full_model_parameters->conventional_to_full( opts[ full_model::twoprime ].resnum_and_chain() ), working_res );
	vector1< Size > const working_fiveprime_cap = working_res_map( full_model_parameters->conventional_to_full( opts[ full_model::fiveprime_cap ].resnum_and_chain() ), working_res );
	vector1< Size > const working_block_stack_above_res = working_res_map( full_model_parameters->conventional_to_full( opts[ full_model::rna::block_stack_above_res ].resnum_and_chain() ), working_res );
	vector1< Size > const working_block_stack_below_res = working_res_map( full_model_parameters->conventional_to_full( opts[ full_model::rna::block_stack_below_res ].resnum_and_chain() ), working_res );
	vector1< Size > const working_virtual_anchor  = working_res_map( full_model_parameters->conventional_to_full( opts[ OptionKeys::rna::denovo::virtual_anchor ].resnum_and_chain() ), working_res );
	vector1< Size > const working_input_res_initialize = working_res_map( input_res_initialize, working_res );
	vector1< vector1< Size > > working_helical_substruct_res;
	for ( Size i=1; i<= helical_substruct_res.size(); ++i ) {
		working_helical_substruct_res.push_back( working_res_map( helical_substruct_res[i] , working_res) );
	}
	vector1< Size > const working_dock_chunks_res = working_res_map( dock_chunks_res, working_res );


	////////////////////
	// Step 10
	////////////////////
	Pose full_pose;
	bool is_rna_and_protein = false;
	initial_pose_setup(
		opts[ OptionKeys::in::path::path ]()[1],
		opts[ OptionKeys::rna::denovo::lores_scorefxn ].user(),
		opts[ OptionKeys::in::file::native ].user(),
		opts[ basic::options::OptionKeys::edensity::mapfile ].user(),
		opts[ OptionKeys::rna::denovo::working_native ].user(),
		opts[ OptionKeys::rna::denovo::working_native ](),
		full_model_parameters, working_res, full_pose, is_rna_and_protein );

	////////////////////
	// Step 11
	////////////////////
	constraint_setup(
		opts[ OptionKeys::constraints::cst_file ].user(),
		opts[ OptionKeys::constraints::cst_file ].user() ? opts[ OptionKeys::constraints::cst_file ](1) : "",
		working_res, full_pose );

	////////////////////
	// Step 12
	////////////////////
	{
		////////////////////////
		// working data
		////////////////////////
		if ( opts[ OptionKeys::rna::data_file].user() ) {
			core::io::rna::RNA_DataReader rna_data_reader( opts[ OptionKeys::in::path::path ]()[1] + opts[ OptionKeys::rna::data_file ]  );
			// note that this actually does look at conventional numbering in a smart way (but not chains yet):
			rna_data_reader.fill_rna_data_info( *pose_ );
		}
	}

	////////////////////
	// Step 13 + 14
	////////////////////
	////////////////////////////
	// working_obligate_pairs.
	// working_obligate_pairs [explicit]
	////////////////////////////
	refine_working_obligate_pairs(
		opts[ OptionKeys::rna::denovo::remove_obligate_pair ].resnum_and_chain(),
		full_model_parameters, obligate_pair, obligate_pair_explicit, working_res, working_obligate_pairs );

	////////////////////
	// Step 15
	////////////////////
	vector1< std::pair< vector1< Size >, vector1< Size > > > working_chain_connections;
	setup_working_chain_connections(
		opts[ OptionKeys::rna::denovo::chain_connection ](),
		full_model_parameters, working_res, working_chain_connections );


	////////////////////
	// Step 16B
	////////////////////
	RNA_BasePairList rna_pairing_list;
	utility::vector1 < utility::vector1 <core::Size > > working_obligate_pair_sets;
	utility::vector1 < utility::vector1 <core::Size > > working_stem_pairing_sets;
	setup_pairing_sets( working_stems, working_obligate_pairs, rna_pairing_list, working_obligate_pair_sets, working_stem_pairing_sets );

	/////////////////
	/// Get the fold tree from a specified silent file
	/// this is a bit of a weird hack for being able refine from an initial structure with density
	/// reading the pose from the silent file directly creates a big mess...
	////////////////
	bool use_fold_tree_from_silent_file = false;
	core::kinematics::FoldTree fold_tree_from_silent_file;
	setup_ft_from_silent(
		opts[ OptionKeys::rna::denovo::get_fold_tree_from_silent_file ].user(),
		opts[ OptionKeys::rna::denovo::get_fold_tree_from_silent_file ](),
		opts[ OptionKeys::rna::denovo::fold_tree_from_silent_file_tag ].user(),
		opts[ OptionKeys::rna::denovo::fold_tree_from_silent_file_tag ](),
		use_fold_tree_from_silent_file,
		fold_tree_from_silent_file
	);

	////////////////////
	// Step 17
	////////////////////
	{
		///////////////////////////////////
		// package above into params
		///////////////////////////////////
		rna_params_ = utility::pointer::make_shared< RNA_DeNovoParameters >();
		rna_params_->set_rna_pairing_list( rna_pairing_list );
		rna_params_->set_obligate_pairing_sets( working_obligate_pair_sets );
		rna_params_->set_stem_pairing_sets( working_stem_pairing_sets );
		rna_params_->set_chain_connections( working_chain_connections );
		rna_params_->set_cutpoints_open( working_cutpoint_open );
		rna_params_->set_cutpoints_closed( working_cutpoint_closed );
		rna_params_->set_cutpoints_cyclize( working_cutpoint_cyclize );
		rna_params_->set_twoprime( working_twoprime );
		rna_params_->set_fiveprime_cap( working_fiveprime_cap );
		rna_params_->set_block_stack_above_res( working_block_stack_above_res );
		rna_params_->set_block_stack_below_res( working_block_stack_below_res );
		rna_params_->set_virtual_anchor_attachment_points( working_virtual_anchor );
		rna_params_->set_rna_and_protein( is_rna_and_protein );
		rna_params_->set_rna_secstruct_legacy( working_secstruct_legacy );
		rna_params_->set_use_fold_tree_from_silent_file( use_fold_tree_from_silent_file );
		rna_params_->set_fold_tree_from_silent_file( fold_tree_from_silent_file );
	}

	////////////////////
	// Step 18
	////////////////////
	options_->set_input_res( working_input_res );
	options_->set_input_res_initialize( working_input_res_initialize );
	options_->set_helical_substruct_res( working_helical_substruct_res );
	options_->set_dock_chunks_res( working_dock_chunks_res );
	setup_final_res_lists(
		working_res,
		cutpoint_open_in_full_model,
		opts[ OptionKeys::rna::denovo::minimize::extra_minimize_res ].resnum_and_chain(),
		opts[ OptionKeys::rna::denovo::minimize::extra_minimize_chi_res ].resnum_and_chain(),
		opts[ OptionKeys::rna::denovo::out::output_jump_res ].resnum_and_chain(),
		opts[ OptionKeys::rna::denovo::minimize_rna ].user(),
		opts[ full_model::rna::force_syn_chi_res_list ].resnum_and_chain(),
		opts[ full_model::rna::force_anti_chi_res_list ].resnum_and_chain(),
		opts[ full_model::rna::block_stack_above_res ].resnum_and_chain(),
		opts[ full_model::rna::block_stack_below_res ].resnum_and_chain(),
		full_model_parameters
	);

	vector1< Size > dummy_domain_map( sequence.size(), 0 );
	full_model_parameters->set_parameter( INPUT_DOMAIN, domain_map /* domain_map */ );
	full_model_parameters->set_parameter( FIXED_DOMAIN, dummy_domain_map /*stepwise::setup::figure_out_fixed_domain_map( domain_map, extra_minimize_res ) */  );

	// Set up FullModelInfo (so we can use info stored there like SYN_CHI_RES)
	FullModelInfoOP full_model_info( new FullModelInfo( full_model_parameters ) );
	full_model_info->set_res_list( working_res );
	TR.Debug << working_res << std::endl;
	// need to include the virtual res if we're using a density map

	set_full_model_info( *pose_, full_model_info );
}

void
RNA_DeNovoProtocolMover::de_novo_setup_from_tag( utility::tag::TagCOP const & tag ) {
	using namespace core::import_pose::options;
	using namespace core::chemical;
	using namespace core::id;
	using namespace core::pose;
	using namespace core::pose::rna;
	using namespace core::pose::full_model_info;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::scoring::func;

	////////////////////
	// Step 1
	////////////////////

	// Sequence setup:
	// FullModelParameters is a nice object that holds sequence, non-standard residues ("Z[Mg]"),
	// and what chains and residue numbers to use.
	FullModelParametersOP full_model_parameters;
	vector1< Size > cutpoint_open_in_full_model;
	initialize_sequence_information(
		tag->getOption< core::Size >( "offset", 0 ),
		tag->hasOption( "edensity_mapfile" ),
		tag->getOption< bool >( "rna_protein_docking", false ),
		tag->hasOption( "virtual_anchor" ),
		tag->hasOption( "cutpoint_open" ),
		utility::get_resnum_and_chain( tag->getOption< std::string >( "cutpoint_open", "" ) ),
		full_model_parameters, cutpoint_open_in_full_model );

	std::string const sequence = full_model_parameters->full_sequence();

	////////////////////
	// Step 3
	////////////////////
	// secondary structure setup.
	RNA_SecStruct secstruct( tag->getOption< std::string >( "secstruct", "" ), tag->getOption< std::string >( "secstruct_file", "" ), sequence );
	// "general" secondary structure includes non-canonical pairs that should be connected by jumps during run; used with -bps_moves.
	RNA_SecStruct secstruct_general( tag->getOption< std::string >( "secstruct_general", "" ), tag->getOption< std::string >( "secstruct_general_file", "" ), sequence );
	check_secstructs(
		tag->hasOption( "edensity_mapfile" ),
		tag->getOption< bool >( "rna_protein_docking", false ),
		sequence, secstruct, secstruct_general );

	vector1< Size > remove_pair =
		full_model_parameters->conventional_to_full(
		utility::get_resnum_and_chain( tag->getOption< std::string >( "remove_pair", "" ) ) );
	runtime_assert( remove_pair.size() % 2 == 0 );
	for ( Size n = 1; n <= remove_pair.size(); n += 2 ) {
		secstruct.remove_pair( std::make_pair( remove_pair[ n ], remove_pair[ n + 1 ] ) );
	}

	////////////////////
	// Step 5
	////////////////////
	// in full_model numbering.
	vector1< Size > input_res;
	vector1< vector1< int > > resnum_list;
	vector1< Size > input_res_initialize;
	vector1< vector1< Size > > helical_substruct_res;
	vector1< Size > dock_chunks_res;
	input_numbering_setup(
		utility::get_resnum_and_chain( tag->getOption< std::string >( "input_res", "") ),
		utility::get_resnum_and_chain( tag->getOption< std::string >( "input_silent_res", "") ),
		sequence, full_model_parameters, input_res, resnum_list, input_res_initialize, helical_substruct_res, dock_chunks_res );

	////////////////////
	// Step 7
	////////////////////
	vector1< Size > obligate_pair =
		full_model_parameters->conventional_to_full(
		utility::get_resnum_and_chain( tag->getOption< std::string >( "obligate_pair", "") ) );

	// empty string => empty array with this function
	vector1< std::string > obligate_pair_explicit =
		utility::string_split_simple( tag->getOption< std::string >( "obligate_pair_explicit", "") );
	vector1< Size > domain_map( sequence.size(), 0 );
	setup_obligate_pair(
		utility::get_resnum_and_chain( tag->getOption< std::string >( "extra_minimize_res", "") ),
		full_model_parameters, resnum_list, secstruct, secstruct_general, cutpoint_open_in_full_model,
		obligate_pair, obligate_pair_explicit, domain_map );

	////////////////////
	// Step 8
	////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Figure out working sequence, secstruct, secstruct_general
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	vector1< Size > working_res  = full_model_parameters->conventional_to_full(
		utility::get_resnum_and_chain( tag->getOption< std::string >( "working_res", "") ) ); //all working stuff
	std::sort( working_res.begin(), working_res.end() ); // some the following depends on correct order.
	if ( working_res.size() == 0 ) {
		for ( Size n = 1; n <= sequence.size(); n++ ) working_res.push_back( n );
	}

	// previously, kept spaces in sequences to mimic user-input with spaces.
	// don't worry about that anymore -- try to handle above in get_sequence_info.
	std::string working_sequence = working_res_map( sequence, working_res );
	RNA_SecStruct working_secstruct = working_res_map( secstruct, working_res );
	RNA_SecStruct working_secstruct_general = working_res_map( secstruct_general, working_res );

	TR << "Sequence:            " << working_sequence << std::endl;
	TR << "Secstruct:           " << working_secstruct.secstruct() << std::endl;
	if ( !secstruct_general.blank() ) TR << "Secstruct [general]: " << working_secstruct_general.secstruct() << std::endl;

	// Step 8B [good ol' secstruct_legacy]
	std::string secstruct_legacy = tag->getOption< std::string >( "secstruct_legacy", "" );
	std::string working_secstruct_legacy = working_res_map( secstruct_legacy, working_res );
	if ( secstruct_legacy.size() > 0  ) TR << "Secstruct [legacy]: " << working_secstruct_legacy << std::endl;

	////////////////////
	// Step 9
	////////////////////
	utility::vector1< core::pose::rna::BasePair > working_obligate_pairs;
	vector1< Size > const working_input_res = working_res_map( input_res, working_res );
	vector1< vector1< std::pair< Size, Size > > > working_stems = working_secstruct.stems();
	if ( options_->fixed_stems() ) {
		update_working_obligate_pairs_with_stems( working_obligate_pairs, working_stems, working_input_res );
	}
	///////////////////////////////////////////////////////////////////
	// initialize variables needed for RNA_DeNovoParams (rna_params_)
	///////////////////////////////////////////////////////////////////
	vector1< Size > working_cutpoint_open   = working_res_map( cutpoint_open_in_full_model, working_res, true /*leave out last residue*/ );
	for ( Size n = 1; n < working_res.size(); n++ ) {
		if ( working_res[ n+1 ] > working_res[n] + 1  && !working_cutpoint_open.has_value( n ) )  working_cutpoint_open.push_back( n );
	}

	vector1< Size > const working_cutpoint_closed = working_res_map( full_model_parameters->conventional_to_full( utility::get_resnum_and_chain( tag->getOption< std::string >( "cutpoint_closed", "") ) ), working_res );
	vector1< Size > const working_cutpoint_cyclize = working_res_map( full_model_parameters->conventional_to_full( utility::get_resnum_and_chain( tag->getOption< std::string >( "cyclize_res", "") ) ), working_res );
	vector1< Size > const working_twoprime = working_res_map( full_model_parameters->conventional_to_full( utility::get_resnum_and_chain( tag->getOption< std::string >( "twoprime_res", "") ) ), working_res );
	vector1< Size > const working_fiveprime_cap = working_res_map( full_model_parameters->conventional_to_full( utility::get_resnum_and_chain( tag->getOption< std::string >( "fiveprime_cap_res", "") ) ), working_res );
	vector1< Size > const working_block_stack_above_res = working_res_map( full_model_parameters->conventional_to_full( utility::get_resnum_and_chain( tag->getOption< std::string >( "block_stack_above_res", "") ) ), working_res );
	vector1< Size > const working_block_stack_below_res = working_res_map( full_model_parameters->conventional_to_full( utility::get_resnum_and_chain( tag->getOption< std::string >( "block_stack_below_res", "") ) ), working_res );
	vector1< Size > const working_virtual_anchor  = working_res_map( full_model_parameters->conventional_to_full( utility::get_resnum_and_chain( tag->getOption< std::string >( "virtual_anchor", "") ) ), working_res );
	vector1< Size > const working_input_res_initialize = working_res_map( input_res_initialize, working_res );
	vector1< vector1< Size > > working_helical_substruct_res;
	for ( Size i=1; i<= helical_substruct_res.size(); ++i ) {
		working_helical_substruct_res.push_back( working_res_map( helical_substruct_res[i] , working_res) );
	}
	vector1< Size > const working_dock_chunks_res = working_res_map( dock_chunks_res, working_res );


	////////////////////
	// Step 10
	////////////////////
	Pose full_pose;
	bool is_rna_and_protein = false;
	initial_pose_setup(
		utility::string_split( tag->getOption< std::string >( "in_path", "" ) )[ 1 ],
		tag->hasOption( "lores_scorefxn" ),
		tag->hasOption( "native" ),
		tag->hasOption( "edensity_mapfile" ),
		tag->hasOption( "working_native" ),
		tag->getOption< std::string >( "working_native", "" ),
		full_model_parameters, working_res, full_pose, is_rna_and_protein );

	////////////////////
	// Step 11
	////////////////////
	constraint_setup(
		tag->hasOption( "cst_file" ),
		tag->getOption< std::string >( "cst_file", "" ),
		working_res, full_pose );

	////////////////////
	// Step 12
	////////////////////
	{
		////////////////////////
		// working data
		////////////////////////
		if ( tag->hasOption( "rna_data_file" ) ) {
			core::io::rna::RNA_DataReader rna_data_reader( tag->getOption< std::string >( "in_path", "" ) + tag->getOption< std::string >( "rna_data_file", "" ) );
			// note that this actually does look at conventional numbering in a smart way (but not chains yet):
			rna_data_reader.fill_rna_data_info( *pose_ );
		}
	}

	////////////////////
	// Step 13 + 14
	////////////////////
	////////////////////////////
	// working_obligate_pairs.
	// working_obligate_pairs [explicit]
	////////////////////////////
	refine_working_obligate_pairs(
		utility::get_resnum_and_chain( tag->getOption< std::string >( "remove_obligate_pair", "" ) ),
		full_model_parameters, obligate_pair, obligate_pair_explicit, working_res, working_obligate_pairs );

	////////////////////
	// Step 15
	////////////////////
	vector1< std::pair< vector1< Size >, vector1< Size > > > working_chain_connections;
	setup_working_chain_connections(
		// need an empty string => empty vector behavior
		utility::string_split_simple( tag->getOption< std::string >( "chain_connection", "" ) ),
		full_model_parameters, working_res, working_chain_connections );


	////////////////////
	// Step 16B
	////////////////////
	RNA_BasePairList rna_pairing_list;
	utility::vector1 < utility::vector1 <core::Size > > working_obligate_pair_sets;
	utility::vector1 < utility::vector1 <core::Size > > working_stem_pairing_sets;
	setup_pairing_sets( working_stems, working_obligate_pairs, rna_pairing_list, working_obligate_pair_sets, working_stem_pairing_sets );

	/////////////////
	/// Get the fold tree from a specified silent file
	/// this is a bit of a weird hack for being able refine from an initial structure with density
	/// reading the pose from the silent file directly creates a big mess...
	////////////////
	bool use_fold_tree_from_silent_file = false;
	core::kinematics::FoldTree fold_tree_from_silent_file;
	setup_ft_from_silent(
		tag->hasOption( "get_fold_tree_from_silent_file" ),
		tag->getOption< std::string >( "get_fold_tree_from_silent_file", "" ),
		tag->hasOption( "fold_tree_from_silent_file_tag" ),
		utility::string_split( tag->getOption< std::string >( "fold_tree_from_silent_file_tag", "" ) ),
		use_fold_tree_from_silent_file,
		fold_tree_from_silent_file
	);

	////////////////////
	// Step 17
	////////////////////
	{
		///////////////////////////////////
		// package above into params
		///////////////////////////////////
		rna_params_ = utility::pointer::make_shared< RNA_DeNovoParameters >();
		rna_params_->set_rna_pairing_list( rna_pairing_list );
		rna_params_->set_obligate_pairing_sets( working_obligate_pair_sets );
		rna_params_->set_stem_pairing_sets( working_stem_pairing_sets );
		rna_params_->set_chain_connections( working_chain_connections );
		rna_params_->set_cutpoints_open( working_cutpoint_open );
		rna_params_->set_cutpoints_closed( working_cutpoint_closed );
		rna_params_->set_cutpoints_cyclize( working_cutpoint_cyclize );
		rna_params_->set_twoprime( working_twoprime );
		rna_params_->set_fiveprime_cap( working_fiveprime_cap );
		rna_params_->set_block_stack_above_res( working_block_stack_above_res );
		rna_params_->set_block_stack_below_res( working_block_stack_below_res );
		rna_params_->set_virtual_anchor_attachment_points( working_virtual_anchor );
		rna_params_->set_rna_and_protein( is_rna_and_protein );
		rna_params_->set_rna_secstruct_legacy( working_secstruct_legacy );
		rna_params_->set_use_fold_tree_from_silent_file( use_fold_tree_from_silent_file );
		rna_params_->set_fold_tree_from_silent_file( fold_tree_from_silent_file );
	}

	////////////////////
	// Step 18
	////////////////////
	options_->set_input_res( working_input_res );
	options_->set_input_res_initialize( working_input_res_initialize );
	options_->set_helical_substruct_res( working_helical_substruct_res );
	options_->set_dock_chunks_res( working_dock_chunks_res );
	setup_final_res_lists(
		working_res,
		cutpoint_open_in_full_model,
		utility::get_resnum_and_chain( tag->getOption< std::string >( "extra_minimize_res", "" ) ),
		utility::get_resnum_and_chain( tag->getOption< std::string >( "extra_minimize_chi_res", "" ) ),
		utility::get_resnum_and_chain( tag->getOption< std::string >( "output_jump_res", "" ) ),
		tag->hasOption( "minimize_rna" ),
		utility::get_resnum_and_chain( tag->getOption< std::string >( "force_syn_chi_res_list", "" ) ),
		utility::get_resnum_and_chain( tag->getOption< std::string >( "force_anti_chi_res_list", "" ) ),
		utility::get_resnum_and_chain( tag->getOption< std::string >( "block_stack_above_res", "" ) ),
		utility::get_resnum_and_chain( tag->getOption< std::string >( "block_stack_below_res", "" ) ),
		full_model_parameters
	);
	//
	vector1< Size > dummy_domain_map( sequence.size(), 0 );
	full_model_parameters->set_parameter( INPUT_DOMAIN, domain_map /* domain_map */ );
	full_model_parameters->set_parameter( FIXED_DOMAIN, dummy_domain_map /*stepwise::setup::figure_out_fixed_domain_map( domain_map, extra_minimize_res ) */  );

	// Set up FullModelInfo (so we can use info stored there like SYN_CHI_RES)
	FullModelInfoOP full_model_info( new FullModelInfo( full_model_parameters ) );
	full_model_info->set_res_list( working_res );
	TR.Debug << working_res << std::endl;
	// need to include the virtual res if we're using a density map

	set_full_model_info( *pose_, full_model_info );
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
// Utility functions
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
vector1< Size >
RNA_DeNovoProtocolMover::working_res_map( vector1< Size > const & vec,
	vector1< Size > const & working_res,
	bool const leave_out_last_working_residue /* = false */ ) const
{
	if ( working_res.size() == 0 ) return vec;
	vector1< Size > working_vec;
	for ( Size const m : vec ) {
		if ( leave_out_last_working_residue && m == working_res[ working_res.size() ]  ) continue;
		if ( working_res.has_value( m ) ) {
			working_vec.push_back( working_res.index( m ) );
		}
	}
	return working_vec;
}

std::string
RNA_DeNovoProtocolMover::working_res_map( std::string const & seq_input,
	vector1< Size > const & working_res,
	bool const annotations_in_brackets /* = true */ ) const
{
	if ( seq_input.size() == 0 ) return "";
	std::string seq( seq_input );
	core::sequence::strip_spacers( seq, annotations_in_brackets );
	if ( working_res.size() == 0 ) return seq;
	std::string working_seq;
	for ( Size m = 1; m <= seq.size(); m++ ) {
		if ( working_res.has_value( m ) ) {
			working_seq += seq[ m - 1 ];
		}
	}
	return working_seq;
}

// Following not handling spacers correctly...
core::pose::rna::RNA_SecStruct
RNA_DeNovoProtocolMover::working_res_map( core::pose::rna::RNA_SecStruct const & rna_secstruct,
	vector1< Size > const & working_res ) const
{
	std::string working_secstruct = working_res_map( rna_secstruct.secstruct(), working_res, false /*annotations_in_brackets*/ );
	return core::pose::rna::RNA_SecStruct( working_secstruct );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Add by FCC: process PDBs and generate reasonable obligate pairs
void
RNA_DeNovoProtocolMover::get_seq_and_resnum( std::string const & pdb,
	std::string & seq,
	vector1< int > & resnum,
	vector1< char > & chain,
	vector1< std::string > & segid ) const
{
	using namespace core::pose;
	using namespace core::import_pose;
	PoseOP pose_op = pose_from_file( pdb );
	Pose & pose = *pose_op;
	PDBInfoOP pdb_info = pose.pdb_info();
	seq = pose.sequence();
	resnum.clear();
	chain.clear();
	for ( Size n = 1; n <= pose.size(); n++ ) {
		resnum.push_back( pdb_info->number( n ) );
		chain.push_back( pdb_info->chain( n ) );
		segid.push_back( pdb_info->segmentID( n ) );
	}
}

std::string
RNA_DeNovoProtocolMover::get_silent_seq( std::string const & silent_file ) const
{
	using namespace core::io::silent;
	SilentFileOptions opts;
	SilentFileData silent_file_data( opts );
	return silent_file_data.get_sequence( silent_file );
}

std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > >
RNA_DeNovoProtocolMover::get_silent_resnum( std::string const & silent_file ) const
{
	using namespace core::io::silent;
	SilentFileOptions opts;
	SilentFileData silent_file_data( opts );
	return silent_file_data.get_resnum( silent_file );
}

bool
RNA_DeNovoProtocolMover::already_listed_in_obligate_pair( vector1< Size > const & new_pair,
	vector1< Size > const & obligate_pair ) const
{
	Size const new_pos1 = new_pair[ 1 ];
	Size const new_pos2 = new_pair[ 2 ];
	for ( Size m = 0; m < obligate_pair.size()/2; m++ ) {
		Size pos1 = obligate_pair[ 2*m+1 ];
		Size pos2 = obligate_pair[ 2*m+2 ];
		if ( pos1 == new_pos1 && pos2 == new_pos2 ) {
			return true;
		}
	}
	return false;
}

bool
RNA_DeNovoProtocolMover::already_listed_in_obligate_pair( vector1< Size > const & new_pair,
	vector1< Size > const & obligate_pair,
	vector1< Size > const & obligate_pair_explicit ) const
{
	vector1< Size > all_pair = obligate_pair;
	//for ( Size n = 1; n <= obligate_pair.size(); n++ ) all_pair.push_back( obligate_pair[ n ] );
	for ( Size const exp : obligate_pair_explicit ) all_pair.push_back( exp );
	return already_listed_in_obligate_pair( new_pair, all_pair );
}

///////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocolMover::update_working_obligate_pairs_with_stems(
	vector1< pose::rna::BasePair > & working_obligate_pairs,
	vector1< vector1< std::pair< Size, Size > > > const & working_stems,
	vector1< Size > const & working_input_res ) const
{
	for ( auto const & working_stem : working_stems ) {
		bool stem_in_input_res( false );

		for ( auto const & pair : working_stem ) {
			if ( working_input_res.has_value( pair.first ) &&
					working_input_res.has_value( pair.second ) ) {
				stem_in_input_res = true; break;
			}
		}
		if ( stem_in_input_res ) continue;

		core::pose::rna::BasePair base_pair( working_stem[1].first, working_stem[1].second, WATSON_CRICK, WATSON_CRICK, ANTIPARALLEL );
		working_obligate_pairs.push_back( base_pair );
	}
}

std::string RNA_DeNovoProtocolMoverCreator::keyname() const {
	return RNA_DeNovoProtocolMover::mover_name();
}

protocols::moves::MoverOP
RNA_DeNovoProtocolMoverCreator::create_mover() const {
	return utility::pointer::make_shared< RNA_DeNovoProtocolMover >();
}

void RNA_DeNovoProtocolMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RNA_DeNovoProtocolMover::provide_xml_schema( xsd );
}



} //movers
} //denovo
} //rna
} //protocols
