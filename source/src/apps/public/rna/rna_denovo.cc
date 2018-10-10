// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Rhiju Das

#include <utility/json_utilities.hh>

#if defined(ZEROMQ)  and  defined(_NLOHMANN_JSON_ENABLED_)
#include <protocols/network/hal.hh>
#include <protocols/network/util.hh>
#include <json.hpp>
#include <protocols/network/ui_mover.hh>
#include <core/import_pose/import_pose.hh>
#endif

// libRosetta headers
#include <core/types.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <protocols/viewer/viewers.hh>
#include <core/pose/Pose.hh>
#include <devel/init.hh>
#include <utility/vector1.hh>

//RNA stuff.
#include <protocols/rna/denovo/RNA_DeNovoProtocol.hh>
#include <core/import_pose/RNA_DeNovoSetup.hh>
#include <core/import_pose/options/RNA_DeNovoProtocolOptions.hh>
#include <protocols/rna/denovo/RNA_FragmentMonteCarlo.hh>
#include <core/fragment/rna/FullAtomRNA_Fragments.hh>
#include <protocols/rna/movers/RNA_LoopCloser.hh>
#include <core/import_pose/RNA_BasePairHandler.hh>
#include <protocols/rna/denovo/movers/RNA_Minimizer.hh>
#include <core/import_pose/RNA_DeNovoParameters.hh>
#include <protocols/rna/denovo/movers/RNA_Relaxer.hh>
#include <protocols/rna/denovo/RNA_DeNovoPoseInitializer.hh>
#include <core/import_pose/libraries/RNA_ChunkLibrary.hh>
#include <protocols/rna/setup/RNA_MonteCarloJobDistributor.hh>
#include <protocols/rna/setup/RNA_CSA_JobDistributor.hh>
#include <core/pose/rna/RNA_BasePairClassifier.hh>
#include <protocols/rna/denovo/util.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_VDW_BinChecker.hh>
#include <protocols/scoring/VDW_CachedRepScreenInfo.hh>


#include <core/conformation/Residue.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/chemical/rna/util.hh>
#include <core/scoring/rna/chemical_shift/RNA_ChemicalShiftPotential.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/io/silent/RNA_SilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/util.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/pose/rna/RNA_BaseDoubletClasses.hh>
#include <core/pose/rna/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>




#include <protocols/jd3/chunk_library/ChunkLibraryJobQueen.hh>
#include <protocols/jd3/chunk_library/MoverAndChunkLibraryJob.hh>
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/JobDistributor.hh>
#include <protocols/jd3/JobDistributorFactory.hh>
#include <protocols/jd3/JobQueen.hh>
#include <core/scoring/Energies.hh>

//#include <protocols/moves/PyMOLMover.hh>

// C++ headers
#include <iostream>
#include <string>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

#include <core/pose/annotated_sequence.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/sequence/Sequence.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/excn/Exceptions.hh>

#include <utility/options/keys/OptionKeyList.hh>


using namespace core;
using namespace pose;
using namespace core::io::silent;
using namespace protocols;
using namespace protocols::rna::denovo;
using namespace core::import_pose;
using namespace core::import_pose::libraries;
using namespace protocols::rna::setup;
using namespace basic::options::OptionKeys;
using namespace core::import_pose;
using namespace core::import_pose::options;
using namespace protocols::jd3;
using namespace protocols::jd3::chunk_library;
using utility::vector1;

#if defined(ZEROMQ)  and  defined(_NLOHMANN_JSON_ENABLED_)
using namespace utility;
using namespace protocols::network;
#endif

static basic::Tracer TR( "apps.public.rna.denovo" );

OPT_1GRP_KEY( Boolean, denovo, use_legacy_job_distributor )

class RNA_DeNovoJobQueen : public ChunkLibraryJobQueen {

	// Known issues (times to use legacy JD):
	// 1. Screws up secstruct input when there is no argument to -s.
	// 2. Screws up cutpoint setup, maybe other issues (edensity scoring? when there is
	// -new_fold_tree_initializer and density.
	// 3. Screws up setup from an explicit params file (to be deprecated in 2018 anyway).
	// 4. RNA-protein: screws up scoring function modification -- I think some option
	// is getting dropped at some point, maybe is_rna_and_protein().


public:

	RNA_DeNovoJobQueen() :
		ChunkLibraryJobQueen()
	{
		using namespace basic::options;
		utility::options::OptionKeyList opts;
		RNA_DeNovoSetup::list_options_read( opts );
		RNA_DeNovoProtocolOptions::list_options_read( opts );
		add_options( opts );
		add_option( OptionKeys::rna::denovo::output_res_num );
	}

	virtual
	JobOP
	complete_larval_job_maturation(
		protocols::jd3::LarvalJobCOP larval_job,
		utility::options::OptionCollectionCOP job_options,
		utility::vector1< JobResultCOP > const & // input_job_results
	) {


		MoverAndChunkLibraryJobOP mature_job( new MoverAndChunkLibraryJob );

		core::pose::PoseOP poseop = pose_for_job( larval_job, *job_options );
		mature_job->pose( poseop );


		RNA_DeNovoSetupOP rna_de_novo_setup( new RNA_DeNovoSetup );
		// AMW TODO: The below calls initialize_inputs_from_options, meaning however hard you try
		// the chunk library input source isn't the "only" source of input. (Not necessarily bad.)
		//rna_de_novo_setup->initialize_from_options( *job_options );
		rna_de_novo_setup->initialize_from_options( basic::options::option );

		Pose & pose = *poseop;

		// First step: extract protocol stuff back out. It does very little now save for app setup.

		//protocols::rna::denovo::RNA_DeNovoProtocol rna_de_novo_protocol( rna_de_novo_setup->options(),
		// rna_de_novo_setup->rna_params() );
		rna_params_ = rna_de_novo_setup->rna_params();
		options_ = rna_de_novo_setup->options();

		if ( rna_params_ == nullptr ) {
			if ( !options_->rna_params_file().empty() ) {
				std::cout << std::endl << options_->rna_params_file() << std::endl;
				rna_params_ = RNA_DeNovoParametersCOP( new RNA_DeNovoParameters( options_->rna_params_file() ) );
			} else {
				rna_params_ = RNA_DeNovoParametersCOP( new RNA_DeNovoParameters );
			}
		}
		native_pose_ = rna_de_novo_setup->native_pose();
		refine_pose_list_ = rna_de_novo_setup->refine_pose_list();

#if defined(ZEROMQ)  and  defined(_NLOHMANN_JSON_ENABLED_)
		protocols::network::AddUIObserver( pose );
#endif
		protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 600, 600 );
		// protocols::moves::AddPyMOLObserver( pose, false, 0.01);





		///////////////////////////////////////////////////////////////////////////
		// A bunch of initialization
		///////////////////////////////////////////////////////////////////////////
		if ( options_->dump_pdb() ) pose.dump_pdb( "init.pdb" );

		// Set up the cached vdw rep screen info in the pose
		// requires further setup later (once user input fragments have been inserted in the pose)
		if ( options_->filter_vdw() ) {
			protocols::scoring::fill_vdw_cached_rep_screen_info_from_options( pose, *job_options );
			vdw_grid_ = protocols::stepwise::modeler::rna::checker::RNA_VDW_BinCheckerOP( new protocols::stepwise::modeler::rna::checker::RNA_VDW_BinChecker( pose ) );
		}

		// RNA score function (both low-res and high-res).
		initialize_scorefxn( pose );

		// Some other silent file setup
		if ( options_->overwrite() ) remove_silent_file_if_it_exists( options_->silent_file() );
		std::string lores_silent_file_;
		if ( options_->output_lores_silent_file() ) {
			lores_silent_file_ = core::io::silent::get_outfile_name_with_tag( options_->silent_file(), "_LORES" );
		}
		auto tag_is_done_ = core::io::silent::initialize_tag_is_done( options_->silent_file() );

		RNA_DeNovoPoseInitializerOP rna_de_novo_pose_initializer( new RNA_DeNovoPoseInitializer( *rna_params_ )  );
		rna_de_novo_pose_initializer->set_bps_moves( options_->bps_moves() );
		rna_de_novo_pose_initializer->set_root_at_first_rigid_body( options_->root_at_first_rigid_body() );
		rna_de_novo_pose_initializer->set_dock_each_chunk( options_->dock_each_chunk() );
		rna_de_novo_pose_initializer->set_dock_each_chunk_per_chain( options_->dock_each_chunk_per_chain() );
		rna_de_novo_pose_initializer->set_center_jumps_in_single_stranded( options_->center_jumps_in_single_stranded() );
		rna_de_novo_pose_initializer->set_new_fold_tree_initializer( options_->new_fold_tree_initializer() );
		rna_de_novo_pose_initializer->set_model_with_density( options_->model_with_density() );
		bool refine_pose( refine_pose_list_.size() > 0 || options_->refine_pose() );
		if ( !refine_pose ) rna_de_novo_pose_initializer->initialize_for_de_novo_protocol( pose, false /*options_->ignore_secstruct()*/ ); // virtualize phosphates, but no chainbreaks -- PUT HIGHER?

		//Keep a copy for resetting after each decoy.
		Pose start_pose = pose;
		// This copy of the pose does not have its grid_vdw stuff set up
		// But it will be set up in RNA_FragmentMonteCarlo before any sampling happens

		///////////////////////////////////////////////////////////////////////////
		// Main Loop.
		///////////////////////////////////////////////////////////////////////////
		Size refine_pose_id( 1 );

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace basic::options::OptionKeys::stepwise::monte_carlo;

		std::string const silent_file = option[ out::file::silent ]();

		RNA_ChunkLibraryOP user_input_chunk_library( new RNA_ChunkLibrary( options_->chunk_pdb_files(), options_->chunk_silent_files(), pose,
			options_->input_res(), rna_params_->allow_insert_res() ) );
		RNA_BasePairHandlerOP rna_base_pair_handler( refine_pose ? new RNA_BasePairHandler( pose ) : new RNA_BasePairHandler( *rna_params_ ) );

		// main loop
		if ( rna_fragment_monte_carlo_ == nullptr ) {
			rna_fragment_monte_carlo_ = RNA_FragmentMonteCarloOP( new RNA_FragmentMonteCarlo( options_ ) );
		}

		rna_fragment_monte_carlo_->set_user_input_chunk_library( user_input_chunk_library );
		rna_fragment_monte_carlo_->set_rna_base_pair_handler( rna_base_pair_handler ); // could later have this look inside pose's sec_struct_info

		if ( options_->initial_structures_provided() ) {
			utility::vector1< std::string > silent_files_empty;
			RNA_ChunkLibraryOP user_input_chunk_initialization_library(
				new RNA_ChunkLibrary( options_->chunk_initialization_pdb_files(),
				silent_files_empty, pose, options_->input_res_initialize(), rna_params_->allow_insert_res() ) );
			rna_fragment_monte_carlo_->set_user_input_chunk_initialization_library( user_input_chunk_initialization_library );
		}

		protocols::rna::denovo::output::RNA_FragmentMonteCarloOutputterOP outputter;

		if ( refine_pose_list_.size() > 0 ) {
			pose = *refine_pose_list_[ refine_pose_id ];
			++refine_pose_id;
		} else {
			pose = start_pose;
		}

		rna_fragment_monte_carlo_->set_native_pose( get_native_pose() );
		rna_fragment_monte_carlo_->set_denovo_scorefxn( denovo_scorefxn_ );
		rna_fragment_monte_carlo_->set_hires_scorefxn( hires_scorefxn_ );
		rna_fragment_monte_carlo_->set_refine_pose( refine_pose );
		rna_fragment_monte_carlo_->set_is_rna_and_protein( rna_params_->is_rna_and_protein() ); // need to know this for the high resolution stuff
		if ( options_->filter_vdw() ) rna_fragment_monte_carlo_->set_vdw_grid( vdw_grid_ );
		if ( !refine_pose ) rna_fragment_monte_carlo_->set_rna_de_novo_pose_initializer( rna_de_novo_pose_initializer ); // only used for resetting fold-tree & cutpoints on each try.
		if ( outputter != 0 ) rna_fragment_monte_carlo_->set_outputter( outputter); // accumulate stats in histogram.

		mature_job->mover( rna_fragment_monte_carlo_ );

		return mature_job;
	}





	void initialize_scorefxn( core::pose::Pose & pose ) {

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

		initialize_constraints( pose );

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

	void initialize_constraints( core::pose::Pose & pose ) {

		using namespace core::scoring;

		if ( pose.constraint_set()->has_constraints() ) {
			denovo_scorefxn_->set_weight( atom_pair_constraint, 1.0 );
			denovo_scorefxn_->set_weight( coordinate_constraint, 1.0 ); // now useable in RNA denovo!
			denovo_scorefxn_->set_weight( base_pair_constraint, 1.0 );
		}
		if ( options_->rmsd_screen() > 0.0 && !denovo_scorefxn_->has_nonzero_weight( coordinate_constraint) ) {
			denovo_scorefxn_->set_weight( coordinate_constraint, 1.0 ); // now useable in RNA denovo!
		}
	}

	void calc_rmsds( core::io::silent::SilentStruct & s, core::pose::Pose & pose,
		std::string const & out_file_tag ) const
	{
		Real const rmsd = rna_fragment_monte_carlo_->get_rmsd_no_superimpose( pose );
		TR << "All atom rmsd: " << rmsd << " for " << out_file_tag << std::endl;
		s.add_energy( "rms",  rmsd );

		Real const rmsd_stems = rna_fragment_monte_carlo_->get_rmsd_stems_no_superimpose( pose );
		if ( rmsd_stems > 0.0 ) {
			TR << "All atom rmsd over stems: " << rmsd_stems << " for " << out_file_tag << std::endl;
			s.add_energy( "rms_stem", rmsd_stems );
		}
	}

	void output_silent_struct(
		core::io::silent::SilentStruct & s, core::io::silent::SilentFileData & silent_file_data,
		std::string const & silent_file, pose::Pose & pose, std::string const & out_file_tag,
		bool const score_only /* = false */ ) const
	{
		using namespace core::io::silent;
		using namespace core::scoring;
		using namespace core::pose::rna;

		if ( get_native_pose() ) calc_rmsds( s, pose, out_file_tag );

		add_number_base_pairs( pose, s );
		if ( get_native_pose() ) {
			add_number_native_base_pairs( pose, *get_native_pose(), s );
		}

		// hopefully these will end up in silent file...
		if ( options_->output_filters() && ( rna_fragment_monte_carlo_ != nullptr ) ) {
			s.add_energy(  "lores_early", rna_fragment_monte_carlo_->lores_score_early() );
			if ( options_->minimize_structure() ) s.add_energy( "lores_final", rna_fragment_monte_carlo_->lores_score_final() );
		}

		TR << "Outputting to silent file: " << silent_file << std::endl;
		silent_file_data.write_silent_struct( s, silent_file, score_only );
	}

	void output_to_silent_file(
		core::pose::Pose & pose,
		std::string const & silent_file,
		std::string const & out_file_tag,
		bool const score_only /* = false */ ) const
	{

		using namespace core::io::silent;
		using namespace core::scoring;

		if ( rna_params_->is_rna_and_protein() && !options_->minimize_structure() ) {
			// convert back to full atom (should give stupid coords... ok for now b/c protein doesn't move)
			// if protein sidechains have moved, then the pose should already be in full atom by now (?!)
			// if the structure is getting minimized, then it should already be converted back to full atom
			core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD, false /* no sloppy match */, true /* only switch protein residues */, true /* keep energies! */ );
			// but as soon as I score again, it tries to recalculate rnp scores (so they get set to 0)
		}

		// Silent file setup?
		SilentFileOptions opts;
		SilentFileData silent_file_data( opts );

		// What is all this rigamarole, making the silent struct data?
		// Why do I need to supply the damn file name? That seems silly.
		TR << "Making silent struct for " << out_file_tag << std::endl;

		SilentStructOP s = ( options_->binary_rna_output() ) ? SilentStructOP( new BinarySilentStruct( opts, pose, out_file_tag ) ) :
			SilentStructOP( new RNA_SilentStruct( opts, pose, out_file_tag ) );

		if ( options_->use_chem_shift_data() ) add_chem_shift_info( *s, pose);

		output_silent_struct( *s, silent_file_data, silent_file, pose, out_file_tag, score_only );
	}


	void align_and_output_to_silent_file( core::pose::Pose & pose, std::string const & silent_file, std::string const & out_file_tag ) const
	{
		rna_fragment_monte_carlo_->align_pose( pose, true /*verbose*/ );
		if ( rna_params_->is_rna_and_protein() && !options_->minimize_structure() ) {
			// copy the pose because we'll switch the residue type set when we output to silent file
			// this is a stupid hack b/c right now the pose has both centroid and full atom residues!
			core::pose::Pose pose_copy =  pose;
			output_to_silent_file( pose_copy, silent_file, out_file_tag, false /*score_only*/ );
		} else {
			output_to_silent_file( pose, silent_file, out_file_tag, false /*score_only*/ );
		}
	}

	void add_chem_shift_info(core::io::silent::SilentStruct & silent_struct,
		core::pose::Pose const & const_pose) const {

		using namespace core::scoring;
		using namespace core::pose;

		runtime_assert( options_->use_chem_shift_data() );

		pose::Pose chem_shift_pose = const_pose; //HARD COPY SLOW!
		core::scoring::ScoreFunctionOP temp_scorefxn( new ScoreFunction );
		temp_scorefxn->set_weight( rna_chem_shift  , 1.00 );
		(*temp_scorefxn)(chem_shift_pose);

		EnergyMap const & energy_map=chem_shift_pose.energies().total_energies();

		Real const rosetta_chem_shift_score= energy_map[ rna_chem_shift ];

		//This statement should be very fast except possibly the 1st call.
		core::scoring::rna::chemical_shift::RNA_ChemicalShiftPotential const &
			rna_chemical_shift_potential( core::scoring::ScoringManager::
			get_instance()->get_RNA_ChemicalShiftPotential() );

		Size const num_chem_shift_data_points=rna_chemical_shift_potential.get_total_exp_chemical_shift_data_points();

		//rosetta_chem_shift_score --> Sum_square chemical_shift deviation.

		Real const chem_shift_RMSD=sqrt( rosetta_chem_shift_score /
			float(num_chem_shift_data_points) );

		silent_struct.add_energy( "chem_shift_RMSD", chem_shift_RMSD);

		silent_struct.add_energy( "num_chem_shift_data",
			float(num_chem_shift_data_points) );

		if ( silent_struct.has_energy("rna_chem_shift")==false ) {
			//If missing this term, then the rna_chem_shift weight is probably
			//zero in the weight_file.
			silent_struct.add_energy( "rna_chem_shift", 0.0);
		}
	}

private:
	core::import_pose::options::RNA_DeNovoProtocolOptionsCOP options_;
	core::import_pose::RNA_DeNovoParametersCOP rna_params_;
	protocols::stepwise::modeler::rna::checker::RNA_VDW_BinCheckerOP vdw_grid_;

	std::string lores_silent_file_;

	RNA_FragmentMonteCarloOP rna_fragment_monte_carlo_;

	std::map< std::string, bool > tag_is_done_;

	core::scoring::ScoreFunctionOP denovo_scorefxn_;
	core::scoring::ScoreFunctionOP hires_scorefxn_;

	utility::vector1<core::pose::PoseOP> refine_pose_list_;
	//core::Size refine_pose_id_ = 1;

	pose::PoseCOP native_pose_;

	pose::PoseCOP get_native_pose() const { return native_pose_; }

	//std::list< core::Real > all_lores_score_final_; // used for filtering.

};



///////////////////////////////////////////////////////////////////////////////
#if defined(ZEROMQ)  and  defined(_NLOHMANN_JSON_ENABLED_)
core::pose::Pose
#else
void
#endif
rna_denovo()
{
	using namespace basic::options;
	if ( option[ OptionKeys::stepwise::superimpose_over_all ].user() ) {
		std::cout << "The use of -superimpose_over_all is deprecated. The behavior in question now defaults to TRUE and is turned off by providing a particular residue that is part of an anchoring input domain as -alignment_anchor_res." << std::endl;
	}

	using namespace core::pose;
	using namespace protocols::rna::denovo;


	protocols::jd3::JobDistributorOP jd = protocols::jd3::JobDistributorFactory::create_job_distributor();
	protocols::jd3::JobQueenOP queen( new RNA_DeNovoJobQueen );
	jd->go( queen );
#if defined(ZEROMQ)  and  defined(_NLOHMANN_JSON_ENABLED_)
	return Pose();
#endif
}

///////////////////////////////////////////////////////////////////////////////
#if defined(ZEROMQ)  and  defined(_NLOHMANN_JSON_ENABLED_)
core::pose::Pose
#else
void
#endif
rna_denovo_legacy()
{
	using namespace basic::options;
	if ( option[ OptionKeys::stepwise::superimpose_over_all ].user() ) {
		std::cout << "The use of -superimpose_over_all is deprecated. The behavior in question now defaults to TRUE and is turned off by providing a particular residue that is part of an anchoring input domain as -alignment_anchor_res." << std::endl;
	}

	using namespace core::pose;
	using namespace protocols::rna::denovo;

	RNA_DeNovoSetupOP rna_de_novo_setup( new RNA_DeNovoSetup );
	rna_de_novo_setup->initialize_inputs_from_options( basic::options::option );
	rna_de_novo_setup->initialize_from_command_line();
	Pose & pose = *( rna_de_novo_setup->pose() );

	protocols::rna::denovo::RNA_DeNovoProtocol rna_de_novo_protocol( rna_de_novo_setup->options(),
		rna_de_novo_setup->rna_params() );
	rna_de_novo_protocol.set_native_pose( rna_de_novo_setup->native_pose() );
	rna_de_novo_protocol.set_refine_pose_list( rna_de_novo_setup->refine_pose_list() );

#if defined(ZEROMQ)  and  defined(_NLOHMANN_JSON_ENABLED_)
	protocols::network::AddUIObserver( pose );
#endif

	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 600, 600 );
	// protocols::moves::AddPyMOLObserver( pose, false, 0.01);

	rna_de_novo_protocol.apply( pose );
#if defined(ZEROMQ)  and  defined(_NLOHMANN_JSON_ENABLED_)
	return pose;
#endif
}

#if defined(ZEROMQ)  and  defined(_NLOHMANN_JSON_ENABLED_)
///////////////////////////////////////////////////////////////
core::pose::Pose
my_main()
{
	using namespace basic::options;

	if ( option[ OptionKeys::denovo::use_legacy_job_distributor ].value() ) {
		return rna_denovo_legacy();
	} else {
		return rna_denovo();
	}
	protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}

#else
///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	using namespace basic::options;

	if ( option[ OptionKeys::denovo::use_legacy_job_distributor ].value() ) {
		rna_denovo_legacy();
	} else {
		rna_denovo();
	}
	protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}
#endif


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;

		NEW_OPT( denovo::use_legacy_job_distributor, "Use the legacy RNA_MonteCarloJobDistributor", false );

		std::cout << std::endl << "Basic usage:  " << argv[0] << "  -fasta <fasta file with sequence>  [ -native <native pdb file> ] " << std::endl;
		std::cout << std::endl << " Type -help for full slate of options." << std::endl << std::endl;

		option.add_relevant( in::file::fasta );
		option.add_relevant( in::file::native );
		option.add_relevant( in::file::input_res );
		option.add_relevant( out::file::silent );
		option.add_relevant( out::nstruct );
		option.add_relevant( out::overwrite );
		option.add_relevant( score::weights );
		option.add_relevant( score::set_weights );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::sequence );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::secstruct );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::secstruct_file );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::minimize_rna );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::rounds );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::fixed_stems );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::obligate_pair );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::obligate_pair_explicit );
		// option.add_relevant( basic::options::OptionKeys::rna::denovo::relax_rna );
		// option.add_relevant( basic::options::OptionKeys::rna::denovo::simple_relax );
		// option.add_relevant( basic::options::OptionKeys::rna::denovo::ignore_secstruct );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::lores_scorefxn );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::set_lores_weights);
		option.add_relevant( basic::options::OptionKeys::rna::denovo::cycles );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::temperature );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::jump_change_frequency );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::close_loops );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::close_loops_after_each_move );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::heat );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::staged_constraints );
		// option.add_relevant( basic::options::OptionKeys::rna::denovo::jump_library_file );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::params_file );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::filter_lores_base_pairs );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::filter_lores_base_pairs_early );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::filter_chain_closure );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::filter_chain_closure_halfway );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::filter_chain_closure_distance );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::output_filters );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::autofilter );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::no_filters );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::vall_torsions );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::use_1jj2_torsions );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::rna_lores_chainbreak_weight );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::rna_lores_linear_chainbreak_weight );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::allow_bulge  );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::allowed_bulge_res );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::allow_consecutive_bulges );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::move_first_rigid_body );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::root_at_first_rigid_body );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::suppress_bp_constraint );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::output_res_num );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::offset );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::tag );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::refine_silent_file );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::refine_native );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::bps_moves );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::out::dump );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::out::output_lores_silent_file );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::out::binary_output );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::minimize::minimizer_use_coordinate_constraints );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::minimize::extra_minimize_res );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::minimize::extra_minimize_chi_res );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::minimize::minimize_bps );
		option.add_relevant( full_model::cyclize );
		option.add_relevant( full_model::cutpoint_closed );
		option.add_relevant( full_model::cutpoint_open );
		option.add_relevant( full_model::rna::block_stack_above_res );
		option.add_relevant( full_model::rna::block_stack_below_res );
		option.add_relevant( basic::options::OptionKeys::rna::vary_geometry );
		option.add_relevant( basic::options::OptionKeys::rna::data_file );

		option.add_relevant( constraints::cst_file );

		// see if we can get density scoring to work
		option.add_relevant( basic::options::OptionKeys::edensity::mapfile );

		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////
		devel::init(argc, argv);

		option[ OptionKeys::chemical::patch_selectors ].push_back( "VIRTUAL_BASE" ); // for chemical mapping.
		option[ OptionKeys::chemical::patch_selectors ].push_back( "VIRTUAL_SIDE_CHAIN" ); // for proteins for vdw screen

		////////////////////////////////////////////////////////////////////////////
		// end of setup
		////////////////////////////////////////////////////////////////////////////
#if defined(ZEROMQ)  and  defined(_NLOHMANN_JSON_ENABLED_)
		{ // creating dummy pose object to trigger database load so later we can create Pose immeditaly
			core::pose::Pose p;
			core::import_pose::pose_from_pdbstring(p, "ATOM     17  N   ILE A   1      16.327  47.509  23.466  1.00  0.00\n");
		}

		//protocols::network::hal(specification, hal_executioner, protocols::network::CommandLineArguments{argc, argv} );

		hal({
			{"main", {my_main, {
			}
			} } },
			protocols::network::CommandLineArguments{argc, argv} );

		return 0;


#else
		protocols::viewer::viewer_main( my_main );
#endif

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

