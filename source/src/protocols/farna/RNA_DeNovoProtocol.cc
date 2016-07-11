// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RNA de novo fragment assembly
/// @brief protocols that are specific to RNA_DeNovoProtocol
/// @details
/// @author Rhiju Das, Parin Sripakdeevong, Fang-Chieh Chou


// Unit headers
#include <protocols/farna/RNA_DeNovoProtocol.hh>
#include <protocols/farna/options/RNA_DeNovoProtocolOptions.hh>
#include <protocols/farna/RNA_FragmentMonteCarlo.hh>
#include <protocols/farna/fragments/FullAtomRNA_Fragments.hh>
#include <protocols/farna/movers/RNA_LoopCloser.hh>
#include <protocols/farna/base_pairs/RNA_BasePairHandler.hh>
#include <protocols/farna/movers/RNA_Minimizer.hh>
#include <protocols/farna/setup/RNA_DeNovoParameters.hh>
#include <protocols/farna/movers/RNA_Relaxer.hh>
#include <protocols/farna/setup/RNA_DeNovoPoseInitializer.hh>
#include <protocols/farna/libraries/RNA_ChunkLibrary.hh>

// Package headers
#include <protocols/toolbox/AtomLevelDomainMap.hh>
#include <core/pose/rna/RNA_BasePairClassifier.hh>
#include <protocols/farna/util.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_VDW_BinChecker.hh>
#include <protocols/stepwise/modeler/rna/checker/VDW_CachedRepScreenInfo.hh>

// Project headers
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
#include <core/io/silent/util.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/pose/rna/RNA_BaseDoubletClasses.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

#include <utility/file/file_sys_util.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>

// External library headers

//C++ headers
#include <iostream>
#include <ctime>

//Auto Headers
#include <protocols/viewer/GraphicsState.hh>
#include <utility/vector1.hh>

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The original protocol for Fragment Assembly of RNA (FARNA), first developed in rosetta++ in 2006.
//
// Refactored in 2015 so that actual monte carlo sampling is encapsulated in *RNA_FragmentMonteCarlo*, along with minimize/relax --
//   this allows call of fragment assembly from within other protocols like stepwise.
//
// Setup of options has moved into RNA_DeNovoProtocolOptions and RNA_FragmentMonteCarloOptions.
//
// So the jobs remaining of RNA_DeNovoProtocol are simply:
//
//   1. Setup of various libraries and movers that might be used by RNA_DeNovoProtocol.
//   2. Compute metrics on final poses (RMSD, etc.)
//   3. Output to silent files
//
// These remaining functionalities could be subsumed into a JobDistributor with appropriate PoseMetrics.
//
//                       -- Rhiju, 2015
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using namespace ObjexxFCL::format; // AUTO USING NS
using namespace core;

namespace protocols {
namespace farna {

static THREAD_LOCAL basic::Tracer TR( "protocols.farna.RNA_DeNovoProtocol" );

RNA_DeNovoProtocol::RNA_DeNovoProtocol( RNA_DeNovoProtocolOptionsCOP options,
	RNA_DeNovoParametersCOP params):
	Mover(),
	options_( options ),
	rna_params_( params )
{
	if ( rna_params_ == 0 ) {
		if ( options_->rna_params_file().size() > 0 ) {
			rna_params_ = RNA_DeNovoParametersCOP( new RNA_DeNovoParameters( options_->rna_params_file() ) );
		} else {
			rna_params_ = RNA_DeNovoParametersCOP( new RNA_DeNovoParameters );
		}
	}
	Mover::type("RNA_DeNovoProtocol");
}

/// @brief Clone this object
protocols::moves::MoverOP RNA_DeNovoProtocol::clone() const {
	return protocols::moves::MoverOP( new RNA_DeNovoProtocol(*this) );
}

//////////////////////////////////////////////////
RNA_DeNovoProtocol::~RNA_DeNovoProtocol() {}

/// @details  Apply the RNA de novo modeling protocol to a pose.
///
void RNA_DeNovoProtocol::apply( core::pose::Pose & pose ) {

	using namespace core::pose;
	using namespace core::scoring;
	using namespace core::io::pdb;
	using namespace core::io::silent;
	using namespace protocols::farna;

	///////////////////////////////////////////////////////////////////////////
	// A bunch of initialization
	///////////////////////////////////////////////////////////////////////////
	if ( options_->dump_pdb() ) pose.dump_pdb( "init.pdb" );

	// Set up the cached vdw rep screen info in the pose
	// requires further setup later (once user input fragments have been inserted in the pose)
	if ( options_->filter_vdw() ) {
		protocols::stepwise::modeler::rna::checker::fill_vdw_cached_rep_screen_info_from_command_line( pose );
		vdw_grid_ = protocols::stepwise::modeler::rna::checker::RNA_VDW_BinCheckerOP( new protocols::stepwise::modeler::rna::checker::RNA_VDW_BinChecker( pose ) );
	}

	// RNA score function (both low-res and high-res).
	initialize_scorefxn( pose );

	// Some other silent file setup
	initialize_lores_silent_file();
	initialize_tag_is_done();

	RNA_DeNovoPoseInitializerOP rna_de_novo_pose_setup( new RNA_DeNovoPoseInitializer( *rna_params_ )  );
	rna_de_novo_pose_setup->set_bps_moves( options_->bps_moves() );
	rna_de_novo_pose_setup->set_root_at_first_rigid_body( options_->root_at_first_rigid_body() );
	bool refine_pose( refine_pose_list_.size() > 0 || options_->refine_pose() );
	if ( !refine_pose ) rna_de_novo_pose_setup->initialize_for_de_novo_protocol( pose, options_->ignore_secstruct() ); // virtualize phosphates, but no chainbreaks -- PUT HIGHER?

	//Keep a copy for resetting after each decoy.
	Pose start_pose = pose;
	// This copy of the pose does not have its grid_vdw stuff set up
	// But it will be set up in RNA_FragmentMonteCarlo before any sampling happens

	///////////////////////////////////////////////////////////////////////////
	// Main Loop.
	///////////////////////////////////////////////////////////////////////////
	Size refine_pose_id( 1 );
	std::list< core::Real > all_lores_score_final; // used for filtering.
	for ( Size n = 1; n <= options_->nstruct(); n++ ) {

		std::string const out_file_tag = "S_"+lead_zero_string_of( n, 6 );
		if ( tag_is_done_[ out_file_tag ] ) continue;
		std::clock_t start_time = clock();

		if ( refine_pose_list_.size() > 0 ) {
			pose = *refine_pose_list_[ refine_pose_id ];
			++refine_pose_id;
			if ( refine_pose_id > refine_pose_list_.size() ) refine_pose_id = 1;
		} else {
			pose = start_pose;
		}

		RNA_ChunkLibraryOP user_input_chunk_library( new RNA_ChunkLibrary( options_->chunk_pdb_files(), options_->chunk_silent_files(), pose,
			options_->input_res(), rna_params_->allow_insert_res() ) );
		RNA_BasePairHandlerOP rna_base_pair_handler( refine_pose ? new RNA_BasePairHandler( pose ) : new RNA_BasePairHandler( *rna_params_ ) );

		rna_fragment_monte_carlo_ = RNA_FragmentMonteCarloOP( new RNA_FragmentMonteCarlo( options_ ) );
		rna_fragment_monte_carlo_->set_out_file_tag( out_file_tag );
		rna_fragment_monte_carlo_->set_native_pose( get_native_pose() );
		rna_fragment_monte_carlo_->set_denovo_scorefxn( denovo_scorefxn_ );
		rna_fragment_monte_carlo_->set_hires_scorefxn( hires_scorefxn_ );
		rna_fragment_monte_carlo_->set_user_input_chunk_library( user_input_chunk_library );
		rna_fragment_monte_carlo_->set_rna_base_pair_handler( rna_base_pair_handler ); // could later have this look inside pose's sec_struct_info
		rna_fragment_monte_carlo_->set_refine_pose( refine_pose );
		// Set the vdw grid here
		if ( options_->filter_vdw() ) {
			rna_fragment_monte_carlo_->set_vdw_grid( vdw_grid_ );
		}
		if ( !refine_pose ) rna_fragment_monte_carlo_->set_rna_de_novo_pose_setup( rna_de_novo_pose_setup ); // only used for resetting fold-tree & cutpoints on each try.
		rna_fragment_monte_carlo_->set_all_lores_score_final( all_lores_score_final );  // used for filtering.

		rna_fragment_monte_carlo_->apply( pose );

		all_lores_score_final = rna_fragment_monte_carlo_->all_lores_score_final(); // might have been updated, user for filtering.
		if ( options_->output_lores_silent_file() ) align_and_output_to_silent_file( *(rna_fragment_monte_carlo_->lores_pose()), lores_silent_file_, out_file_tag );

		std::string const out_file_name = out_file_tag + ".pdb";
		if ( options_->dump_pdb() ) dump_pdb( pose,  out_file_name );

		if ( options_->save_times() ) setPoseExtraScore( pose, "time", static_cast< Real >( clock() - start_time ) / CLOCKS_PER_SEC );
		align_and_output_to_silent_file( pose, options_->silent_file(), out_file_tag );
	} //nstruct
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
std::string
RNA_DeNovoProtocol::get_name() const {
	return "RNA_DeNovoProtocol";
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocol::show(std::ostream & output) const
{
	Mover::show(output);
	output <<   "nstruct:                       " << options_->nstruct()  <<
		"\nUser defined MC cycles:        " << (options_->user_defined_cycles()  ? "True" : "False") <<
		"\nAll RNA fragment file:         " << options_->all_rna_fragments_file() <<
		"\nDump pdb:                      " << (options_->dump_pdb() ? "True" : "False") <<
		"\nMinimize structure:            " << (options_->minimize_structure() ? "True" : "False") <<
		"\nRelax structure:               " << (options_->relax_structure() ? "True" : "False") <<
		"\nIgnore secstruct:              " << (options_->ignore_secstruct() ? "True" : "False") <<
		"\nClose loops at end:            " << (options_->close_loops() ? "True" : "False") <<
		"\nClose loops in last round:     " << (options_->close_loops() ? "True" : "False") <<
		"\nClose loops after each move:   " << (options_->close_loops_after_each_move() ? "True" : "False") <<
		"\nSimple rmsd cutoff relax:      " << (options_->simple_rmsd_cutoff_relax() ? "True" : "False") <<
		"\nAllow bulges:                  " << (options_->allow_bulge() ? "True" : "False") <<
		"\nAllow consecutive bulges:      " << (options_->allow_consecutive_bulges() ? "True" : "False") <<
		"\nUse chem shift data:           " << (options_->use_chem_shift_data() ? "True" : "False") <<
		"\nDefault temperature for MC:    " << options_->temperature() <<
		"\nInput rna params file?:        " << ((options_->rna_params_file() == "" ) ? "No" : "Yes") <<
		"\nJump library file:             " << options_->jump_library_file() <<
		"\nOutput lores silent file:      " << (options_->output_lores_silent_file() ? "True" : "False") <<
		"\nFilter lores base pairs:       " << (options_->filter_lores_base_pairs() ? "True" : "False") <<
		"\nFilter lores base pairs early: " << (options_->filter_lores_base_pairs_early() ? "True" : "False") <<
		"\nFilter chain closure:          " << (options_->filter_chain_closure() ? "True" : "False") <<
		"\nFilter chain closure distance: " << options_->filter_chain_closure_distance() <<
		"\nFilter chain closure halfway:  " << (options_->filter_chain_closure_halfway() ? "True" : "False") <<
		"\nVary bond geometry:            " << (options_->vary_bond_geometry() ? "True" : "False") <<
		"\nBinary RNA output:             " << (options_->binary_rna_output() ? "True" : "False") <<
		"\nStaged constraints:            " << (options_->staged_constraints() ? "True" : "False") <<
		"\nTitrate stack bonus:           " << (options_->titrate_stack_bonus() ? "True" : "False") <<
		"\nMove first rigid body:         " << (options_->move_first_rigid_body() ? "True" : "False") <<
		"\nRoot at first rigid body:      " << (options_->root_at_first_rigid_body() ? "True" : "False") <<
		"\nOutput Filters:                " << (options_->output_filters() ? "True" : "False") <<
		"\nAutofilter:                    " << (options_->autofilter() ? "True" : "False") <<
		"\nAutofilter score quantile:     " << options_->autofilter_score_quantile() <<
		"\nBase pair step moves:          " << (options_->bps_moves() ? "True" : "False");

}

///////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocol::initialize_scorefxn( core::pose::Pose & pose ) {

	using namespace core::scoring;
	using namespace basic::options;

	// RNA low-resolution score function.
	denovo_scorefxn_ = ScoreFunctionFactory::create_score_function( options_->lores_scorefxn() );
	if ( scoring::rna::nonconst_rna_scoring_info_from_pose( pose ).rna_data_info().rna_reactivities().size() > 0 ) {
		denovo_scorefxn_->set_weight( core::scoring::rna_chem_map_lores, option[ OptionKeys::score::rna_chem_map_lores_weight ]() );
	}

	if ( options_->filter_vdw() ) {
		denovo_scorefxn_->set_weight( core::scoring::grid_vdw, options_->grid_vdw_weight() );
	}

	// initial_denovo_scorefxn_ = denovo_scorefxn_->clone();
	if ( options_->chainbreak_weight() > -1.0 ) denovo_scorefxn_->set_weight( chainbreak, options_->chainbreak_weight() );
	if ( options_->linear_chainbreak_weight() > -1.0 ) denovo_scorefxn_->set_weight( linear_chainbreak, options_->linear_chainbreak_weight() );

	initialize_constraints( pose );

	// RNA high-resolution score function.
	hires_scorefxn_ = get_rna_hires_scorefxn();//->clone();

}

///////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocol::initialize_constraints( core::pose::Pose & pose ) {

	using namespace core::scoring;

	if ( pose.constraint_set()->has_constraints() ) {
		denovo_scorefxn_->set_weight( atom_pair_constraint, 1.0 );
		denovo_scorefxn_->set_weight( coordinate_constraint, 1.0 ); // now useable in RNA denovo!
		denovo_scorefxn_->set_weight( base_pair_constraint, 1.0 );
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocol::initialize_tag_is_done()
{
	tag_is_done_ = core::io::silent::initialize_tag_is_done( options_->silent_file() );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocol::initialize_lores_silent_file() {

	if ( !options_->output_lores_silent_file() ) return;

	static std::string const new_prefix( "_LORES.out" );

	std::string::size_type pos = options_->silent_file().find( ".out", 0 );
	if ( pos == std::string::npos ) {
		utility_exit_with_message(  "If you want to output a lores silent file, better use .out suffix ==> " + options_->silent_file() );
	}
	lores_silent_file_ = options_->silent_file();
	lores_silent_file_.replace( pos, new_prefix.length(), new_prefix );
}


//////////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocol::calc_rmsds( core::io::silent::SilentStruct & s, core::pose::Pose & pose,
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

///////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocol::output_silent_struct(
	core::io::silent::SilentStruct & s, core::io::silent::SilentFileData & silent_file_data,
	std::string const & silent_file, pose::Pose & pose, std::string const & out_file_tag,
	bool const score_only /* = false */ ) const
{

	using namespace core::io::silent;
	using namespace core::scoring;

	if ( get_native_pose() ) calc_rmsds( s, pose, out_file_tag );

	// TR << "ADD_NUMBER_BASE_PAIRS" << std::endl;
	add_number_base_pairs( pose, s );
	// TR << "ADD_NUMBER_NATIVE_BASE_PAIRS" << std::endl;
	if ( get_native_pose() ) add_number_native_base_pairs( pose, s );

	// hopefully these will end up in silent file...
	if ( options_->output_filters() && ( rna_fragment_monte_carlo_ != 0 ) ) {
		s.add_energy(  "lores_early", rna_fragment_monte_carlo_->lores_score_early() );
		if ( options_->minimize_structure() ) s.add_energy( "lores_final", rna_fragment_monte_carlo_->lores_score_final() );
	}

	TR << "Outputting to silent file: " << silent_file << std::endl;
	silent_file_data.write_silent_struct( s, silent_file, score_only );

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocol::output_to_silent_file(
	core::pose::Pose & pose,
	std::string const & silent_file,
	std::string const & out_file_tag,
	bool const score_only /* = false */ ) const
{

	using namespace core::io::silent;
	using namespace core::scoring;

	// Silent file setup?
	//static SilentFileData silent_file_data;
	SilentFileData silent_file_data;

	// What is all this rigamarole, making the silent struct data?
	// Why do I need to supply the damn file name? That seems silly.
	TR << "Making silent struct for " << out_file_tag << std::endl;

	SilentStructOP s = ( options_->binary_rna_output() ) ? SilentStructOP( new BinarySilentStruct( pose, out_file_tag ) ) :
		SilentStructOP( new RNA_SilentStruct(   pose, out_file_tag ) );

	if ( options_->use_chem_shift_data() ) add_chem_shift_info( *s, pose);

	output_silent_struct( *s, silent_file_data, silent_file, pose, out_file_tag, score_only );

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocol::align_and_output_to_silent_file( core::pose::Pose & pose, std::string const & silent_file, std::string const & out_file_tag ) const
{
	rna_fragment_monte_carlo_->align_pose( pose, true /*verbose*/ );
	output_to_silent_file( pose, silent_file, out_file_tag, false /*score_only*/ );
}


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
// Following may better fit in a util.cc ,or pose_metrics...
void
RNA_DeNovoProtocol::add_number_base_pairs( pose::Pose const & pose, io::silent::SilentStruct & s ) const
{
	using namespace scoring::rna;
	using namespace pose::rna;
	using namespace conformation;
	using namespace core::chemical::rna;

	utility::vector1< core::pose::rna::BasePair > base_pair_list;
	utility::vector1< bool > is_bulged;
	core::pose::rna::classify_base_pairs( pose, base_pair_list, is_bulged );

	Size N_WC( 0 ), N_NWC( 0 );

	for ( Size n = 1; n <= base_pair_list.size(); n++ ) {

		BasePair const base_pair = base_pair_list[ n ];

		Size const i = base_pair.res1();
		Size const j = base_pair.res2();

		BaseEdge const k = base_pair.edge1();
		BaseEdge const m = base_pair.edge2();

		Residue const & rsd_i( pose.residue( i ) );
		Residue const & rsd_j( pose.residue( j ) );

		if ( ( k == WATSON_CRICK && m == WATSON_CRICK && base_pair.orientation() == ANTIPARALLEL )  &&
				core::chemical::rna::possibly_canonical( rsd_i.aa(), rsd_j.aa() ) )  {
			N_WC++;
		} else {
			N_NWC++;
		}
	}

	s.add_string_value( "N_WC",  ObjexxFCL::format::I( 9, N_WC) );
	s.add_string_value( "N_NWC", ObjexxFCL::format::I( 9, N_NWC ) );
	s.add_string_value( "N_BS",  ObjexxFCL::format::I( 9, core::pose::rna::get_number_base_stacks( pose ) ) );

}

/////////////////////////////////////////////////////////////////////
bool
check_in_base_pair_list( pose::rna::BasePair const & base_pair /*from native*/,
	utility::vector1< core::pose::rna::BasePair > const & base_pair_list /*for decoy*/)
{
	using namespace pose::rna;

	bool in_list( false );

	for ( Size n = 1; n <= base_pair_list.size(); n++ ) {

		BasePair const base_pair2 = base_pair_list[ n ];

		if ( base_pair == base_pair2 ) {
			in_list = true;
			break;
		}

		if ( base_pair.flipped() == base_pair2 ) {
			in_list = true;
			break;
		}

	}

	return in_list;

}

/////////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocol::add_number_native_base_pairs(pose::Pose & pose, io::silent::SilentStruct & s ) const
{
	if ( !get_native_pose() ) return;

	using namespace scoring::rna;
	using namespace chemical::rna;
	using namespace conformation;

	pose::Pose native_pose = *get_native_pose();

	utility::vector1< core::pose::rna::BasePair > base_pair_list;
	utility::vector1< bool > is_bulged;
	core::pose::rna::classify_base_pairs( pose, base_pair_list, is_bulged );

	utility::vector1< core::pose::rna::BasePair > base_pair_list_native;
	utility::vector1< bool > is_bulged_native;
	core::pose::rna::classify_base_pairs( native_pose, base_pair_list_native, is_bulged_native );

	Size N_WC_NATIVE( 0 ), N_NWC_NATIVE( 0 );
	Size N_WC( 0 ), N_NWC( 0 );

	for ( Size n = 1; n <= base_pair_list_native.size(); n++ ) {

		//Real const score = it->first;
		//  Real const SCORE_CUTOFF( -1.0 );
		//  if (score > SCORE_CUTOFF) continue;

		core::pose::rna::BasePair const base_pair = base_pair_list_native[ n ];

		Size const i = base_pair.res1();
		Size const j = base_pair.res2();

		BaseEdge const k = base_pair.edge1();
		BaseEdge const m = base_pair.edge2();

		Residue const & rsd_i( pose.residue( i ) );
		Residue const & rsd_j( pose.residue( j ) );

		//std::cout << " NATIVE BASE PAIR " << i << " " << j << " " << k << " " << m << " " << it->first << std::endl;

		if ( ( k == WATSON_CRICK && m == WATSON_CRICK && base_pair.orientation() == ANTIPARALLEL )  &&
				possibly_canonical( rsd_i.aa(), rsd_j.aa() ) )  {
			N_WC_NATIVE++;
			if ( check_in_base_pair_list( base_pair /*from native*/, base_pair_list /*for decoy*/) ) N_WC++;
		} else {
			N_NWC_NATIVE++;
			if ( check_in_base_pair_list( base_pair /*from native*/, base_pair_list /*for decoy*/) ) {
				N_NWC++;
			} else {
				std::cout << "Missing native base pair " << pose.residue( i ).name1() << i << " " << pose.residue(j).name1() << j << "  " << get_edge_from_num( k ) << " " << get_edge_from_num( m ) << " " << std::endl;
			}
		}
	}

	Real f_natWC( 0.0 ), f_natNWC( 0.0 ), f_natBP( 0.0 );
	if ( N_WC_NATIVE > 0 ) f_natWC = ( N_WC / (1.0 * N_WC_NATIVE) );
	if ( N_NWC_NATIVE > 0 ) f_natNWC = ( N_NWC / (1.0 * N_NWC_NATIVE) );
	if ( (N_WC_NATIVE + N_NWC_NATIVE) > 0 ) f_natBP = ( (N_WC+N_NWC) / (1.0 * (N_WC_NATIVE + N_NWC_NATIVE) ));

	s.add_energy( "f_natWC" , f_natWC );
	s.add_energy( "f_natNWC", f_natNWC );
	s.add_energy( "f_natBP" , f_natBP );

}

void
RNA_DeNovoProtocol::add_chem_shift_info(core::io::silent::SilentStruct & silent_struct,
	core::pose::Pose const & const_pose) const {

	using namespace core::scoring;
	using namespace core::pose;

	runtime_assert( options_->use_chem_shift_data() );

	pose::Pose chem_shift_pose=const_pose; //HARD COPY SLOW!

	core::scoring::ScoreFunctionOP temp_scorefxn( new ScoreFunction );

	temp_scorefxn->set_weight( scoring::rna_chem_shift  , 1.00 );

	(*temp_scorefxn)(chem_shift_pose);

	EnergyMap const & energy_map=chem_shift_pose.energies().total_energies();

	Real const rosetta_chem_shift_score= energy_map[ scoring::rna_chem_shift ];

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

std::ostream & operator<< ( std::ostream &os, RNA_DeNovoProtocol const & mover )
{
	mover.show(os);
	return os;
}

} //farna
} //protocols
