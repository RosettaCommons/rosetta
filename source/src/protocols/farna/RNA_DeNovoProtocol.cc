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
#include <protocols/farna/RNA_DeNovoProtocolOptions.hh>
#include <protocols/farna/RNA_FragmentMonteCarlo.hh>
#include <protocols/farna/FullAtomRNA_Fragments.hh>
#include <protocols/farna/RNA_LoopCloser.hh>
#include <protocols/farna/RNA_Minimizer.hh>
#include <protocols/farna/RNA_Relaxer.hh>
#include <protocols/farna/RNA_StructureParameters.hh>
#include <protocols/farna/RNA_ChunkLibrary.hh>

// Package headers
#include <protocols/toolbox/AllowInsert.hh>
#include <core/pose/rna/RNA_BasePairClassifier.hh>
#include <protocols/stepwise/modeler/align/util.hh> //move this to toolbox/
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/farna/util.hh>

// Project headers
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
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
#include <core/io/pdb/pose_io.hh>
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

// External library headers

//C++ headers
#include <vector>
#include <list>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#ifdef WIN32
#include <ctime>
#endif

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

static THREAD_LOCAL basic::Tracer TR( "protocols.rna.RNA_DeNovoProtocol" );

RNA_DeNovoProtocol::RNA_DeNovoProtocol( RNA_DeNovoProtocolOptionsCOP options ):
	Mover(),
	options_( options )
{
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

	// RNA score function (both low-res and high-res).
	initialize_scorefxn( pose );

	//Keep a copy for resetting after each decoy.
	Pose start_pose = pose;

	// Some other silent file setup
	initialize_lores_silent_file();
	initialize_tag_is_done();

	///////////////////////////////////////////////////////////////////////////
	// Main Loop.
	///////////////////////////////////////////////////////////////////////////
	Size refine_pose_id( 1 );
	std::list< core::Real > all_lores_score_final; // used for filtering.
	for ( Size n = 1; n <= options_->nstruct(); n++ ) {

		std::string const out_file_tag = "S_"+lead_zero_string_of( n, 6 );
		if ( tag_is_done_[ out_file_tag ] ) continue;

		if ( refine_pose_list_.size() > 0 ) {
			pose = *refine_pose_list_[ refine_pose_id ];
			++refine_pose_id;
			if ( refine_pose_id > refine_pose_list_.size() ) refine_pose_id = 1;
		} else {
			pose = start_pose;
		}

		rna_fragment_monte_carlo_ = RNA_FragmentMonteCarloOP( new RNA_FragmentMonteCarlo( options_ ) );
		rna_fragment_monte_carlo_->set_out_file_tag( out_file_tag );
		rna_fragment_monte_carlo_->set_native_pose( get_native_pose() );
		rna_fragment_monte_carlo_->set_denovo_scorefxn( denovo_scorefxn_ );
		rna_fragment_monte_carlo_->set_hires_scorefxn( hires_scorefxn_ );
		rna_fragment_monte_carlo_->set_all_lores_score_final( all_lores_score_final );
		rna_fragment_monte_carlo_->set_refine_pose( refine_pose_list_.size() > 0 || options_->refine_pose() );

		rna_fragment_monte_carlo_->apply( pose );

		all_lores_score_final = rna_fragment_monte_carlo_->all_lores_score_final(); // might have been updated.
		if ( options_->output_lores_silent_file() ) align_and_output_to_silent_file( *(rna_fragment_monte_carlo_->lores_pose()), lores_silent_file_, out_file_tag );

		std::string const out_file_name = out_file_tag + ".pdb";
		if ( options_->dump_pdb() ) dump_pdb( pose,  out_file_name );

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
RNA_DeNovoProtocol::check_for_loop_modeling_case( std::map< core::id::AtomID, core::id::AtomID > & atom_id_map, pose::Pose const & /*pose*/ ) const
{
	// special case -- we only care about the loop(s). Pose has already been aligned to fixed residues.
	// this will be decided in align_and_output_to_silent_file.
	if ( rna_fragment_monte_carlo_ != 0 && rna_fragment_monte_carlo_->rna_chunk_library()->single_user_input_chunk() ) {
		std::map< core::id::AtomID, core::id::AtomID > loop_atom_id_map;
		TR << "In loop modeling mode, since there is a single user-inputted pose" << std::endl;
		for ( std::map< core::id::AtomID, core::id::AtomID >::const_iterator it = atom_id_map.begin(); it != atom_id_map.end(); it++ ) {
			Size domain( rna_fragment_monte_carlo_->rna_chunk_library()->allow_insert()->get_domain( it->second ) );
			if ( domain == 0 || domain == ROSETTA_LIBRARY_DOMAIN ) {
				loop_atom_id_map[ it->first ] = it->second;
				// TR << TR.Cyan << "Loop atom: " << atom_id_to_named_atom_id( it->second, pose ) << TR.Reset << std::endl;
			}
		}
		atom_id_map = loop_atom_id_map;
	}
}

//////////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocol::calc_rmsds( core::io::silent::SilentStruct & s, core::pose::Pose & pose,
	std::string const & out_file_tag ) const
{
	using namespace core::scoring;

	std::map< core::id::AtomID, core::id::AtomID > atom_id_map;
	setup_matching_heavy_atoms( *get_native_pose(), pose, atom_id_map ); // no virtuals, no hydrogens.
	check_for_loop_modeling_case( atom_id_map, pose );

	Real const rmsd = rms_at_corresponding_atoms_no_super( *get_native_pose(), pose, atom_id_map );
	TR << "All atom rmsd: " << rmsd << " for " << out_file_tag << std::endl;
	s.add_energy( "rms", rmsd );

	Real rmsd_stems = 0.0;
	std::list< Size > stem_residues( rna_fragment_monte_carlo_->rna_structure_parameters()->get_stem_residues( pose ) );

	if ( !stem_residues.empty() ) { //size() > 0 ) {
		rmsd_stems = all_atom_rmsd( *get_native_pose(), pose, stem_residues );
		TR << "All atom rmsd over stems: " << rmsd_stems << " for " << out_file_tag << std::endl;
	}
	s.add_energy( "rms_stem", rmsd_stems );

}

///////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoProtocol::output_silent_struct(
	core::io::silent::SilentStruct & s, core::io::silent::SilentFileData & silent_file_data,
	std::string const & silent_file, pose::Pose & pose, std::string const out_file_tag,
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

	bool loop_modeling_into_single_structure( false );

	// if input pdbs were specified with -s or -silent, then automatic alignment to first of these input chunks.
	// otherwise, align to native pose, if specified.
	if ( options_->input_res().size() > 0 ) {
		loop_modeling_into_single_structure = rna_fragment_monte_carlo_->rna_chunk_library()->superimpose_to_single_user_input_chunk( pose );
	}

	if ( !loop_modeling_into_single_structure && get_native_pose() ) {

		Pose const & native_pose = *get_native_pose();

		//realign to native for ease of viewing.
		// check for any fixed domains.
		utility::vector1< Size > superimpose_res; // = get_moving_res( pose, rna_fragment_monte_carlo_->rna_structure_parameters()->allow_insert() );

		// if no fixed domains, just superimpose over all residues.
		if ( superimpose_res.size() == 0 ) {
			for ( Size n = 1; n <= pose.total_residue(); n++ )  superimpose_res.push_back( n );
		}

		id::AtomID_Map< id::AtomID > const & alignment_atom_id_map_native =
			protocols::stepwise::modeler::align::create_alignment_id_map_legacy( pose, native_pose, superimpose_res ); // perhaps this should move to toolbox.

		TR << "Aligning pose to native." << std::endl;

		//pose.dump_pdb( "before_align.pdb");
		//  native_pose.dump_pdb( "native.pdb" );
		core::scoring::superimpose_pose( pose, native_pose, alignment_atom_id_map_native );
		//  pose.dump_pdb( "after_align.pdb");

	}

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
