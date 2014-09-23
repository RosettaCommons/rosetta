// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RNA de novo fragment assembly
/// @brief protocols that are specific to CoarseRNA_DeNovoProtocol
/// @detailed
/// @author Rhiju Das


// Unit headers
#include <protocols/coarse_rna/CoarseRNA_DeNovoProtocol.hh>

// Package headers
#include <protocols/coarse_rna/CoarseRNA_Fragments.hh>

// Project headers
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/toolbox/AllowInsert.hh>
#include <protocols/farna/MultipleDomainMover.hh>
#include <protocols/farna/RNA_ChunkLibrary.hh>
#include <protocols/farna/RNA_FragmentMover.hh>
#include <protocols/farna/RNA_StructureParameters.hh>
#include <protocols/farna/util.hh>
#include <core/io/rna/RNA_DataReader.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/coarse_rna/CoarseRNA_LoopCloser.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/pose/Pose.hh>
#include <basic/database/open.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/ShortestPathInFoldTree.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

#include <utility/file/file_sys_util.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>

#include <numeric/random/random.hh>
// AUTO-REMOVED #include <numeric/conversions.hh>

// External library headers

//C++ headers
#include <vector>
#include <list>
#include <string>
#include <sstream>
#include <fstream>
#ifdef WIN32
#include <ctime>
#endif

//Auto Headers
#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>


using namespace ObjexxFCL::format; // AUTO USING NS

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//  This was originally meant to be separate from RNA_DeNovoProtocol, but its
//  basically converging -- might be a good idea to unify the two classes in the
//  near future -- Rhiju, March 2010.
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

using namespace core;
using basic::T;

namespace protocols {
namespace coarse_rna {


static thread_local basic::Tracer TR( "protocols.coarse_rna.coarse_rna_denovo_protocol" );

CoarseRNA_DeNovoProtocol::CoarseRNA_DeNovoProtocol(
	 Size const nstruct,
	 Size const monte_carlo_cycles,
	 std::string const silent_file ):
    Mover(),
		nstruct_( nstruct ),
		monte_carlo_cycles_( monte_carlo_cycles ),
		rounds_( 10  ),
		silent_file_( silent_file ),
		freeze_domains_( false ),
		dump_pdb_( false ),
		domain_move_frequency_( 0.0 ),
		m_Temperature_( 5.0 ),
		sim_anneal_( true ),
		staged_constraints_( false ),
		frag_size_( 0 ),
		rna_params_file_( "" ),
		rna_data_file_( "" ),
		all_rna_fragments_file_( basic::database::full_name("1jj2_coarse_coords.txt") ),
		jump_library_file_( basic::database::full_name("chemical/rna/1jj2_coarse_jumps.dat" ) ),
		lores_scorefxn_( "farna/coarse_rna.wts" ),
		rna_structure_parameters_( protocols::farna::RNA_StructureParametersOP( new protocols::farna::RNA_StructureParameters ) ),
		rna_loop_closer_( protocols::coarse_rna::CoarseRNA_LoopCloserOP( new protocols::coarse_rna::CoarseRNA_LoopCloser ) ),
		close_loops_( false ),
		choose_best_solution_( false ),
		force_ideal_chainbreak_( false ),
		add_base_pair_constraints_( true ),
		view_monte_carlo_( false )
{
	Mover::type("CoarseRNA_DeNovoProtocol");
}

/// @brief Clone this object
protocols::moves::MoverOP CoarseRNA_DeNovoProtocol::clone() const {
	return protocols::moves::MoverOP( new CoarseRNA_DeNovoProtocol(*this) );
}

//////////////////////////////////////////////////
CoarseRNA_DeNovoProtocol::~CoarseRNA_DeNovoProtocol()
{
}

/// @details  Apply the RNA de novo modeling protocol to a pose.
///
void CoarseRNA_DeNovoProtocol::apply( core::pose::Pose & pose	) {

	using namespace core::pose;
	using namespace core::scoring;
	using namespace core::io::pdb;
	using namespace core::io::silent;

	///////////////////////////////////////////////////////////////////////////
	// A bunch of initialization
	///////////////////////////////////////////////////////////////////////////
	denovo_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( lores_scorefxn_ );

	initialize_tag_is_done();

	rna_structure_parameters_->initialize( pose, rna_params_file_, jump_library_file_, true /*ignore_secstruct*/ );

	rna_data_reader_ = core::io::rna::RNA_DataReaderOP( new core::io::rna::RNA_DataReader( rna_data_file_ ) );
	rna_data_reader_->fill_rna_data_info( pose ); // this seems repeated below? get rid of one instance?

	if( input_res_.size() > 0 )	rna_chunk_library_ = protocols::farna::RNA_ChunkLibraryOP( new protocols::farna::RNA_ChunkLibrary( chunk_silent_files_, pose, input_res_ ) );
	rna_structure_parameters_->set_allow_insert( rna_chunk_library_->allow_insert() );
	rna_structure_parameters_->setup_fold_tree_and_jumps_and_variants( pose );

	if (add_base_pair_constraints_) rna_structure_parameters_->setup_base_pair_constraints( pose );

	rna_chunk_library_->initialize_random_chunks( pose, dump_pdb_ );

	if ( dump_pdb_) std::cout << "Allow insert: " << std::endl;	rna_structure_parameters_->allow_insert()->show();

	protocols::farna::RNA_FragmentsOP rna_fragments( new CoarseRNA_Fragments( all_rna_fragments_file_ ) );
	frag_mover_ = protocols::farna::RNA_FragmentMoverOP( new protocols::farna::RNA_FragmentMover( rna_fragments, rna_structure_parameters_->allow_insert() ) );

	rna_data_reader_->fill_rna_data_info( pose );
	initialize_constraints( pose );

	// "force_ideal_chainbreak" means you can't change angles and dists at cutpoints.
	// This will also disallow torsion angle changes in OVL1, OVL2, OVU1 virtual atoms during fragment closure...
	rna_chunk_library_->allow_insert()->set_force_ideal_chainbreak( force_ideal_chainbreak_ );
	rna_loop_closer_->set_allow_insert( rna_structure_parameters_->allow_insert() );
	if ( choose_best_solution_ ) rna_loop_closer_->choose_best_solution_based_on_score_function( denovo_scorefxn_ );

	multiple_domain_mover_ = protocols::farna::MultipleDomainMoverOP( new protocols::farna::MultipleDomainMover( pose, rna_loop_closer_ ) );
	domain_move_frequency_ = 0.0;
	domain_move_frequency_ =  ( multiple_domain_mover_->num_domains() > 1 && !freeze_domains_ ) ? 0.7: 0.0;

	if ( check_pairing_dists_ ) 	protocols::farna::print_internal_coords( pose );

	if (dump_pdb_) pose.dump_pdb( "start.pdb" );
	std::cout << "FOLD TREE " << pose.fold_tree();

	Pose start_pose = pose;

	monte_carlo_ = protocols::moves::MonteCarloOP( new protocols::moves::MonteCarlo( pose, *denovo_scorefxn_, m_Temperature_ ) );

	if ( view_monte_carlo_ ) protocols::viewer::add_monte_carlo_viewer( *monte_carlo_, "", 400,400 );

	///////////////////////////////////////////////////////////////////////////
	// Main Loop.
	///////////////////////////////////////////////////////////////////////////
	for (Size n = 1; n <= nstruct_; n++ ) {

std::string const out_file_tag = "S_" + ObjexxFCL::lead_zero_string_of( n, 6 );
		if (tag_is_done_[ out_file_tag ] ) continue; // put this in later!

		pose = start_pose;

		if ( domain_move_frequency_ > 0.0 ) multiple_domain_mover_->randomize_pose_rigid_bodies( pose );

		monte_carlo_->reset( pose );

		TR << "Beginning main loop... " << std::endl;

		frag_size_ = 3;

		clock_t const time_start( clock() );

		for (Size r = 1; r <= rounds_; r++ ) {

			monte_carlo_->score_function( *denovo_scorefxn_ );

			pose = monte_carlo_->lowest_score_pose();

			// Introduce constraints in stages.
			update_pose_constraints( r, rounds_, pose );

			update_domain_rot_trans_mag( r, rounds_ );

			monte_carlo_->reset( pose );
			monte_carlo_->set_temperature( get_temperature( r, rounds_ )  );

			//////////////////////
			// This is it ... do the loop.
			//////////////////////
			for( Size i=1; i <= monte_carlo_cycles_/rounds_ ; ++i ) {
				RNA_move_trial( pose );
			}

			if ( get_native_pose() ) {
				Real const rmsd = rms_at_corresponding_heavy_atoms( *get_native_pose(), pose );
				TR << "All atom rmsd: " << rmsd << std::endl;
			}

			monte_carlo_->recover_low( pose );
			monte_carlo_->show_counters();
			monte_carlo_->reset_counters();
		}

		pose = monte_carlo_->lowest_score_pose();
		denovo_scorefxn_->show( std::cout, pose );

		TR << "Finished fragment assembly of " << out_file_tag << " in " << static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC << " seconds. " << std::endl;

		if ( dump_pdb_ ) pose.dump_pdb( out_file_tag + ".pdb" );
		output_to_silent_file( pose, silent_file_, out_file_tag );

	} //nstruct

}

///////////////////////////////////////////////////////////////////////////////////////////////////////
std::string
CoarseRNA_DeNovoProtocol::get_name() const {
	return "CoarseRNA_DeNovoProtocol¯";
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
core::Real
CoarseRNA_DeNovoProtocol::get_temperature( Size const & r, Size const & rounds ) const{

	Real temperature = m_Temperature_;
	if ( sim_anneal_ ) {
		temperature =  m_Temperature_ * ( (rounds - static_cast<Real>(r)) / rounds );
	}

	std::cout << "Will set temperature to: " << temperature << std::endl;
	return temperature;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
CoarseRNA_DeNovoProtocol::initialize_tag_is_done()
{

	using namespace core::io::silent;

	tag_is_done_.clear();

	utility::vector1< std::string > tags_done;

	SilentFileData silent_file_data;
	if ( utility::file::file_exists( silent_file_ ) ) {
		tags_done = silent_file_data.read_tags_fast( silent_file_ );
		for ( utility::vector1< std::string >::const_iterator iter = tags_done.begin(); iter != tags_done.end(); iter++ ) {
			std::cout << "Already done: " << *iter << std::endl;
			tag_is_done_[ *iter ] = true;
		}
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
CoarseRNA_DeNovoProtocol::output_to_silent_file( core::pose::Pose & pose, std::string const & silent_file, std::string const & out_file_tag, bool const score_only /* = false */) const
{

	using namespace core::io::silent;
	using namespace core::scoring;

	static SilentFileData silent_file_data;
	TR << "Making silent struct for " << out_file_tag << std::endl;

	BinarySilentStruct s( pose, out_file_tag );

	if ( get_native_pose() ){

		s.add_energy( "all_rms", all_atom_rmsd( *get_native_pose(), pose ) );

		std::list< Size > stem_residues( rna_structure_parameters_->get_stem_residues( pose ) );
		Real rmsd_stems( 0.0 );
		if ( stem_residues.size() > 0 ) {
			rmsd_stems = all_atom_rmsd( *get_native_pose(), pose, stem_residues );
			TR << "All atom rmsd over stems: " << rmsd_stems << " for " << out_file_tag << std::endl;
		}
		s.add_energy( "rms_stem", rmsd_stems );

	}

	silent_file_data.write_silent_struct( s, silent_file, score_only );

}


////////////////////////////////////////////////////////////////////////////////////////
void
CoarseRNA_DeNovoProtocol::RNA_move_trial( pose::Pose & pose ) {

	if ( numeric::random::rg().uniform() < domain_move_frequency_ ) {
		random_domain_move_trial( pose );
	} else {
		random_fragment_trial( pose );
	}

}

////////////////////////////////////////////////////////////////
void
CoarseRNA_DeNovoProtocol::fill_pairing_dists( pose::Pose & pose ) {

	using core::id::NamedAtomID;

	pairing_dists_.clear();
	std::map< Size, Size>  const & pairs = rna_structure_parameters_->connections();

	for ( 	std::map< Size, Size >::const_iterator it = pairs.begin();
					it != pairs.end(); it ++ ) {
		Size const i = it->first;
		Size const j = it->second;
		pairing_dists_.push_back( ( pose.xyz( NamedAtomID( " CEN", i ) ) -
																pose.xyz( NamedAtomID( " CEN", j ) ) ).length() );
	}

}


////////////////////////////////////////////////////////////////////////////////////////
void
CoarseRNA_DeNovoProtocol::check_new_pairing_dists( pose::Pose & pose, Size const & frag_pos ) {

	using core::id::NamedAtomID;
	Size count( 0 );
	bool OK( true );

	std::map< Size, Size>  const & pairs = rna_structure_parameters_->connections();

	for ( 	std::map< Size, Size >::const_iterator it = pairs.begin();
					it != pairs.end(); it ++ ) {
		Size const i = it->first;
		Size const j = it->second;
		Real const dist_new = ( pose.xyz( NamedAtomID( " CEN", i ) ) -
														pose.xyz( NamedAtomID( " CEN", j ) ) ).length();
		count++;
		if ( std::abs( dist_new - pairing_dists_[ count ] ) > 1e-3 ) {
			std::cout << "AFTER frag insert at " << frag_pos << " CHANGED PAIRING " << i << " " << j << std::endl;
			OK = false;
		}

	}
	if ( !OK ) std::cout << pose.fold_tree() << std::endl;
}


////////////////////////////////////////////////////////////////////////////////////////
void
CoarseRNA_DeNovoProtocol::random_domain_move_trial( pose::Pose & pose ) {

	if ( freeze_domains_ ) return;

	//	std::cout << "BEFORE: " << (*denovo_scorefxn_)( pose ) << std::endl;

	Size const jumpno = multiple_domain_mover_->apply_and_return_jump( pose );

	//	std::cout << "AFTER1: " << (*denovo_scorefxn_)( pose ) << std::endl;

	rna_loop_closer_->apply_after_jump_change( pose, jumpno );

	//	std::cout << "AFTER2: " << (*denovo_scorefxn_)( pose )<< std::endl;

	monte_carlo_->boltzmann( pose, "domain" );

}


////////////////////////////////////////////////////////////////////////////////////////
void
CoarseRNA_DeNovoProtocol::random_fragment_trial( pose::Pose & pose ) {

	// The check_pairing_dists_ is not typically one -- clutters up code. Remove when we're confident?

	if ( check_pairing_dists_ ) fill_pairing_dists( pose );

	Size const frag_pos = frag_mover_->random_fragment_insertion( pose, frag_size_ );

	if ( check_pairing_dists_ ) check_new_pairing_dists( pose, frag_pos );
	if ( check_pairing_dists_ ) fill_pairing_dists( pose );

	if ( close_loops_ ) rna_loop_closer_->apply( pose, frag_pos );

	Size const dummy_pos( 0 );
	if ( check_pairing_dists_ ) check_new_pairing_dists( pose, dummy_pos );

	//	std::cout << "frag_pos " << frag_pos << std::endl;
	//	std::cout << "; after: " << ( pose.xyz( id::NamedAtomID( " Y  ", 6 ) ) -  pose.xyz( id::NamedAtomID( " Y  ", 7 ) ) ).length() << std::endl;

	monte_carlo_->boltzmann( pose, "frag" + SS(frag_size_) );

}

////////////////////////////////////////////////////////////////////////////////////////
Size
CoarseRNA_DeNovoProtocol::figure_out_constraint_separation_cutoff( Size const & r, Size const & rounds, Size const & max_dist )
{

	//Keep score function coarse for early rounds.
	Real const suppress  = ( r )/(rounds - 2.0);

	Size separation_cutoff = static_cast< Size > ( suppress * max_dist ) + 2;
	if ( separation_cutoff > max_dist ) separation_cutoff = max_dist;

	return separation_cutoff;

}


////////////////////////////////////////////////////////////////////////////////////////
void
CoarseRNA_DeNovoProtocol::update_pose_constraints( Size const & r, Size const & rounds, core::pose::Pose & pose )
{
	using namespace core::scoring::constraints;

	if ( !staged_constraints_) return;

	if ( !constraint_set_ ) return;

	ConstraintSetOP cst_set_new( new scoring::constraints::ConstraintSet );

	static core::kinematics::ShortestPathInFoldTree shortest_path_in_fold_tree( pose.fold_tree() );
	Size const separation_cutoff = figure_out_constraint_separation_cutoff( r, rounds, shortest_path_in_fold_tree.max_dist() );
	TR << "ROUND " << r << " out of " << rounds << std::endl;
	TR << "FOLD_TREE CURRENT SEPARATION CUTOFF " << separation_cutoff << " out of " << shortest_path_in_fold_tree.max_dist() << std::endl;

	ConstraintCOPs csts( constraint_set_->get_all_constraints() );

	for ( Size n = 1; n <= csts.size(); n++ ) {

		ConstraintCOP const & cst( csts[n] );

		if ( cst->natoms() == 2 )  { // currently only defined for pairwise distance constraints.
			Size const i = cst->atom( 1 ).rsd();
			Size const j = cst->atom( 2 ).rsd();
			Size const dist( shortest_path_in_fold_tree.dist( i , j ) );
			if ( dist  > separation_cutoff ) continue;
		}

		cst_set_new->add_constraint( cst );
	}

	pose.constraint_set( cst_set_new );

	TR << "NUM CONSTRAINTS " << pose.constraint_set()->get_all_constraints().size() << " out of " <<
		csts.size() << std::endl;

}

////////////////////////////////////////////////////////////////////////////////////////
void
CoarseRNA_DeNovoProtocol::update_domain_rot_trans_mag( Size const & r, Size const & rounds ){

	Real const scale_factor = ( static_cast<Real>( rounds - r + 1) / rounds );
	Real const rot_mag = 5.0 * scale_factor;
	Real const trans_mag = 1.0 * scale_factor;

	multiple_domain_mover_->update_rot_trans_mag( rot_mag, trans_mag );

}

///////////////////////////////////////////////////////////////////////////////////////////////////////
void
CoarseRNA_DeNovoProtocol::initialize_constraints( core::pose::Pose & pose ) {

	using namespace core::scoring;

	if (pose.constraint_set()->has_constraints() )	{
		denovo_scorefxn_->set_weight( atom_pair_constraint, 1.0 );
		constraint_set_ = pose.constraint_set()->clone();
	}

}



} // namespace rna
} // namespace protocols
