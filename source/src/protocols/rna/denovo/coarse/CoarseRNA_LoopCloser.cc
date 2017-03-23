// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file CoarseRNA_LoopCloser
/// @brief protocols that are specific to CoarseRNA_LoopCloser
/// @details
/// @author Rhiju Das

// Unit headers
#include <protocols/rna/denovo/coarse/CoarseRNA_LoopCloser.hh>

// Package headers
#include <protocols/toolbox/AtomLevelDomainMap.hh>
//
#include <core/types.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <basic/Tracer.hh>
#include <core/scoring/ScoreFunction.hh>
#include <numeric/random/random.hh>
#include <numeric/angle.functions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>
#include <numeric/kinematic_closure/bridgeObjects.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

using namespace core;

static THREAD_LOCAL basic::Tracer TR( "protocols.rna.denovo.coarse.CoarseRNA_LoopCloser" );
using ObjexxFCL::format::I;
using ObjexxFCL::format::F;
using core::id::AtomID;
using core::id::NamedAtomID;
using core::id::DOF_ID;
using core::id::BOGUS_DOF_ID;
using numeric::conversions::radians;
using numeric::conversions::degrees;
using numeric::angle_radians;
using numeric::principal_angle;
using numeric::dihedral_radians;

namespace protocols {
namespace rna {
namespace denovo {
namespace coarse {

CoarseRNA_LoopCloser::CoarseRNA_LoopCloser():
	a_little_verbose_( false ),
	verbose_( false ),
	cutpos_( 0 ),
	nsol_( 0 ),
	which_scratch_res_is_cut_( 0 ),
	choose_least_perturb_solution_( true ),
	choose_best_solution_( false ),
	choose_random_solution_( false ),
	save_all_solutions_( false )
{
	Mover::type("CoarseRNA_LoopCloser");
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details  Apply the RNA Loop Closer -- currently specialized for coarse-grained RNA representation.
///
void
CoarseRNA_LoopCloser::apply( core::pose::Pose & ){
	std::cout << "Does nothing for now -- input the residue that you changed." << std::endl;
	// Later make this loop over all cutpoints?
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
CoarseRNA_LoopCloser::apply( core::pose::Pose & pose, Size const & seqpos_moved )
{

	seqpos_moved_ = seqpos_moved;
	if ( a_little_verbose_ || verbose_ ) TR << "Checking move at " << seqpos_moved << std::endl;

	// which cutpoint was affected?
	partition_definition_.dimension( pose.size(), false );
	pose.fold_tree().partition_by_residue( seqpos_moved, partition_definition_ );

	figure_out_which_cutpoints_were_affected( pose );
	// if ( cutpos_list_.size() > 1 ) {
	//  utility_exit_with_message( "Found more than one cutpoint affected by a frag insertion -- should not happen!" );
	// }

	return close_at_all_cutpoints( pose );
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
CoarseRNA_LoopCloser::apply_after_jump_change( core::pose::Pose & pose, Size const & jumpno )
{

	seqpos_moved_ = 0;
	if ( a_little_verbose_ || verbose_ ) TR << "Checking move at jump " << jumpno << std::endl;

	// which cutpoint was affected?
	partition_definition_.dimension( pose.size(), false );
	pose.fold_tree().partition_by_jump( jumpno, partition_definition_ );

	figure_out_which_cutpoints_were_affected( pose );

	return close_at_all_cutpoints( pose );
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
std::string
CoarseRNA_LoopCloser::get_name() const {
	return "CoarseRNA_LoopCloser¯";
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
CoarseRNA_LoopCloser::close_at_all_cutpoints( core::pose::Pose & pose ){

	if ( cutpos_list_.size() == 0 ) return true; // no change to any cutpoints.

	bool success( true );
	for ( Size i = 1; i <= cutpos_list_.size(); i++ ) {

		cutpos_ = cutpos_list_[ i ];

		if ( a_little_verbose_ || verbose_ ) TR << "Cutpoint affected: " << cutpos_ << std::endl;

		bool const success_in_finding_pivots  = figure_out_pivots( pose );
		if ( !success_in_finding_pivots ) {
			success = false;
			continue;
		}

		close_the_loop( pose );
		//  std::cout << nsol_ << std::endl;

		if ( nsol_ == 0 ) success = false;
		//  if ( nsol_ > 0 ) exit( 0 );
	}

	return success;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
CoarseRNA_LoopCloser::figure_out_which_cutpoints_were_affected( core::pose::Pose const & pose ){

	cutpos_list_.clear();

	if ( verbose_ ) {
		for ( Size n = 1; n <= pose.size(); n++ ) std::cout << partition_definition_( n );
		std::cout << std::endl;
	}

	cutpos_ = 0;
	for ( Size n = 1; n < pose.size(); n++ ) {
		if ( pose.fold_tree().is_cutpoint( n ) &&
				partition_definition_( n ) != partition_definition_( n+1 ) &&
				pose.residue_type( n   ).has_variant_type( chemical::CUTPOINT_LOWER ) &&
				pose.residue_type( n+1 ).has_variant_type( chemical::CUTPOINT_UPPER ) ) {
			cutpos_list_.push_back( n );
			if ( verbose_ ) std::cout << "AFFECTED CUTPOINT: " << n << std::endl;
		}
	}


}

void
CoarseRNA_LoopCloser::output_forward_backward_res(){
	std::cout << "BACKWARD_RES ";
	for ( Size n = 1; n <= backward_res_.size(); n++ ) std::cout << ' ' << backward_res_[n];
	std::cout << std::endl;

	std::cout << "FORWARD_RES ";
	for ( Size n = 1; n <= forward_res_.size(); n++ ) std::cout << ' ' << forward_res_[n];
	std::cout << std::endl << std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////
bool
CoarseRNA_LoopCloser::figure_out_pivots( core::pose::Pose const & pose ){

	if ( !atom_level_domain_map_ ) atom_level_domain_map_ = protocols::toolbox::AtomLevelDomainMapOP( new protocols::toolbox::AtomLevelDomainMap( pose ) );

	figure_out_forward_backward_res_by_backtracking( pose );

	// Filter out perturb residue -- if that leaves enough dofs!
	Size const tot_pivot_res =  forward_res_.size() + backward_res_.size();

	if ( tot_pivot_res < 2 ) return false;

	if ( tot_pivot_res > 2 && seqpos_moved_ > 0 ) {
		remove_res( forward_res_, seqpos_moved_+1 /* +1 because pivot atom (phosphate) is on cutpos+1*/ );
		remove_res( backward_res_, seqpos_moved_+1 );
	}

	if ( verbose_ ) output_forward_backward_res();

	figure_out_pivot_res_and_scratch_res();

	return true;
}

////////////////////////////////////////////////////////////////////////////////////////
void
CoarseRNA_LoopCloser::remove_res( utility::vector1< core::Size > & res_vector, core::Size const & res ){

	utility::vector1< core::Size > new_res_vector;

	for ( Size n = 1; n <= res_vector.size(); n++ ) {
		if ( res_vector[ n ] != res ) new_res_vector.push_back( res_vector[ n ] );
	}

	res_vector = new_res_vector;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
// backtrack through atom_tree to figure out which torsions control each end of the cutpoint.
// Keep track of controlling phosphates (potential "pivots")
void
CoarseRNA_LoopCloser::backtrack( core::kinematics::tree::Atom const * current_atom,
	utility::vector1< core::Size > & upstream_res,
	utility::vector1< bool > & is_upstream_res,
	pose::Pose const & pose ){
	while ( current_atom->parent() ) {

		current_atom = ( current_atom->parent().get() );
		AtomID const & atom_id( current_atom->id() );
		if ( pose.residue_type( atom_id.rsd() ).atom_name( atom_id.atomno() ) == " P  " ) {
			upstream_res.push_back( current_atom->id().rsd() );
			is_upstream_res[ current_atom->id().rsd() ] = true;
		}
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
CoarseRNA_LoopCloser::figure_out_forward_backward_res_by_backtracking( pose::Pose const & pose ){

	is_backward_res_.clear();
	is_forward_res_.clear();
	for ( Size i = 1; i <= pose.size(); i++ ) {
		is_backward_res_.push_back( false );
		is_forward_res_.push_back( false );
	}

	backward_res_.clear();
	forward_res_.clear();

	//Backtrack from sugar before chainbreak -- could generalize this to be the atom immediately before chainbreak.
	core::id::AtomID const ref_atom_id( 2, cutpos_ );
	core::kinematics::tree::Atom const * current_atom ( & pose.atom_tree().atom_dont_do_update( ref_atom_id ) );
	backtrack( current_atom, backward_res_, is_backward_res_, pose );

	// Backtrack from phosphate after chainbreak -- again, could generalize.
	AtomID atom_id( named_atom_id_to_atom_id( NamedAtomID( " P  ", cutpos_ + 1 ), pose ) );
	current_atom = & pose.atom_tree().atom_dont_do_update( atom_id );
	backtrack( current_atom, forward_res_, is_forward_res_, pose );

	if ( verbose_ ) output_forward_backward_res();

	// Trick to filter out just controlling pivots -- any upstream residues
	//  shared in the backtrack paths from either end of the chainbreak
	//  actually will have no effect on the relative positions of the chainbreak ends.
	// Also, this is our opportunity to filter out atoms which are not allowed to move
	//  (takes advantage of the atom_level_domain_map_ object).
	filter_path( backward_res_, is_forward_res_, pose );
	filter_path( forward_res_,  is_backward_res_, pose );

	if ( verbose_ ) output_forward_backward_res();

}

/////////////////////////////////////////////////////////////////////////////////////////////////
void
CoarseRNA_LoopCloser::filter_path( utility::vector1< core::Size > & upstream_res,
	utility::vector1< bool > const & is_filter_res,
	pose::Pose const & pose ){

	utility::vector1< core::Size > new_upstream_res;

	for ( Size n = 1; n <= upstream_res.size(); n++ ) {

		Size const i( upstream_res[ n ] );
		// AtomID atom_id = named_atom_id_to_atom_id( NamedAtomID( " P  ", i ), pose ); // Unused variable causes warning.

		//  if ( !is_filter_res[ i ] &&
		//     atom_level_domain_map_->get( atom_id ) )  new_upstream_res.push_back( i );

		if ( !is_filter_res[ i ] &&
				atom_level_domain_map_->get( id::TorsionID( i, id::BB, 1 ), pose.conformation() ) )  new_upstream_res.push_back( i );

	}

	upstream_res = new_upstream_res;

}

/////////////////////////////////////////////////////////////////////////////////////////////////
// We need 3 pivots -- one better involves the cutpoint.
// For kinematic loop closure, need these three pivots -- and also
//  the residue 5' of the cutpoint (this may already be one of the pivots).
void
CoarseRNA_LoopCloser::figure_out_pivot_res_and_scratch_res(){

	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Size total_pivots = forward_res_.size() + backward_res_.size() + 1;

	if ( total_pivots < 3 ) {
		utility_exit_with_message( "Not enough pivots to close chainbreak!" );
	}

	// line up the pivots in order.
	utility::vector1< Size > pivots_in_order;
	// first the "backward" residues, those upstream of the cutpoint. Their order needs to be reversed.
	for ( Size n = backward_res_.size();  n >= 1; n-- ) pivots_in_order.push_back( backward_res_[ n ] );
	// Cutpoint is one of the pivots.
	pivots_in_order.push_back( cutpos_ + 1 );
	// Then the "forward" residues (downstream).
	for ( Size n = 1; n <= forward_res_.size()  ; n++ ) pivots_in_order.push_back( forward_res_[ n ] );

	// Keep track of which pivots have been selected with an array.
	utility::vector1< bool > pivots_selected( total_pivots, false );

	// Need to pick a random subset of these potential pivot points.
	// well, one is the cutpoint pivot for sure.
	Size const cutpoint_pivot = backward_res_.size() + 1;
	pivots_selected[ cutpoint_pivot ] = true;
	total_pivots--;

	// two more to go.
	for ( Size n = 1; n <= 2; n++ ) {

		Size const which_pivot = static_cast<int>( numeric::random::rg().uniform() * total_pivots ) + 1;
		Size count( 0 );

		for ( Size i = 1; i <= pivots_selected.size(); i++ ) {
			if ( !pivots_selected[ i ] ) count++;
			if ( count == which_pivot ) {
				pivots_selected[ i ] = true;
				total_pivots--;
				break;
			}
		}

	}

	////////////////////
	// manual override!
	////////////////////
	// for( Size i=1; i<=pivots_selected.size();i++) pivots_selected[i]= false;
	// pivots_selected[ 2 ] = true;
	// pivots_selected[ 3 ] = true;
	// pivots_selected[ 4 ] = true;
	////////////////////

	////////////////////////////////////////////////////////////////////////////
	// Okey doke, what are the three pivots?
	// Also, will need to carry out kinematic loop closure on a little "scratch pose"
	// that includes these pivots as well as ...
	// the residues before/after the cutpoint (one or both of these is a pivot as well)
	//  [I think there may be a simpler way to write this ]
	is_scratch_res_.clear();
	is_pivot_res_.clear();
	for ( Size i = 1; i <= is_backward_res_.size() /*pose.size()*/ ; i++ ) {
		is_scratch_res_.push_back( false );
		is_pivot_res_.push_back( false );
	}

	if ( verbose_ ) {
		std::cout << "PIVOTS_SELECTED: ";
		for ( Size i = 1; i <= pivots_selected.size(); i++ ) std::cout << pivots_selected[ i ];
		std::cout << std::endl;
	}

	for ( Size i = 1; i <= pivots_selected.size(); i++ ) {
		if ( pivots_selected[ i ] ) {
			is_pivot_res_[   pivots_in_order[ i ]     ]   = true;
			is_scratch_res_[ pivots_in_order[ i ]     ]   = true;
		}
	}

	is_scratch_res_[ cutpos_ ] = true; /* cutpos+1 is already a pivot_res, this is the residue before.*/

	if ( verbose_ ) {
		std::cout << "IS_PIVOT_RES: ";
		for ( Size i = 1; i <= is_pivot_res_.size(); i++ ) std::cout << is_pivot_res_[ i ];
		std::cout << std::endl;

		std::cout << "IS_SCRATCH_RES: ";
		for ( Size i = 1; i <= is_scratch_res_.size(); i++ ) std::cout << is_scratch_res_[ i ];
		std::cout << std::endl;
	}

	scratch_res_.clear();
	pivot_res_.clear();
	pivot_to_scratch_res_.resize( 3, 0 ); /* vector of size 3*/
	which_scratch_res_is_cut_ = 0;

	for ( Size i = 1; i <= is_scratch_res_.size(); i++ ) {
		if ( is_scratch_res_[ i ] ) {
			scratch_res_.push_back( i );
			if ( i == cutpos_ ) which_scratch_res_is_cut_ = scratch_res_.size();
		}
		if ( is_pivot_res_[ i ] ) {
			pivot_res_.push_back( i );
			pivot_to_scratch_res_[ pivot_res_.size() ] = scratch_res_.size();
		}
	}

	debug_assert( which_scratch_res_is_cut_ > 0 );
	debug_assert( scratch_res_.size() >= 3 );
	debug_assert( pivot_res_.size() == 3 );

	if ( verbose_ ) {
		std::cout << "SCRATCH RES: ";
		for ( Size n = 1; n <= scratch_res_.size(); n++ ) std::cout << ' ' << scratch_res_[ n ];
		std::cout << std::endl;
	}

	if ( a_little_verbose_ || verbose_ ) {
		TR << "PIVOT RES: ";
		for ( Size n = 1; n <= pivot_res_.size(); n++ ) TR << ' ' << pivot_res_[ n ];
		TR << std::endl;
	}

	if ( verbose_ ) {
		std::cout << "SCRATCH RES: ";
		std::cout << "PIVOT_TO_SCRATCH_RES: ";
		for ( Size n = 1; n <= pivot_to_scratch_res_.size(); n++ ) std::cout << ' ' << pivot_to_scratch_res_[ n ];
		std::cout << std::endl;

		std::cout << "WHICH_SCRATCH_RES_IS_CUT: " << which_scratch_res_is_cut_ << std::endl;
	}

}

//////////////////////////////////////////////////////////////////////////
void
CoarseRNA_LoopCloser::close_the_loop( pose::Pose & pose ){

	using namespace core::kinematics;
	using namespace numeric::kinematic_closure;

	///// kinematic loop close.
	// Following copied from, e.g., KinematicMover.cc.  Need to elaborate for terminal residues!
	// inputs to loop closure
	utility::vector1<utility::fixedsizearray1<Real,3> > atoms;
	utility::vector1<Size> pivots (3), order (3);
	utility::vector1<Real> dt_ang, db_len, db_ang;

	// doesn't matter.
	order[1]=1;
	order[2]=2;
	order[3]=3;

	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	// This is a terrible hack to make S-P-S triplets look like
	//  N-CA-C for proteins.
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	atom_ids_.clear();
	Size const nscratch = scratch_res_.size();
	// first triplet --dummy for anchoring. Will be replaced by fixed coordinate system
	atom_ids_.push_back( NamedAtomID( " S  ", scratch_res_[ nscratch ] ) );
	atom_ids_.push_back( NamedAtomID( " P  ", scratch_res_[ nscratch ] ) );
	atom_ids_.push_back( NamedAtomID( " S  ", scratch_res_[ nscratch ] ) );
	// scratch res.
	for ( Size i = 1; i <= nscratch; i++ ) {
		atom_ids_.push_back( NamedAtomID( " S  ", scratch_res_[ i ]-1 ) ); // wait won't this be a problem if scratch_res = 1?
		atom_ids_.push_back( NamedAtomID( " P  ", scratch_res_[ i ]   ) );
		atom_ids_.push_back( NamedAtomID( " S  ", scratch_res_[ i ]   ) );
	}
	// last triplet --dummy for anchoring. Will be replaced by fixed coordinate system
	atom_ids_.push_back( NamedAtomID( " S  ", scratch_res_[ 1 ]-1 ) );
	atom_ids_.push_back( NamedAtomID( " P  ", scratch_res_[ 1 ]   ) );
	atom_ids_.push_back( NamedAtomID( " S  ", scratch_res_[ 1 ]   ) );

	for ( Size i = 1; i <= 3; i++ ) {
		pivots[ i ] = 3 * ( pivot_to_scratch_res_[ i ] ) + 2;
	}

	fill_chainTORS( pose, atom_ids_, atoms, dt_ang, db_ang, db_len );

	if ( verbose_ ) output_chainTORS( dt_ang, db_ang, db_len );

	// These atoms_xyz are the same as computed in fill_chainTORS, but I'm
	// having some trouble sending them out (can't clear or copy vector1< Vector > ?)
	utility::vector1< Vector > atoms_xyz;
	for ( Size i = 1; i <= atoms.size(); i++ ) {
		atoms_xyz.push_back( Vector(atoms[i][1],atoms[i][2],atoms[i][3]) );
	}

	//////////////////////////////////////////////
	// Parameter at chainbreak.
	// This looks a bit weird, because of a hack.
	// In fill_chainTORS, there's a check that one S-P-S triplet
	//  may overlap with the next S-P-S triplet. This definitely
	//  happens at the cutpoint. In that case the degeneracy is broken by
	//  nudging the first S of the second S-P-S triplet a little towards the
	//  centroid of the pre-cutpoint residue. So we can ask for a clean
	//  geometry of the P after the cutpoint with respect to the
	//  previous residue's P-S-CEN coordinate system.
	static const Real idl_S_nextP_( 3.838 );
	Real const d_S_nextP = ( pose.xyz( NamedAtomID(" S  ", cutpos_) ) -
		pose.xyz( NamedAtomID("OVL1", cutpos_) ) ).length();
	if ( verbose_ )  std::cout << " D_S_NEXT_P " << d_S_nextP << " " << idl_S_nextP_ << std::endl;
	db_len[ 3 * which_scratch_res_is_cut_ + 4 ] = d_S_nextP;

	static const Real idl_S_Snudge_nextP_(  65.908 );
	// The "180.0 - " is because Snudge is sort of in the opposite direction to CEN (the math works out).
	// Real const theta_Snudge_S_nextP = 180.0 - degrees( angle_radians( pose.xyz( NamedAtomID(" CEN", cutpos_) ),
	//                                 pose.xyz( NamedAtomID(" S  ", cutpos_) ),
	//                                 pose.xyz( NamedAtomID("OVL1", cutpos_) ) ) );
	Real const theta_S_Snudge_nextP = degrees( angle_radians(atoms_xyz[ 3*which_scratch_res_is_cut_+3],
		atoms_xyz[ 3*which_scratch_res_is_cut_+4],
		pose.xyz( NamedAtomID("OVL1", cutpos_) ) ) );
	if ( verbose_ )  std::cout << " THETA_S_Snudge_NEXTP " << theta_S_Snudge_nextP << " " << idl_S_Snudge_nextP_ << std::endl;
	db_ang[ 3 * which_scratch_res_is_cut_ + 4 ] = theta_S_Snudge_nextP;


	static const Real idl_S_nextP_nextS_( 84.947 );
	Real const theta_S_nextP_nextS = degrees( angle_radians(atoms_xyz[ 3*which_scratch_res_is_cut_+4],
		pose.xyz( NamedAtomID("OVL1", cutpos_) ),
		pose.xyz( NamedAtomID("OVL2", cutpos_) ) ) );
	if ( verbose_ )  std::cout << " THETA_S_NEXTP_NEXTS " << theta_S_nextP_nextS << " " << idl_S_nextP_nextS_ << std::endl;
	db_ang[ 3 * which_scratch_res_is_cut_ + 5 ] = theta_S_nextP_nextS;


	static const Real idl_P_S_CEN_nextP( 86.592 );
	// Real const phi_P_S_Snudge_nextP = degrees( dihedral_radians( atoms_xyz[ 3*which_scratch_res_is_cut_+2],
	//                                atoms_xyz[ 3*which_scratch_res_is_cut_+3],
	//                                atoms_xyz[ 3*which_scratch_res_is_cut_+4],
	//                                pose.xyz( NamedAtomID("OVL1", cutpos_) ) ) );
	Real const phi_P_S_Snudge_nextP = degrees( dihedral_radians( atoms_xyz[ 3*which_scratch_res_is_cut_+2],
		atoms_xyz[ 3*which_scratch_res_is_cut_+3],
		atoms_xyz[ 3*which_scratch_res_is_cut_+4],
		pose.xyz( NamedAtomID("OVL1", cutpos_) ) ) );
	if ( verbose_ ) std::cout << " PHI_P_S_Snudge_NEXTP " << phi_P_S_Snudge_nextP  << " " << idl_P_S_CEN_nextP << std::endl;

	dt_ang[ 3 * which_scratch_res_is_cut_ + 3 ] = phi_P_S_Snudge_nextP;


	//Following are not critical but will allow cutpoint virtual atoms (OVL1, OVL2, and OVU1)
	// to swing to the right place.
	Real const phi_cutpivot1 = degrees( dihedral_radians(  atoms_xyz[ 3*which_scratch_res_is_cut_+3],
		atoms_xyz[ 3*which_scratch_res_is_cut_+4],
		pose.xyz( NamedAtomID("OVL1", cutpos_) ),
		pose.xyz( NamedAtomID("OVL2", cutpos_) ) ) );
	dt_ang[ 3 * which_scratch_res_is_cut_ + 4 ] = phi_cutpivot1;
	if ( verbose_ ) std::cout << " PHI_S_Snudge_nextP_nextS " << phi_cutpivot1  << std::endl;

	Real const phi_cutpivot2 = degrees( dihedral_radians( pose.xyz( NamedAtomID("OVU1", cutpos_+1) ),
		atoms_xyz[ 3*which_scratch_res_is_cut_+5],
		atoms_xyz[ 3*which_scratch_res_is_cut_+6],
		atoms_xyz[ 3*which_scratch_res_is_cut_+7] ) );
	dt_ang[ 3 * which_scratch_res_is_cut_ + 5 ] = phi_cutpivot2;


	// Following *should* be unnecessary, but bridgeObjects reads all atom coordinates
	// before first pivot and after last pivot to define boundary conditions for loop.
	// (In principal bridgeObjects just needs the first three and last three atoms!)
	if ( which_scratch_res_is_cut_ == 1 ) {
		Vector ovl1 = pose.xyz( NamedAtomID( "OVL1", cutpos_ ) );
		for ( Size i = 1; i <= 3; i++ ) {
			atoms[ 3*which_scratch_res_is_cut_ + 5 ][i] = ovl1( i );
		}
	}

	if ( verbose_ ) {
		std::cout <<  "after chainbreak geometry fix" << std::endl;
		output_chainTORS( dt_ang, db_ang, db_len );
	}


	///////////////////////////////////
	// Perform loop closure
	///////////////////////////////////
	t_ang_.clear();
	b_ang_.clear();
	b_len_.clear();
	nsol_ = 0;
	bridgeObjects(atoms, dt_ang, db_ang, db_len, pivots, order, t_ang_, b_ang_, b_len_, nsol_);

	// Let's look at the solutions.
	if ( a_little_verbose_ || verbose_ ) TR << "Kinematic loop closure found this many solutions: " << nsol_ << std::endl;
	if ( nsol_ == 0 ) return;

	figure_out_dof_ids_and_offsets( pose, dt_ang );

	apply_solutions( pose );

}

///////////////////////////////////////////////////////////////////////////////////////
void
CoarseRNA_LoopCloser::figure_out_dof_ids_and_offsets( pose::Pose const & pose,
	utility::vector1<Real> const & dt_ang ){
	////////////////////////////////////////////////////////////////////////////////////
	// Note that the torsion angles that we solved for do not directly correspond to
	// torsion angles in the atom-tree. But they are right up to an *offset*, which
	// we pre-calculate now. Also this is a good time to figure out exactly which
	// DOFs need to be changed in the atom-tree.
	//
	// I think this could be nicely generalized to the fullatom RNA case, or even a
	//    general polymer, but currently don't have the time.
	//
	////////////////////////////////////////////////////////////////////////////////////

	offset_save1_.clear();
	offset_save2_.clear();

	dof_ids1_.clear();
	dof_ids2_.clear();

	DOF_ID dof_id1,dof_id2;

	debug_assert( pivot_res_.size() == 3 );

	AtomID id1, id2, id3, id4;
	for ( Size i = 1; i <= 3; i++ ) {

		Size const pivot = pivot_res_[ i ]; // residue number of phosphate.

		/////////////////////////////////////////////////////////////////////////////////////////
		id1 = named_atom_id_to_atom_id( NamedAtomID( " P  ", pivot-1 ), pose );
		id2 = named_atom_id_to_atom_id( NamedAtomID( " S  ", pivot-1 ), pose );
		if ( pose.fold_tree().is_cutpoint( pivot - 1 ) ) {
			debug_assert( pose.residue_type( pivot-1 ).has_variant_type( chemical::CUTPOINT_LOWER ) ); //this better be the case!!!
			id3 = named_atom_id_to_atom_id( NamedAtomID( "OVL1", pivot-1 ), pose );
			id4 = named_atom_id_to_atom_id( NamedAtomID( "OVL2", pivot-1 ), pose );
		} else {
			id3 = named_atom_id_to_atom_id( NamedAtomID( " P  ", pivot   ), pose );
			id4 = named_atom_id_to_atom_id( NamedAtomID( " S  ", pivot   ), pose );
		}

		dof_id1 = pose.atom_tree().torsion_angle_dof_id( id1, id2, id3, id4 );
		//  std::cout << "PIVOT " << pivot << "  DOF_ID1 " << dof_id1 << "  ID1 " << id1<< "  ID2 " << id2 << "  ID3 " << id3 << "   ID4 " << id4 <<  std::endl;
		dof_ids1_.push_back( dof_id1 );
		figure_out_offset( pose, pivot, dof_id1, dt_ang[ 3*pivot_to_scratch_res_[ i ] + 1 ], offset_save1_ );

		/////////////////////////////////////////////////////////////////////////////////////////
		if ( pose.fold_tree().is_cutpoint( pivot - 1 ) ) {
			debug_assert( pose.residue_type( pivot ).has_variant_type( chemical::CUTPOINT_UPPER ) ); //this better be the case!!!
			id1 = named_atom_id_to_atom_id( NamedAtomID( "OVU1", pivot ), pose );
		} else {
			id1 = named_atom_id_to_atom_id( NamedAtomID( " S  ", pivot-1 ), pose );
		}
		id2 = named_atom_id_to_atom_id( NamedAtomID( " P  ", pivot   ), pose );
		id3 = named_atom_id_to_atom_id( NamedAtomID( " S  ", pivot   ), pose );
		if ( pose.fold_tree().is_cutpoint( pivot ) ) {
			id4 = named_atom_id_to_atom_id( NamedAtomID( " CEN", pivot   ), pose );
		} else {
			id4 = named_atom_id_to_atom_id( NamedAtomID( " P  ", pivot+1 ), pose );
		}

		dof_id2 = pose.atom_tree().torsion_angle_dof_id( id1, id2, id3, id4 );
		//  std::cout << "PIVOT " << pivot << "  DOF_ID2 " << dof_id2 << "  ID1 " << id1<< "  ID2 " << id2 << "  ID3 " << id3 << "   ID4 " << id4 <<  std::endl;
		dof_ids2_.push_back( dof_id2 );
		figure_out_offset( pose, pivot, dof_id2, dt_ang[ 3*pivot_to_scratch_res_[ i ] + 2 ], offset_save2_ );

	}

}

////////////////////////////////////////////////////////////////////////////////////
void
CoarseRNA_LoopCloser::figure_out_offset(
	core::pose::Pose const & pose,
	core::Size const & pivot,
	core::id::DOF_ID const & dof_id,
	core::Real const & original_torsion_value,
	utility::vector1< core::Real > & offset_save ){

	if ( dof_id == BOGUS_DOF_ID ) { //expected at cutpoint!
		//   if ( pivot == cutpos_+1 ) {
		//    offset_save.push_back( 0.0 ); //placeholder
		//    if ( verbose_ ){
		//     std::cout << dof_id;
		//     std::cout << "  offset --- " << std::endl;
		//    }
		//   } else {
		std::cout <<  "Problem with DOF_ID "<< pivot << " " << dof_id << std::endl;
		utility_exit_with_message( "Problem with DOF_ID" );
		//  }
	} else {
		offset_save.push_back( pose.dof( dof_id ) - radians( original_torsion_value ) );
		if ( verbose_ ) {
			std::cout << dof_id;
			std::cout << "  offset " << pose.dof( dof_id ) << " " << radians( original_torsion_value )
				<< " " << pose.dof( dof_id ) - radians( original_torsion_value ) << std::endl;
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
void
CoarseRNA_LoopCloser::apply_solutions( core::pose::Pose & pose ){

	debug_assert( t_ang_.size() == Size( nsol_ ) );

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// Finally, ready to check out the solutions
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	if ( nsol_ == 0 ) return;

	if ( choose_least_perturb_solution_ ) {


		if ( verbose_ )  {
			std::cout << "---------------------------------- " << std::endl;
			std::cout << "   start pose " << std::endl;
			std::cout << "---------------------------------- " << std::endl;
			utility::vector1<Real> dt_ang, db_len, db_ang;
			utility::vector1<utility::fixedsizearray1<Real,3> > atoms;
			fill_chainTORS( pose, atom_ids_, atoms, dt_ang, db_ang, db_len );
			output_chainTORS( dt_ang, db_ang, db_len );
			pose.dump_pdb( "before_closed.pdb" );
		}

		Real best_deviation2( 0.0 );
		Size best_sol( 0 );
		// could save time by just looking over a subset of residues. But I don't think this is rate limiting
		utility::vector1< Vector > ref_vectors;
		Size const ref_atom( 1 );
		for ( Size i = 1; i <= pose.size(); i++ ) {
			ref_vectors.push_back( pose.xyz( id::AtomID(ref_atom,i) ) );
		}

		for ( Size n = 1; n <= Size( nsol_ ); n++ ) {
			fill_solution( pose, n );
			Real deviation2( 0.0 );
			for ( Size i = 1; i <= pose.size(); i++ ) {
				deviation2 += ( pose.xyz( id::AtomID(ref_atom,i) ) - ref_vectors[i] ).length_squared();
			}
			if ( n==1 || deviation2 < best_deviation2 ) {
				best_deviation2 = deviation2;
				best_sol = n;
			}
		}

		fill_solution( pose, best_sol );

		if ( verbose_ ) { //consistency check.

			std::cout << "---------------------------------- " << std::endl;
			std::cout << "   solution " << best_sol << std::endl;
			std::cout << "---------------------------------- " << std::endl;
			output_chainTORS( t_ang_[best_sol], b_ang_[best_sol], b_len_[best_sol] );

		}

		fill_solution( pose, best_sol );

		if ( verbose_ )  {
			std::cout << "pose " << best_sol << ": " << std::endl;
			utility::vector1<Real> dt_ang, db_len, db_ang;
			utility::vector1<utility::fixedsizearray1<Real,3> > atoms;
			fill_chainTORS( pose, atom_ids_, atoms, dt_ang, db_ang, db_len );
			output_chainTORS( dt_ang, db_ang, db_len );
			pose.dump_pdb( "closed.pdb" );
		}

	} else if ( choose_best_solution_ ) {

		debug_assert( scorefxn_ != nullptr );

		Real best_score( 0.0 );
		Size best_sol( 0 );
		for ( Size n = 1; n <= Size( nsol_ ); n++ ) {
			fill_solution( pose, n );
			Real const score = (*scorefxn_)( pose );
			if ( score < best_score || n == 1 ) {
				best_score = score;
				best_sol = n;
			}

			if ( verbose_ && n == 2 ) { //consistency check.
				std::cout << "solution " << n << ": " << std::endl;
				output_chainTORS( t_ang_[n], b_ang_[n], b_len_[n] );

				std::cout << "pose " << n << ": " << std::endl;

				utility::vector1<Real> dt_ang, db_len, db_ang;
				utility::vector1<utility::fixedsizearray1<Real,3> > atoms;
				fill_chainTORS( pose, atom_ids_, atoms, dt_ang, db_ang, db_len );
				output_chainTORS( dt_ang, db_ang, db_len );
			}

		}
		fill_solution( pose, best_sol );

	} else {

		debug_assert( choose_random_solution_ );
		Size const n = static_cast<int>( nsol_ * numeric::random::rg().uniform() ) + 1;
		fill_solution( pose, n );

	}

}

//////////////////////////////////////////////////////////////////////////////////////////////
void
CoarseRNA_LoopCloser::get_all_solutions( core::pose::Pose & pose,
	utility::vector1< core::pose::PoseOP > & pose_list ){

	pose_list.clear();

	for ( Size n = 1; n <= Size( nsol_ ); n++ ) {

		fill_solution( pose, n );

		core::pose::PoseOP pose_save( new Pose );
		*pose_save = pose;
		pose_list.push_back( pose_save );

		if ( verbose_ ) {
			pose.dump_pdb( "KIC_"+ ObjexxFCL::string_of( n )+".pdb" );
		}

	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
void
CoarseRNA_LoopCloser::fill_solution( core::pose::Pose & pose,
	Size const n ) const
{

	for ( Size i = 1; i <= dof_ids1_.size(); i++ ) {

		if ( dof_ids1_[ i ] != BOGUS_DOF_ID )  {
			pose.set_dof( dof_ids1_[i], principal_angle( radians( t_ang_[ n ][ 3 * pivot_to_scratch_res_[ i ] + 1 ] ) + offset_save1_[i] ) );
		}

		if ( dof_ids2_[ i ] != BOGUS_DOF_ID )  {
			pose.set_dof( dof_ids2_[i], principal_angle( radians( t_ang_[ n ][ 3 * pivot_to_scratch_res_[ i ] + 2 ] ) + offset_save2_[i] ) );
		}

	}

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
void
CoarseRNA_LoopCloser::choose_best_solution_based_on_score_function( core::scoring::ScoreFunctionOP scorefxn ){
	choose_best_solution_ = true;
	choose_least_perturb_solution_ = false;
	scorefxn_ = scorefxn;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
void
CoarseRNA_LoopCloser::choose_least_perturb_solution(){
	choose_best_solution_ = false;
	choose_least_perturb_solution_ = true;
}

///////////////////////////////////////////////////////////////////////////
void
CoarseRNA_LoopCloser::output_chainTORS( utility::vector1< core::Real > const & dt_ang,
	utility::vector1< core::Real > const & db_ang,
	utility::vector1< core::Real > const & db_len ) const {

	std::cout << "------  chainTORS output ---- " << std::endl;
	for ( Size i = 1; i <= dt_ang.size(); i++ ) {
		std::cout << I( 3, i ) << " ";
		std::cout << "TORSIONS: ";
		std::cout << F(8,3,dt_ang[ i ]) << " ";

		std::cout << "   BOND_ANGLES: ";
		std::cout << F(8,3,db_ang[ i ]) << " ";

		std::cout << "   BOND_LENGTHS: ";
		std::cout << F(8,3,db_len[ i ]) << " ";

		std::cout << std::endl;

	}
}

///////////////////////////////////////////////////////////
void
CoarseRNA_LoopCloser::fill_chainTORS(
	core::pose::Pose const & pose,
	utility::vector1< id::NamedAtomID> const & atom_ids_,
	utility::vector1<utility::fixedsizearray1<Real,3> > & atoms,
	utility::vector1<Real> & dt_ang,
	utility::vector1<Real> & db_ang,
	utility::vector1<Real> & db_len) const {

	using namespace core::kinematics;
	using namespace numeric::kinematic_closure;

	utility::fixedsizearray1<utility::fixedsizearray1<Real,3>,3 > Q0;
	utility::fixedsizearray1<Real,3> R0;

	utility::vector1< Vector > atoms_xyz;
	for ( Size i = 1; i <= atom_ids_.size(); i++ ) {
		//  std::cout << "filling: " << atom_ids_[i].atomno() << " " << atom_ids_[i].rsd() << std::endl;
		atoms_xyz.push_back( pose.xyz( atom_ids_[ i ] ) );
	}

	//replace first and last with coordinate systems?
	atoms_xyz[ 1 ] = Vector( 0.0, 0.0, 0.0 );
	atoms_xyz[ 2 ] = Vector( 0.0, 0.0, 1.0 );
	atoms_xyz[ 3 ] = Vector( 0.0, 1.0, 0.0 );

	atoms_xyz[ atom_ids_.size() - 2 ] = Vector( 1.0, 0.0, 0.0 );
	atoms_xyz[ atom_ids_.size() - 1 ] = Vector( 1.0, 0.0, 1.0 );
	atoms_xyz[ atom_ids_.size()     ] = Vector( 1.0, 1.0, 0.0 );

	// Some of the pivot atoms may not be distinct -- nan.
	// luckily there's a little hack we can do.
	// where we nudge one atom the slightest bit.
	static Real const nudge( 0.000001 );
	for ( Size n = 1; n <= ( (atom_ids_.size()/3) - 3 ) ; n++ ) {
		Size const i = 3 + (n * 3); // Look at S at the end of one triplet that may overlap with starting S of next triplet
		if ( atom_ids_[ i ] == atom_ids_[ i+1 ] ) {
			// This should be the S at the end of one triplet overlapping with
			// the S beginning the next triplet.
			Size const seqpos = atom_ids_[ i ].rsd();
			atoms_xyz[ i+1 ]  = atoms_xyz[ i ] +
				nudge * ( pose.xyz( NamedAtomID( " CEN", seqpos ) ) - atoms_xyz[i] ).normalize();
		}
	}

	// formatting.
	atoms.clear();
	for ( Size i = 1; i <= atom_ids_.size(); i++ ) {
		utility::fixedsizearray1< Real,3 > atom_xyz_vals;
		atom_xyz_vals[1] = atoms_xyz[i].x();
		atom_xyz_vals[2] = atoms_xyz[i].y();
		atom_xyz_vals[3] = atoms_xyz[i].z();
		atoms.push_back( atom_xyz_vals );
	}

	chainTORS(atoms.size(), atoms, dt_ang, db_ang, db_len, R0, Q0);
}


///////////////////////////////////////////////////////////
void
CoarseRNA_LoopCloser::set_atom_level_domain_map( protocols::toolbox::AtomLevelDomainMapOP atom_level_domain_map ) {
	atom_level_domain_map_ = atom_level_domain_map;
}

///////////////////////////////////////////////////////////
protocols::toolbox::AtomLevelDomainMapOP
CoarseRNA_LoopCloser::atom_level_domain_map(){
	return atom_level_domain_map_;
}


} //coarse
} //denovo
} //rna
} //protocols
