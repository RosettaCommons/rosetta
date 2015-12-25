// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_protocols
/// @brief protocols that are specific to RNA_LoopCloser
/// @details
/// @author Rhiju Das

#include <protocols/farna/movers/RNA_LoopCloser.hh>
#include <protocols/toolbox/AtomLevelDomainMap.hh>
#include <protocols/farna/libraries/RNA_ChunkLibrary.hh> // for ROSETTA_LIBRARY_DOMAIN id.
#include <core/conformation/Residue.hh>
#include <core/chemical/rna/util.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/util.hh>


//Minimizer stuff
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <basic/options/keys/OptionKeys.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>

//#include <numeric/random/random.hh>

// External library headers


//C++ headers
#include <vector>
#include <string>
#include <sstream>

//Auto Headers
#include <core/chemical/VariantType.hh>
#include <core/conformation/Conformation.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/pose/util.hh>
#include <utility/vector1.hh>
#include <numeric/conversions.hh>
#include <ObjexxFCL/string.functions.hh>

//Auto using namespaces
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
//Auto using namespaces end


using namespace core;
using basic::T;

static THREAD_LOCAL basic::Tracer TR( "protocols.farna.movers.RNA_LoopCloser" );

typedef  numeric::xyzMatrix< Real > Matrix;

namespace protocols {
namespace farna {
namespace movers {

RNA_LoopCloser::RNA_LoopCloser():
	verbose_( false ),
	NUM_ROUNDS_( 100 ),
	check_tolerance_( false ),
	ccd_tolerance_( 0.000001 ),
	absolute_ccd_tolerance_( 0.01 ),
	attempt_closure_cutoff_( 20.0 ),
	gap_distance_cutoff_( 8.0 ),
	fast_scan_( true )
{
	Mover::type("RNA_LoopCloser");
}

///////////////////////////////////////////////////////////////////////
void RNA_LoopCloser::apply( core::pose::Pose & pose )
{
	std::map< Size, Size > connections_dummy;
	apply( pose, connections_dummy );
}

///////////////////////////////////////////////////////////////////////////////////////////////
/// @details  Apply the RNA loop closer
///
void RNA_LoopCloser::apply( core::pose::Pose & pose, std::map< Size, Size> const & connections )
{
	utility::vector1< Size > cutpoints_closed = get_cutpoints_closed( pose );
	for ( Size n = 1; n <= cutpoints_closed.size(); n++ ) {

		Size const & i = cutpoints_closed[ n ];
		// do_fast_scan will check if this cutpoint has changed since the
		//  last movement, or is the chainbreak is already well closed,
		//  or if the chainbreak is too big to close.
		if ( fast_scan_ && !passes_fast_scan( pose, i ) )  continue;

		//TR << "Trying to close chain near cutpoint " << i << "-" << i+1 << std::endl;

		apply( pose, connections, i );
	}

}

///////////////////////////////////////////////////////////////////////////////////////////////
std::string
RNA_LoopCloser::get_name() const {
	return "RNA_LoopCloser";
}


///////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
RNA_LoopCloser::get_cutpoints_closed( pose::Pose const & pose ) const {

	using namespace core::id;

	utility::vector1< Size > cutpoints_closed;

	// Loop through all residues and look for potential chainbreaks to close --
	// marked by CUTPOINT_LOWER and CUTPOINT_UPPER variants.
	for ( Size i = 1; i < pose.total_residue(); i++ ) {

		if ( !pose.residue( i   ).has_variant_type( chemical::CUTPOINT_LOWER )  ) continue;
		if ( !pose.residue( i+1 ).has_variant_type( chemical::CUTPOINT_UPPER )  ) continue;

		if ( atom_level_domain_map_ != 0 ) {
			Size const domain1( atom_level_domain_map_->get_domain( NamedAtomID( "OVL1", i ), pose ) );
			Size const domain2( atom_level_domain_map_->get_domain( NamedAtomID( "OVU1", i+1 ), pose ) );
			if ( domain1 == domain2 && domain1 > 0 && domain1 != libraries::ROSETTA_LIBRARY_DOMAIN ) {
				//    TR << TR.Red << "Will not close!!!!!!!!!!!!! " << i << " " << domain1 << std::endl;
				continue;
			}
		}

		cutpoints_closed.push_back( i );
	}

	return cutpoints_closed;
}

///////////////////////////////////////////////////////////////////////////////////////////////
Real RNA_LoopCloser::apply( core::pose::Pose & pose, std::map< Size, Size > const & connections, Size const & cutpoint )
{

	if ( verbose_ ) TR << "Closing loop at: " << cutpoint << std::endl;
	// TR << "Closing loop at: " << cutpoint << std::endl;

	// In the future might have a separate option, e.g., for kinematic loop closure.
	return rna_ccd_close( pose, connections, cutpoint );
}

///////////////////////////////////////////////////////////////////////////////////////////////
Real RNA_LoopCloser::apply( core::pose::Pose & pose, Size const & cutpoint )
{

	if ( verbose_ ) TR << "Closing loop at: " << cutpoint << std::endl;
	// TR << "Closing loop at: " << cutpoint << std::endl;

	std::map< Size, Size > connections; /* empty*/

	// In the future might have a separate option, e.g., for kinematic loop closure.
	return rna_ccd_close( pose, connections, cutpoint );
}

///////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_LoopCloser::apply( core::pose::Pose & pose, utility::vector1< Size > const & cutpoints )
{
	for ( Size n = 1; n <= cutpoints.size(); n++ ) {
		apply( pose, cutpoints[ n ] );
	}
}

////////////////////////////////////////////////////////////////////
// returns false if failure.
bool
RNA_LoopCloser::passes_fast_scan( core::pose::Pose & pose, Size const i ) const
{
	//Don't bother if there hasn't been any movement...
	core::kinematics::DomainMap domain_map;
	pose.conformation().update_domain_map( domain_map );
	if ( domain_map( i ) == domain_map( i+1 ) ) {
		if ( verbose_ ) TR << "Skipping " << i << " due to domain map."  << std::endl;
		return false;
	}

	Real const current_dist_err =   get_dist_err( pose, i );
	//Don't bother if chain is really bad...
	if ( current_dist_err > attempt_closure_cutoff_ )  {
		if ( verbose_ ) TR << "Cutpoint " << i << " will be tough to close: " << current_dist_err << std::endl;
		return false;
	}

	Real const current_gap_distance =   get_gap_distance( pose, i );
	//Don't bother if chain is really bad...
	if ( current_gap_distance > gap_distance_cutoff_ )  {
		if ( verbose_ ) TR << "Cutpoint " << i << " will be tough to close since gap is : " << current_gap_distance << std::endl;
		return false;
	}

	//Don't bother if chain already closed...
	if ( current_dist_err < absolute_ccd_tolerance_ )  {
		if ( verbose_ ) TR << "Cutpoint " << i << " already pretty close to closed: " << current_dist_err << std::endl;
		return false;
	}

	return true;
}

////////////////////////////////////////////////////////////////////////////
bool
RNA_LoopCloser::check_closure( core::pose::Pose const & pose, core::Size const i, Real ccd_tolerance )
{
	if ( ccd_tolerance <= 0.0 ) ccd_tolerance = absolute_ccd_tolerance_;

	runtime_assert( pose.residue( i   ).has_variant_type( chemical::CUTPOINT_LOWER )  );
	runtime_assert( pose.residue( i+1 ).has_variant_type( chemical::CUTPOINT_UPPER )  );

	Real const current_dist_err =  get_dist_err( pose, i );
	// TR << "CURRENT_DIST_ERR  " << current_dist_err << " " << ccd_tolerance << std::endl;

	return ( current_dist_err < ccd_tolerance );
}

////////////////////////////////////////////////////////////////////////////
bool
RNA_LoopCloser::check_closure( core::pose::Pose const & pose, Real ccd_tolerance )
{
	if ( ccd_tolerance <= 0.0 ) ccd_tolerance = absolute_ccd_tolerance_;

	utility::vector1< Size > cutpoints_closed = get_cutpoints_closed( pose );
	for ( Size n = 1; n <= cutpoints_closed.size(); n++ ) {
		Size const & i = cutpoints_closed[ n ];
		if ( !check_closure( pose, i, ccd_tolerance ) ) return false;
	}

	return true;
}


///////////////////////////////////////////////////////////////////////////////////////////////
Real
RNA_LoopCloser::rna_ccd_close( core::pose::Pose & input_pose, std::map< Size, Size > const & connections, Size const & cutpoint ) const
{
	using namespace core::id;

	if ( !input_pose.residue( cutpoint ).is_RNA() ||
			!input_pose.residue( cutpoint+1 ).is_RNA() ) {
		utility_exit_with_message( "RNA CCD closure at "+string_of( cutpoint )+" but residues are not RNA?");
	}

	if ( !input_pose.residue( cutpoint ).has_variant_type( chemical::CUTPOINT_LOWER ) ||
			!input_pose.residue( cutpoint+1 ).has_variant_type( chemical::CUTPOINT_UPPER ) ) {
		utility_exit_with_message( "RNA CCD closure at "+string_of( cutpoint )+" but CUTPOINT_LOWER or CUTPOINT_UPPER variants not properly set up." );
	}

	/////////////////////////////////////////////////////////////////////
	//Just to get all the atom tree connectivities right, make a little scratch pose
	// to play with.
	/////////////////////////////////////////////////////////////////////
	// In the future (if refold could be faster?!), we won't need this
	// scratch pose.

	pose::Pose pose;
	pose.append_residue_by_bond( input_pose.residue( cutpoint)   );
	pose.append_residue_by_jump( input_pose.residue( cutpoint+1), 1 );

	// This is a nice extra option -- instead of just tweaking torsions in the residue right before and after the chainbreak,
	// there are some situations in which it makes sense to tweak torsions in their base pairing partners.
	bool close_two_base_pairs = false;
	Size cutpoint_partner( 0 ), cutpoint_next_partner( 0 );

	if ( connections.find( cutpoint ) != connections.end() &&
			connections.find( cutpoint+1 ) != connections.end() ) {

		cutpoint_partner = connections.find( cutpoint )->second;
		cutpoint_next_partner = connections.find( cutpoint+1 )->second;

		if ( cutpoint_partner == cutpoint_next_partner + 1 ) {
			close_two_base_pairs = true;
		}
	}
	if ( close_two_base_pairs ) {
		TR << "Also including some torsions in partners in CCD: " << cutpoint_partner << " " << cutpoint_next_partner << std::endl;
	}

	// This is totally hard-wired and hacky -- should be easy to fix though.
	if ( close_two_base_pairs ) {
		pose.append_residue_by_jump( input_pose.residue( cutpoint_next_partner /* 9 */ ), 1 );
		pose.append_residue_by_jump( input_pose.residue( cutpoint_partner /* 10 */ ), 1 );
	}

	kinematics::FoldTree f( pose.total_residue() );
	if ( close_two_base_pairs ) {
		// This is totally hard-wired and hacky -- should be easy to fix though.
		f.new_jump( 1, 4, 1 );
		f.set_jump_atoms( 1,
			core::chemical::rna::chi1_torsion_atom( pose.residue(1) ),
			core::chemical::rna::chi1_torsion_atom( pose.residue(4) )   );
		f.new_jump( 2, 3, 2 );
		f.set_jump_atoms( 2,
			core::chemical::rna::chi1_torsion_atom( pose.residue(2) ),
			core::chemical::rna::chi1_torsion_atom( pose.residue(3) )   );

	} else {
		f.new_jump( 1, 2, 1 );
		f.set_jump_atoms( 1,
			core::chemical::rna::chi1_torsion_atom( pose.residue(1) ),
			core::chemical::rna::chi1_torsion_atom( pose.residue(2) )   );
	}

	pose.fold_tree( f );

	// pose.dump_pdb( "scratch.pdb" );

	//Vectors of the three atoms upstream or downstream of the cutpoint that need to be matched.
	utility::vector1 <Vector> upstream_xyzs, downstream_xyzs;
	Real mean_dist_err( -1.0 ), mean_dist_err_prev( 9999.9999 );
	mean_dist_err = get_chainbreak_xyz( pose, 1, upstream_xyzs, downstream_xyzs );
	// This get_chainbreak_xyz will get called repeatedly after each torsion change.

	// What torsions are in play? Note that this could be expanded (or contracted) at will.
	utility::vector1< TorsionID > tor_ids, input_tor_ids;

	// epsilon and zeta before the cutpoint.
	for ( Size j = 5; j <=6; ++ j ) {
		tor_ids.push_back( TorsionID( 1 /*cutpoint*/, BB, j ) );
		input_tor_ids.push_back( TorsionID( cutpoint, BB, j ) );
	}

	// alpha, beta, gamma, delta after the cutpoint?
	for ( Size j = 1; j <=3; ++ j ) {
		tor_ids.push_back( TorsionID( 2 /*cutpoint+1*/, BB, j ) );
		input_tor_ids.push_back( TorsionID( cutpoint+1, BB, j ) );
	}

	//Parin May 5, 2009
	// Also include delta (sugar pucker!) from second residue if its a free floater (cutpoints before and after)
	// if ( input_pose.residue( cutpoint+1 ).has_variant_type( chemical::CUTPOINT_LOWER )  ||
	//    input_pose.residue( cutpoint+1 ).has_variant_type( chemical::UPPER_TERMINUS )  ) {
	//  Size const j = 4;
	//  tor_ids.push_back( TorsionID( 2 /*cutpoint+1*/, BB, j ) );
	//  input_tor_ids.push_back( TorsionID( cutpoint+1, BB, j ) );
	// }

	// Also include delta (sugar pucker!) from second residue if its a free floater (cutpoints before and after)
	// if ( input_pose.residue( cutpoint+1 ).has_variant_type( chemical::CUTPOINT_LOWER )  ||
	//    input_pose.residue( cutpoint+1 ).has_variant_type( chemical::UPPER_TERMINUS )  ) {
	//  Size const j = 4;
	//  tor_ids.push_back( TorsionID( 2 /*cutpoint+1*/, BB, j ) );
	//  input_tor_ids.push_back( TorsionID( cutpoint+1, BB, j ) );
	// }


	//Need to make following user-settable.
	Size nrounds( 0 );

	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////

	while ( nrounds++ < NUM_ROUNDS_ ) {

		if ( check_tolerance_ && ( (mean_dist_err_prev - mean_dist_err) < ccd_tolerance_ ) ) {
			if ( verbose_ ) TR << "Reached tolerance of " << ccd_tolerance_ << " after " << nrounds << " rounds. " << std::endl;
			break;
		}
		mean_dist_err_prev = mean_dist_err;

		/////////////////////////////////////////////////////////////////////
		// One round of chain closure.
		// This doesn't have to be sequential ... allow random traversal,
		// and also biased towards "moveable torsions" (i.e., not beta or epsilon)
		// Also to do? -- check on beta/epsilon -- could perhaps make use of
		//   torsion potential as a simple and general screen.
		for ( Size n = 1; n <= tor_ids.size(); n++ ) {
			TorsionID const & tor_id( tor_ids[ n ] );
			AtomID id1,id2,id3,id4;
			pose.conformation().get_torsion_angle_atom_ids( tor_id, id1, id2, id3, id4 );

			//Is this torsion before or after the chainbreak?
			AtomID my_id;
			Real dir( 0.0 );
			if ( Size( tor_id.rsd() ) == 1 /*cutpoint*/ ) {
				my_id = id3;
				dir = 1;
			} else {
				my_id = id2;
				dir = -1;
			}

			core::kinematics::tree::Atom const & current_atom ( pose.atom_tree().atom( my_id ) );

			kinematics::Stub const & stub_i( current_atom.get_stub() );
			Matrix const & M_i( stub_i.M );
			Vector const & x_i = M_i.col_x();

			Real weighted_sine( 0.0 ), weighted_cosine( 0.0 );
			for ( Size m = 1; m <= upstream_xyzs.size(); m++ ) {
				Vector const current_xyz( current_atom.xyz() );

				Vector const r1 = upstream_xyzs[m] - current_xyz;
				Vector const rho1 = r1 - dot( r1, x_i) * x_i;

				Vector const r2 = downstream_xyzs[m] - current_xyz;
				Vector const rho2 = r2 - dot( r2, x_i) * x_i;

				Real const current_sine   = dir * dot( x_i, cross( rho1, rho2 ) );
				Real const current_cosine = dot( rho1, rho2 );
				//    std::cout << "PREFERRED ANGLE: " << numeric::conversions::degrees( std::atan2( current_sine, current_cosine) ) << std::endl;

				weighted_sine += current_sine;
				weighted_cosine += current_cosine;

				mean_dist_err = get_chainbreak_xyz( pose, 1 /*cutpoint*/, upstream_xyzs, downstream_xyzs );
			}

			Real const twist_torsion = numeric::conversions::degrees( std::atan2( weighted_sine, weighted_cosine) );

			//   std::cout << "CHECK: " << twist_torsion << std::endl;

			Real const current_val = pose.torsion( tor_id );
			pose.set_torsion( tor_id, current_val + twist_torsion );


		}

		if ( verbose_ ) std::cout << "Distance error: " << mean_dist_err << std::endl;

	}

	if ( verbose_ ) pose.dump_pdb( "scratch_close.pdb" );

	/////////////////////////////////////////////////////////////////////
	// OK, done with mini_pose ... copy torsions back into main pose.
	for ( Size n = 1; n <= tor_ids.size(); n++ ) {
		input_pose.set_torsion( input_tor_ids[n], pose.torsion( tor_ids[n] ) );
	}

	if ( verbose_ ) input_pose.dump_pdb( "pose_close.pdb" );

	return mean_dist_err;

}


///////////////////////////////////////////////////////////////
Real
RNA_LoopCloser::get_dist_err( pose::Pose const & pose,
	Size const cutpoint
) const
{
	utility::vector1< Vector > upstream_xyzs;
	utility::vector1< Vector > downstream_xyzs;
	return get_chainbreak_xyz( pose, cutpoint, upstream_xyzs, downstream_xyzs );
}

///////////////////////////////////////////////////////////////
Real
RNA_LoopCloser::get_chainbreak_xyz( pose::Pose const & pose,
	Size const cutpoint,
	utility::vector1< Vector > & upstream_xyzs,
	utility::vector1< Vector > & downstream_xyzs
) const
{
	upstream_xyzs.clear();
	downstream_xyzs.clear();

	upstream_xyzs.push_back( pose.residue( cutpoint ).xyz( " O3'" ) );
	upstream_xyzs.push_back( pose.residue( cutpoint ).xyz( "OVL1" ) );
	upstream_xyzs.push_back( pose.residue( cutpoint ).xyz( "OVL2" ) );

	downstream_xyzs.push_back( pose.residue( cutpoint+1 ).xyz( "OVU1" ) );
	downstream_xyzs.push_back( pose.residue( cutpoint+1 ).xyz( " P  " ) );
	downstream_xyzs.push_back( pose.residue( cutpoint+1 ).xyz( " O5'" ) );

	Real mean_dist_err( 0.0 );
	for ( Size m = 1; m <= upstream_xyzs.size(); m++ ) {
		mean_dist_err += ( upstream_xyzs[m]  - downstream_xyzs[m] ).length();
	}
	mean_dist_err /= upstream_xyzs.size();

	return mean_dist_err;
}

///////////////////////////////////////////////////////////////
Real
RNA_LoopCloser::get_gap_distance( pose::Pose & pose,
	Size const cutpoint
) const
{
	return ( pose.residue(cutpoint).xyz( " O3'" ) - pose.residue(cutpoint).xyz(" C5'" ) ).length();
}

///////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
RNA_LoopCloser::get_extra_cutpoints( pose::Pose const & pose ) const
{

	using namespace core::kinematics;
	using namespace core::scoring::constraints;

	FoldTree const & f( pose.fold_tree() );

	utility::vector1< Size > extra_cutpoints;

	for ( Size cutpos = 1; cutpos < pose.total_residue(); cutpos++ ) {

		if ( ! f.is_cutpoint( cutpos ) ) continue;

		if ( !pose.residue( cutpos   ).has_variant_type( chemical::CUTPOINT_LOWER )  ||
				!pose.residue( cutpos+1 ).has_variant_type( chemical::CUTPOINT_UPPER )  ) {

			// go through pose constraints -- was there any constraint holding residue i's O3'  to i+1's P?
			bool cst_found( false );
			ConstraintCOPs csts( pose.constraint_set()->get_all_constraints() );

			for ( Size n = 1; n <= csts.size(); n++ ) {

				ConstraintCOP const & cst( csts[n] );

				if ( cst->natoms() == 2 )  { // currently only defined for pairwise distance constraints.
					Size const i = cst->atom( 1 ).rsd();
					std::string const name1 = pose.residue( i ).atom_name( cst->atom( 1 ).atomno() );
					Size const j = cst->atom( 2 ).rsd();
					std::string const name2 = pose.residue( j ).atom_name( cst->atom( 2 ).atomno() );

					if ( i == cutpos && j == cutpos+1 && name1 == " O3'" && name2 == " P  " ) {
						cst_found = true; break;
					}
					if ( j == cutpos && i == cutpos+1 && name2 == " O3'" && name1 == " P  " ) {
						cst_found = true; break;
					}

				}

			}

			if ( !cst_found ) continue;

			extra_cutpoints.push_back( cutpos );
		}
	}

	return extra_cutpoints;
}

////////////////////////////////////////////////////////////////////////////////////////
void
RNA_LoopCloser::setup_variants_at_extra_cutpoints( pose::Pose & pose, utility::vector1< Size > const & extra_cutpoints ) const
{
	pose::Pose pose_copy( pose );

	for ( Size n = 1; n <= extra_cutpoints.size(); n++ ) {

		Size cutpos( extra_cutpoints[ n ] );
		//std::cout << "ADDING CUTPOINT VARIANTS " << cutpos << std::endl;

		core::pose::correctly_add_cutpoint_variants( pose, cutpos );

		for ( Size i = cutpos; i <= cutpos + 1; i++ ) {
			for ( Size j = 1; j <= chemical::rna::NUM_RNA_MAINCHAIN_TORSIONS; j++ ) {
				id::TorsionID torsion_id( i, id::BB, j );
				pose.set_torsion( torsion_id, pose_copy.torsion( torsion_id ) ) ;
			} // j
		} // i

	} // n
}


///////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_LoopCloser::remove_variants_at_extra_cutpoints( pose::Pose & pose, utility::vector1< Size > const & extra_cutpoints ) const
{

	pose::Pose pose_copy( pose );

	for ( Size n = 1; n <= extra_cutpoints.size(); n++ ) {
		Size cutpos( extra_cutpoints[ n ] );
		core::pose::remove_variant_type_from_pose_residue( pose, chemical::CUTPOINT_LOWER, cutpos   );
		core::pose::remove_variant_type_from_pose_residue( pose, chemical::CUTPOINT_UPPER, cutpos+1 );

		for ( Size i = cutpos; i <= cutpos + 1; i++ ) {
			for ( Size j = 1; j <= chemical::rna::NUM_RNA_MAINCHAIN_TORSIONS; j++ ) {
				id::TorsionID torsion_id( i, id::BB, j );
				pose.set_torsion( torsion_id, pose_copy.torsion( torsion_id ) ) ;
			} // j
		} // i

	} // n

}

///////////////////////////////////////////////////////////////////////////////////////////////
// void
// RNA_LoopCloser::close_loops_carefully(
//   core::pose::Pose & pose,
//   core::scoring::ScoreFunctionOP const & scorefxn,
//  Size close_loops_rounds )
// {

//  for (Size k = 1; k <= close_loops_rounds; k++ ) {
//   std::cout << "ROUND " << k << " of " << close_loops_rounds << ": " << std::endl;
//   scorefxn->show( std::cout, pose );
//   close_loops_carefully_one_round( pose, scorefxn );
//  }

// }

///////////////////////////////////////////////////////////////////////////////////////////////
// void
// RNA_LoopCloser::close_loops_carefully(
//   core::pose::Pose & pose,
//   core::scoring::ScoreFunctionOP const & scorefxn )
// {
//  close_loops_carefully( pose, scorefxn, 1 );
// }

///////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_LoopCloser::close_loops_carefully( core::pose::Pose & pose,  std::map< Size, Size > const & connections )
{
	using namespace core::scoring;

	// ScoreFunctionOP scorefxn_local( scorefxn->clone() );

	// A quick minimization
	//static protocols::farna::RNA_Minimizer rna_minimizer;
	// rna_minimizer.set_score_function( scorefxn_local );
	// rna_minimizer.apply( pose );
	// pose.dump_pdb( "dump1.pdb" );

	// OK, close the loops...
	// Look through constraints -- there may be some cutpoints to close that are not
	// yet flagged... set up chainbreak variants there...
	utility::vector1< Size > const extra_cutpoints = get_extra_cutpoints( pose );
	setup_variants_at_extra_cutpoints( pose, extra_cutpoints );

	// Now locally minimize around chainbreaks -- ignore chainbreak atoms because they're
	// probably in crazy places (wait until CCD).
	//scorefxn_local->set_weight( linear_chainbreak, 0.0 );
	// local_minimize_at_chainbreaks( pose, scorefxn_local );
	// pose.dump_pdb( "dump2.pdb" );

	// CCD loop close.
	// bool fast_scan_save( fast_scan_ );
	// fast_scan_ = false;
	apply( pose, connections );
	// fast_scan_ = fast_scan_save;
	// pose.dump_pdb( "dump3.pdb" );


	// Now minimize with strong coordinate constraints.
	tight_minimize( pose );

	// local_minimize_at_chainbreaks( pose, scorefxn_local );

	// Leave the pose sequence/variants as we received it.
	remove_variants_at_extra_cutpoints( pose, extra_cutpoints );
	// pose.dump_pdb( "dump4.pdb" );

	//In case the pose is about to be output ...
	// pose.constraint_set( cst_set_save );
	// (*scorefxn)( pose );

}

//////////////////////////////////////////////////////////////////////////////////
void
RNA_LoopCloser::tight_minimize( core::pose::Pose & pose ) const
{
	using namespace core::optimization;
	using namespace core::kinematics;
	using namespace core::scoring;

	core::scoring::constraints::ConstraintSetOP const & cst_set_save( pose.constraint_set()->clone() );
	pose.constraint_set( 0 );
	core::scoring::constraints::add_coordinate_constraints( pose );

	static AtomTreeMinimizer minimizer;
	float const dummy_tol( 0.0000025);
	bool const use_nblist( true );
	static MinimizerOptions options( "dfpmin", dummy_tol, use_nblist, false /*deriv_check_*/, false /*deriv_check_*/ );
	static ScoreFunctionOP lores_scorefxn( ScoreFunctionFactory::create_score_function( core::scoring::RNA_LORES_WTS ) );
	lores_scorefxn->set_weight( linear_chainbreak, 5.0 );
	lores_scorefxn->set_weight( coordinate_constraint, 2.0 );

	MoveMap mm;
	mm.set_bb( true );
	mm.set_jump ( false );
	mm.set_chi( false );
	minimizer.run( pose, mm, *lores_scorefxn, options );

	pose.constraint_set( cst_set_save );
}

/////////////////////////////////////////////////////////////////////////////////
void
RNA_LoopCloser::local_minimize_at_chainbreaks(
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionOP & scorefxn_local ) const
{

	using namespace core::optimization;

	AtomTreeMinimizer minimizer;
	float const dummy_tol( 0.0000025);
	bool const use_nblist( true );
	MinimizerOptions options( "dfpmin", dummy_tol, use_nblist, false /*deriv_check_*/, false /*deriv_check_*/ );

	kinematics::MoveMap mm;
	mm.set_bb( false );
	mm.set_chi( false );
	mm.set_jump( false );

	scorefxn_local->set_weight( core::scoring::linear_chainbreak, 1.0 );

	kinematics::FoldTree const & f( pose.fold_tree() );

	for ( Size cutpos = 1; cutpos < pose.total_residue(); cutpos++ ) {

		mm.set_bb( false );

		if ( ! f.is_cutpoint( cutpos ) ) continue;

		if ( !pose.residue( cutpos   ).has_variant_type( chemical::CUTPOINT_LOWER )  ||
				!pose.residue( cutpos+1 ).has_variant_type( chemical::CUTPOINT_UPPER )  ) continue;

		mm.set_bb( cutpos, true );
		mm.set_bb( cutpos+1, true );
		minimizer.run( pose, mm, *scorefxn_local, options );

		TR << "Trying to minimize chain near cutpoint " << cutpos << std::endl;
	}
}


} //movers
} //farna
} //protocols
