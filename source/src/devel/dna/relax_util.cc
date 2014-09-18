// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


// libRosetta headers

#include <devel/dna/relax_util.hh>

#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>

#include <core/scoring/dna/setup.hh>
#include <core/scoring/dna/BasePartner.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.hh>

#include <basic/Tracer.hh>

// Utility Headers
#include <utility/assert.hh>

// Numeric Headers
#include <numeric/random/random.hh>

#include <core/pose/util.hh>
#include <utility/vector1.hh>
#include <numeric/conversions.hh>


namespace devel {
namespace dna {

using namespace core;
using namespace ObjexxFCL;
using utility::vector1;

static thread_local basic::Tracer tt( "devel.dna.relax_util", basic::t_trace );
static thread_local basic::Tracer td( "devel.dna.relax_util", basic::t_debug );
static thread_local basic::Tracer ti( "devel.dna.relax_util", basic::t_info );
static thread_local basic::Tracer tw( "devel.dna.relax_util", basic::t_warning );


/// @details  Add constraints to the pose's constraint set that try to keep the dna chain connected

void
setup_dna_chainbreak_constraints(
	pose::Pose & pose
)
{
	using namespace scoring::constraints;
	using namespace conformation;
	using namespace id;
	using numeric::conversions::radians;

	Real const O3_P_distance( 1.608 );
	Real const O3_angle( 119.8 );
	Real const  P_angle( 103.4 );

	Real const distance_stddev( 0.3 ); // amber is 0.0659
	Real const angle_stddev_degrees( 35 ); // amber is 8.54 (P angle), 5.73 (O3 angle)

	core::scoring::func::FuncOP const distance_func( new core::scoring::func::HarmonicFunc( O3_P_distance, distance_stddev ) );
	core::scoring::func::FuncOP const O3_angle_func( new core::scoring::func::HarmonicFunc( radians( O3_angle ), radians( angle_stddev_degrees ) ) );
	core::scoring::func::FuncOP const  P_angle_func( new core::scoring::func::HarmonicFunc( radians(  P_angle ), radians( angle_stddev_degrees ) ) );

	for ( Size i=1; i< pose.total_residue(); ++i ) {
		Residue const & rsd1( pose.residue( i   ) );
		Residue const & rsd2( pose.residue( i+1 ) );
		if ( rsd1.is_DNA() && !rsd1.is_upper_terminus() && rsd2.is_DNA() && !rsd2.is_lower_terminus() ) {
			tt << "adding dna chainbreak constraint between residues " << i << " and " << i+1 << std::endl;

			AtomID const C3_id( rsd1.atom_index( "C3*" ), i   );
			AtomID const O3_id( rsd1.atom_index( "O3*" ), i   );
			AtomID const  P_id( rsd2.atom_index( "P"   ), i+1 );
			AtomID const O5_id( rsd2.atom_index( "O5*" ), i+1 );

			// distance from O3* to P
			pose.add_constraint( new AtomPairConstraint( O3_id, P_id, distance_func ) );

			// angle at O3*
			pose.add_constraint( new AngleConstraint( C3_id, O3_id, P_id, O3_angle_func ) );

			// angle at P
			pose.add_constraint( new AngleConstraint( O3_id, P_id, O5_id,  P_angle_func ) );
		}
	}

}


/// @details  Choose a random base pair, returns the sequence number of the base-partner that comes first

Size
choose_random_base_pair( pose::Pose const & pose )
{
	using namespace scoring::dna;

	BasePartner const & partner( retrieve_base_partner_from_pose( pose ) );
	vector1< Size > bps;
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( partner[i] > i ) bps.push_back( i );
	}
	assert( bps.size() );
	return bps[ static_cast< int >( numeric::random::rg().uniform() * bps.size() ) + 1 ];
}


/// @details  Choose a random base step jump, returns seqpos of 1st residue, not the jump number

Size
choose_random_base_step_jump( pose::Pose const & pose )
{
	vector1< Size > bs;
	for ( Size i=1; i<= pose.total_residue()-1; ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );
		if ( rsd.is_DNA() && !rsd.is_upper_terminus() && pose.fold_tree().jump_exists( i, i+1 ) ) {
			bs.push_back( i );
		}
	}
	assert( bs.size() );
	return bs[ static_cast< int >( numeric::random::rg().uniform() * bs.size() ) + 1 ];
}


/// @details  Adds jumps into the foldtree to support base-centric kinematics
/// Assumes no intra-dna jumps/cuts to begin with, and NO UNPAIRED BASES
/// (should relax these assumptions).
/// just sets the foldtree's connections, does not fiddle with jump atoms, does not change the pose itself
///

void
add_dna_base_jumps_to_fold_tree(
																pose::Pose const & pose,
																kinematics::FoldTree & f,
																bool const flip
																)
{
	using namespace conformation;
	using namespace scoring::dna;
	using namespace id;
	using namespace kinematics;

	Size const nres( pose.total_residue() );
	assert( f.nres() == nres );

	BasePartner const & partner( retrieve_base_partner_from_pose( pose ) );


	// figure out where dna begins and ends
	Size dna_begin(0), dna_end(0);
	for ( Size i=1; i<= nres; ++i ) {
		Residue const & rsd( pose.residue(i) );
		if ( !dna_begin && rsd.is_DNA() ) dna_begin = i;
		if ( rsd.is_DNA() ) {
			assert( partner[i] );
			if ( i == nres || !pose.residue(i+1).is_DNA() ) { // last dna pos
				dna_end = i;
				break;
			} else {
				assert( !f.is_cutpoint(i) );
			}
		}
	}
	assert( dna_begin && dna_end );

	Size npairs( dna_end - dna_begin + 1 );
	if ( npairs%2 != 0 ) utility_exit_with_message( "bad pose-- unpaired dna?" );
	npairs /= 2;

	for ( Size i=dna_begin; i< dna_begin + npairs; ++i ) {
		Size const ii( i - dna_begin + 1 ); // basepair index
		Size const p_i( partner[i] );
		if ( p_i != dna_end - ii + 1 ) utility_exit_with_message( "bad pose!2" );

		// the basepair jump
		bool const odd( flip ? ( i%2 == 0 ) : ( i%2 == 1 ) );
		int const basepair_cut( odd ? i : p_i - 1 );

		f.new_jump( i, p_i, basepair_cut );

		// the basestep jump
		if ( ii < npairs ) {
			if ( odd ) {
				// from p_i to p_i-1
				f.new_jump( p_i - 1, p_i, p_i - 1 );
			} else {
				f.new_jump( i, i+1, i );
			}
		}
	}

// 	Size root( std::max( Size(1), npairs/2 ) );
// 	if ( flip ) root = partner[ root ];
// 	f.reorder( root );


}

/// @details  Sets the jump-atoms in the foldtree so that the atomtree will have desired connectivity
/// for intra-dna and dna-protein jumps.
/// dna jumps are anchored at the 4th chi1 atom, protein at the C-alpha.
/// Finally, asks the atomtree to update (if necessary, not sure its necessary) the order of the children
/// of the jump anchor atoms so that the jumps have the desired stub, namely:
/// Stub atoms:
/// 1 -- 4th chi1-atom
/// 2 -- 3rd chi1-atom
/// 3 -- 2nd chi1-atom

void
set_dna_jump_atoms( pose::Pose & pose )
{
	using conformation::Residue;
	using namespace id;

	kinematics::FoldTree f( pose.fold_tree() );

	// anchor intra-dna jumps at fourth chi1 atom (out in the base)
	//
	for ( Size i=1; i<= f.num_jump(); ++i ) {
		Residue const & rsd1( pose.residue( f.  upstream_jump_residue( i ) ) );
		Residue const & rsd2( pose.residue( f.downstream_jump_residue( i ) ) );
		if ( rsd1.is_DNA() && rsd2.is_DNA() ) {
			f.set_jump_atoms( i, rsd1.atom_name( rsd1.chi_atoms(1)[4] ), rsd2.atom_name( rsd2.chi_atoms(1)[4] ) );
		} else if ( rsd1.is_DNA() && rsd2.is_protein() ) {
			f.set_jump_atoms( i, rsd1.atom_name( rsd1.chi_atoms(1)[4] ),  "CA" );
		} else if ( rsd2.is_DNA() && rsd1.is_protein() ) {
			f.set_jump_atoms( i, "CA", rsd2.atom_name( rsd2.chi_atoms(1)[4] ) );
		}
	}

	pose.fold_tree( f );

	// tinker with atomorder in atomtree so that we'll get the jump stubs we want
	// useful for graphics -- keeps root base fixed in space
	//
	// I'm not actually sure that this is so important anymore
	// I mean, it seems likely that this will in fact be the stub... but can't remember the details...
	//
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		Residue const & rsd( pose.residue(i) );
		if ( rsd.is_DNA() ) {
			pose.conformation().set_jump_atom_stub_id( StubID ( AtomID( rsd.chi_atoms(1)[4], i ),
																													AtomID( rsd.chi_atoms(1)[3], i ),
																													AtomID( rsd.chi_atoms(1)[2], i ) ) );
		}
	}


}


/// @details  Sets a foldtree for base-centric kinematics in a pose
/// legacy helper code: do not reuse, prefer eg add_dna_base_jumps_to_foldtree and set_dna_base_atoms
void
setup_dna_only_fold_tree( pose::Pose & jump_pose, bool const flip)
{
	using namespace conformation;
	using namespace scoring::dna;
	using namespace id;
	using namespace kinematics;


	BasePartner const & partner( retrieve_base_partner_from_pose( jump_pose ) );
	Size const nres( jump_pose.total_residue() );
	if ( nres%2 != 0 ) utility_exit_with_message( "bad-dna-only-pose!" );
	Size const npairs( nres/ 2 );
	FoldTree f( nres );

	for ( Size i=1; i<= npairs; ++i ) {
		Size const p_i( partner[i] );
		if ( p_i != nres-i+1 ) utility_exit_with_message( "bad-dna-only-pose!2" );

		// the basepair jump
		bool const odd( flip ? ( i%2 == 0 ) : ( i%2 == 1 ) );
		int const basepair_cut( odd ? i : p_i - 1 );

		f.new_jump( i, p_i, basepair_cut );

		// the basestep jump
		if ( i < npairs ) {
			if ( odd ) {
				// from p_i to p_i-1
				f.new_jump( p_i - 1, p_i, p_i - 1 );
			} else {
				f.new_jump( i, i+1, i );
			}
		}

	}

	Size root( std::max( Size(1), npairs/2 ) );
	if ( flip ) root = partner[ root ];
	f.reorder( root );

	for ( Size i=1; i<= f.num_jump(); ++i ) {
		Residue const & rsd1( jump_pose.residue( f.  upstream_jump_residue( i ) ) );
		Residue const & rsd2( jump_pose.residue( f.downstream_jump_residue( i ) ) );
		f.set_jump_atoms( i, rsd1.atom_name( rsd1.chi_atoms(1)[4] ),  rsd2.atom_name( rsd2.chi_atoms(1)[4] ) ); // was C1*, C1*
	}

	jump_pose.fold_tree( f );

	td << "jump_pose new foldtree: " << jump_pose.fold_tree() << std::endl;

	for ( Size i=1; i<= nres; ++i ) {
		Residue const & rsd( jump_pose.residue(i) );
		jump_pose.conformation().set_jump_atom_stub_id( StubID ( AtomID( rsd.chi_atoms(1)[4], i ),
																														 AtomID( rsd.chi_atoms(1)[3], i ),
																														 AtomID( rsd.chi_atoms(1)[2], i ) ) );
	}


}

/// @details  Sets up a DNA-only pose by taking the paired residues in the first DNA chain in start_pose,
/// together with their partners. Sets the intra-base jumps in the foldtree to support base-centric kinematics.


void
setup_dna_only_jump_pose( pose::Pose const & start_pose, pose::Pose & jump_pose )
{
	using namespace conformation;
	using namespace scoring::dna;
	using namespace id;
	using namespace kinematics;

	jump_pose.clear();

	Size const nres( start_pose.total_residue() );

	BasePartner const & partner( retrieve_base_partner_from_pose( start_pose ) );

	std::string const jump_anchor_atom( "C1*" ); // temporary
	for ( Size chain1_begin=1; chain1_begin<= nres; ++chain1_begin ) {
		if ( partner[ chain1_begin ] ) {
			// found first dna residue
			int const chain1( start_pose.chain( chain1_begin ) );

			for ( Size j=chain1_begin; ( j<= nres && partner[j] && start_pose.chain(j) == chain1 ); ++j ) {
				Residue const &         rsd( start_pose.residue(            j ) );
				Residue const & partner_rsd( start_pose.residue( partner[ j ] ) );

				Size const rsd_seqpos( j - chain1_begin + 1 );
				if ( rsd_seqpos == 1 ) {
					jump_pose.append_residue_by_bond( rsd );
				} else {
					jump_pose.insert_residue_by_jump( rsd, rsd_seqpos, rsd_seqpos-1,
																						jump_anchor_atom, // anchor
																						jump_anchor_atom ); // root_atomno
				}

				jump_pose.insert_residue_by_jump( partner_rsd, rsd_seqpos + 1, rsd_seqpos, jump_anchor_atom,
																					jump_anchor_atom );
			}
			break;
		}
	}

	Size const npairs( jump_pose.total_residue() / 2 );
	jump_pose.conformation().insert_chain_ending( npairs );
	core::pose::add_variant_type_to_pose_residue( jump_pose, chemical::LOWER_TERMINUS_VARIANT, 1 );
	core::pose::add_variant_type_to_pose_residue( jump_pose, chemical::LOWER_TERMINUS_VARIANT, npairs+1 );
	core::pose::add_variant_type_to_pose_residue( jump_pose, chemical::UPPER_TERMINUS_VARIANT, npairs );
	core::pose::add_variant_type_to_pose_residue( jump_pose, chemical::UPPER_TERMINUS_VARIANT, jump_pose.total_residue() );

	td << "jump_pose foldtree: " << jump_pose.fold_tree() << std::endl;
	set_base_partner( jump_pose );

	setup_dna_only_fold_tree( jump_pose );

}




/// @details  Delete unpaired DNA bases from the input pose

void
delete_unpaired_bases( pose::Pose & pose )
{
	using namespace kinematics;

	Size const nres( pose.total_residue() );
	vector1< bool > delete_me( nres, false );
	{
		scoring::dna::BasePartner const & partner( scoring::dna::retrieve_base_partner_from_pose( pose ) );
		for ( Size i=1; i<= nres; ++i ) {
			if ( pose.residue(i).is_DNA() && !partner[i] ) {
				delete_me[i] = true;
			}
		}
	}

	// locate nearest undeleted position for rearranging jumps
	vector1< Size > closest( nres, 0 );
	FoldTree const & f( pose.fold_tree() );

	for ( Size i=1; i<= nres; ++i ) {
		if ( delete_me[i] && f.is_jump_point(i) ) {
			// first look backward
			Size j(i);
			while ( !f.is_cutpoint(j-1) ) {
				assert( delete_me[j] );
				--j;
				if ( !delete_me[j] ) break;
			}
			if ( delete_me[j] ) {
				// now try forward
				j = i;
				while ( !f.is_cutpoint(j) ) {
					assert( delete_me[j] );
					++j;
					if ( !delete_me[j] ) break;
				}
			}
			if ( delete_me[j] ) {
				utility_exit_with_message( "remove unpaired bases: bad strand!" );
			}
			assert( !delete_me[j] );
			tt << "sliding deleted jumppoint from " << i << " to " << j << std::endl;
			closest[i] = j;
		}
	}

	FoldTree new_f;
	{
		Size const num_jump( f.num_jump() );
		FArray2D_int jumps(2,num_jump);
		FArray1D_int cuts(num_jump);
		for ( Size i=1; i<= num_jump; ++i ) {
			int pos1( f.jump_point(1,i) );
			int pos2( f.jump_point(2,i) );
			if ( delete_me[pos1] ) pos1 = closest[pos1];
			if ( delete_me[pos2] ) pos2 = closest[pos2];
			jumps(1,i) = pos1;
			jumps(2,i) = pos2;
			cuts(i) = f.cutpoint(i);
		}
		ASSERT_ONLY(bool valid_tree = new_f.tree_from_jumps_and_cuts( nres, num_jump, jumps, cuts );)
		assert( valid_tree );
	}
	pose.fold_tree( new_f );

	// now delete the positions
	for ( Size i=nres; i>=1; --i ) {
		if ( delete_me[i] ) {
			tt << "deleting position: " << i << std::endl;
			pose.delete_polymer_residue(i);
		}
	}

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );
		if ( i>1 && rsd.is_DNA() && pose.chain(i-1) != pose.chain(i) && !rsd.is_lower_terminus() ) {
			tt << "adding lower terminus variant: " << i << std::endl;
			core::pose::add_lower_terminus_type_to_pose_residue( pose, i );
		} else if ( i< pose.total_residue() && rsd.is_DNA() && pose.chain(i) != pose.chain(i+1) && !rsd.is_upper_terminus()){
			tt << "adding upper terminus variant: " << i << std::endl;
			core::pose::add_upper_terminus_type_to_pose_residue( pose, i );
		}
	}

	scoring::dna::set_base_partner( pose );
}

} // ns dna
} // ns devel
