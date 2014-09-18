// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file disulfides.hh
/// @brief Some utility functions for dealing with disulfides.
/// @note These were written before I know all the tricks of mini
/// some have been supersceded by core functions.
/// @author Spencer Bliven <blivens@u.washington.edu>
/// @date Created October 2008


#ifndef INCLUDED_apps_pilot_blivens_disulfides_hh
#define INCLUDED_apps_pilot_blivens_disulfides_hh

#include <core/pose/Pose.hh>

//chemistry
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>
// AUTO-REMOVED #include <core/id/AtomID.hh>
#include <core/kinematics/Stub.hh>
#include <basic/Tracer.hh>
#include <core/pose/disulfide_util.hh>
#include <core/conformation/util.hh>

// AUTO-REMOVED #include <utility/vector1.hh>

//Auto Headers
#include <core/chemical/VariantType.hh>


using namespace core;
using namespace std;
using utility::vector1;

static thread_local basic::Tracer TR_apps_pilot_blivens_disulfides_hh( "apps.pilot.blivens.disulfides" );

/*! @brief Determines the distance between two residues.
 *  @param atom The atom to measure from within each residue, eg "CA"
 *  @return -1 on error (for instance, CB distance of Gly)
 */
Real residue_distance(pose::Pose const& pose,
		Size residueA_pos, string atomA,
		Size residueB_pos, string atomB)
{
	Vector a( pose.residue( residueA_pos ).xyz(atomA) );
	Vector b( pose.residue( residueB_pos ).xyz(atomB) );

	return a.distance(b);
}

Real CA_distance(pose::Pose const& pose, Size posA, Size posB)
{
	return residue_distance(pose,posA,"CA",posB,"CA");
}
Real CB_distance(pose::Pose const& pose, Size posA, Size posB)
{
	return residue_distance(pose,posA,"CB",posB,"CB");
}

/*! @brief Decides whether a disulfide bond could exist between two residues based on backbone orientation
 *
 * @details Looks at the backbone distance and orientation between two residues.
 * Returns whether this backbone would be a reasonable place to put a disulfide
 * bond.
 *
 * Current criteria used:
 *  1. Distance between CA-CA
 * possible enhancements:
 *  - CB distance
 *  - general orientation (pointing towards one another)
 *  - ramachandra plot orientation
 *  - other sidechains which need to be repacked/mutated
 */
	bool
reasonable_disulfide_orientation( pose::Pose const& pose, Size residueA_pos, Size residueB_pos)
{
	//Typical Disulfide Ca-Ca distances:
	const Real ds_CA_length_mean(5.50); //Mean
	const Real ds_CA_length_sd(  0.71); //Standard Deviation


	Real r = CA_distance(pose,residueA_pos,residueB_pos);
	//TR_apps_pilot_blivens_disulfides_hh << "Bond "<<residueA_pos<<"-"<<residueB_pos<<" has length "<< r << endl;

	if( ( r < ds_CA_length_mean - 2*ds_CA_length_sd) ||
			( r > ds_CA_length_mean + 2*ds_CA_length_sd) )
		return false;

	//All filters passed!
	return true;
} //end possible_disulfide



/*! @brief Detect whether a disulfide bond could exist.
 */
	bool
possible_disulfide(pose::Pose const& pose, Size residueA_pos, Size residueB_pos)
{
	return reasonable_disulfide_orientation(pose,residueA_pos,residueB_pos);
}

/*! @brief Detect 'actual' disulfide bonds.
 *
 * @details Since disulfide bonds are generally not annotated, we define a
 * disulfide to exist between a pair of residues if
 *   1) They are both cysteines
 *   2) They are bonded (which is marked by conformer.detect_disulfides)
 */
bool
actual_disulfide( pose::Pose const& pose, Size residueA_pos, Size residueB_pos)
{
  return core::conformation::is_disulfide_bond(pose.conformation(), residueA_pos, residueB_pos);
}

void find_disulfides( pose::Pose const & pose, vector1< pair<Size,Size> > /*out*/ & disulfides )
{
	utility::vector1< bool > is_disulfide;
	is_disulfide.resize( pose.total_residue() );
	std::fill( is_disulfide.begin(), is_disulfide.end(), false );

	disulfides.clear();

	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		conformation::Residue A = pose.residue( ii );
		if ( A.aa() == chemical::aa_cys &&
				A.has_variant_type( chemical::DISULFIDE ) &&
				! is_disulfide[ ii ] ) {
			Size const ii_connect_atom( A.atom_index( "SG" ) );
			Size other_res_ii( 0 );
			for ( Size jj = 1; jj <= A.type().n_residue_connections(); ++jj ) {
				if ( (Size) A.type().residue_connection( jj ).atomno() == ii_connect_atom ) {
					other_res_ii = A.connect_map( jj ).resid();
					break;
				}
			}
			if ( other_res_ii == 0 ) {
				std::cerr << "ERROR: Could not find disulfide partner for residue " << ii << std::endl;
				utility_exit();
			}
			assert( other_res_ii > ii );
			//Can only bond residues of the same residue type set (eg centroid to centroid)
			assert( pose.residue_type(other_res_ii).residue_type_set().name() ==
					pose.residue_type(ii).residue_type_set().name() );

			TR_apps_pilot_blivens_disulfides_hh.Info << "Found disulf between " << ii << " and " << other_res_ii << std::endl;
			is_disulfide[ii] = true;
			is_disulfide[other_res_ii] = true;
			disulfides.push_back( std::pair< Size, Size >( ii, other_res_ii ) );

		}
	}
	TR_apps_pilot_blivens_disulfides_hh.Info << "Found " << disulfides.size() << " DS" << std::endl;
}


/*! @brief Determines whether two residues are on the same secondary structure element
 * @details Before running the correct secondary structure must be assigned to
 * each residue.
 *
 * @par Algorithm
 * Walk from one residue to the next, stopping if you cross a ss border.
 * O(longest ss element).
 */
bool same_secondary_structure(pose::Pose const& pose, Size resA_pos, Size resB_pos)
{
	Size l,u;//lower & upper positions
	if( resA_pos < resB_pos) {
		l = resA_pos;
		u = resB_pos;
	}
	else {
		l = resB_pos;
		u = resA_pos;
	}

	char structure = pose.secstruct(l);

	do {
		++l;
		if(pose.secstruct(l) != structure) //ss changed
			return false;
	}while(l<u);

	return true;
}


#endif //INCLUDED_apps_pilot_blivens_disulfides_HH

