// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RB_geometry - rigid body geometry
/// @brief functions that are needed to do geometric functions for rigid body moves
/// @author Monica Berrondo


#include <protocols/geometry/RB_geometry.hh>

// Rosetta Headers
// AUTO-REMOVED #include <basic/basic.hh>
#include <core/pose/Pose.hh>
#include <protocols/scoring/Interface.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>

// ObjexxFCL Headers

// C++ Headers

//Utility Headers
#include <numeric/conversions.hh>
//#include <numeric/random.functions.hh>
#include <numeric/random/random.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>
//#include <numeric/xyzVector.io.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.geometry.RB_geometry");
static numeric::random::RandomGenerator RG(62451); // <- Magic number, do not change it!!!

using namespace ObjexxFCL;

//     RB_geometry.cc: some supporting functions for rigid body moves

//------------------------------------------------------------------------------
//

namespace protocols {
namespace geometry {

using namespace core;

//////////////////////////////////////////////////////////////////////////
// these functions should probably be moves somewhere else
// copied from rosetta++ jumping_util.cc
numeric::xyzMatrix_double
random_reorientation_matrix(const double phi_range, const double psi_range)
{
	// a genuine rotation matrix which will randomly reorient the coord sys.
	// from Euler theorem
	const double phi( phi_range * RG.uniform() ); // degrees
	const double psi( psi_range * RG.uniform() ); // degrees
	const double theta(
		numeric::conversions::degrees( std::acos(numeric::sin_cos_range( 1.0 - 2.0 *RG.uniform() ) ) )
	); // degrees

	TR << "random_reorientation_matrix phi: " << phi << " psi: " << psi << " theta: " << theta << std::endl;
	return
		numeric::z_rotation_matrix_degrees(  psi   ) *
		numeric::y_rotation_matrix_degrees(  theta ) *
		numeric::z_rotation_matrix_degrees(  phi   );
}


/// @brief Unweighted centroids of all atoms upstream of the jump
///  vs. all atoms downstream of the jump.
/// @details Deliberately includes H -- is this OK?
void
centroids_by_jump(
	core::pose::Pose const & pose,
	core::Size const jump_id,
	core::Vector & upstream_ctrd, //< output
	core::Vector & downstream_ctrd //< output
)
{
	utility::vector1< bool > ok_for_centroid_calculation; //empty is fine.
	centroids_by_jump( pose, jump_id, upstream_ctrd, downstream_ctrd, ok_for_centroid_calculation );
}

/// @brief Unweighted centroids of all atoms upstream of the jump
///  vs. all atoms downstream of the jump.
/// @details Deliberately includes H -- is this OK?
void
centroids_by_jump(
	core::pose::Pose const & pose,
	core::Size const jump_id,
	core::Vector & upstream_ctrd, //< output
	core::Vector & downstream_ctrd, //< output
	utility::vector1< bool > ok_for_centroid_calculation
)
{
	FArray1D_bool is_upstream ( pose.total_residue(), false );
	pose.fold_tree().partition_by_jump( jump_id, is_upstream );

	upstream_ctrd = 0;
	downstream_ctrd = 0;
	int upstream_atoms = 0, downstream_atoms = 0;

	for(int ii = 1, ii_end = pose.total_residue(); ii <= ii_end; ++ii) {
		int & natoms = ( is_upstream(ii) ? upstream_atoms : downstream_atoms );
		core::Vector & ctrd = ( is_upstream(ii) ? upstream_ctrd : downstream_ctrd );
		core::conformation::Residue const & rsd = pose.residue(ii);

		if ( ok_for_centroid_calculation.size() > 0 && !ok_for_centroid_calculation[ ii ] ) {
			//			std::cout << "skipping res  "<< ii << std::endl;
			continue;
		}

		for(int jj = 1, jj_end = rsd.natoms(); jj <= jj_end; ++jj) {
			ctrd += rsd.xyz(jj);
			natoms += 1;
		}
	}

	upstream_ctrd /= upstream_atoms;
	downstream_ctrd /= downstream_atoms;
	//TR << "Upstream:   " << upstream_atoms << " atoms, " << upstream_ctrd << std::endl;
	//TR << "Downstream: " << downstream_atoms << " atoms, " << downstream_ctrd << std::endl;
}

std::pair< core::Vector, core::Vector > centroid_pair_by_jump(
		core::pose::Pose const & pose,
		core::Size jump_id
){
	std::pair<core::Vector, core::Vector > centroids;
	centroids_by_jump(pose, jump_id, centroids.first, centroids.second);
	return centroids;
}

core::Vector upstream_centroid_by_jump(
	core::pose::Pose const & pose,
	core::Size jump_id
){
	return (centroid_pair_by_jump(pose, jump_id)).second;
}

core::Vector downstream_centroid_by_jump(
	core::pose::Pose const & pose,
	core::Size jump_id
){
	return centroid_pair_by_jump(pose, jump_id).second;
}

/// @brief Unweighted centroids of all atoms upstream of the jump
///  vs. all atoms downstream of the jump for interface residues only.
/// needed to calculate rb_centers for fullatom docking - Sid C.
/// @details Deliberately includes H -- is this OK?
void
centroids_by_jump_int(
	core::pose::Pose const & pose,
	core::Size jump_id,
	core::Vector & upstream_ctrd, //< output
	core::Vector & downstream_ctrd //< output
)
{
	FArray1D_bool is_upstream ( pose.total_residue(), false );
	TR.Debug << "fold-tree: " << pose.fold_tree() << std::endl;
	TR.Debug << "partition by jump " << jump_id << std::endl;
	pose.fold_tree().partition_by_jump( jump_id, is_upstream );

	upstream_ctrd = 0;
	downstream_ctrd = 0;
	core::Size upstream_atoms = 0, downstream_atoms = 0;

	protocols::scoring::Interface interface( jump_id );
	//interface.calculate( pose );

	Size int_res_num = 0;
	Real int_distance = 8.0;

  //increments the interface distance until an interface is found, starting from 8A
	for( Size ll = 1; ll <= 5; ll++){
		if (int_res_num == 0){
			interface.distance( int_distance );
			interface.calculate( pose );

			//count interface residues
			for(Size kk = 1; kk <= pose.total_residue(); kk++){
				if (interface.is_interface( kk )) int_res_num++;
			}

			//increment interface distance by 2A
			int_distance = int_distance+2;
		}
	}

	for ( Size ii = 1, ii_end = pose.total_residue(); ii <= ii_end; ++ii ) {
		Size & natoms = ( is_upstream(ii) ? upstream_atoms : downstream_atoms );
		core::Vector & ctrd = ( is_upstream(ii) ? upstream_ctrd : downstream_ctrd );
		core::conformation::Residue const & rsd = pose.residue(ii);
		for ( Size jj = 1, jj_end = rsd.natoms(); jj <= jj_end; ++jj ) {
			if ( interface.is_interface(rsd) ) {
				ctrd += rsd.xyz(jj);
				natoms += 1;
			}
		}
	}

	if ( upstream_atoms == 0 || downstream_atoms == 0 ){
		TR.Warning << "centroids_by_jump_int called but no interface detected!!" << std::endl;
		TR.Warning << "calling centroids_by_jump..." << std::endl;
		centroids_by_jump( pose, jump_id, upstream_ctrd, downstream_ctrd);
	} else {
		upstream_ctrd /= upstream_atoms;
		downstream_ctrd /= downstream_atoms;
	}

	//TR << "Upstream:   " << upstream_atoms << " atoms, " << upstream_ctrd << std::endl;
	//TR << "Downstream: " << downstream_atoms << " atoms, " << downstream_ctrd << std::endl;
}

////////////////////////////////////////////////////////////////////////////////////
/// @begin residue_center_of_mass
///
/// @brief calculates the center of mass of a pose
/// @detailed
///				the start and stop positions (or residues) within the pose are used to
///				find the starting and finishing locations
///
/// @authors Monica Berrondo June 14 2007
///
/// @last_modified June 29 2007
/////////////////////////////////////////////////////////////////////////////////
int
residue_center_of_mass(
	pose::Pose const & pose,
	int const start,
	int const stop
)
{
	Vector center( 0.0 );
	for ( int i=start; i<=stop; ++i ) {
		//if( !pose.residue( i ).is_protein() ) continue;
		if( !pose.residue( i ).is_protein()) {
			Vector ca_pos( pose.residue( i ).nbr_atom_xyz() );
			center += ca_pos;
	 	} else {
			Vector ca_pos( pose.residue( i ).atom( "CA" ).xyz() );
			center += ca_pos;
			}
	}
	center /= (stop-start+1);

	return return_nearest_residue( pose, start, stop, center );
}

////////////////////////////////////////////////////////////////////////////////////
/// @begin return_nearest_residue
///
/// @brief finds the residue nearest some position passed in (normally a
///		center of mass)
/// @detailed
///				the start and stop positions (or residues) within the pose are used to
///				find the starting and finishing locations
///
/// @authors Monica Berrondo June 14 2007
///
/// @last_modified June 29 2007
/////////////////////////////////////////////////////////////////////////////////
int
return_nearest_residue(
	pose::Pose const & pose,
	int const begin,
	int const end,
	Vector center
)
{
	Real min_dist = 9999.9;
	int res = 0;
	for ( int i=begin; i<=end; ++i )
	{
		Vector ca_pos;
		if( !pose.residue( i ).is_protein() ){
			ca_pos = pose.residue( i ).nbr_atom_xyz();
			} else {
		//Vector ca_pos( pose.residue( i ).atom( "CA" ).xyz() );
			ca_pos = pose.residue( i ).atom( "CA" ).xyz() ;
			}

		ca_pos -= center;
		Real tmp_dist( ca_pos.length_squared() );
		if ( tmp_dist < min_dist ) {
			res = i;
			min_dist = tmp_dist;
		}
	}
	return res;
}

} // namespace geometry
} // namespace core
