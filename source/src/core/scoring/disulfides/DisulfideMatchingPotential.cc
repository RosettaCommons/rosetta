// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/disulfides/CentroidDisulfidePotential.cc
/// @brief  Centroid Disulfide Energy Potentials
/// @author rvernon@u.washington.edu
/// @date   02/09/10

// Unit Headers
#include <core/scoring/disulfides/DisulfideMatchingPotential.hh>
#include <core/scoring/disulfides/DisulfideMatchingDatabase.hh>

// Project Headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>
#include <core/kinematics/RT.hh>
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <core/conformation/util.hh>

#include <basic/Tracer.hh>

//Auto Headers
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/kinematics/FoldTree.hh>

//Auto Headers
#include <numeric/xyz.functions.hh>
#include <ObjexxFCL/FArray2A.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option.hh>


static THREAD_LOCAL basic::Tracer TR( "core.scoring.disulfides.CentroidMatchingPotential" );


using namespace core;
using core::scoring::disulfides::DisulfideMatchingPotential;
using core::conformation::Residue;
using std::string;
using utility::vector1;

namespace core {
namespace scoring {
namespace disulfides {

/**
* Constructor
*/
DisulfideMatchingPotential::DisulfideMatchingPotential()
{


}

/**
* Deconstructor
*/
DisulfideMatchingPotential::~DisulfideMatchingPotential() {}

/**
* @brief Calculates scoring terms
*
*/
void
DisulfideMatchingPotential::score_disulfide(
	Residue const & res1,
	Residue const & res2,
	Energy & match_t,
	Energy & match_r,
	Energy & match_rt,
	bool const mirror /*For scoring mirror-image disulfides (DCYS-DCYS)*/
) const
{

	const core::Real probe_radius( basic::options::option[ basic::options::OptionKeys::score::disulf_matching_probe ] );

	//Calculate the distances and angles of this disulfide, mirroring it if so specified
	core::kinematics::RT const scoring_RT( disulfide_RT(res1, res2, mirror) );

	utility::vector1< core::kinematics::RT > db_disulfides( matching_database_.get_all_disulfides() );

	float mt_dist, mr_dist, mrt_dist; //best distance observed
	float r_dist, rt_dist; //current distance being compared

	mr_dist  = 3.0;//2.8; // maybe max should be pi?
	mrt_dist = 3.0+probe_radius;//3.75;

	mt_dist = std::sqrt( scoring_RT.get_translation().distance_squared( db_disulfides[1].get_translation() ));
	r_dist = std::sqrt( scoring_RT.get_rotation().col(1).distance_squared( db_disulfides[1].get_rotation().col(1) ) +
		scoring_RT.get_rotation().col(2).distance_squared( db_disulfides[1].get_rotation().col(2) ) +
		scoring_RT.get_rotation().col(3).distance_squared( db_disulfides[1].get_rotation().col(3) ) );
	rt_dist = core::kinematics::distance( scoring_RT, db_disulfides[1] );

	if ( ( r_dist <= mr_dist ) && ( mt_dist <= probe_radius ) ) mr_dist = r_dist;

	if ( ( rt_dist <= mrt_dist ) && ( mt_dist <= probe_radius ) ) mrt_dist = rt_dist;

	//std::cout << "CHECKING " << db_disulfides.size() << " " << mt_dist << " " << mr_dist << " " << mrt_dist << " " << std::endl;

	for ( Size d = 2; d <= db_disulfides.size(); ++d ) {
		/*
		As far as I can tell, the following four lines were broken before, calculating the same thing for every iteration of this for loop.
		VKM, 17 Aug 2015.
		*/
		float t_dist = std::sqrt( scoring_RT.get_translation().distance_squared( db_disulfides[d].get_translation() ));
		r_dist = std::sqrt( scoring_RT.get_rotation().col(1).distance_squared( db_disulfides[d].get_rotation().col(1) ) +
			scoring_RT.get_rotation().col(2).distance_squared( db_disulfides[d].get_rotation().col(2) ) +
			scoring_RT.get_rotation().col(3).distance_squared( db_disulfides[d].get_rotation().col(3) ) );
		rt_dist = core::kinematics::distance( scoring_RT, db_disulfides[d] );

		//std::cout << "HEYO " << d << " " << mt_dist << " " << mr_dist << " " << mrt_dist << " " << t_dist << " " << r_dist << " " << rt_dist << std::endl;

		if ( t_dist <= mt_dist ) mt_dist = t_dist;
		if ( ( r_dist <= mr_dist ) && ( t_dist <= probe_radius ) ) mr_dist = r_dist;
		if ( ( rt_dist <= mrt_dist ) && ( t_dist <= probe_radius ) ) mrt_dist = rt_dist;
	}

	match_t = mt_dist;
	match_r = mr_dist;
	match_rt = mrt_dist;
}

// Not used by scoring machinery, exists so that other apps can compute the score directly
Energy DisulfideMatchingPotential::compute_matching_energy( pose::Pose const & pose ) const {

	Energy match_RT(0.0);

	utility::vector1< std::pair<Size, Size> > disulfides;
	core::conformation::disulfide_bonds( pose.conformation(), disulfides );

	if ( disulfides.size() > 0 ) {

		for ( Size i = 1; i <= disulfides.size(); ++i ) {

			Energy temp_RT(0.0), junk_rot(0.0), junk_trans(0.0);

			score_disulfide(
				pose.residue(disulfides[i].first),
				pose.residue(disulfides[i].second),
				junk_rot,
				junk_trans,
				temp_RT
			);

			match_RT += temp_RT;
		}
	}

	return match_RT;
}


/**
* @brief calculates some degrees of freedom between two centroid cys residues
*
* If one of the residues is glycine it will be substituted with an idealize
* alanine geometry for the calculations which require a Cb molecule.
*
* centroid_distance requires CEN atoms be defined. If full atom residues
* are specified this function returns centroid_distance of -1.
*
* @param cbcb_distance     The distance between Cbetas squared
* @param centroid_distance The distance between centroids squared
* @param cacbcb_angle_1    The Ca1-Cb1-Cb2 planar angle, in degrees
* @param cacbcb_angle_2    The Ca2-Cb2-Cb1 planar angle, in degrees
* @param cacbcbca_dihedral The Ca1-Cb1-Cb2-Ca2 dihedral angle
* @param backbone_dihedral The N-Ca1-Ca2-C2 dihedral angle
*/
core::kinematics::RT
DisulfideMatchingPotential::disulfide_RT(
	Residue const& res1,
	Residue const& res2,
	bool const mirror
) const {

	core::Real const mirrormult( mirror ? -1.0 : 1.0 );

	Size const MAX_POS( 5 );
	ObjexxFCL::FArray2D_float Epos1(3, MAX_POS), Epos2(3,MAX_POS);

	Epos1(1,2) = res1.atom(res1.atom_index("CA")).xyz().x() * mirrormult;
	Epos1(2,2) = res1.atom(res1.atom_index("CA")).xyz().y();
	Epos1(3,2) = res1.atom(res1.atom_index("CA")).xyz().z();

	Epos1(1,1) = res1.atom(res1.atom_index("N")).xyz().x() * mirrormult;
	Epos1(2,1) = res1.atom(res1.atom_index("N")).xyz().y();
	Epos1(3,1) = res1.atom(res1.atom_index("N")).xyz().z();

	Epos1(1,4) = res1.atom(res1.atom_index("C")).xyz().x() * mirrormult;
	Epos1(2,4) = res1.atom(res1.atom_index("C")).xyz().y();
	Epos1(3,4) = res1.atom(res1.atom_index("C")).xyz().z();

	Epos2(1,2) = res2.atom(res2.atom_index("CA")).xyz().x() * mirrormult;
	Epos2(2,2) = res2.atom(res2.atom_index("CA")).xyz().y();
	Epos2(3,2) = res2.atom(res2.atom_index("CA")).xyz().z();

	Epos2(1,1) = res2.atom(res2.atom_index("N")).xyz().x() * mirrormult;
	Epos2(2,1) = res2.atom(res2.atom_index("N")).xyz().y();
	Epos2(3,1) = res2.atom(res2.atom_index("N")).xyz().z();

	Epos2(1,4) = res2.atom(res2.atom_index("C")).xyz().x() * mirrormult;
	Epos2(2,4) = res2.atom(res2.atom_index("C")).xyz().y();
	Epos2(3,4) = res2.atom(res2.atom_index("C")).xyz().z();

	core::kinematics::RT this_RT(RT_helper::RT_from_epos(Epos1,Epos2));

	return this_RT;
}

void
RT_helper::get_coordinate_system(
	numeric::xyzMatrix_double const & p, //FArray2A_double p, // input
	numeric::xyzMatrix_double & m //FArray2A_double m // output
)
{
	using namespace numeric;

	xyzVector_double a1 = p.col_x() - p.col_y();
	xyzVector_double a2 = p.col_z() - p.col_y();
	a1.normalize();
	xyzVector_double a3 = cross( a1, a2 );
	a3.normalize();
	a2 = cross( a3, a1 );

	m = xyzMatrix_double::cols( a1, a2, a3 );
}

void
RT_helper::get_ncac(
	ObjexxFCL::FArray2A_float pos,
	numeric::xyzMatrix_double & p
)
{
	pos.dimension( 3, 5 );
	using namespace numeric;

	xyzVector_double n( &pos(1,1) );
	xyzVector_double ca( &pos(1,2) );
	xyzVector_double c( &pos(1,4) );

	p = xyzMatrix_double::cols( n, ca, c );

}

numeric::xyzMatrix_double
RT_helper::get_ncac ( ObjexxFCL::FArray2A_float pos )
{
	pos.dimension( 3, 5 );
	using namespace numeric;

	xyzVector_double n( &pos(1,1) );
	xyzVector_double ca( &pos(1,2) );
	xyzVector_double c( &pos(1,4) );

	return xyzMatrix_double::cols( n, ca, c );
}


//helper code to make an RT from two Epos
// does this live somewhere else in mini, haven't found it !
core::kinematics::RT
RT_helper::RT_from_epos( ObjexxFCL::FArray2A_float Epos1, ObjexxFCL::FArray2A_float Epos2)
{
	/// rotation matrix, written in stub1 frame
	core::kinematics::RT::Matrix rotation( 0.0 ); // 3x3
	/// tranlsation vector, written in stub1 frame
	core::kinematics::RT::Vector translation( 0.0 ); // 3

	Size const MAX_POS( 5 ); // param::MAX_POS
	Epos1.dimension(3,MAX_POS);
	Epos2.dimension(3,MAX_POS);

	//bool const local_debug ( false );

	numeric::xyzMatrix_double p1, p2, m1, m2;

	// get coordinate systems from both residues
	get_ncac(Epos1,p1);
	get_ncac(Epos2,p2);
	get_coordinate_system(p1,m1);
	get_coordinate_system(p2,m2);

	// consider this:       |xx xy xz|
	// coordinate frame M = |yx yy yz|
	//                      |zx zy zz|
	// each column is a unit vector written in genuine frame.
	//
	// vector A in frame M can be rewritten as B in genuine frame
	// by the formula B = M x A, thus A = M^T x B
	// a simple example of this would be: the unit vector (1,0,0) in frame M
	// is actually (xx,yx,zx) in genuine frame. mathematically,
	// |xx|   |xx xy xz|   |1|
	// |yx| = |yx yy yz| x |0| ==> B = M x A
	// |zx|   |zx zy zz|   |0|
	//
	// the above formula has another layer of meaning: rotation
	// keeping the genuine frame fixed, a vector can be rotated by applying
	// matrix M onto it, e.g., (1,0,0) rotated to (xx,yx,zx)

	numeric::xyzVector_double e1( &Epos1(1,2) );
	numeric::xyzVector_double e2( &Epos2(1,2) );

	// ( e2 - e1 ) is the vector in genuine frame,
	// translation is the vector in m1 frame. so m1^T is multiplied.
	translation = m1.transposed() * ( e2 - e1 );

	// let's look at the rotation matrix
	// A, B, C are three vectors in genuine frame and are related by rotation
	// B = M1 x A; C = M2 x A;
	// ==> A = M1^T x B = M2^T x C
	// ==> C = M2 x M1^T x B
	// but note that C and B are both in genuine frame and we want a rotation
	// matrix to be applied onto a vector in M1 frame, so comes another step of
	// conversion -- left-multiply M1^T on both sides:
	// M1^T x C = M1^T x M2 x M1^T x B
	// C' = M1^T x M2 x B', as C' and B' are written in M1 frame.
	// so we get the rotation matrix as M1^T x M2.
	// but wait a minute, what Phil orginally got below is M2^T x M1 and it is
	// impossible for that to be wrong, then what happens?

	// It turns out when this rotation matrix is further applied to a vector,
	// it uses Charlies' (col,row) convention (see Dvect_multiply()
	// in RT::make_jump) which means there is one more transpose to do.
	// Now an agreement is reached:
	//  (M2^T x M1)^T = M1^T x (M2^T)^T = M1^T x M2
	// since xyzMatrix uses the normal (row,col) convention, we will switch to
	// rotation = M1^T x M2

	rotation = m1.transposed() * m2;
	/************************Phil's legacy code *********************/
	// rotation(j,*) is the j-th unit vector of 2nd coord sys written in 1st coord-sys
	//  for ( int i=1; i<=3; ++i ) {
	//    for ( int j=1; j<=3; ++j ) {
	//     // DANGER: not sure about the order... ////////////////////////
	//     // should sync with make_jump
	//      rotation(j,i) = Ddotprod( m1(1,i), m2(1,j) ); // 2nd guess
	//      //rotation(j,i) = Ddotprod( m1(1,j), m2(1,i) ); // 1st guess
	//    }
	//   }
	/************************Phil's legacy code ********************/
	core::kinematics::RT rt;
	rt.set_translation( translation );
	rt.set_rotation( rotation );

	return rt;
}


} // disulfides
} // scoring
} // core
