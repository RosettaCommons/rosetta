// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/scoring/ImplicitFastClashCheck.cc
/// @brief  does implicit fast clash checking WRT the provided pose
/// @author Will Sheffler (will@sheffler.me)

#include <protocols/scoring/ImplicitFastClashCheck.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/id/AtomID.hh>
#include <ObjexxFCL/format.hh>
#include <numeric/xyz.functions.hh>

#include <core/kinematics/Stub.hh>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>


namespace protocols {
namespace scoring {

/// @details Auto-generated virtual destructor
ImplicitFastClashCheck::~ImplicitFastClashCheck() {}

using core::pose::Pose;
using core::Size;
using core::Real;
using core::id::AtomID;
using core::kinematics::Stub;
using namespace numeric;
using utility::vector1;
typedef numeric::xyzVector<Real> Vec;
typedef numeric::xyzMatrix<Real> Mat;

ImplicitFastClashCheck::ImplicitFastClashCheck(Pose const & pose_in, Real clash_dis, utility::vector1<core::Size> ignore) {
	pose_ = core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose(pose_in) ) );
	utility::vector1<Pose> poses;
	poses.push_back(pose_in);
	init_clash_check( poses, clash_dis, ignore );
}

ImplicitFastClashCheck::ImplicitFastClashCheck(Pose const & pose_in, Real clash_dis) {
	pose_ = core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose(pose_in) ) );
	utility::vector1<Pose> poses;
	poses.push_back(pose_in);
	utility::vector1<core::Size> ignore;
	init_clash_check( poses, clash_dis, ignore );
}

ImplicitFastClashCheck::ImplicitFastClashCheck(
	utility::vector1<core::pose::Pose> const & poses_in,
	core::Real clash_dis,
	utility::vector1<core::Size> ignore
){
	pose_ = core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose(poses_in[1]) ) );
	init_clash_check(poses_in,clash_dis,ignore);
}


void ImplicitFastClashCheck::init_clash_check(utility::vector1<Pose> const & poses, Real neighbor_cutoff, utility::vector1<core::Size> ignore) {
	using numeric::min;
	using numeric::max;
	using numeric::square;
	typedef  numeric::xyzTriple< Size >  CubeDim; // Cube dimensions
	typedef  numeric::xyzTriple< Size >  CubeKey; // Cube index-triple key

	neighbor_cutoff_ = neighbor_cutoff;
	neighbor_cutoff_sq_ = ( neighbor_cutoff*neighbor_cutoff);
	// points_.reserve((pose_->n_residue()-ignore.size())*5);
	// resno_ .reserve((pose_->n_residue()-ignore.size())*5);
	// atomno_.reserve((pose_->n_residue()-ignore.size())*5);
	for ( utility::vector1<Pose>::const_iterator pi = poses.begin(); pi != poses.end(); ++pi ) {
		for ( Size i = 0; i < pi->n_residue(); ++i ) {
			if ( std::find(ignore.begin(),ignore.end(),i+1) != ignore.end() ) continue;
			//Size const natom = min(5ul,pi->residue(i+1).nheavyatoms());
			Size const natom = pi->residue(i+1).nheavyatoms();
			for ( Size j = 1; j <= natom; ++j ) {
				// TODO could check for point redundance here
				points_.push_back( pi->xyz(AtomID(j,i+1)) );
				resno_ .push_back( i+1 );
				atomno_.push_back( j );
			}
		}
	}

	Vec bbu( points_[ 1 ] ); // Lower and upper corners of bounding box
	bbl_ = bbu;
	for ( Size ii = 2; ii <= points_.size(); ++ii ) { bbl_.min( points_[ ii ] ); bbu.max( points_[ ii ] ); }
	bbl_ -= 10 * std::numeric_limits< Real >::epsilon();
	bbu += 10 * std::numeric_limits< Real >::epsilon();
	Real const side( neighbor_cutoff );
	assert( side > Real( 0 ) );
	side_inv_ = ( Real( 1 ) / side );
	cube_dim_ = CubeDim( // Cube dimensions
		Size( std::ceil( ( bbu.x() - bbl_.x() ) * side_inv_ ) ),             // Test that ceil values == Size values
		Size( std::ceil( ( bbu.y() - bbl_.y() ) * side_inv_ ) ),
		Size( std::ceil( ( bbu.z() - bbl_.z() ) * side_inv_ ) )
	);

	cubes_.dimension( cube_dim_.x(), cube_dim_.y(), cube_dim_.z() );

	for ( Size i = 1; i <= points_.size(); ++i ) {
		Vec const pp( points_[ i ]);
		CubeKey const cube_key(
			Size( ( pp.x() - bbl_.x() ) * side_inv_ ) + 1,
			Size( ( pp.y() - bbl_.y() ) * side_inv_ ) + 1,
			Size( ( pp.z() - bbl_.z() ) * side_inv_ ) + 1
		);
		assert( cube_key.x() <= cube_dim_.x() );
		assert( cube_key.y() <= cube_dim_.y() );
		assert( cube_key.z() <= cube_dim_.z() );
		Size i_index = cubes_.index( cube_key.x(), cube_key.y(), cube_key.z() );
		cubes_[ i_index ].push_back( i );
	}
}

bool ImplicitFastClashCheck::clash_check(Vec const & pp ) const {
	Size const icx( Size( ( pp.x() - bbl_.x() ) * side_inv_ ) + 1 );
	Size const icy( Size( ( pp.y() - bbl_.y() ) * side_inv_ ) + 1 );
	Size const icz( Size( ( pp.z() - bbl_.z() ) * side_inv_ ) + 1 );
	for ( Size ix = max( icx, Size( 2 ) ) - 1,  ixe = min( icx + 1, cube_dim_.x() ); ix <= ixe; ++ix ) {
		for ( Size iy = max( icy, Size( 2 ) ) - 1, iye = min( icy + 1, cube_dim_.y() ); iy <= iye; ++iy ) {
			for ( Size iz = max( icz, Size( 2 ) ) - 1, ize = min( icz + 1, cube_dim_.z() ); iz <= ize; ++iz ) {
				Size cube_index = cubes_.index( ix, iy, iz );
				if ( cubes_[ cube_index ].size() != 0 ) { // Cube exists
					for ( vector1<unsigned int>::const_iterator ia = cubes_[ cube_index ].begin(), iae = cubes_[ cube_index ].end(); ia != iae; ++ia ) {
						Vec const j( points_[*ia] );
						Real const d_sq( pp.distance_squared( j ) );
						if ( d_sq <= neighbor_cutoff_sq_ ) {
							// std::cerr << "CLASHCHECK_FAIL_" << d_sq << "_" << resno_[*ia] << "_" << atomno_[*ia] << " ";
							return false;
						}
					}
				}
			}
		}
	}
	return true;
}

platform::uint ImplicitFastClashCheck::clash_count(Vec const & pp ) const {
	platform::uint count = 0;
	Size const icx( Size( ( pp.x() - bbl_.x() ) * side_inv_ ) + 1 );
	Size const icy( Size( ( pp.y() - bbl_.y() ) * side_inv_ ) + 1 );
	Size const icz( Size( ( pp.z() - bbl_.z() ) * side_inv_ ) + 1 );
	for ( Size ix = max( icx, Size( 2 ) ) - 1,  ixe = min( icx + 1, cube_dim_.x() ); ix <= ixe; ++ix ) {
		for ( Size iy = max( icy, Size( 2 ) ) - 1, iye = min( icy + 1, cube_dim_.y() ); iy <= iye; ++iy ) {
			for ( Size iz = max( icz, Size( 2 ) ) - 1, ize = min( icz + 1, cube_dim_.z() ); iz <= ize; ++iz ) {
				Size cube_index = cubes_.index( ix, iy, iz );
				if ( cubes_[ cube_index ].size() != 0 ) { // Cube exists
					for ( vector1<unsigned int>::const_iterator ia = cubes_[ cube_index ].begin(), iae = cubes_[ cube_index ].end(); ia != iae; ++ia ) {
						Vec const j( points_[*ia] );
						Real const d_sq( pp.distance_squared( j ) );
						if ( d_sq <= neighbor_cutoff_sq_ ) {
							++count;
						}
					}
				}
			}
		}
	}
	return count;
}

bool ImplicitFastClashCheck::clash_check(Vec const & pp, Size resno ) const {
	Size const icx( Size( ( pp.x() - bbl_.x() ) * side_inv_ ) + 1 );
	Size const icy( Size( ( pp.y() - bbl_.y() ) * side_inv_ ) + 1 );
	Size const icz( Size( ( pp.z() - bbl_.z() ) * side_inv_ ) + 1 );
	for ( Size ix = max( icx, Size( 2 ) ) - 1,  ixe = min( icx + 1, cube_dim_.x() ); ix <= ixe; ++ix ) {
		for ( Size iy = max( icy, Size( 2 ) ) - 1, iye = min( icy + 1, cube_dim_.y() ); iy <= iye; ++iy ) {
			for ( Size iz = max( icz, Size( 2 ) ) - 1, ize = min( icz + 1, cube_dim_.z() ); iz <= ize; ++iz ) {
				Size cube_index = cubes_.index( ix, iy, iz );
				if ( cubes_[ cube_index ].size() != 0 ) { // Cube exists
					for ( vector1<unsigned int>::const_iterator ia = cubes_[ cube_index ].begin(), iae = cubes_[ cube_index ].end(); ia != iae; ++ia ) {
						Vec const j( points_[*ia] );
						// TR << "resno " << resno << " resno_[*ia] " << resno_[*ia] << std::endl;
						// if(resno_[*ia]==27 && atomno_[*ia]==1) {
						//  std::cout << "clash_check 27 1 " << j.x() << " " << j.y() << " " << j.z() << " " << pp.x() << " " << pp.y() << " " << pp.z() << std::endl;
						// }
						if ( resno == resno_[*ia] && (atomno_[*ia]==2 || atomno_[*ia]==5) ) {
							// std::cout << "ignoring same res clash " << resno_[*ia] << " " << atomno_[*ia] << std::endl;
							continue; // ignore if same res && CA or CB
						}
						Real const d_sq( pp.distance_squared( j ) );
						// if(resno_[*ia]==27 && atomno_[*ia]==1) {
						//  std::cout << "clash_check 27 1 " << d_sq << std::endl;
						// }
						if ( d_sq <= neighbor_cutoff_sq_ ) {
							return false;
							// if( resno == resno_[*ia] ) TR << "clash from same res with res/atm " << resno_[*ia] << " / " << atomno_[*ia] << std::endl;
							// if( resno != resno_[*ia] ) return false;
							// if(sqrt(d_sq)+0.5 <= neighbor_cutoff_) return false;
						}
					}
				}
			}
		}
	}
	return true;
}

// bool ImplicitFastClashCheck::clash_check(Pose const & pose, Size refrsd) const {
//  Stub stubl(pose_->xyz(AtomID(2,1)),pose_->xyz(AtomID(2,2)),pose_->xyz(AtomID(2,3)));
//  Stub stub1(pose  .xyz(AtomID(2,1)),pose  .xyz(AtomID(2,2)),pose  .xyz(AtomID(2,3)));
//  for(Size i = 9; i <= pose.residue(refrsd).nheavyatoms(); ++i) {
//   if( ! clash_check( stubl.local2global(stub1.global2local(pose.xyz(AtomID(i,refrsd+0*pose_->n_residue())))) ) ) return false;
//  }
//  return true;
// }

// bool ImplicitFastClashCheck::clash_check(Stub const & stub, Vec pos) const {
//  // Stub stubl(pose_->xyz(AtomID(2,1)),pose_->xyz(AtomID(2,2)),pose_->xyz(AtomID(2,3)));
//  // return clash_check( stubl.local2global(stub.global2local(pos)) );
//  return clash_check( stub.local2global(pos) );
// }

bool ImplicitFastClashCheck::clash_check_trimer(Pose const & pose, Size refrsd) const {
	Stub stubl(pose_->xyz(AtomID(2,1)),pose_->xyz(AtomID(2,2)),pose_->xyz(AtomID(2,3)));
	Stub stub1(pose  .xyz(AtomID(2,1)),pose  .xyz(AtomID(2,2)),pose  .xyz(AtomID(2,3)));
	for ( Size i = 1; i <= pose.residue(refrsd).nheavyatoms(); ++i ) {
		if ( i > 9 ) if ( ! clash_check( stubl.local2global(stub1.global2local(pose.xyz(AtomID(i,refrsd+0*pose_->n_residue())))) ) ) return false;
		if ( ! clash_check( stubl.local2global(stub1.global2local(pose.xyz(AtomID(i,refrsd+1*pose_->n_residue())))) ) ) return false;
		if ( ! clash_check( stubl.local2global(stub1.global2local(pose.xyz(AtomID(i,refrsd+2*pose_->n_residue())))) ) ) return false;
	}
	Vec cen = stubl.local2global(stub1.global2local(Vec(0,0,0)));
	Vec axs = stubl.local2global(stub1.global2local(Vec(0,0,1)));
	axs = axs - cen;
	Mat rot = numeric::rotation_matrix_degrees(axs,120.0);
	for ( vector1<Vec>::const_iterator i = points_.begin(); i != points_.end(); ++i ) {
		if ( ! clash_check( rot*(*i-cen)+cen ) ) return false;
	}
	return true;
}

void ImplicitFastClashCheck::dump_debug_pdb( utility::io::ozstream & out, core::kinematics::Stub const & stub, char chain ) const {
	using namespace ObjexxFCL::format;
	for ( Size i = 1; i <= points_.size(); ++i ) {
		numeric::xyzVector<Real> p = stub.global2local(points_[i]);
		std::string rname = pose_->residue(resno_[i]).name3();
		std::string aname = pose_->residue(resno_[i]).atom_name(atomno_[i]);
		out << "HETATM" + I( 5, i ) + " "+aname+" "+rname+" "+chain
			+ I( 4, resno_[i] ) + "    "
			+ F( 8, 3, p.x() )
			+ F( 8, 3, p.y() )
			+ F( 8, 3, p.z() )
			+ F( 6, 2, 0.0 ) + ' '
			+ F( 5, 2, 0.0);
		out << "" << std::endl;
	}
}

void ImplicitFastClashCheck::dump_debug_pdb( std::string const & fname, core::kinematics::Stub const & stub, char chain ) const {
	utility::io::ozstream out(fname);
	dump_debug_pdb(out,stub,chain);
}

bool ImplicitFastClashCheck::clash_check_test( numeric::xyzVector<core::Real> const & pp ) const {

	bool bclash = true;
	utility::vector1<numeric::xyzVector<core::Real> >::const_iterator i;
	for ( i = points_.begin(); i != points_.end(); ++i ) {
		Real const d_sq( pp.distance_squared( *i ) );
		if ( d_sq <= neighbor_cutoff_sq_ ) {
			bclash = false;
			return true;
		}
	}
	return false;

	bool clash = clash_check(pp);

	if ( clash != bclash ) utility_exit_with_message("clash check test fails!!!");

	return clash;
}


} // namespace scoring {
} // namespace
