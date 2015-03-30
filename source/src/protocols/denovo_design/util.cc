/// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/util.cc
/// @brief util functions for denovo design of structures
/// @detailed
/// @author Tom Linsky


//Unit Headers
#include <protocols/denovo_design/util.hh>

//Project Headers

//Protocol Headers
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>

//Core Headers
#include <core/pose/Pose.hh>

//Basic/Utility/Numeric Headers
#include <basic/Tracer.hh>

// Boost/ObjexxFCL Headers

//C++ Headers

static thread_local basic::Tracer TR("protocols.denovo_design.components.util");

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace denovo_design {
////////////////////////////////////////////////////////////////////////////////////////////////////

/// @brief Tells whether the two given poses are identical based on # resides and dihedrals
bool same_pose( core::pose::Pose const & pose1, core::pose::Pose const & pose2 )
{
	if (pose1.total_residue() != pose2.total_residue()) {
		return false;
	}

	for (core::Size i = 1; i <= pose1.total_residue(); ++i) {
		if ( pose1.residue(i).name() != pose2.residue(i).name() )
			return false;
		if ( pose1.residue(i).is_protein() && !pose2.residue(i).is_protein() )
			return false;
		if ( !pose1.residue(i).is_protein() && pose2.residue(i).is_protein() )
			return false;
		if ( !pose1.residue(i).is_protein() && !pose2.residue(i).is_protein() )
			continue;

		if (std::abs(pose1.phi(i) - pose2.phi(i)) > 0.0001) {
			return false;
		}
		if (std::abs(pose1.psi(i) - pose2.psi(i)) > 0.0001) {
			return false;
		}
		if (std::abs(pose1.omega(i) - pose2.omega(i)) > 0.0001) {
			return false;
		}
	}
	return true;
}

/// @brief creates a poly-ala pose where every non-gly, non-cyd, protein residue except those in the given set are converted to alanine
void construct_poly_ala_pose(
	core::pose::Pose & pose,
	bool const keep_disulf,
	std::set< core::Size > const & set1,
	std::set< core::Size > const & set2 )
{
	std::set< core::Size > setunion = set1;
	for ( std::set< core::Size >::const_iterator r=set2.begin(), endr=set2.end(); r!=endr; ++r ) {
		setunion.insert( *r );
	}
	construct_poly_ala_pose( pose, keep_disulf, setunion );
}

/// @brief creates a poly-ala pose where every non-gly, non-cyd, protein residue except those in the given set are converted to alanine
void construct_poly_ala_pose(
	core::pose::Pose & pose,
	bool const keep_disulf,
	std::set< core::Size > const & res_set )
{
	utility::vector1< core::Size > positions;
	for ( core::Size i=1, endi=pose.total_residue(); i<=endi; ++i ) {
		if ( pose.residue(i).is_protein() && ( res_set.find(i) != res_set.end() ) ) {
			positions.push_back(i);
		}
		if ( !keep_disulf && ( pose.residue(i).name() == "CYD" ) ) {
			protocols::simple_moves::MutateResidue mut( i, "ALA" );
			mut.apply( pose );
		}
	}
	protocols::toolbox::pose_manipulation::construct_poly_ala_pose(
			pose, positions,
			false, // bool keep_pro,
			true, // bool keep_gly,
			keep_disulf ); // bool keep_disulfide_cys
}

//////////////////////////////////////////////////////////////////////////
/// Output operators for std classes                                   ///
//////////////////////////////////////////////////////////////////////////

/// @brief outputs a matrix
std::ostream &
operator<<( std::ostream & os, numeric::xyzMatrix< core::Real > const & mat ) {
	os << "[ [" << mat.xx() << ", " << mat.xy() << ", " << mat.xz() << "]" << std::endl;
	os << "  [" << mat.yx() << ", " << mat.yy() << ", " << mat.yz() << "]" << std::endl;
 	os << "  [" << mat.zx() << ", " << mat.zy() << ", " << mat.zz() << "] ]";
	return os;	
}

/// @brief outputs a set
std::ostream &
operator<<( std::ostream & os, std::set< core::Size > const & set ) {
	os << "[ ";
	for ( std::set< core::Size >::const_iterator it=set.begin(); it != set.end(); ++it ) {
		os << *it << " ";
	}
	os << "]";
	return os;
}

/// @brief outputs a list of strings
std::ostream &
operator<<( std::ostream & os, std::list< std::string > const & list ) {
	os << "[ ";
	for ( std::list< std::string >::const_iterator c=list.begin(), end=list.end(); c != end; ++c ) {
		os << *c << " ";
	}
	os << "]";
	return os;
}

/// @brief outputs a set
std::ostream &
operator<<( std::ostream & os, std::set< std::string > const & set ) {
	os << "[ ";
	for ( std::set< std::string >::const_iterator it=set.begin(); it != set.end(); ++it ) {
		os << *it << " ";
	}
	os << "]";
	return os;
}

/// @brief outputs a map
std::ostream &
operator<<( std::ostream & os, std::map< core::Size, core::Size > const & map ) {
	os << "{";
	std::map< core::Size, core::Size >::const_iterator it;
	for ( it = map.begin(); it != map.end(); ++it ) {
		if ( it != map.begin() ) {
			os << ", ";
		}
		os << " " << it->first << ":" << it->second;
	}
	os << " }";
	return os;
}

/// @brief outputs a map
std::ostream &
operator<<( std::ostream & os, std::map< std::string, core::Size > const & map ) {
	os << "{";
	std::map< std::string, core::Size >::const_iterator it;
	for ( it = map.begin(); it != map.end(); ++it ) {
		if ( it != map.begin() ) {
			os << ", ";
		}
		os << " " << it->first << ":" << it->second;
	}
	os << " }";
	return os;
}

/// @brief outputs a map
std::ostream &
operator<<( std::ostream & os, std::map< std::pair< std::string, std::string >, core::Size > const & map ) {
	os << "{";
	std::map< std::pair< std::string, std::string >, core::Size >::const_iterator it;
	for ( it = map.begin(); it != map.end(); ++it ) {
		if ( it != map.begin() ) {
			os << ", ";
		}
		os << " " << it->first.first << "__" << it->first.second << ":" << it->second;
	}
	os << " }";
	return os;
}

/// @brief outputs a map
std::ostream &
operator<<( std::ostream & os, std::map< std::string, core::Real > const & map ) {
	os << "{";
	std::map< std::string, core::Real >::const_iterator it;
	for ( it = map.begin(); it != map.end(); ++it ) {
		if ( it != map.begin() ) {
			os << ", ";
		}
		os << " " << it->first << ":" << it->second;
	}
	os << " }";
	return os;
}

} // namespace denovo_design
} // namespace protocols
