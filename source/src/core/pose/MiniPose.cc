// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/MiniPose.cc
/// @brief  MiniPose class
/// @details minimal class with xyz and fold_tree but not the atom_tree + energies machinery...
/// @author Rhiju Das

// Unit headers
#include <core/pose/Pose.hh>
#include <core/pose/MiniPose.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/FoldTree.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// Utility headers
#include <utility/vector1.hh>

// External headers
#include <ObjexxFCL/format.hh>

// C++ headers
#include <string>

namespace core {
namespace pose {

// @brief Auto-generated virtual destructor
MiniPose::~MiniPose() {}

///////////////////////////////////////////////////////////////////////
MiniPose::MiniPose( core::pose::Pose const & pose )
{
	coords_.clear();
	atom_names_list_.clear();
	variant_types_list_.clear();

	for ( Size i = 1; i <= pose.total_residue(); i++ ) {
		utility::vector1< PointPosition > xyz;
		utility::vector1< std::string > atom_name;
		for ( Size j = 1; j <= pose.residue(i).natoms(); j++ ) {
			xyz.push_back( pose.residue(i).xyz( j ) ) ;
			atom_name.push_back( pose.residue_type(i).atom_name( j ) );
		}
		coords_.push_back( xyz );
		atom_names_list_.push_back(atom_name);
		variant_types_list_.push_back( pose.residue_type( i ).properties().get_list_of_variants() );
	}
	sequence_ = pose.sequence();
	fold_tree_ = pose.fold_tree();
}

//////////////////////////////////////////////////////////////////////////////////
MiniPose::MiniPose( utility::vector1< utility::vector1< PointPosition > > const & coords,
	core::kinematics::FoldTree const & fold_tree,
	std::string const & sequence ):
	coords_( coords ),
	fold_tree_( fold_tree ),
	sequence_( sequence )
{}

//////////////////////////////////////////////////////////////////////////////////
core::kinematics::FoldTree const &
MiniPose::fold_tree() const{
	return fold_tree_;
}

utility::vector1< utility::vector1< PointPosition > > const &
MiniPose::coords() const{
	if ( coords_.size()==0 ) utility_exit_with_message("coords_ is empty!");

	return coords_;
}

utility::vector1< utility::vector1< std::string > > const &
MiniPose::atom_names_list() const{
	if ( atom_names_list_.size()==0 ) utility_exit_with_message("atom_names_list_ is empty!");
	return atom_names_list_;
}

utility::vector1< utility::vector1< std::string > > const &
MiniPose::variant_types_list() const{
	if ( variant_types_list_.size() == 0 ) {
		utility_exit_with_message("variant_types_list_ is empty!");
	}
	return variant_types_list_;
}

Size
MiniPose::size() const {
	return coords_.size();
}

Size
MiniPose::total_residue() const {
	return coords_.size();
}

std::string const &
MiniPose::sequence() const {
	return sequence_;
}

PointPosition const &
MiniPose::xyz( core::id::AtomID atom_id ) const{

	if ( coords_.size()==0 ) utility_exit_with_message("coords_ is empty!");

	if ( coords_.size()<atom_id.rsd() ) {
		std::cout << "atom_id.rsd()= " << atom_id.rsd() << " coords_.size()= " << coords_.size() << std::endl;
		utility_exit_with_message("atom_id.rsd()" +ObjexxFCL::string_of(atom_id.rsd())+ " is out of range!");
	}

	if ( coords_[ atom_id.rsd() ].size()<atom_id.atomno() ) {
		std::cout << "atom_id.atomno()= " << atom_id.atomno() << " coords_[" << atom_id.rsd() << "].size()= " << coords_[ atom_id.rsd() ].size() << std::endl;
		utility_exit_with_message("atom_id.atomno()" +ObjexxFCL::string_of(atom_id.atomno())+ " is out of range!");
	}

	return coords_[ atom_id.rsd() ][ atom_id.atomno() ];
}

std::string const &
MiniPose::atom_name( core::id::AtomID atom_id ) const{
	if ( atom_names_list_.size()==0 ) utility_exit_with_message("atom_names_list_ is empty!");

	if ( atom_names_list_.size()<atom_id.rsd() ) {
		std::cout << "atom_id.rsd()= " << atom_id.rsd() << " atom_names_list_.size()= " << atom_names_list_.size() << std::endl;
		utility_exit_with_message("atom_id.rsd()" +ObjexxFCL::string_of(atom_id.rsd())+ " is out of range!");
	}

	if ( atom_names_list_[ atom_id.rsd() ].size()<atom_id.atomno() ) {
		std::cout << "atom_id.atomno()= " << atom_id.atomno() << " atom_names_list_[" << atom_id.rsd() << "].size()= " << atom_names_list_[ atom_id.rsd() ].size() << std::endl;
		utility_exit_with_message("atom_id.atomno()" +ObjexxFCL::string_of(atom_id.atomno())+ " is out of range!");
	}

	return atom_names_list_[ atom_id.rsd() ][ atom_id.atomno() ];
}

utility::vector1< std::string > const &
MiniPose::variant_types( Size const seq_num ) const
{
	if ( variant_types_list_.size() == 0 ) {
		utility_exit_with_message("variant_types_list_ is empty!");
	}

	if ( variant_types_list_.size() < seq_num ) {
		std::cout << "seq_num= " << seq_num << " variant_types_list_.size()= " << variant_types_list_.size() << std::endl;
		utility_exit_with_message("seq_num" +ObjexxFCL::string_of(seq_num)+ " is out of range!");
	}

	return variant_types_list_[ seq_num ];
}


} // pose
} // core
