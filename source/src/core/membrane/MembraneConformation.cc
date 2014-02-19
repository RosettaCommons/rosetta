// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 	 core/membrane/MembraneConformation.cc
///
/// @brief 	 Membrane Conformation
/// @details The Membrane Conformation is responsible for:
///             - maintaining a correct membrane foldtree
///             - maintaining references to the membrane and embedding residues
///             - providing access to membrane related data
///
/// @note    Last Modified 1/9/14
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_MembraneConformation_cc
#define INCLUDED_core_membrane_MembraneConformation_cc

// Unit Headers
#include <core/membrane/MembraneConformation.hh>

// Project Headers
#include <core/membrane/properties/SpanningTopology.hh>
#include <core/membrane/properties/LipidAccInfo.hh>

// Package Headers
#include <core/conformation/Conformation.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/membrane/util/Exceptions.hh>

#include <core/conformation/util.hh>
#include <core/types.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>

using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.membrane.MembraneConformation");

using namespace core::conformation;

namespace core {
namespace membrane {

//// Constructors //////////////////////////

/// @brief Default Constructor
MembraneConformation::MembraneConformation() :
	Conformation(),
	membrane_(1)
{}

/// @brief Standard Constructor
MembraneConformation::MembraneConformation(
		Conformation const & conf,
		utility::vector1< std::pair< int, int > > embres_map,
		int membrane
) :
		Conformation( conf ),
		embres_map_( embres_map ),
		membrane_( membrane )

{
	build_membrane_foldtree();
	assert( is_membrane() );
}

/// @brief Copy Constructor
MembraneConformation::MembraneConformation( MembraneConformation const & src ) :
		Conformation( src ),
		embres_map_( src.embres_map() ),
		membrane_( src.membrane() ),
		lipid_acc_data_( src.lipid_acc_data() ),
		spanning_topology_( src.spanning_topology() )
{
	build_membrane_foldtree();
	assert( is_membrane() );
}

/// @brief Assignemnt Operator
Conformation &
MembraneConformation::operator=( MembraneConformation const & src )
{
	Conformation::operator=( src );
	membrane_ = src.membrane();
	embres_map_ = src.embres_map();
	spanning_topology_ = src.spanning_topology_;
	lipid_acc_data_ = src.lipid_acc_data_;
	return *this;
}

/// @brief Clone Method
ConformationOP
MembraneConformation::clone() const {
	return new MembraneConformation(*this);
}

/// @brief Return Membrane Embedding Map
utility::vector1< std::pair< int, int > >
MembraneConformation::embres_map() const { return embres_map_; }

/// @brief Return the index of the membrane root
int
MembraneConformation::membrane() const { return membrane_; }

/// @brief Assert Membrane Conformation Invariants
bool MembraneConformation::is_membrane() {

	// Check Class Invariants

	// Check membrane root is within the bounds of the pose
	if ( membrane_ <= 0 || membrane_ > (int) size() ) {
		TR << "Membrane root specified" << membrane_ << "must be between 1 and " << size() << std::endl;
		return false;
	}

	// Check that there is an embedding residue anchored by jump for each protein chain
	if ( embres_map_.size() != num_chains()-1 ) {
		TR << "There must exist an embedding residue for every protein chain in the pose" << std::endl;
		TR << "Specified " << embres_map_.size() << " embedding residues but needs " << num_chains()-1 << std::endl;
		return false;
	}

	// Check that all of the jump anchors are within the bounds of the pose
	for ( core::Size i = 1; i <= num_chains()-1; ++i ) {
		if ( embres_map_[ i ].first <= 0 || embres_map_[ i ].first > (int) size() ) {
			TR << "Embedding residue specified must be within the pose (between 1 and " << size() << ")" << std::endl;
			TR << "Specified " << embres_map_[ i ].first << std::endl;
			return false;
		}
	}

	// Check Fold tree
	if ( !fold_tree().check_fold_tree() ) {
		TR << "Membrane foldtree is an invalid foldtree!" << std::endl;
		return false;
	}

	// Check the membrane embedding residue is the root of the pose
	if ( !fold_tree().is_root( membrane_ ) ) {
		TR << "Membrane root " << membrane_ << " must be the root of the foldtree!" << std::endl;
		return false;
	}

	// Else return true
	return true;
}

/// Manage Derived Spanning and Lipid Acc Data ///////////

void
MembraneConformation::add_topology_by_chain( core::membrane::properties::SpanningTopology sp, core::Size ) {
	spanning_topology_.push_back( sp );
}

utility::vector1< properties::SpanningTopology >
MembraneConformation::spanning_topology() const {
	return spanning_topology_;
}

/// @brief Add Lipid Accessibility Info
void
MembraneConformation::add_lips_by_chain( properties::LipidAccInfo const & sp, core::Size ) {
	lipid_acc_data_.push_back( sp );
}

utility::vector1< properties::LipidAccInfo >
MembraneConformation::lipid_acc_data() const {
	return lipid_acc_data_;
}


/// Access Methods for Membrane Info ///////////

/// @brief Get Membrane Center Coords
core::Vector
MembraneConformation::membrane_center() {
	return residue( membrane_ ).atom( 2 ).xyz();
}

/// @brief Get Membrane Normal Coords
core::Vector
MembraneConformation::membrane_normal() {
	return residue( membrane_ ).atom( 1 ).xyz();
}

/// @brief Get Membrane Thickness Parameter
core::Real
MembraneConformation::membrane_thickness() {
	return residue( membrane_ ).atom( 3 ).xyz().y();
}

/// @brief Get Chain Embedding Center Coords
core::Vector
MembraneConformation::embedding_center( core::Size chain ) {
	core::Size resnum = embres_map_[ chain ].second;
	return residue( resnum ).atom( 2 ).xyz();
}

/// @brief Get Chain Embedding Normal coords
core::Vector
MembraneConformation::embedding_normal( core::Size chain ) {
	core::Size resnum = embres_map_[ chain ].second;
	return residue( resnum ).atom( 1 ).xyz();
}

/// @brief Get Chain Embedding Depth Parameter
core::Real
MembraneConformation::embedding_depth( core::Size chain ) {
	core::Size resnum = embres_map_[ chain ].second;
	return residue( resnum ).atom( 3 ).xyz().y();
}

//// FoldTree ////////////////////////////

/// @brief Build Membrane FoldTree
void
MembraneConformation::build_membrane_foldtree() {

	using namespace core::kinematics;

	// Construct a new foldtree
	FoldTreeOP ft = new FoldTree();

	// Add Initial membrane edge and maintain a jump counter
	int jump_counter = 1;
	ft->add_edge( 1, membrane_, jump_counter );
	jump_counter++;

	// Add peptide edges
	for ( core::Size i = 1; i <= num_chains()-1; ++i ) {

		// Add edge for each chain
		ft->add_edge( chain_begin(i), chain_end(i), -1 );

		// Add a chain-connecting edge if between chains
		if ( i != num_chains()-1 ) {
			ft->add_edge( chain_begin(i), chain_begin(i+1), jump_counter );
			jump_counter++;
		}
	}

	// Set the membrane to the root of the foldtree. This must occur before adding embedding edges
	// because the embedding edges are added to a non-vertex posiiton and will be lost in the reorder arbitrarily
	ft->reorder(membrane_);

	// Add membrane edges
	for ( core::Size i = 1; i <= embres_map_.size(); ++i ) {

		// Add edge between each chain beginning and its corresponding embedding
		ft->add_edge( embres_map_[ i ].first, embres_map_[ i ].second, jump_counter );
		jump_counter++;
	}

	// Setup the new foldtree
	fold_tree( *ft );
}

/// @brief Override Insert Residue
void
MembraneConformation::insert_residue_by_jump(
		 Residue const &,
		 Size const, // desired seqpos of new_rsd
		 Size, // in the current sequence numbering, ie before insertion of seqpos
		 std::string const&, // could be ""
		 std::string const&, // ditto
		 bool
		 )
{
	using namespace core::membrane::util;
	throw new EXCN_NonMembrane("Cannot insert a residue into a membrane topology - throws membrane spanning topology out of sync!");
}

} // membrane
} // core

#endif // INCLUDED_core_membrane_MembraneConformation_cc
