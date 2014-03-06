// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 	 core/membrane/MembraneConformation.hh
///
/// @brief 	 Membrane Conformation
/// @details The Membrane Conformation is responsible for:
///             - maintaining a correct membrane foldtree
///             - maintaining references to the membrane and embedding residues
///             - providing access to membrane related data
///
/// @note    Last Modified 1/9/14
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_MembraneConformation_hh
#define INCLUDED_core_membrane_MembraneConformation_hh

// Unit headers
#include <core/membrane/MembraneConformation.fwd.hh>
#include <core/conformation/Conformation.hh>

// Project Headers
#include <core/membrane/properties/SpanningTopology.hh> 
#include <core/membrane/properties/LipidAccInfo.hh>

// Package Headers
#include <core/kinematics/FoldTree.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/vector1.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>

using namespace core::conformation;

namespace core {
namespace membrane {

/// @brief A Membrane conformation: Additional data for maintaining a mebrane
/// @details Handles membrane proteins
class MembraneConformation : public core::conformation::Conformation {

public: // construction, copying, access, invariants

	//// Constructors //////////////////////////

	/// @brief Default Constructor
	MembraneConformation();

	/// @brief Standard Constructor
	MembraneConformation(
						 Conformation const & conf,
						 utility::vector1< std::pair< int, int > > embres_map,
						 int membrane
						 );

	/// @brief Copy Constructor
	MembraneConformation( MembraneConformation const & src );

	/// @brief Assignemnt Operator
	Conformation &
	operator=( MembraneConformation const & src );

	/// @brief Clone Method
	ConformationOP
	clone() const;

	/// @brief Return Membrane Embedding Map
	utility::vector1< std::pair< int, int > > embres_map() const;

	/// @brief Return Membrane root
	int membrane() const;

	/// @brief Assert Membrane Conformation Invariants
	bool is_membrane();

public: // membrane info

	/// Access Methods for Membrane Info ///////////

	/// @brief Get Membrane Center Coords
	core::Vector membrane_center();

	/// @brief Get Membrane Normal Coords
	core::Vector membrane_normal();

	/// @brief Get Membrane Thickness Parameter
	core::Real membrane_thickness();

	/// @brief Get Chain Embedding Center Coords
	core::Vector embedding_center( core::Size chain );

	/// @brief Get Chain Embedding Normal coords
	core::Vector embedding_normal( core::Size chain );

	/// @brief Get Chain Embedding Depth Parameter
	core::Real embedding_depth( core::Size chain );
	
	/// @brief Return the total number of polymer reisdues in the pose (no mp residues)
	core::Size total_polymer_residue();
	
	/// @brief Return the total number of polymer chains in the pose (no mp chain)
	core::Size num_polymer_chains();

	///// Overrided Methods ///////////
	void insert_residue_by_jump(
		 Residue const & new_rsd_in,
		 Size const seqpos, // desired seqpos of new_rsd
		 Size anchor_pos, // in the current sequence numbering, ie before insertion of seqpos
		 std::string const& anchor_atom, // could be ""
		 std::string const& root_atom, // ditto
		 bool new_chain
								);

public: // membrane conformation data


	/// @brief Add Spanning Topology
	void add_topology_by_chain( properties::SpanningTopology sp, core::Size chain );
	utility::vector1< properties::SpanningTopology > spanning_topology() const;

	/// @brief Add Lipid Accessibility Info
	void add_lips_by_chain( properties::LipidAccInfo const & sp, core::Size chain );
	utility::vector1< properties::LipidAccInfo > lipid_acc_data() const;

private: // methods

	//// FoldTree ////////////////////////////

	/// @brief Build Membrane FoldTree
	void build_membrane_foldtree();

private: // data

	// Store Embedding Residue Map
	utility::vector1< std::pair< int, int > > embres_map_;

	// Store index of the membrane origin
	int membrane_;

	// Lipid accessibility data
	utility::vector1< properties::LipidAccInfo > lipid_acc_data_;
	utility::vector1< properties::SpanningTopology > spanning_topology_;

}; // MembraneConformation

} // membrane
} // core

#endif // INCLUDED_core_membrane_MembraneConformation_hh

