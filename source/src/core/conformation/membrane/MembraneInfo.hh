// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  core/conformation/membrane/MembraneInfo.hh
///
/// @brief  Information about the membrane bilayer and its relationship with the protein(s)
/// @details MembraneInfo is responsible for describing attributes of the
///    membrane bilayer *and* the position & orientation of the bilayer
///    in 3D space. All of the data members work together to accomplish
///    this representation. And the players are:
///     = A pointer to the Pose Conformation, which includes an MEM residue
///     = Topology of transmembrane spans
///     = Per-residue lipophilicity
///     = Membrane thickness & steepness - derived from chemical profiles
///
///    This object is a member of Conformation and should be accessed by
///    pose.conofrmation().membrane_info(). DO NOT access the MEM residue
///    outside of the framework!!!
///
///    Last Updated: 7/23/15
/// @author  Rebecca Faye Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_conformation_membrane_MembraneInfo_hh
#define INCLUDED_core_conformation_membrane_MembraneInfo_hh

// Unit headers
#include <core/conformation/membrane/MembraneInfo.fwd.hh>

// Package Headers
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/membrane/SpanningTopology.fwd.hh>
#include <core/conformation/membrane/LipidAccInfo.fwd.hh>

// Project Headers
#include <core/types.hh>

#include <core/kinematics/FoldTree.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <cstdlib>

namespace core {
namespace conformation {
namespace membrane {

using namespace core::conformation;
using namespace core::kinematics;

/// @brief MembraneInfo describes the membrane bilayer and its relationship with the protein
class MembraneInfo : public utility::pointer::ReferenceCount {

public: // Constructors & Setup

	/// @brief Create a default version of MembraneInfo (DONT USE)
	/// @details Initializes all data members to dummy or empty values
	/// Use the fully specified constructors instead. MembraneInfo is a
	/// data container but is NOT responsible for initialization.
	MembraneInfo();

	/// @brief Create MembraneInfo from initialized data
	/// @details Creates a MembraneInfo object by linking the conformation
	/// to the pose, specify the  membrane residue number, membrane jump number,
	/// spanning topology object and optional lipophilicity data. Thickness and
	/// steepness are currently constants
	MembraneInfo(
		Conformation & conformation,
		core::Size membrane_pos,
		core::SSize membrane_jump,
		SpanningTopologyOP topology
	);

	/// @brief Create MembraneInfo from initialized data with lipophilicity
	/// @details Creates a MembraneInfo object by linking the conformation
	/// to the pose, specify the  membrane residue number, membrane jump number,
	/// spanning topology object and optional lipophilicity data. Thickness and
	/// steepness are currently constants
	MembraneInfo(
		Conformation & conformation,
		core::Size membrane_pos,
		core::SSize membrane_jump,
		LipidAccInfoOP lips,
		SpanningTopologyOP topology
	);

	/// @brief Create a deep copy of all data in this object.
	MembraneInfo( MembraneInfo const & src );

	/// @brief create a deep copy of all data in thsi object upon assignment
	MembraneInfo &
	operator=( MembraneInfo const & src );

	/// @brief Destructor
	~MembraneInfo();

	/// @brief Generate a string representation of information represented by ths MembraneInfo
	virtual void show( std::ostream & output=std::cout ) const;

public: // Chemical information about this membrane

	/// @brief Effective thickness of the membrane
	/// @details For IMM = default is 15.
	core::Real membrane_thickness() const;

	/// @brief Steepness of hydrophobic -> hydrophillic transition
	/// @details For IMM - default is 10
	core::Real membrane_steepness() const;

public: // membrane position & orientation

	/// @brief Membrane center
	/// @details Returns the xyzVector describing the center of the membrane
	/// This is the same as the MPct atom of the membrane (MEM) residue.
	Vector membrane_center() const;

	/// @brief Membrane normal
	/// @details Returns the membrane normal, which describes the membrane
	/// orientation. This is the same as the xyzVector in the MPnm atom
	/// in the membrane residue.
	Vector membrane_normal() const;

	/// @brief Is residue in the membrane? Takes CA coordinate
	/// @details Uses the thickness stored in MembraneInfon and the residue_z_position
	bool in_membrane( core::Size resnum ) const;

	/// @brief Compute residue position relative to membrane normal
	/// @details Calculate the z coordinate of the residue, projected onto
	/// the membrane normal axis. Objective is to maintain correct coordinates
	/// in relative coordinate frame.
	Real
	residue_z_position( core::Size resnum ) const;

	/// @brief Compute atom position relative to membrane normal
	/// @details Calculate the z coordinate of the atom, projected onto
	/// the membrane normal axis. Objective is to maintain correct coordinates
	/// in relative coordinate frame.
	Real
	atom_z_position( core::Size resnum, core::Size atomnum ) const;

	/// @brief Sequence position of the membrane residue
	/// @details Return the residue number of MEM (rsd.seqpos()) in the pose
	core::Size membrane_rsd_num() const;

	/// @brief Indeitifier for the membrane jump
	/// @details Returns an integer (core::Size) denoting the jump number in the foldtree
	/// representing the jump relating the membrane residue to the rest of the molecule
	core::SSize membrane_jump() const;

	/// @brief Allow a protocol to set a new jump number for the membrane jump
	/// @details Set the membrane jump number (core::SSize)
	void
	set_membrane_jump( core::SSize jumpnum );

	/// @brief Somewhat weak check that a membrane foldtree is valid. Use checks in
	/// protocols/membrane/util.hh instead!
	bool check_membrane_fold_tree( FoldTree const & ft_in ) const;

public: // topology of TM spans and lipophilicity

	/// @brief Transmembrane spaning topology
	/// @details Return a SpanningTopology object, which includes a
	/// list of Span objects, describing the start and end sequence
	/// positions of each transmembrane span
	SpanningTopologyOP spanning_topology() const;

	/// @brief Does this MembraneInfo includes lipophilicity information?
	bool
	include_lips() const;

	/// @brief Per-residue lipophilicity (probability of exposure to lipid)
	/// @details Returns a LipidAccInfo describing per residue probability
	/// of exposure to lipid. Data calcualted via the run_lips.pl script
	/// and provided by the user on the commandline if applicable
	LipidAccInfoOP lipid_acc_data() const;

private: // data

	// Keep track of the Pose's conformation
	Conformation& conformation_;

	// Fullatom constants
	core::Real thickness_;
	core::Real steepness_;

	// membrane residue number in the pose
	core::Size membrane_rsd_num_;

	// membrane jump position
	core::SSize membrane_jump_;

	// Lipit Accessibility and Topology Info
	LipidAccInfoOP lipid_acc_data_;
	SpanningTopologyOP spanning_topology_;

}; // MembraneInfo

/// @brief Show MembraneInfo method for pyrosetta
std::ostream & operator << ( std::ostream & os, MembraneInfo const & mem_info );

} // membrane
} // conformation
} // core

#endif // INCLUDED_core_conformation_membrane_MembraneInfo_hh

