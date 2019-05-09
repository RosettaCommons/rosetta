// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file     core/conformation/membrane/MembraneInfo.hh
/// @brief    Data describing the relationship between protein(s) and a membrane environment
///
/// @details  MembraneInfo is a container object that describes membrane-protein relationships
///             1. Coordinates of the membrane
///             2. A pointer to MEM which describes relative orientation (Residue)
///             3. Topology of the transmembrane spans (SpanningTopology)
///             4. Physical and chemical properties of the implicit lipid membrane (ImplicitLipidInfo)
///
/// @note     This object is a member of Conformation and should only be accessed using
///           pose.conformation().membrane_info(). Do not access MEM outside of the framework!
///
/// @author   Rebecca Alford (ralford3@jhu.edu)

#ifndef INCLUDED_core_conformation_membrane_MembraneInfo_hh
#define INCLUDED_core_conformation_membrane_MembraneInfo_hh

// Unit headers
#include <core/conformation/membrane/MembraneInfo.fwd.hh>

// Package Headers
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/membrane/SpanningTopology.fwd.hh>
#include <core/conformation/membrane/ImplicitLipidInfo.fwd.hh>

// Project Headers
#include <core/types.hh>

#include <core/kinematics/FoldTree.hh>

// Utility Headers
#include <numeric/MathMatrix.hh>
#include <utility/pointer/ReferenceCount.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace conformation {
namespace membrane {

/// @brief Data describing the relationship between protein(s) and a membrane environment
class MembraneInfo : public utility::pointer::ReferenceCount {

public: // Constructors & Setup

	/// @brief Create MembraneInfo from initialized data
	/// @details Creates a MembraneInfo object by linking the conformation
	/// to the pose, specify the  membrane residue number, membrane jump number,
	/// spanning topology object and optional lipophilicity data. Thickness and
	/// steepness are currently constants
	MembraneInfo(
		core::Size membrane_pos,
		core::SSize membrane_jump,
		core::Size membrane_core,
		core::Real thickness,
		core::Real steepness,
		SpanningTopologyOP topology
	);

	/// @brief Create MembraneInfo from initialized data
	/// @details Creates a MembraneInfo object by linking the conformation
	/// to the pose, specify the  membrane residue number, membrane jump number,
	/// spanning topology object and optional lipophilicity data. Thickness and
	/// steepness are currently constants
	MembraneInfo(
		core::Size membrane_pos,
		core::SSize membrane_jump,
		core::Real steepness,
		SpanningTopologyOP topology,
		std::string lipid_composition_name,
		core::Real lipid_composition_temp
	);

	/// @brief Create a deep copy of all data in this object.
	MembraneInfo( MembraneInfo const & src );

	/// @brief create a deep copy of all data in thsi object upon assignment
	MembraneInfo &
	operator=( MembraneInfo const & src );

	/// @brief Destructor
	~MembraneInfo();

	/// @brief Generate a string representation of information represented by this MembraneInfo and send it to std::cout
	virtual void show() const;

	/// @brief Generate a string representation of information represented by ths MembraneInfo
	virtual void show( std::ostream & output ) const;

public: // Chemical information about this membrane

	/// @brief Effective thickness of the membrane (default = 15)
	virtual core::Real membrane_thickness() const;

	/// @brief Steepness of hydrophobic -> hydrophillic transition (defualt = 10)
	virtual core::Real membrane_steepness() const;

	/// @brief core membrane thickness
	core::Real membrane_core() const;

public: // membrane position & orientation

	/// @brief Membrane center
	Vector membrane_center( Conformation const & conf ) const;

	/// @brief Membrane normal
	Vector membrane_normal( Conformation const & conf ) const;

	/// @brief Is residue in the membrane?
	bool in_membrane( Conformation const & conf, core::Size resnum ) const;

	/// @brief Compute residue position relative to membrane normal
	Real
	residue_z_position( Conformation const & conf, core::Size resnum ) const;

	/// @brief Compute atom position relative to membrane normal
	Real
	atom_z_position( Conformation const & conf, core::Size resnum, core::Size atomnum ) const;

	/// @brief Sequence position of the membrane residue
	core::Size
	membrane_rsd_num() const;

	/// @brief Indeitifier for the membrane jump
	core::SSize
	membrane_jump() const;

	/// @brief Allow a protocol to set a new jump number for the membrane jump
	void
	set_membrane_jump( core::SSize jumpnum );

	/// @brief Somewhat weak check that a membrane foldtree is valid
	bool check_membrane_fold_tree( core::kinematics::FoldTree const & ft_in ) const;

public: // topology of TM spans

	/// @brief Transmembrane spaning topology
	SpanningTopologyOP spanning_topology() const;

public: // implicit lipid information

	/// @brief Get implicit lipid information
	ImplicitLipidInfoOP implicit_lipids() const;

private: // default constructor

	/// @brief Create a default version of MembraneInfo (DONT USE)
	MembraneInfo();

private: // data

	// membrane thickness & transition steepness
	core::Real thickness_;
	core::Real steepness_;

	core::Real membrane_core_;

	// membrane residue number in the pose
	core::Size membrane_rsd_num_;

	// membrane jump position
	core::SSize membrane_jump_;

	// Spanning Topology information
	SpanningTopologyOP spanning_topology_;

	// Implicit lipid information
	ImplicitLipidInfoOP implicit_lipids_;

#ifdef    SERIALIZATION
public:
	friend class cereal::access;
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // MembraneInfo

/// @brief Show MembraneInfo method for pyrosetta
std::ostream & operator << ( std::ostream & os, MembraneInfo const & mem_info );

} // membrane
} // conformation
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_conformation_membrane_MembraneInfo )
#endif // SERIALIZATION


#endif // INCLUDED_core_conformation_membrane_MembraneInfo_hh



