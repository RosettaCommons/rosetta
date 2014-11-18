// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		core/conformation/membrane/MembraneInfo.hh
///
/// @brief		MembraneInfo - Membrane Pose Definition Object
/// @details	The membrane info object is responsible for storing non-coordinate
///				derived information, extending the traditional definiton of a pose
///				to describe a membrane protein in Rosetta. This information includes:
///				 - the resiue number of the membrane virtual residue containing
///					the posiiton of the membrane
///				 - membrane spanning topology
///				 - membrane lipophilicity info (user-specified)
///
///				This object belongs to the conformation and should be accessed
///				via pose.conformation().membrane_info().
///
///	     		Last Modified: 7/25/14
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_conformation_membrane_MembraneInfo_hh
#define INCLUDED_core_conformation_membrane_MembraneInfo_hh

// Unit headers
#include <core/conformation/membrane/MembraneInfo.fwd.hh>

// Project Headers
#include <core/conformation/Conformation.fwd.hh>

#include <core/conformation/membrane/SpanningTopology.fwd.hh>
#include <core/conformation/membrane/LipidAccInfo.fwd.hh>

// Package Headers
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
using namespace core::conformation::membrane;
using namespace core::kinematics;
	
/// @brief Membrane Info Object - Contains information for describing
/// a membrane protein in Rosetta, including the position of the membrane virtual
/// residue, spanning topology, and lips acc data.
class MembraneInfo : public utility::pointer::ReferenceCount {
		
public:

	////////////////////
	/// Constructors ///
	////////////////////
	
	/// @brief Default Constructor
	/// @details Creates a default copy of the membrane info object
	/// with the membrane virtual at nres = 1, thicnkess = 15A, steepenss
	/// of fullatom membrane transition = 10A and no spanning topology or
	/// lips info defined
	MembraneInfo();
	
	/// @brief Custom Constructor - Membrane pos & topology
	/// @details Creates a default copy of the membrane info object
	/// with the membrane virtual at nres = membrane_pos, thicnkess = 15A, steepenss
	/// of fullatom membrane transition = 10A and spanning topology defined
	/// by API provided vector1 of spanning topology obejcts (by chain)
	MembraneInfo(
		Conformation & conformation,
		core::Size membrane_pos,
		SpanningTopologyOP topology,
		core::SSize membrane_jump = 2
	);
	
	/// @brief Custom Constructor - Membrane pos, topology & lips
	/// @details Creates a default copy of the membrane info object
	/// with the membrane virtual at nres = membrane_pos, thicnkess = 15A, steepenss
	/// of fullatom membrane transition = 10A, spanning topology defined
	/// by API provided vector1 of spanning topology obejcts (by chain), and
	/// lipid accessibility defined by a vector1 of lipid acc objects.
	MembraneInfo(
		Conformation & conformation,
		core::Size membrane_pos,
		SpanningTopologyOP topology,
		LipidAccInfoOP lips,
		core::SSize membrane_jump = 2	);
		
	/// @brief Copy Constructor
	/// @details Create a deep copy of this object
	MembraneInfo( MembraneInfo const & src );
	
	/// @brief Assignment Operator
	/// @details Create a deep copy of this object while overloading the assignemnt
	/// operator "="
	MembraneInfo &
	operator=( MembraneInfo const & src );
	
	/// @brief Destructor
	~MembraneInfo();
	
	/// @brief  Generate string representation of MembraneInfo for debugging purposes.
	virtual void show(std::ostream & output=std::cout) const;
	
	////////////////////////////////////////
	/// Coordinate Derived Membrane Info ///
	////////////////////////////////////////
	
	/// @brief Return the center coordinate of the membrane
	/// @details Return the center xyz coordinate of the membrane described in the
	/// MPct atom of the membrane virtual residue
	Vector membrane_center() const;
	
	/// @brief Returns the normal of the membrane
	/// @details Returns the normal (direction) of the membrane described in the MPnm
	/// atom of the membrane virtual residue.
	Vector membrane_normal() const;
	
	/// @brief Compute Residue Z Position relative to mem
	/// @details Compute the z position of a residue relative to the pre-defined
	/// layers in the membrane. Maintians the relative coordinate frame
	Real
	residue_z_position( core::Size resnum ) const;
	
	/// @brief Compute atom Z Position relative to mem
	/// @details Compute the z position of an atom relative to the pre-defined
	/// layers in the membrane. Maintians the relative coordinate frame
	Real
	atom_z_position( core::Size resnum, core::Size atomnum ) const;

	////////////////////////////////////////////
	/// Non-Coordinate Derived Membrane Data ///
	////////////////////////////////////////////
	
	/// @brief Return position of membrane residue
	/// @details Return the residue number of the membrane virtual
	/// residue in the pose
	core::Size membrane_rsd_num() const;

	/// @brief Return membrane thickness
	/// @details Return the membrane thicnkess used by the
	/// fullatom energy method (centroid is hard coded for now)
	core::Real membrane_thickness() const;
	
	/// @brief Return the membrane transition steepness
	/// @details Return the steepness betwen isotropic and ansitropic
	/// layers of the membrane used by the fullatom energy methods
	core::Real membrane_steepness() const;
	
	/// @brief Return a list of membrane spanning topology objects
	/// @details Return a vector1 of spanning topology objects defining
	/// the starting and ending position of membrane spans per chain.
	// TODO: should be const reference?
	SpanningTopologyOP spanning_topology();
	
	/// @brief Return a list of lipid accessibility objects
	/// @details Return a vector1 of lipid accessibility info objects
	/// describing lipid exposre of individual residues in the pose
	LipidAccInfoOP lipid_acc_data() const;
	
	////////////////////////////////////
	/// Membrane Base Fold Tree Info ///
	////////////////////////////////////
	
	/// @brief Get the number of the membrane jump
	/// @details Get a core::SSize (int) denoting the number of the fold tree jump
	/// relating the membrane residue to the rest of the pose
	core::SSize membrane_jump() const;
		
	/// @brief	 Check membrane fold tree
	/// @details Check that the membrane jump num is a jump point and checking
	///			 for a reasonable fold tree
	bool check_membrane_fold_tree( FoldTree const & ft_in ) const;
	
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
	
} // membrane
} // conformation
} // core

#endif // INCLUDED_core_membrane_MembraneInfo_hh

