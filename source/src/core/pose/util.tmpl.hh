// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/util.hh
/// @brief  Pose utilities
/// @author Phil Bradley
/// @author Modified by Sergey Lyskov

#ifndef INCLUDED_core_pose_util_TMPL_HH
#define INCLUDED_core_pose_util_TMPL_HH

// Unit headers
#include <core/pose/util.hh>

// C/C++ headers
#include <map>
#include <set>

// Utility headers
#include <numeric/xyzVector.hh>
#include <utility/vector1.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/scoring/ScoreType.hh>

// Package headers
#include <core/pose/Pose.hh>
#include <core/pose/MiniPose.fwd.hh>

#include <utility/vector1.hh>

namespace core {
namespace pose {


/// @brief Initialize a DOF_ID_Map for a given Pose using the DOF_ID_Map's current default fill values
template< typename T >
void
initialize_dof_id_map( id::DOF_ID_Map< T > & dof_map, Pose const & pose )
{
	Size const n_res( pose.size() );
	dof_map.clear();
	dof_map.resize( n_res );
	for ( Size i = 1; i <= n_res; ++i ) {
		dof_map.resize( i, pose.residue(i).natoms() );
	}
}


/// @brief Initialize a DOF_ID_Map for a given Pose using a specified fill value
template< typename T >
void
initialize_dof_id_map( id::DOF_ID_Map< T > & dof_map, Pose const & pose, T const & value )
{
	Size const n_res( pose.size() );
	dof_map.clear();
	dof_map.resize( n_res );
	for ( Size i = 1; i <= n_res; ++i ) {
		dof_map.resize( i, pose.residue(i).natoms(), value );
	}
}

/// @brief returns a Distance
core::Real
pose_max_nbr_radius( pose::Pose const & pose );

/// @brief Initialize an AtomID_Map for a given Pose using the AtomID_Map's current default fill values
template< typename T >
void
initialize_atomid_map( id::AtomID_Map< T > & atom_map, pose::Pose const & pose )
{
	Size const n_res( pose.size() );
	atom_map.clear();
	atom_map.resize( n_res );
	for ( Size i = 1; i <= n_res; ++i ) {
		atom_map.resize( i, pose.residue_type(i).natoms() );
	}
}


/// @brief Initialize an AtomID_Map for a given Pose using a specified fill value
template< typename T >
void
initialize_atomid_map( id::AtomID_Map< T > & atom_map, pose::Pose const & pose, T const & value )
{
	Size const n_res( pose.size() );
	atom_map.clear();
	atom_map.resize( n_res );
	for ( Size i = 1; i <= n_res; ++i ) {
		atom_map.resize( i, pose.residue_type(i).natoms(), value );
	}
}


/// @brief Initialize an AtomID_Map for a given Conformation using the AtomID_Map's current default fill values
template< typename T >
void
initialize_atomid_map( id::AtomID_Map< T > & atom_map, conformation::Conformation const & conformation )
{
	Size const n_res( conformation.size() );
	atom_map.clear();
	atom_map.resize( n_res );
	for ( Size i = 1; i <= n_res; ++i ) {
		atom_map.resize( i, conformation.residue_type(i).natoms() );
	}
}


/// @brief Initialize an AtomID_Map for a given Conformation using a specified fill value
template< typename T >
void
initialize_atomid_map( id::AtomID_Map< T > & atom_map, conformation::Conformation const & conformation, T const & value )
{
	Size const n_res( conformation.size() );
	atom_map.clear();
	atom_map.resize( n_res );
	for ( Size i = 1; i <= n_res; ++i ) {
		atom_map.resize( i, conformation.residue_type(i).natoms(), value );
	}
}

/// @brief Initialize an AtomID_Map for a given Pose using the AtomID_Map's current default fill values
template< typename T >
void
initialize_atomid_map_heavy_only( id::AtomID_Map< T > & atom_map, pose::Pose const & pose )
{
	Size const n_res( pose.size() );
	atom_map.clear();
	atom_map.resize( n_res );
	for ( Size i = 1; i <= n_res; ++i ) {
		atom_map.resize( i, pose.residue_type(i).nheavyatoms() );
	}
}


/// @brief Initialize an AtomID_Map for a given Pose using a specified fill value
template< typename T >
void
initialize_atomid_map_heavy_only( id::AtomID_Map< T > & atom_map, pose::Pose const & pose, T const & value )
{
	Size const n_res( pose.size() );
	atom_map.clear();
	atom_map.resize( n_res );
	for ( Size i = 1; i <= n_res; ++i ) {
		atom_map.resize( i, pose.residue_type(i).nheavyatoms(), value );
	}
}


/// @brief Initialize an AtomID_Map for a given Conformation using the AtomID_Map's current default fill values
template< typename T >
void
initialize_atomid_map_heavy_only( id::AtomID_Map< T > & atom_map, conformation::Conformation const & conformation )
{
	Size const n_res( conformation.size() );
	atom_map.clear();
	atom_map.resize( n_res );
	for ( Size i = 1; i <= n_res; ++i ) {
		atom_map.resize( i, conformation.residue_type(i).nheavyatoms() );
	}
}


/// @brief Initialize an AtomID_Map for a given Conformation using a specified fill value
template< typename T >
void
initialize_atomid_map_heavy_only( id::AtomID_Map< T > & atom_map, conformation::Conformation const & conformation, T const & value )
{
	Size const n_res( conformation.size() );
	atom_map.clear();
	atom_map.resize( n_res );
	for ( Size i = 1; i <= n_res; ++i ) {
		atom_map.resize( i, conformation.residue_type(i).nheavyatoms(), value );
	}
}


// PyRosetta concreate version of functions
inline void initialize_atomid_map_AtomID( id::AtomID_Map< id::AtomID > & atom_map, pose::Pose const & pose )
{ initialize_atomid_map(atom_map, pose); }

inline void initialize_atomid_map_AtomID( id::AtomID_Map< id::AtomID > & atom_map, pose::Pose const & pose, id::AtomID const & value )
{ initialize_atomid_map(atom_map, pose, value); }


inline void initialize_atomid_map_AtomID( id::AtomID_Map< id::AtomID > & atom_map, conformation::Conformation const & conformation )
{ initialize_atomid_map(atom_map, conformation); }

inline void initialize_atomid_map_AtomID( id::AtomID_Map< id::AtomID > & atom_map, conformation::Conformation const & conformation, id::AtomID const & value )
{ initialize_atomid_map(atom_map, conformation, value); }

} // pose
} // core

#endif // INCLUDED_core_pose_util_HH
