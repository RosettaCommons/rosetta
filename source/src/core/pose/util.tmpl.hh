// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
// AUTO-REMOVED #include <basic/datacache/BasicDataCache.hh>
#include <numeric/xyzVector.hh>
#include <utility/vector1.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/VariantType.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Map.hh>
// AUTO-REMOVED #include <core/id/DOF_ID_Map.hh>
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

//Auto Headers
#include <platform/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/MMAtomType.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResConnID.fwd.hh>
#include <core/chemical/ResConnID.hh>
#include <core/chemical/ResidueConnection.fwd.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/orbitals/ICoorOrbitalData.hh>
#include <core/chemical/orbitals/OrbitalType.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/orbitals/OrbitalXYZCoords.hh>
#include <core/conformation/signals/ConnectionEvent.fwd.hh>
#include <core/conformation/signals/ConnectionEvent.hh>
#include <core/conformation/signals/GeneralEvent.fwd.hh>
#include <core/conformation/signals/GeneralEvent.hh>
#include <core/conformation/signals/IdentityEvent.fwd.hh>
#include <core/conformation/signals/IdentityEvent.hh>
#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/conformation/signals/LengthEvent.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/conformation/signals/XYZEvent.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Map.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/types.hh>
#include <core/kinematics/AtomPointer.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/AtomWithDOFChange.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/Edge.fwd.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/ResidueCoordinateChangeList.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/pose/PDB_Info.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector0_bool.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/keys/Key2Tuple.fwd.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/keys/Key3Tuple.fwd.hh>
#include <utility/keys/Key3Tuple.hh>
#include <utility/keys/Key4Tuple.fwd.hh>
#include <utility/keys/Key4Tuple.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/signals/BufferedSignalHub.fwd.hh>
#include <utility/signals/BufferedSignalHub.hh>
#include <utility/signals/Link.fwd.hh>
#include <utility/signals/Link.hh>
#include <utility/signals/LinkUnit.fwd.hh>
#include <utility/signals/LinkUnit.hh>
#include <utility/signals/PausableSignalHub.fwd.hh>
#include <utility/signals/PausableSignalHub.hh>
#include <utility/signals/SignalHub.fwd.hh>
#include <utility/signals/SignalHub.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/internal/ColPointers.hh>
#include <numeric/internal/ColVectors.hh>
#include <numeric/internal/ColsPointer.hh>
#include <numeric/internal/RowPointers.hh>
#include <numeric/internal/RowVectors.hh>
#include <numeric/internal/RowsPointer.hh>
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/IndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/InitializerSentinel.hh>
#include <ObjexxFCL/Observer.fwd.hh>
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/ProxySentinel.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>
#include <ObjexxFCL/Star.fwd.hh>
#include <ObjexxFCL/Star.hh>
#include <ObjexxFCL/proxy_const_assert.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <ostream>
#include <string>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>



namespace core {
namespace pose {


/// @brief Initialize a DOF_ID_Map for a given Pose using the DOF_ID_Map's current default fill values
template< typename T >
void
initialize_dof_id_map( id::DOF_ID_Map< T > & dof_map, Pose const & pose )
{
	Size const n_res( pose.n_residue() );
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
	Size const n_res( pose.n_residue() );
	dof_map.clear();
	dof_map.resize( n_res );
	for ( Size i = 1; i <= n_res; ++i ) {
		dof_map.resize( i, pose.residue(i).natoms(), value );
	}
}

///@brief returns a Distance
core::Real
pose_max_nbr_radius( pose::Pose const & pose );

/// @brief Initialize an AtomID_Map for a given Pose using the AtomID_Map's current default fill values
template< typename T >
void
initialize_atomid_map( id::AtomID_Map< T > & atom_map, pose::Pose const & pose )
{
	Size const n_res( pose.n_residue() );
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
	Size const n_res( pose.n_residue() );
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
	Size const n_res( pose.n_residue() );
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
	Size const n_res( pose.n_residue() );
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
