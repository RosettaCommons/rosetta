// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/forge/methods/pose_mod.hh
/// @brief methods for Pose modifications
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_methods_pose_mod_hh
#define INCLUDED_protocols_forge_methods_pose_mod_hh

// type headers
#include <core/types.hh>

// project headers

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/pose/Pose.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/MMAtomType.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResConnID.fwd.hh>
#include <core/chemical/ResConnID.hh>
#include <core/chemical/ResidueConnection.fwd.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/VariantType.fwd.hh>
//XRW_B_T1
//#include <core/coarse/Translator.fwd.hh>
//XRW_E_T1
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
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
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Map.fwd.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/types.hh>
#include <core/kinematics/AtomPointer.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/AtomWithDOFChange.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/Edge.fwd.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/ResidueCoordinateChangeList.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/util.hh>
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
#include <core/pack/dunbrack/RotamerLibrary.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <basic/MetricValue.fwd.hh>
// AUTO-REMOVED #include <basic/OStream.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector0_bool.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
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
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
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
#include <ObjexxFCL/FArray.all.fwd.hh>
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
#include <ObjexxFCL/Observer.fwd.hh>
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>
#include <ObjexxFCL/Star.fwd.hh>
#include <ObjexxFCL/Star.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <map>
#include <ostream>
#include <string>
#include <vector>
#include <boost/bind.hpp>
#include <boost/function.hpp>



namespace protocols {
namespace forge {
namespace methods {


/// @brief grow a series of residues to the left of a position
/// @param[in,out] pose Pose the pose to modify
/// @param[in] anchor the anchor position
/// @param[in] begin iterator that points to the first ResidueOP
/// @param[in] end iterator that points just beyond the last ResidueOP
/// @param[in] correct_terminus re-add lower terminus if found, default true
/// @param[in] use_existing_crd No idealization: place residues with existing
///  coordinates, default false.
/// @return the left endpoint of the new growth or zero if nothing to do
/// @remarks Use a reverse iterator to grow right -> left along vector containers.
template< typename ResidueOPIterator >
core::Size
grow_left_r(
	core::pose::Pose & pose,
	core::Size anchor,
	ResidueOPIterator begin,
	ResidueOPIterator end,
	bool const correct_terminus = true,
	bool const use_existing_crd = false
)
{
	using core::Size;

	if ( begin == end ) {
		return 0;
	}

	// query to see if position we're extending from is a true terminus
	bool const had_lower_terminus = pose.residue( anchor ).is_lower_terminus();

	// grow extension
	Size current_pos = anchor; // tracks the anchor as it moves
	for ( ResidueOPIterator i = begin; i != end; ++i, ++current_pos ) {
		pose.conformation().safely_prepend_polymer_residue_before_seqpos( **i, anchor, !use_existing_crd ); // will remove any terminus
	}

	// add terminus if necessary
	if ( correct_terminus && had_lower_terminus && !pose.residue( anchor ).is_lower_terminus() ) {
		core::pose::add_lower_terminus_type_to_pose_residue( pose, anchor );
	}

	return anchor; // the left endpoint of the new growth
}


/// @brief grow a series of residues to the left of a position
/// @param[in,out] pose Pose the pose to modify
/// @param[in] anchor the anchor position
/// @param[in] begin iterator that points to the first ResidueTypeOP
/// @param[in] end iterator that points just beyond the last ResidueTypeOP
/// @param[in] correct_terminus re-add lower terminus if found, default true
/// @return the left endpoint of the new growth or zero if nothing to do
/// @remarks Use a reverse iterator to grow right -> left along vector containers.
template< typename ResidueTypeOPIterator >
core::Size
grow_left_rtype(
	core::pose::Pose & pose,
	core::Size anchor,
	ResidueTypeOPIterator begin,
	ResidueTypeOPIterator end,
	bool const correct_terminus = true
)
{
	using core::conformation::ResidueOPs;
	using core::conformation::ResidueFactory;

	// create residues
	ResidueOPs residues;
	for ( ResidueTypeOPIterator i = begin; i != end; ++i ) {
		residues.push_back( ResidueFactory::create_residue( **i ) );
	}

	return grow_left_r( pose, anchor, residues.begin(), residues.end(), correct_terminus, false );
}


/// @brief grow a series of residues to the right of a position
/// @param[in,out] pose Pose the pose to modify
/// @param[in] anchor the anchor position, can be 0 if operating on empty
///  Pose
/// @param[in] begin iterator that points to the first ResidueOP
/// @param[in] end iterator that points just beyond the last ResidueOP
/// @param[in] correct_terminus re-add upper terminus if found, default true
/// @param[in] use_existing_crd No idealization: place residues with existing
///  coordinates, default false.
/// @return the right endpoint of the new growth or zero if nothing to do
template< typename ResidueOPIterator >
core::Size
grow_right_r(
	core::pose::Pose & pose,
	core::Size anchor,
	ResidueOPIterator begin,
	ResidueOPIterator end,
	bool const correct_terminus = true,
	bool const use_existing_crd = false
)
{
	using core::Size;

	if ( begin == end ) {
		return 0;
	}

	// query to see if position we're extending from is a true terminus
	bool const had_upper_terminus = anchor > 0 ? pose.residue( anchor ).is_upper_terminus() : false;

	// grow extension
	Size current_pos = anchor; // tracks the right endpoint of the new growth

	ResidueOPIterator i = begin;
	if ( anchor == 0 ) { // special case where pose is empty
		pose.append_residue_by_bond( **i );
		++i;
		++current_pos;
	}

	for ( ; i != end; ++i, ++current_pos ) {
		pose.conformation().safely_append_polymer_residue_after_seqpos( **i, current_pos, !use_existing_crd ); // will remove any terminus
	}

	// add terminus if necessary
	if ( correct_terminus && had_upper_terminus && !pose.residue( current_pos ).is_upper_terminus() ) {
		core::pose::add_upper_terminus_type_to_pose_residue( pose, current_pos );
	}

	return current_pos; // the right endpoint of the new growth
}


/// @brief grow a series of residues to the right of a position
/// @param[in,out] pose Pose the pose to modify
/// @param[in] anchor the anchor position, can be 0 if operating on empty
///  Pose
/// @param[in] begin iterator that points to the first ResidueTypeOP
/// @param[in] end iterator that points just beyond the last ResidueTypeOP
/// @param[in] correct_terminus re-add upper terminus if found, default true
/// @return the right endpoint of the new growth or zero if nothing to do
template< typename ResidueTypeOPIterator >
core::Size
grow_right_rtype(
	core::pose::Pose & pose,
	core::Size anchor,
	ResidueTypeOPIterator begin,
	ResidueTypeOPIterator end,
	bool const correct_terminus = true
)
{
	using core::conformation::ResidueOPs;
	using core::conformation::ResidueFactory;

	// create residues
	ResidueOPs residues;
	for ( ResidueTypeOPIterator i = begin; i != end; ++i ) {
		residues.push_back( ResidueFactory::create_residue( **i ) );
	}

	return grow_right_r( pose, anchor, residues.begin(), residues.end(), correct_terminus, false );
}


/// @brief add cutpoint variants at a specific position
/// @param[in,out] Pose to modify
/// @param[in] position at which to add cutpoint variants
/// @return true if cutpoint variants added, false if position not a topological cutpoint
bool
add_cutpoint_variants(
	core::pose::Pose & pose,
	core::Size const pos
);


/// @brief remove cutpoint variants at a specific position
/// @param[in,out] Pose to modify
/// @param[in] position at which to remove cutpoint variants
/// @return true if cutpoint variants removed, false if no cutpoint variants found
///  or position not a topological cutpoint
bool
remove_cutpoint_variants(
	core::pose::Pose & pose,
	core::Size const pos
);


/// @brief restore residues (i.e. sidechains)
/// @param[in] old2new map indicating residues to be transferred and
///  the mapping from archive_pose position -> pose position
/// @param[in] archive_pose the original Pose to take residues from
/// @param[out] pose the altered Pose that needs residue restoration
void
restore_residues(
	std::map< core::Size, core::Size > const & old2new,
	core::pose::Pose & archive_pose,
	core::pose::Pose & pose
);


/// @brief restore residues (i.e. sidechains)
/// @param[in] archive_pose the original Pose to take residues from
/// @param[out] pose the altered Pose that needs residue restoration
/// @remarks length between two poses must be equal
void
restore_residues(
	core::pose::Pose & archive_pose,
	core::pose::Pose & pose
);


} // methods
} // forge
} // protocols


#endif /* INCLUDED_protocols_forge_methods_pose_mod_HH */
