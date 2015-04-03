// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/rms/rms_util.tmpl.hh
/// @brief  RMS stuff from rosetta++ that requires templates
/// @author James Thompson
/// @author Ian Davis
/// @date   Wed Aug 22 12:10:37 2007

#ifndef INCLUDED_core_scoring_rms_util_tmpl_HH
#define INCLUDED_core_scoring_rms_util_tmpl_HH

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/rms_util.hh>

#include <numeric/model_quality/rms.hh>
#include <numeric/xyzVector.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1A.hh>
#include <ObjexxFCL/FArray2A.hh>
#include <list>

//Auto Headers
#include <platform/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/MMAtomType.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResConnID.fwd.hh>
#include <core/chemical/ResConnID.hh>
#include <core/chemical/ResidueConnection.fwd.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/types.hh>
#include <core/chemical/orbitals/ICoorOrbitalData.hh>
#include <core/chemical/orbitals/OrbitalType.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Conformation.hh>
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
#include <core/conformation/symmetry/SymDof.fwd.hh>
#include <core/conformation/symmetry/SymSlideInfo.fwd.hh>
#include <core/conformation/symmetry/SymSlideInfo.hh>
#include <core/conformation/symmetry/SymmData.fwd.hh>
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Map.fwd.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.hh>
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
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/ResidueCoordinateChangeList.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/pose/PDBInfo.fwd.hh>
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
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/assert.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/stream_util.hh>
#include <utility/string_util.hh>
#include <utility/vector0_bool.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
#include <utility/keys/AutoKey.fwd.hh>
#include <utility/keys/AutoKey.hh>
#include <utility/keys/Key.fwd.hh>
#include <utility/keys/Key.hh>
#include <utility/keys/Key2Tuple.fwd.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/keys/Key3Tuple.fwd.hh>
#include <utility/keys/Key3Tuple.hh>
#include <utility/keys/Key4Tuple.fwd.hh>
#include <utility/keys/Key4Tuple.hh>
#include <utility/keys/KeyLess.fwd.hh>
#include <utility/keys/KeyLookup.fwd.hh>
#include <utility/keys/KeyLookup.hh>
#include <utility/keys/NoClient.fwd.hh>
#include <utility/keys/NoClient.hh>
#include <utility/keys/SmallKeyVector.fwd.hh>
#include <utility/keys/SmallKeyVector.hh>
#include <utility/keys/UserKey.fwd.hh>
#include <utility/keys/VariantKey.fwd.hh>
#include <utility/keys/VariantKey.hh>
#include <utility/options/AnyOption.fwd.hh>
#include <utility/options/AnyOption.hh>
#include <utility/options/AnyVectorOption.fwd.hh>
#include <utility/options/AnyVectorOption.hh>
#include <utility/options/BooleanOption.fwd.hh>
#include <utility/options/BooleanOption.hh>
#include <utility/options/BooleanVectorOption.fwd.hh>
#include <utility/options/BooleanVectorOption.hh>
#include <utility/options/FileOption.fwd.hh>
#include <utility/options/FileOption.hh>
#include <utility/options/FileVectorOption.fwd.hh>
#include <utility/options/FileVectorOption.hh>
#include <utility/options/IntegerOption.fwd.hh>
#include <utility/options/IntegerOption.hh>
#include <utility/options/IntegerVectorOption.fwd.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <utility/options/Option.fwd.hh>
#include <utility/options/Option.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/PathOption.fwd.hh>
#include <utility/options/PathOption.hh>
#include <utility/options/PathVectorOption.fwd.hh>
#include <utility/options/PathVectorOption.hh>
#include <utility/options/RealOption.fwd.hh>
#include <utility/options/RealOption.hh>
#include <utility/options/RealVectorOption.fwd.hh>
#include <utility/options/RealVectorOption.hh>
#include <utility/options/ScalarOption.fwd.hh>
#include <utility/options/ScalarOption.hh>
#include <utility/options/ScalarOption_T_.fwd.hh>
#include <utility/options/ScalarOption_T_.hh>
#include <utility/options/StringOption.fwd.hh>
#include <utility/options/StringOption.hh>
#include <utility/options/StringVectorOption.fwd.hh>
#include <utility/options/StringVectorOption.hh>
#include <utility/options/VariantOption.fwd.hh>
#include <utility/options/VariantOption.hh>
#include <utility/options/VectorOption.fwd.hh>
#include <utility/options/VectorOption.hh>
#include <utility/options/VectorOption_T_.fwd.hh>
#include <utility/options/VectorOption_T_.hh>
#include <utility/options/mpi_stderr.hh>
#include <utility/options/keys/AnyOptionKey.fwd.hh>
#include <utility/options/keys/AnyOptionKey.hh>
#include <utility/options/keys/AnyVectorOptionKey.fwd.hh>
#include <utility/options/keys/AnyVectorOptionKey.hh>
#include <utility/options/keys/BooleanOptionKey.fwd.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/options/keys/BooleanVectorOptionKey.fwd.hh>
#include <utility/options/keys/BooleanVectorOptionKey.hh>
#include <utility/options/keys/FileOptionKey.fwd.hh>
#include <utility/options/keys/FileOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.fwd.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>
#include <utility/options/keys/IntegerOptionKey.fwd.hh>
#include <utility/options/keys/IntegerOptionKey.hh>
#include <utility/options/keys/IntegerVectorOptionKey.fwd.hh>
#include <utility/options/keys/IntegerVectorOptionKey.hh>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/keys/OptionKeys.hh>
#include <utility/options/keys/PathOptionKey.fwd.hh>
#include <utility/options/keys/PathOptionKey.hh>
#include <utility/options/keys/PathVectorOptionKey.fwd.hh>
#include <utility/options/keys/PathVectorOptionKey.hh>
#include <utility/options/keys/RealOptionKey.fwd.hh>
#include <utility/options/keys/RealOptionKey.hh>
#include <utility/options/keys/RealVectorOptionKey.fwd.hh>
#include <utility/options/keys/RealVectorOptionKey.hh>
#include <utility/options/keys/ScalarOptionKey.fwd.hh>
#include <utility/options/keys/ScalarOptionKey.hh>
#include <utility/options/keys/StringOptionKey.fwd.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/options/keys/StringVectorOptionKey.fwd.hh>
#include <utility/options/keys/StringVectorOptionKey.hh>
#include <utility/options/keys/VectorOptionKey.fwd.hh>
#include <utility/options/keys/VectorOptionKey.hh>
#include <utility/options/keys/all.hh>
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
#include <numeric/HomogeneousTransform.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/constants.hh>
#include <numeric/conversions.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/sphericalVector.hh>
#include <numeric/trig.functions.hh>
#include <numeric/types.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyz.functions.hh>
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
#include <ObjexxFCL/FArray.all.fwd.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArray1.all.fwd.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray1A.fwd.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1P.fwd.hh>
#include <ObjexxFCL/FArray1P.hh>
#include <ObjexxFCL/FArray2.all.fwd.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2A.fwd.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2P.fwd.hh>
#include <ObjexxFCL/FArray2P.hh>
#include <ObjexxFCL/FArray3.all.fwd.hh>
#include <ObjexxFCL/FArray3.fwd.hh>
#include <ObjexxFCL/FArray3A.fwd.hh>
#include <ObjexxFCL/FArray3D.fwd.hh>
#include <ObjexxFCL/FArray3P.fwd.hh>
#include <ObjexxFCL/FArray4.all.fwd.hh>
#include <ObjexxFCL/FArray4.fwd.hh>
#include <ObjexxFCL/FArray4A.fwd.hh>
#include <ObjexxFCL/FArray4D.fwd.hh>
#include <ObjexxFCL/FArray4P.fwd.hh>
#include <ObjexxFCL/FArray5.all.fwd.hh>
#include <ObjexxFCL/FArray5.fwd.hh>
#include <ObjexxFCL/FArray5A.fwd.hh>
#include <ObjexxFCL/FArray5D.fwd.hh>
#include <ObjexxFCL/FArray5P.fwd.hh>
#include <ObjexxFCL/FArray6.all.fwd.hh>
#include <ObjexxFCL/FArray6.fwd.hh>
#include <ObjexxFCL/FArray6A.fwd.hh>
#include <ObjexxFCL/FArray6D.fwd.hh>
#include <ObjexxFCL/FArray6P.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/IndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/InitializerSentinel.hh>
#include <ObjexxFCL/KeyFArray1D.fwd.hh>
#include <ObjexxFCL/KeyFArray2D.fwd.hh>
#include <ObjexxFCL/KeyFArray3D.fwd.hh>
#include <ObjexxFCL/KeyFArray4D.fwd.hh>
#include <ObjexxFCL/KeyFArray5D.fwd.hh>
#include <ObjexxFCL/KeyFArray6D.fwd.hh>
#include <ObjexxFCL/Observer.fwd.hh>
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/ProxySentinel.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>
#include <ObjexxFCL/Star.fwd.hh>
#include <ObjexxFCL/Star.hh>
#include <ObjexxFCL/StaticIndexRange.fwd.hh>
#include <ObjexxFCL/StaticIndexRange.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/proxy_const_assert.hh>
#include <ObjexxFCL/string.functions.hh>
#include <algorithm>
#include <utility/assert.hh>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option.hh>
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>
#include <boost/unordered/unordered_map.hpp>

//Auto Headers

namespace core {
namespace scoring {

template< class T >
core::Real
rmsd_with_super(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	std::list< core::Size > const & subset_residues,
	T* predicate
)
{
	ASSERT_ONLY(core::Size const nres1 = pose1.total_residue();)
	ASSERT_ONLY(core::Size const nres2 = pose2.total_residue();)
debug_assert( nres1 == nres2 );

	std::vector< core::Vector > p1_coords;
	std::vector< core::Vector > p2_coords;

	BOOST_FOREACH(Size i, subset_residues){
	    Size num_atoms ( pose1.residue(i).natoms() );
	    if ( predicate == is_ligand_heavyatom ||
		 predicate == is_polymer_heavyatom ||
		 predicate == is_heavyatom ) {
	      //		debug_assert( pose1.residue(i).nheavyatoms() == pose2.residue(i).nheavyatoms() );
	    } else if ( num_atoms > pose2.residue(i).natoms() ){
	      num_atoms = pose2.residue(i).natoms();
	    }
	    for ( core::Size j = 1; j <= num_atoms; ++j ) {
	      if( j <= pose1.residue(i).natoms() &&
		  j <= pose2.residue(i).natoms() ) {
		if ( predicate( pose1, pose2, i, j ) ) {
		  p1_coords.push_back( pose1.residue(i).xyz(j) );
		  p2_coords.push_back( pose2.residue(i).xyz(j) );
		}
	      }
	    }
	  }
debug_assert( p1_coords.size() == p2_coords.size() );

	int const natoms = p1_coords.size();
	ObjexxFCL::FArray2D< core::Real > p1a( 3, natoms );
	ObjexxFCL::FArray2D< core::Real > p2a( 3, natoms );
	for ( int i = 0; i < natoms; ++i ) {
	  for ( int k = 0; k < 3; ++k ) { // k = X, Y and Z
	    p1a(k+1,i+1) = p1_coords[i][k];
	    p2a(k+1,i+1) = p2_coords[i][k];
	  }
	}

	return numeric::model_quality::rms_wrapper( natoms, p1a, p2a );
}

/// @brief Select atoms for RMS via a predicate function/functor.
/// @details Calculates minimal rms, allowing rotation/translation for best fit.
/// Parameter "predicate" should be a function pointer
/// [ or class that defines operator() ] with the following signature:
///
///   bool my_pred_func(Pose const & pose1, Pose const & pose2, Size resno, Size atomno);
///
/// It should return true if the atom should be included and false otherwise.
///
/// Example of use, to calculate C-alpha RMSD:
///   rmsd_with_super(pose1, pose2, is_protein_CA);
template< class T >
core::Real
rmsd_with_super(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	T* predicate
)
{
	std::list< core::Size > all_residues;
	for (core::Size n = 1; n <= pose1.total_residue(); n++ ){
		all_residues.push_back( n );
	}
	return rmsd_with_super( pose1, pose2, all_residues, predicate );
}


/// @brief Select a subset atoms for RMS via a predicate function/functor.
/// @details Calculates minimal rms, allowing rotation/translation for best fit.
///		Same as above function, but allows a subset of residues over which to
///		superposition to be passed in
/// Parameter "predicate" should be a function pointer
/// [ or class that defines operator() ] with the following signature:
///
///   bool my_pred_func(Pose const & pose1, Pose const & pose2, Size resno, Size atomno);
///
/// It should return true if the atom should be included and false otherwise.
///
/// Example of use, to calculate C-alpha RMSD:
///   rmsd_with_super(pose1, pose2, is_protein_CA, subset);
template< class T >
core::Real
rmsd_with_super_subset(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	ObjexxFCL::FArray1D_bool const & subset,
	T* predicate
)
{
	core::Size const nres1 = pose1.total_residue();
	//core::Size const nres2 = pose2.total_residue();
//debug_assert( nres1 == nres2 );

	std::vector< core::Vector > p1_coords;
	std::vector< core::Vector > p2_coords;

	for ( core::Size i = 1; i <= nres1; ++i ) {
		if ( subset( i ) ) {
			Size num_atoms ( pose1.residue(i).natoms() );
			if ( predicate == is_ligand_heavyatom ||
				predicate == is_polymer_heavyatom ||
				predicate == is_heavyatom ) {
			debug_assert( pose1.residue(i).natoms() == pose2.residue(i).natoms() );
			} else if ( num_atoms > pose2.residue(i).natoms() ){
				num_atoms = pose2.residue(i).natoms();
			}
			for ( core::Size j = 1; j <= num_atoms; ++j ) {
				if ( predicate( pose1, pose2, i, j ) ) {
					p1_coords.push_back( pose1.residue(i).xyz(j) );
					p2_coords.push_back( pose2.residue(i).xyz(j) );
				}
			}
		}
	}
debug_assert( p1_coords.size() == p2_coords.size() );

	int const natoms = p1_coords.size();
	ObjexxFCL::FArray2D< core::Real > p1a( 3, natoms );
	ObjexxFCL::FArray2D< core::Real > p2a( 3, natoms );
	for ( int i = 0; i < natoms; ++i ) {
		for ( int k = 0; k < 3; ++k ) { // k = X, Y and Z
			p1a(k+1,i+1) = p1_coords[i][k];
			p2a(k+1,i+1) = p2_coords[i][k];
		}
	}

	return numeric::model_quality::rms_wrapper( natoms, p1a, p2a );
}


/// @brief Select atoms for RMS via a predicate function/functor.
/// @details Calculates minimal rms, NOT allowing rotation/translation --
/// uses current coordinates as-is.
/// Parameter "predicate" should be a function pointer
/// [ or class that defines operator() ] with the following signature:
///
///   bool my_pred_func(Pose const & pose1, Pose const & pose2, Size resno, Size atomno);
///
/// It should return true if the atom should be included and false otherwise.
///
/// Example of use, to calculate C-alpha RMSD:
///   rmsd_no_super(pose1, pose2, is_protein_CA);
template< class T >
core::Real
rmsd_no_super(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	T* predicate
)
{
	core::Size const nres1 = pose1.total_residue();
	ASSERT_ONLY(core::Size const nres2 = pose2.total_residue();)
debug_assert( nres1 == nres2 );

	core::Real sum2( 0.0 );
	core::Size natoms( 0 );
	for ( core::Size i = 1; i <= nres1; ++i ) {
		Size num_atoms ( pose1.residue(i).natoms() );
		if ( predicate == is_ligand_heavyatom ||
			predicate == is_polymer_heavyatom ||
			predicate == is_heavyatom ) {
		debug_assert( pose1.residue(i).natoms() == pose2.residue(i).natoms() );
		} else if ( num_atoms > pose2.residue(i).natoms() ){
			num_atoms = pose2.residue(i).natoms();
		}
		for ( core::Size j = 1; j <= num_atoms; ++j ) {
			if ( predicate( pose1, pose2, i, j ) ) {
				core::Vector diff = pose1.residue(i).xyz(j) - pose2.residue(i).xyz(j);
				sum2 += diff.length_squared();
				natoms += 1;
			}
		}
	}
	return std::sqrt(sum2 / natoms);
}

template< class T >
core::Real
rmsd_no_super(
	core::conformation::ResidueCOPs const & residues1,
	core::conformation::ResidueCOPs const & residues2,
	T* predicate
)
{
	core::Size const nres1 = residues1.size();
	ASSERT_ONLY(core::Size const nres2 = residues2.size();)
debug_assert( nres1 == nres2 );

	core::Real sum2( 0.0 );
	core::Size natoms( 0 );
	for ( core::Size i = 1; i <= nres1; ++i ) {
		Size num_atoms ( residues1[i]->natoms() );
		if (  predicate == is_ligand_heavyatom_residues
			/* predicate == is_polymer_heavyatom || */
			/* predicate == is_heavyatom */) { // for Rotate.cc
		debug_assert( residues1[i]->natoms() == residues2[i]->natoms() );
		} else if ( num_atoms > residues2[i]->natoms() ){
			num_atoms = residues2[i]->natoms();
		}
		for ( core::Size j = 1; j <= num_atoms; ++j ) {
			if ( predicate( *residues1[i], *residues2[i], j ) ) {
				core::Vector diff = residues1[i]->xyz(j) - residues2[i]->xyz(j);
				sum2 += diff.length_squared();
				natoms += 1;
			}
		}
	}
	return std::sqrt(sum2 / natoms);
}

/// @brief Select atoms for RMS via a predicate function/functor.
/// @details Calculates minimal rms over a subset of residues,
/// NOT allowing rotation/translation --
/// uses current coordinates as-is.
/// Parameter "predicate" should be a function pointer
/// [ or class that defines operator() ] with the following signature:
///
///   bool my_pred_func(Pose const & pose1, Pose const & pose2, Size resno, Size atomno);
/// the subset is a vector of all the residues, with true for those
/// over which to calculate the rms
///
/// It should return true if the atom should be included and false otherwise.
///
/// Example of use, to calculate C-alpha RMSD:
///   rmsd_no_super_subset(pose1, pose2, subset, is_protein_CA);
template< class T >
core::Real
rmsd_no_super_subset(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	ObjexxFCL::FArray1D_bool const & subset,
	T* predicate
)
{
	core::Size const nres1 = pose1.total_residue();
	ASSERT_ONLY(core::Size const nres2 = pose2.total_residue();)
debug_assert( nres1 == nres2 );

	core::Real sum2( 0.0 );
	core::Size natoms( 0 );
	for ( core::Size i = 1; i <= nres1; ++i ) {
	if ( subset( i ) ) {
		Size num_atoms ( pose1.residue(i).natoms() );
			if ( predicate == is_ligand_heavyatom ||
				predicate == is_polymer_heavyatom ||
				predicate == is_heavyatom ) {
			debug_assert( pose1.residue(i).natoms() == pose2.residue(i).natoms() );
			} else if ( num_atoms > pose2.residue(i).natoms() ){
				num_atoms = pose2.residue(i).natoms();
			}
			for ( core::Size j = 1; j <= num_atoms; ++j ) {
				if ( predicate( pose1, pose2, i, j ) ) {
					core::Vector diff = pose1.residue(i).xyz(j) - pose2.residue(i).xyz(j);
					sum2 += diff.length_squared();
					natoms += 1;
				}
			}
		}
	}
	return std::sqrt(sum2 / natoms);
}

/// @brief like function above, but uses sequence mapping,
/// i.e. sections of poses of different lengths can be compared
/// at the moment sorta assumes that residues at corresponding
/// positions have the same identity, mainly becaue the predicates
/// are structured that way...
template< class T >
core::Real
rmsd_no_super_subset(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	ObjexxFCL::FArray1D_bool const & subset,
	core::id::SequenceMapping const & seqmap,
	T* predicate
)
{
  core::Size const nres1 = pose1.total_residue();
  core::Real sum2( 0.0 );
  core::Size natoms( 0 );

  for ( core::Size i = 1; i <= nres1; ++i ) {
    if ( subset( i ) ) {
      core::Size res2 = seqmap[i];
      runtime_assert(res2 != 0 );
      Size num_atoms ( pose1.residue(i).natoms() );
      if ( predicate == is_ligand_heavyatom ||
	predicate == is_polymer_heavyatom ||
	predicate == is_heavyatom ) {
debug_assert( pose1.residue(i).natoms() == pose2.residue(res2).natoms() );
      } else if ( num_atoms > pose2.residue(res2).natoms() ){
	num_atoms = pose2.residue(res2).natoms();
      }
      for ( core::Size j = 1; j <= num_atoms; ++j ) {
	if ( predicate( pose1, pose2, i, j ) ) {
		core::Vector diff = pose1.residue(i).xyz(j) - pose2.residue(res2).xyz(j);
		sum2 += diff.length_squared();
		natoms += 1;
	}
      }
    }
  }
  return std::sqrt(sum2 / natoms);
}


/// @brief function to return the biggest deviation between an atom in a pair of poses,
template< class T >
core::Real
biggest_residue_deviation_no_super(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	T* predicate
){

	core::Size const nres1 = pose1.total_residue();
	ASSERT_ONLY(core::Size const nres2 = pose2.total_residue();)
debug_assert( nres1 == nres2 );

	core::Real biggest_dev2( 0.0 );
	for ( core::Size i = 1; i <= nres1; ++i ) {

		Real res_dev2(0.0);
		Size res_count_atoms(0);
		Size num_atoms ( pose1.residue(i).natoms() );
		if ( predicate == is_ligand_heavyatom ||
			predicate == is_polymer_heavyatom ||
			predicate == is_heavyatom ) {
		debug_assert( pose1.residue(i).natoms() == pose2.residue(i).natoms() );
		} else if ( num_atoms > pose2.residue(i).natoms() ){
			num_atoms = pose2.residue(i).natoms();
		}
		for ( core::Size j = 1; j <= num_atoms; ++j ) {
			if ( predicate( pose1, pose2, i, j ) ) {
				core::Vector diff = pose1.residue(i).xyz(j) - pose2.residue(i).xyz(j);
				res_dev2 += diff.length_squared();
				res_count_atoms++;
			}
		}
		res_dev2 /= res_count_atoms;
		if( res_dev2 > biggest_dev2 ) biggest_dev2 = res_dev2;
	}
	return std::sqrt( biggest_dev2 );
}


/// @brief function to return the biggest deviation between an atom in a pair of poses,
/// @brief as specified by the predicate and the subset
template< class T >
core::Real
biggest_residue_deviation_no_super_subset(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	ObjexxFCL::FArray1D_bool const & subset,
	T* predicate
){

	core::Size const nres1 = pose1.total_residue();
	ASSERT_ONLY(core::Size const nres2 = pose2.total_residue());
debug_assert( nres1 == nres2 );

	core::Real biggest_dev2( 0.0 );
	for ( core::Size i = 1; i <= nres1; ++i ) {
	if ( subset( i ) ) {
		Real res_dev2(0.0);
		Size res_count_atoms(0);
		Size num_atoms ( pose1.residue(i).natoms() );
			if ( predicate == is_ligand_heavyatom ||
				predicate == is_polymer_heavyatom ||
				predicate == is_heavyatom ) {
			debug_assert( pose1.residue(i).natoms() == pose2.residue(i).natoms() );
			} else if ( num_atoms > pose2.residue(i).natoms() ){
				num_atoms = pose2.residue(i).natoms();
			}
			for ( core::Size j = 1; j <= num_atoms; ++j ) {
				if ( predicate( pose1, pose2, i, j ) ) {
					core::Vector diff = pose1.residue(i).xyz(j) - pose2.residue(i).xyz(j);
					res_dev2 += diff.length_squared();
					res_count_atoms++;
				}
			}
			res_dev2 /= res_count_atoms;
			if( res_dev2 > biggest_dev2 ) biggest_dev2 = res_dev2;
		}
	}
	return std::sqrt( biggest_dev2 );
}

template< class T >
void
fill_rmsd_coordinates(
	int & natoms,
	ObjexxFCL::FArray2D< core::Real > & p1a,
	ObjexxFCL::FArray2D< core::Real > & p2a,
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	T* predicate
)
{
	core::Size const nres1 = pose1.total_residue();
	core::Size const nres2 = pose2.total_residue();
//debug_assert( nres1 == nres2 );

	std::vector< core::Vector > p1_coords;
	std::vector< core::Vector > p2_coords;

	for ( core::Size i = 1; i <= std::min( nres1, nres2 ); ++i ) {
		if( pose1.residue(i).is_virtual_residue() || pose2.residue(i).is_virtual_residue() ) continue;
	//debug_assert( pose1.residue(i).natoms() == pose2.residue(i).natoms() );
		for ( core::Size j = 1; j <= pose1.residue(i).natoms(); ++j ) {
			if ( (*predicate)( pose1, pose2, i, j ) ) {
				p1_coords.push_back( pose1.residue(i).xyz(j) );
				p2_coords.push_back( pose2.residue(i).xyz(j) );
			}
		}
	}
debug_assert( p1_coords.size() == p2_coords.size() );

	natoms = p1_coords.size();
	p1a.dimension( 3, natoms );
	p2a.dimension( 3, natoms );
	for ( int i = 0; i < natoms; ++i ) {
		for ( int k = 0; k < 3; ++k ) { // k = X, Y and Z
			p1a(k+1,i+1) = p1_coords[i][k];
			p2a(k+1,i+1) = p2_coords[i][k];
		}
	}
}

/// @brief Select a subset atoms for Symmetric RMS via a predicate function/functor.
/// Example of use, to calculate C-alpha RMSD:
///   rmsd_with_super(pose1, pose2, is_protein_CA, subset);
template< class T >
core::Real
sym_rmsd_with_super_subset(
	core::pose::Pose const & native_pose,
	core::pose::Pose const & pose2,
	ObjexxFCL::FArray1D_bool const & subset,
	T* predicate
)
{
	core::Size const nres1 = native_pose.total_residue();
	//core::Size const nres2 = pose2.total_residue();
//debug_assert( nres1 == nres2 );

	std::vector< core::Vector > p1_coords;
	std::vector< core::Vector > p2_coords;

	for ( core::Size i = 1; i <= nres1; ++i ) {
		if ( subset( i ) ) {
			Size num_atoms ( native_pose.residue(i).natoms() );
			if ( predicate == is_ligand_heavyatom ||
				predicate == is_polymer_heavyatom ||
				predicate == is_heavyatom ) {
			debug_assert( native_pose.residue(i).natoms() == pose2.residue(i).natoms() );
			} else if ( num_atoms > pose2.residue(i).natoms() ){
				num_atoms = pose2.residue(i).natoms();
			}
			for ( core::Size j = 1; j <= num_atoms; ++j ) {
				if ( predicate( native_pose, pose2, i, j ) ) {
					p1_coords.push_back( native_pose.residue(i).xyz(j) );
					p2_coords.push_back( pose2.residue(i).xyz(j) );
				}
			}
		}
	}
debug_assert( p1_coords.size() == p2_coords.size() );

	int const natoms = p1_coords.size();
	ObjexxFCL::FArray2D< core::Real > p1a( 3, natoms );
	ObjexxFCL::FArray2D< core::Real > p2a( 3, natoms );
	for ( int i = 0; i < natoms; ++i ) {
		for ( int k = 0; k < 3; ++k ) { // k = X, Y and Z
			p1a(k+1,i+1) = p1_coords[i][k];
			p2a(k+1,i+1) = p2_coords[i][k];
		}
	}


debug_assert( dynamic_cast< core::conformation::symmetry::SymmetricConformation const * >( &pose2.conformation() ) );

	core::conformation::symmetry::SymmetricConformation const & symm_conf (
		dynamic_cast<core::conformation::symmetry::SymmetricConformation const & > ( pose2.conformation()) );
	core::conformation::symmetry::SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

	int const N ( symm_info->subunits() );
	int const nres ( symm_info->num_total_residues_without_pseudo() );
	ObjexxFCL::FArray2D< core::Real > p1a_shuffle( 3, nres );

	core::Real rms = 1e3; //Since fast_rms has not been evaluated yet
	int const num_atoms_subunit (natoms);
	std::vector< std::vector<int> > shuffle_map;
	create_shuffle_map_recursive_rms(std::vector<int>(), N,shuffle_map);
	for (int j=1; j < int (shuffle_map.size()); j++ ){
		for (int i=0; i < N; ++i ) {
			int const begin ( shuffle_map.at(j).at(i)*num_atoms_subunit*3);
			for ( int k = 0; k < num_atoms_subunit*3; ++k ) {
				int const begin_shuffled (i*num_atoms_subunit*3);
					p1a_shuffle[begin_shuffled+k] = p1a[begin+k];
			}
		}
		Real rms_shuffle = numeric::model_quality::rms_wrapper( natoms, p1a_shuffle, p2a );
		if ( rms_shuffle < rms ) {
			rms = rms_shuffle;
		}
	}
	return rms;
}

} // end namespace scoring
} // end namespace core

#endif
