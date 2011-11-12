// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constraints/AtomPairConstraint.cc
///
/// @brief
/// @author Oliver Lange

// Unit Headers
#include <core/scoring/constraints/AtomPairConstraint.hh>

// Package Headers
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/FuncFactory.hh>

// Project Headers
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/kinematics/ShortestPathInFoldTree.hh>
// Utility Headers
#include <basic/Tracer.hh>

// Numeric Headers
#include <numeric/deriv/distance_deriv.hh>

#include <platform/types.hh>
#include <core/types.hh>
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
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/VariantType.fwd.hh>
#include <core/chemical/types.hh>
#include <core/chemical/orbitals/ICoorOrbitalData.hh>
#include <core/chemical/orbitals/OrbitalType.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <core/chemical/sdf/MolData.fwd.hh>
#include <core/chemical/sdf/MolData.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Residue.hh>
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
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Map.fwd.hh>
#include <core/id/DOF_ID_Map.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
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
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/ResidueCoordinateChangeList.fwd.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/pose/MiniPose.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/AtomPairConstraint.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintCreator.fwd.hh>
#include <core/scoring/constraints/ConstraintFactory.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/Func.fwd.hh>
#include <core/scoring/constraints/Func.hh>
#include <core/scoring/constraints/FuncFactory.fwd.hh>
#include <core/scoring/constraints/HarmonicFunc.fwd.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/XYZ_Func.fwd.hh>
#include <core/scoring/constraints/XYZ_Func.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
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
#include <utility/factory/WidgetRegistrator.hh>
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
#include <numeric/xyzVector.hh>
#include <numeric/internal/ColPointers.hh>
#include <numeric/internal/ColVectors.hh>
#include <numeric/internal/ColsPointer.hh>
#include <numeric/internal/RowPointers.hh>
#include <numeric/internal/RowVectors.hh>
#include <numeric/internal/RowsPointer.hh>
#include <numeric/random/random.fwd.hh>
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
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/proxy_const_assert.hh>
#include <ObjexxFCL/string.functions.hh>
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
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableData.fwd.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/DataCache.fwd.hh>
#include <basic/datacache/DataCache.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option.hh>
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>


static basic::Tracer tr("core.io.constraints");


namespace core {
namespace scoring {
namespace constraints {

using core::pose::atom_id_to_named_atom_id;
using core::pose::named_atom_id_to_atom_id;

///
void
AtomPairConstraint::score( XYZ_Func const & xyz, EnergyMap const &, EnergyMap & emap ) const
{
	Real const score_val =  score( xyz( atom1_ ), xyz( atom2_ ) );
	emap[ this->score_type() ] += score_val;
}

///
Real
AtomPairConstraint::score(
													Vector const & xyz1,
													Vector const & xyz2
													) const
{
	//	std::cout << "score " <<  atom1_ << ' ' << atom2_ << ' ' << xyz1.distance( xyz2 ) << " --> " << func( xyz1.distance( xyz2 ) ) << std::endl;;
	return func( xyz1.distance( xyz2 ) );
}

/// @brief Copies the data from this Constraint into a new object and returns an OP
/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
/// if a sequence_mapping is present it is used to map residue numbers .. NULL = identity mapping
/// to the new object. Intended to be implemented by derived classes.
ConstraintOP AtomPairConstraint::remapped_clone( pose::Pose const& src, pose::Pose const& dest, id::SequenceMappingCOP smap ) const {
	id::NamedAtomID atom1( atom_id_to_named_atom_id( atom(1), src ) );
	id::NamedAtomID atom2( atom_id_to_named_atom_id( atom(2), src ) );

	if ( smap ) {
		atom1.rsd() = (*smap)[ atom1_.rsd() ];
		atom2.rsd() = (*smap)[ atom2_.rsd() ];
	}

	//get AtomIDs for target pose
	id::AtomID id1( named_atom_id_to_atom_id( atom1, dest ) );
	id::AtomID id2( named_atom_id_to_atom_id( atom2, dest ) );
	if ( id1.valid() && id2.valid() ) {
		return new AtomPairConstraint( id1, id2, func_, score_type() );
	} else {
		return NULL;
	}
}

bool
AtomPairConstraint::operator == ( Constraint const & other_cst ) const
{
	if( !dynamic_cast< AtomPairConstraint const * > ( &other_cst ) ) return false;

	AtomPairConstraint const & other( static_cast< AtomPairConstraint const & > (other_cst) );

	if( atom1_ != other.atom1_ ) return false;
	if( atom2_ != other.atom2_ ) return false;
	if( func_ != other.func_ ) return false;
	if( this->score_type() != other.score_type() ) return false;

	return true;
}

void AtomPairConstraint::show( std::ostream& out ) const {
	out << "AtomPairConstraint ("
			<< atom1_.atomno() << "," << atom1_.rsd() << "-"
			<< atom2_.atomno() << "," << atom2_.rsd() << ")" << std::endl;
	func_->show( out );
}

void AtomPairConstraint::show_def( std::ostream& out, pose::Pose const& pose ) const {
	out << type() << " " << atom_id_to_named_atom_id( atom1_, pose ) << " " << atom_id_to_named_atom_id( atom2_, pose ) << " ";
	func_->show_definition( out );
}

Real
AtomPairConstraint::dist( pose::Pose const & pose ) const {
	return dist( pose.conformation() );
}

Real AtomPairConstraint::dist( conformation::Conformation  const& conformation ) const {
		conformation::Residue const& res1( conformation.residue( atom1_.rsd() ) );
		conformation::Residue const& res2( conformation.residue( atom2_.rsd() ) );
		bool fail( false );
		if ( atom1_.atomno() == 0 || atom1_.atomno() > res1.natoms() ) {
			std::cerr << "AtomPairConstraint: atom1 out of bounds" << atom1_ << std::endl;
			fail = true;
		}
		if ( atom2_.atomno() == 0 || atom2_.atomno() > res2.natoms() ) {
			std::cerr << "AtomPairConstraint: atom2 out of bounds" << atom2_ << std::endl;
			fail = true;
		}
		assert( conformation.atom_tree().has( atom1_ ) );
		assert( conformation.atom_tree().has( atom2_ ) );

		if ( !conformation.atom_tree().has( atom1_ ) ) {
			std::cerr << "AtomPairConstraint: cannot find atom " << atom1_ << std::endl;
		}
		if ( !conformation.atom_tree().has( atom2_ ) ) {
			std::cerr << "AtomPairConstraint: cannot find atom " << atom2_ << std::endl;
		}
	  assert( !fail );
	Vector const & xyz1( conformation.xyz( atom1_ ) ), xyz2( conformation.xyz( atom2_ ) );
	Vector const f2( xyz1 - xyz2 );
	Real const dist( f2.length() );
	return dist;
}

Real AtomPairConstraint::dist( XYZ_Func const & xyz ) const {
	assert( atom1_.atomno() );
	assert( atom2_.atomno() );
	return xyz( atom1_ ).distance( xyz( atom2_ ) );
}

Size AtomPairConstraint::show_violations(
	std::ostream& out,
	pose::Pose const& pose,
	Size verbose_level,
	Real threshold
) const {

	if ( verbose_level > 80 ) {
		out << "AtomPairConstraint ("
			<< pose.residue_type(atom1_.rsd() ).atom_name( atom1_.atomno() ) << ":"
			<< atom1_.atomno() << "," << atom1_.rsd() << "-"
			<< pose.residue_type(atom2_.rsd() ).atom_name( atom2_.atomno() ) << ":"
			<< atom2_.atomno() << "," << atom2_.rsd() << ") ";
	}
	if ( verbose_level > 120 ) { //don't ask but I had a really weird bug to track down!
		conformation::Conformation const & conformation( pose.conformation() );
		Vector const & xyz1( conformation.xyz( atom1_ ) ), xyz2( conformation.xyz( atom2_ ) );
		out << "\ncoords1: " << xyz1[ 1 ] << " " << xyz1[ 2 ] << " " << xyz1[ 3 ] << " --- ";
		out << "coords1: " << xyz2[ 1 ] << " " << xyz2[ 2 ] << " " << xyz2[ 3 ] << "\n";
	}

	return func_->show_violations( out, dist( pose ), verbose_level, threshold );
}

// atom deriv
void
AtomPairConstraint::fill_f1_f2(
	AtomID const & atom,
	XYZ_Func const & xyz,
	Vector & F1,
	Vector & F2,
	EnergyMap const & weights
) const
{
	AtomID other_atom;
	if ( atom == atom1_ ) other_atom = atom2_;
	else if ( atom == atom2_ ) other_atom = atom1_;
	else {
		// std::cout << "Error in AtomPairConstraint::fill_f1_f2()" << std::endl;
		// std::cout << "Bad AtomID: (" << atom.rsd() << ", " << atom.atomno() << ") -- options: (";
		// std::cout << atom1_.rsd() << ", " << atom1_.atomno() << ") and (";
		// std::cout << atom2_.rsd() << ", " << atom2_.atomno() << ");" << std::endl;
		return;
	}

	//Vector const & xyz1( conformation.xyz( atom ) ), xyz2( conformation.xyz( other_atom ) );

	/*
	Vector const f2( xyz1 - xyz2 );
	Real const dist( f2.length() ), deriv( dfunc( dist ) );
	if ( deriv != 0.0 && dist != 0.0 ) {
	Vector const f1( xyz1.cross( xyz2 ) );
		F1 += ( ( deriv / dist ) * f1 ) * weights[ this->score_type() ];
		F2 += ( ( deriv / dist ) * f2 ) * weights[ this->score_type() ];
	}*/

	Real dist(0.0);
	Vector f1(0.0), f2(0.0);
	numeric::deriv::distance_f1_f2_deriv( xyz( atom ), xyz( other_atom ), dist, f1, f2 );
	Real wderiv( weights[ this->score_type() ] * dfunc( dist ));
	F1 += wderiv * f1;
	F2 += wderiv * f2;

}

ConstraintOP
AtomPairConstraint::remap_resid( core::id::SequenceMapping const &seqmap ) const
{
  if ( seqmap[atom1_.rsd()] != 0 && seqmap[atom2_.rsd()] != 0 ) {
    AtomID remap_a1( atom1_.atomno(), seqmap[atom1_.rsd()] ),
      remap_a2( atom2_.atomno(), seqmap[atom2_.rsd()] );
    return ConstraintOP( new AtomPairConstraint( remap_a1, remap_a2, this->func_ ) );
  } else {
    return NULL;
  }
}

///@details one line definition "AtomPairs atom1 res1 atom2 res2 function_type function_definition"
void
AtomPairConstraint::read_def(
	std::istream & data,
	core::pose::Pose const & pose,
	FuncFactory const & func_factory
) {
	Size res1, res2;
	std::string tempres1, tempres2;
	std::string name1, name2;
	std::string func_type;
	std::string type;

	data
		>> name1 >> tempres1
		>> name2 >> tempres2
		>> func_type;

	ConstraintIO::parse_residue( pose, tempres1, res1 );
	ConstraintIO::parse_residue( pose, tempres2, res2 );

	tr.Debug << "read: " << name1 << " " << name2 << " " << res1 << " " << res2 << " func: " << func_type << std::endl;
	if ( res1 > pose.total_residue() || res2 > pose.total_residue() ) {
		tr.Warning 	<< "ignored constraint (requested residue numbers exceed numbers of residues in pose): " << "Total in Pose: " << pose.total_residue() << " "
			<< name1 << " " << name2 << " " << res1 << " " << res2 << std::endl;
		data.setstate( std::ios_base::failbit );
		return;
	}

	if ( name1 == "H" && res1 == 1 && pose.is_fullatom() ) name1 = "1H";
	if ( name2 == "H" && res2 == 1 && pose.is_fullatom() ) name2 = "1H";

	atom1_ = named_atom_id_to_atom_id( id::NamedAtomID( name1, res1 ), pose );
	atom2_ = named_atom_id_to_atom_id( id::NamedAtomID( name2, res2 ), pose );

	if ( atom1_.atomno() == 0 || atom2_.atomno() == 0 ) {
		tr.Warning << "Error reading atoms: read in atom names("
			<< name1 << "," << name2 << "), "
			<< "and found AtomIDs (" << atom1_ << "," << atom2_ << ")" << std::endl;
			data.setstate( std::ios_base::failbit );
			runtime_assert( false );
			return;
	}

	func_ = func_factory.new_func( func_type );
	func_->read_data( data );

	if ( data.good() ) {
	//chu skip the rest of line since this is a single line defintion.
		while( data.good() && (data.get() != '\n') ) {}
		if ( !data.good() ) data.setstate( std::ios_base::eofbit );
	}

	if ( tr.Debug.visible() ) {
		func_->show_definition( tr.Debug );
		tr.Debug << std::endl;
	}
	runtime_assert( atom1_.valid() && atom2_.valid() );
} // read_def

core::Size
AtomPairConstraint::effective_sequence_separation( core::kinematics::ShortestPathInFoldTree const& sp ) const {
	return sp.dist( atom1_.rsd(), atom2_.rsd() );
}

} // constraints
} // scoring
} // core
