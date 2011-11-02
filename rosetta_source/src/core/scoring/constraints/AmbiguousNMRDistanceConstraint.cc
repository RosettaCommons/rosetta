// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constraints/AmbiguousNMRDistanceConstraint.cc
///
/// @brief
/// @author Oliver Lange

// Unit Headers
#include <core/scoring/constraints/AmbiguousNMRDistanceConstraint.hh>

// Package Headers
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/FuncFactory.hh>
#include <core/scoring/constraints/XYZ_Func.hh>
#include <core/kinematics/ShortestPathInFoldTree.hh>
// Project Headers
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/chemical/AA.hh>
// Utility Headers
#include <basic/Tracer.hh>
#include <ObjexxFCL/string.functions.hh>
#include <basic/prof.hh>
// Numeric Headers
#include <numeric/deriv/distance_deriv.hh>

#include <platform/types.hh>
#include <core/types.hh>
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
#include <core/scoring/constraints/AmbiguousNMRDistanceConstraint.fwd.hh>
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
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <execinfo.h>
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
#include <time.h>
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
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>


static basic::Tracer tr("core.io.constraints");


namespace core {
namespace scoring {
namespace constraints {


/// @brief Copies the data from this Constraint into a new object and returns an OP
/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
/// if a sequence_mapping is present it is used to map residue numbers .. NULL = identity mapping
/// to the new object. Intended to be implemented by derived classes.
ConstraintOP AmbiguousNMRDistanceConstraint::remapped_clone( pose::Pose const& src, pose::Pose const& dest, id::SequenceMappingCOP smap ) const {
	Atoms ids1, ids2;
	for ( Atoms::const_iterator it = atoms1_.begin(); it != atoms1_.end(); ++it ) {
		id::NamedAtomID atom(pose::atom_id_to_named_atom_id( *it, src ) );
		if ( smap ) {
			atom.rsd() = (*smap)[ it->rsd() ];
		}
		id::AtomID id1( core::pose::named_atom_id_to_atom_id(atom, dest ));
		if ( !id1.valid() ) return NULL;
		ids1.push_back( id1 );
	}
	for ( Atoms::const_iterator it = atoms2_.begin(); it != atoms2_.end(); ++it ) {
		id::NamedAtomID atom(atom_id_to_named_atom_id( *it, src ) );
		if ( smap ) {
			atom.rsd() = (*smap)[ it->rsd() ];
		}
		id::AtomID id2( core::pose::named_atom_id_to_atom_id( atom, dest ));
		if ( !id2.valid() ) return NULL;
		ids2.push_back( id2 );
	}
	return new AmbiguousNMRDistanceConstraint( ids1, ids2, func_, score_type() );
}

void AmbiguousNMRDistanceConstraint::show( std::ostream& out ) const {
	out << "AmbiguousNMRDistanceConstraint (";
	for ( Atoms::const_iterator it = atoms1_.begin(); it != atoms1_.end(); ++it ) {
		out << it->atomno() << "," << it->rsd() << "||";
	}
	out << " - ";
	for ( Atoms::const_iterator it = atoms2_.begin(); it != atoms2_.end(); ++it ) {
		out << it->atomno() << "," << it->rsd();
	}
	func_->show( out );
}

inline bool is_aromatic( pose::Pose const& pose, core::Size res ) {
	using namespace core::chemical;
	return pose.residue_type( res ).aa() == aa_phe
		|| pose.residue_type( res ).aa() == aa_tyr
		|| pose.residue_type( res ).aa() == aa_trp ;
}

void parse_NMR_name( std::string name, core::Size res, AmbiguousNMRDistanceConstraint::Atoms& atoms, core::pose::Pose const& pose ) {
	using core::chemical::AA;
	using namespace core::chemical;
	using core::id::NamedAtomID;
	using core::pose::named_atom_id_to_atom_id;
	AA const aa( pose.residue_type( res ).aa() );
	//tr.Debug << "[ERROR]: name is " << name << " res " << res << " and name " << pose.residue( res ).name3() << std::endl;
	//use named_atom_id_to_atom_id because it can throw an exception if atoms are missing instead of hard-exit...

	//put histidine protons onto the the respective Carbon atom... it is just unpredictable if the respective H is present...
	if ( ( name.substr(0,2) == "HD" ) && aa == aa_his ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "CD2", res ), pose ) );
		return;
	}
	if ( ( name.substr(0,2) == "HE" ) && aa == aa_his ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "CE1", res ), pose ) );
		return;
	}

	///now fix all the problems with atomnames:
	if ( name == "H" && res == 1 && pose.is_fullatom() ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1H", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2H", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3H", res ), pose ) );
	} else if ( (name.substr(0,2) == "QA" || name == "HA") && ( aa == aa_gly )) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HA", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HA", res ), pose ) );
	} else if ( ( name == "QB" || name == "HB" ) && (
																								 aa == aa_ala
																								 || aa == aa_leu
																								 || aa == aa_ser
																								 || aa == aa_asn
																								 || aa == aa_gln
																								 || aa == aa_pro
																								 || aa == aa_lys
																								 || aa == aa_cys
																								 || aa == aa_asp
																								 || aa == aa_glu
																								 || aa == aa_arg
																								 || aa == aa_tyr
																								 || aa == aa_phe
																								 || aa == aa_trp
																								 || aa == aa_his
																								 || aa == aa_met ) ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HB", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HB", res ), pose ) );
		if ( aa == aa_ala ) atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HB", res ), pose ) );
	} else if ( ( name == "QD1" ||  name == "HD1" ) && (
																											aa == aa_ile
																											|| aa == aa_leu )) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HD1", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HD1", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HD1", res ), pose ) );
	} else if ( ( name == "QD2" || name == "HD2" ) && (
																										 aa == aa_leu
																										 ||aa == aa_asn )) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HD2", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HD2", res ), pose ) );
		if ( aa == aa_leu ) atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HD2", res ), pose ) );
	} else if ( name == "QQD" && aa == aa_leu ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HD2", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HD2", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HD2", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HD1", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HD1", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HD1", res ), pose ) );
	} else if ( ( name == "QD" || name == "HD" ) && (
																									 aa == aa_pro
																									 || aa == aa_lys
																									 || aa == aa_arg
																									 || aa == aa_tyr
																									 || aa == aa_phe
																									 || aa == aa_his ) ) {
		if ( aa == aa_arg || aa == aa_lys || aa == aa_pro ) {
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HD", res ), pose ) );
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HD", res ), pose ) );
		} else { // tyr, phe, his (his doesn't get down here... )
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "HD1", res ), pose ) );
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "HD2", res ), pose ) );
		}
	} else if ( ( name == "QE" || name == "HE" ) && (
																									 aa == aa_tyr
																									 || aa == aa_phe
																									 || aa == aa_trp
																									 || aa == aa_met
																									 || aa == aa_lys ) ){
		if ( aa == aa_phe || aa == aa_tyr ) {
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "HE1", res ), pose ) );
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "HE2", res ), pose ) );
		}	else if ( aa == aa_trp ) {
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "HE1", res ), pose ) );
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "HE3", res ), pose ) );
		} else if ( aa == aa_met || aa == aa_lys ) { //MET LYS
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HE", res ), pose ) );
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HE", res ), pose ) );
			if ( aa == aa_met ) {
				atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HE", res ), pose ) );
			}
		} //QE MET LYS
	}  else if ( ( name == "QG" || name == "HG" ) && (
																									 aa == aa_met
																									 || aa == aa_gln
																									 || aa == aa_pro
																									 || aa == aa_lys
																									 || aa == aa_glu
																									 || aa == aa_arg )) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HG", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HG", res ), pose ) );
		//		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HE", res ), pose ) );
	} else if (( name == "QG2" || name == "HG2" ) && (
																									 aa == aa_ile
																									 || aa == aa_thr
																									 || aa == aa_val ) ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HG2", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HG2", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HG2", res ), pose ) );
	} else if ( ( name == "QQG" || name == "HG" ) && (
																										aa == aa_ile
																										|| aa == aa_val ) ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HG2", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HG2", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HG2", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HG1", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HG1", res ), pose ) );
		if ( aa != aa_ile ) atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HG1", res ), pose ) );
	} else if (( name == "QG1" || name == "HG1" ) && (
																									 aa == aa_ile
																									 || aa == aa_val ) ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HG1", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HG1", res ), pose ) );
		if ( aa != aa_ile ) atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HG1", res ), pose ) );
	} else if (( name == "QE2" || name == "HE2" || name =="HE" ) && aa == aa_gln ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HE2", res ), pose ) );
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HE2", res ), pose ) );
	} else if (( name == "HZ" || name == "QZ" ) && ( aa == aa_trp || aa == aa_lys ) ) {
		if ( aa == aa_lys ) {
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HZ", res ), pose ) );
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HZ", res ), pose ) );
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HZ", res ), pose ) );
		} else { //aa==trp
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "HZ2", res ), pose ) );
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "HZ3", res ), pose ) );
		}
	} else 	if ( name == "HB1" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HB", res ), pose ) );
	} else 	if ( name == "HB2" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HB", res ), pose ) );
	} else 	if ( name == "HB3" ) {
		if (  aa != aa_ala ) {
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HB", res ), pose ) ); //yeah they call it 2HB and 3HB...
		} else {
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HB", res ), pose ) );
		}
	}	else 	if ( name == "HD1" && !is_aromatic( pose, res ) ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HD", res ), pose ) );
	} else 	if ( name == "HD2" && !is_aromatic( pose, res ) ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HD", res ), pose ) );
	} else 	if ( name == "HD3" ) { //LYS, PRO, ARG  no other has HD3
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HD", res ), pose ) );

	} else 	if ( name == "HG1" && aa != aa_thr ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HG", res ), pose ) );
	} else 	if ( name == "HG2" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HG", res ), pose ) );
	} else 	if ( name == "HG3" ) { //GLU, ARG, GLN, MET
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HG", res ), pose ) );

	} else 	if ( name == "HA1" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HA", res ), pose ) );
	} else 	if ( name == "HA2" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HA", res ), pose ) );
	} else 	if ( name == "HA3" ) { //GLY
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HA", res ), pose ) );
	}

	else if (( name == "HE1" || name == "HE2" || name=="HE3" ) && !is_aromatic( pose, res ) )
		{
			if ( name == "HE3" && aa != aa_met ) {
				atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HE", res ), pose ) ); //e.g. LYS
			} else {
				atoms.push_back( id::AtomID( pose.residue_type(res).atom_index(name.substr(2,1)+name.substr(0,2)), res ) );
			}
		}

	else if ( name == "HD11" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HD1", res ), pose ) );
	} else if ( name == "HD12" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HD1", res ), pose ) );
	} else 	if ( name == "HD13" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HD1", res ), pose ) );
	} else 	if ( name == "HD21" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HD2", res ), pose ) );
	} else if ( name == "HD22" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HD2", res ), pose ) );
	} else 	if ( name == "HD23" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HD2", res ), pose ) );
	} else 	if ( name == "HG11" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HG1", res ), pose ) );
	} else 	if ( name == "HG12" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HG1", res ), pose ) );
	} else if ( name == "HG13" ) {
		if ( aa == aa_ile ) {
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HG1", res ), pose ) );
		} else {
			atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HG1", res ), pose ) );
		}
	} else 	if ( name == "HG21" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HG2", res ), pose ) );
	} else 	if ( name == "HG22" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HG2", res ), pose ) );
	} else if ( name == "HG23" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HG2", res ), pose ) );
	} else if ( name == "HE11" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HE1", res ), pose ) );
	} else if ( name == "HE12" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HE1", res ), pose ) );
	} else if ( name == "HE13" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "3HE1", res ), pose ) );
	} else if ( name == "HE21" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "1HE2", res ), pose ) );
	} else if ( name == "HE22" ) {
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( "2HE2", res ), pose ) );
	} else {
		tr.Trace << "adding " << id::NamedAtomID( name, res ) << " "
						 << pose.residue_type( res ).name3() << std::endl;
		atoms.push_back( named_atom_id_to_atom_id( NamedAtomID( name, res ), pose ) );
		tr.Trace << "as atom: " << atoms.back();
	}

	while ( atoms.size() && !atoms.back().valid() ) atoms.pop_back();
}

bool requires_CB_mapping( AmbiguousNMRDistanceConstraint::Atoms atoms, pose::Pose const& pose ) {
	basic::ProfileThis doit( basic::NOESY_ASSIGN_REQUIRES_CB_MAPPING );
	using core::chemical::AA;
	using namespace core::chemical;
	Size ct_backbone( 0 );
	for ( Size i=1; i<=atoms.size(); ++i ) {
		chemical::ResidueType const& res_type( pose.residue_type( atoms[ i ].rsd() ) );
		ct_backbone += res_type.atom_is_backbone( atoms[ i ].atomno() );
	}
	chemical::ResidueType const& res_type( pose.residue_type( atoms[ 1 ].rsd() ) );
	ct_backbone -= ( ct_backbone == 1 && res_type.aa() != aa_gly && res_type.atom_index( "CB" ) == atoms[ 1 ].atomno() );
	if ( ct_backbone == 0 ) return true;
	if ( ct_backbone == atoms.size() ) return false;
	runtime_assert( 0 ); //this should not happen, as there should only be methyl-Hs combined into one constraint
	return false;
}


	///

void combine_NMR_atom_string( AmbiguousNMRDistanceConstraint::Atoms atoms, std::string &atom_str, pose::Pose const& pose) {
	basic::ProfileThis doit( basic::NOESY_ASSIGN_NMR_STRING );
	atom_str = " BOGUS ";
	if ( atoms.size() == 1 ) {
		atom_str = pose.residue(  atoms[ 1 ].rsd() ).atom_name(  atoms[ 1 ].atomno() );
	} else if ( atoms.size() == 2 && atoms[ 1 ].rsd() == atoms[ 2 ].rsd()
			&& pose.residue(  atoms[ 1 ].rsd() ).atom_name(  atoms[ 1 ].atomno() ).substr( 1, 2 ) == "HB"
		&& pose.residue(  atoms[ 2 ].rsd() ).atom_name(  atoms[ 2 ].atomno() ).substr( 1, 2 ) == "HB" ) {
		atom_str = "QB";
	} else if ( atoms.size() == 3
		&& atoms[ 1 ].rsd() == atoms[ 2 ].rsd() && atoms[ 1 ].rsd() == atoms[ 3 ].rsd() ) {
		std::string methyl = pose.residue(  atoms[ 1 ].rsd() ).atom_name(  atoms[ 1 ].atomno() ).substr( 1, 3 );
		if ( pose.residue(  atoms[ 2 ].rsd() ).atom_name(  atoms[ 2 ].atomno() ).substr( 1, 3 ) == methyl
			&& pose.residue(  atoms[ 3 ].rsd() ).atom_name(  atoms[ 3 ].atomno() ).substr( 1, 3 ) == methyl ) {
			if ( atoms[ 1 ].rsd() == 1 ) { //1H, 2H, 3H Nterminus
				atom_str = "H";
			} else {
				atom_str = "Q"+methyl.substr(1);
			}
		} else {
			atom_str = " untranslatable-3-atom-combi";
			tr.Trace << "could not translate: "
							 << pose.residue(  atoms[ 1 ].rsd() ).atom_name(  atoms[ 1 ].atomno() ) << " "
							 << pose.residue(  atoms[ 2 ].rsd() ).atom_name(  atoms[ 2 ].atomno() ) << " "
							 << pose.residue(  atoms[ 3 ].rsd() ).atom_name(  atoms[ 3 ].atomno() ) << std::endl;
		}
	} else if ( atoms.size() == 2 //QE or QE2 or QD/PHE
		&& atoms[ 1 ].rsd() == atoms[ 2 ].rsd() ) {
		if ( pose.residue(  atoms[ 1 ].rsd() ).atom_name(  atoms[ 1 ].atomno() ).substr(0,2)==" H" ) {
			atom_str = "Q"+pose.residue(  atoms[ 1 ].rsd() ).atom_name(  atoms[ 1 ].atomno() ).substr(2,1);
		} else {
			std::string methyl = pose.residue(  atoms[ 1 ].rsd() ).atom_name(  atoms[ 1 ].atomno() ).substr( 1, 3 );
			if ( pose.residue(  atoms[ 2 ].rsd() ).atom_name(  atoms[ 2 ].atomno() ).substr( 1, 3 ) == methyl ) {
				atom_str = "Q"+methyl.substr(1);
			} else {
				atom_str = " untranslatable-2-methyl-atom-combi ";
				tr.Trace << "could not translate: " << pose.residue(  atoms[ 1 ].rsd() ).atom_name(  atoms[ 1 ].atomno() ) << " " << pose.residue(  atoms[ 2 ].rsd() ).atom_name(  atoms[ 2 ].atomno() ) << std::endl;
			}
		}
	} else if ( atoms.size() == 6
		&& atoms[ 1 ].rsd() == atoms[ 2 ].rsd() && atoms[ 1 ].rsd() == atoms[ 3 ].rsd()
		&& atoms[ 1 ].rsd() == atoms[ 4 ].rsd() && atoms[ 1 ].rsd() == atoms[ 5 ].rsd()
		&& atoms[ 1 ].rsd() == atoms[ 6 ].rsd() ) {
		std::string methyl = pose.residue(  atoms[ 1 ].rsd() ).atom_name(  atoms[ 1 ].atomno() ).substr( 1, 2 );
		if ( pose.residue(  atoms[ 2 ].rsd() ).atom_name(  atoms[ 2 ].atomno() ).substr( 1, 2 ) == methyl
			&& pose.residue(  atoms[ 3 ].rsd() ).atom_name(  atoms[ 3 ].atomno() ).substr( 1, 2 ) == methyl ) {
			atom_str = "QQ"+methyl.substr(1);
		} else {
			atom_str = " untranslatable-6-atom-combi ";
			tr.Trace << "could not translate: "
							 << pose.residue(  atoms[ 1 ].rsd() ).atom_name(  atoms[ 1 ].atomno() ) << " "
							 << pose.residue(  atoms[ 2 ].rsd() ).atom_name(  atoms[ 2 ].atomno() ) << " "
							 << pose.residue(  atoms[ 3 ].rsd() ).atom_name(  atoms[ 3 ].atomno() ) << " "
							 << pose.residue(  atoms[ 4 ].rsd() ).atom_name(  atoms[ 4 ].atomno() ) << " "
							 << pose.residue(  atoms[ 5 ].rsd() ).atom_name(  atoms[ 5 ].atomno() ) << " "
							 << pose.residue(  atoms[ 6 ].rsd() ).atom_name(  atoms[ 6 ].atomno() ) << std::endl;
		}
	}
	std::istringstream eat_white_space( atom_str );
	eat_white_space >> atom_str;
}

void AmbiguousNMRDistanceConstraint::show_def( std::ostream& out, pose::Pose const& pose ) const {
	std::string def_atoms1;
	std::string def_atoms2;
	combine_NMR_atom_string( atoms1_, def_atoms1, pose);
	combine_NMR_atom_string( atoms2_, def_atoms2, pose);
	out << type() << " "
			<< def_atoms1 << " " << atoms1_[ 1 ].rsd() << " "
			<< def_atoms2 << " " << atoms2_[ 1 ].rsd() << " ";
	if ( func_ ) func_->show_definition( out );
	else out << std::endl;
}

ConstraintOP AmbiguousNMRDistanceConstraint::map_to_CB( pose::Pose const& fa_pose, pose::Pose const& centroid, core::Size &mapped ) const {
	bool still_ambiguous( false );
	std::string atom1,atom2;
	basic::ProfileThis doit( basic::NOESY_ASSIGN_MAP2CB );
	using core::id::NamedAtomID;
	using core::pose::named_atom_id_to_atom_id;
	mapped = 2;
	if ( requires_CB_mapping( atoms1_, fa_pose) ) atom1 = "CB";
	else {
		--mapped;
		combine_NMR_atom_string( atoms1_, atom1, fa_pose );
		still_ambiguous |= ( atoms1_.size() > 1 );
	}
	if ( requires_CB_mapping( atoms2_, fa_pose) ) atom2 = "CB";
	else {
		--mapped;
		combine_NMR_atom_string( atoms2_, atom2, fa_pose );
		still_ambiguous |= ( atoms2_.size() > 1 );
	}
	tr.Debug << "map_to_CB:  " << atom1 << " " << resid( 1 ) << " --> " << atom2 << " " << resid( 2 ) << (still_ambiguous ? " ambiguous " : " straight ") << std::endl;
	{ //scope for profile
		basic::ProfileThis doit( basic::NOESY_ASSIGN_MAP2CB_NEW );
		if ( still_ambiguous ) {
			return new AmbiguousNMRDistanceConstraint( id::NamedAtomID( atom1, resid( 1 ) ), id::NamedAtomID( atom2, resid( 2 ) ), centroid, func_, score_type() );
		} else {
			return new AtomPairConstraint( named_atom_id_to_atom_id( NamedAtomID( atom1, resid( 1 ) ), centroid ),
				named_atom_id_to_atom_id( NamedAtomID( atom2, resid( 2 ) ), centroid ), func_, score_type() );
		}
	} //scope
	return NULL; // cannot be reached
}

Real
AmbiguousNMRDistanceConstraint::dist( pose::Pose const & pose ) const {

	//	tr.Trace << "compute distance for  ";
	//	if ( tr.Trace.visible() ) show_def( tr.Trace, pose );
	//	tr.Trace << std::endl;

	//conformation::Conformation const & conformation( pose.conformation() );
	return dist( ConformationXYZ( pose.conformation() ) );
	/*
	Real cum_dist( 0.0 );
	for ( Atoms::const_iterator it1 = atoms1_.begin(); it1 != atoms1_.end(); ++it1 ) {
		for (Atoms::const_iterator it2 = atoms2_.begin(); it2 != atoms2_.end(); ++it2 ) {
			Vector const & xyz1( conformation.xyz( *it1 ) ), xyz2( conformation.xyz( *it2 ) );
			Vector const f2( xyz1 - xyz2 );
			Real const dist( f2.length() );
			Real const inv_dist( 1.0/dist );
			Real const inv_dist2( inv_dist*inv_dist );
			Real const inv_dist6( inv_dist2 * inv_dist2 * inv_dist2 );
			cum_dist += inv_dist6;
			tr.Trace << *it1 << " " << *it2 << " " << dist << " " << pow(cum_dist, -1.0/6) << std::endl;
		}
	}
	tr.Trace << "finished distance" << std::endl;
	return pow(cum_dist, -1.0/6 );
	*/
}

Real
AmbiguousNMRDistanceConstraint::dist(
	XYZ_Func const & xyz
) const
{
	return pow( inv_dist6( xyz ), -1.0/6 );
}

Real
AmbiguousNMRDistanceConstraint::inv_dist6(
	XYZ_Func const & xyz
) const
{
	Real cum_dist( 0.0 );
	for ( Atoms::const_iterator it1 = atoms1_.begin(); it1 != atoms1_.end(); ++it1 ) {
		for (Atoms::const_iterator it2 = atoms2_.begin(); it2 != atoms2_.end(); ++it2 ) {
			Vector const & xyz1( xyz( *it1 ) ), xyz2( xyz( *it2 ) );
			Vector const f2( xyz1 - xyz2 );
			Real const dist( f2.length() );
			Real const inv_dist( 1.0/dist );
			Real const inv_dist2( inv_dist*inv_dist );
			Real const inv_dist6( inv_dist2 * inv_dist2 * inv_dist2 );
			cum_dist += inv_dist6;
			//			tr.Trace << *it1 << " " << *it2 << " " << dist << " " << pow(cum_dist, -1.0/6) << std::endl;
		}
	}
	//	tr.Trace << "finished distance" << std::endl;
	return cum_dist;
}

/*Real
AmbiguousNMRDistanceConstraint::dist( conformation::Conformation const & conformation ) const {
	return pow( inv_dist6( conformation ), -1.0/6 );
}*/

/*Real AmbiguousNMRDistanceConstraint::inv_dist6( XYZ_Func const & xyz ) const {
	Real cum_dist( 0.0 );
	for ( Atoms::const_iterator it1 = atoms1_.begin(); it1 != atoms1_.end(); ++it1 ) {
		for (Atoms::const_iterator it2 = atoms2_.begin(); it2 != atoms2_.end(); ++it2 ) {
			Vector const & xyz1( xyz( *it1 ) ), xyz2( xyz( *it2 ) );
			Vector const f2( xyz1 - xyz2 );
			Real const dist( f2.length() );
			Real const inv_dist( 1.0/dist );
			Real const inv_dist2( inv_dist*inv_dist );
			cum_dist += inv_dist2 * inv_dist2 * inv_dist2;
		}
	}
	return cum_dist;
}*/

void AmbiguousNMRDistanceConstraint::score( XYZ_Func const & xyz, EnergyMap const &, EnergyMap & emap ) const {
	Real eff_dist = pow( inv_dist6( xyz ), -1.0/6 );
	emap[ this->score_type() ] += func( eff_dist );
}

Size AmbiguousNMRDistanceConstraint::show_violations(
	std::ostream& out,
	pose::Pose const& pose,
	Size verbose_level,
	Real threshold
) const {

//  	if ( verbose_level > 80 ) {
//  		out << "AmbiguousNMRDistanceConstraint ("
//  			<< pose.residue_type(atom1_.rsd() ).atom_name( atom1_.atomno() ) << ":"
//  			<< atom1_.atomno() << "," << atom1_.rsd() << "-"
//  			<< pose.residue_type(atom2_.rsd() ).atom_name( atom2_.atomno() ) << ":"
//  			<< atom2_.atomno() << "," << atom2_.rsd() << ") ";
// 	}
// 	if ( verbose_level > 120 ) { //don't ask but I had a really weird bug to track down!
// 		conformation::Conformation const & conformation( pose.conformation() );
// 		Vector const & xyz1( conformation.xyz( atom1_ ) ), xyz2( conformation.xyz( atom2_ ) );
// 		out << "\ncoords1: " << xyz1[ 1 ] << " " << xyz1[ 2 ] << " " << xyz1[ 3 ] << " --- ";
// 		out << "coords1: " << xyz2[ 1 ] << " " << xyz2[ 2 ] << " " << xyz2[ 3 ] << "\n";
// 	}

	return func_->show_violations( out, dist( pose ), verbose_level, threshold );
}

// atom deriv
void
AmbiguousNMRDistanceConstraint::fill_f1_f2(
	AtomID const & atom,
	XYZ_Func const & xyz,
	Vector & F1,
	Vector & F2,
	EnergyMap const & weights
) const
{

	bool not_methyl_1 = std::find( atoms1_.begin() , atoms1_.end() , atom ) == atoms1_.end();
	bool not_methyl_2 = std::find( atoms2_.begin() , atoms2_.end() , atom ) == atoms2_.end();
	if ( not_methyl_1 && not_methyl_2 ) {
		//		tr.Trace << "atom " << atom << " not found " << std::endl;
		return;
	}

// 	if ( !not_methyl_1 && !not_methyl_2 ) {
// 		tr.Error << "atom " << atom << " has been found in both atom arrays of AmbiguousNMRDistanceConstraint " << std::endl;
// 	}

	Real eff_dist = dist( xyz );
	Real out_wderiv( weights[ this->score_type() ] * dfunc( eff_dist ));
	Real in_deriv = -1.0/6.0 * pow( eff_dist, 7.0 );
	//	tr.Trace << "deriv for atom " << atom << eff_dist << " " << out_wderiv << " " << in_deriv << std::endl;

	Atoms const& the_other_atoms( not_methyl_1 ? atoms1_ : atoms2_ );
	//	tr.Trace << "the_other_atoms: " << the_other_atoms.size() << " " << the_other_atoms.front() << std::endl;
	for (Atoms::const_iterator it = the_other_atoms.begin(); it != the_other_atoms.end(); ++it ) {
		AtomID other_atom = *it;
		//		tr.Trace << "contribution from " << other_atom << " to " << atom << std::endl;
		Real rdist(0.0);
		Vector f1(0.0), f2(0.0);
		numeric::deriv::distance_f1_f2_deriv( xyz( atom ), xyz( other_atom ), rdist, f1, f2 );
		Real wderiv = -6.0*pow(rdist,-7.0) * in_deriv * out_wderiv ;
		//		tr.Trace << "wderiv " << wderiv << std::endl;
		F1 += wderiv * f1;
		F2 += wderiv * f2;

	}

}

ConstraintOP
AmbiguousNMRDistanceConstraint::remap_resid( core::id::SequenceMapping const &smap ) const
{
	//	runtime_assert( 0 );
	Atoms ids1, ids2;
	for ( Atoms::const_iterator it = atoms1_.begin(); it != atoms1_.end(); ++it ) {
		id::AtomID atom( it->atomno(), smap[ it->rsd() ] );
		if ( !atom.valid() ) return NULL;
		ids1.push_back( atom );
	}
	for ( Atoms::const_iterator it = atoms2_.begin(); it != atoms2_.end(); ++it ) {
		id::AtomID atom( it->atomno(), smap[ it->rsd() ] );
		if ( !atom.valid() ) return NULL;
		ids2.push_back( atom );
	}
	return new AmbiguousNMRDistanceConstraint( ids1, ids2, func_, score_type() );
}



AmbiguousNMRDistanceConstraint::AmbiguousNMRDistanceConstraint(
    id::NamedAtomID const & a1, //digests names like "QG1"
		id::NamedAtomID const & a2,
		core::pose::Pose const& pose,
	 	FuncOP func,
		ScoreType scoretype
) : Constraint( scoretype ),
		func_( func )
{
	parse_NMR_name( a1.atom(), a1.rsd(), atoms1_, pose );
	parse_NMR_name( a2.atom(), a2.rsd(), atoms2_, pose );
	// } catch ( EXCN_AtomNotFound &excn ) {
// 		tr.Error << "AmbiguousNMRDistanceConstraint cannot parse input atoms: " << excn << std::endl;
// 		atoms1_.clear();
// 		atoms2_.clear();
// 	}

	if ( atoms1_.size() == 0 || atoms2_.size() == 0 ) {
		tr.Warning << "Error constructing from atoms: read in atom names("
							 << a1.atom() << "," << a2.atom() << "), " << std::endl;
	}
}

///@details one line definition "AmbiguousNMRDistances atom1 res1 atom2 res2 function_type function_definition"
void
AmbiguousNMRDistanceConstraint::read_def(
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
		tr.Warning 	<< "ignored constraint (residue number to high for pose: " << pose.total_residue() << " !)"
			<< name1 << " " << name2 << " " << res1 << " " << res2 << std::endl;
		data.setstate( std::ios_base::failbit );
		return;
	}

	parse_NMR_name( name1, res1, atoms1_, pose );
	parse_NMR_name( name2, res2, atoms2_, pose );

	if ( atoms1_.size() == 0 || atoms2_.size() == 0 ) {
		tr.Warning << "Error reading atoms: read in atom names("
							 << name1 << "," << name2 << "), " << std::endl;
		//			<< "and found AtomIDs (" << atom1_ << "," << atom2_ << ")" << std::endl;
			data.setstate( std::ios_base::failbit );
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
} // read_def

core::Size
AmbiguousNMRDistanceConstraint::effective_sequence_separation( core::kinematics::ShortestPathInFoldTree const& sp ) const {
	return sp.dist( resid( 1 ), resid( 2 ) );
}

} // constraints
} // scoring
} // core
