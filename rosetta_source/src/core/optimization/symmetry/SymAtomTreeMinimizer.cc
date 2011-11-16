// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   core/optimization/symmetry/SymAtomTreeMinimizer.cc
/// @brief  High-level atom tree minimizer class for symmetrical minimization
/// @author Ingemar Andre

// Unit headers
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>

// Symmetry headers
// AUTO-REMOVED #include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>


// Package headers
#include <core/optimization/types.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
// AUTO-REMOVED #include <core/optimization/MinimizerMap.hh>
#include <core/optimization/symmetry/SymAtomTreeMultifunc.hh>
#include <core/optimization/symmetry/SymMinimizerMap.hh>

#include <core/optimization/NumericalDerivCheckResult.fwd.hh>

// Project headers
#include <core/id/DOF_ID.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>

#include <basic/Tracer.hh>

#include <ObjexxFCL/format.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
// AUTO-REMOVED #include <core/chemical/AtomType.hh>
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
#include <core/chemical/VariantType.fwd.hh>
// AUTO-REMOVED #include <core/chemical/types.hh>
#include <core/chemical/orbitals/ICoorOrbitalData.hh>
#include <core/chemical/orbitals/OrbitalType.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <core/chemical/sdf/MolData.fwd.hh>
#include <core/chemical/sdf/MolData.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/PointGraph.fwd.hh>
#include <core/conformation/PointGraphData.fwd.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>
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
#include <core/graph/ArrayPool.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/graph/Graph.hh>
#include <core/graph/UpperEdgeGraph.fwd.hh>
#include <core/graph/unordered_object_pool.fwd.hpp>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID_Map.fwd.hh>
#include <core/id/DOF_ID_Map.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/JumpID.fwd.hh>
#include <core/id/JumpID.hh>
#include <core/id/NamedAtomID.fwd.hh>
// AUTO-REMOVED #include <core/id/NamedAtomID.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/TorsionID.hh>
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
#include <core/kinematics/MinimizerMapBase.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/ResidueCoordinateChangeList.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/types.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/optimization/AtomTreeMinimizer.fwd.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/DOF_Node.fwd.hh>
// AUTO-REMOVED #include <core/optimization/DOF_Node.hh>
#include <core/optimization/LineMinimizer.fwd.hh>
// AUTO-REMOVED #include <core/optimization/LineMinimizer.hh>
// AUTO-REMOVED #include <core/optimization/MinimizerMap.fwd.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/optimization/Multifunc.fwd.hh>
#include <core/optimization/Multifunc.hh>
#include <core/optimization/symmetry/SymMinimizerMap.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/ContextGraph.fwd.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/scoring/DerivVectorPair.fwd.hh>
// AUTO-REMOVED #include <core/scoring/DerivVectorPair.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/EnergyGraph.fwd.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/LREnergyContainer.fwd.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/scoring/NeighborList.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/TenANeighborGraph.fwd.hh>
#include <core/scoring/TwelveANeighborGraph.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/Methods.hh>
#include <core/scoring/methods/TwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.fwd.hh>
// AUTO-REMOVED #include <utility/Bound.fwd.hh>
// AUTO-REMOVED #include <utility/Bound.hh>
#include <utility/PyAssert.hh>
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
// AUTO-REMOVED #include <utility/file/FileName.fwd.hh>
// AUTO-REMOVED #include <utility/file/FileName.hh>
// AUTO-REMOVED #include <utility/file/PathName.fwd.hh>
// AUTO-REMOVED #include <utility/file/PathName.hh>
// AUTO-REMOVED #include <utility/keys/AutoKey.fwd.hh>
// AUTO-REMOVED #include <utility/keys/AutoKey.hh>
// AUTO-REMOVED #include <utility/keys/Key.fwd.hh>
// AUTO-REMOVED #include <utility/keys/Key.hh>
#include <utility/keys/Key2Tuple.fwd.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/keys/Key3Tuple.fwd.hh>
#include <utility/keys/Key3Tuple.hh>
#include <utility/keys/Key4Tuple.fwd.hh>
#include <utility/keys/Key4Tuple.hh>
// AUTO-REMOVED #include <utility/keys/KeyLess.fwd.hh>
// AUTO-REMOVED #include <utility/keys/KeyLookup.fwd.hh>
// AUTO-REMOVED #include <utility/keys/KeyLookup.hh>
// AUTO-REMOVED #include <utility/keys/NoClient.fwd.hh>
// AUTO-REMOVED #include <utility/keys/NoClient.hh>
// AUTO-REMOVED #include <utility/keys/SmallKeyVector.fwd.hh>
// AUTO-REMOVED #include <utility/keys/SmallKeyVector.hh>
// AUTO-REMOVED #include <utility/keys/UserKey.fwd.hh>
// AUTO-REMOVED #include <utility/keys/VariantKey.fwd.hh>
// AUTO-REMOVED #include <utility/keys/VariantKey.hh>
// AUTO-REMOVED #include <utility/options/AnyOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/AnyOption.hh>
// AUTO-REMOVED #include <utility/options/AnyVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/AnyVectorOption.hh>
// AUTO-REMOVED #include <utility/options/BooleanOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/BooleanOption.hh>
// AUTO-REMOVED #include <utility/options/BooleanVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/BooleanVectorOption.hh>
// AUTO-REMOVED #include <utility/options/FileOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/FileOption.hh>
// AUTO-REMOVED #include <utility/options/FileVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/FileVectorOption.hh>
// AUTO-REMOVED #include <utility/options/IntegerOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/IntegerOption.hh>
// AUTO-REMOVED #include <utility/options/IntegerVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/IntegerVectorOption.hh>
// AUTO-REMOVED #include <utility/options/Option.fwd.hh>
// AUTO-REMOVED #include <utility/options/Option.hh>
// AUTO-REMOVED #include <utility/options/OptionCollection.fwd.hh>
// AUTO-REMOVED #include <utility/options/OptionCollection.hh>
// AUTO-REMOVED #include <utility/options/PathOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/PathOption.hh>
// AUTO-REMOVED #include <utility/options/PathVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/PathVectorOption.hh>
// AUTO-REMOVED #include <utility/options/RealOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/RealOption.hh>
// AUTO-REMOVED #include <utility/options/RealVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/RealVectorOption.hh>
// AUTO-REMOVED #include <utility/options/ScalarOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/ScalarOption.hh>
// AUTO-REMOVED #include <utility/options/ScalarOption_T_.fwd.hh>
// AUTO-REMOVED #include <utility/options/ScalarOption_T_.hh>
// AUTO-REMOVED #include <utility/options/StringOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/StringOption.hh>
// AUTO-REMOVED #include <utility/options/StringVectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/StringVectorOption.hh>
// AUTO-REMOVED #include <utility/options/VariantOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/VariantOption.hh>
// AUTO-REMOVED #include <utility/options/VectorOption.fwd.hh>
// AUTO-REMOVED #include <utility/options/VectorOption.hh>
// AUTO-REMOVED #include <utility/options/VectorOption_T_.fwd.hh>
// AUTO-REMOVED #include <utility/options/VectorOption_T_.hh>
// AUTO-REMOVED #include <utility/options/mpi_stderr.hh>
// AUTO-REMOVED #include <utility/options/keys/AnyOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/AnyOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/AnyVectorOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/AnyVectorOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/BooleanOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/BooleanOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/BooleanVectorOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/BooleanVectorOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/FileOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/FileOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/FileVectorOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/FileVectorOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/IntegerOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/IntegerOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/IntegerVectorOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/IntegerVectorOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/OptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/OptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/OptionKeys.hh>
// AUTO-REMOVED #include <utility/options/keys/PathOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/PathOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/PathVectorOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/PathVectorOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/RealOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/RealOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/RealVectorOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/RealVectorOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/ScalarOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/ScalarOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/StringOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/StringOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/StringVectorOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/StringVectorOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/VectorOptionKey.fwd.hh>
// AUTO-REMOVED #include <utility/options/keys/VectorOptionKey.hh>
// AUTO-REMOVED #include <utility/options/keys/all.hh>
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
#include <ObjexxFCL/Fstring.fwd.hh>
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
#include <ObjexxFCL/byte.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/proxy_const_assert.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/ubyte.fwd.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <istream>
#include <limits>
#include <list>
#include <map>
#include <ostream>
// AUTO-REMOVED #include <set>
#include <sstream>
#include <string>
#include <typeinfo>
// AUTO-REMOVED #include <utility>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableData.fwd.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/DataCache.fwd.hh>
#include <basic/datacache/DataCache.hh>
// AUTO-REMOVED #include <basic/options/keys/OptionKeys.hh>
// AUTO-REMOVED #include <basic/options/option.hh>
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
// AUTO-REMOVED #include <boost/config.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>
// AUTO-REMOVED #include <boost/pool/detail/mutex.hpp>
#include <boost/pool/poolfwd.hpp>

//Auto Headers



using namespace ObjexxFCL::fmt;

namespace core {
namespace optimization {
namespace symmetry {

///////////////////////////////////////////////////////////////////////////////
Real
SymAtomTreeMinimizer::run(
	pose::Pose & pose,
	kinematics::MoveMap const & move_map,
	scoring::ScoreFunction const & scorefxn,
	MinimizerOptions const & options
)
{
	using namespace core::conformation::symmetry;

	typedef SymmetryInfo::DOF_IDs DOF_IDs;

	bool const use_nblist( options.use_nblist() );

	// it's important that the structure be scored prior to nblist setup
	Real const start_score( scorefxn( pose ) );

	kinematics::MoveMap semisym_move_map;
	make_semisymmetric_movemap( pose, move_map, semisym_move_map );

	SymmetricConformation const & symm_conf ( dynamic_cast<SymmetricConformation const &> ( pose.conformation()) );
	assert( conformation::symmetry::is_symmetric( symm_conf ) );
	SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

	// setup the minimizer map using the semi-symetric min map
	SymMinimizerMap sym_min_map( pose, semisym_move_map, symm_info );
	//kinematics::DomainMap const & dm( sym_min_map.domain_map() );
	//for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
	//	std::cout << "(" << ii << ", " << dm(ii) << "), ";
	//	if ( ii % 10 == 0 ) std::cout << std::endl;
	//}
	//std::cout << std::endl;

	// if we are using the nblist, set it up
	if ( use_nblist ) {
		// setup a mask of the moving dofs
		pose.energies().set_use_nblist( pose, sym_min_map.domain_map(), options.nblist_auto_update() );
	}

	// etable functions must be initialized with asymm domain map
	scorefxn.setup_for_minimizing( pose, sym_min_map );

	// setup the function that we will pass to the low-level minimizer
	SymAtomTreeMultifunc f( pose, sym_min_map, scorefxn,
		options.deriv_check(), options.deriv_check_verbose() );

	// starting position -- "dofs" = Degrees Of Freedom
	Multivec dofs( sym_min_map.nangles() );
	sym_min_map.copy_dofs_from_pose( pose, dofs );

	Real const start_func( f( dofs ) );

	// std::cout << "Start score: " << start_score << " start_func " << start_func << std::endl;
	// APL Note: there might be a discrepancy between start_score and start_func if bb/bb hydrogen bonds
	// are being used.  Why?  Well, during minimization, bb/bb hbonds are calculated in the residue_pair_energy_ext
	// calls, but in regular scoring, bb/bb hbonds are calculated in setup_for_scoring.  This means that
	// bb/bb hbond energies are not stored in the EnergyGraph. Hbonding residue pairs that are not moving wrt
	// each other are not rescored during minimization, so their energies would be stored in the "fixed_energies_"
	// EnergyMap in the MinimizationGraph -- except for the fact that this EnergyMap is filled using energies
	// stored in the EnergyGraph, and the bb/bb hbond energies are not stored there.
	// The delta between start_score and start_func is not worrisome, however, because these energies are constant.
	// The delta has no impact on the minimizer's behavior.

	// now do the optimization with the low-level minimizer function
	Minimizer minimizer( f, options );
	minimizer.run( dofs );

	Real const end_func( f( dofs ) );

	//std::cout << "end_func:" << std::endl;
	//pose.energies().show( std::cout );

	// turn off nblist
	if ( use_nblist ) pose.energies().reset_nblist();

	// if we were doing rigid-body minimization, fold the rotation and
	// translation offsets into the jump transforms
	//
	// also sets rb dofs to 0.0, so in principle func value should be the same
	//
	sym_min_map.reset_jump_rb_deltas( pose, dofs );

	// rescore
	Real const end_score( scorefxn( pose ) );

	//std::cout << "end_score:" << std::endl;
	//pose.energies().show( std::cout );

	// we may not really need all these extra function evaluations
	// good for diagnostics though

	static basic::Tracer core_optimize( "core.optimize",  basic::t_debug);
	core_optimize << "SymAtomTreeMinimizer::run: nangles= " << sym_min_map.nangles() <<
		" start_score: " << F(12,3,start_score) <<
		" start_func: "  << F(12,3,start_func ) <<
		" end_score: "   << F(12,3,end_score  ) <<
		" end_func: "    << F(12,3,end_func   ) << std::endl;


	return end_score;
}

void
SymAtomTreeMinimizer::make_semisymmetric_movemap(
	pose::Pose & pose,
	kinematics::MoveMap const & move_map_sym,
	kinematics::MoveMap & move_map_semisym
)
{
	using namespace core::conformation::symmetry;
	typedef SymmetryInfo::DOF_IDs DOF_IDs;

	SymmetricConformation const & SymmConf (
		dynamic_cast<SymmetricConformation const &> ( pose.conformation()) );
	assert( conformation::symmetry::is_symmetric( SymmConf ) );
	SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

	// copy bb/chi DOFs from symm movemap
	move_map_semisym = move_map_sym;

	for ( Size i=1; i<= pose.conformation().fold_tree().num_jump(); ++i ) {
		id::DOF_ID const & null_id
		( pose.conformation().dof_id_from_torsion_id(id::TorsionID(i,id::JUMP,1)));
		if ( symm_info->get_dof_derivative_weight( null_id , SymmConf ) > 0 ) {
			// if this is not the master get the master
			if ( symm_info->jump_is_independent( i ) ) continue;

			Size master_i = symm_info->jump_follows( i );

			for ( int j=1; j<= 6; ++j ) {
			id::DOF_ID const & id_master
				( pose.conformation().dof_id_from_torsion_id(id::TorsionID(master_i,id::JUMP,j)));
				id::DOF_ID const & id
				( pose.conformation().dof_id_from_torsion_id(id::TorsionID(i,id::JUMP,j)));

				bool allow ( move_map_sym.get( id_master ) );
				move_map_semisym.set(id, allow );
			}
		}
	}

}


void
SymAtomTreeMinimizer::make_assymetric_movemap(
	pose::Pose & pose,
	kinematics::MoveMap const & move_map_sym,
	kinematics::MoveMap & move_map_asym
)
{
	using namespace core::conformation::symmetry;
	typedef SymmetryInfo::DOF_IDs DOF_IDs;

	SymmetricConformation const & SymmConf (
		dynamic_cast<SymmetricConformation const &> ( pose.conformation()) );
	assert( conformation::symmetry::is_symmetric( SymmConf ) );
	SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

	for ( Size i=1; i<= pose.conformation().size(); ++i ) {
		if ( symm_info->bb_is_independent( i ) ) {
			bool bb ( move_map_sym.get_bb(i) );
			bool chi ( move_map_sym.get_chi(i) );
			move_map_asym.set_bb ( i, bb );
			move_map_asym.set_chi( i, chi );
			for ( std::vector< Size>::const_iterator
					clone     = symm_info->bb_clones( i ).begin(),
					clone_end = symm_info->bb_clones( i ).end();
					clone != clone_end; ++clone ){
				move_map_asym.set_bb ( *clone, bb );
				move_map_asym.set_chi( *clone, chi );
			}
		}
	}
	for ( Size i=1; i<= pose.conformation().fold_tree().num_jump(); ++i ) {
		if ( symm_info->jump_is_independent( i ) ) {
			for ( int j=1; j<= 6; ++j ) {
				id::DOF_ID const & id
				( pose.conformation().dof_id_from_torsion_id(id::TorsionID(i,id::JUMP,j)));
				DOF_IDs const & dofs( symm_info->dependent_dofs( id, SymmConf ) );
				bool allow ( move_map_sym.get( id ) );
				move_map_asym.set(id, allow );
				for ( DOF_IDs::const_iterator dof =dofs.begin(), dofe= dofs.end(); dof != dofe; ++dof ) {
					move_map_asym.set( *dof, allow );
				}
			}
		}
	}
}

} // symmetry
} // namespace optimization
} // namespace core
