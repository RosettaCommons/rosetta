// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Oliver Lange
/// @author Mike Tyka

// Unit Headers
#include <protocols/loops/util.hh>

// Package Headers
#include <protocols/loops/Loops.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/kinematics/MoveMap.hh>
#include <protocols/loops/loops_main.hh> //for getting ss from dssp

#include <core/fragment/SecondaryStructure.hh>
// AUTO-REMOVED #include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/util.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>

#include <core/conformation/util.hh> //idealize
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameList.hh>
#include <core/fragment/FragSet.hh>
#ifdef WIN32
#include <core/fragment/FragID.hh>
#endif
#include <core/pose/util.hh>
#include <core/id/AtomID.hh>
#include <core/scoring/rms_util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// Utility headers
#include <basic/Tracer.hh>

//numeric headers

//// C++ headers
#include <list>

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
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/MMAtomType.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResConnID.fwd.hh>
#include <core/chemical/ResConnID.hh>
#include <core/chemical/ResidueConnection.fwd.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueSelector.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
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
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
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
#include <core/fragment/BaseCacheUnit.hh>
#include <core/fragment/FragData.fwd.hh>
// AUTO-REMOVED #include <core/fragment/FragData.hh>
#include <core/fragment/FragID.fwd.hh>
// AUTO-REMOVED #include <core/fragment/FragID_Iterator.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/Frame.fwd.hh>
#include <core/fragment/FrameIterator.fwd.hh>
// AUTO-REMOVED #include <core/fragment/FrameIterator.hh>
// AUTO-REMOVED #include <core/fragment/FrameIteratorWorker_.fwd.hh>
// AUTO-REMOVED #include <core/fragment/FrameIteratorWorker_.hh>
#include <core/fragment/FrameList.fwd.hh>
#include <core/fragment/SecondaryStructure.fwd.hh>
#include <core/fragment/SingleResidueFragData.fwd.hh>
// AUTO-REMOVED #include <core/fragment/SingleResidueFragData.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Map.fwd.hh>
// AUTO-REMOVED #include <core/id/DOF_ID_Map.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/JumpID.fwd.hh>
#include <core/id/JumpID.hh>
#include <core/id/NamedAtomID.fwd.hh>
// AUTO-REMOVED #include <core/id/NamedAtomID.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
// AUTO-REMOVED #include <core/id/SequenceMapping.hh>
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
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/ResidueCoordinateChangeList.fwd.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/types.hh>
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
#include <core/scoring/EnergyGraph.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/LREnergyContainer.fwd.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/Func.fwd.hh>
#include <core/scoring/constraints/Func.hh>
#include <core/scoring/constraints/FuncFactory.fwd.hh>
#include <core/scoring/constraints/HarmonicFunc.fwd.hh>
#include <core/scoring/constraints/XYZ_Func.fwd.hh>
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/TwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.fwd.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.fwd.hh>
#include <protocols/loops/Loop.fwd.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/PyAssert.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
// AUTO-REMOVED #include <utility/stream_util.hh>
// AUTO-REMOVED #include <utility/string_util.hh>
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
// AUTO-REMOVED #include <numeric/NumericTraits.hh>
// AUTO-REMOVED #include <numeric/constants.hh>
// AUTO-REMOVED #include <numeric/conversions.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
// AUTO-REMOVED #include <numeric/sphericalVector.hh>
#include <numeric/trig.functions.hh>
// AUTO-REMOVED #include <numeric/types.hh>
#include <numeric/xyz.functions.fwd.hh>
// AUTO-REMOVED #include <numeric/xyz.functions.hh>
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
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/proxy_const_assert.hh>
#include <ObjexxFCL/ubyte.fwd.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <istream>
#include <limits>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
// AUTO-REMOVED #include <typeinfo>
#include <utility>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
// AUTO-REMOVED #include <basic/datacache/BasicDataCache.hh>
// AUTO-REMOVED #include <basic/datacache/CacheableData.fwd.hh>
// AUTO-REMOVED #include <basic/datacache/CacheableData.hh>
// AUTO-REMOVED #include <basic/datacache/DataCache.fwd.hh>
// AUTO-REMOVED #include <basic/datacache/DataCache.hh>
#include <basic/options/keys/OptionKeys.hh>
// AUTO-REMOVED #include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>
#include <boost/unordered/unordered_map.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <protocols/rosetta_scripts/util.hh>
#include <utility/string_util.hh>

//Auto Headers

// Auto-header: duplicate removed #include <core/pose/util.hh>

//Auto using namespaces
namespace ObjexxFCL { namespace fmt { } } using namespace ObjexxFCL::fmt; // AUTO USING NS

static basic::Tracer TR("protocols.loops.util");

namespace protocols {
namespace loops {

/// TODO(someone) so bad!
using namespace core;
using namespace pose;
using namespace kinematics;

// Numeric constants
static const double EXT_PHI = -150;
static const double EXT_PSI = +150;
static const double EXT_OMG =  180;

void
fix_with_coord_cst( Loops const& rigid, core::pose::Pose& pose, bool bCstAllAtom, utility::vector1< core::Real > &weights ) {
	bool const bReadWeights = ( weights.size() >= pose.total_residue() );
	if ( !bReadWeights ) {
		weights.resize( pose.total_residue() );
	}

  for ( Loops::const_iterator it = rigid.begin(), eit = rigid.end();
	it!=eit; ++it ) {
    for ( Size pos = it->start(); pos <= it->stop(); ++pos ) {
      Size const seq_dist( std::min( (int) pos - it->start(), (int) it->stop() - pos ) + 1);
			Real coord_sdev;
			if ( bReadWeights ) {
				coord_sdev = weights[ pos ];
			} else {
				coord_sdev = ( 1.0/seq_dist ); //or something else?
				weights[ pos ] = coord_sdev;
			}
      conformation::Residue const & rsd( pose.residue( pos ) );
			if ( bCstAllAtom ) {
				for ( Size ii = 1; ii<= rsd.natoms(); ++ii ) {
					pose.add_constraint( new scoring::constraints::CoordinateConstraint(
	      id::AtomID( ii, pos),
	      id::AtomID( 1, pos ) /*this is completely ignored! */,
	      rsd.xyz( ii ),
	      new scoring::constraints::HarmonicFunc( 0.0, coord_sdev )
						) );
				}
      } else {
				id::AtomID atomID( pose.residue_type(pos).atom_index("CA"), pos );
				pose.add_constraint( new scoring::constraints::CoordinateConstraint(
						atomID,
						id::AtomID( 1, pos ) /*this is completely ignored! */,
						rsd.xyz( atomID.atomno() ),
						new scoring::constraints::HarmonicFunc( 0.0, coord_sdev )
					) );
      }
    }
  }
}

///@brief get frags that are fully within the Loop --- shorten(=true/false) frags that are close to the end of loops.
void select_loop_frags(
				 loops::Loops const& loops,
				 core::fragment::FragSet& source,
				 core::fragment::FragSet& loop_frags,
				 Size shorten
) {
	using namespace core::fragment;
	//assuming backbone degrees of freedom are wanted by the loops.
	// Jumps are filtered out
	kinematics::MoveMap movemap;
	movemap.set_bb( false );
	loops.switch_movemap( movemap, id::BB, true );


	InsertMap insert_map;
	InsertSize insert_size;
	source.generate_insert_map( movemap, insert_map, insert_size);

	{ //debug
		Size const total_insert = insert_map.size();
		TR.Trace << "size of insertmap: " << total_insert << " -- ";
		for ( Size i = 1; i<=total_insert; i++ ) TR.Trace << " " << insert_map[ i ];
		TR.Trace << "insert_size: \nRESIDUES: ";
		for ( Size i = 1; i<=insert_map[ total_insert ]; i++ ) TR.Trace << " " << RJ(3,i);
		TR.Trace <<"\nINSERT:   ";
		for ( Size i = 1; i<=insert_map[ total_insert ]; i++ ) TR.Trace << " " << RJ(3,insert_size[ i ]);
		TR.Trace << std::endl;
	}

	Size const total_insert = insert_map.size();

	for ( Size i = 1; i<=total_insert; i++ ) {
		Size const pos ( insert_map[ i ] );
		Size const size ( insert_size[ pos ] );
		FrameList copy_frames;
		source.frames( pos, copy_frames );
		for ( FrameList::iterator it = copy_frames.begin(), eit = copy_frames.end();
					it!=eit; ++it ) {
			TR.Trace << "add frame at pos " << pos << " " << (*it)->length() << " insert_size " << size << std::endl;
			if ( (*it)->length() == size ) loop_frags.add( *it );
			else if ( shorten && size > shorten ) loop_frags.add( (*it)->generate_sub_frame( size ) );
		}
	}
} //select_loop_frags

void safe_set_extended_torsions_and_idealize_loops(const protocols::loops::Loops& loops,
                                                   core::pose::Pose* pose) {
  if (loops.empty())
    return;

  set_extended_torsions_and_idealize_loops(*pose, loops);
}

void set_extended_torsions_and_idealize_loops(core::pose::Pose& pose, loops::Loops loops) {
  using core::Size;
  using core::conformation::Conformation;
  using protocols::loops::Loop;
  using protocols::loops::Loops;

	// if no loops, we want to have extended structure everywhere.
	// it is a by-value parameter -- as intenden the change is kept local
	if (loops.empty())
		loops.add_loop(1, pose.total_residue(), 0);

	TR.Debug << "extend structure for " << loops << std::endl;

  Conformation& conf = pose.conformation();
  for (Loops::const_iterator i = loops.begin(); i != loops.end(); ++i) {
    const Loop& loop = *i;

    for (Size j = loop.start(); j <= loop.stop(); ++j) {
			core::conformation::idealize_position(j, conf);
      pose.set_phi(j, EXT_PHI);
      pose.set_psi(j, EXT_PSI);
      pose.set_omega(j, EXT_OMG);
    }
  }
}

void addScoresForLoopParts(
	core::pose::Pose & pose,
	loops::Loops loops,
	const core::scoring::ScoreFunction & scorefxn,
	core::pose::Pose & native_pose,
	core::Size nloops
) {

	using namespace core;

	Size nres = pose.total_residue();
	utility::vector1< core::Size > all_loop_list;
	for( Size i = 1; i < nres; i ++ ){
		if( loops.is_loop_residue(i) ) all_loop_list.push_back( i );
	}
	scorefxn(pose);
	setPoseExtraScores( pose, "ScoreCore", 	scorefxn.get_sub_score_exclude_res( pose, all_loop_list ) );

	Real score =  scorefxn(pose);

	for( Size l = 1; l <= nloops; l++ ){
		if( l > loops.size() ){
			setPoseExtraScores( pose, "ScoreLoopI" + ObjexxFCL::right_string_of(l,3,'0'), 0 );
			setPoseExtraScores( pose, "ScoreLoopC" + ObjexxFCL::right_string_of(l,3,'0'), 0 );
			setPoseExtraScores( pose, "ScoreLoopL" + ObjexxFCL::right_string_of(l,3,'0'), 0 );
			continue;
		}
		utility::vector1< core::Size > loop_list;
		utility::vector1< core::Size > non_loop_list;
		for( Size i = 1; i < nres; i ++ ){
			if( ( i < loops[l].start() ) || ( i > loops[l].stop() ) ){
				loop_list.push_back( i );
			}else{
				non_loop_list.push_back( i );
			}
		}

		Real loopscore = scorefxn.get_sub_score_exclude_res( pose, loop_list );
		Real nonloopscore = scorefxn.get_sub_score_exclude_res( pose, non_loop_list );
		setPoseExtraScores( pose, "ScoreLoopI" + ObjexxFCL::right_string_of(l,3,'0'), loopscore );
		setPoseExtraScores( pose, "ScoreLoopC" + ObjexxFCL::right_string_of(l,3,'0'), score - loopscore - nonloopscore );
		setPoseExtraScores( pose, "ScoreLoopL" + ObjexxFCL::right_string_of(l,3,'0'), score - nonloopscore );
	}

	// Work out RMS values too

	core::pose::Pose native_pose_super = native_pose;
	id::AtomID_Map< id::AtomID > atom_map;
	core::pose::initialize_atomid_map( atom_map, native_pose_super, core::id::BOGUS_ATOM_ID );
	for ( core::Size ir=1; ir <= native_pose.total_residue(); ++ir ) {
		runtime_assert( ir <=  pose.total_residue() );
		runtime_assert( ir <=  native_pose_super.total_residue() );
		if( ( !loops.is_loop_residue( ir ) ) && pose.residue(ir).is_protein() ) {
			id::AtomID const id1( native_pose_super.residue(ir).atom_index("CA"), ir );
			id::AtomID const id2( pose.residue(ir).atom_index("CA"), ir );
			atom_map.set(id1, id2);
		}
	}
	core::scoring::superimpose_pose( native_pose_super, pose, atom_map );

	int corelength;
	setPoseExtraScores(	pose, "corerms", native_loop_core_CA_rmsd( native_pose, pose, loops, corelength )	);

	for( Size l = 1; l <= nloops; l++ ){
		if( l > loops.size() ){
			setPoseExtraScores( pose, "RMSLoop" + ObjexxFCL::right_string_of(l,3,'0'), 0 );
			continue;
		}
		Loops temploops;
		temploops.add_loop( loops[l] );
		setPoseExtraScores( pose, "RMSLoop" + ObjexxFCL::right_string_of(l,3,'0'),
												loops::loop_rmsd( native_pose_super, pose, temploops, true ) );
	}
}

loops::Loops compute_ss_regions(
	core::Real max_loop_frac,
	core::Size min_length,
	core::fragment::SecondaryStructure const & ss
) {

	using core::Size;

	Size start( 0 );
	Size last( 0 );
	Size max_gap( 2 );
	loops::Loops ss_regions;
	for ( Size pos = 1; pos <= ss.total_residue(); ++pos ) {
		if ( ss.loop_fraction( pos ) <= max_loop_frac ) {
			if ( !start ) {
				start = pos;
				last = pos - 1;
			}
			if ( last + max_gap < pos ) {
				if ( last - start >= min_length ) {
					ss_regions.add_loop( start, last );
				}
				start=0;
			}
			last = pos;
		}
	}
	return ss_regions;
}

core::scoring::ScoreFunctionOP get_cen_scorefxn() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	std::string const weights( option[ OptionKeys::loops::cen_weights ]() ),
		patch( option[ OptionKeys::loops::cen_patch ]() );
	return scoring::ScoreFunctionFactory::create_score_function(
		weights, patch
	);
}

core::scoring::ScoreFunctionOP get_fa_scorefxn() {
	return core::scoring::getScoreFunction();
}


void add_coordinate_constraints_to_pose( core::pose::Pose & pose, const core::pose::Pose &constraint_target_pose,  protocols::loops::Loops &exclude_regions ){
  using namespace conformation;
  using namespace core;
  using namespace core::scoring;
  using namespace core::scoring::constraints;
  using namespace id;
  using namespace scoring::constraints;

  core::Size nnonvrt_cst_target = constraint_target_pose.total_residue();
  core::Size nnonvrt_pose = pose.total_residue();

  while ( pose.residue( nnonvrt_pose ).aa() == core::chemical::aa_vrt ) { nnonvrt_pose--; }
  while ( constraint_target_pose.residue( nnonvrt_cst_target ).aa() == core::chemical::aa_vrt ) { nnonvrt_cst_target--; }

  protocols::loops::Loops coordconstraint_segments;
  coordconstraint_segments = exclude_regions.invert( nnonvrt_cst_target );

  //TR.Info << coordconstraint_segments << std::endl;

  if ( nnonvrt_pose != nnonvrt_cst_target ) {
    TR.Error << "ERROR coord constraint pose length mismatch with input pose: " << nnonvrt_cst_target << " vs. " << nnonvrt_pose << std::endl;
    utility_exit();
  }

  if ( pose.residue( pose.fold_tree().root() ).aa() != core::chemical::aa_vrt ) {
    pose.append_residue_by_jump
      ( *ResidueFactory::create_residue( pose.residue(1).residue_type_set().name_map( "VRT" ) ),
        pose.total_residue()/2 );
  }


  Size nres = pose.total_residue();
  Real const coord_sdev( 0.5 );
  for ( Size i = 1; i<= (Size)nres; ++i ) {
    if ( i==(Size)pose.fold_tree().root() ) continue;
    if( coordconstraint_segments.is_loop_residue( i ) ) {
      Residue const & nat_i_rsd( pose.residue(i) );
      for ( Size ii = 1; ii<= nat_i_rsd.last_backbone_atom(); ++ii ) {
        pose.add_constraint( new CoordinateConstraint(
          AtomID(ii,i), AtomID(1,nres), nat_i_rsd.xyz( ii ),
          new HarmonicFunc( 0.0, coord_sdev ) ) );
      }
    }
  }

}

Loops
loops_from_string( std::string const loop_str, core::pose::Pose const & pose ){
  utility::vector1< std::string > const loops_vec( utility::string_split( loop_str, ',' ) );

// each loop should have the format loop_start:loop_end:cut
// if cut is not set then it's taken to be 0. Residue numbering can follow the
// pdb numbering
  Loops loops_from_tag;
 	foreach( std::string const residue_pair, loops_vec ){
    utility::vector1< std::string > const residues( utility::string_split( residue_pair, ':' ) );
    runtime_assert( residues.size() == 2 || residues.size() == 3 );
    core::Size const loop_start( protocols::rosetta_scripts::parse_resnum( residues[ 1 ], pose ) );
    core::Size const loop_stop( protocols::rosetta_scripts::parse_resnum( residues[ 2 ], pose ) );
    core::Size loop_cut( 0 );
	  if( residues.size() == 3 )
      loop_cut = protocols::rosetta_scripts::parse_resnum( residues[ 3 ], pose );
    runtime_assert( loop_start <= loop_stop );
	  runtime_assert( loop_start >= 1 );
    runtime_assert( loop_stop <= pose.total_residue() );
    loops_from_tag.add_loop( loop_start, loop_stop, loop_cut );
	}
	return( loops_from_tag );
}

} // loops
} // protocols
