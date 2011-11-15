// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_protocols
/// @brief protocols that are specific to RNA_LoopCloser
/// @detailed
/// @author Rhiju Das

#include <protocols/rna/RNA_LoopCloser.hh>
#include <protocols/rna/RNA_Minimizer.hh>
#include <protocols/rna/RNA_ProtocolUtil.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/rna/RNA_Util.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/util.hh>


//Minimizer stuff
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
// AUTO-REMOVED #include <basic/options/util.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>

//#include <numeric/random/random.hh>

// External library headers


//C++ headers
#include <vector>
#include <string>
#include <sstream>
// AUTO-REMOVED #include <fstream>
// AUTO-REMOVED #include <ctime>

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
#include <core/chemical/VariantType.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/types.hh>
#include <core/chemical/orbitals/ICoorOrbitalData.hh>
#include <core/chemical/orbitals/OrbitalType.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <core/chemical/sdf/MolData.fwd.hh>
#include <core/chemical/sdf/MolData.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Conformation.hh>
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
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID_Map.fwd.hh>
#include <core/id/DOF_ID_Map.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/JumpID.fwd.hh>
#include <core/id/JumpID.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/SequenceMapping.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/TorsionID.hh>
#include <core/id/types.hh>
#include <core/kinematics/AtomPointer.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
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
#include <core/optimization/AtomTreeMinimizer.fwd.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/optimization/NumericalDerivCheckResult.fwd.hh>
#include <core/pose/MiniPose.fwd.hh>
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
#include <core/scoring/EnergyGraph.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/LREnergyContainer.fwd.hh>
#include <core/scoring/MinimizationData.fwd.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/Constraints.fwd.hh>
#include <core/scoring/constraints/Constraints.hh>
#include <core/scoring/constraints/DOF_Constraint.fwd.hh>
#include <core/scoring/constraints/DOF_Constraint.hh>
#include <core/scoring/constraints/Func.fwd.hh>
#include <core/scoring/constraints/Func.hh>
#include <core/scoring/constraints/FuncFactory.fwd.hh>
#include <core/scoring/constraints/HarmonicFunc.fwd.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
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
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/jobdist/Jobs.fwd.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/rna/AllowInsert.fwd.hh>
#include <protocols/rna/AllowInsert.hh>
#include <protocols/rna/RNA_MatchType.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
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
#include <utility/file/FileName.fwd.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
#include <utility/io/ozstream.fwd.hh>
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
#include <utility/tag/Tag.fwd.hh>
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
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3.fwd.hh>
#include <ObjexxFCL/FArray3.hh>
#include <ObjexxFCL/FArray3D.fwd.hh>
#include <ObjexxFCL/FArray3D.hh>
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
#include <typeinfo>
#include <utility>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableData.fwd.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/DataCache.fwd.hh>
#include <basic/datacache/DataCache.hh>
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>

//Auto using namespaces
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
//Auto using namespaces end


using namespace core;
using basic::T;

static basic::Tracer TR( "protocols.rna.rna_loop_closer" ) ;

typedef  numeric::xyzMatrix< Real > Matrix;

namespace protocols {
namespace rna {

RNA_LoopCloser::RNA_LoopCloser():
	verbose_( false ),
	NUM_ROUNDS_( 100 ),
	check_tolerance_( false ),
	ccd_tolerance_( 0.000001 ),
	absolute_ccd_tolerance_( 0.01 ),
	attempt_closure_cutoff_( 20.0 ),
	fast_scan_( false )
{
	Mover::type("RNA_LoopCloser");
}



////////////////////////////////////////////
////////////////////////////////////////////
// HEY NEED TO DEFINE A CLONER?
////////////////////////////////////////////
////////////////////////////////////////////

/// @details  Apply the RNA full atom minimizer.
///
void RNA_LoopCloser::apply( core::pose::Pose & pose, std::map< Size, Size> const & connections )
{
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// Loop through all residues and look for potential chainbreaks to close --
	// marked by CUTPOINT_LOWER and CUTPOINT_UPPER variants.
	for (Size i = 1; i <= pose.total_residue(); i++ ) {

		if ( !pose.residue( i   ).has_variant_type( chemical::CUTPOINT_LOWER )  ) continue;
		if ( !pose.residue( i+1 ).has_variant_type( chemical::CUTPOINT_UPPER )  ) continue;

		// do_fast_scan will check if this cutpoint has changed since the
		//  last movement, or is the chainbreak is already well closed,
		//  or if the chainbreak is too big to close.
		if ( fast_scan_ && !passes_fast_scan( pose, i ) )  continue;

		//TR << "Trying to close chain near cutpoint " << i << "-" << i+1 << std::endl;

		apply( pose, connections, i );

	}

}

///////////////////////////////////////////////////////////////////////
void RNA_LoopCloser::apply( core::pose::Pose & pose )
{
	std::map< Size, Size > connections_dummy;
	apply( pose, connections_dummy );
}

std::string
RNA_LoopCloser::get_name() const {
	return "RNA_LoopCloser";
}

///////////////////////////////////////////////////////////////////////////////////////////////
Real RNA_LoopCloser::apply( core::pose::Pose & pose, std::map< Size, Size > const & connections, Size const & cutpoint	)
{

	if (verbose_) TR << "Closing loop at: " << cutpoint << std::endl;
	//	TR << "Closing loop at: " << cutpoint << std::endl;

	// In the future might have a separate option, e.g., for kinematic loop closure.
	return rna_ccd_close( pose, connections, cutpoint );
}

///////////////////////////////////////////////////////////////////////////////////////////////
Real RNA_LoopCloser::apply( core::pose::Pose & pose, Size const & cutpoint	)
{

	if (verbose_) TR << "Closing loop at: " << cutpoint << std::endl;
	//	TR << "Closing loop at: " << cutpoint << std::endl;

	std::map< Size, Size > connections; /* empty*/

	// In the future might have a separate option, e.g., for kinematic loop closure.
	return rna_ccd_close( pose, connections, cutpoint );
}


////////////////////////////////////////////////////////////////////
// returns false if failure.
bool
RNA_LoopCloser::passes_fast_scan( core::pose::Pose & pose, Size const i ) const
{
	//Don't bother if there hasn't been any movement...
	core::kinematics::DomainMap domain_map;
	pose.conformation().update_domain_map( domain_map );
	if ( domain_map( i ) == domain_map( i+1 ) ) {
		if (verbose_) TR << "Skipping " << i << " due to domain map."  << std::endl;
		return false;
	}

	Real const current_dist_err =   get_dist_err( pose, i );
	//Don't bother if chain is really bad...
	if ( current_dist_err > attempt_closure_cutoff_ )		{
		if (verbose_) TR << "Cutpoint " << i << " will be tough to close: " << current_dist_err << std::endl;
		return false;
	}

	//Don't bother if chain already closed...
	if ( current_dist_err < absolute_ccd_tolerance_ )		{
		if (verbose_) TR << "Cutpoint " << i << " already pretty close to closed: " << current_dist_err << std::endl;
		return false;
	}

	return true;
}



///////////////////////////////////////////////////////////////////////////////////////////////
Real
RNA_LoopCloser::rna_ccd_close( core::pose::Pose & input_pose, std::map< Size, Size > const & connections, Size const & cutpoint ) const
{
	using namespace core::scoring::rna;
	using namespace core::id;

	if ( !input_pose.residue( cutpoint ).is_RNA() ||
			 !input_pose.residue( cutpoint+1 ).is_RNA() ) {
		utility_exit_with_message( "RNA CCD closure at "+string_of( cutpoint )+" but residues are not RNA?");
	}

	if ( !input_pose.residue( cutpoint ).has_variant_type( chemical::CUTPOINT_LOWER ) ||
			 !input_pose.residue( cutpoint+1 ).has_variant_type( chemical::CUTPOINT_UPPER ) ) {
		utility_exit_with_message( "RNA CCD closure at "+string_of( cutpoint )+" but CUTPOINT_LOWER or CUTPOINT_UPPER variants not properly set up." );
	}

	/////////////////////////////////////////////////////////////////////
	//Just to get all the atom tree connectivities right, make a little scratch pose
	// to play with.
	/////////////////////////////////////////////////////////////////////
	// In the future (if refold could be faster?!), we won't need this
	// scratch pose.

	pose::Pose pose;
	pose.append_residue_by_bond( input_pose.residue( cutpoint)   );
	pose.append_residue_by_jump( input_pose.residue( cutpoint+1), 1 );

	// This is a nice extra option -- instead of just tweaking torsions in the residue right before and after the chainbreak,
	// there are some situations in which it makes sense to tweak torsions in their base pairing partners.
	bool close_two_base_pairs = false;
	Size cutpoint_partner( 0 ), cutpoint_next_partner( 0 );

	if ( connections.find( cutpoint ) != connections.end() &&
			 connections.find( cutpoint+1 ) != connections.end() ) {

		cutpoint_partner = connections.find( cutpoint )->second;
		cutpoint_next_partner = connections.find( cutpoint+1 )->second;

		if (cutpoint_partner == cutpoint_next_partner + 1 ) {
			close_two_base_pairs = true;
		}
	}
 	if ( close_two_base_pairs ) {
		TR << "Also including some torsions in partners in CCD: " << cutpoint_partner << " " << cutpoint_next_partner << std::endl;
	}

	// This is totally hard-wired and hacky -- should be easy to fix though.
	if (close_two_base_pairs ) {
		pose.append_residue_by_jump( input_pose.residue( cutpoint_next_partner /* 9 */ ), 1 );
		pose.append_residue_by_jump( input_pose.residue( cutpoint_partner /* 10 */ ), 1 );
	}

	kinematics::FoldTree f( pose.total_residue() );
	if ( close_two_base_pairs ) {
	// This is totally hard-wired and hacky -- should be easy to fix though.
		f.new_jump( 1, 4, 1 );
		f.set_jump_atoms( 1,
											core::scoring::rna::chi1_torsion_atom( pose.residue(1) ),
											core::scoring::rna::chi1_torsion_atom( pose.residue(4) )   );
		f.new_jump( 2, 3, 2 );
		f.set_jump_atoms( 2,
											core::scoring::rna::chi1_torsion_atom( pose.residue(2) ),
											core::scoring::rna::chi1_torsion_atom( pose.residue(3) )   );

	} else {
		f.new_jump( 1, 2, 1 );
		f.set_jump_atoms( 1,
											core::scoring::rna::chi1_torsion_atom( pose.residue(1) ),
											core::scoring::rna::chi1_torsion_atom( pose.residue(2) )   );
	}

	pose.fold_tree( f );

//	pose.dump_pdb( "scratch.pdb" );

	//Vectors of the three atoms upstream or downstream of the cutpoint that need to be matched.
	utility::vector1 <Vector> upstream_xyzs, downstream_xyzs;
	Real mean_dist_err( -1.0 ), mean_dist_err_prev( 9999.9999 );
	mean_dist_err = get_chainbreak_xyz( pose, 1, upstream_xyzs, downstream_xyzs );
	// This get_chainbreak_xyz will get called repeatedly after each torsion change.

	// What torsions are in play? Note that this could be expanded (or contracted) at will.
	utility::vector1< TorsionID > tor_ids, input_tor_ids;

	// epsilon and zeta before the cutpoint.
	for (Size j = 5; j <=6; ++ j){
		tor_ids.push_back( TorsionID( 1 /*cutpoint*/, BB, j ) );
		input_tor_ids.push_back( TorsionID( cutpoint, BB, j ) );
	}

	// alpha, beta, gamma, delta after the cutpoint?
	for (Size j = 1; j <=3; ++ j){
		tor_ids.push_back( TorsionID( 2 /*cutpoint+1*/, BB, j ) );
		input_tor_ids.push_back( TorsionID( cutpoint+1, BB, j ) );
	}

	//Parin May 5, 2009
	// Also include delta (sugar pucker!) from second residue if its a free floater (cutpoints before and after)
//	if ( input_pose.residue( cutpoint+1 ).has_variant_type( chemical::CUTPOINT_LOWER )  ||
//			 input_pose.residue( cutpoint+1 ).has_variant_type( chemical::UPPER_TERMINUS )  ) {
//		Size const j = 4;
//		tor_ids.push_back( TorsionID( 2 /*cutpoint+1*/, BB, j ) );
//		input_tor_ids.push_back( TorsionID( cutpoint+1, BB, j ) );
//	}

	// Also include delta (sugar pucker!) from second residue if its a free floater (cutpoints before and after)
	//	if ( input_pose.residue( cutpoint+1 ).has_variant_type( chemical::CUTPOINT_LOWER )  ||
	//			 input_pose.residue( cutpoint+1 ).has_variant_type( chemical::UPPER_TERMINUS )  ) {
	//		Size const j = 4;
	//		tor_ids.push_back( TorsionID( 2 /*cutpoint+1*/, BB, j ) );
	//		input_tor_ids.push_back( TorsionID( cutpoint+1, BB, j ) );
	//	}


	//Need to make following user-settable.
	Size nrounds( 0 );

	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////

	while ( nrounds++ < NUM_ROUNDS_ ) {

		if (check_tolerance_ && ( (mean_dist_err_prev - mean_dist_err) < ccd_tolerance_ ) ) {
			if (verbose_ ) TR << "Reached tolerance of " << ccd_tolerance_ << " after " << nrounds << " rounds. " << std::endl;
			break;
		}
		mean_dist_err_prev = mean_dist_err;

		/////////////////////////////////////////////////////////////////////
		// One round of chain closure.
		// This doesn't have to be sequential ... allow random traversal,
		// and also biased towards "moveable torsions" (i.e., not beta or epsilon)
		// Also to do? -- check on beta/epsilon -- could perhaps make use of
		//   torsion potential as a simple and general screen.
		for (Size n = 1; n <= tor_ids.size(); n++ )
		{
			TorsionID const & tor_id( tor_ids[ n ] );
			AtomID id1,id2,id3,id4;
			pose.conformation().get_torsion_angle_atom_ids( tor_id, id1, id2, id3, id4 );

			//Is this torsion before or after the chainbreak?
			AtomID my_id;
			Real dir( 0.0 );
			if ( Size( tor_id.rsd() ) == 1 /*cutpoint*/ ){
				my_id = id3;
				dir = 1;
			} else {
				my_id = id2;
				dir = -1;
			}

			core::kinematics::tree::Atom const * current_atom ( & pose.atom_tree().atom( my_id ) );

			kinematics::Stub const & stub_i( current_atom->get_stub() );
			Matrix const & M_i( stub_i.M );
			Vector const & x_i = M_i.col_x();

			Real weighted_sine( 0.0 ), weighted_cosine( 0.0 );
			for (Size m = 1; m <= upstream_xyzs.size(); m++ ){
				Vector const current_xyz( current_atom->xyz() );

				Vector const r1 = upstream_xyzs[m] - current_xyz;
				Vector const rho1 = r1 - dot( r1, x_i) * x_i;

				Vector const r2 = downstream_xyzs[m] - current_xyz;
				Vector const rho2 = r2 - dot( r2, x_i) * x_i;

				Real const current_sine   = dir * dot( x_i, cross( rho1, rho2 ) );
				Real const current_cosine = dot( rho1, rho2 );
				//				std::cout << "PREFERRED ANGLE: " << numeric::conversions::degrees( std::atan2( current_sine, current_cosine) ) << std::endl;

				weighted_sine += current_sine;
				weighted_cosine += current_cosine;

				mean_dist_err = get_chainbreak_xyz( pose, 1 /*cutpoint*/, upstream_xyzs, downstream_xyzs );
		}

			Real const twist_torsion = numeric::conversions::degrees( std::atan2( weighted_sine, weighted_cosine) );

			//			std::cout << "CHECK: " << twist_torsion << std::endl;

			Real const current_val = pose.torsion( tor_id );
			pose.set_torsion( tor_id, current_val + twist_torsion );


		}

		if ( verbose_ ) std::cout << "Distance error: " << mean_dist_err << std::endl;

	}

	if (verbose_) pose.dump_pdb( "scratch_close.pdb" );

	/////////////////////////////////////////////////////////////////////
	// OK, done with mini_pose ... copy torsions back into main pose.
	for (Size n = 1; n <= tor_ids.size(); n++ ) {
		input_pose.set_torsion( input_tor_ids[n], pose.torsion( tor_ids[n] ) );
	}

	if (verbose_) input_pose.dump_pdb( "pose_close.pdb" );

	return mean_dist_err;

}



///////////////////////////////////////////////////////////////
Real
RNA_LoopCloser::get_dist_err( pose::Pose & pose,
										Size const cutpoint
										) const
{
	utility::vector1< Vector > upstream_xyzs;
	utility::vector1< Vector > downstream_xyzs;
	return get_chainbreak_xyz( pose, cutpoint, upstream_xyzs, downstream_xyzs );
}

///////////////////////////////////////////////////////////////
Real
RNA_LoopCloser::get_chainbreak_xyz( pose::Pose & pose,
										Size const cutpoint,
										utility::vector1< Vector > & upstream_xyzs,
										utility::vector1< Vector > & downstream_xyzs
										) const
{
	upstream_xyzs.clear();
	downstream_xyzs.clear();

	upstream_xyzs.push_back( pose.residue( cutpoint ).xyz( " O3*" ) );
	upstream_xyzs.push_back( pose.residue( cutpoint ).xyz( "OVL1" ) );
	upstream_xyzs.push_back( pose.residue( cutpoint ).xyz( "OVL2" ) );

	downstream_xyzs.push_back( pose.residue( cutpoint+1 ).xyz( "OVU1" ) );
	downstream_xyzs.push_back( pose.residue( cutpoint+1 ).xyz( " P  " ) );
	downstream_xyzs.push_back( pose.residue( cutpoint+1 ).xyz( " O5*" ) );

	Real mean_dist_err( 0.0 );
	for (Size m = 1; m <= upstream_xyzs.size(); m++ ){
		mean_dist_err += ( upstream_xyzs[m]  - downstream_xyzs[m] ).length();
	}
	mean_dist_err /= upstream_xyzs.size();

	return mean_dist_err;
}

///////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
RNA_LoopCloser::get_extra_cutpoints( pose::Pose const & pose ) const
{

	using namespace core::kinematics;
	using namespace core::scoring::constraints;

	FoldTree const & f( pose.fold_tree() );

	utility::vector1< Size > extra_cutpoints;

	for ( Size cutpos = 1; cutpos < pose.total_residue(); cutpos++ ){

		if ( ! f.is_cutpoint( cutpos ) ) continue;

		if ( !pose.residue( cutpos   ).has_variant_type( chemical::CUTPOINT_LOWER )  ||
				 !pose.residue( cutpos+1 ).has_variant_type( chemical::CUTPOINT_UPPER )  ) {

			// go through pose constraints -- was there any constraint holding residue i's O3*  to i+1's P?
			bool cst_found( false );
			ConstraintCOPs csts( pose.constraint_set()->get_all_constraints() );

			for ( Size n = 1; n <= csts.size(); n++ ) {

				ConstraintCOP const & cst( csts[n] );

				if ( cst->natoms() == 2 )  { // currently only defined for pairwise distance constraints.
					Size const i = cst->atom( 1 ).rsd();
					std::string const name1 = pose.residue( i ).atom_name( cst->atom( 1 ).atomno() );
					Size const j = cst->atom( 2 ).rsd();
					std::string const name2 = pose.residue( j ).atom_name( cst->atom( 2 ).atomno() );

					if ( i == cutpos && j == cutpos+1 && name1 == " O3*" && name2 == " P  " ){
						cst_found = true; break;
					}
					if ( j == cutpos && i == cutpos+1 && name2 == " O3*" && name1 == " P  " ){
						cst_found = true; break;
					}

				}

			}

			if (!cst_found) continue;

			extra_cutpoints.push_back( cutpos );
		}
	}

	return extra_cutpoints;
}

////////////////////////////////////////////////////////////////////////////////////////
void
RNA_LoopCloser::setup_variants_at_extra_cutpoints( pose::Pose & pose, utility::vector1< Size > const & extra_cutpoints ) const
{
	pose::Pose pose_copy( pose );

	for ( Size n = 1; n <= extra_cutpoints.size(); n++ ) {

		Size cutpos( extra_cutpoints[ n ] );
		//std::cout << "ADDING CUTPOINT VARIANTS " << cutpos << std::endl;

		core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, cutpos   );
		core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, cutpos+1 );

		for (Size i = cutpos; i <= cutpos + 1; i++ ){
			for (Size j = 1; j <= scoring::rna::NUM_RNA_MAINCHAIN_TORSIONS; j++ ) {
				id::TorsionID torsion_id( i, id::BB, j );
				pose.set_torsion( torsion_id, pose_copy.torsion( torsion_id ) ) ;
			} // j
		} // i

	} // n
}



///////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_LoopCloser::remove_variants_at_extra_cutpoints( pose::Pose & pose, utility::vector1< Size > const & extra_cutpoints ) const
{

	pose::Pose pose_copy( pose );

	for ( Size n = 1; n <= extra_cutpoints.size(); n++ ) {
		Size cutpos( extra_cutpoints[ n ] );
		core::pose::remove_variant_type_from_pose_residue( pose, chemical::CUTPOINT_LOWER, cutpos   );
		core::pose::remove_variant_type_from_pose_residue( pose, chemical::CUTPOINT_UPPER, cutpos+1 );

		for (Size i = cutpos; i <= cutpos + 1; i++ ){
			for (Size j = 1; j <= scoring::rna::NUM_RNA_MAINCHAIN_TORSIONS; j++ ) {
				id::TorsionID torsion_id( i, id::BB, j );
				pose.set_torsion( torsion_id, pose_copy.torsion( torsion_id ) ) ;
			} // j
		} // i

	} // n

}

///////////////////////////////////////////////////////////////////////////////////////////////
// void
// RNA_LoopCloser::close_loops_carefully(
//   core::pose::Pose & pose,
//   core::scoring::ScoreFunctionOP const & scorefxn,
// 	Size close_loops_rounds )
// {

// 	for (Size k = 1; k <= close_loops_rounds; k++ ) {
// 		std::cout << "ROUND " << k << " of " << close_loops_rounds << ": " << std::endl;
// 		scorefxn->show( std::cout, pose );
// 		close_loops_carefully_one_round( pose, scorefxn );
// 	}

// }

///////////////////////////////////////////////////////////////////////////////////////////////
// void
// RNA_LoopCloser::close_loops_carefully(
//   core::pose::Pose & pose,
//   core::scoring::ScoreFunctionOP const & scorefxn )
// {
// 	close_loops_carefully( pose, scorefxn, 1 );
// }

///////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_LoopCloser::close_loops_carefully( core::pose::Pose & pose,  std::map< Size, Size > const & connections )
{
	using namespace core::scoring;

	//	ScoreFunctionOP scorefxn_local( scorefxn->clone() );

	// A quick minimization
	//static protocols::rna::RNA_Minimizer rna_minimizer;
	//	rna_minimizer.set_score_function( scorefxn_local );
	//	rna_minimizer.apply( pose );
	//	pose.dump_pdb( "dump1.pdb" );

	// OK, close the loops...
	// Look through constraints -- there may be some cutpoints to close that are not
	// yet flagged... set up chainbreak variants there...
	utility::vector1< Size > const extra_cutpoints = get_extra_cutpoints( pose );
	setup_variants_at_extra_cutpoints( pose, extra_cutpoints );

	// Now locally minimize around chainbreaks -- ignore chainbreak atoms because they're
	// probably in crazy places (wait until CCD).
	//scorefxn_local->set_weight( linear_chainbreak, 0.0 );
	//	local_minimize_at_chainbreaks( pose, scorefxn_local );
	//	pose.dump_pdb( "dump2.pdb" );

	// CCD loop close.
	//	bool fast_scan_save( fast_scan_ );
	//	fast_scan_ = false;
	apply( pose, connections );
	//	fast_scan_ = fast_scan_save;
	//	pose.dump_pdb( "dump3.pdb" );


	// Now minimize with strong coordinate constraints.
	tight_minimize( pose );

	//	local_minimize_at_chainbreaks( pose, scorefxn_local );

	// Leave the pose sequence/variants as we received it.
	remove_variants_at_extra_cutpoints( pose, extra_cutpoints );
	//	pose.dump_pdb( "dump4.pdb" );

	//In case the pose is about to be output ...
	//	pose.constraint_set( cst_set_save );
	//	(*scorefxn)( pose );

}

//////////////////////////////////////////////////////////////////////////////////
void
RNA_LoopCloser::tight_minimize( core::pose::Pose & pose	) const
{
	using namespace core::optimization;
	using namespace core::kinematics;
	using namespace core::scoring;

	core::scoring::constraints::ConstraintSetOP const & cst_set_save( pose.constraint_set()->clone() );
	pose.constraint_set( 0 );
	core::scoring::constraints::add_coordinate_constraints( pose );

	static AtomTreeMinimizer minimizer;
	float const dummy_tol( 0.0000025);
	bool const use_nblist( true );
	static MinimizerOptions options( "dfpmin", dummy_tol, use_nblist, false /*deriv_check_*/, false /*deriv_check_*/ );
	static ScoreFunctionOP lores_scorefxn( ScoreFunctionFactory::create_score_function( core::scoring::RNA_LORES_WTS ) );
	lores_scorefxn->set_weight( linear_chainbreak, 5.0 );
	lores_scorefxn->set_weight( coordinate_constraint, 2.0 );

	MoveMap mm;
	mm.set_bb( true );
	mm.set_jump ( false );
	mm.set_chi( false );
	minimizer.run( pose, mm, *lores_scorefxn, options );

	pose.constraint_set( cst_set_save );
}

/////////////////////////////////////////////////////////////////////////////////
void
RNA_LoopCloser::local_minimize_at_chainbreaks(
  core::pose::Pose & pose,
  core::scoring::ScoreFunctionOP & scorefxn_local ) const
{

	using namespace core::optimization;

	AtomTreeMinimizer minimizer;
	float const dummy_tol( 0.0000025);
	bool const use_nblist( true );
	MinimizerOptions options( "dfpmin", dummy_tol, use_nblist, false /*deriv_check_*/, false /*deriv_check_*/ );

	kinematics::MoveMap mm;
	mm.set_bb( false );
	mm.set_chi( false );
	mm.set_jump( false );

	scorefxn_local->set_weight( core::scoring::linear_chainbreak, 1.0 );

	kinematics::FoldTree const & f( pose.fold_tree() );

	for ( Size cutpos = 1; cutpos < pose.total_residue(); cutpos++ ){

		mm.set_bb( false );

		if ( ! f.is_cutpoint( cutpos ) ) continue;

		if ( !pose.residue( cutpos   ).has_variant_type( chemical::CUTPOINT_LOWER )  ||
				 !pose.residue( cutpos+1 ).has_variant_type( chemical::CUTPOINT_UPPER )  ) continue;

		mm.set_bb( cutpos, true );
		mm.set_bb( cutpos+1, true );
		minimizer.run( pose, mm, *scorefxn_local, options );

		TR << "Trying to minimize chain near cutpoint " << cutpos << std::endl;
	}
}


} // namespace rna
} // namespace protocols
