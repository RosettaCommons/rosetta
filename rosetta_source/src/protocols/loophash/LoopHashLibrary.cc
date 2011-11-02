// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/LoopHashLibrary.cc
/// @brief
/// @author Mike Tyka
/// @author Ken Jung

#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LoopHashMap.hh>
#include <protocols/loophash/BackboneDB.hh>
#include <protocols/loops/util.hh>
#include <protocols/loophash/Exceptions.hh>
// AUTO-REMOVED #include <protocols/frag_picker/PdbIdChunkFilter.hh>

#include <core/kinematics/Edge.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
// AUTO-REMOVED #include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <core/kinematics/FoldTree.hh>
#include <basic/options/option.hh>
#include <core/import_pose/pose_stream/util.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pose/annotated_sequence.hh>
// AUTO-REMOVED #include <core/scoring/constraints/util.hh>
// AUTO-REMOVED #include <core/scoring/constraints/CoordinateConstraint.hh>
// AUTO-REMOVED #include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>
#include <core/scoring/rms_util.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/loops/Loops.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>

// AUTO-REMOVED #include <boost/cstdint.hpp>
#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
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
#include <core/fragment/BBTorsionSRFD.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/SecondaryStructure.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Map.fwd.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/JumpID.fwd.hh>
#include <core/id/JumpID.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/TorsionID.hh>
#include <core/id/types.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/io/pdb/file_data.fwd.hh>
#include <core/io/silent/SharedSilentData.hh>
#include <core/io/silent/SilentEnergy.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/io/silent/SilentStructCreator.fwd.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/kinematics/AtomPointer.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/AtomWithDOFChange.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/Edge.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
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
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/optimization/NumericalDerivCheckResult.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/operation/TaskOperation.fwd.hh>
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
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/Func.fwd.hh>
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
#include <core/sequence/Sequence.fwd.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceProfile.fwd.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/checkpoint/CheckPointer.fwd.hh>
#include <protocols/checkpoint/CheckPointer.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/frag_picker/ContactTypes.hh>
#include <protocols/frag_picker/VallChunk.fwd.hh>
#include <protocols/frag_picker/VallChunk.hh>
#include <protocols/frag_picker/VallProvider.fwd.hh>
#include <protocols/frag_picker/VallProvider.hh>
#include <protocols/frag_picker/VallResidue.fwd.hh>
#include <protocols/frag_picker/VallResidue.hh>
#include <protocols/jobdist/Jobs.fwd.hh>
#include <protocols/loophash/LoopHashLibrary.fwd.hh>
#include <protocols/loops/Loop.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/match/Hit.fwd.hh>
#include <protocols/match/SixDHasher.fwd.hh>
#include <protocols/match/SixDHasher.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverCreator.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/relax/RelaxProtocolBase.fwd.hh>
#include <protocols/relax/RelaxProtocolBase.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/PyAssert.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/fixedsizearray1.fwd.hh>
#include <utility/fixedsizearray1.hh>
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
#include <utility/excn/EXCN_Base.hh>
#include <utility/excn/Exceptions.hh>
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
#include <utility/tag/Tag.fwd.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <numeric/geometry/BoundingBox.fwd.hh>
#include <numeric/geometry/BoundingBox.hh>
#include <numeric/internal/ColPointers.hh>
#include <numeric/internal/ColVectors.hh>
#include <numeric/internal/ColsPointer.hh>
#include <numeric/internal/RowPointers.hh>
#include <numeric/internal/RowVectors.hh>
#include <numeric/internal/RowsPointer.hh>
#include <numeric/kdtree/WrappedPrimitive.hh>
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
#include <execinfo.h>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <iterator>
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
#include <basic/datacache/CacheableData.fwd.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>
#include <boost/unordered/unordered_map.hpp>
#include <boost/unordered_map.hpp>

//Auto using namespaces
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
namespace ObjexxFCL { namespace fmt { } } using namespace ObjexxFCL::fmt; // AUTO USING NS
//Auto using namespaces end

#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif

using namespace core::scoring;
using namespace core;
using namespace core::pose;
using namespace conformation;
using namespace kinematics;
using namespace protocols::frag_picker;




namespace protocols {
namespace loophash {

	static basic::Tracer TR("LoopHashLibrary");


	LoopHashLibrary::LoopHashLibrary( const utility::vector1< core::Size > &init_sizes, const core::Size num_partitions, const core::Size assigned_num):
	  options( "dfpmin", 0.2, true , false ),
	  options2( "dfpmin", 0.02,true , false )
	{
		// create the score functions needed for the grafting process
		set_default_score_functions();
		do_sanity_check_ = true ;
		for ( core::Size i=1; i<= init_sizes.size(); ++i) hash_sizes_.push_back( init_sizes[i] );
		setup_hash_maps();
		num_partitions_ = num_partitions;
		assigned_num_ = assigned_num;
		extra_ = true;
		loopdb_range_.first = 0;
		loopdb_range_.second = 0;
		db_path_ = basic::options::option[basic::options::OptionKeys::lh::db_path]();
		assigned_string_ = "";  // This is because I don't know if an initialized string is null or empty
		// we don't want db names like "part0of20" so we increment assigned_num by one to get "part1of20"
		if ( num_partitions_ > 1 ) assigned_string_ = ".part" + utility::to_string( assigned_num + 1 ) + "of" + utility::to_string( num_partitions_);

	}


	void
	LoopHashLibrary::mem_foot_print(){
		for( std::vector< core::Size >::const_iterator it = hash_sizes_.begin(); it != hash_sizes_.end(); ++it ){
			TR << "Hash: " << *it << std::endl;
			hash_[ *it ].mem_foot_print();
		}
		TR << "BackboneDB: " << bbdb_.get_mem_foot_print() << std::endl;
	}


	void
	LoopHashLibrary::setup_hash_maps()
	{
		for( std::vector< core::Size >::const_iterator it = hash_sizes_.begin(); it != hash_sizes_.end(); ++it ){
			TR.Info << "HASHSIZE: " << *it << std::endl;
			LoopHashMap newhashmap( *it );
			hash_[ *it ] = newhashmap;
		}
	}

	LoopHashMap &
	LoopHashLibrary::gethash( core::Size size )
	{
		if( hash_.count( size ) == 1 ) return hash_[ size ];
		// and if that's not true something is wrong
		throw EXCN_Invalid_Hashmap( size );

		// We should never get here -this is just to satisfy the compiler.
		return hash_[ 0 ];
	}

	void
	LoopHashLibrary::sort() {
		for( std::vector< core::Size >::const_iterator it = hash_sizes_.begin(); it != hash_sizes_.end(); ++it ){
			hash_[ *it ].sort();
		}
	}

	void
	LoopHashLibrary::save_db()
	{
		long starttime = time(NULL);
		TR.Info << "Saving bbdb_ (BackboneDatabase) " << assigned_string_ << " with extras" << std::endl;
		bbdb_.write_db( db_path_ + "backbone" + assigned_string_ + ".db" );
		for( std::vector< core::Size >::const_iterator it = hash_sizes_.begin(); it != hash_sizes_.end(); ++it ){
			TR.Info << "Saving loopdb (LoopHashDatabase) " <<  assigned_string_ << " with loop size " << *it << std::endl;
			hash_[ *it ].write_db(db_path_ + "loopdb." + utility::to_string( *it ) + assigned_string_ +  ".db" );
		}
		long endtime = time(NULL);
		TR << "Save LoopHash Library: " << endtime - starttime << " seconds " << std::endl;
	}

	void
	LoopHashLibrary::delete_db()
	{
		long starttime = time(NULL);
		TR.Info << "Deleting database files " << assigned_string_  << std::endl;
		std::string dbstring = db_path_ + "backbone" + assigned_string_ + ".db";
		if ( remove( dbstring.c_str() ) != 0 ) throw EXCN_DB_IO_Failed( dbstring , "delete" );
		TR.Info << "bbdb deletion successful" << std::endl;
		for( std::vector< core::Size >::const_iterator it = hash_sizes_.begin(); it != hash_sizes_.end(); ++it ){
			std::string dbstring = db_path_ + "loopdb." + utility::to_string( *it ) + assigned_string_ +  ".db" ;
			if( remove( dbstring.c_str() ) != 0 ) throw EXCN_DB_IO_Failed( dbstring, "delete" );
			TR.Info << "loopdb size " <<  utility::to_string( *it ) << " deletion successful" << std::endl;
		}
		long endtime = time(NULL);
		TR << "Deleted LoopHash Library: " << endtime - starttime << " seconds " << std::endl;
	}

	void
	LoopHashLibrary::load_db()
	{
		long starttime = time(NULL);
		TR.Info << "Reading bbdb_ (BackboneDatabase) " << assigned_string_ << " with extras" << std::endl;
		bbdb_.read_db( db_path_ + "backbone" + assigned_string_  + ".db", extra_ );
		for( std::vector< core::Size >::const_iterator it = hash_sizes_.begin(); it != hash_sizes_.end(); ++it ){
			TR.Info << "Reading loopdb (LoopHashDatabase) " <<  assigned_string_ << " with loop size " << *it << std::endl;
			hash_[ *it ].read_db( db_path_ + "loopdb." + utility::to_string( *it ) + assigned_string_ + ".db" );
		}
		long endtime = time(NULL);
		TR << "Read LoopHash Library from disk: " << endtime - starttime << " seconds " << std::endl;
	}

	void
	LoopHashLibrary::load_mergeddb()
	{
		// Currently, reads a slice of the backbonedb and whatever proteins
		// are included, the loops from those proteins are loaded.
		long starttime = time(NULL);

		TR.Info << "Reading merged bbdb_ (BackboneDatabase) " << assigned_string_;
		if( extra_ ) TR.Info << " with extras";
		TR.Info << std::endl;
		// Indices of homologs is returned in homolog_map
		std::map< core::Size, bool > homolog_map;
		std::string db_filename = db_path_ + "backbone.db";
		TR.Info << "Reading " <<  db_filename << std::endl;
		bbdb_.read_db( db_filename, extra_, num_partitions_, assigned_num_, loopdb_range_, homolog_map );
		for( std::vector< core::Size >::const_iterator it = hash_sizes_.begin(); it != hash_sizes_.end(); ++it ){
			TR.Info << "Reading loopdb (LoopHashDatabase) " <<  assigned_string_ << " with loop size " << *it << std::endl;
			// pass the range to the loophashmap so it knows which loops to read
			// also pass the map of homologs
			db_filename = db_path_ + "loopdb." + utility::to_string( *it ) + ".db";
			hash_[ *it ].read_db(db_filename, loopdb_range_, homolog_map );
		}
		long endtime = time(NULL);
		TR << "Read MergedLoopHash Library from disk: " << endtime - starttime << " seconds " << std::endl;
	}


	void
	LoopHashLibrary::merge(
		LoopHashLibraryOP second_lib,
		utility::vector1< core::Real> rms_cutoffs )
	{
		long starttime = time(NULL);

		// Might want to split this function into subroutines, its kinda big

		// Concat the entire second_bbdb to the master one
		// Can later add a removal step, where proteins who aren't referenced can be removed
		core::Size index_offset;
		if( extra_ != second_lib->get_extra() ) {
			throw EXCN_bbdb_Merge_Failed( extra_, second_lib->get_extra() );
		}
		if ( ! merge_bbdb( second_lib->backbone_database(), index_offset ) ) {
			throw EXCN_bbdb_Merge_Failed( "bbdb merge failed for unknown reasons" );
		}

		TR.Debug << "BBDB concated" << std::endl;
		core::Size rms_cutoff_counter = 1;  // Because hash_sizes is using an iterator instead of an index
		for( std::vector< core::Size >::const_iterator jt = hash_sizes_.begin(); jt != hash_sizes_.end(); ++jt ){

			core::Size loop_size = *jt;
			core::Real rms_cutoff = rms_cutoffs[ rms_cutoff_counter++ ];

			LoopHashMap &hashmap = gethash( loop_size );
			LoopHashMap &second_hashmap = second_lib->gethash( loop_size );
			TR.Debug << "Hashmaps loaded for frag size " << loop_size <<std::endl;

			// do NOT use bucket interface, boost can't guarantee 1 bucket = 1 key even with rehash
			//iterate through all the members of second_loopdb
			std::pair< BackboneIndexMap::iterator, BackboneIndexMap::iterator > range;
			second_hashmap.bbdb_range( range );

			std::vector < BackboneSegment > bs_vec_;
			std::vector < LeapIndex > leap_vec_;
			boost::uint64_t key = 0;
			for( BackboneIndexMap::iterator it = range.first; it != range.second; it++ ) {
				bool same_as_last = false;
				bool add_this_ = true;

				//Now grab key of that loop
				core::Size bb_index = it->second; //technically it->first == cp.key
				LeapIndex cp = second_hashmap.get_peptide( bb_index );
				// if the key is the same as the last checked loop, keep bs_vec
				if( key == cp.key ) same_as_last = true;
				key = cp.key;

				// lookup the seconddb loop bs
				BackboneSegment bs_;
				second_lib->backbone_database().get_backbone_segment( cp.index, cp.offset, loop_size , bs_ );

				if( rms_cutoff != 0 ) {
					if( !same_as_last ) {
						//Grab loops from the main loopdb that correspond to that key
						std::vector < core::Size > leap_index_equals;
						//hashmap.lookup_withkey( key, leap_index_equals );
						hashmap.radial_lookup_withkey( key, 3, leap_index_equals );

						bs_vec_.clear();
						leap_vec_.clear();
						//Need to generate vector of equal backbone segments of the master lib for the RMS check
						for( std::vector < core::Size >::const_iterator itx = leap_index_equals.begin();
								itx != leap_index_equals.end();
								++itx ){
							core::Size bb_index_equals = *itx;
							LeapIndex cp_equals = hashmap.get_peptide( bb_index_equals );
							BackboneSegment bs_equals;
							bbdb_.get_backbone_segment( cp_equals.index, cp_equals.offset , loop_size , bs_equals );
							bs_vec_.push_back( bs_equals );
							leap_vec_.push_back( cp_equals );
							//if( cp_equals.key != key ) TR.Info<< "These keys don't match: " << cp_equals.key << " " << key << std::endl;
						}
					}
					// Now do an RMS check against every bs in the master lib bucket
					for( core::Size j = 0; j < bs_vec_.size(); j++ ) {
						core::Real BBrms = get_rmsd( bs_vec_[j], bs_ );
						if ( BBrms < rms_cutoff ) {
							// If any bs is within rms_cutoff RMS, then we don't add it
							// and also don't check any others
							add_this_ = false;
							TR.Debug << "RMS too close, skipping this frag" << std::endl;
							//if(cp.key != leap_vec_[j].key )TR.Info << "keys " << cp.key << " " <<leap_vec_[j].key << std::endl;
							continue;
						}
					}
				}
				if( add_this_ ) {
					LeapIndex leap_index;
					// Since the second bbdb was just concatenated onto the end of the master bbdb
					// new index is just the size of the master bbdb + the original ba of the loop
					leap_index.index   = index_offset + cp.index;
					leap_index.offset  = cp.offset;
					leap_index.key     = key;
					hashmap.add_leap( leap_index, key );
					// Might be better to have the slaves do RMS checks when they do their partitions
					// then the following line won't be needed
					bs_vec_.push_back( bs_ );
					leap_vec_.push_back( leap_index );
				}
			}
		}
		long endtime = time(NULL);
		TR << "Merged LoopHash Library: " << endtime - starttime << " seconds " << std::endl;
	}

	bool LoopHashLibrary::merge_bbdb( const BackboneDB & second_bbdb, core::Size & index_offset ) {
		index_offset = bbdb_.size();
		core::Size extra_key_offset = bbdb_.extra_size();
		BBData tmp;
		for(core::Size i = 0; i < second_bbdb.size(); i++ ) {
			second_bbdb.get_protein( i, tmp );
			if ( extra_ ) {
				BBExtraData tmp_extra;
				second_bbdb.get_extra_data( tmp.extra_key, tmp_extra );
				bbdb_.add_extra_data( tmp_extra );
				// need to modify the extra_key in BBData now
				tmp.extra_key = tmp.extra_key + extra_key_offset;
			}
			bbdb_.add_protein( tmp );
		}
		return true;
	}

	void
	LoopHashLibrary::create_db()
	{
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		// Either obtain structural data from a vall file
		if ( option[in::file::vall].user() ){

			// by default read extra data from Vall
			extra_ = true;

			VallProviderOP chunks = new VallProvider();
			// Use partition information generate line numbers
			core::Size vall_nlines = chunks->vallNumLines(option[in::file::vall]()[1]);
			core::Size startline = 2;
			core::Size endline = vall_nlines;

			//Use an overlap of one to avoid off-by-one errors
			if (assigned_num_ != 0 ) startline = vall_nlines * assigned_num_ / num_partitions_ - 1;
			if (assigned_num_ != (num_partitions_ - 1) ) endline = vall_nlines * ( assigned_num_ + 1 ) / num_partitions_ + 1;

			// Read Vall
			chunks->vallChunksFromLibrary(option[in::file::vall]()[1], startline, endline );

			core::Size nchunks = chunks->size();
			for( core::Size i=1; i <= nchunks; ++i ){
				// Now the total number refers to within in partition
				TR.Info << i << "/" << nchunks << " in " << assigned_string_ << std::endl;
				VallChunkOP chunk = chunks->at(i);
				core::pose::PoseOP newpose = chunk->get_pose();
				extract_data_from_pose( *newpose, chunk->size(), chunk );
			}
		}

		// also obtain data from input structures

		core::chemical::ResidueTypeSetCAP rsd_set;
		if ( option[ in::file::fullatom ]() ) {
			rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
		} else {
			rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
		}
		core::import_pose::pose_stream::MetaPoseInputStream input = core::import_pose::pose_stream::streams_from_cmd_line();
		while( input.has_another_pose() ) {
			core::pose::Pose pose;
			input.fill_pose( pose, *rsd_set );
			extract_data_from_pose(  pose  );
		}
	}



	void
	LoopHashLibrary::set_default_score_functions()
	{
		using namespace core::scoring;

		scorefxn_rama_cst.set_weight( coordinate_constraint, 0.5 );
		scorefxn_rama_cst.set_weight( rama, 1.0 );


		scorefxn_cen_cst.set_weight( coordinate_constraint, 0.05 );
		scorefxn_cen_cst.set_weight( env      , 1.0);
		scorefxn_cen_cst.set_weight( pair     , 1.0);
		scorefxn_cen_cst.set_weight( cbeta    , 1.0);
		scorefxn_cen_cst.set_weight( vdw      , 1.0);
		scorefxn_cen_cst.set_weight( rg       , 3.0);
		scorefxn_cen_cst.set_weight( cenpack  , 1.0);
		scorefxn_cen_cst.set_weight( hs_pair  , 1.0);
		scorefxn_cen_cst.set_weight( ss_pair  , 1.0);
		scorefxn_cen_cst.set_weight( rsigma   , 1.0);
		scorefxn_cen_cst.set_weight( sheet    , 1.0);
	}




	void
	LoopHashLibrary::graft_loop(
		const core::pose::Pose& src_pose,
		core::pose::Pose& tgt_pose,
		protocols::loops::Loop myloop
	)
	{

		core::optimization::MinimizerOptions options( "dfpmin", 0.2 , true , false );
		core::optimization::MinimizerOptions options2( "dfpmin", 0.02 ,true , false );


		kinematics::MoveMap final_mm;
		final_mm.set_bb(true);

		core::pose::Pose pose(tgt_pose);

		// Set up contraints
		protocols::loops::Loops exclude_region;
		exclude_region.add_loop( myloop );
		add_coordinate_constraints_to_pose( pose, tgt_pose, exclude_region );

		// copy over stretch of phi/psi/omega angles

		core::pose::transfer_phi_psi( src_pose, pose, myloop.start(), myloop.stop() );

		core::optimization::AtomTreeMinimizer().run( pose, final_mm, scorefxn_rama_cst, options );

		core::Real premin_rms = core::scoring::CA_rmsd( pose, tgt_pose );
		TR.Info << "Graft: Premin RMS: " << premin_rms << std::endl;
		TR.Info << "Graft: Min Score3 " << std::endl;
		//scorefxn_cen_cst.show( TR.Info, *newpose );
		core::optimization::AtomTreeMinimizer().run( pose, final_mm, scorefxn_cen_cst, options2 );

		transfer_phi_psi( pose, tgt_pose );
	}





	void
	LoopHashLibrary::apply( core::pose::Pose& pose )
	{
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		std::string prefix = option[ out::prefix ]();
		core::Size skim_size = option[ lh::skim_size ]();

		for(int round = 0; round < 100; round ++ ){
			//core::Size count;
			static int casecount = 0;
			core::pose::Pose opose = pose;
			std::vector< core::io::silent::SilentStructOP > lib_structs;

			TR.Info << "Loophash apply function ! " << std::endl;

			// fix any shitty backbone angles.

			// Set up contraints
			ScoreFunctionOP fascorefxn = core::scoring::getScoreFunction();
			//protocols::relax::FastRelax *qrelax = new protocols::relax::FastRelax( fascorefxn, 1 );
			protocols::relax::FastRelax *relax = new protocols::relax::FastRelax( fascorefxn,  option[ OptionKeys::relax::sequence_file ]() );

			// convert pose to centroid pose:
			core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID);
			core::pose::set_ss_from_phipsi( pose );

			core::Size starttime2 = time(NULL);
			get_all( pose, lib_structs, 1, 0, 20,1400.0, 0.5, 4.0  );
			core::Size endtime2 = time(NULL);
			TR.Info << "FOUND " << lib_structs.size() << " alternative states in time: " << endtime2 - starttime2 << std::endl;

			std::random_shuffle( lib_structs.begin(), lib_structs.end());

			std::vector< core::io::silent::SilentStructOP > select_lib_structs;

			for( core::Size k=0;k< std::min(skim_size, lib_structs.size() ) ;k++){
				select_lib_structs.push_back( lib_structs[k] );
			}

			core::pose::Pose native_pose;
			core::import_pose::pose_from_pdb( native_pose, option[ in::file::native ]() );

			{ // Save centorids
				core::io::silent::SilentFileData sfd;
				std::string silent_file_ = option[ OptionKeys::out::file::silent ]() + ".centroid.out" ;
				for( core::Size h = 0; h < select_lib_structs.size(); h++){
					core::pose::Pose rpose;
					select_lib_structs[h]->fill_pose( rpose );
					core::Real rms = scoring::CA_rmsd( native_pose, rpose );
					select_lib_structs[h]->add_energy( "round", round, 1.0 );
					select_lib_structs[h]->add_energy( "rms", rms, 1.0 );
					select_lib_structs[h]->set_decoy_tag( "S_" + string_of( round ) + "_" + string_of(  h )  );
					sfd.write_silent_struct( *(select_lib_structs[h]) , silent_file_ );
				}

			}


			core::Real bestscore = MAXIMAL_FLOAT;
			core::Size bestindex = 0;
			// Batch relax the result:

			core::Size starttime = time(NULL);
			relax->batch_apply( select_lib_structs );
			core::Size endtime = time(NULL);
			TR.Info << "Batchrelax time: " << endtime - starttime << " for " << select_lib_structs.size() << " structures " << std::endl;


			for( core::Size h = 0; h < select_lib_structs.size(); h++){
				TR.Info << "DOING: " << h << " / " << select_lib_structs.size() << std::endl;
				core::pose::Pose rpose;

				select_lib_structs[h]->fill_pose( rpose );

				//core::Real score = scoring::CA_rmsd( native_pose, rpose );
				core::Real score = (*fascorefxn)(rpose);
				TR.Info << "score: " << h << "  " << score << std::endl;

				if( score < bestscore ){
					bestscore = score;
					bestindex = h;
					pose = rpose;
				}
			}
			casecount++;
			//test_loop_sample( pose, pose.total_residue() );

			core::Real bestrms = scoring::CA_rmsd( native_pose, pose );
			TR.Info << "BESTSCORE: " << bestscore << "BESTRMS" << bestrms << std::endl;
			//pose.dump_pdb( "lhb_" + prefix + "_" + utility::to_string( round ) + ".pdb" );


			core::io::silent::SilentFileData sfd;
			std::string silent_file_ = option[ OptionKeys::out::file::silent ]();
			for( core::Size h = 0; h < select_lib_structs.size(); h++){

				if( h == bestindex ) {
					core::pose::Pose rpose;
					select_lib_structs[h]->fill_pose( rpose );
					core::Real rms = scoring::CA_rmsd( native_pose, rpose );
					select_lib_structs[h]->add_energy( "round", round, 1.0 );
					select_lib_structs[h]->add_energy( "rms", rms, 1.0 );
					select_lib_structs[h]->set_decoy_tag( "S_" + string_of( round ) + "_" + string_of(  h )  );
					sfd.write_silent_struct( *(select_lib_structs[h]) , silent_file_ );
				}
			}

		}

	}


	void
	LoopHashLibrary::apply_random(
			core::pose::Pose& pose,
			core::Size &fir,
			core::Size &fjr,
			core::Real min_rms,
			core::Real max_rms
	)
	{
		using namespace core;
		using namespace core::pose;
		using namespace conformation;
		using namespace kinematics;
		using namespace protocols::match;

		core::pose::Pose original_pose = pose;

		Size nres = pose.total_residue();
		Size ir, jr;
		//Size newpep_index = 0;

		//core::Size backbone_offset;
		//bbdb_.add_pose( pose, backbone_offset );

		int runcount=0;
		runcount++;

		fir = 0;
		fjr = 0;

		while( runcount++ < 1000 ){

			// pick a random loop length
			core::Size loop_size = hash_sizes_[ rand() % hash_sizes_.size() ];

			// pick a starting residue
			ir = rand() % (nres - loop_size-2) + 2;
			jr = ir + loop_size;
			if ( ir > nres ) continue;
			if ( jr > nres ) continue;

			// find any loops

			BackboneSegment pose_bs;
			pose_bs.read_from_pose( pose, ir, loop_size );

			Real6 t;
			if(!get_rt_over_leap( original_pose, ir, jr, t )) continue;

			LoopHashMap &hashmap = gethash( loop_size );
			std::vector < core::Size > leap_index_bucket;
			hashmap.lookup( t, leap_index_bucket );

			TR.Info << "G: " << runcount << " " << ir << "  " << jr << " " << leap_index_bucket.size() << "  " << t[1] << "  " << t[2] << "  " << t[3] << "  " << t[4] << "  " << t[5] << "  " << t[6] << std::endl;

			if( leap_index_bucket.size() == 0) continue;
			TR.Info << "B: " << leap_index_bucket.size() << std::endl;
			std::vector < core::Size > filter_leap_index_bucket;
			for(  std::vector < core::Size >::const_iterator it = leap_index_bucket.begin();
					it != leap_index_bucket.end();
					++it ){

				//LeapIndex *cp = (LeapIndex*)(*it);
				core::Size retrieve_index = (core::Size) (*it);
				LeapIndex cp = hashmap.get_peptide( retrieve_index );

				// Also retrieve the backbone structures
				BackboneSegment new_bs;
				bbdb_.get_backbone_segment( cp.index, cp.offset, hashmap.get_loop_size() , new_bs );

				core::Real BBrms = get_rmsd( pose_bs, new_bs );
				if( ( BBrms > min_rms ) && ( BBrms < max_rms ) ){
					filter_leap_index_bucket.push_back( *it );
				}
			}

			if( filter_leap_index_bucket.size() == 0) continue;

			core::Size loop_choice = rand() % filter_leap_index_bucket.size();


			// APPLY LOOP and return

			// Also retrieve the backbone structures

			fir = ir;
			fjr = jr;

			core::Size retrieve_index = (core::Size) (filter_leap_index_bucket[loop_choice]);
			LeapIndex cp = hashmap.get_peptide( retrieve_index );


			BackboneSegment new_bs;
			bbdb_.get_backbone_segment( cp.index, cp.offset, hashmap.get_loop_size() , new_bs );


			core::Real BBrms = get_rmsd( pose_bs, new_bs );
			TR.Info << "Applying: " << ir << "  " << jr << "  " << BBrms << nres << loop_size << std::endl;;
			pose_bs.print();

			new_bs.apply_to_pose( pose, ir );
			return;
		}

	}


	void
	LoopHashLibrary::get_all(
			core::pose::Pose& start_pose,
			std::vector< core::io::silent::SilentStructOP > &lib_structs,
			core::Size start_res,
			core::Size stop_res,

			core::Real min_bbrms,
			core::Real max_bbrms,
			core::Real min_rms,
			core::Real max_rms
	)
	{
		using namespace core;
		using namespace core::pose;
		using namespace conformation;
		using namespace kinematics;
		using namespace protocols::match;
		using namespace optimization;
		using namespace id;

		core::pose::Pose original_pose = start_pose;
		core::pose::Pose edit_pose = start_pose;



		kinematics::MoveMap final_mm;
		final_mm.set_bb(true);
		// setup movemap & minimisation

		//				for ( Size ii=ir; ii<= jr; ++ii ) {
		//					final_mm.set_bb( ii, true );
		//					if ( newpose->residue(ii).aa() == chemical::aa_pro ) final_mm.set( TorsionID( phi_torsion, BB, ii ), false );
		//				}

		Size nres = start_pose.total_residue();
		Size ir, jr;
		//Size newpep_index = 0;

		//core::Size backbone_offset;
		//bbdb_.add_pose( pose, backbone_offset );

		int runcount=0;
		runcount++;

		// figure out start and stop residues
		if ( stop_res == 0 ) stop_res = nres;  // to do the whole protein just set stop_res to 0
		start_res = std::min( start_res, (core::Size)2 ); // dont start before 2  - WHY ?

		for( ir = 2; ir < nres; ir ++ ){
			for( core::Size k = 0; k < hash_sizes_.size(); k ++ ){
				core::Size loop_size = hash_sizes_[ k ];

				jr = ir + loop_size;
				if ( ir > nres ) continue;
				if ( jr > nres ) continue;

				// get the rigid body transform for the current segment
				BackboneSegment pose_bs;
				pose_bs.read_from_pose( start_pose, ir, loop_size );
				Real6 t;
				if(!get_rt_over_leap( original_pose, ir, jr, t )) continue;

				// Look up the bin index of that transform in the hash map
				LoopHashMap &hashmap = gethash( loop_size );
				std::vector < core::Size > leap_index_bucket;
				hashmap.lookup( t, leap_index_bucket );

				TR.Info << "G: " << runcount << " " << ir << "  " << jr << " " << leap_index_bucket.size() << "  " << t[1] << "  " << t[2] << "  " << t[3] << "  " << t[4] << "  " << t[5] << "  " << t[6] << std::endl;



				// Now for every hit, get the internal coordinates and make a short list of replacement loops
				// according to the RMS criteria

				if( leap_index_bucket.size() == 0) continue;
				std::vector < core::Size > filter_leap_index_bucket;
				for(  std::vector < core::Size >::const_iterator it = leap_index_bucket.begin();
						it != leap_index_bucket.end();
						++it ){

					// Get the actual strucure index (not just the bin index)
					core::Size retrieve_index = (core::Size) (*it);
					LeapIndex cp = hashmap.get_peptide( retrieve_index );

					// Retrieve the actual backbone structure
					BackboneSegment new_bs;
					bbdb_.get_backbone_segment( cp.index, cp.offset, hashmap.get_loop_size() , new_bs );

					// Check the values against against any RMS limitations imposed by the caller
					core::Real BBrms = get_rmsd( pose_bs, new_bs );
					if( ( BBrms > min_bbrms ) && ( BBrms < max_bbrms ) ){
						filter_leap_index_bucket.push_back( *it );
					}
				}

				// If no loops pass the previous filter - abort
				if( filter_leap_index_bucket.size() == 0) continue;

				// Now go through the chosen loops in random order
				core::Size explore_count = 0;
				std::random_shuffle( filter_leap_index_bucket.begin(), filter_leap_index_bucket.end());
				for(  std::vector < core::Size >::const_iterator it = filter_leap_index_bucket.begin();
						it != filter_leap_index_bucket.end();
						++it ){

					explore_count ++;

					clock_t starttime = clock();


					core::Size retrieve_index = *it;
					LeapIndex cp = hashmap.get_peptide( retrieve_index );

					BackboneSegment new_bs;
					bbdb_.get_backbone_segment( cp.index, cp.offset, hashmap.get_loop_size() , new_bs );

					//core::Real BBrms = get_rmsd( pose_bs, new_bs );
					// Distance measures of end point

					/* no longer applicable since leapindex doesnt have transform info
						 if( TR.Debug.visible() ){
						 core::Real xyzdist = sqrt(sqr(t[1] - cp.vecx) + sqr(t[2] - cp.vecy) + sqr(t[3] - cp.vecz));
						 core::Real ang1 = t[4] - cp.rotx; while( ang1 > 180 ) ang1 -= 360.0; while( ang1 < -180.0 ) ang1 += 360.0;
						 core::Real ang2 = t[5] - cp.roty; while( ang2 > 180 ) ang2 -= 360.0; while( ang2 < -180.0 ) ang2 += 360.0;
						 core::Real ang3 = t[6] - cp.rotz; while( ang3 > 180 ) ang3 -= 360.0; while( ang3 < -180.0 ) ang3 += 360.0;
						 core::Real angdist = sqrt(sqr(ang1)+sqr(ang2)+sqr(ang3) );
						 TR.Info << "  X: " << xyzdist << "  " << angdist << "  " << cp.rotx << "  " << cp.roty << "  " << cp.rotz << std::endl;
						 }*/

					// set newpose
					protocols::loops::Loops exclude_region;
					exclude_region.add_loop( protocols::loops::Loop( ir, jr ) );
					core::pose::Pose newpose( start_pose );
					core::pose::transfer_phi_psi( start_pose, newpose );
					add_coordinate_constraints_to_pose( newpose, original_pose, exclude_region );
					new_bs.apply_to_pose( newpose, ir );
					//scorefxn_rama_cst.show( TR.Info, *newpose );


					// just for comparison with cut!
					//core::pose::PoseOP newpose2( new Pose( original_pose ) );
					//new_bs.apply_to_pose( *newpose2, ir, true );
					//newpose2->dump_pdb("rep_" + utility::to_string( ir ) + "_" + utility::to_string( jr ) + "_" + utility::to_string( int(xyzdist) ) + "_" + utility::to_string( int(angdist) ) + ".cut.pdb" );




					//scorefxn_rama_cst.show( TR.Info, *newpose );
					//newpose->dump_pdb("rep_" + utility::to_string( ir ) + "_" + utility::to_string( jr ) + "_" + utility::to_string( int(xyzdist) ) + "_" + utility::to_string( int(angdist) ) + ".bef.pdb" );
					AtomTreeMinimizer().run( newpose, final_mm, scorefxn_rama_cst, options );
					//scorefxn_rama_cst.show( TR.Info, *newpose );
					//newpose->dump_pdb("rep_" + utility::to_string( ir ) + "_" + utility::to_string( jr ) + "_" + utility::to_string( int(xyzdist) ) + "_" + utility::to_string( int(angdist) ) + ".aft.pdb" );
					//newpose->dump_pdb("rep_" + utility::to_string( ir ) + "_" + utility::to_string( jr ) + "_" + utility::to_string( int(xyzdist) ) + "_" + utility::to_string( int(angdist) ) + ".pdb" );

					core::Real premin_rms = core::scoring::CA_rmsd( newpose, original_pose );
					TR.Info << "Premin RMS: " << premin_rms << std::endl;
					TR.Info << "Min Score3 " << std::endl;
					//scorefxn_cen_cst.show( TR.Info, *newpose );
					AtomTreeMinimizer().run( newpose, final_mm, scorefxn_cen_cst, options2 );
					//scorefxn_cen_cst.show( TR.Info, *newpose );

					// get final RMS

					core::Real final_rms = core::scoring::CA_rmsd( newpose, original_pose );
					TR.Info << "Final RMS: " << final_rms << std::endl;
					if ( ( final_rms < max_rms ) && ( final_rms > min_rms ) ){

						core::pose::Pose mynewpose( start_pose );

						transfer_phi_psi( newpose, mynewpose );

						core::io::silent::SilentStructOP new_struct = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();
						new_struct->fill_struct( mynewpose );
						lib_structs.push_back( new_struct );
					}

					//if ( lib_structs.size() > 2  ) return;

					clock_t endtime = clock();

					TR.Info << "Clocks: " << endtime - starttime << std::endl;

				}
			}
		}

	}


	void
	LoopHashLibrary::extract_data_from_pose( core::pose::Pose& pose ){
		extract_data_from_pose( pose, pose.total_residue() );
	}

	void
	LoopHashLibrary::extract_data_from_pose( core::pose::Pose& pose, core::Size nres, protocols::frag_picker::VallChunkOP chunk )
	{
		using namespace core;
		using namespace core::pose;
		using namespace conformation;
		using namespace kinematics;
		using namespace protocols::match;

		Size ir, jr;
		//Size newpep_index = 0;
		core::Size index;
		bbdb_.add_pose( pose, nres, index, chunk );

		static int runcount=0;
		runcount++;

		for( std::vector< core::Size >::const_iterator it = hash_sizes_.begin(); it != hash_sizes_.end(); ++it ){
			TR.Info << "Setting up hash: Size:  " << *it << std::endl;
			Size loop_size = *it;

			LoopHashMap &hashmap = gethash( loop_size );
			if( loop_size + 2 > nres ) continue;
			for( ir = 2; ir < ( nres - loop_size ); ir ++ ){

				jr = ir+loop_size;

				Real6 t;
				if(!get_rt_over_leap_fast( pose, ir, jr, t )) return;
				LeapIndex leap_index;
				leap_index.index   = index;
				leap_index.offset = (ir-1)*3;


				TR.Debug << "ADD: "
					<< runcount  << " "
					<< ir  << " "
					<< jr  << " "
					<< t[1]	<<	 " "
					<< t[2]	<<	 " "
					<< t[3]	<<	 " "
					<< t[4]	<<	 " "
					<< t[5]	<<	 " "
					<< t[6] << " "
					<< leap_index.index << " "
					<< leap_index.offset;
				TR.Debug << std::endl;
				hashmap.add_leap( leap_index, t );

				BackboneSegment pose_bs;
				pose_bs.read_from_pose( pose, ir, loop_size );
			}
		}



		// reset the fold tree
		FoldTree f;
		f.add_edge( 1, pose.total_residue() , Edge::PEPTIDE );
		if( f.reorder(1) == false ){
			TR.Error << "ERROR During resetting reordering of fold tree - am ignoring this LOOP ! Cannot continue " << std::endl;
			return; // continuing leads to a segfault - instead ignore this loop !
		}
		pose.fold_tree( f );

	}




	void LoopHashLibrary::test_loop_sample( core::pose::Pose& pose, core::Size nres )
	{
		using namespace core;
		using namespace core::pose;
		using namespace conformation;
		using namespace kinematics;
		using namespace protocols::match;

		core::pose::Pose original_pose = pose;

		Size ir, jr;

		core::Size index;
		bbdb_.add_pose( pose, nres, index );

		static int runcount=0;
		runcount++;

		for( std::vector< core::Size >::const_iterator it = hash_sizes_.begin(); it != hash_sizes_.end(); ++it ){
			TR.Info << "Setting up hash: Size:  " << *it << std::endl;
			Size loop_size = *it;

			LoopHashMap &hashmap = gethash( loop_size );

			for( ir = 2; ir < ( nres - loop_size ); ir ++ ){
				jr = ir+loop_size;

				Real6 t;
				if(!get_rt_over_leap( original_pose, ir, jr, t )) continue;
				LeapIndex leap_index;
				leap_index.index = index;
				leap_index.offset = (ir-1)*3;


				TR.Info << "ADD: "
					<< runcount  << " "
					<< ir  << " "
					<< jr  << " "
					<< t[1]	<<	 " "
					<< t[2]	<<	 " "
					<< t[3]	<<	 " "
					<< t[4]	<<	 " "
					<< t[5]	<<	 " "
					<< t[6] << " "
					<< leap_index.index << " "
					<< leap_index.offset;
				TR.Info << std::endl;
				hashmap.add_leap( leap_index, t );

				BackboneSegment pose_bs;
				pose_bs.read_from_pose( pose, ir, loop_size );
				pose_bs.print();

				// sanity check one
				{
					BackboneSegment check_bs;
					bbdb_.get_backbone_segment( leap_index.index, leap_index.offset, hashmap.get_loop_size() , check_bs );
					check_bs.print();


					Pose tmp_pose = original_pose;
					TR.Info << tmp_pose.fold_tree() << std::endl;
					check_bs.apply_to_pose( tmp_pose, ir );
				}


				if( do_sanity_check_ ){
					// Now retrieve everything in that bin: - sanity check:

					std::vector < core::Size > leap_index_bucket;

					// change this to looking up with key instead of transform
					hashmap.lookup_withkey( leap_index.key, leap_index_bucket );

					core::Size sani_count = 0;

					for(  std::vector < core::Size >::const_iterator it = leap_index_bucket.begin();
							it != leap_index_bucket.end();
							++it ){
						sani_count++;

						core::Size retrieve_index = (core::Size) (*it);
						LeapIndex cp = hashmap.get_peptide( retrieve_index );

						// Also retrieve the backbone structures
						BackboneSegment new_bs;
						bbdb_.get_backbone_segment( cp.index, cp.offset, hashmap.get_loop_size() , new_bs );

						core::Real BBrms = get_rmsd( pose_bs, new_bs );
						if( BBrms > 10.0 && leap_index_bucket.size() >= 2 ){

							TR.Info << "RMS: " << BBrms << std::endl;
							TR.Info << "SANI: "
								<< runcount  << " "
								<< ir  << " "
								<< jr  << " "
								<< cp.index         << "  "
								<< cp.offset         << "  "
								<< cp.key         << "  "
								<< std::endl;

							new_bs.print();
							// construct a pose with the alternative
							Pose tmp_pose = original_pose;
							TR.Info << tmp_pose.fold_tree() << std::endl;
							new_bs.apply_to_pose( tmp_pose, ir );

							Real6 t;
							get_rt_over_leap( tmp_pose, ir, jr, t );
							TR.Info << "R6CHECKHERE: " << t[1] << " " << t[2] << " " <<t[3] << " " <<t[4] << " " <<t[5] << " " <<t[6] << std::endl;

						}
						}



					}

				}
			}

			// reset the fold tree
			FoldTree f;
			f.add_edge( 1, pose.total_residue() , Edge::PEPTIDE );
			if( f.reorder(1) == false ){
				TR.Error << "ERROR During reordering of fold tree - am ignoring this LOOP ! I am done. " << std::endl;
				return; // continuing leads to a segfault - instead ignore this loop !
			}
			pose.fold_tree( f );

	}



} // namespace loops
} // namespace protocols




