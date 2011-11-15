// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loophash/MPI_LoopHashRefine_Master.cc
/// @brief
/// @author Mike Tyka

#define TRDEBUG TR.Debug

// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif

#include <protocols/loophash/MPI_LoopHashRefine.hh>
#include <protocols/loophash/MPI_LoopHashRefine_Master.hh>
#include <protocols/loophash/WorkUnit_LoopHash.hh>
#include <protocols/wum/WorkUnit_BatchRelax.hh>
#include <protocols/wum/WorkUnitBase.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>
#include <protocols/wum/SilentStructStore.hh>
#include <ObjexxFCL/format.hh>
/// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>
// Auto-header: duplicate removed #include <ObjexxFCL/format.hh>

#include <numeric/random/random.hh>

#ifndef _WIN32 // REQUIRED FOR WINDOWS
// AUTO-REMOVED #include <unistd.h>
// AUTO-REMOVED #include <cctype>
#endif

#include <fstream>
#include <utility/string_util.hh>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

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
#include <core/chemical/ResidueType.hh>
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
#include <core/conformation/Conformation.hh>
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
#include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymSlideInfo.fwd.hh>
#include <core/conformation/symmetry/SymSlideInfo.hh>
#include <core/conformation/symmetry/SymmData.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/fragment/BBTorsionSRFD.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Map.fwd.hh>
#include <core/id/DOF_ID_Map.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/NamedStubID.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/types.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/io/raw_data/DisulfideFile.hh>
#include <core/io/silent/EnergyNames.hh>
#include <core/io/silent/ProteinSilentStruct.tmpl.hh>
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
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pose/MiniPose.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/annotated_sequence.hh>
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
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/frag_picker/ContactTypes.hh>
#include <protocols/frag_picker/VallChunk.fwd.hh>
#include <protocols/frag_picker/VallChunk.hh>
#include <protocols/frag_picker/VallProvider.fwd.hh>
#include <protocols/frag_picker/VallProvider.hh>
#include <protocols/frag_picker/VallResidue.fwd.hh>
#include <protocols/frag_picker/VallResidue.hh>
#include <protocols/jobdist/Jobs.fwd.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/loophash/BackboneDB.hh>
#include <protocols/loophash/LocalInserter.fwd.hh>
#include <protocols/loophash/LoopHashLibrary.fwd.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LoopHashMap.hh>
#include <protocols/loophash/LoopHashSampler.fwd.hh>
#include <protocols/loops/Loop.fwd.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/match/Hit.fwd.hh>
#include <protocols/match/Hit.hh>
#include <protocols/match/SixDHasher.fwd.hh>
#include <protocols/match/SixDHasher.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverList.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/wum/MPI_WorkUnitManager.fwd.hh>
#include <protocols/wum/MPI_WorkUnitManager.hh>
#include <protocols/wum/SilentStructStore.fwd.hh>
#include <protocols/wum/WorkUnitBase.fwd.hh>
#include <protocols/wum/WorkUnitList.hh>
#include <protocols/wum/WorkUnitManager.fwd.hh>
#include <protocols/wum/WorkUnitManager.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/fixedsizearray1.fwd.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/stream_util.hh>
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
#include <utility/io/mpistream.hh>
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
#include <numeric/geometry/BoundingBox.fwd.hh>
#include <numeric/geometry/BoundingBox.hh>
#include <numeric/internal/ColPointers.hh>
#include <numeric/internal/ColVectors.hh>
#include <numeric/internal/ColsPointer.hh>
#include <numeric/internal/RowPointers.hh>
#include <numeric/internal/RowVectors.hh>
#include <numeric/internal/RowsPointer.hh>
#include <numeric/kdtree/WrappedPrimitive.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.fwd.hh>
#include <numeric/random/uniform.fwd.hh>
#include <numeric/random/uniform.hh>
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray1A.fwd.hh>
#include <ObjexxFCL/FArray1A.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray1P.fwd.hh>
#include <ObjexxFCL/FArray1P.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2A.fwd.hh>
#include <ObjexxFCL/FArray2A.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray2P.fwd.hh>
#include <ObjexxFCL/FArray2P.hh>
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
#include <ObjexxFCL/StaticIndexRange.fwd.hh>
#include <ObjexxFCL/StaticIndexRange.hh>
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
#include <ios>
#include <iosfwd>
#include <iostream>
#include <istream>
#include <iterator>
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
#include <basic/prof.hh>
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>

//Auto using namespaces
namespace std { } using namespace std; // AUTO USING NS
//Auto using namespaces end

using namespace ObjexxFCL;
using namespace ObjexxFCL::fmt;

namespace protocols {
namespace loophash {

using namespace protocols::wum;

static basic::Tracer TR("MPI.LHR.Master");

static numeric::random::RandomGenerator RG(3893251);  // <- Magic number, do not change it (and dont try and use it anywhere else)

void
MPI_LoopHashRefine_Master::set_defaults(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	max_loophash_per_structure_ = option[ OptionKeys::lh::max_loophash_per_structure ]();
	batch_relax_chunks_         = option[ OptionKeys::lh::mpi_batch_relax_chunks ]();
	batch_relax_absolute_max_   = option[ OptionKeys::lh::mpi_batch_relax_absolute_max ]();
	outbound_wu_buffer_size_    = option[ OptionKeys::lh::mpi_outbound_wu_buffer_size ]();
	loophash_split_size_        = option[ OptionKeys::lh::mpi_loophash_split_size     ]();
	library_expiry_time_        = option[ OptionKeys::lh::library_expiry_time ]();
	expire_after_rounds_        = option[ OptionKeys::lh::expire_after_rounds ]();
	mpi_master_save_score_only_ = option[ OptionKeys::lh::mpi_master_save_score_only ]();
}


void
MPI_LoopHashRefine_Master::init(){
	// Are we resuming an old job ?
	if( mpi_resume() != "" ){
		TR << "Resuming job from IDENT:  " <<  mpi_resume() << std::endl;
		load_state( mpi_resume() );
	} else {
		load_structures_from_cmdline_into_library( max_lib_size() * master_rank() );
	}
	
	// sample_weight cannot be initialized until after structures are imported, so we can check size
		load_sample_weight();
	TR << "STARTLIB: " << std::endl;
	print_library();
}



void
MPI_LoopHashRefine_Master::go()
{
	// initialize master (this is a virtual functino call and this function is overloaded by the children of this class)
	TR << "Init Master: " << mpi_rank() << std::endl;
	init();

	TR << "Master Node: Waiting for job requests..." << std::endl;
	while(true){
		// process any incoming messages such as incoming
		TRDEBUG << "Master: processing msgs.." << std::endl;
		process_incoming_msgs();

		TRDEBUG << "Master: process incoming" << std::endl;
		process_inbound_wus();

		TRDEBUG << "Master: process outbound" << std::endl;
		process_outbound_wus();

		// ok, we've done all our work, now wait until we hear from our slaves
		process_incoming_msgs( true );

		print_stats_auto();
	}
}



/// @brief figure out what to do with incoming WUs.
/// Some will be returning WUs that need to be resent others will be finished and will need
/// reintegration into the library
void
MPI_LoopHashRefine_Master::process_inbound_wus(){
	using namespace protocols::loops;

	check_library_expiry_dates();
	TRDEBUG << "Finished checking library dates"<<std::endl;
	if( inbound().size() > 0 ){
		TRDEBUG << "Processing inbound WUs on master.." << std::endl;
	}
	while( inbound().size() > 0 )
	{
		WorkUnitBaseOP  next_wu =  inbound().pop_next();
		runtime_assert( next_wu );

		// skip returning waiting WUs
		if ( next_wu->get_wu_type() == "waitwu" ) continue;

		// Upcast to a StructureModifier WU
		WorkUnit_SilentStructStore* structure_wu = dynamic_cast<  WorkUnit_SilentStructStore * > ( (WorkUnitBase*) (&(*next_wu)) );

		// If upcast was unsuccessful - warn and ignore.
		if ( structure_wu == NULL ){
			TR << "Cannot save structural data for WU: " << std::endl;
			next_wu->print( TR );
			continue;
		}

		// Otherwise extract structures and figure out what to do with them
		TRDEBUG << "Saving decoy store.. " << std::endl;
		SilentStructStore &decoys = structure_wu->decoys();

		if ( structure_wu->get_wu_type() == "loophasher" ){
			totaltime_loophash() += structure_wu->get_run_time();
			TR << "LoopHash return: " << decoys.size() << " structs in " << structure_wu->get_run_time() << "s " << " frm " << structure_wu->last_received_from() << std::endl;
			// Add the node that returned to the blacklist of all WUs with the same ssid and start_ir
			for( WorkUnitQueue::iterator iter = outbound().begin(); iter != outbound().end(); iter++ ) {
					if( (*iter)->get_wu_type() == "loophasher" ) {
			/*			// Upcast to a StructureModifier WU
						// Why use dynamic cast when we're sure of type? (copied from above)
						// So slow!
					TR << "im here1" << std::endl;
						WorkUnit_SilentStructStore* i = dynamic_cast<  WorkUnit_SilentStructStore * > ( (WorkUnitBase*) (&(*(*iter))) );
					TR << "im here2" << std::endl;
							if( ( core::Size )(*i->decoys().begin())->get_energy("ssid") == (core::Size)(*decoys.begin())->get_energy("ssid") ) {
//	i->extra_data_1() == structure_wu->extra_data_1() ) {
//	*/
								if( (*iter)->extra_data_1() == structure_wu->extra_data_1() && (*iter)->extra_data_3() == structure_wu->extra_data_3() ) {
									(*iter)->add_blacklist( structure_wu->last_received_from() );
									TRDEBUG << "Added node " << structure_wu->last_received_from() << " to blacklist of WU " << (*iter)->id() << std::endl;
							}
				  }
		  }
			n_loophash_++;
			//to_be_relaxed_.add( decoys );
			if( decoys.size() > 0 ){
				add_relax_batch( decoys );
				total_structures_ += decoys.size();
			}
		} else
		if ( structure_wu->get_wu_type() == "resultpack" ){
			decoys.all_sort_silent_scores();
			// dump structures
			TR << "Emperor sent: " << decoys.size() << " structs" << std::endl;
			print_library();
			add_structures_to_library( decoys, "add_n_limit" );
			print_library();
			// dont dump structures that came straight from emperor. 
			//dump_structures( decoys, mpi_master_save_score_only_ );
		} else
		if ( structure_wu->get_wu_type() == "batchrelax" ){
			decoys.all_sort_silent_scores();
			decoys.all_add_energy("state", 2 ); // mark structures are just having come thoruhg batchrelax
			totaltime_batchrelax_ += structure_wu->get_run_time();
			n_batchrelax_ ++;
			TR << "BatchRelax return: " << decoys.size() << " structs in " << structure_wu->get_run_time() << "s " << " frm " << structure_wu->last_received_from() << std::endl;
			add_structures_to_library( decoys );
			dump_structures( decoys, mpi_master_save_score_only_ );
		} else {
			TR.Error << "Unknown workunit received. " << std::endl;
		}


	}

	print_stats();
}



void
MPI_LoopHashRefine_Master::process_outbound_wus(){
	TRDEBUG << "Adding loophash WUs if necessary .. " << std::endl;
	if( outbound().size() < outbound_wu_buffer_size_ ){
		if ( library_central().size() == 0 ){
			TR.Error << "FATAL ERROR:  library_central_ is empty! " << std::endl;
      utility_exit_with_message( "FATAL ERROR:  library_central_ is empty! " );
		}
		// pick a random structure from the library

		core::Size finished_structures=0;
		for( SilentStructStore::iterator it = library_central().begin(); it !=  library_central().end(); it ++ ){
			if( max_loophash_per_structure_ > (*it)->get_energy("lhcount"))
			{
					TRDEBUG << "Adding: " << (*it) << "  " << (*it)->get_energy("lhcount") << std::endl;
				 (*it)->add_energy( "lhcount",  (*it)->get_energy("lhcount") + 1.0 );
				create_loophash_WUs( *it );
			}else{
					finished_structures += 1;
					TRDEBUG << "Already done: " << (*it) << "  " << (*it)->get_energy("lhcount") << std::endl;
			}
		}
		TR << "WARNING: " << finished_structures << "  " << library_central().size() << std::endl;
		if ( finished_structures >= library_central().size() ){
			TR << "WARNING: The starting structs exhausted!" << std::endl;
		}
	}

	save_state_auto();
}


void
MPI_LoopHashRefine_Master::create_loophash_WUs( const core::io::silent::SilentStructOP &start_struct ){

		runtime_assert( start_struct );
		core::pose::Pose start_pose;
		start_struct->fill_pose( start_pose );
  	core::util::switch_to_residue_type_set( start_pose, core::chemical::CENTROID);
		core::pose::set_ss_from_phipsi( start_pose );

		//refresh the sampling weight comment, as it may have changed
		// easier to do it here, copy_scores copies comments as well
		core::pose::delete_comment(start_pose,"sample_weight");
		core::pose::add_comment(start_pose,"sample_weight", sample_weight_str_);

		core::io::silent::ProteinSilentStruct pss;
		pss.fill_struct( start_pose );
		pss.copy_scores( *start_struct );

		// first cound up "round" counter - just counts how many times each structure has been
		// thorugh the loop hasher
		core::Size round = (core::Size)pss.get_energy("round");
		round++;
		pss.add_energy("round", round );
		pss.add_energy("masterid", mpi_rank() );
		pss.add_energy("parent_score", pss.get_energy("score") );

		core::Size start_ir = 1;
		core::Size end_ir = 1;
		core::Size ssid = (core::Size)pss.get_energy("ssid");

		core::Size count_wus = 0;
		for( ;start_ir< start_pose.total_residue(); start_ir+=loophash_split_size_ )
		{
			end_ir =  std::min( start_ir + loophash_split_size_ - 1, start_pose.total_residue());
			if( end_ir < start_ir) end_ir = start_ir;
		  if( start_pose.total_residue() - end_ir < loophash_split_size_ ) end_ir = start_pose.total_residue();
			TRDEBUG << "Adding a new loophash WU: " << start_ir << " - " << end_ir << ", ssid = " << ssid << std::endl;
			
			count_wus++;
			WorkUnit_LoopHashOP new_wu = new WorkUnit_LoopHash( start_ir, end_ir, ssid );
			// this is unsatisfying.. why can't i use the template ? grrr C++ thou are limited.
			new_wu->set_wu_type("loophasher");
			new_wu->decoys().add( pss );
			new_wu->clear_serial_data();
			outbound().add( new_wu );
		  if( start_pose.total_residue() - end_ir < loophash_split_size_ ) start_ir = start_pose.total_residue();
		}
		TR << "Added " << count_wus << " loophash WUs to queue. ssid=" << ssid << std::endl;

}


void
MPI_LoopHashRefine_Master::add_relax_batch( SilentStructStore &start_decoys ){
	if( start_decoys.size() == 0 ) return;
	TR << "Adding relax WUs.." << start_decoys.size() << std::endl;

	core::Size count_adds = 0;
	core::Size count_adds_b4_limit = 0;
	core::Size count_wus = 0;

	core::Size chunks  = 1 + core::Size( floor( core::Real(start_decoys.size()) / core::Real( batch_relax_chunks_ ) ) );
	core::Size batchrelax_batchsize_ = (start_decoys.size() / chunks) + 1;
	core::Size dcount=0;
	while( dcount < start_decoys.size() ){
		WorkUnit_BatchRelaxOP new_wu = new WorkUnit_BatchRelax_and_PostRescore();
		new_wu->set_wu_type("batchrelax");
		core::Size lcount=0;

		for(lcount=0; lcount < batchrelax_batchsize_; lcount++ ){
			if ( dcount < start_decoys.size() ){
				core::io::silent::SilentStructOP new_relax_structure =  start_decoys.get_struct( dcount );
				TRDEBUG << "AddRelaxStructure: " << format_silent_struct(new_relax_structure)  << std::endl;
				new_wu->decoys().add( new_relax_structure );
			}
			dcount++;
		}

		// Mix up the order
    std::random_shuffle( new_wu->decoys().begin(), new_wu->decoys().end());

		// make sure the chunk size doesnt exceed batch_relax_absolute_max_
		core::Size chunk_size = new_wu->decoys().size();
		new_wu->decoys().limit( batch_relax_absolute_max_ );

		total_structures_relax_ += new_wu->decoys().size();
		new_wu->clear_serial_data();

		count_adds += new_wu->decoys().size();
		count_adds_b4_limit += chunk_size;
		count_wus ++;
		// Relax work units have a lot of structures and fill up the queue and lead to memory crashes. Thus they get prioritized and added at the beginning of the queue!!
		outbound().push_front( new_wu );
	}

	TR << "Adding " << count_adds << "/" << count_adds_b4_limit << " structs for batchrlx. " << count_wus << " WUs" << std::endl;

}


// this goes through the library and identifies structures that have not managed to get replaced
// for some cutoff amount of time. It will send back this structure and request a new structure with the same ssid from
// the emperor.
void
MPI_LoopHashRefine_Master::check_library_expiry_dates(){
	core::Size current_time = time(NULL);

	SilentStructStore::iterator jt_last = library_central().begin();

	for( SilentStructStore::iterator jt =  library_central().begin();
									jt != library_central().end(); jt ++ )
	{
		TR.Debug << "Checking structure.." << std::endl;
		core::Size struct_time = (core::Size)(*jt)->get_energy("ltime");
		core::Size ssid        = (core::Size)(*jt)->get_energy("ssid");
		core::Size round       = (core::Size)(*jt)->get_energy("round");
		
		bool expired = false;
		// is the structure expired due to time limit ?		
		if( (int(current_time) - int(struct_time)) > (int)library_expiry_time_ ){
			expired = true;
			TR << "Structure: " << ssid << " is expired: " << int(current_time) - int(struct_time) << " > " << (int)library_expiry_time_ <<  std::endl;
		}

		// is the structure expired because it has done too many rounds ?
		if( (expire_after_rounds_ > 0) && ( round >= expire_after_rounds_ ) ){
			expired = true; 
			TR << "Structure: " << ssid << " Round:  is expired: " << round << " >= " << expire_after_rounds_ << std::endl; 
		}

		if( ! expired ){
			jt_last = jt;
			continue;
		}

		// ok, so the structure is expired. send it to the emperor and kill it. wait for a new structure to arrive

		(*jt)->add_energy("expire", (core::Size)(*jt)->get_energy("expire") + 1);
		

		// send the expired structure to the emperor (who will in due time send back a new one)
		WorkUnit_SilentStructStoreOP getnewstruct = new WorkUnit_SilentStructStore( );
		getnewstruct->set_wu_type( "getnewstruct" );
		getnewstruct->decoys().add( (*jt) );
		send_MPI_workunit( getnewstruct, 0 ); // The 0 is the MPI_RANK of the master - constant would be better here!
		
		// clear the queue of loophash WUs from previous struct to avoid false blacklisting
		// assume that false blacklisting from currently processing loophash WU is unlikely 
		
		core::Size erase_count = 0;
		for( WorkUnitQueue::iterator iter = outbound().begin(); iter != outbound().end();) {
				if( (*iter)->get_wu_type() == "loophasher" && ssid == (*iter)->extra_data_3() ) {
					TRDEBUG<<"erasing wu" <<std::endl;
					iter->reset_to_null();
					TRDEBUG<<"erasing wu from list" <<std::endl;
					iter = outbound().erase( iter );
					TRDEBUG<<"erasing done" <<std::endl;
					erase_count ++;
				} else {
						++iter;
				}
		}
		TR << "Erased " << erase_count << " deprecated WUs from outbound queue" << std::endl;
	
		// now delete this expired structure - it is now at the emperor's mercy
		library_central().erase(jt);

		TR << "Reported expired structure to emperor: - waiting for new structure" << std::endl;
		receive_MPI_workunit( 0 ); //receive the reply from master and add it to the normal inbound queue. the 0 here is the emperor's MPIRANK. Better replace with a function or constant
		TR << "Done. Restarting reporting.." << std::endl;	
		break; // only one at a time.
		// reset the iterator to the beginning - we must do that because we could have added the new structure whereever - beginning is the only save iterator 
		jt=library_central().begin();
			
		TRDEBUG << "Library state: " << std::endl;	
		print_library();
	}
	TRDEBUG << "end of check_library_expiry_dates" << std::endl;
}

/// This is a virtual over load of the base class MPI_LoopHashRefine:: add_structure_to_library with an extra behavioural step
/// that reports any successful library add-ons to the emperor. This behaviour is master specific and thus should not be in the base class.

bool
MPI_LoopHashRefine_Master::add_structure_to_library( core::io::silent::ProteinSilentStruct &pss, std::string add_algorithm ){
	bool result = MPI_LoopHashRefine::add_structure_to_library( pss, add_algorithm );
	TR << "MPI_LoopHashRefine_Master::add_structure_to_library: " << std::endl;
	if(result) report_structure_to_emperor( pss );
	return result;
}

void
MPI_LoopHashRefine_Master::report_structure_to_emperor(  core::io::silent::SilentStructOP &ss ) {
	WorkUnit_SilentStructStoreOP resultpack = new WorkUnit_SilentStructStore( );
	resultpack->set_wu_type( "resultpack" );
	resultpack->decoys().add( ss );
	send_MPI_workunit( resultpack, my_emperor() );
	TR << "Reported structure to emperor: " << format_silent_struct( ss ) << std::endl;
}

void
MPI_LoopHashRefine_Master::report_structure_to_emperor(  core::io::silent::ProteinSilentStruct &pss ) {
	WorkUnit_SilentStructStoreOP resultpack = new WorkUnit_SilentStructStore( );
	resultpack->set_wu_type( "resultpack" );
	resultpack->decoys().add( pss );
	send_MPI_workunit( resultpack, my_emperor() );
	TR << "Reported structure to emperor: " << format_silent_struct(pss) << std::endl;
}


void
MPI_LoopHashRefine_Master::load_sample_weight() {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		// This just loads sample weights from a file
		// I assume that the optionkeys sanitizes input
		
		// using ifstream instead of utility::io::izstream because izstream doesn't return success bool 
		if( option[ OptionKeys::lh::sample_weight_file ].active() ) {
				std::string pathtofile = option[ OptionKeys::lh::sample_weight_file ]();
				std::ifstream file( pathtofile.c_str() ); 
				if (!file) utility_exit_with_message( "Failed to open sample_weight file.  Check path." );
				std::string line, tmp;
				while(getline( file, line ) ) {

					boost::trim(line);
					std::vector < std::string > r;
					boost::split(r, line, boost::is_any_of("\t "));
					
				
					int i;
					// Check for correct format
					try {
						i = boost::lexical_cast<int> (r[1] );
					} catch( boost::bad_lexical_cast &) {
							utility_exit_with_message( "Sample weight second column can't be casted to an int.");
					}

					if ( ! (i >= 0 && i <=100 )) {
							utility_exit_with_message( "Sample weight second column is not an int between 0 and 100." );
						} else {
							tmp += r[1] + " ";
						}
				}
				// check for correct length
				boost::trim(tmp);
				std::list < std::string > t;
				t = utility::split_to_list(tmp);
				if( t.size() != (*(library_central().begin()))->nres() )
					utility_exit_with_message( "Sample weight file either improperly formatted or does not have same number of residues as structure." );
				TR << "Sample weight file successfully loaded" << std::endl;
				sample_weight_str_ = tmp;
		} else {
				TR << "Using default sample weight of 50 for every residue" << std::endl;
				std::string t = "50";
				for( int i = 0; i < (*(library_central().begin()))->nres() - 1; i++ ) {
						t += " 50";
				}
				sample_weight_str_ = t;
		}
}






} // namespace loophash
} // namespace protocols


