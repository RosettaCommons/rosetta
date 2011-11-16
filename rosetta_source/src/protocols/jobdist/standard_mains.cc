// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jobdist/standard_mains.cc
///
/// @brief
/// @author Ian W. Davis

#include <utility/io/izstream.hh>
#include <protocols/jobdist/JobDistributors.hh>

#include <protocols/jobdist/standard_mains.hh>

//#include <core/init.hh>
#include <core/types.hh>

// AUTO-REMOVED #include <core/sequence/util.hh>
// AUTO-REMOVED #include <core/io/pdb/file_data.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentStruct.hh>
// AUTO-REMOVED #include <core/import_pose/pose_stream/util.hh>
// AUTO-REMOVED #include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/io/silent/SilentStructFactory.hh>

// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>
#include <basic/datacache/DiagnosticData.hh>
#include <basic/Tracer.hh>

// AUTO-REMOVED #include <protocols/viewer/viewers.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/idealize/IdealizeMover.hh>
// AUTO-REMOVED #include <protocols/evaluation/RmsdEvaluator.hh>
#include <numeric/random/random_permutation.hh>
#include <utility/exit.hh>
#include <utility/file/FileName.hh>
#include <utility/pointer/owning_ptr.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/MembraneTopology.hh>

#include <core/chemical/ChemicalManager.hh>

#include <algorithm>
// AUTO-REMOVED #include <ctime>
#include <map>
#include <string>
#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif
// option key includes

// AUTO-REMOVED #include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

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
#include <core/chemical/ChemicalManager.fwd.hh>
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
// AUTO-REMOVED #include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/SecondaryStructure.fwd.hh>
// AUTO-REMOVED #include <core/fragment/SecondaryStructure.hh>
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
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/atom_tree_diffs/atom_tree_diff.hh>
#include <core/io/pdb/file_data.fwd.hh>
#include <core/io/raw_data/DisulfideFile.hh>
#include <core/io/silent/EnergyNames.hh>
#include <core/io/silent/ProteinSilentStruct.tmpl.hh>
#include <core/io/silent/SharedSilentData.hh>
#include <core/io/silent/SilentEnergy.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/io/silent/SilentStructCreator.fwd.hh>
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
#include <core/scoring/MembraneTopology.fwd.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
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
#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/evaluation/util.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/idealize/IdealizeMover.fwd.hh>
#include <protocols/jobdist/JobDistributors.fwd.hh>
#include <protocols/jobdist/Jobs.fwd.hh>
// AUTO-REMOVED #include <protocols/loops/Loop.fwd.hh>
// AUTO-REMOVED #include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.fwd.hh>
// AUTO-REMOVED #include <protocols/loops/Loops.hh>
#include <protocols/moves/DataMap.fwd.hh>
// AUTO-REMOVED #include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MoverCreator.hh>
#include <protocols/moves/MoverStatus.hh>
// AUTO-REMOVED #include <protocols/viewer/GraphicsState.hh>
// AUTO-REMOVED #include <protocols/viewer/triangle.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/down_cast.hh>
#include <utility/stream_util.hh>
#include <utility/string_util.hh>
#include <utility/vector0.fwd.hh>
#include <utility/vector0.hh>
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
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/file/gzip_util.hh>
#include <utility/io/irstream.fwd.hh>
#include <utility/io/irstream.hh>
#include <utility/io/izstream.fwd.hh>
#include <utility/io/mpistream.hh>
#include <utility/io/orstream.fwd.hh>
#include <utility/io/orstream.hh>
#include <utility/io/ozstream.fwd.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/zipstream.hpp>
#include <utility/io/zipstream.ipp>
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
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.fwd.hh>
#include <numeric/random/random.hh>
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
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/proxy_const_assert.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/ubyte.fwd.hh>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iosfwd>
#include <iostream>
#include <istream>
#include <limits>
#include <list>
#include <ostream>
#include <set>
#include <sstream>
#include <time.h>
#include <typeinfo>
#include <utility>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/datacache/CacheableData.fwd.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/CacheableString.fwd.hh>
#include <basic/datacache/DataCache.fwd.hh>
#include <basic/datacache/DataCache.hh>
#include <basic/datacache/DiagnosticData.fwd.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/prof.hh>
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/unordered/unordered_map.hpp>
#include <zlib/zlib.h>
#include <zlib/zutil.h>

//Auto Headers



namespace protocols {
namespace jobdist {

typedef utility::pointer::owning_ptr< BaseJobDistributor > BaseJobDistributorOP;

basic::Tracer TR("protocols.jobdist.main");
static numeric::random::RandomGenerator RG(32342524);
////////////////////////////////////////////////////////////////////////////////////////////////
std::string extract_tag_from_pose( core::pose::Pose &pose )
{
	//using core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG;
	using basic::datacache::CacheableString;

	if( pose.data().has( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ) ){
			CacheableString *data =  dynamic_cast< CacheableString* > (  (pose.data().get_raw_ptr( ( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG  )) ));
			if( data == NULL ) return std::string("UnknownTag");
			else               return data->str();
	}

	return std::string("UnknownTag");
}

////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< BasicJobOP > load_s_and_l()
{
	using basic::options::option;
	using utility::vector1;
	using utility::file::FileName;
	using namespace basic::options::OptionKeys;

	// concatenate -s and -l flags together to get total list of PDB files
	vector1< FileName > pdb_file_names;
	if ( option[ in::file::s ].active() )
		pdb_file_names = option[ in::file::s ]().vector(); // make a copy (-s)

	vector1< FileName > list_file_names;
	if ( option[ in::file::l ].active() )
		list_file_names = option[ in::file::l ]().vector(); // make a copy (-l)
	if ( option[ in::file::list ].active() ){
		vector1< FileName > better_list_file_names;
		better_list_file_names= option[in::file::list ]().vector(); // make a copy (-list)
		for(vector1< FileName >::iterator i = better_list_file_names.begin(), i_end = better_list_file_names.end(); i != i_end; ++i) {
			list_file_names.push_back(*i); // make a copy (-l)
		}
	}

	for(vector1< FileName >::iterator i = list_file_names.begin(), i_end = list_file_names.end(); i != i_end; ++i) {
		std::string filename( i->name() );
		utility::io::izstream data( filename.c_str() );
		if ( !data.good() ) {
			utility_exit_with_message( "Unable to open file: " + filename + '\n' );
		}
		std::string line;
		while( getline(data, line) ) {
			pdb_file_names.push_back( FileName(line) );
		}
		data.close();
	}

	vector1< BasicJobOP > jobs;
	int const nstruct_flag = option[ out::nstruct ];
	int const nstruct = std::max( 1, nstruct_flag );
	for(vector1< FileName >::iterator i = pdb_file_names.begin(), i_end = pdb_file_names.end(); i != i_end; ++i) {
		std::string native_file (i->base()+".native.pdb");
		//TR << "Looking for native: " << native_file << " ... ";
		if ( option[ in::file::native ].user() ) {
			native_file = option[ in::file::native ]().name();
		} else if ( !utility::file::file_exists(native_file) ) {
			native_file = i->name();
			//TR << " not found!" << std::endl;
		} else {
			//TR << " found!" << std::endl;
		}

		BasicJobOP job = new BasicJob(i->name(), native_file, nstruct);
		jobs.push_back( job );
	}

	return jobs;
}

////////////////////////////////////////////////////////////////////////////////////////////////
///@brief Helper function to safely get current output tag that's cached in Pose.
std::string get_output_tag(core::pose::Pose const & pose)
{
	//using core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG;
	using basic::datacache::CacheableString;

	if( !pose.data().has( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ) ) return "NO_OUTPUT_TAG_CACHED_SORRY";
	runtime_assert( dynamic_cast< CacheableString const *>( &( pose.data().get( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ))));
	return ( static_cast< CacheableString const &>(    pose.data().get( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ))).str();
}

////////////////////////////////////////////////////////////////////////////////////////////////
///@brief Helper function to safely get score_map that's cached in Pose.
std::map < std::string, core::Real > get_score_map(core::pose::Pose const & pose)
{
	//using core::pose::datacache::CacheableDataType::SCORE_MAP;
	using basic::datacache::DiagnosticData;

	if( !pose.data().has( core::pose::datacache::CacheableDataType::SCORE_MAP ) ) {
		std::map < std::string, core::Real > map;
		map[ "NO_OUTPUT_TAG_CACHED_SORRY" ] = 0.0;
		return map;
	}
	runtime_assert( dynamic_cast< DiagnosticData const *>( &( pose.data().get( core::pose::datacache::CacheableDataType::SCORE_MAP ))));
	return ( static_cast< DiagnosticData const &>(    pose.data().get( core::pose::datacache::CacheableDataType::SCORE_MAP ))).data();
}

////////////////////////////////////////////////////////////////////////////////////////////////
///@details Universal IO handler for score, relax and cluster (and hopefully more) executables
///         Handles reading/writing of Silent OR PDB files.
void register_options_universal_main(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
  option.add_relevant( in::file::centroid_input      );
  option.add_relevant( in::file::fullatom            );
  option.add_relevant( in::file::l                   );
  option.add_relevant( in::file::native              );
  option.add_relevant( in::file::s                   );
  option.add_relevant( in::file::silent              );
  option.add_relevant( in::file::silent_struct_type  );
  option.add_relevant( in::file::silent_list         );
  option.add_relevant( in::file::tags                );
  option.add_relevant( out::file::scorefile          );
  option.add_relevant( out::file::silent             );
  option.add_relevant( out::nooutput                 );
  option.add_relevant( out::nstruct                  );
  option.add_relevant( out::prefix                   );
  option.add_relevant( run::repeat                   );
}

int universal_main(
	protocols::moves::Mover & mover,
	float thinout_factor
)
{
	using namespace protocols;
	using namespace protocols::jobdist;
	using namespace protocols::moves;
	using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace utility::file;
	//using core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG;
	//using core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY;
	using basic::datacache::CacheableString;

	// open native pose if it exists
	core::pose::Pose native_pose;
	if (option[ in::file::native ].user()) {
		core::import_pose::pose_from_pdb( native_pose, option[ in::file::native ]() );
		// Set the native pose into the mover
		mover.set_native_pose( new core::pose::Pose(native_pose) );
#ifdef BOINC_GRAPHICS
		// set native for graphics
		boinc::Boinc::set_graphics_native_pose( native_pose );
#endif
	}

	/// ---------- SILENT ----------------------------------
	// Are we reading a silent file here ?
	if( option[ in::file::silent ].user() ||  option[ in::file::silent_list ].user()  ){

		// setup residue types
		core::chemical::ResidueTypeSetCAP rsd_set;
		if ( option[ in::file::fullatom ]() ) {
			rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
		} else if ( option[ in::file::residue_type_set ]() == "rna"  || option[ out::file::residue_type_set ]() == "rna"  ) {
			rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );
		} else {
			rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
		}

		// Very crude & slow way to extract specific tags:
		// Are there specific tags that are required ?

		utility::vector1< std::string > requested_tags;
		if ( option[ in::file::tags ].user() ) {
			requested_tags = option[ in::file::tags ]().vector();
		}

		utility::vector1< std::string > requested_user_tags;
		if ( option[ in::file::user_tags ].user() ) {
			requested_user_tags = option[ in::file::user_tags ]().vector();
		}

		// Grab silent filenames from -in:file:silent option
		utility::vector1< FileName > input_silent_files;
		if ( option[ in::file::silent ].active() )
			input_silent_files = option[ in::file::silent ]();

		// Grab lists of silent filenames from -in:file:silent_list option
		utility::vector1< FileName > list_file_names;
		if ( option[ in::file::silent_list ].active() )
			list_file_names = option[ in::file::silent_list ]();

		for(utility::vector1< FileName >::iterator i = list_file_names.begin(), i_end = list_file_names.end(); i != i_end; ++i) {
			std::string filename( i->name() );
			utility::io::izstream data( filename.c_str() );
			if ( !data.good() ) {	utility_exit_with_message( "Unable to open file: " + filename + '\n' );	}
			std::string line;
			while( getline(data, line) ) {
				input_silent_files.push_back( FileName(line) );
			}
			data.close();
		}

		int const nstruct_flag = option[ out::nstruct ];
		int const nstruct = std::max( 1, nstruct_flag );

//		utility::vector1< FileName >::iterator silent_file_it = input_silent_files.begin(), silent_file_it_end = input_silent_files.end();
//		// Loop over all the silent input files
//		for ( ; silent_file_it != silent_file_it_end; ++silent_file_it ) {

		if( input_silent_files.size() > 1 ){
			utility_exit_with_message( "Sorry - extracitng/reading multiple silent files is disabled.");
		}

		std::string infile  = *input_silent_files.begin();

		//	Read the silent structures
		core::io::silent::SilentFileData sfd;
		sfd.read_file( infile );

		// Grab a prefix if it exists
		std::string outfile = option[ out::prefix ]();

		if( input_silent_files.size() > 1 ) outfile = outfile + infile;

		utility::vector1< BasicJobOP > input_jobs;
		BaseJobDistributorOP jobdist;

		// create joblist
		for ( core::io::silent::SilentFileData::iterator iter = sfd.begin(), end = sfd.end(); iter != end; ++iter ) {
			if( RG.uniform() < thinout_factor ){
				//std::cout << "Thinout: Skipping " << iter->decoy_tag() << std::endl;
				continue; // ignore structures if thinout is required!
			}

			// if the user requested specific tags then check for them
			if ( option[ in::file::tags ].user() ) {
				// Attempt to find the tag in the requested tag list
				bool foundtag = false;
				utility::vector1< std::string >::iterator it = requested_tags.begin(), end = requested_tags.end();
				for ( ; it != end; ++it ) {
					if( iter->decoy_tag() ==  (*it) ){
						foundtag = true; break;
					}
				}
				if( foundtag  ) requested_tags.erase( it );
				if( !foundtag ) continue; // skip this structure if none of the tags match
			}

			// if the user requested specific tags then check for them
			if ( option[ in::file::user_tags ].user() ) {
				// Attempt to find the user_tag in the requested user_tag list
				bool founduser_tag = false;
				utility::vector1< std::string >::iterator it = requested_user_tags.begin(), end = requested_user_tags.end();
				std::cout << iter->get_comment("user_tag") << std::endl;
				for ( ; it != end; ++it ) {
					if( iter->get_comment("user_tag") == (*it) ){
						std::cerr << "FOUND: " << iter->get_comment("user_tag") << std::endl;
						founduser_tag = true; break;
					}
				}
				if( !founduser_tag ) continue; // skip this structure if none of the user_tags match
			}

			BasicJobOP job = new BasicJob( iter->decoy_tag() , "", nstruct);
			job->set_preserve_whole_input_tag( true );
			input_jobs.push_back( job );
		}

		bool const silent_output = option[ out::file::silent ].user();
		if ( silent_output ) {
		#ifdef BOINC
			std::cerr << "Silent Output Mode " << std::endl;
		#endif
			TR << "Silent Output Mode " << std::endl;
			jobdist = new PlainSilentFileJobDistributor(input_jobs);
		} else {
		#ifdef BOINC
			std::cerr << "PDB Output Mode " << std::endl;
		#endif
			TR << "PDB Output Mode " << std::endl;
			jobdist = new PlainPdbJobDistributor(input_jobs, "none");
		}
		if( option[ out::nooutput ]() ){
			jobdist->disable_output();
			jobdist->enable_ignorefinished();
		}



		evaluation::MetaPoseEvaluatorOP evaluator = new evaluation::MetaPoseEvaluator;
		evaluation::read_common_evaluator_options(*evaluator);

		BasicJobOP curr_job, prev_job;
		int curr_nstruct, num_structures_processed = 0;
		#ifdef BOINC
		std::cerr << "Jobdist startup.." << std::endl;
		#endif
		jobdist->startup();
		while( jobdist->next_job(curr_job, curr_nstruct) ) {

			// Now loop over each structure in that silent file
			for ( core::io::silent::SilentFileData::iterator iter = sfd.begin(), end = sfd.end(); iter != end; ++iter ) {

				if( iter->decoy_tag() != curr_job->input_tag() ) continue;

				time_t pdb_start_time = time(NULL);

				std::string curr_job_tag = curr_job->output_tag(curr_nstruct);
				#ifdef BOINC
				#endif
				std::cerr << "Starting work on structure: " << curr_job_tag << " <--- " << curr_job->input_tag() << std::endl;
				// SilentStruct
				core::pose::Pose input_pose;
				iter->fill_pose( input_pose, *rsd_set );
				setPoseExtraScores( input_pose, "silent_score", iter->get_energy( "score" ) );

				std::string user_tag( iter->get_comment( "user_tag" ) );
	    	if ( user_tag == "" ) user_tag = iter->get_comment( "user_ta" );
	    	if ( user_tag == "" ) user_tag = iter->get_comment( "user_t" );

				std::string alignment_id( iter->get_comment( "alignment_id" ) );

				std::cout << "USERTAGS: " << alignment_id << "  " <<  user_tag  << std::endl;

				// are we a centroid scoring function but the input is fullatom ? )
				if (( (option[ OptionKeys::score::weights ]() == "score0") ||
							(option[ OptionKeys::score::weights ]() == "score2") ||
							(option[ OptionKeys::score::weights ]() == "score3") ||
							(option[ OptionKeys::score::weights ]() == "score5") ||
							(option[ OptionKeys::score::weights ]() == "score_membrane"))
			        && ( input_pose.is_fullatom() ) ) {
					std::cout << "switching to centroid" << std::endl;
					core::util::switch_to_residue_type_set( input_pose, core::chemical::CENTROID );
				}

				// are we a centroid scoring function but the input is fullatom ? )
//				if (( (option[ OptionKeys::score::weights ]() == "score12") ||
//							(option[ OptionKeys::score::weights ]() == "standard") )
//			        && ( !input_pose.is_fullatom() ) ) {
//					std::cout << "switching to fullatom" << std::endl;
//					core::util::switch_to_residue_type_set( input_pose, core::chemical::STAMDARD );
//				}



				// Work out the tag. If we are processing more than one silent file, add the file name too!
				std::string tag = iter->decoy_tag();

				mover.set_current_tag( tag );
				TR << "Working on: " << tag << std::endl;

				core::pose::PoseOP the_pose = new core::pose::Pose( input_pose );
				if( the_pose->is_fullatom() ) core::scoring::constraints::add_fa_constraints_from_cmdline_to_pose( *the_pose );
				else                     core::scoring::constraints::add_constraints_from_cmdline_to_pose( *the_pose );
				the_pose->data().set(core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG, new CacheableString(curr_job->output_tag(curr_nstruct)));

				// Membrane protein specific scoring - only centroid score function here:
				if ( option[ OptionKeys::score::weights ]() == "score_membrane" && option[in::file::spanfile].user() && option[ in::file::centroid_input ].user() )	{
					std::string const spanfile = option[ in::file::spanfile ]();
					core::scoring::MembraneTopologyOP topologyOP = new core::scoring::MembraneTopology;
					the_pose->data().set( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY, topologyOP );
					core::scoring::MembraneTopology & topology=*( static_cast< core::scoring::MembraneTopology * >( the_pose->data().get_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY )() ));
					topology.initialize(spanfile);
				}

#ifdef BOINC_GRAPHICS
				// attach boinc graphics pose observer
				protocols::boinc::Boinc::attach_graphics_current_pose_observer( *the_pose );
#endif

				//////////////////////////////////////////////////////////////////////////////////////
				////  Maybe idealize the structure before relax ?
				if ( option[ OptionKeys::run::idealize_before_protocol ].user() ) {
					protocols::idealize::IdealizeMover idealizer;
					idealizer.fast( false );
					idealizer.apply( *the_pose );
				}

				// run the mover !!!!
				for(int repeat = 0; repeat < int(option[ run::repeat ]()); ++repeat ){
					mover.apply( *the_pose );
				}

		    // Statistics
	      core::pose::setPoseExtraScores( *the_pose, "irms",  core::scoring::CA_rmsd( input_pose, *the_pose ) );
	      if ( option[ in::file::native ].user() ){
	        core::pose::setPoseExtraScores( *the_pose, "rms",   core::scoring::native_CA_rmsd( native_pose, *the_pose ) );
	        core::pose::setPoseExtraScores( *the_pose, "srms",  core::scoring::native_CA_rmsd( native_pose, input_pose ) );
				}

				if( ! option[ out::nooutput ]() ){
					// for now, output pdbs
					if (!silent_output) {
						jobdist->dump_pose_and_map( curr_job_tag, *the_pose );    // output PDB
					} else {
						PlainSilentFileJobDistributor *jd =
						 dynamic_cast< PlainSilentFileJobDistributor * > (jobdist());
						core::io::silent::SilentStructOP pss;
						//if (option[ out::file::binary_silentfile ].user()) {
						//	pss = new core::io::silent::BinaryProteinSilentStruct;
						//} else {
						//	pss = new core::io::silent::ProteinSilentStruct;
						//}
						pss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();
						pss->fill_struct( *the_pose, curr_job_tag );
						// run PoseEvaluators
						evaluator->apply( *the_pose, curr_job_tag, *pss );
	    			if ( user_tag != "" ) pss->add_string_value("user_tag", user_tag  );
	    			if ( alignment_id != "" ) pss->add_string_value("alignment_id", alignment_id  );
						jd->dump_silent( curr_nstruct, *pss );
					}
				}

				// output scorefile if specified or not in silent output mode
				if ( option[ out::file::scorefile ].user() || !silent_output ) {
					core::io::silent::SilentStructOP ss( new core::io::silent::ScoreFileSilentStruct );
					ss->fill_struct( *the_pose );
					// run PoseEvaluators
					evaluator->apply( *the_pose, curr_job_tag, *ss );
	    		if ( user_tag != "" ) ss->add_string_value("user_tag", user_tag  );
	    		if ( alignment_id != "" ) ss->add_string_value("alignment_id", alignment_id  );

					core::io::silent::SilentFileData sfd_out;
					sfd_out.write_silent_struct( *ss, option[ out::file::scorefile ]() );
				}

				prev_job = curr_job; // pointer assignment, not a copy op
				num_structures_processed += 1;
				time_t pdb_end_time = time(NULL);
				TR << "Finished " << curr_job_tag << " in " << (long)(pdb_end_time - pdb_start_time) << " seconds." << std::endl;

			} // iterate over silent structures from a silent-file
		}
		jobdist->shutdown();

		// Announce any un-found tags !!
		for ( utility::vector1< std::string >::const_iterator it = requested_tags.begin(),
					end = requested_tags.end(); it != end; ++it ) {
			std::cerr << "WARNING: Cannot find tag: " << (*it) << std::endl;
		}

	} // SILENT done

	/// -------------------------- PDB ----------------------------------------------
	if( option[ in::file::s ].user() ||  option[ in::file::l ].user()  ){

		// are we reading PDBs ?
		time_t overall_start_time = time(NULL);
		utility::vector1< BasicJobOP > input_jobs = load_s_and_l();

		BaseJobDistributorOP jobdist;

		// pick the job distributor type based on a flag
		using basic::options::option;
		using namespace basic::options::OptionKeys;

		bool const silent_output = option[ out::file::silent ].user();
		if ( silent_output ) {
			TR << "Silent Output Mode " << std::endl;
			jobdist = new PlainSilentFileJobDistributor(input_jobs);
		} else {
			TR << "PDB Output Mode " << std::endl;
			jobdist = new PlainPdbJobDistributor(input_jobs, "none");
		}

		if( option[ out::nooutput ]() ){
			jobdist->disable_output();
			jobdist->enable_ignorefinished();
		}

		BasicJobOP curr_job, prev_job;
		int curr_nstruct, num_structures_processed = 0;
		core::pose::PoseOP input_pose; // starts NULL, coords *never* modified!
		jobdist->startup();
		while( jobdist->next_job(curr_job, curr_nstruct) ) {

			if( RG.uniform() < thinout_factor ){
				//std::cout << "Thinout: Skipping " << curr_job->output_tag(curr_nstruct) << std::endl;
				continue; // ignore structures if thinout is required!
			}

			/// Load in input pose.
			if ( utility::file::file_exists( curr_job->native_tag() ) ) {

				if ( option[ in::file::centroid_input ].user() ) {
					core::import_pose::centroid_pose_from_pdb( native_pose, curr_job->native_tag() );
				} else {
					core::import_pose::pose_from_pdb( native_pose, curr_job->native_tag() );
				}
				// Set the native pose into the mover
				mover.set_native_pose( new core::pose::Pose(native_pose) );
#ifdef BOINC_GRAPHICS
				// set native for graphics
				boinc::Boinc::set_graphics_native_pose( native_pose );
#endif
			}

			std::string curr_job_tag = curr_job->output_tag(curr_nstruct);
			mover.set_current_tag( curr_job_tag );
			time_t pdb_start_time = time(NULL);
			//std::cerr << "Starting work on PDB structure: " << curr_job_tag << " <--- " << curr_job->input_tag() << std::endl;
			TR << "Starting " << curr_job->output_tag(curr_nstruct) << " ..." << std::endl;

			// we read each PDB just once to save on disk I/O
			if( curr_job.get() != prev_job.get() || input_pose.get() == NULL ) {
				input_pose = new core::pose::Pose();
				if ( option[ in::file::centroid_input ].user() ) {
					core::import_pose::centroid_pose_from_pdb( *input_pose, curr_job->input_tag() );
				} else {
					core::import_pose::pose_from_pdb( *input_pose, curr_job->input_tag() );
				}
			}

			if( input_pose->total_residue() == 0 ){
				utility_exit_with_message( "Unable to read PDB file: " + curr_job->input_tag() + '\n' );
			}

			// Make a modifiable copy of the pose read from disk
			core::pose::PoseOP the_pose = new core::pose::Pose( *input_pose );
			if( the_pose->is_fullatom() ) core::scoring::constraints::add_fa_constraints_from_cmdline_to_pose( *the_pose );
			else                          core::scoring::constraints::add_constraints_from_cmdline_to_pose( *the_pose );

			the_pose->data().set(core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG, new CacheableString(curr_job->output_tag(curr_nstruct)));

			// Membrane protein specific scoring - only centroid score function here:
			if ( option[ OptionKeys::score::weights ]() == "score_membrane" && option[in::file::spanfile].user() && option[ in::file::centroid_input ].user() )	{
				std::string const spanfile = option[ in::file::spanfile ]();
				core::scoring::MembraneTopologyOP topologyOP = new core::scoring::MembraneTopology;
				the_pose->data().set( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY, topologyOP );
				core::scoring::MembraneTopology & topology=*( static_cast< core::scoring::MembraneTopology * >( the_pose->data().get_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY )() ));
				topology.initialize(spanfile);
			}

#ifdef BOINC_GRAPHICS
			// attach boinc graphics pose observer
			protocols::boinc::Boinc::attach_graphics_current_pose_observer( *the_pose );
#endif
				//////////////////////////////////////////////////////////////////////////////////////
				////  Maybe idealize the structure before relax ?
				if ( option[ OptionKeys::run::idealize_before_protocol ].user() ) {
					protocols::idealize::IdealizeMover idealizer;
					idealizer.fast( false );
					idealizer.apply( *the_pose );
				}

			for(int repeat = 0; repeat < int(option[ run::repeat ]()); ++repeat ){
				mover.apply( *the_pose );
			}

			// Statistics
			core::pose::setPoseExtraScores( *the_pose, "irms",  core::scoring::CA_rmsd( *input_pose, *the_pose ) );
			if ( option[ in::file::native ].user() )
				core::pose::setPoseExtraScores( *the_pose, "rms",   core::scoring::native_CA_rmsd( native_pose, *the_pose ) );

			// for now, output pdbs
			if (!silent_output) {
				jobdist->dump_pose_and_map( curr_job_tag, *the_pose );    // output PDB
			} else {
				PlainSilentFileJobDistributor *jd =
				 dynamic_cast< PlainSilentFileJobDistributor * > (jobdist());
				core::io::silent::SilentStructOP pss;
				//if (option[ out::file::binary_silentfile ].user()) {
				//	pss = new core::io::silent::BinaryProteinSilentStruct;
				//} else {
				//	pss = new core::io::silent::ProteinSilentStruct;
				//}
				pss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();
				pss->fill_struct( *the_pose, curr_job_tag );
				jd->dump_silent( curr_nstruct, *pss );
			}
			// output scorefile if specified or not in silent output mode
			if ( option[ out::file::scorefile ].user() || !silent_output ) {
				core::io::silent::SilentStructOP ss( new core::io::silent::ScoreFileSilentStruct );
				ss->fill_struct( *the_pose );
				core::io::silent::SilentFileData sfd_out;
				sfd_out.write_silent_struct( *ss, option[ out::file::scorefile ]() );
			}

			prev_job = curr_job; // pointer assignment, not a copy op
			num_structures_processed += 1;
			time_t pdb_end_time = time(NULL);
			TR << "Finished " << curr_job->output_tag(curr_nstruct) << " in " << (long)(pdb_end_time - pdb_start_time) << " seconds." << std::endl;
		} // loop over jobs and nstructs
		jobdist->shutdown();

		time_t overall_end_time = time(NULL);
		TR << "Finished all " << num_structures_processed << " structures in " << (long)(overall_end_time - overall_start_time) << " seconds." << std::endl;
		if ( num_structures_processed == 0 )
			basic::Warning() << "No structures processed.  Existing output files may have been skipped, did you mean to delete them or to use the -overwrite flag?" << std::endl;
		return 0;
	} // PDB done

	return 0;
} // universal_main

////////////////////////////////////////////////////////////////////////////////////////////////
///@details Example of how to use the job distributor that will write
/// either a pdb or a raw_data (silent) file depending on an input flag
/// Because PDBs are so big and inefficient, this function is NOT RECOMMENDED
/// for use in production environements.
/// If this function doesn't meet your needs as is, please COPY it into your own main method
/// rather than modifying it in place!
/// The goal is to keep this one as a simple-as-possible EXAMPLE for others,
/// although it will suffice for many protocols as-is.
int main_plain_mover(
	protocols::moves::Mover & mover,
	bool random_permutation
)
{
	//using core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG;
	using basic::datacache::CacheableString;

	time_t overall_start_time = time(NULL);
	utility::vector1< BasicJobOP > input_jobs = load_s_and_l();

#ifndef USEMPI
	// Reduce read contention between processes by randomizing the order in which structures are processed
	// Do not randomize, though, if job distribution is controlled by MPI
	if( random_permutation ) {
		numeric::random::random_permutation( input_jobs, numeric::random::RG );
	}
#endif

	BasicJobOP curr_job, prev_job;
	int curr_nstruct, num_structures_processed = 0;
	core::pose::PoseOP input_pose; // starts NULL, coords *never* modified!
	core::pose::PoseOP native_pose; // starts NULL, coords *never* modified!

	BaseJobDistributorOP jobdist;

	// pick the job distributor type based on a flag
	using basic::options::option;
	using namespace basic::options::OptionKeys;

 	// What on earth is "raw" ? Surely these are called silent files ?
	bool const is_raw = option[ out::file::raw ]();
	bool const silent_output = option[ out::file::silent ].user();
	if ( is_raw || silent_output ) {
		jobdist = new PlainRawJobDistributor(input_jobs, ".out");
	} else {
		jobdist = new PlainPdbJobDistributor(input_jobs, "score");
	}

	if( option[ out::nooutput ]() ){
		jobdist->disable_output();
		jobdist->enable_ignorefinished();
	}

	std::map < std::string, core::Real > score_map;

	jobdist->startup();
	while( jobdist->next_job(curr_job, curr_nstruct) ) {
		time_t pdb_start_time = time(NULL);
		TR << "Starting " << curr_job->output_tag(curr_nstruct) << " ..." << std::endl;
		jobdist->temp_file( curr_job->output_tag(curr_nstruct) );

		// we read each PDB just once to save on disk I/O
		if( curr_job.get() != prev_job.get() || input_pose.get() == NULL ) {
			input_pose = new core::pose::Pose();
			if ( option[ in::file::centroid_input ].user() ) {
				core::import_pose::centroid_pose_from_pdb( *input_pose, curr_job->input_tag() );
				native_pose = new core::pose::Pose();
				core::import_pose::centroid_pose_from_pdb( *native_pose, curr_job->native_tag() );
			} else {
				core::import_pose::pose_from_pdb( *input_pose, curr_job->input_tag() );
				native_pose = new core::pose::Pose();
				core::import_pose::pose_from_pdb( *native_pose, curr_job->native_tag() );
			}
		}
		mover.set_input_pose( input_pose );
		mover.set_native_pose( native_pose );

		// Make a modifiable copy of the pose read from disk
		core::pose::PoseOP the_pose = new core::pose::Pose( *input_pose );
		the_pose->data().set(core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG, new CacheableString(curr_job->output_tag(curr_nstruct)));

		mover.apply( *the_pose );

		prev_job = curr_job; // pointer assignment, not a copy op
		num_structures_processed += 1;
		time_t pdb_end_time = time(NULL);
		TR << "Finished " << curr_job->output_tag(curr_nstruct) << " in " << (long)(pdb_end_time - pdb_start_time) << " seconds." << std::endl;

		score_map = get_score_map( *the_pose );

		if ( option[ run::timer ].user() ){
			score_map["time"] = pdb_end_time - pdb_start_time;
			}

		jobdist->score_map( score_map );
		jobdist->dump_pose_and_map( curr_job->output_tag(curr_nstruct), *the_pose );

	} // loop over jobs and nstructs
	jobdist->shutdown();

	time_t overall_end_time = time(NULL);
	TR << "Finished all " << num_structures_processed << " structures in " << (long)(overall_end_time - overall_start_time) << " seconds." << std::endl;
	if ( num_structures_processed == 0 )
		basic::Warning() << "No structures processed.  Existing output files may have been skipped, did you mean to delete them?" << std::endl;
	return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////
///@details Example of how to use the job distributor with plain old PDB files.
/// Because PDBs are so big and inefficient, this function is NOT RECOMMENDED
/// for use in production environements.
/// If this function doesn't meet your needs as is, please COPY it into your own main method
/// rather than modifying it in place!
/// The goal is to keep this one as a simple-as-possible EXAMPLE for others,
/// although it will suffice for many protocols as-is.
int main_plain_pdb_mover(
	protocols::moves::Mover & mover,
	core::scoring::ScoreFunctionOP scorefxn
)
{
	//using core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG;
	using basic::datacache::CacheableString;

	time_t overall_start_time = time(NULL);
	utility::vector1< BasicJobOP > input_jobs = load_s_and_l();

#ifndef USEMPI
	// Reduce read contention between processes by randomizing the order in which structures are processed
	// Do not randomize, though, if job distribution is controlled by MPI
	numeric::random::random_permutation( input_jobs, numeric::random::RG );
#endif

	BaseJobDistributorOP jobdist;

	// pick the job distributor type based on a flag
	using basic::options::option;
	using namespace basic::options::OptionKeys;

	bool const silent_output = option[ out::file::silent ].user();
	if ( silent_output ) {
		TR << "Silent Output Mode " << std::endl;
		jobdist = new PlainSilentFileJobDistributor(input_jobs);
	} else {
		TR << "PDB Output Mode " << std::endl;
		jobdist = new PlainPdbJobDistributor(input_jobs, "score");
	}

	if( option[ out::nooutput ]() ){
		jobdist->disable_output();
		jobdist->enable_ignorefinished();
	}

	// load native pose (if provided)
	core::pose::Pose native_pose;
	if ( option[ in::file::native ].user() ) {
		core::import_pose::pose_from_pdb( native_pose, option[ in::file::native ]() );
	}

	BasicJobOP curr_job, prev_job;
	int curr_nstruct, num_structures_processed = 0;
	core::pose::PoseOP input_pose; // starts NULL, coords *never* modified!
	jobdist->startup();
	while( jobdist->next_job(curr_job, curr_nstruct) ) {
		time_t pdb_start_time = time(NULL);
		TR << "Starting " << curr_job->output_tag(curr_nstruct) << " ..." << std::endl;

		// we read each PDB just once to save on disk I/O
		if( curr_job.get() != prev_job.get() || input_pose.get() == NULL ) {
			input_pose = new core::pose::Pose();
			if ( option[ in::file::centroid_input ].user() ) {
				core::import_pose::centroid_pose_from_pdb( *input_pose, curr_job->input_tag() );
			} else {
				core::import_pose::pose_from_pdb( *input_pose, curr_job->input_tag() );
			}
		}

		// Make a modifiable copy of the pose read from disk
		core::pose::PoseOP the_pose = new core::pose::Pose( *input_pose );
		the_pose->data().set(core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG, new CacheableString(curr_job->output_tag(curr_nstruct)));

		for(int repeat = 0; repeat < int(option[ run::repeat ]()); ++repeat ){
			mover.apply( *the_pose );
		}

		// Score new structure.  Cached energies (including *residue* energies)
		// must be up-to-date in order to get sensible output.  If you remove these
		// lines, you *must* insert equivalent logic at the end of all apply() methods
		// (or at least for all movers that might be passed to this function).
		(*scorefxn)( *the_pose );
		/// Now handled automatically.  scorefxn->accumulate_residue_total_energies( *the_pose );

		// Statistics
		core::pose::setPoseExtraScores( *the_pose, "irms",  core::scoring::CA_rmsd( *input_pose, *the_pose ) );

		if ( option[ in::file::native ].user() )
			core::pose::setPoseExtraScores( *the_pose, "rms",   core::scoring::native_CA_rmsd( native_pose, *the_pose ) );

		// for now, output pdbs
		if (!silent_output) {
			jobdist->dump_pose_and_map( curr_job->output_tag(curr_nstruct), *the_pose );    // output PDB
		} else {
			protocols::jobdist::PlainSilentFileJobDistributor *jd =
			     dynamic_cast< protocols::jobdist::PlainSilentFileJobDistributor * > (jobdist());
			core::io::silent::ProteinSilentStruct pss;
			pss.fill_struct( *the_pose, curr_job->output_tag(curr_nstruct) );
			jd->dump_silent( curr_nstruct, pss );
		}

		prev_job = curr_job; // pointer assignment, not a copy op
		num_structures_processed += 1;
		time_t pdb_end_time = time(NULL);
		TR << "Finished " << curr_job->output_tag(curr_nstruct) << " in " << (long)(pdb_end_time - pdb_start_time) << " seconds." << std::endl;
	} // loop over jobs and nstructs
	jobdist->shutdown();

	time_t overall_end_time = time(NULL);
	TR << "Finished all " << num_structures_processed << " structures in " << (long)(overall_end_time - overall_start_time) << " seconds." << std::endl;
	if ( num_structures_processed == 0 )
		basic::Warning() << "No structures processed.  Existing output files may have been skipped, did you mean to delete them or to use the -overwrite flag?" << std::endl;
	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////
///@details Example of how to use the job distributor with silent files.
/// If this function doesn't meet your needs as is, please COPY it into your own main method
/// rather than modifying it in place!
/// The goal is to keep this one as a simple-as-possible EXAMPLE for others,
/// although it will suffice for many protocols as-is.
int main_atom_tree_diff_mover(
	protocols::moves::Mover & mover,
	core::scoring::ScoreFunctionOP scorefxn
)
{
	//using core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG;
	using basic::datacache::CacheableString;

	time_t overall_start_time = time(NULL);
	utility::vector1< BasicJobOP > input_jobs = load_s_and_l();
	// Reduce read contention between processes by randomizing the order in which structures are processed
	numeric::random::random_permutation( input_jobs, numeric::random::RG );
	AtomTreeDiffJobDistributor jobdist( input_jobs, "silent.out" );

	BasicJobOP curr_job, prev_job;
	int curr_nstruct, num_structures_processed = 0;
	core::pose::PoseOP input_pose; // starts NULL, coords *never* modified!
	jobdist.startup();
	while( jobdist.next_job(curr_job, curr_nstruct) ) {
		time_t pdb_start_time = time(NULL);
		TR << "Starting " << curr_job->output_tag(curr_nstruct) << " ..." << std::endl;

		// we read each PDB just once to save on disk I/O
		if( curr_job.get() != prev_job.get() || input_pose.get() == NULL ) {
				input_pose = new core::pose::Pose();
				core::import_pose::pose_from_pdb( *input_pose, curr_job->input_tag() );
		}

		// Make a modifiable copy of the pose read from disk
		core::pose::PoseOP the_pose = new core::pose::Pose( *input_pose );
		the_pose->data().set(core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG, new CacheableString(curr_job->output_tag(curr_nstruct)));
		mover.apply( *the_pose );

		// Score new structure and add to silent file
		std::map< std::string, core::Real > scores;
		core::import_pose::atom_tree_diffs::map_of_weighted_scores(*the_pose, *scorefxn, scores);
		jobdist.dump_pose( curr_job->output_tag(curr_nstruct), scores, *input_pose, *the_pose );

		prev_job = curr_job; // pointer assignment, not a copy op
		num_structures_processed += 1;
		time_t pdb_end_time = time(NULL);
		TR << "Finished " << curr_job->output_tag(curr_nstruct) << " in " << (long)(pdb_end_time - pdb_start_time) << " seconds." << std::endl;
	} // loop over jobs and nstructs
	jobdist.shutdown();

	time_t overall_end_time = time(NULL);
	TR << "Finished all " << num_structures_processed << " structures in " << (long)(overall_end_time - overall_start_time) << " seconds." << std::endl;
	if ( num_structures_processed == 0 )
		basic::Warning() << "No structures processed.  Existing output files may have been skipped, did you mean to delete them or use the -overwrite flag?" << std::endl;
	return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////
} // namespace jobdist
} // namespace protocols
