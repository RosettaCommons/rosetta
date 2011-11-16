// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief  read and write XEASY format peak-lists
/// @author Oliver Lange
#include <protocols/noesy_assign/CrossPeakList.hh>
#include <protocols/noesy_assign/ResonanceList.hh>
// AUTO-REMOVED #include <protocols/noesy_assign/PeakAssignmentResidueMap.hh>
// AUTO-REMOVED #include <protocols/noesy_assign/PeakAssignment.hh>
#include <protocols/noesy_assign/PeakAssignmentParameters.hh>
#include <protocols/noesy_assign/NoesyModule.hh>
#include <protocols/noesy_assign/PeakFileFormat.hh>
//#include <devel/NoesyAssign/NoeNetwork.hh>
#include <devel/init.hh>

#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <core/io/silent/SilentFileData.hh>

// AUTO-REMOVED #include <core/fragment/FragSet.hh>
// AUTO-REMOVED #include <core/fragment/FragmentIO.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>

// for switching residue type set to centroid

#include <core/chemical/ChemicalManager.fwd.hh> //required for type-set-name FA_STANDARD


#include <basic/prof.hh>

// AUTO-REMOVED #include <utility/excn/Exceptions.hh>
// AUTO-REMOVED #include <utility/io/izstream.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>
#include <basic/Tracer.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
// AUTO-REMOVED #include <core/chemical/Adduct.fwd.hh>
// AUTO-REMOVED #include <core/chemical/Adduct.hh>
// AUTO-REMOVED #include <core/chemical/AtomICoor.fwd.hh>
// AUTO-REMOVED #include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomType.hh>
// AUTO-REMOVED #include <core/chemical/AtomTypeSet.fwd.hh>
// AUTO-REMOVED #include <core/chemical/ElementSet.fwd.hh>
// AUTO-REMOVED #include <core/chemical/MMAtomType.fwd.hh>
// AUTO-REMOVED #include <core/chemical/MMAtomTypeSet.fwd.hh>
// AUTO-REMOVED #include <core/chemical/ResConnID.fwd.hh>
// AUTO-REMOVED #include <core/chemical/ResConnID.hh>
// AUTO-REMOVED #include <core/chemical/ResidueConnection.fwd.hh>
// AUTO-REMOVED #include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.fwd.hh>
// AUTO-REMOVED #include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
// AUTO-REMOVED #include <core/chemical/VariantType.fwd.hh>
#include <core/chemical/types.hh>
// AUTO-REMOVED #include <core/chemical/orbitals/ICoorOrbitalData.hh>
// AUTO-REMOVED #include <core/chemical/orbitals/OrbitalType.fwd.hh>
// AUTO-REMOVED #include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
// AUTO-REMOVED #include <core/chemical/sdf/MolData.fwd.hh>
// AUTO-REMOVED #include <core/chemical/sdf/MolData.hh>
// AUTO-REMOVED #include <core/conformation/Atom.fwd.hh>
// AUTO-REMOVED #include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
// AUTO-REMOVED #include <core/conformation/Conformation.hh>
// AUTO-REMOVED #include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/conformation/orbitals/OrbitalXYZCoords.hh>
// AUTO-REMOVED #include <core/conformation/signals/ConnectionEvent.fwd.hh>
// AUTO-REMOVED #include <core/conformation/signals/ConnectionEvent.hh>
// AUTO-REMOVED #include <core/conformation/signals/GeneralEvent.fwd.hh>
// AUTO-REMOVED #include <core/conformation/signals/GeneralEvent.hh>
// AUTO-REMOVED #include <core/conformation/signals/IdentityEvent.fwd.hh>
// AUTO-REMOVED #include <core/conformation/signals/IdentityEvent.hh>
// AUTO-REMOVED #include <core/conformation/signals/LengthEvent.fwd.hh>
// AUTO-REMOVED #include <core/conformation/signals/LengthEvent.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
// AUTO-REMOVED #include <core/conformation/signals/XYZEvent.hh>
// AUTO-REMOVED #include <core/fragment/BaseCacheUnit.hh>
// AUTO-REMOVED #include <core/fragment/FragData.fwd.hh>
// AUTO-REMOVED #include <core/fragment/FragData.hh>
// AUTO-REMOVED #include <core/fragment/FragID.fwd.hh>
// AUTO-REMOVED #include <core/fragment/FragID_Iterator.fwd.hh>
// AUTO-REMOVED #include <core/fragment/FragSet.fwd.hh>
// AUTO-REMOVED #include <core/fragment/Frame.fwd.hh>
// AUTO-REMOVED #include <core/fragment/Frame.hh>
// AUTO-REMOVED #include <core/fragment/FrameIterator.fwd.hh>
// AUTO-REMOVED #include <core/fragment/FrameIterator.hh>
// AUTO-REMOVED #include <core/fragment/FrameIteratorWorker_.fwd.hh>
// AUTO-REMOVED #include <core/fragment/FrameIteratorWorker_.hh>
// AUTO-REMOVED #include <core/fragment/FrameList.fwd.hh>
// AUTO-REMOVED #include <core/fragment/FrameList.hh>
// AUTO-REMOVED #include <core/fragment/SingleResidueFragData.fwd.hh>
// AUTO-REMOVED #include <core/fragment/SingleResidueFragData.hh>
#include <core/id/AtomID.fwd.hh>
// AUTO-REMOVED #include <core/id/AtomID.hh>
// AUTO-REMOVED #include <core/id/AtomID_Map.fwd.hh>
// AUTO-REMOVED #include <core/id/AtomID_Map.hh>
// AUTO-REMOVED #include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
// AUTO-REMOVED #include <core/id/DOF_ID.hh>
// AUTO-REMOVED #include <core/id/DOF_ID_Map.fwd.hh>
// AUTO-REMOVED #include <core/id/DOF_ID_Map.hh>
// AUTO-REMOVED #include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
// AUTO-REMOVED #include <core/id/types.hh>
// AUTO-REMOVED #include <core/import_pose/import_pose.hh>
// AUTO-REMOVED #include <core/io/pdb/file_data.fwd.hh>
// AUTO-REMOVED #include <core/io/silent/SharedSilentData.hh>
// AUTO-REMOVED #include <core/io/silent/SilentEnergy.hh>
// AUTO-REMOVED #include <core/io/silent/SilentFileData.fwd.hh>
// AUTO-REMOVED #include <core/io/silent/SilentStruct.fwd.hh>
// AUTO-REMOVED #include <core/io/silent/SilentStruct.hh>
// AUTO-REMOVED #include <core/io/silent/silent.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/AtomPointer.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/AtomTree.hh>
// AUTO-REMOVED #include <core/kinematics/AtomWithDOFChange.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/DomainMap.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/Edge.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Jump.hh>
// AUTO-REMOVED #include <core/kinematics/MinimizerMapBase.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/RT.hh>
// AUTO-REMOVED #include <core/kinematics/ResidueCoordinateChangeList.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Stub.hh>
// AUTO-REMOVED #include <core/kinematics/tree/Atom.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/tree/Atom.hh>
// AUTO-REMOVED #include <core/pose/MiniPose.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/annotated_sequence.hh>
// AUTO-REMOVED #include <core/pose/util.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
// AUTO-REMOVED #include <core/scoring/ScoreType.hh>
// AUTO-REMOVED #include <core/scoring/constraints/AmbiguousNMRConstraint.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/AmbiguousNMRDistanceConstraint.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/Func.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/Func.hh>
#include <core/sequence/DerivedSequenceMapping.fwd.hh>
#include <core/sequence/ScoringScheme.fwd.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <core/sequence/SequenceAlignment.fwd.hh>
// AUTO-REMOVED #include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/noesy_assign/CrossPeak.fwd.hh>
// AUTO-REMOVED #include <protocols/noesy_assign/CrossPeak.hh>
#include <protocols/noesy_assign/CrossPeakInfo.fwd.hh>
#include <protocols/noesy_assign/CrossPeakInfo.hh>
#include <protocols/noesy_assign/CrossPeakList.fwd.hh>
#include <protocols/noesy_assign/NoesyModule.fwd.hh>
#include <protocols/noesy_assign/PeakAssignment.fwd.hh>
#include <protocols/noesy_assign/PeakAssignmentResidueMap.fwd.hh>
#include <protocols/noesy_assign/PeakCalibrator.fwd.hh>
#include <protocols/noesy_assign/PeakFileFormat.fwd.hh>
#include <protocols/noesy_assign/Resonance.hh>
#include <protocols/noesy_assign/ResonanceList.fwd.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
// AUTO-REMOVED #include <utility/stream_util.hh>
// AUTO-REMOVED #include <utility/string_util.hh>
// AUTO-REMOVED #include <utility/vector0_bool.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
// AUTO-REMOVED #include <utility/file/gzip_util.hh>
// AUTO-REMOVED #include <utility/io/irstream.fwd.hh>
// AUTO-REMOVED #include <utility/io/irstream.hh>
// AUTO-REMOVED #include <utility/io/izstream.fwd.hh>
// AUTO-REMOVED #include <utility/io/mpistream.hh>
// AUTO-REMOVED #include <utility/io/orstream.fwd.hh>
// AUTO-REMOVED #include <utility/io/orstream.hh>
// AUTO-REMOVED #include <utility/io/ozstream.fwd.hh>
// AUTO-REMOVED #include <utility/io/zipstream.hpp>
// AUTO-REMOVED #include <utility/io/zipstream.ipp>
#include <utility/keys/AutoKey.fwd.hh>
#include <utility/keys/AutoKey.hh>
#include <utility/keys/Key.fwd.hh>
#include <utility/keys/Key.hh>
// AUTO-REMOVED #include <utility/keys/Key2Tuple.fwd.hh>
// AUTO-REMOVED #include <utility/keys/Key2Tuple.hh>
// AUTO-REMOVED #include <utility/keys/Key3Tuple.fwd.hh>
// AUTO-REMOVED #include <utility/keys/Key3Tuple.hh>
// AUTO-REMOVED #include <utility/keys/Key4Tuple.fwd.hh>
// AUTO-REMOVED #include <utility/keys/Key4Tuple.hh>
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
// AUTO-REMOVED #include <utility/signals/PausableSignalHub.fwd.hh>
// AUTO-REMOVED #include <utility/signals/PausableSignalHub.hh>
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
// AUTO-REMOVED #include <ObjexxFCL/Dimension.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/Dimension.hh>
// AUTO-REMOVED #include <ObjexxFCL/DimensionExpression.hh>
// AUTO-REMOVED #include <ObjexxFCL/DynamicIndexRange.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/DynamicIndexRange.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray1.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray1.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray1D.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray1D.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray2.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray2.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray2D.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray2D.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArrayInitializer.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArrayInitializer.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArraySection.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArraySection.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArrayTraits.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArrayTraits.hh>
// AUTO-REMOVED #include <ObjexxFCL/IndexRange.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/IndexRange.hh>
// AUTO-REMOVED #include <ObjexxFCL/InitializerSentinel.hh>
// AUTO-REMOVED #include <ObjexxFCL/Observer.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/Observer.hh>
// AUTO-REMOVED #include <ObjexxFCL/ObserverMulti.hh>
// AUTO-REMOVED #include <ObjexxFCL/ObserverSingle.hh>
// AUTO-REMOVED #include <ObjexxFCL/ProxySentinel.hh>
// AUTO-REMOVED #include <ObjexxFCL/SetWrapper.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/Star.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/Star.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/char.functions.hh>
// AUTO-REMOVED #include <ObjexxFCL/proxy_const_assert.hh>
#include <ObjexxFCL/string.functions.hh>
#include <algorithm>
#include <bitset>
#include <cassert>
#include <cmath>
#include <cstddef>
// AUTO-REMOVED #include <cstdio>
#include <cstdlib>
// AUTO-REMOVED #include <execinfo.h>
// AUTO-REMOVED #include <fstream>
#include <iomanip>
// AUTO-REMOVED #include <ios>
#include <iosfwd>
#include <iostream>
// AUTO-REMOVED #include <istream>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <time.h>
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
#include <basic/options/after_opts.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>
// AUTO-REMOVED #include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
// AUTO-REMOVED #include <boost/functional/hash.hpp>
// AUTO-REMOVED #include <zlib/zlib.h>
// AUTO-REMOVED #include <zlib/zutil.h>

//Auto Headers



static basic::Tracer tr("main");

using namespace core;
//using namespace protocols;

// OPT_KEY( FileVector, i )
// OPT_KEY( File, o )
// OPT_KEY( File, r )
// OPT_KEY( File, ro )
// OPT_KEY( File, simple_cst )
// //OPT_KEY( File, s )
OPT_2GRP_KEY( File, noesy, out, cst )
// OPT_KEY( String, assign )
// OPT_KEY( Boolean, ignore_assignments )
// OPT_KEY( Boolean, force_symmetry )
// OPT_KEY( Boolean, sequential )

// OPT_KEY( File, frags_for_ss )
// OPT_KEY( Real, hop_penalty )
// OPT_KEY( Integer, max_hops )
// OPT_KEY( Integer, prune_cutoff )
// OPT_KEY( Real, start_weight )
// OPT_KEY( Integer, edge_redundancy )
// OPT_KEY( Boolean, no_decoys )
// OPT_KEY( Boolean, no_network )
// OPT_KEY( Boolean, no_symm )
// OPT_KEY( Boolean, no_cs )
// OPT_KEY( Boolean, no_upper )
// OPT_KEY( Boolean, no_remove_diagonal )
// OPT_KEY( Boolean, no_calibrate )


// Had these "iterative" options in here to enable testing of NOE in script with same flags as archive run...

//OPT_1GRP_KEY( Boolean, iterative, assign_noes )
OPT_2GRP_KEY( Integer, noesy, out, min_seq_sep )
//OPT_1GRP_KEY( Real, iterative, centroid_before_quickrelax_weight )
//OPT_1GRP_KEY( Real, iterative, fullatom_after_quickrelax_weight )

void register_options() {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
	protocols::noesy_assign::NoesyModule::register_options();
  //  Templates::register_options();
  // NEW_OPT( i, "read NOEs from this file", "NOE_in.dat" );
//   NEW_OPT( o, "write NOEs to this file", "NOE_out.dat" );
//   NEW_OPT( r, "read resonance assignments from file", "cs_in.dat" );
//   NEW_OPT( ro, "write resonance assignments to file", "cs_out.dat" );
//   NEW_OPT( simple_cst, "make simple constraints from assigned peaks", "NOE_simple.cst" );
//	OPT( in::file::s );
// 	OPT( in::file::silent );
	OPT( in::file::fasta );
// 	//NEW_OPT( s, "structure", "native.pdb" );
 	NEW_OPT( noesy::out::cst, "write constraints", "assigned.cst" );
	NEW_OPT( noesy::out::min_seq_sep, "do not write constraints from peaks that have any residue pair less than sep. assigned", 2 );
// 	NEW_OPT( assign, "use an assignment strategy [simple, cyana]", "cyana" );
// 	NEW_OPT( ignore_assignments, "don't read existing assignments from peak-file", true );
// 	NEW_OPT( force_symmetry, "assign only symmetric peaks (makes sense, really)", false );
// 	NEW_OPT( sequential, "always assign cross-peak as sequential NOE if possible", false );
// 	NEW_OPT( frags_for_ss, "fragments that define SS structure for NOE Network analysis", "frags3.dat" );
// 	NEW_OPT( hop_penalty, "extra resistance for each hop", 2.0 );
// 	NEW_OPT( max_hops, "how many hops do we evaluate at most", 5);
// 	NEW_OPT( prune_cutoff, "if expected R_inv is at least X times larger, then total R, don't evaluated path", 10 );
// 	NEW_OPT( start_weight, "start with this weight", 0.1 );
// 	NEW_OPT( edge_redundancy, "need for direct NOEs to keep weight ", 4 );
// 	NEW_OPT( no_decoys, "check comp. with decoys", false );
// 	NEW_OPT( no_network, "check comp. with decoys", false );
// 	NEW_OPT( no_symm, "check comp. with decoys", false );
// 	NEW_OPT( no_cs, "check comp. with decoys", false );
// 	NEW_OPT( no_upper, "check upper", false );
// 	NEW_OPT( no_remove_diagonal, "", false );
// 	NEW_OPT( no_calibrate, "don't calibrate the NOE distance bound", false );
//	NEW_OPT( iterative::fullatom_after_quickrelax_weight, "[IGNORED] just to make cyana_test happy", 1.0);
//	NEW_OPT( iterative::centroid_before_quickrelax_weight, "[IGNORED] just to make cyana_test happy", 1.0);
}

void run_old() {
	using namespace protocols::noesy_assign;


// 	std::string fasta_sequence;
// 	pose::PoseOP native_pose = NULL;
// 	if ( basic::options::option[ basic::options::OptionKeys::in::file::fasta ].user() ) {
// 		fasta_sequence = core::sequence::read_fasta_file( basic::options::option[ basic::options::OptionKeys::in::file::fasta ]()[1] )[1]->sequence();
// 		tr.Info << "read fasta sequence: " << fasta_sequence.size() << " residues\n"  << fasta_sequence << std::endl;
// 	} else if ( native_pose ) {
// 		fasta_sequence = native_pose->sequence();
// 		tr.Info << "take sequence from native : " << fasta_sequence << std::endl;
// 	}
// 	ResonanceListOP resonances( new ResonanceList( fasta_sequence ) );
//   { // read Resonances
//     utility::io::izstream input_file(basic::options::option[ basic::options::OptionKeys::r ]() );
//     utility::io::ozstream output_file( basic::options::option[ basic::options::OptionKeys::ro ]() );
//     if ( input_file.good() ) {
//       resonances->read_from_stream( input_file );
//       resonances->write_to_stream( output_file );
//     } else {
//       tr.Error << "cannot read " << basic::options::option[ basic::options::OptionKeys::r ]() << std::endl;
//     }
//   }

// 	CrossPeakList cpl( resonances );
//   { // read CrossPeaks
// 		Size nfiles( basic::options::option[ basic::options::OptionKeys::i ]().size() );
// 		for ( core::Size ifile = 1; ifile <= nfiles; ++ifile ) {
// 			std::string file( basic::options::option[ basic::options::OptionKeys::i ]()[ ifile ] );
// 			utility::io::izstream input_file( file );
// 			if ( input_file.good() ) {
// 				PeakFileFormat_xeasy format;
// 				format.set_ignore_assignments( basic::options::option[ basic::options::OptionKeys::ignore_assignments ]() );
// 				cpl.read_from_stream( input_file, format );
// 			} else {
// 				tr.Error << "cannot read " << file << std::endl;
// 			}
// 		}
//   }



// 	if ( basic::options::option[ basic::options::OptionKeys::assign ].user() ) {
// 		std::string assign_strategy( basic::options::option[ basic::options::OptionKeys::assign ]() );
// 		if ( assign_strategy == "simple" ) {
// 			cpl.find_assignments();
// 		} else if ( assign_strategy == "cyana" ) {
// 			using namespace basic::options;
// 			using namespace basic::options::OptionKeys;
// 			cpl.find_assignments();
// 			if ( !option[ no_remove_diagonal ]() ) cpl.delete_diagonal_peaks();
// 			if ( !option[ no_cs ]() ) cpl.update_chemshiftscore();
// 			if ( !option[ no_symm ]() ) cpl.update_symmetry_score();
// 			if ( !option[ no_upper ]() ) cpl.update_upperdistance_score();
// 			core::io::silent::SilentFileData sfd;
// 			if ( !option[ no_decoys ]() && option[ in::file::silent ].user() ) {
// 				sfd.read_file( basic::options::option[ basic::options::OptionKeys::in::file::silent ]()[ 1 ] );
// 				cpl.update_decoy_compatibility_score( sfd.begin(), sfd.end() );
// 			} else {
// 				cpl.set_trivial_decoy_compatibility_score();
// 			}
// 			if ( !option[ no_network ]() ) cpl.network_analysis();
// 			cpl.update_peak_volumina();
// 			if ( !option[  no_calibrate ] () ) cpl.calibrate( sfd.begin(), sfd.end() );
// 			cpl.eliminate_spurious_peaks();
// 		} else {
// 			tr.Error << "assing strategy " << assign_strategy << " not recognized " << std::endl;
// 		}
// 	}

// 	if ( basic::options::option[ basic::options::OptionKeys::force_symmetry ]() ) {
// 		PeakAssignmentList assignments( resonances );
// 		assignments.add( cpl );
// 		assignments.check_for_symmetric_peaks( cpl );
// 	}

// 	if ( basic::options::option[ basic::options::OptionKeys::sequential ]() ) {
// 		PeakAssignmentList assignments( resonances );
// 		assignments.add( cpl );
// 		assignments.invalidate_competitors_to_sequential_NOE( cpl );
// 	}

// 	core::pose::Pose pose;
// 	core::import_pose::pose_from_pdb( pose, basic::core::options::option[ options::OptionKeys::in::file::s ]()[ 1 ] );

// 	core::scoring::constraints::ConstraintSetOP cstset = cpl.generate_constraints( pose );
// 	core::scoring::constraints::ConstraintIO::write_constraints(  basic::options::option[ basic::options::OptionKeys::cst_out ](), *cstset, pose );

// 	core::scoring::constraints::ConstraintSetOP centroid_cstset = cpl.generate_constraints( pose, true );
// 	core::pose::Pose centroid_pose = pose;
// 	core::util::switch_to_residue_type_set( centroid_pose, core::chemical::CENTROID );
// 	core::scoring::constraints::ConstraintIO::write_constraints(  basic::options::option[ basic::options::OptionKeys::cst_out ]().name()
// 		+ ".centroid", *centroid_cstset, centroid_pose );



// 	utility::io::ozstream output_file( basic::options::option[ basic::options::OptionKeys::o ]() );
// 	PeakFileFormat_xeasy format;
// 	cpl.write_to_stream( output_file, format );

// 	fragment::FragSetOP ss_frags = fragment::FragmentIO().read_data( basic::options::option[ basic::options::OptionKeys::frags_for_ss ]() );
// 	core::conformation::SecondaryStructure ss_def( *ss_frags, false /*no JustUseCentralResidue */ );

// 	utility::vector1< bool > beta_ss( ss_def.total_residue(), false );
// 	utility::vector1< bool > helix_ss( ss_def.total_residue(), false );
// 	utility::vector1< bool > any_ss( ss_def.total_residue(), false );
// 	for ( Size i=1; i<=ss_def.total_residue(); ++i ) {
// 		if ( ss_def.strand_fraction( i )>0.7 ) beta_ss[ i ]=true;
// 		if ( ss_def.helix_fraction( i )>0.7 ) helix_ss[ i ] = true;
// 		if ( helix_ss[ i ] ) tr.Debug << "helix " << i << std::endl;
// 		if ( beta_ss[ i ] ) tr.Debug << "beta: " << i << std::endl;
// 		if ( helix_ss[ i ] || beta_ss[ i ] ) any_ss[ i ]=true;
// 	}

// 	NoeNetwork network(
// 																	resonances,
// 																	basic::options::option[ basic::options::OptionKeys::max_hops ](),
// 																	basic::options::option[ basic::options::OptionKeys::hop_penalty ](),
// 																	basic::options::option[ basic::options::OptionKeys::prune_cutoff ]()
// 	);
// 	network.set_starting_weight( basic::options::option[ basic::options::OptionKeys::start_weight ]() );
// 	network.set_edge_redundancy( basic::options::option[ basic::options::OptionKeys::edge_redundancy ]() );
// 	network.add( cpl );
// 	network.add( ss_def );
// 	Real delta( 1.0 );
// 	while( delta > 10 ) {
// 		tr.Debug << network << std::endl;
// 		delta = network.evaluate_edges();
// 	}
// 	tr.Debug << network << std::endl;
	//pose.dump_pdb("the_reference_pose.pdb" );
	/*
	utility::io::ozstream cst_output_file( basic::options::option[ basic::options::OptionKeys::cst_out ]() );
	//  core::scoring::constraints::ConstraintSetOP cstset( new core::scoring::constraints::ConstraintSet );
  {//make simple constraint list from assigned peaks
    using namespace core::scoring::constraints;
    Size ct( 1 );
    for ( CrossPeakList::CrossPeaks::const_iterator it = cpl.peaks().begin(); it != cpl.peaks().end(); ++it, ++ct ) {
			tr.Debug << " find assignments for peak " << ct << std::endl;
			//			(*it)->find_assignments( resonances );
      if ( (*it)->assigned() && !(*it)->ambiguous() ) {
				PeakAssignment const& assignment( **((*it)->assignments().begin()));
				id::NamedAtomID const& atom1( assignment.atom( *resonances, 1 ) ); //[ (*it)->proton( 1 ).assignment( ind_assigned ) ].atom() );
				id::NamedAtomID const& atom2( assignment.atom( *resonances, 2 ) ); //[ (*it)->proton( 2 ).assignment( ind_assigned ) ].atom() );
				if ( atom1.rsd() != atom2.rsd() ) {
			// 		if ( beta_ss[ atom1.rsd() ] && beta_ss[ atom2.rsd() ]
// 						&& atom1.rsd() > 1 && atom1.rsd() < ss_def.total_residue() && ( beta_ss[ atom1.rsd()-1 ] || beta_ss[ atom1.rsd()+1] )
// 						&& atom2.rsd() > 1 && atom2.rsd() < ss_def.total_residue() && ( beta_ss[ atom2.rsd()-1 ] || beta_ss[ atom2.rsd()+1] )
// 					) {
					cst_output_file << "AmbiguousNMRDistance " << atom1 << " " << atom2 << " BOUNDED 1.5 6 0.5 NOE ; unambiguous NOE " << ct << std::endl;					//}
				}
				tr.Debug << " write assignments as constraints for peak " << ct << std::endl;
				//	cstset->add_constraint( new AmbiguousNMRDistanceConstraint( atom1, atom2, pose, new BoundFunc( 1.5, 5.5, 0.5, "NOE" ) ) );
      } else if ( (*it)->assigned() ) {
				cst_output_file << "AmbiguousConstraint" << std::endl;
				//				Size ind_assigned( 0 );
				for ( CrossPeak::PeakAssignments::const_iterator ait=(*it)->assignments().begin(); ait!=(*it)->assignments().end(); ++ait ) {
					PeakAssignment const& assignment( **ait );
					id::NamedAtomID const& atom1( assignment.atom( *resonances, 1 ) ); //[ (*it)->proton( 1 ).assignment( ind_assigned ) ].atom() );
					id::NamedAtomID const& atom2( assignment.atom( *resonances, 2 ) ); //[ (*it)->proton( 2 ).assignment( ind_assigned ) ].atom() );
					if ( atom1.rsd() != atom2.rsd() ) {
	// 				if ( beta_ss[ atom1.rsd() ] && beta_ss[ atom2.rsd() ]
// 						&& atom1.rsd() > 1 && atom1.rsd() < ss_def.total_residue() && ( beta_ss[ atom1.rsd()-1 ] || beta_ss[ atom1.rsd()+1] )
// 						&& atom2.rsd() > 1 && atom2.rsd() < ss_def.total_residue() && ( beta_ss[ atom2.rsd()-1 ] || beta_ss[ atom2.rsd()+1] )
						//	) {
						cst_output_file << "AmbiguousNMRDistance " << atom1 << " " << atom2 << " BOUNDED 1.5 6 0.5 NOE ; ambiguous NOE " << ct << std::endl;
					}
				}
				cst_output_file << "END_Ambiguous" << std::endl;
      }
    }
  }
	*/
}

void run() {
	using namespace protocols::noesy_assign;

	std::string fasta_sequence;
	if ( basic::options::option[ basic::options::OptionKeys::in::file::fasta ].user() ) {
		fasta_sequence = core::sequence::read_fasta_file( basic::options::option[ basic::options::OptionKeys::in::file::fasta ]()[1] )[1]->sequence();
		tr.Info << "read fasta sequence: " << fasta_sequence.size() << " residues\n"  << fasta_sequence << std::endl;
	}

	PROF_START( basic::NOESY_ASSIGN_TOTAL );

	NoesyModule nm( fasta_sequence );

	nm.assign();

 	core::pose::Pose pose;
// 	core::import_pose::pose_from_pdb( pose, basic::options::option[ options::OptionKeys::in::file::s ]()[ 1 ] );

	core::pose::make_pose_from_sequence(
		pose,
		fasta_sequence,
		chemical::FA_STANDARD
	);

	std::string cst_file( basic::options::option[ basic::options::OptionKeys::noesy::out::cst ]() );
	std::string cst_centroid_file( cst_file + ".centroid");

	nm.generate_constraint_files( pose, cst_file, cst_centroid_file,
		basic::options::option[ basic::options::OptionKeys::noesy::out::min_seq_sep ]() );

	nm.write_assignments();

	PROF_STOP( basic::NOESY_ASSIGN_TOTAL );
	basic::prof_show();
}

int
main( int argc, char * argv [] )
{
	register_options();
	protocols::noesy_assign::PeakAssignmentParameters::register_options();
	protocols::noesy_assign::PeakFileFormat_xeasy::register_options();
	devel::init( argc, argv );
	protocols::noesy_assign::PeakAssignmentParameters::get_instance();
	try{
		run();
	} catch ( utility::excn::EXCN_Base& anExcn ) {
		anExcn.show( std::cerr );
	}

	return 0;
}
