// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file score_aln.cc
/// @brief
/// @author James Thompson

#include <core/types.hh>
#include <devel/init.hh>

#include <basic/Tracer.hh>

#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/scoring/rms_util.hh>

#include <basic/options/option.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/ScoringSchemeFactory.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/Residue.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>

#include <protocols/comparative_modeling/ThreadingMover.hh>

#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <numeric/model_quality/rms.hh>

// C++ headers
#include <map>
#include <iostream>
#include <string>

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

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
#include <core/fragment/FragData.hh>
#include <core/fragment/FragID.fwd.hh>
#include <core/fragment/FragID_Iterator.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.fwd.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.fwd.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FrameIteratorWorker_.fwd.hh>
#include <core/fragment/FrameIteratorWorker_.hh>
#include <core/fragment/FrameList.fwd.hh>
#include <core/fragment/FrameList.hh>
#include <core/fragment/SingleResidueFragData.fwd.hh>
#include <core/fragment/SingleResidueFragData.hh>
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
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/types.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/file_data.fwd.hh>
#include <core/io/silent/SharedSilentData.hh>
#include <core/io/silent/SilentEnergy.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
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
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/sequence/DerivedSequenceMapping.fwd.hh>
#include <core/sequence/ScoringScheme.fwd.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <core/sequence/SequenceAlignment.fwd.hh>
#include <core/sequence/SequenceProfile.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/jobdist/Jobs.fwd.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/down_cast.hh>
#include <utility/stream_util.hh>
#include <utility/vector0_bool.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
#include <utility/file/gzip_util.hh>
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
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/proxy_const_assert.hh>
#include <ObjexxFCL/ubyte.fwd.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <execinfo.h>
#include <fstream>
#include <functional>
#include <iomanip>
#include <ios>
#include <iosfwd>
#include <istream>
#include <limits>
#include <list>
#include <ostream>
#include <set>
#include <sstream>
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
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>
#include <boost/unordered/unordered_map.hpp>
#include <zlib/zlib.h>
#include <zlib/zutil.h>



///////////////////////////////////////////////////////////////////////////////

using std::map;
using std::string;
using core::Size;
using core::Real;
using utility::vector1;
using utility::file::FileName;
using ObjexxFCL::fmt::A;
using ObjexxFCL::fmt::F;
using ObjexxFCL::fmt::I;
using ObjexxFCL::string_of;
using namespace basic;
using namespace core::sequence;
using namespace basic::options;
using namespace basic::options::OptionKeys;

class FastThreadingMover : public protocols::comparative_modeling::ThreadingMover {

public:
	FastThreadingMover(
		core::sequence::SequenceAlignment aln,
		core::pose::Pose templ
	) : ThreadingMover( aln, templ ) {}

	void apply( core::pose::Pose & query_pose ) {
		// don't rebuild loops
		build_loops( false );
		repack_query( false );
		ThreadingMover::apply( query_pose );
	}

private:
	core::pose::Pose native_pose;
};

std::map< std::string, core::pose::Pose >
poses_from_cmd_line(
	utility::vector1< std::string > const & fn_list
) {
	using std::map;
	using std::string;
	using core::pose::Pose;
	using utility::file::file_exists;
	using core::import_pose::pose_from_pdb;
	using namespace core::chemical;

	ResidueTypeSetCAP rsd_set = rsd_set_from_cmd_line();
	map< string, Pose > poses;

	typedef vector1< string >::const_iterator iter;
	for ( iter it = fn_list.begin(), end = fn_list.end(); it != end; ++it ) {
		if ( file_exists(*it) ) {
			Pose pose;
			core::import_pose::pose_from_pdb( pose, *rsd_set, *it );
			string name = utility::file_basename( *it );
			name = name.substr( 0, 5 );
			poses[name] = pose;
		}
	}

	return poses;
}

void print_seq_map(
	std::ostream & out,
	std::map< std::string, core::sequence::SequenceOP > const & seqs
) {
	using std::map;
	using std::string;
	using core::sequence::SequenceOP;
	typedef map< string, SequenceOP >::const_iterator iter;
	for ( iter it = seqs.begin(), end = seqs.end(); it != end; ++it ) {
		out << it->first << " => " << *it->second << std::endl;
	}
}

std::map< std::string, core::sequence::SequenceOP >
sequences_from_cmd_line(
	utility::vector1< std::string > const & fn_list
) {
	using std::map;
	using std::string;
	using utility::file::file_exists;
	using namespace core::sequence;

	map< string, SequenceOP > seqs;

	typedef vector1< string >::const_iterator iter;
	for ( iter it = fn_list.begin(), end = fn_list.end(); it != end; ++it ) {
		SequenceProfileOP prof( new SequenceProfile );
		if ( file_exists(*it) ) {
			prof->read_from_file( *it, 1.0 );
			string name = utility::file_basename( *it );
			name = name.substr( 0, 5 );
			seqs[name] = prof;
		}
	}

	return seqs;
}

void save_per_residue_scores(
	std::string const & fn,
	core::sequence::SequenceAlignment const & aln,
	core::sequence::ScoringSchemeOP ss,
	std::string const & aln_id
) {
	using core::Real;
	using utility::vector1;

	vector1< Real > scores( aln.calculate_per_position_scores(ss) );
	utility::io::ozstream output(fn);
	output << "resi score aln_id" << std::endl;
	for ( Size ii = 1; ii <= aln.length(); ++ii ) {
		if ( !aln.sequence(1)->is_gap(ii) ) {
			Size const resi ( aln.sequence(1)->resnum(ii) );
			Real const score( scores[ii] );
			output << resi << ' ' << score << ' ' << aln_id << std::endl;
		}
	}
	output.close();
}

int
main( int argc, char* argv [] ) {

	// options, random initialization
	devel::init( argc, argv );

	using std::map;
	using std::string;
	using core::Real;
	using core::Size;
	using core::pose::Pose;
	using utility::vector1;
	using core::sequence::SequenceAlignment;
	using core::sequence::SequenceProfile;
	using core::import_pose::pose_from_pdb;
	using namespace core::chemical;
	using namespace core::io::silent;

	basic::Tracer tr( "score_aln" );

	SequenceProfileOP query_prof( new SequenceProfile );
	map< string, SequenceOP > seqs;
	if ( option[ in::file::pssm ].user() ) {
		query_prof->read_from_file( option[ in::file::pssm ]()[1], 1.0 );
		seqs = sequences_from_cmd_line(
			option[ in::file::pssm ]()
		);
	}

	vector1< std::string > align_fns = option[ in::file::alignment ]();

	Pose native_pose;
	bool have_native( false );
	if ( option[ in::file::native ].user() ) {
		core::import_pose::pose_from_pdb(
			native_pose,
			*(rsd_set_from_cmd_line()),
			option[ in::file::native ]()
		);
		have_native = true;
	}

	map< string, Pose > poses;
	if ( have_native ) {
		poses = poses_from_cmd_line(
			option[ in::file::template_pdb ]()
		);
	}

	// scoring scheme for aligning profiles
	std::string const scoring_scheme_type( option[ cm::seq_score ]()[1] );
	ScoringSchemeFactory ssf;
	ScoringSchemeOP ss( ssf.get_scoring_scheme( scoring_scheme_type ) );

	// for the ProfSim scoring scheme, the optimal opening and extension
	// penalties were 2 and 0.2, with a scoring shift of -0.45 applied to
	// all ungapped aligned pairs for database searches.
	Real const gap_open  ( option[ cm::min_gap_open ]() );
	Real const gap_extend( option[ cm::min_gap_extend ]() );
	ss->gap_open  ( gap_open   );
	ss->gap_extend( gap_extend );

	SilentFileData sfd;

	typedef vector1< string >::const_iterator aln_iter;
	for ( aln_iter aln_fn = align_fns.begin(), aln_end = align_fns.end();
				aln_fn != aln_end; ++aln_fn
	) {
		vector1< SequenceAlignment > alns = core::sequence::read_aln(
			option[ cm::aln_format ](), *aln_fn
		);

		for ( vector1< SequenceAlignment >::iterator it = alns.begin(),
				end = alns.end();
				it != end; ++it
		) {
			string const template_id( it->sequence(2)->id().substr(0,5) );
			tr << *it << std::endl;
			tr << "id " << it->sequence(2)->id() << " => " << template_id
				<< std::endl;
			string const ungapped_query( it->sequence(1)->ungapped_sequence() );

			SilentStructOP ss_out( new ScoreFileSilentStruct );
			//bool encountered_error(false);
			map< string, SequenceOP >::iterator seq_it = seqs.find( template_id );
			if ( seq_it == seqs.end() ) {
				//print_seq_map( std::cerr, seqs );
				string msg( "Error: can't find seq (id = " + template_id + ")" );
				//utility_exit_with_message(msg);
				tr.Error << msg << std::endl;
				//continue;
			} else {
				SequenceOP template_prof = seq_it->second;

				utility::vector1< SequenceOP > my_seqs;
				my_seqs.push_back( query_prof->clone() );
				my_seqs.push_back( template_prof->clone() );
				SequenceAlignment rescore_aln = steal_alignment( *it, my_seqs );
				Real const aln_score( rescore_aln.calculate_score_sum_of_pairs( ss ) );
				rescore_aln.score( aln_score );
				//string const ungapped_templ( it->sequence(2)->ungapped_sequence() );
				save_per_residue_scores(
					it->sequence(2)->id() + ".dat", rescore_aln, ss,
					it->sequence(2)->id()
				);

				ss_out->add_energy( "aln_score", aln_score );
				//ss_out->add_energy( "n_ali_templ", ungapped_templ.length() );
				tr.Debug << "score(" << ss_out->decoy_tag() << ") = " << aln_score
					<< std::endl;
			}

			if ( have_native ) {
				// calc rmsd/gdt stats

				Pose query_pose, template_pose;
				core::pose::make_pose_from_sequence(
					query_pose,
					ungapped_query,
					*(rsd_set_from_cmd_line())
				);

				map< string, Pose >::iterator pose_it = poses.find( template_id );
				if ( pose_it == poses.end() ) {
					string msg( "Error: can't find pose (id = "
						+ template_id + ")"
					);
					//utility_exit_with_message(msg);
					tr.Error << msg << std::endl;
					continue;
				} else {
					template_pose = pose_it->second;

					static string const atm( "CA" );
					core::id::SequenceMapping mapping = it->sequence_mapping(1,2);
					using utility::vector1;
					using numeric::xyzVector;
					using numeric::model_quality::calc_rms;
					using numeric::model_quality::rms_wrapper;
					using core::scoring::xyz_gdtmm;
					vector1< xyzVector< Real > > native_coords, template_coords;

					int natoms(0);
					for ( Size ii = 1; ii <= native_pose.total_residue(); ++ii ) {
						Size const templ_ii( mapping[ii] );
						bool skip(
							templ_ii == 0 ||
							ii > native_pose.total_residue() ||
							templ_ii > template_pose.total_residue()
						);
						if ( !skip ) ++natoms;
					}

					ObjexxFCL::FArray2D< core::Real > p1a( 3, natoms );
					ObjexxFCL::FArray2D< core::Real > p2a( 3, natoms );
					Size n_gap(0);
					for ( Size ii = 1; ii <= native_pose.total_residue(); ++ii ) {
						Size const templ_ii( mapping[ii] );
						bool skip(
							templ_ii == 0 ||
							ii > native_pose.total_residue() ||
							templ_ii > template_pose.total_residue()
						);
						if ( skip ) {
							n_gap++;
							continue;
						}

						core::Vector native_xyz  ( native_pose.residue(ii).xyz(atm) );
						core::Vector template_xyz( template_pose.residue(templ_ii).xyz(atm) );
						native_coords.push_back(native_xyz);
						template_coords.push_back(template_xyz);
						for ( Size jj = 1; jj <= 3; ++jj ) {
							p1a(jj,ii - n_gap) = native_xyz  [jj-1];
							p2a(jj,ii - n_gap) = template_xyz[jj-1];
						}
					}
					runtime_assert( native_coords.size() == template_coords.size() );

					//Real const debug_rmsd_ali( calc_rms( native_coords, template_coords ) );
					Real const rmsd_ali ( rms_wrapper( natoms, p1a, p2a ) );
					Real const gdtmm_ali( xyz_gdtmm( p1a, p2a ) );
					Real const coverage(
						(Real) native_coords.size() / (Real) native_pose.total_residue()
					);
					Real const gdtmm_overall( gdtmm_ali * coverage );

					//Real const gdtmm_ali( 1.0 );
					ss_out->add_energy( "coverage", coverage );
					ss_out->add_energy( "rmsd_ali", rmsd_ali );
					ss_out->add_energy( "gdtmm_ali", gdtmm_ali );
					ss_out->add_energy( "gdtmm_overall", gdtmm_overall );
					ss_out->add_energy( "nres_query", native_pose.total_residue() );
					//ss_out->add_energy( "debug_rmsd_ali", debug_rmsd_ali );
				} // template pdb check
			} // have native
			ss_out->scoreline_prefix( "" );
			ss_out->decoy_tag( it->sequence(2)->id() );
			ss_out->add_energy( "n_ali_query", ungapped_query.length() );
			ss_out->add_string_value( "template", template_id );
			sfd.write_silent_struct( *ss_out, option[ out::file::silent ]() );
		} // alns
	} // for ( it in aligns )

	tr.Debug << "finished rescoring alignments." << std::endl;
	tr.flush_all_channels();
} // int main( int argc, char * argv [] )
