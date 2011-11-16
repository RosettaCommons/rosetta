// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file sequence_recovery.cc
/// @brief A protocol which outputs sequence recovery statistics ala the table in the "Native sequences are close to optimal" paper.
/// @author Ron Jacak
/// @author P. Douglas Renfrew (renfrew@unc.edu) ( added rotamer recovery, cleanup )

// Unit headers
#include <devel/init.hh>

//project Headers
#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>
#include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <basic/options/util.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/TaskOperationFactory.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>

#include <protocols/moves/PackRotamersMover.hh>

// Utility Headers
#include <basic/Tracer.hh>
// AUTO-REMOVED #include <basic/MetricValue.hh>
#include <basic/prof.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/ozstream.hh>
#include <utility/vector1.hh>

// Option keys
// AUTO-REMOVED #include <basic/options/keys/optE.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

// Numeric Headers

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray2D.hh>

// C++ headers
#include <sstream>

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
#include <core/conformation/PointGraph.fwd.hh>
#include <core/conformation/PointGraphData.fwd.hh>
#include <core/conformation/PointGraphData.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>
#include <core/conformation/find_neighbors.fwd.hh>
#include <core/conformation/orbitals/OrbitalXYZCoords.hh>
#include <core/conformation/signals/ConnectionEvent.fwd.hh>
#include <core/conformation/signals/IdentityEvent.fwd.hh>
#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/graph/ArrayPool.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/graph/Graph.hh>
#include <core/graph/UpperEdgeGraph.fwd.hh>
#include <core/graph/UpperEdgeGraph.hh>
#include <core/graph/unordered_object_pool.fwd.hpp>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/file_data.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.fwd.hh>
#include <core/pack/rotamer_set/RotamerCouplings.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/pack/task/IGEdgeReweightContainer.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/ResfileReader.fwd.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/operation/ResFilterCreator.fwd.hh>
#include <core/pack/task/operation/ResLvlTaskOperationCreator.fwd.hh>
#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperationCreator.fwd.hh>
#include <core/pack/task/operation/TaskOperationFactory.fwd.hh>
#include <core/pack/task/operation/TaskOperations.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/PDBPoseMap.fwd.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Remarks.fwd.hh>
#include <core/pose/Remarks.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
// AUTO-REMOVED #include <core/pose/metrics/PoseMetricCalculatorBase.fwd.hh>
// AUTO-REMOVED #include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/ContextGraph.fwd.hh>
#include <core/scoring/ContextGraph.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/EnergyGraph.fwd.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/LREnergyContainer.fwd.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/scoring/NeighborList.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
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
#include <core/scoring/symmetry/SymmetricScoreFunction.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/jobdist/Jobs.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/moves/PackRotamersMover.fwd.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/PyAssert.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/stream_util.hh>
#include <utility/string_util.hh>
#include <utility/vector0.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector0_bool.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
#include <utility/file/gzip_util.hh>
#include <utility/io/mpistream.hh>
#include <utility/io/orstream.fwd.hh>
#include <utility/io/orstream.hh>
#include <utility/io/ozstream.fwd.hh>
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
#include <utility/signals/SignalHub.fwd.hh>
#include <utility/signals/SignalHub.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/Tag.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzTriple.fwd.hh>
#include <numeric/xyzTriple.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
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
#include <ObjexxFCL/string.functions.hh>
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
#include <iomanip>
#include <ios>
#include <iosfwd>
#include <iostream>
#include <istream>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <set>
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
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/pool/poolfwd.hpp>
#include <boost/unordered_map.hpp>
#include <zlib/zlib.h>
#include <zlib/zutil.h>

//Auto using namespaces
namespace std { } using namespace std; // AUTO USING NS
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
//Auto using namespaces end



static basic::Tracer TR("sequence_recovery");

using namespace core;
using namespace protocols;
using namespace basic::options;
using namespace basic::options::OptionKeys;

using namespace ObjexxFCL::fmt;

namespace sequence_recovery {
	FileOptionKey const native_pdb_list( "sequence_recovery::native_pdb_list" );
	FileOptionKey const redesign_pdb_list( "sequence_recovery::redesign_pdb_list" );
	FileOptionKey const parse_taskops_file( "sequence_recovery::parse_taskops_file" );
	BooleanOptionKey const rotamer_recovery( "sequence_recovery::rotamer_recovery" );
	StringOptionKey const seq_recov_filename( "sequence_recovery::seq_recov_filename" );
	StringOptionKey const sub_matrix_filename( "sequence_recovery::sub_matrix_filename" );
	IntegerOptionKey const se_cutoff( "sequence_recovery::se_cutoff" );
}

std::string usage_string;

void init_usage_prompt( std::string exe ) {

	// place the prompt up here so that it gets updated easily; global this way, but that's ok
	std::stringstream usage_stream;
	usage_stream << "No files given: Use either -file:s or -file:l to designate a single pdb or a list of pdbs.\n\n"
					<< "Usage: " << exe
					<< "\n\t-database path/to/minidb"
					<< "\n\t-native_pdb_list <list file>"
					<< "\n\t-redesign_pdb_list <list file>"

					<< "\n\t[-seq_recov_filename <file>] file to output sequence recoveries to (default: sequencerecovery.txt)"
					<< "\n\t[-sub_matrix_filename <file>] file to output substitution matrix to (default: submatrix.txt)"
					<< "\n\t[-parse_taskops_file <file>] tagfile which contains task operations to apply before measuring recovery (optional)"
					<< "\n\t[-ignore_unrecognized_res]"

					<< "\n\n";

	usage_string = usage_stream.str();

}


///@brief load custom TaskOperations according to an xml-like utility::tag file
core::pack::task::TaskFactoryOP setup_tf( core::pack::task::TaskFactoryOP task_factory_ ) {

	using namespace core::pack::task::operation;

	if ( option[ sequence_recovery::parse_taskops_file ].user() ) {
		std::string tagfile_name( option[ sequence_recovery::parse_taskops_file ]() );
		TaskOperationFactory::TaskOperationOPs tops;
		TaskOperationFactory::get_instance()->newTaskOperations( tops, tagfile_name );
		for ( TaskOperationFactory::TaskOperationOPs::iterator it( tops.begin() ), itend( tops.end() ); it != itend; ++it ) {
			task_factory_->push_back( *it );
		}
	} else {
		task_factory_->push_back( new pack::task::operation::InitializeFromCommandline );
	}

	return task_factory_;

}


///@brief helper method which uses the tenA nb graph in the pose object to fill a vector with nb counts
void fill_num_neighbors( pose::Pose & pose, utility::vector1< core::Size > & num_nbs ) {

	using core::conformation::PointGraph;
	using core::conformation::PointGraphOP;

	PointGraphOP pg( new PointGraph ); // create graph
	core::conformation::residue_point_graph_from_conformation( pose.conformation(), *pg ); // create vertices
	core::conformation::find_neighbors<core::conformation::PointGraphVertexData,core::conformation::PointGraphEdgeData>( pg, 10.0 /* ten angstrom distance */ ); // create edges

	num_nbs.resize( pose.n_residue(), 0 );
	for ( core::Size ii=1; ii <= pose.total_residue(); ++ii ) {

		// a PointGraph is a typedef of UpperEdgeGraph< PointGraphVertexData, PointGraphEdgeData >
		// so any of the method in UpperEdgeGraph should be avail. here. The UpperEdgeGraph provides access to nodes
		// via a get_vertex() method, and each vertex can report back how many nbs it has.
		// So something that before was really complicated (nb count calculation) is done in <10 lines of code.
		// the assumption we're making here is that a pose residue position ii is the same index as the point graph vertex
		// that is indeed the case if you look at what the function residue_point_graph_from_pose().
		num_nbs[ ii ] = pg->get_vertex(ii).num_neighbors_counting_self();
	}

	return;
}

///@brief return the set of residues that are designable based given pose
std::set< Size > fill_designable_set( pose::Pose & pose, pack::task::TaskFactoryOP & tf ) {

	//we need to score the pose for many of the task operations passed from cmd line
	if ( option[ sequence_recovery::parse_taskops_file ].user() ) {
		scoring::ScoreFunctionOP scorefxn = scoring::getScoreFunction();
		(*scorefxn)( pose );
	}
	std::set< Size > designable_set;
	core::pack::task::PackerTaskOP design_task( tf->create_task_and_apply_taskoperations( pose ) );

#ifndef NDEBUG
	TR<< "Task for " << pose.pdb_info()->name() << " is: \n" << *(design_task)  << std::endl;
#endif

	// iterate over all residues
	for ( Size ii = 1; ii<= design_task->total_residue(); ++ii ) {
		if( design_task->being_designed( ii ) )
			designable_set.insert( ii );
	}

	return designable_set;

}


///@brief iterates over all designed positions and determines identity to native. outputs recoveries to file.
void measure_sequence_recovery( utility::vector1<core::pose::Pose> & native_poses, utility::vector1<core::pose::Pose> & redesign_poses ) {

	// setup main arrays used for calculation
	utility::vector1< core::Size > n_correct( chemical::num_canonical_aas, 0 );
	utility::vector1< core::Size > n_native( chemical::num_canonical_aas, 0 );
	utility::vector1< core::Size > n_designed( chemical::num_canonical_aas, 0 );

	utility::vector1< core::Size > n_correct_core( chemical::num_canonical_aas, 0 );
	utility::vector1< core::Size > n_native_core( chemical::num_canonical_aas, 0 );
	utility::vector1< core::Size > n_designed_core( chemical::num_canonical_aas, 0 );

	utility::vector1< core::Size > n_correct_surface( chemical::num_canonical_aas, 0 );
	utility::vector1< core::Size > n_native_surface( chemical::num_canonical_aas, 0 );
	utility::vector1< core::Size > n_designed_surface( chemical::num_canonical_aas, 0 );

	ObjexxFCL::FArray2D_int sub_matrix( chemical::num_canonical_aas, chemical::num_canonical_aas, 0 );

	Size n_correct_total(0); Size n_total(0);
	Size n_correct_total_core(0); Size n_total_core(0);
	Size n_correct_total_surface(0); Size n_total_surface(0);

	Size surface_exposed_cutoff = option[ sequence_recovery::se_cutoff ];
	Size core_cutoff = 24;

	// iterate through all the structures
	utility::vector1< core::pose::Pose >::iterator native_itr( native_poses.begin() ), native_last( native_poses.end() );
	utility::vector1< core::pose::Pose >::iterator redesign_itr( redesign_poses.begin() ), redesign_last( redesign_poses.end() );

	while( ( native_itr != native_last ) && (redesign_itr != redesign_last ) ) {

		// get local copies of the poses
		core::pose::Pose native_pose( *native_itr );
		core::pose::Pose redesign_pose( *redesign_itr );

		// figure out the task & neighbor info
		core::pack::task::TaskFactoryOP task_factory( new core::pack::task::TaskFactory );
		std::set< Size > design_set;
		utility::vector1< core::Size > num_neighbors;

		// setup what residues we are going to look at...
		setup_tf( task_factory );
		design_set = fill_designable_set( native_pose, task_factory );
		fill_num_neighbors( native_pose, num_neighbors );

		// record native sequence
		// native_sequence vector is sized for the WHOLE pose not just those being designed
		// it doesn't matter because we only iterate over the number of designed positions
		Size const nres( native_pose.total_residue() );
		utility::vector1< chemical::AA > native_sequence( nres );

		// iterate over designable positions
		for ( std::set< core::Size >::const_iterator it = design_set.begin(), end = design_set.end(); it != end; ++it ) {

			if ( ! native_pose.residue(*it).is_protein() ) {
				native_sequence[ *it ] = chemical::aa_unk;
				continue;
			}
			//figure out info about the native pose
			native_sequence[ *it ] = native_pose.residue( *it ).aa();
			n_native[ native_pose.residue(*it).aa() ]++;

			//determine core/surface
			if ( num_neighbors[*it] >= core_cutoff ) {
				n_native_core[ native_pose.residue(*it).aa() ]++;
				n_total_core++;
			}

			if ( num_neighbors[*it] < surface_exposed_cutoff ) {
				n_native_surface[ native_pose.residue(*it).aa() ]++;
				n_total_surface++;
			}

		} // end finding native seq

		/// measure seq recov
		for ( std::set< core::Size >::const_iterator it = design_set.begin(), end = design_set.end(); it != end; ++it ) {

			// don't worry about recovery of non-protein residues
			if ( redesign_pose.residue( *it ).is_protein() ) {
				n_total++;

				// increment the designed count
				n_designed[ redesign_pose.residue(*it).aa() ]++;

				if ( num_neighbors[*it] >= core_cutoff ) { n_designed_core[ redesign_pose.residue(*it).aa() ]++; }
				if ( num_neighbors[*it] < surface_exposed_cutoff ) { n_designed_surface[ redesign_pose.residue(*it).aa() ]++; }

				// then check if it's the same
				if ( native_sequence[ *it ] == redesign_pose.residue(*it).aa() ) {
					n_correct[ redesign_pose.residue(*it).aa() ]++;

					if ( num_neighbors[*it] >= core_cutoff ) {
						n_correct_core[ redesign_pose.residue(*it).aa() ]++;
						n_correct_total_core++;
					}
					if ( num_neighbors[*it] < surface_exposed_cutoff ) {
						n_correct_surface[ redesign_pose.residue(*it).aa() ]++;
						n_correct_total_surface++;
					}
					n_correct_total++;
				}

				// set the substitution matrix for this go round
				sub_matrix( native_pose.residue(*it).aa(), redesign_pose.residue(*it).aa() )++;
			}

		} // end measure seq reovery

		// increment iterators
		native_itr++; redesign_itr++;
	}

	// open sequence recovery file stream
	utility::io::ozstream outputFile( option[ sequence_recovery::seq_recov_filename ].value() ) ;

	// write header
	outputFile << "Residue\tNo.correct core\tNo.native core\tNo.designed core\tNo.correct/ No.native core\tNo.correct/ No.designed core\t"
				 << "No.correct\tNo.native\tNo.designed\tNo.correct/ No.native\tNo.correct/ No.designed\t"
				 << "Residue\tNo.correct surface\tNo.native surface\tNo.designed surface\tNo.correct/ No.native\tNo.correct/ No.designed" << std::endl;

	// write AA data
	for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) {

		outputFile << chemical::name_from_aa( chemical::AA(ii) ) << "\t"
					<< n_correct_core[ ii ] << "\t" << n_native_core[ ii ] << "\t" << n_designed_core[ ii ] << "\t";

		if ( n_native_core[ii] != 0 ) outputFile << F(4,2, (float)n_correct_core[ii]/n_native_core[ii] ) << "\t";
		else outputFile << "---\t";
		if ( n_designed_core[ii] != 0 ) outputFile << F(4,2, (float)n_correct_core[ii]/n_designed_core[ii] ) << "\t";
		else outputFile << "---\t";

		// debug
		//if ( n_native_core[ii] != 0 ) std::cout << F(4,2, (float)n_correct_core[ii]/n_native_core[ii] ) << "\t";
		//if ( n_designed_core[ii] != 0 ) std::cout << F(4,2, (float)n_correct_core[ii]/n_designed_core[ii] ) << "\t";

		outputFile << n_correct[ ii ] << "\t" << n_native[ ii ] << "\t" << n_designed[ ii ] << "\t";
		if ( n_native[ii] != 0 ) outputFile << F(4,2, (float)n_correct[ii]/n_native[ii] ) << "\t";
		else outputFile << "---\t";
		if ( n_designed[ii] != 0 ) outputFile << F(4,2, (float)n_correct[ii]/n_designed[ii] ) << "\t";
		else outputFile << "---\t";

		// debug
		//if ( n_native[ii] != 0 ) std::cout << F(4,2, (float)n_correct[ii]/n_native[ii] ) << "\t";
		//if ( n_designed[ii] != 0 ) std::cout << F(4,2, (float)n_correct[ii]/n_designed[ii] ) << "\t";

		outputFile << chemical::name_from_aa( chemical::AA(ii) ) << "\t"
							 << n_correct_surface[ ii ] << "\t" << n_native_surface[ ii ] << "\t" << n_designed_surface[ ii ] << "\t";

		if ( n_native_surface[ii] != 0 ) outputFile << F(4,2, (float)n_correct_surface[ii]/n_native_surface[ii] ) << "\t";
		else outputFile << "---\t";
		if ( n_designed_surface[ii] != 0 ) outputFile << F(4,2, (float)n_correct_surface[ii]/n_designed_surface[ii] ) << "\t";
		else outputFile << "---\t";

		// debug
		//if ( n_native_surface[ii] != 0 ) std::cout << F(4,2, (float)n_correct_surface[ii]/n_native_surface[ii] ) << "\t";
		//if ( n_designed_surface[ii] != 0 ) std::cout << F(4,2, (float)n_correct_surface[ii]/n_designed_surface[ii] ) << "\t";

		outputFile << std::endl;
	}

	// write totals
	outputFile << "Total\t"
				<< n_correct_total_core << "\t" << n_total_core << "\t\t" << F(5,3, (float)n_correct_total_core/n_total_core ) << "\t\t"
				<< n_correct_total << "\t" << n_total << "\t\t" << F(5,3, (float)n_correct_total/n_total ) << "\t\tTotal\t"
				<< n_correct_total_surface << "\t" << n_total_surface << "\t\t" << F(5,3, (float)n_correct_total_surface/n_total_surface )
				<< std::endl;


	// output the sequence substitution file
	utility::io::ozstream matrixFile( option[ sequence_recovery::sub_matrix_filename ].value() ) ; //defaults to submatrix.txt

	// write the header
	matrixFile << "AA_TYPE" << "\t" ;
	for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) {
		matrixFile << "nat_"<<chemical::name_from_aa( chemical::AA(ii) ) << "\t";
	}
	matrixFile<<std::endl;

	// now write the numbers
	for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) { //redesigns
		matrixFile << "sub_" << chemical::name_from_aa( chemical::AA(ii) );
		for ( Size jj = 1; jj <= chemical::num_canonical_aas; ++jj ) { //natives
			//std::cout<<"Native: "<< jj << " Sub: " << ii << "  Value: "<<sub_matrix( jj, ii ) << std::endl;
			matrixFile<< "\t" << sub_matrix( jj, ii );
		}
		matrixFile << std::endl;
	}

// 	///output the sequence substitution file with percent of native recovered
// 	utility::io::ozstream matrixFileratio( "submatrix.ratio.txt" ) ; //allow naming later
// 	//write the header
// 	matrixFileratio << "AA_TYPE" << "\t" ;
// 	for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) {
// 	  matrixFileratio << "nat_"<<chemical::name_from_aa( chemical::AA(ii) ) << "\t";
// 	}
// 	matrixFileratio<<std::endl;

// 	//now write the numbers
// 	for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) { //redesigns
// 	  matrixFileratio << "sub_" << chemical::name_from_aa( chemical::AA(ii) );
// 	  for ( Size jj = 1; jj <= chemical::num_canonical_aas; ++jj ) { //natives
// 	    //std::cout<<"Native: "<< jj << " Sub: " << ii << "  Value: "<<sub_matrix( jj, ii ) << std::endl;
// 	    matrixFileratio<< "\t"<< F(4,2, (float)sub_matrix( jj, ii )/sub_matrix( jj, jj ) );
// 	  }
// 	  matrixFileratio << std::endl;
// 	}

}


//@brief method which contains logic for calculating rotamer recovery. not implemented.
void measure_rotamer_recovery( utility::vector1<core::pose::Pose> & /*native_poses*/, utility::vector1<core::pose::Pose> & /*redesign_poses*/ ) {}


//@brief main method for the sequence recovery protocol
int main( int argc, char* argv[] ) {

	using utility::file::file_exists;
	using utility::file::FileName;

	option.add( sequence_recovery::native_pdb_list, "List of pdb files of the native structures." );
	option.add( sequence_recovery::redesign_pdb_list, "List of pdb files of the redesigned structures." );
	option.add( sequence_recovery::parse_taskops_file, "XML file which contains task operations to apply before measuring recovery (optional)" );
	option.add( sequence_recovery::rotamer_recovery, "Compare the rotamer recovery instead of sequence recovery." ).def( false );
	option.add( sequence_recovery::seq_recov_filename, "Name of file for sequence recovery output." ).def("sequencerecovery.txt");
	option.add( sequence_recovery::sub_matrix_filename, "Name of file substitution matrix output." ).def("submatrix.txt");
	option.add( sequence_recovery::se_cutoff, "Integer for how many nbs a residue must have less than or equal to to be considered surface exposed." ).def( 16 );


	devel::init( argc, argv );

	// changing this so that native_pdb_list and redesign_pdb_list do not have default values. giving these options can lead
	// to users measuring recovery against the wrong set of PDBs.
	if ( argc == 1 || !option[ sequence_recovery::native_pdb_list ].user() || !option[ sequence_recovery::redesign_pdb_list ].user() ) {
		init_usage_prompt( argv[0] );
		utility_exit_with_message_status( usage_string, 1 );
	}

	// read list file. open the file specified by the flag 'native_pdb_list' and read in all the lines in it
	std::vector< FileName > native_pdb_file_names;
	std::string native_pdb_list_file_name( option[ sequence_recovery::native_pdb_list ].value() );
	std::ifstream native_data( native_pdb_list_file_name.c_str() );
	std::string native_line;
	if ( !native_data.good() ) {
		utility_exit_with_message( "Unable to open file: " + native_pdb_list_file_name + '\n' );
	}
	while( getline( native_data, native_line ) ) {
		native_pdb_file_names.push_back( FileName( native_line ) );
	}

	native_data.close();

	// read list file. open the file specified by the flag 'redesign_pdb_list' and read in all the lines in it
	std::vector< FileName > redesign_pdb_file_names;
	std::string redesign_pdb_list_file_name( option[ sequence_recovery::redesign_pdb_list ].value() );
	std::ifstream redesign_data( redesign_pdb_list_file_name.c_str() );
	std::string redesign_line;
	if ( !redesign_data.good() ) {
		utility_exit_with_message( "Unable to open file: " + redesign_pdb_list_file_name + '\n' );
	}
	while( getline( redesign_data, redesign_line ) ) {
		redesign_pdb_file_names.push_back( FileName( redesign_line ) );
	}
	redesign_data.close();

	// check that the vectors are the same size. if not error out immediately.
	if ( native_pdb_file_names.size() != redesign_pdb_file_names.size() ) {
		utility_exit_with_message( "Size of native pdb list file: " + native_pdb_list_file_name + " does not equal size of redesign pdb list: " + redesign_pdb_list_file_name + "!\n" );
	}

	// iterate over both FileName vector and read in the PDB files
	utility::vector1< pose::Pose > native_poses;
	utility::vector1< pose::Pose > redesign_poses;

	std::vector< FileName >::iterator native_pdb( native_pdb_file_names.begin() ), native_last_pdb(native_pdb_file_names.end());
	std::vector< FileName >::iterator redesign_pdb( redesign_pdb_file_names.begin() ), redesign_last_pdb(redesign_pdb_file_names.end());

	while ( ( native_pdb != native_last_pdb ) && ( redesign_pdb != redesign_last_pdb ) ) {

		// check to make sure the file exists
		if ( !file_exists( *native_pdb ) ) {
			utility_exit_with_message( "Native pdb " + std::string(*native_pdb) + " not found! skipping" );
		}
		if ( !file_exists( *redesign_pdb ) ) {
			utility_exit_with_message( "Redesign pdb " + std::string(*redesign_pdb) + " not found! skipping" );
		}

		TR << "Reading in poses " << *native_pdb << " and " << *redesign_pdb << std::endl;
		core::pose::Pose native_pose, redesign_pose;
		core::import_pose::pose_from_pdb( native_pose, *native_pdb );
		core::import_pose::pose_from_pdb( redesign_pose, *redesign_pdb );

		native_poses.push_back( native_pose ); redesign_poses.push_back( redesign_pose );
		native_pdb++; redesign_pdb++;
	}


	if ( option[ sequence_recovery::rotamer_recovery ].value() ) {
		TR << "Measuring rotamer recovery"  << std::endl;
		measure_rotamer_recovery( native_poses, redesign_poses );
	} else {
		TR << "Measuring sequence recovery" << std::endl;
		measure_sequence_recovery( native_poses, redesign_poses );
	}

}

