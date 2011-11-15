// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file protocols/pmut_scan/point_mut_scan.cc
/// @brief A protocol that tries to find stability enhancing mutations
/// @author Ron Jacak

// Unit headers
#include <protocols/pmut_scan/PointMutScanDriver.hh>
#include <protocols/pmut_scan/Mutant.hh>

//project Headers
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>
#include <basic/options/util.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/init.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/graph/Graph.hh>
#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/import_pose/import_pose.hh>

#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunctionInfo.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <protocols/moves/PackRotamersMover.hh>
#include <protocols/moves/MinMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TaskAwareMinMover.hh>

#include <protocols/toolbox/pose_metric_calculators/NeighborsByDistanceCalculator.hh>
#include <protocols/toolbox/task_operations/RestrictToNeighborhoodOperation.hh>

// Utility Headers
#include <utility/file/FileName.hh>

// Numeric Headers

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>

// C++ headers
#include <iostream>
#include <fstream>
#include <string>

#ifdef USEMPI
/// MPI
#include <mpi.h>
#endif

// option key includes
#include <basic/options/keys/run.OptionKeys.gen.hh>

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
#include <core/graph/Graph.fwd.hh>
#include <core/graph/UpperEdgeGraph.fwd.hh>
#include <core/graph/UpperEdgeGraph.hh>
#include <core/graph/unordered_object_pool.fwd.hpp>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.fwd.hh>
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
#include <core/io/pdb/file_data.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/types.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.fwd.hh>
#include <core/pack/rotamer_set/RotamerCouplings.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/pack/task/IGEdgeReweightContainer.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/ResfileReader.fwd.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <core/pack/task/operation/TaskOperations.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/PDBPoseMap.fwd.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Remarks.fwd.hh>
#include <core/pose/Remarks.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricCalculatorBase.fwd.hh>
#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergyGraph.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/LREnergyContainer.fwd.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/SecondaryStructureWeights.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/hbonds/HBondOptions.fwd.hh>
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
#include <core/scoring/mm/MMBondAngleResidueTypeParamSet.fwd.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/jobdist/Jobs.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/MinMover.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/moves/PackRotamersMover.fwd.hh>
#include <protocols/moves/TaskAwareMinMover.fwd.hh>
#include <protocols/pmut_scan/Mutant.fwd.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborsByDistanceCalculator.fwd.hh>
#include <protocols/toolbox/task_operations/RestrictOperationsBase.fwd.hh>
#include <protocols/toolbox/task_operations/RestrictOperationsBase.hh>
#include <protocols/toolbox/task_operations/RestrictToNeighborhoodOperation.fwd.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/PyAssert.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
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
#include <utility/file/FileName.fwd.hh>
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
#include <utility/signals/SignalHub.fwd.hh>
#include <utility/signals/SignalHub.hh>
#include <utility/tag/Tag.fwd.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzTriple.fwd.hh>
#include <numeric/xyzTriple.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
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
#include <ObjexxFCL/FArray1D.fwd.hh>
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
#include <iomanip>
#include <iosfwd>
#include <istream>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <time.h>
#include <utility>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/prof.hh>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/pool/poolfwd.hpp>
#include <boost/unordered_map.hpp>

//Auto using namespaces
namespace std { } using namespace std; // AUTO USING NS
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
//Auto using namespaces end



using namespace basic::options;
using namespace basic::options::OptionKeys;

using namespace core;
using namespace core::pack::task::operation;
using namespace core::pack::task;

using namespace protocols;
using namespace ObjexxFCL::fmt;
using namespace utility;


namespace protocols {
namespace pmut_scan {


static basic::Tracer TR("protocols.pmut_scan.PointMutScanDriver");


///
/// @begin PointMutScanDriver::PointMutScanDriver
///
/// @brief
/// Main constructor for the class. What all does it do?
///
PointMutScanDriver::PointMutScanDriver( utility::vector1< std::string > & pdb_file_names, bool double_mutant_scan, std::string list_file, bool output_mutant_structures ) :
	pdb_file_names_( pdb_file_names ),
	double_mutant_scan_( double_mutant_scan ),
	mutants_list_file_( list_file ),
	output_mutant_structures_( output_mutant_structures )
{

#ifdef USEMPI
	tag_ = 1; // need to initialize the tag on all nodes to 1 or MPI_Send/_Recv calls start acting funny
#endif

	DDG_CUTOFF_ = -1.0;

	int mpi_rank( 0 ), mpi_nprocs( 1 );
#ifdef USEMPI
	MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );/* get current process id */
	MPI_Comm_size( MPI_COMM_WORLD, &mpi_nprocs );/* get number of processes */
#endif

	MPI_rank_ = (Size)( mpi_rank );
	MPI_nprocs_ = (Size)( mpi_nprocs );

	read_in_structures(); // all processes read in the structures

}


///
/// @begin PointMutScanDriver::~PointMutScanDriver
///
/// @brief
/// Destructor. What all needs to be done here?
///
PointMutScanDriver::~PointMutScanDriver() {}


///
/// @begin PointMutScanDriver::go()
///
/// @brief
/// Entry point for the pmut_scan protocol.  This is function the app calls to do the scan.
///
void PointMutScanDriver::go() {

	clock_t entire_starttime;
	if ( MPI_rank_ == 0 ) {
		// time the protocol, doesn't include the time spent reading in input structures.
		entire_starttime = clock();
	}

	TR << "go(): " << node_name( MPI_rank_ ) << std::endl;

	if ( MPI_rank_ == 0 ) {
		// set up the list of mutations that will be tried.
		// if the user specified a list, then do just those. if not, try all possible combinations of mutants.
		fill_mutations_list();
	}

	barrier(); // do we really want all processes to hold here?
	divide_up_mutations();

	barrier(); // do we really want all processes to hold here?
	make_mutants();

	barrier();
	if ( MPI_rank_ == 0 ) {
		clock_t entire_stoptime = clock();
		TR << "main(): whole protocol took " << ((double)entire_stoptime-entire_starttime) / CLOCKS_PER_SEC << " seconds" << std::endl;
		TR << "go(): DONE with pmut scan." << std::endl;
	}

}

///
/// @begin PointMutScanDriver::node_name()
///
std::string PointMutScanDriver::node_name( Size rank ) {

	if ( rank == 0 ) {
		return "master node";
	} else {
		std::stringstream r;
		r << "slave node " << rank;
		return r.str();
	}
}

///
/// @begin PointMutScanDriver::barrier()
///
/// Make all processes stop and wait here.
///
void PointMutScanDriver::barrier() {

#ifdef USEMPI
	MPI_Barrier( MPI_COMM_WORLD );
#endif
	std::cout.flush();

}


///
/// @begin PointMutScanDriver::read_in_structures()
///
/// @brief
/// Reads in the structure (or list of structures) specified on the command line to a class member variable. Create
/// Pose objects for all of the structures because we'll pass these out to the slave nodes later.
///
/// NOTE: This protocol assumes that if you pass multiple structures, they are all variants of the same structure and you
/// want to use all of them for the ddG calculation.
///
///
void PointMutScanDriver::read_in_structures() {

	//
	// read in all the PDB files into a vector of Pose objects
	//
	utility::vector1< std::string >::iterator input_pdb_filename, last_pdb;
	for ( input_pdb_filename = pdb_file_names_.begin(), last_pdb = pdb_file_names_.end(); input_pdb_filename != last_pdb; ++input_pdb_filename ) {
		pose::Pose pose;
		core::import_pose::pose_from_pdb( pose, *input_pdb_filename );
		input_poses_.push_back( pose );
	}

}


///
/// @begin PointMutScanDriver::fill_mutations_list
///
/// @brief
/// Determines whether the user specified a list of mutants or just wants to do a scan over all possible combinations.
///
/// If we're doing a scan over all possible mutations:
///
/// If we have a two residue protein, 1 and 2, there are 19 possible aa's we can mutate each residue to. Since 1 and 2
/// are independent, the number of possible double mutants is 19*19 = 361. It's easy to enumerate all of the possible
/// double mutants, but we want to come up with an efficient way of making those mutants and calculating the ddGs. The
/// easiest solution is to enumerate all of the possible double mutants, make that a work unit, and then distribute all
/// of the work units out to a large cluster. The downside to this approach is that several processors will end up making
/// the same mutation, at least in part. For example, the double mutant A1C A2C is similar to the mutant A1C A2D. In fact,
/// A1C will have to be paired not only with A2C and A2D, but also A2E, A2F and so on. It would be more efficient to make
/// a pose for A1C, and then go into a second for loop that tries A2C-A2W on that already mutated pose.
///
/// What's the outermost thing this protocol has to do.  For the single mutant scan, you have to try all 19 non-wt aas at
/// every position.  That lends itself to parallelization rather easily.  Each protein position is independent of the others
/// so you can have nres processes each testing the mutations at that position.  At most, each processor will test 19 mutations.
/// With double_mutants, you have to fix one mutation (eg. A1C) and then try all possible other mutations at all other
/// positions.  So, if you have nres processors, each processor will fix some position to 1 of 19 aas and then scan through
/// a mutant at all other positions. Let's assume we have a 10 residue protein.  Position 1 will mutate to 1 of 19 aas.
/// For each of 1 of those 19, we have to test 19 * 9 = 171 other mutations (at the other positions). That results in a
/// grand total of 3249 possibilites.  And that's only residue 1's mutants!  We also have to try to fix the 19 non-wt aas
/// for position 2 and try 19 * 8 = 152 mutations at the other locations for a total of 2888 mutations for just position
/// 2. 3: 19 * 19 * 7 = 2527. 4: 19 * 19 * 6 = 2166.  5: 19 * 19 * 5 = 1805. Continuing on in this fashion leads to a
/// grand grand total of 16245 possible double mutants in a 10 residue protein. Doing the same kind of protocol for a
/// 233 residue protein results in 9,841,221 possible double mutants!
///
/// Testing ~10 million mutants even on 512 cpus could take quite a bit of time. We really need to find a way to prune
/// down the number of possible mutants to just the ones that will be most interesting. I definitely could change it so
/// that if the two mutations are more than some number of Angstroms apart, then don't bother making that mutant and
/// scoring. The question then becomes how often you have a stabilizing first mutant, and then find a stabilizing (better
/// than -0.1) second mutant on the first structure that is more than xAng away. Probably happens often.
///
/// Another problem is that the parallelization is not balanced. Because we have directionality in the approach for
/// testing double mutants - for example, if we've already done 1AC 2AC we don't have to do 2AC 1AC - processor 1 which
/// handles all of the possible mutants at 1 and every other residue has to do way way less
///
///
/// For triple mutants, assuming a 10 residue protein there would be 19 * 19 * 19 * nres(nres+1)/2, or ~377,000, possible
/// mutants. The 233 residue antibody: 186,983,199 possible combinations.
///
///
///
void PointMutScanDriver::fill_mutations_list() {

	if ( !mutants_list_file_.empty() ) {
		read_mutants_list_file( mutants_list_file_ );
		return;
	}

	// otherwise, we're just going to do a scan over all mutations
	// this outer for loop is over all sets of mutations: either single mutants, double mutants, triple mutants, combinations
	// of single, double and triple mutants, etc.
	//utility::vector1< Mutant > stabilizing_mutants;
	//scan_for_mutations( input_poses, scorefxn, stabilizing_mutants, double_mutant_scan_ );

	Size no_double_mutants_possible = 0;
	Size no_double_mutants_excluded_for_distance = 0;

	// use the first structure to determine neighborship for all residues. this neighbor_graph will be used inside the
	// nested for loops to skip mutants that are on opposite sides of the protein.
	utility::vector1< utility::vector1< bool > > neighbors;
	calculate_neighbor_table( input_poses_[1], neighbors );

	pose::Pose & pose = input_poses_[1];
	Size n_residue = pose.n_residue();

	for ( Size resid1 = 1; resid1 <= n_residue; ++resid1 ) {

		// try every type at each position (well, except the native type at this position!)
		for ( Size aa_enum_index_a = 1; aa_enum_index_a <= chemical::num_canonical_aas; ++aa_enum_index_a ) {

//if ( resid1 > 1 ) { break; } // for debugging only

			if ( pose.residue( resid1 ).aa() == chemical::AA( aa_enum_index_a ) ) { continue; }
			if ( !pose.residue_type( resid1 ).is_protein() ) { continue; }

			MutationData md1(
				pose.residue( resid1 ).name1(),
				oneletter_code_from_aa( chemical::AA( aa_enum_index_a ) ),
				resid1,
				pose.pdb_info()->number( resid1 ),
				pose.pdb_info()->icode( resid1 ),
				pose.pdb_info()->chain( resid1 )
			);

			// only do a double mutant scan if the user asked for it
			if ( double_mutant_scan_ ) {

				// only need to iterate over higher indexed residues. can't make two mutations at the same position!
				for ( Size resid2 = resid1 + 1; resid2 <= n_residue; ++resid2 ) {

					// check to see if these residues are neighbors of each other. we don't want to make double mutants
					// where the mutants are on opposite sides of the protein.
					if ( neighbors[ resid1 ][ resid2 ] == false ) {
						no_double_mutants_possible += 19;
						no_double_mutants_excluded_for_distance += 19;
						continue;
					}

					// try every type at each position (well, except the native type at this position!)
					for ( Size aa_enum_index_b = 1; aa_enum_index_b <= chemical::num_canonical_aas; ++aa_enum_index_b ) {

						if ( pose.residue( resid2 ).aa() == chemical::AA( aa_enum_index_b ) ) { continue; }
						if ( !pose.residue_type( resid2 ).is_protein() ) { continue; }

						no_double_mutants_possible++;

						MutationData md2(
							pose.residue( resid2 ).name1(),
							oneletter_code_from_aa( chemical::AA( aa_enum_index_b ) ),
							resid2,
							pose.pdb_info()->number( resid2 ),
							pose.pdb_info()->icode( resid2 ),
							pose.pdb_info()->chain( resid2 )
						);

						Mutant m;
						m.add_mutation( md1 ); // the variable mutations is a vector of vectors!
						m.add_mutation( md2 ); // the variable mutations is a vector of vectors!

						all_mutants_.push_back( m );
						//TR << "fill_mutations_list(): adding mutation: " << m << std::endl;
					}

				} // all residues resid2

			} else {
				Mutant m;
				m.add_mutation( md1 ); // the variable mutations is a vector of vectors!
				all_mutants_.push_back( m );
				//TR << "fill_mutations_list(): adding mutation: " << m << std::endl;
			}
		}

	} // all residues resid1

	if ( MPI_rank_ == 0 ) {
		Size single_possible = 19 * n_residue;
		TR << "fill_mutations_list(): number single mutants possible: " << single_possible << std::endl;

		TR << "fill_mutations_list(): number double mutants possible: " << no_double_mutants_possible << std::endl;
		TR << "fill_mutations_list(): number double mutants excluded for distance: " << no_double_mutants_excluded_for_distance << std::endl;
	}

}


///
/// @begin PointMutScanDriver::read_mutants_list_file()
///
/// @brief
/// If the user specified mutants, it reads the lines in the mutant list file and parses those lines to get mutation
/// data and then saves them all to the class member variable.
/// Needs access to a pose to translate the lines in the mutations_list file to pose numbering
///
void PointMutScanDriver::read_mutants_list_file( std::string & list_file ) {

	std::ifstream data( list_file.c_str() );
	if ( !data.good() ) {
		utility_exit_with_message( "Unable to open mutations file: " + list_file + '\n' );
	}

	// read in all lines in file
	utility::vector1< std::string > mutant_file_lines;
	std::string line;
	while ( getline( data, line ) ) {
		if ( line.size() < 1 || line[0] == '#' ) continue; // skip comment lines
		mutant_file_lines.push_back( line );
	}
	data.close();


	// iterate over all the lines
	for ( Size ii=1; ii <= mutant_file_lines.size(); ++ii ) {
		std::string const & line( mutant_file_lines[ ii ] );
		std::istringstream iss( line );

		char wt_residue, mut_residue, chain;
		std::string position_code;

		Mutant m;

		// there might be more than one mutation per line!
		while ( iss.peek() && !iss.eof() ) {

			iss >> chain >> wt_residue >> position_code >> mut_residue;

			// check to see if an insertion code is present in the position_code string
			// if the string is made of all digits, no icode is present
			Size pdb_resnum; char icode = ' ';
			std::stringstream ss;

			if ( position_code.find_first_not_of("0123456789") == std::string::npos ) {
				icode = ' ';
				ss << position_code;
				ss >> pdb_resnum;

			} else {
				for ( std::string::iterator it = position_code.begin(); it < position_code.end(); ++it ) {
					if ( isdigit(*it) ) {
						ss << (*it);
					} else {
						icode = *it; // assumes that insertion code is only 1-letter!!
					}
				}
				ss >> pdb_resnum; // converts the ss buffer contents to a Size type
			}

			// figure out what the pose residue number for this residue is
			pose::Pose & pose = input_poses_[ 1 ];
			Size pose_resnum = (pose.pdb_info())->pdb2pose( chain, pdb_resnum, icode );

			if ( pose.residue( pose_resnum ).name1() != wt_residue ) {
				TR << "wt_residue: " << wt_residue << ", pdb resnum: " << pdb_resnum << ", pose resnum: " << pose_resnum
					<< ", residue at pose resnum: " << pose.residue( pose_resnum ).name1() << std::endl;
				utility_exit_with_message("Error. Wild-type residue given in mutatons_list file does not match input structure. Please try again.");
			}

			//TR << "Found mutation of " << wt_residue << " to " << mut_residue  << " at position " << pose_resnum << " (pdb chain: '" << chain << "', resnum: '" << pdb_resnum << "', icode: '" << icode << "')" << std::endl;

			MutationData md( wt_residue, mut_residue, pose_resnum, pdb_resnum, icode, chain );
			m.add_mutation( md ); // the variable mutations is a vector of vectors!

		} // done parsing line

		all_mutants_.push_back( m );

	} // end iterating over lines read from input file

}

///
/// @begin PointMutScanDriver::calculate_neighbor_table
///
/// @brief
/// Calculates the 10A neighbor graph using the given pose object and then sets values in a 2D array to indicate which
/// resids are neighbors.
///
void PointMutScanDriver::calculate_neighbor_table( pose::Pose & pose, utility::vector1< utility::vector1< bool > > & neighbors ) {

	// size the neighbors 2D table
	neighbors.resize( pose.n_residue(), utility::vector1< bool >( pose.n_residue(), false ) );

	// PointGraph is a one-way graph, which makes it somewhat annoying for iterating over neighbors of a certain
	// position. Only edges to higher-indexed nodes exist. So instead, make a graph which has all the edges at every
	// node to simplify iterating over all neighboring edges.
	core::conformation::PointGraphOP pg( new core::conformation::PointGraph ); // create graph
	core::conformation::residue_point_graph_from_conformation( pose.conformation(), *pg ); // create vertices
	core::conformation::find_neighbors( pg, 10.0 /* Angstrom cutoff */ ); // create edges

	// actually create the neighbor graph from the point graph
	core::graph::Graph neighbor_graph( pose.n_residue() );
	for ( Size r=1; r <= pose.total_residue(); ++r ) {
		for ( core::conformation::PointGraph::UpperEdgeListConstIter edge_iter = pg->get_vertex(r).upper_edge_list_begin(),
			edge_end_iter = pg->get_vertex(r).upper_edge_list_end(); edge_iter != edge_end_iter; ++edge_iter ) {
				neighbor_graph.add_edge(r, edge_iter->upper_vertex());
		}
	}

	for ( Size ii=1; ii <= pose.n_residue(); ++ii ) {

		conformation::Residue const & ii_rsd( pose.residue( ii ) );
		for ( core::graph::EdgeListConstIterator eli = neighbor_graph.get_node( ii )->const_edge_list_begin(),
			eli_end = neighbor_graph.get_node( ii )->const_edge_list_end(); eli != eli_end; ++eli ) {

			Size nb_resnum = (*eli)->get_other_ind( ii );
			if ( nb_resnum < ii ) { continue; } // only want higher indexed residues

			// check to see if any of the atoms on this neighboring residue "interact" with any atoms on the ii residue.
			// our definition of interact: one sc-sc atom pair within 4.5A (BK's suggestion)
			conformation::Residue const & jj_rsd( pose.residue( nb_resnum ) );

			for ( Size jja = jj_rsd.first_sidechain_atom(); jja <= jj_rsd.nheavyatoms(); ++jja ) {
				conformation::Atom const & jja_atom( jj_rsd.atom( jja ) );
				Vector const & jja_atom_xyz = jja_atom.xyz();

				for ( Size iia = ii_rsd.first_sidechain_atom(); iia <= ii_rsd.nheavyatoms(); ++iia ) {
					conformation::Atom const & iia_atom( ii_rsd.atom( iia ) );
					Vector const & iia_atom_xyz = iia_atom.xyz();

					if ( iia_atom_xyz.distance( jja_atom_xyz ) < 4.5 ) {
						neighbors[ ii ][ nb_resnum ] = true; // only set the upper half of the 2D table; i.e. res1 must always be < res2
						break;
					}

				} // ii rsd atoms

				if ( neighbors[ ii ][ nb_resnum ] ) {
					// already found an atom pair within 4.5A; no point in going through all the rest of jj rsd's atoms!
					break;
				}

			} // jj rsd atoms
		}
	}

}


///
/// @begin PointMutScanDriver::divide_up_mutations
///
/// @brief
/// This function takes the vector of all possible mutants and splits them up as evenly as possible among all the CPUs.
///
void PointMutScanDriver::divide_up_mutations() {

	//TR << "Node " << MPI_rank_ << ", entered method divide_up_mutations()" << std::endl;

	if ( MPI_rank_ == 0 ) {
		//utility::vector1< Mutant > all_mutants_;

		Size const num_mutants_per_cpu = all_mutants_.size() / MPI_nprocs_;
 		Size const nextra = all_mutants_.size() - ( num_mutants_per_cpu * MPI_nprocs_ );

		Size my_njobs = ( nextra >= 1 ? 1 : 0 ) + num_mutants_per_cpu;
		for ( Size ii = 1; ii <= my_njobs; ++ii ) {
			mutants_list_.push_back( all_mutants_[ ii ] );
		}

#ifdef USEMPI
		//TR << "divide_up_mutations(): number of nodes " << MPI_nprocs_ << std::endl;
		Size mutant_offset = my_njobs;

		// send the other nodes their mutations lists so they know what they'll be working on
		for ( Size node_index = 1; node_index < MPI_nprocs_; ++node_index ) {
			Size node_njobs = ( nextra > node_index ? 1 : 0 ) + num_mutants_per_cpu;
			MPI_Send( & node_njobs, 1, MPI_UNSIGNED_LONG, node_index, tag_, MPI_COMM_WORLD );

			for ( Size mutant_index = mutant_offset + 1; mutant_index <= mutant_offset + node_njobs; ++mutant_index ) {
				send_mutant_data_to_node( node_index, all_mutants_[ mutant_index ] );
			}
			mutant_offset += node_njobs;
		}

	} else {
		// slave node. need to receive work order from master node.
		Size my_njobs;
		MPI_Recv( & my_njobs, 1, MPI_UNSIGNED_LONG, 0, tag_, MPI_COMM_WORLD, & stat_ );

		//TR << "divide_up_mutations(): received my_njobs: '" << my_njobs << "'" << std::endl;
		mutants_list_.reserve( my_njobs );
		for ( Size ii = 1; ii <= my_njobs; ++ii ) {
			mutants_list_.push_back( receive_mutant_data_from_node( 0 ) );
		}
#endif
	}

#ifdef USEMPI
	sleep( MPI_rank_ ); // a crude way to order processes...
	for ( Size ii = 1; ii <= mutants_list_.size(); ++ii ) {
		char hostname[256];
		gethostname(hostname, sizeof(hostname));
		TR << "divide_up_pdbs(): mutation '" << mutants_list_[ ii ] << "' assigned to " << hostname << " (rank = " << MPI_rank_ << ")" << std::endl;
	}
#endif

}


#ifdef USEMPI
///
/// @begin PointMutScanDriver::send_mutant_data_to_node
///
/// @brief
/// Takes a Mutant and a destination and constructs the MPI_Send call.
///
void PointMutScanDriver::send_mutant_data_to_node( int destination, const protocols::pmut_scan::Mutant & m ) {

	int tag( 1 );

	// each particular mutant can have one, two or more mutations associated with it, make sure to send all of them!
	Size mutant_num_mutations = m.n_mutations();
	//TR << "sending mutant_num_mutations: " << mutant_num_mutations << " to node " << destination << std::endl;
	MPI_Send( & mutant_num_mutations, 1, MPI_UNSIGNED_LONG, destination, tag, MPI_COMM_WORLD );

	for ( utility::vector1< MutationData >::const_iterator iter = m.mutations_begin(); iter != m.mutations_end(); ++iter ) {

		char wt_residue = iter->wt_residue_;
		char mut_residue = iter->mut_residue_;
		//TR << "sending wt_residue: '" << wt_residue << "' and mut_residue: '" << mut_residue << "' to node " << destination << "." << std::endl;
		MPI_Send( & wt_residue, 1, MPI_CHAR, destination, tag, MPI_COMM_WORLD );
		MPI_Send( & mut_residue, 1, MPI_CHAR, destination, tag, MPI_COMM_WORLD );

		Size pose_resnum = iter->pose_resnum_;
		Size pdb_resnum = iter->pdb_resnum_;
		MPI_Send( & pose_resnum, 1, MPI_UNSIGNED_LONG, destination, tag, MPI_COMM_WORLD );
		MPI_Send( & pdb_resnum, 1, MPI_UNSIGNED_LONG, destination, tag, MPI_COMM_WORLD );

		char icode = iter->icode_;
		char chain = iter->chain_;
		//TR << "sending icode: '" << icode << "' and chain: '" << chain << "' to node " << destination << "." << std::endl;
		MPI_Send( & icode, 1, MPI_CHAR, destination, tag, MPI_COMM_WORLD );
		MPI_Send( & chain, 1, MPI_CHAR, destination, tag, MPI_COMM_WORLD );

	}

}

///
/// @begin PointMutScanDriver::receive_mutant_data_to_node
///
/// @brief
/// Receive mutant data from the master node.  First find out how many mutations are in this mutant and then actually
/// get the mutation data.
///
Mutant PointMutScanDriver::receive_mutant_data_from_node( int source ) {

	int tag( 1 );
	MPI_Status stat;

	Size num_mutations;
	MPI_Recv( & num_mutations, 1, MPI_UNSIGNED_LONG, source, tag, MPI_COMM_WORLD, & stat );
	//TR << "received mutant_num_mutations from node " << source << ": " << num_mutations <<  std::endl;

	Mutant m;
	for ( Size ii = 1; ii <= num_mutations; ++ii ) {

		char wt_residue, mut_residue;
		MPI_Recv( & wt_residue, 1, MPI_CHAR, source, tag, MPI_COMM_WORLD, & stat );
		MPI_Recv( & mut_residue, 1, MPI_CHAR, source, tag, MPI_COMM_WORLD, & stat );
		//TR << "received wt_residue: " << wt_residue << " and mut_residue: " << mut_residue << " from node " << source << "." << std::endl;

		Size pose_resnum = 1, pdb_resnum = 1;
		MPI_Recv( & pose_resnum, 1, MPI_UNSIGNED_LONG, source, tag, MPI_COMM_WORLD, & stat );
		MPI_Recv( & pdb_resnum, 1, MPI_UNSIGNED_LONG, source, tag, MPI_COMM_WORLD, & stat );

		char icode, chain;
		MPI_Recv( & icode, 1, MPI_CHAR, source, tag, MPI_COMM_WORLD, & stat );
		MPI_Recv( & chain, 1, MPI_CHAR, source, tag, MPI_COMM_WORLD, & stat );
		//TR << "received icode: '" << icode << "' and chain: '" << chain << "' from node " << source << "." << std::endl;

		MutationData md( wt_residue, mut_residue, pose_resnum, pdb_resnum, icode, chain );
		m.add_mutation( md );
	}

	//TR << "receive_mutant_data_from_node(): received mutant '" << m << "'" << std::endl;

	return m;
}

#endif


///
/// @begin PointMutScanDriver::make_mutants
///
/// @brief
/// Calls make_specific_mutant on all mutants assigned to this node.
/// Also responsible for creating the score function that's used for all mutants.
///
void PointMutScanDriver::make_mutants() {

	utility::vector1< pose::Pose > mutant_poses( input_poses_.size() ); // this will get set in the function below
	utility::vector1< pose::Pose > native_poses( input_poses_.size() );

	// create a scorefxn that will be used for all the mutants
	// (to enable hpatch scoring, the command line weights file flag will have to be used)
	// decompose bb hbond energies into pair energies
	//
	scoring::ScoreFunctionOP scorefxn = scoring::getScoreFunction();

	scoring::methods::EnergyMethodOptions energymethodoptions( scorefxn->energy_method_options() );
	energymethodoptions.hbond_options().decompose_bb_hb_into_pair_energies( true );
	scorefxn->set_energy_method_options( energymethodoptions );

	// print out a header to the terminal
	if ( MPI_rank_ == 0 ) {
		TR << A( "mutation" ) << X(3) << A( "average ddG" ) << X(3) << A( "average total energy" ) << std::endl;
	}

	for ( Size ii=1; ii <= mutants_list_.size(); ++ii ) {
		Mutant & m = mutants_list_[ ii ];

		//TR << "make_mutants(): making mutant: " << m << std::endl;

		// the make specific mutant function changes both the mutant and native poses. we want to start with our
		// original starting structures each time though. so we have to copy the input poses to some working native
		// and mutant poses vectors.
		for ( Size ii=1; ii <= input_poses_.size(); ++ii ) {
			mutant_poses[ ii ] = input_poses_[ ii ];
			native_poses[ ii ] = input_poses_[ ii ];
		}

		make_specific_mutant( mutant_poses, native_poses, scorefxn, m, "" );
		// this will result in the Mutant object 'm' being modified, and since m is a reference, the original mutants_list_
		// will be modified, as well.
	}

}


///
/// @begin PointMutScanDriver::make_specific_mutant
///
/// @brief
/// Function which takes in a mutation to make (either single, double or more) and calls itself recursively until the desired
/// structure is created. Useful for testing certain combinations of mutations (like putting two double mutants together)
/// without having to run an entire scan protocol that would cover those mutations.
///
void PointMutScanDriver::make_specific_mutant( utility::vector1< pose::Pose > & mutant_poses, utility::vector1< pose::Pose > & native_poses,
	scoring::ScoreFunctionOP scorefxn, Mutant & m, std::string mutation_string ) {

	//TR << "make_specific_mutant() called. mutant_poses.size(): " << mutant_poses.size() << ", native_poses.size(): " << native_poses.size()
	//	<< ", num mutations: " << m.n_mutations() << ", mutation_string: " << mutation_string << std::endl;

	// if the mutants vector has more than element, we have to take out the first element of the vector
	if ( m.n_mutations() > 1 ) {

		// need to make the first mutation and call this function recursively
		MutationData md = m.pop_mutation();

		// make the first mutation on the mutant_poses
		for ( Size ii = 1; ii <= native_poses.size(); ++ii ) {

			// make the specific mutation, but don't do any scoring; the scorefxn is needed for packing
			make_mutant_structure( mutant_poses[ ii ], native_poses[ ii ], md, scorefxn );

		}

		std::stringstream out;
		out << mutation_string;
		if ( mutation_string != "" ) { out << ","; }
		out << md.mutation_string();
		std::string updated_mutation_string = out.str();

		make_specific_mutant( mutant_poses, native_poses, scorefxn, m, updated_mutation_string );

	} else {
		// make the last mutation, calculate the ddG, and print out the results
		MutationData md = m.pop_mutation();

		//TR << "make_specific_mutant(): making final mutation: " << md << std::endl;

		Energy sum_mutant_scores = 0.0;
		Energy average_mutant_score = 0.0;

		Energy sum_native_scores = 0.0;
		Energy average_native_score = 0.0;

		utility::vector1< Real > native_poses_total_energies( native_poses.size() );
		utility::vector1< Real > mutant_poses_total_energies( native_poses.size() );

		for ( Size ii=1; ii <= native_poses.size(); ++ii ) {
			// make the specific mutation, but don't do any scoring; the scorefxn is needed for packing
			// send in the input_pose for the mutant. that way the mutant poses will be "returned" because mutant_poses
			// is actually a reference!
			make_mutant_structure( mutant_poses[ii], native_poses[ii], md, scorefxn );

			// score the created mutant structure
			pose::Pose & mutant_pose = mutant_poses[ ii ];
			Energy mutant_score = (*scorefxn)( mutant_pose );
			mutant_poses_total_energies[ ii ] = mutant_score;
			sum_mutant_scores += mutant_score;

			// score the update native structure
			pose::Pose & native_pose = native_poses[ ii ];
			Energy native_score = (*scorefxn)( native_pose );
			native_poses_total_energies[ ii ] = native_score;
			sum_native_scores += native_score;

			// of course, we want to output structures if the user specified a mutants file!
			if ( output_mutant_structures_ ) {
				std::stringstream out;
				out << mutation_string;
				if ( mutation_string != "" ) { out << ", "; }
				out << md.mutation_string();
				utility::file::FileName fn( pdb_file_names_[ ii ] );
				std::string mutant_filename = fn.base() + "." + out.str() + ".pdb";
				mutant_pose.dump_scored_pdb( mutant_filename, *(scorefxn()) );
			}

		}

		average_mutant_score = sum_mutant_scores / mutant_poses.size();
		average_native_score = sum_native_scores / native_poses.size();

		Real ddG_mutation = average_mutant_score - average_native_score;
		if ( ddG_mutation > DDG_CUTOFF_ ) {
			return;
		}

		std::stringstream out;
		out << mutation_string;
		if ( mutation_string != "" ) { out << ","; }
		out << md.mutation_string();
		std::string final_mutation_string = out.str();

		TR << final_mutation_string << X(3) << F( 9,3,ddG_mutation ) << X(3) << F( 9,2,average_mutant_score ) << std::endl;


		/*TR << "native poses total energies: ";
		for ( Size ii=1; ii <= native_poses_total_energies.size(); ++ii ) {
			TR << native_poses_total_energies[ ii ] << ", ";
		}
		TR << std::endl;
		TR << "mutant poses total energies: ";
		for ( Size ii=1; ii <= mutant_poses_total_energies.size(); ++ii ) {
			TR << mutant_poses_total_energies[ ii ] << ", ";
		}
		TR << std::endl;*/
		TR.flush_all_channels();


	} // end loop over all mutants

}


///
/// @begin PointMutScanDriver::make_mutant_structure
///
/// @brief
/// Given mutant and native pose references and the mutation to make, this function constructs all the necessary PackerTask
/// Operations and Movers to apply the mutation and repacking steps to both the mutant and native poses.
///
void PointMutScanDriver::make_mutant_structure( pose::Pose & mutant_pose, pose::Pose & native_pose, MutationData const & md, scoring::ScoreFunctionOP scorefxn ) {

	Size resid = md.pose_resnum();
	chemical::AA mut_aa = chemical::aa_from_oneletter_code( md.mut_residue() );

	// need to create a neighborhood by distance calculator so we can identify neighbors of the mutated residue
	std::stringstream out;
	out << md.mutation_string() << "_mutant_nb_calculator";
	std::string calculator_name = out.str();

	pose::metrics::PoseMetricCalculatorOP mutant_nb_calculator = new toolbox::pose_metric_calculators::NeighborsByDistanceCalculator( resid );
	pose::metrics::CalculatorFactory::Instance().register_calculator( calculator_name, mutant_nb_calculator );

	basic::MetricValue< std::set< Size > > mv_neighbors;
	mutant_pose.metric( calculator_name, "neighbors", mv_neighbors );
	std::set< Size > const neighbor_set( mv_neighbors.value() );

	//TR << "make_mutant_structure(): neighbor_set: ";
	//for ( std::set< Size >::iterator it = neighbor_set.begin() ; it != neighbor_set.end(); it++ ) {
	//	TR << *it << ", ";
	//}
	//TR << std::endl;

	TaskFactoryOP native_tf = new TaskFactory();
	TaskFactoryOP mutant_tf = new TaskFactory();

	// the restrict operation class (which in the end is just a TaskOperation) takes a calculator during construction. I've already
	// created that calculator above.  This operation will disable repacking and design at all positions except those in the neighborhood
	// of the mutated position.
	TaskOperationCOP nb_op = new toolbox::task_operations::RestrictToNeighborhoodOperation( calculator_name );
	native_tf->push_back( nb_op ); mutant_tf->push_back( nb_op );

	// extra task operations we want to also include
	// the restrict residue to repacking ops are used to make sure that only repacking and not design is done to the residues in the neighborhood
	InitializeFromCommandlineOP init_op = new InitializeFromCommandline();
	native_tf->push_back( init_op ); mutant_tf->push_back( init_op );

	IncludeCurrentOP ic_op = new IncludeCurrent();
	native_tf->push_back( ic_op ); mutant_tf->push_back( ic_op );

	RestrictResidueToRepackingOP mutant_repack_op = new RestrictResidueToRepacking();
	RestrictResidueToRepackingOP wt_repack_op = new RestrictResidueToRepacking(); // will include one extra residue to repack
	for ( Size ii = 1; ii <= mutant_pose.n_residue(); ++ii ) {
		// resid is the position on the original pose. ii is the position on the copy.
		if ( ii == resid ) {
			// do design on this position
			utility::vector1< bool > keep_canonical_aas( chemical::num_canonical_aas, false );
			keep_canonical_aas[ mut_aa ] = true;
			RestrictAbsentCanonicalAASOP restrict_op = new RestrictAbsentCanonicalAAS( ii, keep_canonical_aas );
			mutant_tf->push_back( restrict_op );
			wt_repack_op->include_residue( ii ); // for the wild type, don't design on the mutant resid - but do allow repacking
		} else {
			// make this position repackable only; because of the commutativity of packer task ops, only the residues that are in the neighborhood
			// of the mutant will be allowed to repack. the restrict to neighborhood op will disallow packing at all positions not near the mutant.
			mutant_repack_op->include_residue( ii );
			wt_repack_op->include_residue( ii );
		}
	}
	native_tf->push_back( wt_repack_op );
	mutant_tf->push_back( mutant_repack_op );

	//TR << "Finished creating all TaskOperation's and TaskFactory's. Creating MoveMap." << std::endl;

	kinematics::MoveMapOP movemap = new core::kinematics::MoveMap();
	std::set< core::Size >::const_iterator iter;
	for ( iter = neighbor_set.begin(); iter != neighbor_set.end(); iter++ ) {
		//movemap_->set_bb(i, true); // don't do any backbone minimization
		movemap->set_chi( *iter, true ); // but do minimize the side chains
	}
	//movemap->show( std::cout, mutant_pose.n_residue() );

	//TR << "Movemap created... Beginning repacking/minimization of mutant pose." << std::endl;

	// create an actual PackerTask from the TaskFactory
	pack::task::PackerTaskOP scan_task = mutant_tf->create_task_and_apply_taskoperations( mutant_pose );
	scan_task->num_to_be_packed();
	//TR << "mutant packer task: " << *scan_task << std::endl;  // generates a TON of output

	// now create the movers that will do the repacking and minimization
	moves::PackRotamersMoverOP mutant_repacker_mover = new moves::PackRotamersMover( scorefxn, scan_task, 2 ); // ndruns: 2
	moves::MinMoverOP min_mover = new moves::MinMover( movemap, scorefxn, option[ OptionKeys::run::min_type ].value(), 0.01, true ); // use nb_list: true
	moves::TaskAwareMinMoverOP task_aware_min_mover = new moves::TaskAwareMinMover( min_mover, mutant_tf );
	protocols::moves::SequenceMoverOP seq_mover = new protocols::moves::SequenceMover;
	seq_mover->add_mover( mutant_repacker_mover );
	seq_mover->add_mover( task_aware_min_mover );

	seq_mover->apply( mutant_pose );

	//TR << "Beginning repacking/minimization of wt pose." << std::endl;

	// create an actual PackerTask from the TaskFactory
	pack::task::PackerTaskOP wt_task = native_tf->create_task_and_apply_taskoperations( native_pose );

	// now create the movers that will do the repacking and minimization of the native structure
	moves::PackRotamersMoverOP native_pack_mover = new moves::PackRotamersMover( scorefxn, wt_task, 2 ); // ndruns: 2
	min_mover = new moves::MinMover( movemap, scorefxn, option[ OptionKeys::run::min_type ].value(), 0.01, true ); // use nb_list: true
	task_aware_min_mover = new moves::TaskAwareMinMover( min_mover, native_tf );
	seq_mover = new protocols::moves::SequenceMover;
	seq_mover->add_mover( native_pack_mover );
	seq_mover->add_mover( task_aware_min_mover );

	seq_mover->apply( native_pose );

	// this needs to get recreated each time around
	pose::metrics::CalculatorFactory::Instance().remove_calculator( calculator_name );
	pose::metrics::CalculatorFactory::Instance().clear_calculators();
	mutant_nb_calculator = NULL;


} // done with make_mutant_structure


// setters used by the unit tests only
void PointMutScanDriver::set_ddG_cutoff( Real threshold ) {
	DDG_CUTOFF_ = threshold;
}


///
/// @begin PointMutScanDriver::mutants_begin
///
/// @brief
/// returns a const iterator to the beginning of the Mutant data member variable vector
///
utility::vector1< Mutant >::const_iterator PointMutScanDriver::mutants_begin() const {
	return all_mutants_.begin();
}


///
/// @begin PointMutScanDriver::mutants_end
///
/// @brief
/// returns a const iterator to the end of the Mutant data member variable vector
///
utility::vector1< Mutant >::const_iterator PointMutScanDriver::mutants_end() const {
	return all_mutants_.end();
}


///
/// @begin PointMutScanDriver::n_mutants
///
/// @brief
/// returns the size of the Mutant data member variable vector
///
Size PointMutScanDriver::n_mutants() const {
	return all_mutants_.size();
}


} // namespace pmut_scan
} // namespace protocols

