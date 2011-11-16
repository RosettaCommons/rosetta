// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/JobDistributorFactory
/// @brief  JobDistributorFactory class, part of August 2008 job distributor as planned at RosettaCon08.
/// @author Andrew Leaver-Fay
/// @author Steven Lewis smlewi@gmail.com

// Unit headers
#include <protocols/jd2/JobDistributorFactory.hh>

// Package heaaders
#include <protocols/jd2/FileSystemJobDistributor.hh>
#include <protocols/jd2/ShuffleJobDistributor.hh>
#include <protocols/jd2/BOINCJobDistributor.hh>
#include <protocols/jd2/MPIWorkPartitionJobDistributor.hh>
#include <protocols/jd2/MPIWorkPoolJobDistributor.hh>
#include <protocols/jd2/MPIFileBufJobDistributor.hh>
#include <protocols/jd2/archive/MPIArchiveJobDistributor.hh>

#include <protocols/jd2/PDBJobInputter.hh>
#include <protocols/jd2/GenericJobInputter.hh>
#include <protocols/jd2/SilentFileJobInputter.hh>
#include <protocols/jd2/LazySilentFileJobInputter.hh>
#include <protocols/jd2/ThreadingJobInputter.hh>
#include <protocols/jd2/DatabaseJobInputter.hh>
#include <protocols/jd2/PoseInputStreamJobInputter.hh>

#include <protocols/jd2/PDBJobOutputter.hh>
#include <protocols/jd2/SilentFileJobOutputter.hh>
#include <protocols/jd2/AtomTreeDiffJobOutputter.hh>
#include <protocols/jd2/AtomTreeDiffJobInputter.hh>
#include <protocols/jd2/NoOutputJobOutputter.hh>
#include <protocols/jd2/ScoreOnlyJobOutputter.hh>
#include <protocols/jd2/EnzdesJobOutputter.hh>
#include <protocols/jd2/DatabaseJobOutputter.hh>

#include <protocols/jd2/DockDesignParser.hh>
#include <protocols/protein_interface_design/ParserJobInputter.hh>
#include <protocols/jd2/EnzdesJobInputter.hh>

#include <basic/options/option.hh>
// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>

// Utility headers
#include <utility/exit.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
// AUTO-REMOVED #include <core/conformation/RotamerSetBase.fwd.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/import_pose/atom_tree_diffs/atom_tree_diff.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/io/silent/SharedSilentData.hh>
#include <core/io/silent/SilentEnergy.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
// AUTO-REMOVED #include <core/scoring/EnergyGraph.fwd.hh>
// AUTO-REMOVED #include <core/scoring/EnergyMap.fwd.hh>
// AUTO-REMOVED #include <core/scoring/EnergyMap.hh>
// AUTO-REMOVED #include <core/scoring/LREnergyContainer.fwd.hh>
// AUTO-REMOVED #include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionInfo.fwd.hh>
// AUTO-REMOVED #include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/ContextDependentOneBodyEnergy.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/ContextDependentTwoBodyEnergy.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/ContextIndependentOneBodyEnergy.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/ContextIndependentTwoBodyEnergy.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/EnergyMethod.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/LongRangeTwoBodyEnergy.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/TwoBodyEnergy.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/WholeStructureEnergy.fwd.hh>
#include <core/sequence/SequenceAlignment.fwd.hh>
#include <protocols/enzdes/EnzFilters.fwd.hh>
#include <protocols/enzdes/EnzdesLoopsFile.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.fwd.hh>
// AUTO-REMOVED #include <protocols/features/DatabaseFilters.fwd.hh>
// AUTO-REMOVED #include <protocols/features/FeaturesReporter.fwd.hh>
// AUTO-REMOVED #include <protocols/features/FeaturesReporter.hh>
// AUTO-REMOVED #include <protocols/features/JobDataFeatures.fwd.hh>
// AUTO-REMOVED #include <protocols/features/PdbDataFeatures.fwd.hh>
// AUTO-REMOVED #include <protocols/features/PoseCommentsFeatures.fwd.hh>
// AUTO-REMOVED #include <protocols/features/PoseConformationFeatures.fwd.hh>
// AUTO-REMOVED #include <protocols/features/ProteinResidueConformationFeatures.fwd.hh>
#include <protocols/features/ProteinSilentReport.fwd.hh>
// AUTO-REMOVED #include <protocols/features/ProteinSilentReport.hh>
// AUTO-REMOVED #include <protocols/features/ProtocolFeatures.fwd.hh>
// AUTO-REMOVED #include <protocols/features/Report.fwd.hh>
// AUTO-REMOVED #include <protocols/features/Report.hh>
// AUTO-REMOVED #include <protocols/features/ResidueConformationFeatures.fwd.hh>
// AUTO-REMOVED #include <protocols/features/ScoreTypeFeatures.fwd.hh>
// AUTO-REMOVED #include <protocols/features/StructureFeatures.fwd.hh>
// AUTO-REMOVED #include <protocols/features/StructureScoresFeatures.fwd.hh>
// AUTO-REMOVED #include <protocols/filters/Filter.fwd.hh>
#include <protocols/jd2/AtomTreeDiffJobInputter.fwd.hh>
#include <protocols/jd2/AtomTreeDiffJobOutputter.fwd.hh>
#include <protocols/jd2/DatabaseJobInputter.fwd.hh>
#include <protocols/jd2/DatabaseJobOutputter.fwd.hh>
#include <protocols/jd2/EnzdesJobInputter.fwd.hh>
#include <protocols/jd2/EnzdesJobOutputter.fwd.hh>
#include <protocols/jd2/FileJobOutputter.fwd.hh>
#include <protocols/jd2/FileJobOutputter.hh>
#include <protocols/jd2/FileSystemJobDistributor.fwd.hh>
#include <protocols/jd2/GenericJobInputter.fwd.hh>
#include <protocols/jd2/InnerJob.fwd.hh>
#include <protocols/jd2/Job.fwd.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.fwd.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobDistributorFactory.fwd.hh>
#include <protocols/jd2/JobInputter.fwd.hh>
#include <protocols/jd2/JobInputter.hh>
#include <protocols/jd2/JobOutputter.fwd.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/MPIFileBufJobDistributor.fwd.hh>
#include <protocols/jd2/MPIWorkPartitionJobDistributor.fwd.hh>
#include <protocols/jd2/MPIWorkPoolJobDistributor.fwd.hh>
#include <protocols/jd2/PDBJobInputter.fwd.hh>
#include <protocols/jd2/PDBJobOutputter.fwd.hh>
#include <protocols/jd2/Parser.fwd.hh>
#include <protocols/jd2/Parser.hh>
#include <protocols/jd2/PoseInputStreamJobInputter.fwd.hh>
#include <protocols/jd2/ScoreOnlyJobOutputter.fwd.hh>
#include <protocols/jd2/SilentFileJobInputter.fwd.hh>
#include <protocols/jd2/SilentFileJobOutputter.fwd.hh>
#include <protocols/jd2/ThreadingJobInputter.fwd.hh>
#include <protocols/jd2/archive/MPIArchiveJobDistributor.fwd.hh>
// AUTO-REMOVED #include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MoverFactory.fwd.hh>
#include <protocols/protein_interface_design/ParserJobInputter.fwd.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/down_cast.hh>
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
#include <utility/file/gzip_util.hh>
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
// AUTO-REMOVED #include <utility/sql_database/DatabaseSessionManager.fwd.hh>
// AUTO-REMOVED #include <utility/sql_database/DatabaseSessionManager.hh>
// AUTO-REMOVED #include <utility/tag/Tag.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <deque>
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
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <boost/bind.hpp>
#include <boost/function.hpp>
// AUTO-REMOVED #include <boost/scoped_ptr.hpp>
// AUTO-REMOVED #include <cppdb/frontend.h>
#include <zlib/zlib.h>
#include <zlib/zutil.h>


//Haiku is classy
#ifdef BOINC
#ifdef USEMPI
Compilation stops
MPI or BOINC: pick one
check your build settings
#endif
#endif


namespace protocols {
namespace jd2 {

using namespace basic::options;
using namespace basic::options::OptionKeys;
/// @details All the convoluted logic combining compile-time decisions and
/// run-time decisions for which job distributor to use lives here.
/// As of right now, this logic isn't all that convoluted.
JobDistributor *
JobDistributorFactory::create_job_distributor() {
#ifdef USEMPI
	int npes_;
	MPI_Comm_size( MPI_COMM_WORLD, ( int* )( &npes_ ) );

	if ( npes_ > 3 && option[ OptionKeys::run::archive ] ) return new archive::MPIArchiveJobDistributor;
	if ( npes_ > 2 && option[ OptionKeys::jd2::mpi_work_partition_job_distributor ].value() == true ) {
		return new MPIWorkPartitionJobDistributor;
	} else {

		///NOTE: MPIFileBufJobDistributor has not been tested with PDB-Outputter. No idea if this works.
		/// according to wishes of the community the check in the lines below is turned off.
		/// if you see that MPI and PDB output is not working: you have a couple of options:
		///    1) change to silent-output it is better
		///    2) change to silent-output it is way better
		///    3) dude: change to silent-output
		///    4) debug PDBOutputter to make it work with MPIFileBufJobDistributor --- and remove this note
		///    5) uncomment the test for -out:file:silent below
		///   OL 6/9/09

		///NOTE: PDB output does not work with MPIFileBufJobDistributor as of 6/15/09 (SML)
		///Someone should debug MPIFileBufJobDistributor to follow the parent classes and work with all
		///JobInputters and JobOutputters.  The silent check below can then be removed, since it's dangerous overloading
		///to hack around a bug...

		///NOTE: the option jd2:mpi_filebuf_jobdistributor is now true by default, according to exchange on minirosetta list Nov 2010. (OL)
		if ( npes_ > 2 && option[ out::file::silent ].user() &&
			( option[ OptionKeys::jd2::mpi_file_buf_job_distributor ].value() == true
				|| option[ OptionKeys::jd2::mpi_filebuf_jobdistributor ].value() == true ) ) {
			return new MPIFileBufJobDistributor;
		} else {
			if ( npes_ > 1 ) return new MPIWorkPoolJobDistributor;
		}
	}
#endif

#ifdef BOINC
	return new BOINCJobDistributor;
#endif

	if ( option[  OptionKeys::run::shuffle ]() ){
		 return new ShuffleFileSystemJobDistributor;
	}

	if ( option[ OptionKeys::run::proc_id ].user()
		&& option [ OptionKeys::run::nproc ]() > 1 ) {
		return new MPIWorkPartitionJobDistributor;
	}

// #else
// 	if ( option[ n_worker_threads ].user() ) {
// 		return new MultiThreadedJobDistributor;
// 	} else {
// 		return new FileSystemJobDistributor;
// 	}
	return new FileSystemJobDistributor; //SML override until we have other child classes
}


/// @details All the logic for determining job input type lives here. Logic for
/// other stuff doesn't have to go home, but it can't live here ...
JobInputterOP
JobDistributorFactory::create_job_inputter() {
	if ( option[ basic::options::OptionKeys::jd2::pose_input_stream ]() ) {
		return new PoseInputStreamJobInputter;
	}

	if ( option[ in::file::s ].user() || option[ in::file::l ].user() || option[ in::file::list ].user() ) {
		if ( option[ basic::options::OptionKeys::enzdes::parser_read_cloud_pdb ].user() ) return new EnzdesJobInputter;
		if ( option[ basic::options::OptionKeys::jd2::dd_parser ].user() && option[ basic::options::OptionKeys::parser::patchdock ].user() )
			return new protocols::protein_interface_design::ParserJobInputter;
		else
			return new PDBJobInputter; //SML override until we have other child classes
// 		//if ( option[ in::file::zip ] ) {
// 		//	return new GzippedPDBJobInputter; // unnecessary since the izstream handles this seamlessly
// 		//} else {
// 			return new PDBJobInputter;
// 		//}
 	} else if ( option[ in::file::silent ].user() ) {
		if ( option[ OptionKeys::jd2::lazy_silent_file_reader ].user() ){
			return new LazySilentFileJobInputter;
		} else {
			return new SilentFileJobInputter;
		}
 	} else if (option[in::file::atom_tree_diff].user() ){
 			return new AtomTreeDiffJobInputter;
 	} else if ( option[ in::file::template_pdb ].user() || option[ in::file::template_silent ].user() ) {
		return new ThreadingJobInputter;
	} else if (option[in::use_database].user() ){
		return new DatabaseJobInputter;
	} else {
		// abinitio doesn't start with a pdb or silent file, it starts with a fasta!
		// That's missing the point of the jd2 - code for making a new Pose should
		// in one place! The logic for making a pose from a sequence could live
		// in it's own JobDistributor, so that you could run abinitio fragment
		// assembly on Poses that aren't necessarily extended.
		return new GenericJobInputter; //handles -nstruct alone
	}
}

/// @details this function handles the runtime + compiletime determination of
/// which JobOutputter to use
JobOutputterOP
JobDistributorFactory::create_job_outputter() {

	if ( option[ out::file::silent ].user()  ) {
		return new SilentFileJobOutputter;
	} else if (option[out::file::atom_tree_diff].user() ){
		return new AtomTreeDiffJobOutputter;
	}else if (option[out::file::score_only].user()) {
		return new ScoreOnlyJobOutputter;
	} else if ( option[ basic::options::OptionKeys::jd2::no_output ].value() || option[ out::nooutput ] ){
		return new NoOutputJobOutputter;
	} else if ( option[ basic::options::OptionKeys::jd2::enzdes_out].user() ){
		return new EnzdesJobOutputter;
	} else if ( option[ basic::options::OptionKeys::out::use_database].user() ){
		return new DatabaseJobOutputter;
	}	else {
		return new PDBJobOutputter;
	}
	return new PDBJobOutputter; //SML override until we have other child classes
}

/// @details this function handles the runtime + compiletime determination of
/// which JobOutputter to use
JobOutputterOP
JobDistributorFactory::create_job_outputter( JobOutputterOP default_jobout ) {

	if ( option[ out::file::silent ].user()  ) {
		return new SilentFileJobOutputter;
	} else if (option[out::pdb].user() ){
		return new PDBJobOutputter;
	} else if (option[out::file::atom_tree_diff].user() ){
		return new AtomTreeDiffJobOutputter;
	}else if (option[out::file::score_only].user()) {
		return new ScoreOnlyJobOutputter;
	} else if ( option[ basic::options::OptionKeys::jd2::no_output ].value() || option[ out::nooutput ] ){
		return new NoOutputJobOutputter;
	} else if ( option[ basic::options::OptionKeys::jd2::enzdes_out].user() ){
		return new EnzdesJobOutputter;
	} else if ( option[ basic::options::OptionKeys::out::use_database].user() ){
		return new DatabaseJobOutputter;
	}	else {
		return default_jobout;
	}
	return default_jobout; //SML override until we have other child classes
}

/// @details this function handles the determination of which Parser is required
/// (if any; returning NULL is valid if no parser is desired)
ParserOP
JobDistributorFactory::create_parser()
{
	if ( option[ OptionKeys::jd2::dd_parser ].user() )
		return new DockDesignParser;
	return NULL;
}

} // namespace jd2
} // namespace protocols
