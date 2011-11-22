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
#include <utility/vector1.hh>
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
