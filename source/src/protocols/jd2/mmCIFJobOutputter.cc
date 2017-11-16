// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/mmCIFJobOutputter.cc
/// @brief  mmCIFJobOutputter writes mmCIF formatted outputs from finished Jobs in JD2
/// @author Steven Lewis smlewi@gmail.com

///Unit headers
#include <protocols/jd2/mmCIFJobOutputter.hh>
#include <protocols/jd2/mmCIFJobOutputterCreator.hh>
#include <protocols/jd2/Job.hh>

///Project headers
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/io/util.hh>
#include <core/io/StructFileRep.fwd.hh>
#include <core/io/StructFileRepOptions.hh>
#include <core/io/mmcif/cif_writer.hh>

///Utility headers
#include <basic/Tracer.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh>
//#include <utility/file/FileName.hh>
#include <core/types.hh>
#include <basic/options/option.hh>

///Basic headers
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>
#include <basic/datacache/CacheableStringFloatMap.hh>
#include <basic/datacache/CacheableStringMap.hh>

///C++ headers
#include <string>
#include <map>
#include <sstream>

// option key includes

#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <utility/vector1.hh>


static basic::Tracer TR( "protocols.jd2.mmCIFJobOutputter" );

namespace protocols {
namespace jd2 {

protocols::jd2::mmCIFJobOutputter::mmCIFJobOutputter()
: parent()
{
	using namespace basic::options::OptionKeys;
	using basic::options::option;

	TR.Debug << "Using mmCIFJobOutputter for JobDistributor" << std::endl;

	set_extension(".cif");
	//Not clear we should support mmCIFgz output, given that we don't support its input (I think)
	if ( option[ out::mmCIF_gz ] ) {
		set_extension(".cif.gz");
	}

	if ( option[ out::path::mmCIF ].user() ) {
		set_path(option[ out::path::mmCIF ]().path());
	}
}

protocols::jd2::mmCIFJobOutputter::~mmCIFJobOutputter()= default;

/// @details private function (just prevents code duplication) to fill ozstream
void protocols::jd2::mmCIFJobOutputter::dump_pose(
	JobCOP job,
	core::pose::Pose const & pose,
	utility::io::ozstream & out,
	std::string const & filename /* filename is an optional label in the score data table */
)
{

	//cif_extra_data_separate_file says, DO print the data, but NOT in the mmCIF
	bool const cif_extra_data_separate_file( basic::options::option[ basic::options::OptionKeys::out::file::cif_extra_data_separate_file ] );

	core::io::StructFileRepOptionsOP options =  core::io::StructFileRepOptionsOP( new core::io::StructFileRepOptions );

	//Turn off extra output if we put it in a separate file.
	if ( cif_extra_data_separate_file ) {
		options->set_output_pose_energies_table( false );
		options->set_output_pose_cache_data( false );
	}

	//This actually writes the cif, of course
	core::io::StructFileRepOP sfr = core::io::mmcif::dump_cif(
		pose,
		out,
		options
	);

	//...then, depending on cif_extra_data_separate_file, we may write a bunch more data.
	if ( cif_extra_data_separate_file ) {
		//This call to out.filename() is a stupid way to do this, but it ensures this extra data goes exactly where the cif file did, which otherwise requires burrowing a duplicate copy of a "tag" string through a few layers/classes that don't care about it and duplicating the file name generation calls.  We can rewrite the interfaces to pass tag along instead if necessary.

		sfr->score_table_filename() = filename;
		std::string const & parent_filename(out.filename());
		dump_extra_data_file(job, /*pose,*/ parent_filename);
		dump_energies_file(/*job,*/ sfr, parent_filename);
	}

	return;
}

///@details this function takes "extra data" associated with the job and  writes it to JOBNAME.extradata
void protocols::jd2::mmCIFJobOutputter::dump_extra_data_file(
	JobCOP job,
	//core::pose::Pose const & pose,
	std::string const & parent_filename
)
{

	std::string const extra_data_filename(parent_filename + ".extra_data"); //I assume the compiler can optimize this out
	//possibly should check if user wanted gzip, but I don't expect these files to be large
	utility::io::ozstream out( extra_data_filename );
	out << extract_data_from_Job(job) << std::endl;
	out.close();

	return;
}

///@details this function takes energies from the pose and writes it to JOBNAME.energies
void protocols::jd2::mmCIFJobOutputter::dump_energies_file(
	//JobCOP job,
	core::io::StructFileRepCOP sfr,
	std::string const & parent_filename
)
{
	std::string const energies_filename(parent_filename + ".energies"); //I assume the compiler can optimize this out
	//possibly should check if user wanted gzip, but I don't expect these files to be large
	utility::io::ozstream out( energies_filename );
	//arguably duplicated in that this happens separately in the mmCIF writer
	out << core::io::pose_energies_from_sfr(*sfr) << std::endl;
	out << core::io::pose_data_cache_from_sfr(*sfr) << std::endl;
	out.close();

	return;
}

//CREATOR SECTION
std::string
mmCIFJobOutputterCreator::keyname() const
{
	return "mmCIFJobOutputter";
}

protocols::jd2::JobOutputterOP
mmCIFJobOutputterCreator::create_JobOutputter() const {
	return protocols::jd2::JobOutputterOP( new mmCIFJobOutputter );
}

}//jd2
}//protocols
