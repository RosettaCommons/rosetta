// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/mmtfJobOutputter.cc
/// @brief  mmtfJobOutputter writes mmtf formatted outputs from finished Jobs in JD2
/// @author Danny Farrell danpf@uw.edu

///Unit headers
#include <protocols/jd2/mmtfJobOutputter.hh>
#include <protocols/jd2/mmtfJobOutputterCreator.hh>
#include <protocols/jd2/Job.hh>

///Project headers
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/io/util.hh>
#include <core/io/StructFileRep.fwd.hh>
#include <core/io/StructFileRepOptions.hh>
#include <core/io/mmtf/mmtf_writer.hh>

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


static basic::Tracer TR( "protocols.jd2.mmtfJobOutputter" );

namespace protocols {
namespace jd2 {

protocols::jd2::mmtfJobOutputter::mmtfJobOutputter()
: parent()
{
	using namespace basic::options::OptionKeys;
	using basic::options::option;

	TR.Debug << "Using mmtfJobOutputter for JobDistributor" << std::endl;

	set_extension(".mmtf");
	//Not clear we should support mmtfgz output, given that we don't support its input (I think)
	if ( option[ out::mmtf_gz ] ) {
		set_extension(".mmtf.gz");
	}

	if ( option[ out::path::mmtf ].user() ) {
		set_path(option[ out::path::mmtf ]().path());
	}
}

protocols::jd2::mmtfJobOutputter::~mmtfJobOutputter()= default;

/// @details private function (just prevents code duplication) to fill ozstream
void protocols::jd2::mmtfJobOutputter::dump_pose(
	JobCOP,
	core::pose::Pose const & pose,
	utility::io::ozstream & out,
	std::string const &
)
{
	// Currently dont have a way to write extra data to mmtf
	core::io::StructFileRepOptionsOP options =  utility::pointer::make_shared< core::io::StructFileRepOptions >();
	core::io::mmtf::set_mmtf_default_options(*options);
	core::io::StructFileRepOP sfr = core::io::mmtf::dump_mmtf(
		pose,
		out,
		options
	);

	return;
}

///@details this function takes "extra data" associated with the job and  writes it to JOBNAME.extradata
void protocols::jd2::mmtfJobOutputter::dump_extra_data_file(
	JobCOP job,
	// core::pose::Pose const & pose,
	std::string const & parent_filename
)
{

	std::string const extra_data_filename(parent_filename + ".extra_data");
	// possibly should check if user wanted gzip, but I don't expect these files to be large
	utility::io::ozstream out( extra_data_filename );
	out << extract_data_from_Job(job) << std::endl;
	out.close();

	return;
}

///@details this function takes energies from the pose and writes it to JOBNAME.energies
void protocols::jd2::mmtfJobOutputter::dump_energies_file(
	//JobCOP job,
	core::io::StructFileRepCOP sfr,
	std::string const & parent_filename
)
{
	std::string const energies_filename(parent_filename + ".energies");
	//possibly should check if user wanted gzip, but I don't expect these files to be large
	utility::io::ozstream out( energies_filename );
	//arguably duplicated in that this happens separately in the mmtf writer
	out << core::io::pose_energies_from_sfr(*sfr) << std::endl;
	out << core::io::pose_data_cache_from_sfr(*sfr) << std::endl;
	out.close();

	return;
}

//CREATOR SECTION
std::string
mmtfJobOutputterCreator::keyname() const
{
	return "mmtfJobOutputter";
}

protocols::jd2::JobOutputterOP
mmtfJobOutputterCreator::create_JobOutputter() const {
	return utility::pointer::make_shared< mmtfJobOutputter >();
}

}//jd2
}//protocols
