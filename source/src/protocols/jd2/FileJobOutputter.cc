// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/FileJobOutputter.cc
/// @brief  August 2008 job distributor as planned at RosettaCon08 - FileJobOutputter
/// @author Steven Lewis smlewi@gmail.com

///Unit headers
#include <protocols/jd2/FileJobOutputter.hh>
#include <protocols/jd2/Job.hh>

///Project headers
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/io/raw_data/ScoreFileData.hh>
#include <protocols/jd2/ScoreMap.hh>
#include <core/scoring/Energies.hh>

///Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <utility/io/ozstream.hh>
#include <utility/vector1.hh>

///Basic headers
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>
#include <basic/datacache/CacheableStringFloatMap.hh>
#include <basic/datacache/CacheableStringMap.hh>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>



///C++ headers
//#include <string>
//#include <sstream>

static basic::Tracer TR("protocols.jd2.FileJobOutputter");

namespace protocols {
namespace jd2 {

using namespace basic::options;
using namespace basic::options::OptionKeys;

FileJobOutputter::FileJobOutputter() : parent() {
	set_defaults();
}

void FileJobOutputter::set_defaults() {
	write_scorefile_ = false;
	scorefile_name_ = "";
	TR.Debug << "FileJobOutputter ctor" << std::endl;
	utility::file::FileName default_path( option[ out::path::all ]() );
	if ( option[ out::file::scorefile ].user() ) {
		write_scorefile_ = true;
		scorefile_name_ = option[ out::file::scorefile ]();
		if ( option[ out::path::score ].user() ) {
			scorefile_name_.path( option[ out::path::score ]().path() );
			} else if (! scorefile_name_.absolute() ) {
					scorefile_name_.path( default_path.path() + "/" + scorefile_name_.path() );
			}
	} else if ( option [out::file::score_only].user() ) {
		write_scorefile_ = true;
		scorefile_name_ = option[out::file::score_only]();
		if ( !scorefile_name_.absolute() ) scorefile_name_.path( default_path.path() + "/" + scorefile_name_.path() );
	} else if ( !option[ run::no_scorefile ]() ) {
		write_scorefile_ = true;
		// set up all the information for the scorefile
		utility::file::FileName outfile("score");
		std::ostringstream oss;
		//prefix, suffix
		oss << option[ out::prefix ]() << outfile.base()
		    << option[ out::suffix ]();
		outfile.base( oss.str() );
		//path
		if ( option[ out::path::score ].user() ) {
			outfile.path( option[ out::path::score ]().path() );
		  outfile.vol( option[ out::path::score ]().vol() );
		} else outfile.path( default_path.path() );
		// determine the extension based on fullatom or centroid; this logic is usually wrong but .sc is okay anyway
		if( option[ out::file::fullatom ] ) {
			outfile.ext( ".fasc" );
		} else {
			outfile.ext( ".sc" );
		}
		scorefile_name_ = outfile.name();
	}
}

FileJobOutputter::~FileJobOutputter(){}

void FileJobOutputter::scorefile(
	JobCOP job,
	core::pose::Pose const & pose,
	std::string prefix_tag,
	std::string suffix_tag,
	std::string scorefile
)
{
	TR.Debug << "FileJobOutputter scorefile" << std::endl;
	if (!write_scorefile_) return;

	core::io::raw_data::ScoreFileData sfd((scorefile.empty() ? scorefile_name_.name() : scorefile));
	std::map < std::string, core::Real > score_map;
  std::map < std::string, std::string > string_map;
	ScoreMap::score_map_from_scored_pose( score_map, pose );

	// Adds StringReal job info into the score map for output in the scorefile.
	for( Job::StringRealPairs::const_iterator it(job->output_string_real_pairs_begin()), end(job->output_string_real_pairs_end());
			 it != end;
			 ++it) {
		score_map[it->first] = it->second;
	}

	// Adds StringString job info into a map for output in the scorefile.
	for( Job::StringStringPairs::const_iterator it(job->output_string_string_pairs_begin()), end(job->output_string_string_pairs_end());
			it != end;
			++it) {
		string_map[it->first] = it->second;
	}

	/// Extra Scores from Pose DataCache

	// ARBITRARY_FLOAT_DATA
	if( pose.data().has( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) ) {
		basic::datacache::CacheableStringFloatMapCOP data
			= dynamic_cast< basic::datacache::CacheableStringFloatMap const * >
				( pose.data().get_raw_const_ptr( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) );
		assert( data.get() != NULL );

		for(std::map< std::string, float >::const_iterator it( data->map().begin() ), end( data->map().end() );
		  it != end;
		  ++it) {
			score_map[it->first] = it->second;
		}
	}

	// ARBITRARY_STRING_DATA
	if( pose.data().has( core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA ) ) {
		basic::datacache::CacheableStringMapCOP data
			= dynamic_cast< basic::datacache::CacheableStringMap const * >
				( pose.data().get_raw_const_ptr( core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA ) );
		assert( data.get() != NULL );

		for(std::map< std::string, std::string >::const_iterator it( data->map().begin() ), end( data->map().end() );
		  it != end;
		  ++it) {
			string_map[it->first] = it->second;
		}
	}

	/*
	// Don't write score file entries that only contain one column: total_score = 0
	// 20140504: Disabled by request

	if(score_map.empty()) {
		// No scores at all
		return;
	}
	if(score_map.size() == 1) {
		std::map < std::string, core::Real >::iterator it( score_map.begin() );
		if(it != score_map.end() &&
		  it->first == name_from_score_type(core::scoring::total_score) &&
		  it->second == 0) {
			// Only one score term, and it's the total_score with value 0
			return;
		}
	}
	*/

	sfd.write_pose( pose, score_map, (prefix_tag + output_name(job) + suffix_tag), string_map );
}


///@details this base class implementation will try to _append_ whatever string it gets to a file named after the job, with a user-specified suffix.
void FileJobOutputter::file(
				  JobCOP job,
				  std::string const & data ){
  TR.Debug << "FileJobOutputter::file" << std::endl;
  if (data.empty()) return; //who needs empty files?

  //this almost certainly needs some sort of error handling - what do you think Andrew?
  utility::io::ozstream out;
  std::string const & miscfile_ext(option[ run::jobdist_miscfile_ext ].value());
  out.open_append( output_name(job) + (miscfile_ext[0] == '.' ? "" : ".") + miscfile_ext );
  out << data << std::flush;
  out.close();
}


//void other_scorefile( std::string const & tag, core::pose::Pose const & pose, std::string const & o_scorefile );

}//jd2
}//protocols
