// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jobdist/Jobs.cc
///
/// @brief
/// @author Ian W. Davis


#include <basic/options/option.hh>
#include <protocols/jobdist/Jobs.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/string_util.hh>

#include <iomanip>
#include <sstream>


// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace jobdist {


/// @details Deliberately discards any path information in the input tag
/// as well as any file name extension (since input tags are usually file names).
/// There is some possibility this could lead to non-unique output tags,
/// which deserves further consideration at some point...
std::string BasicJob::output_tag(int struct_n) const
{
	using basic::options::option;
	using namespace basic::options::OptionKeys;

	// Use at least 4 digits in number to match Rosetta++
	int nstruct_width = 0;
	for(int i = 1; i <= nstruct_ || nstruct_width < 4; i *= 10) nstruct_width += 1;
	// Treat tags as file names so that we put the number before the extension.
	// Everything will still work if they're not file names, though.
	utility::vector1<std::string> temp_out_names= utility::split(input_id_);
	utility::file::FileName out_name = utility::file::combine_names( temp_out_names);
//jobs_tracer<< out_name.base()<< std::endl;
	if( option[ run::shuffle ].user() ) out_name = "S_shuffle";
	std::ostringstream oss;

	std::string user_tag("");
	if ( basic::options::option[ basic::options::OptionKeys::out::user_tag ].user() ) {
		user_tag = "_" + basic::options::option[ basic::options::OptionKeys::out::user_tag ];
	}

	if( !preserve_whole_input_tag_ ){
		oss << option[ out::prefix ]() << out_name.base() << option[ out::suffix ]()
				<< user_tag << '_' << std::setfill('0') << std::setw(nstruct_width) << (struct_n);
	}else{
		oss << option[ out::prefix ]() << out_name << option[ out::suffix ]()
		    << user_tag << '_' << std::setfill('0') << std::setw(nstruct_width) << (struct_n);
	}
	out_name.base( oss.str() );
	return out_name.base();
}


} // namespace jobdist
} // namespace protocols
