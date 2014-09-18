// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/outputter/FormatStringOutputter.cc
/// @brief An abstract base class that implements format string parsing
/// @author Ken Jung

// Unit Headers
#include <protocols/outputter/FormatStringOutputter.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <utility/file/FileName.hh>

// tracer
#include <basic/Tracer.hh>

namespace protocols {
namespace outputter {

static thread_local basic::Tracer TR( "protocols.outputter.FormatStringOutputter" );

#ifdef USELUA
void lregister_FormatStringOutputter( lua_State * lstate ) {
	lregister_Outputter(lstate);

	luabind::module(lstate, "protocols")
	[
		luabind::namespace_("outputter")
		[
			luabind::class_<FormatStringOutputter, Outputter>("FormatStringOutputter")
						.def("format_string", (void (FormatStringOutputter::*) (std::string) ) &FormatStringOutputter::format_string)
		]
	];
}
#endif

FormatStringOutputter::FormatStringOutputter(){}
FormatStringOutputter::~FormatStringOutputter(){}

OutputterSP FormatStringOutputter::create() {
	// never should get here
	return OutputterSP();
}

void FormatStringOutputter::write( PipeMap & p ){
	TR << "-------Outputting PipeMap to File--------" << std::endl;
	for( PipeMap::iterator itr = p.begin(); itr != p.end(); itr++ ) {
		PipeSP current_pipe = itr->second;
		filenameparts_["pipe_name"] = itr->first;
		if( current_pipe ) {
			write( *current_pipe );
		}
	}
}

void FormatStringOutputter::write( Pipe & p ) {
	for( core::Size idx = 0; idx < p.size(); idx++ ) {
		PoseSP current_pose = p[idx];
		filenameparts_["pipe_idx"] = idx;
		if( current_pose ) {
			std::string filename;
			core::pose::get_comment(*current_pose, "inputfile", filename );
			utility::file::FileName f(filename);
			filenameparts_["filepath"] = f.path();
			filenameparts_["filename"] = f.base() + f.extension();
			filenameparts_["filebasename"] = f.base();
			core::pose::get_comment(*current_pose, "filemultiplier", filenameparts_["filemultiplier"] );

			write( *current_pose );
		}
	}
}

// walk through format string char by char, look for key matches to filenameparts and replace as needed
void FormatStringOutputter::parse_format_string( boost::unordered_map< std::string, std::string> & filenameparts, std::string const & format_string, std::string & filename ) {
	std::string potential_key = "";
	bool in_potential_key = false;
	for( core::Size i = 0; i < format_string.size() ; i++ ) {
		if( format_string[i] != '%' && ! in_potential_key ) {
			filename.push_back( format_string[i] );
		} else if ( format_string[i] == '%' && in_potential_key ) {
			// the key we though wasn't, so push it onto the filename
			filename.push_back('%');
			filename += potential_key;
			potential_key = "";
		} else if ( format_string[i] == '%' ) {
			in_potential_key = true;
		} else if ( in_potential_key ) {
			potential_key.push_back(format_string[i]);
			boost::unordered_map< std::string, std::string>::iterator itr = filenameparts.find( potential_key );
			if( itr != filenameparts.end() ) {
				// potential key is real key
				filename += itr->second;
				in_potential_key=false;
				potential_key="";
			}
		}
	}
	if( in_potential_key ) 
		// the rest of the input name was a potential key that turned out to be false, so just append it
		filename += '%' + potential_key;
}


#ifdef USELUA
void FormatStringOutputter::parse_def( utility::lua::LuaObject const & def,
		utility::lua::LuaObject const & tasks ) {
	format_string_ = def["format_string"] ? def["format_string"].to<std::string>() : "%filebasename_%filemultiplier.pdb";
}

void FormatStringOutputter::lregister( lua_State * lstate ) {
	lregister_FormatStringOutputter(lstate);
}
#endif

} // outputter
} // protocols
