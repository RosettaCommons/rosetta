// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/parser/LoopsDefinerLoader.cc
/// @brief  Implementation the LoopsDefinerLoader class which implements the DataLoader interface
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <protocols/jd2/parser/LoopsDefinerLoader.hh>
#include <protocols/jd2/parser/StandardLoaderCreators.hh>

// Project Headers
#include <protocols/loops/loops_definers/LoopsDefiner.hh>
#include <protocols/loops/loops_definers/LoopsDefinerFactory.hh>
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <protocols/moves/DataMap.hh>


// Utility headers
#include <utility/tag/Tag.hh>

// Boost Headers
#include <boost/foreach.hpp>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

#define foreach BOOST_FOREACH

using std::string;
using std::endl;
using core::pose::Pose;
using utility::tag::TagPtr;
using protocols::moves::DataMap;
using utility::vector0;
using protocols::loops::loops_definers::LoopsDefinerOP;
using protocols::loops::loops_definers::LoopsDefinerFactory;

namespace protocols {
namespace jd2 {
namespace parser {

static basic::Tracer TR( "protocols.jd2.parser.LoopsDefinerLoader" );

LoopsDefinerLoader::LoopsDefinerLoader() {}
LoopsDefinerLoader::~LoopsDefinerLoader() {}

void LoopsDefinerLoader::load_data(
	Pose const & pose,
	TagPtr const tag,
	DataMap & data
) const
{
	typedef vector0< TagPtr > TagPtrs;

	foreach(TagPtr tag, tag->getTags()){
		string const type( tag->getName() );
		if ( ! tag->hasOption("name") ) {
			utility_exit_with_message( "Can't create unnamed Loops definition (type: " + type + ")" );
		}
		string const name( tag->getOption<string>("name") );
		if ( data.has( "loops_definers", name ) ) {
			TR.Error << "Error LoopsDefiner of name \"" << name
				<< "\" (with type " << type << ") already exists. \n" << tag << endl;
			utility_exit_with_message("Duplicate definition of LoopsDefiner with name " + name);
		}
		LoopsDefinerOP loops_definer( LoopsDefinerFactory::get_instance()->create_loops_definer( type ) );
		loops_definer->parse_my_tag(tag, data, pose);
		data.add("loops_definers", name, loops_definer );
		TR << "Created LoopsDefiner named \"" << name << "\" of type " << type << endl;
	}
	TR.flush();
}

DataLoaderOP
LoopsDefinerLoaderCreator::create_loader() const { return new LoopsDefinerLoader; }

string
LoopsDefinerLoaderCreator::keyname() const { return "LOOP_DEFINITIONS"; }


} //namespace parser
} //namespace jd2
} //namespace protocols
