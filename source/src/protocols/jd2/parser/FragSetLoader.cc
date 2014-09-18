// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/parser/DataLoader.cc
/// @brief  Implementation of the XML parser's DataLoader base class (ctor & dstor)
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit Headers
#include <protocols/jd2/parser/FragSetLoader.hh>
#include <protocols/jd2/parser/StandardLoaderCreators.hh>

// for fragments
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/OrderedFragSet.hh>
#include <protocols/jd2/parser/FragmentReader.hh>

#include <basic/Tracer.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

// Basic headers

// Boost Headers
#include <boost/foreach.hpp>

#include <basic/datacache/DataMap.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

static thread_local basic::Tracer TR( "protocols.jd2.parser.FragSetLoader" );

namespace protocols {
namespace jd2 {
namespace parser {

FragSetLoader::FragSetLoader() {}
FragSetLoader::~FragSetLoader() {}

void FragSetLoader::load_data(
	core::pose::Pose const &,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
) const
{
	using namespace utility::tag;

	using protocols::jd2::parser::FragmentReader;
	using protocols::jd2::parser::FragmentReaderOP;
	typedef std::map< std::string, FragmentReaderOP > FragmentReaderMap;

	FragmentReaderMap frag_readers_map;
	if ( tag->hasTag( "FRAGMENTS" ) ) {
		BOOST_FOREACH(TagCOP tag, tag->getTag( "FRAGMENTS" )->getTags()){
			std::string const name ( tag->getName() ); // this name is used when fragsets are defined later.
			runtime_assert( !name.empty() );
			FragmentReaderOP frop = new FragmentReader( tag );
			frag_readers_map[ name ] = frop;
		}
	} else {
		TR << "No tag of FRAGMENTS" << std::endl;
		runtime_assert( false );
	}

	BOOST_FOREACH( TagCOP tag, tag->getTags() ){
		std::string const name ( tag->getName() );
		if( name == "FRAGMENTS" ) continue;

		std::string const frag_name ( tag->getOption<std::string>( "frag_name", "" ) );
		std::string const output ( tag->getOption<std::string>( "output", "" ) );
		runtime_assert( !name.empty() && frag_name != "" );

		core::fragment::FragSetOP fragset = new core::fragment::OrderedFragSet;
		utility::vector1< std::string > fnames ( utility::string_split( frag_name, ',' ) );
		BOOST_FOREACH(std::string fname, fnames){
			std::map< std::string, FragmentReaderOP >::const_iterator itr;
			itr = frag_readers_map.find( fname );
			if ( itr != frag_readers_map.end() ){
				FragmentReaderOP frop ( frag_readers_map[ fname ] );
				frop->apply( fragset );
			}else{
				TR << "frag_name " << fname << " does not exist." << std::endl;
				runtime_assert( false );
			}
		}
		runtime_assert( fragset->nr_frames() != 0 );
		data.add( "fragsets", name,  fragset );
		// output flagments to fyile
		if( !output.empty() ){
			core::fragment::FragmentIO().write_data( output, *fragset );
		}
	}

}

DataLoaderOP
FragSetLoaderCreator::create_loader() const { return new FragSetLoader; }

std::string
FragSetLoaderCreator::keyname() const { return "FRAGSETS"; }


} //namespace parser
} //namespace jd2
} //namespace protocols
