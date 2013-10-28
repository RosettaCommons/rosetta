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
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <protocols/jd2/parser/ResidueSelectorLoader.hh>
#include <protocols/jd2/parser/ResidueSelectorLoaderCreator.hh>

// Project headers
#include <core/pack/task/residue_selector/ResidueSelector.hh>
#include <core/pack/task/residue_selector/ResidueSelectorFactory.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/mistakes.OptionKeys.gen.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>

namespace protocols {
namespace jd2 {
namespace parser {

static basic::Tracer TR( "protocols.jd2.parser.ResidueSelectorLoader" );

ResidueSelectorLoader::ResidueSelectorLoader() {}
ResidueSelectorLoader::~ResidueSelectorLoader() {}

void ResidueSelectorLoader::load_data(
	core::pose::Pose const &,
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap & datamap
) const
{
	using namespace utility::tag;
	using core::pack::task::residue_selector::ResidueSelectorOP;
	typedef utility::vector0< TagCOP > TagCOPs;

	TagCOPs const & selector_tags( tag->getTags() );
	for ( core::Size ii = 0; ii < selector_tags.size(); ++ii ) {
		TagCOP ii_tag = selector_tags[ ii ];
		ResidueSelectorOP selector = core::pack::task::residue_selector::ResidueSelectorFactory::get_instance()->new_residue_selector(
			ii_tag->getName(),
			ii_tag,
			datamap
		);

		bool const data_add_status = datamap.add( "ResidueSelector" , ii_tag->getName(), selector );
		if( !data_add_status ) {
			utility_exit_with_message( "ResidueSelector '" + ii_tag->getName() + "' already exists in the basic::datacache::DataMap. Please rename." );
		}
	}
	TR.flush();
}

DataLoaderOP
ResidueSelectorLoaderCreator::create_loader() const { return new ResidueSelectorLoader; }

std::string
ResidueSelectorLoaderCreator::keyname() const { return "RESIDUE_SELECTORS"; }


} //namespace parser
} //namespace jd2
} //namespace protocols
