// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/parser/DataLoader.cc
/// @brief  Implementation of the XML parser's DataLoader base class (ctor & dstor)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <protocols/jd2/parser/MonteCarloLoader.hh>
#include <protocols/jd2/parser/StandardLoaderCreators.hh>

// Project Headers
#include <core/scoring/ScoreFunction.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/tag/Tag.hh>

// Boost Headers
#include <boost/foreach.hpp>

#include <basic/datacache/DataMap.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace jd2 {
namespace parser {

static THREAD_LOCAL basic::Tracer TR( "protocols.jd2.parser.MonteCarloLoader" );

MonteCarloLoader::MonteCarloLoader() {}
MonteCarloLoader::~MonteCarloLoader() {}

void MonteCarloLoader::load_data(
	core::pose::Pose const &,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
) const
{
	using namespace utility::tag;
	typedef utility::vector0< TagCOP > TagCOPs;

	TagCOPs const montecarlo_tags( tag->getTags() );

	BOOST_FOREACH ( TagCOP montecarlo_tag, montecarlo_tags ) {
		std::string const mc_name( montecarlo_tag->getName() );
		core::Real const mctemp( montecarlo_tag->getOption< core::Real >( "temperature", 2.0 ));
		core::scoring::ScoreFunctionOP scorefxn =
			rosetta_scripts::parse_score_function( montecarlo_tag, data )->clone();

		protocols::moves::MonteCarloOP mc( new protocols::moves::MonteCarlo( *scorefxn, mctemp ) );
		// add more options for the MonteCarlo object here, e.g.
		// 1. autotemp / quenchtemp
		// 2. heat after cycles
		data.add( "montecarlos" , mc_name, mc );
	}
	TR.flush();
}

DataLoaderOP
MonteCarloLoaderCreator::create_loader() const { return DataLoaderOP( new MonteCarloLoader ); }

std::string
MonteCarloLoaderCreator::keyname() const { return "MONTECARLOS"; }


} //namespace parser
} //namespace jd2
} //namespace protocols
