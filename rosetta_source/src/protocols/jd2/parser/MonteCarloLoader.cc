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
#include <protocols/jd2/parser/MonteCarloLoader.hh>
#include <protocols/jd2/parser/StandardLoaderCreators.hh>

// Project Headers
#include <core/scoring/ScoreFunction.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/tag/Tag.hh>


namespace protocols {
namespace jd2 {
namespace parser {

static basic::Tracer TR( "protocols.jd2.parser.MonteCarloLoader" );

MonteCarloLoader::MonteCarloLoader() {}
MonteCarloLoader::~MonteCarloLoader() {}

void MonteCarloLoader::load_data(
	core::pose::Pose const &,
	utility::tag::TagPtr const tag,
	moves::DataMap & data
) const
{
	using namespace utility::tag;
	typedef utility::vector0< TagPtr > TagPtrs;

	TagPtrs const montecarlo_tags( tag->getTags() );

	for( TagPtrs::const_iterator montecarlo_tag_it=montecarlo_tags.begin(); montecarlo_tag_it!=montecarlo_tags.end();
			++montecarlo_tag_it ) {
		TagPtr const montecarlo_tag = *montecarlo_tag_it;

		std::string const mc_name( montecarlo_tag->getName() );
		core::Real const mctemp( montecarlo_tag->getOption< core::Real >( "temperature", 2.0 ));
		std::string const sfxn_name( montecarlo_tag->getOption< std::string > ( "scorefunction", "score12" ));
		core::scoring::ScoreFunctionOP scorefxn = new core::scoring::ScoreFunction(
			*data.get< core::scoring::ScoreFunction * >( "scorefxns", sfxn_name ));

		protocols::moves::MonteCarloOP mc = new protocols::moves::MonteCarlo( *scorefxn, mctemp );
		// add more options for the MonteCarlo object here, e.g.
		// 1. autotemp / quenchtemp
		// 2. heat after cycles
		data.add( "montecarlos" , mc_name, mc );
	}
	TR.flush();
}

DataLoaderOP
MonteCarloLoaderCreator::create_loader() const { return new MonteCarloLoader; }

std::string
MonteCarloLoaderCreator::keyname() const { return "MONTECARLOS"; }


} //namespace parser
} //namespace jd2
} //namespace protocols
