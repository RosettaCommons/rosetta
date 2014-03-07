// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/evaluation/EvaluatorFactory.cc
/// @brief  Factory for creating Evaluators objects
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <protocols/evaluation/EvaluatorFactory.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/evaluation/EvaluatorCreator.hh>

// Package Headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// Project Headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector0.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/thread/threadsafe_creation.hh>

// C++ Headers
#include <sstream>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/foreach.hpp>

namespace protocols {
namespace evaluation {

using std::endl;
using std::string;
using std::pair;
using std::stringstream;
using core::pose::Pose;
using core::scoring::ScoreFunctionOP;
using protocols::filters::Filters_map;
using basic::datacache::DataMap;
using protocols::moves::Movers_map;
using utility::tag::TagCOP;

static basic::Tracer tr("protocols.evaluator.EvaluatorFactory");

EvaluatorFactory * EvaluatorFactory::instance_( 0 );


#ifdef MULTI_THREADED
#ifdef CXX11

std::mutex EvaluatorFactory::singleton_mutex_;

std::mutex & EvaluatorFactory::singleton_mutex() { return singleton_mutex_; }

#endif
#endif

/// @brief static function to get the instance of ( pointer to) this singleton class
EvaluatorFactory * EvaluatorFactory::get_instance()
{
	boost::function< EvaluatorFactory * () > creator = boost::bind( &EvaluatorFactory::create_singleton_instance );
	utility::thread::safely_create_singleton( creator, instance_ );
	return instance_;
}

EvaluatorFactory *
EvaluatorFactory::create_singleton_instance()
{
	return new EvaluatorFactory;
}


/// @details Private constructor insures correctness of singleton.
EvaluatorFactory::EvaluatorFactory() {}

EvaluatorFactory::EvaluatorFactory(
	const EvaluatorFactory &
) {}

EvaluatorFactory::~EvaluatorFactory() {}

void
EvaluatorFactory::factory_register(
	EvaluatorCreatorCOP creator
) {
	types_.push_back(pair<string,EvaluatorCreatorCOP>(creator->type_name(), creator));
}


void
EvaluatorFactory::add_evaluators(
	string const & type_name,
	MetaPoseEvaluator & eval
) {

	tr.Trace << "generate Evaluator of type " << type_name << std::endl;

	bool found(false);

	for(EvaluatorCreatorMap::const_iterator type = types_.begin(), type_end = types_.end(); type != type_end; ++type) {
		if(type->first == type_name){
			type->second->add_evaluators(eval);
			found=true;
			break;
		}
	}
	if (!found) {
		stringstream error_msg;
		error_msg
			<< "Attempting to create unrecognized Evaluator "
			<< "'" << type_name << "'." << endl
			<< "check spelling or "
			<< "register a new Evaluator in the EvaluatorFactory" << endl
			<< "known Evaluator types are:" << endl;

	for(EvaluatorCreatorMap::const_iterator type = types_.begin(), type_end = types_.end(); type != type_end; ++type) {
			error_msg << "\t" << type->first << endl;
		}
		utility_exit_with_message(error_msg.str());
	}
}

void
EvaluatorFactory::add_all_evaluators(
	MetaPoseEvaluator & eval ) {

	for(EvaluatorCreatorMap::const_iterator type = types_.begin(), type_end = types_.end(); type != type_end; ++type) {
		type->second->add_evaluators(eval);
	}
}

} // namespace
} // namespace
