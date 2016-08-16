// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_filters/PredictedBurialFnEvaluatorCreator.hh
/// @brief  Header for PredictedBurialFnEvaluatorCreator
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_simple_filters_PredictedBurialFnEvaluatorCreator_hh
#define INCLUDED_protocols_simple_filters_PredictedBurialFnEvaluatorCreator_hh

// Unit Headers
#include <protocols/evaluation/EvaluatorCreator.hh>

#include <protocols/evaluation/PoseEvaluator.hh>

#include <core/types.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace simple_filters {

/// @brief creator for the PredictedBurialFnEvaluatorCreator class
class PredictedBurialFnEvaluatorCreator : public evaluation::EvaluatorCreator
{
public:
	PredictedBurialFnEvaluatorCreator() : options_registered_(false) {};
	virtual ~PredictedBurialFnEvaluatorCreator();

	virtual void register_options();

	virtual void add_evaluators( evaluation::MetaPoseEvaluator & eval ) const;

	virtual std::string type_name() const;

private:
	bool options_registered_;
};

} //namespace
} //namespace

#endif
