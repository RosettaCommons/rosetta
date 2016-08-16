// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/feature/EvaluatorCreator.hh
/// @brief  Base class for EvaluatorCreators for the Evaluator load-time factory registration scheme
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_evaluation_EvaluatorCreator_hh
#define INCLUDED_protocols_evaluation_EvaluatorCreator_hh

// Unit Headers
#include <protocols/evaluation/EvaluatorCreator.fwd.hh>

// Package Headers
#include <protocols/evaluation/PoseEvaluator.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// c++ headers
#include <string>

namespace protocols {
namespace evaluation {

/// @brief The Creator class is responsible for creating a particular
/// mover class.
class EvaluatorCreator : public utility::pointer::ReferenceCount
{
public:
	EvaluatorCreator() {}
	virtual ~EvaluatorCreator() {}

	virtual void add_evaluators( MetaPoseEvaluator & ) const = 0;
	virtual std::string type_name() const = 0;

};

} //namespace
} //namespace

#endif
