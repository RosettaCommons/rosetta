// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Zhe Zhang

#ifndef INCLUDED_devel_replica_docking_WrapFilterAsEvaluatorCreator_hh
#define INCLUDED_devel_replica_docking_WrapFilterAsEvaluatorCreator_hh

// Unit Headers
#include <protocols/evaluation/EvaluatorCreator.hh>

#include <protocols/evaluation/PoseEvaluator.hh>

#include <core/types.hh>
#include <utility/vector1.hh>


namespace devel {
namespace replica_docking {

/// @brief creator for the ChemicalShiftsEvaluatorCreator class
class WrapFilterAsEvaluatorCreator : public protocols::evaluation::EvaluatorCreator
{
public:
	WrapFilterAsEvaluatorCreator() : options_registered_(false) {};
	~WrapFilterAsEvaluatorCreator() override;

	virtual void register_options();

	void add_evaluators( protocols::evaluation::MetaPoseEvaluator & eval ) const override;

	std::string type_name() const override;

private:
	bool options_registered_;
};

} //namespace
} //namespace

#endif
