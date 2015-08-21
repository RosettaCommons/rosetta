// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/metrics/CalculatorFactory.cc
/// @brief  CalculatorFactory class
/// @author John Karanicolas

// Unit headers
#include <core/pose/metrics/CalculatorFactory.hh>

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>


// Utility headers
#include <basic/Tracer.hh>
#include <utility/exit.hh>

// C++ Headers
#include <map>
#include <iostream>

#include <utility/vector1.hh>


namespace core {
namespace pose {
namespace metrics {

CalculatorFactory& CalculatorFactory::Instance() {
	static CalculatorFactory singleton;
	return singleton;
}


void CalculatorFactory::register_calculator( std::string const & calculator_name, PoseMetricCalculatorOP const new_calculator ) {
	if ( check_calculator_exists( calculator_name ) ) {
		basic::Error() << "Cannot register a calculator with name: " << calculator_name << std::endl;
		basic::Error() << "This calculator already exists" << std::endl;
		utility_exit();
	}
	calculators_.insert( std::make_pair ( calculator_name, new_calculator ) );
	return;
}


bool CalculatorFactory::check_calculator_exists( std::string const & calculator_name ) {
	std::map< std::string, PoseMetricCalculatorOP >::const_iterator calculator_iter;
	calculator_iter = calculators_.find( calculator_name );
	if ( calculator_iter == calculators_.end() ) {
		// this calculator has not yet been setup
		return false;
	}
	return true;
}


/// @brief remove a calculator from the factory
/// @return true if calculator removed, false if no such calculator
bool CalculatorFactory::remove_calculator( std::string const & calculator_name ) {
	std::map< std::string, PoseMetricCalculatorOP >::iterator i = calculators_.find( calculator_name );

	if ( i != calculators_.end() ) {
		calculators_.erase( i );
		return true;
	}

	return false;
}


/// @brief clear all calculators from factory
/// @return false if no calculators in list, true otherwise
bool CalculatorFactory::clear_calculators() {
	if ( calculators_.empty() ) {
		return false;
	}

	calculators_.clear();
	return true;
}


PoseMetricCalculatorOP CalculatorFactory::retrieve_calculator( std::string const & calculator_name ) {
	std::map< std::string, PoseMetricCalculatorOP >::const_iterator calculator_iter;
	calculator_iter = calculators_.find( calculator_name );
	if ( calculator_iter == calculators_.end() ) {
		// this calculator has not yet been setup
		basic::Error() << "Could not find calculator " << calculator_name << " - need to register it before use" << std::endl;
		utility_exit();
	}
	return calculator_iter->second->clone();
}


} // metrics
} // pose
} // core
