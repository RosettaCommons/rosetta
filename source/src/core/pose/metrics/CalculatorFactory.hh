// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author John Karanicolas


#ifndef INCLUDED_core_pose_metrics_CalculatorFactory_hh
#define INCLUDED_core_pose_metrics_CalculatorFactory_hh

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>


// C++ Headers
#include <map>

#include <utility/vector1.hh>


namespace core {
namespace pose {
namespace metrics {

// Note: CalculatorFactory is a singleton class
class CalculatorFactory {
public:

	static CalculatorFactory& Instance();

	void register_calculator( std::string const & calculator_name, PoseMetricCalculatorOP const new_calculator );

	bool check_calculator_exists( std::string const & calculator_name );

	/// @brief remove a calculator from the factory
	/// @return true if calculator removed, false if no such calculator
	bool remove_calculator( std::string const & calculator_name );

	/// @brief clear all calculators from factory
	/// @return false if no calculators in list, true otherwise
	bool clear_calculators();

	PoseMetricCalculatorOP retrieve_calculator( std::string const & calculator_name );

private:

	CalculatorFactory() {};
	CalculatorFactory( CalculatorFactory const & src );
	CalculatorFactory const & operator=( CalculatorFactory const & src );

	std::map< std::string, PoseMetricCalculatorOP > calculators_;

};


} // namespace metrics
} // namespace pose
} // namespace core


#endif
