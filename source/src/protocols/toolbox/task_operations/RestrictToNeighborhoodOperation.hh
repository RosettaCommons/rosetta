// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/RestrictToNeighborhoodOperation.hh
/// @brief  TaskOperation class that finds a neighborhood and makes it mobile in the PackerTask
/// @author Steven Lewis smlewi@unc.edu

#ifndef INCLUDED_protocols_toolbox_task_operations_RestrictToNeighborhoodOperation_hh
#define INCLUDED_protocols_toolbox_task_operations_RestrictToNeighborhoodOperation_hh

// Unit Headers
#include <protocols/toolbox/task_operations/RestrictToNeighborhoodOperation.fwd.hh>
#include <protocols/toolbox/task_operations/RestrictOperationsBase.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>

#include <protocols/toolbox/pose_metric_calculators/NeighborhoodByDistanceCalculator.fwd.hh>

// Utility Headers
#include <core/types.hh>

// C++ Headers
#include <string>
#include <set>

#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>


namespace protocols {
namespace toolbox {
namespace task_operations {

/// @details this class is a TaskOperation to prevent repacking of residues not near a neighborhood.  Internally it just user NeighborhoodByDistanceCalculator to do all the work.  You are allowed a "calculator name" interface to tell it which you-constructed calculator to use, or you can give it a set of residues and a desired distance cutoff to have it define the calculator itself.  (There is no need to use both interfaces).
class RestrictToNeighborhoodOperation : public RestrictOperationsBase
{
public:
	typedef RestrictOperationsBase parent;
	typedef std::set< core::Size > SizeSet;

	RestrictToNeighborhoodOperation();

	RestrictToNeighborhoodOperation( RestrictToNeighborhoodOperation const & rhs );

	RestrictToNeighborhoodOperation( std::set< core::Size > const & central_residues, core::Real const dist_cutoff );

	/// @brief uses option system default for dist_cutoff
	RestrictToNeighborhoodOperation( std::set< core::Size > const & central_residues );

	RestrictToNeighborhoodOperation( std::string const & calculator );

	virtual ~RestrictToNeighborhoodOperation();

	/// @brief assignment operator
	RestrictToNeighborhoodOperation & operator=( RestrictToNeighborhoodOperation const & rhs );

	virtual core::pack::task::operation::TaskOperationOP clone() const;

	virtual
	void
	apply( core::pose::Pose const &, core::pack::task::PackerTask & ) const;
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "RestrictToNeighborhood"; }

	//getters
	/// @brief this nontrivially checks the underlying calculator
	SizeSet const & get_central_residues() const;
	/// @brief this nontrivially checks the underlying calculator
	core::Real get_distance_cutoff() const;
	/// @brief trivially returns underlying calculator string name
	std::string const & get_calculator_name() const { return calculator_name_; };
	/// @brief look up actual calculator object
	protocols::toolbox::pose_metric_calculators::NeighborhoodByDistanceCalculatorCOP get_calculator() const;

	//setters
	/// @brief reskin of normal make_calculator
	void set_neighborhood_parameters( SizeSet const & central_residues, core::Real dist_cutoff );

	/// @brief reskin of normal make_calculator
	void set_neighborhood_parameters( SizeSet const & central_residues );

	void set_calculator_by_name(std::string const & calculator_name) { calculator_name_ = calculator_name; };

private:
	/// @brief constructor helper function - makes the PoseMetricCalculator
	void make_calculator( std::set< core::Size > const & central_residues, core::Real dist_cutoff );

	/// @brief constructor helper function - makes the PoseMetricCalculator
	void make_calculator( std::set< core::Size > const & central_residues );

	/// @brief constructor helper function - names the PoseMetricCalculator
	void make_name( std::set< core::Size > const & central_residues );

	std::string calculator_name_;
};

} //namespace protocols
} //namespace toolbox
} //namespace task_operations

#endif // INCLUDED_protocols_toolbox_TaskOperations_RestrictToNeighborhoodOperation_HH
