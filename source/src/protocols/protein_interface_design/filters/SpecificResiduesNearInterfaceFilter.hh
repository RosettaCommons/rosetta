// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/SpecificResiduesNearInterfaceFilter.hh
/// @brief Reports the average degree of connectivity of interface residues
/// @author Sarel Fleishman (sarelf@uw.edu)

#ifndef INCLUDED_protocols_protein_interface_design_filters_SpecificResiduesNearInterfaceFilter_hh
#define INCLUDED_protocols_protein_interface_design_filters_SpecificResiduesNearInterfaceFilter_hh


// Project Headers
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/protein_interface_design/filters/SpecificResiduesNearInterfaceFilter.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

#include <utility/vector1.hh>

// Unit headers

namespace protocols {
namespace protein_interface_design {
namespace filters {

class SpecificResiduesNearInterfaceFilter : public protocols::filters::Filter
{
private:
	typedef protocols::filters::Filter parent;
public:
	/// @brief default ctor
	SpecificResiduesNearInterfaceFilter();
	SpecificResiduesNearInterfaceFilter(
		SpecificResiduesNearInterfaceFilter const & src);
	~SpecificResiduesNearInterfaceFilter();

	core::pack::task::TaskFactoryOP task_factory() const;
	void task_factory( core::pack::task::TaskFactoryOP tf );

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );


	/// @brief Constructor with a single target residue
	core::Size compute( core::pose::Pose const & pose ) const;
	bool apply( core::pose::Pose const & pose ) const;
	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	protocols::filters::FilterOP clone() const;
	protocols::filters::FilterOP fresh_instance() const;


private:

	core::pack::task::TaskFactoryOP task_factory_;
	core::Size rb_jump_;
	// TODO make this work
	//core::Real interface_distance_threshold_;
};

} // filters
} //protein_interface_design
} // protocols

#endif //INCLUDED_protocols_Filters_SpecificResiduesNearInterfaceFilter_HH_
