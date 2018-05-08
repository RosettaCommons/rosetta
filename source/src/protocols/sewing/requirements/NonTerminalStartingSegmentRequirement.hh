// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/requirements/NonTerminalStartingSegmentRequirement.hh
/// @brief a Requirement that an Assembly have less than a certain number of clashes
/// @author frankdt (frankdt@email.unc.edu)


#ifndef INCLUDED_protocols_sewing_requirements_NonTerminalStartingSegmentRequirement_hh
#define INCLUDED_protocols_sewing_requirements_NonTerminalStartingSegmentRequirement_hh

#include <protocols/sewing/requirements/NonTerminalStartingSegmentRequirement.fwd.hh>
#include <protocols/sewing/requirements/AssemblyRequirement.hh>
#include <protocols/sewing/data_storage/SmartAssembly.hh>
#include <protocols/sewing/data_storage/SmartSegment.hh>
#include <protocols/sewing/data_storage/SmartSewingResidue.hh>
// Utility headers
#include <core/types.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace sewing {
namespace requirements {

///@brief a Requirement that an Assembly have less than a certain number of clashes
class NonTerminalStartingSegmentRequirement : public AssemblyRequirement {

public:

	NonTerminalStartingSegmentRequirement();
	NonTerminalStartingSegmentRequirement(NonTerminalStartingSegmentRequirement const & src);

	virtual ~NonTerminalStartingSegmentRequirement() override;

	NonTerminalStartingSegmentRequirementOP
	clone() const;

	std::pair<bool,bool>
	test(data_storage::SmartAssemblyOP assembly) override;

	std::string
	get_name() override {
		return type_name();
	}
	void
	set_options_from_tag(
		utility::tag::TagCOP requirement_tag,
		basic::datacache::DataMap& datamap,
		protocols::filters::Filters_map const & filtermap,
		protocols::moves::Movers_map const & movermap,
		core::pose::Pose const & pose) override;

	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & );
	static std::string type_name();
private:
	std::pair<bool,bool> test_results_;
	data_storage::SmartSegmentOP active_segment_;
};


} //protocols
} //sewing
} //requirements



#endif //INCLUDED_protocols_sewing_requirements_NonTerminalStartingSegmentRequirement_hh





