// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/requirements/SizeInSegmentsRequirement.hh
/// @brief a Requirement that an Assembly be within a certain range of sizes
/// @author frankdt (frankdt@email.unc.edu)


#ifndef INCLUDED_protocols_sewing_requirements_SizeInSegmentsRequirement_hh
#define INCLUDED_protocols_sewing_requirements_SizeInSegmentsRequirement_hh

#include <protocols/sewing/requirements/SizeInSegmentsRequirement.fwd.hh>
#include <protocols/sewing/data_storage/SmartAssembly.hh>
#include <core/types.hh>
// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace sewing {
namespace requirements {

///@brief a Requirement that an Assembly be within a certain range of sizes
class SizeInSegmentsRequirement : public AssemblyRequirement {

public:

	SizeInSegmentsRequirement();
	SizeInSegmentsRequirement(core::Size,core::Size);
	SizeInSegmentsRequirement(SizeInSegmentsRequirement const & src);

	virtual ~SizeInSegmentsRequirement() override;

	SizeInSegmentsRequirementOP
	clone() const;

	std::pair<bool,bool>
	test (data_storage::SmartAssemblyOP assembly) override;

	std::string
	get_name() override {
		return type_name();
	}

	//Getters and Setters
	core::Size
	get_minimum_size() const;

	core::Size
	get_maximum_size() const;

	void
	set_minimum_size( core::Size );

	void
	set_maximum_size( core::Size );




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
	core::Size size_;
	core::Size minimum_size_;
	core::Size maximum_size_;
	std::pair<bool,bool> test_results_;
};


} //protocols
} //sewing
} //requirements



#endif //INCLUDED_protocols_sewing_requirements_SizeInSegmentsRequirement_hh





