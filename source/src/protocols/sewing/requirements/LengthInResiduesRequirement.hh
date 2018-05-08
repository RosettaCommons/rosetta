// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/requirements/LengthInResiduesRequirement.hh
/// @brief a Requirement that an Assembly be within a certain range of lengths
/// @author frankdt (frankdt@email.unc.edu)


#ifndef INCLUDED_protocols_sewing_requirements_LengthInResiduesRequirement_hh
#define INCLUDED_protocols_sewing_requirements_LengthInResiduesRequirement_hh

#include <protocols/sewing/requirements/LengthInResiduesRequirement.fwd.hh>
#include <protocols/sewing/data_storage/SmartAssembly.hh>
#include <core/types.hh>
// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace sewing {
namespace requirements {

///@brief a Requirement that an Assembly be within a certain range of lengths
class LengthInResiduesRequirement : public AssemblyRequirement {

public:

	LengthInResiduesRequirement();
	LengthInResiduesRequirement(core::Size,core::Size);
	LengthInResiduesRequirement(LengthInResiduesRequirement const & src);

	virtual ~LengthInResiduesRequirement() override;

	LengthInResiduesRequirementOP
	clone() const;

	std::pair<bool,bool>
	test (data_storage::SmartAssemblyOP assembly) override;

	std::string
	get_name() override {
		return type_name();
	}

	core::Size
	get_minimum_length() const;

	core::Size
	get_maximum_length() const;

	void
	set_minimum_length( core::Size );

	void
	set_maximum_length( core::Size );



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
	core::Size length_;
	core::Size minimum_length_;
	core::Size maximum_length_;
	std::pair<bool,bool> test_results_;
};


} //protocols
} //sewing
} //requirements



#endif //INCLUDED_protocols_sewing_requirements_LengthInResiduesRequirement_hh





