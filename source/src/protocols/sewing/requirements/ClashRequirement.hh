// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/requirements/ClashRequirement.hh
/// @brief a Requirement that an Assembly have less than a certain number of clashes
/// @author frankdt (frankdt@email.unc.edu)


#ifndef INCLUDED_protocols_sewing_requirements_ClashRequirement_hh
#define INCLUDED_protocols_sewing_requirements_ClashRequirement_hh

#include <protocols/sewing/requirements/ClashRequirement.fwd.hh>
#include <protocols/sewing/requirements/AssemblyRequirement.hh>
#include <protocols/sewing/data_storage/SmartAssembly.hh>
#include <protocols/sewing/data_storage/SmartSegment.hh>
#include <protocols/sewing/data_storage/SmartSewingResidue.hh>
#include <core/conformation/Residue.hh>
// Utility headers
#include <core/types.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace sewing {
namespace requirements {

///@brief a Requirement that an Assembly have less than a certain number of clashes
class ClashRequirement : public AssemblyRequirement {

public:

	ClashRequirement();
	ClashRequirement(ClashRequirement const & src);

	virtual ~ClashRequirement() override;

	ClashRequirementOP
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

	core::Size
	get_maximum_clashes_allowed() const;

	core::Real
	get_clash_radius() const;


	void
	set_maximum_clashes_allowed( core::Size );

	void
	set_clash_radius( core::Real );



	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & );

	static utility::tag::AttributeList
	get_xml_attributes();

	static std::string
	type_name();

private:
	core::Size maximum_clashes_allowed_;
	core::Real clash_radius_;
	std::pair<bool,bool> test_results_;
	data_storage::SmartSegmentOP active_segment_;
	data_storage::SmartSegmentOP reference_segment_;
	core::Size active_resnum_;
	core::Size reference_resnum_;
	core::Size partner_resnum_;
	//core::Size active_atom_number_;
	//core::Size reference_atom_number_;
	data_storage::SmartSewingResidueOP active_residue_;
	data_storage::SmartSewingResidueOP reference_residue_;
	core::Size current_clashes_;
};


} //protocols
} //sewing
} //requirements



#endif //INCLUDED_protocols_sewing_requirements_ClashRequirement_hh



