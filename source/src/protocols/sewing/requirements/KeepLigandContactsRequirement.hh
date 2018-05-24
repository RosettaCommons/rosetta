// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/requirements/KeepLigandContactsRequirement.hh
/// @brief a Requirement that an Assembly lose no more than a certain number of ligand contacts during conformer switches
/// @author frankdt (frankdt@email.unc.edu)


#ifndef INCLUDED_protocols_sewing_requirements_KeepLigandContactsRequirement_hh
#define INCLUDED_protocols_sewing_requirements_KeepLigandContactsRequirement_hh

#include <protocols/sewing/requirements/KeepLigandContactsRequirement.fwd.hh>
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
class KeepLigandContactsRequirement : public AssemblyRequirement {

public:

	KeepLigandContactsRequirement();
	KeepLigandContactsRequirement(KeepLigandContactsRequirement const & src);

	virtual ~KeepLigandContactsRequirement() override;

	KeepLigandContactsRequirementOP
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


	core::Real
	get_contact_distance_cutoff() const;

	void
	set_contact_distance_cutoff( core::Real );

	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & );

	static utility::tag::AttributeList
	get_xml_attributes();

	static std::string
	type_name();

private:
	core::Real contact_distance_cutoff_=2.0;
	std::pair<bool,bool> test_results_;
};


} //protocols
} //sewing
} //requirements



#endif //INCLUDED_protocols_sewing_requirements_KeepLigandContactsRequirement_hh



