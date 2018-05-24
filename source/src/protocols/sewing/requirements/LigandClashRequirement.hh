// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/requirements/LigandClashRequirement.hh
/// @brief Assembly requirement that clash checks the assembly backbone with its bound ligands
/// @author Minnie Langlois (minnie@email.unc.edu)


#ifndef INCLUDED_protocols_sewing_requirements_LigandClashRequirement_hh
#define INCLUDED_protocols_sewing_requirements_LigandClashRequirement_hh

//Project headers
//#include <protocols/sewing/requirements/ClashRequirement.hh>
#include <protocols/sewing/requirements/LigandAssemblyRequirement.hh>
#include <protocols/sewing/requirements/LigandClashRequirement.fwd.hh>
#include <protocols/sewing/data_storage/SmartAssembly.fwd.hh>
#include <protocols/sewing/data_storage/SmartSegment.fwd.hh>
#include <protocols/sewing/data_storage/SmartSewingResidue.fwd.hh>
#include <protocols/sewing/data_storage/LigandResidue.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/types.hh>
// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace sewing {
namespace requirements {

///@brief Assembly requirement that clash checks the assembly backbone with its bound ligands
class LigandClashRequirement : public LigandAssemblyRequirement {

public:

	LigandClashRequirement();
	LigandClashRequirement(LigandClashRequirement const & src);

	virtual ~LigandClashRequirement();

	LigandClashRequirementOP
	clone() const;

	std::pair<bool,bool>
	test(data_storage::SmartAssemblyOP assembly) override;

	std::string
	get_name() override{
		return type_name();
	}

	//Getters and Setters
	core::Size
	get_maximum_clashes_allowed() const;

	core::Real
	get_clash_radius() const;

	void
	set_maximum_clashes_allowed( core::Size );

	void
	set_clash_radius( core::Real );


	void
	set_options_from_tag(
		utility::tag::TagCOP requirement_tag,
		basic::datacache::DataMap& datamap,
		protocols::filters::Filters_map const & filtermap,
		protocols::moves::Movers_map const & movermap,
		core::pose::Pose const & pose) override;


	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & );

	static std::string
	type_name();

private:
	core::Size maximum_clashes_allowed_;
	core::Real clash_radius_;
	data_storage::LigandResidueOP active_ligand_residue_;
	data_storage::LigandContactOP active_ligand_contact_;
	core::chemical::AtomTypeSetCAP atom_types_;
	core::Size current_clashes_;
	data_storage::SmartSegmentOP active_segment_;
	core::Size active_resnum_;
	data_storage::SmartSewingResidueOP active_residue_;
	std::pair< bool,bool > test_results_;

};


} //protocols
} //sewing
} //requirements



#endif //INCLUDED_protocols_sewing_requirements_LigandClashRequirement_hh





