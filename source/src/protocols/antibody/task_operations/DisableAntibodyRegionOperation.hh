// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/task_operations/DisableAntibodyRegionOperation.hh
/// @brief Task Operation to Disable a Region of Antibody
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_antibody_task_operations_DisableAntibodyRegionOperation_hh
#define INCLUDED_protocols_antibody_task_operations_DisableAntibodyRegionOperation_hh

#include <protocols/antibody/task_operations/DisableAntibodyRegionOperation.fwd.hh>
#include <protocols/antibody/AntibodyInfo.fwd.hh>
#include <protocols/antibody/AntibodyEnum.hh>

#include <core/pack/task/operation/TaskOperation.hh>


// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace antibody {
namespace task_operations {


/// @brief A TaskOperation that Disables packing +/or design of a particular antibody region.
/// By Default, disables packing and design of the cdr_region.  Make sure to set the region you want disabled.
///
class DisableAntibodyRegionOperation : public core::pack::task::operation::TaskOperation {
public:

	/// @brief Default constructor.  Do not use this.
	DisableAntibodyRegionOperation();

	/// @brief Constructor setting only AntibodyInfo.
	DisableAntibodyRegionOperation(AntibodyInfoCOP ab_info);

	/// @brief Constructor setting AntibodyInfo and the region.
	DisableAntibodyRegionOperation(AntibodyInfoCOP ab_info, AntibodyRegionEnum region);

	/// @brief Constructor setting the region and to fully disable the region.
	/// If disable_packing_and_design is set to False, will only disable design.
	DisableAntibodyRegionOperation(AntibodyInfoCOP ab_info, AntibodyRegionEnum region, bool disable_packing_and_design);

	DisableAntibodyRegionOperation(DisableAntibodyRegionOperation const & src);

	virtual ~DisableAntibodyRegionOperation();

	core::pack::task::operation::TaskOperationOP
	clone() const;

	/// @brief Configure from a RosettaScripts XML tag.
	virtual void
	parse_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & );

	//////////////////////

	virtual
	void
	apply(core::pose::Pose const & pose, core::pack::task::PackerTask & task) const;

	void
	set_region(AntibodyRegionEnum region);

	/// @brief If disable_packing_and_design is set to False, will only disable design.
	void
	set_disable_packing_and_design(bool disable_packing_and_design);


	void
	set_defaults();

	static std::string keyname();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:

	AntibodyInfoCOP ab_info_;
	AntibodyRegionEnum region_;

	bool disable_packing_and_design_;

	///Needed for default and RS constructor.
	AntibodyNumberingSchemeEnum numbering_scheme_;
	CDRDefinitionEnum cdr_definition_;

};

} //task_operations
} //antibody
} //protocols

#endif //INCLUDED_protocols_antibody_task_operations_DisableAntibodyRegionOperation_hh



