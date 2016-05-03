// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file /DisableCDRsOperation.hh
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_antibody_task_operations_DisableCDRsOperation_hh
#define INCLUDED_protocols_antibody_task_operations_DisableCDRsOperation_hh

#include <protocols/antibody/task_operations/DisableCDRsOperation.fwd.hh>
#include <protocols/antibody/AntibodyInfo.fwd.hh>
#include <protocols/antibody/AntibodyEnum.hh>

#include <core/pack/task/operation/TaskOperation.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace antibody {
namespace task_operations {


/// @brief TaskOperation to Disable Packing and/or design of a set of CDRs.
/// By default, disables both Packing and Design of all CDRs.
///
/// @details See options for setting specific CDRs and whether to only disable design.
///
class DisableCDRsOperation : public core::pack::task::operation::TaskOperation {
public:

	/// @brief Default Constructor.
	DisableCDRsOperation();

	/// @brief Regular Constructor
	DisableCDRsOperation(AntibodyInfoCOP ab_info);

	/// @brief Constructor Specifying set of cdrs to use
	DisableCDRsOperation(AntibodyInfoCOP ab_info, utility::vector1<bool> const & cdrs);

	/// @brief Constructor Specifying set of cdrs and what to disable
	DisableCDRsOperation(
		AntibodyInfoCOP ab_info,
		utility::vector1<bool> const & cdrs,
		bool disable_packing_and_design);

	virtual ~DisableCDRsOperation();

	DisableCDRsOperation(DisableCDRsOperation const & src);

	virtual core::pack::task::operation::TaskOperationOP
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

	/// @brief Set the CDRs we will be disabling - 6/8 long vector corresponding to Enum.
	void
	set_cdrs(utility::vector1<bool> const & cdrs);

	/// @brief Set only a single CDR to disable.
	void
	set_cdr_only(CDRNameEnum cdr);

	/// @brief Set to disable packing and design, or only just design.
	void
	set_disable_packing_and_design(bool disable_packing_and_design);

	void
	set_defaults();

	static std::string keyname();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	AntibodyInfoCOP ab_info_;
	utility::vector1<bool> cdrs_; //This is so easy, but debugging sucks.  Wouldn't recommend vectors as storage in future.

	bool disable_packing_and_design_;

	///Needed for default and RS constructor.
	AntibodyNumberingSchemeEnum numbering_scheme_;
	CDRDefinitionEnum cdr_definition_;

};


} //task_operations
} //antibody
} //protocols

#endif //INCLUDED_protocols_antibody_task_operations_DisableCDRsOperation_hh



