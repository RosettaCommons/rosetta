// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/task/operation/FillAUTOTaskOperation
/// @brief fills the AUTO behavior for all residues in Task. Useful if a protocol expects AUTO-style resfile, but no resfile present.
/// @author Steven Lewis (smlewi@gmail.com)


#ifndef INCLUDED_core_pack_task_operation_FillAUTOTaskOperation_hh
#define INCLUDED_core_pack_task_operation_FillAUTOTaskOperation_hh

#include <core/pack/task/operation/FillAUTOTaskOperation.fwd.hh>

#include <core/pack/task/operation/TaskOperation.hh>

// Utility headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace core {
namespace pack {
namespace task {
namespace operation {

///@brief fills the AUTO behavior for all residues in Task. Useful if a protocol expects AUTO-style resfile, but no resfile present.
class FillAUTOTaskOperation: public core::pack::task::operation::TaskOperation {
public:

	FillAUTOTaskOperation();

	FillAUTOTaskOperation(FillAUTOTaskOperation const & src);

	~FillAUTOTaskOperation() override;

	core::pack::task::operation::TaskOperationOP
	clone() const override;

	/// @brief Configure from a RosettaScripts XML tag.
	void
	parse_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & ) override;

	//////////////////////

	void
	apply(core::pose::Pose const & pose, core::pack::task::PackerTask & task) const override;

	/// @brief Return the name used to construct this TaskOperation from an XML file
	static std::string keyname();

	/// @brief Describe the format of XML file used to initialize this TaskOperation
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:


};

} //core
} //pack
} //task
} //operation

#endif //INCLUDED_core/pack/task/operation_FillAUTOTaskOperation_fwd_hh
