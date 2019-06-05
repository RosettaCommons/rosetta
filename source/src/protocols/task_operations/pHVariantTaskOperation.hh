// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/task_operations/pHVariantTaskOperation
/// @brief During repacking, allow interchangeability of protonated and deprotonated variants
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_task_operations_pHVariantTaskOperation_hh
#define INCLUDED_protocols_task_operations_pHVariantTaskOperation_hh

#include <protocols/task_operations/pHVariantTaskOperation.fwd.hh>

#include <core/pack/task/operation/TaskOperation.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Utility headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace task_operations {

///@brief During repacking, allow interchangeability of protonated and deprotonated variants
class pHVariantTaskOperation: public core::pack::task::operation::TaskOperation {

public:

	pHVariantTaskOperation();

	pHVariantTaskOperation(pHVariantTaskOperation const & src);

	~pHVariantTaskOperation() override;

	core::pack::task::operation::TaskOperationOP
	clone() const override;

	/// @brief Configure from a RosettaScripts XML tag.
	void
	parse_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & ) override;

	//////////////////////

	/// @brief A TaskOp that enables sampling of protonation variants during repacking. The ionizable sidechains considered are ASP, GLU, LYS, TYR, and HIS. All other posiitons are restricted to repacking. The -pH_mode commandline option is still currently required to initially read the variants into the database.
	void
	apply(core::pose::Pose const & pose, core::pack::task::PackerTask & task) const override;

	/// @brief Return the name used to construct this TaskOperation from an XML file
	static std::string keyname();

	/// @brief Describe the format of XML file used to initialize this TaskOperation
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:

	// Selectors for ionizable residues
	core::select::residue_selector::ResidueSelectorOP select_asp_;
	core::select::residue_selector::ResidueSelectorOP select_glu_;
	core::select::residue_selector::ResidueSelectorOP select_tyr_;
	core::select::residue_selector::ResidueSelectorOP select_lys_;
	core::select::residue_selector::ResidueSelectorOP select_his_;


};

} //protocols
} //task_operations

#endif //INCLUDED_protocols/task_operations_pHVariantTaskOperation_fwd_hh
