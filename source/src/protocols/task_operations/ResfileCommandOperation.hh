// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/task_operations/ResfileCommandOperation
/// @brief Applies the equivalent of a resfile line (without the resnums) to residues specified in a residue selector.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_task_operations_ResfileCommandOperation_hh
#define INCLUDED_protocols_task_operations_ResfileCommandOperation_hh

#include <protocols/task_operations/ResfileCommandOperation.fwd.hh>

#include <core/pack/task/operation/TaskOperation.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/pack/task/ResfileReader.hh>

// Utility headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace task_operations {

///@brief Applies the equivalent of a resfile line (without the resnums) to residues specified in a residue selector.
class ResfileCommandOperation: public core::pack::task::operation::TaskOperation {
public:

	ResfileCommandOperation();

	ResfileCommandOperation(
		core::select::residue_selector::ResidueSelectorCOP selector,
		std::string const & command
	);

	ResfileCommandOperation(ResfileCommandOperation const & src);

	~ResfileCommandOperation() override;

	core::pack::task::operation::TaskOperationOP
	clone() const override;

	/// @brief Configure from a RosettaScripts XML tag.
	void
	parse_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & ) override;

public:
	//////////////////////

	void
	apply(core::pose::Pose const & pose, core::pack::task::PackerTask & task) const override;

public:

	///@brief A resfile command string without any numbers in the front.
	///EX:
	/// POLAR
	/// EMPTY NC R2 NC T6 NC OP5
	///
	void
	set_command( std::string const & command);

	void
	set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector );

public:

	/// @brief Return the name used to construct this TaskOperation from an XML file
	static std::string keyname();

	/// @brief Describe the format of XML file used to initialize this TaskOperation
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:

	std::map< std::string, core::pack::task::ResfileCommandOP > command_map_;

	core::select::residue_selector::ResidueSelectorCOP selector_ = nullptr;
	std::string command_ = "";

	utility::vector1< core::pack::task::ResfileCommandOP > commands_;


};

} //protocols
} //task_operations

#endif //INCLUDED_protocols/task_operations_ResfileCommandOperation_fwd_hh
