// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/toolbox/task_operations/SequenceMotifTaskOperation
/// @brief A TaskOp that takes a regex-like pattern and turns it into a set of design residues.  The string should identify what to do for each position.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_toolbox_task_operations_SequenceMotifTaskOperation_hh
#define INCLUDED_protocols_toolbox_task_operations_SequenceMotifTaskOperation_hh

#include <protocols/toolbox/task_operations/SequenceMotifTaskOperation.fwd.hh>

#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/ResfileReader.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>



// Utility headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace toolbox {
namespace task_operations {

///@brief A TaskOp that takes a regex-like pattern and turns it into a set of design residues.  The string should identify what to do for each position.
///
///@details
///
///  See the MOTIF option for a description of the syntax used to define a motif.
class SequenceMotifTaskOperation: public core::pack::task::operation::TaskOperation {
public:

	SequenceMotifTaskOperation();

	SequenceMotifTaskOperation(
		core::select::residue_selector::ResidueSelectorCOP selector,
		std::string const & motif );

	SequenceMotifTaskOperation(SequenceMotifTaskOperation const & src);

	~SequenceMotifTaskOperation() override;

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

	///@brief Set a string that tells this operation how to design.
	///@details
	///  This is slightly similar to a regex, but not quite. We are not matching a sequence, we are designing in a motif regardless of the current sequence, anywhere in a protein.
	///
	///   - Each letter corresponds to a position. Using [ ] indicates a more complicated expression for that position.
	///   - An X indicates it can be anything, and that we are designing here.
	///   - An AA Letter, like V, indicates that that position will be designed to a V.
	///   - A - charactor indicates that that position stays with whatever it is currently.  We essentially skip this position.
	///   - An expression like: [^PAV] indicates that we will design anything except Proline, Alanine, and Valine
	///   - An expression like: [NTS] indicates that that position can be Asparigine, Threonine, or Serine and only of these will be enabled during the design.
	///   - RESFILE commands are accepted as well. Resfile commands need a % in front of the whole set.
	///      FOr example [%POLAR] is totally cool.  This is a whole resfile line, so NCs also work:
	///       [%EMPTY NC R2 NC T6 NC OP5]
	///
	/// EXAMPLE:
	///  Glycosylation N-Linked motif design: N[^P][ST]
	///
	void
	set_motif( std::string const & motif );

	void
	set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector );

public:

	/// @brief Return the name used to construct this TaskOperation from an XML file
	static std::string keyname();

	/// @brief Describe the format of XML file used to initialize this TaskOperation
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:

	///@brief Setup our private list of commands.  THis saves a lot of time instead of needing to do this at every apply.
	void
	setup_commands();

private:
	core::select::residue_selector::ResidueSelectorCOP selector_ = nullptr;
	std::string motif_ = "";
	utility::vector1< utility::vector1< core::pack::task::ResfileCommandOP >> commands_;


};

} //protocols
} //toolbox
} //task_operations

#endif //INCLUDED_protocols/toolbox/task_operations_SequenceMotifTaskOperation_fwd_hh
