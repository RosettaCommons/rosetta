// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/task_operations/SequenceMotifTaskOperation.cc
/// @brief  Creator for A TaskOp that takes a regex-like pattern and turns it into a set of design residues.  The string should identify what to do for each position.  An X indicates any residue, and is the same as [A_Z].  Anything other than a charactor should be in [].  The ^ denotes a not.  An example Regex for glycosylation is NX[ST].  This would design an ASN into the first position, skip the second, and allow S and T mutations only at the third position.  N[^P][ST] would denote that at the second position, we do not allow proline.  Sets of charactors using _ can be denoted, even though this doesnt really help us in design.  In the future one can imagine having sets of polars, etc. etc.  This TaskOp is like a simple RESFILE in a string.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_toolbox_task_operations_SequenceMotifTaskOperationCreator_HH
#define INCLUDED_protocols_toolbox_task_operations_SequenceMotifTaskOperationCreator_HH

#include <core/pack/task/operation/TaskOperationCreator.hh>

#include <string>


namespace protocols {
namespace toolbox {
namespace task_operations {

class SequenceMotifTaskOperationCreator : public core::pack::task::operation::TaskOperationCreator {
public:
	core::pack::task::operation::TaskOperationOP create_task_operation() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};



} //protocols
} //toolbox
} //task_operations


#endif //INCLUDED_protocols/toolbox/task_operations_SequenceMotifTaskOperationCreator_HH
