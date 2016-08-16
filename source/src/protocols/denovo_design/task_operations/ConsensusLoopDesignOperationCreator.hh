// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/task_operations/ConsensusLoopDesignOperationCreator.cc
/// @brief Design residues that make up the highest-energy regions in a protein
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_task_operations_ConsensusLoopDesignOperationCreator_hh
#define INCLUDED_protocols_denovo_design_task_operations_ConsensusLoopDesignOperationCreator_hh

// unit headers

// package headers

// project headers
#include <core/pack/task/operation/TaskOperationCreator.hh>

namespace protocols {
namespace denovo_design {
namespace task_operations {

class ConsensusLoopDesignOperationCreator : public core::pack::task::operation::TaskOperationCreator {
public:
	virtual core::pack::task::operation::TaskOperationOP create_task_operation() const;
	virtual std::string keyname() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

}
}
}

#endif
