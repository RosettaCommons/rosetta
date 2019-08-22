// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/task/operation/KeepSequenceSymmetry.cc
/// @brief Use this when you give Rosetta a multimer to design and you want the sequences of the chains to be the same but you don't need strict physical symmetry.
/// @author Jack Maguire, jackmaguire1444@gmail.com

#include <core/pack/task/operation/KeepSequenceSymmetry.hh>
#include <core/pack/task/operation/KeepSequenceSymmetryCreator.hh>

#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/task_op_schemas.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

static basic::Tracer TR( "core.pack.task.operation.KeepSequenceSymmetry" );

namespace core {
namespace pack {
namespace task {
namespace operation {

KeepSequenceSymmetry::KeepSequenceSymmetry():
	TaskOperation(),
	setting_( true )
{}

KeepSequenceSymmetry::~KeepSequenceSymmetry() = default;

TaskOperationOP
KeepSequenceSymmetry::clone() const {
	return utility::pointer::make_shared< KeepSequenceSymmetry >( *this );
}

KeepSequenceSymmetry::KeepSequenceSymmetry(
	KeepSequenceSymmetry const & /*src*/
) = default;

void
KeepSequenceSymmetry::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & )
{
	if ( tag->hasOption( "setting" ) ) {
		setting_ = tag->getOption< bool >( "setting" );
	}
}

void
KeepSequenceSymmetry::apply( core::pose::Pose const &, core::pack::task::PackerTask & task ) const
{
	task.keep_sequence_symmetry( setting_ );
}

std::string
KeepSequenceSymmetry::keyname() {
	return "KeepSequenceSymmetry";
}

void
KeepSequenceSymmetry::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	std::string const description = "If true, Rosetta will activate the SequenceSymmetricAnnealer. Use this when you give Rosetta a multimer to design and you want the sequences of the chains to be the same but you don't need strict physical symmetry.";
	AttributeList attrs;
	attrs + XMLSchemaAttribute( "setting", xsct_rosetta_bool, description );
	task_op_schema_w_attributes( xsd, keyname(), attrs, description );
}

core::pack::task::operation::TaskOperationOP
KeepSequenceSymmetryCreator::create_task_operation() const
{
	return utility::pointer::make_shared< KeepSequenceSymmetry >();
}

std::string
KeepSequenceSymmetryCreator::keyname() const
{
	return KeepSequenceSymmetry::keyname();
}

void
KeepSequenceSymmetryCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	KeepSequenceSymmetry::provide_xml_schema( xsd );
}

} //core
} //pack
} //task
} //operation
