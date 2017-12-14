// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/task/operation/EnableMultiCoolAnnealer.cc
/// @brief Task Operation for turning on the multi-cool annealer
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#include <core/pack/task/operation/EnableMultiCoolAnnealer.hh>
#include <core/pack/task/operation/EnableMultiCoolAnnealerCreator.hh>

#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/task_op_schemas.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

static basic::Tracer TR( "core.pack.task.operation.EnableMultiCoolAnnealer" );

namespace core {
namespace pack {
namespace task {
namespace operation {

EnableMultiCoolAnnealer::EnableMultiCoolAnnealer():
	TaskOperation(),
	increase_history_size_( false ),
	history_size_( 0 )
{}

EnableMultiCoolAnnealer::~EnableMultiCoolAnnealer() = default;

TaskOperationOP
EnableMultiCoolAnnealer::clone() const {
	return TaskOperationOP( new EnableMultiCoolAnnealer( *this ) );
}

EnableMultiCoolAnnealer::EnableMultiCoolAnnealer( EnableMultiCoolAnnealer const & /*src*/ ) = default;

void
EnableMultiCoolAnnealer::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & )
{
	if ( tag->hasOption( "history_size" ) ) {
		increase_history_size_ = true;
		history_size_ = tag->getOption< core::Size >( "history_size" );
	}
}

void
EnableMultiCoolAnnealer::apply( core::pose::Pose const &, core::pack::task::PackerTask & task ) const
{
	task.or_multi_cool_annealer( true );
	if ( increase_history_size_ ) {
		task.increase_multi_cool_annealer_history_size( history_size_ );
	}
}

std::string
EnableMultiCoolAnnealer::keyname() {
	return "EnableMultiCoolAnnealer";
}

void
EnableMultiCoolAnnealer::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attrs;
	attrs + XMLSchemaAttribute( "history_size", xsct_non_negative_integer, "The number of low-energy states that the MultiCoolAnnealer will keep from its initial cooling trajectory that will be explored at very low temperature during the final optimization phase. Default of 10." );
	task_op_schema_w_attributes( xsd, keyname(), attrs, "This task operation will turn on the MultiCoolAnnealer when invoking pack_rotamers; the MultiCoolAnnealer produces much lower energies than the regular annealer by spending more time at lower temperatures, and by quenching from many different low-energy rotamer assignments that it encounters during an initial cooling trajectory. See Leaver-Fay, Jacak, Stranges and Kuhlman, PLoS One 2011" );
}

core::pack::task::operation::TaskOperationOP
EnableMultiCoolAnnealerCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new EnableMultiCoolAnnealer );
}

std::string
EnableMultiCoolAnnealerCreator::keyname() const
{
	return EnableMultiCoolAnnealer::keyname();
}

void
EnableMultiCoolAnnealerCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	EnableMultiCoolAnnealer::provide_xml_schema( xsd );
}

} //core
} //pack
} //task
} //operation
