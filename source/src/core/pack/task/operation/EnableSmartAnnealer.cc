// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/task/operation/EnableSmartAnnealer.cc
/// @author Jack Maguire, jackmaguire1444@gmail.com

#include <core/pack/task/operation/EnableSmartAnnealer.hh>
#include <core/pack/task/operation/EnableSmartAnnealerCreator.hh>

#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/task_op_schemas.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

static basic::Tracer TR( "core.pack.task.operation.EnableSmartAnnealer" );

namespace core {
namespace pack {
namespace task {
namespace operation {

EnableSmartAnnealer::EnableSmartAnnealer():
	TaskOperation()
{}

EnableSmartAnnealer::~EnableSmartAnnealer() = default;

TaskOperationOP
EnableSmartAnnealer::clone() const {
	return utility::pointer::make_shared< EnableSmartAnnealer >( *this );
}

EnableSmartAnnealer::EnableSmartAnnealer( EnableSmartAnnealer const & /*src*/ ) = default;

void
EnableSmartAnnealer::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & )
{
	if ( tag->hasOption( "model" ) ) {
		smart_annealer_model_ = tag->getOption< std::string >( "model" );
	}
	if ( tag->hasOption( "cutoff" ) ) {
		smart_annealer_cutoff_ = tag->getOption< core::Real >( "cutoff" );
	}
	if ( tag->hasOption( "pick_again" ) ) {
		smart_annealer_pick_again_ = tag->getOption< bool >( "pick_again" );
	}
	if ( tag->hasOption( "disable_during_quench" ) ) {
		smart_annealer_disable_during_quench_ = tag->getOption< bool >( "disable_during_quench" );
	}
}

void
EnableSmartAnnealer::apply( core::pose::Pose const &, core::pack::task::PackerTask & task ) const
{
	task.set_smart_annealer( true );
	task.set_smart_annealer_model( smart_annealer_model_ );
	task.set_smart_annealer_cutoff( smart_annealer_cutoff_ );
	task.set_smart_annealer_pick_again( smart_annealer_pick_again_ );
	task.set_smart_annealer_disable_during_quench( smart_annealer_disable_during_quench_ );
}

std::string
EnableSmartAnnealer::keyname() {
	return "EnableSmartAnnealer";
}

void
EnableSmartAnnealer::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attrs;
	attrs + XMLSchemaAttribute::attribute_w_default( "model", xs_string, "Choose which neural network to use for the smart annealer. Look at database/protocol_data/tensorflow_graphs/smart_annealer/ to see the options.", "generation2" );
	attrs + XMLSchemaAttribute::attribute_w_default( "cutoff", xsct_real, "Choose a number from 0 to 1 to tune how aggressive the smart annealer is. Higher numbers are more agressive (risky) but have a potentially greater speedup (speedup requires pick_again=false)", "0.25" );
	attrs + XMLSchemaAttribute::attribute_w_default( "pick_again", xsct_rosetta_bool, "f disabled, the smart annealer just skips unfruitful amino acids. Enabling this option tells the annealer to pick a fruitful rotamer to sample this round instead of skipping the round. Will not give you a speedup but may give you a better final outcome.", "true" );
	attrs + XMLSchemaAttribute::attribute_w_default( "disable_during_quench", xsct_rosetta_bool, "Run the final quenching stage as normal, regardless of how bad an amino acid may be.", "true" );
	task_op_schema_w_attributes( xsd, keyname(), attrs, "The smart annealer uses tensorflow to decrease the sample space of a typical packing run" );
}

core::pack::task::operation::TaskOperationOP
EnableSmartAnnealerCreator::create_task_operation() const
{
	return utility::pointer::make_shared< EnableSmartAnnealer >();
}

std::string
EnableSmartAnnealerCreator::keyname() const
{
	return EnableSmartAnnealer::keyname();
}

void
EnableSmartAnnealerCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	EnableSmartAnnealer::provide_xml_schema( xsd );
}

} //operation
} //task
} //pack
} //core
