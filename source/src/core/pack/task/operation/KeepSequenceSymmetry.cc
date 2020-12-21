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
/// @author Updated by Tim Neary, timdot10@gmail.com


#include <core/pack/task/operation/KeepSequenceSymmetry.hh>
#include <core/pack/task/operation/KeepSequenceSymmetryCreator.hh>

#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/task_op_schemas.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/select/residue_selector/util.hh>

#include <basic/Tracer.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/vector0.hh>
#include <utility/string_util.hh>

// Cpp Headers
#include <string>

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
	if ( tag->hasOption( "name" ) ) {
		set_prefix_name( tag->getOption< std::string >( "name" ) );
	} else {
		utility_exit_with_message( "KeepSequenceSymmetry taskoperation must have a name." );
	}

	if ( tag->hasOption( "setting" ) ) {
		setting_ = tag->getOption< bool >( "setting" );
	}
}

void
KeepSequenceSymmetry::apply( core::pose::Pose const &, core::pack::task::PackerTask & task ) const
{
	task.keep_sequence_symmetry( setting_ );
	task.sequence_symmetric_uid_prefix( setup_magic_name_prefix_ );
}

std::string
KeepSequenceSymmetry::keyname() {
	return "KeepSequenceSymmetry";
}

void
KeepSequenceSymmetry::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	AttributeList attrs;
	attrs + XMLSchemaAttribute::attribute_w_default( "setting", xsct_rosetta_bool,
		"If true, Rosetta will activate the SequenceSymmetricAnnealer."
		"Use this when you give Rosetta a multimer to design and you want the sequences of the chains"
		" to be the same but you don't need strict physical symmetry.",
		"true" );
	/*
	attrs + XMLSchemaAttribute::attribute_w_default( "rotamer_selection_logic", xs_string,
	"Used to define the rotamer selection logic during the annealing process."
	" There are three possible options: \"identity\" which will only allow identical rotamers for"
	" each linked position, \"try_identity\" (default) which will first attempt identical rotamers"
	" but will also allow any rotamers of the same residue type if the initial attempt fails,"
	" \"type_only\" will only attempt to match the residue types for each linked position.",
	"identity" );
	*/

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & complex_type_name_for_task_op )
		.element_name( keyname() )
		.description( "Task operation to enforce sequence symmetry on." )
		.add_attributes( attrs )
		.add_optional_name_attribute()
		.write_complex_type_to_schema( xsd );
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
