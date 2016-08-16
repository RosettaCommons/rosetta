// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/operation/OptH.cc
/// @brief  run optH
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <core/pack/task/operation/OptH.hh>
#include <core/pack/task/operation/OptHCreator.hh>

// project headers
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

// utility headers
#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace task {
namespace operation {


/// @brief default constructor
OptH::OptH() :
	Super(),
	init_from_command_line_( true ),
	include_current_( true ),
	flip_HNQ_( false ),
	use_multicool_annealer_( false )
{}


/// @brief copy constructor
OptH::OptH( OptH const & rval ) :
	Super( rval ),
	disallowed_resids_( rval.disallowed_resids_ ),
	init_from_command_line_( rval.init_from_command_line_ ),
	include_current_( rval.include_current_ ),
	flip_HNQ_( rval.flip_HNQ_ ),
	use_multicool_annealer_( rval.use_multicool_annealer_ )
{}


/// @brief default destructor
OptH::~OptH() {}

/// @brief clone this object
OptH::TaskOperationOP OptH::clone() const {
	return OptH::TaskOperationOP( new OptH( *this ) );
}


/// @brief apply operations to PackerTask
void OptH::apply( Pose const & , PackerTask & task ) const {
	// optH parameters
	if ( init_from_command_line_ ) {
		task.initialize_from_command_line();
	}
	task.or_optimize_h_mode( true );
	task.or_include_current( include_current_ );
	task.or_flip_HNQ( flip_HNQ_ );
	task.or_multi_cool_annealer( use_multicool_annealer_ );

	// lock down disallowed residues
	for ( utility::vector1< Size >::const_iterator i = disallowed_resids_.begin(), ie = disallowed_resids_.end(); i != ie; ++i ) {
		task.nonconst_residue_task( *i ).prevent_repacking();
	}
}


/// @brief prevent a position from being optimized
void OptH::disallow_resid( Size const resid ) {
	disallowed_resids_.push_back( resid );
}


/// @brief init flags from command line? (default true)
void OptH::init_from_comand_line( bool const flag ) {
	init_from_command_line_ = flag;
}


/// @brief include current sidechain rotamers? (default true)
void OptH::include_current( bool const flag ) {
	include_current_ = flag;
}


/// @brief allow sidechain flips of HNQ? (default false)
void OptH::flip_HNQ( bool const flag ) {
	flip_HNQ_ = flag;
}


/// @brief use multicool annealer? (default false)
void OptH::use_multicool_annealer( bool const flag ) {
	use_multicool_annealer_ = flag;
}

std::string OptH::keyname() { return "OptH"; }

void OptH::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	task_op_schema_empty( xsd, keyname() );
}

TaskOperationOP OptHCreator::create_task_operation() const
{
	return TaskOperationOP( new OptH );
}

std::string OptHCreator::keyname() const { return OptH::keyname(); }

void OptHCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	OptH::provide_xml_schema( xsd );
}


} // namespace operation
} // namespace task
} // namespace pack
} // namespace core
