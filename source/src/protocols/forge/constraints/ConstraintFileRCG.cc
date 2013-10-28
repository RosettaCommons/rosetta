// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/forge/constraints/ConstraintFileRCG.cc
///
/// @brief
/// @author Nobuyasu Koga( nobuyasu@uw.edu ) , October 2009

// Unit header
#include <protocols/forge/constraints/ConstraintFileRCG.hh>
#include <protocols/forge/constraints/ConstraintFileCstGeneratorCreator.hh>

// Package headers
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/id/SequenceMapping.hh>

// Project headers
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>


static basic::Tracer tr( "protocols.forge.constraints.ConstraintFileRCG" );

namespace protocols{
namespace forge{
namespace constraints{

std::string
ConstraintFileCstGeneratorCreator::keyname() const
{
	return ConstraintFileCstGeneratorCreator::mover_name();
}

protocols::moves::MoverOP
ConstraintFileCstGeneratorCreator::create_mover() const
{
	return new ConstraintFileRCG();
}

std::string
ConstraintFileCstGeneratorCreator::mover_name()
{
	return "ConstraintFileCstGenerator";
}

/// @brief
ConstraintFileRCG::ConstraintFileRCG():
	RemodelConstraintGenerator(),
	filename_( "" )
{
	/*using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if( option[ constraints::cst_file ].user() ){
		filename_ = option[ constraints::cst_file ].value().at( 1 );
		}*/
}

/// @brief
ConstraintFileRCG::ConstraintFileRCG( String const & filename ):
	RemodelConstraintGenerator(),
	filename_( filename )
{}

/// @brief
ConstraintFileRCG::~ConstraintFileRCG() {}

void
ConstraintFileRCG::parse_my_tag( TagCOP const tag,
												basic::datacache::DataMap & data,
												protocols::filters::Filters_map const & filters,
												protocols::moves::Movers_map const & movers,
												core::pose::Pose const & pose )
{
	RemodelConstraintGenerator::parse_my_tag( tag, data, filters, movers, pose );
	set_cstfile( tag->getOption< std::string >( "filename", filename_ ) );
}

std::string
ConstraintFileRCG::get_name() const
{
	return ConstraintFileCstGeneratorCreator::mover_name();
}

protocols::moves::MoverOP
ConstraintFileRCG::fresh_instance() const
{
	return new ConstraintFileRCG();
}

/// @brief
void
ConstraintFileRCG::set_cstfile( String const & filename )
{
	filename_ = filename;
}

/// @brief
void
ConstraintFileRCG::generate_remodel_constraints( Pose const & pose )
{
	using namespace core::scoring::constraints;
	if ( filename_ == "" ) {
		utility_exit_with_message( "ConstraintFileCstGenerator requires that a constraint filename be specified." );
	}
	ConstraintSetOP constraints = ConstraintIO::get_instance()->read_constraints( filename_, new ConstraintSet, pose );
	this->add_constraints( constraints->get_all_constraints() );
} //generate constraints


} //namespace constraints
} //namespace forge
} //namespace protocols
