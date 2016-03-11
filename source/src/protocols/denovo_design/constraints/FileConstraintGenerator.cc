// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/denovo_design/constraints/FileConstraintGenerator.cc
///
/// @brief
/// @author Nobuyasu Koga( nobuyasu@uw.edu ) , October 2009
/// @author Tom Linsky ( tlinsky at uw dot edu ), Mar 2016

// Unit header
#include <protocols/denovo_design/constraints/FileConstraintGenerator.hh>
#include <protocols/denovo_design/constraints/FileConstraintGeneratorCreator.hh>

// Package headers
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/id/SequenceMapping.hh>

// Project headers
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>


static THREAD_LOCAL basic::Tracer tr( "protocols.denovo_design.constraints.FileConstraintGenerator" );

namespace protocols {
namespace denovo_design {
namespace constraints {

std::string
FileConstraintGeneratorCreator::keyname() const
{
	return FileConstraintGeneratorCreator::mover_name();
}

protocols::moves::MoverOP
FileConstraintGeneratorCreator::create_mover() const
{
	return protocols::moves::MoverOP( new FileConstraintGenerator() );
}

std::string
FileConstraintGeneratorCreator::mover_name()
{
	return "FileConstraintGenerator";
}

/// @brief
FileConstraintGenerator::FileConstraintGenerator():
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
FileConstraintGenerator::FileConstraintGenerator( String const & filename ):
	RemodelConstraintGenerator(),
	filename_( filename )
{}

/// @brief
FileConstraintGenerator::~FileConstraintGenerator() {}

void
FileConstraintGenerator::parse_my_tag( TagCOP const tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose )
{
	RemodelConstraintGenerator::parse_my_tag( tag, data, filters, movers, pose );
	set_cstfile( tag->getOption< std::string >( "filename", filename_ ) );
}

std::string
FileConstraintGenerator::get_name() const
{
	return FileConstraintGeneratorCreator::mover_name();
}

protocols::moves::MoverOP
FileConstraintGenerator::fresh_instance() const
{
	return protocols::moves::MoverOP( new FileConstraintGenerator() );
}

protocols::moves::MoverOP
FileConstraintGenerator::clone() const
{
	return protocols::moves::MoverOP( new FileConstraintGenerator( *this ) );
}

/// @brief
void
FileConstraintGenerator::set_cstfile( String const & filename )
{
	filename_ = filename;
}

/// @brief
core::scoring::constraints::ConstraintCOPs
FileConstraintGenerator::generate_constraints( Pose const & pose )
{
	using namespace core::scoring::constraints;
	if ( filename_ == "" ) {
		utility_exit_with_message( "FileConstraintGenerator requires that a constraint filename be specified." );
	}
	ConstraintSetOP constraints = ConstraintIO::get_instance()->read_constraints( filename_, ConstraintSetOP( new ConstraintSet ), pose );
	return constraints->get_all_constraints();
} //generate constraints


} //namespace constraints
} //namespace denovo_design
} //namespace protocols
