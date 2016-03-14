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
#include <protocols/denovo_design/components/StructureData.hh>

// Core headers
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/id/SequenceMapping.hh>

// Utility headers
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.constraints.FileConstraintGenerator" );

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
FileConstraintGenerator::FileConstraintGenerator( std::string const & filename ):
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
FileConstraintGenerator::set_cstfile( std::string const & filename )
{
	filename_ = filename;
}

/// @brief
core::scoring::constraints::ConstraintCOPs
FileConstraintGenerator::generate_constraints( Pose const & pose )
{
	if ( filename_ == "" ) {
		utility_exit_with_message( "FileConstraintGenerator requires that a constraint filename be specified." );
	}

	// Read constraint file and substitute any StructureData variables
	// File format: %%SEGMENT_NAME#resid%% will be replaced by an integer resid
	TR.Debug << "Opening " << filename_ << std::endl;
	utility::io::izstream infile( filename_ );
	if ( !infile.good() ) {
		throw utility::excn::EXCN_BadInput( "Could not open " + filename_ + " for reading." );
	}
	components::StructureDataOP sd = components::StructureData::create_from_pose( pose, "FileConstraintGenerator" );
	debug_assert( sd );
	std::stringstream cst_stream( clean_constraint_string( sd->substitute_variables( infile ) ) );
	infile.close();
	TR.Debug << "Parsed " << filename_ << " : " << cst_stream.str() << std::endl;

	core::scoring::constraints::ConstraintSetOP constraints( new core::scoring::constraints::ConstraintSet );
	constraints = core::scoring::constraints::ConstraintIO::get_instance()->read_constraints(
		cst_stream,
		constraints,
		pose );
	TR.Debug << "generated " << constraints->get_all_constraints().size() << " constraints" << std::endl;

	return constraints->get_all_constraints();
} //generate constraints

std::string
FileConstraintGenerator::clean_constraint_string( std::string const & cst_str ) const
{
	utility::vector1< std::string > const lines = utility::string_split( cst_str, '\n' );
	std::string newstr = "";
	for ( utility::vector1< std::string >::const_iterator l=lines.begin(); l!=lines.end(); ++l ) {
		if ( l->find_first_not_of( "\t\n\v\f\r " ) != std::string::npos ) newstr += *l + "\n";
	}
	return newstr;
}

} //namespace constraints
} //namespace denovo_design
} //namespace protocols
