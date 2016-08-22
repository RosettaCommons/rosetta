// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <protocols/denovo_design/components/StructureDataFactory.hh>

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
	return FileConstraintGeneratorCreator::constraint_generator_name();
}

protocols::constraint_generator::ConstraintGeneratorOP
FileConstraintGeneratorCreator::create_constraint_generator() const
{
	return protocols::constraint_generator::ConstraintGeneratorOP( new FileConstraintGenerator() );
}

std::string
FileConstraintGeneratorCreator::constraint_generator_name()
{
	return "FileConstraintGenerator";
}

/// @brief
FileConstraintGenerator::FileConstraintGenerator():
	protocols::constraint_generator::ConstraintGenerator( FileConstraintGeneratorCreator::constraint_generator_name() ),
	filename_( "" )
{
}

/// @brief
FileConstraintGenerator::FileConstraintGenerator( std::string const & filename ):
	protocols::constraint_generator::ConstraintGenerator( FileConstraintGeneratorCreator::constraint_generator_name() ),
	filename_( filename )
{
}

/// @brief
FileConstraintGenerator::~FileConstraintGenerator() {}

protocols::constraint_generator::ConstraintGeneratorOP
FileConstraintGenerator::clone() const
{
	return protocols::constraint_generator::ConstraintGeneratorOP( new FileConstraintGenerator( *this ) );
}

void
FileConstraintGenerator::parse_tag(
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap & )
{
	set_cstfile( tag->getOption< std::string >( "filename", filename_ ) );
	if ( filename_.empty() ) {
		throw utility::excn::EXCN_RosettaScriptsOption( "FileConstraintGenerator requires the 'filename' option to be set." );
	}
}

/// @brief
core::scoring::constraints::ConstraintCOPs
FileConstraintGenerator::apply( core::pose::Pose const & pose ) const
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

	// Get stream of cst definitions
	std::stringstream cst_stream;
	components::StructureDataFactory const & factory = *components::StructureDataFactory::get_instance();
	if ( factory.has_cached_data( pose ) ) {
		components::StructureData const & sd = factory.get_from_const_pose( pose );
		cst_stream << clean_constraint_string( sd.substitute_variables( infile ) );
	} else {
		cst_stream << infile.rdbuf();
	}
	infile.close();
	TR.Debug << "Parsed " << filename_ << " : " << cst_stream.str() << std::endl;

	core::scoring::constraints::ConstraintSetOP constraints( new core::scoring::constraints::ConstraintSet );
	constraints = core::scoring::constraints::ConstraintIO::get_instance()->read_constraints(
		cst_stream,
		constraints,
		pose );
	TR.Debug << "generated " << constraints->get_all_constraints().size() << " constraints" << std::endl;

	return constraints->get_all_constraints();
} // apply

/// @brief
void
FileConstraintGenerator::set_cstfile( std::string const & filename )
{
	filename_ = filename;
}

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
