// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/quantum_annealing/ExternalPackerResultLoader.cc
/// @brief A mover that reads the result of an externally-called packer, and applies it to a pose to generate a structure.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Unit headers
#include <protocols/quantum_annealing/ExternalPackerResultLoader.hh>
#include <protocols/quantum_annealing/ExternalPackerResultLoaderCreator.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pack/interaction_graph/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/pointer/memory.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.quantum_annealing.ExternalPackerResultLoader" );

namespace protocols {
namespace quantum_annealing {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
ExternalPackerResultLoader::ExternalPackerResultLoader():
	protocols::moves::Mover( ExternalPackerResultLoader::mover_name() ),
	rotamer_selection_file_data_(""),
	packer_problem_file_data_(""),
	problem_was_loaded_(false),
	solution_was_loaded_(false)
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
ExternalPackerResultLoader::~ExternalPackerResultLoader(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
ExternalPackerResultLoader::apply( core::pose::Pose & pose ){
	std::string const errmsg("Error in protocols::quantum_annealing::ExternalPackerResultLoader::apply(): ");
	runtime_assert_string_msg( problem_was_loaded_ && solution_was_loaded_, errmsg + "A packer problem definition file and a packer solution file must be loaded before the apply() function can be called." );
	core::pack::interaction_graph::set_externally_generated_packer_solution( pose, packer_problem_file_data_, rotamer_selection_file_data_ );
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
ExternalPackerResultLoader::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
/// @details Warning!  This mover reads from disk at parse time!
void
ExternalPackerResultLoader::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	TR << "Parsing options for ExternalPackerResultLoader." << std::endl;
	TR.Warning << TR.Red << "Warning!  The ExternalPackerResultLoader discards the input pose, and builds a new pose based on the description in the packer problem definition file." << TR.Reset << std::endl;
	std::string const problem_file( tag->getOption<std::string>("packer_problem_definition_file") );
	std::string const solution_file( tag->getOption<std::string>("solution_definition_file") );
	read_packer_problem_file(problem_file);
	read_rotamer_selection_file(solution_file);
}
void ExternalPackerResultLoader::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;
	// TO DO: perhaps this is not the right function to call? -- also, delete this comment
	attlist + XMLSchemaAttribute( "packer_problem_definition_file", xs_string, "The file from which a description of the packer problem will be loaded.")
		+ XMLSchemaAttribute( "solution_definition_file", xs_string, "The file from which the solution to the packer problem will be loaded.");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "This mover allows a solution to a packer problem to be imported from an external file.  This allows Rosetta to communicate with external optimization packages, and provides a means of visualizing the result easily.", attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
ExternalPackerResultLoader::fresh_instance() const
{
	return protocols::moves::MoverOP( utility::pointer::make_shared< ExternalPackerResultLoader >() );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
ExternalPackerResultLoader::clone() const
{
	return protocols::moves::MoverOP( utility::pointer::make_shared< ExternalPackerResultLoader >( *this ) );
}

std::string ExternalPackerResultLoader::get_name() const {
	return mover_name();
}

std::string ExternalPackerResultLoader::mover_name() {
	return "ExternalPackerResultLoader";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
ExternalPackerResultLoaderCreator::create_mover() const
{
	return protocols::moves::MoverOP( new ExternalPackerResultLoader );
}

std::string
ExternalPackerResultLoaderCreator::keyname() const
{
	return ExternalPackerResultLoader::mover_name();
}

void ExternalPackerResultLoaderCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ExternalPackerResultLoader::provide_xml_schema( xsd );
}

///////////////////////
/// public methods  ///
///////////////////////

/// @brief Reads a description of the packer problem from disk and caches it.
/// @details WARNING!  READS FROM DISK!
void
ExternalPackerResultLoader::read_packer_problem_file(
	std::string const &filename
) {
	TR << "Reading packer problem from \"" << filename << "\"." << std::endl;
	set_packer_problem_string( utility::file_contents(filename) );
}

/// @brief Reads a description of the packer solution from disk and caches it.
/// @details WARNING!  READS FROM DISK!
void
ExternalPackerResultLoader::read_rotamer_selection_file(
	std::string const &filename
) {
	TR << "Reading packer solution from \"" << filename << "\"." << std::endl;
	set_rotamer_selection_string( utility::file_contents(filename) );
}

/// @brief Set the packer problem file contents string.
/// @details Does not read from disk.  Called by read_packer_problem_file().
void
ExternalPackerResultLoader::set_packer_problem_string(
	std::string const &file_contents
) {
	packer_problem_file_data_ = file_contents;
	problem_was_loaded_ = true;
}

/// @brief Set the packer solution file contents string.
/// @details Does not read from disk.  Called by read_rotamer_selection_file().
void
ExternalPackerResultLoader::set_rotamer_selection_string(
	std::string const &file_contents
) {
	rotamer_selection_file_data_ = file_contents;
	solution_was_loaded_ = true;
}

///////////////////////
/// private methods ///
///////////////////////


std::ostream &
operator<<( std::ostream & os, ExternalPackerResultLoader const & mover )
{
	mover.show(os);
	return os;
}

} //protocols
} //quantum_annealing
