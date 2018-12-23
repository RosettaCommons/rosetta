// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/quantum_annealing/ExternalPackerResultLoader.hh
/// @brief A mover that reads the result of an externally-called packer, and applies it to a pose to generate a structure.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_quantum_annealing_ExternalPackerResultLoader_HH
#define INCLUDED_protocols_quantum_annealing_ExternalPackerResultLoader_HH

// Unit headers
#include <protocols/quantum_annealing/ExternalPackerResultLoader.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover

namespace protocols {
namespace quantum_annealing {

///@brief A mover that reads the result of an externally-called packer, and applies it to a pose to generate a structure.
class ExternalPackerResultLoader : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	ExternalPackerResultLoader();

	/// @brief Copy constructor (not needed unless you need deep copies)
	ExternalPackerResultLoader( ExternalPackerResultLoader const & /*src*/ ) = default;

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~ExternalPackerResultLoader() override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

public:
	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output = std::cout ) const override;

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	/// @details Warning!  This mover reads from disk at parse time!
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	//ExternalPackerResultLoader & operator=( ExternalPackerResultLoader const & src );

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

public: // methods

	/// @brief Reads a description of the packer problem from disk and caches it.
	/// @details WARNING!  READS FROM DISK!
	void read_packer_problem_file( std::string const &filename );

	/// @brief Reads a description of the packer solution from disk and caches it.
	/// @details WARNING!  READS FROM DISK!
	void read_rotamer_selection_file( std::string const &filename );

	/// @brief Set the packer problem file contents string.
	/// @details Does not read from disk.  Called by read_packer_problem_file().
	void set_packer_problem_string( std::string const &file_contents );

	/// @brief Set the packer solution file contents string.
	/// @details Does not read from disk.  Called by read_rotamer_selection_file().
	void set_rotamer_selection_string( std::string const &file_contents );


private: // methods

private: // data

	/// @brief The contents of a data file (or stream) indicating a selection of 1 rotamer at each packable position.
	std::string rotamer_selection_file_data_;

	/// @brief The contents of a data file (or stream) describing the packer problem.
	/// @details Includes full information needed to rebuild the pose.
	std::string packer_problem_file_data_;

	/// @brief Has the problem been loaded from disk?
	bool problem_was_loaded_;

	/// @brief Has the solution been loaded from disk?
	bool solution_was_loaded_;

};

std::ostream &
operator<<( std::ostream & os, ExternalPackerResultLoader const & mover );

} //protocols
} //quantum_annealing

#endif //protocols_quantum_annealing_ExternalPackerResultLoader_HH
