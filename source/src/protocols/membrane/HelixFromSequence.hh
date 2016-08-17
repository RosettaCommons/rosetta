// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/membrane/HelixFromSequence.hh
/// @brief Generates a (transmembrane) helix from a fasta file containing the sequence
/// @author Julia Koehler Leman (julia.koehler1982@gmail.com)


#ifndef INCLUDED_protocols_membrane_HelixFromSequence_hh
#define INCLUDED_protocols_membrane_HelixFromSequence_hh

// Unit headers
#include <protocols/membrane/HelixFromSequence.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>


namespace protocols {
namespace membrane {

///@brief Generates a (transmembrane) helix from a fasta file containing the sequence
class HelixFromSequence : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	HelixFromSequence();

	/// @brief Copy constructor
	HelixFromSequence( HelixFromSequence const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	virtual ~HelixFromSequence();

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Apply the mover
	virtual void
	apply( core::pose::Pose & pose );

	/// @brief Show the contents of the Mover
	virtual void
	show( std::ostream & output=std::cout ) const;

	/// @brief Get the name of the Mover
	std::string
	get_name() const;
	
	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	//HelixFromSequence & operator=( HelixFromSequence const & src );

	/// @brief required in the context of the parser/scripting scheme
	virtual moves::MoverOP
	fresh_instance() const;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const;

private: // methods
	
	/// @brief register options
	void register_options();
	
	/// @brief set defaults
	void set_defaults();

	/// @brief init from commandline
	void init_from_cmd();
	
private: // data

	/// @brief Fasta sequence
	std::string seq_;
	
	/// @brief is it a membrane protein?
	bool mem_;
	
	/// @brief optimize membrane embedding via the scorefunction?
	bool opt_mem_;
	
	/// @brief skip relax run after building the helix?
	bool skip_rlx_;
	
};

std::ostream &operator<< (std::ostream &os, HelixFromSequence const &mover);


} //protocols
} //membrane


#endif //protocols/membrane_HelixFromSequence_hh







