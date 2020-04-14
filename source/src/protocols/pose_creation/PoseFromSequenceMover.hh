// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_creation/PoseFromSequenceMover.hh
/// @brief A class for generating a pose from a sequence/fasta
/// @author Dan Farrell (danpf@uw.edu)

#ifndef INCLUDED_protocols_pose_creation_PoseFromSequenceMover_HH
#define INCLUDED_protocols_pose_creation_PoseFromSequenceMover_HH

// Unit headers
#include <protocols/pose_creation/PoseFromSequenceMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace pose_creation {

///@brief A class for generating a pose from a sequence/fasta
class PoseFromSequenceMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	PoseFromSequenceMover();

	/// @brief Default constructor
	PoseFromSequenceMover(
		std::string const & sequence, std::string const & residue_type_set = "fa_standard",
		bool const extended = false);

	/// @brief Copy constructor (not needed unless you need deep copies)
	PoseFromSequenceMover( PoseFromSequenceMover const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~PoseFromSequenceMover() override;


public:

	/////////////////////
	/// Mover Methods ///
	/////////////////////


	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output = std::cout ) const override;

	void
	set_sequence(std::string const &);

	std::string const &
	get_sequence() const;

	void
	set_residue_type_set(std::string const &);

	std::string const &
	get_residue_type_set() const;

	void
	set_extended(bool const);

	bool
	get_extended() const;


public:

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) override;

	PoseFromSequenceMover & operator=( PoseFromSequenceMover const & src );

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

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

private: // data
	std::string sequence_;
	std::string residue_type_set_;
	bool extended_;
};

std::ostream &
operator<<( std::ostream & os, PoseFromSequenceMover const & mover );

} //pose_creation
} //protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_pose_creation_PoseFromSequenceMover )
#endif // SERIALIZATION

#endif //protocols_pose_creation_PoseFromSequenceMover_HH
