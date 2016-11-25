// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/CopyRotamerMover.hh
/// @brief A mover to copy a rotamer (residue identity and conformation) from one position in a pose to another.
/// @author Vikram K. Mulligan (vmullig@uw.edu), Baker laboratory.

#ifndef INCLUDED_protocols_simple_moves_CopyRotamerMover_hh
#define INCLUDED_protocols_simple_moves_CopyRotamerMover_hh

// Unit headers
#include <protocols/simple_moves/CopyRotamerMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace simple_moves {

///@brief A mover to copy a rotamer (residue identity and conformation) from one position in a pose to another.
class CopyRotamerMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	CopyRotamerMover();

	/// @brief Copy constructor (not needed unless you need deep copies)
	CopyRotamerMover( CopyRotamerMover const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	virtual ~CopyRotamerMover();

	/////////////////////
	/// Mover Methods ///
	/////////////////////

public:

	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	/// @brief Show the contents of the Mover
	// XRW TEMP  static std::string
	// XRW TEMP  class_name();

	void
	show( std::ostream & output = std::cout ) const override;

	/// @brief Get the name of the Mover
	// XRW TEMP  virtual std::string
	// XRW TEMP  get_name() const;

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	) override;

	//CopyRotamerMover & operator=( CopyRotamerMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	///////////////////////////////
	/// Setters          ///
	///////////////////////////////

	/// @brief Set pose index of the residue FROM which we're copying.
	///
	void set_template_res_index( core::Size const index_in );

	/// @brief Set pose index of the residue TO which we're copying.
	///
	void set_target_res_index( core::Size const index_in );

	/// @brief Set whether we're copying the residue identity.
	///
	void set_copy_identity( bool const setting_in );

	/// @brief Set whether we're copying the residue side-chain torsions.
	///
	void set_copy_torsions( bool const setting_in );

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


	// / @brief Set whether we're copying the residue side-chain bond angles.
	// /
	//void set_copy_bondangles( bool const setting_in );

	// / @brief Set whether we're copying the residue side-chain bond lengths.
	// /
	//void set_copy_bondlengths( bool const setting_in );

private: // methods

private: // data

	/// @brief Pose index of the residue FROM which we're copying.
	///
	core::Size template_res_index_;

	/// @brief Pose index of the residue TO which we're copying.
	///
	core::Size target_res_index_;

	/// @brief Are we copying the residue identity?
	/// @details Default true.
	bool copy_identity_;

	/// @brief Are we copying side-chain torsions?
	/// @details Default true.
	bool copy_torsions_;

	// / @brief Are we copying side-chain bond angles?
	// / @details Default true.
	//bool copy_bondangles_;

	// / @brief Are we copying side-chain bond lengths?
	// / @details Default true.
	//bool copy_bondlengths_;

};

std::ostream &
operator<<( std::ostream & os, CopyRotamerMover const & mover );

} //protocols
} //simple_moves

#endif //protocols/simple_moves_CopyRotamerMover_hh
