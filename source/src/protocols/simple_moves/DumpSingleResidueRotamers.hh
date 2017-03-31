// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/DumpSingleResidueRotamers.hh
/// @brief Given a residue index, dump all of the rotamers to individual PDB files within 0-1 sd of the mean
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_simple_moves_DumpSingleResidueRotamers_HH
#define INCLUDED_protocols_simple_moves_DumpSingleResidueRotamers_HH

// Unit headers
#include <protocols/simple_moves/DumpSingleResidueRotamers.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace simple_moves {

///@brief Given a residue index, dump all of the rotamers to individual PDB files within 0-1 sd of the mean
class DumpSingleResidueRotamers : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	DumpSingleResidueRotamers();

	/// @brief Custom constructor
	DumpSingleResidueRotamers( core::Size rsd_index );

	/// @brief Copy constructor (not needed unless you need deep copies)
	DumpSingleResidueRotamers( DumpSingleResidueRotamers const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~DumpSingleResidueRotamers() override;

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
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

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

private: // methods

	void
	init_from_options();

	void
	write_poses_to_pdbs(
		utility::vector1< core::pose::PoseOP > poses,
		core::Size position );

	utility::vector1< core::pose::PoseOP >
	enumerate_aa_rotamer(
		core::pack::dunbrack::SingleResidueDunbrackLibraryCOP rsd_rl,
		core::pose::Pose & pose,
		core::Size rsd_index ) const;

private: // data

	core::Size rsd_index_;
	std::string prefix_;

	bool write_rotamers_to_pdbs_;
	bool all_positions_;

};

std::ostream &
operator<<( std::ostream & os, DumpSingleResidueRotamers const & mover );

} // simple_moves
} // protocols

#endif // protocols_simple_moves_DumpSingleResidueRotamers_HH
