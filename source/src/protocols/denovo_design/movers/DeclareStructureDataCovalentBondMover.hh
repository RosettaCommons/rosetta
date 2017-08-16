// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/movers/DeclareStructureDataCovalentBondMover.hh
/// @brief declares covalent bond
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_movers_DeclareStructureDataCovalentBondMover_hh
#define INCLUDED_protocols_denovo_design_movers_DeclareStructureDataCovalentBondMover_hh

// Unit headers
#include <protocols/denovo_design/movers/DeclareStructureDataCovalentBondMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace denovo_design {
namespace movers {

///@brief declares covalent bond
class DeclareStructureDataCovalentBondMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	DeclareStructureDataCovalentBondMover();

	/// @brief Copy constructor (not needed unless you need deep copies)
	// DeclareStructureDataCovalentBondMover( DeclareStructureDataCovalentBondMover const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	virtual ~DeclareStructureDataCovalentBondMover();

	/////////////////////
	/// Mover Methods ///
	/////////////////////

public:
	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	/// @brief Show the contents of the Mover

	void
	show( std::ostream & output = std::cout ) const override;

	/// @brief Get the name of the Mover

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

	//DeclareStructureDataCovalentBondMover & operator=( DeclareStructureDataCovalentBondMover const & src );

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

private: // data
	std::string atom1_;
	std::string atom2_;

};

std::ostream &
operator<<( std::ostream & os, DeclareStructureDataCovalentBondMover const & mover );

} //protocols
} //denovo_design
} //movers

#endif //protocols/denovo_design/movers_DeclareStructureDataCovalentBondMover_hh
