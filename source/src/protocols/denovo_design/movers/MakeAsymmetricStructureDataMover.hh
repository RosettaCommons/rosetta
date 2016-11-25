// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/movers/MakeAsymmetricStructureDataMover.hh
/// @brief Converts a StructureData for a symmetric pose into an asymmetric representation
/// @author Tom Linsky (tlinsky@gmail.com)

#ifndef INCLUDED_protocols_denovo_design_movers_MakeAsymmetricStructureDataMover_hh
#define INCLUDED_protocols_denovo_design_movers_MakeAsymmetricStructureDataMover_hh

// Unit headers
#include <protocols/denovo_design/movers/MakeAsymmetricStructureDataMover.fwd.hh>
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

///@brief Converts a StructureData for a symmetric pose into an asymmetric representation
class MakeAsymmetricStructureDataMover : public protocols::moves::Mover {

public:

	MakeAsymmetricStructureDataMover();

	// copy constructor (not needed unless you need deep copies)
	//MakeAsymmetricStructureDataMover( MakeAsymmetricStructureDataMover const & src );

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~MakeAsymmetricStructureDataMover();

	// XRW TEMP  static std::string
	// XRW TEMP  class_name();

public:
	// mover virtual API
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output = std::cout ) const override;

	// XRW TEMP  virtual std::string
	// XRW TEMP  get_name() const;

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	//MakeAsymmetricStructureDataMover & operator=( MakeAsymmetricStructureDataMover const & src );

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


private:

};

std::ostream &
operator<<( std::ostream & os, MakeAsymmetricStructureDataMover const & mover );

} //protocols
} //denovo_design
} //movers

#endif //protocols/denovo_design/movers_MakeAsymmetricStructureDataMover_hh
