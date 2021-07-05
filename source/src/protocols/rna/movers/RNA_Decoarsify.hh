// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rna/movers/RNA_Decoarsify.hh
/// @brief Make a coarse RNA pose into a fullatom representation.
/// @author Andy Watkins (andy.watkins2@gmail.com)

#ifndef INCLUDED_protocols_rna_movers_RNA_Decoarsify_HH
#define INCLUDED_protocols_rna_movers_RNA_Decoarsify_HH

// Unit headers
#include <protocols/rna/movers/RNA_Decoarsify.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover

namespace protocols {
namespace rna {
namespace movers {

///@brief Make a coarse RNA pose into a fullatom representation.
class RNA_Decoarsify : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	RNA_Decoarsify();

	/// @brief Copy constructor (not needed unless you need deep copies)
	RNA_Decoarsify( RNA_Decoarsify const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~RNA_Decoarsify() override;


public:

	/////////////////////
	/// Mover Methods ///
	/////////////////////


	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output = std::cout ) const override;


public:

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data ) override;

	//RNA_Decoarsify & operator=( RNA_Decoarsify const & src );

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

};

std::ostream &
operator<<( std::ostream & os, RNA_Decoarsify const & mover );

} //movers
} //rna
} //protocols

#endif //protocols_rna_movers_RNA_Decoarsify_HH
