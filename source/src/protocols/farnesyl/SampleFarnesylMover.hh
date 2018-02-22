// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/farnesyl/SampleFarnesylMover.hh
/// @brief Modifies a free cysteine residue with a branch of 3 DMA residues (the terpene monomer) to create farnesyl-cysteine
/// @author Andy Watkins (andy.watkins2@gmail.com)

#ifndef INCLUDED_protocols_farnesyl_SampleFarnesylMover_HH
#define INCLUDED_protocols_farnesyl_SampleFarnesylMover_HH

// Unit headers
#include <protocols/farnesyl/SampleFarnesylMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover

namespace protocols {
namespace farnesyl {

///@brief Modifies a free cysteine residue with a branch of 3 DMA residues (the terpene monomer) to create farnesyl-cysteine
class SampleFarnesylMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	SampleFarnesylMover();

	/// @brief Copy constructor (not needed unless you need deep copies)
	SampleFarnesylMover( SampleFarnesylMover const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~SampleFarnesylMover() override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

public:
	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	void sample_farnesyl(
		core::pose::Pose & pose,
		Size const cys_idx,
		Size const dma_one_idx,
		Size const dma_two_idx,
		Size const dma_three_idx );

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

	//SampleFarnesylMover & operator=( SampleFarnesylMover const & src );

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

	void set_score_function( core::scoring::ScoreFunctionOP const & sfxn ) { sfxn_ = sfxn; }
	void set_enumerate( bool const enumerate ) { enumerate_ = enumerate; }

private: // methods

private: // data
	core::scoring::ScoreFunctionOP sfxn_ = nullptr;
	bool enumerate_ = true;

};

std::ostream &
operator<<( std::ostream & os, SampleFarnesylMover const & mover );

} //protocols
} //farnesyl

#endif //protocols_farnesyl_SampleFarnesylMover_HH
