// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/farnesyl/InstallFarnesylMover.hh
/// @brief Modifies a free cysteine residue with a branch of 3 DMA residues (the terpene monomer) to create farnesyl-cysteine
/// @author Andy Watkins (andy.watkins2@gmail.com)

#ifndef INCLUDED_protocols_farnesyl_InstallFarnesylMover_HH
#define INCLUDED_protocols_farnesyl_InstallFarnesylMover_HH

// Unit headers
#include <protocols/farnesyl/InstallFarnesylMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover

namespace protocols {
namespace farnesyl {

///@brief Modifies a free cysteine residue with a branch of 3 DMA residues (the terpene monomer) to create farnesyl-cysteine
class InstallFarnesylMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	InstallFarnesylMover();

	/// @brief Copy constructor (not needed unless you need deep copies)
	InstallFarnesylMover( InstallFarnesylMover const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~InstallFarnesylMover() override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

public:
	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	void
	sample_first( core::pose::Pose & pose, Size const cyspos );

	void
	sample_second( core::pose::Pose & pose );

	void
	sample_third( core::pose::Pose & pose );


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

	//InstallFarnesylMover & operator=( InstallFarnesylMover const & src );

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

	void set_selector( core::select::residue_selector::ResidueSelectorCOP const & selector ) {
		std::cout << "In set_selector_ " << std::endl;
		if ( selector ) { std::cout << "passed selector is nonnull " << std::endl; }
		selector_ = selector;
	}

private: // methods

private: // data
	core::select::residue_selector::ResidueSelectorCOP selector_ = nullptr;
	bool sample_per_residue_ = false;

};

std::ostream &
operator<<( std::ostream & os, InstallFarnesylMover const & mover );

} //protocols
} //farnesyl

#endif //protocols_farnesyl_InstallFarnesylMover_HH
