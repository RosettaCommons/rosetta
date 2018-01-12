// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/carbohydrates/GlycanTreeMinMover.hh
/// @brief A class that selects the downstream branch from residues in a movemap/selector, and minimizes those residues if on in the primary glycan movemap. Multiple Applies randomly select a different residue in the movemap/selector
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_carbohydrates_GlycanTreeMinMover_hh
#define INCLUDED_protocols_carbohydrates_GlycanTreeMinMover_hh

// Unit headers
#include <protocols/carbohydrates/GlycanTreeMinMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/minimization_packing/MinMover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace carbohydrates {

///@brief A class that selects the downstream branch from residues in a movemap/selector, and minimizes those residues if on in the primary glycan movemap. Multiple Applies randomly select a different residue in the movemap/selector
///
/// Supports Symmetry
///
///@details
class GlycanTreeMinMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	GlycanTreeMinMover();

	GlycanTreeMinMover(
		core::select::residue_selector::ResidueSelectorCOP selector,
		bool min_rings = false,
		bool min_bb = true,
		bool min_chi = true
	);

	/// @brief Copy constructor (not needed unless you need deep copies)
	GlycanTreeMinMover( GlycanTreeMinMover const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	virtual ~GlycanTreeMinMover();

	/////////////////////
	/// Mover Methods ///
	/////////////////////

public:

	/// @brief Apply the mover
	virtual void
	apply( core::pose::Pose & pose ) override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

public:

	void
	set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector);

	///@brief Minimize Rings? Default False
	void
	set_min_rings( bool min_rings );

	///@brief Minimize Chi? Default True
	void
	set_min_chi( bool min_chi );

	///@brief Minimize BB? Default True
	void
	set_min_bb( bool min_bb );

	///@brief Set a pre-configured MinMover for this class.
	///  Will OVERRIDE movemap settings.
	void
	set_minmover( protocols::minimization_packing::MinMoverCOP min_mover);

public:

	/// @brief Show the contents of the Mover
	static std::string
	class_name();

	/// @brief Get the name of the Mover
	virtual std::string
	get_name() const override;

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	//GlycanTreeMinMover & operator=( GlycanTreeMinMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	clone() const override;

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private: // methods

private:

	bool min_rings_ = false;
	bool min_bb_ = true;
	bool min_chi_ = true;

	core::select::residue_selector::ResidueSelectorCOP selector_ = nullptr;
	protocols::minimization_packing::MinMoverOP min_mover_ = nullptr;

};


} //protocols
} //carbohydrates

#endif //protocols/carbohydrates_GlycanTreeMinMover_hh
