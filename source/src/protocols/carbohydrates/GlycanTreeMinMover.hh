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
#include <protocols/simple_moves/MinMover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace carbohydrates {

///@brief A class that selects the downstream branch from residues in a movemap/selector, and minimizes those residues if on in the primary glycan movemap. Multiple Applies randomly select a different residue in the movemap/selector
///
///@details
class GlycanTreeMinMover : public protocols::simple_moves::MinMover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	GlycanTreeMinMover();
	
	GlycanTreeMinMover(
		core::kinematics::MoveMapOP movemap_in,
		ScoreFunctionCOP scorefxn_in,
		std::string const & min_type_in = "dfpmin_armijo_nonmonotone",
		Real tolerance_in = .01
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

	void
	set_movemap( core::kinematics::MoveMapCOP movemap_in) override;
	
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
	utility::vector1< core::Size > mm_residues_;
};


} //protocols
} //carbohydrates

#endif //protocols/carbohydrates_GlycanTreeMinMover_hh
