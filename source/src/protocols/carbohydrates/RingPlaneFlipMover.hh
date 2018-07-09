// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    protocols/carbohydrates/RingPlaneFlipMover.hh
/// @brief   Declarations and simple accessor/mutator definitions for RingPlaneFlipMover.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_protocols_carbohydrates_RingPlaneFlipMover_HH
#define INCLUDED_protocols_carbohydrates_RingPlaneFlipMover_HH

// Unit headers
#include <protocols/carbohydrates/RingPlaneFlipMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/id/TorsionID.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/select/movemap/MoveMapFactory.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>

// Utility Headers
#include <utility/vector1.hh>


namespace protocols {
namespace carbohydrates {

/// @brief    A Mover class for flipping the plane of a carbohydrate pyranose ring 180 degrees about its anomeric bond.
/// @details  Based on a given ResidueSelector and limited by a MoveMap, this mover selects applicable cyclic residues
/// and performs a 180-degree shearing move in which the anomeric bond and the main-chain bond on the opposite side of
/// the ring are moved in opposite directions.  An "applicable" residue is limited to 1,4-linked aldopyranoses or 2,5-
/// linked ketopyranoses for which both the anomeric bond and the glycosidic linkage bond are equatorial.
class RingPlaneFlipMover: public moves::Mover {
public:  // Standard methods //////////////////////////////////////////////////
	/// @brief  Default constructor
	RingPlaneFlipMover();

	/// @brief  Copy constructor
	RingPlaneFlipMover( RingPlaneFlipMover const & object_to_copy );

	/// @brief  Constructor with MoveMap input option
	RingPlaneFlipMover( core::kinematics::MoveMapOP input_movemap );

	/// @brief  Constructor with ResidueSelector input option
	RingPlaneFlipMover( core::select::residue_selector::ResidueSelectorCOP selector );

	// Assignment operator
	RingPlaneFlipMover & operator=( RingPlaneFlipMover const & object_to_copy );

	// Destructor
	~RingPlaneFlipMover() override = default;


public: // Standard Rosetta methods ///////////////////////////////////////////
	// General methods
	/// @brief  Generate string representation of RingPlaneFlipMover for debugging purposes.
	void show( std::ostream & output=std::cout ) const override;


	// Mover methods
	/// @brief  Register options with the option system.
	static void register_options();

	/// @brief  Return the name of the Mover.
	std::string get_name() const override;

	static std::string mover_name() { return "RingPlaneFlipMover"; }


	protocols::moves::MoverOP clone() const override;

	protocols::moves::MoverOP fresh_instance() const override;


	void parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data,
		Filters_map const & /*filters*/,
		moves::Movers_map const & /*movers*/,
		Pose const & /*pose*/ ) override;

	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


	/// @brief  Apply the corresponding move to <input_pose>.
	void apply( core::pose::Pose & input_pose ) override;


public: // Accessors/Mutators /////////////////////////////////////////////////
	/// @brief  Get the current MoveMap, creating it if needed.
	core::kinematics::MoveMapCOP movemap( core::pose::Pose const & pose ) const;

	/// @brief  Get the current ResidueSelector.
	core::select::residue_selector::ResidueSelectorCOP selector() const { return selector_; }


	/// @brief  Set the MoveMap.
	void movemap( core::kinematics::MoveMapOP new_movemap ) { movemap_ = new_movemap; }

	/// @brief  Set the MoveMapFactory.
	void movemap_factory( core::select::movemap::MoveMapFactoryCOP new_movemap_factory ) {
		movemap_factory_ = new_movemap_factory;
	}

	/// @brief  Set the ResidueSelector.
	void selector( core::select::residue_selector::ResidueSelectorCOP new_selector ) { selector_ = new_selector; }


private:  // Private methods //////////////////////////////////////////////////
	// Initialize data members.
	void init();

	// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
	void copy_data( RingPlaneFlipMover & object_to_copy_to, RingPlaneFlipMover const & object_to_copy_from);


	// Setup list of movable residues from MoveMap and/or ResidueSelector.
	void setup_movable_torsion_pairs( core::pose::Pose const & pose );


private:  // Private data /////////////////////////////////////////////////////
	core::kinematics::MoveMapOP mutable movemap_;
	core::select::movemap::MoveMapFactoryCOP movemap_factory_;
	core::select::residue_selector::ResidueSelectorCOP selector_;

	utility::vector1< std::pair< core::id::TorsionID, core::id::TorsionID > > movable_torsion_pairs_;
};  // class RingPlaneFlipMover


// Insertion operator (overloaded so that RingPlaneFlipMover can be "printed" in PyRosetta).
std::ostream & operator<<( std::ostream & output, RingPlaneFlipMover const & object_to_output );

}  // namespace carbohydrates
}  // namespace protocols

#endif  // INCLUDED_protocols_carbohydrates_RingPlaneFlipMover_HH
