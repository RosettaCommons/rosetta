// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    protocols/carbohydrates/TautomerizeAnomerMover.hh
/// @brief   Declarations and simple accessor/mutator definitions for TautomerizeAnomerMover.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_protocols_carbohydrates_TautomerizeAnomerMover_HH
#define INCLUDED_protocols_carbohydrates_TautomerizeAnomerMover_HH

// Unit headers
#include <protocols/carbohydrates/TautomerizeAnomerMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Utility Headers
#include <utility/vector1.hh>


namespace protocols {
namespace carbohydrates {

/// @brief    A Mover class for tautomerizing from one anomer to another at a reducing end.
/// @details  This carbohydrate-specific Mover randomly selects a free reducing end (not a glycoside) and inverts the
/// stereochemistry, swapping alpha anomers for beta and beta for alpha.  (This could be considered an extremely
/// limited design case; however, reducing ends readily tautomerize in solution, in contrast to other cases, in which
/// residues do not readily mutate into others!)  It is generally not certain which form is preferred (if any) in
/// sugar-binding proteins, and crystal structures sometimes arbitrarily assign one anomer over another when fitting
/// density, so this Mover can assure that each anomer is sampled.
/// If a ResidueSelector is set, the Mover will select from the subset at random; it will not guarantee
/// tautomerization of every Residue in the subset.
/// @note     This Mover assumes that, any time a new saccharide .params file is added to the Rosetta database, its
/// anomer is also added.
class TautomerizeAnomerMover: public moves::Mover {
public:  // Standard methods //////////////////////////////////////////////////
	/// @brief  Default constructor
	TautomerizeAnomerMover();

	/// @brief  Copy constructor
	TautomerizeAnomerMover( TautomerizeAnomerMover const & object_to_copy );

	/// @brief  Constructor with ResidueSelector input option
	TautomerizeAnomerMover( core::select::residue_selector::ResidueSelectorCOP selector );

	// Assignment operator
	TautomerizeAnomerMover & operator=( TautomerizeAnomerMover const & object_to_copy );

	// Destructor
	~TautomerizeAnomerMover() override = default;


public: // Standard Rosetta methods ///////////////////////////////////////////
	// General methods
	/// @brief  Generate string representation of TautomerizeAnomerMover for debugging purposes.
	void show( std::ostream & output=std::cout ) const override;


	// Mover methods
	/// @brief  Register options with the option system.
	static void register_options();

	/// @brief  Return the name of the Mover.
	std::string get_name() const override;

	static std::string mover_name() { return "TautomerizeAnomerMover"; }


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
	/// @brief  Get the current ResidueSelector.
	core::select::residue_selector::ResidueSelectorCOP selector() const { return selector_; }

	/// @brief  Set the ResidueSelector.
	void selector( core::select::residue_selector::ResidueSelectorCOP new_selector ) { selector_ = new_selector; }


private:  // Private methods //////////////////////////////////////////////////
	// Initialize data members.
	void init();

	// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
	void copy_data( TautomerizeAnomerMover & object_to_copy_to, TautomerizeAnomerMover const & object_to_copy_from);


	// Setup list of movable residues from the ResidueSelector.
	void setup_movable_reducing_ends( core::pose::Pose const & pose );


private:  // Private data /////////////////////////////////////////////////////
	core::select::residue_selector::ResidueSelectorCOP selector_;

	utility::vector1< core::uint > movable_reducing_ends_;
};  // class TautomerizeAnomerMover


// Insertion operator (overloaded so that TautomerizeAnomerMover can be "printed" in PyRosetta).
std::ostream & operator<<( std::ostream & output, TautomerizeAnomerMover const & object_to_output );

}  // namespace carbohydrates
}  // namespace protocols

#endif  // INCLUDED_protocols_carbohydrates_TautomerizeAnomerMover_HH
