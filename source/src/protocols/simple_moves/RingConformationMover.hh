// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    protocols/simple_moves/RingConformationMover.hh
/// @brief   Declarations and simple accessor/mutator definitions for RingConformationMover.
/// @author  Labonte <JWLabonte@jhu.edu>

#ifndef INCLUDED_protocols_simple_moves_RingConformationMover_HH
#define INCLUDED_protocols_simple_moves_RingConformationMover_HH

// Unit headers
#include <protocols/simple_moves/RingConformationMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/types.hh>
#include <core/id/types.hh>
#include <core/kinematics/MoveMap.fwd.hh>

// Utility headers
#include <utility/vector1.fwd.hh>

// C++ headers
#include <string>
#include <iostream>


namespace protocols {
namespace simple_moves {

/// @details  Based on a given MoveMap, this mover selects movable cyclic residues and flips their rings to an
/// idealized ring conformer.
class RingConformationMover: public moves::Mover {
public:  // Standard methods //////////////////////////////////////////////////
	/// @brief  Default constructor
	RingConformationMover();

	/// @brief  Copy constructor
	RingConformationMover( RingConformationMover const & object_to_copy );

	/// @brief  Constructor with MoveMap input option
	RingConformationMover( core::kinematics::MoveMapOP input_movemap );

	// Assignment operator
	RingConformationMover & operator=( RingConformationMover const & object_to_copy );

	// Destructor
	~RingConformationMover() override;


public: // Standard Rosetta methods ///////////////////////////////////////////
	// General methods
	/// @brief  Register options with the option system.
	static void register_options();

	/// @brief  Generate string representation of RingConformationMover for debugging purposes.
	void show( std::ostream & output=std::cout ) const override;


	// Mover methods
	/// @brief  Return the name of the Mover.

	protocols::moves::MoverOP clone() const override;

	protocols::moves::MoverOP fresh_instance() const override;

	void parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data,
		Filters_map const & /*filters*/,
		moves::Movers_map const & /*movers*/,
		Pose const & pose ) override;

	/// @brief  Apply the corresponding move to <input_pose>.
	void apply( core::pose::Pose & input_pose ) override;


public: // Accessors/Mutators /////////////////////////////////////////////////
	/// @brief  Get the current MoveMap.
	core::kinematics::MoveMapCOP movemap() const;

	/// @brief  Set the MoveMap.
	void movemap( core::kinematics::MoveMapOP new_movemap );

	/// @brief  Get whether or not this Mover will sample all ring conformers, regardless of energy.
	bool sample_all_conformers() const { return sample_all_conformers_; }

	/// @brief  Set whether or not this Mover will sample all ring conformers, regardless of energy.
	void sample_all_conformers( bool setting ) { sample_all_conformers_ = setting; }

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );



private:  // Private methods //////////////////////////////////////////////////
	// Set command-line options.  (Called by init())
	void set_commandline_options();

	// Initialize data members from arguments.
	void init( core::kinematics::MoveMapOP movemap );

	// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
	void copy_data( RingConformationMover & object_to_copy_to, RingConformationMover const & object_to_copy_from);

	// Setup list of movable cyclic residues from MoveMap.
	void setup_residue_list( core::pose::Pose & pose );


private:  // Private data /////////////////////////////////////////////////////
	core::kinematics::MoveMapOP movemap_;
	utility::vector1<core::Size> residue_list_;  // list of movable cyclic residues by residue number
	bool locked_;  // Is this mover locked from the command line?
	bool sample_all_conformers_;  // Does this Mover sample both energy maxima and minima among ring conformers?

};  // class RingConformationMover

// Insertion operator (overloaded so that RingConformationMover can be "printed" in PyRosetta).
std::ostream & operator<<( std::ostream & output, RingConformationMover const & object_to_output );

}  // namespace simple_moves
}  // namespace protocols

#endif  // INCLUDED_simple_moves_protocols_RingConformationMover_HH
