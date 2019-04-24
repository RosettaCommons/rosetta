// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/enzymatic_movers/KinaseMover.hh
/// @brief  Declarations and simple accessor/mutator definitions for KinaseMover
/// @author Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_protocols_enzymatic_movers_KinaseMover_HH
#define INCLUDED_protocols_enzymatic_movers_KinaseMover_HH

// Unit headers
#include <protocols/enzymatic_movers/KinaseMover.fwd.hh>
#include <protocols/enzymatic_movers/EnzymaticMover.hh>

// Project headers
#include <core/types.hh>


namespace protocols {
namespace enzymatic_movers {

/// @details  This Mover simulates the activity of a virtual kinase enzyme on a Pose by phosphorylating a residue at a
///           biologically relevant sequon position.
class KinaseMover : public EnzymaticMover {
public:  // Standard methods //////////////////////////////////////////////////
	/// @brief  Default constructor
	KinaseMover();

	/// @brief  Copy constructor
	KinaseMover( KinaseMover const & object_to_copy );

	// Assignment operator
	KinaseMover & operator=( KinaseMover const & object_to_copy );

	// Destructor
	~KinaseMover() override = default;


public: // Standard Rosetta methods ///////////////////////////////////////////
	// General methods
	/// @brief  Register options with the option system.
	static void register_options();


	// Mover methods
	/// @brief  Return the name of the Mover.
	std::string get_name() const override;

	static std::string mover_name() { return "KinaseMover"; }


	moves::MoverOP clone() const override;

	moves::MoverOP fresh_instance() const override;


	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


protected:
	void perform_reaction(
		core::pose::Pose & input_pose,
		core::uint const site,
		std::string const & cosubstrate ) override;
};  // class KinaseMover

// Insertion operator (overloaded so that KinaseMover can be "printed" in PyRosetta).
std::ostream & operator<<( std::ostream & output, KinaseMover const & object_to_output );

}  // namespace enzymatic_movers
}  // namespace protocols

#endif  // INCLUDED_protocols_enzymatic_movers_KinaseMover_HH
