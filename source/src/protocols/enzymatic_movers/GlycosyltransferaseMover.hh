// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/enzymatic_movers/GlycosyltransferaseMover.hh
/// @brief  Declarations and simple accessor/mutator definitions for GlycosyltransferaseMover
/// @author Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_protocols_enzymatic_movers_GlycosyltransferaseMover_HH
#define INCLUDED_protocols_enzymatic_movers_GlycosyltransferaseMover_HH

// Unit headers
#include <protocols/enzymatic_movers/GlycosyltransferaseMover.fwd.hh>
#include <protocols/enzymatic_movers/EnzymaticMover.hh>

// Project headers
#include <core/types.hh>


namespace protocols {
namespace enzymatic_movers {

typedef std::pair< core::uint, std::string > GlycosylationSite;

/// @details  WiP
class GlycosyltransferaseMover : public EnzymaticMover {
public:  // Standard methods //////////////////////////////////////////////////
	/// @brief  Default constructor
	GlycosyltransferaseMover();

	/// @brief  Copy constructor
	GlycosyltransferaseMover( GlycosyltransferaseMover const & object_to_copy );

	// Assignment operator
	GlycosyltransferaseMover & operator=( GlycosyltransferaseMover const & object_to_copy );

	// Destructor
	virtual ~GlycosyltransferaseMover();


public: // Standard Rosetta methods ///////////////////////////////////////////
	// General methods
	/// @brief  Register options with the option system.
	static void register_options();


	// Mover methods
	/// @brief  Return the name of the Mover.

	moves::MoverOP clone() const override;

	moves::MoverOP fresh_instance() const override;

	void parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data,
		Filters_map const & /*filters*/,
		moves::Movers_map const & /*movers*/,
		core::pose::Pose const & pose ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


	/// @brief  Apply the corresponding move to <input_pose>.
	//virtual void apply( core::pose::Pose & input_pose );


protected:
	void perform_reaction( core::pose::Pose & input_pose, core::uint const sepos, std::string const & cosubstrate ) override;


private:  // Private methods //////////////////////////////////////////////////
	// Set command-line options.  (Called by init().)
	void set_commandline_options();

	// Initialize data members from arguments.
	void init();

	// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
	void copy_data(
		GlycosyltransferaseMover & object_to_copy_to, GlycosyltransferaseMover const & object_to_copy_from );
};  // class GlycosyltransferaseMover

// Insertion operator (overloaded so that GlycosyltransferaseMover can be "printed" in PyRosetta).
std::ostream & operator<<( std::ostream & output, GlycosyltransferaseMover const & object_to_output );

}  // namespace enzymatic_movers
}  // namespace protocols

#endif  // INCLUDED_protocols_enzymatic_movers_GlycosyltransferaseMover_HH
