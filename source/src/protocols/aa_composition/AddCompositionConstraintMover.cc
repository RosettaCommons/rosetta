// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/aa_composition/AddCompositionConstraintMover.cc
/// @brief Assigns an AACompositionConstraint to a pose, initializing it from a file.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#include <protocols/aa_composition/AddCompositionConstraintMover.hh>
#include <protocols/aa_composition/AddCompositionConstraintMoverCreator.hh>


#include <core/pose/Pose.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/scoring/aa_composition_energy/AACompositionConstraint.hh>

//Auto Headers

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.aa_composition.AddCompositionConstraintMover" );

namespace protocols {
namespace aa_composition {

using namespace core;
using namespace core::scoring;
using namespace constraints;
using namespace utility::tag;

// XRW TEMP std::string
// XRW TEMP AddCompositionConstraintMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return AddCompositionConstraintMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP AddCompositionConstraintMoverCreator::create_mover() const
// XRW TEMP {
// XRW TEMP  return protocols::moves::MoverOP( new AddCompositionConstraintMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP AddCompositionConstraintMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "AddCompositionConstraintMover";
// XRW TEMP }

/// @brief Default Constructor
///
AddCompositionConstraintMover::AddCompositionConstraintMover():
	protocols::moves::Mover( AddCompositionConstraintMover::mover_name() ),
	//TODO initialize variables here
	constraint_()
{
}

/// @brief Copy Constructor
///
AddCompositionConstraintMover::AddCompositionConstraintMover( AddCompositionConstraintMover const &src ):
	protocols::moves::Mover( AddCompositionConstraintMover::mover_name() ),
	//TODO initialize variables here
	constraint_( src.constraint_ )
{
}


/// @brief Destructor
///
AddCompositionConstraintMover::~AddCompositionConstraintMover()= default;

/// @brief Copy this object and return a pointer to the copy.
///
protocols::moves::MoverOP AddCompositionConstraintMover::clone() const { return protocols::moves::MoverOP( new protocols::aa_composition::AddCompositionConstraintMover( *this ) ); }

/// @brief Create a new object of this type and return a pointer to it.
///
protocols::moves::MoverOP AddCompositionConstraintMover::fresh_instance() const { return protocols::moves::MoverOP( new AddCompositionConstraintMover ); }

/// @brief Returns the name of this mover ("AddCompositionConstraintMover").
///
// XRW TEMP std::string
// XRW TEMP AddCompositionConstraintMover::get_name() const {
// XRW TEMP  return AddCompositionConstraintMover::mover_name();
// XRW TEMP }

/// @brief Actually apply the mover to a pose.
///
void
AddCompositionConstraintMover::apply( Pose &pose )
{
	runtime_assert_string_msg( constraint(), "Error in protocols::aa_composition::AddCompositionConstraintMover::apply(): The AACompositionConstraint object was not initialized before the apply() function was called." );
	pose.add_constraint( utility::pointer::dynamic_pointer_cast < core::scoring::constraints::Constraint const >( constraint() ) );
	return;
}

/// @brief Parse RosettaScripts XML tag to set up the mover.
///
void
AddCompositionConstraintMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &datamap,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const &
)
{
	runtime_assert_string_msg( tag->hasOption("filename"), "Error in protocols::aa_composition::AddCompositionConstraintMover::parse_my_tag(): A \"filename\" option MUST be provided." );
	std::string const compfile( tag->getOption<std::string>("filename") );
	if ( TR.visible() ) TR << "Set filename to " << compfile << "." << std::endl;
	create_constraint_from_file( compfile );

	if ( tag->hasOption("selector") ) {
		std::string const selector_name ( tag->getOption< std::string >( "selector" ) );
		if ( TR.visible() ) TR << "Set selector name to " << selector_name << "." << std::endl;
		core::select::residue_selector::ResidueSelectorCOP residue_selector;
		try {
			residue_selector = datamap.get_ptr< core::select::residue_selector::ResidueSelector const >( "ResidueSelector", selector_name );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::string error_message = "Failed to find ResidueSelector named '" + selector_name + "' from the Datamap from AddCompositionConstraintMover::parse_tag()\n" + e.msg();
			throw utility::excn::EXCN_Msg_Exception( error_message );
		}
		runtime_assert( residue_selector );
		add_residue_selector( residue_selector );
	}

	if ( TR.visible() ) TR.flush();
	return;
}

/// @brief Create the AACompositionConstraint object and initialize it from a .comp file.
///
void
AddCompositionConstraintMover::create_constraint_from_file( std::string const &filename ) {
	runtime_assert_string_msg( !constraint_, "Error in protocols::aa_composition::AddCompositionConstraintMover::create_constraint_from_file():  The constraint object already has been created!" );
	constraint_ = core::scoring::aa_composition_energy::AACompositionConstraintOP( new core::scoring::aa_composition_energy::AACompositionConstraint() );
	constraint_->initialize_from_file( filename );
	if ( TR.visible() ) {
		TR << "Initialized AACompositionConstraint object from file " << filename << "." << std::endl;
		TR.flush();
	}
	return;
}

/// @brief Create the AACompositionConstraint object from the data from a .comp file.
/// @details Allows external code to create the constraint object without having it read directly from disk.
void
AddCompositionConstraintMover::create_constraint_from_file_contents(
	std::string const &filecontents
) {
	runtime_assert_string_msg( !constraint_, "Error in protocols::aa_composition::AddCompositionConstraintMover::create_constraint_from_filecontents():  The constraint object already has been created!" );
	constraint_ = core::scoring::aa_composition_energy::AACompositionConstraintOP( new core::scoring::aa_composition_energy::AACompositionConstraint() );
	constraint_->initialize_from_file_contents( filecontents );
	if ( TR.visible() ) {
		TR << "Initialized AACompositionConstraint object from file contents:\n" << filecontents << std::endl;
		TR.flush();
	}
	return;
}

/// @brief Add a ResidueSelector to the constraint to use as a mask.
/// @details The constraint must already have been created with the create_constraint_from_file() function before this function is called.
void
AddCompositionConstraintMover::add_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector_in ) {
	runtime_assert_string_msg( constraint_, "Error in protocols::aa_composition::AddCompositionConstraintMover::add_residue_selector(): The AACompositionConstraint object was not initialized before the add_residue_selector() function was called." );
	constraint_->set_selector( selector_in );
	return;
}

std::string AddCompositionConstraintMover::get_name() const {
	return mover_name();
}

std::string AddCompositionConstraintMover::mover_name() {
	return "AddCompositionConstraintMover";
}

void AddCompositionConstraintMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute("filename", xs_string, "Name of composition constraint file" );
	// I would consider a complex type that ends in .pdb, but I want users to be
	// able to provide cifs! Or name their files whatever!
	attlist + XMLSchemaAttribute("selector", xs_string, "Residue selector named somewhere else in the script" );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Add composition constraints from the provided file to the selected region", attlist );
}

std::string AddCompositionConstraintMoverCreator::keyname() const {
	return AddCompositionConstraintMover::mover_name();
}

protocols::moves::MoverOP
AddCompositionConstraintMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new AddCompositionConstraintMover );
}

void AddCompositionConstraintMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AddCompositionConstraintMover::provide_xml_schema( xsd );
}


} // aa_composition
} // protocols
