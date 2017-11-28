// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/aa_composition/AddNetChargeConstraintMover.cc
/// @brief Assigns an NetChargeConstraint to a pose, initializing it from a file.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#include <protocols/aa_composition/AddNetChargeConstraintMover.hh>
#include <protocols/aa_composition/AddNetChargeConstraintMoverCreator.hh>


#include <core/pose/Pose.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/scoring/netcharge_energy/NetChargeConstraint.hh>

//Auto Headers

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.aa_composition.AddNetChargeConstraintMover" );

namespace protocols {
namespace aa_composition {

using namespace core;
using namespace core::scoring;
using namespace constraints;
using namespace utility::tag;

/// @brief Default Constructor
///
AddNetChargeConstraintMover::AddNetChargeConstraintMover():
	protocols::moves::Mover( AddNetChargeConstraintMover::mover_name() ),
	//TODO initialize variables here
	constraint_()
{
}

/// @brief Copy Constructor
///
AddNetChargeConstraintMover::AddNetChargeConstraintMover( AddNetChargeConstraintMover const &src ):
	protocols::moves::Mover( AddNetChargeConstraintMover::mover_name() ),
	//TODO initialize variables here
	constraint_( src.constraint_ )
{
}

/// @brief Destructor
///
AddNetChargeConstraintMover::~AddNetChargeConstraintMover()= default;

/// @brief Copy this object and return a pointer to the copy.
///
protocols::moves::MoverOP AddNetChargeConstraintMover::clone() const { return protocols::moves::MoverOP( new protocols::aa_composition::AddNetChargeConstraintMover( *this ) ); }

/// @brief Create a new object of this type and return a pointer to it.
///
protocols::moves::MoverOP AddNetChargeConstraintMover::fresh_instance() const { return protocols::moves::MoverOP( new AddNetChargeConstraintMover ); }

/// @brief Actually apply the mover to a pose.
///
void
AddNetChargeConstraintMover::apply( Pose &pose )
{
	runtime_assert_string_msg( constraint(), "Error in protocols::aa_composition::AddNetChargeConstraintMover::apply(): The NetChargeConstraint object was not initialized before the apply() function was called." );
	pose.add_constraint( utility::pointer::dynamic_pointer_cast < core::scoring::constraints::Constraint const >( constraint() ) );
	return;
}

/// @brief Parse RosettaScripts XML tag to set up the mover.
///
void
AddNetChargeConstraintMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &datamap,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const &
)
{
	runtime_assert_string_msg( tag->hasOption("filename"), "Error in protocols::aa_composition::AddNetChargeConstraintMover::parse_my_tag(): A \"filename\" option MUST be provided." );
	std::string const compfile( tag->getOption<std::string>("filename") );
	if ( TR.visible() ) TR << "Set filename to " << compfile << "." << std::endl;
	create_constraint_from_file( compfile );

	if ( tag->hasOption("selector") ) {
		std::string const selector_name ( tag->getOption< std::string >( "selector" ) );
		if ( TR.visible() ) TR << "Set selector name to " << selector_name << "." << std::endl;
		core::select::residue_selector::ResidueSelectorCOP residue_selector;
		try {
			residue_selector = datamap.get_ptr< core::select::residue_selector::ResidueSelector const >( "ResidueSelector", selector_name );
		} catch ( utility::excn::Exception & e ) {
			std::string error_message = "Failed to find ResidueSelector named '" + selector_name + "' from the Datamap from AddNetChargeConstraintMover::parse_tag()\n" + e.msg();
			throw CREATE_EXCEPTION(utility::excn::Exception,  error_message );
		}
		runtime_assert( residue_selector );
		add_residue_selector( residue_selector );
	}

	if ( TR.visible() ) TR.flush();
	return;
}

/// @brief Create the NetChargeConstraint object and initialize it from a .comp file.
///
void
AddNetChargeConstraintMover::create_constraint_from_file( std::string const &filename ) {
	runtime_assert_string_msg( !constraint_, "Error in protocols::aa_composition::AddNetChargeConstraintMover::create_constraint_from_file():  The constraint object already has been created!" );
	constraint_ = core::scoring::netcharge_energy::NetChargeConstraintOP( new core::scoring::netcharge_energy::NetChargeConstraint() );
	constraint_->initialize_from_file( filename );
	if ( TR.visible() ) {
		TR << "Initialized NetChargeConstraint object from file " << filename << "." << std::endl;
		TR.flush();
	}
	return;
}

/// @brief Create the NetChargeConstraint object from the data from a .comp file.
/// @details Allows external code to create the constraint object without having it read directly from disk.
void
AddNetChargeConstraintMover::create_constraint_from_file_contents(
	std::string const &filecontents
) {
	runtime_assert_string_msg( !constraint_, "Error in protocols::aa_composition::AddNetChargeConstraintMover::create_constraint_from_filecontents():  The constraint object already has been created!" );
	constraint_ = core::scoring::netcharge_energy::NetChargeConstraintOP( new core::scoring::netcharge_energy::NetChargeConstraint() );
	constraint_->initialize_from_file_contents( filecontents );
	if ( TR.visible() ) {
		TR << "Initialized NetChargeConstraint object from file contents:\n" << filecontents << std::endl;
		TR.flush();
	}
	return;
}

/// @brief Add a ResidueSelector to the constraint to use as a mask.
/// @details The constraint must already have been created with the create_constraint_from_file() function before this function is called.
void
AddNetChargeConstraintMover::add_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector_in ) {
	runtime_assert_string_msg( constraint_, "Error in protocols::aa_composition::AddNetChargeConstraintMover::add_residue_selector(): The NetChargeConstraint object was not initialized before the add_residue_selector() function was called." );
	constraint_->set_selector( selector_in );
	return;
}

std::string AddNetChargeConstraintMover::get_name() const {
	return mover_name();
}

std::string AddNetChargeConstraintMover::mover_name() {
	return "AddNetChargeConstraintMover";
}

void AddNetChargeConstraintMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute("filename", xs_string, "Name of net charge constraint file.  (These files have a \".charge\" suffix)." );
	// I would consider a complex type that ends in .pdb, but I want users to be
	// able to provide cifs! Or name their files whatever!
	attlist + XMLSchemaAttribute("selector", xs_string, "Residue selector named somewhere else in the script.  If provided, the net charge of the selected region is constrained; otherwise the global net charge is constrained." );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Add net charge constraints from the provided file to the selected region or whole pose.", attlist );
}

std::string AddNetChargeConstraintMoverCreator::keyname() const {
	return AddNetChargeConstraintMover::mover_name();
}

protocols::moves::MoverOP
AddNetChargeConstraintMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new AddNetChargeConstraintMover );
}

void AddNetChargeConstraintMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AddNetChargeConstraintMover::provide_xml_schema( xsd );
}


} // aa_composition
} // protocols
