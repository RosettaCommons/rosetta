// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/aa_composition/AddMHCEpitopeConstraintMover.cc
/// @brief Assigns an MHCEpitopeConstraint to a pose, initializing it from a file.
/// @author Chris Bailey-Kellogg (cbk@cs.dartmouth.edu), based on Vikram Mulligan's MHCEpitopeConstraint

#include <protocols/aa_composition/AddMHCEpitopeConstraintMover.hh>
#include <protocols/aa_composition/AddMHCEpitopeConstraintMoverCreator.hh>


#include <core/pose/Pose.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/scoring/mhc_epitope_energy/MHCEpitopeConstraint.hh>
#include <core/scoring/mhc_epitope_energy/MHCEpitopeEnergySetup.hh>

//Auto Headers

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>
#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.aa_composition.AddMHCEpitopeConstraintMover" );

namespace protocols {
namespace aa_composition {

using namespace core;
using namespace core::scoring;
using namespace constraints;
using namespace utility::tag;

/// @brief Default Constructor
///
AddMHCEpitopeConstraintMover::AddMHCEpitopeConstraintMover():
	protocols::moves::Mover( AddMHCEpitopeConstraintMover::mover_name() ),
	//TODO initialize variables here
	constraint_()
{
}

/// @brief Copy Constructor
///
AddMHCEpitopeConstraintMover::AddMHCEpitopeConstraintMover( AddMHCEpitopeConstraintMover const &src ):
	protocols::moves::Mover( AddMHCEpitopeConstraintMover::mover_name() ),
	//TODO initialize variables here
	constraint_( src.constraint_ )
{
}

/// @brief Destructor
///
AddMHCEpitopeConstraintMover::~AddMHCEpitopeConstraintMover()= default;

/// @brief Copy this object and return a pointer to the copy.
///
protocols::moves::MoverOP AddMHCEpitopeConstraintMover::clone() const { return utility::pointer::make_shared< protocols::aa_composition::AddMHCEpitopeConstraintMover >( *this ); }

/// @brief Create a new object of this type and return a pointer to it.
///
protocols::moves::MoverOP AddMHCEpitopeConstraintMover::fresh_instance() const { return utility::pointer::make_shared< AddMHCEpitopeConstraintMover >(); }

/// @brief Actually apply the mover to a pose.
///
void
AddMHCEpitopeConstraintMover::apply( Pose &pose )
{
	runtime_assert_string_msg( constraint(), "Error in protocols::aa_composition::AddMHCEpitopeConstraintMover::apply(): The MHCEpitopeConstraint object was not initialized before the apply() function was called." );
	TR << "Applying MHCEpitopeConstraint with a weight of " << constraint_->get_cst_weight() << " and the " << constraint_->get_selector_name() << " residue selector." << std::endl;
	TR << "Constraints set with the following settings: " << constraint_->mhc_epitope_energy_setup()->report();
	TR.flush();
	pose.add_constraint( utility::pointer::dynamic_pointer_cast < core::scoring::constraints::Constraint const >( constraint() ) );
	return;
}

/// @brief Parse RosettaScripts XML tag to set up the mover.
///
void
AddMHCEpitopeConstraintMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &datamap,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const &
)
{
	runtime_assert_string_msg( tag->hasOption("filename"), "Error in protocols::aa_composition::AddMHCEpitopeConstraintMover::parse_my_tag(): A \"filename\" option MUST be provided." );
	std::string const mhcfile( tag->getOption<std::string>("filename") );
	if ( TR.visible() ) TR << "Set filename to " << mhcfile << "." << std::endl;
	create_constraint_from_file( mhcfile );

	if ( tag->hasOption("selector") ) {
		std::string const selector_name ( tag->getOption< std::string >( "selector" ) );
		if ( TR.visible() ) TR << "Set selector name to " << selector_name << "." << std::endl;
		core::select::residue_selector::ResidueSelectorCOP residue_selector;
		try {
			residue_selector = datamap.get_ptr< core::select::residue_selector::ResidueSelector const >( "ResidueSelector", selector_name );
		} catch ( utility::excn::Exception & e ) {
			std::string error_message = "Failed to find ResidueSelector named '" + selector_name + "' from the Datamap from AddMHCEpitopeConstraintMover::parse_tag()\n" + e.msg();
			throw CREATE_EXCEPTION(utility::excn::Exception,  error_message );
		}
		runtime_assert( residue_selector );
		add_residue_selector( residue_selector );
		constraint_->set_selector_name( selector_name );
	}

	//Get the cst-specific weight from tag.
	if ( tag->hasOption("weight") ) {
		add_weight( tag->getOption<core::Real>("weight") );
		if ( TR.visible() ) TR << "Set MHC constraint weight to " << tag->getOption<core::Real>("weight") << "." << std::endl;
	}

	if ( TR.visible() ) TR.flush();
	return;
}

/// @brief Create the MHCEpitopeConstraint object and initialize it from a .mhc file.
///
void
AddMHCEpitopeConstraintMover::create_constraint_from_file( std::string const &filename ) {
	runtime_assert_string_msg( !constraint_, "Error in protocols::aa_composition::AddMHCEpitopeConstraintMover::create_constraint_from_file():  The constraint object already has been created!" );
	constraint_ = utility::pointer::make_shared< core::scoring::mhc_epitope_energy::MHCEpitopeConstraint >();
	constraint_->initialize_from_file( filename );
	if ( TR.visible() ) {
		TR << "Initialized MHCEpitopeConstraint object from file " << filename << "." << std::endl;
		TR.flush();
	}
	return;
}

/// @brief Create the MHCEpitopeConstraint object from the data from a .mhc file.
/// @details Allows external code to create the constraint object without having it read directly from disk.
void
AddMHCEpitopeConstraintMover::create_constraint_from_file_contents(
	std::string const &filecontents
) {
	runtime_assert_string_msg( !constraint_, "Error in protocols::aa_composition::AddMHCEpitopeConstraintMover::create_constraint_from_filecontents():  The constraint object already has been created!" );
	constraint_ = utility::pointer::make_shared< core::scoring::mhc_epitope_energy::MHCEpitopeConstraint >();
	constraint_->initialize_from_file_contents( filecontents );
	if ( TR.visible() ) {
		TR << "Initialized MHCEpitopeConstraint object from file contents:\n" << filecontents << std::endl;
		TR.flush();
	}
	return;
}

/// @brief Add a ResidueSelector to the constraint to use as a mask.
/// @details The constraint must already have been created with the create_constraint_from_file() function before this function is called.
void
AddMHCEpitopeConstraintMover::add_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector_in ) {
	runtime_assert_string_msg( constraint_, "Error in protocols::aa_composition::AddMHCEpitopeConstraintMover::add_residue_selector(): The MHCEpitopeConstraint object was not initialized before the add_residue_selector() function was called." );
	constraint_->set_selector( selector_in );
	return;
}

/// @brief Add an individual weight to the constraint.
/// @details The constraint must already have been created with the create_constraint_from_file() function before this function is called.
void
AddMHCEpitopeConstraintMover::add_weight( core::Real cst_weight ) {
	runtime_assert_string_msg( constraint_, "Error in protocols::aa_composition::AddMHCEpitopeConstraintMover::add_residue_selector(): The MHCEpitopeConstraint object was not initialized before the add_weight() function was called." );
	constraint_->set_cst_weight( cst_weight );
	return;
}

std::string AddMHCEpitopeConstraintMover::get_name() const {
	return mover_name();
}

std::string AddMHCEpitopeConstraintMover::mover_name() {
	return "AddMHCEpitopeConstraintMover";
}

void AddMHCEpitopeConstraintMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute("filename", xs_string, "Name of mhc epitope constraint file.  (These files have a \".mhc\" suffix)." );
	attlist + XMLSchemaAttribute("selector", xs_string, "Residue selector named somewhere else in the script.  If provided, the mhc epitope score of the selected region is constrained; otherwise the global score is constrained." );
	attlist + XMLSchemaAttribute::attribute_w_default("weight", xsct_real, "Adjust the mhc_epitope weight for the epitope predictor set up by this constraint mover ONLY.  The weight will be the mhc_epitope weight of the scorefunction * this weight.", "1.0" );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Add mhc epitope constraints from the provided file to the selected region or whole pose.", attlist );
}

std::string AddMHCEpitopeConstraintMoverCreator::keyname() const {
	return AddMHCEpitopeConstraintMover::mover_name();
}

protocols::moves::MoverOP
AddMHCEpitopeConstraintMoverCreator::create_mover() const {
	return utility::pointer::make_shared< AddMHCEpitopeConstraintMover >();
}

void AddMHCEpitopeConstraintMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AddMHCEpitopeConstraintMover::provide_xml_schema( xsd );
}


} // aa_composition
} // protocols
