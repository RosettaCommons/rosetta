// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

static THREAD_LOCAL basic::Tracer TR( "protocols.aa_composition.AddCompositionConstraintMover" );

namespace protocols {
namespace aa_composition {

using namespace core;
using namespace core::scoring;
using namespace constraints;
using namespace utility::tag;

std::string
AddCompositionConstraintMoverCreator::keyname() const
{
	return AddCompositionConstraintMoverCreator::mover_name();
}

protocols::moves::MoverOP
AddCompositionConstraintMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new AddCompositionConstraintMover );
}

std::string
AddCompositionConstraintMoverCreator::mover_name()
{
	return "AddCompositionConstraintMover";
}

/// @brief Default Constructor
///
AddCompositionConstraintMover::AddCompositionConstraintMover():
	protocols::moves::Mover( AddCompositionConstraintMoverCreator::mover_name() ),
	//TODO initialize variables here
	constraint_()
{
}

/// @brief Copy Constructor
///
AddCompositionConstraintMover::AddCompositionConstraintMover( AddCompositionConstraintMover const &src ):
	protocols::moves::Mover( AddCompositionConstraintMoverCreator::mover_name() ),
	//TODO initialize variables here
	constraint_( src.constraint_ )
{
}


/// @brief Destructor
///
AddCompositionConstraintMover::~AddCompositionConstraintMover(){}

/// @brief Copy this object and return a pointer to the copy.
///
protocols::moves::MoverOP AddCompositionConstraintMover::clone() const { return protocols::moves::MoverOP( new protocols::aa_composition::AddCompositionConstraintMover( *this ) ); }

/// @brief Create a new object of this type and return a pointer to it.
///
protocols::moves::MoverOP AddCompositionConstraintMover::fresh_instance() const { return protocols::moves::MoverOP( new AddCompositionConstraintMover ); }

/// @brief Returns the name of this mover ("AddCompositionConstraintMover").
///
std::string
AddCompositionConstraintMover::get_name() const {
	return AddCompositionConstraintMoverCreator::mover_name();
}

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

/// @brief Add a ResidueSelector to the constraint to use as a mask.
/// @details The constraint must already have been created with the create_constraint_from_file() function before this function is called.
void
AddCompositionConstraintMover::add_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector_in ) {
	runtime_assert_string_msg( constraint_, "Error in protocols::aa_composition::AddCompositionConstraintMover::add_residue_selector(): The AACompositionConstraint object was not initialized before the add_residue_selector() function was called." );
	constraint_->set_selector( selector_in );
	return;
}

} // aa_composition
} // protocols
