// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/chemistries/PatchChemistry.cc
/// @brief  A chemistry which applies a patch
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <protocols/chemistries/PatchChemistry.hh>
#include <protocols/chemistries/PatchChemistryCreator.hh>

#include <protocols/chemistries/util.hh>

#include <core/pose/Pose.hh>
#include <core/chemical/AtomRefMapping.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/PatchOperation.hh>
#include <core/chemical/PatchOperationFactory.hh>

#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/Tag.hh>

#include <basic/Tracer.hh>

namespace protocols {
namespace chemistries {

static basic::Tracer TR("protocols.chemistry.PatchChemistry");

ChemistryOP
PatchChemistryCreator::create_chemistry() const {
	return utility::pointer::make_shared< PatchChemistry >();
}

std::string
PatchChemistryCreator::keyname() const {
	return PatchChemistry::class_name();
}

void
PatchChemistryCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	PatchChemistry::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////

void
PatchChemistry::apply( core::chemical::MutableResidueType & restype ) {

	core::chemical::MutableResidueTypeOP new_type;

	if ( ! patch_name_.empty() ) {
		TR.Warning << "Applying a named patch without a Pose context -- assuming the global full atom ResidueTypeSet." << std::endl;
		core::chemical::ResidueTypeSetCOP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FULL_ATOM_t );
		debug_assert( rts != nullptr );
		new_type = apply_patches( restype, rts->get_patches( patch_name_ ) );
	} else if ( ! patch_file_.empty() ) {
		core::chemical::Patch patch{ restype.mode() };
		patch.read_file( patch_file_ );

		new_type = apply_patch( restype, patch );
	} else if ( ! operations_.empty() ) {
		new_type = restype.clone();
	} else {
		TR.Warning << "PatchChemistry hasn't been given any patch information!" << std::endl;
	}

	if ( new_type ) {
		apply_operations( *new_type );
		update_type( restype, *new_type );
	} else {
		TR.Warning << "PatchChemistry was not able to successfully apply patches." << std::endl;
	}
}

void
PatchChemistry::apply( core::chemical::MutableResidueType & restype, core::pose::Pose const & pose ) {

	if ( patch_name_.empty() ) { // We only special case for name in RTS
		apply( restype );
		return;
	}

	core::chemical::ResidueTypeSetCOP rts = pose.residue_type_set_for_pose( restype.mode() );
	debug_assert( rts != nullptr );
	core::chemical::MutableResidueTypeOP new_type = apply_patches( restype, rts->get_patches( patch_name_ ) );

	if ( new_type ) {
		apply_operations( *new_type );
		update_type( restype, *new_type );
	} else {
		TR.Warning << "PatchChemistry was not able to successfully apply patches." << std::endl;
	}
}

core::chemical::VDVDMapping
PatchChemistry::get_mapping() const {

	return mapping_;
}

void
PatchChemistry::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;

	utility::tag::AttributeList attributes;

	attributes
		+ XMLSchemaAttribute( "patch_name", xs_string, "The name of the patch from the ResidueTypeSet to apply.")
		+ XMLSchemaAttribute( "patch_file", xs_string, "The name of a file from which to get the patches.");

	XMLSchemaSimpleSubelementList subelements;
	AttributeList subattlist;
	subattlist + XMLSchemaAttribute::required_attribute("line", xs_string, "The patch file line contents.");

	subelements.add_simple_subelement( "Op", subattlist, "A manually specified patch operation line" );

	xsd_type_definition_w_attributes_and_repeatable_subelements(
		xsd, class_name(),
		"A patch chemistry applies a given patch to the residue type.",
		attributes, subelements );
}

void
PatchChemistry::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &)
{

	patch_name( tag->getOption< std::string >("patch_name", "") );
	patch_file( tag->getOption< std::string >("patch_file", "") );

	for ( utility::tag::TagCOP const & subtag: tag->getTags() ) {
		add_patch_operation_line( subtag->getOption< std::string >("line") );
	}

	TR << "Creating PatchChemistry with name: '" << patch_name_ << "' file: '" << patch_file_ << "' and " << operations_.size() << " operations." << std::endl;
	for ( auto const & op: operations_ ) {
		TR.Debug << "\t" << op->name() << std::endl;
	}
}

std::string
PatchChemistry::class_name() {
	return "PatchChemistry";
}

void
PatchChemistry::add_patch_operation_line( std::string const & line ) {
	std::map< std::string, core::Real > charge_reassign; // Not needed.

	core::chemical::PatchOperationFactory const & pofact = *core::chemical::PatchOperationFactory::get_instance();
	std::istringstream l( line );
	std::string tag;
	l >> tag; // Need to pull the tag to advance the stream to the next location.

	core::chemical::PatchOperationOP op = pofact.newPatchOperation( tag, line, l, charge_reassign );
	if ( op ) {
		operations_.push_back( op );
	}
}

core::chemical::MutableResidueTypeOP
PatchChemistry::apply_patches( core::chemical::MutableResidueType const & restype, utility::vector1< core::chemical::PatchCOP > const & patches ) const {
	for ( core::chemical::PatchCOP const & p: patches ) {
		core::chemical::MutableResidueTypeOP patched = apply_patch( restype, *p );
		if ( patched ) return patched;
	}
	return nullptr;
}

core::chemical::MutableResidueTypeOP
PatchChemistry::apply_patch( core::chemical::MutableResidueType const & restype, core::chemical::Patch const & patch ) const {
	if ( ! patch.applies_to( restype ) ) {
		TR.Warning << "Patch " << patch.name() << " does not apply to restype " << restype.name() << std::endl;
		return nullptr;
	}

	return patch.apply( restype );
}

void
PatchChemistry::apply_operations( core::chemical::MutableResidueType & restype ) const {
	for ( core::chemical::PatchOperationOP const & operation: operations_ ) {
		if ( operation->apply( restype ) ) {
			TR.Warning << "Operation " << operation->name() << " failed on ResidueType " << restype.name() << std::endl;
		} // else we successfully applied
	}
}

void
PatchChemistry::update_type( core::chemical::MutableResidueType & restype, core::chemical::MutableResidueType const & new_restype ) {

	// A name-based mapping isn't the best, but it's the best we've got at the moment.
	mapping_ = core::chemical::VDNameMapping(restype).downstream_combine( core::chemical::NameVDMapping(new_restype) );

	// Now we can overwrite
	restype = new_restype;
}


}
}


