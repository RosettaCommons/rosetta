// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/metal_interface/RemoveMetalConnectionsMover.cc
/// @brief A mover that removes the connections to metals that were added by the SetupMetalsMover
/// or by the -auto_setup_metals flag.
/// @details This mover:
///     - Removes the bonds between metals and metal-binding residues.
///     - Reverts metal-liganding residues back to their pre-bonded types.
///     - Reverts metals back to their pre-bonded types.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Unit headers
#include <protocols/metal_interface/RemoveMetalConnectionsMover.hh>
#include <protocols/metal_interface/RemoveMetalConnectionsMoverCreator.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/ResidueTypeSet.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>
#include <utility/string_util.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// Citation Manager
#include <utility/vector1.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

static basic::Tracer TR( "protocols.metal_interface.RemoveMetalConnectionsMover" );

namespace protocols {
namespace metal_interface {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
RemoveMetalConnectionsMover::RemoveMetalConnectionsMover():
	protocols::moves::Mover( RemoveMetalConnectionsMover::mover_name() )
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
RemoveMetalConnectionsMover::~RemoveMetalConnectionsMover() = default;

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
RemoveMetalConnectionsMover::apply(
	core::pose::Pose & pose
) {
	// Generate the subset of residues to which we are applying this selector:
	core::select::residue_selector::ResidueSubset const selection(
		residue_selector_ == nullptr ?
		core::select::residue_selector::ResidueSubset(pose.total_residue(), true) :
		residue_selector_->apply(pose)
	);

	// Iterate over all residues, and apply to any pairs of metals that we encounter.
	// By doing this backwards, it is more likely that the intermediate residue types that
	// we create have already been created:
	for ( core::Size i(pose.total_residue()); i>=1; --i ) {
		if ( !selection[i] ) continue; //Skip residues that aren't selected.
		if ( pose.residue_type(i).is_metal() ) continue; //Skip metals.  They're processed when we process the metal-liganding residue.
		if ( !pose.residue_type(i).is_metalbinding() ) continue; //Skip residue types that can't bind metal.

		//Iterate over this residue's connections:
		for ( core::Size iconn(1); iconn<=pose.residue_type(i).n_possible_residue_connections(); ++iconn ) {
			//Skip unconnected:
			if ( pose.residue(i).connection_incomplete(iconn) ) continue;

			//Residue connected to this residue at connection iconn:
			core::Size const connected_res_index( pose.residue(i).connected_residue_at_resconn(iconn) );
			if ( !selection[connected_res_index] ) continue; //Skip if we haven't selected the connected residue.

			//Skip if this isn't a connection to a metal:
			core::chemical::ResidueType const & connected_res_type( pose.residue_type(connected_res_index) );
			if ( !connected_res_type.is_metal() ) continue;

			//If we reach here, we have a metal ligand-metal bond, and we should remove it.
			//First, store the atom names of the connection:
			core::Size const connected_res_iconn( pose.residue(i).residue_connection_conn_id(iconn) ); //The connection index on the OTHER residue that connets to THIS residue.
			std::string const atomname_thisres( utility::strip( pose.residue_type(i).atom_name( pose.residue_type(i).residue_connect_atom_index(iconn) ) ) ); //Atom name on THIS residue.
			std::string const atomname_otherres( utility::strip( connected_res_type.atom_name( connected_res_type.residue_connect_atom_index( connected_res_iconn ) ) ) ); //Atom name on OTHER residue.

			//Second, delete the bond (using sever_chemical_bond):
			TR << "Severing chemical bond between residue " << pose.residue_type(i).base_name() << i << ", atom " << atomname_thisres << " and metal " << connected_res_type.base_name() << connected_res_index << ", atom " << atomname_otherres << "." << std::endl;
			pose.conformation().sever_chemical_bond( i, iconn, connected_res_index, connected_res_iconn );

			//Third, try to remove the variant type on this residue:
			remove_variant_types_from_res( pose, i, atomname_thisres, false );

			//Fourth, try to remove the variant type on the other residue:
			remove_variant_types_from_res( pose, connected_res_index, atomname_otherres, true );

			//Decrement iconn since we've removed a connection; next pass increments so that we're trying a
			//new connection ID:
			--iconn;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
RemoveMetalConnectionsMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
RemoveMetalConnectionsMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
) {
	//Parse the residue selector:
	if ( tag->hasOption("residue_selector") ) {
		set_residue_selector( core::select::residue_selector::parse_residue_selector( tag, datamap, "residue_selector" ) );
	}
}
void RemoveMetalConnectionsMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector", "A residue selector for selecting those residues from which to remove metal bonds and variant types.  If a residue AND the metal that it binds are both selected, bonds and variant types are removed from both.  The mover applies to the whole pose if no residue selector is provided." );

	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"A mover that removes bonds between metals and metal-liganding atoms that were added with the SetupMetalsMover or with the -auto_setup_metals commandline flag.  This mover also reverts the metal residue and the metal-liganding residue back to the types that lack the extra bonds.  Note that at the present time it does not remove metal constraints.",
		attlist
	);
}

/// @brief Sets the residue selector for selecting residues from which metal bonds
/// will be removed.
/// @details If a polymer residue AND a metal that it binds are selected, the corresponding
/// metal bonds and variants are removed.
void
RemoveMetalConnectionsMover::set_residue_selector(
	core::select::residue_selector::ResidueSelectorCOP const & selector_in
) {
	residue_selector_ = selector_in;
}

/// @brief Gets the residue selector for selecting residues from which metal bonds
/// will be removed.
/// @details If a polymer residue AND a metal that it binds are selected, the corresponding
/// metal bonds and variants are removed.
core::select::residue_selector::ResidueSelectorCOP
RemoveMetalConnectionsMover::residue_selector() const {
	return residue_selector_;
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
RemoveMetalConnectionsMover::fresh_instance() const
{
	return utility::pointer::make_shared< RemoveMetalConnectionsMover >();
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
RemoveMetalConnectionsMover::clone() const
{
	return utility::pointer::make_shared< RemoveMetalConnectionsMover >( *this );
}

std::string RemoveMetalConnectionsMover::get_name() const {
	return mover_name();
}

std::string RemoveMetalConnectionsMover::mover_name() {
	return "RemoveMetalConnectionsMover";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
RemoveMetalConnectionsMoverCreator::create_mover() const
{
	return utility::pointer::make_shared< RemoveMetalConnectionsMover >();
}

std::string
RemoveMetalConnectionsMoverCreator::keyname() const
{
	return RemoveMetalConnectionsMover::mover_name();
}

void RemoveMetalConnectionsMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RemoveMetalConnectionsMover::provide_xml_schema( xsd );
}

/// @brief Provide citations to the passed CitationCollectionList.
/// @details This mover is unpublished.  It returns Vikram K. Mulligan as its author.
void
RemoveMetalConnectionsMover::provide_citation_info(
	basic::citation_manager::CitationCollectionList & citations
) const {
	using namespace basic::citation_manager;
	citations.add(
		utility::pointer::make_shared< UnpublishedModuleInfo >(
		"RemoveMetalConnectionsMover",
		CitedModuleType::Mover,
		"Vikram K. Mulligan",
		"Center for Computational Biology, Flatiron Institute",
		"vmulligan@flatironinstitute.org",
		"Created this mover."
		)
	);
	if ( residue_selector_ != nullptr ) {
		residue_selector_->provide_citation_info( citations );
	}
}

////////////////////////////////////////////////////////////////////////////////
/// Private methods ///
///////////////////////

/// @brief Remove one or two variant types from a pose residue.
/// @details Based on name lookup, so hopefully should work for on-the-fly variants.
void
RemoveMetalConnectionsMover::remove_variant_types_from_res(
	core::pose::Pose & pose,
	core::Size const resindex,
	std::string const & this_atom,
	bool const this_is_metal
) const {
	//Generate the variant names to remove:
	std::string const variant1(
		this_is_metal ?
		"MP-" + this_atom + "-metal_connect" :
		"MP-" + this_atom + "-connect"
	);
	std::string const variant2(
		this_is_metal ?
		"" :
		"MP-" + this_atom + "-pruneH"
	);

	//Summarize what we're doing for the user:
	TR << "\tAttempting to remove \"" << variant1 << "\" ";
	if ( !this_is_metal ) {
		TR << "and \"" << variant2 << "\" types ";
	} else {
		TR << "type ";
	}
	TR << "from pose residue " << pose.residue_type(resindex).base_name() << resindex << "." << std::endl;

	//The current residue type:
	core::chemical::ResidueType const & this_restype( pose.residue_type(resindex) );

	//Find the new residue type and replace:
	std::string new_resname(this_restype.name());
	if ( !this_is_metal ) {
		new_resname = utility::replace_first_in( new_resname, ":"+variant2, "" );
	}
	new_resname = utility::replace_first_in( new_resname, ":"+variant1, "" );
	core::chemical::ResidueTypeSetCOP typeset( pose.residue_type_set_for_pose( this_restype.mode() )  );
	core::chemical::ResidueTypeCOP new_type( typeset->name_mapOP(new_resname) );
	runtime_assert( new_type != nullptr );
	core::pose::replace_pose_residue_copying_existing_coordinates( pose, resindex, *new_type );

}


std::ostream &
operator<<( std::ostream & os, RemoveMetalConnectionsMover const & mover )
{
	mover.show(os);
	return os;
}


} //metal_interface
} //protocols
