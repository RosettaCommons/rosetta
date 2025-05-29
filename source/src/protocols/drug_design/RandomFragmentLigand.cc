// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/drug_design/RandomFragmentLigand.hh
/// @brief Fragment a ligand to make it smaller
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <protocols/drug_design/RandomFragmentLigand.hh>
#include <protocols/drug_design/RandomFragmentLigandCreator.hh>

#include <core/chemical/ResidueGraphTypes.hh>
#include <core/chemical/MutableResidueType.hh>
#include <core/chemical/Atom.hh>
#include <core/chemical/Bond.hh>
#include <core/chemical/Elements.hh>
#include <core/chemical/Element.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/bond_support.hh>
#include <core/chemical/ChemicalManager.hh>
#include <protocols/chemistries/util.hh>

#include <core/chemical/gasteiger/GasteigerAtomTypeSet.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeData.hh>
#include <core/chemical/gasteiger/GasteigerAtomTyper.hh>
#include <core/chemical/modifications/ValenceHandler.hh>
#include <core/chemical/residue_support.hh>
#include <core/chemical/atomtype_support.hh>

#include <numeric/random/random.hh>

#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/Tracer.hh>

#include <boost/graph/connected_components.hpp>

namespace protocols {
namespace drug_design {

static basic::Tracer TR("protocols.drug_design.RandomFragmentLigand");

//------------------------- Creator -----------------------------

protocols::chemistries::ChemistryOP
RandomFragmentLigandCreator::create_chemistry() const {
	return protocols::chemistries::ChemistryOP( new RandomFragmentLigand );
}

std::string
RandomFragmentLigandCreator::keyname() const {
	return RandomFragmentLigand::class_name();
}

void
RandomFragmentLigandCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	RandomFragmentLigand::provide_xml_schema( xsd );
}



//------------------------- Chemistry -----------------------------

RandomFragmentLigand::RandomFragmentLigand():
	Chemistry(class_name()),
	keep_bigger_( false ),
	ccbond_( false )
{}

void
RandomFragmentLigand::keep_atom( std::string const & keep_atom ) {
	keep_atom_ = keep_atom;
	utility::strip_whitespace( keep_atom_ );
}

void
RandomFragmentLigand::apply( core::chemical::MutableResidueType & restype )
{
	using namespace core::chemical;

	if ( keep_bigger_ && keep_atom_ != "" ) {
		utility_exit_with_message("Cannot set both keep_bigger and keep_atom for RandomFragmentLigand.");
	}

	ResidueGraph const & graph( restype.graph() );

	// Make sure rings are annotated.
	find_bonds_in_rings( restype );

	// Pick a random single bond, not on a ring and not connected to a hydrogen
	utility::vector1<ED> valid_bonds;
	EIter iter, iter_end;
	for ( boost::tie(iter, iter_end) = boost::edges(graph); iter != iter_end; ++iter ) {
		Atom const & a1( graph[ boost::source( *iter, graph ) ] );
		assert( a1.element_type() );
		Atom const & a2( graph[ boost::target( *iter, graph ) ] );
		assert( a2.element_type() );
		if ( graph[ *iter ].ringness() != BondNotInRing ) continue;
		if ( a1.element_type()->element() == element::H ||
				a2.element_type()->element() == element::H ) continue;
		// If ccbond_ is set, skip bonds which aren't carbon-carbon
		if ( ccbond_ &&
				( a1.element_type()->element() != element::C ||
				a2.element_type()->element() != element::C ) ) continue;

		valid_bonds.push_back( *iter );
	}

	if ( valid_bonds.size() == 0 ) {
		TR << "No valid bonds on which to fragment ligand - skipping." << std::endl;
		set_last_status( core::chemical::modifications::FAIL_DO_NOT_RETRY );
		return;
	}

	ED selected_bond = numeric::random::rg().random_element(valid_bonds);
	VD a1( boost::source( selected_bond, graph ) );
	VD a2( boost::target( selected_bond, graph ) );

	TR << "Deleting bond between " << restype.atom_name( a1 ) << " (" << a1 << ") and " << restype.atom_name( a2 ) << " (" << a2 << ")" << std::endl;
	// Delete the bond.
	restype.delete_bond( a1, a2 );

	// Determine the connected components (Do either side of bond, look at connections recursively).
	utility::vector1<VD> side1, side2;

	typedef std::map< core::chemical::VD, core::Size > ComponentMap;
	ComponentMap components;
	core::Size ncomponent;
	boost::associative_property_map< ComponentMap > component_property_map(components);
	ncomponent = boost::connected_components( graph, component_property_map );

	if ( ncomponent != 2 ) {
		utility_exit_with_message("Expected to generate two pieces when fragmenting along non-ring single bond. Found "+utility::to_string(ncomponent)+" instead." );
	}
	VD side1_attach(MutableResidueType::null_vertex), side2_attach(MutableResidueType::null_vertex);
	for ( ComponentMap::const_iterator iter(components.begin()), iter_end(components.end()); iter != iter_end; ++iter ) {
		// The connected_components documentation doesn't mention, but I'm assuming the components are numbered 0 to ncomponent-1
		if ( iter->second == 0 ) {
			side1.push_back( iter->first );
			if ( iter->first == a1 || iter->first == a2 ) {
				side1_attach = iter->first;
			}
		} else if (  iter->second == 1 ) {
			side2.push_back( iter->first );
			if ( iter->first == a1 || iter->first == a2 ) {
				side2_attach = iter->first;
			}
		} else {
			utility_exit_with_message("Found unexpected component "+utility::to_string(iter->second) );
		}
	}

	assert( side1_attach != MutableResidueType::null_vertex && side2_attach != MutableResidueType::null_vertex);

	// Pick which atoms to keep - either the bigger side, the side with a certain atom, or random
	VD keep_atom = side1_attach; // The atom on the bond we're breaking that's on the side we're keeping.

	if ( keep_atom_.size() ) {
		// Search for atom
		bool found1(false), found2(false);
		for ( core::Size ii(1); ii <= side1.size(); ++ii ) {
			if ( utility::stripped_whitespace( graph[ side1[ii] ].name() ) == keep_atom_ ) {
				found1 = true;
				break;
			}
		}
		for ( core::Size ii(1); ii <= side2.size(); ++ii ) {
			if ( utility::stripped_whitespace( graph[ side2[ii] ].name() ) == keep_atom_ ) {
				found2 = true;
				break;
			}
		}
		if ( !found1 && !found2 ) {
			utility_exit_with_message("Cannot find atom named "+keep_atom_+" in residue type "+restype.name());
		}
		if ( !found1 && found2 ) {
			side1.swap(side2); // Put the side to keep first
			keep_atom=side2_attach;
		}
	} else if ( keep_bigger_ ) {
		if ( side2.size() > side1.size() ) {
			side1.swap(side2); // Put the side to keep first
			keep_atom=side2_attach;
		}
	} else { // Randomly pick the side to keep.
		if ( numeric::random::random_range(0,1) == 1 ) {
			side1.swap(side2); // Put the side to keep first
			keep_atom=side2_attach;
		}
	}

	// Delete the atoms which aren't on the side we want

	for ( core::Size ii(1); ii <= side2.size(); ++ii ) { // VD's should stay valid through the deletion
		restype.delete_atom( side2[ii] );
	}
	if ( ! restype.has( keep_atom ) ) {
		TR << "Can't find atom " << keep_atom << std::endl;
		restype.show_all_atom_names(TR);
		runtime_assert( restype.has( keep_atom ) );
	}

	// All atoms now in ResidueTypes should keep the same VDs in the final residue
	mapping_ = core::chemical::VDVDMapping(restype); // Identity map for existing atoms.

	// Add a hydrogen to the atom we just kept.

	// TODO: From the AddHydrogens mover - we probably want to generalize this out.
	if ( restype.gasteiger_atom_typeset() == nullptr ) {
		restype.set_gasteiger_atom_typeset(core::chemical::ChemicalManager::get_instance()->gasteiger_atom_type_set("default"));
	}
	//assign the gastieger atom types, which is an absolute if you want to add hydrogens
	core::chemical::gasteiger::assign_gasteiger_atom_types( restype, restype.gasteiger_atom_typeset(), /*keep_existing=*/ false );

	utility::vector1<numeric::xyzVector<core::Real> > coords( modifications::determine_coordinates(restype, keep_atom));
	for ( core::Size new_hydrogens=1; new_hydrogens <= coords.size(); ++new_hydrogens ) {
		core::Size hydro_number( restype.natoms() + 1 );
		std::string name = "H" + utility::to_string( hydro_number );
		while ( restype.has( name ) ) {
			++hydro_number;
			name = "H" + utility::to_string( hydro_number );
		}
		VD vertex_hydrogen = restype.add_atom(name);
		restype.atom(vertex_hydrogen).ideal_xyz(coords[new_hydrogens]);
		restype.atom(vertex_hydrogen).element_type(restype.element_set().element(element::name_from_elements(element::H))); //dont forget to add the element type!
		restype.add_bond(vertex_hydrogen, keep_atom, SingleBond);
	}

	// Redo the entire atom tree, typing, etc.

	//restype.setup_atom_ordering(); // Needed to redo atom_ordering map, etc. TODO: We need to fix this up to be consistent.

	// TODO: also from AddHydrogens - need to normalize, etc.
	//rename_atoms(res, false); //not sure why, but if you try and add more hydrogens, shit hits the fan. Try and re do this
	rosetta_retype_fullatom(restype, false); //need to do this, fails currently

	rosetta_recharge_fullatom(restype);
	core::chemical::find_bonds_in_rings( restype );

	VD nbr_atom = restype.nbr_vertex();
	if ( nbr_atom == MutableResidueType::null_vertex || ! restype.has( nbr_atom ) || restype.atom(nbr_atom).element_type()->element() == element::H ) {
		TR << "Dumping neighbor atom and recreating." << std::endl;
		nbr_atom = MutableResidueType::null_vertex;
	}
	restype.nbr_radius( find_nbr_dist( restype, nbr_atom ) );
	restype.nbr_atom( nbr_atom );
	restype.assign_internal_coordinates( nbr_atom ); // Also sets atom base. Needs nbr atom assignment
	restype.autodetermine_chi_bonds(); // Must be after internal coordinate setup
	set_last_status( core::chemical::modifications::SUCCESS );
}

core::chemical::VDVDMapping
RandomFragmentLigand::get_mapping() const {
	return mapping_;
}

std::string
RandomFragmentLigand::class_name() {
	return "RandomFragmentLigand";
}

void
RandomFragmentLigand::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default("keep_bigger", xsct_rosetta_bool,
		"If true, keep the fragment which is bigger instead of a random one", "0")
		+ XMLSchemaAttribute("keep_atom", xs_string,
		"If set, keep the fragment which contains the atom of the given name")
		+ XMLSchemaAttribute::attribute_w_default("ccbonds", xsct_rosetta_bool,
		"If true, only fragment on carbon-carbon bonds", "0");

	protocols::chemistries::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"Fragment a ResidueType on a random bond, and discard one half. ",
		attlist );
}

/// @brief Initialize any data members of this instance from an input tag
/// and a DataMap object
void
RandomFragmentLigand::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & ) {

	keep_bigger( tag->getOption<bool>("keep_bigger", false) );
	keep_atom( tag->getOption<std::string>("keep_atom", "") );

	if ( keep_bigger_ && keep_atom_ != "" ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Cannot set both keep_bigger and keep_atom for RandomFragmentLigand.");
	}

	ccbond( tag->getOption<bool>("ccbonds", false) );

}


}
}
