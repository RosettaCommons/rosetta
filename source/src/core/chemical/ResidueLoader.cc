// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/ResidueLoader.fwd.hh
/// @brief implementation of the Residue Loader
/// @author Sam DeLuca

#include <core/chemical/ResidueLoader.hh>
#include <core/chemical/ResidueLoaderCreator.hh>
#include <core/chemical/ResidueLoaderOptions.hh>
#include <core/chemical/ResidueType.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/residue_io.hh>


//utility headers
#include <utility/excn/Exceptions.hh>

namespace core {
namespace chemical {

basic::resource_manager::ResourceLoaderOP ResidueLoaderCreator::create_resource_loader() const
{
	return basic::resource_manager::ResourceLoaderOP( new ResidueLoader );
}

std::string ResidueLoaderCreator::loader_type() const
{
	return "ResidueType";
}

utility::pointer::ReferenceCountOP
ResidueLoader::create_resource(
	basic::resource_manager::ResourceOptions const & options,
	basic::resource_manager::LocatorID const & locator_id,
	std::istream & /*istream*/
) const
{
	if ( ! dynamic_cast< ResidueLoaderOptions const * > ( &options ) ) {
		throw utility::excn::EXCN_Msg_Exception( "ResidueLoader expected to be given a ResidueLoaderOptions object, " \
			"but was given a non-ResidueLoaderOptions object of type '" + options.type() + "', which has the name '" + options.name() + "'." );
	}

	ResidueLoaderOptions const & residue_options = static_cast<ResidueLoaderOptions const & >(options);

	AtomTypeSetCOP atom_type_set = ChemicalManager::get_instance()->atom_type_set(residue_options.atom_type_set_tag());
	ElementSetCOP element_set = ChemicalManager::get_instance()->element_set(residue_options.element_set_tag());
	MMAtomTypeSetCOP mm_atom_type_set = ChemicalManager::get_instance()->mm_atom_type_set(residue_options.mm_atom_type_set_tag());
	orbitals::OrbitalTypeSetCOP orbital_type_set = ChemicalManager::get_instance()->orbital_type_set(residue_options.orbital_set_tag());
	ResidueTypeSetCOP residue_type_set = ChemicalManager::get_instance()->residue_type_set(residue_options.residue_type_set_tag());

	ResidueTypeOP new_residue_type(read_topology_file(locator_id,atom_type_set,element_set,mm_atom_type_set,orbital_type_set,residue_type_set));
	return new_residue_type;

}

basic::resource_manager::ResourceOptionsOP ResidueLoader::default_options() const
{
	return basic::resource_manager::ResourceOptionsOP( new ResidueLoaderOptions() );
}

}
}
