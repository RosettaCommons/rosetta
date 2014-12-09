// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/rmoretti/rosetta_retype_check.cc
/// @brief  Check the rosetta retyping typing scheme in Rosetta
/// @details  Load a params file, do rosetta atom typing, print results with comparison to input types.
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.hh>

#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/residue_io.hh>

#include <core/types.hh>

#include <devel/init.hh>

#include <utility/options/FileVectorOption.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#
#include <basic/Tracer.hh>

#include <fstream>

using namespace basic::options;

static thread_local basic::Tracer TR( "apps.rosetta_retype_check" );

int
main( int argc, char * argv [] )
{

try {

	  using namespace basic::options;
	  using namespace basic::options::OptionKeys;
	  using namespace core::chemical;

	devel::init(argc, argv);

	utility::options::FileVectorOption & fvec
		= basic::options::option[ basic::options::OptionKeys::in::file::extra_res_fa ];

	core::chemical::ChemicalManager * chem_mang = core::chemical::ChemicalManager::get_instance();
	core::chemical::AtomTypeSetCAP atom_types = chem_mang->atom_type_set("fa_standard");
	core::chemical::ElementSetCAP elements = chem_mang->element_set("default");
	core::chemical::MMAtomTypeSetCAP mm_atom_types = chem_mang->mm_atom_type_set("fa_standard");
	core::chemical::orbitals::OrbitalTypeSetCAP orbital_types = chem_mang->orbital_type_set("fa_standard");

	core::chemical::ResidueTypeSet res_set;

	TR << "Loading residue types " << std::endl;
	// We don't need to load all the residue types - just the extra_res_fa ones.
	// Grab each and go.

	for(core::Size i = 1, e = fvec.size(); i <= e; ++i) {
		utility::file::FileName fname = fvec[i];
		std::string filename = fname.name();

		core::chemical::ResidueTypeOP ref( read_topology_file(
				filename, atom_types, elements, mm_atom_types, orbital_types, &res_set ) );

		TR << "Typing residue type " << fname.base()<< std::endl;

		core::chemical::ElementMap emap;
		for( core::Size ii(1); ii <= ref->natoms(); ++ii ) {
			emap[ ref->vd_from_index( ii ) ] = ref->atom_type(ii).element();
		}
		ResidueType rsd( *ref );
		rsd.retype_atoms(emap);
		for( core::Size ii(1); ii <= ref->natoms(); ++ii ) {
			if( ref->atom_type(ii).name() != rsd.atom_type(ii).name() ) {
				TR << "Atom " << ref->atom_name(ii) << " element " << emap[ ref->vd_from_index( ii ) ] << " was " << ref->atom_type(ii).name() << " now " << rsd.atom_type(ii).name() << std::endl;
				TR << "<<<<<<<<<<< Difference >>>>>>>>>>>" << std::endl;
			}
		}


	}

	TR << "Done outputing typeinfo" << std::endl;

} catch ( utility::excn::EXCN_Base const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}

}

