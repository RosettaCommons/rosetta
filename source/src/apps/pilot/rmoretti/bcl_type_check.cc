// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   apps/pilot/rmoretti/gasteiger_type_check.cc
/// @brief  Check the gasteiger typing scheme in Rosetta
/// @details  Load a params file, do BCL atom typing, print results, and then output sdf file for comparison with BCL types.
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <core/chemical/gasteiger/GasteigerAtomTyper.hh>

#include <core/chemical/gasteiger/GasteigerAtomTypeSet.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.hh>

#include <core/chemical/ResidueType.hh>
#include <core/chemical/residue_io.hh>
#include <core/chemical/sdf/mol_writer.hh>

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

static basic::Tracer TR( "apps.gasteiger_type_check" );

int
main( int argc, char * argv [] )
{

	try {

		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		devel::init(argc, argv);

		utility::options::FileVectorOption & fvec
			= basic::options::option[ basic::options::OptionKeys::in::file::extra_res_fa ];

		core::chemical::ChemicalManager * chem_mang = core::chemical::ChemicalManager::get_instance();
		core::chemical::AtomTypeSetCAP atom_types = chem_mang->atom_type_set("fa_standard");
		core::chemical::ElementSetCAP elements = chem_mang->element_set("default");
		core::chemical::MMAtomTypeSetCAP mm_atom_types = chem_mang->mm_atom_type_set("fa_standard");
		core::chemical::orbitals::OrbitalTypeSetCAP orbital_types = chem_mang->orbital_type_set("fa_standard");
		core::chemical::gasteiger::GasteigerAtomTypeSetCOP gasteiger_atom_type_set = chem_mang->gasteiger_atom_type_set("default");

		std::ofstream datafile( "rosetta.types" );
		core::chemical::sdf::MolWriter molwriter;

		TR << "Loading residue types " << std::endl;
		// We don't need to load all the residue types - just the extra_res_fa ones.
		// Grab each and go.

		for ( core::Size i = 1, e = fvec.size(); i <= e; ++i ) {
			utility::file::FileName fname = fvec[i];
			std::string filename = fname.name();

			core::chemical::ResidueTypeOP restype( read_topology_file(
				filename, atom_types, elements, mm_atom_types, orbital_types ) );

			TR << "Typing residue type " << fname.base()<< std::endl;

			core::chemical::gasteiger::assign_gasteiger_atom_types( *restype, gasteiger_atom_type_set, false);

			//TR << "Outputting residue type " << fname.base()<< std::endl;

			for ( core::Size ii(1); ii <= restype->natoms(); ++ii ) {
				if ( ! restype->atom(ii).gasteiger_atom_type() ) {
					TR.Warning << "Untyped atom: " << restype->name() << " " << restype->atom(ii).name() << std::endl;
					continue;
				}
				datafile << "TYPE " << restype->name()
					//<< " " << restype->atom(ii).name()  //The sdf output doesn't include atom names
					<< " " << ii
					<<" : " << restype->atom(ii).gasteiger_atom_type()->get_name() << std::endl;
			}

			std::string outfile = fname.base() + "_post.sdf";
			//TR << "Outputing SDF file " << outfile << std::endl;

			molwriter.output_residue(outfile, restype);

			//TR << "Done with sdf output." << std::endl;
		}

		datafile.close();

		TR << "Done outputing typeinfo" << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

