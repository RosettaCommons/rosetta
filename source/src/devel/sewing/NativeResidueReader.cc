// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file NativeResidueReader.cc
///
/// @brief

/// @author Tim jacobs

//Unit
#include <devel/sewing/NativeResidueReader.hh>

//Core
#include <core/types.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

//Utility
#include <utility/io/izstream.hh>
#include <utility/file/FileName.hh>
#include <utility/string_util.hh>

//Basic
#include <basic/Tracer.hh>

//C++
#include <string>
#include <map>
#include <stdio.h>

namespace devel {
namespace sewing {

static THREAD_LOCAL basic::Tracer TR( "NativeResidueReader" );

std::map<core::Size, utility::vector1<core::conformation::ResidueOP> >
NativeResidueReader::generateResiduesFromFile(utility::file::FileName file){

	utility::io::izstream input_stream(file);

	//generate residue type set from FA_STANDARD
	core::chemical::ResidueTypeSetCOP residue_set
		( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );

	std::map<core::Size, utility::vector1<core::conformation::ResidueOP> > res_map;

	//Read all residue lines here
	int res_num = 0;
	//std::string res_name;
	core::conformation::ResidueOP new_rsd;
	while ( input_stream.good() ) {
		std::string line;
		input_stream.getline(line);

		utility::vector1<std::string> tokens = utility::split(line);

		TR.Debug << "Parsing: " << line << std::endl;
		if ( tokens.size() > 0 ) { //Allow blank lines and comments in the file
			if ( tokens[1]=="RESNUM" ) {
				res_num = utility::string2int(tokens[2]);
				TR << "Creating Residue object for resnum: " << res_num << std::endl;
			} else if ( tokens[1]=="RESIDUE" ) {
				//get ResidueType from NativeResidue file residue name
				core::chemical::ResidueType const & res_type = residue_set->name_map(tokens[2]);

				//create new residue
				new_rsd = core::conformation::ResidueFactory::create_residue(res_type);
				new_rsd->seqpos(res_num);
				res_map[res_num].push_back(new_rsd);

				TR << "Creating default " << res_type.name3() << " to be populated: " << new_rsd->natoms() << std::endl;
			} else if ( tokens[1]=="ATOM" ) {
				double offset = 1e-250; /// coordinates now double, so we can use _really_ small offset.

				core::Size atom_index(utility::string2int(tokens[2]));
				core::Vector xyz(utility::string2float(tokens[3]),utility::string2float(tokens[4]),utility::string2float(tokens[5]));
				new_rsd->atom(atom_index).xyz(xyz+offset);

				TR.Debug << "Modified xyz of atom with index " << atom_index << std::endl;
			}
		}
	}
	input_stream.close();


	for ( std::map<core::Size, utility::vector1<core::conformation::ResidueOP> >::const_iterator map_it = res_map.begin();
			map_it != res_map.end(); ++map_it ) {

		TR.Debug << "Resnum: " << map_it->first << std::endl;
		for ( core::Size p=1; p<=map_it->second.size(); p++ ) {
			TR.Debug << "Residue: " << map_it->second[p]->name3() << std::endl;
			TR.Debug << "Num atoms: " << map_it->second[p]->natoms() << std::endl;
		}
	}

	return res_map;
}

} //sewing namespace
} //devel namespace


//ResidueType &
//    nonconst_name_map( std::string const & name )
//
//ResidueOP new_rsd( ResidueFactory::create_residue( rsd_type ) );
//
////Fill in residue coords
//for ( ResidueCoords::const_iterator iter=xyz.begin(), iter_end=xyz.end(); iter!= iter_end; ++iter ) {
//                if ( new_rsd->has( local_strip_whitespace(iter->first) ) ) {
//                    // offsetting all coordinates by a small constant prevents problems with atoms located
//                    // at position (0,0,0).
//                    // This is a bit of a dirty hack but it fixes the major problem of reading in rosetta
//                    // pdbs which suually start at 0,0,0. However the magnitude of this offset is so small
//                    // that the output pdbs should still match input pdbs. hopefully. yes. aehm.
//
//                    double offset = 1e-250; /// coordinates now double, so we can use _really_ small offset.
//                    new_rsd->atom( local_strip_whitespace(iter->first) ).xyz( iter->second + offset );
//                }
//                //else runtime_assert( iter->first == " H " && rsd_type.is_terminus() ); // special casee
//            }
