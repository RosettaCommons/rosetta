// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file util.cc
///
/// @brief
/// @author Tim Jacobs

//Unit
#include <devel/sewing/util.hh>

//Utility
#include <utility/io/ozstream.hh>

namespace devel {
namespace sewing {

void
dump_native_residue_file(
	NativeRotamersMap native_residue_map,
	std::string filename
){
	utility::io::ozstream native_residue_file;
	native_residue_file.open(filename);
	
	for(NativeRotamersMap::const_iterator pos_it=native_residue_map.begin();
		pos_it != native_residue_map.end(); ++pos_it)
	{
		native_residue_file << "RESNUM " << pos_it->first << "\n";
		for(utility::vector1<core::conformation::ResidueOP>::const_iterator res_it=pos_it->second.begin();
			res_it != pos_it->second.end(); ++res_it)
		{
			native_residue_file << "RESIDUE " << (*res_it)->type().name() << "\n";
			for(core::Size i=1; i<=(*res_it)->natoms(); ++i)
			{
				native_residue_file << "ATOM " << i << " " <<
					(*res_it)->xyz(i).x() << " " <<
					(*res_it)->xyz(i).y() << " " <<
					(*res_it)->xyz(i).z() << " " <<
					"\n";
			}
		}
		native_residue_file << "\n";
	}
	native_residue_file.close();
}

}
}
