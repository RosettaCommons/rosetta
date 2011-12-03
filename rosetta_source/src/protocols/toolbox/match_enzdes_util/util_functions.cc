// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/toolbox/match_enzdes_util/util_functions.cc
/// @brief bunch of utility functions
/// @author Florian Richter, floric@u.washington.edu

//unit headers
#include <protocols/toolbox/match_enzdes_util/util_functions.hh>

// project headers
#include <basic/Tracer.hh>

#include <core/types.hh>

#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>


//utility headers
#include <utility/string_util.hh>

//stl headers
#include <map>

namespace protocols {
namespace toolbox {
namespace match_enzdes_util{

static basic::Tracer tr("protocols.toolbox.match_enzdes_util.util_functions");

void
replace_residue_keeping_all_atom_positions(
	core::pose::Pose & pose,
	core::conformation::Residue new_res,
	core::Size res_pos
)
{

	//have to set the position of the new res to their old values, so we gotta save them now
	std::map< std::string, core::PointPosition > atom_name_to_xyz;

	for( core::Size at_ct = 1; at_ct <= pose.residue(res_pos).natoms(); at_ct++){
		atom_name_to_xyz.insert( 	std::pair< std::string, core::PointPosition > (pose.residue(res_pos).atom_name(at_ct), pose.residue(res_pos).xyz( at_ct ) ) );
	}

	//replacing the residue
	pose.replace_residue( res_pos, new_res, true);

	//and resetting the xyz positions
	for( core::Size at_ct = 1; at_ct <= pose.residue(res_pos).natoms(); at_ct++){

		std::map< std::string, core::PointPosition>::iterator xyz_map_it = atom_name_to_xyz.find( pose.residue(res_pos).atom_name(at_ct) );

		if(xyz_map_it == atom_name_to_xyz.end() ) {
			std::cerr << "ERROR: when trying to make dsflkj constraint covalent, atom " << pose.residue(res_pos).atom_name(at_ct) << " was not found for residue " << pose.residue(res_pos).name3() << " at position " << res_pos << std::endl;
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}
		else{
			pose.set_xyz( core::id::AtomID (at_ct, res_pos), xyz_map_it->second );
		}
	}

} //replace_residues_keeping_positions


std::string
assemble_remark_line(
	std::string chainA,
	std::string resA,
	int seqposA,
	std::string chainB,
	std::string resB,
	int seqposB,
	core::Size cst_block,
	core::Size ex_geom_id
)
{
	std::string posA = utility::to_string( seqposA );
	utility::add_spaces_right_align( posA, 4 );

	std::string posB = utility::to_string( seqposB );
	utility::add_spaces_right_align( posB, 4 );

	return "MATCH TEMPLATE "+ chainA +" "+ resA +" "+ posA +  " MATCH MOTIF "+ chainB + " " + resB + " "+posB + "  " + utility::to_string( cst_block ) + "  " + utility::to_string( ex_geom_id );

} //assemble remark line function


bool
split_up_remark_line(
	std::string line,
	std::string & chainA,
	std::string & resA,
	int & seqposA,
	std::string & chainB,
	std::string & resB,
	int & seqposB,
	core::Size & cst_block,
	core::Size & ex_geom_id
)
{

	std::istringstream line_stream;
	std::string buffer(""), tag("");

	line_stream.clear();
	line_stream.str( line );

	line_stream >> buffer >> tag;
	if( tag == "TEMPLATE"){
		line_stream >> chainA >> resA >> seqposA >> buffer >> buffer;
		line_stream >> chainB >> resB >> seqposB >> cst_block;
		if( resA.size() == 2 ) resA = " " + resA;
		if( resB.size() == 2 ) resB = " " + resB;

		if( !line_stream.good() ){
			tr << "ERROR when trying to split up pdb remark line. Not all fields seem to have been specified." << std::endl;
			return false;
		}

		line_stream >> ex_geom_id;
		if( !line_stream.good() ) ex_geom_id = 1;

		return true;
	}

	return false;
}  //split up remark line function

}
} // enzdes
} //protocols

