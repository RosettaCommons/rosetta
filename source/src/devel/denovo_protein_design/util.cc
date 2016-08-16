// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @details
/// @author Grant Murphy


// Package Headers
#include <core/init/init.hh>

// Project Headers
#include <core/pose/Pose.hh>

#include <core/import_pose/import_pose.hh>

#include <basic/options/option.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Residue.hh>


#include <core/conformation/Residue.fwd.hh>


#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>


#include <core/scoring/dssp/Dssp.hh>// dssp info

//#include <protocols/simple_moves/BackboneMover.hh> //Small/ShearMover
#include <basic/options/keys/DenovoProteinDesign.OptionKeys.gen.hh>
#include <devel/denovo_protein_design/util.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/exit.hh>
#include <numeric/random/random.hh>

// C++ Headers
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>


namespace devel {
namespace denovo_protein_design {

void design_setup( core::pose::Pose & pose, core::pack::task::TaskFactoryOP designtaskfactory ){

	if ( basic::options::option[ basic::options::OptionKeys::DenovoProteinDesign::redesign_complete ].value() == true ) {
		basic::options::option[ basic::options::OptionKeys::DenovoProteinDesign::redesign_core ].value(true);
		basic::options::option[ basic::options::OptionKeys::DenovoProteinDesign::redesign_loops ].value(true);
		basic::options::option[ basic::options::OptionKeys::DenovoProteinDesign::redesign_surface ].value(true);
	}

	// create three vectors, containing the positions that are in the core, in loops, or on the surface
	utility::vector1< core::Size > CorePositions;
	utility::vector1< core::Size > LoopPositions;
	utility::vector1< core::Size > SurfPositions;

	// create two boolean vectors for designable and for repackable
	utility::vector1< bool > designablePositions(pose.n_residue(),false);
	utility::vector1< bool > packablePositions(pose.n_residue(),false); // packable doesn't imply designing

	// get dssp info
	core::scoring::dssp::Dssp dssp( pose );
	dssp.insert_ss_into_pose( pose );

	utility::vector1< core::Size > ind_neighbors; // neighbors at each position
	utility::vector1< utility::vector1< core::Size> > neighbors_; // holds neighbors of all positions


	for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		for ( core::Size jj = 1; jj <= pose.total_residue(); ++jj ) {
			core::Size currentposition = ii;
			core::conformation::Residue const & rsd1 = pose.residue( currentposition );
			core::Size nextneighborposition = jj;
			if ( currentposition == nextneighborposition ) { continue; }
			core::conformation::Residue const & rsd2 = pose.residue( nextneighborposition );
			core::Real distanceBetweenAtoms = rsd1.xyz( rsd1.nbr_atom() ).distance( rsd2.xyz( rsd2.nbr_atom() ) );
			if ( distanceBetweenAtoms <= 10 ) { ind_neighbors.push_back(jj); }
		}
		neighbors_.push_back( ind_neighbors );
		ind_neighbors.clear();
	}

	using namespace basic::options;
	for ( core::Size ii=1; ii <= pose.n_residue(); ++ii ) {
		if ( neighbors_[ii].size() >= 19 && pose.secstruct(ii) != 'L' ) {
			CorePositions.push_back(ii);
			if ( option[ OptionKeys::DenovoProteinDesign::redesign_core ].value() == true ) {
				designablePositions[ii] = true;
			}
		}
		if ( pose.secstruct(ii) == 'L' ) {
			LoopPositions.push_back(ii);
			if ( option[ OptionKeys::DenovoProteinDesign::redesign_loops ].value() == true ) {
				designablePositions[ii] = true;
			}
		}
		if ( neighbors_[ii].size() < 19 && pose.secstruct(ii) != 'L' ) {
			SurfPositions.push_back(ii);
			if ( option[ OptionKeys::DenovoProteinDesign::redesign_surface ].value() == true ) {
				designablePositions[ii] = true;
			}
		}
	}

	// residues that are not designable and are not neighbors of designables are held fixed

	// set true those positions that are packable
	for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		if ( designablePositions[ii] == true ) {
			packablePositions[ii] = true;
			for ( core::Size jj = 1; jj <= neighbors_[ii].size(); ++jj ) {
				packablePositions[ neighbors_[ii][jj] ] = true;
			}
		}
	}

	// designtaskfactory->set_bump_check( true );
	using namespace core::pack::task;
	operation::PreventRepackingOP preventrepacking( new operation::PreventRepacking );
	operation::RestrictResidueToRepackingOP allowrepacking( new operation::RestrictResidueToRepacking );


	for ( core::Size ii = 1; ii <= pose.n_residue(); ++ii ) {

		if ( packablePositions[ii] == false ) {  // if a position is false in the packablePositions vector then it is prevented from repacking

			preventrepacking->include_residue( ii );

		}

		if ( packablePositions[ii] == true && designablePositions[ii] == false ) { // position is only repackable not designable

			allowrepacking->include_residue( ii );

		}

		if ( designablePositions[ii] == true ) {

			//   core::Size aa(pose.residue_type(ii).aa()); // gives the Enum of the amino acid at res pos designPositions[ii]

			utility::vector1< bool > keep_aas( core::chemical::num_canonical_aas, true ); // default is set all amino acids true for this position

			// disallow native aa
			//   if (option[ OptionKeys::DenovoProteinDesign::disallow_native_aa ].value() == true) {
			//    keep_aas[aa] = false; // disallow native AAs at designable positions
			//   }

			operation::RestrictAbsentCanonicalAASOP allowdesign( new operation::RestrictAbsentCanonicalAAS( ii, keep_aas) );
			designtaskfactory->push_back( allowdesign );
		}

	}

	designtaskfactory->push_back(preventrepacking);
	designtaskfactory->push_back(allowrepacking);

}

core::Size numberhelices( core::pose::Pose & pose ){

	core::scoring::dssp::Dssp dssp( pose );
	dssp.insert_ss_into_pose( pose );

	utility::vector1< core::Size > Ind_SS_E;
	utility::vector1< utility::vector1< core::Size > > Total_SS_E;
	Total_SS_E.clear();

	for ( core::Size ii = 1; ii <= pose.n_residue(); ++ii ) {
		if ( ii != pose.n_residue() ) {
			if ( pose.secstruct(ii) == pose.secstruct(ii+1) ) {
				Ind_SS_E.push_back(ii);
			}
			if ( pose.secstruct(ii) != pose.secstruct(ii+1) ) {
				Ind_SS_E.push_back(ii); // this gets the last residue in the ss element into the list
				Total_SS_E.push_back(Ind_SS_E);
				Ind_SS_E.clear();
			}
		}
		if ( ii == pose.n_residue() ) {
			if ( pose.secstruct(ii) == pose.secstruct(ii-1) ) { Ind_SS_E.push_back(ii); Total_SS_E.push_back(Ind_SS_E); }
			if ( pose.secstruct(ii) != pose.secstruct(ii-1) ) {
				Total_SS_E.push_back(Ind_SS_E); // if current_pose.secstruct[i] != current_pose.secstruct[i-1] then the last residue is a new SS Element so pushback previous element
				Ind_SS_E.clear(); // clear the vector of positions of individual SS Elements
				Ind_SS_E.push_back(ii); // add last residue to a new individual SS Elements
				Total_SS_E.push_back(Ind_SS_E); // and push back that new individual SS Element
			}
		}
	}


	core::Size numHELIX(0);
	utility::vector1< core::Size > HELICES;
	HELICES.clear();


	for ( core::Size ii =1; ii<= Total_SS_E.size(); ++ii ) {
		if ( Total_SS_E[ii].size() > 0 ) {
			if ( pose.secstruct(Total_SS_E[ii][1]) == 'H' ) { numHELIX++; HELICES.push_back( Total_SS_E[ii].size() ); }
		}
	}

	return numHELIX;
}

void create_nucleated_sequence_from_template_pdb( std::string & nucleated_sequence , utility::vector1< bool > & KeepNativeTorsions ){
	core::pose::Pose template_pose;
	std::string template_pdb = basic::options::option[ basic::options::OptionKeys::DenovoProteinDesign::create_from_template_pdb ]().name();
	core::import_pose::pose_from_file( template_pose, template_pdb , core::import_pose::PDB_file);

	// get the template poses secondary structure information
	core::scoring::dssp::Dssp dssp( template_pose );
	dssp.insert_ss_into_pose( template_pose);

	// get nucleated sequence from template pose B/P pattern or use template sequence
	if ( basic::options::option[ basic::options::OptionKeys::DenovoProteinDesign::use_template_sequence ].value() == true ) {
		nucleated_sequence = template_pose.sequence();
	} else {
		for ( core::Size pos = 1; pos <= template_pose.n_residue(); pos++ ) {

			if (  (template_pose.residue( pos ).is_polar()) ) {  nucleated_sequence += devel::denovo_protein_design::random_polar_aa();  }
			if ( !(template_pose.residue( pos ).is_polar()) ) {  nucleated_sequence += devel::denovo_protein_design::random_apolar_aa(); }
		}
	}

	if ( basic::options::option[ basic::options::OptionKeys::DenovoProteinDesign::use_template_topology ].value() == true ) {
		KeepNativeTorsions.resize( template_pose.n_residue(), false );
		for ( core::Size seqpos = 1; seqpos <= template_pose.total_residue(); ++seqpos ) {

			if ( template_pose.secstruct( seqpos ) == 'L' || template_pose.secstruct( seqpos - 1 ) == 'L' || template_pose.secstruct( seqpos + 1) == 'L' ) {
				KeepNativeTorsions[seqpos] = true;
			}

		}
	}
}

void read_ss_file_into_vector( std::string const & filename,  utility::vector1< char > & secstructs  ){

	utility::io::izstream data( filename );
	if ( !data ) {
		std::cout << "can not open secondary structure file " << filename << std::endl;
	}

	std::string line;
	char sec;
	getline(data, line); // skip the > line
	while ( getline(data,line) ) {
		std::istringstream line_stream( line );
		while ( line_stream >> sec ) {

			if ( sec != 'H' && sec != 'E' && sec != 'C' ) {
				std::cout << "unrecognized secondary structure element : " << sec << std::endl;
			}
			if ( sec == 'C' ) {
				secstructs.push_back( 'L' );
			} else {
				secstructs.push_back( sec );
			}
		}
	}
	return;
}


void read_hp_file_into_vector( std::string const &filename,  utility::vector1< char > & hydrophobic_polar_sequence_v ){

	utility::io::izstream data( filename );
	if ( !data ) {
		std::cout << "can not open secondary structure file " << filename << std::endl;
	}

	std::string line;
	char sec;
	getline(data, line); // skip the > line
	while ( getline(data,line) ) {
		std::istringstream line_stream( line );
		while ( line_stream >> sec ) {

			if ( sec != 'P' && sec != 'B' ) {
				std::cout << "unrecognized hydrophobic polar element, options are B and P : " << sec << std::endl;
			}
			hydrophobic_polar_sequence_v.push_back( sec );
		}
	}
}

void create_hydrophobic_polar_pattern_from_secstruct_vector( utility::vector1< char > & hydrophobic_polar_sequence_v , utility::vector1< char > const & secstructs ){
	utility::vector1< core::Size > Ind_SS_E;
	utility::vector1< utility::vector1< core::Size > > Total_SS_E;

	for ( core::Size ii = 1; ii <= secstructs.size(); ++ii ) {

		if ( ii != secstructs.size() ) {
			if ( secstructs[ii] == secstructs[ii+1] ) {
				Ind_SS_E.push_back(ii);
			}

			if ( secstructs[ii] != secstructs[ii+1] ) {
				Ind_SS_E.push_back(ii); // this gets the last residue in the ss element into the list
				Total_SS_E.push_back(Ind_SS_E);
				Ind_SS_E.clear();
			}
		}

		if ( ii == secstructs.size() ) {
			if ( secstructs[ii] == secstructs[ii-1] ) { Ind_SS_E.push_back(ii); Total_SS_E.push_back(Ind_SS_E); }
			if ( secstructs[ii] != secstructs[ii-1] ) {
				Total_SS_E.push_back(Ind_SS_E); // if secstructs[i] != secstructs[i-1] then the last residue is a new SS Element so pushback previous element
				Ind_SS_E.clear(); // clear the vector of positions of individual SS Elements
				Ind_SS_E.push_back(ii); // add last residue to a new individual SS Elements
				Total_SS_E.push_back(Ind_SS_E); // and push back that new individual SS Element
			}
		} // this handles the last residue, we need to look back instead of forward( because we will fall off the vector )

	}


	hydrophobic_polar_sequence_v.resize(secstructs.size());

	for ( core::Size ii = 1; ii<=Total_SS_E.size(); ++ii ) {
		for ( core::Size jj = 1; jj <= Total_SS_E[ii].size(); ++jj ) {


			if ( secstructs[Total_SS_E[ii][jj]] == 'L' ) { hydrophobic_polar_sequence_v[Total_SS_E[ii][jj]] = random_helix_position(); }

			if ( secstructs[Total_SS_E[ii][jj]] == 'E' ) {
				if ( jj == 1 ) { hydrophobic_polar_sequence_v[Total_SS_E[ii][jj]] = random_hydrophobic_polar();}
				if ( jj != 1 ) {
					if ( hydrophobic_polar_sequence_v[Total_SS_E[ii][jj]-1] == 'B' ) { hydrophobic_polar_sequence_v[Total_SS_E[ii][jj]] = 'P'; } else { hydrophobic_polar_sequence_v[Total_SS_E[ii][jj]] = 'B'; }
				}
			} // end if == 'E'

			if ( secstructs[Total_SS_E[ii][jj]] == 'H' ) {
				if ( jj == 1 ) { hydrophobic_polar_sequence_v[Total_SS_E[ii][jj]] = random_helix_position();}
				if ( jj != 1 ) {
					if ( hydrophobic_polar_sequence_v[Total_SS_E[ii][jj] - 1 ] == 'a' ) { hydrophobic_polar_sequence_v[Total_SS_E[ii][jj]] = 'b'; }
					if ( hydrophobic_polar_sequence_v[Total_SS_E[ii][jj] - 1 ] == 'b' ) { hydrophobic_polar_sequence_v[Total_SS_E[ii][jj]] = 'c'; }
					if ( hydrophobic_polar_sequence_v[Total_SS_E[ii][jj] - 1 ] == 'c' ) { hydrophobic_polar_sequence_v[Total_SS_E[ii][jj]] = 'd'; }
					if ( hydrophobic_polar_sequence_v[Total_SS_E[ii][jj] - 1 ] == 'd' ) { hydrophobic_polar_sequence_v[Total_SS_E[ii][jj]] = 'e'; }
					if ( hydrophobic_polar_sequence_v[Total_SS_E[ii][jj] - 1 ] == 'e' ) { hydrophobic_polar_sequence_v[Total_SS_E[ii][jj]] = 'f'; }
					if ( hydrophobic_polar_sequence_v[Total_SS_E[ii][jj] - 1 ] == 'f' ) { hydrophobic_polar_sequence_v[Total_SS_E[ii][jj]] = 'g'; }
					if ( hydrophobic_polar_sequence_v[Total_SS_E[ii][jj] - 1 ] == 'g' ) { hydrophobic_polar_sequence_v[Total_SS_E[ii][jj]] = 'a'; }
				}
			} // end if == 'H'
		}
	}


}

void create_sequence_from_hydrophobic_polar_pattern( std::string & sequence, utility::vector1< char > const & hydrophobic_polar_sequence_v){
	for ( core::Size ii = 1; ii <= hydrophobic_polar_sequence_v.size(); ii++ ) {

		if ( hydrophobic_polar_sequence_v[ii] == 'P'  ) { sequence += random_polar_aa();  }
		if ( hydrophobic_polar_sequence_v[ii] == 'B'  ) { sequence += random_apolar_aa(); }
		if ( hydrophobic_polar_sequence_v[ii] == 'a'  ) { sequence += random_apolar_aa(); }
		if ( hydrophobic_polar_sequence_v[ii] == 'b'  ) { sequence += random_polar_aa(); }
		if ( hydrophobic_polar_sequence_v[ii] == 'c'  ) { sequence += random_polar_aa(); }
		if ( hydrophobic_polar_sequence_v[ii] == 'd'  ) { sequence += random_apolar_aa(); }
		if ( hydrophobic_polar_sequence_v[ii] == 'e'  ) { sequence += random_apolar_aa(); }
		if ( hydrophobic_polar_sequence_v[ii] == 'f'  ) { sequence += random_polar_aa(); }
		if ( hydrophobic_polar_sequence_v[ii] == 'g'  ) { sequence += random_apolar_aa(); }
	}


}

char random_hydrophobic_polar(){

	utility::vector1< char > hp;
	hp.clear();
	hp.push_back('B');
	hp.push_back('P');
	core::Size randnum = numeric::random::random_range( 1 , hp.size());
	return hp[randnum];
}

char random_polar_aa(){

	utility::vector1< char > polarAA;
	polarAA.clear();
	polarAA.push_back( 'D' );
	polarAA.push_back( 'E' );
	polarAA.push_back( 'K' );
	polarAA.push_back( 'N' );
	polarAA.push_back( 'Q' );
	polarAA.push_back( 'R' );
	polarAA.push_back( 'S' );
	polarAA.push_back( 'T' );
	polarAA.push_back( 'Y' );
	core::Size randnum = numeric::random::random_range( 1 , polarAA.size());
	return polarAA[randnum];
}

char random_apolar_aa(){
	utility::vector1< char > apolarAA;
	apolarAA.clear();
	apolarAA.push_back( 'A' );
	apolarAA.push_back( 'C' );
	apolarAA.push_back( 'F' );
	apolarAA.push_back( 'H' );
	apolarAA.push_back( 'I' );
	apolarAA.push_back( 'L' );
	apolarAA.push_back( 'M' );
	apolarAA.push_back( 'V' );
	apolarAA.push_back( 'W' );
	core::Size randnum = numeric::random::random_range( 1 , apolarAA.size());
	return apolarAA[randnum];
}

char random_helix_position(){
	utility::vector1< char > helixposition;
	helixposition.clear();
	helixposition.push_back( 'a' );
	helixposition.push_back( 'b' );
	helixposition.push_back( 'c' );
	helixposition.push_back( 'd' );
	helixposition.push_back( 'e' );
	helixposition.push_back( 'f' );
	helixposition.push_back( 'g' );
	core::Size randnum = numeric::random::random_range( 1 , helixposition.size());
	return helixposition[randnum];
}


}
}
