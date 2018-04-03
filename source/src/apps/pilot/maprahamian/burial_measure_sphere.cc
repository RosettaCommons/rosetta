// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// // vi: set ts=2 noet:
// //
// // (c) Copyright Rosetta Commons Member Institutions.
// // (c) This file is part of the Rosetta software suite and is made available under license.
// // (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// // (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// // (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
// /// @file
// /// @brief Determines the number of residue neighbors within a sphere around each residue.
//     [-mid] midpoint of distance cutoff
//     [-exp_dist] distance falloff exponent factor
// /// @author Melanie Aprahamian

#include <iostream>
#include <string>
#include <protocols/jd2/util.hh>
#include <devel/init.hh>
#include <utility/excn/Exceptions.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>

#include <utility/io/izstream.hh>
#include <utility/file/FileName.hh>
#include <utility/vector1.hh>

#include <numeric/NumericTraits.hh>

using namespace core;
using namespace core::scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace pose;

OPT_KEY( Real, dist_midpoint )
OPT_KEY( Real, dist_exponent )
OPT_KEY( Boolean, FPOP_label )

// Main program
int main( int argc, char * argv [] )
{
	try
{

		// get option keys
		NEW_OPT( dist_midpoint, "midpoint of distance falloff", 9.0);
		NEW_OPT( dist_exponent, "exponent factor for distance falloff", 1.0);
		NEW_OPT( FPOP_label, "measure from FPOP label site", false);

		// Initialize Rosetta
		devel::init( argc, argv );

		// Namespaces
		using namespace core;
		using namespace core::scoring;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace pose;

		// check if input PDB file is provided from command line with option -s
		if ( !option[in::file::s].user() ) {
			// exit if no PDB file is found
			utility_exit_with_message("Input PDB file not found");
		}

		// get measure parameters
		core::Real m = option[ dist_midpoint ];
		core::Real n = option[ dist_exponent ];
		bool FPOP = option[ FPOP_label ];

		// read in and import pose
		pose::Pose p;
		std::string pdb_file = option[in::file::s]()[1];
		import_pose::pose_from_file(p,pdb_file);

		// iterate over all residues to determine neighbor counts
		core::Size num_residues (p.total_residue());

		for ( core::Size res_count_target = 1; res_count_target <= num_residues; res_count_target++ ) {
			char residue_name (p.residue(res_count_target).type().name1());

			if ( residue_name == 'C' || residue_name == 'M' || residue_name == 'W' || residue_name == 'Y' || residue_name == 'F' || residue_name == 'H' || residue_name == 'L' || residue_name == 'I' || residue_name == 'R' || residue_name == 'K' || residue_name == 'V' || residue_name == 'V' || residue_name == 'P' || residue_name == 'E' || residue_name == 'Q' ) {

				// define label site
				numeric::xyzVector<core::Real> label_site = p.residue(res_count_target).xyz("CB");

				if ( FPOP == true ) {
					if ( residue_name == 'C' ) {
						label_site = p.residue(res_count_target).xyz("SG");
					}
					if ( residue_name == 'M' ) {
						label_site = p.residue(res_count_target).xyz("SD");
					}
					if ( residue_name == 'W' ) {
						label_site = (p.residue(res_count_target).xyz("CD1") + p.residue(res_count_target).xyz("CE3") + p.residue(res_count_target).xyz("CZ2") + p.residue(res_count_target).xyz("CZ3") + p.residue(res_count_target).xyz("CH2"))/5.0;
					}
					if ( residue_name == 'Y' ) {
						label_site = (p.residue(res_count_target).xyz("CE1") + p.residue(res_count_target).xyz("CE2"))/2.0;
					}
					if ( residue_name == 'F' ) {
						label_site = (p.residue(res_count_target).xyz("CD1") + p.residue(res_count_target).xyz("CD2") + p.residue(res_count_target).xyz("CE1") + p.residue(res_count_target).xyz("CE2") + p.residue(res_count_target).xyz("CZ"))/5.0;
					}
					if ( residue_name == 'H' ) {
						label_site = p.residue(res_count_target).xyz("CE1");
					}
					if ( residue_name == 'L' ) {
						label_site = (p.residue(res_count_target).xyz("CD1") + p.residue(res_count_target).xyz("CD2"))/2.0;
					}
					if ( residue_name == 'I' ) {
						label_site = (p.residue(res_count_target).xyz("CG2") + p.residue(res_count_target).xyz("CD1"))/2.0;
					}
					if ( residue_name == 'R' ) {
						label_site = p.residue(res_count_target).xyz("CD");
					}
					if ( residue_name == 'K' ) {
						label_site = (p.residue(res_count_target).xyz("CG") + p.residue(res_count_target).xyz("CD"))/2.0;
					}
					if ( residue_name == 'V' ) {
						label_site = (p.residue(res_count_target).xyz("CG2") + p.residue(res_count_target).xyz("CG1"))/2.0;
					}
					if ( residue_name == 'P' ) {
						label_site = (p.residue(res_count_target).xyz("CG") + p.residue(res_count_target).xyz("CD") + p.residue(res_count_target).xyz("CB"))/3.0;
					}
					if ( residue_name == 'E' ) {
						label_site = (p.residue(res_count_target).xyz("CG") + p.residue(res_count_target).xyz("CB"))/2.0;
					}
					if ( residue_name == 'Q' ) {
						label_site = (p.residue(res_count_target).xyz("CG") + p.residue(res_count_target).xyz("CB"))/2.0;
					}
				}

				// use any atom for neighbor count
				core::Real neighbor_count (0.0);
				for ( core::Size res_count_neighbor = 1; res_count_neighbor <= num_residues; res_count_neighbor++ ) {
					if ( p.residue(res_count_target).seqpos() != p.residue(res_count_neighbor).seqpos() ) {
						std::string neighbor_atom (p.residue(res_count_neighbor).atom_name(1));
						// find the atom in neighbor residue that is closest to target CB
						for ( core::Size atom_count = 2; atom_count <= p.residue(res_count_neighbor).natoms(); atom_count++ ) {
							if ( p.residue(res_count_neighbor).xyz(atom_count).distance(p.residue(res_count_target).xyz("CB")) < p.residue(res_count_neighbor).xyz(1).distance(p.residue(res_count_target).xyz("CB")) ) {
								neighbor_atom = p.residue(res_count_neighbor).atom_name(atom_count);
							}
						}

						// compute the distance betwee target CA and neighbor atom
						core::Real distance = p.residue(res_count_neighbor).xyz(neighbor_atom).distance(label_site);
						// calculate weighted neighbor count per res
						neighbor_count += 1.0/(1.0 + std::exp(n*(distance-m)));
					}
				}
				std::cout << p.residue(res_count_target).type().name1() << " " << p.residue(res_count_target).seqpos() << " " << neighbor_count << std::endl;
			}
		}
	}
catch ( utility::excn::EXCN_Base const & e )
{
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}
	return 0;
} //main

