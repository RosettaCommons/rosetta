// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief
/// @author Ragul Gowthaman

// Project Headers
#include <devel/init.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/after_opts.hh>
//Auto Headers
#include <core/io/pdb/pose_io.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <protocols/pockets/GenPharmacophore.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/rna/RNA_ResidueType.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

// C++ Headers
#include <iomanip>
#include <fstream>
#include <ostream>
#include <sstream>
#include <iostream>
#include <string>

using namespace std;
using namespace core;
using namespace core::pose::datacache;
using namespace core::optimization;
using namespace core::pose::metrics;
using namespace core::scoring;
using namespace core::scoring::constraints;
using namespace core::id;
using namespace core::chemical;
using namespace core::chemical::rna;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::conformation;

OPT_KEY( String, input_rna )
OPT_KEY( String, input_protein )
OPT_KEY( String, pre_phr )


int main( int argc, char * argv [] ){

	try{

		NEW_OPT( input_rna, "rna file name", "rna.pdb" );
		NEW_OPT( input_protein, "rna protein name", "protein.pdb" );
		NEW_OPT( pre_phr, "pre clustering file name", "" );

		devel::init(argc, argv);

		protocols::pockets::GenPharmacophore rphr;
		std::string const pre_phr_pose = option[ pre_phr ];
		std::string const input_rna_pose = option[ input_rna ];
		std::string const input_protein_pose = option[ input_protein ];

		//cteate 'tag' for output files
		int pfounddir = input_protein_pose.find_last_of("/\\");
		int pfounddot = input_protein_pose.find_last_of(".");
		//get the basename of the protein pdb file
		std::string tag = input_protein_pose.substr((pfounddir+1),(pfounddot-(pfounddir+1)));

		if ( pre_phr_pose.size() > 0 ) {
			std::cout << "Clustering " << pre_phr_pose << std::endl;

			std::ifstream ifs(pre_phr_pose.c_str());
			std::string content( (std::istreambuf_iterator<char>(ifs) ),\
				(std::istreambuf_iterator<char>()    ));

			rphr.cluster_KeyFeatures(content, tag);
			return 0;
		}

		pose::Pose rna_pose, protein_pose;
		core::import_pose::pose_from_pdb( rna_pose, input_rna_pose );
		core::import_pose::pose_from_pdb( protein_pose, input_protein_pose );

		std::string keyFeatures_hbond = rphr.extract_Hbond_atoms_from_protein_rna_complex(protein_pose, rna_pose);
		//rphr.print_string_to_PDBfile(keyFeatures_hbond, "hbond.pdb");
		std::string keyFeatures_rings = rphr.extract_rna_rings_from_protein_rna_complex(protein_pose, rna_pose);
		//rphr.print_string_to_PDBfile(keyFeatures_rings, "rings.pdb");

		//concatenate hbond and rings
		//ring atoms should be placed first (clustering code requires like that!)
		std::string keyFeatures_string1 = keyFeatures_rings + keyFeatures_hbond;

		std::string precluster_PDBfile = tag + "_preclusterPHR.pdb";
		rphr.print_string_to_PDBfile(keyFeatures_string1, precluster_PDBfile);

		std::string keyFeatures_string2 = rphr.make_compatible_with_ROCS_custom_ForceField(keyFeatures_string1);
		rphr.cluster_KeyFeatures(keyFeatures_string2, tag);

	}
catch ( utility::excn::EXCN_Base const & e ) {
	std::cerr << "caught exception " << e.msg() << std::endl;
	return -1;
}

	return 0;

}
