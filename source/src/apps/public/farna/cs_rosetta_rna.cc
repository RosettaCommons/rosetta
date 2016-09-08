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

// libRosetta headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rna/chemical_shift/RNA_ChemicalShiftPotential.hh>
#include <core/sequence/util.hh>

#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/rna/BasePair.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/init/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/sequence/Sequence.hh>

#include <protocols/toolbox/AtomLevelDomainMap.hh>

//RNA stuff.
#include <protocols/farna/movers/RNA_Minimizer.fwd.hh>
#include <protocols/farna/movers/RNA_Minimizer.hh>
#include <protocols/farna/setup/RNA_DeNovoPoseInitializer.hh>
#include <protocols/farna/setup/RNA_DeNovoParameters.hh>
#include <protocols/farna/RNA_DeNovoProtocol.hh>
#include <protocols/farna/movers/RNA_LoopCloser.hh>

// C++ headers
#include <iostream>
#include <string>

// option key includes
#include <basic/database/open.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <utility/file/file_sys_util.hh>
#include <utility/excn/Exceptions.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/viewer/viewers.hh>
#include <utility/vector1.hh>

OPT_KEY( String, pdb )
OPT_KEY( String, mode )
OPT_KEY( String, params_file )

// OPT_KEY( Integer, In_example )
// OPT_KEY( IntegerVector, IV_example )

using namespace basic::options;
using namespace basic::options::OptionKeys;

using namespace core;
using namespace core::chemical;
using namespace core::conformation;
using namespace core::kinematics;
using namespace core::id;
using namespace core::io::silent;
using namespace core::pose;
using namespace core::pose::rna;
using namespace core::scoring;
using namespace core::scoring::rna::chemical_shift;

using namespace protocols;
using namespace protocols::farna;
using namespace protocols::farna::setup;

using namespace ObjexxFCL;

using namespace std;

using utility::vector1;
using io::pdb::dump_pdb;

////////////////////////////////////////////////////////////////////////////////
// Helper functions

// Create the CS-ROSETTA-RNA score function.
// Note
// ----
// This is the score function used in the CS-ROSETTA-RNA study:
//   "Structure determination of noncanonical RNA motifs guided by 1H NMR
//    chemical shifts" (2014).
//
// This looks to have been written by P. Sripakdeevong in 2013-2014;
//  modified in 2015 to work with 'renovated' FARFAR code by R. Das.
//
// NOTE: We should get rid of this, which copies code with rna_denovo,
//  and instead use the 'normal' FARFAR workflow, which now has some nice
//  bells and whistles. See rna_denovo.cc -- and possibly move its rna_denovo_test()
//  function into protocols/farna/setup.cc or something, shared with this function.
//

ScoreFunctionOP
create_cs_rosetta_rna_scorefxn(){

	string score_weight_file = "stepwise/rna/rna_hires_07232011_with_intra_base_phosphate.wts";

	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( score_weight_file );
	//scorefxn->set_weight(rna_chem_shift, 4.0);

	cout << "---------CS-ROSETTA-RNA score function weights----------" << endl;
	scorefxn->show(cout);
	cout << "--------------------------------------------------------" << endl;

	return scorefxn;
}

// Get all obligated base pairs
utility::vector1< BasePair >
get_obligated_rna_pairings(RNA_DeNovoParameters const & rna_struct_params) {

	utility::vector1 < utility::vector1 <core::Size > >
		obligate_pairing_sets = rna_struct_params.obligate_pairing_sets();

	utility::vector1< BasePair >
		all_rna_pairings = rna_struct_params.rna_pairing_list();

	utility::vector1< BasePair > obl_rna_pairings;

	for ( Size i = 1; i <= obligate_pairing_sets.size(); i++ ) {
		for ( Size j = 1; j <= obligate_pairing_sets[i].size(); j++ ) {
			BasePair rna_pairing = all_rna_pairings[ obligate_pairing_sets[i][j] ];
			obl_rna_pairings.push_back(rna_pairing);

			cout << endl;
			cout << "cs_rosetta_rna:: obl_rna_pairing.res1() = " << rna_pairing.res1();
			cout << " obl_rna_pairing.res2() = " << rna_pairing.res2();
			cout << endl;

		}
	}

	return obl_rna_pairings;

}

// Setup RNA_MinimizerOP (for minimizing pose).
// Note
// ----
// Reuses existing functionality of FARNA.
RNA_MinimizerOP
setup_rna_minimizer(ScoreFunctionOP scorefxn,
	RNA_DeNovoParameters const & rna_struct_params,
	core::pose::Pose & pose ){

	RNA_MinimizerOP rna_minimizer( new RNA_Minimizer );
	//  rna_minimizer->vary_bond_geometry( false ); // defaults to false.
	rna_minimizer->set_score_function( scorefxn );
	//return rna_minimizer;

	utility::vector1< BasePair >
		obl_rna_pairings = get_obligated_rna_pairings(rna_struct_params);

	// Determine DOFs to allow minimization.
	toolbox::AtomLevelDomainMapOP allow_minimize ( new toolbox::AtomLevelDomainMap( pose ) );

	if ( obl_rna_pairings.size() > 0 ) {
		allow_minimize->set( false );
		for ( Size res_num = 1; res_num <= pose.size(); res_num++ ) {

			bool is_obl_res = false;
			for ( Size id = 1; id <= obl_rna_pairings.size(); id++ ) {
				if ( res_num == obl_rna_pairings[id].res1() ||
						res_num == obl_rna_pairings[id].res2() ) {
					is_obl_res = true;
				}
			}
			// Allow allow everywhere except obligated rna_pairing nucleotides.
			if ( is_obl_res ) continue;

			allow_minimize->set( res_num, true );
			// new -- make sure loops are moveable (& closeable!) at the 3'-endpoint.
			if ( res_num < pose.size() ) {
				allow_minimize->set_phosphate( res_num+1, pose, true );
			}
		}
	} else {
		allow_minimize->set( true );
	}

	cout << endl;
	cout << "allow_minimize DOFs:" << endl;
	allow_minimize->show();
	cout << endl;

	rna_minimizer->set_atom_level_domain_map( allow_minimize );

	return rna_minimizer;

}

// Compute and return the chem_shift_RMSD value (in ppm unit). Chem_shift_RMSD
// is the root-mean-deviation between the 'back-calculated' and the experimental
// 1H chemical shift. A low chem_shit_RMSD indicates that the RNA 3D structure
// agrees well with the experimental 1H chemical shift data.
Real
compute_chem_shift_RMSD(Pose const & pose) {

	if ( !option[ score::rna_chemical_shift_exp_data ].user() ) {
		utility_exit_with_message("-score::rna_chemical_shift_exp_data option required!");
	}

	Pose chem_shift_pose = pose;

	string score_weight_file = "stepwise/rna/rna_hires_07232011_with_intra_base_phosphate.wts";
	ScoreFunctionOP rna_chem_shift_scorefxn_ = ScoreFunctionFactory::create_score_function( score_weight_file );
	rna_chem_shift_scorefxn_->set_weight(rna_chem_shift, 1.0);

	/* Parin's Original
	ScoreFunctionOP rna_chem_shift_scorefxn_ = new ScoreFunction;
	rna_chem_shift_scorefxn_->set_weight( scoring::rna_chem_shift , 1 );
	*/

	(*rna_chem_shift_scorefxn_)(chem_shift_pose);
	EnergyMap const & energy_map = chem_shift_pose.energies().total_energies();
	Real const rosetta_chem_shift_score= energy_map[ scoring::rna_chem_shift ];

	RNA_ChemicalShiftPotential const & rna_cs_potential( ScoringManager::get_instance()->get_RNA_ChemicalShiftPotential() );
	Size const num_chem_shift_data_points = rna_cs_potential.get_total_exp_chemical_shift_data_points();

	Real const chem_shift_RMSD = sqrt( rosetta_chem_shift_score /
		float(num_chem_shift_data_points) );

	return chem_shift_RMSD;
}

// Apply CCD closer closer to all cutpoint closes positions.
void
close_chain_breaks(Pose & pose) {

	RNA_LoopCloser rna_loop_closer;
	rna_loop_closer.fast_scan( false );
	rna_loop_closer.apply( pose );

}
////////////////////////////////////////////////////////////////////////////////

// This function minimizes the user's specified PDB using the CS-ROSETTA-RNA
// score function. The function then print the PDB's CS-ROSETTA-RNA score
// and outputs the minimized PDB at <pdb>_out
//
// If perform_minimize is false, then the minimization is not performed. The
// function will still print the PDB's CS-ROSETTA-RNA score and outputs the PDB
// at <pdb>_out
//
// Required cmdline options
// ------------------------
// -pdb : Input PDB
// -score:rna_chemical_shift_exp_data : Experimental 1H CS data file
//
// Optional cmdline options
// ------------------------
// -param_file : Param file in FARNA format (specify fold_tree, cutpoints, insert_res)
//               If Param file is specified, then pose will have a simple fold_tree,
//               and no cutpoints. Also will minimize pose over every nucleotides.
void
cs_rosetta_rna_pdb(bool perform_minimize)
{

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	ScoreFunctionOP scorefxn = create_cs_rosetta_rna_scorefxn();

	if ( !option[ pdb ].user() ) {
		utility_exit_with_message("-pdb option required!");
	}
	string pdb_file  = option[ pdb ]();

	if ( !option[ score::rna_chemical_shift_exp_data ].user() ) {
		utility_exit_with_message("-score::rna_chemical_shift_exp_data option required!");
	}

	pose::Pose pose;
	core::import_pose::pose_from_file( pose, *rsd_set, pdb_file , core::import_pose::PDB_file);

	protocols::viewer::add_conformation_viewer( pose.conformation(),
		"current", 600, 600 );
	// Setup fold_tree, jump points, variants and chain breaks.
	RNA_DeNovoPoseInitializer rna_denovo_pose_setup( option[ params_file ]() );
	rna_denovo_pose_setup.initialize_for_de_novo_protocol( pose );
	rna_denovo_pose_setup.setup_fold_tree_and_jumps_and_variants( pose );

	// Setting fold_tree appear to introduce inproperly closed chain breaks.
	// Apply CCD chain closure to these chain break positions.
	close_chain_breaks( pose );

	if ( perform_minimize ) {
		RNA_MinimizerOP rna_minimizer = setup_rna_minimizer( scorefxn, rna_denovo_pose_setup.rna_params(), pose);
		rna_minimizer->apply( pose );
	} else { /*scores PDB*/
		cout << "CS-ROSETTA-RNA score for PDB (" << pdb_file << "):" << endl;
		scorefxn->show(std::cout, pose);
	}

	Real const total_score = (*scorefxn)( pose );

	Real const chem_shift_RMSD = compute_chem_shift_RMSD( pose );

	cout << endl;
	cout << "hybrid_CS-ROSETTA-RNA_all-atom energy: " << total_score << endl;
	cout << endl;
	cout << "chem_shift_RMSD: " << chem_shift_RMSD << endl;
	cout << endl;


	size_t found = pdb_file.find_last_of("/\\");
	string out_pdb_file;

	if ( found != string::npos ) {
		string basename = pdb_file.substr(found+1);
		out_pdb_file = basename + "_out";
	} else {
		out_pdb_file = pdb_file + "_out";
	}
	dump_pdb( pose, out_pdb_file );

}

///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	string user_mode = option[ mode ]();

	if ( user_mode == "score_pdb" ) {
		bool perform_minimize = false;
		cs_rosetta_rna_pdb(perform_minimize);
	} else if ( user_mode == "minimize_pdb" ) {
		bool perform_minimize = true;
		cs_rosetta_rna_pdb(perform_minimize);
	} else {
		utility_exit_with_message("Invalid user specified mode ("+ user_mode +")!");
	}

	protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

		std::cout << std::endl << "Basic usage:  " << argv[0] << "  -mode <mode>  -pdb <pdb_file>  -score:rna_chemical_shift_exp_data <cs_data_file> " << std::endl;
		std::cout << std::endl << " Type -help for full slate of options." << std::endl << std::endl;

		utility::vector1< Size > blank_size_vector;

		NEW_OPT( pdb, "Path to PDB file", "" );
		NEW_OPT( mode, "Mode to run", "minimize_pdb" );
		NEW_OPT( params_file, "Param file (in FARNA format)", "" );
		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////
		core::init::init(argc, argv);

		////////////////////////////////////////////////////////////////////////////
		// end of setup
		////////////////////////////////////////////////////////////////////////////

		protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

