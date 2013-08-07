// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file
/// @brief


// libRosetta headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>

#include <core/init/init.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/util.hh>
#include <basic/options/option_macros.hh>

#include <protocols/viewer/viewers.hh>

#include <utility/excn/Exceptions.hh>


//////////////////////////////////////////////////
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh> //FCC:12/13/10
///////////////////////////////////////////////////

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/types.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
//////////////////////////////////////////////////////////


// C++ headers
#include <iostream>
#include <string>

using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using utility::vector1;

OPT_KEY ( String, force_field_file )

/////////////////////////////////////////////////////////////////////////////
//FCC: Adding Virtual res
void
add_virtual_res ( core::pose::Pose & pose, bool set_res_as_root ) {
	int nres = pose.total_residue();

	// if already rooted on virtual residue , return
	if ( pose.residue ( pose.fold_tree().root() ).aa() == core::chemical::aa_vrt ) {
		std::cout << "addVirtualResAsRoot() called but pose is already rooted on a VRT residue ... continuing." << std::endl;
		return;
	}

	// attach virt res there
	bool fullatom = pose.is_fullatom();
	core::chemical::ResidueTypeSet const & residue_set = pose.residue_type ( 1 ).residue_type_set();
	core::chemical::ResidueTypeCOPs const & rsd_type_list ( residue_set.name3_map ( "VRT" ) );
	core::conformation::ResidueOP new_res ( core::conformation::ResidueFactory::create_residue ( *rsd_type_list[1] ) );
	pose.append_residue_by_jump ( *new_res , 1 );

	// make the virt atom the root
	if ( set_res_as_root ) {
		kinematics::FoldTree newF ( pose.fold_tree() );
		newF.reorder ( nres + 1 );
		pose.fold_tree ( newF );
	}
}

void
pdb_scoring() {
	using namespace core::pose;
	using namespace core::chemical;
	using namespace core::kinematics;
	using namespace core::scoring;
	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->
	          residue_type_set ( RNA );
	// Read in the pose from pdb.
	PoseOP pose;

	if ( option[ in::file::native ].user() ) {
		pose = PoseOP ( new Pose );
		import_pose::pose_from_pdb ( *pose, *rsd_set, option[in::file::native]() );
	}

	add_virtual_res ( *pose, true );
	//Score the pose.
	std::string const force_field_file_option = option[force_field_file];
	core::scoring::ScoreFunctionOP scorefxn =
		ScoreFunctionFactory::create_score_function ( force_field_file_option );
	scorefxn -> show ( std::cout, *pose );
	//pose->dump_pdb("test.pdb");
}

///////////////////////////////////////////////////////////////
void*
my_main ( void* ) {
	pdb_scoring();
	exit ( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main ( int argc, char * argv [] ) {
    try {
	using namespace basic::options;
	NEW_OPT ( force_field_file, "score_file", "rna/rna_hires_elec_dens" );
	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	core::init::init ( argc, argv );
	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////
	protocols::viewer::viewer_main ( my_main );
    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
                                  }
        return 0;
}
