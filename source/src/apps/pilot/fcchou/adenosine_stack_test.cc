// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


// libRosetta headers
#include <core/scoring/rms_util.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/chemical/rna/util.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <utility/io/ozstream.hh>

#include <protocols/idealize/idealize.hh>
#include <protocols/stepwise/StepWiseUtil.cc>

#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <core/init/init.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/EnergyMap.hh> //for EnergyMap
#include <core/scoring/EnergyMap.fwd.hh> //for EnergyMap
#include <core/import_pose/import_pose.hh>

#include <protocols/viewer/viewers.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include <utility/excn/Exceptions.hh>


using namespace core;
using namespace protocols;
using namespace basic::options::OptionKeys;
using namespace basic::options;

using utility::vector1;

typedef  numeric::xyzMatrix< Real > Matrix;

OPT_KEY( Boolean, sample_water )
OPT_KEY( Real, alpha_increment )
OPT_KEY( Real, cosbeta_increment )
OPT_KEY( Real, gamma_increment )

/////////////////////////////////////////////////////////////////////////////
//FCC: Adding Virtual res
void
add_virtual_res ( core::pose::Pose & pose, bool set_res_as_root = true ) {
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
/////////////////////////////////////////////////////////////////////////////////////////////
// Rhiju -- rotate to my favorite frame. Base centroid is now at origin.
//         X points to N1 atom. Z points normal to base. Y is orthonormal and points towards Hoogsteen edge, I think.
void
rotate_into_nucleobase_frame( core::pose::Pose & pose ){

	using namespace core::conformation;
	using namespace core::chemical::rna;
	using namespace core::id;

	// assuming pose has an RNA at residue 1 -- will rotate just that residue.
	Size const base_pos( 1 );
	Residue const & rsd = pose.residue( base_pos );

	Vector centroid = get_rna_base_centroid( rsd, true /*verbose*/ );
	Matrix M = get_rna_base_coordinate_system( rsd, centroid );
	kinematics::Stub stub( M, centroid );

	for (Size i = 1; i <= pose.total_residue(); ++i) {
		Residue const & res = pose.residue(i);
		for (Size j = 1; j <= res.natoms(); j++ ){
			Vector xyz_new = stub.global2local( res.xyz( j ) ); // it is either this or M-inverse.
			pose.set_xyz( AtomID( j, i ), xyz_new );
		}
	}
}

/////////////////////////////////////////////////
void
methane_pair_score_test()
{
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;

	//////////////////////////////////////////////////
	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	// Read in pose with two methane. "Z" = ligand. Note need flag:
	//         -extra_res_fa CH4.params -s two_methane.pdb
	pose::Pose pose;
	std::string infile  = option[ in ::file::s ][1];
	import_pose::pose_from_file( pose, *rsd_set, infile , core::import_pose::PDB_file);

	rotate_into_nucleobase_frame( pose );
	core::chemical::ResidueTypeSet const & residue_set = pose.residue_type ( 1 ).residue_type_set();

	core::conformation::ResidueOP new_res;

	pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", 1 );
	pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_RIBOSE", 1 );
	pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", 2 );
	pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_RIBOSE", 2 );

	pose.dump_pdb( "START.pdb" );
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 800, 800 );


	//////////////////////////////////////////////////////////////////
	// OK, how about a score function?
	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( "farna/rna_hires" );
	//scorefxn->set_weight( fa_elec, 1.0 );

	(*scorefxn)( pose );
	scorefxn->show( std::cout, pose );
}

///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	methane_pair_score_test();

	protocols::viewer::clear_conformation_viewers();

	exit( 0 );

}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
    try {

	NEW_OPT( sample_water, "use a water probe instead of carbon", false );
	NEW_OPT( alpha_increment, "input parameter", 40.0 );
	NEW_OPT( cosbeta_increment, "input parameter", 0.25 );
	NEW_OPT( gamma_increment, "input parameter", 40.0 );

	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	core::init::init(argc, argv);

  protocols::viewer::viewer_main( my_main );
    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }
        return 0;

}
