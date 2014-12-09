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
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomType.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <basic/database/open.hh>
#include <core/init/init.hh>
#include <core/io/pdb/pose_io.hh>

#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.fwd.hh>
//////////////////////////////////////////////////
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
///////////////////////////////////////////////////
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/util.hh>
#include <basic/options/option_macros.hh>
#include <protocols/idealize/idealize.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <protocols/viewer/viewers.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh> //for EnergyMap
#include <core/scoring/EnergyMap.fwd.hh> //for EnergyMap
#include <basic/basic.hh>


#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/conversions.hh>


#include <string>
#include <map>

//////////////////////////////////////////////////////////

#include <protocols/farna/util.hh>
#include <protocols/farna/RNA_BasePairClassifier.hh>


#include <core/scoring/rms_util.tmpl.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>



#include <protocols/farna/RNA_LoopCloser.hh>
#include <protocols/farna/RNA_LoopCloser.fwd.hh>

#include <core/scoring/rna/RNA_BaseDoubletClasses.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/BinarySilentStruct.hh>
// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>  //Test
#include <cctype>
#include <iomanip>
#include <map>
#include <cstdlib>
#include <ctime>

//Added by Parin
//#include <core/scoring/ScoreType.hh>
#include <list>
#include <stdio.h>
#include <math.h>

//silly using/typedef
#include <basic/Tracer.hh>

//#include <core/util/tracer.hh>
using basic::T;
using basic::Error;
using basic::Warning;
using namespace core;
using namespace protocols;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using utility::vector1;
using io::pdb::dump_pdb;

typedef  numeric::xyzMatrix< Real > Matrix;


OPT_KEY( String, 	output_silent_file)


//~/src/mini/bin/convert_pdb_to_silent_file.macosgccrelease  -s 1zih_RNA.pdb -output_silent_file output_test.out -database ~/minirosetta_database/

std::string
get_tag_from_pdb_filename(std::string const pdb_filename){

	std::string tag;

	size_t found=pdb_filename.rfind('/');

	if(found!=std::string::npos){
		tag = pdb_filename.substr(found+1);
	} else {
		tag=pdb_filename;
	}

	size_t found_2=tag.rfind(".pdb");

	if(found_2!=std::string::npos){
		tag = tag.substr(0,tag.size()-4);
	}

	return tag;
}

void
pdb_to_silent_file_simple(){

	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace core::conformation;
	using namespace ObjexxFCL;
	using namespace core::io::silent;
	using namespace core::scoring;
	using namespace core::pose;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );
	SilentFileData silent_file_data;

	pose::Pose viewer_pose;

	if ( !option[ in::file::s ].user() ) utility_exit_with_message( "User must supply in::file::s!" );

	if ( !option[ output_silent_file ].user() ) utility_exit_with_message( "User must supply output_silent_file!" );

	std::string const pdb_file= option[ in::file::s ]()[1];

	std::string const silent_outfile=option[ output_silent_file]();

	std::string tag=get_tag_from_pdb_filename(pdb_file);

	std::cout << "importing pdb_file: " << pdb_file << std::endl;

	pose::Pose pose;

	core::import_pose::pose_from_pdb( pose, *rsd_set, pdb_file );

	// NEW! from rhiju -- pay attention to chain breaks.
	protocols::farna::figure_out_reasonable_rna_fold_tree( pose );

	if(pose.residue(1).atom(1).xyz().length() < 2.0){
	    numeric::xyzMatrix< Real > R( 0.0 );
	    R.to_identity();
	    Vector v( 5.0 );
	    pose.apply_transform_Rx_plus_v( R, v );
	}

	BinarySilentStruct s( pose, tag );

	silent_file_data.write_silent_struct(s, silent_outfile, false);

}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

  using namespace basic::options;

	pdb_to_silent_file_simple();
  exit( 0 );

}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try{
  using namespace basic::options;

	NEW_OPT( output_silent_file, "output_silent_file", "");


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



