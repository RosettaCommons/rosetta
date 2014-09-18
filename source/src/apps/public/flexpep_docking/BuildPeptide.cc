// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:f;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2008 University of Washington
// (C) 199x-2008 University of California Santa Cruz
// (C) 199x-2008 University of California San Francisco
// (C) 199x-2008 Johns Hopkins University
// (C) 199x-2008 University of North Carolina, Chapel Hill
// (C) 199x-2008 Vanderbilt University
// (C) 199x-2008 Hebrew University, Jerusalem
//
/// @file   BuildPeptide.cc
//
/// @brief Application that reads in a peptides sequence file and outputs a linear peptide.
/// @author Barak Raveh, Nir London
/// @date June 01, 2009

//#define GL_GRAPHICS

// Project headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>
// AUTO-REMOVED #include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.hh>
//#include <core/chemical/ResidueSelector.hh>
// AUTO-REMOVED #include <core/chemical/util.hh>
//#include <core/scoring/ScoreFunction.hh>
//#include <core/scoring/ScoreFunctionFactory.hh>
//#include <core/pack/pack_rotamers.hh>
//#include <core/pack/task/PackerTask.hh>
//#include <core/pack/task/TaskFactory.hh>
#include <core/pose/annotated_sequence.hh>
//#include <core/kinematics/MoveMap.hh>
//#include <core/optimization/AtomTreeMinimizer.hh>
//#include <core/optimization/MinimizerOptions.hh>
#include <core/pose/Pose.hh>
//#include <core/options/util.hh>
#include <devel/init.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
//#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
//#include <core/conformation/ResidueFactory.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>

#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers

using basic::T;
using basic::Error;
using basic::Warning;


static thread_local basic::Tracer TR( "BuildPeptide" );

using namespace core;
using namespace basic::options;
using namespace OptionKeys;


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
	using namespace pose;
	using namespace scoring;
	using namespace conformation;
	using namespace core::chemical;

	//setup random numbers and options
	devel::init(argc, argv);

   //create a pose
   pose::Pose pose;

	//protein pose
	//pose::Pose prot_pose;
	//io::pdb::pose_from_pdb( prot_pose, options::start_file() ); // gets filename from -s option

  //read peptides fasta file
  std::string pepSeq = core::sequence::read_fasta_file( basic::options::option[ in::file::fasta ]()[1] )[1]->sequence();
	int seqLen = pepSeq.length();

  make_pose_from_sequence(pose,pepSeq, *ChemicalManager::get_instance()->residue_type_set( "fa_standard" ));

   //make peptide extended
   for (int i=1; i<=seqLen; i++) {
            pose.set_phi(i,-135.0);
            pose.set_psi(i,135.0);
            pose.set_omega(i,180.0);
   }

   //dump pdb to output
	 std::string oFileName= basic::options::option[ out::file::o ]();
   pose.dump_pdb(oFileName);

   } catch ( utility::excn::EXCN_Base const & e ) {
	  std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
   }

   return 0;
}
