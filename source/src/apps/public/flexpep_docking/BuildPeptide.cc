// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
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
#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>
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
//#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
//#include <core/conformation/ResidueFactory.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>

#include <basic/options/option_macros.hh>

#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers

using basic::T;
using basic::Error;
using basic::Warning;


static THREAD_LOCAL basic::Tracer TR( "BuildPeptide" );

using namespace core;
using namespace basic::options;
using namespace OptionKeys;

OPT_KEY( Boolean, helix )
OPT_KEY( Real,  phi )
OPT_KEY( Real,  psi )

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		using namespace pose;
		using namespace scoring;
		using namespace conformation;
		using namespace core::chemical;

		NEW_OPT( helix, "Make helical peptide", false );
		NEW_OPT( phi, "Repeated value for phi", -135.0 );
		NEW_OPT( psi, "Repeated value for psi", 135.0 );

		//setup random numbers and options
		devel::init(argc, argv);

		//create a pose
		pose::Pose pose;

		//protein pose
		//pose::Pose prot_pose;
		//io::pdb::pose_from_file( prot_pose, options::start_file() , core::import_pose::PDB_file); // gets filename from -s option

		//read peptides fasta file
		std::string pepSeq = core::sequence::read_fasta_file( basic::options::option[ in::file::fasta ]()[1] )[1]->sequence();
		int seqLen = pepSeq.length();

		make_pose_from_sequence(pose,pepSeq, *ChemicalManager::get_instance()->residue_type_set( "fa_standard" ));

		//make peptide extended
		for ( int i=1; i<=seqLen; i++ ) {
			if ( basic::options::option[ phi ].user() || basic::options::option[ psi ].user() ) {
				pose.set_phi(i, basic::options::option[ phi ]() );
				pose.set_psi(i, basic::options::option[ psi ]() );
				pose.set_omega(i,180.0);
			} else if ( basic::options::option[ helix ]() ) {
				pose.set_phi(i,-57.8); // geometrically ideal values from http://www.cryst.bbk.ac.uk/PPS2/course/section8/ss-960531_6.html
				pose.set_psi(i,-47.0);
				pose.set_omega(i,180.0);
			} else {
				pose.set_phi(i,-135.0);
				pose.set_psi(i,135.0);
				pose.set_omega(i,180.0);
			}
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
