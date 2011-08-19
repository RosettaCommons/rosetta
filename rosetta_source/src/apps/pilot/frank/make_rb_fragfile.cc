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
#include <devel/init.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/util.hh>
#include <core/fragment/FragmentIO.hh>

#include <core/kinematics/FoldTree.hh>

#include <protocols/abinitio/ClassicAbinitio.hh>
#include <protocols/relax_protocols.hh>
#include <protocols/abinitio/FoldConstraints.hh>
#include <protocols/idealize/idealize.hh>
#include <protocols/jumping/Dssp.hh>

#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/evaluation/RmsdEvaluator.hh>

#include <core/chemical/ChemicalManager.hh>



#include <core/conformation/ResidueFactory.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>


#include <basic/options/option.hh>
//#include <basic/options/OptionKeys.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>
#include <protocols/evaluation/ChemicalShiftEvaluator.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>
// C++ headers
#include <fstream>
#include <iostream>
#include <string>

#include <protocols/viewer/viewers.hh>
#include <basic/Tracer.hh>

#include <core/id/NamedStubID.hh>
#include <core/kinematics/Stub.hh>
#include <core/id/types.hh>
#include <numeric/xyz.functions.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>



int main(int argc, char **argv) {
	using namespace core::fragment;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;

	using core::Size;
	using std::string;

	devel::init( argc,argv );

	core::pose::Pose pose;
	std::vector< utility::file::FileName > pdb_file_names = option[ in::file::s ]().vector();

	for (int i=0; i<pdb_file_names.size(); ++i) {
		core::import_pose::pose_from_pdb( pose, *core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ), pdb_file_names[i] );

		protocols::jumping::Dssp dssp_obj( pose );
		dssp_obj.insert_ss_into_pose( pose );

		core::Size nres = pose.total_residue();

		for (core::Size i=1;i<=nres;++i) {
			//std::sscanf( line+6 , "%1c", &seq);
			char seq = pose.residue(i).name1();
			//std::sscanf( line+8 , "%1c", &ss);
			char ss = pose.secstruct(i);
			//std::sscanf( line+25 , "%9f", &x);
			//std::sscanf( line+34 , "%9f", &y);
			//std::sscanf( line+43 , "%9f", &z);
			numeric::xyzVector<core::Real> x=pose.residue(i).atom("CA").xyz();
			//std::sscanf( line+52, "%9f", &phi);
			//std::sscanf( line+61, "%9f", &psi);
			//std::sscanf( line+70, "%9f", &omega);
			core::Real phi = pose.residue(i).mainchain_torsion( 1 );
			core::Real psi = pose.residue(i).mainchain_torsion( 2 );
			core::Real omega = pose.residue(i).mainchain_torsion( 3 );
			std::printf( "xxxxX %1c %1c  %4i %4i %4i %8.2f %8.2f %8.2f %8.3f %8.3f %8.3f \n", seq,ss,(int)i,(int)(nres-i),(int)i,x[0],x[1],x[2],phi,psi,omega);
		}
	}

	return 0;
}
