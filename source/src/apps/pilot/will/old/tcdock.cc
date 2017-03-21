// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/////////// WISH LIST ///////{
	//		buried unsats
	//			scan for rotamers that cna make H-bonds
	//			detect and penalize missing BB density
	//		strand pairs
	//		centriod score components ?
	//		disulfide-compatible positions
	//		low-res sc? maybe patchdock-like complimentarity measure?
	//		statistically "good" BB xforms? all of by good contacts?
	//	IGNORE
	//		contacts by SS ?? avg. degree probably sufficient
	//	DONE
	//		termini distance
	//		contacts weight by avg. deg.

////////////////// INCLUDS /////////////////////////////////////////////{
	#include <protocols/sic_dock/SICFast.hh>
	#include <protocols/sic_dock/RigidScore.hh>
	#include <protocols/sic_dock/scores/MotifHashRigidScore.hh>
	#include <protocols/sic_dock/util.hh>
	#include <protocols/sic_dock/types.hh>
	#include <basic/options/keys/in.OptionKeys.gen.hh>
	#include <basic/options/keys/out.OptionKeys.gen.hh>
	#include <basic/options/keys/mh.OptionKeys.gen.hh>
	#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
	#include <basic/options/keys/sicdock.OptionKeys.gen.hh>
	#include <basic/options/option.hh>
	#include <basic/options/option_macros.hh>
	#include <basic/Tracer.hh>
	#include <core/chemical/ChemicalManager.hh>
	#include <core/chemical/AtomType.hh>
	#include <core/chemical/ResidueTypeSet.hh>
	#include <core/conformation/Residue.hh>
	#include <core/conformation/symmetry/SymmetricConformation.hh>
	// #include <core/conformation/symmetry/SymDof.hh>
	#include <core/id/AtomID.hh>
	#include <core/import_pose/import_pose.hh>
	#include <core/io/pdb/pdb_writer.hh>
	#include <core/pose/Pose.hh>
	#include <core/pose/symmetry/util.hh>
	#include <core/pose/util.hh>
	#include <core/scoring/dssp/Dssp.hh>
	#include <core/scoring/sasa.hh>
	#include <core/util/SwitchResidueTypeSet.hh>
	#include <devel/init.hh>
	#include <numeric/constants.hh>
	#include <numeric/xyz.functions.hh>
	#include <numeric/xyz.io.hh>
	#include <ObjexxFCL/FArray2D.hh>
	#include <ObjexxFCL/FArray3D.hh>
	#include <ObjexxFCL/format.hh>
	#include <ObjexxFCL/string.functions.hh>
	#include <utility/io/ozstream.hh>
	#include <utility/io/izstream.hh>
	#include <utility/string_util.hh>
	#include <utility/file/file_sys_util.hh>
	#include <utility/vector1.hh>
	#include <core/scoring/methods/RG_Energy_Fast.hh>
	#include <numeric/xyzTransform.hh>
	#include <basic/database/open.hh>
	#include <protocols/sic_dock/xyzStripeHashPose.hh>

	#ifdef USE_OPENMP
	#include <omp.h>
	#endif

////////////////////////// TYPEDEFS ////////////////////////////////////{
	using core::Size;
	using core::Real;
	using core::pose::Pose;
	using std::string;
	using utility::vector1;
	using ObjexxFCL::format::I;
	using ObjexxFCL::format::F;
	using ObjexxFCL::format::RJ;
	using numeric::min;
	using numeric::max;
	using std::cout;
	using std::cerr;
	using std::endl;
	using namespace protocols::sic_dock;


static numeric::Xforms SYMOCT;
static numeric::Xforms SYMTET;
static numeric::Xforms SYMICS;
void initsyms(){
	SYMTET.push_back( Xform( Mat::rows( Vec( 1.000000, 0.000000, 0.000000 ), Vec( 0.000000, -1.000000, 0.000000 ), Vec( 0.000000, -0.000000, -1.000000 ) ), Vec( 0.000000, 0.000000, 0.000000 ) ) );
	SYMTET.push_back( Xform( Mat::rows( Vec( 0.000000, 1.000000, -0.000000 ), Vec( -0.000000, 0.000000, 1.000000 ), Vec( 1.000000, -0.000000, 0.000000 ) ), Vec( 0.000000, 0.000000, 0.000000 ) ) );
	SYMTET.push_back( Xform( Mat::rows( Vec( 1.000000, 0.000000, 0.000000 ), Vec( 0.000000, 1.000000, -0.000000 ), Vec( 0.000000, 0.000000, 1.000000 ) ), Vec( 0.000000, 0.000000, 0.000000 ) ) );
	SYMTET.push_back( Xform( Mat::rows( Vec( 0.000000, -1.000000, 0.000000 ), Vec( -0.000000, -0.000000, -1.000000 ), Vec( 1.000000, 0.000000, -0.000000 ) ), Vec( 0.000000, 0.000000, 0.000000 ) ) );
	SYMTET.push_back( Xform( Mat::rows( Vec( 0.000000, 1.000000, -0.000000 ), Vec( 0.000000, -0.000000, -1.000000 ), Vec( -1.000000, 0.000000, -0.000000 ) ), Vec( 0.000000, 0.000000, 0.000000 ) ) );
	SYMTET.push_back( Xform( Mat::rows( Vec( -0.000000, 0.000000, 1.000000 ), Vec( 1.000000, -0.000000, 0.000000 ), Vec( 0.000000, 1.000000, -0.000000 ) ), Vec( 0.000000, 0.000000, 0.000000 ) ) );
	SYMTET.push_back( Xform( Mat::rows( Vec( 0.000000, -1.000000, 0.000000 ), Vec( 0.000000, 0.000000, 1.000000 ), Vec( -1.000000, -0.000000, 0.000000 ) ), Vec( 0.000000, 0.000000, 0.000000 ) ) );
	SYMTET.push_back( Xform( Mat::rows( Vec( -0.000000, -0.000000, -1.000000 ), Vec( 1.000000, 0.000000, -0.000000 ), Vec( 0.000000, -1.000000, 0.000000 ) ), Vec( 0.000000, 0.000000, 0.000000 ) ) );
	SYMTET.push_back( Xform( Mat::rows( Vec( 0.000000, -0.000000, -1.000000 ), Vec( -1.000000, -0.000000, -0.000000 ), Vec( -0.000000, 1.000000, -0.000000 ) ), Vec( 0.000000, 0.000000, 0.000000 ) ) );
	SYMTET.push_back( Xform( Mat::rows( Vec( -0.000000, 0.000000, 1.000000 ), Vec( -1.000000, 0.000000, -0.000000 ), Vec( -0.000000, -1.000000, 0.000000 ) ), Vec( 0.000000, 0.000000, 0.000000 ) ) );
	SYMTET.push_back( Xform( Mat::rows( Vec( -1.000000, -0.000000, -0.000000 ), Vec( -0.000000, 1.000000, 0.000000 ), Vec( 0.000000, 0.000000, -1.000000 ) ), Vec( 0.000000, 0.000000, 0.000000 ) ) );
	SYMTET.push_back( Xform( Mat::rows( Vec( -1.000000, 0.000000, -0.000000 ), Vec( -0.000000, -1.000000, 0.000000 ), Vec( -0.000000, 0.000000, 1.000000 ) ), Vec( 0.000000, 0.000000, 0.000000 ) ) );
	SYMOCT.push_back( Xform( Mat::rows( Vec(0.000000,1.000000,0.000000), Vec(1.000000,0.000000,-0.000000), Vec(-0.000000,0.000000,-1.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMOCT.push_back( Xform( Mat::rows( Vec(0.000000,-0.000000,1.000000), Vec(1.000000,0.000000,-0.000000), Vec(-0.000000,1.000000,0.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMOCT.push_back( Xform( Mat::rows( Vec(1.000000,0.000000,-0.000000), Vec(0.000000,1.000000,0.000000), Vec(0.000000,-0.000000,1.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMOCT.push_back( Xform( Mat::rows( Vec(1.000000,0.000000,-0.000000), Vec(0.000000,-0.000000,1.000000), Vec(0.000000,-1.000000,-0.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMOCT.push_back( Xform( Mat::rows( Vec(-0.000000,0.000000,-1.000000), Vec(0.000000,1.000000,0.000000), Vec(1.000000,-0.000000,-0.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMOCT.push_back( Xform( Mat::rows( Vec(-0.000000,1.000000,0.000000), Vec(0.000000,-0.000000,1.000000), Vec(1.000000,0.000000,-0.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMOCT.push_back( Xform( Mat::rows( Vec(0.000000,1.000000,0.000000), Vec(-0.000000,0.000000,-1.000000), Vec(-1.000000,0.000000,0.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMOCT.push_back( Xform( Mat::rows( Vec(0.000000,-0.000000,1.000000), Vec(-0.000000,1.000000,0.000000), Vec(-1.000000,-0.000000,0.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMOCT.push_back( Xform( Mat::rows( Vec(0.000000,-1.000000,-0.000000), Vec(1.000000,0.000000,0.000000), Vec(0.000000,-0.000000,1.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMOCT.push_back( Xform( Mat::rows( Vec(1.000000,-0.000000,-0.000000), Vec(-0.000000,0.000000,-1.000000), Vec(0.000000,1.000000,0.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMOCT.push_back( Xform( Mat::rows( Vec(1.000000,0.000000,0.000000), Vec(0.000000,-1.000000,-0.000000), Vec(-0.000000,0.000000,-1.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMOCT.push_back( Xform( Mat::rows( Vec(-0.000000,0.000000,-1.000000), Vec(1.000000,-0.000000,-0.000000), Vec(-0.000000,-1.000000,-0.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMOCT.push_back( Xform( Mat::rows( Vec(-1.000000,0.000000,0.000000), Vec(0.000000,1.000000,-0.000000), Vec(-0.000000,0.000000,-1.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMOCT.push_back( Xform( Mat::rows( Vec(-1.000000,-0.000000,0.000000), Vec(0.000000,0.000000,1.000000), Vec(-0.000000,1.000000,0.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMOCT.push_back( Xform( Mat::rows( Vec(0.000000,-0.000000,1.000000), Vec(0.000000,-1.000000,-0.000000), Vec(1.000000,0.000000,0.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMOCT.push_back( Xform( Mat::rows( Vec(0.000000,1.000000,-0.000000), Vec(-1.000000,0.000000,0.000000), Vec(0.000000,-0.000000,1.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMOCT.push_back( Xform( Mat::rows( Vec(0.000000,0.000000,1.000000), Vec(-1.000000,-0.000000,0.000000), Vec(0.000000,-1.000000,-0.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMOCT.push_back( Xform( Mat::rows( Vec(0.000000,-1.000000,-0.000000), Vec(-0.000000,-0.000000,1.000000), Vec(-1.000000,-0.000000,-0.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMOCT.push_back( Xform( Mat::rows( Vec(-0.000000,-1.000000,-0.000000), Vec(-0.000000,0.000000,-1.000000), Vec(1.000000,-0.000000,-0.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMOCT.push_back( Xform( Mat::rows( Vec(-0.000000,0.000000,-1.000000), Vec(-1.000000,0.000000,0.000000), Vec(0.000000,1.000000,-0.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMOCT.push_back( Xform( Mat::rows( Vec(0.000000,0.000000,-1.000000), Vec(-0.000000,-1.000000,-0.000000), Vec(-1.000000,0.000000,0.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMOCT.push_back( Xform( Mat::rows( Vec(-1.000000,0.000000,0.000000), Vec(-0.000000,0.000000,-1.000000), Vec(-0.000000,-1.000000,0.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMOCT.push_back( Xform( Mat::rows( Vec(-1.000000,-0.000000,-0.000000), Vec(0.000000,-1.000000,-0.000000), Vec(-0.000000,-0.000000,1.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMOCT.push_back( Xform( Mat::rows( Vec(0.000000,-1.000000,-0.000000), Vec(-1.000000,-0.000000,-0.000000), Vec(0.000000,0.000000,-1.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(1.000000,0.000000,0.000000), Vec(0.000000,-1.000000,-0.000000), Vec(0.000000,0.000000,-1.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.809017,-0.309017,0.500000), Vec(0.309017,-0.500000,-0.809017), Vec(0.500000,0.809017,-0.309017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(1.000000,0.000000,0.000000), Vec(0.000000,1.000000,0.000000), Vec(0.000000,-0.000000,1.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.809017,-0.309017,0.500000), Vec(-0.309017,0.500000,0.809017), Vec(-0.500000,-0.809017,0.309017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.809017,0.309017,-0.500000), Vec(0.309017,0.500000,0.809017), Vec(0.500000,-0.809017,0.309017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.809017,0.309017,0.500000), Vec(-0.309017,-0.500000,0.809017), Vec(0.500000,-0.809017,-0.309017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.809017,0.309017,-0.500000), Vec(-0.309017,-0.500000,-0.809017), Vec(-0.500000,0.809017,-0.309017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.809017,0.309017,0.500000), Vec(0.309017,0.500000,-0.809017), Vec(-0.500000,0.809017,0.309017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.500000,-0.809017,0.309017), Vec(0.809017,0.309017,-0.500000), Vec(0.309017,0.500000,0.809017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.809017,-0.309017,-0.500000), Vec(-0.309017,0.500000,-0.809017), Vec(0.500000,0.809017,0.309017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.500000,-0.809017,0.309017), Vec(-0.809017,-0.309017,0.500000), Vec(-0.309017,-0.500000,-0.809017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.809017,-0.309017,-0.500000), Vec(0.309017,-0.500000,0.809017), Vec(-0.500000,-0.809017,-0.309017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.500000,0.809017,-0.309017), Vec(0.809017,-0.309017,0.500000), Vec(0.309017,-0.500000,-0.809017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.309017,0.500000,0.809017), Vec(0.500000,-0.809017,0.309017), Vec(0.809017,0.309017,-0.500000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.309017,-0.500000,0.809017), Vec(-0.500000,-0.809017,-0.309017), Vec(0.809017,-0.309017,-0.500000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.500000,0.809017,-0.309017), Vec(-0.809017,0.309017,-0.500000), Vec(-0.309017,0.500000,0.809017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.309017,0.500000,0.809017), Vec(-0.500000,0.809017,-0.309017), Vec(-0.809017,-0.309017,0.500000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.309017,-0.500000,0.809017), Vec(0.500000,0.809017,0.309017), Vec(-0.809017,0.309017,0.500000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.500000,-0.809017,-0.309017), Vec(0.809017,0.309017,0.500000), Vec(-0.309017,-0.500000,0.809017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.309017,-0.500000,-0.809017), Vec(0.500000,0.809017,-0.309017), Vec(0.809017,-0.309017,0.500000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.309017,0.500000,-0.809017), Vec(-0.500000,0.809017,0.309017), Vec(0.809017,0.309017,0.500000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.500000,0.809017,0.309017), Vec(-0.809017,0.309017,0.500000), Vec(0.309017,-0.500000,0.809017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.500000,-0.809017,-0.309017), Vec(-0.809017,-0.309017,-0.500000), Vec(0.309017,0.500000,-0.809017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.309017,-0.500000,-0.809017), Vec(-0.500000,-0.809017,0.309017), Vec(-0.809017,0.309017,-0.500000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.309017,0.500000,-0.809017), Vec(0.500000,-0.809017,-0.309017), Vec(-0.809017,-0.309017,-0.500000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.500000,0.809017,0.309017), Vec(0.809017,-0.309017,-0.500000), Vec(-0.309017,0.500000,-0.809017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.000000,0.000000,1.000000), Vec(1.000000,-0.000000,-0.000000), Vec(0.000000,1.000000,-0.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(-0.309017,-0.500000,0.809017), Vec(0.500000,-0.809017,-0.309017), Vec(0.809017,0.309017,0.500000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.000000,-1.000000,0.000000), Vec(-0.000000,-0.000000,-1.000000), Vec(1.000000,0.000000,-0.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.000000,0.000000,1.000000), Vec(-1.000000,0.000000,0.000000), Vec(-0.000000,-1.000000,0.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(-0.309017,-0.500000,0.809017), Vec(-0.500000,0.809017,0.309017), Vec(-0.809017,-0.309017,-0.500000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.000000,-1.000000,0.000000), Vec(0.000000,0.000000,1.000000), Vec(-1.000000,-0.000000,0.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.000000,-0.000000,-1.000000), Vec(1.000000,0.000000,0.000000), Vec(0.000000,-1.000000,0.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(-0.309017,0.500000,-0.809017), Vec(0.500000,0.809017,0.309017), Vec(0.809017,-0.309017,-0.500000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.000000,1.000000,0.000000), Vec(0.000000,-0.000000,1.000000), Vec(1.000000,-0.000000,-0.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(-0.309017,0.500000,0.809017), Vec(-0.500000,-0.809017,0.309017), Vec(0.809017,-0.309017,0.500000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.000000,-0.000000,-1.000000), Vec(-1.000000,-0.000000,-0.000000), Vec(-0.000000,1.000000,-0.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(-0.309017,0.500000,-0.809017), Vec(-0.500000,-0.809017,-0.309017), Vec(-0.809017,0.309017,0.500000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(0.000000,1.000000,0.000000), Vec(-0.000000,0.000000,-1.000000), Vec(-1.000000,0.000000,0.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(-0.309017,0.500000,0.809017), Vec(0.500000,0.809017,-0.309017), Vec(-0.809017,0.309017,-0.500000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(-0.500000,-0.809017,0.309017), Vec(0.809017,-0.309017,0.500000), Vec(-0.309017,0.500000,0.809017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(-0.500000,-0.809017,-0.309017), Vec(0.809017,-0.309017,-0.500000), Vec(0.309017,-0.500000,0.809017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(-0.309017,-0.500000,-0.809017), Vec(-0.500000,0.809017,-0.309017), Vec(0.809017,0.309017,-0.500000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(-0.500000,-0.809017,0.309017), Vec(-0.809017,0.309017,-0.500000), Vec(0.309017,-0.500000,-0.809017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(-0.500000,-0.809017,-0.309017), Vec(-0.809017,0.309017,0.500000), Vec(-0.309017,0.500000,-0.809017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(-0.309017,-0.500000,-0.809017), Vec(0.500000,-0.809017,0.309017), Vec(-0.809017,-0.309017,0.500000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(-0.500000,0.809017,-0.309017), Vec(0.809017,0.309017,-0.500000), Vec(-0.309017,-0.500000,-0.809017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(-0.500000,0.809017,0.309017), Vec(0.809017,0.309017,0.500000), Vec(0.309017,0.500000,-0.809017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(-0.809017,0.309017,0.500000), Vec(0.309017,-0.500000,0.809017), Vec(0.500000,0.809017,0.309017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(-0.809017,-0.309017,0.500000), Vec(-0.309017,-0.500000,-0.809017), Vec(0.500000,-0.809017,0.309017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(-0.500000,0.809017,-0.309017), Vec(-0.809017,-0.309017,0.500000), Vec(0.309017,0.500000,0.809017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(-0.500000,0.809017,0.309017), Vec(-0.809017,-0.309017,-0.500000), Vec(-0.309017,-0.500000,0.809017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(-0.809017,0.309017,0.500000), Vec(-0.309017,0.500000,-0.809017), Vec(-0.500000,-0.809017,-0.309017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(-0.809017,-0.309017,0.500000), Vec(0.309017,0.500000,0.809017), Vec(-0.500000,0.809017,-0.309017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(-0.809017,-0.309017,-0.500000), Vec(0.309017,0.500000,-0.809017), Vec(0.500000,-0.809017,-0.309017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(-0.809017,0.309017,-0.500000), Vec(-0.309017,0.500000,0.809017), Vec(0.500000,0.809017,-0.309017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(-0.809017,-0.309017,-0.500000), Vec(-0.309017,-0.500000,0.809017), Vec(-0.500000,0.809017,0.309017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(-0.809017,0.309017,-0.500000), Vec(0.309017,-0.500000,-0.809017), Vec(-0.500000,-0.809017,0.309017) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(-1.000000,-0.000000,0.000000), Vec(0.000000,-1.000000,-0.000000), Vec(0.000000,0.000000,1.000000) ), Vec(0.000000,0.000000,0.000000) ) );
	SYMICS.push_back( Xform( Mat::rows( Vec(-1.000000,-0.000000,0.000000), Vec(-0.000000,1.000000,0.000000), Vec(-0.000000,-0.000000,-1.000000) ), Vec(0.000000,0.000000,0.000000) ) );
}

///////////////////////////// ARCH /////////////////////////////////////{
	class Arch : public utility::pointer::ReferenceCount {
	public:
		string name_;
		int    nfold1_;
		int    nfold2_;
		Vec    axis1_;
		Vec    axis2_;
		Vec    move1_;
		Vec    move2_;
		Xform  intra_xform1_;
		Xform  intra_xform2_;
		Vec    intra_slide1_;
		Vec    intra_slide2_;
		Vec    intra_move1A_;
		Vec    intra_move1B_;
		Vec    intra_move2A_;
		Vec    intra_move2B_;
		Real   alpha_;
		Real   sin_alpha_;
		Real   tan_alpha_;
		int    nangle1_;
		int    nangle2_;
		Real   angle1_;
		Real   angle2_;
	public:
		string name()         const { return name_; }
		int    nfold1()       const { return nfold1_; }
		int    nfold2()       const { return nfold2_; }
		Vec    axis1()        const { return axis1_; }
		Vec    axis2()        const { return axis2_; }
		Vec    move1()        const { return move1_; }
		Vec    move2()        const { return move2_; }
		Xform  intra_xform1() const { return intra_xform1_; }
		Xform  intra_xform2() const { return intra_xform2_; }
		Vec    intra_slide1() const { return intra_slide1_; }
		Vec    intra_slide2() const { return intra_slide2_; }
		Vec    intra_move1A() const { return intra_move1A_; }
		Vec    intra_move1B() const { return intra_move1B_; }
		Vec    intra_move2A() const { return intra_move2A_; }
		Vec    intra_move2B() const { return intra_move2B_; }
		int    nangle1()      const { return nangle1_; }
		int    nangle2()      const { return nangle2_; }
		Real   angle1()       const { return angle1_; }
		Real   angle2()       const { return angle2_; }
		virtual int max_ori() const { return 360; }
		virtual bool dihedral () const { return false; }
		virtual bool dihedral1() const { return false; }
		virtual bool dihedral2() const { return false; }
		Xform xform1( Real const & a1, Real const & r1 ) const { return Xform(rotation_matrix_degrees(axis1_,a1),r1*move1_); }
		Xform xform2( Real const & a2, Real const & r2 ) const { return Xform(rotation_matrix_degrees(axis2_,a2),r2*move2_); }

		virtual Xform generator1(Real /*a1*/, Real /*a2*/, Real /*d1*/, Real /*d2*/) const {
			return Xform::rot_deg(axis1_,angle1_,Vec(0,0,0));
		}
		virtual Xform generator2(Real /*a1*/, Real /*a2*/, Real /*d1*/, Real /*d2*/) const {
			return Xform::rot_deg(axis2_,angle2_,Vec(0,0,0));
		}

		// virtual Real disp_dof1(Real const & omega, Real const & slidedis) const {
		// 	Real const gamma = numeric::conversions::radians(omega-alpha_/2.0);
		// 	Real const sin_gamma=sin(gamma);
		// 	// Real const cos_gamma=cos(gamma);
		// 	Real const x = slidedis*sin_gamma;
		// 	// Real const y = slidedis*cos_gamma;
		// 	Real const w = x/sin_alpha_;
		// 	// Real const z = x/tan_alpha_;
		// 	// dcmp1 = w;
		// 	// dcmp2 = y+z;
		// 	return w;
		// }
		// virtual Real disp_dof2(Real const & omega, Real const & slidedis) const {
		// 	Real const gamma = numeric::conversions::radians(omega-alpha_/2.0);
		// 	Real const sin_gamma=sin(gamma);
		// 	Real const cos_gamma=cos(gamma);
		// 	Real const x = slidedis*sin_gamma;
		// 	Real const y = slidedis*cos_gamma;
		// 	// Real const w = x/sin_alpha_;
		// 	Real const z = x/tan_alpha_;
		// 	// dcmp1 = w;
		// 	// dcmp2 = y+z;
		// 	return y+z;
		// }
		virtual void get_disp_dofs(Real const & omega, Real const & slidedis, Real & dcmp1_out, Real & dcmp2_out) const {
			Real const gamma = numeric::conversions::radians(omega-alpha_/2.0);
			Real sin_gamma=sin(gamma);
			Real cos_gamma=cos(gamma);
			Real const x = slidedis*sin_gamma;
			Real const y = slidedis*cos_gamma;
			Real const w = x/sin_alpha_;
			Real const z = x/tan_alpha_;
			dcmp1_out = w;   //disp_dof1(omega,slidedis);
			dcmp2_out = y+z; //disp_dof2(omega,slidedis);
			if( dcmp1_out * dcmp2_out < 0.0 ) utility_exit_with_message("dcmp1 & dcmp2 have opposite signs!");
		}
		virtual bool ori_is_allowed(Real const omega) const {
			if(omega <= alpha_ || omega >= 360.0 || (180.0 <= omega && omega <= 180.0+alpha_) ) return false;
			return true;

		}

		virtual Real intra_check_dof1(Real const &   d1  , Real const & /*d2*/) const { return d1; }
		virtual Real intra_check_dof2(Real const & /*d1*/, Real const &   d2  ) const { return d2; }
		virtual Real intra_slidedis_to_component_disp(int /*slidecomp*/, int dispcomp, Real d) const {
			if(1==dispcomp){
				return d/2.0/sin( angle_radians(axis1_,Vec(0,0,0),intra_xform1_*axis1_)/2.0 );
			} else {
				return d/2.0/sin( angle_radians(axis2_,Vec(0,0,0),intra_xform2_*axis2_)/2.0 );
			}
			utility_exit_with_message("bad input");
			return 0.0;
		}

		virtual Vec slide_axis(Real omega) const {
			return rotation_matrix_degrees( move1_.cross(move2_),omega)*(axis1_+axis2_).normalized();
		}

		void init( string _name, int _nfold1, Vec _axis1, int _nfold2, Vec _axis2 ){
			name_         = _name;
			nfold1_       = _nfold1;
			nfold2_       = _nfold2;
			axis1_        = _axis1;
			axis2_        = _axis2;
			move1_        = _axis1;
			move2_        = _axis2;
			nangle1_      = 360/nfold1_;
			nangle2_      = 360/nfold2_;
			angle1_       = 360.0/(Real)nfold1_;
			angle2_       = 360.0/(Real)nfold2_;
			intra_xform1_ = Xform( rotation_matrix_degrees(axis2_,angle2_), Vec(0,0,0) );
			intra_xform2_ = Xform( rotation_matrix_degrees(axis1_,angle1_), Vec(0,0,0) );
			intra_slide1_ = ( intra_xform1().R*axis1_ - axis1_ ).normalized();
			intra_slide2_ = ( intra_xform2().R*axis2_ - axis2_ ).normalized();
			intra_move1A_ =                 move1_;
			intra_move1B_ = intra_xform1_ * move1_;
			intra_move2A_ =                 move2_;
			intra_move2B_ = intra_xform2_ * move2_;
			alpha_        = angle_degrees(move1_,Vec(0,0,0),move2_);
			sin_alpha_    = sin(numeric::conversions::radians(alpha_));
			tan_alpha_    = tan(numeric::conversions::radians(alpha_));
		}
	};
	typedef utility::pointer::shared_ptr<Arch      > ArchOP;
	typedef utility::pointer::shared_ptr<Arch const> ArchCOP;

	class ArchT32 : public Arch { public: ArchT32(){ init( "T32",
		3, Vec( 1.00000000000000, 1.00000000000000, 1.00000000000000 ).normalized()   ,
	    2, Vec( 1.00000000000000, 0.00000000000000, 0.00000000000000 ).normalized()   ); } };
	class ArchT33 : public Arch { public: ArchT33(){ init( "T33",
		3, Vec( 1.00000000000000, 1.00000000000000, 1.00000000000000 ).normalized()   ,
	    3, Vec( 1.00000000000000, 1.00000000000000,-1.00000000000000 ).normalized()   ); } };
	class ArchO32 : public Arch { public: ArchO32(){ init( "O32",
		3, Vec( 1.00000000000000, 1.00000000000000, 1.00000000000000 ).normalized()   ,
	    2, Vec( 1.00000000000000, 1.00000000000000, 0.00000000000000 ).normalized()   ); } };
	class ArchO42 : public Arch { public: ArchO42(){ init( "O42",
		4, Vec( 1.00000000000000, 0.00000000000000, 0.00000000000000 ).normalized()   ,
	    2, Vec( 1.00000000000000, 1.00000000000000, 0.00000000000000 ).normalized()   ); } };
	class ArchO43 : public Arch { public: ArchO43(){ init( "O43",
		4, Vec( 1.00000000000000, 0.00000000000000, 0.00000000000000 ).normalized()   ,
	    3, Vec( 1.00000000000000, 1.00000000000000, 1.00000000000000 ).normalized()   ); } };
	class ArchI32 : public Arch { public: ArchI32(){ init( "I32",
		3, Vec( 0.93417235896272, 0.00000000000000, 0.35682208977309 ).normalized()   ,
	    2, Vec( 1.00000000000000, 0.00000000000000, 0.00000000000000 ).normalized()   ); } };
	class ArchI52 : public Arch { public: ArchI52(){ init( "I52",
		5, Vec( 0.85065080835204, 0.52573111211914, 0.00000000000000 ).normalized()   ,
	    2, Vec( 1.00000000000000, 0.00000000000000, 0.00000000000000 ).normalized()   ); } };
	class ArchI53 : public Arch { public: ArchI53(){ init( "I53",
		5, Vec( 0.85065080835204, 0.52573111211914, 0.00000000000000 ).normalized()   ,
	    3, Vec( 0.93417235896272, 0.00000000000000, 0.35682208977309 ).normalized()   ); } };


	// move1 slides along comp1 3fold
	// move2 slides perp to axes
	class ArchP6_32 : public Arch {
	public:
		ArchP6_32(){
			name_         = "P6_32";
			nfold1_       = 3;
			nfold2_       = 2;
			nangle1_      = 360/nfold1_;
			nangle2_      = 360/nfold2_;
			angle1_       = 360.0/(Real)nfold1_;
			angle2_       = 360.0/(Real)nfold2_;
			axis1_        = Vec(0,0,1);
			axis2_        = Vec(0,0,1);
			move1_        = Vec(0,0,1);
			move2_        = Vec(-1,0,0);
			intra_move1A_ =                 move2_;
			intra_move1B_ = intra_xform1_ * move2_;
			intra_move2A_ =                 move2_;
			intra_move2B_ = intra_xform2_ * move2_;

			alpha_        = angle_degrees(move1_,Vec(0,0,0),move2_);
			sin_alpha_    = sin(numeric::conversions::radians(alpha_));
			tan_alpha_    = tan(numeric::conversions::radians(alpha_));
			intra_xform1_ = Xform( rotation_matrix_degrees(axis2_,angle2_), Vec(0,0,0) );
			intra_xform2_ = Xform( rotation_matrix_degrees(axis1_,angle1_), Vec(0,0,0) );
		    intra_slide1_ = (intra_xform1_*move2_-move2_).normalized(); // move1 != move2!!!
		    intra_slide2_ = (intra_xform2_*move2_-move2_).normalized(); // same
		}
		virtual Real intra_slidedis_to_component_disp(int slidecomp, int /*dispcomp*/, Real d) const {
			if(1==slidecomp){
				return d / 2.0;
			} else {
				return d / 2.0 * 1.1547005383792517; // == 1/sin(pi/3)
			}
			utility_exit_with_message("bad input");
			return 0.0;
		}
		// virtual Real disp_dof1(Real const & omega, Real const & slidedis) const {
		// 	Real const gamma = numeric::conversions::radians(omega-90.0);
		// 	Real const sin_gamma=sin(gamma);
		// 	// Real const cos_gamma=cos(gamma);
		// 	Real const x = slidedis*sin_gamma;
		// 	// Real const y = slidedis*cos_gamma;
		// 	Real const w = x/sin_alpha_;
		// 	// Real const z = x/tan_alpha_;
		// 	// dcmp1 = w;
		// 	// dcmp2 = y+z;
		// 	return w;
		// }
		// virtual Real disp_dof2(Real const & omega, Real const & slidedis) const {
		// 	Real const gamma = numeric::conversions::radians(omega-90.0);
		// 	Real const sin_gamma=sin(gamma);
		// 	Real const cos_gamma=cos(gamma);
		// 	Real const x = slidedis*sin_gamma;
		// 	Real const y = slidedis*cos_gamma;
		// 	// Real const w = x/sin_alpha_;
		// 	Real const z = x/tan_alpha_;
		// 	// dcmp1 = w;
		// 	// dcmp2 = y+z;
		// 	return y+z;
		// }
		virtual void get_disp_dofs(Real const & omega, Real const & slidedis, Real & dcmp1_out, Real & dcmp2_out) const {
			Real const gamma = numeric::conversions::radians(omega-90.0);
			Real const sin_gamma=sin(gamma);
			Real const cos_gamma=cos(gamma);
			Real const x = slidedis*sin_gamma;
			Real const y = slidedis*cos_gamma;
			Real const w = x/sin_alpha_;
			Real const z = x/tan_alpha_;
			dcmp1_out = w;   //disp_dof1(omega,slidedis);
			dcmp2_out = y+z; //disp_dof2(omega,slidedis);
		}

		virtual Real intra_check_dof1(Real const & /*d1*/, Real const & d2) const { return d2; }
		virtual Real intra_check_dof2(Real const & /*d1*/, Real const & d2) const { return d2; }
		virtual Xform generator1(Real /*a1*/, Real /*a2*/, Real /*d1*/, Real /*d2*/) const { return Xform::rot_deg(axis1_,angle1_,Vec(0,0,0)); }
		virtual Xform generator2(Real /*a1*/, Real /*a2*/, Real /*d1*/, Real   d2  ) const { return Xform::rot_deg(axis2_,angle2_,move2_*d2); }

		virtual int max_ori() const { return 180; }
	};


	// move1 slides along comp1 3fold
	// move2 slides perp to axes
	class ArchP6m_32d : public ArchP6_32 {
	public:
		ArchP6m_32d(){
			name_         = "P6m_32d";
			nfold1_       = 3;
			nfold2_       = 2;
			nangle1_      = 360/nfold1_;
			nangle2_      = 360/nfold2_;
			angle1_       = 360.0/(Real)nfold1_;
			angle2_       = 360.0/(Real)nfold2_;
			axis1_        = Vec(0,0,1);
			axis2_        = Vec(0,0,1);
			move1_        = Vec(0,0,1);
			move2_        = Vec(1,0,0);
			intra_move1A_ =                 move2_;
			intra_move1B_ = intra_xform1_ * move2_;
			intra_move2A_ =                 move2_;
			intra_move2B_ = intra_xform2_ * move2_;

			alpha_        = angle_degrees(move1_,Vec(0,0,0),move2_);
			sin_alpha_    = sin(numeric::conversions::radians(alpha_));
			tan_alpha_    = tan(numeric::conversions::radians(alpha_));
			intra_xform1_ = Xform( rotation_matrix_degrees(axis2_,angle2_), Vec(0,0,0) );
			intra_xform2_ = Xform( rotation_matrix_degrees(axis1_,angle1_), Vec(0,0,0) );
		    intra_slide1_ = (intra_xform1_*move2_-move2_).normalized(); // move1 != move2!!!
		    intra_slide2_ = (intra_xform2_*move2_-move2_).normalized(); // same
		}
		virtual bool dihedral () const { return true; }
		virtual bool dihedral2() const { return true; }

	};


	// move1 slides along comp1 3fold
	// move2 slides perp to axes
	class ArchP4_42 : public Arch {
	public:
		ArchP4_42(){
			name_         = "P4_42";
			nfold1_       = 4;
			nfold2_       = 2;
			nangle1_      = 360/nfold1_;
			nangle2_      = 360/nfold2_;
			angle1_       = 360.0/(Real)nfold1_;
			angle2_       = 360.0/(Real)nfold2_;
			axis1_        = Vec(0,0,1);
			axis2_        = Vec(0,0,1);
			move1_        = Vec(0,0,1);
			move2_        = Vec(0,1,0);
			intra_move1A_ =                 move2_;
			intra_move1B_ = intra_xform1_ * move2_;
			intra_move2A_ =                 move2_;
			intra_move2B_ = intra_xform2_ * move2_;

			alpha_        = angle_degrees(move1_,Vec(0,0,0),move2_);
			sin_alpha_    = sin(numeric::conversions::radians(alpha_));
			tan_alpha_    = tan(numeric::conversions::radians(alpha_));
			intra_xform1_ = Xform( rotation_matrix_degrees(axis2_,angle2_), Vec(0,0,0) );
			intra_xform2_ = Xform( rotation_matrix_degrees(axis1_,angle1_), Vec(0,0,0) );
		    intra_slide1_ = (intra_xform1_*move2_-move2_).normalized(); // move1 != move2!!!
		    intra_slide2_ = (intra_xform2_*move2_-move2_).normalized(); // same
		}
		virtual Real intra_slidedis_to_component_disp(int slidecomp, int /*dispcomp*/, Real d) const {
			if(1==slidecomp){
				return d / 2.0;
			} else {
				return d / 2.0 * 1.414214; // == 1/sin(pi/4)
			}
			utility_exit_with_message("bad input");
			return 0.0;
		}
		// virtual Real disp_dof1(Real const & omega, Real const & slidedis) const {
		// 	Real const gamma = numeric::conversions::radians(omega-90.0);
		// 	Real const sin_gamma=sin(gamma);
		// 	// Real const cos_gamma=cos(gamma);
		// 	Real const x = slidedis*sin_gamma;
		// 	// Real const y = slidedis*cos_gamma;
		// 	Real const w = x/sin_alpha_;
		// 	// Real const z = x/tan_alpha_;
		// 	// dcmp1 = w;
		// 	// dcmp2 = y+z;
		// 	return w;
		// }
		// virtual Real disp_dof2(Real const & omega, Real const & slidedis) const {
		// 	Real const gamma = numeric::conversions::radians(omega-90.0);
		// 	Real const sin_gamma=sin(gamma);
		// 	Real const cos_gamma=cos(gamma);
		// 	Real const x = slidedis*sin_gamma;
		// 	Real const y = slidedis*cos_gamma;
		// 	// Real const w = x/sin_alpha_;
		// 	Real const z = x/tan_alpha_;
		// 	// dcmp1 = w;
		// 	// dcmp2 = y+z;
		// 	return y+z;
		// }
		virtual void get_disp_dofs(Real const & omega, Real const & slidedis, Real & dcmp1_out, Real & dcmp2_out) const {
			Real const gamma = numeric::conversions::radians(omega-90.0);
			Real const sin_gamma=sin(gamma);
			Real const cos_gamma=cos(gamma);
			Real const x = slidedis*sin_gamma;
			Real const y = slidedis*cos_gamma;
			Real const w = x/sin_alpha_;
			Real const z = x/tan_alpha_;
			dcmp1_out = w;   //disp_dof1(omega,slidedis);
			dcmp2_out = y+z; //disp_dof2(omega,slidedis);
		}
		virtual Real intra_check_dof1(Real const & /*d1*/, Real const & d2) const { return d2; }
		virtual Real intra_check_dof2(Real const & /*d1*/, Real const & d2) const { return d2; }
		virtual Xform generator1(Real /*a1*/, Real /*a2*/, Real /*d1*/, Real /*d2*/) const { return Xform::rot_deg(axis1_,angle1_,Vec(0,0,0)); }
		virtual Xform generator2(Real /*a1*/, Real /*a2*/, Real /*d1*/, Real   d2  ) const { return Xform::rot_deg(axis2_,angle2_,move2_*d2); }

		virtual int max_ori() const { return 180; }
	};

	class ArchP4m_42d : public ArchP4_42 {
	public:
		ArchP4m_42d(){
			name_         = "P4m_42d";
			nfold1_       = 4;
			nfold2_       = 2;
			nangle1_      = 360/nfold1_;
			nangle2_      = 360/nfold2_;
			angle1_       = 360.0/(Real)nfold1_;
			angle2_       = 360.0/(Real)nfold2_;
			axis1_        = Vec(0,0,1);
			axis2_        = Vec(0,0,1);
			move1_        = Vec(0,0,1);
			move2_        = Vec(0,1,0);
			intra_move1A_ =                 move2_;
			intra_move1B_ = intra_xform1_ * move2_;
			intra_move2A_ =                 move2_;
			intra_move2B_ = intra_xform2_ * move2_;

			alpha_        = angle_degrees(move1_,Vec(0,0,0),move2_);
			sin_alpha_    = sin(numeric::conversions::radians(alpha_));
			tan_alpha_    = tan(numeric::conversions::radians(alpha_));
			intra_xform1_ = Xform( rotation_matrix_degrees(axis2_,angle2_), Vec(0,0,0) );
			intra_xform2_ = Xform( rotation_matrix_degrees(axis1_,angle1_), Vec(0,0,0) );
		    intra_slide1_ = (intra_xform1_*move2_-move2_).normalized(); // move1 != move2!!!
		    intra_slide2_ = (intra_xform2_*move2_-move2_).normalized(); // same
		}
		virtual bool dihedral () const { return true; }
		virtual bool dihedral2() const { return true; }

	};


	class ArchP4g_42d : public ArchP4_42 {
	public:
		ArchP4g_42d(){
			name_         = "P4g_42d";
			nfold1_       = 4;
			nfold2_       = 2;
			nangle1_      = 360/nfold1_;
			nangle2_      = 360/nfold2_;
			angle1_       = 360.0/(Real)nfold1_;
			angle2_       = 360.0/(Real)nfold2_;
			axis1_        = Vec(0,0,1);
			axis2_        = Vec(0,0,1);
			move1_        = Vec(0,0,1);
			move2_        = Vec(1,1,0);
			intra_move1A_ =                 move2_;
			intra_move1B_ = intra_xform1_ * move2_;
			intra_move2A_ =                 move2_;
			intra_move2B_ = intra_xform2_ * move2_;

			alpha_        = angle_degrees(move1_,Vec(0,0,0),move2_);
			sin_alpha_    = sin(numeric::conversions::radians(alpha_));
			tan_alpha_    = tan(numeric::conversions::radians(alpha_));
			intra_xform1_ = Xform( rotation_matrix_degrees(axis2_,angle2_), Vec(0,0,0) );
			intra_xform2_ = Xform( rotation_matrix_degrees(axis1_,angle1_), Vec(0,0,0) );
		    intra_slide1_ = (intra_xform1_*move2_-move2_).normalized(); // move1 != move2!!!
		    intra_slide2_ = (intra_xform2_*move2_-move2_).normalized(); // same
		}
		virtual bool dihedral () const { return true; }
		virtual bool dihedral2() const { return true; }

	};

/////////////////////////////////// OPTIONS ////////////////////////////{
	OPT_1GRP_KEY( Real , tcdock, intra )
	OPT_1GRP_KEY( Real , tcdock, intra1 )
	OPT_1GRP_KEY( Real , tcdock, intra2 )
	OPT_1GRP_KEY( Real , tcdock, termini_weight )
	OPT_1GRP_KEY( Real , tcdock, termini_cutoff )
	OPT_1GRP_KEY( Real , tcdock, termini_cutoff_short )
	OPT_1GRP_KEY( Real , tcdock, termini_dist_filter )
	OPT_1GRP_KEY( Real , tcdock, min_score_filter )
	OPT_1GRP_KEY( Integer , tcdock, termini_trim )
	OPT_1GRP_KEY( Integer , tcdock, nsamp1 )
	OPT_1GRP_KEY( Integer , tcdock, topx )
	OPT_1GRP_KEY( Integer , tcdock, peak_grid_size)
	OPT_1GRP_KEY( Integer , tcdock, peak_grid_smooth)
	OPT_1GRP_KEY( Integer , tcdock, dump_pdb )
	OPT_1GRP_KEY( IntegerVector , tcdock, dump_pdbs )
	OPT_1GRP_KEY( Boolean , tcdock, dump_pdb_generated_symmetry_no_symdef )

	OPT_1GRP_KEY( Boolean , tcdock, reverse )

	OPT_1GRP_KEY( Integer , tcdock, dump_pdb_generated_depth )
	OPT_1GRP_KEY( Real    , tcdock, dump_pdb_generated_radius )
	OPT_1GRP_KEY( Boolean , tcdock, dump_gz )
	OPT_1GRP_KEY( Integer , tcdock, dump_pdb_primary_sub )
	OPT_1GRP_KEY( IntegerVector , tcdock, dump_peak_grids )
	OPT_1GRP_KEY( String , tcdock, architecture )
	OPT_1GRP_KEY( FileVector, tcdock, component1 )
	OPT_1GRP_KEY( FileVector, tcdock, component2 )
	OPT_1GRP_KEY( Boolean , tcdock, usen1 )
	OPT_1GRP_KEY( Boolean , tcdock, usec1 )
	OPT_1GRP_KEY( Boolean , tcdock, require_exposed_termini )
	OPT_1GRP_KEY( Boolean , tcdock, dry_run )
	OPT_1GRP_KEY( Real , tcdock, term_min_expose )
	OPT_1GRP_KEY( Real , tcdock, term_max_angle )
	OPT_1GRP_KEY( String , tcdock, clash_atoms )
	OPT_1GRP_KEY( Real , tcdock, redundant_angle )
	// OPT_1GRP_KEY( Real , tcdock, redundant_dist )
	OPT_1GRP_KEY( Boolean , tcdock, fast_stage_one )
	OPT_1GRP_KEY( Integer , tcdock, max_linker_len )
	OPT_1GRP_KEY( Integer , tcdock, linker_lookup_radius )
	OPT_1GRP_KEY( Integer , tcdock, max_res )
	OPT_1GRP_KEY( Boolean , tcdock, cb_weight_secstruct )
	OPT_1GRP_KEY( Boolean , tcdock, cb_weight_average_degree )
	OPT_1GRP_KEY( Real    , tcdock, trim_floppy_termini )
	OPT_1GRP_KEY( Boolean , tcdock, debug )

	OPT_1GRP_KEY( Boolean , tcdock, require_inter_contact )
	OPT_1GRP_KEY( Boolean , tcdock, print_dumped_motifs )

	OPT_1GRP_KEY( File , tcdock, dump_transforms )

	OPT_1GRP_KEY( Real , tcdock, motif_hash_weight )
	OPT_1GRP_KEY( Real , tcdock, hash_speed_filter_hack )

	void register_options() {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		OPT( in::file::s );

		NEW_OPT( tcdock::architecture     ,"target arch"        , ""     );
		NEW_OPT( tcdock::component1            ,"file(s) for component 1"        , ""     );
		NEW_OPT( tcdock::component2            ,"file(s) for component 2"       , ""     );

		NEW_OPT( tcdock::intra            ,"weight for intra contacts"           ,  1.0   );
		NEW_OPT( tcdock::intra1           ,"weight for comp1 intra contacts"     ,  1.0   );
		NEW_OPT( tcdock::intra2           ,"weight for conp2 intra contacts"     ,  1.0   );
		NEW_OPT( tcdock::nsamp1           ,"check around X top lowres hits"      , 20000   );
		NEW_OPT( tcdock::topx             ,"output top X hits"                   , 10     );
		NEW_OPT( tcdock::peak_grid_size   ,"peak detect grid size (2*N+1)"       , 24     );
		NEW_OPT( tcdock::peak_grid_smooth ,"peak detect grid smooth (0+)"        ,  1     );
		NEW_OPT( tcdock::dump_pdb         ,"dump pdbs up to #"                       , 0  );
		NEW_OPT( tcdock::dump_pdbs        ,"dump particular pdbs"                    , -1 );
		NEW_OPT( tcdock::dump_pdb_generated_symmetry_no_symdef         ,"dump autogen sym pdbs" , false );
		NEW_OPT( tcdock::dump_pdb_generated_depth                      ,"dump autogen sym depth" , 8 );
		NEW_OPT( tcdock::dump_pdb_generated_radius                     ,"dump autogen sym radius" , 100.0 );
		NEW_OPT( tcdock::dump_gz          ,"dump gz pdbs"                        , true   );

		NEW_OPT( tcdock::reverse         ,""                     , false   );

		NEW_OPT( tcdock::dump_pdb_primary_sub ,"dump pdb primary subunit only"   , 0  );
		NEW_OPT( tcdock::dump_peak_grids  ,"dump specific grids grids"           , -1     );
		NEW_OPT( tcdock::termini_cutoff   ,"tscore = w*max(0,cut-x)"             , 20.0   );
		NEW_OPT( tcdock::termini_cutoff_short   ,"tscore = w*max(0,cut-x)"       ,  0.0   );
		NEW_OPT( tcdock::termini_weight   ,"tscore = w*max(0,cut-x)"             ,  0.0   );
		NEW_OPT( tcdock::termini_trim     ,"trim termini up to for termini score",  0     );
		NEW_OPT( tcdock::termini_dist_filter ,"trim termini up to for termini score",  9e9     );
		NEW_OPT( tcdock::min_score_filter ,"min total score",  0.0001     );
		NEW_OPT( tcdock::usec1            ,"use comp1 cterm"                     , true   );
		NEW_OPT( tcdock::usen1            ,"use comp1 nterm"                     , true   );
		NEW_OPT( tcdock::require_exposed_termini,""         , false );
		NEW_OPT( tcdock::dry_run          ,"no calculations, just setup"         , false  );
		NEW_OPT( tcdock::term_max_angle ,""             , 45.0  );
		NEW_OPT( tcdock::term_min_expose ,""            , 0.1  );
		NEW_OPT( tcdock::clash_atoms     ,""            , "BB"  );
		NEW_OPT( tcdock::redundant_angle ,""            , 0  );
		// NEW_OPT( tcdock::redundant_dist  ,""            , 0  );
		NEW_OPT( tcdock::fast_stage_one  ,"faster stage one, may miss some"   , false  );
		NEW_OPT( tcdock::max_linker_len  ,""   ,  10    );
		NEW_OPT( tcdock::linker_lookup_radius  ,""   , 1 );
		NEW_OPT( tcdock::max_res  ,""   , 99999 );
		NEW_OPT( tcdock::cb_weight_secstruct     ,""         , true );
		NEW_OPT( tcdock::cb_weight_average_degree,""         , true );
		NEW_OPT( tcdock::trim_floppy_termini ,"", 0.0 );
		NEW_OPT( tcdock::debug,""         , false );

		NEW_OPT( tcdock::require_inter_contact,""         , false );
		NEW_OPT( tcdock::print_dumped_motifs,""         , false );

		NEW_OPT( tcdock::dump_transforms,""         ,"" );

		NEW_OPT( tcdock::motif_hash_weight,""         , 1.0 );
		NEW_OPT( tcdock::hash_speed_filter_hack,""    , 10.0 );

	}

//////////////////////////////// UTILITY CRAP //////////////////////////{
static THREAD_LOCAL basic::Tracer TR( "symdock_enum" );


/////////////////////////// TCDOCK MISC ////////////////////////////////{
	void prune_cb_pairs(vector1<Vec> & cba, vector1<Vec> & cbb, vector1<Real> & wa_in, vector1<Real> & wb_in, Real CTD2) {
		vector1<Vec> a,b;
		vector1<Real> wa,wb;
		vector1<Real>::const_iterator iwa = wa_in.begin();
		for(vector1<Vec>::const_iterator ia = cba.begin(); ia != cba.end(); ++ia,++iwa) {
			vector1<Real>::const_iterator iwb = wb_in.begin();
			for(vector1<Vec>::const_iterator ib = cbb.begin(); ib != cbb.end(); ++ib,++iwb) {
				if( ib->distance_squared(*ia) < CTD2 ) {
					a.push_back(*ia);
					b.push_back(*ib);
					wa.push_back(*iwa);
					wb.push_back(*iwb);
				}
			}
		}
		cba = a;
		cbb = b;
		wa_in = wa;
		wb_in = wb;
	}

	struct LMAX {
		Real score,radius;
		int icmp2,icmp1,iori;
		bool reverse;
		LMAX() :	score(0),radius(0),icmp2(0),icmp1(0),iori(0) {}
		LMAX(Real _score, Real _radius, int _icmp2, int _icmp1, int _iori, bool _reverse) :
			score(_score),radius(_radius),icmp2(_icmp2),icmp1(_icmp1),iori(_iori),reverse(_reverse) {}
	};
	int compareLMAX(const LMAX a,const LMAX b) {
		return a.score > b.score;
	}


Real sphere_coverage(vector1<Vec> points, vector1<Xform> const & xforms){
	static vector1<Vec> sphere;
	if(sphere.size()==0){
		utility::io::izstream is;
		basic::database::open(is,"sampling/spheres/sphere_4092.dat");
		for(int i = 1; i <= 4092; ++i) {
			Real x,y,z;
			is >> x >> y >> z;
			sphere.push_back(Vec(x,y,z));
		}
		is.close();
	}
	Real tot = 0;
	for(vector1<Vec  >::const_iterator is = sphere.begin(); is != sphere.end(); ++is){
		Vec const & s(*is);
		bool covered = false;
		for(vector1<Vec>  ::const_iterator ip = points.begin(); ip != points.end(); ++ip){
			Vec const p( *ip / ip->length() );
			for(vector1<Xform>::const_iterator ix = xforms.begin(); ix != xforms.end(); ++ix){
				Xform const & x(*ix);
				if( s.distance_squared(x*p) < 0.002 ){
					covered = true;
					goto done;
				}
			}
		}
		done:
		if(covered) tot += 1.0;
	}
	return tot / 4092.0;
}

struct TCDock {
	private: ArchCOP arch_; public: Arch const & arch() const { return *arch_; }

	vector1<Real> cmp2mnpos_,cmp1mnpos_,cmp2mnneg_,cmp1mnneg_;
	ObjexxFCL::FArray2D<Real> cmp2cbpos_,cmp1cbpos_,cmp2cbneg_,cmp1cbneg_;

	ObjexxFCL::FArray3D<Real> gradii,gscore;

	Real cmp1diapos_,cmp1dianeg_,cmp2diapos_,cmp2dianeg_;

	Real cmp1scale_,cmp2scale_;

	core::pose::Pose cmp1olig_,cmp2olig_;
	core::pose::Pose const & cmp1_,cmp2_;
	vector1<Vec> cmp1pts_,cmp1cbs_,cmp2pts_,cmp2cbs_;
	vector1<Real> cmp1wts_,cmp2wts_;
	core::id::AtomID_Map<Real> cmp1wtmap_,cmp2wtmap_;

	string cmp1name_,cmp2name_;

	core::id::AtomID_Map<core::Real> clashmap1_,clashmap2_;

	protocols::sic_dock::SICFast sic_;
	protocols::sic_dock::RigidScoreCOP rigid_sfxn_;
	protocols::sic_dock::RigidScoreCOP final_sfxn_;
	protocols::sic_dock::LinkerScoreCOP lnscore_;
	protocols::sic_dock::scores::MotifHashRigidScoreCOP mhscore_;

	bool abort_;
	Size conf_count_;
	Real start_time_;

	utility::io::ozstream * xform_out_;

	TCDock(
		ArchCOP _arch,
		Pose const & _pose1,
		Pose const & _pose2,
		string _cmp1pdb,
		string _cmp2pdb
	):
		arch_(_arch),
		cmp1_(_pose1),
		cmp2_(_pose2),
		rigid_sfxn_(/* NULL */),
		final_sfxn_(/* NULL */),
		lnscore_(/* NULL */),
		mhscore_(/* NULL */),
		abort_(false),
		conf_count_(0),
		xform_out_(NULL)
	{
		#ifdef USE_OPENMP
		cout << "OMP info: " << num_threads() << " " << thread_num() << endl;
		#endif
		using basic::options::option;
		using namespace basic::options::OptionKeys;

		cout << "tcdock run " << arch().name() << " nang1 " << arch().nangle1() << " nang2 " << arch().nangle2() << endl;

		if(option[tcdock::dump_transforms].user()){
			string xfname = option[out::file::o]()+"/"+ObjexxFCL::string_of(option[tcdock::dump_transforms]())+"_"+arch().name()+"_"+utility::file_basename(_cmp1pdb)+"_"+utility::file_basename(_cmp2pdb)+".docks";
			cout << "dump xforms to " << xfname << endl;
			xform_out_ = new utility::io::ozstream(xfname);
			if(xform_out_) (*xform_out_) << arch().name() << " " << _cmp1pdb << " " << _cmp2pdb << endl;
		}

		cmp1olig_ = _pose1;
		cmp2olig_ = _pose2;
		// core::chemical::ResidueTypeSetCAP crs=core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD);
		// core::import_pose::pose_from_file(cmp1olig_,*crs,_cmp1pdb, core::import_pose::PDB_file);
		// core::import_pose::pose_from_file(cmp2olig_,*crs,_cmp2pdb, core::import_pose::PDB_file);

		cmp1diapos_=0.0,cmp1dianeg_=0.0,cmp2diapos_=0.0,cmp2dianeg_=0.0;
		for(Size i = 1; i <= cmp1olig_.size(); ++i) {
			int const natom = (cmp1olig_.residue(i).name3()=="GLY") ? 4 : 5;
			for(int j = 1; j <= natom; ++j) {
				Vec const & v( cmp1olig_.residue(i).xyz(j) );
				cmp1diapos_ = max(cmp1diapos_, v.z());
				cmp1dianeg_ = max(cmp1dianeg_,-v.z());
			}
		}
		for(Size i = 1; i <= cmp2olig_.size(); ++i) {
			int const natom = (cmp2olig_.residue(i).name3()=="GLY") ? 4 : 5;
			for(int j = 1; j <= natom; ++j) {
				Vec const & v( cmp2olig_.residue(i).xyz(j) );
				cmp2diapos_ = max(cmp2diapos_, v.z());
				cmp2dianeg_ = max(cmp2dianeg_,-v.z());
			}
		}


		// set up structure
		make_Cx(cmp1olig_,arch().nfold1());
		make_Cx(cmp2olig_,arch().nfold2());

		if(arch().dihedral1()){
			utility_exit_with_message("nto implemented!!");
			cout << "DIHEDRAL component 1" << endl;
			make_Cx(cmp1olig_,2,Vec(0,1,0));
		}
		if(arch().dihedral2()){
			cout << "DIHEDRAL component 2" << endl;
			make_Cx(cmp2olig_,2,Vec(0,1,0));
		}

		cmp1name_ = utility::file_basename(_cmp1pdb);
		cmp2name_ = utility::file_basename(_cmp2pdb);

	}
	void init(bool reverse){
		using basic::options::option;
		using namespace basic::options::OptionKeys;

		if(reverse) rot_pose(cmp1olig_,Vec(0,1,0),180.0);

		core::scoring::dssp::Dssp(cmp1olig_).insert_ss_into_pose(cmp1olig_);
		core::scoring::dssp::Dssp(cmp2olig_).insert_ss_into_pose(cmp2olig_);
		if(option[tcdock::debug]()) cmp1olig_.dump_pdb("cmp1olig.pdb");
		if(option[tcdock::debug]()) cmp2olig_.dump_pdb("cmp2olig.pdb");

		bool nt1good=1,nt2good=1,ct1good=1,ct2good=1;
		// protocols::sic_dock::termini_exposed(cmp1olig_,nt1good,ct1good);
		// protocols::sic_dock::termini_exposed(cmp2olig_,nt2good,ct2good);
		TR << cmp1name_ << " nterm 1: " << nt1good << endl;
		TR << cmp1name_ << " cterm 1: " << ct1good << endl;
		TR << cmp2name_ << " nterm 2: " << nt2good << endl;
		TR << cmp2name_ << " cterm 2: " << ct2good << endl;
		if(!option[tcdock::usec1].user()) {
			if( ct1good && nt2good ) option[tcdock::usec1](true);
			else                     option[tcdock::usec1](false);
		}
		if(!option[tcdock::usen1].user()) {
			if( ct2good && nt1good ) option[tcdock::usen1](true);
			else                     option[tcdock::usen1](false);
		}
		cout << "TERMINI " << cmp1name_ << " " << cmp2name_ << " " << option[tcdock::usec1]() << " " << option[tcdock::usen1]() << endl;
		if( !option[tcdock::usec1]() && !option[tcdock::usen1]() ) {
			if(option[tcdock::require_exposed_termini]()) {
				std::cout << "no compatible good termini" << std::endl;
				abort_ = true;
				return;
			}
		}

		if(cmp1name_.substr(cmp1name_.size()-3)==".gz" ) cmp1name_ = cmp1name_.substr(0,cmp1name_.size()-3);
		if(cmp2name_.substr(cmp2name_.size()-3)==".gz" ) cmp2name_ = cmp2name_.substr(0,cmp2name_.size()-3);
		if(cmp1name_.substr(cmp1name_.size()-4)==".pdb") cmp1name_ = cmp1name_.substr(0,cmp1name_.size()-4);
		if(cmp2name_.substr(cmp2name_.size()-4)==".pdb") cmp2name_ = cmp2name_.substr(0,cmp2name_.size()-4);

		alignaxis(cmp1olig_,arch().axis1(),Vec(0,0,1),Vec(0,0,0));
		alignaxis(cmp2olig_,arch().axis2(),Vec(0,0,1),Vec(0,0,0));

		// COMPUTE WEIGHTS BB, AND CBS
		using core::id::AtomID;
		core::pose::initialize_atomid_map(cmp1wtmap_,cmp1olig_, 0.0);
		core::pose::initialize_atomid_map(cmp2wtmap_,cmp2olig_, 0.0);
		core::pose::initialize_atomid_map(clashmap1_,cmp1olig_,-1.0);
		core::pose::initialize_atomid_map(clashmap2_,cmp2olig_,-1.0);
		TR.Warning << "skip O atoms!!!!" << endl;
		for(Size i12 = 0; i12 < 2; ++i12){
			Pose const & ptmp( i12?cmp1olig_:cmp2olig_ );
			for(Size i = 1; i <= ptmp.size(); ++i) {
				if(ptmp.residue(i).has("CB")) {
					Real wt = 1.0;
					if( option[tcdock::cb_weight_average_degree]()                          ) wt *= min(1.0,(Real)neighbor_count(ptmp,i)/20.0);
					if( option[tcdock::cb_weight_secstruct     ]() && ptmp.secstruct(i)=='L') wt /= 3.0; //TODO make option somehow
					(i12?cmp1cbs_:cmp2cbs_).push_back(Vec(ptmp.residue(i).xyz("CB")));
					(i12?cmp1wts_:cmp2wts_).push_back(wt);
					(i12?cmp1wtmap_:cmp2wtmap_)[AtomID(ptmp.residue(i).atom_index("CB"),i)] = wt;
				}
				int natom=0;
				if( "BB"    == option[tcdock::clash_atoms]() ) natom = (ptmp.residue(i).name3()=="GLY") ? 4 : 5;
				if( "HEAVY" == option[tcdock::clash_atoms]() ) natom = ptmp.residue(i).nheavyatoms();
				for(int j = 1; j <= natom; ++j) {
					if(j==4) continue; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					Vec const & v( ptmp.residue(i).xyz(j) );
					(i12?clashmap1_:clashmap2_)[AtomID(j,i)] = ptmp.residue(i).atom_type(j).lj_radius();
					(i12?cmp1pts_:cmp1pts_).push_back(v);
				}
			}
		}

		sic_.init(cmp1olig_,cmp2olig_,clashmap1_,clashmap2_);
		Real CTD(basic::options::option[basic::options::OptionKeys::sicdock::contact_dis]());
		Real CLD(basic::options::option[basic::options::OptionKeys::sicdock::  clash_dis]());

		TR << "creating score functions" << std::endl;
		protocols::sic_dock::RigidScoreCOP cbscore( protocols::sic_dock::RigidScoreOP( new protocols::sic_dock::CBScore(cmp1olig_,cmp2olig_,CLD,CTD,cmp1wtmap_,cmp2wtmap_) ) );
		protocols::sic_dock::JointScoreOP jtscore( new protocols::sic_dock::JointScore );
		protocols::sic_dock::JointScoreOP fnscore( new protocols::sic_dock::JointScore );
		rigid_sfxn_ = jtscore;
		final_sfxn_ = fnscore;

		if(option[basic::options::OptionKeys::mh::xform_score_data].user()){
			mhscore_ = protocols::sic_dock::scores::MotifHashRigidScoreCOP( protocols::sic_dock::scores::MotifHashRigidScoreOP( new protocols::sic_dock::scores::MotifHashRigidScore(cmp1olig_,cmp2olig_) ) );
			jtscore->add_score(mhscore_,option[tcdock::motif_hash_weight]());
			// Real const approx_surf = std::pow((Real)cmp1olig_.size(),0.666666) * std::pow((Real)cmp2olig_.size(),0.666666);
			// Real const hackcut = approx_surf/25.0 * option[tcdock::hash_speed_filter_hack]();
			// cout << "MotifHashRigidScore / JointScore Speed filter hack scorecut: " << hackcut << endl;
		} else {
			jtscore->add_score(cbscore ,1.0);
		}
		if(option[tcdock::max_linker_len]() > 0){
			string lnscore_tag = cmp1name_+"_"+cmp2name_+arch().name()+"_";
			lnscore_ = protocols::sic_dock::LinkerScoreCOP( protocols::sic_dock::LinkerScoreOP( new protocols::sic_dock::LinkerScore(cmp1olig_,cmp2olig_,option[tcdock::max_linker_len](),option[tcdock::linker_lookup_radius](),lnscore_tag) ) );
			jtscore->add_score(lnscore_,1.0);
			// fnscore->add_score(lnscore_,1.0);
		}


		cmp1mnpos_.resize   (arch().nangle1(),    0.0); cmp2mnpos_.resize   (arch().nangle2(),    0.0);
		cmp1cbpos_.dimension(arch().nangle1(),200,0.0); cmp2cbpos_.dimension(arch().nangle2(),200,0.0);
		cmp1mnneg_.resize   (arch().nangle1(),    0.0); cmp2mnneg_.resize   (arch().nangle2(),    0.0);
		cmp1cbneg_.dimension(arch().nangle1(),200,0.0); cmp2cbneg_.dimension(arch().nangle2(),200,0.0);
		gradii.dimension(arch().nangle2(),arch().nangle1(),arch().max_ori(),-9e9);
		gscore.dimension(arch().nangle2(),arch().nangle1(),arch().max_ori(),-9e9);

		core::scoring::methods::RG_Energy_Fast rgcalc;
		Real tmp1 = rgcalc.calculate_rg_score(cmp1olig_);
		Real tmp2 = rgcalc.calculate_rg_score(cmp2olig_);
		cmp1scale_ = 2.0*tmp1/(tmp1+tmp2);
		cmp2scale_ = 2.8*tmp2/(tmp1+tmp2);
	 }
	virtual ~TCDock() {
		if(xform_out_) delete xform_out_;
	 }
	int num_threads() {
		#ifdef USE_OPENMP
			return omp_get_max_threads();
		#else
			return 1;
		#endif
	 }
	int thread_num() {
		#ifdef USE_OPENMP
			return omp_get_thread_num();
		#else
			return 0;
		#endif
	 }
	void precompute_intra() {
		// return;
		using basic::options::option;
		using namespace basic::options::OptionKeys;
		Real const CONTACT_D  = basic::options::option[basic::options::OptionKeys::sicdock::contact_dis]();
		Real const CLASH_D    = basic::options::option[basic::options::OptionKeys::sicdock::  clash_dis]();
		Real const CONTACT_D2 = sqr(CONTACT_D);
		// Real const CLASH_D2	= sqr(CLASH_D);
		// compute high/low min dis for pent and cmp1 here, input to sicfast and don't allow any below
		cout << "precomputing one-component interactions every 1°" << endl;
		for(int i12 = 0; i12 < 2; ++i12) {
			Vec const axis = (i12?arch().axis1():arch().axis2());
			vector1<Real>             & cmpwts  ( i12?cmp1wts_  :cmp2wts_   );
			vector1<Real>             & cmpmnpos( i12?cmp1mnpos_:cmp2mnpos_ );
			vector1<Real>             & cmpmnneg( i12?cmp1mnneg_:cmp2mnneg_ );
			ObjexxFCL::FArray2D<Real> & cmpcbpos( i12?cmp1cbpos_:cmp2cbpos_ );
			ObjexxFCL::FArray2D<Real> & cmpcbneg( i12?cmp1cbneg_:cmp2cbneg_ );

			protocols::sic_dock::SICFast sic;
			sic.init( i12? cmp1olig_   :cmp2olig_,
				      i12? cmp1olig_   :cmp2olig_,
				      i12? clashmap1_:clashmap2_,
				      i12? clashmap1_:clashmap2_);
			protocols::sic_dock::CBScore sfxn(i12?cmp1olig_:cmp2olig_,i12?cmp1olig_:cmp2olig_,CLASH_D,CONTACT_D);

			for(int ipn = 0; ipn < 2; ++ipn) {
				#ifdef USE_OPENMP
				#pragma omp parallel for schedule(dynamic,1)
				#endif
				for(int icmp = 0; icmp < (i12?arch().nangle1():arch().nangle2()); icmp+=1) {
					Xform const x_to_related_comp =        ( i12 ? arch().intra_xform1() : arch().intra_xform2() );
					Vec   const sicaxis = (ipn?1.0:-1.0) * ( i12 ? arch().intra_slide1() : arch().intra_slide2() );
					vector1<Vec> ptsA,cbA;
					if(i12) get_cmp1(icmp,ptsA,cbA); else get_cmp2(icmp,ptsA,cbA);
					vector1<Vec> ptsB(ptsA),cbB(cbA);
					assert( cbA.size() == cmpwts.size() && cbB.size() == cmpwts.size() );
					Xform const Xintra( i12? arch().intra_xform1() : arch().intra_xform2() );
					std::transform(ptsB.begin(),ptsB.end(),ptsB.begin(),Xintra);
					std::transform( cbB.begin(), cbB.end(), cbB.begin(),Xintra);
					// for(vector1<Vec>::iterator i = ptsB.begin(); i != ptsB.end(); ++i) *i = Xintra*(*i);
					// for(vector1<Vec>::iterator i =  cbB.begin(); i !=  cbB.end(); ++i) *i = Xintra*(*i);
					Vec const TMP_axis2 = x_to_related_comp.R * (i12?arch().axis1():arch().axis2());


					Real score = 0;
					Xform xa,xb;
					xa.R = rotation_matrix_degrees(TMP_axis2,(Real)icmp) * x_to_related_comp.R; /// FIXME
					xb.R = rotation_matrix_degrees(    axis ,(Real)icmp)                      ;
					Real const d = -protocols::sic_dock::slide_into_contact_and_score(sic,sfxn,xa,xb,sicaxis,score);
					if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for cmppos! "+ObjexxFCL::string_of(icmp));
					if(fabs(d) > 9e8){
						utility_exit_with_message("precompute_intra error");
					}

					Real mindcomp1 = (ipn?-1.0:1.0) * arch().intra_slidedis_to_component_disp(i12?1:2,1,d);
					Real mindcomp2 = (ipn?-1.0:1.0) * arch().intra_slidedis_to_component_disp(i12?1:2,2,d);
					if(i12)	(ipn?cmpmnpos:cmpmnneg)[icmp+1] = mindcomp1;
					else    (ipn?cmpmnpos:cmpmnneg)[icmp+1] = mindcomp2;

					// if( i12==0 && ipn==0 && icmp==33 ){
					// 	cout << i12 <<" "<<ipn<<" "<<icmp<<" slide axis: " << (i12?"tri ":"dim ") << sicaxis << endl;
					// 	cout << i12 <<" "<<ipn<<" "<<icmp<<" slide dist: " << d << endl;
					// 	cout << i12 <<" "<<ipn<<" "<<icmp<<" dcomp1 " << mindcomp1 << endl;
					// 	cout << i12 <<" "<<ipn<<" "<<icmp<<" dcomp2 " << mindcomp2 << endl;

					// 	// Pose tmp1(cmp1olig_);
					// 	// rot_pose(tmp1,arch().axis1(),icmp);
					// 	// tmp1.dump_pdb("test1A.pdb");
					// 	// protocols::sic_dock::xform_pose(tmp1,Xintra);
					// 	// trans_pose(tmp1,d*sicaxis);
					// 	// tmp1.dump_pdb("test1B.pdb");

					// 	// Pose tmp2(cmp2olig_);
					// 	// rot_pose(tmp2,arch().axis2(),icmp);
					// 	// tmp2.dump_pdb("test2A.pdb");
					// 	// protocols::sic_dock::xform_pose(tmp2,Xintra);
					// 	// trans_pose(tmp2,d*sicaxis);
					// 	// tmp2.dump_pdb("test2B.pdb");


					// 	Pose tmp1(cmp1olig_);
					// 	Pose tmp2(cmp2olig_);
					// 	  rot_pose(tmp1,arch().axis1(),icmp);
					// 	  rot_pose(tmp2,arch().axis2(),icmp);
					// 	trans_pose(tmp1,mindcomp1*arch().move2());
					// 	trans_pose(tmp2,mindcomp2*arch().move2());
					// 	if(i12==1) tmp1.dump_pdb("test1A.pdb");
					// 	if(i12==0) tmp2.dump_pdb("test2A.pdb");
					// 	protocols::sic_dock::xform_pose(tmp1,Xintra);
					// 	protocols::sic_dock::xform_pose(tmp2,Xintra);
					// 	if(i12==1) tmp1.dump_pdb("test1B.pdb");
					// 	if(i12==0) tmp2.dump_pdb("test2B.pdb");

					// 	utility_exit_with_message("test precomp intra");
					// }

					Vec move_intra1 = i12 ? arch().intra_move1A() : arch().intra_move2A();
					Vec move_intra2 = i12 ? arch().intra_move1B() : arch().intra_move2B();

					// move 2nd CBs into contact
					std::transform(cbB.begin(),cbB.end(),cbB.begin(), Xform(Mat::identity(),-d*sicaxis) );

					vector1<Real> wA(cmpwts),wB(cmpwts);
					prune_cb_pairs(cbA,cbB,wA,wB,CONTACT_D2); // mutates inputs!!!!!!!!!!
					// Real lastcbc = 9e9;
					for(int i = 1; i <= 200; ++i) {
						Real cbc = 0.0;
						vector1<Real>::const_iterator iwa=wA.begin(), iwb=wB.begin();
						for( vector1<Vec>::const_iterator ia=cbA.begin(), ib=cbB.begin(); ia != cbA.end(); ++ia,++ib,++iwa,++iwb) {
							cbc += CBScore_dist_score( ia->distance_squared(*ib) , CLASH_D, CONTACT_D ) * (*iwa) * (*iwa);
						}
						// assert(lastcbc >= cbc);
						// lastcbc = cbc;
						if(cbc < 0.0) utility_exit_with_message("negative score!");
						(ipn?cmpcbpos:cmpcbneg)(icmp+1,i) = cbc;
						if(cbc==0.0){
							for(int jj=i; jj <= 200; ++jj) (ipn?cmpcbpos:cmpcbneg)(icmp+1,jj) = 0.0;
							break;
						}
						for(vector1<Vec>::iterator iv = cbA.begin(); iv != cbA.end(); ++iv) *iv = (*iv) + (ipn?0.1:-0.1)*move_intra1;
						for(vector1<Vec>::iterator iv = cbB.begin(); iv != cbB.end(); ++iv) *iv = (*iv) + (ipn?0.1:-0.1)*move_intra2;
					}
				}
			}
		}
	 }
	void dump_onecomp() {
		using	namespace	basic::options;
		using	namespace	basic::options::OptionKeys;
		using utility::file_basename;
		cout << "dumping 1D stats: " << option[out::file::o]()+"/"+cmp2name_+"_POS_1D.dat" << endl;
		{ utility::io::ozstream out(option[out::file::o]()+"/"+cmp2name_+"_POS_1D.dat");
			for(int i = 1; i <= arch().nangle2(); ++i) out << i << " " << cmp2cbpos_(i,1) << " " << cmp2cbpos_(i,2) << " " << cmp2cbpos_(i,3) << " " << cmp2cbpos_(i,4) << endl;
			out.close(); }
		{ utility::io::ozstream out(option[out::file::o]()+"/"+cmp2name_+"_NEG_1D.dat");
			for(int i = 1; i <= arch().nangle2(); ++i) out << i << " " << cmp2cbneg_(i,1) << " " << cmp2cbneg_(i,2) << " " << cmp2cbneg_(i,3) << " " << cmp2cbneg_(i,4) << endl;
			out.close(); }
		{ utility::io::ozstream out(option[out::file::o]()+"/"+cmp1name_+"_POS_1D.dat");
			for(int i = 1; i <= arch().nangle1(); ++i) out << i << " " << cmp1cbpos_(i,1) << " " << cmp1cbpos_(i,2) << " " << cmp1cbpos_(i,3) << " " << cmp1cbpos_(i,4) << endl;
			out.close(); }
		{ utility::io::ozstream out(option[out::file::o]()+"/"+cmp1name_+"_NEG_1D.dat");
			for(int i = 1; i <= arch().nangle1(); ++i) out << i << " " << cmp1cbneg_(i,1) << " " << cmp1cbneg_(i,2) << " " << cmp1cbneg_(i,3) << " " << cmp1cbneg_(i,4) << endl;
			out.close(); }
		utility_exit_with_message("1COMP!");
	 }
	void get_cmp1(Real acmp1, vector1<Vec> & pts, vector1<Vec> & cbs ) {
		pts.resize(cmp1pts_.size()); std::transform(cmp1pts_.begin(),cmp1pts_.end(),pts.begin(), arch().xform1(acmp1,0.0) );
		cbs.resize(cmp1cbs_.size()); std::transform(cmp1cbs_.begin(),cmp1cbs_.end(),cbs.begin(), arch().xform1(acmp1,0.0) );
	 }
	void get_cmp2(Real acmp2, vector1<Vec> & pts, vector1<Vec> & cbs ){
		pts.resize(cmp2pts_.size()); std::transform(cmp2pts_.begin(),cmp2pts_.end(),pts.begin(), arch().xform2(acmp2,0.0) );
		cbs.resize(cmp2cbs_.size()); std::transform(cmp2cbs_.begin(),cmp2cbs_.end(),cbs.begin(), arch().xform2(acmp2,0.0) );
	 }
	void get_best_sub1_contact_delta_rotations( Real icmp1, Real dcmp1, Real icmp2, Real dcmp2, Real & out_delta_ang1,  Real & out_delta_ang2){
		Real const contact_dis = basic::options::option[basic::options::OptionKeys::sicdock::contact_dis]();
		Real const contact_dis2 = contact_dis*contact_dis;
		Xform x1 = arch().xform1(icmp1,dcmp1);
		Xform x2 = arch().xform2(icmp2,dcmp2);
		int best1=0,best2=0,bestcount=0;
		for(int ir1 = 0; ir1 < arch().nfold1(); ++ir1){
			for(int ir2 = 0; ir2 < arch().nfold2(); ++ir2){
				int ccount = 0;
				for(Size ir=1; ir<=cmp1cbs_.size()/arch().nfold1(); ++ir){
					Vec const cb1( x1*(cmp1cbs_[ir]) );
					for(Size jr=1; jr<=cmp2cbs_.size()/arch().nfold2(); ++jr){
						Vec const cb2( x2*(cmp2cbs_[jr]) );
						if( cb1.distance_squared(cb2) < contact_dis2 ) ccount++;
					}
				}
				if(ccount > bestcount){
					bestcount = ccount;
					best1 = ir1;
					best2 = ir2;
				}
				x2.R = rotation_matrix_degrees(arch().axis2(),(Real)arch().nangle2()) * x2.R;
			}
			x1.R = rotation_matrix_degrees(arch().axis1(),(Real)arch().nangle1()) * x1.R;

		}
		out_delta_ang1 = (Real)best1*(Real)arch().nangle1();
		out_delta_ang2 = (Real)best2*(Real)arch().nangle2();
	 }
	void get_contacts_by_component(Real const & ang1, Real const & ang2, Real const & dis1, Real const & dis2,
		                          int  & nresiface, int  & nresiface0, int  & nresiface1, int  & nresiface2,
		                          int  & ncontacts, int  & ncontacts0, int  & ncontacts1, int  & ncontacts2,
		                          Real & nweighted, Real & nweighted0, Real & nweighted1, Real & nweighted2
	 ){
		Real contact_dis2 = sqr(basic::options::option[basic::options::OptionKeys::sicdock::contact_dis]());
		vector1<Vec> cb1 = cmp1cbs_; std::transform(cb1.begin(),cb1.end(),cb1.begin(),arch().xform1(ang1,dis1));
		vector1<Vec> cb2 = cmp2cbs_; std::transform(cb2.begin(),cb2.end(),cb2.begin(),arch().xform2(ang2,dis2));
		vector1<Vec> cb1intra = cb1; std::transform(cb1intra.begin(),cb1intra.end(),cb1intra.begin(),arch().intra_xform1());
		vector1<Vec> cb2intra = cb2; std::transform(cb2intra.begin(),cb2intra.end(),cb2intra.begin(),arch().intra_xform2());
		ncontacts0=0; ncontacts1=0; ncontacts2=0; nweighted0=0; nweighted1=0; nweighted2=0;

		std::set<vector1<Vec>::const_iterator > s,s0,s1,s2;
		Real CTD(basic::options::option[basic::options::OptionKeys::sicdock::contact_dis]());
		Real CLD(basic::options::option[basic::options::OptionKeys::sicdock::  clash_dis]());


		for(vector1<Vec>::const_iterator i = cb1.begin(); i != cb1.end(); ++i){
		for(vector1<Vec>::const_iterator j = cb2.begin(); j != cb2.end(); ++j){
			Real d2 = i->distance_squared(*j);
			if(d2 < contact_dis2){
				ncontacts0 += 1;
				nweighted0 += CBScore_dist_score(d2,CLD,CTD);
				s.insert(i); s.insert(j); s0.insert(i); s0.insert(j);
			}
		}}
		for(vector1<Vec>::const_iterator i = cb1     .begin(); i != cb1     .end(); ++i){
		for(vector1<Vec>::const_iterator j = cb1intra.begin(); j != cb1intra.end(); ++j){
			Real d2 = i->distance_squared(*j);
			if(d2 < contact_dis2){
				ncontacts1 += 1;
				nweighted1 += CBScore_dist_score(d2,CLD,CTD);
				s.insert(i); s.insert(j); s1.insert(i); s1.insert(j);
			}
		}}
		for(vector1<Vec>::const_iterator i = cb2     .begin(); i != cb2     .end(); ++i){
		for(vector1<Vec>::const_iterator j = cb2intra.begin(); j != cb2intra.end(); ++j){
			Real d2 = i->distance_squared(*j);
			if(d2 < contact_dis2){
				ncontacts2 += 1;
				nweighted2 += CBScore_dist_score(d2,CLD,CTD);
				s.insert(i); s.insert(j); s2.insert(i); s2.insert(j);
			}
		}}
		ncontacts = ncontacts0 + ncontacts1 + ncontacts2;
		nweighted = nweighted0 + nweighted1 + nweighted2;
		nresiface  = s .size();
		nresiface0 = s0.size();
		nresiface1 = s1.size();
		nresiface2 = s2.size();
	 }
	void dump_pdb(int icmp2, int icmp1, int iori, string fname, bool reverse, bool dumpsym) {
		using basic::options::option;
		using namespace basic::options::OptionKeys;
		using ObjexxFCL::string_of;
		Real d=0,dcmp2=0,dcmp1=0,rscore=0,cmp1cbc=0,cmp2cbc=0;
		dock_get_geom(icmp2,icmp1,iori,d,dcmp2,dcmp1,rscore,cmp2cbc,cmp1cbc);
		Real da1=0,da2=0;
		get_best_sub1_contact_delta_rotations(icmp1,dcmp1,icmp2,dcmp2,da1,da2);


		Pose p1,p2,symm;
		#ifdef USE_OPENMP
		#pragma omp critical
		#endif
		{


			// // residues in starting pose... silly
			// Size Ncmp1 = ((arch().dihedral1()) ? cmp1olig_.size()/2 : cmp1olig_.size()) / arch().nfold1();
			// Size Ncmp2 = ((arch().dihedral2()) ? cmp2olig_.size()/2 : cmp2olig_.size()) / arch().nfold2();
			// p1.append_residue_by_jump(cmp1olig_.residue(1),1,"","",false);
			// for(Size i = 2; i <= Ncmp1; ++i) {
			// 	if(p1.residue(i-1).is_terminus()||p1.residue(i-1).is_ligand())
			// 	     p1.append_residue_by_jump(cmp1olig_.residue(i),1);
			// 	else p1.append_residue_by_bond(cmp1olig_.residue(i));
			// }
			// if(true && !p1.residue(       1      ).is_lower_terminus()) add_lower_terminus_type_to_pose_residue(p1,      1       );
			// if(true && !p1.residue(p1.size()).is_upper_terminus()) add_upper_terminus_type_to_pose_residue(p1,p1.size());
			// p2.append_residue_by_jump(cmp2olig_.residue(1),1,"","", true ); // other chain iff not dumping sym complex (too many chains)
			// for(Size i = 2; i <= Ncmp2; ++i) {
			// 	if(p2.residue(p2.size()).is_terminus()||p2.residue(p2.size()).is_ligand())
			// 	     p2.append_residue_by_jump(cmp2olig_.residue(i),1);
			// 	else p2.append_residue_by_bond(cmp2olig_.residue(i));
			// }
			// if(true && !p2.residue(       1      ).is_lower_terminus()) add_lower_terminus_type_to_pose_residue(p2,      1       );
			// if(true && !p2.residue(p2.size()).is_upper_terminus()) add_upper_terminus_type_to_pose_residue(p2,p2.size());
			// p1.dump_pdb("testB1.pdb");
			// p2.dump_pdb("testB2.pdb");

			p1 = cmp1_;
			p2 = cmp2_;
			if(reverse) rot_pose(p1,Vec(0,1,0),180.0);
			alignaxis(p1,arch().axis1(),Vec(0,0,1),Vec(0,0,0));
			alignaxis(p2,arch().axis2(),Vec(0,0,1),Vec(0,0,0));
			// p1.dump_pdb("testA1.pdb");
			// p2.dump_pdb("testA2.pdb");
			// utility_exit_with_message("foo");
 		}


		protocols::sic_dock::xform_pose(p1,arch().xform1(icmp1,dcmp1));
		protocols::sic_dock::xform_pose(p2,arch().xform2(icmp2,dcmp2));


		// make joint structure
		{
			symm.append_residue_by_jump(p1.residue(1),1,"","",false);
			for(Size i = 2; i <= p1.size(); ++i) {
				if(p1.residue(i-1).is_upper_terminus()||
				   p1.residue(i  ).is_lower_terminus()||
				   p1.residue(i-1).is_ligand()||
				   p1.residue(i  ).is_ligand()
				){
					symm.append_residue_by_jump(p1.residue(i),1);
				} else {
					symm.append_residue_by_bond(p1.residue(i));
				}
			}
			symm.append_residue_by_jump(p2.residue(1),1,"","", true ); // other chain iff not dumping sym complex (too many chains)
			for(Size i = 2; i <= p2.size(); ++i) {
				if(p2.residue(i-1).is_upper_terminus()||
				   p2.residue(i  ).is_lower_terminus()||
				   p2.residue(i-1).is_ligand()||
				   p2.residue(i  ).is_ligand()
				){
					symm.append_residue_by_jump(p2.residue(i),1);
				} else {
					symm.append_residue_by_bond(p2.residue(i));
				}
			}
		}

		if(mhscore_){
			utility::io::ozstream out(option[out::file::o]()+"/"+fname+"_motifs.pdb");
			Pose a(cmp1olig_),b(cmp2olig_);
			xform_pose(a,arch().xform1(icmp1,dcmp1));
			xform_pose(b,arch().xform2(icmp2,dcmp2));
			int tmp;
			mhscore_->dump_matching_motifs(a,b,out,tmp,NULL,option[tcdock::print_dumped_motifs]());
			out.close();
		}

		// symm.dump_pdb("test_"+string(dumpsym?"presym":"nosym")+".pdb");
		if(dumpsym) {
				option[ basic::options::OptionKeys::symmetry::symmetry_definition]("input/"+arch().name()+"_tcdock.sym");
			core::pose::symmetry::make_symmetric_pose(symm);
			// symm.dump_pdb("test_sym.pdb");

			if("P6_32"==arch().name()){

				using namespace core::conformation::symmetry;
				SymmetricConformation const & symm_conf ( dynamic_cast<SymmetricConformation const & > ( symm.conformation() ) );
				SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );
				std::map< Size, SymDof > dofs ( symm_info->get_dofs() );
				std::map< Size, SymDof >::iterator it;
				std::map< Size, SymDof >::iterator it_begin = dofs.begin();
				std::map< Size, SymDof >::iterator it_end = dofs.end();
				for(it=it_begin; it != it_end; ++it){
				 	core::kinematics::Jump j(symm.jump(it->first));
				 	Vec    t = j.get_translation();
				 	if( t.length() > 0.0001 ){
					 	Vec newt = Vec(0,-dcmp2,0);
					 	j.set_translation(newt);
					 	symm.set_jump(it->first,j);
					 	Vec delta = p2.xyz(core::id::AtomID(1,1))-symm.xyz(core::id::AtomID(1,p1.size()+1));
					 	trans_pose(symm, delta ,p1.size()+1,p1.size()+p2.size());
					 	// cout << "trans jump! " << j << " " << it->first << " " << dcmp1 << " " << dcmp2 << " " << newt-t << " " << delta << endl;
					}
				}
			}
		}

		{
			utility::io::ozstream out(option[out::file::o]()+"/"+fname);
			core::io::pdb::dump_pdb(symm,out);
			out.close();
		}

		if(
			!dumpsym && option[tcdock::dump_pdb_generated_symmetry_no_symdef]() ){
			string end = ".pdb";
			if(fname.substr(fname.size()-3,3)==".gz") end = ".pdb.gz";
			utility::io::ozstream out(option[out::file::o]()+"/"+fname+"_autosym"+end);
			Xform xA = arch().generator1(icmp1,icmp2,dcmp1,dcmp2);
			Xform xB = arch().generator2(icmp1,icmp2,dcmp1,dcmp2);
			vector1<Xform> xforms;
			std::back_insert_iterator<vector1<Xform> > back_it(xforms);
			numeric::expand_xforms(back_it,xA,xB,option[tcdock::dump_pdb_generated_depth](),option[tcdock::dump_pdb_generated_radius]());
			int xcount = 0;
			Xform xlast;
			Xform Xdihedral( rotation_matrix_degrees(Vec(0,1,0),180.0),Vec(0,0,0) );
			for(vector1<Xform>::const_iterator ix = xforms.begin(); ix != xforms.end(); ++ix){
				bool is_dihedral = arch().dihedral();
				for(int idh = 0; idh <= is_dihedral; ++idh){
					Xform x = *ix;
					x = x * (idh? Xdihedral : Xform() );
					protocols::sic_dock::xform_pose(symm, x * ~xlast);
					xlast = x;
					out<<"MODEL S"+ObjexxFCL::string_of(++xcount)<<endl;
					core::io::pdb::dump_pdb(symm,out);
					out<<"ENDMDL"<<endl;
					// cout << xcount << " " << *ix << endl;
				}
			}
			out.close();
		}

		if(lnscore_){
			lnscore_->dump_linkers( arch().xform1(icmp1,dcmp1), arch().xform2(icmp2,dcmp2), option[out::file::o]()+"/"+fname );
		}

		// find contacts
	 }
	Real __dock_base__(int icmp2,int icmp1,int iori,Real&dori,Real&dcmp2,Real&dcmp1,Real&rscore,Real&cmp2cbc,Real&cmp1cbc,bool cache=true){
		if(!cache || gradii(icmp2+1,icmp1+1,iori+1)==-9e9 ) {
			using basic::options::option;
			using namespace basic::options::OptionKeys;

			if(rscore!=-12345.0) gscore(icmp2+1,icmp1+1,iori+1) = 0.0;
			gradii(icmp2+1,icmp1+1,iori+1) = 9e9;

			rscore = 0.0;
 			cmp1cbc = 0.0;
			cmp2cbc = 0.0;

			icmp1 = ( icmp1 % arch().nangle1() + arch().nangle1() ) % arch().nangle1();
			icmp2 = ( icmp2 % arch().nangle2() + arch().nangle2() ) % arch().nangle2();
			iori  = ( iori  % arch().max_ori() + arch().max_ori() ) % arch().max_ori();

			Vec sicaxis = arch().slide_axis(iori);
			Xform xa,xb;
			xa.R = rotation_matrix_degrees(arch().axis1(),(Real)icmp1);
			xb.R = rotation_matrix_degrees(arch().axis2(),(Real)icmp2);

			if( ! arch().ori_is_allowed( (Real)iori + arch().alpha_/2.0 ) ) return 0.0;

			Real contact = -sic_.slide_into_contact(xa,xb,sicaxis);
			if(fabs(contact) > 9e8 || contact > 0 )	return 0.0;

			arch().get_disp_dofs((Real)iori,contact,dcmp1,dcmp2);
			Real dintra1 = arch().intra_check_dof1(dcmp1,dcmp2);
			Real dintra2 = arch().intra_check_dof2(dcmp1,dcmp2);
			gradii(icmp2+1,icmp1+1,iori+1) = contact;

			int ipn = dintra1 > 0 ? 1 : -1;
			Real const cmp1mn = (ipn>0?cmp1mnpos_:cmp1mnneg_)[icmp1+1];
			Real const cmp2mn = (ipn>0?cmp2mnpos_:cmp2mnneg_)[icmp2+1];
			if( ipn*dintra1 < ipn*cmp1mn || ipn*dintra2 < ipn*cmp2mn ) {
				Real const dmncmp1 = cmp1mn / dintra1 * contact;
				Real const dmncmp2 = cmp2mn / dintra2 * contact;
				if( dmncmp1 > 0.0 || dmncmp2 > 0.0 ) utility_exit_with_message("dmncmp1/dmncmp2 sign mismatch "+ObjexxFCL::string_of(dmncmp1)+" "+ObjexxFCL::string_of(dmncmp2)+" "+ObjexxFCL::string_of(cmp1mn)+" "+ObjexxFCL::string_of(cmp2mn)+" "+ObjexxFCL::string_of(dintra1)+" "+ObjexxFCL::string_of(dintra2));
				gradii(icmp2+1,icmp1+1,iori+1) = min(dmncmp2,dmncmp1);
				Real scale = min(dmncmp2,dmncmp1)/contact;
				dintra1  *= scale;
				dintra2  *= scale;
				dcmp1    *= scale;
				dcmp2    *= scale;
				contact  *= scale;
			}
			dori = contact;
			if(rscore==-12345.0) return rscore;

			for(Real slidedis = contact; slidedis > contact - 3.0; slidedis -= 0.2){
				xa.t = -slidedis*sicaxis; // move into slid position
				Real score = rigid_sfxn_->score(xa,xb);

				// ++count;
				// if(score > 30 && score > old) cout << iori << " " << slidedis << " " << score << " " << old << " " << count << endl;
				// old = score;

				if(option[tcdock::require_inter_contact]() && score < 3.0) continue;
				int didx2 = ipn*(int)(dintra2-cmp2mn)*10+1; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				int didx1 = ipn*(int)(dintra1-cmp1mn)*10+1; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				if(didx1 < 1 || didx2 < 1) utility_exit_with_message("bad");
				if( 0 < didx1 && didx1 <= 200 ) score += option[tcdock::intra]()*option[tcdock::intra1]()*(ipn>0?cmp1cbpos_:cmp1cbneg_)(icmp1+1,didx1);
				if( 0 < didx2 && didx2 <= 200 ) score += option[tcdock::intra]()*option[tcdock::intra2]()*(ipn>0?cmp2cbpos_:cmp2cbneg_)(icmp2+1,didx2);
				Real tdis, tscore( termini_score(dcmp1,icmp1,dcmp2,icmp2,tdis) );
				if(option[tcdock::termini_dist_filter]() < tdis) continue;
				score += tscore;

				if(arch().dihedral1()){
					utility_exit_with_message("-tcdock:dihedral only tested with D2 -- will be comp. 2");
					bool dihedral_too_close = false;
					vector1<Vec> cbsA(cmp1cbs_.size()); std::transform(cmp1cbs_.begin(),cmp1cbs_.end(),cbsA.begin(), arch().xform1(icmp1,dcmp1) );
					vector1<Vec> cbsB(cmp1cbs_.size()); std::transform(    cbsA.begin(),    cbsA.end(),cbsB.begin(), Xform(rotation_matrix_degrees(Vec(0,1,0),180.0),Vec(0,0,0)) );
					for(vector1<Vec>::const_iterator i = cbsA.begin(); i != cbsA.end(); ++i){
						for(vector1<Vec>::const_iterator j = cbsB.begin(); j != cbsB.end(); ++j){
							dihedral_too_close |= ( i->distance_squared(*j) < 36.0 );
						}
					}
					if(dihedral_too_close) continue;
				}
				if(arch().dihedral2()){
					bool dihedral_too_close = false;
					vector1<Vec> cbsA(cmp1cbs_.size()); std::transform(cmp1cbs_.begin(),cmp1cbs_.end(),cbsA.begin(), arch().xform1(icmp1,dcmp1) );
					vector1<Vec> cbsB(cmp1cbs_.size()); std::transform(    cbsA.begin(),    cbsA.end(),cbsB.begin(), Xform(rotation_matrix_degrees(Vec(0,1,0),180.0),Vec(0,0,0)) );
					for(vector1<Vec>::const_iterator i = cbsA.begin(); i != cbsA.end(); ++i){
						for(vector1<Vec>::const_iterator j = cbsB.begin(); j != cbsB.end(); ++j){
							dihedral_too_close |= ( i->distance_squared(*j) < 36.0 );
						}
					}
					if(dihedral_too_close) continue;
				}
				if(score > rscore){
					rscore = score;
					dori = slidedis;
					arch().get_disp_dofs((Real)iori,slidedis,dcmp1,dcmp2);
					if( 0 < didx1 && didx1 <= 200 ) cmp1cbc = (ipn>0?cmp1cbpos_:cmp1cbneg_)(icmp1+1,didx1);
					if( 0 < didx2 && didx2 <= 200 ) cmp2cbc = (ipn>0?cmp2cbpos_:cmp2cbneg_)(icmp2+1,didx2);
				}
			}
			if(xform_out_){
				Xform a1 = arch().xform1(icmp1,dcmp1)*alignaxis_xform(arch().axis1(),Vec(0,0,1),Vec(0,0,0));
				Xform a2 = arch().xform2(icmp2,dcmp2)*alignaxis_xform(arch().axis2(),Vec(0,0,1),Vec(0,0,0));
				(*xform_out_) << a1 << endl << a2 << endl;
			}

			gscore(icmp2+1,icmp1+1,iori+1) = rscore;

		}


		return gscore(icmp2+1,icmp1+1,iori+1);
	 }
	void dock_no_score(int icmp2,int icmp1,int iori){
		Real dori=0,dcmp2=0,dcmp1=0,rscore=0,cmp2cbc=0,cmp1cbc=0;
		rscore = -12345.0;
		__dock_base__(icmp2,icmp1,iori,dori,dcmp2,dcmp1,rscore,cmp2cbc,cmp1cbc,true);
	 }
	Real dock_score(int icmp2,int icmp1,int iori,Real&rscore,Real&cmp2cbc,Real&cmp1cbc){
		Real dori,dcmp2,dcmp1;
		return __dock_base__(icmp2,icmp1,iori,dori,dcmp2,dcmp1,rscore,cmp2cbc,cmp1cbc,false);
	 }
	Real dock_score(int icmp2,int icmp1,int iori){
		Real dori=0,dcmp2=0,dcmp1=0,rscore=0,cmp2cbc=0,cmp1cbc=0;
		return __dock_base__(icmp2,icmp1,iori,dori,dcmp2,dcmp1,rscore,cmp2cbc,cmp1cbc,false);
	 }
	Real dock_get_geom(int icmp2,int icmp1,int iori,Real&dori,Real&dcmp2,Real&dcmp1,Real&rscore,Real&cmp2cbc,Real&cmp1cbc){
		return __dock_base__(icmp2,icmp1,iori,dori,dcmp2,dcmp1,rscore,cmp2cbc,cmp1cbc,false);
	 }
	Real min_termini_dis(Real d1, Real a1, Real d2, Real a2){
		using basic::options::option;
		using namespace basic::options::OptionKeys;
		bool usen1 = option[tcdock::usen1]();
		bool usec1 = option[tcdock::usec1]();
		Real td = 9e9;
		for(int t1 = 0; t1 <= option[tcdock::termini_trim](); ++t1) {
			Vec n1_0 = cmp1olig_.xyz(core::id::AtomID(1,1+t1));
			for(int t2 = 0; t2 <= option[tcdock::termini_trim](); ++t2) {
				Vec n2_0 = cmp2olig_.xyz(core::id::AtomID(1,1+t2));
				for(int t3 = 0; t3 <= option[tcdock::termini_trim](); ++t3) {
					Vec c1_0 = cmp1olig_.xyz(core::id::AtomID(3,cmp1olig_.size()-t3));
					for(int t4 = 0; t4 <= option[tcdock::termini_trim](); ++t4) {
						Vec c2_0 = cmp2olig_.xyz(core::id::AtomID(3,cmp2olig_.size()-t4));
						for(int i1 = 0; i1 < arch().nfold1(); ++i1){
							Vec n1 = d1*arch().move1() + rotation_matrix_degrees(arch().axis1(),a1 + 360.0/(Real)arch().nfold1()*(Real)i1) * n1_0;
							Vec c1 = d1*arch().move1() + rotation_matrix_degrees(arch().axis1(),a1 + 360.0/(Real)arch().nfold1()*(Real)i1) * c1_0;
							for(int i2 = 0; i2 < arch().nfold2(); ++i2){
								Vec n2 = d2*arch().move2() + rotation_matrix_degrees(arch().axis2(),a2 + 360.0/(Real)arch().nfold2()*(Real)i2) * n2_0;
								Vec c2 = d2*arch().move2() + rotation_matrix_degrees(arch().axis2(),a2 + 360.0/(Real)arch().nfold2()*(Real)i2) * c2_0;
								if(usen1) td = min( n1.distance(c2), td );
								if(usec1) td = min( c1.distance(n2), td );
							}
						}
					}
				}
			}
		}
		return td>9e8? -1: td;
	 }
	Real min_termini_proj(Real a1, Real a2, Vec sicaxis){
		using basic::options::option;
		using namespace basic::options::OptionKeys;
		//		return 0.0;
		bool usen1 = option[tcdock::usen1]();
		bool usec1 = option[tcdock::usec1]();
		Real td = 9e9;
		int t1=1;
		Vec n1_0 = cmp1olig_.xyz(core::id::AtomID(1,1+t1));
		int t2=1;
		Vec n2_0 = cmp2olig_.xyz(core::id::AtomID(1,1+t2));
		int t3=0;
		Vec c1_0 = cmp1olig_.xyz(core::id::AtomID(3,cmp1olig_.size()-t3));
		int t4=0;
		Vec c2_0 = cmp2olig_.xyz(core::id::AtomID(3,cmp2olig_.size()-t4));
		for(int i1 = 0; i1 < arch().nfold1(); ++i1){
			Vec n1 = rotation_matrix_degrees(arch().axis1(),a1 + 360.0/(Real)arch().nfold1()*(Real)i1) * n1_0;
			Vec c1 = rotation_matrix_degrees(arch().axis1(),a1 + 360.0/(Real)arch().nfold1()*(Real)i1) * c1_0;
			for(int i2 = 0; i2 < arch().nfold2(); ++i2){
				Vec n2 = rotation_matrix_degrees(arch().axis2(),a2 + 360.0/(Real)arch().nfold2()*(Real)i2) * n2_0;
				Vec c2 = rotation_matrix_degrees(arch().axis2(),a2 + 360.0/(Real)arch().nfold2()*(Real)i2) * c2_0;
				if(usen1) td = min( projperp(sicaxis,n1-c2).length(), td );
				if(usec1) td = min( projperp(sicaxis,c1-n2).length(), td );
			}
		}
		//std::cerr << td << std::endl;
		return td;
	 }
	Real termini_score(Real d1, Real a1, Real d2, Real a2, Real & td){
		using basic::options::option;
		using namespace basic::options::OptionKeys;
		td = min_termini_dis(d1,a1,d2,a2);
		return option[tcdock::termini_weight]() * max(0.0,option[tcdock::termini_cutoff]()-max(option[tcdock::termini_cutoff_short](),td));
	 }
	void run() {
		vector1<LMAX> local_maxima;

		using basic::options::option;
		using namespace basic::options::OptionKeys;

		if( ! option[tcdock::reverse]() ){
			// init(true);
			init(false);
			gather_hits(local_maxima,false);
		} else {
			// init(false);
			init(true);
			// rot_pose(cmp1olig_,Vec(0,1,0),180.0);
			gather_hits(local_maxima,true);
		}
		std::sort(local_maxima.begin(),local_maxima.end(),compareLMAX);
		dump_top_hits(local_maxima);
	 }
	void get_moments(vector1<Vec> const & pts, Real & m1, Real & m2, Real & m3, Real & m4){
		m1=0; m2=0; m3=0; m4=0;
		for(vector1<Vec>::const_iterator i = pts.begin(); i != pts.end(); ++i){
			m1 += i->length();
		}
		m1 /= (Real)pts.size();
		for(vector1<Vec>::const_iterator i = pts.begin(); i != pts.end(); ++i){
			Real const r(i->length()-m1);
			m2 += r*r;
			m3 += r*r*r;
			m4 += r*r*r*r;
		}
		m2 /= (Real)pts.size();
		m3 /= (Real)pts.size();
		m4 /= (Real)pts.size();
		m2 = std::pow(m2,0.5);
		m3 = std::pow(fabs(m3),1.0/3.0) * (m3>0?1.0:-1.0);
		m4 = std::pow(m4,0.25);
	 }
	void gather_hits(vector1<LMAX> & local_maxima, bool const & reverse) {
		using basic::options::option;
		using namespace basic::options::OptionKeys;
		using namespace core::id;
		using numeric::conversions::radians;
		if(abort_) {
			cout << "abort run" << endl;
			return;
		}
		// Real const CONTACT_D  = basic::options::option[basic::options::OptionKeys::sicdock::contact_dis]();
		// Real const CLASH_D    = basic::options::option[basic::options::OptionKeys::sicdock::	clash_dis]();
		// Real const CONTACT_D2 = sqr(CONTACT_D);
		// Real const CLASH_D2   = sqr(CLASH_D);
		// Pose const cmp1init(cmp1olig_);
		// Pose const cmp2init(cmp2olig_);
		{                                                         // 3deg loop
			// Real max1=0;
			// dump_onecomp();
			// cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << endl;
			// cout << " !!!!!!!!!!!!!!!!no precomp intra! !!!!!!!!!!!!!!!!!!!!!!!!!!! " << endl;
			// cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << endl;
			// if(arch().name()=="P6_32"){
			// 	cout << "NOT CONSIDERING INTRA CONTACTS FOR " << arch().name() << endl;
			// } else {
			precompute_intra();
			// }

			// start_time_ = time_highres();
			// conf_count_ = 0;

			// utility_exit_with_message("TESTING P6");

			///////////////////////////// STAGE 1 ////////////////////////////////////////
			cout << "main loop 1 over icmp2, icmp1, iori every 3 degrees" << endl;
			Real max_score = 0;
			int max_i1=0,max_i2=0,max_io=0;
			if(option[tcdock::fast_stage_one]()){
				for(int icmp1 = 0; icmp1 < arch().nangle1(); icmp1+=3) {
					if(icmp1%15==0 && icmp1!=0){
						// Real rate = Real(conf_count_) / Real(time_highres()-start_time_);
						cout<<" lowres dock "
						    <<cmp1name_<<" "
						    <<cmp2name_<<" "
						    << (reverse?"R ":"F ")
						    <<I(2,100*icmp1/arch().nangle1())
						    <<"% done, max_score: "
						    <<F(10,3,max_score)
						    // <<" rate: "<<rate<<" confs/sec"
						    // <<" "<<rate/Real(num_threads())<<" confs/sec/thread"
						    <<endl;
					}
					#ifdef USE_OPENMP
					#pragma omp parallel for schedule(dynamic,1)
					#endif
					for(int icmp2 = 0; icmp2 < arch().nangle2(); icmp2+=3) {
						vector1<Vec> pb,cbb; get_cmp2(icmp2,pb,cbb);
						int iori = -1, stg = 1, fs1_count=0;	bool newstage = true;
						while(stg < 5) {
							if(newstage) {
								fs1_count = 0;
								iori = (int)(((stg>2)?270.0:90.0)+angle_degrees(arch().axis1(),Vec(0,0,0),arch().axis2())/2.0);
								iori = (iori / 3) * 3; // round to closest multiple of angle incr
								if(stg==2||stg==4) iori -= 3;
								newstage = false;
								// cout << "NEWSTAGE" << endl;
							} else {
								iori += (stg%2==0) ? -3 : 3;
							}
							iori = iori % 360;
							if(iori > arch().max_ori()) break;
							Real const score = dock_score(icmp2,icmp1,iori);
							#ifdef USE_OPENMP
							#pragma omp critical
							#endif
							{
								if(score > max_score) max_score = score;
							}
							if( -12345==score || ++fs1_count == 30 ){ stg++; newstage=true; continue; }
						}
					}
				}
			} else {
				//	JACOB!!!!
				vector1<numeric::xyzVector<int> > tasks;
				int div = 1;
				if(option[tcdock::dump_transforms].user()) div = 5;
				for(int icmp1 = 0; icmp1 < arch().nangle1()/div; icmp1+=3) {
					if( arch().dihedral1() && icmp1 != 0 && icmp1!=arch().nangle1()/2) continue;
					for(int icmp2 = 0; icmp2 < arch().nangle2()/div; icmp2+=3) {
						if( arch().dihedral2() && icmp2 != 0 && icmp2!=arch().nangle2()/2) continue;
						for(int iori = 0; iori < arch().max_ori(); iori+=3) {
							tasks.push_back( numeric::xyzVector<int>(icmp1,icmp2,iori) );
						}
					}
				}
				#ifdef USE_OPENMP
				#pragma omp parallel for schedule(dynamic,1)
				#endif
				for(int i = 1; i <= (int)tasks.size(); ++i){
					int icmp1 = tasks[i].x();
					int icmp2 = tasks[i].y();
					int iori  = tasks[i].z();
					if(i%28800==0){
						// Real rate = Real(conf_count_) / Real(time_highres()-start_time_);
						cout<<" lowres dock "
						    <<cmp1name_<<" "
						    <<cmp2name_<<" "
						    <<I(2,100*icmp1/arch().nangle1())
						    <<"% done, max_score: "
						    <<F(10,3,max_score) << endl;
						    // <<" rate: "<<rate<<" confs/sec"
						    // <<" "<<rate/num_threads()<<" confs/sec/thread"<<" " <<num_threads()
						if(mhscore_) cout<< " hash lookups:  " << F(10,3,(Real)mhscore_->nhashlookups()/1000.0) << "K";
						cout<< " confs sampled: " << F(10,3,conf_count_/1000.0) << "K"
						    << endl;
					}

					Real const score = dock_score(icmp2,icmp1,iori);
					#ifdef USE_OPENMP
					#pragma omp critical
					#endif
					{
						++conf_count_;
						if(score > max_score){
							max_score = score;
							// max_i1 = icmp1;
							// max_i2 = icmp2;
							// max_io = iori;
						}
					}
				}
			}
			if(max_score<0.00001){
				TR << "max_score<0.00001, cmp1 or cmp2 too large -> no contacts?";
				return;
			}
			cout << "MAX3 " << max_score << " " << max_i1 << " " << max_i2 << " " << max_io << endl;
		}

		if(xform_out_) xform_out_->close();
		if(xform_out_) delete xform_out_;
		xform_out_ = NULL;


		///////////////////////////////// STAGE 2 ////////////////////////////////////
		utility::vector1<vector1<int> > cmp2lmx,cmp1lmx,orilmx; { // set up work for main loop 2
			Real topX_3 = -9e6;
			vector1<Real> cbtmp;
			for(Size i = 0; i < gscore.size(); ++i) if(gscore[i] > 0) cbtmp.push_back(gscore[i]);
			std::sort(cbtmp.begin(),cbtmp.end());
			topX_3 = cbtmp[max(1,(int)cbtmp.size()-option[tcdock::nsamp1]()+1)];
			cout << "scanning top "<<min(option[tcdock::nsamp1](),(int)cbtmp.size())<<" with score3 >= " << topX_3 << endl;
			for(int icmp2 = 0; icmp2 < arch().nangle2(); icmp2+=3) {
				if( arch().dihedral2() && icmp2 != 0 && icmp2!=arch().nangle2()/2) continue;
				for(int icmp1 = 0; icmp1 < arch().nangle1(); icmp1+=3) {
					if( arch().dihedral1() && icmp1 != 0 && icmp1!=arch().nangle1()/2) continue;
					for(int iori = 0; iori < arch().max_ori(); iori+=3) {
						if( gscore(icmp2+1,icmp1+1,iori+1) >= topX_3) {
							vector1<int> cmp2,cmp1,ori;
							for(int i = -1; i <= 1; ++i) cmp2.push_back( (icmp2+i+arch().nangle2())%arch().nangle2() );
							for(int j = -1; j <= 1; ++j) cmp1.push_back( (icmp1+j+arch().nangle1())%arch().nangle1() );
							for(int k = -1; k <= 1; ++k) ori.push_back( (iori+k+arch().max_ori())%arch().max_ori() );
							cmp2lmx.push_back(cmp2);
							cmp1lmx.push_back(cmp1);
							orilmx.push_back(ori);
						}
					}
				}
			}
		}
		{                                                         //main loop 2
			#ifdef USE_OPENMP
			#pragma omp parallel for schedule(dynamic,1)
			#endif
			for(int ilmx = 1; ilmx <= (int)cmp2lmx.size(); ++ilmx)  {       //  MAIN LOOP 2
				if( (ilmx-1)%(option[tcdock::nsamp1]()/10)==0 && ilmx!=1){
					cout<<" highres dock "<<cmp1name_<<" "<<(Real(ilmx-1)/(option[tcdock::nsamp1]()/100))<<"% done"<<endl;
				}
				for(vector1<int>::const_iterator picmp2 = cmp2lmx[ilmx].begin(); picmp2 != cmp2lmx[ilmx].end(); ++picmp2) {
					int icmp2 = *picmp2;
					if( arch().dihedral2() && icmp2 != 0 && icmp2!=arch().nangle2()/2) continue;
					for(vector1<int>::const_iterator picmp1 = cmp1lmx[ilmx].begin(); picmp1 != cmp1lmx[ilmx].end(); ++picmp1) {
						int icmp1 = *picmp1;
						if( arch().dihedral1() && icmp1 != 0 && icmp1!=arch().nangle1()/2) continue;
						for(vector1<int>::const_iterator piori = orilmx[ilmx].begin(); piori != orilmx[ilmx].end(); ++piori) {
							int iori = *piori;
							dock_score(icmp2,icmp1,iori);
							#ifdef USE_OPENMP
							#pragma omp critical
							#endif
							++conf_count_;
						}
					}
				}
			}
		}
		cout << "gather local minima" << endl;
		{                             // get local radial disp maxima (minima, really)
			Real highscore = -9e9;
			for(int i = 1; i <=(int) gradii.size1(); ++i){
				for(int j = 1; j <= (int)gradii.size2(); ++j){
					for(int k = 1; k <= (int)gradii.size3(); ++k){
						Real const val = gradii(i,j,k);
						if( val < -9e8 ) continue;

						// Real nbmax = -9e9, scoremax = -9e9;
						// if( option[tcdock::geometric_minima_only]() ){
						// 	// int nedge = 0;
						// 	for(int di = -1; di <= 1; ++di){
						// 		for(int dj = -1; dj <= 1; ++dj){
						// 			for(int dk = -1; dk <= 1; ++dk){
						// 				if(di==0 && dj==0 && dk==0) continue;
						// 				int i2 = (i+di+gradii.size1()-1)%gradii.size1()+1;
						// 				int j2 = (j+dj+gradii.size2()-1)%gradii.size2()+1;
						// 				int k2 = (k+dk+gradii.size3()-1)%gradii.size3()+1;
						// 				Real const nbval = gradii(i2,j2,k2);
						// 				nbmax = max(nbmax,nbval);
						// 			}
						// 		}
						// 	}
						// }
						// if( option[tcdock::score_minima_only]() ){
						// 	// int nedge = 0;
						// 	for(int di = -1; di <= 1; ++di){
						// 		for(int dj = -1; dj <= 1; ++dj){
						// 			for(int dk = -1; dk <= 1; ++dk){
						// 				if(di==0 && dj==0 && dk==0) continue;
						// 				int i2 = (i+di+gradii.size1()-1)%gradii.size1()+1;
						// 				int j2 = (j+dj+gradii.size2()-1)%gradii.size2()+1;
						// 				int k2 = (k+dk+gradii.size3()-1)%gradii.size3()+1;
						// 				Real const nbval = gradii(i2,j2,k2);
						// 				nbmax = max(nbmax,nbval);
						// 			}
						// 		}
						// 	}
						// }

						// if( /*nbmax != -9e9 &&*/ val >= nbmax ) {
							Real score = gscore(i,j,k); //dock_score(i-1,j-1,k-1);
							if(score <= option[tcdock::min_score_filter]()) continue;
							local_maxima.push_back( LMAX(score,gradii(i,j,k),i-1,j-1,k-1,reverse) );
							highscore = max(score,highscore);
						// }
					}
				}
			}
		}
		if(mhscore_) cout << " hash lookups:  " << F(10,3,(Real)mhscore_->nhashlookups()/1000.0) << "K confs sampled: " << F(10,3,conf_count_/1000.0) << "K" << endl;

	 }
	void dump_top_hits(vector1<LMAX> & local_maxima){
		using basic::options::option;
		using namespace basic::options::OptionKeys;
		using namespace core::id;
		using numeric::conversions::radians;

		int npad = (cmp1name_+"_"+cmp2name_).size()+10;
		if( local_maxima.size() == 0 ){
			TR << endl << "no maxima!" << endl;
			return;
		}
		cout << "N maxima: " << local_maxima.size() << ", best score: " << local_maxima.front().score << endl;
		string pad = ""; for(int i = 1; i <= npad; ++i) pad += " ";
		cout << "  no  tag   " << pad << "    score   sc/nc sc/nwtd     rsc   nr   nc    nwtd    rsc0  nr0  nc0   nwtd0    cb1  nr1  nc1   nwtd1    cb2  nr2  nc2   nwtd2   diam   cover   tdis   n1  a1       r1   n2  a2       r2 ori   mom1    mom2    mom3    mom4  v0.4  v0.8  v1.2  v1.6  v2.0  v2.4  v2.8  v3.2  v3.6  v4.0 ";
		rigid_sfxn_->show(cout);
		final_sfxn_->show(cout);
		cout << "  no" << endl;
		vector1<Vec> dumpedit;
		int ilm = 0, nout = 0;
		Real redundant_thresh2 = option[tcdock::redundant_angle]() * option[tcdock::redundant_angle]();
		while( ++ilm <= (int)local_maxima.size() && nout < option[tcdock::topx]() ){
		// for(Size ilm = 1; ilm <= min(local_maxima.size(),(Size)option[tcdock::topx]()); ++ilm) { // dump top hit info
			LMAX const & h(local_maxima[ilm]);

			// redundency check
			bool toosimilar = false;
			for(vector1<Vec>::const_iterator i = dumpedit.begin(); i != dumpedit.end(); ++i){
				Real x = cmp1scale_*(Real)h.icmp1;
				Real y = cmp2scale_*(Real)h.icmp2;
				Real z = (cmp1scale_+cmp2scale_)/2.0*h.iori;
				if( i->distance_squared( Vec(x,y,z) ) < redundant_thresh2 ){
					toosimilar = true;
				}
			}
			if(toosimilar) continue;

			Real d=0,dcmp2=0,dcmp1=0,rscore=0,cmp1cbc=0,cmp2cbc=0;
			int N = option[tcdock::peak_grid_size]();
			ObjexxFCL::FArray3D<Real> grid(2*N+1,2*N+1,2*N+1,0.0);
			#ifdef USE_OPENMP
			#pragma omp parallel for schedule(dynamic,1)
			#endif
			for(int di = -N; di <= N; ++di) {
				for(int dj = -N; dj <= N; ++dj) {
					for(int dk = -N; dk <= N; ++dk) {
						if( Vec(di,dj,dk).length() > (Real)N+0.5 ) {
							grid(di+N+1,dj+N+1,dk+N+1) = -9e9;
						} else {
							int i = (h.icmp2+di+gscore.size1())%gscore.size1();
							int j = (h.icmp1+dj+gscore.size2())%gscore.size2();
							int k = (h.iori +dk+gscore.size3())%gscore.size3();
							if( gradii(i+1,j+1,k+1) < 0 ) dock_no_score(i,j,k);
							grid(di+N+1,dj+N+1,dk+N+1) = gradii(i+1,j+1,k+1);
						}
					}
				}
			}
			vector1<Real> ffhist(10,0);
			int const Nedge = option[tcdock::peak_grid_smooth]();
			// #ifdef USE_OPENMP
			// #pragma omp parallel for schedule(dynamic,1)
			// #endif
			for(int ifh = 0; ifh < (int)ffhist.size(); ++ifh) {
				protocols::sic_dock::flood_fill3D(N+1,N+1,N+1, grid, grid(N+1,N+1,N+1)-0.4 );
				Real count = 0;
				// ObjexxFCL::FArray3D<Real> grid2(grid);
				for(int i = 1+Nedge; i <= (int)grid.size1()-Nedge; ++i) {
					for(int j = 1+Nedge; j <= (int)grid.size2()-Nedge; ++j) {
						for(int k = 1+Nedge; k <= (int)grid.size3()-Nedge; ++k) {
							if( grid(i,j,k)!=grid(N+1,N+1,N+1) ) continue;
							int ninside = 0;
							for(int di = -Nedge; di <= Nedge; ++di) {
								for(int dj = -Nedge; dj <= Nedge; ++dj) {
									for(int dk = -Nedge; dk <= Nedge; ++dk) {
										ninside += grid(i+di,j+dj,k+dk)==grid(N+1,N+1,N+1);
									}
								}
							}
							Real w = max(0.0, 1.0 - Vec(N+1-i,N+1-j,N+1-k).length() / (Real)N );
							count += w*((Real)ninside/(Real)((2*Nedge+1)*(2*Nedge+1)*(2*Nedge+1)));
							// cerr << F(7,3,((Real)ninside/(Real)((2*Nedge+1)*(2*Nedge+1)*(2*Nedge+1)))) << " " << F(5,3,w) << " " << I(2,ninside) << endl;
							// grid2(i,j,k) += allgood;
						}
					}
				}
				ffhist[ifh+1] = pow(count,1.0/3.0);
				// vector1<int> dumpg = option[tcdock::dump_peak_grids]();
				// #ifdef USE_OPENMP
				// #pragma omp critical
				// #endif
				// if( std::find(dumpg.begin(),dumpg.end(),ilm)!=dumpg.end() ) {
					// utility::io::ozstream o(("grid_"+ObjexxFCL::string_of(ilm)+"_"+ObjexxFCL::string_of(ifh)+".dat.gz"));
					// for(int i = 1; i <= (int)grid.size1(); ++i) {
					// 	for(int j = 1; j <= (int)grid.size2(); ++j) {
					// 		for(int k = 1; k <= (int)grid.size3(); ++k) {
					// 			o << grid(i,j,k) << endl;
					// 		}
					// 	}
					// }
					// o.close();
					// utility_exit_with_message("test_grid");
				// }
			}
			dock_get_geom(h.icmp2,h.icmp1,h.iori,d,dcmp2,dcmp1,rscore,cmp2cbc,cmp1cbc);
			Real da1=0,da2=0;
			get_best_sub1_contact_delta_rotations(h.icmp1,dcmp1,h.icmp2,dcmp2,da1,da2);
			Xform x1(rotation_matrix_degrees(arch().axis1(),(Real)h.icmp1),dcmp1*arch().move1());
			Xform x2(rotation_matrix_degrees(arch().axis2(),(Real)h.icmp2),dcmp2*arch().move2());

			int nresiface,nresiface0,nresiface1,nresiface2;
			int ncontacts,ncontacts0,ncontacts1,ncontacts2;
			Real nweighted,nweighted0,nweighted1,nweighted2;
			get_contacts_by_component(h.icmp1,h.icmp2,dcmp1,dcmp2,
				nresiface,nresiface0,nresiface1,nresiface2,
				ncontacts,ncontacts0,ncontacts1,ncontacts2,
				nweighted,nweighted0,nweighted1,nweighted2
			);
			vector1<Vec> pts;
			vector1<Vec> cb1 = cmp1cbs_; std::transform(cb1.begin(),cb1.end(),cb1.begin(),x1);
			vector1<Vec> cb2 = cmp2cbs_; std::transform(cb2.begin(),cb2.end(),cb2.begin(),x2);
			for(Size i = 1; i <= cb1.size()/3; ++i) pts.push_back(cb1[i]);
			for(Size i = 1; i <= cb2.size()/3; ++i) pts.push_back(cb2[i]);
			Real coverage = 0;
			if(arch().name()[0]=='I') coverage = sphere_coverage(pts,SYMICS);
			if(arch().name()[0]=='O') coverage = sphere_coverage(pts,SYMOCT);
			if(arch().name()[0]=='T') coverage = sphere_coverage(pts,SYMTET);
			Real moment1,moment2,moment3,moment4; get_moments(pts,moment1,moment2,moment3,moment4);

			string fn = cmp1name_+"_"+cmp2name_+"_"+arch().name()+(h.reverse?"R":"F")+"_"+ObjexxFCL::string_of(ilm);
			cout << "| " ;
			cout << I(3,nout+1) << " ";
			cout << ObjexxFCL::format::LJ(npad+7,fn) << " ";
			cout << F(7,3,rscore) << " " << F(7,4,rscore/ncontacts) <<" "<< F(7,4,rscore/nweighted) << " ";
            cout << F(7,3, rscore ) <<" "<< I(4,nresiface )<<" "<<I(4,ncontacts )<<" "<<F(7,3,nweighted )<<" ";
            Real tmp = option[tcdock::intra]()*(option[tcdock::intra1]()*cmp1cbc + option[tcdock::intra2]()*cmp2cbc);
            cout << F(7,3,   rscore  - tmp      ) <<" "<< I(4,nresiface0)<<" "<<I(4,ncontacts0)<<" "<<F(7,3,nweighted0)<<" ";
			cout << F(6,2,   cmp1cbc            ) <<" "<< I(4,nresiface1)<<" "<<I(4,ncontacts1)<<" "<<F(7,3,nweighted1)<<" ";
			cout << F(6,2,   cmp2cbc            ) <<" "<< I(4,nresiface2)<<" "<<I(4,ncontacts2)<<" "<<F(7,3,nweighted2)<<" ";
			cout << F(6,2,2*max(fabs(dcmp1)+(dcmp1>0?cmp1diapos_:cmp1dianeg_),fabs(dcmp2)+(dcmp2>0?cmp2diapos_:cmp2dianeg_))) << " "
			     << F(7,5,coverage)<<" "
			     << F(6,2, min_termini_dis(dcmp1,h.icmp1,dcmp2,h.icmp2) ) << " "
			     << I(4,cmp1olig_.size()) << " "
			     << I(3,h.icmp1+(int)da1) << " "
			     << F(8,3,dcmp1) << " "
			     << I(4,cmp2olig_.size()) << " "
			     << I(3,h.icmp2+(int)da2) << " "
			     << F(8,3,dcmp2) << " "
				 << I(3,h.iori);
			cout << F(7,2,moment1) << " " << F(7,2,moment2) << " " << F(7,2,moment3) << " " << F(7,2,moment4);
			if(lnscore_){
				cout << "   " << F(6,2,lnscore_->score(Xforms(1,x1),Xforms(1,x2))) << "  ";
			}
			for(int ifh = 1; ifh <= (int)ffhist.size(); ++ifh) {
				cout << " " << F(5,2,ffhist[ifh]);
			}
			cout << " ";
			if(mhscore_) mhscore_->show(cout,Xforms(1,x1),Xforms(1,x2));
			final_sfxn_->show(cout,x1,x2);

			cout << " " << I(3,nout+1);
			cout << endl;
			vector1<Size> const & dumplist(option[tcdock::dump_pdbs]());
			bool dumpthis = (std::find(dumplist.begin(),dumplist.end(),Size(nout+1)) != dumplist.end() );
			if(option[tcdock::dump_pdb_primary_sub]() > nout ) dump_pdb(h.icmp2+(int)da2,h.icmp1+(int)da1,h.iori,fn+"_sub1.pdb"+(option[tcdock::dump_gz]()?".gz":""), h.reverse, false );
			if(dumpthis || option[tcdock::dump_pdb]() > nout ) dump_pdb(h.icmp2+(int)da2,h.icmp1+(int)da1,h.iori,fn+"_full.pdb"+(option[tcdock::dump_gz]()?".gz":""), h.reverse, true  );
			#ifdef USE_OPENMP
			#pragma omp critical
			#endif
			{
				nout++;
				Real x = cmp1scale_*(Real)h.icmp1;
				Real y = cmp2scale_*(Real)h.icmp2;
				Real z = (cmp1scale_+cmp2scale_)/2.0*h.iori;
				// this is kinda dumb... to handle periodicity
				dumpedit.push_back(Vec(x                 ,y                 ,z));
				dumpedit.push_back(Vec(x+arch().nangle2(),y                 ,z));
				dumpedit.push_back(Vec(x                 ,y+arch().nangle2(),z));
				dumpedit.push_back(Vec(x+arch().nangle2(),y+arch().nangle2(),z));
				dumpedit.push_back(Vec(x-arch().nangle2(),y                 ,z));
				dumpedit.push_back(Vec(x                 ,y-arch().nangle2(),z));
				dumpedit.push_back(Vec(x-arch().nangle2(),y-arch().nangle2(),z));
				dumpedit.push_back(Vec(x+arch().nangle2(),y-arch().nangle2(),z));
				dumpedit.push_back(Vec(x-arch().nangle2(),y+arch().nangle2(),z));
			}

		}
		cout << "DONE symdock_enum_3_1 " << cmp1name_+" "+cmp2name_+" "+arch().name() << endl;
	 }
};


int main (int argc, char *argv[]) {
	try{

	register_options();
	devel::init(argc,argv);
	initsyms();
	using basic::options::option;
	using namespace basic::options::OptionKeys;

	string archname = option[tcdock::architecture]();
	vector1<string> comp1pdbs(option[tcdock::component1]()), comp2pdbs(option[tcdock::component2]());
	ArchCOP arch;
	     if( "I53"  ==archname ) arch = ArchCOP( ArchOP( new ArchI53 ) );
	else if( "I52"  ==archname ) arch = ArchCOP( ArchOP( new ArchI52 ) );
	else if( "I32"  ==archname ) arch = ArchCOP( ArchOP( new ArchI32 ) );
	else if( "O43"  ==archname ) arch = ArchCOP( ArchOP( new ArchO43 ) );
	else if( "O42"  ==archname ) arch = ArchCOP( ArchOP( new ArchO42 ) );
	else if( "O32"  ==archname ) arch = ArchCOP( ArchOP( new ArchO32 ) );
	else if( "T32"  ==archname ) arch = ArchCOP( ArchOP( new ArchT32 ) );
	else if( "T33"  ==archname ) arch = ArchCOP( ArchOP( new ArchT33 ) );
	else if( "P4_42"==archname ) arch = ArchCOP( ArchOP( new ArchP4_42 ) );
	else if( "P4m_42d"==archname ) arch = ArchCOP( ArchOP( new ArchP4m_42d ) );
	else if( "P4g_42d"==archname ) arch = ArchCOP( ArchOP( new ArchP4g_42d ) );
	else if( "P6_32"==archname ) arch = ArchCOP( ArchOP( new ArchP6_32 ) );
	else if( "P6m_32d"==archname ) arch = ArchCOP( ArchOP( new ArchP6m_32d ) );
	else utility_exit_with_message("don't know what architecture you're asking for!");

	TR << "tcdock: " << archname << ", tcdock.cc build date: " << __DATE__ << " " <<  __TIME__ << endl;


	if(option[tcdock::fast_stage_one]() && "T33"==archname ){
		TR << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		TR << "fast_stage_one seems to cause problems with T33, TELL WILL TO FIX IT! (ignoring it for now)" << endl;
		TR << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		option[tcdock::fast_stage_one](false);
		cerr << "WILL TODO: fix dcmp1/dcmp2 sign issue with T33" << endl;
	}
	if(arch->dihedral()){
		TR.Warning << "doing dihedral, disabling fast_stage_one!" << endl;
		option[tcdock::fast_stage_one](false);
	}
	if( "P6_32"==archname && option[tcdock::intra]()>0 && ( option[tcdock::intra1]()>0 || option[tcdock::intra2]()>0 ) ){
		utility_exit_with_message("2D wallpapers don't work right with intra contact scoring yet, intra/intra1/intra2 should be 0");
	}
	// if( "P6_32"==archname && option[tcdock::termini_weight]()>0 ){
	// 	utility_exit_with_message("2D wallpapers don't work right with termini dis yet, -termini_weight should be 0");
	// }
	if(!utility::file::create_directory_recursive(option[out::file::o]())){
		utility_exit_with_message("can't make out::file::o dir"+string(option[out::file::o]()));
	}

	// Real ttrim_cut = option[tcdock::trim_floppy_termini]();
	for(Size i = 1; i <= comp1pdbs.size(); ++i) {
		Pose pose1;
		core::import_pose::pose_from_file(pose1,comp1pdbs[i], core::import_pose::PDB_file);
		// if(ttrim_cut > 0.001){
		// 	TR << "triming tails on " << comp1pdbs[i] << std::endl;
		// 	protocols::sic_dock::auto_trim_floppy_termini(pose1,ttrim_cut,compnfold[1]);
		// }
		if( pose1.size() > (Size)option[tcdock::max_res]() ){
			TR << "skip " << comp1pdbs[i] << " " << pose1.size() << std::endl;
			continue;
		}
		for(Size j = 1; j <= comp2pdbs.size(); ++j) {
			Pose pose2;
			core::import_pose::pose_from_file(pose2,comp2pdbs[j], core::import_pose::PDB_file);
			// if(ttrim_cut > 0.001){
			// 	TR << "triming tails on " << comp2pdbs[j] << std::endl;
			// 	protocols::sic_dock::auto_trim_floppy_termini(pose2,ttrim_cut,compnfold[2]);
			// }
			if( pose2.size() > (Size)option[tcdock::max_res]() ){
				TR << "skip " << comp2pdbs[j] << " " << pose2.size() << std::endl;
				continue;
			}
			TCDock tcd(arch,pose1,pose2,comp1pdbs[i],comp2pdbs[j]);
			try{
				cout << "RUN " << i << " " << j << endl;
				// cout << "RUN " << compkind[1] << compkind[2][1] << " " << arch->name() <<" "<< i << " " << j << endl;
				tcd.run();
			} catch(int e) {
				continue;
			}
		}
	}
	cout << "DONE symdock_enum_3_1" << endl;
    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
    return 0;
}


