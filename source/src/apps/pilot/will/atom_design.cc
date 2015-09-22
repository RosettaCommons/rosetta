// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/apps/pilat/will/dubois_graft.cc
/// @brief graft peptide structs onto dubois catalyst

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/conformation/Residue.hh>
#include <core/pose/symmetry/util.hh>

#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/MoveMap.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <core/pack/optimizeH.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/moves/Mover.hh>
#include <sstream>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>

#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <protocols/moves/MoverStatistics.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>
#include <apps/pilot/will/will_util.ihh>


using numeric::conversions::radians;

static THREAD_LOCAL basic::Tracer TR( "atom_design" );

using core::Size;
using core::Real;
using core::id::AtomID;
typedef numeric::xyzVector<Real> Vec;
typedef utility::vector1<Vec>    Vecs;
typedef numeric::xyzMatrix<Real> Mat;
using protocols::moves::MoverOP;
using core::scoring::ScoreFunctionOP;
using numeric::random::uniform;


void adesign(core::pose::Pose & pose, ScoreFunctionOP sf, core::chemical::ResidueTypeSet const & rs) {

	std::map<std::string,utility::vector1<std::string> > adtypemap;
	utility::vector1<std::string> tgrp1; tgrp1.push_back("CH5"); tgrp1.push_back("CR5"); tgrp1.push_back("NH5");
	for(utility::vector1<std::string>::iterator i = tgrp1.begin(); i != tgrp1.end(); ++i) adtypemap[*i] = tgrp1;

	using namespace core::pack::task;
	PackerTaskOP task = TaskFactory::create_packer_task(pose);
	utility::vector1< bool > aas(20,false);
	for(Size i = 1; i <= pose.n_residue(); ++i) {
		if( adtypemap.find(pose.residue(i).name()) != adtypemap.end() ) {
			for(Size j = 1; j <= adtypemap.size(); ++j) {
				task->nonconst_residue_task(i).allow_noncanonical_aa(adtypemap[pose.residue(i).name()][j],rs);
			}
		} else {
			// task->nonconst_residue_task(i).prevent_repacking();
			task->nonconst_residue_task(i).restrict_to_repacking();
			task->nonconst_residue_task(i).or_include_current(true);
		}
	}
	TR << *task << std::endl;
	protocols::simple_moves::PackRotamersMover repack( sf, task );
	repack.apply(pose);

}

void mydumppdb(core::pose::Pose & pose,std::string fname) {
	utility::io::ozstream out(fname);
	core::id::AtomID_Map<Real> bfac;
	core::pose::initialize_atomid_map(bfac,pose);
	core::io::pdb::dump_bfactor_pdb(pose,bfac,out);
	out.close();
}

/////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {


	using namespace core;
	using namespace chemical;
	using namespace pose;

	devel::init(argc,argv);

	core::chemical::ResidueTypeSet & rs = core::chemical::ChemicalManager::get_instance()->nonconst_residue_type_set( core::chemical::FA_STANDARD );
	// for(ResidueTypeSet::const_residue_iterator i = rs->all_residues_begin(); i != rs->all_residues_end(); ++i) {
	// 	TR << i->first << " " << i->second->name() << std::endl;
	// }
	// rs.nonconst_name_map("CH5").set_RotamerLibraryName("input/CH5.pdb");
	// rs.nonconst_name_map("NH5").set_RotamerLibraryName("input/NH5.pdb");
	// rs.nonconst_name_map("CR5").set_RotamerLibraryName("input/CR5.pdb");

	ScoreFunctionOP sf = core::scoring::get_score_function();

	Pose pose;
	import_pose::pose_from_pdb(pose,rs,basic::options::option[basic::options::OptionKeys::in::file::s]()[1]);
	pose.dump_pdb("init.pdb");

	// replace_pose_residue_copying_existing_coordinates(pose,1,rs.name_map("NH5"));
	// replace_pose_residue_copying_existing_coordinates(pose,2,rs->name_map("CR5"));

	core::pack::optimizeH(pose,*sf);
	mydumppdb(pose,"pre_design.pdb");
	utility::vector1<std::string> tgrp1; tgrp1.push_back("CH5"); /*tgrp1.push_back("CR5");*/ /*tgrp1.push_back("NH5");*/

	Size start_pos = 0;
	for(Size i = 1; i <= 2; ++i) {
		for(Size j = 1; j <= tgrp1.size(); ++j) {
			replace_pose_residue_copying_existing_coordinates(pose,start_pos+i,rs.name_map(tgrp1[j]));
			// core::pack::optimizeH(pose,*sf);
			sf->score(pose);
			std::cout << tgrp1[j] << " " << std::endl;
			pose.energies().show(std::cout,start_pos+i);
		}
	}
	// utility_exit_with_message("DBEUG");

	for(Size i = 1; i <= 2; ++i) {
		sf->score(pose);
		pose.energies().show(std::cout,start_pos+i);
	}


	sf->show(pose);
	adesign(pose,sf,rs);
	sf->show(pose);

	mydumppdb(pose,"post_design.pdb");

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
