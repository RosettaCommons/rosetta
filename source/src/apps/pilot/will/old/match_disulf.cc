// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /src/apps/pilat/will/spiro.cc
/// @brief design around the proposed spiro photocatalyst for h2 production

#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/izstream.hh>
#include <sstream>
#include <numeric/model_quality/rms.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <iostream>
#include <core/kinematics/FoldTree.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <protocols/relax/FastMultiRelax.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


class DisulfEpos : public utility::pointer::ReferenceCount {
public:
	DisulfEpos() {}
	core::Real rms( ObjexxFCL::FArray2D<core::Real> & farrayin ) {
		return numeric::min(numeric::model_quality::rms_wrapper(6,coords,farrayin),
		                    numeric::model_quality::rms_wrapper(6,coords2,farrayin));
	}
	core::Real rms( core::conformation::Residue const & res1, core::conformation::Residue const & res2) {
		ObjexxFCL::FArray2D<core::Real> f = make_farray(res1.xyz(2),res1.xyz(1),res1.xyz(3),res2.xyz(2),res2.xyz(1),res2.xyz(3));
		return rms(f);
	}
	core::Real rms( core::pose::Pose const & pose, core::Size res1, core::Size res2 ) {
		return rms(pose.residue(res1),pose.residue(res2));
	}
	friend std::istream & operator >> (std::istream & in, DisulfEpos & ds) {
		//./ideal/1FCQA_0001.pdb 13 297 CYS CYS H H 284
		std::string stmp;
		in >> stmp; assert(stmp=="DISULF");
		in >> ds.pdb >> ds.pdbres1 >> ds.pdbres2;
		in >> stmp; assert(stmp=="CYS");
		in >> stmp;	assert(stmp=="CYS");
		in >> ds.ss1 >> ds.ss2 >> ds.sep;
		in >> stmp; assert(stmp=="EPOS");
		in >> ds.ca1 >> ds.n1 >> ds.c1 >> ds.ca2 >> ds.n2 >> ds.c2;
		ds.coords  = ds.make_farray(ds.ca1,ds.n1,ds.c1,ds.ca2,ds.n2,ds.c2);
		ds.coords2 = ds.make_farray(ds.ca2,ds.n2,ds.c2,ds.ca1,ds.n1,ds.c1);
		return in;
	}
	friend std::ostream & operator << (std::ostream & out, DisulfEpos const & ds) {
		out << "DISULF " << ds.pdb << " " << ds.pdbres1 << " " << ds.pdbres2 << " "
		    << "CYS CYS " << ds.ss1 << " " << ds.ss2 << " " << ds.sep << " EPOS "
			 << ds.ca1 << ds.n1 << ds.c1 << ds.ca2 << ds.n2 << ds.c2 << std::endl;
		return out;
	}
	ObjexxFCL::FArray2D<core::Real> make_farray(numeric::xyzVector<core::Real> ca1,numeric::xyzVector<core::Real> n1,numeric::xyzVector<core::Real> c1,numeric::xyzVector<core::Real> ca2,numeric::xyzVector<core::Real> n2,numeric::xyzVector<core::Real> c2) {
		ObjexxFCL::FArray2D<core::Real> f(3,6);
		f(1,1) = ca1.x();	f(2,1) = ca1.y(); f(3,1) = ca1.z();
		f(1,2) =  n1.x();	f(2,2) =  n1.y(); f(3,2) =  n1.z();
		f(1,3) =  c1.x();	f(2,3) =  c1.y(); f(3,3) =  c1.z();
		f(1,4) = ca2.x();	f(2,4) = ca2.y(); f(3,4) = ca2.z();
		f(1,5) =  n2.x();	f(2,5) =  n2.y(); f(3,5) =  n2.z();
		f(1,6) =  c2.x();	f(2,6) =  c2.y(); f(3,6) =  c2.z();
		return f;
	}
// private:
	std::string pdb;
	core::Size pdbres1,pdbres2,sep;
	char ss1,ss2;
	numeric::xyzVector<core::Real> ca1,n1,c1,ca2,n2,c2;
	ObjexxFCL::FArray2D<core::Real> coords,coords2;
};

void addcc(core::pose::Pose & pose, core::id::AtomID aid, core::id::AtomID anchor, core::Real mult = 1.0 ) {
	core::scoring::constraints::ConstraintOP cc = new core::scoring::constraints::CoordinateConstraint(
		aid, anchor, pose.xyz(aid), new core::scoring::constraints::HarmonicFunc(0,
		 mult * basic::options::option[basic::options::OptionKeys::docking::ligand::harmonic_Calphas]()
		) );
	pose.add_constraint(cc);
}

class DisulfEposDatabase : public utility::pointer::ReferenceCount {
public:
	DisulfEposDatabase() {}
	void read_data_file( std::string fname ) {
		utility::io::izstream in(fname);
		DisulfEpos ds;
		char buf[999];
		std::istringstream iss;
		while( in.getline(buf,999) ) {
			iss.str(buf);
			iss >> ds;
			disulfs.push_back(ds);
			// std::cerr << disulfs.size() << " " << ds << std::endl;
		}
		std::cerr << "done " << disulfs.size() << std::endl;
	}
	core::Real rms( core::conformation::Residue const & res1, core::conformation::Residue const & res2) {
		ObjexxFCL::FArray2D<core::Real> f = disulfs[1].make_farray(res1.xyz(2),res1.xyz(1),res1.xyz(3),res2.xyz(2),res2.xyz(1),res2.xyz(3));
		core::Real min = 9e9;
		for(utility::vector1<DisulfEpos>::iterator i = disulfs.begin(); i != disulfs.end(); ++i ) {
			// if( i->ss1!='L' && i->ss2!='L' ) continue;
			core::Real tmp = i->rms(f);
			if(tmp < min) min = tmp;
		}
		return min;
	}
	core::Real rms( core::pose::Pose const & pose, core::Size res1, core::Size res2 ) {
		return rms(pose.residue(res1),pose.residue(res2));
	}
	DisulfEpos & disulf( core::Size i ) { return disulfs[i]; }
private:
	utility::vector1<DisulfEpos> disulfs;
};

int main(int argc, char *argv[]) {

	try {

	using namespace core;
	using namespace chemical;
	using namespace scoring;

	devel::init(argc,argv);

	DisulfEposDatabase ds;
	std::cerr << "reading disulf lib..." ;
	ds.read_data_file("input/DISULF_EPOS.txt");

	pose::Pose pose;
	std::cerr << "reading pdb input/if_x.pdb" << std::endl;
	core::import_pose::pose_from_file(pose,"input/if_x.pdb", core::import_pose::PDB_file);
	using core::conformation::ResidueFactory;
	ResidueTypeSetCAP residue_set( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
	pose.append_residue_by_jump(*ResidueFactory::create_residue(residue_set->name_map("VRT")),1);
	pose.append_residue_by_jump(*ResidueFactory::create_residue(residue_set->name_map("VRT")),140);

	// std::cerr << "searching..." << std::endl;
	// Size mni=0,mnj=0;
	// Real mn = 9e9;
	// for(Size i = 1; i <= 270; ++i) {
	// 	for(Size j = 271; j <= 540; ++j) {
	// 		Real d = pose.residue(i).xyz(2).distance(pose.residue(j).xyz(2));
	// 		if( d > 8.4 ) continue;
	// 		Real r = ds.rms(pose,i,j);
	// 		if( r < 0.2 ) {
	// 			std::cerr << i << " " << j << " " << d << " " << r << std::endl;
	// 		}
	// 		if( r < mn ) {
	// 			mni = i; mnj = j; mn = r;
	// 		}
	// 	}
	// }
	// std::cerr << "min rms " << mn	<< " " << mni << " " << mnj << std::endl;
	std::cout << "=================================== FOLD TREE ================================" << std::endl;
	std::cout << pose.fold_tree() << std::endl;

	// set fold tree?
	ObjexxFCL::FArray2D_int jumps(2,3);
	ObjexxFCL::FArray1D_int cuts(3);
	cuts(1) = 270;
	cuts(2) = 271;
	cuts(3) = 140;
	jumps(1,1) =  271; jumps(2,1) = 272;
	jumps(1,2) =  1  ; jumps(2,2) = 271;
	jumps(1,3) =  141; jumps(2,3) = 272;
	core::kinematics::FoldTree ft = pose.fold_tree();
	ft.tree_from_jumps_and_cuts( pose.n_residue(), 3, jumps, cuts, 271 );
	ft.reorder(271);
	pose.fold_tree(ft);
	std::cerr << "FT " << pose.fold_tree() << std::endl;
	std::cout << "==============================================================================" << std::endl;

	core::kinematics::MoveMapOP mm = new core::kinematics::MoveMap;
	mm->init_from_file( basic::options::option[basic::options::OptionKeys::in::file::movemap]() );
	for( Size i = 1; i <= pose.n_residue(); ++i ) {
		core::id::AtomID anchor(1,271);
		if(i > 140) anchor.rsd() = 272;
		if( ! mm->get_bb(i) ) addcc(pose,core::id::AtomID(2,i),anchor,100.0);
	}
	mm->set_jump(false);
	// mm->set(core::id::RB1,true);
	mm->set_jump(2,false);
	mm->set_jump(3,false);
	// protocols::relax::FastRelax relax(core::scoring::get_score_function(),"");
	protocols::relax::FastMultiRelax relax(core::scoring::get_score_function());
	relax.set_movemap(mm);
	relax.apply(pose);

	// Real mxdis = 0;
	// for(Size i = 1; i <= 3099; ++i ) {
	// 	mxdis = numeric::max(mxdis,ds.disulf(i).ca1.distance(ds.disulf(i).ca2));
	//
	// // 	for(Size j = 1; j < i; ++j) {
	// // 		std::cerr << i << " " << j << " " << ds.disulf(i).rms(ds.disulf(j).coords) << std::endl;
	// // 	}
	// }
	// std::cerr << "mxdis " << mxdis << std::endl;

	pose.dump_pdb(basic::options::option[basic::options::OptionKeys::out::file::o]());

	return 0;


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
