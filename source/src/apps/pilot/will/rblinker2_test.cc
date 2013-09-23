// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

#include <basic/basic.hh>
#include <basic/database/open.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/keys/rblinker.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/simple_moves/BBGaussianMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/scoring/methods/ImplicitClashEnergy.hh>
#include <protocols/viewer/viewers.hh>
#include <sstream>
#include <set>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <apps/pilot/will/mynamespaces.ihh>

using core::kinematics::Stub;
using core::conformation::Residue;
using core::Real;

static basic::Tracer TR("rblinker2");

// static numeric::random::RandomGenerator RG(8334046);
// 
// inline void xform_pose( core::pose::Pose & pose, Stub const & s ) {
// 	for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
// 		for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
// 			core::id::AtomID const aid(core::id::AtomID(ia,ir));
// 			pose.set_xyz( aid, s.local2global(pose.xyz(aid)) );
// 		}
// 	}
// }
// 
// 
// 
// // petf   37   42   45   75
// // hyd   101  156  334  338
// // psI  1492 1536 1489 1495
// inline Real hyd_petf_sf4_dis(Pose const & hyd, Stub const & shyd, Pose const & petf, Stub const & spetf ) {
// 	Vec sf4A = ( hyd.xyz(AtomID(5, 101))+ hyd.xyz(AtomID(5, 156))+ hyd.xyz(AtomID(5, 334))+ hyd.xyz(AtomID(5, 338)))/4.0;
// 	Vec sf4B = (petf.xyz(AtomID(5,  37))+petf.xyz(AtomID(5,  42))+petf.xyz(AtomID(5,  45))+petf.xyz(AtomID(5,  75)))/4.0;
// 	return shyd.local2global(sf4A).distance(spetf.local2global(sf4B));
// }
// inline Real psI_petf_sf4_dis(Pose const & psI, Pose const & petf, Stub const & spetf ) {
// 	Vec sf4A = ( psI.xyz(AtomID(5,1492))+ psI.xyz(AtomID(5,1536))+ psI.xyz(AtomID(5,1489))+ psI.xyz(AtomID(5,1495)))/4.0;
// 	Vec sf4B = (petf.xyz(AtomID(5,  37))+petf.xyz(AtomID(5,  42))+petf.xyz(AtomID(5,  45))+petf.xyz(AtomID(5,  75)))/4.0;
// 	return sf4A.distance(spetf.local2global(sf4B));
// }
// 
// // get stup that aligns r1 to r2
// Stub getxform(Residue const & r1, Residue const & r2) {
// 	Stub s;
// 	s.M = alignVectorSets(r1.xyz(1)-r1.xyz(2),r1.xyz(3)-r1.xyz(2),r2.xyz(1)-r2.xyz(2),r2.xyz(3)-r2.xyz(2));
// 	s.v = r2.xyz(2)-s.M*r1.xyz(2);
// 	return s;
// }



inline void xform_pose( core::pose::Pose & pose, core::kinematics::Stub const & s, Size sres=1, Size eres=0 ) {
	if(eres==0) eres = pose.n_residue();
	for(Size ir = sres; ir <= eres; ++ir) {
		for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
			core::id::AtomID const aid(core::id::AtomID(ia,ir));
			pose.set_xyz( aid, s.local2global(pose.xyz(aid)) );
		}
	}
}
inline void xform_pose_rev( core::pose::Pose & pose, core::kinematics::Stub const & s ) {
	for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
		for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
			core::id::AtomID const aid(core::id::AtomID(ia,ir));
			pose.set_xyz( aid, s.global2local(pose.xyz(aid)) );
		}
	}
}




class SimpleBBMover : public protocols::moves::Mover {
	Size start_,stop_;
	Real mag_;
public:
	SimpleBBMover(Size start, Size stop, Real mag) : start_(start),stop_(stop),mag_(mag) {}
	Real magnitude(        ) { return mag_; }
	void magnitude(Real mag) { mag_ = mag; }	
	void apply(core::pose::Pose & pose) {
		Size i = start_-1 + std::ceil(uniform()*(stop_-start_+1));
		if(uniform()<0.5) pose.set_phi(i,pose.phi(i)+gaussian()*mag_);
		else              pose.set_psi(i,pose.psi(i)+gaussian()*mag_);
	}
	std::string get_name() const { return "SimpleBBMover"; }
};
typedef utility::pointer::owning_ptr<SimpleBBMover> SimpleBBMoverOP;

std::string printbits(unsigned long nn) {
	std::string s = "";
	for(Size i = 0; i < 8*sizeof(nn); ++i) {
		if(i%8==0) s += " ";
		if(nn & 1UL << (8*sizeof(nn)-1-i)) s += "1";
		else                               s += "0";
	}
	return s;
}

std::string bin2string(unsigned long bin, Size nres) {
	std::string s = "";
	for(Size i = 0; i < nres; ++i) {
		for(Size j = 0; j < 2; ++j) {
			int tmp = ( bin >> 8*i+2*j ) % 4;
			s += ObjexxFCL::lead_zero_string_of(tmp,1);
		}
		if(i+1 < nres) s += "-";
	}
	return s;
}

unsigned long pose2bin(core::pose::Pose const & pose) {
	using namespace ObjexxFCL::format;
	unsigned long bin = 0;
	for(int i = 0; i < min( 16, (int)pose.n_residue() ); ++i) {		
		Real phid = pose.phi(i+1);
		Real psid = pose.psi(i+1);
		// TR << phid << " " << psid << std::endl;
		unsigned long phi = basic::unsigned_periodic_range(phid,360.0) / 90.0;
		unsigned long psi = basic::unsigned_periodic_range(psid,360.0) / 90.0;
		phi = phi << 4*i;
		psi = psi << 4*i+2;
		bin += phi;
		bin += psi;
	}
	return bin;
}

Size get_aln_resi(core::pose::Pose const & pose, core::Size const aln_chain, char const aln_termi) {
	Size resi = 0, nnt = 0, nct = 0;
	for(Size i = 1; i <= pose.n_residue(); ++i) {
		if(pose.residue(i).is_lower_terminus()) {
			nnt++;
			if(aln_termi=='N' && nnt==aln_chain) {
				resi = i;
				break;
			}
		}
		if(pose.residue(i).is_upper_terminus()) {
			nct++;
			if(aln_termi=='C' && nnt==aln_chain) {
				resi = i;
				break;
			}
		}
	}
	if( 0 == resi ) utility_exit_with_message("ImplicitClashData: could not find aln resi in input pose chain: "+ObjexxFCL::string_of(aln_chain)+" term: "+aln_termi);
	return resi;
}
core::kinematics::Stub getxform(core::conformation::Residue const & move_resi, core::conformation::Residue const & fixd_resi) {
	core::kinematics::Stub s;
	s.M = alignVectorSets(move_resi.xyz(1)-move_resi.xyz(2),move_resi.xyz(3)-move_resi.xyz(2),fixd_resi.xyz(1)-fixd_resi.xyz(2),fixd_resi.xyz(3)-fixd_resi.xyz(2));
	s.v = fixd_resi.xyz(2)-s.M*move_resi.xyz(2);
	return s;
}
core::kinematics::Stub getxform(core::conformation::ResidueCOP      move_resi, core::conformation::Residue const & fixd_resi) { return getxform(*move_resi, fixd_resi); }
core::kinematics::Stub getxform(core::conformation::Residue const & move_resi, core::conformation::ResidueCOP      fixd_resi) { return getxform( move_resi,*fixd_resi); }
core::kinematics::Stub getxform(core::conformation::ResidueCOP      move_resi, core::conformation::ResidueCOP      fixd_resi) { return getxform(*move_resi,*fixd_resi); }


core::pose::Pose build_algned_linker(core::pose::Pose const & alnpose, Size len, core::chemical::ResidueTypeSetCAP resset) {
	Pose lnk;
	string ggs = "GGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGS";
	make_pose_from_sequence( lnk , "A"+ggs.substr(0,len)+"A" , *resset , true);
	for(Size i = 1; i <= lnk.n_residue(); ++i) {
		lnk.set_phi  (i,-60);
		lnk.set_psi  (i,-45);
		lnk.set_omega(i,180);
	}
	FoldTree ft = lnk.fold_tree();
	if(ft.num_jump() != 0) utility_exit_with_message("why is there a jump in the fold tree alreaduy?!?!?!");
	ObjexxFCL::FArray2D_int jump_point(2,1);
	ObjexxFCL::FArray1D_int cuts(1);
	cuts(1) = lnk.n_residue()/2;
	core::pose::add_variant_type_to_pose_residue(lnk,"CUTPOINT_LOWER",cuts(1)+0);		
	core::pose::add_variant_type_to_pose_residue(lnk,"CUTPOINT_UPPER",cuts(1)+1);
	jump_point(1,1) = 1;
	jump_point(2,1) = lnk.n_residue();
	ft.tree_from_jumps_and_cuts( lnk.n_residue(), 1, jump_point, cuts );
	ft.set_jump_atoms(1,"N","C");
	lnk.fold_tree(ft);		

	xform_pose(lnk, getxform(lnk.residue(     1         ),alnpose.residue(1)) );
	xform_pose(lnk, getxform(lnk.residue(lnk.n_residue()),alnpose.residue(2)), cuts(1)+1 ); // only after cutpoint
	
	return lnk;
}


void* doit(void*) {
	using namespace core;
	using namespace chemical;
	using namespace conformation;
	using namespace pose;
	using namespace protocols;
	using namespace moves;
	using namespace ObjexxFCL::format;
	using numeric::random::uniform;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	ResidueTypeSetCAP cenresset = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	Pose lnk;
	Size linklen1 = option[rblinker::linker1_size]();
	Size linklen2 = option[rblinker::linker2_size]();
	string seq1, seq2;
	string ggs = "GGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGS";
	string lll = "***************************************************************************************";
	seq1 = "A" + ggs.substr(0,linklen1) + "A";
	seq2 = "A" + ggs.substr(0,linklen2) + "A";
	string ss1 = "_" + lll.substr(0,linklen1) + "_";
	string ss2 = "_" + lll.substr(0,linklen2) + "_";

	TR << "linker1 size " << linklen1 << " seq1: " << seq1 << std::endl;
	TR << "linker2 size " << linklen2 << " seq2: " << seq2 << std::endl;

	make_pose_from_sequence(lnk,seq1,*cenresset,true);
	for(Size i = 1; i <= lnk.n_residue(); ++i) {
		lnk.set_phi  (i,-60);
		lnk.set_psi  (i,-45);
		lnk.set_omega(i,180);
	}

	core::scoring::ScoreFunctionOP sf = new core::scoring::ScoreFunction;
	sf->set_weight(core::scoring::vdw  ,5.0);
	sf->set_weight(core::scoring::pair ,1.0);
	sf->set_weight(core::scoring::cbeta,1.0);
	sf->set_weight(core::scoring::rama ,1.0);
	sf = core::scoring::getScoreFunction();
		
	Pose test,hyd;
	core::import_pose::pose_from_pdb(test,*cenresset,"input/psI_0001_strip_0001.pdb");
	core::import_pose::pose_from_pdb(hyd ,*cenresset,"input/hyda1_0001_strip_0001.pdb");
	Size st = test.n_residue();
	
	core::pose::remove_upper_terminus_type_from_pose_residue(test,test.n_residue());
	core::pose::remove_lower_terminus_type_from_pose_residue(lnk,1);
	for(Size i = 1; i <= lnk.n_residue(); ++i) test.append_residue_by_bond(lnk.residue(i),true);
	core::pose::remove_upper_terminus_type_from_pose_residue(test,test.n_residue());
	core::pose::remove_lower_terminus_type_from_pose_residue(hyd,1);
	for(Size i = 1; i <= hyd.n_residue(); ++i) test.append_residue_by_bond(hyd.residue(i),true);
	
	test = lnk;
	
	// MoverOP bbmove = new SimpleBBMover(st+2,st+lnk.n_residue(),30.0); 	
	MoverOP bbmove = new SimpleBBMover(1,lnk.n_residue(),30.0); 	
	MonteCarloOP mc1 = new MonteCarlo( test, *sf, 2.0 );
	mc1->set_autotemp( true, 2.0 ); mc1->set_temperature( 2.0 );
	TrialMover trial(bbmove,mc1);

	for(int j = 1; j <= option[out::nstruct](); j++) {
		// TR << "trial " << j << std::endl;
		trial.apply(test);
	}
	test.dump_pdb("test.pdb");
	TR << "done" << std::endl;

	return NULL;
}



int main( int argc, char * argv [] ) {

	try {

	devel::init(argc,argv);

	void* (*func)(void*) = &doit;
	if (option[ basic::options::OptionKeys::parser::view ]()) {
		protocols::viewer::viewer_main( func );
	} else {
		func(NULL);
	}


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

}

/*

lnk
1560 runs OK
1691 runs OK
1692 runs OK
1756 runs OK

lnk2
1560 runs OK
1691 runs OK
1692 runs OK
1756 runs OK


*/
