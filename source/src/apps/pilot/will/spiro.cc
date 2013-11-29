// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file /src/apps/pilat/will/spiro.cc
/// @brief design around the proposed spiro photocatalyst for h2 production

//moved to top to circumvent compiler errors on some version of GCC - these need to be after core/init/init.hh if using namespace core is in scope.
#include <utility/io/izstream.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>


#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/VariantType.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymmData.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>

// AUTO-REMOVED #include <core/fragment/FragmentIO.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/FragSet.hh>
#include <devel/init.hh>
#include <basic/database/open.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/util.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
// AUTO-REMOVED #include <core/scoring/constraints/AmbiguousConstraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/AngleConstraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/AtomPairConstraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/CoordinateConstraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
// AUTO-REMOVED #include <core/scoring/constraints/util.hh>
#include <core/scoring/Energies.hh>
// AUTO-REMOVED #include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <basic/Tracer.hh>
// AUTO-REMOVED #include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
// AUTO-REMOVED #include <numeric/xyz.functions.hh>
// AUTO-REMOVED #include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/abinitio/FragmentMover.hh>
#include <protocols/flxbb/FlxbbDesign.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
// AUTO-REMOVED #include <protocols/symmetric_docking/SymDockingLowRes.hh>
#include <sstream>

//Auto Headers
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <protocols/toolbox/SwitchResidueTypeSet.hh>
#include <utility/io/mpistream.hh>
#include <ObjexxFCL/format.hh>

//Auto using namespaces
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end


using core::Size;
using core::Real;
using numeric::random::uniform;
using core::id::AtomID;
using core::id::DOF_ID;
using core::scoring::ScoreFunctionOP;

static basic::Tracer TR("spiro");

struct PoseWrap {
	PoseWrap() : hascst(false) {}
	core::pose::Pose pose;
	Size nsub,attach,nres;
	bool hascst;
	void set_phi(Size seqpos, Real setting) {
		if(seqpos==attach) {
			pose.set_chi(3,1,setting);
		} else {
			pose.set_phi(seqpos,setting);
		}
	}
};

void switch_to_fa(PoseWrap & pw) {
	protocols::toolbox::switch_to_residue_type_set( pw.pose, core::chemical::FA_STANDARD );
}

void switch_to_cen(PoseWrap & pw) {
	protocols::toolbox::switch_to_residue_type_set( pw.pose, core::chemical::CENTROID );
}

void minimize(PoseWrap & pw, ScoreFunctionOP sf) {
	core::pose::Pose & pose(pw.pose);
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	movemap->set_jump(false);
	movemap->set_chi(true);
	movemap->set_bb(true);
	// movemap->set_bb(1,false);
	protocols::simple_moves::symmetry::SymMinMover m( movemap, sf, "dfpmin_armijo_nonmonotone", 1e-4, true );
	m.apply(pose);
}


void read_fragdata( utility::vector1< core::fragment::FragDataOP > & fds, utility::io::izstream & in, bool /*design = false*/ ) {
	using namespace core::fragment;
	Size n,count=0;
	while( in >> n ) {
	 	std::string pdb;
		char buf[999];
		FragDataOP fd = new FragData;
		for( Size i = 1; i <= n; ++i ) {
			utility::pointer::owning_ptr<SingleResidueFragData> srfd;
			srfd = new BBTorsionSRFD;
			in >> pdb;
			in.getline(buf,999);
			std::istringstream iss(buf);
			iss >> *srfd;
			fd->add_residue(srfd);
		}
		fd->set_valid(true);
		fds.push_back(fd);
		count++;
		// if( count >= nread ) break;
	}
	in.close();
}

std::map<std::string, utility::vector1<core::fragment::FragDataOP> >
get_frags_map( bool design = false ) {
	using namespace core::fragment;
	TR << "reading frags" << std::endl;
	utility::io::izstream in;
	std::map<std::string,utility::vector1<FragDataOP> > fds;
	// basic::database::open(in,"sampling/ss_fragfiles/EEE.fragfile"); read_fragdata(fds["EEE"],in,design);
	// basic::database::open(in,"sampling/ss_fragfiles/EEH.fragfile"); read_fragdata(fds["EEH"],in,design);
	// basic::database::open(in,"sampling/ss_fragfiles/EEL.fragfile"); read_fragdata(fds["EEL"],in,design);
	// basic::database::open(in,"sampling/ss_fragfiles/EHH.fragfile"); read_fragdata(fds["EHH"],in,design);
	// basic::database::open(in,"sampling/ss_fragfiles/ELE.fragfile"); read_fragdata(fds["ELE"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/ELH.fragfile"); read_fragdata(fds["ELH"],in,design);
	// basic::database::open(in,"sampling/ss_fragfiles/ELL.fragfile"); read_fragdata(fds["ELL"],in,design);
	// basic::database::open(in,"sampling/ss_fragfiles/HEE.fragfile"); read_fragdata(fds["HEE"],in,design);
	// basic::database::open(in,"sampling/ss_fragfiles/HEH.fragfile"); read_fragdata(fds["HEH"],in,design);
	// basic::database::open(in,"sampling/ss_fragfiles/HEL.fragfile"); read_fragdata(fds["HEL"],in,design);
	// basic::database::open(in,"sampling/ss_fragfiles/HHE.fragfile"); read_fragdata(fds["HHE"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/HHH.fragfile"); read_fragdata(fds["HHH"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/HHL.fragfile"); read_fragdata(fds["HHL"],in,design);
	// basic::database::open(in,"sampling/ss_fragfiles/HLE.fragfile"); read_fragdata(fds["HLE"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/HLH.fragfile"); read_fragdata(fds["HLH"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/HLL.fragfile"); read_fragdata(fds["HLL"],in,design);
	// basic::database::open(in,"sampling/ss_fragfiles/LEE.fragfile"); read_fragdata(fds["LEE"],in,design);
	// basic::database::open(in,"sampling/ss_fragfiles/LEH.fragfile"); read_fragdata(fds["LEH"],in,design);
	// basic::database::open(in,"sampling/ss_fragfiles/LEL.fragfile"); read_fragdata(fds["LEL"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/LHH.fragfile"); read_fragdata(fds["LHH"],in,design);
	// basic::database::open(in,"sampling/ss_fragfiles/LLE.fragfile"); read_fragdata(fds["LLE"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/LLH.fragfile"); read_fragdata(fds["LLH"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/LLL.fragfile"); read_fragdata(fds["LLL"],in,design);
	return fds;
}


core::fragment::FragSetOP make_frag_set(Size start, std::string ss, std::map<std::string, utility::vector1<core::fragment::FragDataOP> > fds) {
	using namespace core::fragment;
	FragSetOP frags = new ConstantLengthFragSet();
	Size const stop = ss.size() + start - 3;
	if(start >= stop) return NULL;
	for( Size i = start; i <= stop; ++i ) {
		FrameOP frame = new Frame(i,3);
		utility::vector1<FragDataOP>::iterator beg = fds[ss.substr(i-start,3)].begin();
		utility::vector1<FragDataOP>::iterator end = fds[ss.substr(i-start,3)].end();
		for( utility::vector1<FragDataOP>::iterator fi = beg; fi != end; ++fi ) {
			frame->add_fragment(*fi);
		}
		frags->add(frame);
	}
	return frags;
}


PoseWrap make_pose(std::string seq) {
	using namespace core;
	using namespace chemical;
	using namespace conformation;
	PoseWrap pw;
	pw.attach = 2;
	pw.nsub = 2;
	pw.nres = seq.size();
	ResidueTypeSetCAP residue_set( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
	core::import_pose::pose_from_pdb(pw.pose,"input/start_gly.pdb" ,residue_set);
	// pw.pose.dump_pdb("init0.pdb");
	// core::pose::remove_lower_terminus_type_from_pose_residue(pw.pose,2);
	core::pose::remove_upper_terminus_type_from_pose_residue(pw.pose,pw.pose.n_residue());
	// pw.pose.dump_pdb("init_noterm.pdb");
	// std::string name3 = name_from_aa(aa_from_oneletter_code(seq[0]));
	// pw.pose.replace_residue(2,*ResidueFactory::create_residue(residue_set->name_map(name3)),true);
	for( Size i = 3; i <= pw.nres; ++i ) {
		std::string name3 = name_from_aa(aa_from_oneletter_code(seq[i-2]));
		pw.pose.append_residue_by_bond(*ResidueFactory::create_residue(residue_set->name_map(name3)),true);
	}
	core::pose::add_upper_terminus_type_to_pose_residue(pw.pose,pw.pose.n_residue());
	// pw.pose.dump_pdb("init_allres.pdb");
	kinematics::FoldTree ft = pw.pose.fold_tree();
	ft.jump_edge(1).start() = 1;
	ft.jump_edge(1).stop () = pw.attach;
	ft.jump_edge(1).start_atom() = "VN";
	ft.jump_edge(1).stop_atom()  =  "N";
	// ft.reorder(1);
	pw.pose.fold_tree(ft);
	// std::cout << ft << std::endl;

	core::pose::symmetry::make_symmetric_pose( pw.pose );
	// pw.pose.dump_pdb("make_pose_symm.pdb");
	switch_to_cen(pw);
	// pw.pose.dump_pdb("make_pose_symm_cen.pdb");
	ft = pw.pose.fold_tree();
	for( Size i = 1; i <= pw.nsub; ++i ) {
		ft.jump_edge(i).start_atom() = "VN";
		ft.jump_edge(i).stop_atom()  =  "N";
	}
	pw.pose.fold_tree(ft);
	for( Size i = 2; i <= pw.nres; ++i ) {
		pw.pose.set_phi  (i,-60.16731);
		pw.pose.set_psi  (i,-45.19451);
		pw.pose.set_omega(i,180.00000);
	}
	// pw.pose.dump_pdb("test1.pdb");
	// pw.pose.set_chi(1,1,pw.pose.chi(1,1)+10);
	// pw.pose.set_chi(2,1,pw.pose.chi(2,1)+10);
	// pw.pose.set_chi(3,1,pw.pose.chi(3,1)+10);
	// pw.pose.dump_pdb("test2.pdb");
	// pw.pose.set_chi(1,1,pw.pose.chi(1,1)+10);
	// pw.pose.set_chi(2,1,pw.pose.chi(2,1)+10);
	// pw.pose.set_chi(3,1,pw.pose.chi(3,1)+10);
	// pw.pose.dump_pdb("test3.pdb");
	// pw.pose.set_chi(1,1,pw.pose.chi(1,1)+10);
	// pw.pose.set_chi(2,1,pw.pose.chi(2,1)+10);
	// pw.pose.set_chi(3,1,pw.pose.chi(3,1)+10);
	// pw.pose.dump_pdb("test4.pdb");
	// pw.pose.set_chi(1,1,pw.pose.chi(1,1)+10);
	// pw.pose.set_chi(2,1,pw.pose.chi(2,1)+10);
	// pw.pose.set_chi(3,1,pw.pose.chi(3,1)+10);
	// pw.pose.dump_pdb("test5.pdb");
	// pw.pose.set_chi(1,1,pw.pose.chi(1,1)+10);
	// pw.pose.set_chi(2,1,pw.pose.chi(2,1)+10);
	// pw.pose.set_chi(3,1,pw.pose.chi(3,1)+10);
	// pw.pose.dump_pdb("test6.pdb");
	return pw;
}

void repack(PoseWrap & pw, ScoreFunctionOP sf) {
	core::pose::Pose & pose(pw.pose);
	using namespace core::pack::task;
	PackerTaskOP task = TaskFactory::create_packer_task(pose);
	task->restrict_to_repacking();
	task->or_include_current(true);
	task->nonconst_residue_task(1).prevent_repacking();
	protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
	repack.apply(pose);
}


void design(PoseWrap & pw, ScoreFunctionOP sf) {
	core::pose::Pose & pose(pw.pose);
	using namespace core::pack::task;
	utility::vector1< bool > aas(20,true);
	aas[core::chemical::aa_cys] = false;
	aas[core::chemical::aa_lys] = false;
	// aas[core::chemical::aa_glu] = false;
	PackerTaskOP task = TaskFactory::create_packer_task(pose);
	task->or_include_current(true);
	task->nonconst_residue_task(1).prevent_repacking();
	// task->nonconst_residue_task(2).prevent_repacking();
	for(Size i = 2; i <= task->total_residue(); ++i) {
		// if( pose.residue(i).name3() == "PRO" || pose.residue(i).name3() == "GLY" ) {
		// 	task->restrict_to_repacking();
		// } else {
			task->nonconst_residue_task(i).restrict_absent_canonical_aas(aas);
		// }
	}
	protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
	repack.apply(pose);


}

void flxbb_nobu(PoseWrap & pw, ScoreFunctionOP sfd, ScoreFunctionOP sfr) {
	protocols::flxbb::FlxbbDesign flxbb(sfd,sfr);
	flxbb.apply(pw.pose);
}

class	LigChiMover : public protocols::moves::Mover {
	void apply( core::pose::Pose & pose ) {
		Real a =  5.0*numeric::random::gaussian();
		Real b = 30.0*numeric::random::gaussian();
		Real c = 30.0*numeric::random::uniform();
		Real d = 30.0*numeric::random::uniform();
		pose.set_chi(1,1, pose.chi(1,1) + a );
		pose.set_chi(2,1, pose.chi(2,1) + b );
		pose.set_chi(3,1, pose.chi(3,1) + c );
		pose.set_chi(4,1, pose.chi(4,1) + d );
		// pose.set_psi(2,pose.psi(2)+10.0*numeric::random::gaussian());
	}
	void
	parse_my_tag(
		utility::tag::TagCOP const,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const &
	) {}
};

core::scoring::ScoreFunctionOP
cen_fold(PoseWrap & pw, Size NCYCLES, core::fragment::FragSetOP frags, Real temp=2.0 ) {
	using namespace core;
	using namespace scoring;
	using namespace protocols::moves;
	core::pose::Pose & pose(pw.pose);

	RandomMoverOP random_mover = new RandomMover;
	random_mover->add_mover(new LigChiMover,0.1);
	if( frags() != NULL ) {
		protocols::moves::MoverOP fragins = new protocols::abinitio::ClassicFragmentMover(frags);
		random_mover->add_mover(fragins,0.9);
	}

	scoring::ScoreFunctionOP sf0 = scoring::ScoreFunctionFactory::create_score_function( "score0" );
	scoring::ScoreFunctionOP sf1 = scoring::ScoreFunctionFactory::create_score_function( "score1" );
	scoring::ScoreFunctionOP sf2 = scoring::ScoreFunctionFactory::create_score_function( "score2" );
	scoring::ScoreFunctionOP sf3 = scoring::ScoreFunctionFactory::create_score_function( "score3" );
	scoring::ScoreFunctionOP sf5 = scoring::ScoreFunctionFactory::create_score_function( "score5" );
	sf0 = new scoring::symmetry::SymmetricScoreFunction(*sf0);
	sf1 = new scoring::symmetry::SymmetricScoreFunction(*sf1);
	sf2 = new scoring::symmetry::SymmetricScoreFunction(*sf2);
	sf3 = new scoring::symmetry::SymmetricScoreFunction(*sf3);
	sf5 = new scoring::symmetry::SymmetricScoreFunction(*sf5);
	sf0->set_weight(scoring::rg,0.0);
	sf1->set_weight(scoring::rg,0.0);
	sf2->set_weight(scoring::rg,0.0);
	sf3->set_weight(scoring::rg,0.0);
	sf5->set_weight(scoring::rg,0.0);
	sf0->set_weight(scoring::sym_lig,10.0);
	sf1->set_weight(scoring::sym_lig,10.0);
	sf2->set_weight(scoring::sym_lig,10.0);
	sf3->set_weight(scoring::sym_lig,10.0);
	sf5->set_weight(scoring::sym_lig,10.0);

	core::pose::Pose init = pose;

	pose = init;

	MonteCarloOP mc = new MonteCarlo( pose, *sf0, temp );
	mc->set_autotemp( true, temp );
	mc->set_temperature( temp );

	TR << "stage 1" << std::endl;
	RepeatMover( new TrialMover( random_mover, mc ), NCYCLES ).apply( pose );
	mc->reset( pose );
	// pose.dump_pdb("stage1.pdb");

	TR << "stage 2" << std::endl;
	mc = new MonteCarlo( pose, *sf1, 2.0 );
	RepeatMover( new TrialMover( random_mover, mc ), NCYCLES ).apply( pose );
	mc->reset( pose );
	// design(pose,sf1);
	// pose.dump_pdb("stage2.pdb");

	// for( Size i = 1; i <= 5; ++i ) {
	// 	TR << "stage 3 " << i << std::endl;
	// 	mc = new MonteCarlo( pose, *sf2, 2.0 );
	// 	RepeatMover( new TrialMover( random_mover, mc ), NCYCLES ).apply( pose );
	// 	mc->reset( pose );
	// 	// design(pose,sf2);
	// 	mc = new MonteCarlo( pose, *sf5, 2.0 );
	// 	RepeatMover( new TrialMover( random_mover, mc ), NCYCLES ).apply( pose );
	// 	mc->reset( pose );
	// 	// design(pose,sf5);
	// }
	// // pose.dump_pdb("stage3.pdb");

	for( Size i = 1; i <= 5; ++i ) {
		TR << "stage 3" << i << std::endl;
		mc = new MonteCarlo( pose, *sf3, 2.0 );
		RepeatMover( new TrialMover( random_mover, mc ), NCYCLES ).apply( pose );
		mc->reset( pose );
	}

	for( Size i = 1; i <= 5; ++i ) {
		TR << "stage 4" << i << std::endl;
		mc = new MonteCarlo( pose, *sf3, 2.0 );
		RepeatMover( new TrialMover( random_mover, mc ), NCYCLES ).apply( pose );
		mc->reset( pose );
		// design(pose,sf3);
	}

	return sf3;
}

core::scoring::ScoreFunctionOP
fa_refine_and_design(PoseWrap & pw, Size /*NCYCLE*/) {
	using namespace core;
	using namespace scoring;
	switch_to_fa(pw);
	core::pose::Pose & pose(pw.pose);
	// Size nres_mono = (pose.n_residue()-4)/4;
	ScoreFunctionOP sf1,sf2,sf3,sf4;
	sf1 = getScoreFunction();
	sf2 = getScoreFunction();
	sf3 = getScoreFunction();
	sf4 = getScoreFunction();
	sf1->set_weight(fa_rep,0.025);
	sf2->set_weight(fa_rep,0.050);
	sf3->set_weight(fa_rep,0.100);
	sf4->set_weight(fa_rep,0.200);
	// sf1->set_weight(fa_dun,0.100);
	// sf2->set_weight(fa_dun,0.150);
	// sf3->set_weight(fa_dun,0.200);
	// sf4->set_weight(fa_dun,0.300);
	// sf1->set_weight(fa_sol,0.800);
	// sf2->set_weight(fa_sol,0.800);
	// sf3->set_weight(fa_sol,0.800);
	// sf4->set_weight(fa_sol,0.800);
	sf1->set_weight(scoring::sym_lig,10.0);
	sf2->set_weight(scoring::sym_lig,10.0);
	sf3->set_weight(scoring::sym_lig,10.0);
	sf4->set_weight(scoring::sym_lig,10.0);
	sf1 = new symmetry::SymmetricScoreFunction(*sf1);
	sf2 = new symmetry::SymmetricScoreFunction(*sf2);
	sf3 = new symmetry::SymmetricScoreFunction(*sf3);
	sf4 = new symmetry::SymmetricScoreFunction(*sf4);
	core::pose::Pose best = pose;

	// the old way
	// for(Size i = 1; i <= NCYCLE; ++i ) {
	// 	for(Size k = 1; k <= 4; ++k) {
	// 		repack(pw,sf1); minimize(pw,sf1);
	// 		repack(pw,sf2); minimize(pw,sf2);
	// 		repack(pw,sf3); minimize(pw,sf3);
	// 		repack(pw,sf4); minimize(pw,sf4);
	// 		TR << "repack/min1234    " << I(2,k) << " " << F(10,3,(*sf4)(pose)) << " " << pose.sequence().substr(0,nres_mono) << std::endl;
	// 		if( (*sf4)(best) >= (*sf4)(pose) ) best = pose;
	// 	}
	// 	pose = best;
	// 	for(Size k = 1; k <= 1; ++k) {
	// 		// design(pose,sf1); minimize(pose,sf1);
	// 		//design(pose,sf2); minimize(pose,sf2);
	// 		design(pw,sf3); minimize(pw,sf3);
	// 		design(pw,sf4); minimize(pw,sf4);
	// 		TR << "deisgn/min1234   " << I(2,k) << " " << F(10,3,(*sf4)(pose)) << " " << pose.sequence().substr(0,nres_mono) << std::endl;
	// 		if( (*sf4)(best) >= (*sf4)(pose) ) best = pose;
	// 	}
	// 	pose = best;
	// }

	// nobu's way
	flxbb_nobu(pw,sf3,sf4);

	// switch_to_fa(pw);
	// repack(pw,sf4); minimize(pw,sf4);


	return sf4;
}

void report( std::string tag, PoseWrap & pw, core::pose::Pose & cenpose, ScoreFunctionOP sf_fa, std::ostringstream & oss, Real censcore ) {
	using namespace core;
	core::pose::Pose & pose(pw.pose);
	Size nres_mono = pw.nres;//(pose.n_residue()-4)/4;
	Real ref_rep  =  2827.2 + 2392.02; // pre-calculated

	core::pose::Pose posebefore = pose;
	sf_fa->set_weight(scoring::holes_min,0.005);
	// std::cout << "BEFORE HOLES MIN: " << (*sf4)(pose) << std::endl;
	(*sf_fa)(pose);
	core::Real before = pw.pose.energies().total_energies()[core::scoring::holes_min];
	minimize(pw,sf_fa);
	(*sf_fa)(pose);
	TR << "holes_min " << before << " " << pw.pose.energies().total_energies()[core::scoring::holes_min] << std::endl;
	// sf4->show(pose);

	// for( Size i = 0; i < 4; ++i ) {
	// 	Size nres_mono = (pose.n_residue()-4)/4;
	// 	ref_rep += pose.energies().residue_total_energies(i*nres_mono+1)[scoring::fa_rep];
	//    	ref_rep += pose.energies().residue_total_energies(i*nres_mono+2)[scoring::fa_rep];
	// }
	Real sym_lig      = pose.energies().total_energies()[ scoring::sym_lig      ];
	Real fa_atr       = pose.energies().total_energies()[ scoring::fa_atr       ];
 	Real fa_rep       = pose.energies().total_energies()[ scoring::fa_rep       ] - ref_rep;
	Real fa_sol       = pose.energies().total_energies()[ scoring::fa_sol       ];
	Real fa_intra_rep = pose.energies().total_energies()[ scoring::fa_intra_rep ];
	Real pro_close    = pose.energies().total_energies()[ scoring::pro_close    ];
	Real fa_pair      = pose.energies().total_energies()[ scoring::fa_pair      ];
	Real hbond_sr_bb  = pose.energies().total_energies()[ scoring::hbond_sr_bb  ];
	Real hbond_lr_bb  = pose.energies().total_energies()[ scoring::hbond_lr_bb  ];
	Real hbond_bb_sc  = pose.energies().total_energies()[ scoring::hbond_bb_sc  ];
	Real hbond_sc     = pose.energies().total_energies()[ scoring::hbond_sc     ];
	Real fa_dun       = pose.energies().total_energies()[ scoring::fa_dun       ];
	Real p_aa_pp      = pose.energies().total_energies()[ scoring::p_aa_pp      ];
	Real ref          = pose.energies().total_energies()[ scoring::ref          ];
	Real surface      = pose.energies().total_energies()[ scoring::surface      ];
	Real rholes       = pose.energies().total_energies()[ scoring::holes_min    ];
   ref_rep *= sf_fa->get_weight(scoring::fa_rep);
	Real score = (*sf_fa)(pose) - ref_rep;
	oss << LJ(20,tag)          << " "
	    << F(8,3,score/nres_mono)<< " "
		 << I(3,nres_mono)      << " "
	    << F(8,3,score)        << " "
	    << F(8,3,censcore)     << " "
		 << F(8,3,sym_lig     ) << " "
       << F(8,3,fa_atr      ) << " "
       << F(8,3,fa_rep      ) << " "
       << F(8,3,fa_sol      ) << " "
       << F(8,3,fa_intra_rep) << " "
       << F(8,3,pro_close   ) << " "
       << F(8,3,fa_pair     ) << " "
       << F(8,3,hbond_sr_bb ) << " "
       << F(8,3,hbond_lr_bb ) << " "
       << F(8,3,hbond_bb_sc ) << " "
       << F(8,3,hbond_sc    ) << " "
       << F(8,3,fa_dun      ) << " "
       << F(8,3,p_aa_pp     ) << " "
       << F(8,3,ref         ) << " "
       << F(8,3,surface     ) << " "
       << F(8,3,rholes      ) << " "
	    << pose.sequence().substr(0,pw.nres) << std::endl;



	if( score/nres_mono <= basic::options::option[basic::options::OptionKeys::in::file::silent_energy_cut]() ) {
		using namespace basic::options;
		cenpose.dump_pdb(std::string(option[OptionKeys::out::file::o]())+"/spiro"+tag+"_cen.pdb.gz");
		posebefore .dump_pdb(std::string(option[OptionKeys::out::file::o]())+"/spiro"+tag+"_fa.pdb.gz");
		pose   .dump_pdb(std::string(option[OptionKeys::out::file::o]())+"/spiro"+tag+"_holes.pdb.gz");
		std::cout << oss.str();
		std::cout.flush();
		oss.clear();
		oss.str("");
	}

}

std::string make_rand_ss() {
	bool twohelix = uniform() < 0.8;
	Size len1=(Size)(uniform()*15+7), len2=(Size)(uniform()*20+7);
	std::string SS = "L";
	if( uniform() < 0.5 ) SS += "L";
	for(Size i=0;i<len1;i++) SS += "H";
	if( !twohelix ) return SS;
	SS += "LL";
	if( uniform() < 0.5 ) SS += "L"; // 50/50 3 res loop
	for(Size i=0;i<len2;i++) SS += "H";
	return SS;
}


int
main( int argc, char * argv [] )
{

	try {


	using namespace core;
	using namespace pose;
	using namespace protocols;
	using namespace moves;
	using namespace ObjexxFCL::format;
	using numeric::random::uniform;

	devel::init(argc,argv);
	std::ostringstream oss;
	std::map<std::string, utility::vector1<core::fragment::FragDataOP> > fds = get_frags_map( true );

	while(true) {

		std::string SS = make_rand_ss();
		// SS = "LHHHHHHHH";
		Size NRES = SS.size()+1;
		std::string SEQ = "ZG";
		while( SEQ.size() < NRES ) SEQ += "V";
		// SS  = std::string("ZG");
		// SEQ = std::string("ZG");
		// NRES = SS.size();

// #ifdef NDEBUG
		Size NCEN = 30*NRES, NFA = 10;
// #else
		// Size NCEN = 3*NRES, NFA = 1;
// #endif

		PoseWrap pw = make_pose(SEQ);
		TR << "made pose, size: " << pw.pose.n_residue() << std::endl;

		std::string tag = string_of(uniform());

		core::fragment::FragSetOP frags = make_frag_set(2,SS,fds);
		ScoreFunctionOP sf_cen = cen_fold(pw,NCEN,frags);
		Real censcore = (*sf_cen)(pw.pose);
		Pose cenpose = pw.pose;
		// pw.pose.dump_pdb("cen_after.pdb");

		ScoreFunctionOP sf_fa = fa_refine_and_design(pw,NFA);

		// std::cerr << "ref_rep ";
		// for( Size i = 0; i < 2; ++i ) {
		// 	std::cerr << pw.pose.energies().residue_total_energies(i*NRES+1)[scoring::fa_rep] << " ";
		// }
		// std::cerr << std::endl;

		report( tag, pw, cenpose, sf_fa, oss, censcore );
	}


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

}
