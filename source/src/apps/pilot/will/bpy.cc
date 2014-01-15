// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file /src/apps/pilat/will/coiled_coil.cc
/// @brief samples coiled coils

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>

#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/FragSet.hh>
#include <devel/init.hh>
#include <devel/init.hh>
#include <basic/database/open.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/smhybrid.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/optimizeH.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/util.hh>
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/packstat/compute_sasa.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <basic/Tracer.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/flxbb/FlxbbDesign.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/symmetric_docking/SymDockingLowRes.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>


using core::Size;
using core::Real;
using numeric::random::uniform;
using core::id::AtomID;
using core::id::DOF_ID;
using core::scoring::ScoreFunctionOP;


static basic::Tracer TR("bpy");

static numeric::random::RandomGenerator RG(6054212);

struct PoseWrap {
	PoseWrap() : hascst(false) {}
	core::pose::Pose pose;
	Size nsub,attach,nres;
	bool hascst;
	void dump_pdb(std::string fname) {
		utility::io::ozstream out(fname);
		pose.dump_pdb(out);
		if(basic::options::option[basic::options::OptionKeys::smhybrid::add_metal_at_0]())
			out << "HETATM 9999 ZN    ZN A 999       0.000   0.000   0.000  1.00  0.00              " << std::endl;
		if(basic::options::option[basic::options::OptionKeys::smhybrid::add_cavities]())
			core::scoring::packstat::output_packstat_pdb(pose,out);

	}
};

void switch_to_fa(PoseWrap & pw) {
	core::pose::Pose const cen( pw.pose );
	core::util::switch_to_residue_type_set( pw.pose, core::chemical::FA_STANDARD );
	core::id::AtomID_Map<bool> missing;
	core::pose::initalize_atomid_map(missing,pw.pose,false);
	bool anymissing = false;
	for(Size attach = pw.attach; attach < pw.nsub*pw.nres; attach += pw.nres) {
		for(Size i = 1; i <= pw.pose.residue(attach).natoms(); ++i) {
			std::string aname = pw.pose.residue(attach).atom_name(i);
			if( cen.residue(attach).type().has(aname) ) {
				// std::cerr << "setting coords for res " << attach << " atom " << aname << std::endl;
				pw.pose.set_xyz(core::id::AtomID(i,attach),cen.residue(attach).xyz(aname));
			} else {
				// std::cerr << "switch_to_fa: res " << attach << " missing " << aname << std::endl;
				missing[core::id::AtomID(i,attach)] = true;
				anymissing = true;
			}
		}
	}
	if(anymissing) pw.pose.conformation().fill_missing_atoms( missing );
}

void switch_to_cen(PoseWrap & pw) {
	core::pose::Pose const fa( pw.pose );
	protocols::toolbox::switch_to_residue_type_set( pw.pose, core::chemical::CENTROID );
	core::id::AtomID_Map<bool> missing;
	core::pose::initalize_atomid_map(missing,pw.pose,false);
	bool anymissing = false;
	for(Size attach = pw.attach; attach < pw.nsub*pw.nres; attach += pw.nres) {
		for(Size i = 1; i <= pw.pose.residue(attach).natoms(); ++i) {
			std::string aname = pw.pose.residue(attach).atom_name(i);
			if( fa.residue(attach).type().has(aname) ) {
				// std::cerr << "setting coords for res " << attach << " atom " << aname << std::endl;
				pw.pose.set_xyz(core::id::AtomID(i,attach),fa.residue(attach).xyz(aname));
			} else {
				// std::cerr << "switch_to_cen: res " << attach << " missing " << aname << std::endl;
				missing[core::id::AtomID(i,attach)] = true;
				anymissing = true;
			}
		}
	}
	if(anymissing) pw.pose.conformation().fill_missing_atoms( missing );
}


void add_apc(core::pose::Pose & pose, core::id::AtomID aid1, core::id::AtomID aid2, core::Real mean, core::Real sd) {
	core::scoring::constraints::ConstraintOP cc = new core::scoring::constraints::AtomPairConstraint(
		aid1, aid2, new core::scoring::constraints::HarmonicFunc(mean,sd) );
	pose.add_constraint(cc);
}


void read_fragdata( FragDataCOPs& fds, utility::io::izstream & in, bool /*design = false*/ ) {
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

std::map<std::string, FragDataCOPs >
get_frags_map( bool design = false ) {
	using namespace core::fragment;
	TR << "reading frags" << std::endl;
	utility::io::izstream in;
	std::map<std::string, FragDataCOPs > fds;
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


core::fragment::FragSetOP make_frag_set(Size start, std::string ss, std::map<std::string, FragDataCOPs > const& fds) {
	using namespace core::fragment;
	FragSetOP frags = new ConstantLengthFragSet();
	Size const stop = ss.size() + start - 3;
	if(start >= stop) return NULL;
	for( Size i = start; i <= stop; ++i ) {
		FrameOP frame = new Frame(i,3);
		FragDataCOPs::const_iterator beg = fds[ss.substr(i-start,3)].begin();
		FragDataCOPs::const_iterator end = fds[ss.substr(i-start,3)].end();
		for( FragDataCOPs::const_iterator fi = beg; fi != end; ++fi ) {
			frame->add_fragment(*fi);
		}
		frags->add(frame);
	}
	return frags;
}

class ChiMover : public protocols::moves::Mover {
	Size residue_;
	core::pack::dunbrack::SingleResidueRotamerLibraryCAP lib_;
	core::pack::dunbrack::RotamerLibraryScratchSpace scratch_;
public:
	ChiMover(core::pose::Pose const & pose, Size residue) : residue_(residue) {
		using namespace core::pack::dunbrack;
		using namespace core::chemical;
		std::string name3 = pose.residue(residue).name3();
		if(name3=="BPY") name3 = "TRP";
		if(name3=="TIA") name3 = "LYS";
		ResidueTypeSetCAP rs( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
		lib_ = core::pack::dunbrack::RotamerLibrary::get_rsd_library( rs->name_map(name3) );
		assert(lib_);
	}
	void apply( core::pose::Pose & pose ) {
		using namespace core::pack::dunbrack;
		ChiVector chis;
		lib_->assign_random_rotamer_with_bias(pose.residue(residue_),scratch_,RG,chis,true);
		for(size_t i = 1; i <= chis.size(); ++i)	{
			pose.set_chi(i,residue_,chis[i]);
		}
		// Real a =  5.0*numeric::random::gaussian();
		// Real b = 30.0*numeric::random::gaussian();
		// Real c = 30.0*numeric::random::uniform();
		// Real d = 30.0*numeric::random::uniform();
		// pose.set_chi(1,1, pose.chi(1,1) + a );
		// pose.set_chi(2,1, pose.chi(2,1) + b );
		// pose.set_chi(3,1, pose.chi(3,1) + c );
		// pose.set_chi(4,1, pose.chi(4,1) + d );
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

void design(PoseWrap & pw, ScoreFunctionOP sf) {
	core::pose::Pose & pose(pw.pose);
	using namespace core::pack::task;
	utility::vector1< bool > aas(20,true);
	aas[core::chemical::aa_cys] = false;
	PackerTaskOP task = TaskFactory::create_packer_task(pose);
	task->or_include_current(true);
	task->nonconst_residue_task(pw.attach).prevent_repacking();
	for(Size i = 1; i <= task->total_residue(); ++i) {
		task->nonconst_residue_task(i).restrict_absent_canonical_aas(aas);
	}
	protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
	repack.apply(pose);


}

void minimize(PoseWrap & pw, ScoreFunctionOP sf) {
	core::pose::Pose & pose(pw.pose);
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	movemap->set_chi(true);
	movemap->set_bb(true);
	protocols::simple_moves::symmetry::SymMinMover m( movemap, sf, "dfpmin_armijo_nonmonotone", 1e-8, true );
	m.apply(pose);
}

core::scoring::ScoreFunctionOP
cen_fold(PoseWrap & pw, Size NCYCLES, core::fragment::FragSetOP frags, Real temp=2.0 ) {
	using namespace core;
	using namespace scoring;
	using namespace protocols::moves;
	switch_to_cen(pw);
	core::pose::Pose & pose(pw.pose);

	RandomMoverOP random_mover = new RandomMover;
	random_mover->add_mover(new ChiMover(pw.pose,pw.attach),0.1);
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
	// sf0->set_weight(scoring::rg,5.0);
	// sf0->set_weight(scoring::rg,5.0);
	// sf1->set_weight(scoring::rg,5.0);
	// sf3->set_weight(scoring::rg,5.0);
	// sf5->set_weight(scoring::rg,5.0);

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

	// pose.dump_pdb("stage2.pdb");

	// NO SHEETS AT THE MOMENT....
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

	for( Size i = 1; i <= 10; ++i ) {
		TR << "stage 3" << i << std::endl;
		mc = new MonteCarlo( pose, *sf3, 2.0 );
		RepeatMover( new TrialMover( random_mover, mc ), NCYCLES/5 ).apply( pose );
		mc->reset( pose );
	}
	if(basic::options::option[basic::options::OptionKeys::smhybrid::abinitio_design]()) {
		design(pw,sf3);
	}
	for( Size i = 1; i <= 10; ++i ) {
		TR << "stage 3" << i << std::endl;
		mc = new MonteCarlo( pose, *sf3, 2.0 );
		RepeatMover( new TrialMover( random_mover, mc ), NCYCLES/5 ).apply( pose );
		mc->reset( pose );
	}


	return sf3;
}



PoseWrap make_pose(Size Nprev, Size Npost) {
	using namespace core;
	using namespace chemical;
	using namespace conformation;
	PoseWrap pw;
	pw.attach = Nprev+1;
	pw.nsub = 3;
	pw.nres = 1 + Nprev + Npost;

	core::import_pose::pose_from_pdb(pw.pose,options::option[options::OptionKeys::in::file::s]()[1]);

	core::pose::remove_lower_terminus_type_from_pose_residue(pw.pose,1);
	core::pose::remove_upper_terminus_type_from_pose_residue(pw.pose,pw.pose.n_residue());
	// pw.pose.dump_pdb("make_pose_init.pdb");


	// string seq = seq1 + "X[BPY]" + seq2''
	// core::pose::make_pose_from_sequence(pose,seq,"fa_standard",true);


	ResidueTypeSetCAP residue_set( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
	for( Size i = 1; i <= Nprev; ++i ) {
		std::string name3 = name_from_aa(aa_from_oneletter_code('L'));
		pw.pose.prepend_polymer_residue_before_seqpos(*ResidueFactory::create_residue(residue_set->name_map(name3)),1,true);
	}
	core::pose::add_lower_terminus_type_to_pose_residue(pw.pose,1);
	// pw.pose.dump_pdb("init1.pdb");
	for( Size i = 1; i <= Npost; ++i ) {
		std::string name3 = name_from_aa(aa_from_oneletter_code('L'));
		pw.pose.append_residue_by_bond(*ResidueFactory::create_residue(residue_set->name_map(name3)),true);
	}
	core::pose::add_upper_terminus_type_to_pose_residue(pw.pose,pw.pose.n_residue());

	// set right fold tree... just re-rooting... assumes symmetry machinery will respect this
	kinematics::FoldTree ft = pw.pose.fold_tree();
	ft.reorder(pw.attach);
	ft.set_root_atomno(basic::options::option[basic::options::OptionKeys::smhybrid::anchor_atomno]()); // should be "zinc"
	pw.pose.fold_tree(ft);


	// pw.pose.dump_pdb("make_pose_asym.pdb");
	// for( Size i = 1; i <= pw.nres; ++i ) {
	// 	pw.pose.set_phi  (i,-60.16731);
	// 	pw.pose.set_psi  (i,-45.19451);
	// 	pw.pose.set_omega(i,180.00000);
	// }
	// { ChiMover cm(pw.pose,pw.attach); cm.apply(pw.pose); }
	// pw.pose.dump_pdb("make_pose_helix.pdb");



	core::pose::symmetry::make_symmetric_pose( pw.pose );



	// pw.pose.dump_pdb("make_pose_symm.pdb");
	//
	// { ChiMover cm(pw.pose,pw.attach); cm.apply(pw.pose); }
	// pw.pose.dump_pdb("make_pose_symm1.pdb");
	// std::cerr << pw.pose.fold_tree() << std::endl;
	// pw.pose.set_chi(1,pw.attach,pw.pose.chi(1,pw.attach) + uniform()*10.-5. );
	// pw.pose.set_chi(2,pw.attach,pw.pose.chi(2,pw.attach) + uniform()*10.-5. );
	// for( Size i = 1; i <= pw.nres; ++i ) {
	// 	pw.pose.set_phi  (i,pw.pose.phi(i)   + uniform()*10.- 5.);
	// 	pw.pose.set_psi  (i,pw.pose.psi(i)   + uniform()*10.- 5.);
	// 	pw.pose.set_omega(i,pw.pose.omega(i) + uniform()* 1.-.5 );
	// }

	// switch_to_cen(pw);
	// { ChiMover cm(pw.pose,pw.attach); cm.apply(pw.pose); }
	// pw.pose.dump_pdb("make_pose_symm2.pdb");
	// switch_to_fa(pw);
	// { ChiMover cm(pw.pose,pw.attach); cm.apply(pw.pose); }
	// pw.pose.dump_pdb("make_pose_symm3.pdb");
	//
	// std::exit(-1);

	return pw;
}

ScoreFunctionOP flxbb_nobu(PoseWrap & pw) {
	using namespace core::scoring;
	ScoreFunctionOP sf3,sf4;
	sf3 = getScoreFunction();
	sf4 = getScoreFunction();
	sf3->set_weight(fa_rep,0.100);
	sf4->set_weight(fa_rep,0.400);
	sf3 = new symmetry::SymmetricScoreFunction(*sf3);
	sf4 = new symmetry::SymmetricScoreFunction(*sf4);

	protocols::flxbb::FlxbbDesign flxbb(sf3,sf4);
	flxbb.apply(pw.pose);

	return sf4;
}

void report( PoseWrap & pw, ScoreFunctionOP sf_fa, std::ostringstream & oss, Real censcore ) {
	using namespace core;
	using namespace ObjexxFCL::format;
	core::pose::Pose & pose(pw.pose);
	Size nres_mono = pw.nres;//(pose.n_residue()-4)/4;
	std::string tag = string_of(uniform());

	// core::pose::Pose posebefore = pose;
	// sf_fa->set_weight(scoring::holes_min,0.0005);
	// std::cout << "BEFORE HOLES MIN: " << (*sf4)(pose) << std::endl;
	(*sf_fa)(pose);
	// core::Real before = pw.pose.energies().total_energies()[core::scoring::holes_min];
	// minimize(pw,sf_fa);
	// (*sf_fa)(pose);
	// TR << "holes_min " << before << " " << pw.pose.energies().total_energies()[core::scoring::holes_min] << std::endl;
	// sf4->show(pose);

	// for( Size i = 0; i < 4; ++i ) {
	// 	Size nres_mono = (pose.n_residue()-4)/4;
	// 	ref_rep += pose.energies().residue_total_energies(i*nres_mono+1)[scoring::fa_rep];
	//    	ref_rep += pose.energies().residue_total_energies(i*nres_mono+2)[scoring::fa_rep];
	// }
	Real sym_lig      = pose.energies().total_energies()[ scoring::sym_lig      ];
	Real fa_atr       = pose.energies().total_energies()[ scoring::fa_atr       ];
	Real fa_rep       = pose.energies().total_energies()[ scoring::fa_rep       ];
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
	Real score = (*sf_fa)(pose);
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
		// cenpose.dump_pdb(std::string(option[OptionKeys::out::file::o]())+"/bpy"+tag+"_cen.pdb.gz");
		// posebefore .dump_pdb(std::string(option[OptionKeys::out::file::o]())+"/bpy"+tag+"_fa.pdb.gz");
		pw.dump_pdb(std::string(option[OptionKeys::out::file::o]())+"/bpy"+tag+"_holes.pdb.gz");

		std::cout << oss.str();
		std::cout.flush();
		oss.clear();
		oss.str("");
	}

}




std::string make_rand_ss(core::Size len) {
	std::string ss;
	for(int i = 1; i <= uniform()*15+5; ++i) ss += "H";
	while(ss.size()<len) {
		for(int i = 1; i <= uniform()*6   ; ++i) ss += "L";
		for(int i = 1; i <= uniform()*15+5; ++i) ss += "H";
	}
	return ss;
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
	std::map<std::string, FragDataCOPs > fds = get_frags_map( true );
	std::ostringstream oss;

	while(true) {
		// make topology
		Size N = basic::options::option[basic::options::OptionKeys::smhybrid::nres_mono]();
		Size Ncycle = basic::options::option[basic::options::OptionKeys::smhybrid::abinitio_cycles]();
		std::string SS = make_rand_ss(10);
		Size attach = basic::options::option[basic::options::OptionKeys::smhybrid::attach_rsd]();
		if( 0 == attach ) attach = (Size)(uniform()*((Real)(SS.size()-min(20.,2.*N/3.))))+min((Size)10,(Size)(N/3.));
		TR << "attaching at " << attach << std::endl;
		PoseWrap pw = make_pose(attach-1,SS.size()-attach);

		// fold up
		core::fragment::FragSetOP frags = make_frag_set(1,SS,fds);
		scoring::ScoreFunctionOP sf_cen = cen_fold(pw,Ncycle,frags);
		Real censcore = (*sf_cen)(pw.pose);

		// design and refine using Nobu's protocol
		switch_to_fa(pw);
		ScoreFunctionOP sf_fa = flxbb_nobu(pw);

		report(pw,sf_fa,oss,censcore);

	}


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

}
