// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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

#include <core/fragment/FragmentIO.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/FragSet.hh>
#include <devel/init.hh>
#include <basic/database/open.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/kinematics/MoveMap.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <basic/Tracer.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/abinitio/FragmentMover.hh>
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
#include <core/pose/util.hh>
#include <protocols/toolbox/SwitchResidueTypeSet.hh>


using numeric::conversions::radians;

static basic::Tracer TR( "dubois" );

using core::Size;
using core::Real;
using core::id::AtomID;
typedef numeric::xyzVector<Real> Vec;
typedef utility::vector1<Vec>    Vecs;
typedef numeric::xyzMatrix<Real> Mat;
using protocols::moves::MoverOP;
using core::scoring::ScoreFunctionOP;
using numeric::random::uniform;

void switch_to_fa(core::pose::Pose & pose) {
	protocols::toolbox::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD );
	pose.set_dof(core::id::DOF_ID(core::id::AtomID(6,2),core::id::PHI  ),0);
	pose.set_dof(core::id::DOF_ID(core::id::AtomID(6,2),core::id::THETA),PI);
}

void switch_to_cen(core::pose::Pose & pose) {
	protocols::toolbox::switch_to_residue_type_set( pose, core::chemical::CENTROID );
	pose.set_dof(core::id::DOF_ID(core::id::AtomID(6,2),core::id::PHI  ),0);
	pose.set_dof(core::id::DOF_ID(core::id::AtomID(6,2),core::id::THETA),PI);
}

void minimize(core::pose::Pose & pose, ScoreFunctionOP sf) {
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	movemap->set_jump(false);
	movemap->set_chi(true);
	movemap->set_bb(true);
	// movemap->set_bb(1,false);
	protocols::simple_moves::symmetry::SymMinMover m( movemap, sf, "lbfgs_armijo_nonmonotone", 1e-2, true );
	m.apply(pose);
}


class SeqBBTorSRFD : public core::fragment::BBTorsionSRFD {
public:
	bool apply( core::pose::Pose& pose, Size seqpos ) const {
		if ( seqpos >= 3 && 0.001 > pose.energies().residue_total_energies_unsafe(seqpos)[core::scoring::sym_lig] ) {
			pose.conformation().replace_centroid(seqpos,sequence_);
		}
		return BBTorsionSRFD::apply( pose, seqpos );
	}
};

void read_fragdata( utility::vector1< core::fragment::FragDataOP > & fds, utility::io::izstream & in, bool design ) {
	using namespace core::fragment;
	Size n,count=0;
	while ( in >> n ) {
		std::string pdb;
		char buf[999];
		FragDataOP fd = new FragData;
		for ( Size i = 1; i <= n; ++i ) {
			utility::pointer::owning_ptr<SingleResidueFragData> srfd;
			if ( design ) srfd = new SeqBBTorSRFD;
			else       srfd = new BBTorsionSRFD;
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
get_frags_map( bool design = true ) {
	using namespace core::fragment;
	TR << "reading frags" << std::endl;
	utility::io::izstream in;
	std::map<std::string,utility::vector1<FragDataOP> > fds;
	basic::database::open(in,"sampling/ss_fragfiles/EEE.fragfile"); read_fragdata(fds["EEE"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/EEH.fragfile"); read_fragdata(fds["EEH"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/EEL.fragfile"); read_fragdata(fds["EEL"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/EHH.fragfile"); read_fragdata(fds["EHH"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/ELE.fragfile"); read_fragdata(fds["ELE"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/ELH.fragfile"); read_fragdata(fds["ELH"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/ELL.fragfile"); read_fragdata(fds["ELL"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/HEE.fragfile"); read_fragdata(fds["HEE"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/HEH.fragfile"); read_fragdata(fds["HEH"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/HEL.fragfile"); read_fragdata(fds["HEL"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/HHE.fragfile"); read_fragdata(fds["HHE"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/HHH.fragfile"); read_fragdata(fds["HHH"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/HHL.fragfile"); read_fragdata(fds["HHL"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/HLE.fragfile"); read_fragdata(fds["HLE"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/HLH.fragfile"); read_fragdata(fds["HLH"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/HLL.fragfile"); read_fragdata(fds["HLL"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/LEE.fragfile"); read_fragdata(fds["LEE"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/LEH.fragfile"); read_fragdata(fds["LEH"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/LEL.fragfile"); read_fragdata(fds["LEL"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/LHH.fragfile"); read_fragdata(fds["LHH"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/LLE.fragfile"); read_fragdata(fds["LLE"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/LLH.fragfile"); read_fragdata(fds["LLH"],in,design);
	basic::database::open(in,"sampling/ss_fragfiles/LLL.fragfile"); read_fragdata(fds["LLL"],in,design);
	return fds;
}


core::fragment::FragSetOP make_frag_set(Size start, std::string ss, std::map<std::string, utility::vector1<core::fragment::FragDataOP> > fds) {
	using namespace core::fragment;
	FragSetOP frags = new ConstantLengthFragSet();
	Size const stop = ss.size() + start - 3;
	for ( Size i = start; i <= stop; ++i ) {
		FrameOP frame = new Frame(i,3);
		utility::vector1<FragDataOP>::iterator beg = fds[ss.substr(i-start,3)].begin();
		utility::vector1<FragDataOP>::iterator end = fds[ss.substr(i-start,3)].end();
		for ( utility::vector1<FragDataOP>::iterator fi = beg; fi != end; ++fi ) {
			frame->add_fragment(*fi);
		}
		frags->add(frame);
	}
	return frags;
}

void test( core::pose::Pose & pose ) {
	using namespace core;
	using namespace scoring;


	std::cerr << "FT after symm" << std::endl;
	for ( Size i = 1; i <= pose.num_jump(); ++i ) {
		std::cerr << pose.fold_tree().jump_edge(i) << std::endl;
	}

	pose.dump_pdb("debug_sym0.pdb");
	for ( Size i = 1; i <= 10; ++i ) {
		// pose.set_chi( 1, 1, pose.chi(1,1)+20.0*uniform() );
		pose.set_chi( 2, 1, pose.chi(2,1)+20.0*uniform() );
		// pose.set_phi(    2, pose.phi(  2)+20.0+uniform() );
		// pose.set_psi(    2, pose.psi(  2)+20.0+uniform() );
		pose.dump_pdb("rot"+string_of(i)+".pdb");
	}

	pose.dump_pdb("min0.pdb");
	ScoreFunctionOP sf = get_score_function_legacy( PRE_TALARIS_2013_STANDARD_WTS );
	sf->set_weight(atom_pair_constraint,1.0);

	sf->show(pose);
	minimize(pose,sf);
	sf->show(pose);

	pose.dump_pdb("min1.pdb");

}

core::pose::Pose make_pose(std::string seq) {
	using namespace core;
	using namespace chemical;
	using namespace conformation;
	Size nres = seq.size();
	ResidueTypeSetCAP residue_set( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
	pose::Pose pose;
	core::import_pose::pose_from_file(pose,"input/nipphgly.pdb",residue_set, core::import_pose::PDB_file);
	core::pose::remove_upper_terminus_type_from_pose_residue(pose,2);
	core::pose::add_variant_type_to_pose_residue(pose,"VIRTUAL_NTERM",2);
	pose.set_dof(id::DOF_ID(id::AtomID(6,2),id::PHI  ),0);
	pose.set_dof(id::DOF_ID(id::AtomID(6,2),id::THETA),PI);
	for ( Size i = 3; i <= nres; ++i ) {
		std::string name3 = name_from_aa(aa_from_oneletter_code(seq[i-1]));
		pose.append_residue_by_bond(*ResidueFactory::create_residue(residue_set->name_map(name3)),true);
	}
	core::pose::symmetry::make_symmetric_pose( pose );
	switch_to_cen(pose);
	kinematics::FoldTree ft = pose.fold_tree();
	for ( Size i = 1; i <= 4; ++i ) {
		ft.jump_edge(i).start_atom() = "VN";
		ft.jump_edge(i).stop_atom() = "N";
	}
	pose.fold_tree(ft);
	pose.set_dof(id::DOF_ID(id::AtomID(6,2),id::PHI  ),0);
	pose.set_dof(id::DOF_ID(id::AtomID(6,2),id::THETA),PI);
	for ( Size i = 2; i <= seq.size(); i++ ) {
		pose.set_phi  (i,-60.16731);
		pose.set_psi  (i,-45.19451);
		pose.set_omega(i,180.00000);
	}
	return pose;
}

void repack(core::pose::Pose & pose, ScoreFunctionOP sf) {
	using namespace core::pack::task;
	PackerTaskOP task = TaskFactory::create_packer_task(pose);
	task->restrict_to_repacking();
	task->or_include_current(true);
	task->nonconst_residue_task(1).prevent_repacking();
	task->nonconst_residue_task(2).prevent_repacking();
	protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
	repack.apply(pose);
}


void design(core::pose::Pose & pose, ScoreFunctionOP sf) {
	using namespace core::pack::task;
	utility::vector1< bool > aas(20,true);
	aas[core::chemical::aa_cys] = false;
	aas[core::chemical::aa_asp] = false;
	aas[core::chemical::aa_glu] = false;
	PackerTaskOP task = TaskFactory::create_packer_task(pose);
	task->or_include_current(true);
	task->nonconst_residue_task(1).prevent_repacking();
	task->nonconst_residue_task(2).prevent_repacking();
	for ( Size i = 3; i <= task->size(); ++i ) {
		if ( pose.residue(i).name3() == "PRO" || pose.residue(i).name3() == "GLY" ) {
			task->restrict_to_repacking();
		} else {
			task->nonconst_residue_task(i).restrict_absent_canonical_aas(aas);
		}
	}
	protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
	repack.apply(pose);


}

class LigChiMover : public protocols::moves::Mover {
	void apply( core::pose::Pose & pose ) {
		Real a =  10.0*numeric::random::gaussian();
		Real b =  10.0*numeric::random::gaussian();
		Real n = 180.0*numeric::random::uniform();
		pose.set_chi(1,1, pose.chi(2,1) + b     );
		pose.set_chi(2,1, pose.chi(1,1) + a - n );
		pose.set_chi(3,1, pose.chi(1,1) + n     );
		// pose.set_psi(2,pose.psi(2)+10.0*numeric::random::gaussian());
	}
};

core::scoring::ScoreFunctionOP
cen_fold(core::pose::Pose & pose, Size NCYCLES, core::fragment::FragSetOP frags, Real temp=2.0 ) {
	using namespace core;
	using namespace scoring;
	using namespace protocols::moves;

	protocols::moves::MoverOP fragins = new protocols::abinitio::ClassicFragmentMover(frags);
	RandomMoverOP random_mover = new RandomMover;
	random_mover->add_mover(fragins,0.9);
	random_mover->add_mover(new LigChiMover,0.1);

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
	//  TR << "stage 3 " << i << std::endl;
	//  mc = new MonteCarlo( pose, *sf2, 2.0 );
	//  RepeatMover( new TrialMover( random_mover, mc ), NCYCLES ).apply( pose );
	//  mc->reset( pose );
	//  // design(pose,sf2);
	//  mc = new MonteCarlo( pose, *sf5, 2.0 );
	//  RepeatMover( new TrialMover( random_mover, mc ), NCYCLES ).apply( pose );
	//  mc->reset( pose );
	//  // design(pose,sf5);
	// }
	// // pose.dump_pdb("stage3.pdb");

	for ( Size i = 1; i <= 5; ++i ) {
		TR << "stage 3" << i << std::endl;
		mc = new MonteCarlo( pose, *sf3, 2.0 );
		RepeatMover( new TrialMover( random_mover, mc ), NCYCLES ).apply( pose );
		mc->reset( pose );
	}

	for ( Size i = 1; i <= 5; ++i ) {
		TR << "stage 4" << i << std::endl;
		mc = new MonteCarlo( pose, *sf3, 2.0 );
		RepeatMover( new TrialMover( random_mover, mc ), NCYCLES ).apply( pose );
		mc->reset( pose );
		// design(pose,sf3);
	}

	return sf3;
}

core::scoring::ScoreFunctionOP
fa_refine_and_design(core::pose::Pose & pose, Size NCYCLE) {
	using namespace core;
	using namespace scoring;
	Size nres_mono = (pose.size()-4)/4;
	ScoreFunctionOP sf1,sf2,sf3,sf4;
	sf1 = get_score_function_legacy( PRE_TALARIS_2013_STANDARD_WTS );
	sf2 = get_score_function_legacy( PRE_TALARIS_2013_STANDARD_WTS );
	sf3 = get_score_function_legacy( PRE_TALARIS_2013_STANDARD_WTS );
	sf4 = get_score_function_legacy( PRE_TALARIS_2013_STANDARD_WTS );
	sf1->set_weight(fa_rep,0.025);
	sf2->set_weight(fa_rep,0.050);
	sf3->set_weight(fa_rep,0.100);
	sf4->set_weight(fa_rep,0.200);
	sf1->set_weight(fa_dun,0.100);
	sf2->set_weight(fa_dun,0.150);
	sf3->set_weight(fa_dun,0.200);
	sf4->set_weight(fa_dun,0.300);
	sf1->set_weight(fa_sol,0.800);
	sf2->set_weight(fa_sol,0.800);
	sf3->set_weight(fa_sol,0.800);
	sf4->set_weight(fa_sol,0.800);
	sf1->set_weight(scoring::sym_lig,10.0);
	sf2->set_weight(scoring::sym_lig,10.0);
	sf3->set_weight(scoring::sym_lig,10.0);
	sf4->set_weight(scoring::sym_lig,10.0);
	sf1 = new symmetry::SymmetricScoreFunction(*sf1);
	sf2 = new symmetry::SymmetricScoreFunction(*sf2);
	sf3 = new symmetry::SymmetricScoreFunction(*sf3);
	sf4 = new symmetry::SymmetricScoreFunction(*sf4);
	core::pose::Pose best = pose;
	for ( Size i = 1; i <= NCYCLE; ++i ) {
		for ( Size k = 1; k <= 4; ++k ) {
			repack(pose,sf1); minimize(pose,sf1);
			repack(pose,sf2); minimize(pose,sf2);
			repack(pose,sf3); minimize(pose,sf3);
			repack(pose,sf4); minimize(pose,sf4);
			TR << "repack/min1234    " << I(2,k) << " " << F(10,3,(*sf4)(pose)) << " " << pose.sequence().substr(0,nres_mono) << std::endl;
			if ( (*sf4)(best) >= (*sf4)(pose) ) best = pose;
		}
		pose = best;
		// for(Size k = 1; k <= 1; ++k) {
		//  // design(pose,sf1); minimize(pose,sf1);
		//  //design(pose,sf2); minimize(pose,sf2);
		//  design(pose,sf3); minimize(pose,sf3);
		//  design(pose,sf4); minimize(pose,sf4);
		//  TR << "deisgn/min1234   " << I(2,k) << " " << F(10,3,(*sf4)(pose)) << " " << pose.sequence().substr(0,nres_mono) << std::endl;
		//  if( (*sf4)(best) >= (*sf4)(pose) ) best = pose;
		// }
		// pose = best;
	}
	return sf4;
}

void report( std::string tag, core::pose::Pose & pose, core::pose::Pose & cenpose, ScoreFunctionOP sf_fa, std::ostringstream & oss, Real censcore ) {
	using namespace core;
	Size nres_mono = (pose.size()-4)/4;
	Real ref_rep  =  2260.948; // pre-calculated
	// for( Size i = 0; i < 4; ++i ) {
	//  Size nres_mono = (pose.size()-4)/4;
	//  ref_rep += pose.energies().residue_total_energies(i*nres_mono+1)[scoring::fa_rep];
	//     ref_rep += pose.energies().residue_total_energies(i*nres_mono+2)[scoring::fa_rep];
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
		<< pose.sequence().substr(0,nres_mono) << std::endl;

	if ( score/nres_mono <= basic::options::option[basic::options::OptionKeys::in::file::silent_energy_cut]() ) {
		using namespace basic::options;
		cenpose.dump_pdb(std::string(option[OptionKeys::out::file::o]())+"/dubois"+tag+"_cen.pdb.gz");
		pose   .dump_pdb(std::string(option[OptionKeys::out::file::o]())+"/dubois"+tag+"_fa.pdb.gz");
		std::cout << oss.str();
		std::cout.flush();
		oss.clear();
		oss.str("");
	}

}

std::string make_rand_ss() {
	bool twohelix = uniform() < 0.8;
	Size len1=(uniform()*15+7), len2=(uniform()*20+7);
	std::string SS = "L";
	if ( uniform() < 0.5 ) SS += "L";
	for ( Size i=0; i<len1; i++ ) SS += "H";
	if ( !twohelix ) return SS;
	SS += "LL";
	if ( uniform() < 0.5 ) SS += "L"; // 50/50 3 res loop
	for ( Size i=0; i<len2; i++ ) SS += "H";
	return SS;
}
/////////////////////////////////////////////////////////////////////////////
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

		const Real PI = numeric::NumericTraits<Real>::pi();

		while ( true ) {

			std::string SS = make_rand_ss();
			Size NRES = SS.size()+1;
			std::string SEQ = "ZGSH";
			while ( SEQ.size() < NRES ) SEQ += "V";

#ifdef NDEBUG
		Size NCEN = 30*NRES, NFA = 10;
#else
			Size NCEN = 3*NRES, NFA = 1;
#endif

			Pose init = make_pose(SEQ);

			std::string tag = string_of(uniform());
			Pose pose = init;

			core::fragment::FragSetOP frags = make_frag_set(2,SS,fds);
			ScoreFunctionOP sf_cen = cen_fold(pose,NCEN,frags);
			Real censcore = (*sf_cen)(pose);
			Pose cenpose = pose;

			switch_to_fa(pose);
			pose.set_dof(id::DOF_ID(id::AtomID(6,2),id::PHI  ),0);
			pose.set_dof(id::DOF_ID(id::AtomID(6,2),id::THETA),PI);
			core::pose::add_variant_type_to_pose_residue(pose,"VIRTUAL_NTERM",2);

			ScoreFunctionOP sf_fa = fa_refine_and_design(pose,NFA);

			//std::cerr << "ref_rep ";
			//for( Size i = 0; i < 4; ++i ) {
			//  std::cerr << pose.energies().residue_total_energies(i*NRES+1)[scoring::fa_rep] << " ";
			//  std::cerr << pose.energies().residue_total_energies(i*NRES+2)[scoring::fa_rep] << " ";
			//}
			//std::cerr << std::endl;

			report( tag, pose, cenpose, sf_fa, oss, censcore );

		}


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
