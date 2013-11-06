#include <basic/database/open.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/smhybrid.OptionKeys.gen.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/util.hh>
#include <basic/Tracer.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/conformation/symmetry/VirtualCoordinate.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
#include <devel/init.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/optimizeH.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/annotated_sequence.hh>b
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/CircularHarmonicFunc.hh>
#include <core/scoring/constraints/MultiConstraint.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/XYZ_Func.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/electron_density/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/packstat/compute_sasa.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/abinitio/AbrelaxMover.hh>
#include <protocols/abinitio/ClassicAbinitio.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/electron_density/util.hh>
#include <protocols/flxbb/FlxbbDesign.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/simple_moves/symmetry/SymDockingInitialPerturbation.hh>
#include <protocols/symmetric_docking/SymDockingLowRes.hh>
#include <protocols/viewer/viewers.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
// #include <core/scoring/constraints/LocalCoordinateConstraint.hh>
// #include <devel/init.hh>
#include <apps/pilot/will/will_util.ihh>

#include <apps/pilot/will/mynamespaces.ihh>
using core::conformation::symmetry::SymmData;
using core::conformation::symmetry::SymmDataOP;
using core::conformation::symmetry::SymmetryInfo;
using core::conformation::symmetry::SymmetryInfoOP;
using core::pose::symmetry::make_symmetric_pose;
using protocols::simple_moves::symmetry::SymMinMover;


OPT_KEY( String, hub_sequence )

void register_options() {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  NEW_OPT( hub_sequence ,"the alignment file ", "alignment.hhpred");
}



static core::io::silent::SilentFileData sfd;


// struct HubAbinitioMover : protocols::moves::Mover {
// 	core::pose::Pose init;
// 	ScoreFunctionOP sf0,sf1,sf2,sf3,sf5,sfsym,sfasym;
// 	protocols::moves::MoverOP fragmove9,fragmove3;
// 	HubAbinitioMover() {
// 		using namespace core::scoring;
// 		string seq = option[OptionKeys::hub_sequence]();
// 		if( seq[0] != 'Z' ) utility_exit_with_message("first residue must be Z!!");

// 		sfsym  = getScoreFunction();
// 		sfasym = core::scoring::symmetry::asymmetrize_scorefunction(sfsym);
// 		sf0 = new symmetry::SymmetricScoreFunction(ScoreFunctionFactory::create_score_function("score0"));
// 		sf1 = new symmetry::SymmetricScoreFunction(ScoreFunctionFactory::create_score_function("score1"));
// 		sf2 = new symmetry::SymmetricScoreFunction(ScoreFunctionFactory::create_score_function("score2"));
// 		sf3 = new symmetry::SymmetricScoreFunction(ScoreFunctionFactory::create_score_function("score3"));
// 		sf5 = new symmetry::SymmetricScoreFunction(ScoreFunctionFactory::create_score_function("score5"));
// 		sf1->set_weight(core::scoring::atom_pair_constraint,1.0);
// 		sf1->set_weight(core::scoring::    angle_constraint,1.0);
// 		sf2->set_weight(core::scoring::atom_pair_constraint,1.0);
// 		sf2->set_weight(core::scoring::    angle_constraint,1.0);
// 		sf5->set_weight(core::scoring::atom_pair_constraint,1.0);
// 		sf5->set_weight(core::scoring::    angle_constraint,1.0);
// 		sf3->set_weight(core::scoring::atom_pair_constraint,1.0);
// 		sf3->set_weight(core::scoring::    angle_constraint,1.0);

// 		core::chemical::ResidueTypeSetCAP rtsfa = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
// 		rtsfa->name_map("CHC");

// 		Pose p( *core::import_pose::pose_from_pdb(*rtsfa,"input/CHC_HUB_FA.pdb") );
// 		core::pose::remove_upper_terminus_type_from_pose_residue(p,1);
// 		Size nres = seq.size();
// 		for(Size ir = 2; ir <= nres; ++ir) {
// 			core::conformation::ResidueOP tmp = core::conformation::ResidueFactory::create_residue(
// 				*p.residue(1).residue_type_set().aa_map(core::chemical::aa_from_oneletter_code(seq[ir-1]))[1] );
// 			tmp->seqpos(ir);
// 			tmp->chain(1);
// 			p.append_residue_by_bond( *tmp, true );
// 			p.set_omega(ir-1,180.0);
// 			p.set_omega(ir  ,180.0);
// 		}
// 		core::pose::add_upper_terminus_type_to_pose_residue(p,p.n_residue());
// 		//sfsym->show(p);

// 		make_symmetric_pose(p);
// 		p.conformation().detect_bonds();
// 		FoldTree ft = p.fold_tree();
// 		//for(Size i = 1; i <= ft.num_jump(); ++i) cout << i << " " << ft.jump_edge(i) << endl;
// 		// cout << ft << endl;
// 		ft.set_jump_atoms( 1, "ORIG", "CA");
// 		ft.set_jump_atoms( 2, "ORIG", "CA");
// 		ft.set_jump_atoms( 3, "ORIG", "CA");
// 		p.fold_tree(ft);
// 		// cout << ft << endl;

// 		core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
// 	  	movemap->set_jump(false);
// 	  	movemap->set_bb(true);
// 	  	movemap->set_chi(true);
// 	  	core::pose::symmetry::make_symmetric_movemap(p,*movemap);
// 	//	protocols::simple_moves::symmetry::SymMinMover( movemap, sfsym, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false ).apply(p);


// 		// get frags
// 		core::fragment::FragSetOP frags9 = core::fragment::FragmentIO(
// 			option[ OptionKeys::abinitio::number_9mer_frags ](),
// 			option[ OptionKeys::frags::nr_large_copies ](),
// 			option[ OptionKeys::frags::annotate ]()
// 		).read_data( option[ OptionKeys::in::file::frag9 ] );
// 		core::fragment::FragSetOP frags3 = core::fragment::FragmentIO(
// 			option[ OptionKeys::abinitio::number_3mer_frags ],
// 			1, //nr_copies
// 			option[ OptionKeys::frags::annotate ]
// 		).read_data( option[ OptionKeys::in::file::frag3 ] );
// 		fragmove9 = new protocols::simple_moves::ClassicFragmentMover(frags9,movemap);
// 		fragmove3 = new protocols::simple_moves::ClassicFragmentMover(frags3,movemap);


// 		init = p;
// 		core::util::switch_to_residue_type_set(init,"centroid");

// 	}

// 	void apply(core::pose::Pose & p) {

// 		core::util::switch_to_residue_type_set(init,"centroid");
// 		//cout << iter << " fold" << endl;
// 		p = init;
// 		{
// 			Real temp = 2.0;
// 			protocols::moves::MonteCarloOP mc = new protocols::moves::MonteCarlo( p, *sf0, temp );
// 			mc->set_autotemp( true, temp );
// 			mc->set_temperature( temp );
// 			protocols::moves::RepeatMover( new protocols::moves::TrialMover( fragmove9, mc ), 2000.0 ).apply( p );
// 			mc->reset( p );
// 		}{
// 			Real temp = 2.0;
// 			protocols::moves::MonteCarloOP mc = new protocols::moves::MonteCarlo( p, *sf1, temp );
// 			mc->set_autotemp( true, temp );
// 			mc->set_temperature( temp );
// 			protocols::moves::RepeatMover( new protocols::moves::TrialMover( fragmove9, mc ), 2000.0 ).apply( p );
// 			mc->reset( p );
// 		}{
// 			Real temp = 2.0;
// 			protocols::moves::MonteCarloOP mc = new protocols::moves::MonteCarlo( p, *sf2, temp );
// 			mc->set_autotemp( true, temp );
// 			mc->set_temperature( temp );
// 			protocols::moves::RepeatMover( new protocols::moves::TrialMover( fragmove9, mc ), 2000.0 ).apply( p );
// 			mc->reset( p );
// 		}{
// 			Real temp = 2.0;
// 			protocols::moves::MonteCarloOP mc = new protocols::moves::MonteCarlo( p, *sf5, temp );
// 			mc->set_autotemp( true, temp );
// 			mc->set_temperature( temp );
// 			protocols::moves::RepeatMover( new protocols::moves::TrialMover( fragmove9, mc ), 2000.0 ).apply( p );
// 			mc->reset( p );
// 		}{
// 			Real temp = 2.0;
// 			protocols::moves::MonteCarloOP mc = new protocols::moves::MonteCarlo( p, *sf3, temp );
// 			mc->set_autotemp( true, temp );
// 			mc->set_temperature( temp );
// 			protocols::moves::RepeatMover( new protocols::moves::TrialMover( fragmove3, mc ), 2000.0 ).apply( p );
// 			mc->reset( p );
// 		}

// 		// protocols::abinitio::ClassicAbinitio abi(frags3,frags9,movemap);
// 		// abi.init(p);
// 		// abi.apply(p);

// 		// p.dump_pdb("prerelax.pdb");
// 		protocols::relax::FastRelax rlx(sfsym);
// 		rlx.apply(p);

// 	}
// 	std::string get_name() const { return "HubAbinitioMover"; }
// };

Vec center_ca(Pose const & p, Size st = 1, Size nres = 0) {
	if(0==nres) nres = p.n_residue();
	Vec com(0,0,0);
	Size n = 0;
	for(Size i = st; i <= nres; ++i) {
		if(p.residue(i).has("CA")) {
			com += p.residue(i).xyz("CA");
			n++;
		}
	}
	return com / (Real)n;
}

Vec center_heavy(Pose const & p, Size st = 1, Size nres = 0) {
	if(0==nres) nres = p.n_residue();
	Vec com(0,0,0);
	Size n = 0;
	for(Size i = st; i <= nres; ++i) {
		for(Size ia = 1; ia <= p.residue(i).nheavyatoms(); ++ia) {
			com += p.xyz(AtomID(ia,i));
			n++;
		}
	}
	return com / (Real)n;
}

void calc_c3_rmsd(Size const nres, Pose p, Pose const & native, Vec const & natcom, Vec const & natcomca, Real & carmsd, Real & aarmsd) {
	if(native.n_residue() == 0) return;
	if(native.n_residue()+1 != nres) {
		cout << p.sequence() << endl;
		cout << " " << native.sequence() << endl;
		utility_exit_with_message("native should be 1 chain w/o hub");
	}
	core::pose::symmetry::make_asymmetric_pose(p);
	Real aarmsd1 = 0.0; {
		Vec com   = center_heavy(p,2,nres);
		Real ang   = dihedral_degrees(com  ,Vec(0,0,0),Vec(0,0,1),natcom  );
		trans_pose(p,Vec(0,0,natcom.z()-com.z()));
		rot_pose(p,Vec(0,0,1),ang);
		Size naa = 0;
		for(Size ir = 1; ir < nres; ++ir) {
			if( native.residue(ir).name() != p.residue(ir+1).name() ||
			    native.residue(ir).nheavyatoms() != p.residue(ir+1).nheavyatoms()) {
				native.dump_pdb("debug_native.pdb");
				p.dump_pdb("debug_p.pdb");
				cout << "nat " << ir   << " " << native.residue(ir  ).nheavyatoms() << " " << native.residue(ir  ).name()
				     << "  p " << ir+1 << " " <<      p.residue(ir+1).nheavyatoms() << " " <<      p.residue(ir+1).name() << endl;
				utility_exit_with_message("mismatched native residue "+str(ir));
			}
			for(Size ia = 1; ia <= native.residue(ir).nheavyatoms(); ++ia) {
				aarmsd1 += native.xyz(AtomID(ia,ir)).distance_squared(p.xyz(AtomID(ia,ir+1)));
				naa++;
			}
		}
		aarmsd1 = sqrt(aarmsd1/(Real)naa);
	}
	Real carmsd1 = 0.0; {
		Vec com   = center_ca(p,2,nres);
		Real ang  = dihedral_degrees(com  ,Vec(0,0,0),Vec(0,0,1),natcomca  );
		trans_pose(p,Vec(0,0,natcomca.z()-com.z()));
		rot_pose(p,Vec(0,0,1),ang);
		Size naa = 0;
		for(Size ir = 1; ir < nres; ++ir) {
			carmsd1 += native.xyz(AtomID(2,ir)).distance_squared(p.xyz(AtomID(2,ir+1)));
			naa++;
		}
		carmsd1 = sqrt(carmsd1/(Real)naa);
	}
	rot_pose(p,Vec(1,0,0),180.0,Vec(0,0,0));
	Real aarmsd2 = 0.0; {
		Vec com   = center_heavy(p,2,nres);
		Real ang   = dihedral_degrees(com  ,Vec(0,0,0),Vec(0,0,1),natcom  );
		trans_pose(p,Vec(0,0,natcom.z()-com.z()));
		rot_pose(p,Vec(0,0,1),ang);
		Size naa = 0;
		for(Size ir = 1; ir < nres; ++ir) {
			if(native.residue(ir).nheavyatoms() != p.residue(ir+1).nheavyatoms()) utility_exit_with_message("mismatched native residue "+str(ir));
			for(Size ia = 1; ia <= native.residue(ir).nheavyatoms(); ++ia) {
				aarmsd2 += native.xyz(AtomID(ia,ir)).distance_squared(p.xyz(AtomID(ia,ir+1)));
				naa++;
			}
		}
		aarmsd2 = sqrt(aarmsd2/(Real)naa);
	}
	Real carmsd2 = 0.0; {
		Vec com   = center_ca(p,2,nres);
		Real ang  = dihedral_degrees(com  ,Vec(0,0,0),Vec(0,0,1),natcomca  );
		trans_pose(p,Vec(0,0,natcomca.z()-com.z()));
		rot_pose(p,Vec(0,0,1),ang);
		Size naa = 0;
		for(Size ir = 1; ir < nres; ++ir) {
			carmsd2 += native.xyz(AtomID(2,ir)).distance_squared(p.xyz(AtomID(2,ir+1)));
			naa++;
		}
		carmsd2 = sqrt(carmsd2/(Real)naa);
	}
	if(carmsd1 < carmsd2) {
		carmsd = carmsd1;
		aarmsd = aarmsd1;
	} else {
		carmsd = carmsd2;
		aarmsd = aarmsd2;
	}
}

int main(int argc, char *argv[]) {

	try {

	register_options();
	devel::init(argc,argv);
	using namespace core::scoring;

	string seq = option[OptionKeys::hub_sequence]();
	if( seq[0] != 'Z' ) utility_exit_with_message("first residue must be Z!!");

	Pose native;
	string tagpref = "";
	Vec natcomca,natcom;
	if(option[OptionKeys::in::file::native].user()) {
		native = *core::import_pose::pose_from_pdb(option[OptionKeys::in::file::native]());
		core::pose::remove_lower_terminus_type_from_pose_residue(native,1);
		tagpref = utility::file_basename(option[OptionKeys::in::file::native]());
		tagpref = tagpref.substr(0,tagpref.size()-4);
		natcomca = center_ca   (native);
		natcom   = center_heavy(native);
		Pose tmp(native);
		core::pose::symmetry::make_symmetric_pose(tmp);
		tmp.dump_pdb(option[OptionKeys::out::file::o]() + "/native_trimer.pdb");
	}

	ScoreFunctionOP sfsym  = getScoreFunction();
	ScoreFunctionOP sfasym = new ScoreFunction(*sfsym);

	ScoreFunctionOP sf0 = new symmetry::SymmetricScoreFunction(ScoreFunctionFactory::create_score_function("score0"));
	ScoreFunctionOP sf1 = new symmetry::SymmetricScoreFunction(ScoreFunctionFactory::create_score_function("score1"));
	ScoreFunctionOP sf2 = new symmetry::SymmetricScoreFunction(ScoreFunctionFactory::create_score_function("score2"));
	ScoreFunctionOP sf3 = new symmetry::SymmetricScoreFunction(ScoreFunctionFactory::create_score_function("score3"));
	ScoreFunctionOP sf5 = new symmetry::SymmetricScoreFunction(ScoreFunctionFactory::create_score_function("score5"));
	sf1->set_weight(core::scoring::atom_pair_constraint,1.0);
	sf1->set_weight(core::scoring::    angle_constraint,1.0);
	sf2->set_weight(core::scoring::atom_pair_constraint,1.0);
	sf2->set_weight(core::scoring::    angle_constraint,1.0);
	sf5->set_weight(core::scoring::atom_pair_constraint,1.0);
	sf5->set_weight(core::scoring::    angle_constraint,1.0);
	sf3->set_weight(core::scoring::atom_pair_constraint,1.0);
	sf3->set_weight(core::scoring::    angle_constraint,1.0);

	core::chemical::ResidueTypeSetCAP rtsfa = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
	rtsfa->name_map("CHC");

	Pose p( *core::import_pose::pose_from_pdb(*rtsfa,"input/CHC_HUB_FA.pdb") );
	core::pose::remove_upper_terminus_type_from_pose_residue(p,1);
	Size nres = seq.size();
	for(Size ir = 2; ir <= nres; ++ir) {
		core::conformation::ResidueOP tmp = core::conformation::ResidueFactory::create_residue(
			*p.residue(1).residue_type_set().aa_map(core::chemical::aa_from_oneletter_code(seq[ir-1]))[1] );
		tmp->seqpos(ir);
		tmp->chain(1);
		p.append_residue_by_bond( *tmp, true );
		p.set_omega(ir-1,180.0);
		p.set_omega(ir  ,180.0);
	}
	core::pose::add_upper_terminus_type_to_pose_residue(p,p.n_residue());
	//sfsym->show(p);

	make_symmetric_pose(p);
	p.conformation().detect_bonds();
	FoldTree ft = p.fold_tree();
//	for(Size i = 1; i <= ft.num_jump(); ++i) cout << i << " " << ft.jump_edge(i) << endl;
	// cout << ft << endl;
	ft.set_jump_atoms( 1, "ORIG", "CA");
	ft.set_jump_atoms( 2, "ORIG", "CA");
	ft.set_jump_atoms( 3, "ORIG", "CA");
	p.fold_tree(ft);
	// cout << ft << endl;

	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
  	movemap->set_jump(false);
  	movemap->set_bb(true);
  	movemap->set_chi(true);
  	core::pose::symmetry::make_symmetric_movemap(p,*movemap);
//	protocols::simple_moves::symmetry::SymMinMover( movemap, sfsym, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false ).apply(p);


	// get frags
	core::fragment::FragSetOP frags9 = core::fragment::FragmentIO(
		option[ OptionKeys::abinitio::number_9mer_frags ](),
		option[ OptionKeys::frags::nr_large_copies ](),
		option[ OptionKeys::frags::annotate ]()
	).read_data( option[ OptionKeys::in::file::frag9 ] );
	core::fragment::FragSetOP frags3 = core::fragment::FragmentIO(
		option[ OptionKeys::abinitio::number_3mer_frags ],
		1, //nr_copies
		option[ OptionKeys::frags::annotate ]
	).read_data( option[ OptionKeys::in::file::frag3 ] );
	protocols::moves::MoverOP fragmove9 = new protocols::simple_moves::ClassicFragmentMover(frags9,movemap);
	protocols::moves::MoverOP fragmove3 = new protocols::simple_moves::ClassicFragmentMover(frags3,movemap);


	Pose init(p);
	core::util::switch_to_residue_type_set(init,"centroid");
	for(Size iter = 1; iter <= option[OptionKeys::out::nstruct](); ++iter) {
		// cout << iter << " fold" << endl;
		p = init;
		{
			Real temp = 2.0;
			protocols::moves::MonteCarloOP mc = new protocols::moves::MonteCarlo( p, *sf0, temp );
			mc->set_autotemp( true, temp );
			mc->set_temperature( temp );
			protocols::moves::RepeatMover( new protocols::moves::TrialMover( fragmove9, mc ), 2000.0 ).apply( p );
			mc->reset( p );
		}{
			Real temp = 2.0;
			protocols::moves::MonteCarloOP mc = new protocols::moves::MonteCarlo( p, *sf1, temp );
			mc->set_autotemp( true, temp );
			mc->set_temperature( temp );
			protocols::moves::RepeatMover( new protocols::moves::TrialMover( fragmove9, mc ), 2000.0 ).apply( p );
			mc->reset( p );
		}{
			Real temp = 2.0;
			protocols::moves::MonteCarloOP mc = new protocols::moves::MonteCarlo( p, *sf2, temp );
			mc->set_autotemp( true, temp );
			mc->set_temperature( temp );
			protocols::moves::RepeatMover( new protocols::moves::TrialMover( fragmove9, mc ), 2000.0 ).apply( p );
			mc->reset( p );
		}{
			Real temp = 2.0;
			protocols::moves::MonteCarloOP mc = new protocols::moves::MonteCarlo( p, *sf5, temp );
			mc->set_autotemp( true, temp );
			mc->set_temperature( temp );
			protocols::moves::RepeatMover( new protocols::moves::TrialMover( fragmove9, mc ), 2000.0 ).apply( p );
			mc->reset( p );
		}{
			Real temp = 2.0;
			protocols::moves::MonteCarloOP mc = new protocols::moves::MonteCarlo( p, *sf3, temp );
			mc->set_autotemp( true, temp );
			mc->set_temperature( temp );
			protocols::moves::RepeatMover( new protocols::moves::TrialMover( fragmove3, mc ), 2000.0 ).apply( p );
			mc->reset( p );
		}

		// protocols::abinitio::ClassicAbinitio abi(frags3,frags9,movemap);
		// abi.init(p);
		// abi.apply(p);

		// p.dump_pdb("prerelax.pdb");
		// cout << iter << " relax" << endl;
		protocols::relax::FastRelax rlx(sfsym);
		rlx.apply(p);
		//p.dump_pdb("relaxed_"+lzs(iter,7)+".pdb");

		string fn = option[OptionKeys::out::file::o]() + "/" + tagpref +"_"+ str(uniform()).substr(2,8) + ".pdb.gz";

		Real carmsd = -1.0, aarmsd = -1.0;
		calc_c3_rmsd(nres,p,native,natcom,natcomca,carmsd,aarmsd);

		core::io::silent::SilentStructOP ss_out( new core::io::silent::ScoreFileSilentStruct );
	  	ss_out->fill_struct(p,fn);
	  	ss_out->add_energy( "carmsd", carmsd );
	  	ss_out->add_energy( "aarmsd", aarmsd );
	  	sfd.write_silent_struct( *ss_out, option[OptionKeys::out::file::o]() + "/" + option[ OptionKeys::out::file::silent ]() );
		p.dump_pdb(fn);

	}

	return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

}









	// for(Size ippo = 1; ippo <= 3; ++ippo ) {
	// 	for(Size ir = 1; ir <= nres; ++ir ) {
	// 		Pose tmp(p);
	// 	  	for(Size itrial = 1; itrial <= 10; ++itrial) {
	// 			// for(Size i = 1; i <= nres; ++i) {
	// 			// 	p.set_phi(1,uniform()*360.0);
	// 			// 	p.set_psi(1,uniform()*360.0);
	// 			// 	p.set_omega(1,170.0+uniform()*20.0);
	// 			// }
	// 		 	//  protocols::simple_moves::symmetry::SymMinMover( movemap, sfsym, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false ).apply(p);
	// 			if(1==ippo) tmp.set_phi  (ir,tmp.phi  (ir)+10.0);
	// 			if(2==ippo) tmp.set_psi  (ir,tmp.psi  (ir)+10.0);
	// 			if(3==ippo) tmp.set_omega(ir,tmp.omega(ir)+10.0);
	// 			core::util::switch_to_residue_type_set(tmp,"fa_standard");
	// 			tmp.dump_pdb("test_"+str(ippo)+"_"+str(ir)+"_"+str(itrial) +".pdb");
	// 		}
	// 	}
	// }
