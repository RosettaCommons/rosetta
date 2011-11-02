#include <basic/database/open.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/smhybrid.OptionKeys.gen.hh>
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
#include <core/init.hh>
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
#include <protocols/basic_moves/FragmentMover.hh>
#include <protocols/electron_density/util.hh>
#include <protocols/flxbb/FlxbbDesign.hh>
#include <protocols/jobdist/standard_mains.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/RigidBodyMover.hh>
#include <protocols/moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/moves/symmetry/SymMinMover.hh>
#include <protocols/moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/symmetric_docking/SymDockingInitialPerturbation.hh>
#include <protocols/symmetric_docking/SymDockingLowRes.hh>
#include <protocols/viewer/viewers.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
// #include <core/scoring/constraints/LocalCoordinateConstraint.hh>
// #include <devel/init.hh>
#include <apps/pilot/will/will_util.hh>

#include <apps/pilot/will/mynamespaces.hh>
using core::conformation::symmetry::SymmData;
using core::conformation::symmetry::SymmDataOP;
using core::conformation::symmetry::SymmetryInfo;
using core::conformation::symmetry::SymmetryInfoOP;
using core::pose::symmetry::make_symmetric_pose;
using protocols::moves::symmetry::SymMinMover;


OPT_KEY( String, hub_ss )
OPT_KEY( IntegerVector, hub_ho_cst )


static core::io::silent::SilentFileData sfd;

void read_fragdata( vector1< core::fragment::FragDataOP > & fds, std::istream & in ) {
	using namespace core::fragment;
	Size n,count=0;
	while( in >> n ) {
	 	string pdb;
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
}

std::map<string, vector1<core::fragment::FragDataOP> >
get_frags_map( ) {
	using namespace core::fragment;
	cout << "reading frags" << std::endl;
	utility::io::izstream in;
	std::map<string,vector1<FragDataOP> > fds;
	basic::database::open(in,"ss_fragfiles/EEE.fragfile"); read_fragdata(fds["EEE"],in); in.close();
	basic::database::open(in,"ss_fragfiles/EEH.fragfile"); read_fragdata(fds["EEH"],in); in.close();
	basic::database::open(in,"ss_fragfiles/EEL.fragfile"); read_fragdata(fds["EEL"],in); in.close();
	basic::database::open(in,"ss_fragfiles/EHH.fragfile"); read_fragdata(fds["EHH"],in); in.close();
	basic::database::open(in,"ss_fragfiles/ELE.fragfile"); read_fragdata(fds["ELE"],in); in.close();
	basic::database::open(in,"ss_fragfiles/ELH.fragfile"); read_fragdata(fds["ELH"],in); in.close();
	basic::database::open(in,"ss_fragfiles/ELL.fragfile"); read_fragdata(fds["ELL"],in); in.close();
	basic::database::open(in,"ss_fragfiles/HEE.fragfile"); read_fragdata(fds["HEE"],in); in.close();
	basic::database::open(in,"ss_fragfiles/HEH.fragfile"); read_fragdata(fds["HEH"],in); in.close();
	basic::database::open(in,"ss_fragfiles/HEL.fragfile"); read_fragdata(fds["HEL"],in); in.close();
	basic::database::open(in,"ss_fragfiles/HHE.fragfile"); read_fragdata(fds["HHE"],in); in.close();
	basic::database::open(in,"ss_fragfiles/HHH.fragfile"); read_fragdata(fds["HHH"],in); in.close();
	basic::database::open(in,"ss_fragfiles/HHL.fragfile"); read_fragdata(fds["HHL"],in); in.close();
	basic::database::open(in,"ss_fragfiles/HLE.fragfile"); read_fragdata(fds["HLE"],in); in.close();
	basic::database::open(in,"ss_fragfiles/HLH.fragfile"); read_fragdata(fds["HLH"],in); in.close();
	basic::database::open(in,"ss_fragfiles/HLL.fragfile"); read_fragdata(fds["HLL"],in); in.close();
	basic::database::open(in,"ss_fragfiles/LEE.fragfile"); read_fragdata(fds["LEE"],in); in.close();
	basic::database::open(in,"ss_fragfiles/LEH.fragfile"); read_fragdata(fds["LEH"],in); in.close();
	basic::database::open(in,"ss_fragfiles/LEL.fragfile"); read_fragdata(fds["LEL"],in); in.close();
	basic::database::open(in,"ss_fragfiles/LHH.fragfile"); read_fragdata(fds["LHH"],in); in.close();
	basic::database::open(in,"ss_fragfiles/LLE.fragfile"); read_fragdata(fds["LLE"],in); in.close();
	basic::database::open(in,"ss_fragfiles/LLH.fragfile"); read_fragdata(fds["LLH"],in); in.close();
	basic::database::open(in,"ss_fragfiles/LLL.fragfile"); read_fragdata(fds["LLL"],in); in.close();
	return fds;
}


core::fragment::FragSetOP make_frag_set(std::string ss, std::map<string, vector1<core::fragment::FragDataOP> > & fds) {
	using namespace core::fragment;
	FragSetOP frags = new ConstantLengthFragSet();
	int const stop = ss.size() + 1 - 3;
	if((int)1 >= stop) return NULL;
	for( Size i = 1; i <= (Size)stop; ++i ) {
		string ss3 = ss.substr(i-1,3);
		bool mkframe = true;
		for(Size j = 0; j < ss3.size(); ++j) if(ss3[j]!='H'&&ss3[j]!='E'&&ss3[j]!='L'&&ss3[j]!='*') mkframe = false;
		if(!mkframe) utility_exit_with_message("oaisnrt");
		FrameOP frame = new Frame(i,3);
		vector1<char> ss0,ss1,ss2;
		if('*'==ss3[0]) { ss0.push_back('H'); ss0.push_back('E'); ss0.push_back('L'); } else ss0.push_back(ss3[0]);
		if('*'==ss3[1]) { ss1.push_back('H'); ss1.push_back('E'); ss1.push_back('L'); } else ss1.push_back(ss3[1]);
		if('*'==ss3[2]) { ss2.push_back('H'); ss2.push_back('E'); ss2.push_back('L'); } else ss2.push_back(ss3[2]);				
		for( Size j = 1; j <= ss0.size(); ++j ) {
		for( Size k = 1; k <= ss1.size(); ++k ) {
		for( Size l = 1; l <= ss2.size(); ++l ) {
			string ss=""; ss+=ss0[j]; ss+=ss1[k]; ss+=ss2[l];
			//cout << "adding ss " << ss << " '" << ss0[j] << "' '" << ss1[k] << "' '" << ss2[l] << "'" << std::endl;
			vector1<FragDataOP>::iterator beg = fds[ss].begin();
			vector1<FragDataOP>::iterator end = fds[ss].end();
			for( vector1<FragDataOP>::iterator fi = beg; fi != end; ++fi ) {
				frame->add_fragment(*fi);
			}
		}}}
		frags->add(frame);
		cout << "make frag " << i << ": " << ss3 << std::endl;
	}
	if(frags->size() == 0) return NULL;
	return frags;
}

core::fragment::FragSetOP make_frag_set_9mers(Size nres) {
	using namespace core::fragment;

	vector1<core::fragment::FragDataOP> fds9;
	std::ifstream in("input/loop_helix.9mers");
	read_fragdata(fds9,in);

	FragSetOP frags = new ConstantLengthFragSet();
	for( Size i = 1; i <= nres-8; ++i ) {
		FrameOP frame = new Frame(i,9);
		for( vector1<FragDataOP>::iterator fi = fds9.begin(); fi != fds9.end(); ++fi ) {
			frame->add_fragment(*fi);
		}
		frags->add(frame);
	}
	if(frags->size() == 0) return NULL;
	return frags;
}



int main(int argc, char *argv[]) {
	core::init(argc,argv);
	using namespace core::scoring;

	string ss = option[OptionKeys::hub_ss]();
	vector1<Size> hocst = option[OptionKeys::hub_ho_cst]();

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
	Size nres = ss.size();
	for(Size ir = 2; ir <= nres; ++ir) {
		core::conformation::ResidueOP tmp = core::conformation::ResidueFactory::create_residue(p.residue(1).residue_type_set().name_map("ALA"));
		tmp->seqpos(ir);
		tmp->chain(1);
		p.append_residue_by_bond( *tmp, true );
		p.set_omega(ir-1,180.0);
		p.set_omega(ir  ,180.0);		
	}
	core::pose::add_upper_terminus_type_to_pose_residue(p,p.n_residue());
	sfsym->show(p);

	make_symmetric_pose(p);
	FoldTree ft = p.fold_tree();
	for(Size i = 1; i <= ft.num_jump(); ++i) cout << i << " " << ft.jump_edge(i) << endl;
	// cout << ft << endl;
	ft.set_jump_atoms( 1, "ORIG", "CA");
	ft.set_jump_atoms( 2, "ORIG", "CA");
	ft.set_jump_atoms( 3, "ORIG", "CA");		
	p.fold_tree(ft);
	// cout << ft << endl;

	p.conformation().detect_bonds();

	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
  	movemap->set_jump(false);
  	movemap->set_bb(true);
  	movemap->set_chi(true);
	protocols::moves::symmetry::SymMinMover( movemap, sfsym, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false ).apply(p);

	core::util::switch_to_residue_type_set(p,"centroid");

	sf3->set_weight(core::scoring::ss_pair,10.0);
	sf3->set_weight(core::scoring::sheet,10.0);

	std::map<string, vector1<core::fragment::FragDataOP> > fds( get_frags_map() );
	core::fragment::FragSetOP frags3 = make_frag_set(ss,fds);
	protocols::moves::MoverOP fragins = new protocols::basic_moves::ClassicFragmentMover(frags3);
	Pose init(p);
	for(int iter = 1; iter < 99999999; ++iter) {
		Pose tmp(init);
		
		for(vector1<Size>::const_iterator i = hocst.begin(); i != hocst.end(); ++i) {
			using namespace core::scoring::constraints;
			Size ir = *(  i);
			Size jr = *(++i);
			tmp.add_constraint( new AtomPairConstraint( AtomID(tmp.residue(ir).atom_index("H"),ir),
				                                        AtomID(tmp.residue(ir).atom_index("O"),jr),	
														new HarmonicFunc(1.8,0.5)));
			tmp.add_constraint( new AngleConstraint( AtomID(tmp.residue(ir).atom_index("N"),ir),	
					                                 AtomID(tmp.residue(ir).atom_index("H"),ir),
				                                     AtomID(tmp.residue(ir).atom_index("O"),jr),	
													 new CircularHarmonicFunc(3.14159,0.5)));
		}
		

		Real temp = 2.0;
		protocols::moves::MonteCarloOP mc = new protocols::moves::MonteCarlo( tmp, *sf3, temp );
		mc->set_autotemp( true, temp );
		mc->set_temperature( temp );
		protocols::moves::RepeatMover( new protocols::moves::TrialMover( fragins, mc ), 10000.0 ).apply( tmp );
		mc->reset( tmp );

		string dumbtag = ss+str(uniform());
		tmp.dump_pdb(dumbtag+".pdb");

		core::io::silent::SilentStructOP ss_out( new core::io::silent::ScoreFileSilentStruct );
	  	ss_out->fill_struct(tmp,dumbtag+".pdb");

	  	// // Write the scorefile
	  	sfd.write_silent_struct( *ss_out, option[OptionKeys::out::file::o]() + "/" + option[ OptionKeys::out::file::silent ]() );


	}

	return 0;
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
	// 		 	//  protocols::moves::symmetry::SymMinMover( movemap, sfsym, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false ).apply(p);
	// 			if(1==ippo) tmp.set_phi  (ir,tmp.phi  (ir)+10.0);
	// 			if(2==ippo) tmp.set_psi  (ir,tmp.psi  (ir)+10.0);
	// 			if(3==ippo) tmp.set_omega(ir,tmp.omega(ir)+10.0);
	// 			core::util::switch_to_residue_type_set(tmp,"fa_standard");
	// 			tmp.dump_pdb("test_"+str(ippo)+"_"+str(ir)+"_"+str(itrial) +".pdb");
	// 		}
	// 	}
	// }
