#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/sicdock.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <devel/init.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <utility/io/ozstream.hh>
#include <sstream>
#include <basic/options/option_macros.hh>
#include <utility/io/izstream.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <utility/file/file_sys_util.hh>

// // sicdock stuff
// #include <protocols/sicdock/Assay.hh>
// #include <protocols/sicdock/Bouquet.hh>
// #include <protocols/sicdock/util.hh>
// #include <protocols/sicdock/app/SicDockApp.hh>
// #include <protocols/sicdock/assay/MotifHashAssay.hh>

// motif stuff
#include <core/scoring/motif/motif_hash_stuff.hh>
#include <core/pose/motif/reference_frames.hh>

// rosetta design stuff
#include <core/util/SwitchResidueTypeSet.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/minimization_packing/symmetry/SymPackRotamersMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/minimization_packing/symmetry/SymMinMover.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/kinematics/MoveMap.hh>

static basic::Tracer TR("scheme_score");

bool DUMP = false;

using core::Size;
using core::Real;
using std::cout;
using std::endl;
typedef numeric::xyzTransform<Real> Xform;

OPT_1GRP_KEY( Integer, scheme, nfold )
OPT_1GRP_KEY( Boolean, scheme, debug )
OPT_1GRP_KEY( Boolean, scheme, design )
OPT_1GRP_KEY( Boolean, scheme, logscore )
OPT_1GRP_KEY( Boolean, scheme, overwrite )
OPT_1GRP_KEY( String, scheme, dokfile )
void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( scheme::nfold  , "C syms to test", 1 );
	NEW_OPT( scheme::debug  , "", false );
	NEW_OPT( scheme::design  , "", false );
	NEW_OPT( scheme::logscore  , "", false );
	NEW_OPT( scheme::overwrite  , "", false );
	NEW_OPT( scheme::dokfile  , "", "DEFAULT.dok" );
}


Real compute_bb_motif_score(
	char const & ss1,
	char const & ss2,
	char const & aa1,
	char const & aa2,
	Xform const & x
){
	core::scoring::motif::MotifHashManager & mman = *core::scoring::motif::MotifHashManager::get_instance();
	Real bb_motif = 0;
	// cout << ss1 << " " << ss2 << " " << aa1 << " " << aa2 << endl;
	core::scoring::motif::XformScoreCOP xs_bb_fxn1 = mman.get_xform_score_BB_BB(ss1,ss2,aa1,aa2);
	core::scoring::motif::XformScoreCOP xs_bb_fxn2 = mman.get_xform_score_BB_BB(ss2,ss1,aa2,aa1);
	if ( xs_bb_fxn1 ) bb_motif += xs_bb_fxn1->score_of_bin(x          .rt6());
	if ( xs_bb_fxn2 ) bb_motif += xs_bb_fxn2->score_of_bin(x.inverse().rt6());
	if ( basic::options::option[basic::options::OptionKeys::scheme::logscore]() ) bb_motif = log(bb_motif);
	// cout << bb_motif << endl;
	return bb_motif;
}

Real motif_score_pose(core::pose::Pose const & pose, bool interchain_only=true){
	Size nres1 = pose.size(), nres2=nres1;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		nres1 = core::pose::symmetry::symmetry_info(pose)->num_independent_residues();
		nres2 = core::pose::symmetry::symmetry_info(pose)->num_total_residues_without_pseudo();
	}
	Real score = 0;
	for ( Size ir =    1; ir <= nres1; ++ir ) {
		Xform const xb1 = core::pose::motif::get_backbone_reference_frame(pose,ir);
		for ( Size jr = ir+1; jr <= nres2; ++jr ) {
			if ( interchain_only && pose.chain(ir)==pose.chain(jr) ) continue;
			Xform const xb2 = core::pose::motif::get_backbone_reference_frame(pose,jr);
			if ( xb1.t.distance_squared(xb2.t) > 100.0 ) continue;
			Xform const xbb = xb1.inverse() * xb2;
			Real const symwt = jr > nres1 ? 0.5 : 1.0;
			score += symwt * compute_bb_motif_score( pose.secstruct(ir), pose.secstruct(jr), pose.residue(ir).name1(), pose.residue(jr).name1(), xbb );
		}
	}
	return score;
}

struct ScoreBreakdown {
	Real total;
	Real interchain;
	Real intrachain;
	Real onebody;
	ScoreBreakdown() : total(0),interchain(0),intrachain(0),onebody(0) {}
};
struct DesignOpts {
	bool minimize_bb;
	bool minimize_chi;
	bool minimize_rb;
	bool softrep;
	DesignOpts():
		minimize_bb(false),
		minimize_chi(false),
		minimize_rb(false),
		softrep(false)
	{}
};

void
get_motif_hits_common(
	core::pose::Pose const & pose,
	core::scoring::motif::MotifHits & hits
){
	core::scoring::motif::ResPairMotifQuery query(pose);
	query.clash_check() = true;
	query.interface_only() = false;
	core::scoring::motif::MotifHashManager::get_instance()->get_matching_motifs(query,hits);
}

core::pack::task::PackerTaskOP
make_motif_task(
	core::pose::Pose & pose
){
	using core::pose::Pose;
	core::pack::task::PackerTaskOP task( core::pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_extra_rotamer_flags_from_command_line();
	std::string resfile;
	try{

		// q.match_radius() = motif_match_radius_;
		core::scoring::motif::MotifHits hits;
		get_motif_hits_common(pose,hits);
		if ( hits.size()==0 ) {
			cout << "no motifs!" << endl;
			return NULL;
			utility_exit_with_message("no motifs!!!");
		}
		if ( DUMP ) hits.dump_motifs_pdb("test_motifs.pdb");

		std::set<core::Size> resi_in_resfile;
		resfile = hits.get_resfile(false,resi_in_resfile);
		TR << "RESFILE: " << endl << resfile << endl;
		core::pack::task::parse_resfile_string(pose, *task, resfile);


		utility::vector1<bool> aas(20,false);
		aas[core::chemical::aa_gly] = true;
		// aas[core::chemical::aa_ala] = true;
		// aas[core::chemical::aa_ile] = true;
		// aas[core::chemical::aa_leu] = true;
		// aas[core::chemical::aa_met] = true;
		// aas[core::chemical::aa_val] = true;
		// now make all non-motif residues ala or gly
		for ( Size ir = 1; ir <= pose.size(); ++ir ) {
			// task->nonconst_residue_task(ir).restrict_absent_canonical_aas(aas);
			using namespace core::chemical;
			// if(!core::pose::symmetry::residue_is_independent(pose,ir)) continue; // in in asym unit
			if ( !pose.residue(ir).is_protein() ) continue; // not an aa
			task->nonconst_residue_task(ir).or_include_current(false);
			if ( resi_in_resfile.find(ir)==resi_in_resfile.end() ) {
				// task->nonconst_residue_task(ir).restrict_to_repacking();
				task->nonconst_residue_task(ir).restrict_absent_canonical_aas(aas);
			} else {
				// task->nonconst_residue_task(ir).restrict_absent_canonical_aas(aas);
				// task->nonconst_residue_task(ir).or_ex1_sample_level( static_cast<core::pack::task::ExtraRotSample>(motif_ex1_ );
				// task->nonconst_residue_task(ir).or_ex2_sample_level( static_cast<core::pack::task::ExtraRotSample>(motif_ex2_ );
				// task->nonconst_residue_task(ir).or_ex3_sample_level( static_cast<core::pack::task::ExtraRotSample>(motif_ex3_ );
				// task->nonconst_residue_task(ir).or_ex4_sample_level( static_cast<core::pack::task::ExtraRotSample>(motif_ex4_ );
			}
		}
	} catch(core::pack::task::ResfileReaderException e){
		std::stringstream error_message;
		error_message
			<< "Failed to process resfile" << endl
			<< "RESFILE:" << endl
			<< resfile << endl;
		throw CREATE_EXCEPTION(utility::excn::Exception, error_message.str());
	}
	return task;

}

void
design_pose_motifs_only(
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionCOP sf,
	DesignOpts opts
){
	using namespace core::pack::task;

	PackerTaskOP task = make_motif_task(pose);
	if ( !task ) return;

	if ( core::pose::symmetry::is_symmetric(pose) ) {
		protocols::minimization_packing::symmetry::SymPackRotamersMover(sf,task).apply(pose);
	} else {
		protocols::minimization_packing::PackRotamersMover(sf,task).apply(pose);
	}

	// // pose.dump_pdb("test1.pdb");

	if ( opts.minimize_bb || opts.minimize_chi || opts.minimize_rb ) {
		core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
		movemap->set_jump(opts.minimize_rb);
		movemap->set_bb(opts.minimize_bb);
		movemap->set_chi(opts.minimize_chi);
		// core::pose::symmetry::make_symmetric_movemap( pose, *movemap );
		if ( core::pose::symmetry::is_symmetric(pose) ) {
			protocols::minimization_packing::symmetry::SymMinMover( movemap, sf, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false ).apply(pose);
		} else {
			protocols::minimization_packing::MinMover( movemap, sf, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false ).apply(pose);
		}
	}

}

ScoreBreakdown
extract_scores(
	core::pose::Pose const & pose,
	core::scoring::EnergyMap const & weights
){
	ScoreBreakdown scores;
	using namespace core::scoring;
	Energies    const & energies     ( pose.energies() );
	EnergyGraph const & energy_graph ( energies.energy_graph() );
	Size nres = pose.size();
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		nres = core::pose::symmetry::symmetry_info(pose)->num_independent_residues();
	}
	scores.total = energies.total_energy();
	for ( Size ir = 1; ir <= nres; ++ir ) {
		if ( !pose.residue(ir).is_protein() ) { continue; }
		scores.onebody += energies.onebody_energies(ir).dot(weights);
		for ( utility::graph::Graph::EdgeListConstIter
				iru  = energy_graph.get_node(ir)->const_upper_edge_list_begin(),
				irue = energy_graph.get_node(ir)->const_upper_edge_list_end();
				iru != irue; ++iru
				) {
			EnergyEdge const & edge( static_cast< EnergyEdge const & > (**iru) );
			Size const jr( edge.get_second_node_ind() );
			if ( ir >= jr ) utility_exit_with_message("upper edge fail!");
			if ( pose.chain(ir) == pose.chain(jr) ) {
				scores.intrachain += edge.dot(weights);
			} else {
				scores.interchain += edge.dot(weights);
			}
		}
	}
	return scores;
}

void
centroid_scores_destroys_pose(
	core::pose::Pose & pose,
	ScoreBreakdown & cenvalsc,
	ScoreBreakdown & cenleusc,
	ScoreBreakdown & cencheatsc
){
	core::scoring::ScoreFunctionOP sfcen = core::scoring::ScoreFunctionFactory::create_score_function("cen_std");
	core::util::switch_to_residue_type_set(pose,"centroid");
	if ( DUMP ) pose.dump_pdb("cencheat.pdb");
	sfcen->score(pose);
	cencheatsc = extract_scores(pose,sfcen->weights());
	core::Size nres = pose.size();
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::pose::symmetry::symmetry_info(pose)->num_independent_residues();
	}
	for ( Size ir = 1; ir <= nres; ++ir ) {
		core::pose::replace_pose_residue_copying_existing_coordinates( pose, ir, pose.residue(1).residue_type_set().name_map("VAL") );
	}
	if ( DUMP ) pose.dump_pdb("cenval.pdb");
	sfcen->score(pose);
	cenvalsc = extract_scores(pose,sfcen->weights());
	for ( Size ir = 1; ir <= nres; ++ir ) {
		core::pose::replace_pose_residue_copying_existing_coordinates( pose, ir, pose.residue(1).residue_type_set().name_map("LEU") );
	}
	if ( DUMP ) pose.dump_pdb("cenleu.pdb");
	sfcen->score(pose);
	cenleusc = extract_scores(pose,sfcen->weights());
}

void scheme_score(){
	using namespace basic::options;
	using namespace OptionKeys;
	using core::id::AtomID;

	// int nfold = option[sicdock::nfold]();
	DUMP = option[scheme::debug]();

	utility::io::ozstream dokout(option[basic::options::OptionKeys::scheme::dokfile]());

	core::scoring::ScoreFunctionOP sf = core::scoring::get_score_function(true);
	cout << "score functions FA:  " <<  core::scoring::get_score_functionName(true) << endl;
	core::scoring::methods::EnergyMethodOptions myopt = sf->energy_method_options();
	myopt.hbond_options().decompose_bb_hb_into_pair_energies(true);
	sf->set_energy_method_options(myopt);

	// using namespace protocols::sicdock::scores;
	// using namespace protocols::sicdock;
	// Assay assay;
	// MotifHashAssayOP mhscore_ = new MotifHashAssay;
	// assay.add_method_twobody(mhscore_);

	utility::vector1<std::string> filenames = option[in::file::s]();
	if ( option[in::file::l].user() ) {
		cout << "using list at -l, appending to -s" << endl;
		runtime_assert(option[in::file::l]().size()==1);
		utility::io::izstream infile(option[in::file::l]()[1]);
		std::string tmp;
		while ( infile>>tmp ) filenames.push_back(tmp);
		infile.close();
	}
	cout << "scoring total of " << filenames.size() << " files" << endl;

	// Size count = 0;
	for ( std::string const & path : filenames ) {
		// cout << "progress: " << ((Real)++count)/(Real)filenames.size() << endl;


		// pose->dump_pdb("testala.pdb");

		Real scheme_score = 0; {
			// Bouquet b;
			// Pose tmp;
			// core::pose::PoseOP sicpose = core::import_pose::pose_from_file(path, core::import_pose::PDB_file);
			// for(unsigned ichain = 1; ichain <= sicpose->conformation().num_chains(); ++ichain){
			//  core::pose::extract_pose_chain(*sicpose,ichain,tmp);
			//  std::string rosepath = path;
			//  if( sicpose->conformation().num_chains() > 1 ){
			//   rosepath += ",chain="+ObjexxFCL::string_of(ichain);
			//  }
			//  Rose r(rosepath);
			//  if(nfold==1) b.add_rose( r );
			//  else   b.add_rose( r, new FixedAxisOper(nfold) );
			// }
			// b.arrange_symmetry();
			// // b.dump_pdb("test_b.pdb");
			// scheme_score = assay.score(b);
		}

		core::pose::PoseOP pose = core::import_pose::pose_from_file(path, core::import_pose::PDB_file);
		if ( basic::options::option[basic::options::OptionKeys::symmetry::symmetry_definition].user() ) {
			core::pose::symmetry::make_symmetric_pose(*pose,"");
		}
		core::scoring::dssp::Dssp(*pose).insert_ss_into_pose(*pose,false);
		scheme_score = motif_score_pose(*pose,true);

		// std::cout << "PING /Users/sheffler/rosetta/scheme/source/src/apps/pilot/will/scheme_score.cc/scheme_score.cc:276" << std::endl;
		if ( option[scheme::design]() ) {
			bool exists = utility::file::file_exists(path+"_scheme_score.pdb.gz");
			if ( !exists || option[scheme::overwrite]() ) {
				if ( DUMP ) pose->dump_pdb("test0.pdb");
				DesignOpts opts;
				opts.minimize_chi = true;
				design_pose_motifs_only(*pose,sf,opts);
				// dump asym unit
				core::pose::Pose asymunit;
				core::pose::symmetry::extract_asymmetric_unit(*pose,asymunit,false);
				asymunit.dump_pdb(path+"_scheme_score.pdb.gz");
				if ( DUMP ) pose->dump_pdb("test2.pdb");
			} else {
				cout << "design file exists: " << path+"_scheme_score.pdb.gz" << endl;
				continue;
			}
		}

		// sf->set_weight(core::scoring::scheme_bb_bb,1.0);
		// sf->score(*pose);
		// cout << "scheme_sc_sc " << pose->energies().total_energies()[core::scoring::scheme_sc_sc] << endl;
		// cout << "scheme_sc_bb " << pose->energies().total_energies()[core::scoring::scheme_sc_bb] << endl;
		// cout << "scheme_sc_ph " << pose->energies().total_energies()[core::scoring::scheme_sc_ph] << endl;
		// cout << "scheme_sc_po " << pose->energies().total_energies()[core::scoring::scheme_sc_po] << endl;
		// cout << "scheme_bb_bb " << pose->energies().total_energies()[core::scoring::scheme_bb_bb] << endl;
		// cout << "scheme_bb_ph " << pose->energies().total_energies()[core::scoring::scheme_bb_ph] << endl;
		// cout << "scheme_bb_po " << pose->energies().total_energies()[core::scoring::scheme_bb_po] << endl;
		// cout << "scheme_ph_po_short " << pose->energies().total_energies()[core::scoring::scheme_ph_po_short] << endl;
		// cout << "scheme_ph_po_long  " << pose->energies().total_energies()[core::scoring::scheme_ph_po_long] << endl;
		// sf->set_weight(core::scoring::scheme_bb_bb,0.0);
		// core::scoring::motif::MotifHits hits;
		// get_motif_hits_common(*pose,hits);
		// pose->dump_pdb("test.pdb");
		// hits.dump_motifs_pdb("test_motifs.pdb");

		// sf->score(*pose);
		// ScoreBreakdown fasc = extract_scores(*pose,sf->weights());

		// ScoreBreakdown cenvalsc,cenleusc,cencheatsc;
		// centroid_scores_destroys_pose(*pose,cenvalsc,cenleusc,cencheatsc);

		cout << "SCHEME_SCORE "
			<< path << " "
			// << fasc.interchain << " "
			// << fasc.total << " "
			<< scheme_score << " "
			// << cenvalsc.interchain << " "
			// << cenvalsc.total << " "
			// << cenleusc.interchain << " "
			// << cenleusc.total << " "
			// << cencheatsc.interchain << " "
			// << cencheatsc.total << " "
			<< endl;

		// assay.reset();
		// std::string fname = utility::file_basename(path);
		// protocols::sicdock::app::generic_print_results( std::cout, dokout, fname, assay, mhscore_ );
	}
}


int main(int argc, char *argv[]) {
	try {
		register_options();
		devel::init(argc,argv);
		scheme_score();
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}
}



