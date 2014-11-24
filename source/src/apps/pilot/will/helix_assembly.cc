// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:

// includes
	#include <ObjexxFCL/FArray2D.hh>
	#include <ObjexxFCL/format.hh>
	#include <ObjexxFCL/string.functions.hh>
	#include <basic/Tracer.hh>
	#include <basic/database/open.hh>
	#include <basic/options/keys/in.OptionKeys.gen.hh>
	#include <basic/options/keys/out.OptionKeys.gen.hh>
	#include <basic/options/option_macros.hh>
	#include <core/chemical/AtomType.hh>
	#include <core/chemical/ChemicalManager.hh>
	#include <core/conformation/symmetry/util.hh>
	#include <core/import_pose/import_pose.hh>
	#include <core/io/silent/SilentFileData.hh>
	#include <core/pose/PDBInfo.hh>
	#include <core/pose/Pose.hh>
	#include <core/pose/annotated_sequence.hh>
	#include <core/pose/util.hh>
	#include <core/scoring/Energies.hh>
	#include <core/scoring/EnergyGraph.hh>
	#include <core/scoring/ScoreFunctionFactory.hh>
	#include <core/scoring/ScoreTypeManager.hh>
	#include <core/scoring/dssp/Dssp.hh>
	#include <core/scoring/hbonds/HBondOptions.hh>
	#include <core/scoring/methods/EnergyMethodOptions.hh>
	#include <core/scoring/packing/compute_holes_score.hh>
	#include <core/scoring/rms_util.hh>
	#include <core/scoring/sasa.hh>
	#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
	#include <devel/init.hh>
	#include <numeric/conversions.hh>
	#include <numeric/model_quality/rms.hh>
	#include <numeric/random/random.hh>
	#include <numeric/xyz.functions.hh>
	#include <numeric/xyzTransform.hh>
	#include <numeric/xyz.io.hh>
	#include <numeric/xyzVector.hh>
	#include <protocols/sic_dock/RigidScore.hh>
	#include <protocols/sic_dock/SICFast.hh>
	#include <protocols/sic_dock/util.hh>
	#include <protocols/sic_dock/read_biounit.hh>
	#include <utility/io/izstream.hh>
	#include <utility/io/ozstream.hh>
	#include <utility/fixedsizearray1.hh>
	// #include <protocols/sic_dock/designability_score.hh>

	#include <apps/pilot/will/will_util.ihh>

	#ifdef USE_OPENMP
	#include <omp.h>
	#endif

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

////////////////////////////////////////////////

static std::string const chr_chains("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!@#$&.<>?]{}|-_\\~=%zyxwvutsrqponmlkjihgfedcbaABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!@#$&.<>?]{}|-_\\~=%zyxwvutsrqponmlkjihgfedcbaABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!@#$&.<>?]{}|-_\\~=%zyxwvutsrqponmlkjihgfedcbaABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!@#$&.<>?]{}|-_\\~=%zyxwvutsrqponmlkjihgfedcbaABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!@#$&.<>?]{}|-_\\~=%zyxwvutsrqponmlkjihgfedcbaABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!@#$&.<>?]{}|-_\\~=%zyxwvutsrqponmlkjihgfedcbaABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!@#$&.<>?]{}|-_\\~=%zyxwvutsrqponmlkjihgfedcbaABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!@#$&.<>?]{}|-_\\~=%zyxwvutsrqponmlkjihgfedcba");

#define HELIX_ROT 99.3646
#define HELIX_TRANS 1.55642

// types

	typedef numeric::xyzVector<core::Real> Vec;
	typedef numeric::xyzMatrix<core::Real> Mat;
	typedef numeric::xyzTransform<core::Real> Xform;

	using core::id::AtomID;
	using basic::options::option;
	using core::pose::Pose;
	using core::Real;
	using core::scoring::ScoreFunctionOP;
	using core::Size;
	using numeric::max;
	using numeric::min;
	using numeric::random::gaussian;
	using numeric::random::uniform;
	using numeric::rotation_matrix_degrees;
	using numeric::conversions::radians;
	using numeric::conversions::degrees;
	using ObjexxFCL::format::F;
	using ObjexxFCL::format::I;
	using ObjexxFCL::string_of;
	using std::cerr;
	using std::cout;
	using std::endl;
	using std::ostream;
	using std::istream;
	using std::string;
	using utility::io::izstream;
	using utility::io::ozstream;
	using utility::file_basename;
	using utility::vector1;
	using std::endl;
	using core::import_pose::pose_from_pdb;
	using numeric::geometry::hashing::Real3;
	using numeric::geometry::hashing::Real6;

static thread_local basic::Tracer TR( "helix_assembly" );
// options
	// OPT_1GRP_KEY( Integer    , ha, ideal_helix_size                             )
	// OPT_1GRP_KEY( Real       , ha, helix_match_rmsd                             )
	OPT_1GRP_KEY( FileVector , ha, catalog_helices                                 )
	OPT_1GRP_KEY( Integer    , ha, max_res                                         )
	OPT_1GRP_KEY( Integer    , ha, max_res_chain_avg                               )
	OPT_1GRP_KEY( Integer    , ha, helix_max_len                                   )
	OPT_1GRP_KEY( Integer    , ha, helix_min_len                                   )
	OPT_1GRP_KEY( Integer    , ha, tail_max_len                                    )
	OPT_1GRP_KEY( Integer    , ha, tail_max_len_helix                              )
	OPT_1GRP_KEY( Integer    , ha, tail_max_len_sheet                              )
	// OPT_1GRP_KEY( Real       , ha, tail_max_ddg                                )
	// OPT_1GRP_KEY( Real       , ha, tail_max_bsa                                )
	OPT_1GRP_KEY( Integer    , ha, helix_min_offset                                )
	OPT_1GRP_KEY( Integer    , ha, helix_max_offset                                )
	OPT_1GRP_KEY( Integer    , ha, max_num_chains                                  )
	OPT_1GRP_KEY( Real       , ha, tail_penalty                                    )
	OPT_1GRP_KEY( Real       , ha, helix_rms_score_cut                             )
	OPT_1GRP_KEY( Boolean    , ha, debug                                           )
	OPT_1GRP_KEY( Boolean    , ha, pymol_cmd                                       )
	OPT_1GRP_KEY( Real       , ha, redundant_dis_cut                               )
	OPT_1GRP_KEY( Real       , ha, redundant_ang_cut                               )

	void register_options() {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		// reading / merging
		NEW_OPT( ha::catalog_helices , "pdbs to read", ""                             );
		// NEW_OPT( ha::ideal_helix_size , "size of helix align reference", 10        );
		// NEW_OPT( ha::helix_match_rmsd , "size of helix align reference", 0.2       );
		NEW_OPT( ha::max_res , "max total res", 600                                   );
		NEW_OPT( ha::max_res_chain_avg , "max res / nchain", 300                      );
		NEW_OPT( ha::helix_max_len      , "TODO DOC", 30                              );
		NEW_OPT( ha::helix_min_len      , "TODO DOC", 12                              );
		NEW_OPT( ha::tail_max_len       , "TODO DOC", 10                              );
		NEW_OPT( ha::tail_max_len_helix , "TODO DOC",  3                              );
		NEW_OPT( ha::tail_max_len_sheet , "TODO DOC",  2                              );
		// NEW_OPT( ha::tail_max_ddg       , "TODO DOC",  2.0                            );
		// NEW_OPT( ha::tail_max_bsa       , "TODO DOC", 100.0                           );
		NEW_OPT( ha::helix_min_offset   , "TODO DOC", -3                              );
		NEW_OPT( ha::helix_max_offset   , "TODO DOC",  3                              );
		NEW_OPT( ha::max_num_chains     , "TODO DOC",  12                             );
		NEW_OPT( ha::tail_penalty       , "TODO DOC",  0.3                            );
		NEW_OPT( ha::helix_rms_score_cut, "TODO DOC",  0.4                            );
		NEW_OPT( ha::debug , "", false                                                );
		NEW_OPT( ha::pymol_cmd , "", false                                            );

		NEW_OPT( ha::redundant_dis_cut , "",  3.0                                      );
		NEW_OPT( ha::redundant_ang_cut , "", 10.0                                      );

		// int MAX_HELIX_LEN = 30;
		// int MIN_HELIX_LEN = 12;
		// int MAX_TAIL_LEN  = 10;
		// int MAX_TAIL_LEN_HELIX = 3;
		// Real TAIL_PENALTY  = 0.3;
		// Real HELIX_RMS_SCORE_CUT = 0.4; //???????????

	}

void get_pose_chains(Pose const & pose, int nres, vector1<int> & chainlower, vector1<int> & chainupper,/* vector1<char> & pdbchain,*/ bool const /*DEBUG8*/){
	chainlower.push_back(1);
	// pdbchain.push_back(pose.pdb_info()->chain(1));
	if(!pose.residue(1).is_protein()) utility_exit_with_message("pose starts with non-protein residue!");
	int prevchain = pose.chain(1);
	int prevresi = 1;
	// if(DEBUG) cerr<<1 <<" "<< pose.chain(1) <<" "<< pose.pdb_info()->chain(1)<<endl;
	for(int ir = 2; ir <= nres; ++ir){
		if(!pose.residue(ir).is_protein()) continue;
		// if(DEBUG) cerr<<ir <<" "<< pose.chain(ir) /*<<" "<< pose.pdb_info()->chain(ir)*/<<endl;
		if(pose.chain(ir)!=prevchain){
			// if( pose.pdb_info()->chain(ir) == pose.pdb_info()->chain(prevresi) )
			// 	utility_exit_with_message("pdbchain == at boundary!");
			chainlower.push_back(ir);
			chainupper.push_back(prevresi);
			// pdbchain.push_back(pose.pdb_info()->chain(ir));
			prevchain = pose.chain(ir);
		}// else {
		// 	if( pose.pdb_info()->chain(ir) != pose.pdb_info()->chain(prevresi) )
		// 		utility_exit_with_message("pdbchain != not at boundary!");
		// }
		prevresi = ir;
	}
	chainupper.push_back(prevresi);
	TR<<"FOUND_CHAINS";
	for(int i = 1; i <= (int)chainlower.size(); ++i) TR <<" "<</*pdbchain[i]<<" "<<*/chainlower[i]<<" "<<chainupper[i];
	TR<<endl;
}

Real alignhelix(Pose & helix, Pose const & pose, int ir, int jr, int & cen){
	using namespace core::id;
	AtomID_Map<AtomID> alignmap;
	core::pose::initialize_atomid_map(alignmap,helix,core::id::BOGUS_ATOM_ID);
	int helixoffset = ((int)helix.n_residue()-(jr-ir))/2 + 1 - ir;
	int helixcenres = (1+helix.n_residue())/2.0;
	cen = helixcenres-helixoffset;
	for(int k = ir; k <= jr; ++k){
		// TR<<ir <<" "<< jr <<" "<< k <<" "<< helixoffset+k <<" "<< helix.n_residue()<<endl;
		// TR<<helixoffset-ir+k <<" "<< helix.n_residue()<<"    "<<ir+k<< " "<<pose.n_residue()<<endl;
		if(pose.secstruct(k)!='H'){ return 9e9; }
		int natoms = pose.residue(k).aa()==core::chemical::aa_gly ? 4 : 5;
		// assert(lower <= k && k <= upper);
		assert(1 <= helixoffset+k && helixoffset+k <= (int)helix.n_residue());
		for(int ia = 1; ia <= natoms; ++ia ){
			alignmap[AtomID(ia,helixoffset+k)] = AtomID(ia,k);
		}
	}
	core::scoring::superimpose_pose(helix,pose,alignmap); // returns garbage????
	Real rms = 0.0;
	int rmsdatoms = 0;
	for(int k = ir; k <= jr; ++k){
		int natoms = pose.residue(k).aa()==core::chemical::aa_gly ? 4 : 5;
		rmsdatoms += natoms;
		for(int ia = 1; ia <= natoms; ++ia ) {
			rms += helix.xyz(AtomID(ia,helixoffset+k)).distance_squared(pose.xyz(AtomID(ia,k)));
		}
	}
	return sqrt(rms/(Real)rmsdatoms);

}
Real alignhelix(Pose & helix, Pose const & pose, int ir, int jr){
	int tmp;
	return alignhelix(helix,pose,ir,jr,tmp);
}

Vec proj(    Vec u, Vec v){ return numeric::projection_matrix(u)*v; }
// Vec projperp(Vec u, Vec v){	return v - proj(u,v); } // in will_util.hh

string fname2pdbcode(string const & fname){
	std::string pdb = utility::file_basename(fname);
	if( pdb.substr(pdb.size()-3) == ".gz" ) pdb = pdb.substr(0,pdb.size()-3);
	if( pdb.substr(pdb.size()-4) == ".pdb" ) pdb = pdb.substr(0,pdb.size()-4);
	if( pdb.size() != 9 ){
		return pdb;
		utility_exit_with_message("pdb biounit name "+pdb);
		//TR.Warning<<"pdb file does not have a biounit style name! "<<_fname<<endl;
		//pdb_return[1] = '_'; pdb_return[2] = '_'; pdb_return[3] = '_'; pdb_return[4] = '_'; pdb_return[5] = '_';
	}
	if( pdb.substr(4,4) != ".pdb" ) {
		return pdb;
		utility_exit_with_message("pdb biounit name "+pdb);
	}
	int pdb1 = std::atoi(pdb.substr(pdb.size()-1).c_str());
	if(pdb1 > 9) {
		return pdb;
		utility_exit_with_message("bad pdb biounit name "+pdb);
	}
	string pdb_return(pdb.substr(0,5));
	pdb_return[4] = pdb[pdb.size()-1];
	return pdb_return;
}

struct HelixHit{
	string pdb;
	Vec cen,axis,ori;
	int chain,startres,cenres,stopres,chainlen,totallen,nchains,pdbres;
	char pdbchain;
	HelixHit(
		string const & _pdb,
		Vec    const & _cen,
		Vec    const & _axis,
		Vec    const & _ori,
		int    const & _chain,
		int    const & _startres,
		int    const & _cenres,
		int    const & _stopres,
		int const & _chainlen,
		int const & _totallen,
		int const & _nchains,
		int const & _pdbres,
		char const & _pdbchain
	):
		pdb(_pdb),
		cen(_cen),
		axis(_axis),
		ori(_ori),
		chain(_chain),
		startres(_startres),
		cenres(_cenres),
		stopres(_stopres),
		chainlen(_chainlen),
		totallen(_totallen),
		nchains(_nchains),
		pdbres(_pdbres),
		pdbchain(_pdbchain)
	{}
	Xform frame() const {
		Xform x;
		Vec const & X = axis;
		Vec const & Y = ori;
		Vec const Z = X.cross(Y);
		assert(fabs(X.length()-1)<0.00001);
		assert(fabs(Y.length()-1)<0.00001);
		assert(fabs(Z.length()-1)<0.00001);
		x.R.col_x(X).col_y(Y).col_z(Z);
		x.t = cen;
		return x;
	}
	void header(ostream & out){
		out<<"T_H_HIT_HEADER       ID pdb chn  pres  chn Nchn clen nres  stres cenres endres                   CEN_X                   CEN_Y                   CEN_Z                   AXS_X                   AXS_Y                   AXS_Z                   ORI_X                   ORI_Y                   ORI_Z" <<endl;
		//    TERM_HELIX_HIT 3mon1_2N 3mon1 B    18    2    2   50   94     55     62     69       42.88561227252754       16.15025914896387       35.10579309541517      0.1993102565655764     -0.5799247858610918     -0.7899130739339177     -0.3170937507180697     -0.8008850493472015      0.5079711517278093
	}
};
ostream & operator<<(ostream & out, HelixHit const & h){
	using namespace ObjexxFCL::format;
	out<<h.pdb<<" "
	   <<h.pdbchain<<" "
	   <<I(5,h.pdbres)<<" "
	   <<I(4,h.chain)<<" "
	   <<I(4,h.nchains)<<" "
	   <<I(4,h.chainlen)<<" "
	   <<I(4,h.totallen)<<" "
	   <<I(6,h.startres)<<" "
	   <<I(6,h.cenres)<<" "
	   <<I(6,h.stopres)<<" "
	   <<h.cen<<" "
	   <<h.axis<<" "
	   <<h.ori;
	return out;
}
// istream & operator<<(istream & in, HelixHit & h){
// 	in >> h.pdb;
// 	in >> h.chain;
// 	in >> h.startres;
// 	in >> h.stopres;
// 	in >> h.cen .x() >> h.cen .y() >> h.cen .z();
// 	in >> h.axis.x() >> h.axis.y() >> h.axis.z();
// 	in >> h.ori .x() >> h.ori .y() >> h.ori .z();
// 	return in;
// }

int main(int argc, char *argv[]) {
	register_options();
	devel::init(argc,argv);
	using namespace basic::options::OptionKeys;
	using namespace core::scoring;
	using namespace core::id;

	int  const MAX_HELIX_LEN       = option[ha::helix_max_len        ]();
	int  const MIN_HELIX_LEN       = option[ha::helix_min_len        ]();
	int  const MAX_TAIL_LEN        = option[ha::tail_max_len         ]();
	int  const MAX_TAIL_LEN_HELIX  = option[ha::tail_max_len_helix   ]();
	int  const MAX_TAIL_LEN_SHEET  = option[ha::tail_max_len_sheet   ]();
	Real const TAIL_PENALTY        = option[ha::tail_penalty         ]();
	Real const HELIX_RMS_SCORE_CUT = option[ha::helix_rms_score_cut  ]();
	bool const DEBUG               = option[ha::debug                ]();
	bool const PYMOL_CMD           = option[ha::pymol_cmd            ]();
	int  const HELIX_MIN_OFFSET    = option[ha::helix_min_offset     ]();
	int  const HELIX_MAX_OFFSET    = option[ha::helix_max_offset     ]();
	Real const REDUNDANT_DSQ_CUT   = option[ha::redundant_dis_cut]()*option[ha::redundant_dis_cut]();
	Real dotcut = cos(numeric::conversions::radians(option[ha::redundant_ang_cut]()));
	if(option[ha::redundant_ang_cut]()>180.0) dotcut = -1.0001;
	Real const REDUNDANT_DOT_CUT   = dotcut;

	Pose helix; Size idealres=0; Vec idealcen,idealaxs,idealori;
	{
		string seq = ""; for(int i=-14; i<MAX_HELIX_LEN; ++i) seq+="A";
		core::pose::make_pose_from_sequence( helix, seq, "fa_standard", false );
		for(int i = 1; i <= (int)helix.n_residue(); ++i){
			helix.set_phi  (i,-57);
			helix.set_psi  (i,-47);
			helix.set_omega(i,180);
			helix.set_secstruct(i,'H');
		}
		Vec com =  center_of_geom(helix,                  2,helix.n_residue()-1);
		Vec comN = center_of_geom(helix,                  2,                  8);
		Vec comC = center_of_geom(helix,helix.n_residue()-7,helix.n_residue()-1);
		idealcen = com;
		// TR<<helix.n_residue()-7 <<" "<< helix.n_residue()-1<<endl;
		idealaxs = (comC-comN).normalized();
		idealres = (1+helix.n_residue())/2.0;
		idealcen = idealcen+proj(idealaxs,helix.xyz(AtomID(2,idealres))-idealcen);
		idealori = projperp(idealaxs,helix.xyz(AtomID(2,idealres))-idealcen).normalized();
		helix.set_xyz(AtomID(6,1),idealcen    ); // com
		helix.set_xyz(AtomID(7,1),idealcen+idealaxs); // axis

		// cout<<"showvecfrompoint( 50*Vec("+string_of(idealaxs.x())+","+string_of(idealaxs.y())+","+string_of(idealaxs.z())+"),Vec("+string_of(idealcen.x())+","+string_of(idealcen.y())+","+string_of(idealcen.z())+"))"<<endl;
		// cout<<"showvecfrompoint( 10*Vec("+string_of(idealori.x())+","+string_of(idealori.y())+","+string_of(idealori.z())+"),Vec("+string_of(idealcen.x())+","+string_of(idealcen.y())+","+string_of(idealcen.z())+"))"<<endl;
		// helix.dump_pdb("test.pdb");
		// utility_exit_with_message("helix");

		// Real a(0),d(0);
		// for(int i = 8; i <= 35; ++i){
		// 	a += dihedral_degrees(helix.xyz(AtomID(2,i)),com,com+idealaxs,helix.xyz(AtomID(2,i+1)));
		// 	d += proj(idealaxs,helix.xyz(AtomID(2,i+1))-helix.xyz(AtomID(2,i))).length();
		// }
		// cout<<a/28.0<<" "<<d/28.0<<endl;
		// utility_exit_with_message("ast");
		// helix.dump_pdb("ideal_helix.pdb");
	}

	vector1<string> pdbs = option[ha::catalog_helices]();
	for(int ipdb = 1; ipdb <= (int)pdbs.size(); ++ipdb){
		string const & pdb(pdbs[ipdb]);
		cout<<"processing "<<pdb<<endl;
		string pdbcode = fname2pdbcode(pdb);
		Pose pose;
		vector1<Real> bfac,occ;
		vector1<int> pdbres;
		std::map<int,char> pdbchain;
		int nresmodel1=0;//12345;

		bool goodpdb = protocols::sic_dock::read_biounit(pdb,pose,bfac,occ,pdbres,pdbchain,nresmodel1,option[ha::max_res](),DEBUG);
		if(!goodpdb){cerr<<"ERROR reading in "<<pdb<<", skipping"<<endl; continue; }
		core::scoring::dssp::Dssp dssp(pose);
		dssp.insert_ss_into_pose(pose);

		// core::import_pose::pose_from_pdb(pose,pdb);
		// core::scoring::dssp::Dssp dssp(pose);
		// dssp.insert_ss_into_pose(pose);
		// if(DEBUG) pose.dump_pdb("test0.pdb");

		if(DEBUG) cerr<<"main: GET CHAINS"<<endl;
		vector1<int> chainlower,chainupper;
		// vector1<char> pdbchain;
		get_pose_chains(pose,pose.n_residue(),chainlower,chainupper/*,pdbchain*/,DEBUG);
		if((int)pose.n_residue()/(int)chainlower.size() > option[ha::max_res_chain_avg]()) { cout<<"SKIP max_res_chain_avg "<<pdb<<endl; continue; }
		if((int)chainlower.size()>option[ha::max_num_chains]()){ cout<<"SKIP max_num_chains "<<pdb<<endl; continue; }
		vector1<bool> hit_on_chain(chainupper.size(),false);
		Size nhits=0;
		// if(chainlower.size()!=2) { cout<<"SKIP no 2 chains "<<pdb<<endl; continue; }

		vector1<HelixHit> hits;

		for(int ic = 1; ic <= (int)chainlower.size(); ++ic){
			int lower = chainlower[ic];
			int upper = chainupper[ic];
			if(upper-lower+1<10) { cout<<"SKIP has tiny chain "<<pdb<<endl; continue; }

			/// !!!!!!!!!!!!!!!!!!!!!!!!!!
			for(int Nterm_or_Cterm = 1; Nterm_or_Cterm >= 0; --Nterm_or_Cterm){ // ==1 -> Nterm
				if(DEBUG) cout<<"CHECKING CHAIN "<<ic<<(Nterm_or_Cterm?" Nterm ":" Cterm ")<<pose.chain(lower) /*<<" "<< pose.pdb_info()->chain(lower)*/ <<" "<< lower <<" "<< upper<<endl;

				Real bestsc = 9e9;
				int bestir=0,bestjr=0,bestcen=0;
				if(Nterm_or_Cterm){
					for(    int ir = lower             ; ir <= min(upper,lower+MAX_TAIL_LEN ); ++ir ){
						for(int jr = ir+MIN_HELIX_LEN-1; jr <= min(upper,ir+MAX_HELIX_LEN-1); ++jr ){
							if(ir-lower > MAX_TAIL_LEN) continue;
							int nhelixtail=0; for(int tmp=max(ir-7,lower); tmp <= ir-1; ++tmp) if(pose.secstruct(tmp)=='H') ++nhelixtail;
							if( nhelixtail > MAX_TAIL_LEN_HELIX) continue;
							int nsheettail=0;
							for(int tmp=lower; tmp<ir;++tmp){
								if(pose.secstruct(tmp)!='E') continue;
								bool sheet = false;
								for(int tmp2=1; tmp2 <= (int)pose.n_residue(); ++tmp2)
									if(pose.secstruct(tmp2)=='E' && dssp.bb_pair_score(tmp,tmp2)) sheet=true;
								if(sheet) ++nsheettail;
							}
							if( nsheettail > MAX_TAIL_LEN_SHEET ) continue;
							int cen;
							Real rms = alignhelix(helix,pose,ir,jr,cen);
							Real hscore = (rms+Real(ir-lower)/Real(MAX_TAIL_LEN)*TAIL_PENALTY)/sqrt(Real(jr-ir+1)); // add up to
							if(hscore < bestsc){ bestsc = hscore; bestir = ir; bestjr = jr; bestcen = cen; }
							// cout<<"HIT "<<pdbchain[ic]<<(Nterm_or_Cterm?" Nterm ":" Cterm ")<<ir <<" "<< jr <<" "<< rms <<" "<< (Real(ir-lower)/Real(MAX_TAIL_LEN)*TAIL_PENALTY)/sqrt(Real(jr-ir+1))  <<" "<< hscore <<" "<< endl;
						}
					}
				} else {
					for(    int ir = max(lower,upper-MAX_HELIX_LEN-MAX_TAIL_LEN); ir <= upper-MIN_HELIX_LEN+1; ++ir ){
						for(int jr = ir+MIN_HELIX_LEN-1                         ; jr <= upper                ; ++jr ){
							if(upper-jr > MAX_TAIL_LEN) continue;
							int nhelixtail=0; for(int tmp=jr+1; tmp <= min(jr+7,upper);++tmp) if(pose.secstruct(tmp)=='H') ++nhelixtail;
							// cout<<ir <<" "<< jr <<" "<< nhelixtail<<endl;
							if( nhelixtail > MAX_TAIL_LEN_HELIX) continue;
							int nsheettail = 0;
							for(int tmp=jr+1; tmp<=upper;++tmp){
								if(pose.secstruct(tmp)!='E') continue;
								bool sheet = false;
								for(int tmp2=1; tmp2 <= (int)pose.n_residue(); ++tmp2)
									if(	pose.secstruct(tmp)!='E' && dssp.bb_pair_score(tmp,tmp2)) sheet=true;
								if(sheet) ++nsheettail;
							}
							if( nsheettail > MAX_TAIL_LEN_SHEET ) continue;

							int cen;
							Real rms = alignhelix(helix,pose,ir,jr,cen);
							Real hscore = (rms+Real(upper-jr)/Real(MAX_TAIL_LEN)*TAIL_PENALTY)/sqrt(Real(jr-ir+1)); // add up to
							if(hscore < bestsc){ bestsc = hscore; bestir = ir; bestjr = jr; bestcen = cen; }
						}
					}
				}
				if(bestsc > HELIX_RMS_SCORE_CUT) continue;
				/*Real rms =*/ alignhelix(helix,pose,bestir,bestjr);
				Vec cen = helix.xyz(AtomID(6,1));
				Vec axs = (helix.xyz(AtomID(7,1)) - cen).normalized();
				      cen = cen+proj(axs,pose.xyz(AtomID(2,(bestir+bestjr)/2))-cen);
				Vec lbcen = cen+proj(axs,pose.xyz(AtomID(2,bestir))-cen);
				Vec ubcen = cen+proj(axs,pose.xyz(AtomID(2,bestjr))-cen);
				Vec   ori = projperp(axs,pose.xyz(AtomID(2,(bestir+bestjr)/2))-cen).normalized();
				Vec lbori = projperp(axs,pose.xyz(AtomID(2,        bestir   ))-cen).normalized();
				Vec ubori = projperp(axs,pose.xyz(AtomID(2,        bestjr   ))-cen).normalized();
				int chain = pose.chain(lower);
				if( !Nterm_or_Cterm ) chain = -chain;
				hits.push_back(
					HelixHit(
						pdbcode,
						cen,
						axs,
						ori,
						chain,
						bestir,
						bestcen,
						bestjr,
						upper-lower+1,
						pose.n_residue(),
						chainlower.size(),
						pdbres[bestcen],
						pdbchain[pose.chain(bestcen)]
					)
				);
				hit_on_chain[ic] = true;
				++nhits;
				if(DEBUG){
				// cout<<"RMSD "<<bestir <<" "<< bestjr <<" "<< rms <<" " <<Real(bestir-lower)/Real(MAX_subl TAIL_LEN)*TAIL_PENALTY<<" "<< bestsc <<" "<< endl;
					//	string const & _pdb,
					//	Vec    const & _cen,
					//	Vec    const & _axis,
					//	Vec    const & _ori,
					//	int    const & _chain,
					//	int    const & _startres,
					//	int    const & _stopres
					// cout<<"TERM_HELIX_HIT "
					//     <<pose.chain(lower)<<" "//pdbchain[ic]
					//     <<(Nterm_or_Cterm?" Nterm ":" Cterm ")
					//     <<bestir <<" "
					//     <<bestjr <<" "
					//     <<rms <<" "
					//     <<(Real(Nterm_or_Cterm?(bestir-lower):(upper-bestjr))/Real(MAX_TAIL_LEN)*TAIL_PENALTY)/sqrt(Real(bestjr-bestir+1)) <<" "
					//     <<bestsc <<" "
					//     <<lbcen.distance(ubcen)<<" "
					//     <<numeric::dihedral_degrees(lbori,lbcen,ubcen,ubori)<<"     "
					//     <<lbcen<<"     "
					//     <<lbori<<"     "
					//     <<ubcen<<"     "
					//     <<ubori<<" "
					//     <<endl;
					// 	int helixoffset = ((int)helix.n_residue()-(bestjr-bestir)+1)/2 - bestir;
					// 	// cout<<endl<<endl<<endl;
					// 	pymol_commands<<"color orange , test0 and resi "<<            bestir<<"-"<<            bestjr<< endl;
					// 	pymol_commands<<"color yellow, "<<"test1_"+string_of(pose.chain(lower))+(Nterm_or_Cterm?"_Nterm":"_Cterm") <<" and resi "<<helixoffset+bestir<<"-"<<helixoffset+bestjr<< endl;
					// 	// if(Nterm_or_Cterm) pymol_commands<<"showvecfrompoint(-50*Vec("+string_of(axs.x())+","+string_of(axs.y())+","+string_of(axs.z())+"),Vec("+string_of(cen.x())+","+string_of(cen.y())+","+string_of(cen.z())+"))"<<endl;
					// 	// else               pymol_commands<<"showvecfrompoint( 50*Vec("+string_of(axs.x())+","+string_of(axs.y())+","+string_of(axs.z())+"),Vec("+string_of(cen.x())+","+string_of(cen.y())+","+string_of(cen.z())+"))"<<endl;
					// 	pymol_commands<<"showvecfrompoint( 50*Vec("+string_of(axs.x())+","+string_of(axs.y())+","+string_of(axs.z())+"),Vec("+string_of(cen.x())+","+string_of(cen.y())+","+string_of(cen.z())+"))"<<endl;
					// 	pymol_commands<<"showvecfrompoint( 10*Vec("+string_of(ori.x())+","+string_of(ori.y())+","+string_of(ori.z())+"),Vec("+string_of(cen.x())+","+string_of(cen.y())+","+string_of(cen.z())+"))"<<endl;						// cout<<endl<<endl<<endl;
						helix.dump_pdb(utility::file_basename(pdb)+"_"+string_of(pose.chain(lower))+(Nterm_or_Cterm?"_Nterm.pdb":"_Cterm.pdb"));
					// }
				}// else {
				// 	if(DEBUG) cout<<"BAD "<<pdbchain[ic]<<(Nterm_or_Cterm?" Nterm ":" Cterm ")<<bestir <<" "<< bestjr <<" "<< rms <<" " <<Real(bestir-lower)/Real(MAX_TAIL_LEN)*TAIL_PENALTY<<" "<< bestsc <<" "<< endl;
				// }
			}
		}

		// int chains_with_hits=0; for(int ichain=1; ichain<=(int)chainupper.size();++ichain) if(hit_on_chain[ichain]) ++chains_with_hits;
		// if(chains_with_hits < (int)chainlower.size()) { cout <<"SKIP chains_with_hits < nchains"<<endl; continue; }

		// pick longest helix-ish catalog_helices
		// get termini from chains


		if(hits.size()) hits[1].header(cout);
		for(int ihit = 1; ihit <= (int)hits.size(); ++ihit){
			cout<<"TERM_HELIX_HIT "<<hits[ihit].pdb<<"_"<<std::abs(hits[ihit].chain)<<(hits[ihit].chain<0?"C":"N")<<" "<<hits[ihit]<<endl;
		}

		std::map<string,vector1<std::pair<int,Xform> > > seenit;
		vector1<bool> chainpair_redundant((chainlower.size()+2)*(chainlower.size()+2),false);

		vector1<std::pair<int,int> > hitpairs;
		for(int phase = 1; phase <= 4; ++phase){
			for(int ihit = 1; ihit <= (int)hits.size(); ++ihit){
				HelixHit const & a(hits[ihit]);
				if( a.cenres > nresmodel1 ) continue; // at least one of pair in asym unit
				for(int jhit = 1; jhit <= (int)hits.size(); ++jhit){
					if( ihit==jhit ) continue;
					if( std::abs(hits[ihit].cenres-hits[jhit].cenres) < 10 ) continue;
					HelixHit const & b(hits[jhit]);
					if(phase==1 && a.chain>0 && b.chain>0 )	hitpairs.push_back(std::make_pair(ihit,jhit));
					if(phase==2 && a.chain<0 && b.chain<0 )	hitpairs.push_back(std::make_pair(ihit,jhit));
					if(phase==3 && a.chain>0 && b.chain<0 )	hitpairs.push_back(std::make_pair(ihit,jhit));
					if(phase==4 && a.chain<0 && b.chain>0 )	hitpairs.push_back(std::make_pair(ihit,jhit));
				}
			}
		}
		for(int ipair = 1; ipair <= (int)hitpairs.size(); ++ipair){
			int const ihit = hitpairs[ipair].first;
			int const jhit = hitpairs[ipair].second;
			HelixHit const & a(hits[ihit]);
			HelixHit const & b(hits[jhit]);
			int const chain_pair_id = min( std::abs(a.chain),std::abs(b.chain))*chainlower.size()+max( std::abs(a.chain),std::abs(b.chain));
			Xform const xa0 = a.frame();
			Xform const xb0 = b.frame();

			if(DEBUG)
			for(int ii = 1; ii <= (int)hits.size(); ++ii){
				if(hits[ii].cenres>nresmodel1) break;
				for(int jj = 1; jj <= (int)hits.size(); ++jj){
					int tmp = min( std::abs(hits[ii].chain),std::abs(hits[jj].chain))*chainlower.size()+max( std::abs(hits[ii].chain),std::abs(hits[jj].chain));
					cout <<" "<< (chainpair_redundant[tmp] ? string_of(jj) : "-");
				}
				cout << endl;
			}

			// if NC or CN and chainpair has redundency, skip it -- hopefully keep NC/CN hits in asym unit
			if( a.chain*b.chain < 0 && chainpair_redundant[chain_pair_id] ) {
				if(DEBUG) cout << "SKIPHIT " << ihit << " " << jhit << ", inferred redundant" << endl;
				continue;
			}

			for(int Aoffset = 0; Aoffset <= 0; ++Aoffset){
				Real const Ashift = a.chain<0?Aoffset:-Aoffset;
				Xform const xah = Ashift*a.axis*HELIX_TRANS  +  Xform::rot_deg(a.axis,Ashift*HELIX_ROT,a.cen);
				Xform const xa = xah*xa0;
				int const Ares = a.cenres+(int)Ashift;
				int const Aseqres = a.cenres;//max(min(a.stopres-3,Ares),a.startres+3);
				string Aseq = a.chain<0?"C":"N";
					for(int ir=Aseqres-3;ir<=Aseqres+3;++ir) Aseq += pose.residue(ir).name1();

				for(int Boffset = HELIX_MIN_OFFSET; Boffset <= HELIX_MAX_OFFSET; ++Boffset){
					Real const Bshift = b.chain<0?Boffset:-Boffset;
					int Bres = b.cenres+(int)Bshift;
					int const Bseqres = b.cenres;//max(min(b.stopres-3,Bres),b.startres+3);
					string Bseq = b.chain<0?"C":"N";
						for(int ir=Bseqres-3;ir<=Bseqres+3;++ir) Bseq += pose.residue(ir).name1();

					// rotate around helix axis * translate -- xform by Bshift residue helix
					Xform const xbh = Bshift*b.axis*HELIX_TRANS  +  Xform::rot_deg(b.axis,Bshift*HELIX_ROT,b.cen);
					Xform const xb = xbh*xb0;

					Xform const x  = xb * ~xa;
					Real ang; Vec raxis = numeric::rotation_axis(x.R,ang); ang = numeric::conversions::degrees(ang);
					Real talonga = raxis.dot(x.t);
					int nf = 360.0/ang+0.5;

					bool issym = ( fabs(talonga) < 3.0 && fabs(360.0/(Real)nf-ang) < 5.0 && fabs(360.0/ang-(Real)nf) < 0.1 );

					// cout << "DEBUG "
					//      << I(3,ihit)<<" "
					//      << I(3,jhit)<<" "
					//      << I(3,Aoffset)<<" "
					//      << I(3,Boffset)<<" "
					//      << I(5,a.cenres)<<" "
					//      << I(5,b.cenres)<<" "
					//      << I(5,a.stopres-a.startres+1)<<" "
					//      << I(5,b.stopres-b.startres+1)<<" "
					//      << I(5,a.cenres+(int)Ashift)<<" "
					//      << I(5,b.cenres+(int)Bshift)<<" "
					//      << x.t
					//      << F(10,5,ang) << " "
					//      << F(10,5,fabs(talonga)) << " "
					//      << raxis<<" "
					//      << endl;
					// continue;

					Vec rotcen(0,0,0);
					if(issym){
						Vec tmp(a.cen);
						for(int ii = 1; ii <= nf; ++ii){
							tmp = x*tmp;
							rotcen += tmp;
						}
						rotcen /= (Real)nf;
					}
					// cout <<" hist pair "<<I(3,ihit)<<" "<<I(3,jhit)<<" "<<I(5,a.cenres)<<" "<<I(5,b.cenres+(int)Bshift)<<" "<<" Boffset "<<I(3,Boffset)<<" "<<I(3,Bshift)
					//      <<" trans B by " << F(7,3,Bshift*HELIX_TRANS) << " rot "<< F(7,3,Bshift*HELIX_ROT)
					//      <<F(12,6,ang)        <<" "
					//      <<F(12,6,talonga)    <<"     "
					//      << endl;

					cout << (~xa).R*a.ori << " " << x.R*a.ori << " " << b.ori << endl;

					assert(                             (      x.R*a.axis - b.axis).length() < 0.000001);
					assert( Aoffset!=0 || Boffset!=0 || (      x.R*a.ori  - b.ori ).length() < 0.000001);
					assert( Aoffset!=0 || Boffset!=0 || (      x  *a.cen  - b.cen ).length() < 0.000001);
					assert( Aoffset!=0               || (Vec(1,0,0)-(~xa).R*a.axis).length() < 0.000001);
					assert(               Boffset!=0 || (Vec(0,1,0)-(~xa).R*a.ori ).length() < 0.000001);
					assert( Aoffset!=0               || (Vec(1,0,0)-(~xb).R*b.axis).length() < 0.000001);
					assert(               Boffset!=0 || (Vec(0,1,0)-(~xb).R*b.ori ).length() < 0.000001);

					bool redundant = false;
					for(vector1<std::pair<int,Xform> >::const_iterator ix = seenit[Aseq+Bseq].begin(); ix != seenit[Aseq+Bseq].end(); ++ix){
						// cout << "DSQ " << x.t << endl;
						// cout << "DSQ " << ix->second.t    << endl;
						// cout << "DSQ " << ix->second.t.distance_squared(x.t)    << endl;
						// cout << "DOT " << x.R.col_x().dot(ix->second.R.col_x()) << endl;
						// cout << "DOT " << x.R.col_y().dot(ix->second.R.col_y()) << endl;
						// cout << "DOT " << x.R.col_z().dot(ix->second.R.col_z()) << endl;
						if( ix->second.t.distance_squared(x.t)    > REDUNDANT_DSQ_CUT ) continue;
						if( x.R.col_x().dot(ix->second.R.col_x()) < REDUNDANT_DOT_CUT ) continue;
						if( x.R.col_y().dot(ix->second.R.col_y()) < REDUNDANT_DOT_CUT ) continue;
						if( x.R.col_z().dot(ix->second.R.col_z()) < REDUNDANT_DOT_CUT ) continue;
						redundant = true;
						if( std::abs(a.chain)!= std::abs(b.chain)){
							chainpair_redundant[chain_pair_id] = true;
							chainpair_redundant[ix->first    ] = true;
						}
					}

					// cout
					// 	<<hits[ihit].pdb<<"_"<<std::abs(hits[ihit].chain)<<(hits[ihit].chain<0?"C":"N")<<" "
					// 	<<hits[jhit].pdb<<"_"<<std::abs(hits[jhit].chain)<<(hits[jhit].chain<0?"C":"N")<<" "
					// 	<< Aseq+Bseq << " " << Bseq+Aseq << endl;

					if(redundant) {
						if(DEBUG) cout<<"REDUNDANT      "
						<<hits[ihit].pdb<<"_"<<std::abs(hits[ihit].chain)<<(hits[ihit].chain<0?"C":"N")<<" "
						<<hits[jhit].pdb<<"_"<<std::abs(hits[jhit].chain)<<(hits[jhit].chain<0?"C":"N")<<" "
						<< endl;
						continue;
					}

					cout<<"HELIX_JUNCTION "
						<<hits[ihit].pdb<<"_"<<std::abs(hits[ihit].chain)<<(hits[ihit].chain<0?"C":"N")<<" "
						<<hits[jhit].pdb<<"_"<<std::abs(hits[jhit].chain)<<(hits[jhit].chain<0?"C":"N")<<" "
						// <<I(3,ihit)<<" "
						// <<I(3,jhit)<<" "
					    <<I(5,Ares) << " "
					    // <<I(3,a.stopres-a.startres+1) << " "
   					    // <<(a.chain<0?"C":"N")<<(b.chain<0?"C":"N") <<" "
					    <<I(5,Bres) << " "
					    // <<I(3,b.stopres-b.startres+1) << " "
					    // <<I(3,Aoffset)<<" "
					    // <<I(3,Boffset)<<" "
					    <<I(3,Ashift)<<" "
					    <<I(3,Bshift)<<" "
					    <<F(12,6,ang)        <<" "
					    // <<F(12,6,talonga)    <<"     "
					    <<F(12,6,raxis.x())  <<" "
					    <<F(12,6,raxis.y())  <<" "
					    <<F(12,6,raxis.z())  <<" "
					    // <<F(12,6,rotcen.x()) <<" "
					    // <<F(12,6,rotcen.y()) <<" "
					    // <<F(12,6,rotcen.z()) <<"     "
					    <<F(12,6,x.t.x())    <<" "
					    <<F(12,6,x.t.y())    <<" "
					    <<F(12,6,x.t.z())
						// << "                 " << Aseq+Bseq << " " << Bseq+Aseq
					    <<endl;

					// cout << "SEENIT " << Aseq+Bseq << " " <<  x << endl;
					// cout << "SEENIT " << Bseq+Aseq << " " << ~x << endl;
					seenit[Aseq+Bseq].push_back(std::make_pair(chain_pair_id, x));
					seenit[Bseq+Aseq].push_back(std::make_pair(chain_pair_id,~x));

					if(PYMOL_CMD){
						Vec Aaxs = xa.R*Vec(1,0,0);
						Vec Baxs = xb.R*Vec(1,0,0);
						Vec Aori = xa.R*Vec(0,1,0);
						Vec Bori = xb.R*Vec(0,1,0);
						Vec Acen = xa  *Vec(0,0,0);
						Vec Bcen = xb  *Vec(0,0,0);
						cout<<"PYMOL_CMD "<<ihit<<"_"<<jhit<<" showvecfrompoint(Vec("+string_of((Bcen-Acen).x())+","+string_of((Bcen-Acen).y())+","+string_of((Bcen-Acen).z())+"),Vec("+string_of(Acen.x())+","+string_of(Acen.y())+","+string_of(Acen.z())+"))"<<endl;
					    if(issym) cout<<"PYMOL_CMD "<<ihit<<"_"<<jhit<<" a=Vec("+string_of(raxis.x())+","+string_of(raxis.y())+","+string_of(raxis.z())+")"<<endl;
					    if(issym) cout<<"PYMOL_CMD "<<ihit<<"_"<<jhit<<" c=Vec("+string_of(rotcen.x())+","+string_of(rotcen.y())+","+string_of(rotcen.z())+")"<<endl;
						if(issym) cout<<"PYMOL_CMD "<<ihit<<"_"<<jhit<<" showvecfrompoint(100*a,c)"<<endl;
						if(Boffset==HELIX_MIN_OFFSET) cout<<"PYMOL_CMD "<<ihit<<"_"<<jhit<<" showvecfrompoint( 50*Vec("+string_of(Aaxs.x())+","+string_of(Aaxs.y())+","+string_of(Aaxs.z())+"),Vec("+string_of(Acen.x())+","+string_of(Acen.y())+","+string_of(Acen.z())+"))"<<endl;
						if(Boffset==HELIX_MIN_OFFSET) cout<<"PYMOL_CMD "<<ihit<<"_"<<jhit<<" showvecfrompoint( 10*Vec("+string_of(Aori.x())+","+string_of(Aori.y())+","+string_of(Aori.z())+"),Vec("+string_of(Acen.x())+","+string_of(Acen.y())+","+string_of(Acen.z())+"))"<<endl;
						if(Aoffset==HELIX_MIN_OFFSET) cout<<"PYMOL_CMD "<<ihit<<"_"<<jhit<<" showvecfrompoint( 50*Vec("+string_of(Baxs.x())+","+string_of(Baxs.y())+","+string_of(Baxs.z())+"),Vec("+string_of(Bcen.x())+","+string_of(Bcen.y())+","+string_of(Bcen.z())+"))"<<endl;
						if(Aoffset==HELIX_MIN_OFFSET) cout<<"PYMOL_CMD "<<ihit<<"_"<<jhit<<" showvecfrompoint( 10*Vec("+string_of(Bori.x())+","+string_of(Bori.y())+","+string_of(Bori.z())+"),Vec("+string_of(Bcen.x())+","+string_of(Bcen.y())+","+string_of(Bcen.z())+"))"<<endl;
						if(issym) cout<<"PYMOL_CMD "<<ihit<<"_"<<jhit<<" delete tmp*"<<endl;
						if(issym) cout<<"PYMOL_CMD "<<ihit<<"_"<<jhit<<" for i in range("+string_of(nf)+"): cmd.create('tmp%i'%i,'*_ha')"<<endl;
						if(issym) cout<<"PYMOL_CMD "<<ihit<<"_"<<jhit<<" for i in range("+string_of(nf)+"):        rot('tmp%i'%i,a, "+string_of(360.0/(Real)nf)+"*i,c)"<<endl;
						// utility_exit_with_message("rst");
					}
				}
				if(PYMOL_CMD){
					cout<<"PYMOL_CMD "<<ihit<<"_"<<jhit<<" select HH=( (chain "<<chr_chains[std::abs(a.chain)-1]<< " and resi "<<a.startres<<"-"<<a.stopres<<") or (chain "<<chr_chains[std::abs(b.chain)-1]<< " and resi "<<b.startres<<"-"<<b.stopres<<") )"<<endl;
					cout<<"PYMOL_CMD "<<ihit<<"_"<<jhit<<" select HC=( (chain "<<chr_chains[std::abs(a.chain)-1]<< " and resi "<<a.cenres<<") or (chain "<<chr_chains[std::abs(b.chain)-1]<< " and resi "<<b.cenres<<") )"<<endl;
					cout<<"PYMOL_CMD "<<ihit<<"_"<<jhit<<" select none" << endl;
					cout<<"PYMOL_CMD "<<ihit<<"_"<<jhit<<" hide lines; show rib;"<<endl;
					cout<<"PYMOL_CMD "<<ihit<<"_"<<jhit<<" show sti, HH and name n+ca+c+o+cb+h"<<endl;
					cout<<"PYMOL_CMD "<<ihit<<"_"<<jhit<<" show car, HH"<<endl;
					cout<<"PYMOL_CMD "<<ihit<<"_"<<jhit<<" show sph, HC and name CA"<<endl;
					cout<<"PYMOL_CMD "<<ihit<<"_"<<jhit<<" cbow"<<endl;
					cout<<"PYMOL_CMD "<<ihit<<"_"<<jhit<<" util.cnc"<<endl;
					cout<<"PYMOL_CMD "<<ihit<<"_"<<jhit<<" zo vis"<<endl;
				}
			}
		}
		// utility_exit_with_message("TESTING ONE PDB!!!");
		cout << "DONE " << pdb << endl;
	}
	cout << "DONE" << endl;
}





// crap

		// for(    int ir =           1; ir <= (int)pose.n_residue()-refsize-1; ++ir){
		// 	for(int jr =ir+refsize  ; jr <= (int)pose.n_residue()          ; ++jr){
		// 		if( jr-ir > (int)helix.n_residue()-4 ) continue;
		// 		if( jr-ir > 30 ) continue;
		// 		int helixoffset = ((int)helix.n_residue()-ir+jr)/2;
		// 		AtomID_Map<AtomID> alignmap;
		// 		core::pose::initialize_atomid_map(alignmap,helix,core::id::BOGUS_ATOM_ID);
		// 		Real rms = 0.0;
		// 		for(int k = ir; k <= jr; ++k){
		// 			// TR<<helixoffset-ir+k <<" "<< ir+k<<" ";
		// 			// TR<<helixoffset-ir+k <<" "<< helix.n_residue()<<"    "<<ir+k<< " "<<pose.n_residue()<<endl;
		// 			alignmap[        AtomID(1,helixoffset-ir+k)] =                         AtomID(1,k);
		// 			alignmap[        AtomID(2,helixoffset-ir+k)] =                         AtomID(2,k);
		// 			alignmap[        AtomID(3,helixoffset-ir+k)] =                         AtomID(3,k);
		// 			alignmap[        AtomID(4,helixoffset-ir+k)] =                         AtomID(4,k);
		// 		}
		// 		core::scoring::superimpose_pose(helix,pose,alignmap); // returns garbage????
		// 		for(int k = ir; k <= jr; ++k){
		// 			rms += helix.xyz(AtomID(1,helixoffset-ir+k)).distance_squared(pose.xyz(AtomID(1,k)));
		// 			rms += helix.xyz(AtomID(2,helixoffset-ir+k)).distance_squared(pose.xyz(AtomID(2,k)));
		// 			rms += helix.xyz(AtomID(3,helixoffset-ir+k)).distance_squared(pose.xyz(AtomID(3,k)));
		// 			rms += helix.xyz(AtomID(4,helixoffset-ir+k)).distance_squared(pose.xyz(AtomID(4,k)));
		// 		}
		// 		rms = sqrt(rms/(Real)(jr-ir)/4.0);
		// 		cout<<"RMSD "<<(ir+jr)/2 <<" "<< jr-ir <<" "<< rms <<" "<< rms/Real(jr-ir) <<" "<< endl;
		// 		Vec cen = helix.xyz(AtomID(1,6));
		// 		Vec axs = helix.xyz(AtomID(1,7)) - cen;
		// 		// if(uniform() < 0.01){
		// cout<<"showvecfrompoint(Vec("+string_of(axs.x())+","+string_of(axs.y())+","+string_of(axs.z())+"),Vec("+string_of(axs.x())+","+string_of(axs.y())+","+string_of(axs.z())+"))"<<endl;
		// 		// 	pose .dump_pdb("test0.pdb");
		// 		// 	helix.dump_pdb("test1.pdb");
		// 		// 	utility_exit_with_message("arst");
		// 		// }
		// 	}
		// }
