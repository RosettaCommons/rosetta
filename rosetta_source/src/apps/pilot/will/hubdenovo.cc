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
#include <core/init.hh>
#include <core/import_pose/import_pose.hh>
#include <core/id/NamedAtomID.hh>
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
#include <core/pose/annotated_sequence.hh>
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
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/electron_density/util.hh>
#include <protocols/flxbb/FlxbbDesign.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/jobdist/standard_mains.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/moves/TrialMover.hh>
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
#include <apps/pilot/will/frag_util.hh>

using core::conformation::symmetry::SymmData;
using core::conformation::symmetry::SymmDataOP;
using core::conformation::symmetry::SymmetryInfo;
using core::conformation::symmetry::SymmetryInfoOP;
using core::pose::symmetry::make_symmetric_pose;
using protocols::simple_moves::symmetry::SymMinMover;
using protocols::moves::MoverOP;
using core::scoring::constraints::ConstraintOP;
using core::id::NamedAtomID;
typedef utility::vector1<core::Size> Sizes;
typedef utility::vector1<std::string> Strings;

OPT_KEY( File  , hub_pdb )
OPT_KEY( String, hub_ss )
OPT_KEY( String, hub_sequence )
OPT_KEY( File, hub_cst_cfg )
OPT_KEY( IntegerVector, hub_ho_cst )
OPT_KEY( Real, hub_cen_energy_cut )
OPT_KEY( Boolean, hub_no_fa )
OPT_KEY( Boolean, hub_graphics )

void register_options() {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  NEW_OPT( hub_pdb      ,"", "input/CHC_HUB_FA.pdb" );
  NEW_OPT( hub_ss       ,"", "" );
  NEW_OPT( hub_sequence ,"", "" );
  NEW_OPT( hub_cst_cfg  ,"", "" );
  NEW_OPT( hub_ho_cst   ,"", utility::vector1< Size >() );
  NEW_OPT( hub_cen_energy_cut   ,"", 200.0 );
  NEW_OPT( hub_no_fa   ,"", false );
  NEW_OPT( hub_graphics,"", false );
}

static core::io::silent::SilentFileData sfd;

static basic::Tracer tr("hubdenovo");


/*
Real calc_c3_carmsd(Size const nres, Pose p, Pose const & native) {
	if(native.n_residue() == 0) return;
	if(native.n_residue()+1 != nres) {
		tr << p.sequence() << endl;
		tr << " " << native.sequence() << endl;
		utility_exit_with_message("native should be 1 chain w/o hub");
	}
	core::pose::symmetry::make_asymmetric_pose(p);
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
	return min(carmsd1,carmsd2);
}
*/



// bool hocstcmp (std::pair<Size,Size> i, std::pair<Size,Size> j) {
// 	return abs(i.first-i.second) < abs(j.first-j.second);
// }
string peek(std::istream & in) {
	std::streampos i = in.tellg();
	string r;
	in >> r;
	in.seekg(i);
	return r;
}
template<typename T>
bool read_ignore_comments(std::istream & in, T & s) {
	char buf[9999];
	while(peek(in)[0] == '#') {
		in.getline(buf,9999);
	}
	return (in >> s);
}
string read_ignore_comments(std::istream & in) {
	char buf[9999];
	while(peek(in)[0] == '#') {
		in.getline(buf,9999);
	}
	string s;
	if(!(in >> s)) utility_exit_with_message("unexpected EOF");
	return s;
}
Strings getline(std::istream & in) {
	string buf;
	while(buf.size()==0 && in.good()) {
		std::getline(in,buf);
		buf = buf.substr(0,buf.find('#'));
	}
	std::istringstream iss(buf);
	Strings v;
	string s;
	while(iss >> s) v.push_back(s);
	return v;
}
struct DCST {
	// sc cst
	string atm1,atm2;    	Size rsd1,rsd2;    	Real d,sd;   bool active;
	DCST() : atm1(""),atm2(""),rsd1(0),rsd2(0),d(0),sd(0),active(false) {}
	DCST( string _atm1, Size _rsd1, string _atm2, Size _rsd2, Real _d, Real _sd ) : atm1(_atm1), rsd1(_rsd1), atm2(_atm2), rsd2(_rsd2), d(_d), sd(_sd),active(false) {}
};
struct CST {
	// bb cst
	int tplt,grp,dres1,tres1,dres2,tres2;     bool active;
	CST() : tplt(0),grp(0),dres1(0),tres1(0),dres2(0),tres2(0),active(false) {}
	CST(int _tplt, int _grp, int a, int b, int c, int d) : tplt(_tplt),grp(_grp),dres1(a),tres1(b),dres2(c),tres2(d),active(false) {}
};
std::ostream & operator<<(std::ostream & out, CST const & c) {
	out << c.tplt << " " << c.grp << " " << c.dres1 << " " << c.tres1 << " " << c.dres2 << " " << c.tres2 << endl;
	return out;
}
typedef utility::vector1<CST> CSTs;

// fixed SS
struct ConstraintConfig {
	Size nres,nsub,nhub;
	Real CSTSDMULT;
	string fname;
	vector1<char> ss;
	std::map<string,Sizes> ssmap;
	std::map<Size,Strings> seq;
	Strings templates_fname;
	vector1<PoseOP> templates_fa,templates_cen;
	vector1<std::map<Size,Size> > template_map;
	vector1<Size> template_sc,template_sc_resi;
	vector1<CST> cst_sc,cst_bb;
	vector1<DCST> dcst;
	Size ncst,nintracst,nintercst;
	core::chemical::ResidueTypeSetCAP crs,frs;
	virtual ~ConstraintConfig() {}
	ConstraintConfig()
	: nres(0),nsub(0),nhub(0),CSTSDMULT(1),fname("NONE"),ss(0),templates_fname(0),templates_fa(0),templates_cen(0),template_sc(0),template_sc_resi(0),cst_sc(0),cst_bb(0),dcst(0),crs(NULL),frs(NULL)
	{
		init();
	}
	ConstraintConfig(string cfgfile)
 	: nres(0),nsub(0),nhub(0),CSTSDMULT(1),fname("NONE"),ss(0),templates_fname(0),templates_fa(0),templates_cen(0),template_sc(0),template_sc_resi(0),cst_sc(0),cst_bb(0),dcst(0),crs(NULL),frs(NULL)
	{
		tr << "reading config file: " << cfgfile << endl;
		fname = cfgfile;
		init();
		parse_config_file(cfgfile);
		show_cst_grids(tr);
		show_sequence(tr);
	}
	void init() {
		frs = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
		crs = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
	}
	std::string ssstr() const {
		string s = "";
		for(Size i = 1; i <= nres; ++i) s += ss[i];
		return s;
	}
	void parse_config_file(string infname) {
		utility::io::izstream in(infname);
		if(!in.good()) utility_exit_with_message("err opening "+infname);
		parse_config_file(in);
	}
	Strings residue_names_from_sequence(string s) const {
		core::chemical::ResidueTypeCAPs r = core::pose::residue_types_from_sequence(s,*frs,false);
		Strings rn;
		for(Size i = 1; i <= r.size(); ++i) {
			rn.push_back( str( r[i]->name1() ) );
		}
		return rn;
	}
	void show_sequence(std::ostream & out = tr) const {
		for(Size i = 1; i <= nres; ++i) {
			tr << I(4,i) << " " << ss[i];
			if(seq.count(i)==1) {
				for(Size j = 1; j <= seq.find(i)->second.size(); ++j) tr << " " << seq.find(i)->second[j];
			} else {
				tr << " default";
			}
			tr << endl;
		}
	}
	void show_cst_grids(std::ostream & out = tr) const {
		string alpha = "0abcdefghijklmnop";
		string ALPHA = "0ABCDEFGHIJKLMNOP";
		vector1<ObjexxFCL::FArray2D<char> > grids;
		for(Size i = 1; i <= templates_fname.size(); ++i) {
			ObjexxFCL::FArray2D<char> grid(templates_fa[i]->n_residue(),ss.size(),'.');
			int X=0,Y=0;
			for(CSTs::const_iterator ic = cst_bb.begin(); ic != cst_bb.end(); ++ic) {
				if(ic->tplt==i) {
					grid(ic->tres1,ic->dres2) = alpha[ic->grp];
					//grid(ic->tres2,ic->dres1) = alpha[ics];
					grid(ic->tres1,ic->dres1) = '-';
					grid(ic->tres2,ic->dres2) = '-';
					X = max(X,ic->tres1+3);
					X = max(X,ic->tres2+3);
					Y = max(Y,ic->dres1+3);
					Y = max(Y,ic->dres2+3);
				}
			}
			for(CSTs::const_iterator ic = cst_sc.begin(); ic != cst_sc.end(); ++ic) {
				if(ic->tplt==i) {
					//grid(ic->tres1,ic->dres2) = str(ics)[0];
					grid(ic->tres2,ic->dres1) = str(ic->grp)[0];
					grid(ic->tres1,ic->dres1) = '+';
					grid(ic->tres2,ic->dres2) = '+';
					X = max(X,ic->tres1+3);
					X = max(X,ic->tres2+3);
					Y = max(Y,ic->dres1+3);
					Y = max(Y,ic->dres2+3);
				}
			}
			grids.push_back(grid);
			X = min(X,(int)(templates_cen[i]->n_residue()));
			Y = min(Y,(int)(nsub*nres));

			out << templates_fname[i] << endl << "        ";
			for(Size it = 1; it <= X; ++it) if(it%2==1) out << lss(it,4);
			out << endl << "          ";
			for(Size it = 1; it <= X; ++it) if(it%2==0) out << lss(it,4);
			out << endl << "          ";
			for(Size it = 1; it <= X; ++it) out << " " << templates_fa[i]->residue(it).name1(); out << endl;
			for(Size id = 1; id <= Y; ++id) {
				if((id-1)%nres==0) out << "--- ";
				else               out << "    ";
				out << lss(id,3) << " " << ss[id] << " " ;
				char c = '.';
				for(Size it = 1; it <= X; ++it) out << " " << grid(it,id); out << endl;
			}
		}
	}
	Sizes parse_residues(string s, bool aliases=true) const {
		Sizes r;
		if(aliases && ssmap.count(s)>0) {
			r = ssmap.find(s)->second;
		} else {
			bool isnum=true,isdash=true;
			for(Size i = 0; i < s.size(); ++i) {
				isnum  &= (s[i] >= '0' && s[i] <= '9');
				isdash &= (s[i] >= '0' && s[i] <= '9' || s[i]=='-');
			}
			if(isnum) {
				r.push_back(atoi(s.c_str()));
			} else if(isdash) {
				if(s.find('-')!=s.rfind('-')) utility_exit_with_message("one dash only!");
				string s1=s.substr(0,s.find('-'));
				string s2=s.substr(s.find('-')+1);
				Size lb = atoi(s1.c_str()), ub = atoi(s2.c_str());
				for(Size i = lb; i <= ub; ++i) r.push_back(i);
			} else {
				utility_exit_with_message("parse_residues "+s);
			}
		}
		// tr << "parse_residues " << s << " " << r << endl;
		return r;
	}
	Sizes parse_residues(Strings s, bool aliases=true) const {
		Sizes res;
		for(Strings::const_iterator i = s.begin(); i != s.end(); ++i) {
			Sizes r = parse_residues(*i,aliases);
			res.insert(res.end(),r.begin(),r.end());
		}
		//tr << "parse_residues '" << s << "'' " << res << endl;
		return res;
	}
	Size resi_to_primary(Size resi) const {
		if(resi < 1) utility_exit_with_message("aroist");
		if(resi > nsub*nres) utility_exit_with_message("aroistarst");
		return (resi-1) % nres + 1;
	}
	vector1<Strings > parse_sequence(string s) const {
		vector1<Strings > r;
		r.push_back(residue_names_from_sequence(s));
		return r;
	}

	void add_sym_cst(Pose & p, AtomID id1, AtomID id2, Real d, Real sd) const {
//		if(nsub > 11) utility_exit_with_message("residue mapping unclear!!!!!!!!");
		//tr << "add sym cst " << id1 << " " << id2 << endl;
		using namespace core::scoring::constraints;
		if(id1.rsd() > id2.rsd()) {
			AtomID tmp(id1);
			id1 = id2;
			id2 = tmp;
		}
		if(id2.rsd() > nsub*nres) utility_exit_with_message("2nd constraint rsd "+str(id2.rsd())+" outside of nres*nsub");
		if(id1.rsd() > nres     ) return;//utility_exit_with_message("1st constraint rsd "+str(id1.rsd())+" outside of primary subunit");
		//tr << "SYMCST " << id1.rsd() << "-" << id2.rsd() << endl;
		p.add_constraint( new AtomPairConstraint( id1 , id2 , new HarmonicFunc(d,sd) ) );
		int sub2 = (id2.rsd()-1)/nres + 1;
		if(sub2 > 1 && sub2 <= nhub) {
			AtomID id1B( id2.atomno(), id2.rsd() - nres * (sub2-1)             );
			AtomID id2B( id1.atomno(), id1.rsd() - nres * (sub2-1) + nhub*nres );
			//tr << "SYMCST " << id1.rsd() << "-" << id2.rsd() << " " << id1B.rsd() << "-" << id2B.rsd() << endl;
			p.add_constraint( new AtomPairConstraint( id1B, id2B, new HarmonicFunc(d,sd) ) );
		}
		if(sub2 > nhub) { // !!!!!!!!!!!!!! assuming dimer cst on higher sym!
			AtomID id1B( id2.atomno(), id2.rsd() - nres * (sub2-1) );
			AtomID id2B( id1.atomno(), id1.rsd() + nres * (sub2-1) );
			//tr << "SYMCST " << id1.rsd() << "-" << id2.rsd() << " " << id1B.rsd() << "-" << id2B.rsd() << endl;
			p.add_constraint( new AtomPairConstraint( id1B, id2B, new HarmonicFunc(d,sd) ) );
		}
	}
	int hub_seq_sep(int r1, int r2) const {
		if(r1 >  r2) { int tmp = r1; r1 = r2; r2 = tmp;	}		
		if(r1 >  nres) return nsub*nres;
		if(r1 <= nres) return r2-r1;
		if(r2 <= nhub*nres) return (r2-1)%nres+1 + (r1-1)%nres+1;
		if(r2 >  nhub*nres) return nhub*nres+1;
	}
	void parse_config_file(std::istream & in) {
		string op,op2,val;
		while(true) {
			tr << "ROOT peak: " << (char)in.peek() << endl;
			if(!(read_ignore_comments(in,op))) break;
			if("SUBUNITS"==op) {
				if(!read_ignore_comments(in,nhub)) utility_exit_with_message("error in SUBUNITS HERE");
				if(!read_ignore_comments(in,nsub)) utility_exit_with_message("error in SUBUNITS "+str(nhub)+" HERE");
				tr << "SUBUNITS: " << nsub << endl;
			} else
			if("CSTSDMULT"==op) {
				if(!read_ignore_comments(in,CSTSDMULT)) utility_exit_with_message("error in CSTSDMULT HERE");
				tr << "CSTSDMULT: " << CSTSDMULT << endl;
			} else
			if("SECSTRUCT"==op) {
				if(!nsub) utility_exit_with_message("SUBUNITS not set before SECSTRUCT");
				string ssstr;
				read_ignore_comments(in,ssstr);
				for(Size is = 0; is < ssstr.size(); ++is) {
					if(ssstr[is]!='H'&&ssstr[is]!='E'&&ssstr[is]!='L'&&ssstr[is]!='_') {
						utility_exit_with_message("parse_config_file: bad SECSTRUCT: "+ssstr.substr(is,1)+", "+ssstr);
					}
				}
				nres = ssstr.size();
				template_sc.resize(nres*nsub);
				template_sc_resi.resize(nres*nsub);
				for(Size i = 1; i <= nsub; ++i) {
					for(Size j = 0; j < ssstr.size(); ++j) ss.push_back(ssstr[j]);
				}
				tr << "SECSTRUCT " << ss.size() << endl;
				{ // make ssmap
					std::map<char,Size> nss; nss['H']=0; nss['E']=0; nss['L']=0; nss['_']=0;
					char prev = '@';
					Size ssst = 0;
					string symbol;
					string alphabet = "0ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijkhijklmnopqrstuvwxyz";
					for(Size i = 1; i <= ss.size(); ++i) {
						if(prev != ss[i]) {
							if(ssst>0) { // add whole previous elment
								Sizes v;
								for(Size is = ssst; is < i; ++is) v.push_back(is);
								ssmap[symbol] = v;
							}
							prev = ss[i];
							nss[ss[i]]++;
							symbol = ss[i] + alphabet.substr(nss[ss[i]],1);
							ssst = i;
						}
						ssmap[symbol+str(i-ssst+1)] = Sizes(1,i);
					}
					Sizes v;
					for(Size is = ssst; is <= ss.size(); ++is) v.push_back(is);
					ssmap[symbol] = v;
					for(std::map<string,Sizes >::const_iterator i = ssmap.begin(); i != ssmap.end(); ++i) {
						// tr << "ssmap['" << i->first << "'] = ";
						char ssc = i->first[0];
						for(Sizes::const_iterator iv = i->second.begin(); iv != i->second.end(); ++iv) {
							// tr << " " << *iv;
							if(ss[*iv] != ssc) utility_exit_with_message("SECSTRUCT parse issue: "+i->first+" "+str(*iv));
						}
						// tr << endl;
					}
					vector1<std::pair<string,Sizes> > toadd;
					for(std::map<string,Sizes >::const_iterator i = ssmap.begin(); i != ssmap.end(); ++i) {
						if(i->second.size() == 0) utility_exit_with_message("ssmap of 0 size! "+i->first);
						if(i->second.size() == 1) continue;
						int fr = i->second.front();
						int bk = i->second.back();
						for(int j = 1; j < 10; ++j) {
							string symbol = i->first + str(-j);
							if( fr - j >= 1         ) toadd.push_back(std::pair<string,Sizes>(symbol,Sizes(1,fr-j)));
							symbol = i->first + "+" + str( j);
							if( bk + j <= nres*nsub ) toadd.push_back(std::pair<string,Sizes>(symbol,Sizes(1,bk+j)));
						}
					}
					for(Size i = 1; i <= toadd.size(); ++i) ssmap[toadd[i].first] = toadd[i].second;
					// for(std::map<string,Sizes >::CONST_iterator i = ssmap.begin(); i != ssmap.end(); ++i) {
					// 	tr << "SSMAP " << i->first << " " << i->second.front() << "-" << i->second.back() << endl;
					// }
					// utility_exit_with_message("debug");
				}
			} else
			if("SEQUENCE"==op) {
				if(ssmap.size()==0) utility_exit_with_message("SECSTRUCT must come before SEQUENCE");
				while(true) {
					string buf;
					while(peek(in).size() && peek(in)[0]=='#') std::getline(in,buf);
					if("TEMPLATE" ==peek(in)) break;
					if("SECSTRUCT"==peek(in)) break;
					if("CONSTRAINT"==peek(in)) break;
					if("SEQUENCE" ==peek(in)) break;
					if("SUBUNITS" ==peek(in)) utility_exit_with_message("BAD9679867");
					if(!(read_ignore_comments(in,op2))) break;
					string tmp = read_ignore_comments(in);
					Sizes r = parse_residues(op2);
					vector1<Strings > s = parse_sequence(tmp);
					if(s.size()==1 && r.size()==s[1].size()) {
						Strings tmp = s[1];
						s.clear();
						for(Size i = 1; i <= tmp.size(); ++i) s.push_back(Strings(1,tmp[i]));
					}
					if(r.size()!=s.size()) utility_exit_with_message("SEQUENCE size mismatch: "+op2+" vs. "+tmp);
					for(Size i = 1; i <= r.size(); ++i) {
						if(r[i]==1) utility_exit_with_message("cannot mutate resi 1!!! "+op2+" "+tmp);
						if(r[i] > nres) utility_exit_with_message("AAs can only be specified in primary subunit "+str(r[i])+" > "+str(nres)+", "+op2+" "+tmp);
						if(seq.count(r[i])) utility_exit_with_message("multiple AA sets at position "+str(r[i])+", "+op2+" "+tmp);
						seq[r[i]] = s[i];
					}
					// tr << "sequence: " << op2 << " " << tmp << endl;
					// for(Size i = 1; i <= r.size(); ++i) {
					// 	tr << "   " << r[i];
					// 	for(Size j = 1; j <= s[i].size(); ++j) tr << " " << s[i][j];
					// 	tr << endl;
					// }
				}
				tr << "SEQUENCE" << endl;
			} else
			if("CONSTRAINT"==op) {
				if(ssmap.size()==0) utility_exit_with_message("SECSTRUCT must come before CONSTRAINT");
				string atm1,rsd1,atm2,rsd2;
				Real d,sd;
				in >> atm1 >> rsd1 >> atm2 >> rsd2 >> d >> sd;
				Sizes tmp1 = parse_residues(rsd1);
				if( tmp1.size() != 1 ) utility_exit_with_message("bad resi for DCST: "+rsd1);
				Sizes tmp2 = parse_residues(rsd2);
				if( tmp2.size() != 1 ) utility_exit_with_message("bad resi for DCST: "+rsd2);
				dcst.push_back(DCST(atm1,tmp1[1],atm2,tmp2[1],d,sd));
			} else
			if("TEMPLATE"==op) {
				string tpdb,rsd;
				int scgrp=0,bbgrp=0,exgrp=0;
				char buf[9999];
				read_ignore_comments(in,tpdb);
				if(tpdb.substr(tpdb.size()-4,4)!=".pdb"&&tpdb.substr(tpdb.size()-7,7)!=".pdb.gz") utility_exit_with_message("bad pdb: "+tpdb);
				templates_fname.push_back(tpdb);
				templates_fa .push_back( pose_from_pdb(*frs,tpdb) );
				templates_cen.push_back( pose_from_pdb(*crs,tpdb) );
				Size wcount = 0;
				while(true) {
					if("TEMPLATE"==peek(in)) break;
					if(!(read_ignore_comments(in,op2))) break;
					while(true) {
						++wcount;
						if(wcount > 100) utility_exit_with_message("TEMPLATE "+tpdb+" parse error on "+op2);
						if("MAPPING"==peek(in)||"SIDECHAINS"==peek(in)||"BACKBONE"==peek(in)||"DISTANCES"==peek(in)||"TEMPLATE"==peek(in)||""==peek(in)) break;
						if("MAPPING"==op2) {
							tr << "MAPPING" << endl;
							Sizes dres = parse_residues(getline(in));
							Sizes tres = parse_residues(getline(in));
							if(dres.size()!=tres.size()) utility_exit_with_message("MAPPING res no. mismatch in template "+str(tpdb)+" MAPPING");
							std::map<Size,Size> m;
							for(Size i = 1; i <= dres.size(); ++i) m[tres[i]] = dres[i];
							template_map.push_back(m);
							tr << "MAPPING {t:d} " << tpdb << " " << m << endl;
							break;
						} else {
							if(template_map.size()!=templates_fa.size()) utility_exit_with_message("MAPPING must be defined first in TEMPLATE");
							if("SIDECHAIN"==op2) {
								Sizes tres = parse_residues(getline(in),false);
								Sizes dres;
								for(Sizes::const_iterator i = tres.begin(); i != tres.end(); ++i) {
									Size tresi = *i;
									if(template_map.back().count(tresi)==0) utility_exit_with_message("SIDECHAIN: template "+tpdb+" residue "+str(tresi)+" not mapped to design residue!");
									Size dresi = template_map.back()[tresi];
									dres.push_back(dresi);
									if(template_sc[resi_to_primary(dresi)]!=0 && template_sc[resi_to_primary(dresi)]!=templates_fname.size())
										utility_exit_with_message("template "+tpdb+": design AA "+str(dresi)+" already mapped to template "+templates_fname.back());
									template_sc[resi_to_primary(dresi)] = templates_fname.size();
									template_sc[dresi] = templates_fname.size();
									template_sc_resi[dresi] = tresi;
									if(seq.count(dresi)) {
										if(seq[dresi].size()!=1 || seq[dresi][1]!=str(templates_cen.back()->residue(tresi).name1())) {
											 utility_exit_with_message("incompatible sequence for residue "+str(dresi));
										}
									}
									seq[dresi].resize(1);
									seq[dresi][1] = templates_cen.back()->residue(tresi).name1();
								}
								scgrp++;
								for(Size i = 1; i <= tres.size(); ++i) {
									for(Size j = 1; j < i; ++j) {
										tr << "add cst " << templates_fname.size() << " " << dres[i] << " " << tres[i] << " " << dres[j] << " " << tres[j] << endl;
										cst_sc.push_back(CST(templates_fname.size(),scgrp,dres[i],tres[i],dres[j],tres[j]));
									}
								}
							} else if("BACKBONE"==op2) {
								Sizes tres = parse_residues(getline(in),false);
								Sizes dres;
								for(Sizes::const_iterator i = tres.begin(); i != tres.end(); ++i) {
									if(template_map.back().count(*i)==0) utility_exit_with_message("BACKBONE: template "+tpdb+" residue "+str(*i)+" not mapped to design residue!");
									Size dri = template_map.back()[*i];
									dres.push_back(dri);
								}
								bbgrp++;
								for(Size i = 1; i <= tres.size(); ++i) {
									for(Size j = 1; j < i; ++j) {
										cst_bb.push_back(CST(templates_fname.size(),bbgrp,dres[i],tres[i],dres[j],tres[j]));
									}
								}
							} else {
								utility_exit_with_message("parse_config_file: don't understand TEMPLATE CMD: '"+op2+"'");
							}
						}
						if(!in.good()) break;
					}
				}
			} else {
					utility_exit_with_message("parse_config_file: don't understand CMD: '"+op+"'");
					break;
			}
		}
		for(Size i = 1; i <= cst_bb.size(); ++i) {
			Size r1=cst_bb[i].dres1;
			Size r2=cst_bb[i].dres2;
			if(r2<r1) {
				r2=cst_bb[i].dres1;
				r1=cst_bb[i].dres2;
			}
			if(r1 > nres ) continue;
			if(r2 > nhub*nres) nintercst++;
			else               nintracst++;
			ncst++;
		}
		for(Size i = 1; i <= cst_sc.size(); ++i) {
			Size r1=cst_sc[i].dres1;
			Size r2=cst_sc[i].dres2;
			if(r2<r1) {
				r2=cst_sc[i].dres1;
				r1=cst_sc[i].dres2;
			}
			if(r1 > nres ) continue;
			if(r2 > nhub*nres) nintercst++;
			else               nintracst++;
			ncst++;
		}
	}

	void apply_bb_csts(Pose & p, Size ssep=0, bool earlycst = false) {
		using namespace core::scoring::constraints;
		vector1<string> anames; anames.push_back("N"); anames.push_back("CA"); anames.push_back("C"); anames.push_back("O"); /*anames.push_back("CB");*/ anames.push_back("H");
		for(CSTs::iterator i = cst_bb.begin(); i != cst_bb.end(); ++i) {
			if(ssep && (hub_seq_sep(i->dres1,i->dres2) != ssep) && (!earlycst || numeric::random::uniform() > Real(nres)/6.0/Real(nintercst)) ) continue;
			if(i->active) continue; else i->active = true;
			if(hub_seq_sep(i->dres1,i->dres2) != ssep) tr << "early cst!!!" << endl;
			Pose & t(*templates_cen[i->tplt]);
			for(Size an1 = 1; an1 <= anames.size(); ++an1) {
				if(!p.residue(i->dres1).has(anames[an1]) || !t.residue(i->tres1).has(anames[an1])) continue;
				for(Size an2 = 1; an2 <= anames.size(); ++an2) {
					if(!p.residue(i->dres2).has(anames[an2]) || !t.residue(i->tres2).has(anames[an2])) continue;
					Real d = t.residue(i->tres1).xyz(anames[an1]).distance( t.residue(i->tres2).xyz(anames[an2]) );
					AtomID id1( p.residue(i->dres1).atom_index(anames[an1]), i->dres1 );
					AtomID id2( p.residue(i->dres2).atom_index(anames[an2]), i->dres2 );
					add_sym_cst( p, id1, id2, d, CSTSDMULT*sqrt(d) );
					// p.add_constraint( new AtomPairConstraint( id1, id2, new HarmonicFunc(d,CSTSDMULT*sqrt(d)) ) );
				}
			}
		}
	}
	void apply_dcsts(Pose & p, Size ssep=0) {
		using namespace core::scoring::constraints;
		for(vector1<DCST>::iterator i = dcst.begin(); i != dcst.end(); ++i) {
			if(ssep && (hub_seq_sep(i->rsd1,i->rsd2) != ssep) ) continue;
			if(i->active) continue; else i->active = true;
			if( !p.residue(i->rsd1).has(i->atm1) || !p.residue(i->rsd2).has(i->atm2) ) {
				utility_exit_with_message("CONSTRAINT missing atom "+i->atm1+" "+str(i->rsd1)+" "+i->atm2+" "+str(i->rsd2)+" "+str(i->d)+" "+str(i->sd) );
			}
			AtomID id1( p.residue(i->rsd1).atom_index(i->atm1), i->rsd1 );
			AtomID id2( p.residue(i->rsd2).atom_index(i->atm2), i->rsd2 );
			//tr << "DCST " << i->atm1 << " " << i->rsd1 << " " << i->atm2 << " " << i->rsd2 << " " << i->d << " " << i->sd << endl;
			add_sym_cst( p, id1, id2, i->d, i->sd );
		}
	}
	void apply_cen_csts(Pose & p, Size ssep=0, bool earlycst = false) {
		using namespace core::scoring::constraints;
		vector1<string> anames; anames.push_back("CB"); anames.push_back("CEN");
		if(p.is_fullatom()) utility_exit_with_message("apply_cen_csts called on fa pose!");
		for(CSTs::iterator i = cst_sc.begin(); i != cst_sc.end(); ++i) {
			if(ssep && (hub_seq_sep(i->dres1,i->dres2) != ssep) ) continue;
			if(i->active) continue; else i->active = true;
			Pose & t(*templates_cen[i->tplt]);
			for(Size an1 = 1; an1 <= anames.size(); ++an1) {
				if(!p.residue(i->dres1).has(anames[an1]) || !t.residue(i->tres1).has(anames[an1])) continue;
				for(Size an2 = 1; an2 <= anames.size(); ++an2) {
					if(!p.residue(i->dres2).has(anames[an2]) || !t.residue(i->tres2).has(anames[an2])) continue;
					Real d = t.residue(i->tres1).xyz(anames[an1]).distance( t.residue(i->tres2).xyz(anames[an2]) );
					AtomID id1( p.residue(i->dres1).atom_index(anames[an1]), i->dres1 );
					AtomID id2( p.residue(i->dres2).atom_index(anames[an2]), i->dres2 );
					add_sym_cst( p, id1, id2, d, CSTSDMULT*sqrt(d) );
				}
			}
			// Real d = t.residue(i->tres1).xyz("CEN").distance( t.residue(i->tres2).xyz("CEN") );
			// AtomID id1( p.residue(i->dres1).atom_index("CEN"), i->dres1 );
			// AtomID id2( p.residue(i->dres2).atom_index("CEN"), i->dres2 );
			// //tr << "constraint: " << ssep << " " << i->dres1 << ",CEN " << i->dres2 << ",CEN " << d << endl;
			// p.add_constraint( new AtomPairConstraint( id1, id2, new HarmonicFunc(d,CSTSDMULT*sqrt(d)) ) );
		}
		apply_bb_csts(p,ssep,earlycst);
		apply_dcsts(p,ssep);
	}
	void apply_fa_csts(Pose & p, Size ssep=0) {
		using namespace core::scoring::constraints;
		if(p.is_centroid()) utility_exit_with_message("apply_fa_csts called on cen pose!");
		for(CSTs::iterator i = cst_sc.begin(); i != cst_sc.end(); ++i) {
			if(ssep && (hub_seq_sep(i->dres1,i->dres2) != ssep) ) continue;
			if(i->active) continue; else i->active = true;
			Pose & t(*templates_fa[i->tplt]);
			for(Size a1 = 6; a1 <= p.residue(i->dres1).nheavyatoms(); ++a1) {
				string aname1 = p.residue(i->dres1).atom_name(a1);
				if(!t.residue(i->tres1).has(aname1)) {
					for(Size ia = 1; ia <= t.residue(i->tres1).nheavyatoms(); ++ia) tr << t.residue(i->tres1).atom_name(ia) << endl;
					utility_exit_with_message("can't find atom "+aname1+" in template residue "+str(i->tres1)+" t.aa "+t.residue(i->tres1).name()+" d.aa "+p.residue(i->dres1).name());
				}
				for(Size a2 = 6; a2 <= p.residue(i->dres2).nheavyatoms(); ++a2) {
					string aname2 = p.residue(i->dres2).atom_name(a2);
					if(!t.residue(i->tres2).has(aname2)) utility_exit_with_message("can't find atom "+aname2+" in template residue "+str(i->tres2));
					Real d = t.residue(i->tres1).xyz(aname1).distance( t.residue(i->tres2).xyz(aname2) );
					AtomID id1( p.residue(i->dres1).atom_index(aname1), i->dres1 );
					AtomID id2( p.residue(i->dres2).atom_index(aname2), i->dres2 );
					//tr << "constraint: " << ssep << " " << i->dres1 << "," << aname1 << " " << i->dres2 << "," << aname2 << " " << d << endl;
					add_sym_cst( p, id1, id2, d, CSTSDMULT/2.0*sqrt(d) );
					p.add_constraint( new AtomPairConstraint( id1, id2, new HarmonicFunc(d,CSTSDMULT/2.0*sqrt(d)) ) );
				}
			}
		}
		apply_bb_csts(p,ssep);
		apply_dcsts(p,ssep);
	}
	void apply_csts(Pose & p, Size ssep=0, bool earlycst = false) {
		if(  p.is_centroid() &&  p.is_fullatom() ) utility_exit_with_message("MADNESS!");
		if( !p.is_centroid() && !p.is_fullatom() ) utility_exit_with_message("MADNESS!");
		if(p.is_centroid()) apply_cen_csts(p,ssep,earlycst);
		if(p.is_fullatom()) apply_fa_csts(p,ssep);
	}
	std::string pick_sequence() const {
		string s = "Z";
		for(Size i = 2; i <= nres; ++i) {
			if(seq.count(i)) {
				Strings const & choices( seq.find(i)->second );
				string tmp = choices[ (Size)std::ceil( numeric::random::uniform()*((Real)choices.size()) ) ];
				//tr << i << " " << tmp << endl;
				if(tmp.size() != 1) utility_exit_with_message("sequence is bad: '"+tmp+"'");
				s += tmp;
			} else {
				if(ss[i]=='_') s += "G";
				else if(ss[i]=='H') s += "L";
				else if(ss[i]=='E') s += "V";
				else if(ss[i]=='L') s += "G";
				else utility_exit_with_message("bad ss "+ss[i]);
				//tr << i << " " << "LVG" << endl;
			}
		}
		if( s.size() != nres ) {
			utility_exit_with_message("particular sequence doesn't match length of ss");
		}
		tr << "pick sequence: " << s << endl;
		return s;
	}
	int get_highest_intrahub_seqsep() {
		int mx=1;
		for(CSTs::const_iterator i = cst_bb.begin(); i != cst_bb.end(); ++i) {
			if( i->dres1 <= nhub*nres && i->dres1 > mx) mx = i->dres1;
			if( i->dres2 <= nhub*nres && i->dres2 > mx) mx = i->dres2;
		}
		for(CSTs::const_iterator i = cst_sc.begin(); i != cst_sc.end(); ++i) {
			if( i->dres1 <= nhub*nres && i->dres1 > mx) mx = i->dres1;
			if( i->dres2 <= nhub*nres && i->dres2 > mx) mx = i->dres2;
		}
		return mx;
	}	
};


struct HubDenovo {
	Pose hub;
	ConstraintConfig cfg;
	ScoreFunctionOP sf3,sfsym,sfasym,sfsymnocst;
	protocols::moves::MoverOP rlxcst,rlxnocst,des,fragins3,fraginsL,fragins,cenmin;
	virtual ~HubDenovo() {}
	HubDenovo(std::string cstcfgfile) : cfg(cstcfgfile),sf3(NULL),sfsym(NULL),sfasym(NULL),sfsymnocst(NULL),rlxcst(NULL),rlxnocst(NULL),
		                                  des(NULL),fragins3(NULL),fraginsL(NULL),fragins(NULL),cenmin(NULL)
	{
		using namespace core::scoring;
		sf3 = new symmetry::SymmetricScoreFunction(ScoreFunctionFactory::create_score_function("score4_smooth"));
		Real cstwt = Real(cfg.nres)/Real(cfg.cst_bb.size());
		sf3->set_weight(atom_pair_constraint,cstwt);
		sf3->set_weight(vdw,2.0);
		sfsym  = getScoreFunction();
		sfsymnocst  = getScoreFunction();
		sfasym = new ScoreFunction(*sfsym);
		sfsym->set_weight(atom_pair_constraint,cstwt);
		core::chemical::ResidueTypeSetCAP rtsfa = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");

		std::map<string, vector1<core::fragment::FragDataOP> > fds( get_frags_map() );
		core::fragment::FragSetOP frags3 = make_frag_set(cfg.ss,fds);
		core::fragment::FragSetOP fragsl = make_frag_set(cfg.ssstr().substr(0,6),fds);

		fragins3  = new protocols::simple_moves::ClassicFragmentMover(frags3);
		fraginsL  = new protocols::simple_moves::ClassicFragmentMover(fragsl);
		{
		  	protocols::moves::RandomMoverOP tmp = new protocols::moves::RandomMover;
		  	tmp->add_mover(fragins3,0.5);
		  	tmp->add_mover(fraginsL,0.5);
		  	fragins = tmp;
		}
		des      = new protocols::flxbb::FlxbbDesign( sfsym, sfsym );
		rlxcst   = new protocols::relax::FastRelax (sfsym);
		rlxnocst = new protocols::relax::FastRelax(sfsymnocst);

		hub = *core::import_pose::pose_from_pdb(*rtsfa, option[OptionKeys::hub_pdb]() );

	}
	Pose make_start_pose() {
		Pose p(hub);
		core::pose::remove_upper_terminus_type_from_pose_residue(p,1);

		string tmpseq = cfg.pick_sequence();
		// tr << "start sequence " << seq << endl;
		for(Size ir = 2; ir <= cfg.nres; ++ir) {
			core::conformation::ResidueOP tmp = core::conformation::ResidueFactory::create_residue(
				*p.residue(1).residue_type_set().aa_map(core::chemical::aa_from_oneletter_code(tmpseq[ir-1]))[1] );
			tmp->seqpos(ir);
			tmp->chain(1);
			p.append_residue_by_bond( *tmp, true );
			p.set_omega(ir-1,180.0);
			p.set_omega(ir  ,180.0);
		}
		core::pose::add_upper_terminus_type_to_pose_residue(p,p.n_residue());
		//sfsym->show(p);

		make_symmetric_pose(p);
		FoldTree ft = p.fold_tree();
		//for(Size i = 1; i <= ft.num_jump(); ++i) tr << i << " " << ft.jump_edge(i) << endl;
		// tr << ft << endl;
		ft.set_jump_atoms( 1, "ORIG", "CA");
		ft.set_jump_atoms( 2, "ORIG", "CA");
		ft.set_jump_atoms( 3, "ORIG", "CA");
		p.fold_tree(ft);
		// tr << ft << endl;

		p.conformation().detect_bonds();

		core::util::switch_to_residue_type_set(p,"centroid");

		return p;

	}

	bool cen_fold(Pose & p) {
		
		// set up movemap and min mover
		using core::conformation::symmetry::SymDof;
		core::conformation::symmetry::SymmetryInfoCOP si = core::pose::symmetry::symmetry_info(p);
		std::map<Size,SymDof> const & dofs = si->get_dofs();
		core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	  movemap->set_jump(false);
	  movemap->set_bb(true);
	  movemap->set_chi(true);
		for(std::map<Size,SymDof>::const_iterator i = dofs.begin(); i != dofs.end(); ++i) {
			cout << "sym dof jump: " << i->first << endl;
	  	movemap->set_jump(i->first,true);
		}
		cenmin = new protocols::simple_moves::symmetry::SymMinMover( movemap, sf3, "dfpmin_armijo_nonmonotone", 1e-3, true, false, false );
		
		Size STOP = cfg.get_highest_intrahub_seqsep() + 4;
		tr << "rnd 1 ssep STOP " << STOP << endl;
		
		Real temp = 2.0;
		Pose last_cor_ori = p;
		for(Size icst = 1; icst <= STOP; ++icst) {
			cfg.apply_csts(p,icst,icst > 15);
			protocols::moves::MonteCarloOP mc = new protocols::moves::MonteCarlo( p, *sf3, temp );
			mc->set_autotemp( true, temp );	mc->set_temperature( temp );
			protocols::moves::RepeatMover( new protocols::moves::TrialMover( fragins, mc ), 4000/cfg.nres ).apply( p );
			mc->reset( p );
			if(icst%5==0) {
				tr << "cen min" << endl;
				cenmin->apply(p);
				sf3->show(p);
				//if( option[OptionKeys::hub_cen_energy_cut]()*4.0 < sf3->score(p) ) return false;
			}
			if( center_of_geom(p).z() < 0.0 ) p = last_cor_ori;
			else                              last_cor_ori = p;
			tr << icst << " " << sf3->score(p) << endl;
		}
		for(Size icst = STOP+1; icst <= cfg.nsub*cfg.nres; ++icst) {
			cfg.apply_csts(p,icst);
		}
		Real cstwt = sf3->get_weight(core::scoring::atom_pair_constraint);
		for(Size i = 1; i < 5; ++i) {
			sf3->set_weight(core::scoring::atom_pair_constraint,cstwt/Real(i*i));
			tr << i << " fin " << sf3->score(p) << endl;
			protocols::moves::MonteCarloOP mc = new protocols::moves::MonteCarlo( p, *sf3, temp );
			mc->set_autotemp( true, temp );	mc->set_temperature( temp );
			protocols::moves::RepeatMover( new protocols::moves::TrialMover( fragins, mc ), 400 ).apply( p );
			mc->reset( p );
			cenmin->apply(p);
			sf3->show(p);			
			//if( option[OptionKeys::hub_cen_energy_cut]()*4.0 < sf3->score(p) ) return false;
		}
		return true;
	}

	void min_as_poly_ala(Pose & p, ScoreFunctionOP sf) {
		for(Size i = 2; i <= cfg.nres; ++i) {
			if(p.residue(i).name3()=="GLY") continue;
			core::pose::replace_pose_residue_copying_existing_coordinates(p,i,p.residue(i).residue_type_set().name_map("ALA"));
		}
				// reapply csts
		p.remove_constraints();
		cfg.apply_csts(p);

		sf->show(p);
		sf->set_weight(core::scoring::atom_pair_constraint,1.0);
		sf->set_weight(core::scoring::    angle_constraint,1.0);
		core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	  	movemap->set_jump(false);
	  	movemap->set_bb(true);
	  	movemap->set_chi(true);
		protocols::simple_moves::symmetry::SymMinMover( movemap, sf, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false ).apply(p);
		sf->set_weight(core::scoring::atom_pair_constraint,1.0);
		sf->set_weight(core::scoring::    angle_constraint,1.0);
	}

	void run(Size NITER = 9999999999) {
		using namespace core::scoring;
		for(int iter = 1; iter < NITER; ++iter) {
			Pose tmp = make_start_pose();
			
			if(basic::options::option[basic::options::OptionKeys::hub_graphics]()) {
				protocols::viewer::add_conformation_viewer(tmp.conformation(),"test",1150,1150);
			}
			
			string fn = option[OptionKeys::out::file::o]() + "/" + utility::file_basename(cfg.fname) +"_"+ str(uniform()).substr(2,8) + ".pdb.gz";

			if(!cen_fold(tmp)) continue;

			{
				core::io::silent::SilentStructOP ss_out( new core::io::silent::ScoreFileSilentStruct );
			  	ss_out->fill_struct(tmp,fn);
			  	sfd.write_silent_struct( *ss_out, option[OptionKeys::out::file::o]() + "/" + option[ OptionKeys::out::file::silent ]()+"_cen.sc" );
			}

			if(tmp.energies().total_energy() > option[OptionKeys::hub_cen_energy_cut]() ) {
				tr << "CEN fail " << tmp.energies().total_energy() << endl;
				continue;
			}

			tr << "HIT " << tmp.energies().total_energy() << " " << fn << endl;

			if(!option[OptionKeys::hub_no_fa]()) {
				tmp.dump_pdb(fn+"_cen.pdb.gz");

				core::util::switch_to_residue_type_set(tmp,"fa_standard");
				// reapply csts
				tmp.remove_constraints();
				cfg.apply_csts(tmp);

				{
					sfsym->set_weight(core::scoring::atom_pair_constraint,10.0);
					sfsym->set_weight(core::scoring::    angle_constraint,10.0);
					core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
				  	movemap->set_jump(false); movemap->set_bb(true); movemap->set_chi(true);
					rlxcst->apply(tmp); // no cst
					// protocols::simple_moves::symmetry::SymMinMover( movemap, sfsym, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false ).apply(tmp);
					sfsym->set_weight(core::scoring::atom_pair_constraint,1.0);
					sfsym->set_weight(core::scoring::    angle_constraint,1.0);
					// min_as_poly_ala(tmp,sfsym,hocsts,cfg.nres);
					sfsym->show(tmp);
					//tmp.dump_pdb(fn+"_min.pdb.gz");
					tmp.constraint_set()->show_violations(tr,tmp,1,1.0);
				}
				// Real cstsc = tmp.energies().total_energies()[core::scoring::atom_pair_constraint];
				// tmp.remove_constraints();
				tmp.dump_pdb(fn+"_cst.pdb.gz");

				tr << "rlx nocst" << endl;
				rlxnocst->apply(tmp); // no cst

				sfsym->score(tmp); // rescore with cst
			}

			tmp.dump_pdb(fn);
			core::io::silent::SilentStructOP ss_out( new core::io::silent::ScoreFileSilentStruct );
			// ss_out->add_energy("cstsc",cstsc);
		  	ss_out->fill_struct(tmp,fn);

		  	sfd.write_silent_struct( *ss_out, option[OptionKeys::out::file::o]() + "/" + option[ OptionKeys::out::file::silent ]() );

		}

	}

};





void * run(void *) {
	string cstcfg = option[OptionKeys::hub_cst_cfg]();
	HubDenovo hd(cstcfg);
	hd.run();
}



int main(int argc, char *argv[]) {
	register_options();
	core::init(argc,argv);
	using namespace core::scoring;

	if(basic::options::option[basic::options::OptionKeys::hub_graphics]()) {
		protocols::viewer::viewer_main(&run);
	} else {
		run(NULL);
	}

	return 0;
}



