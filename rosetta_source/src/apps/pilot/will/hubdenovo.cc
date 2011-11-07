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
#include <protocols/relax/FastRelax.hh>
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
#include <apps/pilot/will/will_util.ihh>
#include <apps/pilot/will/mynamespaces.ihh>
#include <apps/pilot/will/frag_util.hh>

using core::conformation::symmetry::SymmData;
using core::conformation::symmetry::SymmDataOP;
using core::conformation::symmetry::SymmetryInfo;
using core::conformation::symmetry::SymmetryInfoOP;
using core::pose::symmetry::make_symmetric_pose;
using protocols::moves::symmetry::SymMinMover;
using protocols::moves::MoverOP;
using core::scoring::constraints::ConstraintOP;
using core::id::NamedAtomID;
typedef utility::vector1<core::Size> Sizes;
typedef utility::vector1<std::string> Strings;

OPT_KEY( String, hub_ss )
OPT_KEY( String, hub_sequence )
OPT_KEY( IntegerVector, hub_ho_cst )

void register_options() {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  NEW_OPT( hub_ss       ,"", "" );
  NEW_OPT( hub_sequence ,"", "" );
  NEW_OPT( hub_ho_cst   ,"", utility::vector1< Size >() );
}

static core::io::silent::SilentFileData sfd;

bool hocstcmp (std::pair<Size,Size> i, std::pair<Size,Size> j) { 
	return abs(i.first-i.second) < abs(j.first-j.second);
}

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
	}
	std::istringstream iss(buf);
	Strings v;
	string s;
	while(iss >> s) v.push_back(s);
	return v;
}


// fixed SS
struct HubDenovoConfig {
	Size nres;
	Size nsub;
	vector1<char> ss;
	std::map<string,Sizes > ssmap;
	std::map<Size,Strings > seq;
	Strings templates_fname;
	vector1<PoseOP> templates_fa,templates_cen;
	vector1<std::map<Size,Size> > template_map;
	vector1<char> sctemplate;
	vector1<std::pair<Size,vector1<std::pair<Size,Size> > > >   bbtemplate; //[ ( tplt, [ ( dres, tres ) ] ) ]
	vector1<std::pair<Size,vector1<std::pair<Size,Size> > > > disttemplate; //[ ( tplt, [ ( dres, tres ) ] ) ]
	core::chemical::ResidueTypeSetCAP crs,frs;
	HubDenovoConfig(string cfgfile) {
		nsub = 0;
		frs = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
		crs = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
		parse_config_file(cfgfile);
	}
	void parse_config_file(string infname) {
		utility::io::izstream in(infname);
		if(!in.good()) utility_exit_with_message("err opening "+infname);
		parse_config_file(in);
		utility_exit_with_message("debug parse_config_file");
	}
	Strings residue_names_from_sequence(string s) {
		core::chemical::ResidueTypeCAPs r = core::pose::residue_types_from_sequence(s,*frs,false);
		Strings rn;
		for(Size i = 1; i <= r.size(); ++i) {
			rn.push_back(r[i]->name());
		}
		return rn;
	}
	Sizes parse_residues(string s, bool aliases=true) {
		Sizes r;
		if(aliases && ssmap.count(s)>0) {
			r = ssmap[s];
		} else {
			bool isnum=true,isdash=true;
			for(Size i = 0; i < s.size(); ++i) {
				isnum &= (s[i] >= '0' && s[i] <= '9');
				isnum &= (s[i] >= '0' && s[i] <= '9' || s[i]=='-');
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
				utility_exit_with_message("parse_residues");	
			}
		}
		cout << "parse_residues " << s << " " << r << endl;
		return r;
	}
	Sizes parse_residues(Strings s, bool aliases=true) {
		cout << "parse_residues " << s << endl;
		Sizes res;
		for(Strings::const_iterator i = s.begin(); i != s.end(); ++i) {
			Sizes r = parse_residues(*i,aliases);
			res.insert(res.end(),r.begin(),r.end());
		}	
		return res;
	}
	vector1<Strings > parse_sequence(string s) {
		vector1<Strings > r;
		r.push_back(residue_names_from_sequence(s));
		return r;
	}
	void parse_config_file(std::istream & in) {
		string op,op2,val;
		while(true) {
			//cout << "ROOT peak: " << (char)in.peek() << endl;
			if(!(read_ignore_comments(in,op))) break;
			if("SUBUNITS"==op) {
				read_ignore_comments(in,nsub);
				//cout << "SUBUNITS: " << nsub << endl;
			} else
			if("SECSTRUCT"==op) {
				if(!nsub) utility_exit_with_message("SUBUNITS not set before SECSTRUCT");
				string ssstr;
				read_ignore_comments(in,ssstr);
				for(Size is = 0; is < ssstr.size(); ++is) {
					if(ssstr[is]!='H'&&ssstr[is]!='E'&&ssstr[is]!='L'&&ssstr[is]!='*') {
						utility_exit_with_message("parse_config_file: bad SECSTRUCT: "+ssstr.substr(is,1)+", "+ssstr);
					}
				}
				nres = ssstr.size();
				for(Size i = 1; i <= nsub; ++i) {
					for(Size j = 0; j < ssstr.size(); ++j) ss.push_back(ssstr[j]);
				}
				cout << "SECSTRUCT " << ss.size() << endl;
				{ // make ssmap
					std::map<char,Size> nss; nss['H']=0; nss['E']=0; nss['L']=0; nss['*']=0;
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
						// cout << "ssmap['" << i->first << "'] = ";
						char ssc = i->first[0];
						for(Sizes::const_iterator iv = i->second.begin(); iv != i->second.end(); ++iv) {
							// cout << " " << *iv;
							if(ss[*iv] != ssc) utility_exit_with_message("SECSTRUCT parse issue: "+i->first+" "+str(*iv));
						}
						// cout << endl;
					}
				}
			} else
			if("SEQUENCE"==op) {
				if(ssmap.size()==0) utility_exit_with_message("SECSTRUCT must come before SEQUENCE");
				while(true) {
					if("TEMPLATE" ==peek(in)) break;
					if("SECSTRUCT"==peek(in)) break;
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
					// cout << "sequence: " << op2 << " " << tmp << endl;
					// for(Size i = 1; i <= r.size(); ++i) {
					// 	cout << "   " << r[i];
					// 	for(Size j = 1; j <= s[i].size(); ++j) cout << " " << s[i][j];
					// 	cout << endl;
					// }
				}
				cout << "SEQUENCE" << endl;
			} else
			if("TEMPLATE"==op) {
				string tpdb,rsd;
				char buf[9999];
				read_ignore_comments(in,tpdb);
				if(tpdb.substr(tpdb.size()-4,4)!=".pdb"&&tpdb.substr(tpdb.size()-7,7)!=".pdb.gz") utility_exit_with_message("bad pdb: "+tpdb);
				templates_fname.push_back(tpdb);
				templates_fa.push_back( pose_from_pdb(tpdb) );
				PoseOP tmpcen = new Pose(*templates_fa.back());
				core::util::switch_to_residue_type_set(*tmpcen,"centroid");
				templates_cen.push_back(tmpcen);
				Size wcount = 0;
				while(true) {
					if(!(read_ignore_comments(in,op2))) break;
					while(true) {					
						++wcount;
						if(wcount > 100) utility_exit_with_message("TEMPLATE "+tpdb+" parse error on "+op2);
						if("MAPPING"==peek(in)||"SIDECHAINS"==peek(in)||"BACKBONE"==peek(in)||"DISTANCES"==peek(in)) break;
						if("MAPPING"==op2) {
							cout << "MAPPING" << endl;
							Sizes dres = parse_residues(getline(in));
							Sizes tres = parse_residues(getline(in));
							if(dres.size()!=tres.size()) utility_exit_with_message("MAPPING res no. mismatch");
							std::map<Size,Size> m;
							for(Size i = 1; i <= dres.size(); ++i) m[tres[i]] = dres[i];
							template_map.push_back(m);
							break;
						} else {
							if(template_map.size()!=templates_fa.size()) utility_exit_with_message("MAPPING must be defined first in TEMPLATE");
							Sizes tres = parse_residues(getline(in),false);
							if("SIDECHAIN"==op2) {
						
							} else if("BACKBONE"==op2) {

							} else if("DISTANCES"==op2) {

							} else {
							utility_exit_with_message("parse_config_file: don't understand TEMPLATE CMD: '"+op2+"'");
							}
						}
					}
				}
			} else {
					utility_exit_with_message("parse_config_file: don't understand CMD: '"+op+"'");
					break;
			}
		}
		for(Size i = 1; i <= templates_fa.size(); ++i) {
			templates_fa[i]->dump_pdb("TEMPLATE_FA__"+str(i)+".pdb");
			templates_fa[i]->dump_pdb("TEMPLATE_CEN_"+str(i)+".pdb");			
		}
	}

	void apply_cen_csts(Pose & p) const {
		if(p.is_fullatom()) utility_exit_with_message("apply_cen_csts called on fa pose!");
	}
	void apply_fa_csts(Pose & p) const {
		if(p.is_centroid()) utility_exit_with_message("apply_fa_csts called on cen pose!");		
	}

};


struct HubDenovo {
	string seq,ss;
	Size nres;
	Pose init;
	vector1<std::pair<Size,Size> > hocsts;
	ScoreFunctionOP sf3,sfsym,sfasym,sfsymnocst;
	MoverOP fragmove3;
	protocols::moves::MoverOP rlxcst,rlxnocst,des,fragins;
	HubDenovo(string _seq, string _ss, vector1<std::pair<Size,Size> > _hocsts): seq(_seq),ss(_ss),hocsts(_hocsts)
	{
		using namespace core::scoring;
		sf3 = new symmetry::SymmetricScoreFunction(ScoreFunctionFactory::create_score_function("score3"));
		sf3->set_weight(atom_pair_constraint,1.0);
		sf3->set_weight(    angle_constraint,1.0);
		sfsym  = getScoreFunction();
		sfsymnocst  = getScoreFunction();
		sfasym = new ScoreFunction(*sfsym);
		sfsym->set_weight(atom_pair_constraint,1.0);
		sfsym->set_weight(    angle_constraint,1.0);
		if( seq[0] != 'Z' ) utility_exit_with_message("first residue must be Z!!");
		core::chemical::ResidueTypeSetCAP rtsfa = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");

		Pose p( *core::import_pose::pose_from_pdb(*rtsfa,"input/CHC_HUB_FA.pdb") );
		core::pose::remove_upper_terminus_type_from_pose_residue(p,1);
		Size nres = ss.size();
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
		FoldTree ft = p.fold_tree();
		for(Size i = 1; i <= ft.num_jump(); ++i) cout << i << " " << ft.jump_edge(i) << endl;
		// cout << ft << endl;
		ft.set_jump_atoms( 1, "ORIG", "CA");
		ft.set_jump_atoms( 2, "ORIG", "CA");
		ft.set_jump_atoms( 3, "ORIG", "CA");		
		p.fold_tree(ft);
		// cout << ft << endl;

		p.conformation().detect_bonds();

		std::map<string, vector1<core::fragment::FragDataOP> > fds( get_frags_map() );
		core::fragment::FragSetOP frags3 = make_frag_set(ss,fds);

		fragins  = new protocols::basic_moves::ClassicFragmentMover(frags3);
		des      = new protocols::flxbb::FlxbbDesign( sfsym, sfsym );
		rlxcst   = new protocols::relax::FastRelax (sfsym);
		rlxnocst = new protocols::relax::FastRelax(sfsymnocst);
		core::util::switch_to_residue_type_set(p,"centroid");

		init = p;

	}

	void add_ho_csts(Pose & p, Size ir, Size jr) {
		// Size ir = hocsts[i].first;
		// Size jr = hocsts[i].second;
		using namespace core::scoring::constraints;
		cout << "bond (name H and resi " << ir << "), (name O and resi " << jr << ")" << endl;
		p.add_constraint( new AtomPairConstraint( AtomID(p.residue(ir).atom_index("H"),ir),
			                                      AtomID(p.residue(jr).atom_index("O"),jr),	
												  new HarmonicFunc(2.0,0.5)) );
		p.add_constraint( new AngleConstraint( AtomID(p.residue(ir).atom_index("N"),ir),	
				                               AtomID(p.residue(ir).atom_index("H"),ir),
			                                   AtomID(p.residue(jr).atom_index("O"),jr),	
										       new CircularHarmonicFunc(3.14159,1.0)) );
		p.add_constraint( new AngleConstraint( AtomID(p.residue(ir).atom_index("H"),ir),
				                               AtomID(p.residue(jr).atom_index("O"),jr),
			                                   AtomID(p.residue(jr).atom_index("C"),jr),	
										       new CircularHarmonicFunc(3.14159,1.0)) );
		if( ir > nres || jr > nres ) {
			cout << "sym cst!" << endl;
			if(ir > nres) ir -=   nres;
			else          ir += 2*nres;
			if(jr > nres) jr -=   nres;
			else          jr += 2*nres;
			cout << "bond (name H and resi " << ir << "), (name O and resi " << jr << ")" << endl;			
			p.add_constraint( new AtomPairConstraint( AtomID(p.residue(ir).atom_index("H"),ir),
				                                      AtomID(p.residue(jr).atom_index("O"),jr),	
													  new HarmonicFunc(2.0,0.5)) );
			p.add_constraint( new AngleConstraint( AtomID(p.residue(ir).atom_index("N"),ir),
					                               AtomID(p.residue(ir).atom_index("H"),ir),
				                                   AtomID(p.residue(jr).atom_index("O"),jr),	
											       new CircularHarmonicFunc(3.14159,1.0)) );
			p.add_constraint( new AngleConstraint( AtomID(p.residue(ir).atom_index("H"),ir),	
					                               AtomID(p.residue(ir).atom_index("O"),ir),
				                                   AtomID(p.residue(jr).atom_index("C"),jr),	
											       new CircularHarmonicFunc(3.14159,1.0)) );
		}

	}

	void cen_fold(Pose & p) {

		for(vector1<std::pair<Size,Size> >::const_iterator i = hocsts.begin(); i != hocsts.end(); ++i) {
			Size ir = i->first; Size jr = i->second;
			add_ho_csts(p,ir,jr);
			Real temp = 2.0;
			protocols::moves::MonteCarloOP mc = new protocols::moves::MonteCarlo( p, *sf3, temp );
			mc->set_autotemp( true, temp );	mc->set_temperature( temp );
			protocols::moves::RepeatMover( new protocols::moves::TrialMover( fragmove3, mc ), 10000/hocsts.size() ).apply( p );
			mc->reset( p );
		}
		Real temp = 2.0;
		protocols::moves::MonteCarloOP mc = new protocols::moves::MonteCarlo( p, *sf3, temp );
		mc->set_autotemp( true, temp );	mc->set_temperature( temp );
		protocols::moves::RepeatMover( new protocols::moves::TrialMover( fragmove3, mc ), 10000 ).apply( p );
		mc->reset( p );
		//sf3->show(p);
	}
	
	void min_as_poly_ala(Pose & p, ScoreFunctionOP sf) {
		for(Size i = 2; i <= nres; ++i) {
			if(p.residue(i).name3()=="GLY") continue;
			core::pose::replace_pose_residue_copying_existing_coordinates(p,i,p.residue(i).residue_type_set().name_map("ALA"));
		}
				// reapply csts
			p.remove_constraints();
			for(vector1<std::pair<Size,Size> >::const_iterator i = hocsts.begin(); i != hocsts.end(); ++i) 
				add_ho_csts(p,i->first,i->second);

		sf->show(p);
		sf->set_weight(core::scoring::atom_pair_constraint,1.0);
		sf->set_weight(core::scoring::    angle_constraint,1.0);
		core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	  	movemap->set_jump(false);
	  	movemap->set_bb(true);
	  	movemap->set_chi(true);
		protocols::moves::symmetry::SymMinMover( movemap, sf, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false ).apply(p);
		sf->set_weight(core::scoring::atom_pair_constraint,1.0);
		sf->set_weight(core::scoring::    angle_constraint,1.0);
	}

	void run(Size NITER = 9999999999) {
		using namespace core::scoring;
		for(int iter = 1; iter < NITER; ++iter) {
			Pose tmp(init);
			
			cen_fold(tmp);

			if(tmp.energies().total_energies()[atom_pair_constraint] > 60.0) {
				cout << "cst fail " << tmp.energies().total_energies()[atom_pair_constraint] << endl;
				continue;
			}

			string fn = option[OptionKeys::out::file::o]() + "/" + ss +"_"+ str(uniform()).substr(2,8) + ".pdb.gz";
			cout << "HIT " << fn << endl;

			tmp.dump_pdb(fn+"_cen.pdb.gz");
			utility_exit_with_message("foo");
			core::util::switch_to_residue_type_set(tmp,"fa_standard");
			// reapply csts
			tmp.remove_constraints();
			for(vector1<std::pair<Size,Size> >::const_iterator i = hocsts.begin(); i != hocsts.end(); ++i) {
				add_ho_csts(tmp,i->first,i->second);
			}

			{
				sfsym->set_weight(core::scoring::atom_pair_constraint,10.0);
				sfsym->set_weight(core::scoring::    angle_constraint,10.0);
				core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
			  	movemap->set_jump(false); movemap->set_bb(true); movemap->set_chi(true);
				rlxcst->apply(tmp); // no cst
				// protocols::moves::symmetry::SymMinMover( movemap, sfsym, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false ).apply(tmp);
				sfsym->set_weight(core::scoring::atom_pair_constraint,1.0);
				sfsym->set_weight(core::scoring::    angle_constraint,1.0);
				// min_as_poly_ala(tmp,sfsym,hocsts,nres);
				sfsym->show(tmp);
				//tmp.dump_pdb(fn+"_min.pdb.gz");
				tmp.constraint_set()->show_violations(cout,tmp,1,1.0);
			}
			tmp.remove_constraints();

			// cout << "flxbb" << endl;
			// des.apply(tmp); // with cst

			cout << "rlx nocst" << endl;
			rlxnocst->apply(tmp); // no cst
			sfsym->show(tmp); // rescore with cst

			tmp.dump_pdb(fn);
			core::io::silent::SilentStructOP ss_out( new core::io::silent::ScoreFileSilentStruct );
		  	ss_out->fill_struct(tmp,fn);

		  	sfd.write_silent_struct( *ss_out, option[OptionKeys::out::file::o]() + "/" + option[ OptionKeys::out::file::silent ]() );

		}
		
	}

};









int main(int argc, char *argv[]) {
	register_options();
	core::init(argc,argv);
	using namespace core::scoring;

	string seq = option[OptionKeys::hub_sequence]();
	string ss = option[OptionKeys::hub_ss]();
	vector1<std::pair<Size,Size> > hocsts; {
		Sizes tmp = option[OptionKeys::hub_ho_cst]();
		for(Size i = 1; i < tmp.size(); i += 2) hocsts.push_back(std::pair<Size,Size>(tmp[i],tmp[i+1]));
		std::sort(hocsts.begin(),hocsts.end(),hocstcmp);
	}
	for(Size i = 1; i <= hocsts.size(); ++i) cout << hocsts[i].first << " " << hocsts[i].second << endl;

	HubDenovoConfig hdc("input/run/cylindrin_hairpin.hubdenovo");

	HubDenovo hd(seq,ss,hocsts);

	hd.run();

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
