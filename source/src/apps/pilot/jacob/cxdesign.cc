#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/smhybrid.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/matdes.OptionKeys.gen.hh>
//#include <basic/options/keys/willmatch.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
#include <core/pack/optimizeH.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/packing/compute_holes_score.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <numeric/conversions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/scoring/ImplicitFastClashCheck.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
// #include <devel/init.hh>

// #include <core/scoring/constraints/LocalCoordinateConstraint.hh>
#include <apps/pilot/will/will_util.ihh>
#include <apps/pilot/will/mynamespaces.ihh>

#include <sys/stat.h>

#include <utility/excn/Exceptions.hh>


OPT_1GRP_KEY( Integer      , cxdock, sphere       ) // 12872 32672 78032 8192
OPT_1GRP_KEY( Real         , cxdock, clash_dis    )
OPT_1GRP_KEY( Real         , cxdock, contact_dis  )
//OPT_1GRP_KEY( Real         , cxdock, num_contacts )
OPT_1GRP_KEY( Boolean      , cxdock, dump  )

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	OPT( in::file::s );
	//NEW_OPT( cxdock::syms  ,"CX symmitries", utility::vector1< Size >() );
	NEW_OPT( cxdock::clash_dis   ,"min acceptable contact dis", 3.6 );
	NEW_OPT( cxdock::contact_dis ,"max acceptable contact dis", 6.0 );
	//NEW_OPT( cxdock::num_contacts ,"required no. contacts", 20.0 );
	NEW_OPT( cxdock::sphere       ,"sph points", 8192 );
	NEW_OPT( cxdock::dump         ,"", false );
}


typedef numeric::xyzVector<Real> Vecf;
typedef numeric::xyzMatrix<Real> Matf;


using core::kinematics::Stub;
using protocols::scoring::ImplicitFastClashCheck;
using core::pose::Pose;
using core::conformation::ResidueOP;

static basic::Tracer TR( "pentcb" );
static core::io::silent::SilentFileData sfd;


#include <apps/pilot/will/sicfast.ihh>


struct Hit {
	int iss,irt,cbc,sym;
	core::kinematics::Stub s1,s2;
	Hit(int is, int ir, int cb, int sm) : iss(is),irt(ir),cbc(cb),sym(sm) {}
};
bool cmp(Hit i,Hit j) { return i.cbc > j.cbc; }


// Pose needs to be scored before this will work.
void new_sc(core::pose::Pose &pose, int nmono, Real& int_area, Real& sc) {
	using namespace core;
	core::scoring::sc::ShapeComplementarityCalculator scc;
	scc.Init();
	for ( Size i=1; i<=nmono; ++i ) scc.AddResidue(0, pose.residue(i));
	for ( Size i=1; i<=nmono; ++i ) scc.AddResidue(1, pose.residue(i+nmono));
	if ( scc.Calc() ) {
		sc = scc.GetResults().sc;
		int_area = scc.GetResults().surface[2].trimmedArea;
	}
}


Real
average_degree (Pose const &pose, Size nmono, vector1<Size> mutalyze_pos, Real distance_threshold=8.0)
{

	Size count_neighbors( 0 );

	for ( Size i = 1; i <= mutalyze_pos.size(); ++i ) {
		Size ires = mutalyze_pos[i];
		core::conformation::Residue const resi( pose.conformation().residue( ires ) );
		Size resi_neighbors( 0 );
		for ( Size jres = 1; jres <= nmono; ++jres ) {
			core::conformation::Residue const resj( pose.residue( jres ) );
			Real const distance( resi.xyz( resi.nbr_atom() ).distance( resj.xyz( resj.nbr_atom() ) ) );
			if ( distance <= distance_threshold ) {
				++count_neighbors;
				++resi_neighbors;
			}
		}
		//TR << "avg_deg of " << resi.name3() << ires << " = " << resi_neighbors << std::endl;
	}

	return( (Real) count_neighbors / mutalyze_pos.size() );

}


vector1<Size> design(Pose & p, Size nmono, Size nsym) {
	using namespace core::pack::task;
	ScoreFunctionOP sf = core::scoring::get_score_function();
	PackerTaskOP task = TaskFactory::create_packer_task(p);
	// task->initialize_extra_rotamer_flags_from_command_line();
	vector1< bool > aac(20,false);
	aac[core::chemical::aa_ala] = true;
	//aac[core::chemical::aa_cys] = true;
	//aac[core::chemical::aa_asp] = true;
	//aac[core::chemical::aa_glu] = true;
	aac[core::chemical::aa_phe] = true;
	//aac[core::chemical::aa_gly] = true;
	//aac[core::chemical::aa_his] = true;
	aac[core::chemical::aa_ile] = true;
	//aac[core::chemical::aa_lys] = true;
	aac[core::chemical::aa_leu] = true;
	aac[core::chemical::aa_met] = true;
	//aac[core::chemical::aa_asn] = true;
	//aac[core::chemical::aa_pro] = true;
	//aac[core::chemical::aa_gln] = true;
	//aac[core::chemical::aa_arg] = true;
	//aac[core::chemical::aa_ser] = true;
	//aac[core::chemical::aa_thr] = true;
	aac[core::chemical::aa_val] = true;
	aac[core::chemical::aa_trp] = true;
	aac[core::chemical::aa_tyr] = true;

	vector1<Size> iface(nmono,0);
	for ( Size ir = 1; ir <= nmono; ++ir ) {
		if ( p.residue(ir).name3()=="CYS" || p.residue(ir).name3()=="GLY" || p.residue(ir).name3()=="PRO" ) {
			iface[ir] = 0;
			continue;
		}
		if ( !p.residue(ir).is_protein()  ) continue;
		Real closestcb = 9e9;
		for ( Size jr = nmono+1; jr <= nsym*nmono; ++jr ) {
			if ( p.residue(jr).name3()=="GLY" ) continue;
			if ( !p.residue(jr).is_protein()  ) continue;
			Real d = p.xyz(AtomID(5,ir)).distance( p.xyz(AtomID(5,jr)) );
			closestcb = min(closestcb,d);
		}
		if ( closestcb < 9.0 ) {
			iface[ir] = 2;
		} else if ( closestcb < 12.0 ) {
			iface[ir] = 1;
		}
	}

	for ( Size i = 1; i <= nmono; ++i ) {
		if       ( iface[i] == 2 ) {
			bool tmp = aac[p.residue(i).aa()];
			aac[p.residue(i).aa()] = true;
			task->nonconst_residue_task(i).restrict_absent_canonical_aas(aac);
			task->nonconst_residue_task(i).initialize_extra_rotamer_flags_from_command_line();
			aac[p.residue(i).aa()] = tmp;
		} else if ( iface[i] == 1 ) {
			task->nonconst_residue_task(i).restrict_to_repacking();
			task->nonconst_residue_task(i).initialize_extra_rotamer_flags_from_command_line();
		} else if ( iface[i] == 0 ) {
			task->nonconst_residue_task(i).prevent_repacking();
		}

	}

	//Real rorig = sf->get_weight(core::scoring::fa_rep);
	//sf->set_weight(core::scoring::fa_rep,rorig/4.0);
	Real worig = sf->get_weight(core::scoring::res_type_constraint);
	if ( worig == 0.0 ) sf->set_weight(core::scoring::res_type_constraint,1.0);
	utility::vector1< core::scoring::constraints::ConstraintCOP > res_cst = add_favor_native_cst(p);
	p.add_constraints( res_cst );

	core::pack::make_symmetric_PackerTask_by_truncation(p,task);
	protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
	repack.apply(p);

	// cleanup 2
	p.remove_constraints( res_cst );
	sf->set_weight(core::scoring::res_type_constraint,worig);

	int acnt = 0;
	for ( int i = 1; i <= nmono; ++i ) {
		if ( iface[i]==2 && p.residue(i).name3()=="ALA" ) acnt++;
	}

	vector1<Size> ifa;
	for ( Size i = 1; i <= iface.size(); ++i ) if ( iface[i]==2 ) ifa.push_back(i);
	return ifa;

}


void repack(Pose & p, Size nmono, vector1<Size> const & iface) {
	using namespace core::pack::task;
	ScoreFunctionOP sf = core::scoring::get_score_function();
	PackerTaskOP task = TaskFactory::create_packer_task(p);
	for ( Size i = 1; i <= nmono; ++i ) {
		if       ( iface[i] == 2 ) {
			task->nonconst_residue_task(i).restrict_to_repacking();
			task->nonconst_residue_task(i).initialize_extra_rotamer_flags_from_command_line();
		} else if ( iface[i] == 1 ) {
			task->nonconst_residue_task(i).restrict_to_repacking();
			task->nonconst_residue_task(i).initialize_extra_rotamer_flags_from_command_line();
		} else if ( iface[i] == 0 ) {
			task->nonconst_residue_task(i).prevent_repacking();
		}
	}
	core::pack::make_symmetric_PackerTask_by_truncation(p,task);
	protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
	repack.apply(p);
}


void
minimize(Pose & pose,  utility::vector1<Size> const & iface) {
	ScoreFunctionOP sf = core::scoring::get_score_function();

	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	movemap->set_jump(true);
	movemap->set_chi(false);
	movemap->set_bb(false);
	for ( utility::vector1<Size>::const_iterator i = iface.begin(); i != iface.end(); i++ ) {
		movemap->set_bb (*i,true);
		movemap->set_chi(*i,true);
	}
	core::pose::symmetry::make_symmetric_movemap( pose, *movemap );
	protocols::simple_moves::symmetry::SymMinMover m( movemap, sf, "lbfgs_armijo_nonmonotone", 1e-5, true, false, false );
	m.apply(pose);
}


Real get_ddg(Pose p, Size nmono, vector1<Size> const & iface) {
	ScoreFunctionOP sf = core::scoring::get_score_function();
	repack(p,nmono,iface);
	minimize(p,iface);
	Real s1 = sf->score(p);
	core::conformation::symmetry::SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(p);
	std::map<Size,core::conformation::symmetry::SymDof> dofs = sym_info->get_dofs();
	int sym_jump = 0;
	for ( std::map<Size,core::conformation::symmetry::SymDof>::iterator i = dofs.begin(); i != dofs.end(); i++ ) {
		Size jump_num = i->first;
		if ( sym_jump == 0 ) sym_jump = jump_num;
		else utility_exit_with_message("Can only handle one subunit!");
	}
	if ( sym_jump == 0 ) utility_exit_with_message("No jump defined!");

	core::kinematics::Jump j = p.jump(sym_jump);
	j.set_translation( j.get_translation()*500.0 );
	p.set_jump(sym_jump,j);

	repack(p,nmono,iface);
	minimize(p,iface);
	Real s0 = sf->score(p);

	return s1-s0;
}


void cxdock_design(Pose const init, std::string const & fn, vector1<xyzVector<double> > const & ssamp, int iss, int irt, int ic, Real cbc) {
	using namespace basic::options;

	vector1<double> asamp; for ( Real i = 0; i < 180; ++i ) asamp.push_back(i);

	vector1<Vecf> bb0tmp,cb0tmp;
	for ( int ir = 1; ir <= init.size(); ++ir ) {
		if ( !init.residue(ir).is_protein() ) continue;
		for ( int ia = 1; ia <= ((init.residue(ir).has("CB"))?5:4); ++ia ) {
			bb0tmp.push_back(init.xyz(AtomID(ia,ir)));
		}
		if ( init.secstruct(ir)=='H' ) {
			if ( init.residue(ir).has("CB") ) cb0tmp.push_back(init.xyz(AtomID(5,ir)));
			else                           cb0tmp.push_back(init.xyz(AtomID(4,ir)));
		}
	}
	vector1<Vecf> const bb0(bb0tmp);
	vector1<Vecf> const cb0(cb0tmp);

	Matf Rsym = rotation_matrix_degrees(Vec(1,0,0),360.0/(Real)ic);

	if ( iss > ssamp.size() ) utility_exit_with_message("bad iss");
	if ( irt > asamp.size() ) utility_exit_with_message("bad irt");
	Vec axs = ssamp[iss];
	Real const rot = asamp[irt];
	Matf const R = rotation_matrix_degrees(axs,rot);

	vector1<Vecf> bb1 = bb0;
	vector1<Vecf> cb1 = cb0;
	for ( vector1<Vecf>::iterator i = bb1.begin(); i != bb1.end(); ++i ) *i = R*(*i);
	for ( vector1<Vecf>::iterator i = cb1.begin(); i != cb1.end(); ++i ) *i = R*(*i);
	vector1<Vecf> bb2 = bb1;
	vector1<Vecf> cb2 = cb1;
	for ( vector1<Vecf>::iterator i = bb2.begin(); i != bb2.end(); ++i ) *i = Rsym*(*i);
	for ( vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i ) *i = Rsym*(*i);
	Real cbctmp;
	Real t = sicfast(bb1,bb2,cb1,cb2,cbctmp);

	Pose p(init);
	rot_pose(p,R);
	//Pose q(p);
	//rot_pose(q,Rsym);
	trans_pose(p,Vec(0,0,t));
	Vec cen = Vec(0,t/2.0/tan( numeric::constants::d::pi / (Real)ic ),t/2.0);
	trans_pose(p,-cen);
	//trans_pose(q,-cen);

	core::pose::symmetry::make_symmetric_pose(p);

	if ( option[OptionKeys::cxdock::dump]() ) {
		p.dump_pdb(option[OptionKeys::out::file::o]()+"/"+utility::file_basename(fn)+"_"+str(ic)+"_"+str(iss)+"_"+str(irt)+".pdb.gz");
	}
	if ( fabs(cbc-cbctmp) > 0.1 ) utility_exit_with_message("cb counts disagree!!! "+str(cbctmp)+" should be "+str(cbc));
	if ( option[OptionKeys::cxdock::dump]() ) {
		return;
	}

	for ( int ir = 1; ir <= 1*init.size(); ++ir ) {
		if ( !p.residue(ir).is_protein() ) continue;
		for ( int ia = 1; ia <= ((p.residue(ir).has("CB"))?5:4); ++ia ) {
			Vec const pi = p.xyz(AtomID(ia,ir));
			for ( int jr = init.size()+1; jr <= 2*init.size(); ++jr ) {
				if ( !p.residue(jr).is_protein() ) continue;
				for ( int ja = 1; ja <= ((p.residue(jr).has("CB"))?5:4); ++ja ) {
					if ( pi.distance_squared(p.xyz(AtomID(ja,jr))) < 9.0 ) {
						return;
					}
				}
			}
		}
	}

	vector1<Size> iface = design(p,init.size(),ic);

	// Calculate the surface area and surface complementarity for the interface
	Real int_area = 0; Real sc = 0;
	new_sc(p,init.size(),int_area,sc);
	Real avg_deg = average_degree(p,init.size(),iface);

	Real avg_nchi = 0.0;
	for ( Size i = 1; i <= iface.size(); ++i ) avg_nchi += p.residue(iface[i]).nchi();
	avg_nchi /= iface.size();

	core::scoring::get_score_function()->score(p);

	string tag = utility::file_basename(fn)+"_"+str(ic)+"_"+str(iss)+"_"+str(irt);

	// Create a scorefile struct, add custom metrics to it
	core::io::silent::SilentStructOP ss_out( new core::io::silent::ScoreFileSilentStruct );
	ss_out->fill_struct(p,tag+".pdb.gz");
	ss_out->add_energy("avg_deg",avg_deg);
	ss_out->add_energy("sc_int_area",int_area);
	ss_out->add_energy("sc",sc);
	ss_out->add_energy("avg_nchi",avg_nchi);

	Real ddg = get_ddg(p,init.size(),iface);

	ss_out->add_energy("ddg",ddg);

	// // Write the scorefile
	sfd.write_silent_struct( *ss_out, option[OptionKeys::out::file::o]() + "/" + option[ OptionKeys::out::file::silent ]() );

	//p.dump_pdb(option[OptionKeys::out::file::o]() + "/" + tag+".pdb.gz");
	minimize(p,iface);
	p.dump_pdb(option[OptionKeys::out::file::o]() + "/" + tag+".pdb.gz");

}


int main(int argc, char *argv[]) {
	try {

		using namespace basic::options;

		register_options();
		devel::init(argc,argv);

		Size const NSS = basic::options::option[basic::options::OptionKeys::cxdock::sphere]();
		vector1<xyzVector<double> > ssamp(NSS); {
			izstream is;
			basic::database::open(is,"geometry/sphere_"+str(NSS)+".dat");
			for ( int i = 1; i <= NSS; ++i ) {
				double x,y,z;
				is >> x >> y >> z;
				ssamp[i] = xyzVector<double>(x,y,z);
			}
			is.close();
		}

		utility::io::izstream res(option[OptionKeys::in::file::l]()[1]);
		string dummy,fn;
		int iss,ic,nss,irt;
		Real cbc;
		while ( res >> dummy >> fn >> ic >> iss >> nss >> irt >> cbc ) {
			TR <<" "<< fn <<" "<< iss <<" "<< ic <<" "<< nss <<" "<< irt <<" "<< cbc << endl;
			option[OptionKeys::symmetry::symmetry_definition]("input/sym/C"+str(ic)+".sym");
			if ( NSS != nss ) utility_exit_with_message("wrong ssamp!!!");
			Pose pnat;
			pose_from_file(pnat,fn, core::import_pose::PDB_file);
			trans_pose(pnat,-center_of_geom(pnat,1,pnat.size()));
			core::scoring::dssp::Dssp dssp(pnat);
			dssp.insert_ss_into_pose(pnat);
			Size nhelix=0;
			for ( Size ir = 2; ir <= pnat.size()-1; ++ir ) {
				if ( pnat.residue(ir).is_lower_terminus() ) remove_lower_terminus_type_from_pose_residue(pnat,ir);//goto cont1;
				if ( pnat.residue(ir).is_upper_terminus() ) remove_upper_terminus_type_from_pose_residue(pnat,ir);//goto cont1;
			} goto done1; cont1: TR << "skipping " << fn << std::endl; continue; done1:

			cxdock_design(pnat,fn,ssamp,iss,irt,ic,cbc);
		}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}

