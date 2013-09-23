// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

#include <basic/basic.hh>
#include <basic/database/open.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/keys/rblinker.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/simple_moves/BBGaussianMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/scoring/methods/ImplicitClashEnergy.hh>
#include <protocols/viewer/viewers.hh>
#include <sstream>
#include <set>
//#include <utility/io/izstream.hh>
//#include <utility/io/ozstream.hh>

#include <apps/pilot/will/mynamespaces.ihh>

using core::kinematics::Stub;
using core::conformation::Residue;
using core::Real;

static basic::Tracer TR("rblinker2");

// static numeric::random::RandomGenerator RG(8334046);
// 
// inline void xform_pose( core::pose::Pose & pose, Stub const & s ) {
// 	for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
// 		for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
// 			core::id::AtomID const aid(core::id::AtomID(ia,ir));
// 			pose.set_xyz( aid, s.local2global(pose.xyz(aid)) );
// 		}
// 	}
// }
// 
// 
// 
// // petf   37   42   45   75
// // hyd   101  156  334  338
// // psI  1492 1536 1489 1495
// inline Real hyd_petf_sf4_dis(Pose const & hyd, Stub const & shyd, Pose const & petf, Stub const & spetf ) {
// 	Vec sf4A = ( hyd.xyz(AtomID(5, 101))+ hyd.xyz(AtomID(5, 156))+ hyd.xyz(AtomID(5, 334))+ hyd.xyz(AtomID(5, 338)))/4.0;
// 	Vec sf4B = (petf.xyz(AtomID(5,  37))+petf.xyz(AtomID(5,  42))+petf.xyz(AtomID(5,  45))+petf.xyz(AtomID(5,  75)))/4.0;
// 	return shyd.local2global(sf4A).distance(spetf.local2global(sf4B));
// }
// inline Real psI_petf_sf4_dis(Pose const & psI, Pose const & petf, Stub const & spetf ) {
// 	Vec sf4A = ( psI.xyz(AtomID(5,1492))+ psI.xyz(AtomID(5,1536))+ psI.xyz(AtomID(5,1489))+ psI.xyz(AtomID(5,1495)))/4.0;
// 	Vec sf4B = (petf.xyz(AtomID(5,  37))+petf.xyz(AtomID(5,  42))+petf.xyz(AtomID(5,  45))+petf.xyz(AtomID(5,  75)))/4.0;
// 	return sf4A.distance(spetf.local2global(sf4B));
// }
// 
// // get stup that aligns r1 to r2
// Stub getxform(Residue const & r1, Residue const & r2) {
// 	Stub s;
// 	s.M = alignVectorSets(r1.xyz(1)-r1.xyz(2),r1.xyz(3)-r1.xyz(2),r2.xyz(1)-r2.xyz(2),r2.xyz(3)-r2.xyz(2));
// 	s.v = r2.xyz(2)-s.M*r1.xyz(2);
// 	return s;
// }



inline void xform_pose( core::pose::Pose & pose, core::kinematics::Stub const & s, Size sres=1, Size eres=0 ) {
	if(eres==0) eres = pose.n_residue();
	for(Size ir = sres; ir <= eres; ++ir) {
		for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
			core::id::AtomID const aid(core::id::AtomID(ia,ir));
			pose.set_xyz( aid, s.local2global(pose.xyz(aid)) );
		}
	}
}
inline void xform_pose_rev( core::pose::Pose & pose, core::kinematics::Stub const & s ) {
	for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
		for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
			core::id::AtomID const aid(core::id::AtomID(ia,ir));
			pose.set_xyz( aid, s.global2local(pose.xyz(aid)) );
		}
	}
}




class SimpleBBMover : public protocols::moves::Mover {
	Size start_,stop_;
	Real mag_;
public:
	SimpleBBMover(Size start, Size stop, Real mag) : start_(start),stop_(stop),mag_(mag) {}
	Real magnitude(        ) { return mag_; }
	void magnitude(Real mag) { mag_ = mag; }	
	void apply(core::pose::Pose & pose) {
		Size i = start_-1 + (Size)std::ceil(uniform()*(stop_-start_+1));
		if(uniform()<0.5) pose.set_phi(i,pose.phi(i)+gaussian()*mag_);
		else              pose.set_psi(i,pose.psi(i)+gaussian()*mag_);
	}
	std::string get_name() const { return "SimpleBBMover"; }
};
typedef utility::pointer::owning_ptr<SimpleBBMover> SimpleBBMoverOP;

std::string printbits(unsigned long nn) {
	std::string s = "";
	for(Size i = 0; i < 8*sizeof(nn); ++i) {
		if(i%8==0) s += " ";
		if(nn & 1UL << (8*sizeof(nn)-1-i)) s += "1";
		else                               s += "0";
	}
	return s;
}

std::string bin2string(unsigned long bin, Size nres) {
	std::string s = "";
	for(Size i = 0; i < nres; ++i) {
		for(Size j = 0; j < 2; ++j) {
			int tmp = ( bin >> 8*i+2*j ) % 4;
			s += ObjexxFCL::lead_zero_string_of(tmp,1);
		}
		if(i+1 < nres) s += "-";
	}
	return s;
}

unsigned long pose2bin(core::pose::Pose const & pose) {
	using namespace ObjexxFCL::format;
	unsigned long bin = 0;
	for(int i = 0; i < min( 16, (int)pose.n_residue() ); ++i) {		
		Real phid = pose.phi(i+1);
		Real psid = pose.psi(i+1);
		// TR << phid << " " << psid << std::endl;
		unsigned long phi = (unsigned long)(basic::unsigned_periodic_range(phid,360.0) / 90.0);
		unsigned long psi = (unsigned long)(basic::unsigned_periodic_range(psid,360.0) / 90.0);
		phi = phi << 4*i;
		psi = psi << 4*i+2;
		bin += phi;
		bin += psi;
	}
	return bin;
}

Size get_aln_resi(core::pose::Pose const & pose, core::Size const aln_chain, char const aln_termi) {
	Size resi = 0, nnt = 0, nct = 0;
	for(Size i = 1; i <= pose.n_residue(); ++i) {
		if(pose.residue(i).is_lower_terminus()) {
			nnt++;
			if(aln_termi=='N' && nnt==aln_chain) {
				resi = i;
				break;
			}
		}
		if(pose.residue(i).is_upper_terminus()) {
			nct++;
			if(aln_termi=='C' && nnt==aln_chain) {
				resi = i;
				break;
			}
		}
	}
	if( 0 == resi ) utility_exit_with_message("ImplicitClashData: could not find aln resi in input pose chain: "+ObjexxFCL::string_of(aln_chain)+" term: "+aln_termi);
	return resi;
}
core::kinematics::Stub getxform(core::conformation::Residue const & move_resi, core::conformation::Residue const & fixd_resi) {
	core::kinematics::Stub s;
	s.M = alignVectorSets(move_resi.xyz(1)-move_resi.xyz(2),move_resi.xyz(3)-move_resi.xyz(2),fixd_resi.xyz(1)-fixd_resi.xyz(2),fixd_resi.xyz(3)-fixd_resi.xyz(2));
	s.v = fixd_resi.xyz(2)-s.M*move_resi.xyz(2);
	return s;
}
core::kinematics::Stub getxform(core::conformation::ResidueCOP      move_resi, core::conformation::Residue const & fixd_resi) { return getxform(*move_resi, fixd_resi); }
core::kinematics::Stub getxform(core::conformation::Residue const & move_resi, core::conformation::ResidueCOP      fixd_resi) { return getxform( move_resi,*fixd_resi); }
core::kinematics::Stub getxform(core::conformation::ResidueCOP      move_resi, core::conformation::ResidueCOP      fixd_resi) { return getxform(*move_resi,*fixd_resi); }


core::pose::Pose build_algned_linker(core::pose::Pose const & alnpose, Size len, core::chemical::ResidueTypeSetCAP resset, Pose const & oldlnk) {
	Pose lnk;
	string ggs = "GGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGS";
	make_pose_from_sequence( lnk , "A"+ggs.substr(0,len)+"A" , *resset , true);
	for(Size i = 1; i <= lnk.n_residue(); ++i) {
		lnk.set_phi  (i,-60);
		lnk.set_psi  (i,-45);
		lnk.set_omega(i,180);
	}
	FoldTree ft = lnk.fold_tree();
	if(ft.num_jump() != 0) utility_exit_with_message("why is there a jump in the fold tree alreaduy?!?!?!");
	ObjexxFCL::FArray2D_int jump_point(2,1);
	ObjexxFCL::FArray1D_int cuts(1);
	cuts(1) = lnk.n_residue()/2;
	core::pose::add_variant_type_to_pose_residue(lnk,"CUTPOINT_LOWER",cuts(1)+0);		
	core::pose::add_variant_type_to_pose_residue(lnk,"CUTPOINT_UPPER",cuts(1)+1);
	jump_point(1,1) = 1;
	jump_point(2,1) = lnk.n_residue();
	ft.tree_from_jumps_and_cuts( lnk.n_residue(), 1, jump_point, cuts );
	ft.set_jump_atoms(1,"N","C");
	lnk.fold_tree(ft);		

	xform_pose(lnk, getxform(lnk.residue(     1         ),alnpose.residue(1)) );
	xform_pose(lnk, getxform(lnk.residue(lnk.n_residue()),alnpose.residue(2)), cuts(1)+1 ); // only after cutpoint

	if( oldlnk.n_residue() != 0 && oldlnk.n_residue()+1 == lnk.n_residue()) {
		Size c = oldlnk.fold_tree().cutpoint(1);
		for(Size i = 1; i <= c; ++i) {
			lnk.set_phi(i, oldlnk.phi(i) );
			lnk.set_psi(i, oldlnk.psi(i) );
		}
		for(Size i = 0; i < oldlnk.n_residue()-c; ++i) {
			lnk.set_phi(lnk.n_residue()-i, oldlnk.phi(oldlnk.n_residue()-i) );
			lnk.set_psi(lnk.n_residue()-i, oldlnk.psi(oldlnk.n_residue()-i) );
		}
		// lnk.dump_pdb("lnk.pdb");
		// oldlnk.dump_pdb("oldlnk.pdb");	
	}	
	return lnk;
}


void* doit(void*) {
	using namespace core;
	using namespace chemical;
	using namespace conformation;
	using namespace pose;
	using namespace protocols;
	using namespace moves;
	using namespace ObjexxFCL::format;
	using numeric::random::uniform;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	ResidueTypeSetCAP cenresset = ChemicalManager::get_instance()->residue_type_set( CENTROID );

	Pose lnk;
	Size linklen1 = option[rblinker::linker1_size]();
	Size linklen2 = option[rblinker::linker2_size]();
	string seq1, seq2;
	string ggs = "GGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGSGGS";
	string lll = "***************************************************************************************";
	seq1 = "A" + ggs.substr(0,linklen1) + "A";
	seq2 = "A" + ggs.substr(0,linklen2) + "A";
	string ss1 = "_" + lll.substr(0,linklen1) + "_";
	string ss2 = "_" + lll.substr(0,linklen2) + "_";

	TR << "linker1 size " << linklen1 << " seq1: " << seq1 << std::endl;
	TR << "linker2 size " << linklen2 << " seq2: " << seq2 << std::endl;

	make_pose_from_sequence(lnk,seq1,*cenresset,true);
	for(Size i = 1; i <= lnk.n_residue(); ++i) {
		lnk.set_phi  (i,-60);
		lnk.set_psi  (i,-45);
		lnk.set_omega(i,180);
	}

	// if(linklen2 > 0) {
	// 	Pose tmp;
	// 	make_pose_from_sequence(tmp,seq2,*cenresset,true);
	// 	lnk.append_residue_by_jump(tmp.residue(1),lnk.n_residue(),"C","N"); // is this working????
	// 	for(Size i = 2; i <= tmp.n_residue(); ++i) {
	// 		lnk.append_residue_by_bond(tmp.residue(i));
	// 	}
	// 	for(Size i = 1; i <= lnk.n_residue(); ++i) {
	// 		lnk.set_phi  (i,-60);
	// 		lnk.set_psi  (i,-45);
	// 		lnk.set_omega(i,180);
	// 	}
	// 	
	// 	utility::io::izstream in( basic::options::option[basic::options::OptionKeys::score::fastclash::implicit_clash_aln_config]() );
	// 	Size c1,c2,c3,c4;
	// 	char t1,t2,t3,t4;
	// 	std::string fn;
	// 	Size count = 0;
	// 	while( in >> c1 >> t1 >> c2 >> t2 >> fn >> c3 >> t3 >> c4 >> t4 ) {			
	// 		count++;
	// 	}
	// 	if(count!=1) utility_exit_with_message("require exactly one line in fastclash aln config");
	// 	if( c1 != 1 || t1 != 'C' || t2 != 'N' || t3 != 'C' || t4 != 'N' || c4 != 2 ) utility_exit_with_message("improper fastclash aln config for 2 linker system!");
	// 	Pose alnpose;
	// 	core::import_pose::pose_from_pdb(alnpose,fn);
	// 	Size aln_resi_lnk_ref = get_aln_resi(lnk,c1,t1);
	// 	Size aln_resi_lnk_mov = get_aln_resi(lnk,c4,t4);
	// 	Size aln_resi_other_ref = get_aln_resi(alnpose,c3,t3);
	// 	Size aln_resi_other_mov = get_aln_resi(alnpose,c2,t2);
	// 	TR << "aligning two linkers lnkref " << aln_resi_lnk_ref << " lnkmov " << aln_resi_lnk_mov << " oref " << aln_resi_other_ref << " omov " << aln_resi_other_mov << std::endl;
	// 	Stub aln_alnpose = getxform( alnpose.residue(aln_resi_other_mov), lnk.residue(aln_resi_lnk_ref) );
	// 	xform_pose(alnpose,aln_alnpose);
	// 	if(option[rblinker::debug]()) alnpose.dump_pdb("init_alnpose.pdb");
	// 	Stub aln_chain2 = getxform( lnk.residue(aln_resi_lnk_mov), alnpose.residue(aln_resi_other_ref) );
	// 	{
	// 		for(Size ir = aln_resi_lnk_mov; ir <= lnk.n_residue(); ++ir) {
	// 			for(Size ia = 1; ia <= lnk.residue_type(ir).natoms(); ++ia) {
	// 				core::id::AtomID const aid(core::id::AtomID(ia,ir));
	// 				lnk.set_xyz( aid, aln_chain2.local2global(lnk.xyz(aid)) );
	// 			}
	// 		}
	// 	}
	// 	
	// 	// init test...
	// 	protocols::scoring::methods::ImplicitClashEnergy::dump_implicit_structures("imp0.pdb",lnk);	lnk.dump_pdb("pose0.pdb"); lnk.set_phi(10,uniform()*360.0);
	// 	protocols::scoring::methods::ImplicitClashEnergy::dump_implicit_structures("imp1.pdb",lnk);	lnk.dump_pdb("pose1.pdb"); lnk.set_phi(30,uniform()*360.0); 
	// 	protocols::scoring::methods::ImplicitClashEnergy::dump_implicit_structures("imp2.pdb",lnk);	lnk.dump_pdb("pose2.pdb"); lnk.set_phi(22,uniform()*360.0); 
	// 	protocols::scoring::methods::ImplicitClashEnergy::dump_implicit_structures("imp3.pdb",lnk); lnk.dump_pdb("pose3.pdb"); lnk.set_psi(22,uniform()*360.0); 
	// 	protocols::scoring::methods::ImplicitClashEnergy::dump_implicit_structures("imp4.pdb",lnk);	lnk.dump_pdb("pose4.pdb"); lnk.set_phi(23,uniform()*360.0); 
	// 	protocols::scoring::methods::ImplicitClashEnergy::dump_implicit_structures("imp5.pdb",lnk);	lnk.dump_pdb("pose5.pdb"); lnk.set_psi(23,uniform()*360.0); 
	// 	protocols::scoring::methods::ImplicitClashEnergy::dump_implicit_structures("imp6.pdb",lnk);	lnk.dump_pdb("pose6.pdb"); 
	// 	utility_exit_with_message("3 RB 2 linker systems not working properly yet!");
	// }

	core::scoring::ScoreFunctionOP sf = new core::scoring::ScoreFunction;
	sf->set_weight(core::scoring::vdw  ,5.0);
	sf->set_weight(core::scoring::pair ,1.0);
	sf->set_weight(core::scoring::cbeta,1.0);
	sf->set_weight(core::scoring::rama ,1.0);
	sf->set_weight(core::scoring::fastclash,1.0);
	sf->set_weight(core::scoring::implicitdock,option[rblinker::linearscore_wt]());

	bool CHAINBREAK = false;
	if(option[in::file::s].user()) {
		utility::vector1<std::string> fns = option[in::file::s]();		
		sf->set_weight(core::scoring::linear_chainbreak,20.0);
		for(Size ifn = 1; ifn <= fns.size(); ++ifn) {
			std::string fn = fns[ifn];
			// if(fns.size() != 1) utility_exit_with_message("rblinker2 currently only accepts 1 in:file:s input");
			Pose alnpose;
			core::import_pose::pose_from_pdb(alnpose,fn);
			if(alnpose.n_residue() != 2) utility_exit_with_message("reference pose must have only two residues");
			
			try {			
				Pose start;
				for(Size LNKLEN = 2; LNKLEN < 53; LNKLEN++) {
					start = build_algned_linker(alnpose,LNKLEN,cenresset,start);
					Pose best = start;
					bool linkable = false;
					for(Size CBTRIAL = 1; CBTRIAL <= 10; CBTRIAL++) {
						lnk = start;			
						if(option[rblinker::debug]()) TR << "try to find starting config without clash...." << std::endl;
						MoverOP bigbbmove = new SimpleBBMover(1,lnk.n_residue(),9e6);
						for(int precount = 1; precount < option[rblinker::nclashtrials](); precount++ ) {
							bigbbmove->apply(lnk);
							if( sf->score(lnk) < 9e5) break;
						}
						Pose cbtrial_best = lnk;
						if(lnk.energies().total_energies()[core::scoring::fastclash] > 9e5) {
							throw std::string("EVER-CLASH");
							utility_exit_with_message("failed to find starting pos w/o clash! nclashtrial "+string_of(option[rblinker::nclashtrials]()));							
						}
						SimpleBBMoverOP bbmove = new SimpleBBMover(1,lnk.n_residue(),option[rblinker::ntrials]()*2-20);
						sf->score(cbtrial_best);
						for(int ITER = 1; ITER <= option[rblinker::ntrials](); ITER++) {
							lnk = cbtrial_best;
							MonteCarloOP mc1 = new MonteCarlo( lnk, *sf, 2.0 );
							mc1->set_autotemp( true, 2.0 ); mc1->set_temperature( 2.0 );
							TrialMover trial(bbmove,mc1);
							for(int j = 1; j <= option[rblinker::nsubtrials](); j++) {
								trial.apply(lnk);
								if( lnk.energies().total_energies()[core::scoring::linear_chainbreak] < cbtrial_best.energies().total_energies()[core::scoring::linear_chainbreak] ) {
									cbtrial_best = lnk;
									best = lnk;
									if( lnk.energies().total_energies()[core::scoring::linear_chainbreak] < 1.0 ) {
										linkable = true;
										break;
									}
								}
							}
							if(linkable) {
								break;
							} else if( ITER > 10 && best.energies().total_energies()[core::scoring::linear_chainbreak] > 10.0 ) {
								break; // hopeless... don't keep trying at this linker length
							}
							Real m = max( 1.0, bbmove->magnitude()-1.0 );
							if(option[rblinker::debug]()) TR << "len " << LNKLEN << " setting bbmove mag to " << m << std::endl;
							bbmove->magnitude( m );
							if(option[rblinker::debug]()) TR << "len " << LNKLEN << " cb: " << lnk.energies().total_energies()[core::scoring::linear_chainbreak] << " " << cbtrial_best.energies().total_energies()[core::scoring::linear_chainbreak] << std::endl;
						}
						if(linkable) {
							lnk = best;
							break;
						} else if( best.energies().total_energies()[core::scoring::linear_chainbreak] > 3.0 ) {
							break; // hopeless... don't keep trying at this linker length
						}
					}
					if(linkable) {
						TR << "input configuration " << fn << " is linkable with linker size " << lnk.n_residue()-2 << std::endl;
						std::string fnout = option[out::file::o]()+"/"+utility::file_basename(fn)+"_linker_"+ObjexxFCL::lead_zero_string_of(lnk.n_residue()-2,3)+".pdb.gz";
						protocols::scoring::methods::ImplicitClashEnergy::dump_implicit_structures(fnout,lnk);
						break;
					} else {
						TR << "input configuration " << fn << " seems unlinkable with linker size " << lnk.n_residue()-2 << std::endl;
					}
				}
			} catch(std::string s) {
				if( s == "EVER-CLASH" ) {
					TR << "EVER-CLASH on file " << fn << std::endl;
				} else {
					utility_exit_with_message(s);
				}
			}
		}
	} else {
		protocols::moves::MoverOP bbmove;
		string const bb_samp_method = option[ rblinker::bb_samp_method ]();
		/**/ if( "simplegaussian1"  == bb_samp_method ) { bbmove = new SimpleBBMover(1,lnk.n_residue(), 1.0); TR << "bb_samp_method bbg1"    << std::endl; }
		else if( "simplegaussian5"  == bb_samp_method ) { bbmove = new SimpleBBMover(1,lnk.n_residue(), 5.0); TR << "bb_samp_method bbg5"    << std::endl; }
		else if( "simplegaussian10" == bb_samp_method ) { bbmove = new SimpleBBMover(1,lnk.n_residue(),10.0); TR << "bb_samp_method bbg10"   << std::endl; }
		else if( "simplegaussian20" == bb_samp_method ) { bbmove = new SimpleBBMover(1,lnk.n_residue(),20.0); TR << "bb_samp_method bbg20"   << std::endl; }
		else if( "simplegaussian30" == bb_samp_method ) { bbmove = new SimpleBBMover(1,lnk.n_residue(),30.0); TR << "bb_samp_method bbg30"   << std::endl; }
		else if( "simplegaussian60" == bb_samp_method ) { bbmove = new SimpleBBMover(1,lnk.n_residue(),60.0); TR << "bb_samp_method bbg60"   << std::endl; }
		else if( "uniform"          == bb_samp_method ) { bbmove = new SimpleBBMover(1,lnk.n_residue(), 9e6); TR << "bb_samp_method uniform" << std::endl; }
		else if( "bbg8t3a"          == bb_samp_method ) { bbmove = new protocols::simple_moves::BBG8T3AMover();      TR << "bb_samp_method bbg8t3a" << std::endl; }
		else if( "jump"             == bb_samp_method ) {
			TR << "bb_samp_method jump" << std::endl;
			FoldTree ft = lnk.fold_tree();
			if(ft.num_jump() == 0) {
				ObjexxFCL::FArray2D_int jump_point(2,1);
				ObjexxFCL::FArray1D_int cuts(1);
				cuts(1) = 1;
				jump_point(1,1) = 1;
				jump_point(2,1) = lnk.n_residue();
				ft.tree_from_jumps_and_cuts( lnk.n_residue(), 1, jump_point, cuts );
				lnk.fold_tree(ft);
			}
			bbmove = new rigid::RigidBodyPerturbMover(1, 10.0, 1.0);
		}
		else utility_exit_with_message("unknown option value for -rblinker:bb_samp_methd");

		TR << "try to find starting config without clash...." << std::endl;
		MoverOP bigbbmove = new SimpleBBMover(1,lnk.n_residue(),9e6);
		if( "jump" == bb_samp_method ) bigbbmove = bbmove;
		for(int precount = 1; precount < option[rblinker::nclashtrials](); precount++ ) {
			bigbbmove->apply(lnk);
			if(option[rblinker::debug]()) lnk.dump_pdb("test_"+ObjexxFCL::lead_zero_string_of(precount,7)+"_0.pdb");
			if( sf->score(lnk) < 9e5) break;
			if(option[rblinker::debug]()) lnk.dump_pdb("test_"+ObjexxFCL::lead_zero_string_of(precount,7)+"_1.pdb");
		}
		if(lnk.energies().total_energies()[core::scoring::fastclash] > 9e5) utility_exit_with_message("failed to find starting pos w/o clash! nclashtrial "+string_of(option[rblinker::nclashtrials]()));
		TR << "found starting config without clash...." << std::endl;

		utility::vector1<Size> hist(100,0);
		std::set<unsigned long> coverage;
		Size accepts = 0;
		for(int ITER = 1; ITER <= option[rblinker::ntrials](); ITER++) {
			// TR << "ITER " << ITER << std::endl;	
	
			MonteCarloOP mc1 = new MonteCarlo( lnk, *sf, 2.0 );
			mc1->set_autotemp( true, 2.0 ); mc1->set_temperature( 2.0 );
			TrialMover trial(bbmove,mc1);
		
			for(int j = 1; j <= option[rblinker::nsubtrials](); j++) {
				trial.apply(lnk);
				coverage.insert( pose2bin(lnk) );
				Real d = lnk.energies().total_energies()[core::scoring::implicitdock];
				hist[numeric::max(numeric::min(100,int(d)),1)]++;
			}
			accepts += trial.num_accepts();
		
			if(lnk.energies().total_energies()[core::scoring::fastclash] > 9e5) utility_exit_with_message("has clash!!!! exiting!");
		
			TR << "ITER: " << ITER*10000 << " accepts: " << accepts << " coverage: " << coverage.size() << " HIST ";
			for(Size i = 1; i <= hist.size(); ++i) {
				TR << hist[i] << " ";
			}
			TR << "HISTEND" << std::endl;
		
		
			Real d = lnk.energies().total_energies()[core::scoring::implicitdock];
			std::string  fn = "rblinker2_out_"+ObjexxFCL::lead_zero_string_of(ITER,9)+"_unbnd.pdb";
			if(d < 15.0) {
				int tmp = (int)	d;
				fn = option[out::file::o]()+"/"+"rblinker2_out_"+ObjexxFCL::lead_zero_string_of(ITER,9)+"_bound"+ObjexxFCL::lead_zero_string_of(tmp,2)+".pdb.gz";
			}
			protocols::scoring::methods::ImplicitClashEnergy::dump_implicit_structures(fn,lnk);

		}
	}	

	return NULL;
}



int main( int argc, char * argv [] ) {

	try {

	devel::init(argc,argv);

	void* (*func)(void*) = &doit;
	if (option[ basic::options::OptionKeys::parser::view ]()) {
		protocols::viewer::viewer_main( func );
	} else {
		func(NULL);
	}


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

}

/*

lnk
1560 runs OK
1691 runs OK
1692 runs OK
1756 runs OK

lnk2
1560 runs OK
1691 runs OK
1692 runs OK
1756 runs OK


*/
