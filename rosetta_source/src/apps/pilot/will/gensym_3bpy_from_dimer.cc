// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file /src/apps/pilat/will/genmatch.cc
/// @brief ???

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/smhybrid.OptionKeys.gen.hh>
#include <basic/options/keys/willmatch.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
#include <core/pack/optimizeH.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.hh>
#include <numeric/conversions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/moves/symmetry/SymMinMover.hh>
#include <protocols/moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/scoring/ImplicitFastClashCheck.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
// #include <devel/init.hh>

// #include <core/scoring/constraints/LocalCoordinateConstraint.hh>
#include "apps/pilot/will/will_util.hh"
#include "mynamespaces.hh"

using core::kinematics::Stub;
using protocols::scoring::ImplicitFastClashCheck;

static basic::Tracer TR("gensym_3bpy_from_dimer");
static core::io::silent::SilentFileData sfd;


inline Real const sqr(Real const r) { return r*r; }
inline Real sigmoidish_neighbor( Real const & sqdist ) {
  if( sqdist > 9.*9. ) {
    return 0.0;
  } else if( sqdist < 6.*6. ) {
    return 1.0;
  } else {
    Real dist = sqrt( sqdist );
    return sqr(1.0  - sqr( (dist - 6.) / (9. - 6.) ) );
  }
}


Real iface_check_c3(Pose & pose, Size nres, vector1<Size> const & iface_candidates) {
  Real num = 0;
  for(vector1<Size>::const_iterator i=iface_candidates.begin(),ie=iface_candidates.end(); i != ie; ++i) {
    for(vector1<Size>::const_iterator j=iface_candidates.begin(),je=iface_candidates.end(); j != je; ++j) {
      num += sigmoidish_neighbor(pose.residue(*i).xyz(5).distance_squared(pose.residue(*j+1*nres).xyz(5)));
      num += sigmoidish_neighbor(pose.residue(*i).xyz(5).distance_squared(pose.residue(*j+2*nres).xyz(5)));
    }
  }
  return num;
}

vector1<Size> read_res_list(string fn) {
  vector1<Size> l;
  if(fn=="") return l;
  if(fn=="_") return l;
  if(fn.size()==1 && fn[0]==(char)0) return l;
  izstream in(fn);
  if(!in.good()) {
    utility_exit_with_message("can't open res list file '"+fn+"'");
  }
  Size r;
  while( in >> r ) l.push_back(r);
  return l;
}

vector1<Vec> line_cone_intersection(Vec p, Vec d, Vec v, Vec a, Real t) {
  vector1<Vec> sol;
  t = numeric::conversions::radians(t);
  Mat M = numeric::outer_product(a,a) - cos(t)*cos(t)*Mat::identity();
  Vec D = p-v;
  Real c2 =   d.dot(M*d);
  Real c1 = 2*d.dot(M*D);
  Real c0 =   D.dot(M*D);
  Real disc = c1*c1 - 4*c0*c2;
  if( disc == 0) sol.push_back( p + (-c1)/(2.0*c2)*d );
  else if( disc > 0) {
    disc = sqrt(disc);
    sol.push_back(p+(-c1+disc)/(2.0*c2)*d);
    sol.push_back(p+(-c1-disc)/(2.0*c2)*d);
  }
  return sol;
}

// Vec projperp(Vec u, Vec v) {
//  return v - projection_matrix(u)*v;
// }

vector1<std::pair<Vec,Vec> > intersecting_bpy_axes(Vec CB, Vec CG, Vec FE, Vec symmaxis, Vec symmcen = Vec(0,0,0)) {
  vector1<std::pair<Vec,Vec> > sol;
  Vec  p = symmcen;
  Vec  d = symmaxis.normalized();
  Vec  a = (CG-CB).normalized();
  Vec  v = CG + projection_matrix(a)*(FE-CG);
  Real t = 35.2643434495;
  Real l = v.distance(FE);
  vector1<Vec> X = line_cone_intersection(p,d,v,a,t);
  for(Size i = 1; i <= X.size(); ++i) {
    Vec x = X[i];
    Vec o = projperp(a,x-v);
    Real L = o.length();
    o = o.normalized() * l;
    Real ang = 90.0 - numeric::conversions::degrees( atan(l/L) );
    Vec o1 = rotation_matrix_degrees(a, ang) * o;
    Vec o2 = rotation_matrix_degrees(a,-ang) * o;
    sol.push_back(std::pair<Vec,Vec>(x,v+o1));
    sol.push_back(std::pair<Vec,Vec>(x,v+o2));
  }
  return sol;
}



void fixH(core::pose::Pose & pose) {
  for(Size i = 1; i <= pose.n_residue(); ++i) {
    if(!pose.residue(i).has("H")) continue;
    numeric::xyzVector<Real> n  = pose.residue(i).xyz("N");
    numeric::xyzVector<Real> ca = pose.residue(i).xyz("CA");
    Size in = i-1;
    if(in == 0) in = pose.n_residue();
    numeric::xyzVector<Real> c  = pose.residue(in).xyz("C");
    numeric::xyzVector<Real> h  = n + (n-(ca+c)/2.0).normalized()*1.01;
    pose.set_xyz(AtomID(pose.residue(i).atom_index("H"),i), h );
  }
}

bool iface_check_C2Z(Pose const & p, Size rsd) {
  Mat Rc2 = rotation_matrix_degrees(Vec(0,0,1),180.0);
  Vec CA = p.xyz(AtomID(2,rsd));
  for (Size ir; ir < p.n_residue(); ++ir) {
    if(p.residue(ir).is_protein())
      if( CA.distance_squared( Rc2*p.xyz(AtomID(2,ir))) < 100.0 )
        return true;
  }
  return false;
}


void refine(Pose & pose, Size ibpy, Size dsub) {
	using namespace core::chemical;
	using namespace core::conformation::symmetry;
	using namespace core::pack::task;
	using namespace core::scoring::constraints;
	ScoreFunctionOP sf = core::scoring::getScoreFunction();
	sf->set_weight(core::scoring::atom_pair_constraint,1.0);
	sf->set_weight(core::scoring::coordinate_constraint,1.0);

  // Get the symmetry info and make the packer task
  SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);

	pose.add_constraint( new AtomPairConstraint( AtomID(pose.residue(ibpy).atom_index("ZN"),ibpy),
																							 AtomID(pose.residue(ibpy).atom_index("ZN"),ibpy+1*sym_info->num_independent_residues()),
																							 new HarmonicFunc(0,0.02) ) );
	pose.add_constraint( new AtomPairConstraint( AtomID(pose.residue(ibpy).atom_index("ZN"),ibpy),
																							 AtomID(pose.residue(ibpy).atom_index("ZN"),ibpy+2*sym_info->num_independent_residues()),
																							 new HarmonicFunc(0,0.02) ) );

	AtomID REF(1,sym_info->num_total_residues_without_pseudo()+1);
	for(Size ir = 1; ir <= sym_info->num_independent_residues(); ++ir) {
		if(ir==ibpy) continue;
		for(Size ia = 1; ia <= pose.residue(ir).nheavyatoms(); ++ia) {
			pose.add_constraint(new CoordinateConstraint( AtomID(ia,ir), REF, pose.xyz(AtomID(ia,ir)), new HarmonicFunc(0.0,1.0) ) );
		}
	}
  core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	movemap->set_jump(true);
  movemap->set_bb(true);
  movemap->set_chi(true);
	core::pose::symmetry::make_symmetric_movemap(pose,*movemap);
  protocols::moves::symmetry::SymMinMover( movemap, sf, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false ).apply(pose);

	pose.remove_constraints();

  // Set allowed AAs.
  vector1<bool> allowed_aas(20,false);
  allowed_aas[aa_ala] = true;
  allowed_aas[aa_asp] = true;
  allowed_aas[aa_phe] = true;
  allowed_aas[aa_ile] = true;
  allowed_aas[aa_leu] = true;
  allowed_aas[aa_met] = true;
  allowed_aas[aa_asn] = true;
  allowed_aas[aa_ser] = true;
  allowed_aas[aa_thr] = true;
  allowed_aas[aa_val] = true;
  allowed_aas[aa_trp] = true;
  allowed_aas[aa_tyr] = true;

  PackerTaskOP task( TaskFactory::create_packer_task( pose ));


	vector1<Size> design_pos;
  for( Size i=1; i<=pose.n_residue(); i++) {
		if(i==ibpy) continue;
		if(!sym_info->bb_is_independent(i)) continue;
		for( Size j=1; j<=pose.n_residue(); j++) {
			if(sym_info->bb_is_independent(j)) continue;
			if(sym_info->subunit_index(j) == dsub ) continue;
			if( pose.residue(i).nbr_atom_xyz().distance_squared( pose.residue(j).nbr_atom_xyz() ) < 100.0 ) {
				//TR << "design_pos " << i << " " << j << std::endl;
				if( find(design_pos.begin(), design_pos.end(), i) == design_pos.end()) design_pos.push_back(i);
			}
		}
	}

  // Set which residues can be designed
  for( Size i=1; i<=pose.n_residue(); i++) {
    if (!sym_info->bb_is_independent(i)) {
      task->nonconst_residue_task(i).prevent_repacking();
    } else if( pose.residue(i).name3() == "PRO" || pose.residue(i).name3() == "GLY") {
      // Don't mess with Pros or Glys at the interfaces
      task->nonconst_residue_task(i).prevent_repacking();
    } else if (find(design_pos.begin(), design_pos.end(), i) == design_pos.end()) {
      task->nonconst_residue_task(i).prevent_repacking();
    } else {
      bool temp = allowed_aas[pose.residue(i).aa()];
      allowed_aas[pose.residue(i).aa()] = true;
      task->nonconst_residue_task(i).restrict_absent_canonical_aas(allowed_aas);
      task->nonconst_residue_task(i).or_include_current(true);
      task->nonconst_residue_task(i).initialize_from_command_line();
      allowed_aas[pose.residue(i).aa()] = temp;
    }
  }

  // Actually perform design.
	core::pack::make_symmetric_PackerTask(pose, task);
  protocols::moves::MoverOP packer = new protocols::moves::symmetry::SymPackRotamersMover(sf, task);
  packer->apply(pose);

}


void repack(Pose & pose, Size ibpy, Size dsub) {
	using namespace core::chemical;
	using namespace core::conformation::symmetry;
	using namespace core::pack::task;
	using namespace core::scoring::constraints;
	ScoreFunctionOP sf = core::scoring::getScoreFunction();

  // Get the symmetry info and make the packer task
  SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);

  PackerTaskOP task( TaskFactory::create_packer_task( pose ));

	vector1<Size> design_pos;
  for( Size i=1; i<=pose.n_residue(); i++) {
		if(i==ibpy) continue;
		if(!sym_info->bb_is_independent(i)) continue;
		for( Size j=1; j<=pose.n_residue(); j++) {
			if(sym_info->bb_is_independent(j)) continue;
			if(sym_info->subunit_index(j) == dsub ) continue;
			if( pose.residue(i).nbr_atom_xyz().distance_squared( pose.residue(j).nbr_atom_xyz() ) < 100.0 ) {
				//TR << "design_pos " << i << " " << j << std::endl;
				if( find(design_pos.begin(), design_pos.end(), i) == design_pos.end()) design_pos.push_back(i);
			}
		}
	}

  // Set which residues can be designed
  for( Size i=1; i<=pose.n_residue(); i++) {
    if (!sym_info->bb_is_independent(i)) {
      task->nonconst_residue_task(i).prevent_repacking();
    } else if( pose.residue(i).name3() == "PRO" || pose.residue(i).name3() == "GLY") {
      // Don't mess with Pros or Glys at the interfaces
      task->nonconst_residue_task(i).prevent_repacking();
    } else if (find(design_pos.begin(), design_pos.end(), i) == design_pos.end()) {
      task->nonconst_residue_task(i).prevent_repacking();
    } else {
      task->nonconst_residue_task(i).restrict_to_repacking();
    }
  }

  // Actually perform design.
	core::pack::make_symmetric_PackerTask(pose, task);
  protocols::moves::MoverOP packer = new protocols::moves::symmetry::SymPackRotamersMover(sf, task);
  packer->apply(pose);

}

// Pose needs to be scored before this will work.
void
new_sc(Pose &pose, utility::vector1<Size> intra_subs, Real& int_area, Real& sc) {

  using namespace core;

  core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);
  core::scoring::sc::ShapeComplementarityCalculator scc;
  scc.Init();

  // Figure out which chains touch chain A, and add the residues from those chains
  // into the sc surface objects
  Size nres_monomer = symm_info->num_independent_residues();
  for (Size i=1; i<=nres_monomer; ++i) {
    scc.AddResidue(0, pose.residue(i));
  }
  for (Size i=1; i<=symm_info->subunits(); ++i) {
    if (std::find(intra_subs.begin(), intra_subs.end(), i) != intra_subs.end()) continue;
    bool contact = false;
    Size start = (i-1)*nres_monomer;
    for (Size ir=1; ir<=nres_monomer; ir++) {
      if (pose.energies().residue_total_energies(ir+start)[core::scoring::fa_atr] < 0) {
        contact = true;
        break;
      }
    }
    if (contact) {
      for (Size ir=1; ir<=nres_monomer; ir++) {
        scc.AddResidue(1, pose.residue(ir+start));
      }
    }
  }
  if (scc.Calc()) {
    sc = scc.GetResults().sc;
    int_area = scc.GetResults().surface[2].trimmedArea;
  }
}


Real packing_score(Pose const & psym) {
	
}


#define ATET 54.735610317245360079 // asin(sr2/sr3)
#define AOCT 35.264389682754668343 // asin(sr1/sr3)
#define AICS 20.89774264557        // asin(G/2/sr3)
#define ANGTH 1.0

void run() {

  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  ScoreFunctionOP sf = core::scoring::getScoreFunction();
  ScoreFunctionOP sfrep = new core::scoring::ScoreFunction;
  sfrep->set_weight(core::scoring::fa_rep,1.0);

  core::chemical::ResidueTypeSetCAP rs = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
  Pose bpy,ala;
  core::import_pose::pose_from_pdb(bpy ,*rs,"input/bpy_ideal.pdb");
  core::pose::remove_lower_terminus_type_from_pose_residue(bpy,1);
  core::pose::remove_upper_terminus_type_from_pose_residue(bpy,1);
  make_pose_from_sequence(ala,"A","fa_standard",false);
  remove_lower_terminus_type_from_pose_residue(ala,1);
  remove_upper_terminus_type_from_pose_residue(ala,1);

  Vec a2f1 = Vec(0,0,1);
  Vec c2f1 = Vec(0,0,0);
  Mat Rc2 = rotation_matrix_degrees(a2f1,180.0);

  vector1<string> const & fnames(option[OptionKeys::in::file::s]());
  TR << "gensym_3bpy_from_dimer " << fnames.size() << " files" << std::endl;
  for(Size ifile = 1; ifile <= fnames.size(); ++ifile) {
    string fname = fnames[ifile];
    string infile = utility::file_basename(fnames[ifile]);
    Pose nat;
    core::import_pose::pose_from_pdb(nat,*rs,fname);
    Pose base(nat);

		Size cyscnt = 0;
		for(Size ir = 1; ir <= base.n_residue(); ++ir) {
			if(base.residue(ir).name3()=="CYS") cyscnt++;
		}

		if(cyscnt > 3) continue;
    //TR << "gensym_3bpy_from_dimer " << ifile << " " << fnames[ifile] << " " << base.n_residue() << " residues" << " " << cyscnt << std::endl;
    for(Size ibpy = 1; ibpy <= base.n_residue(); ++ibpy) {
      if(!base.residue(ibpy).is_protein()) continue;
      if(iface_check_C2Z(base,ibpy)) continue;
      if(base.residue(ibpy).is_lower_terminus()) continue;
      if(base.residue(ibpy).is_upper_terminus()) continue;
      base = nat;
      base.replace_residue(ibpy,bpy.residue(1),true);
      Real chi1_incr = option[willmatch::chi1_increment]();
      for(Real bch1 = 0.0; bch1 <= 360; bch1 += chi1_incr) {
        base.set_chi(1,ibpy,bch1);
        for(Size ir = 1; ir <= base.n_residue(); ++ir) {
          Size natom = (ir==ibpy) ? 5 : base.residue(ir).nheavyatoms();
          for(Size ia = 1; ia <= natom; ia++) {
            if( base.xyz(AtomID(ia,ir)).distance_squared(base.residue(ibpy).xyz("CZ")) < 9.0 ) goto clash2;
            if( base.xyz(AtomID(ia,ir)).distance_squared(base.residue(ibpy).xyz("CP")) < 9.0 ) goto clash2;
            if( base.xyz(AtomID(ia,ir)).distance_squared(base.residue(ibpy).xyz("CM")) < 9.0 ) goto clash2;
          }
        }
        goto noclash2;  clash2: continue; noclash2:
        Vec const CB = base.residue(ibpy).xyz("CB");
        Vec const CG = base.residue(ibpy).xyz("CG");
        Vec const FE = base.residue(ibpy).xyz("ZN");
        vector1<std::pair<Vec,Vec> > baxes = intersecting_bpy_axes(CB,CG,FE,a2f1,c2f1);
        for(Size jaxs = 1; jaxs <= baxes.size(); ++jaxs) {
          Vec const isct = baxes[jaxs].first;
          Vec const c3f1 = baxes[jaxs].second;
          if( isct.distance_squared(c3f1) < 25.0 ) continue;
          if( isct.distance_squared(c2f1) < 25.0 ) continue;
          Vec const a3f1 = (isct-c3f1).normalized();
          Real const orig_ang = angle_degrees( c2f1, isct, c3f1);
          Real const ang = (orig_ang > 90.0) ? 180.0-orig_ang : orig_ang;
          if( fabs(ang-ATET) > ANGTH && fabs(ang-AOCT) > ANGTH && fabs(ang-AICS) > ANGTH ) continue;
          // aln bpy
          Vec v = CG + projection_matrix(CB-CG)*(FE-CG);
          Real const dang = dihedral_degrees( c3f1,v,CG,FE );
          base.set_chi(2,ibpy, base.chi(2,ibpy) + dang );
          for(Size ir = 1; ir <= base.n_residue(); ++ir) {
            if(ir==ibpy) continue;
            for(Size ia = 1; ia <= base.residue(ir).nheavyatoms(); ia++) {
              for(Size ja = 7; ja <= base.residue(ibpy).nheavyatoms(); ja++) {
                if( base.xyz(AtomID(ia,ir)).distance_squared(base.xyz(AtomID(ja,ibpy))) < 9.0 ) goto clash3;
              }
            }
          }
          goto noclash3;  clash3: continue; noclash3:

          Mat const R1 = rotation_matrix_degrees(a3f1,120.0);
          Mat const R2 = rotation_matrix_degrees(a3f1,240.0);
          for(Size ir = 1; ir <= base.n_residue(); ++ir) {
            for(Size ia = 1; ia <= min(5ul,base.residue(ir).nheavyatoms()); ia++) {
              if(ir==ibpy) {
                if(base.residue(ir).atom_name(ia)=="NE1" || base.residue(ir).atom_name(ia)=="NN1") continue;
              }
              Vec const x1 = R1*(base.xyz(AtomID(ia,ir))-c3f1)+c3f1;
              Vec const x2 = R2*(base.xyz(AtomID(ia,ir))-c3f1)+c3f1;
              for(Size jr = 1; jr <= base.n_residue(); ++jr) {
                for(Size ja = 1; ja <= min(5ul,base.residue(jr).nheavyatoms()); ja++) {
                  if( x1.distance_squared(base.xyz(AtomID(ja,jr))) < 9.0 ) goto clash4;
                  if( x1.distance_squared( Rc2*(base.xyz(AtomID(ja,jr))-c2f1)+c2f1 ) < 9.0 ) goto clash4;
                  if( x2.distance_squared( Rc2*(base.xyz(AtomID(ja,jr))-c2f1)+c2f1 ) < 9.0 ) goto clash4;

                }
              }
            }
          }
          goto noclash4;  clash4: continue; noclash4:

          Real const baserep = sfrep->score(base);
          Pose psym = base;
          trans_pose(psym,-isct);
          Pose psym_bare = psym;
          Mat Rsymm;
          string symtag;
          Vec f1=(c3f1-isct).normalized(),f2=(c2f1-isct).normalized(),t1=Vec(0,0,1);
          Size dimersub = 0;
          if       ( fabs(ang-ATET) <= ANGTH ) {
            symtag = "TET";
            option[OptionKeys::symmetry::symmetry_definition]("input/sym/tetra.sym");
            Rsymm = alignVectorSets(f1,f2,t1, (orig_ang<90.0?1.0:-1.0)*Vec(0.8164965743782284,0.0,0.5773502784520137) );
            dimersub =  4;
          } else if( fabs(ang-AOCT) <= ANGTH ) {
            symtag = "OCT";
            option[OptionKeys::symmetry::symmetry_definition]("input/sym/octa.sym");
            Rsymm = alignVectorSets(f1,f2,t1, (orig_ang<90.0?1.0:-1.0)*Vec(0.4082482904638630,0.4082482904638626,0.8164965809277260));
            dimersub =  4;
          } else if( fabs(ang-AICS) <= ANGTH ) {
            symtag = "ICS";
            option[OptionKeys::symmetry::symmetry_definition]("input/sym/icosa.sym");
            Rsymm = alignVectorSets(f1,f2,t1, (orig_ang<90.0?1.0:-1.0)*rotation_matrix_degrees(Vec(0,0,1),120.0)*Vec(0.35670090519235864157,0.0,0.93421863834701557305));
            dimersub =  11;
          } else {
						continue;
            TR << "closed symm from ang " << ang << " not yet supported" << std::endl;
          }

          rot_pose(psym     ,Rsymm); core::pose::symmetry::make_symmetric_pose(psym     );
          rot_pose(psym_bare,Rsymm); core::pose::symmetry::make_symmetric_pose(psym_bare);

          string fname = symtag+"_"+infile+"_B"+lead_zero_string_of(ibpy,3)+"-"+lead_zero_string_of((Size)(bch1),3)+"_"+lead_zero_string_of(jaxs,1)+".pdb";

          sf->score(psym_bare);
          // TR << "REP " << 2*baserep << " " << psym_bare.energies().total_energies()[core::scoring::fa_rep] << std::endl;
          if(! psym_bare.energies().total_energies()[core::scoring::fa_rep] - 2*baserep > 1000.0) continue;

          //refine(psym,ibpy,dimersub);
          psym.dump_pdb(option[OptionKeys::out::file::o]()+"/"+fname);

					Real rms = core::scoring::CA_rmsd(psym,psym_bare);
					Real sc,int_area;
					vector1<Size> intra_subs(1,4);
					sf->score(psym);
					//new_sc(psym,intra_subs,int_area,sc);


					//repack(psym,ibpy,dimersub);
					Real s0 = sf->score(psym);
					for(Size ir = 1; ir <= base.n_residue(); ++ir) {
						for(Size ia = 1; ia <= psym.residue_type(ir).natoms(); ++ia) {
							core::id::AtomID const aid(core::id::AtomID(ia,ir));
							if("ICS"==symtag) psym.set_xyz( aid, psym.xyz(aid) + 20*Vec(-0.178411, 0.309017, 0.934172 ) );
							if("OCT"==symtag)	psym.set_xyz( aid, psym.xyz(aid) + 20*Vec( 0.408248, 0.408248, 0.816497 ) );
							if("TET"==symtag)	psym.set_xyz( aid, psym.xyz(aid) + 20*Vec( 0.816497, 0.000000, 0.577350 ) );
						}
					}
					//repack(psym,ibpy,dimersub);
					Real s1 = sf->score(psym);
          psym.dump_pdb(option[OptionKeys::out::file::o]()+"/"+fname+"_DIMER.pdb");

          TR << "HIT " << fname << " " <<  ang << " " << ibpy << " " << bch1 << " " << jaxs << " ddG: " << s0-s1 << " sc: " << sc << " " << int_area <<std::endl;

          core::io::silent::SilentStructOP ss_out( new core::io::silent::ScoreFileSilentStruct );
          ss_out->fill_struct(psym_bare,fname);
					ss_out->add_energy("rms",rms);
					ss_out->add_energy("ddg",s0-s1);
					ss_out->add_energy("sc",sc);
					ss_out->add_energy("int_area",int_area);
          sfd.write_silent_struct( *ss_out, option[OptionKeys::out::file::o]()+"/"+option[basic::options::OptionKeys::out::file::silent ]() );

        }
      }
    }
  }
  /*

    string fname = infile+"_"+lead_zero_string_of(irsd,3)+"_"+lead_zero_string_of((Size)(chi1+180.0),3)+"_"+lead_zero_string_of((Size)(chi2+180.0),3)+"_"+lead_zero_string_of(jrsd,3)+"_"+lead_zero_string_of((Size)(totrot),3)+"_"+lead_zero_string_of(iaxs,3)+".pdb";

    sf->score(psym);
    if(psym.energies().total_energies()[core::scoring::rg] > 14.0) {
    TR << "DANGER DANGER DANGER!!!" << std::endl;
    continue;
    }
    psym.dump_pdb(fname);

    core::pose::replace_pose_residue_copying_existing_coordinates(psym,jrsd,rs->name_map("ALA"));
    sf->score(psym);
    core::io::silent::SilentStructOP ss_out( new core::io::silent::ScoreFileSilentStruct );
    ss_out->fill_struct(psym,fname);
    sfd.write_silent_struct( *ss_out, option[ basic::options::OptionKeys::out::file::silent ]() );

    }

    }



    }
    // utility_exit_with_message("debug SG");
    } else {
    //matches?;
    }


    }
    }
  */
}

int main (int argc, char *argv[]) {

  // Vec p(-6.456746, 5.922204, -0.982538);
  // Vec d(0.393718,  0.677101,  0.621707);
  // Vec v(0.000000,  0.000000,  0.000000);
  // Vec a(0.998233, -0.003844, -0.059301);
  // Real t = 45;
  // vector1<Vec> X = line_cone_intersection(p,d,v,a,t);
  // for(Size i = 1; i <= X.size(); ++i) {
  //  TR << "X " << i << " " << X[i] << std::endl;
  // }
  // utility_exit_with_message("debug line_cone_intersection");

  core::init(argc,argv);
  run();
}




//
//








