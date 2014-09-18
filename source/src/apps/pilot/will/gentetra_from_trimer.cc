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
#include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/util.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/util.hh>
// AUTO-REMOVED #include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymDof.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymmData.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymmetryInfo.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
// AUTO-REMOVED #include <core/pack/optimizeH.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
// AUTO-REMOVED #include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
// AUTO-REMOVED #include <core/scoring/Energies.hh>
// AUTO-REMOVED #include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
// AUTO-REMOVED #include <core/scoring/ScoringManager.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <numeric/conversions.hh>
// AUTO-REMOVED #include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
// AUTO-REMOVED #include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/scoring/ImplicitFastClashCheck.hh>
#include <sstream>
#include <utility/io/izstream.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>
// #include <devel/init.hh>

// #include <core/scoring/constraints/LocalCoordinateConstraint.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>
#include <apps/pilot/will/mynamespaces.ihh>
#include <apps/pilot/will/will_util.ihh>


using core::kinematics::Stub;
using protocols::scoring::ImplicitFastClashCheck;

static thread_local basic::Tracer TR( "gentetra" );
static core::io::silent::SilentFileData sfd;


inline Real sqr(Real const r) { return r*r; }
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

void repack(Pose & pose, Size nres, ScoreFunctionOP sf) {
  using namespace core::pack::task;
  PackerTaskOP task = TaskFactory::create_packer_task(pose);
  task->initialize_extra_rotamer_flags_from_command_line();
  for(Size i = 1; i <= nres; ++i) {
    if(pose.residue(i).name3()=="BPY") {
      task->nonconst_residue_task(i).prevent_repacking();
    } else {
      task->nonconst_residue_task(i).restrict_to_repacking();
    }
  }
  // TR << *task << std::endl;
  protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
  repack.apply(pose);
}

void design(Pose & pose, Size nres, ScoreFunctionOP sf) {
  core::id::AtomID_Map< bool > atom_map;
  core::pose::initialize_atomid_map( atom_map, pose, false );
  for ( Size ir = 1; ir <= pose.total_residue(); ++ir ) {
    atom_map.set(AtomID(2,ir) , true );
    atom_map.set(AtomID(3,ir) , true );
    atom_map.set(AtomID(5,ir) , true );
  }
  core::id::AtomID_Map<Real> atom_sasa; utility::vector1<Real> sasa;
  core::scoring::calc_per_atom_sasa( pose, atom_sasa, sasa, 2.3, false, atom_map );
  for(Size i = 1; i <= sasa.size(); ++i) if( atom_sasa.n_atom(i) > 4 ) sasa[i] = atom_sasa[AtomID(5,i)];

  using namespace core::pack::task;
  PackerTaskOP task = TaskFactory::create_packer_task(pose);
  vector1< bool > aas(20,true);
  aas[core::chemical::aa_cys] = false;
  aas[core::chemical::aa_his] = false;
  // aas[core::chemical::aa_met] = false;
  aas[core::chemical::aa_pro] = false;
  aas[core::chemical::aa_gly] = false;
  aas[core::chemical::aa_gly] = false;
  if(option[basic::options::OptionKeys::willmatch::exclude_ala]()) aas[core::chemical::aa_ala] = false;
  if(option[basic::options::OptionKeys::smhybrid::design_hydrophobic]()) {
    aas[core::chemical::aa_ser] = false;
    aas[core::chemical::aa_thr] = false;
    aas[core::chemical::aa_asp] = false;
    aas[core::chemical::aa_glu] = false;
    aas[core::chemical::aa_lys] = false;
    aas[core::chemical::aa_arg] = false;
    aas[core::chemical::aa_asn] = false;
    aas[core::chemical::aa_gln] = false;
  }

  vector1<Size> fixed;
  if(option[basic::options::OptionKeys::willmatch::fixed_res].user()) {
    utility::io::izstream in(option[basic::options::OptionKeys::willmatch::fixed_res]());
    Size tmp;
    while(in>>tmp) fixed.push_back(tmp);
    in.close();
  }

  vector1<Size> interface;
  if(option[basic::options::OptionKeys::willmatch::design_interface]()) {
    for(Size i = 1; i <= nres; ++i) {
      if( sasa.size() >= i && sasa[i] > 15.0 ) continue;
      AtomID aid(5,i);
      if(pose.residue(i).nheavyatoms() < 5) aid.atomno() = 2;
      for(Size j = nres+1; j <= 3*nres; ++j) {
        AtomID aid2(5,j);
        if(pose.residue(j).nheavyatoms() < 5) aid.atomno() = 2;
        if(pose.xyz(aid).distance_squared(pose.xyz(aid2)) < 49) {
          interface.push_back(i);
        }
      }
    }
    for(Size i = 1; i <= nres; ++i) {
      Vec xyz = pose.xyz(AtomID(5,i));
      xyz.z() = 0;
      if( xyz.length() < 8.0 ) interface.push_back(i);
    }
  }
  for(Size i = 1; i <= nres; ++i) {
    if(pose.residue(i).name3()=="BPY") {
      task->nonconst_residue_task(i).prevent_repacking();
      // std::exit(-1);
    } else if(std::find(fixed.begin(),fixed.end(),i)!=fixed.end()){
      task->nonconst_residue_task(i).restrict_to_repacking();
      task->nonconst_residue_task(i).or_ex1_sample_level( core::pack::task::EX_ONE_STDDEV );
      task->nonconst_residue_task(i).or_ex2_sample_level( core::pack::task::EX_ONE_STDDEV );
    } else if(std::find(interface.begin(),interface.end(),i)!=interface.end()){
      task->nonconst_residue_task(i).restrict_absent_canonical_aas(aas);
      task->nonconst_residue_task(i).initialize_extra_rotamer_flags_from_command_line();
    } else {
      task->nonconst_residue_task(i).restrict_to_repacking();
      task->nonconst_residue_task(i).or_ex1_sample_level( core::pack::task::EX_ONE_STDDEV );
      task->nonconst_residue_task(i).or_ex2_sample_level( core::pack::task::EX_ONE_STDDEV );
    }
  }
  for(Size i = nres+1; i <= pose.n_residue(); ++i) {
    task->nonconst_residue_task(i).prevent_repacking();
  }
  TR << *task << std::endl;

  protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
  repack.apply(pose);
}

vector1<Vec> line_cone_intersection(Vec p, Vec d, Vec v, Vec a, Real t) {
  using namespace numeric;
  vector1<Vec> sol;
  t = conversions::radians(t);
  Mat M = outer_product(a,a) - cos(t)*cos(t)*Mat::identity();
  Vec D = p-v;
  core::Real c2 =   d.dot(M*d);
  core::Real c1 = 2*d.dot(M*D);
  core::Real c0 =   D.dot(M*D);
  core::Real disc = c1*c1 - 4*c0*c2;
  if( disc == 0) sol.push_back( p + (-c1)/(2.0*c2)*d );
  else if( disc > 0) {
    disc = sqrt(disc);
    sol.push_back(p+(-c1+disc)/(2.0*c2)*d);
    sol.push_back(p+(-c1-disc)/(2.0*c2)*d);
  }
  return sol;
}

vector1<std::pair<Vec,Vec> > intersecting_disulfide_axes(Vec CB, Vec SG, Vec symmaxis, Vec symmcen = Vec(0,0,0)) {
  vector1<std::pair<Vec,Vec> > sol;
  Vec p = symmcen;
  Vec d = symmaxis.normalized();
  Vec v = SG;
  Vec a = (SG-CB).normalized();
  Real t = 45.0; // dihedral should be 90 = 2*45
  vector1<Vec> X = line_cone_intersection(p,d,v,a,t);
  for(Size i = 1; i <= X.size(); ++i) {
    Vec x = X[i];
    Vec o = projperp(a,x-v);
    Real L = o.length();
    o = o.normalized();
    Real ang = 90.0 - numeric::conversions::degrees( atan(1.0/L) );
    Vec o1 = rotation_matrix_degrees(a, ang) * o;
    Vec o2 = rotation_matrix_degrees(a,-ang) * o;
    sol.push_back(std::pair<Vec,Vec>(x,v+o1));
    sol.push_back(std::pair<Vec,Vec>(x,v+o2));
  }
  return sol;
}



void minimize(Pose & pose, Size nres, Size , ScoreFunctionOP sf, int bb=0) {
  core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
  // core::pose::symmetry::make_symmetric_movemap(pose,*movemap);
  movemap->set_chi(true);
  movemap->set_bb(false);
  movemap->set_jump(false);
  if(bb==1) for(Size i = 1; i <= nres; ++i) if(pose.secstruct(i)=='L') movemap->set_bb(i,true);
  if(bb>=2) for(Size i = 1; i <= nres; ++i) movemap->set_bb(i,true);
  // movemap->set_chi(bpyres,false);

  core::pose::symmetry::make_symmetric_movemap( pose, *movemap );

  protocols::simple_moves::symmetry::SymMinMover m( movemap, sf, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false );

  m.apply(pose);

}

std::pair<Size,Size> makesplitwork_3bpy(Size total) {
  using namespace basic::options::OptionKeys;
  Size size1 = 1;
  Size size2 = total;
  if( option[willmatch::splitwork].user() ) {
    Size part   = option[willmatch::splitwork]()[1];
    Size nparts = option[willmatch::splitwork]()[2];
    size1 = (part-1)*(Size)std::ceil(((Real)total)/(Real)nparts)+1;
    size2 = (part  )*(Size)std::ceil(((Real)total)/(Real)nparts);
    if( option[in::file::s]().size() == 1 ) {
      Real frac1 = ((Real)part-1)/(Real)nparts;
      Real frac2 = ((Real)part  )/(Real)nparts;
      frac1 = 1.0 - sqrt(1.0-frac1);
      frac2 = 1.0 - sqrt(1.0-frac2);
      size1 = (Size)std::ceil(frac1*(Real)total)+1;
      size2 = (Size)std::ceil(frac2*(Real)total);
    }
  }
  TR << "SIZE " << size1 << " " << size2 << std::endl;
  return std::pair<Size,Size>(size1,size2);
}



void run() {

  using namespace basic::options::OptionKeys;

  core::chemical::ResidueTypeSetCAP  rs = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

  utility::vector1<utility::file::FileName> const & fns( option[in::file::s]() );
  for(Size ifn = 1; ifn <= fns.size(); ++ifn) {
    string fn = fns[ifn];

    Pose init;
    core::import_pose::pose_from_pdb(init,*rs,fn);

    Size nres = init.n_residue();
    // ScoreFunctionOP sf = core::scoring::get_score_function();
    ScoreFunctionOP sf = new core::scoring::symmetry::SymmetricScoreFunction(core::scoring::get_score_function());

    Pose pose = init;

    Mat R1 = rotation_matrix_degrees(Vec(0,0,1), 120.0);
    Mat R2 = rotation_matrix_degrees(Vec(0,0,1),-120.0);

    Size irsd = 0;
    for(Size i = 1; i <= pose.n_residue(); ++i) {
      if(pose.residue(i).name3()=="BPY") irsd = i;
    }

    vector1<Size> scanres;
    if(option[willmatch::residues].user()) {
      TR << "input scanres!!!!!!" << std::endl;
      scanres = option[willmatch::residues]();
    } else {
      for(Size i = 1; i <= pose.n_residue(); ++i) {
        if(!pose.residue(i).has("N" )) { continue; }
        if(!pose.residue(i).has("CA")) { continue; }
        if(!pose.residue(i).has("C" )) { continue; }
        if(!pose.residue(i).has("O" )) { continue; }
        if(!pose.residue(i).has("CB")) { continue; }
        if(pose.residue(i).name3()=="PRO") { continue; }
        scanres.push_back(i);
      }
    }

    for(vector1<Size>::const_iterator jiter = scanres.begin(); jiter != scanres.end(); ++jiter) {
      Size jrsd = *jiter;
      if(jrsd==irsd) continue;
      if( pose.residue(jrsd).name3()=="GLY" || pose.residue(jrsd).name3()=="PRO" ) continue;
      Vec CA = pose.residue(jrsd).xyz("CA");
      Vec CB = pose.residue(jrsd).xyz("CB");
      Vec SG = CB + (CA-CB).normalized()*1.8;
      SG = rotation_matrix_degrees( (CA-CB).cross(Vec(1,1,1)) ,117.5) * (SG-CB) + CB;
      Real rotang = 1.0, totrot = 0.0;
      Mat Rchi = rotation_matrix_degrees(CB-CA,rotang);
      while(totrot < 360.0) {
        // do stuff

        totrot += rotang;
        SG = Rchi * (SG-CB) + CB;

        vector1<std::pair<Vec,Vec> > axes = intersecting_disulfide_axes(CB,SG,Vec(0,0,1));
        for(Size iaxs = 1; iaxs <= axes.size(); ++iaxs) {
          Vec isct = axes[iaxs].first;
          Vec c2f1 = axes[iaxs].second;
          Real ang = angle_degrees( c2f1, isct, Vec(0,0,0) );
          if(ang > 90.0) ang -= 90.0;
          if( fabs(ang-54.7356563997) > 1.0 ) continue;
          Mat R2f = rotation_matrix_degrees(isct-c2f1,180.0);
          Vec c3f1(0,0,0);
          Vec c3f2 = R2f*(c3f1-c2f1) + c2f1;

          // clash check
          for(Size ir = 1; ir <= nres; ++ir) {
            for(Size ia = 1; ia <= /*pose.residue(ir).nheavyatoms()*/5; ++ia) {
              if(pose.residue(ir).atom_name(ia) == "ZN" && pose.residue(ir).atom_name(ia) == "NE1" && pose.residue(ir).atom_name(ia) == "NN1") continue;
              numeric::xyzVector<Real> X2 = R2f * (pose.xyz(AtomID(ia,ir))-c2f1) + c2f1;
              for(Size jr = 1; jr <= nres; ++jr) {
                for(Size ja = 1; ja <= /*pose.residue(jr).nheavyatoms()*/5; ++ja) {
                  if( X2.distance_squared(    pose.xyz(AtomID(ja,jr)) ) < 9.0 ||
                      X2.distance_squared( R1*pose.xyz(AtomID(ja,jr)) ) < 9.0 ||
                      X2.distance_squared( R2*pose.xyz(AtomID(ja,jr)) ) < 9.0) {
                    if( pose.residue(jr).atom_name(ja) != "ZN" && pose.residue(jr).atom_name(ja) != "NE1" && pose.residue(jr).atom_name(ja) != "NN1" )
                      goto cont2;
                  }
                }
              }
            }
          }
          goto done2;
        cont2:          continue;
        done2:          TR << "found nonclashing trimer/dimer: " << jrsd << " " << totrot << " " << iaxs << std::endl;
          // in my tetrahedral sym file, first 2fold should be along X and pos Z (threefold is olong Z)


          // TR << "HIT " << irsd << " " << chi1 << " " << chi2 << " " << jrsd << " " << totrot << " " << ang3f << "       " << c2f1 << "      " << isct << "       " << c3f2 << "         " << (isct-c2f1).normalized() << std::endl;
          // utility::io::ozstream out("test.pdb");
          // pose.dump_pdb(out);
          // Vec viz(0,0,0);
          // viz = viz; out<<"HETATM"<<I(5,9999)<<' '<<"ZN  "<<' '<<" ZN"<<' '<<"Z"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
          // viz = CB;  out<<"HETATM"<<I(5,9999)<<' '<<"ZN  "<<' '<<" ZN"<<' '<<"Y"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
          // viz = SG;  out<<"HETATM"<<I(5,9999)<<' '<<"ZN  "<<' '<<" ZN"<<' '<<"W"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
          // viz = isct;out<<"HETATM"<<I(5,9999)<<' '<<"ZN  "<<' '<<" ZN"<<' '<<"X"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
          // viz = c2f1;out<<"HETATM"<<I(5,9999)<<' '<<"ZN  "<<' '<<" ZN"<<' '<<"V"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
          // out.close();

          Pose symm = pose;
          core::pose::replace_pose_residue_copying_existing_coordinates(symm,jrsd,rs->name_map("CYS"));
          symm.set_xyz( AtomID(symm.residue(jrsd).atom_index("SG"),jrsd), SG );
          symm.set_xyz( AtomID(symm.residue(jrsd).atom_index("HG"),jrsd), c2f1 );
          trans_pose(symm,-isct);
          if( (c2f1-isct).z() < 0 ) {
            rot_pose(symm,Vec(1,0,0),180.0);
          }
          ang = dihedral_degrees( Vec(1,0,0), Vec(0,0,0), Vec(0,0,1), symm.residue(jrsd).xyz("HG") );
          rot_pose(symm,Vec(0,0,1),-ang);
          // TR << "HG " << symm.residue(jrsd).xyz("HG") << std::endl;

          core::pose::symmetry::make_symmetric_pose( symm );

          string fname = lead_zero_string_of(jrsd,3)+"_"+lead_zero_string_of((Size)(totrot),3)+"_"+lead_zero_string_of(iaxs,3)+"_"+utility::file_basename(fn)+".pdb";

          symm.dump_pdb(fname);

          core::pose::replace_pose_residue_copying_existing_coordinates(symm,jrsd,rs->name_map("ALA"));
          sf->score(symm);
          core::io::silent::SilentStructOP ss_out( new core::io::silent::ScoreFileSilentStruct );
          ss_out->fill_struct(symm,fname);
          sfd.write_silent_struct( *ss_out, option[ basic::options::OptionKeys::out::file::silent ]() );

        }
      }
    }
  }
}

int main (int argc, char *argv[]) {

	try {


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

  devel::init(argc,argv);
  run();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}




//
//








