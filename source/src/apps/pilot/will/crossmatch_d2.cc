// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/apps/pilat/will/crossmatch.cc
/// @brief crosses matches fast

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/smhybrid.OptionKeys.gen.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>
#include <basic/options/keys/crossmatch.OptionKeys.gen.hh>
#include <basic/options/keys/willmatch.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
#include <core/pack/optimizeH.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/packstat/compute_sasa.hh>
#include <core/scoring/packstat/sasa_dot_locations.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <numeric/NumericTraits.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/toolbox/SwitchResidueTypeSet.hh>
#include <sstream>
#include <utility/excn/Exceptions.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>

#include <apps/pilot/will/will_util.ihh>

// #include <devel/init.hh>
// #include <core/scoring/constraints/LocalCoordinateConstraint.hh>

using core::id::AtomID;
using basic::options::option;
using core::kinematics::Stub;
using core::kinematics::FoldTree;
using core::pose::Pose;
using core::pose::PoseAP;
using core::pose::PoseCAP;
using core::pose::PoseCOP;
using core::pose::PoseOP;
using core::Real;
using core::Size;
using numeric::max;
using numeric::min;
using numeric::random::gaussian;
using numeric::random::uniform;
using numeric::conversions::radians;
using numeric::constants::r::pi;
using numeric::rotation_matrix_degrees;
using numeric::xyzVector;
using numeric::xyzMatrix;
using numeric::x_rotation_matrix_degrees;
using numeric::y_rotation_matrix_degrees;
using numeric::z_rotation_matrix_degrees;
using ObjexxFCL::format::F;
using ObjexxFCL::format::I;
using ObjexxFCL::format::LJ;
using ObjexxFCL::format::SS;
using ObjexxFCL::string_of;
using std::cerr;
using std::cout;
using std::string;
using std::pair;
using utility::io::izstream;
using utility::io::ozstream;
using utility::pointer::owning_ptr;
using utility::pointer::access_ptr;
using utility::pointer::ReferenceCount;
using utility::vector1;
using std::endl;
using core::import_pose::pose_from_pdb;
typedef numeric::xyzMatrix<Real> Mat;
typedef numeric::xyzVector<Real> Vec;
typedef utility::vector1<Vec>    Vecs;


using core::pose::Pose;
using core::scoring::ScoreFunctionOP;

static thread_local basic::Tracer TR( "crossmatch" );

static core::io::silent::SilentFileData sfd;

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

inline Vec randvec() {
  Vec v(uniform(),uniform(),uniform());
  while( v.length() > 1.0 ) v = Vec(uniform(),uniform(),uniform());
  return v.normalized();
}


// in will_util now
// inline void xform_pose( core::pose::Pose & pose, Stub const & s ) {
//  for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
//    for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
//      core::id::AtomID const aid(core::id::AtomID(ia,ir));
//      pose.set_xyz( aid, s.local2global(pose.xyz(aid)) );
//    }
//  }
// }
// inline void xform_pose_rev( core::pose::Pose & pose, Stub const & s ) {
//  for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
//    for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
//      core::id::AtomID const aid(core::id::AtomID(ia,ir));
//      pose.set_xyz( aid, s.global2local(pose.xyz(aid)) );
//    }
//  }
// }

// inline void trans_pose( core::pose::Pose & pose, Vec const & trans ) {
//  for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
//    for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
//      core::id::AtomID const aid(core::id::AtomID(ia,ir));
//      pose.set_xyz( aid, pose.xyz(aid) + trans );
//    }
//  }
// }
//
// inline void rot_pose( core::pose::Pose & pose, Mat const & rot ) {
//  for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
//    for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
//      core::id::AtomID const aid(core::id::AtomID(ia,ir));
//      pose.set_xyz( aid, rot * pose.xyz(aid) );
//    }
//  }
// }
//
// inline void rot_pose( core::pose::Pose & pose, Mat const & rot, Vec const & cen ) {
//  trans_pose(pose,-cen);
//  rot_pose(pose,rot);
//  trans_pose(pose,cen);
// }
//
// inline void rot_pose( core::pose::Pose & pose, Vec const & axis, Real const & ang ) {
//  rot_pose(pose,rotation_matrix_degrees(axis,ang));
// }
//
// inline void rot_pose( core::pose::Pose & pose, Vec const & axis, Real const & ang, Vec const & cen ) {
//  rot_pose(pose,rotation_matrix_degrees(axis,ang),cen);
// }


void minimize(Pose & pose, ScoreFunctionOP sf, vector1<Size> matchres) {
  core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
  movemap->set_chi(true);
  movemap->set_bb(false);
  movemap->set_jump(false);
  for(Size i = 1; i <= matchres.size(); ++i) {
    movemap->set_chi(matchres[i],false);
  }

  if(core::pose::symmetry::is_symmetric(pose)) {
    core::pose::symmetry::make_symmetric_movemap(pose,*movemap);
    for(Size i = 1; i <= pose.fold_tree().num_jump(); ++i) {
      if((Size)pose.fold_tree().upstream_jump_residue(i)==pose.n_residue()-1 && (Size)pose.fold_tree().downstream_jump_residue(i) < pose.n_residue()/2) {
        TR << "set jump true " << i << std::endl;
        movemap->set_jump(i,true);
      }
    }
    protocols::simple_moves::symmetry::SymMinMover m( movemap, sf, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false );
    m.apply(pose);
  } else {
    protocols::simple_moves::MinMover m( movemap, sf, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false );
    m.apply(pose);
  }

}

void myoptH(Pose & pose, ScoreFunctionOP sf) {
  add_lower_terminus_type_to_pose_residue(pose,1);
  add_upper_terminus_type_to_pose_residue(pose,2);
  core::pack::optimizeH(pose,*sf);
  remove_lower_terminus_type_from_pose_residue(pose,1);
  remove_upper_terminus_type_from_pose_residue(pose,2);
}

void design(Pose & pose, ScoreFunctionOP sf, Size end_of_prot_1, vector1<Size> const & matchres, Size end_of_prot_2=0, bool designall = false ){
  using namespace core::pack::task;
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


  Real worig = sf->get_weight(core::scoring::res_type_constraint);
  if( worig == 0.0 ) sf->set_weight(core::scoring::res_type_constraint,1.0);
  utility::vector1< core::scoring::constraints::ConstraintCOP > res_cst = add_favor_native_cst(pose);

  Size nres = pose.n_residue();
  PackerTaskOP task = TaskFactory::create_packer_task(pose);
  // task->initialize_extra_rotamer_flags_from_command_line();
  vector1< bool > aas(20,true);
  aas[core::chemical::aa_cys] = false;
  aas[core::chemical::aa_his] = false;
  // aas[core::chemical::aa_met] = false;
  aas[core::chemical::aa_pro] = false;
  aas[core::chemical::aa_gly] = false;
  if(option[basic::options::OptionKeys::willmatch::exclude_ala]()) aas[core::chemical::aa_ala] = false;
  if(!designall) {
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
  }

  vector1<Size> fixed;
  if(option[basic::options::OptionKeys::willmatch::fixed_res].user()) {
    utility::io::izstream in(option[basic::options::OptionKeys::willmatch::fixed_res]());
    Size tmp;
    while(in>>tmp) fixed.push_back(tmp);
    in.close();
  }

  vector1<Size> interface;
  if( option[basic::options::OptionKeys::willmatch::design_interface]()) {
    for(Size i = 1; i <= end_of_prot_2; ++i) {
      if( sasa.size() >= i && sasa[i] > option[basic::options::OptionKeys::willmatch::cb_sasa_thresh]() ) continue;
      AtomID aid(5,i);
      if(pose.residue(i).nheavyatoms() < 5) continue; // don't design GLY or VIRT
      for(Size j = 1; j <= nres; ++j) {
        if( nres > end_of_prot_2 && sasa[i] > option[basic::options::OptionKeys::willmatch::cb_sasa_thresh]()/2.0 ) continue;
        if( i <= end_of_prot_1 && j <= end_of_prot_1 ) continue;
        if( end_of_prot_1 < i && i <= end_of_prot_2 && end_of_prot_1 < j && j <= end_of_prot_2 ) continue;
        AtomID aid2(5,j);
        if(pose.residue(j).nheavyatoms() < 5) continue; // don't design GLY or VIRT
        if(pose.xyz(aid).distance_squared(pose.xyz(aid2)) < 64) {
          interface.push_back(i);
        }
      }
    }
  }
  // TR << "INTERFACE ";
  for(Size i = 1; i <= nres; ++i) {
    if(std::find(matchres.begin(),matchres.end(),i)!=matchres.end()) {
      task->nonconst_residue_task(i).prevent_repacking();
    } else if(std::find(fixed.begin(),fixed.end(),i)!=fixed.end()){
      task->nonconst_residue_task(i).restrict_to_repacking();
    } else if(std::find(interface.begin(),interface.end(),i)!=interface.end()){
      bool tmp = aas[pose.residue(i).aa()];
      aas[pose.residue(i).aa()] = true;
      task->nonconst_residue_task(i).restrict_absent_canonical_aas(aas);
      aas[pose.residue(i).aa()] = tmp;
      task->nonconst_residue_task(i).initialize_extra_rotamer_flags_from_command_line();
    } else {
      task->nonconst_residue_task(i).restrict_to_repacking();
    }
  }
  // task->or_include_current(true);
  // TR << std::endl;
  // TR << *task << std::endl;
  // pose.dump_pdb("test.pdb");
  // if(uniform() > 0.2) std::exit(-1);
  if(core::pose::symmetry::is_symmetric(pose)) {
    protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
    repack.apply(pose);
  } else {
    protocols::simple_moves::PackRotamersMover repack( sf, task );
    repack.apply(pose);
  }

  // cleanup 2
  pose.remove_constraints( res_cst );
  sf->set_weight(core::scoring::res_type_constraint,worig);

}
void repack(Pose & pose, ScoreFunctionOP sf, Size /*end_of_prot_1*/, vector1<Size> const & matchres, Size end_of_prot_2=0 ){
  using namespace core::pack::task;
  PackerTaskOP task = TaskFactory::create_packer_task(pose);
  task->initialize_extra_rotamer_flags_from_command_line();
  if(end_of_prot_2 == 0) end_of_prot_2 = pose.n_residue();
  for(Size i = 1; i <= end_of_prot_2; ++i) {
    if(std::find(matchres.begin(),matchres.end(),i)!=matchres.end()) {
      task->nonconst_residue_task(i).prevent_repacking();
    } else {
      task->nonconst_residue_task(i).restrict_to_repacking();
    }
  }
  if(core::pose::symmetry::is_symmetric(pose)) {
    protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
    repack.apply(pose);
  } else {
    protocols::simple_moves::PackRotamersMover repack( sf, task );
    repack.apply(pose);
  }
}

void design_homodimer(Pose & pose, ScoreFunctionOP sf, vector1<Size> const & matchres, bool designall = false ){
  using namespace core::pack::task;

  Real worig = sf->get_weight(core::scoring::res_type_constraint);
  if( worig == 0.0 ) sf->set_weight(core::scoring::res_type_constraint,1.0);
  utility::vector1< core::scoring::constraints::ConstraintCOP > res_cst = add_favor_native_cst(pose);

  Size nres = pose.n_residue();
  PackerTaskOP task = TaskFactory::create_packer_task(pose);
  // task->initialize_extra_rotamer_flags_from_command_line();
  vector1< bool > aas(20,true);
  aas[core::chemical::aa_cys] = false;
  aas[core::chemical::aa_his] = false;
  // aas[core::chemical::aa_met] = false;
  aas[core::chemical::aa_pro] = false;
  aas[core::chemical::aa_gly] = false;
  if(option[basic::options::OptionKeys::willmatch::exclude_ala]()) aas[core::chemical::aa_ala] = false;
  if(!designall) {
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
  }

  vector1<Size> fixed;
  if(option[basic::options::OptionKeys::willmatch::fixed_res].user()) {
    utility::io::izstream in(option[basic::options::OptionKeys::willmatch::fixed_res]());
    Size tmp;
    while(in>>tmp) fixed.push_back(tmp);
    in.close();
  }

  core::conformation::symmetry::SymmetryInfoCOP symminfo = core::pose::symmetry::symmetry_info(pose);
  Size end_of_prot_1 = symminfo->get_nres_subunit();
  Size end_of_prot_2 = 2*end_of_prot_1;
  vector1<Size> interface;
  for(Size i = 1; i <= end_of_prot_2; ++i) {
    // if( sasa.size() >= i && sasa[i] > option[basic::options::OptionKeys::willmatch::cb_sasa_thresh]() ) continue;
    AtomID aid(5,i);
    if(pose.residue(i).nheavyatoms() < 5) continue; // don't design GLY or VIRT
    for(Size j = 1; j <= nres; ++j) {
      // if( nres > end_of_prot_2 && sasa[i] > option[basic::options::OptionKeys::willmatch::cb_sasa_thresh]()/2.0 ) continue;
      if( i <= end_of_prot_1 && j <= end_of_prot_1 ) continue;
      if( end_of_prot_1 < i && i <= end_of_prot_2 && end_of_prot_1 < j && j <= end_of_prot_2 ) continue;
      AtomID aid2(5,j);
      if(pose.residue(j).nheavyatoms() < 5) continue; // don't design GLY or VIRT
      if(pose.xyz(aid).distance_squared(pose.xyz(aid2)) < 100.0) {
        interface.push_back(i);
      }
    }
  }

  // TR << "INTERFACE ";
  for(Size i = 1; i <= nres; ++i) {
    if(std::find(matchres.begin(),matchres.end(),i)!=matchres.end()) {
      task->nonconst_residue_task(i).prevent_repacking();
    } else if(std::find(fixed.begin(),fixed.end(),i)!=fixed.end()){
      task->nonconst_residue_task(i).restrict_to_repacking();
    } else if(std::find(interface.begin(),interface.end(),i)!=interface.end()){
      bool tmp = aas[pose.residue(i).aa()];
      aas[pose.residue(i).aa()] = true;
      task->nonconst_residue_task(i).restrict_absent_canonical_aas(aas);
      aas[pose.residue(i).aa()] = tmp;
      task->nonconst_residue_task(i).initialize_extra_rotamer_flags_from_command_line();
    } else {
      task->nonconst_residue_task(i).prevent_repacking();
    }
  }
  // task->or_include_current(true);
  // TR << std::endl;
  // TR << *task << std::endl;
  // pose.dump_pdb("test.pdb");
  // if(uniform() > 0.2) std::exit(-1);
  if(core::pose::symmetry::is_symmetric(pose)) {
    protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
    repack.apply(pose);
  } else {
    protocols::simple_moves::PackRotamersMover repack( sf, task );
    repack.apply(pose);
  }

  // cleanup 2
  pose.remove_constraints( res_cst );
  sf->set_weight(core::scoring::res_type_constraint,worig);

}


struct Tags {
  string tag1;
  string tag2;
  Size rot;
  Tags(string a, string b, Size r) : tag1(a), tag2(b), rot(r) {}
};

struct Filter {
  virtual bool filter(string tag1, string tag2) const = 0;
  virtual bool filter(string tag1, string tag2, Size rot) const = 0;

};
struct NullFilter : Filter {
  bool filter(string /*tag1*/, string /*tag2*/) const { return true; }
  bool filter(string /*tag1*/, string /*tag2*/, Size /*rot*/) const { return true; }
};
struct ListFilter : Filter {
  vector1<Tags> tagslist;
  std::set<string> tags1;
  std::set<string> tags2;
  ListFilter(vector1<Tags> tagslist_in) : tagslist(tagslist_in) {
    for(Size i = 1; i<tagslist.size();++i) {
      tags1.insert(tagslist[i].tag1);
      tags2.insert(tagslist[i].tag2);
    }
  }
  bool filter(string tag1, string tag2) const {
    if( tags1.find(tag1) == tags1.end() ) return false;
    if( tags2.find(tag2) == tags2.end() ) return false;
    for(Size i = 1; i <= tagslist.size(); ++i) {
      if( tagslist[i].tag1 == tag1 && tagslist[i].tag2 == tag2) return true;
    }
    return false;
  }
  bool filter(string tag1, string tag2, Size rot) const {
    if( tags1.find(tag1) == tags1.end() ) return false;
    if( tags2.find(tag2) == tags2.end() ) return false;
    for(Size i = 1; i <= tagslist.size(); ++i) {
      if( tagslist[i].tag1 == tag1 && tagslist[i].tag2 == tag2 && tagslist[i].rot == rot) return true;
    }
    return false;
  }
};


struct MatchAlignInfo {
  Vec cen,axis,ortho;
  bool sqp;
  Size r1d,r2d;
  MatchAlignInfo() : sqp(false),r1d(0),r2d(0) {}
  // friend inline bool operator <( MatchAlignInfo const & a, MatchAlignInfo const & b ) {
  //  return ( a.cen < b.cen ) && ( a.axis < b.axis ) && ( a.ortho < b.ortho );
  // }
  friend inline bool operator ==( MatchAlignInfo const & a, MatchAlignInfo const & b ) {
    Real const TOL = basic::options::option[basic::options::OptionKeys::willmatch::identical_match_dis]();
    if( fabs(a.cen.x()  -b.cen.x()  ) > TOL ) return false;
    if( fabs(a.cen.y()  -b.cen.y()  ) > TOL ) return false;
    if( fabs(a.cen.z()  -b.cen.z()  ) > TOL ) return false;
    if( fabs(a.axis.x() -b.axis.x() ) > TOL ) return false;
    if( fabs(a.axis.y() -b.axis.y() ) > TOL ) return false;
    if( fabs(a.axis.z() -b.axis.z() ) > TOL ) return false;
    if( fabs(a.ortho.x()-b.ortho.x()) > TOL ) return false;
    if( fabs(a.ortho.y()-b.ortho.y()) > TOL ) return false;
    if( fabs(a.ortho.z()-b.ortho.z()) > TOL ) return false;
    return true;
  }
  void dump_pdb(std::ostream & out) {
    Vec axs = axis+cen;
    Vec ort = ortho+cen;
    out << "HETATM" << I(5,9994) << ' ' << " CEN" << ' ' << "MAI" << ' ' << "A" << I(4,994) << "    " << F(8,3,cen.x()) << F(8,3,cen.y()) << F(8,3,cen.z()) << F(6,2,1.0) << F(6,2,1.0) << std::endl;
    out << "HETATM" << I(5,9995) << ' ' << "AXIS" << ' ' << "MAI" << ' ' << "A" << I(4,995) << "    " << F(8,3,axs.x()) << F(8,3,axs.y()) << F(8,3,axs.z()) << F(6,2,1.0) << F(6,2,1.0) << std::endl;
    out << "HETATM" << I(5,9996) << ' ' << "ORTH" << ' ' << "MAI" << ' ' << "A" << I(4,996) << "    " << F(8,3,ort.x()) << F(8,3,ort.y()) << F(8,3,ort.z()) << F(6,2,1.0) << F(6,2,1.0) << std::endl;
  }
  void dump_pdb(std::string fname) {
    ozstream out(fname);
    dump_pdb(out);
    out.close();
  }
};
std::ostream & operator<< (std::ostream & out, MatchAlignInfo & ai) {
  out << ai.cen << " " << ai.axis << " " << ai.ortho << std::endl;
  return out;

}

// struct ltmai
// {
//   bool operator()(const MatchAlignInfo & a, const MatchAlignInfo & b) const
//   {
//    return ( a.cen < b.cen ) && ( a.axis < b.axis ) && ( a.ortho < b.ortho );
//   }
// };

void setup_favor_native(core::pose::Pose & pose, core::pose::Pose const & native) {
  using namespace basic::options;
  using namespace core::scoring::constraints;
  core::Real bonus = option[OptionKeys::enzdes::favor_native_res].value();
  TR.Info << "favor_native_res: adding a bonus of " << bonus << " for native residues to pose." << std::endl;
  // should remove exising residue type constraints???
  utility::vector1< core::scoring::constraints::ConstraintCOP > favor_native_constraints;
  for( core::Size i = 1; i <= pose.total_residue(); ++i){
    // if( task->design_residue(i) ){
    ConstraintOP resconstraint = new ResidueTypeConstraint( native, i, bonus );
    favor_native_constraints.push_back( resconstraint );
    // }
  }
  pose.add_constraints( favor_native_constraints );
}

struct MatchAligner : public utility::pointer::ReferenceCount {
  virtual MatchAlignInfo align_info(Pose const & lig ) = 0;
  virtual vector1<Stub> align_rot( MatchAlignInfo const & mi1, MatchAlignInfo const & mi2 ) = 0;
  virtual bool checkalign( Pose const & a, Pose const & b, numeric::xyzVector<Real> & c ) = 0;
};
typedef utility::pointer::owning_ptr<MatchAligner> MatchAlignerOP;

struct Tet4HMatchAligner : public MatchAligner {
  Size bondedN(Pose const & lig, Size i) {
    Vec m = lig.residue(3).xyz(1);
    Size nd1 = lig.residue(i).atom_index("ND1");
    Size ne2 = lig.residue(i).atom_index("NE2");
    if( lig.xyz(AtomID(nd1,i)).distance(m) < lig.xyz(AtomID(ne2,i)).distance(m) ) {
      return nd1;
    } else {
      return ne2;
    }
  }
  MatchAlignInfo align_info(Pose const & lig ) {
    const core::Real PI = numeric::NumericTraits<Real>::pi();
    MatchAlignInfo mai;
    Size in1 = bondedN(lig,1);
    Size in2 = bondedN(lig,2);
    mai.r1d = in1;
    mai.r2d = in2;
    Vec  n1 = lig.residue(1).xyz(in1);
    Vec  n2 = lig.residue(2).xyz(in2);
    // Vec  h1,h2;
    // if( in1 == lig.residue(1).atom_index("NE2") ) h1 = lig.residue(1).xyz("HE2");
    // else                                          h1 = lig.residue(1).xyz("HD1");
    // if( in2 == lig.residue(2).atom_index("NE2") ) h2 = lig.residue(2).xyz("HE2");
    // else                                          h2 = lig.residue(2).xyz("HD1");
    // Vec cd1 = lig.residue(1).xyz("CD2");
    // Vec ce1 = lig.residue(1).xyz("CE1");
    // Vec cd2 = lig.residue(2).xyz("CD2");
    // Vec ce2 = lig.residue(2).xyz("CE1");
    Vec   m = lig.residue(3).xyz(1);
    // Vec m1 = n1 + 2.0*( n1 - (cd1+ce1)/2.0 ).normalized();
    // Vec m2 = n2 + 2.0*( n2 - (cd2+ce2)/2.0 ).normalized();
    // if( numeric::conversions::degrees(numeric::angle_radians(n1,m,n2)) < 100.0 ) {
    //  mai.sqp = true;
    // }
    // mai.cen   = ( m1+m2 )/ 2.0;
    mai.cen = m;
    mai.axis = (mai.cen - (n1+n2)/2.0).normalized();
    mai.ortho = (n1-n2).cross(mai.axis).normalized();
    Real dany = basic::options::option[basic::options::OptionKeys::willmatch::max_dis_any]();
    Real dall = basic::options::option[basic::options::OptionKeys::willmatch::max_dis_all]();
    // Real dmet = basic::options::option[basic::options::OptionKeys::willmatch::max_dis_metal]();
    if( (n1.distance(mai.cen) > dall && n2.distance(mai.cen) > dall) ||
        n1.distance(mai.cen) > dany ||
        n2.distance(mai.cen) > dany /*||*/
        /*m1.distance(m2) > dmet*/ )
      // if( false )
      {
        mai.cen = 0;
        mai.axis = 0;
        mai.ortho = 0;
        // TR << "align_info FAILED FOR MATCH!!!!" << std::endl;
      }
    assert(fabs(mai.axis.dot(mai.ortho))<0.001);
    return mai;
  }
  vector1<Stub> align_rot( MatchAlignInfo const & mi1, MatchAlignInfo const & mi2 ) {
    vector1<Stub> rots;
    if(mi1.sqp != mi2.sqp) return rots;
    Vec axis = mi2.axis.cross(mi1.axis);
    while( axis.length() < 0.001 ) {
      // axis = mi2.axis.cross(Vec(uniform(),uniform(),uniform()));
      axis = mi2.axis.cross(Vec(1,1,1).normalized());
    }
    axis = axis.normalized();
    Real ang1 = 180.0-acos(max(-1.0,min(1.0,mi2.axis .dot(                                   mi1.axis ))))*180.0/PI;
    Real ang2 =      -acos(max(-1.0,min(1.0,mi2.ortho.dot(rotation_matrix_degrees(axis,ang1)*mi1.ortho))))*180.0/PI;
    Vec axis2 = mi2.ortho.cross(rotation_matrix_degrees(axis,ang1)*mi1.ortho);
    // TR << "ANG2 " << ang2 << std::endl;
    // TR << "AXIS1 " <<     axis << " " << ang1 << std::endl;
    // TR << "AXIS2 " << mi2.axis << " " << ang2 << std::endl;
    // ang2 = 0;
    Mat rot1;
    Mat rot2;
    if(mi1.sqp) {
      // TR << "HHGEOM SQP" << std::endl;
      rot1 = rotation_matrix_degrees(axis2,ang2+  0.0) * rotation_matrix_degrees(axis,ang1);
      rot2 = rotation_matrix_degrees(axis2,ang2+180.0) * rotation_matrix_degrees(axis,ang1);Vec trans1 = mi2.cen - rot1*mi1.cen;
    } else {
      // TR << "HHGEOM TET" << std::endl;
      rot1 = rotation_matrix_degrees(axis2,ang2+ 90.0) * rotation_matrix_degrees(axis,ang1);
      rot2 = rotation_matrix_degrees(axis2,ang2- 90.0) * rotation_matrix_degrees(axis,ang1);Vec trans1 = mi2.cen - rot1*mi1.cen;
    }
    Vec trans1 = mi2.cen - rot1*mi1.cen;
    Vec trans2 = mi2.cen - rot2*mi1.cen;
    rots.push_back( Stub(rot1,trans1) );
    rots.push_back( Stub(rot2,trans2) );

    return rots;
  }
  bool checkalign( Pose const & a, Pose const & b, numeric::xyzVector<Real> & c ) {
    /////!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Size  na1 = bondedN(a,1);
    Size  na2 = bondedN(a,2);
    Size  na3 = bondedN(b,1);
    Size  na4 = bondedN(b,2);
    numeric::xyzVector<Real> n1 = a.xyz(AtomID(na1,1));
    numeric::xyzVector<Real> n2 = a.xyz(AtomID(na2,2));
    numeric::xyzVector<Real> n3 = b.xyz(AtomID(na3,1));
    numeric::xyzVector<Real> n4 = b.xyz(AtomID(na4,2));
    c = (n1+n2+n3+n4)/4.0;
    Real zndis1 = c.distance(n1);
    Real zndis2 = c.distance(n2);
    Real zndis3 = c.distance(n3);
    Real zndis4 = c.distance(n4);
    Real znang12 = angle_degrees(n1,c,n2);
    Real znang13 = angle_degrees(n1,c,n3);
    Real znang14 = angle_degrees(n1,c,n4);
    Real znang23 = angle_degrees(n2,c,n3);
    Real znang24 = angle_degrees(n2,c,n4);
    Real znang34 = angle_degrees(n3,c,n4);
    TR << "ZNSITE GEOM:" << std::endl;
    TR << "zndis1 " << zndis1 << std::endl;
    TR << "zndis3 " << zndis3 << std::endl;
    TR << "zndis2 " << zndis2 << std::endl;
    TR << "zndis4 " << zndis4 << std::endl;
    TR << "znang12 " << znang12 << std::endl;
    TR << "znang34 " << znang34 << std::endl;
    TR << "znang13 " << znang13 << std::endl;
    TR << "znang24 " << znang24 << std::endl;
    TR << "znang14 " << znang14 << std::endl;
    TR << "znang23 " << znang23 << std::endl;
    if( zndis1 > 2.4 || zndis2 > 2.4 || zndis3 > 2.4 || zndis4 > 2.4 ) return false;
    if( zndis1 < 1.8 || zndis1 < 1.8 || zndis1 < 1.8 || zndis1 < 1.8 ) return false;
    Real MXA = 120.0 , MNA= 97.0;
    if( znang12 > MXA || znang13 > MXA || znang14 > MXA || znang23 > MXA || znang24 > MXA || znang34 > MXA ) return false;
    if( znang12 < MNA || znang13 < MNA || znang14 < MNA || znang23 < MNA || znang24 < MNA || znang34 < MNA ) return false;
    Real d = dihedral_degrees(a.xyz(AtomID(na1,1)),a.xyz(AtomID(na2,2)),b.xyz(AtomID(na3,1)),b.xyz(AtomID(na4,2)));
    return true;
    // TR << "D " << d << std::endl;
    // if(sqp) {
    if( fabs(d) < 5.0 ) return true;
    if( fabs(fabs(d)-180.0) < 10.0 ) return true;
    // } else {
    if( fabs(fabs(d)- 75.0) < 10.0 ) return true;
    // if( fabs(fabs(d)-270.0) < 5.0 ) return true;
    // }
    return false;
  }

};

struct MatchBase : public utility::pointer::ReferenceCount {
  Pose pose;
  // utility::vector1< core::id::AtomID > bb_bur_uns;
  MatchBase() {}
  void init(Pose & ipose) {
    // bb_bur_uns = count_bb_bur_uns(ipose);
    Pose tpose = ipose;
    core::chemical::ResidueTypeSetCAP cen_residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID );
    core::chemical::ResidueType const & ala( cen_residue_set->name_map("ALA") );
    vector1<Size> positions;
    Size start = 1,i;
    while( start >= tpose.n_residue() && !tpose.residue(start).is_protein() ) start++;
    for( i = start; i <= tpose.n_residue(); ++i) {
      if(!tpose.residue(i).is_protein()) break;
      core::pose::replace_pose_residue_copying_existing_coordinates(tpose,i,ala);
      positions.push_back(i);
    }
    core::kinematics::FoldTree f(positions.size());
    core::pose::create_subpose(tpose,positions,f,pose);
    core::pose::add_lower_terminus_type_to_pose_residue(pose,1);
    core::pose::add_upper_terminus_type_to_pose_residue(pose,pose.n_residue());

  }
};
typedef utility::pointer::owning_ptr<MatchBase> MatchBaseOP;
typedef utility::pointer::owning_ptr<const MatchBase> MatchBaseCOP;

struct MatchLig : public utility::pointer::ReferenceCount {
  Pose pose;
  Size rsd1,rsd2;
  string aa1,aa2,tag;
  MatchAlignInfo align_info;
  MatchLig(Pose & ipose, Size irsd1, Size irsd2, string iaa1, string iaa2, string itag, MatchAlignerOP mop)
    : rsd1(irsd1), rsd2(irsd2), aa1(iaa1), aa2(iaa2), tag(itag)
  {
    if(ipose.n_residue()>3) {
      assert(aa1[0]==ipose.residue(rsd1).name1());
      assert(aa2[0]==ipose.residue(rsd2).name1());
      core::conformation::Residue r1 = ipose.residue(rsd1);
      core::conformation::Residue r2 = ipose.residue(rsd2);
      pose.append_residue_by_jump(r1,1,"","",true);
      pose.append_residue_by_jump(r2,1,"","",true);

      Size i = 1;
      while(i <= ipose.n_residue() && ipose.residue(i).is_protein()) ++i;
      pose.append_residue_by_jump(ipose.residue(i),1,"","",true);
      // for(Size j = 1; j <= pose.residue(3).natoms(); ++j) {
      //  pose.set_xyz(AtomID(j,3),ipose.residue(i).xyz(j));
      // }
    } else {
      pose = ipose;
    }
    for(Size i = 1; i <= 2; ++i) {
      if(pose.residue(i).is_lower_terminus()) core::pose::remove_lower_terminus_type_from_pose_residue(pose,i);
      if(pose.residue(i).is_upper_terminus()) core::pose::remove_upper_terminus_type_from_pose_residue(pose,i);
    }
    align_info = mop->align_info(pose);

  }
  bool isbad() { return align_info.axis.length() < 0.01; }
};

struct MatchSet {
  MatchBase base;
  Pose native;
  vector1<MatchLig> ligs;
  MatchAlignerOP mop;
  utility::vector1<Vec> points_;
  ObjexxFCL::FArray3D< vector1<Vec> > cubes_;
  Vec bbl_;
  numeric::xyzTriple< platform::Size > cube_dim_;
  Real side_inv_, neighbor_cutoff_sq_;
  vector1<MatchAlignInfo> xforms_seen_;
  vector1<Size> iface_candidates;
  MatchSet() {}
  MatchSet(vector1<Pose> & poses, MatchAlignerOP imop, Pose & native_in, utility::vector1<Size> allowed_res = utility::vector1<Size>() ) {
    init(poses,imop,native_in,allowed_res);
  }
  void init(vector1<Pose> & poses, MatchAlignerOP imop, Pose & native_in, utility::vector1<Size> allowed_res) {
    native = native_in;
    mop = imop;
    base.init(poses[1]);
    // for(Size i = 1; i <= base.pose.n_residue(); ++i) {
    //  if(base.pose.residue(i).is_lower_terminus()) core::pose::remove_lower_terminus_type_from_pose_residue(base.pose,i);
    //  if(base.pose.residue(i).is_upper_terminus()) core::pose::remove_upper_terminus_type_from_pose_residue(base.pose,i);
    // }
    ScoreFunctionOP sflig = core::scoring::get_score_function();
    sflig->set_weight(core::scoring::fa_dun,0.01);
    sflig->set_weight(core::scoring::fa_intra_rep,0.44);
    sflig->set_weight(core::scoring::ref,0.0);
    Size numbad = 0, numskip=0, scorefail=0, num_disallowed=0;
    for(Size idx = 2; idx <= poses.size(); ++idx) {
      // TR << "read lig " << idx << std::endl;
      {
        Pose tmp;
        tmp.append_residue_by_jump(poses[idx].residue(1),1);
        tmp.append_residue_by_jump(poses[idx].residue(2),1);
        myoptH(tmp,sflig);
        if(tmp.energies().total_energy() > basic::options::option[basic::options::OptionKeys::crossmatch::his_score_cut]() ) {
          // tmp.energies().show();
          // tmp.dump_pdb("his_his_score_fail_"+string_of(tmp.energies().total_energy())+".pdb");
          scorefail++;
          continue;
        }
      }
      // sflig->show(tmp);
      // TR << "CROSSMATCHLIG lig" << idx << ".pdb " << tmp.energies().total_energy() << std::endl;
      // tmp.dump_pdb("lig"+string_of(idx)+".pdb");

      std::string s = poses[idx].pdb_info()->modeltag();
      // TR << "'" << s << "'" << std::endl;
      s = utility::file_basename(s);
      if(s == "BASE.pdb") {
        TR << "unexpeted BASE.pdb!" << std::endl;
        std::exit(-1);
      }
      string aa1,aa2,ri1,ri2,geom="";
      if(s.substr(0,7)=="HIS HIS") {
        aa1="H"; aa2="H";
        std::istringstream iss(s);
        string tmp;
        iss >> tmp >> tmp >> tmp >> geom >> ri1 >> tmp >> tmp >> ri2;
        // TR << "S " << ri1 << " " << ri2 << std::endl;
        if( poses[idx].residue(1).name3()!="HIS" ||
            poses[idx].residue(2).name3()!="HIS" ||
            poses[idx].residue(3).name3()!=" ZN" ) {
          TR << "LIG out of order! " << s << std::endl;
          TR << "'" << poses[idx].residue(1).name() << "'" << std::endl;
          TR << "'" << poses[idx].residue(2).name() << "'" << std::endl;
          TR << "'" << poses[idx].residue(3).name() << "'" << std::endl;
          std::exit(-1);
        }
      } else {
        s = s.substr(3);
        s = s.substr( s.find("_")+1 );
        s = s.substr(0, s.find("_"));
        Size i = 0;
        while( i < s.size() && 64 < (int)s[i] && (int)s[i] < 91 ) aa1 += s[i++]; // A-Z
        while( i < s.size() && 47 < (int)s[i] && (int)s[i] < 58 ) ri1 += s[i++]; // 0-9
        while( i < s.size() && 64 < (int)s[i] && (int)s[i] < 91 ) aa2 += s[i++];
        while( i < s.size() && 47 < (int)s[i] && (int)s[i] < 58 ) ri2 += s[i++];
      }
      // TR << "s " << poses[idx].pdb_info()->modeltag() << " " << s << " " << ri1 << " " << ri2 << std::endl;
      Size rsd1,rsd2;
      rsd1 = atoi(ri1.c_str());
      rsd2 = atoi(ri2.c_str());

      if( allowed_res.size() > 0 ) { // skip if not in allowed res
        if( std::find(allowed_res.begin(),allowed_res.end(),rsd1) == allowed_res.end() ||
            std::find(allowed_res.begin(),allowed_res.end(),rsd2) == allowed_res.end() ) {
          num_disallowed++;
          continue;
        }
      }
      // why skip ends?
      // if(    rsd1 == 1 || rsd2 == 1 || rsd1 == base.pose.n_residue() || rsd2 == base.pose.n_residue()) {
      //  TR << "skipping end " << rsd1 << " " << rsd2 << std::endl;
      //  continue;
      // }
      // assert(rsd1  > 1 && rsd2  > 1 && rsd1 <= base.pose.n_residue() && rsd2 <= base.pose.n_residue());

      // TR << "make MatchLig: " << s << std::endl;
      MatchLig ml = MatchLig(poses[idx],rsd1,rsd2,aa1,aa2,utility::file_basename(poses[idx].pdb_info()->modeltag()),mop);
      if(geom!="") {
        if     (geom=="tet") ml.align_info.sqp = false;
        else if(geom=="sqp") ml.align_info.sqp = false; // ALWAYS TET FOR NOW!!!
        else {
          utility_exit_with_message("unknown geom '"+geom+"'");
        }
      }
      if(!ml.isbad()) {
        if( std::find(xforms_seen_.begin(),xforms_seen_.end(),ml.align_info) == xforms_seen_.end() ) {
          // TR << "good AI " << ml.align_info << std::endl;
          xforms_seen_.push_back(ml.align_info);
          ligs.push_back(ml);
        } else {
          // TR << "skipping match with already-seen align-info" << std::endl;
          numskip++;
        }
      } else {
        numbad++;
        // TR << "lig bad" << std::endl;
      }
    }
    TR << "MatchSet total " << poses.size()-1 << " " << numskip << " redundant, " << numbad << " badgeom, " << scorefail << " badscore, " << num_disallowed << " disallowed residues " << ligs.size() << " good!" << std::endl;
    init_clash_check( basic::options::option[basic::options::OptionKeys::crossmatch::clash_dis]() );
  }
  void write_to_file(string fname) {
    bool ov = basic::options::option[basic::options::OptionKeys::out::file::output_virtual]();
    basic::options::option[basic::options::OptionKeys::out::file::output_virtual](true);
    utility::io::ozstream out(fname);
    TR << "writing " << fname << " BASE" << std::endl;
    out << "MODEL BASE" << endl;
    core::io::pdb::dump_pdb(base.pose,out);
    out << "ENDMDL" << endl;
    for(Size i = 1; i <= ligs.size(); ++i) {
      TR << "writing " << fname << " " << ligs[i].tag << std::endl;
      out << "MODEL " << ligs[i].tag << endl;
      Pose tmp = ligs[i].pose;
      if(!tmp.residue(1).is_lower_terminus()) core::pose::add_lower_terminus_type_to_pose_residue(tmp,1);
      if(!tmp.residue(1).is_upper_terminus()) core::pose::add_upper_terminus_type_to_pose_residue(tmp,1);
      if(!tmp.residue(1).is_lower_terminus()) core::pose::add_lower_terminus_type_to_pose_residue(tmp,2);
      if(!tmp.residue(1).is_upper_terminus()) core::pose::add_upper_terminus_type_to_pose_residue(tmp,2);
      core::io::pdb::dump_pdb(tmp,out);
      out << "ENDMDL" << endl;
    }
    out.close();
    basic::options::option[basic::options::OptionKeys::out::file::output_virtual](ov);
  }
  void init_clash_check(Real neighbor_cutoff) {
    using numeric::min;
    using numeric::max;
    using numeric::square;
    typedef  numeric::xyzTriple< platform::Size >  CubeDim; // Cube dimensions
    typedef  numeric::xyzTriple< platform::Size >  CubeKey; // Cube index-triple key

    neighbor_cutoff_sq_ = ( neighbor_cutoff*neighbor_cutoff);
    points_.resize(base.pose.n_residue()*5);
    for(Size i = 0; i < base.pose.n_residue(); ++i) {
      for(Size j = 1; j <= 5; ++j) points_[5*i+j] = base.pose.xyz(AtomID(j,i+1));
    }

    Vec bbu( points_[ 1 ] ); // Lower and upper corners of bounding box
    bbl_ = bbu;
    for ( Size ii = 2; ii <= points_.size(); ++ii ) { bbl_.min( points_[ ii ] ); bbu.max( points_[ ii ] ); }
    bbl_ -= 10 * std::numeric_limits< Real >::epsilon();
    bbu += 10 * std::numeric_limits< Real >::epsilon();
    Real const side( neighbor_cutoff );
    assert( side > platform::Real( 0 ) );
    side_inv_ = ( platform::Real( 1 ) / side );
    cube_dim_ = CubeDim( // Cube dimensions
                        platform::Size( std::ceil( ( bbu.x() - bbl_.x() ) * side_inv_ ) ),             // Test that ceil values == platform::Size values
                        platform::Size( std::ceil( ( bbu.y() - bbl_.y() ) * side_inv_ ) ),
                        platform::Size( std::ceil( ( bbu.z() - bbl_.z() ) * side_inv_ ) )
                         );

    cubes_.dimension( cube_dim_.x(), cube_dim_.y(), cube_dim_.z() );

    for ( Size i = 1; i <= points_.size(); ++i ) {
      Vec const pp( points_[ i ]);
      CubeKey const cube_key(
                             platform::Size( ( pp.x() - bbl_.x() ) * side_inv_ ) + 1,
                             platform::Size( ( pp.y() - bbl_.y() ) * side_inv_ ) + 1,
                             platform::Size( ( pp.z() - bbl_.z() ) * side_inv_ ) + 1
                             );
      assert( cube_key.x() <= cube_dim_.x() );
      assert( cube_key.y() <= cube_dim_.y() );
      assert( cube_key.z() <= cube_dim_.z() );
      Size i_index = cubes_.index( cube_key.x(), cube_key.y(), cube_key.z() );
      cubes_[ i_index ].push_back( pp );
    }
  }
  inline bool clash_check(Vec const & pp) const {
    Size const icx( Size( ( pp.x() - bbl_.x() ) * side_inv_ ) + 1 );
    Size const icy( Size( ( pp.y() - bbl_.y() ) * side_inv_ ) + 1 );
    Size const icz( Size( ( pp.z() - bbl_.z() ) * side_inv_ ) + 1 );
    for ( Size ix = max( icx, Size( 2 ) ) - 1,  ixe = min( icx + 1, cube_dim_.x() ); ix <= ixe; ++ix ) {
      for ( Size iy = max( icy, Size( 2 ) ) - 1, iye = min( icy + 1, cube_dim_.y() ); iy <= iye; ++iy ) {
        for ( Size iz = max( icz, Size( 2 ) ) - 1, ize = min( icz + 1, cube_dim_.z() ); iz <= ize; ++iz ) {
          Size cube_index = cubes_.index( ix, iy, iz );
          if ( cubes_[ cube_index ].size() != 0 ) { // Cube exists
            for ( vector1<Vec>::const_iterator ia = cubes_[ cube_index ].begin(), iae = cubes_[ cube_index ].end(); ia != iae; ++ia ) {
              Vec const j( *ia );
              Real const d_sq( distance_squared( pp, j ) );
              if ( d_sq <= neighbor_cutoff_sq_ ) {
                return false;
              }
            }
          }
        }
      }
    }
    return true;
  }
  inline bool clash_check(Stub const & stub) const {
    for(vector1<Vec>::const_iterator i=points_.begin(),ie=points_.end(); i!=ie; ++i) {
      if(!clash_check( stub.local2global(*i))) return false;
    }
    return true;
  }
  inline bool clash_check(Vec const & pp, Stub const & stub) const {
    return clash_check( stub.local2global(pp) );
  }
  inline bool clash_check_inv(Vec const & pp, Stub const & stub) const {
    return clash_check( stub.global2local(pp) );
  }
  inline bool clash_check(vector1<Vec> const & pps) const {
    for ( vector1<Vec>::const_iterator ia = pps.begin(), iae = pps.end(); ia != iae; ++ia ) {
      if(!clash_check(*ia)) return false;
    }
    return true;
  }
  inline bool clash_check(Pose const & pose) const {
    for(Size i = 1; i <= pose.n_residue(); ++i) {
      for(Size j = 1; j <= 5; ++j) {
        if(!clash_check(pose.xyz(AtomID(j,i)))) return false;
      }
    }
    return true;
  }
  inline bool clash_check(Pose const & pose, Stub const & stub) const {
    for(Size i = 1; i <= pose.n_residue(); ++i) {
      for(Size j = 1; j <= 5; ++j) {
        if(!clash_check(stub.local2global(pose.xyz(AtomID(j,i))))) return false;
      }
    }
    return true;
  }
  inline bool c2_clash_check(Stub const & hdstub, MatchSet const & ms2, Pose const & pose, Stub const & c2stub, vector1<Vec> const & extra = vector1<Vec>() ) const {
    for(Size i = 1; i <= pose.n_residue(); ++i) for(Size j = 1; j <= 5; ++j) if(!clash_check(c2stub.M*pose.xyz(AtomID(j,i))+c2stub.v)) return false;
    for(Size i = 1; i <= extra.size(); ++i) if(!clash_check(c2stub.M*extra[i]+c2stub.v)) return false;
    Mat const rot   = hdstub.M.transposed() * c2stub.M;
    Vec const trans = hdstub.M.transposed() * (c2stub.v-hdstub.v) ;
    for(Size i = 1; i <= pose.n_residue(); ++i) for(Size j = 1; j <= 5; ++j) if(!ms2.clash_check(rot*pose.xyz(AtomID(j,i))+trans)) return false;
    for(Size i = 1; i <= extra.size(); ++i) if(!clash_check(rot*extra[i]+trans)) return false;
    return true;
  }
  inline bool clash_check_naive(Pose & pose) const {
    for(Size i = 1; i <= base.pose.n_residue(); ++i) {
      for(Size j = 1; j <= 5; ++j) {
        Vec const xyz1( base.pose.xyz(AtomID(j,i)) );
        for(Size i2 = 1; i2 <= pose.n_residue(); ++i2) {
          for(Size j2 = 1; j2 <= 5; ++j2) {
            Vec const xyz2( pose.xyz(AtomID(j2,i2)) );
            Real const d_sq( distance_squared( xyz1, xyz2 ) );
            if ( d_sq <= neighbor_cutoff_sq_ ) {
              return false;
            }
          }
        }
      }
    }
    return true;
  }
  inline Real sqr(Real const & r) const { return r*r; }
  inline Real sigmoidish_neighbor( Real const & sqdist ) const
  {
    if( sqdist > 9.*9. ) {
      return 0.0;
    } else if( sqdist < 6.*6. ) {
      return 1.0;
    } else {
      Real dist = sqrt( sqdist );
      return sqr(1.0  - sqr( (dist - 6.) / (9. - 6.) ) );
    }
  }
  inline Real iface_check(Pose const & pose, vector1<Size> const & include_other) const {
    Real num = 0;
    for(vector1<Size>::const_iterator i=include_other.begin(),ie=include_other.end(); i != ie; ++i) {
      for(vector1<Size>::const_iterator j=iface_candidates.begin(),je=iface_candidates.end(); j != je; ++j) {
        num += sigmoidish_neighbor( pose.residue(*i).xyz(5).distance_squared(base.pose.residue(*j).xyz(5)) );
      }
    }
    return num;
  }
  inline Real c2_iface_check(Pose const & pose, Stub const & stub, vector1<Size> const & include) const {
    Real num = 0;
    for(vector1<Size>::const_iterator i=include.begin(),ie=include.end(); i != ie; ++i) {
      for(vector1<Size>::const_iterator j=include.begin(),je=include.end(); j != je; ++j) {
        num += sigmoidish_neighbor(pose.residue(*i).xyz(5).distance_squared( (stub.M*(pose.residue(*j).xyz(5)))+stub.v ));
      }
    }
    return num;
  }
  inline Size c2_linker_check_dist(Pose const & outpose, MatchSet const & b, core::kinematics::Stub const & s) {
    Vec ct1 = outpose.residue(b.base.pose.n_residue()  ).xyz("C");
    Vec nt1 = outpose.residue(                        1).xyz("N");
    Vec ct2 = outpose.residue(    outpose.n_residue()  ).xyz("C");
    Vec nt2 = outpose.residue(b.base.pose.n_residue()+1).xyz("N");
    ct2 = s.M*ct2+s.v;
    nt2 = s.M*nt2+s.v;
    // if( ct1.distance_squared(nt2) < thresh*thresh || ct2.distance_squared(nt1) < thresh*thresh ) {
    //  TR << "linker pass " << thresh << std::endl;
    // }
    return numeric::min( ct1.distance_squared(nt2), ct2.distance_squared(nt1) );

  }
  inline Vec com(Pose const & pose) const {
    Vec com(0,0,0);
    for(Size i = 1; i <= pose.n_residue(); ++i) com += pose.xyz(AtomID(2,i));
    return com / pose.n_residue();
  }
  inline bool get_contacting_stub(Pose const & outpose, Vec & trans, Stub & s, Stub const & hdstub,
                                  MatchSet const & b,vector1<Vec> const & his_atoms, Real const & LINK_CST,
                                  Real step = 16.0) {
    Vec com1 = com(outpose);
    s.v = com1 - s.M*com1;
    bool linkerfail = false;
    while(fabs(step) > 0.2) {
      s.v += trans*step;
      Size count = 0;
      while( step * (2.0*c2_clash_check(hdstub,b,outpose,s,his_atoms)-1.0) < 0 ) {
        s.v += trans*step;
        if(++count > 10) utility_exit_with_message("C2 contact search fail");
      }
      if( c2_linker_check_dist(outpose,b,s) > (LINK_CST+fabs(step))*(LINK_CST+fabs(step)) ) linkerfail = true;
      if(linkerfail) break;
      step *= -0.5;
    }
    if(linkerfail) return false;
    return true;
  }
  void selfcross_c2(Size idx1, Size itrans) {
    vector1<Stub> rots = mop->align_rot(ligs[idx1].align_info,ligs[idx1].align_info);
    if(rots.size() < itrans) return;
    if(!clash_check(base.pose,rots[itrans])) return;

    utility_exit_with_message("!!!!!!!!!!!!!!!!!!!!!!!!!");
    // Pose outpose = base.pose;
    // rot_pose  (outpose,rots[itrans].M);
    // trans_pose(outpose,rots[itrans].v);
    // Pose move     = b.ligs[idx2].pose;
    // rot_pose  (move,rots[itrans].M);
    // trans_pose(move,rots[itrans].v);
    // // if(!clash_check(move)) return;
    //
    // vector1<Size> include_res = b.iface_candidates;
    // for(vector1<Size>::const_iterator i=iface_candidates.begin(),ie=iface_candidates.end(); i != ie; ++i) {
    //  include_res.push_back(*i+b.base.pose.n_residue());
    // }
    // Real iface = iface_check(outpose,b.iface_candidates);
    // if(iface < (Real)basic::options::option[basic::options::OptionKeys::willmatch::interface_size]()) return;


  }
  void cross_homodimer(Size idx1,Size itrans) {
    // TODO how do I build a new residue on a terminal spot? get upper/lower connect error
    // need lower/upper patch on types?
    {
      Size r1 = ligs[idx1].rsd1;
      Size r2 = ligs[idx1].rsd2;
      if( r1 == 1 ) return;
      if( r2 == 1 ) return;
      if( r1 == base.pose.n_residue() ) return;
      if( r2 == base.pose.n_residue() ) return;
      if( base.pose.chain(r1)!=base.pose.chain(r1-1) ) return;
      if( base.pose.chain(r1)!=base.pose.chain(r1+1) ) return;
      if( base.pose.chain(r2)!=base.pose.chain(r2-1) ) return;
      if( base.pose.chain(r2)!=base.pose.chain(r2+1) ) return;
    }

    vector1<Stub> rots = mop->align_rot(ligs[idx1].align_info,ligs[idx1].align_info);
    if(rots.size() < itrans) return;

    // main clash check
    if(!clash_check(base.pose,rots[itrans])) return;

    Pose fixd = ligs[idx1].pose;
    Pose move = ligs[idx1].pose;
    rot_pose  (move,rots[itrans].M);
    trans_pose(move,rots[itrans].v);

    // clash check base (self) vs move (other HISs)
    for(Size ia = 1; ia <= move.residue(1).nheavyatoms(); ++ia) {
      if( !this->clash_check(move.residue(1).xyz(ia)) ||
          !this->clash_check(move.residue(2).xyz(ia)) ) {
        // TR << "clash_check on move HIS HIS fails!" << std::endl;
        // base.pose.dump_pdb("base.pdb");
        // move.dump_pdb("move.pdb");
        // utility_exit_with_message("lasdfgj");
        return;
      }
    }

    // check his vs his
    for(Size ir = 1; ir <= 2; ++ir) {
      for(Size ia = 1; ia <= ligs[idx1].pose.residue(ir).nheavyatoms(); ++ia) {
        if( 7==ia || 10==ia ) continue; //sheffler RISK assuming ND NE are 7 and 10!!!!
        Vec const xyz1 = ligs[idx1].pose.residue(ir).xyz(ia);
        for(Size jr = 1; jr <= 2; ++jr) {
          for(Size ja = 1; ja <= move.residue(jr).nheavyatoms(); ++ja) {
            if( 7==ja || 10==ja ) continue; //sheffler RISK assuming ND NE are 7 and 10!!!!
            if( xyz1.distance_squared( move.residue(jr).xyz(ja) ) <= 10.0 ) {
              // TR << "clash_check HIS pair of pairs! " << ir << "/" << ia << " " << jr << "/" << ja << std::endl;
              // outpose.dump_pdb("b.base.pdb");
              // ligs[idx1].pose.dump_pdb("fixd.pdb");
              // base.pose.dump_pdb("base.pdb");
              // move.dump_pdb("move.pdb");
              // utility_exit_with_message("lasdfgj");
              return;
            }
          }
        }
      }
    }

    TR << "found non-clashing homodimer match: " << ligs[idx1].tag+"___"+string_of(itrans) << std::endl;

    Size end1 = base.pose.n_residue();
    Size end2 = 2*end1; // homodimer case

    vector1<Size> matchres;
    matchres.push_back(ligs[idx1].rsd1);
    matchres.push_back(ligs[idx1].rsd2);
    numeric::xyzVector<Real> zncen;
    if(!mop->checkalign(fixd,move,zncen)) {
      TR << "BAD ALIGNMENT HOMODIMER, (NOT) DUMPING " << "align_"+ligs[idx1].tag+"___"+string_of(itrans)+".pdb" << std::endl;
      return;
      ozstream out("BAD_ALIGN.pdb");
      fixd.dump_pdb(out);
      move.dump_pdb(out);
      out.close();
    }
    // fixd.dump_pdb("fixd.pdb");
    // move.dump_pdb("move.pdb");

    ScoreFunctionOP sf = new core::scoring::symmetry::SymmetricScoreFunction(core::scoring::get_score_function());

    // fix crappy H placement ( needed for alignment)
    // fixd.dump_pdb("0_lig_fixd.pdb");
    // move.dump_pdb("0_lig_move.pdb");
    myoptH(move,sf);
    myoptH(fixd,sf);
    // fixd.dump_pdb("1_lig_fixd.pdb");
    // move.dump_pdb("1_lig_move.pdb");


    Pose tmppose = base.pose;
    rot_pose  (tmppose,rots[itrans].M);
    trans_pose(tmppose,rots[itrans].v);
    tmppose.append_residue_by_jump(base.pose.residue(1),1,"","",true);
    for(Size l = 2; l <= base.pose.n_residue(); ++l) tmppose.append_residue_by_bond(base.pose.residue(l));

    core::conformation::symmetry::SymmDataOP sd = core::pose::symmetry::symm_data_from_asym_pose(tmppose,end1);
    Pose sympose = base.pose;
    core::pose::symmetry::make_symmetric_pose(sympose,*sd);

    if(false) {
      sympose.replace_residue(ligs[idx1].rsd1,fixd.residue(1),true);
      sympose.replace_residue(ligs[idx1].rsd2,fixd.residue(2),true);
      utility::io::ozstream out("sympose.pdb");
      sympose.dump_pdb(out);
      out << "HETATM" << I(5,9999) << "  ZN   ZN Z" << I(4,994) << "    " << F(8,3,  zncen.x()) << F(8,3,  zncen.y()) << F(8,3,  zncen.z()) << F(6,2,1.0) << F(6,2,1.0) << std::endl;
      out.close();
      utility_exit_with_message("testing....");
    }

    for(Size j = 1; j <= native.n_residue(); ++j) sympose.replace_residue(j,native.residue(j),true);
    sympose.replace_residue(ligs[idx1].rsd1,fixd.residue(1),true);
    sympose.replace_residue(ligs[idx1].rsd2,fixd.residue(2),true);
    // sympose.dump_pdb("symtest.pdb");
    // utility_exit_with_message("testing...");

    // core::conformation::symmetry::SymmetryInfoCOP symminfo = core::pose::symmetry::symmetry_info(sympose);
    // Size end1 = symminfo->get_nres_subunit();
    // Size end2 = 2*end1;
    Pose outpose = sympose;

    Pose native = outpose;

    // // buried uns stuff
    // utility::vector1<std::pair<Size,Size> > sections;
    // sections.push_back( std::pair<Size,Size>(     1,end1) );
    // sections.push_back( std::pair<Size,Size>(end1+1,end2) );
    // utility::vector1< core::id::AtomID > bb_bur_uns_iface = get_bb_bur_uns_iface(outpose,sections,2.2);
    // if( bb_bur_uns_iface.size() > (Size)option[basic::options::OptionKeys::crossmatch::max_bur_uns]() ) {
    //  // std::string fn = "bur_uns_fail_"+ObjexxFCL::string_of(idx1)+"_"+ObjexxFCL::string_of(idx2)+"_"+ObjexxFCL::string_of(itrans)+"_"+".pdb";
    //  // TR << "bur_uns_fail_iface " << fn << " ";
    //  // for(Size i = 1; i <= bb_bur_uns_iface.size(); ++i) {
    //  //  TR << "( name " << outpose.residue(bb_bur_uns_iface[i].rsd()).atom_name(bb_bur_uns_iface[i].atomno()) << " and resi " << bb_bur_uns_iface[i].rsd() << ") or ";
    //  // }
    //  // TR << std::endl;
    //  // outpose.dump_pdb(fn);
    //  TR << "buried uns fail " << bb_bur_uns_iface.size() << std::endl;
    //  return;
    // }

    // if(true) { // within lines design and minimization
    // TR << "adding zn site csts" << std::endl;
    // { // add constraints on ZN site
    //  vector1<Size> ha(4), hr(4);
    //  ha[1] = b.ligs[idx2].align_info.r1d;
    //  ha[2] = b.ligs[idx2].align_info.r2d;
    //  ha[3] =   ligs[idx1].align_info.r1d;
    //  ha[4] =   ligs[idx1].align_info.r2d;
    //  hr[1] = b.ligs[idx2].rsd1;
    //  hr[2] = b.ligs[idx2].rsd2;
    //  hr[3] =   ligs[idx1].rsd1+end1;
    //  hr[4] =   ligs[idx1].rsd2+end1;
    //  Size ce1 = outpose.residue(hr[1]).atom_index("CE1");
    //  for(Size i = 1; i <= 4; ++i) {
    //    for(Size j = i+1; j <= 4; ++j) {
    //      using namespace core::scoring::constraints;
    //      Real d = 4.0;
    //      if( ligs[idx1].align_info.sqp == false ) d = 3.266566;
    //      outpose.add_constraint( new AtomPairConstraint( AtomID(ha[i],hr[i]), AtomID(ha[j],hr[j]), new HarmonicFunc(d,0.01) ) );
    //
    //      Real ang1 = angle_radians(      outpose.xyz( AtomID(ce1,hr[i])),outpose.xyz(AtomID(ha[i],hr[i])), outpose.xyz(AtomID(ha[j],hr[j])) );
    //      outpose.add_constraint( new AngleConstraint( AtomID(ce1,hr[i]),             AtomID(ha[i],hr[i]),              AtomID(ha[j],hr[j]), new HarmonicFunc(ang1,0.1) ) );
    //
    //      Real ang2 = angle_radians(      outpose.xyz( AtomID(ha[i],hr[i])), outpose.xyz(AtomID(ha[j],hr[j])), outpose.xyz(AtomID(ce1,hr[j])) );
    //      outpose.add_constraint( new AngleConstraint( AtomID(ha[i],hr[i]),              AtomID(ha[j],hr[j])   ,           AtomID(ce1,hr[j]), new HarmonicFunc(ang2,0.1) ) );
    //
    //      Real dih = dihedral_radians(        outpose.xyz(AtomID(ce1,hr[i])), outpose.xyz(AtomID(ha[i],hr[i])), outpose.xyz(AtomID(ha[j],hr[j])), outpose.xyz(AtomID(ce1,hr[j])) );
    //      outpose.add_constraint( new DihedralConstraint( AtomID(ce1,hr[i]) ,             AtomID(ha[i],hr[i]),              AtomID(ha[j],hr[j]),              AtomID(ce1,hr[j]), new CircularHarmonicFunc(dih,0.1) ) );
    //    }
    //  }
    // }

    // outpose.dump_pdb("3_premin.pdb");

    Real orig_rep = sf->get_weight(core::scoring::fa_rep);
    Size MIN_STEPS = 1;
    TR << "design/min " << MIN_STEPS << " rounds" << std::endl;
    // sf->set_weight(core::scoring::atom_pair_constraint, 1.0 );
    // sf->set_weight(core::scoring::    angle_constraint, 1.0 );
    // sf->set_weight(core::scoring:: dihedral_constraint, 1.0 );
    Size imin = 1;
    // for(Size imin = 1; imin <= MIN_STEPS; imin++) {
    sf->set_weight(core::scoring::fa_rep, orig_rep /((Real)MIN_STEPS*MIN_STEPS-(Real)(imin*imin)+1.0));
    design_homodimer(outpose,sf,matchres,false);
    // outpose.dump_pdb(string_of(imin+3)+"_design.pdb");
    // minimize(outpose,sf,matchres);
    // outpose.dump_pdb(string_of(imin+3)+"_min.pdb");
    // Size tmpnmut = 0; for(Size i = 1; i <= end2; ++i) if(outpose.residue(i).name3()!=native.residue(i).name3()) tmpnmut++;
    // TR << "design/min " << imin << " " << orig_rep /((Real)MIN_STEPS*MIN_STEPS-(Real)(imin*imin)+1.0) << " nmut: " << tmpnmut << std::endl;
    // }
    // design(outpose,sf,end1,matchres,end2,false);
    // outpose.dump_pdb("5_postmin.pdb");


    Size nmut  = 0; for(Size i = 1; i <= end2; ++i) if(outpose.residue(i).name3()!=native.residue(i).name3()) nmut++;
    Size nhmut = 0; for(Size i = 1; i <= end2; ++i) {
      if(outpose.residue(i).name3()==native.residue(i).name3()) continue;
      if(outpose.residue(i).name3()=="PHE" || outpose.residue(i).name3()=="ILE" || outpose.residue(i).name3()=="VAL" ||
         outpose.residue(i).name3()=="TRP" || outpose.residue(i).name3()=="MET" || outpose.residue(i).name3()=="TYR" ||
         outpose.residue(i).name3()=="LEU"  )
        nhmut++;
    }



    sf->set_weight(core::scoring::atom_pair_constraint, 0.0 );
    sf->set_weight(core::scoring::    angle_constraint, 0.0 );
    sf->set_weight(core::scoring:: dihedral_constraint, 0.0 );
    sf->score(outpose);

    // Real rholes = core::scoring::packstat::compute_packing_score( outpose , 0 );
    string fname = "ALIGN_"+ligs[idx1].tag+"___"+string_of(itrans)+".pdb";
    for(Size i = 0; i <= fname.size(); ++i) if(fname[i]==' ') fname[i] = '_';
    // TR << "CROSSMATCH_SCORE " << iface << " " << nhmut << " " << (*sf)(outpose) << " " << outpose.energies().total_energies()[core::scoring::fa_atr] << " " << rholes << " " << fname << std::endl;
    utility::io::ozstream out(option[basic::options::OptionKeys::out::file::o]()+"/"+fname);
    TR << "dumping pdb " << fname << std::endl;
    outpose.dump_pdb(out);
    Vec cen = ligs[idx1].align_info.cen;
    out << "HETATM" << I(5,9999) << "  ZN   ZN Z" << I(4,994) << "    " << F(8,3,  cen.x()) << F(8,3,  cen.y()) << F(8,3,  cen.z()) << F(6,2,1.0) << F(6,2,1.0) << std::endl;
    if(basic::options::option[basic::options::OptionKeys::smhybrid::add_cavities]()) {
      core::scoring::packstat::output_packstat_pdb(outpose,out);
    }
    out.close();



    core::io::silent::SilentStructOP ss_out( new core::io::silent::ScoreFileSilentStruct );
    ss_out->fill_struct(outpose,fname);
    // ss_out->add_energy( "iface", bb_bur_uns_iface.size() );
    // ss_out->add_energy( "bbburuns", bb_bur_uns_iface.size() );
    ss_out->add_energy( "nmut", nmut );
    ss_out->add_energy( "nhmut", nhmut );
    //TOOD count buried unsatisfied polars
    sfd.write_silent_struct( *ss_out, option[ basic::options::OptionKeys::out::file::silent ]() );

    // outpose.energies().show();

    // utility_exit_with_message("DEBUG");

    // std::ofstream out(("ALIGN_"+b.ligs[idx1].tag+"___"+b.ligs[idx2].tag+"___"+string_of(itrans)+".pdb").c_str());
    // outpose.dump_pdb(out);
    // {
    // Vec cen   = b.ligs[idx1].align_info.cen;
    // Vec axis  = cen + b.ligs[idx1].align_info.axis;
    // Vec ortho = cen + b.ligs[idx1].align_info.ortho;
    // out << "HETATM" << I(5,9994) << ' ' << "XXXX" << ' ' <<  "XXX" << ' ' << "A" << I(4,994) << "    " << F(8,3,  cen.x()) << F(8,3,  cen.y()) << F(8,3,  cen.z()) << F(6,2,1.0) << F(6,2,1.0) << '\n';
    // out << "HETATM" << I(5,9995) << ' ' << "XXXX" << ' ' <<  "XXX" << ' ' << "A" << I(4,995) << "    " << F(8,3, axis.x()) << F(8,3, axis.y()) << F(8,3, axis.z()) << F(6,2,1.0) << F(6,2,1.0) << '\n';
    // out << "HETATM" << I(5,9996) << ' ' << "XXXX" << ' ' << "XXX" << ' ' << "A" << I(4,996) << "    " << F(8,3,ortho.x()) << F(8,3,ortho.y()) << F(8,3,ortho.z()) << F(6,2,1.0) << F(6,2,1.0) << '\n';
    // }
    // {
    // Vec cen   = rots[itrans].M * b.ligs[idx2].align_info.cen + rots[itrans].v;
    // Vec axis  = cen + rots[itrans].M * b.ligs[idx2].align_info.axis;
    // Vec ortho = cen + rots[itrans].M * b.ligs[idx2].align_info.ortho;
    // out << "HETATM" << I(5,9997) << ' ' << "XXXX" << ' ' <<  "XXX" << ' ' << "A" << I(4,997) << "    " << F(8,3,  cen.x()) << F(8,3,  cen.y()) << F(8,3,  cen.z()) << F(6,2,1.0) << F(6,2,1.0) << '\n';
    // out << "HETATM" << I(5,9998) << ' ' << "XXXX" << ' ' <<  "XXX" << ' ' << "A" << I(4,998) << "    " << F(8,3, axis.x()) << F(8,3, axis.y()) << F(8,3, axis.z()) << F(6,2,1.0) << F(6,2,1.0) << '\n';
    // out << "HETATM" << I(5,9999) << ' ' << "XXXX" << ' ' << "XXX" << ' ' << "A" << I(4,999) << "    " << F(8,3,ortho.x()) << F(8,3,ortho.y()) << F(8,3,ortho.z()) << F(6,2,1.0) << F(6,2,1.0) << '\n';
    //      }
    // out.close();


  }
  void cross(MatchSet & b, Size idx1, Size idx2, Size itrans, Filter & filter) {
    if(!filter.filter(ligs[idx1].tag,b.ligs[idx2].tag,itrans)) {
      TR << "skipping " << ligs[idx1].tag << " " << b.ligs[idx2].tag << " " << itrans << std::endl;
      return;
    } else {
      // TR << "checking " << ligs[idx1].tag << " " << b.ligs[idx2].tag << " " << itrans << std::endl;
    }

    // TODO how do I build a new residue on a terminal spot? get upper/lower connect error
    // need lower/upper patch on types?
    {
      Size r1 = ligs[idx1].rsd1;
      Size r2 = ligs[idx1].rsd2;
      if( r1 == 1 ) return;
      if( r2 == 1 ) return;
      if( r1 == base.pose.n_residue() ) return;
      if( r2 == base.pose.n_residue() ) return;
      if( base.pose.chain(r1)!=base.pose.chain(r1-1) ) return;
      if( base.pose.chain(r1)!=base.pose.chain(r1+1) ) return;
      if( base.pose.chain(r2)!=base.pose.chain(r2-1) ) return;
      if( base.pose.chain(r2)!=base.pose.chain(r2+1) ) return;
    }
    {
      Size r1 = b.ligs[idx2].rsd1;
      Size r2 = b.ligs[idx2].rsd2;
      if( r1 == 1 ) return;
      if( r2 == 1 ) return;
      if( r1 == b.base.pose.n_residue() ) return;
      if( r2 == b.base.pose.n_residue() ) return;
      if( b.base.pose.chain(r1)!=b.base.pose.chain(r1-1) ) return;
      if( b.base.pose.chain(r1)!=b.base.pose.chain(r1+1) ) return;
      if( b.base.pose.chain(r2)!=b.base.pose.chain(r2-1) ) return;
      if( b.base.pose.chain(r2)!=b.base.pose.chain(r2+1) ) return;
    }
    // if(  base.pose.residue(  ligs[idx1].rsd1).is_lower_terminus() ||   base.pose.residue(  ligs[idx1].rsd1).is_upper_terminus()) return;
    // if(  base.pose.residue(  ligs[idx1].rsd2).is_lower_terminus() ||   base.pose.residue(  ligs[idx1].rsd2).is_upper_terminus()) return;
    // if(b.base.pose.residue(b.ligs[idx2].rsd1).is_lower_terminus() || b.base.pose.residue(b.ligs[idx2].rsd1).is_upper_terminus()) return;
    // if(b.base.pose.residue(b.ligs[idx2].rsd2).is_lower_terminus() || b.base.pose.residue(b.ligs[idx2].rsd2).is_upper_terminus()) return;

    vector1<Stub> rots = mop->align_rot(b.ligs[idx2].align_info,ligs[idx1].align_info);
    if(rots.size() < itrans) return;

    // if( b.base.pose.n_residue()!=407 || base.pose.n_residue()!=1956 ) {
    //  TR << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! NOT checking problem-specific props! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    // } else {
    // Vec hydsf4(0,0,0), psIsf4(0,0,0);
    // hydsf4 += rots[itrans].local2global(b.base.pose.xyz(AtomID(5, 101)))/4.0;
    // hydsf4 += rots[itrans].local2global(b.base.pose.xyz(AtomID(5, 156)))/4.0;
    // hydsf4 += rots[itrans].local2global(b.base.pose.xyz(AtomID(5, 334)))/4.0;
    // hydsf4 += rots[itrans].local2global(b.base.pose.xyz(AtomID(5, 338)))/4.0;
    // psIsf4 +=                             base.pose.xyz(AtomID(5,1305)) /4.0;
    // psIsf4 +=                             base.pose.xyz(AtomID(5, 563)) /4.0;
    // psIsf4 +=                             base.pose.xyz(AtomID(5,1314)) /4.0;
    // psIsf4 +=                             base.pose.xyz(AtomID(5, 572)) /4.0;
    // if(hydsf4.distance(psIsf4) > 23.0) {
    //  // TR << "found SF4 " << " " << hydsf4.distance(psIsf4) << std::endl;
    //  return;
    // } else {
    //  TR << "found match wiht SF4 within 30A " << " " << hydsf4.distance(psIsf4) << std::endl;
    // }
    // }

    if(!clash_check(b.base.pose,rots[itrans])) return;
    Pose outpose = b.base.pose;
    rot_pose  (outpose,rots[itrans].M);
    trans_pose(outpose,rots[itrans].v);
    Pose move     = b.ligs[idx2].pose;
    rot_pose  (move,rots[itrans].M);
    trans_pose(move,rots[itrans].v);

    // clash check base (self) vs move (other HISs)
    for(Size ia = 1; ia <= move.residue(1).nheavyatoms(); ++ia) {
      if( !this->clash_check(move.residue(1).xyz(ia)) ||
          !this->clash_check(move.residue(2).xyz(ia)) ) {
        // TR << "clash_check on move HIS HIS fails!" << std::endl;
        // base.pose.dump_pdb("base.pdb");
        // move.dump_pdb("move.pdb");
        // utility_exit_with_message("lasdfgj");
        return;
      }
    }

    // clash check other base vs fixd (our HISs)
    for(Size ir = 1; ir <= 2; ++ir) {
      for(Size ia = 1; ia <= ligs[idx1].pose.residue(ir).nheavyatoms(); ++ia) {
        if( !b.clash_check_inv(ligs[idx1].pose.residue(ir).xyz(ia),rots[itrans]) ) {
          // TR << "clash_check on fixd HIS HIS fails! " << ir << " " << ia << std::endl;
          // outpose.dump_pdb("b.base.pdb");
          // ligs[idx1].pose.dump_pdb("fixd.pdb");
          // base.pose.dump_pdb("base.pdb");
          // move.dump_pdb("move.pdb");
          // utility_exit_with_message("lasdfgj");
          return;
        }
      }
    }

    for(Size ir = 1; ir <= 2; ++ir) {
      for(Size ia = 1; ia <= ligs[idx1].pose.residue(ir).nheavyatoms(); ++ia) {
        if( 7==ia || 10==ia ) continue; //sheffler RISK assuming ND NE are 7 and 10!!!!
        Vec const xyz1 = ligs[idx1].pose.residue(ir).xyz(ia);
        for(Size jr = 1; jr <= 2; ++jr) {
          for(Size ja = 1; ja <= move.residue(jr).nheavyatoms(); ++ja) {
            if( 7==ja || 10==ja ) continue; //sheffler RISK assuming ND NE are 7 and 10!!!!
            if( xyz1.distance_squared( move.residue(jr).xyz(ja) ) <= 10.0 ) {
              // TR << "clash_check HIS pair of pairs! " << ir << "/" << ia << " " << jr << "/" << ja << std::endl;
              // outpose.dump_pdb("b.base.pdb");
              // ligs[idx1].pose.dump_pdb("fixd.pdb");
              // base.pose.dump_pdb("base.pdb");
              // move.dump_pdb("move.pdb");
              // utility_exit_with_message("lasdfgj");
              return;
            }
          }
        }
      }
    }


    vector1<Size> include_res = b.iface_candidates;
    for(vector1<Size>::const_iterator i=iface_candidates.begin(),ie=iface_candidates.end(); i != ie; ++i) {
      include_res.push_back(*i+b.base.pose.n_residue());
    }
    Real iface = iface_check(outpose,b.iface_candidates);
    if(iface < (Real)basic::options::option[basic::options::OptionKeys::willmatch::interface_size]()) return;

    TR << "found non-clashing matches: " << ligs[idx1].tag+"___"+b.ligs[idx2].tag+"___"+string_of(itrans) << std::endl;

    TR << "build joint pose" << std::endl;
    Pose fixd     = ligs[idx1].pose;
    Pose fixdbase = base.pose;
    outpose.append_residue_by_jump(fixdbase.residue(1),1,"","",true);
    for(Size l = 2; l <= fixdbase.n_residue(); ++l) outpose.append_residue_by_bond(fixdbase.residue(l));

    Size end1 = b.base.pose.n_residue();
    vector1<Size> matchres;
    matchres.push_back(b.ligs[idx2].rsd1);
    matchres.push_back(b.ligs[idx2].rsd2);
    matchres.push_back(ligs[idx1].rsd1+end1);
    matchres.push_back(ligs[idx1].rsd2+end1);
    numeric::xyzVector<Real> zncen;
    if(!mop->checkalign(fixd,move,zncen)) {
      TR << "BAD ALIGNMENT, DUMPING " << "align_"+b.ligs[idx1].tag+"___"+b.ligs[idx2].tag+"_"+string_of(itrans)+".pdb" << std::endl;
      ozstream out("BAD_ALIGN.pdb");
      fixd.dump_pdb(out);
      move.dump_pdb(out);
      out.close();
    }
    using core::kinematics::Stub;
    Real LINK_CST = option[basic::options::OptionKeys::willmatch::c2_linker_dist]();
    Real const c2_incr = option[basic::options::OptionKeys::willmatch::c2_symm_increment]();

    if(option[basic::options::OptionKeys::willmatch::symmetry_c2_dock]()) {
      utility_exit_with_message("symmetry_c2_dock may not work!!!! check code!!!");
      vector1<Vec> sdots = core::scoring::packstat::get_sasa_dot_locations();
      Real best_iface = 0;
      Stub best_stub;
      bool found = false;
      vector1<Vec> his_atoms;
      for(Size i = 6; i <= fixd.residue(1).nheavyatoms(); ++i) {
        his_atoms.push_back(fixd.xyz(AtomID(i,1)));
        his_atoms.push_back(fixd.xyz(AtomID(i,2)));
        his_atoms.push_back(rots[itrans].M*move.xyz(AtomID(i,1))+rots[itrans].v);
        his_atoms.push_back(rots[itrans].M*move.xyz(AtomID(i,2))+rots[itrans].v);
      }
      for(Size idot = 1; idot <= sdots.size(); ++idot) {
        Stub s( rotation_matrix_degrees( sdots[idot], 180.0 ), Vec(0,0,0) );
        Vec trans0 = sdots[idot].cross(Vec(0.77183184,0.05846532,0.54337440)).normalized();
        for(Real ic2 = 0; ic2 < 360; ic2 += c2_incr ) {
          Vec trans = rotation_matrix_degrees(sdots[idot],ic2) * trans0;
          // TR << "checking c2 " << idot << " " << ic2 << std::endl;
          if(!get_contacting_stub(outpose,trans,s,rots[itrans],b,his_atoms,LINK_CST)) continue;
          Real iface = c2_iface_check(outpose,s,include_res);// - c2_linker_check_dist(outpose,b,s)/3.0;
          if( iface > best_iface ) {
            best_stub = s;
            best_iface = iface;
            found = true;
          }
          // TR << "IFACE " << iface << std::endl;
        }
      }
      if(!found) {
        TR << "COULDNT FIND C2 MATCH" << std::endl;
        return;
      }
      TR << "best_iface BEFORE perterb " << best_iface << std::endl;
      Vec com1 = com(outpose);
      for(Size i = 1; i <= 1000; ++i) {
        Real scale = (Real)i/100. - 5.0;
        Stub s = best_stub;
        Real theta;
        Vec axis = rotation_axis(s.M,theta);
        assert((fabs(theta)-180) < 0.001);
        Vec trans = (s.v-(com1-s.M*com1)).normalized();
        // TR << "axis " << axis << " trans " << trans << " " << axis.dot(trans) << std::endl;
        axis = rotation_matrix_degrees(axis.cross(randvec()).normalized(),uniform()*scale) * axis;
        s.M = rotation_matrix_degrees(axis,180.0);
        trans = rotation_matrix_degrees(trans.cross(randvec()).normalized(),uniform()*scale/2.0) * trans;
        trans = (trans-axis.dot(trans)*axis).normalized();
        // TR << "axis " << axis << " trans " << trans << std::endl;
        if(!get_contacting_stub(outpose,trans,s,rots[itrans],b,his_atoms,9999.0)) continue;
        Real iface = c2_iface_check(outpose,s,include_res);// - c2_linker_check_dist(outpose,b,s)/3.0;
        // TR << "best_iface WHILE  " << iface << std::endl;
        if(iface > best_iface) {
          best_stub = s;
          best_iface = iface;
        }
      }
      TR << "best_iface AFTER perterb " << best_iface << std::endl;
      // std::exit(-1);
      // Pose primary = outpose;
      // Pose partner = outpose;
      // rot_pose(partner,best_stub.M);
      // trans_pose(partner,best_stub.v);
      // primary.append_residue_by_jump(partner.residue(1),1,"","",true);
      // for(Size l = 2; l <= b.base.pose.n_residue(); ++l) {
      //  primary.append_residue_by_bond(partner.residue(l));
      // }
      // primary.append_residue_by_jump(partner.residue(b.base.pose.n_residue()+1),1,"","",true);
      // for(Size l = b.base.pose.n_residue()+2; l <= partner.n_residue(); ++l) {
      //  primary.append_residue_by_bond(partner.residue(l));
      // }
      // ozstream out("test.pdb");
      // primary.dump_pdb(out);
      // out.close();

      Real ang;
      Vec symaxis = rotation_axis(best_stub.M,ang);
      assert(fabs(ang-radians(180.0))<0.0001);
      trans_pose(outpose,-best_stub.v/2.0);
      if(symaxis != Vec(1,0,0)) {
        rot_pose(outpose,rotation_matrix(symaxis.cross(Vec(1,0,0)),acos(symaxis.x())));
      }
      core::pose::symmetry::make_symmetric_pose(outpose);
      // // rot_pose(outpose,rotation_matrix(symaxis.cross(Vec(1,0,0)),-acos(symaxis.x())));
      // // trans_pose(outpose, best_stub.v/2.0);
      // outpose.dump_pdb("symm.pdb");
      // std::exit(-1);
    }

    ScoreFunctionOP sf = core::scoring::get_score_function();

    Size end2 = base.pose.n_residue()+b.base.pose.n_residue();

    // fix crappy H placement ( needed for alignment)
    // fixd.dump_pdb("0_lig_fixd.pdb");
    // move.dump_pdb("0_lig_move.pdb");
    myoptH(move,sf);
    myoptH(fixd,sf);
    // fixd.dump_pdb("1_lig_fixd.pdb");
    // move.dump_pdb("1_lig_move.pdb");

    TR << "checking sasa and buried unsat polar" << std::endl;
    protocols::toolbox::switch_to_residue_type_set( outpose, "FA_STANDARD" );
    repack(outpose,sf,end1,utility::vector1<Size>(),end2);
    outpose.replace_residue(   b.ligs[idx2].rsd1,move.residue(1),true);
    outpose.replace_residue(   b.ligs[idx2].rsd2,move.residue(2),true);
    outpose.replace_residue(end1+ligs[idx1].rsd1,fixd.residue(1),true);
    outpose.replace_residue(end1+ligs[idx1].rsd2,fixd.residue(2),true);

    // buried uns stuff
    utility::vector1<std::pair<Size,Size> > sections;
    sections.push_back( std::pair<Size,Size>(     1,end1) );
    sections.push_back( std::pair<Size,Size>(end1+1,end2) );
    utility::vector1< core::id::AtomID > bb_bur_uns_iface = get_bb_bur_uns_iface(outpose,sections,2.2);
    if( bb_bur_uns_iface.size() > (Size)option[basic::options::OptionKeys::crossmatch::max_bur_uns]() ) {
      // std::string fn = "bur_uns_fail_"+ObjexxFCL::string_of(idx1)+"_"+ObjexxFCL::string_of(idx2)+"_"+ObjexxFCL::string_of(itrans)+"_"+".pdb";
      // TR << "bur_uns_fail_iface " << fn << " ";
      // for(Size i = 1; i <= bb_bur_uns_iface.size(); ++i) {
      //  TR << "( name " << outpose.residue(bb_bur_uns_iface[i].rsd()).atom_name(bb_bur_uns_iface[i].atomno()) << " and resi " << bb_bur_uns_iface[i].rsd() << ") or ";
      // }
      // TR << std::endl;
      // outpose.dump_pdb(fn);
      TR << "buried uns fail " << bb_bur_uns_iface.size() << std::endl;
      return;
    }
    // outpose.dump_pdb("bur_uns_pass_"+ObjexxFCL::string_of(idx1)+"_"+ObjexxFCL::string_of(idx2)+"_"+ObjexxFCL::string_of(itrans)+"_"+".pdb");
    // return;

    TR << "placing HIS and optH" << std::endl;
    for(Size j = 1; j <= b.native.n_residue(); ++j) outpose.replace_residue(j     ,b.native.residue(j),true);
    for(Size j = 1; j <=   native.n_residue(); ++j) outpose.replace_residue(j+end1,  native.residue(j),true);
    outpose.replace_residue(   b.ligs[idx2].rsd1,move.residue(1),true);
    outpose.replace_residue(   b.ligs[idx2].rsd2,move.residue(2),true);
    outpose.replace_residue(end1+ligs[idx1].rsd1,fixd.residue(1),true);
    outpose.replace_residue(end1+ligs[idx1].rsd2,fixd.residue(2),true);



    Pose native = outpose;

    // if(true) { // within lines design and minimization
    // TR << "adding zn site csts" << std::endl;
    // { // add constraints on ZN site
    //  vector1<Size> ha(4), hr(4);
    //  ha[1] = b.ligs[idx2].align_info.r1d;
    //  ha[2] = b.ligs[idx2].align_info.r2d;
    //  ha[3] =   ligs[idx1].align_info.r1d;
    //  ha[4] =   ligs[idx1].align_info.r2d;
    //  hr[1] = b.ligs[idx2].rsd1;
    //  hr[2] = b.ligs[idx2].rsd2;
    //  hr[3] =   ligs[idx1].rsd1+end1;
    //  hr[4] =   ligs[idx1].rsd2+end1;
    //  Size ce1 = outpose.residue(hr[1]).atom_index("CE1");
    //  for(Size i = 1; i <= 4; ++i) {
    //    for(Size j = i+1; j <= 4; ++j) {
    //      using namespace core::scoring::constraints;
    //      Real d = 4.0;
    //      if( ligs[idx1].align_info.sqp == false ) d = 3.266566;
    //      outpose.add_constraint( new AtomPairConstraint( AtomID(ha[i],hr[i]), AtomID(ha[j],hr[j]), new HarmonicFunc(d,0.01) ) );
    //
    //      Real ang1 = angle_radians(      outpose.xyz( AtomID(ce1,hr[i])),outpose.xyz(AtomID(ha[i],hr[i])), outpose.xyz(AtomID(ha[j],hr[j])) );
    //      outpose.add_constraint( new AngleConstraint( AtomID(ce1,hr[i]),             AtomID(ha[i],hr[i]),              AtomID(ha[j],hr[j]), new HarmonicFunc(ang1,0.1) ) );
    //
    //      Real ang2 = angle_radians(      outpose.xyz( AtomID(ha[i],hr[i])), outpose.xyz(AtomID(ha[j],hr[j])), outpose.xyz(AtomID(ce1,hr[j])) );
    //      outpose.add_constraint( new AngleConstraint( AtomID(ha[i],hr[i]),              AtomID(ha[j],hr[j])   ,           AtomID(ce1,hr[j]), new HarmonicFunc(ang2,0.1) ) );
    //
    //      Real dih = dihedral_radians(        outpose.xyz(AtomID(ce1,hr[i])), outpose.xyz(AtomID(ha[i],hr[i])), outpose.xyz(AtomID(ha[j],hr[j])), outpose.xyz(AtomID(ce1,hr[j])) );
    //      outpose.add_constraint( new DihedralConstraint( AtomID(ce1,hr[i]) ,             AtomID(ha[i],hr[i]),              AtomID(ha[j],hr[j]),              AtomID(ce1,hr[j]), new CircularHarmonicFunc(dih,0.1) ) );
    //    }
    //  }
    // }

    // outpose.dump_pdb("3_premin.pdb");

    Real orig_rep = sf->get_weight(core::scoring::fa_rep);
    Size MIN_STEPS = 1;
    TR << "design/min " << MIN_STEPS << " rounds" << std::endl;
    // sf->set_weight(core::scoring::atom_pair_constraint, 1.0 );
    // sf->set_weight(core::scoring::    angle_constraint, 1.0 );
    // sf->set_weight(core::scoring:: dihedral_constraint, 1.0 );
    Size imin = 1;
    // for(Size imin = 1; imin <= MIN_STEPS; imin++) {
    sf->set_weight(core::scoring::fa_rep, orig_rep /((Real)MIN_STEPS*MIN_STEPS-(Real)(imin*imin)+1.0));
    design(outpose,sf,end1,matchres,end2,false);
    // outpose.dump_pdb(string_of(imin+3)+"_design.pdb");
    // minimize(outpose,sf,matchres);
    // outpose.dump_pdb(string_of(imin+3)+"_min.pdb");
    Size tmpnmut = 0; for(Size i = 1; i <= end2; ++i) if(outpose.residue(i).name3()!=native.residue(i).name3()) tmpnmut++;
    TR << "design/min " << imin << " " << orig_rep /((Real)MIN_STEPS*MIN_STEPS-(Real)(imin*imin)+1.0) << " nmut: " << tmpnmut << std::endl;
    // }
    // design(outpose,sf,end1,matchres,end2,false);
    // outpose.dump_pdb("5_postmin.pdb");


    Size nmut  = 0; for(Size i = 1; i <= end2; ++i) if(outpose.residue(i).name3()!=native.residue(i).name3()) nmut++;
    Size nhmut = 0; for(Size i = 1; i <= end2; ++i) {
      if(outpose.residue(i).name3()==native.residue(i).name3()) continue;
      if(outpose.residue(i).name3()=="PHE" || outpose.residue(i).name3()=="ILE" || outpose.residue(i).name3()=="VAL" ||
         outpose.residue(i).name3()=="TRP" || outpose.residue(i).name3()=="MET" || outpose.residue(i).name3()=="TYR" ||
         outpose.residue(i).name3()=="LEU"  )
        nhmut++;
    }

    bool reset_surface_mutations = false;
    Size revert_iter = 0;
    while(reset_surface_mutations) {
      TR << "reverting surface mutations" << std::endl;
      vector1<Size> reset_res;
      core::id::AtomID_Map< bool > atom_map;
      core::pose::initialize_atomid_map( atom_map, outpose, false );
      for ( Size ir = 1; ir <= outpose.n_residue(); ++ir ) {
        for(Size ia = 1; ia <= outpose.residue(ir).nheavyatoms(); ++ia) {
          Size t = outpose.residue(ir).atom_type_index(ia);
          if( t==3 || t==4 || t==5 || t==6 || t==19 || t==20 ) atom_map.set(AtomID(ia,ir) , true );
        }
      }
      core::id::AtomID_Map<Real> atom_sasa; utility::vector1<Real> sasa;
      core::scoring::calc_per_atom_sasa( outpose, atom_sasa, sasa, 2.0, false, atom_map );
      for(Size i = 1; i <= end2; ++i) {
        if(outpose.residue(i).name3()=="PHE" && sasa[i] > 15.0 ) { /*TR << "SASA PHE " << sasa[i] << std::endl;*/ reset_res.push_back(i); }
        if(outpose.residue(i).name3()=="TRP" && sasa[i] > 20.0 ) { /*TR << "SASA TRP " << sasa[i] << std::endl;*/ reset_res.push_back(i); }
        if(outpose.residue(i).name3()=="TYR" && sasa[i] > 20.0 ) { /*TR << "SASA TYR " << sasa[i] << std::endl;*/ reset_res.push_back(i); }
        if(outpose.residue(i).name3()=="ILE" && sasa[i] > 15.0 ) { /*TR << "SASA ILE " << sasa[i] << std::endl;*/ reset_res.push_back(i); }
        if(outpose.residue(i).name3()=="LEU" && sasa[i] > 15.0 ) { /*TR << "SASA LEU " << sasa[i] << std::endl;*/ reset_res.push_back(i); }
        if(outpose.residue(i).name3()=="VAL" && sasa[i] > 15.0 ) { /*TR << "SASA VAL " << sasa[i] << std::endl;*/ reset_res.push_back(i); }
        if(outpose.residue(i).name3()=="MET" && sasa[i] > 15.0 ) { /*TR << "SASA MET " << sasa[i] << std::endl;*/ reset_res.push_back(i); }
        if(outpose.residue(i).name3()=="ALA" && sasa[i] >  5.0 ) { /*TR << "SASA ALA " << sasa[i] << std::endl;*/ reset_res.push_back(i); }
      }
      Size reset = 0;
      for(Size i = 1; i <= reset_res.size(); ++i) {
        if(outpose.residue(reset_res[i]).name3()=="HIS" || outpose.residue(reset_res[i]).name3()==native.residue(reset_res[i]).name3()) continue;
        TR << "reverting res " << reset_res[i] << " " << outpose.residue(reset_res[i]).name3() << std::endl;
        core::pose::replace_pose_residue_copying_existing_coordinates(outpose,reset_res[i],native.residue(reset_res[i]).type());
        reset++;
      }
      if(reset==0) break;
      repack(outpose,sf,end1,matchres,end2);
      minimize(outpose,sf,matchres);
      revert_iter++;
      // outpose.dump_pdb("7_revert_"+string_of(revert_iter)+".pdb");
    }
    nmut  = 0; for(Size i = 1; i <= end2; ++i) if(outpose.residue(i).name3()!=native.residue(i).name3()) nmut++;
    nhmut = 0; for(Size i = 1; i <= end2; ++i) {
      if(outpose.residue(i).name3()==native.residue(i).name3()) continue;
      if(outpose.residue(i).name3()=="PHE" || outpose.residue(i).name3()=="ILE" || outpose.residue(i).name3()=="VAL" ||
         outpose.residue(i).name3()=="TRP" || outpose.residue(i).name3()=="MET" || outpose.residue(i).name3()=="TYR" ||
         outpose.residue(i).name3()=="LEU"  )
        nhmut++;
    }
    TR << "done, num mutations " << nmut << " hydrophbic " << nhmut << std::endl;
    outpose.dump_pdb("9_done.pdb");
    // }
    // if(nhmut < 5) {
    // TR << "not enough useful mutations" << std::endl;
    //  return;
    // }

    sf->set_weight(core::scoring::atom_pair_constraint, 0.0 );
    sf->set_weight(core::scoring::    angle_constraint, 0.0 );
    sf->set_weight(core::scoring:: dihedral_constraint, 0.0 );
    sf->score(outpose);

    // Real rholes = core::scoring::packstat::compute_packing_score( outpose , 0 );
    string fname = "ALIGN_"+ligs[idx1].tag+"___"+b.ligs[idx2].tag+"___"+string_of(itrans)+".pdb";
    for(Size i = 0; i <= fname.size(); ++i) if(fname[i]==' ') fname[i] = '_';
    // TR << "CROSSMATCH_SCORE " << iface << " " << nhmut << " " << (*sf)(outpose) << " " << outpose.energies().total_energies()[core::scoring::fa_atr] << " " << rholes << " " << fname << std::endl;
    utility::io::ozstream out(option[basic::options::OptionKeys::out::file::o]()+"/"+fname);
    TR << "dumping pdb " << fname << std::endl;
    outpose.dump_pdb(out);
    Vec cen = ((ligs[idx1].align_info.cen) + (rots[itrans].M * b.ligs[idx2].align_info.cen + rots[itrans].v) ) / 2.0;
    out << "HETATM" << I(5,9999) << "  ZN   ZN Z" << I(4,994) << "    " << F(8,3,  cen.x()) << F(8,3,  cen.y()) << F(8,3,  cen.z()) << F(6,2,1.0) << F(6,2,1.0) << std::endl;
    if(basic::options::option[basic::options::OptionKeys::smhybrid::add_cavities]()) {
      core::scoring::packstat::output_packstat_pdb(outpose,out);
    }
    out.close();

    core::io::silent::SilentStructOP ss_out( new core::io::silent::ScoreFileSilentStruct );
    ss_out->fill_struct(outpose,fname);
    ss_out->add_energy( "nmut", nmut );
    ss_out->add_energy( "nhmut", nhmut );
    //TOOD count buried unsatisfied polars
    sfd.write_silent_struct( *ss_out, option[ basic::options::OptionKeys::out::file::silent ]() );

    outpose.energies().show();

    // utility_exit_with_message("DEBUG");

    // std::ofstream out(("ALIGN_"+b.ligs[idx1].tag+"___"+b.ligs[idx2].tag+"___"+string_of(itrans)+".pdb").c_str());
    // outpose.dump_pdb(out);
    // {
    // Vec cen   = b.ligs[idx1].align_info.cen;
    // Vec axis  = cen + b.ligs[idx1].align_info.axis;
    // Vec ortho = cen + b.ligs[idx1].align_info.ortho;
    // out << "HETATM" << I(5,9994) << ' ' << "XXXX" << ' ' <<  "XXX" << ' ' << "A" << I(4,994) << "    " << F(8,3,  cen.x()) << F(8,3,  cen.y()) << F(8,3,  cen.z()) << F(6,2,1.0) << F(6,2,1.0) << '\n';
    // out << "HETATM" << I(5,9995) << ' ' << "XXXX" << ' ' <<  "XXX" << ' ' << "A" << I(4,995) << "    " << F(8,3, axis.x()) << F(8,3, axis.y()) << F(8,3, axis.z()) << F(6,2,1.0) << F(6,2,1.0) << '\n';
    // out << "HETATM" << I(5,9996) << ' ' << "XXXX" << ' ' << "XXX" << ' ' << "A" << I(4,996) << "    " << F(8,3,ortho.x()) << F(8,3,ortho.y()) << F(8,3,ortho.z()) << F(6,2,1.0) << F(6,2,1.0) << '\n';
    // }
    // {
    // Vec cen   = rots[itrans].M * b.ligs[idx2].align_info.cen + rots[itrans].v;
    // Vec axis  = cen + rots[itrans].M * b.ligs[idx2].align_info.axis;
    // Vec ortho = cen + rots[itrans].M * b.ligs[idx2].align_info.ortho;
    // out << "HETATM" << I(5,9997) << ' ' << "XXXX" << ' ' <<  "XXX" << ' ' << "A" << I(4,997) << "    " << F(8,3,  cen.x()) << F(8,3,  cen.y()) << F(8,3,  cen.z()) << F(6,2,1.0) << F(6,2,1.0) << '\n';
    // out << "HETATM" << I(5,9998) << ' ' << "XXXX" << ' ' <<  "XXX" << ' ' << "A" << I(4,998) << "    " << F(8,3, axis.x()) << F(8,3, axis.y()) << F(8,3, axis.z()) << F(6,2,1.0) << F(6,2,1.0) << '\n';
    // out << "HETATM" << I(5,9999) << ' ' << "XXXX" << ' ' << "XXX" << ' ' << "A" << I(4,999) << "    " << F(8,3,ortho.x()) << F(8,3,ortho.y()) << F(8,3,ortho.z()) << F(6,2,1.0) << F(6,2,1.0) << '\n';
    //      }
    // out.close();

  }
};

std::pair<vector1<Size>,vector1<Size> > makesplitwork(Size total, Size total2 = 0) {
  using namespace basic::options::OptionKeys;
  vector1<Size> idx1;
  vector1<Size> idx2;
  Size part   = 1;
  Size nparts = 1;
  if( option[willmatch::splitwork].user() ) {
    part   = option[willmatch::splitwork]()[1];
    nparts = option[willmatch::splitwork]()[2];
  }
  TR << "makesplitwork " << total << " " << total2 << " " << part << " " << nparts << std::endl;
  if( total2 == 0 ){
    for( Size i = part; i <= total; i += nparts ) {
      idx1.push_back(i);
    }
  } else {
    if( option[in::file::s]().size() == 1 ) {
      for( Size i = part; i <= total; i += nparts ) {
        for( Size j = i; j <= total; j += 1 ) {
          idx1.push_back(i);
          idx2.push_back(j);
        }
      }
    } else {
      for( Size i = part; i <= total; i += nparts ) {
        for( Size j = 1; j <= total2; j += 1 ) {
          idx1.push_back(i);
          idx2.push_back(j);
        }
      }
    }
  }
  return std::pair<vector1<Size>,vector1<Size> >(idx1,idx2);
}

int main (int argc, char *argv[])
{

	try {

  using namespace basic::options::OptionKeys;

  devel::init(argc,argv);

  MatchAlignerOP mop = new Tet4HMatchAligner;
  // MatchSet ms2(option[in::file::s],mop);

  if(basic::options::option[basic::options::OptionKeys::willmatch::write_reduced_matchset].user()) {
    //
    vector1<string> v = basic::options::option[basic::options::OptionKeys::willmatch::write_reduced_matchset]();
    string fname = v[1];
    vector1<Pose> poses;
    poses.resize(v.size()-1);
    for(Size i = 2; i <= v.size(); ++i){
      TR << "reading pose " << v[i] << std::endl;
      core::import_pose::pose_from_pdb(poses[i-1],v[i]);
    }
    MatchSet ms1(poses,mop,poses[1]);
    ms1.write_to_file(fname);
  } else
    if(basic::options::option[basic::options::OptionKeys::willmatch::symmetry_d2]()) {
      core::chemical::ResidueTypeSetCAP  fa_residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
      vector1< Pose > poses;
      Pose native;
      core::import_pose::pose_from_pdb(poses,option[in::file::s][1]);
      core::import_pose::pose_from_pdb(native,option[willmatch::native1]());
      MatchSet ms(poses,mop,native);
      std::pair<vector1<Size>,vector1<Size> > splitwork = makesplitwork(ms.ligs.size());
      vector1<Size> ILIG1 = splitwork.first;
      vector1<Size> ILIG2 = splitwork.second;
      assert(ILIG2.size()==ILIG1.size());
      vector1<Pose> his(2);
      core::pose::make_pose_from_sequence(his[1],"H[HIS_D]",*fa_residue_set,false);
      core::pose::make_pose_from_sequence(his[2],"H[HIS]"  ,*fa_residue_set,false);
      core::pose::remove_lower_terminus_type_from_pose_residue(his[1],1);
      core::pose::remove_lower_terminus_type_from_pose_residue(his[2],1);
      core::pose::remove_upper_terminus_type_from_pose_residue(his[1],1);
      core::pose::remove_upper_terminus_type_from_pose_residue(his[2],1);
      his[1].set_dof(core::id::DOF_ID(AtomID(his[1].residue(1).atom_index("HD1"),1),core::id::D),2.1);
      his[2].set_dof(core::id::DOF_ID(AtomID(his[2].residue(1).atom_index("HE2"),1),core::id::D),2.1);
      core::conformation::Residue vrt(fa_residue_set->name_map("VRT"),true);
      his[1].append_residue_by_jump(vrt,1,"HD1","ORIG");
      his[2].append_residue_by_jump(vrt,1,"HE2","ORIG");
      FoldTree ft = his[1].fold_tree(); ft.reorder(2); his[1].fold_tree(ft);
      ft = his[2].fold_tree();  ft.reorder(2); his[2].fold_tree(ft);
      Real chi1incr = option[willmatch::chi1_increment]();
      Real chi2incr = option[willmatch::chi2_increment]();
      vector1<Real> CHI1,CHI2;
      for(Real i = 0; i < 360; i+= chi1incr) CHI1.push_back(i);
      for(Real i = 0; i < 360; i+= chi2incr) CHI2.push_back(i);
      Real const HIS_CLASH_DIS2 = 10.0;
      // ObjexxFCL::FArray3D<Vec>   chi2axis(2,CHI1.size(),CHI2.size());
      // ObjexxFCL::FArray3D<Vec>  chi2trans(2,CHI1.size(),CHI2.size());
      // ObjexxFCL::FArray3D<Real> chi2angle(2,CHI1.size(),CHI2.size());
      // Stub stub(his[2].xyz(AtomID(1,1)),his[2].xyz(AtomID(2,1)),his[2].xyz(AtomID(3,1)));
      // Real angle;
      // Vec axis = rotation_axis(stub.M,angle);
      // for(Size i = 1; i <= CHI1.size(); ++i) {
      //  his[1].set_chi(1,1,CHI1[i]);
      //  his[2].set_chi(1,1,CHI1[i]);
      //  for(Size j = 1; j <= CHI2.size(); ++j) {
      //    his[2].set_chi(2,1,CHI2[j]);
      //    his[1].set_chi(2,1,CHI2[j]);
      //    Stub stubd(his[1].residue(1).xyz("HD1"),his[2].residue(1).xyz("ND1"),his[2].residue(1).xyz("CE1"));
      //    Stub stube(his[2].residue(1).xyz("HE2"),his[2].residue(1).xyz("NE2"),his[2].residue(1).xyz("CE1"));
      //    chi2axis (1,i,j) = stubd.global2local( axis );
      //    chi2axis (2,i,j) = stube.global2local( axis );
      //    chi2angle(1,i,j) = angle;
      //    chi2angle(2,i,j) = angle;
      //    chi2trans(1,i,j) = stubd.global2local( stub.v );
      //    chi2trans(2,i,j) = stube.global2local( stub.v );
      //  }
      // }
      for(Size iligidx = 1; iligidx <= ILIG1.size(); ++iligidx) {
        Size ilig = ILIG1[iligidx];
        TR << "symm d2 with match " << ilig << std::endl;
        // if(ilig < 113) {
        //  TR << "!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << std::endl;
        //  continue;
        // }
        MatchAlignInfo & mai(ms.ligs[ilig].align_info);
        Pose tmppose = ms.base.pose;
        Size rsd1 = ms.ligs[ilig].rsd1, rsd2 = ms.ligs[ilig].rsd2;
        tmppose.replace_residue(rsd1,ms.ligs[ilig].pose.residue(1),true);
        tmppose.replace_residue(rsd2,ms.ligs[ilig].pose.residue(2),true);
        core::pose::remove_lower_terminus_type_from_pose_residue(tmppose,1);
        core::pose::remove_upper_terminus_type_from_pose_residue(tmppose,tmppose.n_residue());
        // tmppose.set_xyz(AtomID(tmppose.residue(rsd1).atom_index("HE2"),rsd1),mai.cen);
        // tmppose.set_xyz(AtomID(tmppose.residue(rsd2).atom_index("HE2"),rsd2),mai.cen);
        Pose const pose(tmppose);
        // pose.dump_pdb("init.pdb");
        // mai.dump_pdb("mai.pdb");
        for(Size ihis = 2; ihis <= 2; ihis++) {
          Pose h = his[ihis];
          // h.set_chi(1,1,10);
          // h.dump_pdb("test1.pdb");
          // h.set_chi(1,1,20);
          // h.dump_pdb("test2.pdb");
          // utility_exit_with_message("HTEST");
          Size hatm,natm;
          if( ihis == 1 ) {
            hatm = h.residue(1).atom_index("HD1");
            natm = h.residue(1).atom_index("ND1");
          } else {
            hatm = h.residue(1).atom_index("HE2");
            natm = h.residue(1).atom_index("NE2");
          }
          if(mai.sqp) {
            // TR << "SQP" << std::endl;
            Vec v = (h.residue(1).xyz(hatm)-h.residue(1).xyz(natm)).normalized();
            Mat r1 = rotation_matrix( v.cross(-mai.axis), acos(v.dot(-mai.axis)) );
            Mat r2 = rotation_matrix( mai.ortho, numeric::conversions::radians(45.0) );
            rot_pose(h,r2*r1);
            trans_pose(h,mai.cen-h.residue(1).xyz(hatm));
          } else {
            // TR << "TET" << std::endl;
            Vec v = (h.residue(1).xyz(hatm)-h.residue(1).xyz(natm)).normalized();
            Mat r1 = rotation_matrix( v.cross(-mai.axis), acos(v.dot(-mai.axis)) );
            Mat r2 = rotation_matrix( mai.ortho, numeric::conversions::radians(55.25) );
            Mat r3 = rotation_matrix( mai.axis , numeric::conversions::radians(90.0) );
            rot_pose(h,r3*r2*r1);
            trans_pose(h,mai.cen-h.residue(1).xyz(hatm));
          }



          Vec rot_cen  = mai.cen;
          Vec rot_axis = (h.xyz(AtomID(natm,1))-mai.cen).normalized();
          vector1<Vec> FXD1,FXD2,CA;
          for(Size i = 1; i <= pose.n_residue(); ++i) {
            Vec fxd1 = (pose.residue(i).xyz("CB")-pose.residue(i).xyz("CA")).normalized();
            Vec fxd2 = (pose.residue(i).xyz( "N")-pose.residue(i).xyz("CA")).normalized();
            fxd2 = (fxd2-fxd1.dot(fxd2)*fxd1).normalized();
            FXD1.push_back(fxd1);
            FXD2.push_back(fxd2);
            CA.push_back(pose.residue(i).xyz("CA"));
          }
          // h.dump_pdb("h3.pdb");



          for(Size ich1 = 1; ich1 <= CHI1.size(); ++ich1) {
            for(Size ich2 = 1; ich2 <= CHI2.size(); ++ich2) {
              if( fabs(CHI1[ich1]-188.3) > 3.0 ) continue;
              if( fabs(CHI2[ich2]-274.7) > 3.0 ) continue;
              h.set_chi(1,1,CHI1[ich1]);  h.set_chi(2,1,CHI2[ich2]);
              TR << "checking " << CHI1[ich1] << " " << CHI2[ich2] << std::endl;
              Vec const mov10 = (h.residue(1).xyz("CB")-h.residue(1).xyz("CA")).normalized();
              Vec const tmp  = (h.residue(1).xyz( "N")-h.residue(1).xyz("CA")).normalized();
              Vec const mov20 = (tmp-mov10.dot(tmp)*mov10).normalized();
              for(Size irsd = 1; irsd <= pose.n_residue(); ++irsd) {
                if(irsd==rsd1 || irsd==rsd2) continue;
                if(irsd != 59) continue;
                Stub symm_stub1;
                Vec symm_axis1, symm_cen1;
                {
                  // Vec fxd1 = (pose.residue(irsd).xyz("CB")-pose.residue(irsd).xyz("CA")).normalized();
                  // Vec fxd2 = (pose.residue(irsd).xyz( "N")-pose.residue(irsd).xyz("CA")).normalized();
                  // fxd2 = (fxd2-fxd1.dot(fxd2)*fxd1).normalized();
                  Vec mov1 = (h.residue(1).xyz("CB")-h.residue(1).xyz("CA")).normalized();
                  Vec mov2 = (h.residue(1).xyz( "N")-h.residue(1).xyz("CA")).normalized();
                  mov2 = (mov2-mov1.dot(mov2)*mov1).normalized();
                  Vec fxd1 = FXD1[irsd];
                  Vec fxd2 = FXD2[irsd];
                  // Vec mov1 = mov10;
                  // Vec mov2 = mov20;

                  // Vec m1=mov1, m2=mov2, f1=fxd1, f2=fxd2;
                  Mat fxd2mov = rotation_matrix(mov1.cross(        fxd1),-acos(mov1.dot(        fxd1)));
                  fxd2mov     = rotation_matrix(mov2.cross(fxd2mov*fxd2),-acos(mov2.dot(fxd2mov*fxd2))) * fxd2mov;
                  Mat alignmov = rotation_matrix( rot_axis.cross(mov1), -acos(rot_axis.dot(mov1)) );
                  fxd1 = fxd2mov.transposed()*alignmov*fxd2mov*fxd1;
                  fxd2 = fxd2mov.transposed()*alignmov*fxd2mov*fxd2;
                  mov1 = alignmov*mov1;
                  mov2 = alignmov*mov2;
                  symm_axis1 = (fxd1.normalized()+mov1.normalized()).normalized();
                  Vec tgt = rotation_matrix( fxd2.cross(symm_axis1), 2.0*acos(fxd2.dot(symm_axis1)) ) * fxd2;
                  tgt      = (tgt - tgt.dot(rot_axis)*rot_axis).normalized();
                  Vec tgt0 = (mov2-mov2.dot(rot_axis)*rot_axis).normalized();
                  Real ang = acos(tgt.dot(tgt0));
                  if( tgt0.cross(tgt).dot(rot_axis) < 0 ) ang = 2.0*pi-ang;

                  rot_pose(h,rotation_matrix(rot_axis,ang),rot_cen);
                  symm_cen1 = (h.residue(1).xyz("CA")+CA[irsd])/2.0;
                  symm_stub1.M = rotation_matrix(rot_axis,ang)*fxd2mov;
                  symm_stub1.v = h.residue(1).xyz("CA") - symm_stub1.M*CA[irsd];
                  if( (symm_stub1.v.dot(symm_axis1)*symm_axis1).length() > 0.5 ) {
                    TR << "symm1 axis fail" << std::endl;
                    continue;
                  }
                  symm_stub1.v = symm_stub1.v - symm_stub1.v.dot(symm_axis1)*symm_axis1;
                }

                bool clash = false;
                for(Size i = 1; i <= h.residue(1).nheavyatoms(); ++i) {
                  for(Size j = 1; j <= 5; ++j) {
                    if( pose.xyz(AtomID(i,rsd1)).distance_squared(h.xyz(AtomID(j,1))) < HIS_CLASH_DIS2 ) clash = true;
                    if( pose.xyz(AtomID(i,rsd2)).distance_squared(h.xyz(AtomID(j,1))) < HIS_CLASH_DIS2 ) clash = true;
                    if(clash) break;
                  }
                  if(i > 5) if( !ms.clash_check(pose.xyz(AtomID(i,rsd1)),symm_stub1) ) clash = true;
                  if(i > 5) if( !ms.clash_check(pose.xyz(AtomID(i,rsd2)),symm_stub1) ) clash = true;
                  if(i > 5) if( !ms.clash_check(   h.xyz(AtomID(i,   1))           ) ) clash = true;
                  if(clash) break;
                }
                if(clash) {
                  TR << "symm1 his3 clash" << std::endl;
                  continue;
                }
                for(Size i = 5; i <= h.residue(1).nheavyatoms(); ++i) {
                  for(Size j = 5; j <= h.residue(1).nheavyatoms(); ++j) {
                    if( pose.xyz(AtomID(i,rsd1)).distance_squared( symm_stub1.local2global(pose.xyz(AtomID(j,rsd1))) ) < 1.0*HIS_CLASH_DIS2 ) clash = true;
                    if( pose.xyz(AtomID(i,rsd1)).distance_squared( symm_stub1.local2global(pose.xyz(AtomID(j,rsd2))) ) < 1.0*HIS_CLASH_DIS2 ) clash = true;
                    if( pose.xyz(AtomID(i,rsd2)).distance_squared( symm_stub1.local2global(pose.xyz(AtomID(j,rsd1))) ) < 1.0*HIS_CLASH_DIS2 ) clash = true;
                    if( pose.xyz(AtomID(i,rsd2)).distance_squared( symm_stub1.local2global(pose.xyz(AtomID(j,rsd2))) ) < 1.0*HIS_CLASH_DIS2 ) clash = true;
                    // if( symm_stub1.global2local(h.xyz(AtomID(i,1))).distance_squared( symm_stub1.local2global(pose.xyz(AtomID(j,rsd1))) ) < 1.0*HIS_CLASH_DIS2 ) clash = true;
                    // if( symm_stub1.global2local(h.xyz(AtomID(i,1))).distance_squared( symm_stub1.local2global(pose.xyz(AtomID(j,rsd2))) ) < 1.0*HIS_CLASH_DIS2 ) clash = true;
                    // if( symm_stub1.global2local(h.xyz(AtomID(i,1))).distance_squared(                            h.xyz(AtomID(j,   1))  ) < 1.0*HIS_CLASH_DIS2 ) clash = true;
                    if( h.xyz(AtomID(j,1)).distance_squared( mai.cen ) < 8.0 ) continue;
                    if( pose.xyz(AtomID(i,rsd1)).distance_squared(                            h.xyz(AtomID(j,   1))  ) < 1.0*HIS_CLASH_DIS2 ) clash = true;
                    if( pose.xyz(AtomID(i,rsd2)).distance_squared(                            h.xyz(AtomID(j,   1))  ) < 1.0*HIS_CLASH_DIS2 ) clash = true;
                  }
                  if(clash) break;
                }
                if(clash) {
                  TR << "symm1 symmetric his clash fail" << std::endl;
                  continue;
                }
                if( !ms.clash_check(symm_stub1) ) {
                  TR << "base clash fail" << std::endl;
                  continue;
                }

                Pose self = pose;
                replace_pose_residue_copying_existing_coordinates(self,irsd,h.residue_type(1));
                self.set_chi(1,irsd,CHI1[ich1]); self.set_chi(2,irsd,CHI2[ich2]);
                // Pose other = self;
                // replace_pose_residue_copying_existing_coordinates(other,irsd,h.residue_type(1));
                // other.set_chi(1,irsd,CHI1[ich1]); other.set_chi(2,irsd,CHI2[ich2]);
                // xform_pose(other,symm_stub1);
                // self.dump_pdb("self.pdb");
                // other.dump_pdb("other.pdb");
                // h.dump_pdb("h.pdb");
                // utility_exit_with_message("DEBUG");


                symm_cen1 = (self.residue(1).xyz("CA")+symm_stub1.local2global(self.residue(1).xyz("CA")))/2.0;

                if(mai.sqp) {
                  // TR << "SQP" << std::endl;
                  Vec v = (h.residue(1).xyz(hatm)-h.residue(1).xyz(natm)).normalized();
                  Mat r1 = rotation_matrix( v.cross(-mai.axis), acos(v.dot(-mai.axis)) );
                  Mat r2 = rotation_matrix( mai.ortho, numeric::conversions::radians(-45.0) );
                  rot_pose(h,r2*r1);
                  trans_pose(h,mai.cen-h.residue(1).xyz(hatm));
                } else {
                  // TR << "TET" << std::endl;
                  Vec v = (h.residue(1).xyz(hatm)-h.residue(1).xyz(natm)).normalized();
                  Mat r1 = rotation_matrix( v.cross(-mai.axis), acos(v.dot(-mai.axis)) );
                  Mat r2 = rotation_matrix( mai.ortho, numeric::conversions::radians(55.25) );
                  Mat r3 = rotation_matrix( mai.axis , numeric::conversions::radians(-90.0) );
                  rot_pose(h,r3*r2*r1);
                  trans_pose(h,mai.cen-h.residue(1).xyz(hatm));
                }
                Vec rot_cen  = mai.cen;
                Vec rot_axis = (h.xyz(AtomID(natm,1))-mai.cen).normalized();

                TR << "MATCH partner 1" << std::endl;
                // rot_pose(self,rotation_matrix(symm_axis1.cross(Vec(1,0,0)),-acos(symm_axis1.dot(Vec(1,0,0)))));
                // rot_pose(other,rotation_matrix(symm_axis1.cross(Vec(1,0,0)),-acos(symm_axis1.dot(Vec(1,0,0)))));
                // rot_pose(h,rotation_matrix(symm_axis1.cross(Vec(1,0,0)),-acos(symm_axis1.dot(Vec(1,0,0)))));

                // Real ymin=9e9,ymax=-9e9,zmin=9e9,zmax=-9e9,rmin=9e9,rmax=-9e9;
                clash = false;
                for(Real jrot = 0; jrot < 360; jrot += 5.0) {
                  TR << "JROT " << jrot << std::endl;
                  rot_pose( h, rotation_matrix_degrees(rot_axis,5.0), rot_cen );
                  for(Size jch1 = 1; jch1 <= CHI1.size(); ++jch1) {
                    for(Size jch2 = 1; jch2 <= CHI2.size(); ++jch2) {
                      h.set_chi(1,1,CHI1[jch1]);  h.set_chi(2,1,CHI2[jch2]);
                      for(Size i = 1; i <= 7; ++i) {
                        for(Size j = 7; j <= h.residue(1).nheavyatoms(); ++j) {
                          if( self.xyz(AtomID(j,rsd1)).distance_squared( h.xyz(AtomID(i,1)) ) < HIS_CLASH_DIS2 ) clash = true;
                          if( self.xyz(AtomID(j,rsd2)).distance_squared( h.xyz(AtomID(i,1)) ) < HIS_CLASH_DIS2 ) clash = true;
                          if( self.xyz(AtomID(j,irsd)).distance_squared( h.xyz(AtomID(i,1)) ) < HIS_CLASH_DIS2 ) clash = true;
                        }
                        if(clash) break;
                        if( !ms.clash_check( h.xyz(AtomID(i,1)) ) ) clash = true;
                        if( !ms.clash_check( symm_stub1.local2global(h.xyz(AtomID(i,1))) ) ) clash = true;
                        if(clash) break;
                      }
                      if(clash) continue;

                      Mat rot = rotation_matrix(symm_axis1.cross(Vec(1,0,0)),acos(symm_axis1.dot(Vec(1,0,0))));
                      Stub align( rot, -(rot*symm_cen1) );
                      xform_pose(h,align);
                      xform_pose(self,align);

                      Vec const ca = h.residue(1).xyz(2);
                      Real const r2 = ca.y()*ca.y()+ca.z()*ca.z();
                      Vec hcacb = (h.xyz(AtomID(2,1))-h.xyz(AtomID(5,1))).normalized(); // 6 = CB?? removed termini... 5 should be ok
                      for(Size jrsd = 1; jrsd <= self.n_residue(); ++jrsd) {
                        if( jrsd==irsd || jrsd==rsd1 || jrsd==rsd2 ) continue;
                        Vec ca2 = self.residue(jrsd).xyz(2);
                        Real const r2j = ca2.y()*ca2.y()+ca2.z()*ca2.z();
                        if( fabs(r2j-r2) > 20.0 ) continue;
                        if( (sqrt(r2)-sqrt(r2j)) > 0.9 ) continue;

                        Mat h4rot = rotation_matrix_degrees(Vec(0,1,0),180.0);

                        Vec hca =          h.residue(   1).xyz("CA"); hca.x() = 0.0;
                        Vec jca = h4rot*self.residue(jrsd).xyz("CA"); jca.x() = 0.0;
                        Real hjang = acos(hca.normalized().dot(jca.normalized()));
                        Mat hjrot = rotation_matrix(hca.cross(jca),-hjang);

                        Vec jcacb = (self.xyz(AtomID(2,jrsd))-self.xyz(AtomID(5,jrsd))).normalized();
                        if( (hjrot*h4rot*jcacb).dot(hcacb) < 0.90 ) continue;
                        // TR << "checking jrsd " << jrsd << std::endl;
                        Vec dca = ca - hjrot*h4rot*self.xyz(AtomID(2,jrsd));
                        Stub symm_stub2( hjrot*h4rot, Vec(dca.x(),0,0) );
                        if( !ms.clash_check(symm_stub2) ) continue;


                        Pose self2 = self;
                        replace_pose_residue_copying_existing_coordinates(self2,irsd,h.residue_type(1));
                        self2.set_chi(1,irsd,CHI1[ich1]); self2.set_chi(2,irsd,CHI2[ich2]);
                        replace_pose_residue_copying_existing_coordinates(self2,jrsd,h.residue_type(1));
                        Real best = 9e9, bestch1=0, bestch2=0;
                        for(Real kch1 = 0.0; kch1 < 360.0; kch1 += 2.0 ) {
                          self2.set_chi(1,jrsd,kch1);
                          for(Real kch2 = 0.0; kch2 < 360.0; kch2 += 2.0 ) {
                            self2.set_chi(2,jrsd,kch2);
                            if( h.xyz(AtomID(natm,1)).distance_squared(self2.xyz(AtomID(natm,jrsd))) < best ) {
                              best = h.xyz(AtomID(natm,1)).distance_squared(self2.xyz(AtomID(natm,jrsd)));
                              bestch1 = kch1;
                              bestch2 = kch2;
                            }
                          }
                        }
                        if(best > 0.25) continue;
                        self2.set_chi(1,jrsd,bestch1); self2.set_chi(2,jrsd,bestch2);

                        Pose other = self2;
                        xform_pose_rev(other,align);
                        xform_pose(other,symm_stub1);
                        xform_pose(other,align);

                        Pose other2 = self2;
                        // xform_pose(other2,align); // = self, already aligned to x axis
                        xform_pose(other2,symm_stub2);

                        Pose other3 = other2;
                        rot_pose(other3,x_rotation_matrix_degrees(180.0));


                        TR << jrsd << " " << hjang*180.0/pi << " " << jca << " " << hca << std::endl;
                        string tag = string_of(uniform());
                        self.dump_pdb(tag+"self.pdb");
                        other.dump_pdb(tag+"other.pdb");
                        other2.dump_pdb(tag+"other2.pdb");
                        other3.dump_pdb(tag+"other3.pdb");
                        h.dump_pdb(tag+"h4.pdb");
                        h.dump_pdb(tag+"h4B.pdb");
                        utility_exit_with_message("DEBUG");

                      }
                      if( (symm_axis1-Vec(1,0,0)).length() > 0.00001 ) {
                        Mat rot = rotation_matrix(symm_axis1.cross(Vec(1,0,0)),acos(symm_axis1.dot(Vec(1,0,0))));
                        Stub align( rot, -(rot*symm_cen1) );
                        xform_pose_rev(h,align);
                        xform_pose_rev(self,align);
                      }
                    }
                  }
                }



                // Mat rot = rotation_matrix( rot_axis, ang );
                // m1 = rot*m1;
                // m2 = rot*m2;
                // ozstream out2("rot.pdb");
                // out2 << "HETATM" << I(5,1) << ' ' << "XXXX" << ' ' <<  "FX1" << ' ' << "A" << I(4,1) << "    " << F(8,3, f1.x()*10.) << F(8,3, f1.y()*10.) << F(8,3, f1.z()*10.) << F(6,2,1.0) << F(6,2,1.0) << std::endl;
                // out2 << "HETATM" << I(5,2) << ' ' << "XXXX" << ' ' <<  "SYM" << ' ' << "A" << I(4,2) << "    " << F(8,3, symm_axis1.x()*10.) << F(8,3, symm_axis1.y()*10.) << F(8,3, symm_axis1.z()*10.) << F(6,2,1.0) << F(6,2,1.0) << std::endl;
                // out2 << "HETATM" << I(5,3) << ' ' << "XXXX" << ' ' <<  "FX2" << ' ' << "A" << I(4,3) << "    " << F(8,3, f2.x()*10.) << F(8,3, f2.y()*10.) << F(8,3, f2.z()*10.) << F(6,2,1.0) << F(6,2,1.0) << std::endl;
                // out2 << "HETATM" << I(5,4) << ' ' << "XXXX" << ' ' <<  "MV1" << ' ' << "A" << I(4,4) << "    " << F(8,3, m1.x()*10.) << F(8,3, m1.y()*10.) << F(8,3, m1.z()*10.) << F(6,2,1.0) << F(6,2,1.0) << std::endl;
                // out2 << "HETATM" << I(5,5) << ' ' << "XXXX" << ' ' <<  "MV2" << ' ' << "A" << I(4,5) << "    " << F(8,3, m2.x()*10.) << F(8,3, m2.y()*10.) << F(8,3, m2.z()*10.) << F(6,2,1.0) << F(6,2,1.0) << std::endl;
                // out2 << "HETATM" << I(5,6) << ' ' << "XXXX" << ' ' <<  "ORI" << ' ' << "A" << I(4,6) << "    " << F(8,3,         0.0) << F(8,3,         0.0) << F(8,3,         0.0) << F(6,2,1.0) << F(6,2,1.0) << std::endl;
                // out2 << "HETATM" << I(5,6) << ' ' << "XXXX" << ' ' <<  "ROT" << ' ' << "A" << I(4,6) << "    " << F(8,3,rot_axis.x()*10) << F(8,3,rot_axis.y()*10) << F(8,3,rot_axis.z()*10) << F(6,2,1.0) << F(6,2,1.0) << std::endl;
                // out2.close();

              }
            }
          }
        }
      }
    }
  if(!basic::options::option[basic::options::OptionKeys::willmatch::symmetry_d2]()) {
    vector1< Pose > poses1,poses2;
    Pose native1,native2;
    core::import_pose::pose_from_pdb(poses1,option[in::file::s][1]);
    core::import_pose::pose_from_pdb(native1,option[willmatch::native1]());
    vector1<Size> allowed_res1,allowed_res2;
    if( option[crossmatch::allowed_res1].user() ) allowed_res1 = option[crossmatch::allowed_res1]();
    if( option[crossmatch::allowed_res2].user() ) allowed_res2 = option[crossmatch::allowed_res2]();

    vector1<Size> iface_candidates1,iface_candidates2;
    vector1<Size> tmp = read_res_list(basic::options::option[basic::options::OptionKeys::willmatch::exclude_res1]());
    for(Size i = 1; i <= native1.n_residue(); ++i) if(std::find(tmp.begin(),tmp.end(),i)==tmp.end()) iface_candidates1.push_back(i);
    MatchSet ms1(poses1,mop,native1,allowed_res1);
    ms1.iface_candidates = iface_candidates1;
    MatchSet ms2;
    if(option[in::file::s].size()==1) {
      poses2 = poses1;
      native2 = native1;
      ms2 = ms1;
    } else {
      core::import_pose::pose_from_pdb(poses2,option[in::file::s][2]);
      core::import_pose::pose_from_pdb(native2,option[willmatch::native2]());
      ms2.init(poses2,mop,native2,allowed_res2);
      tmp = read_res_list(basic::options::option[basic::options::OptionKeys::willmatch::exclude_res2]());
      for(Size i = 1; i <= native2.n_residue(); ++i) if(std::find(tmp.begin(),tmp.end(),i)==tmp.end()) iface_candidates2.push_back(i);
      ms2.iface_candidates = iface_candidates2;
    }
    Filter *filter;
    std::pair<vector1<Size>,vector1<Size> > splitwork = makesplitwork(ms1.ligs.size(),ms2.ligs.size());
    vector1<Size> ILIG1 = splitwork.first;
    vector1<Size> ILIG2 = splitwork.second;
    assert(ILIG2.size()==ILIG1.size());
    if(basic::options::option[basic::options::OptionKeys::willmatch::taglist].user()) {
      vector1<Tags> tagslist;
      utility::io::izstream in(basic::options::option[basic::options::OptionKeys::willmatch::taglist]());
      string s1,s2;
      Size rot;
      while(in >> s1 >> s2 >> rot) {
        tagslist.push_back(Tags(s1,s2,rot));
      }
      filter = new ListFilter(tagslist);
    } else {
      filter = new NullFilter;
    }
    for(Size iligidx = 1; iligidx <= ILIG1.size(); ++iligidx) {
      Size ilig = ILIG1[iligidx];
      Size jlig = ILIG2[iligidx];
      // TR << "SPLITWORK " << ilig << "-" << jlig << std::endl;
      // continue;
      // // just to speed things up?
      // std::cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
      if( basic::options::option[basic::options::OptionKeys::willmatch::homodimer]() && ilig!=jlig) continue;
      if(iligidx%100==0) TR << "cross " << ilig << " " << jlig << std::endl;
      for(Size i = 1; i <= 2; ++i) {
        try {
          // TR << "try cross " << ilig << " " << jlig << " " << i << std::endl;
          if( basic::options::option[basic::options::OptionKeys::willmatch::homodimer]()) {
            ms1.cross_homodimer(ilig,i);
          } else {
            ms1.cross(ms2,ilig,jlig,i,*filter);
          }
        } catch(...) {
          TR << "EXCEPTION " << std::endl;
          continue;
        }
      }
    }
    delete filter;
    // ms1.write_to_file("test.pdb");
  }

  return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}










