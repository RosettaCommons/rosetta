// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/apps/pilat/will/genmatch.cc
/// @brief ???

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
//#include <basic/options/keys/willmatch.OptionKeys.gen.hh>
#include <basic/options/option.hh>
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
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/scoring/ImplicitFastClashCheck.hh>
#include <sstream>
#include <utility/io/ozstream.hh>
#include <utility/string_constants.hh>
// #include <devel/init.hh>

// #include <core/scoring/constraints/LocalCoordinateConstraint.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/pose/util.tmpl.hh>
#include <apps/pilot/will/mynamespaces.ihh>
#include <apps/pilot/will/will_util.ihh>


using core::kinematics::Stub;
using protocols::scoring::ImplicitFastClashCheck;
using core::pose::Pose;
using core::conformation::ResidueOP;

static THREAD_LOCAL basic::Tracer TR( "genI213" );
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

vector1<Size> get_scanres(Pose const & pose) {
  vector1<Size> scanres;
  //if(basic::options::option[basic::options::OptionKeys::willmatch::residues].user()) {
	//TR << "input scanres!!!!!!" << std::endl;
	//scanres = basic::options::option[basic::options::OptionKeys::willmatch::residues]();
	//} else {
    for(Size i = 1; i <= pose.n_residue(); ++i) {
      if(!pose.residue(i).has("N" )) { continue; }
      if(!pose.residue(i).has("CA")) { continue; }
      if(!pose.residue(i).has("C" )) { continue; }
      if(!pose.residue(i).has("O" )) { continue; }
      if(!pose.residue(i).has("CB")) { continue; }
      if(pose.residue(i).name3()=="PRO") { continue; }
      scanres.push_back(i);
    }
		//}
  return scanres;
}


void dumpsym(Pose const & pose, Mat R2, Mat R3a, Mat R3b, Vec cen2, string fname) {
  vector1<Vec> seenit;
  Mat R3[3];
  R3[0] = Mat::identity();
  R3[1] = R3a;
  R3[2] = R3b;
  TR << "output" << std::endl;
	vector1<string> ANAME(3);
	ANAME[1] = " N  ";
	ANAME[2] = " CA ";
	ANAME[3] = " C  ";
	string const & CHAIN( utility::UPPERCASE_LETTERS );
  ozstream out( fname );
  Size acount=0,rcount=0,ccount=0;
  for(Size i3a = 0; i3a < 3; i3a++) {
    for(Size i2a = 0; i2a < 2; i2a++) {
      for(Size j3a = 0; j3a < 3; j3a++) {
        for(Size j2a = 0; j2a < 2; j2a++) {
          for(Size k3a = 0; k3a < 3; k3a++) {
            for(Size k2a = 0; k2a < 2; k2a++) {
              for(Size l3a = 0; l3a < 3; l3a++) {
                for(Size l2a = 0; l2a < 2; l2a++) {
                  for(Size m3a = 0; m3a < 2; m3a++) {
                    for(Size m2a = 0; m2a < 2; m2a++) {
                      for(Size n3a = 0; n3a < 2; n3a++) {
                        for(Size n2a = 0; n2a < 2; n2a++) {

                          Vec chk( pose.xyz(AtomID(1,1)) );
													chk = R3[i3a]*chk; if(i2a) chk = R2*(chk-cen2)+cen2;
													chk = R3[j3a]*chk; if(j2a) chk = R2*(chk-cen2)+cen2;
													chk = R3[k3a]*chk; if(k2a) chk = R2*(chk-cen2)+cen2;
													chk = R3[l3a]*chk; if(l2a) chk = R2*(chk-cen2)+cen2;
													chk = R3[m3a]*chk; if(m2a) chk = R2*(chk-cen2)+cen2;
													chk = R3[n3a]*chk; if(n2a) chk = R2*(chk-cen2)+cen2;
                          for(vector1<Vec>::const_iterator i = seenit.begin(); i != seenit.end(); ++i) {
                            if( i->distance_squared(chk) < 1.0 ) goto cont2;
                          }
                          goto done2; cont2: continue; done2:
                          seenit.push_back(chk);

													char chain = CHAIN[ccount];
													if( (i2a+j2a+k2a+l2a+m2a+n2a) % 2 == 1 ) chain = 'B';
                          for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
														if( rcount >= 9999) {
															rcount = 0;
															ccount++;
														}
                            Size rn = ++rcount;
														for(Size ia = 1; ia <= 3; ia++) {
															Vec tmp(pose.residue(ir).xyz(ia));
															tmp = R3[i3a]*tmp; if(i2a) tmp = R2*(tmp-cen2)+cen2;
															tmp = R3[j3a]*tmp; if(j2a) tmp = R2*(tmp-cen2)+cen2;
															tmp = R3[k3a]*tmp; if(k2a) tmp = R2*(tmp-cen2)+cen2;
															tmp = R3[l3a]*tmp; if(l2a) tmp = R2*(tmp-cen2)+cen2;
															tmp = R3[m3a]*tmp; if(m2a) tmp = R2*(tmp-cen2)+cen2;
															tmp = R3[n3a]*tmp; if(n2a) tmp = R2*(tmp-cen2)+cen2;
															string X = F(8,3,tmp.x());
															string Y = F(8,3,tmp.y());
															string Z = F(8,3,tmp.z());
															out<<"ATOM  "<<I(5,++acount)<<' '<<ANAME[ia]<<' '<<"ALA"<<' '<<chain<<I(4,rn)<<"    "<<X<<Y<<Z<<F(6,2,1.0)<<F(6,2,0.0)<<'\n';
														}
                          }
													out << "TER" << std::endl;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  out.close();
}

void dumpsym6(Pose const & pose, Mat R2, Mat R3a, Mat R3b, Vec HG, string fname, Size irsd, Size ch1, Size ch2, Real ANG) {
  vector1<Vec> seenit;
  Mat R3[3];
  R3[0] = Mat::identity();
  R3[1] = R3a;
  R3[2] = R3b;
  TR << "output" << std::endl;
  for(Size i3a = 0; i3a < 3; i3a++) {
    for(Size i2a = 0; i2a < 2; i2a++) {
      for(Size j3a = 0; j3a < 3; j3a++) {
        for(Size j2a = 0; j2a < 2; j2a++) {
          for(Size k3a = 0; k3a < 3; k3a++) {
            for(Size k2a = 0; k2a < 2; k2a++) {
              for(Size l3a = 0; l3a < 3; l3a++) {
                for(Size l2a = 0; l2a < 2; l2a++) {
                  for(Size m3a = 0; m3a < 2; m3a++) {
                    for(Size m2a = 0; m2a < 2; m2a++) {
                      for(Size n3a = 0; n3a < 2; n3a++) {
                        for(Size n2a = 0; n2a < 2; n2a++) {
                          Pose tmp(pose);
                          rot_pose(tmp,R3[i3a]); if(i2a) rot_pose(tmp,R2,HG);
                          rot_pose(tmp,R3[j3a]); if(j2a) rot_pose(tmp,R2,HG);
                          rot_pose(tmp,R3[k3a]); if(k2a) rot_pose(tmp,R2,HG);
                          rot_pose(tmp,R3[l3a]); if(l2a) rot_pose(tmp,R2,HG);
                          rot_pose(tmp,R3[m3a]); if(m2a) rot_pose(tmp,R2,HG);
                          rot_pose(tmp,R3[n3a]); if(n2a) rot_pose(tmp,R2,HG);
                          Vec chk( tmp.xyz(AtomID(1,1)) );
                          for(vector1<Vec>::const_iterator i = seenit.begin(); i != seenit.end(); ++i) {
                            if( i->distance_squared(chk) < 1.0 ) goto cont2;
                          }
                          goto done2; cont2: continue; done2:
                          seenit.push_back(chk);
                          tmp.dump_pdb( utility::file_basename(fname)
                                        +"_"+  F(4,1,ANG)+"_"+ lzs(irsd,3)+"_"+ lzs((Size)ch1,3)+"_"+ lzs((Size)ch2,3)+"_"+
                                        lzs(i3a,1)+"_"+lzs(i2a,1)+"_"+lzs(j3a,1)+"_"+lzs(j2a,1)+"_"+lzs(k3a,1)+"_"+lzs(k2a,1)+"_"+
                                        lzs(l3a,1)+"_"+lzs(l2a,1)+"_"+lzs(m3a,1)+"_"+lzs(m2a,1)+"_"+lzs(n3a,1)+"_"+lzs(n2a,1)+".pdb");
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

void dumpsym2(Pose const & pose, Mat R2, Mat R3a, Mat R3b, Vec HG, string fname, Size irsd, Size ch1, Size ch2, Real ANG) {
  vector1<Vec> seenit;
  Mat R3[3];
  R3[0] = Mat::identity();
  R3[1] = R3a;
  R3[2] = R3b;
  TR << "output" << std::endl;
  for(Size i3a = 0; i3a < 3; i3a++) {
    for(Size i2a = 0; i2a < 2; i2a++) {
      for(Size j3a = 0; j3a < 3; j3a++) {
        for(Size j2a = 0; j2a < 2; j2a++) {
          Pose tmp(pose);
          rot_pose(tmp,R3[i3a]); if(i2a) rot_pose(tmp,R2,HG);
          rot_pose(tmp,R3[j3a]); if(j2a) rot_pose(tmp,R2,HG);
          Vec chk( tmp.xyz(AtomID(1,1)) );
          for(vector1<Vec>::const_iterator i = seenit.begin(); i != seenit.end(); ++i) {
            if( i->distance_squared(chk) < 1.0 ) goto cont2;
          }
          goto done2; cont2: continue; done2:
          seenit.push_back(chk);
          tmp.dump_pdb( utility::file_basename(fname)
                        +"_"+  F(4,1,ANG)+"_"+ lzs(irsd,3)+"_"+ lzs((Size)ch1,3)+"_"+ lzs((Size)ch2,3)+"_"+
                        lzs(i3a,1)+"_"+lzs(i2a,1)+"_"+lzs(j3a,1)+"_"+lzs(j2a,1)+".pdb");
        }
      }
    }
  }
}

struct Hit {
  Size rsd,cbc;
  Real chi1,chi2;
  Vec axs,cen;
  bool sym;
};


#define ATET 54.735610317245360079 // asin(sr2/sr3)
#define AOCT 35.264389682754668343 // asin(sr1/sr3)


std::pair<vector1<Hit>,vector1<Hit> >
dock(Pose & init, string fname) {

  using namespace basic::options::OptionKeys;

  core::chemical::ResidueTypeSetCAP  rs = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
  Pose cys;

  make_pose_from_sequence(cys,"C",core::chemical::FA_STANDARD,false);
  remove_lower_terminus_type_from_pose_residue(cys,1);
  remove_upper_terminus_type_from_pose_residue(cys,1);
  //  add_variant_type_to_pose_residue(cys,"DISULF_PARTNER",1);

  core::scoring::dssp::Dssp dssp(init);
  dssp.insert_ss_into_pose(init);

  Vec com(0,0,0);
  for(Size ir = 1; ir <= init.n_residue(); ++ir) {
    init.replace_residue(ir,cys.residue(1),true);
    replace_pose_residue_copying_existing_coordinates(init,ir,init.residue(ir).residue_type_set().name_map("ALA"));
    com += init.xyz(AtomID(2,ir));
  }
  com /= init.n_residue();

  protocols::scoring::ImplicitFastClashCheck ifc3(init,3.5);
  protocols::scoring::ImplicitFastClashCheck ifc2(init,2.8);

  Size nres = init.n_residue();
  // ScoreFunctionOP sf = core::scoring::get_score_function();
  ScoreFunctionOP sf = new core::scoring::symmetry::SymmetricScoreFunction(core::scoring::get_score_function());

  Pose pose = init;

  Mat R3a = rotation_matrix_degrees(Vec(0,0,1), 120.0);
  Mat R3b = rotation_matrix_degrees(Vec(0,0,1),-120.0);

  vector1<Size> scanres = get_scanres(pose);

  core::pack::rotamers::SingleResidueRotamerLibraryCAP dunlib = core::pack::dunbrack::RotamerLibrary::get_instance()->get_rsd_library( rs->name_map("CYS") );
  core::pack::dunbrack::RotamerLibraryScratchSpace scratch;

  vector1<Hit> hits,allhits;

  for(vector1<Size>::const_iterator iiter = scanres.begin(); iiter != scanres.end(); ++iiter) {
    Size irsd = *iiter;
    // ssq
    for( int iss = int(irsd)-2; iss <= int(irsd)+2; ++iss) {
      if( iss < 1 || iss > pose.n_residue() ) goto contss;
      if( pose.secstruct(iss) != 'L' ) goto doness;
    } goto doness; contss: continue; doness:

		//    if(irsd != 42) continue;

    //TR << fname << " RES " << irsd << std::endl;
    //if( pose.residue(irsd).name3()=="GLY" || pose.residue(irsd).name3()=="PRO" ) continue;
    ResidueOP rprev = pose.residue(irsd).clone();
    pose.replace_residue(irsd,cys.residue(1),true);
    Vec CA = pose.residue(irsd).xyz("CA");
    Vec CB = pose.residue(irsd).xyz("CB");
    vector1<Vec> seen_axs,seen_cen;
    vector1<Vec> asym_axs,asym_cen;

    bool foundhit = false;
    for(Size stage=0; stage<=1; ++stage) {

      for(Real ch1 = 0.0; ch1 < 360.0; ch1+=2.0) {
        pose.set_chi(1,irsd,ch1);
        Vec SG = pose.residue(irsd).xyz("SG");
        if( ! ifc2.clash_check(SG,irsd) ) continue;

        for(Real ch2 = 0.0; ch2 < 360.0; ch2+=3.0) {
          pose.set_chi(2,irsd,ch2);

          Real dun = dunlib->rotamer_energy( pose.residue(irsd), scratch );
          //if( dun > 4.0 ) continue;                                                // DUNBRACK

          Vec HG = pose.residue(irsd).xyz("HG");
          if( ! ifc2.clash_check(HG,irsd) ) continue;

          for(Size iaxs = 0; iaxs <= 1; iaxs++) {
            for(Size isg = 1; isg <= 2; isg++) {
              Real const ANG(  isg==1 ? 54.7 : 35.3  );

              Vec axs = rotation_matrix_degrees(HG-SG, iaxs?45.0:-45.0 ) * projperp(HG-SG,SG-CB).normalized();
              Real a1 = fabs(angle_degrees( axs,Vec(0,0,0), Vec(0,0,1) ) - ANG);
              Real a2 = fabs(angle_degrees( axs,Vec(0,0,0),-Vec(0,0,1) ) - ANG);
              Vec axs1;

              if(stage==1 && foundhit) { ////////////////////////////////////////////////////////////////////////////////

                axs1 = axs;
                Mat R2 = rotation_matrix_degrees(axs1,180.0);
                if( ! ifc2.clash_check( R2*(SG-HG)+HG ,irsd) ) continue;

                Vec cen = R2*(com-HG)+HG; // filter redundency at ~0.7Á and 5deg
                for(Size i = 1; i <= asym_cen.size(); ++i) {
                  if( asym_cen[i].distance_squared(cen) < 0.25 && asym_axs[i].dot(axs1) > 0.9962 ) goto cont0;
                } goto done0; cont0: continue; done0:

                for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
                  for(Size ia = 1; ia <= 5; ++ia) {
                    if( ! ifc3.clash_check( ( R2*((pose.xyz(AtomID(ia,ir)))-HG) + HG ) ) ) goto cont1;
                  }
                } goto done1; cont1: continue; done1:

                Size cbc = 0;
                for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
                  Vec CB12 = R2 * (       (pose.xyz(AtomID(5,ir)))-HG) + HG;
                  Vec CB22 = R2 * ( R3a * (pose.xyz(AtomID(5,ir)))-HG) + HG;
                  Vec CB32 = R2 * ( R3b * (pose.xyz(AtomID(5,ir)))-HG) + HG;
                  for(Size jr = 1; jr <= pose.n_residue(); ++jr) {
                    Vec CB11 =       pose.xyz(AtomID(5,jr));
                    Vec CB21 = R3a * pose.xyz(AtomID(5,jr));
                    Vec CB31 = R3b * pose.xyz(AtomID(5,jr));
                    if( CB11.distance_squared(CB12) < 49.0 ) cbc++;
                    if( CB21.distance_squared(CB12) < 49.0 ) cbc++;
                    if( CB31.distance_squared(CB12) < 49.0 ) cbc++;
                    if( CB11.distance_squared(CB22) < 49.0 ) cbc++;
                    if( CB21.distance_squared(CB22) < 49.0 ) cbc++;
                    if( CB31.distance_squared(CB22) < 49.0 ) cbc++;
                    if( CB11.distance_squared(CB32) < 49.0 ) cbc++;
                    if( CB21.distance_squared(CB32) < 49.0 ) cbc++;
                    if( CB31.distance_squared(CB32) < 49.0 ) cbc++;
                  }
                }
                if( cbc < 20 ) continue;
								//						TR << "4" << std::endl;
                for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
                  for(Size ia = 1; ia <= 5; ++ia) {
                    if( ! ifc3.clash_check(       ( R2 * (       (pose.xyz(AtomID(ia,ir)))-HG) + HG ) ) ) goto cont2;
                    if( ! ifc3.clash_check(       ( R2 * ( R3a * (pose.xyz(AtomID(ia,ir)))-HG) + HG ) ) ) goto cont2;
                    if( ! ifc3.clash_check(       ( R2 * ( R3b * (pose.xyz(AtomID(ia,ir)))-HG) + HG ) ) ) goto cont2;
                    if( ! ifc3.clash_check( R3a * ( R2 * (       (pose.xyz(AtomID(ia,ir)))-HG) + HG ) ) ) goto cont2;
                    if( ! ifc3.clash_check( R3a * ( R2 * ( R3a * (pose.xyz(AtomID(ia,ir)))-HG) + HG ) ) ) goto cont2;
                    if( ! ifc3.clash_check( R3a * ( R2 * ( R3b * (pose.xyz(AtomID(ia,ir)))-HG) + HG ) ) ) goto cont2;
                    if( ! ifc3.clash_check( R3b * ( R2 * (       (pose.xyz(AtomID(ia,ir)))-HG) + HG ) ) ) goto cont2;
                    if( ! ifc3.clash_check( R3b * ( R2 * ( R3a * (pose.xyz(AtomID(ia,ir)))-HG) + HG ) ) ) goto cont2;
                    if( ! ifc3.clash_check( R3b * ( R2 * ( R3b * (pose.xyz(AtomID(ia,ir)))-HG) + HG ) ) ) goto cont2;
                  }
                } goto done2; cont2: continue; done2:
								//TR << "5" << std::endl;
                asym_cen.push_back(R2*(com-HG)+HG);
                asym_axs.push_back(axs1);

                Hit h;
                h.rsd = irsd;
                h.cbc = cbc;
                h.axs = axs1;
                h.cen = HG;
                h.chi1 = ch1;
                h.chi2 = ch2;
                h.sym = false;
                allhits.push_back(h);

              }

							if(stage != 0) continue;
              if( dun > 2.0 ) continue;                                                // DUNBRACK
              if( a1 > 15.0 && a2 > 15.0 ) continue; ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

              if( a1 < a2 ) axs1 = rotation_matrix_degrees( Vec(0,0, 1).cross(axs), ANG ) * Vec(0,0, 1);
              else          axs1 = rotation_matrix_degrees( Vec(0,0,-1).cross(axs), ANG ) * Vec(0,0,-1);
              if( angle_degrees(axs1,Vec(0,0,0),axs) > 15.0 ) {
                TR << "VEC " << axs << " " << axs1 << " " << a1 << " " << a2 << "   " << angle_degrees(axs1,Vec(0,0,0),axs) << std::endl;
                utility_exit_with_message("axs alignment issue");
              }
              if(axs1.z() < 0.0) axs1 = -axs1;

              Mat R2 = rotation_matrix_degrees(axs1,180.0);
              if( ! ifc2.clash_check( R2*(SG-HG)+HG ,irsd) ) continue;

              Vec cen = R2*(com-HG)+HG; // filter redundency at 0.5Á and 5deg
              for(Size i = 1; i <= seen_cen.size(); ++i) {
                if( seen_cen[i].distance_squared(cen) < 0.25 && seen_axs[i].dot(axs1) > 0.9962 ) goto cont3;
              } goto done3; cont3: continue; done3:

              for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
                for(Size ia = 1; ia <= 5; ++ia) {
                  if( ! ifc3.clash_check( ( R2*((pose.xyz(AtomID(ia,ir)))-HG) + HG ) ) ) goto cont4;
                }
              } goto done4; cont4: continue; done4:

              Size cbc = 0;
              for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
                Vec CB12 = R2 * (       (pose.xyz(AtomID(5,ir)))-HG) + HG;
                Vec CB22 = R2 * ( R3a * (pose.xyz(AtomID(5,ir)))-HG) + HG;
                Vec CB32 = R2 * ( R3b * (pose.xyz(AtomID(5,ir)))-HG) + HG;
                for(Size jr = 1; jr <= pose.n_residue(); ++jr) {
                  Vec CB11 =       pose.xyz(AtomID(5,jr));
                  Vec CB21 = R3a * pose.xyz(AtomID(5,jr));
                  Vec CB31 = R3b * pose.xyz(AtomID(5,jr));
                  if( CB11.distance_squared(CB12) < 49.0 ) cbc++;
                  if( CB21.distance_squared(CB12) < 49.0 ) cbc++;
                  if( CB31.distance_squared(CB12) < 49.0 ) cbc++;
                  if( CB11.distance_squared(CB22) < 49.0 ) cbc++;
                  if( CB21.distance_squared(CB22) < 49.0 ) cbc++;
                  if( CB31.distance_squared(CB22) < 49.0 ) cbc++;
                  if( CB11.distance_squared(CB32) < 49.0 ) cbc++;
                  if( CB21.distance_squared(CB32) < 49.0 ) cbc++;
                  if( CB31.distance_squared(CB32) < 49.0 ) cbc++;
                }
              }
              if( cbc < 30 ) continue;

              for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
                for(Size ia = 1; ia <= 5; ++ia) {
                  if( ! ifc3.clash_check(       ( R2 * (       (pose.xyz(AtomID(ia,ir)))-HG) + HG ) ) ) goto cont5;
                  if( ! ifc3.clash_check(       ( R2 * ( R3a * (pose.xyz(AtomID(ia,ir)))-HG) + HG ) ) ) goto cont5;
                  if( ! ifc3.clash_check(       ( R2 * ( R3b * (pose.xyz(AtomID(ia,ir)))-HG) + HG ) ) ) goto cont5;
                  if( ! ifc3.clash_check( R3a * ( R2 * (       (pose.xyz(AtomID(ia,ir)))-HG) + HG ) ) ) goto cont5;
                  if( ! ifc3.clash_check( R3a * ( R2 * ( R3a * (pose.xyz(AtomID(ia,ir)))-HG) + HG ) ) ) goto cont5;
                  if( ! ifc3.clash_check( R3a * ( R2 * ( R3b * (pose.xyz(AtomID(ia,ir)))-HG) + HG ) ) ) goto cont5;
                  if( ! ifc3.clash_check( R3b * ( R2 * (       (pose.xyz(AtomID(ia,ir)))-HG) + HG ) ) ) goto cont5;
                  if( ! ifc3.clash_check( R3b * ( R2 * ( R3a * (pose.xyz(AtomID(ia,ir)))-HG) + HG ) ) ) goto cont5;
                  if( ! ifc3.clash_check( R3b * ( R2 * ( R3b * (pose.xyz(AtomID(ia,ir)))-HG) + HG ) ) ) goto cont5;
                }
              } goto done5; cont5: continue; done5:

              Mat R3[3];
              R3[0] = Mat::identity();
              R3[1] = R3a;
              R3[2] = R3b;
              {
                vector1<Vec> clash_olap;
                {
                  Vec chk( pose.xyz(AtomID(1,1)) );
                  clash_olap.push_back(    chk);
                  clash_olap.push_back(R3a*chk);
                  clash_olap.push_back(R3b*chk);
                }
                for(Size i3a = 0; i3a < 3; i3a++) {
                for(Size i2a = 0; i2a < 2; i2a++) {
                for(Size j3a = 0; j3a < 3; j3a++) {
                for(Size j2a = 0; j2a < 2; j2a++) {
                for(Size k3a = 0; k3a < 3; k3a++) {
                for(Size k2a = 0; k2a < 2; k2a++) {
                for(Size l3a = 0; l3a < 3; l3a++) {
                for(Size l2a = 0; l2a < 2; l2a++) {
                for(Size m3a = 0; m3a < 3; m3a++) {
                for(Size m2a = 0; m2a < 2; m2a++) {
                for(Size n3a = 0; n3a < 3; n3a++) {
                for(Size n2a = 0; n2a < 2; n2a++) {
									Vec chk( pose.xyz(AtomID(1,1)) ); Vec chk0 = chk;
									chk = R3[i3a] * chk; if(i2a) chk = R2*(chk-HG)+HG;
									chk = R3[j3a] * chk; if(j2a) chk = R2*(chk-HG)+HG;
									chk = R3[k3a] * chk; if(k2a) chk = R2*(chk-HG)+HG;
									chk = R3[l3a] * chk; if(l2a) chk = R2*(chk-HG)+HG;
									chk = R3[m3a] * chk; if(m2a) chk = R2*(chk-HG)+HG;
									chk = R3[n3a] * chk; if(n2a) chk = R2*(chk-HG)+HG;
									for(vector1<Vec>::const_iterator i = clash_olap.begin(); i != clash_olap.end(); ++i) {
										if( i->distance_squared(chk) < 1.0 ) goto cont6;
									} goto done6; cont6: continue; done6:
									clash_olap.push_back(chk);
									if( chk.distance_squared(chk0) > 30000 ) continue;
									for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
										for(Size ia = 1; ia <= 5; ++ia) {
											Vec chk( pose.xyz(AtomID(ia,ir)) );
											chk = R3[i3a] * chk; if(i2a) chk = R2*(chk-HG)+HG;
											chk = R3[j3a] * chk; if(j2a) chk = R2*(chk-HG)+HG;
											chk = R3[k3a] * chk; if(k2a) chk = R2*(chk-HG)+HG;
											chk = R3[l3a] * chk; if(l2a) chk = R2*(chk-HG)+HG;
											chk = R3[m3a] * chk; if(m2a) chk = R2*(chk-HG)+HG;
											chk = R3[n3a] * chk; if(n2a) chk = R2*(chk-HG)+HG;
											if( ! ifc3.clash_check(chk) ) goto cont7;
										}
									}
								}
								}
								}
								}
								}
								}
								}
								}
								}
								}
								}
                } goto done7; cont7: continue; done7:
                if(clash_olap.size() < 25) continue;
              }

              TR << "HIT " << irsd << " " << ch1 << " " << ch2 << " " << iaxs << " " << ANG << std::endl;
              seen_cen.push_back(R2*(com-HG)+HG);
              seen_axs.push_back(axs1);
              dumpsym2( pose, R2, R3a, R3b, HG, fname, irsd, ch1, ch2, ANG);

              Hit h;
              h.rsd = irsd;
              h.cbc = cbc;
              h.axs = axs1;
              h.cen = HG;
              h.chi1 = ch1;
              h.chi2 = ch2;
              h.sym = true;
              hits.push_back(h);
							foundhit = true;

              // std::cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
              // return std::pair<vector1<Hit>,vector1<Hit> > (hits,allhits);

              //pose.set_xyz( AtomID(pose.residue(irsd).atom_index("HG"),irsd), SG+axs1 );

              // for(Size i3a = 0; i3a < 3; i3a++) {
              //   for(Size i2a = 0; i2a < 2; i2a++) {
              //     Pose tmp(pose);
              //     rot_pose(tmp,R3[i3a]); if(i2a) rot_pose(tmp,R2,HG);
              //     string fn = utility::file_basename(fname) +"_"+  F(4,1,ANG)+"_"+ lzs(irsd,3)+"_"
              //       + lzs((Size)ch1,3)+"_"+ lzs((Size)ch2,3)+"_"+ lzs(i3a,1)+"_"+lzs(i2a,1)+".pdb";
              //     tmp.dump_pdb( fn );
              //   }
              // }
              //utility_exit_with_message("foo");

            }
          }
        }
      }
    }
    pose.replace_residue(irsd,*rprev,false);
  }

  return std::pair<vector1<Hit>,vector1<Hit> > (hits,allhits);
}

int main (int argc, char *argv[]) {

	try {

  devel::init(argc,argv);
  using namespace basic::options::OptionKeys;
  for(Size ifn = 1; ifn <= option[in::file::s]().size(); ++ifn) {
    string fn = option[in::file::s]()[ifn];
    Pose pnat;
    core::import_pose::pose_from_file(pnat,fn, core::import_pose::PDB_file);
    // if( pnat.n_residue() > 300 ) continue;
    for(Size ir = 2; ir <= pnat.n_residue()-1; ++ir) {
      if(!pnat.residue(ir).is_protein()) goto cont1;
      if(pnat.residue(ir).is_lower_terminus()) core::pose::remove_lower_terminus_type_from_pose_residue(pnat,ir);
      if(pnat.residue(ir).is_upper_terminus()) core::pose::remove_upper_terminus_type_from_pose_residue(pnat,ir);
      // if(pnat.residue(ir).name3()=="CYS") goto cont1;
    } goto done1; cont1: TR << "skipping " << fn << std::endl; continue; done1:
		//continue;
    Pose pala(pnat);
    std::pair<vector1<Hit>,vector1<Hit> > hts = dock(pala,fn);
    // design_hits(pnat,fn,hts.first,hts.second);
  }

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

