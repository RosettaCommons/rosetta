// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
/*
 * DisulfideBondEnergy.cc
 *
 *  Created on: Mar 9, 2009
 *      Author: dgront
 */

#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <protocols/Protocol.hh>
#include <protocols/abinitio/AbrelaxApplication.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option_macros.hh>

#include <core/id/NamedStubID.hh>
#include <core/kinematics/Stub.hh>
#include <core/id/types.hh>
#include <numeric/xyz.functions.hh>
#include <core/kinematics/types.hh>

#include <utility/excn/Exceptions.hh>

// Numeric headers
#include <numeric/xyzMatrix.hh>

using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;

OPT_1GRP_KEY( Integer, in, cysid )

void register_options() {

  OPT( in::file::native );
  NEW_OPT( in::cysid, "id of a CYS residue in a disulfide bond to be tested" , 1 );
  option.add_relevant(in::cysid);
}

class DisulfideBondEnergy {

  public:
    static double energyBondCutoff;
    static double pseudocountsFraction;

    DisulfideBondEnergy(double, double, double);
    double probability(double ca_ca_sg, double cb_sg_sg_cb, double sg_cb_ca, double d);
    double probability(core::PointPosition cb1, core::PointPosition sg1, core::PointPosition sg2, core::PointPosition cb2);
    bool detectDisulfideBond(pose::Pose &pose, int cysId1, int cysId2);

    double evaluate(pose::Pose &pose, int cysId1, int cysId2);
    bool rebuildAndDetectDisulfideBond(pose::Pose &pose, int cysId1, double chi1, int cysId2, double chi2);
    bool rebuildAndDetectDisulfideBond(pose::Pose &pose, int cysId1, int cysId2);
    double tabulate(pose::Pose &pose, int cysId1, int cysId2);

  private:
    static double theta;
    static double dist;
    static double meanPos[];
    double* stdevPos;
    static double meanNeg[];
    double* stdevNeg;
    static double maxCBDistance;
    static double sqrtTwoPi;
    double* twoVarPos;
    double* twoVarNeg;
};

double DisulfideBondEnergy::energyBondCutoff = -0.3;
double DisulfideBondEnergy::pseudocountsFraction = 0.000000001;

double DisulfideBondEnergy::meanPos[] = { 104.4, 93.6, 104.4, 2.055 };
double DisulfideBondEnergy::meanNeg[] = { 103.8, -85.89, 103.9, 2.055 };
double DisulfideBondEnergy::theta = 1.19426;
double DisulfideBondEnergy::dist = 1.808;
double DisulfideBondEnergy::maxCBDistance = 4.5;
double DisulfideBondEnergy::sqrtTwoPi = sqrt(M_PI * 2.0);

DisulfideBondEnergy::DisulfideBondEnergy(double planarWidthScaling, double torsionWidthScaling, double distanceWidthScaling) {

  stdevPos = new double[4];
  stdevNeg = new double[4];
  stdevPos[0] = 2.3 * planarWidthScaling;
  stdevPos[1] = 11.16 * torsionWidthScaling;
  stdevPos[2] = 2.3 * planarWidthScaling;
  stdevPos[3] = 0.048 * distanceWidthScaling;
  stdevNeg[0] = 2.14 * planarWidthScaling;
  stdevNeg[1] = 8.43 * torsionWidthScaling;
  stdevNeg[2] = 2.0 * planarWidthScaling;
  stdevNeg[3] = 0.048 * distanceWidthScaling;

  twoVarPos = new double[4];
  twoVarNeg = new double[4];
  twoVarPos[0] = 2 * stdevPos[0] * stdevPos[0];
  twoVarPos[1] = 2 * stdevPos[1] * stdevPos[1];
  twoVarPos[2] = 2 * stdevPos[2] * stdevPos[2];
  twoVarPos[3] = 2 * stdevPos[3] * stdevPos[3];
  twoVarNeg[0] = 2 * stdevNeg[0] * stdevNeg[0];
  twoVarNeg[1] = 2 * stdevNeg[1] * stdevNeg[1];
  twoVarNeg[2] = 2 * stdevNeg[2] * stdevNeg[2];
  twoVarNeg[3] = 2 * stdevNeg[3] * stdevNeg[3];
}
double DisulfideBondEnergy::probability(double cb_sg_sg, double cb_sg_sg_cb, double sg_sg_cb, double d) {

  double x[] = { cb_sg_sg, cb_sg_sg_cb, sg_sg_cb, d };
  double val = 1.0;
  if (cb_sg_sg_cb > 0) {
    for (int i = 0; i < 4; i++) {
      double t = x[i] - meanPos[i];
      val *= exp(-t * t / (twoVarPos[i])) / (stdevPos[i] * sqrtTwoPi);
      //      std::cerr << x[i] << " " << meanPos[i] << " " << stdevPos[i] << " " << exp(-t * t / (twoVarPos[i])) / (stdevPos[i]
      //          * sqrtTwoPi) << "\n";
    }
  } else {
    for (int i = 0; i < 4; i++) {
      double t = x[i] - meanNeg[i];
      val *= exp(-t * t / (twoVarNeg[i])) / (stdevNeg[i] * sqrtTwoPi);
      //      std::cerr << x[i] << " " << meanNeg[i] << " " << stdevNeg[i] << " " << exp(-t * t / (twoVarNeg[i])) / (stdevNeg[i]
      //          * sqrtTwoPi) << "\n";
    }
  }

  return val;
}

double DisulfideBondEnergy::probability(core::PointPosition cb1, core::PointPosition sg1, core::PointPosition sg2,
    core::PointPosition cb2) {

  double cb_sg_sg_cb = dihedral_degrees(cb1, sg1, sg2, cb2);
  double cb_sg_sg = angle_radians(cb1, sg1, sg2) * 180.0 / M_PI;
  double sg_sg_cb = angle_radians(sg1, sg2, cb2) * 180.0 / M_PI;

  //  std::cerr << cb_sg_sg << " " << cb_sg_sg_cb << " " << sg_sg_cb << " " << distance(sg1, sg2) << "\n";

  return probability(cb_sg_sg, cb_sg_sg_cb, sg_sg_cb, distance(sg1, sg2));
}

double DisulfideBondEnergy::evaluate(pose::Pose &pose, int cysId1, int cysId2) {

  Size const seqpos1(cysId1); //the first CYS
  id::NamedStubID CB_stub_id1("CB", "CA", "N", seqpos1);
  kinematics::Stub stub1(pose.stub_from_id(CB_stub_id1));
  core::PointPosition n1 = pose.xyz(id::NamedAtomID("N", seqpos1));
  core::PointPosition ca1 = pose.xyz(id::NamedAtomID("CA", seqpos1));
  core::PointPosition cb1 = pose.xyz(id::NamedAtomID("CB", seqpos1));
  core::PointPosition sg1 = pose.xyz(id::NamedAtomID("SG", seqpos1));

  Size const seqpos2(cysId2); //the second CYS
  id::NamedStubID CB_stub_id2("CB", "CA", "N", seqpos2);
  kinematics::Stub stub2(pose.stub_from_id(CB_stub_id2));
  core::PointPosition n2 = pose.xyz(id::NamedAtomID("N", seqpos2));
  core::PointPosition ca2 = pose.xyz(id::NamedAtomID("CA", seqpos2));
  core::PointPosition cb2 = pose.xyz(id::NamedAtomID("CB", seqpos2));
  core::PointPosition sg2 = pose.xyz(id::NamedAtomID("SG", seqpos2));

  id::TorsionID chi1_id(seqpos1, id::CHI, 1);
  std::cout << pose.torsion(chi1_id) << std::endl;
  id::TorsionID chi2_id(seqpos2, id::CHI, 1);
  std::cout << pose.torsion(chi2_id) << std::endl;

  double p = probability(cb1, sg1, sg2, cb2);

  std::cout << sg1[0] << " " << sg1[1] << " " << sg1[2] << "\n";
  std::cout << sg2[0] << " " << sg2[1] << " " << sg2[2] << "\n";

  return (-log(pseudocountsFraction + p) + log(pseudocountsFraction));
}

bool DisulfideBondEnergy::rebuildAndDetectDisulfideBond(pose::Pose &pose, int cysId1, int cysId2) {

  //  double chiCrude[] = { -63, -176, 63 };
  double chiCrude[] = { -72.029, -69.2231, 0 };

  int nChi = 3;
  Size const seqpos1(cysId1); //the first CYS
  id::NamedStubID CB_stub_id1("CB", "CA", "N", seqpos1);
  kinematics::Stub stub1(pose.stub_from_id(CB_stub_id1));
  core::PointPosition cb1 = pose.xyz(id::NamedAtomID("CB", seqpos1));

  Size const seqpos2(cysId2); //the second CYS
  id::NamedStubID CB_stub_id2("CB", "CA", "N", seqpos2);
  kinematics::Stub stub2(pose.stub_from_id(CB_stub_id2));
  core::PointPosition cb2 = pose.xyz(id::NamedAtomID("CB", seqpos2));

  if (distance(cb1, cb2) > maxCBDistance)
    return false;

  for (int i1 = 0; i1 < nChi; i1++) {
    double chi1 = chiCrude[i1];
    for (int i2 = 0; i2 < nChi; i2++) {
      double chi2 = chiCrude[i2];
      kinematics::Stub new_stub1(stub1.M * numeric::x_rotation_matrix_degrees(chi1) * numeric::z_rotation_matrix_radians(theta),
          stub1.v);
      new_stub1.v += dist * new_stub1.M.col_x();
      kinematics::Stub new_stub2(stub2.M * numeric::x_rotation_matrix_degrees(chi2) * numeric::z_rotation_matrix_radians(theta),
          stub2.v);
      new_stub2.v += dist * new_stub2.M.col_x();
      double d = distance(new_stub1.v, new_stub2.v);
      if (d > 3.0)
        continue;
      for (double eps1 = -10; eps1 <= 10.01; eps1 += 1.5) {
        for (double eps2 = -10; eps2 <= 10.01; eps2 += 1.5) {
          kinematics::Stub new_stub1E(stub1.M * numeric::x_rotation_matrix_degrees(chi1 + eps1)
              * numeric::z_rotation_matrix_radians(theta), stub1.v);
          new_stub1E.v += dist * new_stub1E.M.col_x();
          kinematics::Stub new_stub2E(stub2.M * numeric::x_rotation_matrix_degrees(chi2 + eps2)
              * numeric::z_rotation_matrix_radians(theta), stub2.v);
          new_stub2E.v += dist * new_stub2E.M.col_x();
          if (distance(new_stub1E.v, new_stub2E.v) > 2.5)
            continue;
          double p = probability(cb1, new_stub1E.v, new_stub2E.v, cb2);
//          std::cout << (chi1 + eps1) << " " << (chi2 + eps2) << " " << p << " " << distance(new_stub1E.v, new_stub2E.v) << " "
//              << (-log(pseudocountsFraction + p) + log(pseudocountsFraction)) << "\n";
          if (-log(pseudocountsFraction + p) + log(pseudocountsFraction) < energyBondCutoff)
            return true;
        }
      }
    }
  }

  return false;
}

bool DisulfideBondEnergy::rebuildAndDetectDisulfideBond(pose::Pose &pose, int cysId1, double chi1, int cysId2, double chi2) {

  Size const seqpos1(cysId1); //the first CYS
  id::NamedStubID CB_stub_id1("CB", "CA", "N", seqpos1);
  kinematics::Stub stub1(pose.stub_from_id(CB_stub_id1));
  core::PointPosition cb1 = pose.xyz(id::NamedAtomID("CB", seqpos1));

  Size const seqpos2(cysId2); //the second CYS
  id::NamedStubID CB_stub_id2("CB", "CA", "N", seqpos2);
  kinematics::Stub stub2(pose.stub_from_id(CB_stub_id2));
  core::PointPosition cb2 = pose.xyz(id::NamedAtomID("CB", seqpos2));

  kinematics::Stub new_stub1(stub1.M * numeric::x_rotation_matrix_degrees(chi1) * numeric::z_rotation_matrix_radians(theta),
      stub1.v);
  new_stub1.v += dist * new_stub1.M.col_x();
  kinematics::Stub new_stub2(stub2.M * numeric::x_rotation_matrix_degrees(chi2) * numeric::z_rotation_matrix_radians(theta),
      stub2.v);
  new_stub2.v += dist * new_stub2.M.col_x();

  double p = probability(cb1, new_stub1.v, new_stub2.v, cb2);

  std::cout << new_stub1.v[0] << " " << new_stub1.v[1] << " " << new_stub1.v[2] << "\n";
  std::cout << new_stub2.v[0] << " " << new_stub2.v[1] << " " << new_stub2.v[2] << "\n";
  std::cout << (-log(pseudocountsFraction + p) + log(pseudocountsFraction)) << "\n";

  if (-log(pseudocountsFraction + p) + log(pseudocountsFraction) < energyBondCutoff)
    return true;
  return false;
}

double DisulfideBondEnergy::tabulate(pose::Pose &pose, int cysId1, int cysId2) {

  Size const seqpos1(cysId1); //the first CYS
  id::NamedStubID CB_stub_id1("CB", "CA", "N", seqpos1);
  kinematics::Stub stub1(pose.stub_from_id(CB_stub_id1));
  core::PointPosition cb1 = pose.xyz(id::NamedAtomID("CB", seqpos1));

  Size const seqpos2(cysId2); //the second CYS
  id::NamedStubID CB_stub_id2("CB", "CA", "N", seqpos2);
  kinematics::Stub stub2(pose.stub_from_id(CB_stub_id2));
  core::PointPosition cb2 = pose.xyz(id::NamedAtomID("CB", seqpos2));

  for (double chi1 = -180.0; chi1 < 180.0; chi1 += 2.5) {
    for (double chi2 = -180.0; chi2 < 180.0; chi2 += 2.5) {
      kinematics::Stub new_stub1(stub1.M * numeric::x_rotation_matrix_degrees(chi1) * numeric::z_rotation_matrix_radians(theta),
          stub1.v);
      new_stub1.v += dist * new_stub1.M.col_x();
      kinematics::Stub new_stub2(stub2.M * numeric::x_rotation_matrix_degrees(chi2) * numeric::z_rotation_matrix_radians(theta),
          stub2.v);
      new_stub2.v += dist * new_stub2.M.col_x();
      double p = probability(cb1, new_stub1.v, new_stub2.v, cb2);
      std::cout << chi1 << " " << chi2 << " " << p << " " << -log(pseudocountsFraction + p) + log(pseudocountsFraction)
          << std::endl;
    }
  }
  return 0.0;
}

int main(int argc, char * argv[]) {

  try {
  devel::init(argc, argv);
  register_options();

  pose::Pose init_pose;

  int id1 = 131;
  int id2 = 138;
  //  if (option[in::cysid].user()) {
  //    id = option[in::cysid]();
  //  }

  DisulfideBondEnergy* ssEn = new DisulfideBondEnergy(1.5, 1.5, 1.5);

  if (option[in::file::native].user()) {
    core::import_pose::pose_from_file(init_pose, option[in::file::native](), core::import_pose::PDB_file);
    std::cout << ssEn->evaluate(init_pose, id1, id2) << std::endl;
    std::cout << ssEn->rebuildAndDetectDisulfideBond(init_pose, id1, id2) << std::endl;
    std::cout << ssEn->rebuildAndDetectDisulfideBond(init_pose, id1, -72.029, id2, -69.2231) << std::endl;
    //    tabulate(init_pose, id1, id2);
  }
  } catch (utility::excn::Exception const & e ) {
                            std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                }
      return 0;

}
