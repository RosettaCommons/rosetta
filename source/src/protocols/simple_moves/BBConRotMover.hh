// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @notes Monte Carlo backbone sampling for polypeptides with variable bond angles
/// and dihedral angles using concerted rotations and a Gaussian bias
/// Jakob P. Ulmschneider, JCP, VOl 118, 9

#ifndef INCLUDED_protocols_simple_moves_BBConRotMover_hh
#define INCLUDED_protocols_simple_moves_BBConRotMover_hh

#include <protocols/simple_moves/BBConRotMover.fwd.hh>
//core
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
// AUTO-REMOVED #include <core/chemical/ResidueType.fwd.hh>

//protocols
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/BBGaussianMover.hh>
//movemap
#include <core/kinematics/MoveMap.fwd.hh>
//
#include <utility/vector1.hh>

//std
#include <string>

namespace protocols {
namespace simple_moves {

class BBConRotMover : public BBGaussianMover
{
public:
    //typedef core::Vector Vector;
    //typedef core::Real Real;
    //typedef core::Size Size;
    //typedef core::pose::Pose Pose;
    //typedef core::pose::PoseOP PoseOP;
    //typedef core::pose::PoseCOP PoseCOP;
    //typedef core::conformation::Residue Residue;
    //typedef utility::vector1< Vector > VdRdPhi;
    //typedef utility::vector1< VdRdPhi > VMatrix;
    //typedef utility::vector1< utility::vector1<Real> > Matrix;
    //typedef numeric::xyzMatrix<Real> xyzMatrix;
    typedef  numeric::xyzMatrix< core::Real > xyzMatrix;

public:
    BBConRotMover();
    ~BBConRotMover();

    void apply(Pose &);
    virtual std::string get_name() const;
    void factorA( core::Real const fA );
    void factorB( core::Real const fB );
    void factorC( core::Real const fC );

protected:
		bool make_move(Pose &);
    void get_VdRdPhi(Pose const &);
    void get_G();
    void get_A();
    Real get_L_move(Pose &);
    Real get_L_prime();

private:
    /// @brief get xyz of d, given reference atoms a,b,c (a-b-c-d)
    /// and
    /// distance of c-d
    /// bond angle theta of d-c-b
    /// dihedral angle phi of d-c-b-a
    void get_xyz(Vector const &a, Vector const &b, Vector const &c, Vector &d,
                    Real distance, Real theta, Real phi);

    //key function, from Jakob
    Real calc_jacobian_cartesians(Vector const &, Vector const &,
        Vector const &, Vector const &, Vector const &);

    bool get_determinant(double a[5][5], int n, double &det);

    bool closure(
        //before closure
        Vector &r0,
        Vector &r1,
        Vector &r2,
        Vector &r3,
        Vector &r4,
        Vector &r5,
        Vector &r6,
        //after closure
        Vector &p1,
        Vector &p2,
        Vector &p3,
        //old angle/dih
        Real const theta1_old,
        Real const theta2_old,
        Real const theta3_old,
        Real const theta4_old,
        Real const alpha1_old,
        Real const alpha2_old,
        Real const alpha3_old,
        //new angle/dih
        Real &theta1,
        Real &theta2,
        Real &theta3,
        Real &theta4,
        Real &alpha1,
        Real &alpha2,
        Real &alpha3
    );

private:
    utility::vector1< Real > dphi;
    utility::vector1< Real > oldphi;
    Real factorA_;
    Real factorB_;
    Real factorC_;

    //no sidechain and branch atom
    //core::chemical::ResidueType const &cen_gly_restyp;
};

}//moves
}//protocols

#endif




