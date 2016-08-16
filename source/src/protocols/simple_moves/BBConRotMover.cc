// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


#include <protocols/simple_moves/BBConRotMover.hh>
//core
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueType.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>

#include <core/pose/Pose.hh>
#include <basic/basic.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/bbg.OptionKeys.gen.hh>

#include <core/kinematics/Stub.hh>

//util
#include <utility/exit.hh>
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/constants.hh>
#include <numeric/internal/RowVectors.hh>

#include <iostream>
#include <sstream>

#include <utility/vector1.hh>


using namespace std;
using namespace core;
using namespace core::pose;
using namespace utility;
using namespace numeric;
using namespace basic::options;
using namespace basic::options::OptionKeys;

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.BBConRotMover" );

namespace protocols {
namespace simple_moves {

void BBConRotMover::factorA( core::Real const fA )
{
	factorA_ = fA;
}

void BBConRotMover::factorB( core::Real const fB )
{
	factorB_ = fB;
}

void BBConRotMover::factorC( core::Real const fC )
{
	factorC_ = fC;
}

BBConRotMover::BBConRotMover()
:BBGaussianMover(1,15,5),
	dphi(Vector(n_dof_angle_)),
	oldphi(Vector(n_dof_angle_))
{
	//using numeric::constants::d::pi;

	protocols::moves::Mover::type("BBConRotMover");

	//in bbg, set_phi(degree), 1 (A/2) ~ 57.3/2
	//here, set_dof(rad)
	//paper param: 50 (C1) ~ 0.87 -> A~1.74

	//init the ABC factor
	//factorA_ = option[ bbg::factorA ]*180.0/pi;
	factorA_ = option[ bbg::factorA ]*2.0; //the g() is (0,1), not (0,1/sqrt(2))
	factorB_ = option[ bbg::factorB ];
	factorC_ = 20.0;
}

BBConRotMover::~BBConRotMover(){}

std::string BBConRotMover::get_name() const
{
	return "BBConRotMover";
}

void BBConRotMover::apply(Pose &pose)
{
	Size iter=0;
	while ( make_move(pose) )
			{
		if ( iter++ > 100 ) break;
	}
	TR.Debug << "apply: iter=" << iter << std::endl;
}

bool BBConRotMover::make_move(Pose &pose)
{
	using basic::periodic_range;
	using basic::unsigned_periodic_range;
	using numeric::constants::d::pi;
	//using numeric::constants::d::pi_2;

	setup_list(pose);
	int ndx=static_cast< int >( numeric::random::rg().uniform()*available_seg_list_.size()+1 );

	Size left = available_seg_list_[ ndx ].first;
	resnum_ = available_seg_list_[ ndx ].second;
	TR.Debug << "Pick: " << left << " <--> " << resnum_ << std::endl;
	if ( resnum_-left+1 < n_pert_res_ ) {
		//do random
		TR.Debug << "Do random rot ... " << std::endl;
		pivot_range_randomly(pose, left, resnum_);
		return false;
	}
	assert(resnum_-left == n_pert_res_-1);

	Size nres(pose.n_residue());
	xyzVector oldv(pose.residue(nres).atom("CA").xyz());

	//using whole pose
	xyzVector r0(pose.residue(resnum_).atom("C").xyz());
	xyzVector r1(pose.residue(resnum_).atom("CA").xyz());
	xyzVector r2(pose.residue(resnum_).atom("N").xyz());
	xyzVector r3(pose.residue(resnum_-1).atom("C").xyz());
	xyzVector r4(pose.residue(resnum_-1).atom("CA").xyz());
	xyzVector r5(pose.residue(resnum_-1).atom("N").xyz());
	xyzVector r6(pose.residue(resnum_-2).atom("C").xyz());
	xyzVector p1((r2-r1).length(),0.0,0.0);
	xyzVector p2((r3-r2).length(),0.0,0.0);
	xyzVector p3((r4-r3).length(),0.0,0.0);

	Real jac_new0 = calc_jacobian_cartesians(r5, r4, r3, r2, r1);
	TR.Debug << "Jac_new0=" << jac_new0 << std::endl;

	bool failed = true;
	//check jac>0
	if ( jac_new0>0 ) failed = false;

	if ( !failed ) {

		id::DOF_ID dih12( id::AtomID(pose.residue(resnum_  ).atom_index("C" ), resnum_  ), id::PHI );
		id::DOF_ID dih11( id::AtomID(pose.residue(resnum_  ).atom_index("CA"), resnum_  ), id::PHI );
		id::DOF_ID dih10( id::AtomID(pose.residue(resnum_  ).atom_index("N" ), resnum_  ), id::PHI );
		id::DOF_ID dih9 ( id::AtomID(pose.residue(resnum_-1).atom_index("C" ), resnum_-1), id::PHI );
		id::DOF_ID ang11( id::AtomID(pose.residue(resnum_  ).atom_index("CA"), resnum_  ), id::THETA );
		id::DOF_ID ang10( id::AtomID(pose.residue(resnum_  ).atom_index("N" ), resnum_  ), id::THETA );
		id::DOF_ID ang9 ( id::AtomID(pose.residue(resnum_-1).atom_index("C" ), resnum_-1), id::THETA );
		Real theta1_old = pose.dof(dih12); //dih_12
		Real theta2_old = pose.dof(dih11);; //dih_11
		Real theta3_old = pose.dof(dih10); //dih_10
		Real theta4_old = pose.dof(dih9); //dih_9
		Real alpha1_old = pi - pose.dof(ang11); //ang_11
		Real alpha2_old = pi - pose.dof(ang10); //ang_10
		Real alpha3_old = pi - pose.dof(ang9); //ang_9

		//perturb forward
		get_VdRdPhi(pose);
		get_G();
		get_A();
		//TR.Debug << "perturb..." << std::endl;
		Real W_old = get_L_move(pose);
		//backward
		get_VdRdPhi(pose);
		get_G();
		get_A();
		Real W_new = get_L_prime();

		//proposal density
		last_proposal_density_ratio_ = W_new / W_old;
		TR.Debug << "W_old=" << W_old << " W_new=" << W_new << std::endl;
		TR.Debug << "ratio=" << last_proposal_density_ratio_ << std::endl;

		//new r6, r5, r4
		r4 = pose.residue(resnum_-1).atom("CA").xyz();
		r5 = pose.residue(resnum_-1).atom("N").xyz();
		r6 = pose.residue(resnum_-2).atom("C").xyz();

		//closure
		Real theta1,theta2,theta3,theta4;
		Real alpha1,alpha2,alpha3;


		//TR<<"Closure..." << std::endl;
		failed=true;

		if ( !closure(
				//Vector
				r0,r1,r2,r3,r4,r5,r6,p1,p2,p3,
				//old angle/dih
				theta1_old,theta2_old,theta3_old,theta4_old,
				alpha1_old,alpha2_old,alpha3_old,
				//new angle/dih
				theta1,theta2,theta3,theta4,
				alpha1,alpha2,alpha3) ) {
			//successfully close
			Real jac_new1 = calc_jacobian_cartesians(r5, r4, r3, r2, r1);
			TR.Debug << "Jac_new1=" << jac_new1 << std::endl;

			if ( jac_new1>0 ) {
				last_proposal_density_ratio_ *= jac_new1 / jac_new0;
				TR.Debug << "last_proposal=" << last_proposal_density_ratio_ << std::endl;

				//using whole pose
				pose.set_dof(dih12, theta1); //dih_12
				pose.set_dof(dih11, theta2); //dih_11
				pose.set_dof(dih10, theta3); //dih_10
				pose.set_dof(dih9, theta4); //dih_9
				pose.set_dof(ang11, pi-alpha1); //ang_11
				pose.set_dof(ang10, pi-alpha2); //ang_10
				pose.set_dof(ang9, pi-alpha3); //ang_9

				//make sure the downstream didn't flip
				xyzVector dd(pose.residue(nres).atom("CA").xyz()-oldv);
				if ( dd.length_squared()<1.0e-6 ) {
					failed=false;
					TR.Debug << "Closure Success!" << std::endl;
				}
			}
		}

	}//if (!failed)

	if ( failed ) {
		//false
		//back
		id::DOF_ID dih0(   id::AtomID(pose.residue(resnum_-4).atom_index("C" ), resnum_-4), id::PHI );
		id::DOF_ID dih1(   id::AtomID(pose.residue(resnum_-3).atom_index("N" ), resnum_-3), id::PHI );
		id::DOF_ID dih3(   id::AtomID(pose.residue(resnum_-3).atom_index("C" ), resnum_-3), id::PHI );
		id::DOF_ID dih4(   id::AtomID(pose.residue(resnum_-2).atom_index("N" ), resnum_-2), id::PHI );
		id::DOF_ID dih6(   id::AtomID(pose.residue(resnum_-2).atom_index("C" ), resnum_-2), id::PHI );
		id::DOF_ID dih7(   id::AtomID(pose.residue(resnum_-1).atom_index("N" ), resnum_-1), id::PHI );
		id::DOF_ID angle1( id::AtomID(pose.residue(resnum_-4).atom_index("C" ), resnum_-4), id::THETA );
		id::DOF_ID angle2( id::AtomID(pose.residue(resnum_-3).atom_index("N" ), resnum_-3), id::THETA );
		id::DOF_ID angle3( id::AtomID(pose.residue(resnum_-3).atom_index("CA"), resnum_-3), id::THETA );
		id::DOF_ID angle4( id::AtomID(pose.residue(resnum_-3).atom_index("C" ), resnum_-3), id::THETA );
		id::DOF_ID angle5( id::AtomID(pose.residue(resnum_-2).atom_index("N" ), resnum_-2), id::THETA );
		id::DOF_ID angle6( id::AtomID(pose.residue(resnum_-2).atom_index("CA"), resnum_-2), id::THETA );
		id::DOF_ID angle7( id::AtomID(pose.residue(resnum_-2).atom_index("C" ), resnum_-2), id::THETA );
		id::DOF_ID angle8( id::AtomID(pose.residue(resnum_-1).atom_index("N" ), resnum_-1), id::THETA );
		id::DOF_ID angle9( id::AtomID(pose.residue(resnum_-1).atom_index("CA"), resnum_-1), id::THETA );

		//TR<< "Rebuild ..." << std::endl;
		pose.set_dof(dih0, oldphi[1]);
		pose.set_dof(dih1, oldphi[2]);
		pose.set_dof(dih3, oldphi[3]);
		pose.set_dof(dih4, oldphi[4]);
		pose.set_dof(dih6, oldphi[5]);
		pose.set_dof(dih7, oldphi[6]);
		pose.set_dof(angle1, oldphi[7]);
		pose.set_dof(angle2, oldphi[8]);
		pose.set_dof(angle3, oldphi[9]);
		pose.set_dof(angle4, oldphi[10]);
		pose.set_dof(angle5, oldphi[11]);
		pose.set_dof(angle6, oldphi[12]);
		pose.set_dof(angle7, oldphi[13]);
		pose.set_dof(angle8, oldphi[14]);
		pose.set_dof(angle9, oldphi[15]);

		last_proposal_density_ratio_ = 1.0;
		TR.Debug << "Closure Failed!" << std::endl;
	}

	return failed;
}

void BBConRotMover::get_VdRdPhi(Pose const &segment)
{
	Size nres=resnum_; //using the whole pose
	conformation::Residue const & rsd0(segment.residue(nres-4));
	conformation::Residue const & rsd1(segment.residue(nres-3));
	conformation::Residue const & rsd2(segment.residue(nres-2));
	conformation::Residue const & rsd3(segment.residue(nres-1));

	//the only end
	xyzVector end_xyz = rsd3.atom("CA").xyz();

	///////////////////
	//dihedral x 6
	///////////////////
	matrix_dRdPhi[1][1] = get_dRdPhi(
		rsd0.atom("N").xyz(),
		rsd0.atom("CA").xyz(),
		end_xyz);

	matrix_dRdPhi[1][2] = get_dRdPhi(
		rsd0.atom("CA").xyz(),
		rsd0.atom("C").xyz(),
		end_xyz);

	matrix_dRdPhi[1][3] = get_dRdPhi(
		rsd1.atom("N").xyz(),
		rsd1.atom("CA").xyz(),
		end_xyz);

	matrix_dRdPhi[1][5] = get_dRdPhi(
		rsd1.atom("CA").xyz(),
		rsd1.atom("C").xyz(),
		end_xyz);

	matrix_dRdPhi[1][5] = get_dRdPhi(
		rsd2.atom("N").xyz(),
		rsd2.atom("CA").xyz(),
		end_xyz);

	matrix_dRdPhi[1][6] = get_dRdPhi(
		rsd2.atom("CA").xyz(),
		rsd2.atom("C").xyz(),
		end_xyz);

	////////////////////////
	//bond angle x 9
	////////////////////////
	matrix_dRdPhi[1][7] = get_dRdTheta(
		rsd0.atom("N").xyz(),
		rsd0.atom("CA").xyz(),
		rsd0.atom("C").xyz(),
		end_xyz);

	matrix_dRdPhi[1][8] = get_dRdTheta(
		rsd0.atom("CA").xyz(),
		rsd0.atom("C").xyz(),
		rsd1.atom("N").xyz(),
		end_xyz);

	matrix_dRdPhi[1][9] = get_dRdTheta(
		rsd0.atom("C").xyz(),
		rsd1.atom("N").xyz(),
		rsd1.atom("CA").xyz(),
		end_xyz);

	matrix_dRdPhi[1][10] = get_dRdTheta(
		rsd1.atom("N").xyz(),
		rsd1.atom("CA").xyz(),
		rsd1.atom("C").xyz(),
		end_xyz);

	matrix_dRdPhi[1][11] = get_dRdTheta(
		rsd1.atom("CA").xyz(),
		rsd1.atom("C").xyz(),
		rsd2.atom("N").xyz(),
		end_xyz);

	matrix_dRdPhi[1][12] = get_dRdTheta(
		rsd1.atom("C").xyz(),
		rsd2.atom("N").xyz(),
		rsd2.atom("CA").xyz(),
		end_xyz);

	matrix_dRdPhi[1][13] = get_dRdTheta(
		rsd2.atom("N").xyz(),
		rsd2.atom("CA").xyz(),
		rsd2.atom("C").xyz(),
		end_xyz);

	matrix_dRdPhi[1][14] = get_dRdTheta(
		rsd2.atom("CA").xyz(),
		rsd2.atom("C").xyz(),
		rsd3.atom("N").xyz(),
		end_xyz);

	matrix_dRdPhi[1][15] = get_dRdTheta(
		rsd2.atom("C").xyz(),
		rsd3.atom("N").xyz(),
		rsd3.atom("CA").xyz(),
		end_xyz);
}

void BBConRotMover::get_G()
{
	for ( Size i=1; i<=n_dof_angle_; i++ ) {
		for ( Size j=i; j<=n_dof_angle_; j++ ) {
			matrix_G[i][j] = matrix_dRdPhi[1][i].dot(matrix_dRdPhi[1][j]);
			if ( i<j ) matrix_G[j][i]=matrix_G[i][j];
		}
	}
}

void BBConRotMover::get_A()
{
	//A = a(1+bG)
	//A = L^-1 * L
	//L(angle) -> c * L(angle)
	for ( Size i=1; i<=n_dof_angle_; i++ ) {
		for ( Size j=i; j<=n_dof_angle_; j++ ) {
			matrix_A[i][j] = factorB_ * matrix_G[i][j];
			if ( i==j ) matrix_A[i][j] += 1.0;
			matrix_A[i][j] *= factorA_;

			if ( i<j ) matrix_A[j][i] = matrix_A[i][j];
		}
	}
}

core::Real BBConRotMover::get_L_move(Pose &segment)
{
	using basic::periodic_range;
	using basic::unsigned_periodic_range;
	using numeric::constants::d::pi;
	using numeric::constants::d::pi_2;

	Size nres=resnum_; //using the whole pose
	//Size nres=6; //using copy segment

	//gerate a Gaussian dx vector
	Vector delta(n_dof_angle_);
	for ( Size i=1; i<=n_dof_angle_; i++ ) delta[i]=numeric::random::rg().gaussian();
	//Debug: no angle changes
	//for (Size i=7; i<=n_dof_angle_; i++) delta[i]=0;

	Real d2=0.0;
	for ( Size i=1; i<=n_dof_angle_; i++ ) d2+=delta[i]*delta[i];

	//cholesky, get L^t, L^-1
	Real detL = cholesky_fw(matrix_A, n_dof_angle_, delta, dphi, 7, 15, factorC_);

	//should i use factorC here or in get_A?
	Real W_old = detL*exp(-d2/2.0);

	//set the new phi, psi, theta
	//res total: n_pert_res_+2, from 2 to n_pert_res_+1
	id::DOF_ID dih0(   id::AtomID(segment.residue(nres-4).atom_index("C" ), nres-4), id::PHI );
	id::DOF_ID dih1(   id::AtomID(segment.residue(nres-3).atom_index("N" ), nres-3), id::PHI );
	id::DOF_ID dih3(   id::AtomID(segment.residue(nres-3).atom_index("C" ), nres-3), id::PHI );
	id::DOF_ID dih4(   id::AtomID(segment.residue(nres-2).atom_index("N" ), nres-2), id::PHI );
	id::DOF_ID dih6(   id::AtomID(segment.residue(nres-2).atom_index("C" ), nres-2), id::PHI );
	id::DOF_ID dih7(   id::AtomID(segment.residue(nres-1).atom_index("N" ), nres-1), id::PHI );
	id::DOF_ID angle1( id::AtomID(segment.residue(nres-4).atom_index("C" ), nres-4), id::THETA );
	id::DOF_ID angle2( id::AtomID(segment.residue(nres-3).atom_index("N" ), nres-3), id::THETA );
	id::DOF_ID angle3( id::AtomID(segment.residue(nres-3).atom_index("CA"), nres-3), id::THETA );
	id::DOF_ID angle4( id::AtomID(segment.residue(nres-3).atom_index("C" ), nres-3), id::THETA );
	id::DOF_ID angle5( id::AtomID(segment.residue(nres-2).atom_index("N" ), nres-2), id::THETA );
	id::DOF_ID angle6( id::AtomID(segment.residue(nres-2).atom_index("CA"), nres-2), id::THETA );
	id::DOF_ID angle7( id::AtomID(segment.residue(nres-2).atom_index("C" ), nres-2), id::THETA );
	id::DOF_ID angle8( id::AtomID(segment.residue(nres-1).atom_index("N" ), nres-1), id::THETA );
	id::DOF_ID angle9( id::AtomID(segment.residue(nres-1).atom_index("CA"), nres-1), id::THETA );

	oldphi[1]=segment.dof(dih0);
	oldphi[2]=segment.dof(dih1);
	oldphi[3]=segment.dof(dih3);
	oldphi[4]=segment.dof(dih4);
	oldphi[5]=segment.dof(dih6);
	oldphi[6]=segment.dof(dih7);
	oldphi[7]=segment.dof(angle1);
	oldphi[8]=segment.dof(angle2);
	oldphi[9]=segment.dof(angle3);
	oldphi[10]=segment.dof(angle4);
	oldphi[11]=segment.dof(angle5);
	oldphi[12]=segment.dof(angle6);
	oldphi[13]=segment.dof(angle7);
	oldphi[14]=segment.dof(angle8);
	oldphi[15]=segment.dof(angle9);

	segment.set_dof(dih0, periodic_range(segment.dof(dih0)+dphi[1],pi_2));
	segment.set_dof(dih1, periodic_range(segment.dof(dih1)+dphi[2],pi_2));
	segment.set_dof(dih3, periodic_range(segment.dof(dih3)+dphi[3],pi_2));
	segment.set_dof(dih4, periodic_range(segment.dof(dih4)+dphi[4],pi_2));
	segment.set_dof(dih6, periodic_range(segment.dof(dih6)+dphi[5],pi_2));
	segment.set_dof(dih7, periodic_range(segment.dof(dih7)+dphi[6],pi_2));
	segment.set_dof(angle1, unsigned_periodic_range(segment.dof(angle1)+dphi[7],pi));
	segment.set_dof(angle2, unsigned_periodic_range(segment.dof(angle2)+dphi[8],pi));
	segment.set_dof(angle3, unsigned_periodic_range(segment.dof(angle3)+dphi[9],pi));
	segment.set_dof(angle4, unsigned_periodic_range(segment.dof(angle4)+dphi[10],pi));
	segment.set_dof(angle5, unsigned_periodic_range(segment.dof(angle5)+dphi[11],pi));
	segment.set_dof(angle6, unsigned_periodic_range(segment.dof(angle6)+dphi[12],pi));
	segment.set_dof(angle7, unsigned_periodic_range(segment.dof(angle7)+dphi[13],pi));
	segment.set_dof(angle8, unsigned_periodic_range(segment.dof(angle8)+dphi[14],pi));
	segment.set_dof(angle9, unsigned_periodic_range(segment.dof(angle9)+dphi[15],pi));

	return W_old;
}

core::Real BBConRotMover::get_L_prime()
{
	Vector delta(n_dof_angle_);
	Real detL = cholesky_bw(matrix_A, n_dof_angle_, dphi, delta, 7, 15, factorC_);
	Real d2=0.0;
	for ( Size i=1; i<=n_dof_angle_; i++ ) d2+=delta[i]*delta[i];
	return detL*exp(-d2/2.0);
}

void BBConRotMover::get_xyz(
	xyzVector const &a,
	xyzVector const &b,
	xyzVector const &c,
	xyzVector &d,
	Real distance,
	Real theta,
	Real phi
)
{
	using numeric::x_rotation_matrix_radians;
	using numeric::z_rotation_matrix_radians;
	using numeric::constants::d::pi;
	using namespace core::kinematics;

	Stub stub(c,b,a);
	xyzMatrix M(stub.M * x_rotation_matrix_radians( phi ));
	M *= z_rotation_matrix_radians( pi - theta );
	//M *= z_rotation_matrix_radians( theta );
	d = stub.v + distance * M.col_x();
}

core::Real BBConRotMover::calc_jacobian_cartesians(
	xyzVector const &v6,
	xyzVector const &v7,
	xyzVector const &v8,
	xyzVector const &v9,
	xyzVector const &v10
)
{
	double A[5][5];

	xyzVector u1(v7-v6);
	u1.normalize();
	xyzVector u2(v8-v7);
	u2.normalize();
	xyzVector u3(v9-v8);
	u3.normalize();
	xyzVector u4(v10-v9);
	u4.normalize();

	xyzVector s1(u2.cross(u1));
	s1.normalize();
	xyzVector s2(u3.cross(u2));
	s2.normalize();
	xyzVector s3(u4.cross(u3));
	s3.normalize();

	xyzVector r42(v9-v7);
	xyzVector r43(v9-v8);

	//build matrix
	xyzVector b;

	b = (u1.cross(r42));
	A[0][0] = b.x();
	A[1][0] = b.y();
	A[2][0] = b.z();
	b = (u2.cross(r43));
	A[0][1] = b.x();
	A[1][1] = b.y();
	A[2][1] = b.z();
	b = (s1.cross(r42));
	A[0][2] = b.x();
	A[1][2] = b.y();
	A[2][2] = b.z();
	b = (s2.cross(r43));
	A[0][3] = b.x();
	A[1][3] = b.y();
	A[2][3] = b.z();
	// last element is 0
	A[0][4] = 0;
	A[1][4] = 0;
	A[2][4] = 0;

	// lines 4, 5 of matrix
	b = (u1.cross(u4));
	A[3][0] = b.x();
	A[4][0] = b.y();
	b = (u2.cross(u4));
	A[3][1] = b.x();
	A[4][1] = b.y();
	b = (s1.cross(u4));
	A[3][2] = b.x();
	A[4][2] = b.y();
	b = (s2.cross(u4));
	A[3][3] = b.x();
	A[4][3] = b.y();
	b= (s3.cross(u4));
	A[3][4] = b.x();
	A[4][4] = b.y();

	double det;
	Real jacobian;
	if ( get_determinant( A, 5, det ) ) {
		jacobian = 1.0 / std::fabs(det);
	} else {
		jacobian = -1.0;
	}
	return jacobian;
}

bool BBConRotMover::closure(
	//before closure
	xyzVector &r0,
	xyzVector &r1,
	xyzVector &r2,
	xyzVector &r3,
	xyzVector &r4,
	xyzVector &r5,
	xyzVector &r6,
	//after closure
	xyzVector &p1,
	xyzVector &p2,
	xyzVector &p3,
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
)
{
	//0 [ 1 2 3 4 (5) ] 6
	//move 'C' of res4, close 4=5
	using numeric::x_rotation_matrix_radians;
	using numeric::z_rotation_matrix_radians;
	using basic::periodic_range;
	using numeric::constants::d::pi;
	using numeric::constants::d::pi_2;

	Real sinA1, cosA1, sinA2, cosA2;
	Real sinO1, cosO1, /*sinO2,*/ cosO2;
	//Real a[3],b[3],c[3],d[3],k[3],j[3],n[3],v[3];
	xyzVector a, b, c, d, k, j, n, v;
	Real k_2, j_2, v_2, va;
	Real w, w_2, h, h_2, g, p1_2, p2_2, p3_2;
	//Real R2[3][3], R2t[3][3];
	//Real T0[3][3], T0t[3][3], T1[3][3], T1t[3][3], T2[3][3], T2t[3][3];
	//xyzMatrix R2, R2t;
	//xyzMatrix T0, T0t, T1, T1t, T2, T2t;
	static Real const MAXDIH2 = 50.0*pi/180.0;
	static Real const MAXANG2 = 20.0*pi/180.0;

	bool crfailed = false;

	p1_2 = p1.x()*p1.x();
	p2_2 = p2.x()*p2.x();
	p3_2 = p3.x()*p3.x();

	/* construct T0 */
	a = (r2 - r1).normalize();
	b = (r1 - r0).normalize();
	c = a.cross(b);
	g = a.dot(b);
	g = std::sqrt( 1.0 - g*g );
	c *= 1.0/g;
	d = c.cross(a);
	xyzMatrix T0;
	T0.col_x(a).col_y(d).col_z(c);
	xyzMatrix T0t(T0);
	T0t.transpose();

	k = r4 - r2;  k_2 = k.length_squared();
	j = r4 - r1;  j_2 = j.length_squared();

	cosA2 = (k_2 - p2_2 - p3_2) / (2.0*p2.x()*p3.x());
	if ( std::fabs(cosA2)>1.0 ) {
		//bad case
		TR.Debug <<"can't close" << std::endl;
		return true;
	}
	sinA2 = -std::sqrt( 1.0 - cosA2*cosA2 );
	//set_rotation_matrix_Z( cosA2, sinA2, T2, T2t );
	//sin -> minus?
	//TR << "cos: " << cosA2 << std::endl;
	//TR << "sin: " << sinA2 << std::endl;
	//TR << "atan2: " << std::atan2(sinA2, cosA2)*180/pi << std::endl;
	//TR << "acos:" << std::acos(cosA2)*180/pi << std::endl;
	xyzMatrix T2(z_rotation_matrix_radians( std::atan2(sinA2, cosA2) ));
	a = T2 * p3;
	k = a + p2;
	cosO2 = cos( theta2_old );
	//sinO2 = sin( theta2_old );  // set but never used ~Labonte
	//set_rotation_matrix_X( cosO2, sinO2, R2, R2t );
	xyzMatrix R2(x_rotation_matrix_radians(theta2_old));
	v = R2 * k;
	v_2 = v.length_squared();

	w   = (j_2 - p1_2 - v_2) / (2.0 * p1.x());
	w_2 = w*w;
	va  = v.x()*v.x() + v.y()*v.y();
	h_2 = va - w_2;

	if ( h_2 < 0 ) {
		TR.Debug << "h_2=" << h_2 << std::endl;
		crfailed = true;

		TR.Debug <<"can't close" << std::endl;
	} else {
		// do 1st branch
		h = std::sqrt(h_2);
		cosA1 =  ( w * v.x() - v.y() * h ) / va;
		sinA1 = -( w * v.y() + v.x() * h ) / va;
		//set_rotation_matrix_Z( cosA1, sinA1, T1, T1t );
		xyzMatrix T1(z_rotation_matrix_radians( std::atan2(sinA1, cosA1) ));

		a = T1 * v;
		j = a + p1;

		a = r4 - r1;
		n = T0t * a;

		cosO1 = ( j.z()*n.z() + j.y()*n.y() ) / (j.y()*j.y() + j.z()*j.z() );
		sinO1 = ( j.y()*n.z() - j.z()*n.y() ) / (j.y()*j.y() + j.z()*j.z() );

		alpha1 = std::atan2( sinA1, cosA1 );
		alpha1 += alpha1<0 ? pi : 0;
		theta1 = std::atan2( sinO1, cosO1 );
		alpha2 = std::atan2( sinA2, cosA2 );
		alpha2 += alpha2<0 ? pi : 0;
		theta2 = theta2_old;

		if ( alpha1 != alpha1 || theta1 != theta1 ) {
			//something wrong of param calculation
			TR.Debug << "ERROR COS SIN" << std::endl;
			TR.Debug << "DB: A1 " << sinA1 << " " << cosA1 << std::endl;
			TR.Debug << "DB: A2 " << sinA2 << " " << cosA2 << std::endl;
			TR.Debug << "DB: O1 " << sinO1 << " " << cosO2 << std::endl;
		}

		TR.Debug << "alpha1=" << alpha1*180.0/pi << " alpha2=" << alpha2*180.0/pi << std::endl;
		//calc_cartesian2( p2[0], *alpha1, *theta1, r0, r1, r2, r3 );
		get_xyz(r0, r1, r2, r3, p2.x(), alpha1, theta1);
		alpha3 = numeric::angle_radians( r3, r4, r5 );
		theta3 = numeric::dihedral_radians( r2, r3, r4, r5 );
		theta4 = numeric::dihedral_radians( r3, r4, r5, r6 );

		//TR << "dihs: " << theta1 << " " << theta2 << " " << theta3 << " " << theta4 << std::endl;
		//TR << "angs: " << alpha1 << " " << alpha2 << " " << alpha3 << std::endl;

		if (
				(std::fabs(periodic_range( theta1_old-theta1, pi_2 )) > MAXDIH2) ||
				(std::fabs(periodic_range( theta2_old-theta2, pi_2 )) > MAXDIH2) ||
				(std::fabs(periodic_range( theta3_old-theta3, pi_2 )) > MAXDIH2) ||
				(std::fabs(periodic_range( theta4_old-theta4, pi_2 )) > MAXDIH2) ||
				(std::fabs(periodic_range( alpha1_old-alpha1, pi_2 )) > MAXANG2) ||
				(std::fabs(periodic_range( alpha2_old-alpha2, pi_2 )) > MAXANG2) ||
				(std::fabs(periodic_range( alpha3_old-alpha3, pi_2 )) > MAXANG2) ||
				alpha1<0.017 || alpha1>3.124
				) {
			//printf("in 2nd branch\n");
			// do 2nd branch
			h = std::sqrt(h_2);
			cosA1 =  ( w * v.x() + v.y() * h ) / va;
			sinA1 = -( w * v.y() - v.x() * h ) / va; //+/-

			//set_rotation_matrix_Z( cosA1, sinA1, T1, T1t );
			xyzMatrix T1(z_rotation_matrix_radians( std::atan2(sinA1, cosA1) ));
			a = T1 * v;
			j = a + p1;

			a = r4 - r1;
			n = T0t * a;
			cosO1 = ( j.z()*n.z() + j.y()*n.y() ) / (j.y()*j.y() + j.z()*j.z() );
			sinO1 = ( j.y()*n.z() - j.z()*n.y() ) / (j.y()*j.y() + j.z()*j.z() );

			alpha1 = std::atan2( sinA1, cosA1 );
			alpha1 += alpha1<0 ? pi : 0;
			theta1 = std::atan2( sinO1, cosO1 );
			alpha2 = std::atan2( sinA2, cosA2 );
			alpha2 += alpha2<0 ? pi : 0;
			theta2 = theta2_old;

			TR.Debug << "alpha1=" << alpha1*180.0/pi << " alpha2=" << alpha2*180.0/pi << std::endl;

			//calc_cartesian2( p2.x(), *alpha1, *theta1, r0, r1, r2, r3 );
			get_xyz(r0, r1, r2, r3, p2.x(), alpha1, theta1);
			alpha3 = numeric::angle_radians( r3, r4, r5 );
			theta3 = numeric::dihedral_radians( r2, r3, r4, r5 );
			theta4 = numeric::dihedral_radians( r3, r4, r5, r6 );

			//TR << "dihs: " << theta1 << " " << theta2 << " " << theta3 << " " << theta4 << std::endl;
			//TR << "angs: " << alpha1 << " " << alpha2 << " " << alpha3 << std::endl;

			if (
					(std::fabs(periodic_range( theta1_old-theta1, pi_2 )) > MAXDIH2) ||
					(std::fabs(periodic_range( theta2_old-theta2, pi_2 )) > MAXDIH2) ||
					(std::fabs(periodic_range( theta3_old-theta3, pi_2 )) > MAXDIH2) ||
					(std::fabs(periodic_range( theta4_old-theta4, pi_2 )) > MAXDIH2) ||
					(std::fabs(periodic_range( alpha1_old-alpha1, pi_2 )) > MAXANG2) ||
					(std::fabs(periodic_range( alpha2_old-alpha2, pi_2 )) > MAXANG2) ||
					(std::fabs(periodic_range( alpha3_old-alpha3, pi_2 )) > MAXANG2) ||
					alpha1<0.017 || alpha1>3.124
					) {
				crfailed = true;
			}
		}
	}

	return crfailed;
}

bool BBConRotMover::get_determinant(double a[5][5], int n, double &d)
{
	int i,imax=0,j,k;
	double big,dum,sum,temp;
	double vv[10];
	//int indx[10];

	d=1.0;
	for ( i=0; i<n; i++ ) {
		big=0.0;
		for ( j=0; j<n; j++ ) {
			if ( (temp=fabs(a[i][j])) > big ) big=temp;
		}
		if ( big == 0.0 ) {
			//printf("Singular matrix in routine get_determinant\n");
			//exit(1);
			TR.Warning << "Singular matrix in routine get_determinant" << std::endl;
			return false;
		}
		vv[i]=1.0/big;
	}
	for ( j=0; j<n; j++ ) {
		for ( i=0; i<j; i++ ) {
			sum=a[i][j];
			for ( k=0; k<i; k++ ) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for ( i=j; i<n; i++ ) {
			sum=a[i][j];
			for ( k=0; k<j; k++ ) {
				sum -= a[i][k]*a[k][j];
			}
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big ) {
				big=dum;
				imax=i;
			}
		}
		if ( j != imax ) {
			for ( k=0; k<n; k++ ) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			d = -(d);
			vv[imax]=vv[j];
		}
		//indx[j]=imax;  // set but never used ~Labonte
		if ( a[j][j] == 0.0 ) a[j][j]=1.0e-20;
		if ( j != (n-1) ) {
			dum=1.0/(a[j][j]);
			for ( i=j+1; i<n; i++ ) a[i][j] *= dum;
		}
	}

	for ( j = 0; j < n; j++ )  d *= a[j][j];

	return(true);
}

}
}

