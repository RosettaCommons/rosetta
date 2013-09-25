/// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_moves/BBGaussianMover.fwd.hh
/// @brief  Gaussian Perturbation to backbone
/// @author Yuan Liu (wendao@u.washington.edu)

#ifndef INCLUDED_protocols_simple_moves_BBGaussianMover_hh
#define INCLUDED_protocols_simple_moves_BBGaussianMover_hh

#include <protocols/simple_moves/BBGaussianMover.fwd.hh>
//core
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
//protocols
#include <protocols/moves/Mover.hh>
//movemap
#include <core/kinematics/MoveMap.fwd.hh>
//
#include <utility/vector1.hh>

//std
#include <string>

#include <core/conformation/Residue.fwd.hh>
#include <numeric/xyzVector.hh>


using namespace std;

namespace protocols {
namespace simple_moves {

/// @brief the basic class for implementing Gaussian Perturbation of bb
/// @note "Monte Carlo update for chain molecules: Biased Gaussian steps in torsional space"
/// "The Journal of Chemical Physics, Vol. 114, No. 18. (2001), pp. 8154-8158."
/// Two steps for perturbing the backbone and keeping the geometry constrain
/// Step 1: Gaussian Biased steps in torsional space:
/// the end atoms of the moving segment should be choosen by user
/// and to keep the geometry constrain(6), the DOF of the moving segment > 6
/// Step 2: "pivot" update the bb conformation, or followed by a chainclosure method (BBConRot)
class BBGaussianMover : public moves::Mover
{
public:
    typedef core::Vector Vector;
    typedef core::Real Real;
    typedef core::Size Size;
    typedef core::pose::Pose Pose;
    typedef core::pose::PoseOP PoseOP;
    typedef core::pose::PoseCOP PoseCOP;
    typedef core::conformation::Residue Residue;
    typedef utility::vector1< Vector > VdRdPhi;
    typedef utility::vector1< VdRdPhi > VMatrix;
    typedef utility::vector1< utility::vector1<Real> > Matrix;
    typedef numeric::xyzMatrix<Real> xyzMatrix;

public:
	BBGaussianMover( Size, Size, Size);
	~BBGaussianMover();

	//go through all the residue, count how many movable residue passed
	//count > n_pert_res (depends on the n_dof_angle)
	//stop by frozen res and cut point
	//TODO: apply smallmover to the end segment if their length is not
	//satisfy the n_pert_res_, dof: L-(n-1) => L+(n+1)

	core::kinematics::MoveMapCOP movemap();
	void movemap(core::kinematics::MoveMapCOP new_movemap);

	using moves::Mover::last_proposal_density_ratio;

	/// @brief get the ratio of proposal densities for the last move
	Real last_proposal_density_ratio() const;

protected:
    void setup_list(Pose const &);
	//r1, r2, r: r rotate around axis r1->r2
	Vector get_dRdPhi(Vector const &r1, Vector const &r2, Vector const &r)
    {
        //dr/dphi = r x axis
        Vector axis((r2-r1).normalize());
        Vector r_local(r-r1);
        return r_local.cross(axis);
    }

	//r1, r2, r3, r: r rotate around axis (r3-r2)x(r2xr1)
	Vector get_dRdTheta(Vector const &r1, Vector const &r2, Vector const &r3, Vector const &r)
    {
        Vector axis((r3-r2).cross(r2-r1).normalize());
        Vector r_local(r-r2);
        return r_local.cross(axis);
    }

	//interface
	virtual void get_VdRdPhi(Pose const &)=0;
	virtual void get_G()=0;
	virtual void get_A()=0;
	//calculate the L matrix and update pose
	virtual Real get_L_move(Pose &)=0;
	//calculate the L', update last_proposal_density
	virtual Real get_L_prime()=0;

	//foreward
	Real cholesky_fw(Matrix &, Size, utility::vector1<Real> &, utility::vector1<Real> &, Size from=1, Size to=0, Real scale=1.0);
	//backward
	Real cholesky_bw(Matrix &, Size, utility::vector1<Real> &, utility::vector1<Real> &, Size from=1, Size to=0, Real scale=1.0);

	void pivot_range_randomly(Pose &, Size, Size);

protected:
	Size n_end_atom_;
	Size n_dof_angle_;
	Size n_pert_res_;
	Size resnum_;
	//utility::vector1< Size > available_res_list_;
	utility::vector1< std::pair<Size,Size> > available_seg_list_;
	/// @note {dr/dphi_i * dr/dphi_j} -- n_dof_angle^2
	Matrix matrix_G;
	/// @note perturbation from G
	Matrix matrix_A;
	/// @note {dr_i/dphi_j} -- n_dof_angle*n_end_atom
	VMatrix matrix_dRdPhi;

	//to avoid -0.0
	//static const Real ZERO;

	// proposal density ratio
	Real last_proposal_density_ratio_;

	//movemap
	core::kinematics::MoveMapCOP movemap_;
};

/// @brief a particular gaussian mover from the original paper
/// @note using 8 torsion angles as DOF, 3 atoms (Ca,C,O) as end
class BBG8T3AMover : public BBGaussianMover
{
public:
	virtual protocols::moves::MoverOP clone() const;

public:
	BBG8T3AMover();
	~BBG8T3AMover();


	//static void register_options();

	void apply(Pose &);
	virtual std::string get_name() const;
	void factorA( core::Real const fA );
	void factorB( core::Real const fB );
	//Real get_last_delta_square();
	
	//interface for rosetta_scripts
	virtual void parse_my_tag(
		TagPtr const,
		protocols::moves::DataMap &,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & );

protected:
	void get_VdRdPhi(Pose const &);
	void get_G();
	void get_A();

	Real get_L_move(Pose &);
	Real get_L_prime();

private:
	utility::vector1< string > end_atom_list_;
	utility::vector1< Real > dphi;
	Real factorA_;
	Real factorB_;
	Real last_delta_square_;
};

}//namespace simple_moves
}//namespace protocols

#endif

