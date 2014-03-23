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

#include <protocols/simple_moves/BBGaussianMover.hh>
#include <protocols/simple_moves/BBGaussianMoverCreator.hh>

//core
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <basic/basic.hh>

#include <core/kinematics/MoveMap.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/bbg.OptionKeys.gen.hh>

//Rosetta Scripts
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <protocols/rosetta_scripts/util.hh>

//util
#include <utility/exit.hh>
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>

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

static basic::Tracer TR("protocols.simple_moves.BBGaussianMover");
static numeric::random::RandomGenerator RG(6233); //Magic Number

namespace protocols {
namespace simple_moves {

std::string
BBGaussianMoverCreator::keyname() const
{
	return BBGaussianMoverCreator::mover_name();
}

protocols::moves::MoverOP
BBGaussianMoverCreator::create_mover() const {
	return new BBG8T3AMover;
}

std::string
BBGaussianMoverCreator::mover_name()
{
	return "BBG8T3AMover";
}

///////////////////////////////////////////////////
BBGaussianMover::BBGaussianMover( Size n_end_atom, Size n_dof_angle, Size n_pert_res)
        :Mover(),
        n_end_atom_(n_end_atom),
        n_dof_angle_(n_dof_angle),
        n_pert_res_(n_pert_res),
        matrix_G(n_dof_angle_,utility::vector1<Real>(n_dof_angle_)),
        matrix_A(n_dof_angle_,utility::vector1<Real>(n_dof_angle_)),
        matrix_dRdPhi(n_end_atom_,utility::vector1<Vector>(n_dof_angle_)),
        last_proposal_density_ratio_(1.0)
{
    protocols::moves::Mover::type("BBGaussianMover");
		core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
		movemap->set_bb(true);
		movemap_=movemap;
		//    setup_list(pose);
}

BBGaussianMover::~BBGaussianMover(){}

core::Real BBGaussianMover::last_proposal_density_ratio() const
{
    return last_proposal_density_ratio_;
}

void BBGaussianMover::setup_list(Pose const &pose)
{
    using namespace id;
    using namespace basic::options;
    using namespace basic::options::OptionKeys;

    //!!! can not handle cutpoint
    available_seg_list_.erase(available_seg_list_.begin(),available_seg_list_.end());
    int seg_len = 0;
    for (int n=1,end=static_cast< int >(pose.n_residue()); n<=end; n++)
    {
        conformation::Residue const & rsd( pose.residue( n ) );
        if ( rsd.is_protein()
            && movemap_->get( TorsionID( n, BB, phi_torsion ) )
            && movemap_->get( TorsionID( n, BB, psi_torsion ) )
            && n<end
            && (!option[bbg::ignore_improper_res] || pose.residue(n).name1()!='P')  )
        {
            //skip Proline
            seg_len++;
        }
        else
        {
            //TR << "No: " << n << std::endl;
            if (seg_len>0)
            {
                //put them in
                int first(n-seg_len);
                int last;
                if( n==end &&
                	 pose.residue( n ).is_protein()
                	 && movemap_->get( TorsionID( n, BB, phi_torsion ) )
                	 && movemap_->get( TorsionID( n, BB, psi_torsion ) )
                	 && (!option[bbg::ignore_improper_res] || pose.residue(n).name1()!='P') )
                {
                	last = end;
                }
                else
                {
                	last = n-1;
                }
                //TR << "first=" << first << " last=" << last << std::endl;

                int r_num(first);
                int l_num(first + 1 - static_cast< int >(n_pert_res_));
                //TR << "l=" << l_num << " r=" << r_num << std::endl;
                for (;l_num<=last;r_num++,l_num++)
                {
                    //TR << "l=" << l_num << " r=" << r_num << std::endl;
                    //in the middle of the scaffold, no end fix case
                    int part1, part2;

                    if (option[bbg::fix_short_segment] && l_num<first) continue;
                    else part1 = l_num<first?first:l_num;
                    if (option[bbg::fix_short_segment] && r_num>last) continue;
                    else part2 = r_num>last?last:r_num;

                    //get one
                    std::pair<Size, Size> seg(part1,part2);
                    available_seg_list_.push_back(seg);
                }
            }
            seg_len = 0;
        }
    }
}

core::kinematics::MoveMapCOP BBGaussianMover::movemap()
{
		return movemap_;
}

void BBGaussianMover::movemap(core::kinematics::MoveMapCOP new_movemap)
{
		movemap_=new_movemap;
		available_seg_list_.clear();
}

core::Real BBGaussianMover::cholesky_fw(Matrix &a, Size n, utility::vector1<Real> &delta, utility::vector1<Real> &dphi, Size from, Size to, Real scale)
{
    Size i,j,k;
    Real sum;
    utility::vector1<Real> p(n);

    for (i=1;i<=n;i++)
    {
        for (j=i;j<=n;j++)
        {
            for (sum=a[i][j],k=i-1;k>=1;k--) sum -= a[i][k]*a[j][k];
            if (i == j)
            {
                runtime_assert(sum>0.0);
                p[i]=sqrt(sum);
            }
            else
            {
                a[j][i]=sum/p[i];
            }
        }
    }

    for(i=from; i<=to; i++)
    {
        p[i]*=scale;
        for(j=from; j<i; j++)
        {
            a[i][j] *= scale;
        }
    }

    for (i=n;i>0;i--)
    {
        for (sum=delta[i],k=i+1;k<=n;k++) sum-=a[k][i]*dphi[k];
        dphi[i]=sum/p[i];
    }

    Real detL = 1.0;
    for (i=1;i<=n;i++) detL*=p[i];
    return detL;
}

core::Real BBGaussianMover::cholesky_bw(Matrix &a, Size n, utility::vector1<Real> &dphi, utility::vector1<Real> &delta, Size from, Size to, Real scale)
{
    Size i,j,k;
    Real sum;
    utility::vector1<Real> p(n);

    for (i=1;i<=n;i++)
    {
        for (j=i;j<=n;j++)
        {
            for (sum=a[i][j],k=i-1;k>=1;k--) sum -= a[i][k]*a[j][k];
            if (i == j)
            {
                runtime_assert(sum>0.0);
                p[i]=sqrt(sum);
            }
            else
            {
                a[j][i]=sum/p[i];
            }
        }
    }

    for(i=from; i<=to; i++)
    {
        p[i]*=scale;
        for(j=from; j<i; j++)
        {
            a[i][j] *= scale;
        }
    }

    for (i=1;i<=n;i++)
    {
        delta[i]=p[i]*dphi[i];
        for (j=i+1;j<=n;j++) delta[i]+=a[j][i]*dphi[j];
    }

    Real detL = 1.0;
    for (i=1;i<=n;i++) detL*=p[i];
    return detL;
}

/// @brief randomly rotate the dih angle in this range for avoiding the fixed ends
void BBGaussianMover::pivot_range_randomly(Pose &pose, Size i, Size to)
{
    //copy from Backbone mover, make_move()
    static Real big_angle_= 12.0;//loop parameter

    for (;i<=to;i++)
    {
        Real old_phi_ = pose.phi(i);
        Real new_phi_ = basic::periodic_range( old_phi_ + RG.gaussian() * big_angle_, 360.0 );
        //Real new_phi_ = basic::periodic_range( old_phi_ + RG.uniform() * big_angle_, 360.0 );

        Real old_psi_ = pose.psi(i);
        Real new_psi_ = basic::periodic_range( old_psi_ + RG.gaussian() * big_angle_, 360.0 );
        //Real new_psi_ = basic::periodic_range( old_psi_ + RG.uniform() * big_angle_, 360.0 );

        pose.set_phi( i, new_phi_ );
        pose.set_psi( i, new_psi_ );
    }
}

////////////////////////////////////////////
// BBG8T3AMover
////////////////////////////////////////////

protocols::moves::MoverOP
BBG8T3AMover::clone() const {
		protocols::simple_moves::BBG8T3AMoverOP mp = new BBG8T3AMover();
		mp->factorA(factorA_);
		mp->factorB(factorB_);
		mp->movemap(movemap_);
		return static_cast< protocols::moves::MoverOP >(mp);
}

void
BBG8T3AMover::factorA( core::Real const fA ){
	factorA_ = fA;
}

void
BBG8T3AMover::factorB( core::Real const fB ){
	factorB_ = fB;
}

BBG8T3AMover::BBG8T3AMover()
        :BBGaussianMover(3,8,4),
        dphi(utility::vector1<Real>(n_dof_angle_))
{
    protocols::moves::Mover::type("BBG8T3AMover");
    //build end atom list
    end_atom_list_.push_back("CA");
    end_atom_list_.push_back("C");
    end_atom_list_.push_back("O");

    //init the A/C factor
    factorA_ = option[ bbg::factorA ];
    factorB_ = option[ bbg::factorB ];
}

BBG8T3AMover::~BBG8T3AMover(){}

std::string
BBG8T3AMover::get_name() const {
	return "BBG8T3AMover";
}

void BBG8T3AMover::get_VdRdPhi(Pose const &pose)
{
    conformation::Residue const & rsd4( pose.residue( resnum_ ) );
    conformation::Residue const & rsd3( pose.residue( resnum_-1 ) );
    conformation::Residue const & rsd2( pose.residue( resnum_-2 ) );
    conformation::Residue const & rsd1( pose.residue( resnum_-3 ) );

    for (Size i=1;i<=end_atom_list_.size();i++)
    {
        //for each end atom
        Vector end_xyz = rsd4.atom(end_atom_list_[i]).xyz();
        //TR << "Phi: 8" << endl;
        matrix_dRdPhi[i][8] = get_dRdPhi(rsd4.atom("CA").xyz(),
                                         rsd4.atom("C").xyz(),
                                         end_xyz);
        //TR << "Phi: 7" << endl;
        matrix_dRdPhi[i][7] = get_dRdPhi(rsd4.atom("N").xyz(),
                                         rsd4.atom("CA").xyz(),
                                         end_xyz);
        //TR << "Phi: 6" << endl;
        matrix_dRdPhi[i][6] = get_dRdPhi(rsd3.atom("CA").xyz(),
                                         rsd3.atom("C").xyz(),
                                         end_xyz);
        //TR << "Phi: 5" << endl;
        matrix_dRdPhi[i][5] = get_dRdPhi(rsd3.atom("N").xyz(),
                                         rsd3.atom("CA").xyz(),
                                         end_xyz);
        //TR << "Phi: 4" << endl;
        matrix_dRdPhi[i][4] = get_dRdPhi(rsd2.atom("CA").xyz(),
                                         rsd2.atom("C").xyz(),
                                         end_xyz);
        //TR << "Phi: 3" << endl;
        matrix_dRdPhi[i][3] = get_dRdPhi(rsd2.atom("N").xyz(),
                                         rsd2.atom("CA").xyz(),
                                         end_xyz);
        //TR << "Phi: 2" << endl;
        matrix_dRdPhi[i][2] = get_dRdPhi(rsd1.atom("CA").xyz(),
                                         rsd1.atom("C").xyz(),
                                         end_xyz);
        //TR << "Phi: 1" << endl;
        matrix_dRdPhi[i][1] = get_dRdPhi(rsd1.atom("N").xyz(),
                                         rsd1.atom("CA").xyz(),
                                         end_xyz);
    }
}

void BBG8T3AMover::get_G()
{
    for (Size i=1; i<=n_dof_angle_; i++)
    {
        for (Size j=i; j<=n_dof_angle_; j++)
        {
            matrix_G[i][j] = 0.0;
            for (Size n=1; n<=end_atom_list_.size();n++)
            {
                matrix_G[i][j] += matrix_dRdPhi[n][i].dot(matrix_dRdPhi[n][j]);
            }
            //if (matrix_G[i][j]>-ZERO && matrix_G[i][j]<ZERO) matrix_G[i][j] = 0.0;
            if (i<j) matrix_G[j][i]=matrix_G[i][j];

        }
    }
}

void BBG8T3AMover::get_A()
{
    for (Size i=1; i<=n_dof_angle_; i++)
    {
        for (Size j=i; j<=n_dof_angle_; j++)
        {
            matrix_A[i][j] = factorB_ * matrix_G[i][j];
            if (i==j) matrix_A[i][j] += 1.0;
            matrix_A[i][j] *= factorA_ / 2.0;
            if (i<j) matrix_A[j][i] = matrix_A[i][j];
        }
    }
}

core::Real BBG8T3AMover::get_L_move(Pose &pose)
{
    //gerate a Gaussian dx vector
    utility::vector1<Real> delta(n_dof_angle_);
    for (Size i=1; i<=n_dof_angle_; i++) delta[i]=RG.gaussian();

    //calculate d^2 = delta^2
    Real d2=0.0;
    for (Size i=1; i<=n_dof_angle_; i++) d2+=delta[i]*delta[i];

    //cholesky, get L^t, L^-1
    Real detL = cholesky_fw(matrix_A, n_dof_angle_, delta, dphi);

    //W_old *= exp(-d^2)
    Real W_old = detL*exp(-d2/2.0);

    //set the new phi,psi (above all called phi, actually 4 phi, 4 psi)
    pose.set_psi(resnum_, basic::periodic_range( pose.psi(resnum_)+dphi[8], 360.0 ) );
    pose.set_phi(resnum_, basic::periodic_range( pose.phi(resnum_)+dphi[7], 360.0 ) );
    pose.set_psi(resnum_-1, basic::periodic_range( pose.psi(resnum_-1)+dphi[6], 360.0 ) ) ;
    pose.set_phi(resnum_-1, basic::periodic_range( pose.phi(resnum_-1)+dphi[5], 360.0 ) );
    pose.set_psi(resnum_-2, basic::periodic_range( pose.psi(resnum_-2)+dphi[4], 360.0 ) );
    pose.set_phi(resnum_-2, basic::periodic_range( pose.phi(resnum_-2)+dphi[3], 360.0 ) );
    pose.set_psi(resnum_-3, basic::periodic_range( pose.psi(resnum_-3)+dphi[2], 360.0 ) );
    pose.set_phi(resnum_-3, basic::periodic_range( pose.phi(resnum_-3)+dphi[1], 360.0 ) );

    return W_old;
}

core::Real BBG8T3AMover::get_L_prime()
{
    utility::vector1<Real> delta(n_dof_angle_);
    //get L
    Real detL = cholesky_bw(matrix_A, n_dof_angle_, dphi, delta);
    
    //calculate d^2 = delta^2
    Real d2=0.0;
    for (Size i=1; i<=n_dof_angle_; i++)d2+=delta[i]*delta[i];

    Real W_new = detL*exp(-d2/2.0);

    return W_new;
}

void BBG8T3AMover::apply(Pose &pose)
{
    //setup_list(pose);
    if(available_seg_list_.size()==0)setup_list(pose);
    if(available_seg_list_.size()==0)return;
    //randomly select a residue in the list
    int ndx=static_cast< int >( RG.uniform()*available_seg_list_.size()+1 );

    Size left = available_seg_list_[ ndx ].first;
    resnum_ = available_seg_list_[ ndx ].second;

    if (resnum_-left+1 < n_pert_res_)
    {
        //do random
        pivot_range_randomly(pose, left, resnum_);
        return;
    }

    get_VdRdPhi(pose);
    get_G();
    get_A();
    Real W_old = get_L_move(pose);

    //TR << endl;
    //TR << "------new-------" << endl;
    get_VdRdPhi(pose);
    get_G();
    get_A();
    Real W_new = get_L_prime();

    last_proposal_density_ratio_ = W_new / W_old;
}

void BBG8T3AMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const & pose )
{
		factorA_ = tag->getOption< Real >("factorA", 1.0);
		factorB_ = tag->getOption< Real >("factorB", 10.0);

		core::kinematics::MoveMapOP mm = new core::kinematics::MoveMap;
		mm->set_bb(true);
		protocols::rosetta_scripts::parse_movemap( tag, pose, mm, data, false );
		movemap(mm);
}

}//namespace simple_moves
}//namespace protocols

