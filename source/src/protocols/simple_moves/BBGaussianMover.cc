// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_moves/BBGaussianMover.fwd.hh
/// @brief  Gaussian Perturbation to backbone
/// @author Yuan Liu (wendao@u.washington.edu)

// Unit headers
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

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.BBGaussianMover" );

namespace protocols {
namespace simple_moves {

std::string
BBGaussianMoverCreator::keyname() const
{
	return BBGaussianMoverCreator::mover_name();
}

protocols::moves::MoverOP
BBGaussianMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new BBGaussianMover );
}

std::string
BBGaussianMoverCreator::mover_name()
{
	return "BBGaussian";
}

///////////////////////////////////////////////////
BBGaussianMover::BBGaussianMover()
:protocols::canonical_sampling::ThermodynamicMover()
{
	init();
}

BBGaussianMover::BBGaussianMover(Size n_end_atom, Size n_dof_angle, Size n_pert_res)
:protocols::canonical_sampling::ThermodynamicMover(),
	n_end_atom_(n_end_atom),
	n_dof_angle_(n_dof_angle),
	n_pert_res_(n_pert_res),
	matrix_G(n_dof_angle_, Vector(n_dof_angle_)),
	matrix_A(n_dof_angle_, Vector(n_dof_angle_)),
	matrix_dRdPhi(n_end_atom_, VdRdPhi(n_dof_angle_)),
	last_proposal_density_ratio_(1.0),
	dphi(n_dof_angle_)
{
	init();
	//no need to do it here, may assign movemap later
	//setup_list(pose);
}

BBGaussianMover::~BBGaussianMover(){
	TR << "Final factorA=" << factorA_ << std::endl;
}

protocols::moves::MoverOP
BBGaussianMover::clone() const {
	protocols::simple_moves::BBGaussianMoverOP mp( new BBGaussianMover() );
	mp->factorA(factorA_);
	mp->factorB(factorB_);
	mp->movemap(movemap_);
	return static_cast< protocols::moves::MoverOP >(mp);
}

std::string BBGaussianMover::get_name() const {
	return "BBGaussianMover";
}

void BBGaussianMover::resize(Size n_end_atom, Size n_dof_angle, Size n_pert_res)
{
	n_end_atom_ = n_end_atom;
	n_dof_angle_ = n_dof_angle;
	n_pert_res_ = n_pert_res;

	matrix_G.resize(n_dof_angle, Vector(n_dof_angle));
	matrix_A.resize(n_dof_angle, Vector(n_dof_angle));
	matrix_dRdPhi.resize(n_end_atom, VdRdPhi(n_dof_angle));
	dphi.resize(n_dof_angle);
}

void BBGaussianMover::init()
{
	protocols::moves::Mover::type("BBGaussianMover");
	//default fix everything
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
	movemap_=movemap;

	//build end atom list
	end_atom_list_.push_back(std::pair<Size, std::string>(0,"CA"));
	end_atom_list_.push_back(std::pair<Size, std::string>(0,"C"));
	end_atom_list_.push_back(std::pair<Size, std::string>(0,"O"));
	n_end_atom_ = 3;

	//init the A/C factor
	factorA_ = option[ bbg::factorA ];
	factorB_ = option[ bbg::factorB ];

	last_proposal_density_ratio_ = 1.0;
	preserve_detailed_balance_ = true;

	fix_short_segment_ = option[bbg::fix_short_segment];
	use_all_pivot_res_ = false;
	shrink_frag_ends_ = false;
	auto_adjust_factorA_ = false;

	N_auto_all = 0;
	N_auto_small = 0;
}

void BBGaussianMover::setup_list(Pose const &pose)
{
	using namespace id;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//!!! can not handle cutpoint
	available_seg_list_.erase(available_seg_list_.begin(),available_seg_list_.end());

	int seg_len = 0;
	if ( use_all_pivot_res_ ) {
		for ( int n=1,end=static_cast< int >(pose.size()); n<=end; n++ ) {
			conformation::Residue const & rsd( pose.residue( n ) );
			if ( rsd.is_protein() && movemap_->get_bb(n) && n<end
					&& (!option[bbg::ignore_improper_res] || pose.residue(n).name1()!='P')  ) {
				seg_len++;
			} else {
				if ( seg_len>0 ) {
					int first(n-seg_len);
					int last;
					if ( n==end &&
							pose.residue( n ).is_protein()
							&& movemap_->get_bb( n )
							&& (!option[bbg::ignore_improper_res] || pose.residue(n).name1()!='P') ) {
						last = end;
					} else {
						last = n-1;
					}
					std::pair<Size, Size> seg(first, last);
					available_seg_list_.push_back(seg);
				}
				seg_len = 0;
			}
		}
	} else {
		for ( int n=1,end=static_cast< int >(pose.size()); n<=end; n++ ) {
			conformation::Residue const & rsd( pose.residue( n ) );
			if ( rsd.is_protein() && movemap_->get_bb(n) && n<end
					&& (!option[bbg::ignore_improper_res] || pose.residue(n).name1()!='P')  ) {
				//skip Proline
				seg_len++;
			} else {
				//TR << "No: " << n << std::endl;
				if ( seg_len>0 ) {
					//put them in
					int first(n-seg_len);
					int last;
					if ( n==end &&
							pose.residue( n ).is_protein()
							&& movemap_->get_bb( n )
							&& (!option[bbg::ignore_improper_res] || pose.residue(n).name1()!='P') ) {
						last = end;
					} else {
						last = n-1;
					}
					//TR << "first=" << first << " last=" << last << std::endl;

					int r_num(first);
					int l_num(first + 1 - static_cast< int >(n_pert_res_));
					//TR << "l=" << l_num << " r=" << r_num << std::endl;
					for ( ; l_num<=last; r_num++,l_num++ ) {
						//TR << "l=" << l_num << " r=" << r_num << std::endl;
						//in the middle of the scaffold, no end fix case
						int part1, part2;

						if ( fix_short_segment_ && l_num<first ) {
							//short frag
							continue;
						} else {
							if ( l_num<first ) {
								//short
								if ( first==1 ) {
									//only for tail
									part1 = first;
								} else {
									continue;
								}
							} else {
								//normal
								part1 = l_num;
							}
						}
						if ( fix_short_segment_ && r_num>last ) {
							continue;
						} else {
							if ( r_num>last ) {
								if ( last==end ) {
									part2 = last;
								} else {
									continue;
								}
							} else {
								part2 = r_num;
							}
						}

						//get one
						std::pair<Size, Size> seg(part1,part2);
						available_seg_list_.push_back(seg);
					}
				}
				seg_len = 0;
			}
		}
	}

	Size nseg = available_seg_list_.size();
	//check logic
	if ( use_all_pivot_res_ ) {
		if ( n_dof_angle_ != (n_pert_res_-nseg)*2 + (shrink_frag_ends_ ? 0 : 2*nseg) ) {
			TR.Debug << "DOF=" << n_dof_angle_ << " N_RES=" << n_pert_res_ << " N_SEG=" << nseg << std::endl;
			throw utility::excn::EXCN_RosettaScriptsOption( " BBG: dof and res number mismatch" );
		}
	}
	//debug info
	for ( Size i=1; i<=nseg; i++ ) {
		TR.Debug << "Seg " << i << ": " << available_seg_list_[i].first << "<-->" << available_seg_list_[i].second << std::endl;
	}
}

core::kinematics::MoveMapCOP BBGaussianMover::movemap()
{
	return movemap_;
}

void BBGaussianMover::movemap(core::kinematics::MoveMapCOP new_movemap)
{
	movemap_=new_movemap;
	available_seg_list_.erase(available_seg_list_.begin(),available_seg_list_.end());
}

core::Real BBGaussianMover::cholesky_fw(Matrix &a, Size n, Vector &delta, Vector &dphi, Size from, Size to, Real scale)
{
	Size i,j,k;
	Real sum;
	Vector p(n);

	for ( i=1; i<=n; i++ ) {
		for ( j=i; j<=n; j++ ) {
			for ( sum=a[i][j],k=i-1; k>=1; k-- ) sum -= a[i][k]*a[j][k];
			if ( i == j ) {
				if ( sum>0.0 ) {
					p[i]=sqrt(sum);
				} else {
					throw utility::excn::EXCN_RosettaScriptsOption( " BBG: Cholesky decomposition failed, may be a wrong structure" );
				}
			} else {
				a[j][i]=sum/p[i];
			}
		}
	}

	for ( i=from; i<=to; i++ ) {
		p[i]*=scale;
		for ( j=from; j<i; j++ ) {
			a[i][j] *= scale;
		}
	}

	for ( i=n; i>0; i-- ) {
		for ( sum=delta[i],k=i+1; k<=n; k++ ) sum-=a[k][i]*dphi[k];
		dphi[i]=sum/p[i];
	}

	Real detL = 1.0;
	for ( i=1; i<=n; i++ ) detL*=p[i];
	return detL;
}

core::Real BBGaussianMover::cholesky_bw(Matrix &a, Size n, Vector &dphi, Vector &delta, Size from, Size to, Real scale)
{
	Size i,j,k;
	Real sum;
	Vector p(n);

	for ( i=1; i<=n; i++ ) {
		for ( j=i; j<=n; j++ ) {
			for ( sum=a[i][j],k=i-1; k>=1; k-- ) sum -= a[i][k]*a[j][k];
			if ( i == j ) {
				if ( sum>0.0 ) {
					p[i]=sqrt(sum);
				} else {
					throw utility::excn::EXCN_RosettaScriptsOption( " BBG: Cholesky decomposition failed, may be a wrong structure" );
				}
			} else {
				a[j][i]=sum/p[i];
			}
		}
	}

	for ( i=from; i<=to; i++ ) {
		p[i]*=scale;
		for ( j=from; j<i; j++ ) {
			a[i][j] *= scale;
		}
	}

	for ( i=1; i<=n; i++ ) {
		delta[i]=p[i]*dphi[i];
		for ( j=i+1; j<=n; j++ ) delta[i]+=a[j][i]*dphi[j];
	}

	Real detL = 1.0;
	for ( i=1; i<=n; i++ ) detL*=p[i];
	return detL;
}

/// @brief randomly rotate the dih angle in this range for avoiding the fixed ends
void BBGaussianMover::pivot_range_randomly(Pose &pose, Size i, Size to)
{
	//copy from Backbone mover, make_move()
	static Real big_angle_= 12.0;//loop parameter

	for ( ; i<=to; i++ ) {
		Real old_phi_ = pose.phi(i);
		Real new_phi_ = basic::periodic_range( old_phi_ + numeric::random::rg().gaussian() * big_angle_, 360.0 );
		//Real new_phi_ = basic::periodic_range( old_phi_ + numeric::random::rg().uniform() * big_angle_, 360.0 );

		Real old_psi_ = pose.psi(i);
		Real new_psi_ = basic::periodic_range( old_psi_ + numeric::random::rg().gaussian() * big_angle_, 360.0 );
		//Real new_psi_ = basic::periodic_range( old_psi_ + numeric::random::rg().uniform() * big_angle_, 360.0 );

		pose.set_phi( i, new_phi_ );
		pose.set_psi( i, new_psi_ );
	}
}

void
BBGaussianMover::factorA( core::Real const fA ){
	factorA_ = fA;
}

void
BBGaussianMover::factorB( core::Real const fB ){
	factorB_ = fB;
}

void BBGaussianMover::get_G()
{
	for ( Size i=1; i<=n_dof_angle_; i++ ) {
		for ( Size j=i; j<=n_dof_angle_; j++ ) {
			matrix_G[i][j] = 0.0;
			for ( Size n=1; n<=end_atom_list_.size(); n++ ) {
				matrix_G[i][j] += matrix_dRdPhi[n][i].dot(matrix_dRdPhi[n][j]);
			}
			//if (matrix_G[i][j]>-ZERO && matrix_G[i][j]<ZERO) matrix_G[i][j] = 0.0;
			if ( i<j ) matrix_G[j][i]=matrix_G[i][j];

		}
	}
}

void BBGaussianMover::get_A()
{
	for ( Size i=1; i<=n_dof_angle_; i++ ) {
		for ( Size j=i; j<=n_dof_angle_; j++ ) {
			matrix_A[i][j] = factorB_ * matrix_G[i][j];
			if ( i==j ) matrix_A[i][j] += 1.0;
			matrix_A[i][j] *= factorA_ / 2.0;
			if ( i<j ) matrix_A[j][i] = matrix_A[i][j];
		}
	}
}

core::Real BBGaussianMover::get_L_prime()
{
	Vector delta(n_dof_angle_);
	//get L
	Real detL = cholesky_bw(matrix_A, n_dof_angle_, dphi, delta);

	//calculate d^2 = delta^2
	Real d2=0.0;
	for ( Size i=1; i<=n_dof_angle_; i++ ) d2+=delta[i]*delta[i];

	Real W_new = detL*exp(-d2/2.0);

	return W_new;
}

core::Real BBGaussianMover::get_L_move(Pose &pose)
{
	//gerate a Gaussian dx vector
	Vector delta(n_dof_angle_);
	for ( Size i=1; i<=n_dof_angle_; i++ ) delta[i] = numeric::random::rg().gaussian();

	//calculate d^2 = delta^2
	Real d2=0.0;
	for ( Size i=1; i<=n_dof_angle_; i++ ) d2 += delta[i]*delta[i];

	//cholesky, get L^t, L^-1
	Real detL = cholesky_fw(matrix_A, n_dof_angle_, delta, dphi);

	//W_old *= exp(-d^2)
	Real W_old = detL*exp(-d2/2.0);

	Size n_dof=n_dof_angle_;
	if ( use_all_pivot_res_ ) {
		//go through each segment
		for ( Size seg=1, eseg=available_seg_list_.size(); seg<=eseg; seg++ ) {
			Size left = available_seg_list_[seg].first;
			Size right = available_seg_list_[seg].second;
			for ( Size j=right; j>=left; j-- ) {
				//check count
				runtime_assert(n_dof>0);
				//conformation::Residue const & rsd( pose.residue( j ) );
				if ( j<right || (!shrink_frag_ends_) ) {
					//skip j==right, unless shrink==false
					pose.set_psi(j, basic::periodic_range( pose.psi(j)+dphi[n_dof--], 360.0 ) );
				}
				if ( j>left || (!shrink_frag_ends_) ) {
					//skip j==left, unless shring==false
					pose.set_phi(j, basic::periodic_range( pose.phi(j)+dphi[n_dof--], 360.0 ) );
				}
			}
		}
	} else {
		//set the new phi,psi (above all called phi, actually 4 phi, 4 psi)
		for ( Size j=0; j<n_pert_res_; j++ ) {
			Size ndx = resnum_-j;
			pose.set_psi(ndx, basic::periodic_range( pose.psi(ndx)+dphi[n_dof--], 360.0 ) );
			pose.set_phi(ndx, basic::periodic_range( pose.phi(ndx)+dphi[n_dof--], 360.0 ) );
		}
	}
	runtime_assert(n_dof==0);

	return W_old;
}

void BBGaussianMover::get_VdRdPhi(Pose const &pose)
{
	for ( Size i=1; i<=end_atom_list_.size(); i++ ) {
		//for each end atom
		Size lock_res = end_atom_list_[i].first;
		xyzVector end_xyz;
		if ( lock_res==0 ) {
			runtime_assert(resnum_>0);
			conformation::Residue const & rsd_end( pose.residue( resnum_ ) );
			if ( rsd_end.has(end_atom_list_[i].second) ) {
				end_xyz = rsd_end.atom(end_atom_list_[i].second).xyz();
			} else {
				throw utility::excn::EXCN_RosettaScriptsOption( " BBG: No lock atom found on res (end)" );
			}
		} else {
			conformation::Residue const & rsd_end( pose.residue( end_atom_list_[i].first ) );
			if ( rsd_end.has(end_atom_list_[i].second) ) {
				end_xyz = rsd_end.atom(end_atom_list_[i].second).xyz();
			} else {
				throw utility::excn::EXCN_RosettaScriptsOption( " BBG: No lock atom found on res (lock)" );
			}
		}

		Size n_dof = n_dof_angle_;
		if ( use_all_pivot_res_ ) {
			//use_all and fix_tail
			//go through each segment
			for ( Size seg=1, eseg=available_seg_list_.size(); seg<=eseg; seg++ ) {
				Size left = available_seg_list_[seg].first;
				Size right = available_seg_list_[seg].second;
				for ( Size j=right; j>=left; j-- ) {
					//check count
					runtime_assert(n_dof>0);
					conformation::Residue const & rsd( pose.residue( j ) );
					if ( j<right || (!shrink_frag_ends_) ) {
						matrix_dRdPhi[i][n_dof--] = get_dRdPhi(rsd.atom("CA").xyz(), rsd.atom("C").xyz(), end_xyz);
					}
					if ( j>left || (!shrink_frag_ends_) ) {
						matrix_dRdPhi[i][n_dof--] = get_dRdPhi(rsd.atom("N").xyz(), rsd.atom("CA").xyz(), end_xyz);
					}
				}
			}
		} else {
			//use N, should be continous
			// old logic:
			// go through continues n
			for ( Size j=0; j<n_pert_res_; j++ ) {
				conformation::Residue const & rsd( pose.residue( resnum_ - j ) );
				matrix_dRdPhi[i][n_dof--] = get_dRdPhi(rsd.atom("CA").xyz(), rsd.atom("C").xyz(), end_xyz);
				matrix_dRdPhi[i][n_dof--] = get_dRdPhi(rsd.atom("N").xyz(), rsd.atom("CA").xyz(), end_xyz);
			}
		}
		runtime_assert(n_dof==0);
	}
}

void BBGaussianMover::apply(Pose &pose)
{
	//setup_list(pose);
	if ( available_seg_list_.size()==0 ) setup_list(pose);
	if ( available_seg_list_.size()==0 ) return;
	//randomly select a residue in the list

	// use all or not?
	// lock or not?

	if ( use_all_pivot_res_ ) {
		//left = available_seg_list_[1].first;
		//last res of the last frag
		resnum_ = available_seg_list_.back().second;
	} else {
		Size left=0;
		Size ndx = static_cast< int >( numeric::random::rg().uniform()*available_seg_list_.size()+1 );
		left = available_seg_list_[ ndx ].first;
		resnum_ = available_seg_list_[ ndx ].second;

		//for short tail
		if ( resnum_-left+1 < n_pert_res_ ) {
			//do random
			pivot_range_randomly(pose, left, resnum_);
			return;
		}
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

void BBGaussianMover::init_kic_loop(Size looplength, core::kinematics::MoveMapCOP mm)
{
	//restrain the last residue
	end_atom_list_.erase(end_atom_list_.begin(), end_atom_list_.end());
	end_atom_list_.push_back(std::pair<Size, std::string>(looplength, "N"));
	end_atom_list_.push_back(std::pair<Size, std::string>(looplength, "CA"));
	//end_atom_list_.push_back(std::pair<Size, std::string>(looplength-1, "C"));

	n_end_atom_ = 2;

	//setup dof
	//n_pert_res_ = looplength-2; // the first and the last res are anchor
	//n_dof_angle_ = n_pert_res_ * 2 - 2; // two angles for each mobile res
	n_pert_res_ = 0;
	n_dof_angle_ = 0;
	for ( Size i=1; i<=looplength; i++ ) {
		if ( mm->get_bb(i) ) {
			++n_pert_res_;
			n_dof_angle_+=2;
		}
	}
	//n_dof_angle_ -= 2;

	use_all_pivot_res_ = true; // use phi/psi for all res
	shrink_frag_ends_ = false; // don't shrink
	auto_adjust_factorA_ = true; // need to be auto
	fix_short_segment_ = true; // don't randomly pertub loop ends

	TR.Debug << "DOF=" << n_dof_angle_ << " N_RES=" << n_pert_res_ << " END=" << n_end_atom_ << std::endl;
	TR.Debug << "factorA=" << factorA_ << " factorB=" << factorB_ << std::endl;

	movemap(mm);
	//allocate the vector and matrix
	resize(end_atom_list_.size(), n_dof_angle_, n_pert_res_);
}

void BBGaussianMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const & pose )
{
	factorA_ = tag->getOption< Real >("factorA", 1.0);
	factorB_ = tag->getOption< Real >("factorB", 10.0);

	//if end_atom_list is defined by user, change it
	//otherwise it's the default 3 atoms
	if ( tag->hasOption("end_atoms") ) {
		end_atom_list_.erase(end_atom_list_.begin(), end_atom_list_.end());
		std::string const end_atoms_string( tag->getOption<std::string>("end_atoms") );
		utility::vector1<std::string> end_atoms = utility::string_split( end_atoms_string, ',' );
		Size n_lock = end_atoms.size();
		runtime_assert(n_lock%2==0);
		n_lock /= 2;
		for ( Size i=1; i<=n_lock; i++ ) {
			Size lock;
			// readin new end
			std::stringstream ss(end_atoms[i*2-1]);
			ss >> lock;
			end_atom_list_.push_back(std::pair<Size, std::string>(lock, end_atoms[i*2]));
			TR.Debug << end_atom_list_[i].first << " " << end_atom_list_[i].second << std::endl;
		}
		n_end_atom_ = n_lock;
	}

	//default value is for 4 continous res, and 8 assosiated dih
	n_dof_angle_ = tag->getOption< Size >("dof", 8);
	n_pert_res_ = tag->getOption< Size >("pivot", 4);

	use_all_pivot_res_ = tag->getOption< bool >("use_all", false);

	//the logic here is:
	//by default, shrink should be false
	//if use all res, set shrink to be true
	//unless manually set it to be false
	//a little bit confuse but no need to change the old scripts
	if ( use_all_pivot_res_ ) {
		shrink_frag_ends_ = true;
	} else {
		shrink_frag_ends_ = false;
	}

	if ( tag->hasOption("shrink") ) {
		//
		shrink_frag_ends_ = tag->getOption< bool >("shrink", false);
	}

	auto_adjust_factorA_ = tag->getOption< bool >("autoA", false);
	fix_short_segment_ = tag->getOption< bool >("fix_tail", fix_short_segment_);

	//debug info
	TR.Debug << "DOF=" << n_dof_angle_ << " N_RES=" << n_pert_res_ << " END=" << n_end_atom_ << std::endl;
	TR.Debug << "factorA=" << factorA_ << " factorB=" << factorB_ << std::endl;

	//check logic
	if ( use_all_pivot_res_ ) {
		if ( !fix_short_segment_ ) {
			throw utility::excn::EXCN_RosettaScriptsOption( " BBG: use_all should go with fix_tail" );
		}
	} else {
		if ( n_dof_angle_ != (n_pert_res_*2) ) {
			throw utility::excn::EXCN_RosettaScriptsOption( " BBG: dof and res number mismatch" );
		}
	}

	//reset movemap, default is empty (no move)
	core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
	protocols::rosetta_scripts::parse_movemap( tag, pose, mm, data, false );
	movemap(mm);

	//allocate the vector and matrix
	resize(end_atom_list_.size(), n_dof_angle_, n_pert_res_);
}

////////////// Thermo /////////////
core::Real BBGaussianMover::last_proposal_density_ratio()
{
	//TR << last_proposal_density_ratio_ << std::endl;

	if ( auto_adjust_factorA_ ) {
		update_counting_last_PDR(last_proposal_density_ratio_);
	}

	return last_proposal_density_ratio_;
}

void BBGaussianMover::update_counting_last_PDR(Real lastP)
{
	N_auto_all++;
	N_auto_small += lastP<0.6 ? 1 : 0;

	if ( N_auto_all > 200 ) {
		Real p = (Real)N_auto_small/N_auto_all;
		TR.Debug << "Psmall=" << p << std::endl;
		if ( p>0.4 ) {
			//factorA too small
			factorA_ *= 2;
			TR.Debug << "Auto increase factorA to " << factorA_ << std::endl;
		} else if ( p < 0.2 ) {
			//factorA too big
			factorA_ *= 0.667;
			TR.Debug << "Auto decrease factorA to " << factorA_ << std::endl;
		}

		N_auto_all = 0;
		N_auto_small = 0;
	}
}

void BBGaussianMover::set_preserve_detailed_balance( bool preserve_detailed_balance )
{
	preserve_detailed_balance_ = preserve_detailed_balance;
	if ( !preserve_detailed_balance ) {
		TR.Warning << "This mover should use preserve_detailed_balance=true !" << std::endl;
	}
}

bool BBGaussianMover::preserve_detailed_balance() const
{
	return preserve_detailed_balance_;
}

////////////////////////////////////////////
// BBG8T3AMover
////////////////////////////////////////////
BBG8T3AMover::BBG8T3AMover():BBGaussianMover()
{
	protocols::moves::Mover::type("BBG8T3AMover");
	core::kinematics::MoveMapOP mymm( new core::kinematics::MoveMap );
	mymm->set_bb(true); //use all bb dof
	movemap(mymm);
	resize(3,8,4);
}

protocols::moves::MoverOP
BBG8T3AMover::clone() const {
	protocols::simple_moves::BBG8T3AMoverOP mp( new BBG8T3AMover() );
	mp->factorA(factorA_);
	mp->factorB(factorB_);
	mp->movemap(movemap_);
	return static_cast< protocols::moves::MoverOP >(mp);
}

BBG8T3AMover::~BBG8T3AMover()= default;

std::string BBG8T3AMover::get_name() const {
	return "BBG8T3AMover";
}

}//namespace simple_moves
}//namespace protocols

