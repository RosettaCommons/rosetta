// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file test_bbmc.cc
/// @brief run Monte Carlo for sampling protein conformation
/// @author Yuan Liu
/// @detailed
/// Modified from Colin's backrub.cc, put all backbone algorithm(backrub, bbg, conrot)
/// and all sidechain algorithm(sc, scmc), and MonteCarlo/ReplicaExchange together

// Protocols Headers
#include <protocols/branch_angle/BranchAngleOptimizer.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/standard_mains.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/backrub/BackrubMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/ReplicaExchangeMC.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMCMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/canonical_sampling/TrajectoryRecorder.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/moves/PyMolMover.hh>
#include <protocols/canonical_sampling/mc_convergence_checks/Pool_ConvergenceCheck.hh>

// Core Headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <devel/init.hh>
#include <core/pose/PDBInfo.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/scoring/mm/MMBondAngleResidueTypeParamSet.hh>
#include <core/types.hh>

#include <basic/Tracer.hh>
#include <basic/basic.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <core/scoring/rms_util.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/RG_Energy_Fast.hh>

//Backbone Gaussian Mover
#include <protocols/simple_moves/BBGaussianMover.hh>
#include <protocols/simple_moves/BBConRotMover.hh>

#include <core/optimization/Minimizer.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

// Utility Headers
#include <utility/io/ozstream.hh>
#include <utility/vector1.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/file/FileName.hh>
#include <utility/exit.hh>

// Numeric Headers
#include <numeric/random/random.hh>
#include <numeric/MultiDimensionalHistogram.hh>

// Platform Headers
#include <platform/types.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/backrub.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/mc.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <protocols/moves/MoverStatistics.hh>
#include <utility/io/mpistream.hh>
#include <fstream>
#include <string>
#include <ObjexxFCL/Fmath.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <basic/prof.hh>

static numeric::random::RandomGenerator RG(62331900);
basic::Tracer TR("pilot.wendao.bbmc");

//params for all
OPT_1GRP_KEY(Integer, mc, ntrials) //how many steps
OPT_1GRP_KEY(Real, mc, sm_prob) //prob of smallmover
OPT_1GRP_KEY(Real, mc, sm_angle_max)
OPT_1GRP_KEY(Real, mc, backrub_prob) //prob of backrubmover
OPT_1GRP_KEY(Real, mc, conrot_prob) //prob of conrotmover
OPT_1GRP_KEY(Real, mc, kt) //temperatrue
OPT_1GRP_KEY(Real, mc, mm_bend_weight) //mm bend energy
OPT_1GRP_KEY(Boolean, mc, detailed_balance) //detailed balance correction
OPT_1GRP_KEY(Boolean, mc, initial_pack) //for packer
OPT_1GRP_KEY(File, mc, movemap) //movemap
OPT_1GRP_KEY(Integer, mc, rmsd_region_start) //beginning of the rmsd region
OPT_1GRP_KEY(Integer, mc, rmsd_region_stop) //end of the rmsd region
OPT_1GRP_KEY(File, mc, minimize_movemap) //movemap
OPT_1GRP_KEY(RealVector, mc, trajectory_tlist)
OPT_1GRP_KEY(Integer, mc, score_stride) //score only silent file
OPT_1GRP_KEY(Integer, mc, trajectory_stride)
OPT_1GRP_KEY(Integer, mc, output_stride)
OPT_1GRP_KEY(Integer, mc, cluster_ndx)
OPT_1GRP_KEY(Boolean, mc, centroid)
OPT_1GRP_KEY(Boolean, mc, noscore)
OPT_1GRP_KEY(Boolean, mc, xyzcst)

//replica
OPT_1GRP_KEY(Boolean, mc, replica)
OPT_1GRP_KEY(String, mc, re_pdb_prefix)
OPT_1GRP_KEY(Boolean, mc, re_pdb_suffix)
OPT_1GRP_KEY(Integer, mc, re_ninterval)
OPT_1GRP_KEY(RealVector, mc, re_tlist)

//params for sidechain
OPT_1GRP_KEY(Real, mc, sc_prob)
OPT_1GRP_KEY(Real, mc, sc_prob_uniform)
OPT_1GRP_KEY(Real, mc, sc_prob_withinrot)
OPT_1GRP_KEY(Real, mc, sc_prob_random_pert_current)
OPT_1GRP_KEY(Boolean, mc, fast_sc)
OPT_1GRP_KEY(Real, mc, fast_sc_prob)//two ways for calling fastsc
OPT_1GRP_KEY(Boolean, mc, sc_strategy2)
OPT_1GRP_KEY(Boolean, mc, fast_sc_strategy2)
OPT_1GRP_KEY(Integer, mc, sc_ntrials)
OPT_1GRP_KEY(IntegerVector, mc, sc_statistic)
OPT_1GRP_KEY(Integer, mc, bb_dih_statistic)
OPT_1GRP_KEY(String,mc,restart_from_silent)

OPT_1GRP_KEY(Real, mc, near_native_threshold )
OPT_1GRP_KEY(Boolean, mc, follow_classic_naming_convention )
OPT_1GRP_KEY(String,mc, movable_segment)

using namespace core;

void *my_main( void* );
void get_resmap( pose::Pose const &pose, pose::Pose const &ref_pose, std::map< Size, Size > &resmap );

///////////////////////////////////////
// test new mover
///////////////////////////////////////

using namespace protocols::simple_moves;

// class BBG8T3A_Jump_Mover : public BBGaussianMover
// {
// public:

// public:
// 	BBG8T3A_Jump_Mover():BBGaussianMover(3,8,4),
// 	    dphi(utility::vector1<Real>(n_dof_angle_))
// 	{
// 	    protocols::moves::Mover::type("BBG8T3A_Jump_Mover");
// 	    //build end atom list
// 	    end_atom_list_.push_back("CA");
// 	    end_atom_list_.push_back("C");
// 	    end_atom_list_.push_back("O");

// 	    //init the A/C factor
// 	    factorA_ = 0.6;
// 	    factorB_ = 10;
// 	}

// 	~BBG8T3A_Jump_Mover(){}

// 	void apply(Pose &pose)
// 	{
// 	    //if(available_seg_list_.size()==0)setup_list(pose);
// 	    //if(available_seg_list_.size()==0)return;

// 	    //select four hinge
// 	    //for lysozyme, 104-117
// 	    Real rA, rB, rC, rD;
// 	    if (RG.uniform()<0.5) {
// 	    	Size Ns(RG.uniform()*2+1);
// 		    Real randn = RG.uniform();
// 		    rA = Size(randn*(15.0-3.0*Ns))+104;
// 		    rB = rA+Ns;
// 		    rC = rB+Ns;
// 		    rD = rC+Ns;
// 	    }
// 	    else {
// 	    	rA = 105;
// 	    	rB = 106;
// 	    	rC = 116;
// 	    	rD = 117;
// 	    }

// 	    std::cout << rA << "," << rB << "," << rC << "," << rD << std::endl;

// 	    get_VdRdPhi(pose, rA, rB, rC, rD);
// 	    get_G();
// 	    get_A();
// 	    Real W_old = get_L_move(pose, rA, rB, rC, rD);

// 	    get_VdRdPhi(pose, rA, rB, rC, rD);
// 	    get_G();
// 	    get_A();
// 	    Real W_new = get_L_prime();

// 	    last_proposal_density_ratio_ = W_new / W_old;
// 	}

// 	virtual std::string get_name() const
// 	{
// 		return "BBG8T3A_Jump_Mover";
// 	}

// protected:
// 	void get_VdRdPhi(Pose const &pose, Size rA, Size rB, Size rC, Size rD)
// 	{
// 	    conformation::Residue const & rsd4( pose.residue( rD ) );
// 	    conformation::Residue const & rsd3( pose.residue( rC ) );
// 	    conformation::Residue const & rsd2( pose.residue( rB ) );
// 	    conformation::Residue const & rsd1( pose.residue( rA ) );

// 	    for (Size i=1;i<=end_atom_list_.size();i++)
// 	    {
// 	        //for each end atom
// 	        xyzVector end_xyz = rsd4.atom(end_atom_list_[i]).xyz();
// 	        //TR << "Phi: 8" << endl;
// 	        matrix_dRdPhi[i][8] = get_dRdPhi(rsd4.atom("CA").xyz(),
// 	                                         rsd4.atom("C").xyz(),
// 	                                         end_xyz);
// 	        //TR << "Phi: 7" << endl;
// 	        matrix_dRdPhi[i][7] = get_dRdPhi(rsd4.atom("N").xyz(),
// 	                                         rsd4.atom("CA").xyz(),
// 	                                         end_xyz);
// 	        //TR << "Phi: 6" << endl;
// 	        matrix_dRdPhi[i][6] = get_dRdPhi(rsd3.atom("CA").xyz(),
// 	                                         rsd3.atom("C").xyz(),
// 	                                         end_xyz);
// 	        //TR << "Phi: 5" << endl;
// 	        matrix_dRdPhi[i][5] = get_dRdPhi(rsd3.atom("N").xyz(),
// 	                                         rsd3.atom("CA").xyz(),
// 	                                         end_xyz);
// 	        //TR << "Phi: 4" << endl;
// 	        matrix_dRdPhi[i][4] = get_dRdPhi(rsd2.atom("CA").xyz(),
// 	                                         rsd2.atom("C").xyz(),
// 	                                         end_xyz);
// 	        //TR << "Phi: 3" << endl;
// 	        matrix_dRdPhi[i][3] = get_dRdPhi(rsd2.atom("N").xyz(),
// 	                                         rsd2.atom("CA").xyz(),
// 	                                         end_xyz);
// 	        //TR << "Phi: 2" << endl;
// 	        matrix_dRdPhi[i][2] = get_dRdPhi(rsd1.atom("CA").xyz(),
// 	                                         rsd1.atom("C").xyz(),
// 	                                         end_xyz);
// 	        //TR << "Phi: 1" << endl;
// 	        matrix_dRdPhi[i][1] = get_dRdPhi(rsd1.atom("N").xyz(),
// 	                                         rsd1.atom("CA").xyz(),
// 	                                         end_xyz);
// 	    }
// 	}

// 	void get_G()
// 	{
// 		for (Size i=1; i<=n_dof_angle_; i++)
// 	    {
// 	        for (Size j=i; j<=n_dof_angle_; j++)
// 	        {
// 	            matrix_G[i][j] = 0.0;
// 	            for (Size n=1; n<=end_atom_list_.size();n++)
// 	            {
// 	                matrix_G[i][j] += matrix_dRdPhi[n][i].dot(matrix_dRdPhi[n][j]);
// 	            }
// 	            //if (matrix_G[i][j]>-ZERO && matrix_G[i][j]<ZERO) matrix_G[i][j] = 0.0;
// 	            if (i<j) matrix_G[j][i]=matrix_G[i][j];

// 	        }
// 	    }
// 	}

// 	void get_A()
// 	{
// 	    for (Size i=1; i<=n_dof_angle_; i++)
// 	    {
// 	        for (Size j=i; j<=n_dof_angle_; j++)
// 	        {
// 	            matrix_A[i][j] = factorB_ * matrix_G[i][j];
// 	            if (i==j) matrix_A[i][j] += 1.0;
// 	            matrix_A[i][j] *= factorA_ / 2.0;
// 	            if (i<j) matrix_A[j][i] = matrix_A[i][j];
// 	        }
// 	    }
// 	}

// 	void get_VdRdPhi(Pose const &){}
// 	//Real get_L_move(Pose &){}

// 	Real get_L_move(Pose &pose, Size rA, Size rB, Size rC, Size rD)
// 	{
// 	    //gerate a Gaussian dx vector
// 	    utility::vector1<Real> delta(n_dof_angle_);
// 	    for (Size i=1; i<=n_dof_angle_; i++) delta[i]=RG.gaussian();
// 	    //calculate d^2 = delta^2
// 	    Real d2=0.0;
// 	    for (Size i=1; i<=n_dof_angle_; i++) d2+=delta[i]*delta[i];
// 	    //cholesky, get L^t, L^-1
// 	    Real detL = cholesky_fw(matrix_A, n_dof_angle_, delta, dphi);

// 	    //W_old *= exp(-d^2)
// 	    Real W_old = detL*exp(-d2/2.0);
// 	    //set the new phi,psi (above all called phi, actually 4 phi, 4 psi)
// 	    pose.set_psi(rD, basic::periodic_range( pose.psi(rD)+dphi[8], 360.0 ) );
// 	    pose.set_phi(rD, basic::periodic_range( pose.phi(rD)+dphi[7], 360.0 ) );
// 	    pose.set_psi(rC, basic::periodic_range( pose.psi(rC)+dphi[6], 360.0 ) ) ;
// 	    pose.set_phi(rC, basic::periodic_range( pose.phi(rC)+dphi[5], 360.0 ) );
// 	    pose.set_psi(rB, basic::periodic_range( pose.psi(rB)+dphi[4], 360.0 ) );
// 	    pose.set_phi(rB, basic::periodic_range( pose.phi(rB)+dphi[3], 360.0 ) );
// 	    pose.set_psi(rA, basic::periodic_range( pose.psi(rA)+dphi[2], 360.0 ) );
// 	    pose.set_phi(rA, basic::periodic_range( pose.phi(rA)+dphi[1], 360.0 ) );

// 	    return W_old;
// 	}

// 	Real get_L_prime()
// 	{
// 		utility::vector1<Real> delta(n_dof_angle_);
// 	    //get L
// 	    Real detL = cholesky_bw(matrix_A, n_dof_angle_, dphi, delta);
// 	    //delta = L^t * dphi
// 	    //calculate d^2 = delta^2
// 	    Real d2=0.0;
// 	    for (Size i=1; i<=n_dof_angle_; i++)d2+=delta[i]*delta[i];
// 	    Real W_new = detL*exp(-d2/2.0);
// 	    return W_new;
// 	}

// private:
// 	utility::vector1< std::string > end_atom_list_;
// 	utility::vector1< Real > dphi;
// 	Real factorA_;
// 	Real factorB_;
// 	Real last_delta_square_;
// };

//////////////////////////////
//////////////////////////////


core::Real
periodic_range( core::Real a, core::Real x ){
  using namespace ObjexxFCL;
  core::Real const halfx = 0.5f * x;
  return ( ( a >= halfx || a < -halfx ) ? mod( mod( a, x ) + ( x + halfx ), x ) - halfx : a );
}

std::string get_ABGEO_string( core::pose::Pose & p, core::Size start, core::Size stop ) {
	std::string ABGEO_assignment = "";
	for( core::Size ii = start; ii <= stop; ii++ ){
		core::Real phi = p.phi(ii);
		core::Real psi = p.psi(ii);
		core::Real omega = p.omega(ii);
		periodic_range( phi  , 360.0 );  //does this get applied to phi??
		periodic_range( psi  , 360.0 );
		periodic_range( omega, 360.0 );
		std::string position_assignment="";
		if ( std::abs( omega ) < 90 ) {
			position_assignment= "O";
		} else if ( phi >= 0.0 ) {
			if ( -100 < psi && psi <= 100 ) {
				position_assignment= "G"; // alpha-L
			} else {
				position_assignment= "E"; // E
			}
		} else {
			if ( -125 < psi && psi <= 50 ) {
				position_assignment= "A"; // helical
			} else {
				position_assignment= "B"; // extended
			}
		}
		ABGEO_assignment = ABGEO_assignment + position_assignment;
	}

  return ABGEO_assignment;
}

std::string get_tag( core::Size itrial, core::Real kT, int rank ) {
	using namespace ObjexxFCL;
	return   "S_" + lead_zero_string_of( itrial, 8 ) + "_"
    + lead_zero_string_of( kT, 5 ).substr(0,8) + "_"
    + lead_zero_string_of( rank, 5);
}

std::string get_filename(std::string const suffix, core::Real t) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	int index=0;
	core::Size end=option[mc::re_tlist]().size();

	if (end>0) {
		for(core::Size i=1; i<=end; i++) {
			if (std::fabs(option[mc::re_tlist]()[i]-t)<1.0e-4) index=static_cast<int>(i-1);
		}
		runtime_assert(index>=0 && index < static_cast<int>(end));
	}

	std::ostringstream inputfn;
	if (option[ mc::re_pdb_prefix ].user()) {
		inputfn << option[ mc::re_pdb_prefix ]();
		inputfn << "_" << index;
	}
	else {
		//using -s
		inputfn << utility::file::FileName(option[in::file::s]().vector()[0]).base() << "_" << 0;
	}

	return inputfn.str() + suffix;
}

int main( int argc, char * argv [] )
{

	try {

	//all
	NEW_OPT(mc::ntrials, "number of Monte Carlo trials to run", 1000);
	NEW_OPT(mc::kt, "value of kT for Monte Carlo", 0.56);
	NEW_OPT(mc::mm_bend_weight, "weight of mm_bend bond angle energy term", 0);
	NEW_OPT(mc::detailed_balance, "preserve detailed balance", true);
	NEW_OPT(mc::initial_pack, "force a repack at the beginning regardless of whether mutations are set in the resfile", false);
	NEW_OPT(mc::movemap, "specify degrees of freedom for simulation", "");
	NEW_OPT(mc::minimize_movemap, "specify degrees of freedom for minimization", "");
	NEW_OPT(mc::trajectory_tlist, "record trajectory for specified T", utility::vector1<core::Real>());
	NEW_OPT(mc::trajectory_stride, "write out a trajectory frame every N steps", 0);
	NEW_OPT(mc::score_stride, "write out a score only silent file every N steps", 0);
	NEW_OPT(mc::output_stride, "write out a statistical info every N steps", 100);
	NEW_OPT(mc::cluster_ndx, "index of the center of clusters in silent file", 0);
	NEW_OPT(mc::centroid, "using centroid mode", false);
	NEW_OPT(mc::noscore, "in absence of score", false);
	NEW_OPT(mc::xyzcst, "coordinate constraints", false);
	NEW_OPT(mc::rmsd_region_start, "beginning of rmsd region", 0);
	NEW_OPT(mc::rmsd_region_stop, "end of rmsd region", 0);
	//bb
	NEW_OPT(mc::sm_prob, "probability of making a small move", 0);
	NEW_OPT(mc::sm_angle_max, "small move maximum band of angluar perturbation", 6);
	NEW_OPT(mc::backrub_prob, "probability of making a backrub move", 0);
	NEW_OPT(mc::conrot_prob, "probability of making a conrot move", 0);
	//sc
	NEW_OPT(mc::sc_prob, "probability of making a side chain move", 0);
	NEW_OPT(mc::sc_prob_uniform, "probability of uniformly sampling chi angles", 0.1);
	NEW_OPT(mc::sc_prob_withinrot, "probability of sampling within the current rotamer", 0.0);
	NEW_OPT(mc::sc_prob_random_pert_current, "probability of sampling within the current rotamer", 0.0);
	NEW_OPT(mc::fast_sc, "using fast sidechainmover", false);
	NEW_OPT(mc::fast_sc_prob, "probability of making a fast side chain move", 0);
	NEW_OPT(mc::fast_sc_strategy2, "using fast sidechainmover, strategy II", false);
	NEW_OPT(mc::sc_strategy2, "using sidechainmover, strategy II", false);
	NEW_OPT(mc::sc_ntrials, "fast sidechainmover(scmc)'s internal sc move", 100 );
	NEW_OPT(mc::sc_statistic, "specify residue id which is gona output chis", 0 );
	NEW_OPT(mc::bb_dih_statistic, "how many dihs to be stat", 0 );
	//replica exchange
	NEW_OPT(mc::replica, "using replica exchange -- mpi only", false);
	NEW_OPT(mc::re_pdb_prefix, "load seperate pdb", "default.pdb");
	NEW_OPT(mc::re_pdb_suffix, "use suffix or not", false);
	NEW_OPT(mc::re_ninterval, "exchange interval of RE", 100);
	NEW_OPT(mc::re_tlist, "temperature list of RE", utility::vector1<core::Real>(1,0));
	NEW_OPT(mc::restart_from_silent,"restart using a specified silent-file","in.out");
	NEW_OPT(mc::near_native_threshold,"threshold with which to regard a structure as near-native",2.5);
	NEW_OPT(mc::follow_classic_naming_convention,"use old (yuan's) naming convention or new (which allows you to track trajectories)",false);
	NEW_OPT(mc::movable_segment,"","movables");
	devel::init(argc, argv);
	protocols::viewer::viewer_main( my_main );

	return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

void
get_resmap( pose::Pose const &pose,
				pose::Pose const &ref_pose,
				std::map< Size, Size > &resmap
		)
{
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
				Size ii_pdb( pose.pdb_info()->number( ii ) );
				//pose_resmap[ii_pdb] = ii;

				for ( Size jj = 1; jj <= ref_pose.total_residue(); ++jj ) {
						Size jj_pdb( ref_pose.pdb_info()->number( jj ) );

						if( ii_pdb == jj_pdb ){
								id::AtomID id1( pose.residue(ii).atom_index( "CA" ), ii );
								id::AtomID id2( ref_pose.residue(jj).atom_index( "CA" ), jj );

								resmap[ii] = jj;
								TR << "Map: " << ii << " " << ii_pdb << " mapped to ";
								TR << jj << " " << jj_pdb << std::endl;
								break;
						}
				}
		}
}

void *
my_main( void* )
{
	using namespace core;
	using namespace core::id;
	using namespace core::io::silent;
	using namespace core::io::pdb;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	using namespace core::pose;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace protocols::moves;

	// create a TaskFactory with the resfile
	using namespace core::pack::task;
	TaskFactoryOP main_task_factory = new TaskFactory;
	main_task_factory->push_back( new operation::InitializeFromCommandline );
	if ( option[ packing::resfile ].user() ) {
		main_task_factory->push_back( new operation::ReadResfile );
	}
	else {
		operation::RestrictToRepackingOP rtrop = new operation::RestrictToRepacking;
		main_task_factory->push_back( rtrop );
	}
	// C-beta atoms should not be altered during packing because branching atoms are optimized
	main_task_factory->push_back( new operation::PreserveCBeta );

	//setup score function
	core::scoring::ScoreFunctionOP score_fxn;
	if ( option[mc::noscore] ) {
		score_fxn = new core::scoring::ScoreFunction();
	}
	else if ( option[mc::centroid] ) {
		if ( option[score::weights].user() ) {
			score_fxn = core::scoring::get_score_function();
		}
		else {
			score_fxn = core::scoring::ScoreFunctionFactory::create_score_function("cen_std");
		}
	}
	else {
		score_fxn = core::scoring::get_score_function();
	}

	//bend score(for conrot and backrub)
	if ( !(option[mc::centroid]) && option[ mc::mm_bend_weight ]>0) {
		score_fxn->set_weight(core::scoring::mm_bend, option[ mc::mm_bend_weight ]);
	}

	core::scoring::methods::EnergyMethodOptions energymethodoptions(score_fxn->energy_method_options());
	energymethodoptions.hbond_options().decompose_bb_hb_into_pair_energies(true);
	energymethodoptions.bond_angle_central_atoms_to_score(option[ backrub::pivot_atoms ]);
	score_fxn->set_energy_method_options(energymethodoptions);
	if ( option[ in::file::centroid_input ].user() ) {
		core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn(*score_fxn);
	}
	else {
		core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn);
	}


	//set up the BackrubMover
	protocols::backrub::BackrubMover backrubmover;
	//read known and unknown optimization parameters from the database
	backrubmover.branchopt().read_database();
	// tell the branch angle optimizer about the score function MMBondAngleResidueTypeParamSet, if any
	if (energymethodoptions.bond_angle_residue_type_param_set()) {
		backrubmover.branchopt().bond_angle_residue_type_param_set(energymethodoptions.bond_angle_residue_type_param_set());
	}
	backrubmover.set_preserve_detailed_balance(option[ mc::detailed_balance ]);

	//setup the SmallMover
	protocols::simple_moves::SmallMover smallmover;
	smallmover.nmoves(1);
	smallmover.set_preserve_detailed_balance(option[ mc::detailed_balance ]);
	if ( option[ mc::sm_angle_max ].user() ) smallmover.angle_max(option[ mc::sm_angle_max ]);

	//setup the BBGMover
	protocols::simple_moves::BBG8T3AMover bbgmover;
	//BBG8T3A_Jump_Mover bbgmover;

	//setup the ConRotMover
	protocols::simple_moves::BBConRotMover bbcrmover;

	//setup Movemap
	if ( option[ mc::movemap ].user() ) {
		core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
		movemap->init_from_file(option[ mc::movemap ]);
		smallmover.movemap(movemap);
		bbgmover.movemap(movemap);
		bbcrmover.movemap(movemap);
	}

	// set up the SidechainMover
	protocols::simple_moves::sidechain_moves::SidechainMover sidechainmover;
	sidechainmover.set_task_factory(main_task_factory);
	sidechainmover.set_prob_uniform(option[ mc::sc_prob_uniform ]);
	sidechainmover.set_prob_withinrot(option[ mc::sc_prob_withinrot ]);
	sidechainmover.set_prob_random_pert_current( option[ mc::sc_prob_random_pert_current ] );
	sidechainmover.set_preserve_detailed_balance(option[ mc::detailed_balance ]);

	//setup switch mover
	protocols::simple_moves::SwitchResidueTypeSetMover to_centroid("centroid");
	protocols::simple_moves::SwitchResidueTypeSetMover to_fullatom("fa_standard");

	//setup MPI
	int rank=0;
	//int size=0;
#ifdef USEMPI
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	//MPI_Comm_size( MPI_COMM_WORLD, &size );
#endif

	/////////////////////////
	///code for recorder
	/////////////////////////
	//RG
	core::scoring::methods::RG_Energy_Fast rge;
	//score
	std::ostringstream inputfn;
	if (option[ mc::re_pdb_prefix ].user()) {
		inputfn << option[ mc::re_pdb_prefix ]();
		inputfn << "_" << rank;
	}
	else {
		//using -s
		inputfn << utility::file::FileName(option[in::file::s]().vector()[0]).base() << "_0";
	}

	//setup native pose for rmsd
	ResidueTypeSetCAP rsd_set = ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	PoseOP native_pose;
	if ( option[in::file::native].user() ) {
		native_pose = new Pose();
		core::import_pose::pose_from_pdb( *native_pose, *rsd_set, option[ in::file::native ]() );
	}

	//setup starting pose
	PoseOP pose = new Pose();
	Pose &p(*pose);
	if ( option[ in::file::s ].user() ) {
		core::import_pose::pose_from_pdb( p, *rsd_set, option[ in::file::s ]().vector()[ 0 ] );
	}
	else if ( option[ mc::re_pdb_prefix ].user() ) {
		std::ostringstream infn;
		infn << option[ mc::re_pdb_prefix ]();
		if (option[mc::re_pdb_suffix])
		{
			infn << "_" << rank;
		}
		infn << ".pdb";
		std::cout << "Loading " << infn.str() << std::endl;
		core::import_pose::pose_from_pdb( p, *rsd_set, infn.str());
	}
	else {
		std::cerr << "User did not specify the pdb file!" << std::endl;
		exit( EXIT_FAILURE );
	}

	/////////////////////////////////
	// setup coordinate cst
	kinematics::MoveMap cstmm;
	cstmm.set_bb(false);
	cstmm.set_chi(false);
	cstmm.set_jump(true);
	core::optimization::AtomTreeMinimizer cstmin;
	core::optimization::MinimizerOptions cstoptions( "lbfgs_armijo_nonmonotone", 1e-2, true, false, false );
	cstoptions.max_iter(5);
	if ( option[mc::xyzcst] ) {
		// virtual atom
		core::pose::addVirtualResAsRoot(p);
		core::Size nmonomerres = p.total_residue()-1;

		// add cst to pose
		for ( Size i = 1; i<=nmonomerres; ++i ) {
			if ( !p.residue(i).is_polymer() ) continue;
				core::conformation::Residue const & nat_i_rsd( p.residue(i) );

				//add constraints on CA
				Size CA_i = nat_i_rsd.atom_index("CA");
				p.add_constraint( new core::scoring::constraints::CoordinateConstraint(
					AtomID(CA_i,i), AtomID(1, p.fold_tree().root()), nat_i_rsd.xyz( CA_i ),
					new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) ) );

		}
	}

	//move it here, yuan
	Pose start_pose(p); //save the starting structure so we can compare by rmsd later ek 1-21-2011

	//switch to centroid
	if ( option[ mc::centroid ] ) {
		to_centroid.apply(p);
	}

	//if constraints are specified, add them!
	if( option[ basic::options::OptionKeys::constraints::cst_file ].user() ) {
		core::scoring::constraints::add_constraints_from_cmdline( p, *score_fxn);
	}

	// setup resmap for gdt
	std::map<Size, Size> resmap;
	if ( option[in::file::native].user() ) {
		get_resmap(*native_pose, p, resmap);
	}

	//call score_fxn after pose loaded
	TR << "Score After PDB Load:" << std::endl;
	score_fxn->show(TR, p);
	TR.flush();

	//ss, score only silent file
	core::io::silent::SilentStructOP ss;
	core::io::silent::SilentFileData sfd;
	if (option[mc::score_stride]) {
		//dump score only
		ss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct( "score" );
		ss->fill_struct(p, inputfn.str()+"_"+ObjexxFCL::lead_zero_string_of(0, 6));
		ss->add_energy( "temperature", 0.0 );
		ss->add_energy( "GDTMM", 0.0 );
		ss->add_energy( "GDTHA", 0.0 );
		ss->add_energy( "rmsd", 0.0 );
		ss->add_energy( "srmsd", 0.0 );

		std::string score_fn(inputfn.str() + ".out");
		if (!utility::file::file_exists(score_fn)) {
		  //new file create header
		  std::ofstream scoreos(score_fn.c_str());
  	  ss->print_header( scoreos );
		}
	}

	Size i=1; //hack, for restarting a run
	Size maxstep=0;
	std::string maxtag;

	//trajactory
	core::io::silent::SilentStructOP trajss;
	core::io::silent::SilentFileData trajsfd;
	if (option[mc::trajectory_stride]>0 ) {
		trajss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out( p );
		if (rank==0) { //only one proc do this
			for (core::Size ndx_traj=1, num_traj=option[mc::trajectory_tlist]().size(); ndx_traj<=num_traj; ++ndx_traj) {
				if( option[mc::follow_classic_naming_convention ]() ) {
					trajss->fill_struct(p, inputfn.str()+"_"+ObjexxFCL::lead_zero_string_of(0, 6));
				}else {
					trajss->fill_struct(p, get_tag( 0, option[mc::trajectory_tlist]()[ndx_traj], rank ) );
				}
				trajss->add_energy( "temperature", 0.0 );
				trajss->add_energy( "GDTMM", 0.0 );
				trajss->add_energy( "GDTHA", 0.0 );
				trajss->add_energy( "rmsd", 0.0  );
				trajss->add_energy( "srmsd", 0.0 );
				std::string traj_fn = get_filename("_traj.out",option[mc::trajectory_tlist]()[ndx_traj]);
				if (utility::file::file_exists(traj_fn)) {
					//need to restart from a traj
					std::cout << "Restarting ..." << std::endl;
					SilentFileData resfd(traj_fn, false, true, option[in::file::silent_struct_type]());
					utility::vector1< std::string > tags;
					resfd.read_tags_fast( traj_fn, tags );
					resfd.read_file( traj_fn );
					//tags = resfd.tags();
					std::cout << "number of tags read in: " << tags.size() << " " << resfd.size() << std::endl;
					for( core::Size ii = 1; ii <= tags.size(); ii++ ) {
							std::stringstream stepstr(tags[ii].substr(2,8));
							Size step; stepstr >> step;
							if (step>maxstep) {
									maxstep = step;
									maxtag = tags[ii];
							}
					}
					//find last pose
					if (maxstep>0) {
							SilentStructOP ress = resfd[ maxtag ];
							ress->fill_pose(p);
							core::pose::clearPoseExtraScores(p);
							i = maxstep+1;
					}
				}
				else {
					//new file create init
					std::ofstream trajos(traj_fn.c_str());
  					trajss->print_header( trajos );
				}
			}
		}
	}

	//nresidue
	core::Size n_res = p.n_residue();
	core::Size mres = static_cast<core::Size>(p.n_residue()*0.5);
	core::Size nsegment = static_cast<core::Size>(option[mc::bb_dih_statistic]);
	core::Size l_offset=mres+1, r_offset=mres;
	for (core::Size ndx=1; ndx<=nsegment; ndx++) {
		if (ndx % 2 == 0) {
			if (r_offset<n_res) r_offset++;
		}
		else {
			if (l_offset>1) l_offset--;
		}
	}


	//debug
	//std::cout << "mres=" << mres << " half_seg=" << half_segment << std::endl;
	//rmsd region
	core::Size rmsd_start=1;
	core::Size rmsd_stop=p.n_residue();
	if (option[mc::xyzcst]) rmsd_stop--; //skip virt
	if (option[mc::rmsd_region_start].user()) rmsd_start=option[mc::rmsd_region_start];
	if (option[mc::rmsd_region_stop].user()) rmsd_stop=option[mc::rmsd_region_stop];

	if (option[mc::bb_dih_statistic].user()) {
		if (option[mc::rmsd_region_start].user()) l_offset = rmsd_start;
		if (option[mc::rmsd_region_stop].user()) r_offset = rmsd_stop;
	}

	//init backrub
	backrubmover.clear_segments();
	backrubmover.set_input_pose(pose);
	if( !option[ mc::movable_segment ].user() ) {
		backrubmover.add_mainchain_segments_from_options();
	} else {
		//backrubmover.add_mainchain_segments_from_options();
		TR.Debug << "using special backrub segment specification" << std::endl;
		backrubmover.clear_segments();
		core::Size min_atoms = 3;
		core::Size max_atoms = 66;
		std::string segments = option[ mc::movable_segment ]();
		utility::vector1<core::id::AtomID> atomids;
		std::ifstream in;
		in.open( segments.c_str());
		core::Size segment_i = 0;
		while( !in.eof() ) {
			core::Size mobilestart,rigidstart,rigidend,mobileend;
			segment_i++;
			in >> mobilestart >> rigidstart >> rigidend >> mobileend;
			TR.Debug << "segment-i: " << segment_i << " mobile-start " << mobilestart << " rigid parts: " << rigidstart << " <--> " << rigidend << " mobile-end: " << mobileend << std::endl;
			for( core::Size ii = mobilestart; ii <= mobileend; ii++ ) {
				if( (ii < rigidstart || ii > rigidend) && ii > 0 ) { //0 means that there is no rigid segment
					atomids.push_back(core::id::AtomID(pose->residue(ii).atom_index("CA"),ii));
				}
			}
			backrubmover.add_mainchain_segments(atomids,min_atoms,max_atoms);
			atomids.clear();
		}
	}

	//optimize
	if (!option[mc::centroid]) {
		sidechainmover.idealize_sidechains(p);
		backrubmover.optimize_branch_angles(p);
	}
	//TR << "Score After Branch Angle Optimization/Side Chain Idealization:" << endl;
	//score_fxn->show(TR, p);
	//TR.flush();

	if ( !option[mc::centroid] && false ) { //add another option
		protocols::simple_moves::PackRotamersMover packrotamersmover;
		packrotamersmover.task_factory(main_task_factory);
		packrotamersmover.score_function(score_fxn);
		packrotamersmover.apply(p);

		// if a minimization movemap was specified, go through a series of minimizations
		if ( option[ mc::minimize_movemap ].user() ) {
			// setup the MoveMaps
			core::kinematics::MoveMapOP minimize_movemap = new core::kinematics::MoveMap;
			minimize_movemap->init_from_file(option[ mc::minimize_movemap ]);
			core::kinematics::MoveMapOP minimize_movemap_progressive = new core::kinematics::MoveMap;

			// setup the MinMover
			protocols::simple_moves::MinMover minmover;
			minmover.score_function(score_fxn);
			minmover.min_type("dfpmin");

			// first minimize just the side chains
			for (core::kinematics::MoveMap::MoveMapTorsionID_Map::const_iterator iter = minimize_movemap->movemap_torsion_id_begin();
						   	iter != minimize_movemap->movemap_torsion_id_end(); ++iter) {
				if (iter->first.second == core::id::CHI) minimize_movemap_progressive->set(iter->first, iter->second);
			}
			minmover.movemap(minimize_movemap_progressive);
			minmover.apply(p);
			//pose->dump_pdb(input_jobs[jobnum]->output_tag(structnum) + "_postminchi.pdb");

			// next minimize the side chains and backbone
			for (core::kinematics::MoveMap::MoveMapTorsionID_Map::const_iterator iter = minimize_movemap->movemap_torsion_id_begin();
					iter != minimize_movemap->movemap_torsion_id_end(); ++iter)
			{
				if (iter->first.second == core::id::BB) minimize_movemap_progressive->set(iter->first, iter->second);
			}
			minmover.movemap(minimize_movemap_progressive);
			minmover.apply(p);
			//pose->dump_pdb(input_jobs[jobnum]->output_tag(structnum) + "_postminbb.pdb");

			// finally minimize everything
			for (core::kinematics::MoveMap::MoveMapTorsionID_Map::const_iterator iter = minimize_movemap->movemap_torsion_id_begin();
				iter != minimize_movemap->movemap_torsion_id_end(); ++iter)
			{
				if (iter->first.second == core::id::JUMP) minimize_movemap_progressive->set(iter->first, iter->second);
			}
			minmover.movemap(minimize_movemap_progressive);
			minmover.apply(p);
			//pose->dump_pdb(input_jobs[jobnum]->output_tag(structnum) + "_postminjump.pdb");
		}
	}

	//setup Monte Carlo
	MonteCarloOP mc;
	core::Real kT = option[ mc::kt ];
	if ( option[ mc::replica ] ) {
		//setup RE
		std::cout << "Seting up Replica Exchange ..." << std::endl;
		mc = new ReplicaExchangeMC(p,*score_fxn, option[mc::re_tlist](), option[mc::re_ninterval]);
		std::cout << "Rank: " <<  rank << " Done!" << std::endl;
		//get the right kT
		kT = option[mc::re_tlist]()[rank+1];
	}
	else {
		//setup normal mc
		mc = new MonteCarlo(p, *score_fxn, kT);
	}

	//setup Sicechainmover correction
	if (option[mc::sc_strategy2]) sidechainmover.set_sampling_temperature(kT);

	//setup SidechainMC mover
	protocols::simple_moves::sidechain_moves::SidechainMCMover scmc;
	if ( option[mc::fast_sc] || option[mc::fast_sc_prob]>0 ) {
		pack::task::PackerTaskOP pt = core::pack::task::TaskFactory::create_packer_task( p );
		scmc.set_task( pt );
		pt->restrict_to_repacking();
		scmc.init_task( p );
		scmc.set_ntrials( option[mc::sc_ntrials] );
		TR << "sc_ntrials are " << scmc.ntrials() << std::endl;
		scmc.set_prob_uniform( option[ mc::sc_prob_uniform ] );
		scmc.set_prob_withinrot( option[ mc::sc_prob_withinrot ] );
		scmc.set_prob_random_pert_current( option[ mc::sc_prob_random_pert_current ] );
		scmc.set_preserve_detailed_balance(option[ mc::detailed_balance ]);
		scmc.set_temperature( kT ); //only for intra mc criteria
		core::scoring::ScoreFunctionOP scfxn = score_fxn->clone();

		if (option[mc::fast_sc_strategy2]) {
			scfxn->set_weight(core::scoring::fa_dun, 0); //turn off the dunbrack term
			scmc.set_preserve_detailed_balance( false ); //don't use detailed balance correction
			scmc.set_sampling_temperature( kT ); //using temperature correction
		}
		//else {
			//fa_dun should be in the original score term
			//unless you know what you are doing
			//scmc.set_preserve_detailed_balance( true );
		//}
		scmc.set_scorefunction( *scfxn );
		scmc.setup( scfxn );
	}

	//viewer
	if (!option[mc::replica]) protocols::viewer::add_monte_carlo_viewer(*mc, "Gaussian", 600, 600);
	//pymol viewer
	//protocols::moves::AddPyMolLink(p, false);

	TR.flush();
#ifdef USEMPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	core::Real sm_prob = option[ mc::sm_prob ];
	core::Real backrub_prob = option[ mc::backrub_prob ];
	core::Real conrot_prob = option[ mc::conrot_prob ];
	//only if in fullatom model
	core::Real sc_prob = option[ mc::sc_prob ];
	core::Real fast_sc_prob = option[ mc::fast_sc_prob ];
	if (option[mc::centroid]) sc_prob=0.0;


	// if this is a continuation of a previous rep-exch run, then read in the
	// structure that corresponds to this temperature
	TR << "reading in structure corresponding to this node's temperature: " << kT << std::endl;
	if( option[mc::restart_from_silent].user() ){
		std::string silent_input = option[mc::restart_from_silent]();
		//read in silent-file
		//read tags and pick temperature that matches
		//throw an error if you cannot find the structure corresponding to the right temperature...
		SilentFileData sfd(silent_input,false,true,option[in::file::silent_struct_type]());
		utility::vector1< std::string > tags;
		sfd.read_tags_fast( silent_input, tags );
		sfd.read_file( silent_input );
		tags = sfd.tags();
		bool decoy_found = false;
		std::cout << "number of tags read in: " << tags.size() << " " << sfd.size() << std::endl;
		for( core::Size ii = 1; ii <= tags.size(); ii++ ) {
			TR << "tag of structure: |" << tags[ ii ] << "|" << std::endl;
			SilentStructOP ss = sfd[ tags[ ii ] ];
			core::Real temp = ss->get_energy( "temperature" );
			TR << "temperature of structure with tag: " << tags[ii] << " is: " << temp << " comparing the kT " << kT << std::endl;
			if( fabs(temp - kT) <= 1.0e-4 ) {
				runtime_assert( ss );
				TR << "found the structure!" << std::endl;
				//runtime_assert(p);
				ss->fill_pose(p);
				TR << "total residue: " << p.total_residue() << std::endl;
				TR << "finished filling pose " << std::endl;
				TR << "temperature of structure with tag: " << tags[ii] << " corresponds to this node's temperature: " << kT << " filled pdb has " << p.total_residue() << std::endl;
				decoy_found = true;
				break;
			}
			//
		}
		if( !decoy_found ) {
			utility_exit_with_message("we did not find structure corresponding to temperature:  you cannot use option:  mc::restart_from_silent, exiting");
		} else {
			//TR << "starting from decoy "  << ss->decoy_tag() << " corresponding to a kT of " << kT << std::endl;
			std::ostringstream os;
			os << "input_" << rank << "_" << kT << ".pdb";
			p.dump_pdb(os.str());
		}
	}

	Size ntrials = option[ mc::ntrials ];
	std::cout << "Job is working ..." << std::endl;
	basic::prof_reset();
	//"i", init from 1 or traj
	for (; i <= ntrials; ++i) {
		//init
		std::string move_type("fake");
		Real proposal_density_ratio=1.0;
		//random number
		core::Real prob = RG.uniform();

		//choose one
		if ( prob > sm_prob + backrub_prob + conrot_prob + sc_prob + fast_sc_prob ) { //bbg
			bbgmover.apply(p);
			move_type = bbgmover.type();
			proposal_density_ratio = bbgmover.last_proposal_density_ratio();

			//if xyzcst, refit
			if (option[mc::xyzcst]) {
				cstmin.run( p, cstmm, *score_fxn, cstoptions );
			}

			mc->boltzmann(p, move_type, proposal_density_ratio);
		}
		else if ( prob > backrub_prob+conrot_prob+sc_prob+ fast_sc_prob ) { //small
			smallmover.apply(p);
			move_type = smallmover.type();
			proposal_density_ratio = smallmover.last_proposal_density_ratio();
			mc->boltzmann(p, move_type, proposal_density_ratio);
		}
		else if ( prob > conrot_prob+sc_prob+fast_sc_prob ) { //backrub
			backrubmover.apply(p);
			move_type = backrubmover.type();
			proposal_density_ratio = backrubmover.last_proposal_density_ratio();
			mc->boltzmann(p, move_type, proposal_density_ratio);
		}
		else if ( prob > sc_prob+fast_sc_prob ) { //conrot
			bbcrmover.apply(p);
			move_type = bbcrmover.type();
			proposal_density_ratio = bbcrmover.last_proposal_density_ratio();
			mc->boltzmann(p, move_type, proposal_density_ratio);
		}
		else if ( prob > fast_sc_prob && !option[mc::centroid] ) { //sidechain
			//TR << "probabilities are: (SC) " << sidechainmover.prob_uniform() << " " << sidechainmover.prob_withinrot() << " " << sidechainmover.prob_random_pert_current() << std::endl;
			//	exit(1);
			sidechainmover.apply(p);
			move_type = sidechainmover.type();
			proposal_density_ratio = sidechainmover.last_proposal_density_ratio();
			mc->boltzmann(p, move_type, proposal_density_ratio);
		}
		else if (!option[mc::centroid]) { //fast sidechain
			//TR << "probabilities are (FAST-SC): " << scmc.prob_uniform() << " " << scmc.prob_withinrot() << " " << scmc.prob_random_pert_current() << std::endl;
			//exit(1);
			scmc.apply(p);
			mc->set_last_accepted_pose(p);
		}

		//fast sidechain
		if (option[ mc::fast_sc ] && !option[mc::centroid]) {
			//only if using fast_sc in fullatom model
			scmc.apply(p);
			//save the last accepted pose, but we will lost the lowest energy pose
			mc->set_last_accepted_pose(p);
		}

		//sync temperature
		if (!option[mc::centroid]) {
			//reset sidechainmover's temperature
			if (option[mc::sc_strategy2]) sidechainmover.set_sampling_temperature(mc->temperature());
			if (option[mc::fast_sc_strategy2]) scmc.set_sampling_temperature(mc->temperature());
			scmc.set_temperature(mc->temperature()); //should be
		}

		if ( i%option[mc::output_stride] == 0 ) {
			//output
			TR << "STEP=" << i;
			TR << " T=" << mc->temperature() << " " << get_ABGEO_string(p, rmsd_start, rmsd_stop) ;
			if (mres>1) TR << " RG=" << rge.calculate_rg_score(p);
			core::Real rmsd = 10000;
			if (native_pose) {
				rmsd = core::scoring::CA_rmsd(p,*native_pose, rmsd_start, rmsd_stop);
				TR << " RMSD=" << rmsd;
			}
			TR << " SCORE=" << (*score_fxn)(p);
			TR << std::endl;

			if( rmsd <= option[ mc::near_native_threshold ]() ){
				TR << "near native structure produced! " << rmsd << std::endl;
				SilentFileData sfd_nn("near_natives.out",false,false,option[ out::file::silent_struct_type ]());
				core::io::silent::SilentStructOP ss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();
				std::ostringstream os;
				//os << "nn_" << rank << "_" << i << "_" << kT;
				ss->fill_struct( p, get_tag( i, kT, rank ) );
				sfd_nn.write_silent_struct( *ss, "near_natives.out", false );
			}

			//side chain stat
			if (option[mc::sc_statistic].user()) {
				for (core::Size ndx=1, nr=option[mc::sc_statistic]().size(); ndx<=nr; ++ndx) {
					int nres = option[mc::sc_statistic]()[ndx];
					TR << "STAT_RES_" << nres << "_CHI: ";
					utility::vector1<core::Real> const &chis(p.residue(nres).chi());
					for (core::Size j=1; j<=chis.size(); j++) {
						TR << chis[j] << " ";
					}
					TR << std::endl;
				}
			}

			//middle dihs stat
			if (l_offset<=mres && mres>1) {
				TR << "STAT_DIHS: ";
				for (core::Size j=l_offset; j<=r_offset; j++)
				{
					TR << p.phi(j) << " " << p.psi(j) << " ";
				}
				TR << std::endl;
			}

#ifdef USEMPI
			if (rank==0) std::cout << "STEP: " << i << std::endl;
#endif
		}

		//save sore only silent file
		if ( option[mc::score_stride]>0 ) {
			if ( i%option[mc::score_stride]==0 ) {
				//dump silent file
				int ndump = static_cast<int>(i/option[mc::score_stride]);
				if(option[mc::follow_classic_naming_convention]()) {
					ss->fill_struct(p, inputfn.str()+"_"+ObjexxFCL::lead_zero_string_of(ndump, 6));
				}else{
					ss->fill_struct(p, get_tag(i,mc->temperature(),rank) );
				}
				ss->add_energy( "temperature", mc->temperature() );
				if( option[ in::file::native ].user() ) {
				  ss->add_energy( "GDTMM", core::scoring::CA_gdtmm(*native_pose, p, resmap) );
				  ss->add_energy( "GDTHA", core::scoring::gdtha(*native_pose, p, resmap) );
					ss->add_energy( "rmsd", core::scoring::CA_rmsd( *native_pose, p) ); //ek add rmsd to native in silent-structure header information
				}
				ss->add_energy( "srmsd", core::scoring::CA_rmsd( start_pose, p ) ); //ek add rmsd to starting structure in silent-structure header information
				sfd.write_silent_struct(*ss, get_filename(".out", mc->temperature()), true );
			}
		}

		//save trajctory for specified temperature
		if ( option[mc::trajectory_stride]>0 ) {
			if ( i%option[mc::trajectory_stride]==0 ) {
				for (core::Size ndx_traj=1, num_traj=option[mc::trajectory_tlist]().size(); ndx_traj<=num_traj; ++ndx_traj) {
					if (fabs(mc->temperature() - option[mc::trajectory_tlist]()[ndx_traj])>1.0e-4) continue;
					//dump silent file
					int ndump = static_cast<int>(i/option[mc::trajectory_stride]);
					if( option[mc::follow_classic_naming_convention]() ){
						trajss->fill_struct(p, inputfn.str()+"_"+ObjexxFCL::lead_zero_string_of(ndump, 6));
					} else {
						trajss->fill_struct(p, get_tag(i,mc->temperature(),rank) );
					}
					trajss->add_energy( "temperature", mc->temperature() );
					if( option[ in::file::native ].user() ) {
					trajss->add_energy( "GDTMM", core::scoring::CA_gdtmm(*native_pose, p, resmap) );
					trajss->add_energy( "GDTHA", core::scoring::gdtha(*native_pose, p, resmap) );
						trajss->add_energy( "rmsd", core::scoring::CA_rmsd( *native_pose, p) ); //ek add rmsd to native in silent-structure header information
					}
					trajss->add_energy( "srmsd", core::scoring::CA_rmsd( start_pose, p ) ); //ek add rmsd to starting structure in silent-structure header information
					trajsfd.write_silent_struct(*trajss, get_filename("_traj.out", mc->temperature()), false );
				}
			}
		}
	} // for i = 1:ntrials
	basic::prof_show();
	//finish
	std::string const outfn_last(inputfn.str() + "_last.pdb");
	std::string const outfn_lowest(inputfn.str() + "_lowest.pdb");
	std::cout << "Saving " << outfn_last << " and " << outfn_lowest << std::endl;

	mc->last_accepted_pose().dump_pdb(outfn_last);
	mc->lowest_score_pose().dump_pdb(outfn_lowest);
	mc->show_counters();

	backrubmover.branchopt().write_database();

#ifdef USEMPI
  MPI_Barrier( MPI_COMM_WORLD );
  MPI_Finalize();
#endif

	return 0;
}
