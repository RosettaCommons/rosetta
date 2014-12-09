// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

//std
#include <iostream>
#include <map>
#include <string>
#include <sstream>

#include <devel/init.hh>
#include <core/types.hh>
#include <core/io/pdb/file_data.hh>
#include <core/io/pdb/pose_io.hh>

#include <core/pose/Pose.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/AtomVDW.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

#include <core/scoring/packstat/compute_sasa.hh>

#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/standard_mains.hh>

#include <utility/string_util.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

//from comput_sasa.cc
#include <basic/options/keys/packstat.OptionKeys.gen.hh>
#include <ObjexxFCL/format.hh>
#include <utility/exit.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <numeric/xyz.functions.hh>

//from
#include <core/graph/Graph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/etable/Etable.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace std;
using namespace core;
using namespace core::io::pdb;
using namespace core::conformation;
using namespace core::chemical;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::scoring;
using namespace core::scoring::packstat;

using namespace ObjexxFCL::format;
using namespace numeric;
using namespace utility;

static thread_local basic::Tracer TR( "apps.pilot.wendao.gen_cav_ball" );

#define MAXCAVN 9
#define MAXRESD 10.0
#define T_FACTOR 2.0
#define THRESHOLD 10.0

int main( int argc, char * argv [] )
{

	try {

    //normal init
    devel::init( argc, argv );

    pose::Pose pose;
    TR << "[get pdb] Get the pose from pdb ..." << endl;
    if ( option[ in::file::s ].user() )
    {
        core::import_pose::pose_from_pdb( pose, option[ in::file::s ]().vector()[ 0 ] );
    }
    else
    {
        utility_exit_with_message("User did not specify the pdb file!");
    }
    TR << "[get pdb] Done!" << endl;

	//////////////////////////
	//calculate the cavities
	//////////////////////////

	//get user options
	bool surface_accessibility = option[ OptionKeys::packstat::surface_accessibility ]();
	core::Real burial_radius   = option[ OptionKeys::packstat::cavity_burial_probe_radius ]();

	//get the atom spheres and center of the residues
	PosePackData pd = pose_to_pack_data(pose);
	//std::cerr << "spheres len: " << pd.spheres.size() << std::endl;
	//std::cerr << "centers len: " << pd.centers.size() << std::endl;

	//set up the options
	SasaOptions opts;
	opts.prune_max_iters = 999;
	opts.prune_max_delta = 0;
	opts.num_cav_ball_layers = 10;
	opts.frac_cav_ball_required_exposed = 0.00;
	opts.area_cav_ball_required_exposed = 0.00;
	opts.surrounding_sasa_smoothing_window = 3;
	opts.min_cav_ball_radius = option[ OptionKeys::packstat::min_cav_ball_radius ]();
	opts.min_cluster_overlap = option[ OptionKeys::packstat::min_cluster_overlap ]();
	opts.cluster_min_volume = option[ OptionKeys::packstat::cluster_min_volume ]();
	for( PackstatReal pr = 3.0; pr >  2.0; pr -= 0.2 ) opts.probe_radii.push_back(pr);
	for( PackstatReal pr = 2.0; pr >= 0.7; pr -= 0.1 ) opts.probe_radii.push_back(pr);
	opts.prune_cavity_burial_probe_radii.push_back( burial_radius );
	//if swicth is on, do more calculation
	if( surface_accessibility ) {
		for( PackstatReal pr = burial_radius-0.1; pr >= 0.1; pr -= 0.1 ) {
			opts.prune_cavity_burial_probe_radii.push_back(pr);
		}
	}

	TR << "compute MSAs" << std::endl;
	SasaResultOP sr = compute_sasa( pd.spheres, opts );

	////////////////////////////////////////////////////////////////////////////////////////////////
	CavBalls cavballs = sr->cavballs;
	TR << "pruning hidden cav balls " << cavballs.size() << std::endl;
	cavballs = prune_hidden_cavity_balls( cavballs, opts );

	TR << "pruning exposed cav balls " << cavballs.size() << std::endl;
	cavballs = prune_cavity_balls( pd.spheres, cavballs, opts );

	TR << "compute cav ball volumes	" << cavballs.size() << std::endl;
	compute_cav_ball_volumes( cavballs, opts );
	//judge the number of cavballs
	if (cavballs.size()<1)
	{
		//search failed, there is no one left
		utility_exit_with_message("No cavity ball found!");
	}

	//prepare to calculate the lj energy
	//scoring function
	ScoreFunctionOP scorefunc = new ScoreFunction();
	(*scorefunc)( pose );
	//atom type
	chemical::AtomTypeSetCAP atom_set(ChemicalManager::get_instance()->atom_type_set(core::chemical::FA_STANDARD));
	core::Size n_atomtypes = atom_set->n_atomtypes();

	//for pair energy
	etable::EtableCAP et = ScoringManager::get_instance()->etable( "FA_STANDARD_DEFAULT");
	core::Real min_dis = 0.01;
	core::Real min_dis2 = min_dis*min_dis;
	core::Real max_dis = 6.0;
	core::Real max_dis2 = max_dis*max_dis;
	core::Real Wradius = 1.0;
	core::Real lj_switch_dis2sigma = 0.6;
	core::Real lj_slope_intercept = 0.0;
	core::Real lk_min_dis2sigma = 0.89;
	utility::vector1< core::Real > lj_radius( n_atomtypes );
	utility::vector1< core::Real > lj_wdepth( n_atomtypes );
	utility::vector1< core::Real > lj_sigma( n_atomtypes );
	utility::vector1< core::Real > lj_r6_coeff( n_atomtypes );
	utility::vector1< core::Real > lj_r12_coeff( n_atomtypes );
	utility::vector1< core::Real > lj_switch_intercept( n_atomtypes );
	utility::vector1< core::Real > lj_switch_slope( n_atomtypes );
	utility::vector1< core::Real > lk_lambda( n_atomtypes );
	utility::vector1< core::Real > lk_dgfree( n_atomtypes );
	utility::vector1< core::Real > lk_volume( n_atomtypes );
	utility::vector1< core::Real > lk_inv_lambda2( n_atomtypes );
	utility::vector1< core::Real > lk_coeff( n_atomtypes );
	utility::vector1< core::Real > lk_min_dis2sigma_value( n_atomtypes );

	//calculate the coeff
	for( core::Size id=1; id<=n_atomtypes; id++)
	{
		lj_radius[id] = et->lj_radius(id);
		lj_wdepth[id] = et->lj_wdepth(id);
		lk_lambda[id] = et->lk_lambda(id);
		lk_dgfree[id] = et->lk_dgfree(id);
		lk_volume[id] = et->lk_volume(id);
	}

	//make_pairenergy
	core::Real lj_switch_sigma2dis = 1.0/lj_switch_dis2sigma;
	// ctsa - value of the lennard-jones potential at the linear
	//   switch point divided by the wdepth (note that this
	//   coefficient is independent of atomtype)
	core::Real lj_switch_value2wdepth = std::pow( lj_switch_sigma2dis, 12 ) -
	 2.0 * std::pow( lj_switch_sigma2dis, 6 );
	// ctsa - slope of the lennard-jones potential at the linear
	//   switch point times sigma divided by wdepth (note that this
	//   coefficient is independent of atomtype)
	core::Real lj_switch_slope_sigma2wdepth = -12.0 * (
	 std::pow( lj_switch_sigma2dis, 13 ) -
	 std::pow( lj_switch_sigma2dis, 7 ) );

	core::Real sigma,sigma6,sigma12,wdepth;
	core::Real inv_lambda;
	utility::vector1< core::Real > lk_coeff_tmp( n_atomtypes );
	core::Real thresh_dis,inv_thresh_dis2,x_thresh;
	core::Real const inv_neg2_tms_pi_sqrt_pi = { -0.089793561062583294 };

	for( core::Size i=1; i<=n_atomtypes; i++)
	{
		inv_lambda = 1.0/lk_lambda[i];
		lk_inv_lambda2[i] = inv_lambda * inv_lambda;
		lk_coeff_tmp[i] = inv_neg2_tms_pi_sqrt_pi *
		lk_dgfree[i] * inv_lambda;
	}

	for ( core::Size i = 1; i <= n_atomtypes; ++i )
	{
		sigma = Wradius * ( lj_radius[i] + lj_radius[i] );
		//jjh temporary fix to prevent division by zero below
		sigma = ( sigma < 1.0e-9 ? 1.0e-9 : sigma );

		// (bk) modify sigma for hbond donors and acceptors
		//   ctsa -- should these scale down by Wradius as well?

		// pb specific sigma correction for pairs between charged oxygen acceptors (15)
		// pb and hydroxyl oxygen donors (13). sigma correction for all polar H and charged oxygen
		// pb acceptor. Combinations of these corrections allow better prediction of both
		// pb hydroxyl O donor/charged O acceptor and charged NH donor/charged O acceptor
		// pb distances.

		sigma6  = std::pow( sigma, 6 );
		sigma12 = sigma6 * sigma6;
		wdepth = lj_wdepth[i];

		lj_sigma[i] = sigma;
		lj_r6_coeff[i] = -2. * wdepth * sigma6;
		lj_r12_coeff[i] = wdepth * sigma12;

		// ctsa - create coefficients for linear projection of lj repulsive used
		//  for low distance values
		bool lj_use_lj_deriv_slope=true;
		if ( lj_use_lj_deriv_slope )
	       	{
			// ctsa - use the slope of the true lj repulsive at the
			//  linear switch point to create a linear projection of
			//  lj for low distances

			//  slope = wdepth/sigma *
			//          (slope@switch_point*sigma/wdepth)
			lj_switch_slope[i] = (wdepth/sigma)*
			lj_switch_slope_sigma2wdepth;

			// intercept = wdepth*(lj@switch_point/wdepth)
			//             - slope*switch_point_distance
			lj_switch_intercept[i] = wdepth*lj_switch_value2wdepth -
			lj_switch_slope[i]*sigma*lj_switch_dis2sigma;
		}
		else
	       	{
			// ctsa - create a linear projection of lj for low distances which
			//  is defined by a constant y intercept and the true lj repulsive
			//  value at the linear switch point
			lj_switch_slope[i] = -(1./sigma)*lj_switch_sigma2dis*
			(lj_slope_intercept-wdepth*lj_switch_value2wdepth);
			lj_switch_intercept[i] = lj_slope_intercept;
		}

		// ctsa - precalculated lk solvation coefficients
		lk_coeff[i] = lk_coeff_tmp[i] * lk_volume[i];

		// ctsa - when dis/sigma drops below lk_min_dis2sigma,
		//   a constant lk solvation value equal to the value at the
		//   switchover point is used. That switchover-point value
		//   is calculated here and stored in lk_min_dis2sigma_value
		thresh_dis = lk_min_dis2sigma*sigma;
		inv_thresh_dis2 = 1./( thresh_dis * thresh_dis );
		core::Real dis_rad = thresh_dis - lj_radius[i];
		x_thresh = ( dis_rad * dis_rad ) * lk_inv_lambda2[i];
		lk_min_dis2sigma_value[i] = std::exp(-x_thresh) * lk_coeff[i] *
		inv_thresh_dis2;

		dis_rad = thresh_dis - lj_radius[i];
		x_thresh = ( dis_rad * dis_rad ) * lk_inv_lambda2[i];
		lk_min_dis2sigma_value[i] = std::exp(-x_thresh) * lk_coeff[i] *
		inv_thresh_dis2;
	}

	//go through all cav ball
	utility::vector1< core::Real > freeenergys(cavballs.size());
	for( core::Size i=1, l=cavballs.size(); i<=l; i++)
	{
		core::Real score = 0.0;
		core::Real factor = T_FACTOR;
		//go through all residues
		for( core::Size j=1, m=pose.n_residue(); j<=m; j++)
		{
			core::conformation::ResidueCOP res( pose.residue(j).get_self_ptr() );
			//go through all atoms
			for( core::Size k=1, n=res->atoms().size(); k<=n; k++)
			{
				//do
				//for each neighbor residue
				double dis2 = cavballs[i].xyz().distance_squared(res->atoms()[k].xyz());
				int atypid = res->atom_type_index(k);

				// locals
				core::Real ljE,d_ljE,x1,x2;
				core::Real dis;
				core::Real inv_dis,inv_dis2,inv_dis6,inv_dis7,inv_dis12,inv_dis13;
				core::Real dis2sigma;

				// include after local variables to allow data statements to initialize
				core::Real atrE = 0.;
				core::Real d_atrE = 0.;
				core::Real repE = 0.;
				core::Real d_repE = 0.;
				core::Real solvE1 = 0.;
				core::Real solvE2 = 0.;
				core::Real dsolvE1 = 0.;
				core::Real dsolvE2 = 0.;

				//  ctsa - epsilon allows final bin value to be calculated
				if ( dis2 > max_dis2 ) continue;

				if ( dis2 < min_dis2 ) dis2 = min_dis2;

				dis = std::sqrt(dis2);
				inv_dis = 1.0/dis;
				inv_dis2 = inv_dis * inv_dis;

				//lj_sigma
				dis2sigma = dis / lj_sigma[atypid];

				if ( dis2sigma < lj_switch_dis2sigma )
				{
					//  ctsa - use linear ramp instead of lj when the dis/sigma
					//    ratio drops below theshold
					d_ljE = lj_switch_slope[atypid];
					ljE = dis*d_ljE + lj_switch_intercept[atypid];
				}
				else
				{
					//  ctsa - calc regular lennard-jones
					inv_dis6  = inv_dis2 * inv_dis2 * inv_dis2;
					inv_dis7  = inv_dis6 * inv_dis;
					inv_dis12 = inv_dis6 * inv_dis6;
					inv_dis13 = inv_dis12 * inv_dis;

					ljE = lj_r12_coeff[atypid] * inv_dis12 +
					lj_r6_coeff[atypid] * inv_dis6;

					d_ljE = -12.*lj_r12_coeff[atypid] * inv_dis13-6. *
					lj_r6_coeff[atypid] * inv_dis7;
				}

				if ( ljE < 0. )
				{
					atrE = ljE;
					d_atrE = d_ljE;
				}
				else
				{
					repE = ljE;
					d_repE = d_ljE;
				}

				// ctsa - calc lk
				if ( dis2sigma < lk_min_dis2sigma )
				{
					// ctsa - solvation is constant when the dis/sigma ratio
					//   falls below minimum threshold
					solvE1 = lk_min_dis2sigma_value[atypid];
					solvE2 = lk_min_dis2sigma_value[atypid];
					dsolvE1 = dsolvE2 = 0.0;
				}
				else
				{
					core::Real dis_rad = dis - lj_radius[atypid];
					x1 = ( dis_rad * dis_rad ) * lk_inv_lambda2[atypid];
					dis_rad = dis - lj_radius[atypid];
					x2 = ( dis_rad * dis_rad ) * lk_inv_lambda2[atypid];

					solvE1 = std::exp(-x1) * lk_coeff[atypid] * inv_dis2;
					solvE2 = std::exp(-x2) * lk_coeff[atypid] * inv_dis2;

					// ctsa - get d(lk_E)/dr
					dsolvE1 = -2.0 * solvE1 *
					(((dis-lj_radius[atypid])*lk_inv_lambda2[atypid])+inv_dis);
					dsolvE2 = -2.0 * solvE2 *
						(((dis-lj_radius[atypid])*lk_inv_lambda2[atypid])+inv_dis);
				}

				//save the energy
				score+=atrE+repE;
			}//atom
		}//res

		freeenergys[i] = score - factor*log(cavballs[i].radius());

		//save the energy
		cavballs[i].evdw = score;
	}

	ostringstream fnball;
	fnball << "out_ball.pdb";
	ofstream fpball(fnball.str().c_str());
	for( core::Size i=1, l=cavballs.size(); i<=l; i++)
	{
		cavballs[i].cluster_=1;
		if (cavballs[i].radius()>1.0)
		fpball << cavballs[i].hetero_atom_line(	i, i, 0.0 )
		//<< F(8, 2, cavballs[i].evdw)
		//<< F(6, 2, cavballs[i].radius())
	       	<< endl;
	}


	/*
	CavBalls goodballs;
	for( core::Size i=1, l=cavballs.size(); i<=l; i++)
	{
		if(freeenergys[i]<THRESHOLD)goodballs.push_back(cavballs[i]);
	}

	//cluster
	vector1< CavityBallCluster > clusters = compute_cav_ball_clusters( goodballs, opts );
	if (clusters.size()<1)
	{
		//no cluster big enough
		utility_exit_with_message("No cluster big enough!");
	}

	//out put the clusters of cavballs
	for( core::Size i = 1; i <= clusters.size() && i<=MAXCAVN ; i++ )
	{
		ostringstream fnball;
		ostringstream fnpara;
		fnball << "out_ball_" << i << ".pdb";
		fnpara << "out_para_" << i << ".param";
		ofstream fpball(fnball.str().c_str());
		ofstream fppara(fnpara.str().c_str());

		TR << "cluster " << i << ": vol=" << clusters[i].volume << endl;

		//volume,surface_area,surface_accessibility;
		fppara << "volume " << clusters[i].volume << endl;
		fppara << "surface_area " << clusters[i].surface_area << endl;
		fppara << "surface_accessibiliy " << clusters[i].surface_accessibility << endl;

		//statistics of residue type
		core::Size n_res=0;
		//neighbor number of the active residues
		utility::vector1< core::Size > active_res_nn;
		utility::vector1< core::Size > active_res_id;

		//for each residue
		for( core::Size j=1; j<=pose.n_residue(); j++)
		{
			core::conformation::ResidueCOP res( &pose.residue(j) );
			bool flag=true;

			//for each atom
			for( core::Size k=1; k<=res->atoms().size() && flag; k++)
			{
				Vector const & xyzatom = res->atoms()[k].xyz();
				//for each cav ball
				for( core::Size n=1; n<=clusters[i].cavballs.size() && flag; n++ )
				{
					Vector const d = xyzatom - clusters[i].cavballs[n].xyz();

					//cutoff 5.0A
					if ( d.length_squared() < 25.0) flag=false;
				}
			}

			//save the neighbor number of this residue if it is the active res
			if (!flag)
			{
				n_res++;
				active_res_id.push_back(j);
				active_res_nn.push_back( pose.energies().tenA_neighbor_graph()\
							.get_node( j )->num_neighbors_counting_self());
			}
		}
		fppara << "n_active_res " << n_res << endl;

		for( core::Size j=1; j<clusters[i].cavballs.size(); j++)
		{
			//output the ball
			if( clusters[i].cavballs[j].radius() > 0.7 )
			{
				fpball << clusters[i].cavballs[j].hetero_atom_line( \
				pose.total_residue()+i, i, clusters[i].cavballs[j].radius() )
					//	energy
					//<< F(6, 2, score)
					<< F(6, 2, clusters[i].cavballs[j].radius())
				       	<< endl;
			}
		}

		//fppara << "lowest_freeenergy " << min_freeenergy << endl;
		fpball.close();
		fppara.close();
	}//cluster
*/
	return EXIT_SUCCESS;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}


