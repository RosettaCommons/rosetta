#include <core/pose/Pose.hh>
#include <core/scoring/packstat/CavityBall.hh>
#include <core/scoring/packstat/compute_sasa.hh>


utility::vector1< core::scoring::packstat::CavityBallCluster >
sheffler_voids( core::pose::Pose & pose ) {
	using namespace core::scoring::packstat;
 	using namespace std;
	using namespace core;
	using namespace numeric;
	using namespace utility;

	PosePackData pd = pose_to_pack_data(pose);

	core::Real burial_radius = 1.6;
	SasaOptions opts;
	opts.prune_max_iters = 999;
	opts.prune_max_delta = 0;
	opts.num_cav_ball_layers = 10;
	opts.frac_cav_ball_required_exposed = 0.00;
	opts.area_cav_ball_required_exposed = 0.00;
	opts.surrounding_sasa_smoothing_window = 3;
	opts.min_cav_ball_radius = 0.8;
	opts.min_cluster_overlap = 0.1;
	opts.cluster_min_volume = 20.0;
	for( PackstatReal pr = 3.0; pr >  2.0; pr -= 0.2 ) opts.probe_radii.push_back(pr);
	for( PackstatReal pr = 2.0; pr >= opts.min_cav_ball_radius; pr -= 0.1 ) opts.probe_radii.push_back(pr);
	opts.prune_cavity_burial_probe_radii.push_back( burial_radius );

	//TRps << "compute MSAs" << std::endl;
	SasaResultOP sr = compute_sasa( pd.spheres, opts );

	////////////////////////////////////////////////////////////////////////////////////////////////
	CavBalls cavballs = sr->cavballs;
	//TRps << "pruning hidden cav balls " << cavballs.size() << std::endl;
	cavballs = prune_hidden_cavity_balls( cavballs, opts );

	//TRps << "pruning exposed cav balls " << cavballs.size() << std::endl;
	cavballs = prune_cavity_balls( pd.spheres, cavballs, opts );

	//TRps << "compute cav ball volumes	" << cavballs.size() << std::endl;
	compute_cav_ball_volumes( cavballs, opts );

	return compute_cav_ball_clusters( cavballs, opts );

}

#include <devel/init.hh>
#include <iostream>

#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>


int main (int argc, char *argv[])
{

	try {

	devel::init(argc,argv);

	core::pose::Pose test;
	core::import_pose::pose_from_file(test,"test.pdb", core::import_pose::PDB_file);

	utility::vector1< core::scoring::packstat::CavityBallCluster > cavs = sheffler_voids(test);

	for( core::Size i = 1; i <= cavs.size(); i++ ) {
		for( core::Size j = 1; j <= cavs[i].cavballs.size(); j++ ) {
			numeric::xyzVector<core::Real> xyz = cavs[i].cavballs[j].xyz();
			std::cout << xyz.x() << " " << xyz.y() << " " << xyz.z() << " " << cavs[i].cavballs[j].radius() << std::endl;
		}
	}


	return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
