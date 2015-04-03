// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/vip/VIP_Utils.cc
/// @brief Utility functions for VIP mover

#include <core/types.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/packstat/types.hh>
#include <core/pose/Pose.hh>
#include <protocols/simple_filters/PackStatFilter.hh>
#include <basic/options/option.hh>

namespace protocols {
namespace vip {

using core::Real;
using namespace core::scoring::packstat;

std::string base_name(const std::string& str) {
  size_t begin = 0;
  size_t end = str.length();

  for (size_t i=0; i<str.length(); ++i) {
    if (str[i] == '/') begin = i+1;}
		return str.substr(begin,end);}

std::string get_out_tag(std::string fname) {
	std::string base = base_name(fname);
	std::transform( base.begin(), base.end(), base.begin(), tolower );
	if (system( ("mkdir -p out/" + base.substr(1,2)).c_str() ) == -1) {
		return "Error! Shell process failed!";
	}
	std::string OUT_TAG = "out/" + base.substr(1,2) + "/" + base;
	return OUT_TAG;}


core::Real output_packstat( core::pose::Pose & pose ) {

	using core::Size;
	using namespace core::scoring::packstat;
  	using namespace std;
	using namespace basic;
	using namespace options;
	using namespace OptionKeys;
//	using namespace ObjexxFCL::format;
	using namespace numeric;
	using namespace utility;
/*
	core::Size oversample = option[ OptionKeys::packstat::oversample ]();
	bool include_water  = option[ OptionKeys::packstat::include_water ]();
	bool residue_scores = option[ OptionKeys::packstat::residue_scores ]();
	bool packstat_pdb   = option[ OptionKeys::packstat::packstat_pdb ]();
	bool raw_stats      = option[ OptionKeys::packstat::raw_stats ]();

	AtomRadiusMap arm;
	SimplePDB pdb;
	utility::io::izstream in(fname.c_str());
	in >> pdb;

	Spheres spheres;
	spheres = pdb.get_spheres(arm);
	vector1< xyzVector<PackstatReal> > centers( pdb.get_res_centers() );

	PosePackData pd;
	pd.spheres = spheres;
	pd.centers = centers;
	core::Real packing_score = compute_packing_score( pd, oversample );

	utility::vector1<core::Real> res_scores; // needed if output res scores or pdb
	if( raw_stats ) {	// stupid duplicate code....
		assert( pd.spheres.size() > 0 );
		assert( pd.centers.size() > 0 );

		SasaOptions opts;
		opts.surrounding_sasa_smoothing_window = 1+2*oversample;
		opts.num_surrounding_sasa_bins = 7;
		for( core::Size ipr = 1; ipr <= 31; ++ipr ) {
			PackstatReal pr = 3.0 - ((double)(ipr-1))/10.0;
			PackstatReal ostep = 0.1 / (oversample*2.0+1.0);
			for( core::Size i = 1; i <= oversample; ++i )	opts.probe_radii.push_back( pr + i*ostep );
			opts.probe_radii.push_back( pr );
			for( core::Size i = 1; i <= oversample; ++i )	opts.probe_radii.push_back( pr - i*ostep );}
		SasaResultOP result = compute_sasa( pd.spheres, opts );
		for( core::Size i = 1; i <= pd.centers.size(); ++i ) {
			// std::cout << i << " ";
			PackingScoreResDataCOP dat( compute_surrounding_sasa( pd.centers[i], pd.spheres, result, opts ) );
			for( Size i =1; i <= dat->nrad(); ++i ) {
				for( Size j =1; j <= dat->npr(); ++j ) {
				}}
			}}
*/

	using namespace protocols;
	using namespace simple_filters;

	PackStatFilter();
	PackStatFilter psfilter;

	core::Real packing_score = psfilter.compute( pose );


	return packing_score;}


	void
	set_local_movemap(
		core::pose::Pose & pose,
		core::Size position,
		core::kinematics::MoveMapOP mmap ) {

		// This helper function sets up a movemap to only allow freedom for
		// a position and its neighbors

		core::scoring::EnergyGraph const & energy_graph( pose.energies().energy_graph() );

		// Defaulting to immobile
		mmap->set_chi( false );
		mmap->set_bb( false );
		mmap->set_jump( false );

    for ( core::graph::Graph::EdgeListConstIter
        iru = energy_graph.get_node( position )->const_upper_edge_list_begin(),
        irue = energy_graph.get_node( position )->const_upper_edge_list_end();
        iru != irue; ++iru ) {
      int const nbr_position( (*iru)->get_second_node_ind() );
			//std::cout << "Neighbors for " << position << " include " << nbr_position << std::endl;

			mmap->set_chi( nbr_position, true );
			mmap->set_bb( nbr_position, true );
		}
	}


}}


