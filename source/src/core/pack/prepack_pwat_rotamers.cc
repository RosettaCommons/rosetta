// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/pack_rotamers.cc
/// @brief  pack rotamers module
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/pack/prepack_pwat_rotamers.hh>

// Package Headers
#include <core/pack/packer_neighbors.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSetsBase.hh>
#include <core/pack/task/IGEdgeReweightContainer.hh>

#include <core/pack/annealer/AnnealerFactory.hh>
#include <core/pack/annealer/SimAnnealerBase.hh>
#include <core/pack/annealer/FixbbPwatSimAnnealer.hh>
#include <core/pack/interaction_graph/InteractionGraphFactory.hh>
#include <core/pack/interaction_graph/AnnealableGraphBase.hh>

#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <utility/graph/Graph.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/lkball/LK_BallInfo.hh>

// util
#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>

// option key includes
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>


using namespace ObjexxFCL;

namespace core {
namespace pack {

using core::conformation::symmetry::SymmetryInfoCOP;
using core::conformation::symmetry::SymmetricConformation;

static basic::Tracer TR("core.pack.prepack_pwat_rotamers",basic::t_info );

void
get_clash_and_neighbors(
	pose::Pose const & pose,
	task::PackerTaskCOP task,
	utility::vector1< core::Vector > & clash_atoms,
	utility::vector1< core::Vector > & neighbor_atoms
)
{
	Size nres = pose.total_residue();
	// get bump check atom positions
	for ( Size i=1; i <= nres; ++i ) {
		core::conformation::Residue const &res = pose.residue(i);
		if ( res.is_water() ) continue;
		if ( res.is_protein() ) {
			// get all atoms for buriedness check
			if ( res.name3() != "GLY" ) {
				neighbor_atoms.push_back( res.atom("CB").xyz() );
			} else {
				neighbor_atoms.push_back( res.atom("CA").xyz() );
			}
			// get all backbone atoms
			for ( int j=1; j<=(int)res.last_backbone_atom(); ++j ) {
				clash_atoms.push_back( res.atom(j).xyz() );
			}
			// if gly is not a target type, add CB to the clash list
			bool gly_allowed = false;
			core::pack::task::ResidueLevelTask const &task_i = task->residue_task(i);
			for ( std::list< core::chemical::ResidueTypeCOP >::const_iterator type( task_i.allowed_residue_types_begin() );
					type != task_i.allowed_residue_types_end() && !gly_allowed; ++type ) {
				gly_allowed = ((**type).aa() == core::chemical::aa_gly);
			}
			if ( !gly_allowed && res.aa() != core::chemical::aa_gly ) {
				clash_atoms.push_back( res.atom("CB").xyz() );
			}
		} else {
			for ( Size j=1; j<= res.nheavyatoms(); ++j ) {
				clash_atoms.push_back( res.atom(j).xyz() );
			}
		}
	}
}


// return centroid from input points
Vector
centroid(
	utility::vector1< Vector > cluster
)
{
	Vector cent(0);
	for ( Size i = 1; i <= cluster.size(); ++i ) {
		cent += cluster[i];
	}
	return (cent /= cluster.size());
}


// bottom-up clustering based on distance from centroid
void
cluster_sites(
	utility::vector1< Vector > overlap_waters,
	utility::vector1< utility::vector1< Vector > > & overlap_clusters,
	Real clust_rad
)
{
	// assign first overlap_water to first cluster
	utility::vector1< Vector > newcluster;
	newcluster.push_back( overlap_waters[1] );
	overlap_clusters.push_back( newcluster );
	// check rest of overlap water sites
	bool added;
	for ( Size i=2; i<=overlap_waters.size(); ++i ) {
		added = false;
		for ( Size j=1; j<=overlap_clusters.size(); ++j ) {
			Vector clust_cent = centroid(overlap_clusters[j]);
			if ( overlap_waters[i].distance_squared(clust_cent) <= clust_rad*clust_rad ) {
				overlap_clusters[j].push_back( overlap_waters[i] );
				added = true;
				break;
			}
		}
		if ( added ) continue;
		utility::vector1< Vector > newcluster;
		newcluster.push_back( overlap_waters[i] );
		overlap_clusters.push_back( newcluster );
	}
}


void
get_all_possible_water_sites(
	pose::Pose & pose,
	task::PackerTaskCOP task,
	rotamer_set::RotamerSetsOP const & rotsets,
	utility::vector1< Vector > const & clash_atoms,
	utility::vector1< Vector > const & neighbor_atoms,
	utility::vector1< utility::vector1< Vector > > & lkb_wats,
	utility::vector1< utility::vector1< Vector > > & bb_wats,
	utility::vector1< Vector > & non_pack_waters,
	utility::vector1< Vector > & all_bb_pwat,
	bool exclude_exposed
)
{

	// backbone atom_pairs for superposition when cycling through rotamers
	utility::vector1< std::pair< std::string, std::string > > atom_pairs;
	atom_pairs.push_back(std::pair<std::string, std::string>("C","C") );
	atom_pairs.push_back(std::pair<std::string, std::string>("CA","CA") );
	atom_pairs.push_back(std::pair<std::string, std::string>("N","N") );

	// get pwat and lkball positions from packable side chains
	Size nexposed = 0;
	Size n_total_lkb = 0;
	time_t const get_sites_start = time(NULL);
	TR << "collecting all possible solvation sites..." << std::endl;
	for ( Size i=1; i <= task->total_residue(); ++i ) {

		// get sidechain lkball positions for packable / designable sidechains
		if ( task->design_residue( i ) || task->pack_residue( i ) ) {
			// get single residue rotamer_set from rotsets
			using namespace core::pack::rotamer_set;
			RotamerSetOP res_rotset = rotsets->rotamer_set_for_residue( i );

			// iterate over all possible rotamers,
			// apply to pose, and get lkball waters for each
			utility::vector1< Vector > all_waters_for_single_res;
			for ( Size ii=1; ii <= res_rotset->num_rotamers(); ++ii ) {
				core::conformation::Residue res_rot = *res_rotset->rotamer( ii );

				// for pwat solvation sites
				if ( (res_rot.name() == "PWAT") || (res_rot.name() == "PWAT_V") ) {
					pose.replace_residue( i, res_rot, false );
					core::Vector O ( res_rot.xyz("O") );

					// perform clash check
					bool isclash = false;
					for ( Size kk=1; kk <= clash_atoms.size() && !isclash; ++kk ) {
						core::Real dist2 = O.distance_squared(clash_atoms[kk]);
						isclash = (dist2 < 2.5*2.5);
					}
					if ( isclash ) continue;

					// perform exposure check
					bool isexposed = true;
					Size nneighbors = 0;
					for ( Size kk=1; kk <= neighbor_atoms.size() && isexposed; ++kk ) {
						core::Real dist2 = O.distance_squared(neighbor_atoms[kk]);
						if ( dist2 <= 10.0*10.0 ) ++nneighbors;
						if ( nneighbors > 16 ) isexposed = false;
					}
					if ( isexposed && exclude_exposed ) {
						++nexposed;
						continue;
					}
					all_bb_pwat.push_back( O ); // add to background for check with lkball overlap
					all_waters_for_single_res.push_back( O ); // add to single_res for pwat/pwat overlap check

				} else {
					// lkball sites
					pose.replace_residue( i, res_rot, atom_pairs );

					core::scoring::lkball::LKB_ResidueInfoOP lkb_resinfo( new core::scoring::lkball::LKB_ResidueInfo( res_rot ) );
					utility::vector1< utility::fixedsizearray1< Vector, 4 > > res_waters = lkb_resinfo->waters();
					//          utility::vector1< utility::vector1< Vector > > res_waters = lkb_resinfo->waters();
					// extract lkball->waters()
					for ( Size j=1; j <= res_waters.size(); ++j ) {
						if ( res_waters[j].size() > 0 ) {
							for ( Size k=1; k <= res_waters[j].size(); ++k ) {

								// perform clash check before keeping water
								bool isclash = false;
								for ( Size kk=1; kk <= clash_atoms.size() && !isclash; ++kk ) {
									core::Real dist2 = res_waters[j][k].distance_squared(clash_atoms[kk]);
									isclash = (dist2 < 2.5*2.5);
								}
								if ( isclash ) continue;

								// perform exposure check
								bool isexposed = true;
								Size nneighbors = 0;
								for ( Size kk=1; kk <= neighbor_atoms.size() && isexposed; ++kk ) {
									core::Real dist2 = res_waters[j][k].distance_squared(neighbor_atoms[kk]);
									if ( dist2 <= 10.0*10.0 ) ++nneighbors;
									if ( nneighbors > 16 ) isexposed = false;
								}
								if ( isexposed && exclude_exposed ) {
									++nexposed;
									continue;
								}
								// passed clash and exposure check; keep water
								all_waters_for_single_res.push_back( res_waters[j][k] );
							}
						}
					}
				}
			} // for res_rotset

			if ( (pose.residue(i).name() == "PWAT") || (pose.residue(i).name() == "PWAT_V") ) {
				if ( all_waters_for_single_res.size() > 0 ) {
					bb_wats.push_back( all_waters_for_single_res );
				}
			} else {
				// remove duplicate backbone lkball sites
				if ( all_waters_for_single_res.size() > 0 ) {
					utility::vector1< Vector > nonredundant_single_res;
					for ( Size k=1; k<=all_waters_for_single_res.size(); ++k ) {
						if ( !nonredundant_single_res.has_value( all_waters_for_single_res[k] ) ) {
							nonredundant_single_res.push_back( all_waters_for_single_res[k] );
						}
					}
					lkb_wats.push_back( nonredundant_single_res );
					n_total_lkb += nonredundant_single_res.size();
				}
			}
		} else {  // residue not being packed/designed
			// get lkball positions for fixed residues
			core::scoring::lkball::LKB_ResidueInfoOP lkb_resinfo( new core::scoring::lkball::LKB_ResidueInfo( pose.residue(i) ) );
			utility::vector1< utility::fixedsizearray1< Vector, 4 > > res_waters = lkb_resinfo->waters();
			//utility::vector1< utility::vector1< Vector > > res_waters = lkb_resinfo->waters();

			for ( Size j=1; j <= res_waters.size(); ++j ) {
				if ( res_waters[j].size() > 0 ) {
					for ( Size k=1; k <= res_waters[j].size(); ++k ) {
						non_pack_waters.push_back( res_waters[j][k] );
					}
				}
			}
		}
	}  // loop over nres

	TR << "total number of " << n_total_lkb << " lkb sites on " << lkb_wats.size() << " different side chains." << std::endl;
	TR << "number of exposed positions excluded = " << nexposed << std::endl;
	TR << "total number of lkb sites from fixed side chains = " << non_pack_waters.size() << std::endl;
	TR << "total number of bb pwat sites = " << all_bb_pwat.size() << std::endl << std::endl;

	TR << "time spent collecting all water sites = " << time(NULL) - get_sites_start << " seconds" << std::endl << std::endl;
}


//
//
void
find_overlap(
	utility::vector1< utility::vector1< Vector > > lkb_wats,
	utility::vector1< utility::vector1< Vector > > bb_wats,
	utility::vector1< Vector > non_pack_waters,
	utility::vector1< Vector > & overlap_waters,
	bool use_average
)
{
	using namespace basic::options;

	time_t const compare_start = time(NULL);
	TR << "seaching for overlap of possible solvation sites" << std::endl;
	Real overlap_dist = option[ OptionKeys::corrections::water::lkb_overlap_distance ].value();

	utility::vector1< Vector > all_bb_pwat;
	for ( Size i=1; i <= bb_wats.size(); ++i ) {
		for ( Size j=1; j <= bb_wats[i].size(); ++j ) {
			all_bb_pwat.push_back(bb_wats[i][j]);
		}
	}
	TR << "total number of bb_pwat = " << all_bb_pwat.size() << std::endl;

	//
	// find overlap of lkball sites
	//
	Size end_val = ( use_average ) ? lkb_wats.size()-1 : lkb_wats.size();
	for ( Size i=1; i <= end_val; ++i ) {
		Vector closest_lkb;
		for ( Size ii=1; ii <= lkb_wats[i].size(); ++ii ) {

			Real min_dist2 = 1e6;

			// compare current lkb position to fixed sidechain lkb positions
			for ( Size k=1; k <= non_pack_waters.size(); ++k ) {
				Real dist2 = lkb_wats[i][ii].distance_squared(non_pack_waters[k]);
				if ( dist2 < min_dist2 ) {
					min_dist2 = dist2;
					closest_lkb = non_pack_waters[k];
				}
			}

			// compare to all bb_pwat positions
			for ( Size k=1; k <= all_bb_pwat.size(); ++k ) {
				Real dist2 = lkb_wats[i][ii].distance_squared(all_bb_pwat[k]);
				if ( dist2 < min_dist2 ) {
					min_dist2 = dist2;
					closest_lkb = all_bb_pwat[k];
				}
			}

			// compare current lkb position to all other lkb positions
			Size start_val = ( use_average ) ? ( i+1 ) : 1;
			for ( Size j=start_val; j <= lkb_wats.size(); ++j ) {
				if ( j == i ) continue;
				for ( Size jj=1; jj <= lkb_wats[j].size(); ++jj ) {
					Real dist2 = lkb_wats[i][ii].distance_squared(lkb_wats[j][jj]);
					if ( dist2 < min_dist2 ) {
						min_dist2 = dist2;
						closest_lkb = lkb_wats[j][jj];
					}
				}
			}

			// fill overlap_waters vector based on overlap_dist cutoff
			if ( min_dist2 <= overlap_dist*overlap_dist ) {
				if ( use_average ) {
					overlap_waters.push_back( (lkb_wats[i][ii]+closest_lkb)/2.0 );
				} else {
					overlap_waters.push_back( lkb_wats[i][ii] );
				}
			}
		}
	}

	//
	// find overlap of backbone pwat sites
	//
	end_val = ( use_average ) ? bb_wats.size()-1 : bb_wats.size();
	for ( Size i=1; i <= end_val; ++i ) {
		Vector closest_bbwat;
		for ( Size ii=1; ii <= bb_wats[i].size(); ++ii ) {

			Real min_dist2 = 1e6;

			// compare current bb_pwat position to fixed sidechain lkb positions and pwat positions
			for ( Size k=1; k <= non_pack_waters.size(); ++k ) {
				Real dist2 = bb_wats[i][ii].distance_squared(non_pack_waters[k]);
				if ( dist2 < min_dist2 ) {
					min_dist2 = dist2;
					closest_bbwat = non_pack_waters[k];
				}
			}

			// compare current bb_pwat position to all other bb_pwat positions
			Size start_val = ( use_average ) ? ( i+1 ) : 1;
			for ( Size j=start_val; j <= bb_wats.size(); ++j ) {
				if ( j == i ) continue;
				for ( Size jj=1; jj <= bb_wats[j].size(); ++jj ) {
					Real dist2 = bb_wats[i][ii].distance_squared(bb_wats[j][jj]);
					if ( dist2 < min_dist2 ) {
						min_dist2 = dist2;
						closest_bbwat = bb_wats[j][jj];
					}
				}
			}

			if ( min_dist2 <= overlap_dist*overlap_dist ) {
				if ( use_average ) {
					overlap_waters.push_back( (bb_wats[i][ii]+closest_bbwat)/2.0 );
				} else {
					overlap_waters.push_back( bb_wats[i][ii] );
				}
			}
		}
	}

	TR << "number of lkball and pwat waters within cutoff of " << overlap_dist << " = " << overlap_waters.size() << std::endl;
	TR << "time spent comparing lkball overlap = " << time(NULL) - compare_start << " seconds" << std::endl << std::endl;

}


void
lkb_pwat_rotamers_setup(
	pose::Pose & pose,
	scoring::ScoreFunction const & scfxn,
	task::PackerTaskCOP task,
	rotamer_set::RotamerSetsOP rotsets,
	SetofSets & new_pwat_rotsets,
	bool exclude_exposed,
	bool use_average
)
{
	using namespace basic::options;
	using namespace interaction_graph;
	pack_scorefxn_pose_handshake( pose, scfxn );

	pose.update_residue_neighbors();

	scfxn.setup_for_packing( pose, task->repacking_residues(), task->designing_residues() );

	utility::graph::GraphOP packer_neighbor_graph = create_packer_graph( pose, scfxn, task );

	rotsets->set_task( task );
	rotsets->initialize_pose_for_rotsets_creation(pose);
	rotsets->build_rotamers( pose, scfxn, packer_neighbor_graph );
	rotsets->prepare_sets_for_packing( pose, scfxn );   /// probably only need to do this once rotsets is finalized

	// keep lkb and backbone_wats separate so can exclude self when looking for overlap
	utility::vector1< utility::vector1< Vector > > lkb_wats, bb_wats;
	utility::vector1< Vector > clash_atoms, neighbor_atoms, non_pack_waters, all_bb_pwat;

	// get all relevant bb_pwat and lkball sites
	get_clash_and_neighbors( pose, task, clash_atoms, neighbor_atoms );
	get_all_possible_water_sites( pose, task, rotsets, clash_atoms, neighbor_atoms, lkb_wats, bb_wats, non_pack_waters, all_bb_pwat, exclude_exposed );

	// find intersection of bb_pwat and lkball sites
	utility::vector1< Vector > overlap_waters;
	find_overlap( lkb_wats, bb_wats, non_pack_waters, overlap_waters, use_average );

	// cluster the intersection points
	utility::vector1< utility::vector1< Vector > > overlap_clusters;
	Real clust_rad = option[ OptionKeys::corrections::water::lkb_cluster_radius ].value();
	cluster_sites( overlap_waters, overlap_clusters, clust_rad );

	// get the centroids for each cluster
	utility::vector1< Vector > all_centroids;
	for ( Size i=1; i <= overlap_clusters.size(); ++i ) {
		all_centroids.push_back( centroid(overlap_clusters[i]) );
	}

	TR << "for " << overlap_waters.size() << " initial pwat positions: " << all_centroids.size() << " clusters at radius " << clust_rad << std::endl;

	// build new rotsets out of the lkball/pwat intersection cluster centroids
	Real rotset_cutoff = 3.0;
	cluster_sites( all_centroids, new_pwat_rotsets, rotset_cutoff );
	TR << " using cluster radius of " << rotset_cutoff << ", " << all_centroids.size() << " new pwat sites are placed in " << new_pwat_rotsets.size() << " rotamer clouds" << std::endl;

}


// @begin pack_rotamers_run
// @details as simple as possible -- runs simulated annealing and then places
// the optimal rotamers onto the backbone of the input pose.
Real
pack_pwat_rotamers_run(
	pose::Pose & pose,
	rotamer_set::FixbbRotamerSetsCOP rotsets,
	interaction_graph::AnnealableGraphBaseOP ig,
	utility::vector0< int > rot_to_pack, // defaults to an empty vector (no effect)
	utility::vector1< PointDwell > & all_rot
)
{
	using namespace ObjexxFCL::format;

	FArray1D_int bestrotamer_at_seqpos( pose.total_residue() );
	core::PackerEnergy bestenergy( 0.0 );

	pack_pwat_rotamers_run( pose, rotsets, ig, rot_to_pack, bestrotamer_at_seqpos, bestenergy, all_rot );

	// place new rotamers on input pose
	for ( uint ii = 1; ii <= rotsets->nmoltenres(); ++ii ) {
		uint iiresid = rotsets->moltenres_2_resid( ii );
		uint iibestrot = rotsets->rotid_on_moltenresidue( bestrotamer_at_seqpos( iiresid ) );
		conformation::ResidueCOP bestrot( rotsets->rotamer_set_for_moltenresidue( ii )->rotamer( iibestrot ) );
		conformation::ResidueOP newresidue( bestrot->create_residue() );
		pose.replace_residue ( iiresid, *newresidue, false );
	}
	return bestenergy;
}

/// @brief Runs simulated annealing and returns the
void
pack_pwat_rotamers_run(
	pose::Pose const & pose,
	rotamer_set::FixbbRotamerSetsCOP rotsets,
	interaction_graph::AnnealableGraphBaseOP ig,
	utility::vector0< int > rot_to_pack,
	ObjexxFCL::FArray1D_int & bestrotamer_at_seqpos,
	core::PackerEnergy & bestenergy,
	utility::vector1< PointDwell > & all_rot
)
{
	using namespace annealer;

	bool start_with_current = false;
	FArray1D_int current_rot_index( pose.total_residue(), 0 );
	bool calc_rot_freq = false;
	FArray1D< core::PackerEnergy > rot_freq( ig->get_num_total_states(), 0.0 );

	SimAnnealerBaseOP annealer = SimAnnealerBaseOP( new FixbbPwatSimAnnealer(
		rot_to_pack, bestrotamer_at_seqpos, bestenergy, start_with_current, ig,
		rotsets, current_rot_index, calc_rot_freq, rot_freq, all_rot ) );

	time_t const pwatanneal_start = time(NULL);

	PROF_START( basic::SIMANNEALING );
	annealer->run();
	PROF_STOP( basic::SIMANNEALING );

	TR << "time spent packing pwat = " << time(NULL) - pwatanneal_start << " seconds -- runtime" << std::endl;
}

} // namespace pack
} // namespace core
