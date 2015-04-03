// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/cluster/APCluster.hh
///
/// @brief
/// @author Ian W. Davis


#ifndef INCLUDED_protocols_cluster_APCluster_hh
#define INCLUDED_protocols_cluster_APCluster_hh

#include <protocols/cluster/APCluster.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <core/types.hh>

#include <utility/vector1_bool.hh>


namespace protocols {
namespace cluster {


/// @brief Data structure for one similarity measurement (s_ik) for affinity propagation clustering.
///
/// @details There will be one instance of this class for each (finite) similarity between two input points,
/// up to a maximum of N*N instances if the similarity matrix is fully populated.
///
class Exemplar //: public utility::pointer::ReferenceCount
{
public:

	Exemplar(core::Size k_in, core::Real s_ik_in):
		k(k_in),
		s_ik(s_ik_in),
		r_ik(0),
		a_ik(0)
	{}

	~Exemplar() {}

	/// @brief "Less than" (actually greater than) comparator for making a heap of exemplars
	static bool min_heap(Exemplar a, Exemplar b);

public:
	core::Size k; //< label for this data point, 1 .. N
	core::Real s_ik; //< similarity: higher value = point k is a better exemplar for point i
	core::Real r_ik; //< responsibility: higher value = point k is a better exemplar for point i relative to other possible exemplars for i
	core::Real a_ik; //< availability: higher value = point k would be a good exemplar for many points i

}; // Exemplar


/// @brief Data structure for one input data point for affinity propagation clustering.
///
/// @details There should be one instance of this class for each input point.
/// Fields are public because it's a glorified struct -- clients shouldn't use this directly.
///
class DataPoint //: public utility::pointer::ReferenceCount
{
public:

	DataPoint(core::Size i_in):
		i(i_in),
		s_kk(0),
		r_kk(0),
		a_kk(0),
		sum(0),
		curr_exemplar(i_in), // != 0 used as a flag that s_kk hasn't been set yet
		best_exemplar(0),
		candidates(),
		candidate_for()
	{}

	~DataPoint() {}

	/// @brief Set similarity s(i,k), the suitability of point k to be an exemplar for this point.
	void add_similarity(core::Size k, core::Real s_ik, core::Size max_sims);

	bool is_set_s_kk() const { return curr_exemplar == 0; }

public:
	core::Size i; //< label for this data point, 1 .. N
	core::Real s_kk; //< self-similarity, aka "preference".  Higher values = more likely to be a cluster center (exemplar).
	core::Real r_kk; //< self-responsibility
	core::Real a_kk; //< self-availability
	core::Real sum; //< = sum(e.r_ik for e in candidate_for if e.r_ik > 0); a.k.a. the update value for a(k,k); cached for efficiency
	core::Size curr_exemplar; //< label of point that was best exemplar for this in last round
	core::Size best_exemplar; //< exemplar from the best round (largest net_sim) in case we want to go back to it
	utility::vector1< Exemplar > candidates; //< candidate exemplars: s(i,k) exists with i = this.  Organized as a max-heap using STL functions.
	utility::vector1< Exemplar* > candidate_for; //< reverse lookup for s(i,k), with k = this

}; // DataPoint


/// @brief Public interface for doing affinity propagation clustering.
///
/// @details Based on Frey and Dueck, "Clustering by Passing Messages Between Data Points", Science 315 (2007).
/// Useful for choosing a set of representative data points (exemplars)
/// out of a large set (e.g. all decoys from a large Rosetta run)
/// given a measure of similarity (e.g. RMSD, maxsub, GDT, ...).
///
/// As I understand it, this procedures tries to maximize the sum of similarities between
/// each data point and its exemplar, while balancing that with total number of clusters.
/// Reasonable measures of similarity would be negative RMSD, log-likelihoods,
/// or squared distance (i.e. squared error), depending on what the points represent.
/// Note there is no requirement for symmetry:  s(i,j) need not equal s(j,i).
/// The self-similarity s(k,k) ("preference") for each point controls the likelihood it will be selected as an exemplar,
/// and thus indirectly controls the total number of clusters.
/// There is no way to directly specify a specific number of clusters.
/// The authors suggest that using the median of all other similarities will give a moderate number of clusters,
/// and using the minimum of the other similaries will give a small number of clusters.
///
/// This implementation is designed for clustering very large numbers of data points
/// with sparse similarity [ s(i,k) = -Inf for most i,k ].
/// Similarities for each input point are kept in a heap so that you can limit to only the L highest for each.
/// (This scheme is quite likely to break symmetry, as some points will have more close neighbors than others.)
/// Alternately, you may choose to do your own pre-filtering and only enter the G globally highest similarities
/// between any points in the data set.
/// Run time (per cycle) is linear in the number of similarities, or O(N^2) in the limit of a dense similarity matrix.
///
/// I follow the conventions of the original paper, where "i" is the index of some generic data point,
/// and "k" is the index of a data point being considered as an exemplar (cluster center).
///
class APCluster : public utility::pointer::ReferenceCount
{
public:

	/// @brief Create new clustering class for total_pts input data points.
	/// Optionally, set a limit on the number of similarities stored per point.
	APCluster(core::Size total_pts, core::Size max_sims_per_pt = 0);
	virtual ~APCluster();

	/// @brief How appropriate is k as an exemplar for i?
	virtual void set_sim(core::Size i, core::Size k, core::Real sim);
	/// @brief Run the actual clustering algorithm.
	virtual bool cluster(core::Size maxits, core::Size convits, core::Real lambda);

	virtual core::Size num_pts() const { return pts_.size(); }
	/// @brief Return the index of the point that is the exemplar for point i.
	virtual core::Size get_exemplar_for(core::Size i) const { return pts_[i].curr_exemplar; }
	/// @brief The number of exemplars selected (number of clusters).
	/// Monotonically related to the self-preferences s(k,k).
	virtual core::Size get_num_exemplars() const;
	/// @brief Return the indices of data points chosen as exemplars (cluster centers).
	virtual void get_all_exemplars(utility::vector1< core::Size > & exemplars) const;
	/// @brief Returns the indices of points with the specified exemplar k.
	/// Note that k is the index of an (input) data point that was chosen as an exemplar,
	/// not some "cluster index" between 1 and get_num_exemplars().
	virtual void get_cluster_for(core::Size k, utility::vector1< core::Size > & cluster) const;
	/// @brief The sum of similarities s(i,k) between every data point i and its exemplar k,
	/// plus the self preferences of the exemplars.
	/// The algorithm *should* minimize this value -- if it dips and climbs again, increase lambda.
	virtual core::Real get_net_sim() const;

	/// @brief Saves the (sparse) similarity matrix and current cluster assignments (if any),
	/// but not the accumulated evidence from the last clustering [ r(i,k) and a(i,k) ].
	/// File format is custom binary and is not portable (host endian-ness).
	virtual bool save_binary(std::string const & filename) const;
	/// @brief Wipes all currently held data and reads in similarity values and cluster assignments.
	/// Afterwards, points may be re-clustered with different parameters if desired.
	/// File format is custom binary and is not portable (host endian-ness).
	virtual bool load_binary(std::string const & filename);

protected:
	virtual void freeze();
	virtual void reinitialize();
	virtual void update_r_ik(core::Real lambda);
	virtual void update_a_ik(core::Real lambda);
	virtual core::Size assign_exemplars();
	virtual void save_best_exemplars();
	virtual void restore_best_exemplars();

private:
	utility::vector1< DataPoint > pts_; //< the data points to be clustered
	core::Size max_sims_per_pt_; //< if more than this many similarities for some point i, discard the lowest ones
	bool is_frozen_; //< have the DataPoint.candidate_for vectors been filled in yet?

}; // APCluster


} // namespace cluster
} // namespace protocols

#endif // INCLUDED_protocols_cluster_APCluster_HH
