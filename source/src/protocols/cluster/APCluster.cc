// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/cluster/APCluster.cc
///
/// @brief
/// @author Ian W. Davis


#include <protocols/cluster/APCluster.hh>

#include <basic/Tracer.hh>

#include <utility/exit.hh>

#include <cstdio>

static thread_local basic::Tracer TR( "protocols.cluster.APCluster" );

namespace protocols {
namespace cluster {


/// @details Less-than comparator works backwards on similarity values,
/// to effectively implement a min-heap instead of a max-heap.
bool Exemplar::min_heap(Exemplar a, Exemplar b)
{
	return a.s_ik > b.s_ik;
}


/// @details If this point has more than max_sims similarities already stored,
/// the lowest one (most negative) will be discarded.
///
/// There is currently no protection against adding s(i,k) twice,
/// which will not be caught and will screw up the computation royally.
void DataPoint::add_similarity(core::Size k, core::Real s_ik, core::Size max_sims)
{
	if ( k == i ) {
		s_kk = s_ik;
		curr_exemplar = 0; // used as a flag that s_kk has been set
	} else {
		Exemplar ex(k, s_ik);
		candidates.push_back(ex); // add element
		std::push_heap( candidates.begin(), candidates.end(), Exemplar::min_heap ); // sort smallest s_ik to front
		while ( candidates.size() > max_sims ) {
			std::pop_heap( candidates.begin(), candidates.end(), Exemplar::min_heap ); // sort smallest s_ik to back
			candidates.pop_back(); // remove element
		}
	}
}


APCluster::APCluster(core::Size total_pts, core::Size max_sims_per_pt /*= 0*/):
	utility::pointer::ReferenceCount(),
	pts_(),
	is_frozen_(false)
{
	if ( max_sims_per_pt > 0 ) max_sims_per_pt_ = max_sims_per_pt;
	else max_sims_per_pt_ = total_pts;

	for ( core::Size i = 1; i <= total_pts; ++i ) {
		pts_.push_back( DataPoint(i) );
	}
}


APCluster::~APCluster() {}


/// @details Adding s(i,j) is not the same as adding s(j,i) -- you must do both if you want symmetry.
///
/// There is currently no protection against adding s(i,k) twice,
/// which will not be caught and will screw up the computation royally.
void APCluster::set_sim(core::Size i, core::Size k, core::Real sim)
{
	// If we've already filled in the DataPoint.candidate_for arrays,
	// we'll have to go back and recompute them before we can cluster again.
	// As a special case, we can change the self preferences without trouble.
	if ( i != k && is_frozen_ ) is_frozen_ = false;
	pts_[i].add_similarity(k, sim, max_sims_per_pt_);
}


/// @details
/// @param maxits   maximum number of iterations, period; 100 - 4000 reasonable.
/// @param convits  terminate after clusters don't change for this many iterations, 10 - 400 reasonable.
/// @param lambda   damping factor, 0.50 - 0.95 is a reasonable range
/// @return true if the clustering converges and false otherwise.
/// Failure to converge is not necessarily a problem -- you may be close to a good solution.
/// Try increasing maxits and/or making lambda closer to 1 if you want convergence.
bool APCluster::cluster(core::Size maxits, core::Size convits, core::Real lambda)
{
	runtime_assert( 0 <= lambda && lambda < 1 );
	if ( !is_frozen_ ) freeze();
	reinitialize();

	bool is_converged = false;
	core::Size no_change_its = 0;
	//core::Real last_net_sim = -1e99;
	core::Real max_net_sim = -1e99;
	for ( core::Size itr = 1; itr <= maxits; ++itr ) {
		update_r_ik(lambda);
		update_a_ik(lambda);
		core::Size const changes = assign_exemplars();
		// Check for early termination if exemplar assignments haven't changed for several cycles.
		// I find that sometimes the run gets stuck with a few unassigned points (which count as changes)
		// and if you continue on, you shoot off to a less-optimal solution.  Here, "a few" is < 5%.
		// This heuristic seems to hacky though, so instead I'm just going to track net_sim
		// and return the solution with the lowest value of that.
		core::Real const net_sim = get_net_sim();
		if ( net_sim > max_net_sim ) {
			save_best_exemplars();
			max_net_sim = net_sim;
		}
		if ( changes == 0 /*|| (net_sim == last_net_sim && changes <= pts_.size()/20 )*/ ) {
			no_change_its += 1;
			if ( no_change_its >= convits ) {
				is_converged = true;
				break;
			}
		} else {
			no_change_its = 0;
		}
		TR.Debug << "Iteration " << itr << ", " << get_num_exemplars() << " clusters, net_sim " << net_sim << ", "
			<< changes << " changes, stable for " << no_change_its << std::endl;
		//last_net_sim = net_sim;
	}
	restore_best_exemplars();
	core::Real const net_sim = get_net_sim();
	TR << "Finished, " << get_num_exemplars() << " clusters, net_sim " << net_sim << std::endl;
	return is_converged;
}


core::Size APCluster::get_num_exemplars() const
{
	core::Size count = 0;
	for ( core::Size i = 1, N = pts_.size(); i <= N; ++i ) {
		DataPoint const & p = pts_[i];
		if ( p.i == p.curr_exemplar ) count += 1;
	}
	return count;
}


void APCluster::get_all_exemplars(utility::vector1< core::Size > & exemplars) const
{
	exemplars.clear();
	for ( core::Size i = 1, N = pts_.size(); i <= N; ++i ) {
		DataPoint const & p = pts_[i];
		if ( p.i == p.curr_exemplar ) exemplars.push_back(i);
	}
}


void APCluster::get_cluster_for(core::Size k, utility::vector1< core::Size > & cluster) const
{
	cluster.clear();
	for ( core::Size i = 1, N = pts_.size(); i <= N; ++i ) {
		DataPoint const & p = pts_[i];
		if ( k == p.curr_exemplar ) cluster.push_back(i);
	}
}


core::Real APCluster::get_net_sim() const
{
	core::Real sum = 0;
	for ( core::Size i = 1, N = pts_.size(); i <= N; ++i ) {
		DataPoint const & p = pts_[i];
		if ( p.curr_exemplar == i ) {
			sum += p.s_kk; // Should this be included?  Yes.
			continue;
		}
		for ( core::Size j = 1, j_end = p.candidates.size(); j <= j_end; ++j ) {
			Exemplar const & e = p.candidates[j];
			if ( p.curr_exemplar == e.k ) {
				sum += e.s_ik;
				break;
			}
		}
	}
	return sum;
}


void APCluster::freeze()
{
	if ( is_frozen_ ) return;

	for ( core::Size i = 1, N = pts_.size(); i <= N; ++i ) {
		DataPoint & p = pts_[i];
		p.candidate_for.clear();
	}
	for ( core::Size i = 1, N = pts_.size(); i <= N; ++i ) {
		DataPoint & p = pts_[i];
		for ( core::Size j = 1, j_end = p.candidates.size(); j <= j_end; ++j ) {
			// Tell point j that s(i,j) exists -- point i already knows by traversing candidates
			Exemplar & e = p.candidates[j];
			pts_[e.k].candidate_for.push_back(&e);
		}
	}

	is_frozen_ = true;
}


/// @details Prepares the data points for another clustering run.
/// Does not erase similarity data, etc.
void APCluster::reinitialize()
{
	for ( core::Size i = 1, N = pts_.size(); i <= N; ++i ) {
		DataPoint & p = pts_[i];
		p.r_kk = 0;
		p.a_kk = 0;
		p.curr_exemplar = 0;
		p.best_exemplar = 0;
		for ( core::Size j = 1, j_end = p.candidates.size(); j <= j_end; ++j ) {
			Exemplar & e = p.candidates[j];
			e.r_ik = 0;
			e.a_ik = 0;
		}
	}
}


void APCluster::update_r_ik(core::Real lambda)
{
	for ( core::Size i = 1, N = pts_.size(); i <= N; ++i ) {
		DataPoint & p = pts_[i];

		// Compute the max and runner-up just once per point.
		// Reduces time complexity from O(N^3) to O(N^2)
		core::Real max1 = p.a_kk + p.s_kk, max2 = p.a_kk + p.s_kk;
		core::Size max1_k = i/*, max2_k = i*/;
		for ( core::Size kk = 1, kk_end = p.candidates.size(); kk <= kk_end; ++kk ) {
			Exemplar const & e_kk = p.candidates[kk];
			core::Real const as = e_kk.a_ik + e_kk.s_ik;
			if ( as >= max2 ) {
				if ( as >= max1 ) {
					max2 = max1;
					//max2_k = max1_k;  // set but never used ~Labonte
					max1 = as;
					max1_k = e_kk.k;
				} else {
					max2 = as;
					//max2_k = e_kk.k;  // set but never used ~Labonte
				}
			}
		}

		// General update rule for r(i,k)
		// r(i,k) = s(i,k) - max{ a(i,k') + s(i,k') }, k' != k
		for ( core::Size k = 1, k_end = p.candidates.size(); k <= k_end; ++k ) {
			Exemplar & e = p.candidates[k];

			//core::Real max_as = p.a_kk + p.s_kk;
			//for(core::Size kk = 1, kk_end = p.candidates.size(); kk <= kk_end; ++kk) {
			// Exemplar const & e_kk = p.candidates[kk];
			// core::Real const as = e_kk.a_ik + e_kk.s_ik;
			// if( as >= max_as && k != kk ) max_as = as;
			//}
			// Global max unless that was from this k, then runner up.
			core::Real const max_as = ( max1_k == e.k ? max2 : max1 );

			core::Real const new_r_ik = e.s_ik - max_as;
			e.r_ik = (lambda)*e.r_ik + (1-lambda)*new_r_ik;
		}
		// Special update rule for r(k,k)
		{
			//core::Real max_as = -1e99;
			//for(core::Size kk = 1, kk_end = p.candidates.size(); kk <= kk_end; ++kk) {
			// Exemplar const & e_kk = p.candidates[kk];
			// core::Real const as = e_kk.a_ik + e_kk.s_ik;
			// if( as >= max_as ) max_as = as;
			//}
			core::Real const max_as = ( max1_k == i ? max2 : max1 );

			core::Real const new_r_kk = p.s_kk - max_as;
			p.r_kk = (lambda)*p.r_kk + (1-lambda)*new_r_kk;
		}
	}
}


void APCluster::update_a_ik(core::Real lambda)
{
	// Precalculate the sum term to reduce complexity from O(N^3) to O(N^2)
	for ( core::Size i = 1, N = pts_.size(); i <= N; ++i ) {
		DataPoint & p = pts_[i];
		p.sum = 0;
		for ( core::Size ii = 1, ii_end = p.candidate_for.size(); ii <= ii_end; ++ii ) {
			Exemplar const & e_ii = *( p.candidate_for[ii] );
			// Can't ever have k == i because the diagonal elements have their own vars
			if ( e_ii.r_ik > 0 /*&& e_ii.k != i*/ ) p.sum += e_ii.r_ik;
		}
		runtime_assert( p.sum >= 0 );
	}
	for ( core::Size i = 1, N = pts_.size(); i <= N; ++i ) {
		DataPoint & p = pts_[i];
		// General update rule for a(i,k)
		// a(i,k) = min{ 0, r(k,k) + sum{ max{0, r(i',k)} } }, i' != i and i' != k
		for ( core::Size k = 1, k_end = p.candidates.size(); k <= k_end; ++k ) {
			Exemplar & e = p.candidates[k];
			DataPoint const & p_k = pts_[e.k];

			//core::Real sum = 0;
			//for(core::Size ii = 1, ii_end = p_k.candidate_for.size(); ii <= ii_end; ++ii) {
			// Exemplar const & e_ii = *( p_k.candidate_for[ii] );
			// if( e_ii.r_ik > 0 && e_ii.k != i && e_ii.k != e.k ) sum += e_ii.r_ik;
			//}
			core::Real sum = p_k.sum;
			// We know i' == k has been excluded from the sum because that would be r(k,k),
			// and diagonal terms are stored in their own vars in the DataPoints rather than in Exemplars.
			// If i' == i, then r(i',k) -> r(i,k), which we have access to as e.r_ik.
			if ( e.r_ik > 0 ) sum -= e.r_ik;
			runtime_assert( sum >= 0 );

			core::Real new_a_ik = p_k.r_kk + sum;
			if ( new_a_ik > 0 ) new_a_ik = 0;
			e.a_ik = (lambda)*e.a_ik + (1-lambda)*new_a_ik;
		}
		// Special update rule for a(k,k)
		{
			core::Real const new_a_kk = p.sum;
			p.a_kk = (lambda)*p.a_kk + (1-lambda)*new_a_kk;
		}
		// It appears that a(i,k) <= 0 while a(k,k) >= 0 ... strange.
	}
}


core::Size APCluster::assign_exemplars()
{
	core::Size changes = 0;

	// The paper claims that the exemplar for each point i is the point k
	// that maximizes a(i,k)+r(i,k), where i and k may be equal.
	// HOWEVER, this can lead to the case where a point k is the exemplar
	// for other point and yet not for itself.
	// Also, it's at odds with what's actually in their MATLAB code -- see below.

	//for(core::Size i = 1, N = pts_.size(); i <= N; ++i) {
	// DataPoint & p = pts_[i];
	// core::Size new_exemplar = i;
	// core::Real max_ar = p.a_kk + p.r_kk;
	// for(core::Size k = 1, k_end = p.candidates.size(); k <= k_end; ++k) {
	//  Exemplar const & e = p.candidates[k];
	//  core::Real const ar = e.a_ik + e.r_ik;
	//  if( ar > max_ar ) {
	//   new_exemplar = e.k;
	//   max_ar = ar;
	//  }
	// }
	// if( new_exemplar != p.curr_exemplar ) changes += 1;
	// p.curr_exemplar = new_exemplar;
	//}

	// The MATLAB code says that exemplars are points for which a(k,k)+r(k,k) > 0,
	// and all other points i are assigned to the exemplar they are closest to,
	// i.e. the exemplar k that maximizes s(i,k)

	// 1. Find exemplars
	for ( core::Size i = 1, N = pts_.size(); i <= N; ++i ) {
		DataPoint & p = pts_[i];
		if ( p.a_kk + p.r_kk > 0 ) {
			if ( p.curr_exemplar != i ) changes += 1;
			p.curr_exemplar = i;
		} else if ( p.curr_exemplar == i ) {
			changes += 1;
			p.curr_exemplar = 0; // to avoid confusion below if this point was an exemplar last time
		}
	}
	// 2. Assign other points to nearest exemplar.
	for ( core::Size i = 1, N = pts_.size(); i <= N; ++i ) {
		DataPoint & p = pts_[i];
		if ( p.curr_exemplar == i ) continue; // already an exemplar itself
		core::Size new_exemplar = 0;
		core::Real max_s = -1e99;
		for ( core::Size k = 1, k_end = p.candidates.size(); k <= k_end; ++k ) {
			Exemplar const & e = p.candidates[k];
			DataPoint const & p_k = pts_[e.k];
			if ( p_k.curr_exemplar != p_k.i ) continue; // not an exemplar
			if ( e.s_ik > max_s ) {
				max_s = e.s_ik;
				new_exemplar = e.k; // == p_k.i
			}
		}
		if ( new_exemplar != p.curr_exemplar && p.curr_exemplar != 0 ) changes += 1;
		p.curr_exemplar = new_exemplar;
	}
	// 3. Mop up unassigned points by putting them in their own cluster.
	// This is possible in sparse cases; shouldn't be in dense cases.
	// Due to the logic above, these points are ALWAYS counted as changes,
	// even if they ended up in this predicament during the last cycle too.
	for ( core::Size i = 1, N = pts_.size(); i <= N; ++i ) {
		DataPoint & p = pts_[i];
		if ( p.curr_exemplar != 0 ) continue; // already has an exemplar
		//changes += 1; already counted in other steps
		p.curr_exemplar = i;
	}


	return changes;
}


void APCluster::save_best_exemplars()
{
	for ( core::Size i = 1, N = pts_.size(); i <= N; ++i ) {
		DataPoint & p = pts_[i];
		p.best_exemplar = p.curr_exemplar;
	}
}


void APCluster::restore_best_exemplars()
{
	for ( core::Size i = 1, N = pts_.size(); i <= N; ++i ) {
		DataPoint & p = pts_[i];
		p.curr_exemplar = p.best_exemplar;
	}
}


// Helper functions for IO
// I'm not too worried about efficiency as files should be block-buffered by the OS by default.
template<typename T> inline void write1(FILE /*const*/ * f, T t)
{ if ( fwrite(&t, sizeof(T), 1, f) != 1 ) { assert(false); } }
template<typename T> inline void read1(FILE /*const*/ * f, T & t)
{ if ( fread(&t, sizeof(T), 1, f) != 1 ) { assert(false); } }


bool APCluster::save_binary(std::string const & filename) const
{
	using namespace std;
	FILE* f = fopen( filename.c_str(), "wb" );
	if ( f == NULL ) return false;
	core::Size const N = pts_.size(); write1(f, N);
	for ( core::Size i = 1; i <= N; ++i ) {
		DataPoint const & p = pts_[i];
		write1(f, p.i);
		write1(f, p.s_kk);
		write1(f, p.curr_exemplar);
		core::Size const k_end = p.candidates.size(); write1(f, k_end);
		for ( core::Size k = 1; k <= k_end; ++k ) {
			Exemplar const & e = p.candidates[k];
			write1(f, e.k);
			write1(f, e.s_ik);
		}
	}
	write1(f, max_sims_per_pt_);
	fclose(f);
	return true;
}


bool APCluster::load_binary(std::string const & filename)
{
	using namespace std;
	FILE* f = fopen( filename.c_str(), "rb" );
	if ( f == NULL ) return false;
	pts_.clear();
	core::Size N; read1(f, N);
	for ( core::Size i = 1; i <= N; ++i ) {
		pts_.push_back( DataPoint(i) );
		DataPoint & p = pts_[i];
		read1(f, p.i);
		read1(f, p.s_kk);
		read1(f, p.curr_exemplar);
		core::Size k_end; read1(f, k_end);
		for ( core::Size k = 1; k <= k_end; ++k ) {
			Exemplar e(0,0);
			read1(f, e.k);
			read1(f, e.s_ik);
			p.candidates.push_back(e);
		}
	}
	read1(f, max_sims_per_pt_);
	is_frozen_ = false;
	fclose(f);
	return true;
}


} // namespace cluster
} // namespace protocols
