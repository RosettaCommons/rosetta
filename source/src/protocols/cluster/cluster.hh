// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
///
/// @author Mike Tyka

#ifndef INCLUDED_protocols_cluster_cluster_hh
#define INCLUDED_protocols_cluster_cluster_hh

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/rosetta_scripts/PosePropertyReporter.fwd.hh>
#include <protocols/rosetta_scripts/PosePropertyReporter.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray2D.hh>

// C++ headers
#include <iostream>
#include <string>
#include <deque>
#include <vector>
#include <algorithm>

#include <utility/vector1.hh>


namespace protocols {
namespace cluster {

class GatherPosesMover;
typedef utility::pointer::shared_ptr< GatherPosesMover > GatherPosesMoverOP;
typedef utility::pointer::shared_ptr< GatherPosesMover const > GatherPosesMoverCOP;

class GatherPosesMover: public moves::Mover {
public:
	GatherPosesMover();

	void set_score_function( core::scoring::ScoreFunctionOP sfxn );
	void set_filter( core::Real filter );
	core::Real get_filter() const;

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	bool check_tag( const std::string &query_tag );

	void push_back( const core::pose::Pose & pose ) {
		poselist.push_back( pose );
	}

	void clear() {
		poselist.clear();
	}

	void set_cluster_by_protein_backbone( bool const & setting ) {  cluster_by_protein_backbone_ = setting; }

	void set_cluster_by_all_atom( bool const & setting ) {  cluster_by_all_atom_ = setting; }

	/// @brief Return the number of poses gathered by this mover
	/// This is the post-filtered number of poses.
	core::Size nposes() const { return poselist.size(); }

	std::vector< core::pose::Pose > & get_pose_list() { return poselist; }

	std::vector< std::string > const & get_tag_list() const{ return tag_list; }

	virtual core::Real
	get_distance_measure(
		const core::pose::Pose & pose1,
		const core::pose::Pose & pose2
	) const;

protected:
	std::vector< core::pose::Pose > poselist;

private:
	core::scoring::ScoreFunctionOP sfxn_;
	core::Real filter_;
	std::vector< std::string > tag_list;
	bool cluster_by_protein_backbone_;
	bool cluster_by_all_atom_;

	std::map< std::string, core::Real > template_scores;

};


// Main base class for making constraints out of groups of structures

class EnsembleConstraints;
typedef utility::pointer::shared_ptr< EnsembleConstraints > EnsembleConstraintsOP;
typedef utility::pointer::shared_ptr< EnsembleConstraints const > EnsembleConstraintsCOP;

class EnsembleConstraints: public protocols::cluster::GatherPosesMover {
public:
	EnsembleConstraints(): GatherPosesMover() {};
#ifndef BOINC // gives windows build error
	//EnsembleConstraints* clone() const = 0;
#endif
	virtual void createConstraints( std::ostream &out ) = 0;
	virtual std::string get_name() const;
};


// A super simple implementation of the above - jsut create bounds for close CA atoms.
class EnsembleConstraints_Simple;
typedef utility::pointer::shared_ptr< EnsembleConstraints_Simple > EnsembleConstraints_SimpleOP;
typedef utility::pointer::shared_ptr< EnsembleConstraints_Simple const > EnsembleConstraints_SimpleCOP;

class EnsembleConstraints_Simple: public protocols::cluster::EnsembleConstraints {
public:
	EnsembleConstraints_Simple( core::Real minimum_width = 0): EnsembleConstraints() {
		minimum_width_ = minimum_width;
	};
#ifndef BOINC // gives windows build error
	protocols::moves::MoverOP clone() const { return protocols::moves::MoverOP( new EnsembleConstraints_Simple( *this ) ) ; }
#endif
	virtual void createConstraints( std::ostream &out);
	virtual std::string get_name() const;

protected:
	core::Real minimum_width_;
};


class ClusterBase;
typedef utility::pointer::shared_ptr< ClusterBase > ClusterBaseOP;
typedef utility::pointer::shared_ptr< ClusterBase const > ClusterBaseCOP;

struct Cluster {
public:
	Cluster():
		cluster_center_( 1 ),
		group_size_(0)
	{

	}

	Cluster( int new_cluster_center ):
		cluster_center_( new_cluster_center ),
		group_size_(0)
	{
		add_member( new_cluster_center );
	}

public:

	int 							cluster_center_;
	std::deque< int > member;
	core::Size        group_size_;       // Log this seperately, because sometimes members are prunned, but its still nice to know what the original size would have been

public:

	int  get_cluster_center() const { return cluster_center_; }
	//void set_cluster_center( int new_cluster_center ) { cluster_center_ = new_cluster_center; }

	// Actual removing or addition of members
	void add_member( int new_member ){
		member.push_back( new_member );
		group_size_ ++;
	}
	void add_member_front( int new_member ){
		member.push_front( new_member );
		group_size_ ++;
	}
	void remove_member( int old_member ){
		erase( old_member );
		group_size_ --;
	}

	// Editing operations (they are more crude and dont affect the group_size count )
	void push_back( int new_member ){ member.push_back( new_member ); }
	void push_front( int new_member ){ member.push_front( new_member ); }
	int&  operator [] ( int index ){ return member[index]; }
	int  operator [] ( int index ) const { return member[index]; }
	core::Size size() const { return member.size(); }
	void clear(){ member.clear(); }
	void erase( core::Size j ){
				std::deque< int >::iterator it = member.begin();
				it += j;
				member.erase(it);
	}
	void shuffle();

	core::Size group_size(){ return group_size_; }
};

class ClusterBase: public GatherPosesMover {
public:
	ClusterBase();
	virtual std::string get_name() const;
	void set_cluster_radius( core::Real cluster_radius );
	void set_population_weight( core::Real population_weight );

	core::Real get_cluster_radius();
	core::Real get_median_rms(){ return median_rms_; }

	void calculate_distance_matrix();

	void add_structure( core::pose::Pose & pose );


	// PostProcessing ---------------------------------------------------------

	void remove_highest_energy_member_of_each_group();
	void sort_each_group_by_energy( );
	void sort_groups_by_energy( );
	void remove_singletons();
	void export_only_low( bool value = false ){ export_only_low_ = value; };
	void limit_groupsize( int limit = -1);
	void limit_groupsize( core::Real percent_limit );
	void random_limit_groupsize( core::Real percent_limit );
	void limit_groups( int limit = -1);
	void limit_total_structures( int limit = -1);
	void clean_pose_store();

	// Printing --------------------------------------------------------------

	void print_summary();
	void print_raw_numbers();
	void print_cluster_assignment();
	void print_cluster_PDBs( std::string prefix );
	void print_clusters_silentfile( std::string prefix );
	std::vector< core::pose::PoseOP > return_lowest_poses_in_clusters();
	std::vector< core::pose::PoseOP > return_top_poses_in_clusters(core::Size count);
	void create_constraints( std::string prefix, EnsembleConstraints &constraint_maker);

  std::vector < Cluster >  const & get_cluster_list() const{ return clusterlist; }

protected:
	ObjexxFCL::FArray2D< core::Real >    distance_matrix;
	std::vector < Cluster >    clusterlist;

	bool  export_only_low_;
	core::Real median_rms_;
	core::Real population_weight_;
	core::Real cluster_radius_;
};


class ClusterPhilStyle;
typedef utility::pointer::shared_ptr< ClusterPhilStyle > ClusterPhilStyleOP;
typedef utility::pointer::shared_ptr< ClusterPhilStyle const > ClusterPhilStyleCOP;

class ClusterPhilStyle: public ClusterBase {
public:
	ClusterPhilStyle();
	virtual ~ClusterPhilStyle() {};
	protocols::moves::MoverOP clone() const { return protocols::moves::MoverOP( new ClusterPhilStyle( *this ) ) ; }
	virtual std::string get_name() const;
	virtual void do_clustering( core::Size max_total_cluster );

	// this ensures every structure is in the cluster to who's cluster center it is most similar too
	void do_redistribution();
};

class ClusterPhilStyle_Loop: public ClusterPhilStyle {
public:
	ClusterPhilStyle_Loop( protocols::loops::Loops loop_def )
		: loop_def_(loop_def)
	{}

	virtual ~ClusterPhilStyle_Loop() {}
	virtual std::string get_name() const;
	protocols::moves::MoverOP clone() const {
		return protocols::moves::MoverOP( new ClusterPhilStyle_Loop( *this ) );
	}

	virtual core::Real
	get_distance_measure(
		const core::pose::Pose & pose1,
		const core::pose::Pose & pose2
	) const;

private:
	protocols::loops::Loops loop_def_;
};

class ClusterPhilStyle_PoseReporter : public ClusterPhilStyle {
public:
	ClusterPhilStyle_PoseReporter( protocols::rosetta_scripts::PosePropertyReporterOP reporter )
		: reporter_(reporter)
	{}

	virtual ~ClusterPhilStyle_PoseReporter() {}
	virtual std::string get_name() const;
	protocols::moves::MoverOP clone() const {
		return protocols::moves::MoverOP( new ClusterPhilStyle_PoseReporter( *this ) );
	}

	virtual core::Real
	get_distance_measure(
		const core::pose::Pose & pose1,
		const core::pose::Pose & pose2
	) const;

private:
	protocols::rosetta_scripts::PosePropertyReporterOP reporter_;
};

class AssignToClustersMover;
typedef utility::pointer::shared_ptr< AssignToClustersMover > AssignToClustersMoverOP;
typedef utility::pointer::shared_ptr< AssignToClustersMover const > AssignToClustersMoverCOP;

class AssignToClustersMover: public GatherPosesMover {
public:
	AssignToClustersMover( ClusterBaseOP cluster_base );
#ifndef BOINC // gives windows build error
	protocols::moves::MoverOP clone() const {
		return protocols::moves::MoverOP( new AssignToClustersMover( *this ) );
	}
#endif
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;


protected:
	ClusterBaseOP cluster_base_;
};

} //namespace cluster
} //namespace protocols

#endif
