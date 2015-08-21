// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author

#ifndef INCLUDED_protocols_toolbox_Cluster_hh
#define INCLUDED_protocols_toolbox_Cluster_hh


#include <core/types.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <deque>
#include <protocols/toolbox/DecoySetEvaluation.fwd.hh>

#include <utility/vector1_bool.hh>


namespace protocols {
namespace toolbox {

class ClusterBase {
public:
	ClusterBase( core::Size dim ) :
		distance_( dim, dim, 0.0 ),
		dim_( dim )
	{};

	core::Real dist( core::Size i, core::Size j ) const {
		return distance_( i, j );
	}

	core::Real & dist( core::Size i, core::Size j ) {
		return distance_( i, j);
	}

	ObjexxFCL::FArray2D< core::Real > & distance_matrix() { return distance_; }
	ObjexxFCL::FArray2D< core::Real > const& distance_matrix() const { return distance_; }

	core::Size dim() const { return dim_; }

	void print_cluster_assignment( std::ostream& out ) const;

	typedef std::deque< core::Size > Cluster;
	typedef utility::vector1 < Cluster > ClusterList;
	typedef Cluster::const_iterator IntraClusterIterator;
	typedef ClusterList::const_iterator ClusterIterator;

	core::Size size() const {
		return clusterlist_.size();
	}

	Cluster const& cluster( core::Size i ) const {
		return clusterlist_[ i ];
	}

	ClusterList const& clusterlist() const {
		return clusterlist_;
	}

	void sort_each_group_by_energy( utility::vector1< core::Real > all_energies, bool keep_center = false );
	void limit_groupsize( core::Size limit );
	void print_summary( utility::vector1< std::string > tags, utility::vector1< core::Real > all_energies );

	void show( std::ostream& out ) const;
	void read( std::istream& in );

protected:
	ClusterList clusterlist_;
	ObjexxFCL::FArray2D< core::Real > distance_;
	core::Size dim_;

};

class ClusterPhilStyle : public ClusterBase {
public:
	ClusterPhilStyle( core::Size dim, core::Real rad = 1.0 ) :
		ClusterBase( dim ),
		cluster_radius_ ( rad ),
		n_max_cluster_( 10 )
	{};

	void compute();
	void do_clustering() { compute(); };

	void set_radius( core::Real setting ) {
		cluster_radius_ = setting;
	}

	void set_n_max_cluster( core::Size setting ) {
		n_max_cluster_ = setting;
	}
private:
	core::Real cluster_radius_;
	core::Size n_max_cluster_;
};


extern std::ostream& operator<< ( std::ostream& out, ClusterBase const& cl);
extern std::ostream& operator<< ( std::ostream& out, ClusterBase::ClusterList const& cl);
extern std::istream& operator>> ( std::istream& in, ClusterBase& cl);
extern std::istream& operator>> ( std::istream& in, ClusterBase::ClusterList& cl);

class ClusterOptions {
public:
	ClusterOptions( bool assign_new_cluster_tag = false ); //set options from cmd-line

	bool assign_new_cluster_tag;
	core::Size limit_cluster_size;
	core::Real cluster_radius;
	bool keep_center;
};


// use with any iterator that points to a SilentStruct& :
// (1) SilentFileData sfd;
// cluster_silent_structs( sfd.size(), sfd.begin(), sfd.end(), ... )
//
// or
// (2) typedef vector1< SilentStructOP > SilentStructs;
// typedef toolbox::OP_const_iterator< SilentStructs::const_iterator, SilentStructs::value_type > Deref;
// SilentStructs list;
// /* fill list with decoys ... */
// cluster_silent_struct( list.size(), Deref( list.begin() ), Deref( list.end() ), ... );
//
template< typename SilentStructIterator, typename StructureContainer >
void cluster_silent_structs(
	core::Size n_decoys,
	SilentStructIterator input_decoys_begin,
	SilentStructIterator input_decoys_end,
	StructureContainer& new_structs, //provides a "push_back" method
	ClusterOptions opts
);

template< typename SilentStructIterator, typename StructureContainer >
void cluster_silent_structs( DecoySetEvaluation const& CA_set,
	SilentStructIterator input_decoys_begin,
	SilentStructIterator input_decoys_end,
	StructureContainer& new_structs,
	ClusterOptions opts
);

}
}

#endif
