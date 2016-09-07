// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/simple_filters/PackerNeighborGraphFilter.hh
/// @brief header file for packer neighbor graph based filter
/// @details
///
///
///
/// @author Florian Richter, floric@u.washington.edu (march 09 )

#ifndef INCLUDED_protocols_simple_filters_PackerNeighborGraphFilter_hh
#define INCLUDED_protocols_simple_filters_PackerNeighborGraphFilter_hh

// Unit Headers
#include <protocols/filters/Filter.hh>

// Project Headers
//#include <core/graph/Graph.fwd.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>


// Utility headers
//#include <utility/pointer/ReferenceCount.hh>

//// C++ headers
#include <string>
#include <set>
#include <map>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_filters {


/// @details helper class for PackerNeighborGraphFilter
class RegionalConnections {

public:

	RegionalConnections(
		std::set< core::Size >  reg1,
		std::set< core::Size >  reg2,
		core::Size required_cons
	) : region1_(std::move( reg1 )), region2_(std::move( reg2 )), required_connections_( required_cons ), num_cons_( 0 )
	{}

	std::set< core::Size > const &
	region1() const {
		return region1_; }

	std::set< core::Size > const &
	region2() const {
		return region2_; }

	bool
	enough_connections() const {
		return (required_connections_ <= num_cons_); }

	//core::Size
	//num_connections() const {
	// return num_connections_; }

	void
	reset_num_connections() const{
		num_cons_ = 0; }

	void
	check_if_connected_residues_belong_to_regions(
		core::Size res1,
		core::Size res2
	) const;

private:
	std::set< core::Size > region1_;
	std::set< core::Size > region2_;

	core::Size required_connections_;

	mutable core::Size num_cons_;

};

/// @details filter that creates a packer neighbor graph of the pose
/// @details in every apply function and returns true if this graph
/// @details satisfies a specified connectivity
class PackerNeighborGraphFilter : public protocols::filters::Filter {

public:

	PackerNeighborGraphFilter(
		core::pack::task::PackerTaskCOP task,
		core::scoring::ScoreFunctionCOP sfxn
	) : task_(std::move( task )), sfxn_(std::move( sfxn )), task_invalidated_( false ) {}


	PackerNeighborGraphFilter();

	~PackerNeighborGraphFilter() override;


	filters::FilterOP clone() const override {
		return filters::FilterOP( new PackerNeighborGraphFilter( *this ) ); }

	filters::FilterOP fresh_instance() const override {
		return filters::FilterOP( new PackerNeighborGraphFilter() ); }


	/// @brief Returns true if the given pose passes the filter, false otherwise.
	bool apply( core::pose::Pose const & pose ) const override;

	
	std::string name() const override { return "PackerNeighborGraphFilter"; }

	void
	set_task( core::pack::task::PackerTaskCOP task ){
		task_ = task; }

	/// @brief note: will overwrite if information for this residue has already been entered
	void
	set_required_connections_for_residue(
		core::Size residue,
		core::Size required_connections
	);

	/// @brief note: will increase required connections for this residue by 1
	void
	add_required_connection_for_residue(
		core::Size residue
	);

	void
	add_required_connections_between_regions(
		std::set< core::Size > const & region1,
		std::set< core::Size > const & region2,
		core::Size required_connections
	);


private:

	core::pack::task::PackerTaskCOP task_;
	core::scoring::ScoreFunctionCOP sfxn_;

	std::map< core::Size, core::Size > required_connections_per_res_;

	utility::vector1< RegionalConnections > required_connections_between_regions_;

	bool task_invalidated_;

};


typedef utility::pointer::shared_ptr< PackerNeighborGraphFilter >  PackerNeighborGraphFilterOP;
typedef utility::pointer::shared_ptr< PackerNeighborGraphFilter const >  PackerNeighborGraphFilterCOP;

} // filters
} // protocols

#endif
