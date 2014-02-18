// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/ResidueBurialFilter.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#include <protocols/simple_filters/ResidueBurialFilter.hh>
#include <protocols/simple_filters/ResidueBurialFilterCreator.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
#include <core/pack/task/TaskFactory.hh>

namespace protocols{
namespace simple_filters {

static basic::Tracer residue_burial_filter_tracer( "protocols.simple_filters.ResidueBurialFilter" );

ResidueBurialFilter::ResidueBurialFilter() :
	Filter( "ResidueBurial" ),
	residue_( "" ),
	neighbors_( 1 ),
	distance_threshold_( 8.0 ),
	task_factory_( NULL ) {}

/*SJF does anyone use this constructor? If so, we need to change target_residue to a string to make it work
ResidueBurialFilter::ResidueBurialFilter( core::Size const target_residue, core::Size const neighbors, core::Real const distance_threshold ) :
	    Filter( "ResidueBurial" ),
			residue_( target_residue ),
			neighbors_( neighbors ),
			distance_threshold_( distance_threshold ),
			task_factory_( NULL ) {}
*/
protocols::filters::FilterOP
ResidueBurialFilterCreator::create_filter() const { return new ResidueBurialFilter; }

std::string
ResidueBurialFilterCreator::keyname() const { return "ResidueBurial"; }

ResidueBurialFilter::~ResidueBurialFilter(){}

void
ResidueBurialFilter::task_factory( core::pack::task::TaskFactoryOP tf ){
	task_factory_ = tf;
}

core::pack::task::TaskFactoryOP
ResidueBurialFilter::task_factory() const{ return task_factory_; }

bool
ResidueBurialFilter::apply( core::pose::Pose const & pose ) const {
	if( task_factory() ){//taskfactory is on, iterate over all designable residues and check whether any one passes the filter
/// This looks like recursion but is really quite limited, b/c the call below to rbf.apply uses the functionality where task_factory() is off, so it goes by a different route, and could never go to a depth of more than 2
		utility::vector1< core::Size > const target_residues( protocols::rosetta_scripts::residue_packer_states( pose, task_factory(), true/*designable*/, false/*packable*/ ) );
		core::Size const total_designable( target_residues.size() );
		residue_burial_filter_tracer<<"total target residues: "<<total_designable<<". At least "<<(core::Real)total_designable * (core::Real)residue_fraction_buried()<<" should be buried to pass."<<std::endl;
		core::Size count_buried( 0 );
		foreach( core::Size const resi, target_residues ){
			ResidueBurialFilter rbf;
			rbf.residue( utility::to_string( resi ) );
			rbf.neighbors( neighbors() );
			rbf.distance_threshold( distance_threshold() );
			if( rbf.apply( pose ) ){
				count_buried++;
				if( (core::Real ) ( (core::Real) count_buried / (core::Real) total_designable ) >= residue_fraction_buried() )
					return true;
			}
		}
		return false;
	}

	core::Size const count_neighbors( compute( pose ) );

	core::Size const residue_num( core::pose::parse_resnum( residue(), pose ) );
	residue_burial_filter_tracer<<"Number of interface neighbors of residue "<<pose.residue( residue_num  ).name3()<<residue_num<<" is "<<count_neighbors<<std::endl;
	return( count_neighbors >= neighbors_ );
}

void
ResidueBurialFilter::parse_my_tag( utility::tag::TagCOP const tag, basic::datacache::DataMap & data,
		filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	residue( tag->getOption< std::string >( "pdb_num", "" ) );
	distance_threshold_ = tag->getOption<core::Real>( "distance", 8.0 );
	neighbors_ = tag->getOption<core::Size>( "neighbors", 1 );
	residue_fraction_buried( tag->getOption< core::Real >( "residue_fraction_buried", 0.0001 ) );
	if( tag->hasOption( "task_operations" ) )
		task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );

	residue_burial_filter_tracer << "ResidueBurialFilter with distance threshold of " << distance_threshold_ <<
			" around residue " << residue() << " residue_fraction_buried " << residue_fraction_buried() << " with " <<
			neighbors_ << " neighbors." << std::endl;
}

void
ResidueBurialFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	if( !task_factory() ){
		core::Size const count_neighbors( compute( pose ) );

		core::Size const residue_num( core::pose::parse_resnum( residue(), pose ) );
		out<<"Number of interface neighbors of residue "<<pose.residue( residue_num ).name3()<<residue_num<<" is "<<count_neighbors<<'\n';
	}
}

core::Real
ResidueBurialFilter::report_sm( core::pose::Pose const & pose ) const {
	if( !task_factory() ){
		core::Size const count_neighbors( compute( pose ) );

		return( count_neighbors );
	}
	else return( 0 );
}

/// @details counts the number of residues to target_residue_ across all chains in the pose, other than the one containing target_residue_
core::Size
ResidueBurialFilter::compute( core::pose::Pose const & pose ) const {
	core::Size const residue_num( core::pose::parse_resnum( residue(), pose ) );
	residue_burial_filter_tracer<<"Residue: "<<residue()<<" is serialized to: "<<residue_num<<std::endl;
	core::Size chain( 1 );
	for( ; chain <= pose.conformation().num_chains(); ++chain )
		if( pose.conformation().chain_begin( chain ) <= residue_num && pose.conformation().chain_end( chain ) >= residue_num ) break;

	core::Size const chain_begin( pose.conformation().chain_begin( chain ) );
	core::Size const chain_end( pose.conformation().chain_end( chain ) );

	residue_burial_filter_tracer<<"chain span "<<chain_begin<< " "<<chain_end<<std::endl;
	core::Size count_neighbors( 0 );
	core::conformation::Residue const res_target( pose.conformation().residue( residue_num ) );
	for( core::Size i=1; i<=pose.total_residue(); ++i ){
		if( i>=chain_begin && i<=chain_end ) continue;
		core::conformation::Residue const resi( pose.residue( i ) );
		core::Real const distance( resi.xyz( resi.nbr_atom() ).distance( res_target.xyz( res_target.nbr_atom() ) ) );
		if( distance <= distance_threshold_ ) ++count_neighbors;
	}
	return( count_neighbors);
}

filters::FilterOP
ResidueBurialFilter::clone() const {
    return new ResidueBurialFilter( *this );
}

filters::FilterOP
ResidueBurialFilter::fresh_instance() const{
    return new ResidueBurialFilter();
}

}
}
