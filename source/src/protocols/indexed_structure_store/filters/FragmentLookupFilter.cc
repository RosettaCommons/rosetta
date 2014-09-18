// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ./src/protocols/indexed_structure_store/filters/FragmentLookupFilter.cc
/// @brief FragmentLookupFilter class implementation.
/// @detailed
/// @author Alex Ford (fordas@uw.edu)
//
#include <basic/Tracer.hh>

#include <utility/tag/Tag.hh>
#include <utility/json_spirit/json_spirit.h>

#include <core/pose/Pose.hh>

#include <core/indexed_structure_store/FragmentLookup.hh>
#include <core/indexed_structure_store/StructureStoreManager.hh>

#include <protocols/indexed_structure_store/filters/FragmentLookupFilter.hh>
#include <protocols/indexed_structure_store/filters/FragmentLookupFilterCreator.hh>

#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/Filter.hh>

#include <iterator>

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

#include <time.h>

namespace protocols {
namespace indexed_structure_store {
namespace filters {

static thread_local basic::Tracer TR( "protocols.indexed_structure_store.filters.FragmentLookupFilter" );

using namespace core::indexed_structure_store;

FragmentLookupFilter::FragmentLookupFilter() :
  target_lookup_(NULL)
{
}

FragmentLookupFilter::FragmentLookupFilter(
			std::string lookup_name,
			LookupMode mode) :
			target_lookup_(StructureStoreManager::get_instance()->load_fragment_lookup(lookup_name)),
			lookup_mode_(mode),
			target_chain_(0),
			threshold_(0),
			b_target_chain_(false)
{
}

FragmentLookupFilter::FragmentLookupFilter(
			std::string lookup_name,
			std::string store_path,
			LookupMode mode,
			core::Size target_chain,
			core::Size threshold,
			bool b_target_chain) :
			target_lookup_(StructureStoreManager::get_instance()->load_fragment_lookup(lookup_name, store_path)),
			lookup_mode_(mode),
			target_chain_(target_chain),
			threshold_(threshold),
			b_target_chain_(b_target_chain)
{
}

FragmentLookupFilter::FragmentLookupFilter( FragmentLookupFilter const & rval ) :
	protocols::filters::Filter(rval),
  target_lookup_(rval.target_lookup_)
{
}

core::Size FragmentLookupFilter::compute( Pose const & pose ) const
{
	  using namespace core::indexed_structure_store;


	clock_t start_t = clock();

	std::vector<FragmentLookupResult> lookup_result;
	std::vector<core::Size> lookup_residue;

	core::pose::Pose target_pose = pose;
	TR << boost::format("apply( pose=<%s residues> )") % target_pose.n_residue() << std::endl;
	//Apply to a particular chain

	if (b_target_chain_){
		TR <<  "Target chain: " << target_chain_ << "in a pose with #Chains: " << pose.conformation().num_chains() << std::endl;
		if(target_chain_ > pose.conformation().num_chains() ){
			utility_exit_with_message(" FragmentLookupFilter invalid chain" );
		}
		target_pose = *pose.split_by_chain(target_chain_);
		TR << boost::format("apply mod by chain! (Now pose=<%s residues> )") % target_pose.n_residue() << std::endl;
	}

	if(lookup_mode_ == First)
	{
		target_lookup_->lookup_pose_fragments(target_pose, std::back_inserter(lookup_result), std::back_inserter(lookup_residue));
	}
	else
	{
		target_lookup_->lookup_closest_pose_fragments(target_pose, std::back_inserter(lookup_result), std::back_inserter(lookup_residue));
	}

	float elapsed_seconds = (clock() - start_t) / CLOCKS_PER_SEC;

	TR << boost::format("lookup returned results: %s seconds: %s") % lookup_result.size() % elapsed_seconds << std::endl;

	if ( TR.Debug.visible() )
	{
		TR.Debug << "\n";

		for (core::Size i = 0; i < lookup_result.size(); i++)
		{
			TR.Debug << boost::format("residue: %s result: %s") % lookup_residue[i] % lookup_result[i].found_match << "\n";
		}

		TR.Debug << std::endl;
	}

	cached_lookup_result_.clear();

	core::Size num_failed_matches=0;
	for (core::Size i = 0; i < lookup_result.size(); i++)
	{
		if(!lookup_result[i].found_match){
			num_failed_matches+=1;
		}
		cached_lookup_result_[lookup_residue[i]] = lookup_result[i];
	}

	return num_failed_matches;
}



void FragmentLookupFilter::report( std::ostream & os, core::pose::Pose const & ) const
{
	using namespace core::indexed_structure_store;
	typedef std::map<core::Size, FragmentLookupResult> result_t;
	using utility::json_spirit::Value;
	using utility::json_spirit::Object;
	using utility::json_spirit::Pair;

	Object result_container;
	core::Size num_failed_matches=0;
	BOOST_FOREACH(result_t::value_type &r, cached_lookup_result_)
	{
		Object result_object;

		result_object.push_back(Pair("found_match", r.second.found_match));

		if(!r.second.found_match){
			num_failed_matches+=1;
		}

		if (r.second.found_match)
		{
			result_object.push_back(Pair("match_score", r.second.match_score));
			result_object.push_back(Pair("match_index", static_cast<boost::uint64_t>(r.second.match_index)));
			result_object.push_back(Pair("match_rmsd_threshold", r.second.match_rmsd_threshold));
			result_object.push_back(Pair("match_rmsd", r.second.match_rmsd));
		}

		result_container.push_back(Pair( boost::lexical_cast<std::string>(r.first), Value(result_object) ));
	}
	os << "\n";
	utility::json_spirit::write(result_container, os, utility::json_spirit::pretty_print);
	os << std::endl;
}

core::Real FragmentLookupFilter::report_sm( core::pose::Pose const & pose ) const {
        //Get pose size just to avoid warning regarding pose
	pose.n_residue();
	using namespace core::indexed_structure_store;
	typedef std::map<core::Size, FragmentLookupResult> result_t;
	using utility::json_spirit::Value;
	using utility::json_spirit::Object;
	using utility::json_spirit::Pair;

	Object result_container;
	core::Size num_failed_matches=0;
	BOOST_FOREACH(result_t::value_type &r, cached_lookup_result_)
	{
		Object result_object;

		result_object.push_back(Pair("found_match", r.second.found_match));

		if(!r.second.found_match){
			num_failed_matches+=1;
		}

		if (r.second.found_match)
		{
			result_object.push_back(Pair("match_score", r.second.match_score));
			result_object.push_back(Pair("match_index", static_cast<boost::uint64_t>(r.second.match_index)));
			result_object.push_back(Pair("match_rmsd_threshold", r.second.match_rmsd_threshold));
			result_object.push_back(Pair("match_rmsd", r.second.match_rmsd));
		}

		result_container.push_back(Pair( boost::lexical_cast<std::string>(r.first), Value(result_object) ));
	}
	TR << "\n";
	utility::json_spirit::write(result_container, TR, utility::json_spirit::pretty_print);
	TR << std::endl;

	return ((core::Real)num_failed_matches);
}

bool FragmentLookupFilter::apply( Pose const & pose ) const
{
	TR << "Calculating Fragments Ideality: " ;
	core::Size const num_failed_matches( compute( pose ));
	//If # bad fragments are lower than the threshold
	if (num_failed_matches <= threshold_ ){
		return true;
	}
	//else by default
	return false;
}

void FragmentLookupFilter::parse_my_tag( utility::tag::TagCOP tag,
                           DataMap &,
                           Filters_map const &,
                           Movers_map const &,
                           Pose const &)
{
  std::string lookup_name = tag->getOption< std::string >( "lookup_name", "" );
  std::string store_path = tag->getOption< std::string >( "store_path", "");
  std::string lookup_mode = tag->getOption< std::string >( "lookup_mode", "");

  if (lookup_name == "")
  {
    utility_exit_with_message("FragmentLookupFilter tag without mandatory parameter: lookup_name");
  }



  if (store_path != "")
  {
    target_lookup_ = StructureStoreManager::get_instance()->load_fragment_lookup(lookup_name, store_path);
  }
  else
  {
    target_lookup_ = StructureStoreManager::get_instance()->load_fragment_lookup(lookup_name);
  }



  if (lookup_mode == "" || lookup_mode == "first")
  {
		lookup_mode_ = First;
  }
	else if(lookup_mode == "closest")
	{
		lookup_mode_ = Closest;
	}
	else
	{
    utility_exit_with_message("FragmentLookupFilter tag with invalid lookup_mode: " + lookup_mode);
	}

	if( tag->hasOption("chain") ){
		target_chain_ = tag->getOption<core::Size>( "chain" );
		b_target_chain_ = true;
	}
	if( tag->hasOption("threshold") ){
		threshold_ = tag->getOption<core::Size>( "threshold" );
	}else{
		TR.Warning << "No threshold option selected, using default = 0" << std::endl;
		threshold_ = 0;
	}
}

FragmentSpecification const & FragmentLookupFilter::fragment_specification()
{
	return target_lookup_->fragment_specification();
}

protocols::filters::FilterOP
FragmentLookupFilterCreator::create_filter() const { return new FragmentLookupFilter; }

std::string
FragmentLookupFilterCreator::keyname() const { return "FragmentLookupFilter"; }

}
}
}
