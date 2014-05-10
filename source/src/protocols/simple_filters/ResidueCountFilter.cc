// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/ResidueCountFilter.cc
/// @brief Filter on the total number of residues in the structure
/// @author Matthew O'Meara (mattjomeara@gmail.com)
/// @author Tom Linsky (tlinsky@gmail.com) [residue type option]
/// @author Doo	Nam	Kim (doonam.kim@gmail.com) [enables	count/filter by percentage, not by raw counted number]

//Unit Headers
#include <protocols/simple_filters/ResidueCountFilter.hh>
#include <protocols/simple_filters/ResidueCountFilterCreator.hh>
#include <utility/tag/Tag.hh>
#include <core/pose/Pose.hh>
#include <protocols/filters/Filter.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
//Project Headers
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <basic/Tracer.hh>
//Utility Headers
#include <utility/string_util.hh>
#include <protocols/rosetta_scripts/util.hh>

namespace protocols{
namespace simple_filters {

using namespace core;
using namespace core::scoring;

static basic::Tracer TR( "protocols.simple_filters.ResidueCountFilter" );

protocols::filters::FilterOP
ResidueCountFilterCreator::create_filter() const { return new ResidueCountFilter; }

std::string
ResidueCountFilterCreator::keyname() const { return "ResidueCount"; }

//default ctor
ResidueCountFilter::ResidueCountFilter() :
	protocols::filters::Filter( "ResidueCount" ),
	max_residue_count_(0),
	enable_max_residue_count_(false),
	min_residue_count_(0),
	enable_min_residue_count_(false),
	count_as_percentage_(false),	// for a user who does not use rosetta_scripts, count_as_percentage_ = false by default here
  packable_( 0 ),
 	task_factory_( NULL )
{}

ResidueCountFilter::ResidueCountFilter(
	ResidueCountFilter const & src
) :
	protocols::filters::Filter( "ResidueCount" ),
	max_residue_count_(src.max_residue_count_),
	enable_max_residue_count_(src.enable_max_residue_count_),
	min_residue_count_(src.min_residue_count_),
	enable_min_residue_count_(src.enable_min_residue_count_),
	count_as_percentage_(src.count_as_percentage_),	
	res_types_( src.res_types_ ),
  packable_( src.packable_ ),
  task_factory_( src.task_factory_ )
{}


ResidueCountFilter::~ResidueCountFilter() {}


filters::FilterOP
ResidueCountFilter::clone() const {
	return new ResidueCountFilter( *this );
}

filters::FilterOP
ResidueCountFilter::fresh_instance() const {
	return new ResidueCountFilter();
}


Real
ResidueCountFilter::round_to_Real(
	Real x)	const
{
	Real rounded = floor((x * 10) +	0.5)	/	10;
	return rounded;
} //round_to_Real

void
ResidueCountFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	filters::Filters_map const &,
	moves::Movers_map const &,
	core::pose::Pose const & pose
) {
	if(tag->hasOption("max_residue_count")){
		enable_max_residue_count(true);
		max_residue_count(tag->getOption< core::Size >("max_residue_count"));
	}

	if(tag->hasOption("min_residue_count")){
		enable_min_residue_count(true);
		min_residue_count(tag->getOption< core::Size >("min_residue_count"));
	}

	count_as_percentage_ = tag->getOption<bool>("count_as_percentage", false);
		//	definition: if true, count residues as percentage (=100*raw_number_of_specified_residue/total_residue)	instead of counting raw number of it
		//				max_residue_count/min_residue_count	are also	assumed to be entered as percentage

	if(tag->hasOption("residue_types")){
		std::string const res_type_str( tag->getOption< std::string >("residue_types", "") );
		utility::vector1< std::string > const res_type_vec( utility::string_split( res_type_str, ',' ) );
		TR << "Residue types specified: " << res_type_vec << std::endl;
		// get the residue type set from the first residue of the input pose
		runtime_assert( pose.total_residue() >= 1 );
		core::chemical::ResidueTypeSet const & res_type_set( pose.residue( 1 ).residue_type_set() );
		// loop through residue types, check to see if they are valid, and add them if they are valid
		for ( core::Size i=1; i<=res_type_vec.size(); ++i ) {
			if ( ! add_residue_type_by_name( res_type_set, res_type_vec[i] ) ) {
				// try to add the residue type -- if it doesn't exist in the residue type set, inform the user and exit
				utility_exit_with_message( "An invalid residue type (" + res_type_vec[i] + ") was specified to the ResidueCount filter." );
			}
		} // for all res type specified
	}
  task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
  packable( tag->getOption< bool >( "packable", 0 ) );
}

bool
ResidueCountFilter::apply(
	core::pose::Pose const & pose
) const {
	if(enable_max_residue_count() && compute(pose) > max_residue_count()){
		return false;
	}

	if(enable_min_residue_count() && compute(pose) < min_residue_count()){
		return false;
	}

	return true;
}

void
ResidueCountFilter::report(
	std::ostream & out,
	core::pose::Pose const & pose
) const {
	out << "Residue Count: " << compute( pose ) <<std::endl;

	//commented out below since it breaks an existing unit test (Doonam Kim)
//	if	(count_as_percentage_)
//	{
//		out << "Residue Count as percentage: " << compute( pose ) <<std::endl;
//	}
//	else //	(!count_as_percentage_)
//	{
//		out << "Residue Count as raw numer: " << compute( pose ) <<std::endl;
//	}
}

core::Real
ResidueCountFilter::report_sm(
	core::pose::Pose const & pose
) const {
	return compute( pose );
}

core::Real
ResidueCountFilter::compute(
	core::pose::Pose const & pose
) const {

	utility::vector1< core::Size > target_res;
	target_res.clear();
	if( !task_factory() ){
		for( core::Size i = 1; i <= pose.total_residue(); ++i )
			target_res.push_back( i );
	} else {
		if( packable_ ) {
			target_res = protocols::rosetta_scripts::residue_packer_states( pose, task_factory(), true/*designable*/, true/*packable*/ );
		} else {
			target_res = protocols::rosetta_scripts::residue_packer_states( pose, task_factory(), true/*designable*/, false/*packable*/ );
		} 
	}

  core::Size count( 0 );
  for ( core::Size i=1; i<=target_res.size(); ++i ) {
  	if( res_types_.size() > 0 ) {
  		for ( core::Size j=1; j<=res_types_.size(); ++j ) {
  			if ( pose.residue( target_res[i] ).name3() == res_types_[j] || pose.residue( target_res[i] ).name() == res_types_[j] ) {
  				++count;
  			} // if
  		} // for all res types specified
  	} else {
  		++count; // for all packable/designable residues if task specified or all residues if no task specified
  	}
	}
	if	(count_as_percentage_) {
		Real	count_as_percentage	=	round_to_Real((static_cast<Real>(count)*100)/static_cast<Real>(target_res.size()));
			//without	static_cast<Real>,	division returns only Size (integer-like)
		return count_as_percentage;
	} else {
		return count;
	}
}

core::Size
ResidueCountFilter::max_residue_count() const {
	return max_residue_count_;
}

void
ResidueCountFilter::max_residue_count(
	core::Size value
) {
	max_residue_count_ = value;
}

utility::vector1< std::string >
ResidueCountFilter::res_types() const {
	return res_types_;
}

void
ResidueCountFilter::res_types(
  utility::vector1< std::string > const & res_types
) {
	res_types_.clear();
	for ( core::Size i=1; i<=res_types.size(); ++i ) {
		res_types_.push_back( res_types[i] );
	}
}

bool
ResidueCountFilter::enable_max_residue_count() const {
	return enable_max_residue_count_;
}

void
ResidueCountFilter::enable_max_residue_count(
	bool value
) {
	enable_max_residue_count_ = value;
}

core::Size
ResidueCountFilter::min_residue_count() const {
	return min_residue_count_;
}

void
ResidueCountFilter::min_residue_count(
	core::Size value
) {
	min_residue_count_ = value;
}

bool
ResidueCountFilter::enable_min_residue_count() const {
	return enable_min_residue_count_;
}

void
ResidueCountFilter::enable_min_residue_count(
	bool value
) {
	enable_min_residue_count_ = value;
}

core::pack::task::TaskFactoryOP
ResidueCountFilter::task_factory() const
{
	return task_factory_;
}

void
ResidueCountFilter::task_factory( core::pack::task::TaskFactoryOP task_factory )
{
	task_factory_ = task_factory;
}

bool
ResidueCountFilter::packable() const
{
	return( packable_ );
}

void
ResidueCountFilter::packable( bool const pack )
{
	packable_ = pack;
}


///@brief Checks whether a residue type is present in the provided residue type set, and if so, adds it to res_types_
///@detail given user specified residue type string, look for residue type name and add it to the list of types to count if it is present in the specified ResidueTypeSet.
///@input res_type_set, the residue type set of the input structure
///@input res_type_input, the user specified residue type name
///@return false if res_type_input doesn't match any residue type names, true otherwise
bool
ResidueCountFilter::add_residue_type_by_name(
  core::chemical::ResidueTypeSet const & res_type_set,
  std::string const & res_type_input
) {
	using namespace core::chemical;

	// check for the name in the residue type set
	if ( !res_type_set.has_name( res_type_input ) ) {
		return false;
	}

	res_types_.push_back( res_type_input );
	return true;
}

} // namespace
} // namespace
