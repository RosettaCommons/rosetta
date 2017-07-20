// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/InnerLarvalJob.cc
/// @brief  Method definitions for the InnerLarvalJob class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/jd3/InnerLarvalJob.hh>

// Project headers
#include <protocols/jd3/InputSource.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.hh>

// Basic headers
#include <basic/datacache/ConstDataMap.hh>
#include <basic/resource_manager/JobOptions.hh>

//C++ headers
#include <string>
#include <sstream>

// Utility headers
#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {


InnerLarvalJob::InnerLarvalJob() :
	input_sources_(),
	nstruct_( 1 ),
	bad_( false ),
	const_data_cache_( new basic::datacache::ConstDataMap )
{}

InnerLarvalJob::InnerLarvalJob( InnerLarvalJob const & src ) :
	input_tag_                ( src.input_tag_                ),
	job_tag_                  ( src.job_tag_                  ),
	outputter_                ( src.outputter_                ),
	input_sources_            ( src.input_sources_            ),
	nstruct_                  ( src.nstruct_                  ),
	input_job_result_indices_ ( src.input_job_result_indices_ ),
	bad_                      ( src.bad_                      ),
	const_data_cache_         ( src.const_data_cache_         ),
	jobdef_tag_               ( src.jobdef_tag_               )
{}


InnerLarvalJob::InnerLarvalJob( core::Size nstruct ) :
	input_sources_(),
	nstruct_( nstruct ),
	bad_( false ),
	const_data_cache_( new basic::datacache::ConstDataMap )
{}

InnerLarvalJob::~InnerLarvalJob() = default;

/// @brief Return true if either pointers point at the same object or
/// if they point at objects that are themselves equal.  Requires that
/// class T implements operator ==
template < class T >
bool
pointer_equals(
	utility::pointer::shared_ptr< T > const & first,
	utility::pointer::shared_ptr< T > const & second
)
{
	return first == second || ( first != 0 && second != 0 && *first == *second );
}

bool InnerLarvalJob::sources_same( InnerLarvalJob const & other ) const
{
	if ( input_sources_.size() != other.input_sources_.size() ) return false;
	for ( core::Size ii = 1; ii <= input_sources_.size(); ++ii ) {
		if ( *input_sources_[ii] != *other.input_sources_[ii] ) return false;
	}
	return true;
}

/// @brief Mutual comparison of this inner job to the other inner job
/// so that if either one thinks it's not the same as the other, then
/// it returns false.  Invokes the same() function on both this and other
bool
InnerLarvalJob::operator == ( InnerLarvalJob const & other ) const
{
	return input_tag_ == other.input_tag_ &&
		job_tag_ == other.job_tag_ &&
		outputter_ == other.outputter_ &&
		nstruct_ == other.nstruct_ &&
		// bad_ == other.bad_ ? No. Two jobs can be equal even if one has been marked bad and the other not yet marked
		sources_same( other ) &&
		*const_data_cache_ == *other.const_data_cache_ &&
		// pointer_equals( job_options_, other.job_options_ ) &&
		same_type( other ) && // make sure the dynamic type of other is equivalent to that of this
		other.same_type( *this ); // make sure the dynamic type of this is equivalent to other
}

bool
InnerLarvalJob::operator != (InnerLarvalJob const & other ) const
{
	return ! (*this == other);
}

/// @details Since this is the base-class function, then by construction
/// other is equivalent to this.
/// @note classes derived from InnerLarvalJob must perform dynamic casts
/// to ensure the other InnerLarvalJob has the same type as them
bool
InnerLarvalJob::same_type( InnerLarvalJob const & ) const
{
	return true;
}

void
InnerLarvalJob::show( std::ostream & out ) const
{
	out << "InnerLarvalJob::show stubbed out";
}

std::ostream &
operator<< ( std::ostream & out, const InnerLarvalJob & inner_job )
{
	inner_job.show( out );
	return out;
}

void InnerLarvalJob::input_source( InputSourceCOP setting )
{
	input_sources_.clear();
	input_sources_.push_back( setting );
}

void
InnerLarvalJob::clear_input_sources()
{
	input_sources_.clear();
}

void
InnerLarvalJob::append_input_source( InputSourceCOP setting )
{
	input_sources_.push_back( setting );
}

std::string
InnerLarvalJob::input_tag() const
{
	if ( input_tag_ == "" ) {
		return concatenated_input_source_names();
	} else {
		return input_tag_;
	}
}

void InnerLarvalJob::input_tag( std::string const & setting )
{
	input_tag_ = setting;
}


std::string InnerLarvalJob::job_tag() const
{
	if ( job_tag_ == "" ) { return input_tag(); }
	else { return job_tag_; }
}


void InnerLarvalJob::job_tag( std::string const & setting )
{
	job_tag_ = setting;
}

std::string InnerLarvalJob::outputter() const
{
	return outputter_;
}

void InnerLarvalJob::outputter( std::string const & setting )
{
	outputter_ = setting;
}

core::Size InnerLarvalJob::nstruct_max() const
{
	return nstruct_;
}

basic::datacache::ConstDataMap const &
InnerLarvalJob::const_data_map() const
{
	return *const_data_cache_;
}

utility::tag::TagCOP
InnerLarvalJob::jobdef_tag() const
{
	return jobdef_tag_;
}

void
InnerLarvalJob::jobdef_tag( utility::tag::TagCOP setting )
{
	jobdef_tag_ = setting;
}

// basic::resource_manager::JobOptions const &
// InnerLarvalJob::job_options() const
// {
//  return *job_options_;
// }

void
InnerLarvalJob::nstruct_max( core::Size setting )
{
	nstruct_ = setting;
}

core::Size InnerLarvalJob::n_input_sources() const
{
	return input_sources_.size();
}

/// @brief Read access to the InputSource object that describes how the Pose for
/// this job should be constructed.
InputSource const &
InnerLarvalJob::input_source() const
{
	debug_assert( input_sources_.size() == 1 );
	return *input_sources_[1];
}

InputSource const &
InnerLarvalJob::input_source( Size index ) const
{
	return *input_sources_[ index ];
}

utility::vector1< JobOutputID > const &
InnerLarvalJob::input_job_result_indices() const
{
	return input_job_result_indices_;
}

void InnerLarvalJob::add_input_job_result_index( JobOutputID job_index )
{
	input_job_result_indices_.push_back( job_index );
}


/// @brief Store the fact that this job is bad, for some reason or another, e.g. it has
/// malformed inputs.
void InnerLarvalJob::set_bad( bool setting )
{
	bad_ = setting;
}

/// @brief Has this job been labeled as bad?
bool InnerLarvalJob::bad() const
{
	return bad_;
}

/// @brief Write access to the ConstDataMap that the %InnerLarvalJob holds.  The objects pointed
/// to by the const data map can be shared between two %InnerLarvalJobs, but the data maps
/// themselves should not be shared.
basic::datacache::ConstDataMap &
InnerLarvalJob::const_data_map()
{
	return *const_data_cache_;
}

// @brief Set the JobOptions object held by the %InnerLarvalJob -- two %InnerLarvalJobs may point to
// the same JobOptions object.
//void InnerLarvalJob::job_options( basic::resource_manager::JobOptionsCOP setting )
//{
// job_options_ = setting;
//}

/// @brief construct a string from the input_sources_
std::string
InnerLarvalJob::concatenated_input_source_names() const
{
	std::ostringstream oss;
	if ( input_sources_.size() ) {
		for ( core::Size ii = 1; ii < input_sources_.size(); ++ii ) {
			oss << input_sources_[ii]->input_tag() << "_";
		}
		oss << input_sources_[ input_sources_.size() ]->input_tag();
		return oss.str();
	} else {
		return "unspecified";
	}
}

} // namespace jd3
} // namespace protocols

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::jd3::InnerLarvalJob::save( Archive & arc ) const {
	arc( CEREAL_NVP( input_tag_ ) ); // std::string
	arc( CEREAL_NVP( job_tag_ ) ); // std::string
	arc( CEREAL_NVP( outputter_ ) ); // std::string
	arc( CEREAL_NVP( input_sources_ ) ); // utility::vector1<InputSourceCOP>
	arc( CEREAL_NVP( nstruct_ ) ); // core::Size
	arc( CEREAL_NVP( input_job_result_indices_ ) ); // utility::vector1<core::Size>
	arc( CEREAL_NVP( bad_ ) ); // _Bool
	arc( CEREAL_NVP( const_data_cache_ ) ); // basic::datacache::ConstDataMapOP
	arc( CEREAL_NVP( jobdef_tag_ ) ); // utility::tag::TagCOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::jd3::InnerLarvalJob::load( Archive & arc ) {
	arc( input_tag_ ); // std::string
	arc( job_tag_ ); // std::string
	arc( outputter_ ); // std::string
	utility::vector1< std::shared_ptr< protocols::jd3::InputSource > > local_input_sources;
	arc( local_input_sources ); // utility::vector1<InputSourceCOP>
	input_sources_ = local_input_sources; // copy the non-const pointer(s) into the const pointer(s)
	arc( nstruct_ ); // core::Size
	arc( input_job_result_indices_ ); // utility::vector1<core::Size>
	arc( bad_ ); // _Bool
	arc( const_data_cache_ ); // basic::datacache::ConstDataMapOP
	utility::tag::TagOP local_jobdef_tag;
	arc( local_jobdef_tag ); // utility::tag::TagOP
	jobdef_tag_ = local_jobdef_tag;
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::InnerLarvalJob );
CEREAL_REGISTER_TYPE( protocols::jd3::InnerLarvalJob )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_jd3_InnerLarvalJob )
#endif // SERIALIZATION
