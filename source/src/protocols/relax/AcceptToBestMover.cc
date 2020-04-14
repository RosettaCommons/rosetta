// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/relax/AcceptToBestMover.cc
/// @brief Mover equivalent of the accept_to_best command
/// @author Jack Maguire, jackmaguire1444@gmail.com

#include <protocols/relax/AcceptToBestMover.hh>
#include <protocols/relax/AcceptToBestMoverCreator.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/util.hh>
#include <numeric/random/random.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/CacheableStringMap.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/DiagnosticData.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <core/scoring/ScoreFunction.hh>


#ifdef SERIALIZATION
// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/list.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/details/helpers.hpp>

#endif

static  basic::Tracer TR( "protocols.relax.AcceptToBestMover" );

namespace protocols {
namespace relax {

////////////////////////////////////////////////////////////////////////////////////////////////////

AcceptToBestMover::AcceptToBestMover() :
	protocols::moves::Mover( AcceptToBestMover::mover_name() )
{}

AcceptToBestMover::AcceptToBestMover( AcceptToBestMover const & ) = default;

AcceptToBestMover::~AcceptToBestMover() = default;

moves::MoverOP
AcceptToBestMover::fresh_instance() const {
	return utility::pointer::make_shared< AcceptToBestMover >();
}

moves::MoverOP
AcceptToBestMover::clone() const {
	return utility::pointer::make_shared< AcceptToBestMover >( *this );
}

std::string
AcceptToBestMover::get_name() const {
	return AcceptToBestMoverCreator::mover_name();
}

#ifdef SERIALIZATION
std::string
name_for_cacheable_data(){
	return "AcceptToBestMover_cached_data";
}

core::Real
get_cached_energy( core::pose::Pose const & pose ){
	auto const cacheable_data_type =
		core::pose::datacache::CacheableDataType::SCORE_MAP;

	if ( ! pose.data().has( cacheable_data_type ) ) {
		return 0;
	}

	basic::datacache::DiagnosticDataCOP data =
		utility::pointer::dynamic_pointer_cast< basic::datacache::DiagnosticData const > (
		pose.data().get_const_ptr( cacheable_data_type )
	);
	debug_assert( data.get() != nullptr );

	auto iter = data->data().find( name_for_cacheable_data() );
	if ( iter == data->data().end() ) return 0;

	core::Real const score = iter->second;
	return score;
}

void
set_cached_energy( core::pose::Pose & pose, core::Real const energy ){
	using namespace core::pose::datacache;
	using namespace basic::datacache;

	auto const cacheable_data_type =
		core::pose::datacache::CacheableDataType::SCORE_MAP;

	if ( ! pose.data().has( cacheable_data_type ) ) {
		std::map < std::string, double > dummy;
		pose.data().set(
			core::pose::datacache::CacheableDataType::SCORE_MAP,
			utility::pointer::make_shared< basic::datacache::DiagnosticData >( dummy )
		);
	}

	basic::datacache::DiagnosticDataOP data =
		utility::pointer::dynamic_pointer_cast< basic::datacache::DiagnosticData > (
		pose.data().get_ptr( cacheable_data_type )
	);
	debug_assert( data.get() != nullptr );

	data->data()[ name_for_cacheable_data() ] = energy;
}

void
detached_copy_from_string(
	std::string const & serialized_pose_to_copy,
	core::pose::Pose & pose //this will become the new copy
){
	core::pose::Pose deserialized_pose;
	std::istringstream iss( serialized_pose_to_copy );
	cereal::BinaryInputArchive arc( iss );
	arc( deserialized_pose );
	pose.detached_copy( deserialized_pose );
	//can we just do arc( pose ) here and forget about deserialized_pose?
}

std::string
get_cached_pose_str( core::pose::Pose const & pose ){
	using namespace core::pose::datacache;
	using namespace basic::datacache;

	auto const cacheable_data_type =
		core::pose::datacache::CacheableDataType::STRING_MAP; //SCORE_MAP

	if ( ! pose.data().has( cacheable_data_type ) ) {
		return "";
	}

	CacheableStringMapCOP data =
		utility::pointer::dynamic_pointer_cast< CacheableStringMap const > (
		pose.data().get_const_ptr( cacheable_data_type )
	);
	debug_assert( data.get() != nullptr );

	auto iter = data->map().find( name_for_cacheable_data() );
	if ( iter == data->map().end() ) return "";

	std::string const pose_string = iter->second;
	return pose_string;

}

std::string
serialize_pose( core::pose::Pose const & pose ){
	std::ostringstream oss;
	cereal::BinaryOutputArchive arc( oss );
	arc( pose );
	return oss.str();
}

void
set_cached_pose_str_with_self( core::pose::Pose & pose ){
	using namespace core::pose::datacache;
	using namespace basic::datacache;

	auto const cacheable_data_type =
		core::pose::datacache::CacheableDataType::STRING_MAP;

	if ( ! pose.data().has( cacheable_data_type ) ) {
		using basic::datacache::CacheableStringMap;
		pose.data().set(
			core::pose::datacache::CacheableDataType::STRING_MAP,
			utility::pointer::make_shared< CacheableStringMap >()
		);
	}

	CacheableStringMapOP data =
		utility::pointer::dynamic_pointer_cast< CacheableStringMap > (
		pose.data().get_ptr( cacheable_data_type )
	);
	debug_assert( data.get() != nullptr );

	data->map()[ name_for_cacheable_data() ] = serialize_pose( pose );
}

void
clear_cached_pose_str( core::pose::Pose & pose ){
	using namespace core::pose::datacache;
	using namespace basic::datacache;

	auto const cacheable_data_type =
		core::pose::datacache::CacheableDataType::STRING_MAP;

	if ( ! pose.data().has( cacheable_data_type ) ) {
		using basic::datacache::CacheableStringMap;
		pose.data().set(
			core::pose::datacache::CacheableDataType::STRING_MAP,
			utility::pointer::make_shared< CacheableStringMap >()
		);
	}

	CacheableStringMapOP data =
		utility::pointer::dynamic_pointer_cast< CacheableStringMap > (
		pose.data().get_ptr( cacheable_data_type )
	);
	debug_assert( data.get() != nullptr );

	data->map()[ name_for_cacheable_data() ] = serialize_pose( pose );
}


#endif

void
AcceptToBestMover::apply( core::pose::Pose & pose ) {

#ifdef SERIALIZATION

	runtime_assert( sfxn_ != nullptr );
	core::Real const score = (*sfxn_)(pose);

	std::string const & cached_pose_string = get_cached_pose_str( pose );
	core::Real const score_for_best_pose = get_cached_energy( pose );

	if ( cached_pose_string.empty() ) {
		//This is the first pose considered and thus the best, cache it
		set_cached_energy( pose, score );
		TR << "First Pose Observed" << std::endl;

	} else if ( score <= score_for_best_pose ) {
		//The new pose is "better", update cache
		set_cached_energy( pose, score );
		TR << "Change Accepted" << std::endl;

	} else {
		//The cached pose is "better", revert to it
		detached_copy_from_string( cached_pose_string, pose );
		TR << "Change Rejected, Reverting Back" << std::endl;
	}

	//No matter what, the current pose is the "best" one seen so far
	//Now, clear the cached pose (to prevent recursive caching) and cache this pose inside itself
	clear_cached_pose_str( pose );
	set_cached_pose_str_with_self( pose );

#else
	TR << pose.size() << std::endl;//Making sure variable doesn't go unused. Silly in this case, right?
	utility_exit_with_message( "AcceptToBestMover requires you to build with extras=serialization until someone smarter than Jack Maguire comes along to improve it." );
#endif
}


void
AcceptToBestMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
) {
	sfxn_ = rosetta_scripts::parse_score_function( tag, data );
	runtime_assert( sfxn_ != nullptr );
}


std::string AcceptToBestMover::mover_name() {
	return "AcceptToBestMover";
}

void AcceptToBestMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	rosetta_scripts::attributes_for_parse_score_function( attlist );
	moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Mover equivalent of the accept_to_best command", attlist );
}

std::string AcceptToBestMoverCreator::mover_name() {
	return "AcceptToBestMover";
}

std::string AcceptToBestMoverCreator::keyname() const {
	return AcceptToBestMover::mover_name();
}

moves::MoverOP
AcceptToBestMoverCreator::create_mover() const {
	return utility::pointer::make_shared< AcceptToBestMover >();
}

void
AcceptToBestMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AcceptToBestMover::provide_xml_schema( xsd );
}

} // namespace relax
} // namespace protocols
