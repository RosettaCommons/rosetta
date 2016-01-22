// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/components/NamedMover.cc
/// @brief Named segment of residues
/// @details
/// @author Tom Linsky

//Unit Headers
#include <protocols/denovo_design/components/NamedMover.hh>

//Project Headers
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/util.hh>

//Protocols Headers

//Core Headers

//Basic/Utility/Numeric Headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

//Boost
#include <boost/algorithm/string.hpp>

//C++ Headers

static basic::Tracer TR("protocols.denovo_design.components.NamedMover");

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace denovo_design {
namespace components {

NamedMover::NamedMover() :
	protocols::moves::Mover(),
	id_( "" ),
	parent_id_( "" ),
	segment_names_()
{
	segment_names_.clear();
}

NamedMover::NamedMover( std::string const & idval, std::string const & parent_id ):
	protocols::moves::Mover(),
	id_( idval ),
	parent_id_( parent_id ),
	segment_names_()
{
	segment_names_.clear();
}

/// @brief setup the parameters via an xml tag
void
NamedMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	if ( tag->hasOption( "name" ) ) {
		set_id( tag->getOption< std::string >( "name" ) );
		if ( segment_names().empty() ) add_segment_name( id() );
	} else {
		throw utility::excn::EXCN_RosettaScriptsOption( "Name is required for " + tag->getName() );
	}
	if ( tag->hasOption( "parent" ) ) {
		set_parent_id( tag->getOption< std::string >( "parent" ) );
	}
}

/// @brief performs setup and applies loop building
/// @details Steps to apply():
/// 1. Pulls StructureData from the pose
/// 2. setup_permutation() stores data about this connection which is not
///    static for every apply() call
/// 3. apply_permutation() uses the StructureData object to build the loop
/// 4. check() checks the built loop
void
NamedMover::apply( core::pose::Pose & pose )
{
	if ( id().empty() ) {
		std::stringstream err;
		err << "A component of type " << get_name() << " has no name -- all NamedMover derivatives must have names." << std::endl;
		throw utility::excn::EXCN_Msg_Exception( err.str() );
	}

	// preserve header always should be on
	if ( !basic::options::option[basic::options::OptionKeys::run::preserve_header].value() ) {
		throw utility::excn::EXCN_BadInput( "To use NamedMover derived classes properly, the -run:preserve_header option must be specified." );
	}

	set_last_move_status( protocols::moves::MS_SUCCESS );

	//StructureDataOP sd = StructureData::create_from_pose( pose, id() );
	StructureDataOP sd( new MultiChainStructureData( id() ) );
	debug_assert( sd );

	// try 100 times to produce a valid permutation
	bool permutation_is_valid = false;
	for ( core::Size i = 1; i <= 100; ++i ) {
		try {
			setup_permutation( *sd );
		} catch ( EXCN_Setup const & e ) {
			TR.Error << id() << ": setup_permutation failed, error=";
			e.show( TR.Error );
			TR.Error << std::endl << "sd=" << *sd << std::endl;
			continue;
		}
		try {
			check_permutation( *sd );
		} catch ( EXCN_PreFilterFailed const & e ) {
			e.show( TR.Debug );
			continue;
		}
		permutation_is_valid = true;
		break;
	}
	if ( !permutation_is_valid ) {
		set_last_move_status( protocols::moves::FAIL_RETRY );
		TR.Error << "Failed to generate a valid permutation from the user input." << std::endl;
		return;
	}

	apply_permutation( *sd );
	protocols::moves::MoverStatus const applystatus = get_last_move_status();
	if ( applystatus != protocols::moves::MS_SUCCESS ) {
		TR.Error << id() << ": apply_permutation failed, sd=" << *sd << std::endl;
		set_last_move_status( applystatus );
		return;
	}

	if ( ! check( *sd ) ) {
		TR << id() << ": structure failed checks. Setting status to failure and not altering the pose." << std::endl;
		TR.Debug << *sd << std::endl;
		set_last_move_status( protocols::moves::FAIL_RETRY );
		return;
	}

	debug_assert( sd->pose() );
	pose.clear();
	pose = *(sd->pose());
}

/// @brief adds prefix if necessary, returns result
std::string
NamedMover::add_parent_prefix( std::string const & s ) const
{
	debug_assert( s != "" );
	std::string prefix = parent_id() + PARENT_DELIMETER;
	if ( !parent_id().empty() && !boost::starts_with(s, prefix) ) {
		return prefix + s;
	} else {
		return s;
	}
}

std::string const &
NamedMover::id() const
{
	return id_;
}

std::string const &
NamedMover::parent_id() const
{
	return parent_id_;
}

void
NamedMover::set_id( std::string const & idval )
{
	id_ = add_parent_prefix( idval );
}

void
NamedMover::set_parent_id( std::string const & parentid )
{
	parent_id_ = parentid;
	for ( StringList::iterator s=segment_names_.begin(); s!=segment_names_.end(); ++s ) {
		utility::vector1< std::string > const names = utility::string_split( *s, PARENT_DELIMETER );
		std::string const & basename = names[ names.size() ];
		*s = parent_id_ + PARENT_DELIMETER + basename;
	}
	if ( ! id_.empty() ) {
		set_id( id() );
	}
}

} // namespace components
} // namespace denovo_design
} // namespace protocols
