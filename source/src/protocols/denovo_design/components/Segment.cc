// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/denovo_design/components/Segment.cc
/// @brief Named segment of residues
/// @details
/// @author Tom Linsky

//Unit Headers
#include <protocols/denovo_design/components/Segment.hh>

//Protocol Headers
#include <protocols/denovo_design/util.hh>
#include <protocols/grafting/simple_movers/DeleteRegionMover.hh>

//Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/variant_util.hh>

//Basic/Utility/Numeric Headers
#include <basic/Tracer.hh>
#include <utility>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

//C++ Headers

static basic::Tracer TR("protocols.denovo_design.components.Segment");

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace denovo_design {
namespace components {

#ifdef    SERIALIZATION
Segment::Segment() :
	id_( "" ),
	posestart_( 1 ),
	movable_group_( 1 ),
	saferes_( 0 ),
	cutpoint_( 0 ),
	ss_( "" ),
	abego_( "" ),
	nterm_included_( false ),
	cterm_included_( false ),
	lower_segment_( "" ),
	upper_segment_( "" ),
	template_pose_(),
	lower_dihedrals_(),
	upper_dihedrals_(),
	lower_residue_(),
	upper_residue_()
{
}
#endif // SERIALIZATION

Segment::Segment( std::string const & id_val ) :
	id_( id_val ),
	posestart_( 1 ),
	movable_group_( 1 ),
	saferes_( 0 ),
	cutpoint_( 0 ),
	ss_( "" ),
	abego_( "" ),
	nterm_included_( false ),
	cterm_included_( false ),
	lower_segment_( "" ),
	upper_segment_( "" ),
	template_pose_(),
	lower_dihedrals_(),
	upper_dihedrals_(),
	lower_residue_(),
	upper_residue_()
{
}

Segment::Segment(
	std::string const & id_val,
	std::string const & ss_val,
	std::string const & abego_val,
	bool const start_inc,
	bool const stop_inc ):
	id_( id_val ),
	posestart_( 1 ),
	movable_group_( 1 ),
	saferes_( 0 ),
	cutpoint_( 0 ),
	ss_( ss_val ),
	abego_( abego_val ),
	nterm_included_( start_inc ),
	cterm_included_( stop_inc ),
	lower_segment_( "" ),
	upper_segment_( "" ),
	template_pose_(),
	lower_dihedrals_(),
	upper_dihedrals_(),
	lower_residue_(),
	upper_residue_()
{
	debug_assert( abego_.size() == ss_.size() );
}

SegmentOP
Segment::clone() const
{
	return SegmentOP( new Segment( *this ) );
}

std::string const &
Segment::id() const
{
	return id_;
}

void
Segment::set_id( std::string const & id_val )
{
	id_ = id_val;
}

void
Segment::clear()
{
	ss_ = "";
	set_abego( "" );
	template_pose_ = core::pose::PoseOP();
}

void
Segment::set_pose_start( core::Size const pose_resid )
{
	posestart_ = pose_resid;
}

core::Size
Segment::movable_group() const
{
	return movable_group_;
}

void
Segment::set_movable_group( core::Size const mg )
{
	movable_group_ = mg;
}

core::Size
Segment::lower_local() const
{
	if ( length() == 0 ) {
		return 0;
	} else {
		return 1;
	}
}

core::Size
Segment::upper_local() const
{
	return length();
}

core::Size
Segment::start_local() const
{
	if ( length() == 0 ) return 0;

	if ( nterm_included_ ) {
		return 1;
	} else {
		if ( length() > 1 ) {
			return 2;
		} else {
			return 0;
		}
	}
}

core::Size
Segment::stop_local() const
{
	if ( cterm_included_ ) {
		return length();
	} else {
		if ( length() > 0 ) {
			return length() - 1;
		} else {
			return 0;
		}
	}
}

core::Size
Segment::lower_padding() const
{
	return start_local() - lower_local();
}

core::Size
Segment::upper_padding() const
{
	return length() - stop_local();
}

core::Size
Segment::length() const
{
	return ss_.size();
}

core::Size
Segment::elem_length() const
{
	if ( stop_local() == 0 ) return 0;
	return stop_local() - start_local() + 1;
}

void
Segment::set_lower_segment( std::string const & lower_seg )
{
	lower_segment_ = lower_seg;
}

void
Segment::set_upper_segment( std::string const & upper_seg )
{
	upper_segment_ = upper_seg;
}

bool
Segment::contains( core::Size const pose_resid ) const
{
	return ( ( lower() <= pose_resid ) && ( pose_resid <= upper() ) );
}

/// @brief construct from xml tag
void
Segment::parse_tag( utility::tag::TagCOP tag )
{
	debug_assert( tag->getName() == "ResidueRange" );

	// check for required options
	if ( !tag->hasOption( "start" ) ) {
		TR << "start must be specified as an option to ResidueRange xml tag!!" << *tag << std::endl;
		debug_assert( tag->hasOption( "start" ) );
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "start must be specified as an option to ResidueRange xml tag!!" );
	}
	if ( !tag->hasOption( "ss" ) ) {
		TR << "ss must be specified as an option to ResidueRange xml tag!!" << *tag << std::endl;
		debug_assert( tag->hasOption( "ss" ) );
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "ss must be specified as an option to ResidueRange xml tag!!" );
	}
	if ( !tag->hasOption( "abego" ) ) {
		TR << "abego must be specified as an option to ResidueRange xml tag!!" << *tag << std::endl;
		debug_assert( tag->hasOption( "abego" ) );
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "abego must be specified as an option to ResidueRange xml tag!!" );
	}

	// sort out termini first
	nterm_included_ = tag->getOption< bool >( "nterm", nterm_included_ );
	cterm_included_ = tag->getOption< bool >( "cterm", cterm_included_ );

	posestart_ = tag->getOption< core::Size >( "start" );
	if ( !nterm_included_ ) {
		--posestart_;
	}

	ss_ = tag->getOption< std::string >( "ss" );
	set_abego( tag->getOption< std::string >( "abego" ) );
	if ( abego_.size() != length() ) {
		std::stringstream msg;
		msg << "Segment::parse_tag(): Invalid abego length ("
			<< abego_.size() << ") vs secondary structure length ("
			<< length() << ") when parsing " << *tag << std::endl;
		msg << "The abego length must match the secondary strcuture length." << std::endl;
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  msg.str() );
	}

	set_safe( tag->getOption< core::Size >( "safe", saferes_ ) );
	if ( saferes_ > elem_length() ) {
		std::stringstream msg;
		msg << "Segment::parse_tag(): Safe res (" << saferes_
			<< ") is larger than the length of the element (" << elem_length()
			<< ") encountered while parsing " << *tag << std::endl;
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  msg.str() );
	}

	set_cutpoint( tag->getOption< core::Size >( "cutpoint", cutpoint_ ) );
	if ( cutpoint_ > elem_length() ) {
		std::stringstream msg;
		msg << "Segment::parse_tag(): Cutpoint res (" << saferes_
			<< ") is larger than the length of the element (" << elem_length()
			<< ") encountered while parsing " << *tag << std::endl;
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  msg.str() );
	}

	set_lower_segment( tag->getOption< std::string >( "lower_segment", "" ) );
	set_upper_segment( tag->getOption< std::string >( "upper_segment", "" ) );

	movable_group_ = tag->getOption< core::Size >( "mgroup", movable_group_ );
	movable_group_ = tag->getOption< core::Size >( "group", movable_group_ );
}

void
Segment::parse_motif( std::string const & motif_str )
{
	clear();
	extend( "L", "X" );

	std::string secstruct = "";
	std::string abego = "";
	parse_motif_string( motif_str, secstruct, abego );
	extend( secstruct, abego );

	extend( "L", "X" );
}

/// @brief converts a pose resid to a segment resid
///        Numbering starts at start_res
///        start() --> N - posestart - start_local() + 1
///        start() + 1 --> N - posestart - start_local() + 2
SegmentResid
Segment::pose_to_segment( core::Size const pose_resid ) const
{
	core::Size const local_resid = pose_to_local( pose_resid );
	return local_to_segment( local_resid );
}

/// @brief converts a local resid to a pose resid
///        1 --> start_res - 1 + 1
///        2 --> start_res - 1 + 1
///        N --> start_res - 1 + N
core::Size
Segment::segment_to_pose( SegmentResid const segment_resid ) const
{
	core::Size const local_resid = segment_to_local( segment_resid );
	return local_to_pose( local_resid );
}

/// @brief converts an internal local resid to a segment resid
SegmentResid
Segment::local_to_segment( core::Size const local_resid ) const
{
	return local_resid - start_local() + 1;
}

/// @brief converts a pose resid to an internal local resid
core::Size
Segment::pose_to_local( core::Size const pose_resid ) const
{
	return pose_resid - posestart_ + 1;
}

/// @brief converts an internal local resid to a pose resid
///        1 --> posestart - 1 + 1
///        2 --> posestart - 1 + 2
///        N --> posestart - 1 + N
core::Size
Segment::local_to_pose( core::Size const local_resid ) const
{
	debug_assert( posestart_ );
	debug_assert( local_resid <= length() );
	return posestart_ - 1 + local_resid;
}

/// @brief converts a segment resid into a local resid
core::Size
Segment::segment_to_local( SegmentResid const segment_resid ) const
{
	debug_assert( segment_resid );
	if ( static_cast< core::Size >( std::abs( segment_resid ) ) > length() ) {
		std::stringstream msg;
		msg << "Segment::segment_to_local(): Given segment resid (" << segment_resid
			<< ") is larger than the elem_length of the segment (" << elem_length()
			<< "), segment=" << *this << std::endl;
		utility_exit_with_message( msg.str() );
	}
	if ( segment_resid < 0 ) {
		return stop_local() + segment_resid + 1;
	} else {
		return start_local() + segment_resid - 1;
	}
}

core::Size
Segment::start() const
{
	return local_to_pose( start_local() );
}

core::Size
Segment::stop() const
{
	return local_to_pose( stop_local() );
}

core::Size
Segment::safe() const
{
	if ( saferes_ ) return local_to_pose( saferes_ );

	core::Size const comp_safe = ( elem_length() + 1 )/ 2;
	if ( comp_safe ) {
		return segment_to_pose( comp_safe );
	} else {
		return 0; // should only happen in 0-elem-length segments
	}
}

SegmentResid
Segment::safe_segment() const
{
	if ( saferes_ == 0 ) return saferes_;

	return local_to_segment( saferes_ );
}

core::Size
Segment::cutpoint() const
{
	if ( cutpoint_ ) {
		return local_to_pose( cutpoint_ );
	} else {
		return cutpoint_; // 0
	}
}

SegmentResid
Segment::cutpoint_segment() const
{
	if ( cutpoint_ == 0 ) return cutpoint_;

	return local_to_segment( cutpoint_ );
}

/// @brief sets safe residue for this segment to be the ith residue in the segment
///        safe = segment_start - 1 + cut_res
void
Segment::set_safe( SegmentResid const segment_resid )
{
	if ( static_cast< core::Size >( std::abs( segment_resid ) ) > elem_length() ) {
		std::stringstream msg;
		msg << "Segment::set_safe(): The given safe residue value (" << segment_resid
			<< ") is larger than the length of the segment (" << elem_length()
			<< ")!" << std::endl;
		utility_exit_with_message( msg.str() );
	}
	if ( segment_resid == 0 ) {
		saferes_ = 0;
	} else {
		saferes_ = segment_to_local( segment_resid );
	}
}

/// @brief sets cutpoint for this segment to be the ith residue in the segment
///        cut = segment_start - 1 + cut_res
void
Segment::set_cutpoint( SegmentResid const segment_resid )
{
	if ( static_cast< core::Size >( std::abs( segment_resid ) ) > elem_length() ) {
		std::stringstream msg;
		msg << "Segment::set_cutpoint(): The given cutpoint value (" << segment_resid
			<< ") is larger than the length of the segment (" << elem_length()
			<< ")!" << std::endl;
		utility_exit_with_message( msg.str() );
	}
	if ( segment_resid == 0 ) {
		cutpoint_ = 0;
	} else {
		cutpoint_ = segment_to_local( segment_resid );
	}
}

/// @brief returns the n-terminal residue of this segment
core::Size
Segment::lower() const
{
	return local_to_pose( lower_local() );
}

/// @brief returns the n-terminal residue of this segment
core::Size
Segment::upper() const
{
	return local_to_pose( upper_local() );
}

core::Size
Segment::n_residues_before_cutpoint() const
{
	return cutpoint_;
}

core::Size
Segment::n_residues_after_cutpoint() const
{
	return length() - cutpoint_;
}

void
Segment::add_lower_padding()
{
	//debug_assert( !template_pose_ );
	if ( !nterm_included_ ) return;
	ss_ = 'L' + ss_;
	set_abego( 'X' + abego_ );
	if ( saferes_ ) saferes_ += 1;
	if ( cutpoint_ ) cutpoint_ += 1;
	nterm_included_ = false;
}

void
Segment::add_upper_padding()
{
	//debug_assert( !template_pose_ );
	if ( !cterm_included_ ) return;
	ss_ = ss_ + 'L';
	set_abego( abego_ + 'X' );
	cterm_included_ = false;
}

void
Segment::delete_lower_padding()
{
	// nothing to do
	if ( nterm_included_ ) {
		return;
	}

	core::Size const pad = lower_padding();
	ss_ = ss_.substr( pad, std::string::npos );
	set_abego( abego_.substr( pad, std::string::npos ) );
	if ( saferes_ ) {
		debug_assert( saferes_ > pad );
		saferes_ -= pad;
	}
	if ( cutpoint_ ) {
		debug_assert( cutpoint_ > pad );
		cutpoint_ -= pad;
	}

	//if ( template_pose_ ) {
	// core::kinematics::FoldTree ft;
	// ft.add_edge( template_pose_->size(), 1, -1 );
	// template_pose_->fold_tree( ft );
	// template_pose_->delete_polymer_residue( 1 );
	//}
	nterm_included_ = true;
}

void
Segment::delete_upper_padding()
{
	// nothing to do
	if ( cterm_included_ ) {
		return;
	}

	core::Size const stop_residue = stop_local();
	debug_assert( stop_residue > 0 );
	ss_ = ss_.substr( 0, stop_residue );
	set_abego( abego_.substr( 0, stop_residue ) );

	//if ( template_pose_ ) {
	// core::kinematics::FoldTree ft;
	// ft.add_edge( 1, template_pose_->size(), -1 );
	// template_pose_->fold_tree( ft );
	// template_pose_->delete_polymer_residue( template_pose_->size() );
	//}
	cterm_included_ = true;
}

/// @brief expands this residue set to include the dummy trailing residues
//void
//Segment::engulf_lower_padding()
//{
// nterm_included_ = true;
// // abego and ss stay the same
//}
//
///// @brief expands this residue set to include the dummy trailing residues
//void
//Segment::engulf_upper_padding()
//{
// cterm_included_ = true;
// // abego and ss stay the same
//}

core::Size
Segment::template_resid( SegmentResid const segment_resid ) const
{
	core::Size template_resid;
	if ( segment_resid < 0 ) {
		template_resid = elem_length() + 1 + segment_resid;
	} else {
		template_resid = static_cast< core::Size >( segment_resid );
	}
	return template_resid;
}

/// @brief given a segment residue number, delete that residue. Resid for start_local() == 1
void
Segment::delete_residue( SegmentResid const segment_resid )
{
	core::Size const local_resid = segment_to_local( segment_resid );
	if ( cutpoint_ && ( local_resid <= cutpoint_ ) ) {
		--cutpoint_;
	}
	if ( saferes_ && ( local_resid <= saferes_ ) ) {
		--saferes_;
	}

	// fix secondary structure
	std::stringstream newss;
	newss << ss_.substr( 0, local_resid - 1 ) << ss_.substr( local_resid, std::string::npos );
	ss_ = newss.str();

	// fix abego
	std::stringstream newabego;
	newabego << abego_.substr( 0, local_resid - 1 ) << abego_.substr( local_resid, std::string::npos );
	set_abego( newabego.str() );

	// delete residue in template pose
	if ( template_pose_ ) {
		core::Size const resid = template_resid( segment_resid );
		protocols::grafting::simple_movers::DeleteRegionMover del( resid, resid );
		del.apply( *template_pose_ );
	}
}

/// @brief given a residue number range local to this 1=start, length=end, delete the residue
void
Segment::delete_residues( core::Size const local_resnum_start, core::Size const local_resnum_stop )
{
	debug_assert( local_resnum_start >= 1 );
	debug_assert( local_resnum_stop <= stop_local() );
	debug_assert( local_resnum_start <= local_resnum_stop );

	core::Size const len = local_resnum_stop - local_resnum_start + 1;
	debug_assert( len < stop_local() );

	if ( local_resnum_stop <= cutpoint_ ) {
		cutpoint_ -= len;
		debug_assert( cutpoint_ >= start_local() );
	} else if ( local_resnum_start <= cutpoint_ ) {
		cutpoint_ = local_resnum_start;
	}
	if ( local_resnum_stop <= saferes_ ) {
		saferes_ -= len;
		debug_assert( saferes_ >= start_local() );
	} else if ( local_resnum_start <= saferes_ ) {
		saferes_ = local_resnum_start;
	}

	// fix secondary structure
	std::string newss = ss_.substr( 0, local_resnum_start-1 );
	newss += ss_.substr( local_resnum_stop, std::string::npos );
	ss_ = newss;

	// fix abego
	std::string newab = abego_.substr( 0, local_resnum_start-1 );
	newab += abego_.substr( local_resnum_stop, std::string::npos );
	set_abego( newab );
}

/// @brief Returns whether or not this segment has a template pose
bool
Segment::has_template_pose() const
{
	if ( template_pose_ ) {
		return true;
	} else {
		return false;
	}
}

/// @brief returns template pose
core::pose::PoseCOP
Segment::template_pose() const
{
	return template_pose_;
}

/// @brief sets template pose
void
Segment::set_template_pose(
	core::pose::Pose const & full_template_pose,
	core::Size const start_resid,
	core::Size const stop_resid )
{
	// this adds terminal variants that I may not want
	core::pose::PoseOP subpose( new core::pose::Pose( full_template_pose, start_resid, stop_resid ) );

	// if it does, remove them
	if ( !full_template_pose.residue( start_resid ).is_lower_terminus() ) {
		core::pose::remove_lower_terminus_type_from_pose_residue( *subpose, 1 );

		// reset locations of terminal hydrogen
		if ( subpose->residue(1).type().has( "H" ) ) {
			core::id::AtomID const id( subpose->residue(1).type().atom_index( "H" ), 1 );
			subpose->set_xyz( id, full_template_pose.residue( start_resid ).xyz( "H" ) );
			TR << "Set position of H from " << start_resid << std::endl;
		}
	}

	if ( !full_template_pose.residue( stop_resid ).is_upper_terminus() ) {
		core::pose::remove_upper_terminus_type_from_pose_residue( *subpose, subpose->size() );
		if ( subpose->residue( subpose->size() ).has( "O" ) ) {
			core::id::AtomID const id( subpose->residue(subpose->size()).type().atom_index( "O" ), subpose->size() );
			subpose->set_xyz( id, full_template_pose.residue( stop_resid ).xyz( "O" ) );
			TR << "Set position of O from " << stop_resid << std::endl;
		}
	}

	template_pose_ = subpose;
	if ( template_pose_->size() != elem_length() ) {
		std::stringstream msg;
		msg << "Segment::set_template_pose(): Template pose's length (" << template_pose_->size()
			<< ") does not match the segment length (" << elem_length() << "). Quitting." << std::endl;
		msg << " Segment definition: " << *this << std::endl;
		msg << "start=" << start_resid << " stop=" << stop_resid << " template_len=" <<
			full_template_pose.size() << " subpose_len=" << subpose->size() << std::endl;
		std::string const template_dump = "template_pose.pdb";
		std::string const subpose_dump = "template_subpose.pdb";
		msg << "Full template pose dumped to " << template_dump << ", extracted subpose dumped to"
			<< subpose_dump << std::endl;
		full_template_pose.dump_pdb( template_dump );
		subpose->dump_pdb( subpose_dump );
		utility_exit_with_message( msg.str() );
	}
	// TODO: The is-protein() requirement here is to work around a bug in carboyhydrate code. It should be eventually removed
	if ( ( start_resid - 1 > 0 ) &&
			( full_template_pose.residue( start_resid - 1 ).is_protein() ) &&
			( full_template_pose.residue( start_resid ).is_protein() ) ) {
		lower_dihedrals_ = ResidueDihedrals( full_template_pose, start_resid - 1 );
		lower_residue_ = full_template_pose.residue( start_resid - 1 ).clone();
	} else {
		lower_dihedrals_ = ResidueDihedrals();
	}
	// TODO: The is-protein() requirement here is to work around a bug in carboyhydrate code.  It should be eventually removed
	if ( ( stop_resid + 1 <= full_template_pose.size() ) &&
			( full_template_pose.residue( stop_resid + 1 ).is_protein() ) &&
			( full_template_pose.residue( stop_resid ).is_protein() ) ) {
		upper_dihedrals_ = ResidueDihedrals( full_template_pose, stop_resid );
		upper_residue_ = full_template_pose.residue( stop_resid + 1 ).clone();
	} else {
		upper_dihedrals_ = ResidueDihedrals();
	}
}

ResidueDihedrals const &
Segment::lower_dihedrals() const
{
	return lower_dihedrals_;
}

ResidueDihedrals const &
Segment::upper_dihedrals() const
{
	return upper_dihedrals_;
}

bool
Segment::has_lower_residue() const
{
	if ( lower_residue_ ) return true;
	else return false;
}

core::conformation::Residue const &
Segment::lower_residue() const
{
	return *lower_residue_;
}

bool
Segment::has_upper_residue() const
{
	if ( upper_residue_ ) return true;
	else return false;
}

core::conformation::Residue const &
Segment::upper_residue() const
{
	return *upper_residue_;
}

void
Segment::extend( std::string const & secstruct, std::string const & abego )
{
	ss_ += secstruct;
	abego_ += abego;
}

std::string const &
Segment::ss() const
{
	return ss_;
}

void
Segment::set_ss( SegmentResid const segment_resid, char const ss_type )
{
	ss_[ segment_to_local( segment_resid ) - 1 ] = ss_type;
}

std::string const &
Segment::abego() const
{
	return abego_;
}

char
Segment::abego( SegmentResid const segment_resid ) const
{
	return abego_[ segment_to_local( segment_resid ) - 1 ];
}

void
Segment::set_abego( std::string const & abego_str )
{
	if ( abego_str.size() != length() ) {
		std::stringstream msg;
		msg << "Segment::set_abego(): Abego string's length (" << abego_str.size()
			<< ") does not match the segment length (" << length() << "). Quitting." << std::endl;
		msg << " Segment definition: " << *this << std::endl;
		utility_exit_with_message( msg.str() );
	}
	abego_ = abego_str;
}

void
Segment::set_abego( utility::vector1< std::string > const & abego_vec )
{
	set_abego( abego_str( abego_vec ) );
}

/////////////// ResidueDihedrals Definitions ///////////////////

ResidueDihedrals::ResidueDihedrals():
	lower_phi_( 42.0 ),
	psi_( 42.0 ),
	omega_( 42.0 ),
	phi_( 42.0 ),
	upper_psi_( 42.0 ),
	upper_omega_( 42.0 )
{}

ResidueDihedrals::ResidueDihedrals( core::pose::Pose const & input, core::Size const lower_resid ):
	lower_phi_( input.phi( lower_resid ) ),
	psi_( input.psi( lower_resid ) ),
	omega_( input.omega( lower_resid ) ),
	phi_( input.phi( lower_resid + 1 ) ),
	upper_psi_( input.psi( lower_resid + 1 ) ),
	upper_omega_( input.omega( lower_resid + 1 ) )
{}

core::Real
ResidueDihedrals::phi() const
{
	return phi_;
}

core::Real
ResidueDihedrals::psi() const
{
	return psi_;
}

core::Real
ResidueDihedrals::omega() const
{
	return omega_;
}

core::Real
ResidueDihedrals::lower_phi() const
{
	return lower_phi_;
}

core::Real
ResidueDihedrals::upper_psi() const
{
	return upper_psi_;
}

core::Real
ResidueDihedrals::upper_omega() const
{
	return upper_omega_;
}

bool
residue_is_compatible( core::conformation::Residue const & rsd )
{
	if ( rsd.is_ligand() ) return false;
	if ( rsd.is_carbohydrate() ) return false;
	return true;
}

void
ResidueDihedrals::set_in_pose( core::pose::Pose & pose, core::Size const lower_resid ) const
{
	if ( residue_is_compatible( pose.residue( lower_resid ) ) ) {
		pose.set_psi( lower_resid, psi_ );
		pose.set_omega( lower_resid, omega_ );
		TR.Debug << "Set psi for " << lower_resid << " from " << pose.psi( lower_resid ) << " to " << psi_ << std::endl;
		TR.Debug << "Set omega for " << lower_resid << " from " << pose.omega( lower_resid ) << " to " << omega_ << std::endl;
	}
	if ( residue_is_compatible( pose.residue( lower_resid + 1 ) ) ) {
		pose.set_phi( lower_resid + 1, phi_ );
		pose.set_psi( lower_resid + 1, upper_psi_ );
		TR.Debug << "Set phi for " << lower_resid + 1 << " from " << pose.phi( lower_resid + 1 ) << " to " << phi_ << std::endl;
		TR.Debug << "Set psi for " << lower_resid + 1 << " from " << pose.psi( lower_resid + 1 ) << " to " << upper_psi_ << std::endl;
	}
}

/// output residueinfo
std::ostream &
operator<<( std::ostream & os, Segment const & res )
{
	os << "<ResidueRange name=\"" << res.id()
		<< "\" start=\"" << res.start()
		<< "\" stop=\"" << res.stop()
		<< "\" safe=\"" << res.safe_segment()
		<< "\" nterm=\"" << res.nterm_included_
		<< "\" cterm=\"" << res.cterm_included_
		<< "\" group=\"" << res.movable_group_
		<< "\" cutpoint=\"" << res.cutpoint_segment()
		<< "\" ss=\"" << res.ss()
		<< "\" abego=\"" << res.abego();

	if ( res.lower_segment() != "" ) {
		os << "\" lower_segment=\"" << res.lower_segment();
	}
	if ( res.upper_segment() != "" ) {
		os << "\" upper_segment=\"" << res.upper_segment();
	}

	os << "\" />";
	return os;
}

} // namespace components
} // namespace denovo_design
} // namespace protocols

#ifdef    SERIALIZATION

/// @brief Serialization method for ResidueDihedrals
template< class Archive >
void
protocols::denovo_design::components::ResidueDihedrals::save( Archive & arc ) const
{
	//arc( cereal::base_class< utility::pointer::ReferenceCount >( this ) );
	arc( CEREAL_NVP( lower_phi_ ) );
	arc( CEREAL_NVP( psi_ ) );
	arc( CEREAL_NVP( omega_ ) );
	arc( CEREAL_NVP( phi_ ) );
	arc( CEREAL_NVP( upper_psi_ ) );
	arc( CEREAL_NVP( upper_omega_ ) );
}

/// @brief Deserialization method for ResidueDihedrals
template< class Archive >
void
protocols::denovo_design::components::ResidueDihedrals::load( Archive & arc )
{
	//arc( cereal::base_class< utility::pointer::ReferenceCount >( this ) );
	arc( lower_phi_ );
	arc( psi_ );
	arc( omega_ );
	arc( phi_ );
	arc( upper_psi_ );
	arc( upper_omega_ );
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::denovo_design::components::ResidueDihedrals );
CEREAL_REGISTER_TYPE( protocols::denovo_design::components::ResidueDihedrals )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_denovo_design_components_ResidueDihedrals )

/// @brief Serialization method for Segment
template< class Archive >
void
protocols::denovo_design::components::Segment::save( Archive & arc ) const
{
	arc( CEREAL_NVP( id_ ) );
	arc( CEREAL_NVP( posestart_ ) );
	arc( CEREAL_NVP( movable_group_ ) );
	arc( CEREAL_NVP( saferes_ ) );
	arc( CEREAL_NVP( cutpoint_ ) );
	arc( CEREAL_NVP( ss_ ) );
	arc( CEREAL_NVP( abego_ ) );
	arc( CEREAL_NVP( nterm_included_ ) );
	arc( CEREAL_NVP( cterm_included_ ) );
	arc( CEREAL_NVP( lower_segment_ ) );
	arc( CEREAL_NVP( upper_segment_ ) );
	arc( CEREAL_NVP( template_pose_ ) );
	arc( CEREAL_NVP( lower_dihedrals_ ) );
	arc( CEREAL_NVP( upper_dihedrals_ ) );
	arc( CEREAL_NVP( lower_residue_ ) );
	arc( CEREAL_NVP( upper_residue_ ) );
}

/// @brief Deserialization method for Segment
template< class Archive >
void
protocols::denovo_design::components::Segment::load( Archive & arc )
{
	arc( id_ );
	arc( posestart_ );
	arc( movable_group_ );
	arc( saferes_ );
	arc( cutpoint_ );
	arc( ss_ );
	arc( abego_ );
	arc( nterm_included_ );
	arc( cterm_included_ );
	arc( lower_segment_ );
	arc( upper_segment_ );
	arc( template_pose_ );
	arc( lower_dihedrals_ );
	arc( upper_dihedrals_ );
	core::conformation::ResidueOP lower_res;
	arc( lower_res );
	lower_residue_ = lower_res;

	core::conformation::ResidueOP upper_res;
	arc( upper_res );
	upper_residue_ = upper_res;
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::denovo_design::components::Segment );
CEREAL_REGISTER_TYPE( protocols::denovo_design::components::Segment )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_denovo_design_components_Segment )

#endif // SERIALIZATION

