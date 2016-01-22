// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/components/Segment.cc
/// @brief Named segment of residues
/// @details
/// @author Tom Linsky

//Unit Headers
#include <protocols/denovo_design/components/Segment.hh>

//Project Headers
#include <protocols/denovo_design/util.hh>

//Protocols Headers

//Core Headers
#include <core/pose/Pose.hh>

//Basic/Utility/Numeric Headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

//C++ Headers

static basic::Tracer TR("protocols.denovo_design.components.Segment");

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace denovo_design {
namespace components {

Segment::Segment() :
	movable_group( 0 ),
	is_loop( false ),
	posestart_( 0 ),
	start_( 0 ),
	stop_( 0 ),
	saferes_( 0 ),
	cutpoint_( 0 ),
	ss_( "" ),
	nterm_included_( false ),
	cterm_included_( false ),
	lower_segment_( "" ),
	upper_segment_( "" )
{
	abego_.clear();
}

Segment::Segment(
	core::Size const pose_length_val,
	core::Size const local_saferes,
	core::Size const local_cutpoint,
	core::Size const movable_group_val,
	bool const is_loop_val,
	bool const start_inc,
	bool const stop_inc,
	std::string const & lower,
	std::string const & upper,
	std::string const & ss_val,
	utility::vector1< std::string > const & abego_val ) :
	movable_group( movable_group_val ),
	is_loop( is_loop_val ),
	posestart_( 1 ),
	start_( 0 ),
	stop_( 0 ),
	saferes_( local_saferes ),
	cutpoint_( local_cutpoint ),
	ss_( ss_val ),
	abego_( abego_val ),
	nterm_included_( start_inc ),
	cterm_included_( stop_inc ),
	lower_segment_( lower ),
	upper_segment_( upper )
{
	debug_assert( pose_length_val );
	debug_assert( local_cutpoint <= pose_length_val );
	if ( start_inc ) {
		start_ = 1;
	} else {
		debug_assert( pose_length_val > 1 );
		start_ = 2;
	}
	if ( stop_inc ) {
		stop_ = pose_length_val;
	} else {
		debug_assert( pose_length_val > 1 );
		stop_ = pose_length_val - 1;
	}
	debug_assert( local_saferes >= start_ );
	debug_assert( local_saferes <= stop_ );
	debug_assert( ss_.size() == length() );
	debug_assert( ss_.size() == abego_.size() );
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
		throw utility::excn::EXCN_RosettaScriptsOption( "start must be specified as an option to ResidueRange xml tag!!" );
	}
	if ( !tag->hasOption( "stop" ) ) {
		TR << "stop must be specified as an option to ResidueRange xml tag!!" << *tag << std::endl;
		debug_assert( tag->hasOption( "stop" ) );
		throw utility::excn::EXCN_RosettaScriptsOption( "stop must be specified as an option to ResidueRange xml tag!!" );
	}
	if ( !tag->hasOption( "ss" ) ) {
		TR << "ss must be specified as an option to ResidueRange xml tag!!" << *tag << std::endl;
		debug_assert( tag->hasOption( "ss" ) );
		throw utility::excn::EXCN_RosettaScriptsOption( "ss must be specified as an option to ResidueRange xml tag!!" );
	}
	if ( !tag->hasOption( "abego" ) ) {
		TR << "abego must be specified as an option to ResidueRange xml tag!!" << *tag << std::endl;
		debug_assert( tag->hasOption( "abego" ) );
		throw utility::excn::EXCN_RosettaScriptsOption( "abego must be specified as an option to ResidueRange xml tag!!" );
	}

	// sort out termini first
	nterm_included_ = tag->getOption< bool >( "nterm", nterm_included_ );
	cterm_included_ = tag->getOption< bool >( "cterm", cterm_included_ );

	posestart_ = tag->getOption< core::Size >( "start" );
	if ( !nterm_included_ ) {
		--posestart_;
	}

	if ( nterm_included_ ) {
		start_ = 1;
	} else {
		start_ = 2;
	}

	stop_ = tag->getOption< core::Size >( "stop" );
	debug_assert( stop_ >= posestart_ );
	stop_ = stop_ - posestart_ + 1;

	if ( tag->hasOption( "safe" ) ) {
		saferes_ = tag->getOption< core::Size >( "safe" );
		debug_assert( saferes_ >= posestart_ );
		saferes_ = saferes_ - posestart_ + 1;
	} else {
		saferes_ = ( stop_ + start_ )/ 2;
	}

	if ( ( saferes_ > stop_ ) || ( saferes_ < start_ ) ) {
		TR.Error << "Invalid safe residue in " << *tag << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( "Invalid safe residue" );
	}

	movable_group = tag->getOption< core::Size >( "mgroup", movable_group );
	is_loop = tag->getOption< bool >( "loop", is_loop );
	cutpoint_ = tag->getOption< core::Size >( "cutpoint", cutpoint_ );
	lower_segment_ = tag->getOption< std::string >( "lower_segment", "" );
	upper_segment_ = tag->getOption< std::string >( "upper_segment", "" );
	ss_ = tag->getOption< std::string >( "ss" );
	if ( ss_.size() != length() ) {
		TR.Error << "Invalid ss length " << ss_.size() << " vs segment length " << length() << " in " << *tag << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( "Invalid ss length specified to ResidueSegment" );
	}

	std::string const abego_str = tag->getOption< std::string >( "abego" );
	abego_.clear();
	abego_ = abego_vector( abego_str );
	if ( abego_.size() != length() ) {
		TR.Error << "Invalid abego length " << abego_.size() << " vs segment length " << length() << " in " << *tag << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( "Invalid abego length specified to ResidueSegment" );
	}
}

core::Size
Segment::resid( core::Size const local_resnum ) const
{
	debug_assert( local_resnum );
	debug_assert( local_resnum <= stop_ );
	return start() + local_resnum - 1;
}

core::Size
Segment::start() const
{
	return posestart_ + start_ - 1;
}

core::Size
Segment::stop() const
{
	return posestart_ + stop_ - 1;
}

core::Size
Segment::safe() const
{
	return posestart_ + saferes_ - 1;
}

core::Size
Segment::cutpoint() const
{
	if ( cutpoint_ ) {
		return posestart_ + cutpoint_ - 1;
	} else {
		return 0;
	}
}

/// @brief returns the n-terminal residue of this segment
core::Size
Segment::nterm_resi() const
{
	if ( nterm_included_ ) {
		return start();
	} else {
		debug_assert( start() > 1 );
		return start() - 1;
	}
}

core::Size
Segment::cterm_resi() const
{
	if ( cterm_included_ ) {
		return stop();
	} else {
		return stop() + 1;
	}
}

std::string
Segment::serialize() const
{
	std::string ab = "";
	for ( core::Size i=1; i<=abego_.size(); ++i ) {
		ab += abego_[i][0];
	}
	return boost::lexical_cast< std::string >( start() ) + ":" +
		boost::lexical_cast< std::string >( stop() ) + ":" +
		boost::lexical_cast< std::string >( safe() ) + ":" +
		boost::lexical_cast< std::string >( movable_group ) + ":" +
		boost::lexical_cast< std::string >( is_loop ) + ":" +
		boost::lexical_cast< std::string >( nterm_included_ ) + ":" +
		boost::lexical_cast< std::string >( cterm_included_ ) + ":" +
		lower_segment_ + ":" + upper_segment_ + ":" + ss_ + ":" + ab;
}

void
Segment::set_nterm_included( bool const ntermval )
{
	if ( nterm_included_ && ( ntermval == false ) ) {
		ss_ = 'L' + ss_;
		utility::vector1< std::string > newabego;
		newabego.push_back( "X" );
		for ( int i=1, endi=abego_.size(); i<=endi; ++i ) {
			newabego.push_back( abego_[i] );
		}
		debug_assert( newabego.size() == abego_.size() + 1 );
		abego_ = newabego;
		start_ += 1;
		saferes_ += 1;
		stop_ += 1;
		if ( cutpoint_ ) {
			cutpoint_ += 1;
		}
	}
	nterm_included_ = ntermval;
}

void
Segment::set_cterm_included( bool const ctermval )
{
	if ( cterm_included_ && ( ctermval == false ) ) {
		ss_ = ss_ + 'L';
		abego_.push_back( "X" );
	}
	cterm_included_ = ctermval;
}

void
Segment::delete_leading_residues()
{
	// nothing to do
	if ( nterm_included_ ) {
		return;
	}

	core::Size const pad = nterm_pad();
	std::string const newss = ss_.substr(pad,std::string::npos);
	ss_ = newss;
	utility::vector1< std::string > newabego( abego_.size() - pad );
	std::copy( abego_.begin() + pad, abego_.end(), newabego.begin() );
	abego_ = newabego;
	nterm_included_ = true;
	start_ = 1;
	debug_assert( stop_ > pad );
	stop_ -= pad;
	debug_assert( saferes_ > pad );
	saferes_ -= pad;
	if ( cutpoint_ ) {
		debug_assert( cutpoint_ > pad );
		cutpoint_ -= pad;
	}
}

void
Segment::delete_trailing_residues()
{
	// nothing to do
	if ( cterm_included_ ) {
		return;
	}

	core::Size const pad = cterm_pad();
	debug_assert( length()-pad > 0 );
	std::string const newss = ss_.substr(0,length()-pad);
	ss_ = newss;
	utility::vector1< std::string > newabego( abego_.size() - pad );
	std::copy( abego_.begin(), abego_.end()-pad, newabego.begin() );
	abego_ = newabego;
	cterm_included_ = true;
}

/// @brief expands this residue set to include the dummy trailing residues
void
Segment::engulf_leading_residues()
{
	start_ = 1;
	nterm_included_ = true;
	// abego and ss stay the same
}

/// @brief expands this residue set to include the dummy trailing residues
void
Segment::engulf_trailing_residues()
{
	core::Size const pad = cterm_pad();
	stop_ += pad;
	cterm_included_ = true;
	// abego and ss stay the same
}

/// @brief given a residue number range local to this 1=start, length=end, delete the residue
void
Segment::delete_residues( core::Size const local_resnum_start, core::Size const local_resnum_stop )
{
	debug_assert( local_resnum_start >= 1 );
	debug_assert( local_resnum_stop <= stop_ );
	debug_assert( local_resnum_start <= local_resnum_stop );

	core::Size const len = local_resnum_stop - local_resnum_start + 1;
	debug_assert( len < stop_ );
	stop_ -= len;
	debug_assert( stop_ >= start_ );

	if ( local_resnum_stop <= cutpoint_ ) {
		cutpoint_ -= len;
		debug_assert( cutpoint_ >= start_ );
	} else if ( local_resnum_start <= cutpoint_ ) {
		cutpoint_ = local_resnum_start;
	}
	if ( local_resnum_stop <= saferes_ ) {
		saferes_ -= len;
		debug_assert( saferes_ >= start_ );
	} else if ( local_resnum_start <= saferes_ ) {
		saferes_ = local_resnum_start;
	}

	// fix secondary structure
	std::string newss = ss_.substr(0,local_resnum_start-1);
	newss += ss_.substr( local_resnum_stop, std::string::npos );
	ss_ = newss;

	// fix abego
	utility::vector1< std::string > newabego;
	for ( core::Size i=1; i<local_resnum_start; ++i ) {
		newabego.push_back( abego_[i] );
	}
	for ( core::Size i=local_resnum_stop+1, end=length()+len; i<=end; ++i ) {
		newabego.push_back( abego_[i] );
	}
	abego_ = newabego;
}

/// output residueinfo
std::ostream &
operator<<( std::ostream & os, Segment const & res )
{
	std::string ab_str = "";
	for ( core::Size i=1, end=res.abego().size(); i<=end; ++i ) {
		ab_str += res.abego()[i][0];
	}
	os << "start=\"" << res.start()
		<< "\" stop=\"" << res.stop()
		<< "\" safe=\"" << res.safe()
		<< "\" nterm=\"" << res.nterm_included_
		<< "\" cterm=\"" << res.cterm_included_
		<< "\" loop=\"" << res.is_loop
		<< "\" mgroup=\"" << res.movable_group
		<< "\" cutpoint=\"" << res.cutpoint()
		<< "\" ss=\"" << res.ss()
		<< "\" abego=\"" << ab_str;

	if ( res.lower_segment() != "" ) {
		os << "\" lower_segment=\"" << res.lower_segment();
	}
	if ( res.upper_segment() != "" ) {
		os << "\" upper_segment=\"" << res.upper_segment();
	}

	os << "\"";
	return os;
}

std::ostream &
operator<<( std::ostream & os, NamedSegment const & resis )
{
	os << "<ResidueRange name=\"" << resis.first << "\" " << resis.second << " />";
	return os;
}

} // namespace components
} // namespace denovo_design
} // namespace protocols
