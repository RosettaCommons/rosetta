// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/devel/denovo_design/components/StructureData.cc
/// @brief StructureData of a segment -- basically interfaces between segment and pose
/// @details
/// @author Tom Linsky

//Unit Headers
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataCreator.hh>

//Project Headers
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/SegmentPairing.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/denovo_design/util.hh>

//Protocol Headers

//Core Headers
#include <core/chemical/VariantType.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/io/Remarks.hh>
#include <core/pose/symmetry/util.hh>

//Basic/Utility/Numeric Headers
#include <basic/Tracer.hh>
#include <utility>
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>

// Boost/ObjexxFCL Headers
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

//Auto Headers
#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/list.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

// C++ Headers
#include <algorithm>

#include <stack> // AUTO IWYU For stack
#include <utility/stream_util.hh> // AUTO IWYU For operator<<

static basic::Tracer TR( "protocols.denovo_design.components.StructureData" );

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace denovo_design {
namespace components {
using core::io::Remarks;

char const StructureData::DATA_DELIMETER = '#';

StructureData::StructureData():
	id_( "" ),
	ss_( "" ),
	abego_( "" ),
	pose_length_( 0 ),
	length_( 0 ),
	segments_(),
	aliases_(),
	covalent_bonds_(),
	segment_order_(),
	block_signals_( false )
{
}

StructureData::StructureData( std::string const & id_val ) :
	id_( id_val ),
	ss_( "" ),
	pose_length_( 0 ),
	length_( 0 ),
	segments_(),
	aliases_(),
	covalent_bonds_(),
	segment_order_(),
	block_signals_( false )
{
}

/// @brief destructors
StructureData::~StructureData() = default;

std::string
StructureData::class_name()
{
	return "StructureData";
}

StructureData::CacheableDataOP
StructureData::clone() const
{
	return utility::pointer::make_shared< StructureData >( *this );
}

std::string
StructureData::datatype() const
{
	return class_name();
}

void
StructureData::write( std::ostream & os ) const
{
	std::stringstream xml;
	xml << *this;
	std::string xml_string = xml.str();
	clean_for_storage( xml_string );
	os << class_name() << " " << xml_string;
}

void
StructureData::parse_tag( utility::tag::TagCOP tag )
{
	std::string const id_str = tag->getOption< std::string >( "name", "" );
	if ( ! id_str.empty() ) {
		set_id( id_str );
	}
	for ( const auto & t : tag->getTags() ) {
		parse_subtag( t );
	}

	auto const length_check = tag->getOption< core::Size >( "pose_length", 0 );
	if ( length_check && ( pose_length_ != length_check ) ) {
		std::stringstream msg;
		msg << "StructureData::parse_tag(): Pose length in tag (" << length_check
			<< ") does not match length computed from Segment lengths ("
			<< pose_length_ << ")" << std::endl;
		msg << *this << std::endl;
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  msg.str() );
	}

	auto const elem_length_check = tag->getOption< core::Size >( "length", 0 );
	if ( elem_length_check && ( length_ != elem_length_check ) ) {
		std::stringstream msg;
		msg << "StructureData::parse_tag(): Element length in tag (" << elem_length_check
			<< ") does not match length computed from Segment lengths ("
			<< length_ << ")" << std::endl;
		msg << *this << std::endl;
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  msg.str() );
	}
}

void
StructureData::parse_subtag( utility::tag::TagCOP tag )
{
	if ( tag->getName() == "ResidueRange" ) {
		debug_assert( tag->hasOption( "name" ) );
		Segment newresis( tag->getOption< std::string >( "name" ) );
		newresis.parse_tag( tag );
		add_segment( newresis );
	} else if ( tag->getName() == "Int" ) {
		debug_assert( tag->hasOption( "name" ) );
		debug_assert( tag->hasOption( "value" ) );
		set_data_int( tag->getOption< std::string >( "name" ), tag->getOption< int >( "value" ) );
	} else if ( tag->getName() == "Real" ) {
		debug_assert( tag->hasOption( "name" ) );
		debug_assert( tag->hasOption( "value" ) );
		set_data_real( tag->getOption< std::string >( "name" ), tag->getOption< core::Real >( "value" ) );
	} else if ( tag->getName() == "Str" ) {
		debug_assert( tag->hasOption( "name" ) );
		if ( tag->hasOption( "value" ) ) {
			set_data_str( tag->getOption< std::string >( "name" ), tag->getOption< std::string >( "value" ) );
		} else {
			set_data_str( tag->getOption< std::string >( "name" ), "" );
		}
	} else if ( tag->getName() == "Alias" ) { debug_assert( tag->hasOption( "name" ) );
		debug_assert( tag->hasOption( "segment" ) );
		debug_assert( tag->hasOption( "res" ) );
		set_alias( tag->getOption< std::string >( "name" ),
			tag->getOption< std::string >( "segment" ),
			tag->getOption< core::Size >( "res" ) );
	} else if ( tag->getName() == "CovalentBond" ) {
		debug_assert( tag->hasOption( "segment1" ) );
		debug_assert( tag->hasOption( "segment2" ) );
		debug_assert( tag->hasOption( "residue1" ) );
		debug_assert( tag->hasOption( "residue2" ) );
		debug_assert( tag->hasOption( "atom1" ) );
		debug_assert( tag->hasOption( "atom2" ) );
		add_covalent_bond(
			tag->getOption< std::string >( "segment1" ),
			tag->getOption< core::Size >( "residue1" ),
			tag->getOption< std::string >( "atom1" ),
			tag->getOption< std::string >( "segment2" ),
			tag->getOption< core::Size >( "residue2" ),
			tag->getOption< std::string >( "atom2" )
		);
	} else if ( tag->getName() == SegmentPairing::TAG_NAME ) {
		SegmentPairingOP newpairing = SegmentPairing::create( tag->getOption< std::string >( "type" ) );
		newpairing->parse_my_tag( *tag );
		add_pairing( *newpairing );
	} else {
		tag->write( TR.Error );
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "Unknown tag in permutation: " + tag->getName() );
	}
}

/// @brief returns the id of this permutation
std::string const &
StructureData::id() const
{
	return id_;
}

/// @brief returns the length of this permutation
core::Size
StructureData::length() const
{
	return length_;
}

/// @brief returns the total length of this permutation, including n-, c-terminal loop residues which are basically for show
core::Size
StructureData::pose_length() const
{
	return pose_length_;
}

/// @brief sets id name
void
StructureData::set_id( std::string const & id_str )
{
	id_ = id_str;
}

std::string
StructureData::specific_enzdes_header(
	core::pose::Pose const & pose,
	std::string const & generic_remark_str ) const
{
	utility::vector1< std::string > const fields = utility::string_split_simple( generic_remark_str, ' ' );

	// compute chains/resids
	std::stringstream resid1_stream( fields[5] );
	auto const resid1 = boost::lexical_cast< core::Size >( substitute_variables( resid1_stream ) );
	std::string chain1 = "A";
	if ( resid1 ) {
		chain1 = core::conformation::canonical_chain_letter_for_chain_number( resid1 );
		if ( pose.pdb_info() ) {
			std::string pdbchain = pose.pdb_info()->chain( resid1 );
			if ( pdbchain != "^" ) chain1 = pdbchain;
		}
	}

	std::stringstream resid2_stream( fields[10] );
	auto const resid2 = boost::lexical_cast< core::Size >( substitute_variables( resid2_stream ) );
	std::string chain2 = "A";
	if ( resid2 ) {
		chain2 = core::conformation::canonical_chain_letter_for_chain_number( resid2 );
		if ( pose.pdb_info() ) {
			std::string pdbchain = pose.pdb_info()->chain( resid2 );
			if ( pdbchain != "^" ) chain2 = pdbchain;
		}
	}

	std::stringstream ss;
	ss << fields[1] << " " << fields[2] << " ";
	if ( resid1 ) ss << chain1;
	else ss << fields[3];
	ss << " " << fields[4] << " " << std::setw( 4 ) << resid1 << " ";
	ss << fields[6] << " " << fields[7] << " ";
	if ( resid2 ) ss << chain2;
	else ss << fields[8];
	ss << " " << fields[9] << " " << std::setw( 4 ) << resid2;
	for ( core::Size field_idx=11; field_idx<=fields.size(); ++field_idx ) {
		ss << " " << std::setw( 2 ) << fields[ field_idx ];
	}
	ss << "     ";

	return ss.str();
}

std::string
StructureData::generic_enzdes_header( std::string const & remark_str ) const
{
	//enzdes header
	utility::vector1< std::string > fields = utility::string_split_simple( remark_str, ' ' );
	auto const resid1 = boost::lexical_cast< core::Size >( fields[5] );
	std::stringstream ss;
	ss << fields[1] << " " << fields[2] << " " << fields[3] << " " << fields[4] << " ";
	if ( resid1 ) {
		std::string const seg = segment_name( resid1 );
		debug_assert( resid1 >= segment(seg).start() );
		debug_assert( resid1 <= segment(seg).stop() );
		core::Size const localres = segment(seg).pose_to_segment( resid1 );
		ss << "%%" << seg << DATA_DELIMETER << localres << "%% ";
	} else {
		ss << fields[5] << " ";
	}
	ss << fields[6] << " " << fields[7] << " " << fields[8] << " " << fields[9] << " ";

	auto const resid2 = boost::lexical_cast< core::Size >( fields[10] );
	if ( resid2 ) {
		std::string const seg2 = segment_name( resid2 );
		debug_assert( resid2 >= segment(seg2).start() );
		debug_assert( resid2 <= segment(seg2).stop() );
		core::Size const localres2 = segment(seg2).pose_to_segment( resid2 );
		ss << "%%" << seg2 << DATA_DELIMETER << localres2 << "%% ";
	} else {
		ss << fields[10] << " ";
	}
	ss << fields[11] << " " << fields[12];
	return ss.str();
}

/// @brief Saves remarks of the given pose into the pose's datacache -- changes enzdes residues to segment name/number
void
StructureData::save_remarks( core::io::Remarks const & remarks )
{
	core::Size remcount = 1;
	block_signals();
	for ( const auto & remark : remarks ) {
		// skip remarks related to StructureData
		if ( remark.num == StructureDataFactory::REMARK_NUM ) {
			continue;
		}

		std::stringstream remark_tag;
		remark_tag << "REMARK_" << remark.num;

		if ( remark.num == 666 ) {
			TR << "Processing enzdes header " << remark.value << std::endl;
			set_data_str( remark_tag.str(), boost::lexical_cast< std::string >( remcount ), generic_enzdes_header( remark.value ) );
		} else {
			set_data_str( remark_tag.str(), boost::lexical_cast< std::string >( remcount ), remark.value );
		}
		++remcount;
	}
	unblock_signals();
	changed();
}

/// @brief retrieves cached remarks specific to the given pose
core::io::Remarks
StructureData::retrieve_remarks( core::pose::Pose const & pose ) const
{
	using RemarkStrings = utility::vector1<std::string>;
	RemarkStrings numbers;
	RemarkStrings remarks;
	for ( const auto & s : data_str_ ) {
		if ( ! boost::starts_with( s.first, "REMARK_" ) ) continue;

		utility::vector1< std::string > const fields = utility::string_split( s.first, '_' );
		if ( fields.size() < 2 ) {
			std::stringstream msg;
			msg << "StructureData::retrieve_remarks(): invalid remark name " << s.first
				<< " remark value = " << s.second << std::endl;
			msg << *this << std::endl;
			utility_exit_with_message( msg.str() );
		}

		numbers.push_back( fields[2] );
		remarks.push_back( s.second );
	}

	debug_assert( numbers.size() == remarks.size() );

	using RemarkInfoOP = utility::pointer::shared_ptr<core::io::RemarkInfo>;
	using RemarkInfos = utility::vector1<RemarkInfoOP>;
	RemarkInfos ordered_remarks( numbers.size() );
	for ( RemarkStrings::const_iterator n=numbers.begin(), r=remarks.begin(); (n!=numbers.end()) && (r!=remarks.end()); ++n, ++r ) {
		utility::vector1< std::string > const fields = utility::string_split( *n, DATA_DELIMETER );
		if ( fields.size() < 2 ) {
			std::stringstream msg;
			msg << "StructureData::retrieve_remarks(): invalid remark number " << *n
				<< " remark value = " << *r << std::endl;
			msg << *this << std::endl;
			utility_exit_with_message( msg.str() );
		}
		auto const remark_num = boost::lexical_cast< core::Size >( fields[1] );
		auto const remark_order = boost::lexical_cast< core::Size >( fields[2] );

		if ( remark_order > ordered_remarks.size() ) {
			std::stringstream msg;
			msg << "StructureData::retrieve_remarks(): Remark order number (" << remark_order
				<< ") is large than the number of remarks found (" << ordered_remarks.size()
				<< " -- remarks must be stored with contiguous numbers starting with 1." << std::endl;
			msg << *this << std::endl;
			utility_exit_with_message( msg.str() );
		}

		TR.Debug << "Cached line = " << *r << std::endl;
		RemarkInfoOP me( new core::io::RemarkInfo );
		me->num = remark_num;
		if ( me->num == 666 ) { // enzdes header
			me->value = specific_enzdes_header( pose, *r );
		} else {
			me->value = *r;
		}

		if ( ordered_remarks[ remark_order ] ) {
			std::stringstream msg;
			msg << "StructureData::retrieve_remarks(): Remark order number (" << remark_order
				<< ") exists twice.  Each remark order number can only be used once." << std::endl;
			msg << *this << std::endl;
			utility_exit_with_message( msg.str() );
		}

		ordered_remarks[ remark_order ] = me;
	}

	core::io::Remarks retval;
	for ( RemarkInfos::const_iterator r=ordered_remarks.begin(); r!=ordered_remarks.end(); ++r ) {
		if ( ! *r ) {
			std::stringstream msg;
			msg << "StructureData::retrieve_remarks(): A remark order number between 1 and " << ordered_remarks.size()
				<< " has not been used. Each remark order number can only be used once." << std::endl;
			msg << *this << std::endl;
			utility_exit_with_message( msg.str() );
		}
		retval.push_back( **r );
	}
	return retval;
}

/// @brief tells whether a non-polymer bond exists between the given segments
bool
StructureData::non_polymer_bond_exists( std::string const & seg1, std::string const & seg2 ) const
{
	return ( find_non_polymer_bond( seg1, seg2 ) != covalent_bonds_.end() );
}

BondInfo const &
StructureData::non_polymer_bond( std::string const & seg1, std::string const & seg2 ) const
{
	auto bi = find_non_polymer_bond( seg1, seg2 );
	if ( bi == covalent_bonds_end() ) {
		std::stringstream msg;
		msg << "StructureData::non_peptidic_bond(): A non-polymer bond between segments "
			<< seg1 << " and " << seg2 << " could not be found." << std::endl;
		msg << *this << std::endl;
		utility_exit_with_message( msg.str() );
	}

	return *bi;
}

BondInfos::const_iterator
StructureData::find_non_polymer_bond( std::string const & seg1, std::string const & seg2 ) const
{
	for ( auto bi=covalent_bonds_begin(); bi!=covalent_bonds_end(); ++bi ) {
		if ( ( ( seg1 == bi->seg1 ) && ( seg2 == bi->seg2 ) ) ||
				( ( seg2 == bi->seg1 ) && ( seg1 == bi->seg2 ) ) ) {
			return bi;
		}
	}
	return covalent_bonds_end();
}

BondInfos::const_iterator
StructureData::covalent_bonds_begin() const
{
	return covalent_bonds_.begin();
}

BondInfos::const_iterator
StructureData::covalent_bonds_end() const
{
	return covalent_bonds_.end();
}

std::ostream & operator<<( std::ostream & os, BondInfo const & b ) {
	os << "\t<CovalentBond segment1=\"" << b.seg1 << "\" segment2=\"" << b.seg2
		<< "\" residue1=\"" << b.res1 << "\" residue2=\"" << b.res2
		<< "\" atom1=\"" << b.atom1 << "\" atom2=\"" << b.atom2 << "\" />";
	return os;
}


/// @brief given an input stream, substitute all variables
/// @details variables are of the form: %%SEGMENTNAME#residue%%
/// SEGMENTNAME = name of the segment
/// residue = local residue number within the segment
/// The substituted value will be an core::Size corresponding to the pose residue
std::string
StructureData::substitute_variables( std::istream & input ) const
{
	std::stringstream sub_str;
	core::Size linecount = 0;
	// File format: %%SEGMENT_NAME#resid%% will give an integer resid
	while ( input.good() ) {
		std::string line = "";
		std::getline( input, line );
		if ( line[0] == '#' ) continue;
		++linecount;
		TR.Debug << "Line=" << line << std::endl;
		core::Size next_sub = line.find("%%");
		while ( next_sub != std::string::npos ) {
			core::Size second_sub = line.find("%%", next_sub+1);
			if ( second_sub == std::string::npos ) {
				throw CREATE_EXCEPTION(utility::excn::BadInput,  "Malformed line in constraint file : " + line );
			}
			debug_assert( second_sub - next_sub >= 5 );
			std::string const variable = line.substr( next_sub+2, second_sub-next_sub-2 );
			utility::vector1< std::string > fields = utility::string_split( variable, DATA_DELIMETER );
			core::Size local_resid = 0;
			bool from_start = true;
			if ( fields.size() > 1 ) {
				std::string res_str = fields[ 2 ];
				if ( *(res_str.begin()) == '-' ) {
					from_start = false;
					res_str = res_str.substr( 1, std::string::npos );
				}
				local_resid = boost::lexical_cast< core::Size >( res_str );
			}
			core::Size actual_resid = 0;
			if ( from_start ) {
				actual_resid = pose_residue( fields[ 1 ], local_resid );
			} else {
				actual_resid = segment( fields[ 1 ] ).stop() - local_resid + 1;
				debug_assert( actual_resid >= segment( fields[ 1 ] ).start() );
			}

			line = line.substr(0,next_sub) + boost::lexical_cast< std::string >(actual_resid) + line.substr(second_sub+2,std::string::npos);
			next_sub = line.find("%%");
		}
		TR.Debug << "New line=" << line << std::endl;
		if ( linecount > 1 ) {
			sub_str << std::endl;
		}
		if ( !line.empty() ) {
			sub_str << line;
		}
	}
	return sub_str.str();
}

/// @brief finds a segment in the segment_order list and returns an iterator to it
/// @throws utility_exit if not found
SegmentNameList::iterator
StructureData::find_segment_name( std::string const & segname )
{
	auto c = std::find( segment_order_.begin(), segment_order_.end(), segname );
	if ( c == segment_order_.end() ) {
		std::stringstream msg;
		msg << "StructureData::find_segment_name(): Could not find a segment named " << segname
			<< ". Available segments are: " << segment_order_ << std::endl;
		msg << *this << std::endl;
		utility_exit_with_message( msg.str() );
	}
	return c;
}

/// @brief finds a segment in the segment_order list and returns a const iterator to it
/// @throws utility_exit if not found
SegmentNameList::const_iterator
StructureData::find_segment_name( std::string const & segname ) const
{
	auto c = std::find( segment_order_.begin(), segment_order_.end(), segname );
	if ( c == segment_order_.end() ) {
		std::stringstream msg;
		msg << "StructureData::find_segment_name(): Could not find a segment named " << segname
			<< ". Available segments are: " << segment_order_ << std::endl;
		msg << *this << std::endl;
		utility_exit_with_message( msg.str() );
	}
	return c;
}

/// @brief finds a segment in the segments map and returns an iterator to it
/// @throws utility_exit if not found
SegmentMap::iterator
StructureData::find_segment( std::string const & segname )
{
	auto s = segments_.find( segname );
	if ( s == segments_.end() ) {
		std::stringstream msg;
		msg << "StructureData::find_segment(): Could not find a segment named " << segname
			<< ". Available segments are: " << segment_order_ << std::endl;
		msg << *this << std::endl;
		utility_exit_with_message( msg.str() );
	}
	return s;
}

/// @brief finds a segment in the segments map and returns an iterator to it
/// @throws utility_exit if not found
SegmentMap::const_iterator
StructureData::find_segment( std::string const & segname ) const
{
	auto c = segments_.find( segname );
	if ( c == segments_.end() ) {
		std::stringstream msg;
		msg << "StructureData::find_segment(): Could not find a segment named " << segname
			<< ". Available segments are: " << segment_order_ << std::endl;
		msg << *this << std::endl;
		utility_exit_with_message( msg.str() );
	}
	return c;
}

SegmentNameList
StructureData::connected_segments( std::string const & seg, bool const stop_at_cutpoint ) const
{
	SegmentNameList segmentlist;
	std::set< std::string > visited;
	std::stack< std::string > nodes;

	// go forward from segment
	nodes.push( seg );
	while ( nodes.size() ) {
		std::string const cur = nodes.top();
		nodes.pop();
		if ( visited.find( cur ) != visited.end() ) {
			continue;
		}
		visited.insert( cur );
		segmentlist.push_back( cur );
		Segment const & res = segment( cur );
		if ( !res.has_free_upper_terminus() ) {
			if ( !stop_at_cutpoint || !res.cutpoint() ) {
				nodes.push( res.upper_segment() );
			}
		}
	}

	// now go backward from segment
	if ( !segment( seg ).has_free_lower_terminus() ) {
		nodes.push( segment( seg ).lower_segment() );
	}
	while ( nodes.size() ) {
		std::string const cur = nodes.top();
		nodes.pop();
		if ( visited.find( cur ) != visited.end() ) {
			continue;
		}
		visited.insert( cur );
		segmentlist.push_front( cur );
		Segment const & res = segment( cur );
		if ( !res.has_free_lower_terminus() ) {
			if ( !stop_at_cutpoint || !res.cutpoint() ) {
				nodes.push( res.lower_segment() );
			}
		}
	}
	return segmentlist;
}

/// @brief marks the given segments as covanlently connected
void
StructureData::mark_connected(
	std::string const & lower_seg,
	std::string const & upper_seg )
{
	auto low = segments_.find( lower_seg );
	debug_assert( low != segments_.end() );
	auto up = segments_.find( upper_seg );
	debug_assert( up != segments_.end() );
	debug_assert( low->second.upper_segment() == "" );
	low->second.set_upper_segment( upper_seg );
	debug_assert( up->second.lower_segment() == "" );
	up->second.set_lower_segment( lower_seg );
	changed();
}

/// @brief removes jump and cutpoint between the two segments to create a single polymer chain
void
StructureData::delete_jump_and_intervening_cutpoint( std::string const & segment1, std::string const & segment2 )
{
	if ( segment(segment1).cutpoint() ) {
		segment_nonconst( segment1 ).set_cutpoint( 0 );
	} else {
		if ( segment(segment2).cutpoint() ) {
			segment_nonconst( segment2 ).set_cutpoint( 0 );
		}
	}
	changed();
}

/// @brief return secondary structure string
std::string const &
StructureData::ss() const
{
	return ss_;
}

/// @brief return secondary structure at residue resid
char
StructureData::ss( core::Size const resid ) const
{
	return ss_[ resid - 1 ];
}

/// @brief sets secondary structure for residue resid
void
StructureData::set_ss( core::Size const resid, char const ss_type )
{
	std::string const seg_name = segment_name( resid );
	core::Size const seg_resid = segment( seg_name ).pose_to_segment( resid );
	TR << "Found residue " << resid << " : " << seg_name << "#" << seg_resid << std::endl;
	segment_nonconst( seg_name ).set_ss( seg_resid, ss_type );
	update_numbering();
}

/// @brief return abego string
std::string const &
StructureData::abego() const
{
	return abego_;
}

/// @brief return secondary structure at residue resid
char
StructureData::abego( core::Size const resid ) const
{
	return abego_[ resid - 1 ];
}

void
StructureData::set_abego( std::string const & segment, std::string const & abego )
{
	segment_nonconst( segment ).set_abego( abego );
	update_numbering();
}

void
StructureData::set_abego( std::string const & segment, utility::vector1< std::string > const & abego )
{
	set_abego( segment, abego_str( abego ) );
}

void
StructureData::set_abego( std::string const & abego )
{
	for ( auto s=segments_begin(); s!=segments_end(); ++s ) {
		Segment & seg = segment_nonconst( *s );
		std::string seg_abego = "";
		for ( core::Size resid=seg.lower(); resid<=seg.upper(); ++resid ) {
			seg_abego += abego[ resid - 1 ];
		}
		seg.set_abego( seg_abego );
	}
}

void
StructureData::set_abego( utility::vector1< std::string > const & abego )
{
	set_abego( abego_str( abego ) );
}

/// @brief re-arranges residues such that segment 2 follows segment 1 in sequence
void
StructureData::move_segment( std::string const & segment1, std::string const & segment2 )
{
	Segment const & s1 = segment( segment1 );
	Segment const & s2 = segment( segment2 );

	// don't do anything if the segments are already moved
	if ( s1.upper() + 1 == s2.lower() ) {
		return;
	}

	// can't move a segment to after itself
	if ( segment1 == segment2 ) {
		std::stringstream msg;
		msg << "StructureData::move_segment(): The segment1 (" << segment1
			<< ") that you are trying to move is the same as the segment2 (" << segment2 << ")" << std::endl;
		msg << "Moving a segment to after itself cannot be done." << std::endl;
		msg << *this << std::endl;
		utility_exit_with_message( msg.str() );
	}

	// check connectivity
	if ( !s1.has_free_upper_terminus() && ( s1.upper_segment() != segment2 ) ) {
		std::stringstream msg;
		msg << "StructureData::move_segment(): The segment1 you are trying to move ("
			<< segment1 << ") has a different upper connected segment listed ("
			<< s1.upper_segment() << ") than the one being moved after it ("
			<< segment2 << ")." << std::endl;
		msg << *this << std::endl;
		utility_exit_with_message( msg.str() );
	}

	// check connections
	if ( !s2.has_free_lower_terminus() && ( s2.lower_segment() != segment1 ) ) {
		std::stringstream msg;
		msg << "StructureData::move_segment(): The segment2 you are trying to move ("
			<< segment2 << ") has a different lower connected segment listed ("
			<< s2.upper_segment() << ") than the one it is being moved to follow ("
			<< segment1 << ")." << std::endl;
		msg << *this << std::endl;
		utility_exit_with_message( msg.str() );
	}

	// list of segments connected to segment2 to be moved
	SegmentNameList const segments_to_move = connected_segments( segment2, false );

	// ensure segment2 is the first segment in the list
	if ( *segments_to_move.begin() != segment2 ) {
		std::stringstream msg;
		msg << "StructureData::move_segment(): The segment2 you are trying to move ("
			<< segment2 << ") is not the first segment in its chain. "
			<< "Segments in the same chain : " << segments_to_move << std::endl;
		msg << *this << std::endl;
		utility_exit_with_message( msg.str() );
	}

	move_segments( segment1, *segments_to_move.begin(), *segments_to_move.rbegin() );
}

void
StructureData::move_segments(
	std::string const & segment1,
	std::string const & segment2_lower,
	std::string const & segment2_upper )
{
	// don't do anything if the segments are already moved
	if ( segment( segment1 ).upper() + 1 == segment( segment2_lower ).lower() ) {
		return;
	}

	// get segment range to move
	auto segment2_segment_begin = find_segment_name( segment2_lower );
	auto segment2_segment_end = find_segment_name( segment2_upper );
	debug_assert( segment2_segment_begin != segment_order_.end() );
	debug_assert( segment2_segment_end != segment_order_.end() );
	++segment2_segment_end;

	// Store segments to be moved
	SegmentNameList const segments_to_move( segment2_segment_begin, segment2_segment_end );

	// remove segments in segment2
	debug_assert( segment_order_.size() == segments_.size() );
	segment_order_.erase( segment2_segment_begin, segment2_segment_end );

	// find insertion position
	auto insert_pos = find_segment_name( segment1 );
	debug_assert( insert_pos != segment_order_.end() );
	++insert_pos;

	// insert the segment2 segments after segment1
	segment_order_.insert( insert_pos, segments_to_move.begin(), segments_to_move.end() );
	debug_assert( segment_order_.size() == segments_.size() );
	update_numbering();
}

/// @brief deletes the residues between the segment N terminus and the N anchor point
void
StructureData::delete_leading_residues( std::string const & seg )
{
	segment_nonconst( seg ).delete_lower_padding();
	update_numbering();
}

/// @brief deletes the residues between the segment C terminus and the C anchor point
void
StructureData::delete_trailing_residues( std::string const & seg )
{
	segment_nonconst( seg ).delete_upper_padding();
	update_numbering();
}

void
StructureData::delete_segment( std::string const & seg_val )
{
	TR.Debug << "StructureData::delete_segment(" << seg_val << ")" << std::endl;

	auto r = find_segment( seg_val );
	std::string const segment1_c = r->second.lower_segment();
	std::string const segment2_n = r->second.upper_segment();

	// Track size of deleted segment
	core::Size start_del = r->second.lower();
	core::Size stop_del = r->second.upper();

	if ( !r->second.has_free_lower_terminus() ) {
		auto r2 = find_segment( segment1_c );
		debug_assert( r2 != segments_.end() );
		if ( start_del <= stop_del ) {
			++start_del;
			r2->second.add_upper_padding();
		}
		r2->second.set_upper_segment( "" );
	}

	if ( !r->second.has_free_upper_terminus() ) {
		auto r2 = find_segment( segment2_n );
		debug_assert( r2 != segments_.end() );
		if ( start_del <= stop_del ) {
			--stop_del;
			r2->second.add_lower_padding();
		}
		r2->second.set_lower_segment( "" );
	}

	auto remove_me = find_segment_name( seg_val );
	debug_assert( remove_me != segment_order_.end() );
	segment_order_.erase( remove_me );
	segments_.erase( r );

	debug_assert( segment_order_.size() == segments_.size() );
	update_numbering();
}

/// @brief the total number of free termini in the permutation
core::Size
StructureData::num_chains() const
{
	core::Size lower_count = 0;
	core::Size upper_count = 0;
	for ( const auto & segment : segments_ ) {
		if ( segment.second.has_free_lower_terminus() ) {
			++lower_count;
		}
		if ( segment.second.has_free_upper_terminus() ) {
			++upper_count;
		}
	}

	if ( lower_count != upper_count ) {
		TR.Error << "Error in counting number of chains..." << lower_count << " vs. " << upper_count << std::endl;
		TR.Error << *this << std::endl;
	}
	debug_assert( lower_count == upper_count );
	return lower_count;
}

/// @brief returns the actual residue number of the given name and res #
core::Size
StructureData::pose_residue( std::string const & segment_name, core::Size const local_res ) const
{
	// first look in residues
	auto r = segments_.find( segment_name );
	if ( r == segments_.end() ) {
		if ( has_alias( segment_name ) ) {
			return alias( segment_name );
		}
	}
	if ( r == segments_.end() ) {
		std::stringstream err;
		err << "Can't resolve the segment name " << segment_name
			<< "into a valid segment in the permutation. Valid segments are: " << segments_ << std::endl;
		throw CREATE_EXCEPTION(utility::excn::BadInput,  err.str() );
	}
	return r->second.segment_to_pose(local_res);
}

/// @brief returns segments which have free lower termini in sequence order
SegmentNames
StructureData::available_lower_termini() const
{
	SegmentNames retval;
	for ( auto r=segments_begin(); r!=segments_end(); ++r ) {
		if ( segment( *r ).has_free_lower_terminus() ) {
			TR.Debug << *r << " has free lower terminus: " << segment( *r ).lower_segment() << std::endl;
			retval.push_back( *r );
		}
	}
	return retval;
}

/// @brief returns segments which have free upper termini in sequence order
SegmentNames
StructureData::available_upper_termini() const
{
	SegmentNames retval;
	for ( auto r=segments_begin(); r!=segments_end(); ++r ) {
		if ( segment( *r ).has_free_upper_terminus() ) {
			TR.Debug << *r << " has free upper terminus: " << segment( *r ).upper_segment() << std::endl;
			retval.push_back( *r );
		}
	}
	return retval;
}

/// @brief replace segment named seg_name with the given new segment
void
StructureData::replace_segment( std::string const & seg_name, Segment const & segment )
{
	auto s = find_segment( seg_name );
	if ( s == segments_.end() ) {
		std::stringstream msg;
		msg << "StructureData::replace_segment(): Segment named " << seg_name
			<< " was not found in the StructureData.  SD=" << *this << std::endl;
		utility_exit_with_message( msg.str() );
	}
	s->second = segment;
	update_numbering();
	changed();
}

/// @brief adds a residues segment to the end of the list
void
StructureData::add_segment( Segment const & resis )
{
	TR.Debug << "Adding " << resis << " to end of list" << std::endl;
	add_segment( resis, segment_order_.end() );
}

/// @brief adds a residues segment -- final will be inserted before segment named "insert_segment"
/// if "after_seg" is empty, or not found in the list of segments, the new segment will be inserted at the end
void
StructureData::add_segment( Segment const & resis, std::string const & insert_before_segment )
{
	add_segment( resis, find_segment_name( insert_before_segment ) );
}

/// @brief adds a residues segment -- final will be inserted before the iterator given as insert_pose
void
StructureData::add_segment(
	Segment const & resis,
	SegmentNameList::iterator insert_pos )
{
	if ( insert_pos != segment_order_.end() ) {
		TR.Debug << "Adding " << resis << " to after " << *insert_pos << std::endl;
	} else {
		TR.Debug << "Adding " << resis << " to end of list" << std::endl;
	}
	debug_assert( segment_order_.size() == segments_.size() );

	std::string const & id_val = resis.id();
	auto s = segments_.find( id_val );
	if ( s == segments_.end() ) {
		if ( !segments_.insert( std::make_pair( id_val, resis ) ).second ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "failed to insert segment " + id_val );
		}
		segment_order_.insert( insert_pos, id_val );
		debug_assert( segments_.size() == segment_order_.size() );
	} else {
		s->second = resis;
		debug_assert( std::count( segment_order_.begin(), segment_order_.end(), id_val ) == 1 );
	}
	TR.Debug << "Segment order is now " << segment_order_ << std::endl;
	update_numbering();
}

/// @brief returns n and c terminal segments of the chain which includes res
std::pair< std::string, std::string >
StructureData::termini( std::string const & seg ) const
{
	// go forward
	std::set< std::string > visited;
	std::string lower_segment = find_lower_terminus( visited, seg );
	TR.Debug << "Lower segment is " << lower_segment << " visited=" << visited << std::endl;
	visited.clear();
	std::string upper_segment = find_upper_terminus( visited, seg );
	TR.Debug << "Upper segment is " << upper_segment << " visited=" << visited << std::endl;
	return std::make_pair( lower_segment, upper_segment );
}

/// @brief performs dfs in lower direction looking for termini
std::string
StructureData::find_lower_terminus( std::set< std::string > & visited, std::string const & seg ) const
{
	// if we've already seen this segment, we are cyclic
	if ( visited.find(seg) != visited.end() ) {
		return "";
	}
	visited.insert(seg);
	if ( segment(seg).has_free_lower_terminus() ) {
		return seg;
	}
	debug_assert( segment(seg).lower_segment() != "" );
	return find_lower_terminus( visited, segment(seg).lower_segment() );
}

/// @brief performs dfs in upper direction looking for termini
std::string
StructureData::find_upper_terminus( std::set< std::string > & visited, std::string const & seg ) const
{
	// if we've already seen this segment, we are cyclic
	if ( visited.find(seg) != visited.end() ) {
		return "";
	}
	visited.insert(seg);
	if ( segment(seg).has_free_upper_terminus() ) {
		return seg;
	}
	debug_assert( segment(seg).upper_segment() != "" );
	return find_upper_terminus( visited, segment(seg).upper_segment() );
}

/// @brief sets an "alias" for a particular residue inside a segment which allows for it to be easily accessed
void
StructureData::set_alias(
	std::string const & alias_name,
	std::string const & segment_name,
	core::Size const resi )
{
	if ( has_segment( segment_name ) ) {
		debug_assert( resi );
		debug_assert( resi <= segment(segment_name).elem_length() );
		aliases_[ alias_name ] = Alias( segment_name, resi );
	} else if ( has_alias( segment_name ) ) {
		aliases_[ alias_name ] = aliases_[ segment_name ];
	} else {
		std::stringstream ss;
		ss << "Segment named " << segment_name
			<< " does not exist in the permutation as an alias or residue segment. Perm="
			<< *this << std::endl;
		throw CREATE_EXCEPTION(utility::excn::BadInput,  ss.str() );
	}
	changed();
}

/// @brief sets an "alias" for a particular residue which allows for it to be easily accessed
void
StructureData::set_alias(
	std::string const & alias_name,
	core::Size const resi )
{
	debug_assert( resi );
	debug_assert( resi <= pose_length() );
	std::string const segmentname = segment_name( resi );
	TR.Debug << "Going to set alias to " << resi << " with segment = " << segmentname << std::endl;
	core::Size const localresi = segment( segmentname ).pose_to_segment( resi );
	set_alias( alias_name, segmentname, localresi );
}

/// @brief given a residue alias, returns a pose residue number
core::Size
StructureData::alias( std::string const & alias ) const
{
	auto alias_it = aliases_.find( alias );
	if ( alias_it == aliases_.end() ) {
		TR.Warning << " alias " << alias << " not found!!  returning 0" << std::endl;
	}

	std::string const & segment_name = alias_it->second.first;
	core::Size local_res = alias_it->second.second;
	TR.Debug << "Looking up Alias " << segment_name << " # " << local_res << std::endl;
	debug_assert( has_segment( segment_name ) );
	return segment( segment_name ).segment_to_pose( local_res );
}

/// @brief renames a residue segment and updates all connections. If the prefix is already there for a given segment name, it does nothing
void
StructureData::add_prefix_to_segments( std::string const & prefix )
{
	add_prefix_to_segments( prefix, PARENT_DELIMETER );
}

/// @brief renames a residue segment and updates all connections. If the prefix is already there for a given segment name, it does nothing
void
StructureData::add_prefix_to_segments( std::string const & prefix, char const delimeter )
{
	std::string const prepend_str = prefix + delimeter;
	SegmentMap newmap;
	for ( auto & segment : segments_ ) {
		if ( ! segment.second.has_free_lower_terminus() ) {
			if ( !boost::starts_with( segment.second.lower_segment(), prepend_str ) ) {
				segment.second.set_lower_segment( prepend_str + segment.second.lower_segment() );
			}
		}
		if ( ! segment.second.has_free_upper_terminus() ) {
			if ( !boost::starts_with( segment.second.upper_segment(), prepend_str ) ) {
				segment.second.set_upper_segment( prepend_str + segment.second.upper_segment() );
			}
		}
		std::string newname;
		if ( boost::starts_with( segment.first, prepend_str ) ) {
			newname = segment.first;
		} else {
			newname = prepend_str + segment.first;
		}
		Segment newseg = segment.second;
		newseg.set_id( newname );
		newmap.insert( std::make_pair( newname, newseg ) );
		if ( segment.first != newname ) {
			auto c = find_segment_name( segment.first );
			debug_assert( c != segment_order_.end() );
			*c = newname;
		}
	}
	// fix aliases too
	for ( auto & aliase : aliases_ ) {
		if ( !boost::starts_with( aliase.second.first, prepend_str ) ) {
			aliase.second.first = prepend_str + aliase.second.first;
		}
	}
	segments_ = newmap;
	changed();
}

/// @brief renames a residue segment and updates all connections
void
StructureData::rename_segment( std::string const & old_name, std::string const & new_name )
{
	// find and rename the original segment
	auto r = find_segment( old_name );
	Segment resis = r->second;
	resis.set_id( new_name );
	segments_.erase( r );
	segments_.insert( std::make_pair( new_name, resis ) );

	auto c = find_segment_name( old_name );
	*c = new_name;

	for ( auto & segment : segments_ ) {
		if ( segment.second.lower_segment() == old_name ) {
			segment.second.set_lower_segment( new_name );
		}
		if ( segment.second.upper_segment() == old_name ) {
			segment.second.set_upper_segment( new_name );
		}
	}

	// fix aliases too
	for ( auto & aliase : aliases_ ) {
		if ( aliase.second.first == old_name ) aliase.second.first = new_name;
	}
	changed();
}

/// @brief true if this permutation contains a residue segment named seg
bool
StructureData::has_segment( std::string const & seg ) const
{
	return ( ( segments_.find(seg) != segments_.end() ) ||
		( segments_.find( id() + PARENT_DELIMETER + seg ) != segments_.end() ) );
}

/// @brief start of segments list
SegmentNameList::const_iterator
StructureData::segments_begin() const
{
	return segment_order_.begin();
}

/// @brief end of segment list
SegmentNameList::const_iterator
StructureData::segments_end() const
{
	return segment_order_.end();
}

core::Size
StructureData::choose_new_movable_group() const
{
	std::set< core::Size > const mgs( movable_groups_.begin(), movable_groups_.end() );
	core::Size mg = 1;
	while ( mgs.find( mg ) != mgs.end() ) {
		++mg;
	}
	return mg;
}

/// @brief Returns a StructureData containing only the given segments
///        The resulting StructureData will contain all relevant data from the current
StructureData
StructureData::slice( SegmentNameSet const & segmentset, bool const force_padding ) const
{
	StructureData sd( id() );

	// Add desired segments
	for ( auto s=segments_begin(); s!=segments_end(); ++s ) {
		if ( segmentset.find( *s ) == segmentset.end() ) continue;
		Segment newseg = segment( *s );

		if ( segmentset.find( newseg.lower_segment() ) == segmentset.end() ) {
			if ( !newseg.lower_segment().empty() ) {
				newseg.set_lower_segment( "" );
				if ( force_padding || !newseg.template_pose() ) newseg.add_lower_padding();
			}
		}
		if ( segmentset.find( newseg.upper_segment() ) == segmentset.end() ) {
			if ( !newseg.upper_segment().empty() ) {
				newseg.set_upper_segment( "" );
				if ( force_padding || !newseg.template_pose() ) newseg.add_upper_padding();
			}
		}
		sd.add_segment( newseg );
	}

	sd.copy_data( *this );
	return sd;
}

/// @brief merge all data and segments from "other" into this StructureData
void
StructureData::merge( StructureData const & other )
{
	if ( other.segments_begin() == other.segments_end() ) {
		std::stringstream msg;
		msg << "Not merging because the StructureData to be merged does not contain any segments!"
			<< " SD=" << *this << " Other SD=" << other << std::endl;
		utility_exit_with_message( msg.str() );
	}
	merge( other, SegmentNames( other.segments_begin(), other.segments_end() ) );
}


/// @brief merge given data and segments from "other" into this StructureData
void
StructureData::merge( StructureData const & other, SegmentNames const & segments )
{
	core::Size movable_group_template_pose = 0;

	for ( const auto & segment : segments ) {
		Segment newseg = other.segment( segment );

		// if there is no template pose, no pose exists to provide orientation.
		// in this case, the segment should be placed into its own movable group
		MovableGroup group;
		if ( !newseg.template_pose() ) {
			group = choose_new_movable_group();
		} else {
			if ( movable_group_template_pose == 0 ) {
				movable_group_template_pose = choose_new_movable_group();
			}
			group = movable_group_template_pose;
		}
		newseg.set_movable_group( group );

		add_segment( newseg, segment_order_.end() );
		TR.Debug << "Added " << newseg << std::endl;
	}

	copy_data( other );
}

void
StructureData::merge_before( StructureData const & other, std::string const & position, SegmentNames const & segments )
{
	core::Size movable_group_template_pose = 0;
	for ( const auto & segment : segments ) {
		if ( has_segment( segment ) ) {
			std::stringstream msg;
			msg << "StructureData::merge_before(): StructureData contains a segment named "
				<< segment << ", but one of the segments being merged has the same name. "
				<< "Segment names must be unique." << std::endl
				<< "Current SD = " << *this << std::endl
				<< "Other SD = " << other << std::endl;
			utility_exit_with_message( msg.str() );
		}
		Segment newseg = other.segment( segment );

		// if there is no template pose, no pose exists to provide orientation.
		// in this case, the segment should be placed into its own movable group
		MovableGroup group;
		if ( !newseg.template_pose() ) {
			group = choose_new_movable_group();
		} else {
			if ( movable_group_template_pose == 0 ) {
				movable_group_template_pose = choose_new_movable_group();
			}
			group = movable_group_template_pose;
		}
		newseg.set_movable_group( group );

		add_segment( newseg, find_segment_name( position ) );
		TR.Debug << "Added " << newseg << std::endl;
	}

	copy_data( other );
}

void
StructureData::merge_before( StructureData const & other, std::string const & position )
{
	merge_before( other, position, SegmentNames( other.segments_begin(), other.segments_end() ) );
}

/// @brief computes and returns a set of segments which are in the given movable group
SegmentNames
StructureData::segments_in_movable_group( core::Size const group ) const
{
	SegmentNames segments;
	for ( const auto & segment : segments_ ) {
		if ( segment.second.movable_group() == group ) {
			segments.push_back( segment.first );
		}
	}
	return segments;
}

/// @brief computes and returns a set of movable groups
MovableGroups const &
StructureData::movable_groups() const
{
	return movable_groups_;
}

void
StructureData::clear_pairings()
{
	pairings_.clear();
}

void
StructureData::add_pairing( SegmentPairing const & pairing )
{
	pairings_.push_back( pairing.clone() );
}

SegmentPairingCOPs::const_iterator
StructureData::pairings_begin() const
{
	return pairings_.begin();
}

SegmentPairingCOPs::const_iterator
StructureData::pairings_end() const
{
	return pairings_.end();
}

/// @brief declares a covalent bond between the specified atoms
void
StructureData::declare_covalent_bond(
	std::string const & seg1, core::Size const res1, std::string const & atom1,
	std::string const & seg2, core::Size const res2, std::string const & atom2 )
{
	declare_covalent_bond( pose_residue( seg1, res1 ), atom1, pose_residue( seg2, res2 ), atom2 );
}

/// @brief declares a covalent bond using pose residues
void
StructureData::declare_covalent_bond(
	core::Size const res1, std::string const & atom1,
	core::Size const res2, std::string const & atom2 )
{
	TR.Debug << "Creating covalent bond between " << res1 << " and " << res2 << std::endl;
	if ( ( res1 != res2 + 1 ) && ( res1 + 1 != res2 ) ) {
		add_covalent_bond( res1, atom1, res2, atom2 );
	}
}

void
StructureData::add_covalent_bond(
	core::Size const res1, std::string const & atom1,
	core::Size const res2, std::string const & atom2 )
{
	std::string const & seg1 = segment_name( res1 );
	core::Size const localres1 = res1 - segment( seg1 ).start() + 1;
	std::string const & seg2 = segment_name( res2 );
	core::Size const localres2 = res2 - segment( seg2 ).start() + 1;
	if ( localres1 && localres2 ) {
		add_covalent_bond( seg1, localres1, atom1, seg2, localres2, atom2 );
	} else {
		TR.Warning << "connection between residues " << res1 << " and " << res2 << " may be lost since one/both of the residues are present as \"padding\" residues." << std::endl;
	}
}

void
StructureData::add_covalent_bond(
	std::string const & seg1, core::Size const res1, std::string const & atom1,
	std::string const & seg2, core::Size const res2, std::string const & atom2 )
{
	add_covalent_bond( BondInfo( seg1, seg2, res1, res2, atom1, atom2 ) );
}

void
StructureData::add_covalent_bond( BondInfo const & bi )
{
	// look for this bond, stop if it's found
	for ( utility::vector1< BondInfo >::const_iterator bi_it=covalent_bonds_.begin(); bi_it != covalent_bonds_.end(); ++bi_it ) {
		if ( bi == *bi_it ) {
			TR.Debug << "Skipping adding existing bond info " << *bi_it << std::endl;
			return;
		}
	}
	covalent_bonds_.push_back( bi );
	changed();
}

/// @brief returns a list of valid cutpoints, in N-->C order
Cutpoints
StructureData::cutpoints() const
{
	Cutpoints cutpoints;
	for ( auto s=segments_begin(); s!=segments_end(); ++s ) {
		core::Size const cut = segment( *s ).cutpoint();
		if ( cut ) cutpoints.push_back( cut );
	}
	return cutpoints;
}

void
StructureData::set_cutpoint( core::Size const resid )
{
	std::string const seg_name = segment_name( resid );
	SegmentResid const segment_resid = segment( seg_name ).pose_to_segment( resid );
	segment_nonconst( seg_name ).set_cutpoint(  segment_resid );
}

/// @brief marks the resi-th residue of seg as a cutpoint
void
StructureData::set_cutpoint( std::string const & seg, SegmentResid const resi )
{
	segment_nonconst( seg ).set_cutpoint( resi );
	changed();
}

/// @brief connects the given chains together, doesn't update anything -- don't call this on its own unless you know what you're doing.
void
StructureData::connect_segments(
	std::string const & segment1_c,
	std::string const & segment2_n )
{
	// connected segments
	if ( has_free_upper_terminus( segment1_c ) && has_free_lower_terminus( segment2_n ) ) {
		mark_connected( segment1_c, segment2_n );
	}
	debug_assert( segment(segment1_c).upper_segment() == segment2_n );
	debug_assert( segment(segment2_n).lower_segment() == segment1_c );

	TR.Debug << "Marked " << segment1_c << " and " << segment2_n << " as connected." << std::endl;
	changed();
}

void
StructureData::mark_disconnected(
	std::string const & seg1,
	std::string const & seg2 )
{
	segment_nonconst( seg1 ).set_upper_segment( "" );
	segment_nonconst( seg2 ).set_lower_segment( "" );

	TR.Debug << "Marked " << seg1 << " and " << seg2 << " as disconnected." << std::endl;
	changed();
}

/// @brief disconnects the given chains doesn't update anything -- don't call this on its own unless you know what you're doing.
void
StructureData::disconnect_segments(
	std::string const & segment1_c,
	std::string const & segment2_n )
{
	mark_disconnected( segment1_c, segment2_n );
}

/// @brief merges two segments into one that has the name new_name. They must be next to each other in sequence.
void
StructureData::merge_segments(
	std::string const & segment1,
	std::string const & segment2,
	std::string const & new_name )
{
	Segment const & c1 = segment(segment1);
	Segment const & c2 = segment(segment2);

	SegmentNames c1_grp = segments_in_movable_group( c1.movable_group() );
	SegmentNames c2_grp = segments_in_movable_group( c2.movable_group() );
	debug_assert( ( c1_grp.size() <= 1 ) || ( c2_grp.size() <= 1 ) );

	core::Size new_group = c1.movable_group();
	if ( c2_grp.size() > c1_grp.size() ) {
		new_group = c2.movable_group();
	}

	// segments should be contiguous
	debug_assert( c1.upper()+1 == c2.lower() );

	std::string const newss = segment(segment1).ss() + segment(segment2).ss();
	std::string const new_abego = segment(segment1).abego() + segment(segment2).abego();
	debug_assert( newss.size() == c1.length() + c2.length() );
	debug_assert( new_abego.size() == c1.length() + c2.length() );

	// handle merging of residues segments
	std::string const low_seg = segment( segment1 ).lower_segment();
	std::string const upp_seg = segment( segment2 ).upper_segment();
	Segment resis( new_name, newss, new_abego, c1.nterm_included(), c2.cterm_included() );
	resis.set_movable_group( new_group );

	TR << "Removing old segments " << *segments_.find( segment1 ) << " and " << *segments_.find( segment2 )
		<< " and adding new merged segment " << resis << std::endl;
	SegmentNameList::const_iterator after_me = find_segment_name( segment2 );
	if ( after_me != segment_order_.end() ) ++after_me;
	std::string const after_str = *after_me;
	delete_segment( segment2 );
	core::Size const num_groups = movable_groups().size();
	delete_segment( segment1 );
	if ( movable_groups().size() < num_groups ) {
		resis.set_movable_group( movable_groups().size() + 1 );
	}
	add_segment( resis, find_segment_name( after_str ) );
	if ( !low_seg.empty() ) {
		connect_segments( low_seg, new_name );
	}
	if ( !upp_seg.empty() ) {
		connect_segments( new_name, upp_seg );
	}
}

/// @brief checks pose vs. StructureData info
/// @throw EXCN_BadInput if things don't match
void
StructureData::check_pose( core::pose::Pose const & pose ) const
{
	TR.Debug << "Checking pose properties..." << std::endl;
	if ( ! pose.size() ) return;

	// subtract only those virtual residues added to the end of the pose
	core::Size total_residue = pose.size();
	for ( core::Size resid=pose.size(); resid>=1; --resid ) {
		if ( pose.residue( resid ).aa() == core::chemical::aa_vrt ) --total_residue;
		else break;
	}

	if ( total_residue != pose_length() ) {
		std::stringstream msg;
		msg << id() << ": pose length does not match StructureData.  Pose length = "
			<< total_residue << " SD length = " << pose_length() << std::endl;
		msg << *this << std::endl;
		throw CREATE_EXCEPTION(EXCN_PoseInconsistent, msg.str() );
	}
	if ( total_residue != ss().size()  ) {
		std::stringstream msg;
		msg << id() << ": pose length does not match StructureData secstruct length. Pose length = "
			<< total_residue << " SS length = " << ss().size() << " SD: " << *this << std::endl;
		throw CREATE_EXCEPTION(EXCN_PoseInconsistent, msg.str() );
	}
}

/// @brief returns n-terminal residue of the chain represented by given string
core::Size
StructureData::lower_anchor( std::string const & id_val ) const
{
	return segment(id_val).start();
}

/// @brief returns c-terminal residue of the chain represented by given string
core::Size
StructureData::upper_anchor( std::string const & id_val ) const
{
	return segment(id_val).stop();
}

/// @brief returns non-const residue range of the segment represented by given string
Segment &
StructureData::segment_nonconst( std::string const & id_val )
{
	auto it = segments_.find( id_val );
	if ( it == segments_.end() ) {
		it = segments_.find( id() + PARENT_DELIMETER + id_val );
	}
	if ( it == segments_.end() ) {
		std::stringstream err;
		err << id() << ": Segment not found in residue lists! ";
		err << "Search term is: " << id_val << "; Segment map is: " << segments_ << std::endl;
		throw CREATE_EXCEPTION(utility::excn::Exception,  err.str() );
	}
	return it->second;
}

/// @brief returns residue range of the segment represented by given string
Segment const &
StructureData::segment( std::string const & id_val ) const
{
	auto it = segments_.find( id_val );
	if ( it == segments_.end() ) {
		it = segments_.find( id() + PARENT_DELIMETER + id_val );
	}
	if ( it == segments_.end() ) {
		std::stringstream err;
		err << id() << ": Segment not found in residue lists! ";
		err << "Search term is: " << id_val << "; Segment map is: " << segments_ << std::endl;
		throw CREATE_EXCEPTION(utility::excn::Exception,  err.str() );
	}
	return it->second;
}

bool
StructureData::has_data_int( std::string const & segment_id, std::string const & data_name ) const
{
	return has_data_int( segment_id + DATA_DELIMETER + data_name );
}

/// @brief check for real number data
bool
StructureData::has_data_real( std::string const & segment_id, std::string const & data_name ) const
{
	return has_data_real( segment_id + DATA_DELIMETER + data_name );
}

/// @brief gets real number data
bool
StructureData::has_data_str( std::string const & segment_id, std::string const & data_name ) const
{
	return has_data_str( segment_id + DATA_DELIMETER + data_name );
}

/// @brief sets integer number data
void
StructureData::set_data_int( std::string const & segment_id, std::string const & data_name, int const val )
{
	set_data_int( segment_id + DATA_DELIMETER + data_name, val );
}

/// @brief sets integer number data
void
StructureData::set_data_int( std::string const & data_name, int const val )
{
	auto it = data_int_.find(data_name);
	if ( it == data_int_.end() ) {
		data_int_[data_name] = val;
	} else {
		it->second = val;
	}
	changed();
}

/// @brief sets real number data
void
StructureData::set_data_real( std::string const & segment_id, std::string const & data_name, core::Real const val )
{
	set_data_real( segment_id + DATA_DELIMETER + data_name, val );
}

/// @brief sets real number data
void
StructureData::set_data_real( std::string const & data_name, core::Real const val )
{
	auto it = data_real_.find(data_name);
	if ( it == data_real_.end() ) {
		data_real_[data_name] = val;
	} else {
		it->second = val;
	}
	changed();
}

/// @brief sets string data
void
StructureData::set_data_str( std::string const & segment_id, std::string const & data_name, std::string const & val )
{
	return set_data_str( segment_id + DATA_DELIMETER + data_name, val );
}

/// @brief sets string data
void
StructureData::set_data_str( std::string const & data_name, std::string const & val )
{
	auto it = data_str_.find(data_name);
	if ( it == data_str_.end() ) {
		data_str_[data_name] = val;
	} else {
		it->second = val;
	}
	changed();
}

/// @brief gets int number data
int
StructureData::get_data_int( std::string const & segment_id, std::string const & data_name ) const
{
	return get_data_int( segment_id + DATA_DELIMETER + data_name );
}

/// @brief gets real number data
int
StructureData::get_data_int( std::string const & data_name ) const
{
	auto it = data_int_.find(data_name);
	if ( it == data_int_.end() ) {
		std::stringstream err( "In StructureData: " );
		err << id() << ": " << data_name << " not found in data map." << std::endl;
		for ( it = data_int_.begin(); it != data_int_.end(); ++it ) {
			err << it->first << " : " << it->second << std::endl;
		}
		throw CREATE_EXCEPTION(utility::excn::Exception,  err.str() );
	}
	return it->second;
}

/// @brief gets real number data
core::Real
StructureData::get_data_real( std::string const & segment_id, std::string const & data_name ) const
{
	return get_data_real( segment_id + DATA_DELIMETER + data_name );
}

/// @brief gets real number data
core::Real
StructureData::get_data_real( std::string const & data_name ) const
{
	auto it = data_real_.find(data_name);
	if ( it == data_real_.end() ) {
		std::stringstream err( "In StructureData: " );
		err << id() << ": " << data_name << " not found in data map." << std::endl;
		for ( it = data_real_.begin(); it != data_real_.end(); ++it ) {
			err << it->first << " : " << it->second << std::endl;
		}
		throw CREATE_EXCEPTION(utility::excn::Exception,  err.str() );
	}
	return it->second;
}

/// @brief gets string data
std::string const &
StructureData::get_data_str( std::string const & segment_id, std::string const & data_name ) const
{
	return get_data_str( segment_id + DATA_DELIMETER + data_name );
}

/// @brief gets string data
std::string const &
StructureData::get_data_str( std::string const & data_name ) const
{
	auto it = data_str_.find(data_name);
	if ( it == data_str_.end() ) {
		std::stringstream err( "In StructureData: " );
		err << id() << ": " << data_name << " not found in data map." << std::endl;
		for ( it = data_str_.begin(); it != data_str_.end(); ++it ) {
			err << it->first << " : " << it->second << std::endl;
		}
		throw CREATE_EXCEPTION(utility::excn::Exception,  err.str() );
	}
	return it->second;
}

/// @brief copies user data fields from one permutation to this one -- existing data is overridden
void
StructureData::copy_data( StructureData const & perm )
{
	copy_data( perm, true );
}

template< class T >
void
set_map_data(
	std::map< std::string, T > & map,
	std::string const & key,
	T const & value,
	bool const overwrite )
{
	if ( overwrite ) {
		map[ key ] = value;
	} else {
		if ( map.insert( std::make_pair( key, value ) ).second == false ) {
			TR.Debug << "skipping new map data with different value for key " << key
				<< " : " << map.find( key )->second << " exists vs. new value " << value << std::endl;
		}
	}
}

/// @brief copies user data fields from one permutation to this one -- gives the option to not overwrite data
void
StructureData::copy_data( StructureData const & perm, bool const overwrite )
{
	for ( const auto & d : perm.data_real() ) {
		set_map_data( data_real_, d.first, d.second, overwrite );
	}
	for ( const auto & d : perm.data_int() ) {
		set_map_data( data_int_, d.first, d.second, overwrite );
	}
	for ( const auto & d : perm.data_str() ) {
		set_map_data( data_str_, d.first, d.second, overwrite );
	}
	for ( const auto & a : perm.aliases() ) {
		set_map_data( aliases_, a.first, a.second, overwrite );
	}
	for ( auto bi=perm.covalent_bonds_begin(); bi!=perm.covalent_bonds_end(); ++bi ) {
		add_covalent_bond( *bi );
	}
	for ( auto p=perm.pairings_begin(); p!=perm.pairings_end(); ++p ) {
		// Add pairings where all segments are in the slice
		bool missing_a_segment = false;
		for ( const auto & s : (*p)->segments() ) {
			if ( !has_segment( s ) ) {
				missing_a_segment = true;
				break;
			}
		}
		TR.Debug << "Checking pairing segments " << (*p)->segments()
			<< " vs existing segments " << segment_order_ << " segment missing: "
			<< missing_a_segment << std::endl;
		if ( !missing_a_segment ) add_pairing( **p );
	}
	changed();
}

/// @brief returns segment which includes residue number res
std::string const &
StructureData::segment_name( core::Size const res ) const
{
	for ( const auto & segment : segments_ ) {
		if ( ( res >= segment.second.lower() ) && ( res <= segment.second.upper() ) ) {
			return segment.first;
		}
	}
	std::stringstream err;
	err << " Residue " << res << " was not found in the residues map!" << std::endl;
	throw CREATE_EXCEPTION(utility::excn::Exception,  err.str() );
}

/// @brief updates numbering based on the saved order of Segment objects
void
StructureData::update_numbering()
{
	movable_groups_.clear();
	debug_assert( segment_order_.size() == segments_.size() );
	core::Size cur_num = 1;
	core::Size non_dummy_count = 0;
	std::string new_ss = "";
	std::string new_abego = "";

	std::set< core::Size > mgs;
	for ( auto c=segments_begin(); c!=segments_end(); ++c ) {
		auto r = segments_.find( *c );
		debug_assert( r != segments_.end() );
		r->second.set_pose_start( cur_num );
		cur_num += r->second.length();
		non_dummy_count += ( r->second.stop() - r->second.start() + 1 );
		new_ss += r->second.ss();
		new_abego += r->second.abego();
		mgs.insert( r->second.movable_group() );
	}
	pose_length_= cur_num-1;
	length_ = non_dummy_count;
	if ( new_ss.size() != pose_length_ ) {
		std::stringstream err;
		err << id() << ": StructureData pose size doesn't match secondary structure string size.  this is probably an internal bug that needs to be fixed." << std::endl;
		err << *this << std::endl;
		err << "new ss= " << new_ss << std::endl;
		throw CREATE_EXCEPTION(utility::excn::BadInput,  err.str() );
	}
	debug_assert( new_ss.size() == pose_length_ );
	ss_ = new_ss;
	debug_assert( new_abego.size() == pose_length_ );
	abego_ = new_abego;
	TR.Debug << "Numbering updated - new pose length = " << pose_length_ << " new ss = " << new_ss << std::endl;
	movable_groups_ = MovableGroups( mgs.begin(), mgs.end() );
}

/// @brief attaches a template pose to the given segment
void
StructureData::set_template_pose(
	std::string const & seg_id,
	core::pose::Pose const & template_pose,
	core::Size const start_resid,
	core::Size const stop_resid )
{
	auto r = segments_.find( seg_id );
	if ( r == segments_.end() ) {
		std::stringstream msg;
		msg << "StructureData::set_template_pose(): Segment named " << seg_id
			<< " could not be found in the StructureData. Valid segments are: " <<
			segments_ << std::endl;
		utility_exit_with_message( msg.str() );
	}
	r->second.set_template_pose( template_pose, start_resid, stop_resid );
}

/// @brief overridden by derived classes if they need to do anything when the SD changes
void
StructureData::changed()
{
	if ( !block_signals_ ) {
		on_change();
	}
}

void
StructureData::on_change()
{
	TR.Debug << "StructureData::on_change()" << std::endl;
}

/// @brief renumbers movable group "oldg" to have new number "newg"
void
StructureData::renumber_movable_group( core::Size const oldg, core::Size const newg )
{
	for ( auto & segment : segments_ ) {
		if ( segment.second.movable_group() == oldg ) {
			segment.second.set_movable_group( newg );
		}
	}
	changed();
}

/// @brief returns true if this object has a group of segments with the given name
bool
StructureData::has_segment_group( std::string const & sname ) const
{
	if ( has_segment( sname ) ) return true;
	std::string const match_str = sname + PARENT_DELIMETER;
	for ( auto c=segments_begin(); c!=segments_end(); ++c ) {
		if ( boost::starts_with( *c, match_str ) ) {
			return true;
		}
	}
	return false;
}

/// @brief returns true if this object has a group of segments with the given name
SegmentNames
StructureData::segment_group( std::string const & sname ) const
{
	std::string const match_str = sname + PARENT_DELIMETER;
	SegmentNames retval;

	// check for perfect match first
	if ( has_segment( sname ) ) {
		retval.push_back( sname );
		return retval;
	}

	for ( auto c=segments_begin(); c!=segments_end(); ++c ) {
		if ( boost::starts_with( *c, match_str ) ) {
			retval.push_back( *c );
		}
	}
	return retval;
}

/// @brief checks consistency of the data
/// @throws EXCN_PoseInconsistent if there is a problem
void
StructureData::check_consistency() const
{
	check_residues();
	check_movable_groups();
}

void
StructureData::check_pose_consistency( core::pose::Pose const & pose ) const
{
	using core::pose::symmetry::extract_asymmetric_unit;

	check_residues();
	check_movable_groups();

	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		TR.Debug << "Extracting asymm unit " << std::endl;
		core::pose::Pose posecopy( pose );
		posecopy.remove_constraints();

		core::pose::Pose asymm_pose;
		extract_asymmetric_unit( posecopy, asymm_pose, false );
		check_pose( asymm_pose );
		check_chains( asymm_pose );
		check_improper_termini( asymm_pose );
	} else {
		check_pose( pose );
		check_chains( pose );
		check_improper_termini( pose );
	}
}

void
StructureData::check_residues() const
{
	TR.Debug << "Checking residues..." << std::endl;
	std::stringstream msg;
	utility::vector1< bool > accounted_for( pose_length(), false );
	for ( auto r=segments_.begin(); r!=segments_.end(); ++r ) {
		if ( r->second.upper() > pose_length() ) {
			msg << r->second << " has a upper terminal residue ("
				<< r->second.upper() << ") that is large than the pose length." << std::endl;
			throw CREATE_EXCEPTION(EXCN_PoseInconsistent, msg.str() );
		}

		for ( core::Size i=r->second.lower(); i<=r->second.upper(); ++i ) {
			if ( !accounted_for[i] ) {
				accounted_for[i] = true;
			} else {
				msg << r->second << " overlaps with something else at position " << i << std::endl << *this << std::endl;
				throw CREATE_EXCEPTION(EXCN_PoseInconsistent, msg.str() );
			}
		}
		if ( r->second.lower_segment() != "" ) {
			auto r2 = segments_.find( r->second.lower_segment() );
			if ( r2 == segments_.end() ) {
				msg << "Lower segment of " << r->second << " does not exist. SD=" << *this << std::endl;
				throw CREATE_EXCEPTION(EXCN_PoseInconsistent, msg.str() );
			}
		}
		if ( r->second.upper_segment() != "" ) {
			auto r2 = segments_.find( r->second.upper_segment() );
			if ( r2 == segments_.end() ) {
				msg << "Upper segment of " << r->second << " does not exist. SD=" << *this << std::endl;
				throw CREATE_EXCEPTION(EXCN_PoseInconsistent, msg.str() );
			}
		}

	}
	for ( core::Size i=1; i<=pose_length(); ++i ) {
		if ( !accounted_for[i] ) {
			msg << " Residue " << i << " is not accounted for by the permutation. " << *this << std::endl;
			throw CREATE_EXCEPTION(EXCN_PoseInconsistent, msg.str() );
		}
	}
}

void
StructureData::check_improper_termini( core::pose::Pose const & pose ) const
{
	TR.Debug << "Checking for improper termini..." << std::endl;
	for ( const auto & segment : segments_ ) {
		for ( core::Size i=segment.second.lower()+1; i<segment.second.upper(); ++i ) {
			if ( pose.residue(i).is_terminus() ) {
				std::stringstream msg;
				msg << " Residue " << i << " has a terminal variant but is inside segment " << segment.second << "." << std::endl;
				throw CREATE_EXCEPTION(EXCN_PoseInconsistent, msg.str() );
			}
		}
	}
}

utility::vector1< core::Size >
compute_cutpoints( core::pose::Pose const & pose )
{
	utility::vector1< core::Size > cutpoints;
	for ( core::Size resid=1; resid!=pose.size(); ++resid ) {
		if ( pose.residue( resid ).has_variant_type( core::chemical::CUTPOINT_LOWER ) ) {
			cutpoints.push_back( resid );
			if ( ( resid < pose.size() ) && !pose.residue( resid + 1 ).has_variant_type( core::chemical::CUTPOINT_UPPER ) ) {
				std::stringstream msg;
				msg << "StructureData.cc:compute_cutpoints(): Residue " << resid << " has a lower cutpoint variant, but residue "
					<< resid + 1 << " does not have an upper cutpoint variant." << std::endl;
				utility_exit_with_message( msg.str() );
			}
		}
	}
	return cutpoints;
}

utility::vector1< core::Size >
compute_chain_beginnings( core::pose::Pose const & pose )
{
	utility::vector1< core::Size > beginnings;
	utility::vector1< core::Size > endings = pose.conformation().chain_endings();
	beginnings.push_back( 1 );
	for ( utility::vector1< core::Size >::const_iterator r=endings.begin(); r!=endings.end(); ++r ) {
		beginnings.push_back( *r + 1 );
	}
	return beginnings;
}

void
advance_past_virtual_residues(
	core::pose::Pose const & pose,
	utility::vector1< core::Size >::const_iterator & cur_resid,
	utility::vector1< core::Size >::const_iterator const & end_resids )
{
	for ( ; cur_resid != end_resids; ++cur_resid ) {
		if ( pose.residue( *cur_resid ).aa() != core::chemical::aa_vrt ) break;
	}
}

void
StructureData::check_chains( core::pose::Pose const & pose ) const
{
	TR.Debug << "Checking chains..." << std::endl;

	// get ordered sets of chain endings in the pose
	utility::vector1< core::Size > chain_beginnings = compute_chain_beginnings( pose );
	utility::vector1< core::Size > chain_endings = pose.conformation().chain_endings();
	utility::vector1< core::Size > cutpoints = compute_cutpoints( pose );

	if ( pose.size() ) {
		chain_endings.push_back( pose.size() );
	}

	TR.Trace << "Beginnings: " << chain_beginnings << std::endl;
	TR.Trace << "Endings: " << chain_endings << std::endl;
	TR.Trace << "Cutpoints: " << cutpoints << std::endl;

	debug_assert( chain_beginnings.size() == chain_endings.size() );

	utility::vector1< core::Size >::const_iterator cur_beginning = chain_beginnings.begin();
	utility::vector1< core::Size >::const_iterator cur_ending = chain_endings.begin();
	utility::vector1< core::Size >::const_iterator cur_cutpoint = cutpoints.begin();
	for ( auto s=segments_begin(); s!=segments_end(); ++s ) {
		TR.Debug << "Checking for an ending in this segment: " << *s << std::endl;
		// Handle cases for this segment being a chain end:
		//  1. Cutpoints count as endings
		if ( ( cur_cutpoint != cutpoints.end() ) && ( segment( *s ).cutpoint() == *cur_cutpoint ) ) {
			++cur_cutpoint;
		} else if ( segment( *s ).cutpoint() ) {
			std::stringstream msg;
			msg << "StructureData::check_chains(): StructureData has a cutpoint at position "
				<< segment( *s ).cutpoint() << ", but no cutpoint was found in the Pose. Pose cutpoints: "
				<< cutpoints << std::endl;
			msg << *this << std::endl;
			throw CREATE_EXCEPTION(EXCN_PoseInconsistent, msg.str() );
		}

		// 2. Lower termini count as chain beginnings
		if ( ( cur_beginning != chain_beginnings.end() ) && ( segment( *s ).lower() == *cur_beginning ) ) {
			if ( segment( *s ).has_free_lower_terminus() ) {
				++cur_beginning;
				if ( ( cur_cutpoint != cutpoints.end() ) && ( *cur_cutpoint+1 == segment( *s ).lower() ) ) {
					++cur_cutpoint;
				}
			} else {
				std::stringstream msg;
				msg << "StructureData::check_chains(): StructureData has no chain beginning at position "
					<< segment( *s ).lower() << ", but that residue is a chain beginning in the pose. Pose chain beginnings: "
					<< chain_beginnings << std::endl;
				msg << *this << std::endl;
				throw CREATE_EXCEPTION(EXCN_PoseInconsistent, msg.str() );
			}
		}

		// 3. Upper termini count as chain endings
		if ( ( cur_ending != chain_endings.end() ) && ( segment( *s ).upper() == *cur_ending ) ) {
			if ( segment( *s ).has_free_upper_terminus() ) {
				++cur_ending;
				if ( ( cur_cutpoint != cutpoints.end() ) && ( *cur_cutpoint == segment( *s ).upper() ) ) {
					++cur_cutpoint;
				}
			} else {
				std::stringstream msg;
				msg << "StructureData::check_chains(): StructureData has no chain ending at position "
					<< segment( *s ).upper() << ", but that residue is a chain ending in the pose. Pose chain endings: "
					<< chain_endings << " current ending: " << *cur_ending << std::endl;
				msg << *this << std::endl;
				throw CREATE_EXCEPTION(EXCN_PoseInconsistent, msg.str() );
			}
		}
	}

	advance_past_virtual_residues( pose, cur_cutpoint, cutpoints.end() );
	advance_past_virtual_residues( pose, cur_beginning, chain_beginnings.end() );
	advance_past_virtual_residues( pose, cur_ending, chain_endings.end() );

	if ( cur_cutpoint != cutpoints.end() ) {
		std::stringstream msg;
		msg << "StructureData::check_chains(): Pose has a cutpoint at position "
			<< *cur_cutpoint << ", but no cutpoint was set in the StructureData.  Pose cutpoints: "
			<< cutpoints << std::endl;
		msg << *this << std::endl;
		throw CREATE_EXCEPTION(EXCN_PoseInconsistent, msg.str() );
	}

	if ( cur_beginning != chain_beginnings.end() ) {
		std::stringstream msg;
		msg << "StructureData::check_chains(): Pose has a chain beginning at position "
			<< *cur_beginning << ", but that residue is not a chain beginning in the StructureData. Pose chain beginnings: "
			<< chain_beginnings << std::endl;
		msg << *this << std::endl;
		throw CREATE_EXCEPTION(EXCN_PoseInconsistent, msg.str() );
	}

	if ( cur_ending != chain_endings.end() ) {
		std::stringstream msg;
		msg << "StructureData::check_chains(): Pose has a chain ending at position "
			<< *cur_ending << ", but that residue is not a chain ending in the StructureData. Pose chain endings: "
			<< chain_endings << std::endl;
		msg << *this << std::endl;
		throw CREATE_EXCEPTION(EXCN_PoseInconsistent, msg.str() );
	}

	TR.Debug << "Done checking chain endings " << std::endl;
}

void
StructureData::check_movable_groups() const
{
	// check movable groups
	TR.Debug << "Movable groups set is " << movable_groups_ << std::endl;
	for ( unsigned long movable_group : movable_groups_ ) {
		if ( movable_group <= 0 ) {
			std::stringstream msg;
			msg << " StructureData has a movable group <= 0 : " << *this << std::endl;
			throw CREATE_EXCEPTION(EXCN_PoseInconsistent, msg.str() );
		}
	}
}

/// @brief delete residue at resid
void
StructureData::delete_residue( core::Size const pose_resid )
{
	std::string const segname = segment_name( pose_resid );
	core::Size const seg_resid = segment( segname ).pose_to_segment( pose_resid );
	TR.Debug << "Deleting residue " << pose_resid << " from StructureData - " << segname << " " << seg_resid << std::endl;
	if ( seg_resid < 1 ) {
		segment_nonconst( segname ).delete_lower_padding();
	} else if ( seg_resid > segment( segname ).elem_length() ) {
		segment_nonconst( segname ).delete_upper_padding();
	} else if ( segment( segname ).elem_length() == 1 ) {
		delete_segment( segname );
	} else {
		segment_nonconst( segname ).delete_residue( seg_resid );
	}
	update_numbering();
	changed();
}

core::Size
StructureData::movable_group( core::Size const resid ) const
{
	return segment( segment_name( resid ) ).movable_group();
}

/// @brief sets movable group of a segment
void
StructureData::set_movable_group( std::string const & segid, core::Size const mg )
{
	auto s = segments_.find( segid );
	if ( s == segments_.end() ) {
		std::stringstream err;
		err << "Error in StructureData::set_movable_group( " << segid << ", " << mg << "):"
			<< "segment not found! perm=" << *this << std::endl;
		throw CREATE_EXCEPTION(utility::excn::Exception, err.str() );
	}
	s->second.set_movable_group( mg );
	changed();
}

/// @brief tells if the segment given has an available lower terminus
bool
StructureData::has_free_lower_terminus( std::string const & id_val ) const
{
	if ( segments_.find( id_val ) == segments_.end() ) {
		return false;
	}
	return segment( id_val ).has_free_lower_terminus();
}

/// @brief tells if the segment given has an available lower terminus
bool
StructureData::has_free_upper_terminus( std::string const & id_val ) const
{
	if ( segments_.find( id_val ) == segments_.end() ) {
		return false;
	}
	return segment( id_val ).has_free_upper_terminus();
}

SegmentPairingCOP
StructureData::pairing( SegmentNames const & segments ) const
{
	SegmentPairingCOP retval;
	for ( SegmentPairingCOP const & pair : pairings_ ) {
		if ( pair->segments() == segments ) {
			if ( retval ) {
				std::stringstream msg;
				msg << "StructureData::pairing(): More than one pairing was found containing segments "
					<< segments << std::endl;
				msg << "SD=" << *this << std::endl;
				utility_exit_with_message( msg.str() );
			}
			retval = pair;
		}
	}
	return retval;
}

/// @brief non-const iterator to end of segment names list
SegmentNameList::iterator
StructureData::segments_end_nonconst()
{
	return segment_order_.end();
}

void
StructureData::block_signals()
{
	block_signals_ = true;
}

void
StructureData::unblock_signals()
{
	block_signals_ = false;
}

/// output
std::ostream &
operator<<( std::ostream & os, StructureData const & perm )
{
	os << "<StructureData name=\"" << perm.id()
		<< "\" length=\"" << perm.length()
		<< "\" pose_length=\"" << perm.pose_length()
		<< "\" >" << std::endl;

	// residues
	for ( auto c=perm.segment_order_.begin(); c!=perm.segment_order_.end(); ++c ) {
		os << "\t" << perm.segment( *c ) << std::endl;
	}
	// int data
	std::map< std::string, int >::const_iterator dat;
	for ( dat = perm.data_int_.begin(); dat != perm.data_int_.end(); ++dat ) {
		os << "\t<Int name=\"" << dat->first << "\" value=\"" << dat->second << "\" />" << std::endl;
	}
	// real data
	std::map< std::string, core::Real >::const_iterator datr;
	for ( datr = perm.data_real_.begin(); datr != perm.data_real_.end(); ++datr ) {
		os << "\t<Real name=\"" << datr->first << "\" value=\"" << datr->second << "\" />" << std::endl;
	}
	// string data
	std::map< std::string, std::string >::const_iterator dats;
	for ( dats = perm.data_str_.begin(); dats != perm.data_str_.end(); ++dats ) {
		os << "\t<Str name=\"" << dats->first << "\" value=\"" << dats->second << "\" />" << std::endl;
	}
	// aliases
	for ( const auto & aliase : perm.aliases_ ) {
		os << "\t<Alias name=\"" << aliase.first << "\" segment=\"" << aliase.second.first
			<< "\" res=\"" << aliase.second.second << "\" />" << std::endl;
	}
	// covalent bonds
	for ( auto b=perm.covalent_bonds_begin(); b!=perm.covalent_bonds_end(); ++b ) {
		os << *b << std::endl;
	}
	// pairings
	for ( auto p=perm.pairings_begin(); p!=perm.pairings_end(); ++p ) {
		os << **p << std::endl;
	}
	os << "</StructureData>";
	return os;
}

/*
/// dump contents of residues map
std::ostream &
operator<<( std::ostream & os, SegmentMap const & resmap )
{
for ( SegmentMap::const_iterator r=resmap.begin(); r!=resmap.end(); ++r ) {
os << r->second;
}
return os;
}
*/

std::ostream &
operator<<( std::ostream & os, Alias const & alias )
{
	os << "segment=\"" << alias.first << "\" res=\"" << alias.second << "\"";
	return os;
}

//////////////////////////////////////////////////////////////
///  StructureDataCreator
//////////////////////////////////////////////////////////////
basic::datacache::WriteableCacheableDataOP
StructureDataCreator::create_data( std::istream & in ) const
{
	StructureDataFactory const & factory = *StructureDataFactory::get_instance();
	return utility::pointer::make_shared< StructureData >( factory.create_from_cacheable_data( in ) );
}

std::string
StructureDataCreator::keyname() const
{
	return StructureData::class_name();
}

} // namespace components
} // namespace denovo_design
} // namespace protocols

#ifdef    SERIALIZATION

/// @brief Serialization method for BondInfo
template< class Archive >
void
protocols::denovo_design::components::BondInfo::save( Archive & arc ) const
{
	arc( CEREAL_NVP( seg1 ) );
	arc( CEREAL_NVP( seg2 ) );
	arc( CEREAL_NVP( res1 ) );
	arc( CEREAL_NVP( res2 ) );
	arc( CEREAL_NVP( atom1 ) );
	arc( CEREAL_NVP( atom2 ) );
}

/// @brief Deserialization method for BondInfo
template< class Archive >
void
protocols::denovo_design::components::BondInfo::load( Archive & arc )
{
	arc( seg1 );
	arc( seg2 );
	arc( res1 );
	arc( res2 );
	arc( atom1 );
	arc( atom2 );
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::denovo_design::components::BondInfo );
CEREAL_REGISTER_TYPE( protocols::denovo_design::components::BondInfo )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_denovo_design_components_BondInfo )

/// @brief Serialization method for StructureData
/// @details This does not serialize either the length_event_ nor length_event_link_
/// data members as these hold raw pointers and would be invalid following deserialization.
template< class Archive >
void
protocols::denovo_design::components::StructureData::save( Archive & arc ) const
{
	arc( CEREAL_NVP( id_ ) );
	arc( CEREAL_NVP( ss_ ) );
	arc( CEREAL_NVP( abego_ ) );
	arc( CEREAL_NVP( pose_length_ ) );
	arc( CEREAL_NVP( length_ ) );
	arc( CEREAL_NVP( data_int_ ) );
	arc( CEREAL_NVP( data_real_ ) );
	arc( CEREAL_NVP( data_str_ ) );
	arc( CEREAL_NVP( segments_ ) );
	arc( CEREAL_NVP( aliases_ ) );
	arc( CEREAL_NVP( covalent_bonds_ ) );
	arc( CEREAL_NVP( segment_order_ ) );
	arc( CEREAL_NVP( movable_groups_ ) );
	arc( CEREAL_NVP( pairings_ ) );
	arc( CEREAL_NVP( block_signals_ ) );
}

/// @brief Deserialization method
/// @brief This does not deserialize either the length_event_ vector nor the length_event_link_
/// data member because both of these classes hold raw pointers and would be invalid following
/// deserialization.
template< class Archive >
void
protocols::denovo_design::components::StructureData::load( Archive & arc )
{
	arc( id_ );
	arc( ss_ );
	arc( abego_ );
	arc( pose_length_ );
	arc( length_ );
	arc( data_int_ );
	arc( data_real_ );
	arc( data_str_ );
	arc( segments_ );
	arc( aliases_ );
	arc( covalent_bonds_ );
	arc( segment_order_ );
	arc( movable_groups_ );
	arc( pairings_ );
	arc( block_signals_ );
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::denovo_design::components::StructureData );
CEREAL_REGISTER_TYPE( protocols::denovo_design::components::StructureData )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_denovo_design_components_StructureData )
#endif // SERIALIZATION
