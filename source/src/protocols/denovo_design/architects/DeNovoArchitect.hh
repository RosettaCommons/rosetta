// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/denovo_design/architects/StructureArchitect.hh
/// @brief Designs topologies
/// @author Tom Linsky (tlinsky@uw.edu)
/// @note   This is interface: it has no fields, and only
///         pure virtual methods.  No further constructors should
///         be defined.

#ifndef INCLUDED_protocols_denovo_design_architects_DeNovoArchitect_hh
#define INCLUDED_protocols_denovo_design_architects_DeNovoArchitect_hh

// Unit headers
#include <protocols/denovo_design/architects/DeNovoArchitect.fwd.hh>
#include <protocols/denovo_design/architects/StructureArchitect.hh>

// Protocol headers
#include <protocols/denovo_design/components/Segment.fwd.hh>
#include <protocols/denovo_design/components/StructureData.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.fwd.hh>

// Boost headers
#include <boost/lexical_cast.hpp>

// C++ headers
#include <set>

namespace protocols {
namespace denovo_design {
namespace architects {

/// @brief for planning ideal pieces of structures
/// @details Derived classes must still implment their own parse_tag() and type() functions.
///          This base handles the apply() virtual, though. In parse_tag(), derived classes
///          MUST be sure to set the motif list.
class DeNovoArchitect : public StructureArchitect {
public:
	DeNovoArchitect( std::string const & id );

	virtual
	~DeNovoArchitect();

	// pure virtual API
public:
	virtual std::string
	type() const = 0;

	virtual DeNovoArchitectOP
	clone() const = 0;

	virtual components::StructureDataOP
	design( core::pose::Pose const & pose, core::Real & random ) const = 0;

protected:
	virtual void
	parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data ) = 0;

public:
	components::StructureDataOP
	apply( core::pose::Pose const & pose ) const;

}; // DeNovoArchitect

/// @brief for planning arbitrary motifs
class DeNovoMotifArchitect : public DeNovoArchitect {
public:
	typedef components::Segment Motif;
	typedef components::SegmentOP MotifOP;
	typedef components::SegmentCOP MotifCOP;
	typedef utility::vector1< MotifCOP > MotifCOPs;

public:
	DeNovoMotifArchitect( std::string const & id );

	virtual ~DeNovoMotifArchitect();

	virtual DeNovoArchitectOP
	clone() const;

	static std::string
	architect_name() { return "DeNovoMotif"; }

	virtual std::string
	type() const;

	virtual components::StructureDataOP
	design( core::pose::Pose const & pose, core::Real & random ) const;

protected:
	virtual void
	parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data );

public:
	MotifCOPs::const_iterator
	motifs_begin() const;

	MotifCOPs::const_iterator
	motifs_end() const;

	/*u
	void
	set_lengths( std::string const & length_str );

	void
	set_lengths( Lengths const & lengths );
	*/

	void
	set_motifs( std::string const & motif_str );

	void
	set_motifs( MotifCOPs const & motifs );

private:
	MotifCOPs motifs_;
};

template< class T >
utility::vector1< T >
parse_length_str( std::string const & len_str )
{
	std::set< T > retval;
	utility::vector1< std::string > const str_residues( utility::string_split( len_str , ',' ) );
	for ( core::Size i = 1; i <= str_residues.size(); ++i ) {
		if ( str_residues[i] == "" ) continue;
		utility::vector1< std::string > const ranges( utility::string_split( str_residues[i], ':' ) );
		if ( ranges.size() == 1 ) {
			retval.insert( boost::lexical_cast< T >( ranges[1] ) );
		} else if ( ranges.size() == 2 ) {
			core::Size const start( boost::lexical_cast< T >( ranges[1] ) );
			core::Size const end( boost::lexical_cast< T >( ranges[2] ) );
			for ( core::Size i=start; i<=end; ++i ) {
				retval.insert( i );
			}
		} else {
			utility_exit_with_message( "Invalid length input: " + len_str );
		}
	}
	return Lengths( retval.begin(), retval.end() );
}

SecStructInfo
generate_secstruct_for_length(
	char const ss_char,
	std::string const & abego,
	core::Size const len );

} //protocols
} //denovo_design
} //architects

#endif //INCLUDED_protocols_denovo_design_architects_DeNovoArchitect_hh
