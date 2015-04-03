// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Oliver Lange

#ifndef INCLUDED_protocols_topology_broker_SequenceNumberResolver_hh
#define INCLUDED_protocols_topology_broker_SequenceNumberResolver_hh

// Unit Headers
#include <protocols/topology_broker/SequenceNumberResolver.fwd.hh>

// Package Headers
#include <utility/pointer/ReferenceCount.hh>
#include <protocols/topology_broker/Exceptions.hh>


// Project Headers
#include <core/types.hh>

// C/C++ headers
#include <map>
#include <iterator>
#include <stdexcept>

namespace protocols {
namespace topology_broker {

class SequenceNumberResolver : public utility::pointer::ReferenceCount {

public:
	SequenceNumberResolver() {}

	SequenceNumberResolver( const SequenceNumberResolver& );

	//@throws EXCN_BadInput if label or offset is already taken
	void register_label_offset( std::string const& label, core::Size offset );

	//Returns global position of element at position resid in sequence claim with the specified label.
	core::Size find_global_pose_number( std::string const& label, core::Size resid) const;

	core::Size find_global_pose_number( std::pair< std::string, core::Size> pos_pair ) const;


	//Returns global position of first element of sequence claim with the specified label.
	core::Size find_global_pose_number( std::string const& label) const;

	//Returns label of sequence claim that corresponds to the global sequence position <pose_number>
	std::string find_label( core::Size pose_number ) const;

	//Returns position of global <pose_number> in corresponding sequence claim.
	core::Size find_local_pose_number( core::Size pose_number ) const;

	//Returns offset of a given label.
	core::Size offset( std::string const& label ) const {
		OffsetMap::const_iterator p = offset_map_.find(label);

		if( p != offset_map_.end() ) return p->second;
		else throw utility::excn::EXCN_BadInput( "SequenceNumberResolver asked to resolve SequenceClaim label '"
												 + label + "', which does not match any SequenceClaim labels." );
	}

	// returns the map element with the largest offset = SequenceClaim at the end of the sequence
	std::pair <core::Size, std::string> terminal_pair() const;


private:

	// map sequence labels to offsets
	typedef std::map< std::string, core::Size > OffsetMap;
	OffsetMap offset_map_;

	// Inverse map, mapping offsets to labels
	std::map< core::Size, std::string > offset_map_reversed_;

	//Searches for entry in reversed map with highest offset that is smaller than the given pose_number.
	//This entry contains label and offset corresponding to the given pose_number.
	std::map<core::Size, std::string>::const_iterator search_reversed_map( core::Size pose_number ) const;

}; //class SequenceNumberResolver

}
}

#endif  // INCLUDED_protocols_topology_broker_SequenceNumberResolver_hh
