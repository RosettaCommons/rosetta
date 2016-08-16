// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/integer_mapping.hh
/// @brief  A set of useful classes to map between two enumerations.  So far, only a subset mapping is implemented.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_utility_integer_mapping_hh
#define INCLUDED_utility_integer_mapping_hh

#include <utility/integer_mapping.fwd.hh>

#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <platform/types.hh>

namespace utility {

/// @brief This class handles the bookeeping to map between a set of
/// integer ids in the "source" enumeration to a subset of those
/// ids -- the destination enumartion.  Elements in the source enumeration
/// that do not map to elements in the destination enumeration
/// are represented by the value 0.  Both enumerations should
/// count from 1.  Once the class has been initialized, this class
/// offers O(1) mapping between elements in the enumerations.
class subset_mapping : public utility::pointer::ReferenceCount {
public:
	typedef utility::pointer::ReferenceCount parent;
	static platform::Size const UNMAPPED;

public:
	subset_mapping();
	subset_mapping( platform::Size source_enumeration_size );
	subset_mapping( subset_mapping const & src );
	subset_mapping & operator = ( subset_mapping const & rhs );

	virtual ~subset_mapping();

	/// @brief Required before the first call to set_next_correspondence may be called.
	/// The size of the source enumeration must be known before the mapping may begin
	void set_source_size( platform::Size );

	/// @brief If you know the size of the destination enumeration, then you can save
	/// some under-the-hood vector resizing operations by informing the subset_mapping
	/// its size up front.  This call must proceed the first call to set_next_correspondence
	void reserve_destination_size( platform::Size );

	/// @brief Inform the mapping of the next source-enumeration id that should be mapped
	/// to a destination-enumeration id.  This will increase the size of the destination
	/// enumeration by one.  It is not essential that the source-enumeration ids appear
	/// in sorted order, however, by construction, the destination ids will be in sorted order.
	void set_next_correspondence( platform::Size source_id );

	/// @brief The number of elements in the source enumeration
	platform::Size source_size() const;

	/// @brief The number of elements in the destination enumeration -- this
	/// represents the number of calls that have been made to set_next_correspondence()
	platform::Size destination_size() const;

	/// @brief Map from the id of an element in  source enumeration to an id in the
	/// the destination enumeration, which may in fact be UNMAPPED.
	platform::Size s2d( platform::Size source_id ) const;

	/// @brief Map from the id of an element in the destination enumeration to an id
	/// in the source enumeration.  This is guaranteed to return a non-zero value.
	platform::Size d2s( platform::Size destination_id ) const;

	bool source_id_is_mapped( platform::Size source_id ) const;

private:
	utility::vector1< platform::Size > src_2_dst_;
	utility::vector1< platform::Size > dst_2_src_;

};

}

#endif
