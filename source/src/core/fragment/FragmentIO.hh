// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief  some utilities for fragments
/// @author Oliver Lange (olange@u.washington.edu)
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_core_fragment_FragmentIO_HH
#define INCLUDED_core_fragment_FragmentIO_HH

// Project Headers
#include <core/types.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/Frame.fwd.hh>
#include <core/fragment/FrameList.fwd.hh>
#include <core/fragment/SingleResidueFragData.fwd.hh>

// Utility headers
#include <utility/vector1.fwd.hh>

// C/C++
#include <map>

//Auto Headers
namespace core {
namespace fragment {

class FragFactory {
	typedef std::map< std::string, SingleResidueFragDataOP > SRFD_Types;
	typedef std::map< std::string, FrameOP > FrameTypes;

public:
	FragFactory(void);

	void add_frag_type(std::string const& type_name,
		SingleResidueFragDataOP frag_type);

	void add_frame_type( std::string const& type_name, FrameOP new_frag );

	FrameOP frame( std::string const& frame_name ) const;
	SingleResidueFragDataOP frag_type( std::string const& frag_name ) const;

private:
	SRFD_Types frag_types_;
	FrameTypes frame_types_;
}; // class FragFactory


class FragmentIO {
	typedef std::map< std::string, FragSetOP > FragFileCache;

public:
	FragmentIO() : top_( 0 ), ncopies_( 1 ), bAnnotate_( true ) {};
	FragmentIO( Size top, Size ncopies = 1, bool bAnnotate = true ) :
		top_( top ),
		ncopies_( ncopies ),
		bAnnotate_( bAnnotate )
	{};

	/// @brief read a FragSet... note that this function caches the fragment set.
	/// i.e., if you read the same set from multiple positions in the code you get
	/// the same set. if you read a new file ... we'll check the cache for stale
	/// references and delete them...
	FragSetOP read_data( std::string const & filename );

	void read_data( std::string const& filename, FrameList& );

	void read_data( std::istream& data, FrameList& );

	void write_data( std::string const& file, FragSet const& frags );

	FragFactory& get_frag_factory();

	/// @brief Updates the number of distinct fragments to keep
	void set_top_frag_num( Size setting ) {
		top_ = setting;
	}

	/// @brief Updates the number of copies of each fragment to keep (default 1)
	void set_ncopies( Size setting ) {
		ncopies_ = setting;
	}

	/// @brief Toggles between reading annotated (true) and non-annotated (false)
	// fragment files
	void set_read_annotation( bool setting = true ) {
		bAnnotate_ = setting;
	}

	/// @brief remove all FragSets that are not referenced outside the cache.
	void clean_frag_cache();

private:
	void read_next_frames(std::istream& data,
		std::string& next_line,
		FrameList &next_frames);

	void read_frag_data(std::istream& data,
		std::string& next_line,
		FrameList &next_frames );

	static FragFactory frag_factory_;
	static FragFileCache frag_cache_;

	// Number of distinct fragments to keep
	Size top_;

	// Number of times each of the <top_> distinct fragments is included in the
	// result. At the conclusion of file I/O, there will be an equal number of
	// distinct fragments. The total size of the list is <top_> * <copies_>.
	Size ncopies_;

	// Toggles between reading fragment files in annotated format (includes PDB id
	// and chain) and non-annotated format
	bool bAnnotate_;
};

} // fragment
} // core

#endif
