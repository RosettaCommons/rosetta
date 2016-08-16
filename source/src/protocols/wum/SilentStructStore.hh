// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/wum/SilentStructStore.hh
/// @brief
/// @author Mike Tyka


#ifndef INCLUDED_protocols_wum_SilentStructStore_hh
#define INCLUDED_protocols_wum_SilentStructStore_hh

#include <protocols/wum/SilentStructStore.fwd.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/io/silent/SilentStruct.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <string>
#include <vector>

//Auto Headers
namespace protocols {
namespace wum {

/// @brief this little class is a predicate for finding silent structures in say a vector of silent structures
class find_SilentStructOPs
{
public:
	find_SilentStructOPs(std::string field, core::Real value ): field_(field), value_(value) {}

	bool operator () (const core::io::silent::SilentStructOP& check);

private:
	std::string field_;
	core::Real value_;
};

class SilentStructStore : public utility::pointer::ReferenceCount  {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~SilentStructStore();
	SilentStructStore()
	{

	}

public:

	// some wrapper stuff to allow iteration of the structure store
	typedef std::vector < core::io::silent::SilentStructOP >::iterator iterator;
	typedef std::vector < core::io::silent::SilentStructOP >::const_iterator const_iterator;

	/// @brief Returns an iterator to the start of the members of this container.
	iterator begin() { return store_.begin(); }

	/// @brief Returns an iterator to the start of the members of this container.
	const_iterator begin() const { return store_.begin(); }

	/// @brief Returns an iterator to the end of the members of this container.
	iterator end()   { return store_.end(); }

	/// @brief Returns an iterator to the end of the members of this container.
	const_iterator end() const { return store_.end(); }

	void sort_by( std::string field = "score" );

public:

	/// @brief Remove all structures
	void clear();

	/// @brief add a pose
	void add( const core::pose::Pose &pose );

	/// @brief Add a core::io::silent::SilentStruct
	void add( core::io::silent::SilentStructOP new_struct );

	/// @brief Add a core::io::silent::SilentStruct
	void add( const core::io::silent::SilentStruct &new_struct );

	/// @brief Add structures froma silent file data object
	void add( core::io::silent::SilentFileData const& sfd ) ;

	/// @brief Add the contents of another SilentStructStore
	void add( SilentStructStore &mergestore ) ;

	/// @brief THis uses the pose stream to read in everything from -l, -s and -in:file:silent into this store.
	void read_from_cmd_line( );

	/// @brief read from string
	void read_from_string( const std::string & input );

	/// @brief read from silent file
	void read_from_stream( std::istream & input );

	/// @brief read from string
	void read_from_file( const std::string & filename );

	/// @brief Obtain a new pose from a given index. must provide a template pose though!
	void get_pose( core::Size index,  core::pose::Pose &pose ) const ;

	/// @brief How many structures
	core::Size size(){ return store_.size(); }

	/// @brief Get a structure with a certain index
	core::io::silent::SilentStructCOP get_struct( core::Size index ) const {
		if ( index >= store_.size() ) runtime_assert( false );
		return (core::io::silent::SilentStructCOP) store_[ index ];
	}

	/// @brief Get a structure with a certain index
	core::io::silent::SilentStructOP& get_struct( core::Size index )  {
		if ( index >= store_.size() ) runtime_assert( false );
		return store_[ index ];
	}

	/// @brief Get a random structure
	core::io::silent::SilentStructCOP get_struct_random() const ;

	/// @brief Print silent file
	void serialize( std::ostream & out ) const ;

	/// @brief Print silent file
	void serialize( std::string & out ) const ;

	/// @brief Print silent file
	void serialize_to_file( const std::string &file ) const ;

	/// @brief Print silent file
	void print(  std::ostream & out ) const ;

	std::vector < core::io::silent::SilentStructOP > &  store() { return store_; }

	void limit( core::Size limit_size ){
		if ( store().size() > limit_size ) store().resize( limit_size );
	}

	/// @brief Return memory usage
	virtual core::Size mem_footprint() const;

	/// @brief return numner of structures
	core::Size size() const { return store_.size(); }

	/// @brief return numner of structures
	void erase( iterator it ) { store_.erase(it); }

	/// Manipulators:
	void all_add_energy( std::string scorename, core::Real value, core::Real weight = 1.0 );

	/// Manipulators:
	void all_sort_silent_scores( );


private:
	std::vector < core::io::silent::SilentStructOP > store_;
};


// Some additional tools!

std::string generate_unique_structure_id();


}
}

#endif


