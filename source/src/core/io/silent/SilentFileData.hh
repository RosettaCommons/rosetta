// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/silent/SilentFileData.hh
///
/// @brief silent input file reader for mini
/// @author James Thompson

#ifndef INCLUDED_core_io_silent_SilentFileData_hh
#define INCLUDED_core_io_silent_SilentFileData_hh

// mini headers
#include <core/types.hh>

#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SharedSilentData.hh>
#include <core/pose/full_model_info/FullModelParameters.hh>

#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <string>
#include <map>

#include <utility/vector1.hh>


namespace core {
namespace io {
namespace silent {

/// @brief Abstract base class for classes that read and write different types of
/// silent-files. Silent-files can contain SilentStruct objects which are expected,
/// to be uniquely identified by some sort of string-based tag inside the file.
class SilentFileData : public utility::pointer::ReferenceCount {

private:
	///
	/// @brief mapping from tags to structure data pointers
	Structure_Map structure_map_;


	/// @brief mapping from index in a silent file to data pointers, for situations in which the ordering in the file matters.
	utility::vector1< SilentStructOP > structure_list_;

	utility::vector1< std::string > comment_lines_;

	mutable std::map< SharedSilentDataType, SharedSilentDataOP > shared_silent_data_;
	std::string filename_; // filename of the last file that we wrote a structure to
	bool store_argv_in_file_;
	bool strict_column_mode_;
	bool record_source_;
	std::string silent_struct_type_;
	bool verbose_;
	core::pose::full_model_info::FullModelParametersOP full_model_parameters_;

	SilentFileOptions options_;

public:
	///////////////////////////////////////////////////////////////////////////
	// constructor

	SilentFileData( SilentFileOptions const & options );
	SilentFileData( std::string const& filename, SilentFileOptions const & options );

	SilentFileData(
		const std::string &filename,
		bool  store_argv_in_file,
		bool  strict_column_mode,
		const std::string & silent_struct_type,
		SilentFileOptions const & options
	);

	/// @brief Read in the SilentStruct objects contained in the given filename.
	/// this version will throw an exception if things go wrong (boolean return value is thus always true)
	bool read_file(
		std::string const & filename
	);

	/// @brief Read in the SilentStruct objects contained in the given filename.
	/// this version returns with a boolean to tell you about success
	bool _read_file(
		std::string const & filename,
		bool throw_exception_on_bad_structs = false
	);


	/// @brief Read in the SilentStruct objects contained in the given filename.
	/// Ignore any SilentStruct with a tag not in the tags vector.
	/// throw an exception if things go wrong (returned boolean always true)
	bool read_file(
		std::string const & filename,
		utility::vector1< std::string > const & tags
	);

	bool read_stream(
		std::istream & data,
		utility::vector1< std::string > const & tags,
		bool throw_exception_on_bad_structs, /*default false*/
		std::string filename="read_from_stream" /** for error reporting **/
	);

	/// @brief Read in the SilentStruct objects contained in the given filename.
	/// Ignore any SilentStruct with a tag not in the tags vector.
	/// returns with a boolean to tell about success
	bool _read_file(
		std::string const & filename,
		utility::vector1< std::string > const & tags,
		bool throw_exception_on_bad_structs = false
	);


	/// @brief Returns true if this silent-file plans on storing
	/// option.get_argv() as a comment in the silent-file, false
	/// otherwise.
	bool store_argv_in_file() const {
		return store_argv_in_file_;
	}

	/// @brief Sets a boolean variable governing if this silent-file
	/// plans on storing option.get_argv() as a comment in the
	/// silent-file.
	void store_argv_in_file( bool new_setting ) {
		store_argv_in_file_ = new_setting;
	}


	/// @brief Write the given silent-struct to the given outfile.
	bool write_silent_struct(
		SilentStruct & s,
		std::string const & filename,
		bool bWriteScoreOnly = false
	) const;

	/// @brief Write the given silent-struct, including the header, to the given
	/// ostream.
	bool _write_silent_struct(
		SilentStruct & s,
		std::ostream & out,
		bool bWriteScoreOnly = false
	) const;

	/// @brief Write the given silent-struct, minus the header, to the given
	/// ostream.
	bool write_silent_struct(
		SilentStruct & s,
		std::ostream & out,
		bool bWriteScoreOnly = false
	) const;

	/// @brief Write a comment to the given output stream.
	void write_comment(
		std::ostream & out,
		std::string const & line
	) const;

	/// @brief Returns the number of structures contained in this container.
	Size size() const { return structure_map_.size(); }

	/// @brief Returns the number of residues in the first structure in
	/// this object. Not guaranteed to be fixed for all structures in
	/// this container.
	int nres() const;

	/// @brief Sets the filename that this SilentFileData object will
	/// write to.
	void set_filename( std::string filename ) {
		filename_ = filename;
	}

	/// @brief Sets the filename that this SilentFileData object will
	/// write to.
	void set_verbose( bool const setting ) { verbose_ = setting; }

	SilentStructOP operator[] (std::string tag);

	/// @brief Gets the filename that this SilentFileData object will
	/// write to.
	std::string const& filename() const {
		return filename_;
	}
	/// @brief Return all tags in this container.
	utility::vector1< std::string > tags() const;

	/// @brief quickly read a list of tags from a silent-input file. Only checks
	/// lines beginning with SCORE: strings.
	utility::vector1< std::string > read_tags_fast(
		std::string const & filename
	) const;

	/// @brief Function to access the vector of silent structure owning pointers
	/// ordered as they were in the input file.
	utility::vector1 <SilentStructOP> structure_list() { return structure_list_; }

	std::string
	get_sequence( std::string const & filename );

	std::pair< utility::vector1< int >, utility::vector1< char > >
	get_resnum( std::string const & filename );

	/// @brief quickly read a list of tags from a silent-input file. Only checks
	/// lines beginning with SCORE: strings.
	bool read_tags_fast(
		std::string const & filename, utility::vector1< std::string >&
	) const;

	/// @brief return mode=first,last,all matched tags -- currently matching 'expression*' to tags in file, boost::regexp possible later
	bool matched_tags(
		std::string const& expression,
		std::string const& mode,
		utility::vector1< std::string >& tags_in_file
		//  utility::vector1< SilentStructOP >& decoys_in_file,
		//  bool ignore_decoys = false
	) const;

	/// @brief Returns a boolean indicating whether or not the strict_column_mode
	/// is turned on when printing scores.
	/// @details If strict_column_mode() is true, then the first SilentStruct
	/// printed to this SilentFileData object sets the EnergyNames that will be
	/// printed for all other SilentStruct objects. Extra EnergyNames in
	/// subsequent SilentStruct objects are ignored. If new objects are missing
	/// energies that should be printed in strict_column_mode, missing energy
	/// values are set to zero. In !strict_column_mode(), when each SilentStruct
	/// is printed, a new SCORE: header is printed if energies differ from last
	/// printed SilentStruct.
	bool strict_column_mode() const {
		return strict_column_mode_;
	}

	/// @brief Sets value for strict column mode. See strict_column_mode() for
	/// more information.
	void strict_column_mode( bool new_mode ) {
		strict_column_mode_ = new_mode;
	}

	////////////////////////////////////////////////////
	void set_record_source( bool const & new_mode ) {
		record_source_ = new_mode;
	}

	/// @brief Returns true if we have a SilentStruct matching the given tag,
	/// false otherwise.
	bool has_tag( std::string const & tag ) const {
		return ( structure_map_.count( tag ) > 0 );
	}

	/// @brief Returns a vector1 of the comment lines read in from
	/// a silent-file. Comment lines are any lines in the silent-file
	/// beginning with a # character.
	utility::vector1< std::string > comment_lines() {
		return comment_lines_;
	}

	/// @brief Adds a comment-line that will be printed in write_all method.
	/// Comment lines are simply lines in the silent-file that begin with the #
	/// character, and are printed in the write_all method.
	void comment_line( std::string const & line ) {
		comment_lines_.push_back( line );
	}

	/// @brief Removes the worst ( 1 - score_fraction ) percent of the decoys by
	/// score. The value of score_fraction should be between 0 and 1.0.
	void
	score_filter(
		Real const score_fraction
	);

	void
	reverse_score_filter(
		Real const score_fraction
	);

	/// @brief Orders silent structs by energy.
	void
	order_by_energy();

	// This was basically the old add_structure...
	void add_structure_replace_tag_if_necessary( SilentStructOP & new_struct );

	/// @brief Adds a SilentStructOP to the structure_map_. If the SilentStruct's
	/// tag already exists in the structure_map_, a new tag is assigned. Careful
	/// with this method, as it stores an owning pointer. If you change the
	/// SilentStruct later, it will change your already stored structures.
	void add_structure( SilentStructOP const & new_struct );

	/// @brief In stepwise modeling, sometimes poses carry refs to 'other poses'
	bool add_as_other_struct_if_relevant( SilentStructOP const & new_struct );

	/// @brief push_back to provide compatibility with other std containers.
	void push_back( SilentStructOP const & new_struct ) {
		add_structure( new_struct );
	}

	/// @brief Saves a copy of the silent struct. This method is:
	/// - SAFE! in the sense that it actually copies the SilentStruct object, not
	/// just the pointer to the object.
	/// - SLOW! in the sense that copying the object takes a small amount of time.
	void add_structure( SilentStruct const & new_struct );

	/// @brief Return a SilentStruct referred to by the given tag. Assumes that
	/// we have checked the tag!!
	SilentStruct const &
	get_structure (
		const std::string & tag
	) const {
		return *( structure_map_.find(tag)->second );
	}

	/// @brief Remove all of the SilentStruct objects from this object.
	void clear_structure_map() {
		structure_map_.clear();
		structure_list_.clear();
	}

	/// @brief Clears all of the data associated with this object.
	void clear() {
		clear_structure_map();
		comment_lines_.clear();

		shared_silent_data_.clear();
		filename_.clear();
		store_argv_in_file_ = false;
		strict_column_mode_ = false;
	}

	/// @brief Clears all of the data associated with this object.
	void clear_shared_silent_data() {
		shared_silent_data_.clear();
	}

	/// @brief Destructor.
	virtual ~SilentFileData();

	/// @brief write all SilentStruct objects in the structure_map_ to the given
	/// filename.
	void write_all(
		std::string const & filename, bool bWriteScoreOnly = false
	) const;

	/// @brief renumber all of the decoys in this SilentFileData object. This
	/// obliterates decoy tag associated with every SilentStruct object, and tries
	/// to sensibly number the decoys with similar and increasingly numbered decoy
	/// tags.
	void renumber_all_decoys();

	SilentStructOP create_SilentStructOP();

	/// @brief SharedSilentData methods
	SharedSilentDataOP get_shared_silent_data( SharedSilentDataType ssdt ) const;
	void set_shared_silent_data(
		SharedSilentDataType ssdt, SharedSilentDataOP ssd_op
	) const;
	bool has_shared_silent_data( SharedSilentDataType ssdt ) const;

private:
	//some utility function for the reading process --- returns bool if new type
	bool read_silent_struct_type_from_remark( std::string const& line, bool const header = false /*make true if this is one of the first 3 lines*/ );

	bool read_full_model_parameters_from_remark( std::string const& line, bool const header = false );

	bool check_if_rna_from_sequence_line( std::string const& sequence_line );

public:
	/// @brief Iterator class for SilentFileData container.
	class iterator {

		friend class const_iterator;

	public:
		typedef SilentStructOP  value_type;
		typedef SilentStructOP* pointer;
		typedef SilentStructOP& reference;
		typedef std::ptrdiff_t  difference_type;
		typedef std::forward_iterator_tag iterator_category;
		/// @brief empty constructor
		iterator() {}

		/// @brief Constructor, given an iterator into the Structure_Map.
		iterator( Structure_Map::iterator s_iter ) :
			it_( s_iter)
		{}

		~iterator() {}

		iterator& operator=( const iterator& src ) {
			it_ = src.it_;
			return (*this);
		}

		bool operator==( const iterator& other ) const {
			return ( it_ == other.it_ );
		}

		bool operator!=( const iterator& other ) const {
			return ( it_ != other.it_ );
		}

		iterator& operator++() {
			++it_;
			return (*this);
		}

		iterator& operator--() {
			--it_;
			return (*this);
		}

		SilentStructOP operator->() const {
			return it_->second;
		}

		SilentStructOP operator*() const {
			return it_->second;
		}

	protected:
		Structure_Map::iterator it_; // keep track of my place in a Structure_Map
	}; // class iterator

	/// @brief const_iterator class for SilentFileData container.
	class const_iterator {
	public:
		typedef SilentStructOP  value_type;
		typedef SilentStructOP* pointer;
		typedef SilentStructOP& reference;
		typedef std::ptrdiff_t  difference_type;
		typedef std::bidirectional_iterator_tag iterator_category;

		/// @brief empty constructor
		const_iterator() {}

		/// @brief Constructor, given an iterator into the Structure_Map.
		const_iterator( Structure_Map::const_iterator s_iter ) :
			it_( s_iter )
		{}

		~const_iterator() {}

		const_iterator& operator=( const const_iterator& src ) {
			it_ = src.it_;
			return (*this);
		}

		bool operator==( const const_iterator& other ) {
			return ( it_ == other.it_ );
		}

		bool operator!=( const const_iterator& other ) {
			return ( it_ != other.it_ );
		}

		const_iterator& operator++() {
			++it_;
			return (*this);
		}

		const_iterator& operator--() {
			--it_;
			return (*this);
		}

		SilentStructOP operator->() const {
			return it_->second;
		}

		SilentStructOP operator*() const {
			return it_->second;
		}

	private:
		Structure_Map::const_iterator it_; // keeps track of my place in a Structure_Map
	}; // class iterator

	//  void open_for_writing( utility::io::ozstream&, std::string const& filename, std::stringstream& ) const; //open silent-file and write header if first

	/// @brief Returns an iterator to the start of the members of this container.
	///
	iterator begin() { return ( iterator( structure_map_.begin() ) ); }

	/// @brief Returns an iterator to the start of the members of this container.
	///
	const_iterator begin() const { return ( const_iterator( structure_map_.begin() ) ); }

	/// @brief Returns an iterator to the end of the members of this container.
	///
	iterator end()   { return ( iterator( structure_map_.end()   ) ); }

	/// @brief Returns an iterator to the end of the members of this container.
	///
	const_iterator end() const { return ( const_iterator( structure_map_.end()   ) ); }

	iterator
	get_iterator_for_tag (
		const std::string & tag
	) {
		return iterator( structure_map_.find(tag) );
	}

	const_iterator
	get_iterator_for_tag (
		const std::string & tag
	) const {
		return const_iterator( structure_map_.find(tag) );
	}

	//const_iterator begin_const() const { return ( const_iterator( structure_map_.begin() ) ); }
	// const_iterator end_const()   const { return ( const_iterator( structure_map_.end()   ) ); }

}; // class SilentFileData

} // namespace silent
} // namespace io
} // namespace core

#endif
