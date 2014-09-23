// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/raw_data/RawFileData.hh
///
/// @brief FileData base class
/// @author James Thompson, Monica Berrondo

#ifndef INCLUDED_core_io_raw_data_RawFileData_hh
#define INCLUDED_core_io_raw_data_RawFileData_hh

// mini headers
#include <core/types.hh>

#include <core/io/raw_data/Raw.fwd.hh>
#include <core/io/raw_data/RawStruct.hh>

// AUTO-REMOVED #include <utility/file/file_sys_util.hh>

// C++ Headers
// AUTO-REMOVED #include <string>

#include <utility/vector1.hh>


namespace core {
namespace io {
namespace raw_data {
	/// @brief Abstract base class for classes that writes different types of
	/// silent-files that contain a mixture of Struct objects which are expected
	/// to be uniquely identified by some sort of string-based tag.
	class RawFileData {
	public:
		///////////////////////////////////////////////////////////////////////////
		// constructor

		RawFileData() { }

		/// @brief Destructor.
		virtual ~RawFileData() {
			clear_structure_map();
		}

		/// @brief Returns the number of structures contained in this container.
		int size() const { return structure_map_.size(); }

		/// @brief Returns the number of residues in the first structure in this object. Not
		/// guaranteed to be fixed for all structures in this container.
		int nres() const {
			return begin_const()->nres();
		}

		/// @brief Returns the sequence of the first structure in this object. Not
		/// guaranteed to be fixed for all structures in this container.
		std::string sequence() const {
			return sequence_;
		}

		/// @brief Remove all of the RawStruct objects from this object.
		void clear_structure_map()
		{
			structure_map_.clear();
		}

		/// @brief quickly read a list of tags from a silent-input file. Only checks lines beginning
		/// with SCORE: strings.
		utility::vector1< std::string > read_tags_fast( std::string const filename ) const;

		/// @brief write all RawStruct objects in the structure_map_ to the given filename.
		void write_all( const std::string filename, std::map < std::string, core::Real > const & score_map );

	protected:
		// mapping from tags to structure data pointers
		StructureMap structure_map_;
		std::string sequence_;

	public:
		class const_iterator;

		/// @brief Iterator class for RawFileData container.
		class iterator {

			friend class const_iterator;

		public:
			/// @brief empty constructor
			iterator() {}

			/// @brief Constructor, given an iterator into the StructureMap.
			iterator( StructureMap::iterator s_iter ) {
				it_ = s_iter;
			}

			~iterator() {}

			iterator& operator=( const iterator& src ) {
				it_ = src.it_;
				return (*this);
			}

			bool operator==( const iterator& other ) {
				return ( it_ == other.it_ );
			}

			bool operator!=( const iterator& other ) {
				return ( it_ != other.it_ );
			}

			iterator& operator++() {
				it_++;
				return (*this);
			}

			RawStructOP operator->() {
				return it_->second;
			}

			RawStructOP operator*() const {
				return it_->second;
			}

		private:
			StructureMap::iterator it_; // keep track of my place in a StructureMap
		}; // class iterator

		/// @brief const_iterator class for RawFileData container.
		class const_iterator {

		public:
			/// @brief empty constructor
			const_iterator() {}

			/// @brief Constructor, given an iterator into the StructureMap.
			const_iterator( StructureMap::const_iterator s_iter ) {
				it_ = s_iter;
			}

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
				it_++;
				return (*this);
			}

			RawStructOP operator->() const {
				return it_->second;
			}

			RawStructOP operator*() const {
				return it_->second;
			}

		private:
			StructureMap::const_iterator it_; // keep track of my place in a StructureMap
		}; // class iterator

		/// @brief Returns an iterator to the start of the members of this container.
		iterator begin() { return ( iterator( structure_map_.begin() ) ); }
		/// @brief Returns an iterator to the end of the members of this container.
		iterator end()   { return ( iterator( structure_map_.end()   ) ); }

		const_iterator begin_const() const { return ( const_iterator( structure_map_.begin() ) ); }
		const_iterator end_const()   const { return ( const_iterator( structure_map_.end()   ) ); }

	}; // class RawFileData
} // namespace silent
} // namespace io
} // namespace core

#endif
