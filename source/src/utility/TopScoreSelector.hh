// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/TopScoreSelector.hh
/// @brief  Class for keeping the N top-scoring objects
/// @author  Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_utility_TopScoreSelector_hh
#define INCLUDED_utility_TopScoreSelector_hh

// Unit headers
#include <utility/TopScoreSelector.fwd.hh>

// Package headers
#include <utility/vector1.hh>

// Project headers
#include <platform/types.hh>

// C++ headers
#include <utility/assert.hh>
#include <iostream>

namespace utility {

/// @brief The TopScoreSelector keeps an in-place, sorted, singly-linked list
/// of the N top-scoring objects and their scores that it encounters.  In the
/// case of ties, the first object(s) encountered are kept.
/// Class T must have a default constructor, a copy constructor, and must
/// be assignable. Class S must be assignable and comparable, both by a
/// < operator and by a > operator.
template < class T, class S >
class TopScoreSelector {
public:
	typedef 	platform::Size Size;

public:

	TopScoreSelector() :
		n_to_keep_( 0 ),
		low_is_better_( true ),
		worst_ind_( 0 ),
		is_sorted_( true )
	{}

	void n_to_keep( Size setting ) {
		n_to_keep_ = setting;
		top_scorers_.reserve( n_to_keep_ );
		next_ptr_.reserve( n_to_keep_ );
	}

	Size n_to_keep() const { return n_to_keep_; }

	void low_is_better( bool setting ) {
		low_is_better_ = setting;
	}

	bool low_is_better() const {
		return low_is_better_;
	}

	bool better( S s1, S s2 ) const {
		if ( low_is_better_ ) {
			return s1 < s2;
		} else {
			return s1 > s2;
		}
	}

	/// @brief Keep a particular element if it is better than the worst
	/// element, or if we have not yet stored n_to_keep_ elements.
	void insert( T const & obj, S score ) {
		if ( top_scorers_.size() < n_to_keep_ || better( score, top_scorers_[ worst_ind_ ].second  )) {
			add_to_set( obj, score );
		}
	}

	Size size() const {
		return top_scorers_.size();
	}

	T const &
	operator [] ( Size ind ) const {
		if ( ! is_sorted_ ) sort_top_scores();
		return sorted_top_scorers_[ ind ].first;
	}

	S
	top_score( Size ind ) const {
		if ( ! is_sorted_ ) sort_top_scores();
		return sorted_top_scorers_[ ind ].second;
	}

private:

	void add_to_set( T const & obj, S score ) {
		is_sorted_ = false;

		Size new_index( 0 );
		if ( top_scorers_.size() < n_to_keep_ ) {
			top_scorers_.push_back( std::make_pair( obj, score ) );
			next_ptr_.push_back( 0 );
			new_index = top_scorers_.size();

			if ( worst_ind_ == 0 ) {
				worst_ind_ = new_index;
				return;
			}

		} else {
			/// replace the previous worst entry as the new worst entry
			new_index = worst_ind_;
			top_scorers_[ worst_ind_ ].first = obj;
			top_scorers_[ worst_ind_ ].second = score;
			worst_ind_ = next_ptr_[ worst_ind_ ];
		}

		Size ind_next( worst_ind_ ), last( 0 );
		while ( true ) {
			if ( ind_next == 0 ) {
				next_ptr_[ last ] = new_index;
				next_ptr_[ new_index ] = 0;
				break;
			}
			//std::cout << "add_to_set: " << ind_next << " " << new_index << " " << score << " " << top_scorers_[ ind_next ].second << std::endl;
			if ( better( score, top_scorers_[ ind_next ].second ) ) {
				last = ind_next;
				ind_next = next_ptr_[ ind_next ];
			} else {
				next_ptr_[ new_index ] = ind_next;
				if ( ind_next == worst_ind_ ) {
					worst_ind_ = new_index;
				} else {
					next_ptr_[ last ] = new_index;
				}
				break;
			}
		}

		/*Size next = worst_ind_;
		while ( next != 0 ) {
			std::cout << " " << next << " " << top_scorers_[ next ].first << " " << top_scorers_[ next ].second << std::endl;
			next = next_ptr_[ next ];
		}*/

	}

	/// @brief Copy the top_scorers_ data into the sorted_top_scorers_ array.
	void sort_top_scores() const {
		is_sorted_ = true;
		sorted_top_scorers_.resize( top_scorers_.size() );
		Size count_down = top_scorers_.size();
		Size next_ind = worst_ind_;
		while ( count_down > 0 ) {
		debug_assert( next_ind != 0 );
			sorted_top_scorers_[ count_down ] = top_scorers_[ next_ind ];
			--count_down;
			next_ind = next_ptr_[ next_ind ];
		}
	}

private:

	Size n_to_keep_;
	bool low_is_better_;
	Size worst_ind_;
	vector1< std::pair< T, S > > top_scorers_;
	vector1< Size > next_ptr_;

	mutable bool is_sorted_;
	mutable vector1< std::pair< T, S > > sorted_top_scorers_;
};

}

#endif
