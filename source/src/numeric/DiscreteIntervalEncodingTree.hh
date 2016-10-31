// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/DiscreteIntervalEncodingTree.hh
/// @brief  definition of a discrete interval encoding tree ("DIET") as described by Martin Erwig
///         https://web.engr.oregonstate.edu/~erwig/papers/Diet_JFP98.pdf
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_numeric_DiscreteIntervalEncodingTree_HH
#define INCLUDED_numeric_DiscreteIntervalEncodingTree_HH

// Unit headers
// #include <numeric/DiscreteIntervalEncodingTree.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>

// C++ headers
#include <list>

namespace numeric {

template < class T >
class DietNode : public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< DietNode< T > >
{
public:
	typedef T Value;
	typedef typename utility::pointer::shared_ptr< DietNode< Value > > DietNodeOP;

	struct split_tuple {
		split_tuple( Value ext, DietNodeOP par ) : extrema( ext ), parent( par ) {}

		Value extrema;
		DietNodeOP parent;
	};

	typedef typename std::list< std::pair< Value, Value > > RangeList;

public:

	DietNode( Value value ) : lower_( value ), upper_( value ) {}

	//DietNode( Value lower, Value upper, DietNodeOP left, DietNodeOP right );

	/// @brief Return whether a particular value is contained by any of the ranges
	/// in the subtree rooted at this node.
	bool member( Value val ) const {
		if ( val < lower_ ) {
			if ( left_ ) return left_->member( val );
			return false;
		}
		if ( val > upper_ ) {
			if ( right_ ) return right_->member( val );
			return false;
		}
		// val is somewhere in the inclusive range between lower_ and upper_
		return true;
	}

	/// @brief Return the tuple of the largest value of the left child and the
	/// node that is the parent of that largest value.
	split_tuple
	split_max( DietNodeOP parent )
	{
		if ( ! right_ ) {
			return split_tuple( upper_, parent );
		} else {
			return right_->split_max( this->shared_from_this() );
		}
	}

	/// @brief Return the tuple of the smallest value of the right child and the
	/// node that is the parent of that smallest value.
	split_tuple
	split_min( DietNodeOP parent )
	{
		if ( ! left_ ) {
			return split_tuple( lower_, parent );
		} else {
			return left_->split_min( this->shared_from_this() );
		}
	}


	/// @brief Used in insert( val ) where val + 1 == lower_; this may reorder
	/// nodes in the tree.
	void
	join_left( Value new_lower ) {
		assert( new_lower + 1 == lower_ );
		if ( ! left_ ) {
			lower_ = new_lower;
		} else {
			split_tuple tup = left_->split_max( this->shared_from_this() );
			if ( tup.extrema + 1 == new_lower ) {
				if ( tup.parent.get() == this ) {
					lower_ = left_->lower_;
					left_ = left_->left_;
				} else {
					// remove the pointer to tup.parent's right child
					// and merge this child into the range covered by this node
					lower_ = tup.parent->right_->lower_;
					tup.parent->right_ = 0;
				}
			} else {
				lower_ = new_lower;
			}
		}
	}


	void
	join_right( Value new_upper )
	{
		assert( new_upper == upper_ + 1 );
		if ( ! right_ ) {
			upper_ = new_upper;
		} else {
			split_tuple tup = right_->split_min( this->shared_from_this() );
			if ( tup.extrema == new_upper + 1 ) {
				if ( tup.parent.get() == this ) {
					// merge the right node into this node; the right node then
					// disappears
					upper_ = right_->upper_;
					right_ = right_->right_;
				} else {
					// remove the pointer to tup.parent's left child
					// and merge this child into the range covered by this node
					upper_ = tup.parent->left_->upper_;
					tup.parent->left_ = 0;
				}
			} else {
				upper_ = new_upper;
			}
		}
	}

	void insert( Value new_val )
	{
		if ( new_val < lower_ ) {
			if ( new_val + 1 == lower_ ) {
				join_left( new_val );
			} else if ( ! left_ ) {
				left_ = DietNodeOP( new DietNode( new_val ) );
			} else {
				left_->insert( new_val );
			}
		}

		if ( upper_ < new_val ) {
			if ( upper_ + 1 == new_val ) {
				join_right( new_val );
			} else if ( ! right_ ) {
				right_ = DietNodeOP( new DietNode( new_val ) );
			} else {
				right_->insert( new_val );
			}
		}

		// else new_val is inside the range between lower_ and upper_
	}

	int
	size() const {
		return 1 + ( left_ ? left_->size() : 0 ) + ( right_ ? right_->size() : 0 );
	}

	void
	inorder_range_list( RangeList & ranges ) const {
		if ( left_ ) left_->inorder_range_list( ranges );
		ranges.push_back( std::make_pair( lower_, upper_ ) );
		if ( right_ ) right_->inorder_range_list( ranges );
	}

	bool
	correct( RangeList & ranges ) const {
		if ( left_ ) {
			if ( ! left_->correct( ranges ) ) {
				return false;
			}
		}
		if ( ! ranges.empty() ) {
			if ( ranges.back().second + 1 >= lower_ ) {
				return false;
			}
		}
		if ( lower_ > upper_ ) return false;
		ranges.push_back( std::make_pair( lower_, upper_ ) );
		if ( right_ ) {
			return right_->correct( ranges );
		}
		return true;
	}

private:

	Value lower_;
	Value upper_;

	DietNodeOP left_;
	DietNodeOP right_;

};

template < class T >
class DiscreteIntervalEncodingTree
{
public:
	typedef T Value;
	typedef typename utility::pointer::shared_ptr< DietNode< Value > > DietNodeOP;
	typedef typename std::list< std::pair< Value, Value > > RangeList;

public:

	bool
	member( Value val ) const {
		if ( ! root_ ) return false;
		return root_->member( val );
	}

	void
	insert( Value val ) {
		if ( ! root_ ) root_ = DietNodeOP( new DietNode< Value >( val ) );
		root_->insert( val );
	}

	int
	size() const {
		return root_ ? root_->size() : 0;
	}

	RangeList
	ranges() const {
		RangeList ranges;
		if ( root_ ) { root_->inorder_range_list( ranges ); }
		return ranges;
	}

	bool
	correct() const {
		if ( root_ ) {
			RangeList ranges;
			return root_->correct( ranges );
		}
		return true;
	}

private:
	DietNodeOP root_;

};

} // namespace numeric

#endif

