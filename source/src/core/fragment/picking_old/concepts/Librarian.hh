// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragment/picking_old/concepts/Librarian.hh
/// @brief  Librarian template for sorting through and extracting desired fragments
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_fragment_picking_old_concepts_Librarian_hh
#define INCLUDED_core_fragment_picking_old_concepts_Librarian_hh

// unit headers
#include <core/fragment/picking_old/concepts/Librarian.fwd.hh>

#include <core/fragment/picking_old/vall/gen/VallFragmentGen.hh>
#include <core/fragment/picking_old/vall/eval/VallFragmentEval.hh>
#include <core/fragment/picking_old/vall/scores/VallFragmentScore.hh>
#include <core/fragment/picking_old/vall/VallLibrary.hh>

// type headers
#include <core/types.hh>

// utility headers

// C++ headers
#include <algorithm>
#include <functional>

#include <utility/vector1.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>

namespace core {
namespace fragment {
namespace picking_old {
namespace concepts {


/// @brief Librarian template for sorting through and extracting desired fragments
template< typename Bookmark, typename ExtentEvaluator, typename ExtentGenerator, typename Library >
class Librarian {


public: // typedefs


	typedef core::Size Size;

	typedef typename Library::BookIterator BookIterator;
	typedef typename Library::BookConstIterator BookConstIterator;
	typedef typename Library::Book Book;

	typedef typename Book::PageConstIterator PageConstIterator;
	typedef typename Book::PageIterator PageIterator;
	typedef typename Book::Page Page;

	typedef typename ExtentGenerator::Extent Extent;

	typedef utility::pointer::shared_ptr< ExtentEvaluator > ExtentEvalOP;
	typedef utility::pointer::shared_ptr< ExtentEvaluator const > ExtentEvalCOP;
	typedef utility::pointer::shared_ptr< ExtentGenerator > ExtentGenOP;
	typedef utility::pointer::shared_ptr< ExtentGenerator const > ExtentGenCOP;

	typedef utility::vector1< Bookmark > Bookmarks;
	typedef typename Bookmarks::const_iterator BookmarkConstIterator;
	typedef typename Bookmarks::iterator BookmarkIterator;


protected: // typedefs


	typedef utility::vector1< ExtentGenOP > ExtentGenOPs;
	typedef utility::vector1< ExtentEvalOP > ExtentEvalOPs;


public: // construct/destruct


	/// @brief default constructor
	inline
	Librarian() {}


	/// @brief default destructor
	inline
	virtual
	~Librarian() {}


private: // disallow copy


	/// @brief disallow copy constructor
	// NOTE: If implementing copy in the future, remember to clone OPs.
	Librarian( Librarian const & rval );


	/// @brief disallow copy assignment
	// NOTE: If implementing copy in the future, remember to clone OPs.
	Librarian & operator =( Librarian const & rval );


public: // library operations


	/// @brief create sorted list corresponding to fragments in Library
	/// @details uses Bookmark < for evaluation
	/// @return true if creation successful, false otherwise (e.g. no ExtentEvaluators or ExtentGenerators found)
	inline
	bool catalog( Library const & library ) {
		return catalog( library, std::less< Bookmark >() );
	}


	/// @brief create sorted list corresponding to fragments in Library
	/// @tparam LessThan predicate <tt> Pr( left, right ) </tt> evaluating <tt> left < right <tt> for Bookmarks
	/// @return true if creation successful, false otherwise (e.g. no ExtentEvaluators or ExtentGenerators found)
	template< typename LessThan >
	bool catalog( Library const & library, LessThan const & lt ) {
		using std::push_heap;
		using std::sort_heap;

		if ( egen_.size() == 0 || eeval_.size() == 0 ) {
			return false;
		}

		bookmarks_.clear();

		for ( BookConstIterator book = library.begin(), end_of_library = library.end(); book != end_of_library; ++book ) {
			for ( PageConstIterator page = book->begin(), end_of_book = book->end(); page != end_of_book; ++page ) {

				// run through extents generated from this page
				for ( typename ExtentGenOPs::const_iterator e = egen_.begin(), ee = egen_.end(); e != ee; ++e ) {
					Extent const extent = (**e)( page, end_of_book );

					if ( extent.valid ) {
						Bookmark mark;
						if ( evaluate_extent( extent, mark ) ) {
							bookmarks_.push_back( mark );
							push_heap( bookmarks_.begin(), bookmarks_.end(), lt ); // maintain heap property
						}
					}
				} // foreach extent

			} // foreach page
		} // foreach book

		// sort heap with ordering '<'
		sort_heap( bookmarks_.begin(), bookmarks_.end(), lt );

		return true;
	}


public: // extent management


	/// @brief add an extent generator
	inline
	void add_extent_gen( ExtentGenCOP const & gen ) {
		egen_.push_back( gen->clone() );
	}


	/// @brief clear list of generators
	inline
	void clear_extent_gen() {
		egen_.clear();
	}


public: // evaluator management


	/// @brief add extent evaluator
	inline
	void add_extent_eval( ExtentEvalCOP val ) {
		eeval_.push_back( val->clone() );
	}


	/// @brief clear list of evaluators
	inline
	void clear_extent_eval() {
		eeval_.clear();
	}


protected: // methods


	/// @brief evaluate a fragment starting from Page at iterator
	/// @return Bookmark containing scores for fragment
	/// @remarks at least one helper needs to score the fragment extent
	inline
	bool evaluate_extent( Extent const & extent, Bookmark & mark ) {
		bool extent_allowed = true;

		for ( typename ExtentEvalOPs::iterator i = eeval_.begin(), ie = eeval_.end(); extent_allowed && i != ie; ++i ) {
			extent_allowed = extent_allowed && (**i)( extent, mark );
		}

		return extent_allowed;
	}


	/// @brief get the current bookmark heap
	inline
	Bookmarks const & bookmarks() const {
		return bookmarks_;
	}


	/// @brief get the current bookmark heap
	inline
	Bookmarks & bookmarks() {
		return bookmarks_;
	}


	/// @brief the list of extent generators
	inline
	ExtentGenOPs const & extent_gen() const {
		return egen_;
	}


	/// @brief the list of extent generators
	inline
	ExtentGenOPs & extent_gen() {
		return egen_;
	}


	/// @brief the list of extent evaluators
	inline
	ExtentEvalOPs const & extent_eval() const {
		return eeval_;
	}


	/// @brief the list of extent evaluators
	inline
	ExtentEvalOPs & extent_eval() {
		return eeval_;
	}


private: // data


	/// @brief heap of bookmarks
	/// @brief after catalog(), will be sorted w/respect to operator < of Bookmark
	Bookmarks bookmarks_;


	/// @brief generators for page extents
	ExtentGenOPs egen_;


	/// @brief list of page extent evaluators
	ExtentEvalOPs eeval_;


};

class Librarian_VallFragmentScore_VallFragmentEval_VallFragmentGen_VallLibrary: public
	Librarian<core::fragment::picking_old::vall::scores::VallFragmentScore,
	core::fragment::picking_old::vall::eval::VallFragmentEval,
	core::fragment::picking_old::vall::gen::VallFragmentGen,
	core::fragment::picking_old::vall::VallLibrary> {};

} // concepts
} // picking_old
} // fragment
} // core


#endif /* INCLUDED_core_fragment_picking_old_concepts_Librarian_HH */
