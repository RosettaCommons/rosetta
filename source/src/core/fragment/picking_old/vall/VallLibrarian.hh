// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragment/picking_old/vall/VallLibrarian.hh
/// @brief  Librarian that picks fragments from the Vall
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_fragment_picking_old_vall_VallLibrarian_hh
#define INCLUDED_core_fragment_picking_old_vall_VallLibrarian_hh


// unit headers
#include <core/fragment/picking_old/vall/VallLibrarian.fwd.hh>

// type headers
#include <core/types.hh>

// package headers
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/FragData.fwd.hh>
#include <core/fragment/picking_old/concepts/Librarian.hh>
#include <core/fragment/picking_old/vall/VallLibrary.hh>
#include <core/fragment/picking_old/vall/VallResidue.hh>
#include <core/fragment/picking_old/vall/eval/VallFragmentEval.hh>
#include <core/fragment/picking_old/vall/gen/VallFragmentGen.hh>

// project headers
#include <basic/Tracer.hh>

// utility headers
#include <utility/pointer/owning_ptr.hh>

// C++ headers
#include <functional>

#include <utility/vector1.hh>

#ifdef WIN32
#include <ctime>
#endif


namespace core {
namespace fragment {
namespace picking_old {
namespace vall {


/// @brief Librarian that picks fragments from the Vall
class VallLibrarian : public core::fragment::picking_old::concepts::Librarian< scores::VallFragmentScore, eval::VallFragmentEval, gen::VallFragmentGen, VallLibrary > {


private: // typedefs


	typedef core::fragment::picking_old::concepts::Librarian< scores::VallFragmentScore, eval::VallFragmentEval, gen::VallFragmentGen, VallLibrary > Super;


public: // typedefs


	typedef core::Size Size;

	typedef core::fragment::FragDataOPs FragDataOPs;
	typedef core::fragment::FragData FragData;
	typedef core::fragment::FragDataOP FragDataOP;

	typedef eval::VallFragmentEval VallFragmentEval;
	typedef scores::VallFragmentScore VallFragmentScore;


public: // concept typedefs


	typedef Super::BookIterator BookIterator;
	typedef Super::BookConstIterator BookConstIterator;
	typedef Super::PageConstIterator PageConstIterator;
	typedef Super::PageIterator PageIterator;
	typedef Super::ExtentEvalOP ExtentEvalOP;
	typedef Super::ExtentEvalCOP ExtentEvalCOP;
	typedef Super::ExtentGenOP ExtentGenOP;
	typedef Super::ExtentGenCOP ExtentGenCOP;
	typedef Super::Bookmarks Bookmarks;
	typedef Super::BookmarkConstIterator BookmarkConstIterator;
	typedef Super::BookmarkIterator BookmarkIterator;

	typedef Super::ExtentGenOPs ExtentGenOPs;
	typedef Super::ExtentEvalOPs ExtentEvalOPs;


public: // concept translation typedefs


	typedef ExtentEvalOP VallFragmentEvalOP;
	typedef ExtentEvalCOP VallFragmentEvalCOP;
	typedef ExtentGenOP VallFragmentGenOP;
	typedef ExtentGenCOP VallFragmentGenCOP;
	typedef Bookmarks Scores;
	typedef BookmarkConstIterator ScoreConstIterator;
	typedef BookmarkIterator ScoreIterator;


protected: // concept translation typedefs


	typedef BookConstIterator VallSectionConstIterator;
	typedef BookIterator VallSectionIterator;
	typedef PageConstIterator VallResidueConstIterator;
	typedef PageIterator VallResidueIterator;
	typedef ExtentGenOPs VallFragmentGenOPs;
	typedef ExtentEvalOPs VallFragmentEvalOPs;


public: // construct/destruct


	/// @brief default constructor
	VallLibrarian();


	/// @brief default destructor
	virtual
	~VallLibrarian();


private: // disallow copy


	/// @brief disallow copy constructor
	// NOTE: If implementing copy in the future, remember to clone OPs.
	VallLibrarian( VallLibrarian const & rval );


	/// @brief disallow copy assignment
	// NOTE: If implementing copy in the future, remember to clone OPs.
	VallLibrarian & operator =( VallLibrarian const & rval );


public: // accessors


	/// @brief preallocate scores container prior to catalog() to attempt
	///  speedup?, default true
	inline
	bool preallocate() const {
		return preallocate_;
	}


public: // mutators


	/// @brief set flag to preallocate scores container prior to catalog() to
	///  attempt speedup
	inline
	void preallocate( bool const flag ) {
		preallocate_ = flag;
	}


public: // generator management


	/// @brief add a fragment generator (aka extent generator)
	inline
	void add_fragment_gen( VallFragmentGenCOP const & gen ) {
		Super::add_extent_gen( gen );
	}


	/// @brief clear list of generators
	inline
	void clear_fragment_gen() {
		Super::clear_extent_gen();
	}


	/// @brief the number of currently defined fragment generators
	inline
	Size n_fragment_gen() const {
		return Super::extent_gen().size();
	}


public: // evaluator management


	/// @brief add a fragment evaluator (aka extent evaluator)
	inline
	void add_fragment_eval( VallFragmentEvalCOP eval ) {
		Super::add_extent_eval( eval );
	}


	/// @brief clear list of evaluators
	inline
	void clear_fragment_eval() {
		Super::clear_extent_eval();
	}


	/// @brief the number of currently defined fragment evaluators
	inline
	Size n_fragment_eval() const {
		return Super::extent_eval().size();
	}


public: // library operations


	/// @brief create sorted list corresponding to fragments in Library
	/// @details uses Score's '<' for evaluation
	/// @return true if creation successful, false otherwise (e.g. no VallFragmentEval or VallFragmentGen found)
	bool catalog( VallLibrary const & library ) {
		return catalog( library, std::less< VallFragmentScore >() );
	}


	/// @brief create sorted list corresponding to fragments in Library
	/// @tparam LessThan predicate <tt> Pr( left, right ) </tt> evaluating <tt> left < right <tt> for Scores (aka Bookmarks)
	/// @return true if creation successful, false otherwise (e.g. no VallFragmentEval or VallFragmentGen found)
	template< typename LessThan >
	bool catalog( VallLibrary const & library, LessThan const & lt ) {
		basic::Tracer TR( "core.fragment.picking_old.vall.VallLibrarian" );
		pre_catalog_ops( library );

		if ( TR.Debug.visible() ) {
			TR.Debug << "Cataloging " << library.size() << " residues in fragment library..." << std::endl;
		}

		time_t time_start = time( NULL );
		bool const status = Super::catalog( library, lt );
		time_t time_end = time( NULL );

		if ( TR.Debug.visible() ) {
			TR.Debug << "... done.  " << scores().size() << " scores filed.  Time elapsed: " << ( time_end - time_start ) << " seconds." << std::endl;
		}

		post_catalog_ops( library );

		return status;
	}


	/// @brief number of scores currently filed
	inline
	Size n_scores() const {
		return scores().size();
	}


public: // fragment extraction


	/// @brief get top 'N' fragments from prior catalog()
	/// @param n The number of fragments to get.
	/// @param srfd_type The BBTorsionSRFD type to use.
	FragDataOPs top_fragments(
		Size const n,
		BBTorsionSRFD const & srfd_type = BBTorsionSRFD()
	) const;


	/// @brief get fragments from prior catalog() [from, to]
	/// @param from index of the first fragment in the list, indexing starts from '1'
	/// @param to index of the last fragment in the list (inclusive)
	/// @param srfd_type The BBTorsionSRFD type to use.
	/// @return filled FragDataOPs if sort() was called successfully, otherwise empty FragDataOPs
	FragDataOPs fragments(
		Size from,
		Size to,
		BBTorsionSRFD const & srfd_type = BBTorsionSRFD()
	) const;


public: // concept translation


	/// @brief return scores container
	inline
	Scores const & scores() const {
		return Super::bookmarks();
	}


	/// @brief return scores container
	inline
	Scores & scores() {
		return Super::bookmarks();
	}


private: // additional catalog methods

	/// @brief this function runs before main routine in catalog() starts
	void pre_catalog_ops( VallLibrary const & library );

	/// @brief this function runs after main routine catalog() finishes
	void post_catalog_ops( VallLibrary const & library );

private: // data


	/// @brief flag controls preallocation of score container, default true
	bool preallocate_;

};


} // namespace vall
} // namespace picking_old
} // namespace fragment
} // namespace core


#endif /* INCLUDED_core_fragment_picking_old_vall_VallLibrarian_HH */
