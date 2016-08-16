// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragment/picking_old/vall/eval/IdentityEval.hh
/// @brief  scores a fragment based on secondary structure identity and sequence identity
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_fragment_picking_old_vall_eval_IdentityEval_hh
#define INCLUDED_core_fragment_picking_old_vall_eval_IdentityEval_hh


// unit headers
#include <core/fragment/picking_old/vall/eval/IdentityEval.fwd.hh>

// type headers
#include <core/types.hh>

// package headers
#include <core/fragment/picking_old/vall/eval/VallFragmentEval.hh>
#include <core/fragment/picking_old/vall/VallLibrary.fwd.hh>

#include <utility/vector1.hh>


// C++ headers


namespace core {
namespace fragment {
namespace picking_old {
namespace vall {
namespace eval {


/// @brief scores a fragment based on sum of secondary structure identity and sequence identity
class IdentityEval : public VallFragmentEval {


private: // typedefs


	typedef VallFragmentEval Super;


public: // typedefs


	typedef Super::VallFragmentScore VallFragmentScore;

	typedef std::string String;
	typedef core::Real Real;


public: // concept typedefs


	/// @brief typedef for ExtentEvaluator concept
	typedef Super::PageConstIterator PageConstIterator;


	/// @brief typedef for ExtentEvaluator concept
	typedef Super::Extent Extent;


public: // concept translation typedefs


	typedef PageConstIterator VallResidueConstIterator;


public: // construct/destruct


	/// @brief default constructor
	IdentityEval();


	/// @brief full values constructor
	/// @param ss secondary structure string to match against
	/// @param aa amino acid structure string to match against
	IdentityEval(
		String const & ss,
		String const & aa,
		Real const ss_penalty = 1.0,
		Real const aa_penalty = 1.0,
		bool const randomize = true
	);


	/// @brief secondary structure constructor
	IdentityEval(
		String const & ss,
		Real const ss_penalty = 1.0,
		bool const randomize = true
	);


	/// @brief default copy constructor
	IdentityEval( IdentityEval const & rval );


	/// @brief default destructor
	virtual
	~IdentityEval();


public: // copy assignment


	/// @brief copy assignment
	IdentityEval & operator =( IdentityEval const & rval );


public: // virtual constructors


	/// @brief clone this object
	virtual
	VallFragmentEvalOP clone() const;


public: // virtual evaluation methods


	/// @brief for a fragment extent, evaluate and store results in a VallFragmentScore
	/// @return true, so score is always stored during VallLibrarian::catalog()
	virtual
	bool eval_impl(
		Extent const & extent,
		VallFragmentScore & fs
	);


public: // accessor/mutators


	/// @brief get secondary structure string
	inline
	String const & ss_str() const {
		return ss_;
	}


	/// @brief set secondary structure string
	inline
	void ss_str( String const & ss ) {
		ss_ = ss;
	}


	/// @brief get amino acid string
	inline
	String const & aa_str() const {
		return aa_;
	}


	/// @brief set amino acid string
	inline
	void aa_str( String const & aa ) {
		aa_ = aa;
	}


	/// @brief get secondary structure penalty
	inline
	Real ss_penalty() const {
		return ss_penalty_;
	}


	/// @brief set secondary structure penalty
	inline
	void ss_penalty( Real const penalty ) {
		ss_penalty_ = penalty;
	}


	/// @brief get amino acid penalty
	inline
	Real aa_penalty() const {
		return aa_penalty_;
	}


	/// @brief set amino acid penalty
	inline
	void aa_penalty( Real const penalty ) {
		aa_penalty_ = penalty;
	}


	/// @brief adding random noise to score?
	inline
	bool randomize() const {
		return randomize_;
	}


	/// @brief set flag to add random noise between [0, 0.001) to score
	inline
	void randomize( bool const flag ) {
		randomize_ = flag;
	}


public: // additional hooks


	/// @brief operation to be perform before catalog() starts
	virtual
	void pre_catalog_op( VallLibrary const & );


private: // data


	/// @brief secondary structure string to match against
	String ss_;


	/// @brief amino acid string to match against
	String aa_;


	/// @brief secondary structure penalty if non-matching
	Real ss_penalty_;


	/// @brief amino acid penalty if non-matching
	Real aa_penalty_;


	/// @brief flag to add random noise between [0, 0.001) into score
	bool randomize_;


};


} // namespace eval
} // namespace vall
} // namespace picking_old
} // namespace fragment
} // namespace core


#endif /* INCLUDED_core_fragment_picking_old_vall_eval_IdentityEval_HH */
