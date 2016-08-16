// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragment/picking_old/vall/eval/ABEGOEval.hh
/// @brief  scores a fragment based on ABEGO
/// @author Nobuyasu Koga (nobuyasu@u.washington.edu)

#ifndef INCLUDED_core_fragment_picking_old_vall_eval_ABEGOEval_hh
#define INCLUDED_core_fragment_picking_old_vall_eval_ABEGOEval_hh


// unit headers
#include <core/fragment/picking_old/vall/eval/ABEGOEval.fwd.hh>
#include <core/sequence/ABEGOManager.fwd.hh>

// type headers
#include <core/types.hh>

// package headers
#include <core/fragment/picking_old/vall/eval/VallFragmentEval.hh>
#include <core/fragment/picking_old/vall/VallLibrary.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace fragment {
namespace picking_old {
namespace vall {
namespace eval {


/// @brief scores a fragment based on sum of secondary structure identity and sequence identity
class ABEGOEval : public VallFragmentEval {


private: // typedefs


	typedef VallFragmentEval Super;


public: // typedefs


	typedef Super::VallFragmentScore VallFragmentScore;

	typedef std::string String;
	typedef core::Real Real;

	typedef core::sequence::ABEGOManager ABEGOManager;
	typedef core::sequence::ABEGOManagerOP ABEGOManagerOP;

public: // concept typedefs


	/// @brief typedef for ExtentEvaluator concept
	typedef Super::PageConstIterator PageConstIterator;


	/// @brief typedef for ExtentEvaluator concept
	typedef Super::Extent Extent;


public: // concept translation typedefs


	typedef PageConstIterator VallResidueConstIterator;


public: // construct/destruct


	/// @brief default constructor
	ABEGOEval();


	/// @brief full values constructor
	/// @param ss secondary structure string to match against
	/// @param aa amino acid structure string to match against
	ABEGOEval(
		utility::vector1< String > const & input,
		Real const penalty = 1.0,
		bool const randomize = true
	);


	/// @brief default copy constructor
	ABEGOEval( ABEGOEval const & rval );


	/// @brief default destructor
	virtual
	~ABEGOEval();


public: // copy assignment


	/// @brief copy assignment
	ABEGOEval & operator =( ABEGOEval const & rval );


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
	utility::vector1< String > const & abego_str() const {
		return abego_;
	}

	/// @brief set abego vector of string
	inline
	void abego_str( utility::vector1< String > const & abego ) {
		abego_ = abego;
	}

	/// @brief get secondary structure penalty
	inline
	Real penalty() const {
		return penalty_;
	}

	/// @brief adding random noise to score?
	inline
	bool randomize() const {
		return randomize_;
	}


public: // additional hooks


	/// @brief operation to be perform before catalog() starts
	virtual
	void pre_catalog_op( VallLibrary const & );


private: // data


	/// @brief abego string to match against
	utility::vector1< String > abego_;

	/// @brief abego penalty if non-matching
	Real penalty_;

	/// @brief flag to add random noise between [0, 0.001) into score
	bool randomize_;

	ABEGOManagerOP am_;


};


} // namespace eval
} // namespace vall
} // namespace picking_old
} // namespace fragment
} // namespace core


#endif /* INCLUDED_core_fragment_picking_old_vall_eval_ABEGOEval_HH */
