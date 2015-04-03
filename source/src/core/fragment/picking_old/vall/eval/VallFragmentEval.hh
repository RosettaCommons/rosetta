// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragment/picking_old/vall/eval/VallFragmentEval.hh
/// @brief  base class for Vall ExtentEvaluator
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_fragment_picking_old_vall_eval_VallFragmentEval_hh
#define INCLUDED_core_fragment_picking_old_vall_eval_VallFragmentEval_hh

// unit headers
#include <core/fragment/picking_old/vall/eval/VallFragmentEval.fwd.hh>

// package headers
#include <core/fragment/picking_old/vall/VallSection.hh>
#include <core/fragment/picking_old/vall/scores/VallFragmentScore.hh>

// utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <core/fragment/picking_old/concepts/Extent.fwd.hh>
#include <core/fragment/picking_old/vall/VallLibrary.fwd.hh>
#include <utility/vector1.hh>


namespace core {
namespace fragment {
namespace picking_old {
namespace vall {
namespace eval {


/// @brief  base class for Vall ExtentEvaluator
class VallFragmentEval : public utility::pointer::ReferenceCount {


private: // typedefs


	typedef utility::pointer::ReferenceCount Super;


public: // typedefs


	typedef core::fragment::picking_old::vall::scores::VallFragmentScore VallFragmentScore;


public: // concept typedefs


	/// @brief typedef for ExtentEvaluator concept
	typedef VallFragmentScore::PageConstIterator PageConstIterator;


	/// @brief typedef for ExtentEvaluator concept
	typedef core::fragment::picking_old::concepts::Extent< VallSection::PageConstIterator > Extent;


public: // concept translation typedefs


	typedef PageConstIterator VallResidueConstIterator;


public: // construct/destruct


	/// @brief default constructor
	VallFragmentEval();


	/// @brief default copy constructor
	VallFragmentEval( VallFragmentEval const & rval );


	/// @brief default destructor
	virtual
	~VallFragmentEval();


public: // copy assignment


	/// @brief copy assignment
	VallFragmentEval & operator =( VallFragmentEval const & rval );


public: // virtual constructors


	/// @brief clone this object
	virtual
	VallFragmentEvalOP clone() const = 0;


public: // operators


	/// @brief called by VallLibrarian: for a fragment extent, evaluate and store
	///  results in a VallFragmentScore
	/// @return true if score should be stored, false otherwise
	bool operator ()(
		Extent const & extent,
		VallFragmentScore & fs
	);


public: // virtual evaluation methods


	/// @brief do the actual work of fragment evaluation
	/// @return true if score should be stored, false otherwise
	virtual
	bool eval_impl(
		Extent const & extent,
		VallFragmentScore & fs
	) = 0;


public: // additional hooks


	/// @brief operation to be perform before catalog() starts
	virtual
	void pre_catalog_op( VallLibrary const & )
	{}


	/// @brief operation to be performed after catalog() finished
	virtual
	void post_catalog_op( VallLibrary const & )
	{}


};


} // namespace eval
} // namespace vall
} // namespace picking_old
} // namespace fragment
} // namespace core


#endif /* INCLUDED_core_fragment_picking_old_vall_eval_VallFragmentEval_HH */
