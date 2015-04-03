// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragment/picking_old/vall/eval/EnergyEval.hh
/// @brief  scores a fragment by inserting its backbone angles into a Pose
///         and evaluating its energy using a given ScoreFunction
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_fragment_picking_old_vall_eval_EnergyEval_hh
#define INCLUDED_core_fragment_picking_old_vall_eval_EnergyEval_hh


// unit headers
#include <core/fragment/picking_old/vall/eval/EnergyEval.fwd.hh>

// type headers
#include <core/types.hh>

// package headers
#include <core/fragment/picking_old/vall/eval/VallFragmentEval.hh>
#include <core/fragment/picking_old/vall/VallLibrary.fwd.hh>

// project headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>


// C++ headers
#include <iostream>
#include <string>

#include <utility/vector1.hh>


namespace core {
namespace fragment {
namespace picking_old {
namespace vall {
namespace eval {


/// @brief scores a fragment by inserting its backbone angles into a Pose and
///  evaluating its energy using a given ScoreFunction
/// @remarks Useful for e.g. evaluating a chainbreak score, constraint score, etc.
/// @warning Class currently assumes insertion of backbone angles into protein
///  residues.
class EnergyEval : public VallFragmentEval {


private: // typedefs


	typedef VallFragmentEval Super;


public: // typedefs


	typedef Super::VallFragmentScore VallFragmentScore;

	typedef core::pose::Pose Pose;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef core::scoring::ScoreFunctionCOP ScoreFunctionCOP;


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
	EnergyEval();


	/// @brief constructor
	/// @param[in] pose insert backbone angles using a copy of this Pose
	/// @param[in] insert_position insert backbone angles starting from this
	///  position in the Pose
	/// @param[in] score_function evaluate the Pose using a copy of this
	///  ScoreFunction
	/// @param[in] randomize flags that indicates whether a small amount
	///  of noise between [0, 0.000001) will be added to the energy
	EnergyEval(
		Pose const & pose,
		Size const insert_position,
		ScoreFunction const & score_function,
		bool const randomize = false
	);


	/// @brief default copy constructor
	EnergyEval( EnergyEval const & rval );


	/// @brief default destructor
	virtual
	~EnergyEval();


public: // copy assignment


	/// @brief copy assignment
	EnergyEval & operator =( EnergyEval const & rval );


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


public: // accessor


	/// @brief the Pose to evaluate
	inline
	Pose const & pose() const {
		return pose_;
	}


	/// @brief insert angles into the Pose starting from this position
	inline
	Size insert_position() const {
		return insert_position_;
	}


	/// @brief the ScoreFunction used to evaluate the pose
	inline
	ScoreFunctionCOP score_function() const {
		return score_function_;
	}


	/// @brief adding random noise between [0, 0.000001) to the energy?
	inline
	bool randomize() const {
		return randomize_;
	}


public: // mutators


	/// @brief the Pose to evaluate
	inline
	void pose( Pose const & p ) {
		pose_ = p;
	}


	/// @brief insert angles into the Pose starting from this position
	inline
	void insert_position( Size const position ) {
		insert_position_ = position;
	}


	/// @brief the ScoreFunction used to evaluate the pose
	inline
	void score_function( ScoreFunction const & fx ) {
		score_function_ = fx.clone();
	}


	/// @brief set flag to add random noise between [0, 0.000001) to the energy
	inline
	void randomize( bool const flag ) {
		randomize_ = flag;
	}


public: // additional hooks


	/// @brief operation to be perform before catalog() starts
	virtual
	void pre_catalog_op( VallLibrary const & );


private: // data


	/// @brief insert backbone angles into this Pose
	Pose pose_;


	/// @brief insert backbone angles into pose starting from this position
	Size insert_position_;


	/// @brief evaluate the Pose using this ScoreFunction
	ScoreFunctionOP score_function_;


	/// @brief flag to add random noise between [0, 0.000001) into the energy
	bool randomize_;


};


} // namespace eval
} // namespace vall
} // namespace picking_old
} // namespace fragment
} // namespace core


#endif /* INCLUDED_core_fragment_picking_old_vall_eval_EnergyEval_HH */
