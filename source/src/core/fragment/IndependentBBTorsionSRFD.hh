// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragment/IndependentBBTorsionSRFD.hh
/// @brief  A version of BBTorsionSRFD that considers each torsion independently
///         during is_applicable() and apply() calls when passed a MoveMap (vs
///         the all-torsions-must-be-moveable-or-nothing-is behavior in the
///         original BBTorsionSRFD).
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_fragment_IndependentBBTorsionSRFD_hh
#define INCLUDED_core_fragment_IndependentBBTorsionSRFD_hh


// unit headers
#include <core/fragment/IndependentBBTorsionSRFD.fwd.hh>

// package headers
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/SecstructSRFD.hh>

#include <core/types.hh>

#include <utility/vector1.hh>


namespace core {
namespace fragment {


/// @brief A version of BBTorsionSRFD that considers each torsion independently
///  during is_applicable() and apply() calls when passed a MoveMap (vs the
///  all-torsions-must-be-moveable-or-nothing-is behavior in the original
///  BBTorsionSRFD).
class IndependentBBTorsionSRFD : public BBTorsionSRFD {


private: // typedefs


	typedef BBTorsionSRFD Super;
	typedef SecstructSRFD Super2;


public: // typedefs


	typedef core::Size Size;
	typedef core::kinematics::MoveMap MoveMap;
	typedef core::pose::Pose Pose;


public: // construct/destruct


	/// @brief default constructor
	IndependentBBTorsionSRFD();


	/// @brief constructor
	/// @param[in] n_bbtorsions Number of backbone torsions.
	/// @param[in] secstruct The single character secondary structure type.
	/// @param[in] sequence The single character sequence type.
	IndependentBBTorsionSRFD(
		Size const n_bbtorsions,
		char const secstruct,
		char const sequence
	);


	/// @brief copy constructor
	IndependentBBTorsionSRFD( IndependentBBTorsionSRFD const & rval );


	/// @brief default destructor
	
	~IndependentBBTorsionSRFD() override;


public: // assignment


	/// @brief copy assignment
	IndependentBBTorsionSRFD & operator =( IndependentBBTorsionSRFD const & rval );


public: // virtual constructors


	/// @brief clone this object
	
	SingleResidueFragDataOP clone() const override;


	/// @brief create a new instance of this object
	
	SingleResidueFragDataOP create() const override;


public: // methods


	/// @brief apply only torsions in this fragment marked as moveable in the given
	///  MoveMap
	/// @param[in] movemap Check for moveable torsions in this MoveMap.
	/// @param[in,out] pose The Pose to modify.
	/// @param[in] seqpos Insert at this sequence position.
	/// @return True if at least one torsion inserted and second level superclass
	///  <tt>SecstructSRFD::apply()</tt> succeeded, otherwise false.
	
	bool apply(
		MoveMap const & movemap,
		Pose & pose,
		Size const seqpos
	) const override;


	/// @brief is at least one torsion marked as moveable in the given MoveMap?
	/// @param[in] movemap Check for moveable torsions in this MoveMap.
	/// @param[in] seqpos Check at this sequence position.
	/// @return True if at least one torsion moveable and second level superclass
	///  <tt>SecstructSRFD::is_applicable()</tt>, otherwise False.
	
	bool is_applicable(
		MoveMap const & movemap,
		Size seqpos
	) const override;


};


} // namespace fragment
} // namespace core


#endif /* INCLUDED_core_fragment_IndependentBBTorsionSRFD_HH */
