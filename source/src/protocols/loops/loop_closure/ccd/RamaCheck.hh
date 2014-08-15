// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/loops/loop_closure/ccd/RamaCheck.hh
/// @brief  Header for RamaCheck classes (RamaCheckBase, RamaCheck1B, RamaCheck2B)
/// @author Brian D. Weitzner


#ifndef INCLUDED_protocols_loops_loop_closure_ccd_RamaCheck_HH
#define INCLUDED_protocols_loops_loop_closure_ccd_RamaCheck_HH

// Unit Headers
#include <protocols/loops/loop_closure/ccd/RamaCheck.fwd.hh>

// Project Headers
#include <core/id/TorsionID.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/Ramachandran.fwd.hh>
#include <core/scoring/Ramachandran2B.fwd.hh>

#include <core/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

// C++ Headers
#include <iostream>

namespace protocols {
namespace loops {
namespace loop_closure {
namespace ccd {

class RamaCheckBase : public utility::pointer::ReferenceCount {

typedef utility::pointer::ReferenceCount Parent;
typedef utility::vector1< core::Real > RamaScoreVector;

public:
	/// @brief constructor
	RamaCheckBase();

	/// @brief  Copy constructor
	RamaCheckBase( RamaCheckBase const & object_to_copy );

	// Assignment operator
	RamaCheckBase & operator=( RamaCheckBase const & object_to_copy );

	// destructor
	virtual ~RamaCheckBase();

	/// @brief  Generate a string representation of RamaCheck for debugging purposes.
	virtual void show( std::ostream & output=std::cout ) const;

	virtual void parse_my_tag( utility::tag::TagCOP tag );

	/// @brief  Register options with the option system.
	static void register_options();

	/// @brief Return the name of the RamaCheck class being used
	virtual std::string name() const = 0;

	/// @brief Return a pointer to a new, fully configured copy of RamaCheckBase
	/// @note  The clone method should be used to copy an instance of RamaCheck if the derived class type is unknown.
	virtual RamaCheckBaseOP
	clone() const =0;
	
public: // accessors and mutators
	/// @brief  Get the "temperature" used for rama score-checking with the Metropolis criterion.
	/// @note   This value will only be used if the option for checking the rama score is turned on.
	core::Real temperature() const { return temperature_; }

	/// @brief  Set the "temperature" used for rama score-checking with the Metropolis criterion.
	/// @note   This value will only be used if the option for checking the rama score is turned on.
	void temperature( core::Real input_temperature ) { temperature_ = input_temperature; }

	/// @brief  Get the maximum rama score increase allowed during rama score-checking with the Metropolis criterion.
	/// @note   This value will only be used if the option for checking the rama score is turned on.
	core::Real
	max_rama_score_increase() const
	{
		return max_rama_score_increase_;
	}

	/// @brief  Set the maximum rama score increase allowed during rama score-checking with the Metropolis criterion.
	/// @note   This value will only be used if the option for checking the rama score is turned on.
	void
	max_rama_score_increase( core::Real input_max_rama_score_increase )
	{
		max_rama_score_increase_ = input_max_rama_score_increase;
	}

public: // methods used to compute Ramachandran checks
	/// @brief Store the Ramachandran scores of each residue in the supplied pose.
	void initialize_starting_rama_scores( core::pose::Pose const & pose ) const;

	/// @brief Determine whether or not a candidate conformation should be accepted based on the Ramachandran score.
	/// @return true if the move should be accepted, false if it should be rejected.
	bool
	accept_new_conformation(
		core::pose::Pose const & pose,
		core::id::TorsionID const & torsion_id,
		core::Angle const alpha
	) const;

	/// @brief Compute the total net change in Ramachandran score between the initial pose and the current pose from
	/// <first_res> to <last_res>.
	/// @return The total net change in Ramachandran score.
	core::Real
	total_net_change_in_rama_score_over_range(
		core::pose::Pose const & pose,
		core::uint const first_res,
		core::uint const last_res
	) const;

	/// @brief Compute the average change in Ramachandran score between the initial pose and the current pose from
	/// <first_res> to <last_res>.
	/// @return The average change in Ramachandran score.
	core::Real
	average_change_in_rama_score_over_range(
		core::pose::Pose const & pose,
		core::uint const first_res,
		core::uint const last_res
	) const;

	/// @brief Compute the Ramachandran score of residue <seqpos> in <pose> with a hypothetical conformation <phi>, <psi>
	virtual core::Real
	compute_rama_score(
		core::pose::Pose const & pose,
		core::uint const seqpos,
		core::Real const phi,
		core::Real const psi
	) const = 0;

private: // methods
	// Copy all data members from <from> to <to>.
	void copy_data( RamaCheckBase & to, RamaCheckBase const & from ) const;

	// setup private data
	void init();

	// Initialize data members from option system.
	void init_options();

	// acceptance with scores supplied
	bool
	accept_new_conformation(
		core::Real const current_score,
		core::Real const proposed_score,
		core::Real const starting_score
	) const;

private: // data
	core::Real temperature_;
	core::Real max_rama_score_increase_;
	mutable RamaScoreVector starting_rama_scores_;

private: // static constants
	static core::Real const BAD_SCORE;

}; // RamaCheckBase

class RamaCheck1B : public RamaCheckBase {
public:
	/// @brief constructor
	RamaCheck1B();

	/// @brief  Copy constructor
	RamaCheck1B( RamaCheck1B const & object_to_copy );

	// Assignment operator
	RamaCheck1B & operator=( RamaCheck1B const & object_to_copy );

	// destructor
	virtual ~RamaCheck1B();

	/// @brief Return "RamaCheck1B"
	virtual std::string name() const;

	/// @brief Return a pointer to a new, fully configured copy of RamaCheck1B
	virtual RamaCheckBaseOP
	clone() const;

	/// @brief Compute the Ramachandran score of residue <seqpos> in <pose> with a hypothetical conformation <phi>, <psi>.
	/// The score is independent of the identity of the neighboring residues.
	virtual core::Real
	compute_rama_score(
		core::pose::Pose const & pose,
		core::uint const seqpos,
		core::Real const phi,
		core::Real const psi
	) const;

private: //methods
	// Copy all data members from <from> to <to>.
	void copy_data( RamaCheck1B & to, RamaCheck1B const & from ) const;

private: // data
	mutable core::scoring::RamachandranOP rama_;

}; // RamaCheck1B

class RamaCheck2B : public RamaCheckBase {
public:
	/// @brief constructor
	RamaCheck2B();

	/// @brief  Copy constructor
	RamaCheck2B( RamaCheck2B const & object_to_copy );

	// Assignment operator
	RamaCheck2B & operator=( RamaCheck2B const & object_to_copy );

	// destructor
	virtual ~RamaCheck2B();

	/// @brief Return "RamaCheck2B"
	virtual std::string name() const;

	/// @brief Return a pointer to a new, fully configured copy of RamaCheck2B
	virtual RamaCheckBaseOP
	clone() const;

	/// @brief Compute the Ramachandran score of residue <seqpos> in <pose> with a hypothetical conformation <phi>, <psi>.
	/// The score depends on the identity of the neighboring residues.
	virtual core::Real
	compute_rama_score(
		core::pose::Pose const & pose,
		core::uint const seqpos,
		core::Real const phi,
		core::Real const psi
	) const;

private: // methods
	// Copy all data members from <from> to <to>.
	void copy_data( RamaCheck2B & to, RamaCheck2B const & from ) const;

private: // data
	mutable core::scoring::Ramachandran2BOP rama_;

}; // RamaCheck2B

// Insertion operator (overloaded so that RamaCheck can be "printed" in PyRosetta).
std::ostream & operator<< ( std::ostream & os, RamaCheckBase const & rama_check );

} // namespace ccd
} // namespace loop_closure
} // namespace loops
} // namespace protocols

#endif // INCLUDED_protocols_loops_loop_closure_ccd_RamaCheck_HH
