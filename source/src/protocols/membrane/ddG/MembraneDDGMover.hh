// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/MembraneDDGMover.hh
///
/// @brief      Compute ddG Scores for Membrane Protein
/// @details	Initialize a membrane pose, compute an initial membrane position,
///				compute a native score, make mutation, repack sidechains, and score new
///				structure. Uses a user-provided resfile for repacking/mutations
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (7/7/14)

#ifndef INCLUDED_protocols_membrane_ddG_MembraneDDGMover_hh
#define INCLUDED_protocols_membrane_ddG_MembraneDDGMover_hh

// Unit Headers
#include <protocols/membrane/ddG/MembraneDDGMover.fwd.hh>

// Package headers
#include <protocols/moves/Mover.hh>
#include <protocols/membrane/ddG/Mutation.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>

// C++ Headers
#include <cstdlib>
#include <set>

namespace protocols {
namespace membrane {
namespace ddG {

using namespace core;
using namespace core::scoring;
using namespace protocols::moves;

/// @brief Mutate Residue, Repack within pack_radius and compute
/// ddG score of mutation
class MembraneDDGMover : public Mover {
	
public:
	
	////////////////////
	/// Constructors ///
	////////////////////

	/// @brief Default Constructor
	/// @details Construct a defauilt version of this mover: pack_radius = 0.0, membrane
	/// env smooth sfxn
	MembraneDDGMover();
	
	/// @brief Custom Constructor
	/// @details Specify a pack radius, scorefunction to use, etc.
	MembraneDDGMover(
		Real pack_radius_in,
		ScoreFunctionOP sfxn_in,
		utility::vector1< MutationOP > mutations_in
	);
	
	/// @brief Copy Constructor
	/// @details Create a deep copy of this object
	MembraneDDGMover( MembraneDDGMover const & src );
	
	/// @brief Assignment Operator
	/// @details Create a deep copy of this object overloading the assignment operator
	MembraneDDGMover &
	operator=( MembraneDDGMover const & src );
	
	/// @brief Destructor
	~MembraneDDGMover();
	
	///////////////////////////////
	/// Rosetta Scripts Methods ///
	///////////////////////////////
	
	/// @brief Create a Clone of this mover
	virtual protocols::moves::MoverOP clone() const;
	
	/// @brief Create a Fresh Instance of this Mover
	virtual protocols::moves::MoverOP fresh_instance() const;
	
	/// @brief Pase Rosetta Scripts Options for this Mover
	void parse_my_tag(
	  utility::tag::TagCOP tag,
	  basic::datacache::DataMap &,
	  protocols::filters::Filters_map const &,
	  protocols::moves::Movers_map const &,
	  core::pose::Pose const &
	  );
	
	/////////////////////
	/// Mover Methods ///
	/////////////////////
	
	/// @brief Get the name of this Mover (MembraneDDGMover)
	virtual std::string get_name() const;
	
	/// @brief Compute ddG Scores for Membranes
	virtual void apply( Pose & pose );
	
private:

/// @brief Compute ddG Score
/// @details Compute ddG from a copy of the native pose, native score
/// and provided resfile. Doing this PyRosetta Style - thanks to Evan Baugh's Mutate.py for
/// instructions
core::Real
compute_ddG_score(
	Pose & pose,
	AA aa,
	Size position,
	core::Real native_score
	);
	
/// @brief List of AAs to test
/// @details Return a vector 1 of characters representing the 20 AAs to substitute. This is
/// possibly too specific to my current task as well
utility::vector1< char >
designed_amino_acids();

/// @brief Access Neighbors within Radius
/// @details Determine the number of neighbors within
/// the user provided radius
std::set< Size >
get_residue_neighbors( Pose & pose, core::Size position, core::Real radius );

private:

	// Pack Radius
	Real pack_radius_;

	// ScoreFunction
	ScoreFunctionOP sfxn_;
	
	// List of mutations
	utility::vector1< MutationOP > mutations_;
	
};


} // ddG
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_ddG_MembraneDDGMover_hh
