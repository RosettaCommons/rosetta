// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    mp_mutate_relax.cc
/// @brief   Mutate a residue, then do quick relax for a membrane protein
/// @author  JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_MPMutateRelaxMover_hh
#define INCLUDED_protocols_membrane_MPMutateRelaxMover_hh

// Unit Headers
#include <protocols/relax/membrane/MPMutateRelaxMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project Headers

// Package Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh> // AUTO IWYU For ScoreFunctionOP

namespace protocols {
namespace membrane {

class MPMutateRelaxMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default Constructor for calls from JD2
	/// @details Mutations will be read in by commandline
	MPMutateRelaxMover();

	/// @brief Copy Constructor
	MPMutateRelaxMover( MPMutateRelaxMover const & src );

	/// @brief Assignment Operator
	MPMutateRelaxMover & operator = ( MPMutateRelaxMover const & src );

	/// @brief Destructor
	~MPMutateRelaxMover() override;

	///////////////////////////////
	/// Rosetta Scripts Methods ///
	///////////////////////////////

	/// @brief Create a Clone of this mover
	protocols::moves::MoverOP clone() const override;

	/// @brief Create a Fresh Instance of this Mover
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief Pase Rosetta Scripts Options for this Mover
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &
	) override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Get the name of this Mover (MPMutateRelaxMover)
	std::string get_name() const override;

	/// @brief Mutate residue and then quick relax the membrane protein
	void apply( core::pose::Pose & pose ) override;

private: // methods

	/////////////////////
	/// Setup Methods ///
	/////////////////////

	/// @brief Register Options from Command Line
	void register_options();

	/// @brief Initialize from commandline
	void init_from_cmd();

	/// @brief Get repack residues
	void get_repack_residues( Pose & pose );

	/// @brief Finalize setup
	void finalize_setup( Pose & pose );

	/// @brief Make mutations
	std::string make_mutations( Pose & pose, core::Size num_construct );

	////////////////////////////////////////////////////////////////////////////////
	/*
	THIS IS HOW THE INPUT FORMAT OF THE FILE LOOKS LIKE:
	= each line belongs to a single construct, i.e. a single sequence
	= a single entry (format A163F) is a single mutation, multiple mutations per construct are possible
	= example input:

	A163F
	W4N R27G K94E L45P
	G32V P34N

	= this means the first run is carried out for the single point mutation A163F
	= the second run of the mover is carried out with a quadrupel mutation
	= the third run is carried out with a double mutation
	... and so on.

	*/
	////////////////////////////////////////////////////////////////////////////////

	/// @brief Create output filename
	std::string output_filename( std::string mutation_tag, core::Size counter );

private: // data

	/// @brief Scorefunction
	core::scoring::ScoreFunctionOP sfxn_;

	/// @brief Input file containing desired mutants to make
	std::string mutant_file_;

	/// @brief Wildtype residue
	/// @details Outer vector is different constructs, inner vector
	///   is multiple residues per construct
	utility::vector1< utility::vector1< char > > wt_res_;

	/// @brief Pose residue numbers to mutate
	/// @details Outer vector is different constructs, inner vector
	///   is multiple residues per construct
	utility::vector1< utility::vector1< core::Size > > resn_;

	/// @brief Residue identities to mutate into, three letter code
	/// @details Outer vector is different constructs, inner vector
	///   is multiple residues per construct
	utility::vector1< utility::vector1< char > > mut_res_;

	/// @brief Number of iterations the mover runs on the inside, i.e. number of
	///   output models (Mover doesn't run in JD2!)
	core::Size nstruct_;

	/// @brief Protein name for dumping PDBs
	std::string protein_;

	/// @brief Repack?
	bool repack_mutation_only_;
	core::Real repack_radius_;
	utility::vector1< utility::vector1< bool > > repack_residues_;

	/// @brief Full relax?
	bool relax_;

};

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_MPMutateRelaxMover_hh
