// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    mp_mutate_relax.cc
/// @brief   Mutate a residue, then do quick relax for a membrane protein
/// @author  JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_MPMutateRelaxMover_hh
#define INCLUDED_protocols_membrane_MPMutateRelaxMover_hh

// Unit Headers
#include <protocols/membrane/MPMutateRelaxMover.fwd.hh>
#include <protocols/membrane/MPMutateRelaxMoverCreator.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/kinematics/MoveMap.fwd.hh>

// Package Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.fwd.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/Tracer.fwd.hh>

namespace protocols {
namespace membrane {

using namespace core;
using namespace core::pose;
using namespace protocols::moves;
using namespace core::kinematics;

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
	virtual ~MPMutateRelaxMover();

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

	/// @brief Get the name of this Mover (MPMutateRelaxMover)
	virtual std::string get_name() const;

	/// @brief Mutate residue and then quick relax the membrane protein
	virtual void apply( Pose & pose );

private: // methods

	/////////////////////
	/// Setup Methods ///
	/////////////////////

	/// @brief Register Options from Command Line
	void register_options();

	/// @brief Initialize from commandline
	void init_from_cmd();

	////////////////////////////////////////////////////////////////////////////////
	/*
	THIS IS HOW THE INPUT FORMAT OF THE FILE LOOKS LIKE:
	= each line belongs to a single construct, i.e. a single run
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
	/// @brief Read mutant file
	void read_mutant_file();

	/// @brief Add mutants to private data: A163F into vectors
	/// @details This is an entire line of the mutant input file
	void add_mutant_to_vectors( std::string mutations );

	/// @brief Check mutant file for errors
	/// @details If Rosetta doesn't start crying, you're good to go
	void check_mutant_file( Pose & pose );

	/// @brief Convert one-letter code to three-letter code
	std::string one2three( std::string one );

	/// @brief Create output filename
	std::string output_filename( std::string mutation_tag, Size counter );

private: // data

	/// @brief Input file containing desired mutants to make
	std::string mutant_file_;

	/// @brief Wildtype residue
	utility::vector1< utility::vector1< std::string > > wt_res_;

	/// @brief Pose residue numbers to mutate
	utility::vector1< utility::vector1< Size > > resn_;

	/// @brief Residue identities to mutate into, three letter code
	utility::vector1< utility::vector1< std::string > > new_res_;

	/// @brief Number of iterations the mover runs on the inside
	Size iter_;

	/// @brief Protein name for dumping PDBs
	std::string protein_;

};

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_MPMutateRelaxMover_hh
