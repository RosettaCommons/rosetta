// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pmut_scan/PointMutScanDriver.hh
/// @brief A protocol that tries to find stability enhancing mutations
/// @author Ron Jacak (ron.jacak@gmail.com)

#ifndef INCLUDED_protocols_pmut_scan_PointMutScanDriver_HH
#define INCLUDED_protocols_pmut_scan_PointMutScanDriver_HH

// MPI Headers have to be #included first
#ifdef USEMPI
#include <mpi.h>
#endif

// Project Headers
#include <protocols/pmut_scan/Mutant.fwd.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

// ObjexxFCL header

// C++
#include <string>

namespace protocols {
namespace pmut_scan {

class PointMutScanDriver {

public:
	PointMutScanDriver( utility::vector1< std::string > & pdb_file_names, bool double_mutant_scan, std::string list_file, bool output_mutant_structures );
	virtual ~PointMutScanDriver();

	//This protocol is its own Job Distributor - this fires it off from application-level
	void go();

	//mutant-scanning related functions

	void read_in_structures();
	void fill_mutations_list();

	void calculate_neighbor_table( core::pose::Pose & pose, utility::vector1< utility::vector1< bool > > & neighbors );

	void make_specific_mutant( utility::vector1< core::pose::Pose > & mutant_poses, utility::vector1< core::pose::Pose > & native_poses, protocols::pmut_scan::Mutant & m, std::string mutation_string = "", std::string mutation_string_PDB_numbering = "" );

	/// @brief score the pose for the purposes of determining if a mutation is "good" or not.  In the base implementation, it's just a scorefunction call, but in child implementations it may be fancier (for example, calculating a binding energy instead)
	virtual core::Energy score(core::pose::Pose & pose);

	/// @brief accessor for scorefxn_ now that it is private member data
	core::scoring::ScoreFunctionCOP get_scorefxn() const;// { return scorefxn_; }

	/// @brief offers a chance for child classes to inject mutant selection logic
	virtual bool reject_mutant( Mutant const & /*m*/, core::pose::Pose const & /*pose*/ ) { return false; }

	//unit test utilities
	void set_ddG_cutoff( core::Real threshold );

	utility::vector1< Mutant >::const_iterator mutants_begin() const; // used only by unit tests
	utility::vector1< Mutant >::const_iterator mutants_end() const;   // used only by unit tests

	core::Size n_mutants() const; // used only by unit tests

private:

	void make_mutants();

	void make_mutant_structure( core::pose::Pose & mutant_pose, core::pose::Pose & native_pose, protocols::pmut_scan::MutationData const & md );

private: //mutant scanning data
	bool double_mutant_scan_;
	std::string mutants_list_file_;
	bool output_mutant_structures_;

	utility::vector1< core::pose::Pose > input_poses_;

	utility::vector1< Mutant > all_mutants_;
	utility::vector1< Mutant > mutants_list_;

	utility::vector1< std::string > pdb_file_names_;

	core::Real DDG_cutoff_;

	core::scoring::ScoreFunctionOP scorefxn_;

private: //Job Distribution related functions
	void barrier();
	std::string node_name( core::Size rank );

	void read_mutants_list_file( std::string & list_file );
	void divide_up_mutations();

#ifdef USEMPI
	static void send_mutant_data_to_node( int destination, protocols::pmut_scan::Mutant const & m );
	static protocols::pmut_scan::Mutant receive_mutant_data_from_node( int source );
#endif

private: //Job Distribution related data

#ifdef USEMPI
	MPI_Status stat_;
	int tag_;
#endif

	core::Size MPI_rank_;
	core::Size MPI_nprocs_;

}; // class PointMutScanDriver

} // namespace pmut_scan
} // namespace protocols

#endif //INCLUDED_protocols_pmut_scan_PointMutScanDriver_HH
