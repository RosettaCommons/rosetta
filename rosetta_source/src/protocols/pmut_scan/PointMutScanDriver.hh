// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file protocols/stability/point_mut_scan.hh
/// @brief A protocol that tries to find stability enhancing mutations
/// @author Ron Jacak

#ifndef INCLUDED_protocols_pmut_scan_PointMutScanDriver_HH
#define INCLUDED_protocols_pmut_scan_PointMutScanDriver_HH

// MPI Headers have to be #included first
#ifdef USEMPI
#include <mpi.h>
#endif

// Project Headers
#include <protocols/pmut_scan/Mutant.fwd.hh>

#include <core/graph/Graph.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>

#include <basic/options/util.hh>

// Utility headers

// ObjexxFCL header

// C++
#include <string>
#include <iostream>

namespace protocols {
namespace pmut_scan {

class PointMutScanDriver {

public:
	PointMutScanDriver( utility::vector1< std::string > & pdb_file_names, bool double_mutant_scan, std::string list_file, bool output_mutant_structures );
	~PointMutScanDriver();

	void go();

	void read_in_structures();
	void fill_mutations_list();

	void calculate_neighbor_table( core::pose::Pose & pose, utility::vector1< utility::vector1< bool > > & neighbors );

	void make_specific_mutant( utility::vector1< core::pose::Pose > & mutant_poses, utility::vector1< core::pose::Pose > & native_poses,
		core::scoring::ScoreFunctionOP scorefxn, protocols::pmut_scan::Mutant & m, std::string mutation_string = "" );

	void set_ddG_cutoff( core::Real threshold );
	
	utility::vector1< Mutant >::const_iterator mutants_begin() const; // used only by unit tests
	utility::vector1< Mutant >::const_iterator mutants_end() const;   // used only by unit tests
	
	core::Size n_mutants() const; // used only by unit tests


private:
	void barrier();
	std::string node_name( core::Size rank );
	
	void read_mutants_list_file( std::string & list_file );
	void divide_up_mutations();

#ifdef USEMPI
	static void send_mutant_data_to_node( int destination, protocols::pmut_scan::Mutant const & m );
	static protocols::pmut_scan::Mutant receive_mutant_data_from_node( int source );
#endif

	void make_mutants();

	void make_mutant_structure( core::pose::Pose & mutant_pose, core::pose::Pose & native_pose, protocols::pmut_scan::MutationData const & md, core::scoring::ScoreFunctionOP scorefxn );


private:
	bool double_mutant_scan_;
	std::string mutants_list_file_;
	bool output_mutant_structures_;

#ifdef USEMPI
	MPI_Status stat_;
	int tag_;
#endif

	core::Size MPI_rank_;
	core::Size MPI_nprocs_;

	utility::vector1< core::pose::Pose > input_poses_;

	utility::vector1< Mutant > all_mutants_;
	utility::vector1< Mutant > mutants_list_;

	utility::vector1< std::string > pdb_file_names_;

	core::Real DDG_CUTOFF_;


}; // class PointMutScanDriver

} // namespace pmut_scan
} // namespace protocols


#endif
