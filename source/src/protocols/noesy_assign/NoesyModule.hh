// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file NoesyModule.hh
/// @brief main hook-up for the automatic NOESY assignment module
/// @details
///   handling of input-output options
///       provides class NoesyModule
///
/// @author Oliver Lange

#ifndef INCLUDED_protocols_noesy_assign_NoesyModule_HH
#define INCLUDED_protocols_noesy_assign_NoesyModule_HH


// Unit Header
#include <protocols/noesy_assign/NoesyModule.fwd.hh>

// Package Headers
//#include <devel/noesy_assign/ResonanceList.fwd.hh>

// Project Headers
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <core/pose/Pose.fwd.hh>
#include <protocols/noesy_assign/CrossPeakList.fwd.hh>
#include <protocols/noesy_assign/ResonanceList.fwd.hh>
#include <string>


//// C++ headers
// #include <cstdlib>
// #include <string>
// #include <list>
// #include <map>

namespace protocols {
namespace noesy_assign {

/// WARNING WARNING WARNING THREAD UNSAFE FOR USING THE COVALENTCOMPLIANCE CLASS IN A NON-CONST WAY
class NoesyModule : public utility::pointer::ReferenceCount {

private:
	static bool options_registered_;

public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~NoesyModule();
	/// @brief register options
	static void register_options();

public:

	/// @brief constructor  -- initialize with fast-sequence
	NoesyModule( std::string const& fasta_sequence );

	/// @brief assign NOE data, use models provided by DecoyIterator for scoring and restraint exclusion, if cycle=0 read cycle from cmd-line
	template< class DecoyIterator >
	void assign( DecoyIterator const& begin, DecoyIterator const& end );

	/// @brief same as above, but decoy file will be determined from commandline and read directly
	void assign ();

	/// @brief after assignment --> create appropriate constraints
	void generate_constraint_files(
		core::pose::Pose const& pose,
		std::string const& cst_fa_file,
		std::string const& cst_centroid_file,
		core::Size min_seq_separation,
		core::Size min_quali = 0,
		core::Size max_quali = 4,
		bool ignore_elimination_candidates = true,
		bool elimination_candidates = false
	) const;

	/// @brief write assignments into peak-file...
	void write_assignments( std::string file = "use_cmd_line" );

	/// @brief reset assignments...  -- call before re-assigning peaks
	void reset();

	/// @brief returns true if -noesy::in::resonances and -noesy::in::peaks are in cmd-line
	static bool cmdline_options_activated();


	/// @brief return (cross)peak-list (peak-file)
	CrossPeakList const& crosspeaks() const { return *crosspeaks_; }

	/// @brief return resonance assignments (prot-file)
	//  ResonanceList const& resonances() const { return *main_resonances_; }

	void add_dist_viol_to_assignments( core::pose::Pose native_pose);

	std::string const& sequence() {
		return sequence_;
	}

private:

	/// @brief return all input files
	void read_input_files();

	//  bool skip_network_analysis_; //moved to PeakAssignmentParameters
	/// @brief private data, peak-list and master-resonances (sometimes different resonances for different peak-lists, thus the name)
	CrossPeakListOP crosspeaks_;
	//ResonanceListOP main_resonances_;
	std::string sequence_;


};

}
}

#endif
