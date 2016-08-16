// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file CrossPeakList.hh
/// @author Oliver Lange

#ifndef INCLUDED_protocols_noesy_assign_PeakAssignmentResidueMap_HH
#define INCLUDED_protocols_noesy_assign_PeakAssignmentResidueMap_HH


// Unit Headers
#include <protocols/noesy_assign/PeakAssignmentResidueMap.fwd.hh>
#include <protocols/noesy_assign/CrossPeakList.fwd.hh>
#include <core/types.hh>

// Package Headers
#include <protocols/noesy_assign/PeakAssignment.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>

// Project Headers

// Utility headers
// #include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
// #include <numeric/numeric.functions.hh>
// #include <core/util/prof.hh>
//#include <core/util/Tracer.hh>
// #include <core/options/option.hh>
// #include <core/options/keys/abinitio.OptionKeys.gen.hh>
// #include <core/options/keys/run.OptionKeys.gen.hh>
//#include <core/options/keys/templates.OptionKeys.gen.hh>

//// C++ headers
#include <list>
#include <set>
#include <map>

#include <protocols/noesy_assign/ResonanceList.fwd.hh>


namespace protocols {
namespace noesy_assign {


/// @brief fast access to assignments by residue number
class PeakAssignmentResidueMap : public utility::pointer::ReferenceCount {
public:
	typedef std::list< PeakAssignmentOP > PeakAssignments;
	typedef std::map< core::Size, PeakAssignments > PeakAssignmentMap;
	typedef utility::vector1< PeakAssignmentMap > ResidueList;

	PeakAssignmentResidueMap();
	virtual ~PeakAssignmentResidueMap();

	/// @brief add all PeakAssignments in all Crosspeaks of list
	void add( CrossPeakList const& );

	/// @brief add individual PeakAssignment
	void add( PeakAssignmentOP const& );


	/// @brief add all resonances for backward compatibility in covalent-part of network-analysis
	void add_all_atoms( ResonanceList const& );

	/// @brief remove individual PeakAssignment
	void remove( PeakAssignment const& );

	/// @brief invalidate non symmetric peaks
	void check_for_symmetric_peaks( CrossPeakList&, bool accumulate_symmetry );

	/// @brief remove all ambiguous assignments to i,i+1 CrossPeaks.

	//commented out 9/10/12 because the only relevant action within the subroutine is commented out already
	//void invalidate_competitors_to_sequential_NOE( CrossPeakList& );

	void network_analysis( Size n_total_assignments );
	void network_analysis2(); // ResonanceList const& resonances );
	/// @brief get list of PeakAssignments for pair of residues --- throws Exception
	PeakAssignments const& assignments( core::Size resi, core::Size resj ) const;
	PeakAssignments& assignments( core::Size resi, core::Size resj );

	/// @brief add assignments found between resi and resj to collector
	void assignments( core::Size resi, core::Size resj, PeakAssignments& collector ) const;

	/// @brief has some (valid or invalid) assignments between residue pair
	bool has( core::Size resi, core::Size resj );

	core::Size total_residue() const {
		return residues_.size();
	}

private:
	/// @brief same as "assignments()" but returns BOGUS_ASSIGNMENTS if not found
	PeakAssignments const& _assignments( core::Size resi, core::Size resj ) const;
	PeakAssignments& _assignments( core::Size resi, core::Size resj );

	/// @brief subroutine to compute Nk for alpha->gamma->beta path.
	core::Real compute_Nk(
		PeakAssignment const& alpha_beta,
		core::id::NamedAtomID const& gamma_atom,
		bool connect_in_i,
		bool connect_in_j,
		bool sequential,
		PeakAssignments const& close_to_i_assignments,
		PeakAssignments const& close_to_j_assignments,
		core::Real longrange_peak_volume
		/*either the a-g or g-b peak_volume (corresponding to the "false" in one of the connect_in_i/j */
	) const;

	/// @brief subroutine to collect putative gammas that need to be queried due to covalent structure
	void fill_covalent_gammas( Size, std::map< core::id::NamedAtomID, bool >& collector ) const;


	//find out if we have an initial assignment to gamma
	/// map with resid of atom(2) as search key.
	///per residue - a PeakAssignmentMap --- all residues that are connected by an initial assignment
	ResidueList residues_;
	PeakAssignments BOGUS_ASSIGNMENTS; //emtpy list

	// a little redundancy in the residue number
	// but saves lot of time later we can just reference these objects.
	typedef std::set< core::id::NamedAtomID > AtomList;
	typedef utility::vector1< AtomList > AtomByResList;

	//for reach residue list of NamedAtomIDs
	AtomByResList atoms_;
};

}
}

#endif
