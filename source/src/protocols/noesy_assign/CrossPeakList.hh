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

#ifndef INCLUDED_protocols_noesy_assign_CrossPeakList_HH
#define INCLUDED_protocols_noesy_assign_CrossPeakList_HH


// Unit Header
#include <protocols/noesy_assign/CrossPeakList.fwd.hh>

// Package Headers
#include <protocols/noesy_assign/PeakFileFormat.fwd.hh>
#include <protocols/noesy_assign/PeakAssignmentResidueMap.fwd.hh>
#include <protocols/noesy_assign/PeakCalibrator.fwd.hh>
//#include <devel/noesy_assign/PeakAssignment.hh>
//#include <devel/noesy_assign/ResonanceList.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>

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

#include <protocols/noesy_assign/CrossPeak.fwd.hh>
#include <protocols/noesy_assign/ResonanceList.fwd.hh>


namespace protocols {
namespace noesy_assign {

class CrossPeakList : public utility::pointer::ReferenceCount {
public:
	typedef std::list< CrossPeakOP > CrossPeaks;
	typedef CrossPeaks::const_iterator const_iterator;
	typedef CrossPeaks::iterator iterator;
	CrossPeakList();
	~CrossPeakList() override;

	void read_from_stream( std::istream&, PeakFileFormat& input_adaptor, ResonanceListOP resonances );
	void write_to_stream( std::ostream&, PeakFileFormat& output_adaptor ) const;
	void write_peak_files( std::string const& prefix, PeakFileFormat& output_adaptor ) const;
	void find_assignments();
	void update_chemshiftscore();
	void update_symmetry_score();
	void update_upperdistance_score();

	template < class DecoyIterator >
	void update_decoy_compatibility_score( DecoyIterator const& begin, DecoyIterator const& end ); //core::io::silent::SilentFileData const& sfd

	void eliminate_spurious_peaks();

	template < class DecoyIterator >
	void calibrate( DecoyIterator const& begin, DecoyIterator const& end );
	//  void calibrate( core::io::silent::SilentFileData const& decoys );

#if 0
	core::scoring::constraints::ConstraintSetOP generate_constraints( core::pose::Pose const& pose, bool centroid = false, core::Size min_seq_separation = 2 ) const;
#endif

	void generate_fa_and_cen_constraints(
		core::scoring::constraints::ConstraintSetOP fa_set,
		core::scoring::constraints::ConstraintSetOP cen_set,
		core::pose::Pose const& pose,
		core::pose::Pose const& centroid_pose,
		core::Size min_seq_separation,
		core::Size min_quali,
		core::Size max_quali,
		core::Real padding = 0.0,
		bool ignore_elimination_candidates = true,
		bool elimination_candidates = false
	) const;

	PeakAssignmentResidueMap const& assignments() const {
		return *assignments_;
	}

	PeakAssignmentResidueMap& assignments() {
		return *assignments_;
	}

	//   ResonanceList const& resonances() const {
	//     return *resonances_;
	//   }

	core::Size count_assignments() const;
	void delete_diagonal_peaks();
	void update_peak_volumina();
	void network_analysis();// ResonanceList const& );
	void set_trivial_decoy_compatibility_score();
	CrossPeaks const& peaks() const { return peaks_; }
	const_iterator begin() const { return peaks_.begin(); }
	const_iterator end() const { return peaks_.end(); }
	iterator begin() { return peaks_.begin(); }
	core::Size size() const { return peaks_.size(); }
private:
	void update_assignment_list(); //after find_assignments or delete_diagonal_peaks()

	/// @brief return average upper distance bound
	core::Real calibrate( PeakCalibrator const& calibrator );
	CrossPeaks peaks_;
	//probably good to have all things in this class?
	//  ResonanceListOP resonances_;
	PeakAssignmentResidueMapOP assignments_;
};

}
}

#endif
