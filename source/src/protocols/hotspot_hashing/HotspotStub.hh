// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author John Karanicolas, Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_hotspot_hashing_HotspotStub_hh
#define INCLUDED_protocols_hotspot_hashing_HotspotStub_hh

// Unit headers
#include <protocols/hotspot_hashing/HotspotStub.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <protocols/filters/Filter.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace hotspot_hashing {

enum StubStatus {
	reject=1, // stub is not compatible with scaffold position[i]
	accept,  // stub is compatible with scaffold position[i]
	unchecked // stub has not been checked yet
};

class HotspotStub : public utility::pointer::ReferenceCount {
public:
	/// @brief no-argument constructor REQUIRED FOR WINDOWS
	HotspotStub() {}

	//SJF do we need to instantiate the filter at construction? At that point, we don't actually know what filter we would
	// use...
	HotspotStub( core::conformation::ResidueCOP const & residue,
		core::Real const bonus_value,
		core::pose::PoseOP pose,
		core::Size chain_to_design,
		protocols::filters::FilterCOP filter );

	HotspotStub( HotspotStub const & src );

	virtual ~HotspotStub();

#ifndef BOINC // gives windows build error
#ifndef WIN32 // gives windows build error
	HotspotStub & operator=( HotspotStub const & src );
	bool operator<( HotspotStub const & right ) const ;
#endif
#endif

	/// @brief Return this potential hotspot's bonus value
	core::Real bonus_value() const;

	/// @brief Return the residue associated with this potential hotspot
	core::conformation::ResidueCOP residue() const;

	/// @brief Set status at position
	void set_scaffold_status( core::Size const seqpos, StubStatus const status );

	/// @brief Get status at position (setting status if not set)
	bool get_scaffold_status( core::Size const seqpos );

	void set_filter( protocols::filters::FilterCOP filter );

	/// @brief find the residue that's closest to the stub on chain_to_design_
	core::Size get_nearest_residue( core::pose::Pose const & pose ) const;

private:
	friend class HotspotStubSet;
	core::conformation::ResidueCOP residue_;
	core::Real bonus_value_;
	core::pose::PoseOP pose_; // SJF DANGEROUS to make this an OP, but copying the pose for each operation is too costly
	protocols::filters::FilterCOP filter_;
	core::Size chain_to_design_;

	// 0 = stub is not compatible with scaffold position[i]
	// 1 = stub is compatible with scaffold position[i]
	// 2 = stub has not been checked yet
	// nonconst, since it can be modified by the containing stubset
	std::vector<StubStatus> scaffold_status_;

	/// @brief Associate this stub with its set's scaffold for design
	void pair_with_scaffold( );
	void pair_with_scaffold( core::pose::PoseOP pose, protocols::filters::FilterCOP filter, core::Size chain_to_design );

	/// @brief Check the stub's match to scaffold position
	bool scaffold_match( core::Size const seqpos );
};

} // namespace hotspot_hashing
} // namespace protocols

#endif
