// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/precomputed/PrecomputedLibraryMover.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_modeler_precomputed_PrecomputedLibraryMover_HH
#define INCLUDED_protocols_stepwise_modeler_precomputed_PrecomputedLibraryMover_HH

#include <protocols/moves/Mover.hh>
#include <protocols/stepwise/modeler/precomputed/PrecomputedLibraryMover.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace stepwise {
namespace modeler {
namespace precomputed {

class PrecomputedLibraryMover: public protocols::moves::Mover {

public:

	//constructor
	PrecomputedLibraryMover();

	//destructor
	~PrecomputedLibraryMover();

public:

	/// @brief Apply the minimizer to one pose
	using protocols::moves::Mover::apply;
	virtual void apply( core::pose::Pose & pose_to_visualize );

	void apply( core::pose::Pose & pose_to_visualize ) const;

	virtual std::string get_name() const { return "PrecomputedLibraryMover"; }

	bool
	has_precomputed_move( core::pose::Pose const & pose ) const;

private:
	void
	initialize_from_directory( std::string const & dir_name );

private:

	std::map< std::string, core::io::silent::SilentFileDataOP > library_map_;

};

} //precomputed
} //modeler
} //stepwise
} //protocols

#endif
