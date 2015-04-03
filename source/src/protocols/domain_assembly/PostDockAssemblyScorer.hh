// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/domain_assembly/PostDockAssemblyScorer.fwd.hh
/// @brief  Computes crmsd of the assembly to the staring point
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_domain_assembly_PostDockAssemblyScorer_hh
#define INCLUDED_protocols_domain_assembly_PostDockAssemblyScorer_hh

#include <protocols/domain_assembly/PostDockAssemblyScorer.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace domain_assembly {

class PostDockAssemblyScorer : public protocols::moves::Mover {
public:

	PostDockAssemblyScorer() :
		Mover("PostDockAssemblyScorer"),
		score_prefix_("rebuild_dist")
	{}
	PostDockAssemblyScorer( std::string const & prefix ) :
		Mover("PostDockAssemblyScorer"),
		score_prefix_(prefix)
	{}

	virtual std::string get_name() const {
		return "PostDockAssemblyScorer";
	}

	void apply( core::pose::Pose & pose );

private:
	std::string const score_prefix_;
};


}//namespace domain_assembly
}//namespace protocols

#endif

