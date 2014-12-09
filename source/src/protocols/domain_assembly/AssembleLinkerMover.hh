// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/domain_assembly/AssembleLinkerMover.fwd.hh
/// @brief  A mover for assembling domains
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_domain_assembly_AssembleLinkerMover_hh
#define INCLUDED_protocols_domain_assembly_AssembleLinkerMover_hh

#include <protocols/domain_assembly/AssembleLinkerMover.fwd.hh>

#include <protocols/moves/Mover.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>

#include <core/fragment/FragSet.hh>
#include <utility/vector1.hh>
#include <core/types.hh>

#include <string>

namespace protocols {
namespace domain_assembly {

class AssembleLinkerMover : public protocols::moves::Mover {
public:

	AssembleLinkerMover(
		std::string const & loop_mover_name,
		Size const min_loop_size,
		utility::vector1< core::fragment::FragSetOP > frags
	) :
		Mover("AssembleLinkerMover"),
		loop_mover_name_(loop_mover_name),
		min_loop_size_(min_loop_size),
		frag_libs_(frags)
	{}

	virtual std::string get_name() const {
		return "AssembleLinkerMover";
	}

	void apply( core::pose::Pose & pose );
private:
	// fragments
	std::string const loop_mover_name_;
	core::Size const min_loop_size_;
	utility::vector1< core::fragment::FragSetOP > frag_libs_;
};

}//namespace domain_assembly
}//namespace protocols

#endif

