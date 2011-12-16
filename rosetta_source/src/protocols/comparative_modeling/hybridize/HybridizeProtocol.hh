// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Add constraints to the current pose conformation.
/// @author Yifan Song

#ifndef INCLUDED_protocols_moves_HybridizeProtocol_hh
#define INCLUDED_protocols_moves_HybridizeProtocol_hh

#include <protocols/moves/Mover.hh>

#include <core/pack/task/TaskFactory.fwd.hh>


namespace protocols {
namespace moves {

class HybridizeProtocol : public Mover {

public:
	HybridizeProtocol();
	virtual ~HybridizeProtocol();
		
	virtual void apply( Pose & );
	virtual std::string get_name() const;

	virtual MoverOP clone() const;
	virtual MoverOP fresh_instance() const;
	
	//virtual void
	//parse_my_tag( TagPtr const, DataMap &, Filters_map const &, Movers_map const &, Pose const & );

	
private:
	utility::vector1 < core::pose::PoseCOP > templates_;
	utility::vector1 < protocols::loops::Loops > template_chunks_;
	utility::vector1 < protocols::loops::Loops > template_contigs_;
	core::fragment::FragSetOP fragments9_, fragments3_; // abinitio frag9,frag3 flags

	// native pose, aln
	core::pose::PoseOP native_;
	core::sequence::SequenceAlignmentOP aln_;
};

} // moves
} // protocols

#endif
