// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/RotamerSetOperations/SpecialRotamerRotSetOps.hh
/// @brief  classes for rigid body movement during rotamer packing
/// @author Florian Richter, floric@u.washington.edu, sep 2009

#ifndef INCLUDED_protocols_toolbox_rotamer_set_operations_SpecialRotamerRotSetOps_hh
#define INCLUDED_protocols_toolbox_rotamer_set_operations_SpecialRotamerRotSetOps_hh

// Unit Headers
#include <protocols/toolbox/rotamer_set_operations/SpecialRotamerRotSetOps.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.hh>

//Project headers
#include <core/conformation/Residue.fwd.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>

#ifdef WIN32
#include <core/graph/Graph.hh>
#endif


namespace protocols {
namespace toolbox {
namespace rotamer_set_operations {

class SpecialRotamerRSO : public core::pack::rotamer_set::RotamerSetOperation
{
public:
	typedef core::pack::rotamer_set::RotamerSetOperation parent;
	typedef core::Real Real;
	typedef core::Size Size;

	SpecialRotamerRSO( core::Size seqpos );
	SpecialRotamerRSO( SpecialRotamerRSO const & src );
	~SpecialRotamerRSO();

	virtual
	core::pack::rotamer_set::RotamerSetOperationOP
	clone() const;

	virtual
	void
	alter_rotamer_set(
		core::pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		core::pack::task::PackerTask const & ptask,
		core::graph::GraphCOP packer_neighbor_graph,
		core::pack::rotamer_set::RotamerSet & rotamer_set
	);

	virtual //virtual or not?
	void
	set_new_rots(
		core::pack::rotamer_set::Rotamers & new_rots
	);

private:
	core::Size seqpos_;
	core::pack::rotamer_set::Rotamers new_rots_;

};

} //namespace rotamer_set_operations
} //namespace toolbox
} //namespace protocols

#endif
