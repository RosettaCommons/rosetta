// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/RotamerSetOperations/SpecialRotamerRotSetOps.hh
/// @brief  class to add rotamers that interact especially well with
///         desired target residues
/// @author Florian Richter, florian.richter.1@hu-berlin.de, may 2014

#ifndef INCLUDED_protocols_toolbox_rotamer_set_operations_AddGood2BPairEnergyRotamers_hh
#define INCLUDED_protocols_toolbox_rotamer_set_operations_AddGood2BPairEnergyRotamers_hh

// Unit Headers
#include <protocols/toolbox/rotamer_set_operations/AddGood2BPairEnergyRotamers.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.hh>

//Project headers
#include <core/conformation/Residue.fwd.hh>
#include <utility/graph/Graph.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>

#ifdef WIN32
#include <utility/graph/Graph.hh>
#endif


namespace protocols {
namespace toolbox {
namespace rotamer_set_operations {

/// @brief implementation of rosetta 2.x style rotamer explosion
/// builds a whole bunch of extra rotamers, but only adds them to the
/// set if they score better than a cutoff with a given set of target
/// residues
class AddGood2BPairEnergyRotamers : public core::pack::rotamer_set::RotamerSetOperation
{
public:
	typedef core::pack::rotamer_set::RotamerSetOperation parent;
	typedef core::Real Real;
	typedef core::Size Size;

	AddGood2BPairEnergyRotamers(
		Size seqpos,
		Size ex_level,
		utility::vector1< core::Size > const & target_seqpos,
		Real score_cutoff,
		bool drop_rots_above_cutoff
	);

	AddGood2BPairEnergyRotamers( AddGood2BPairEnergyRotamers const & src );
	~AddGood2BPairEnergyRotamers();

	virtual
	core::pack::rotamer_set::RotamerSetOperationOP
	clone() const;

	virtual
	void
	alter_rotamer_set(
		core::pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		core::pack::task::PackerTask const & ptask,
		utility::graph::GraphCOP packer_neighbor_graph,
		core::pack::rotamer_set::RotamerSet & rotamer_set
	);


	/// @brief helper function
	Real
	get_res_res_score(
		core::conformation::Residue const & rsd1,
		core::conformation::Residue const & rsd2,
		core::pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn
	) const;

	void
	set_debug( bool setting ){
		debug_ = setting; }

	/// @brief helper function to check if the passed in candi
	bool
	rotamer_set_contains_rotamer(
		core::pack::rotamer_set::RotamerSet const & rotamer_set,
		core::conformation::Residue const & candidte_rot ) const;

private:
	core::Size seqpos_;
	core::Size ex_level_;
	utility::vector1< core::Size > target_seqpos_;
	core::Real score_cutoff_;
	bool drop_rots_above_cutoff_;
	bool debug_;

	bool disabled_;

};

} //namespace rotamer_set_operations
} //namespace toolbox
} //namespace protocols

#endif
