// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file
/// @brief
/// @author Ingemar Andre

#ifndef INCLUDED_protocols_moves_symmetry_SymPackRotamersMover_hh
#define INCLUDED_protocols_moves_symmetry_SymPackRotamersMover_hh

// Unit headers
#include <protocols/moves/symmetry/SymPackRotamersMover.fwd.hh>
#include <protocols/moves/PackRotamersMover.hh>

// Project headers
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.fwd.hh>

#include <utility/vector0.hh>

namespace protocols {
namespace moves {
namespace symmetry {

class SymPackRotamersMover : public PackRotamersMover {

public:
	// default constructor
	SymPackRotamersMover();

	SymPackRotamersMover(
		core::scoring::ScoreFunctionCOP scorefxn,
		core::pack::task::PackerTaskCOP task = 0,
		core::Size nloop = 1
	);

	// destructor (important for properly forward-declaring smart-pointer members)
	~SymPackRotamersMover();

	// copy constructor
	SymPackRotamersMover( PackRotamersMover const & other );

//	virtual void apply( core::pose::Pose & pose );

	void
	make_symmetric_task(
		core::pose::Pose & pose,
		core::pack::task::PackerTaskOP task
	);
	virtual std::string get_name() const;

private:

	// to be used/redefined by derived classes
	virtual void setup( core::pose::Pose & pose );
	// need a more elegant rot_to_pack implementation than this
	virtual core::PackerEnergy run(
		core::pose::Pose & pose,
		utility::vector0< int > rot_to_pack = utility::vector0<int>()
	) const;

private:
	// pointers to data that are passed in
	core::pack::rotamer_set::symmetry::SymmetricRotamerSetsOP sym_rotamer_sets_;
	core::pack::task::PackerTaskOP symmetric_task_;
	InteractionGraphBaseOP symmetric_ig_;
};

} // symmetry
} // moves
} // protocols

#endif
