// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/protocols/simple_moves/Rotamerize.hh
/// @brief class declaration for RotamerizeMover
/// @author Jim Havranek

#ifndef INCLUDED_protocols_simple_moves_RotamerizeMover_HH
#define INCLUDED_protocols_simple_moves_RotamerizeMover_HH

// Unit headers
#include <protocols/simple_moves/RotamerizeMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/types.hh>
//#include <core/pack/types.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

//#include <utility/Tag/Tag.fwd.hh>
//#include <utility/vector0.hh>
//#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

/// @brief a mover that replaces the repackable parts - sidechains, bases, etc. - based
/// purely on geometric similarity to the starting structure.  The purpose is to generate
/// the best case output of a repacking calculation as a positive control for benchmarking
/// and parameter fitting.  Thus, the soft rep potential was originally derived as the
/// score function that gave the best fit between repacked and rotamerized structures for
/// a large test set.
class RotamerizeMover : public protocols::moves::Mover {

public:
	typedef core::pack::rotamer_set::RotamerSetsOP RotamerSetsOP;
	typedef core::pack::rotamer_set::RotamerSetsCOP RotamerSetsCOP;
	typedef core::pack::task::PackerTaskOP PackerTaskOP;
	typedef core::pack::task::PackerTaskCOP PackerTaskCOP;
	typedef core::pack::task::TaskFactoryCOP TaskFactoryCOP;

public:
	///@brief default constructor
	RotamerizeMover();
	///@brief constructor with typename
	RotamerizeMover( std::string const & );

	RotamerizeMover(
			PackerTaskOP task
	);

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~RotamerizeMover();

	// copy constructor
	RotamerizeMover( RotamerizeMover const & other );

	// methods
	virtual void apply( Pose & pose );
	virtual std::string get_name() const;
	bool task_is_valid( Pose const & pose ) const; // should this be virtual?

	// setters
	void task_factory( TaskFactoryCOP tf );
	void task( PackerTaskOP t );

	// accessors
	PackerTaskOP task() const;
	TaskFactoryCOP task_factory() const;
	RotamerSetsCOP rotamer_sets() const;

protected:
	///@brief get rotamers, energies. Also performs lazy initialization of ScoreFunction, PackerTask.
	virtual void setup( Pose & pose );

private:
	// pointers to data that are passed in
	PackerTaskOP task_;
	TaskFactoryCOP task_factory_;

	// 'really private:' packer data, actually created and owned by this class
	RotamerSetsOP rotamer_sets_;

};

// note: it is better to create new files, instead of adding additional classes here

} // simple_moves
} // protocols

#endif
