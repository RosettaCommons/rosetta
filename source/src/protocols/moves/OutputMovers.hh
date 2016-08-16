// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file OutputMovers.hh
/// @brief File to contain classes that deal with output and pdb dumping
/// @details
/// @author Monica Berrondo


#ifndef INCLUDED_protocols_moves_OutputMovers_hh
#define INCLUDED_protocols_moves_OutputMovers_hh

// Unit headers

#include <protocols/moves/OutputMovers.fwd.hh>

#include <core/pose/Pose.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>


//#include <basic/Tracer.hh>

// ObjexxFCL Headers

// C++ Headers
#include <map>
#include <string>

#include <utility/vector1.hh>


// Utility Headers

namespace protocols {
namespace moves {

//using basic::T;
//using basic::Error;
//using basic::Warning;

// a mover to dump pdbs within cycles using movers (see DockingHighRes and DockingLowRes for examples on usage
class PDBDumpMover : public Mover
{
public:
	//constructor
	//PDBDumpMover(std::string name_in, basic::Tracer const & tr=core::io::pdb::TR_dump_pdb_dummy);
	PDBDumpMover(std::string name_in);

	//destructor
	~PDBDumpMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	void name( std::string name_in );

	void clear();

private:
	std::string name_;
	int num_;

	//basic::Tracer tracer_;

}; // class PDBDumpMover

// allows something profiler output to be printed to the screen during a move cycle
class ProfilerMover : public Mover
{
public:

	//constructor
	ProfilerMover();

	//destructor
	~ProfilerMover();

	virtual void apply( core::pose::Pose & /* pose*/ );
	virtual std::string get_name() const;
};

// allows mc.show_scores to be used inside a cycle of movers
class MCShowMover : public Mover
{
public:

	//constructors
	MCShowMover( MonteCarloOP mc_in );

	//destructor
	~MCShowMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

private:
	// the monte carlo that needs to be shown
	MonteCarloOP mc_;
};
} // moves
} // protocols

#endif //INCLUDED_protocols_moves_OutputMovers_HH
