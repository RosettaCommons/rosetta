// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file OutputMovers.cc
/// @brief File to contain classes that deal with output and pdb dumping
/// @details
/// @author Monica Berrondo

#include <protocols/moves/OutputMovers.hh>

// Rosetta Headers
#include <core/pose/Pose.hh> // REMOVE THIS ASAP! -- this have very obviously been cut and pasted...

#include <protocols/moves/OutputMovers.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>

#include <basic/Tracer.hh>
#include <basic/prof.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ Headers
#include <map>
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace moves {

using basic::T;
using basic::Error;
using basic::Warning;

//constructor
PDBDumpMover::PDBDumpMover(std::string name_in)
	 : Mover(), name_(name_in), num_(0)
{
		//tracer_.init(tr);

		Mover::type("PDBDump");
}

//destructor
PDBDumpMover::~PDBDumpMover() {}

void PDBDumpMover::apply( core::pose::Pose & pose )
{
		num_+=1;
		std::string filename ( name_+ObjexxFCL::right_string_of(num_,2,'0')+".pdb" );
		//core::io::pdb::traced_dump_pdb( tracer_, pose, filename );
		core::io::pdb::dump_pdb( pose, filename );
}

std::string
PDBDumpMover::get_name() const {
	return "PDBDumpMover";
}

void PDBDumpMover::name( std::string name_in )
{
		name_ = name_in;
	  clear();
}

void PDBDumpMover::clear() { num_=0; }

//constructor
ProfilerMover::ProfilerMover() :
	Mover()
{
		Mover::type( "Profiler" );
}

//destructor
ProfilerMover::~ProfilerMover() {}

void ProfilerMover::apply( core::pose::Pose & /*pose*/ )
{
		basic::prof_show();
}

std::string
ProfilerMover::get_name() const {
	return "ProfilerMover";
}

//constructor
MCShowMover::MCShowMover( MonteCarloOP mc_in ) :
	Mover(), mc_( mc_in )
{
		Mover::type( "MCShow" );
}

//destructor
MCShowMover::~MCShowMover() {}

void MCShowMover::apply( core::pose::Pose & pose )
{
		using namespace ObjexxFCL::format;
		mc_->show_scores();
		mc_->score_function()( pose );
		/// Now handled automatically.  mc_->score_function().accumulate_residue_total_energies( pose );
		for ( core::scoring::EnergyMap::const_iterator it=pose.energies().total_energies().begin(),
						it_end = pose.energies().total_energies().end(); it != it_end; ++it ) {
			if ( *it != 0.0 ) {
				std::cout << "total_energy " << core::scoring::ScoreType( it - pose.energies().total_energies().begin() + 1 ) << ' ' << F(12,3,*it) << std::endl;
			}
		}
}

std::string
MCShowMover::get_name() const {
	return "MCShowMover";
}

} // moves
} // protocols
