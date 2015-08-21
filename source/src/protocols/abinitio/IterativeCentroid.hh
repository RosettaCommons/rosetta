// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file AbrelaxMover
/// @brief  this class will be handled to a SampleProtocol as a control instance
/// @details responsibilities:
///           know which chainbreaks to penalize and close
///           know which jumps to use during sampling, which (if any) to keep after loop-closing
///           supply a JumpMover if jumps should be moved
///           supply a MoveMap
///           supply a "StrictMoveMap": the protocol should not move anything that is dissallowed in strict_movemap(),
///                      it should try to move just stuff in movemap()
/// should this class also know how to ramp score terms ?
/// handle the titration of constraints ?
/// @author Oliver Lange


#ifndef INCLUDED_protocols_abinitio_IterativeCentroid_hh
#define INCLUDED_protocols_abinitio_IterativeCentroid_hh

// Unit Headers
//#include <protocols/abinitio/IterativeCentroid.fwd.hh>

// Package Headers
#include <protocols/abinitio/IterativeBase.hh>
#include <protocols/abinitio/IterativeFullatom.hh>

// Project Headers
#include <protocols/jd2/archive/ArchiveManager.fwd.hh>

#include <core/io/silent/SilentStruct.hh>

// ObjexxFCL Headers

// Utility headers

//// C++ headers
//#include <cstdlib>
//#include <string>

#include <utility/vector1.hh>

//Auto Headers


namespace protocols {
namespace abinitio {

class IterativeCentroid : public IterativeBase {
	typedef IterativeBase Parent;
public:
	IterativeCentroid(IterativeFullatom* fullatom_pool_ptr ) :
		IterativeBase( "centroid_pool" ),
		fullatom_pool_ptr_( fullatom_pool_ptr ) {};

	virtual void gen_diversity_pool( jd2::archive::Batch& batch, bool fullatom = false );

	virtual void update_noesy_filter_files(
		std::string const& current,
		bool fullatom
	);

	/// @brief save and restore archive to file-system
	virtual void save_to_file( std::string suffix = "" );
	virtual bool restore_from_file();

protected:


	virtual void erase_decoy(
		std::string const& tag
	);

	/// @brief call to insert structure at position given by iterator
	virtual void add_structure_at_position (
		SilentStructs::iterator iss,
		core::io::silent::SilentStructOP new_decoy,
		core::io::silent::SilentStructOP alternative_decoy
	);

	virtual void collect_alternative_decoys( SilentStructs primary_decoys, std::string alternative_decoy_file, SilentStructVector& output_decoys );

private:
	/// @brief also have to keep the stage2 decoys for the stage2 resampling (stages IV and VI)
	SilentStructs stage2_decoys_;

	IterativeFullatom* fullatom_pool_ptr_;
};


}
}

#endif
