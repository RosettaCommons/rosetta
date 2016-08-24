// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    helix_from_sequence.cc
/// @details Creates a helix pose from the sequence only, works for soluble and
///    membrane proteins; the latter will be transformed into the membrane
/// @author  JKLeman (julia.koehler1982@gmail.com)

// App headers
#include <devel/init.hh>

// Project Headers
#include <core/conformation/membrane/SpanningTopology.hh>

#include <protocols/moves/Mover.hh>

#include <protocols/membrane/HelixFromSequence.hh>

#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/OptimizeProteinEmbeddingMover.hh>
#include <protocols/membrane/TransformIntoMembraneMover.hh>
#include <protocols/relax/membrane/MPRangeRelaxMover.hh>
#include <protocols/simple_moves/ScoreMover.hh>

#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>

// Package Headers
#include <apps/benchmark/performance/init_util.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/types.hh>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>
#include <string>

using basic::Error;
using basic::Warning;

using namespace basic::options;
using namespace core;
using namespace core::pose;
using namespace core::conformation::membrane;
using namespace protocols::membrane;

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.jkleman.helix_from_sequence" );

////////////////////////////////////////////////////////////////////////////

// read the fasta file
std::string read_fasta() {

	std::string sequence;

	if ( option[OptionKeys::in::file::fasta].user() ) {

		sequence = core::sequence::read_fasta_file( option[ OptionKeys::in::file::fasta ]()[1] )[1]->sequence();
		TR << "Read in fasta file" << option[OptionKeys::in::file::fasta]()[1] << std::endl;
	} else {
		throw new utility::excn::EXCN_Msg_Exception("Please provide fasta file!");
	}

	return sequence;

} // read fasta file

////////////////////////////////////////////////////////////////////////////

// transform pose into membrane if it is a membrane pose
bool read_membrane() {

	bool membrane( false );

	if ( option[OptionKeys::mp::setup::transform_into_membrane].user() ) {
		membrane = option[OptionKeys::mp::setup::transform_into_membrane]();

		if ( membrane == true ) {
			TR << "Pose is a membrane protein and will be transformed into the membrane" <<  std::endl;
		}
	}

	return membrane;

} // read membrane

////////////////////////////////////////////////////////////////////////////

// optimize embedding?
bool optimize_embedding() {

	bool optimize( false );

	if ( option[OptionKeys::mp::transform::optimize_embedding].user() ) {
		optimize = option[OptionKeys::mp::transform::optimize_embedding]();

		if ( optimize == true ) {
			TR << "Protein embedding in the membrane will be optimized" <<  std::endl;
		}
	}

	return optimize;

} // optimize embedding?

////////////////////////////////////////////////////////////////////////////

// skip relax?
bool skipping_relax() {

	bool skip_relax( false );

	if ( option[OptionKeys::relax::range::skip_relax].user() ) {
		skip_relax = option[OptionKeys::relax::range::skip_relax]();

		if ( skip_relax == true ) {
			TR << "Skipping refinement step..." << std::endl;
		}
	}

	return skip_relax;

} // skip relax?

////////////////////////////////////////////////////////////////////////////

// create an ideal helix from the fasta and dump it
void helix_from_sequence() {

	using namespace protocols::relax::membrane;
	using namespace protocols::simple_moves;

	// read fasta file
	std::string seq = read_fasta();

	// get input options
	bool mem = read_membrane();
	bool opt = optimize_embedding();
	bool skip_relax = skipping_relax();

	// initialize pose
	core::pose::Pose pose;

	// create ideal helix pose from the sequence:
	// 1. create pose from sequence
	make_pose_from_sequence( pose, seq, core::chemical::FA_STANDARD );
	TR << "pose:total_residue: " << pose.total_residue() << std::endl;

	// 2. need to set the PDBInfo object in the pose, because
	// make_pose_from_sequence doesn't take care of that!
	PDBInfoOP pdbinfo( new PDBInfo( pose ) );
	pose.pdb_info( pdbinfo );

	// 3. create ideal helices from SSEs
	for ( Size i = 1; i <= pose.total_residue(); ++i ) {
		pose.set_phi(   i, -62 );
		pose.set_psi(   i, -41 );
		pose.set_omega( i, 180 );
	}

	// if a membrane protein: transform helix into membrane:
	if ( mem == true ) {

		// create topology object from first to last residue
		SpanningTopologyOP topo( new SpanningTopology() );
		topo->add_span( 1, pose.total_residue() );

		if ( opt == true ) {

			// run AddMembraneMover with topology
			AddMembraneMoverOP addmem( new AddMembraneMover( topo, 1, 0 ));
			addmem->apply( pose );

			topo->show();

			// transform into membrane and optimize embedding
			// runs TransformIntoMembrane underneath
			OptimizeProteinEmbeddingMoverOP opt( new OptimizeProteinEmbeddingMover() );
			opt->apply( pose );

			// run MPRangeRelax with default values (nres and 0.1 angle max)
			if ( ! skip_relax ) {
				MPRangeRelaxMoverOP relax( new MPRangeRelaxMover() );
				relax->apply( pose );
			}

		} else {

			// run AddMembraneMover with topology
			AddMembraneMoverOP addmem( new AddMembraneMover( topo, 1, 0 ));
			addmem->apply( pose );

			// transform into membrane
			TransformIntoMembraneMoverOP transform( new TransformIntoMembraneMover() );
			transform->apply( pose );

			// run MPRangeRelax with default values (nres and 0.1 angle max)
			if ( ! skip_relax ) {
				MPRangeRelaxMoverOP relax( new MPRangeRelaxMover() );
				relax->apply( pose );
			}
		}
	}

	// score pose for scorefile
	ScoreMoverOP score( new ScoreMover( "mpframework_smooth_fa_2012.wts" ) );
	score->apply( pose );

	// dump PDB
	pose.dump_pdb("helix_from_sequence.pdb");

} // helix_from_sequence

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// MAIN ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {

		using namespace protocols::jd2;
		using namespace protocols::membrane;

		// initialize option system, RNG, and all factory-registrators
		devel::init(argc, argv);

		HelixFromSequenceOP helix( new HelixFromSequence() );
		JobDistributor::get_instance()->go( helix );

		// call my function
		//  helix_from_sequence();

	}
catch ( utility::excn::EXCN_Base const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}
}
