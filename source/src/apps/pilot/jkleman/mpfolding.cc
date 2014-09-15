// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    mpfolding.cc
/// @details last Modified: 4/4/14
/// @author  JKLeman (julia.koehler1982@gmail.com)

// App headers
#include <devel/init.hh>

// Project Headers
#include <protocols/moves/Mover.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>


// Package Headers
#include <apps/benchmark/performance/init_util.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <basic/Tracer.hh>

#include <protocols/viewer/viewers.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>


///////////////////////////////////////////////////////////////////////////////
// HEADERS FROM PROTOCOL

	// Unit Headers
#include <protocols/moves/Mover.hh>

	// Project Headers
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/docking/DockMCMProtocol.hh>
#include <protocols/docking/DockingProtocol.hh>
#include <protocols/moves/MoverContainer.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMolMover.hh>

	// Package Headers
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/membrane/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/sequence/util.hh>
#include <core/types.hh>

	// Utility Headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/membrane_new.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

	// C++ Headers
#include <cstdlib>

using basic::Error;
using basic::Warning;

using namespace core;
using namespace core::pose;
using namespace core::conformation;
using namespace core::conformation::membrane;

static basic::Tracer TR( "apps.pilot.jkleman.MPFolding" );

////////////////////////////////////////////////////////////////////////////////
//////////////////////////// HEADER ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

class MPFoldingMover : public protocols::moves::Mover {
	
public:
	
	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default Constructor
	MPFoldingMover();
	
	/// @brief Copy Constructor
	MPFoldingMover( MPFoldingMover const & src );
	
	/// @brief Destructor
	virtual ~MPFoldingMover();
	
public: // methods
	
	/// @brief Create a Clone of this mover
	virtual protocols::moves::MoverOP clone() const;
	
	/// @brief Create a Fresh Instance of this Mover
	virtual protocols::moves::MoverOP fresh_instance() const;
	
	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Get the name of this Mover (MPFoldingMover)
	virtual std::string get_name() const;
	
	/// @brief Fold MP
	virtual void apply( Pose & pose );
	
private: // methods
	
	// setup the foldtree with cutpoints according to topology
//	void setup_foldtree();
	
	
	
private: // data
	
	// topology object storing SSEs
	SpanningTopologyOP SSE_topo_;
	
	// topology object storing loops
	SpanningTopologyOP loops_;
	
	// foldtree
	FoldTreeOP foldtree_;
	
	// add membrane mover
	protocols::membrane::AddMembraneMoverOP add_membrane_mover_;
	
	// docking protocols protocol
	protocols::docking::DockingProtocolOP docking_protocol_;
	
	// sequence mover
	protocols::moves::RandomMoverOP random_mover_;
	
	// scorefunction
	core::scoring::ScoreFunctionOP lowres_scorefxn_;
	core::scoring::ScoreFunctionOP highres_scorefxn_;
	
	// kT for MCM protocol
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// Real kT_;
	
};

////////////////////////////////////////////////////////////////////////////////
//////////////////////////// IMPLEMENTATION/////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default Constructor
/// @details Create a membrane pose setting the membrane center
/// at center=(0, 0, 0), normal=(0, 0, 1) and loads in spans
/// and lips from the command line interface.
MPFoldingMover::MPFoldingMover() :
protocols::moves::Mover()
{}

/// @brief Copy Constructor
/// @details Create a deep copy of this mover
MPFoldingMover::MPFoldingMover( MPFoldingMover const & src ) :
protocols::moves::Mover( src )
{}

/// @brief Destructor
MPFoldingMover::~MPFoldingMover() {}

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
MPFoldingMover::clone() const {
	return ( new MPFoldingMover( *this ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
MPFoldingMover::fresh_instance() const {
	return new MPFoldingMover();
}

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this Mover (MPFoldingMover)
std::string
MPFoldingMover::get_name() const {
	return "MPFoldingMover";
}


/// @brief Add Membrane Components to Pose
void
MPFoldingMover::apply( Pose & pose ) {
	
	using namespace basic;
	using namespace basic::options;
	using namespace core;
	using namespace core::pose;
	using namespace core::sequence;
	using namespace core::conformation::membrane;

	// register options
	option.add_relevant( OptionKeys::in::file::fasta );

	// check if fasta there
	if ( ! option[OptionKeys::in::file::fasta].user() ){
		utility_exit_with_message("Please provide a fasta file!");
	}

	// read fasta
	SequenceOP sequence = read_fasta_file(
					option[ OptionKeys::in::file::fasta ]()[1] )[1];

	// get number of residues
	std::string seq = sequence->to_string();
	//Size nres = sequence->length();

	// create pose from sequence
//	Pose pose;
	make_pose_from_sequence( pose, seq, "centroid" );
	
	// create topology from spanfile
	SpanningTopologyOP topo = new SpanningTopology( spanfile_name() );

	
	// create ideal helices from SSEs
	
	
	// setup foldtree
	// 1) add membrane at root
	

	// 2) add jumps to all COMs of SSEs


	// dock each SSE starting from center of sequence
	
	
	// filter via loop length constraint
	
	
	// consider other regular constraints
	
	
	// switch fragments and minimize































	
}

////////////////////////////////////////////////////////////////////////////////

// read spanfiles
std::string read_spanfile(){
	
	using namespace basic;
	using namespace basic::options;

	// cry if spanfiles not given
	if ( ! option[OptionKeys::membrane_new::setup::spanfiles].user() ){
		utility_exit_with_message("Please provide a single spanfiles!");
	}
	
	// get filenames from Optionsystem
	std::string spanfile = option[OptionKeys::membrane_new::setup::spanfiles]()[1];
	
	return spanfile;
}// read spanfile


////////////////////////////////////////////////////////////////////////////////

typedef utility::pointer::owning_ptr< MPFoldingMover > MPFoldingMoverOP;

/////////////////////////////////////// MAIN ///////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {
//		using namespace protocols::docking::membrane;
		using namespace protocols::jd2;
		
		// initialize option system, random number generators, and all factory-registrators
		devel::init(argc, argv);
		//protocols::init(argc, argv);
		
		MPFoldingMoverOP mpfm = new MPFoldingMover();
		JobDistributor::get_instance()->go(mpfm);
	}
	catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
