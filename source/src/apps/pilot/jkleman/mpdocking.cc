// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    mpdocking.cc
/// @brief   Dock two membrane proteins in the membrane
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

static basic::Tracer TR( "apps.pilot.jkleman.mpdocking" );

////////////////////////////////////////////////////////////////////////////////
//////////////////////////// HEADER ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

class MPDockingMover : public protocols::moves::Mover {
	
public:
	
	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default Constructor
	/// @details Create a membrane pose setting the membrane center
	/// at center=(0, 0, 0), normal=(0, 0, 1) and loads in spans
	/// and lips from the command line interface.
	MPDockingMover();
	
	/// @brief Copy Constructor
	/// @details Create a deep copy of this mover
	MPDockingMover( MPDockingMover const & src );
	
	/// @brief Destructor
	virtual ~MPDockingMover();
	
public: // methods
	
	/// @brief Create a Clone of this mover
	virtual protocols::moves::MoverOP clone() const;
	
	/// @brief Create a Fresh Instance of this Mover
	virtual protocols::moves::MoverOP fresh_instance() const;
	
	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Get the name of this Mover (MPDockingMover)
	virtual std::string get_name() const;
	
	/// @brief Add Membrane Components to Pose
	/// @details Add membrane components to pose which includes
	///	spanning topology, lips info, embeddings, and a membrane
	/// virtual residue describing the membrane position
	virtual void apply( Pose & pose );
	
private: // methods
	
	// setup
	void setup();
	
private: // data
	
	// add membrane mover
	protocols::membrane::AddMembraneMoverOP add_membrane_mover_;

	// docking protocols protocol
	protocols::docking::DockMCMProtocolOP dock_mcm_protocol_;
	protocols::docking::DockingProtocolOP docking_protocol_;

	// sequence mover
	protocols::moves::RandomMoverOP random_mover_;

	// scorefunction
	core::scoring::ScoreFunctionOP lowres_scorefxn_;
	core::scoring::ScoreFunctionOP highres_scorefxn_;

	// kT
	Real kT_;

	// Membrane Center/Normal pair used for setup
	Vector center_;
	Vector normal_;

	// SpanningTopology
	std::string spanfile_;
	
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
MPDockingMover::MPDockingMover() :
protocols::moves::Mover(),
center_(0, 0, 0),
normal_(0, 0, 1)
{}

/// @brief Copy Constructor
/// @details Create a deep copy of this mover
MPDockingMover::MPDockingMover( MPDockingMover const & src ) :
protocols::moves::Mover( src ),
center_( src.center_ ),
normal_( src.normal_ )
{}

/// @brief Destructor
MPDockingMover::~MPDockingMover() {}

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
MPDockingMover::clone() const {
	return ( new MPDockingMover( *this ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
MPDockingMover::fresh_instance() const {
	return new MPDockingMover();
}

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this Mover (MPDockingMover)
std::string
MPDockingMover::get_name() const {
	return "MPDockingMover";
}


/// @brief Add Membrane Components to Pose
/// @details Add membrane components to pose which includes
///	spanning topology, lips info, embeddings, and a membrane
/// virtual residue describing the membrane position
void
MPDockingMover::apply( Pose & pose ) {
	
	using namespace core::conformation::membrane;
	
	TR << "calling setup" << std::endl;
	setup();
	
	// assuming that protein 1 is fixed in the membrane!!!
	// add membrane VRT, call AddMembraneMover
	TR << "adding MEM" << std::endl;
	add_membrane_mover_->apply( pose );
	
	// foldtree
	TR << "creating foldtree from pose: " << std::endl;
	pose.fold_tree().show(std::cout);
	core::kinematics::FoldTree foldtree = pose.fold_tree();

	// reorder only reorders, but does not rename jump edges
	foldtree.reorder( pose.conformation().membrane()->membrane_rsd_num() );
	pose.fold_tree( foldtree );
	
	// add a jump from protein 1 (fixed) to protein 2 (flexible)
	//	foldtree.add_edge( protein2_start, protein2_end, number );

	// show foldtree
	TR << "foldtree reordered" << std::endl;
	pose.fold_tree().show(std::cout);
	
	// attach Pymol observer
	TR << "test print" << std::endl;
	TR << "attach Pymol observer" << std::endl;
	protocols::moves::AddPyMolObserver( pose );

	// create MC-object
//	TR << "create MC object" << std::endl;
//	protocols::moves::MonteCarloOP montecarlo = new protocols::moves::MonteCarlo( pose, *lowres_scorefxn_, kT_);
		
	// dock MCM protocol - high-res only
//	TR << "calling dock_MCM_protocol" << std::endl;
//	dock_mcm_protocol_->apply( pose );
	
	// regular docking protocol for low-res
	TR << "calling docking protocol" << std::endl;
	docking_protocol_->apply( pose );
	
	// score
//	TR << "calling Metropolis" << std::endl;
//	montecarlo->boltzmann( pose );
	
}

////////////////////////////////////////////////////////////////////////////////

void MPDockingMover::setup(){

	using namespace protocols::membrane;
	using namespace protocols::docking;
	using namespace core::scoring;

	add_membrane_mover_ = new AddMembraneMover();
//	low_res_scorefxn_ = getScoreFunction();

	// create scorefunctions for lowres and highres
//	ScoreFunctionOP lowres_scorefxn_ = ScoreFunctionFactory::create_score_function( "cen_membrane_2014.wts" );
	ScoreFunctionOP lowres_scorefxn_ = ScoreFunctionFactory::create_score_function( "mpdocking_cen_14-6-25.wts" );
	ScoreFunctionOP highres_scorefxn_ = ScoreFunctionFactory::create_score_function( "mpdocking_fa_14-6-25.wts" );
	
	dock_mcm_protocol_ = new DockMCMProtocol( 1, lowres_scorefxn_, highres_scorefxn_);

	// low res only (jump#, bool low-res only, bool local refine only, bool autofoldtree, low scorefxn, high scorefxn)
//	docking_protocol_ = new DockingProtocol( 1, true, false, false, lowres_scorefxn_, highres_scorefxn_ );
	
	// high res only
//	docking_protocol_ = new DockingProtocol( 1, false, true, false, lowres_scorefxn_, highres_scorefxn_ );
	
	// both low res and high res
	docking_protocol_ = new DockingProtocol( 1, true, false, false, lowres_scorefxn_, highres_scorefxn_ );
	
	random_mover_ = new protocols::moves::RandomMover();

	kT_ = 1;
	
}

typedef utility::pointer::owning_ptr< MPDockingMover > MPDockingMoverOP;


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
		
		MPDockingMoverOP mpdm = new MPDockingMover();
		JobDistributor::get_instance()->go(mpdm);
	}
	catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
