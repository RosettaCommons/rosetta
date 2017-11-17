// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/apps/pilot/ralford/sample_tilt_angles.cc
/// @brief Calculates the energy at all possible tilt angles (0->180 degrees)
/// @author Rebecca Alford (rfalford12@gmail.com)

// Note - going to do this with three energy functions... will make more general later (TODO)
// Note - it might be worthwhile to also add a low resolution component to this analysis

// Unit Headers
#include <devel/init.hh>

// Package Headers
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/util.hh>

#include <protocols/membrane/TransformIntoMembraneMover.hh>
#include <protocols/membrane/geometry/EmbeddingDef.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>

// Project Headers
#include <protocols/rigid/RigidBodyMover.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/Job.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/Conformation.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/types.hh>

// Utility Headers
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>

#include <numeric/xyzVector.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>
#include <ostream>
#include <fstream>

static basic::Tracer TR( "apps.pilot.ralford.SampleTiltProtocol" );

class SampleTiltProtocol : public protocols::moves::Mover {

public:

	SampleTiltProtocol()
	{
		register_options();
		init_options();
	}

	virtual ~SampleTiltProtocol() {};

	virtual
	void
	apply( core::pose::Pose & pose ) {

		using namespace core;
		using namespace core::scoring;
		using namespace protocols::jd2;
		using namespace protocols::membrane;
		using namespace protocols::membrane::geometry;
		using namespace protocols::rigid;
		using namespace numeric;

		TR << "Adding a membrane residue to the pose at the residue COM" << std::endl;
		// Initialize the membrane framework via AddMembraneMover
		// Anchor the membrane jump to the residue center of mass
		core::Size rsd_com( residue_center_of_mass( pose, 1, pose.size()-1 ) );
		AddMembraneMoverOP add_memb( new AddMembraneMover( rsd_com, 0 ) );
		add_memb->apply( pose );

		// Determine the membrane jump number (pre-determined in the framework)
		Size membrane_jump( pose.conformation().membrane_info()->membrane_jump() );
		Vector center( pose.conformation().membrane_info()->membrane_center( pose.conformation() ) );
		TR << "Setting up an initial membrane orientation by aligning the peptide to the current membrane normal" << std::endl;
		// Calculate the current helix axes based on the COM of the transmembrane
		// span start and end posiitons
		numeric::xyzVector< core::Real > helix_axis = calc_helix_axis( pose, 1 /* only for single span poses */ );

		// Calculate the angle between the current helix axis and the membrane normal
		Vector spin_axis(1, 0, 0);

		// Get ready to score
		ScoreFunctionOP talaris = get_score_function();
		ScoreFunctionOP membrane07 = ScoreFunctionFactory::create_score_function( "mpframework_fa_2007" );
		ScoreFunctionOP membrane12 = ScoreFunctionFactory::create_score_function( "mpframework_smooth_fa_2012" );

		// Rotate the helix to align with the membrane normal (representing a tilt angle of 0 degrees)
		EmbeddingDefOP embedding = compute_structure_based_embedding( pose );
		embedding->set_normal( helix_axis );

		// Transform the helix into a membrane coordinate frame such that the helix axis is parallel
		// to the membrane normal.
		TransformIntoMembraneMoverOP transform_into_memb( new TransformIntoMembraneMover( embedding ) );
		transform_into_memb->apply( pose );

		// Intermediate checks!
		TR << calc_helix_tilt_angle( pose, 1 /* only working with 1 helix poses for now */ ) << std::endl;

		core::Real talaris_score( talaris->score( pose ) );
		core::Real membrane07_score( membrane07->score( pose ) );
		core::Real membrane12_score( membrane12->score( pose ) );

		// Setup score & angle arrays
		utility::vector1< core::Real > angles;
		utility::vector1< core::Real > talaris_scores;
		utility::vector1< core::Real > membrane07_scores;
		utility::vector1< core::Real > membrane12_scores;

		// Add the first scores and angles
		angles.push_back( 0 ); // angle should be zero or close to it
		talaris_scores.push_back( talaris_score );
		membrane07_scores.push_back( membrane07_score );
		membrane12_scores.push_back( membrane12_score );

		// Reset the deterministic step size to 1
		// Rotate the helix to align with the membrane normal (representing a tilt angle of 0 degrees)
		core::Real step_size( 1.0 );
		RigidBodyDeterministicSpinMoverOP deterministic_z( new RigidBodyDeterministicSpinMover( membrane_jump, spin_axis, center, step_size ) );

		// For each tilt angle (sample every degree for a smooth surface), score the protein and add to arrays
		for ( core::Size angle = 1; angle <= 180; ++angle ) {

			deterministic_z->apply( pose );

			// Intermediate checks!
			TR << calc_helix_tilt_angle( pose, 1 /* only working with 1 helix poses for now */ ) << std::endl;

			// Get the current state of the pose
			talaris_score = talaris->score( pose );
			membrane07_score = membrane07->score( pose );
			membrane12_score = membrane12->score( pose );

			// Add state to arrays
			angles.push_back( angle );
			talaris_scores.push_back( talaris_score );
			membrane07_scores.push_back( membrane07_score );
			membrane12_scores.push_back( membrane12_score );
		}

		// Write scores to output
		write_score_to_outfiles( angles, talaris_scores, membrane07_scores, membrane12_scores );

		TR << "Sample Tilt Anlge Protocol Complete!" << std::endl;
	}

	void
	write_score_to_outfiles(
		utility::vector1< core::Real > angles,
		utility::vector1< core::Real > talaris_scores,
		utility::vector1< core::Real > membrane07_scores,
		utility::vector1< core::Real > membrane12_scores
	) {

		using namespace std;

		// Open three output files for writing - one for each energy function type
		std::ofstream talaris_out( "M13-coat_talaris.out" );
		std::ofstream membrane07_out( "M13-coat_membrane07.out" );
		std::ofstream membrane12_out( "M13-coat_membrane12.out" );

		// Write the score/angle pairs to an output file
		talaris_out << "TiltAngle Score" << std::endl;
		membrane07_out << "TiltAngle Score" << std::endl;
		membrane12_out << "TiltAngle Score" << std::endl;

		// Check that angles and score lists are of equal length
		assert( angles.size() == talaris_scores.size() );
		assert( angles.size() == membrane07_scores.size() );
		assert( angles.size() == membrane12_scores.size() );

		for ( core::Size i = 1; i <= angles.size(); ++i ) {
			talaris_out << angles[i] << " " << talaris_scores[i] << std::endl;
			membrane07_out << angles[i] << " " << membrane07_scores[i] << std::endl;
			membrane12_out << angles[i] << " " << membrane12_scores[i] << std::endl;
		}

		talaris_out.close();
		membrane07_out.close();
		membrane12_out.close();
	}

	std::string get_name() const { return "SampleTiltProtocol"; }

	void
	register_options() {

		using namespace basic::options;

		option.add_relevant( OptionKeys::score::weights );
		option.add_relevant( OptionKeys::mp::benchmark::tilt_angle::output );

	}

	void
	init_options() {

		using namespace basic::options;

		// Read score function weights
		/**
		if ( option[ OptionKeys::score::weights ].user() ) {
		sfxn_weights_ = option[ OptionKeys::score::weights ]();
		}
		**/

		// Read user specified output file
		if ( option[ OptionKeys::mp::benchmark::tilt_angle::output ].user() ) {
			outfile_ = option[ OptionKeys::mp::benchmark::tilt_angle::output ]();
		}

	}

private:

	std::string outfile_;
	//std::string sfxn_weights_;

};


typedef utility::pointer::shared_ptr< SampleTiltProtocol > SampleTiltProtocolOP;
typedef utility::pointer::shared_ptr< SampleTiltProtocol const  > SampleTiltProtocolCOP;

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try{

		devel::init(argc, argv);
		SampleTiltProtocolOP sample_tilt( new SampleTiltProtocol );
		protocols::jd2::JobDistributor::get_instance()->go(sample_tilt);

	} catch (utility::excn::Exception const & e ) {
		std::cout << "Caught Exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
