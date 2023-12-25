// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/membrane_benchmark/MembraneEnergyLandscapeSampler.cc
/// @brief Sample a protein energy landscape as a function of orientation in the membrane
/// @details Orientation is sampled as a function of three degrees of freedom: helix tilt relative to the membrane
/// normal, azimuthal angle about  its own axis and distance between the membrane center and peptide residue center of mass.
/// @author Rebecca F. Alford (rfalford12@gmail.com)
/// The code is modified on Sept 26th by Rituparna Samanta (rituparna@utexas.edu)
/// An additional degree of freedom is included while sampling the
/// orientation: azimuthal angle beside helix tilt
/// relative to the membrane normal, distance between the membrane center and peptide residue center of mass.
/// Unit headers
///In this version we have first aligned the peptide to the center of the membrane and along the z-axis.
///The definition of center is the average of all CA atoms or the center of the axis joining the upper disk and lower disk.
//It is different from the previous codes, in terms of the center definition and alignment.
//Unit headers
#include <protocols/membrane_benchmark/MembraneEnergyLandscapeSampler.hh>
#include <protocols/membrane_benchmark/MembraneEnergyLandscapeSamplerCreator.hh>

// Package headers
#include <protocols/membrane/TransformIntoMembraneMover.hh>
#include <protocols/membrane/TranslationRotationMover.hh>
#include <protocols/membrane/util.hh>
//#include <protocols/membrane/scoring/MembranepHEnergy.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_moves/MutateResidue.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <core/pose/PDBInfo.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/Energies.hh>
#include <core/energy_methods/pHEnergy.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Jump.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/palette/DefaultPackerPalette.hh>
//#include <core/pack/palette/NoDesignPackerPalette.hh>
//#include <protocols/task_operations/pHVariantTaskOperation.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/Span.hh>

#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/MembranePotential.hh>
#include <core/scoring/MembraneTopology.hh>
#include <core/conformation/util.hh>
#include <core/chemical/AA.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/io/ozstream.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <basic/options/keys/pH.OptionKeys.gen.hh>

#include <numeric/conversions.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/constants.hh>
#include <numeric/HomogeneousTransform.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/random/random.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/down_cast.hh>
#include <utility/vector1.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// C++ Headers
#include <ostream>

#include <protocols/rosetta_scripts/util.hh> // AUTO IWYU For attributes_for_parse_score_...
#include <utility/string_util.hh> // AUTO IWYU For string_split

static basic::Tracer TR( "protocols.membrane_benchmark.MembraneEnergyLandscapeSampler" );

namespace protocols {
namespace membrane_benchmark {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
MembraneEnergyLandscapeSampler::MembraneEnergyLandscapeSampler():
	protocols::moves::Mover( MembraneEnergyLandscapeSampler::mover_name() ),
	sfxn_(),
	sfxn_weights_(),
	rotation_type_( "YZ" ),
	interface_( false ),
	start_z_( -60 ),
	end_z_( 60 ),
	flag_axis_( 1.0 ),
	azimuthal_delta_( 30 ),
	repack_( false ),
	pH_mode_( false )
{
	init_from_options();
}

/// @brief Non-defualt cnstructor for landscape sampling
MembraneEnergyLandscapeSampler::MembraneEnergyLandscapeSampler(
	std::string sfxn_weights,
	std::string rotation_type,
	bool interface,
	core::Real start_z,
	core::Real end_z,
	core::Real flag_axis,
	core::Real azimuthal_delta,
	bool repack ):
	protocols::moves::Mover( MembraneEnergyLandscapeSampler::mover_name() ),
	sfxn_( core::scoring::ScoreFunctionFactory::create_score_function( sfxn_weights ) ),
	sfxn_weights_( sfxn_weights ),
	rotation_type_( rotation_type ),
	interface_( interface ),
	start_z_( start_z ),
	end_z_( end_z ),
	flag_axis_( flag_axis ),
	azimuthal_delta_( azimuthal_delta ),
	repack_( repack ),
	pH_mode_( false )
{

	init_from_options();
}

/// @brief Copy constructor
MembraneEnergyLandscapeSampler::MembraneEnergyLandscapeSampler( MembraneEnergyLandscapeSampler const & src ):
	protocols::moves::Mover( src ),
	sfxn_( src.sfxn_ ),
	sfxn_weights_( src.sfxn_weights_ ),
	rotation_type_( src.rotation_type_ ),
	interface_( src.interface_ ),
	start_z_( src.start_z_ ),
	end_z_( src.end_z_ ),
	flag_axis_( src.flag_axis_ ),
	azimuthal_delta_( src.azimuthal_delta_ ),
	repack_( src.repack_ ),
	pH_mode_( src.pH_mode_ )
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
MembraneEnergyLandscapeSampler::~MembraneEnergyLandscapeSampler() {}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
MembraneEnergyLandscapeSampler::apply( core::pose::Pose & pose ) {

	using namespace numeric;
	using namespace core;
	using namespace protocols;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace protocols::membrane;
	using namespace protocols::membrane::geometry;
	using namespace core::conformation::membrane;
	using namespace numeric::conversions;
	using namespace protocols::minimization_packing;

	// Check that the pose is a membrane protein before continuing
	if ( !pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Pose is not a membrane protein, and I cannot score the membrane energy landscape of a non membrane protein!");
	}

	// Abort if the membrane is not fixed. Leads to significant rounding errors
	if ( !is_membrane_fixed( pose ) ) {
		utility_exit_with_message( "The pose is moveable and the membrane is not fixed. Exiting..." );
	}

	// Perform an initial transformation of the pose into the membrane
	TransformIntoMembraneMoverOP transform_into_memb( utility::pointer::make_shared< TransformIntoMembraneMover >() );
	transform_into_memb->apply( pose );

	/*-----------
	std::string filename8( tempstr1+"_inmem_nopack"+".pdb" );
	pose.dump_pdb( filename8 );
	--------------*/
	core::Vector peptide_center( 0,0,0 );
	peptide_center = getcenter( pose, flag_axis_ );

	TR <<"center after transforming into membrane" << peptide_center.x() << "|" << peptide_center.y() << "|" << peptide_center.z() <<std::endl;

	/*The center matches the center of the membrane only for alpha-helices. Since the com of
	peptides were taken as the average of the starting and ending CA of the span. I am
	replacing here with the com of all alpha atoms.*/

	TR << "Initialize energy landscape sampling with rotation type " << rotation_type_ << std::endl;

	// Setup energy function based on user specified weights file
	TR << "Creating an energy function with weights " << sfxn_weights_ << std::endl;

	// Setup output filename and Header
	TR << "Configuring output file with 3D lanscape data" << std::endl;
	utility::vector1< std::string > temp( utility::string_split( pose.pdb_info()->name(), '/') );
	std::string tempstr = temp[ temp.size() ].substr(0, temp[ temp.size() ].size()-4 );
	std::string filename( tempstr +  "_" + sfxn_->get_name() + "_" + std::to_string( end_z_ ) + "_" + std::to_string( start_z_ ) + "_landscape.dat" );
	utility::io::ozstream output( filename );

	std::string filename2( tempstr +  "_" + sfxn_->get_name() + "_" + std::to_string( end_z_ ) + "_" + std::to_string( start_z_ ) + "_energybreakdown.dat" );
	utility::io::ozstream output2( filename2 );
	output2<<"zcoord angle azimuthal res_type res# total_energy fa_atr fa_rep fa_sol fa_water_to_bilayer f_elec_lipidlayer fa_imm_elec fa_elec"<<std::endl;

	// Write header, including score types available in the given energy function
	output << "loop zcoord angle azimuthal total_score ";
	ScoreTypes const nonzero_types( sfxn_->get_nonzero_weighted_scoretypes() );
	for ( core::Size ii = 1; ii <= nonzero_types.size(); ++ii ) {
		output << name_from_score_type( nonzero_types[ii] ) << " ";
	}
	output << std::endl;


	PackRotamersMoverOP pack_mover = get_pH_aware_packer( sfxn_ );

	std::string filename5( tempstr+"_posetomembrane_"+std::to_string( end_z_ ) + "_" + std::to_string( start_z_ ) +".pdb" );
	pose.dump_pdb( filename5 );
	//Side chains are packed once it is inside the membrane.

	// Calculate the protein length to determine bounds
	TR << "Computing helix length to determine landscape boundaries" << std::endl;
	//core::Real start_z( pose.residue( 1 ).xyz( 2 ).z() );
	//core::Real end_z( pose.residue( pose.total_residue() ).xyz( 2 ).z() );

	core::Real length( std::abs( pose.residue( 1 ).xyz( 2 ).z()  - pose.residue( pose.total_residue() ).xyz( 2 ).z() ) );
	core::Real limit1( 0 );
	core::Real limit2( 0 );

	limit1 = start_z_;
	limit2 = end_z_;

	//it is better to do directly depth, instead of making it a length dependent task.
	TR << "limit 2 is:~" << limit2 << std::endl;
	TR << "The length is:~" << length << std::endl;

	// get the membrane jump
	core::Size membrane_jump( pose.conformation().membrane_info()->membrane_jump() );

	// Apply an initial translation of the membrane center
	Vector initial_move( 0, 0, limit1 );
	TranslationMoverOP initialize( utility::pointer::make_shared< TranslationMover >( initial_move, membrane_jump ) );
	initialize->apply( pose );

	/*std::string filename7( tempstr+"_poseafterinitimove"+".pdb" );
	pose.dump_pdb( filename7 );*/

	Vector init_center( pose.conformation().membrane_info()->membrane_center( pose.conformation() ) );
	TR << "Shifting the protein to a new center position of (" << initial_move.x() << "," << initial_move.y() << "," << initial_move.z() << ")" << std::endl;

	peptide_center = getcenter( pose, flag_axis_ );

	//TR <<"center while transfering" << peptide_center.x() << "|" << peptide_center.y() << "|" << peptide_center.z() <<std::endl;


	core::Vector difference_in_center = initial_move - peptide_center;
	TranslationMoverOP adjusting_center( utility::pointer::make_shared< TranslationMover >( difference_in_center, membrane_jump ) );
	adjusting_center->apply( pose );

	peptide_center = getcenter( pose, flag_axis_ );


	// Set up deltas for rotations in normal and azimuthal angle & translations and spin axes
	// We densely sample the grid to more easily inspect for discontinuities

	core::Real rotation_delta( 1 ); // Degrees
	//core::Real azimuthal_delta( 30 ); //Degrees
	core::Vector translation_delta( 0, 0, 1 ); // Angstroms
	TR << "Setting a translation delta of: " << translation_delta.z() << std::endl;
	TR << "Setting a rotation delta of " << rotation_delta << std::endl;

	// Set up a constant translation mover along the z axis
	TranslationMoverOP translate_memb( utility::pointer::make_shared< TranslationMover >( translation_delta, membrane_jump ) );

	// Setup the axis of rotation (formerly the membrane normal
	core::Vector axis( 0,0,0 );
	if ( rotation_type_ == "YZ" || rotation_type_ == "XZ" ) {
		axis.z() = 1;
	} else if ( rotation_type_ == "XY" ) {
		axis.x() = 1;
	} else {
		utility_exit_with_message( "Unknown rotation type " + rotation_type_ );
	}


	//Align the protein with the z-axis as the initial position
	core::Vector peptide_normal( 0,0,0 );
	peptide_normal = getaxis( pose, flag_axis_ );
	TR << " Peptide normal is ::" << peptide_normal.x() << "i+ " << peptide_normal.y() << "j+" << peptide_normal.z() << "k" << std::endl;
	//rotate the peptide to the axis
	core::Real angle = numeric::conversions::degrees( angle_of( peptide_normal, axis ) );
	TR << "Angle with z-axis is: " << angle << "degrees" << std::endl;

	RotationMoverOP align_axis( utility::pointer::make_shared< RotationMover >( peptide_normal, axis, peptide_center, membrane_jump ) );
	align_axis->apply( pose );
	Vector rot_center = pose_tm_com( pose );
	peptide_normal = getaxis( pose, flag_axis_ );
	peptide_center = getcenter( pose, flag_axis_ );

	/*std::string filename_aa( tempstr + "_after_align_start" + std::to_string( end_z_ )  + "_end" + std::to_string( start_z_ ) + ".pdb" );
	pose.dump_pdb( filename_aa );*/

	TR << " peptide normal after alignment is ::" << peptide_normal.x() << "i+ " << peptide_normal.y() << "j +" << peptide_normal.z() << "k" << std::endl;
	TR << " peptide center after alignment is ::" << peptide_center.x() << "i+ " << peptide_center.y() << "j +" << peptide_center.z() << "k" << std::endl;


	// If testing interface sampling, rotate the pose by 90 degrees and reset the rotation axis
	if ( interface_ ) {

		TR << "Interface Mode: Rotate the pose by 90 degrees in the YZ plane so it lies horizontal to the membrane plane" << std::endl;

		// this step gets the peptide onto its side
		core::Vector rot_center = pose_tm_com( pose );
		core::Real delta_degrees( 90 );
		core::Real rotation_delta( numeric::conversions::radians( delta_degrees ) );
		Vector new_normal( 0, axis.y()*cos(rotation_delta) + axis.z()*sin(rotation_delta), -axis.y()*sin(rotation_delta) + axis.z()*cos(rotation_delta) );
		RotationMoverOP rotate_onto_side( utility::pointer::make_shared< RotationMover >( axis, new_normal, peptide_center, membrane_jump ) );
		rotate_onto_side->apply( pose );

		peptide_normal = getaxis( pose, flag_axis_ );
		peptide_center = getcenter( pose, flag_axis_ );
		TR << " peptide normal after rotation is ::" << peptide_normal.x() << "i+ " << peptide_normal.y() << "j +" << peptide_normal.z() << "k" << std::endl;
		TR << " peptide center after rotation is ::" << peptide_center.x() << "i+ " << peptide_center.y() << "j +" << peptide_center.z() << "k" << std::endl;

	}


	// Translate along the z axis in very small steps
	for ( core::Real z_coord = limit1; z_coord < limit2 ; z_coord += translation_delta.z() ) {

		// Define the rotation center from the transmembrane center of mass
		Vector rot_center = pose_tm_com( pose );

		peptide_center = getcenter( pose, flag_axis_ );
		//TR <<"center while in the loop" << peptide_center.x() << "|" << peptide_center.y() << "|" << peptide_center.z() <<std::endl;

		// Rotate from zero to 360 degrees in very small steps
		for ( core::Real normal_angle = 0; normal_angle <= 180; normal_angle += rotation_delta ) {

			// Verify the membrane information
			TR << "Membrane center: " << pose.conformation().membrane_info()->membrane_center( pose.conformation() ) << std::endl;


			for ( core::Real azimuthal_angle = 0; azimuthal_angle <= 360; azimuthal_angle += azimuthal_delta_ ) {

				// Make a copy of the pose to avoid rounding errors
				core::pose::PoseOP pose_copy = pose.clone();

				TR << " azimuthal angle " << azimuthal_angle << std::endl;
				// Verify the membrane information
				TR << "Membrane center:" << pose.conformation().membrane_info()->membrane_center( pose.conformation() ) << std::endl;

				using namespace core::conformation::membrane;


				/*Azimuthal angle rotation about the axis of minimum moment of Inertia
				The axis is calculated in the RotationMover*/
				/*if flag==2 :: eigen vector
				flag == 1 :: vector joining the coordinates of CA*/

				RotationMoverOP rotate_azim( utility::pointer::make_shared< RotationMover >( peptide_normal, peptide_normal, peptide_center, membrane_jump, azimuthal_angle ) );
				rotate_azim->apply( *pose_copy );

				/*Rotation about the Normal angle*/

				RotationMoverOP rotate( get_rotation( normal_angle, peptide_normal, peptide_center, membrane_jump, rotation_type_ ));
				rotate->apply( *pose_copy );

				if ( (z_coord == 0.0 || z_coord == 40.0) && (std::fmod(normal_angle,30.0) == 0.0) ) {

					std::string filename4( tempstr + "limit_" + std::to_string(z_coord) + "_norm" + std::to_string( normal_angle ) + "_azrot" + std::to_string( azimuthal_angle ) +".pdb" );
					pose_copy->dump_pdb( filename4 );
				}


				if ( pH_mode_ ) {
					TR << "=========================================" << std::endl;
					TR << "the pH mode implementation is under constaruction. It will be added in future versions. In this form the calculations are done at neutral pH." << std::endl;
					TR << "=========================================" << std::endl;

				} else {
					core::Real nloop( 0.0 );

					pack_mover->apply( *pose_copy );
					// Write data to output file
					output << nloop << " " << z_coord << " " << normal_angle << " " << azimuthal_angle;
					output << " " << sfxn_->score( *pose_copy );
					for ( core::Size ii = 1; ii <= nonzero_types.size(); ++ii ) {
						output << " " << sfxn_->score_by_scoretype( *pose_copy, nonzero_types[ii] );
					} // finish writing scores
					output <<std::endl;

					if ( z_coord == limit1 ) {
						if ( (normal_angle == 0 || normal_angle == 90) && (azimuthal_angle == 0 || azimuthal_angle == 90) ) {
							utility::vector1< bool > is_scoringres( pose.size(), false);
							//is_scoringres[ resno ] = true;

							// Energies E( pose_copy->energies() );
							for ( core::Size jj=1; jj<pose.total_residue(); jj++ ) {
								is_scoringres[ jj ] = true;
								core::Real score = pose_copy->energies().residue_total_energies(jj)[total_score];
								core::Real res_score_1 = pose_copy->energies().residue_total_energies(jj)[fa_water_to_bilayer];
								core::Real res_score_2 = pose_copy->energies().residue_total_energies(jj)[f_elec_lipidlayer];
								core::Real res_score_3 = pose_copy->energies().residue_total_energies(jj)[fa_imm_elec];
								core::Real res_score_5 = pose_copy->energies().residue_total_energies(jj)[fa_elec];
								core::Real res_score_6 = pose_copy->energies().residue_total_energies(jj)[fa_atr];
								core::Real res_score_7 = pose_copy->energies().residue_total_energies(jj)[fa_rep];
								core::Real res_score_8 = pose_copy->energies().residue_total_energies(jj)[fa_sol];


								output2 << z_coord << " " << normal_angle << " " << azimuthal_angle << " ";
								output2 << pose_copy->residue(jj).name() << " ";
								output2 << jj << " " << score << " " << res_score_6 << " " << res_score_7 << " " << res_score_8 << " " << res_score_1 << " " << res_score_2 << " " << res_score_3 << " " << res_score_5 << std::endl;
								is_scoringres[ jj ] = false;
							}
							std::string filename4( tempstr + "limit_" + std::to_string(z_coord) + "_norm" + std::to_string( normal_angle ) + "_azrot" + std::to_string( azimuthal_angle )  + + "_nloop" + std::to_string( nloop ) + "afterdesign.pdb" );
							pose_copy->dump_pdb( filename4 );
							//    E.show(jj);

						}
					}

				}
			}// end of azimuthal angle loop


		}// end of normal angle loop

		// increment the translation
		translate_memb->apply( pose );


	} // end translation loop

} // apply

/// @brief initialize the options
void
MembraneEnergyLandscapeSampler::init_from_options() {

	using namespace basic::options;
	if ( option[ OptionKeys::pH::pH_mode ].user() ) {
		pH_mode_ = option[ OptionKeys::pH::pH_mode ].value();
	}
}

/// @brief get user input pH value
core::Real
MembraneEnergyLandscapeSampler::get_pH_value() {

	using namespace basic::options;
	if ( pH_mode_ ) {
		return option[ OptionKeys::pH::value_pH ].value();
	} else {
		utility_exit_with_message( " pH mode is not enabled. calculations being done at neutral pH. " );
	}
}

/// @brief count the number of residue "res" in the pose "pose"
core::Real
MembraneEnergyLandscapeSampler::count_res(std::string res, core::pose::Pose const & pose){
	core::Real count(0.0);
	for ( core::Size ii=1; ii<pose.total_residue(); ++ii ) {
		if ( pose.residue(ii).name() == res ) {
			count = count+1.0;
		}
	}

	return count;
}
/// @brief count the number of occurences the 3 letter name of residue is different in pose1 and pose2
core::Real
MembraneEnergyLandscapeSampler::count_diff(core::pose::Pose const & pose1, core::pose::Pose const & pose2){
	core::Real count(0.0);
	for ( core::Size ii=1; ii<pose1.total_residue(); ++ii ) {
		//TR << "residues: " << pose2.residue(ii).name3() << std::endl;
		if ( pose1.residue(ii).name3() != pose2.residue(ii).name3() ) {
			count = count+1.0;
		}
	}

	return count;
}



core::Vector
MembraneEnergyLandscapeSampler::get_rotation_axis() {

	// Setup the axis of rotation (formerly the membrane normal
	core::Vector axis(0,0,0);
	if ( rotation_type_ == "YZ" || rotation_type_ == "XZ" ) {
		axis.z() = 1;
	} else if ( rotation_type_ == "XY" ) {
		axis.x() = 1;
	} else {
		utility_exit_with_message( "Unknown rotation type " + rotation_type_ );
	}

	return axis;
}


protocols::minimization_packing::PackRotamersMoverOP
MembraneEnergyLandscapeSampler::get_pH_aware_packer( core::scoring::ScoreFunctionOP sfxn ) const {

	//using namespace protocols::task_operations;
	using namespace core::pack::task::operation;
	using namespace protocols::minimization_packing;
	using namespace core::pack::task;
	using namespace core::pack::palette;

	//NoDesignPackerPaletteOP pp( new NoDesignPackerPalette() );
	//this was totally stopping from design
	DefaultPackerPaletteOP pp( utility::pointer::make_shared< DefaultPackerPalette >() );
	TaskFactoryOP tf( new TaskFactory() );
	tf->set_packer_palette(pp);

	PackRotamersMoverOP packer( utility::pointer::make_shared< PackRotamersMover >() );

	if ( pH_mode_ ) {

		// tf->push_back( utility::pointer::make_shared< pHVariantTaskOperation >() );
		tf->push_back( utility::pointer::make_shared< core::pack::task::operation::RestrictToRepacking>() );
		packer->score_function(sfxn);
		packer->task_factory(tf);
		packer->nloop(3);

	} else {

		tf->push_back( utility::pointer::make_shared< core::pack::task::operation::RestrictToRepacking>() );
		PackRotamersMoverOP packer( utility::pointer::make_shared< PackRotamersMover >() );
		packer->score_function(sfxn);
		packer->task_factory(tf);

	}
	//This packertask after every move is different from franklin2019;

	return packer;


}

protocols::membrane::RotationMoverOP
MembraneEnergyLandscapeSampler::get_rotation( core::Real normal_angle, core::Vector axis, core::Vector rot_center, core::Size membrane_jump, std::string rotation_type ) {

	using namespace protocols::membrane;
	using namespace numeric::conversions;

	// Pick the rotation delta from the normal angle
	core::Real delta( radians( normal_angle ) );


	// Calculate a new normal based on an XY or XZ rotation
	core::Vector updated_axis( 0, 0, 0 );

	if ( rotation_type == "XZ" ) {
		updated_axis.x() = axis.x() * cos(delta) + axis.z() * sin(delta);
		updated_axis.y() = 0;
		updated_axis.z() = -axis.x() * sin(delta) + axis.z() * cos(delta);
	} else if ( rotation_type == "YZ" ) {
		updated_axis.x() = 0;
		updated_axis.y() = axis.y() * cos(delta) + axis.z() * sin(delta);
		updated_axis.z() = -axis.y() * sin(delta) + axis.z() * cos(delta);
	} else if ( rotation_type == "XY" ) {
		updated_axis.x() = axis.x() * cos(delta) + axis.y() * sin(delta);
		updated_axis.y() = -axis.x() * sin(delta) + axis.y() * cos(delta);
		updated_axis.z() = 0;
	}

	RotationMoverOP rotation( utility::pointer::make_shared< RotationMover >( axis, updated_axis, rot_center, membrane_jump ) );
	return rotation;

}

core::Vector
MembraneEnergyLandscapeSampler::getaxis( core::pose::Pose & pose, core::Real flag ){

	using namespace protocols::membrane;
	core::Real lastres( pose.total_residue());
	core::Vector new_normal( 0,0,0 );
	core::Vector rot_center = pose_tm_com( pose );
	core::Vector rot_center2( 0, 0, 0 );
	/*if flag==2 :: eigen vector
	flag == 0 :: vector joining the coordinates of CA
	flag == 1 :: average of all TM span vectors. A TM span
	vector is one joining the start and end CA atoms (from
	bottom to top) of a helix or beta sheet spanning TM */

	if ( flag == 2 ) {

		numeric::xyzMatrix<core::Real> Mat1( numeric::xyzMatrix<core::Real>::identity() );
		core::Vector xyz_w_w;
		numeric::xyzMatrix<core::Real> xyz_eVec( numeric::xyzMatrix<core::Real>::identity() );

		core::Real const tol( 0.000001 );
		core::Real minimumelement( 0 );
		core::Real counter_min( 0 );

		Mat1.xx() = 0; Mat1.xy() = 0; Mat1.xz() = 0;
		Mat1.yx() = 0; Mat1.yy() = 0; Mat1.yz() = 0;
		Mat1.zx() = 0; Mat1.zy() = 0; Mat1.zz() = 0;

		//core::Real count_tm( 0 );
		/*lastres is the membrane*/
		for ( core::Size r=1; r<lastres; r++ ) {

			rot_center2 = rot_center2 + pose.residue(r).xyz(2);
			//  TR << r << "th ele" << pose.residue(r).xyz(2).x()<<"|" << pose.residue(r).xyz(2).y()<<"|" <<pose.residue(r).xyz(2).z()<< std::endl;
		}
		rot_center2 = rot_center2/( lastres-1 );

		TR << "center of CA : " << rot_center2.x() << "i+ " << rot_center2.y() << "j+ " << rot_center2.z() << "k" << std::endl;
		core::Vector avg_axis ( 0, 0, 0 );


		for ( core::Size r=1; r<lastres; r++ ) {

			core::Vector pose_centered ( 0, 0, 0 );


			//only the alpha_atoms
			pose_centered = pose.residue(r).xyz(2) - rot_center2;

			//if( abs( pose.residue(r).xyz(t).z() ) < 15 ){
			Mat1.xx() = Mat1.xx() + ( pose_centered.y()*pose_centered.y() + pose_centered.z()*pose_centered.z() );
			Mat1.yy() = Mat1.yy() + ( pose_centered.x()*pose_centered.x() + pose_centered.z()*pose_centered.z() );
			Mat1.zz() = Mat1.zz() + ( pose_centered.x()*pose_centered.x() + pose_centered.y()*pose_centered.y() );
			Mat1.xy() = Mat1.xy() - ( pose_centered.x()*pose_centered.y() );
			Mat1.xz() = Mat1.xz() - ( pose_centered.x()*pose_centered.z() );
			Mat1.yz() = Mat1.yz() - ( pose_centered.y()*pose_centered.z() );

			avg_axis.x() = avg_axis.x() + pose_centered.x();
			avg_axis.y() = avg_axis.y() + pose_centered.y();
			avg_axis.z() = avg_axis.z() + pose_centered.z();

			//   }
			//  }

		}
		//avg_axis.normalize();

		Mat1.yx() = Mat1.xy();
		Mat1.zx() = Mat1.xz();
		Mat1.zy() = Mat1.yz();

		/*Finding the minimum Eigen value*/
		xyz_w_w = numeric::eigenvector_jacobi( Mat1, tol, xyz_eVec);
		minimumelement = xyz_w_w[0];
		for ( core::Size r=0; r<=2; r++ ) {

			if ( minimumelement > xyz_w_w[r] ) {
				minimumelement = xyz_w_w[r];
				counter_min = r;
			}

		}

		if ( counter_min == 0 ) {
			TR << "Eigen vector" << xyz_eVec.xx() << "|" << xyz_eVec.yx() << "|" << xyz_eVec.zx() <<std::endl;
			new_normal.x() = xyz_eVec.xx();
			new_normal.y() = xyz_eVec.yx();
			new_normal.z() = xyz_eVec.zx();
		} else if ( counter_min == 1 ) {
			TR << "Eigen vector" << xyz_eVec.xy() << "|" << xyz_eVec.yy() << "|" << xyz_eVec.zy() <<std::endl;
			new_normal.x() = xyz_eVec.xy();
			new_normal.y() = xyz_eVec.yy();
			new_normal.z() = xyz_eVec.zy();
		} else {
			TR << "Eigen vector" << xyz_eVec.xz() << "|" << xyz_eVec.yz() << "|" << xyz_eVec.zz() <<std::endl;
			new_normal.x() = xyz_eVec.xz();
			new_normal.y() = xyz_eVec.yz();
			new_normal.z() = xyz_eVec.zz();
		}



	} else if ( flag == 0 ) {


		core::Real modulus( 0 );
		core::Vector center_top( 0, 0, 0 );
		core::Vector center_bottom( 0, 0, 0 );

		/*idea2:: instead we take the line joining diagonally oposite points at the top layer and bottom layer*/
		/*few proteins had the last few alpha atoms not on the helix*/
		center_top = 0.5*( pose.residue(lastres-4).xyz(2) + pose.residue(lastres-6).xyz(2) );
		center_bottom = 0.5*( pose.residue(7).xyz(2) + pose.residue(5).xyz(2) );

		new_normal = center_top - center_bottom;
		//rot_center2 = 0.5* ( center_top + center_bottom );

		modulus = sqrt( dot( new_normal, new_normal ) );
		new_normal = new_normal/modulus;

		//TR << " new_normal is ::" << new_normal.x() << "i+ " << new_normal.y() << "j +" << new_normal.z() << "k" << std::endl;
		//TR << "new_center is ::" << rot_center2.x() << "i+" << rot_center2.y() << "j +" << rot_center2.z() << "k" << std::endl;
	} else if ( flag == 1 ) {
		/*trying the average of tm helices*/

		using namespace core::conformation::membrane;

		// get topology from MembraneInfo
		SpanningTopologyOP topo( pose.membrane_info()->spanning_topology() );

		// initialize vector
		core::Vector tm_vector( 0, 0, 0 );
		core::Vector avg_vector( 0, 0, 0 );

		// go through topology and avg coords of start and end residues
		for ( core::Size i = 1; i <= topo->nspans(); ++i ) {

			// get CA coords

			core::Vector start_coord = pose.residue( topo->span(i)->start() ).xyz("CA");
			core::Vector end_coord = pose.residue( topo->span(i)->end() ).xyz("CA");
			//TR << "start res:" << topo->span(i)->start() <<std::endl;
			//TR << "end res:" << topo->span(i)->end() <<std::endl;

			if ( ( start_coord.z() - end_coord.z() ) > 0.0 ) {

				tm_vector = start_coord - end_coord ;

			} else {

				tm_vector = end_coord - start_coord ;

			}

			avg_vector = avg_vector + tm_vector ;


		}
		avg_vector = avg_vector/( topo->nspans() );
		avg_vector.normalize();
		new_normal = avg_vector;
		TR << " new_normal is ::" << new_normal.x() << "i+ " << new_normal.y() << "j +" << new_normal.z() << "k" << std::endl;
		//TR << "new_center is ::" << rot_center2.x() << "i+" << rot_center2.y() << "j +" << rot_center2.z() << "k" << std::endl;

	} else {
		TR <<" wrong choice of flag" << std::endl;

	}

	return( new_normal );
}

core::Vector
MembraneEnergyLandscapeSampler::getcenter( core::pose::Pose & pose, core::Real flag ){

	core::Real lastres( pose.total_residue() );
	core::Vector rot_center2( 0, 0, 0 );
	/*if flag==2 :: eigen vector
	flag == 0 :: vector joining the coordinates of CA
	flag == 1 :: membrane normal-axis; used as the default option*/

	if ( flag == 2 ) {

		for ( core::Size r=1; r<lastres; r++ ) {

			rot_center2 = rot_center2 + pose.residue(r).xyz(2);

		}
		rot_center2 = rot_center2/( lastres-1 );

	} else if ( flag == 0 ) {

		/*if we take the vector joining the first and last CA*/
		core::Vector center_top( 0, 0, 0 );
		core::Vector center_bottom( 0, 0, 0 );
		center_top = 0.5*( pose.residue(lastres-4).xyz(2) + pose.residue(lastres-6).xyz(2) );
		center_bottom = 0.5*( pose.residue(7).xyz(2) + pose.residue(5).xyz(2) );

		rot_center2 = 0.5* ( center_top + center_bottom );
	} else if ( flag == 1 ) {

		rot_center2 = 0.0;
		for ( core::Size r=1; r<lastres; r++ ) {

			rot_center2 = rot_center2 + pose.residue(r).xyz(2);

		}
		rot_center2 = rot_center2/( lastres-1 );

	}

	return( rot_center2 );


}


////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
MembraneEnergyLandscapeSampler::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
// declaring an overload
void
MembraneEnergyLandscapeSampler::show() const
{
	show( TR );
}
////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
MembraneEnergyLandscapeSampler::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
) {

	// Get scoring function and set the weights
	using namespace protocols::rosetta_scripts;
	using namespace core::scoring;
	ScoreFunctionOP candidate_sfxn( parse_score_function( tag, "scorefxn", data ) );
	if ( candidate_sfxn == nullptr ) {
		utility_exit_with_message("No energy function specified!" );
	} else {
		sfxn_ = candidate_sfxn;
		sfxn_weights_ = sfxn_->get_name();
		TR << "Setting up a new scoring function with weights " << sfxn_->get_name() << std::endl;
	}

	// Specify rotation type
	if ( tag->hasOption( "rotation_type" ) ) {
		rotation_type_ = tag->getOption< std::string >( "rotation_type" );
	}

	// Should I treat this like an interface landscape problem?
	if ( tag->hasOption( "interface" ) ) {
		interface_ = tag->getOption< bool >( "interface" );
	}

	//What is the z_lowerlimit
	if ( tag->hasOption( "start_z" ) ) {
		start_z_ = tag->getOption< core::Real >( "start_z");

	}

	//What is the z_upperlimit
	if ( tag->hasOption( "end_z" ) ) {
		end_z_ = tag->getOption< core::Real >( "end_z");

	}

	//For azimuthal_angle rotation, is the axis center-to-center
	//or axis of minimum moment of inertia
	if ( tag->hasOption( "flag_axis" ) ) {
		flag_axis_ = tag->getOption< core::Real >( "flag_axis" );
	}

	//The frequency of azimuthal or rotation angle. For peptides, we use the default option ~30 or 45;
	//for multipass proteins; we must use this as 5.
	if ( tag->hasOption( "azimuthal_delta" ) ) {
		azimuthal_delta_ = tag->getOption< core::Real >( "azimuthal_delta" );
	}

	// Should I repack each pose prior to scoring
	if ( tag->hasOption( "repack" ) ) {
		repack_ = tag->getOption< bool >( "repack" );
	}

	// Should I incorporate protonation variants during repacking
	if ( tag->hasOption( "pH_mode" ) ) {
		pH_mode_ = tag->getOption< bool >( "pH_mode" );
	}

}

////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
MembraneEnergyLandscapeSampler::fresh_instance() const
{
	return utility::pointer::make_shared< MembraneEnergyLandscapeSampler >();
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
MembraneEnergyLandscapeSampler::clone() const
{
	return utility::pointer::make_shared< MembraneEnergyLandscapeSampler >(*this);
}

std::string MembraneEnergyLandscapeSampler::get_name() const {
	return mover_name();
}

std::string MembraneEnergyLandscapeSampler::mover_name() {
	return "MembraneEnergyLandscapeSampler";
}

void MembraneEnergyLandscapeSampler::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;
	rosetta_scripts::attributes_for_parse_score_function( attlist, "scorefxn" );
	attlist
		+ XMLSchemaAttribute( "sfxn_weights", xs_string, "Energy function weights file" )
		+ XMLSchemaAttribute( "rotation_type", xs_string, "Rotation axis: XY, XZ, or YZ" )
		+ XMLSchemaAttribute( "interface", xsct_rosetta_bool, "Should I treat this like an interface landscape scanning problem?")
		+ XMLSchemaAttribute( "start_z", xsct_real, "lower limit of z-coord" )
		+ XMLSchemaAttribute( "end_z", xsct_real, "upper limit of z-coord" )
		+ XMLSchemaAttribute( "flag_axis", xsct_real, "select the axis of rotation" )
		+ XMLSchemaAttribute( "azimuthal_delta", xsct_real, "select the frequency of azimuthal/rotation angle" )
		+ XMLSchemaAttribute( "repack", xsct_rosetta_bool, "Should I repack each pose prior to scoring?" )
		+ XMLSchemaAttribute( "pH_mode", xsct_rosetta_bool, "Should I include protonation variants during packing?" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Sample the membrane energy function landscape using a given energy function", attlist );

}

/////////////// Creator ///////////////

protocols::moves::MoverOP
MembraneEnergyLandscapeSamplerCreator::create_mover() const
{
	return utility::pointer::make_shared< MembraneEnergyLandscapeSampler >();
}

std::string
MembraneEnergyLandscapeSamplerCreator::keyname() const
{
	return MembraneEnergyLandscapeSampler::mover_name();
}

void MembraneEnergyLandscapeSamplerCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MembraneEnergyLandscapeSampler::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////

std::ostream &
operator<<( std::ostream & os, MembraneEnergyLandscapeSampler const & mover )
{
	mover.show(os);
	return os;
}

} //protocols
} //membrane_benchmark
