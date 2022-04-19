// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW CoMotion, email: license@uw.edu.

//This app should optimize the parameters of either a bicelle, vesicle, or double vesicle membrane geometry.
//
//
#include <iostream>
#include <cmath>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <devel/init.hh>
#include <utility/vector1.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <utility/pointer/owning_ptr.hh>
#include <core/pose/variant_util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Conformation.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>


#include <numeric/xyzVector.hh>
#include <core/id/AtomID.hh>

//membrane specific
#include <protocols/membrane/AddMembraneMover.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/MembraneGeometry.hh>
#include <core/conformation/membrane/membrane_geometry/Bicelle.hh>
#include <core/conformation/membrane/membrane_geometry/DoubleVesicle.hh>
#include <core/conformation/membrane/membrane_geometry/Vesicle.hh>
#include <protocols/membrane/OptimizeProteinEmbeddingMover.hh>
#include <protocols/membrane/OptimizeMembranePositionMover.hh>

#include <utility/io/ozstream.hh>


//Updates the MembraneGeometry in pose to be a Vesicle with the given radius
//called in opt_vesicle_params to update the raidus
void update_vesicle_radius( core::pose::Pose & pose, core::conformation::membrane::membrane_geometry::VesicleCOP vesicle, core::Real radius ) {
	//create new Vesicle with the specified radius
	core::conformation::membrane::membrane_geometry::VesicleCOP mp_ves_new( new core::conformation::membrane::membrane_geometry::Vesicle( vesicle->membrane_steepness(), vesicle->membrane_thickness(), radius ));
	//convert VesicleCOP to MembraneGeometryCOP, so it can be used to replace the MembraneGeometryCOP in MembraneInfo
	core::conformation::membrane::MembraneGeometryCOP mp_geo_new = utility::pointer::dynamic_pointer_cast< core::conformation::membrane::MembraneGeometry const > (mp_ves_new );
	//replace old membrane geometry, with the new one with updated radius
	pose.conformation().membrane_info()->set_membrane_geometry( mp_geo_new );
}

//Updates the MembraneGeometry in pose to be a Bicelle with the given inner_radius
////called in opt_bicelle_params to update the bicelle inner_radius
void update_bicelle_radius( core::pose::Pose & pose, core::conformation::membrane::membrane_geometry::BicelleCOP bicelle, core::Real inner_radius ) {
	//create new Bicelle with the specified radius
	core::conformation::membrane::membrane_geometry::BicelleCOP mp_bic_new( new core::conformation::membrane::membrane_geometry::Bicelle( bicelle->membrane_steepness(), bicelle->membrane_thickness(), inner_radius ));
	//convert BicelleCOP to MembraneGeometryCOP, so it can be used to replace the MembraneGeometryCOP in MembraneInfo
	core::conformation::membrane::MembraneGeometryCOP mp_geo_new = utility::pointer::dynamic_pointer_cast< core::conformation::membrane::MembraneGeometry const > (mp_bic_new );
	//replace old membrane geometry, with the new one with updated radius
	pose.conformation().membrane_info()->set_membrane_geometry( mp_geo_new );
}

//Updates the MembraneGeometry in pose to be a DoubleVesicle with the given distance between the two vesicle membranes and the outer_radius
//called in opt_double_vesicle_params to update either the distance or outer_radius parameters
void update_double_vesicle( core::pose::Pose & pose, core::conformation::membrane::membrane_geometry::DoubleVesicleCOP dvesicle, core::Real distance, core::Real outer_r ) {
	//create new DoubelVesicle with the specified radius
	core::conformation::membrane::membrane_geometry::DoubleVesicleCOP mp_dves_new( new core::conformation::membrane::membrane_geometry::DoubleVesicle( dvesicle->membrane_steepness(), dvesicle->membrane_thickness(), outer_r, distance ));
	//convert DoubleVesicleCOP to MembraneGeometryCOP, so it can be used to replace the MembraneGeometryCOP in MembraneInfo
	core::conformation::membrane::MembraneGeometryCOP mp_geo_new = utility::pointer::dynamic_pointer_cast< core::conformation::membrane::MembraneGeometry const > (mp_dves_new );
	//replace old membrane geometry, with the new one with updated parameters
	pose.conformation().membrane_info()->set_membrane_geometry( mp_geo_new );
}

//called in main to find vesicle radius that gives lowest energy for the current pose
void opt_vesicle_params( core::pose::Pose & pose, core::conformation::membrane::MembraneGeometryCOP mp_geometry, core::scoring::ScoreFunctionOP sfxn ) {
	core::conformation::membrane::membrane_geometry::VesicleCOP mp_vesicle = utility::pointer::dynamic_pointer_cast< core::conformation::membrane::membrane_geometry::Vesicle const > ( mp_geometry );

	core::Real lowest_score( sfxn->score( pose ) );
	core::Real best_r(mp_vesicle->get_radius());
	core::Real best_r_coarse( best_r );

	core::Real big_step_size( 10 );
	core::Real small_step_size( 1 );

	//starting at radius = 20, increase by big_step_size until 600 A, keep lowest scoring
	for ( core::Real r=20; r<=600; r+=big_step_size ) {
		update_vesicle_radius( pose, mp_vesicle, r );
		core::Real temp_score = sfxn->score( pose);
		// save best radius
		if ( temp_score < lowest_score ) {
			best_r_coarse = r;
			lowest_score = temp_score;
		}
	}

	//set radius for the lowest scoring value
	//starting at current radius minus big_step_size + small_step_size, loop through radius with small step size until at current_radius + big_step_size - small_step_size
	core::Real r_min = best_r_coarse-big_step_size+small_step_size;
	for ( core::Real r=r_min; r < best_r_coarse+big_step_size; r+= small_step_size ) {
		update_vesicle_radius( pose, mp_vesicle, r );
		//opt->apply( pose );
		core::Real temp_score = sfxn->score( pose);
		// save best radius
		if ( temp_score <= lowest_score ) {
			best_r = r;
			lowest_score = temp_score;
		}
	}


	update_vesicle_radius( pose, mp_vesicle, best_r );
	core::Real new_score( sfxn->score( pose ) );
	std::cout << "new radius: " << best_r << std::endl;
	std::cout << "new score: " << new_score << std::endl;

	utility::vector1< std::string > temp( utility::string_split( pose.pdb_info()->name(), '/') );
	std::string tempstr = temp[ temp.size() ].substr(0, temp[ temp.size() ].size()-4 );
	std::string filename( tempstr + "_vesicle_opt_params.txt" );
	utility::io::ozstream output( filename );
	output << tempstr << std::endl;
	output << "Geometry: Vesicle" << std::endl;
	output << "Optimized Radius: " << best_r << std::endl;
}

//called in main to find the radius of the outer vesicle and distance to the inner vesicle membrane that gives lowest energy for the current pose
void opt_double_vesicle_params( core::pose::Pose & pose, core::conformation::membrane::MembraneGeometryCOP mp_geometry, core::scoring::ScoreFunctionOP sfxn ) {
	core::conformation::membrane::membrane_geometry::DoubleVesicleCOP mp_dvesicle = utility::pointer::dynamic_pointer_cast< core::conformation::membrane::membrane_geometry::DoubleVesicle const > ( mp_geometry );

	core::Real lowest_score( sfxn->score( pose ) );
	core::Real best_dis( mp_dvesicle->get_distance() );
	core::Real best_outer_r( mp_dvesicle->get_outer_radius() );
	core::Real best_outer_r_coarse( best_outer_r );


	//optimize distance first
	//distance between outer edge of inner membrane and inner edge of outer membrane
	for ( core::Real d = 10; d<100; d++ ) {
		update_double_vesicle( pose, mp_dvesicle, d, best_outer_r); //best outer_r here is just the starting outer radius because it hasn't been updated
		core::Real temp_score = sfxn->score( pose );
		//save best distance
		if ( temp_score < lowest_score ) {
			best_dis = d;
			lowest_score = temp_score;
		}
	}

	update_double_vesicle( pose, mp_dvesicle, best_dis, best_outer_r );
	core::Real min_outer_r = 20 + best_dis + 2*mp_dvesicle->membrane_thickness(); //if min inner_radius = 20
	core::Real max_outer_r = 600 + best_dis + 2*mp_dvesicle->membrane_thickness(); //if max inner_radius = 600

	core::Real big_step_size( 10 );
	core::Real small_step_size( 1 ) ;

	//optimize outer_radius coarse (inner_radius depends on outer radius and distance, so it will be updated)
	for ( core::Real r = min_outer_r; r< max_outer_r; r+=big_step_size ) {
		update_double_vesicle( pose, mp_dvesicle, best_dis, r );
		core::Real temp_score = sfxn->score( pose );
		//save best outer radius
		if ( temp_score <= lowest_score ) {
			best_outer_r_coarse = r;
			lowest_score = temp_score;
		}
	}

	//optimize outer_radius
	for ( core::Real r=best_outer_r_coarse-big_step_size+small_step_size; r < best_outer_r_coarse; r+= small_step_size ) {
		update_double_vesicle( pose, mp_dvesicle, best_dis, r );
		core::Real temp_score = sfxn->score( pose );
		//save best outer radius
		if ( temp_score < lowest_score ) {
			best_outer_r = r;
			lowest_score = temp_score;
		}
	}

	update_double_vesicle( pose, mp_dvesicle, best_dis, best_outer_r );
	core::Real new_score( sfxn->score( pose ) );
	std::cout << "new distance: " << best_dis << std::endl;
	std::cout << "new outer radius: " << best_outer_r << std::endl;
	std::cout << "new score: " << new_score << std::endl;

	utility::vector1< std::string > temp( utility::string_split( pose.pdb_info()->name(), '/') );
	std::string tempstr = temp[ temp.size() ].substr(0, temp[ temp.size() ].size()-4 );
	std::string filename( tempstr + "_double_vesicle_opt_params.txt" );
	utility::io::ozstream output( filename );
	output << tempstr << std::endl;
	output << "Geometry: Double Vesicle" << std::endl;
	output << "Optimized Outer Radius: " << best_outer_r << std::endl;
	output << "Optimized distance: " << best_dis << std::endl;

}

//called in main to find bicelle inner_radius that gives lowest energy for the current pose
void opt_bicelle_params( core::pose::Pose & pose, core::conformation::membrane::MembraneGeometryCOP mp_geometry, core::scoring::ScoreFunctionOP sfxn ) {
	core::conformation::membrane::membrane_geometry::BicelleCOP mp_bicelle = utility::pointer::dynamic_pointer_cast< core::conformation::membrane::membrane_geometry::Bicelle const > ( mp_geometry );

	core::Real lowest_score( sfxn->score( pose ) );
	core::Real best_r( mp_bicelle->bicelle_inner_radius() );

	core::Real step_size( 0.5 );
	core::Real protein_slice_diameter( mp_bicelle->protein_slice_diameter() );

	//starting at inner_radius = 0, increase by step_size until twice the protein_slice_diameter, keep lowest scoring
	for ( core::Real r=0; r<=2*protein_slice_diameter; r+=step_size ) {
		update_bicelle_radius( pose, mp_bicelle, r );
		core::Real temp_score = sfxn->score( pose);
		// save best radius
		if ( temp_score < lowest_score ) {
			best_r = r;
			lowest_score = temp_score;
		}
	}


	update_bicelle_radius( pose, mp_bicelle, best_r );
	core::Real new_score( sfxn->score( pose ) );
	std::cout << "new radius: " << best_r << std::endl;
	std::cout << "new score: " << new_score << std::endl;

	utility::vector1< std::string > temp( utility::string_split( pose.pdb_info()->name(), '/') );
	std::string tempstr = temp[ temp.size() ].substr(0, temp[ temp.size() ].size()-4 );
	std::string filename( tempstr + "_bicelle_opt_params.txt" );
	utility::io::ozstream output( filename );
	output << tempstr << std::endl;
	output << "Geometry: Bicelle" << std::endl;
	output << "Optimized Inner Radius: " << best_r << std::endl;
}


int main( int argc, char ** argv ) {
	try {
		devel::init( argc, argv );

		//read in pdb provided in command line
		utility::vector1< std::string > filenames = basic::options::option[basic::options::OptionKeys::in::file::s ].value();
		if ( filenames.size() > 0 ) {
			std::cout << "You entered: " << filenames[ 1 ] << " as the PDB file to be read" << std::endl;
		} else {
			std::cout << "You didn’t provide a PDB file with the -in::file::s option" << std::endl;
			return 1;
		}

		//construct a pose object from pdb file
		core::pose::PoseOP pose_op = core::import_pose::pose_from_file( filenames[1] );
		core::pose::Pose& mypose = *pose_op;

		//add MEM virtual residue to pose
		protocols::membrane::AddMembraneMoverOP addmem( new protocols::membrane::AddMembraneMover() );
		addmem->apply( mypose );

		protocols::membrane::OptimizeProteinEmbeddingMoverOP opt( new protocols::membrane::OptimizeProteinEmbeddingMover() );
		//protocols::membrane::OptimizeMembranePositionMoverOP opt( new protocols::membrane::OptimizeMembranePositionMover() );
		opt->apply( mypose );

		//create conformation object from pose
		const core::conformation::Conformation& myconf = mypose.conformation();
		core::conformation::membrane::MembraneGeometryCOP mp_geometry( myconf.membrane_info()->membrane_geometry() );

		//get score funciton
		core::scoring::ScoreFunctionOP sfxn = core::scoring::ScoreFunctionFactory::create_score_function("mpframework_smooth_fa_2012.wts");

		//Identify membrane geometry
		//Membrane geometry is set in AddMembraneMover using command line option flags -mp:geometry bicelle, vesicle, or double_vesicle
		core::conformation::membrane::MP_GEOMETRY_TRANSITION mp_geo_enum( mypose.conformation().membrane_info()->membrane_geometry()->geometry_enum() );


		//different parameters depending on geometry
		switch ( mp_geo_enum )
				{
				case core::conformation::membrane::MP_GEOMETRY_TRANSITION::BICELLE :
					std::cout << "bicelle" << std::endl;
					opt_bicelle_params( mypose, mp_geometry, sfxn );
					break;
				case core::conformation::membrane::MP_GEOMETRY_TRANSITION::VESICLE :
					std::cout << "vesicle" << std::endl;
					opt_vesicle_params( mypose, mp_geometry, sfxn );
					break;
				case core::conformation::membrane::MP_GEOMETRY_TRANSITION::DOUBLE_VESICLE :
					std::cout << "double_vesicle" << std::endl;
					opt_double_vesicle_params( mypose, mp_geometry, sfxn);
					break;
				default :
					std::cout << "This application optimizes the membrane geometry parameters for a bicelle, vesicle, or double_vesicle. This doesn't run with the default slab. Use the -mp:geometry option to specify geometry." << std::endl;
				}



	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}
