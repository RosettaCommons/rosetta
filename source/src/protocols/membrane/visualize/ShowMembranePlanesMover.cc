// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		protocols/membrane/visualize/ShowMembranePlanesMover.cc
///
/// @brief 		Add Anchor Residues for Membrane Planes
/// @details    Add 6 virtual membrane residues to the pose as an additional
///				chain to the protein. Anchor residues are attached by jump to the membrane center
///				and are pointing downstream to allow movement of the planes. This plug in
///				works directly with the PyMol Mover
///
/// @author 	Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified: 6/27/14

// Unit Headers
#include <protocols/membrane/visualize/ShowMembranePlanesMover.hh>
#include <protocols/membrane/visualize/ShowMembranePlanesMoverCreator.hh>
#include <protocols/moves/Mover.hh>

// Package headers
#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/methods/RG_Energy_Fast.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

// Utility Headers
#include <numeric/conversions.hh>
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>

#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/membrane_new.OptionKeys.gen.hh>

// Utility Headers
#include <numeric/xyzVector.hh>
#include <utility/tag/Tag.hh>

static basic::Tracer TR( "protocols.membrane.visualize.ShowMembranePlanesMover" );

namespace protocols {
namespace membrane {
namespace visualize {

using namespace core;
using namespace core::conformation;
using namespace protocols::moves;

/////////////////////
//// Constructors ///
/////////////////////

/// @brief Default Constructor
/// @details Construct a Show membrane planes mover - uses the membrane framework
/// normal and center parameters to compute 6 residues defining an upper/lower plane
/// as the membrane planes
ShowMembranePlanesMover::ShowMembranePlanesMover() :
	Mover(),
	thickness_( 12.5 ),
	npoints_( 4 )
{
	register_options();
	init_from_cmd();
}


/// @brief Constructor
/// @details Construct a Show membrane planes mover - uses the membrane framework
/// normal and center parameters to compute 6 residues defining an upper/lower plane
/// as the membrane planes. Can specify membrane thicnkess and/or raidus of gyration to define
/// the radii of the planes (triangle shaped)
ShowMembranePlanesMover::ShowMembranePlanesMover(
	core::Real thickness,
	core::Size npoints
	) :
	thickness_( thickness ),
	npoints_( npoints )
{}

/// @brief Copy Constructor
/// @details Create a deep copy of the show membrane planes mover
ShowMembranePlanesMover::ShowMembranePlanesMover( ShowMembranePlanesMover const & src ) :
	Mover( src ),
	thickness_( src.thickness_ ),
	npoints_( src.npoints_ )
{}

/// @brief Assignemnt Operator
/// @details Create a deep copy of the Show membrane planes mover while
/// overloading the assignment operator
ShowMembranePlanesMover &
ShowMembranePlanesMover::operator=( ShowMembranePlanesMover const & src ) {
	
	// Abort self-assignment.
	if (this == &src) {
		return *this;
	}
	
	// Otherwise, create a new object
	return *( new ShowMembranePlanesMover( *this ) );
	
}


/// @brief Destructor - Get rid of this class
ShowMembranePlanesMover::~ShowMembranePlanesMover() {}

//////////////////////
//// Mover Methods ///
//////////////////////

/// @brief Get the name of this mover (ShowMembranePlanesMover)
std::string
ShowMembranePlanesMover::get_name() const {
	return "ShowMembranePlanesMover";
}

/// @brief Add Anchor Residues to Pose
/// @details Add 6 anchor residues to pose defining the bounds of the membrane planes.
/// Width of membrane planes are determined by 2x RG of the pose.
void
ShowMembranePlanesMover::apply( Pose & pose ) {
	
	using namespace core::scoring::methods;
	using namespace protocols::membrane;
	
	// Check that my pose is a membrane pose
	if (! pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Cannot show membrane planes on a non membrane pose!" );
	}
	
	TR << "Adding virtual residues defining the bounds of membrane planes to the pose for visualization in pymol" << std::endl;
	
	// Grab the residue typeset of the pose
	bool fullatom = pose.is_fullatom();
	
	// Grab the membrane center & normal parameters
	Vector center( pose.conformation().membrane_info()->membrane_center() );
	Vector normal( pose.conformation().membrane_info()->membrane_normal() );
	normal.normalize();
	
	// Grab a list of relevant residues and go
	// Compute radius of gyration of the pose
	utility::vector1< bool > relevant_residues;
	relevant_residues.resize( pose.total_residue() );
	
	// I thought there was a fill method or better templated method, but couldn't find it
	for ( Size i = 1; i < relevant_residues.size(); ++i ) {
		relevant_residues[i] = true;
	}
	
	// Define radius of gyration
	RG_Energy_Fast rg_method;
	core::Real rg =  2*rg_method.calculate_rg_score( pose, relevant_residues );

	
	// Determine the plane points by projecting onto the normal vector
	Vector upper = center + normal*thickness_;
	Vector lower = center - normal*thickness_;
	
	// Compute Upper points
	utility::vector1< Vector > upper_points = select_plane_points( upper, normal, npoints_, rg );
	utility::vector1< Vector > lower_points = select_plane_points( lower, normal, npoints_, rg );
	
	// Attach plane anchor residues by jump to the membrane virtual residue
	core::Size membrane_rsd = pose.conformation().membrane_info()->membrane_rsd_num();
	
	// Attach and track residues
	utility::vector1< Size > top_points;
	utility::vector1< Size > bottom_points;
	
	// Attach top plane points
	for ( Size i = 1; i <= npoints_; ++i ) {
		
		ResidueOP upper = create_membrane_virtual( upper_points[i], fullatom );
		
		if ( i == 1 ) {
			pose.append_residue_by_jump( *upper, membrane_rsd, "", "", true );
			top_points.push_back( pose.total_residue() );
		} else {
			pose.append_residue_by_jump( *upper, membrane_rsd, "", "", false );
			top_points.push_back( pose.total_residue() );
		}
	}
	
	// Attach bottom points
	for ( Size i = 1; i <= npoints_; ++i ) {
		
		ResidueOP lower = create_membrane_virtual( lower_points[i], fullatom );
		
		if ( i == 1 ) {
			pose.append_residue_by_jump( *lower, membrane_rsd, "", "", true );
			bottom_points.push_back( pose.total_residue() );
		} else {
			pose.append_residue_by_jump( *lower, membrane_rsd, "", "", false );
			bottom_points.push_back( pose.total_residue() );
		}
	}
	
	// Track plane residues in membrane info
	pose.conformation().membrane_info()->setup_plane_visualization( top_points, bottom_points );
	
}

////////////////////////////////
//// Rosetta Scripts Methods ///
////////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
ShowMembranePlanesMover::clone() const {
	return new ShowMembranePlanesMover( *this );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
ShowMembranePlanesMover::fresh_instance() const {
	return new ShowMembranePlanesMover();
}

/// @brief Pass Rosetta Scripts Options for this Mover
void
ShowMembranePlanesMover::parse_my_tag(
  utility::tag::TagCOP tag,
  basic::datacache::DataMap &,
  protocols::filters::Filters_map const &,
  protocols::moves::Movers_map const &,
  core::pose::Pose const &
  ) {
	
	if ( tag->hasOption( "thickness" ) ) {
		thickness_ = tag->getOption< Real >( "thickness" );
	}
	

	if ( tag->hasOption( "num_points" ) ) {
		npoints_ = tag->getOption< Size >( "num_points" );
		
		if ( npoints_ < 3 ) {
			utility_exit_with_message( "Cannot define planes with less than 3 points!" );
		}
	}
	
}


/// @brief Create a new copy of this mover
protocols::moves::MoverOP
ShowMembranePlanesMoverCreator::create_mover() const {
	return new ShowMembranePlanesMover();
}

/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
ShowMembranePlanesMoverCreator::keyname() const {
	return ShowMembranePlanesMoverCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
ShowMembranePlanesMoverCreator::mover_name() {
	return "ShowMembranePlanesMover";
}

//////////////////////
/// Helper Methods ///
//////////////////////

/// @brief Register Options with JD2
void
ShowMembranePlanesMover::register_options() {
	
	using namespace basic::options;
	
	option.add_relevant( OptionKeys::membrane_new::viewer::thickness );
	option.add_relevant( OptionKeys::membrane_new::viewer::num_points );
	
}

/// @brief Initialize Options from the Command Line
/// @details Options allowed are vrt spacing and plane width
void
ShowMembranePlanesMover::init_from_cmd() {
	
	using namespace basic::options;
	
	// Specify membrane thickness option
	if ( option[ OptionKeys::membrane_new::visualize::thickness ].user() ) {
		thickness_ = option[ OptionKeys::membrane_new::visualize::thickness ]();
	}
		
	// Specify number of points to define the membrane plane
	if ( option[ OptionKeys::membrane_new::viewer::num_points ].user() ) {
		npoints_ = option[ OptionKeys::membrane_new::viewer::num_points ]();
		
		if ( npoints_ < 3 ) {
			utility_exit_with_message( "Cannot define planes with less than 3 points!" );
		}
	}
	
}

/// @brief Create a Membrane Residue
/// @details Given a centered position and residue typeset, return
/// a ResidueOP with the xyz coordinate pos, type MEM, from typeset given
ResidueOP
ShowMembranePlanesMover::create_membrane_virtual( Vector pos, bool fullatom ) {
	
	using namespace core::conformation;
	using namespace core::chemical;
	
	// Center position
	int center = 2;
	
	// Grab the current residue typeset and create a new residue
	ResidueTypeSetCAP const & residue_set(
										  core::chemical::ChemicalManager::get_instance()->residue_type_set( fullatom ? core::chemical::FA_STANDARD : core::chemical::CENTROID )
										  );
	
	// Create a new Residue from rsd typeset of type MEM
	ResidueTypeCOPs const & rsd_type_list( residue_set->name3_map( "MEM" ) );
	ResidueType const & membrane( *rsd_type_list[1] );
	ResidueOP rsd( ResidueFactory::create_residue( membrane ) );
	
	// Fill the Residue with normal/center info
	rsd->set_xyz( center, pos );
	return rsd;
	
}

/// @brief Select Points to define the upper and lower planes
utility::vector1< Vector >
ShowMembranePlanesMover::select_plane_points(
											 core::Vector center,
											 core::Vector normal,
											 core::Size n_points,
											 core::Real radius
											 ) {
	
	using namespace numeric::constants;
	using namespace numeric::conversions;
	
	// Pick point angles based on number of points to define the CGO plane
	utility::vector1< core::Real > angles;
	for ( core::Size i = 0; i <= (n_points-1); ++i ) {
		Real angle = ( (i + 1) * 2 * f::pi )/n_points;
		angles.push_back( angle );
	}
	
	// Pick an arbitrary orthogonal Unit Vector
	core::Real tolerance = pow(10, -7);
	core::Vector p( 0, 0, 0 );
	if ( abs(normal.x() + normal.y()) < tolerance ) {
		p.x() = -normal.y() - normal.z();
		p.y() = normal.x();
		p.z() = normal.x();
	} else {
		p.x() = normal.z();
		p.y() = normal.z();
		p.z() = -normal.x() - normal.y();
	}
	
	// Scale point by radius of gyration
	p = 2 * radius * p;
	
	// For the remaining angles, rotate p around the normal vector
	utility::vector1< Vector > points;
	for ( Size i = 1; i <= angles.size(); ++i ) {
		
		numeric::xyzMatrix< Real > rot_matrix = rotation_matrix_radians( normal, angles[i] );
		Vector new_point = center + rot_matrix*p;
		points.push_back( new_point );
		
	}
	
	return points;
}



} // visualize
} // membrane
} // protocols


