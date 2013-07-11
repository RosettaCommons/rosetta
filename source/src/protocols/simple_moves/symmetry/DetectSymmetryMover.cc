// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/symmetry/DetectSymmetryMover.cc
/// @brief Automatical detection and setup of the symmetry machinery from an assymetric pose made of symmetric chains. Only works with cyclic simmetries.
/// @author Javier Castellanos ( javiercv@uw.edu ) 

#include <protocols/simple_moves/symmetry/DetectSymmetryMover.hh>
#include <protocols/simple_moves/symmetry/DetectSymmetryMoverCreator.hh>

#include <core/pose/Pose.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <basic/database/open.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/rms_util.hh>
#include <utility/tag/Tag.hh>

// Utility headers
#include <boost/lexical_cast.hpp>
#include <utility/exit.hh>
#include <basic/Tracer.hh>


#include <basic/options/option.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>


namespace protocols {
namespace simple_moves {
namespace symmetry {

static basic::Tracer TR("protocols.simple_moves.symmetry.DetectSymmetry");

// creators
std::string
DetectSymmetryMoverCreator::keyname() const {
	return DetectSymmetryMoverCreator::mover_name();
}

protocols::moves::MoverOP
DetectSymmetryMoverCreator::create_mover() const {
	return new DetectSymmetry;
}

std::string
DetectSymmetryMoverCreator::mover_name() {
	return "DetectSymmetry";
}

////////////////////
DetectSymmetry::DetectSymmetry():
	subunit_tolerance_( 0.01),
	plane_tolerance_( 1e-3 )
{ }

DetectSymmetry::DetectSymmetry(core::Real subunit_tolerance, core::Real plane_tolerance):
	subunit_tolerance_( subunit_tolerance ),
	plane_tolerance_( plane_tolerance )
{
}

void 
DetectSymmetry::apply(Pose & pose) {
	Size n_jumps = pose.num_jump();
	if( n_jumps == 0 ) utility_exit_with_message("Only one chain! no posible symmetry.");
	Size symmetric_type = n_jumps + 1;
	TR << symmetric_type << " number of subunits found" << std::endl;
	//check that the chains have the same sequence
	std::string seq1 = pose.chain_sequence(1);
	const Pose ref_pose( pose, 1, seq1.size() );
	for(Size i = 1; i < symmetric_type; ++i){
		if( seq1 != pose.chain_sequence(i) ) utility_exit_with_message("subunits have different sequence! structure can't be symmetric");
		Pose test_pose( pose, i*seq1.size()+1 , (i + 1)*seq1.size() );
		core::Real rms = core::scoring::CA_rmsd(ref_pose, test_pose);
		TR.Debug << "rms chain " << i + 1 << " " << rms << std::endl;
		if( rms > subunit_tolerance_ ) utility_exit_with_message("rmsd between subunits higher than subunit tolerance");

	}
	TR << seq1.size() << " residues per subunit" << std::endl;

	// Translate the center of the mass of the pose to the origin
  xyzMatrix id_rot_mat = numeric::xyzMatrix< core::Real >::identity();
  xyzVector cm_pose = protocols::geometry::center_of_mass(pose, 1, pose.total_residue());
	while( cm_pose[2] > plane_tolerance_) {
		pose.apply_transform_Rx_plus_v(id_rot_mat, -1*cm_pose);
  	cm_pose = protocols::geometry::center_of_mass(pose, 1, pose.total_residue());
	}

	// align the center of mass of chain A in the Y axis
	// rotate around x to align the center of mass of chain a to the xy plane
  xyzVector cm_chain_A = protocols::geometry::center_of_mass(pose, 1, seq1.size());
	core::Real angle_cm_orig_x = numeric::angle_degrees(cm_chain_A, xyzVector(cm_chain_A[0],0,0), xyzVector(cm_chain_A[0],cm_chain_A[1], 0));
	TR.Debug << "Angle between center of mass of chain A and x axis: " << angle_cm_orig_x << std::endl;
	TR.Debug << "start - center of mass chain A " << cm_chain_A[0] << " " << cm_chain_A[1] << " " << cm_chain_A[2] << " " << std::endl;
	xyzMatrix x_rot = numeric::x_rotation_matrix_degrees(angle_cm_orig_x * ((cm_chain_A[2] < 0.0)? -1 : 1) );
	pose.apply_transform_Rx_plus_v(x_rot, xyzVector(0,0,0));
  cm_chain_A = protocols::geometry::center_of_mass(pose, 1, seq1.size());
	angle_cm_orig_x = numeric::angle_degrees(cm_chain_A, xyzVector(cm_chain_A[0],0,0), xyzVector(cm_chain_A[0],cm_chain_A[1], 0));
	TR.Debug << "t1 - Angle between center of mass of chain A and x axis: " << angle_cm_orig_x << std::endl;
	TR.Debug << "t1 - center of mass chain A " << cm_chain_A[0] << " " << cm_chain_A[1] << " " << cm_chain_A[2] << " " << std::endl;

	// rotate around z to move the center of mass of chain A to Y axis
	core::Real angle_cm_orig_z = numeric::angle_degrees(cm_chain_A, xyzVector(0,0,0), xyzVector(0,1,0));
	xyzMatrix z_rot = numeric::z_rotation_matrix_degrees(angle_cm_orig_z * ((cm_chain_A[0] < 0.0) ? -1 : 1) );
	pose.apply_transform_Rx_plus_v(z_rot, xyzVector(0,0,0));
  cm_chain_A = protocols::geometry::center_of_mass(pose, 1, seq1.size());
	TR.Debug << "t2 - Angle between center of mass of chain A and y axis: " << angle_cm_orig_z << std::endl;
	TR.Debug << "t2 - center of mass chain A " << cm_chain_A[0] << " " << cm_chain_A[1] << " " << cm_chain_A[2] << " " << std::endl;

    // rotate around y again to put the center of mass of the other subunits in the xy-plane
    xyzVector cm_chain_B = protocols::geometry::center_of_mass(pose, seq1.size(), 2* seq1.size());
	core::Real angle_cm_orig_y = numeric::angle_degrees(cm_chain_B, xyzVector(0,cm_chain_B[1],0), xyzVector(cm_chain_B[0],cm_chain_B[1], 0));
	if(  cm_chain_B[ 0 ] < 0.0 )
	{
		if( cm_chain_B[2] < 0.0 ) {
			angle_cm_orig_y = 180 + (180.0 + angle_cm_orig_y);
            //	TR << "adjusting for com in -x -z" << std::endl;
		} else {
			angle_cm_orig_y = 180.0 - angle_cm_orig_y;
            //		TR << "adjusting for com in -x +z" << std::endl;
		}
	} else if( cm_chain_B[2] > 0.0 ) {
			angle_cm_orig_y = 180 + ( 180.0 - angle_cm_orig_y );
	}
        
  xyzMatrix y_rot = numeric::y_rotation_matrix_degrees(angle_cm_orig_y * ((cm_chain_B[2] > 0.0) ? -1 : 1));
	TR.Debug << "t4 - Angle between center of mass of chain B and y axis: " << angle_cm_orig_y << std::endl;
	pose.apply_transform_Rx_plus_v(y_rot, xyzVector(0,0,0));
  cm_chain_B = protocols::geometry::center_of_mass(pose, seq1.size(), 2* seq1.size());
	TR.Debug << "t4 - center of mass chain B " << cm_chain_B[0] << " " << cm_chain_B[1] << " " << cm_chain_B[2] << " " << std::endl;
    
    for( Size i = 0; i < symmetric_type; i++)
    {
        xyzVector cm_chain = protocols::geometry::center_of_mass(pose, i*seq1.size()+1, i * seq1.size() + seq1.size());
        runtime_assert_msg(cm_chain[2] > -0.0001 & cm_chain[2] < 0.0001, "com of chain " + boost::lexical_cast< std::string >( i ) +  " is not properly aligned to the x-y plane");
    }

    // if C2 check that the vector formed my the middle point between a point and the symmetric counterpart is aligned to the z-axis
    // if the aligment ia correct the middle point between any atom and its symmetric copy has to be aligned to the z-axis.
    if( symmetric_type == 2 )
    {
      xyzVector v1 = ( pose.residue( 1 ).xyz("CA") + pose.residue( seq1.size() + 1 ).xyz("CA") ) / 2;
      core::Real angle_rot_axis_z_axis = numeric::angle_degrees(v1, xyzVector(0,0,0), xyzVector(0,0,1));
      TR << v1[0] << " " << v1[1] << " " << v1[2] <<std::endl;
      TR << angle_rot_axis_z_axis << std::endl;
      xyzMatrix rot = numeric::y_rotation_matrix_degrees(  -1 * angle_rot_axis_z_axis );
      pose.apply_transform_Rx_plus_v(rot, xyzVector(0,0,0));
      for( Size i = 1, j = seq1.size()+1; i <= seq1.size(); i++,j++ )
      {
        xyzVector v = ( pose.residue( i ).xyz("CA") + pose.residue( j ).xyz("CA") )/ 2;
        runtime_assert_string_msg(v[0] < 0.0001 && v[0] > -0.0001 && v[1] < 0.0001 && v[1] > -0.0001, "rotation axis not properly aligned!");
      }
    }
  
  // Now that chain A is properly oriented around the Z-axis copy it to a new pose.
	Pose new_pose( pose, 1, seq1.size() );
  

  // set up for symmetry
	std::string db_file = "symmetry/cyclic/C" + boost::lexical_cast< std::string >( symmetric_type ) + "_Z.sym";
	std::string path_to_symdef = basic::database::full_name(db_file);
	core::conformation::symmetry::SymmData symmdef;
	symmdef.read_symmetry_data_from_file( path_to_symdef );
	// Turn symmetry hacks on
	basic::options::option[basic::options::OptionKeys::symmetry::symmetry_definition].value( db_file );
	core::pose::symmetry::make_symmetric_pose( new_pose, symmdef );
	pose = new_pose;
	
}


void
DetectSymmetry::parse_my_tag(
	TagPtr const tag,
	protocols::moves::DataMap &,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const &
)
{
	subunit_tolerance_ = tag->getOption< core::Real >("subunit_tolerance", 0.01);
	plane_tolerance_ = tag->getOption< core::Real >("plane_tolerance",1e-3);
}


} // symmetry
} // simple_moves
} // protocols
