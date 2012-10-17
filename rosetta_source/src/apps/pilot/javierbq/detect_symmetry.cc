
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   
/// @brief  
/// @author 


#include <core/pose/Pose.hh>
#include <protocols/moves/Mover.hh>

#include <protocols/rigid/RB_geometry.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>

// Utility headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>
#include <ObjexxFCL/string.functions.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <devel/init.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

using namespace basic::options;
using namespace basic::options::OptionKeys;

static basic::Tracer TR("main");


class ThisApplication  {
public:
	ThisApplication();
	static void register_options();
private:
};

ThisApplication::ThisApplication()
{}

//OPT_1GRP_KEY( File, cluster, out )

void ThisApplication::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
  OPT(in::file::s);
}

class DetectSymmetry : public protocols::moves::Mover {
public:
	typedef core::pose::Pose Pose;
	typedef core::pose::PoseOP PoseOP;
	typedef numeric::xyzMatrix< core::Real > xyzMatrix;
	typedef numeric::xyzVector< core::Real > xyzVector;
public:
	DetectSymmetry(){ }

	virtual void apply(Pose & pose) {
		Size n_jumps = pose.num_jump();
		if( n_jumps == 0 ) return; // NO SYMMETRIC CHAINS
		Size symmetric_type = n_jumps + 1;
		TR << symmetric_type << " number of subunits found" << std::endl;
		//check that the chains have the same sequence
		std::string seq1 = pose.chain_sequence(1);
		for(int i = 2; i <= symmetric_type; i++){
			if( seq1 != pose.chain_sequence(i) ) return;
			// Perhaps check that chains are the same?? (rms == 0.0 between superposed chains)
		 }

		// Translate the center of the mass of the pose to the origin
	  xyzMatrix id_rot_mat = numeric::xyzMatrix< core::Real >::identity();
	  xyzVector cm_pose = protocols::geometry::center_of_mass(pose, 1, pose.total_residue());
		pose.apply_transform_Rx_plus_v(id_rot_mat, -1*cm_pose);
		// align the center of mass of chain A in the X axis
		// rotate around z
	  xyzVector cm_chain_A = protocols::geometry::center_of_mass(pose, 1, seq1.size());
		core::Real angle_cm_orig_x = numeric::angle_degrees(xyzVector(cm_chain_A[1] ,cm_chain_A[2] ,0), xyzVector(0,0,0), xyzVector(1,0,0));
		TR << "Angle between center of mass of chain A and x axis: " << angle_cm_orig_x << std::endl;
		xyzMatrix z_rot = numeric::z_rotation_matrix_degrees( angle_cm_orig_x);
		pose.apply_transform_Rx_plus_v(z_rot, xyzVector(0,0,0));

		// rotate around y
		core::Real angle_cm_orig_z = numeric::angle_degrees(xyzVector(cm_chain_A[1],0 ,cm_chain_A[3]), xyzVector(0,0,0), xyzVector(1,0,0));
		xyzMatrix y_rot = numeric::y_rotation_matrix_degrees( angle_cm_orig_z);
		TR << "Angle between center of mass of chain A and z axis: " << angle_cm_orig_x << std::endl;
		pose.apply_transform_Rx_plus_v(y_rot, xyzVector(0,0,0));
		
	}

	virtual std::string get_name() const {return "DetectSymmetry";}

};

using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;

int main( int argc, char** argv ) {
	ThisApplication::register_options();
	devel::init( argc, argv );
	// mover
	protocols::moves::MoverOP protocol;
	protocol = new DetectSymmetry( );

	// run
	protocols::jd2::JobDistributor::get_instance()->go( protocol );
	
}
