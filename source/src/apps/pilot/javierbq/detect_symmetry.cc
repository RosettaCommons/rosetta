// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author


#include <core/pose/Pose.hh>
#include <protocols/moves/Mover.hh>

#include <protocols/rigid/RB_geometry.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <basic/database/open.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>

#include <boost/lexical_cast.hpp>
#include <utility/exit.hh>
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

#include <utility/excn/Exceptions.hh>

using namespace basic::options;
using namespace basic::options::OptionKeys;

static THREAD_LOCAL basic::Tracer TR( "main" );


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
	DetectSymmetry():
		subunit_tolerance_( 0.01)
	{ }

	virtual void apply(Pose & pose) {
		Size n_jumps = pose.num_jump();
		if( n_jumps == 0 ) return; // NO SYMMETRIC CHAINS
		Size symmetric_type = n_jumps + 1;
		TR << symmetric_type << " number of subunits found" << std::endl;
		//check that the chains have the same sequence
		std::string seq1 = pose.chain_sequence(1);
		const Pose ref_pose( pose, 1, seq1.size() );
		for(int i = 1; i < symmetric_type; i++){
			if( seq1 != pose.chain_sequence(i) ) utility_exit_with_message("subunits have different sequence! structure can't be symmetric");
			Pose test_pose( pose, i*seq1.size()+1 , (i + 1)*seq1.size() );
			core::Real rms = core::scoring::CA_rmsd(ref_pose, test_pose);
			TR.Debug << "rms chain " << i + 1 << " " << rms << std::endl;
			if( rms > subunit_tolerance_ ) utility_exit_with_message("rmsd between subunits higher than subunit tolerance");

		}
		TR << seq1.size() << " residues per subunit" << std::endl;

		// Translate the center of the mass of the pose to the origin
	  xyzMatrix id_rot_mat = numeric::xyzMatrix< core::Real >::identity();
	  xyzVector cm_pose = core::pose::center_of_mass(pose, 1, pose.total_residue());
		pose.apply_transform_Rx_plus_v(id_rot_mat, -1*cm_pose);
		// align the center of mass of chain A in the Y axis
		// rotate around x to align the center of mass of chain a to the xy plane
	  xyzVector cm_chain_A = core::pose::center_of_mass(pose, 1, seq1.size());
		core::Real angle_cm_orig_x = numeric::angle_degrees(cm_chain_A, xyzVector(cm_chain_A[0],0,0), xyzVector(cm_chain_A[0],cm_chain_A[1], 0));
		TR.Debug << "Angle between center of mass of chain A and x axis: " << angle_cm_orig_x << std::endl;
		TR.Debug << "start - center of mass chain A " << cm_chain_A[0] << " " << cm_chain_A[1] << " " << cm_chain_A[2] << " " << std::endl;
		xyzMatrix x_rot = numeric::x_rotation_matrix_degrees(angle_cm_orig_x * ((cm_chain_A[2] < 0.0)? -1 : 1) );
		pose.apply_transform_Rx_plus_v(x_rot, xyzVector(0,0,0));
	  cm_chain_A = core::pose::center_of_mass(pose, 1, seq1.size());
		angle_cm_orig_x = numeric::angle_degrees(cm_chain_A, xyzVector(cm_chain_A[0],0,0), xyzVector(cm_chain_A[0],cm_chain_A[1], 0));
		TR.Debug << "t1 - Angle between center of mass of chain A and x axis: " << angle_cm_orig_x << std::endl;
		TR.Debug << "t1 - center of mass chain A " << cm_chain_A[0] << " " << cm_chain_A[1] << " " << cm_chain_A[2] << " " << std::endl;

		// rotate around z to move the center of mass of chain A to Y axis
		core::Real angle_cm_orig_z = numeric::angle_degrees(cm_chain_A, xyzVector(0,0,0), xyzVector(0,1,0));
		xyzMatrix z_rot = numeric::z_rotation_matrix_degrees(angle_cm_orig_z * ((cm_chain_A[0] < 0.0) ? -1 : 1) );
		pose.apply_transform_Rx_plus_v(z_rot, xyzVector(0,0,0));
	  cm_chain_A = core::pose::center_of_mass(pose, 1, seq1.size());
		TR.Debug << "t2 - Angle between center of mass of chain A and y axis: " << angle_cm_orig_z << std::endl;
		TR.Debug << "t2 - center of mass chain A " << cm_chain_A[0] << " " << cm_chain_A[1] << " " << cm_chain_A[2] << " " << std::endl;

		// rotate around y again to put the center of mass of the other subunits in the xy-plane
	  xyzVector cm_chain_B = core::pose::center_of_mass(pose, seq1.size(), 2* seq1.size());
		TR.Debug << "t4 - center of mass chain B " << cm_chain_B[0] << " " << cm_chain_B[1] << " " << cm_chain_B[2] << " " << std::endl;
		core::Real angle_cm_orig_y = numeric::angle_degrees(cm_chain_B, xyzVector(0,cm_chain_B[1],0), xyzVector(cm_chain_B[0],cm_chain_B[1], 0));
		xyzMatrix y_rot = numeric::y_rotation_matrix_degrees(angle_cm_orig_y * ((cm_chain_B[2] > 0.0) ? -1 : 1));
		TR.Debug << "t4 - Angle between center of mass of chain B and y axis: " << angle_cm_orig_y << std::endl;
		pose.apply_transform_Rx_plus_v(y_rot, xyzVector(0,0,0));
	  cm_chain_B = core::pose::center_of_mass(pose, seq1.size(), 2* seq1.size());
		TR.Debug << "t4 - center of mass chain B " << cm_chain_B[0] << " " << cm_chain_B[1] << " " << cm_chain_B[2] << " " << std::endl;

		// Now that chain A is properly oriented around the Z-axis copy it to a new pose.
		Pose new_pose( pose, 1, seq1.size() );

		// set up for symmetry
		std::string db_file = "symmetry/cyclic/C" + boost::lexical_cast< std::string >( symmetric_type ) + "_Z.sym";
		std::string path_to_symdef = basic::database::full_name(db_file);
		core::conformation::symmetry::SymmData symmdef;
		symmdef.read_symmetry_data_from_file( path_to_symdef );
		core::pose::symmetry::make_symmetric_pose( new_pose, symmdef );
		pose = new_pose;

	}

	virtual std::string get_name() const {return "DetectSymmetry";}
private:
	core::Real subunit_tolerance_;

};

using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;

int main( int argc, char** argv ) {
	try {

	ThisApplication::register_options();
	devel::init( argc, argv );
	// mover
	protocols::moves::MoverOP protocol;
	protocol = new DetectSymmetry( );

	// run
	protocols::jd2::JobDistributor::get_instance()->go( protocol );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
