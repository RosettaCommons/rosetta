// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/apps/pilot/nmarze/tilt_angle.cc
/// @brief pilot app for extracting VL/VH packing angle from pdb/list of pdbs
/// @author Nick Marze (nickmarze@gmail.com)


#include <protocols/jd2/util.hh>
#include <utility/excn/Exceptions.hh>
//#include <utility/exit.hh>

#include <devel/init.hh>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/exit.hh>
#include <protocols/antibody/metrics.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <core/pose/PDBInfo.hh>
#include <protocols/moves/Mover.hh>
#include <utility/vector1.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobOutputter.hh>

static thread_local basic::Tracer TR( "apps.pilot.nmarze.packing_angle" );


class PackingAngle : public protocols::moves::Mover {
public:
	// default constructor
	PackingAngle(){};
	// destructor
	virtual ~PackingAngle(){};

virtual void apply( core::pose::Pose & pose_in )
{
	TR << "Applying Packing Angle Calculator" << std::endl;

	using namespace core;
	using namespace protocols;
	using namespace core::pose;
	using namespace protocols::moves;
	using namespace protocols::antibody;

	protocols::jd2::JobOP job( protocols::jd2::JobDistributor::get_instance()->current_job() );

	antibody::AntibodyInfoCOP ab_info_ = antibody::AntibodyInfoCOP ( new AntibodyInfo( pose_in ) );
	PoseCOP new_pose = PoseOP( new Pose( pose_in ) );

	vector1< Real > orientation_coords_ = vl_vh_orientation_coords( *new_pose , *ab_info_ );

	job->add_string_real_pair( "VL_VH_distance", orientation_coords_[1] );
	job->add_string_real_pair( "VL_VH_opening_angle", orientation_coords_[2] );
	job->add_string_real_pair( "VL_VH_opposite_opening_angle", orientation_coords_[3] );
	job->add_string_real_pair( "VL_VH_packing_angle", orientation_coords_[4] );
	job->add_string_real_pair( "H1_length", ab_info_->get_CDR_length( h1 ) );
	job->add_string_real_pair( "H2_length", ab_info_->get_CDR_length( h2 ) );
	job->add_string_real_pair( "H3_length", ab_info_->get_CDR_length( h3 ) );
	job->add_string_real_pair( "L1_length", ab_info_->get_CDR_length( l1 ) );
	job->add_string_real_pair( "L2_length", ab_info_->get_CDR_length( l2 ) );
	job->add_string_real_pair( "L3_length", ab_info_->get_CDR_length( l3 ) );

	TR << "Finished applying Packing Angle Calculator" << std::endl;

	return;
} // PackingAngle::apply()

std::string get_name() const { return "PackingAngle"; }

  virtual
  protocols::moves::MoverOP
  fresh_instance() const {
	  return protocols::moves::MoverOP( new PackingAngle() );
  }

  virtual
  bool
  reinitialize_for_each_job() const { return false; }

  virtual
  bool
  reinitialize_for_new_input() const { return false; }

private:

};

typedef utility::pointer::shared_ptr< PackingAngle > PackingAngleOP;


/////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
    try {

		protocols::jd2::register_options();

		// initialize core
		devel::init(argc, argv);

		PackingAngleOP packing_angle = PackingAngleOP( new PackingAngle );
		protocols::jd2::JobDistributor::get_instance()->go( packing_angle );

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
				return -1;
    }
    return 0;
}
