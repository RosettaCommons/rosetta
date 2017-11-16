// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/apps/pilot/nmarze/tilt_angle.cc
/// @brief pilot app for extracting VL/VH packing angle from pdb/list of pdbs
/// @author Nick Marze (nickmarze@gmail.com)


#include <protocols/jd2/util.hh>
#include <protocols/jd2/internal_util.hh>
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
#include <core/pose/PDBInfo.hh>
#include <protocols/moves/Mover.hh>
#include <utility/vector1.hh>
#include <core/pose/symmetry/util.hh>
#include <protocols/rigid/RB_geometry.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobOutputter.hh>

static basic::Tracer TR( "apps.pilot.nmarze.hemagg_trimer" );


class HemaggTrimer : public protocols::moves::Mover {
public:
	// default constructor
	HemaggTrimer(){};
	// destructor
	virtual ~HemaggTrimer(){};

	virtual void apply( core::pose::Pose & pose_in )
	{
		TR << "Applying HemaggTrimer" << std::endl;

		using namespace core;
		using namespace protocols;
		using namespace core::pose;
		using namespace protocols::moves;

		protocols::jd2::JobOP job( protocols::jd2::JobDistributor::get_instance()->current_job() );

		core::Size com = protocols::geometry::residue_center_of_mass( pose_in, 1, pose_in.size() );

		core::pose::symmetry::make_symmetric_pose( pose_in );

		TR << "Center of mass residue = " << com << std::endl;

		TR << "Finished applying HemaggTrimer" << std::endl;

		return;
	} // HemaggTrimer::apply()

	std::string get_name() const { return "PackingAngle"; }

	virtual
	protocols::moves::MoverOP
	fresh_instance() const {
		return new HemaggTrimer;
	}

	virtual
	bool
	reinitialize_for_each_job() const { return false; }

	virtual
	bool
	reinitialize_for_new_input() const { return false; }

private:

};

typedef utility::pointer::owning_ptr< HemaggTrimer > HemaggTrimerOP;


/////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

		protocols::jd2::register_options();

		// initialize core
		devel::init(argc, argv);

		HemaggTrimerOP hemagg_trimer = new HemaggTrimer;
		protocols::jd2::JobDistributor::get_instance()->go( hemagg_trimer );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

