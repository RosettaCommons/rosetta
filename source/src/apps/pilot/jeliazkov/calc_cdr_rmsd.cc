// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/apps/pilot/jeliazkov/calc_cdr_rmsd.cc
/// @brief pilot app for two cdr rms + LHOC metric calculation
/// @author Jeliazko Jeliazkov

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/pose/Pose.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

#include <devel/init.hh>

#include <protocols/antibody/metrics.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/moves/Mover.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.fwd.hh>
#include <protocols/jd2/internal_util.hh>

#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>


static basic::Tracer TR( "apps.pilot.jeliazkov.calc_cdr_rmsd" );


class CalcCdrRms : public protocols::moves::Mover {
public:
	CalcCdrRms() : Mover() {
		using namespace basic::options;
		using namespace core;
		// set native otherwise its meaningless to run this code
		// if native is not given, what will the error be?
		import_pose::pose_from_file( native_pose_, option[ OptionKeys::in::file::native ](), import_pose::PDB_file );
		// get antibody info as well
		native_abi_ = protocols::antibody::AntibodyInfoOP ( new protocols::antibody::AntibodyInfo( native_pose_ ) );
	}
	// destructor
	~CalcCdrRms() override= default;

	void apply( core::pose::Pose & pose_in ) override
	{
		TR << "Comparing CDR RMS using old (AntibodyInfo) and new (JRJ) approach..." << std::endl;

		using namespace core;
		using namespace utility;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace protocols::jd2;
		using namespace protocols::moves;
		using namespace protocols::antibody;

		JobOP job( JobDistributor::get_instance()->current_job() );

		AntibodyInfoOP ab_info = AntibodyInfoOP ( new AntibodyInfo( pose_in ) );

		// calculate rmsds with my new methods
		vector1< Real > results =
			cdr_backbone_rmsds( pose_in, native_pose_, ab_info, native_abi_ );

		setPoseExtraScore( pose_in, "OCD",  results[1]);
		setPoseExtraScore( pose_in, "FRH",  results[2]);
		setPoseExtraScore( pose_in, "H1_new",  results[3]);
		setPoseExtraScore( pose_in, "H2_new",  results[4]);
		setPoseExtraScore( pose_in, "H3_new",  results[5]);
		setPoseExtraScore( pose_in, "FRL",  results[6]);
		setPoseExtraScore( pose_in, "L1_new",  results[7]);
		setPoseExtraScore( pose_in, "L2_new",  results[8]);
		setPoseExtraScore( pose_in, "L3_new",  results[9]);

		// repeat with other method
		/*
		align_to_native( pose_in, native_pose_, ab_info, native_abi_, "H" );
		for ( auto cdr : {h1, h2, h3} ) {
		protocols::loops::LoopsOP loops = ab_info->get_CDR_in_loopsop(cdr);
		Real rms = global_loop_rmsd( pose_in, native_pose_, loops );
		setPoseExtraScore( pose_in, ab_info->get_CDR_name( cdr )  + "_old", rms );
		}

		align_to_native( pose_in, native_pose_, ab_info, native_abi_, "L" );
		for ( auto cdr : {l1, l2, l3} ) {
		protocols::loops::LoopsOP loops = ab_info->get_CDR_in_loopsop( cdr );
		Real rms = global_loop_rmsd( pose_in, native_pose_, loops );
		setPoseExtraScore( pose_in, ab_info->get_CDR_name( cdr ) + "_old", rms );
		}
		*/
		ScoreFunctionOP sf( get_score_function() );
		( *sf )( pose_in );

		TR << "Finished!" << std::endl;

		return;
	} // PackingAngle::apply()

	std::string get_name() const override { return "CalcCdrRms"; }


	protocols::moves::MoverOP
	fresh_instance() const override {
		return utility::pointer::make_shared< CalcCdrRms >();
	}


	bool
	reinitialize_for_each_job() const override { return false; }


	bool
	reinitialize_for_new_input() const override { return false; }

private:

	core::pose::Pose native_pose_;
	protocols::antibody::AntibodyInfoOP native_abi_;

};

using CalcCdrRmsOP = utility::pointer::shared_ptr< CalcCdrRms >;


/////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

		using namespace protocols::jd2;

		register_options();

		// initialize core
		devel::init(argc, argv);

		CalcCdrRmsOP CalcCdrRmsInst = utility::pointer::make_shared< CalcCdrRms >();
		JobDistributor::get_instance()->go( CalcCdrRmsInst );

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}
