// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/smlewis/sasa_qd.cc
/// @brief Q&D protocol to run SasaCalculatorLegacy as protocol (why is this not already there?)
/// @author Steven Lewis

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
//#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>

#include <basic/MetricValue.hh>
#include <core/pose/metrics/CalculatorFactory.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("apps.pilot.smlewis.sasa_qd");

class sasa_qdMover : public protocols::moves::Mover {
public:
	sasa_qdMover() : Sasa_("sasaqd") {

		using core::pose::metrics::CalculatorFactory;
		//create the SasaCalculatorLegacy in the constructor (to ensure it will always exist)
		if( CalculatorFactory::Instance().check_calculator_exists( Sasa_ ) ){
			Warning() << "In sasa_qd, calculator " << Sasa_
			<< " already exists, this is hopefully correct for your purposes" << std::endl;
		} else {
			CalculatorFactory::Instance().register_calculator( Sasa_, new core::pose::metrics::simple_calculators::SasaCalculatorLegacy);
		}

	}

	virtual
	void
	apply(core::pose::Pose & pose ){

		//apply the calculator and extract the values
		basic::MetricValue< utility::vector1< core::Real > > by_residue_sasa_m;
		pose.metric( Sasa_, "residue_sasa" /*keystring for per-residue sasa vector*/, by_residue_sasa_m);
		utility::vector1< core::Real > const by_residue_sasa(by_residue_sasa_m.value());

		//print the values
		core::Size const nres(pose.total_residue());
		runtime_assert(nres == by_residue_sasa.size());
		//protocols::jd2::JobOP job_me( JobDistributor::get_instance()->current_job() );
		//std::string me(JobDistributor::get_instance()->job_outputter()->output_name( job_me ) );
		std::string const jobname(protocols::jd2::JobDistributor::get_instance()->current_job()->input_tag());

		TR << jobname << " PDBres PDBchain residue restype SASA" << std::endl;
		for(core::Size i(1); i<=nres; ++i){
			TR
				<< jobname << " "
				<< pose.pdb_info()->pose2pdb(i) << " "
				//<< pose.pdb_info()->chain(i) << " " //not needed - pose2pdb does both
				<< i << " "
				<< pose.residue_type(i).name1() << " "
				<< by_residue_sasa[i] << std::endl;
		}

		set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY); //hacky way to prevent output (we already output to tracer)
		return;
	}

	virtual
	std::string
	get_name() const { return "sasa_qdMover"; }

	std::string const Sasa_;

};

int
main( int argc, char* argv[] )
{

	try {

	devel::init(argc, argv);

	protocols::jd2::JobDistributor::get_instance()->go(new sasa_qdMover());

	TR << "************************d**o**n**e**************************************" << std::endl;

	return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

}
