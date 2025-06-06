// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file David E Kim
/// @brief


// libRosetta headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobDistributorFactory.hh>
#include <protocols/jd2/internal_util.hh>
#include <protocols/jd2/PDBJobOutputter.hh>

#include <protocols/moves/Mover.hh>

#include <basic/Tracer.hh>
#include <devel/init.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.fwd.hh>

// Project headers

// Utility headers
#include <basic/options/keys/OptionKeys.hh>
#include <utility/excn/Exceptions.hh>



static basic::Tracer TR( "main" );

using namespace protocols::moves;
using namespace core::scoring;

class MyMover : public Mover {
public:
	MyMover();

	void apply( core::pose::Pose& pose ) override;
	std::string get_name() const override { return "ConvertToCentroidMover"; }

	MoverOP clone() const override {
		return utility::pointer::make_shared< MyMover >( *this );
	}

	MoverOP fresh_instance() const override {
		return utility::pointer::make_shared< MyMover >();
	}

};

MyMover::MyMover()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;

}

void MyMover::apply( core::pose::Pose& pose ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;

	core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID );

}

int
main( int argc, char * argv [] )
{
	try {
		using namespace protocols;
		using namespace protocols::jd2;

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core;

		jd2::register_options();

		// initialize core
		devel::init(argc, argv);

		MoverOP mymover( new MyMover );

		using namespace protocols::jd2;

		// Set up a job outputter that writes a scorefile and no PDBs and no Silent Files.
		PDBJobOutputterOP jobout( new PDBJobOutputter );

		// If the user chooses something else, then so be it, but by default score(_jd2) should only create a score
		// file and nothing else.
		protocols::jd2::JobDistributor::get_instance()->set_job_outputter( JobDistributorFactory::create_job_outputter( jobout ));

		JobDistributor::get_instance()->go( mymover );
	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}

