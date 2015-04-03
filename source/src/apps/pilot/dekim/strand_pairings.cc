// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file David E Kim
/// @brief


// libRosetta headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobDistributorFactory.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/PDBJobOutputter.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>

#include <basic/Tracer.hh>
#include <devel/init.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/Pose.hh>

#include <core/scoring/dssp/StrandPairing.hh>
#include <core/scoring/dssp/PairingsList.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/VariantType.hh>

// Project headers

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/exit.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>

#include <sstream>
#include <fstream>


static thread_local basic::Tracer TR( "main" );

using namespace protocols::moves;
using namespace core::scoring;

class MyMover : public Mover {
public:
	MyMover();

	virtual void apply( core::pose::Pose& pose );
	std::string get_name() const { return "StrandPairingsMover"; }

	virtual MoverOP clone() const {
		return MoverOP( new MyMover( *this ) );
	}

	virtual	MoverOP	fresh_instance() const {
		return MoverOP( new MyMover );
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


	core::scoring::dssp::StrandPairingSet strand_pairings( pose );
	TR.Info << strand_pairings << std::endl;

  for ( core::scoring::dssp::StrandPairingSet::const_iterator it = strand_pairings.begin();
	       it != strand_pairings.end();  it++ ) {
		core::scoring::dssp::PairingList beta_pairs;

		it->get_beta_pairs( beta_pairs );
		for ( core::scoring::dssp::PairingList::const_iterator iti = beta_pairs.begin(), eiti = beta_pairs.end();
		        iti != eiti; ++iti ) {
			TR.Info << *iti << std::endl;
		}
	}

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
	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cerr << "Exception: " << std::endl;
		excn.show( std::cerr );
		std::cout << "Exception: " << std::endl;
		excn.show( std::cout ); //so its also seen in a >LOG file
		return -1;
	}
	return 0;
}

