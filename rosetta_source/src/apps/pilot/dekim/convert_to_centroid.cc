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
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/scoring/ResidualDipolarCoupling.hh>

// Project headers
#include <core/scoring/rms_util.hh>
#include <numeric/model_quality/rms.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/RelaxProtocolBase.hh>
#include <protocols/relax/util.hh>
#include <protocols/protein_interface_design/dock_design_filters.hh>
#include <protocols/simple_filters/RmsdEvaluator.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/exit.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/string_util.hh>

#include <sstream>
#include <fstream>



static basic::Tracer TR("main");

using namespace protocols::moves;
using namespace core::scoring;
// using namespace basic::options;
// using namespace basic::options::OptionKeys;



class MyScoreMover : public Mover {
public:
	MyScoreMover();

	virtual void apply( core::pose::Pose& pose );
	std::string get_name() const { return "NonLocalFragsScoreMover"; }

	virtual MoverOP clone() const {
		return new MyScoreMover( *this );
	}

	virtual	MoverOP	fresh_instance() const {
		return new MyScoreMover;
	}

private:
	core::scoring::ScoreFunctionOP sfxn_;

};

MyScoreMover::MyScoreMover()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;

}

void MyScoreMover::apply( core::pose::Pose& pose ) {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
	using namespace core;

  core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID );

}

int
main( int argc, char * argv [] )
{
	using namespace protocols;
	using namespace protocols::jd2;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;

	jd2::register_options();

	// initialize core
	devel::init(argc, argv);

	//MyScoreMover* scoremover = new MyScoreMover;
	MoverOP scoremover = new MyScoreMover;

	using namespace protocols::jd2;

	// Make sure the default JobOutputter is SilentJobOutputter to ensure that when score_jd2
	// is called with default arguments is prints a proper scorefile and not the hacky thing that
	// the  JobOutputter scorefile() function produces (which for example skips Evaluators!!)

	// Set up a job outputter that writes a scorefile and no PDBs and no Silent Files.
	PDBJobOutputterOP jobout = new PDBJobOutputter;

	// If the user chooses something else, then so be it, but by default score(_jd2) should only create a score
	// file and nothing else.
	protocols::jd2::JobDistributor::get_instance()->set_job_outputter( JobDistributorFactory::create_job_outputter( jobout ));

	try{
		JobDistributor::get_instance()->go( scoremover );
	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cerr << "Exception: " << std::endl;
		excn.show( std::cerr );
		std::cout << "Exception: " << std::endl;
		excn.show( std::cout ); //so its also seen in a >LOG file
	}
	return 0;
}

