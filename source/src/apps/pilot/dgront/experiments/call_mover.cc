#include <iostream>
#include <string>
#include <sstream>

#include <numeric/random/random.hh>

#include <protocols/jobdist/standard_mains.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/relax/ClassicRelax.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/backrub/BackrubMover.hh>

#include <devel/init.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/backrub.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/canonical_sampling.OptionKeys.gen.hh>

// Core Headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <utility/excn/Exceptions.hh>


using basic::Error;
using basic::Warning;

using namespace core;
using namespace protocols;
using namespace protocols::jobdist;
using namespace protocols::moves;

using namespace basic::options;
using namespace basic::options::OptionKeys;

void use_backrub(core::pose::PoseOP & pose,core::scoring::ScoreFunctionOP scorefxn) {

	core::pose::Pose &p(*pose);

	// set up BackrubMover and read from the database
	protocols::backrub::BackrubMover backrubmover;
	backrubmover.branchopt().read_database();

	// Setup MC
	core::Real const mc_kT( option[basic::options::OptionKeys::canonical_sampling::sampling::mc_kt] );
	protocols::moves::MonteCarlo mc(p, *scorefxn, mc_kT );
	mc.reset(p);

	// Where will it happed?
	utility::vector1<core::Size> pivot_residues;
	for ( core::Size i = 1; i <= option[ basic::options::OptionKeys::backrub::pivot_residues ].size(); ++i ) {
		if ( option[ OptionKeys::backrub::pivot_residues ][i] >= 1 ) {
			pivot_residues.push_back(option[ basic::options::OptionKeys::backrub::pivot_residues ][i]);
		}
	}
	backrubmover.set_pivot_residues(pivot_residues);

	//clear segments and set the input pose
	backrubmover.clear_segments();
	backrubmover.set_input_pose( pose );

	(*scorefxn)(p);
	Size const nstruct( option[ out::nstruct] );
	std::string move_type = backrubmover.type();
	for ( core::Size trial = 1; trial <= nstruct; ++trial ) {
		backrubmover.apply(p);
		mc.boltzmann(pose, move_type);
		std::ostringstream oss;
		oss << "p" << trial<<".pdb";
		std::cerr << (*scorefxn)(p) << "\n";
		p.dump_pdb(oss.str());
	}
}

int
main( int argc, char * argv [] )
{
	try {

		// initialize core
		devel::init(argc, argv);

		core::scoring::ScoreFunctionOP scorefxn;
		scorefxn = core::scoring::get_score_function();

		// Load a pose from a PDB file
		core::chemical::ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
		core::pose::PoseOP pose = new core::pose::Pose();
		core::pose::Pose &p(*pose);
		core::import_pose::pose_from_file( p, *rsd_set, option[ in::file::s ]().vector()[ 0 ] , core::import_pose::PDB_file);

		use_backrub(pose,scorefxn);

		std::cout << "Done! -------------------------------\n";
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

