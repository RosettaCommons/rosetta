#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/relax/util.hh>
#include <utility/excn/Exceptions.hh>

#include <devel/init.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/relax/RelaxProtocolBase.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/simple_moves/BBGaussianMover.hh>

//options
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/bbg.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

//
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>



static thread_local basic::Tracer TR( "fpd_bbg" );

OPT_1GRP_KEY(Integer, bbg, ntrials)
OPT_1GRP_KEY(Integer, bbg, nstride)
OPT_1GRP_KEY(Real, bbg, kT)


///////////////////////////////////////////////////////////////////////////////

class BBGWrapperMover : public protocols::moves::Mover {
public:
	BBGWrapperMover() {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		// options[]
    	ntrials = option[bbg::ntrials];
		kT = option[bbg::kT];

		// should be options[]
		scorefxn_low = core::scoring::ScoreFunctionFactory::create_score_function("cen_std");
		scorefxn_high = core::scoring::get_score_function();
	}

	virtual std::string get_name() const { return "BBGWrapperMover";}


	void apply( Pose & pose ) {
		protocols::simple_moves::SwitchResidueTypeSetMover to_centroid("centroid");
		core::Size n_res = pose.n_residue();

		//mc
		TR << "MC mover" << std::endl;
		to_centroid.apply(pose);
		protocols::moves::MonteCarlo mc(pose, *scorefxn_low, kT);

		protocols::viewer::add_monte_carlo_viewer(mc, "Gaussian", 600, 600);

		protocols::simple_moves::BBG8T3AMover bbg8t3amover;

		//ref pose
		core::pose::Pose ref_pose(pose);

		std::string move_type = bbg8t3amover.type();
		core::Real proposal_density_ratio=1.0;
		for (Size i = 1; i <= ntrials; ++i) {
			bbg8t3amover.apply(pose);
			proposal_density_ratio = bbg8t3amover.last_proposal_density_ratio();

			//TR << "proposal density ratio = " << proposal_density_ratio << endl;
			mc.boltzmann(pose, move_type, proposal_density_ratio);
		}

		mc.show_counters();
	}

private:
    core::Size ntrials;
	core::Real kT;
    core::scoring::ScoreFunctionOP scorefxn_low;
    core::scoring::ScoreFunctionOP scorefxn_high;
};

///////////////////////////////////////////////////////////////////////////////

void*
my_main( void* ) {
	using namespace protocols::moves;

	SequenceMoverOP seq( new SequenceMover() );
	//seq->add_mover( new symmetry::SetupForSymmetryMover() );
	seq->add_mover( new BBGWrapperMover() );
	seq->add_mover( new protocols::simple_moves::SwitchResidueTypeSetMover("fa_standard") );
	protocols::relax::RelaxProtocolBaseOP fa_rlx( protocols::relax::generate_relax_from_cmd() );
	seq->add_mover( fa_rlx );

	try{
		protocols::jd2::JobDistributor::get_instance()->go( seq );
	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cerr << "Exception: " << std::endl;
		excn.show( std::cerr );
	}

	return 0;
}


int main(int argc, char *argv[])
{
    try {

    NEW_OPT(bbg::ntrials, "number of Monte Carlo trials to run", 1000);
    NEW_OPT(bbg::kT, "temperature of MC mover", 0.6);

    devel::init(argc, argv);
    protocols::viewer::viewer_main( my_main );
    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }
    return 0;
}

