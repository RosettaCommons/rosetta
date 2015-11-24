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
#include <numeric/random/random_permutation.hh>
#include <numeric/random/random.hh>

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


#include <iostream>
#include <string>
#include <fstream>
#include <sstream>


static THREAD_LOCAL basic::Tracer TR( "fpd_bbg" );

OPT_1GRP_KEY(Real, perturb, mutrate)


///////////////////////////////////////////////////////////////////////////////

class PerturbStruct : public protocols::moves::Mover {
private:
	core::Real mutrate_;
    core::scoring::ScoreFunctionOP scorefxn_;

public:
	PerturbStruct() {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

    	mutrate_ = option[perturb::mutrate];
		scorefxn_ = core::scoring::get_score_function();
	}

	virtual std::string get_name() const { return "PerturbStruct";}

	//////
	void do_mutate( core::pose::Pose & pose, core::Size nmut ) {
		utility::vector1< core::Size > positions;
		for (core::Size i=1; i<=pose.total_residue(); ++i) {
			core::chemical::ResidueType const & cur_restype = pose.residue_type( i );
			if (!cur_restype.is_protein()) continue;
			if (cur_restype.aa() == core::chemical::aa_cys && cur_restype.has_variant_type( core::chemical::DISULFIDE )) continue;
			positions.push_back(i);
		}
		numeric::random::random_permutation(positions, numeric::random::rg());

		for (core::Size i=1; i<=std::min(nmut,pose.total_residue()); ++i) {
			core::Size res_to_mut = positions[i];

			//
			core::chemical::AA to_mutate = (core::chemical::AA) numeric::random::rg().random_range(1,20);
			while (to_mutate == pose.residue(res_to_mut).aa()) {
				to_mutate = (core::chemical::AA) numeric::random::rg().random_range(1,20);
			}

			// random AA
			core::chemical::ResidueTypeSet const & rsd_set( pose.residue(res_to_mut).residue_type_set() );
			core::chemical::ResidueTypeCOP rtype = rsd_set.get_representative_type_aa(to_mutate);
			core::conformation::Residue replace_res( *rtype, true );

			pose.replace_residue( res_to_mut, replace_res, true);
		}
	}

	//////
	void revert_muations( core::pose::Pose & pose, core::pose::Pose const & ref ) {
		for (core::Size i=1; i<=pose.total_residue(); ++i) {
			if (pose.residue(i).aa() != ref.residue(i).aa()) {
				pose.replace_residue( i, ref.residue(i), true);
			}
		}
	}


	void apply( core::pose::Pose & pose ) {
		core::pose::Pose ref_pose(pose);

		// random mutations
		core::Size nmut = (core::Size) std::floor( mutrate_ * pose.total_residue() + 0.5 );

		protocols::relax::RelaxProtocolBaseOP fa_rlx( protocols::relax::generate_relax_from_cmd() );

		do_mutate( pose, nmut );
		fa_rlx->apply( pose );
		revert_muations( pose, ref_pose );
		fa_rlx->apply( pose );
	}

};

///////////////////////////////////////////////////////////////////////////////

void*
my_main( void* ) {
	using namespace protocols::moves;

	MoverOP seq( new PerturbStruct() );

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

    NEW_OPT(perturb::mutrate, "mutrate", 0.1);

    devel::init(argc, argv);
    protocols::viewer::viewer_main( my_main );
    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }
    return 0;
}

