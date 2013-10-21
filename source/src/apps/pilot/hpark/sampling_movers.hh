#ifndef INCLUDED_apps_pilot_hpark_sampling_movers_hh
#define INCLUDED_apps_pilot_hpark_sampling_movers_hh

#include <apps/pilot/hpark/sampling_utils.hh>
#include <apps/pilot/hpark/sampling_movers.hh>

#include <protocols/normalmode/NormalModeRelaxMover.hh>
#include <protocols/normalmode/NormalModeRelaxMover.fwd.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/md/CartesianMD.hh>
#include <protocols/simple_moves/CombinePoseMover.hh>
#include <protocols/loophash/LoopHashSampler.fwd.hh>
#include <protocols/loophash/LocalInserter.fwd.hh>
#include <protocols/loophash/LoopHashLibrary.fwd.hh>

#include <core/pose/Pose.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/io/silent/SilentStruct.fwd.hh>

#include <core/optimization/Minimizer.fwd.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>

#include <numeric/random/random.hh>
#include <utility/vector1.hh>

namespace myspace {

using namespace core;

utility::vector1< pose::Pose >
test_NMrelaxer( pose::Pose const pose,
		scoring::ScoreFunctionCOP sfxn,
		optimization::MinimizerOptionsCOP minoption,
		Size const nstruct,
		bool const cartesian,
		bool const repack
		);

utility::vector1< pose::Pose >
test_MD( pose::Pose const pose,
	 scoring::ScoreFunctionCOP sfxn,
	 scoring::ScoreFunctionCOP sfxn_obj,
	 optimization::MinimizerOptionsCOP minoption,
	 Size const nstruct,
	 Size const nper_trj
	 );

utility::vector1< pose::Pose >
test_relax( pose::Pose const pose,
	    scoring::ScoreFunctionCOP sfxn,
	    Size const nstruct,
	    bool const cartesian
	    );

// full Repack & Cartmin
utility::vector1< pose::Pose >
test_recombine( pose::Pose const pose,
		pose::Pose const pose2,
		scoring::ScoreFunctionCOP sfxn,
		optimization::MinimizerOptionsCOP minoption,
		Size const nstruct,
		bool const repack
		);

// by default run torsion minimization
utility::vector1< pose::Pose >
test_bbgauss( pose::Pose const pose,
	      scoring::ScoreFunctionCOP sfxn,
	      optimization::MinimizerOptionsCOP minoption,
	      Size const nstruct
	      );

// Repack & Cartmin
utility::vector1< pose::Pose >
test_loophash( pose::Pose const pose,
	       scoring::ScoreFunctionCOP sfxn,
	       optimization::MinimizerOptionsCOP minoptions,
	       Size const nstruct,
	       bool const local
	       );

}

#endif
