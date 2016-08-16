// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief


// libRosetta headers
#include <core/types.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/pose/selection.hh>

#include <core/scoring/electron_density/util.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>

#include <basic/basic.hh>
#include <basic/database/open.hh>
#include <devel/init.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>

#include <protocols/hybridization/CartesianSampler.hh>

#include <protocols/moves/MoverContainer.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>


#include <core/io/pdb/pdb_writer.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>

#include <core/scoring/constraints/util.hh>

#include <utility/excn/Exceptions.hh>

#include <set>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//silly using/typedef


//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <basic/Tracer.hh>

using namespace core;
using namespace protocols;

using basic::T;
using basic::Error;
using basic::Warning;
using utility::vector1;



OPT_1GRP_KEY(String, cartrefine, residues)
OPT_1GRP_KEY(Boolean, cartrefine, fullatom)
OPT_1GRP_KEY(Boolean, cartrefine, bbmove)


class CartRefineWrapperMover : public protocols::moves::Mover {
public:
	CartRefineWrapperMover(){}
	void apply( core::pose::Pose & pose) {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace protocols::moves;
		using namespace protocols::hybridization;

		// parse positions
		std::set<core::Size> user_pos;
		std::string residues = option[cartrefine::residues]();
		if (residues.length()==0) {
			for (int i=1; i<=pose.total_residue(); ++i)
				user_pos.insert( i );
		} else {
			user_pos = core::pose::get_resnum_list( residues , pose );
		}

		// frag sizes
		utility::vector1< core::Size > frag_sizes;
		frag_sizes.push_back( 17 );

		core::scoring::ScoreFunctionOP fascorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "soft_rep" );
		fascorefxn->set_weight( core::scoring::cart_bonded, 0.1 );

		bool fullatom = option[cartrefine::fullatom]();
		bool bbmove = option[cartrefine::bbmove]() || option[cartrefine::fullatom]();

		CartesianSampler sampler;
		sampler.set_userpos( user_pos );
		sampler.set_strategy( "user" );
		sampler.set_rms_cutoff(3.0);
		sampler.set_fullatom(fullatom);
		sampler.set_bbmove(bbmove);
		sampler.set_nminsteps(fullatom ? 25 : 10);
		sampler.set_restore_csts(false);
		sampler.set_fa_scorefunction( fascorefxn );
		sampler.set_frag_sizes(frag_sizes);
		sampler.set_ncycles(bbmove ? 100 : 200);

		sampler.apply( pose );

		// relax (with CSTS)
		core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
		scorefxn->set_weight( core::scoring::atom_pair_constraint, 0.5 );
		scorefxn->set_weight( core::scoring::cart_bonded, 0.5 );
		protocols::relax::FastRelax relax_prot( scorefxn, 4 );
		relax_prot.min_type("lbfgs_armijo_nonmonotone");
		relax_prot.apply(pose);
	}

	virtual std::string get_name() const {
		return "CartRefineWrapperMover";
	}
};

void*
my_main( void* ) {
	using namespace protocols::moves;
	using namespace protocols::simple_moves::symmetry;

	try{
		protocols::jd2::JobDistributor::get_instance()->go( new CartRefineWrapperMover() );
	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cerr << "Exception: " << std::endl;
		excn.show( std::cerr );
	}

	return 0;
}


///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {
		NEW_OPT(cartrefine::residues, "residues to rebuild", "");
		NEW_OPT(cartrefine::fullatom, "in fullatom?", false);
		NEW_OPT(cartrefine::bbmove, "min frag backbones before graft?", false);

		devel::init(argc, argv);

		protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
