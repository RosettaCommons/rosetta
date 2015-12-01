// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/LK_BallEnergy.hh>
#include <core/scoring/sasa.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.fwd.hh>

#include <core/scoring/electron_density/util.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>

#include <basic/basic.hh>
#include <basic/database/open.hh>
#include <devel/init.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>


#include <core/io/pdb/pose_io.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>

#include <core/scoring/constraints/util.hh>

#include <utility/excn/Exceptions.hh>


#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <basic/Tracer.hh>
using basic::T;
using basic::Error;
using basic::Warning;

using namespace core;
using namespace core::scoring;
using namespace core::chemical;
using namespace core::scoring::etable;
using namespace core::scoring::etable::count_pair;
using namespace core::scoring::methods;

static basic::Tracer TR("fasol_refit");

OPT_1GRP_KEY(Real, fasol, del_ddg)
OPT_1GRP_KEY(Real, fasol, max_sasa)
OPT_1GRP_KEY(Boolean, fasol, revert)
OPT_1GRP_KEY(Boolean, fasol, classic)

class FaSolReporter : public protocols::moves::Mover {
	public:
	FaSolReporter(){}
		void apply( core::pose::Pose & pose) {
		using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::scoring::ScoreFunctionOP sf = core::scoring::get_score_function();
	(*sf)(pose);

	// get per-atom SASA
	core::Real probe_radius=1.4;
	core::id::AtomID_Map< Real > atom_sasa;
	utility::vector1< Real > rsd_sasa;
	core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probe_radius );

	EnergyGraph const & energy_graph( pose.energies().energy_graph() );
	methods::EnergyMethodOptions e_opts = sf->energy_method_options();
	Etable etable = *( ScoringManager::get_instance()->etable( e_opts ).lock() ); // copy

	core::scoring::methods::LK_BallEnergy lkb(e_opts);
	lkb.setup_for_scoring(pose, *sf);
	LKB_PoseInfo const & lkbposeinfo
		( static_cast< LKB_PoseInfo const & >( pose.data().get( pose::datacache::CacheableDataType::LK_BALL_POSE_INFO ) ) );

	// atomtype used for burial calcs
	int bur_type=1;
	while ( etable.lk_dgfree( bur_type ) == 0 || etable.lk_lambda( bur_type ) != 3.5 ) {
		bur_type++;
	}

	for ( Size ires=1; ires<= pose.total_residue(); ++ires ) {
		core::conformation::Residue const & rsd1( pose.residue(ires) );
		if ( !rsd1.is_protein() ) continue;
		if ( rsd_sasa[ires] != 0 ) continue;

		LKB_ResidueInfo const &lkbinfo1 = lkbposeinfo[ires];
		utility::vector1< utility::vector1< numeric::xyzVector<Real> > > const & rsd1_waters( lkbinfo1.waters() );
		utility::vector1< utility::vector1< Real > > const & rsd1_atom_wts( lkbinfo1.atom_weights() );

		core::Real fa_sol_i = 0.0, lk_ball_i = 0.0, fa_atr_i=0, fa_rep_i=0;
		core::Real burial_i = 0.0;

		for ( Size iatm=rsd1.first_sidechain_atom(); iatm<= rsd1.natoms(); ++iatm ) {
			if ( iatm > rsd1.nheavyatoms() && iatm < rsd1.first_sidechain_hydrogen() ) continue;

			core::conformation::Atom a_i = rsd1.atom(iatm);
			a_i.type(bur_type);

			utility::vector1< numeric::xyzVector< core::Real > > const & atom1_waters( rsd1_waters[ iatm ] );
			//numeric::xyzVector< core::Real > const & atom1_xyz( rsd1.xyz( iatm ) );
			utility::vector1< core::Real > const & atom1_wts( rsd1_atom_wts[iatm] );

			Real const sasa_this_atom( atom_sasa[ core::id::AtomID( iatm, ires ) ] );
			if ( sasa_this_atom > option[fasol::max_sasa] ) continue;

			for ( graph::Graph::EdgeListConstIter
					iru  = energy_graph.get_node(ires)->const_edge_list_begin(),
					irue = energy_graph.get_node(ires)->const_edge_list_end();
					iru != irue; ++iru ) {
				EnergyEdge const * edge( static_cast< EnergyEdge const *> (*iru) );
				int jres=(int)edge->get_other_ind(ires);
				core::conformation::Residue const &rsd2( pose.residue(jres) );
				if ( !rsd2.is_protein() ) continue;

				//LKB_ResidueInfo const &lkbinfo2 = lkbposeinfo[jres];
				//utility::vector1< utility::vector1< numeric::xyzVector<Real> > > const & rsd2_waters( lkbinfo2.waters() );
				//utility::vector1< utility::vector1< Real > > const & rsd2_atom_wts( lkbinfo2.atom_weights() );

				// count pair
				CountPairFunctionOP cpfxn =
					CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

				for ( Size jatm=1; jatm<= rsd2.natoms(); ++jatm ) {
					numeric::xyzVector< core::Real > const & atom2_xyz( rsd2.xyz( iatm ) );

					core::Real weight=sf->get_weight( core::scoring::fa_sol );
					core::Real weightA=sf->get_weight( core::scoring::fa_atr );
					core::Real weightR=sf->get_weight( core::scoring::fa_rep );
					core::Real weightLKb=sf->get_weight( core::scoring::lk_ball_wtd );

					core::Size path_dist;
					if ( cpfxn->count( iatm, jatm, weight, path_dist ) ) {
						core::Real fasol1,fasol2, ljatr, ljrep, lkjunk, dis2;

						etable.analytic_etable_evaluation( rsd1.atom(iatm), rsd2.atom(jatm), ljatr, ljrep, lkjunk, dis2);
						fa_atr_i += 0.5*(weightA*ljatr);
						fa_rep_i += 0.5*(weightR*ljrep);

						if ( iatm <= rsd1.nheavyatoms() && jatm <= rsd2.nheavyatoms() ) {
							etable.analytic_lk_energy(rsd1.atom(iatm), rsd2.atom(jatm), fasol1,fasol2 );
							fa_sol_i += weight*fasol1;

							if ( !atom1_waters.empty() ) {
								Real const fasol1_lkball =
									fasol1 * lkb.get_lk_fractional_contribution( atom2_xyz, rsd2.atom(jatm).type(), atom1_waters );
								lk_ball_i += weightLKb * ( atom1_wts[1] * fasol1 + atom1_wts[2] * fasol1_lkball );
							}

							etable.analytic_lk_energy(a_i, rsd1.atom(jatm), fasol1,fasol2 );
							burial_i += weight*fasol1; ///etable.lk_dgfree( bur_type );
						}
					}
				}
			}

			// add fa_intra_* contribution
			CountPairFunctionOP cpfxn =
				CountPairFactory::create_intrares_count_pair_function( rsd1, CP_CROSSOVER_4 );
			for ( Size jatm=1; jatm<= rsd1.natoms(); ++jatm ) {
				core::Real weight=sf->get_weight( core::scoring::fa_intra_sol_xover4 );
				core::Real weightA=sf->get_weight( core::scoring::fa_intra_atr_xover4 );
				core::Real weightR=sf->get_weight( core::scoring::fa_intra_rep_xover4 );
				core::Size path_dist;
				if ( cpfxn->count( iatm, jatm, weight, path_dist ) ) {
					core::Real fasol1,fasol2, ljatr, ljrep, lkjunk, dis2;

					etable.analytic_etable_evaluation( rsd1.atom(iatm), rsd1.atom(jatm), ljatr, ljrep, lkjunk, dis2);
					fa_atr_i += 0.5*(weightA*ljatr);
					fa_rep_i += 0.5*(weightR*ljrep);

					if ( iatm <= rsd1.nheavyatoms() && jatm <= rsd1.nheavyatoms() ) {
						etable.analytic_lk_energy(rsd1.atom(iatm), rsd1.atom(jatm), fasol1,fasol2 );
						fa_sol_i += weight*fasol1;
						etable.analytic_lk_energy(a_i, rsd1.atom(jatm), fasol1,fasol2 );
						burial_i += weight*fasol1; ///etable.lk_dgfree( bur_type );
					}
				}
			}
		}
		std::string base_name = protocols::jd2::JobDistributor::get_instance()->current_job()->input_tag();

		TR << base_name << " " << pose.pdb_info()->number(ires) << " " << pose.total_residue()
			<< " " << pose.residue(ires).name1() << " " << fa_sol_i+lk_ball_i << " " << fa_atr_i+fa_rep_i << " " << burial_i << std::endl;
	}

}
virtual std::string get_name() const {
	return "FaSolReporter";
}
};

///////
///////
void*
my_main( void* ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	try{
		protocols::jd2::JobDistributor::get_instance()->go( protocols::moves::MoverOP( new FaSolReporter() ) );
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
		NEW_OPT(fasol::del_ddg, "del_ddg value", -0.2);
		NEW_OPT(fasol::max_sasa, "max_sasa", 0.0);
		NEW_OPT(fasol::revert, "revert?", false);
		NEW_OPT(fasol::classic, "use classic parameters", false);

		devel::init(argc, argv);

		protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
