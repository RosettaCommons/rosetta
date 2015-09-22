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
#include <core/kinematics/FoldTree.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>

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


#include <core/io/pdb/pose_io.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>

#include <core/scoring/constraints/util.hh>

#include <utility/excn/Exceptions.hh>


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
using basic::T;
using basic::Error;
using basic::Warning;
using namespace core;
using namespace protocols;
using utility::vector1;
using io::pdb::dump_pdb;


OPT_1GRP_KEY(Boolean, min, debug)
OPT_1GRP_KEY(Boolean, min, debug_verbose)
OPT_1GRP_KEY(Boolean, min, cartesian)
OPT_1GRP_KEY(Boolean, min, pack)
OPT_1GRP_KEY(Integer, min, repeats)
OPT_1GRP_KEY(String, min, minimizer)

static THREAD_LOCAL basic::Tracer TR( "min_test" );

bool rama_list_pred( const std::pair < core::Size, core::Real > &left, const std::pair < core::Size, core::Real > &right ) {
	return left.second > right.second;
}


class MinTestMover : public protocols::moves::Mover {
public:
	MinTestMover(){}

	// TO DO make symm-friendly
	void add_coordinate_constraints_to_pose( core::pose::Pose & pose, const core::pose::Pose &constraint_pose ){
		using namespace core;
		using namespace conformation;
		using namespace pose;
		using namespace scoring;
		using namespace constraints;
		using namespace id;
		using namespace kinematics;
		using namespace moves;

		core::Size nnonvrt_cst_target = constraint_pose.total_residue();
		core::Size nnonvrt_pose = pose.total_residue();

		while ( pose.residue( nnonvrt_pose ).aa() == core::chemical::aa_vrt ) { nnonvrt_pose--; }
		while ( constraint_pose.residue( nnonvrt_cst_target ).aa() == core::chemical::aa_vrt ) { nnonvrt_cst_target--; }

		if ( nnonvrt_pose != nnonvrt_cst_target ) {
			TR << "ERROR coord constraint pose length mismatch with input pose: " << nnonvrt_cst_target << " vs. " << nnonvrt_pose << std::endl;
			utility_exit();
		}

		if ( pose.residue( pose.fold_tree().root() ).aa() != core::chemical::aa_vrt ) {
			core::pose::addVirtualResAsRoot( pose );
		}

		Size nres = pose.total_residue();
		Real const coord_sdev( 0.5 );
		for ( Size i = 1; i<= (Size)nres; ++i ) {
			if ( i==(Size)pose.fold_tree().root() ) continue;
			Residue const & nat_i_rsd( constraint_pose.residue(i) );
			for ( Size ii = 1; ii<= 3; ++ii ) {  // N/CA/C only
				func::FuncOP fx( new core::scoring::func::HarmonicFunc( 0.0, coord_sdev ) );
				pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new CoordinateConstraint(
					AtomID(ii,i), AtomID(1,nres), nat_i_rsd.xyz( ii ),
					fx ) ) ) );
			}
		}
	}

	// TO DO make symm-friendly
	void fix_worst_bad_ramas( core::pose::Pose & pose, core::Real limit_RMS=1.0, core::Real limit_rama=2.0 ){
		using namespace core;
		using namespace id;
		using namespace optimization;
		using namespace protocols::moves;
		using namespace core::scoring;
		using namespace core::pose;
		using namespace conformation;
		using namespace kinematics;

		const core::Real limit_rama_min = limit_rama;
		const core::Real stepsize = 360.0;

		// Original RAFT set
		//-140  153 180 0.135 B
		// -72  145 180 0.155 B
		//-122  117 180 0.073 B
		// -82  -14 180 0.122 A
		// -61  -41 180 0.497 A
		//  57   39 180 0.018 L

		core::Real ok_phi[] =  { -140, -72, -122, -82, -61, 57 } ;
		core::Real ok_psi[] =  {  153, 145,  117, -14, -41, 39 } ;
		core::pose::Pose original_pose = pose;
		core::pose::Pose temp_pose, best_pose;

		add_coordinate_constraints_to_pose( pose, original_pose ) ;

		core::scoring::ScoreFunction rama_scorefxn;
		rama_scorefxn.set_weight( coordinate_constraint, 0.2 );
		rama_scorefxn.set_weight( rama, 1.0 );
		rama_scorefxn.set_weight( omega, 0.2 );

		Energies & energies( pose.energies() );
		rama_scorefxn(pose); //apply score

		core::optimization::AtomTreeMinimizer minimizer;
		MinimizerOptions options( "lbfgs_armijo", 0.1, true /*use_nblist*/, false /*deriv_check*/ );
		options.max_iter(25);  //?
		kinematics::MoveMap final_mm;
		final_mm.set_bb(true);

		std::vector < std::pair < core::Size, core::Real > > rama_list;

		TR << "INITIAL: " << std::endl;
		for ( Size g=0; g<rama_list.size(); g++ ) {
			TR << "RAMALIST: " << rama_list[g].first << "  " << rama_list[g].second << std::endl;
		}

		kinematics::MoveMap my_mm;
		my_mm.set_bb(true);
		for ( Size ii=1; ii<= pose.total_residue(); ii ++ ) {
			my_mm.set( TorsionID( omega_torsion, BB, ii), false );
			// disallow proline PHI
			if ( pose.residue(ii).aa() == chemical::aa_pro ) my_mm.set( TorsionID( phi_torsion, BB, ii), false );
		}

		Size nrounds = 5; // rama_list.size()
		for ( Size g=0; g< nrounds; g++ ) {
			// find bad ramas
			for ( Size j=1; j<= pose.total_residue(); ++j ) {
				EnergyMap & emap( energies.onebody_energies( j ) );
				if (  emap[ rama ] > limit_rama_min ) {
					rama_list.push_back( std::make_pair( j, emap[ rama ] ) );
				}
			}
			if ( rama_list.size() == 0 ) return;
			std::sort(rama_list.begin(), rama_list.end(), rama_list_pred);


			core::Size i = rama_list[g].first;
			EnergyMap & emap( energies.onebody_energies( i ) );
			//core::Real rama_e =  emap[ rama ];

			core::Real curphi = pose.phi(i);
			core::Real curpsi = pose.psi(i);

			// save the original angles
			temp_pose = pose;

			core::Real bestscore = 1e6;
			for ( Size a=0; a<6; ++a ) {
				core::Real diffphi =  ok_phi[a]  - curphi; while ( diffphi > 180 ) diffphi-=360.0;  while ( diffphi < -180 ) diffphi += 360.0;
				core::Real diffpsi =  ok_psi[a]  - curpsi; while ( diffpsi > 180 ) diffpsi-=360.0;  while ( diffpsi < -180 ) diffpsi += 360.0;
				core::Real dist = sqrt( diffphi*diffphi + diffpsi * diffpsi );
				core::Size nsteps = (core::Size) std::ceil ( dist / stepsize );

				core::Real stepphi = diffphi / nsteps;
				core::Real steppsi = diffpsi / nsteps;

				for ( Size n=1; n<=nsteps; ++n ) {
					pose.set_phi( i, curphi + n*stepphi );
					pose.set_psi( i, curpsi + n*steppsi );
					minimizer.run( pose, my_mm, rama_scorefxn, options );
				}

				// check score
				core::Real tgt = rama_scorefxn(pose); //apply score
				core::Real rama_final =  emap[ rama ];
				core::Real rms_final = core::scoring::CA_rmsd( original_pose, pose );
				core::Real score_final = tgt + 10*rms_final;

				TR << "TRY " << g << "." << a << "  " << nsteps << "  " << score_final << "  " << rama_final << "  " << rms_final << std::endl;

				if ( score_final<bestscore && rms_final<limit_RMS ) {
					bestscore = score_final;
					best_pose = pose;
				}

				// restore
				pose = temp_pose;
			}
			TR << "BEST " << bestscore << std::endl;

			pose = best_pose;
		}

		// find bad ramas
		rama_list.clear();
		for ( Size j=1; j<= pose.total_residue(); ++j ) {
			EnergyMap & emap( energies.onebody_energies( j ) );
			if (  emap[ rama ] > limit_rama_min ) {
				rama_list.push_back( std::make_pair( j, emap[ rama ] ) );
			}
		}
		if ( rama_list.size() == 0 ) return;
		std::sort(rama_list.begin(), rama_list.end(), rama_list_pred);

		TR << "FINAL: " << std::endl;
		for ( Size g=0; g<rama_list.size(); g++ ) {
			TR << "RAMALIST: " << rama_list[g].first << "  " << rama_list[g].second << std::endl;
		}
	}

	void apply( core::pose::Pose & pose) {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::pack;
		using namespace core::pack::task;
		using namespace core::pack::task::operation;

		using namespace protocols::moves;
		using namespace scoring;

		// steal relax flags
		scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
		kinematics::MoveMap mm;
		mm.set_bb  ( true );
		mm.set_chi ( true );
		mm.set_jump( true );
		mm.set( core::id::THETA, option[ OptionKeys::relax::minimize_bond_angles ]() );
		mm.set( core::id::D, option[ OptionKeys::relax::minimize_bond_lengths ]() );

		core::Size repeats = option[ OptionKeys::min::repeats ]();

		if ( option[ OptionKeys::relax::jump_move ].user() ) {
			mm.set_jump( option[ OptionKeys::relax::jump_move ]() );
		}
		if ( option[ OptionKeys::relax::bb_move ].user() ) {
			mm.set_bb( option[ OptionKeys::relax::bb_move ]() );
		}
		if ( option[ OptionKeys::relax::chi_move ].user() ) {
			mm.set_chi( option[ OptionKeys::relax::chi_move ]() );
		}


		if ( option[ OptionKeys::symmetry::symmetry_definition ].user() )  {
			protocols::simple_moves::symmetry::SetupForSymmetryMoverOP symm( new protocols::simple_moves::symmetry::SetupForSymmetryMover );
			symm->apply( pose );
			core::pose::symmetry::make_symmetric_movemap( pose, mm );
		}

		// csts
		if ( option[ OptionKeys::constraints::cst_fa_file ].user() ) {
			protocols::simple_moves::ConstraintSetMoverOP loadCsts( new protocols::simple_moves::ConstraintSetMover );
			loadCsts->constraint_file( core::scoring::constraints::get_cst_fa_file_option() );
			loadCsts->apply(pose);
		}

		// now add density scores from cmd line
		if ( option[ edensity::mapfile ].user() ) {
			core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *scorefxn );
		}

		// set pose for density scoring if a map was input
		//   + (potentially) dock map into density
		if ( option[ edensity::mapfile ].user() ) {
			protocols::electron_density::SetupForDensityScoringMoverOP edens( new protocols::electron_density::SetupForDensityScoringMover );
			edens->apply( pose );
		}

		pose::PoseOP start_pose( new pose::Pose(pose) );
		(*scorefxn)(pose);
		scorefxn->show(TR, pose);

		// ramady
		if ( option[ OptionKeys::relax::ramady ]() )  {
			fix_worst_bad_ramas( pose );
		}

		// repack
		if ( option[ OptionKeys::min::pack ]() || option[ OptionKeys::relax::ramady ]() )  {
			TaskFactoryOP local_tf( new TaskFactory() );
			local_tf->push_back(TaskOperationCOP( new InitializeFromCommandline() ));
			local_tf->push_back(TaskOperationCOP( new RestrictToRepacking() ));
			local_tf->push_back(TaskOperationCOP( new IncludeCurrent() ));

			protocols::simple_moves::PackRotamersMoverOP pack_full_repack( new protocols::simple_moves::PackRotamersMover( scorefxn ) );
			if ( core::pose::symmetry::is_symmetric( pose ) )  {
				pack_full_repack = protocols::simple_moves::PackRotamersMoverOP( new simple_moves::symmetry::SymPackRotamersMover( scorefxn ) );
			}
			pack_full_repack->task_factory(local_tf);
			pack_full_repack->apply( pose );

			(*scorefxn)(pose);
			scorefxn->show(TR, pose);
		}


		bool debug_verbose = option[ OptionKeys::min::debug_verbose ]();
		bool debug_derivs = option[ OptionKeys::min::debug ]() | debug_verbose;
		std::string minimizer_name = option[ OptionKeys::min::minimizer ]();

		// setup the options
		if ( !option[ OptionKeys::min::cartesian ]() )  {
			if ( option[ OptionKeys::symmetry::symmetry_definition ].user() )  {
				core::optimization::MinimizerOptions options( minimizer_name, 0.00001, true, debug_derivs, debug_derivs );
				core::optimization::symmetry::SymAtomTreeMinimizer minimizer;
				TR << "SYMTORSION MINTEST: " << std::endl;
				TR << "start score: " << (*scorefxn)(pose)  << std::endl;
				long t1=clock();
				for ( int i=1; i<=(int)repeats; ++i ) minimizer.run( pose, mm, *scorefxn, options );
				long t2=clock();
				double time = ((double)t2 - t1) / CLOCKS_PER_SEC;
				TR << "end score: " << (*scorefxn)(pose)  << std::endl;
				TR << "MIN TIME: " << time << " sec " << std::endl;
			} else {
				core::optimization::MinimizerOptions options( minimizer_name, 0.00001, true, debug_derivs, debug_verbose );
				core::optimization::AtomTreeMinimizer minimizer;
				TR << "TORSION MINTEST: "  << std::endl;
				TR << "start score: " << (*scorefxn)(pose) << std::endl;
				long t1=clock();
				for ( int i=1; i<=(int)repeats; ++i ) minimizer.run( pose, mm, *scorefxn, options );
				long t2=clock();
				double time = ((double)t2 - t1) / CLOCKS_PER_SEC;
				TR << "end score: " << (*scorefxn)(pose) << std::endl;
				TR << "MIN TIME: " << time << " sec" << std::endl;
				(*scorefxn)(pose);
				scorefxn->show(TR, pose);
			}
		} else {
			core::optimization::MinimizerOptions options( minimizer_name, 0.00001, true, debug_derivs, debug_verbose );
			core::optimization::CartesianMinimizer minimizer;
			TR << "CART MINTEST: "  << std::endl;
			TR << "start score: " << (*scorefxn)(pose)  << std::endl;
			long t1=clock();
			for ( int i=1; i<=(int)repeats; ++i ) minimizer.run( pose, mm, *scorefxn, options );
			long t2=clock();
			double time = ((double)t2 - t1) / CLOCKS_PER_SEC;
			TR << "end score: " << (*scorefxn)(pose)  << std::endl;
			TR << "MIN TIME: " << time << " sec " << std::endl;
			(*scorefxn)(pose);
			scorefxn->show(TR, pose);
		}
	}
	virtual std::string get_name() const {
		return "MinTestMover";
	}
};

void*
my_main( void* ) {
	using namespace protocols::moves;
	using namespace protocols::simple_moves::symmetry;

	try{
		protocols::jd2::JobDistributor::get_instance()->go( protocols::moves::MoverOP( new MinTestMover() ) );
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
		NEW_OPT(min::debug, "debug derivs?", false);
		NEW_OPT(min::debug_verbose, "debug derivs verbose?", false);
		NEW_OPT(min::cartesian, "cartesian minimization?", false);
		NEW_OPT(min::repeats, "#repeats", 1);
		NEW_OPT(min::pack, "pack first?", false);
		NEW_OPT(min::minimizer, "minimizer?", "lbfgs_armijo_nonmonotone");

		devel::init(argc, argv);

		protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
