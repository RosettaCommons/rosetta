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
#include <core/kinematics/FoldTree.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
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
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/symmetry/SymPackRotamersMover.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/pose/util.hh>

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
#include <protocols/jd2/JobDistributorFactory.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/SilentFileJobOutputter.hh>

#include <core/io/pdb/pdb_writer.hh>
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
using basic::Error;
using basic::Warning;
using namespace core;
using namespace protocols;
using utility::vector1;



OPT_1GRP_KEY(Boolean, min, scoreonly)
OPT_1GRP_KEY(Boolean, min, debug)
OPT_1GRP_KEY(Boolean, min, debug_verbose)
OPT_1GRP_KEY(Boolean, min, cartesian)
OPT_1GRP_KEY(Boolean, min, pack)
OPT_1GRP_KEY(RealVector, min, ramp)
OPT_1GRP_KEY(String, min, minimizer)
OPT_1GRP_KEY(StringVector, min, fix_chains)

static basic::Tracer TR( "min_test" );

bool rama_list_pred( const std::pair < core::Size, core::Real > &left, const std::pair < core::Size, core::Real > &right ) {
	return left.second > right.second;
}


class MinTestMover : public protocols::moves::Mover {
public:
	MinTestMover()= default;

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

		core::Size nnonvrt_cst_target = constraint_pose.size();
		core::Size nnonvrt_pose = pose.size();

		while ( pose.residue( nnonvrt_pose ).aa() == core::chemical::aa_vrt ) { nnonvrt_pose--; }
		while ( constraint_pose.residue( nnonvrt_cst_target ).aa() == core::chemical::aa_vrt ) { nnonvrt_cst_target--; }

		if ( nnonvrt_pose != nnonvrt_cst_target ) {
			TR << "ERROR coord constraint pose length mismatch with input pose: " << nnonvrt_cst_target << " vs. " << nnonvrt_pose << std::endl;
			utility_exit();
		}

		if ( pose.residue( pose.fold_tree().root() ).aa() != core::chemical::aa_vrt ) {
			core::pose::addVirtualResAsRoot( pose );
		}

		Size nres = pose.size();
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
		for ( auto & g : rama_list ) {
			TR << "RAMALIST: " << g.first << "  " << g.second << std::endl;
		}

		kinematics::MoveMap my_mm;
		my_mm.set_bb(true);
		for ( Size ii=1; ii<= pose.size(); ii ++ ) {
			my_mm.set( TorsionID( omega_torsion, BB, ii), false );
			// disallow proline PHI
			if ( pose.residue(ii).aa() == chemical::aa_pro ) my_mm.set( TorsionID( phi_torsion, BB, ii), false );
		}

		Size nrounds = 5; // rama_list.size()
		for ( Size g=0; g< nrounds; g++ ) {
			// find bad ramas
			for ( Size j=1; j<= pose.size(); ++j ) {
				EnergyMap & emap( energies.onebody_energies( j ) );
				if (  emap[ rama ] > limit_rama_min ) {
					rama_list.emplace_back( j, emap[ rama ] );
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
				auto nsteps = (core::Size) std::ceil ( dist / stepsize );

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
		for ( Size j=1; j<= pose.size(); ++j ) {
			EnergyMap & emap( energies.onebody_energies( j ) );
			if (  emap[ rama ] > limit_rama_min ) {
				rama_list.emplace_back( j, emap[ rama ] );
			}
		}
		if ( rama_list.size() == 0 ) return;
		std::sort(rama_list.begin(), rama_list.end(), rama_list_pred);

		TR << "FINAL: " << std::endl;
		for ( auto & g : rama_list ) {
			TR << "RAMALIST: " << g.first << "  " << g.second << std::endl;
		}
	}

	void set_foldtree_for_variable_movement(core::pose::Pose & pose) {
		core::pose::addVirtualResAsRoot( pose );

		core::kinematics::FoldTree const &f_in = core::conformation::symmetry::get_asymm_unit_fold_tree( pose.conformation() );
		Size totres = f_in.nres();
		Size nres = totres - 1;
		core::Size vrtid = nres+1;

		utility::vector1< core::Size > cuts;
		utility::vector1< std::pair<core::Size,core::Size> > jumps;
		utility::vector1< Size > cuts_in = f_in.cutpoints();
		std::sort( cuts_in.begin(), cuts_in.end() );
		for ( core::Size i=1; i<=cuts_in.size(); ++i ) {
			core::Size seg_start = (i==1) ? 1 : cuts_in[i-1]+1;
			core::Size seg_end = cuts_in[i];
			core::Size jump_end = seg_start + (seg_end-seg_start)/2;
			cuts.push_back( cuts_in[i] );
			jumps.push_back( std::pair<core::Size,core::Size>( vrtid, jump_end ) );
		}

		ObjexxFCL::FArray2D< Size > fjumps( 2, jumps.size() );
		ObjexxFCL::FArray1D< Size > fcuts ( cuts.size() );
		for ( Size i=1; i<=jumps.size(); ++i ) {
			fjumps(1,i) = std::min( jumps[i].first , jumps[i].second );
			fjumps(2,i) = std::max( jumps[i].first , jumps[i].second );
		}
		for ( Size i = 1; i<=cuts.size(); ++i ) {
			fcuts(i) = cuts[i];
		}

		kinematics::FoldTree f;
		bool valid_tree = f.tree_from_jumps_and_cuts( nres+1, jumps.size(), fjumps, fcuts );
		runtime_assert( valid_tree );
		f.reorder( vrtid );

		TR << "New (asu) fold tree: " << f << std::endl;
		core::pose::symmetry::set_asymm_unit_fold_tree( pose , f );
	}

	void apply( core::pose::Pose & pose) override {
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
		mm.set_branches( true );
		mm.set_nu  ( true );
		mm.set_chi ( true );
		mm.set_jump( true );
		mm.set( core::id::THETA, option[ OptionKeys::relax::minimize_bond_angles ]() );
		mm.set( core::id::D, option[ OptionKeys::relax::minimize_bond_lengths ]() );

		utility::vector1< core::Real > ramp =  option[ OptionKeys::min::ramp ]();
		core::Size repeats = ramp.size();

		if ( option[ OptionKeys::relax::jump_move ].user() ) {
			mm.set_jump( option[ OptionKeys::relax::jump_move ]() );
		}
		if ( option[ OptionKeys::relax::bb_move ].user() ) {
			mm.set_bb( option[ OptionKeys::relax::bb_move ]() );
		}
		if ( option[ OptionKeys::relax::chi_move ].user() ) {
			mm.set_chi( option[ OptionKeys::relax::chi_move ]() );
		}

		if ( option[ OptionKeys::min::fix_chains ].user() ) {
			set_foldtree_for_variable_movement(pose);

			utility::vector1<std::string> chains_to_fix = option[ OptionKeys::min::fix_chains ]();
			for ( core::Size i=1; i<=chains_to_fix.size(); ++i ) {
				runtime_assert( chains_to_fix[i].length() == 1);
				for ( core::Size j=1; j<=pose.size(); ++j ) {
					if ( pose.pdb_info()->chain(j) == chains_to_fix[i][0] ) {
						mm.set_bb  ( j, false );
						mm.set_chi ( j, false );
					}
				}

				for ( core::Size j=1; j<=pose.num_jump(); ++j ) {
					id::AtomID atm_j = pose.conformation().jump_atom_id( j );
					if ( pose.pdb_info()->chain(atm_j.rsd()) == chains_to_fix[i][0] ) {
						mm.set_jump ( j, false );
					}
				}
			}
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

			if ( scorefxn->get_weight( core::scoring::atom_pair_constraint ) == 0 ) {
				scorefxn->set_weight( core::scoring::atom_pair_constraint, option[ OptionKeys::constraints::cst_fa_weight ] );
			}
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

		long t1=clock();
		TR << "start score: " << (*scorefxn)(pose) << std::endl;
		for ( int i=1; i<=(int)repeats; ++i ) {
			core::scoring::ScoreFunctionOP local_sf = scorefxn->clone();
			local_sf->set_weight( core::scoring::fa_rep,
				ramp[i]*scorefxn->get_weight( core::scoring::fa_rep )
			);

			// repack
			if ( option[ OptionKeys::min::pack ]() || option[ OptionKeys::relax::ramady ]() )  {
				TaskFactoryOP local_tf( new TaskFactory() );
				local_tf->push_back(TaskOperationCOP( new InitializeFromCommandline() ));
				local_tf->push_back(TaskOperationCOP( new RestrictToRepacking() ));
				local_tf->push_back(TaskOperationCOP( new IncludeCurrent() ));

				// mask by movemap
				bool repack = true;
				if ( basic::options::option[ basic::options::OptionKeys::relax::chi_move].user() ) {
					repack = basic::options::option[ basic::options::OptionKeys::relax::chi_move]();
				}

				PreventRepackingOP prevent_some( new PreventRepacking() );
				for ( Size i = 1; i<= pose.size() ; ++i ) {
					if ( !mm.get_chi(i) || !repack ) {
						prevent_some->include_residue(i);
					}
				}
				local_tf->push_back( prevent_some );

				protocols::minimization_packing::PackRotamersMoverOP pack_full_repack( new protocols::minimization_packing::PackRotamersMover( scorefxn ) );
				if ( core::pose::symmetry::is_symmetric( pose ) )  {
					pack_full_repack = protocols::minimization_packing::PackRotamersMoverOP( new minimization_packing::symmetry::SymPackRotamersMover( scorefxn ) );
				}
				pack_full_repack->task_factory(local_tf);
				pack_full_repack->apply( pose );

				(*scorefxn)(pose);
			}

			bool debug_verbose = option[ OptionKeys::min::debug_verbose ]();
			bool debug_derivs = option[ OptionKeys::min::debug ]() | debug_verbose;
			std::string minimizer_name = option[ OptionKeys::min::minimizer ]();

			// setup the options
			if ( !option[ OptionKeys::min::cartesian ]() )  {
				if ( option[ OptionKeys::symmetry::symmetry_definition ].user() )  {
					core::optimization::MinimizerOptions options( minimizer_name, 0.00001, true, debug_derivs, debug_verbose );
					core::optimization::symmetry::SymAtomTreeMinimizer minimizer;
					minimizer.run( pose, mm, *scorefxn, options );
				} else {
					core::optimization::MinimizerOptions options( minimizer_name, 0.00001, true, debug_derivs, debug_verbose );
					core::optimization::AtomTreeMinimizer minimizer;
					minimizer.run( pose, mm, *scorefxn, options );
				}
			} else {
				core::optimization::MinimizerOptions options( minimizer_name, 0.00001, true, debug_derivs, debug_verbose );
				core::optimization::CartesianMinimizer minimizer;
				minimizer.run( pose, mm, *scorefxn, options );
			}
		}
		scorefxn->show(TR, pose);
		long t2=clock();
		double time = ((double)t2 - t1) / CLOCKS_PER_SEC;
		TR << "end score: " << (*scorefxn)(pose) << std::endl;
		TR << "MIN TIME: " << time << " sec" << std::endl;

		core::scoring::dssp::Dssp my_dssp( pose );
		//core::scoring::dssp::DsspOP my_dssp_OP = new core::scoring::dssp::Dssp( pose );
		my_dssp.insert_ss_into_pose( pose );

		// add counts of ss elts to scorefile
		core::Size nL=0,nH=0,nE=0;
		for ( core::Size i=1; i<=pose.size(); ++i ) {
			if ( !pose.residue(i).is_protein() ) continue;
			if ( pose.secstruct(i) == 'L' ) nL++;
			if ( pose.secstruct(i) == 'H' ) nH++;
			if ( pose.secstruct(i) == 'E' ) nE++;
		}
		core::pose::setPoseExtraScore( pose, "nL", nL);
		core::pose::setPoseExtraScore( pose, "nH", nH);
		core::pose::setPoseExtraScore( pose, "nE", nE);
	}

	std::string get_name() const override {
		return "MinTestMover";
	}
};

void*
my_main( void* ) {
	using namespace protocols::moves;
	using namespace protocols::jd2;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	try{
		// Set up a job outputter that writes a scorefile calling evaluators
		if ( option[min::scoreonly]() ) {
			SilentFileJobOutputterOP jobout( new SilentFileJobOutputter );
			jobout->set_write_no_structures();
			jobout->set_write_separate_scorefile(true);
			protocols::jd2::JobDistributor::get_instance()->set_job_outputter( JobDistributorFactory::create_job_outputter( jobout ));
		}

		protocols::jd2::JobDistributor::get_instance()->go( protocols::moves::MoverOP( new MinTestMover() ) );
	} catch (utility::excn::Exception& excn ) {
		std::cerr << "Exception: " << std::endl;
		excn.show( std::cerr );
	}

	return nullptr;
}


///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {
		NEW_OPT(min::scoreonly, "scoresonly?", false);
		NEW_OPT(min::debug, "debug derivs?", false);
		NEW_OPT(min::debug_verbose, "debug derivs verbose?", false);
		NEW_OPT(min::cartesian, "cartesian minimization?", false);
		NEW_OPT(min::ramp, "ramp", utility::vector1<core::Real>(1,1.0));
		NEW_OPT(min::pack, "pack first?", false);
		NEW_OPT(min::minimizer, "minimizer?", "lbfgs_armijo_nonmonotone");
		NEW_OPT(min::fix_chains, "fix chains", utility::vector1<std::string>());

		devel::init(argc, argv);

		protocols::viewer::viewer_main( my_main );
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
