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
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
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

#include <protocols/membrane/AddMembraneMover.hh>

#include <core/scoring/electron_density/util.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>

#include <basic/basic.hh>
#include <basic/database/open.hh>
#include <devel/init.hh>
#include <protocols/symmetry/SetupForSymmetryMover.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MoverContainer.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobDistributorFactory.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/SilentFileJobOutputter.hh>

#include <core/io/pdb/pdb_writer.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <protocols/constraint_movers/ConstraintSetMover.hh>

#include <core/scoring/constraints/util.hh>

#include <utility/excn/Exceptions.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>

#include <devel/dna/relax_util.hh>
#include <core/scoring/dna/setup.hh>
#include <core/scoring/dna/BasePartner.hh>

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
OPT_1GRP_KEY(Boolean, min, dna)
OPT_1GRP_KEY(Boolean, min, dualspace)
OPT_1GRP_KEY(RealVector, min, ramp)
OPT_1GRP_KEY(RealVector, min, ramp_cst)
OPT_1GRP_KEY(BooleanVector, min, cart)
OPT_1GRP_KEY(String, min, minimizer)
OPT_1GRP_KEY(String, min, ref)
OPT_1GRP_KEY(StringVector, min, fix_chains)

static basic::Tracer TR( "min_test" );

bool rama_list_pred( const std::pair < core::Size, core::Real > &left, const std::pair < core::Size, core::Real > &right ) {
	return left.second > right.second;
}


class MinTestMover : public protocols::moves::Mover {
public:
	MinTestMover()= default;

	// TO DO make symm-friendly
	void add_coordinate_constraints_to_pose( core::pose::Pose & pose, const core::pose::Pose &constraint_pose, bool allatom=false ){
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

		Real const coord_sdev( 0.5 );
		for ( Size i = 1; i<=nnonvrt_pose; ++i ) {
			Residue const & nat_i_rsd( constraint_pose.residue(i) );
			if ( allatom ) {
				for ( Size ii = 1; ii<=nat_i_rsd.nheavyatoms(); ++ii ) {
					func::FuncOP fx( new core::scoring::func::HarmonicFunc( 0.0, coord_sdev ) );
					pose.add_constraint( scoring::constraints::ConstraintCOP( utility::pointer::make_shared< CoordinateConstraint >(
						AtomID(ii,i), AtomID(1,nnonvrt_pose+1), nat_i_rsd.xyz( ii ),
						fx ) ) );
				}
			} else {
				core::Size atm_index = 0;
				if ( nat_i_rsd.is_protein() ) {
					atm_index = nat_i_rsd.atom_index("CA");
				} else if ( nat_i_rsd.is_NA() ) {
					atm_index = nat_i_rsd.atom_index("P");
				}
				if ( atm_index != 0 ) {
					func::FuncOP fx( new core::scoring::func::HarmonicFunc( 0.0, coord_sdev ) );
					pose.add_constraint( scoring::constraints::ConstraintCOP( utility::pointer::make_shared< CoordinateConstraint >(
						AtomID(atm_index,i), AtomID(1,nnonvrt_pose+1), nat_i_rsd.xyz( atm_index ),
						fx ) ) );
				}
			}
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

		if (  option[ OptionKeys::in::membrane ].user() ) {
			using namespace protocols::membrane;
			AddMembraneMoverOP add_memb = AddMembraneMoverOP( new AddMembraneMover() );
			add_memb->apply( pose );
		}

		scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
		kinematics::MoveMap mm;
		mm.set_bb  ( true );
		mm.set_branches( true );
		mm.set_nu  ( true );
		mm.set_chi ( true );
		mm.set_jump( true );

		utility::vector1< core::Real > ramp =  option[ OptionKeys::min::ramp ]();
		core::Size repeats = ramp.size();

		utility::vector1< core::Real > ramp_cst =  option[ OptionKeys::min::ramp_cst ]();
		runtime_assert( ramp_cst.size() == 0 || ramp_cst.size() == repeats );

		utility::vector1< bool > cart =  option[ OptionKeys::min::cart ]();
		runtime_assert( cart.size() == 0 || cart.size() == repeats );


		// steal relax movemap flags
		mm.set( core::id::THETA, option[ OptionKeys::relax::minimize_bond_angles ]() );
		mm.set( core::id::D, option[ OptionKeys::relax::minimize_bond_lengths ]() );
		if ( option[ OptionKeys::relax::jump_move ].user() ) {
			mm.set_jump( option[ OptionKeys::relax::jump_move ]() );
		}
		if ( option[ OptionKeys::relax::bb_move ].user() ) {
			mm.set_bb( option[ OptionKeys::relax::bb_move ]() );
		}
		if ( option[ OptionKeys::relax::chi_move ].user() ) {
			mm.set_chi( option[ OptionKeys::relax::chi_move ]() );
		}

		// optionally freeze chains based on chain ID
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

		// optionally constrain to a reference pose
		if ( option[ OptionKeys::min::ref ].user() ) {
			core::pose::Pose reference_pose;
			core::import_pose::pose_from_file(reference_pose, option[ OptionKeys::min::ref ](), core::import_pose::PDB_file);
			add_coordinate_constraints_to_pose( pose, reference_pose, true );

			if ( scorefxn->get_weight( core::scoring::coordinate_constraint ) == 0 && ramp_cst.size() == 0 ) {
				scorefxn->set_weight( core::scoring::coordinate_constraint, option[ OptionKeys::constraints::cst_fa_weight ] );
			}
		}

		// symmetry
		if ( option[ OptionKeys::symmetry::symmetry_definition ].user() )  {
			protocols::symmetry::SetupForSymmetryMoverOP symm( new protocols::symmetry::SetupForSymmetryMover );
			symm->apply( pose );
			core::pose::symmetry::make_symmetric_movemap( pose, mm );
		}

		// constraints
		if ( option[ OptionKeys::constraints::cst_fa_file ].user() ) {
			protocols::constraint_movers::ConstraintSetMoverOP loadCsts( new protocols::constraint_movers::ConstraintSetMover );
			loadCsts->constraint_file( core::scoring::constraints::get_cst_fa_file_option() );
			loadCsts->apply(pose);

			if ( scorefxn->get_weight( core::scoring::atom_pair_constraint ) == 0 ) {
				scorefxn->set_weight( core::scoring::atom_pair_constraint, option[ OptionKeys::constraints::cst_fa_weight ] );
				scorefxn->set_weight( core::scoring::angle_constraint, option[ OptionKeys::constraints::cst_fa_weight ] );
				scorefxn->set_weight( core::scoring::dihedral_constraint, option[ OptionKeys::constraints::cst_fa_weight ] );
				scorefxn->set_weight( core::scoring::coordinate_constraint, option[ OptionKeys::constraints::cst_fa_weight ] );
			}
		}

		// add density scores from cmd line
		if ( option[ edensity::mapfile ].user() ) {
			core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *scorefxn );
		}

		// set pose for density scoring if a map was input
		if ( option[ edensity::mapfile ].user() ) {
			protocols::electron_density::SetupForDensityScoringMoverOP edens( new protocols::electron_density::SetupForDensityScoringMover );
			edens->apply( pose );
		}

		// start!  score pose and report pre-min energies
		pose::PoseOP start_pose( new pose::Pose(pose) );
		(*scorefxn)(pose);
		scorefxn->show(TR, pose);
		TR << "start score: " << (*scorefxn)(pose) << std::endl;

		// set foldtree for DNA
		if (  option[ min::dna ]() ) {
			kinematics::FoldTree ft (pose.size());
			core::scoring::dna::set_base_partner( pose );
			devel::dna::add_dna_base_jumps_to_fold_tree( pose, ft, false );
			devel::dna::set_dna_jump_atoms_in_fold_tree( pose, false, false, ft );
			pose.fold_tree( ft );
			TR << ft << std::endl;
		}

		// optimization loop
		long t1=clock();
		for ( core::Size i=1; i<=repeats; ++i ) {
			core::scoring::ScoreFunctionOP local_sf = scorefxn->clone();
			local_sf->set_weight( core::scoring::fa_rep,
				ramp[i]*scorefxn->get_weight( core::scoring::fa_rep )
			);
			if ( ramp_cst.size() != 0 ) {
				local_sf->set_weight( core::scoring::atom_pair_constraint, ramp_cst[i]);
				local_sf->set_weight( core::scoring::angle_constraint, ramp_cst[i]);
				local_sf->set_weight( core::scoring::dihedral_constraint, ramp_cst[i]);
				local_sf->set_weight( core::scoring::coordinate_constraint, ramp_cst[i]);
			}

			// repack
			if ( option[ OptionKeys::min::pack ]() )  {
				TaskFactoryOP local_tf( new TaskFactory() );
				local_tf->push_back(utility::pointer::make_shared< InitializeFromCommandline >());
				local_tf->push_back(utility::pointer::make_shared< RestrictToRepacking >());
				local_tf->push_back(utility::pointer::make_shared< IncludeCurrent >());

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

				protocols::minimization_packing::PackRotamersMoverOP pack_full_repack( new protocols::minimization_packing::PackRotamersMover( local_sf ) );
				if ( core::pose::symmetry::is_symmetric( pose ) )  {
					pack_full_repack = utility::pointer::make_shared< minimization_packing::symmetry::SymPackRotamersMover >( local_sf );
				}
				pack_full_repack->task_factory(local_tf);
				pack_full_repack->apply( pose );

				(*scorefxn)(pose);
			}

			bool debug_verbose = option[ OptionKeys::min::debug_verbose ]();
			bool debug_derivs = option[ OptionKeys::min::debug ]() | debug_verbose;
			std::string minimizer_name = option[ OptionKeys::min::minimizer ]();

			// setup the options
			bool cart_this_cycle = option[ OptionKeys::min::cartesian ]();
			if ( cart.size() > 0 ) {
				cart_this_cycle = cart[i];
			}

			if ( !cart_this_cycle )  {
				if ( option[ OptionKeys::symmetry::symmetry_definition ].user() )  {
					core::optimization::MinimizerOptions options( minimizer_name, 0.00001, true, debug_derivs, debug_verbose );
					core::optimization::symmetry::SymAtomTreeMinimizer minimizer;
					minimizer.run( pose, mm, *local_sf, options );
				} else {
					core::optimization::MinimizerOptions options( minimizer_name, 0.00001, true, debug_derivs, debug_verbose );
					core::optimization::AtomTreeMinimizer minimizer;
					minimizer.run( pose, mm, *local_sf, options );
				}
			} else {
				core::optimization::MinimizerOptions options( minimizer_name, 0.00001, true, debug_derivs, debug_verbose );
				core::optimization::CartesianMinimizer minimizer;
				minimizer.run( pose, mm, *local_sf, options );
			}

			if ( i != repeats ) {
				local_sf->show(TR, pose);
			}
		}

		// done!  report final energies
		long t2=clock();
		scorefxn->show(TR, pose);
		double time = ((double)t2 - t1) / CLOCKS_PER_SEC;
		TR << "end score: " << (*scorefxn)(pose) << std::endl;
		TR << "MIN TIME: " << time << " sec" << std::endl;

		//core::pose::setPoseExtraScore( pose, "nL", nL);
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
		protocols::jd2::JobDistributor::get_instance()->go( utility::pointer::make_shared< MinTestMover >() );
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
		NEW_OPT(min::ramp_cst, "ramp_cst", utility::vector1<core::Real>(0));
		NEW_OPT(min::cart, "cart at cycles", utility::vector1<bool>(0));
		NEW_OPT(min::dna, "dna mode?", false);
		NEW_OPT(min::dualspace, "dualspace mode?", false);
		NEW_OPT(min::pack, "pack first?", false);
		NEW_OPT(min::minimizer, "minimizer?", "lbfgs_armijo_nonmonotone");
		NEW_OPT(min::ref, "reference structure", "");
		NEW_OPT(min::fix_chains, "fix chains", utility::vector1<std::string>());

		devel::init(argc, argv);

		protocols::viewer::viewer_main( my_main );
	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}
