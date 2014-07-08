/// @file
/// @brief


#include <devel/init.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AA.hh>

#include <core/optimization/Multifunc.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/NelderMeadSimplex.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/elec/FA_ElecEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/util/cryst_util.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/pack/rtmin.hh>
#include <core/graph/Graph.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <utility/vector1.hh>
#include <numeric/random/random.hh>

#include <boost/thread.hpp>
#include <boost/bind.hpp>

using namespace core;
using namespace conformation;
using namespace chemical;
using namespace pack::rotamer_set;
using namespace pack::scmin;
using namespace pose;
using namespace scoring;
using namespace scoring::methods;
using namespace optimization;
using namespace graph;
using namespace basic::options;


OPT_1GRP_KEY(Integer, rropt, j)
OPT_1GRP_KEY(Integer, rropt, ncyc)
OPT_1GRP_KEY(Integer, rropt, maxiter)
OPT_1GRP_KEY(Real, rropt, H)
OPT_1GRP_KEY(Real, rropt, sigma)
OPT_1GRP_KEY(Boolean, rropt, allwts)
OPT_1GRP_KEY(Boolean, rropt, bbmove)
OPT_1GRP_KEY(Boolean, rropt, ideal)
OPT_1GRP_KEY(Boolean, rropt, etable)
OPT_1GRP_KEY(Boolean, rropt, hackelec)
OPT_1GRP_KEY(String, rropt, sflo)
OPT_1GRP_KEY(String, rropt, sfhi)


utility::vector1<core::pose::PoseOP>
load_models () {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	utility::vector1<core::pose::PoseOP> models = core::import_pose::poseOPs_from_pdbs( option[OptionKeys::in::file::s]());
	return models;
}

utility::vector1< utility::vector1<bool> >
select_residues_for_evaluation(utility::vector1<core::pose::PoseOP> const &models) {
	utility::vector1< utility::vector1<bool> > retval( models.size() );

	for (int i=1; i<=models.size(); ++i) {
		retval[i].resize( models[i]->total_residue() );
		core::Size accepted_residues = 0;

		for (int j=1; j<=models[i]->total_residue(); ++j) {
			core::Real Bsum=0.0, Bcount=0.0;
			core::conformation::Residue const &rsd_i = models[i]->residue(j);

			for (Size k=rsd_i.first_sidechain_atom(); k<=rsd_i.nheavyatoms(); ++k) {
				if (rsd_i.is_virtual(k)) continue;
				Real B = models[i]->pdb_info()->temperature( j, k );
				Bsum += B;
				Bcount += 1.0;
			}
			Bsum /= Bcount;

			retval[i][j] = (Bsum<=30 && rsd_i.aa() != core::chemical::aa_gly && rsd_i.aa() != core::chemical::aa_ala);
			if (retval[i][j]) accepted_residues++;
		}
		std::cerr << "model " << i << " repacking " << accepted_residues << " / " << models[i]->total_residue() << std::endl;
	}
	return retval;
}


void
sf2multivec( core::scoring::ScoreFunctionCOP scorefxn, Multivec &y, bool islow=false, bool ishigh=false) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if (!option[rropt::ideal]()) {
		y.push_back( scorefxn->get_weight( core::scoring::cart_bonded_length ) );
		y.push_back( scorefxn->get_weight( core::scoring::cart_bonded_angle ) );
		y.push_back( scorefxn->get_weight( core::scoring::cart_bonded_torsion ) );
	}
	if (option[rropt::allwts]()) {
		y.push_back( scorefxn->get_weight( core::scoring::fa_dun ) );
		y.push_back( scorefxn->get_weight( core::scoring::fa_atr ) );
		y.push_back( scorefxn->get_weight( core::scoring::fa_rep ) );
		y.push_back( scorefxn->get_weight( core::scoring::fa_sol ) );
		y.push_back( scorefxn->get_weight( core::scoring::fa_intra_rep ) );
		y.push_back( scorefxn->get_weight( core::scoring::fa_intra_sol ) );
		y.push_back( scorefxn->get_weight( core::scoring::fa_intra_atr ) );
		y.push_back( scorefxn->get_weight( core::scoring::fa_elec ) );
		y.push_back( scorefxn->get_weight( core::scoring::hbond_sr_bb ) );  // other hbond 'riding'
		if (option[rropt::bbmove]()) {
			y.push_back( scorefxn->get_weight( core::scoring::rama ) );
			y.push_back( scorefxn->get_weight( core::scoring::omega ) );
		}
	}
}

void
multivec2sf( Multivec const &vars, core::scoring::ScoreFunctionOP scorefxn ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	core::Size counter=1;

	if (!option[rropt::ideal]()) {
		scorefxn->set_weight( core::scoring::cart_bonded_length       , vars[counter++]) ;
		scorefxn->set_weight( core::scoring::cart_bonded_angle        , vars[counter++]) ;
		scorefxn->set_weight( core::scoring::cart_bonded_torsion      , vars[counter++]) ;
	}

	if (option[rropt::allwts]()) {
		scorefxn->set_weight( core::scoring::fa_dun , vars[counter++]) ;
		scorefxn->set_weight( core::scoring::fa_atr , vars[counter++]) ;
		scorefxn->set_weight( core::scoring::fa_rep , vars[counter++]) ;
		scorefxn->set_weight( core::scoring::fa_sol , vars[counter++]) ;
		scorefxn->set_weight( core::scoring::fa_intra_rep , vars[counter++]) ;
		scorefxn->set_weight( core::scoring::fa_intra_sol , vars[counter++]) ;
		scorefxn->set_weight( core::scoring::fa_intra_atr , vars[counter++]) ;
		scorefxn->set_weight( core::scoring::fa_elec , vars[counter++]) ;
		scorefxn->set_weight( core::scoring::hbond_sr_bb , vars[counter]) ;
		scorefxn->set_weight( core::scoring::hbond_lr_bb , vars[counter]) ;
		scorefxn->set_weight( core::scoring::hbond_bb_sc , vars[counter]) ;
		scorefxn->set_weight( core::scoring::hbond_sc , vars[counter++]*1.1/1.17) ;
		if (option[rropt::bbmove]()) {
			scorefxn->set_weight( core::scoring::rama , vars[counter++]) ;
			scorefxn->set_weight( core::scoring::omega, vars[counter++]) ;
		}
	}
}


void
dump_multivec( Multivec const &vars ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::cerr << "V[" << vars[1];
	for (int i=2; i<=vars.size(); ++i) {
		std::cerr << "," << vars[i];
	}
	std::cerr << "]";
}

/////////////////

void
do_rtmin(
		core::scoring::ScoreFunction const &sf,
		core::pack::task::PackerTask const &packer_task,
		core::pose::Pose const &pose, core::Real &rr_score,
		core::Real &abs_rr_score ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using core::chemical::aa_unk;

	core::Size npacked=0, nrotrec=0;

	core::pack::task::PackerTaskOP one_res_task( packer_task.clone() );
	core::graph::GraphOP packer_neighbor_graph = core::pack::create_packer_graph( pose, sf, one_res_task );

	core::pack::RTMin rtmin;
	rtmin.set_cartesian(!option[rropt::ideal]());
	rtmin.set_nonideal(option[rropt::bbmove]());

	one_res_task->temporarily_fix_everything();
	core::pose::Pose working_pose = pose;

	rr_score=0;
	abs_rr_score=0;
	for( core::Size ii = 1; ii <= pose.total_residue(); ++ii ){
		if ( !packer_task.pack_residue(ii) ) continue;
		if ( !pose.residue(ii).is_polymer() ) continue;

		one_res_task->temporarily_set_pack_residue( ii, true );
	    if ( ! packer_task.include_current( ii ) ) {
			core::pack::rotamer_set::RotamerSetFactory rsf;
			core::pack::rotamer_set::RotamerSetOP rotset( rsf.create_rotamer_set( working_pose.residue( ii ) ));
			rotset->set_resid( ii );
			rotset->build_rotamers( working_pose, sf, *one_res_task, packer_neighbor_graph );
			if ( rotset->num_rotamers() > 0 ) {
				working_pose.replace_residue( ii, *rotset->rotamer(1), false );
			}
		} else {
			utility_exit_with_message( "Including current!" );
		}

		rtmin.rtmin( working_pose, sf, one_res_task );
/*
		utility::vector1 <core::Real> &scores = rtmin.get_scores(1);
		utility::vector1 <utility::vector1 <core::Real> > &rots = rtmin.get_rotamers(1);

		core::Size nrots = scores.size();

		if (nrots == 0) continue;  // needed?

		core::Real minscore=1e30;
		for (core::Size jj=1; jj<=nrots; ++jj) {
			minscore = std::min( scores[jj], minscore );
		}

		// find closest rotamer
		core::Real prob_all=0.0, prob_right=0.0;
		core::Real best_right=1e30, best_wrong=1e30;
		core::Real best_all=1e30;
		core::Real best_chidiff=0;
		core::conformation::Residue const &res_ref = pose.residue(ii);
		for (core::Size jj=1; jj<=nrots; ++jj) {
			utility::vector1 <core::Real> const &res_rt = rots[jj];
			core::Real max_chidiff_i=0.0;
			for ( core::Size chi_index=1; chi_index <= res_ref.nchi(); ++chi_index ) {
				if ( res_ref.type().chi_2_proton_chi( chi_index ) != 0 ) continue;
				core::Real chidiff = std::abs(
						core::pack::dunbrack::subtract_chi_angles( res_ref.chi(chi_index), res_rt[chi_index], res_ref.aa(), chi_index )
					);
				max_chidiff_i = std::max( max_chidiff_i, chidiff );
			}

			// find low-energy rotamer _excluding_ closest
			//    anything within 20 deg is "right"
			if (max_chidiff_i<20) {
				prob_right += exp( -scores[jj]+minscore );
				best_right = std::min( scores[jj], best_right );
			} else {
				best_wrong = std::min( scores[jj], best_wrong );
			}
			prob_all += exp( -scores[jj]+minscore );

			if (scores[jj] < best_all) {
				best_all = scores[jj];
				best_chidiff = max_chidiff_i;
			}
		}

		//rr_score += prob_right / prob_all;
		core::Real k = (best_wrong - best_right);
		rr_score += 1.0/(1.0+exp(-option[rropt::sigma]()*k));
		if (best_right<best_wrong) abs_rr_score+= 1.0;

		working_pose.replace_residue( ii, pose.residue( ii ), false ); // restore the original conformation.
		one_res_task->temporarily_set_pack_residue( ii, false );
		*/
	}
}




//
class RTminOptMultifunc : public core::optimization::Multifunc {
public:
	RTminOptMultifunc(
		utility::vector1< core::pose::PoseOP > const &poses_in,
		utility::vector1<std::string> const &names_in,
		utility::vector1< utility::vector1<bool> > const &masks_in,
		core::scoring::ScoreFunction const & scorefxn_in
	) {
		poses_ = poses_in;
		names_ = names_in;
		masks_ = masks_in;
		scorefxn_ = new core::scoring::ScoreFunction( scorefxn_in );
		task_factory_ = new core::pack::task::TaskFactory;
		//task_factory_->push_back( new core::pack::task::operation::InitializeFromCommandline );
		task_factory_->push_back( new core::pack::task::operation::RestrictToRepacking );
	}

	virtual ~RTminOptMultifunc() {}

	// function
	virtual core::Real
	operator ()( core::optimization::Multivec const & vars ) const {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		multivec2sf( vars, scorefxn_);
		core::Real retval=0;


		utility::vector1< core::Real > thread_scores( poses_.size(), 0.0 );
		utility::vector1< core::Real > abs_thread_scores( poses_.size(), 0.0 );
		utility::vector1< core::pack::task::PackerTaskOP > packer_tasks( poses_.size() );

		// just to be safe give each thread their own copy
		utility::vector1< core::scoring::ScoreFunctionOP > scorefunctions( poses_.size() );

#ifdef BOOST_THREAD
		boost::thread_group threads;
		for (core::Size i=1; i<= poses_.size(); ++i) {
			packer_tasks[i] = task_factory_->create_task_and_apply_taskoperations( *poses_[i] );
			packer_tasks[i]->restrict_to_residues( masks_[i] );

			//core::scoring::ScoreFunction const &sf = *scorefxn_;
			scorefunctions[i] = scorefxn_->clone();
			core::scoring::ScoreFunction const &sf = *scorefunctions[i];
			(sf)(*poses_[i]);

			threads.create_thread(
				boost::bind(
					&do_rtmin, boost::cref(*scorefxn_), boost::cref(*packer_tasks[i]), boost::cref(*poses_[i]), boost::ref(thread_scores[i]) , boost::ref(abs_thread_scores[i])
				) );
		}
		threads.join_all();
#else
		std::cerr << "compile with extras=boost_thread!" << std::endl;
		exit(1);
#endif

		core::Real absretval=0;
		for (core::Size i=1; i<= poses_.size(); ++i) {
			retval += thread_scores[i];
			absretval += abs_thread_scores[i];
		}

		dump_multivec( vars );
		std::cerr << " = " << retval << "  " << absretval << std::endl;
		return -retval;
	}

	// derivatives
	virtual void
	dfunc( core::optimization::Multivec const & vars, core::optimization::Multivec & dE_dvars ) const {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		core::Real H=option[rropt::H]();
		core::Real score = (*this)(vars);

		for (int i=1; i<=vars.size(); ++i) {
			utility::vector1< core::Real > varsp = vars; varsp[ i ] += H;
			core::Real offset_score = (*this)(varsp);
			dE_dvars[i] = ( offset_score - score ) / (H);
		}
	}

private:
	utility::vector1< core::pose::PoseOP > poses_;
	utility::vector1< std::string > names_;
	utility::vector1< utility::vector1<bool> > masks_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::pack::task::TaskFactoryOP task_factory_;
}; // GradientOptMultifunc



int main( int argc, char * argv [] ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	NEW_OPT(rropt::j, "j", 1);
	NEW_OPT(rropt::ncyc, "ncyc", 50);
	NEW_OPT(rropt::maxiter, "maxiter", 1);
	NEW_OPT(rropt::H, "H", 0.01);
	NEW_OPT(rropt::sigma, "sigma", 6.0);
	NEW_OPT(rropt::allwts, "allwts?", false);
	NEW_OPT(rropt::bbmove, "bbmove?", false);
	NEW_OPT(rropt::ideal, "ideal?", false);
	NEW_OPT(rropt::sfhi, "sfhi", "dummy");

	devel::init(argc, argv);

	runtime_assert( !( option[rropt::bbmove]() && option[rropt::ideal]() ) );

	utility::vector1<core::pose::PoseOP> models = load_models();
	utility::vector1< utility::vector1<bool> > masks = select_residues_for_evaluation(models);
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
	core::scoring::ScoreFunctionOP start_scorefxn = scorefxn->clone();

	core::scoring::ScoreFunctionOP scorefxn_hi = core::scoring::ScoreFunctionFactory::create_score_function(option[rropt::sfhi]());

	// starting point
	core::optimization::Multivec y, yHI;
	sf2multivec( scorefxn, y );
	sf2multivec( scorefxn_hi, yHI, false, true );

	core::optimization::MinimizerOptions options( "lbfgs_armijo_nonmonotone", 1e-12, true, false, false );
	options.max_iter(option[rropt::maxiter]());

	for (int i=1; i<=option[rropt::ncyc](); ++i) {
		utility::vector1<core::pose::PoseOP> models_i;
		utility::vector1<std::string> names_i;
		utility::vector1< utility::vector1<bool> > masks_i;

		if (option[rropt::j]() < models.size()) {
			for (int j=0; j< option[rropt::j](); ++j) {
				core::Size modelnum = numeric::random::random_range(1, models.size());
				core::pose::PoseOP pose_i = new core::pose::Pose( *models[ modelnum ] );
				models_i.push_back( pose_i );
				names_i.push_back( option[OptionKeys::in::file::s]()[ modelnum ] );
				masks_i.push_back( masks[ modelnum ] );
			}
		} else {
			models_i = models;
			names_i = option[OptionKeys::in::file::s]();
			masks_i = masks;
		}
		RTminOptMultifunc f( models_i, names_i , masks_i,  *start_scorefxn );

		//core::optimization::Minimizer minimizer( f, options );
		//minimizer.run( y );

		core::optimization::NelderMeadSimplex minimizer( f, options );
		minimizer.run( y, yHI );
	}
}


