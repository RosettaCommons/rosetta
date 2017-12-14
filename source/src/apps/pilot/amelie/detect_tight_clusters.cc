// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Noah Ollikainen (original version, ca. 2008)
/// @author Amelie Stein (amelie.stein@ucsf.edu) -- edit/update/extra filters and features, August 2012

// Protocols Headers
//#include <protocols/rigid/RigidBodyMover.hh>

// Core Headers
// AS -- it's quite possible that a lot of the includes aren't needed...
#include <devel/init.hh>
//#include <core/io/pdb/pdb_writer.hh>
#include <core/import_pose/import_pose.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
//#include <core/id/AtomID_Map.Pose.hh>
//#include <core/util/MetricValue.hh>
#include <core/id/types.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
// #include <core/optimization/AtomTreeMinimizer.hh> // needed?
#include <protocols/simple_moves/MinMover.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>
#include <core/types.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <protocols/relax/FastRelax.hh>
//#include <protocols/simple_moves/MinPackMover.hh> // for MinPacker
#include <core/pack/min_pack.hh> // core/unwrapped MinPacker
#include <core/pack/rotamer_trials.hh>
#include <core/pack/rtmin.hh> // rotamer trials with minimization
//#include <core/pose/metrics/CalculatorFactory.hh>
//#include <core/pack/dunbrack/RotamerLibrary.hh> // AS -- for symmetry-respecting subtract_chi_angles function
#include <protocols/rotamer_recovery/RRComparer.hh>
#include <protocols/backrub/BackrubMover.hh>
//#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <protocols/moves/MonteCarlo.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/excn/Exceptions.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

#include <ObjexxFCL/format.hh>

using namespace core;

// Platform Headers

// option key includes
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>


#include <fstream>
#include <stack>
#include <map>

// constants / parameters that probably won't change much:
//core::Size min_num_neighbors = 9;
const core::Size buried_threshold = 14; // >14 neighbors <=> buried cluster
const core::Size neighbor_dist = 8;
//const core::Size b_factor_threshold = 30;
//core::Real heavy_atom_dist = 4.5;
// const core::Size chi_deviation_threshold = 40; // for now we're only reporting really bad repacking, deviations of 40+ in chi1 or chi2 -- might be interesting to control this via a flag though
const core::Real repack_sphere_radius = 8;
const core::Size rot_trials_iterations = 10;
const core::Size env_quality_check_dist = 6; // no zero-occupancy or Rosetta-rebuilt atoms within this distance of any cluster atom
const core::Size backrub_iterations = 10;

static basic::Tracer TR( "apps.pilot.amelie.detect_tight_clusters" );

OPT_1GRP_KEY(Boolean, detect_tight_clusters, generate_output_structures) // (no semicolon here!) -- this would generate huge amounts of output, so by default it runs silently, only reporting which residues are involved in clusters
OPT_1GRP_KEY(Boolean, detect_tight_clusters, require_renumbered_structures) // requires numbering to match length of structures -- hack for making the output work with external tools like OSCARstar, but has some drawbacks (see below)
OPT_1GRP_KEY(Boolean, detect_tight_clusters, debug)
OPT_1GRP_KEY(Boolean, detect_tight_clusters, repack_8A_sphere) // if not, just repack the 4 residues in the cluster
//OPT_1GRP_KEY(Integer, detect_tight_clusters, chi_deviation_threshold)
OPT_1GRP_KEY( IntegerVector, detect_tight_clusters, chi_deviation_thresholds ) // keep list of chi_dev_thresholds and iterate over those

OPT_1GRP_KEY(Boolean, detect_tight_clusters, min_pack)
//OPT_1GRP_KEY(Boolean, detect_tight_clusters, stochastic_pack)
OPT_1GRP_KEY(Boolean, detect_tight_clusters, rot_trials_min)
OPT_1GRP_KEY(Boolean, detect_tight_clusters, rot_trials)
OPT_1GRP_KEY(Boolean, detect_tight_clusters, sidechain_fastrelax)
OPT_1GRP_KEY(Boolean, detect_tight_clusters, cartmin)
OPT_1GRP_KEY(Boolean, detect_tight_clusters, soft_repack_min_hard_repack_min)
OPT_1GRP_KEY(Boolean, detect_tight_clusters, backrub)
OPT_1GRP_KEY(Real, detect_tight_clusters, backrub_sc_prob) // probability of making a side chain move within the backrub option -- default 0.25

OPT_1GRP_KEY(Boolean, detect_tight_clusters, detect_interchain_clusters) // standard detection works within proteins only, i.e. all members of the cluster have to be in the same chain -- if this flag is switched on, at least one residue must be in another chain -- may require adaptation of the burial threshold etc., in particular min_num_neighbors
OPT_1GRP_KEY(Real, detect_tight_clusters, min_num_neighbors) // within a chain we only want buried or intermediate clusters (9+) -- in interfaces this is not realistic though, may be adjusted

OPT_1GRP_KEY(Real, detect_tight_clusters, heavy_atom_dist) // 4.5A is the value Noah used -- can be increased slightly if we get too few clusters, in particular for interchain clusters
OPT_1GRP_KEY(Real, detect_tight_clusters, b_factor_threshold)  // 30 worked well for intrachain clusters

/*
the following set of options os for providing an external cluster (e.e., from a baseline run)
that can either be repacked and evaluated, or the user provides an externally repacked structure,
which needs to come with a reference structure via -in:file:native

a native structure can also be provided e.g. for repacking after backbone relaxation
-- in this case we may first want to check whether chi1/2 of the cluster residues are
still the same as in the native structure, but this isn't implemented yet
*/

OPT_1GRP_KEY(StringVector, detect_tight_clusters, external_cluster) // PDB numbering -- positions will be casted to int, probably won't be able to deal well with insertion codes etc.... but when starting the initial (cluster detection) run from a renumbered structure this should be fine
OPT_1GRP_KEY(Boolean, detect_tight_clusters, evaluate_externally_repacked_structure)
OPT_1GRP_KEY(Boolean, detect_tight_clusters, repack_external_cluster)

bool sort_min(const std::pair<core::Size, core::Real> & left, const std::pair<core::Size, core::Real> & right) {
	return left.second < right.second;
}


// simple function to detect all residues within a certain distance of the starting position
void detect_neighbors(
	const core::pose::Pose & p,
	const core::Size pos,
	const core::Real radius,
	std::map<core::Size, bool> & neighbor_map) {

	for ( core::Size i = 1; i <= p.size(); i++ ) {
		neighbor_map[i] = false;
		if ( pos != i ) {
			for ( core::Size a = 1; a <= p.residue(pos).natoms(); a++ ) {
				for ( core::Size b = 1; b <= p.residue(i).natoms(); b++ ) {
					core::Real dist = p.residue(pos).xyz(a).distance(p.residue(i).xyz(b));
					if ( dist <= radius ) {
						neighbor_map[i] = true;
						break;
					}
				}
			}
		}
	}
	return;
}


// modularized to simplify comparing structures from non-Rosetta-repacking
// WARNING: if the native structure doesn't match the repacked structure in size etc., we have a problem --> assertion will fail, program will quit
void compare_chi1_2_angles(
	const core::pose::Pose & native_p,
	const core::pose::Pose & repacked_p,
	const std::set<core::Size> & cluster,
	utility::vector1<int> & chi_dev_lst,
	utility::vector1<std::string> & chi_dev_details

)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ detect_tight_clusters::debug ]() ) {
		TR << native_p << std::endl;
		TR << repacked_p << std::endl;
	}

	runtime_assert( native_p.size() == repacked_p.size() );
	protocols::rotamer_recovery::RRComparerChiDiff rrc_chi_diff;
	rrc_chi_diff.set_max_chi_considered( 2 ); // only evaluate chi1 and chi2

	// make sure that chi_dev_details has the correct size
	chi_dev_details.resize(cluster.size());

	for ( auto it = option[ detect_tight_clusters::chi_deviation_thresholds ]().begin(), end = option[ detect_tight_clusters::chi_deviation_thresholds ]().end(); it != end; it++ ) {

		//const core::Size chi_deviation_threshold ( *it );
		rrc_chi_diff.set_recovery_threshold( *it );
		if ( option[ detect_tight_clusters::debug ]() ) {
			TR << " checking for chi dev threshold " << *it << std::endl;
		}

		// also count how many residues have chi1 and/or chi2 off beyond our threshold
		// iterate over cluster, not all positions, as this may include the environment if we repack more
		std::map <core::Size, int> pos_with_chi_dev;

		utility::vector1<core::Size> individual_chi_dev_vector(4); // this could be a vector of bools, but I think Size prints more reliably -- AS_DEBUG -- still needed?

		int i = 0; // size isn't random-access, so I'll need a counter to keep track of where we are -- this already looks error-prone
		for ( unsigned long iter : cluster ) {
			i++; // working with vector1
			core::conformation::Residue const & refres( native_p.residue( iter ) );
			core::conformation::Residue const & rp_res( repacked_p.residue( iter ) );
			bool within_chi_threshold = false; // the logic of the reporter is inverted vs. what I did before -- it reports true if the chi delta is within the threshold (!)
			core::Real max_chi_diff; // we're not reporting this at the moment, but the function wants it
			bool measured = rrc_chi_diff.measure_rotamer_recovery(native_p, repacked_p, refres, rp_res, max_chi_diff, within_chi_threshold);
			bool chi_threshold_exceeded = !within_chi_threshold;

			// record residue identity or fetch existing data
			std::string res_info;
			if ( chi_dev_details[i] != "" ) { // AS_DEBUG -- this might not work
				res_info = chi_dev_details[i];
			} else {
				res_info = refres.name1();
			}

			if ( !measured ) {
				TR << "WARNING -- chi diff could not be assessed for residues at position " << iter << std::endl;
			} else {
				if ( chi_threshold_exceeded ) {
					pos_with_chi_dev[iter]++;
				}
				//pos_with_chi_dev[*iter] = chi_threshold_exceeded; // I'm not sure if this is always 1 2 3 4... check -- this probably was replaced with res_info now -- AS_DEBUG remove?
				//TR << " checking individual chi thresholds: " << i << " " << chi_threshold_exceeded << std::endl; // AS_DEBUG
				res_info += ":" + utility::to_string(chi_threshold_exceeded);
			}
			chi_dev_details[i] = res_info; // store (partial) information about current residue
		}
		chi_dev_lst.push_back(pos_with_chi_dev.size());
	}
}

core::kinematics::MoveMapOP derive_MoveMap_from_cluster_lst(
	core::pose::Pose const & pose,
	utility::vector1< bool > const & is_flexible,
	bool allow_bb_move = false /* does this work if I don't have a separate declaration? */
)
{
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );

	movemap->set_bb(allow_bb_move);
	movemap->set_chi(false);
	movemap->set_jump(false);

	for ( core::Size i=1; i<= pose.size(); ++i ) {
		if ( is_flexible[i] ) {
			movemap->set_chi( i, true );
		}
	}
	return movemap;
}


void
sidechain_fastrelax(
	core::scoring::ScoreFunctionOP scorefxn,
	utility::vector1< bool > const & is_flexible,
	bool const cartesian_min,
	core::pose::Pose & pose
)
{
	protocols::relax::FastRelax fastrelax( scorefxn, 0 );
	fastrelax.cartesian( cartesian_min );

	core::kinematics::MoveMapOP movemap = derive_MoveMap_from_cluster_lst(pose, is_flexible, false); // should be modified to respect a flag for whether to also move the backbone

	fastrelax.set_movemap( movemap );
	fastrelax.apply( pose );
}


void repack_cluster(
	const core::pose::Pose & p, // this is the native/reference structure
	core::pose::Pose repacked, // this one will be modified (well, repacked)
	core::scoring::ScoreFunctionOP score_fxn,
	std::set<core::Size> cluster,
	const std::string cluster_acc, // buried or intermediate?
	std::string p_name,
	utility::io::ozstream & outfile) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	utility::vector1<core::Size> positions;
	utility::vector1<std::string> cluster_pos_in_PDB_numbering; // for use by external packers & evaluation
	utility::vector1<std::string> repacked_pos_in_PDB_numbering;
	// perform actual repacking, or just wrap around RMSD and chi diff calculations?
	if ( !option[ detect_tight_clusters::evaluate_externally_repacked_structure ]() ) {

		core::pack::task::PackerTaskOP base_packer_task( core::pack::task::TaskFactory::create_packer_task( p ));
		// we now need to tell the TaskFactory which parts are allowed to move

		base_packer_task->set_bump_check( true );
		base_packer_task->initialize_from_command_line();
		base_packer_task->or_include_current( false );

		core::pack::task::PackerTaskOP repack_packer_task( base_packer_task->clone() );
		utility::vector1<bool> allow_repacked( p.size(), false );

		utility::vector1<bool> ala_aalist( chemical::num_canonical_aas, false );
		ala_aalist[ chemical::aa_ala ] = true;


		for ( unsigned long iter : cluster ) {
			if ( !allow_repacked.at(iter) ) { // prevent duplicates -- in 8A repacking, this residue could already be in the list
				allow_repacked.at(iter) = true;
				positions.push_back(iter);
				repack_packer_task->nonconst_residue_task(iter).restrict_to_repacking();
				//TR << *iter << " " << p.pdb_info()->chain(*iter) << p.pdb_info()->number(*iter) << std::endl;

				repacked_pos_in_PDB_numbering.push_back(utility::to_string(p.pdb_info()->chain(iter)) + ":" + utility::to_string(p.pdb_info()->number(iter)));
			}

			cluster_pos_in_PDB_numbering.push_back(utility::to_string(p.pdb_info()->chain(iter)) + ":" + utility::to_string(p.pdb_info()->number(iter)));

			// in case we repack the entire 8A shell, we also need to determine all those residues and set their packer task appropriately
			if ( basic::options::option[ basic::options::OptionKeys::detect_tight_clusters::repack_8A_sphere ]() ) {
				std::map<core::Size, bool> neighbor_map;
				detect_neighbors(p, iter, repack_sphere_radius, neighbor_map);
				for ( std::map<core::Size, bool>::const_iterator m_iter = neighbor_map.begin(); m_iter != neighbor_map.end(); m_iter++ ) {
					if ( m_iter->second && !allow_repacked.at(m_iter->first) ) { // try to avoid duplicates
						allow_repacked.at((m_iter->first)) = true;
						positions.push_back(m_iter->first);
						repack_packer_task->nonconst_residue_task(m_iter->first).restrict_to_repacking();
						repacked_pos_in_PDB_numbering.push_back(utility::to_string(p.pdb_info()->chain(m_iter->first)) + ":" + utility::to_string(p.pdb_info()->number(m_iter->first)));

					}
				}
			}
		}

		if ( option[ detect_tight_clusters::debug ]() ) {
			TR << "repacked positions: ";
			for ( core::Size pos = 1; pos <= positions.size(); pos++ ) {
				TR << positions[pos] << ", ";
			}
			TR << std::endl;
		}

		repack_packer_task->restrict_to_residues( allow_repacked );

		// first operation is always standard packing -- only affected by the command line flags, and by which rotamers get to move
		// this sets the baseline and ensures that native side chains are discarded (except when running with -use_input_sc)
		core::pack::pack_rotamers( repacked, *score_fxn, repack_packer_task );

		// making cartmin available to all protocols
		bool const cartmin( option[ detect_tight_clusters::cartmin ] ); /// option: use cartesian-space min
		// Phil: you should also allow setting of just the cart_bonded subterms (angle/length/etc)
		if ( cartmin ) runtime_assert( score_fxn->get_weight( core::scoring::cart_bonded ) > 1e-3 );


		if ( option[ detect_tight_clusters::min_pack ] ) {
			core::pack::min_pack( repacked, *score_fxn, repack_packer_task );

		} else if ( option[ detect_tight_clusters::rot_trials ]() ) {
			for ( core::Size ri = 0; ri < rot_trials_iterations; ri++ ) {
				core::pack::rotamer_trials( repacked, *score_fxn, repack_packer_task );
			}
		} else if ( option[ detect_tight_clusters::rot_trials_min ]() ) {

			core::pack::RTMin rtmin;
			for ( core::Size ri = 0; ri < rot_trials_iterations; ri++ ) {
				rtmin.rtmin( repacked, *score_fxn, repack_packer_task );
			}

		} else if ( option[ detect_tight_clusters::sidechain_fastrelax ] ) {
			/// use Mike's fast-relax protocol, keeping the backbone fixed
			sidechain_fastrelax( score_fxn, allow_repacked, cartmin, repacked );

		} else if ( option[ detect_tight_clusters::soft_repack_min_hard_repack_min ]() ) { // soft-repack+hard-min+hard-repack+hard-min, as suggested by Frank

			// make soft scorefxn

			core::scoring::ScoreFunctionOP soft_sfxn;
			if ( option[ score::soft_wts ].user() ) {
				soft_sfxn = core::scoring::ScoreFunctionFactory::create_score_function( option[ score::soft_wts ]() );
			} else {
				soft_sfxn = core::scoring::ScoreFunctionFactory::create_score_function( "soft_rep_design" );
			}

			// repack with soft
			core::pack::pack_rotamers( repacked, *soft_sfxn, repack_packer_task );

			// minimize with hard -- requires movemap that must be derived from PackerTask
			core::kinematics::MoveMapOP mm = derive_MoveMap_from_cluster_lst(repacked, allow_repacked);
			protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover() ); // no symmetry support for now
			//const std::string min_type = "lbfgs_armijo_nonmonotone"; // use as many defaults as possible -- not sure if we want to fix this, actually
			//min_mover->min_type( min_type );
			min_mover->score_function( score_fxn );
			min_mover->cartesian( cartmin );
			min_mover->movemap( mm );
			min_mover->apply( repacked );

			// repack with hard -- probably requires IncludeCurrent
			core::pack::task::PackerTaskOP after_soft_packer_task( core::pack::task::TaskFactory::create_packer_task( p ));
			// we now need to tell the TaskFactory which parts are allowed to move

			after_soft_packer_task->set_bump_check( true );
			after_soft_packer_task->initialize_from_command_line();
			after_soft_packer_task->or_include_current( true ); // I think we want this, otherwise how to benefit from the previous step?
			after_soft_packer_task->restrict_to_repacking();
			after_soft_packer_task->restrict_to_residues( allow_repacked );
			core::pack::pack_rotamers( repacked, *score_fxn, after_soft_packer_task );

			// minimize with hard, same as above
			min_mover->apply( repacked );


		} else if ( option[ detect_tight_clusters::backrub ]() ) { // note that backrub only makes sense for 8A repacking

			core::pose::Pose best_backrub(repacked);
			protocols::moves::MonteCarlo mc(best_backrub, *score_fxn, 0.6/*option[ backrub::mc_kt ]*/);
			//core::pack::task::TaskFactoryOP backrub_task_factory = new core::pack::task::TaskFactory;
			//backrub_task_factory->push_back( new core::pack::task::operation::InitializeFromCommandline );
			//backrub_task_factory->modify_task(*after_backrub, repack_packer_task); // maybe this needs to be defined on the respective structure, even though the object is just a copy? will this work even though repack_packer_task was defined on another TF?

			for ( core::Size bi = 1; bi <= backrub_iterations; bi++ ) {

				TR << " -- pose setup -- " << std::endl;


				core::pose::PoseOP after_backrub( new core::pose::Pose(repacked) ); // freshly initialize each time
				//main_task_factory->modify_task(*after_backrub, repack_packer_task);
				mc.reset(*after_backrub); // reset MonteCarlo object


				TR << " -- backrub and sc mover setup -- " << std::endl;

				protocols::backrub::BackrubMover backrubmover;// setup -- copied from Colin's backrub_pilot.cc
				// read known and unknown optimization parameters from the database -- note that this is still for mm_bend, thus we'll keep the known issues with bond angle distribution...
				backrubmover.branchopt().read_database();
				/*
				if (energymethodoptions.bond_angle_residue_type_param_set()) {
				backrubmover.branchopt().bond_angle_residue_type_param_set(energymethodoptions.bond_angle_residue_type_param_set());
				}
				*/
				// do anything about detailed balance?

				// looks like the side chain mover needs to be set up separately
				protocols::simple_moves::sidechain_moves::SidechainMover sidechainmover;

				// dummy task factory that allows repacking everything
				//core::pack::task::TaskFactoryOP backrub_task_factory = new core::pack::task::TaskFactory;
				//backrub_task_factory->push_back( new core::pack::task::operation::InitializeFromCommandline );
				//core::pack::task::PackerTaskOP backrub_packer_task( main_task_factory->create_packer_task(repacked) );
				//backrub_packer_task->restrict_to_repacking();
				//backrub_task_factory->push_back( new core::pack::task::operation::RestrictToRepacking );


				//sidechainmover.set_task_factory(main_task_factory);
				sidechainmover.set_task(repack_packer_task);

				sidechainmover.set_prob_uniform( 0.1/* option[ backrub::sc_prob_uniform ]*/ );
				sidechainmover.set_prob_withinrot( 0.0/* option[ backrub::sc_prob_withinrot ] */);
				sidechainmover.set_preserve_detailed_balance(false/* option[ backrub::detailed_balance ]*/ );  // not for now -- would be difficult for Noah's stuff

				//sidechainmover.init_task(); // maybe this is what we're missing? -- no, this is supposed to happen inside the apply function

				// wait -- the side chain mover also needs to be restricted -- or is this done through the main_task_factory?


				// initial setup...
				backrubmover.clear_segments();
				backrubmover.set_input_pose(after_backrub);
				backrubmover.set_pivot_residues(positions);
				backrubmover.add_mainchain_segments(*after_backrub); // _from_options();

				// initialize structure (!?) -- check with Tanja...
				backrubmover.optimize_branch_angles(*after_backrub);
				sidechainmover.idealize_sidechains(*after_backrub);

				// do actual backrub...
				for ( int i = 1; i <= 1000; /* option[ backrub::ntrials ]; */ ++i ) {


					//TR << bi << " backrub/sc iteration " << i << std::endl; // debug


					std::string move_type;
					core::Real proposal_density_ratio(1); // only relevant for detailed balance, I think

					// could use random mover for this...
					core::Real move_prob = numeric::random::rg().uniform();
					if ( move_prob > option[ detect_tight_clusters::backrub_sc_prob ] ) {
						//TR << " backrub mover! " << move_prob << std::endl;

						// AS debug
						// can we get some info about the next segment here?
						//protocols::backrub::BackrubSegment const & segment(backrubmover.segment(backrubmover.next_segment_id()));
						//TR << "Next segment: " << segment.start_atomid() << segment.end_atomid() << std::endl;


						backrubmover.apply(*after_backrub);
						move_type = backrubmover.type();

						//protocols::backrub::BackrubSegment const & last_segment(backrubmover.segment(backrubmover.last_segment_id()));
						//TR << "Last segment was: " << last_segment.start_atomid() << last_segment.end_atomid() << std::endl;


						/*
						if (option[ backrub::backrub_sc_prob ] && segment.size() == 7) {
						core::Size middle_resnum(static_cast<core::Size>((segment.start_atomid().rsd() + segment.end_atomid().rsd())*.5));
						//TR << "Simultaneous move for segment: " << segment.start_atomid() << segment.end_atomid() << " middle: " << middle_resnum << std::endl;
						if (sidechainmover.residue_packed()[middle_resnum]) {
						if (numeric::random::rg().uniform() < option[ backrub::backrub_sc_prob ]) {
						sidechainmover.next_resnum(middle_resnum);
						//TR << "next_resnum before: " << sidechainmover.next_resnum() << std::endl;
						sidechainmover.apply(*pose);
						//TR << "next_resnum after: " << sidechainmover.next_resnum() << std::endl;
						move_type += "_" + sidechainmover.type();
						if (option[ backrub::detailed_balance ]) {
						proposal_density_ratio = sidechainmover.last_proposal_density_ratio();
						}
						//TR << "proposal density: " << proposal_density_ratio << std::endl;
						}
						}
						}
						*/
					} else {

						//TR << " side chain mover! " << move_prob << std::endl;

						sidechainmover.apply(*after_backrub);
						move_type = sidechainmover.type();
						/*
						if (option[ backrub::detailed_balance ]) {
						proposal_density_ratio = sidechainmover.last_proposal_density_ratio();
						}
						*/
					}

					//mdhist_proposed.record(pose->residue(1).chi());

					//TR << "boltzmann check " ;
					/*bool accept =*/ mc.boltzmann(*after_backrub, move_type, proposal_density_ratio); //  ?
					//TR << accept << std::endl;

				}
				//mc.show_counters();


				best_backrub = mc.lowest_score_pose(); // actually get lowest-scoring pose -- note that there's a risk that we never take anything other than the plain repacked pose (if everything else scores worse) -- is that what we want?
				score_fxn->show(TR, *after_backrub);
				TR.flush();

				if ( option[ detect_tight_clusters::debug ]() ) {
					TR << bi << " -- current best backrub model " << best_backrub << std::endl; // debug
				}


				/*

				maybe we don't need this -- maybe MC can do this for us?

				// accept based on energy, except for the first one -- otherwise we might end up without any actual backrub moves
				if (bi == 1 || (*score_fxn)(*after_backrub) < (*score_fxn)(best_backrub)) {

				best_backrub.clear();
				best_backrub = *after_backrub; //new core::pose::Pose(*after_backrub); // check if this works as a copy constructor, or if there's a problem once we exit the loop
				}


				*/

			}
			repacked.clear();
			repacked = best_backrub;

			if ( option[ detect_tight_clusters::debug ]() ) {
				TR << "backrubbed pose is now " << repacked << (*score_fxn)(repacked) << std::endl;
				repacked.dump_pdb( p_name+"_"+"backrub.pdb" );
			}
		}

	} // done with repacking (if applicable) -- the rest is evaluation

	core::pose::Pose starting_cluster;
	core::pose::Pose repacked_cluster;
	utility::vector1<core::Real> individual_pos_RMSDs; // all-atom RMSDs for each cluster position individually

	std::string cluster_name = "";
	std::map<core::Size, core::Size> all_neighbors; // make sure we don't append residues multiple times to a pose, that would probably be a mess
	for ( unsigned long iter : cluster ) {
		cluster_name += utility::to_string(p.pdb_info()->number(iter)) + "_";
		starting_cluster.append_residue_by_jump(p.residue(iter), starting_cluster.size());
		repacked_cluster.append_residue_by_jump(repacked.residue(iter), repacked_cluster.size());

		// get all-atom RMSD for this particular residue
		core::pose::Pose i_native; // empty pose, then just append the single residue - will this work?
		core::pose::Pose i_repacked;
		i_native.append_residue_by_jump(p.residue(iter), i_native.size());
		i_repacked.append_residue_by_jump(repacked.residue(iter), i_repacked.size());
		individual_pos_RMSDs.push_back(core::scoring::rmsd_with_super(i_native, i_repacked, core::scoring::is_heavyatom));


		// we'll also want these residues in the output set -- note that this will affect the RMSD calculation (!)
		if ( basic::options::option[ basic::options::OptionKeys::detect_tight_clusters::repack_8A_sphere ]() ) {
			std::map<core::Size, bool> neighbor_map;
			detect_neighbors(p, iter, repack_sphere_radius, neighbor_map);
			for ( std::map<core::Size, bool>::const_iterator m_iter = neighbor_map.begin(); m_iter != neighbor_map.end(); m_iter++ ) {
				if ( m_iter->second ) {
					if ( all_neighbors.find(m_iter->first) == all_neighbors.end() ) {
						starting_cluster.append_residue_by_jump(p.residue(m_iter->first), starting_cluster.size());
						repacked_cluster.append_residue_by_jump(repacked.residue(m_iter->first), repacked_cluster.size());
						all_neighbors[m_iter->first]++;
					}
				}
			}
		}
	}

	core::Real motif_starting_energy = 0.0;

	for ( core::Size i = 1; i <= positions.size(); i++ ) {
		for ( core::Size j = i+1; j <= positions.size(); j++ ) {
			core::scoring::EMapVector energymap;
			score_fxn->eval_ci_2b_sc_sc(p.residue(positions[i]), p.residue(positions[j]), p, energymap);
			score_fxn->eval_cd_2b_sc_sc(p.residue(positions[i]), p.residue(positions[j]), p, energymap);
			core::Real two_body_energy = score_fxn->weights().dot(energymap);
			motif_starting_energy += two_body_energy;
		}
	}

	core::Real motif_repacked_energy = 0.0;

	for ( core::Size i = 1; i <= positions.size(); i++ ) {
		for ( core::Size j = i+1; j <= positions.size(); j++ ) {
			core::scoring::EMapVector energymap;
			score_fxn->eval_ci_2b_sc_sc(repacked.residue(positions[i]), repacked.residue(positions[j]), repacked, energymap);
			score_fxn->eval_cd_2b_sc_sc(repacked.residue(positions[i]), repacked.residue(positions[j]), repacked, energymap);
			core::Real two_body_energy = score_fxn->weights().dot(energymap);
			motif_repacked_energy += two_body_energy;
		}
	}

	if ( basic::options::option[ basic::options::OptionKeys::detect_tight_clusters::generate_output_structures ]() ) {
		//starting_cluster.dump_pdb( p_name+"_"+cluster_name+"starting.pdb" );
		//repacked_cluster.dump_pdb( option[ out::path::pdb ]().path() + p_name+"_"+cluster_name+"repacked.pdb" );
		repacked_cluster.dump_pdb( p_name+"_"+cluster_name+"repacked.pdb" );// just the repacked residues
		//repacked.dump_pdb( p_name+"_"+cluster_name+"repacked.pdb" ); // full structure
	}

	core::Real motif_rmsd = core::scoring::rmsd_with_super(starting_cluster, repacked_cluster, core::scoring::is_heavyatom);

	utility::vector1<int> chi_dev_lst;
	utility::vector1<std::string> chi_dev_details;
	compare_chi1_2_angles(p, repacked, cluster, chi_dev_lst, chi_dev_details);

	outfile << p_name << "\t"
		<< cluster_name << "\t"
		<< motif_rmsd;

	for ( core::Size i = 1; i <= chi_dev_lst.size(); i++ ) {
		outfile << "\t" << chi_dev_lst[i];
	}

	// report failures per residue
	for ( core::Size i = 1; i <= chi_dev_details.size(); i++ ) {
		outfile << "\t" << chi_dev_details[i];
	}

	// report RMSD per residue -- note that I think these are not symmetry-corrected
	for ( core::Size i = 1; i <= individual_pos_RMSDs.size(); i++ ) {
		outfile << "\t" << individual_pos_RMSDs[i];
	}

	outfile << "\t" << cluster_acc << "\t" << p.energies().total_energy() << "\t"
		<< repacked.energies().total_energy() << "\t"
		<< motif_starting_energy << "\t"
		<< motif_repacked_energy << "\t"
		<< repacked_cluster.sequence() << "\t";

	for ( core::Size i = 1; i <= cluster_pos_in_PDB_numbering.size(); i++ ) {
		outfile << cluster_pos_in_PDB_numbering[i] << ",";
	}

	outfile << "\t";

	for ( core::Size i = 1; i <= repacked_pos_in_PDB_numbering.size(); i++ ) {
		outfile << repacked_pos_in_PDB_numbering[i] << ",";
	}

	outfile << std::endl;
}


// simple function to check for b-factor and buried-ness
// - also make sure that there are no Rosetta-rebuilt or occupancy-zero residues within 6A
// - also recored #neighbors, for exposure classification -- per cluster, keep min #neighbors over the 4 pos
bool passes_quality_check(
	const core::pose::Pose & p,
	core::Size i,
	core::Size & num_neighbors) {
	bool acceptable = true;
	core::Real b_factor_threshold = basic::options::option[ basic::options::OptionKeys::detect_tight_clusters::b_factor_threshold]();
	core::conformation::Residue res_i(p.residue(i));
	for ( core::Size ii=1; ii <= res_i.natoms(); ii++ ) {
		// also remove cases with a b-factor of 0 -- these are usually rebuilt side-chains, so they're not suitable for
		// a packing test -- and those with occupancy 0 -- however, virtual atoms are allowed to have 0 occ
		// added parentheses to avoid logic warning; hope I got it correct ~Labonte -- yes, thanks, looks fine -- AS
		if ( (p.pdb_info()->temperature(i, ii) > b_factor_threshold) ||
				(p.residue(i).atom_type(ii).is_heavyatom() &&
				!res_i.atom_type(ii).is_virtual() &&
				(p.pdb_info()->temperature(i, ii) == 0 || p.pdb_info()->occupancy(i, ii) == 0)) ) {
			if ( basic::options::option[ basic::options::OptionKeys::detect_tight_clusters::debug ]() ) {
				TR << " INITIAL FILTER -- discarding " << i;
				TR << " because its b-factor in atom " << p.residue(i).atom_type(ii).name();
				TR << " is too high (or 0): " << p.pdb_info()->temperature((i), ii) << std::endl;
			}
			acceptable = false;
			break;
		}
	}
	// filter for surface exposure -- each CB must have at least 9 other CB within an 8A radius
	// AS Feb 8, 2013: adapting for complexes with DNA or RNA
	//core::Size num_cb_neighbors = 0;
	std::map<core::Size, bool> cb_neighbors; // to make sure I only count each of them once
	numeric::xyzVector<double> cb_i;
	if ( !(p.residue(i).is_protein()) ) { // the residues in the cluster have to be protein residues for now -- could be changed later to also allow DNA in clusters
		if ( basic::options::option[ basic::options::OptionKeys::detect_tight_clusters::debug ]() ) {
			TR << " INITIAL FILTER -- discarding " << i << " because it isn't a protein residue: " << p.residue(i) << std::endl;
		}
		acceptable = false;
		//break;
	}
	if ( p.residue(i).name1() == 'G' ) { // assumption: we've checked b-factors above
		cb_i = p.residue(i).xyz(" CA ");
	} else {
		cb_i = p.residue(i).xyz(" CB ");
	}
	// iterate over all residues in this protein
	for ( core::Size pos = 1; pos <= p.size(); pos++ ) {
		// non-protein "residues" don't necessarily have a CB, which would cause the program to die -- for DNA, use CA
		if ( ( pos != i && p.residue(pos).is_protein() ) ||
				p.residue(pos).is_DNA() ||
				p.residue(pos).is_RNA() ) {

			// test for any closeby residue with 0 occupancy or 0 b-factor (which indicates Rosetta rebuilding)
			for ( core::Size ii=1; ii <= res_i.natoms(); ii++ ) {
				if ( res_i.atom_type(ii).is_heavyatom() && !res_i.atom_type(ii).is_virtual() && p.pdb_info()->temperature(i, ii) != 0 && p.pdb_info()->occupancy(i, ii) != 0 ) { // hydrogens will sometimes be rebuilt
					for ( core::Size jj=1; jj <= p.residue(pos).natoms(); jj++ ) {
						if ( res_i.xyz(ii).distance(p.residue(pos).xyz(jj)) < env_quality_check_dist  &&  p.residue(pos).atom_type(jj).is_heavyatom() && !p.residue(pos).atom_type(jj).is_virtual() && (p.pdb_info()->temperature(pos, jj) == 0 || p.pdb_info()->occupancy(pos, jj) == 0) ) {
							if ( basic::options::option[ basic::options::OptionKeys::detect_tight_clusters::debug ]() ) {
								TR << " Cluster environment filter for position " << i << " -- discarding " << pos << " because its b-factor or occupancy in atom " << p.residue(pos).atom_type(jj).name() << " are 0: " << p.pdb_info()->temperature(pos, jj) << " / " << p.pdb_info()->occupancy(pos, jj) << std::endl;
							}
							acceptable = false;
							break;
						}
					}
				}
			}

			// count CB neighbors
			numeric::xyzVector<double> cb_pos;
			if ( p.residue(pos).name1() == 'G' || p.residue(pos).is_DNA() || p.residue(pos).is_RNA() ) {
				cb_pos = p.residue(pos).xyz(" CA ");
				// check the temperature of this atom
				if ( p.pdb_info()->temperature((pos), p.residue(pos).atom_index(" CA ")) == 0 || p.pdb_info()->occupancy((pos), p.residue(pos).atom_index(" CA ")) == 0 ) {
					if ( basic::options::option[ basic::options::OptionKeys::detect_tight_clusters::debug ]() ) {
						TR << "INITIAL FILTER -- WARNING -- neighbor candidate " << pos << " has a Rosetta-rebuilt CA or occupancy, skipping... " << std::endl;
					}
					continue;
				}
			} else {
				cb_pos = p.residue(pos).xyz(" CB ");

				// check the temperature and occupancy of this atom
				if ( p.pdb_info()->temperature((pos), p.residue(pos).atom_index(" CB ")) == 0 || p.pdb_info()->occupancy((pos), p.residue(pos).atom_index(" CB ")) == 0 ) {
					if ( basic::options::option[ basic::options::OptionKeys::detect_tight_clusters::debug ]() ) {
						TR << "INITIAL FILTER -- WARNING -- neighbor candidate " << pos << " has a Rosetta-rebuilt or zero-occupancy CB, skipping... " << std::endl;
					}
					continue;
				}

			}

			if ( cb_i.distance(cb_pos) <= neighbor_dist ) {
				cb_neighbors[pos] = true;
			}
		}
	}
	num_neighbors = 0; // just to make sure...
	for ( std::map<core::Size, bool>::const_iterator mi = cb_neighbors.begin(); mi != cb_neighbors.end(); mi++ ) {
		if ( mi->second && mi->first != i ) { // mustn't count self in neighbor list
			num_neighbors++;
		}
	}
	if ( num_neighbors < basic::options::option[ basic::options::OptionKeys::detect_tight_clusters::min_num_neighbors]() ) {
		if ( basic::options::option[ basic::options::OptionKeys::detect_tight_clusters::debug ]() ) {
			TR << "INITIAL FILTER -- discarding " << i << " because it isn't buried enough: " << num_neighbors << std::endl;
		}

		acceptable = false; // surface cluster -- not buried enough
	}
	//num_neighbors = num_cb_neighbors; // will be returned by reference -- for later classification
	return acceptable;
}


void find_clusters(
	core::pose::Pose & p, // can't be const because of the TenANeighborGraph calculation? Or p.energies()? -- actually this shouldn't matter any more, as apparently we're not using the energygraph...
	core::pack::task::PackerTask & input_packer_task,
	std::map<std::set<core::Size>, std::string > & filtered_clusters) {

	using namespace basic::options;

	utility::vector1<std::set<core::Size> > neighbors;
	neighbors.resize(p.size());

	std::map<core::Size, core::Size> num_neighbors_by_pos;
	core::Real heavy_atom_dist = basic::options::option[ basic::options::OptionKeys::detect_tight_clusters::heavy_atom_dist]();
	bool find_interchain_clusters = basic::options::option[ basic::options::OptionKeys::detect_tight_clusters::detect_interchain_clusters]();

	//core::scoring::TenANeighborGraph & energygraph(p.energies().tenA_neighbor_graph());
	for ( core::Size i = 1; i <= p.size(); i++ ) {
		if ( p.residue(i).is_protein() ) {

			// check if this residue is allowed to be in a cluster
			core::pack::task::ResidueLevelTask const & res_i_task(input_packer_task.residue_task(i)); // doesn't seem to work...
			//TR << res_i_task.command_string() << std::endl;
			if ( !res_i_task.being_packed() ) {
				if ( basic::options::option[ basic::options::OptionKeys::detect_tight_clusters::debug ]() ) {
					TR << "Residue " << i << " was prevented from being in the cluster by the input resfile" << std::endl;
				}
				continue;
			}


			core::conformation::Residue res_i(p.residue(i));

			core::Size num_neighbors_i = 0; // will be set by passes_quality_check
			if ( !passes_quality_check(p, i, num_neighbors_i) ) {
				if ( basic::options::option[ basic::options::OptionKeys::detect_tight_clusters::debug ]() ) {
					TR << "Residue " << i << " didn't pass the quality check" << std::endl;
				}
				continue;
			}

			// TODO: for each cluster, store the minimum number of neighbors, so that we can use it for evaluation purposes -- maybe directly map to "intermediate" and "buried"
			num_neighbors_by_pos[i] = num_neighbors_i;

			//TR << " actually considering residue " << i << " for a cluster... " << std::endl;

			std::string name_i = utility::to_string(p.pdb_info()->number(i));
			for ( core::Size j = i+1; j <= p.size(); j++ ) {
				if ( p.residue(j).is_protein() ) {
					core::Size num_neighbors_j = 0;
					if ( !passes_quality_check(p, j, num_neighbors_j) ) {
						if ( basic::options::option[ basic::options::OptionKeys::detect_tight_clusters::debug ]() ) {
							TR << "Residue " << j << " didn't pass the quality check" << std::endl;
						}
						continue;
					}

					num_neighbors_by_pos[j] = num_neighbors_j;
					std::string name_j = utility::to_string(p.pdb_info()->number(j));
					core::conformation::Residue res_j(p.residue(j));
					//core::Real min_dist = 0.0;
					if ( find_interchain_clusters || res_j.chain() == res_i.chain() ) { // AS March 06, 2013: when looking for interchain clusters we don't care about chain comparison here
						bool is_neighbor = false;
						for ( core::Size ii=1; ii <= res_i.natoms(); ii++ ) {
							if ( !res_i.atom_is_backbone(ii) && !res_i.atom_is_hydrogen(ii) ) {
								for ( core::Size jj=1; jj <= res_j.natoms(); jj++ ) {
									if ( !res_j.atom_is_backbone(jj) && !res_j.atom_is_hydrogen(jj) ) { // maybe replace by .is_heavyatom() if we still see differences between the different runs?
										core::Real const distance(res_i.xyz(ii).distance(res_j.xyz(jj)));
										if ( distance <= heavy_atom_dist ) {
											is_neighbor = true;
											//min_dist = distance;
											break;
										}
									}
								}
							}
							if ( is_neighbor ) break;
						}
						if ( is_neighbor ) {
							neighbors[i].insert(j);
						}
					}
				}
			}
		}
	}

	//TR << "candidates (neigbors): " << neighbors.size() << std::endl; // debug

	std::map<std::set<core::Size>, std::string > clusters;

	//  for(core::Size i = 1; i <= p.size(); i++) {
	//   for(std::set<core::Size>::iterator iter_j = neighbors[i].begin(); iter_j != neighbors[i].end(); iter_j++) {
	//    core::Size j = (*iter_j);
	//   }
	//  }

	//find clusters of four interacting residues
	for ( core::Size i = 1; i <= p.size(); i++ ) {
		// check if there is shared neighbor between positions i and j
		for ( auto iter_j = neighbors[i].begin(); iter_j != neighbors[i].end(); iter_j++ ) {
			core::Size j = (*iter_j);
			for ( auto iter_k = neighbors[j].begin(); iter_k != neighbors[j].end(); iter_k++ ) {
				core::Size k = (*iter_k);
				if ( neighbors[i].find(k) != neighbors[i].end() ) {
					for ( auto iter_m = neighbors[k].begin(); iter_m != neighbors[k].end(); iter_m++ ) {
						core::Size m = (*iter_m);
						if ( ( !find_interchain_clusters && p.residue(i).chain() == p.residue(j).chain() && p.residue(i).chain() == p.residue(k).chain() && p.residue(i).chain() == p.residue(m).chain() ) ||
								( find_interchain_clusters && ( p.residue(i).chain() != p.residue(j).chain() || p.residue(i).chain() != p.residue(k).chain() || p.residue(i).chain() != p.residue(m).chain() || p.residue(j).chain() != p.residue(k).chain() || p.residue(j).chain() != p.residue(m).chain() || p.residue(k).chain() != p.residue(m).chain() ) ) ) { // default: only looking for intrachain clusters -- if detect_interchain_clusters(), some chains must (!) differ


							// add debug statement here indicating the different positions and their chain?
							//TR << "candidate positions: " << i << " " << j << " " << k << " " << m << std::endl;

							if ( neighbors[i].find(m) != neighbors[i].end() && neighbors[j].find(m) != neighbors[j].end() ) { // to fulfill the criteria of both neighbors being close enough, we can only check for m here, but not k

								//found a shared neighbor!

								//TR << "found a new cluster for consideration! chains are " << p.residue(i).chain() << p.residue(j).chain() << p.residue(k).chain() << p.residue(m).chain()  << std::endl; // debug


								// check if this is a new cluster
								bool already_found = false;
								//for(core::Size c=1; c <= clusters.size(); c++) {
								for ( std::map<std::set<core::Size>, std::string>::const_iterator mi = clusters.begin(); mi != clusters.end(); mi++ ) {
									if ( (mi->first).find(i) != (mi->first).end()
											&& (mi->first).find(j) !=(mi->first).end()
											&& (mi->first).find(k) !=(mi->first).end()
											&& (mi->first).find(m) !=(mi->first).end() ) {
										already_found = true;
										break;
									}
								}
								if ( !already_found ) {
									std::set<core::Size> new_cluster;
									new_cluster.insert(i);
									new_cluster.insert(j);
									new_cluster.insert(k);
									new_cluster.insert(m);

									// determine buried-ness
									// this may need to be adapted for the interchain clusters
									if ( num_neighbors_by_pos[i] > buried_threshold && num_neighbors_by_pos[j] > buried_threshold && num_neighbors_by_pos[k] > buried_threshold && num_neighbors_by_pos[m] > buried_threshold ) {
										clusters[new_cluster] = "buried";
									} else {
										// if (num_neighbors_by_pos[i] >= min_num_neighbors && num_neighbors_by_pos[j] >= min_num_neighbors && num_neighbors_by_pos[k] >= min_num_neighbors && num_neighbors_by_pos[m] >= min_num_neighbors) {
										//TR << " -- this should be an intermediate cluster: " << i << "/" << num_neighbors_by_pos[i] << " " << j << "/" << num_neighbors_by_pos[j] << " " << k << "/" << num_neighbors_by_pos[k] << " " << m << "/" << num_neighbors_by_pos[m] << std::endl;

										clusters[new_cluster] = "intermediate";
									}
									//clusters.push_back(new_cluster);
								}
							}
						}
					}
				}
			}
		}
	}


	// TR << "Initial number of clusters: " << clusters.size() << std::endl;
	filtered_clusters = clusters; // all filtering is already done in initial step -- "forbidden" residues aren't considered
	TR << "Filtered number of clusters: " << filtered_clusters.size() << std::endl;

}


int main( int argc, char * argv [] )
{
	try{

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using core::pack::task::operation::TaskOperationCOP;

		// register specific options
		OPT(packing::resfile);
		OPT(in::file::native); // for comparison against externally repacked structures
		NEW_OPT(detect_tight_clusters::generate_output_structures, "generate output structures (default: no)", false);
		NEW_OPT(detect_tight_clusters::require_renumbered_structures, "require that input structures be renumbered sequentially, to match Rosetta numbering? Hint: necessary for compatibility with external tools that rely on this numbering. (default: no)", false);
		NEW_OPT(detect_tight_clusters::debug, "debug mode? (tons of extra output)", false);
		NEW_OPT(detect_tight_clusters::repack_8A_sphere, "repack all residues in an 8A sphere around the cluster, instead of just the 4 cluster residues", false);
		// NEW_OPT(detect_tight_clusters::cluster_file_suffix, "extension for the filename with the cluster results", "test");
		NEW_OPT(detect_tight_clusters::chi_deviation_thresholds, "list of chi1/2 deviation thresholds", utility::vector1<core::Size>()); // -- haven't figured out how to give this default weights: 40, 20, 10));
		NEW_OPT(detect_tight_clusters::min_pack, "use MinPacker", false);
		//NEW_OPT(detect_tight_clusters::stochastic_pack, "use stochastic pack", false);
		NEW_OPT(detect_tight_clusters::rot_trials_min, "use rotamer trials with minimization after initial repack (10 trials)", false);
		NEW_OPT(detect_tight_clusters::rot_trials, "use rotamer trials after initial repack", false);
		NEW_OPT(detect_tight_clusters::sidechain_fastrelax, "use sidechain fast-relax protocol", false);
		NEW_OPT(detect_tight_clusters::cartmin, "Use cartesian space minimization. Only for sidechain-fast-relax protocol right now. Score function must have non-zero weight for cart_bonded term.", false);
		NEW_OPT(detect_tight_clusters::backrub, "use backrub to move both backbone and side chain", false);
		NEW_OPT(detect_tight_clusters::backrub_sc_prob, "probability of making a side chain move within backrub (default 0.25)", 0.25);
		NEW_OPT(detect_tight_clusters::external_cluster, "provide all cluster positions in PDB numbering (chain pos vector)", utility::vector1<std::string>()); // half of these will be casted to chars, the other half to ints...
		NEW_OPT(detect_tight_clusters::evaluate_externally_repacked_structure, "evaluate chi1/2 of an externally repacked structure -- requires -in:file:native for the reference structure", false);
		NEW_OPT(detect_tight_clusters::repack_external_cluster, "repack and evaluate an externally provided cluster, e.g. from a baseline run", false);
		NEW_OPT(detect_tight_clusters::detect_interchain_clusters, "detect clusters between different protein chains (default only detects clusters within a single chain)", false);
		NEW_OPT(detect_tight_clusters::min_num_neighbors, "minimum number of neighbors each residue in a cluster must have for the cluster to be accepted (default 9 -- good for within-protein clusters, but too high for interfaces)", 9);
		NEW_OPT(detect_tight_clusters::b_factor_threshold, "maximum allowed b-factor in any cluster residue atom", 30);
		NEW_OPT(detect_tight_clusters::heavy_atom_dist, "maximum heavy atom distance from each cluster residue to at least 2 other residues in the cluster", 4.5);
		NEW_OPT(detect_tight_clusters::soft_repack_min_hard_repack_min, "(hard repack to normalize, as for all runs, then) soft repack, hard minimize, hard repack w/ include_current, hard min", false);


		// initialize Rosetta
		devel::init(argc, argv);

		std::stringstream cluster_filename_suffix;
		if ( basic::options::option[ basic::options::OptionKeys::detect_tight_clusters::repack_8A_sphere ]() ) {
			cluster_filename_suffix << ".8A_sphere";
		}

		/*
		if (option[ detect_tight_clusters::cluster_file_suffix ].user()) {
		cluster_filename_suffix << "." << option[ detect_tight_clusters::cluster_file_suffix ]();
		TR << "setting cluster file suffix: " << cluster_filename_suffix << std::endl;
		}
		*/
		//cluster_filename_suffix << "." << option[ detect_tight_clusters::chi_deviation_threshold ]();
		for ( int it : option[ detect_tight_clusters::chi_deviation_thresholds ]() ) {
			cluster_filename_suffix << "." << it;
		}

		// get scoring function
		core::scoring::ScoreFunctionOP score_fxn = core::scoring::get_score_function();

		utility::vector1<std::string> pdbs = basic::options::option[ in::file::s ]();
		core::pose::Pose p;
		std::string p_name = pdbs[1]; // warning -- this just handles one PDB...

		core::import_pose::pose_from_file( p, p_name , core::import_pose::PDB_file);

		std::string outfile_core = p_name;


		//p.dump_pdb( p_name+"_native.pdb" ); // only for debugging the odd rebuild-neighbor-count issues


		if ( option[ out::path::pdb ].user() ) { // ?
			size_t si = outfile_core.rfind("/"); // last occurrence of /
			if ( si != std::string::npos ) {
				outfile_core = outfile_core.substr(si+1); // // we only want the name, but leave out the input path
			}
			outfile_core = option[ out::path::pdb ]().path() + outfile_core;
		} else if ( option[ out::path::all ].user() ) {
			size_t si = outfile_core.rfind("/"); // last occurrence of /
			if ( si != std::string::npos ) {
				outfile_core = outfile_core.substr(si+1); // // we only want the name, but leave out the input path
			}
			outfile_core = option[ out::path::all ]().path() + outfile_core;
		}

		// later: add option to compare against in::file::native, in case we want to (fast?)relax structures first [ideally restricting side chain conformations, I think there was some protocol for this...] and use those as starting structures, but then compare against the actual native/crystal

		(*score_fxn)(p);
		score_fxn->show(TR, p);
		TR << std::endl;

		// read in a resfile, if specified -- with this the user may disable regions from cluster detection, e.g. those that are too close to a HETATM (which may not be visible in runs with -ignore_unrecognized_res)
		core::pack::task::TaskFactoryOP input_task_factory( new core::pack::task::TaskFactory );
		input_task_factory->push_back( TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ) );
		if ( option[ packing::resfile ].user() ) {
			input_task_factory->push_back( TaskOperationCOP( new core::pack::task::operation::ReadResfile ) );
		} else {
			core::pack::task::operation::RestrictToRepackingOP rtrop( new core::pack::task::operation::RestrictToRepacking ); // make sure that by default everything is allowed to be repacked
			input_task_factory->push_back( rtrop );
		}

		//utility::vector1<std::set<core::Size> > filtered_clusters;
		std::map< std::set<core::Size>, std::string > filtered_clusters; // also store whether the cluster is intermediate or exposed

		core::pack::task::PackerTaskOP input_task(input_task_factory->create_task_and_apply_taskoperations(p)); // this seems to be a really inefficient way do to this...

		utility::io::ozstream outfile( outfile_core + ".cluster" + cluster_filename_suffix.str()); // this should be path-adapted


		if ( option[ detect_tight_clusters::debug ]() ) {
			p.dump_pdb( p_name+"_input.pdb" ); // to make sure Rosetta doesn't do anything funny to OSCAR's structure -- in particular, look at the Histidines!
		}

		pose::Pose native_p(p); // actual native if provided, otherwise this is a copy of the input structure
		(*score_fxn)(native_p);

		if ( option[ detect_tight_clusters::external_cluster ].user() ) {
			if ( option[ detect_tight_clusters::external_cluster ]().size() % 2 != 0 ) {
				TR << "error -- -external_cluster parameters must always be pairs of chain and position!" << std::endl;
				return 2;
			}


			if ( option [ detect_tight_clusters::evaluate_externally_repacked_structure ]() ) { // requires native for reference
				// it might be easier to directly call comparison on the given cluster
				if ( ! option[ in::file::native ].user() ) {
					TR << "error -- native pose required for comparison with externally repacked structures!" << std::endl;
					return 2;
				}
				native_p.clear();
				core::import_pose::pose_from_file(native_p, option[ in::file::native ](), core::import_pose::PDB_file); // does this work?
				(*score_fxn)(native_p); // initial scoring is required to get actual energies
			}

			std::set<core::Size> ext_cluster;
			std::string cluster_name = "";
			// also extract the actual 4 residues for getting the RMSD
			core::pose::Pose native_cluster;
			core::pose::Pose repacked_cluster;

			for ( core::Size c_pos = 1; c_pos <= option[ detect_tight_clusters::external_cluster ]().size(); c_pos += 2 ) {

				core::Size mapped_pos = p.pdb_info()->pdb2pose(option[ detect_tight_clusters::external_cluster ][c_pos][0], atoi(option[ detect_tight_clusters::external_cluster ][c_pos+1].c_str()));

				ext_cluster.insert(mapped_pos); // we have ensured above that there's an even number of entries, so this shouldn't produce segfaults -- the casting from strings to chars and ints is a terrible hack though
				cluster_name += utility::to_string(p.pdb_info()->number(mapped_pos)) + "_";

				native_cluster.append_residue_by_jump(native_p.residue(mapped_pos), native_cluster.size());
				repacked_cluster.append_residue_by_jump(p.residue(mapped_pos), repacked_cluster.size());


				if ( option[ detect_tight_clusters::debug ]() ) {
					TR << " *** evaluating external cluster *** " << option[ detect_tight_clusters::external_cluster ][c_pos]
						<< " " << option[ detect_tight_clusters::external_cluster ][c_pos+1] << " "
						<< mapped_pos << std::endl;
				}


			}

			if ( option[ detect_tight_clusters::evaluate_externally_repacked_structure ]() || option[ detect_tight_clusters::repack_external_cluster ]() ) {

				// append to filtered_clusters in either case -- repacking or just evaluation will be checked inside the function
				filtered_clusters[ext_cluster] = "external"; // we don't know whether this one is buried or intermediate
			} else {
				TR << " -- error -- you provided an external cluster but did not state whether to repack or just evaluate -- exiting" << std::endl;
			}

		} else {
			if ( option[ detect_tight_clusters::require_renumbered_structures]() ) { // switch on for compatibility with external tools -- there must be a better way to check for missing and thus discarded backbones though. The problem with this approach is that the renumbering affects which residues are excluded due to proximity to a ligand (via resfiles), so it may actually lead to exclusion of the wrong residues.

				// check here whether p.size() is the same as the last number in this PDB -- with -remember_unrecognized_residues this might show us which cases had "partial" backbone information that is discarded by Rosetta and would thus confuse numbering when used with other tools, such as OSCAR-star
				core::Size num_res = p.size();
				if ( num_res != core::Size(p.pdb_info()->number(num_res)) ) { // throws a warning if not casted..
					TR << p << std::endl;
					TR << "Structure " << p_name << " seems to have missing backbone data that was discarded -- skipping" << std::endl; // mainly for compatibility with OSCAR*
					return 0;
				}
			}

			find_clusters(p, *input_task, filtered_clusters);
		}


		// repack all the filtered clusters (whether externally provided or detected in this run)
		for ( std::map<std::set<core::Size>, std::string>::const_iterator mi = filtered_clusters.begin(); mi != filtered_clusters.end(); mi++ ) {
			repack_cluster(native_p, p, score_fxn, mi->first, mi->second, outfile_core, outfile);
		}

		outfile.close();

		//score_fxn->show(TR, p);


	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
