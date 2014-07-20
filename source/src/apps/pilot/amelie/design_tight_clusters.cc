// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Amelie Stein (amelie.stein@ucsf.edu)
/// @brief Given 4 residues and a repack_radius, design those 4 residues and report whether the native position was recovered

// Protocols Headers
//#include <protocols/rigid/RigidBodyMover.hh>

// Core Headers
// AS -- it's quite possible that a lot of the includes aren't needed...
#include <devel/init.hh>
#include <core/import_pose/import_pose.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <core/id/types.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDB_Info.hh>
#include <core/conformation/Residue.hh>
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
#include <protocols/rotamer_recovery/RRComparer.hh>
#include <protocols/backrub/BackrubMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/protein_interface_design/design_utils.hh> // for FavorNativeResidue

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
//const core::Size b_factor_threshold = 30;
//core::Real heavy_atom_dist = 4.5;
// const core::Size chi_deviation_threshold = 40; // for now we're only reporting really bad repacking, deviations of 40+ in chi1 or chi2 -- might be interesting to control this via a flag though
const core::Size rot_trials_iterations = 10;

basic::Tracer TR("apps.pilot.amelie.design_tight_clusters");
static numeric::random::RandomGenerator RG(45367418);

OPT_1GRP_KEY(Boolean, design_tight_clusters, generate_output_structures) // (no semicolon here!) -- this would generate huge amounts of output, so by default it runs silently, only reporting which residues are involved in clusters
OPT_1GRP_KEY(Boolean, design_tight_clusters, require_renumbered_structures) // requires numbering to match length of structures -- hack for making the output work with external tools like OSCARstar, but has some drawbacks (see below)
OPT_1GRP_KEY(Boolean, design_tight_clusters, debug)
OPT_1GRP_KEY(Real, design_tight_clusters, repack_radius) // something between 0 and 8 (more is probably just too much / more noise than signal)
OPT_1GRP_KEY(Integer, design_tight_clusters, num_cycles) // number of design/repack cycles -- note that some packing options cannot do design (e.g. sc_fr)
OPT_1GRP_KEY(Real, design_tight_clusters, favor_native) // energy bonus for native AA

//OPT_1GRP_KEY(Integer, design_tight_clusters, chi_deviation_threshold)
OPT_1GRP_KEY( IntegerVector, design_tight_clusters, chi_deviation_thresholds ) // keep list of chi_dev_thresholds and iterate over those -- note that these will only be relevant here if we rediscovered the native residue

OPT_1GRP_KEY(Boolean, design_tight_clusters, min_pack)
//OPT_1GRP_KEY(Boolean, design_tight_clusters, stochastic_pack)
OPT_1GRP_KEY(Boolean, design_tight_clusters, rot_trials_min)
OPT_1GRP_KEY(Boolean, design_tight_clusters, rot_trials)
OPT_1GRP_KEY(Boolean, design_tight_clusters, sidechain_fastrelax)
OPT_1GRP_KEY(Boolean, design_tight_clusters, cartmin)
OPT_1GRP_KEY(Boolean, design_tight_clusters, soft_repack_min_hard_repack_min)
OPT_1GRP_KEY(Boolean, design_tight_clusters, backrub)
OPT_1GRP_KEY(Real, design_tight_clusters, backrub_sc_prob) // probability of making a side chain move within the backrub option -- default 0.25
/*
 the following set of options os for providing an external cluster (e.e., from a baseline run)
 that can either be repacked and evaluated, or the user provides an externally repacked structure,
 which needs to come with a reference structure via -in:file:native

 a native structure can also be provided e.g. for repacking after backbone relaxation
 -- in this case we may first want to check whether chi1/2 of the cluster residues are
 still the same as in the native structure, but this isn't implemented yet
 */

OPT_1GRP_KEY(StringVector, design_tight_clusters, external_cluster) // PDB numbering -- positions will be casted to int, probably won't be able to deal well with insertion codes etc.... but when starting the initial (cluster detection) run from a renumbered structure this should be fine
//OPT_1GRP_KEY(Boolean, design_tight_clusters, evaluate_externally_repacked_structure)

bool sort_min(const std::pair<core::Size, core::Real> & left, const std::pair<core::Size, core::Real> & right) {
	return left.second < right.second;
}


// simple function to detect all residues within a certain distance of the starting position
void detect_neighbors(
					  const core::pose::Pose & p,
					  const core::Size pos,
					  const core::Real radius,
					  std::map<core::Size, bool> & neighbor_map) {

	for (core::Size i = 1; i <= p.total_residue(); i++) {
		neighbor_map[i] = false;
		if (pos != i) {
			for (core::Size a = 1; a <= p.residue(pos).natoms(); a++) {
				for (core::Size b = 1; b <= p.residue(i).natoms(); b++) {
					core::Real dist = p.residue(pos).xyz(a).distance(p.residue(i).xyz(b));
					if (dist <= radius) {
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
// NOTE: this now also checks whether the residue identity has changed or not, as we're testing design -- report both the native and the target residue, so that we can plot them in a heat map
void compare_residues_and_chi1_2_angles(
										const core::pose::Pose & native_p,
										const core::pose::Pose & repacked_p,
										const std::set<core::Size> & cluster,
										utility::vector1<std::string> & mutations,
										utility::vector1<std::string> & chi_dev_details

										)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ design_tight_clusters::debug ]() ) {
		TR << native_p << std::endl;
		TR << repacked_p << std::endl;
	}

	runtime_assert( native_p.total_residue() == repacked_p.total_residue() );
	protocols::rotamer_recovery::RRComparerChiDiff rrc_chi_diff;
	rrc_chi_diff.set_max_chi_considered( 2 ); // only evaluate chi1 and chi2

	// make sure that chi_dev_details has the correct size
	chi_dev_details.resize(cluster.size());
	mutations.clear();
	mutations.resize(cluster.size()); // AS_DEBUG -- for some reason this variable contains garbage, trying to get rid of it...

	for (IntegerVectorOption::const_iterator it = option[ design_tight_clusters::chi_deviation_thresholds ]().begin(), end = option[ design_tight_clusters::chi_deviation_thresholds ]().end(); it != end; it++) {

		//const core::Size chi_deviation_threshold ( *it );
		rrc_chi_diff.set_recovery_threshold( *it );
		if ( option[ design_tight_clusters::debug ]() ) {
			TR << " checking for chi dev threshold " << *it << std::endl;
		}

		// also count how many residues have chi1 and/or chi2 off beyond our threshold
		// iterate over cluster, not all positions, as this may include the environment if we repack more
		std::map <core::Size, int> pos_with_chi_dev;

		utility::vector1<core::Size> individual_chi_dev_vector(4); // this could be a vector of bools, but I think Size prints more reliably -- AS_DEBUG -- still needed?

		int i = 0; // size isn't random-access, so I'll need a counter to keep track of where we are -- this already looks error-prone
		for(std::set<core::Size>::iterator iter = cluster.begin(); iter != cluster.end(); iter++) {
		  i++; // working with vector1
		  core::conformation::Residue const & refres( native_p.residue( *iter ) );
		  core::conformation::Residue const & rp_res( repacked_p.residue( *iter ) );
		  //	  mutations.push_back(refres.name1() + ":" + rp_res.name1()); // AS_DEBUG -- this seems to access some weird part of the code instead of .name1()...
		  mutations[i] = utility::to_string(refres.name1()) + ":" + utility::to_string(rp_res.name1()); // AS_DEBUG -- this seems to access some weird part of the code instead of .name1()...

		  //  TR << " -- debug info within function compare_residues_and_chi12_angles -- " << i << " " << refres.name1() << "->" << rp_res.name1() << " *" << mutations[mutations.size()] << "* " << mutations.size() << std::endl; // AS_DEBUG


		  // record residue identity or fetch existing data
		  std::string res_info;
		  if (chi_dev_details[i] != "") // AS_DEBUG -- this might not work
		    res_info = chi_dev_details[i];
		  else
		    res_info = refres.name1();




		  if (refres.name1() == rp_res.name1()) {
		    bool within_chi_threshold = false; // the logic of the reporter is inverted vs. what I did before -- it reports true if the chi delta is within the threshold (!)
		    core::Real max_chi_diff; // we're not reporting this at the moment, but the function wants it
		    bool measured = rrc_chi_diff.measure_rotamer_recovery(native_p, repacked_p, refres, rp_res, max_chi_diff, within_chi_threshold);
		    bool chi_threshold_exceeded = !within_chi_threshold;

		    if (!measured) {
		      TR << "WARNING -- chi diff could not be assessed for residues at position " << *iter << std::endl;
		    } else {
		      if ( chi_threshold_exceeded ) {
			pos_with_chi_dev[*iter]++;
		      }
		      //pos_with_chi_dev[*iter] = chi_threshold_exceeded; // I'm not sure if this is always 1 2 3 4... check -- this probably was replaced with res_info now -- AS_DEBUG remove?
		      //TR << " checking individual chi thresholds: " << i << " " << chi_threshold_exceeded << std::endl; // AS_DEBUG
		      res_info += ":" + utility::to_string(chi_threshold_exceeded);
		    }
		  } else {
		    res_info += "?"; // for mutations we cannot determine the chi differences
		  }


		  //TR << res_info << std::endl; // AS_DEBUG

		  chi_dev_details[i] = res_info; // store (partial) information about current residue
		}
	}
}

core::kinematics::MoveMapOP derive_MoveMap_from_cluster_lst(
															core::pose::Pose const & pose,
															utility::vector1< bool > const & is_flexible,
															bool allow_bb_move = false
															)
{
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;

	movemap->set_bb(allow_bb_move);
	movemap->set_chi(false);
	movemap->set_jump(false);

	for ( core::Size i=1; i<= pose.total_residue(); ++i ) {
		if ( is_flexible[i] ) {
			movemap->set_chi( i, true );
		}
	}
	return movemap;
}


/// @brief use fastrelax to optimize packing
/// @details note that fastrelax doesn't seem to allow design, so we'll need separate design and packing steps for this
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


/// @brief as some of our packing functions cannot do design (fastrelax, backrub with restrictions) we need to have the option to do packing independent of / interleaved with design
/// @details assumption: the packer_task has already been set up correctly and knows which residues may be designed and which should only be packed
void repack_step (
				  core::pose::Pose & repacked, // will be repacked and thus modified
				  const core::pack::task::PackerTaskOP repack_packer_task,
				  const core::scoring::ScoreFunctionOP score_fxn,
				  utility::vector1<bool> allow_repacked
				  )
{
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  bool const cartmin( option[ design_tight_clusters::cartmin ] );

	if ( option[ design_tight_clusters::min_pack ] ) {
		core::pack::min_pack( repacked, *score_fxn, repack_packer_task );

	} else if ( option[ design_tight_clusters::rot_trials ]() ) {
		for (core::Size ri = 0; ri < rot_trials_iterations; ri++) {
			core::pack::rotamer_trials( repacked, *score_fxn, repack_packer_task );
		}
	} else if ( option[ design_tight_clusters::rot_trials_min ]() ) {

		core::pack::RTMin rtmin;
		for (core::Size ri = 0; ri < rot_trials_iterations; ri++) {
			rtmin.rtmin( repacked, *score_fxn, repack_packer_task );
		}

	} else if ( option[ design_tight_clusters::sidechain_fastrelax ] ) {
		/// use Mike's fast-relax protocol, keeping the backbone fixed
		sidechain_fastrelax( score_fxn, allow_repacked, cartmin, repacked );

	} else if ( option[ design_tight_clusters::soft_repack_min_hard_repack_min ]() ) { // soft-repack+hard-min+hard-repack+hard-min, as suggested by Frank

		// make soft scorefxn

		core::scoring::ScoreFunctionOP soft_sfxn;
		if ( option[ score::soft_wts ].user() )
			soft_sfxn = core::scoring::ScoreFunctionFactory::create_score_function( option[ score::soft_wts ]() );
		else
			soft_sfxn = core::scoring::ScoreFunctionFactory::create_score_function( "soft_rep_design" );

		// repack with soft
		core::pack::pack_rotamers( repacked, *soft_sfxn, repack_packer_task );

		// minimize with hard -- requires movemap that must be derived from PackerTask
		core::kinematics::MoveMapOP mm = derive_MoveMap_from_cluster_lst(repacked, allow_repacked);
		protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover(); // no symmetry support for now
		//const std::string min_type = "dfpmin"; // use as many defaults as possible -- not sure if we want to fix this, actually
		//min_mover->min_type( min_type );
		min_mover->score_function( score_fxn );
		min_mover->cartesian( cartmin );
		min_mover->movemap( mm );
		min_mover->apply( repacked );

		// repack with hard -- probably requires IncludeCurrent -- note that here, unlike in the standard repacking catastrophes, include_current is already true, so we don't need an extra packer task
		core::pack::pack_rotamers( repacked, *score_fxn, repack_packer_task );

		// minimize with hard, same as above
		min_mover->apply( repacked );



	} else if ( option[ design_tight_clusters::backrub ]() ) {// note that backrub only makes sense for 8A repacking

	  // removed for now -- we're not using it and it would make splitting the two steps more difficut (can all be fixed though)

	}
}


void design_cluster(
					const core::pose::Pose & p, // this is the native/reference structure
					core::pose::Pose repacked, // this one will be modified (well, repacked)
					core::scoring::ScoreFunctionOP score_fxn,
					std::set<core::Size> cluster,
					const std::string cluster_acc, // information about the cluster sf_acc -- not very useful at the moment though
					std::string p_name,
					utility::io::ozstream & outfile) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	const Real repack_sphere_radius = option[ OptionKeys::design_tight_clusters::repack_radius ]();

	utility::vector1<core::Size> positions;
	utility::vector1<std::string> cluster_pos_in_PDB_numbering; // for use by external packers & evaluation
	utility::vector1<std::string> repacked_pos_in_PDB_numbering;
	// perform actual repacking, or just wrap around RMSD and chi diff calculations?

	// set up packer task
	core::pack::task::PackerTaskOP base_packer_task( core::pack::task::TaskFactory::create_packer_task( p ));
	// we now need to tell the TaskFactory which parts are allowed to move and/or be designed

	base_packer_task->set_bump_check( true );
	base_packer_task->initialize_from_command_line();
	// base_packer_task->or_include_current( false ); // we'll want to allow include_current on all but the initial packer_task, so if they're derived from this one we need to set or_include_current( false ) after the split

	core::pack::task::PackerTaskOP initial_packer_task( base_packer_task->clone() );
	initial_packer_task->or_include_current( false ); // the initial repack trial must not include the current (native!) side chains
	core::pack::task::PackerTaskOP repack_packer_task( base_packer_task->clone() );
	repack_packer_task->or_include_current( true ); // all later repacking steps may make use of previous results
	utility::vector1<bool> allow_repacked( p.total_residue(), false );
	utility::vector1<bool> allow_design( p.total_residue(), false );

	utility::vector1<bool> ala_aalist( chemical::num_canonical_aas, false );
	ala_aalist[ chemical::aa_ala ] = true;


	// now we first need to iterate over all cluster positions once, to make sure we know about all that should be designed and hence NOT restricted to repacking
	for(std::set<core::Size>::iterator iter = cluster.begin(); iter != cluster.end(); iter++) {
		if (!allow_design.at(*iter)) { // prevent duplicates -- shouldn't matter here as only the cluster residues are allowed to be designed though
			allow_design.at((*iter)) = true;
			allow_repacked.at((*iter)) = true; // don't think we need this, design should be "above" packing if there is a hierarchy
			positions.push_back(*iter);
			//repack_packer_task->nonconst_residue_task((*iter)).restrict_to_repacking();
			//TR << *iter << " " << p.pdb_info()->chain(*iter) << p.pdb_info()->number(*iter) << std::endl;

			repacked_pos_in_PDB_numbering.push_back(utility::to_string(p.pdb_info()->chain(*iter)) + ":" + utility::to_string(p.pdb_info()->number(*iter)));
		}

		cluster_pos_in_PDB_numbering.push_back(utility::to_string(p.pdb_info()->chain(*iter)) + ":" + utility::to_string(p.pdb_info()->number(*iter)));

	} // first iteration done -- now we know which residues may be designed


	for(std::set<core::Size>::iterator iter = cluster.begin(); iter != cluster.end(); iter++) { // determine repack shell
		std::map<core::Size, bool> neighbor_map;
		detect_neighbors(p, *iter, repack_sphere_radius, neighbor_map);
		for (std::map<core::Size, bool>::const_iterator m_iter = neighbor_map.begin(); m_iter != neighbor_map.end(); m_iter++) {
			if (m_iter->second && !allow_repacked.at(m_iter->first)) { // try to avoid duplicates
				allow_repacked.at((m_iter->first)) = true;
				positions.push_back(m_iter->first);
				if (!allow_design.at(m_iter->first)) { // we must not restrict the cluster residues
				  initial_packer_task->nonconst_residue_task(m_iter->first).restrict_to_repacking(); // residues in the shell will only be repacked, but not designed
				  repack_packer_task->nonconst_residue_task(m_iter->first).restrict_to_repacking();
				}
				repacked_pos_in_PDB_numbering.push_back(utility::to_string(p.pdb_info()->chain(m_iter->first)) + ":" + utility::to_string(p.pdb_info()->number(m_iter->first)));

			}
		}
	}

	initial_packer_task->restrict_to_residues( allow_repacked );
	repack_packer_task->restrict_to_residues( allow_repacked );

	if ( option[ design_tight_clusters::debug ]() ) {
		TR << "repacked positions: ";
		for (core::Size pos = 1; pos <= positions.size(); pos++) {
			TR << positions[pos] << ", ";
		}
		TR << std::endl;
		TR << *repack_packer_task << std::endl; // AS_DEBUG

	}

	// favor native residues?
	if ( option[ OptionKeys::design_tight_clusters::favor_native ].user() ) {
		protocols::protein_interface_design::FavorNativeResidue fnr(repacked, option[ OptionKeys::design_tight_clusters::favor_native ]());

		// this is a type of constraint -- and apparently we need to set the constraint weight specifically
		if( score_fxn->get_weight( core::scoring::res_type_constraint ) == 0.0 ){
			//	TR << " -- constraint weight was 0 so far... " << std::endl; // AS_DEBUG
			score_fxn->set_weight( core::scoring::res_type_constraint, option[ OptionKeys::design_tight_clusters::favor_native ]() );
			//TR << *(repacked.constraint_set()) << std::endl;
		}
	}

	//TR << repacked << std::endl;


	// making cartmin available to all protocols
	bool const cartmin( option[ design_tight_clusters::cartmin ] ); /// option: use cartesian-space min
	// Phil: you should also allow setting of just the cart_bonded subterms (angle/length/etc)
	if ( cartmin ) runtime_assert( score_fxn->get_weight( core::scoring::cart_bonded ) > 1e-3 );

	// initial repack & design step
	core::pack::pack_rotamers( repacked, *score_fxn, initial_packer_task );

	for (int c = 1; c <= option[ OptionKeys::design_tight_clusters::num_cycles ](); c++) {

	  TR << " cycle " << c << std::endl;
		// first operation is always standard packing -- only affected by the command line flags, and by which rotamers get to move
		// this sets the baseline and ensures that native side chains are discarded (except when running with -use_input_sc)
		// NOTE: now this also includes a first round of design (!)
	  if (c > 1)
	    core::pack::pack_rotamers( repacked, *score_fxn, repack_packer_task );

	  // next: repack -- and possibly design, depending on the mover -- the pose
	  repack_step(repacked, repack_packer_task, score_fxn, allow_repacked);
	}

	// evaluation

	core::pose::Pose starting_cluster;
	core::pose::Pose repacked_cluster;
	utility::vector1<core::Real> individual_pos_RMSDs; // all-atom RMSDs for each cluster position individually

	std::string cluster_name = "";
	std::map<core::Size, core::Size> all_neighbors; // make sure we don't append residues multiple times to a pose, that would probably be a mess
	for(std::set<core::Size>::iterator iter = cluster.begin(); iter != cluster.end(); iter++) {
		cluster_name += utility::to_string(p.pdb_info()->number(*iter)) + "_";
		starting_cluster.append_residue_by_jump(p.residue(*iter), starting_cluster.total_residue());
		repacked_cluster.append_residue_by_jump(repacked.residue(*iter), repacked_cluster.total_residue());

		// get all-atom RMSD for this particular residue
		core::pose::Pose i_native; // empty pose, then just append the single residue - will this work?
		core::pose::Pose i_repacked;
		i_native.append_residue_by_jump(p.residue(*iter), i_native.total_residue());
		i_repacked.append_residue_by_jump(repacked.residue(*iter), i_repacked.total_residue());
		individual_pos_RMSDs.push_back(core::scoring::rmsd_with_super(i_native, i_repacked, core::scoring::is_heavyatom));


		// we'll also want the repacked residues in the output set -- note that this will affect the RMSD calculation (!)
		std::map<core::Size, bool> neighbor_map;
		detect_neighbors(p, *iter, repack_sphere_radius, neighbor_map);
		for (std::map<core::Size, bool>::const_iterator m_iter = neighbor_map.begin(); m_iter != neighbor_map.end(); m_iter++) {
			if (m_iter->second) {
				if (all_neighbors.find(m_iter->first) == all_neighbors.end()) {
					starting_cluster.append_residue_by_jump(p.residue(m_iter->first), starting_cluster.total_residue());
					repacked_cluster.append_residue_by_jump(repacked.residue(m_iter->first), repacked_cluster.total_residue());
					all_neighbors[m_iter->first]++;
				}
			}
		}
	}

	core::Real motif_starting_energy = 0.0;

	for(core::Size i = 1; i <= positions.size(); i++) {
		for(core::Size j = i+1; j <= positions.size(); j++) {
			core::scoring::EMapVector energymap;
			score_fxn->eval_ci_2b_sc_sc(p.residue(positions[i]), p.residue(positions[j]), p, energymap);
			score_fxn->eval_cd_2b_sc_sc(p.residue(positions[i]), p.residue(positions[j]), p, energymap);
			core::Real two_body_energy = score_fxn->weights().dot(energymap);
			motif_starting_energy += two_body_energy;
		}
	}

	core::Real motif_repacked_energy = 0.0;

	for(core::Size i = 1; i <= positions.size(); i++) {
		for(core::Size j = i+1; j <= positions.size(); j++) {
			core::scoring::EMapVector energymap;
			score_fxn->eval_ci_2b_sc_sc(repacked.residue(positions[i]), repacked.residue(positions[j]), repacked, energymap);
			score_fxn->eval_cd_2b_sc_sc(repacked.residue(positions[i]), repacked.residue(positions[j]), repacked, energymap);
			core::Real two_body_energy = score_fxn->weights().dot(energymap);
			motif_repacked_energy += two_body_energy;
		}
	}

	if (basic::options::option[ basic::options::OptionKeys::design_tight_clusters::generate_output_structures ]() ) {
		//starting_cluster.dump_pdb( p_name+"_"+cluster_name+"starting.pdb" );
		//repacked_cluster.dump_pdb( option[ out::path::pdb ]().path() + p_name+"_"+cluster_name+"repacked.pdb" );
		repacked_cluster.dump_pdb( p_name+"_"+cluster_name+"repacked.pdb" );// just the repacked residues
		//repacked.dump_pdb( p_name+"_"+cluster_name+"repacked.pdb" ); // full structure
	}

	core::Real motif_rmsd = core::scoring::rmsd_with_super(starting_cluster, repacked_cluster, core::scoring::is_heavyatom);

	utility::vector1<std::string> mutations;
	utility::vector1<std::string> chi_dev_details;
	compare_residues_and_chi1_2_angles(p, repacked, cluster, mutations, chi_dev_details);

	outfile << p_name << "\t"
	<< cluster_name << "\t"
	<< motif_rmsd;

	for (core::Size i = 1; i <= mutations.size(); i++) {
		outfile << "\t" << mutations[i];
	}

	// report failures per residue
	for (core::Size i = 1; i <= chi_dev_details.size(); i++) {
		outfile << "\t" << chi_dev_details[i];
	}

	// report RMSD per residue -- note that I think these are not symmetry-corrected
	for (core::Size i = 1; i <= individual_pos_RMSDs.size(); i++) {
		outfile << "\t" << individual_pos_RMSDs[i];
	}

	outfile << "\t" << cluster_acc << "\t" << p.energies().total_energy() << "\t"
	<< repacked.energies().total_energy() << "\t"
	<< motif_starting_energy << "\t"
	<< motif_repacked_energy << "\t"
	<< repacked_cluster.sequence() << "\t";

	for (core::Size i = 1; i <= cluster_pos_in_PDB_numbering.size(); i++) {
		outfile << cluster_pos_in_PDB_numbering[i] << ",";
	}

	outfile << "\t";

	for (core::Size i = 1; i <= repacked_pos_in_PDB_numbering.size(); i++) {
		outfile << repacked_pos_in_PDB_numbering[i] << ",";
	}

	outfile << std::endl;
}



int main( int argc, char * argv [] )
{
	try{

		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		// register specific options
		OPT(packing::resfile);
		OPT(in::file::native); // for comparison against externally repacked structures
		NEW_OPT(design_tight_clusters::generate_output_structures, "generate output structures (default: no)", false);
		NEW_OPT(design_tight_clusters::require_renumbered_structures, "require that input structures be renumbered sequentially, to match Rosetta numbering? Hint: necessary for compatibility with external tools that rely on this numbering. (default: no)", false);
		NEW_OPT(design_tight_clusters::debug, "debug mode? (tons of extra output)", false);
		NEW_OPT(design_tight_clusters::repack_radius, "radius around the 4 cluster residues that should be repacked (but not designed)", 4);
		NEW_OPT(design_tight_clusters::num_cycles, "number of pack/design cycles to be performed", 2);
		NEW_OPT(design_tight_clusters::favor_native, "energy benefit for native residues", 0.0);

		//	NEW_OPT(design_tight_clusters::cluster_file_suffix, "extension for the filename with the cluster results", "test");
		NEW_OPT(design_tight_clusters::chi_deviation_thresholds, "list of chi1/2 deviation thresholds", utility::vector1<core::Size>()); // -- haven't figured out how to give this default weights: 40, 20, 10));
		NEW_OPT(design_tight_clusters::min_pack, "use MinPacker", false);
		//NEW_OPT(design_tight_clusters::stochastic_pack, "use stochastic pack", false);
		NEW_OPT(design_tight_clusters::rot_trials_min, "use rotamer trials with minimization after initial repack (10 trials)", false);
		NEW_OPT(design_tight_clusters::rot_trials, "use rotamer trials after initial repack", false);
		NEW_OPT(design_tight_clusters::sidechain_fastrelax, "use sidechain fast-relax protocol", false);
		NEW_OPT(design_tight_clusters::cartmin, "Use cartesian space minimization. Only for sidechain-fast-relax protocol right now. Score function must have non-zero weight for cart_bonded term.", false);
		NEW_OPT(design_tight_clusters::backrub, "use backrub to move both backbone and side chain", false);
		NEW_OPT(design_tight_clusters::backrub_sc_prob, "probability of making a side chain move within backrub (default 0.25)", 0.25);
		NEW_OPT(design_tight_clusters::external_cluster, "provide all cluster positions in PDB numbering (chain pos vector)", utility::vector1<std::string>()); // half of these will be casted to chars, the other half to ints...
		//	NEW_OPT(design_tight_clusters::evaluate_externally_repacked_structure, "evaluate chi1/2 of an externally repacked structure -- requires -in:file:native for the reference structure", false);
		NEW_OPT(design_tight_clusters::soft_repack_min_hard_repack_min, "(hard repack to normalize, as for all runs, then) soft repack, hard minimize, hard repack w/ include_current, hard min", false);



		// initialize Rosetta
		devel::init(argc, argv);

		std::stringstream cluster_filename_suffix;
		if (basic::options::option[ basic::options::OptionKeys::design_tight_clusters::repack_radius ]()) {
			cluster_filename_suffix << "." << basic::options::option[ basic::options::OptionKeys::design_tight_clusters::repack_radius ]() << "A_repack";
		}

		/*
		 if (option[ design_tight_clusters::cluster_file_suffix ].user()) {
		 cluster_filename_suffix << "." << option[ design_tight_clusters::cluster_file_suffix ]();
		 TR << "setting cluster file suffix: " << cluster_filename_suffix << std::endl;
		 }
		 */
		//cluster_filename_suffix << "." << option[ design_tight_clusters::chi_deviation_threshold ]();
		for (IntegerVectorOption::const_iterator it = option[ design_tight_clusters::chi_deviation_thresholds ]().begin(), end = option[ design_tight_clusters::chi_deviation_thresholds ]().end(); it != end; it++) {
			cluster_filename_suffix << "." << *it;
		}

		// get scoring function
		core::scoring::ScoreFunctionOP score_fxn = core::scoring::get_score_function();

		utility::vector1<std::string> pdbs = basic::options::option[ in::file::s ]();
		core::pose::Pose p;
		std::string p_name = pdbs[1]; // warning -- this just handles one PDB...

		core::import_pose::pose_from_pdb( p, p_name );

		std::string outfile_core = p_name;


		//p.dump_pdb( p_name+"_native.pdb" ); // only for debugging the odd rebuild-neighbor-count issues



		if ( option[ out::path::pdb ].user() ) {
			size_t si = outfile_core.rfind("/"); // last occurrence of /
			if (si != std::string::npos)
				outfile_core = outfile_core.substr(si+1); // // we only want the name, but leave out the input path
			outfile_core = option[ out::path::pdb ]().path() + outfile_core;
		} else if ( option[ out::path::all ].user() ) {
			size_t si = outfile_core.rfind("/"); // last occurrence of /
			if (si != std::string::npos)
				outfile_core = outfile_core.substr(si+1); // // we only want the name, but leave out the input path
			outfile_core = option[ out::path::all ]().path() + outfile_core;
		}

		// later: add option to compare against in::file::native, in case we want to (fast?)relax structures first [ideally restricting side chain conformations, I think there was some protocol for this...] and use those as starting structures, but then compare against the actual native/crystal

		(*score_fxn)(p);
		score_fxn->show(TR, p);
		TR << std::endl;

		// read in a resfile, if specified -- with this the user may disable regions from cluster detection, e.g. those that are too close to a HETATM (which may not be visible in runs with -ignore_unrecognized_res)
		core::pack::task::TaskFactoryOP input_task_factory = new core::pack::task::TaskFactory;
		input_task_factory->push_back( new core::pack::task::operation::InitializeFromCommandline );
		if ( option[ packing::resfile ].user() ) {
			input_task_factory->push_back( new core::pack::task::operation::ReadResfile );
		} else {
			core::pack::task::operation::RestrictToRepackingOP rtrop = new core::pack::task::operation::RestrictToRepacking; // make sure that by default everything is allowed to be repacked
			input_task_factory->push_back( rtrop );
		}

		//utility::vector1<std::set<core::Size> > filtered_clusters;
		std::map< std::set<core::Size>, std::string > filtered_clusters; // also store whether the cluster is intermediate or exposed

		core::pack::task::PackerTaskOP input_task(input_task_factory->create_task_and_apply_taskoperations(p)); // this seems to be a really inefficient way do to this...

		utility::io::ozstream outfile( outfile_core + ".cluster" + cluster_filename_suffix.str()); // this should be path-adapted


		if ( option[ design_tight_clusters::debug ]() ) {
			p.dump_pdb( p_name+"_input.pdb" ); // to make sure Rosetta doesn't do anything funny to OSCAR's structure -- in particular, look at the Histidines!
		}

		pose::Pose native_p(p); // actual native if provided, otherwise this is a copy of the input structure
		(*score_fxn)(native_p);

		if ( option[ design_tight_clusters::external_cluster ]().size() % 2 != 0) {
			TR << "error -- -external_cluster parameters must always be pairs of chain and position!" << std::endl;
			return 2;
		}


		/* --- not implemented yet ---
		 if ( option [ design_tight_clusters::evaluate_externally_repacked_structure ]() ) { // requires native for reference
		 // it might be easier to directly call comparison on the given cluster
		 if (! option[ in::file::native ].user() ) {
		 TR << "error -- native pose required for comparison with externally repacked structures!" << std::endl;
		 return 2;
		 }
		 native_p.clear();
		 core::import_pose::pose_from_pdb(native_p, option[ in::file::native ]()); // does this work?
		 (*score_fxn)(native_p); // initial scoring is required to get actual energies
		 }
		 */

		std::set<core::Size> ext_cluster;
		std::string cluster_name = "";
		// also extract the actual 4 residues for getting the RMSD
		core::pose::Pose native_cluster;
		core::pose::Pose repacked_cluster;

		for (core::Size c_pos = 1; c_pos <= option[ design_tight_clusters::external_cluster ]().size(); c_pos += 2) {

			core::Size mapped_pos = p.pdb_info()->pdb2pose(option[ design_tight_clusters::external_cluster ][c_pos][0], atoi(option[ design_tight_clusters::external_cluster ][c_pos+1].c_str()));

			ext_cluster.insert(mapped_pos); // we have ensured above that there's an even number of entries, so this shouldn't produce segfaults -- the casting from strings to chars and ints is a terrible hack though
			cluster_name += utility::to_string(p.pdb_info()->number(mapped_pos)) + "_";

			native_cluster.append_residue_by_jump(native_p.residue(mapped_pos), native_cluster.total_residue());
			repacked_cluster.append_residue_by_jump(p.residue(mapped_pos), repacked_cluster.total_residue());


			if ( option[ design_tight_clusters::debug ]() ) {
				TR << " *** evaluating external cluster *** " << option[ design_tight_clusters::external_cluster ][c_pos]
				<< " " << option[ design_tight_clusters::external_cluster ][c_pos+1] << " "
				<< mapped_pos << std::endl;
			}


		}


		// append to filtered_clusters in either case -- repacking or just evaluation will be checked inside the function
		filtered_clusters[ext_cluster] = "external"; // we don't know whether this one is buried or intermediate



		// repack all the filtered clusters (whether externally provided or detected in this run)
		for (std::map<std::set<core::Size>, std::string>::const_iterator mi = filtered_clusters.begin(); mi != filtered_clusters.end(); mi++) {
			design_cluster(native_p, p, score_fxn, mi->first, mi->second, outfile_core, outfile);
		}

		outfile.close();

		//score_fxn->show(TR, p);


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
