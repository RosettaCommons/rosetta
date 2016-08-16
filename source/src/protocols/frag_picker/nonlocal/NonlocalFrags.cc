// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/frag_picker/nonlocal/NonlocalFrags.cc
/// @author David Kim (dekim@u.washington.edu)

#ifdef BOINC
#include <protocols/boinc/boinc.hh>
#endif


// Unit headers
#include <protocols/frag_picker/nonlocal/NonlocalFrags.hh>

// C/C++ headers

// Package headers

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/RelaxProtocolBase.hh>
#include <protocols/relax/util.hh>
#include <protocols/simple_filters/DdgFilter.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/Job.fwd.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/JobOutputter.fwd.hh>
#include <protocols/idealize/IdealizeMover.hh>
#include <protocols/simple_filters/RmsdEvaluator.hh>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.hh>
#include <utility/exit.hh>
#include <utility/file/file_sys_util.hh>
#include <numeric/random/random.hh>

#include <sstream>
#include <fstream>


namespace protocols {
namespace frag_picker {
namespace nonlocal {

using namespace core;

static THREAD_LOCAL basic::Tracer TR( "protocols.frag_picker.nonlocal.NonlocalFrags" );

void NonlocalFrags::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// input
	option.add_relevant( in::file::s );   // PDB
	option.add_relevant( in::file::l );   // PDB list
	option.add_relevant( in::file::list ); // PDB list
	option.add_relevant( frags::nonlocal::single_chain ); // pairs from same chain only
	// relax options
	option.add_relevant( frags::nonlocal::relax_input );
	option.add_relevant( frags::nonlocal::relax_input_with_coordinate_constraints );
	option.add_relevant( frags::nonlocal::relax_frags_repeats );
	option.add_relevant( frags::contacts::min_seq_sep );
	option.add_relevant( frags::contacts::dist_cutoffs );
	// non-local contact definition options
	option.add_relevant( frags::nonlocal::min_contacts_per_res );
	option.add_relevant( frags::nonlocal::max_rmsd_after_relax );

	// output
	option.add_relevant( frags::nonlocal::output_idealized );
	option.add_relevant( frags::nonlocal::output_frags_pdbs );
	// fragment sizes to output
	option.add_relevant( frags::frag_sizes );
}

NonlocalFrags::NonlocalFrags() :
	single_chain_( false ),
	relax_input_( false ),
	relax_input_with_coordinate_constraints_( false ),
	relax_frags_repeats_( 1 ),
	checkpointfile_( "nonlocalfrags.checkpoint" ),
	min_seq_sep_( 12 ),
	ca_dist_squared_( 100.0 ),
	min_contacts_per_res_( 1 ),
	max_rmsd_after_relax_( 1.5 ),
	max_ddg_score_( -4.0 ),
	output_frags_pdbs_( false ),
	output_idealized_( false )
{
	initialize();
}

void NonlocalFrags::initialize() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( !option[ in::file::s ].user() && !option[ in::file::l ].user() && !option[ in::file::list ].user() ) {
		utility_exit_with_message( "Error: in:file:s or in:file:l option required!" );
	}
	if ( option[ frags::nonlocal::single_chain ].user() ) {
		single_chain_ = option[ frags::nonlocal::single_chain ]();
	}
	if ( option[ frags::nonlocal::relax_input ].user() ) {
		relax_input_ = option[ frags::nonlocal::relax_input ]();
	}
	if ( option[ frags::nonlocal::relax_input_with_coordinate_constraints ].user() ) {
		relax_input_with_coordinate_constraints_ = option[ frags::nonlocal::relax_input_with_coordinate_constraints ]();
	}
	if ( option[ frags::nonlocal::relax_frags_repeats ].user() ) {
		relax_frags_repeats_ = option[ frags::nonlocal::relax_frags_repeats ]();
	}
	// non-local contact definition
	if ( option[ frags::contacts::min_seq_sep ].user() ) {
		min_seq_sep_ = option[ frags::contacts::min_seq_sep ]();
	}
	if ( option[ frags::contacts::dist_cutoffs ].user() ) {
		utility::vector1<Real> dist_cutoffs = option[ frags::contacts::dist_cutoffs ]();
		if ( dist_cutoffs.size() != 1 ) {
			utility_exit_with_message( "Error: only one frags::contacts::dist_cutoffs value is allowed!" );
		}
		Real min_dist = dist_cutoffs[1];
		ca_dist_squared_ = min_dist*min_dist;
	}
	if ( option[ frags::nonlocal::min_contacts_per_res ].user() ) {
		min_contacts_per_res_ = (Size)option[ frags::nonlocal::min_contacts_per_res ]();
	}
	// for valid interacting pair
	if ( option[ frags::nonlocal::max_ddg_score ].user() ) {
		max_ddg_score_ = option[ frags::nonlocal::max_ddg_score ]();
	}
	if ( option[ frags::nonlocal::max_rmsd_after_relax ].user() ) {
		max_rmsd_after_relax_ = option[ frags::nonlocal::max_rmsd_after_relax ]();
	}
	// fragment sizes
	if ( option[ frags::frag_sizes ].user() ) {
		frag_sizes_ = option[ frags::frag_sizes ]();
	} else {
		frag_sizes_.clear(); // default to no frag sizes
	}
	if ( option[ frags::nonlocal::output_idealized ].user() ) {
		output_idealized_ = option[ frags::nonlocal::output_idealized ]();
	}

	if ( option[ frags::nonlocal::output_frags_pdbs ].user() ) {
		output_frags_pdbs_ = option[ frags::nonlocal::output_frags_pdbs ]();
	}

	// make sure the silent file output is binary
	option[ out::file::silent_struct_type ].value("binary");

	read_checkpoint_file();
}

bool NonlocalFrags::recover_checkpoint( std::string const & tag, pose::Pose& pose ) {
	std::map<std::string, std::string>::iterator it;
	it=checkpoints_map_.find(tag);
	if ( it != checkpoints_map_.end() ) {
		// add recovered data to pose as a comment for output
		if ( checkpoints_map_[ tag ].size() ) pose::add_comment( pose, tag, checkpoints_map_[ tag ] );
		TR.Info << "Recovered checkpoint " << tag << std::endl;
		return true;
	}
	return false;
}

void NonlocalFrags::write_checkpoint( std::string const & tag, std::string const & data ) {
	std::ofstream checkpoint( checkpointfile_.c_str(), std::ios_base::app );
	if ( !checkpoint.good() ) {
		TR.Debug << "Warning: cannot open checkpoint file for writing!" << std::endl;
		return;
	}
	// check data

#ifdef BOINC
	boinc_begin_critical_section();
#endif
	if ( data.size() ) {
		Size pdbinfonumberj, pdbinfonumberk, contact_cnt;
		Real sub_pose_score, relaxed_rmsd, relaxed_score, relaxed_ddg_score;
		std::string chainj, chaink;
		std::istringstream line_stream(data);
		line_stream >> pdbinfonumberj >> chainj >> pdbinfonumberk >> chaink >>
			contact_cnt >> sub_pose_score >> relaxed_rmsd >> relaxed_score >> relaxed_ddg_score;
		if ( line_stream.fail() ) {
			TR.Debug << "Warning: cannot parse checkpoint data!" << std::endl;
			return;
		}
		checkpoint << tag << " " << data << std::endl;
		checkpoints_map_[ tag ] = data;
	} else {
		checkpoint << tag << std::endl;
		checkpoints_map_[ tag ] = "";
	}
#ifdef BOINC
	boinc_end_critical_section();
#endif
	TR.Info << "Write checkpoint " << tag << std::endl;
}

void NonlocalFrags::read_checkpoint_file() {
	if ( !utility::file::file_exists( checkpointfile_ ) ) return;
	std::string line;
	std::ifstream checkpoint( checkpointfile_.c_str() );
	while ( getline(checkpoint,line) ) {
		Size pdbinfonumberj, pdbinfonumberk, contact_cnt;
		Real sub_pose_score, relaxed_rmsd, relaxed_score, relaxed_ddg_score;
		std::string tag, chainj, chaink;

		std::istringstream line_stream(line);
		line_stream >> tag;
		if ( line_stream.eof() ) {
			checkpoints_map_[ tag ] = "";
		} else {
			line_stream >> pdbinfonumberj >> chainj >> pdbinfonumberk >> chaink >>
				contact_cnt >> sub_pose_score >> relaxed_rmsd >> relaxed_score >> relaxed_ddg_score;
			if ( line_stream.fail() ) utility_exit_with_message( "Error: checkpoint file parse error!" );
			std::stringstream data;
			data << pdbinfonumberj << " " << chainj << " " << pdbinfonumberk << " " << chaink <<
				" " << contact_cnt << " " << sub_pose_score << " " << relaxed_rmsd << " " << relaxed_score <<
				" " << relaxed_ddg_score;
			checkpoints_map_[ tag ] = data.str();
		}
		TR.Info << "Read checkpoint " << tag << std::endl;
	}
	checkpoint.close();
}

void NonlocalFrags::apply(pose::Pose& pose) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	pose::Pose unmodified_pose = pose;
	pose::Pose relaxed_pose;

	// set input pose as native
	// Should not use in:file:native because the minirosetta app will add it to BOINC graphics and we will be using sub-poses for the fragment pairs
	protocols::jd2::JobDistributor* jd = protocols::jd2::JobDistributor::get_instance();
	// clear the evaluators or else different input pdbs may cause a runtime error
	// when doing rmsd evaluations
	jd->job_outputter()->clear_evaluators();
	jd->job_outputter()->add_evaluation( evaluation::PoseEvaluatorOP( new simple_filters::RmsdEvaluator( pose::PoseOP( new pose::Pose( unmodified_pose ) ), "" ) ) );

	//scoring::ScoreFunctionOP scorefxn = scoring::get_score_function();
	scoring::ScoreFunctionOP scorefxn = scoring::get_score_function();
	Real unmodified_pose_score = (*scorefxn)( pose );

	std::string output_name = protocols::jd2::current_output_name();

#ifdef BOINC_GRAPHICS
	Size step_cnt = 0;
	std::string input_tag = protocols::jd2::get_current_job()->input_tag();
	protocols::boinc::Boinc::attach_graphics_current_pose_observer( pose );
#endif

	// relax input pose before picking non-local fragment pairs?
	if ( frag_sizes_.size() > 0 && (relax_input_ || relax_input_with_coordinate_constraints_) ) {

#ifdef BOINC_GRAPHICS
		// for some reason this is necessary to display the structures initially, i.e. the relax protocol may not
		protocols::boinc::Boinc::update_graphics_low_energy( pose, unmodified_pose_score );
		protocols::boinc::Boinc::update_graphics_last_accepted( pose, unmodified_pose_score );
		protocols::boinc::Boinc::update_mc_trial_info( step_cnt++, "Relax_" + input_tag );
#endif
		if ( relax_input_with_coordinate_constraints_ ) {
			// a hack to turn on coordinate constraints (we don't want to use them for relaxing the fragment pairs later)
			option[OptionKeys::relax::constrain_relax_to_start_coords].value(true);
		}

		// a hack to checkpoint this relax
		bool origval = option[run::delete_checkpoints]();
		option[run::delete_checkpoints].value(false);
		protocols::relax::RelaxProtocolBaseOP relax_protocol = protocols::relax::generate_relax_from_cmd();
		relax_protocol->set_current_tag( output_name );
		relax_protocol->apply( pose );
		// restore original option value
		option[run::delete_checkpoints].value(origval);

		if ( relax_input_with_coordinate_constraints_ ) {
			// turn off coordinate constraints
			option[OptionKeys::relax::constrain_relax_to_start_coords].value(false);
			scorefxn->set_weight(scoring::coordinate_constraint, 0.0 );
			scorefxn->set_weight(scoring::atom_pair_constraint, 0.0 );
			scorefxn->set_weight(scoring::angle_constraint, 0.0 );
			scorefxn->set_weight(scoring::dihedral_constraint, 0.0 );
			scorefxn->set_weight(scoring::res_type_constraint, 0.0);
		}
		// remove virtual residue
		pose::Pose temp_pose = unmodified_pose;
		for ( Size ii = 1; ii <= unmodified_pose.total_residue(); ++ii ) {
			temp_pose.replace_residue( ii, pose.residue( ii ), false );
		}
		pose = temp_pose;
		relaxed_pose = pose;
	}

	// add input score comment for output
	std::stringstream scorestr;
	scorestr << unmodified_pose_score;
	pose::add_comment( pose, output_name + "_input_score", scorestr.str() );

	// setup relax protocol for sub pose (frag pair)
	protocols::relax::RelaxProtocolBaseOP sub_pose_relax_protocol( new protocols::relax::FastRelax( scorefxn, relax_frags_repeats_ ) );
	kinematics::MoveMapOP mm = sub_pose_relax_protocol->get_movemap();
	mm->set_jump(true); // set jumps movable

	// iterate each frag size
	// note: shall this be random for BOINC?
	Size total_residue = pose.total_residue();

	// Size frag_len = numeric::random::random_element( frag_sizes_ );  //  = frag_sizes_[ i ]; // frag size
	for ( Size i=1; i <= frag_sizes_.size(); ++i ) {   // frag_sizes_ loop

		Size frag_len = frag_sizes_[ i ]; // frag size


		Size total_len = frag_len*2;
		if ( total_residue >= frag_len+min_seq_sep_+frag_len ) {    // continue; // skip if pose is too small

			Size min_contacts = min_contacts_per_res_*frag_len; // min non-local contacts

#ifdef BOINC_GRAPHICS
			pose::Pose best_energy_pose;
			Real min_relaxed_score = 100000000.0;
#endif

			// find non-local pairs
			Size contact_cnt = 0;
			bool continue_k = false;
			for ( Size j=1; j<=total_residue-frag_len-min_seq_sep_; ++j ) {
				for ( Size k=j+frag_len+min_seq_sep_; k<=total_residue-frag_len; ++k ) {
					if ( single_chain_ && (pose.residue(j).chain() != pose.residue(k).chain()) ) continue;

					// get tag for this frag pair
					std::stringstream sstream;
					sstream << "_" << frag_len << "_" << j << "_" << k;
					std::string tag = sstream.str();

					// recover checkpoint
					if ( recover_checkpoint( output_name + tag, pose ) ) {
#ifdef BOINC_GRAPHICS
						protocols::boinc::Boinc::update_mc_trial_info( step_cnt++, "NonLocalFrags_" + input_tag + tag );
#endif
						continue;
					}

					pose::PDBInfoCOP pdbinfo = pose.pdb_info();

					// j and k are starting positions for fragment pairs
					// count CA contacts
					contact_cnt = 0;
					continue_k = false;
					Size prev_chain_m = 0;
					Size prev_chain_n = 0;
					Size prev_resnum_m = 0;
					Size prev_resnum_n = 0;
					numeric::xyzVector< Real> prev_ca_xyz_m(0.0), prev_ca_xyz_n(0.0);
					for ( Size m=0; m<frag_len; ++m ) {
						conformation::Residue const & rsd_m = pose.residue(j+m);
						if ( !rsd_m.has("CA") || !rsd_m.is_protein() ) { continue_k = true; break; }
						Size chain_m = rsd_m.chain();
						Size resnum_m = pdbinfo->number(j+m);
						numeric::xyzVector< Real> ca_xyz_m = rsd_m.xyz("CA");
						// make sure the fragment residues are sequential in the orig PDB and from the same chain!
						if ( m>0 && (chain_m != prev_chain_m || resnum_m != prev_resnum_m+1 ||
								ca_xyz_m.distance_squared(prev_ca_xyz_m) > 25) ) { // CA-CA distance must be within 5 angstroms (i.e. no chainbreak)
							continue_k = true; break;
						}
						prev_chain_m = chain_m;
						prev_resnum_m = resnum_m;
						prev_ca_xyz_m = ca_xyz_m;
						for ( Size n=0; n<frag_len; ++n ) {
							conformation::Residue const & rsd_n = pose.residue(k+n);
							Size chain_n = rsd_n.chain();
							if ( !rsd_n.has("CA") || !rsd_n.is_protein() || (single_chain_ && chain_m != chain_n) ) { continue_k = true; break; }
							Size resnum_n = pdbinfo->number(k+n);
							numeric::xyzVector< Real> ca_xyz_n = rsd_n.xyz("CA");
							// make sure the fragment residues are sequential in the orig PDB and from the same chain!
							if ( n>0 && (chain_n != prev_chain_n || resnum_n != prev_resnum_n+1 ||
									ca_xyz_n.distance_squared(prev_ca_xyz_n) > 25) ) {  // CA-CA distance must be within 5 angstroms (i.e. no chainbreak)
								continue_k = true; break;
							}
							prev_chain_n = chain_n;
							prev_resnum_n = resnum_n;
							prev_ca_xyz_n = ca_xyz_n;
							if ( rsd_m.xyz("CA").distance_squared(rsd_n.xyz("CA")) < ca_dist_squared_ ) contact_cnt++;
						}
						if ( continue_k ) break;
					}

					if ( !continue_k && contact_cnt >= min_contacts ) {
						// At this point, fragment pair has enough non-local contacts, no missing CA, and is a protein.

#ifdef BOINC_GRAPHICS
						protocols::boinc::Boinc::update_mc_trial_info( step_cnt++, "NonLocalFrags_" + input_tag + tag );
#endif
						// Create a sub pose representing the pair
						Size midpoint = static_cast<Size>(ceil(frag_len / 2.0));
						kinematics::FoldTree tree(total_len);
						int jump_id = tree.new_jump( midpoint, frag_len+midpoint, frag_len );
						tree.set_jump_atoms(jump_id, "CA", "CA");
						assert(tree.check_fold_tree());
						//   Create the sub pose
						pose::Pose sub_pose;
						utility::vector1<Size> positions;
						for ( Size m=0; m<frag_len; ++m ) positions.push_back( j+m );
						for ( Size m=0; m<frag_len; ++m ) positions.push_back( k+m );
						pose::create_subpose( pose, positions, tree, sub_pose );

						//   Detect disulfides
						sub_pose.conformation().detect_disulfides();

						//   Treat the fragments as separate chains
						sub_pose.conformation().insert_chain_ending( frag_len );

						// RELAX non-local fragment pair
						pose::Pose relax_sub_pose = sub_pose;
						Real sub_pose_score = (*scorefxn)( sub_pose );
#ifdef BOINC_GRAPHICS
						protocols::boinc::Boinc::attach_graphics_current_pose_observer( relax_sub_pose );
						protocols::boinc::Boinc::update_graphics_low_energy( relax_sub_pose, sub_pose_score );
#endif
						sub_pose_relax_protocol->set_current_tag( output_name + tag );
						sub_pose_relax_protocol->apply( relax_sub_pose );

						// check relaxed pose
						Real relaxed_rmsd = scoring::CA_rmsd( relax_sub_pose, sub_pose );
						if ( relaxed_rmsd > max_rmsd_after_relax_ ) {
							write_checkpoint( output_name + tag, "" );
							continue;
						}
						Real relaxed_score = (*scorefxn)( relax_sub_pose );

						// DDG of relaxed fragment pair
						protocols::simple_filters::DdgFilter ddg = protocols::simple_filters::DdgFilter( 1000, scorefxn, 1, 5);
						Real relaxed_ddg_score = ddg.compute( relax_sub_pose );
						if ( relaxed_ddg_score >= max_ddg_score_ ) {  // DDG FILTER VALUE
							write_checkpoint( output_name + tag, "" );
							continue;
						}

#ifdef BOINC_GRAPHICS
						if ( relaxed_score < min_relaxed_score ) {
							min_relaxed_score = relaxed_score;
							best_energy_pose = relax_sub_pose;
							protocols::boinc::Boinc::update_graphics_low_energy( best_energy_pose, min_relaxed_score, protocols::boinc::Boinc::PERSIST );
						}
#endif

						// add comment in pose for output, this is the main data for non-local fragment pairs
						// make sure this data format matches the checkpoint methods
						//seqj = sub_pose.sequence().substr(0,frag_len)
						//seqk = sub_pose.sequence().substr(frag_len)
						std::stringstream pair_data;
						pair_data << pdbinfo->number(j) << " " << pdbinfo->chain(j) << " " << pdbinfo->number(k) << " " << pdbinfo->chain(k) <<
							" " << contact_cnt << " " << sub_pose_score << " " << relaxed_rmsd << " " << relaxed_score << " " << relaxed_ddg_score;
						pose::add_comment( pose, output_name + tag, pair_data.str() );
						// should we also use this info?
						//  pose.pdb_info()->occupancy( seqpos, ii );
						//  pose.pdb_info()->temperature( seqpos, ii );

						// DEBUG
						TR.Debug << output_name << tag << " " << pair_data.str() << std::endl;
						if ( output_frags_pdbs_ ) sub_pose.dump_pdb( output_name + tag + ".pdb" );

						write_checkpoint( output_name + tag, pair_data.str() );

#ifdef BOINC_GRAPHICS
						protocols::boinc::Boinc::update_graphics_last_accepted( relax_sub_pose, relaxed_score );
#endif
					}
					//else {  // don't bother checkpointing at this point
					// write_checkpoint( output_name + tag, "" );
					//}
				}
			}

		}   // frag_sizes_ loop


	}

	if ( output_idealized_ ) {

		pose::Pose idealize_pose = unmodified_pose;

#ifdef BOINC_GRAPHICS
		Real score = (*scorefxn)(idealize_pose);
		protocols::boinc::Boinc::attach_graphics_current_pose_observer( idealize_pose );
		protocols::boinc::Boinc::update_graphics_last_accepted( idealize_pose, score );
		protocols::boinc::Boinc::update_graphics_low_energy( idealize_pose, score, protocols::boinc::Boinc::RESET );
		protocols::boinc::Boinc::update_mc_trial_info( step_cnt++, "Idealize_" + input_tag );
#endif

		// idealize
		TR.Info << "Idealize pose" << std::endl;
		protocols::idealize::IdealizeMover idealizer;
		idealizer.fast( false );
		idealizer.chainbreaks( true );
		idealizer.apply( idealize_pose );
		Real idealize_score = (*scorefxn)(idealize_pose);
		// add input score comment for output
		std::stringstream scorestr;
		scorestr << idealize_score;
		pose::add_comment( pose, output_name + "_idealize_score", scorestr.str() );

#ifdef BOINC_GRAPHICS
		protocols::boinc::Boinc::update_mc_trial_info( step_cnt++, "Relax_" + input_tag );
#endif

		// relax idealized
		if ( relax_input_with_coordinate_constraints_ ) {
			// a hack to turn on coordinate constraints (we don't want to use them for relaxing the fragment pairs later)
			option[OptionKeys::relax::constrain_relax_to_start_coords].value(true);
		}

		// a hack to checkpoint this relax
		bool origval = option[run::delete_checkpoints]();
		option[run::delete_checkpoints].value(false);
		protocols::relax::RelaxProtocolBaseOP relax_protocol = protocols::relax::generate_relax_from_cmd();
		relax_protocol->set_current_tag( output_name + "_ideal" );

		relax_protocol->apply( idealize_pose );
		// restore original option value
		option[run::delete_checkpoints].value(origval);

		if ( relax_input_with_coordinate_constraints_ ) {
			// turn off coordinate constraints
			option[OptionKeys::relax::constrain_relax_to_start_coords].value(false);
			scorefxn->set_weight(scoring::coordinate_constraint, 0.0 );
			scorefxn->set_weight(scoring::atom_pair_constraint, 0.0 );
			scorefxn->set_weight(scoring::angle_constraint, 0.0 );
			scorefxn->set_weight(scoring::dihedral_constraint, 0.0 );
			scorefxn->set_weight(scoring::res_type_constraint, 0.0);
		}
		// remove virtual residue
		pose::Pose temp_pose = pose;
		for ( Size ii = 1; ii <= unmodified_pose.total_residue(); ++ii ) {
			temp_pose.replace_residue( ii, idealize_pose.residue( ii ), false );
		}
		pose = temp_pose;

		(*scorefxn)(pose);

	}

#ifdef BOINC
	boinc_begin_critical_section();  // for output safety
#endif

	// output relaxed pose?
	if ( relax_input_ || relax_input_with_coordinate_constraints_ ) {
		// first check if it has already been done (from a preempted run)
		// a hack since JobOutputter::affixed_numbered_name( JobCOP job ) is used to create the silent file
		// output name and uses the out::prefix option
		std::string origprefix = option[ out::prefix ].value();
		option[ out::prefix ].value("relaxed_");
		if ( !jd->job_outputter()->job_has_completed( jd->current_job() ) && relaxed_pose.total_residue() ) {
			(*scorefxn)(relaxed_pose);
			jd->job_outputter()->final_pose( jd->current_job(), relaxed_pose );
		}
		option[ out::prefix ].value(origprefix);
	}

}

std::string NonlocalFrags::get_name() const {
	return "NonlocalFrags";
}

protocols::moves::MoverOP NonlocalFrags::clone() const {
	return protocols::moves::MoverOP( new NonlocalFrags(*this) );
}

protocols::moves::MoverOP NonlocalFrags::fresh_instance() const {
	return protocols::moves::MoverOP( new NonlocalFrags() );
}


}  // namespace nonlocal
}  // namespace frag_picker
}  // namespace protocols
