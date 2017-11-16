// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SewingHasher.cc
///
/// @brief An MPI-enabled application that reads in a Hasher hashtable, scores each model in that table against all others, and
/// generates a SewGraph file.
///
/// @author Tim Jacobs

//Package headers
#include <devel/init.hh>
#include <protocols/sewing/hashing/Hasher.hh>
#include <protocols/sewing/conformation/Model.hh>
#include <protocols/sewing/util/io.hh>

//Protocol headers
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>

#include <core/pack/task/operation/TaskOperations.hh>

#include <core/conformation/util.hh>

#include <protocols/relax/FastRelax.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>

#include <utility/io/ozstream.hh>

//Utility headers
#include <utility/io/izstream.hh>
#include <utility/mpi_util.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/file/FileName.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <basic/database/sql_utils.hh>

//C++ headers
#include <map>
#include <set>
#include <iterator>

static basic::Tracer TR("SewingHasher");

int
main( int argc, char * argv [] ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::sewing;

	// initialize core and read options
	devel::init(argc, argv);
	core::Size min_score = option[sewing::min_hash_score].def(10);
	core::Size max_clash_score = option[sewing::max_clash_score].def(0);

	//////////////////////////// MODEL GENERATION ///////////////////////////////////
	std::map< int, Model > models;
	utility::vector1<utility::file::FileName> pdb_library;
	if ( option[ in::file::l ].user() ) {
		utility::vector1<utility::file::FileName> input_lists( option[ in::file::l ]() );
		for ( core::Size i = 1; i <= input_lists.size(); i++ ) {
			utility::io::izstream current_input_list( input_lists[i] );
			if ( !current_input_list.good() ) {
				utility_exit_with_message("unable to open input file file: "+input_lists[i].name()+"\n");
			}
			while ( current_input_list.good() ) {
				std::string name;
				current_input_list.getline(name);
				if ( current_input_list.good() ) pdb_library.push_back( utility::file::FileName(name) );
			}
		}

		for ( core::Size i=1; i<=pdb_library.size(); ++i ) {
			core::pose::Pose pose;
			core::import_pose::pose_from_file(pose, pdb_library[i], core::import_pose::PDB_file);
			if ( pose.size() <= 1 ) { continue; }
			utility::vector1< std::pair<core::Size,core::Size> > segments;
			segments.push_back(std::make_pair(1, pose.size()));
			Model pdb_model = create_model_from_pose(pose, segments, (int)i);
			models.insert(std::make_pair(i, pdb_model));
		}
	} else {
		utility_exit_with_message("Must specify -l flag for motif database");
	}

	Model query_model;
	if ( !option[ in::file::s ].user() ) {
		utility_exit_with_message("Must specify -s flag for query structure");
	}
	utility::file::FileName query_pdb = option[ in::file::s ]()[1];
	core::pose::Pose query_pose;
	core::import_pose::pose_from_file(query_pose, query_pdb, core::import_pose::PDB_file);
	utility::vector1< std::pair<core::Size,core::Size> > segments;
	segments.push_back(std::make_pair(1, query_pose.size()));
	query_model = create_model_from_pose(query_pose, segments, -1);

	std::map< int, Model >::const_iterator it = models.begin();
	std::map< int, Model >::const_iterator it_end = models.end();
	Hasher hasher;
	for ( ; it != it_end; ++it ) {
		hasher.insert(it->second);
	}
	ScoreResults group_scores = hasher.score(query_model, 1, min_score, max_clash_score, true);
	std::cout << "Number of hits with score >=" << min_score << " " << group_scores.size() << std::endl;

	ScoreResults::const_iterator score_it = group_scores.begin();
	ScoreResults::const_iterator score_it_end = group_scores.end();
	int counter = 0;
	for ( ; score_it != score_it_end; ++score_it ) {
		++counter;
		BasisPair bp = score_it->first;
		core::Size query_resnum = bp.first.resnum;
		int hit_model_id = bp.second.model_id;
		core::Size hit_resnum = bp.second.resnum;

		//Check to see if we have hits for every residue in the model
		core::Size n_residues = models[hit_model_id].segments_[1].residues_.size();
		AtomMap matched_atom_map = score_it->second.segment_matches.begin()->second;
		AtomMap::const_iterator atom_it = matched_atom_map.begin();
		AtomMap::const_iterator atom_it_end = matched_atom_map.end();
		std::set<core::Size> resnums;
		for ( ; atom_it != atom_it_end; ++atom_it ) {
			resnums.insert(atom_it->second.rsd());
		}
		core::Size matched_residues = resnums.size();
		if ( matched_residues < n_residues-1 ) { continue; }

		core::pose::Pose hit_pose;
		core::import_pose::pose_from_file(hit_pose, pdb_library[hit_model_id], core::import_pose::PDB_file);
		core::conformation::remove_upper_terminus_type_from_conformation_residue(hit_pose.conformation(), 1);
		core::conformation::remove_lower_terminus_type_from_conformation_residue(hit_pose.conformation(), 1);
		core::conformation::remove_upper_terminus_type_from_conformation_residue(hit_pose.conformation(), hit_pose.size());
		core::conformation::remove_lower_terminus_type_from_conformation_residue(hit_pose.conformation(), hit_pose.size());

		runtime_assert(score_it->second.segment_matches.size() == 1);

		core::id::AtomID_Map< core::id::AtomID > atom_map( core::id::AtomID::BOGUS_ATOM_ID() );
		core::pose::initialize_atomid_map( atom_map, hit_pose, core::id::AtomID::BOGUS_ATOM_ID() );

		//First just superimpose the two matching residues that make up the frame
		atom_map.set( core::id::AtomID(hit_pose.residue( hit_resnum ).atom_index("CA"), hit_resnum ), core::id::AtomID( query_pose.residue(query_resnum).atom_index("CA"), query_resnum ) );
		atom_map.set( core::id::AtomID(hit_pose.residue( hit_resnum ).atom_index("N"), hit_resnum ), core::id::AtomID( query_pose.residue(query_resnum).atom_index("N"), query_resnum ) );
		atom_map.set( core::id::AtomID(hit_pose.residue( hit_resnum ).atom_index("C"), hit_resnum ), core::id::AtomID( query_pose.residue(query_resnum).atom_index("C"), query_resnum ) );
		atom_map.set( core::id::AtomID(hit_pose.residue( hit_resnum ).atom_index("O"), hit_resnum ), core::id::AtomID( query_pose.residue(query_resnum).atom_index("O"), query_resnum ) );

		TR << "Hit " << counter << " mapping (" << hit_model_id << "):" << std::endl;
		TR << query_resnum << " " << query_pose.residue(query_resnum).name3() << " -> " << hit_resnum << " " << hit_pose.residue(hit_resnum).name3() << std::endl;

		//Then superimpose all mapped atoms
		AtomMap atom_mapping = score_it->second.segment_matches.begin()->second;
		AtomMap::const_iterator it = atom_mapping.begin();
		AtomMap::const_iterator it_end = atom_mapping.end();
		std::map<core::Size, core::Size> residue_map;
		for ( ; it != it_end; ++it ) {
			TR << it->first.rsd() << "-" << it->first.atomno() << " : " << it->second.rsd() << "-" << it->second.atomno() << std::endl;
			atom_map.set( it->second, it->first );
			residue_map.insert(std::make_pair(it->first.rsd(), it->second.rsd()));
		}
		core::scoring::superimpose_pose( hit_pose, query_pose, atom_map );
		hit_pose.dump_pdb("hit_pose_"+utility::to_string(counter)+"_"+utility::to_string(hit_model_id)+"_"+utility::to_string(atom_mapping.size())+".pdb");

		//Replace residues in the simplest way possible and dump (this will leave a bunch of nasty bond angles/lengths)
		core::conformation::Residue new_frame_res(hit_pose.residue(hit_resnum));
		new_frame_res.clear_residue_connections();
		query_pose.replace_residue(query_resnum, new_frame_res, false/*orient backbone*/);

		std::map<core::Size, core::Size>::const_iterator rsd_it = residue_map.begin();
		std::map<core::Size, core::Size>::const_iterator rsd_it_end = residue_map.end();
		for ( ; rsd_it != rsd_it_end; ++rsd_it ) {
			TR << "Replacing residue " << rsd_it->first << " in target with residue " << rsd_it->second << " from thermophilic motif" << std::endl;
			core::conformation::Residue new_res(hit_pose.residue(rsd_it->second));
			new_res.clear_residue_connections();
			query_pose.replace_residue(rsd_it->first, new_res, false/*orient backbone*/);
		}
		query_pose.dump_pdb("query_pose_"+utility::to_string(counter)+"_"+utility::to_string(hit_model_id)+"_"+utility::to_string(atom_mapping.size())+"_replaced.pdb");


		//Thread the sequence on using the packer and dump along with a coordinate constraints file
		core::pack::task::TaskFactoryOP task_factory = new core::pack::task::TaskFactory();
		core::pack::task::operation::RestrictResidueToRepackingOP res_repack = new core::pack::task::operation::RestrictResidueToRepacking();
		for ( core::Size i=1; i<=query_pose.size(); ++i ) {
			//If this isn't a mapped residue, don't mutate it, otherwise do
			if ( residue_map.find(i) == residue_map.end() ) {
				res_repack->include_residue(i);
			} else {
				core::pack::task::operation::RestrictAbsentCanonicalAASOP des = new core::pack::task::operation::RestrictAbsentCanonicalAAS();
				des->include_residue(i);
				des->keep_aas(std::string(1, hit_pose.residue(residue_map[i]).type().name1()));
				task_factory->push_back(des);
			}
		}
		task_factory->push_back(res_repack);

		protocols::simple_moves::PackRotamersMoverOP packer = new protocols::simple_moves::PackRotamersMover();
		packer->task_factory(task_factory);
		packer->apply(query_pose);
		core::pose::addVirtualResAsRoot(query_pose);
		query_pose.dump_pdb("query_pose_"+utility::to_string(counter)+"_"+utility::to_string(hit_model_id)+"_"+utility::to_string(atom_mapping.size())+"_repacked.pdb");

		//Setup coordinate constraints
		utility::io::ozstream constraint_file;
		constraint_file.open("constraint_"+utility::to_string(counter)+"_"+utility::to_string(hit_model_id)+"_"+utility::to_string(atom_mapping.size())+".cst");
		core::id::AtomID virt_root(1,query_pose.size());
		rsd_it = residue_map.begin();
		for ( ; rsd_it != rsd_it_end; ++rsd_it ) {
			for ( core::Size i = 1; i<=query_pose.residue(rsd_it->first).atoms().size(); ++i ) {

				constraint_file << "CoordinateConstraint ";

				constraint_file << query_pose.residue(rsd_it->first).type().atom_name(i) << " ";
				constraint_file << rsd_it->first << " ";
				constraint_file << query_pose.residue(query_pose.size()).type().atom_name(1) << " ";
				constraint_file << query_pose.size() << " ";

				numeric::xyzVector<core::Real> xyz = hit_pose.residue(rsd_it->second).atom(i).xyz();
				constraint_file << xyz.x() << " ";
				constraint_file << xyz.y() << " ";
				constraint_file << xyz.z() << " ";

				constraint_file << "HARMONIC 0.0 0.5" << std::endl;
			}
		}
		constraint_file.close();

		//  core::Real cst_std_dev = 0.5;
		//  query_pose.add_constraint(new core::scoring::constraints::CoordinateConstraint(it->first, virt_root, target_coords, func));
		//  core::scoring::func::HarmonicFuncOP func = new core::scoring::func::HarmonicFunc(0.0,cst_std_dev);
		//  rsd_it = residue_map.begin();

		//  core::scoring::ScoreFunctionOP scfxn = core::scoring::ScoreFunctionFactory::create_score_function("talaris2013_cart");
		//  protocols::relax::FastRelaxOP relax = new protocols::relax::FastRelax(1);
		//  relax->cartesian(true);
		//  relax->constrain_relax_to_start_coords(true);
		//  relax->min_type("lbfgs_armijo_nonmonotone");
		//  relax->set_scorefxn(scfxn);
		//  relax->apply(query_pose);
		//  query_pose.dump_pdb("query_pose_"+utility::to_string(counter)+"_"+utility::to_string(hit_model_id)+"_"+utility::to_string(atom_mapping.size())+"_relaxed.pdb");
		//  std::exit(1);

		//core::pack::task::TaskOperationOP op =
		//task_factory->push_back())


	}
}
