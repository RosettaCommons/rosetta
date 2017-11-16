// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DatabaseEntryWorkUnit.cc
///
/// @brief A work unit that runs a database query, processes the results, and returns a string (presumably a database insert statement)

/// @author Tim Jacobs

#include <core/pose/util.hh>
#include <core/types.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <basic/database/sql_utils.hh>
#include <basic/Tracer.hh>

#include <devel/init.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>

#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>

//Utility headers
#include <utility/vector1.hh>

//Protocol headers
#include <protocols/features/ProteinSilentReport.hh>
#include <protocols/sewing/hashing/Hasher.hh>
#include <core/scoring/rms_util.hh>

//C++ headers
#include <map>

static basic::Tracer TR("SmotifHasher");

namespace SmotifHasher {
basic::options::IntegerOptionKey const num_bundles( "num_bundles" );
basic::options::IntegerOptionKey const mode( "mode" ); //1 is hash all smotifs, 2 is score smotifs
basic::options::IntegerOptionKey const query_id( "query_id" );
basic::options::IntegerOptionKey const min_score( "min_score" );
}

std::map< core::Size, utility::vector1<protocols::sewing::ResidueHash> >
get_models(
	utility::sql_database::sessionOP db_session
){
	std::string coords_select_string =
		"SELECT s.smotif_id, coords.seqpos, coords.atomno, coords.x, coords.y, coords.z\n"
		"FROM smotifs s\n"
		"JOIN secondary_structure_segments ss ON\n"
		"       s.struct_id = ss.struct_id AND\n"
		"       (s.secondary_struct_segment_id_1 = ss.segment_id OR\n"
		"       s.secondary_struct_segment_id_2 = ss.segment_id)\n"
		"JOIN residue_atom_coords coords ON\n"
		"       s.struct_id = coords.struct_id AND\n"
		"       coords.seqpos BETWEEN ss.residue_begin AND ss.residue_end\n"
		"WHERE coords.atomno IN (1,2,3)\n"
		"ORDER BY s.smotif_id, coords.seqpos, coords.atomno;";

	cppdb::statement coords_select_stmt=basic::database::safely_prepare_statement(coords_select_string, db_session);
	cppdb::result coords_res=basic::database::safely_read_from_database(coords_select_stmt);


	core::Size smotif_id, seqpos, atomno;
	core::Real x, y, z;

	std::map< core::Size, utility::vector1<protocols::sewing::ResidueHash> > models_to_hash;
	core::Size last_smotif, last_seqpos=0;

	protocols::sewing::ResidueHash cur_residue;
	while ( coords_res.next() ) {
		coords_res >> smotif_id >> seqpos >> atomno >> x >> y >> z;

		if ( seqpos != last_seqpos && last_seqpos != 0 ) {
			runtime_assert(cur_residue.atoms.size() == 3);
			models_to_hash[last_smotif].push_back(cur_residue);
			cur_residue = protocols::sewing::ResidueHash();
		}

		cur_residue.resnum = seqpos;

		protocols::sewing::AtomHash atom;
		atom.atomno = atomno;
		atom.coords = numeric::xyzVector<core::Real>(x,y,z);
		cur_residue.atoms.push_back(atom);

		last_seqpos = seqpos;
		last_smotif = smotif_id;
	}
	models_to_hash[last_smotif].push_back(cur_residue);

	TR << "Generated " << models_to_hash.size() << " models from smotifs table" << std::endl;
	return models_to_hash;
}

void
trim_pdb(core::pose::Pose & pose, std::set<core::Size> const & residue_numbers) {
	core::Size total_res=pose.size();
	int num_removed_residues=0;
	for ( core::Size i=1; i<=total_res; ++i ) {
		//if the residue wasn't specified, delete it
		if ( residue_numbers.find(i)==residue_numbers.end() ) {
			pose.conformation().delete_residue_slow(i-num_removed_residues);
			++num_removed_residues;
		}
	}
}

void
superimpose_smotifs(
	utility::sql_database::sessionOP db_session,
	std::map< std::pair<protocols::sewing::HashValue, protocols::sewing::HashValue>, core::Size > const & scores,
	core::Size score_cutoff
) {
	using namespace protocols::sewing;

	std::string coords_select_string =
		"SELECT s.struct_id, s.smotif_id, coords.seqpos\n"
		"FROM smotifs s\n"
		"JOIN secondary_structure_segments ss ON\n"
		"       s.struct_id = ss.struct_id AND\n"
		"       (s.secondary_struct_segment_id_1 = ss.segment_id OR\n"
		"       s.secondary_struct_segment_id_2 = ss.segment_id)\n"
		"JOIN residue_atom_coords coords ON\n"
		"       s.struct_id = coords.struct_id AND\n"
		"       coords.seqpos BETWEEN ss.residue_begin AND ss.residue_end\n"
		"WHERE s.smotif_id=?;";
	cppdb::statement coords_select_stmt=basic::database::safely_prepare_statement(coords_select_string, db_session);

	protocols::features::ProteinSilentReportOP protein_silent_report = new protocols::features::ProteinSilentReport();

	//First dump the query smotifs
	coords_select_stmt.bind(1, scores.begin()->first.first.model_id);
	cppdb::result query_coords_res=basic::database::safely_read_from_database(coords_select_stmt);

	std::set<core::Size> query_resnums;
	protocols::features::StructureID struct_id;
	core::Size smotif_id, seqpos;
	while ( query_coords_res.next() ) {
		query_coords_res >> struct_id >> smotif_id >> seqpos;
		query_resnums.insert(seqpos);
	}
	core::pose::Pose query_pose;
	protein_silent_report->load_pose(db_session, struct_id, query_pose);

	core::pose::Pose trimmed_query(query_pose);
	trim_pdb(trimmed_query, query_resnums);
	trimmed_query.dump_pdb("query_" + utility::to_string(smotif_id) + ".pdb");

	for ( std::map< std::pair< HashValue, HashValue >, core::Size>::const_iterator it=scores.begin(); it!=scores.end(); ++it ) {
		if ( it->second >= score_cutoff ) {
			HashValue const & query = it->first.first;
			HashValue const & hit = it->first.second;
			TR << "Superimposing query " << query.model_id << "(resnum " << query.resnum << ") with hit " << hit.model_id << " (resnum " << hit.resnum << "). Score is " << it->second << std::endl;

			coords_select_stmt.bind(1, hit.model_id);
			cppdb::result coords_res=basic::database::safely_read_from_database(coords_select_stmt);

			std::set<core::Size> resnums;
			while ( coords_res.next() ) {
				coords_res >> struct_id >> smotif_id >> seqpos;
				resnums.insert(seqpos);
			}
			core::pose::Pose hit_pose;
			protein_silent_report->load_pose(db_session, struct_id, hit_pose);

			core::id::AtomID_Map< core::id::AtomID > atom_map( core::id::AtomID::BOGUS_ATOM_ID() );
			core::pose::initialize_atomid_map( atom_map, hit_pose, core::id::AtomID::BOGUS_ATOM_ID() );
			atom_map.set( core::id::AtomID(hit_pose.residue( hit.resnum ).atom_index("N"), hit.resnum ), core::id::AtomID( query_pose.residue(query.resnum).atom_index("N"), query.resnum ) );
			atom_map.set( core::id::AtomID(hit_pose.residue( hit.resnum ).atom_index("CA"), hit.resnum ), core::id::AtomID( query_pose.residue(query.resnum).atom_index("CA"), query.resnum ) );
			atom_map.set( core::id::AtomID(hit_pose.residue( hit.resnum ).atom_index("C"), hit.resnum ), core::id::AtomID( query_pose.residue(query.resnum).atom_index("C"), query.resnum ) );
			core::scoring::superimpose_pose(hit_pose, query_pose, atom_map);

			trim_pdb(hit_pose, resnums);
			hit_pose.dump_pdb("hit_" + utility::to_string(smotif_id) + ".pdb");

		}
	}

}

int
main( int argc, char * argv [] ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	option.add( SmotifHasher::mode, "The mode");
	option.add( SmotifHasher::query_id, "The smotif to query");
	option.add( SmotifHasher::min_score, "minimum score to print");

	// initialize core and read options
	devel::init(argc, argv);
	core::Size mode = option[SmotifHasher::mode].def(1);
	core::Size query_id = option[SmotifHasher::query_id].def(1);
	core::Size min_score = option[SmotifHasher::min_score].def(10);

	// Initialize DB
	utility::sql_database::sessionOP db_session( basic::database::get_db_session() );

	std::map<core::Size, utility::vector1< protocols::sewing::ResidueHash > > models = get_models(db_session);
	protocols::sewing::Hasher hasher;

	TR << "Smotif Hasher options:" << std::endl;
	TR << "\tMode: " << mode << std::endl;
	TR << "\tQuery ID: " << query_id << std::endl;
	TR << "\tMinimum Score: " << min_score << std::endl;

	if ( mode == 1 ) {
		for ( std::map<core::Size, utility::vector1< protocols::sewing::ResidueHash > >::const_iterator it = models.begin(); it!=models.end(); ++it ) {
			//   TR << "About to hash model: " << it->first << std::endl;
			//   for(core::Size i=1; i<=it->second.size(); ++i) {
			//    TR << "Residue " << it->second[i].resnum << std::endl;
			//   }
			hasher.insert(it->first, it->second);
		}
		hasher.write_to_disk("smotifs.hashtable");
	} else if ( mode ==2 ) {
		hasher.read_from_disk("smotifs.hashtable");
		TR << "Finished reading hash table from disk" << std::endl;

		std::map< std::pair<protocols::sewing::HashValue, protocols::sewing::HashValue>, core::Size > scores;

		core::Size counter=0;
		core::Size ten_percent=models.size()/10;
		for ( std::map<core::Size, utility::vector1< protocols::sewing::ResidueHash > >::const_iterator it = models.begin(); it!=models.end(); ++it ) {
			if ( query_id==0 ) {
				++counter;
				if ( counter % ten_percent == 0 ) {
					TR << "Completed " << counter << " of " << models.size() << std::endl;;
				}
				std::map< std::pair<protocols::sewing::HashValue, protocols::sewing::HashValue>, core::Size > temp_scores = hasher.score(it->first, it->second);
				scores.insert(temp_scores.begin(), temp_scores.end());
			} else if ( it->first == query_id ) {
				scores = hasher.score(it->first, it->second);
			}
		}

		superimpose_smotifs(db_session, scores, min_score);
	} else {
		utility_exit_with_message("Unrecognized mode");
	}

	//superimpose_smotifs(coordinates.begin()->first, best_hit.first.id, best_hit.first.bs);

}
