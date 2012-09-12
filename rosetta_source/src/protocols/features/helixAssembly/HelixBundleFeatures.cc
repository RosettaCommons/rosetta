// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :notabs=false:tabSize=4:indentsize=4:
//
// (c) copyright rosetta commons member institutions.
// (c) this file is part of the rosetta software suite and is made available under license.
// (c) the rosetta software is developed by the contributing members of the rosetta commons.
// (c) for more information, see http://www.rosettacommons.org. questions about this can be
// (c) addressed to university of washington uw techtransfer, email: license@u.washington.edu.

/// @file protocols/features/helixAssembly/HelixBundleFeatures.cc
/// @brief Search through a pose for sets of 3 helices and generate RMSD calculations between all pairs of them
/// @author Tim Jacobs

//Core
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>
#include <core/scoring/sasa.hh>
#include <core/pose/util.hh>

//External
#include <boost/uuid/uuid.hpp>

//Devel
#include <protocols/features/helixAssembly/HelixBundleFeatures.hh>
#include <protocols/features/helixAssembly/HelicalFragment.hh>

//Utility
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/tag/Tag.hh>

//C++
#include <string>
#include <math.h>

//External Headers
#include <cppdb/frontend.h>

//Basic
#include <basic/database/sql_utils.hh>
#include <basic/Tracer.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/helixAssembly.OptionKeys.gen.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/Constraint.hh>

static basic::Tracer TR("protocols.features.helixAssembly.HelixBundleFeatures");

namespace protocols {
namespace features {
namespace helixAssembly {

using namespace std;
using namespace core;
using core::pose::Pose;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using cppdb::result;

HelixBundleFeatures::HelixBundleFeatures() :
bundle_size_(3),
helix_size_(14),
helix_cap_dist_cutoff_(12.0)
{}

utility::vector1<std::string>
HelixBundleFeatures::features_reporter_dependencies() const
{
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	dependencies.push_back("ProteinResidueConformationFeatures");
	dependencies.push_back("ResidueSecondaryStructureFeatures");
	return dependencies;
}

void
HelixBundleFeatures::write_schema_to_db(utility::sql_database::sessionOP db_session) const
{
	using namespace basic::database::schema_generator;

	/******helix_bundles******/
	Column bundle_id("bundle_id", new DbInteger(), false /*not null*/, true /*autoincrement*/);
	Column struct_id("struct_id", new DbUUID(), false /*not null*/, false /*don't autoincrement*/);
	Column bundle_size("bundle_size", new DbInteger(), false, false);
	Column helix_size("helix_size", new DbInteger(), false, false);
	
	Schema helix_bundles("helix_bundles", PrimaryKey(bundle_id));
	helix_bundles.add_column(bundle_size);
	helix_bundles.add_column(helix_size);
	helix_bundles.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true /*defer*/));

	helix_bundles.write(db_session);

	/******bundle_helices******/
	Column helix_id("helix_id", new DbInteger(), false /*not null*/, true /*autoincrement*/);
	Column flipped("flipped", new DbInteger());
	Column sasa("sasa", new DbReal());
	Schema bundle_helices("bundle_helices", PrimaryKey(helix_id));

	bundle_helices.add_column(struct_id);
	bundle_helices.add_column(flipped);
	bundle_helices.add_column(sasa);

	Column residue_begin("residue_begin", new DbInteger());
	Column residue_end("residue_end", new DbInteger());


	utility::vector1<std::string> fkey_reference_cols;
	fkey_reference_cols.push_back("struct_id");
	fkey_reference_cols.push_back("resNum");

	utility::vector1<Column> fkey_cols_begin;
	fkey_cols_begin.push_back(struct_id);
	fkey_cols_begin.push_back(residue_begin);

	utility::vector1<Column> fkey_cols_end;
	fkey_cols_end.push_back(struct_id);
	fkey_cols_end.push_back(residue_end);

	ForeignKey bundle_id_fkey(Column("bundle_id", new DbInteger(),false), "helix_bundles", "bundle_id");

	bundle_helices.add_foreign_key(bundle_id_fkey);
	bundle_helices.add_foreign_key(ForeignKey(fkey_cols_begin, "residues", fkey_reference_cols, true /*defer*/));
	bundle_helices.add_foreign_key(ForeignKey(fkey_cols_end, "residues", fkey_reference_cols, true /*defer*/));

	bundle_helices.write(db_session);
}

//Select all helical segments reported by the ResidueSecondaryStructureFeatures and save them in a vector
utility::vector1<HelicalFragmentOP> HelixBundleFeatures::get_helix_fragments(boost::uuids::uuid struct_id, sessionOP db_session)
{
	std::string select_string =
	"SELECT\n"
	"	helices.helix_id,\n"
	"	helices.residue_begin,\n"
	"	helices.residue_end\n"
	"FROM\n"
	"	helix_segments as helices\n"
	"WHERE\n"
	"	helices.struct_id = ?\n"
	"ORDER BY residue_begin;";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	result res(basic::database::safely_read_from_database(select_statement));

	utility::vector1<HelicalFragmentOP> all_helices;
	
	core::Size prev_residue_begin;
	core::Size prev_residue_end;
	while(res.next())
	{
		Size helix_id, residue_begin, residue_end;
		res >> helix_id >> residue_begin >> residue_end;
		
		if(residue_begin - prev_residue_end == 2 && all_helices.size() > 0)
		{
			all_helices.pop_back();
			all_helices.push_back(new HelicalFragment(prev_residue_begin, residue_end));
			TR  << "combining helix segments: " << prev_residue_begin << "-" << prev_residue_end
				<< " and " << residue_begin << "-" << residue_end << std::endl;
		}
		else
		{
			all_helices.push_back(new HelicalFragment(residue_begin, residue_end));
		}
		prev_residue_begin=residue_begin;
		prev_residue_end=residue_end;
	}

	TR << "number of uninterrupted helix_segments: " << all_helices.size() << std::endl;
	utility::vector1<HelicalFragmentOP> all_helix_fragments;
	for(core::Size i=1; i<=all_helices.size(); ++i)
	{
		HelicalFragmentOP cur_helix=all_helices[i];
		
		if(cur_helix->end() >= helix_size_)
		{
			for(core::Size start_res=cur_helix->start();
				start_res<=cur_helix->end()-helix_size_+1;
				++start_res)
			{
				
				core::Size end_res = start_res+helix_size_-1;
				
				HelicalFragmentOP frag_1 = new HelicalFragment(start_res, end_res);//normal direction
				HelicalFragmentOP frag_2 = new HelicalFragment(end_res, start_res);//reversed
													
				all_helix_fragments.push_back(frag_1);
				all_helix_fragments.push_back(frag_2);
			}
		}
	}

	return all_helix_fragments;
}

void
HelixBundleFeatures::record_helix_sasas(
	core::pose::Pose const & pose,
	utility::vector1<HelicalFragmentOP> & bundle_fragments
){
	utility::vector1<core::Size> positions;
	utility::vector1< std::pair<core::Size, core::Size> > helix_ends;
	core::Size prev_end=0;
	for(core::Size i=1; i<=bundle_fragments.size(); ++i)
	{
		//Record residues from the starting pose for this bundle
		for(core::Size cur_res=bundle_fragments[i]->seq_start();
			cur_res<=bundle_fragments[i]->seq_end();
			++cur_res)
		{
			positions.push_back(cur_res);
		}
		
		//Get the ends of each helix in the new bundle
		core::Size new_start=prev_end+1;
		core::Size new_end=new_start+helix_size_-1;
		helix_ends.push_back(std::make_pair(new_start, new_end));
		prev_end=new_end;
	}
	
	assert(helix_ends.size()==bundle_fragments.size());
	
	core::pose::Pose bundle_pose;
	for(core::Size j=1; j<=helix_ends.size(); ++j)
	{
		kinematics::FoldTree ft;
		core::Size jump_counter=1;
		if(j!=1)
		{
			//peptide edge from 1->residue before helix to be separated
			ft.add_edge(1, helix_ends[j-1].second, kinematics::Edge::PEPTIDE);
		
			//jump from 1 to the edge to be separated
			ft.add_edge(1, helix_ends[j].first, jump_counter);
			jump_counter++;
		}
		
		//peptide edge for the helix to be separated
		ft.add_edge(helix_ends[j].first, helix_ends[j].second, kinematics::Edge::PEPTIDE);
		
		//add jump and peptide edge for the remainder of the bundle if the helix to be separated isn't the last one
		if(j!=helix_ends.size())
		{
			ft.add_edge(1, helix_ends[j+1].first, jump_counter);
			ft.add_edge(helix_ends[j+1].first, helix_ends[helix_ends.size()].second, kinematics::Edge::PEPTIDE);
		}
		
		ft.reorder(1);
		TR << "Fold tree " << j << ": " << ft << std::endl;
		
		//Only make the pose the first time, otherwise just change the fold tree
		if(j==1)
		{
			core::pose::create_subpose(pose, positions, ft, bundle_pose);
		}
		else
		{
			bundle_pose.fold_tree(ft);
		}
				
		sasa_filter_.jump(1);
		core::Real dSasa = sasa_filter_.compute(bundle_pose);
		bundle_fragments[j]->sasa(dSasa);
	}
}

void
HelixBundleFeatures::generate_comparison_list(
	utility::vector1<HelicalFragmentOP> all_helix_fragments,
	core::Size prev_index,
	utility::vector1<HelicalFragmentOP> fragment_list
){
	if(fragment_list.size()==bundle_size_)
	{
		comparison_list_.push_back(fragment_list);
	}
	else
	{
		for(core::Size i=prev_index+1; i<=all_helix_fragments.size(); ++i)
		{
			bool overlapping=false;
			for(core::Size j=1; j<=fragment_list.size(); ++j){
				if((all_helix_fragments[i]->seq_start() >= fragment_list[j]->seq_start() &&
					all_helix_fragments[i]->seq_start() <= fragment_list[j]->seq_end()) ||
				   (all_helix_fragments[i]->seq_end() >= fragment_list[j]->seq_start() &&
					all_helix_fragments[i]->seq_end() <= fragment_list[j]->seq_end()))
				{
					overlapping=true;
					break;
				}
					 
			}
			
			//don't allow fragment lists with two overlapping helices, or fragment
			//lists that start with a reversed helix (this prevents double reporting)
			if( !overlapping &&
				(fragment_list.size()>0 || !all_helix_fragments[i]->reversed()))
			{
				fragment_list.push_back(all_helix_fragments[i]);
				generate_comparison_list(all_helix_fragments, i, fragment_list);
				fragment_list.pop_back();
			}
		}
	}
}

bool
HelixBundleFeatures::check_cap_distances(
	core::pose::Pose const & pose,
	utility::vector1<HelicalFragmentOP> frag_set
){
	core::Real dist_sq_cutoff = pow(helix_cap_dist_cutoff_, 2);
	for(core::Size i=1; i<=frag_set.size(); ++i){
			
		for(core::Size j=i+1; j<=frag_set.size(); ++j){
			
			core::Real helix_start_dist_sq = pose.residue(frag_set[i]->start()).atom("CA").xyz().distance_squared(
				pose.residue(frag_set[j]->start()).atom("CA").xyz());
			
			core::Real helix_end_dist_sq = pose.residue(frag_set[i]->end()).atom("CA").xyz().distance_squared(
				pose.residue(frag_set[j]->end()).atom("CA").xyz());
				
			if(helix_start_dist_sq > dist_sq_cutoff || helix_end_dist_sq > dist_sq_cutoff){
				return false;
			}
		}
	}
	return true;
}


void
HelixBundleFeatures::parse_my_tag(
	utility::tag::TagPtr const tag,
	protocols::moves::DataMap & data,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & pose
){
	runtime_assert(tag->getOption<string>("name") == type_name());
	
	if(tag->hasOption("bundle_size")){
		bundle_size_ = tag->getOption<core::Size>("bundle_size");
	}
	if(tag->hasOption("helix_size")){
		helix_size_ = tag->getOption<core::Size>("helix_size");
	}
	if(tag->hasOption("helix_cap_dist_cutoff")){
		helix_cap_dist_cutoff_ = tag->getOption<core::Real>("helix_cap_dist_cutoff");
	}
}
	
core::Size
HelixBundleFeatures::report_features(
	core::pose::Pose const & pose,
	utility::vector1<bool> const &,
	boost::uuids::uuid struct_id,
	utility::sql_database::sessionOP db_session
){
	utility::vector1<HelicalFragmentOP> all_helix_fragments = get_helix_fragments(struct_id, db_session);
	TR << "Total helical fragments of size " << helix_size_ << ": " << all_helix_fragments.size() << std::endl;
	
	utility::vector1<HelicalFragmentOP> empty_frag_list;
	generate_comparison_list(all_helix_fragments, 0, empty_frag_list);
	
	TR << "Comparison list size: " << comparison_list_.size() << std::endl;
	
	for(core::Size cur_compare_index=1;
		cur_compare_index<=comparison_list_.size();
		++cur_compare_index)
	{
		
		utility::vector1<HelicalFragmentOP> cur_frag_set = comparison_list_[cur_compare_index];
		if(check_cap_distances(pose, cur_frag_set))
		{
			//record the dSasa for each helix
			record_helix_sasas(pose, cur_frag_set);
			
			string bundle_insert =  "INSERT INTO helix_bundles (struct_id, bundle_size, helix_size)  VALUES (?, ?, ?);";
			statement bundle_insert_stmt(basic::database::safely_prepare_statement(bundle_insert,db_session));
			bundle_insert_stmt.bind(1,struct_id);
			bundle_insert_stmt.bind(2,bundle_size_);
			bundle_insert_stmt.bind(3,helix_size_);
			basic::database::safely_write_to_database(bundle_insert_stmt);
			
			//Get bundle primary key
			core::Size bundle_id(bundle_insert_stmt.sequence_last("helix_bundles_bundle_id_seq"));
			
			string helix_insert =  "INSERT INTO bundle_helices (bundle_id, struct_id, residue_begin, residue_end, flipped, sasa) VALUES (?,?,?,?,?,?);";
			statement helix_insert_stmt(basic::database::safely_prepare_statement(helix_insert,db_session));
			TR << "Saving helices" << endl;
			
			for(core::Size frag_index=1;
				frag_index <= cur_frag_set.size();
				++frag_index)
			{
				helix_insert_stmt.bind(1,bundle_id);
				helix_insert_stmt.bind(2,struct_id);
				helix_insert_stmt.bind(3,cur_frag_set[frag_index]->seq_start());
				helix_insert_stmt.bind(4,cur_frag_set[frag_index]->seq_end());
				helix_insert_stmt.bind(5,cur_frag_set[frag_index]->reversed());
				helix_insert_stmt.bind(6,cur_frag_set[frag_index]->sasa());
				basic::database::safely_write_to_database(helix_insert_stmt);
			}
			TR << "Done saving helices" << endl;			
		}
	}
	return 0;
}

} //namespace helixAssembly
} //namespace features
} //namespace protocols
