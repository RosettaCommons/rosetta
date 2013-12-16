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
#include <core/pose/symmetry/util.hh>

//Devel
#include <protocols/features/helixAssembly/HelixBundleFeatures.hh>
#include <protocols/features/helixAssembly/HelicalFragment.hh>

//Utility
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/tag/Tag.hh>
#include <utility/LexicographicalIterator.hh>

//Core
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreType.hh>
#include <core/graph/Graph.hh>

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

//Numeric
#include <numeric/PCA.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.hh>

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
		helix_cap_dist_cutoff_(12.0),
		bundle_size_(3),
		helix_size_(14)
{
	scorefxn_ = core::scoring::getScoreFunction();
}

utility::vector1<std::string>
HelixBundleFeatures::features_reporter_dependencies() const
{
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	dependencies.push_back("ProteinResidueConformationFeatures");
	dependencies.push_back("ResidueSecondaryStructureFeatures");
	dependencies.push_back("SecondaryStructureSegmentFeatures");
	return dependencies;
}

void
HelixBundleFeatures::write_schema_to_db(utility::sql_database::sessionOP db_session) const
{
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", new DbBigInt(), false /*not null*/, false /*don't autoincrement*/);

	/******helix_bundles******/
	Column bundle_id("bundle_id", new DbInteger(), false /*not null*/, false /*autoincrement*/);
//	Column bundle_id("bundle_id", new DbInteger(), false /*not null*/, true /*autoincrement*/);
	Column bundle_size("bundle_size", new DbInteger(), false, false);
	Column helix_size("helix_size", new DbInteger(), false, false);
	
	utility::vector1<Column> helix_bundle_pkey_cols;
	helix_bundle_pkey_cols.push_back(struct_id);
	helix_bundle_pkey_cols.push_back(bundle_id);
	
	Schema helix_bundles("helix_bundles", helix_bundle_pkey_cols);
	helix_bundles.add_column(bundle_size);
	helix_bundles.add_column(helix_size);
	helix_bundles.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true /*defer*/));

	helix_bundles.write(db_session);
	
	/******bundle_helices******/
	Column helix_id("helix_id", new DbInteger(), false /*not null*/, false /*autoincrement*/);
//	Column helix_id("helix_id", new DbInteger(), false /*not null*/, true /*autoincrement*/);
	Column flipped("flipped", new DbInteger());
	Column sasa("sasa", new DbReal());

	Column residue_begin("residue_begin", new DbInteger());
	Column residue_end("residue_end", new DbInteger());
	
	utility::vector1<Column> bundle_helix_pkey_cols;
	bundle_helix_pkey_cols.push_back(struct_id);
	bundle_helix_pkey_cols.push_back(helix_id);

	utility::vector1<std::string> fkey_res_reference_cols;
	fkey_res_reference_cols.push_back("struct_id");
	fkey_res_reference_cols.push_back("resNum");

	utility::vector1<Column> fkey_cols_begin;
	fkey_cols_begin.push_back(struct_id);
	fkey_cols_begin.push_back(residue_begin);

	utility::vector1<Column> fkey_cols_end;
	fkey_cols_end.push_back(struct_id);
	fkey_cols_end.push_back(residue_end);
	
	utility::vector1<std::string> fkey_bundle_reference_cols;
	fkey_bundle_reference_cols.push_back("struct_id");
	fkey_bundle_reference_cols.push_back("bundle_id");
	
	utility::vector1<Column> fkey_bundle_cols;
	fkey_bundle_cols.push_back(struct_id);
	fkey_bundle_cols.push_back(bundle_id);


	Schema bundle_helices("bundle_helices", bundle_helix_pkey_cols);
	bundle_helices.add_foreign_key(ForeignKey(fkey_bundle_cols, "helix_bundles", fkey_bundle_reference_cols, true /*defer*/));
	bundle_helices.add_foreign_key(ForeignKey(fkey_cols_begin, "residues", fkey_res_reference_cols, true /*defer*/));
	bundle_helices.add_foreign_key(ForeignKey(fkey_cols_end, "residues", fkey_res_reference_cols, true /*defer*/));
	bundle_helices.add_column(struct_id);
	bundle_helices.add_column(flipped);
	bundle_helices.add_column(sasa);
	
	bundle_helices.write(db_session);
	
	/******helix_pairs******/
	Column pair_id("pair_id", new DbInteger(), false /*not null*/, false /*autoincrement*/);
//	Column pair_id("pair_id", new DbInteger(), false /*not null*/, true /*autoincrement*/);
	Column helix_id_1("helix_id_1", new DbInteger(), false, false);
	Column helix_id_2("helix_id_2", new DbInteger(), false, false);
	Column shared_fa_atr("shared_fa_atr", new DbReal(), false, false);
	Column interacting_fraction("interacting_fraction", new DbReal(), false, false);
	Column crossing_angle("crossing_angle", new DbReal(), false, false);
	Column end_1_distance("end_1_distance", new DbReal(), false, false);
	Column end_2_distance("end_2_distance", new DbReal(), false, false);
	
	utility::vector1<Column> helix_pair_pkey_cols;
	helix_pair_pkey_cols.push_back(struct_id);
	helix_pair_pkey_cols.push_back(pair_id);
	
	utility::vector1<std::string> fkey_ref_helix;
	fkey_ref_helix.push_back("struct_id");
	fkey_ref_helix.push_back("helix_id");
	
	utility::vector1<Column> fkey_cols_helix_1;
	fkey_cols_helix_1.push_back(struct_id);
	fkey_cols_helix_1.push_back(helix_id_1);
	
	utility::vector1<Column> fkey_cols_helix_2;
	fkey_cols_helix_2.push_back(struct_id);
	fkey_cols_helix_2.push_back(helix_id_2);
	
	Schema helix_pairs("helix_pairs", helix_pair_pkey_cols);
	helix_pairs.add_column(helix_id_1);
	helix_pairs.add_column(helix_id_2);
	helix_pairs.add_column(shared_fa_atr);
	helix_pairs.add_column(interacting_fraction);
	helix_pairs.add_column(crossing_angle);
	helix_pairs.add_column(end_1_distance);
	helix_pairs.add_column(end_2_distance);
	helix_pairs.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true /*defer*/));
	helix_pairs.add_foreign_key(ForeignKey(fkey_bundle_cols, "helix_bundles", fkey_bundle_reference_cols, true /*defer*/));
	helix_pairs.add_foreign_key(ForeignKey(fkey_cols_helix_1, "bundle_helices", fkey_ref_helix, true /*defer*/));
	helix_pairs.add_foreign_key(ForeignKey(fkey_cols_helix_2, "bundle_helices", fkey_ref_helix, true /*defer*/));
	
	helix_pairs.write(db_session);
}

//Select all helical segments reported by the ResidueSecondaryStructureFeatures and save them in a vector
utility::vector1<HelicalFragmentOP> HelixBundleFeatures::get_helix_fragments(StructureID struct_id, sessionOP db_session)
{
	std::string select_string =
		"SELECT\n"
		"	residue_begin,\n"
		"	residue_end\n"
		"FROM\n"
		"	secondary_structure_segments\n"
		"WHERE\n"
		"	struct_id = ?\n AND "
		"	dssp = 'H'\n"
		"ORDER BY residue_begin;";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	result res(basic::database::safely_read_from_database(select_statement));

	utility::vector1<HelicalFragmentOP> all_helices;
	
	core::Size prev_residue_begin = 0;
	core::Size prev_residue_end = 0;
	while(res.next())
	{
		Size residue_begin, residue_end;
		res >> residue_begin >> residue_end;
		
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

PairMap
HelixBundleFeatures::get_helix_pairs(
	core::pose::Pose const & pose,
	utility::vector1<HelicalFragmentOP> helix_fragments
){
	core::Real dist_sq_cutoff = pow(helix_cap_dist_cutoff_, 2);
	
	PairMap pair_map;
	for(core::Size i=1; i<=helix_fragments.size(); ++i)
	{
		for(core::Size j=i+1; j<=helix_fragments.size(); ++j)
		{
			//ensure pairs don't contain overlapping fragments
			if((helix_fragments[i]->seq_start() >= helix_fragments[j]->seq_start() &&
				helix_fragments[i]->seq_start() <= helix_fragments[j]->seq_end()) ||
				
			   (helix_fragments[i]->seq_end() >= helix_fragments[j]->seq_start() &&
				helix_fragments[i]->seq_end() <= helix_fragments[j]->seq_end()))
			{
				continue;
			}
			
			FragmentPair helix_pair(helix_fragments[i], helix_fragments[j]);
			
			//check to see if com's of the two fragments are within cutoff
//			core::Real com_dist_sq = helix_pair.fragment_1->com().distance_squared(
//				helix_pair.fragment_2->com());
//			if(com_dist_sq > dist_sq_cutoff)
//			{
//				continue;
//			}
//			helix_pair.com_distance = std::sqrt(com_dist_sq);

			core::Real end_1_dist_sq = pose.residue(helix_pair.fragment_1->start()).atom("CA").xyz().distance_squared(
				pose.residue(helix_pair.fragment_2->start()).atom("CA").xyz());
			
			core::Real end_2_dist_sq = pose.residue(helix_pair.fragment_1->end()).atom("CA").xyz().distance_squared(
				pose.residue(helix_pair.fragment_2->end()).atom("CA").xyz());
				
			if(end_1_dist_sq > dist_sq_cutoff || end_2_dist_sq > dist_sq_cutoff)
			{
				continue;
			}
			helix_pair.end_1_distance = std::sqrt(end_1_dist_sq);
			helix_pair.end_2_distance = std::sqrt(end_2_dist_sq);
			
//			TR << "end 1 distance (" << helix_pair.fragment_1->start() << "," << helix_pair.fragment_2->start()
//				<< "): " << helix_pair.end_1_distance << std::endl;
//				
//			TR << "end 2 distance (" << helix_pair.fragment_1->end() << "," << helix_pair.fragment_2->end()
//				<< "): " << helix_pair.end_2_distance << std::endl;

			//check to see per-residue fa_attr and percentage of interacting residues are within cutoff
			calc_fa_energy(pose, helix_pair);
			if( helix_pair.fa_attr > -(min_per_residue_fa_attr_) ||
			    helix_pair.fa_fraction < min_interacting_set_fraction_
			){
				continue;
			}
			
			//ensure crossing angle is within cutoff
			calc_crossing_angles(pose, helix_pair);
//			if((helix_pair.crossing_angle < 0 && helix_pair.crossing_angle > (-180+max_degrees_off_parallel_) ) ||
//				helix_pair.crossing_angle > max_degrees_off_parallel_)
//			{
//				continue;
//			}
			
			std::pair<core::Size, core::Size> pair_index = std::make_pair(i, j);
			pair_map.insert(std::make_pair(pair_index, helix_pair));
		}
	}
	return pair_map;
}

///@brief calculate the shared fa_attr for each pair of helices
/// in the bundle
void
HelixBundleFeatures::calc_fa_energy(
	core::pose::Pose const & pose,
	FragmentPair & helix_pair
){
	HelicalFragmentOP frag_1=helix_pair.fragment_1;
	HelicalFragmentOP frag_2=helix_pair.fragment_2;
	
	core::Real helix_pair_fa_attr=0;
	std::set<core::Size> interacting_residues; //residues that make favorable inter-fragment contacts (negative fa_attr)
	for(core::Size ii=frag_1->seq_start(); ii<=frag_1->seq_end(); ++ii)
	{
		//iterate through interaction graph of current residue and add any fa_attr
		//between this residue and any of the residues from fragment j
		core::scoring::EnergyGraph const & eg =
		pose.energies().energy_graph();
		for (core::graph::Node::EdgeListConstIter
			 iter = eg.get_node(ii)->const_upper_edge_list_begin(),
			 iter_end = eg.get_node(ii)->const_upper_edge_list_end();
			 iter != iter_end; ++iter )
		{
			core::Size jj = (*iter)->get_other_ind(ii);
			if(jj >= frag_2->seq_start() && jj <= frag_2->seq_end())
			{
				const core::scoring::EnergyEdge* ee =
				dynamic_cast< const core::scoring::EnergyEdge* > (*iter);
				helix_pair_fa_attr += (*ee)[core::scoring::fa_atr];
				if((*ee)[core::scoring::fa_atr] < 0)
				{
					interacting_residues.insert(ii);
					interacting_residues.insert(jj);
				}
			}
		}
	}
	
	//normalize fragment pair fa attr by number of residues and check for per-residue threshold
	helix_pair.fa_fraction = interacting_residues.size()/(helix_size_*2.0);
	helix_pair.fa_attr = helix_pair_fa_attr/(helix_size_*2.0);
	
//	TR << "residue-normalized FA attractive between: " << (*frag_1)
//		<< " " << (*frag_2) << ": " << per_residue_fa_attr << std::endl;
//	
//	TR << "interacting set fraction: " << (*frag_1)
//		<< " " << (*frag_2) << ": " << interacting_set_fraction << std::endl;
}


///@brief calculate the crossing angles of the helix fragment in the bundle set
void
HelixBundleFeatures::calc_crossing_angles(
	core::pose::Pose const & pose,
	FragmentPair & helix_pair
){
	HelicalFragmentOP frag_1=helix_pair.fragment_1;
	HelicalFragmentOP frag_2=helix_pair.fragment_2;
	
	//populate coms and principal components
	calc_pc_and_com(pose, frag_1);
	calc_pc_and_com(pose, frag_2);
	
	//Find the closest residues between fragments i and j
//	std::pair<core::Size,core::Size> closest_residues;
//	core::Real min_dist_sq=1000000;
//	for(core::Size ii=frag_1->seq_start(); ii<=frag_1->seq_end(); ++ii)
//	{
//		for(core::Size jj=frag_2->seq_start(); jj<=frag_2->seq_end(); ++jj)
//		{
//			core::Real dist_sq = pose.residue(ii).atom("CA").xyz().distance_squared(
//				pose.residue(jj).atom("CA").xyz());
//			if(dist_sq < min_dist_sq)
//			{
//				min_dist_sq=dist_sq;
//				closest_residues=std::make_pair(ii,jj);
//			}
//		}
//	}

//	numeric::xyzVector<core::Real> p0 =
//		closest_point_on_line(frag_1->com(), frag_1->principal_component(),
//		pose.residue(frag_1->seq_start()).atom("N").xyz());
//
//	numeric::xyzVector<core::Real> p1 =
//		closest_point_on_line(frag_1->com(), frag_1->principal_component(),
//		pose.residue(closest_residues.first).atom("CA").xyz());
//		
//	numeric::xyzVector<core::Real> p2 =
//		closest_point_on_line(frag_2->com(), frag_2->principal_component(),
//		pose.residue(closest_residues.second).atom("CA").xyz());
//		
//	numeric::xyzVector<core::Real> p3 =
//		closest_point_on_line(frag_2->com(), frag_2->principal_component(),
//		pose.residue(frag_2->seq_start()).atom("N").xyz());
		
	numeric::xyzVector<core::Real> p0 =
		closest_point_on_line(frag_1->com(), frag_1->principal_component(),
		pose.residue(frag_1->start()).atom("N").xyz());
		
//	TR << "P0: " << p0.x() << "," << p0.y() << "," << p0.z() << endl;

	numeric::xyzVector<core::Real> p1 =
		closest_point_on_line(frag_1->com(), frag_1->principal_component(),
		pose.residue(frag_1->end()).atom("CA").xyz());
		
//	TR << "P1: " << p1.x() << "," << p1.y() << "," << p1.z() << endl;

	numeric::xyzVector<core::Real> p2 =
		closest_point_on_line(frag_2->com(), frag_2->principal_component(),
		pose.residue(frag_2->end()).atom("CA").xyz());
		
//	TR << "P2: " << p2.x() << "," << p2.y() << "," << p2.z() << endl;

	numeric::xyzVector<core::Real> p3 =
		closest_point_on_line(frag_2->com(), frag_2->principal_component(),
		pose.residue(frag_2->start()).atom("N").xyz());
		
//	TR << "P3: " << p3.x() << "," << p3.y() << "," << p3.z() << endl;

	helix_pair.crossing_angle =
		numeric::dihedral_degrees(p0,p1,p2,p3);

//	TR << "crossing angle between: " << (*frag_1)
//		<< " " << (*frag_2) << ": " << crossing_angle << std::endl;
}
	

///@brief create a bundle-pose from the combination of fragments
/// and record the "interface" SASA for each helix against the
/// rest of the bundle
void
HelixBundleFeatures::record_helix_sasas(
	core::pose::Pose const & /*pose*/,
	std::set<HelicalFragmentOP> const & /*frag_set*/
){
//	utility::vector1<core::Size> positions;
//	utility::vector1< std::pair<core::Size, core::Size> > helix_ends;
//	core::Size prev_end=0;
//	for(core::Size i=1; i<=frag_set.size(); ++i)
//	{
//		//Record residues from the starting pose for this bundle
//		for(core::Size cur_res=frag_set[i]->seq_start();
//			cur_res<=frag_set[i]->seq_end();
//			++cur_res)
//		{
//			positions.push_back(cur_res);
//		}
//		
//		//Get the ends of each helix in the new bundle
//		core::Size new_start=prev_end+1;
//		core::Size new_end=new_start+helix_size_-1;
//		helix_ends.push_back(std::make_pair(new_start, new_end));
//		prev_end=new_end;
//	}
//	
//	assert(helix_ends.size()==frag_set.size());
//	
//	core::pose::Pose bundle_pose;
//	for(core::Size j=1; j<=helix_ends.size(); ++j)
//	{
//		kinematics::FoldTree ft;
//		core::Size jump_counter=1;
//		if(j!=1)
//		{
//			//peptide edge from 1->residue before helix to be separated
//			ft.add_edge(1, helix_ends[j-1].second, kinematics::Edge::PEPTIDE);
//		
//			//jump from 1 to the edge to be separated
//			ft.add_edge(1, helix_ends[j].first, jump_counter);
//			jump_counter++;
//		}
//		
//		//peptide edge for the helix to be separated
//		ft.add_edge(helix_ends[j].first, helix_ends[j].second, kinematics::Edge::PEPTIDE);
//		
//		//add jump and peptide edge for the remainder of the bundle if the helix to be separated isn't the last one
//		if(j!=helix_ends.size())
//		{
//			ft.add_edge(1, helix_ends[j+1].first, jump_counter);
//			ft.add_edge(helix_ends[j+1].first, helix_ends[helix_ends.size()].second, kinematics::Edge::PEPTIDE);
//		}
//		
//		ft.reorder(1);
//		TR << "Fold tree " << j << ": " << ft << std::endl;
//		
//		//Only make the pose the first time, otherwise just change the fold tree
//		if(j==1)
//		{
//			core::pose::create_subpose(pose, positions, ft, bundle_pose);
//		}
//		else
//		{
//			bundle_pose.fold_tree(ft);
//		}
//				
//		sasa_filter_.jump(1);
//		core::Real dSasa = sasa_filter_.compute(bundle_pose);
//		frag_set[j]->sasa(dSasa);
//	}
}

void
HelixBundleFeatures::calc_pc_and_com(
	core::pose::Pose const & pose,
	HelicalFragmentOP fragment
){
	utility::vector1< numeric::xyzVector< core::Real> > fragment_coords;
	for(core::Size cur_res=fragment->seq_start(); cur_res<=fragment->seq_end(); ++cur_res)
	{
		fragment_coords.push_back(pose.residue(cur_res).atom("N").xyz());
		fragment_coords.push_back(pose.residue(cur_res).atom("CA").xyz());
		fragment_coords.push_back(pose.residue(cur_res).atom("C").xyz());
		fragment_coords.push_back(pose.residue(cur_res).atom("O").xyz());
	}
	
	numeric::xyzVector<core::Real> com = numeric::center_of_mass(fragment_coords);
	fragment->com(com);
	fragment->principal_component(numeric::first_principal_component(fragment_coords) + com);
}

void
HelixBundleFeatures::parse_my_tag(
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/
){
	runtime_assert(tag->getOption<string>("name") == type_name());
	
	bundle_size_ = tag->getOption<core::Size>("bundle_size", 3);
	helix_size_ = tag->getOption<core::Size>("helix_size", 12);
	helix_cap_dist_cutoff_ = tag->getOption<core::Real>("helix_cap_dist_cutoff", 18);
	
	min_per_residue_fa_attr_ = tag->getOption<core::Real>("min_per_residue_fa_attr", 0.2);
	min_interacting_set_fraction_ = tag->getOption<core::Real>("min_interacting_set_fraction", 0.4);
	max_degrees_off_parallel_ = tag->getOption<core::Real>("max_degrees_off_parallel", 40);

}
	
core::Size
HelixBundleFeatures::report_features(
	core::pose::Pose const & pose,
	utility::vector1<bool> const &,
	StructureID struct_id,
	utility::sql_database::sessionOP db_session
){
	utility::vector1<HelicalFragmentOP> all_helix_fragments = get_helix_fragments(struct_id, db_session);
	if(all_helix_fragments.size()<bundle_size_){ return 0; }
	TR << "Total helical fragments of size " << helix_size_ << ": " << all_helix_fragments.size() << std::endl;
	
	//Non-const pose to score
	core::pose::Pose pose_copy(pose);
	core::pose::symmetry::make_score_function_consistent_with_symmetric_state_of_pose(pose_copy, scorefxn_);
	scorefxn_->score(pose_copy);
	PairMap pair_map = get_helix_pairs(pose_copy, all_helix_fragments);
	if(pair_map.size()==0){ return 0; };
	TR << "Total valid fragment pairs: " << pair_map.size() << std::endl;
	
	//Insert statements
	string bundle_insert =
//		"INSERT INTO helix_bundles (struct_id, bundle_size, helix_size)  VALUES (?, ?, ?);";
		"INSERT INTO helix_bundles (bundle_id, struct_id, bundle_size, helix_size)  VALUES (?, ?, ?, ?);";
	statement bundle_insert_stmt(basic::database::safely_prepare_statement(bundle_insert,db_session));
	
	string helix_insert =
//		"INSERT INTO bundle_helices (bundle_id, struct_id, residue_begin, residue_end, flipped, sasa) VALUES (?,?,?,?,?,?);";
		"INSERT INTO bundle_helices (helix_id, bundle_id, struct_id, residue_begin, residue_end, flipped, sasa) VALUES (?,?,?,?,?,?,?);";
	statement helix_insert_stmt(basic::database::safely_prepare_statement(helix_insert,db_session));
	
	string pair_insert =
//		"INSERT INTO helix_pairs (struct_id, bundle_id, helix_id_1, helix_id_2, shared_fa_atr, interacting_fraction, crossing_angle, end_1_distance, end_2_distance)\n"
//		"VALUES (?,?,?,?,?,?,?,?,?)";
		"INSERT INTO helix_pairs (pair_id, struct_id, bundle_id, helix_id_1, helix_id_2, shared_fa_atr, interacting_fraction, crossing_angle, end_1_distance, end_2_distance)\n"
		"VALUES (?,?,?,?,?,?,?,?,?,?)";
	statement pair_insert_stmt(basic::database::safely_prepare_statement(pair_insert, db_session));
	
	utility::vector1<core::Size> dim_sizes(bundle_size_, all_helix_fragments.size());
	bool done=false;
	core::Size bundle_counter=0;
	core::Size helix_counter=0;
	core::Size pair_counter=0;
	core::Size tracker=1;
	for ( utility::LexicographicalIterator lex( dim_sizes ); ! lex.at_end(); ++lex )
	{
		if(lex[1] != tracker){
			TR << "first lexico iterator is at: " << lex[1] << " of " << all_helix_fragments.size() << std::endl;
			tracker=lex[1];
		}
		//Only evaluate sets of fragments that are non-overlapping. and non duplicates (anytime lex[i] > lex[j] it's a duplicate)
		for(core::Size i=1; i<=bundle_size_; ++i){
			for(core::Size j=i+1; j<=bundle_size_; ++j){
				while( lex[i] > lex[j] || overlapping(all_helix_fragments[lex[i]], all_helix_fragments[lex[j]]) ){
					if(lex.at_end()){
						done=true;
						break;
					}
					lex.continue_at_dimension(j);
					i=1;
					j=i+1;
				}
			}
		}
		if (!done) {
			//check to see if each helix is involved in at least two pairs
			utility::vector1<core::Size> num_partners(bundle_size_,0);
			utility::vector1<FragmentPair> fragment_pairs;
			for(core::Size i=1; i<=bundle_size_; ++i){
				for(core::Size j=i+1; j<=bundle_size_; ++j){
					PairMap::const_iterator frag_pair_it = pair_map.find(std::make_pair(lex[i], lex[j]));
					if(frag_pair_it != pair_map.end()){
//						++bundle_counter;
						num_partners[i]++;
						num_partners[j]++;
						fragment_pairs.push_back(frag_pair_it->second);
					}
				}
			}
			bool valid = true;
			for(core::Size i=1; i<=bundle_size_; ++i){
				if(num_partners[i] < 2){
					valid = false;
				}
			}
			if(valid){
				//********************//
				//****Write Bundle****//
				//********************//
				++bundle_counter;
				//record the dSasa for each helix
//				record_helix_sasas(pose, cur_frag_set);
				
				bundle_insert_stmt.bind(1,bundle_counter);
				bundle_insert_stmt.bind(2,struct_id);
				bundle_insert_stmt.bind(3,bundle_size_);
				bundle_insert_stmt.bind(4,helix_size_);

//				bundle_insert_stmt.bind(1,struct_id);
//				bundle_insert_stmt.bind(2,bundle_size_);
//				bundle_insert_stmt.bind(3,helix_size_);
				basic::database::safely_write_to_database(bundle_insert_stmt);
//				//core::Size bundle_id(bundle_insert_stmt.sequence_last("helix_bundles_bundle_id_seq"));
				
				//*********************//
				//****Write Helices****//
				//*********************//
				std::map<HelicalFragmentOP, core::Size> helix_ids;
				for(core::Size i=1; i<=bundle_size_; ++i)
				{
					HelicalFragmentOP cur_fragment = all_helix_fragments[lex[i]];
					++helix_counter;
					helix_insert_stmt.bind(1,helix_counter);
					helix_insert_stmt.bind(2,bundle_counter);
					helix_insert_stmt.bind(3,struct_id);
					helix_insert_stmt.bind(4,cur_fragment->seq_start());
					helix_insert_stmt.bind(5,cur_fragment->seq_end());
					helix_insert_stmt.bind(6,cur_fragment->reversed());
					helix_insert_stmt.bind(7,cur_fragment->sasa());

//					helix_insert_stmt.bind(1,bundle_id);
//					helix_insert_stmt.bind(2,struct_id);
//					helix_insert_stmt.bind(3,cur_fragment->seq_start());
//					helix_insert_stmt.bind(4,cur_fragment->seq_end());
//					helix_insert_stmt.bind(5,cur_fragment->reversed());
//					helix_insert_stmt.bind(6,cur_fragment->sasa());
					basic::database::safely_write_to_database(helix_insert_stmt);
					
//					core::Size helix_id(bundle_insert_stmt.sequence_last("bundle_helices_helix_id_seq"));
					helix_ids.insert(std::make_pair(cur_fragment, helix_counter));
				}
				
				//********************//
				//****Write Pairs****//
				//********************//
				for(core::Size i=1; i<=fragment_pairs.size(); ++i){
					++pair_counter;
					
					FragmentPair const & cur_pair = fragment_pairs[i];
					
					core::Size frag_1_id = helix_ids[cur_pair.fragment_1];
					core::Size frag_2_id = helix_ids[cur_pair.fragment_2];
					
					pair_insert_stmt.bind(1,pair_counter);
					pair_insert_stmt.bind(2,struct_id);
					pair_insert_stmt.bind(3,bundle_counter);
					pair_insert_stmt.bind(4,frag_1_id);
					pair_insert_stmt.bind(5,frag_2_id);
					pair_insert_stmt.bind(6,cur_pair.fa_attr);
					pair_insert_stmt.bind(7,cur_pair.fa_fraction);
					pair_insert_stmt.bind(8,cur_pair.crossing_angle);
					pair_insert_stmt.bind(9,cur_pair.end_1_distance);
					pair_insert_stmt.bind(10,cur_pair.end_2_distance);
					
	//				pair_insert_stmt.bind(1,struct_id);
	//				pair_insert_stmt.bind(2,bundle_id);
	//				pair_insert_stmt.bind(3,frag_1_id);
	//				pair_insert_stmt.bind(4,frag_2_id);
	//				pair_insert_stmt.bind(5,cur_pair.fa_attr);
	//				pair_insert_stmt.bind(6,cur_pair.fa_fraction);
	//				pair_insert_stmt.bind(7,cur_pair.crossing_angle);
	//				pair_insert_stmt.bind(8,cur_pair.end_1_distance);
	//				pair_insert_stmt.bind(9,cur_pair.end_2_distance);
					basic::database::safely_write_to_database(pair_insert_stmt);
				}
			}
		}
	}
	TR << "Found " << bundle_counter << " total bundles" << std::endl;
	return 0;
}

bool
HelixBundleFeatures::overlapping(
	HelicalFragmentOP const & fragment_1,
	HelicalFragmentOP const & fragment_2
){
	if((fragment_1->seq_start() >= fragment_2->seq_start() &&
		fragment_1->seq_start() <= fragment_2->seq_end()) ||
		
	   (fragment_1->seq_end() >= fragment_2->seq_start() &&
		fragment_1->seq_end() <= fragment_2->seq_end()))
	{
		return true;
	}
	return false;
}
	
} //namespace helixAssembly
} //namespace features
} //namespace protocols
