// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :notabs=false:tabSize=4:indentsize=4:
//
// (c) copyright rosetta commons member institutions.
// (c) this file is part of the rosetta software suite and is made available under license.
// (c) the rosetta software is developed by the contributing members of the rosetta commons.
// (c) for more information, see http://www.rosettacommons.org. questions about this can be
// (c) addressed to university of washington uw techtransfer, email: license@u.washington.edu.

/// @file protocols/features/helixAssembly/ModelFeatures.cc
/// @brief Search through a pose for sets of 3 helices and generate RMSD calculations between all pairs of them
/// @author Tim Jacobs

//Unit Headers
#include <protocols/features/ModelFeatures.hh>
#include <protocols/features/ModelFeaturesCreator.hh>

//Protocol Headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/PoseMetricCalculatorBase.fwd.hh>
#include <core/select/util/interface_vector_calculate.hh>
#include <core/pose/PDBInfo.hh>
#include <core/graph/Graph.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>

#include <protocols/toolbox/pose_metric_calculators/NeighborhoodByDistanceCalculator.hh>

//Utility
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/tag/Tag.hh>
#include <utility/LexicographicalIterator.hh>

//Utility Headers
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/Constraint.hh>
#include <basic/database/schema_generator/DbDataType.hh>
#include <basic/database/sql_utils.hh>
#include <basic/MetricValue.hh>
#include <basic/options/util.hh>
#include <basic/Tracer.hh>

#include <numeric/PCA.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>

//C++
#include <string>
#include <math.h>

//External Headers
#include <cppdb/frontend.h>
#include <boost/graph/bron_kerbosch_all_cliques.hpp>

static basic::Tracer TR("devel.sewing.ModelFeatures");

namespace protocols {
namespace features {

///////////////////////////////////////////////////////////
/////////////  Creator functions   ////////////////////////
///////////////////////////////////////////////////////////

ModelFeaturesCreator::ModelFeaturesCreator() {}
ModelFeaturesCreator::~ModelFeaturesCreator() {}
protocols::features::FeaturesReporterOP ModelFeaturesCreator::create_features_reporter() const {
	return ModelFeaturesOP( new ModelFeatures );
}

std::string ModelFeaturesCreator::type_name() const {
	return "ModelFeatures";
}

///////////////////////////////////////////////////////////
///////////////  Class functions   ////////////////////////
///////////////////////////////////////////////////////////

ModelFeatures::ModelFeatures():
	min_helix_size_(12),
	min_strand_size_(6),
	min_ss_cluster_size_(3),
	max_ss_cluster_size_(5),
	distance_cutoff_(14.0),
	angle_cutoff_(45.0)
{}

utility::vector1<std::string>
ModelFeatures::features_reporter_dependencies() const
{
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	dependencies.push_back("ProteinResidueConformationFeatures");
	dependencies.push_back("ResidueSecondaryStructureFeatures");
	dependencies.push_back("SecondaryStructureSegmentFeatures");
	return dependencies;
}

void
ModelFeatures::write_schema_to_db(utility::sql_database::sessionOP db_session) const
{
	using namespace basic::database::schema_generator;

	////////////////////////////////////
	///////Sewing Models Table//////////
	////////////////////////////////////
	Column struct_id("struct_id", DbDataTypeOP (new DbBigInt()), false /*not null*/, false /*don't autoincrement*/);
	Column model_id("model_id", DbDataTypeOP(new DbInteger()), false /*not null*/, true /*autoincrement*/);
	Column num_segments("num_segments", DbDataTypeOP(new DbInteger()), false /*not null*/, false /*autoincrement*/);

	Schema sewing_models("sewing_models", model_id);
	sewing_models.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", false /*defer*/));
	sewing_models.add_column(num_segments);
	sewing_models.write(db_session);

	////////////////////////////////////
	//////Model Segments Table//////////
	////////////////////////////////////
	Column model_id_fkey("model_id", DbDataTypeOP(new DbInteger()), false /*not null*/, false /*autoincrement*/);
	Column segment_id_fkey("segment_id", DbDataTypeOP(new DbInteger()), false /*not null*/, false /*autoincrement*/);

	utility::vector1<Column> pkey_cols;
	pkey_cols.push_back(model_id_fkey);
	pkey_cols.push_back(segment_id_fkey);

	Schema model_segments("model_segments", pkey_cols);
	model_segments.add_foreign_key(ForeignKey(model_id_fkey, "sewing_models", "model_id", false /*defer*/));
	model_segments.add_foreign_key(ForeignKey(segment_id_fkey, "secondary_structure_segments", "segment_id", false /*defer*/));
	model_segments.write(db_session);
}


utility::vector1<Segment>
ModelFeatures::find_segments(
	protocols::features::StructureID struct_id,
	utility::sql_database::sessionOP db_session
) {
	using namespace cppdb;

	std::string select_string =
		"SELECT\n"
		"\tsegment_id,\n"
		" dssp, \n"
		"\tresidue_begin,\n"
		"\tresidue_end\n"
		"FROM\n"
		"\tsecondary_structure_segments\n"
		"WHERE\n"
		"\tstruct_id = ?\n AND "
		"\tdssp IN('H','E')\n"
		"ORDER BY residue_begin;";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	result res(basic::database::safely_read_from_database(select_statement));

	utility::vector1<Segment> segments;
	while ( res.next() )
			{
		std::string dssp;
		Size segment_id, residue_begin, residue_end;
		res >> segment_id >> dssp >> residue_begin >> residue_end;

		Segment segment;
		segment.id = segment_id;
		segment.dssp = dssp;
		segment.begin = residue_begin;
		segment.end = residue_end;

		if ( dssp == "H" && residue_end - residue_begin + 1 >= min_helix_size_ ) {
			segments.push_back(segment);
		} else if ( dssp == "E" && residue_end - residue_begin + 1 >= min_strand_size_ ) {
			segments.push_back(segment);
		}
	}

	TR << "Found " << segments.size() << " secondary structure segments." << std::endl;
	return segments;
}

///@details
core::Size
ModelFeatures::report_features(
	core::pose::Pose const & pose,
	utility::vector1<bool> const &,
	protocols::features::StructureID struct_id,
	utility::sql_database::sessionOP db_session
){
	std::string pdb_name = utility::string_split(pose.pdb_info()->name(), '/').back();
	utility::trim( pdb_name, ".pdb");

	utility::vector1<Segment> segments =
		find_segments(struct_id, db_session);

	SegmentGraph graph;
	std::map< Segment, SegmentVertex > vertices;
	for ( core::Size seg_i=1; seg_i<=segments.size(); ++seg_i ) {

		//Map this segment to a vertex (this seems ridiculous to have to do)
		SegmentVertex v = graph.add_vertex(segments[seg_i].id);
		vertices[segments[seg_i]] = v;

		utility::vector1< numeric::xyzVector< core::Real> > bb_coords;
		for ( core::Size res_i=segments[seg_i].begin; res_i <= segments[seg_i].end; ++res_i ) {
			bb_coords.push_back(pose.residue(res_i).atom("CA").xyz());
		}
		segments[seg_i].first_pc = numeric::first_principal_component(bb_coords);
		TR << "Segment " << segments[seg_i].id << " pc: " << segments[seg_i].first_pc << std::endl;
	}

	std::map< core::Size, Segment > residue_segments;
	for ( core::Size seg_i=1; seg_i<=segments.size(); ++seg_i ) {
		for ( core::Size res_i=segments[seg_i].begin; res_i <= segments[seg_i].end; ++res_i ) {
			residue_segments[res_i] = segments[seg_i];
		}
	}


	//iterate through every residue to find unique sets of secondary structure clusters
	std::set< std::set< Segment > > ss_groups;
	for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
		std::set<core::Size> central_residues;
		central_residues.insert(i);

		std::string const nb_calc("neighbor_calculator_" + utility::to_string(i));
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( nb_calc,
			core::pose::metrics::PoseMetricCalculatorOP (new protocols::toolbox::pose_metric_calculators::NeighborhoodByDistanceCalculator( central_residues, distance_cutoff_ ) ) );

		//Get the neighbors of the current residue
		basic::MetricValue< std::set< core::Size > > neighbor_mv;
		pose.metric( nb_calc, "neighbors", neighbor_mv);
		std::set<core::Size> neighbors = neighbor_mv.value();
		core::pose::metrics::CalculatorFactory::Instance().remove_calculator(nb_calc);

		std::set<Segment> ss_group;
		for ( std::set<core::Size>::const_iterator it = neighbors.begin(); it != neighbors.end(); ++it ) {
			if ( residue_segments.find(*it) != residue_segments.end() ) {
				ss_group.insert(residue_segments.find(*it)->second);
			}
		}
		if ( ss_group.size() >= 3 ) {
			ss_groups.insert(ss_group);
		}
	}

	TR << "Found " << ss_groups.size() << " unique secondary structure groups of at least 3 segments." << std::endl;
	core::pose::Pose score_pose(pose);
	for ( std::set< std::set< Segment > >::iterator groups_it = ss_groups.begin(); groups_it != ss_groups.end(); ++groups_it ) {
		std::set< Segment > current_group = *groups_it;

		core::scoring::ScoreFunctionOP scorefxn ( new core::scoring::ScoreFunction );
		scorefxn->set_weight(core::scoring::fa_atr, 1.0);
		scorefxn->score(score_pose);

		core::scoring::EnergyGraph const & eg =
			score_pose.energies().energy_graph();

		for ( std::set< Segment >::const_iterator it1 = current_group.begin(); it1 != current_group.end(); ++it1 ) {
			for ( std::set< Segment >::const_iterator it2 = it1; it2!= current_group.end(); ++it2 ) {
				if ( it1 == it2 ) { continue; }

				//Principal component angle (crude check for parallelness)
				core::Real pc_angle = numeric::angle_degrees(it1->first_pc, numeric::xyzVector<core::Real>(0.0,0.0,0.0), it2->first_pc);
				TR << "PC angle for segments " << it1->id << " and " << it2->id << ": " << pc_angle << std::endl;

				if ( pc_angle > angle_cutoff_ && pc_angle < (180-angle_cutoff_) ) {
					continue;
				}

				//Cap distance
				core::Size seg_1_atr_start = 100000000;
				core::Size seg_1_atr_end = 0;
				core::Size seg_2_atr_start = 100000000;
				core::Size seg_2_atr_end = 0;
				for ( core::Size ii = it1->begin; ii <= it1->end; ++ii ) {
					for ( core::Size jj = it2->begin; jj <= it2->end; ++jj ) {
						core::Real distance_sq = pose.residue(ii).atom("CA").xyz().distance_squared(pose.residue(jj).atom("CA").xyz());
						if ( distance_sq < (distance_cutoff_*distance_cutoff_) ) {
							seg_1_atr_start = std::min(seg_1_atr_start, ii);
							seg_1_atr_end = std::max(seg_1_atr_end, ii);

							seg_2_atr_start = std::min(seg_2_atr_start, jj);
							seg_2_atr_end = std::max(seg_2_atr_end, jj);
						}
					}
				}
				if ( it1->dssp == "H"  && seg_1_atr_end-seg_1_atr_start+1 < min_helix_size_ ) {
					continue;
				} else if ( it1->dssp == "E"  && seg_2_atr_end-seg_2_atr_start+1 < min_strand_size_ ) {
					continue;
				}

				bool attraction = false;
				for ( core::Size ii = it1->begin; ii <= it1->end; ++ii ) {
					for ( core::graph::Node::EdgeListConstIter
							iter = eg.get_node(ii)->const_upper_edge_list_begin(),
							iter_end = eg.get_node(ii)->const_upper_edge_list_end();
							iter != iter_end; ++iter ) {
						core::Size jj = (*iter)->get_other_ind(ii);
						if ( jj >= it2->begin && jj <= it2->end ) {
							const core::scoring::EnergyEdge* ee =
								dynamic_cast< const core::scoring::EnergyEdge* > (*iter);

							core::Real atr_energy = (*ee)[core::scoring::fa_atr];
							if ( atr_energy < 0 ) {
								//TR << "Attractive energy between " << *it1 << " and " << *it2 << " residues " << ii << " and " << jj << ": " << atr_energy << std::endl;
								attraction = true;
								break;
							}
						}
					}
				}
				if ( attraction ) {
					graph.add_edge(vertices[*it1], vertices[*it2]);
				}
			}
		}

	}

	// Instantiate the visitor for printing cliques
	std::set< std::set<core::Size> > * cliques = new std::set< std::set<core::Size> >;
	clique_saver vis(cliques);

	// Use the Bron-Kerbosch algorithm to find all cliques, printing them
	// as they are found.
	bron_kerbosch_all_cliques(graph, vis);
	TR << "Found " << cliques->size() << " cliques" << std::endl;
	core::Size clique_counter = 0;
	for ( std::set< std::set<core::Size> >::const_iterator it= cliques->begin(); it != cliques->end(); ++it ) {
		if ( it->size() < min_ss_cluster_size_ || it->size() > max_ss_cluster_size_ ) {
			continue;
		}
		++clique_counter;
		std::set<core::Size> current_clique = *it;
		//
		//  core::pose::Pose clique_pose = pose;
		//  std::set<core::Size> resnums;
		//  for(core::Size i=1; i<=segments.size(); ++i) {
		//   if(current_clique.find(segments[i].id) != current_clique.end()) {
		//    for(core::Size j=segments[i].begin; j<=segments[i].end; ++j) {
		//     resnums.insert(j);
		//    }
		//   }
		//  }
		//  trim_pose(clique_pose, resnums);
		//  clique_pose.dump_pdb(pdb_name + "_clique_" + utility::to_string(clique_counter) + ".pdb");
		write_clique_to_db(current_clique, struct_id, db_session);
	}
	delete cliques;

	return 0;
}

void
ModelFeatures::write_clique_to_db(
	std::set<core::Size> const & clique,
	protocols::features::StructureID struct_id,
	utility::sql_database::sessionOP db_session
) const {

	using namespace cppdb;

	std::string model_insert =
		"INSERT INTO sewing_models (struct_id, num_segments)  VALUES (?,?);";
	statement model_insert_stmt(basic::database::safely_prepare_statement(model_insert,db_session));
	model_insert_stmt.bind(1, struct_id);
	model_insert_stmt.bind(2, clique.size());
	basic::database::safely_write_to_database( model_insert_stmt );
	core::Size model_id = model_insert_stmt.sequence_last("sewing_models_model_id_seq");

	std::string model_segment_insert =
		"INSERT INTO model_segments (model_id, segment_id)  VALUES (?,?);";
	statement model_segment_insert_stmt(basic::database::safely_prepare_statement(model_segment_insert,db_session));
	for ( std::set<core::Size>::const_iterator it=clique.begin(); it != clique.end(); ++it ) {
		model_segment_insert_stmt.bind(1, model_id);
		model_segment_insert_stmt.bind(2, *it);
		basic::database::safely_write_to_database( model_segment_insert_stmt );
	}

}

void
ModelFeatures::trim_pose(
	core::pose::Pose & pose,
	std::set<core::Size> resnums
) const {
	core::Size num_removed_residues=0;
	core::Size total_res = pose.total_residue();
	for ( core::Size i=1; i<=total_res; ++i ) {
		if ( resnums.find(i) == resnums.end() ) {
			pose.conformation().delete_residue_slow(i-num_removed_residues);
			++num_removed_residues;
		}
	}
}


void
ModelFeatures::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/
){
	runtime_assert(tag->getOption<std::string>("name") == type_name());

	if ( tag->hasOption("min_helix_size") ) {
		min_helix_size_ = tag->getOption<core::Size>("min_helix_size");
	}

	if ( tag->hasOption("min_strand_size") ) {
		min_strand_size_ = tag->getOption<core::Size>("min_strand_size");
	}

	if ( tag->hasOption("min_ss_cluster_size") ) {
		min_ss_cluster_size_ = tag->getOption<core::Size>("min_ss_cluster_size");
	}

	if ( tag->hasOption("max_ss_cluster_size") ) {
		max_ss_cluster_size_ = tag->getOption<core::Size>("max_ss_cluster_size");
	}

	if ( tag->hasOption("distance_cutoff") ) {
		distance_cutoff_ = tag->getOption<core::Real>("distance_cutoff");
	}

	if ( tag->hasOption("angle_cutoff") ) {
		angle_cutoff_ = tag->getOption<core::Real>("angle_cutoff");
	}
}

} //namespace sewing
} //namespace devel
