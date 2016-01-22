// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :notabs=false:tabSize=4:indentsize=4:
//
// (c) copyright rosetta commons member institutions.
// (c) this file is part of the rosetta software suite and is made available under license.
// (c) the rosetta software is developed by the contributing members of the rosetta commons.
// (c) for more information, see http://www.rosettacommons.org. questions about this can be
// (c) addressed to university of washington uw techtransfer, email: license@u.washington.edu.

/// @file ModelFeatures.hh
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_features_ModelFeatures_hh
#define INCLUDED_protocols_features_ModelFeatures_hh

//Unit
#include <protocols/features/ModelFeatures.fwd.hh>

//Protocols
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/simple_filters/InterfaceSasaFilter.hh>

//Devel
#include <protocols/features/helixAssembly/HelicalFragment.hh>

//Utility
#include <utility/vector1.hh>

//External
#include <boost/graph/undirected_graph.hpp>

namespace protocols {
namespace features {

typedef boost::undirected_graph<core::Size> SegmentGraph;
//typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS> SegmentGraph;
typedef boost::graph_traits<SegmentGraph>::vertex_descriptor SegmentVertex;
typedef boost::graph_traits<SegmentGraph>::edge_descriptor InteractionEdge;

struct clique_saver
{
	clique_saver(std::set< std::set<core::Size> >* cliques):
		cliques_(cliques)
	{ }

	template <typename Clique, typename SegmentGraph>
	void clique(const Clique& c, const SegmentGraph& g) {
		std::set<core::Size> clique;
		typename Clique::const_iterator i, end = c.end();
		for ( i = c.begin(); i != end; ++i ) {
			clique.insert(g[*i]);
		}
		cliques_->insert(clique);
	}
	std::set< std::set<core::Size> >* cliques_;
};

struct Segment
{
	core::Size id;
	core::Size begin;
	core::Size end;
	std::string dssp;

	numeric::xyzVector< core::Real > first_pc;

	friend
	bool
	operator<(Segment const & a, Segment const &b) {
		return a.id < b.id;
	}
};

class ModelFeatures : public protocols::features::FeaturesReporter
{

public:

	ModelFeatures();

	virtual
	std::string
	type_name() const  {
		return "ModelFeatures";
	};


	/// @brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const;

	/// @brief generate the table schemas and write them to the database
	virtual
	void
	write_schema_to_db(
		utility::sql_database::sessionOP db_session
	) const;

	/// @brief get the vector of beginnings/ends of secondary structure segments
	utility::vector1<Segment>
	find_segments(
		protocols::features::StructureID struct_id,
		utility::sql_database::sessionOP db_session
	);

	/// @brief for testing purposes only
	void
	trim_pose(
		core::pose::Pose & pose,
		std::set<core::Size> resnums
	) const;

	/// @brief collect all the feature data for the pose
	virtual
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1<bool> const & relevant_residues,
		protocols::features::StructureID struct_id,
		utility::sql_database::sessionOP db_session
	);

	void
	write_clique_to_db(
		std::set<core::Size> const & clique,
		protocols::features::StructureID struct_id,
		utility::sql_database::sessionOP db_session
	) const;


	virtual
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & pose
	);

private:

	core::Size min_helix_size_;
	core::Size min_strand_size_;

	core::Size min_ss_cluster_size_;
	core::Size max_ss_cluster_size_;

	core::Real distance_cutoff_;
	core::Real angle_cutoff_;

};

} //namespace sewing
} //namespace devel

#endif
