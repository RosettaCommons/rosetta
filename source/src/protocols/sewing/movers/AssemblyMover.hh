// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/movers/AssemblyMover.hh
/// @brief an interface for making Movers that deal with Assemblies
/// @author frankdt (frankdt@email.unc.edu)

#ifndef INCLUDED_protocols_sewing_movers_AssemblyMover_hh
#define INCLUDED_protocols_sewing_movers_AssemblyMover_hh

// Unit headers
#include <protocols/sewing/movers/AssemblyMover.fwd.hh>
#include <protocols/moves/Mover.hh>

//Project headers
#include <protocols/sewing/data_storage/SmartAssembly.hh>
#include <protocols/sewing/data_storage/HashedSmartAssembly.hh>
#include <protocols/sewing/hashing/BasisMapGenerator.hh>
#include <protocols/sewing/hashing/EdgeMapGenerator.hh>
#include <protocols/sewing/scoring/AssemblyScorer.hh>
#include <protocols/sewing/requirements/AssemblyRequirement.hh>
#include <protocols/sewing/hashing/ModelFileReader.hh>
#include <core/pose/Pose.hh>
// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <numeric/random/random.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
namespace protocols {
namespace sewing {
namespace movers {

///@brief an interface for making Movers that deal with Assemblies

// typedef std::map< std::string,  std::pair<core::Real, scoring::AssemblyScorerOP> > ScorerMap;

class AssemblyMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	AssemblyMover();

	/// @brief Copy constructor (not needed unless you need deep copies)
	AssemblyMover( AssemblyMover const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	virtual ~AssemblyMover()=default;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

public:
	/// @brief Apply the mover
	virtual void
	apply( core::pose::Pose & pose ) override;

	virtual void
	set_up_hashing();

	virtual void
	generate_assembly(data_storage::SmartAssemblyOP, core::pose::Pose &);


	/// @brief Show the contents of the Mover
	static std::string
	class_name();

	virtual void
	show( std::ostream & output = std::cout ) const override;

	/// @brief Get the name of the Mover
	virtual std::string
	get_name() const override;

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	virtual void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	//AssemblyMover & operator=( AssemblyMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	clone() const override;

	///@brief return the weighted sum of all score terms for the assembly.
	virtual core::Real
	score( data_storage::SmartAssemblyOP );

	///@brief resets scorers to their score before the move that is now being reverted
	void
	revert_score();

	///@brief add an AssemblyScorer and its weight to the ScorerMap
	virtual core::Size
	add_scorer( scoring::AssemblyScorerOP scorer );
	// add_scorer( scoring::AssemblyScorerOP scorer, core::Real weight );

	///@brief get the unweighted value for a particular score term
	virtual core::Real
	get_unweighted_score_term( data_storage::SmartAssemblyOP, core::Size scorer_index );
	// get_unweighted_score_term( data_storage::SmartAssemblyOP, std::string scoretype );

	virtual void
	add_motif_scorers_to_score_file( core::pose::Pose & pose, data_storage::SmartAssemblyOP & assembly );

	virtual void
	print_statistics( data_storage::SmartAssemblyOP ) const;

	virtual data_storage::SmartAssemblyOP
	set_up_assembly( core::pose::Pose & pose);

	bool
	assembly_meets_requirements( data_storage::SmartAssemblyOP assembly );


	//Helper functions for parse_my_tag

	///@brief Parses the Requirements subtag of AssemblyMover.
	void
	parse_requirements(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap& datamap,
		protocols::filters::Filters_map const & filtermap,
		protocols::moves::Movers_map const & movermap,
		core::pose::Pose const & pose);

	///@brief Parses the AssemblyScorers subtag of AssemblyMover
	void
	parse_assembly_scorers(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap& datamap,
		protocols::filters::Filters_map const & filtermap,
		protocols::moves::Movers_map const & movermap,
		core::pose::Pose const & pose);

	///@brief Set default score function for this AssemblyMover. Called if no AssemblyScorers are specified through RosettaScripts.
	virtual void
	set_default_assembly_scorers();


	///@brief Set default requirements for this AssemblyMover. Called if no requirements are specified through RosettaScripts.
	virtual void
	set_default_requirements();

	void
	add_requirement( requirements::AssemblyRequirementOP requirement );


	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & );


	static std::string
	assembly_mover_subtag_ct_namer( std::string );

	static std::string mover_name();

	static void
	define_generic_assembly_mover_attributes(utility::tag::AttributeList & );

	//Getters
	core::Size
	get_window_width() const;

	hashing::SegmentVectorCOP
	get_segment_vector() const;

	hashing::BasisMapGeneratorOP
	get_basis_map_generator() const;

	hashing::EdgeMapGeneratorOP
	get_edge_map_generator() const;

	std::string
	get_model_file_name() const;

	std::string
	get_edge_file_name() const;

	// ScorerMap
	utility::vector1< scoring::AssemblyScorerOP >
	get_assembly_scorers() const;

	utility::vector1< requirements::AssemblyRequirementOP >
	get_requirements() const;

	bool
	get_hashed() const;

	core::Real
	get_start_temperature() const;

	core::Real
	get_end_temperature() const;

	core::Size
	get_max_num_segments() const;

	core::Size
	get_max_segment_length() const;

	core::Size
	get_min_cycles() const;

	core::Size
	get_max_cycles() const;

	core::Real
	get_add_probability() const;

	core::Real
	get_delete_probability() const;

	core::Real
	get_conformer_switch_probability() const;

	bool
	get_output_pose_per_move() const;

	core::Size
	get_output_pose_count() const;

	core::Real
	get_lowest_score() const;

	//Setters
	void
	set_window_width( core::Size width );

	void
	set_segment_vector( hashing::SegmentVectorCOP segvec );

	void
	set_basis_map_generator( hashing::BasisMapGeneratorOP bmg );

	void
	set_edge_map_generator( hashing::EdgeMapGeneratorOP emg );

	void
	set_model_file_name( std::string models );

	void
	set_edge_file_name( std::string edges );

	void
	set_assembly_scorers( utility::vector1< scoring::AssemblyScorerOP > assembly_scorers );
	// set_assembly_scorers( ScorerMap assembly_scorers );

	void
	set_requirements( utility::vector1 <requirements::AssemblyRequirementOP > requirements );

	void
	set_hashed( bool hashed );

	void
	set_start_temperature( core::Real temp );

	void
	set_end_temperature( core::Real temp );

	void
	set_max_num_segments( core::Size max_segs );

	void
	set_max_segment_length( core::Size max_seg_length );

	void
	set_min_cycles( core::Size min_cycles );

	void
	set_max_cycles( core::Size max_cycles );

	void
	set_add_probability( core::Real add_probability );

	void
	set_delete_probability( core::Real delete_probability );

	void
	set_conformer_switch_probability( core::Real prob );

	void
	set_output_pose_per_move( bool output_pose_per_move);

	void
	output_pose( data_storage::SmartAssemblyOP assembly , core::pose::Pose &);

	void
	set_lowest_score ( core::Real lowest_score );

	void
	set_recover_lowest_assembly(bool);

	bool
	get_recover_lowest_assembly() const;


private:
	hashing::EdgeMapGeneratorOP edge_map_generator_;
	hashing::BasisMapGeneratorOP basis_map_generator_;
	std::string model_file_name_;
	std::string edge_file_name_;
	// //Map of scorers
	// ScorerMap assembly_scorers_;
	utility::vector1< scoring::AssemblyScorerOP > assembly_scorers_;
	//Requirements (still a RequirementSet? Or maybe a map like the Scorers?)
	utility::vector1< requirements::AssemblyRequirementOP > requirements_;
	core::Size max_num_segments_;
	core::Size max_segment_length_;
	core::Size min_cycles_;
	core::Size max_cycles_;
	core::Real current_add_probability_;
	core::Real add_probability_;
	core::Real delete_probability_;
	core::Real conformer_switch_probability_;
	core::Real random_action_;
	core::Real score_;
	core::Real old_score_;
	core::Real lowest_score_;
	core::Real start_temperature_;
	core::Real end_temperature_;
	bool output_pose_per_move_;
	core::Size output_pose_count_;
	bool hashed_;
	hashing::SegmentVectorCOP segment_vector_;
	utility::vector1< data_storage::BasisPair > lowest_scoring_assembly_;
	data_storage::SmartAssemblyOP best_assembly_;
	core::Size window_width_;
	bool recover_lowest_assembly_;
};

std::ostream &
operator<<( std::ostream & os, AssemblyMover const & mover );

} //protocols
} //sewing
} //movers

#endif //protocols/sewing/movers_AssemblyMover_hh
