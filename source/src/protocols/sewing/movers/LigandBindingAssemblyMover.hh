// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/movers/LigandBindingAssemblyMover.hh
/// @brief an AssemblyMover for adding to existing poses
/// @author frankdt (frankdt@email.unc.edu)

#ifndef INCLUDED_protocols_sewing_movers_LigandBindingAssemblyMover_hh
#define INCLUDED_protocols_sewing_movers_LigandBindingAssemblyMover_hh

// Unit headers
//#include <protocols/sewing/hashing/AlignmentFileGeneratorMover.hh>
#include <protocols/sewing/hashing/LigandBindingResPlacer.hh>
#include <protocols/sewing/movers/AppendAssemblyMover.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace sewing {
namespace movers {

///@brief an AssemblyMover for adding to existing poses
class LigandBindingAssemblyMover : public AppendAssemblyMover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	LigandBindingAssemblyMover();

	/// @brief Copy constructor (not needed unless you need deep copies)
	LigandBindingAssemblyMover( LigandBindingAssemblyMover const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	virtual ~LigandBindingAssemblyMover()=default;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

public:
	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	data_storage::SmartAssemblyOP
	set_up_assembly(core::pose::Pose & pose) override;

	void
	generate_assembly(data_storage::SmartAssemblyOP, core::pose::Pose &) override;

	/// @brief Show the contents of the Mover
	static std::string
	class_name() ;

	void
	show( std::ostream & output = std::cout ) const override;

	/// @brief Get the name of the Mover
	std::string
	get_name() const override;

	core::Real
	get_min_distance( data_storage::SmartSegmentOP, data_storage::LigandResidueOP );

	void
	add_motif_scorers_to_score_file( core::pose::Pose & pose, data_storage::SmartAssemblyOP & assembly ) override;

	//Getters

	//std::map< std::string, core::Real >
	//get_geometry_score_thresholds() const;

	core::Real
	get_distance_cutoff() const;

	core::Size
	get_segment_distance_cutoff() const;

	core::Size
	get_binding_cycles() const;

	utility::vector1< requirements::AssemblyRequirementOP >
	get_ligand_requirements() const;

	utility::vector1< requirements::AssemblyRequirementOP >
	get_non_ligand_requirements() const;

	bool
	get_build_site_only() const;

	//std::map< std::string, utility::pointer::shared_ptr< std::map< char, hashing::LigandHashMap > > >
	//get_ligand_coords() const;


	//Setters
	//void
	//set_ligand_coords( std::map< std::string, utility::pointer::shared_ptr< std::map< char, hashing::LigandHashMap > > > coords );

	//void
	//set_geometry_score_thresholds( std::map< std::string, core::Real > );

	void
	set_distance_cutoff( core::Real cutoff );

	void
	set_segment_distance_cutoff( core::Size cutoff );

	void
	set_binding_cycles( core::Size cycles );

	void
	set_ligand_requirements( utility::vector1 <requirements::AssemblyRequirementOP > requirements );

	void
	set_non_ligand_requirements( utility::vector1 <requirements::AssemblyRequirementOP > requirements );

	void
	set_build_site_only( bool build_only );
	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	//LigandBindingAssemblyMover & operator=( LigandBindingAssemblyMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & );


private: // data


	//These data are exclusive to LigandBindingAssemblyMover (not needed in context of AppendAssemblyMover)
	//Adjust parse_my_tag and provide_xml_schema to reflect this fact


	core::Real distance_cutoff_;
	core::Size segment_distance_cutoff_;
	core::Size binding_cycles_;
	//This is flushed every time we call generate_assembly, so no need to include in copy constructor
	std::map< core::Size, hashing::LigandBindingResPlacerOP > binder_finders_;
	utility::vector1< data_storage::LigandDescription > expanded_ligands_;
	utility::vector1< requirements::AssemblyRequirementOP > ligand_requirements_;
	utility::vector1< requirements::AssemblyRequirementOP > non_ligand_requirements_;
	bool build_site_only_ = false;

};

std::ostream &
operator<<( std::ostream & os, LigandBindingAssemblyMover const & mover );

} //protocols
} //sewing
} //movers

#endif //protocols/sewing/movers_LigandBindingAssemblyMover_hh
