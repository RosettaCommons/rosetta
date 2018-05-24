// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/movers/AppendAssemblyMover.hh
/// @brief an AssemblyMover for adding to existing poses
/// @author frankdt (frankdt@email.unc.edu)

#ifndef INCLUDED_protocols_sewing_movers_AppendAssemblyMover_hh
#define INCLUDED_protocols_sewing_movers_AppendAssemblyMover_hh

// Unit headers
#include <protocols/sewing/hashing/AlignmentFileGeneratorMover.hh>
#include <protocols/sewing/movers/AssemblyMover.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace sewing {
namespace movers {

///@brief an AssemblyMover for adding to existing poses
class AppendAssemblyMover : public AssemblyMover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	AppendAssemblyMover();

	/// @brief Copy constructor (not needed unless you need deep copies)
	AppendAssemblyMover( AppendAssemblyMover const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	virtual ~AppendAssemblyMover()=default;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

public:
	/// @brief Apply the mover
	virtual void
	apply( core::pose::Pose & pose ) override;

	virtual data_storage::SmartAssemblyOP
	set_up_assembly( core::pose::Pose &) override;

	virtual data_storage::SmartSegmentOP
	find_starting_segment( data_storage::SmartAssemblyOP assembly );

	/// @brief Show the contents of the Mover
	static std::string
	class_name();

	virtual void
	show( std::ostream & output = std::cout ) const override;

	/// @brief Get the name of the Mover
	virtual std::string
	get_name() const override;

	void
	set_partner_pdb( core::pose::PoseOP );

	void
	//set_required_starting_residues( utility::vector1< core::Size > );
	set_required_resnums( std::string );

	void
	set_required_selector( core::select::residue_selector::ResidueSelectorCOP );

	void
	set_segments_from_dssp( bool );

	bool
	get_extend_mode() const;

	void
	set_extend_mode(bool);
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

	//AppendAssemblyMover & operator=( AppendAssemblyMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	clone() const override;

	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & );

	static void
	attributes_for_append_assembly_mover( utility::tag::AttributeList & attributes );

	//Getters
	core::pose::PoseOP
	get_partner_pdb() const;

	//utility::vector1< core::Size >
	std::string
	get_required_resnums() const;

	core::select::residue_selector::ResidueSelectorCOP
	get_required_selector() const;

	bool
	get_segments_from_dssp() const;

	void
	set_modifiable_terminus(char);

	char
	get_modifiable_terminus() const;

	void
	set_start_node_vital_segments(std::string);

	std::string
	get_start_node_vital_segments() const;

	void
	set_output_partner(bool);

	bool
	get_output_partner() const;

	std::string
	get_pose_segment_starts_string() const;

	std::string
	get_pose_segment_ends_string() const;

	void
	set_pose_segment_starts_string( std::string const & );

	void
	set_pose_segment_ends_string( std::string const & );

	void set_pose_segment_dssp( std::string const & );

	std::string
	get_pose_segment_dssp() const;

	bool
	get_strict_dssp_changes() const;

	void
	set_strict_dssp_changes( bool );

	void
	set_ligands( utility::vector1< data_storage::LigandDescription > );

	utility::vector1< data_storage::LigandDescription >
	get_ligands() const;

	utility::vector1< data_storage::LigandDescription > &
	get_nonconst_ligands();


private: // data

	core::pose::PoseOP partner_pdb_;
	std::string required_resnums_;
	core::select::residue_selector::ResidueSelectorCOP required_selector_;
	utility::vector1< data_storage::LigandDescription > ligands_;
	bool segments_from_dssp_ = true;
	char modifiable_terminus_ = 'B';
	std::string start_node_vital_segments_;
	bool output_partner_ = true;
	std::string  pose_segment_starts_string_;
	std::string pose_segment_ends_string_;
	std::string pose_segment_dssp_;
	bool strict_dssp_changes_ = false;
	bool extend_mode_ = false;

};
std::ostream &
operator<<( std::ostream & os, AppendAssemblyMover const & mover );

} //protocols
} //sewing
} //movers

#endif //protocols/sewing/movers_AppendAssemblyMover_hh
