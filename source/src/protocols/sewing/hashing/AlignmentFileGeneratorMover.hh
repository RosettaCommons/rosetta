// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/sewing/hashing/AlignmentFileGeneratorMover.hh
/// @brief Given a model file, edge file ,and one or more input structures, generates alignment files for use with AppendAssemblyMover.
/// @author guffysl (guffy@email.unc.edu)
/// @author minniel (minnie@email.unc.edu)

#ifndef INCLUDED_protocols_sewing_hashing_AlignmentFileGeneratorMover_hh
#define INCLUDED_protocols_sewing_hashing_AlignmentFileGeneratorMover_hh

// Unit headers
#include <protocols/sewing/hashing/AlignmentFileGeneratorMover.fwd.hh>
#include <protocols/moves/Mover.hh>
// Project Headers
#include <protocols/sewing/hashing/hasher_data.hh>
#include <protocols/sewing/hashing/BasisMapGenerator.fwd.hh>
#include <protocols/sewing/hashing/EdgeMapGenerator.fwd.hh>
#include <protocols/sewing/data_storage/LigandResidue.fwd.hh> //Need for LigandDescription
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
//#include <protocols/sewing/hashing/AlignmentGenerator.hh>


namespace protocols {
namespace sewing {
namespace hashing {



///@brief Given a model file, edge file ,and one or more input structures, generates alignment files for use with AppendAssemblyMover.
class AlignmentFileGeneratorMover : public protocols::moves::Mover {

public:

	AlignmentFileGeneratorMover();

	// copy constructor
	AlignmentFileGeneratorMover( AlignmentFileGeneratorMover const & src );


	AlignmentFileGeneratorMover(
		EdgeMapGeneratorOP edge_file_reader,
		std::string model_file_name,
		std::string alignment_file_name
	);

	AlignmentFileGeneratorMover(
		EdgeMapGeneratorOP edge_file_reader,
		std::string model_file_name
	);

	AlignmentFileGeneratorMover(
		std::string model_file_name,
		std::string edge_file_name,
		utility::vector1< core::Size > match_segments,
		utility::vector1< core::Size > pose_segment_starts,
		utility::vector1< core::Size > pose_segment_ends,
		core::Size recursive_depth
	);

	AlignmentFileGeneratorMover(
		BasisMapGeneratorOP bmg
	);


	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~AlignmentFileGeneratorMover();

	//Helper methods for constructors
	/*
	void initialize_from_options();
	void edge_file_from_options();
	void model_file_from_options();
	void alignment_settings_from_options();
	*/
	virtual void
	apply( core::pose::Pose & pose );
	//Each apply will have its own Hasher to avoid conflict b/w threads


public:
	virtual void
	show( std::ostream & output=std::cout ) const;

	std::string
	get_name() const;

	static std::string class_name();

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	//AlignmentFileGeneratorMover & operator=( AlignmentFileGeneratorMover const & src );

	//Getters
	BasisMapGeneratorOP basis_map_generator() const;

	std::string
	get_required_resnums() const;
	void
	set_required_resnums( std::string );

	core::select::residue_selector::ResidueSelectorCOP
	get_required_selector() const;

	void
	set_required_selector(  core::select::residue_selector::ResidueSelectorCOP );

	//Setters
	void basis_map_generator( BasisMapGeneratorOP bmg );

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	fresh_instance() const;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const;

	static core::Size add_pose_segments_to_segment_vector(
		core::pose::Pose const & pose,
		core::pose::PoseOP partner_pose,
		SegmentVectorCOP segvec,
		std::map< core::Size, data_storage::SmartSegmentOP > & pdb_segments,
		std::string pose_segment_starts_string,
		std::string pose_segment_ends_string,
		std::string pose_segment_dssp,
		utility::vector1< data_storage::LigandDescription > & ligands,
		std::map< core::Size, data_storage::LigandResidueCOP > & partner_ligands, //This will be filled in the function
		utility::vector1< data_storage::LigandDescription > & expanded_ligands,
		std::string required_resnums,
		core::select::residue_selector::ResidueSelectorCOP required_selector,
		bool strict_dssp_changes
	);
	static core::Size add_pose_segments_to_segment_vector(
		core::pose::Pose const & pose,
		core::pose::PoseOP partner_pose,
		BasisMapGeneratorOP bmg,
		utility::vector1< data_storage::LigandDescription > & ligands,
		std::map< core::Size, data_storage::LigandResidueCOP > & partner_ligands, //This will be filled in the function
		utility::vector1< data_storage::LigandDescription > & expanded_ligands,
		std::string required_resnums,
		core::select::residue_selector::ResidueSelectorCOP required_selector,
		bool set_segments_from_dssp
	);

	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static void append_ligands_subelement( utility::tag::XMLSchemaSimpleSubelementList & subs, utility::tag::XMLSchemaDefinition & xsd, bool include_coord );

	static void
	parse_ligands_tag(
		utility::tag::TagCOP ligands_tag,
		basic::datacache::DataMap const & data,
		utility::vector1< data_storage::LigandDescription > & ligands
	);

	static std::string
	ligands_subtag_ct_namer( std::string );

	static void
	alignment_settings_from_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap&,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &,
		BasisMapGeneratorOP bmg);

	static void
	parse_ligand_conformers(
		utility::vector1< data_storage::LigandDescription > & ligands
	);

	/*
	static std::map< core::Size, std::string >
	map_ligand_id_to_resnum_string(
	core::pose::Pose const & pose,
	utility::vector1< data_storage::LigandDescription > & ligands
	);
	*/

private:
	BasisMapGeneratorOP basis_map_generator_;
	utility::vector1< data_storage::LigandDescription > ligands_;
	std::string required_resnums_;
	core::select::residue_selector::ResidueSelectorCOP required_selector_;
	bool set_segments_from_dssp_;
};

std::ostream &operator<< (std::ostream &os, AlignmentFileGeneratorMover const &mover);


} //hashing
} //sewing
} //protocols


#endif //protocols_sewing_hashing_AlignmentFileGeneratorMover_hh
