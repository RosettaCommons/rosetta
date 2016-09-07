// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/enzdes/DesignProteinLigandInterface.hh
/// @brief a collection of Task operations used in enzyme design
/// @author Florian Richter, Sinisa Bjelic (sbjelic@u.washington.edu), Rocco Moretti (rmoretti@u.washington.edu)


#ifndef INCLUDED_protocols_enzdes_EnzdesTaskOperations_hh
#define INCLUDED_protocols_enzdes_EnzdesTaskOperations_hh

#include <protocols/enzdes/EnzdesTaskOperations.fwd.hh>
#include <protocols/motifs/LigandMotifSearch.fwd.hh>

#include <core/pack/task/operation/TaskOperation.hh>
#include <core/types.hh>

#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>
#include <set>
#include <string>


namespace protocols {
namespace enzdes {


/// @brief queries the pose cst cache for the catalytic residues,
/// and sets the behavior for them as specified in the cstfile
class SetCatalyticResPackBehavior : public core::pack::task::operation::TaskOperation
{

public:
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
	typedef core::pose::Pose Pose;
	typedef core::pack::task::PackerTask PackerTask;

public:

	SetCatalyticResPackBehavior();
	~SetCatalyticResPackBehavior() override;

	TaskOperationOP clone() const override;

	/// @brief Change a packer task in some way.  The input pose is the one to which the input
	/// task will be later applied.
	void apply( Pose const & pose, PackerTask & task ) const override;

	void parse_tag( TagCOP tag , DataMap & ) override;
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "SetCatalyticResPackBehavior"; }

	void
	set_fix_catalytic_aa( bool setting ){
		fix_catalytic_aa_ = setting; }

private:
	bool fix_catalytic_aa_;
	std::string behavior_non_catalytic_; //this string defaults to empty, but if set through the tag, the noncatalytic residues will be changed accordingly

};


/// @brief Given a set of cut1/cut2/cut3/cut4 distance specifications, alter a
///packing task to set residues within alpha carbons within cut1 of a ligand
///(or within cut2 with beta carbons facing inwards) to redesign, and within cut3
///(or cut4 facing inwards) to repack, and all others to fixed. If a resfile is provided,
///only do the detection for those residues set to AUTO in the resfile.
class DetectProteinLigandInterface: public core::pack::task::operation::TaskOperation
{

public:
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
	typedef core::pose::Pose Pose;
	typedef core::pack::task::PackerTask PackerTask;

public:

	DetectProteinLigandInterface();
	~DetectProteinLigandInterface() override;

	TaskOperationOP clone() const override;

	/// @brief Initialize the class based on the command line options.
	void init_from_options();

	/// @brief Change a packer task in some way.  The input pose is the one to which the input
	/// task will be later applied.
	void apply( Pose const & pose, PackerTask & task) const override;

	void parse_tag( TagCOP, DataMap & ) override;
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "DetectProteinLigandInterface"; }

	void
	find_design_interface(
		core::pose::Pose const & pose,
		std::set< core::Size > const & interface_target_res,
		core::Real cut1,
		core::Real cut2,
		core::Real cut3,
		core::Real cut4,
		utility::vector1< bool > & repack_res,
		utility::vector1< bool > & design_res
	) const;

	void
	find_design_interface_arg_sweep(
		core::pose::Pose const & pose,
		std::set< core::Size > const & interface_target_res,
		core::Real cut1,
		core::Real cut2,
		core::Real cut3,
		core::Real cut4,
		core::Real arg_sweep_cutoff,
		utility::vector1< bool > & repack_res,
		utility::vector1< bool > & design_res
	) const;

	static void register_options();

	bool get_design() const { return design_;}
	void set_design(bool const design_in) {design_ = design_in;}

	bool get_no_design_cys() const { return no_design_cys_;}
	void set_no_design_cys(bool const no_design_cys_in) {no_design_cys_ = no_design_cys_in;}

	void
	set_design_target_res( std::set< core::Size > const & target_res )
	{ design_target_res_ = target_res; }

	static
	void
	add_observer_cache_segments_to_set(
		core::pose::Pose const & pose,
		std::set< core::Size > & set
	);

private:
	/// Control whether the detect design interface algorithm is run.
	bool detect_design_interface_;

	/// Use protein-DNA interface like arginine rotamer sweep to identify designable positions
	bool arg_sweep_interface_;
	/// Dicitate the max distance from arginine to ligand for position to be counted as designable
	core::Real arg_sweep_cutoff_;

	/// interface can be declared as all catalytic res
	bool catalytic_res_part_of_interface_;
	/// Depending on the design_ variable setting the cut1 and cut2 can be turned off and no design will take place.
	bool design_;
	/// Turn of design (repack_only_), or repacking and design (score_only_)
	// repack_only_ & score_only_ are "global" options - design_ is more temporary/local
	bool repack_only_, score_only_;
	/// cut1_ through cut4_ define the shells over which the interface detection takes place
	core::Real cut1_,cut2_,cut3_,cut4_;
	/// If resfilename_ is not empty, load it and obey the AUTO directives therein
	std::string resfilename_;
	std::set< core::Size > design_target_res_;
	bool add_observer_cache_segs_to_interface_;
	/// Should we prohibit designing to non disulfide cys?
	bool no_design_cys_;
	bool catres_only_;
	bool use_cstid_list_;
	std::string cstid_list_;
};

/// @brief Class to alter a packer task to speficially upweight the protein-ligand interaction energies
class ProteinLigandInterfaceUpweighter: public core::pack::task::operation::TaskOperation
{

public:
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
	typedef core::pose::Pose Pose;
	typedef core::pack::task::PackerTask PackerTask;

public:

	ProteinLigandInterfaceUpweighter();
	~ProteinLigandInterfaceUpweighter() override;

	TaskOperationOP clone() const override;

	/// @brief Initialize the class based on the command line options.
	void init_from_options();

	/// @brief Change a packer task in some way.  The input pose is the one to which the input
	/// task will be later applied.
	void apply( Pose const & pose, PackerTask & task) const override;

	void parse_tag( TagCOP, DataMap & ) override;
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "ProteinLigandInterfaceUpweighter"; }

	static void register_options();

	bool get_weight() const { return lig_packer_weight_;}
	void set_weight(bool const weight_in) {lig_packer_weight_ = weight_in;}

private:
	/// Reweight protein-ligand interaction by a factor of lig_packer_weight_.
	core::Real lig_packer_weight_;
	core::Real catres_packer_weight_;
};

class AddRigidBodyLigandConfs : public core::pack::task::operation::TaskOperation
{

public:
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
	typedef core::pose::Pose Pose;
	typedef core::pack::task::PackerTask PackerTask;

public:

	AddRigidBodyLigandConfs();
	~AddRigidBodyLigandConfs() override;

	TaskOperationOP clone() const override;

	void parse_tag( TagCOP tag , DataMap & ) override;
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "AddRigidBodyLigandConfs"; }

	/// @brief Change a packer task in some way.  The input pose is the one to which the input
	/// task will be later applied.
	void apply( Pose const &, PackerTask & ) const override;

};

//Add ligand motif rotamers
class AddLigandMotifRotamers : public core::pack::task::operation::TaskOperation
{

public:
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
	typedef core::pose::Pose Pose;
	typedef core::pack::task::PackerTask PackerTask;

private:
	protocols::motifs::LigandMotifSearchOP motif_search_;


public:

	AddLigandMotifRotamers();
	~AddLigandMotifRotamers() override;

	TaskOperationOP clone() const override;

	static void register_options();

	void parse_tag( TagCOP, DataMap & ) override;
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "AddLigandMotifRotamers"; }

	/// @brief Change a packer task in some way.  The input pose is the one to which the input
	/// task will be later applied.
	void apply( Pose const &, PackerTask & ) const override;

};

}//namespace protocols
}//namespace protocols

#endif
