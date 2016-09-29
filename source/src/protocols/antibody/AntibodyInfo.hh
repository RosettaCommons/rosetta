// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file protocols/antibody/AntibodyInfo.hh
/// @brief Class for getting antibody-specific objects and information
/// @author Jianqing Xu (xubest@gmail.com)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_AntibodyInfo_hh
#define INCLUDED_protocols_antibody_AntibodyInfo_hh

#include <protocols/antibody/AntibodyInfo.fwd.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/AntibodyEnumManager.hh>
#include <protocols/antibody/AntibodyNumberingParser.hh>
#include <protocols/antibody/clusters/CDRClusterEnum.hh>
#include <protocols/antibody/clusters/CDRClusterEnumManager.hh>
#include <protocols/antibody/clusters/CDRClusterSet.hh>
#include <protocols/antibody/clusters/CDRCluster.hh>

// Rosetta Headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <utility/vector1.hh>
#include <protocols/docking/types.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/kinematics/FoldTree.hh>


///////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace antibody {


struct FrameWork {
	Size    start;
	Size    stop;
	std::string chain_name;
};


/// @brief This class is used to get all relevant information you would need when dealing with an antibody.
///  @details It mainly holds numbering information, but passes out a variety of Rosetta specific objects like movemaps, Loops, taskfactories, etc.
///  as well as other information such as CDR cluster type, camelid, H3 types, etc.
class AntibodyInfo : public utility::pointer::ReferenceCount {

public:

	/// @brief Constructor that Loads default numbering scheme and cdr definition from options (Chothia/Aroop)
	AntibodyInfo( core::pose::Pose const & pose, bool const cdr_pdb_numbered = true);

	/// @brief Constructor that uses the numbering scheme and cdr_definition given.
	AntibodyInfo( core::pose::Pose const & pose,
		AntibodyNumberingSchemeEnum const & numbering_scheme,
		CDRDefinitionEnum const cdr_definition,
		bool const cdr_pdb_numbered = true);

	/// @brief Constructor using just the cdr definition
	AntibodyInfo( core::pose::Pose const & pose,
		CDRDefinitionEnum const cdr_definition,
		bool const cdr_pdb_numbered = true);

	/// @brief Copy Constructor
	AntibodyInfo(AntibodyInfo const & src);

	/// @brief Clone the pointer.
	AntibodyInfoOP
	clone() const;

	//Default destructor
	~AntibodyInfo() override;

public:

	/// @brief: get the current numbering scheme being used
	std::string
	get_current_AntibodyNumberingScheme() const ;

	std::string
	get_current_CDRDefinition()  const;

	/// @brief get the length of the cdr upon construction of AntibodyInfo
	///
	/// @details Does not work for proto CDR4 (DE loop)
	core::Size
	get_CDR_length(CDRNameEnum const cdr_name) const;

	/// @breif On-the-fly CDR length.
	core::Size
	get_CDR_length(CDRNameEnum const cdr_name, core::pose::Pose const & pose) const;

	core::Size
	get_CDR_length(CDRNameEnum const cdr_name, core::pose::Pose const & pose, CDRDefinitionEnum const transform) const;

	std::string
	get_CDR_name(CDRNameEnum const cdr_name) const;

	CDRNameEnum
	get_CDR_name_enum(std::string const & cdr_name) const;

	char
	get_CDR_chain(CDRNameEnum const cdr_name) const {
		return chains_for_cdrs_[cdr_name];
	}

	/// @brief return this antibody is camelid or not
	bool
	is_camelid()  const {
		return is_camelid_;
	}

	/// @brief Return if antibody is lambda, kappa or unknown.  Type set via cmd-line option for now
	std::string
	get_light_chain_type() const {
		return enum_manager_->light_chain_type_enum_to_string(light_chain_type_);
	}

	/// @brief Return if antibody is lambda, kappa or unknown.  Type set via cmd-line option for now
	LightChainTypeEnum
	get_light_chain_type_enum() const {
		return light_chain_type_;
	}

	/// @brief return whether this pose has antigen or not
	bool
	antigen_present() const {
		return has_antigen_;
	}

	/// @brief return whether a cdr contacts antigen.  If no antigen is present, returns false.
	/// @details Considered 'in_contact' if > 5 atoms of antigen are within 5 A of any atom of the CDR.
	//bool
	//contacts_antigen(pose::Pose const & pose, CDRNameEnum const cdr_name) const;

	/// @brief return num of cdr loops, 3 (nanobody) or 6 (regular antibody)
	CDRNameEnum
	get_total_num_CDRs(bool include_proto_cdr4 = false) const;

	/// @brief Does the antibody have this CDR?
	/// See also:
	///                   is_camelid())
	///                   get_total_num_CDRs())
	///
	bool
	has_CDR(CDRNameEnum const cdr_name) const;

	/// @brief Return PDB residue number of CDR start residue
	core::Size
	get_CDR_start_PDB_num(CDRNameEnum const cdr_name) const {
		return numbering_info_.cdr_numbering[cdr_name][cdr_start]->resnum();
	}

	/// @brief Return PDB residue number of CDR end residue
	core::Size
	get_CDR_end_PDB_num(CDRNameEnum const cdr_name) const {
		return numbering_info_.cdr_numbering[cdr_name][cdr_end]->resnum(); // PDB numbering
	}


	/// @brief Return pose number of CDR start residue
	core::Size
	get_CDR_start(CDRNameEnum const cdr_name, core::pose::Pose const & pose) const;

	// @brief Return pose number of CDR start residue in the definition of the numbering scheme of transform.
	core::Size
	get_CDR_start(CDRNameEnum const cdr_name, core::pose::Pose const & pose, CDRDefinitionEnum const & transform) const;

	/// @brief Return pose number of CDR end residue in the definition of the numbering scheme of transform.
	core::Size
	get_CDR_end(CDRNameEnum const cdr_name, core::pose::Pose const & pose) const; // pose numbering

	core::Size
	get_CDR_end(CDRNameEnum const cdr_name, core::pose::Pose const & pose, CDRDefinitionEnum const & transform) const;

	/// @brief Get the region of the resnum - aka - antigen_region, cdr_region, or framework_region
	///  By default, count the DE region (CDR4) as framework.
	AntibodyRegionEnum
	get_region_of_residue(core::pose::Pose const & pose, core::Size resnum, bool de_as_framework = true) const;

	/// @brief return the framework numbering information
	utility::vector1< utility::vector1<FrameWork> >
	get_AntibodyFrameworkInfo() const ;

	/// @brief get H3 cterminal kink/extended conformation (predicted by constructor)
	H3BaseTypeEnum
	get_H3_kink_type() const ;

	std::string
	get_H3_kink_type_name() const ;

	/// @brief get residues used to calculate VL/VH packing angle
	utility::vector1< Size >
	get_PackingAngleResidues() const {
		return packing_angle_residues_;
	}

	/////////////////////////////Antigen and Antibody Chains ////////////////////////////////////////

	/// @brief gets all non-LH chains.  Empty vector if no antigen present.
	utility::vector1< char >
	get_antigen_chains() const {
		return chains_for_antigen_;
	}
	utility::vector1<core::Size>
	get_antigen_chain_ids(const core::pose::Pose & pose ) const;

	/// @brief Return the antigen chains as a string
	std::string
	get_antigen_chain_string() const;


	utility::vector1< char >
	get_antibody_chains() const;

	utility::vector1<core::Size>
	get_antibody_chain_ids(const core::pose::Pose & pose ) const;

	/// @brief Returns H or LH depeding on camelid.
	std::string
	get_antibody_chain_string() const;


	/// @brief return pose residue number of the first residue of the H3 kink
	core::Size
	kink_begin(core::pose::Pose const & pose) const {
		return get_CDR_loop(h3, pose, Aroop).stop() - 2;
	}

	/// @brief return pose residue number of the last residue of the H3 kink
	core::Size
	kink_end(core::pose::Pose const & pose) const {
		return get_CDR_loop(h3, pose, Aroop).stop() + 1;
	}

	/// @brief return pose residue number of the last residue of the H3 kink
	core::Size
	kink_trp(core::pose::Pose const & pose) const {
		return get_CDR_loop(h3, pose, Aroop).stop() + 1;
	}

	/// @brief return pose residue number of the kink 'anion' (typically Asp/Glu) residue in the kink bulge HBond
	core::Size
	kink_anion_residue(core::pose::Pose const & pose) const {
		return get_CDR_loop(h3, pose, Aroop).stop() - 1;
	}

	/// @brief return pose residue number of the kink 'cation' (typically Arg/Lys) residue in the kink bulge HBond
	Size
	kink_cation_residue(core::pose::Pose const & pose) const {
		return get_CDR_loop(h3, pose, Aroop).start() - 1;
	}

	/// @brief return side chain anion atoms (typically Asp/Glu oxygens) in the kink bulge HBond
	std::vector<core::Vector>
	kink_anion_atoms(const core::pose::Pose & pose) const;

	/// @brief return side chain cation atoms (typically Lys/His nitrogens) in the kink bulge HBond
	std::vector<core::Vector>
	kink_cation_atoms(const core::pose::Pose & pose) const;

public:
	///////////////////////////////////////////////////
	//Landmark Transform
	//
	//


	/// @brief Used to get a residue number of a particular place in the framework or conserved residue in the CDR.
	/// Use this instead of PDBInfo!!
	///
	/// @details If the current numbering scheme is not 'scheme', will return the equivalent position of the current numbering scheme using the transform scheme file in the database.
	/// Should not be used for residues within CDRs since numbering schemes vary greatly in their within cdr alignments and numbering.  Use get_CDR_start/end/loop functions with relative positions for this purpose.
	/// Returns 0 if residue is not found in pose
	core::Size
	get_landmark_resnum(core::pose::Pose const & pose,
		AntibodyNumberingSchemeEnum const scheme,
		char const chain,
		core::Size const pdb_resnum,
		char const insertion_code=' ',
		bool fail_on_missing_resnum = true) const;

	///@brief Get the full antibody numbering struct for use in anything.
	AntibodyNumbering const &
	get_antibody_numbering_info() const {
		return numbering_info_;
	}

public:
	//////////////////////////////////////////////////
	//CDR Clusters
	//

	/// @brief setup the clustering information for each CDR to totalCDRLoops
	///  Optionally, attempt to set the CDR cluster from the Pose data cache if any data is there.
	void
	setup_CDR_clusters(core::pose::Pose const & pose, bool attempt_set_from_pose = false);

	/// @brief setup the clustering information for each CDR according to boolean vector.
	///  Optionally, attempt to set the CDR cluster from the Pose data cache if any data is there.
	void
	setup_CDR_clusters(core::pose::Pose const & pose, utility::vector1<bool> const & cdrs, bool attempt_set_from_pose=false);

	/// @brief setup the clustering information for a particular CDR.
	///  Optionally, attempt to set the CDR cluster from the Pose data cache if any data is there.
	void
	setup_CDR_cluster(core::pose::Pose const & pose, CDRNameEnum cdr, bool attempt_set_from_pose=false);

	/// @brief Manually set the CDR cluster for a particular CDR
	void
	set_CDR_cluster(CDRNameEnum cdr, clusters::CDRClusterCOP cluster);

	/// @brief get the cdr's cluster identity and distance to cluster using it's structure as a CDRCluster object.
	/// @details See North, B., A. Lehmann, et al. (2011). JMB 406(2): 228-256.
	///  Must use setup_CDR_clusters first.
	/// @note If this AbInfo does not have all 6 chains, CDRClusterOP may be NULL if the CDR is not present to begin with.
	///
	clusters::CDRClusterCOP
	get_CDR_cluster(CDRNameEnum const cdr_name) const;

	clusters::CDRClusterSetOP
	get_CDR_cluster_set() const;

	/// @brief Check to make sure AbInfo has a cluster object for this CDR.  In that all 6 CDRs not nessessarily present for each antibody.
	bool
	has_cluster_for_cdr(CDRNameEnum const cdr_name) const;

	std::string
	get_cluster_name(clusters::CDRClusterEnum const cluster) const;

	clusters::CDRClusterEnum
	get_cluster_enum(std::string const & cluster) const;

	CDRNameEnum
	get_cdr_enum_for_cluster(clusters::CDRClusterEnum const cluster) const;

	core::Size
	get_cluster_length(clusters::CDRClusterEnum const cluster) const;

public:
	///////////////////////////////////////////////////
	//Sequence
	//
	//
	/// @brief return the sequence of a particular CDR loop
	///  Deprecated in favor of on the fly sequence.
	std::string
	get_CDR_sequence_with_stem( CDRNameEnum const cdr_name,
		Size left_stem = 0,
		Size right_stem = 0) const;

	///@brief On-The-Fly CDR sequence
	std::string
	get_CDR_sequence_with_stem( CDRNameEnum const cdr_name,
		core::pose::Pose const & pose,
		Size left_stem = 0,
		Size right_stem = 0) const;
	
	///@brief On-The-Fly CDR sequence
	std::string
	get_CDR_sequence_with_stem( CDRNameEnum const cdr_name,
		core::pose::Pose const & pose,
		Size left_stem,
		Size right_stem,
		CDRDefinitionEnum const & transform) const;

	/// @brief return the antibody sequence of LH or just H for camelid
	std::string
	get_antibody_sequence() const;


public:
	///////////////////////////////////////////////////
	//Loops
	//


	/// @brief return the loop of a certain loop type
	loops::Loop
	get_CDR_loop( CDRNameEnum const cdr_name ) const;

	/// @brief return the loop of a certain loop type on the fly
	loops::Loop
	get_CDR_loop( CDRNameEnum const cdr_name, core::pose::Pose const & pose, core::Size overhang = 0) const;

	/// @brief return the loop of a certain loop type in definitions of the numbering scheme transform.
	loops::Loop
	get_CDR_loop(
		CDRNameEnum const cdr_name,
		core::pose::Pose const & pose,
		CDRDefinitionEnum const transform,
		core::Size overhang = 0) const;



	/// @brief return the loop of a certain loop type
	// Please don't use this.  It is being deprecated.
	loops::LoopsOP
	get_CDR_in_loopsop( CDRNameEnum const cdr_name ) const;

	/// @brief return a LoopsOP object, initialized upon class construction.
	// Please don't use this.  It is being deprecated.
	loops::LoopsOP
	get_AllCDRs_in_loopsop() const;

	/// @brief On-the-fly CDR LoopsOP
	loops::LoopsOP
	get_CDR_loops(core::pose::Pose const & pose, core::Size overhang = 0) const;

	loops::LoopsOP
	get_CDR_loops( core::pose::Pose const & pose, const utility::vector1<bool> & cdrs, core::Size overhang = 0) const;

public:
	///////////////////////////////////////////////////
	//FoldTrees
	//
	//
	// FoldTrees //TODO: find a way to remove setup_simple_fold_tree
	core::kinematics::FoldTreeCOP
	setup_simple_fold_tree(core::Size const & jumppoint1,
		Size const & cutpoint,
		Size const & jumppoint2,
		core::pose::Pose const & pose ) const;

	core::kinematics::FoldTreeCOP
	get_FoldTree_AllCDRs_LHDock( core::pose::Pose const & pose) const;

	core::kinematics::FoldTreeCOP
	get_FoldTree_AllCDRs(core::pose::Pose const & pose) const;

	/// @brief SnugDock foldtrees
	core::kinematics::FoldTree
	get_FoldTree_LH_A( core::pose::Pose const & pose ) const;

	core::kinematics::FoldTree
	get_FoldTree_L_HA( core::pose::Pose const & pose ) const;

	core::kinematics::FoldTree
	get_FoldTree_LA_H( core::pose::Pose const & pose ) const;


public:
	///////////////////////////////////////////////////
	//MoveMaps
	//
	//
	/// TODO: this should be a standard utility for loops?
	/// @brief get a movemap for loops.
	core::kinematics::MoveMap
	get_MoveMap_for_Loops(core::pose::Pose const & pose,
		loops::Loops const & the_loops,
		bool const & bb_only = false,
		bool const & include_nb_sc = false,
		core::Real const & nb_dist = 10.0) const;

	core::kinematics::MoveMap
	get_MoveMap_for_AllCDRsSideChains_and_H3backbone(core::pose::Pose const & pose,
		bool const & include_nb_sc = false,
		core::Real const & nb_dist = 10.0) const;

	/// @brief Add CDR flexibility to a movemap. Uses pose for numbering.
	void
	add_CDR_to_MoveMap(core::pose::Pose const & pose,
		core::kinematics::MoveMapOP movemap,
		CDRNameEnum const cdr_name,
		bool const & bb_only = false,
		bool const & include_nb_sc = false,
		core::Real const & nb_dist = 10.0) const;

	/// @brief get a movemap for loops and set the first jump be true
	core::kinematics::MoveMap
	get_MoveMap_for_LoopsandDock(core::pose::Pose const & pose,
		loops::Loops const & the_loops,
		bool const & bb_only = false,
		bool const & include_nb_sc = false,
		core::Real const & nb_dist = 10.0) const;


public:
	///////////////////////////////////////////////////
	//TaskFactories
	//
	//
	/// @brief TaskFactory
	core::pack::task::TaskFactoryOP
	get_TaskFactory_AllCDRs(core::pose::Pose &  pose) const;

	core::pack::task::TaskFactoryOP
	get_TaskFactory_OneCDR(core::pose::Pose & pose, CDRNameEnum const cdr_name) const;


public:

	/// @brief use the H3 cterm coordinates in the pose to calculate the cterminal type
	//std::string calculate_H3_base_by_coordinates(pose::Pose const & pose) const;
	AntibodyEnumManagerCOP
	get_antibody_enum_manager() const;

	clusters::CDRClusterEnumManagerCOP
	get_cdr_cluster_enum_manager() const;

	void show( std::ostream & out=std::cout ) const;
	friend std::ostream & operator<<(std::ostream& out, const AntibodyInfo & ab_info ) ;


private:
	/////////////////////////////////////////////////////////////////////////////////////////
	///Setters for private AntibodyInfo variables           ///
	/////////////////////////////////////////////////////////////////////////////////////////
	///

	///////////////////////////////////////////////////
	// Identification
	//
	//

	void set_default();

	/// @brief check the input pose is nanobody, antibody or wrong.  Sets sequence.  Sets Antigen chains.
	void identify_antibody(core::pose::Pose const & pose);

	/// @brief place-holder for identification of light chain type: lambda/kappa/unknown.
	void identify_light_chain(core::pose::Pose const & pose);

	void init(core::pose::Pose const & pose);

	/// @brief setup the CDR loops objects based on the input numbering scheme
	void setup_CDRsInfo( core::pose::Pose const & pose );

	/// @brief setup the framework information based on the input numbering scheme
	void setup_FrameWorkInfo(core::pose::Pose const & pose );

	/// @brief setup the residues used to calculate VL/VH packing angle
	void setup_VL_VH_packing_angle( core::pose::Pose const & pose );


	/// @brief predict H3 cterminus base as Kinked or Extended
	void predict_H3_base_type( core::pose::Pose const & pose ) ;
	void detect_and_set_camelid_CDR_H3_stem_type( core::pose::Pose const & pose );
	void detect_and_set_regular_CDR_H3_stem_type( core::pose::Pose const & pose );
	void detect_and_set_regular_CDR_H3_stem_type_new_rule( core::pose::Pose const & pose );

	/// @brief identify CDRs on L or H sequence
	void identify_CDR_from_a_sequence(std::string const & querychain);

	/////////////////////////////////////////////////////////////////////////////////////////
	//  Numbering
	//
	//

	/// @brief Setup the Internal AntibodyNumbering variables
	//  Examples:
	//           cdr_numbering_[h1][start]
	//           packing_numbering_[VL_sheet_1][stop]
	void setup_numbering_info_for_scheme(AntibodyNumberingSchemeEnum const & numbering_scheme, CDRDefinitionEnum const cdr_definition);

	/// @brief Gets transform for numbering scheme.  Allows const qualification for other functions.
	utility::vector1< utility::vector1< PDBLandmarkOP > >
	get_cdr_definition_transform(CDRDefinitionEnum const cdr_definition) const;

	/// @brief Is the transform from our numbering scheme to cdr_definition defined?
	bool
	cdr_definition_transform_present(CDRDefinitionEnum const cdr_definition) const;

	utility::vector1< PDBLandmarkOP >
	get_numbering_scheme_landmarks(AntibodyNumberingSchemeEnum const numbering_scheme) const;

	bool
	numbering_scheme_transform_present(AntibodyNumberingSchemeEnum const numbering_scheme) const;


	/// @brief Make sure there are no large chainbreaks - aka missing residues - in the CDR loops or really bad peptide bonds
	/// Controlled via cmd-line flags
	void check_cdr_quality( core::pose::Pose const & pose) const ;

	///                       ///
	/////////////////////////////////////////////////////////////////////////////////////////


	/// @brief copy
	//void init_for_equal_operator_and_copy_constructor( AntibodyInfo & lhs, AntibodyInfo const & rhs);

	/// @brief all the static functions
	static core::scoring::ScoreFunctionCOP get_Pack_ScoreFxn(void);
	static core::scoring::ScoreFunctionCOP get_Dock_ScoreFxn(void);
	static core::scoring::ScoreFunctionCOP get_LoopCentral_ScoreFxn(void);
	static core::scoring::ScoreFunctionCOP get_LoopHighRes_ScoreFxn(void);

	/// private data
private:

	/// the information of the antibody pose
	bool is_camelid_;
	bool has_antigen_;
	bool cdr_pdb_numbered_;

	/// the CDR and Framework information
	utility::vector1<loops::LoopsOP> vector1_loopsop_having_cdr_; // each Loops contains one CDR
	loops::LoopsOP loopsop_having_allcdrs_;  // one Loops object containing the set of Loop objects for all CDRs
	utility::vector1< utility::vector1<FrameWork> > framework_info_ ;
	utility::vector1<char> ab_sequence_;
	std::map<core::Size, char> sequence_map_; //Residue number sequence map
	utility::vector1< core::Size > packing_angle_residues_;
	utility::vector1<utility::vector1< core::Size > > packing_angle_numbering_; //In Chothia_Scheme numbering.
	/// Antibody properties
	H3BaseTypeEnum predicted_H3_base_type_;
	CDRNameEnum total_cdr_loops_;
	utility::vector1<char> chains_for_cdrs_;
	utility::vector1<char> chains_for_antigen_;

	core::Size L_chain_;
	core::Size H_chain_;

	LightChainTypeEnum light_chain_type_;

	clusters::CDRClusterSetOP cdr_cluster_set_;

	///Internal Per-Scheme Numbering
	AntibodyNumbering numbering_info_;
	utility::vector1< utility::vector1< core::Size> > current_transform_;

	// AntibodyNumbering framework_numbering_;

	//Enum Managers
	AntibodyEnumManagerOP enum_manager_;
	clusters::CDRClusterEnumManagerOP cdr_cluster_manager_;

	AntibodyNumberingParserOP numbering_parser_;

};


} //namespace antibody
} //namespace protocols


#endif //INCLUDED_protocols_loops_AntibodyInfo_HH
