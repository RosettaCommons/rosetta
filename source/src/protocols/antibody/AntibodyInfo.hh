// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
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
#include <protocols/antibody/CDRClusterEnum.hh>
#include <protocols/antibody/CDRClusterEnumManager.hh>

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
using namespace core;
using namespace utility;
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
class AntibodyInfo : public pointer::ReferenceCount {

public:
	AntibodyInfo( pose::Pose const & pose,
	              AntibodyNumberingSchemeEnum const & numbering_scheme = Aroop,
	              bool const & cdr_pdb_numbered = true);

	//Default destructor
	virtual ~AntibodyInfo();

public:

	/// @brief: get the current numbering scheme being used
	std::string
	get_Current_AntibodyNumberingScheme();

	/// @brief get the length of the cdr
	Size
	get_CDR_length(CDRNameEnum const & cdr_name) const;
	
	Size
	get_CDR_length(CDRNameEnum const & cdr_name, pose::Pose const & pose);
	
	std::string
	get_CDR_Name(CDRNameEnum const & cdr_name) const;

	char
	get_CDR_chain(CDRNameEnum const & cdr_name) const {
		return Chain_IDs_for_CDRs_[cdr_name];
	}

	/// @brief return this antibody is camelid or not
	bool
	is_camelid()  const {
		return is_camelid_;
	}

	/// @brief return whether this pose has antigen or not
	bool
	get_pose_has_antigen() const {
		return InputPose_has_antigen_;
	}

	/// @brief return num of cdr loops, 3 (nanobody) or 6 (regular antibody)
	CDRNameEnum
	get_total_num_CDRs() const {
		return total_cdr_loops_;
	}

	/// @brief Return PDB residue number of CDR start residue
	Size
	get_CDR_start_PDB_num(CDRNameEnum const & cdr_name) const {
		return cdr_numbering_[cdr_name][start];
	}

	/// @brief Return PDB residue number of CDR end residue
	Size
	get_CDR_end_PDB_num(CDRNameEnum const & cdr_name) const {
		return cdr_numbering_[cdr_name][stop]; // PDB numbering
	}

	/// @brief Return pose number of CDR start residue
	Size
	get_CDR_start(CDRNameEnum const & cdr_name, pose::Pose const & pose) const;

	/// @brief Return pose number of CDR end residue
	Size
	get_CDR_end(CDRNameEnum const & cdr_name, pose::Pose const & pose) const; // pose numbering

	/// @brief return the framework numbering information
	vector1< vector1<FrameWork> >
	get_AntibodyFrameworkInfo() const ;

	/// @brief get H3 cterminal kink/extended conformation (predicted by constructor)
	H3BaseTypeEnum
	get_Predicted_H3BaseType() const ;

	/// @brief get residues used to calculate VL/VH packing angle
	vector1< Size >
	get_PackingAngleResidues() const {
		return packing_angle_residues_;
	}
	
	/// @brief gets all non-LH chains.  Empty vector if no antigen present.
	vector1< char >
	get_antigen_chains() const {
		return Chain_IDs_for_antigen_;
	}

	/// @brief return pose residue number of the first residue of the H3 kink
	 Size
	kink_begin() const {
	return get_CDR_loop(h3).stop() - 2;
	}

	/// @brief return pose residue number of the last residue of the H3 kink
	Size
	kink_end() const {
		return get_CDR_loop(h3).stop() + 1;
	}

	/// @brief return pose residue number of the last residue of the H3 kink
	Size
	kink_trp() const {
		return get_CDR_loop(h3).stop() + 1;
	}

	  /// @brief return pose residue number of the kink 'anion' (typically Asp/Glu) residue in the kink bulge HBond
	Size
	kink_anion_residue() const {
		return get_CDR_loop(h3).stop() - 1;
	}

	/// @brief return pose residue number of the kink 'cation' (typically Arg/Lys) residue in the kink bulge HBond
	Size
	kink_cation_residue() const {
		return get_CDR_loop(h3).start() - 1;
	}
    
	/// @brief return side chain anion atoms (typically Asp/Glu oxygens) in the kink bulge HBond
	std::vector<Vector>
	kink_anion_atoms(const core::pose::Pose & pose) const;
    
	/// @brief return side chain cation atoms (typically Lys/His nitrogens) in the kink bulge HBond
	std::vector<Vector>
	kink_cation_atoms(const core::pose::Pose & pose) const;
    

public:
	//////////////////////////////////////////////////
	//CDR Clusters
	//
	//

	/// @brief setup the clustering information for each CDR
	void
	setup_CDR_clusters(pose::Pose const & pose);

	///@brief Check if the Cluster information is already set.  Useful for protocols.
	bool
	clusters_setup() const;
	
	/// @brief get the cdr's cluster identity and distance to cluster using it's structure
	/// @details See North, B., A. Lehmann, et al. (2011). JMB 406(2): 228-256.
	///  Must use setup_CDR_clusters first.
	std::pair < CDRClusterEnum, Real >
	get_CDR_cluster(CDRNameEnum const cdr_name) const;
	
	std::string
	get_cluster_name(CDRClusterEnum const cluster) const;
	
	CDRClusterEnum
	get_cluster_enum(std::string const cluster) const;
	
	CDRNameEnum
	get_cdr_enum_for_cluster(CDRClusterEnum const cluster) const;

	core::Size
	get_cluster_length(CDRClusterEnum const cluster) const;
	
public:
	///////////////////////////////////////////////////
	//Sequence
	//
	//
	/// @brief return the sequence of a particular CDR loop
	vector1<char>
	get_CDR_sequence_with_stem( CDRNameEnum const & cdr_name,
	                            Size left_stem = 0,
	                            Size right_stem = 0) const;

	/// @brief return the antibody sequence of LH or just H for camelid
	vector1<char> const &
	get_antibody_sequence() const {
		return ab_sequence_;
	}


public:
	///////////////////////////////////////////////////
	//Loops
	//
	//

	
	/// @brief return the loop of a certain loop type
	loops::Loop
	get_CDR_loop( CDRNameEnum const & cdr_name ) const {
		return (*get_CDR_in_loopsop(cdr_name))[1];
	}
	
	///@brief return the loop of a certain loop type on the fly
	loops::Loop
	get_CDR_loop( CDRNameEnum const & cdr_name, pose::Pose const & pose) const;

	/// @brief return the loop of a certain loop type
	loops::LoopsOP
	get_CDR_in_loopsop( CDRNameEnum const & cdr_name ) const {
		return vector1_loopsop_having_cdr_[cdr_name];
	}
	
	/// @brief return a LoopsOP object, initialized upon class construction.
	loops::LoopsOP
	get_AllCDRs_in_loopsop() const {
		return loopsop_having_allcdrs_;
	}
	
	///@brief On-the-fly CDR LoopsOP
	loops::LoopsOP
	get_CDR_loops(pose::Pose const & pose) const;
	
public:
	///////////////////////////////////////////////////
	//FoldTrees
	//
	//
	// FoldTrees //TODO: find a way to remove setup_simple_fold_tree
	kinematics::FoldTreeCOP
	setup_simple_fold_tree(Size const & jumppoint1,
	                       Size const & cutpoint,
	                       Size const & jumppoint2,
	                       pose::Pose const & pose ) const;

	kinematics::FoldTreeCOP
	get_FoldTree_AllCDRs_LHDock( pose::Pose const & pose) const;

	kinematics::FoldTreeCOP
	get_FoldTree_AllCDRs(pose::Pose const & pose) const;

	/// @brief SnugDock foldtrees
	kinematics::FoldTree
	get_FoldTree_LH_A( pose::Pose const & pose ) const;

	kinematics::FoldTree
	get_FoldTree_L_HA( pose::Pose const & pose ) const;

	kinematics::FoldTree
	get_FoldTree_LA_H( pose::Pose const & pose ) const;


public:
	///////////////////////////////////////////////////
	//MoveMaps
	//
	//
	/// TODO: this should be a standard utility for loops?
	/// @brief get a movemap for loops.
	kinematics::MoveMap
	get_MoveMap_for_Loops(pose::Pose const & pose,
	                      loops::Loops const & the_loops,
	                      bool const & bb_only = false,
	                      bool const & include_nb_sc = false,
	                      Real const & nb_dist = 10.0) const;

	kinematics::MoveMap
	get_MoveMap_for_AllCDRsSideChains_and_H3backbone(pose::Pose const & pose,
	        bool const & include_nb_sc = false,
	        Real const & nb_dist = 10.0) const;

	/// @brief Add CDR flexibility to a movemap. Uses pose for numbering.
	void
	add_CDR_to_MoveMap(pose::Pose const & pose,
	                   kinematics::MoveMapOP movemap,
	                   CDRNameEnum const & cdr_name,
	                   bool const & bb_only = false,
	                   bool const & include_nb_sc = false,
	                   Real const & nb_dist = 10.0) const;

	/// @brief get a movemap for loops and set the first jump be true
	kinematics::MoveMap
	get_MoveMap_for_LoopsandDock(pose::Pose const & pose,
	                             loops::Loops const & the_loops,
	                             bool const & bb_only = false,
	                             bool const & include_nb_sc = false,
	                             Real const & nb_dist = 10.0) const;


public:
	///////////////////////////////////////////////////
	//TaskFactories
	//
	//
	/// @brief TaskFactory
	pack::task::TaskFactoryOP
	get_TaskFactory_AllCDRs(pose::Pose &  pose) const;

	pack::task::TaskFactoryOP
	get_TaskFactory_OneCDR(pose::Pose & pose, CDRNameEnum const & cdr_name) const;
    
	
public:

	/// @brief use the H3 cterm coordinates in the pose to calculate the cterminal type
	//std::string calculate_H3_base_by_coordinates(pose::Pose const & pose) const;
	AntibodyEnumManagerOP
	get_antibody_enum_manager() const;

	CDRClusterEnumManagerOP
	get_cdr_cluster_enum_manager() const;

	void show( std::ostream & out=std::cout );
	friend std::ostream & operator<<(std::ostream& out, const AntibodyInfo & ab_info ) ;


private:
	/////////////////////////////////////////////////////////////////////////////////////////
	///Setters for private AntibodyInfo variables									  ///
	/////////////////////////////////////////////////////////////////////////////////////////
	///
	void set_default();
    
	/// @brief check the input pose is nanobody, antibody or wrong.  Sets sequence.  Sets Antigen chains.
	void identify_antibody(pose::Pose const & pose);

	void init(pose::Pose const & pose);

	/// @brief Setup the Internal AntibodyNumbering variables
	//  Examples:
	//           cdr_numbering_[h1][start]
	//           packing_numbering_[VL_sheet_1][stop]
	void setup_numbering_info_for_scheme(AntibodyNumberingSchemeEnum const & numbering_scheme);
	
	/// @brief setup the CDR loops objects based on the input numbering scheme
	void setup_CDRsInfo( pose::Pose const & pose );

	/// @brief setup the framework information based on the input numbering scheme
	void setup_FrameWorkInfo(pose::Pose const & pose );

	/// @brief setup the residues used to calculate VL/VH packing angle
	void setup_VL_VH_packing_angle( pose::Pose const & pose );

	/// @brief calculate the clustering for each CDR
	void set_cdr_cluster(pose::Pose const & pose, CDRNameEnum const & cdr_name);

	/// @brief predict H3 cterminus base as Kinked or Extended
	void predict_H3_base_type( pose::Pose const & pose ) ;
	void detect_and_set_camelid_CDR_H3_stem_type( pose::Pose const & pose );
	void detect_and_set_regular_CDR_H3_stem_type( pose::Pose const & pose );
	void detect_and_set_regular_CDR_H3_stem_type_new_rule( pose::Pose const & pose );

	/// @brief identify CDRs on L or H sequence
	void identify_CDR_from_a_sequence(std::string const & querychain);
	///																					  ///
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
	bool InputPose_has_antigen_;
	bool cdr_pdb_numbered_;

	/// the CDR and Framework information
	vector1<loops::LoopsOP> vector1_loopsop_having_cdr_; // each Loops contains one CDR
	loops::LoopsOP loopsop_having_allcdrs_;  // one Loops object containing the set of Loop objects for all CDRs
	vector1< vector1<FrameWork> > framework_info_ ;
	vector1<char> ab_sequence_;
	std::map<Size, char> sequence_map_; //Residue number sequence map
	vector1< Size > packing_angle_residues_;

	// Clusters - not set on construction.
	vector1< CDRClusterEnum > cdr_clusters_;
	vector1< core::Real > cdr_cluster_distances_;

	/// Antibody properties
	AntibodyNumberingSchemeEnum numbering_scheme_;
	H3BaseTypeEnum predicted_H3_base_type_;
	CDRNameEnum total_cdr_loops_;
	vector1<char> Chain_IDs_for_CDRs_;
	vector1<char> Chain_IDs_for_antigen_;
	
	Size L_chain_;
	Size H_chain_;

	///Internal Per-Scheme Numbering
	AntibodyNumbering cdr_numbering_;
	AntibodyNumbering packing_angle_numbering_;
	// AntibodyNumbering framework_numbering_;

	//Enum Managers
	AntibodyEnumManagerOP antibody_manager_;
	CDRClusterEnumManagerOP cdr_cluster_manager_;

};


} //namespace antibody
} //namespace protocols


#endif //INCLUDED_protocols_loops_AntibodyInfo_HH
