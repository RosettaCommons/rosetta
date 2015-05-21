// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/design/AntibodyCDRSetInstructions.hh
/// @brief Create and hold AntibodyCDRSetInstructions for the graft design step
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_antibody_database_CDRSetOptions_hh
#define INCLUDED_protocols_antibody_database_CDRSetOptions_hh

#include <protocols/antibody/database/CDRSetOptions.fwd.hh>
#include <protocols/antibody/clusters/CDRClusterEnum.hh>
#include <protocols/antibody/AntibodyEnum.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace antibody {
	
//typedef std::map<CDRNameEnum, CDRSetInstructions> AntibodyCDRSetOptions;
	typedef utility::vector1<CDRSetOptionsOP> AntibodyCDRSetOptions;
	
/// @brief Class that holds instructions for a single CDR for loading from the antibody database.
/// Default is instructions to load all CDRs from the database
class CDRSetOptions : public utility::pointer::ReferenceCount {
	
public:

	CDRSetOptions(bool load=true);
	CDRSetOptions(CDRNameEnum cdr, bool load=true);
	
	CDRSetOptions(CDRSetOptions const & src);
	
	//CDRSetOptions(std::string instructions);
	
	~CDRSetOptions();
	
	
public:
	////// General settings //////////////
	
	void
	set_cdr(CDRNameEnum cdr);
	
	CDRNameEnum
	cdr() const {
		return cdr_;
	}
	
	void
	load(bool load);
	
	bool
	load() const {
		return load_;
	}
	
	
	
	/////////// Length Types ////////////
	
	/// @brief Set to only sample with clusters of the given type for this CDR.
	void
	length_type(core::Size const type, bool const setting);
	
	utility::vector1<bool>
	length_type() const {
		return length_types_;
	}
	
	
	////////////// CDR Lengths /////////////////
	
	/// @brief set the minimum cdr length to sample.  Nothing shorter then this will be used during graft.
	void
	min_length(core::Size length);
	
	core::Size
	min_length() const {
		return min_length_;
	}
	
	
	/// @brief set the maximum cdr length to sample.  Nothing longer then this will be used.  Useful for H3.
	void
	max_length(core::Size length);
	
	core::Size
	max_length() const {
		return max_length_;
	}
	
	
	////////////// PDBIds ///////////////////
	
	void
	exclude_pdbs(utility::vector1<std::string> exclude_pdbids);
	
	void
	exclude_pdbs_add(std::string pdbid);
	
	void
	exclude_pdbs_clear();
	
	utility::vector1<std::string>
	exclude_pdbs()  const;
	
	
	void
	include_only_pdbs(utility::vector1<std::string> include_pdbids);
	
	void
	include_only_pdbs_add(std::string pdbid);
	
	void
	include_only_pdbs_clear();
	
	utility::vector1<std::string>
	include_only_pdbs() const;
	
	////////////// Clusters /////////////////
	
	void
	exclude_clusters(utility::vector1< clusters::CDRClusterEnum > exclude);
	
	void
	exclude_clusters_add(clusters::CDRClusterEnum exclude);
	
	void
	exclude_clusters_clear();
	
	utility::vector1<clusters::CDRClusterEnum>
	exclude_clusters() const;
	
	
	void
	include_only_clusters(utility::vector1<clusters::CDRClusterEnum> include_only);
	
	void
	include_only_clusters_add(clusters::CDRClusterEnum include_only);
	
	void
	include_only_clusters_clear();
	
	utility::vector1<clusters::CDRClusterEnum>
	include_only_clusters() const;
	
	///////////// Species ////////////////////
	
	void
	exclude_species(utility::vector1<std::string> exclude);
	
	void
	exclude_species_add(std::string exclude);
	
	void
	exclude_species_clear();
	
	utility::vector1<std::string>
	exclude_species() const;
	
	void
	include_only_species(utility::vector1<std::string> include_only);
	
	void
	include_only_species_add(std::string include_only);
	
	void
	include_only_species_clear();
	
	utility::vector1<std::string>
	include_only_species() const;
	
	
	/////////// Germlines ////////////////////
	
	void
	exclude_germlines(utility::vector1<std::string> exclude);
	
	void
	exclude_germlines_add(std::string exclude);
	
	void
	exclude_germlines_clear();
	
	utility::vector1<std::string>
	exclude_germlines() const;
	
	void
	include_only_germlines(utility::vector1<std::string> include_only);
	
	void
	include_only_germlines_add(std::string include_only);
	
	void
	include_only_germlines_clear();
	
	utility::vector1<std::string>
	include_only_germlines() const;
	
	/////// Current Cluster ////////
	
	/// @brief
	/// If antibody CDRs are loaded in relation to the current PDB in whatever algorithm or app is using this,
	///  Should we only pull clusters of the same type as this PDB?
	/// @details
	/// May or may not be used.  In AntibodyGraftDesign, this is used.  
	///
	void
	include_only_current_cluster(bool only_current_cluster);
	
	bool
	include_only_current_cluster() const;
	
	/// @brief
	/// Do we only include the center members of the clusters in the CDRSet?
	void
	include_only_center_clusters(bool centers_only);
	
	bool
	include_only_center_clusters() const;
	
public:
	CDRSetOptionsOP
	clone() const;
	
	void
	set_defaults();
	
private:
	CDRNameEnum cdr_;
	bool load_;
	
	bool only_current_cluster_;
	bool only_center_clusters_;
	
	utility::vector1< bool>  length_types_; //Only sample within cluster types (Cluster types 1,2,3)
	utility::vector1< std::string > exclude_pdb_ids_;
	utility::vector1< std::string > include_only_pdb_ids_;
	utility::vector1< clusters::CDRClusterEnum > exclude_clusters_;//Leave out these clusters.  
	utility::vector1< clusters::CDRClusterEnum > include_only_clusters_; //Include only these clusters.
	core::Size min_length_;
	core::Size max_length_;
		
	//Additional tweaking parameters:
	utility::vector1< std::string > exclude_species_;
	utility::vector1< std::string > include_only_species_;
		
	utility::vector1< std::string > exclude_germlines_;
	utility::vector1< std::string > include_only_germlines_;
	
};


			
}
}

#endif

