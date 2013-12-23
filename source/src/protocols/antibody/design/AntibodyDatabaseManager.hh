// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/design/AntibodyDatabaseManager.hh
/// @brief Handles all loading of CDR, Framework, and cluster/dmap info from external database file.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_design_AntibodyDatabaseManager_hh
#define INCLUDED_protocols_antibody_design_AntibodyDatabaseManager_hh


#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/design/AntibodyDatabaseManager.fwd.hh>
#include <protocols/antibody/design/AntibodyDesignModeler.hh>
#include <protocols/antibody/design/AntibodyGraftDesigner.hh>
#include <protocols/antibody/design/AntibodyCDRDesigner.hh>
#include <protocols/antibody/AntibodyInfo.hh>

//Core Headers
#include <core/pose/Pose.hh>

//C++ Headers
#include <cppdb/frontend.h>
#include <string>
#include <map>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

namespace protocols {
namespace antibody {
namespace design{
	
	using namespace utility;
	using namespace core::pose;
	using namespace protocols::antibody;
	using namespace protocols::antibody::clusters;
	
typedef std::map< CDRNameEnum, CDRGraftInstructions > GraftInstructions;
typedef std::map< CDRNameEnum, CDRDesignInstructions > DesignInstructions;
typedef std::map< CDRNameEnum, vector1<PoseOP> > CDRSet;
typedef std::map< CDRNameEnum, vector1< CDRClusterEnum > > CDRClusterMap;
typedef std::map< CDRNameEnum, vector1< std::string > > PDBMap;

typedef std::map< core::chemical::AA, core::Real > AAProbabilities;

///@brief Class responsible for loading information from the AntibodyDatabase.  
class AntibodyDatabaseManager : public utility::pointer::ReferenceCount{
public: 
	
	AntibodyDatabaseManager();

	AntibodyDatabaseManager(std::string const database_path);
	
	virtual ~AntibodyDatabaseManager();
        
	///@brief Load data needed for GraftDesigner - CDR structures + their corresponding cluster.
	std::pair< CDRSet, CDRClusterMap >
	load_cdrs_for_grafting(AntibodyInfoCOP ab_info, GraftInstructions & instructions, PDBMap & pdbmap, core::Size overhang=3);
	
	///@brief Load probability data for CDR Design.  Returns CDRs where data could not be loaded. Cutoff indicates number of total sequences needed to use the data.
	vector1< CDRNameEnum >
	load_cdr_design_data(AntibodyInfoCOP ab_info, core::pose::Pose const & pose, std::map< core::Size, AAProbabilities > & prob_set, core::Size const cutoff, DesignInstructions & instructions);
	
	//std::map< core::Size, AAProbabilities >
	//load_framework_design_data(AntibodyInfoCOP ab_info, core::pose::Pose const & pose);
	
private:
	
	///@brief Gets database session.
	void
	start_database_session(std::string const database_path);
	
	
	///@brief  Checks to make sure the instructions make sense before trying to create the statement for the db..  
	void
	check_for_graft_instruction_inconsistencies(AntibodyInfoCOP ab_info, GraftInstructions & instructions);
	
	//protocols::features::ProteinSilentReportOP protein_silent_report_;
	std::string db_path_;
	utility::sql_database::sessionOP db_session_;
	
};
	
} 
}
}

#endif //INCLUDED_protocols_antibody_design_AntibodyDatabaseManager_hh
