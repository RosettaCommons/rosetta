// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/H3CterInsert.hh
/// @brief Build a homology model of an antibody
/// @details
///
///
/// @author Jianqing Xu ( xubest@gmail.com )


#ifndef INCLUDED_protocols_antibody_H3CterInsert_hh
#define INCLUDED_protocols_antibody_H3CterInsert_hh


#include <protocols/antibody/H3CterInsert.fwd.hh>
#include <core/fragment/FragData.hh>
#include <protocols/antibody/AntibodyInfo.fwd.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/moves/Mover.hh>


namespace protocols {
namespace antibody {


//////////////////////////////////////////////////////////////////////////
/// @brief H3 CDR, Fragment Insertion and CCD
/// @details
class H3CterInsert : public protocols::moves::Mover {

public:
	/// @brief default constructor
	H3CterInsert();

	/// @brief constructor with arguments
	H3CterInsert(antibody::AntibodyInfoOP antibody_info, bool camelid );


	/// @brief default destructor
	~H3CterInsert() override;

	void set_default();

	void apply(core::pose::Pose & pose ) override;
	std::string get_name() const override;

	// read CDR H3 C-terminal fragments (size: 4)
	void read_H3_cter_fragment();


private:

	// CDR H3 C-terminal fragments
	utility::vector1< core::fragment::FragData > H3_base_library_;

	AntibodyInfoOP ab_info_;

	bool user_defined_;

	/// @brief benchmark flag
	bool benchmark_;

	/// @brief is camelid antibody without light chain
	bool is_camelid_;


	void init(AntibodyInfoOP antibody_info, bool camelid, bool benchmark);
	//    void setup_objects();
	//    void finalize_setup( core::pose::Pose & pose );


	std::string H3_ter_library_filename_;

};


}//antibody
}//protocols

#endif

