// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RNA_VDW_BinChecker.hh
/// @brief
/// @detailed
///
/// @author Parin Sripkaddevong


#ifndef INCLUDED_protocols_stepwise_rna_RNA_VDW_BinChecker_hh
#define INCLUDED_protocols_stepwise_rna_RNA_VDW_BinChecker_hh



#include <protocols/stepwise/modeler/rna/StepWiseRNA_Classes.hh> /*For core::Size and Torsion_Info*/
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <protocols/stepwise/modeler/rna/checker/VDW_RepScreenInfo.hh>
#include <protocols/stepwise/modeler/rna/checker/VDW_CachedRepScreenInfo.fwd.hh>
#include <protocols/stepwise/modeler/rna/checker/VDW_Grid.hh>



#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <utility/vector1.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <string>
#include <map>
#include <core/chemical/AA.hh>

#include <core/conformation/Residue.hh>
#include <set>
#include <fstream>
typedef  numeric::xyzMatrix< core::Real > Matrix;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace checker {


// Atom_Bin is now declared in VDW_Grid.hh
// struct Atom_Bin{ int x; int y; int z; };

class RNA_VDW_BinChecker: public utility::pointer::ReferenceCount {

public:

	//constructor that is initialized with pose!
	RNA_VDW_BinChecker( core::pose::Pose const & pose );

	//constructor!
	RNA_VDW_BinChecker();

	//destructor -- necessary?
	virtual ~RNA_VDW_BinChecker();

	void
	FARFAR_setup_using_user_input_VDW_pose( utility::vector1< std::string > const & VDW_rep_screen_pose_info, core::pose::Pose const & const_working_pose );


	void
	setup_using_user_input_VDW_pose( utility::vector1< std::string > const & All_VDW_rep_screen_pose_info, core::pose::Pose const & const_working_pose, working_parameters::StepWiseWorkingParametersCOP const & working_parameters );


	void
	setup_using_working_pose( core::pose::Pose const & const_working_pose, working_parameters::StepWiseWorkingParametersCOP const & working_parameters );


	void
	update_VDW_screen_bin( core::pose::Pose const & pose,
		utility::vector1< core::Size > const & ignore_res_list,
		bool const is_prepend,  //associated with ignore_res_list
		std::ofstream & outfile_act );

	void
	create_VDW_screen_bin( core::pose::Pose const & pose,
		utility::vector1< core::Size > const & ignore_res_list,
		bool const is_prepend,  //associated with ignore_res_list
		numeric::xyzVector< core::Real > const & reference_xyz,
		bool const verbose = false );
	void
	create_VDW_screen_bin( utility::vector1< VDW_RepScreenInfo > const & VDW_rep_screen_info_list,
		numeric::xyzVector< core::Real > const & reference_xyz,
		bool const verbose );


	void
	create_VDW_screen_bin( utility::vector1< core::pose::Pose > const & pose_list,
		utility::vector1< utility::vector1< core::Size > > const & list_of_ignore_res_list,
		utility::vector1< bool > const list_of_is_prepend,
		numeric::xyzVector< core::Real > const & reference_xyz,
		bool const verbose );



	//fast version
	bool
	VDW_rep_screen( core::pose::Pose const & screening_pose, //Warning..this pose coordinate is not update...use here for VIRTUAL atom screening.
		core::Size const & moving_res,
		core::conformation::Residue const & rsd_at_origin,
		core::kinematics::Stub const & moving_res_base_stub ); //the actual updated coordiate

	//Slow version (in the sense that position of screening_pose had to be updated before this function is called)..
	bool
	VDW_rep_screen( core::pose::Pose const & screening_pose, core::Size const & moving_res );

	bool
	VDW_rep_screen_with_act_pose( core::pose::Pose const & screening_pose, utility::vector1< core::Size > const & moving_res_list, bool const local_verbose = true );


	bool
	user_inputted_VDW_screen_pose() const;

	void
	reference_xyz_consistency_check( numeric::xyzVector< core::Real > const & inputted_reference_xyz ) const;

	void
	set_VDW_rep_alignment_RMSD_CUTOFF( core::Real const & setting ){ VDW_rep_alignment_RMSD_CUTOFF_ = setting ; }

	void
	set_VDW_rep_delete_matching_res( utility::vector1< std::string > const & setting ){ VDW_rep_ignore_matching_res_ = setting ; }

	//void
	//align_to_first_working_pose(core::pose::Pose & pose, std::string const & tag) const;

	void
	set_physical_pose_clash_dist_cutoff( core::Real const & setting ){ physical_pose_clash_dist_cutoff_ = setting ; }

	void
	set_output_pdb( bool const setting ){ output_pdb_ = setting; }

private:

	void
	check_VDW_screen_bin_is_setup() const;

	bool
	is_atom_bin_in_range( Atom_Bin const & atom_pos_bin ) const;

	bool
	check_atom_bin_in_range( Atom_Bin const & atom_pos_bin );

	core::Vector
	get_reference_xyz( core::pose::Pose const & pose, core::Size const reference_res/*,bool const verbose = false*/ );

	core::Vector
	get_reference_xyz_average( core::pose::Pose const & pose );

	void
	set_reference_xyz( numeric::xyzVector< core::Real > const & reference_xyz );

	Atom_Bin
	get_atom_bin( numeric::xyzVector< core::Real > const & atom_pos ) const;

	numeric::xyzVector< core::Real >
	get_atom_pos( Atom_Bin const & atom_bin ) const;

	void
	output_atom_bin( std::string const filename ) const;


	/*
	void
	delete_matching_res_in_VDW_rep_screen_pose( core::pose::Pose & VDW_rep_screen_pose,
	core::pose::Pose const & working_pose,
	utility::vector1< core::Size > const & VDW_rep_screen_align_res,
	utility::vector1< core::Size > const & working_align_res,
	std::map< core::Size, core::Size > & full_to_sub,
	bool const verbose ) const;
	*/

	//replacement for delete_matching_res_in_VDW_rep_screen_pose()!
	utility::vector1< core::Size >
	get_matching_res_in_VDW_rep_screen_pose( core::pose::Pose const & VDW_rep_screen_pose,
		core::pose::Pose const & working_pose,
		utility::vector1< core::Size > const & VDW_rep_screen_align_res,
		utility::vector1< core::Size > const & working_align_res,
		std::map< core::Size, core::Size > & full_to_sub ) const;


	void
	align_VDW_rep_screen_pose( core::pose::Pose & VDW_rep_screen_pose,
		core::pose::Pose const & working_pose,
		utility::vector1< core::Size > const & VDW_rep_screen_align_res,
		utility::vector1< core::Size > const & working_align_res,
		bool const verbose ) const;

	void
	read_in_VDW_rep_screen_pose( VDW_RepScreenInfo & VDW_rep_screen_info ) const;


	void
	create_VDW_rep_screen_pose( VDW_RepScreenInfo & VDW_rep_screen_info, //This function update this class!
		core::pose::Pose const & working_pose,
		std::map< core::Size, core::Size > & full_to_sub,
		bool const verbose ) const;


private:


	core::Real max_distance_;
	core::Real atom_bin_size_;

	int const bin_min_;
	int const bin_max_;
	int const bin_offset_;
	core::Size const num_clash_atom_cutoff_;
	bool const write_to_file_;
	bool is_reference_xyz_setup_;
	bool is_VDW_screen_bin_setup_;
	bool user_inputted_VDW_screen_pose_;

	VDW_GridCOP VDW_screen_bin_;
	utility::vector1< Atom_Bin > occupied_xyz_bins_;
	numeric::xyzVector< core::Real > reference_xyz_;

	core::Real VDW_rep_alignment_RMSD_CUTOFF_;
	bool tolerate_off_range_atom_bin_;
	int num_atom_pos_bin_out_of_range_message_outputted_;

	bool VDW_rep_screen_with_physical_pose_verbose_;
	core::Real physical_pose_clash_dist_cutoff_;

	utility::vector1< std::string > VDW_rep_ignore_matching_res_;
	bool use_VDW_rep_pose_for_screening_;

	utility::vector1< VDW_RepScreenInfo > VDW_rep_screen_info_list_;

	bool output_pdb_;
	bool optimize_memory_usage_;
	bool optimize_speed_;
	bool verbose_;


};

} //checker
} //rna
} //modeler
} //stepwise
} //protocols

#endif
