// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/rna/sugar/VirtualSugarSamplerFromStringList.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_rna_VirtualSugarSamplerFromStringList_HH
#define INCLUDED_protocols_stepwise_rna_VirtualSugarSamplerFromStringList_HH

#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <protocols/stepwise/modeler/rna/sugar/VirtualSugarSamplerFromStringList.fwd.hh>
#include <protocols/stepwise/modeler/rna/sugar/SugarModeling.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

using namespace core;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace sugar {

class VirtualSugarSamplerFromStringList: public protocols::moves::Mover {

public:

	//constructor
	VirtualSugarSamplerFromStringList( working_parameters::StepWiseWorkingParametersCOP & working_parameters,
		utility::vector1< std::string > const sample_virtual_sugar_string_list);

	//destructor
	~VirtualSugarSamplerFromStringList();

	virtual void apply( pose::Pose & pose_to_visualize );

	virtual std::string get_name() const;

public:

	void set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn );

	void set_use_phenix_geo( bool const & setting ) { use_phenix_geo_ = setting; }

	void set_legacy_mode( bool const & setting ) { legacy_mode_ = setting; }

	void set_choose_random( bool const & setting ) { choose_random_ = setting; }

	void set_integration_test_mode( bool const & setting ){ integration_test_mode_ = setting; }

	void set_tag( std::string const & setting ){ tag_ = setting; }

	void set_silent_file_out( std::string const & setting ){ silent_file_out_ = setting; }

private:

	utility::vector1< SugarModeling >
	setup_sugar_modeling_list( pose::Pose const & pose ) const;

	bool
	empty_sugar_modeling_list( utility::vector1< SugarModeling > const &  sugar_modeling_list );

	bool
	empty_pose_data_list( utility::vector1< pose::PoseOP > const & pose_list, Size const n, std::string tag ) ;

	void
	output_pose_data_list( utility::vector1< pose::PoseOP > & pose_data_list );

private:

	working_parameters::StepWiseWorkingParametersCOP working_parameters_; //need to use the full_to_sub map...should convert to const style.. Parin Feb 28, 2010
	utility::vector1< std::string > const sample_virtual_sugar_string_list_;

	std::string silent_file_out_;
	bool use_phenix_geo_;
	bool legacy_mode_;
	bool choose_random_;
	bool integration_test_mode_;
	std::string tag_;

	scoring::ScoreFunctionOP scorefxn_;

};

} //sugar
} //rna
} //modeler
} //stepwise
} //protocols

#endif
