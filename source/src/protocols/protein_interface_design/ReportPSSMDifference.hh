// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file devel/protein_interface_design/ReportPSSMDifference.hh
/// @brief calculation of the difference in PSSM score between mutated and native pose
/// @author Hermann Zellner (hermann1.zellner@biologie.uni-regensburg.de)


#ifndef INCLUDED_protocols_protein_interface_design_ReportPSSMDifference_hh
#define INCLUDED_protocols_protein_interface_design_ReportPSSMDifference_hh

#include <core/types.hh>

#include <utility/vector1.hh>

#include <map>

#include <core/chemical/AA.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.fwd.hh>


namespace protocols {
namespace protein_interface_design {


class ReportPSSMDifferences
{
public:
	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;
public:
	ReportPSSMDifferences( ) {};

	ReportPSSMDifferences( ReportPSSMDifferences const & init ) : // copy constructor
		res_name1_( init.res_name1_ ),
		pssm_data_( init.pssm_data_ )
	{}
	core::Real calculate( Pose const & pose1, Pose const & pose2, core::pack::task::PackerTaskCOP const & task );

	std::map< Size, std::string > const & res_name1() const { return res_name1_; }

	utility::vector1< std::pair< core::chemical::AA, utility::vector1< Real > > > const & pssm_() const { return pssm_data_; }

	bool load_pssm_data(std::string const & native_filename);

	virtual ~ReportPSSMDifferences() = default;
private:
	std::map< Size, std::string > res_name1_;
	utility::vector1< std::pair< core::chemical::AA, utility::vector1< Real > > > pssm_data_;
};

} //protocols
} //protein_interface_design

#endif /* INCLUDED_protocols_protein_interface_design_REPORTPSSMDIFFERENCE_HH_ */
