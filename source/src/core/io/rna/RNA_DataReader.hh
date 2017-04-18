// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_rna_RNA_DataReader_hh
#define INCLUDED_protocols_rna_RNA_DataReader_hh

#include <core/io/rna/RNA_DataReader.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/full_model_info/FullModelParameters.fwd.hh>
#include <core/pose/rna/RNA_DataInfo.hh>
#include <core/io/rna/RDAT.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>

// C++ Headers
#include <string>
#include <iostream>

#include <utility/vector1.hh>


namespace core {
namespace io {
namespace rna {


//////////////////////////////////////////////////////////////////////////////////////////////
class RNA_DataReader : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~RNA_DataReader();

	//constructor
	RNA_DataReader( std::string const & rna_data_file );

	void
	initialize( std::string const & rna_data_file );

	void
	fill_rna_data_info( core::pose::Pose & pose );

	bool
	has_reactivities(){ return ( rna_data_info_with_conventional_numbering_ != 0 &&
		rna_data_info_with_conventional_numbering_->rna_reactivities().size() > 0 ); }

private:

	void
	read_backbone_info( std::istringstream & line_stream,
		utility::vector1< Size > & backbone_res );

	void
	read_data_info( std::istringstream & line_stream );

	void
	read_data_from_file( std::string const & rna_data_file );

	void
	read_reactivity_info( std::istringstream & line_stream, pose::rna::RNA_ReactivityType const type );

	void
	get_reactivity_from_rdat( core::io::rna::RDAT const & rdat,
		pose::rna::RNA_ReactivityType const & type,
		std::string const & modifier_name );
	void
	read_data_from_rdat( std::string const & filename );

	ObjexxFCL::FArray1D< bool >
	fill_backbone_array( utility::vector1< Size > const & backbone_res,
		core::pose::Pose const & pose );

private:

	pose::rna::RNA_DataInfoOP rna_data_info_with_conventional_numbering_;
	utility::vector1< Size > backbone_burial_res_, backbone_exposed_res_;

};

pose::rna::RNA_DataInfo const &
get_rna_data_info( core::pose::Pose & pose, std::string const & rna_data_file,
	core::scoring::ScoreFunctionOP scorefxn = 0 );


} //rna
} //io
} //core

#endif
