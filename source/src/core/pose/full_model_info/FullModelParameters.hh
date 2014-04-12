// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/pose/full_model_info/FullModelParameters.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_pose_full_model_info_FullModelParameters_HH
#define INCLUDED_core_pose_full_model_info_FullModelParameters_HH

#include <utility/pointer/ReferenceCount.hh>
#include <core/pose/full_model_info/FullModelParameterType.hh>
#include <core/pose/full_model_info/FullModelParameters.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <map>
#include <string>

using namespace core;

////////////////////////////////////////////////////////////////////////////
// This object is held in FullModelInfo, and stores all information
//  related to the eventual full-length model.
//
// These include full_sequence, and 'conventional' numbering/chain scheme,
//  and cutpoint_open_in_full_model, etc.
// See FullModelParameterType for full list of variables.
//
// These really should be 'permanent' -- parameters that won't
//  change during a run.
//
// -- rhiju, 2014
////////////////////////////////////////////////////////////////////////////

namespace core {
namespace pose {
namespace full_model_info {

	class FullModelParameters: public utility::pointer::ReferenceCount {

	public:

		//constructor
		FullModelParameters( std::string const full_sequence );

		FullModelParameters( std::string const full_sequence,
												 utility::vector1< Size > const & cutpoint_open_in_full_model,
												 utility::vector1< Size > const & res_numbers_in_pose	);

		FullModelParameters( pose::Pose const & pose,
												 utility::vector1< Size > const & res_list );

		FullModelParameters( FullModelParameters const & src );

		//destructor
		~FullModelParameters();

	public:

		FullModelParametersOP
		clone() const
		{
			return new FullModelParameters( *this );
		}

		std::string const & full_sequence() const { return full_sequence_;}
		utility::vector1< int > const & conventional_numbering() const { return conventional_numbering_;}

		// this is res_at_value
		void
		set_parameter( FullModelParameterType const type,
									 utility::vector1< Size > const & setting );

		void
		set_parameter_as_res_list( FullModelParameterType const type,
															 utility::vector1< Size > const & setting );

		utility::vector1< Size > const &
		get_res_list( FullModelParameterType const type, Size const value ) const;

		utility::vector1< Size > const &
		get_res_list( FullModelParameterType const type ) const { return get_res_list( type, 1 ); }

		utility::vector1< Size > const &
		get_parameter( FullModelParameterType const type ) const;

		Size
		get_value_at_res( Size const n, FullModelParameterType const type );

		Size size() const { return full_sequence_.size(); }

	private:

		std::map< Size, utility::vector1< Size > >
		convert_to_res_lists_by_value( utility::vector1< Size > const & parameter_values_at_res );

		utility::vector1< Size >
		convert_to_parameter_values_at_res( utility::vector1< Size > const & res_list );

		void
		get_sequence_with_gaps_filled_with_n( pose::Pose const & pose,
																					std::string & sequence,
																					utility::vector1< int > & full_numbering,
																					utility::vector1< Size > const & res_list ) const;

		utility::vector1< Size >
		get_cutpoint_open_from_pdb_info( pose::Pose const & pose ) const;

		utility::vector1< Size >
		chains_in_full_model() const;

		void
		keep_chain_and_cutpoint_open_matched( FullModelParameterType const & type );

	private:

		std::map< FullModelParameterType, utility::vector1< Size > > parameter_values_at_res_;
		// this is set at the same time as above.
		std::map< FullModelParameterType, std::map< Size, utility::vector1< Size > > > parameter_values_as_res_lists_;

		std::string full_sequence_;
		utility::vector1< int > conventional_numbering_; // permits user to start numbering without beginning at 1.
		// std::string conventional_chains_; // not defined yet!

	};

} //full_model_info
} //pose
} //core

#endif
