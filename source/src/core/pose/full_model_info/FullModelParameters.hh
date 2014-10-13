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
// These include
//
//   full_sequence
//   'conventional' numbering/chain scheme,
//   cutpoint_open_in_full_model, etc.
//
// See FullModelParameterType for full list of variables.
//
// Note that there are is no information here on what subset of
//  residues a specific pose contains (thats in FullModelInfo).
//
// The variables in here really should be 'permanent' -- parameters that won't
//  change during a run.
//
// Note that integer lists are stored in two ways, for convenience:
//
//  parameter_values_at_res
//    [ 0, 0, 1, 1, 0, 0, 2, 2 ] (has same size as full_sequence)
//
//  parameter_values_as_res_lists -- same info as above, different format.
//    { 0:[1, 2, 5, 6],  1:[3, 4], 2:[7, 8] }
//
// -- rhiju, 2014
//
////////////////////////////////////////////////////////////////////////////

namespace core {
namespace pose {
namespace full_model_info {

	class FullModelParameters: public utility::pointer::ReferenceCount {

	public:

		//constructor
		FullModelParameters();

		FullModelParameters( std::string const full_sequence );

		FullModelParameters( std::string const full_sequence,
												 utility::vector1< Size > const & cutpoint_open_in_full_model,
												 utility::vector1< Size > const & res_numbers_in_pose	);

		FullModelParameters( pose::Pose const & pose,
												 utility::vector1< Size > & res_list /*will be updated*/ );

		FullModelParameters( FullModelParameters const & src );

		//destructor
		~FullModelParameters();

	public:

		FullModelParametersOP
		clone() const
		{
			return FullModelParametersOP( new FullModelParameters( *this ) );
		}

		std::string const & full_sequence() const { return full_sequence_;}

		utility::vector1< int >  const & conventional_numbering() const { return conventional_numbering_;}
		utility::vector1< char > const & conventional_chains() const { return conventional_chains_;}

		void set_conventional_numbering( utility::vector1< int > const & setting ) { conventional_numbering_  = setting; }
		void set_conventional_chains( utility::vector1< char > const & setting ) { conventional_chains_  = setting; }

		// this is res_at_value
		void
		set_parameter( FullModelParameterType const type,
									 utility::vector1< Size > const & setting );

		void
		set_parameter_as_res_list( FullModelParameterType const type,
															 utility::vector1< Size > const & setting );

		void
		set_parameter_as_res_lists( FullModelParameterType const type,
																std::map< Size, utility::vector1< Size > > const & setting );

		utility::vector1< Size > const &
		get_res_list( FullModelParameterType const type, Size const value ) const;

		utility::vector1< Size > const &
		get_res_list( FullModelParameterType const type ) const { return get_res_list( type, 1 ); }

		utility::vector1< Size > const &
		get_parameter( FullModelParameterType const type ) const;

		std::map< Size, utility::vector1< Size > > const &
		get_parameter_as_res_lists( FullModelParameterType const type ) const;

		utility::vector1< Size >
		conventional_to_full( utility::vector1< int > const & res_list ) const;

		utility::vector1< Size >
		conventional_to_full( std::pair< utility::vector1< int >, utility::vector1< char > > const & resnum_and_chain ) const;

		bool
		has_conventional_residue( int const res_num ) const;

		bool
		has_conventional_residue( int const res_num, char const chain ) const;

		Size
		conventional_to_full( int const res_num ) const;

		Size
		conventional_to_full( int const res_num, char const chain ) const;

		utility::vector1< int >
		full_to_conventional( utility::vector1< Size > const & res_list ) const;

		int
		full_to_conventional( Size const res_num ) const;

		utility::vector1< Size >
		chains_in_full_model() const;

		Size size() const { return full_sequence_.size(); }

	private:

		void
		fill_parameter_values( utility::vector1< Size > & parameter_values_at_res,
													 Size const idx, utility::vector1< Size > const & res_list ) const;

		std::map< Size, utility::vector1< Size > >
		convert_to_res_lists_by_value( utility::vector1< Size > const & parameter_values_at_res );

		utility::vector1< Size >
		convert_to_parameter_values_at_res( utility::vector1< Size > const & res_list );

		utility::vector1< Size >
		convert_to_parameter_values_at_res( std::map< Size, utility::vector1< Size > > const & res_lists );

		void
		get_sequence_with_gaps_filled_with_n( pose::Pose const & pose,
																					std::string & sequence,
																					utility::vector1< int  > & conventional_numbering,
																					utility::vector1< char > & conventional_chains,
																					utility::vector1< Size > & res_list ) const;

		utility::vector1< Size >
		get_cutpoint_open_from_pdb_info( pose::Pose const & pose ) const;

		void
		keep_chain_and_cutpoint_open_matched( FullModelParameterType const & type );

		/// @brief input operator
		friend std::istream & operator >>(std::istream & is, FullModelParameters & t);

		/// @brief output operator
		friend std::ostream & operator <<(std::ostream & os, FullModelParameters const & t);

		/// @brief equal to operator
		friend
		bool
		operator==(
							 FullModelParameters const & a,
							 FullModelParameters const & b
							 );

		/// @brief not equal to operator
		friend
		bool
		operator!=(
							 FullModelParameters const & a,
							 FullModelParameters const & b
							 );

	private:

		std::string full_sequence_;
		utility::vector1< int >  conventional_numbering_; // permits user to use numbering other than 1, 2, 3...
		utility::vector1< char > conventional_chains_;    // permits user to use chains other than A, B, C, ..
		std::map< Size, std::string > non_standard_residues_; // for DNA, non-natural protein/RNA, ligands, ions, etc.

		std::map< FullModelParameterType, utility::vector1< Size > > parameter_values_at_res_;
		// this is set at the same time as above.
		std::map< FullModelParameterType, std::map< Size, utility::vector1< Size > > > parameter_values_as_res_lists_;

	};

} //full_model_info
} //pose
} //core

#endif
