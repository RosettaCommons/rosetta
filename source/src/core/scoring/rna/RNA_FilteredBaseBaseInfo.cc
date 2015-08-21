// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RNA_RawBaseBaseInfo.cc
/// @brief  Statistically derived RNA potential
/// @author Rhiju Das


// Unit headers
#include <core/scoring/rna/RNA_FilteredBaseBaseInfo.hh>
#include <core/scoring/rna/RNA_RawBaseBaseInfo.hh>
#include <core/pose/rna/RNA_BaseDoubletClasses.hh>
#include <core/scoring/rna/data/RNA_DataInfo.hh>
#include <core/chemical/rna/util.hh>

// Package headers

// Project headers
#include <core/chemical/AA.hh>

// Utility headers
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>


// C++

///////////////////////////////////////////////////////
// Keep track of some base geometry that is
// useful for RNA scoring.
///////////////////////////////////////////////////////
using namespace core::chemical::rna;

namespace core {
namespace scoring {
namespace rna {

using namespace ObjexxFCL::format;

RNA_FilteredBaseBaseInfo::RNA_FilteredBaseBaseInfo():
	total_base_pair_score_ ( 0.0 ),
	total_base_axis_score_ ( 0.0 ),
	total_base_stagger_score_ ( 0.0 ),
	total_base_stack_score_ ( 0.0 ),
	total_base_stack_axis_score_ ( 0.0 ),
	scale_axis_stagger_( true ),
	basepair_axis_stagger_scaling_( 0.1 ),
	basestack_axis_scaling_( 1.0 ),
	include_neighbor_base_stacks_( basic::options::option[ basic::options::OptionKeys::score::include_neighbor_base_stacks ]() ),
	calculated_( false ),
	rna_verbose_( false )
{}

/// @details Copy constructors must copy all data, not just some...
RNA_FilteredBaseBaseInfo::RNA_FilteredBaseBaseInfo( RNA_FilteredBaseBaseInfo const & src ) :
	CacheableData(),
	scale_axis_stagger_( src.scale_axis_stagger_ ),
	basepair_axis_stagger_scaling_( src.basepair_axis_stagger_scaling_ ),
	basestack_axis_scaling_( src.basestack_axis_scaling_ ),
	include_neighbor_base_stacks_( src.include_neighbor_base_stacks_ )
{
	filtered_base_pair_array_ = src.filtered_base_pair_array_;
	filtered_base_axis_array_ = src.filtered_base_axis_array_;
	filtered_base_stagger_array_ = src.filtered_base_stagger_array_;
	filtered_base_stack_array_ = src.filtered_base_stack_array_;
	filtered_base_stack_axis_array_ = src.filtered_base_stack_axis_array_;
	total_base_pair_score_ = src.total_base_pair_score_;
	total_base_axis_score_ = src.total_base_axis_score_;
	total_base_stagger_score_ = src.total_base_stagger_score_;
	total_base_stack_score_ = src.total_base_stack_score_;
	total_base_stack_axis_score_ = src.total_base_stack_axis_score_;
	calculated_ = src.calculated_;
	rna_verbose_ = src.rna_verbose_;
}

void
RNA_FilteredBaseBaseInfo::resize( Size const & total_residue )
{
	filtered_base_pair_array_.dimension( total_residue, total_residue );
	filtered_base_stagger_array_.dimension( total_residue, total_residue );
	filtered_base_axis_array_.dimension( total_residue, total_residue );
	filtered_base_stack_array_.dimension( total_residue, total_residue );
	filtered_base_stack_axis_array_.dimension( total_residue, total_residue );
}


//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Sort list of base pairs by energy, and go down list -- don't allow
// any base pairs that are mutually exclusive!
///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
void
RNA_FilteredBaseBaseInfo::carry_out_filtering( RNA_RawBaseBaseInfo const & raw_base_base_info ){

	// Actually not a good idea -- sometimes the Monte Carlo
	// wants a tally of the energy, but no recalculation of pair terms.
	//debug_assert( raw_base_base_info.calculated() );

	Size const total_residue = raw_base_base_info.size();
	resize( total_residue );

	figure_out_rna_base_pairs_to_score(  raw_base_base_info );

	figure_out_rna_base_stacks_to_score( raw_base_base_info );

	set_calculated( true );
}

void
RNA_FilteredBaseBaseInfo::figure_out_rna_base_pairs_to_score(
	RNA_RawBaseBaseInfo const & raw_base_base_info
)
{
	ObjexxFCL::FArray3D < Real > raw_base_pair_array( raw_base_base_info.base_pair_array() );
	ObjexxFCL::FArray3D < Real > raw_base_axis_array( raw_base_base_info.base_axis_array() );
	ObjexxFCL::FArray3D < Real > raw_base_stagger_array( raw_base_base_info.base_stagger_array() );

	ObjexxFCL::FArray2D < Real > raw_base_geometry_orientation_array( raw_base_base_info.base_geometry_orientation_array() );

	//A rigorous list of base pairs, for scoring.
	pose::rna::EnergyBasePairList energy_base_pair_list;

	scored_base_pair_list_.clear();

	filtered_base_pair_array_ = 0.0;
	filtered_base_axis_array_ = 0.0;
	filtered_base_stagger_array_ = 0.0;

	total_base_pair_score_ = 0.0;
	total_base_axis_score_ = 0.0;
	total_base_stagger_score_ = 0.0;

	Real const SCORE_CUTOFF = -0.001;

	Size const total_residue = raw_base_base_info.size();
	ObjexxFCL::FArray2D_bool edge_is_base_pairing( total_residue, NUM_EDGES, false );

	for ( Size i = 1; i <= total_residue; i++ ) {
		for ( Size k = 1; k <= NUM_EDGES; k++ ) {

			for ( Size j = i + 1; j <= total_residue; j++ ) {
				//A base pair is only real if each partner called the other one a true partner.
				if ( raw_base_pair_array( i, j, k ) < SCORE_CUTOFF ) {

					Size found_match( 0 );
					Real tmp_energy = 0.0;
					for ( Size m = 1; m <= NUM_EDGES; m++ ) {
						if ( raw_base_pair_array( j, i, m ) < tmp_energy ) {
							found_match = m;
							tmp_energy = raw_base_pair_array( j, i, m );
						}
					} //m
					if ( found_match == 0 ) continue;

					Real const total_base_pair_energy =
						raw_base_pair_array( i, j, k ) + raw_base_pair_array( j, i, found_match );

					pose::rna::BasePair base_pair;
					base_pair.res1 = i;
					base_pair.edge1 = k;

					base_pair.res2 = j;
					base_pair.edge2 = found_match;

					//orientations are cos( theta ) and should be symmetric!
					debug_assert( std::abs( raw_base_geometry_orientation_array( i, j ) - raw_base_geometry_orientation_array( j, i ) ) < 1.0e-2 );
					base_pair.orientation = ( raw_base_geometry_orientation_array( i, j ) + raw_base_geometry_orientation_array( j, i ) < 0.0 ? 1 : 2 );

					energy_base_pair_list.push_back( std::make_pair( total_base_pair_energy, base_pair )  );

				} //is it a basepair?
			} //j

		} //k
	} //i

	energy_base_pair_list.sort(); //Start with the lowest energy base pairs.

	// static bool const scale_axis_stagger_by_xy_score = true;

	for ( pose::rna::EnergyBasePairList::const_iterator it = energy_base_pair_list.begin();
			it != energy_base_pair_list.end(); ++it ) {
		Real const energy = it->first;
		pose::rna::BasePair const base_pair = it->second;

		Size const i = base_pair.res1;
		Size const k = base_pair.edge1;

		Size const j = base_pair.res2;
		Size const m = base_pair.edge2;

		if ( edge_is_base_pairing( i, k ) ) continue;
		if ( edge_is_base_pairing( j, m ) ) continue;

		edge_is_base_pairing( i, k ) = true;
		edge_is_base_pairing( j, m ) = true;

		Real scalefactor ( 1.0f ) ;
		if ( scale_axis_stagger_ ) scalefactor = -1.0 * basepair_axis_stagger_scaling_ * energy;

		Real const scaled_axis_energy = scalefactor * ( raw_base_axis_array( i, j, k ) + raw_base_axis_array( j, i, m ) );
		Real const scaled_stagger_energy = scalefactor * ( raw_base_stagger_array( i, j, k ) + raw_base_stagger_array( j, i, m ) );
		//rna_axis_score    += scaled_axis_energy;
		//rna_stagger_score += scaled_stagger_energy;

		if ( rna_verbose_ ) {

			std::cout << "BASE PAIR: " << I( 3, i ) << " " << I( 3, j ) << " "
				<< get_edge_from_num( k ) << " "
				<< get_edge_from_num( m ) << " "
				<< get_orientation_from_num( base_pair.orientation )
				<< "  ==  > "
				<< F( 5, 3, raw_base_pair_array( i, j, k ) ) << " " << F( 5, 3, raw_base_pair_array( j, i, m ) )  << " "
				//        << scaled_axis_energy
				<< std::endl;

		}

		//if ( user_defined_basepair(i,j,k,m) > 0.1 ) rna_contact_score += -1.0;

		filtered_base_pair_array_( i, j ) = energy;
		filtered_base_pair_array_( j, i ) = energy;
		total_base_pair_score_ += energy;

		filtered_base_axis_array_( i, j ) = scaled_axis_energy;
		filtered_base_axis_array_( j, i ) = scaled_axis_energy;
		total_base_axis_score_ += scaled_axis_energy;

		filtered_base_stagger_array_( i, j ) = scaled_stagger_energy;
		filtered_base_stagger_array_( j, i ) = scaled_stagger_energy;
		total_base_stagger_score_ += scaled_stagger_energy;

		scored_base_pair_list_.push_back( *it );
	}

}

///////////////////////////////////////////////////////////////////////////////
void
RNA_FilteredBaseBaseInfo::figure_out_rna_base_stacks_to_score(
	RNA_RawBaseBaseInfo const & raw_base_base_info )
{

	ObjexxFCL::FArray2D < Real > raw_base_stack_array( raw_base_base_info.base_stack_array() );
	ObjexxFCL::FArray2D < Real > raw_base_stack_axis_array( raw_base_base_info.base_stack_axis_array() );
	ObjexxFCL::FArray2D < Real > raw_base_geometry_orientation_array( raw_base_base_info.base_geometry_orientation_array() );
	ObjexxFCL::FArray2D < Real > raw_base_geometry_height_array( raw_base_base_info.base_geometry_height_array() );

	//  static bool const include_neighbor_base_stacks = false;//truefalseoption("include_neighbor_base_stacks");

	Size const total_residue = raw_base_base_info.size();

	scored_base_stack_list_.clear();

	filtered_base_stack_array_ = 0.0;
	filtered_base_stack_axis_array_ = 0.0;

	total_base_stack_score_ = 0.0;
	total_base_stack_axis_score_ = 0.0;

	for ( Size i = 1; i <= total_residue; i++ ) {
		for ( Size j = i + 1; j <= total_residue; j++ ) {


			//Base pairs and base stacks are mutually exclusive!
			if ( filtered_base_pair_array_( i, j ) < 0.0 ) continue;

			//Both partners should think they are stacked onto each other.
			if ( raw_base_stack_array( i, j ) < 0.0 &&
					raw_base_stack_array( j, i ) < 0.0 ) {
				if ( rna_verbose_ ) std::cout << "BASE STACK: " << i << " " << j << std::endl;

				pose::rna::BaseStack base_stack;
				base_stack.res1 = i;
				base_stack.res2 = j;

				//orientations are cos( theta ) and should be symmetric!
				debug_assert( std::abs( raw_base_geometry_orientation_array( i, j ) - raw_base_geometry_orientation_array( j, i ) ) < 1.0e-2 );
				base_stack.orientation = ( raw_base_geometry_orientation_array( i, j ) + raw_base_geometry_orientation_array( j, i ) ) < 0.0 ? 1 : 2;

				// height is not necessarily (anti)-symmetric if the planes of the two bases aren't co-planar.
				base_stack.which_side = ( raw_base_geometry_height_array( i, j ) > 0.0 ) ? 1 : 2;

				Real const total_base_stack_energy = raw_base_stack_array( i, j ) + raw_base_stack_array( j, i );
				scored_base_stack_list_.push_back( std::make_pair( total_base_stack_energy,  base_stack ) );

				//By default, don't count stacks between neighboring nucleotides, since that
				// interaction is captured by fragments.
				if ( !include_neighbor_base_stacks_  &&  j == i + 1 ) continue;

				filtered_base_stack_array_( i, j ) = total_base_stack_energy;
				filtered_base_stack_array_( j, i ) = total_base_stack_energy;
				total_base_stack_score_ += total_base_stack_energy;

				Real total_base_stack_axis_energy =  raw_base_stack_axis_array( i, j ) + raw_base_stack_axis_array( j, i ) ;

				// This scaling is actually unity, unless we're fading near the boundaries:
				if ( scale_axis_stagger_ ) total_base_stack_axis_energy *= -1 * basestack_axis_scaling_ * total_base_stack_energy;

				filtered_base_stack_axis_array_( i, j ) = total_base_stack_axis_energy;
				filtered_base_stack_axis_array_( j, i ) = total_base_stack_axis_energy;
				total_base_stack_axis_score_ += total_base_stack_axis_energy;

			}
		}
	}

}

/////////////////////////////////////////////////////////////////////
Real RNA_FilteredBaseBaseInfo::get_data_score( data::RNA_DataInfo const & rna_data_info ) const
{
	Real rna_data_score( 0.0 );

	utility::vector1< data::RNA_Datum > const & rna_data( rna_data_info.rna_data() );

	for ( Size n = 1; n <= rna_data.size(); n++ ) {

		Size const seqpos( rna_data[n].position() );

		// This could be accelerated if needed.
		for ( pose::rna::EnergyBasePairList::const_iterator it = scored_base_pair_list_.begin();
				it != scored_base_pair_list_.end(); ++it ) {

			pose::rna::BasePair const base_pair = it->second;

			Size const i = base_pair.res1;
			Size const j = base_pair.res2;

			Size k( 0 );
			if ( i == seqpos ) {
				k = base_pair.edge1;
			} else if ( j == seqpos ) {
				k = base_pair.edge2;
			} else {
				continue;
			}

			if ( k == rna_data[ n ].edge() ) {
				rna_data_score += rna_data[n].weight();
			}

		}
	}

	return rna_data_score;
}

} //rna
} //scoring
} //core
