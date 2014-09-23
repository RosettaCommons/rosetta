// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/upstream/ProteinUpstreamBuilder.cc
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

// Unit headers
#include <protocols/match/upstream/ProteinUpstreamBuilder.hh>

// Package headers
#include <protocols/match/BumpGrid.hh>
#include <protocols/match/downstream/DownstreamAlgorithm.hh>
// AUTO-REMOVED #include <protocols/match/downstream/DownstreamBuilder.hh>
#include <protocols/match/upstream/ProteinSCSampler.hh>
#include <protocols/match/upstream/OriginalScaffoldBuildPoint.hh>
#include <protocols/match/upstream/UpstreamResTypeGeometry.hh>

// Project headers
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/dunbrack/SingleResidueRotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/LexicographicalIterator.hh>

#include <protocols/match/Hit.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace upstream {

/// @details Auto-generated virtual destructor
FullChiSampleSet::~FullChiSampleSet() {}

static thread_local basic::Tracer TR( "protocols.match.upstream.ProteinUpstreamBuilder" );

/// dummy return value
BuildSet b;


/*enum ChiStrategy {
	rotameric_chi_follow_EX_flags,
	rotameric_chi_mimic_EX_flags,
	rotameric_chi_step_by_value,
	rotameric_chi_step_wi_sd_range,
	rotameric_chi_partition_sd_range,
	nonrotameric_chi_sample_wi_nrchi_bin,
	nonrotameric_chi_sample_wi_nrchi_bin_to_lower_boundary
};*/

SampleStrategyData::SampleStrategyData() :
	strategy_( follow_EX_flags ),
	sample_level_( core::pack::task::NO_EXTRA_CHI_SAMPLES ),
	step_size_( 0.0 ),
	sd_range_( 0.0 ),
	n_samples_wi_sd_range_( 0 ),
	nrchi_prob_minimum_for_extra_samples_( 0.2 ),
	n_samples_per_side_of_nrchi_bin_( 0 )
{}


SampleStrategyData::~SampleStrategyData() {}

void
SampleStrategyData::set_strategy( ChiStrategy setting )
{
	strategy_ = setting;
}

void
SampleStrategyData::set_sample_level( ExtraRotSample setting )
{
	assert( strategy_ == rotameric_chi_mimic_EX_flags );
	sample_level_ = setting;
}

void
SampleStrategyData::set_step_size( Real setting )
{
	assert(
		strategy_ == rotameric_chi_step_by_value ||
		strategy_ == rotameric_chi_step_wi_sd_range );

	step_size_ = setting;
}

void
SampleStrategyData::set_sd_range(  Real setting ) {
	assert(
		strategy_ == rotameric_chi_step_wi_sd_range ||
		strategy_ == rotameric_chi_partition_sd_range );
	sd_range_ = setting;
}

void
SampleStrategyData::set_n_samples_wi_sd_range( Size setting )
{
	assert( strategy_ == rotameric_chi_partition_sd_range );
	n_samples_wi_sd_range_ = setting;
}

void
SampleStrategyData::set_nrchi_prob_minimum_for_extra_samples( Real setting )
{
	assert(
		strategy_ == follow_EX_flags ||
		strategy_ == nonrotameric_chi_sample_wi_nrchi_bin ||
		strategy_ == nonrotameric_chi_sample_wi_nrchi_bin_to_lower_boundary );
	nrchi_prob_minimum_for_extra_samples_ = setting;
}

void
SampleStrategyData::set_n_samples_per_side_of_nrchi_bin( Size setting ) {
	assert(
		strategy_ == nonrotameric_chi_sample_wi_nrchi_bin ||
		strategy_ == nonrotameric_chi_sample_wi_nrchi_bin_to_lower_boundary );
	n_samples_per_side_of_nrchi_bin_ = setting;
}

ChiStrategy
SampleStrategyData::strategy() const { return strategy_; }

SampleStrategyData::ExtraRotSample
SampleStrategyData::sample_level() const
{
	assert( strategy_ == rotameric_chi_mimic_EX_flags );
	return sample_level_;
}


SampleStrategyData::Real
SampleStrategyData::step_size() const
{
	assert(
		strategy_ == rotameric_chi_step_by_value ||
		strategy_ == rotameric_chi_step_wi_sd_range );
	return step_size_;
}

SampleStrategyData::Real
SampleStrategyData::sd_range() const
{
	assert(
		strategy_ == rotameric_chi_step_wi_sd_range ||
		strategy_ == rotameric_chi_partition_sd_range );
	return sd_range_;
}

SampleStrategyData::Size
SampleStrategyData::n_samples_wi_sd_range() const
{
	assert( strategy_ == rotameric_chi_partition_sd_range );
	return n_samples_wi_sd_range_;
}

SampleStrategyData::Real
SampleStrategyData::nrchi_prob_minimum_for_extra_samples() const
{
	assert(
		strategy_ == follow_EX_flags ||
		strategy_ == nonrotameric_chi_sample_wi_nrchi_bin ||
		strategy_ == nonrotameric_chi_sample_wi_nrchi_bin_to_lower_boundary );
	return nrchi_prob_minimum_for_extra_samples_;
}

SampleStrategyData::Size
SampleStrategyData::n_samples_per_side_of_nrchi_bin() const
{
	assert(
		strategy_ == nonrotameric_chi_sample_wi_nrchi_bin ||
		strategy_ == nonrotameric_chi_sample_wi_nrchi_bin_to_lower_boundary );
	return n_samples_per_side_of_nrchi_bin_;
}

/*UpstreamResTypeGeometry::UpstreamResTypeGeometry() :
	restype_name_( "UNINITIALIZED" )
{}

UpstreamResTypeGeometry::UpstreamResTypeGeometry( core::chemical::ResidueType const & res ) :
	restype_name_( res.name() )
{
	initialize_from_residue_type( res );
}

void
UpstreamResTypeGeometry::initialize_from_residue_type(
	core::chemical::ResidueType const & res
)
{
	if ( restype_name_ != res.name() ) {
		restype_name_ = res.name();
	}

	/// 1. Resize arrays that depend on the number of atoms
	Size const n_atoms = res.natoms();

	controlling_chi_for_atom_ =  res.last_controlling_chi();
	which_point_for_atom_.resize( n_atoms );
	std::fill( which_point_for_atom_.begin(), which_point_for_atom_.end(), 0 );

	/// 2. Resize arrays that depend on the number of chi
	Size const n_chi = res.nchi();

	chitip_atoms_.resize( n_chi );
	std::fill( chitip_atoms_.begin(), chitip_atoms_.end(), 0 );

	ht_for_chitip_atoms_.resize( n_chi );
	for ( Size ii = 1; ii <= n_chi; ++ii ) ht_for_chitip_atoms_[ ii ].set_identity();

	nonchitip_atoms_.resize( res.nchi() );
	for ( Size ii = 1; ii <= n_chi; ++ii ) nonchitip_atoms_[ ii ].clear();

	points_for_nonchitip_atoms_.resize( res.nchi() );
	for ( Size ii = 1; ii <= n_chi; ++ii ) points_for_nonchitip_atoms_[ ii ].clear();

	if ( nchi() == 0 ) return; // match from gly? can't see why you'd want to!


	for ( Size ii = 1; ii <= n_chi; ++ii ) {
		assert( res.chi_atoms( ii ).size() == 4 );

		Size const
			chiat2( res.chi_atoms( ii )[ 2 ] ),
			chiat3( res.chi_atoms( ii )[ 3 ] ),
			chiat4( res.chi_atoms( ii )[ 4 ] );

		chitip_atoms_[ ii ] = chiat4;


		ht_for_chitip_atoms_[ ii ].set_xaxis_rotation_rad( -1 * res.icoor( chiat4 ).theta() );
		ht_for_chitip_atoms_[ ii ].walk_along_z( res.icoor( chiat4 ).d() );

		HTReal chi_tip_frame(
			res.xyz( chiat2 ),
			res.xyz( chiat3 ),
			res.xyz( chiat4 ) );

		Size const n_nontip_ats_for_chi = res.atoms_last_controlled_by_chi( ii ).size() - 1;

		nonchitip_atoms_[ ii ].reserve( n_nontip_ats_for_chi );
		points_for_nonchitip_atoms_[ ii ].reserve( n_nontip_ats_for_chi );

		for ( Size jj = 1; jj <= res.atoms_last_controlled_by_chi( ii ).size(); ++jj ) {
			Size const jjatom = res.atoms_last_controlled_by_chi( ii )[ jj ];
			if ( jjatom == chiat4 ) continue;

			Vector jjloc_in_chitip_frame = chi_tip_frame.to_local_coordinate( res.xyz( jjatom ) );
			nonchitip_atoms_[ ii ].push_back( jjatom );
			points_for_nonchitip_atoms_[ ii ].push_back( jjloc_in_chitip_frame );
			which_point_for_atom_[ jjatom ] = points_for_nonchitip_atoms_[ ii ].size();
		}
	}
}*/


BuildSet::BuildSet() :
	backbone_only_( false ),
	rot_prob_accumulation_limit_( 0.98 ),
	check_fa_dun_(false),
	fa_dun_cutoff_(0.0)
{}

BuildSet::~BuildSet() {}

BuildSet::BuildSet( BuildSet const & other ) : parent() {
	(*this) = other;
}


BuildSet const &
BuildSet::operator = ( BuildSet const & rhs )
{
	if ( this != &rhs ) {

		restype_ = rhs.restype_;
		backbone_only_ = rhs.backbone_only_;

		restype_geom_ = rhs.restype_geom_;
		rot_prob_accumulation_limit_ = rhs.rot_prob_accumulation_limit_;
		sample_strategy_for_chi_ = rhs.sample_strategy_for_chi_;
		algorithm_ = rhs.algorithm_; // SHALLOW COPY

		nbonds_from_bb_atom_ = rhs.nbonds_from_bb_atom_;
		atom_radii_ = rhs.atom_radii_;

		check_fa_dun_ = rhs.check_fa_dun_;
		fa_dun_cutoff_ = rhs.fa_dun_cutoff_;
	}
	return *this;
}


void BuildSet::set_residue_type(
	core::chemical::ResidueTypeCOP restype,
	bool backbone_only
)
{
	restype_ = restype;
	backbone_only_ = backbone_only;

	restype_geom_ = UpstreamResTypeGeometryOP( new UpstreamResTypeGeometry( *restype_ ) );

	sample_strategy_for_chi_.clear();
	sample_strategy_for_chi_.resize( restype->nchi() );

	nbonds_from_bb_atom_.resize( restype_->natoms() );
	atom_radii_.resize( restype_->natoms() );

	for ( Size ii = 1; ii <= nbonds_from_bb_atom_.size(); ++ii ) {
		if ( restype_->atom_is_backbone( ii )) {
			nbonds_from_bb_atom_[ ii ] = 0;
		} else {
			Size min_path_dist = 6; /// assume that we only care about bond distances less than 6; anything greater than 5 is "infinite"
			for ( Size jj = 1; jj <= nbonds_from_bb_atom_.size(); ++jj ) {
				if ( restype_->atom_is_backbone( jj ) && (Size) restype_->path_distance( ii, jj ) < min_path_dist ) {
					min_path_dist = restype_->path_distance( ii, jj );
				}
			}
			nbonds_from_bb_atom_[ ii ] = min_path_dist;
		}
	}

	for ( Size ii = 1; ii <= atom_radii_.size(); ++ii ) {
		atom_radii_[ ii ] = probe_radius_for_atom_type( restype_->atom( ii ).atom_type_index() );
	}
}

void BuildSet::set_sample_strategy_for_chi(
	Size chi,
	SampleStrategyData const & data
)
{
	sample_strategy_for_chi_[ chi ] = data;
}


void BuildSet::set_downstream_algorithm(
	downstream::DownstreamAlgorithmOP algorithm
)
{
	algorithm_ = algorithm;
}

void
BuildSet::set_fa_dun_cutoff(
	core::Real cutoff )
{
	fa_dun_cutoff_ = cutoff;
	if( cutoff == 0.0 ) check_fa_dun_ = false;
	else check_fa_dun_ = true;
}


FullChiSampleSet::FullChiSampleSet(
	BuildSet const & build_set,
	core::pack::dunbrack::DunbrackRotamerSampleData const & sample,
	bool dry_run
) :
	dry_run_( dry_run ),
	num_chi_samples_total_( 1 )
{
	Size nchi = build_set.restype().nchi();
	if ( ! dry_run ) {
		if ( nchi == 0 ) {
			n_samples_per_chi_.resize( 1, 1 ); /// pretend we have one sample and one chi
		} else {
			n_samples_per_chi_.resize( nchi );
			std::fill( n_samples_per_chi_.begin(), n_samples_per_chi_.end(), 0 );
			chi_samples_.resize( nchi );
			frames_.resize( nchi );
		}
	}

	for ( Size ii = 1; ii <= nchi; ++ii ) {

		if ( ii > sample.nchi() ) {
			expand_non_dunbrack_chi( ii, build_set );
		} else {
			switch ( build_set.sample_strategy_for_chi( ii ).strategy() ) {
				case no_samples:
					if ( !dry_run ) {
						n_samples_per_chi_[ ii ] = 1;
						chi_samples_[ ii ].resize( 1 );
						frames_[ ii ].resize( 1 );
						chi_samples_[ ii ][ 1 ] = numeric::dihedral_degrees(
							build_set.restype().atom( build_set.restype().chi_atoms( ii )[ 1 ] ).ideal_xyz(),
							build_set.restype().atom( build_set.restype().chi_atoms( ii )[ 2 ] ).ideal_xyz(),
							build_set.restype().atom( build_set.restype().chi_atoms( ii )[ 3 ] ).ideal_xyz(),
							build_set.restype().atom( build_set.restype().chi_atoms( ii )[ 4 ] ).ideal_xyz() );
					}
					break;
				case follow_EX_flags :
					expand_samples_by_ex_behavior( ii, ex_level_from_flags( ii ), sample );
				break;

				case rotameric_chi_mimic_EX_flags :
					expand_samples_by_ex_behavior( ii, build_set.sample_strategy_for_chi( ii ).sample_level(), sample );
				break;

				case rotameric_chi_step_wi_sd_range :
					expand_samples_by_steps_wi_sdrange( ii, build_set.sample_strategy_for_chi( ii ), sample );
				break;

				/// distinction without a difference; when would you ever want to sample the boundary point twice!?
				case nonrotameric_chi_sample_wi_nrchi_bin :
				case nonrotameric_chi_sample_wi_nrchi_bin_to_lower_boundary :
					expand_samples_for_nrchi_wi_nrchi_bin( ii, build_set.sample_strategy_for_chi( ii ), sample );
				break;

				/// UNIMPLEMENTED BELOW:
				case rotameric_chi_step_by_value :
				//break;
				case rotameric_chi_partition_sd_range :
				//break;
				//break;
				utility_exit_with_message( "unimplemented rotamer building strategy" );
				break;
			}
		}
		if ( ! dry_run_ ) {
			create_hts_for_chi( ii );
		}
	}

}

void
FullChiSampleSet::expand_non_dunbrack_chi( Size chi, BuildSet const & build_set )
{
	using namespace core::pack::task;

	if ( build_set.restype().is_proton_chi( chi )  ) {
		Size const proton_chi = build_set.restype().chi_2_proton_chi( chi );
		bool extra_proton_chi( ex_level_from_flags( chi ) > NO_EXTRA_CHI_SAMPLES );

		Size nsamples( 1 );

		if ( build_set.sample_strategy_for_chi( chi ).strategy() == no_samples ) {
			/// noop
			nsamples = 1;
		} else if ( extra_proton_chi ) {
			nsamples =  build_set.restype().proton_chi_samples( proton_chi ).size() *
				( 1 + 2*build_set.restype().proton_chi_extra_samples( proton_chi ).size() );
		} else {
			nsamples = build_set.restype().proton_chi_samples( proton_chi ).size();
		}

		if ( ! dry_run_ ) {
			chi_samples_[ chi ].reserve( nsamples );
			if ( build_set.sample_strategy_for_chi( chi ).strategy() == no_samples ) {
				chi_samples_[ chi ].push_back( build_set.restype().proton_chi_samples( proton_chi )[ 1 ] );
			} else {
				for ( Size ii = 1; ii <= build_set.restype().proton_chi_samples( proton_chi ).size(); ++ii ) {
					if ( extra_proton_chi ) {
						for ( Size jj = 1; jj <= build_set.restype().proton_chi_extra_samples( proton_chi ).size(); ++jj ) {
							chi_samples_[ chi ].push_back( build_set.restype().proton_chi_samples( proton_chi )[ ii ] -
								build_set.restype().proton_chi_extra_samples( proton_chi )[ jj ] );
						}
					}
					chi_samples_[ chi ].push_back( build_set.restype().proton_chi_samples( proton_chi )[ ii ] );
					if ( extra_proton_chi ) {
						for ( Size jj = 1; jj <= build_set.restype().proton_chi_extra_samples( proton_chi ).size(); ++jj ) {
							chi_samples_[ chi ].push_back( build_set.restype().proton_chi_samples( proton_chi )[ ii ] +
								build_set.restype().proton_chi_extra_samples( proton_chi )[ jj ] );
						}
					}
				}
			}
			n_samples_per_chi_[ chi ] = nsamples;
		}
		num_chi_samples_total_ *= nsamples;
	} else {
		utility_exit_with_message( "Unrecognized chi type in ProteinUpstreamBuilder for " +
			build_set.restype().name() + " chi# " + utility::to_string( chi ) );
	}
}

void
FullChiSampleSet::expand_samples_by_ex_behavior(
	Size chi,
	ExtraRotSample behavior,
	core::pack::dunbrack::DunbrackRotamerSampleData const & sample
)
{
	using namespace core::pack::task;

	Size nsamps( 1 );

	if ( sample.chi_is_nonrotameric( chi ) ) {
		if ( behavior != NO_EXTRA_CHI_SAMPLES ) {
			nsamps = 3;
			if ( ! dry_run_ ) {
				n_samples_per_chi_[ chi ] = nsamps;
				chi_samples_[ chi ].reserve( nsamps );
				chi_samples_[ chi ].push_back( sample.nrchi_lower_boundary() + 0.5 * ( sample.chi_mean()[ chi ] - sample.nrchi_lower_boundary() ));
				chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] );
				chi_samples_[ chi ].push_back( sample.nrchi_upper_boundary() + 0.5 * ( sample.nrchi_upper_boundary() - sample.chi_mean()[ chi ] ));
			}
		} else {
			nsamps = 1;
			if ( ! dry_run_ ) {
				n_samples_per_chi_[ chi ] = nsamps;
				chi_samples_[ chi ].reserve( nsamps );
				chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] );
			}
		}
	} else { /// rotameric chi

		switch ( behavior ) {
			case NO_EXTRA_CHI_SAMPLES :
				if ( ! dry_run_ ) {
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] );
					n_samples_per_chi_[ chi ] = 1;
				}
			break;
			case EX_ONE_STDDEV :
				nsamps = 3;
				if ( ! dry_run_ ) {
					chi_samples_[ chi ].reserve( nsamps );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] - sample.chi_sd()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] + sample.chi_sd()[ chi ] );
					n_samples_per_chi_[ chi ] = nsamps;
				}
			break;
			case EX_ONE_HALF_STEP_STDDEV :
				nsamps = 3;
				if ( ! dry_run_ ) {
					chi_samples_[ chi ].reserve( nsamps );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] - 0.5 * sample.chi_sd()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] + 0.5 * sample.chi_sd()[ chi ] );
					n_samples_per_chi_[ chi ] = nsamps;
				}
			break;
			case EX_TWO_FULL_STEP_STDDEVS :
				nsamps = 5;
				if ( ! dry_run_ ) {
					chi_samples_[ chi ].reserve( nsamps );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] - 2 * sample.chi_sd()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] - 1 * sample.chi_sd()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] + 1 * sample.chi_sd()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] + 2 * sample.chi_sd()[ chi ] );
					n_samples_per_chi_[ chi ] = nsamps;
				}
			break;
			case EX_TWO_HALF_STEP_STDDEVS :
				nsamps = 5;
				if ( ! dry_run_ ) {
					chi_samples_[ chi ].reserve( nsamps );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] - 1.0 * sample.chi_sd()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] - 0.5 * sample.chi_sd()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] + 0.5 * sample.chi_sd()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] + 1.0 * sample.chi_sd()[ chi ] );
					n_samples_per_chi_[ chi ] = nsamps;
				}
			break;
			case EX_FOUR_HALF_STEP_STDDEVS :
				nsamps = 9;
				if ( ! dry_run_ ) {
					chi_samples_[ chi ].reserve( nsamps );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] - 2.0 * sample.chi_sd()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] - 1.5 * sample.chi_sd()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] - 1.0 * sample.chi_sd()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] - 0.5 * sample.chi_sd()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] + 0.5 * sample.chi_sd()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] + 1.0 * sample.chi_sd()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] + 1.5 * sample.chi_sd()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] + 2.0 * sample.chi_sd()[ chi ] );
					n_samples_per_chi_[ chi ] = nsamps;
				}
			break;
			case EX_THREE_THIRD_STEP_STDDEVS :
				nsamps = 7;
				if ( ! dry_run_ ) {
					chi_samples_[ chi ].reserve( nsamps );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] - 1.0 * sample.chi_sd()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] - 0.67 * sample.chi_sd()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] - 0.33 * sample.chi_sd()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] + 0.33 * sample.chi_sd()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] + 0.67 * sample.chi_sd()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] + 1.0 * sample.chi_sd()[ chi ] );
					n_samples_per_chi_[ chi ] = nsamps;
				}
			break;
			case EX_SIX_QUARTER_STEP_STDDEVS :
				nsamps = 13;
				if ( ! dry_run_ ) {
					chi_samples_[ chi ].reserve( nsamps );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] - 1.50 * sample.chi_sd()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] - 1.25 * sample.chi_sd()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] - 1.00 * sample.chi_sd()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] - 0.75 * sample.chi_sd()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] - 0.50 * sample.chi_sd()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] - 0.25 * sample.chi_sd()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] + 0.25 * sample.chi_sd()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] + 0.50 * sample.chi_sd()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] + 0.75 * sample.chi_sd()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] + 1.00 * sample.chi_sd()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] + 1.25 * sample.chi_sd()[ chi ] );
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] + 1.50 * sample.chi_sd()[ chi ] );
					n_samples_per_chi_[ chi ] = nsamps;
				}
			break;
			default :
				utility_exit_with_message( "Requested ex sample level not yet supported" );
			break;
		}
	}
	num_chi_samples_total_ *= nsamps;
}

void
FullChiSampleSet::expand_samples_by_steps_wi_sdrange(
	Size chi,
	SampleStrategyData const & stratdat,
	core::pack::dunbrack::DunbrackRotamerSampleData const & sample
)
{
	Size nsamps( 1 );

	runtime_assert( stratdat.strategy() == rotameric_chi_step_wi_sd_range );
	runtime_assert( stratdat.sd_range()  > 0 );
	runtime_assert( stratdat.step_size() > 0 );
	Size nsteps = static_cast< int > ( (sample.chi_sd()[ chi ] *  stratdat.sd_range()) / stratdat.step_size() );
	nsamps = 1 + 2 * nsteps;

	if ( ! dry_run_ ) {
		n_samples_per_chi_[ chi ] = nsamps;
		chi_samples_[ chi ].reserve( nsamps );
		for ( Size ii = nsteps; ii >= 1; --ii ) {
			chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] - ii * stratdat.step_size() );
		}
		chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] );
		for ( Size ii = 1; ii <= nsteps; ++ii ) {
			chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] + ii * stratdat.step_size() );
		}
	}

	num_chi_samples_total_ *= nsamps;
}

void
FullChiSampleSet::expand_samples_for_nrchi_wi_nrchi_bin(
	Size chi,
	SampleStrategyData const & stratdat,
	core::pack::dunbrack::DunbrackRotamerSampleData const & sample
)
{
	Size nsamps( 1 );
	if ( sample.chi_is_nonrotameric( chi ) ) {
		if ( sample.probability() > stratdat.nrchi_prob_minimum_for_extra_samples() ) {
			/// expand!
			nsamps = 1 + ( stratdat.n_samples_per_side_of_nrchi_bin() == 0 ?
				0 : 2 * stratdat.n_samples_per_side_of_nrchi_bin() - 1 );
			if ( ! dry_run_ ) {
				n_samples_per_chi_[ chi ] = nsamps;
				chi_samples_[ chi ].reserve( nsamps );
				{// lower
				Real lower_edge = sample.nrchi_lower_boundary();
				Real step = ( sample.chi_mean()[ chi ] - lower_edge ) / stratdat.n_samples_per_side_of_nrchi_bin();
				for ( Size ii = 0; ii < stratdat.n_samples_per_side_of_nrchi_bin(); ++ii ) {
					chi_samples_[ chi ].push_back( lower_edge + ii * step  );
				}
				}
				chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] );
				{// upper
				Real upper_edge = sample.nrchi_upper_boundary();
				Real step = ( upper_edge - sample.chi_mean()[ chi ] ) / stratdat.n_samples_per_side_of_nrchi_bin();
				for ( Size ii = 1; ii < stratdat.n_samples_per_side_of_nrchi_bin(); ++ii ) {
					chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] + ii * step  );
				}
				}
			}

		} else {
			nsamps = 1;
			if ( ! dry_run_ ) {
				n_samples_per_chi_[ chi ] = nsamps;
				chi_samples_[ chi ].reserve( nsamps );
				chi_samples_[ chi ].push_back( sample.chi_mean()[ chi ] );
			}
		}
	} else {
		utility_exit_with_message( "Cannot treat chi " + utility::to_string( chi ) + " as non-rotameric; did you forget the -dun10 flag?" );
	}
	num_chi_samples_total_ *= nsamps;
}


void
FullChiSampleSet::create_hts_for_chi( Size chi )
{
	assert( ! dry_run_ );
	Size const nsamples = chi_samples_[ chi ].size();
	runtime_assert( nsamples == n_samples_per_chi_[ chi ] );
	frames_[ chi ].resize( nsamples );
	for ( Size ii = 1; ii <= nsamples; ++ii ) {
		frames_[ chi ][ ii ].set_zaxis_rotation_deg( chi_samples_[ chi ][ ii ] );
	}
}


FullChiSampleSet::ExtraRotSample
FullChiSampleSet::ex_level_from_flags( Size chi )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys::packing;
	using namespace core::pack::task;

	ExtraRotSample sample_level( NO_EXTRA_CHI_SAMPLES );

	switch ( chi ) {
		case 1:
			if ( option[ ex1::ex1 ] ) sample_level = EX_ONE_STDDEV;
			if ( option[ ex1::level ].user() && option[ ex1::level ] > sample_level ) sample_level = ExtraRotSample( option[ ex1::level ]() );
		break;
		case 2 :
			if ( option[ ex2::ex2 ] ) sample_level = EX_ONE_STDDEV;
			if ( option[ ex2::level ].user() && option[ ex2::level ] > sample_level ) sample_level = ExtraRotSample( option[ ex2::level ]() );
		break;
		case 3 :
			if ( option[ ex3::ex3 ] ) sample_level = EX_ONE_STDDEV;
			if ( option[ ex3::level ].user() && option[ ex3::level ] > sample_level ) sample_level = ExtraRotSample( option[ ex3::level ]() );
		break;
		case 4 :
			if ( option[ ex4::ex4 ] ) sample_level = EX_ONE_STDDEV;
			if ( option[ ex4::level ].user() && option[ ex4::level ] > sample_level ) sample_level = ExtraRotSample( option[ ex4::level ]() );
		break;
	}
	return sample_level;
}

//Kui 112509 Native. Initial avoid_building_any_rotamers_dueto_native_ to false
ProteinUpstreamBuilder::ProteinUpstreamBuilder() : use_input_sc_( false ), avoid_building_any_rotamers_dueto_native_( false )
{}

ProteinUpstreamBuilder::~ProteinUpstreamBuilder() {}

UpstreamBuilderOP
ProteinUpstreamBuilder::clone() const {
	return UpstreamBuilderOP( new ProteinUpstreamBuilder( *this ) );
}

//Kui 110609 Native
void ProteinUpstreamBuilder::set_native_flag (bool native){
	avoid_building_any_rotamers_dueto_native_ = native;
}


/// @details
/// This is the main workhorse function for the upstream builder.  There are lots of ways
/// to customize the way the loops in this function work, but custom code should be written
/// in subroutines called by this function and should not change the function itself.
/// It is crucial here that the rotamers are created in the same way they were
/// counted in insert().  Make sure that if you change build()'s behavior, that insert()'s
/// behavior changes, too!
///
/// DO NOT LET THIS FUNCTION EXCEED 150 LINES.  CLEAR CODE IS SHORT.
std::list< Hit >
ProteinUpstreamBuilder::build(
	ScaffoldBuildPoint const & build_point
) const
{
	// std::cout << "APL DEBUG ProteinUpstreamBuilder::build begin" << std::endl;
	/// Algorithm in brief:
	/// declare a local hit list; this list will be returned at the end of this method.
	/// Iterate across residue types to build
	///    Get the list of base chi samples.
	///    Build CB coordinate frame.  If CB is outside the found-hits accessibility volume, continue.
	///    for each base chi sample (that's sufficiently probable), enumerate the expanded chi samples,
	///    and precompute the HT's for each chi sample
	///       Enumerate the expanded chi combinations (with a LexicographicalIterator )
	///         If, when building atoms for chi i, an atom collides with the background, continue.
	///         If the chitip atom for chi i is outside of the found-hits accessibility volume, continue.
	///         Fill the ResidueCoordinates object with the coordinates for this chi
	///         Iterate across the ExternalGeomSamplers, calling parent::build_from_three_coords() method
	///            splice the returned hit list of returned hits into the local hit list.
	/// return the local list of hits

	using namespace utility;

	std::list< Hit > local_hit_list;

	Size upstream_state = 1;
	Size n_possible_hits = 0;
	for ( Size ii = 1; ii <= build_sets_.size(); ++ii ) {
		core::conformation::Residue rescoords( build_sets_[ ii ].restype(), false );
		//TR << "Building residue type: " << rescoords.type().name() << std::endl;

		core::pack::dunbrack::RotamerLibraryScratchSpace dunscratch; //in case we're checking fa_dun for rotamer
		// Warning: rotlib is assumed to be always valid, yet it throws a bad_weak_ptr exception in the match integration test -- what should it be?
		// Do we expect the rotamer library to be always there? Was:
		// core::pack::dunbrack::SingleResidueRotamerLibraryCOP rotlib(  core::pack::dunbrack::RotamerLibrary::get_instance().get_rsd_library( rescoords.type() ) );
		core::pack::dunbrack::SingleResidueRotamerLibraryCOP rotlib = core::pack::dunbrack::RotamerLibrary::get_instance().get_rsd_library( rescoords.type() ).lock();
		bool check_fa_dun( build_sets_[ii].check_fa_dun() );
		core::Real fa_dun_cutoff( build_sets_[ii].fa_dun_cutoff() );
		if( check_fa_dun ){
				ProteinBackboneBuildPoint const & bb(
			static_cast< ProteinBackboneBuildPoint const & >
			( build_point ));
				rescoords.mainchain_torsions()[1] = bb.phi();
				rescoords.mainchain_torsions()[2] = bb.psi();
		}

		Size const ii_nchi = build_sets_[ ii ].restype().nchi(); // may be different from the number of chi that the Dunbrack library defines
		UpstreamResTypeGeometry const & geom( build_sets_[ ii ].restype_geometry() );

		vector0< HTReal > chitip_frames( ii_nchi + 1 );
		Real accumulated_probability( 0.0 );

		/// Query the rotamer sampler for the list of samples to take
		DunbrackRotamerSampleDataVector rotamer_samples = sampler_->samples( build_point, build_sets_[ ii ].restype() );

		/// Load the Residue with the build point coordinates, and compute the CBeta frame.
		/// Frame 0 defines the frame at CBeta -- to build the chi-tip atom for chi i,
		/// build from the frame for i-1.  For chi 1, start with the frame for "chi 0"
		chitip_frames[ 0 ] = initialize_rescoords( ii, rescoords, build_point );
		bool const avoid_building_any_rotamers = atom_coordinate_unacceptable(
			ii, rescoords, build_sets_[ ii ].restype_geometry().CB_atom_id() );


		Size ii_nrotamers( 0 );

		for ( Size jj = 1; jj <= rotamer_samples.size(); ++jj ) {

			/// Compute the set of extra samples to take for this rotamer.
			/// The rules for building extra samples are not known to the upstream builder,
			/// so if the upstream builder is to maintain an accurate "state" id for
			/// later reconstructing the rotamer described by a hit, the upstream builder
			/// cannot skip the logic for counting how many extra rotamers there are, and how
			/// many extra states are being skipped.
			FullChiSampleSet additional_chi_samples(
				build_sets_[ ii ], rotamer_samples[ jj ], avoid_building_any_rotamers );

			//Kui 110509 Native codes
			// avoid_building_any_rotamers_dueto_native_ is a native residue flag.
			// If user specific this option, he/she want to use native residue only.
			if ( avoid_building_any_rotamers || avoid_building_any_rotamers_dueto_native_ ) {
				/// With a count of how many extra states there are for this rotamer, the upstream
				/// state may safely advance to the next rotamer.
				upstream_state += additional_chi_samples.num_chi_samples_total();
			} else {
				utility::LexicographicalIterator lex( additional_chi_samples.n_samples_per_chi() );

				/// Iterate across all combinations of chi dihedrals for this rotamer.

				Size n_chi_needing_update = ii_nchi;
				while ( ! lex.at_end() ) {
					assert( ( n_chi_needing_update >= 1 || ii_nchi == 0 ) && n_chi_needing_update <= ii_nchi );
					/// Some rotamer placements may be skipped because they either collide with the
					/// background, or they would produce downstream geometries that are too distant
					/// from any of the hits generated in previous stages (and therefore would be
					/// unable to produce a hit that matched the other hits).
					bool rotamer_acceptable( true );

					/// Iterate outwards from the closest chi to the most distal chi, starting
					/// from the chi which last changed: all coordinate frames before that last chi
					/// are still valid since (by construction) the chi upstream has not changed from the
					/// previous iteration.
					for ( Size kk = ii_nchi - n_chi_needing_update + 1; kk <= ii_nchi; ++kk ) {
						chitip_frames[ kk ] = chitip_frames[ kk - 1 ] *
							additional_chi_samples.frame( kk, lex[ kk ] ) *
							geom.ht_for_chitip_atom( kk ); // matrix * matrix * matrix
						rescoords.set_xyz( geom.chitip_atom( kk ), chitip_frames[ kk ].point());
						rescoords.chi()[ kk ] = additional_chi_samples.chi_sample( kk, lex[ kk ] );
						if ( atom_coordinate_unacceptable( ii, rescoords, geom.chitip_atom( kk ) ) ) {
							rotamer_acceptable = false;
						}
						if ( rotamer_acceptable ) {
							// if the chitip atom is unacceptible, then don't bother computing the coordinates
							// for the other atoms controlled by this chi.
							for ( Size ll = 1; ll <= geom.n_nonchitip_atoms_for_chi( kk ); ++ll ) {
								Size const ll_atno = geom.nonchitip_atom( kk, ll );
								// matrix * vector
								rescoords.set_xyz( ll_atno, chitip_frames[ kk ] * geom.points_for_nonchitip_atoms( kk )[ ll ]);
								if ( atom_coordinate_unacceptable( ii, rescoords, ll_atno )) {
									rotamer_acceptable = false;
									break;
								}
							}
						}

						//dunbrack energy check
						//apparently we have to specifically
						//set each chi in the residue object
						//this check can only be performed after we have specified all the chi angles, so it should only
						//be performed in the very last iteration through this loop.
						if ( kk == ii_nchi && rotamer_acceptable && check_fa_dun ) {
							for ( core::Size res_chi(1); res_chi<= ii_nchi; ++res_chi) {
								rescoords.chi()[res_chi] = additional_chi_samples.chi_sample( res_chi, lex[ res_chi ] );
							}
							if( rotlib->rotamer_energy( rescoords, dunscratch ) >= fa_dun_cutoff ){
								rotamer_acceptable = false;
							}
						}

						if ( ! rotamer_acceptable ) {
							/// Bail out.  We have a collision for an atom whose coordinates depend on chi kk.
							/// or we've determined that chi kk cannot place the downstream partner sufficiently close
							/// to any of the hits generated in previous rounds.
							/// Advance the lexicographical iterator for dimension kk.  If this is the last chi sample
							/// for chi kk, then chi kk - 1 will be advanced (or kk - 2 if kk - 1 is at it's last
							/// chi sample, etc.).  The coordinates for everything downstream of
							/// ii_nchi - n_chi_needing_update + 1 will need updating in the next round (if there is one).
							Size n_skipped = 1;
							for ( Size ll = kk+1; ll <= ii_nchi; ++ll ) {
								n_skipped *= lex.dimsize( ll );
							}
							upstream_state += n_skipped;
							ii_nrotamers += n_skipped;
							n_chi_needing_update = lex.continue_at_dimension( kk );
							break;
						}
					} // end update coordinates for residue
					if ( ! rotamer_acceptable ) continue; // lex iter has already been updated, and upstream_state already advanced.

					//std::cout << "ScaffID: " << build_point.index() << " built rotamer # " << upstream_state << " " << rescoords.name();
					//for ( Size chi = 1; chi <= ii_nchi; ++chi ) std::cout << " " << additional_chi_samples.chi_sample( chi, lex[ chi ] ); std::cout << std::endl;

					/// Excellent!  We've arrived at a conformation for this rotamer that's acceptible (collision free).
					/// Hand off work to the downstream hit construction algorithm.
					std::list< Hit > hits = build_sets_[ ii ].algorithm().build(
						build_point.index(),
						upstream_state,
						rescoords );
					local_hit_list.splice( local_hit_list.end(), hits );

					/// proceed to the next chi sample -- track how many chi changed and build
					/// from the first chi that did not change -- that chi's coordinate frame is still
					/// up to date.
					n_chi_needing_update = ++lex;
					++upstream_state;
					++ii_nrotamers;
				}

			}

			/// build rotamers that are sufficiently probable and bail out once we reach
			/// the dregs.  The user controls the dreg limit.  Note, however, that matching
			/// does nothing to point out that low probability rotamers are poor!  Include
			/// low probability rotamers at your own risk.
			accumulated_probability += rotamer_samples[ jj ].probability();
			if ( accumulated_probability > build_sets_[ ii ].probability_limit() ) break;
		} //jj loop over rotamer_samples.size()

		/// Finally, build rotamers for the input side chain, if the original scaffold build point
		/// is of the appropriate type
		/// Kui 110609 Native
		/// avoid_building_any_rotamers_dueto_native_ specific usage of native residue only.
		if ( use_input_sc_ || avoid_building_any_rotamers_dueto_native_) {

			OriginalBackboneBuildPoint const * orig =
				dynamic_cast< OriginalBackboneBuildPoint const * > ( & build_point );
			if ( orig && & (orig->input_conformation().type()) == & (build_sets_[ ii ].restype()) ) {
				/// restype match
				std::list< Hit > hits = build_sets_[ ii ].algorithm().build(
					build_point.index(),
					upstream_state,
					orig->input_conformation() );
				local_hit_list.splice( local_hit_list.end(), hits );
				++upstream_state;
				n_possible_hits += build_sets_[ ii ].algorithm().n_possible_hits_per_upstream_conformation();
				//Kui 110609 Native debug
				//TR << "use_input_sc_:" << use_input_sc_ << " avoid_building_any_rotamers_dueto_native_:"
				//	<< avoid_building_any_rotamers_dueto_native_ << std::endl;
			}
		}

		n_possible_hits += ii_nrotamers * build_sets_[ ii ].algorithm().n_possible_hits_per_upstream_conformation();
		//TR << "Finished building residue type: " << rescoords.type().name() << std::endl;
	} //ii loop over build sets

	assert( dynamic_cast< ProteinBackboneBuildPoint const * > ( & build_point ) );
	ProteinBackboneBuildPoint const & bb( static_cast< ProteinBackboneBuildPoint const & >
		( build_point ) );

	TR << "Considered " << n_possible_hits << " downstream conformations at residue " << bb.original_insertion_point()  << " and found " << local_hit_list.size() << " hits." << std::endl;
	//std::cout << "APL DEBUG ProteinUpstreamBuilder::build end" << std::endl;
	return local_hit_list;
}



/// @brief Replace residue for the amino acid corresponding to the rotamer indicated in the
/// hit at the Pose position seqpos_to_insert_at.
/// @details It is crucial here that the rotamers are counted in the same way they were
/// iterated across in build.  Make sure that if build() changes its behavior, that this
/// function also changes its behavior!
void
ProteinUpstreamBuilder::recover_hit(
	Hit const & hit,
	ScaffoldBuildPoint const & build_point,
	UpstreamResidueProcessor & processor
) const
{
	std::list< Hit > hitlist;
	hitlist.push_back( hit );
	recover_hits( hitlist.begin(), hitlist.end(), build_point, processor );
}

void
ProteinUpstreamBuilder::recover_hits(
	std::list< Hit >::const_iterator hit_iter,
	std::list< Hit >::const_iterator hits_end,
	ScaffoldBuildPoint const & build_point,
	UpstreamResidueProcessor & processor
) const
{
//	assert( hit_iter == hits_end || hit_iter->scaffold_build_id() == build_point.id() );
  assert( hit_iter == hits_end || hit_iter->scaffold_build_id() == build_point.index() );
	//Hit hit = *hit_iter;
	//std::cout << "ProteinUpstreamBuilder::recover hit " << hit.first()[ 1 ] << " " << hit.first()[ 2 ] << std::endl;
	/// 1. Figure out which amino acid it is that we're inserting.
	Size upstream_state = 1;
	for ( Size ii = 1; ii <= build_sets_.size(); ++ii ) {
		Size const ii_nchi = build_sets_[ ii ].restype().nchi(); // may be different from the number of chi that the Dunbrack library defines
		Real accumulated_probability( 0.0 );

		/// Query the rotamer sampler for the list of samples to take
		DunbrackRotamerSampleDataVector rotamer_samples = sampler_->samples( build_point, build_sets_[ ii ].restype() );

		for ( Size jj = 1; jj <= rotamer_samples.size(); ++jj ) {

			/// Determine how many extra samples to take for this rotamer using the "dry_run" flag of the FullChiSampleSet
			FullChiSampleSet dry_samples(
				build_sets_[ ii ], rotamer_samples[ jj ], /*dry_run*/ true );

			//std::cout << " jj " << jj << " us_state: " << upstream_state << " nsamps: " << dry_samples.num_chi_samples_total() << std::endl;
			//Kui Native 110809 - need more testing for native 2nd matching
			if ( upstream_state + dry_samples.num_chi_samples_total() <= hit_iter->first()[ 2 ]
					|| avoid_building_any_rotamers_dueto_native_ ) {
				/// With a count of how many extra states there are for this rotamer, the upstream
				/// state may safely advance to the next rotamer.
				upstream_state += dry_samples.num_chi_samples_total();
			} else {
				UpstreamResTypeGeometry const & geom( build_sets_[ ii ].restype_geometry() );

				/// Determine the actual chi to be used.
				FullChiSampleSet additional_chi_samples(
					build_sets_[ ii ], rotamer_samples[ jj ], /*dry_run*/ false );

				runtime_assert( dry_samples.num_chi_samples_total() == additional_chi_samples.num_chi_samples_total() );

				utility::LexicographicalIterator lex( additional_chi_samples.n_samples_per_chi() );

				/// Reuse the additional_chi_samples and lex if we have multiple hits that come from
				/// the same base-rotamer.  The exit conditions for this loop are:
				/// the hit list is exhausted, or, the next hit isn't within the set of subrotamers
				/// for this base rotamer.
				while ( true ) {
					lex.set_position_from_index( hit_iter->first()[ 2 ] - upstream_state + 1 );

					utility::vector0< HTReal > chitip_frames( ii_nchi + 1 );
					core::conformation::Residue rescoords( build_sets_[ ii ].restype(), false );
					chitip_frames[ 0 ] = initialize_rescoords( ii, rescoords, build_point );

					for ( Size kk = 1; kk <= ii_nchi; ++kk ) {
						chitip_frames[ kk ] = chitip_frames[ kk - 1 ] *
							additional_chi_samples.frame( kk, lex[ kk ] ) *
							geom.ht_for_chitip_atom( kk ); // matrix * matrix * matrix
						rescoords.set_xyz( geom.chitip_atom( kk ), chitip_frames[ kk ].point());
						rescoords.chi()[ kk ] = additional_chi_samples.chi_sample( kk, lex[ kk ] );
						for ( Size ll = 1; ll <= geom.n_nonchitip_atoms_for_chi( kk ); ++ll ) {
							Size const ll_atno = geom.nonchitip_atom( kk, ll );
							// matrix * vector
							rescoords.set_xyz( ll_atno, chitip_frames[ kk ] * geom.points_for_nonchitip_atoms( kk )[ ll ]);
						}
					}

					//std::cout << "ScaffID: " << build_point.index() << " recover rotamer # " << hit_iter->first()[ 2 ] << " " << rescoords.name();
					//for ( Size chi = 1; chi <= ii_nchi; ++chi ) std::cout << " " << additional_chi_samples.chi_sample( chi, lex[ chi ] ); std::cout << std::endl;

					/// Hand off the coordinates of the Residue to the UpstreamResidueProcessor
					/// Maybe this writes a kinemage file, or inserts the residue into a Pose.
					processor.process_hit( *hit_iter, rescoords );

					std::list< Hit >::const_iterator last_hit = hit_iter;
					++hit_iter;
					if ( hit_iter == hits_end ) return;
					runtime_assert( hit_iter->scaffold_build_id() == last_hit->scaffold_build_id() );
					runtime_assert( hit_iter->upstream_conf_id()  >= last_hit->upstream_conf_id() ); // sorted order!
					if ( upstream_state + dry_samples.num_chi_samples_total() <= hit_iter->first()[ 2 ] ) {
						upstream_state += dry_samples.num_chi_samples_total();
						break;
					}
				}
			}

			/// build rotamers that are sufficiently probable and bail out once we reach
			/// the dregs.  The user controls the dreg limit.  Note, however, that matching
			/// does nothing to point out that low probability rotamers are poor!  Include
			/// low probability rotamers at your own risk.
			accumulated_probability += rotamer_samples[ jj ].probability();
			if ( accumulated_probability > build_sets_[ ii ].probability_limit() ) break;
		}

		/// Is the hit from the input side chain?
		//Kui Native 110809 - need more testing for native 2nd matching
		if ( use_input_sc_ || avoid_building_any_rotamers_dueto_native_ ) {
			OriginalBackboneBuildPoint const * orig =
				dynamic_cast< OriginalBackboneBuildPoint const * > ( & build_point );
			if ( orig && & (orig->input_conformation().type()) == & (build_sets_[ ii ].restype()) ) {
				/// restype match
				if ( hit_iter->first()[ 2 ] == upstream_state ) {
					processor.process_hit( *hit_iter, orig->input_conformation() );
					++hit_iter;
				}
				if ( hit_iter == hits_end ) return;
				++upstream_state;
			}
		}
	}

	// All rotamers should have been found by now; something is wrong if we reach the end
	// of this function instead of leaving trough one of the exit opportunities above.
	utility_exit_with_message("ProteinUpstreamBuilder was unable to recover rotamer " +
		utility::to_string( hit_iter->first()[ 2 ] ) + " from scaffold build position " +
		utility::to_string( hit_iter->first()[ 1 ] ) );
}

ProteinUpstreamBuilder::Size
ProteinUpstreamBuilder::n_restypes_to_build() const
{
	return build_sets_.size();
}

core::chemical::ResidueTypeCOP
ProteinUpstreamBuilder::restype( Size which_restype ) const
{
	return build_sets_[ which_restype ].restype().get_self_ptr();
}


/// @brief The side chain for two
bool
ProteinUpstreamBuilder::compatible(
	Hit const & my_hit,
	ScaffoldBuildPoint const & build_point_mine,
	UpstreamBuilder const & other,
	Hit const & other_hit,
	ScaffoldBuildPoint const & build_point_other,
	bool
) const
{
	return other.compatible( other_hit, build_point_other, *this, my_hit, build_point_mine, false );

}

bool ProteinUpstreamBuilder::compatible(
	Hit const & my_hit,
	ScaffoldBuildPoint const & build_point_mine,
	ProteinUpstreamBuilder const & other_builder,
	Hit const & other_hit,
	ScaffoldBuildPoint const & build_point_other,
	bool
) const
{
	/// This sidechain can only be in one position per residue.  Declare these
	/// SC builders incompatible if they're building from the same build point.
	if ( build_point_mine.index() != build_point_other.index() ) {
		return true;
	} else {
		if ( my_hit.first()[ 2 ] == 1 && build_sets_[ 1 ].backbone_only() ) {
			return true;
		} else if ( other_hit.first()[ 2 ] == 1 && other_builder.build_sets_[ 1 ].backbone_only() ) {
			return true;
		} else {
			return false;
		}
	}
}

void
ProteinUpstreamBuilder::add_build_set(
	BuildSet const & build_set
)
{
	if ( build_set.backbone_only() ) {
		//build_sets_.push_front( build_set );
		runtime_assert( build_sets_.size() == 0 || ! build_sets_[ 1 ].backbone_only() ); /// only one backbone-only build set allowed
		build_sets_.resize( build_sets_.size() + 1 );
		for ( Size ii = build_sets_.size(); ii > 1; --ii ) {
			build_sets_[ ii ] = build_sets_[ ii - 1 ];
		}
		build_sets_[ 1 ] = build_set;
	} else {
		build_sets_.push_back( build_set );
	}
}

BuildSet const &
ProteinUpstreamBuilder::build_set( core::chemical::ResidueTypeCOP restype ) const
{
	return const_cast< ProteinUpstreamBuilder * > (this)->build_set( restype );
}

BuildSet &
ProteinUpstreamBuilder::build_set( core::chemical::ResidueTypeCOP restype )
{
	for ( Size ii = 1; ii <= build_sets_.size(); ++ii ) {
		if ( &build_sets_[ ii ].restype() == &(*restype) ) {
			return build_sets_[ ii ];
		}
	}
	utility_exit_with_message( "Could not retrieve a build set for residue type " + restype->name() );

	// appease the compiler
	return b;
}


void
ProteinUpstreamBuilder::set_sampler(
	ProteinSCSamplerCOP sampler
)
{
	sampler_ = sampler;
}

void ProteinUpstreamBuilder::set_use_input_sidechain( bool setting )
{
	//std::cout << "Setting use_input_sidechain" << setting << std::endl;
	use_input_sc_ = setting;
}



/// @brief Copy the coordinates from the build_point object into the the rescoords
/// object, compute the coordinate frame at CBeta, copy the CB coordinate into the
/// rescoords object, and return the CBeta frame.
ProteinUpstreamBuilder::HTReal
ProteinUpstreamBuilder::initialize_rescoords(
	Size build_set_id,
	core::conformation::Residue & rescoords,
	ScaffoldBuildPoint const & build_point
) const
{
	{ /// scope to test that we have a protein backbone residue

	if ( ! dynamic_cast< ProteinBackboneBuildPoint const * > ( & build_point ) ) {
		utility_exit_with_message( "Input to ProteinUpstreamBuilder not castable to ProteinBackboneBuildPoint *" );
	}

	}

	ProteinBackboneBuildPoint const & bb( static_cast< ProteinBackboneBuildPoint const & >
		( build_point ) );

	UpstreamResTypeGeometry const & geom( build_sets_[ build_set_id ].restype_geometry() );

	rescoords.seqpos( bb.original_insertion_point() );
	rescoords.set_xyz( geom.N_atom_id(), bb.N_pos());
	rescoords.set_xyz( geom.CA_atom_id(),bb.CA_pos());
	rescoords.set_xyz( geom.C_atom_id(), bb.C_pos());
	rescoords.set_xyz( geom.O_atom_id(), bb.O_pos());

	if ( geom.has_H_atom() )  rescoords.set_xyz( geom.H_atom_id(),  bb.H_pos()  );
	if ( geom.has_HA_atom() ) rescoords.set_xyz( geom.HA_atom_id(), bb.HA_pos() );

	if ( geom.has_CB_atom() ) {
		HTReal cb_frame = compute_cb_frame( build_set_id, bb );
		rescoords.set_xyz( geom.CB_atom_id(), cb_frame.point());
		return cb_frame;
	} else {
		/// GLYCINE
		HTReal input_frame( bb.N_pos(), 0.5 * ( bb.N_pos() + bb.C_pos() ), bb.CA_pos() );
		runtime_assert( build_sets_[ build_set_id ].restype().has( "1HA" ) );
		Size oneHA_index( build_sets_[ build_set_id ].restype().atom_index( "1HA" ));
		rescoords.set_xyz( oneHA_index,
			input_frame * geom.coordinate_for_nonchi_atom_in_ideal_frame( oneHA_index ));
		return input_frame;
	}
}

/// @details Follow Phil's convention for placing CBeta based on halfway point of
/// the ideal N and C positions and the halfway point from the the input N and
/// C positions.  In fact, this is barely different from Phil's code except that
/// it uses HT's instead of Stubs.  By "Phil's code", I mean
/// core/conformation/Residue.cc::orient_onto_position()
ProteinUpstreamBuilder::HTReal
ProteinUpstreamBuilder::compute_cb_frame(
	Size build_set_id,
	ProteinBackboneBuildPoint const & build_point
) const
{
	UpstreamResTypeGeometry const & geom( build_sets_[ build_set_id ].restype_geometry()  );

	Size CB_no = geom.CB_atom_id();
	runtime_assert( CB_no != 0 ); // there's no reason that we should be building the CBeta frame if there's no cbeta.

	Vector halfpoint_input = 0.5 * ( build_point.N_pos() + build_point.C_pos() );
	HTReal input_frame( build_point.N_pos(), halfpoint_input, build_point.CA_pos() );
	Vector CBpos = input_frame * geom.coordinate_for_nonchi_atom_in_ideal_frame( CB_no );

	return HTReal( build_point.N_pos(), build_point.CA_pos(), CBpos );
}

/// @details Collision detection against the background goes here.  "Background"
/// is tricky for atoms near the backbone
bool
ProteinUpstreamBuilder::atom_coordinate_unacceptable(
	Size build_set_id,
	core::conformation::Residue const & rescoords,
	Size atomno
) const
{
	/// do not run bump-grid collision detection on any atom that's within
	/// 4 chemical bonds of CA or you run the danger of wrongfully-detecting
	/// a collision with CA

	if ( atomno == 0 ) return false; // glycine CBeta


	Size const MIN_CHEMBOND_SEP_FROM_CA = 4;
	if ( build_sets_[ build_set_id ].nbonds_from_bb_atom( atomno ) > MIN_CHEMBOND_SEP_FROM_CA ) {
		return bbgrid().occupied( build_sets_[ build_set_id ].atom_radius( atomno ), rescoords.xyz( atomno ));
	} else if ( build_sets_[ build_set_id ].nbonds_from_bb_atom( atomno ) == MIN_CHEMBOND_SEP_FROM_CA ) {
		/// check to see if we're not colliding with CA but we are colliding with something else:
		/// ASSUMPTION.  This atom connects to the backbone at a single point -- CAlpha.
		/// This assumes that there are no residues vaguly similar to proline where it would be possible
		/// to be 4 bonds from N but futher than four bonds from CA:
		///  N -- CA
		/// CD      CB
		///  CD1     CG
		///  CD2   CG2
		///   CD3-CG3  <--- CD3 is 4 bonds from N.  From my above assumption, I would incorrectly dismiss.
		///                 some conformations where CD3 collided with N

		Size capos = build_sets_[ build_set_id ].restype_geometry().CA_atom_id();
		Real const min_d = bbgrid().required_separation_distance(
			build_sets_[ build_set_id ].atom_radius( capos ),
			build_sets_[ build_set_id ].atom_radius( atomno ) );
		Real min_d2 = min_d * min_d;
		if ( rescoords.xyz( atomno ).distance_squared( rescoords.xyz( capos )) > min_d2 ) {
			return bbgrid().occupied( build_sets_[ build_set_id ].atom_radius( atomno ), rescoords.xyz( atomno ));
		}
	}


	return false;

}

}
}
}
