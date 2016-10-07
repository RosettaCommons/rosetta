// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/hbonds/HBondOptions.cc
/// @brief  HBondOptions class, holds the options for the hbond energy function
/// @author Matthew O'Meara

/// @detail
/// To add an additional option for hydrogen bonds do the following:
///
/// In HBondOptions.hh:
/// 1) add it to the default constructor
/// 2) add it to the copy constructor
/// 3) add a getter and a setter
/// 4) add it to operator==
/// 5) add it to the private data
/// 6) add it to HBondOptions::show

// Unit Headers
#include <core/scoring/hbonds/HBondOptions.hh>

// Package Headers
#include <basic/options/option.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh> //pba
#include <basic/options/keys/dna.OptionKeys.gen.hh> //sthyme
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/membrane.OptionKeys.gen.hh> //pba
#include <basic/options/keys/mp.OptionKeys.gen.hh>

// C++ Headers
#include <string>
#include <iostream>

#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/vector1.hh>
#include <utility/options/keys/OptionKeyList.hh>


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace hbonds {

/// @details delegating constructor passing in the global option collection
HBondOptions::HBondOptions( std::string params_db_tag ) :
	HBondOptions( basic::options::option, params_db_tag )
{}

HBondOptions::HBondOptions(
	utility::options::OptionCollection const & options,
	std::string params_db_tag
):
	exclude_DNA_DNA_( true ),
	exclude_intra_res_protein_( true ),
	exclude_intra_res_RNA_( false ),
	put_intra_into_total_( false ), // means that intra will go into 'hbond' (all hbond terms)
	exclude_self_hbonds_( true ), // does not seem to be active --rhiju
	use_hb_env_dep_ ( true ),
	use_hb_env_dep_DNA_( true ),
	smooth_hb_env_dep_( true ),
	bb_donor_acceptor_check_( true ),
	decompose_bb_hb_into_pair_energies_( false ),
	params_database_tag_(params_db_tag),
	use_sp2_chi_penalty_( false ),
	sp2_BAH180_rise_( 0.75 ),
	sp2_outer_width_( 0.357 ),
	measure_sp3acc_BAH_from_hvy_( false ),
	fade_energy_( false ),
	Mbhbond_( false ), //pba
	mphbond_( false ), // membrane framework hbond option
	hbond_energy_shift_( 0.0 ),
	length_dependent_srbb_( false ),
	ldsrbb_low_scale_( 0.5 ),
	ldsrbb_high_scale_( 2.0 ),
	ldsrbb_minlength_( 4 ),
	ldsrbb_maxlength_( 17 ),
	use_hb_env_dep_new_ ( false ),  //fpd new envdep HB
	hb_env_dep_new_low_scale_ ( 0.2 ),  //fpd new envdep HB
	hb_env_dep_new_low_nneigh_ ( 10 ),  //fpd new envdep HB
	hb_env_dep_new_high_nneigh_ ( 20 )  //fpd new envdep HB
{
	initialize_from_options( options );
}

void HBondOptions::initialize_from_options()
{
	initialize_from_options( basic::options::option );
}

void HBondOptions::initialize_from_options( utility::options::OptionCollection const & options )
{
	using namespace basic::options;
	if ( params_database_tag_.size() == 0 ) {
		params_database_tag_ = options[OptionKeys::score::hbond_params];
	}

	if ( options.has(OptionKeys::membrane::Mhbond_depth) &&
			options[OptionKeys::membrane::Mhbond_depth].user() ) {
		Mbhbond_ = options[OptionKeys::membrane::Mhbond_depth];//pba
	}

	if ( options.has( OptionKeys::mp::scoring::hbond ) ) {
		mphbond_ = options[ OptionKeys::mp::scoring::hbond ];
	}

	exclude_DNA_DNA_ = options[OptionKeys::dna::specificity::exclude_dna_dna]; // adding because this parameter should absolutely be false for any structure with DNA in it and it doesn't seem to be read in via the weights file method, so now it's an option - sthyme
	decompose_bb_hb_into_pair_energies_ = options[ OptionKeys::score::hbond_bb_per_residue_energy ];
	use_sp2_chi_penalty_ = options[OptionKeys::corrections::score::hb_sp2_chipen ];
	bb_donor_acceptor_check_ = ! options[ OptionKeys::score::hbond_disable_bbsc_exclusion_rule ];
	sp2_BAH180_rise_ = options[ OptionKeys::corrections::score::hb_sp2_BAH180_rise ];
	put_intra_into_total_ = options[basic::options::OptionKeys::score::put_intra_into_total]();

	if ( options[ OptionKeys::score::length_dep_srbb ].user() ) {
		length_dependent_srbb_ = options[ OptionKeys::score::length_dep_srbb ];
	}
	if ( options[ OptionKeys::score::ldsrbb_low_scale ].user() ) {
		ldsrbb_low_scale_ = options[ OptionKeys::score::ldsrbb_low_scale ];
	}
	if ( options[ OptionKeys::score::ldsrbb_high_scale ].user() ) {
		ldsrbb_high_scale_ = options[ OptionKeys::score::ldsrbb_high_scale ];
	}
	if ( options[ OptionKeys::score::ldsrbb_minlength ].user() ) {
		ldsrbb_minlength_ = options[ OptionKeys::score::ldsrbb_minlength ];
	}
	if ( options[ OptionKeys::score::ldsrbb_maxlength ].user() ) {
		ldsrbb_maxlength_ = options[ OptionKeys::score::ldsrbb_maxlength ];
	}

	if ( options[ OptionKeys::score::hb_env_dep_new ].user() ) {
		use_hb_env_dep_new_ = options[ OptionKeys::score::hb_env_dep_new ]();
	}
	if ( options[ OptionKeys::score::hb_env_dep_new_low_scale ].user() ) {
		hb_env_dep_new_low_scale_ = options[ OptionKeys::score::hb_env_dep_new_low_scale ]();
	}
	if ( options[ OptionKeys::score::hb_env_dep_new_low_nneigh ].user() ) {
		hb_env_dep_new_low_nneigh_ = options[ OptionKeys::score::hb_env_dep_new_low_nneigh ]();
	}
	if ( options[ OptionKeys::score::hb_env_dep_new_high_nneigh ].user() ) {
		hb_env_dep_new_high_nneigh_ = options[ OptionKeys::score::hb_env_dep_new_high_nneigh ]();
	}

	if ( options.has(OptionKeys::corrections::score::hb_sp2_outer_width) ) {
		sp2_outer_width_ = options[ OptionKeys::corrections::score::hb_sp2_outer_width ];
	}

	measure_sp3acc_BAH_from_hvy_ = options[ OptionKeys::corrections::score::hbond_measure_sp3acc_BAH_from_hvy ];
	fade_energy_ = options[ OptionKeys::corrections::score::hb_fade_energy ];

	if ( options.has( OptionKeys::corrections::score::hbond_energy_shift) ) {
		hbond_energy_shift_ = options[ OptionKeys::corrections::score::hbond_energy_shift ];
	}
}

void
HBondOptions::list_options_read( utility::options::OptionKeyList & option_list ) {
	using namespace basic::options;
	option_list
		+ OptionKeys::corrections::score::hb_fade_energy
		+ OptionKeys::corrections::score::hbond_energy_shift
		+ OptionKeys::corrections::score::hbond_measure_sp3acc_BAH_from_hvy
		+ OptionKeys::corrections::score::hb_sp2_BAH180_rise
		+ OptionKeys::corrections::score::hb_sp2_chipen
		+ OptionKeys::corrections::score::hb_sp2_outer_width
		+ OptionKeys::dna::specificity::exclude_dna_dna
		+ OptionKeys::membrane::Mhbond_depth
		+ OptionKeys::mp::scoring::hbond
		+ OptionKeys::score::hb_env_dep_new
		+ OptionKeys::score::hb_env_dep_new_high_nneigh
		+ OptionKeys::score::hb_env_dep_new_low_nneigh
		+ OptionKeys::score::hb_env_dep_new_low_scale
		+ OptionKeys::score::hbond_bb_per_residue_energy
		+ OptionKeys::score::hbond_disable_bbsc_exclusion_rule
		+ OptionKeys::score::hbond_params
		+ OptionKeys::score::ldsrbb_high_scale
		+ OptionKeys::score::ldsrbb_low_scale
		+ OptionKeys::score::ldsrbb_maxlength
		+ OptionKeys::score::ldsrbb_minlength
		+ OptionKeys::score::length_dep_srbb;
}

/// copy constructor
HBondOptions::HBondOptions( HBondOptions const & src ):
	ReferenceCount( src )
{
	*this = src;
}

HBondOptions::~HBondOptions(){}

/// copy operator
HBondOptions &
HBondOptions::operator=( HBondOptions const & src )
{
	exclude_DNA_DNA_ = src.exclude_DNA_DNA_;
	exclude_intra_res_protein_ = src.exclude_intra_res_protein_;
	exclude_intra_res_RNA_ = src.exclude_intra_res_RNA_;
	put_intra_into_total_ = src.put_intra_into_total_;
	exclude_self_hbonds_ = src.exclude_self_hbonds_;
	use_hb_env_dep_ = src.use_hb_env_dep_;
	use_hb_env_dep_DNA_ = src.use_hb_env_dep_DNA_;
	smooth_hb_env_dep_ = src.smooth_hb_env_dep_;
	bb_donor_acceptor_check_ = src.bb_donor_acceptor_check_;
	decompose_bb_hb_into_pair_energies_ = src.decompose_bb_hb_into_pair_energies_;
	params_database_tag_ = src.params_database_tag_;
	use_sp2_chi_penalty_ = src.use_sp2_chi_penalty_;
	sp2_BAH180_rise_ = src.sp2_BAH180_rise_;
	sp2_outer_width_ = src.sp2_outer_width_;
	measure_sp3acc_BAH_from_hvy_ = src.measure_sp3acc_BAH_from_hvy_;
	fade_energy_ = src.fade_energy_;
	Mbhbond_ = src.Mbhbond_; //pba
	mphbond_ = src.mphbond_;
	hbond_energy_shift_ = src.hbond_energy_shift_;
	length_dependent_srbb_ = src.length_dependent_srbb_;
	ldsrbb_low_scale_ = src.ldsrbb_low_scale_;
	ldsrbb_high_scale_ = src.ldsrbb_high_scale_;
	ldsrbb_minlength_ = src.ldsrbb_minlength_;
	ldsrbb_maxlength_ = src.ldsrbb_maxlength_;
	use_hb_env_dep_new_ = src.use_hb_env_dep_new_;
	hb_env_dep_new_low_scale_ = src.hb_env_dep_new_low_scale_;
	hb_env_dep_new_low_nneigh_ = src.hb_env_dep_new_low_nneigh_;
	hb_env_dep_new_high_nneigh_ = src.hb_env_dep_new_high_nneigh_;

	return *this;
}


void
HBondOptions::parse_my_tag(
	utility::tag::TagCOP tag
) {
	// hbond options

	// DEPRECATE
	if ( tag->hasOption( "exclude_DNA_DNA_hbond" ) ) {
		exclude_DNA_DNA( tag->getOption<bool>( "exclude_DNA_DNA_hbond" ) );
	}
	if ( tag->hasOption( "use_hb_env_dep_DNA" ) ) {
		use_hb_env_dep_DNA( tag->getOption<bool>( "use_hb_env_dep_DNA" ) );
	}
	if ( tag->hasOption( "use_hb_env_dep" ) ) {
		use_hb_env_dep( tag->getOption<bool>( "use_hb_env_dep" ) );
	}
	if ( tag->hasOption( "smooth_hb_env_dep" ) ) {
		smooth_hb_env_dep( tag->getOption<bool>( "smooth_hb_env_dep" ) );
	}
	if ( tag->hasOption( "decompose_bb_hb_into_pair_energies" ) ) {
		decompose_bb_hb_into_pair_energies( tag->getOption<bool>( "decompose_bb_hb_into_pair_energies" ) );
	}

	if ( tag->hasOption( "hbonds__exclude_DNA_DNA_hbond" ) ) {
		exclude_DNA_DNA( tag->getOption<bool>( "hbonds__exclude_DNA_DNA_hbond" ) );
	}
	if ( tag->hasOption( "hbonds__use_hb_env_dep_DNA" ) ) {
		use_hb_env_dep_DNA( tag->getOption<bool>( "hbonds__use_hb_env_dep_DNA" ) );
	}
	if ( tag->hasOption( "hbonds__put_intra_into_total" ) ) {
		put_intra_into_total( tag->getOption<bool>( "hbonds__put_intra_into_total" ) );
	}
	if ( tag->hasOption( "hbonds__exclude_self_hbonds" ) ) {
		exclude_self_hbonds( tag->getOption<bool>( "hbonds__exclude_self_hbonds" ) );
	}
	if ( tag->hasOption( "hbonds__use_hb_env_dep" ) ) {
		use_hb_env_dep( tag->getOption<bool>( "hbonds__use_hb_env_dep" ) );
	}
	if ( tag->hasOption( "hbonds__smooth_hb_env_dep" ) ) {
		smooth_hb_env_dep( tag->getOption<bool>( "hbonds__smooth_hb_env_dep" ) );
	}
	if ( tag->hasOption( "hbonds__decompose_bb_hb_into_pair_energies" ) ) {
		decompose_bb_hb_into_pair_energies( tag->getOption<bool>( "hbonds__decompose_bb_hb_into_pair_energies" ) );
	}
	if ( tag->hasOption( "hbonds__exclude_intra_res_protein" ) ) {
		exclude_intra_res_protein( tag->getOption<bool>( "hbonds__exclude_intra_res_protein" ) );
	}
	if ( tag->hasOption( "hbonds__exclude_intra_res_RNA" ) ) {
		exclude_intra_res_RNA( tag->getOption<bool>( "hbonds__exclude_intra_res_RNA" ) );
	}
	//if ( tag->hasOption( "hbonds__put_intra_into_total" ) ) {
	// put_intra_into_total( tag->getOption<bool>( "hbonds__put_intra_into_total" ) );
	//}
	if ( tag->hasOption( "hbonds__bb_donor_acceptor_check" ) ) {
		bb_donor_acceptor_check( tag->getOption<bool>( "hbonds__bb_donor_acceptor_check" ) );
	}
	if ( tag->hasOption( "hbonds__params_database_tag" ) ) {
		params_database_tag( tag->getOption<std::string>( "hbonds__params_database_tag" ) );
	}
	if ( tag->hasOption( "hbonds__use_sp2_chi_penalty" ) ) {
		use_sp2_chi_penalty( tag->getOption<bool>( "hbonds__use_sp2_chi_penalty" ) );
	}
	if ( tag->hasOption( "hbonds__sp2_BAH180_rise" ) ) {
		sp2_BAH180_rise( tag->getOption<Real>( "hbonds__sp2_BAH180_rise" ) );
	}
	if ( tag->hasOption( "hbonds__sp2_outer_width" ) ) {
		sp2_outer_width( tag->getOption<Real>( "hbonds__sp2_outer_width" ) );
	}
	if ( tag->hasOption( "hbonds__measure_sp3acc_BAH_from_hvy" ) ) {
		measure_sp3acc_BAH_from_hvy( tag->getOption<bool>( "hbonds__measure_sp3acc_BAH_from_hvy" ) );
	}
	if ( tag->hasOption( "hbonds__fade_energy" ) ) {
		fade_energy( tag->getOption<bool>( "hbonds__fade_energy" ) );
	}
	if ( tag->hasOption( "hbonds__Mbhbond" ) ) {
		bb_donor_acceptor_check( tag->getOption<bool>( "hbonds__Mbhbond" ) );
	}
	if ( tag->hasOption( "hbonds__mphbond" ) ) {
		bb_donor_acceptor_check( tag->getOption<bool>( "hbonds__mphbond" ) );
	}
}

void
HBondOptions::append_schema_attributes( utility::tag::AttributeList & attributes )
{
	using namespace utility::tag;
	attributes
		+ XMLSchemaAttribute( "exclude_DNA_DNA_hbond", xsct_rosetta_bool )
		+ XMLSchemaAttribute( "use_hb_env_dep_DNA", xsct_rosetta_bool )
		+ XMLSchemaAttribute( "use_hb_env_dep", xsct_rosetta_bool )
		+ XMLSchemaAttribute( "smooth_hb_env_dep", xsct_rosetta_bool )
		+ XMLSchemaAttribute( "decompose_bb_hb_into_pair_energies", xsct_rosetta_bool )
		+ XMLSchemaAttribute( "hbonds__exclude_DNA_DNA_hbond", xsct_rosetta_bool )
		+ XMLSchemaAttribute( "hbonds__use_hb_env_dep_DNA", xsct_rosetta_bool )
		+ XMLSchemaAttribute( "hbonds__put_intra_into_total", xsct_rosetta_bool )
		+ XMLSchemaAttribute( "hbonds__exclude_self_hbonds", xsct_rosetta_bool )
		+ XMLSchemaAttribute( "hbonds__use_hb_env_dep", xsct_rosetta_bool )
		+ XMLSchemaAttribute( "hbonds__smooth_hb_env_dep", xsct_rosetta_bool )
		+ XMLSchemaAttribute( "hbonds__decompose_bb_hb_into_pair_energies", xsct_rosetta_bool )
		+ XMLSchemaAttribute( "hbonds__exclude_intra_res_protein", xsct_rosetta_bool )
		+ XMLSchemaAttribute( "hbonds__exclude_intra_res_RNA", xsct_rosetta_bool )
		//+ XMLSchemaAttribute( "hbonds__put_intra_into_total", xsct_rosetta_bool )
		+ XMLSchemaAttribute( "hbonds__bb_donor_acceptor_check", xsct_rosetta_bool )
		+ XMLSchemaAttribute( "hbonds__params_database_tag", xs_string )
		+ XMLSchemaAttribute( "hbonds__use_sp2_chi_penalty", xsct_rosetta_bool )
		+ XMLSchemaAttribute( "hbonds__sp2_BAH180_rise", xs_decimal )
		+ XMLSchemaAttribute( "hbonds__sp2_outer_width", xs_decimal )
		+ XMLSchemaAttribute( "hbonds__measure_sp3acc_BAH_from_hvy", xsct_rosetta_bool )
		+ XMLSchemaAttribute( "hbonds__fade_energy", xsct_rosetta_bool )
		+ XMLSchemaAttribute( "hbonds__Mbhbond", xsct_rosetta_bool )
		+ XMLSchemaAttribute( "hbonds__mphbond", xsct_rosetta_bool );

}


bool
HBondOptions::exclude_DNA_DNA() const
{
	return exclude_DNA_DNA_;
}


void
HBondOptions::exclude_DNA_DNA( bool const setting )
{
	exclude_DNA_DNA_ = setting;
}


bool
HBondOptions::exclude_intra_res_protein() const
{
	return exclude_intra_res_protein_;
}


void
HBondOptions::exclude_intra_res_protein( bool const setting )
{
	exclude_intra_res_protein_ = setting;
}


bool
HBondOptions::exclude_intra_res_RNA() const
{
	return exclude_intra_res_RNA_;
}


void
HBondOptions::exclude_intra_res_RNA( bool const setting )
{
	exclude_intra_res_RNA_ = setting;
}


bool
HBondOptions::put_intra_into_total() const
{
	return put_intra_into_total_;
}


void
HBondOptions::put_intra_into_total( bool const setting )
{
	put_intra_into_total_ = setting;
}


bool
HBondOptions::exclude_self_hbonds() const
{
	return exclude_self_hbonds_;
}


void
HBondOptions::exclude_self_hbonds( bool const setting )
{
	exclude_self_hbonds_ = setting;
}


bool
HBondOptions::use_hb_env_dep_DNA() const
{
	return use_hb_env_dep_DNA_;
}


void
HBondOptions::use_hb_env_dep_DNA( bool const setting )
{
	use_hb_env_dep_DNA_ = setting;
}


bool
HBondOptions::use_hb_env_dep() const
{
	return use_hb_env_dep_;
}


void
HBondOptions::use_hb_env_dep( bool const setting )
{
	use_hb_env_dep_ = setting;
}


bool
HBondOptions::smooth_hb_env_dep() const
{
	return smooth_hb_env_dep_;
}


void
HBondOptions::smooth_hb_env_dep( bool const setting )
{
	smooth_hb_env_dep_ = setting;
}


bool
HBondOptions::decompose_bb_hb_into_pair_energies() const
{
	return decompose_bb_hb_into_pair_energies_;
}


void
HBondOptions::decompose_bb_hb_into_pair_energies( bool const setting )
{
	decompose_bb_hb_into_pair_energies_ = setting;
}


bool
HBondOptions::bb_donor_acceptor_check() const
{
	return bb_donor_acceptor_check_;
}


void
HBondOptions::bb_donor_acceptor_check( bool const setting )
{
	bb_donor_acceptor_check_ = setting;
}


std::string const &
HBondOptions::params_database_tag() const
{
	return params_database_tag_;
}


void
HBondOptions::params_database_tag( std::string const & setting )
{
	params_database_tag_ = setting;
}

///pba
bool
HBondOptions::Mbhbond() const
{
	return Mbhbond_;
}


void
HBondOptions::Mbhbond( bool const setting )
{
	Mbhbond_ = setting;
}


bool
HBondOptions::mphbond() const
{
	return mphbond_;
}


void
HBondOptions::mphbond( bool const setting )
{
	mphbond_ = setting;
}

bool HBondOptions::use_sp2_chi_penalty() const
{
	return use_sp2_chi_penalty_;
}

void HBondOptions::use_sp2_chi_penalty( bool setting )
{
	use_sp2_chi_penalty_ = setting;
}

Real HBondOptions::sp2_BAH180_rise() const { return sp2_BAH180_rise_; }
void HBondOptions::sp2_BAH180_rise( Real setting ) { sp2_BAH180_rise_ = setting; }

Real HBondOptions::sp2_outer_width() const { return sp2_outer_width_; }
void HBondOptions::sp2_outer_width( Real setting ) { sp2_outer_width_ = setting; }

bool HBondOptions::measure_sp3acc_BAH_from_hvy() const { return measure_sp3acc_BAH_from_hvy_; }
void HBondOptions::measure_sp3acc_BAH_from_hvy( bool setting ) { measure_sp3acc_BAH_from_hvy_ = setting; }

bool HBondOptions::fade_energy() const { return fade_energy_; }
void HBondOptions::fade_energy( bool setting ) { fade_energy_ = setting; }

Real HBondOptions::hbond_energy_shift() const { return hbond_energy_shift_; }
void HBondOptions::hbond_energy_shift( Real setting ) { hbond_energy_shift_ = setting; }

bool HBondOptions::length_dependent_srbb() const { return length_dependent_srbb_; }
void HBondOptions::length_dependent_srbb( bool setting ) { length_dependent_srbb_ = setting; }

Real HBondOptions::length_dependent_srbb_lowscale() const { return ldsrbb_low_scale_; }
void HBondOptions::length_dependent_srbb_lowscale( Real setting ) { ldsrbb_low_scale_ = setting; }

Real HBondOptions::length_dependent_srbb_highscale() const { return ldsrbb_high_scale_; }
void HBondOptions::length_dependent_srbb_highscale( Real setting ) { ldsrbb_high_scale_ = setting; }

Size HBondOptions::length_dependent_srbb_minlength() const { return ldsrbb_minlength_; }
void HBondOptions::length_dependent_srbb_minlength( Size setting ) { ldsrbb_minlength_ = setting; }

Size HBondOptions::length_dependent_srbb_maxlength() const { return ldsrbb_maxlength_; }
void HBondOptions::length_dependent_srbb_maxlength( Size setting ) { ldsrbb_maxlength_ = setting; }

bool HBondOptions::use_hb_env_dep_new() const { return use_hb_env_dep_new_; }
void HBondOptions::use_hb_env_dep_new(bool val) { use_hb_env_dep_new_=val; }

core::Real HBondOptions::hb_env_dep_new_low_scale() const { return hb_env_dep_new_low_scale_; }
void HBondOptions::hb_env_dep_new_low_scale(core::Real val) { hb_env_dep_new_low_scale_=val; }

core::Real HBondOptions::hb_env_dep_new_low_nneigh() const { return hb_env_dep_new_low_nneigh_; }
void HBondOptions::hb_env_dep_new_low_nneigh(core::Real val) { hb_env_dep_new_low_nneigh_=val; }

core::Real HBondOptions::hb_env_dep_new_high_nneigh() const { return hb_env_dep_new_high_nneigh_; }
void HBondOptions::hb_env_dep_new_high_nneigh(core::Real val) { hb_env_dep_new_high_nneigh_=val; }

bool
operator==( HBondOptions const & a, HBondOptions const & b )
{
	return ( ( a.exclude_DNA_DNA_ == b.exclude_DNA_DNA_ ) &&
		( a.exclude_intra_res_protein_ == b.exclude_intra_res_protein_  ) &&
		( a.exclude_intra_res_RNA_ == b.exclude_intra_res_RNA_  ) &&
		( a.put_intra_into_total_ == b.put_intra_into_total_  ) &&
		( a.exclude_self_hbonds_ == b.exclude_self_hbonds_ ) &&
		( a.use_hb_env_dep_ == b.use_hb_env_dep_ ) &&
		( a.use_hb_env_dep_DNA_ == b.use_hb_env_dep_DNA_ ) &&
		( a.smooth_hb_env_dep_ == b.smooth_hb_env_dep_ ) &&
		( a.bb_donor_acceptor_check_ == b.bb_donor_acceptor_check_ ) &&
		( a.decompose_bb_hb_into_pair_energies_ == b.decompose_bb_hb_into_pair_energies_ ) &&
		( a.params_database_tag_ == b.params_database_tag_ ) &&
		( a.use_sp2_chi_penalty_ == b.use_sp2_chi_penalty_ ) &&
		( a.sp2_BAH180_rise_ == b.sp2_BAH180_rise_ ) &&
		( a.sp2_outer_width_ == b.sp2_outer_width_ ) &&
		( a.measure_sp3acc_BAH_from_hvy_ == b.measure_sp3acc_BAH_from_hvy_ ) &&
		( a.fade_energy_ == b.fade_energy_ ) &&
		( a.Mbhbond_ == b.Mbhbond_ ) && //pba
		( a.mphbond_ == b.mphbond_ ) &&
		( a.hbond_energy_shift_ == b.hbond_energy_shift_ ) &&
		( a.length_dependent_srbb_ == b.length_dependent_srbb_) &&
		( a.ldsrbb_low_scale_ == b.ldsrbb_low_scale_) &&
		( a.ldsrbb_high_scale_ == b.ldsrbb_high_scale_) &&
		( a.ldsrbb_minlength_ == b.ldsrbb_minlength_) &&
		( a.ldsrbb_maxlength_ == b.ldsrbb_maxlength_) &&
		( a.use_hb_env_dep_new_ == b.use_hb_env_dep_new_) &&
		( a.hb_env_dep_new_low_scale_ == b.hb_env_dep_new_low_scale_) &&
		( a.hb_env_dep_new_low_nneigh_ == b.hb_env_dep_new_low_nneigh_) &&
		( a.hb_env_dep_new_high_nneigh_ == b.hb_env_dep_new_high_nneigh_) );
}

bool
operator!=( HBondOptions const & a, HBondOptions const & b )
{
	return !( a == b );
}

std::ostream &
operator<< ( std::ostream & out, const HBondOptions & options ){
	options.show( out );
	return out;
}


void
HBondOptions::show( std::ostream & out ) const
{
	out <<"HBondOptions::show: exclude_DNA_DNA: "
		<<( exclude_DNA_DNA_ ? "true" : "false" ) << std::endl;
	out <<"HBondOptions::show: exclude_intra_res_protein_: "
		<<( exclude_intra_res_protein_ ? "true" : "false" ) << std::endl;
	out <<"HBondOptions::show: exclude_intra_res_RNA_: "
		<<( exclude_intra_res_RNA_ ? "true" : "false" ) << std::endl;
	out <<"HBondOptions::show: put_intra_into_total_: "
		<<( put_intra_into_total_ ? "true" : "false" ) << std::endl;
	out <<"HBondOptions::show: exclude_self_hbonds: "
		<<( exclude_self_hbonds_ ? "true" : "false" ) << std::endl;
	out <<"HBondOptions::show: use_hb_env_dep: "
		<<( use_hb_env_dep_ ? "true" : "false" ) << std::endl;
	out <<"HBondOptions::show: use_hb_env_dep_DNA: "
		<<( use_hb_env_dep_DNA_ ? "true" : "false" ) << std::endl;
	out <<"HBondOptions::show: smooth_hb_env_dep: "
		<<( smooth_hb_env_dep_ ? "true" : "false " ) << std::endl;
	out <<"HBondOptions::show: bb_donor_acceptor_check: "
		<<( bb_donor_acceptor_check_ ? "true" : "false " ) << std::endl;
	out <<"HBondOptions::show: decompose_bb_hb_into_pair_energies: "
		<<( decompose_bb_hb_into_pair_energies_ ? "true" : "false" ) << std::endl;
	out <<"HBondOptions::show: params_database_tag_: "
		<< params_database_tag_ << std::endl;
	out <<"HBondOptions::show: use_sp2_chi_penalty_: "
		<<( use_sp2_chi_penalty_ ? "true" : "false" ) << std::endl;
	out <<"HBondOptions::show: sp2_BAH180_rise_: "
		<< sp2_BAH180_rise_ << std::endl;
	out <<"HBondOptions::show: sp2_outer_width_: "
		<< sp2_outer_width_ << std::endl;
	out <<"HBondOptions::show: measure_sp3acc_BAH_from_hvy_: "
		<<( measure_sp3acc_BAH_from_hvy_ ? "true" : "false" ) << std::endl;
	out <<"HBondOptions::show: fade_energy_: "
		<< fade_energy_ << std::endl;
	out <<"HBondOptions::show: Mbhbond: "
		<<( Mbhbond_ ? "true" : "false " ) << std::endl; //pba
	out << "HbondOptions::show: mphbond: "
		<<( mphbond_ ? "true" : "false" ) << std::endl;
	out <<"HBondOptions::show: hbond_energy_shift: " << hbond_energy_shift_ << std::endl;
}

}
}
}

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::hbonds::HBondOptions::save( Archive & arc ) const {
	arc( CEREAL_NVP( exclude_DNA_DNA_ ) ); // _Bool
	arc( CEREAL_NVP( exclude_intra_res_protein_ ) ); // _Bool
	arc( CEREAL_NVP( exclude_intra_res_RNA_ ) ); // _Bool
	arc( CEREAL_NVP( put_intra_into_total_ ) ); // _Bool
	arc( CEREAL_NVP( exclude_self_hbonds_ ) ); // _Bool
	arc( CEREAL_NVP( use_hb_env_dep_ ) ); // _Bool
	arc( CEREAL_NVP( use_hb_env_dep_DNA_ ) ); // _Bool
	arc( CEREAL_NVP( smooth_hb_env_dep_ ) ); // _Bool
	arc( CEREAL_NVP( bb_donor_acceptor_check_ ) ); // _Bool
	arc( CEREAL_NVP( decompose_bb_hb_into_pair_energies_ ) ); // _Bool
	arc( CEREAL_NVP( params_database_tag_ ) ); // std::string
	arc( CEREAL_NVP( use_sp2_chi_penalty_ ) ); // _Bool
	arc( CEREAL_NVP( sp2_BAH180_rise_ ) ); // Real
	arc( CEREAL_NVP( sp2_outer_width_ ) ); // Real
	arc( CEREAL_NVP( measure_sp3acc_BAH_from_hvy_ ) ); // _Bool
	arc( CEREAL_NVP( fade_energy_ ) ); // _Bool
	arc( CEREAL_NVP( Mbhbond_ ) ); // _Bool
	arc( CEREAL_NVP( mphbond_ ) ); // _Bool
	arc( CEREAL_NVP( hbond_energy_shift_ ) ); // Real
	arc( CEREAL_NVP( length_dependent_srbb_ ) ); // _Bool
	arc( CEREAL_NVP( ldsrbb_low_scale_ ) ); // Real
	arc( CEREAL_NVP( ldsrbb_high_scale_ ) ); // Real
	arc( CEREAL_NVP( ldsrbb_minlength_ ) ); // Size
	arc( CEREAL_NVP( ldsrbb_maxlength_ ) ); // Size
	arc( CEREAL_NVP( use_hb_env_dep_new_ ) ); // _Bool
	arc( CEREAL_NVP( hb_env_dep_new_low_scale_ ) ); // core::Real
	arc( CEREAL_NVP( hb_env_dep_new_low_nneigh_ ) ); // core::Real
	arc( CEREAL_NVP( hb_env_dep_new_high_nneigh_ ) ); // core::Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::hbonds::HBondOptions::load( Archive & arc ) {
	arc( exclude_DNA_DNA_ ); // _Bool
	arc( exclude_intra_res_protein_ ); // _Bool
	arc( exclude_intra_res_RNA_ ); // _Bool
	arc( put_intra_into_total_ ); // _Bool
	arc( exclude_self_hbonds_ ); // _Bool
	arc( use_hb_env_dep_ ); // _Bool
	arc( use_hb_env_dep_DNA_ ); // _Bool
	arc( smooth_hb_env_dep_ ); // _Bool
	arc( bb_donor_acceptor_check_ ); // _Bool
	arc( decompose_bb_hb_into_pair_energies_ ); // _Bool
	arc( params_database_tag_ ); // std::string
	arc( use_sp2_chi_penalty_ ); // _Bool
	arc( sp2_BAH180_rise_ ); // Real
	arc( sp2_outer_width_ ); // Real
	arc( measure_sp3acc_BAH_from_hvy_ ); // _Bool
	arc( fade_energy_ ); // _Bool
	arc( Mbhbond_ ); // _Bool
	arc( mphbond_ ); // _Bool
	arc( hbond_energy_shift_ ); // Real
	arc( length_dependent_srbb_ ); // _Bool
	arc( ldsrbb_low_scale_ ); // Real
	arc( ldsrbb_high_scale_ ); // Real
	arc( ldsrbb_minlength_ ); // Size
	arc( ldsrbb_maxlength_ ); // Size
	arc( use_hb_env_dep_new_ ); // _Bool
	arc( hb_env_dep_new_low_scale_ ); // core::Real
	arc( hb_env_dep_new_low_nneigh_ ); // core::Real
	arc( hb_env_dep_new_high_nneigh_ ); // core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::hbonds::HBondOptions );
CEREAL_REGISTER_TYPE( core::scoring::hbonds::HBondOptions )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_hbonds_HBondOptions )
#endif // SERIALIZATION
