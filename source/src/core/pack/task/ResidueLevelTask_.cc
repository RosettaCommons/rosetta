// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/PackerTask_.hh
/// @brief  Implementation class for task class to describe packer's behavior header
/// Almost all of rosetta needs to use packer tasks, but very little of rosetta needs
/// to see how it behaves internally.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Steven Lewis (smlewi@gmail.com)


//Unit Headers
#include <core/pack/task/ResidueLevelTask_.hh>

// Package Headers
#include <core/pack/rotamer_set/RotamerSetOperation.hh>
#include <core/pack/rotamer_set/RotamerCouplings.hh>
#include <core/pack/rotamer_set/RotamerLinks.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/task/IGEdgeReweightContainer.hh>
#include <core/pack/task/rna/RNA_ResidueLevelTask.hh>

//Project Headers
#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/pose/Pose.hh>
#include <basic/options/option.hh>
#include <core/pose/PDBInfo.hh>
#include <core/id/SequenceMapping.hh>
#include <core/pose/util.hh>

// Objexx headers
#include <ObjexxFCL/format.hh>

//Utility Headers
#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/options/keys/OptionKeyList.hh>

// basic headers
#include <basic/Tracer.hh>

// c++ headers
#include <iostream>

// option key includes
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>

using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

#ifdef    SERIALIZATION
// Project serialization headers
#include <core/chemical/ResidueType.srlz.hh>

// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/list.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#endif // SERIALIZATION

namespace core {
namespace pack {
namespace task {

static THREAD_LOCAL basic::Tracer T( "core.pack.task", basic::t_info );

/// @details ResidueLevelTask constructor has following defaults:
///all ex set to false with zero sample level
///position is packable and designable to any of the canonical aa's,
///with variant matching (for termini, etc)
///current rotamer is not included for packer.
///bump check is deactivated by default.
ResidueLevelTask_::ResidueLevelTask_(
	conformation::Residue const & original_residue,
	pose::Pose const & pose
)
:
	include_current_( false ),
	adducts_( true ),
	original_residue_type_set_( pose.residue_type_set_for_pose( original_residue.type().mode() ) ),
	original_residue_type_( original_residue.type().get_self_weak_ptr() ),
	target_residue_type_(/* 0 */),
	designing_( original_residue.is_protein() || original_residue.is_peptoid() ), // default -- design at all protein residues
	repacking_( true ),
	optimize_H_mode_( false ),
	preserve_c_beta_( false ),
	flip_HNQ_( false ),
	fix_his_tautomer_( false ),
	include_virtual_side_chain_( false ),
	disabled_( false ),
	design_disabled_( false ),
	sample_proton_chi_( true ),
	ex1_( false ),
	ex2_( false ),
	ex3_( false ),
	ex4_( false ),
	ex1aro_( false ),
	ex2aro_( false ),
	ex1aro_exposed_( false ),
	ex2aro_exposed_( false ),
	ex1_sample_level_( NO_EXTRA_CHI_SAMPLES ),
	ex2_sample_level_( NO_EXTRA_CHI_SAMPLES ),
	ex3_sample_level_( NO_EXTRA_CHI_SAMPLES ),
	ex4_sample_level_( NO_EXTRA_CHI_SAMPLES ),
	ex1aro_sample_level_( NO_EXTRA_CHI_SAMPLES ),
	ex2aro_sample_level_( NO_EXTRA_CHI_SAMPLES ),
	ex1aro_exposed_sample_level_( NO_EXTRA_CHI_SAMPLES ),
	ex2aro_exposed_sample_level_( NO_EXTRA_CHI_SAMPLES ),
	exdna_sample_level_( NO_EXTRA_CHI_SAMPLES ),
	extrachi_cutoff_( EXTRACHI_CUTOFF_LIMIT ),
	operate_on_ex1_( false ),
	operate_on_ex2_( false ),
	operate_on_ex3_( false ),
	operate_on_ex4_( false )
{
	using namespace chemical;

	if ( original_residue.is_protein() || original_residue.is_peptoid() ) {

		//default: all amino acids at all positions -- additional "and" operations will remove
		// amino acids from the list of allowed ones
		//no rule yet to treat chemically modified aa's differently
		ResidueType const & match_residue_type( original_residue_type_set_->get_residue_type_with_variant_removed( original_residue.type(), chemical::VIRTUAL_SIDE_CHAIN ) );
		for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) {
			for ( Size jj = 1; jj <= match_residue_type.variant_types().size(); ++jj ) {}
			ResidueTypeCOPs const & aas( original_residue_type_set_->get_all_types_with_variants_aa( AA( ii ), match_residue_type.variant_types(), pH_mode_exceptions() ) );
			for ( ResidueTypeCOPs::const_iterator
					aas_iter = aas.begin(),
					aas_end = aas.end(); aas_iter != aas_end; ++aas_iter ) {
				allowed_residue_types_.push_back( *aas_iter );
			}
		}
		// allow noncanonical AAs and D-amino acids to be repacked
		if ( original_residue.aa() == aa_unk || core::chemical::is_canonical_D_aa( original_residue.aa() ) ) {
			allowed_residue_types_.push_back( original_residue.type().get_self_ptr() );
		}

		// allow metapatched residues to at least repack???
		if ( allowed_residue_types_.empty() ) {
			allowed_residue_types_.push_back( original_residue.type().get_self_ptr() );
		}

	} else if ( original_residue.is_DNA() ) {

		ResidueTypeCOPs dna_types = ResidueTypeFinder( *original_residue_type_set_ ).variants( original_residue.type().variant_types() ).
			base_property( DNA ).variant_exceptions( utility::tools::make_vector1( ADDUCT_VARIANT ) ).get_all_possible_residue_types();
		for ( Size n = 1; n <= dna_types.size(); n++ ) allowed_residue_types_.push_back( dna_types[ n ] );

	} else if ( original_residue.is_RNA() ) {

		ResidueType const & match_residue_type(
			original_residue_type_set_->get_residue_type_with_variant_removed( original_residue.type(), chemical::VIRTUAL_O2PRIME_HYDROGEN ) );
		allowed_residue_types_.push_back( match_residue_type.get_self_ptr() );

	} else {
		// for non-amino acids, default is to include only the existing residuetype
		allowed_residue_types_.push_back( original_residue.type().get_self_ptr() );
	}
	if ( original_residue.is_RNA() ) rna_task_ = rna::RNA_ResidueLevelTaskOP( new rna::RNA_ResidueLevelTask );
	// The intention is for all ResidueTasks to *start off* as repackable.
	// Some, like protein AAs and DNA, also start off designable.
	determine_if_designing();
	determine_if_repacking();
	// assert packable by default
	runtime_assert( being_packed() );
	// This is the same assertion, really, but coded directly
	runtime_assert( ! allowed_residue_types_.empty() );
}

ResidueLevelTask_::ResidueLevelTask_() {}

ResidueLevelTask_::~ResidueLevelTask_() {}

ExtraRotSample
ResidueLevelTask_::extrachi_sample_level(
	bool buried,
	int chi,
	chemical::ResidueType const & concrete_residue
) const
{

	if ( concrete_residue.is_DNA() ) {
		runtime_assert( chi == 1 );
		return exdna_sample_level_;
	}
	if ( concrete_residue.is_aromatic()  && chi <= 2 ) {
		if ( chi == 1 ) {
			if (  buried ) {
				return ex1aro_sample_level_;
			} else {
				return ex1aro_exposed_sample_level_;
			}
		} else {
			if ( buried ) {
				return ex2aro_sample_level_;
			} else {
				return ex2aro_exposed_sample_level_;
			}
		}
	} else {
		if (  buried && chi <= 4 ) {
			if ( chi == 1 ) {
				return ex1_sample_level_;
			} else if ( chi == 2 ) {
				return ex2_sample_level_;
			} else if ( chi == 3 ) {
				return ex3_sample_level_;
			} else {
				return ex4_sample_level_;
			}
		}
	}
	return NO_EXTRA_CHI_SAMPLES;
}

void
ResidueLevelTask_::initialize_from_command_line()
{
	initialize_from_options( basic::options::option );
}

void ResidueLevelTask_::initialize_from_options( utility::options::OptionCollection const & options )
{
	using namespace basic::options;
	using namespace OptionKeys::packing;

	initialize_extra_rotamer_flags_from_options( options );

	if ( options[ use_input_sc ] ) or_include_current( true );
	if ( options[ OptionKeys::packing::preserve_c_beta ].user() ) or_preserve_c_beta( true );
	if ( options[ OptionKeys::packing::prevent_repacking ].value() ) prevent_repacking();

}

void
ResidueLevelTask_::initialize_extra_rotamer_flags_from_command_line()
{
	initialize_extra_rotamer_flags_from_options( basic::options::option );
}

void
ResidueLevelTask_::initialize_extra_rotamer_flags_from_options( utility::options::OptionCollection const & options )
{
	using namespace basic::options;
	using namespace OptionKeys::packing;

	if ( options[ ex1::ex1 ].user() ) { or_ex1( options[ ex1::ex1 ].value() ); }
	if ( options[ ex2::ex2 ].user() ) { or_ex2( options[ ex2::ex2 ].value() ); }
	if ( options[ ex3::ex3 ].user() ) { or_ex3( options[ ex3::ex3 ].value() ); }
	if ( options[ ex4::ex4 ].user() ) { or_ex4( options[ ex4::ex4 ].value() ); }

	if ( options[ ex1aro::ex1aro ].user() ) {
		or_ex1aro( options[ ex1aro::ex1aro ].value() );
	}

	if ( options[ ex2aro::ex2aro ].user() ) {
		or_ex2aro( options[ ex2aro::ex2aro ].value() );
	}

	if ( options[ ex1aro_exposed::ex1aro_exposed ].user() ) {
		or_ex1aro_exposed( options[ ex1aro_exposed::ex1aro_exposed ].value() );
	}

	if ( options[ ex2aro_exposed::ex2aro_exposed ].user() ) {
		or_ex2aro_exposed( options[ ex2aro_exposed::ex2aro_exposed ].value() );
	}

	if ( options[ ex1::level ].user() ) {
		or_ex1_sample_level( static_cast<ExtraRotSample>( options[ ex1::level ].value() ));
	}

	if ( options[ ex2::level ].user() ) {
		or_ex2_sample_level( static_cast<ExtraRotSample>( options[ ex2::level ].value() ));
	}

	if ( options[ ex3::level ].user() ) {
		or_ex3_sample_level( static_cast<ExtraRotSample>( options[ ex3::level ].value() ));
	}

	if ( options[ ex4::level ].user() ) {
		or_ex4_sample_level( static_cast<ExtraRotSample>( options[ ex4::level ].value() ));
	}

	if ( options[ ex1aro::level ].user() ) {
		or_ex1aro_sample_level( static_cast<ExtraRotSample>( options[ ex1aro::level ].value() ));
	}

	if ( options[ ex2aro::level ].user() ) {
		or_ex2aro_sample_level( static_cast<ExtraRotSample>( options[ ex2aro::level ].value() ));
	}

	if ( options[ ex1aro_exposed::level ].user() ) {
		or_ex1aro_exposed_sample_level( static_cast<ExtraRotSample>( options[ ex1aro_exposed::level ].value() ));
	}

	if ( options[ ex2aro_exposed::level ].user() ) {
		or_ex2aro_exposed_sample_level( static_cast<ExtraRotSample>( options[ ex2aro_exposed::level ].value() ));
	}

	if ( options[ ex1::operate ].user() ) { or_operate_on_ex1( options[ ex1::operate ]() ); }
	if ( options[ ex2::operate ].user() ) { or_operate_on_ex2( options[ ex2::operate ]() ); }
	if ( options[ ex3::operate ].user() ) { or_operate_on_ex3( options[ ex3::operate ]() ); }
	if ( options[ ex4::operate ].user() ) { or_operate_on_ex4( options[ ex4::operate ]() ); }

	// options to control dna rotamer building
	if ( options[ exdna::exdna ].user() && options[ exdna::level ].user() ) {
		if ( options[ exdna::exdna ].value() && options[ exdna::level ].value() ) {
			or_exdna_sample_level( static_cast<ExtraRotSample>( options[ exdna::level ].value() ) );
		}
	} else if ( options[ exdna::exdna ].user() && !options[ exdna::level ].user() ) {
		if ( options[ exdna::exdna ].value() ) {
			or_exdna_sample_level( EX_ONE_STDDEV );
		}
	} else if ( !options[ exdna::exdna ].user() && options[ exdna::level ].user() ) {
		if ( options[ exdna::level ].value() ) {
			or_exdna_sample_level( static_cast<ExtraRotSample>( options[ exdna::level ].value() ) );
		}
	}

	// extra chi cutoff
	if ( options[ OptionKeys::packing::extrachi_cutoff ].user() ) {
		and_extrachi_cutoff( options[ OptionKeys::packing::extrachi_cutoff ] );
	}
}

void
ResidueLevelTask::list_options_read( utility::options::OptionKeyList & read_options )
{
	using namespace basic::options::OptionKeys;
	read_options
		+ packing::use_input_sc
		+ packing::preserve_c_beta
		+ packing::prevent_repacking
		+ packing::ex1::ex1
		+ packing::ex2::ex2
		+ packing::ex3::ex3
		+ packing::ex4::ex4
		+ packing::ex1aro::ex1aro
		+ packing::ex2aro::ex2aro
		+ packing::ex1aro_exposed::ex1aro_exposed
		+ packing::ex2aro_exposed::ex2aro_exposed
		+ packing::ex1::level
		+ packing::ex2::level
		+ packing::ex3::level
		+ packing::ex4::level
		+ packing::ex1aro::level
		+ packing::ex2aro::level
		+ packing::ex1aro_exposed::level
		+ packing::ex2aro_exposed::level
		+ packing::ex1::operate
		+ packing::ex2::operate
		+ packing::ex3::operate
		+ packing::ex4::operate
		+ packing::exdna::exdna
		+ packing::exdna::level
		+ packing::extrachi_cutoff;
}


void ResidueLevelTask_::or_include_current( bool include_current )
{
	include_current_ |= include_current;
}

bool ResidueLevelTask_::include_current() const
{
	return include_current_;
}

void ResidueLevelTask_::add_behavior( std::string const & behavior )
{
	// probably not worth caring about (rare) redundancies
	behaviors_.push_back( behavior );
}
bool ResidueLevelTask_::has_behavior( std::string const & behavior ) const
{
	return std::find( behaviors_.begin(), behaviors_.end(), behavior ) != behaviors_.end();
}
bool ResidueLevelTask_::has_behavior() const
{
	return behaviors_.empty() ? false : true;
}

void ResidueLevelTask_::target_type( chemical::ResidueTypeCOP type ) {
	bool allowed( std::find( allowed_residue_types_.begin(), allowed_residue_types_.end(), type ) !=
		allowed_residue_types_.end() );
	if ( !allowed ) {
		T.Error << "Target type " << type->name() << " is not an allowed type!" << std::endl;
		utility_exit();
		return;
	}
	target_residue_type_ = type; /// non-commutative if multiple target residue types are set.
}
void ResidueLevelTask_::target_type( chemical::AA aa ) {
	target_type( original_residue_type_set_->get_representative_type_aa( aa ) );
}
void ResidueLevelTask_::target_type( std::string name ) {
	target_type( original_residue_type_set_->name_map( name ).get_self_ptr() );
}

void ResidueLevelTask_::or_adducts( bool setting )
{
	if ( setting ) return;

	for ( ResidueTypeCOPListIter type_iter( allowed_residue_types_.begin() ), type_end(  allowed_residue_types_.end() ); type_iter != type_end; ) {
		if ( (*type_iter)->has_variant_type( chemical::ADDUCT_VARIANT ) ) {
			type_iter = allowed_residue_types_.erase( type_iter );
		} else {
			++type_iter;
		}
	}

	adducts_ = false;
	determine_if_designing();
	determine_if_repacking();
}
bool ResidueLevelTask_::adducts() const { return adducts_; }

void ResidueLevelTask_::or_ex1( bool ex1 )
{
	ex1_ |= ex1;
	refresh_ex1_sample_levels();
}

void ResidueLevelTask_::or_ex2( bool ex2 )
{
	ex2_ |= ex2;
	refresh_ex2_sample_levels();
}

void ResidueLevelTask_::or_ex3( bool ex3 )
{
	ex3_ |= ex3;
	refresh_ex3_sample_levels();
}

void ResidueLevelTask_::or_ex4( bool ex4 )
{
	ex4_ |= ex4;
	refresh_ex4_sample_levels();
}

void ResidueLevelTask_::or_ex1_sample_level( ExtraRotSample ex1_sample_level )
{
	if ( ex1_sample_level_ < ex1_sample_level ) {
		ex1_sample_level_ = ex1_sample_level;
		refresh_ex1_sample_levels();
	}
}

void ResidueLevelTask_::or_exdna_sample_level( ExtraRotSample exdna_sample_level )
{
	if ( exdna_sample_level_ < exdna_sample_level ) {
		exdna_sample_level_ = exdna_sample_level;
		refresh_ex1_sample_levels();
	}
}

void ResidueLevelTask_::or_ex2_sample_level( ExtraRotSample ex2_sample_level )
{
	if ( ex2_sample_level_ < ex2_sample_level ) {
		ex2_sample_level_ = ex2_sample_level;
		refresh_ex2_sample_levels();
	}
}
void ResidueLevelTask_::or_ex3_sample_level( ExtraRotSample ex3_sample_level )
{
	if ( ex3_sample_level_ < ex3_sample_level ) {
		ex3_sample_level_ = ex3_sample_level;
		refresh_ex3_sample_levels();
	}
}
void ResidueLevelTask_::or_ex4_sample_level( ExtraRotSample ex4_sample_level )
{
	if ( ex4_sample_level_ < ex4_sample_level ) {
		ex4_sample_level_ = ex4_sample_level;
		refresh_ex4_sample_levels();
	}
}

void ResidueLevelTask_::or_ex1aro( bool ex1aro )
{
	ex1aro_ |= ex1aro;
	refresh_ex1_sample_levels();
}

void ResidueLevelTask_::or_ex2aro( bool ex2aro )
{
	ex2aro_ |= ex2aro;
	refresh_ex2_sample_levels();
}

void ResidueLevelTask_::or_ex1aro_exposed( bool ex1aro_exposed )
{
	ex1aro_exposed_ |= ex1aro_exposed;
	refresh_ex1_sample_levels();
}

void ResidueLevelTask_::or_ex2aro_exposed( bool ex2aro_exposed )
{
	ex2aro_exposed_ |= ex2aro_exposed;
	refresh_ex2_sample_levels();
}

void ResidueLevelTask_::or_ex1aro_sample_level( ExtraRotSample ex1aro_sample_level )
{
	if ( ex1aro_sample_level_ < ex1aro_sample_level ) {
		ex1aro_sample_level_ = ex1aro_sample_level;
		refresh_ex1_sample_levels();
	}
}

void ResidueLevelTask_::or_ex2aro_sample_level( ExtraRotSample ex2aro_sample_level )
{
	if ( ex2aro_sample_level_ < ex2aro_sample_level ) {
		ex2aro_sample_level_ = ex2aro_sample_level;
		refresh_ex2_sample_levels();
	}
}

void ResidueLevelTask_::or_ex1aro_exposed_sample_level(
	ExtraRotSample ex1aro_exposed_sample_level
)
{
	if ( ex1aro_exposed_sample_level_ < ex1aro_exposed_sample_level ) {
		ex1aro_exposed_sample_level_ = ex1aro_exposed_sample_level;
		refresh_ex1_sample_levels();
	}
}

void ResidueLevelTask_::or_ex2aro_exposed_sample_level(
	ExtraRotSample ex2aro_exposed_sample_level
)
{
	if ( ex2aro_exposed_sample_level_ < ex2aro_exposed_sample_level ) {
		ex2aro_exposed_sample_level_ = ex2aro_exposed_sample_level;
		refresh_ex2_sample_levels();
	}
}

void ResidueLevelTask_::or_operate_on_ex1( bool operate )
{ operate_on_ex1_ |= operate; }
void ResidueLevelTask_::or_operate_on_ex2( bool operate )
{ operate_on_ex2_ |= operate; }
void ResidueLevelTask_::or_operate_on_ex3( bool operate )
{ operate_on_ex3_ |= operate; }
void ResidueLevelTask_::or_operate_on_ex4( bool operate )
{ operate_on_ex4_ |= operate; }


void ResidueLevelTask_::sample_proton_chi( bool setting )
{
	sample_proton_chi_ = setting;
}

bool ResidueLevelTask_::sample_proton_chi() const
{
	return sample_proton_chi_;
}

bool ResidueLevelTask_::ex1() const
{
	return ex1_;
}
bool ResidueLevelTask_::ex2() const
{
	return ex2_;
}
bool ResidueLevelTask_::ex3() const
{
	return ex3_;
}
bool ResidueLevelTask_::ex4() const
{
	return ex4_;
}

ExtraRotSample ResidueLevelTask_::ex1_sample_level() const
{
	return ex1_sample_level_;
}
ExtraRotSample ResidueLevelTask_::ex2_sample_level() const
{
	return ex2_sample_level_;
}
ExtraRotSample ResidueLevelTask_::ex3_sample_level() const
{
	return ex3_sample_level_;
}
ExtraRotSample ResidueLevelTask_::ex4_sample_level() const
{
	return ex4_sample_level_;
}

bool ResidueLevelTask_::ex1aro() const
{
	return ex1aro_;
}
bool ResidueLevelTask_::ex2aro() const
{
	return ex2aro_;
}
bool ResidueLevelTask_::ex1aro_exposed() const
{
	return ex1aro_exposed_;
}
bool ResidueLevelTask_::ex2aro_exposed() const
{
	return ex2aro_exposed_;
}

ExtraRotSample ResidueLevelTask_::ex1aro_sample_level() const
{
	return ex1aro_sample_level_;
}
ExtraRotSample ResidueLevelTask_::ex2aro_sample_level() const
{
	return ex2aro_sample_level_;
}
ExtraRotSample ResidueLevelTask_::ex1aro_exposed_sample_level() const
{
	return ex1aro_exposed_sample_level_;
}
ExtraRotSample ResidueLevelTask_::ex2aro_exposed_sample_level() const
{
	return ex2aro_exposed_sample_level_;
}

ExtraRotSample ResidueLevelTask_::exdna_sample_level() const
{
	return exdna_sample_level_;
}

bool ResidueLevelTask_::operate_on_ex1() const
{
	return operate_on_ex1_;
}
bool ResidueLevelTask_::operate_on_ex2() const
{
	return operate_on_ex2_;
}
bool ResidueLevelTask_::operate_on_ex3() const
{
	return operate_on_ex3_;
}
bool ResidueLevelTask_::operate_on_ex4() const
{
	return operate_on_ex4_;
}

void ResidueLevelTask_::or_optimize_h( bool setting )
{
	optimize_H_mode_ |= setting;
	if ( optimize_H_mode_ ) {
		// remove all amino acids that do not match the native amino acid
		// there is no design in optimize_H_mode_
		restrict_to_repacking();
	}
}

bool ResidueLevelTask_::optimize_h() const
{
	return optimize_H_mode_;
}

void ResidueLevelTask_::or_preserve_c_beta( bool setting )
{
	preserve_c_beta_ |= setting;
}

bool ResidueLevelTask_::preserve_c_beta() const
{
	return preserve_c_beta_;
}

void ResidueLevelTask_::or_flip_HNQ( bool setting )
{
	flip_HNQ_ |= setting;
	or_optimize_h( setting );
}

bool ResidueLevelTask_::flip_HNQ() const
{
	return flip_HNQ_;
}

/// @details this function forces a fixed histidine tautomer by removing the alternate tautomer from the ResidueTypesCAPList.  The fix_his_tautomer_ boolean is maintained for reference that this has been done, but the boolean is not the source of the effect.
void ResidueLevelTask_::or_fix_his_tautomer( bool setting )
{
	//TODO: Modify this function to allow D-amino acids (particularly D-his)!
	if ( !setting || original_residue_type_->aa() != chemical::aa_his ) return;
	fix_his_tautomer_ |= setting;

	//logic: iterate through ResidueTypeCOPList.  Remove all ResidueTypeCOPs that ARE histidine but are not of the same tautomer as this residue's original_residue_type_.
	//but how do we decide tautomer?  ResidueType offers two options: string comparisons of the names, and checking which atoms are present.  We're going with the latter.

	//which hydrogen atoms are present in the original residue?  We'll check both and vainly hope this won't break if somone tries protonation variants.
	bool HIS(original_residue_type_->has(" HE2")); //HIS
	bool HIS_D(original_residue_type_->has(" HD1")); //HIS_D

	for ( ResidueTypeCOPListIter
			allowed_iter = allowed_residue_types_.begin(),
			iter_next = allowed_residue_types_.begin(),
			allowed_end = allowed_residue_types_.end();
			allowed_iter != allowed_end;  /* no increment: deletion + iterator incrementing = segfault! */ ) {
		iter_next = allowed_iter;
		++iter_next;

		if ( (*allowed_iter)->aa() == chemical::aa_his ) { //we only want to look at histidines
			//which hydrogen atoms are present in the residue we're checking?
			bool comp_HIS((*allowed_iter)->has(" HE2")); //HIS
			bool comp_HIS_D((*allowed_iter)->has(" HD1")); //HIS_D

			if ( (HIS != comp_HIS) || (HIS_D != comp_HIS_D) ) {
				allowed_residue_types_.erase( allowed_iter );
			}
		}
		allowed_iter = iter_next;
	}

}

bool ResidueLevelTask_::fix_his_tautomer() const
{
	return fix_his_tautomer_;
}


/// @details this function samples a protein conformation with a virtualized side chain.
void ResidueLevelTask_::or_include_virtual_side_chain( bool setting )
{
	include_virtual_side_chain_ |= setting;
}

bool ResidueLevelTask_::include_virtual_side_chain() const
{
	return include_virtual_side_chain_;
}

/// @details and operation -- min -- move only toward a lower cutoff for #neighbors w/i 10A that
///qualify a residue to be considered buried.
void ResidueLevelTask_::and_extrachi_cutoff( Size num_neighbors_to_be_called_buried )
{
	if ( extrachi_cutoff_ > num_neighbors_to_be_called_buried ) {
		extrachi_cutoff_ = num_neighbors_to_be_called_buried;
	}
}

Size ResidueLevelTask_::extrachi_cutoff() const
{
	return extrachi_cutoff_;
}

// remove all ResidueTypes from the list of allowed residue types
void ResidueLevelTask_::prevent_repacking()
{
	disabled_ = true;
	design_disabled_ = true;
	allowed_residue_types_.clear();
	determine_if_designing();
	determine_if_repacking();
	mode_tokens_.push_back("NATRO");
}


/// @details contract (and) the list of available aas for canonical aa's
///if an amino acid is not present (false) in the boolean vector, then do not allow it at this position
///boolean vector is based on the aa enum; see another example with PIKAA.
///The boolean vector is a 20-length vector in alphabetical order by one-letter code.
void
ResidueLevelTask_::restrict_absent_canonical_aas( utility::vector1< bool > const & allowed_aas )
{
	Size num_allowed = std::count( allowed_aas.begin(), allowed_aas.end(), true);
	std::ostringstream aas;
	if ( num_allowed == 0 ) {
		mode_tokens_.push_back("EMPTY");
	} else if ( num_allowed == chemical::num_canonical_aas ) {
		// this doesn't constrain anything...
	} else if ( num_allowed < chemical::num_canonical_aas/2 ) {
		mode_tokens_.push_back("PIKAA");
		for ( Size i = 1; i <= chemical::num_canonical_aas; ++i ) {
			if ( allowed_aas[ i ] ) {
				aas << chemical::oneletter_code_from_aa( static_cast<chemical::AA>(i) );
			}
		}
	} else {
		mode_tokens_.push_back("NOTAA");
		for ( Size i = 1; i <= chemical::num_canonical_aas; ++i ) {
			if ( !allowed_aas[ i ] ) {
				aas << chemical::oneletter_code_from_aa( static_cast<chemical::AA>(i) );
			}
		}
	}
	mode_tokens_.push_back(aas.str() );
	do_restrict_absent_canonical_aas( allowed_aas );
}

void
ResidueLevelTask_::restrict_absent_canonical_aas( utility::vector1< bool > const & allowed_aas, std::string const & mode )
{
	mode_tokens_.push_back( mode );
	do_restrict_absent_canonical_aas( allowed_aas );
}

//The boolean vector is a 20-length vector in alphabetical order by one-letter code.
void
ResidueLevelTask_::do_restrict_absent_canonical_aas( utility::vector1< bool > const & allowed_aas )
{
	runtime_assert( allowed_aas.size() == chemical::num_canonical_aas );

	for ( ResidueTypeCOPListIter
			allowed_iter( allowed_residue_types_.begin() ),
			allowed_end( allowed_residue_types_.end() );
			allowed_iter != allowed_end; ) {

		if ( ( (*allowed_iter)->aa() <= chemical::num_canonical_aas ) && ( ! allowed_aas[ (*allowed_iter)->aa() ] ) ) {
			allowed_iter = allowed_residue_types_.erase( allowed_iter );
		} else {
			++allowed_iter;
		}
	}

	determine_if_designing();
	determine_if_repacking();
}

//@details Same behavior as restrict_absent_canonical_aas except that it always allows the native aa at a position even if it is not included in the allowed residues
void
ResidueLevelTask_::restrict_nonnative_canonical_aas( utility::vector1< bool > const & allowed_aas)
{
	runtime_assert( allowed_aas.size() == chemical::num_canonical_aas );

	for ( ResidueTypeCOPListIter
			allowed_iter = allowed_residue_types_.begin(),
			allowed_end = allowed_residue_types_.end();
			allowed_iter != allowed_end; ) {

		//checks if not in the allowed list and not the original type
		if ( ( (*allowed_iter)->aa() <= chemical::num_canonical_aas ) && ( ! allowed_aas[ (*allowed_iter)->aa() ]) && ( ! is_original_type( *allowed_iter ) ) ) {
			allowed_iter = allowed_residue_types_.erase( allowed_iter );
		} else {
			++allowed_iter;
		}
	}

	determine_if_designing();
	determine_if_repacking();
}


/// @details contract (and) the list of available nas for canonical na's
///if a nucleic acid is not present in the vector, then do not allow it at this position
void
ResidueLevelTask_::restrict_absent_nas(
	utility::vector1< chemical::AA > const & keep_nas
)
{
	typedef utility::vector1< chemical::AA > AAs;
	for ( ResidueTypeCOPListIter
			type_itr = allowed_residue_types_.begin(), next_itr = allowed_residue_types_.begin(),
			end_itr = allowed_residue_types_.end(); type_itr != end_itr;
			/* no increment: deletion + iterator incrementing = segfault! */ ) {
		next_itr = type_itr;
		++next_itr;
		bool keep(false);
		for ( AAs::const_iterator na( keep_nas.begin() ); na != keep_nas.end(); ++na ) {
			if ( (*type_itr)->aa() == *na ) { keep = true; break; }
		}
		if ( ! keep ) allowed_residue_types_.erase( type_itr );
		type_itr = next_itr;
	}
	determine_if_designing();
	determine_if_repacking();


	std::ostringstream nas;
	for ( utility::vector1 < chemical::AA >::const_iterator na = keep_nas.begin(); na != keep_nas.end(); ++na ) {
		// illegal for DNA: single-letter codes redundant
		//nas << chemical::oneletter_code_from_aa( *na );
		nas << chemical::name_from_aa( *na );
	}
	mode_tokens_.push_back( "NA" );
	mode_tokens_.push_back( nas.str() );

}

/// @details removes all residues from the allowed residue types list, except the one that matches
///the original residue; this means only rotameric and not sequence changes are allowed
///if the original residue type has been disabled elsewhere, this function will prevent repacking at
///that residue.
void
ResidueLevelTask_::restrict_to_repacking()
{
	for ( ResidueTypeCOPListIter
			allowed_iter = allowed_residue_types_.begin(),
			allowed_end = allowed_residue_types_.end();
			allowed_iter != allowed_end; ) {

		if ( ! is_original_type( *allowed_iter ) ) {
			allowed_iter = allowed_residue_types_.erase( allowed_iter );
		} else {
			++allowed_iter;
		}
	}

	design_disabled_ = true;
	determine_if_designing();
	determine_if_repacking();
	mode_tokens_.push_back("NATAA");
}

bool ResidueLevelTask_::is_original_type( chemical::ResidueTypeCOP type ) const
{
	if ( original_residue_type_->aa() == chemical::aa_unk ) {
		// unknown aa; go with interchangeability_group
		return ( type->interchangeability_group() == original_residue_type_->interchangeability_group() );
	} else if ( fix_his_tautomer_ && (original_residue_type_->aa() == chemical::aa_his) ) {
		// the only way to distinguish HIS and HIS_D is by their full name
		// this will still fail for cases where variants mismatch but the tautomer does not
		return( type->name() == original_residue_type_->name() );
	} else {
		return ( type->aa() == original_residue_type_->aa() );
	}
}

chemical::ResidueTypeSetCOP
ResidueLevelTask_::get_original_residue_set() const {
	return original_residue_type_set_;
}

chemical::AA const &
ResidueLevelTask_::get_original_residue() const {
	return original_residue_type_->aa();
}

/// @details expand (or) the list of available aa's for non-cannonicals
/// looking for ResidueTypes that share the input "interchangeability_group" id
void ResidueLevelTask_::allow_noncanonical_aa(
	std::string const & interchangeability_group,
	chemical::ResidueTypeSet const & residue_set
)
{
	if ( disabled_ || design_disabled_ ) return;

	using namespace core::chemical;

	for ( Size ii = chemical::num_canonical_aas + 1; ii <= chemical::num_aa_types; ++ii ) {
		ResidueTypeCOPs const & aas( ResidueTypeFinder( residue_set ).aa( AA( ii ) ).variants( original_residue_type_->variant_types() ).interchangeability_group( interchangeability_group ).get_all_possible_residue_types() );
		for ( ResidueTypeCOPs::const_iterator
				aas_iter = aas.begin(),
				aas_end = aas.end(); aas_iter != aas_end; ++aas_iter ) {
			if ( std::find( allowed_residue_types_.begin(), allowed_residue_types_.end(), *aas_iter ) == allowed_residue_types_.end() /* haven't already added it */ ) {
				allowed_residue_types_.push_back( *aas_iter );
			}
		}
	}

	mode_tokens_.push_back("NC");
	mode_tokens_.push_back( interchangeability_group );

	determine_if_designing();
	determine_if_repacking();
}

/// @details Calls the overloaded allow_noncanonical_aas method using the same ResidueTypeSet as original_residue_type_
void ResidueLevelTask_::allow_noncanonical_aa( std::string const & interchangeability_group )
{
	allow_noncanonical_aa( interchangeability_group, *original_residue_type_set_ );
}

/// @details Calls the overloaded allow_noncanonical_aas method using the same ResidueTypeSet as original_residue_type_
/// taking as the interchangeability group what is returned by the "name_from_aa" function.
void ResidueLevelTask_::allow_noncanonical_aa( chemical::AA aa )
{
	if ( aa <= chemical::num_canonical_aas ) {
		T.Warning << "Canonical amino acid, " << aa << ", given in a call to allow_noncanonical_aa has no effect." << std::endl;
		return;
	}
	allow_noncanonical_aa( name_from_aa(aa) );
}

/// @details expand (or) the list of available aa's for non-cannonicals
/// looking for ResidueTypes that share the input "interchangeability_group" id
void ResidueLevelTask_::disallow_noncanonical_aas()
{
	for ( ResidueTypeCOPList::iterator
			aas_iter( allowed_residue_types_.begin() ),
			aas_end( allowed_residue_types_.end() );
			aas_iter != aas_end; ) {

		if ( (*aas_iter)->aa() > chemical::num_canonical_aas ) {
			aas_iter = allowed_residue_types_.erase( aas_iter );
		} else {
			++aas_iter;
		}
	}

	mode_tokens_.push_back("EMPTY");

	determine_if_designing();
	determine_if_repacking();
}

void
ResidueLevelTask_::allow_aa(
	chemical::AA const & aa
)
{
	if ( disabled_ || design_disabled_ ) return;
	disabled_ = false;
	design_disabled_ = false;

	chemical::ResidueTypeCOPs const aas( original_residue_type_set_->get_all_types_with_variants_aa( aa, original_residue_type_->variant_types() ) );

	for ( chemical::ResidueTypeCOPs::const_iterator
			aas_iter = aas.begin(), aas_end = aas.end(); aas_iter != aas_end; ++aas_iter ) {
		if ( std::find( allowed_residue_types_.begin(), allowed_residue_types_.end(), *aas_iter ) ==
				allowed_residue_types_.end() /* haven't already added it */ ) {
			allowed_residue_types_.push_back( *aas_iter );
		}
	}

	determine_if_designing();
	determine_if_repacking();
}

ResidueLevelTask::ResidueTypeCOPList const &
ResidueLevelTask_::allowed_residue_types() const
{
	return allowed_residue_types_;
}

ResidueLevelTask::ResidueTypeCOPListConstIter
ResidueLevelTask_::allowed_residue_types_begin() const
{
	return allowed_residue_types_.begin();
}

ResidueLevelTask::ResidueTypeCOPListConstIter
ResidueLevelTask_::allowed_residue_types_end() const
{
	return allowed_residue_types_.end();
}

chemical::ResidueTypeCOP
ResidueLevelTask_::target_type() const {

	bool allowed( std::find( allowed_residue_types_.begin(), allowed_residue_types_.end(),
		target_residue_type_ ) != allowed_residue_types_.end() );
	if ( !allowed ) return 0;
	return target_residue_type_;
}

void
ResidueLevelTask_::print_allowed_types( std::ostream & os ) const
{
	for ( ResidueTypeCOPListConstIter type( allowed_residue_types_begin() );
			type != allowed_residue_types_end(); ++type ) {
		os << '\t' << (**type).name() << '\n';
	}
}

bool ResidueLevelTask_::being_designed() const { return designing_; } // is this residue up for design?
bool ResidueLevelTask_::being_packed() const { return repacking_; } // is this residue being modified at all by the packer

/// @brief ONLY for the RESET command in resfiles: completely reset this position.
/// @details This does several things.  It:
/// - Removes all noncanonicals allowed at this position.
/// - Resets the list of allowed canonicals to the 20 standard canonicals.
/// - Resets the designability of this position (design allowed).
/// - Resets the repacking of this position (repacking allowed).
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void ResidueLevelTask_::reset() {
	using namespace core::chemical;

	mode_tokens_.push_back( "RESET" );
	disallow_noncanonical_aas();

	ResidueTypeSet const & residue_set( *get_original_residue_set() );
	allowed_residue_types_.clear();
	ResidueType const & match_residue_type( residue_set.get_residue_type_with_variant_removed( *original_residue_type_, chemical::VIRTUAL_SIDE_CHAIN ) );
	for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) {
		ResidueTypeCOPs const & aas( residue_set.get_all_types_with_variants_aa( AA( ii ), match_residue_type.variant_types(), pH_mode_exceptions() ) );
		for ( ResidueTypeCOPs::const_iterator
				aas_iter = aas.begin(),
				aas_end = aas.end(); aas_iter != aas_end; ++aas_iter ) {
			allowed_residue_types_.push_back( *aas_iter );
		}
	}

	//Re-enable design and repacking, if they've been disabled:
	repacking_ = true;
	designing_ = true;
	disabled_ = false;
	design_disabled_ = false;
	determine_if_designing(); // This will ensure that disulfides aren't designed away.
	determine_if_repacking();

	return;
}

/// @details  bookkeeping - increases to EX_ONE_STDDEV if boolean is on, but sample level is zero (AS IT SHOULD!)
void ResidueLevelTask_::refresh_ex1_sample_levels()
{
	if ( ex1_ ) {
		if ( ex1_sample_level_ < EX_ONE_STDDEV ) ex1_sample_level_ = EX_ONE_STDDEV;
	}
	if ( ex1aro_ ) {
		if ( ex1aro_sample_level_ < EX_ONE_STDDEV ) ex1aro_sample_level_ = EX_ONE_STDDEV;
	}
	if ( ex1aro_exposed_ ) {
		if ( ex1aro_exposed_sample_level_ < EX_ONE_STDDEV ) ex1aro_exposed_sample_level_ = EX_ONE_STDDEV;
	}

	if ( ex1aro_sample_level_ < ex1_sample_level_ ) {
		ex1aro_sample_level_ = ex1_sample_level_;
	}
	if ( ex1aro_sample_level_ < ex1aro_exposed_sample_level_ ) ex1aro_sample_level_ = ex1aro_exposed_sample_level_;
}

/// @details  bookkeeping - increases to EX_ONE_STDDEV if boolean is on, but sample level is zero
void ResidueLevelTask_::refresh_ex2_sample_levels()
{
	if ( ex2_ ) {
		if ( ex2_sample_level_ < EX_ONE_STDDEV ) ex2_sample_level_ = EX_ONE_STDDEV;
	}
	if ( ex2aro_ ) {
		if ( ex2aro_sample_level_ < EX_ONE_STDDEV ) ex2aro_sample_level_ = EX_ONE_STDDEV;
	}
	if ( ex2aro_exposed_ ) {
		if ( ex2aro_exposed_sample_level_ < EX_ONE_STDDEV ) ex2aro_exposed_sample_level_ = EX_ONE_STDDEV;
	}

	if ( ex2aro_sample_level_ < ex2_sample_level_ ) ex2aro_sample_level_ = ex2_sample_level_;
	if ( ex2aro_sample_level_ < ex2aro_exposed_sample_level_ ) ex2aro_sample_level_ = ex2aro_exposed_sample_level_;
}

/// @details  bookkeeping - increases to EX_ONE_STDDEV if boolean is on, but sample level is zero
void ResidueLevelTask_::refresh_ex3_sample_levels()
{
	if ( ex3_ ) {
		if ( ex3_sample_level_ < EX_ONE_STDDEV ) ex3_sample_level_ = EX_ONE_STDDEV;
	}
}

/// @details  bookkeeping - increases to EX_ONE_STDDEV if boolean is on, but sample level is zero
void ResidueLevelTask_::refresh_ex4_sample_levels()
{
	if ( ex4_ ) {
		if ( ex4_sample_level_ < EX_ONE_STDDEV ) ex4_sample_level_ = EX_ONE_STDDEV;
	}
}


void
ResidueLevelTask_::determine_if_designing()
{
	designing_ = false;
	bool found_aa_difference = false;
	for ( ResidueTypeCOPListConstIter
			allowed_iter = allowed_residue_types_.begin(),
			allowed_end = allowed_residue_types_.end();
			allowed_iter != allowed_end;  ++allowed_iter ) {
		if ( (*allowed_iter)->interchangeability_group() != original_residue_type_->interchangeability_group() ) {
			designing_ = true;
			found_aa_difference = true;
			break;
		}
	}
	if ( design_disabled_ || disabled_ ) {
		runtime_assert( ! found_aa_difference );
		designing_ = false;
	}
}

void ResidueLevelTask_::determine_if_repacking()
{
	if ( ! allowed_residue_types_.empty() ) {
		repacking_ = true;
	} else {
		repacking_ = false;
	}
}

rotamer_set::RotamerOperations const &
ResidueLevelTask_::rotamer_operations() const
{
	return rotamer_operations_;
}

void
ResidueLevelTask_::append_rotamer_operation(
	rotamer_set::RotamerOperationOP rotop
)
{
	rotamer_operations_.push_back( rotop );
}

rotamer_set::RotSetOperationListIterator
ResidueLevelTask_::rotamer_set_operation_begin() const
{
	return rotsetops_.begin();
}

rotamer_set::RotSetOperationListIterator
ResidueLevelTask_::rotamer_set_operation_end() const
{
	return rotsetops_.end();
}

void
ResidueLevelTask_::append_rotamerset_operation(
	rotamer_set::RotamerSetOperationOP rotsetop
)
{
	rotsetops_.push_back( rotsetop );
}

rna::RNA_ResidueLevelTask const &
ResidueLevelTask_::rna_task() const{ return *rna_task_; }

rna::RNA_ResidueLevelTask &
ResidueLevelTask_::nonconst_rna_task() { return *rna_task_; }

void
ResidueLevelTask_::update_union( ResidueLevelTask const & t )
{
	//T << "ResidueLevelTask_::update_union" << std::endl;
	ResidueLevelTask_ const & o(dynamic_cast<ResidueLevelTask_ const &>(t));

	for ( utility::vector1<std::string>::const_iterator i = o.behaviors_.begin(); i != o.behaviors_.end(); ++i ) {
		if ( std::find(behaviors_.begin(),behaviors_.end(),*i) == behaviors_.end() ) {
			//std::cout << "behaviors_ add " << *i;
			behaviors_.push_back(*i);
		}
	}
	for ( rotamer_set::RotamerOperations::const_iterator i = o.rotamer_operations_.begin(); i != o.rotamer_operations_.end(); ++i ) {
		if ( std::find(rotamer_operations_.begin(),rotamer_operations_.end(),*i) == rotamer_operations_.end() ) {
			//std::cout << "rotamer_operations_ add " << *i;
			rotamer_operations_.push_back(*i);
		}
	}
	for ( rotamer_set::RotSetOperationList::const_iterator i = o.rotsetops_.begin(); i != o.rotsetops_.end(); ++i ) {
		if ( std::find(rotsetops_.begin(),rotsetops_.end(),*i) == rotsetops_.end() ) {
			//std::cout << "rotsetops_ add " << *i;
			rotsetops_.push_back(*i);
		}
	}
	for ( std::vector<std::string>::const_iterator i = o.mode_tokens_.begin(); i != o.mode_tokens_.end(); ++i ) {
		if ( std::find(mode_tokens_.begin(),mode_tokens_.end(),*i) == mode_tokens_.end() ) {
			//std::cout << "mode_tokens_ add " << *i;
			mode_tokens_.push_back(*i);
		}
	}
	for ( ResidueTypeCOPList::const_iterator i = o.allowed_residue_types_.begin(); i != o.allowed_residue_types_.end(); ++i ) {
		if ( std::find(allowed_residue_types_.begin(),allowed_residue_types_.end(),*i) == allowed_residue_types_.end() ) {
			//std::cout << "allowed_residue_types_ add " << (*i)->name();
			allowed_residue_types_.push_back(*i);
			//std::cout << " set disabled false ";

		}
	}

	determine_if_designing();
	determine_if_repacking();
	design_disabled_ = !designing_;
	disabled_  = !repacking_;

	include_current_             |= o.include_current_;
	adducts_                     |= o.adducts_;
	optimize_H_mode_             |= o.optimize_H_mode_;
	preserve_c_beta_             |= o.preserve_c_beta_;
	flip_HNQ_                    |= o.flip_HNQ_;
	fix_his_tautomer_            |= o.fix_his_tautomer_;
	include_virtual_side_chain_  |= o.include_virtual_side_chain_;
	sample_proton_chi_           |= o.sample_proton_chi_;
	ex1_                         |= o.ex1_;
	ex2_                         |= o.ex2_;
	ex3_                         |= o.ex3_;
	ex4_                         |= o.ex4_;
	ex1aro_                      |= o.ex1aro_;
	ex2aro_                      |= o.ex2aro_;
	ex1aro_exposed_              |= o.ex1aro_exposed_;
	ex2aro_exposed_              |= o.ex2aro_exposed_;
	operate_on_ex1_              |= o.operate_on_ex1_;
	operate_on_ex2_              |= o.operate_on_ex2_;
	operate_on_ex3_              |= o.operate_on_ex3_;
	operate_on_ex4_              |= o.operate_on_ex4_;
	ex1_sample_level_             = (ExtraRotSample)std::max((int)ex1_sample_level_,(int)o.ex1_sample_level_);
	ex2_sample_level_             = (ExtraRotSample)std::max((int)ex2_sample_level_,(int)o.ex2_sample_level_);
	ex3_sample_level_             = (ExtraRotSample)std::max((int)ex3_sample_level_,(int)o.ex3_sample_level_);
	ex4_sample_level_             = (ExtraRotSample)std::max((int)ex4_sample_level_,(int)o.ex4_sample_level_);
	ex1aro_sample_level_          = (ExtraRotSample)std::max((int)ex1aro_sample_level_,(int)o.ex1aro_sample_level_);
	ex2aro_sample_level_          = (ExtraRotSample)std::max((int)ex2aro_sample_level_,(int)o.ex2aro_sample_level_);
	ex1aro_exposed_sample_level_  = (ExtraRotSample)std::max((int)ex1aro_exposed_sample_level_,(int)o.ex1aro_exposed_sample_level_);
	ex2aro_exposed_sample_level_  = (ExtraRotSample)std::max((int)ex2aro_exposed_sample_level_,(int)o.ex2aro_exposed_sample_level_);
	exdna_sample_level_           = (ExtraRotSample)std::max((int)exdna_sample_level_,(int)o.exdna_sample_level_);
	extrachi_cutoff_              = std::min(o.extrachi_cutoff_,extrachi_cutoff_);

}


void
ResidueLevelTask_::update_intersection( ResidueLevelTask const & t )
{
	//T << "ResidueLevelTask_::update_union" << std::endl;
	ResidueLevelTask_ const & o(dynamic_cast<ResidueLevelTask_ const &>(t));

	utility::vector1<std::string> new_behaviors;
	for ( utility::vector1<std::string>::const_iterator i = o.behaviors_.begin(); i != o.behaviors_.end(); ++i ) {
		if ( std::find(  behaviors_.begin(),  behaviors_.end(),*i) !=   behaviors_.end() &&
				std::find(o.behaviors_.begin(),o.behaviors_.end(),*i) != o.behaviors_.end() ) {
			new_behaviors.push_back(*i);
		}
	}
	behaviors_ = new_behaviors;

	rotamer_set::RotamerOperations new_rotamer_operations;
	for ( rotamer_set::RotamerOperations::const_iterator i = o.rotamer_operations_.begin(); i != o.rotamer_operations_.end(); ++i ) {
		if ( std::find(  rotamer_operations_.begin(),  rotamer_operations_.end(),*i) !=   rotamer_operations_.end() &&
				std::find(o.rotamer_operations_.begin(),o.rotamer_operations_.end(),*i) != o.rotamer_operations_.end() ) {
			new_rotamer_operations.push_back(*i);
		}
	}
	rotamer_operations_ = new_rotamer_operations;


	rotamer_set::RotSetOperationList new_rotsetops;
	for ( rotamer_set::RotSetOperationList::const_iterator i = o.rotsetops_.begin(); i != o.rotsetops_.end(); ++i ) {
		if ( std::find(  rotsetops_.begin(),  rotsetops_.end(),*i) !=   rotsetops_.end() &&
				std::find(o.rotsetops_.begin(),o.rotsetops_.end(),*i) != o.rotsetops_.end() ) {
			new_rotsetops.push_back(*i);
		}
	}
	rotsetops_ = new_rotsetops;

	std::vector<std::string> new_mode_tokens;
	for ( std::vector<std::string>::const_iterator i = o.mode_tokens_.begin(); i != o.mode_tokens_.end(); ++i ) {
		if ( std::find(  mode_tokens_.begin(),  mode_tokens_.end(),*i) !=   mode_tokens_.end() &&
				std::find(o.mode_tokens_.begin(),o.mode_tokens_.end(),*i) != o.mode_tokens_.end() ) {
			new_mode_tokens.push_back(*i);
		}
	}
	mode_tokens_ = new_mode_tokens;

	ResidueTypeCOPList new_allowed_residue_types;
	for ( ResidueTypeCOPList::const_iterator i = o.allowed_residue_types_.begin(); i != o.allowed_residue_types_.end(); ++i ) {
		if ( std::find(  allowed_residue_types_.begin(),  allowed_residue_types_.end(),*i) !=   allowed_residue_types_.end() &&
				std::find(o.allowed_residue_types_.begin(),o.allowed_residue_types_.end(),*i) != o.allowed_residue_types_.end() ) {
			new_allowed_residue_types.push_back(*i);
		}
	}
	allowed_residue_types_ = new_allowed_residue_types;

	determine_if_designing();
	determine_if_repacking();
	design_disabled_ = !designing_;
	disabled_  = !repacking_;

	include_current_             &= o.include_current_;
	adducts_                     &= o.adducts_;
	optimize_H_mode_             &= o.optimize_H_mode_;
	preserve_c_beta_             &= o.preserve_c_beta_;
	flip_HNQ_                    &= o.flip_HNQ_;
	fix_his_tautomer_            &= o.fix_his_tautomer_;
	include_virtual_side_chain_  &= o.include_virtual_side_chain_;
	sample_proton_chi_           &= o.sample_proton_chi_;
	ex1_                         &= o.ex1_;
	ex2_                         &= o.ex2_;
	ex3_                         &= o.ex3_;
	ex4_                         &= o.ex4_;
	ex1aro_                      &= o.ex1aro_;
	ex2aro_                      &= o.ex2aro_;
	ex1aro_exposed_              &= o.ex1aro_exposed_;
	ex2aro_exposed_              &= o.ex2aro_exposed_;
	operate_on_ex1_              &= o.operate_on_ex1_;
	operate_on_ex2_              &= o.operate_on_ex2_;
	operate_on_ex3_              &= o.operate_on_ex3_;
	operate_on_ex4_              &= o.operate_on_ex4_;
	ex1_sample_level_             = (ExtraRotSample)std::min((int)ex1_sample_level_,(int)o.ex1_sample_level_);
	ex2_sample_level_             = (ExtraRotSample)std::min((int)ex2_sample_level_,(int)o.ex2_sample_level_);
	ex3_sample_level_             = (ExtraRotSample)std::min((int)ex3_sample_level_,(int)o.ex3_sample_level_);
	ex4_sample_level_             = (ExtraRotSample)std::min((int)ex4_sample_level_,(int)o.ex4_sample_level_);
	ex1aro_sample_level_          = (ExtraRotSample)std::min((int)ex1aro_sample_level_,(int)o.ex1aro_sample_level_);
	ex2aro_sample_level_          = (ExtraRotSample)std::min((int)ex2aro_sample_level_,(int)o.ex2aro_sample_level_);
	ex1aro_exposed_sample_level_  = (ExtraRotSample)std::min((int)ex1aro_exposed_sample_level_,(int)o.ex1aro_exposed_sample_level_);
	ex2aro_exposed_sample_level_  = (ExtraRotSample)std::min((int)ex2aro_exposed_sample_level_,(int)o.ex2aro_exposed_sample_level_);
	exdna_sample_level_           = (ExtraRotSample)std::min((int)exdna_sample_level_,(int)o.exdna_sample_level_);
	extrachi_cutoff_              = std::max(o.extrachi_cutoff_,extrachi_cutoff_);

}

void
ResidueLevelTask_::update_commutative(
	ResidueLevelTask const & t
){
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	ResidueLevelTask_ const & o(dynamic_cast<ResidueLevelTask_ const &>(t));

	include_current_              |= o.include_current_;
	adducts_                      |= o.adducts_;
	//designing_                    = o.designing_; // need to call determine_if_designing
	//repacking_                    = o.repacking_; // need to call determine_if_repacking
	optimize_H_mode_              |= o.optimize_H_mode_;
	preserve_c_beta_              |= o.preserve_c_beta_;
	flip_HNQ_                     |= o.flip_HNQ_;
	fix_his_tautomer_             |= o.fix_his_tautomer_;
	include_virtual_side_chain_   |= o.include_virtual_side_chain_;
	//disabled_                     |= o.disabled_; // need to call determine_if_repacking()
	//design_disabled_              |= o.design_disabled_; // need to call determine_if_designing()
	sample_proton_chi_            = o.sample_proton_chi_; // <--- apparently sample_proton_chi is not commutatively assigned
	ex1_                          |= o.ex1_;
	ex2_                          |= o.ex2_;
	ex3_                          |= o.ex3_;
	ex4_                          |= o.ex4_;
	ex1aro_                       |= o.ex1aro_;
	ex2aro_                       |= o.ex2aro_;
	ex1aro_exposed_               |= o.ex1aro_exposed_;
	ex2aro_exposed_               |= o.ex2aro_exposed_;
	ex1_sample_level_             = std::max( ex1_sample_level_, o.ex1_sample_level_ );
	ex2_sample_level_             = std::max( ex2_sample_level_, o.ex2_sample_level_ );
	ex3_sample_level_             = std::max( ex3_sample_level_, o.ex3_sample_level_ );
	ex4_sample_level_             = std::max( ex4_sample_level_, o.ex4_sample_level_ );
	ex1aro_sample_level_          = std::max( ex1aro_sample_level_, o.ex1aro_sample_level_ );
	ex2aro_sample_level_          = std::max( ex2aro_sample_level_, o.ex2aro_sample_level_ );
	ex1aro_exposed_sample_level_  = std::max( ex1aro_exposed_sample_level_, o.ex1aro_exposed_sample_level_ );
	ex2aro_exposed_sample_level_  = std::max( ex2aro_exposed_sample_level_, o.ex2aro_exposed_sample_level_ );
	exdna_sample_level_           = std::max( exdna_sample_level_, o.exdna_sample_level_ );
	extrachi_cutoff_              = std::min( extrachi_cutoff_, o.extrachi_cutoff_ );
	operate_on_ex1_               |= o.operate_on_ex1_;
	operate_on_ex2_               |= o.operate_on_ex2_;
	operate_on_ex3_               |= o.operate_on_ex3_;
	operate_on_ex4_               |= o.operate_on_ex4_;

	target_residue_type_ = o.target_residue_type_; // <-- there can be only one, so, this is obviously a non-commutative member

	// merge the allowed residue types to find the set that overlaps
	ResidueTypeCOPList my_allowed_residue_types( allowed_residue_types_ );
	ResidueTypeCOPList o_allowed_residue_types( o.allowed_residue_types_ );
	my_allowed_residue_types.sort();
	o_allowed_residue_types.sort();
	ResidueTypeCOPList common;
	for ( ResidueTypeCOPListConstIter
			myiter = my_allowed_residue_types.begin(),
			myend = my_allowed_residue_types.end(),
			oiter = o_allowed_residue_types.begin(),
			oend = o_allowed_residue_types.end();
			myiter != myend && oiter != oend; /* no increment */ ) {
		//std::cout << " myiter: " << *myiter << " " << (*myiter)->name() << " oiter: " << *oiter << " " << (*oiter)->name() << std::endl;
		if ( *myiter == *oiter ) {
			common.push_back( *myiter );
			//std::cout << "Common! " << (*myiter)->name() << std::endl;
			++myiter;
			++oiter;
		} else if ( *myiter < *oiter ) {
			++myiter;
		} else {
			++oiter;
		}
	}
	// Now insert the elements of common in their original order into allowed_residue_types_
	ResidueTypeCOPList my_allowed_residue_types2;
	my_allowed_residue_types2.swap( allowed_residue_types_ );
	for ( ResidueTypeCOPListConstIter
			myiter = my_allowed_residue_types2.begin(),
			myend = my_allowed_residue_types2.end();
			myiter != myend; ++myiter ) {
		if ( std::find( common.begin(), common.end(), *myiter ) != common.end() ) {
			allowed_residue_types_.push_back( *myiter );
		}
	}
	determine_if_repacking();
	determine_if_designing();
	disabled_                      = !o.repacking_;
	design_disabled_               = !o.designing_;

	/// Form a union of the following sets
	rotamer_operations_.insert(rotamer_operations_.begin(),o.rotamer_operations_.begin(),o.rotamer_operations_.end());
	rotsetops_         .insert(rotsetops_         .begin(),o.rotsetops_         .begin(),o.rotsetops_         .end());
	behaviors_         .insert(behaviors_         .begin(),o.behaviors_         .begin(),o.behaviors_         .end());
	mode_tokens_       .insert(mode_tokens_       .begin(),o.mode_tokens_       .begin(),o.mode_tokens_       .end());
}

} //namespace task
} //namespace pack
} //namespace core


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pack::task::ResidueLevelTask_::save( Archive & arc ) const {
	arc( cereal::base_class< core::pack::task::ResidueLevelTask >( this ) );
	arc( CEREAL_NVP( include_current_ ) ); // _Bool
	arc( CEREAL_NVP( behaviors_ ) ); // utility::vector1<std::string>
	arc( CEREAL_NVP( adducts_ ) ); // _Bool
	arc( CEREAL_NVP( original_residue_type_set_ ) ); // chemical::ResidueTypeSetCOP
	arc( CEREAL_NVP( allowed_residue_types_ ) );
	arc( CEREAL_NVP( original_residue_type_ ) );
	arc( CEREAL_NVP( target_residue_type_ ) );
	arc( CEREAL_NVP( designing_ ) ); // _Bool
	arc( CEREAL_NVP( repacking_ ) ); // _Bool
	arc( CEREAL_NVP( optimize_H_mode_ ) ); // _Bool
	arc( CEREAL_NVP( preserve_c_beta_ ) ); // _Bool
	arc( CEREAL_NVP( flip_HNQ_ ) ); // _Bool
	arc( CEREAL_NVP( fix_his_tautomer_ ) ); // _Bool
	arc( CEREAL_NVP( include_virtual_side_chain_ ) ); // _Bool
	arc( CEREAL_NVP( disabled_ ) ); // _Bool
	arc( CEREAL_NVP( design_disabled_ ) ); // _Bool
	arc( CEREAL_NVP( sample_proton_chi_ ) ); // _Bool
	arc( CEREAL_NVP( ex1_ ) ); // _Bool
	arc( CEREAL_NVP( ex2_ ) ); // _Bool
	arc( CEREAL_NVP( ex3_ ) ); // _Bool
	arc( CEREAL_NVP( ex4_ ) ); // _Bool
	arc( CEREAL_NVP( ex1aro_ ) ); // _Bool
	arc( CEREAL_NVP( ex2aro_ ) ); // _Bool
	arc( CEREAL_NVP( ex1aro_exposed_ ) ); // _Bool
	arc( CEREAL_NVP( ex2aro_exposed_ ) ); // _Bool
	arc( CEREAL_NVP( ex1_sample_level_ ) ); // enum core::pack::task::ExtraRotSample
	arc( CEREAL_NVP( ex2_sample_level_ ) ); // enum core::pack::task::ExtraRotSample
	arc( CEREAL_NVP( ex3_sample_level_ ) ); // enum core::pack::task::ExtraRotSample
	arc( CEREAL_NVP( ex4_sample_level_ ) ); // enum core::pack::task::ExtraRotSample
	arc( CEREAL_NVP( ex1aro_sample_level_ ) ); // enum core::pack::task::ExtraRotSample
	arc( CEREAL_NVP( ex2aro_sample_level_ ) ); // enum core::pack::task::ExtraRotSample
	arc( CEREAL_NVP( ex1aro_exposed_sample_level_ ) ); // enum core::pack::task::ExtraRotSample
	arc( CEREAL_NVP( ex2aro_exposed_sample_level_ ) ); // enum core::pack::task::ExtraRotSample
	arc( CEREAL_NVP( exdna_sample_level_ ) ); // enum core::pack::task::ExtraRotSample
	arc( CEREAL_NVP( extrachi_cutoff_ ) ); // Size
	arc( CEREAL_NVP( operate_on_ex1_ ) ); // _Bool
	arc( CEREAL_NVP( operate_on_ex2_ ) ); // _Bool
	arc( CEREAL_NVP( operate_on_ex3_ ) ); // _Bool
	arc( CEREAL_NVP( operate_on_ex4_ ) ); // _Bool
	arc( CEREAL_NVP( rotamer_operations_ ) ); // rotamer_set::RotamerOperations
	arc( CEREAL_NVP( rotsetops_ ) ); // rotamer_set::RotSetOperationList
	arc( CEREAL_NVP( mode_tokens_ ) ); // std::vector<std::string>
	arc( CEREAL_NVP( rna_task_ ) ); // rna::RNA_ResidueLevelTaskOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pack::task::ResidueLevelTask_::load( Archive & arc ) {
	arc( cereal::base_class< core::pack::task::ResidueLevelTask >( this ) );
	arc( include_current_ ); // _Bool
	arc( behaviors_ ); // utility::vector1<std::string>
	arc( adducts_ ); // _Bool
	arc( original_residue_type_set_ ); // chemical::ResidueTypeSetCOP
	arc( allowed_residue_types_ );
	arc( original_residue_type_ );
	arc( target_residue_type_ );
	arc( designing_ ); // _Bool
	arc( repacking_ ); // _Bool
	arc( optimize_H_mode_ ); // _Bool
	arc( preserve_c_beta_ ); // _Bool
	arc( flip_HNQ_ ); // _Bool
	arc( fix_his_tautomer_ ); // _Bool
	arc( include_virtual_side_chain_ ); // _Bool
	arc( disabled_ ); // _Bool
	arc( design_disabled_ ); // _Bool
	arc( sample_proton_chi_ ); // _Bool
	arc( ex1_ ); // _Bool
	arc( ex2_ ); // _Bool
	arc( ex3_ ); // _Bool
	arc( ex4_ ); // _Bool
	arc( ex1aro_ ); // _Bool
	arc( ex2aro_ ); // _Bool
	arc( ex1aro_exposed_ ); // _Bool
	arc( ex2aro_exposed_ ); // _Bool
	arc( ex1_sample_level_ ); // enum core::pack::task::ExtraRotSample
	arc( ex2_sample_level_ ); // enum core::pack::task::ExtraRotSample
	arc( ex3_sample_level_ ); // enum core::pack::task::ExtraRotSample
	arc( ex4_sample_level_ ); // enum core::pack::task::ExtraRotSample
	arc( ex1aro_sample_level_ ); // enum core::pack::task::ExtraRotSample
	arc( ex2aro_sample_level_ ); // enum core::pack::task::ExtraRotSample
	arc( ex1aro_exposed_sample_level_ ); // enum core::pack::task::ExtraRotSample
	arc( ex2aro_exposed_sample_level_ ); // enum core::pack::task::ExtraRotSample
	arc( exdna_sample_level_ ); // enum core::pack::task::ExtraRotSample
	arc( extrachi_cutoff_ ); // Size
	arc( operate_on_ex1_ ); // _Bool
	arc( operate_on_ex2_ ); // _Bool
	arc( operate_on_ex3_ ); // _Bool
	arc( operate_on_ex4_ ); // _Bool
	arc( rotamer_operations_ ); // rotamer_set::RotamerOperations
	arc( rotsetops_ ); // rotamer_set::RotSetOperationList
	arc( mode_tokens_ ); // std::vector<std::string>
	arc( rna_task_ ); // rna::RNA_ResidueLevelTaskOP
}

SAVE_AND_LOAD_SERIALIZABLE( core::pack::task::ResidueLevelTask_ );
CEREAL_REGISTER_TYPE( core::pack::task::ResidueLevelTask_ )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_task_ResidueLevelTask_ )
#endif // SERIALIZATION
