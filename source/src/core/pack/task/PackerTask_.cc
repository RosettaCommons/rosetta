// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/PackerTask_.hh
/// @brief  Implementation class for task class to describe packer's behavior header
/// Almost all of rosetta needs to use packer tasks, but very little of rosetta needs
/// to see how it behaves internally.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Steven Lewis (smlewi@unc.edu)


//Unit Headers
#include <core/pack/task/PackerTask_.hh>

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
#include <core/chemical/ResidueSelector.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/pose/Pose.hh>
#include <basic/options/option.hh>
#include <core/pose/PDBInfo.hh>
#include <core/id/SequenceMapping.hh>

// Objexx headers
#include <ObjexxFCL/format.hh>

//Utility Headers
#include <basic/Tracer.hh>

#include <iostream>

// option key includes

// AUTO-REMOVED #include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <utility/vector1.hh>

using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

namespace core {
namespace pack {
namespace task {

static thread_local basic::Tracer T( "core.pack.task", basic::t_info );

///@details ResidueLevelTask constructor has following defaults:
///all ex set to false with zero sample level
///position is packable and designable to any of the canonical aa's,
///with variant matching (for termini, etc)
///current rotamer is not included for packer.
///bump check is deactivated by default.
ResidueLevelTask_::ResidueLevelTask_(
	conformation::Residue const & original_residue
)
:
	include_current_( false ),
	adducts_( true ),
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
	// pb -- get the residue_set from the current residue
	ResidueTypeSet const & residue_set( original_residue.residue_type_set() );

	if ( original_residue.is_protein() || original_residue.is_peptoid() ) {
		//default: all amino acids at all positions -- additional "and" operations will remove
		// amino acids from the list of allowed ones
		//no rule yet to treat chemically modified aa's differently
		ResidueType const & match_residue_type( residue_set.get_residue_type_with_variant_removed( original_residue.type(), chemical::VIRTUAL_SIDE_CHAIN ) );
		for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) {
			ResidueTypeCOPs const & aas( residue_set.aa_map( chemical::AA( ii )));
			for ( ResidueTypeCOPs::const_iterator
					aas_iter = aas.begin(),
					aas_end = aas.end(); aas_iter != aas_end; ++aas_iter ) {
				if ( variants_match( match_residue_type, **aas_iter ) ) {
					allowed_residue_types_.push_back( *aas_iter );
				}
			}
		}
		// allow noncanonical AAs and D-amino acids to be repacked
		if (original_residue.aa() == aa_unk || core::chemical::is_canonical_D_aa( original_residue.aa() ) )
			allowed_residue_types_.push_back( original_residue.type().get_self_ptr() );
	} else if ( original_residue.is_DNA() ) {
		// default: all canonical DNA types w/ adducts
	  ResidueTypeCOPs dna_types( ResidueSelector().set_property("DNA").select( residue_set ) );
		for ( ResidueTypeCOPs::const_iterator type( dna_types.begin() ), end( dna_types.end() );
		      type != end; ++type ) {
			if ( nonadduct_variants_match( original_residue.type(), **type ) ) {
				allowed_residue_types_.push_back( *type );
			}
		}
	} else if ( original_residue.is_RNA() ) {
		ResidueType const & match_residue_type(
				residue_set.get_residue_type_with_variant_removed( original_residue.type(), chemical::VIRTUAL_O2PRIME_HYDROGEN ) );
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
			if (buried ) {
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
	using namespace basic::options;
	using namespace OptionKeys::packing;

	initialize_extra_rotamer_flags_from_command_line();

	if ( option[ use_input_sc ] ) or_include_current( true );
	if ( option[ OptionKeys::packing::preserve_c_beta ].user() ) or_preserve_c_beta( true );
	if ( option[ OptionKeys::packing::prevent_repacking ].value() ) prevent_repacking();

}

void
ResidueLevelTask_::initialize_extra_rotamer_flags_from_command_line()
{
	using namespace basic::options;
	using namespace OptionKeys::packing;

	if ( option[ ex1::ex1 ].user() ) { or_ex1( option[ ex1::ex1 ].value() ); }
	if ( option[ ex2::ex2 ].user() ) { or_ex2( option[ ex2::ex2 ].value() ); }
	if ( option[ ex3::ex3 ].user() ) { or_ex3( option[ ex3::ex3 ].value() ); }
	if ( option[ ex4::ex4 ].user() ) { or_ex4( option[ ex4::ex4 ].value() ); }

	if ( option[ ex1aro::ex1aro ].user() )
	{ or_ex1aro( option[ ex1aro::ex1aro ].value() ); }

	if ( option[ ex2aro::ex2aro ].user() )
	{ or_ex2aro( option[ ex2aro::ex2aro ].value() ); }

	if ( option[ ex1aro_exposed::ex1aro_exposed ].user() )
	{ or_ex1aro_exposed( option[ ex1aro_exposed::ex1aro_exposed ].value() ); }

	if ( option[ ex2aro_exposed::ex2aro_exposed ].user() )
	{ or_ex2aro_exposed( option[ ex2aro_exposed::ex2aro_exposed ].value() ); }

	if ( option[ ex1::level ].user() )
	{ or_ex1_sample_level( static_cast<ExtraRotSample>( option[ ex1::level ].value() )); }

	if ( option[ ex2::level ].user() )
	{ or_ex2_sample_level( static_cast<ExtraRotSample>( option[ ex2::level ].value() )); }

	if ( option[ ex3::level ].user() )
	{ or_ex3_sample_level( static_cast<ExtraRotSample>( option[ ex3::level ].value() )); }

	if ( option[ ex4::level ].user() )
	{ or_ex4_sample_level( static_cast<ExtraRotSample>( option[ ex4::level ].value() )); }

	if ( option[ ex1aro::level ].user() )
	{ or_ex1aro_sample_level( static_cast<ExtraRotSample>( option[ ex1aro::level ].value() )); }

	if ( option[ ex2aro::level ].user() )
	{ or_ex2aro_sample_level( static_cast<ExtraRotSample>( option[ ex2aro::level ].value() )); }

	if ( option[ ex1aro_exposed::level ].user() )
	{ or_ex1aro_exposed_sample_level( static_cast<ExtraRotSample>( option[ ex1aro_exposed::level ].value() )); }

	if ( option[ ex2aro_exposed::level ].user() )
	{ or_ex2aro_exposed_sample_level( static_cast<ExtraRotSample>( option[ ex2aro_exposed::level ].value() )); }

	if ( option[ ex1::operate ].user() ) { or_operate_on_ex1( option[ ex1::operate ]() ); }
	if ( option[ ex2::operate ].user() ) { or_operate_on_ex2( option[ ex2::operate ]() ); }
	if ( option[ ex3::operate ].user() ) { or_operate_on_ex3( option[ ex3::operate ]() ); }
	if ( option[ ex4::operate ].user() ) { or_operate_on_ex4( option[ ex4::operate ]() ); }

	// options to control dna rotamer building
	if ( option[ exdna::exdna ].user() && option[ exdna::level ].user() ) {
		if ( option[ exdna::exdna ].value() && option[ exdna::level ].value() ) {
			or_exdna_sample_level( static_cast<ExtraRotSample>( option[ exdna::level ].value() ) );
		}
	}
	else if ( option[ exdna::exdna ].user() && !option[ exdna::level ].user() ) {
		if ( option[ exdna::exdna ].value() ) {
			or_exdna_sample_level( EX_ONE_STDDEV );
		}
	}
	else if ( !option[ exdna::exdna ].user() && option[ exdna::level ].user() ) {
		if ( option[ exdna::level ].value() ) {
			or_exdna_sample_level( static_cast<ExtraRotSample>( option[ exdna::level ].value() ) );
		}
	}

	// extra chi cutoff
	if ( option[ OptionKeys::packing::extrachi_cutoff ].user() ) {
		and_extrachi_cutoff( option[ OptionKeys::packing::extrachi_cutoff ] );
	}
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
	target_type( original_residue_type_->residue_type_set().aa_map( aa ).front() );
}
void ResidueLevelTask_::target_type( std::string name ) {
	target_type( original_residue_type_->residue_type_set().name_map( name ).get_self_ptr() );
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

///@details this function forces a fixed histidine tautomer by removing the alternate tautomer from the ResidueTypesCAPList.  The fix_his_tautomer_ boolean is maintained for reference that this has been done, but the boolean is not the source of the effect.
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

		if ((*allowed_iter)->aa() == chemical::aa_his) { //we only want to look at histidines
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


///@details this function samples a protein conformation with a virtualized side chain.
void ResidueLevelTask_::or_include_virtual_side_chain( bool setting )
{
	include_virtual_side_chain_ |= setting;
}

bool ResidueLevelTask_::include_virtual_side_chain() const
{
	return include_virtual_side_chain_;
}

///@details and operation -- min -- move only toward a lower cutoff for #neighbors w/i 10A that
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


///@details contract (and) the list of available aas for canonical aa's
///if an amino acid is not present (false) in the boolean vector, then do not allow it at this position
///boolean vector is based on the aa enum; see another example with PIKAA.
///The boolean vector is a 20-length vector in alphabetical order by one-letter code.
void
ResidueLevelTask_::restrict_absent_canonical_aas( utility::vector1< bool > const & allowed_aas )
{
	Size num_allowed = std::count( allowed_aas.begin(), allowed_aas.end(), true);
	std::ostringstream aas;
	if (num_allowed == 0){
		mode_tokens_.push_back("EMPTY");
	} else if ( num_allowed == chemical::num_canonical_aas ){
		// this doesn't constrain anything...
	} else if (num_allowed < chemical::num_canonical_aas/2){
		mode_tokens_.push_back("PIKAA");
		for( Size i = 1; i <= chemical::num_canonical_aas; ++i ){
			if( allowed_aas[ i ] ){
				aas << chemical::oneletter_code_from_aa( static_cast<chemical::AA>(i) );
			}
		}
	} else {
		mode_tokens_.push_back("NOTAA");
		for( Size i = 1; i <= chemical::num_canonical_aas; ++i ){
			if( !allowed_aas[ i ] ){
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


///@details contract (and) the list of available nas for canonical na's
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
	for( utility::vector1 < chemical::AA >::const_iterator na = keep_nas.begin(); na != keep_nas.end(); ++na ){
		// illegal for DNA: single-letter codes redundant
		//nas << chemical::oneletter_code_from_aa( *na );
		nas << chemical::name_from_aa( *na );
	}
	mode_tokens_.push_back( "NA" );
	mode_tokens_.push_back( nas.str() );

}

///@details removes all residues from the allowed residue types list, except the one that matches
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

chemical::ResidueTypeSet const & ResidueLevelTask_::get_original_residue_set() const {
  return original_residue_type_->residue_type_set();
}

chemical::AA const & ResidueLevelTask_::get_original_residue() const {
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

	// get ResidueTypeCOPs vector
	chemical::ResidueTypeCOPs const & aas( residue_set.interchangeability_group_map( interchangeability_group ) );

	for ( chemical::ResidueTypeCOPs::const_iterator	aas_iter = aas.begin(), aas_end = aas.end(); aas_iter != aas_end; ++aas_iter ) {
		if ( variants_match( *original_residue_type_, **aas_iter ) &&
				 (*aas_iter)->aa() >= chemical::num_canonical_aas && // cannot be used to add canonical amino acids
				std::find( allowed_residue_types_.begin(), allowed_residue_types_.end(), *aas_iter ) ==	allowed_residue_types_.end() /* haven't already added it */)
		{
			allowed_residue_types_.push_back( *aas_iter );
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
	allow_noncanonical_aa( interchangeability_group, original_residue_type_->residue_type_set() );
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

	chemical::ResidueTypeSet const & residue_set( original_residue_type_->residue_type_set() );
	chemical::ResidueTypeCOPs const & aas( residue_set.aa_map( aa ) );

	for ( chemical::ResidueTypeCOPs::const_iterator
					aas_iter = aas.begin(), aas_end = aas.end(); aas_iter != aas_end; ++aas_iter ) {
		if ( variants_match( *original_residue_type_, **aas_iter ) &&
				 ( std::find( allowed_residue_types_.begin(), allowed_residue_types_.end(), *aas_iter ) ==
					 allowed_residue_types_.end() ) /* haven't already added it */ ) {
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

///@details  bookkeeping - increases to EX_ONE_STDDEV if boolean is on, but sample level is zero (AS IT SHOULD!)
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

///@details  bookkeeping - increases to EX_ONE_STDDEV if boolean is on, but sample level is zero
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

///@details  bookkeeping - increases to EX_ONE_STDDEV if boolean is on, but sample level is zero
void ResidueLevelTask_::refresh_ex3_sample_levels()
{
	if ( ex3_ ) {
		if ( ex3_sample_level_ < EX_ONE_STDDEV ) ex3_sample_level_ = EX_ONE_STDDEV;
	}
}

///@details  bookkeeping - increases to EX_ONE_STDDEV if boolean is on, but sample level is zero
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

	for(utility::vector1<std::string>::const_iterator i = o.behaviors_.begin(); i != o.behaviors_.end(); ++i){
		if( std::find(behaviors_.begin(),behaviors_.end(),*i) == behaviors_.end() ) {
			//std::cout << "behaviors_ add " << *i;
			behaviors_.push_back(*i);
		}
	}
	for(rotamer_set::RotamerOperations::const_iterator i = o.rotamer_operations_.begin(); i != o.rotamer_operations_.end(); ++i){
		if( std::find(rotamer_operations_.begin(),rotamer_operations_.end(),*i) == rotamer_operations_.end() ) {
			//std::cout << "rotamer_operations_ add " << *i;
			rotamer_operations_.push_back(*i);
		}
	}
	for(rotamer_set::RotSetOperationList::const_iterator i = o.rotsetops_.begin(); i != o.rotsetops_.end(); ++i){
		if( std::find(rotsetops_.begin(),rotsetops_.end(),*i) == rotsetops_.end() ) {
			//std::cout << "rotsetops_ add " << *i;
			rotsetops_.push_back(*i);
		}
	}
	for(std::vector<std::string>::const_iterator i = o.mode_tokens_.begin(); i != o.mode_tokens_.end(); ++i){
		if( std::find(mode_tokens_.begin(),mode_tokens_.end(),*i) == mode_tokens_.end() ) {
			//std::cout << "mode_tokens_ add " << *i;
			mode_tokens_.push_back(*i);
		}
	}
	for(ResidueTypeCOPList::const_iterator i = o.allowed_residue_types_.begin(); i != o.allowed_residue_types_.end(); ++i){
		if( std::find(allowed_residue_types_.begin(),allowed_residue_types_.end(),*i) == allowed_residue_types_.end() ) {
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
	for(utility::vector1<std::string>::const_iterator i = o.behaviors_.begin(); i != o.behaviors_.end(); ++i){
		if( std::find(  behaviors_.begin(),  behaviors_.end(),*i) !=   behaviors_.end() &&
			std::find(o.behaviors_.begin(),o.behaviors_.end(),*i) != o.behaviors_.end() ){
			new_behaviors.push_back(*i);
		}
	}
	behaviors_ = new_behaviors;

	rotamer_set::RotamerOperations new_rotamer_operations;
	for(rotamer_set::RotamerOperations::const_iterator i = o.rotamer_operations_.begin(); i != o.rotamer_operations_.end(); ++i){
		if( std::find(  rotamer_operations_.begin(),  rotamer_operations_.end(),*i) !=   rotamer_operations_.end() &&
			std::find(o.rotamer_operations_.begin(),o.rotamer_operations_.end(),*i) != o.rotamer_operations_.end() ){
			new_rotamer_operations.push_back(*i);
		}
	}
	rotamer_operations_ = new_rotamer_operations;



	rotamer_set::RotSetOperationList new_rotsetops;
	for(rotamer_set::RotSetOperationList::const_iterator i = o.rotsetops_.begin(); i != o.rotsetops_.end(); ++i){
		if( std::find(  rotsetops_.begin(),  rotsetops_.end(),*i) !=   rotsetops_.end() &&
			std::find(o.rotsetops_.begin(),o.rotsetops_.end(),*i) != o.rotsetops_.end() ){
			new_rotsetops.push_back(*i);
		}
	}
	rotsetops_ = new_rotsetops;

	std::vector<std::string> new_mode_tokens;
	for(std::vector<std::string>::const_iterator i = o.mode_tokens_.begin(); i != o.mode_tokens_.end(); ++i){
		if( std::find(  mode_tokens_.begin(),  mode_tokens_.end(),*i) !=   mode_tokens_.end() &&
			std::find(o.mode_tokens_.begin(),o.mode_tokens_.end(),*i) != o.mode_tokens_.end() ){
			new_mode_tokens.push_back(*i);
		}
	}
	mode_tokens_ = new_mode_tokens;

	ResidueTypeCOPList new_allowed_residue_types;
	for(ResidueTypeCOPList::const_iterator i = o.allowed_residue_types_.begin(); i != o.allowed_residue_types_.end(); ++i){
		if( std::find(  allowed_residue_types_.begin(),  allowed_residue_types_.end(),*i) !=   allowed_residue_types_.end() &&
			std::find(o.allowed_residue_types_.begin(),o.allowed_residue_types_.end(),*i) != o.allowed_residue_types_.end() ){
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
	disabled_                     |= o.disabled_;
	design_disabled_              |= o.design_disabled_;
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

	/// Form a union of the following sets
	rotamer_operations_.insert(rotamer_operations_.begin(),o.rotamer_operations_.begin(),o.rotamer_operations_.end());
	rotsetops_         .insert(rotsetops_         .begin(),o.rotsetops_         .begin(),o.rotsetops_         .end());
	behaviors_         .insert(behaviors_         .begin(),o.behaviors_         .begin(),o.behaviors_         .end());
	mode_tokens_       .insert(mode_tokens_       .begin(),o.mode_tokens_       .begin(),o.mode_tokens_       .end());
}

void
PackerTask_::update_residue_union(
	Size resid,
	ResidueLevelTask const & t
){
	ResidueLevelTask_ const & o(dynamic_cast<ResidueLevelTask_ const &>(t));
	residue_tasks_[resid].update_union(o);
	pack_residue_[resid] = residue_tasks_[resid].being_packed();
}

void
PackerTask_::update_residue_intersection(
	Size resid,
	ResidueLevelTask const & t
){
	ResidueLevelTask_ const & o(dynamic_cast<ResidueLevelTask_ const &>(t));
	residue_tasks_[resid].update_intersection(o);
	pack_residue_[resid] = residue_tasks_[resid].being_packed();
}

void
PackerTask_::update_residue_commutative(
	Size resid,
	ResidueLevelTask const & t
){
	ResidueLevelTask_ const & o(dynamic_cast<ResidueLevelTask_ const &>(t));
	residue_tasks_[resid].update_commutative(o);
	pack_residue_[resid] = residue_tasks_[resid].being_packed();
}

void
PackerTask_::update_commutative(
	PackerTask const & t
){

	/*for(int i = 0; i < 9; ++i){
		std::cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		std::cerr << "!!!!!!!!!!!!!!  update update_commutative temporarially is copy! !!!!!!!!!!!!!!" << std::endl;
		std::cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		}*/

	PackerTask_ const & o(dynamic_cast<PackerTask_ const &>(t));

	if( nres_ != o.nres_ ) utility_exit_with_message("unmatching nres_");

	/// leave the state of the pack_residue array as true everywhere.  This is an array
	/// used only inside rotamer_trials and rtmin, so this should be a good idea.
	for(Size i = 1; i <= pack_residue_.size(); ++i) pack_residue_[i] = true;

	for(Size i = 1; i <= residue_tasks_.size(); ++i) residue_tasks_[i].update_commutative(o.residue_tasks_[i]);

	n_to_be_packed_ = o.n_to_be_packed_;/// <-- this derived data will need to be updated
	n_to_be_packed_up_to_date_ = false; /// <-- signal that n_to_be_packed will need to be updated

	linmem_ig_ |= o.linmem_ig_;
	if ( linmem_ig_history_size_at_default_ ) {
		linmem_ig_history_size_at_default_ = o.linmem_ig_history_size_at_default_;
		linmem_ig_history_size_ = o.linmem_ig_history_size_;
	} else {
		if ( linmem_ig_history_size_ > o.linmem_ig_history_size_ ) {
			linmem_ig_history_size_ = o.linmem_ig_history_size_;
		}
	}
	lazy_ig_ |= o.lazy_ig_;

	double_lazy_ig_ |= o.double_lazy_ig_;
	dlig_mem_limit_ = std::min( dlig_mem_limit_, o.dlig_mem_limit_ );

	multi_cool_annealer_ |= o.multi_cool_annealer_;
	mca_history_size_ = std::max( mca_history_size_, o.mca_history_size_ );

	optimize_H_ |= o.optimize_H_;

	bump_check_ = o.bump_check_; /// <-- apparently bump check setting is not commutative!
	max_rotbump_energy_ = std::min( max_rotbump_energy_, o.max_rotbump_energy_ );

	rotamer_couplings_ = o.rotamer_couplings_; /// <--- apparently rotamer couplings assignment is not commutative

	rotamer_links_ = o.rotamer_links_; /// <--- apparently rotamer links assignment is not commutative

	low_temp_ = o.low_temp_;
	high_temp_ = o.high_temp_;
	disallow_quench_ = o.disallow_quench_;


	if(o.IG_edge_reweights_) {
		IG_edge_reweights_ = IGEdgeReweightContainerOP( new IGEdgeReweightContainer(nres_) );
		for(utility::vector1< IGEdgeReweighterOP >::const_iterator i = o.IG_edge_reweights_->reweighters_begin();
			i != o.IG_edge_reweights_->reweighters_end(); ++i) {
			IG_edge_reweights_->add_reweighter(*i);
		}
	} else {
		IG_edge_reweights_ = NULL;
	}

}

void
PackerTask_::request_symmetrize_by_intersection(){
	if( symmetry_status_ == NO_SYMMETRIZATION_REQUEST ){
		symmetry_status_ = REQUEST_SYMMETRIZE_BY_INTERSECTION;
	}
	if( symmetry_status_ != REQUEST_SYMMETRIZE_BY_INTERSECTION ) {
		utility_exit_with_message("PackerTask can be symmetrized by union or by intersection but not both.");
	}
}

void
PackerTask_::request_symmetrize_by_union(){
	if( symmetry_status_ == NO_SYMMETRIZATION_REQUEST ){
		symmetry_status_ = REQUEST_SYMMETRIZE_BY_UNION;
	}
	if( symmetry_status_ != REQUEST_SYMMETRIZE_BY_UNION ) {
		utility_exit_with_message("PackerTask can be symmetrized by union or by intersection but not both.");
	}
}

bool
PackerTask_::symmetrize_by_union() const {
	return symmetry_status_ == REQUEST_SYMMETRIZE_BY_UNION;
}

bool
PackerTask_::symmetrize_by_intersection() const {
	return symmetry_status_ == REQUEST_SYMMETRIZE_BY_INTERSECTION;
}

PackerTask_::PackerTask_(
)
:
	nres_(0 ),
	pack_residue_( nres_, true ),
	n_to_be_packed_( nres_ ),
	n_to_be_packed_up_to_date_( true ),
	linmem_ig_( false ),
	linmem_ig_history_size_at_default_( true ),
	linmem_ig_history_size_( 10 ),
	lazy_ig_( false ),
	double_lazy_ig_( false ),
	dlig_mem_limit_( 0 ),
	multi_cool_annealer_( false ),
	mca_history_size_( 1 ),
	optimize_H_( false ),
	bump_check_( true ),
	max_rotbump_energy_( 5.0 ),
	low_temp_( -1.0 ), //default --> let annealer pick
	high_temp_( -1.0 ), //default --> let annealer pick
	disallow_quench_( false ),
	symmetry_status_(NO_SYMMETRIZATION_REQUEST)
{
	IG_edge_reweights_ = NULL; //default stays empty, no reweighting
}

///@details constructor requires a pose.  most settings are in ResidueLevelTask
///nres_ is copied from the pose, all residues are set to be packable by default, and bump_check is true
///the constructor reads NEITHER the command line flags NOR a resfile; this must be done after creation!
PackerTask_::PackerTask_(
	pose::Pose const & pose
)
:
	nres_( pose.total_residue() ),
	pack_residue_( nres_, true ),
	n_to_be_packed_( nres_ ),
	n_to_be_packed_up_to_date_( true ),
	linmem_ig_( false ),
	linmem_ig_history_size_at_default_( true ),
	linmem_ig_history_size_( 10 ),
	lazy_ig_( false ),
	double_lazy_ig_( false ),
	dlig_mem_limit_( 0 ),
	multi_cool_annealer_( false ),
	mca_history_size_( 1 ),
	optimize_H_( false ),
	bump_check_( true ),
	max_rotbump_energy_( 5.0 ),
	low_temp_( -1.0 ), //default --> let annealer pick
	high_temp_( -1.0 ), //default --> let annealer pick
	disallow_quench_( false ),
	symmetry_status_(NO_SYMMETRIZATION_REQUEST)
{
	//create residue-level tasks
	residue_tasks_.reserve( nres_ );
	for ( Size ii = 1; ii <= nres_; ++ii ) {
		residue_tasks_.push_back( ResidueLevelTask_( pose.residue(ii) ));
	}

	IG_edge_reweights_ = NULL; //default stays empty, no reweighting
}

PackerTask_::~PackerTask_() {}

///@details uses compiler-generated copy ctor
PackerTaskOP
PackerTask_::clone() const
{
	// rely on compiler-generated copy ctor
	PackerTaskOP newtask( new PackerTask_( *this ) );
	return newtask;
}

Size PackerTask_::total_residue() const
{
	return nres_;
}

void
PackerTask_::clean_residue_task( conformation::Residue const & original_residue, Size const seqpos)
{
	if ( !pack_residue_[ seqpos ] ){
		--n_to_be_packed_;
		pack_residue_[ seqpos ] = true;
	}
	residue_tasks_[seqpos] = ResidueLevelTask_( original_residue );
}

///@details turn off packing at all positions.
///This does not affect underlying ResidueLevelTasks, but at the moment there is no method for reversing
void PackerTask_::temporarily_fix_everything()
{
	for ( Size ii = 1; ii <= nres_; ++ii )
	{
		pack_residue_[ ii ] = false;
	}
	n_to_be_packed_ = 0;
	n_to_be_packed_up_to_date_ = true;
}

///@details arbitrarily set the packer mutability for a position
///reverse with same function, opposite bool input
void PackerTask_::temporarily_set_pack_residue( int resid, bool setting )
{
	if ( ! setting && pack_residue_[ resid ]) {
		--n_to_be_packed_;
	} else if ( setting && ! pack_residue_[ resid ] ) {
		++n_to_be_packed_;
	}
	pack_residue_[ resid ] = setting;
	if ( ! n_to_be_packed_up_to_date_ ) update_n_to_be_packed();
}

Size PackerTask_::num_to_be_packed() const
{
	if ( ! n_to_be_packed_up_to_date_ ) update_n_to_be_packed();
	return n_to_be_packed_;
}

bool PackerTask_::pack_residue( int resid ) const
{
	//	if ( resid < 0 || resid > pack_residue_.size() || resid > residue_tasks_.size() ){
	//	T("core.pack.task.PackerTask") << "Resid " << resid << " out of range." << std::endl;
	//}
	return pack_residue_[ resid ] && residue_tasks_[ resid ].being_packed();
}


bool PackerTask_::design_residue( int resid ) const
{
	return ( pack_residue_[ resid ] && residue_tasks_[ resid ].being_designed() );
}

bool PackerTask_::design_any() const
{
	for ( Size ii = 1; ii <= nres_; ++ii ) {
		if ( pack_residue_[ ii ] && residue_tasks_[ ii ].being_designed() ) {
			return true;
		}
	}
	return false;
}

void PackerTask_::set_bump_check( bool setting ) {
	bump_check_ = setting;
}
bool PackerTask_::bump_check() const { return bump_check_; }

void PackerTask_::and_max_rotbump_energy( Real setting )
{
	if ( setting < 0 ) return;
	if ( setting < max_rotbump_energy_ ) {
		max_rotbump_energy_ = setting;
	}
}

Real PackerTask_::max_rotbump_energy() const
{
	return max_rotbump_energy_;
}


void PackerTask_::or_include_current( bool setting )
{
	if ( ! setting ) return; // short circuit
	for ( Size ii = 1; ii <= nres_; ++ii ) {
		residue_tasks_[ ii ].or_include_current( setting );
	}
}

void PackerTask_::or_include_current( bool setting, Size resid )
{
	residue_tasks_[ resid ].or_include_current( setting );
}

bool PackerTask_::include_current( Size resid ) const
{
	return residue_tasks_[ resid ].include_current();
}

void PackerTask_::add_behavior( std::string const & behavior )
{
	for ( Size ii = 1; ii <= nres_; ++ii ) {
		residue_tasks_[ ii ].add_behavior( behavior );
	}
}
void PackerTask_::add_behavior( std::string const & behavior, Size resid )
{
	residue_tasks_[ resid ].add_behavior( behavior );
}
bool PackerTask_::has_behavior( std::string const & behavior, Size resid ) const
{
	return residue_tasks_[ resid ].has_behavior( behavior );
}
bool PackerTask_::has_behavior( Size resid ) const
{
	return residue_tasks_[ resid ].has_behavior();
}


chemical::ResidueTypeCOP
PackerTask_::target_type( Size resid ) const
	{ return residue_tasks_[ resid ].target_type(); }

// adducts
void PackerTask_::or_adducts( bool setting )
{
	if ( ! setting ) return; // short circuit
	for ( Size ii = 1; ii <= nres_; ++ii ) {
		residue_tasks_[ ii ].or_adducts( setting );
	}
}
void PackerTask_::or_adducts( bool setting, Size resid ){
	residue_tasks_[ resid ].or_adducts( setting );
}
bool PackerTask_::adducts( Size resid ) const {
	return residue_tasks_[ resid ].adducts();
}


void PackerTask_::or_optimize_h_mode( bool setting )
{
	if ( ! setting ) return; // short circuit

	optimize_H_ = setting;
	if ( linmem_ig_ ) linmem_ig_ = false; // Don't bother with lmig in optH.

	for ( Size ii = 1; ii <= nres_; ++ii ) {
		residue_tasks_[ ii ].or_optimize_h( setting );
	}
}

void PackerTask_::or_preserve_c_beta( bool setting )
{
	if ( ! setting ) return; // short circuit
	for ( Size ii = 1; ii <= nres_; ++ii ) {
		residue_tasks_[ ii ].or_preserve_c_beta( setting );
	}
}

void PackerTask_::or_flip_HNQ( bool setting )
{
	if ( ! setting ) return; // short circuit
	for ( Size ii = 1; ii <= nres_; ++ii ) {
		residue_tasks_[ ii ].or_flip_HNQ( setting );
	}

}

void PackerTask_::or_fix_his_tautomer( utility::vector1<int> const & positions, bool setting )
{
	if ( ! setting ) return;
	if ( positions.size() == 0 ) { // no positions defined, set for all residue_tasks_
		for (Size ii=1; ii<=nres_; ++ii) {
			residue_tasks_[ii].or_fix_his_tautomer( setting );
		}
	} else { // set for specific positions
		for ( Size ii=1; ii<=positions.size(); ++ii ) {
			residue_tasks_[ positions[ii] ].or_fix_his_tautomer( setting );
		}
	}
}

void PackerTask_::or_linmem_ig( bool setting )
{
	linmem_ig_ |= setting;
}

bool PackerTask_::linmem_ig() const
{
	return linmem_ig_;
}

void PackerTask_::decrease_linmem_ig_history_size( Size setting )
{
	if ( linmem_ig_history_size_at_default_ ) {
		linmem_ig_history_size_at_default_ = false;
		linmem_ig_history_size_ = setting;
	} else {
		if ( setting < linmem_ig_history_size_ ) {
			linmem_ig_history_size_ = setting;
		}
	}
}

Size PackerTask_::linmem_ig_history_size() const {
	return linmem_ig_history_size_;
}

void PackerTask_::or_lazy_ig( bool setting )
{
	lazy_ig_ |= setting;
}

bool PackerTask_::lazy_ig() const
{
	return lazy_ig_;
}

void PackerTask_::or_double_lazy_ig( bool setting )
{
	double_lazy_ig_ = setting;
}

bool PackerTask_::double_lazy_ig() const
{
	return double_lazy_ig_;
}

void PackerTask_::decrease_double_lazy_ig_memlimit( Size nbytes_for_rpes )
{
	if ( nbytes_for_rpes != 0 ) {
		if ( nbytes_for_rpes < dlig_mem_limit_ || dlig_mem_limit_ == 0 ) {
			dlig_mem_limit_ = nbytes_for_rpes;
		}
	}
}

Size PackerTask_::double_lazy_ig_memlimit() const
{
	return dlig_mem_limit_;
}

///@brief read only the command line options for extra rotamer building;
PackerTask &
PackerTask_::initialize_extra_rotamer_flags_from_command_line()
{
	for ( Size ii = 1; ii <= nres_; ++ii ) {
		residue_tasks_[ ii ].initialize_extra_rotamer_flags_from_command_line();
	}
	return *this;
}

void PackerTask_::or_multi_cool_annealer( bool setting )
{
	if ( ! (rotamer_couplings_ || rotamer_links_) ) {
		multi_cool_annealer_ |= setting;
	}
}

bool PackerTask_::multi_cool_annealer() const
{
	return multi_cool_annealer_;
}

void PackerTask_::increase_multi_cool_annealer_history_size( Size setting )
{
	if ( setting > mca_history_size_ ) mca_history_size_ = setting;
}

Size PackerTask_::multi_cool_annealer_history_size() const
{
	return mca_history_size_;
}

void
PackerTask_::show( std::ostream & out ) const {

	out << "#Packer_Task" << std::endl;
	out <<                   std::endl;


	out << "resid\tpack?\tdesign?\tallowed_aas" << std::endl;

	for ( Size i=1, it_end = total_residue();
				i <= it_end; ++i){

		out << i;
		out << "\t" << (pack_residue( i ) ? "TRUE" : "FALSE");
		out << "\t" << (design_residue( i ) ? "TRUE" : "FALSE");

		out << "\t";
		for ( ResidueLevelTask::ResidueTypeCOPListConstIter
						allowed_iter(   residue_task( i ).allowed_residue_types_begin() );
						allowed_iter != residue_task( i ).allowed_residue_types_end();
					++allowed_iter ) {
			out << ((allowed_iter == residue_task( i ).allowed_residue_types_begin()) ? "" : ",") <<
					(*allowed_iter)->name();
		}
		out << std::endl;
	}

	//sml pymol-style selection, great for debugging
	if( basic::options::option[ basic::options::OptionKeys::packing::print_pymol_selection ].value() ){
		for ( Size i=1, it_end = total_residue();	i != it_end; ++i)	if( pack_residue(i) ) out << i << "+";
		out << std::endl;
	}
}

void
PackerTask_::show() const {
	show(std::cout);
}

void
PackerTask_::show_residue_task( std::ostream & out, Size resid ) const {
	out << "Residue " << resid << ": " << residue_task(resid).get_original_residue() << "\n";
	out << " \tpack?\t\t\t\t" << (pack_residue(resid) ? "TRUE" : "FALSE") <<
			"\n" << " \tdesign?\t\t\t\t" << (design_residue(resid) ? "TRUE" : "FALSE") <<
			"\n" << " \tinclude current rotamer?\t" << (residue_task(resid).include_current() ? "TRUE" : "FALSE") <<
			"\n" << " \thas protocol-level behavior?\t" << (residue_task(resid).has_behavior() ? "TRUE" : "FALSE") <<
			"\n" << " \tinclude adducts?\t\t" << (residue_task(resid).adducts() ? "TRUE" : "FALSE") <<
			"\n" << " \toptimize H placements?\t\t" << (residue_task(resid).optimize_h() ? "TRUE" : "FALSE") <<
			"\n" << " \tpreserve C beta?\t\t" << (residue_task(resid).preserve_c_beta() ? "TRUE" : "FALSE") <<
			"\n" << " \tflip HIS/ASN/GLN?\t\t" << (residue_task(resid).flip_HNQ() ? "TRUE" : "FALSE")<<
			"\n" << " \tfix HIS tautamer?\t\t" << (residue_task(resid).fix_his_tautomer() ? "TRUE" : "FALSE") <<
			"\n" << " \tsample proton chi?\t\t" << (residue_task(resid).sample_proton_chi() ? "TRUE" : "FALSE") <<
			"\n" << " \textra chi cutoff:\t\t" << residue_task(resid).extrachi_cutoff() <<
			"\n" << " \textra rotamer for chi 1\t\t" << (residue_task(resid).ex1() ? "TRUE" : "FALSE") <<
			"\n" << " \textra rotamer for chi 2\t\t" << (residue_task(resid).ex2() ? "TRUE" : "FALSE") <<
			"\n" << " \textra rotamer for chi 3\t\t" << (residue_task(resid).ex3() ? "TRUE" : "FALSE") <<
			"\n" << " \textra rotamer for chi 4\t\t" << (residue_task(resid).ex4() ? "TRUE" : "FALSE") <<
			"\n" << " \textra rotamer for chi 1 (aromatics)\t\t" << (residue_task(resid).ex1aro() ? "TRUE" : "FALSE") <<
			"\n" << " \textra rotamer for chi 2 (aromatics)\t\t" << (residue_task(resid).ex2aro() ? "TRUE" : "FALSE") <<
			"\n" << " \textra rotamer for chi 1 (exposed aromatics)\t" <<
			(residue_task(resid).ex1aro_exposed() ? "TRUE" : "FALSE") <<
			"\n" << " \textra rotamer for chi 2 (exposed aromatics)\t" <<
			(residue_task(resid).ex2aro_exposed() ? "TRUE" : "FALSE") << "\n";
}

void
PackerTask_::show_residue_task( Size resid ) const {
	show_residue_task(std::cout, resid);
}

void
PackerTask_::show_all_residue_tasks( std::ostream & out ) const {
	for ( Size i=1, it_end = total_residue(); i <= it_end; ++i){
		show_residue_task( out, i );
	}
}

void
PackerTask_::show_all_residue_tasks() const {
	for ( Size i=1, it_end = total_residue(); i <= it_end; ++i){
		show_residue_task( i );
	}
}

PackerTask &
PackerTask_::initialize_from_command_line()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	T << "Packer task: initialize from command line() " << std::endl;

	if ( option[ packing::linmem_ig ].user() && ! optimize_H_ ) {
		or_linmem_ig( true );
		decrease_linmem_ig_history_size( option[ packing::linmem_ig ] );
	}
	if ( option[ packing::lazy_ig ] && ! optimize_H_ ) {
		or_lazy_ig( true );
	}
	if ( option[ packing::double_lazy_ig ] && ! optimize_H_ ) {
		or_double_lazy_ig( true );
	}
	if ( option[ packing::multi_cool_annealer ].user()  ) {
		or_multi_cool_annealer( true );
		increase_multi_cool_annealer_history_size( option[ packing::multi_cool_annealer ] );
	}

	for ( Size ii = 1; ii <= nres_; ++ii ) {
		residue_tasks_[ ii ].initialize_from_command_line();
	}

	if ( option[ packing::fix_his_tautomer ].user() ) {
		utility::vector1< int > positions_to_fix( option[ packing::fix_his_tautomer ] );
		if ( positions_to_fix.size() == 1 && positions_to_fix[ 1 ] == 0 ) {
			utility::vector1< Size > empty;
			or_fix_his_tautomer( empty, true );
		} else {
			utility::vector1< Size > good_positions;
			good_positions.reserve( positions_to_fix.size() );
			for ( Size ii = 1; ii <= positions_to_fix.size(); ++ii ) {
				if ( positions_to_fix[ ii ] > 0 && positions_to_fix[ ii ] <= int(nres_) ) {
					good_positions.push_back( positions_to_fix[ ii ] );
				} else {
					T.Warning << "Ignoring position " << positions_to_fix[ ii ] << " given on the command line for the fix_his_tautomer flag\n";
					if ( positions_to_fix[ii] == 0 ) {
						T.Warning << "The given value of 0 can only be given if there is only a single value given.\n";
						T.Warning << "If multiple values are given, they must all be greater than 0 and less than the total number of residues in the protein" << std::endl;
					} else if ( positions_to_fix[ii] < 0 ) {
						T.Warning << "A negative value was given, but this flag accepts only positive values" << std::endl;
					} else {
						T.Warning << "There are only " << nres_ << " residues in the Pose at the time this flag is read" << std::endl;
					}
				}
			}
			or_fix_his_tautomer( good_positions, true );
		}
	}

	if ( option[ packing::repack_only ].value() ) {
		restrict_to_repacking();
	}

	and_max_rotbump_energy( option[ packing::max_rotbump_energy ] );

	update_n_to_be_packed();
	return *this;
}

///@details vector boolean is based on residue position, disables packing at false positions
///does nothing to true positions.  Cannot turn on packing.
PackerTask &
PackerTask_::restrict_to_residues(
	utility::vector1< bool > const & residues_allowed_to_be_packed
)
{
	for ( Size ii = 1; ii <= nres_; ++ii ) {
		if ( ! residues_allowed_to_be_packed[ ii ] ) {
			residue_tasks_[ ii ].prevent_repacking(); // permenent disabling of the residue
		}
	}
	update_n_to_be_packed();
	return *this;
}

///@details vector boolean is based on residue position, disables packing at false positions
///does nothing to true positions.  Cannot turn on packing.  Will prevent packing at false positions
///if the original residue type has been otherwise disallowed.
PackerTask &
PackerTask_::restrict_to_repacking()
{
	for ( Size ii = 1; ii <= nres_; ++ii ) {
		residue_tasks_[ ii ].restrict_to_repacking();
	}
	update_n_to_be_packed();
	return *this;
}

ResidueLevelTask const &
PackerTask_::residue_task( Size resid ) const
{
	// packer has no need to read residue-level task for a residue that has been disabled
	//This makes the assumption that ONLY the packer is allowed to use this function
	//assert( pack_residue_[ resid ] );
	return residue_tasks_[ resid ];
}

ResidueLevelTask &
PackerTask_::nonconst_residue_task( Size resid )
{
	n_to_be_packed_up_to_date_ = false;
	return residue_tasks_[ resid ];
}

utility::vector1< bool >
PackerTask_::repacking_residues() const
{
	utility::vector1< bool > being_packed( nres_ );
	for ( Size ii = 1; ii <= nres_; ++ii ) { being_packed[ ii ] = residue_tasks_[ ii ].being_packed(); }
	return being_packed;
}

utility::vector1< bool >
PackerTask_::designing_residues() const
{
	utility::vector1< bool > being_designed( nres_ );
	for ( Size ii = 1; ii <= nres_; ++ii ) { being_designed[ ii ] = residue_tasks_[ ii ].being_designed(); }
	return being_designed;
}

///@brief is there at RotamerCouplings object to worry about? (for DNA GC AT pairing, etc)
bool
PackerTask_::rotamer_couplings_exist() const
{
	return ( rotamer_couplings_ != 0 );
}

///@brief const accessor for the RotamerCouplings object
PackerTask::RotamerCouplingsCOP
PackerTask_::rotamer_couplings() const
{
	return rotamer_couplings_;
}

///@brief setter for the RotamerCouplings object
void
PackerTask_::rotamer_couplings( PackerTask::RotamerCouplingsCOP setting )
{
	rotamer_couplings_ = setting;
	if ( setting ) {
		multi_cool_annealer_ = false;
	}
}

///@brief is there at RotamerLinks object to worry about? (for repeat linking
//equivalent positions, etc)
bool
PackerTask_::rotamer_links_exist() const
{
	return ( rotamer_links_ != 0 );
}

///@brief const accessor for the RotamerCouplings object
PackerTask::RotamerLinksCOP
PackerTask_::rotamer_links() const
{
	return rotamer_links_;
}

///@brief setter for the RotamerCouplings object
void
PackerTask_::rotamer_links( PackerTask::RotamerLinksCOP setting )
{
	rotamer_links_ = setting;
	if ( setting ) {
		multi_cool_annealer_ = false;
	}
}

IGEdgeReweightContainerCOP
PackerTask_::IGEdgeReweights() const{
	return IG_edge_reweights_;
}

IGEdgeReweightContainerOP
PackerTask_::set_IGEdgeReweights()
{
	if( !IG_edge_reweights_ ){
		IG_edge_reweights_ = IGEdgeReweightContainerOP( new IGEdgeReweightContainer( nres_ ) );
	}
	return IG_edge_reweights_;

}

void
PackerTask_::append_rotamer_operation(
	rotamer_set::RotamerOperationOP rotop
)
{
	for ( Size ii = 1; ii <= nres_; ++ii ) {
		residue_tasks_[ ii ].append_rotamer_operation( rotop );
	}
}

void
PackerTask_::append_rotamerset_operation(
	rotamer_set::RotamerSetOperationOP rotsetop
)
{
	for ( Size ii = 1; ii <= nres_; ++ii ) {
		residue_tasks_[ ii ].append_rotamerset_operation( rotsetop );
	}
}

void
PackerTask_::update_n_to_be_packed() const
{
	n_to_be_packed_ = 0;
	for ( Size ii = 1; ii <= nres_; ++ii ) {
		if ( pack_residue_[ ii ] && residue_tasks_[ ii ].being_packed() ) {
			++n_to_be_packed_;
		}
	}
	n_to_be_packed_up_to_date_ = true;
}

//some get-set functions for annealer options
void
PackerTask_::low_temp( Real const & low_temp ){
	low_temp_ = low_temp;
}

void
PackerTask_::high_temp( Real const & high_temp ){
	high_temp_ = high_temp;
}

void
PackerTask_::disallow_quench( bool const & disallow_quench ){
	disallow_quench_ = disallow_quench;
}

Real
PackerTask_::low_temp() const { return low_temp_; }

Real
PackerTask_::high_temp() const { return high_temp_; }

bool
PackerTask_::disallow_quench() const { return disallow_quench_;}

std::string ResidueLevelTask_::task_mode() const{
	std::ostringstream modes;
	for( std::vector<std::string>::const_iterator i = mode_tokens_.begin(); i < mode_tokens_.end(); ++i){
		modes << *i << " ";
	}
	return modes.str();
}

std::string ResidueLevelTask_::get_ex_flags(
	Size const chiid,
	Size const exaro_sample_level,
	Size const ex_sample_level
) const {
	std::ostringstream ex_flags;

	if (exaro_sample_level == 1){
		ex_flags << "EX ARO " << chiid;
	} else if (exaro_sample_level > 1) {
		ex_flags << "EX ARO " << chiid << " LEVEL " << exaro_sample_level;
	} else if (ex_sample_level == 1){
		ex_flags << "EX " << chiid;
	} else if (ex_sample_level > 1){
		ex_flags << "EX " << chiid << " LEVEL " << ex_sample_level;
	} else {
		// use the default
	}
	return ex_flags.str();
}

std::string	ResidueLevelTask_::command_string() const{

	std::ostringstream command_string;
	std::string task_name( task_mode() );
	if (task_name != "NATRO")	command_string << " " << task_name;

	std::string ex1 = get_ex_flags( 1, ex1aro_sample_level_, ex1_sample_level_);
	if (ex1 != "")  command_string << " " << ex1;

	std::string ex2 = get_ex_flags( 2, ex2aro_sample_level_, ex2_sample_level_);
	if (ex2 != "") command_string << " " << ex2;

	std::string ex3 = get_ex_flags( 3, 0, ex3_sample_level_);
	if (ex3 != "") command_string << " " << ex3;

	std::string ex4 = get_ex_flags( 4, 0, ex4_sample_level_);
	if (ex4 != "") command_string << " " << ex4;

	if( extrachi_cutoff_ != 18) command_string << " EX_CUTOFF " << extrachi_cutoff_;

	if( include_current()) command_string << " USE_INPUT_SC";

	//SCAN, AUTO, TARGET etc
	for (utility::vector1<std::string>::const_iterator i = behaviors_.begin(); i < behaviors_.end(); ++i){
		if (*i == "TARGET"){
			command_string << " TARGET ";
			if (target_type() == 0){
				command_string << "_";
			} else {
				command_string << target_type()->name();
			}
		} else {
			command_string << " " << *i;
		}
	}

	return command_string.str();
}

std::string PackerTask_::task_string( pose::Pose const & pose ) const{
	std::ostringstream task_string;
	task_string << "start" << std::endl;

	for (Size i = 1; i <= total_residue(); ++i){
		std::string residue_task_string = residue_tasks_[ i ].command_string( );
		if (residue_task_string != "" && residue_task_string != " ") {
			//std::cout << "residue_task_string -> '" << residue_task_string << "'." << std::endl;
            task_string << pose.pdb_info()->number(i);
			task_string << " " << pose.pdb_info()->chain(i);
			task_string << " " << residue_task_string << std::endl;
		}
	}
	return task_string.str();
}


///@brief output highlights of the internal data of PackerTask_: for
///each residue whether it is to be packed, designed and which amino
///acids are allowed.

std::ostream &
operator << ( std::ostream & os, PackerTask const & task )
{
	// (Matt O'Meara) Note that PackerTask is the BASE pure virtual
	// class while PackerTask_ is the DERIVED class.
	// I have put [this function] in this file (PackerTask_.cc)
	// because it exposed the internal data of a PackerTask_ object.
	// Although PackerTask_ is the only derived class of PackerTask, if
	// there were others, the << operator would probably have to be
	// different too.
	//sml If you have a derived child of PackerTask, and this function isn't working for you, create a PackerTask.print() pure virtual, move these contents to PackerTask_.print(), and have this function call PackerTask.print() so it will virtual-lookup the right function.  You'd need to write your own print too of course.  Google for "virtual friend function idiom".

	task.show( os );
	return os;
}

void PackerTask_::remap_residue_level_tasks(
	core::id::SequenceMappingCOP seqmap,
	core::pose::Pose const & pose
){
	utility::vector1< bool > remapped_pack_residue;
	utility::vector1< ResidueLevelTask_ > remapped_residue_tasks;

	core::id::SequenceMapping reverse_seqmap = core::id::SequenceMapping( *seqmap );
	reverse_seqmap.reverse();

	for( Size ii = 1; ii <= pose.total_residue(); ++ii ){
		//core::Size new_pos =  (*seqmap)[ii];
		core::Size old_pos = reverse_seqmap[ii];

  		if( old_pos == 0 ){
        //find insertions and remapped residues old res index is 0
        remapped_residue_tasks.push_back( ResidueLevelTask_( pose.residue( ii ) ));
        remapped_residue_tasks[ii].initialize_from_command_line();
        remapped_pack_residue.push_back( true );
      }else if( old_pos != 0 ){
        // no change
        remapped_residue_tasks.push_back( residue_tasks_[ old_pos ]);
        remapped_pack_residue.push_back( pack_residue_[ old_pos ] );
      }
	}

	residue_tasks_ = remapped_residue_tasks;
	pack_residue_ = remapped_pack_residue;
	nres_ = pose.total_residue();
	update_n_to_be_packed();

}

PackerTask &
PackerTask_::operator=(PackerTask const &){
	utility_exit_with_message("illegal");
	return *this;
}

} //namespace task
} //namespace pack
} //namespace core

