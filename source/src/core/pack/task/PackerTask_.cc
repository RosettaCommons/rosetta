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
#include <core/pack/task/PackerTask_.hh>
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

// C++ headers
#include <iostream>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/options/keys/OptionKeyList.hh>

using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace pack {
namespace task {

static basic::Tracer T( "core.pack.task", basic::t_info );


void
PackerTask_::update_residue_union(
	Size resid,
	ResidueLevelTask const & t
){
	auto const & o(dynamic_cast<ResidueLevelTask_ const &>(t));
	residue_tasks_[resid].update_union(o);
	pack_residue_[resid] = residue_tasks_[resid].being_packed();
}

void
PackerTask_::update_residue_intersection(
	Size resid,
	ResidueLevelTask const & t
){
	auto const & o(dynamic_cast<ResidueLevelTask_ const &>(t));
	residue_tasks_[resid].update_intersection(o);
	pack_residue_[resid] = residue_tasks_[resid].being_packed();
}

void
PackerTask_::update_residue_commutative(
	Size resid,
	ResidueLevelTask const & t
){
	auto const & o(dynamic_cast<ResidueLevelTask_ const &>(t));
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

	auto const & o(dynamic_cast<PackerTask_ const &>(t));

	if ( nres_ != o.nres_ ) utility_exit_with_message("unmatching nres_");

	/// leave the state of the pack_residue array as true everywhere.  This is an array
	/// used only inside rotamer_trials and rtmin, so this should be a good idea.
	for ( Size i = 1; i <= pack_residue_.size(); ++i ) pack_residue_[i] = true;

	for ( Size i = 1; i <= residue_tasks_.size(); ++i ) residue_tasks_[i].update_commutative(o.residue_tasks_[i]);

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


	if ( o.IG_edge_reweights_ ) {
		IG_edge_reweights_ = IGEdgeReweightContainerOP( new IGEdgeReweightContainer(nres_) );
		for ( auto i = o.IG_edge_reweights_->reweighters_begin();
				i != o.IG_edge_reweights_->reweighters_end(); ++i ) {
			IG_edge_reweights_->add_reweighter(*i);
		}
	} else {
		IG_edge_reweights_ = nullptr;
	}

}

void
PackerTask_::request_symmetrize_by_intersection(){
	if ( symmetry_status_ == NO_SYMMETRIZATION_REQUEST ) {
		symmetry_status_ = REQUEST_SYMMETRIZE_BY_INTERSECTION;
	}
	if ( symmetry_status_ != REQUEST_SYMMETRIZE_BY_INTERSECTION ) {
		utility_exit_with_message("PackerTask can be symmetrized by union or by intersection but not both.");
	}
}

void
PackerTask_::request_symmetrize_by_union(){
	if ( symmetry_status_ == NO_SYMMETRIZATION_REQUEST ) {
		symmetry_status_ = REQUEST_SYMMETRIZE_BY_UNION;
	}
	if ( symmetry_status_ != REQUEST_SYMMETRIZE_BY_UNION ) {
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
	IG_edge_reweights_ = nullptr; //default stays empty, no reweighting
}

/// @details constructor requires a pose.  most settings are in ResidueLevelTask
///nres_ is copied from the pose, all residues are set to be packable by default, and bump_check is true
///the constructor reads NEITHER the command line flags NOR a resfile; this must be done after creation!
PackerTask_::PackerTask_(
	pose::Pose const & pose
)
:
	nres_( pose.size() ),
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
		residue_tasks_.push_back( ResidueLevelTask_( pose.residue(ii), pose ));
	}

	IG_edge_reweights_ = nullptr; //default stays empty, no reweighting
}

PackerTask_::~PackerTask_() = default;

/// @details uses compiler-generated copy ctor
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
PackerTask_::clean_residue_task( conformation::Residue const & original_residue, Size const seqpos, core::pose::Pose const & pose )
{
	if ( !pack_residue_[ seqpos ] ) {
		--n_to_be_packed_;
		pack_residue_[ seqpos ] = true;
	}
	residue_tasks_[seqpos] = ResidueLevelTask_( original_residue, pose );
}

/// @details turn off packing at all positions.
///This does not affect underlying ResidueLevelTasks, but at the moment there is no method for reversing
void PackerTask_::temporarily_fix_everything()
{
	for ( Size ii = 1; ii <= nres_; ++ii ) {
		pack_residue_[ ii ] = false;
	}
	n_to_be_packed_ = 0;
	n_to_be_packed_up_to_date_ = true;
}

/// @details arbitrarily set the packer mutability for a position
///reverse with same function, opposite bool input
void PackerTask_::temporarily_set_pack_residue( int resid, bool setting )
{
	if ( ! setting && pack_residue_[ resid ] ) {
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
	// if ( resid < 0 || resid > pack_residue_.size() || resid > residue_tasks_.size() ){
	// T("core.pack.task.PackerTask") << "Resid " << resid << " out of range." << std::endl;
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
		for ( Size ii=1; ii<=nres_; ++ii ) {
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

void PackerTask_::and_linmem_ig( bool setting )
{
	linmem_ig_ &= setting;
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

/// @brief read only the command line options for extra rotamer building;
PackerTask &
PackerTask_::initialize_extra_rotamer_flags_from_command_line()
{
	return initialize_extra_rotamer_flags_from_options( basic::options::option );
}

/// @brief read only the command line options for extra rotamer building;
PackerTask &
PackerTask_::initialize_extra_rotamer_flags_from_options(
	utility::options::OptionCollection const & options
)
{
	for ( Size ii = 1; ii <= nres_; ++ii ) {
		residue_tasks_[ ii ].initialize_extra_rotamer_flags_from_options( options );
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
			i <= it_end; ++i ) {

		out << i;
		out << "\t" << (pack_residue( i ) ? "TRUE" : "FALSE");
		out << "\t" << (design_residue( i ) ? "TRUE" : "FALSE");

		out << "\t";
		for ( auto
				allowed_iter(   residue_task( i ).allowed_residue_types_begin() );
				allowed_iter != residue_task( i ).allowed_residue_types_end();
				++allowed_iter ) {
			out << ((allowed_iter == residue_task( i ).allowed_residue_types_begin()) ? "" : ",") <<
				(*allowed_iter)->name();
		}
		out << std::endl;
	}

	//sml pymol-style selection, great for debugging
	if ( basic::options::option[ basic::options::OptionKeys::packing::print_pymol_selection ].value() ) {
		for ( Size i=1, it_end = total_residue(); i != it_end; ++i ) if ( pack_residue(i) ) out << i << "+";
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
	for ( Size i=1, it_end = total_residue(); i <= it_end; ++i ) {
		show_residue_task( out, i );
	}
}

void
PackerTask_::show_all_residue_tasks() const {
	for ( Size i=1, it_end = total_residue(); i <= it_end; ++i ) {
		show_residue_task( i );
	}
}

PackerTask &
PackerTask_::initialize_from_command_line()
{
	return initialize_from_options( basic::options::option );
}

PackerTask &
PackerTask_::initialize_from_options( utility::options::OptionCollection const & options )
{
	//using namespace basic::options;
	using namespace basic::options::OptionKeys;

	T << "Packer task: initialize from command line() " << std::endl;

	if ( options[ packing::linmem_ig ].user() && ! optimize_H_ ) {
		or_linmem_ig( true );
		decrease_linmem_ig_history_size( options[ packing::linmem_ig ] );
	}
	if ( options[ packing::lazy_ig ] && ! optimize_H_ ) {
		or_lazy_ig( true );
	}
	if ( options[ packing::double_lazy_ig ] && ! optimize_H_ ) {
		or_double_lazy_ig( true );
	}
	if ( options[ packing::multi_cool_annealer ].user()  ) {
		or_multi_cool_annealer( true );
		increase_multi_cool_annealer_history_size( options[ packing::multi_cool_annealer ] );
	}

	for ( Size ii = 1; ii <= nres_; ++ii ) {
		residue_tasks_[ ii ].initialize_from_options( options );
	}

	if ( options[ packing::fix_his_tautomer ].user() ) {
		utility::vector1< int > positions_to_fix( options[ packing::fix_his_tautomer ] );
		if ( positions_to_fix.size() == 1 && positions_to_fix[ 1 ] == 0 ) {
			utility::vector1< Size > empty;
			or_fix_his_tautomer( empty, true );
		} else {
			utility::vector1< Size > good_positions;
			good_positions.reserve( positions_to_fix.size() );
			for ( auto const position_to_fix : positions_to_fix ) {
				if ( position_to_fix > 0 && position_to_fix <= int(nres_) ) {
					good_positions.push_back( position_to_fix );
				} else {
					T.Warning << "Ignoring position " << position_to_fix << " given on the command line for the fix_his_tautomer flag\n";
					if ( position_to_fix == 0 ) {
						T.Warning << "The given value of 0 can only be given if there is only a single value given.\n";
						T.Warning << "If multiple values are given, they must all be greater than 0 and less than the total number of residues in the protein" << std::endl;
					} else if ( position_to_fix < 0 ) {
						T.Warning << "A negative value was given, but this flag accepts only positive values" << std::endl;
					} else {
						T.Warning << "There are only " << nres_ << " residues in the Pose at the time this flag is read" << std::endl;
					}
				}
			}
			or_fix_his_tautomer( good_positions, true );
		}
	}

	if ( options[ packing::repack_only ].value() ) {
		restrict_to_repacking();
	}

	and_max_rotbump_energy( options[ packing::max_rotbump_energy ] );

	update_n_to_be_packed();
	return *this;
}


void
PackerTask::list_options_read( utility::options::OptionKeyList & read_options )
{
	using namespace basic::options::OptionKeys;
	read_options
		+ packing::linmem_ig
		+ packing::lazy_ig
		+ packing::double_lazy_ig
		+ packing::multi_cool_annealer
		+ packing::fix_his_tautomer
		+ packing::repack_only
		+ packing::max_rotbump_energy;

	ResidueLevelTask::list_options_read( read_options );
}


/// @details vector boolean is based on residue position, disables packing at false positions
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

/// @details vector boolean is based on residue position, disables packing at false positions
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
	//debug_assert( pack_residue_[ resid ] );
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

/// @brief is there at RotamerCouplings object to worry about? (for DNA GC AT pairing, etc)
bool
PackerTask_::rotamer_couplings_exist() const
{
	return ( rotamer_couplings_ != nullptr );
}

/// @brief const accessor for the RotamerCouplings object
PackerTask::RotamerCouplingsCOP
PackerTask_::rotamer_couplings() const
{
	return rotamer_couplings_;
}

/// @brief setter for the RotamerCouplings object
void
PackerTask_::rotamer_couplings( PackerTask::RotamerCouplingsCOP setting )
{
	rotamer_couplings_ = setting;
	if ( setting ) {
		multi_cool_annealer_ = false;
	}
}

/// @brief is there at RotamerLinks object to worry about? (for repeat linking
//equivalent positions, etc)
bool
PackerTask_::rotamer_links_exist() const
{
	return ( rotamer_links_ != nullptr );
}

/// @brief const accessor for the RotamerCouplings object
PackerTask::RotamerLinksCOP
PackerTask_::rotamer_links() const
{
	return rotamer_links_;
}

/// @brief setter for the RotamerCouplings object
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
	if ( !IG_edge_reweights_ ) {
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
	for ( auto i = mode_tokens_.begin(); i < mode_tokens_.end(); ++i ) {
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

	if ( exaro_sample_level == 1 ) {
		ex_flags << "EX ARO " << chiid;
	} else if ( exaro_sample_level > 1 ) {
		ex_flags << "EX ARO " << chiid << " LEVEL " << exaro_sample_level;
	} else if ( ex_sample_level == 1 ) {
		ex_flags << "EX " << chiid;
	} else if ( ex_sample_level > 1 ) {
		ex_flags << "EX " << chiid << " LEVEL " << ex_sample_level;
	} else {
		// use the default
	}
	return ex_flags.str();
}

std::string ResidueLevelTask_::command_string() const{

	std::ostringstream command_string;
	std::string task_name( task_mode() );
	if ( task_name != "NATRO" ) command_string << " " << task_name;

	std::string ex1 = get_ex_flags( 1, ex1aro_sample_level_, ex1_sample_level_);
	if ( ex1 != "" )  command_string << " " << ex1;

	std::string ex2 = get_ex_flags( 2, ex2aro_sample_level_, ex2_sample_level_);
	if ( ex2 != "" ) command_string << " " << ex2;

	std::string ex3 = get_ex_flags( 3, 0, ex3_sample_level_);
	if ( ex3 != "" ) command_string << " " << ex3;

	std::string ex4 = get_ex_flags( 4, 0, ex4_sample_level_);
	if ( ex4 != "" ) command_string << " " << ex4;

	if ( extrachi_cutoff_ != 18 ) command_string << " EX_CUTOFF " << extrachi_cutoff_;

	if ( include_current() ) command_string << " USE_INPUT_SC";

	//SCAN, AUTO, TARGET etc
	for ( std::string const & elem : behaviors_ ) {
		if ( elem == "TARGET" ) {
			command_string << " TARGET ";
			if ( target_type() == nullptr ) {
				command_string << "_";
			} else {
				command_string << target_type()->name();
			}
		} else {
			command_string << " " << elem;
		}
	}

	return command_string.str();
}

std::string PackerTask_::task_string( pose::Pose const & pose ) const{
	std::ostringstream task_string;
	task_string << "start" << std::endl;

	for ( Size i = 1; i <= total_residue(); ++i ) {
		std::string residue_task_string = residue_tasks_[ i ].command_string( );
		if ( residue_task_string != "" && residue_task_string != " " ) {
			//std::cout << "residue_task_string -> '" << residue_task_string << "'." << std::endl;
			task_string << pose.pdb_info()->number(i);
			task_string << " " << pose.pdb_info()->chain(i);
			task_string << " " << residue_task_string << std::endl;
		}
	}
	return task_string.str();
}


/// @brief output highlights of the internal data of PackerTask_: for
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

	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		//core::Size new_pos =  (*seqmap)[ii];
		core::Size old_pos = reverse_seqmap[ii];

		if ( old_pos == 0 ) {
			//find insertions and remapped residues old res index is 0
			remapped_residue_tasks.push_back( ResidueLevelTask_( pose.residue( ii ), pose ));
			remapped_residue_tasks[ii].initialize_from_command_line();
			remapped_pack_residue.push_back( true );
		} else if ( old_pos != 0 ) {
			// no change
			remapped_residue_tasks.push_back( residue_tasks_[ old_pos ]);
			remapped_pack_residue.push_back( pack_residue_[ old_pos ] );
		}
	}

	residue_tasks_ = remapped_residue_tasks;
	pack_residue_ = remapped_pack_residue;
	nres_ = pose.size();
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


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pack::task::PackerTask_::save( Archive & arc ) const {
	arc( cereal::base_class< core::pack::task::PackerTask >( this ) );
	arc( CEREAL_NVP( nres_ ) ); // Size
	arc( CEREAL_NVP( pack_residue_ ) ); // utility::vector1<_Bool>
	arc( CEREAL_NVP( residue_tasks_ ) ); // utility::vector1<ResidueLevelTask_>
	arc( CEREAL_NVP( n_to_be_packed_ ) ); // Size
	arc( CEREAL_NVP( n_to_be_packed_up_to_date_ ) ); // _Bool
	arc( CEREAL_NVP( linmem_ig_ ) ); // _Bool
	arc( CEREAL_NVP( linmem_ig_history_size_at_default_ ) ); // _Bool
	arc( CEREAL_NVP( linmem_ig_history_size_ ) ); // Size
	arc( CEREAL_NVP( lazy_ig_ ) ); // _Bool
	arc( CEREAL_NVP( double_lazy_ig_ ) ); // _Bool
	arc( CEREAL_NVP( dlig_mem_limit_ ) ); // Size
	arc( CEREAL_NVP( multi_cool_annealer_ ) ); // _Bool
	arc( CEREAL_NVP( mca_history_size_ ) ); // Size
	arc( CEREAL_NVP( optimize_H_ ) ); // _Bool
	arc( CEREAL_NVP( bump_check_ ) ); // _Bool
	arc( CEREAL_NVP( max_rotbump_energy_ ) ); // Real
	arc( CEREAL_NVP( rotamer_couplings_ ) ); // RotamerCouplingsCOP
	arc( CEREAL_NVP( rotamer_links_ ) ); // RotamerLinksCOP
	arc( CEREAL_NVP( low_temp_ ) ); // Real
	arc( CEREAL_NVP( high_temp_ ) ); // Real
	arc( CEREAL_NVP( disallow_quench_ ) ); // _Bool
	arc( CEREAL_NVP( IG_edge_reweights_ ) ); // IGEdgeReweightContainerOP
	arc( CEREAL_NVP( symmetry_status_ ) ); // enum core::pack::task::PackerTaskSymmetryStatus
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pack::task::PackerTask_::load( Archive & arc ) {
	arc( cereal::base_class< core::pack::task::PackerTask >( this ) );
	arc( nres_ ); // Size
	arc( pack_residue_ ); // utility::vector1<_Bool>
	arc( residue_tasks_ ); // utility::vector1<ResidueLevelTask_>
	arc( n_to_be_packed_ ); // Size
	arc( n_to_be_packed_up_to_date_ ); // _Bool
	arc( linmem_ig_ ); // _Bool
	arc( linmem_ig_history_size_at_default_ ); // _Bool
	arc( linmem_ig_history_size_ ); // Size
	arc( lazy_ig_ ); // _Bool
	arc( double_lazy_ig_ ); // _Bool
	arc( dlig_mem_limit_ ); // Size
	arc( multi_cool_annealer_ ); // _Bool
	arc( mca_history_size_ ); // Size
	arc( optimize_H_ ); // _Bool
	arc( bump_check_ ); // _Bool
	arc( max_rotbump_energy_ ); // Real
	std::shared_ptr< core::pack::rotamer_set::RotamerCouplings > local_rotamer_couplings;
	arc( local_rotamer_couplings ); // RotamerCouplingsCOP
	rotamer_couplings_ = local_rotamer_couplings; // copy the non-const pointer(s) into the const pointer(s)
	std::shared_ptr< core::pack::rotamer_set::RotamerLinks > local_rotamer_links;
	arc( local_rotamer_links ); // RotamerLinksCOP
	rotamer_links_ = local_rotamer_links; // copy the non-const pointer(s) into the const pointer(s)
	arc( low_temp_ ); // Real
	arc( high_temp_ ); // Real
	arc( disallow_quench_ ); // _Bool
	arc( IG_edge_reweights_ ); // IGEdgeReweightContainerOP
	arc( symmetry_status_ ); // enum core::pack::task::PackerTaskSymmetryStatus
}

SAVE_AND_LOAD_SERIALIZABLE( core::pack::task::PackerTask_ );
CEREAL_REGISTER_TYPE( core::pack::task::PackerTask_ )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_task_PackerTask_ )
#endif // SERIALIZATION
