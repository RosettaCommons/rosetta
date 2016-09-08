// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/swa/FullModelInfo.cc
/// @brief  Mapping from a working pose into a bigger pose, for swa monte carlo stuff.
/// @author Rhiju Das

// Unit headers
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelInfo.fwd.hh>

// Package headers
#include <core/pose/full_model_info/util.hh>
#include <core/pose/full_model_info/FullModelParameters.hh>
#include <core/pose/full_model_info/FullModelParameterType.hh>
#include <core/pose/full_model_info/SubMotifInfo.hh>


// Project headers
#include <core/conformation/Residue.hh>
#include <core/chemical/VariantType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/full_model_info/FullModelParameters.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <basic/datacache/BasicDataCache.hh>

// Utility headers
#include <utility/vector1.hh>

// C++
#include <string>
#include <map>

/////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This object holds information on the 'full model' during stepwise modeling.
//
//  It must be held in the pose, because it is needed in scoring, e.g. to figure out
//   what loops haven't been created yet (necessary to estimate loop_close costs), where there
//   are free groups (like virtualized phosphates) that should be given bonuses,
//   and whether there are other 'sister' poses that need scoring.
//
//  It is also used in stepwise monte carlo to determine what moves are allowed next for the pose.
//
//  -- Rhiju, 2014
//
// A little more detail:
//
//  There is a FullModelParameters object that includes the full sequence that will be modeled,
//   conventional chain/numbering that goes with that full sequence, cutpoints, and any residues
//   that are supposed to have properties like syn chi torsions (for RNA). FullModelParameters is not
//   supposed to change during the run and may be shared between sister poses; that's why its a
//   constant owning pointer and its not explicitly cloned.
//
//  There are also two variables that *do* change during stepwise monte carlo runs:
//
//    res_list
//       which gives the numbering of the pose residues in the context of the full
//       sequence. Note that this numbering is such that the full_sequence has numbering 1,2,...
//       [If you want to specify that the final PDBs uses a different convention, set conventional_numbering_]
//
//    other_pose_list
//       which is a list of pointers to sister poses. Note that the sister poses should not carry lists back to
//       this pose -- I originally had things tha way and started to have major issues
//       in cloning during monte carlo. (So they're formally not really sisters, they're daughters.)
//       In principle, the current framework could also allow the easy setup of 'pose trees' where each model
//       is subdivided into modules that are independent aside from connecting uninstantiated loops, an appealing
//       idea for some RNA structure applications, but I haven't followed this up.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////

using namespace core;
using namespace core::pose::datacache;
using namespace basic::datacache;

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace pose {
namespace full_model_info {


//Constructor
FullModelInfo::FullModelInfo(
	std::string const & full_sequence,
	utility::vector1< Size > const & cutpoint_open_in_full_model,
	utility::vector1< Size > const & res_numbers_in_pose ):
	CacheableData(),
	res_list_( res_numbers_in_pose ),
	full_model_parameters_( FullModelParametersCOP( FullModelParametersOP( new FullModelParameters( full_sequence, cutpoint_open_in_full_model, res_numbers_in_pose ) ) )  )
{
}

//Constructor
FullModelInfo::FullModelInfo( FullModelParametersCOP full_model_parameters ):
	CacheableData(),
	full_model_parameters_( full_model_parameters->clone() )
{
}

//Constructor
FullModelInfo::FullModelInfo( pose::Pose & pose ) :
	CacheableData()
{
	full_model_parameters_ = FullModelParametersCOP( FullModelParametersOP( new FullModelParameters( pose, res_list_ ) ) );
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details Copy constructors must copy all data, not just some...
FullModelInfo::FullModelInfo( FullModelInfo const & src ) :
	CacheableData(),
	res_list_( src.res_list_ ),
	full_model_parameters_( src.full_model_parameters_ ),
	submotif_info_list_( src.submotif_info_list_ )
{
	// we have to tell our daughters in the pose tree that its time to get cloned.
	for ( Size n = 1; n <= src.other_pose_list_.size(); n++ ) {
		other_pose_list_.push_back( src.other_pose_list_[ n ]->clone() );
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
FullModelInfo::~FullModelInfo(){}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// properties of full model.
FullModelParametersCOP
FullModelInfo::full_model_parameters() const {
	return full_model_parameters_;
}
void
FullModelInfo::set_full_model_parameters( FullModelParametersCOP setting ){
	full_model_parameters_ = setting;
}
std::string const &
FullModelInfo::full_sequence() const {
	return full_model_parameters_->full_sequence();
}
utility::vector1< int > const &
FullModelInfo::conventional_numbering() const {
	return full_model_parameters_->conventional_numbering();
}
utility::vector1< char > const &
FullModelInfo::conventional_chains() const {
	return full_model_parameters_->conventional_chains();
}
utility::vector1< Size > const &
FullModelInfo::cutpoint_open_in_full_model() const {
	return full_model_parameters_->get_res_list( CUTPOINT_OPEN );
}
utility::vector1< Size > const &
FullModelInfo::dock_domain_map() const {
	return full_model_parameters_->get_parameter( DOCK_DOMAIN );
}
utility::vector1< Size > const &
FullModelInfo::fixed_domain_map() const {
	return full_model_parameters_->get_parameter( FIXED_DOMAIN );
}
utility::vector1< Size > const &
FullModelInfo::input_domain_map() const {
	return full_model_parameters_->get_parameter( INPUT_DOMAIN );
}
utility::vector1< Size > const &
FullModelInfo::extra_minimize_res() const {
	return full_model_parameters_->get_res_list( EXTRA_MINIMIZE );
}
utility::vector1< Size > const &
FullModelInfo::preferred_root_res() const {
	return full_model_parameters_->get_res_list( PREFERRED_ROOT );
}
utility::vector1< Size > const &
FullModelInfo::sample_res() const {
	return full_model_parameters_->get_res_list( SAMPLE );
}
utility::vector1< Size > const &
FullModelInfo::working_res() const {
	return full_model_parameters_->get_res_list( WORKING );
}
utility::vector1< Size > const &
FullModelInfo::calc_rms_res() const {
	return full_model_parameters_->get_res_list( CALC_RMS );
}
utility::vector1< Size > const &
FullModelInfo::rna_terminal_res() const {
	return full_model_parameters_->get_res_list( RNA_TERMINAL );
}
utility::vector1< Size > const &
FullModelInfo::rna_block_stack_above_res() const {
	return full_model_parameters_->get_res_list( RNA_BLOCK_STACK_ABOVE );
}
utility::vector1< Size > const &
FullModelInfo::rna_block_stack_below_res() const {
	return full_model_parameters_->get_res_list( RNA_BLOCK_STACK_BELOW );
}
utility::vector1< Size > const &
FullModelInfo::rna_bulge_res() const {
	return full_model_parameters_->get_res_list( RNA_BULGE );
}
utility::vector1< Size > const &
FullModelInfo::rna_sample_sugar_res() const {
	return full_model_parameters_->get_res_list( RNA_SAMPLE_SUGAR );
}
utility::vector1< Size > const &
FullModelInfo::rna_syn_chi_res() const {
	return full_model_parameters_->get_res_list( RNA_SYN_CHI );
}
utility::vector1< Size > const &
FullModelInfo::rna_anti_chi_res() const {
	return full_model_parameters_->get_res_list( RNA_ANTI_CHI );
}
utility::vector1< std::pair< Size, Size > >
FullModelInfo::jump_pairs() const {
	return full_model_parameters_->get_res_list_as_pairs( JUMP );
}
utility::vector1< Size > const &
FullModelInfo::jump_res_map() const {
	return full_model_parameters_->get_parameter( JUMP );
}
Size
FullModelInfo::size() const {
	return full_model_parameters_->size();
}
utility::vector1< Size >
FullModelInfo::chains_in_full_model() const {
	return full_model_parameters_->chains_in_full_model();
}

///////////////////////////////////////////////////////////////////////////////////////
std::map< Size, Size >
FullModelInfo::full_to_sub() const {
	std::map< Size, Size > full_to_sub;
	for ( Size n = 1; n <= size(); n++ ) {
		if ( res_list_.has_value( n ) ) full_to_sub[ n ] = res_list_.index( n );
	}
	return full_to_sub;
}

///////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
FullModelInfo::full_to_sub( utility::vector1< Size > const & res_in_full_model_numbering ) const {
	utility::vector1< Size > res;
	for ( Size n = 1; n <= res_in_full_model_numbering.size(); n++ ) {
		Size const i = res_list_.index( res_in_full_model_numbering[ n ] );
		if ( i > 0 ) res.push_back( i );
	}
	return res;
}

///////////////////////////////////////////////////////////////////////////////////////
Size
FullModelInfo::full_to_sub( Size const res_in_full_model_numbering ) const {
	runtime_assert( res_list_.index( res_in_full_model_numbering ) );
	return res_list_.index( res_in_full_model_numbering );
}

///////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
FullModelInfo::sub_to_full( utility::vector1< Size > const & res ) const {
	utility::vector1< Size > res_in_full_model_numbering;
	for ( Size n = 1; n <= res.size(); n++ ) {
		Size const & i = res[n];
		runtime_assert( i >= 1 && i <= res_list_.size() );
		res_in_full_model_numbering.push_back( res_list_[ i ] );
	}
	return res_in_full_model_numbering;
}

///////////////////////////////////////////////////////////////////////////////////////
Size
FullModelInfo::sub_to_full( Size const & res ) const {
	if ( res == 0 ) return 0;
	runtime_assert( res >= 1 && res <= res_list_.size() );
	return res_list_[ res ];
}

///////////////////////////////////////////////////////////////////////////////////////
Size
FullModelInfo::find_index_in_other_pose_list( pose::Pose const & pose) const {

	for ( Size n = 1; n <= other_pose_list_.size(); n++ ) {
		if ( other_pose_list_[ n ].get() == & pose ) return n;
	}
	return 0;
}

///////////////////////////////////////////////////////////////////////////////////////
void
FullModelInfo::clear_other_pose_list() {
	other_pose_list_.clear();
}

///////////////////////////////////////////////////////////////////////////////////////
void
FullModelInfo::clear_res_list() {
	res_list_.clear();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief return index of pose in other_pose_list
Size
FullModelInfo::get_idx_for_other_pose_with_residues( utility::vector1< Size > const & input_res ) const {
	Size const & idx = get_idx_for_other_pose_with_residue( input_res[ 1 ] );
	for ( Size n = 2; n <= input_res.size(); ++n ) {
		if ( idx != get_idx_for_other_pose_with_residue( input_res[ n ] ) ) {
			return 0;
		}
	}
	return idx;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Size
FullModelInfo::get_idx_for_other_pose_with_residue( Size const input_res ) const {
	Size idx( 0 );
	for ( Size i = 1; i <= other_pose_list_.size(); i++ ) {
		utility::vector1< Size > const & other_pose_res_list = const_full_model_info( *other_pose_list_[i] ).res_list();
		if ( other_pose_res_list.has_value( input_res ) ) {
			runtime_assert( idx == 0 ); // should be at most only one other pose with this residue number
			idx = i;
		}
	}
	return idx;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Size
FullModelInfo::get_idx_for_other_pose( pose::Pose const & pose ) const {
	Size idx( 0 );
	if ( pose.size() > 0 ) {
		Size const resnum = get_res_list_from_full_model_info_const( pose )[ 1 ];
		idx = get_idx_for_other_pose_with_residue( resnum );
	} else {
		for ( Size i = 1; i <= other_pose_list_.size(); i++ ) {
			if ( other_pose_list_[i]->size() == 0 ) {
				runtime_assert( idx == 0 ); // only one blank pose allowed in other_pose_list.
				idx = i;
			}
		}
	}
	runtime_assert( idx > 0 );
	return idx;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
FullModelInfo::add_other_pose( core::pose::PoseOP pose )
{
	other_pose_list_.push_back( pose );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
FullModelInfo::set_other_pose_list( utility::vector1< pose::PoseOP > const & setting ){
	other_pose_list_ = setting;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
FullModelInfo::remove_other_pose_at_idx( Size const idx ){

	runtime_assert( idx <= other_pose_list_.size() );
	utility::vector1< core::pose::PoseOP > other_pose_list_new;

	for ( Size i = 1; i <= other_pose_list_.size(); i++ ) {
		if ( i == idx ) continue;
		other_pose_list_new.push_back( other_pose_list_[ i ] );
	}

	other_pose_list_ = other_pose_list_new;
}


//////////////////////////////////////////////////////////////////////////////
/// @brief Return true if res_list_ matches a submotif exactly
bool
FullModelInfo::is_a_submotif() const
{
	return ( submotif_info_list_.size() > 0 && is_a_submotif( res_list_ ) );
}


//////////////////////////////////////////////////////////////////////////////
/// @brief Return true if res_list_ matches a submotif_seed exactly
bool
FullModelInfo::is_a_submotif_seed() const
{
	return ( submotif_info_list_.size() > 0 && is_a_submotif_seed( res_list_ ) );
}


//////////////////////////////////////////////////////////////////////////////
/// @brief Return true if residues match a submotif exactly
bool
FullModelInfo::is_a_submotif( utility::vector1< Size > const & res_list, bool const & check_other_poses /*= false*/ ) const
{
	// The !! converts pointer to boolean. For some reason, some compilers are not happy checking if the SubMotifInfoOP == 0.
	return ( !!submotif_info( res_list, check_other_poses ) );
}


//////////////////////////////////////////////////////////////////////////////
/// @brief Return true if residues match a submotif_seed exactly
bool
FullModelInfo::is_a_submotif_seed( utility::vector1< Size > const & res_list, bool const & check_other_poses /*= false*/ ) const
{
	return ( is_a_submotif( res_list, check_other_poses ) && submotif_info( res_list, check_other_poses )->seed() );
}


//////////////////////////////////////////////////////////////////////////////
/// @brief Return a submotif_info_op with matching residues (in full model
/// numbering).
SubMotifInfoOP
FullModelInfo::submotif_info( utility::vector1< Size > const & res_list, bool const & check_other_poses /*= false*/ ) const
{
	utility::vector1< Size > check_res_list = res_list;
	std::sort( check_res_list.begin(), check_res_list.end() );
	utility::vector1< SubMotifInfoOP >::iterator itr;
	for ( itr = submotif_info_list_.begin(); itr != submotif_info_list_.end(); ++itr ) {
		if ( (*itr)->sorted_res_list() == check_res_list ) return (*itr);
	}
	if ( check_other_poses ) {
		Size const & idx = get_idx_for_other_pose_with_residue( res_list[ 1 ] );
		if ( idx ) return const_full_model_info( *other_pose_list_[ idx ] ).submotif_info( res_list );
	}
	return NULL;
}


//////////////////////////////////////////////////////////////////////////////
/// @brief Return true if residues are found in a submotif
bool
FullModelInfo::in_a_submotif( utility::vector1< Size > const & res_list, bool const & check_other_poses /*= false*/ ) const
{
	// The !! converts pointer to boolean. For some reason, some compilers are not happy checking if the SubMotifInfoOP == 0.
	return ( !!submotif_info_containing_residues( res_list, check_other_poses ) );
}


//////////////////////////////////////////////////////////////////////////////
/// @brief Return true if residues are found in a submotif_seed
bool
FullModelInfo::in_a_submotif_seed( utility::vector1< Size > const & res_list, bool const & check_other_poses /*= false*/ ) const
{
	return ( in_a_submotif( res_list, check_other_poses ) && submotif_info_containing_residues( res_list, check_other_poses )->seed() );
}


//////////////////////////////////////////////////////////////////////////////
/// @brief Return the submotif_info_op that contains (but not limited to) input residues.
SubMotifInfoOP
FullModelInfo::submotif_info_containing_residues( utility::vector1< Size > const & res_list, bool const & check_other_poses /*= false*/ ) const
{
	utility::vector1< SubMotifInfoOP >::iterator itr;
	for ( itr = submotif_info_list_.begin(); itr != submotif_info_list_.end(); ++itr ) {
		bool found_all_residues( true );
		for ( Size n = 1; n <= res_list.size(); n++ ) {
			if ( (*itr)->res_list().has_value( res_list[ n ] ) ) continue;
			found_all_residues = false;
			break;
		}
		if ( !found_all_residues ) continue;
		return (*itr);
	}
	if ( check_other_poses ) {
		Size const & idx = get_idx_for_other_pose_with_residue( res_list[ 1 ] );
		if ( idx ) return const_full_model_info( *other_pose_list_[ idx ] ).submotif_info_containing_residues( res_list );
	}
	return NULL;
}


//////////////////////////////////////////////////////////////////////////////
/// @brief Delete submotif_infos with residues that are not found in res_list_.
void
FullModelInfo::update_submotif_info_list()
{
	utility::vector1< SubMotifInfoOP >::iterator itr;
	utility::vector1< SubMotifInfoOP > submotif_info_to_delete;
	for ( itr = submotif_info_list_.begin(); itr != submotif_info_list_.end(); ++itr ) {
		for ( Size n = 1; n <= (*itr)->res_list().size(); ++n ) {
			if ( !res_list_.has_value( (*itr)->res_list( n ) ) ) submotif_info_to_delete.push_back( *itr );
		}
	}

	// have to delete in a second stage; otherwise iteration above can get tripped up it
	// an iterator disapears in the middle of the loop.
	for ( Size n = 1; n <= submotif_info_to_delete.size(); n++ ) delete_submotif_info( submotif_info_to_delete[ n ] );

}


//////////////////////////////////////////////////////////////////////////////
/// @brief Create new SubMotifInfoOP and add it to the submotif_info_list_.
void
FullModelInfo::add_submotif_info(
	utility::vector1< Size > const & res_list,
	std::string const & tag,
	bool const & seed /*= false*/  )
{
	add_submotif_info( SubMotifInfoOP( new SubMotifInfo( res_list, tag, seed ) ) );
}


//////////////////////////////////////////////////////////////////////////////
/// @brief Add a list of SubMotifInfoOPs to the submotif_info_list_.
void
FullModelInfo::add_submotif_info( utility::vector1< SubMotifInfoOP > submotif_info_list ) {
	utility::vector1< SubMotifInfoOP >::iterator itr;
	for ( itr = submotif_info_list.begin(); itr != submotif_info_list.end(); ++itr ) {
		add_submotif_info( (*itr) );
	}
}


//////////////////////////////////////////////////////////////////////////////
/// @brief Add a SubMotifInfoOP to the submotif_info_list_.
void
FullModelInfo::add_submotif_info( SubMotifInfoOP submotif_info_op )
{
	// could try .index()/.has_value() -- but only if overloaded operator== is used
	utility::vector1< SubMotifInfoOP >::iterator itr;
	for ( itr = submotif_info_list_.begin(); itr != submotif_info_list_.end(); ++itr ) {
		if ( (*itr) == submotif_info_op ) return;
	}
	submotif_info_list_.push_back( submotif_info_op->clone() );
}


//////////////////////////////////////////////////////////////////////////////
/// @brief Delete SubMotifInfoOP (with equal values) from submotif_info_list_
void
FullModelInfo::delete_submotif_info(
	utility::vector1< Size > const & res_list,
	std::string const & tag,
	bool const & seed /*= false*/  )
{
	delete_submotif_info( SubMotifInfoOP( new SubMotifInfo( res_list, tag, seed ) ) );
}


//////////////////////////////////////////////////////////////////////////////
/// @brief Delete a list of SubMotifInfoOPs from submotif_info_list_
void
FullModelInfo::delete_submotif_info( utility::vector1< SubMotifInfoOP > submotif_info_list ) {
	utility::vector1< SubMotifInfoOP >::iterator itr;
	for ( itr = submotif_info_list.begin(); itr != submotif_info_list.end(); ++itr ) {
		delete_submotif_info( (*itr) );
	}
}


//////////////////////////////////////////////////////////////////////////////
/// @brief Delete SubMotifInfoOP in submotif_info_list_ if values are equal
/// to submotif_info. Should submotif_info deleted if seed() == true?
void
FullModelInfo::delete_submotif_info( SubMotifInfoOP submotif_info_op )
{
	utility::vector1< SubMotifInfoOP >::iterator itr;
	for ( itr = submotif_info_list_.begin(); itr != submotif_info_list_.end(); ++itr ) {
		if ( (*itr) == submotif_info_op ) {
			submotif_info_list_.erase( itr );
			break;
		}
	}
}


//////////////////////////////////////////////////////////////////////////////
/// @brief Delete SubMotifInfoOP in submotif_info_list_ if values are equal
/// to submotif_info. Should submotif_info deleted if seed() == true?
void
FullModelInfo::show_submotif_info_list( bool const & for_all_poses /*= false*/ ) const
{
	if ( !submotif_info_list_.size() ) return;
	utility::vector1< SubMotifInfoOP >::iterator itr;
	std::cout << "SUBMOTIF_INFO FOR POSE" << std::endl;
	for ( itr = submotif_info_list_.begin(); itr != submotif_info_list_.end(); ++itr ) {
		std::cout << *itr << std::endl;
	}
	if ( for_all_poses ) {
		for ( Size pose_idx = 1; pose_idx <= other_pose_list_.size(); ++pose_idx ) {
			std::cout << "SUBMOTIF_INFO FOR OTHER POSE [IDX: " << pose_idx << "]" << std::endl;
			const_full_model_info( *other_pose_list_[ pose_idx ] ).show_submotif_info_list();
		}
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// @details Pose must already contain a full_model_info object or this method will fail.
FullModelInfo const &
const_full_model_info( pose::Pose const & pose )
{
	debug_assert( pose.data().has( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO ) );
	return *( utility::pointer::static_pointer_cast< core::pose::full_model_info::FullModelInfo const > ( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO) ) );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details Either returns a non-const reference to the FullModelInfo object already stored
/// in the pose, or creates a new FullModelInfo object, places it in the pose, and returns
/// a non-const reference to it.
FullModelInfo &
nonconst_full_model_info( pose::Pose & pose )
{

	if ( pose.data().has( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO ) ) {
		// following takes time -- so we shouldn't call this nonconst function a lot.
		runtime_assert ( check_full_model_info_OK( pose ) ); // later remove this for speed.
		return *( utility::pointer::static_pointer_cast< core::pose::full_model_info::FullModelInfo > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO ) ));
	}

	FullModelInfoOP full_model_info( new FullModelInfo( pose ) );

	pose.data().set( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO, full_model_info );
	return *full_model_info;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
full_model_info_defined( pose::Pose const & pose ){
	return pose.data().has( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO );
}

//////////////////////////////////////////////////////////////////////////////
// basically an alias
FullModelInfo const &
make_sure_full_model_info_is_setup( pose::Pose & pose )
{
	return nonconst_full_model_info( pose );
}

//////////////////////////////////////////////////////////////////////////////
void
set_full_model_info( pose::Pose & pose, FullModelInfoOP & full_model_info ){
	pose.data().set( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO, full_model_info );
	update_pdb_info_from_full_model_info( pose );
}

//////////////////////////////////////////////////////////////////////////////
void
update_full_model_info_from_pose( pose::Pose & pose ){
	FullModelInfoOP full_model_info( new FullModelInfo( pose ) );
	set_full_model_info( pose, full_model_info );
}


} //full_model_info
} //pose
} //core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::pose::full_model_info::FullModelInfo::FullModelInfo() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pose::full_model_info::FullModelInfo::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( res_list_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( other_pose_list_ ) ); // utility::vector1<core::pose::PoseOP>
	arc( CEREAL_NVP( full_model_parameters_ ) ); // FullModelParametersCOP
	arc( CEREAL_NVP( submotif_info_list_ ) ); // utility::vector1<SubMotifInfoOP>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pose::full_model_info::FullModelInfo::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( res_list_ ); // utility::vector1<Size>
	arc( other_pose_list_ ); // utility::vector1<core::pose::PoseOP>

	std::shared_ptr< core::pose::full_model_info::FullModelParameters > local_full_model_parameters;
	arc( local_full_model_parameters ); // FullModelParametersCOP
	full_model_parameters_ = local_full_model_parameters; // copy the non-const pointer(s) into the const pointer(s)

	arc( CEREAL_NVP( submotif_info_list_ ) ); // utility::vector1<SubMotifInfoOP>
}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::full_model_info::FullModelInfo );
CEREAL_REGISTER_TYPE( core::pose::full_model_info::FullModelInfo )

CEREAL_REGISTER_DYNAMIC_INIT( core_pose_full_model_info_FullModelInfo )
#endif // SERIALIZATION
