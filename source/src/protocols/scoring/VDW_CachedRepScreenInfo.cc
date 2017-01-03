// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/scoring/VDW_CachedRepScreenInfo.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


///#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/scoring/VDW_CachedRepScreenInfo.hh>
#include <core/pose/rna/VDW_RepScreenInfo.hh>
#include <core/pose/rna/VDW_Grid.hh>
//#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.hh> ???
///#include <protocols/farna/util.hh>
#include <core/pose/rna/util.hh>

#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/rna/util.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <utility/vector1.hh>

#include <basic/Tracer.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>

static basic::Tracer TR( "protocols.stepwise.modeler.rna.checker.VDW_CachedRepScreenInfo" );

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace scoring {

// @brief Constructor
VDW_CachedRepScreenInfo::VDW_CachedRepScreenInfo() :
	CacheableData(),
	VDW_screen_bin_( core::pose::rna::VDW_GridCOP( new core::pose::rna::VDW_Grid() ) )
{
	read_in_VDW_rep_screen_pose_from_command_line();
}

// @brief Constructor
VDW_CachedRepScreenInfo::VDW_CachedRepScreenInfo( core::pose::Pose const & pose ) :
	CacheableData()
{
	VDW_rep_screen_info_list_ = const_vdw_cached_rep_screen_info_from_pose( pose ).VDW_rep_screen_info_list();
	VDW_screen_bin_ = const_vdw_cached_rep_screen_info_from_pose( pose ).VDW_screen_bin();
}

/// @details Copy constructors must copy all data, not just some...
VDW_CachedRepScreenInfo::VDW_CachedRepScreenInfo( VDW_CachedRepScreenInfo const & src ) :
	CacheableData(),
	VDW_rep_screen_info_list_( src.VDW_rep_screen_info_list_ ),
	VDW_screen_bin_( src.VDW_screen_bin_ ) // note, we are not cloning -- this is 'scratch' space that takes a while to initialize.
{
	// we will not copy over VDW_screen_bin, which is a huge 3D grid, and
	// most of the time is zero.
	// only instantiate that when we need it.
}


// @brief default destructor
VDW_CachedRepScreenInfo::~VDW_CachedRepScreenInfo()
{}


basic::datacache::CacheableDataOP
VDW_CachedRepScreenInfo::clone() const
{
	return basic::datacache::CacheableDataOP( new VDW_CachedRepScreenInfo( *this ) );
}


utility::vector1< core::pose::rna::VDW_RepScreenInfo > &
VDW_CachedRepScreenInfo::VDW_rep_screen_info_list() const
{
	return VDW_rep_screen_info_list_;
}


core::pose::rna::VDW_GridCOP
VDW_CachedRepScreenInfo::VDW_screen_bin() const
{
	return VDW_screen_bin_;
}


void
VDW_CachedRepScreenInfo::read_in_VDW_rep_screen_pose( core::pose::rna::VDW_RepScreenInfo & VDW_rep_screen_info ) const
{
	using namespace core::chemical;

	TR.Debug << "importing VDW_rep_screen_pose: " << VDW_rep_screen_info.pose_name << std::endl;
	ResidueTypeSetCOP rsd_set( /*core::chemical::*/ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );

	if ( VDW_rep_screen_info.pose_name == "" ) utility_exit_with_message( VDW_rep_screen_info.pose_name == "" );
	VDW_rep_screen_info.VDW_pose = core::pose::PoseOP( new core::pose::Pose );

	core::import_pose::pose_from_file( *VDW_rep_screen_info.VDW_pose, *rsd_set, VDW_rep_screen_info.pose_name, core::import_pose::PDB_file );
	core::pose::rna::make_phosphate_nomenclature_matches_mini( *VDW_rep_screen_info.VDW_pose );
	core::pose::rna::add_virtual_O2Prime_hydrogen( *VDW_rep_screen_info.VDW_pose );
}


void
VDW_CachedRepScreenInfo::read_in_VDW_rep_screen_pose_from_command_line() const
{
	using namespace basic::options;

	if ( !option[ OptionKeys::stepwise::rna::VDW_rep_screen_info ].user() ) return;
	utility::vector1< std::string > const & All_VDW_rep_screen_pose_info = option[ OptionKeys::stepwise::rna::VDW_rep_screen_info ]();

	bool align_res_specified = false;
	if ( All_VDW_rep_screen_pose_info.size() > 1 ) {
		std::string str = All_VDW_rep_screen_pose_info[2];
		align_res_specified = ( str.find(".pdb", str.size()-4) == std::string::npos );
	}

	for ( Size n = 1; n <= All_VDW_rep_screen_pose_info.size(); n++ ) {

		if ( All_VDW_rep_screen_pose_info[n].find(".pdb", All_VDW_rep_screen_pose_info[n].size()-4) == std::string::npos ) continue; // must have '.pdb'

		core::pose::rna::VDW_RepScreenInfo VDW_rep_screen_info = core::pose::rna::VDW_RepScreenInfo();
		VDW_rep_screen_info.pose_name = All_VDW_rep_screen_pose_info[n];
		read_in_VDW_rep_screen_pose( VDW_rep_screen_info );

		// Need to set the VDW_align_res and the working_align_res

		utility::vector1< std::string > VDW_align_res_string, working_align_res_string;
		if ( align_res_specified ) {
			VDW_align_res_string = core::pose::rna::tokenize( All_VDW_rep_screen_pose_info[n + 1], "-");
			working_align_res_string = core::pose::rna::tokenize( All_VDW_rep_screen_pose_info[n + 2], "-");
			for ( Size i = 1; i <= VDW_align_res_string.size(); ++i ) {
				//VDW_rep_screen_info.VDW_align_res.push_back( protocols::stepwise::modeler::rna::string_to_int( VDW_align_res_string[i] ) );
				VDW_rep_screen_info.VDW_align_res.push_back( core::pose::rna::string_to_int( VDW_align_res_string[i] ) );
			}
			for ( Size i = 1; i <= working_align_res_string.size(); ++i ) {
				VDW_rep_screen_info.working_align_res.push_back( core::pose::rna::string_to_int( working_align_res_string[i] ) );
			}
		}

		VDW_rep_screen_info_list_.push_back( VDW_rep_screen_info );
	}
}


/// @details Pose must already contain a vdw_cached_rep_screen_info object or this method will fail.
VDW_CachedRepScreenInfo const &
const_vdw_cached_rep_screen_info_from_pose( core::pose::Pose const & pose )
{
	assert( pose.data().has( core::pose::datacache::CacheableDataType::VDW_REP_SCREEN_INFO ) );
	return *( utility::pointer::static_pointer_cast< VDW_CachedRepScreenInfo const > ( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::VDW_REP_SCREEN_INFO ) ) );
}


/// @details Either returns a non-const reference to the vdw_cached_rep_screen_info object already stored
/// in the pose, or creates a new vdw_cached_rep_screen_info object, places it in the pose, and returns
/// a non-const reference to it.
VDW_CachedRepScreenInfo &
nonconst_vdw_cached_rep_screen_info_from_pose( core::pose::Pose & pose )
{
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::VDW_REP_SCREEN_INFO ) ) {
		pose.data().set( core::pose::datacache::CacheableDataType::VDW_REP_SCREEN_INFO, VDW_CachedRepScreenInfoOP( new VDW_CachedRepScreenInfo() ) );
	}
	assert( pose.data().has( core::pose::datacache::CacheableDataType::VDW_REP_SCREEN_INFO ) );
	return *( utility::pointer::static_pointer_cast< VDW_CachedRepScreenInfo > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::VDW_REP_SCREEN_INFO ) ) );
}


//////////////////////////////////////////////////////////////////////////////
// basically an alias
VDW_CachedRepScreenInfo const &
make_sure_vdw_cached_rep_screen_info_is_setup( core::pose::Pose & pose )
{
	return nonconst_vdw_cached_rep_screen_info_from_pose( pose );
}


//////////////////////////////////////////////////////////////////////////////
bool
vdw_cached_rep_screen_info_is_setup( core::pose::Pose const & pose )
{
	return pose.data().has( core::pose::datacache::CacheableDataType::VDW_REP_SCREEN_INFO );
}


//////////////////////////////////////////////////////////////////////////////
bool
option_vdw_rep_screen_info_user()
{
	return basic::options::option[ basic::options::OptionKeys::stepwise::rna::VDW_rep_screen_info ].user();
}


//////////////////////////////////////////////////////////////////////////////
void
set_vdw_cached_rep_screen_info( core::pose::Pose & pose, VDW_CachedRepScreenInfoOP & vdw_cached_rep_screen_info ){
	pose.data().set( core::pose::datacache::CacheableDataType::VDW_REP_SCREEN_INFO, vdw_cached_rep_screen_info );
}


//////////////////////////////////////////////////////////////////////////////
void
set_vdw_cached_rep_screen_info_from_pose( core::pose::Pose & new_pose, core::pose::Pose const & pose ){
	if ( !option_vdw_rep_screen_info_user() ) return;
	VDW_CachedRepScreenInfoOP vdw_cached_rep_screen_info( new VDW_CachedRepScreenInfo( pose ) );
	set_vdw_cached_rep_screen_info( new_pose, vdw_cached_rep_screen_info );
}


///////////////////////////////////////////////////////////////////////////////////////
void
fill_vdw_cached_rep_screen_info_from_command_line( core::pose::Pose & pose ) {
	if ( !option_vdw_rep_screen_info_user() ) return;
	make_sure_vdw_cached_rep_screen_info_is_setup( pose );
}


///////////////////////////////////////////////////////////////////////////////////////
void
fill_vdw_cached_rep_screen_info_from_command_line( utility::vector1< core::pose::Pose * > & pose_pointers ) {
	if ( !option_vdw_rep_screen_info_user() ) return;
	for ( auto const & poseop : pose_pointers ) {
		fill_vdw_cached_rep_screen_info_from_command_line( *poseop );
	}
}

} //scoring
} //protocols

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::scoring::VDW_CachedRepScreenInfo::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( VDW_rep_screen_info_list_ ) ); // utility::vector1<core::pose::rna::VDW_RepScreenInfo>
	arc( CEREAL_NVP( VDW_screen_bin_ ) ); // core::pose::rna::VDW_GridCOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
 protocols::scoring::VDW_CachedRepScreenInfo::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( VDW_rep_screen_info_list_ ); // utility::vector1<core::pose::rna::VDW_RepScreenInfo>
	std::shared_ptr< core::pose::rna::VDW_Grid > local_VDW_screen_bin;
	arc( local_VDW_screen_bin ); // core::pose::rna::VDW_GridCOP
	VDW_screen_bin_ = local_VDW_screen_bin; // copy the non-const pointer(s) into the const pointer(s)
}

SAVE_AND_LOAD_SERIALIZABLE(  protocols::scoring::VDW_CachedRepScreenInfo );
CEREAL_REGISTER_TYPE(  protocols::scoring::VDW_CachedRepScreenInfo )

CEREAL_REGISTER_DYNAMIC_INIT(  protocols_scoring_VDW_CachedRepScreenInfo )
#endif // SERIALIZATION
