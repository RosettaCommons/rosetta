// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/operation/TaskOperations.cc
/// @brief
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Unit Headers
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/TaskOperationCreators.hh>

// Project Headers
#include <core/id/SequenceMapping.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/rotamer_set/RotamerCouplings.hh>
#include <core/pack/rotamer_set/RotamerLinks.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/selection.hh>
#include <core/pose/datacache/CacheableObserverType.hh>
#include <core/pose/datacache/cacheable_observers.hh>
#include <core/pose/datacache/ObserverCache.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <boost/foreach.hpp>

// Utility Headers
#include <utility/exit.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>

// basic headers
#include <basic/resource_manager/ResourceManager.hh>

// option key includes
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace task {
namespace operation {

using basic::t_warning;
using basic::t_info;
using basic::t_debug;
static thread_local basic::Tracer TR( "core.pack.task.operation.TaskOperations", t_info );
using namespace utility::tag;

/// BEGIN RestrictToRepacking

RestrictToRepacking::~RestrictToRepacking() {}

TaskOperationOP RestrictToRepackingCreator::create_task_operation() const
{
	return TaskOperationOP( new RestrictToRepacking );
}

TaskOperationOP RestrictToRepacking::clone() const
{
	return TaskOperationOP( new RestrictToRepacking( *this ) );
}

void
RestrictToRepacking::apply( pose::Pose const &, PackerTask & task ) const
{
	task.restrict_to_repacking();
}

void
RestrictToRepacking::parse_tag( TagCOP, DataMap & )
{}

void
RestrictToRepacking::parse_def( utility::lua::LuaObject const & ) {}

/// BEGIN RestrictResidueToRepacking
RestrictResidueToRepacking::~RestrictResidueToRepacking() {}

TaskOperationOP RestrictResidueToRepackingCreator::create_task_operation() const
{
	return TaskOperationOP( new RestrictResidueToRepacking );
}

TaskOperationOP RestrictResidueToRepacking::clone() const
{
	return TaskOperationOP( new RestrictResidueToRepacking( *this ) );
}

void
RestrictResidueToRepacking::apply( pose::Pose const &, PackerTask & task ) const
{
	for(utility::vector1< core::Size >::const_iterator it(residues_to_restrict_to_repacking_.begin()), end(residues_to_restrict_to_repacking_.end());
			it != end; ++it)
		{
		//debug_assert( *it ); //begin/end can't return NULL anyway?
			task.nonconst_residue_task(*it).restrict_to_repacking();
		}
	return;
}

void
RestrictResidueToRepacking::include_residue( core::Size resid ) { residues_to_restrict_to_repacking_.push_back(resid); }

void
RestrictResidueToRepacking::clear() { residues_to_restrict_to_repacking_.clear(); }

void
RestrictResidueToRepacking::parse_tag( TagCOP tag , DataMap & )
{
	include_residue( tag->getOption< core::Size >( "resnum", 0 ) );
}

/// BEGIN RestrictAbsentCanonicalAAS
RestrictAbsentCanonicalAAS::RestrictAbsentCanonicalAAS()
	:	parent(),
		resid_(1)
{
	keep_aas_.assign( chemical::num_canonical_aas, false );
}

RestrictAbsentCanonicalAAS::RestrictAbsentCanonicalAAS( core::Size resid, utility::vector1< bool > keep  )
	:
	parent()
{
	keep_aas( keep );
	include_residue( resid );
}

RestrictAbsentCanonicalAAS::~RestrictAbsentCanonicalAAS(){}

TaskOperationOP RestrictAbsentCanonicalAASCreator::create_task_operation() const
{
	return TaskOperationOP( new RestrictAbsentCanonicalAAS );
}

TaskOperationOP RestrictAbsentCanonicalAAS::clone() const
{
	return TaskOperationOP( new RestrictAbsentCanonicalAAS( *this ) );
}

void
RestrictAbsentCanonicalAAS::apply( pose::Pose const &, PackerTask & task ) const
{
	if( resid_ == 0 ){// restrict all residues
		for( core::Size i( 1 ); i <= task.total_residue(); ++i ) {
			task.nonconst_residue_task( i ).restrict_absent_canonical_aas( keep_aas_ );
		}
	}
	else
		task.nonconst_residue_task( resid_ ).restrict_absent_canonical_aas( keep_aas_ );
}

void
RestrictAbsentCanonicalAAS::keep_aas( std::string const keep )
{
	using namespace chemical;
	runtime_assert( keep_aas_.size() == num_canonical_aas );
	utility::vector1< bool > canonical_aas_to_keep( num_canonical_aas, false );
	BOOST_FOREACH( char const c, keep ){
		if ( oneletter_code_specifies_aa( c ) ) {
			//std::cout << "Keeping amino acid " << c << std::endl;
			canonical_aas_to_keep[ aa_from_oneletter_code( c ) ] = true;
		} else {
			TR << "aa letter " << c << " does not not correspond to a canonical AA"<<std::endl;
			utility_exit();
		}
	}
	keep_aas( canonical_aas_to_keep );
}

// if an amino acid is not present (false) in the boolean vector, then do not allow it at this position.  The boolean vector is a 20-length vector in alphabetical order by one-letter code. 
void
RestrictAbsentCanonicalAAS::keep_aas( utility::vector1< bool > const keep )
{
	runtime_assert( keep.size() == chemical::num_canonical_aas );
	keep_aas_ = keep;
	//for ( Size ii = 1; ii <= keep_aas_.size(); ++ii ) {
	//	std::cout << " keeping " << ii << " " << " ? " << keep_aas_[ ii ]  << std::endl;
	//}
}

void
RestrictAbsentCanonicalAAS::include_residue( core::Size const resid )
{
	resid_ = resid;
}

void
RestrictAbsentCanonicalAAS::parse_tag( TagCOP tag , DataMap & )
{
	include_residue( tag->getOption< core::Size >( "resnum", 0 ) );
	keep_aas( tag->getOption< std::string >( "keep_aas" ) );
}


//BEGIN DisallowIfNonnative
DisallowIfNonnative::DisallowIfNonnative():
	disallowed_aas_ ( chemical::num_canonical_aas, false ),
	allowed_aas_( invert_vector( disallowed_aas_ ) )
{}

DisallowIfNonnative::DisallowIfNonnative( utility::vector1< bool > disallowed_aas):
	disallowed_aas_( disallowed_aas ),
	allowed_aas_( invert_vector(disallowed_aas) )
{}

DisallowIfNonnative::DisallowIfNonnative( utility::vector1< bool > disallowed_aas, utility::vector1<core::Size> res_selection):
	residue_selection_( res_selection ),
	disallowed_aas_( disallowed_aas ),
	allowed_aas_( invert_vector(disallowed_aas) )
{}

DisallowIfNonnative::~DisallowIfNonnative(){}

TaskOperationOP DisallowIfNonnativeCreator::create_task_operation() const
{
	return TaskOperationOP( new DisallowIfNonnative );
}

TaskOperationOP DisallowIfNonnative::clone() const
{
	return TaskOperationOP( new DisallowIfNonnative( *this ) );
}

void DisallowIfNonnative::clear(){
	allowed_aas_.clear();
	disallowed_aas_.clear();
	residue_selection_.clear();
}

//private function to invert disallowed aas into allowed aas
utility::vector1< bool >
DisallowIfNonnative::invert_vector( utility::vector1< bool > disallowed_aas){
	utility::vector1< bool > inverted_vec;
	for(core::Size ii=1; ii<=disallowed_aas_.size(); ii++ ){
		inverted_vec.push_back( ! disallowed_aas[ii] );
	}
	return inverted_vec;
}

void DisallowIfNonnative::apply( pose::Pose const &, PackerTask & task ) const
{
 	runtime_assert( allowed_aas_.size() == chemical::num_canonical_aas );
	if ( residue_selection_.empty() ){ //if no residue defined then do all residues
		for (core::Size ii = 1; ii<= task.total_residue(); ii++)
			task.nonconst_residue_task( ii ).restrict_nonnative_canonical_aas( allowed_aas_ );
	}
	else{  //if a residue is defined then do only on that residue
		for(core::Size jj = 1; jj <= residue_selection_.size(); jj++ )
			task.nonconst_residue_task( residue_selection_[jj] ).restrict_nonnative_canonical_aas( allowed_aas_ );
	}
}

//helper functions for DisallowIfNonnative
void
DisallowIfNonnative::disallow_aas( utility::vector1< bool > const & cannonical_disallowed ){
	runtime_assert( cannonical_disallowed.size() == chemical::num_canonical_aas );
	disallowed_aas_ = cannonical_disallowed;
	allowed_aas_ = invert_vector( disallowed_aas_ );
}
void DisallowIfNonnative::disallow_aas( std::string const & aa_string ){
	using namespace chemical;
	utility::vector1< bool > aa_vector ( chemical::num_canonical_aas, false );
	for ( std::string::const_iterator it( aa_string.begin() ), end( aa_string.end() );
				it != end; ++it ) {
		if ( oneletter_code_specifies_aa( *it ) ) {
			aa_vector[ aa_from_oneletter_code( *it ) ] = true;
		} else {
			std::ostringstream os;
			os << "aa letter " << *it << " does not not correspond to a canonical AA";
			utility_exit_with_message( os.str() );
		}
	}
	disallowed_aas_ = aa_vector;
	allowed_aas_ = invert_vector( disallowed_aas_ );
}
	//functions to restrict what residues are looked at by operation
	//selections are additive
void DisallowIfNonnative::restrict_to_residue( core::Size const & resid){
	//chrisk: default to all residues if bogus seqpos 0 passed (e.g. through rosetta_scripts)
	if( resid != 0 ) residue_selection_.push_back( resid );
}
void DisallowIfNonnative::restrict_to_residue( utility::vector1< core::Size > const & residues){
	for(core::Size ii=1; ii<=residues.size(); ii++)
		residue_selection_.push_back( residues[ii] );
}

void DisallowIfNonnative::parse_tag( TagCOP tag , DataMap & )
{
	restrict_to_residue( tag->getOption< core::Size >( "resnum", 0 ) );
	disallow_aas( tag->getOption< std::string >( "disallow_aas" ) );
}
void DisallowIfNonnative::parse_def( utility::lua::LuaObject const & def) {
	restrict_to_residue( def["resnum"] ? def["resnum"].to<core::Size>() : 0 );
	disallow_aas( def["disallow_aas"].to<std::string>());
}

//BEGIN RotamerExplosion
RotamerExplosion::RotamerExplosion(){}

RotamerExplosion::RotamerExplosion( core::Size const resid, ExtraRotSample const sample_level, core::Size const chi ) :
	resid_( resid ),
	chi_( chi ),
	sample_level_( sample_level )
{}

RotamerExplosion::~RotamerExplosion() {}

TaskOperationOP RotamerExplosionCreator::create_task_operation() const
{
	return TaskOperationOP( new RotamerExplosion );
}

TaskOperationOP RotamerExplosion::clone() const
{
	return TaskOperationOP( new RotamerExplosion( *this ) );
}

void
RotamerExplosion::apply( core::pose::Pose const &, PackerTask & task ) const
{
	ResidueLevelTask & restask( task.nonconst_residue_task( resid_ ) );
  if( chi_ > 0 ) restask.or_ex1_sample_level( sample_level_ );
	if( chi_ > 1 ) restask.or_ex2_sample_level( sample_level_ );
	if( chi_ > 2 ) restask.or_ex3_sample_level( sample_level_ );
	if( chi_ > 3 ) restask.or_ex4_sample_level( sample_level_ );
	restask.or_include_current( false ); //not that this call does anything...
}

void
RotamerExplosion::parse_tag( TagCOP tag , DataMap & )
{
	resid( tag->getOption< core::Size >( "resnum" ) );
	chi( tag->getOption< core::Size >( "chi" ) );
	sample_level( EX_THREE_THIRD_STEP_STDDEVS );// hardcoded and ugly, but this is probably too much to expect the user to set
}

void
RotamerExplosion::resid( core::Size const r )
{
	resid_ = r;
}

void
RotamerExplosion::chi( core::Size const c )
{
	chi_ = c;
}

void
RotamerExplosion::sample_level( ExtraRotSample const s )
{
	sample_level_ = s;
}


/// BEGIN InitializeFromCommandline

InitializeFromCommandline::~InitializeFromCommandline() {}

TaskOperationOP InitializeFromCommandlineCreator::create_task_operation() const
{
	return TaskOperationOP( new InitializeFromCommandline );
}

TaskOperationOP InitializeFromCommandline::clone() const
{
	return TaskOperationOP( new InitializeFromCommandline( *this ) );
}

void
InitializeFromCommandline::apply( pose::Pose const &, PackerTask & task ) const
{
	task.initialize_from_command_line();
}

void
InitializeFromCommandline::parse_tag( TagCOP, DataMap & )
{
}

void
InitializeFromCommandline::parse_def( utility::lua::LuaObject const &) {}

/// BEGIN InitializeFromCommandline

InitializeExtraRotsFromCommandline::~InitializeExtraRotsFromCommandline() {}

TaskOperationOP InitializeExtraRotsFromCommandlineCreator::create_task_operation() const
{
	return TaskOperationOP( new InitializeExtraRotsFromCommandline );
}

TaskOperationOP InitializeExtraRotsFromCommandline::clone() const
{
	return TaskOperationOP( new InitializeExtraRotsFromCommandline( *this ) );
}

void
InitializeExtraRotsFromCommandline::apply( pose::Pose const &, PackerTask & task ) const
{
	task.initialize_extra_rotamer_flags_from_command_line();
}


/// BEGIN IncludeCurrent

IncludeCurrent::~IncludeCurrent() {}

TaskOperationOP IncludeCurrentCreator::create_task_operation() const
{
	return TaskOperationOP( new IncludeCurrent );
}

TaskOperationOP IncludeCurrent::clone() const
{
	return TaskOperationOP( new IncludeCurrent( *this ) );
}

void
IncludeCurrent::apply( pose::Pose const &, PackerTask & task ) const
{
	task.or_include_current(true);
}

void
IncludeCurrent::parse_def( utility::lua::LuaObject const & ) {}

/// BEGIN ExtraRotamersGeneric

ExtraRotamersGeneric::ExtraRotamersGeneric() :
	ex1_(false),
	ex2_(false),
	ex3_(false),
	ex4_(false),
	ex1aro_(false),
	ex2aro_(false),
	ex1aro_exposed_(false),
	ex2aro_exposed_(false),
	ex1_sample_level_( NO_EXTRA_CHI_SAMPLES ),
	ex2_sample_level_( NO_EXTRA_CHI_SAMPLES ),
	ex3_sample_level_( NO_EXTRA_CHI_SAMPLES ),
	ex4_sample_level_( NO_EXTRA_CHI_SAMPLES ),
	ex1aro_sample_level_( NO_EXTRA_CHI_SAMPLES ),
	ex2aro_sample_level_( NO_EXTRA_CHI_SAMPLES ),
	ex1aro_exposed_sample_level_( NO_EXTRA_CHI_SAMPLES ),
	ex2aro_exposed_sample_level_( NO_EXTRA_CHI_SAMPLES ),
	exdna_sample_level_( NO_EXTRA_CHI_SAMPLES ),
	extrachi_cutoff_( EXTRACHI_CUTOFF_LIMIT )
{}

ExtraRotamersGeneric::~ExtraRotamersGeneric() {}

TaskOperationOP ExtraRotamersGenericCreator::create_task_operation() const
{
	return TaskOperationOP( new ExtraRotamersGeneric );
}

TaskOperationOP ExtraRotamersGeneric::clone() const
{
	return TaskOperationOP( new ExtraRotamersGeneric( *this ) );
}

void
ExtraRotamersGeneric::parse_tag( TagCOP tag , DataMap & )
{
	ex1_ = tag->getOption<bool>("ex1", false);
	ex2_ = tag->getOption<bool>("ex2", false);
	ex3_ = tag->getOption<bool>("ex3", false);
	ex4_ = tag->getOption<bool>("ex4", false);
	ex1aro_ = tag->getOption<bool>("ex1aro", false);
	ex2aro_ = tag->getOption<bool>("ex2aro", false);
	ex1aro_exposed_ = tag->getOption<bool>("ex1aro_exposed", false);
	ex2aro_exposed_ = tag->getOption<bool>("ex2aro_exposed", false);

	ex1_sample_level_ = static_cast<ExtraRotSample>(tag->getOption<Size>("ex1_sample_level", NO_EXTRA_CHI_SAMPLES));
	ex2_sample_level_ = static_cast<ExtraRotSample>(tag->getOption<Size>("ex2_sample_level", NO_EXTRA_CHI_SAMPLES));
	ex3_sample_level_ = static_cast<ExtraRotSample>(tag->getOption<Size>("ex3_sample_level", NO_EXTRA_CHI_SAMPLES));
	ex4_sample_level_ = static_cast<ExtraRotSample>(tag->getOption<Size>("ex4_sample_level", NO_EXTRA_CHI_SAMPLES));
	ex1aro_sample_level_ = static_cast<ExtraRotSample>(tag->getOption<Size>("ex1aro_sample_level", NO_EXTRA_CHI_SAMPLES));
	ex2aro_sample_level_ = static_cast<ExtraRotSample>(tag->getOption<Size>("ex2aro_sample_level", NO_EXTRA_CHI_SAMPLES));
	ex1aro_exposed_sample_level_ = static_cast<ExtraRotSample>(tag->getOption<Size>("ex1aro_exposed_sample_level", NO_EXTRA_CHI_SAMPLES));
	ex2aro_exposed_sample_level_ = static_cast<ExtraRotSample>(tag->getOption<Size>("ex2aro_exposed_sample_level", NO_EXTRA_CHI_SAMPLES));
	exdna_sample_level_ = static_cast<ExtraRotSample>(tag->getOption<Size>("exdna_sample_level", NO_EXTRA_CHI_SAMPLES));

	extrachi_cutoff_ = tag->getOption<Size>("extrachi_cutoff", EXTRACHI_CUTOFF_LIMIT);
}


void
ExtraRotamersGeneric::apply( pose::Pose const &, PackerTask & task ) const
{
	for(Size i=1; i <= task.total_residue(); ++i){
		ResidueLevelTask & res_task(task.nonconst_residue_task(i));
		res_task.or_ex1(ex1_);
		res_task.or_ex2(ex2_);
		res_task.or_ex3(ex3_);
		res_task.or_ex4(ex4_);
		res_task.or_ex1aro(ex1aro_);
		res_task.or_ex2aro(ex2aro_);
		res_task.or_ex1aro_exposed(ex1aro_exposed_);
		res_task.or_ex2aro_exposed(ex2aro_exposed_);
		res_task.or_ex1_sample_level(ex1_sample_level_);
		res_task.or_ex2_sample_level(ex2_sample_level_);
		res_task.or_ex3_sample_level(ex3_sample_level_);
		res_task.or_ex4_sample_level(ex4_sample_level_);
		res_task.or_ex1aro_sample_level(ex1aro_sample_level_);
		res_task.or_ex2aro_sample_level(ex2aro_sample_level_);
		res_task.or_ex1aro_exposed_sample_level(ex1aro_exposed_sample_level_);
		res_task.or_ex2aro_exposed_sample_level(ex2aro_exposed_sample_level_);
		res_task.or_exdna_sample_level(exdna_sample_level_);
		res_task.and_extrachi_cutoff(extrachi_cutoff_);
	}
}

void ExtraRotamersGeneric::ex1( bool value ) {
	ex1_ = value;
}
void ExtraRotamersGeneric::ex2( bool value ) {
	ex2_ = value;
}
void ExtraRotamersGeneric::ex3( bool value ) {
	ex3_ = value;
}
void ExtraRotamersGeneric::ex4( bool value ) {
	ex4_ = value;
}
void ExtraRotamersGeneric::ex1aro( bool value ) {
	ex1aro_ = value;
}
void ExtraRotamersGeneric::ex2aro( bool value ) {
	ex2aro_ = value;
}
void ExtraRotamersGeneric::ex1aro_exposed( bool value ) {
	ex1aro_exposed_ = value;
}
void ExtraRotamersGeneric::ex2aro_exposed( bool value ) {
	ex2aro_exposed_ = value;
}
void ExtraRotamersGeneric::ex1_sample_level( ExtraRotSample value ) {
	ex1_sample_level_ = value;
}
void ExtraRotamersGeneric::ex2_sample_level( ExtraRotSample value ) {
	ex2_sample_level_ = value;
}
void ExtraRotamersGeneric::ex3_sample_level( ExtraRotSample value ) {
	ex3_sample_level_ = value;
}
void ExtraRotamersGeneric::ex4_sample_level( ExtraRotSample value ) {
	ex4_sample_level_ = value;
}
void ExtraRotamersGeneric::ex1aro_sample_level( ExtraRotSample value ) {
	ex1aro_sample_level_ = value;
}
void ExtraRotamersGeneric::ex2aro_sample_level( ExtraRotSample value ) {
	ex2aro_sample_level_ = value;
}
void ExtraRotamersGeneric::ex1aro_exposed_sample_level( ExtraRotSample value ) {
	ex1aro_exposed_sample_level_ = value;
}
void ExtraRotamersGeneric::ex2aro_exposed_sample_level( ExtraRotSample value ) {
	ex2aro_exposed_sample_level_ = value;
}
void ExtraRotamersGeneric::exdna_sample_level( ExtraRotSample value ) {
	exdna_sample_level_ = value;
}
void ExtraRotamersGeneric::extrachi_cutoff( Size value ) {
	extrachi_cutoff_ = value;
}

/// BEGIN ReadResfile

ReadResfile::ReadResfile() : parent()
{
}

ReadResfile::ReadResfile( std::string const & filename ) :
	parent(),
	resfile_filename_( filename )
{}

ReadResfile::~ReadResfile() {}

TaskOperationOP ReadResfileCreator::create_task_operation() const
{
	return TaskOperationOP( new ReadResfile );
}

TaskOperationOP ReadResfile::clone() const
{
	return TaskOperationOP( new ReadResfile( *this ) );
}

void
ReadResfile::apply( pose::Pose const & pose, PackerTask & task ) const
{
	using namespace basic::resource_manager;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( resfile_filename_.empty() ) {
		/// only apply the read-resfile command if a resfile has been supplied, either through the
		/// resource manager, or through the command line
		if ( ResourceManager::get_instance()->has_option( packing::resfile ) ||  option[ packing::resfile ].user() ) {
			parse_resfile(pose, task, ResourceManager::get_instance()->get_option( packing::resfile )[ 1 ] );
		} /// else -- do not change the input PackerTask at all. Noop.
	} else {
		parse_resfile(pose, task, resfile_filename_ );
	}
}

/// @brief Assign the filename from the ResourceManager, if a resfile has been assigned for the
/// current job, and fall back on the options system, if a resfile has not been assigned.
void
ReadResfile::default_filename()
{
	using namespace basic::resource_manager;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( ResourceManager::get_instance()->has_option( packing::resfile ) ||  option[ packing::resfile ].user() ) {
		resfile_filename_ = ResourceManager::get_instance()->get_option( packing::resfile )[ 1 ];
	} else {
		resfile_filename_ = "";
	}
}

void
ReadResfile::filename( std::string const & filename )
{
	resfile_filename_ = filename;
}

std::string const & ReadResfile::filename() const
{
	return resfile_filename_;
}

void
ReadResfile::parse_tag( TagCOP tag , DataMap & )
{
	if ( tag->hasOption("filename") ) resfile_filename_ = tag->getOption<std::string>("filename");
	// special case: if "COMMANDLINE" string specified, use commandline option setting.
	// This is wholy unneccessary of course.  In the absence of a specified filename, the command line
	// will be read from, anyways.
	// if no filename is given, then the ReadResfile command will read either from the ResourceManager
	// or from the packing::resfile option on the command line.
	if ( resfile_filename_ == "COMMANDLINE" ) default_filename();
}

void
ReadResfile::parse_def( utility::lua::LuaObject const & def) {
	if( def["filename"] ) resfile_filename_ = def["filename"].to<std::string>();
	// special case: if "COMMANDLINE" string specified, use commandline option setting
	if ( resfile_filename_ == "COMMANDLINE" ) default_filename();
}

/// BEGIN ReadResfileAndObeyLengthEvents
ReadResfileAndObeyLengthEvents::ReadResfileAndObeyLengthEvents() :
	parent(),
	apply_default_commands_to_inserts_(false)
{}

ReadResfileAndObeyLengthEvents::ReadResfileAndObeyLengthEvents( std::string const & filename ) :
	parent( filename ),
	apply_default_commands_to_inserts_(false)
 {}

ReadResfileAndObeyLengthEvents::~ReadResfileAndObeyLengthEvents(){}

TaskOperationOP ReadResfileAndObeyLengthEventsCreator::create_task_operation() const
{
	return TaskOperationOP( new ReadResfileAndObeyLengthEvents );
}

TaskOperationOP ReadResfileAndObeyLengthEvents::clone() const
{
	return TaskOperationOP( new ReadResfileAndObeyLengthEvents( *this ) );
}


/// @details IMPORTANT: only use this if any length changes are
/// not reflected in the pose's pdb info, such as seems to be the case after vlb
/// not quite certain on the ideal approach yet. it's prolly best
/// to parse the resfile, and then apply the ResfileCommands to the
/// remapped residues. this necessitates getting the ResfileContents.
/// the code under 2. here is some duplication of ResfileReader::parse_resfile/parse_resfile_string,
/// ideally this and ResfileReader should be refactored a bit
void
ReadResfileAndObeyLengthEvents::apply(
	pose::Pose const & pose,
	PackerTask & ptask ) const
{

	//1. get the length change that the pose was exposed to
	//safeguard
	if( !pose.observer_cache().has( pose::datacache::CacheableObserverType::LENGTH_EVENT_COLLECTOR ) ){
		parent::apply( pose, ptask );
		return;
	}
	pose::datacache::CacheableObserverCOP len_obs = pose.observer_cache().get_const_ptr( pose::datacache::CacheableObserverType::LENGTH_EVENT_COLLECTOR );
	pose::datacache::LengthEventCollectorCOP lencollect( utility::pointer::static_pointer_cast< pose::datacache::LengthEventCollector const >( len_obs ) );

	utility::vector1< core::conformation::signals::LengthEvent > const & events( lencollect->events() );
	utility::vector1< core::id::SequenceMapping > smaps;
	for( Size i =1; i <= events.size(); ++i ){
		smaps.push_back( core::id::SequenceMapping( events[i] ) );
	}
	core::id::SequenceMappingOP fullsmap( core::id::combine_sequence_mappings( smaps ) );
	fullsmap->reverse();

	//2. get the resfile contents, we should probably
	//refactor ResfileReader a little bit to replace the
	//following block by one call
	std::string resfile_string;
	utility::io::izstream file( this->filename() );
	if (!file ) utility_exit_with_message( "Cannot open file " + this->filename() );
	utility::slurp( file, resfile_string );
	std::istringstream resfile(resfile_string);
	ResfileContents contents( pose, resfile );


	//3. apply the ResfileCommands to the remapped seqpos
	for ( Size ii = 1; ii <= ptask.total_residue(); ++ii ) {
		Size ii_resfile_seqpos( (*fullsmap)[ii] );

		if( !ii_resfile_seqpos && !apply_default_commands_to_inserts_ ) continue;

		std::list< ResfileCommandCOP > const & ii_command_list( this->resfile_commands( ii_resfile_seqpos, contents, ptask) );

		for ( std::list< ResfileCommandCOP >::const_iterator
					iter = ii_command_list.begin(), iter_end = ii_command_list.end();
					iter != iter_end; ++iter ) {
			(*iter)->residue_action( ptask, ii );
		}
	} //loop over all task residues

} //ReadResfileAndObeyLengthEvents::apply(

void
ReadResfileAndObeyLengthEvents::parse_tag( TagCOP tag , DataMap & datamap )
{
	parent::parse_tag( tag, datamap );
	if ( tag->hasOption("default_commands_for_inserts") )
		apply_default_commands_to_inserts_ = tag->getOption<bool>("default_commands_for_inserts",1);
}

/// @details note: this function will return default commands for
/// resfile_seqpos == 0
std::list< ResfileCommandCOP > const &
ReadResfileAndObeyLengthEvents::resfile_commands(
	core::Size const resfile_seqpos,
	ResfileContents const & contents,
	PackerTask const & ptask ) const
{
	if( (resfile_seqpos == 0) || ( resfile_seqpos > ptask.total_residue()) ){
		return contents.default_commands();
	}
	return (contents.specialized_commands_exist_for_residue( resfile_seqpos ) ?
		contents.commands_for_residue( resfile_seqpos ) : contents.default_commands() );
}


/// BEGIN SetRotamerCouplings

SetRotamerCouplings::SetRotamerCouplings()
{}

SetRotamerCouplings::~SetRotamerCouplings()
{}

SetRotamerCouplings::SetRotamerCouplings( SetRotamerCouplings const & src )
:
	parent(),
	rotamer_couplings_( src.rotamer_couplings_ )
{}

SetRotamerCouplings const &
SetRotamerCouplings::operator = ( SetRotamerCouplings const & rhs )
{
	rotamer_couplings_ = rhs.rotamer_couplings_;
	return *this;
}

TaskOperationOP SetRotamerCouplingsCreator::create_task_operation() const
{
	return TaskOperationOP( new SetRotamerCouplings );
}

TaskOperationOP SetRotamerCouplings::clone() const
{
	return TaskOperationOP( new SetRotamerCouplings( *this ) );
}

void
SetRotamerCouplings::apply( pose::Pose const &, PackerTask & task ) const
{
	task.rotamer_couplings( rotamer_couplings_ );
}

void
SetRotamerCouplings::set_couplings( rotamer_set::RotamerCouplingsOP couplings )
{
	rotamer_couplings_ = couplings;
}

/// BEGIN SetRotamerLinks

SetRotamerLinks::SetRotamerLinks()
{}

SetRotamerLinks::~SetRotamerLinks()
{}

SetRotamerLinks::SetRotamerLinks( SetRotamerLinks const & src )
:
	parent(),
	rotamer_links_( src.rotamer_links_ )
{}

SetRotamerLinks const &
SetRotamerLinks::operator = ( SetRotamerLinks const & rhs )
{
	rotamer_links_ = rhs.rotamer_links_;
	return *this;
}

TaskOperationOP SetRotamerLinksCreator::create_task_operation() const
{
	return TaskOperationOP( new SetRotamerLinks );
}

TaskOperationOP SetRotamerLinks::clone() const
{
	return TaskOperationOP( new SetRotamerLinks( *this ) );
}

void
SetRotamerLinks::apply( pose::Pose const &, PackerTask & task ) const
{
	task.rotamer_links( rotamer_links_ );
}

void
SetRotamerLinks::set_links( rotamer_set::RotamerLinksOP links )
{
	rotamer_links_ = links;
}

/// BEGIN AppendRotamer

AppendRotamer::AppendRotamer()
	: rotamer_operation_(/* 0 */)
{}

AppendRotamer::~AppendRotamer()
{}

AppendRotamer::AppendRotamer( rotamer_set::RotamerOperationOP rotamer_operation )
 : rotamer_operation_( rotamer_operation )
{}

AppendRotamer::AppendRotamer( AppendRotamer const & src )
: parent(), rotamer_operation_( src.rotamer_operation_ )
{}

TaskOperationOP AppendRotamerCreator::create_task_operation() const
{
	return TaskOperationOP( new AppendRotamer );
}

TaskOperationOP AppendRotamer::clone() const
{
	return TaskOperationOP( new AppendRotamer( *this ) );
}

void
AppendRotamer::apply( pose::Pose const &, PackerTask & task ) const
{
	task.append_rotamer_operation( rotamer_operation_ );
}

void
AppendRotamer::set_rotamer_operation(
	rotamer_set::RotamerOperationOP rotamer_operation
)
{
	rotamer_operation_ = rotamer_operation;
}


/// BEGIN AppendRotamerSet

AppendRotamerSet::AppendRotamerSet()
	: rotamer_set_operation_(/* 0 */)
{}

AppendRotamerSet::~AppendRotamerSet()
{}

AppendRotamerSet::AppendRotamerSet( rotamer_set::RotamerSetOperationOP rotamer_set_operation )
 : rotamer_set_operation_( rotamer_set_operation )
{}

AppendRotamerSet::AppendRotamerSet( AppendRotamerSet const & src )
: parent(), rotamer_set_operation_( src.rotamer_set_operation_ )
{}

TaskOperationOP AppendRotamerSetCreator::create_task_operation() const
{
	return TaskOperationOP( new AppendRotamerSet );
}

TaskOperationOP AppendRotamerSet::clone() const
{
	return TaskOperationOP( new AppendRotamerSet( *this ) );
}

void
AppendRotamerSet::apply( pose::Pose const &, PackerTask & task ) const
{
	task.append_rotamerset_operation( rotamer_set_operation_ );
}

void
AppendRotamerSet::set_rotamer_set_operation(
	rotamer_set::RotamerSetOperationOP rotamer_set_operation
)
{
	rotamer_set_operation_ = rotamer_set_operation;
}


/// BEGIN AppendResidueRotamerSet

AppendResidueRotamerSet::AppendResidueRotamerSet()
	: parent(),
	  resnum_(0),
	  rotamer_set_operation_(/* 0 */)
{}

AppendResidueRotamerSet::~AppendResidueRotamerSet()
{}

AppendResidueRotamerSet::AppendResidueRotamerSet( core::Size resnum,
		rotamer_set::RotamerSetOperationOP rotamer_set_operation )
 : resnum_(resnum),
   rotamer_set_operation_( rotamer_set_operation )
{}

AppendResidueRotamerSet::AppendResidueRotamerSet( AppendResidueRotamerSet const & src )
	: parent(src),
		resnum_( src.resnum_ ),
		rotamer_set_operation_( src.rotamer_set_operation_ )
{}

TaskOperationOP AppendResidueRotamerSetCreator::create_task_operation() const
{
	return TaskOperationOP( new AppendResidueRotamerSet );
}

TaskOperationOP AppendResidueRotamerSet::clone() const
{
	return TaskOperationOP( new AppendResidueRotamerSet( *this ) );
}

void
AppendResidueRotamerSet::apply( pose::Pose const &, PackerTask & task ) const
{
	task.nonconst_residue_task(resnum_).append_rotamerset_operation(rotamer_set_operation_);
}

void
AppendResidueRotamerSet::set_resnum( core::Size resnum )
{
	resnum_ = resnum;
}

void
AppendResidueRotamerSet::set_rotamer_set_operation( rotamer_set::RotamerSetOperationOP rotamer_set_operation )
{
	rotamer_set_operation_ = rotamer_set_operation;
}


/// BEGIN PreserveCBeta

PreserveCBeta::~PreserveCBeta() {}

TaskOperationOP PreserveCBetaCreator::create_task_operation() const
{
	return TaskOperationOP( new PreserveCBeta );
}

TaskOperationOP PreserveCBeta::clone() const
{
	return TaskOperationOP( new PreserveCBeta( *this ) );
}

void
PreserveCBeta::apply( pose::Pose const &, PackerTask & task ) const
{
	task.or_preserve_c_beta( true );
}

/// BEGIN PreventRepacking
PreventRepacking::~PreventRepacking() {}

TaskOperationOP PreventRepackingCreator::create_task_operation() const
{
	return TaskOperationOP( new PreventRepacking );
}

TaskOperationOP PreventRepacking::clone() const
{
	return TaskOperationOP( new PreventRepacking( *this ) );
}

void
PreventRepacking::apply( pose::Pose const & pose, PackerTask & task ) const
{
	utility::vector1<core::Size> const res = core::pose::get_resnum_list_ordered( residue_selection_, pose );
	utility::vector1<core::Size> residues_to_prevent = residues_to_prevent_;
	residues_to_prevent.insert(residues_to_prevent.end(),res.begin(),res.end());

	for(utility::vector1< core::Size >::const_iterator it(residues_to_prevent.begin()), end(residues_to_prevent.end());
			it != end; ++it)
		{task.nonconst_residue_task(*it).prevent_repacking();}
	return;
}

void
PreventRepacking::include_residue( core::Size resid ) { residues_to_prevent_.push_back(resid); }

void
PreventRepacking::clear() { residues_to_prevent_.clear(); }

void
PreventRepacking::parse_tag( TagCOP tag , DataMap & )
{
	residue_selection_ = tag->getOption<std::string>("resnum","0");
}


// BEGIN RestrictYSDesign
RestrictYSDesign::~RestrictYSDesign() {}
RestrictYSDesign::RestrictYSDesign() : gly_switch_( false )  {}
RestrictYSDesign::RestrictYSDesign( RestrictYSDesign const & src ) :
	parent(), YSresids_( src.YSresids_ ), gly_switch_( src.gly_switch_ )
{}
RestrictYSDesign::RestrictYSDesign( utility::vector1< core::Size > const & resids ) {
	gly_switch_ = false;
	YSresids_ = resids;
}

void
RestrictYSDesign::apply( pose::Pose const &, PackerTask & task ) const {
	utility::vector1<bool> restrict_to_aa( 20, false );
	for( utility::vector1<core::Size>::const_iterator res_it=YSresids_.begin(); res_it!=YSresids_.end(); ++res_it ) {
		if( gly_switch_ ) restrict_to_aa[chemical::aa_from_name( "GLY" )] = true;
		restrict_to_aa[chemical::aa_from_name( "TYR" )] = true;
		restrict_to_aa[chemical::aa_from_name( "SER" )] = true;
		task.nonconst_residue_task(*res_it).restrict_absent_canonical_aas( restrict_to_aa );
	}
}

TaskOperationOP RestrictYSDesignCreator::create_task_operation() const
{
	return TaskOperationOP( new RestrictYSDesign );
}

TaskOperationOP RestrictYSDesign::clone() const {
	return TaskOperationOP( new RestrictYSDesign( *this ) );
}

void
RestrictYSDesign::include_resid( core::Size const resid ) { YSresids_.push_back(resid); }

void
RestrictYSDesign::include_gly( bool const gly ) { gly_switch_ = gly; }


//////////////////////////////////////////////////////
// This class could easily be expanded to handle sample_level, etc.
// Someone has probably already written this, and should replace this
// dinky thing.
ExtraRotamers::ExtraRotamers() :
	resid_( 0 ), // default apply to all residues
	chi_( 0 ),
	level_( 0 )
{}

ExtraRotamers::ExtraRotamers( core::Size const resid, core::Size const chi, core::Size level ):
	resid_( resid ),
	chi_( chi ),
	level_( level )
{}

ExtraRotamers::~ExtraRotamers() {}

TaskOperationOP ExtraRotamersCreator::create_task_operation() const
{
	return TaskOperationOP( new ExtraRotamers );
}

TaskOperationOP ExtraRotamers::clone() const
{
	return TaskOperationOP( new ExtraRotamers( *this ) );
}

void
ExtraRotamers::apply( core::pose::Pose const & p, PackerTask & task ) const
{
	if ( resid_ != 0 ) {
		ResidueLevelTask & restask( task.nonconst_residue_task( resid_ ) );
		if( chi_ == 1 ) restask.or_ex1( ExtraRotSample( level_ ) );
		if( chi_ == 2 ) restask.or_ex2( ExtraRotSample( level_ ) );
		if( chi_ == 3 ) restask.or_ex3( ExtraRotSample( level_ ) );
		if( chi_ == 4 ) restask.or_ex4( ExtraRotSample( level_ ) );
	} else {
		// apply to all residues
		TR << "Enabling extra rotamers for chi " << chi_ << " at all positions" << std::endl;
		for ( Size ii = 1; ii <= p.total_residue(); ++ii ) {
			ResidueLevelTask & restask( task.nonconst_residue_task( ii ) );
			if( chi_ == 1 ) restask.or_ex1( ExtraRotSample( level_ ) );
			if( chi_ == 2 ) restask.or_ex2( ExtraRotSample( level_ ) );
			if( chi_ == 3 ) restask.or_ex3( ExtraRotSample( level_ ) );
			if( chi_ == 4 ) restask.or_ex4( ExtraRotSample( level_ ) );
		}
	}
}

void ExtraRotamers::parse_tag( TagCOP tag , DataMap & )
{
	if ( tag->hasOption("resid") ) {
		resid_ = tag->getOption< core::Size >( "resid" );
	} else {
		resid_ = 0; // apply to all residues
	}

	if ( ! tag->hasOption("chi")  ) {
		utility_exit_with_message("ExtraRotamers Task Operation requires the chi option");
	}
	chi_ = tag->getOption< core::Size >("chi");
	if ( chi_ > 4 ) {
		utility_exit_with_message("ExtraRotamers Task Operation given a value for chi outside of legal range 1-4.  Given value of " +  utility::to_string(chi_) );
	}
	if ( tag->hasOption("level") ) {
		level_ = tag->getOption< core::Size >("level");
		if ( level_ > core::Size( ExtraRotSampleCardinality ) ) {
			utility_exit_with_message( "ExtraRotamers Task Operation gien a value for level outside of legal range 1-" + utility::to_string( core::Size( ExtraRotSampleCardinality ) ) + ".  Given value of " + utility::to_string( level_ ) );

		}
	}
}

//////////////////////////////////////////////////////
// This class could easily be expanded ...
// Someone has probably already written this, and should replace this
// dinky thing.
ExtraChiCutoff::ExtraChiCutoff() {}

ExtraChiCutoff::ExtraChiCutoff( core::Size const resid, core::Size const extrachi_cutoff):
	resid_( resid ),
	extrachi_cutoff_( extrachi_cutoff )
{}

ExtraChiCutoff::~ExtraChiCutoff() {}

TaskOperationOP ExtraChiCutoffCreator::create_task_operation() const
{
	return TaskOperationOP( new ExtraChiCutoff );
}

TaskOperationOP ExtraChiCutoff::clone() const
{
	return TaskOperationOP( new ExtraChiCutoff( *this ) );
}

void
ExtraChiCutoff::apply( core::pose::Pose const & p, PackerTask & task ) const
{
	if ( resid_ != 0 ) {
		task.nonconst_residue_task( resid_ ).and_extrachi_cutoff( extrachi_cutoff_ );
	} else {
		//appy to all residues
		TR << "Enabling extrachi cutoff at all positions" << std::endl;
		for ( Size ii = 1; ii <= p.total_residue(); ++ii ) {
			task.nonconst_residue_task( ii ).and_extrachi_cutoff( extrachi_cutoff_ );
		}
	}
}

void ExtraChiCutoff::parse_tag( TagCOP tag , DataMap & )
{
	if ( tag->hasOption("resid") ) {
		resid_ = tag->getOption< core::Size >( "resid" );
	} else {
		resid_ = 0; // apply to all residues
	}

	if ( ! tag->hasOption("extrachi_cutoff")  ) {
		utility_exit_with_message("ExtraChiCutoff Task Operation requires the extrachi_cutoff option");
	}
	extrachi_cutoff_ = tag->getOption< core::Size >("extrachi_cutoff");
}


} //namespace operation
} //namespace task
} //namespace pack
} //namespace core
