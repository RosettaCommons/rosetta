// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/operation/TaskOperations.cc
/// @brief
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Unit Headers
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/TaskOperationCreators.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

// Project Headers
#include <core/id/SequenceMapping.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/rotamer_set/RotamerCouplings.hh>
#include <core/pack/rotamer_set/RotamerLinks.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/selection.hh>
#include <core/pose/datacache/CacheableObserverType.hh>
#include <core/pose/datacache/cacheable_observers.hh>
#include <core/pose/datacache/ObserverCache.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility>
#include <utility/exit.hh>
#include <utility/io/izstream.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/keys/OptionKeyList.hh>
#include <utility/keys/VariantKey.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

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

static basic::Tracer TR( "core.pack.task.operation.TaskOperations", t_info );

using namespace utility::tag;

/// BEGIN RestrictToRepacking

RestrictToRepacking::~RestrictToRepacking() = default;

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

std::string RestrictToRepacking::keyname() { return "RestrictToRepacking"; }

void RestrictToRepacking::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	task_op_schema_empty( xsd, keyname(), "Only allow residues to repack. No design at any position." );
}

TaskOperationOP RestrictToRepackingCreator::create_task_operation() const
{
	return TaskOperationOP( new RestrictToRepacking );
}

void RestrictToRepackingCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RestrictToRepacking::provide_xml_schema( xsd );
}

std::string RestrictToRepackingCreator::keyname() const {
	return RestrictToRepacking::keyname();
}

/// BEGIN RestrictResidueToRepacking
RestrictResidueToRepacking::~RestrictResidueToRepacking() = default;

TaskOperationOP RestrictResidueToRepacking::clone() const
{
	return TaskOperationOP( new RestrictResidueToRepacking( *this ) );
}

void
RestrictResidueToRepacking::apply( pose::Pose const &, PackerTask & task ) const
{
	for ( unsigned long it : residues_to_restrict_to_repacking_ ) {
		//debug_assert( *it ); //begin/end can't return NULL anyway?
		task.nonconst_residue_task(it).restrict_to_repacking();
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

std::string RestrictResidueToRepacking::keyname() { return "RestrictResidueToRepacking"; }

void RestrictResidueToRepacking::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attributes;
	attributes + XMLSchemaAttribute::required_attribute(  "resnum", xsct_positive_integer, "Single residue number." );
	task_op_schema_w_attributes( xsd, keyname(), attributes , "Restrict a single residue to repacking. No design." );
}

TaskOperationOP RestrictResidueToRepackingCreator::create_task_operation() const
{
	return TaskOperationOP( new RestrictResidueToRepacking );
}

std::string RestrictResidueToRepackingCreator::keyname() const {
	return RestrictResidueToRepacking::keyname();
}

void RestrictResidueToRepackingCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	return RestrictResidueToRepacking::provide_xml_schema( xsd );
}


/// BEGIN RestrictAbsentCanonicalAAS
RestrictAbsentCanonicalAAS::RestrictAbsentCanonicalAAS()
: parent(),
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

RestrictAbsentCanonicalAAS::~RestrictAbsentCanonicalAAS()= default;


TaskOperationOP RestrictAbsentCanonicalAAS::clone() const
{
	return TaskOperationOP( new RestrictAbsentCanonicalAAS( *this ) );
}

void
RestrictAbsentCanonicalAAS::apply( pose::Pose const &, PackerTask & task ) const
{
	if ( resid_ == 0 ) { // restrict all residues
		for ( core::Size i( 1 ); i <= task.total_residue(); ++i ) {
			task.nonconst_residue_task( i ).restrict_absent_canonical_aas( keep_aas_ );
		}
	} else {
		task.nonconst_residue_task( resid_ ).restrict_absent_canonical_aas( keep_aas_ );
	}
}

void
RestrictAbsentCanonicalAAS::keep_aas( std::string const & keep )
{
	using namespace chemical;
	runtime_assert( keep_aas_.size() == num_canonical_aas );
	utility::vector1< bool > canonical_aas_to_keep( num_canonical_aas, false );
	for ( char const c : keep ) {
		if ( oneletter_code_specifies_aa( c ) ) {
			//std::cout << "Keeping amino acid " << c << std::endl;
			canonical_aas_to_keep[ aa_from_oneletter_code( c ) ] = true;
		} else {
			TR.Fatal << "aa letter " << c << " does not not correspond to a canonical AA"<<std::endl;
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
	// std::cout << " keeping " << ii << " " << " ? " << keep_aas_[ ii ]  << std::endl;
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

std::string RestrictAbsentCanonicalAAS::keyname() { return "RestrictAbsentCanonicalAAS"; }

void RestrictAbsentCanonicalAAS::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::attribute_w_default(  "resnum", xsct_non_negative_integer, "Restrict design to user-specified residues. If =0, all residues.",  "0"  )
		+ XMLSchemaAttribute( "keep_aas", utility::tag::xs_string , "Canonical amino acids to keep." );
	task_op_schema_w_attributes( xsd, keyname(), attributes, "Restrict design to user-specified residues. If resnum is left as 0, the restriction will apply throughout the pose." );
}

TaskOperationOP RestrictAbsentCanonicalAASCreator::create_task_operation() const
{
	return TaskOperationOP( new RestrictAbsentCanonicalAAS );
}

std::string RestrictAbsentCanonicalAASCreator::keyname() const {
	return RestrictAbsentCanonicalAAS::keyname();
}

void RestrictAbsentCanonicalAASCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RestrictAbsentCanonicalAAS::provide_xml_schema( xsd );
}



//BEGIN DisallowIfNonnative
DisallowIfNonnative::DisallowIfNonnative():
	disallowed_aas_ ( chemical::num_canonical_aas, false ),
	allowed_aas_( invert_vector( disallowed_aas_ ) )
{}

DisallowIfNonnative::DisallowIfNonnative( utility::vector1< bool > const & disallowed_aas):
	disallowed_aas_( disallowed_aas ),
	allowed_aas_( invert_vector(disallowed_aas) )
{}

DisallowIfNonnative::DisallowIfNonnative( utility::vector1< bool > const & disallowed_aas, utility::vector1<core::Size> const & res_selection):
	residue_selection_( res_selection ),
	disallowed_aas_( disallowed_aas ),
	allowed_aas_( invert_vector(disallowed_aas) )
{}

DisallowIfNonnative::~DisallowIfNonnative()= default;

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
	for ( core::Size ii=1; ii<=disallowed_aas_.size(); ii++ ) {
		inverted_vec.push_back( ! disallowed_aas[ii] );
	}
	return inverted_vec;
}

void DisallowIfNonnative::apply( pose::Pose const &, PackerTask & task ) const
{
	runtime_assert( allowed_aas_.size() == chemical::num_canonical_aas );
	if ( residue_selection_.empty() ) { //if no residue defined then do all residues
		for ( core::Size ii = 1; ii<= task.total_residue(); ii++ ) {
			task.nonconst_residue_task( ii ).restrict_nonnative_canonical_aas( allowed_aas_ );
		}
	} else {  //if a residue is defined then do only on that residue
		for ( core::Size jj = 1; jj <= residue_selection_.size(); jj++ ) {
			task.nonconst_residue_task( residue_selection_[jj] ).restrict_nonnative_canonical_aas( allowed_aas_ );
		}
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
	for ( char it : aa_string ) {
		if ( oneletter_code_specifies_aa( it ) ) {
			aa_vector[ aa_from_oneletter_code( it ) ] = true;
		} else {
			std::ostringstream os;
			os << "aa letter " << it << " does not not correspond to a canonical AA";
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
	if ( resid != 0 ) residue_selection_.push_back( resid );
}
void DisallowIfNonnative::restrict_to_residue( utility::vector1< core::Size > const & residues){
	for ( core::Size ii=1; ii<=residues.size(); ii++ ) {
		residue_selection_.push_back( residues[ii] );
	}
}

void DisallowIfNonnative::parse_tag( TagCOP tag , DataMap & )
{
	restrict_to_residue( tag->getOption< core::Size >( "resnum", 0 ) );
	disallow_aas( tag->getOption< std::string >( "disallow_aas" ) );
}

std::string DisallowIfNonnative::keyname() { return "DisallowIfNonnative"; }

void DisallowIfNonnative::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::attribute_w_default(  "resnum", xsct_non_negative_integer, "If resnum is left as 0, the restriction will apply throughout the pose.",  "0"  )
		+ XMLSchemaAttribute( "disallow_aas", xs_string , "Disallow_aas takes a string of one letter amino acid codes, no separation needed. For example disallow_aas=GCP would prevent Gly, Cys, and Pro from being designed unless they were the native amino acid at a position." );
	task_op_schema_w_attributes( xsd, keyname(), attributes, "Restrict design to not include a residue as an possibility in the task at a position unless it is the starting residue." );
}

TaskOperationOP DisallowIfNonnativeCreator::create_task_operation() const
{
	return TaskOperationOP( new DisallowIfNonnative );
}

std::string DisallowIfNonnativeCreator::keyname() const {
	return DisallowIfNonnative::keyname();
}

void DisallowIfNonnativeCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DisallowIfNonnative::provide_xml_schema( xsd );
}


//BEGIN RotamerExplosion
RotamerExplosion::RotamerExplosion()= default;

RotamerExplosion::RotamerExplosion( core::Size const resid, ExtraRotSample const sample_level, core::Size const chi ) :
	resid_( resid ),
	chi_( chi ),
	sample_level_( sample_level )
{}

RotamerExplosion::~RotamerExplosion() = default;

TaskOperationOP RotamerExplosion::clone() const
{
	return TaskOperationOP( new RotamerExplosion( *this ) );
}

void
RotamerExplosion::apply( core::pose::Pose const &, PackerTask & task ) const
{
	ResidueLevelTask & restask( task.nonconst_residue_task( resid_ ) );
	if ( chi_ > 0 ) restask.or_ex1_sample_level( sample_level_ );
	if ( chi_ > 1 ) restask.or_ex2_sample_level( sample_level_ );
	if ( chi_ > 2 ) restask.or_ex3_sample_level( sample_level_ );
	if ( chi_ > 3 ) restask.or_ex4_sample_level( sample_level_ );
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

std::string RotamerExplosion::keyname() { return "RotamerExplosionCreator"; }

void RotamerExplosion::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::required_attribute( "resnum", xsct_positive_integer , "Residue number." )
		+ XMLSchemaAttribute::required_attribute( "chi",    xsct_non_negative_integer , "Chi level 0-3?" );
	task_op_schema_w_attributes( xsd, keyname(), attributes, "Sample residue chi angles much more finely during packing. Currently hardcoded to use three 1/3 step standard deviation." );
}

TaskOperationOP RotamerExplosionCreator::create_task_operation() const
{
	return TaskOperationOP( new RotamerExplosion );
}

std::string RotamerExplosionCreator::keyname() const { return RotamerExplosion::keyname(); }

void RotamerExplosionCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RotamerExplosion::provide_xml_schema( xsd );
}



/// BEGIN InitializeFromCommandline

InitializeFromCommandline::~InitializeFromCommandline() = default;

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

std::string InitializeFromCommandline::keyname() { return "InitializeFromCommandline"; }

void InitializeFromCommandline::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	task_op_schema_empty( xsd, keyname(), "Reads commandline options. For example, -ex1 -ex2 (does not read resfile from command line options) This taskoperation will complain about an unimplemented method, but you can safely ignore the message." );
}

TaskOperationOP InitializeFromCommandlineCreator::create_task_operation() const
{
	return TaskOperationOP( new InitializeFromCommandline );
}

std::string InitializeFromCommandlineCreator::keyname() const {
	return InitializeFromCommandline::keyname();
}

void InitializeFromCommandlineCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	InitializeFromCommandline::provide_xml_schema( xsd );
}

/// BEGIN UseMultiCoolAnnealer

UseMultiCoolAnnealer::~UseMultiCoolAnnealer() = default;

UseMultiCoolAnnealer::UseMultiCoolAnnealer()
: parent(),
	states_(10)
{
}

UseMultiCoolAnnealer::UseMultiCoolAnnealer( core::Size states )
: parent()
{
	set_states( states );
}

TaskOperationOP UseMultiCoolAnnealer::clone() const
{
	return TaskOperationOP( new UseMultiCoolAnnealer( *this ) );
}

void
UseMultiCoolAnnealer::apply( pose::Pose const &, PackerTask & task ) const
{
	task.or_multi_cool_annealer( true );
	task.increase_multi_cool_annealer_history_size( states_ );
}

void
UseMultiCoolAnnealer::set_states( core::Size states ) { states_ = states; }

void
UseMultiCoolAnnealer::parse_tag( TagCOP tag, DataMap & )
{
	set_states( tag->getOption< core::Size >( "states", 10 ) );
}

std::string UseMultiCoolAnnealer::keyname() { return "UseMultiCoolAnnealer"; }

void UseMultiCoolAnnealer::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::required_attribute( "states", xsct_non_negative_integer, "Number of states" );
	task_op_schema_w_attributes( xsd, keyname(), attributes );
}

TaskOperationOP UseMultiCoolAnnealerCreator::create_task_operation() const
{
	return TaskOperationOP( new UseMultiCoolAnnealer );
}

std::string UseMultiCoolAnnealerCreator::keyname() const {
	return UseMultiCoolAnnealer::keyname();
}

void UseMultiCoolAnnealerCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	UseMultiCoolAnnealer::provide_xml_schema( xsd );
}


/// BEGIN InitializeFromOptionCollection

InitializeFromOptionCollection::InitializeFromOptionCollection() = default;

InitializeFromOptionCollection::InitializeFromOptionCollection( utility::options::OptionCollectionCOP options )
{
	options_ = options;
}

InitializeFromOptionCollection::InitializeFromOptionCollection(
	InitializeFromOptionCollection const & src
) :
	options_( src.options_ )
{}

InitializeFromOptionCollection::~InitializeFromOptionCollection() = default;

TaskOperationOP InitializeFromOptionCollection::clone() const
{
	return TaskOperationOP( new InitializeFromOptionCollection( *this ) );
}

void
InitializeFromOptionCollection::apply( pose::Pose const &, PackerTask & task ) const
{
	runtime_assert( options_ );
	task.initialize_from_options( *options_ );
}

void
InitializeFromOptionCollection::parse_tag( TagCOP tag, DataMap & datamap )
{
	//runtime_assert( datamap.has( "options"));
	std::string which_options = tag->getOption< std::string >( "option_collection", "job_options" );
	if ( ! datamap.has( "options", which_options ) ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Failed to find options named \"" + which_options + "\" in the InitializeFromOptionCollection task operation" );
	}
	options_ = datamap.get_ptr< utility::options::OptionCollection const >( "options", which_options );
}

void
InitializeFromOptionCollection::options( utility::options::OptionCollectionCOP options )
{
	options_ = options;
}


std::string InitializeFromOptionCollection::keyname() { return "InitializeFromOptionCollection"; }

void InitializeFromOptionCollection::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attributes;
	attributes + XMLSchemaAttribute::attribute_w_default(  "option_collection", xs_string, "Options in datamap.",  "job_options"  );
	task_op_schema_w_attributes( xsd, keyname(), attributes, "Initializes task operations from an option collection. Similar to InitializeFromCommandLine, but uses a specified option collection object instead of basic::option::options." );
}

TaskOperationOP InitializeFromOptionCollectionCreator::create_task_operation() const
{
	return TaskOperationOP( new InitializeFromOptionCollection );
}

std::string InitializeFromOptionCollectionCreator::keyname() const {
	return InitializeFromOptionCollection::keyname();
}

void InitializeFromOptionCollectionCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	InitializeFromOptionCollection::provide_xml_schema( xsd );
}

/// BEGIN InitializeExtraRotsFromCommandline

InitializeExtraRotsFromCommandline::~InitializeExtraRotsFromCommandline() = default;

TaskOperationOP InitializeExtraRotsFromCommandline::clone() const
{
	return TaskOperationOP( new InitializeExtraRotsFromCommandline( *this ) );
}

void
InitializeExtraRotsFromCommandline::apply( pose::Pose const &, PackerTask & task ) const
{
	task.initialize_extra_rotamer_flags_from_command_line();
}

std::string InitializeExtraRotsFromCommandline::keyname() { return "InitializeExtraRotsFromCommandline"; }

void InitializeExtraRotsFromCommandline::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	task_op_schema_empty( xsd, keyname(), "Creates extra rotamer operations from command line options. Under development and untested. Use at your own risk." );
}

TaskOperationOP InitializeExtraRotsFromCommandlineCreator::create_task_operation() const
{
	return TaskOperationOP( new InitializeExtraRotsFromCommandline );
}

std::string InitializeExtraRotsFromCommandlineCreator::keyname() const {
	return InitializeExtraRotsFromCommandline::keyname();
}


void InitializeExtraRotsFromCommandlineCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	InitializeExtraRotsFromCommandline::provide_xml_schema( xsd );
}


/// BEGIN IncludeCurrent

IncludeCurrent::~IncludeCurrent() = default;


TaskOperationOP IncludeCurrent::clone() const
{
	return TaskOperationOP( new IncludeCurrent( *this ) );
}

void
IncludeCurrent::apply( pose::Pose const &, PackerTask & task ) const
{
	task.or_include_current(true);
}

std::string IncludeCurrent::keyname() { return "IncludeCurrent"; }

void IncludeCurrent::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	task_op_schema_empty( xsd, keyname(), "Includes current rotamers (eg - from input pdb) in the rotamer set. These rotamers will be lost after a packing run, so they are only effective upon initial loading of a pdb!" );
}

TaskOperationOP IncludeCurrentCreator::create_task_operation() const
{
	return TaskOperationOP( new IncludeCurrent );
}

std::string IncludeCurrentCreator::keyname() const {
	return IncludeCurrent::keyname();
}

void IncludeCurrentCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	IncludeCurrent::provide_xml_schema( xsd );
}

/// BEGIN ExtraRotamersGeneric
ExtraRotamerSamplingData::ExtraRotamerSamplingData() :
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

ExtraRotamersGeneric::ExtraRotamersGeneric() :
	sampling_data_()
{}

ExtraRotamersGeneric::~ExtraRotamersGeneric() = default;


TaskOperationOP ExtraRotamersGeneric::clone() const
{
	return TaskOperationOP( new ExtraRotamersGeneric( *this ) );
}

void
ExtraRotamersGeneric::parse_tag( TagCOP tag , DataMap & )
{
	parse_rotamer_sampling_data( tag, sampling_data_ );
}

void
ExtraRotamersGeneric::apply( pose::Pose const &, PackerTask & task ) const
{
	for ( Size i=1; i <= task.total_residue(); ++i ) {
		ResidueLevelTask & res_task(task.nonconst_residue_task(i));
		set_rotamer_sampling_data_for_RLT( sampling_data_, res_task );
	}
}

void ExtraRotamersGeneric::ex1( bool value ) {
	sampling_data_.ex1_ = value;
}
void ExtraRotamersGeneric::ex2( bool value ) {
	sampling_data_.ex2_ = value;
}
void ExtraRotamersGeneric::ex3( bool value ) {
	sampling_data_.ex3_ = value;
}
void ExtraRotamersGeneric::ex4( bool value ) {
	sampling_data_.ex4_ = value;
}
void ExtraRotamersGeneric::ex1aro( bool value ) {
	sampling_data_.ex1aro_ = value;
}
void ExtraRotamersGeneric::ex2aro( bool value ) {
	sampling_data_.ex2aro_ = value;
}
void ExtraRotamersGeneric::ex1aro_exposed( bool value ) {
	sampling_data_.ex1aro_exposed_ = value;
}
void ExtraRotamersGeneric::ex2aro_exposed( bool value ) {
	sampling_data_.ex2aro_exposed_ = value;
}
void ExtraRotamersGeneric::ex1_sample_level( ExtraRotSample value ) {
	sampling_data_.ex1_sample_level_ = value;
}
void ExtraRotamersGeneric::ex2_sample_level( ExtraRotSample value ) {
	sampling_data_.ex2_sample_level_ = value;
}
void ExtraRotamersGeneric::ex3_sample_level( ExtraRotSample value ) {
	sampling_data_.ex3_sample_level_ = value;
}
void ExtraRotamersGeneric::ex4_sample_level( ExtraRotSample value ) {
	sampling_data_.ex4_sample_level_ = value;
}
void ExtraRotamersGeneric::ex1aro_sample_level( ExtraRotSample value ) {
	sampling_data_.ex1aro_sample_level_ = value;
}
void ExtraRotamersGeneric::ex2aro_sample_level( ExtraRotSample value ) {
	sampling_data_.ex2aro_sample_level_ = value;
}
void ExtraRotamersGeneric::ex1aro_exposed_sample_level( ExtraRotSample value ) {
	sampling_data_.ex1aro_exposed_sample_level_ = value;
}
void ExtraRotamersGeneric::ex2aro_exposed_sample_level( ExtraRotSample value ) {
	sampling_data_.ex2aro_exposed_sample_level_ = value;
}
void ExtraRotamersGeneric::exdna_sample_level( ExtraRotSample value ) {
	sampling_data_.exdna_sample_level_ = value;
}
void ExtraRotamersGeneric::extrachi_cutoff( Size value ) {
	sampling_data_.extrachi_cutoff_ = value;
}

ExtraRotamerSamplingData const &
ExtraRotamersGeneric::sampling_data() const
{
	return sampling_data_;
}

std::string ExtraRotamersGeneric::keyname() { return "ExtraRotamersGeneric"; }

void ExtraRotamersGeneric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	AttributeList attributes = rotamer_sampling_data_xml_schema_attributes( xsd );
	task_op_schema_w_attributes( xsd, keyname(), attributes, "During packing, extra rotamers can be used to increase sampling. Use this TaskOperation to specify for all residues at once what extra rotamers should be used." );
}

TaskOperationOP ExtraRotamersGenericCreator::create_task_operation() const
{
	return TaskOperationOP( new ExtraRotamersGeneric );
}

std::string ExtraRotamersGenericCreator::keyname() const {
	return ExtraRotamersGeneric::keyname();
}

void ExtraRotamersGenericCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ExtraRotamersGeneric::provide_xml_schema( xsd );
}

void parse_rotamer_sampling_data(
	utility::tag::TagCOP tag,
	ExtraRotamerSamplingData & sampling_data
)
{
	sampling_data.ex1_ = tag->getOption<bool>("ex1", false);
	sampling_data.ex2_ = tag->getOption<bool>("ex2", false);
	sampling_data.ex3_ = tag->getOption<bool>("ex3", false);
	sampling_data.ex4_ = tag->getOption<bool>("ex4", false);
	sampling_data.ex1aro_ = tag->getOption<bool>("ex1aro", false);
	sampling_data.ex2aro_ = tag->getOption<bool>("ex2aro", false);
	sampling_data.ex1aro_exposed_ = tag->getOption<bool>("ex1aro_exposed", false);
	sampling_data.ex2aro_exposed_ = tag->getOption<bool>("ex2aro_exposed", false);

	sampling_data.ex1_sample_level_ = static_cast<ExtraRotSample>(tag->getOption<Size>("ex1_sample_level", NO_EXTRA_CHI_SAMPLES));
	sampling_data.ex2_sample_level_ = static_cast<ExtraRotSample>(tag->getOption<Size>("ex2_sample_level", NO_EXTRA_CHI_SAMPLES));
	sampling_data.ex3_sample_level_ = static_cast<ExtraRotSample>(tag->getOption<Size>("ex3_sample_level", NO_EXTRA_CHI_SAMPLES));
	sampling_data.ex4_sample_level_ = static_cast<ExtraRotSample>(tag->getOption<Size>("ex4_sample_level", NO_EXTRA_CHI_SAMPLES));
	sampling_data.ex1aro_sample_level_ = static_cast<ExtraRotSample>(tag->getOption<Size>("ex1aro_sample_level", NO_EXTRA_CHI_SAMPLES));
	sampling_data.ex2aro_sample_level_ = static_cast<ExtraRotSample>(tag->getOption<Size>("ex2aro_sample_level", NO_EXTRA_CHI_SAMPLES));
	sampling_data.ex1aro_exposed_sample_level_ = static_cast<ExtraRotSample>(tag->getOption<Size>("ex1aro_exposed_sample_level", NO_EXTRA_CHI_SAMPLES));
	sampling_data.ex2aro_exposed_sample_level_ = static_cast<ExtraRotSample>(tag->getOption<Size>("ex2aro_exposed_sample_level", NO_EXTRA_CHI_SAMPLES));
	sampling_data.exdna_sample_level_ = static_cast<ExtraRotSample>(tag->getOption<Size>("exdna_sample_level", NO_EXTRA_CHI_SAMPLES));

	sampling_data.extrachi_cutoff_ = tag->getOption<Size>("extrachi_cutoff", EXTRACHI_CUTOFF_LIMIT);
}

void
define_extra_rotamers_sampling_level_restriction( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	XMLSchemaRestriction exchi_sample_level = integer_range_restriction( "exchi_sample_level", NO_EXTRA_CHI_SAMPLES, EX_SIX_QUARTER_STEP_STDDEVS ) ;
	xsd.add_top_level_element( exchi_sample_level );
}

AttributeList
rotamer_sampling_data_xml_schema_attributes( XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	define_extra_rotamers_sampling_level_restriction( xsd );

	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::attribute_w_default(  "ex1", xsct_rosetta_bool, "Use extra chi1 sub-rotamers for all residues that pass the extrachi_cutof",  "0"  )
		+ XMLSchemaAttribute::attribute_w_default(  "ex2", xsct_rosetta_bool, "Use extra chi2 sub-rotamers for all residues that pass the extrachi_cutoff",  "0"  )
		+ XMLSchemaAttribute::attribute_w_default(  "ex3", xsct_rosetta_bool, "Use extra chi3 sub-rotamers for all residues that pass the extrachi_cutoff",  "0"  )
		+ XMLSchemaAttribute::attribute_w_default(  "ex4", xsct_rosetta_bool, "Use extra chi4 sub-rotamers for all residues that pass the extrachi_cutoff",  "0"  )
		+ XMLSchemaAttribute::attribute_w_default(  "ex1aro", xsct_rosetta_bool, "Use extra chi1 sub-rotamers for aromatic residues that pass the extrachi_cutoff",  "0"  )
		+ XMLSchemaAttribute::attribute_w_default(  "ex2aro", xsct_rosetta_bool, "Use extra chi2 sub-rotamers for aromatic residues that pass the extrachi_cutoff",  "0"  )
		+ XMLSchemaAttribute::attribute_w_default(  "ex1aro_exposed", xsct_rosetta_bool, "Use extra chi1 sub-rotamers for all aromatic residues",  "0"  )
		+ XMLSchemaAttribute::attribute_w_default(  "ex2aro_exposed", xsct_rosetta_bool, "Use extra chi2 sub-rotamers for all aromatic residues",  "0"  )

		+ XMLSchemaAttribute::attribute_w_default(  "ex1_sample_level", "exchi_sample_level", "Sampling level for chi1 sampling; all include the current chi dihedral. There are currently 8 options: NO_EXTRA_CHI_SAMPLES (0) (original only); EX_ONE_STDDEV (1) (default, +- one sd, 3 samples); EX_ONE_HALF_STEP_STDDEV (2)(+- 0.5 sd, 3 samples); EX_TWO_FULL_STEP_STDDEVS(3) ( +- 1 and 2 sd, 5 samples); EX_TWO_HALF_STEP_STDDEVS(4) ( +- 0.5 and +- 1 sd, 5 samples); EX_FOUR_HALF_STEP_STDDEVS(5) ( +- 0.5, 1, 1.5, 2 sd, 9 samples); EX_THREE_THIRD_STEP_STDDEVS (6) ( +- 0.33, 0.67, 1 sd; 7 samples); EX_SIX_QUARTER_STEP_STDDEVS (7) ( += 0.25, 0.5, 0.75, 1.25, 1.5 sd; 13 samples)",  "1"  )
		+ XMLSchemaAttribute::attribute_w_default(  "ex2_sample_level", "exchi_sample_level", "Sampling level for chi2 sampling; see ex1_sample_level for available options.",  "0"  )
		+ XMLSchemaAttribute::attribute_w_default(  "ex3_sample_level", "exchi_sample_level", "Sampling level for chi3 sampling; see ex1_sample_level for available options",  "1"  )
		+ XMLSchemaAttribute::attribute_w_default(  "ex4_sample_level", "exchi_sample_level", "Sampling level for chi4 sampling; see ex1_sample_level for available options",  "1"  )
		+ XMLSchemaAttribute::attribute_w_default(  "ex1aro_sample_level", "exchi_sample_level", "Sampling level for aromatic chi1 sampling; see ex1_sample_level for available options",  "1"  )
		+ XMLSchemaAttribute::attribute_w_default(  "ex2aro_sample_level", "exchi_sample_level", "Sampling level for aromatic chi2 sampling; see ex1_sample_level for available options",  "0"  )
		+ XMLSchemaAttribute::attribute_w_default(  "ex1aro_exposed_sample_level", "exchi_sample_level", "Sampling level for aromatic chi1 sampling; see ex1_sample_level for available options",  "1"  )
		+ XMLSchemaAttribute::attribute_w_default(  "ex2aro_exposed_sample_level", "exchi_sample_level",  "Sampling level for aromatic chi1 sampling; see ex1_sample_level for available options",  "1"  )

		+ XMLSchemaAttribute::attribute_w_default(  "extrachi_cutoff", xsct_non_negative_integer, "Number of neighbors a residue must have before extra rotamers are used. Default 18",  utility::to_string( EXTRACHI_CUTOFF_LIMIT  ));

	return attributes;
}

void set_rotamer_sampling_data_for_RLT(
	ExtraRotamerSamplingData const & sampling_data,
	ResidueLevelTask & res_task
)
{
	res_task.or_ex1(sampling_data.ex1_);
	res_task.or_ex2(sampling_data.ex2_);
	res_task.or_ex3(sampling_data.ex3_);
	res_task.or_ex4(sampling_data.ex4_);
	res_task.or_ex1aro(sampling_data.ex1aro_);
	res_task.or_ex2aro(sampling_data.ex2aro_);
	res_task.or_ex1aro_exposed(sampling_data.ex1aro_exposed_);
	res_task.or_ex2aro_exposed(sampling_data.ex2aro_exposed_);
	res_task.or_ex1_sample_level(sampling_data.ex1_sample_level_);
	res_task.or_ex2_sample_level(sampling_data.ex2_sample_level_);
	res_task.or_ex3_sample_level(sampling_data.ex3_sample_level_);
	res_task.or_ex4_sample_level(sampling_data.ex4_sample_level_);
	res_task.or_ex1aro_sample_level(sampling_data.ex1aro_sample_level_);
	res_task.or_ex2aro_sample_level(sampling_data.ex2aro_sample_level_);
	res_task.or_ex1aro_exposed_sample_level(sampling_data.ex1aro_exposed_sample_level_);
	res_task.or_ex2aro_exposed_sample_level(sampling_data.ex2aro_exposed_sample_level_);
	res_task.or_exdna_sample_level(sampling_data.exdna_sample_level_);
	res_task.and_extrachi_cutoff(sampling_data.extrachi_cutoff_);
}


/// BEGIN ReadResfile

ReadResfile::ReadResfile() :
	parent(),
	resfile_filename_(""),
	file_was_read_(false),
	resfile_cache_(""),
	residue_selector_()
	//TODO -- UPDATE THIS WHEN NEW PRIVATE MEMBER VARIABLES ARE ADDED
{
	cache_resfile();
}

ReadResfile::ReadResfile( utility::options::OptionCollection const & options ) :
	parent(),
	resfile_filename_(""),
	file_was_read_(false),
	resfile_cache_(""),
	residue_selector_()
{
	if ( options[ basic::options::OptionKeys::packing::resfile ].user() ) {
		filename( options[ basic::options::OptionKeys::packing::resfile ]()[1] );
	} else {
		cache_resfile();
	}
}

ReadResfile::ReadResfile( std::string const & filename ) :
	parent(),
	resfile_filename_( filename ),
	file_was_read_(false),
	resfile_cache_(""),
	residue_selector_()
	//TODO -- UPDATE THIS WHEN NEW PRIVATE MEMBER VARIABLES ARE ADDED
{
	cache_resfile(); //Read in the file.
}

/// @brief Copy constructor.
/// @details Needed if a ResidueSelector is used.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
ReadResfile::ReadResfile( ReadResfile const &src ) :
	parent(),
	resfile_filename_( src.resfile_filename_ ),
	file_was_read_( src.file_was_read_ ),
	resfile_cache_( src.resfile_cache_ ),
	residue_selector_() //Copied (cloned) below
	//TODO -- UPDATE THIS WHEN NEW PRIVATE MEMBER VARIABLES ARE ADDED
{
	if ( src.residue_selector_ ) {
		residue_selector_ = src.residue_selector_->clone(); //Clone the residue selector(s) of the source.
	}
}


ReadResfile::~ReadResfile() = default;


TaskOperationOP ReadResfile::clone() const
{
	return TaskOperationOP( new ReadResfile( *this ) );
}

void
ReadResfile::apply( pose::Pose const & pose, PackerTask & task ) const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( !file_was_read_ ) {
		if ( option[ packing::resfile ].user() ) {
			parse_resfile(pose, task, option[ packing::resfile ][ 1 ] ); //Note that in this case, the selector is not applied
		} // else, if no file was provided directly or through the options system -- do not change the input PackerTask at all. Noop.
		return; //Do nothing if no resfile was read in and nothing is provided via the options system.
	}

	//At this point, we're dealing with a cached resfile.
	if ( residue_selector_ ) { //If there's a residue selector provided, use it to make a mask, then apply this TaskOperation only to unmasked (selected) residues:
		core::select::residue_selector::ResidueSubset const mask( residue_selector_->apply(pose) );
		parse_resfile_string(pose, task, resfile_filename_, resfile_cache_, mask );
	} else { //Otherwise, apply this TaskOperation to all residues:
		parse_resfile_string(pose, task, resfile_filename_, resfile_cache_ );
	}

	return;
}

/// @brief Set the residue selector.
/// @details The input selector is cloned and the clone is stored.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
ReadResfile::set_residue_selector(
	core::select::residue_selector::ResidueSelectorCOP selector_in
) {
	runtime_assert_string_msg( selector_in, "Error in core::pack::task::operation::ReadResfile::set_residue_selector(): A null pointer was passed to this function.  This shouldn't happen.  Consult a developer." );
	residue_selector_ = selector_in->clone();
}

/// @brief Get the residue selector, if one exists.  (Const-access owning pointer).
/// @details Returns NULL pointer if one does not.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
core::select::residue_selector::ResidueSelectorCOP
ReadResfile::residue_selector() const { return residue_selector_; }

/// @brief Assign the filename from options system, if the command-line flag has been set,
/// or to the empty string otherwise.
void
ReadResfile::default_filename()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( option[ packing::resfile ].user() ) {
		resfile_filename_ = option[ packing::resfile ][ 1 ];
	} else {
		resfile_filename_ = "";
	}
}

void
ReadResfile::filename( std::string const & filename )
{
	resfile_filename_ = filename;
	cache_resfile();
}

std::string const & ReadResfile::filename() const
{
	return resfile_filename_;
}

void
ReadResfile::parse_tag( TagCOP tag , DataMap &datamap )
{
	if ( tag->hasOption("filename") ) resfile_filename_ = tag->getOption<std::string>("filename");
	// special case: if "COMMANDLINE" string specified, use commandline option setting.
	// This is wholy unneccessary of course.  In the absence of a specified filename, the command line
	// will be read from, anyways.
	// if no filename is given, then the ReadResfile command will read from
	// the packing::resfile option on the command line.
	if ( resfile_filename_ == "COMMANDLINE" ) default_filename();

	if ( tag->hasOption( "selector" ) ) {
		std::string const selector_name ( tag->getOption< std::string >( "selector" ) );
		try {
			set_residue_selector( datamap.get_ptr< core::select::residue_selector::ResidueSelector const >( "ResidueSelector", selector_name ) );
		} catch ( utility::excn::Exception & e ) {
			std::string error_message = "Failed to find ResidueSelector named '" + selector_name + "' from the Datamap from ReadResfile::parse_tag()\n" + e.msg();
			throw CREATE_EXCEPTION(utility::excn::Exception,  error_message );
		}
		debug_assert(residue_selector_);
	}

	// Read in the resfile and store it:
	cache_resfile();
}

/// @brief Read in the resfile and store it, so that it
/// doesn't have to be read over and over again at apply time.
void
ReadResfile::cache_resfile() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( resfile_filename_ == "" ) {
		if ( option[ packing::resfile ].user() ) {
			resfile_filename_ = option[ packing::resfile ][ 1 ];
		}
	}

	if ( resfile_filename_ == "" ) {
		file_was_read_ = false;
		resfile_cache_ = "";
		return; //Do nothing if we have no filename to read.
	}

	utility::io::izstream file( resfile_filename_ );
	if ( !file ) {
		TR.Error << "File:" << resfile_filename_ << " not found!\n";
		TR.Error.flush();
		utility_exit_with_message( "Cannot open file " + resfile_filename_ );
	}
	resfile_cache_=""; //Clear current cache, if it exists.
	utility::slurp( file, resfile_cache_ ); //Cache the resfile in memory.
	file.close();

	file_was_read_ = true;
	return;
}

/// @brief Allows code to provide resfile contents, so that this TaskOperation doesn't directly have to
/// handle file i/o.  Handly on large systems (e.g. Blue Gene), where one might only want the
/// master process to read a file.
void
ReadResfile::set_cached_resfile(
	std::string const &file_contents
) {
	resfile_cache_ = file_contents;
	file_was_read_ = true;
}


std::string ReadResfile::keyname() { return "ReadResfile"; }

utility::tag::AttributeList
ReadResfile::xml_schema_attributes() {
	utility::tag::AttributeList attributes;
	attributes
		+ utility::tag::XMLSchemaAttribute( "filename", xs_string , "If a filename is given, read from that file. Otherwise, read the file specified on the commandline with -packing:resfile." )
		+ utility::tag::XMLSchemaAttribute( "selector", xs_string , "Optionally, a previously-defined ResidueSelector may be specified using the selector=(some string) option. If this is used, then the ResidueSelector is used as a mask, and the ReadResfile TaskOperation is applied only to those residues selected by the ResidueSelector, even if the resfile lists other residues as well." );
	return attributes;
}

void ReadResfile::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	utility::tag::AttributeList attributes = xml_schema_attributes();
	task_op_schema_w_attributes( xsd, keyname(), attributes, "Read a resfile." );
}

void ReadResfile::list_options_read( utility::options::OptionKeyList & options )
{
	options + basic::options::OptionKeys::packing::resfile;
}

TaskOperationOP ReadResfileCreator::create_task_operation() const
{
	return TaskOperationOP( new ReadResfile );
}

std::string ReadResfileCreator::keyname() const { return ReadResfile::keyname(); }

void ReadResfileCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ReadResfile::provide_xml_schema( xsd );
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

ReadResfileAndObeyLengthEvents::~ReadResfileAndObeyLengthEvents()= default;

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
	if ( !pose.observer_cache().has( pose::datacache::CacheableObserverType::LENGTH_EVENT_COLLECTOR ) ) {
		parent::apply( pose, ptask );
		return;
	}
	pose::datacache::CacheableObserverCOP len_obs = pose.observer_cache().get_const_ptr( pose::datacache::CacheableObserverType::LENGTH_EVENT_COLLECTOR );
	pose::datacache::LengthEventCollectorCOP lencollect( utility::pointer::static_pointer_cast< pose::datacache::LengthEventCollector const >( len_obs ) );

	utility::vector1< core::conformation::signals::LengthEvent > const & events( lencollect->events() );
	utility::vector1< core::id::SequenceMapping > smaps;
	TR << "ReadResfileAndObeyLengthEvents: processing " << events.size() << " length events." << std::endl;
	for ( Size i =1; i <= events.size(); ++i ) {
		TR.Debug <<  events[i] << std::endl;
		smaps.push_back( core::id::SequenceMapping( events[i] ) );
	}
	core::id::SequenceMappingOP fullsmap( core::id::combine_sequence_mappings( smaps ) );
	fullsmap->reverse();

	//2. get the resfile contents, we should probably
	//refactor ResfileReader a little bit to replace the
	//following block by one call
	std::string resfile_string;
	utility::io::izstream file( filename() );
	if ( !file ) utility_exit_with_message( "Cannot open file " + filename() );
	utility::slurp( file, resfile_string );
	std::istringstream resfile(resfile_string);
	ResfileContents contents( pose, filename(), resfile );


	//3. apply the ResfileCommands to the remapped seqpos
	for ( Size ii = 1; ii <= ptask.total_residue(); ++ii ) {
		Size ii_resfile_seqpos( (*fullsmap)[ii] );

		if ( !ii_resfile_seqpos && !apply_default_commands_to_inserts_ ) continue;

		std::list< ResfileCommandCOP > const & ii_command_list( this->resfile_commands( ii_resfile_seqpos, contents, ptask) );

		for ( auto const & command : ii_command_list ) {
			command->residue_action( ptask, ii );
		}
	} //loop over all task residues

} //ReadResfileAndObeyLengthEvents::apply(

void
ReadResfileAndObeyLengthEvents::parse_tag( TagCOP tag , DataMap & datamap )
{
	parent::parse_tag( tag, datamap );
	if ( tag->hasOption("default_commands_for_inserts") ) {
		apply_default_commands_to_inserts_ = tag->getOption<bool>("default_commands_for_inserts",true);
	}
}

/// @details note: this function will return default commands for
/// resfile_seqpos == 0
std::list< ResfileCommandCOP > const &
ReadResfileAndObeyLengthEvents::resfile_commands(
	core::Size const resfile_seqpos,
	ResfileContents const & contents,
	PackerTask const & ptask ) const
{
	if ( (resfile_seqpos == 0) || ( resfile_seqpos > ptask.total_residue()) ) {
		return contents.default_commands();
	}
	return (contents.specialized_commands_exist_for_residue( resfile_seqpos ) ?
		contents.commands_for_residue( resfile_seqpos ) : contents.default_commands() );
}

std::string ReadResfileAndObeyLengthEvents::keyname() { return "ReadResfileAndObeyLengthEvents"; }

void ReadResfileAndObeyLengthEvents::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	utility::tag::AttributeList attributes = parent::xml_schema_attributes();
	attributes + utility::tag::XMLSchemaAttribute::attribute_w_default(  "default_commands_for_inserts", xsct_rosetta_bool , "apply the ResfileCommands to the remapped seqpos.",  "1"  );
	task_op_schema_w_attributes( xsd, keyname(), attributes, "Read resfile and obey length events. IMPORTANT: only use this if any length changes are not reflected in the pose's pdb info, such as seems to be the case after vlb not quite certain on the ideal approach yet. it's prolly best to parse the resfile, and then apply the ResfileCommands to the remapped residues. this necessitates getting the ResfileContents. the code under 2. here is some duplication of ResfileReader::parse_resfile/parse_resfile_string, ideally this and ResfileReader should be refactored a bit." );
}

TaskOperationOP ReadResfileAndObeyLengthEventsCreator::create_task_operation() const
{
	return TaskOperationOP( new ReadResfileAndObeyLengthEvents );
}

std::string ReadResfileAndObeyLengthEventsCreator::keyname() const {
	return ReadResfileAndObeyLengthEvents::keyname();
}

void ReadResfileAndObeyLengthEventsCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ReadResfileAndObeyLengthEvents::provide_xml_schema( xsd );
}

/// BEGIN SetRotamerCouplings

SetRotamerCouplings::SetRotamerCouplings() = default;

SetRotamerCouplings::~SetRotamerCouplings() = default;

SetRotamerCouplings::SetRotamerCouplings( SetRotamerCouplings const & src ) :
	parent(),
	rotamer_couplings_( src.rotamer_couplings_ )
{}

SetRotamerCouplings &
SetRotamerCouplings::operator = ( SetRotamerCouplings const & rhs )
{
	rotamer_couplings_ = rhs.rotamer_couplings_;
	return *this;
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

std::string SetRotamerCouplings::keyname() { return "SetRotamerCouplings"; }

void SetRotamerCouplings::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	task_op_schema_empty( xsd, keyname(), "Rotamer couplings cannot currently be provided through the RosettaScripts interface. Under development and untested. Use at your own risk." );
}

TaskOperationOP SetRotamerCouplingsCreator::create_task_operation() const
{
	return TaskOperationOP( new SetRotamerCouplings );
}

std::string SetRotamerCouplingsCreator::keyname() const {
	return SetRotamerCouplings::keyname();
}


void SetRotamerCouplingsCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SetRotamerCouplings::provide_xml_schema( xsd );
}

/// BEGIN SetRotamerLinks

SetRotamerLinks::SetRotamerLinks() = default;

SetRotamerLinks::~SetRotamerLinks() = default;

SetRotamerLinks::SetRotamerLinks( SetRotamerLinks const & src )
:
	parent(),
	rotamer_links_( src.rotamer_links_ )
{}

SetRotamerLinks &
SetRotamerLinks::operator = ( SetRotamerLinks const & rhs )
{
	rotamer_links_ = rhs.rotamer_links_;
	return *this;
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

std::string SetRotamerLinks::keyname() { return "SetRotamerLinks"; }

void SetRotamerLinks::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	task_op_schema_empty( xsd, keyname(), "Rotamer links cannot currently be specified through the RosettaScripts interface (C++/PyRosetta only). Sets rotamer links" );
}

TaskOperationOP SetRotamerLinksCreator::create_task_operation() const
{
	return TaskOperationOP( new SetRotamerLinks );
}

std::string SetRotamerLinksCreator::keyname() const {
	return SetRotamerLinks::keyname();
}

void SetRotamerLinksCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SetRotamerLinks::provide_xml_schema( xsd );
}

/// BEGIN AppendRotamer

AppendRotamer::AppendRotamer()
: rotamer_operation_(/* 0 */)
{}

AppendRotamer::~AppendRotamer() = default;

AppendRotamer::AppendRotamer( rotamer_set::RotamerOperationOP rotamer_operation )
: rotamer_operation_(std::move( rotamer_operation ))
{}

AppendRotamer::AppendRotamer( AppendRotamer const & src )
: parent(), rotamer_operation_( src.rotamer_operation_ )
{}

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

std::string AppendRotamer::keyname() { return "AppendRotamer"; }

void AppendRotamer::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	task_op_schema_empty( xsd, keyname(), "Appends a rotamer opteration to the packer task (cannot currently be specified through RosettaScripts). Under development and untested. Use at your own risk." );
}

TaskOperationOP AppendRotamerCreator::create_task_operation() const
{
	return TaskOperationOP( new AppendRotamer );
}

std::string AppendRotamerCreator::keyname() const {
	return AppendRotamer::keyname();
}

void AppendRotamerCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AppendRotamer::provide_xml_schema( xsd );
}

/// BEGIN AppendRotamerSet

AppendRotamerSet::AppendRotamerSet()
: rotamer_set_operation_(/* 0 */)
{}

AppendRotamerSet::~AppendRotamerSet() = default;

AppendRotamerSet::AppendRotamerSet( rotamer_set::RotamerSetOperationOP rotamer_set_operation )
: rotamer_set_operation_(std::move( rotamer_set_operation ))
{}

AppendRotamerSet::AppendRotamerSet( AppendRotamerSet const & src )
: parent(), rotamer_set_operation_( src.rotamer_set_operation_ )
{}

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

std::string AppendRotamerSet::keyname() { return "AppendRotamerSet"; }

void AppendRotamerSet::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	task_op_schema_empty( xsd, keyname(), "Appends a rotamer set operation to the packer task (cannot currently be specified through RosettaScripts).  Under development and untested. Use at your own risk." );
}

TaskOperationOP AppendRotamerSetCreator::create_task_operation() const
{
	return TaskOperationOP( new AppendRotamerSet );
}

std::string AppendRotamerSetCreator::keyname() const {
	return AppendRotamerSet::keyname();
}

void AppendRotamerSetCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AppendRotamerSet::provide_xml_schema( xsd );
}


/// BEGIN AppendResidueRotamerSet

AppendResidueRotamerSet::AppendResidueRotamerSet()
:
	parent(),
	resnum_(0),
	rotamer_set_operation_(/* 0 */)
{}

AppendResidueRotamerSet::~AppendResidueRotamerSet() = default;

AppendResidueRotamerSet::AppendResidueRotamerSet( core::Size resnum,
	rotamer_set::RotamerSetOperationOP rotamer_set_operation )
: resnum_(resnum),
	rotamer_set_operation_(std::move( rotamer_set_operation ))
{}

AppendResidueRotamerSet::AppendResidueRotamerSet( AppendResidueRotamerSet const & /*src*/ ) = default;

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

std::string AppendResidueRotamerSet::keyname() { return "AppendResidueRotamerSet"; }

void AppendResidueRotamerSet::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	task_op_schema_empty( xsd, keyname(), "Appends a rotamer set operation to the packer task for the specified residue. Cannot currently be specified through RosettaScripts." );
}

TaskOperationOP AppendResidueRotamerSetCreator::create_task_operation() const
{
	return TaskOperationOP( new AppendResidueRotamerSet );
}

std::string AppendResidueRotamerSetCreator::keyname() const {
	return AppendResidueRotamerSet::keyname();
}

void AppendResidueRotamerSetCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AppendResidueRotamerSet::provide_xml_schema( xsd );
}


/// BEGIN PreserveCBeta

PreserveCBeta::~PreserveCBeta() = default;

TaskOperationOP PreserveCBeta::clone() const
{
	return TaskOperationOP( new PreserveCBeta( *this ) );
}

void
PreserveCBeta::apply( pose::Pose const &, PackerTask & task ) const
{
	task.or_preserve_c_beta( true );
}

std::string PreserveCBeta::keyname() { return "PreserveCBeta"; }

void PreserveCBeta::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	task_op_schema_empty( xsd, keyname(), "Under development and untested. Use at your own risk." );
}

TaskOperationOP PreserveCBetaCreator::create_task_operation() const
{
	return TaskOperationOP( new PreserveCBeta );
}

std::string PreserveCBetaCreator::keyname() const {
	return PreserveCBeta::keyname();
}

void PreserveCBetaCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PreserveCBeta::provide_xml_schema( xsd );
}

/// BEGIN PreventRepacking
PreventRepacking::~PreventRepacking() = default;

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

	for ( utility::vector1< core::Size >::const_iterator it(residues_to_prevent.begin()), end(residues_to_prevent.end());
			it != end; ++it ) {
		task.nonconst_residue_task(*it).prevent_repacking();
	}
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

std::string PreventRepacking::keyname() { return "PreventRepacking"; }

void PreventRepacking::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	utility::tag::AttributeList attributes;
	core::pose::attributes_for_get_resnum_selector( attributes, xsd, "resnum" );
	task_op_schema_w_attributes( xsd, keyname(), attributes , "Do not allow repacking at all for the specified residues. Freezes residues." );
}

TaskOperationOP PreventRepackingCreator::create_task_operation() const
{
	return TaskOperationOP( new PreventRepacking );
}

std::string
PreventRepackingCreator::keyname() const { return PreventRepacking::keyname(); }

void PreventRepackingCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PreventRepacking::provide_xml_schema( xsd );
}

// BEGIN RestrictYSDesign
RestrictYSDesign::~RestrictYSDesign() = default;
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
	for ( unsigned long YSresid : YSresids_ ) {
		if ( gly_switch_ ) restrict_to_aa[chemical::aa_from_name( "GLY" )] = true;
		restrict_to_aa[chemical::aa_from_name( "TYR" )] = true;
		restrict_to_aa[chemical::aa_from_name( "SER" )] = true;
		task.nonconst_residue_task(YSresid).restrict_absent_canonical_aas( restrict_to_aa );
	}
}

TaskOperationOP RestrictYSDesign::clone() const {
	return TaskOperationOP( new RestrictYSDesign( *this ) );
}

void
RestrictYSDesign::include_resid( core::Size const resid ) { YSresids_.push_back(resid); }

void
RestrictYSDesign::include_gly( bool const gly ) { gly_switch_ = gly; }


std::string RestrictYSDesign::keyname() { return "RestrictYSDesign"; }

void RestrictYSDesign::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	task_op_schema_empty( xsd, keyname() , "Restrict amino acid choices during design to Tyr and Ser. This is similar to the restricted YS alphabet used by Sidhu's group during in vitro evolution experiments. Under development and untested. Use at your own risk." );
}

TaskOperationOP RestrictYSDesignCreator::create_task_operation() const
{
	return TaskOperationOP( new RestrictYSDesign );
}

std::string RestrictYSDesignCreator::keyname() const { return RestrictYSDesign::keyname(); }

void RestrictYSDesignCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RestrictYSDesign::provide_xml_schema( xsd );
}

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

ExtraRotamers::~ExtraRotamers() = default;

TaskOperationOP ExtraRotamers::clone() const
{
	return TaskOperationOP( new ExtraRotamers( *this ) );
}

void
ExtraRotamers::apply( core::pose::Pose const & p, PackerTask & task ) const
{
	if ( resid_ != 0 ) {
		ResidueLevelTask & restask( task.nonconst_residue_task( resid_ ) );
		if ( chi_ == 1 ) restask.or_ex1( ExtraRotSample( level_ ) );
		if ( chi_ == 2 ) restask.or_ex2( ExtraRotSample( level_ ) );
		if ( chi_ == 3 ) restask.or_ex3( ExtraRotSample( level_ ) );
		if ( chi_ == 4 ) restask.or_ex4( ExtraRotSample( level_ ) );
	} else {
		// apply to all residues
		TR << "Enabling extra rotamers for chi " << chi_ << " at all positions" << std::endl;
		for ( Size ii = 1; ii <= p.size(); ++ii ) {
			ResidueLevelTask & restask( task.nonconst_residue_task( ii ) );
			if ( chi_ == 1 ) restask.or_ex1( ExtraRotSample( level_ ) );
			if ( chi_ == 2 ) restask.or_ex2( ExtraRotSample( level_ ) );
			if ( chi_ == 3 ) restask.or_ex3( ExtraRotSample( level_ ) );
			if ( chi_ == 4 ) restask.or_ex4( ExtraRotSample( level_ ) );
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

std::string ExtraRotamers::keyname() { return "ExtraRotamers"; }

void ExtraRotamers::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	define_extra_rotamers_sampling_level_restriction( xsd );

	XMLSchemaRestriction one_to_four = utility::tag::integer_range_restriction( "one_to_four", 1, 4 );
	xsd.add_top_level_element( one_to_four );

	utility::tag::AttributeList attributes;
	attributes
		+ utility::tag::XMLSchemaAttribute::attribute_w_default(  "resid", xsct_non_negative_integer, "Apply to specific residues, or =0, apply to all residues.",  "0"  )
		+ utility::tag::XMLSchemaAttribute::required_attribute( "chi", "one_to_four" , "Chi value between 1 and 4." )
		+ utility::tag::XMLSchemaAttribute::attribute_w_default(  "level", "exchi_sample_level", "The level of extra rotamers you want.",  "0"  );
	task_op_schema_w_attributes( xsd, keyname(), attributes, "The ability to add various levels of extra rotamers to specific or all residues." );
}

TaskOperationOP ExtraRotamersCreator::create_task_operation() const
{
	return TaskOperationOP( new ExtraRotamers );
}

std::string ExtraRotamersCreator::keyname() const {
	return ExtraRotamers::keyname();
}

void ExtraRotamersCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ExtraRotamers::provide_xml_schema( xsd );
}


//////////////////////////////////////////////////////
// This class could easily be expanded ...
// Someone has probably already written this, and should replace this
// dinky thing.
ExtraChiCutoff::ExtraChiCutoff() = default;

ExtraChiCutoff::ExtraChiCutoff( core::Size const resid, core::Size const extrachi_cutoff):
	resid_( resid ),
	extrachi_cutoff_( extrachi_cutoff )
{}

ExtraChiCutoff::~ExtraChiCutoff() = default;

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
		for ( Size ii = 1; ii <= p.size(); ++ii ) {
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

std::string ExtraChiCutoff::keyname() { return "ExtraChiCutoff"; }

void ExtraChiCutoff::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	utility::tag::AttributeList attributes;
	attributes
		+ utility::tag::XMLSchemaAttribute::attribute_w_default(  "resid", xsct_non_negative_integer, "Residue ID of specific residues or =0, all residues.",  "0"  )
		+ utility::tag::XMLSchemaAttribute::required_attribute( "extrachi_cutoff", xsct_non_negative_integer , "Extra rotamers chi cutoff." );
	task_op_schema_w_attributes( xsd, keyname(), attributes, "Assign a chi cutoff to specific or all residues." );
}

TaskOperationOP ExtraChiCutoffCreator::create_task_operation() const
{
	return TaskOperationOP( new ExtraChiCutoff );
}

/*std::string ExtraChiCutoffCreator::keyname() const {
return ExtraChiCutoff::keyname();
}*/

void ExtraChiCutoffCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ExtraChiCutoff::provide_xml_schema( xsd );
}

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core
