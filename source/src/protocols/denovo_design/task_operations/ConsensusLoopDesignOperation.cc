// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/task_operations/ConsensusLoopDesignOperation.cc
/// @brief Restrict loop residue identities to those commonly found in nature
/// @author Tom Linsky (tlinsky@uw.edu)

// unit headers
#include <protocols/denovo_design/task_operations/ConsensusLoopDesignOperation.hh>
#include <protocols/denovo_design/task_operations/ConsensusLoopDesignOperationCreator.hh>

// package headers
#include <protocols/denovo_design/util.hh>
#include <protocols/jd2/parser/BluePrint.hh>

// core headers
#include <core/chemical/ResidueType.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/ResidueLevelTask.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/select/residue_selector/SecondaryStructureSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/sequence/ABEGOManager.hh>

// utility headers
#include <basic/database/open.hh>
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

// boost headers
#include <boost/assign.hpp>

// c++ headers

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.task_operations.ConsensusLoopDesignOperation" );

namespace protocols {
namespace denovo_design {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace utility::tag;

TaskOperationOP
ConsensusLoopDesignOperationCreator::create_task_operation() const
{
	return TaskOperationOP( new ConsensusLoopDesignOperation );
}

void
ConsensusLoopDesignOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ConsensusLoopDesignOperation::provide_xml_schema( xsd );
}

std::string
ConsensusLoopDesignOperationCreator::keyname() const
{
	return ConsensusLoopDesignOperation::class_name();
}

// default constructor
ConsensusLoopDesignOperation::ConsensusLoopDesignOperation() :
	core::pack::task::operation::TaskOperation(),
	secstruct_( "" ),
	include_adjacent_residues_( false ),
	use_dssp_( true ),
	enrichment_threshold_( -2.0 ),
	selector_()
{
}

// destructor
ConsensusLoopDesignOperation::~ConsensusLoopDesignOperation()
{}

/// @brief make clone
core::pack::task::operation::TaskOperationOP
ConsensusLoopDesignOperation::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new ConsensusLoopDesignOperation( *this ) );
}

std::string
ConsensusLoopDesignOperation::class_name()
{
	return "ConsensusLoopDesign";
}

std::string
ConsensusLoopDesignOperation::get_name() const
{
	return ConsensusLoopDesignOperation::class_name();
}

/// @brief apply
void
ConsensusLoopDesignOperation::apply(
	core::pose::Pose const & pose,
	core::pack::task::PackerTask & task ) const
{
	LoopInfoVec info = get_loop_info( pose );
	for ( LoopInfoVec::const_iterator l=info.begin(); l!=info.end(); ++l ) {
		TR << "Restricting AAs in loop: " << *l << std::endl;
		disallow_aas( task, *l );
	}
	task.show( TR );
	TR.flush();
}

void
ConsensusLoopDesignOperation::parse_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data )
{
	set_secstruct( tag->getOption< std::string >( "secstruct", secstruct_ ) );

	if ( tag->hasOption( "blueprint" ) ) {
		set_secstruct_from_blueprint( tag->getOption< std::string >( "blueprint" ) );
		if ( secstruct_.empty() ) {
			std::stringstream msg;
			msg << "ConsensusLoopDesign: Error getting secondary structure from blueprint file" << std::endl;
			throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
		}
	}

	using core::select::residue_selector::parse_residue_selector;
	using core::select::residue_selector::ResidueSelectorCOP;
	ResidueSelectorCOP selector = parse_residue_selector( tag, data );
	if ( selector ) set_residue_selector( *selector );

	set_include_adjacent_residues( tag->getOption< bool >( "include_adjacent_residues", include_adjacent_residues_ ) );
	enrichment_threshold_ = tag->getOption< core::Real >( "threshold", enrichment_threshold_ );
	use_dssp_ = tag->getOption< bool >( "use_dssp", use_dssp_ );
}

/// @brief Returns a residue selector to be used to choose loops
/// @param[in] secstruct  Secondary structure to be used, if default selector is used
/// @details If selector_ is provided, simply return that.  If not,
///          create a default secondary structure selector using the given secstruct
core::select::residue_selector::ResidueSelectorCOP
ConsensusLoopDesignOperation::residue_selector( std::string const & secstruct ) const
{
	if ( selector_ ) return selector_;

	using core::select::residue_selector::SecondaryStructureSelector;
	using core::select::residue_selector::SecondaryStructureSelectorOP;

	SecondaryStructureSelectorOP sel( new SecondaryStructureSelector( "L" ) );
	sel->set_include_terminal_loops( false );
	sel->set_pose_secstruct( secstruct );
	return sel;
}

void
ConsensusLoopDesignOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;

	attributes
		+ XMLSchemaAttribute( "blueprint", xs_string )
		+ XMLSchemaAttribute( "residue_selector", xs_string )
		+ XMLSchemaAttribute( "include_adjacent_residues", xs_boolean )
		+ XMLSchemaAttribute( "enrichment_threshold", xs_string );

	task_op_schema_w_attributes( xsd, class_name(), attributes );
}

void
ConsensusLoopDesignOperation::set_secstruct_from_blueprint( std::string const & bpfile )
{
	using namespace core::select::residue_selector;

	protocols::jd2::parser::BluePrint bp( bpfile );
	set_secstruct( bp.secstruct() );

}

void
ConsensusLoopDesignOperation::set_secstruct( std::string const & ss )
{
	secstruct_ = ss;
}

void
ConsensusLoopDesignOperation::set_include_adjacent_residues( bool const include_res )
{
	include_adjacent_residues_ = include_res;
}

void
ConsensusLoopDesignOperation::set_enrichment_threshold( core::Real const threshold )
{
	enrichment_threshold_ = threshold;
}

AAFrequencies const &
ConsensusLoopDesignOperation::aa_frequencies(
	LoopInfo const & info,
	core::Size const loop_resid ) const
{
	ConsensusLoopDatabase const & db = *ConsensusLoopDatabase::get_instance();

	return db.frequencies( info.ss_around, info.abego, loop_resid );
}

ConsensusLoopDesignOperation::AAs
ConsensusLoopDesignOperation::forbidden_aas( AAFrequencies const & frequencies ) const
{
	AAs aas = "";
	for ( AAFrequencies::AAFrequencyMap::const_iterator f=frequencies.begin(); f!=frequencies.end(); ++f ) {
		if ( f->second->enrichment() >= enrichment_threshold_ ) continue;
		aas += f->first;
	}

	return aas;
}

utility::vector1< bool >
make_aa_bitmap_from_forbidden_aas( std::string const & forbidden_aas )
{
	utility::vector1< bool > bitmap( core::chemical::num_canonical_aas, true );
	for ( std::string::const_iterator aa=forbidden_aas.begin(); aa!=forbidden_aas.end(); ++aa ) {
		bitmap[ core::chemical::aa_from_oneletter_code( *aa ) ] = false;
	}
	return bitmap;
}

utility::vector1< bool >
make_aa_bitmap_from_allowed_aas( std::string const & allowed_aas )
{
	utility::vector1< bool > bitmap( core::chemical::num_canonical_aas, false );
	for ( std::string::const_iterator aa=allowed_aas.begin(); aa!=allowed_aas.end(); ++aa ) {
		bitmap[ core::chemical::aa_from_oneletter_code( *aa ) ] = true;
	}
	return bitmap;
}

ConsensusLoopDesignOperation::AAs
ConsensusLoopDesignOperation::compute_aas_after_disallowing(
	AAs const & aas,
	core::pack::task::ResidueLevelTask const & task ) const
{
	AAs remaining_aas = "";
	for ( core::pack::task::ResidueLevelTask::ResidueTypeCOPListConstIter t=task.allowed_residue_types_begin(); t!=task.allowed_residue_types_end(); ++t ) {
		if ( std::find( aas.begin(), aas.end(), (*t)->name1() ) == aas.end() ) {
			remaining_aas += (*t)->name1();
		}
	}
	return remaining_aas;
}

ConsensusLoopDesignOperation::AAs
ConsensusLoopDesignOperation::compute_best_allowed_aas(
	AAFrequencies const & aa_freqs,
	core::pack::task::ResidueLevelTask const & task ) const
{
	std::pair< char, char > best_aas = std::make_pair( 0, 0 );
	std::pair< core::Real, core::Real > best_scores = std::make_pair( -1200, -1200 );

	utility::vector1< char > no_freq;
	for ( core::pack::task::ResidueLevelTask::ResidueTypeCOPListConstIter t=task.allowed_residue_types_begin(); t!=task.allowed_residue_types_end(); ++t ) {
		if ( !aa_freqs.has_frequency( (*t)->name1() ) ) {
			no_freq.push_back( (*t)->name1() );
			continue;
		}
		core::Real const enrichment = aa_freqs.frequency( (*t)->name1() ).enrichment();
		if ( best_aas.first == 0 ) {
			best_aas.first = (*t)->name1();
			best_scores.first = enrichment;
		} else if ( enrichment > best_scores.first ) {
			best_aas.second = best_aas.first;
			best_scores.second = best_scores.first;
			best_aas.first = (*t)->name1();
			best_scores.first = enrichment;
		} else if ( best_aas.second == 0 ) {
			best_aas.second = (*t)->name1();
			best_scores.second = enrichment;
		} else if ( enrichment > best_scores.second ) {
			best_aas.second = (*t)->name1();
			best_scores.second = enrichment;
		}
	}
	std::stringstream aas;
	aas << best_aas.first << best_aas.second;
	TR.Debug << "Best aas are " << aas.str() << " best enrichments are: " << best_scores.first << " " << best_scores.second << std::endl;

	return aas.str();
}

void
ConsensusLoopDesignOperation::disallow_aas(
	core::pack::task::PackerTask & task,
	LoopInfo const & loop ) const
{
	core::Size start = 1;
	if ( !include_adjacent_residues_ ) start += 1;

	core::Size stop = loop.abego.size();
	if ( !include_adjacent_residues_ ) stop -= 1;

	core::Size pose_resid = loop.startres;
	if ( include_adjacent_residues_ ) pose_resid -= 1;

	if ( ! ConsensusLoopDatabase::get_instance()->has_frequencies( loop.ss_around, loop.abego ) ) {
		TR << "Loop " << loop << " not found in database" << std::endl;
		return;
	}

	for ( core::Size resid=start; resid<=stop; ++resid, ++pose_resid ) {
		AAFrequencies const & aa_freqs = aa_frequencies( loop, resid );
		AAs const aas = forbidden_aas( aa_freqs );
		if ( aas.empty() ) continue;
		TR << "Residue: " << pose_resid << "; forbidden aas: " << aas << std::endl;
		debug_assert( pose_resid <= task.total_residue() );

		AAs const remaining_aas = compute_aas_after_disallowing( aas, task.residue_task( pose_resid ) );

		utility::vector1< bool > aa_bitmap;
		if ( remaining_aas.empty() ) {
			AAs const next_best_aas = compute_best_allowed_aas( aa_freqs, task.residue_task( pose_resid ) );
			TR.Warning << "WARNING: ConsensusLoopDesign would disallow all amino acids at position "
				<< pose_resid << ". Allowing the best two currently-allowed amino acids (" << next_best_aas
				<< ") even though their frequency scores are below the user-specified threshold of "
				<< enrichment_threshold_ << std::endl;
			aa_bitmap = make_aa_bitmap_from_allowed_aas( next_best_aas );
		} else {
			aa_bitmap = make_aa_bitmap_from_forbidden_aas( aas );
		}
		task.nonconst_residue_task( pose_resid ).restrict_absent_canonical_aas( aa_bitmap );
	}
}

/// @brief Gets secondary structure to be used for determining what is surrounding
///        the loops
/// @param[in] pose input pose
/// @returns Secondary structure string according to the following rules:
///          1. If secstruct_ is set, return that.
///          2. If use_dssp_ is set, return string computed from DSSP.
///          3. return pose.secstruct()
std::string
ConsensusLoopDesignOperation::get_secstruct( core::pose::Pose const & pose ) const
{
	if ( !secstruct_.empty() ) return secstruct_;

	if ( use_dssp_ ) {
		core::scoring::dssp::Dssp dssp( pose );
		std::string const dssp_ss = dssp.get_dssp_secstruct();
		TR.Debug << "Secondary structure determined by DSSP: " << dssp_ss << std::endl;
		return dssp_ss;
	}

	TR.Debug << "Secondary structure determined by pose.secstruct(): "
		<< pose.secstruct() << std::endl;

	return pose.secstruct();
}

LoopInfoVec
ConsensusLoopDesignOperation::get_loop_info( core::pose::Pose const & pose ) const
{
	std::string const ss = get_secstruct( pose );
	if ( ss.size() != pose.size() ) {
		std::stringstream msg;
		msg << class_name() << "::get_loop_info(): Number of residues in secondary structure ("
			<< ss.size() << ") does not match pose size (" << pose.size()
			<< ")" << std::endl;
		utility_exit_with_message( msg.str() );
	}

	core::select::residue_selector::ResidueSubset const subset = residue_selector( ss )->apply( pose );

	return loop_info_from_subset( pose, ss, subset );
}

LoopInfoVec
ConsensusLoopDesignOperation::loop_info_from_subset(
	core::pose::Pose const & pose,
	std::string const & ss,
	core::select::residue_selector::ResidueSubset const & subset ) const
{
	utility::vector1< std::string > abego = core::sequence::get_abego( pose, 1 );
	debug_assert( subset.size() == ss.size() );
	debug_assert( subset.size() == abego.size() );

	LoopInfoVec retval;
	bool inloop = false;
	LoopInfo info;
	for ( core::Size res=1; res<=subset.size(); ++res ) {
		if ( !inloop && subset[res] && ( res > 1 ) ) {
			inloop = true;
			info.startres = res;
			info.ss_around.before = ss[ res - 2 ];
			info.abego += abego[ res - 1 ];
		}

		if ( inloop && !subset[res] ) {
			info.ss_around.after = ss[ res - 1 ];
			info.abego += abego[ res ];
			retval.push_back( info );
			info = LoopInfo();
			inloop = false;
		}

		if ( inloop ) {
			debug_assert( abego[ res ].size() == 1 );
			info.abego += abego[ res ];
		}
	}
	return retval;
}

void
ConsensusLoopDesignOperation::set_residue_selector( core::select::residue_selector::ResidueSelector const & selector_val )
{
	selector_ = selector_val.clone();
}

/////////////////// AAFrequency /////////////////////////////////////////

AAFrequency::AAFrequency( core::Real const freq, core::Real const enrich ):
	utility::pointer::ReferenceCount(),
	frequency_( freq ),
	enrichment_( enrich )
{}

core::Real
AAFrequency::enrichment() const
{
	return enrichment_;
}

core::Real
AAFrequency::frequency() const
{
	return frequency_;
}

/////////////////// AAFrequencies ///////////////////////////////////////

AAFrequencies::AAFrequencies():
	utility::pointer::ReferenceCount(),
	aa_freq_()
{}

bool
AAFrequencies::has_frequency( char const aa ) const
{
	return ( aa_freq_.find( aa ) != aa_freq_.end() );
}

void
AAFrequencies::set_frequency( char const aa, AAFrequencyCOP frequency )
{
	if ( has_frequency( aa ) ) {
		std::stringstream msg;
		msg << "Error: overwriting frequency data for " << aa;
		utility_exit_with_message( msg.str() );
	}
	aa_freq_.insert( std::make_pair( aa, frequency ) );
}

AAFrequency const &
AAFrequencies::frequency( char const aa ) const
{
	AAFrequencyMap::const_iterator aa_f = aa_freq_.find( aa );
	if ( aa_f == aa_freq_.end() ) {
		std::stringstream msg;
		msg << "AAFrequencies: could not find " << aa << " in the frequency map." << std::endl;
		utility_exit_with_message( msg.str() );
	}
	if ( ! aa_f->second ) {
		std::stringstream msg;
		msg << "AAFrequencies: data for aa " << aa << " is NULL in the frequency map." << std::endl;
		utility_exit_with_message( msg.str() );
	}
	return *aa_f->second;
}

AAFrequencies::AAFrequencyMap::const_iterator
AAFrequencies::begin() const
{
	return aa_freq_.begin();
}

AAFrequencies::AAFrequencyMap::const_iterator
AAFrequencies::end() const
{
	return aa_freq_.end();
}

/////////////////// ConsensusLoopDatabase ///////////////////////////////

ConsensusLoopDatabase::ConsensusLoopDatabase():
	ss_map_()
{
	read_db();
}

ConsensusLoopDatabase::~ConsensusLoopDatabase()
{}

bool
ConsensusLoopDatabase::has_frequencies(
	SurroundingSS const & surrounding,
	Abego const & loop_abego ) const
{
	SurroundingSSMap::const_iterator ss = ss_map_.find( surrounding );
	if ( ss == ss_map_.end() ) {
		return false;
	}

	AbegoToLoopAAsMap::const_iterator ab = ss->second.find( loop_abego );
	if ( ab == ss->second.end() ) {
		return false;
	}
	return true;
}

bool
ConsensusLoopDatabase::has_frequencies(
	SurroundingSS const & surrounding,
	Abego const & loop_abego,
	core::Size const loop_resid ) const
{
	SurroundingSSMap::const_iterator ss = ss_map_.find( surrounding );
	if ( ss == ss_map_.end() ) {
		return false;
	}

	AbegoToLoopAAsMap::const_iterator ab = ss->second.find( loop_abego );
	if ( ab == ss->second.end() ) {
		return false;
	}

	if ( ab->second.size() < loop_resid ) {
		return false;
	}

	return true;
}

bool
ConsensusLoopDatabase::has_frequency(
	SurroundingSS const & surrounding,
	Abego const & loop_abego,
	core::Size const loop_resid,
	char const aa ) const
{
	if ( !has_frequencies( surrounding, loop_abego, loop_resid ) ) {
		return false;
	}

	return frequencies( surrounding, loop_abego, loop_resid ).has_frequency( aa );
}

AAFrequencies const &
ConsensusLoopDatabase::frequencies(
	SurroundingSS const & surrounding,
	Abego const & loop_abego,
	core::Size const loop_resid ) const
{
	SurroundingSSMap::const_iterator ss = ss_map_.find( surrounding );
	if ( ss == ss_map_.end() ) {
		std::stringstream msg;
		msg << "ConsensusLoopDatabase: could not find data for loops starting/ending with"
			<< surrounding << "." << std::endl;
		utility_exit_with_message( msg.str() );
	}

	AbegoToLoopAAsMap::const_iterator ab = ss->second.find( loop_abego );
	if ( ab == ss->second.end() ) {
		std::stringstream msg;
		msg << "ConsensusLoopDatabase: could not find data for loops starting/ending with "
			<< surrounding << " with abego " << loop_abego << "." << std::endl;
		utility_exit_with_message( msg.str() );
	}

	if ( ab->second.size() < loop_resid ) {
		std::stringstream msg;
		msg << "ConsensusLoopDatabase: could not find data for loops starting/ending with "
			<< surrounding << " with abego " << loop_abego << " for residue " << loop_resid << "." << std::endl;
		utility_exit_with_message( msg.str() );
	}

	return ab->second[ loop_resid ];
}

AAFrequency const &
ConsensusLoopDatabase::frequency(
	SurroundingSS const & surrounding,
	Abego const & loop_abego,
	core::Size const loop_resid,
	char const aa ) const
{
	return frequencies( surrounding, loop_abego, loop_resid ).frequency( aa );
}

void
ConsensusLoopDatabase::set_frequency(
	SurroundingSS const & surrounding,
	std::string const & loop_abego,
	core::Size const loop_resid,
	char const aa,
	AAFrequencyCOP frequency )
{
	SurroundingSSMap::iterator ss_data = ss_map_.find( surrounding );
	if ( ss_data == ss_map_.end() ) {
		ss_data = ss_map_.insert( std::make_pair( surrounding, AbegoToLoopAAsMap() ) ).first;
	}
	debug_assert( ss_data != ss_map_.end() );

	AbegoToLoopAAsMap::iterator seq = ss_data->second.find( loop_abego );
	if ( seq == ss_data->second.end() ) {
		seq = ss_data->second.insert( std::make_pair( loop_abego, LoopAAs() ) ).first;
	}
	debug_assert( seq != ss_data->second.end() );

	LoopAAs & loop_aas = seq->second;
	while ( loop_aas.size() < loop_resid ) {
		loop_aas.push_back( AAFrequencies() );
	}
	loop_aas[ loop_resid ].set_frequency( aa, frequency );
}

void
ConsensusLoopDatabase::read_db()
{
	// constants
	static std::string const db_file = "protocol_data/denovo_design/aa_abego_frequencies.gz";
	static std::map< char, char > const ss_to_abego = boost::assign::map_list_of ('E', 'B')('H', 'A');

	utility::io::izstream infile;
	basic::database::open( infile, db_file );

	// database format is: <Before SS> <After SS> <Loop Abego> <Allowed AAs1> <Allowed AAs2> ... <Allowed AAsN>
	if ( !infile.good() ) {
		throw utility::excn::EXCN_Msg_Exception( "ConsensusLoopDesign: Could not open database file " + db_file + "\n" );
	}

	while ( infile.good() ) {
		char beforess, afterss, aa;
		core::Size loop_resid, n_examples;
		std::string loop_abego;
		core::Real frequency, bg_frequency, score;

		infile >> beforess >> afterss
			>> loop_abego
			>> loop_resid
			>> n_examples
			>> aa
			>> frequency
			>> bg_frequency
			>> score;

		std::map< char, char >::const_iterator preabego = ss_to_abego.find( beforess );
		std::map< char, char >::const_iterator postabego = ss_to_abego.find( afterss );
		debug_assert( preabego != ss_to_abego.end() );
		debug_assert( postabego != ss_to_abego.end() );

		// construct full loop abego
		std::stringstream abegostream;
		abegostream << preabego->second << loop_abego << postabego->second;

		// set in map
		set_frequency( SurroundingSS( beforess, afterss ),
			abegostream.str(),
			loop_resid,
			aa,
			AAFrequencyOP( new AAFrequency( frequency, score ) ) );

		if ( !infile.good() ) {
			break;
		}
	}
}

bool
SurroundingSS::operator<( SurroundingSS const & ab2 ) const
{
	if ( before < ab2.before ) {
		return true;
	}
	if ( after < ab2.after ) {
		return true;
	}
	return false;
}

std::ostream &
operator<<( std::ostream & os, SurroundingSS const & ss )
{
	os << ss.before << ss.after;
	return os;
}

std::ostream &
operator<<( std::ostream & os, LoopInfo const & info )
{
	os << "Start: " << info.startres << " Abego: " << info.abego << " Before: " << info.ss_around.before << " After: " << info.ss_around.after;
	return os;
}

} // task_operations
} // denovo_design
} // protocols
