// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/denovo_design/task_operations/ConsensusLoopDesign.cc
/// @brief Restrict loop residue identities to those commonly found in nature
/// @author Tom Linsky (tlinsky@uw.edu)

// unit headers
#include <protocols/denovo_design/task_operations/ConsensusLoopDesign.hh>
#include <protocols/denovo_design/task_operations/ConsensusLoopDesignCreator.hh>

// package headers
#include <protocols/denovo_design/util.hh>

// core headers
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/residue_selector/SecondaryStructureSelector.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/sequence/ABEGOManager.hh>

// utility headers
#include <basic/database/open.hh>
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>
#include <utility/tag/Tag.hh>

// boost headers
#include <boost/assign.hpp>

// c++ headers

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.task_operations.ConsensusLoopDesign" );

namespace protocols {
namespace denovo_design {
namespace task_operations {

core::pack::task::operation::TaskOperationOP
ConsensusLoopDesignOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new ConsensusLoopDesignOperation );
}

std::string
ConsensusLoopDesignOperationCreator::keyname() const
{
	return "ConsensusLoopDesign";
}

// default constructor
ConsensusLoopDesignOperation::ConsensusLoopDesignOperation() :
	core::pack::task::operation::TaskOperation(),
	selector_()
{
	using namespace core::pack::task::residue_selector;
	SecondaryStructureSelectorOP sel( new SecondaryStructureSelector( "L" ) );
	sel->set_include_terminal_loops( false );
	set_selector( sel );
	read_db();
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
ConsensusLoopDesignOperation::get_name() const
{
	return "ConsensusLoopDesign";
}

/// @brief apply
void
ConsensusLoopDesignOperation::apply(
	core::pose::Pose const & pose,
	core::pack::task::PackerTask & task ) const
{
	LoopInfoVec info = get_loop_info( pose );
	for ( LoopInfoVec::const_iterator l = info.begin(), endl = info.end(); l != endl; ++l ) {
		TR << "Restricting AAs in loop: " << *l << std::endl;
		disallow_aas( task, *l );
	}
	task.show(TR);
}

void
ConsensusLoopDesignOperation::parse_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data )
{
	if ( tag->hasOption( "residue_selector" ) ) {
		set_selector( get_residue_selector( data, tag->getOption< std::string >( "residue_selector" ) ) );
	}
}

LoopAAs const &
ConsensusLoopDesignOperation::allowed_aas(
	SurroundingSS const & surrounding,
	std::string const & loop_abego ) const
{
	static LoopAAs const emptylist;

	ConsensusSequenceTable::const_iterator ss_data = seqtable_.find( surrounding );
	if ( ss_data == seqtable_.end() ) {
		TR << "ConsensusLoopDesign: no data found for loops between " << surrounding.before << " and " << surrounding.after << std::endl;
		return emptylist;
	}
	debug_assert( ss_data != seqtable_.end() );

	AbegoToSequenceMap::const_iterator abegotoseq = ss_data->second.find( loop_abego );
	if ( abegotoseq == ss_data->second.end() ) {
		TR << "ConsensusLoopDesign: no data found for loop with abego " << loop_abego << " between secondary structure elements " << surrounding.before << " and " << surrounding.after << std::endl;
		return emptylist;
	}
	debug_assert( abegotoseq != ss_data->second.end() );
	return abegotoseq->second;
}

void
ConsensusLoopDesignOperation::set_allowed_aas(
	SurroundingSS const & surrounding,
	std::string const & loop_abego,
	LoopAAs const & allowed_aa )
{
	ConsensusSequenceTable::iterator ss_data = seqtable_.find( surrounding );
	if ( ss_data == seqtable_.end() ) {
		ss_data = seqtable_.insert( std::make_pair( surrounding, AbegoToSequenceMap() ) ).first;
	}
	debug_assert( ss_data != seqtable_.end() );

	AbegoToSequenceMap::iterator seq = ss_data->second.find( loop_abego );
	if ( seq == ss_data->second.end() ) {
		seq = ss_data->second.insert( std::make_pair( loop_abego, allowed_aa ) ).first;
	} else {
		ss_data->second[ loop_abego ] = allowed_aa;
	}
}

void
ConsensusLoopDesignOperation::read_db()
{
	static std::string const db_file = "protocol_data/denovo_design/aa_abego_db";
	utility::io::izstream infile;
	basic::database::open( infile, db_file );

	// database format is: <Before SS> <After SS> <Loop Abego> <Allowed AAs1> <Allowed AAs2> ... <Allowed AAsN>
	if ( !infile.good() ) {
		throw utility::excn::EXCN_Msg_Exception( "ConsensusLoopDesign: Could not open database file " + db_file + "\n" );
	}
	static std::map< char, char > const ss_to_abego = boost::assign::map_list_of ('E', 'B')('H', 'A');

	while ( infile.good() ) {
		char beforess, afterss;
		std::string loop_abego;
		LoopAAs allowed_aa;
		infile >> beforess >> afterss >> loop_abego;

		std::stringstream abegostream;
		std::map< char, char >::const_iterator preabego = ss_to_abego.find( beforess );
		std::map< char, char >::const_iterator postabego = ss_to_abego.find( afterss );
		debug_assert( preabego != ss_to_abego.end() );
		debug_assert( postabego != ss_to_abego.end() );
		abegostream << preabego->second << loop_abego << postabego->second;
		for ( core::Size i=1, endi=abegostream.str().size(); i<=endi; ++i ) {
			std::string aas;
			infile >> aas;
			allowed_aa.push_back( aas );
		}
		if ( !infile.good() ) {
			break;
		}
		set_allowed_aas( SurroundingSS( beforess, afterss ), abegostream.str(), allowed_aa );
	}
}

utility::vector1< bool >
make_aa_bitmap( std::string const & aas )
{
	utility::vector1< bool > bitmap( core::chemical::num_canonical_aas, false );
	for ( core::Size i=1, endi=aas.size(); i<=endi; ++i ) {
		bitmap[ core::chemical::aa_from_oneletter_code( aas[ i-1 ] ) ] = true;
	}
	return bitmap;
}

void
ConsensusLoopDesignOperation::disallow_aas(
	core::pack::task::PackerTask & task,
	LoopInfo const & loop ) const
{
	LoopAAs aas = allowed_aas( loop.ss_around, loop.abego );
	if ( aas.empty() ) {
		return;
	}
	core::Size lres = loop.startres;
	for ( LoopAAs::const_iterator a = aas.begin(), enda = aas.end(); a != enda; ++a, ++lres ) {
		debug_assert( lres <= task.total_residue() );
		TR << "Residue: " << lres << "; allowed aas: " << *a << std::endl;
		utility::vector1< bool > aa_bitmap = make_aa_bitmap( *a );
		task.nonconst_residue_task( lres ).restrict_absent_canonical_aas( aa_bitmap );
	}
}

LoopInfoVec
ConsensusLoopDesignOperation::get_loop_info( core::pose::Pose const & pose ) const
{
	debug_assert( selector_ );
	core::pack::task::residue_selector::ResidueSubset subset = selector_->apply( pose );

	core::scoring::dssp::Dssp dssp( pose );
	std::string const ss = dssp.get_dssp_secstruct();

	return loop_info_from_subset( pose, ss, subset );
}

LoopInfoVec
ConsensusLoopDesignOperation::loop_info_from_subset(
	core::pose::Pose const & pose,
	std::string const & ss,
	core::pack::task::residue_selector::ResidueSubset const & subset ) const
{
	utility::vector1< std::string > abego = core::sequence::get_abego( pose, 1 );
	debug_assert( subset.size() == ss.size() );
	debug_assert( subset.size() == abego.size() );

	LoopInfoVec retval;
	bool inloop = false;
	LoopInfo info;
	for ( core::Size res=1, endr=subset.size(); res<=endr; ++res ) {

		if ( !inloop && subset[res] && ( res > 1 ) ) {
			inloop = true;
			info.startres = res - 1;
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
ConsensusLoopDesignOperation::set_selector( core::pack::task::residue_selector::ResidueSelectorCOP selector_val )
{
	debug_assert( selector_val );
	selector_ = selector_val;
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
operator<<( std::ostream & os, LoopInfo const & info )
{
	os << "Start: " << info.startres << " Abego: " << info.abego << " Before: " << info.ss_around.before << " After: " << info.ss_around.after;
	return os;
}

//namespaces
} // task_operations
} // denovo_design
} // protocols
