// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/PDBInfo.cc
/// @brief  class to hold PDB information so it's not loose in the pose
/// @author Steven Lewis
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// Unit headers
#include <core/pose/PDBInfo.hh>

// Project headers
#include <core/chemical/types.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/signals/IdentityEvent.hh>
#include <core/conformation/signals/LengthEvent.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBPoseMap.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/py/PyAssert.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>

// C++ headers
#include <algorithm>

#include <utility/vector0.hh>
#include <ObjexxFCL/format.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/conformation/signals/ConnectionEvent.hh>
//using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

using core::chemical::chr_chains;

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace core {
namespace pose {

static THREAD_LOCAL basic::Tracer TR( "core.pose.PDBInfo" );

core::pose::UnrecognizedAtomRecord::UnrecognizedAtomRecord() {}



/// @brief default constructor, obsolete is *true*
PDBInfo::PDBInfo() :
	Super(),
	obsolete_( true ),
	name_( "" ),
	modeltag_( "" ),
	conf_( /* NULL */ ),
	num_unrecognized_res_(0),
	num_unrecognized_atoms_(0)
{}


/// @brief size constructor (ensure space for 'n' residue records),
///  obsolete is *true*
PDBInfo::PDBInfo( Size const n ) :
	Super(),
	obsolete_( true ),
	name_( "" ),
	modeltag_( "" ),
	residue_rec_( n ),
	conf_( /* NULL */ ),
	num_unrecognized_res_(0),
	num_unrecognized_atoms_(0)
{}


/// @brief Pose constructor (ensures space for residue and atom records
///  relative to Pose)
/// @param[in] pose  Pose
/// @param[in] init  if true (default), then residue records are initialized
///  and obsolete set to false, otherwise obsolete is true
///  using Pose residue numbering and chains of the Residues in the Conformation
PDBInfo::PDBInfo(
	Pose const & pose,
	bool init
) :
	Super(),
	obsolete_( !init ),
	name_( "" ),
	modeltag_( "" ),
	residue_rec_( pose.size() ),
	conf_( /* NULL */ ),
	num_unrecognized_res_(0),
	num_unrecognized_atoms_(0)
{
	resize_atom_records( pose );

	// initialize using Pose residue numbering and chains of the Residues
	// in the Pose's Conformation
	if ( init ) {
		// do this manually to save on method calls
		for ( Size i = 1, ie = pose.size(); i <= ie; ++i ) {
			ResidueRecord & rr = residue_rec_[ i ];
			// fpd: wrap around if > 62 chains.
			// rhiju: note -- got rid of _ as first char of chr_chains.
			rr.chainID = chr_chains[ (pose.residue( i ).chain() - 1)  % chr_chains.size() ];
			rr.resSeq = static_cast< int >( i );
		}
		rebuild_pdb2pose();
	}
}


/// @brief copy constructor
PDBInfo::PDBInfo( PDBInfo const & info ) :
	Super( info ),
	obsolete_( info.obsolete_ ),
	name_( info.name_ ),
	modeltag_( info.modeltag_ ),
	header_information_( info.header_information_ ),
	remarks_( info.remarks_ ),
	residue_rec_( info.residue_rec_ ),
	pdb2pose_( info.pdb2pose_ ),
	conf_( /* NULL */ ), // no observation on copy
	unrecognized_atoms_       ( info.unrecognized_atoms_ ),
	unrecognized_res_num2name_( info.unrecognized_res_num2name_ ),
	unrecognized_res_size_    ( info.unrecognized_res_size_ ),
	num_unrecognized_res_     ( info.num_unrecognized_res_ ),
	num_unrecognized_atoms_   ( info.num_unrecognized_atoms_ ),
	crystinfo_( info.crystinfo_ )
{}


/// @brief default destructor
PDBInfo::~PDBInfo()
{
	detach_from(); // stop observing Conformation
}


/// @brief copy assignment
/// @details any Conformation already being observed stays constant, there is no
///  re-assignment
PDBInfo &
PDBInfo::operator =( PDBInfo const & info )
{
	if ( this != &info ) {
		Super::operator =( info );

		obsolete_ = info.obsolete_;

		name_ = info.name_;
		modeltag_ = info.modeltag_;
		remarks_ = info.remarks_;
		if ( info.header_information_ ) {
			header_information_ = io::HeaderInformationOP( new io::HeaderInformation(*info.header_information_) );
		}

		residue_rec_ = info.residue_rec_;

		pdb2pose_ = info.pdb2pose_;

		unrecognized_atoms_ = info.unrecognized_atoms_;
		unrecognized_res_num2name_ = info.unrecognized_res_num2name_;
		unrecognized_res_size_ = info.unrecognized_res_size_;
		num_unrecognized_res_ = info.num_unrecognized_res_;
		num_unrecognized_atoms_ = info.num_unrecognized_atoms_;
		crystinfo_ = info.crystinfo_;
	}
	return *this;
}


/// @brief is this PDBInfo currently observing a conformation?
/// @return the Conformation being observed, otherwise NULL
core::conformation::ConformationCAP
PDBInfo::is_observing() {
	return conf_;
}


/// @brief attach to Conformation and begin observation
void
PDBInfo::attach_to( core::conformation::Conformation & conf ) {
	// detach from prior Conformation if necessary
	if ( !conf_.expired() ) {
		detach_from();
	}

	conf.attach_connection_obs( &PDBInfo::on_connection_change, this );
	conf.attach_identity_obs( &PDBInfo::on_identity_change, this );
	conf.attach_length_obs( &PDBInfo::on_length_change, this );

	conf_ = conf.get_self_weak_ptr();
}


/// @brief detach from Conformation and stop observation
/// @remarks takes no arguments because PDBInfo can only observe one
///  Conformation at a time
void
PDBInfo::detach_from() {
	core::conformation::ConformationCOP conf = conf_.lock();
	if ( conf ) {
		conf->detach_connection_obs( &PDBInfo::on_connection_change, this );
		conf->detach_identity_obs( &PDBInfo::on_identity_change, this );
		conf->detach_length_obs( &PDBInfo::on_length_change, this );
	}

	conf_.reset(); // to NULL
}


/// @brief update when connection to Conformation is changed
void
PDBInfo::on_connection_change( core::conformation::signals::ConnectionEvent const & event ) {
	using core::conformation::signals::ConnectionEvent;

	switch ( event.tag ) {
	case ConnectionEvent::DISCONNECT :
		obsolete( true );
		detach_from();
		break;
	case ConnectionEvent::TRANSFER :
		// Disconnect -- PDBInfo does not honor TRANSFER tag.
		break;
	default : // do nothing
		break;
	}

}


/// @brief update atom records when residue identity changes in Conformation
void
PDBInfo::on_identity_change( core::conformation::signals::IdentityEvent const & event ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using core::conformation::signals::IdentityEvent;

	switch( event.tag ) {
	case IdentityEvent::INVALIDATE : {
		obsolete( true );
		// detach for now, in the future consider whether it's appropriate to
		// rebuild all the data
		detach_from();
		break;
	}
	case IdentityEvent::RESIDUE : {
		// replace and zero the records
		if ( option[ in::preserve_crystinfo ]() ) {
			replace_res_remap_bfactors( event.position, *event.residue ); //fpd  maintain reasonable bfactors
		} else {
			replace_res( event.position, event.residue->natoms() );
		}
	}
	default : {
		// do nothing, fall through
		break;
	}
	}
}


/// @brief update residue and atom records when length changes in Conformation,
///  obsoletes PDBInfo
/// @details Residue and atoms records in PDBInfo will be automatically resized
///  and the obsolete flag will be set to true, causing PDBInfo to not be used
///  on pdb output.  Existing data inside the PDBInfo is not touched.  After
///  filling in the data for the newly updated residues, setting
///  PDBInfo::obsolete( false ) will re-enable usage of PDBInfo on pdb output.
///  If PDBInfo receives a LengthEvent::INVALIDATE, it will obsolete and then
///  detach itself for safety.
void
PDBInfo::on_length_change( core::conformation::signals::LengthEvent const & event ) {
	using core::conformation::signals::LengthEvent;

	switch( event.tag ) {
	case LengthEvent::INVALIDATE : {
		obsolete( true );
		detach_from(); // if we don't detach, an array-out-of-bounds may occur on future signals
		break;
	}
	case LengthEvent::RESIDUE_APPEND : {
		append_res( event.position, event.residue->natoms(), std::abs( event.length_change ) );
		obsolete( true );
		break;
	}
	case LengthEvent::RESIDUE_PREPEND : {
		prepend_res( event.position, event.residue->natoms(), std::abs( event.length_change ) );
		obsolete( true );
		break;
	}
	case LengthEvent::RESIDUE_DELETE : {
		delete_res( event.position, std::abs( event.length_change ) );
		//Shouldn't obsolete PDBInfo on delete - every residue that still exists has valid information
		//obsolete( true );
		break;
	}
	default : {
		// do nothing, fall through
		break;
	}
	}
}


/// @brief resize for 'n' residue records
/// @details Leaves atom record state inconsistent.  Atom records for
///  remaining residues are untouched while new residues have no atom
///  records, so make sure and call one of resize_atom_records()
///  afterwards if necessary.
/// @warning Do not use this method for ins/del of residues, as it leaves
///  the data state inconsistent.  See append_res/prepend_res/delete_res
///  for that type of functionality.
void
PDBInfo::resize_residue_records( Size const n )
{
	residue_rec_.resize( n );
	rebuild_pdb2pose();
}


/// @brief ensure 'n' available atom records for particular residue
/// @param[in] res  residue
/// @param[in] n  number of atoms
/// @param[in] zero  if true, zero the atom records for this residue
void
PDBInfo::resize_atom_records(
	Size const res,
	Size const n,
	bool const zero
)
{
	if ( zero ) {
		residue_rec_[ res ].atomRec = AtomRecords( n );
	} else {
		residue_rec_[ res ].atomRec.resize( n );
	}
}


/// @brief ensure 'n' available atom records for every residue
/// @param[in] res  residue
/// @param[in] n  number of atoms
/// @param[in] zero  if true, zero the atom records
void
PDBInfo::resize_atom_records(
	Size const n,
	bool const zero
)
{
	if ( zero ) {
		for ( auto & i : residue_rec_ ) {
			i.atomRec = AtomRecords( n );
		}
	} else {
		for ( auto & i : residue_rec_ ) {
			i.atomRec.resize( n );
		}
	}
}


/// @brief update number of atom records with respect to atoms in Pose
/// @details Number of internally available atom records will be adjusted
///  to match number of atoms within each residue in Pose.  Only newly
///  created records will be zeroed, any existing records are untouched.
void
PDBInfo::resize_atom_records( Pose const & pose )
{
	using core::conformation::Residue;

	debug_assert( residue_rec_.size() == pose.size() );

	for ( Size r = 1, re = pose.size(); r <= re; ++r ) {
		residue_rec_[ r ].atomRec.resize( pose.residue( r ).natoms() );
	}
}


/// @brief tighten memory usage
void
PDBInfo::tighten_memory()
{
	// tighten remarks vector
	if ( remarks_.capacity() > remarks_.size() ) {
		Remarks( remarks_ ).swap( remarks_ );
	}

	// tighten residue vector
	if ( residue_rec_.capacity() > residue_rec_.size() ) {
		ResidueRecords( residue_rec_ ).swap( residue_rec_ );
	}

	// tighten each of the atom vectors
	for ( auto & i : residue_rec_ ) {
		if ( i.atomRec.capacity() > i.atomRec.size() ) {
			AtomRecords( i.atomRec ).swap( i.atomRec );
		}
	}
}


/// @brief translates the pose number to pdb numbering string
/// for use in PyRosetta.
/// @param[in] res pose residue number
/// @return pdb string containing chainID and number
PDBInfo::String
PDBInfo::pose2pdb( Size const res ) const
{
	PyAssert((res > 0) && (res <= residue_rec_.size()), "PDBInfo::pose2pdb( Size const res ): res is not in this PDBInfo!" );
	std::stringstream pdb_num, pdb_chain;
	pdb_chain << residue_rec_[res].chainID;
	pdb_num << residue_rec_[res].resSeq;
	return pdb_num.str() + " " + pdb_chain.str() + " ";
}

/// @brief returns for pose( resnumber) the label associated to the residue
/// @note the retrun string is a concatenation of all the strings inside of the vector label<>
/// @param[in] res  residue in pose numbering
utility::vector1 < std::string >
PDBInfo::get_reslabels( Size const res ) const
{
	PyAssert((res > 0) && (res <= residue_rec_.size()), "PDBInfo::get_label( Size const res ): res is not in this PDBInfo!" );
	utility::vector1< std::string > pdb_label;
	for ( core::Size i=1; i<= residue_rec_[res].label.size(); ++i ) {
		pdb_label.push_back(residue_rec_[res].label[i]);
	}
	return pdb_label;
}

/// @brief set chain id for residue
/// @remarks chain id should not be the empty record character, currently '^'
void
PDBInfo::chain(
	Size const res,
	char const chain_id
)
{
	PyAssert((res > 0) && (res <= residue_rec_.size()), "PDBInfo::chain( Size const res, char const chain_id ): res is not in this PDBInfo!" );
	ResidueRecord & rr = residue_rec_[ res ];

	// sync map
	pdb2pose_.conditional_erase( rr.chainID, rr.resSeq, rr.iCode, rr.segmentID, res );
	pdb2pose_.insert( chain_id, rr.resSeq, rr.iCode, rr.segmentID, res );

	// set new residue info
	rr.chainID = chain_id;
}

/// @brief set pdb residue sequence number
void
PDBInfo::number(
	Size const res,
	int const pdb_res
)
{
	PyAssert((res > 0) && (res <= residue_rec_.size()), "PDBInfo::number( Size const res, int const pdb_res ): res is not in this PDBInfo!" );
	ResidueRecord & rr = residue_rec_[ res ];

	// sync map
	pdb2pose_.conditional_erase( rr.chainID, rr.resSeq, rr.iCode, rr.segmentID, res );
	pdb2pose_.insert( rr.chainID, pdb_res, rr.iCode, rr.segmentID, res );

	// set new residue info
	rr.resSeq = pdb_res;
}

/// @brief set insertion code for residue
void
PDBInfo::icode(
	Size const res,
	char const ins_code
)
{
	PyAssert((res > 0) && (res <= residue_rec_.size()), "PDBInfo::icode( Size const res, ins_code ): res is not in this PDBInfo!" );
	ResidueRecord & rr = residue_rec_[ res ];

	// sync map
	pdb2pose_.conditional_erase( rr.chainID, rr.resSeq, rr.iCode, rr.segmentID, res );
	pdb2pose_.insert( rr.chainID, rr.resSeq, ins_code, rr.segmentID, res );

	// set new residue info
	rr.iCode = ins_code;
}

void
PDBInfo::segmentID(
	Size const res,
	std::string const & segmentID
) {
	PyAssert((res > 0) && (res <= residue_rec_.size()), "PDBInfo::icode( Size const res, ins_code ): res is not in this PDBInfo!" );
	ResidueRecord & rr = residue_rec_[ res ];

	// no map sync needed, pray

	// set new residue info
	rr.segmentID = segmentID;

	runtime_assert( segmentID.size() == 4 );
}

/// @brief Displays the PDB info by expressing continuous chain segments
void
PDBInfo::show(
	std::ostream & out
) const
{
	// counters - initialize to the first residue
	Size current_res_tot = 1; // Total residues in this segment
	Size current_atm_tot = natoms(1); // Total atoms in this segment
	Size atm_tot = current_atm_tot; // Overall number of atoms
	Size current_start_pose = 1; // Starting pose number for this segment
	// first line
	out << "PDB file name: " << name() << std::endl;
	out << " Pose Range  Chain    PDB Range  |   #Residues         #Atoms\n" << std::endl;
	for ( Size i=2 ; i <= nres() ; i++ ) {
		// loop through residue records
		if ( icode(i)==' ' && icode(i-1)==' ' && (number(i)-number(i-1))==1 && chain(i)==chain(i-1) ) {
			// Same block - extend it
			current_res_tot += 1;
			current_atm_tot += natoms(i);
			atm_tot += natoms(i);
		} else {
			// resi i is in new segment - print out old segment (to i-1)
			out << I(4,4,current_start_pose) << " -- " << I(4,4,i-1) << A(5,chain(i-1)) <<
				' ' << I(4,4,number(current_start_pose)) << icode(current_start_pose) <<
				" -- " << I(4,4,number(i-1)) << icode(i-1) <<
				" | " << I(6,4,current_res_tot) << " residues; " << I(8,5,current_atm_tot) << " atoms" << std::endl;
			// Reset segment values
			current_res_tot = 1;
			current_atm_tot = natoms(i);
			atm_tot += natoms(i);
			current_start_pose = i;
		}
	}
	//Print out the last segment
	out << I(4,4,current_start_pose) << " -- " << I(4,4,nres()) << A(5,chain(nres())) <<
		' ' << I(4,4,number(current_start_pose)) << icode(current_start_pose) <<
		" -- " << I(4,4,number(nres())) << icode(nres()) <<
		" | " << I(6,4,current_res_tot) << " residues; " << I(8,5,current_atm_tot) << " atoms" << std::endl;

	// last line
	out << "                           TOTAL | " << I(6,4,nres()) << " residues; " << I(8,5,atm_tot) << " atoms" << std::endl;
}

std::string
PDBInfo::short_desc() const {
	utility::vector1< int > res_vector;
	utility::vector1< char > chain_vector;
	utility::vector1< std::string > segid_vector;

	for ( core::Size ii(1); ii <= nres(); ++ii ) {
		res_vector.push_back( number(ii) );
		chain_vector.push_back( chain(ii) );
		segid_vector.push_back( segmentID(ii) );
	}

	return utility::make_tag_with_dashes( res_vector, chain_vector, segid_vector );
}


/// @brief set chain/pdb/insertion code for residue simultaneously
/// @note convenience method; more efficient than doing each individually
///  due to map updates
/// @param[in] res  residue in pose numbering
/// @param[in] chain_id  pdb chain id
/// @param[in] pdb_res  residue in pdb numbering
/// @param[in] ins_code  pdb insertion code
void
PDBInfo::set_resinfo(
	Size const res,
	char const chain_id,
	int const pdb_res,
	char const ins_code,
	std::string const & segmentID
)
{
	ResidueRecord & rr = residue_rec_[ res ];

	// sync map
	pdb2pose_.conditional_erase( rr.chainID, rr.resSeq, rr.iCode, rr.segmentID, res );
	pdb2pose_.insert( chain_id, pdb_res, ins_code, segmentID, res );

	// set new residue info
	rr.chainID = chain_id;
	rr.resSeq = pdb_res;
	rr.iCode = ins_code;
	rr.segmentID = segmentID;
	rr.label.clear();
}


/// @brief adds a label associated to a pose resid.
/// @param[in] res  residue in pose numbering
/// @param[in] label  string that is the "label"
void
PDBInfo::add_reslabel(
	Size res,
	std::string const & label
)
{
	PyAssert((res > 0) && (res <= residue_rec_.size()), "PDBInfo::icode( Size const res, ins_code ): res is not in this PDBInfo!" );
	ResidueRecord & rr = residue_rec_[ res ];
	//Avoid label duplication
	if ( !(std::find( residue_rec_[res].label.begin(), residue_rec_[res].label.end(), label ) != residue_rec_[res].label.end()) ) {
		rr.label.push_back(label);
	}
}

/// @brief clean all the label(s) associated to a pose resid.
/// @param[in] res  residue in pose numbering
void
PDBInfo::clear_reslabel(
	Size const res
)
{
	PyAssert((res > 0) && (res <= residue_rec_.size()), "PDBInfo::icode( Size const res, ins_code ): res is not in this PDBInfo!" );
	ResidueRecord & rr = residue_rec_[ res ];
	rr.label.clear();
}

/// @brief set all residue chain IDs to a single character
void
PDBInfo::set_chains( char const id ) {
	for ( auto i = residue_rec_.begin(), ie = residue_rec_.end(); i < ie; ++i ) {
		i->chainID = id;
	}

	rebuild_pdb2pose();
}


/// @brief copy a section from another PDBInfo
/// @param[in] input_info the PDBInfo to copy from
/// @param[in] copy_from the first residue position in input_info to copy
/// @param[in] copy_to the final residue position in input_info to copy
/// @param[in] start_from the first residue position in this PDBInfo to
///  copy into
void
PDBInfo::copy(
	PDBInfo const & input_info,
	Size const copy_from,
	Size const copy_to,
	Size const start_from
)
{
	debug_assert( copy_from <= input_info.residue_rec_.size() );
	debug_assert( copy_to <= input_info.residue_rec_.size() );
	debug_assert( start_from <= residue_rec_.size() );

	// force erase data from map
	Size idx = start_from;
	ResidueRecords::const_iterator const begin = residue_rec_.begin() + ( start_from - 1 );
	ResidueRecords::const_iterator const end = begin + ( copy_to - copy_from + 1 );
	for ( ResidueRecords::const_iterator i = begin; i < end; ++i, ++idx ) {
		pdb2pose_.erase( i->chainID, i->resSeq, i->iCode, i->segmentID );
	}

	// do the copy
	std::copy(
		input_info.residue_rec_.begin() + ( copy_from - 1 ),
		input_info.residue_rec_.begin() + copy_to,
		residue_rec_.begin() + ( start_from - 1 )
	);

	// add new data to map
	Size count = start_from;
	for ( ResidueRecords::const_iterator i = begin; i < end; ++i, ++count ) {
		pdb2pose_.insert( i->chainID, i->resSeq, i->iCode, i->segmentID, count );
	}
}


/// @brief append residue records after given residue number
/// @param[in] res  residue to append after (in internal/pose numbering)
/// @param[in] natoms  number of atoms in type of appended residue
/// @param[in] n    number of residue records to append
void
PDBInfo::append_res(
	Size const res,
	Size const natoms,
	Size const n
)
{
	ResidueRecord rr;
	rr.atomRec.resize( natoms );
	residue_rec_.insert( residue_rec_.begin() + res, n, rr );

	// sync map -- rebuild from res+2 onwards
	Size count = res + 2;
	ResidueRecords::const_iterator i = residue_rec_.begin() + ( res + 1 );
	for ( ResidueRecords::const_iterator ie = residue_rec_.end(); i < ie; ++i, ++count ) {
		pdb2pose_.insert( i->chainID, i->resSeq, i->iCode, i->segmentID, count );
	}
}


/// @brief prepend residue records before given residue number
/// @param[in] res  residue to prepend before (in internal/pose numbering)
/// @param[in] natoms  number of atoms in type of appended residue
/// @param[in] n    number of residue records to prepend
void
PDBInfo::prepend_res(
	Size const res,
	Size const natoms,
	Size const n
)
{
	ResidueRecord rr;
	rr.atomRec.resize( natoms );
	residue_rec_.insert( residue_rec_.begin() + ( res - 1 ), n, rr );

	// sync map -- rebuild from res+1 onwards
	Size count = res + 1;
	ResidueRecords::const_iterator i = residue_rec_.begin() + res;
	for ( ResidueRecords::const_iterator ie = residue_rec_.end(); i < ie; ++i, ++count ) {
		pdb2pose_.insert( i->chainID, i->resSeq, i->iCode, i->segmentID, count );
	}
}


/// @brief "replace" residue record for given residue number
/// @details Leaves information in residue record untouched, but resizes
///  and zeroes atom records for the residue.
/// @param[in] res residue to replace
/// @param[in] natoms number of atoms in type of residue
void
PDBInfo::replace_res(
	Size const res,
	Size const natoms
)
{
	AtomRecords ar( natoms );
	residue_rec_[ res ].atomRec = ar;

	// no need to sync map, data should stay constant
}

void
PDBInfo::replace_res_remap_bfactors(
	Size const res,
	conformation::Residue const & tgt
)
{
	AtomRecords ar( tgt.natoms() );

	//fpd unfortunately, we don't have access to the source residue at this point.
	//    we'll have to make a reasonable guess at a mapping given the # atoms and the target residue
	if ( residue_rec_[ res ].atomRec.size() == tgt.natoms() ) {
		// #atoms doesn't change, assume repack and copy B factors
		for ( Size i=1; i<=tgt.natoms(); ++i ) {
			ar[i].temperature = residue_rec_[ res ].atomRec[i].temperature;
		}
	} else {
		if ( tgt.is_protein() ) {
			// redesign or cen <-> fa: copy backbone B's, average the rest
			for ( Size i=1; i<=tgt.last_backbone_atom(); ++i ) {
				ar[i].temperature = residue_rec_[ res ].atomRec[i].temperature;
			}

			// set sidechain B factors to the mean
			Real Bsum = 0.0, Bcount = 0.0;
			for ( Size i=tgt.first_sidechain_atom(); i<=residue_rec_[ res ].atomRec.size(); ++i ) {
				Bsum +=  residue_rec_[ res ].atomRec[i].temperature;
				Bcount+=1.0;
			}
			if ( Bcount>0 ) Bsum /= Bcount;
			else Bsum = ar[2].temperature;
			for ( Size i=tgt.first_sidechain_atom(); i<=tgt.nheavyatoms(); ++i ) {
				ar[i].temperature = Bsum;
			}
		} else {
			// ligand (or some other non-protein residue): average all B's
			Real Bsum = 0.0, Bcount = 0.0;
			for ( Size i=1; i<=residue_rec_[ res ].atomRec.size(); ++i ) {
				Bsum +=  residue_rec_[ res ].atomRec[i].temperature;
				Bcount+=1.0;
			}
			if ( Bcount>0 ) Bsum /= Bcount;
			else Bsum = 30.0;
			for ( Size i=tgt.first_sidechain_atom(); i<=tgt.nheavyatoms(); ++i ) {
				ar[i].temperature = Bsum;
			}
		}

		// set H to 1.2 x attached heavy atom
		for ( Size i = 1; i <= tgt.nheavyatoms(); ++i ) {
			Real hB = 1.2*ar[i].temperature;
			for ( Size hid = tgt.attached_H_begin(i), hid_end = tgt.attached_H_end(i);
					hid <= hid_end; ++hid ) {
				ar[hid].temperature = hB;
			}
		}
	}

	residue_rec_[ res ].atomRec = ar;
}


/// @brief delete 'n' residue records starting from given residue
/// @param[in] res  residue to start deleting from (in internal/pose numbering)
/// @param[in] n    number of residue records to delete
void
PDBInfo::delete_res(
	Size const res,
	Size const n
)
{
	auto start = residue_rec_.begin() + ( res - 1 );

	// sync map first (force erase)
	Size idx = res;
	for ( auto i = start, ie = start + n; i < ie; ++i, ++idx ) {
		pdb2pose_.erase( i->chainID, i->resSeq, i->iCode, i->segmentID );
	}

	// delete residue records
	residue_rec_.erase( start, start + n );

	// sync map -- renumber existing records
	for ( Size i = res, ie = residue_rec_.size(); i <= ie; ++i ) {
		ResidueRecord const & rr = residue_rec_[ i ];
		pdb2pose_.insert( rr.chainID, rr.resSeq, rr.iCode, rr.segmentID, i );
	}
}


/// @brief remembers info about all atoms not read into the pose
void
PDBInfo::add_unrecognized_atoms(
	utility::vector1< UnrecognizedAtomRecord > UAs
) {
	unrecognized_atoms_ = UAs;
	num_unrecognized_atoms_ = unrecognized_atoms_.size();
	for ( Size ii = 1; ii <= num_unrecognized_atoms_; ++ii ) {
		unrecognized_res_num2name_[ UAs[ ii ].res_num() ] = UAs[ ii ].res_name();
		if ( unrecognized_res_size_.find( UAs[ ii ].res_num() ) == unrecognized_res_size_.end() ) {
			unrecognized_res_size_[ UAs[ ii ].res_num() ] = 0;
			num_unrecognized_res_++;
		}
		unrecognized_res_size_[ UAs[ ii ].res_num() ] += 1;
	}
}


/// @brief if size of residue records != passed value, fail fast
/// @note This is meant to be used only for en masse methods, not individual
///  residue/atom methods
void
PDBInfo::check_residue_records_size( Size const size ) const
{
	if ( residue_rec_.size() != size ) {
		basic::Error() << "PDBInfo::check_residue_records_size() failed.  An en masse action attempted to access or set more residues than exists in PDBInfo."
			<< std::endl;
		utility_exit();
	}
}


/// @brief rebuilds PDBPoseMap from scratch
void
PDBInfo::rebuild_pdb2pose()
{
	pdb2pose_.clear();
	pdb2pose_.fill( *this );
}

std::ostream & operator << (std::ostream & os, PDBInfo const & info)
{
	//os << "PDB file name: " << info.name() << std::endl;
	info.show(os);
	//os << "Number of Residues: " << info.nres() << std::endl;
	// perhaps add per residue information: pose_num, pdb_num, chain, icode, occupancy, temperature, etc.
	return os;
}

std::ostream & operator << ( std::ostream & os, core::pose::UnrecognizedAtomRecord const & uar )
{
	/* Size res_num,
	std::string res_name,
	std::string atom_name,
	numeric::xyzVector<Real> coords,
	Real temp
	*/
	os << "Res Number " << uar.res_num() << " Atom Name " <<  uar.atom_name() << " Coords " << uar.coords().to_string() << " Temp " << uar.temp() << std::endl;
	return os;
}



} // pose
} // core


#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
//core::pose::UnrecognizedAtomRecord::UnrecognizedAtomRecord() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pose::UnrecognizedAtomRecord::save( Archive & arc ) const {
	arc( CEREAL_NVP( res_num_ ) ); // Size
	arc( CEREAL_NVP( res_name_ ) ); // std::string
	arc( CEREAL_NVP( atom_name_ ) ); // std::string
	arc( CEREAL_NVP( coords_ ) ); // numeric::xyzVector<Real>
	arc( CEREAL_NVP( temp_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pose::UnrecognizedAtomRecord::load( Archive & arc ) {
	arc( res_num_ ); // Size
	arc( res_name_ ); // std::string
	arc( atom_name_ ); // std::string
	arc( coords_ ); // numeric::xyzVector<Real>
	arc( temp_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::UnrecognizedAtomRecord );

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pose::PDBInfo::ResidueRecord::save( Archive & arc ) const {
	arc( CEREAL_NVP( chainID ) ); // char
	arc( CEREAL_NVP( resSeq ) ); // int
	arc( CEREAL_NVP( iCode ) ); // char
	arc( CEREAL_NVP( segmentID ) ); // std::string
	arc( CEREAL_NVP( atomRec ) ); // AtomRecords
	arc( CEREAL_NVP( label ) ); // utility::vector1<std::string>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pose::PDBInfo::ResidueRecord::load( Archive & arc ) {
	arc( chainID ); // char
	arc( resSeq ); // int
	arc( iCode ); // char
	arc( segmentID ); // std::string
	arc( atomRec ); // AtomRecords
	arc( label ); // utility::vector1<std::string>
}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::PDBInfo::ResidueRecord );

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pose::PDBInfo::AtomRecord::save( Archive & arc ) const {
	arc( CEREAL_NVP( isHet ) ); // _Bool
	arc( CEREAL_NVP( altLoc ) ); // char
	arc( CEREAL_NVP( occupancy ) ); // Real
	arc( CEREAL_NVP( temperature ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pose::PDBInfo::AtomRecord::load( Archive & arc ) {
	arc( isHet ); // _Bool
	arc( altLoc ); // char
	arc( occupancy ); // Real
	arc( temperature ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::PDBInfo::AtomRecord );

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pose::PDBInfo::save( Archive & arc ) const {
	arc( CEREAL_NVP( obsolete_ ) ); // _Bool
	arc( CEREAL_NVP( name_ ) ); // String
	arc( CEREAL_NVP( modeltag_ ) ); // String
	arc( CEREAL_NVP( header_information_ ) ); // io::HeaderInformationOP
	arc( CEREAL_NVP( remarks_ ) ); // Remarks
	arc( CEREAL_NVP( residue_rec_ ) ); // ResidueRecords
	arc( CEREAL_NVP( pdb2pose_ ) ); // class core::pose::PDBPoseMap
	arc( CEREAL_NVP( conf_ ) ); // core::conformation::ConformationCAP
	arc( CEREAL_NVP( unrecognized_atoms_ ) ); // utility::vector1<UnrecognizedAtomRecord>
	arc( CEREAL_NVP( unrecognized_res_num2name_ ) ); // std::map<core::Size, std::string>
	arc( CEREAL_NVP( unrecognized_res_size_ ) ); // std::map<core::Size, core::Size>
	arc( CEREAL_NVP( num_unrecognized_res_ ) ); // core::Size
	arc( CEREAL_NVP( num_unrecognized_atoms_ ) ); // core::Size
	arc( CEREAL_NVP( crystinfo_ ) ); // class core::io::CrystInfo
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pose::PDBInfo::load( Archive & arc ) {
	arc( obsolete_ ); // _Bool
	arc( name_ ); // String
	arc( modeltag_ ); // String
	arc( header_information_ ); // io::HeaderInformationOP
	arc( remarks_ ); // Remarks
	arc( residue_rec_ ); // ResidueRecords
	arc( pdb2pose_ ); // class core::pose::PDBPoseMap

	core::conformation::ConformationAP local_conf;
	arc( local_conf ); // core::conformation::ConformationCAP
	conf_ = local_conf;

	arc( unrecognized_atoms_ ); // utility::vector1<UnrecognizedAtomRecord>
	arc( unrecognized_res_num2name_ ); // std::map<core::Size, std::string>
	arc( unrecognized_res_size_ ); // std::map<core::Size, core::Size>
	arc( num_unrecognized_res_ ); // core::Size
	arc( num_unrecognized_atoms_ ); // core::Size
	arc( crystinfo_ ); // class core::io::CrystInfo
}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::PDBInfo );
CEREAL_REGISTER_TYPE( core::pose::PDBInfo )

CEREAL_REGISTER_DYNAMIC_INIT( core_pose_PDBInfo )
#endif // SERIALIZATION


