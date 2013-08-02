// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/PDBInfo.cc
/// @brief  class to hold PDB information so it's not loose in the pose
/// @author Steven Lewis
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// Unit headers
#include <core/pose/PDBInfo.hh>

// Project headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/signals/IdentityEvent.hh>
#include <core/conformation/signals/LengthEvent.hh>
// AUTO-REMOVED #include <core/io/pdb/file_data.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBPoseMap.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/PyAssert.hh>
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
using namespace ObjexxFCL::fmt;


namespace core {
namespace pose {

static basic::Tracer TR("core.pose.PDBInfo");

/// @brief default constructor, obsolete is *true*
PDBInfo::PDBInfo() :
	Super(),
	obsolete_( true ),
	name_( "" ),
	modeltag_( "" ),
	conf_( NULL ),
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
	conf_( NULL ),
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
	residue_rec_( pose.total_residue() ),
	conf_( NULL ),
	num_unrecognized_res_(0),
	num_unrecognized_atoms_(0)
{
	resize_atom_records( pose );

	// initialize using Pose residue numbering and chains of the Residues
	// in the Pose's Conformation
	if ( init ) {
		// do this manually to save on method calls
		for ( Size i = 1, ie = pose.total_residue(); i <= ie; ++i ) {
			ResidueRecord & rr = residue_rec_[ i ];
			//assert( pose.residue( i ).chain() < chr_chains.size() );
			rr.chainID = chr_chains[ pose.residue( i ).chain() % chr_chains.size() ];  //fpd wrap around if > 62 chains
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
	conf_( NULL ), // no observation on copy
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
		if(info.header_information_){
			header_information_ =
				new io::pdb::HeaderInformation(*info.header_information_);
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
core::conformation::Conformation const *
PDBInfo::is_observing() {
	return conf_;
}


/// @brief attach to Conformation and begin observation
void
PDBInfo::attach_to( core::conformation::Conformation & conf ) {
	// detach from prior Conformation if necessary
	if ( conf_ ) {
		detach_from();
	}

	conf.attach_connection_obs( &PDBInfo::on_connection_change, this );
	conf.attach_identity_obs( &PDBInfo::on_identity_change, this );
	conf.attach_length_obs( &PDBInfo::on_length_change, this );

	conf_ = &conf;
}


/// @brief detach from Conformation and stop observation
/// @remarks takes no arguments because PDBInfo can only observe one
///  Conformation at a time
void
PDBInfo::detach_from() {
	if ( conf_ ) {
		conf_->detach_connection_obs( &PDBInfo::on_connection_change, this );
		conf_->detach_identity_obs( &PDBInfo::on_identity_change, this );
		conf_->detach_length_obs( &PDBInfo::on_length_change, this );
	}

	conf_ = NULL;
}


/// @brief update when connection to Conformation is changed
void
PDBInfo::on_connection_change( core::conformation::signals::ConnectionEvent const & event ) {
	using core::conformation::signals::ConnectionEvent;

	switch ( event.tag ) {
		case ConnectionEvent::DISCONNECT:
			obsolete( true );
			detach_from();
			break;
		case ConnectionEvent::TRANSFER:
			// Disconnect -- PDBInfo does not honor TRANSFER tag.
			break;
		default: // do nothing
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
		case IdentityEvent::INVALIDATE: {
			obsolete( true );
			// detach for now, in the future consider whether it's appropriate to
			// rebuild all the data
			detach_from();
			break;
		}
		case IdentityEvent::RESIDUE: {
			// replace and zero the records
			if ( option[ in::preserve_crystinfo ]() ) {
				replace_res_remap_bfactors( event.position, *event.residue ); //fpd  maintain reasonable bfactors
			} else {
				replace_res( event.position, event.residue->natoms() );
			}
		}
		default: {
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
		case LengthEvent::INVALIDATE: {
			obsolete( true );
			detach_from(); // if we don't detach, an array-out-of-bounds may occur on future signals
			break;
		}
		case LengthEvent::RESIDUE_APPEND: {
			append_res( event.position, event.residue->natoms(), abs( event.length_change ) );
			obsolete( true );
			break;
		}
		case LengthEvent::RESIDUE_PREPEND: {
			prepend_res( event.position, event.residue->natoms(), abs( event.length_change ) );
			obsolete( true );
			break;
		}
		case LengthEvent::RESIDUE_DELETE: {
			delete_res( event.position, abs( event.length_change ) );
			//Shouldn't obsolete PDBInfo on delete - every residue that still exists has valid information
			//obsolete( true );
			break;
		}
		default: {
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
		for ( ResidueRecords::iterator i = residue_rec_.begin(), ie = residue_rec_.end(); i != ie; ++i ) {
			i->atomRec = AtomRecords( n );
		}
	} else {
		for ( ResidueRecords::iterator i = residue_rec_.begin(), ie = residue_rec_.end(); i != ie; ++i ) {
			i->atomRec.resize( n );
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

	assert( residue_rec_.size() == pose.n_residue() );

	for ( Size r = 1, re = pose.n_residue(); r <= re; ++r ) {
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
	for ( ResidueRecords::iterator i = residue_rec_.begin(), ie = residue_rec_.end(); i != ie; ++i ) {
		if ( i->atomRec.capacity() > i->atomRec.size() ) {
			AtomRecords( i->atomRec ).swap( i->atomRec );
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
	pdb_num << residue_rec_[res].chainID;
	pdb_chain << residue_rec_[res].resSeq;
	return pdb_chain.str() + " " + pdb_num.str() + " ";
}

/// @brief returns for pose( resnumber) the label associated to the residue
/// @note the retrun string is a concatenation of all the strings inside of the vector label<>
/// @param[in] res  residue in pose numbering
utility::vector1 < std::string >
PDBInfo::get_reslabels( Size const res ) const
{
	PyAssert((res > 0) && (res <= residue_rec_.size()), "PDBInfo::get_label( Size const res ): res is not in this PDBInfo!" );
	utility::vector1< std::string > pdb_label;
	for (core::Size i=1; i<= residue_rec_[res].label.size(); ++i ){
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
	pdb2pose_.conditional_erase( rr.chainID, rr.resSeq, rr.iCode, res );
	pdb2pose_.insert( chain_id, rr.resSeq, rr.iCode, res );

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
	pdb2pose_.conditional_erase( rr.chainID, rr.resSeq, rr.iCode, res );
	pdb2pose_.insert( rr.chainID, pdb_res, rr.iCode, res );

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
	pdb2pose_.conditional_erase( rr.chainID, rr.resSeq, rr.iCode, res );
	pdb2pose_.insert( rr.chainID, rr.resSeq, ins_code, res );

	// set new residue info
	rr.iCode = ins_code;
}

/// @brief Displays the PDB info by expressing continuous chain segments
void
PDBInfo::show(
	std::ostream & out
) const
{
	// counters
	Size current_res_tot = 0;
	Size current_atm_tot = 0;
	Size atm_tot = 0;
	Size current_start_pose = 0;
	char icode_start_pdb(' ');
	// first line
	out << "PDB file name: " << name() << std::endl;
	out << " Pose Range  Chain    PDB Range  | #Resi        #Atoms\n" << std::endl;
	for (Size i=1 ; i <= nres() ; i++){
		// loop through residue records
		if ( i<nres() && icode(i)==' ' && (number(i+1)-number(i))==1 && chain(i)==chain(i+1) ){
			if (current_start_pose<1){	// its empty
				current_res_tot = 1;
				current_atm_tot = natoms(i);
				current_start_pose = i;
				icode_start_pdb = icode(i);
			} else {	// its started, update info
				current_res_tot += 1;
				current_atm_tot += natoms(i);
			}
		} else {	// print it
		out << I(4,4,current_start_pose) << " -- " << I(4,4,i) << A(5,chain(i)) <<
			' ' << I(4,4,number(i)-current_res_tot) << icode_start_pdb << " -- " << I(4,4,number(i)) << icode(i) <<
			" | " << I(6,4,current_res_tot+1) << " residues; " << I(8,5,current_atm_tot) << " atoms\n";
		// then update/blank data
		atm_tot += current_atm_tot;
		current_start_pose = 0;
		}
	}
	// last line
	out << "\t\t\t   TOTAL | " << I(6,4,nres()) << " residues; " << I(8,5,atm_tot) << std::endl;
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
	char const ins_code
)
{
	ResidueRecord & rr = residue_rec_[ res ];

	// sync map
	pdb2pose_.conditional_erase( rr.chainID, rr.resSeq, rr.iCode, res );
	pdb2pose_.insert( chain_id, pdb_res, ins_code, res );

	// set new residue info
	rr.chainID = chain_id;
	rr.resSeq = pdb_res;
	rr.iCode = ins_code;
	rr.label.clear();
}


/// @brief adds a label associated to a pose resid.
/// @param[in] res  residue in pose numbering
/// @param[in] label  string that is the "label"
void
PDBInfo::add_reslabel( 
	Size const res,
	std::string const label 
)
{
	PyAssert((res > 0) && (res <= residue_rec_.size()), "PDBInfo::icode( Size const res, ins_code ): res is not in this PDBInfo!" );
	ResidueRecord & rr = residue_rec_[ res ];
	//Avoid label duplication
	if (!(std::find( residue_rec_[res].label.begin(), residue_rec_[res].label.end(), label ) != residue_rec_[res].label.end())){
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
	for ( ResidueRecords::iterator i = residue_rec_.begin(), ie = residue_rec_.end(); i < ie; ++i ) {
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
	assert( copy_from <= input_info.residue_rec_.size() );
	assert( copy_to <= input_info.residue_rec_.size() );
	assert( start_from <= residue_rec_.size() );

	// force erase data from map
	Size idx = start_from;
	ResidueRecords::const_iterator const begin = residue_rec_.begin() + ( start_from - 1 );
	ResidueRecords::const_iterator const end = begin + ( copy_to - copy_from + 1 );
	for ( ResidueRecords::const_iterator i = begin; i < end; ++i, ++idx ) {
		pdb2pose_.erase( i->chainID, i->resSeq, i->iCode );
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
		pdb2pose_.insert( i->chainID, i->resSeq, i->iCode, count );
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
		pdb2pose_.insert( i->chainID, i->resSeq, i->iCode, count );
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
		pdb2pose_.insert( i->chainID, i->resSeq, i->iCode, count );
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
		for (Size i=1; i<=tgt.natoms(); ++i) {
			ar[i].temperature = residue_rec_[ res ].atomRec[i].temperature;
		}
	} else {
		if (tgt.is_protein()) {
			// redesign or cen <-> fa: copy backbone B's, average the rest
			for (Size i=1; i<=tgt.last_backbone_atom(); ++i)
				ar[i].temperature = residue_rec_[ res ].atomRec[i].temperature;

			// set sidechain B factors to the mean
			Real Bsum = 0.0, Bcount = 0.0;
			for (Size i=tgt.first_sidechain_atom(); i<=residue_rec_[ res ].atomRec.size(); ++i) {
				Bsum +=  residue_rec_[ res ].atomRec[i].temperature;
				Bcount+=1.0;
			}
			if (Bcount>0) Bsum /= Bcount;
			else Bsum = ar[2].temperature;
			for (Size i=tgt.first_sidechain_atom(); i<=tgt.nheavyatoms(); ++i)
				ar[i].temperature = Bsum;
		} else {
			// ligand (or some other non-protein residue): average all B's
			Real Bsum = 0.0, Bcount = 0.0;
			for (Size i=1; i<=residue_rec_[ res ].atomRec.size(); ++i) {
				Bsum +=  residue_rec_[ res ].atomRec[i].temperature;
				Bcount+=1.0;
			}
			if (Bcount>0) Bsum /= Bcount;
			else Bsum = 30.0;
			for (Size i=tgt.first_sidechain_atom(); i<=tgt.nheavyatoms(); ++i)
				ar[i].temperature = Bsum;
		}

		// set H to 1.2 x attached heavy atom
		for (Size i = 1; i <= tgt.nheavyatoms(); ++i) {
			Real hB = 1.2*ar[i].temperature;
			for (Size hid = tgt.attached_H_begin(i), hid_end = tgt.attached_H_end(i);
			          hid <= hid_end; ++hid) {
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
	ResidueRecords::iterator start = residue_rec_.begin() + ( res - 1 );

	// sync map first (force erase)
	Size idx = res;
	for ( ResidueRecords::iterator i = start, ie = start + n; i < ie; ++i, ++idx ) {
		pdb2pose_.erase( i->chainID, i->resSeq, i->iCode );
	}

	// delete residue records
	residue_rec_.erase( start, start + n );

	// sync map -- renumber existing records
	for ( Size i = res, ie = residue_rec_.size(); i <= ie; ++i ) {
		ResidueRecord const & rr = residue_rec_[ i ];
		pdb2pose_.insert( rr.chainID, rr.resSeq, rr.iCode, i );
	}
}


// added by sheffler
/// @brief remembers info about atoms not read into the pose
void
PDBInfo::add_unrecognized_atom(
	Size resnum,
	std::string resname,
	std::string atomname,
	numeric::xyzVector<Real> coords,
	Real temp
) {
	unrecognized_atoms_.push_back( UnrecognizedAtomRecord(resnum,resname,atomname,coords,temp) );
	unrecognized_res_num2name_[resnum] = resname;
	if( unrecognized_res_size_.find(resnum) == unrecognized_res_size_.end() ) {
		unrecognized_res_size_[resnum] = 0;
		num_unrecognized_res_++;
	}
	num_unrecognized_atoms_++;
	unrecognized_res_size_[resnum] += 1;
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

} // pose
} // core
