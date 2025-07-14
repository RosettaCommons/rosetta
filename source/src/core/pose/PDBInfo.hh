// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/PDBInfo.hh
/// @brief  pose information so it's not loose in the pose
/// @author Steven Lewis
/// @author Yih-En Andrew Ban (yab@u.washington.edu)


#ifndef INCLUDED_core_pose_PDBInfo_hh
#define INCLUDED_core_pose_PDBInfo_hh

#include <core/pose/PDBInfo.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

// Unit headers
#include <core/pose/Pose.fwd.hh>

// type headers
#include <core/types.hh>

// Project headers
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/signals/ConnectionEvent.fwd.hh>
#include <core/conformation/signals/IdentityEvent.fwd.hh>
#include <core/conformation/signals/LengthEvent.fwd.hh>

#include <core/pose/PDBPoseMap.hh>

#include <core/io/CrystInfo.hh>
#include <core/io/HeaderInformation.fwd.hh>
#include <core/io/Remarks.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>

// C++ Headers
#include <algorithm>
#include <map>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pose {

/// @brief info about an atom in a unrecognized res (not in pose, but we want to remember it)
class UnrecognizedAtomRecord {

public:
	UnrecognizedAtomRecord();

	UnrecognizedAtomRecord(
		Size res_num,
		std::string const & res_name,
		std::string const & atom_name,
		numeric::xyzVector<Real> coords,
		Real temp
	) : res_num_(res_num),
		res_name_(res_name),
		atom_name_(atom_name),
		coords_(coords),
		temp_(temp)
	{}

	inline core::Size               const & res_num()   const { return res_num_; }
	inline std::string              const & res_name()  const { return res_name_; }
	inline std::string              const & atom_name() const { return atom_name_; }
	inline numeric::xyzVector<Real> const & coords()    const { return coords_; }
	inline Real                     const & temp()      const { return temp_; }

private:

	Size res_num_;
	std::string res_name_, atom_name_;
	numeric::xyzVector<Real> coords_;
	Real temp_;

#ifdef    SERIALIZATION
	//public:
	/// @brief Default constructor for the purpose of serialization only. Do not use.
	// UnrecognizedAtomRecord();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


/// @brief maintains pdb residue & atom information inside a Pose
/// @details Upon creation of new residue records, e.g. when calling the
///  constructors without 'init' or appending/prepending residues, the
///  chain letter for the new records will be set to a character, currently
///  "^", denoting "empty record".  This value may be looked up by
///  calling the static method PDBInfo::empty_record().
/// @remarks Class implementation is biased towards simplicity and fast lookup.
///  Residue/atom information are kept in vectors.  An internally maintained
///  PDBPoseMap provides mapping from pdb -> pose residue numbering.  This
///  causes residue mutators to be a bit more expensive due to map updates,
///  but this is ok because they are typically called sparingly.  Accessors
///  and mutators have overloaded method convention, while special mutators
///  use .set_* convention.
class PDBInfo : public utility::VirtualBase {


public: // typedefs


	typedef utility::VirtualBase Super;

	typedef core::Size Size;
	typedef core::Real Real;
	typedef std::string String;
	typedef core::io::Remarks Remarks;
	typedef core::io::CrystInfo CrystInfo;


public: // forward declarations


	struct AtomRecord;
	struct ResidueRecord;


public: // typedefs


	typedef utility::vector1< AtomRecord > AtomRecords;
	typedef utility::vector1< ResidueRecord > ResidueRecords;


public: // structs


	/// @brief internal struct for storing PDB atom related information
	struct AtomRecord {
		/// @brief default constructor
		AtomRecord() :
			isHet( false ),
			altLoc( ' ' ),
			occupancy( 1.0 ),
			temperature( 0.0 )
		{}

		/// @brief Returns true if there is a heteroatom in the pdb
		bool isHet;
		/// @brief alternate location
		char altLoc;
		/// @brief occupancy
		Real occupancy;
		/// @brief temperature factor
		Real temperature;

#ifdef    SERIALIZATION
		template< class Archive > void save( Archive & arc ) const;
		template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION
	};


	/// @brief internal struct for storing PDB residue related information
	struct ResidueRecord {
		/// @brief default constructor
		ResidueRecord () :
			chainID( PDBInfo::empty_record() ),
			resSeq( 0 ),
			iCode( ' '),
			segmentID( "    " )
		{}

		/// @brief chain -- Normally a letter, but can be longer
		std::string chainID;
		/// @brief residue sequence number
		int resSeq;
		/// @brief insertion code
		char iCode;
		/// @brief segment ID
		std::string segmentID;
		/// @brief vector of AtomRecord
		/// @details sized the same as number of atoms for a given instance of core::conformation::Residue
		AtomRecords atomRec;
		/// @brief Added in 2013 by dadriano. A vector of labels in the AtomRecord that can be used to store
		/// residue-based information that you want/can-use to communicate movers with moverts or afthermath for task-opperations
		utility::vector1< std::string >  label;


#ifdef    SERIALIZATION
		template< class Archive > void save( Archive & arc ) const;
		template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION
	};


public: // constructors


	/// @brief default constructor, obsolete is *true*
	PDBInfo();


	/// @brief size constructor (ensure space for 'n' residue records),
	///  obsolete is *true*
	PDBInfo( Size const n );


	/// @brief Pose constructor (ensures space for residue and atom records
	///  relative to Pose)
	/// @param[in] pose  Pose
	/// @param[in] init  if true (default), then residue records are initialized
	///  and obsolete set to false, otherwise obsolete is true
	///  using Pose residue numbering and chains of the Residues in the Conformation
	PDBInfo(
		Pose const & pose,
		bool init = true
	);


	/// @brief copy constructor
	PDBInfo( PDBInfo const & info );


public: // destructor


	/// @brief default destructor
	~PDBInfo() override;


public: // assignment


	/// @brief copy assignment
	PDBInfo &
	operator =( PDBInfo const & info );


public: // observer interface


	/// @brief Returns the Conformation if this PDBInfo is currently observing
	/// a conformation, otherwise return NULL
	///
	/// example(s):
	///     pose.pdb_info().is_observing()
	/// See also:
	///     Pose
	///     PDBInfo
	core::conformation::ConformationCAP
	is_observing();


	/// @brief Attaches the Conformation  <conf>  and begins observation
	///
	/// example(s):
	///
	/// See also:
	///     Pose
	///     PDBInfo
	void
	attach_to( core::conformation::Conformation & conf );


	/// @brief Detaches the Conformation and stops observation
	/// @remarks takes no arguments because PDBInfo can only observe one
	/// Conformation at a time
	///
	/// example(s):
	///     pose.pdb_info().detach_from()
	/// See also:
	///     Pose
	///     PDBInfo
	void
	detach_from();


	/// @brief Updates when connection to Conformation is changed
	void
	on_connection_change( core::conformation::signals::ConnectionEvent const & event );


	/// @brief Updates atom records when residue identity changes in Conformation
	void
	on_identity_change( core::conformation::signals::IdentityEvent const & event );


	/// @brief Updates residue and atom records when length changes in Conformation,
	/// obsoletes PDBInfo
	void
	on_length_change( core::conformation::signals::LengthEvent const & event );


public: // obsolescence


	/// @brief Returns true if PDBInfo is obsolete and needs updating
	/// @details This flag is currently not used within the class and
	/// is provided for user convenience.  Setting this will
	/// forcibly turn off pdb numbering when dumping pdbs.
	///
	/// example(s):
	///     pose.pdb_info().obsolete()
	/// See also:
	///     Pose
	///     PDBInfo
	inline
	bool
	obsolete() const
	{
		return obsolete_;
	}


	/// @brief Sets the obsolete state to  <flag>
	/// @details this flag is currently not used within the class and
	/// is provided for user convenience.  Setting this will
	/// forcibly turn off pdb numbering when dumping pdbs.
	inline
	void
	obsolete( bool flag )
	{
		obsolete_ = flag;
	}


public: // state


	/// @brief Returns the number of residues represented in PDBInfo
	///
	/// example(s):
	///     pose.pdb_info().nres()
	/// See also:
	///     Pose
	///     PDBInfo
	inline
	Size
	nres() const
	{
		return residue_rec_.size();
	}


	/// @brief Returns the number of atoms represented for the residue  <res>
	/// @note: residue  <res>  must be in pose numbering
	///
	/// example(s):
	///     pose.pdb_info().natoms(3)
	/// See also:
	///     Pose
	///     PDBInfo
	inline
	Size
	natoms( Size const res ) const
	{
		//PyAssert( (res>0) && (res<=residue_rec_.size()), "PDBInfo::natoms( Size const res): res is not in this PDBInfo!" );
		return residue_rec_[ res ].atomRec.size();
	}


	/// @brief Resizes for  <n>  residue records
	/// @details Leaves atom record state inconsistent.  Atom records for
	/// remaining residues are untouched while new residues have no atom
	/// records, so make sure and call one of resize_atom_records()
	/// afterwards if necessary.
	/// @warning Do not use this method for ins/del of residues, as it leaves
	/// the data state inconsistent.  See append_res/prepend_res/delete_res
	/// for that type of functionality.
	void
	resize_residue_records( Size const n );


	/// @brief Ensures  <n>  atom records for residue  <res>  are available
	/// @note: if  <zero>  is true, zero the atom records for this residue
	void
	resize_atom_records(
		Size const res,
		Size const n,
		bool const zero = true
	);


	/// @brief Ensures  <n>  atom records for residue  <res>  are available
	/// @note: if  <zero>  is true, zero the atom records for this residue
	void
	resize_atom_records(
		Size const n,
		bool const zero = true
	);


	/// @brief Updates the number of atom records with respect to atoms in  <pose>
	/// @details Number of internally available atom records will be adjusted
	/// to match number of atoms within each residue in Pose.  Only newly
	/// created records will be zeroed, any existing records are untouched.
	void
	resize_atom_records( Pose const & pose );


	/// @brief Tightens memory usage
	void
	tighten_memory();


	/// @brief Returns the chain letter character specifying "empty record",
	/// currently '^'
	inline
	static
	std::string empty_record() {
		return "^";
	}


public: // pdb-wide accessors/mutators


	/// @brief Returns the pdb name
	///
	/// example(s):
	///     pose.pdb_info().name()
	/// See also:
	///     Pose
	///     PDBInfo
	inline
	String const &
	name() const
	{
		return name_;
	}


	/// @brief Sets the pdb name
	///
	/// example(s):
	///     pose.pdb_info().name('MyPDB')
	/// See also:
	///     Pose
	///     PDBInfo
	inline
	void
	name( String const & s )
	{
		name_ = s;
	}


	/// @brief Returns the model tag for a multi-model pdb
	///
	/// example(s):
	///     pose.pdb_info().modeltag()
	/// See also:
	///     Pose
	///     PDBInfo
	inline
	String const &
	modeltag() const
	{
		return modeltag_;
	}


	/// @brief Sets the model tag for a multi-model pdb to  <tag>
	///
	/// example(s):
	///     pose.pdb_info().modeltag('Number1')
	/// See also:
	///     Pose
	///     PDBInfo
	inline
	void
	modeltag( String const & tag )
	{
		modeltag_ = tag;
	}


	/// @brief Returns the pdb remarks (const)
	///
	/// example(s):
	///     pose.pdb_info().remarks()
	/// See also:
	///     Pose
	///     PDBInfo
	inline
	Remarks const &
	remarks() const
	{
		return remarks_;
	}


	/// @brief Returns the pdb remarks (mutable)
	/// @note we allow direct access to the remarks vector because its
	/// state is independent of the rest of PDBInfo and it's much more
	/// convenient for the user
	///
	/// example(s):
	///     pose.pdb_info().remarks()
	/// See also:
	///     Pose
	///     PDBInfo
	inline
	Remarks &
	remarks()
	{
		return remarks_;
	}


	/// @brief Sets the pdb remarks to  <in>
	///
	/// example(s):
	///
	/// See also:
	///     Pose
	///     PDBInfo
	inline
	void
	remarks( Remarks const & in )
	{
		remarks_ = in;
	}

	/// @brief Returns the pdb crystinfo
	inline
	CrystInfo
	crystinfo() const
	{
		return crystinfo_;
	}

	/// @brief Sets the pdb crystinfo
	inline
	void
	set_crystinfo(CrystInfo crystinfoin)
	{
		crystinfo_ = crystinfoin;
	}


	/// @brief For structures deposited into the protein databank, the
	/// header information record stores the classfication, deposition
	/// date and the 4 character identification code.

	/// For now only allow initalizing it with copy and accessing it
	/// with a constant owning pointer.

	/// Note: The header information is only initialized if it is needed.
	/// Generally this requires using the -run:preserve_header options flag.
	inline
	void
	header_information(io::HeaderInformationOP header_information){
		header_information_ = header_information;
	}

	inline
	io::HeaderInformationCOP
	header_information() const {
		return header_information_;
	}

	inline
	io::HeaderInformationOP
	header_information() {
		return header_information_;
	}


	/// @brief For structures deposited into the protein databank, the
	/// chain sequences stores the SEQRES records. A listing of the
	/// consecutive chemical components covalently linked in a linear
	/// fashion to form a polymer. In short, the sequence (as a vector of
	/// three letter codes) that the researcher performed the experiment on.

	/// Note: The chain sequences information is only initialized if it is needed.
	/// Generally this requires using the -run:preserve_header options flag.

	inline
	void
	chain_sequences( std::map< char, utility::vector1< std::string > > chain_sequences ) {
		chain_sequences_ = chain_sequences;
	}

	inline
	utility::vector1< std::string >
	chain_sequences( char chain ) {
		return chain_sequences_[ chain ];
	}


	inline
	std::map< char, utility::vector1< std::string > >
	chain_sequences() {
		return chain_sequences_;
	}

public: // single residue accessors


	/// @brief Returns the chain letter for pose residue  <res>
	///
	/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
	///
	/// example(s):
	///     pose.pdb_info().chain(3)
	/// See also:
	///     Pose
	///     PDBInfo
	///     PDBInfo.pdb2pose
	///     PDBInfo.pose2pdb
	inline
	std::string const &
	chain( Size const res ) const
	{
		//PyAssert( (res>0) && (res<=residue_rec_.size()), "PDBInfo::chain( Size const res): res is not in this PDBInfo!" );
		return residue_rec_[ res ].chainID;
	}


	/// @brief Returns the pdb sequence number of pose residue  <res>
	///
	/// example(s):
	///     pose.pdb_info().number(3)
	/// See also:
	///     Pose
	///     PDBInfo
	///     PDBInfo.pdb2pose
	///     PDBInfo.pose2pdb
	inline
	int const &
	number( Size const res ) const
	{
		//PyAssert( (res>0) && (res<=residue_rec_.size()), "PDBInfo::number( Size const res): res is not in this PDBInfo!" );
		return residue_rec_[ res ].resSeq;
	}


	/// @brief Returns the pdb insertion code of residue  <res>
	///
	/// example(s):
	///     pose.pdb_info().icode(3)
	/// See also:
	///     Pose
	///     PDBInfo
	///     PDBInfo.pdb2pose
	///     PDBInfo.pose2pdb
	inline
	char const &
	icode( Size const res ) const
	{
		//PyAssert( (res>0) && (res<=residue_rec_.size()), "PDBInfo::icode( Size const res): res is not in this PDBInfo!" );
		return residue_rec_[ res ].iCode;
	}

	/// @brief Returns the segment ID code of residue  <res>
	///
	/// example(s):
	///     pose.pdb_info().icode(3)
	inline
	std::string const &
	segmentID( Size const res ) const
	{
		return residue_rec_[ res ].segmentID;
	}


	/// @brief Returns a tag for residue  <res> like C:1235 (or C:QA:1235 if segment ID is present)
	///
	/// example(s):
	///     pose.pdb_info().icode(3)
	std::string
	tag( Size const res ) const;

	/// @brief Returns the pose numbering of the pdb residue with chain letter <chain>,
	/// pdb residue  <res>, and insertion code  <icode>
	/// @note:  <icode>  is not required, returns 0 if not found,
	/// pose numbering is sequential ie. starts at 1 and goes up
	/// pdb numbering can be anything
	///
	/// example(s):
	///     pose.pdb_info().pdb2pose("B",5)
	/// See also:
	///     Pose
	///     PDBInfo
	///     PDBInfo.pose2pdb
	inline
	Size
	pdb2pose(
		std::string const & chain,
		int const res,
		char const icode = ' ',
		std::string const & segmentID = "    "
	) const
	{
		return pdb2pose_.find( chain, res, icode, segmentID );
	}


	/// @brief Returns the pdb numbering string of pose residue  <res>
	/// @note: pdb string contains the chain letter and number
	/// pose numbering is sequential ie. starts at 1 and goes up
	/// pdb numbering can be anything
	///
	/// example(s):
	///     pose.pdb_info().pose2pdb(25)
	/// See also:
	///     Pose
	///     PDBInfo
	///     PDBInfo.pdb2pose
	String
	pose2pdb( Size const res, std::string const & conn = " " ) const;

	/// @brief Returns the total number of (chain letter) chains in the pose
	///
	/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
	///
	/// example(s):
	///     pose.pdb_info().num_chains()
	/// see also:
	///     pose
	///     pose.annotated_sequence
	///     pose.fold_tree
	///     pose.num_jump
	///   pose.pdb_info().chain
	///   pose.pdb_info().pdb2pose
	int
	num_chains() const;

	/// @brief returns the pose(number) labels associated to the residue
	/// @param[in] res  residue in pose numbering
	utility::vector1 < std::string >
	get_reslabels( Size const res ) const;

public: // single residue mutators


	/// @brief Sets the chain id of pose residue  <res>  to  <chain_id>
	/// @remarks chain_id should not be the empty record character, currently "^"
	///
	/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
	///
	/// example(s):
	///     pose.pdb_info().chain(3,"R")
	/// See also:
	///     Pose
	///     PDBInfo
	///     PDBInfo.pdb2pose
	///     PDBInfo.pose2pdb
	void
	chain(
		Size const res,
		std::string const & chain_id
	);


	/// @brief Sets the pdb sequence residue number of pose residue  <res>  to  <pdb_res>
	/// @brief Returns the pdb insertion code of residue  <res>
	///
	/// example(s):
	///     pose.pdb_info().number(3,81)
	/// See also:
	///     Pose
	///     PDBInfo
	///     PDBInfo.pdb2pose
	///     PDBInfo.pose2pdb
	void
	number(
		Size const res,
		int const pdb_res
	);


	/// @brief Sets the insertion code of pose residue  <res>  to  <ins_code>
	/// @brief Returns the pdb insertion code of residue  <res>
	///
	/// example(s):
	///
	/// See also:
	///     Pose
	///     PDBInfo
	///     PDBInfo.pdb2pose
	///     PDBInfo.pose2pdb
	void
	icode(
		Size const res,
		char const ins_code
	);

	/// @brief Sets the insertion code of pose residue  <res>  to  <ins_code>
	/// @brief Returns the pdb insertion code of residue  <res>
	///
	/// example(s):
	///
	/// See also:
	///     Pose
	///     PDBInfo
	///     PDBInfo.pdb2pose
	///     PDBInfo.pose2pdb
	void
	segmentID(
		Size const res,
		std::string const & segmentID
	);

	/// @brief Sets the chain letter, pdb sequence residue numbering, and insertion
	/// code for pose residue  <res>  to  <chain_id>  ,  <pdb_res>  , and
	/// <ins_code>  respectfully
	/// @note: convenience method; more efficient than doing each individually
	/// due to map updates
	/// @brief Returns the pdb insertion code of residue  <res>
	///
	/// example(s):
	///     pose.pdb_info().icode(3)
	/// See also:
	///     Pose
	///     PDBInfo
	///     PDBInfo.chain
	///     PDBInfo.icode
	///     PDBInfo.number
	///     PDBInfo.pdb2pose
	///     PDBInfo.pose2pdb
	void
	set_resinfo(
		Size const res,
		std::string const & chain_id,
		int const pdb_res,
		char const ins_code = ' ',
		std::string const & segmentID = "    ",
		bool clear_labels = true
	);

	/// @brief adds a label associated to a pose resid.
	/// @param[in] res  residue in pose numbering
	/// @param[in] label  string that is the "label"
	void
	add_reslabel(
		Size res,
		std::string const & label
	);

	/// @brief clean all the label(s) associated to a pose resid.
	/// @param[in] res  residue in pose numbering
	void
	clear_reslabel(
		Size const res
	);

	///@brief Clean a specific label of a specific residue if that label contains
	/// the string.
	void
	clear_reslabel(
		Size const res,
		std::string const & match_label
	);


	/// @brief uses std iterators to check if a residue has a label associated to it
	/// @param[in] res  residue in pose numbering
	/// @param[in] target_label string to look for inside the labes associated to the residue
	bool
	res_haslabel(
		Size res,
		std::string const & target_label
	) const
	{
		return ( std::find( residue_rec_[res].label.begin(), residue_rec_[res].label.end(), target_label ) != residue_rec_[res].label.end() );
	}

	/// @brief parse PDBInfo-LABEL: lines from a collection of lines
	/// @details Lines not containing PDBInfo-LABEL are ignored
	/// @param[in] lines the lines to parse
	void
	parse_pdbinfo_labels( utility::vector1 < std::string > const & lines, Pose const & pose );

	/// @brief write REMARK PDBInfo-LABEL: lines for output to pose
	/// @param[out] remark_lines the output container
	void
	write_pdbinfo_labels( utility::vector1 < std::string > & remark_lines ) const;

public: // atom accessors


	/// @brief Returns true if the  <atom_index>  atom of pose residue  <res>
	/// is a heteratom
	/// @note: atom index within instance of core::conformation::Residue,
	/// currently Rosetta's pdb file output treats the value of .is_het() as an
	/// override -- if this is set to true, the record will be hetatm, otherwise it
	/// will decide by Residue type (the usual default)
	inline
	bool const &
	is_het(
		Size const res,
		Size const atom_index
	) const
	{
		return residue_rec_[ res ].atomRec[ atom_index ].isHet;
	}


	/// @brief Returns the alternate location for the  <atom_index>  atom of pose
	/// residue  <res>
	///
	/// example(s):
	///     pose.pdb_info().alt_loc(1,1)
	/// See also:
	///     Pose
	///     PDBInfo
	///     PDBInfo.pdb2pose
	///     PDBInfo.pose2pdb
	inline
	char const &
	alt_loc(
		Size const res,
		Size const atom_index
	) const
	{
		return residue_rec_[ res ].atomRec[ atom_index ].altLoc;
	}


	/// @brief Returns the occupancy for the  <atom_index>  atom  of pose residue  <res>
	///
	/// example(s):
	///     pose.pdb_info().occupancy(1,1)
	/// See also:
	///     Pose
	///     PDBInfo
	///     PDBInfo.pdb2pose
	///     PDBInfo.pose2pdb
	inline
	Real const &
	occupancy(
		Size const res,
		Size const atom_index
	) const
	{
		return residue_rec_[ res ].atomRec[ atom_index ].occupancy;
	}


	/// @brief Returns the temperature for the  <atom_index>  atom of pose residue  <res>
	///
	/// example(s):
	///     pose.pdb_info().temperature(1,1)
	/// See also:
	///     Pose
	///     PDBInfo
	///     PDBInfo.pdb2pose
	///     PDBInfo.pose2pdb
	inline
	Real const &
	temperature(
		Size const res,
		Size const atom_index
	) const
	{
		return residue_rec_[ res ].atomRec[ atom_index ].temperature;
	}

	/// @brief Returns the temperature for the  <atom_index>  atom of pose residue  <res>
	///
	/// example(s):
	///     pose.pdb_info().temperature(1,1)
	/// See also:
	///     Pose
	///     PDBInfo
	///     PDBInfo.pdb2pose
	///     PDBInfo.pose2pdb
	inline
	Real const &
	bfactor(
		Size const res,
		Size const atom_index
	) const
	{
		return residue_rec_[ res ].atomRec[ atom_index ].temperature;
	}


public: // atom mutators


	/// @brief Sets the heteroatom flag of the  <atom_index>  atom of pose
	/// residue  <res>  to  <flag>
	/// @note currently Rosetta's pdb file output treats the value of .is_het() as an
	/// override -- if this is set to true, the record will be hetatm, otherwise it
	/// will decide by Residue type (the usual default)
	inline
	void
	is_het(
		Size const res,
		Size const atom_index,
		bool const flag
	)
	{
		residue_rec_[ res ].atomRec[ atom_index ].isHet = flag;
	}


	/// @brief Sets the alternate location of the  <atom_index>  atom of pose
	/// residue  <res>  to   <loc>
	/// @note:  <loc> is an alternate location character
	inline
	void
	alt_loc(
		Size const res,
		Size const atom_index,
		char const loc
	)
	{
		residue_rec_[ res ].atomRec[ atom_index ].altLoc = loc;
	}


	/// @brief Sets the occupancy of the  <atom_index>  atom of pose residue
	/// <res>  to  <occ>
	inline
	void
	occupancy(
		Size const res,
		Size const atom_index,
		Real const occ
	)
	{
		residue_rec_[ res ].atomRec[ atom_index ].occupancy = occ;
	}


	/// @brief Sets the temperature of the  <atom_index>  atom of pose residue
	/// <res>  to  <t>
	inline
	void
	temperature(
		Size const res,
		Size const atom_index,
		Real const t
	)
	{
		residue_rec_[ res ].atomRec[ atom_index ].temperature = t;
	}

	/// @brief Sets the temperature of the  <atom_index>  atom of pose residue
	/// <res>  to  <t>
	inline
	void
	bfactor(
		Size const res,
		Size const atom_index,
		Real const t
	)
	{
		residue_rec_[ res ].atomRec[ atom_index ].temperature = t;
	}

public: // residue accessors en masse


	/// @brief Returns the internally maintained PDBPoseMap
	///
	/// example(s):
	///
	/// See also:
	///     Pose
	///     PDBInfo
	///     PDBInfo.pdb2pose
	///     PDBInfo.pose2pdb
	inline
	PDBPoseMap const &
	pdb2pose() const
	{
		return pdb2pose_;
	}

	/// @brief Displays segments of PDB information, segments may or may not
	/// be entire chains
	///
	/// example(s);
	///     pose.pdb_info().show()
	///
	/// See Also:
	///     PDBInfo
	///     PDBInfo.chain
	///     PDBInfo.icode
	///     PDBInfo.nres
	///     PDBInfo.number
	///     Pose
	///     Pose.pdb_info
	void
	show( std::ostream & out ) const;

	/// @brief Returns a short description of the contents of the PDBInfo.
	/// The current format of this is something like "A:1-3 B:9-11"
	///
	/// example(s);
	///     pose.pdb_info().short_desc()
	///
	/// See Also:
	///     PDBInfo
	///     PDBInfo.show
	///     Pose.pdb_info
	std::string
	short_desc() const;

public: // residue mutators en masse


/// @brief Sets the residue chain letters from some string container, iterator version
/// @warning This function does not check if number of elements within span
/// of iterators exceeds number of residues currently defined in object.
/// User must ensure this independently.
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
template< typename StringIterator >
inline
void
set_chains(
	StringIterator const & begin,
	StringIterator const & end
)
{
	auto rr = residue_rec_.begin();

	for ( StringIterator i = begin; i < end; ++i, ++rr ) {
		rr->chainID = *i;
		debug_assert( rr < residue_rec_.end() );
	}

	rebuild_pdb2pose();
}



/// @brief Sets the residue chain letters from some string container, e.g. utility::vector1
/// @details Container must have const .size(), .begin() and .end() iterator
/// access methods.
/// @warning This function checks for size first and will cause failure if
/// size of container does not match number of residues currently defined
/// in object.
///
/// See the documentation of Pose::num_chains() for details about chain numbers, chain letters and jumps.
///
template< typename StringContainer >
inline
void
set_chains( StringContainer const & c )
{
	debug_assert( residue_rec_.size() == c.size() );
	check_residue_records_size( c.size() ); // run-time check

	set_chains( c.begin(), c.end() );
}


	/// @brief Sets the pdb sequence residue numbering from some int container, iterator version
	/// @warning This function does not check if number of elements within span
	/// of iterators exceeds number of residues currently defined in object.
	/// User must ensure this independently.
	template< typename IntIterator >
	inline
	void
	set_numbering(
		IntIterator const & begin,
		IntIterator const & end
	)
	{
		auto rr = residue_rec_.begin();

		for ( IntIterator i = begin; i < end; ++i, ++rr ) {
			rr->resSeq = *i;
			debug_assert( rr < residue_rec_.end() );
		}

		rebuild_pdb2pose();
	}


	/// @brief Sets the pdb sequence residue numbering from some int container, e.g. utility::vector1
	/// @details Container must have const .size(), .begin() and .end() iterator
	/// access methods.
	/// @warning This function checks for size first and will cause failure if
	/// size of container does not match number of residues currently defined
	/// in object.
	template< typename IntContainer >
	inline
	void
	set_numbering( IntContainer const & c )
	{
		debug_assert( residue_rec_.size() == c.size() );
		check_residue_records_size( c.size() ); // run-time check

		set_numbering( c.begin(), c.end() );
	}


	/// @brief Sets the insertion codes from some char container, iterator version
	/// @warning This function does not check if number of elements within span
	/// of iterators exceeds number of residues currently defined in object.
	/// User must ensure this independently.
	template< typename CharIterator >
	inline
	void
	set_icodes(
		CharIterator const & begin,
		CharIterator const & end
	)
	{
		auto rr = residue_rec_.begin();

		for ( CharIterator i = begin; i < end; ++i, ++rr ) {
			rr->iCode = *i;
			debug_assert( rr < residue_rec_.end() );
		}

		rebuild_pdb2pose();
	}


	/// @brief Sets the insertion codes from some char container, e.g. utility::vector1
	/// @details Container must have const .size(), .begin() and .end() iterator
	/// access methods.
	/// @warning This function checks for size first and will cause failure if
	/// size of container does not match number of residues currently defined
	/// in object.
	template< typename CharContainer >
	inline
	void
	set_icodes( CharContainer const & c )
	{
		debug_assert( residue_rec_.size() == c.size() );
		check_residue_records_size( c.size() ); // run-time check

		set_icodes( c.begin(), c.end() );
	}

	template< typename StringIterator >
	inline
	void
	set_segment_ids(
		StringIterator const & begin,
		StringIterator const & end
	) {
		auto rr = residue_rec_.begin();

		for ( StringIterator i = begin; i < end; ++i, ++rr ) {
			rr->segmentID = *i;
			runtime_assert( rr->segmentID.size() == 4 );
			debug_assert( rr < residue_rec_.end() );
		}

		rebuild_pdb2pose();
	}

	template< typename StringContainer >
	inline
	void
	set_segment_ids( StringContainer const & s )
	{
		debug_assert( residue_rec_.size() == s.size() );
		check_residue_records_size( s.size() ); // run-time check

		set_segment_ids( s.begin(), s.end() );
	}

	/// @brief Copyies a section from PDBInfo  <input_info>
	/// @param[in] input_info the PDBInfo to copy from
	/// @param[in] copy_from the first residue position in input_info to copy
	/// @param[in] copy_to the final residue position in input_info to copy
	/// @param[in] start_from the first residue position in this PDBInfo to
	///  copy into
	void
	copy(
		PDBInfo const & input_info,
		Size const copy_from,
		Size const copy_to,
		Size const start_from
	);


public: // residue insertion/deletion


	/// @brief Appends residue records given a pose residue number  <res>
	/// @param[in] res  residue to append after (in internal/pose numbering)
	/// @param[in] natoms  number of atoms in type of appended residue
	/// @param[in] n    number of residue records to append
	void
	append_res(
		Size const res,
		Size const natoms,
		Size const n = 1
	);


	/// @brief Prepends residue records before given pose residue number  <res>
	/// @param[in] res  residue to prepend before (in internal/pose numbering)
	/// @param[in] natoms  number of atoms in type of appended residue
	/// @param[in] n    number of residue records to prepend
	void
	prepend_res(
		Size const res,
		Size const natoms,
		Size const n = 1
	);


	/// @brief "Replaces" residue record for pose residue  <res>
	/// @details Leaves information in residue record untouched, but resizes
	/// and zeroes atom records for the residue.
	/// @param[in] res residue to replace
	/// @param[in] natoms number of atoms in type of residue
	void
	replace_res(
		Size const res,
		Size const natoms
	);

	/// @brief same as replace_res BUT remaps B factors from src to tgt
	void
	replace_res_remap_bfactors(
		Size const res,
		conformation::Residue const & tgt );


	/// @brief Deletes  <n>  residue records starting from pose residue  <res>
	/// @param[in] res  residue to start deleting from (in internal/pose numbering)
	/// @param[in] n    number of residue records to delete
	void
	delete_res(
		Size const res,
		Size const n = 1
	);


	// added by sheffler
	/// @brief remembers info about atoms not read into the pose
	inline
	utility::vector1< UnrecognizedAtomRecord > const &
	get_unrecognized_atoms() const {
		return unrecognized_atoms_;
	}


	core::Size const &
	get_num_unrecognized_atoms() const {
		return num_unrecognized_atoms_;
	}


	inline
	core::Size const &
	get_num_unrecognized_res() const {
		return num_unrecognized_res_;
	}


	inline
	std::string const &
	get_unrecognized_res_name( core::Size const & i ) const {
		return unrecognized_res_num2name_.find(i)->second;
	}


	inline
	core::Size const &
	get_unrecognized_res_size( core::Size const & i ) const {
		return unrecognized_res_size_.find(i)->second;
	}


	// added by sheffler
	/// @brief remembers info about atoms not read into the pose
	void
	add_unrecognized_atoms(
		utility::vector1< UnrecognizedAtomRecord > UAs
	);

	/// @brief rebuilds PDBPoseMap from scratch
	//fpd  making this public because methods in this class don't keep the pdb2pose mapping in sync
	void
	rebuild_pdb2pose();

private: // methods


	/// @brief if size of residue records != passed value, fail fast
	/// @note This is meant to be used only for en masse methods, not individual
	///  residue/atom methods
	void
	check_residue_records_size( Size const size ) const;


private: // data


	/// @brief indicates object is out of sync with reference (e.g. parent Pose)
	/// @details control boolean prevents PDB emitter access if info out of sync
	bool obsolete_;


	/// @brief name of pdb/structure
	String name_;


	/// @brief model tag for multi-model pdbs
	String modeltag_;


	/// @brief header information
	io::HeaderInformationOP header_information_;

	/// @brief chain sequences
	std::map< char, utility::vector1< std::string > > chain_sequences_;

	/// @brief pdb remarks
	Remarks remarks_;


	/// @brief residue records in internal rosetta numbering from 1 .. n
	ResidueRecords residue_rec_;


	/// @brief maps PDB chain,residue -> internal residue numbering
	PDBPoseMap pdb2pose_;


	/// @brief Conformation being observed, NULL if not attached
	core::conformation::ConformationCAP conf_;


	// added by sheffler
	/// @brief information about unrecognized residues
	utility::vector1< UnrecognizedAtomRecord > unrecognized_atoms_;
	std::map< core::Size, std::string > unrecognized_res_num2name_;
	std::map< core::Size, core::Size >  unrecognized_res_size_;
	core::Size num_unrecognized_res_;
	core::Size num_unrecognized_atoms_;

	// fpd spacegroup and crystal parameters
	CrystInfo crystinfo_;


#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //end class PDBInfo

// for Python bindings
std::ostream & operator << ( std::ostream & os, PDBInfo const & info);
std::ostream & operator << ( std::ostream & os, UnrecognizedAtomRecord const & uar );

} // namespace pose
} // namespace core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pose_PDBInfo )
#endif // SERIALIZATION


#endif //INCLUDED_core_pose_PDBInfo_HH
