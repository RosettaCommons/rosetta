// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pdb/file_data.cc
/// @brief  Method definitions for FileData and related classes.
/// @author Sergey Lyskov

// Note: AVOID ACCESSING THE OPTIONS SYSTEM DIRECTLY IN THIS FILE, ESPECIALLY FOR PDB INPUT!
// Doing so will mean the Resource Manager may not work properly.
// Instead, modify FileDataOptions to include the option.


// Unit headers
#include <core/io/pdb/Field.hh>
#include <core/io/pdb/HeaderInformation.hh>
#include <core/io/pdb/file_data.hh>

// Package headers
#include <core/io/pdb/pose_io.hh>
#include <core/io/pdb/file_data_options.hh>
#include <core/io/pdb/file_data_fixup.hh>
#include <core/io/pdb/pdb_dynamic_reader.hh>
#include <core/io/pdb/pdb_dynamic_reader_options.hh>
#include <core/io/pdb/NomenclatureManager.hh>

// Project headers
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/NamedAtomID_Map.hh>
#include <core/io/raw_data/DisulfideFile.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/chemical/types.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/parametric/Parameters.hh>
#include <core/conformation/parametric/ParametersSet.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/pose/PDBInfo.hh>

#include <core/pose/util.hh>
#include <core/pose/ncbb/util.hh>
#include <core/pose/util.tmpl.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/cryst/util.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/exit.hh>

// External headers
#include <ObjexxFCL/format.hh>
#include <boost/lexical_cast.hpp>

// C++ headers
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <utility>


namespace core {
namespace io {
namespace pdb {

using core::Size;
using core::SSize;

using core::chemical::chr_chains;

using basic::T;
using basic::Error;
using basic::Warning;

using ObjexxFCL::strip_whitespace;
using ObjexxFCL::stripped_whitespace;
using ObjexxFCL::rstripped_whitespace;
using namespace ObjexxFCL::format;

using std::string;
using std::iostream;

// Tracer instance for this file
static thread_local basic::Tracer TR( "core.io.pdb.file_data" );


ResidueInformation::ResidueInformation() :
	// resid( "" ),
	resName( "" ),
	chainID( ' ' ),
	resSeq( 0 ),
	iCode( ' ' ),
	terCount( 0 ),
	atoms(),
	xyz(),
	temps()
{}

ResidueInformation::ResidueInformation(
	AtomInformation const & ai) :
	// resid( "" ),
	resName( ai.resName ),
	chainID( ai.chainID ),
	resSeq( ai.resSeq ),
	iCode( ai.iCode ),
	terCount( ai.terCount ),
	atoms(),
	xyz(),
	temps()
{}

bool
ResidueInformation::operator==(
	ResidueInformation const & that) const {
	return
		resName == that.resName &&
		chainID == that.chainID &&
		resSeq == that.resSeq &&
		iCode == that.iCode &&
		terCount == that.terCount;
}

bool
ResidueInformation::operator!=(
	ResidueInformation const & that) const {
	return !(*this == that);
}

String ResidueInformation::resid() const {
	String buf;
	buf.resize(1024);
	// This is horribly hacky. Is this necessary?
	sprintf(&buf[0], "%4d%c%c", this->resSeq, this->iCode, this->chainID);
	buf.resize(6);
	return buf;
}

///////////////////////////////////////////////////////////////////////////////
FileData::FileData() :
	header(/* 0 */),
	remarks(pose::RemarksOP( new pose::Remarks ))
{}

FileData::~FileData()
{}


// Header Information methods /////////////////////////////////////////////////
/// @details prepare the HeaderInformation data structure;
void
FileData::initialize_header_information() {
	header = HeaderInformationOP( new HeaderInformation() );
}

HeaderInformationOP
FileData::header_information() const {
	return header;
}

/// @details Store information in the header record into the HeaderInformation
/// @remarks HeaderInformation must be created explicitly before it can be filled!
void
FileData::store_header_record(Record & R) {
	header->store_record(R);
}

/// @details Populate the header records from the data in the HeaderInformation
/// @remarks HeaderInformation must be created explicitly before it can be filled!
void
FileData::fill_header_records(
	std::vector<Record> & VR
) const {
	header->fill_records(VR);
}

/// @details finalize storing records from the data in the HeaderInformation
/// @remarks HeaderInformation must be created explicitly before it can be filled!
void
FileData::finalize_header_information() {
	header->finalize_parse();
}


// Store (non-standard) polymer linkages in a map. ////////////////////////////
/// @author Labonte
void
FileData::store_link_record( Record & record )
{
	using namespace std;
	using namespace utility;

	LinkInformation link;
	vector1< LinkInformation > links;

	// Extract values from record fields.
	link.name1 = record[ "name1" ].value;  // 1st atom name
	link.resName1 = record[ "resName1" ].value;
	link.chainID1 = record[ "chainID1" ].value[ 0 ];
	link.resSeq1 = atof( record[ "resSeq1" ].value.c_str() );
	link.iCode1 = record[ "iCode1" ].value[ 0 ];

	link.resID1 = record[ "resSeq1" ].value + record[ "iCode1" ].value + record[ "chainID1" ].value;

	link.name2 = record[ "name2" ].value;  // 2nd atom name
	link.resName2 = record[ "resName2" ].value;
	link.chainID2 = record[ "chainID2" ].value[ 0 ];
	link.resSeq2 = atof( record[ "resSeq2" ].value.c_str() );
	link.iCode2 = record[ "iCode2" ].value[ 0 ];

	link.resID2 = record[ "resSeq2" ].value + record[ "iCode2" ].value + record[ "chainID2" ].value;

	link.length = atof( record[ "length" ].value.c_str() );  // bond length

	// If key is found in the links map, add this new linkage information to the links already keyed to this residue.
	if ( link_map.count( link.resID1 ) ) {
		links = link_map[ link.resID1 ];
	}
	links.push_back( link );

	link_map[ link.resID1 ] = links;

	if ( TR.Debug.visible() ) {
		TR.Debug << "LINK record information stored successfully." << std::endl;
	}
}


// Heterogen Information methods //////////////////////////////////////////////
// Store heterogen name information in map.
/// @remarks  heterogen "names" for carbohydrates (from "Rosetta-ready" PDB files) instead have the name field parsed
/// to extract the base (non-variant) ResidueType needed for a particular residue.
/// @author   Labonte
void
FileData::store_heterogen_names(std::string const & hetID, std::string & text)
{
	using namespace std;
	using namespace core::chemical::carbohydrates;

	if ( hetID.empty() ) {
		TR.Warning << "PDB HETNAM record is missing an heterogen ID field." << endl;
		return;
	}
	if ( text.empty() ) {
		TR.Warning << "PDB HETNAM chemical name field is an empty string." << endl;
		return;
	}

	// If the hetID is found in the map of Rosetta-allowed carbohydrate 3-letter codes....
	if ( CarbohydrateInfo::code_to_root_map().count(hetID) ) {
		strip_whitespace(text);
		parse_heterogen_name_for_carbohydrate_residues(text);
	} else {
		// Search through current list of HETNAM records: append or create records as needed.
		bool record_found = false;
		Size const n_heterogen_names = heterogen_names.size();
		for ( uint i = 1; i <= n_heterogen_names; ++i ) {
			// If a record already exists with this hetID, this is a continuation line; append.
			if ( hetID == heterogen_names[i].first ) {
				heterogen_names[i].second.append(rstripped_whitespace(text));
				record_found = true;
				break;
			}
		}
		if ( !record_found ) {
			// Non-carbohydrate heterogen names are simply stored in the standard PDB way.
			strip_whitespace(text);
			heterogen_names.push_back(make_pair(hetID, text));
		}
	}
}

// Parse heterogen name data for a given carbohydrate and save the particular base (non-variant) ResidueType needed in
// a map.
/// @details  The standard PDB HETNAM record is insufficient for indicating the type of carbohydrate residue found at
/// each position of the sequence.  The PDB format only allows one heterogen name per 3-letter code.  To work around
/// this, a PDB file containing carbohydrates will need to be made "Rosetta-ready".  3-letter codes will need to be
/// converted from, e.g., GLC, which implies the vague "ALPHA-D-GLUCOSE", to Glc.  When Rosetta reads in a "sentence
/// case" code such as this, it will check a resID-to-ResidueType map to determine which type of alpha-D-glucose to
/// use, e.g., ->4)-alpha-D-glucopyranosyl or ->6)-alpha-D-glucopyranosyl or ->4)-alpha-D-glucofuranosyl, etc.  This
/// function fills that map from a "Rosetta-ready" HETNAM text field, which includes the resID information (in the
/// same order as in an ATOM or HETATOM record) followed by a space and the base (non-variant) ResidueType.
/// @author Labonte
void
FileData::parse_heterogen_name_for_carbohydrate_residues(std::string const & text)
{
	using namespace std;

	string chainID = string(text.begin(), text.begin() + 1);  // 1 character for chainID
	string resSeq = string(text.begin() + 1, text.begin() + 5);  // 4 characters for resSeq
	string iCode = string(text.begin() + 5, text.begin() + 6);  // 1 character for iCode
	string key = resSeq + iCode + chainID;  // a resID, as defined elsewhere in FileData

	string needed_residue_type_base_name = string(text.begin() + 7, text.end());  // name starts after 7th character

	residue_type_base_names[key] = needed_residue_type_base_name;
}


// Return the PDB resName, chainID, resSeq, and iCode for the given Rosetta sequence position.
ResidueInformation
FileData::get_residue_information(core::pose::Pose const & pose, core::uint const seqpos,
	bool use_PDB, bool renumber_chains) const
{
	using namespace std;
	using namespace utility;
	using namespace core;
	using namespace core::pose;

	ResidueInformation res_info;

	res_info.resName = pose.residue(seqpos).name3();

	// Use PDB-specific information?
	if ( use_PDB ) {
		PDBInfoCOP pdb_info = pose.pdb_info();

		res_info.chainID = pdb_info->chain(seqpos);
		if ( res_info.chainID == PDBInfo::empty_record() ) {  // safety
			TR.Warning << "PDBInfo chain ID was left as character '" << PDBInfo::empty_record()
				<< "', denoting an empty record; for convenience, replacing with space." << endl;
			res_info.chainID = ' ';
		}
		res_info.resSeq = pdb_info->number(seqpos);
		res_info.iCode = pdb_info->icode(seqpos);

		// ...or not?
	} else {
		uint chain_num = pose.chain(seqpos);
		runtime_assert(chain_num > 0);

		res_info.chainID = chr_chains[(chain_num - 1) % chr_chains.size()];
		res_info.resSeq = seqpos;
		res_info.iCode = ' ';

		// If option is specified, renumber per-chain.
		if ( renumber_chains ) {
			vector1<uint> const & chn_ends = pose.conformation().chain_endings();
			for ( uint i = 1; i <= chn_ends.size(); ++i ) {
				if ( chn_ends[i] < seqpos ) {
					res_info.resSeq = seqpos - chn_ends[i];
				}
			}
		}

		// Fix for >10k residues.
		res_info.resSeq %= 10000;
	}

	return res_info;
}

// Pose to FileData methods ///////////////////////////////////////////////////
// Append pdb information to FileData for a single residue.
void
FileData::append_residue(
	core::conformation::Residue const & rsd,
	core::Size & atom_index,
	// for pdb numbering and chains, could change to PDBInfo if necessary (but casting here is perhaps best)
	core::pose::Pose const & pose,
	bool /*preserve_crystinfo*/
)
{
	using namespace core;
	using namespace basic::options;
	using namespace utility;

	// Extract PDBInfo pointer.
	pose::PDBInfoCOP pdb_info = pose.pdb_info();

	// Setup options.
	bool use_PDB(false);
	if (
			pdb_info
			&& !(pdb_info->obsolete())
			&& rsd.seqpos() <= pdb_info->nres()
			&& !(option[ OptionKeys::out::file::renumber_pdb ].value()) ) {
		use_PDB = true;
	}

	bool renumber_chains(false);
	if ( option[ OptionKeys::out::file::per_chain_renumbering ].value() ) {
		renumber_chains = true;
	}


	// Determine residue identifier information.
	ResidueInformation res_info = get_residue_information(pose, rsd.seqpos(), use_PDB, renumber_chains);

	// Generate HETNAM data, if applicable.
	// TODO: For now, only output HETNAM records for saccharide residues, but in the future, outputting HETNAM records
	// for any HETATM residues could be done.
	if ( rsd.is_carbohydrate() ) {
		string const & hetID = rsd.name3();
		string resSeq = print_i("%4d", res_info.resSeq);
		resSeq.resize(4);
		string const resID = string(1, res_info.chainID) + resSeq + string(1, res_info.iCode);

		string const text = resID + " " + residue_type_base_name(rsd.type());

		heterogen_names.push_back(make_pair(hetID, text));
	}


	// Loop through each atom in the residue and generate ATOM or HETATM data.
	for ( Size j = 1; j <= rsd.natoms(); ++j ) {
		//skip outputting virtual atom unless specified
		if ( !option[ OptionKeys::out::file::output_virtual ]() &&
				rsd.atom_type(j).is_virtual() ) continue;

		//fpd optionally don't output centroids
		if ( option[ OptionKeys::out::file::no_output_cen ]() &&
				rsd.atom_name(j) == " CEN" ) continue;

		// skip outputting zero occupancy atoms if specified
		if ( use_PDB && option[ OptionKeys::out::file::suppress_zero_occ_pdb_output ]() &&
				( rsd.seqpos() <= pdb_info->nres() ) ) {
			if ( pdb_info->occupancy( rsd.seqpos(), j ) < 0.0001 ) continue;
		}

		conformation::Atom const & atom( rsd.atom(j) );

		++atom_index;

		AtomInformation ai;
		AtomInformation orb;  //have to initialize this out here.

		ai.isHet = (!rsd.is_polymer() || rsd.is_ligand());
		ai.chainID = res_info.chainID;
		ai.resSeq = res_info.resSeq;
		ai.iCode = res_info.iCode;
		ai.serial = atom_index;
		ai.name = rsd.atom_name(j);
		ai.resName = rsd.name3();
		ai.x = atom.xyz()(1);
		ai.y = atom.xyz()(2);
		ai.z = atom.xyz()(3);
		ai.occupancy = 1.0; // dummy occupancy, can be overridden by PDBInfo

		// Output with pdb-specific info if possible.
		if ( use_PDB ) {
			if ( pdb_info->is_het( rsd.seqpos(), j ) ) { // override standard het only if .is_het() is true
				ai.isHet = true;
			}
			ai.altLoc = pdb_info->alt_loc( rsd.seqpos(), j );
			ai.occupancy = pdb_info->occupancy( rsd.seqpos(), j );
			ai.temperature = pdb_info->temperature( rsd.seqpos(), j );
		}

		// Element
		// (written by fpd; moved here by Labonte)
		core::chemical::AtomTypeSet const &ats = rsd.type().atom_type_set();
		ai.element = ats[atom.type()].element();
		if ( ai.element.length() == 1 ) ai.element = " "+ai.element;

		// 'chains' is member data
		if ( chains.size() < Size(rsd.chain() + 1) ) chains.resize( rsd.chain() + 1 );
		AtomChain & AC(chains[rsd.chain()]);
		AC.push_back(ai);
	}
}

/// @details Read atoms/residue information from Pose object and put it in FileData object.
void
FileData::init_from_pose(core::pose::Pose const & pose)
{
	FileDataOptions options;
	init_from_pose( pose, options );
}

/// @details Read atoms/residue information from Pose object and put it in FileData object using options defined in
/// FileDataOptions.
void
FileData::init_from_pose(core::pose::Pose const & pose, FileDataOptions const & options)
{
	using namespace core;
	using core::pose::PDBInfo;

	// Get Title Section information.
	if ( (options.preserve_header() == true || options.preserve_crystinfo() == true ) && pose.pdb_info() ) {
		*remarks = pose.pdb_info()->remarks();  // Get OP to PDBInfo object for remarks.
		if ( pose.pdb_info()->header_information() ) {
			header = HeaderInformationOP( new HeaderInformation(*(pose.pdb_info()->header_information())) );
		} else {
			header = HeaderInformationOP( new HeaderInformation() );
		}
	}

	// Get parametric information
	if ( options.write_pdb_parametric_info() ) {
		get_parametric_info(remarks, pose);
	}

	// Get Connectivity Annotation Section information.
	if ( options.write_pdb_link_records() ) {
		using namespace utility;
		using namespace id;
		using namespace kinematics;
		using namespace conformation;

		FoldTree const & ft = pose.fold_tree();

		FoldTree::const_iterator end_of_tree = ft.end();
		for ( FoldTree::const_iterator edge = ft.begin(); edge != end_of_tree; ++edge ) {
			if ( edge->is_chemical_bond() ) {
				string const & start_atom = edge->start_atom();
				string const & stop_atom = edge->stop_atom();
				uint const start_num = edge->start();
				uint const stop_num = edge->stop();

				// TODO: Don't assume use of PDBInfo.
				ResidueInformation start_res = get_residue_information(pose, start_num, pose.pdb_info() != NULL );
				ResidueInformation stop_res = get_residue_information(pose, stop_num, pose.pdb_info() != NULL );

				// Fill LinkInformation
				LinkInformation link;
				vector1<LinkInformation> links;

				link.name1 = start_atom;
				link.resName1 = start_res.resName;
				link.chainID1 = start_res.chainID;
				link.resSeq1 = start_res.resSeq;
				link.iCode1 = start_res.iCode;
				link.resID1 = start_res.resid();

				link.name2 = stop_atom;
				link.resName2 = stop_res.resName;
				link.chainID2 = stop_res.chainID;
				link.resSeq2 = stop_res.resSeq;
				link.iCode2 = stop_res.iCode;
				link.resID2 = stop_res.resid();

				// Calculate bond distance.
				uint start_atom_index = pose.residue(start_num).atom_index(start_atom);
				uint stop_atom_index = pose.residue(stop_num).atom_index(stop_atom);
				link.length = pose.conformation().bond_length(
					AtomID(start_atom_index, start_num),
					AtomID(stop_atom_index, stop_num));

				// If key is found in the links map, add this new linkage information to the links already keyed to
				// this residue.
				if ( link_map.count(link.resID1) ) {
					links = link_map[link.resID1];
				}
				links.push_back(link);

				link_map[link.resID1] = links;
			}
		}
	}

	// Get Crystallographic and Coordinate Transformation Section information.
	if ( options.preserve_crystinfo() && pose.pdb_info() ) {
		crystinfo = pose.pdb_info()->crystinfo();
	}

	// Get Coordinate Section information, as well as Heterogen Section information.
	Size const nres( pose.total_residue() );
	Size atom_index(0);

	chains.resize(0);

	for ( Size i=1; i<= nres; ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );
		append_residue( rsd, atom_index, pose, options.preserve_crystinfo() );
	}
}

/// @brief Get parametric information from the Pose object and add it to the PDB remarks.
///
void FileData::get_parametric_info(pose::RemarksOP remarks, core::pose::Pose const & pose) {

	using namespace core::conformation::parametric;

	core::Size const nsets(pose.conformation().n_parameters_sets()); //How many ParametersSet objects are there in the pose?
	if ( nsets==0 ) return; //No need to proceed if this isn't a parametric conformation.

	for ( core::Size iset=1; iset<=nsets; ++iset ) { //Loop through all of the ParametersSet objects.
		ParametersSetCOP curset( pose.conformation().parameters_set(iset) );
		std::stringstream curset_summary;
		curset->get_pdb_remark( curset_summary );

		int cur_remark_number(1); //int instead of core::Size, to match the Remarks class.
		if ( remarks->size() > 0 ) cur_remark_number = (*remarks)[remarks->size()-1].num + 1;

		std::string cur_remark_str;
		while ( std::getline(curset_summary,cur_remark_str) ) {
			pose::RemarkInfo cur_remark;
			cur_remark.value = cur_remark_str; // Ugh.  The RemarkInfo class provides no getters or setters -- only public access to its members.
			cur_remark.num=cur_remark_number;
			remarks->push_back( cur_remark );
			cur_remark_str.clear();
		}
	}

	return;
}

/// @details A lightweight, direct way of limiting pose pdb output to a subset of residues.
/// @note The alternative of constructing new subposes for output only would be unnecessary/less efficient (?)
void
FileData::init_from_pose(
	core::pose::Pose const & pose,
	utility::vector1< core::Size > const & residue_indices
)
{
	using namespace core;
	Size const nres( pose.total_residue() );
	Size atom_index(0);

	chains.resize(0); // 'chains' is member data
	for ( utility::vector1< Size >::const_iterator index( residue_indices.begin() ),
			end( residue_indices.end() ); index != end; ++index ) {
		if ( *index < 1 || *index > nres ) { runtime_assert(false); continue; }
		append_residue( pose.residue( *index ), atom_index, pose );
	}
}


/// @details Convert given Pose object in to PDB format and send it to the given stream.
void
FileData::dump_pdb(
	core::pose::Pose const & pose,
	std::ostream & out,
	string const & /* tag */,
	bool write_fold_tree
)
{
	string data;
	FileData fd;
	fd.init_from_pose(pose);

	data = PDB_DReader::createPDBData(fd);
	out.write( data.c_str(), data.size() );

	write_additional_pdb_data( out, pose, fd, write_fold_tree );
}

/// @details Convert given Pose object into PDB format and save it to 'file_name' file.
/// @return true if operation was completed without error, false otherwise.
bool
FileData::dump_pdb(
	core::pose::Pose const & pose,
	string const & file_name,
	string const & tag,
	bool write_fold_tree)
{
	utility::io::ozstream file(file_name.c_str(), std::ios::out | std::ios::binary);
	if ( !file ) {
		Error() << "FileData::dump_pdb: Unable to open file:" << file_name << " for writing!!!" << std::endl;
		return false;
	}
	dump_pdb(pose, file, tag, write_fold_tree);

	file.close();

	return true;
}

/// @details Convert given Pose object into PDB format and send it to the given stream.
/// only the residues corresponding to indices in 'residue_indices' will be output
void
FileData::dump_pdb(
	core::pose::Pose const & pose,
	std::ostream & out,
	utility::vector1< core::Size > const & residue_indices,
	string const & /* tag */
)
{
	FileData fd;
	string data;
	fd.init_from_pose( pose, residue_indices );

	data = PDB_DReader::createPDBData(fd);
	out.write( data.c_str(), data.size() );

	write_additional_pdb_data( out, pose, fd );
}


/// @details Debug/Info function.
/// Output FileData object to TR like stream in human readable format.
std::ostream&
operator <<(std::ostream &os, FileData const & fd)
{
	os << "<FileData>{";
	for ( Size i=0; i<fd.chains.size(); i++ ) {
		os << "Chain<" << i << ">";
		for ( Size j=0; j<fd.chains[i].size(); j++ ) {
			os << "[" << j << ":" << fd.chains[i][j] << "]" << "\n";
		}
	}
	os << "}";
	return os;
}


/// @details Convert FileData in to set of residues, sequences, coordinates.
/// this is a convenience function, no magic done here.
/// Well, maybe a little.
void
FileData::create_working_data(
	utility::vector1< ResidueInformation > & rinfo
)
{
	FileDataOptions options;
	create_working_data( rinfo, options );
}


/// @details Convert FileData in to set of residues, sequences, coordinates.
/// this is a convenience function, no magic done here.
/// Well, maybe a little.
void
FileData::create_working_data(
	utility::vector1< ResidueInformation > & rinfo,
	FileDataOptions const & options
)
{
	rinfo.clear();

	for ( Size ch=0; ch<chains.size(); ++ch ) {
		for ( Size i=0; i<chains[ch].size(); ++i ) {
			AtomInformation & ai( chains[ch][i] );
			// we should make a copy instead of taking a reference if "fixing" the names causes problems
			std::string const  res_name( convert_res_name( ai.resName ) );
			std::string const atom_name( convert_atom_name( res_name, ai.name ) );
			ai.resName = res_name;
			ai.name = atom_name;

			bool const ok = update_atom_information_based_on_occupancy( ai, options );
			if ( !ok ) continue;

			ResidueInformation new_res( ai );
			if ( rinfo.size() == 0 || rinfo.back() != new_res ) rinfo.push_back(new_res);
			ResidueInformation & curr_res = rinfo.back();
			// Only insert atoms once, so we capture just the first alt conf.
			// Would be nice in the future to take the highest occupancy instead...
			if ( curr_res.xyz.count(ai.name) == 0 ) {
				curr_res.atoms.push_back(ai); // this *does* make a copy
				Vector coords( ai.x, ai.y, ai.z );
				curr_res.xyz[ ai.name ] = coords;
				curr_res.temps[ ai.name ] = ai.temperature;
			}
		}
	}

}


// what to do if occupancy is 0.0?
//chu modify the logic how atoms are treated with zero or negative occupancy field.
bool
FileData::update_atom_information_based_on_occupancy( AtomInformation & ai, FileDataOptions const & options ) const {

	if ( ai.occupancy == 0.0 ) {
		if ( options.randomize_missing_coords() ) {
			randomize_missing_coords( ai );
		} else if ( !options.ignore_zero_occupancy() ) {
			// do nothing and keep this atom as it is
		} else {
			//When flag default changes from true to false, change to TR.Debug and remove second line
			TR.Warning << "PDB reader is ignoring atom " << ai.name << " in residue " << ai.resSeq << ai.iCode << ai.chainID
				<< ".  Pass flag -ignore_zero_occupancy false to change this behavior" << std::endl;
			return false; // skip this atom with zero occ by default
		}
	} else if ( ai.occupancy < 0.0 ) { // always randomize coords for atoms with negative occ
		randomize_missing_coords( ai );
	} else {
		// do nothing for normal atoms with positive occ
	}
	return true;
}


// Helper Functions ///////////////////////////////////////////////////////////

void fill_name_map( core::io::pdb::NameBimap & name_map,
	ResidueInformation const & rinfo,
	chemical::ResidueType const & rsd_type,
	FileDataOptions const & options) {
	//Reset name map
	bool rename = rsd_type.remap_pdb_atom_names();
	for ( core::Size ii(1); ii <= options.residues_for_atom_name_remapping().size(); ++ii ) {
		rename = rename || (options.residues_for_atom_name_remapping()[ ii ] == rsd_type.name3() );
	}
	if ( rename ) {
		// Remap names according to bonding pattern and elements
		// Will reset whatever is in name_map
		remap_names_on_geometry( name_map, rinfo, rsd_type );
	} else {
		name_map.clear();
		// Using names as-is, except for canonicalizing
		for ( utility::vector1< AtomInformation >::const_iterator iter=rinfo.atoms.begin(), iter_end=rinfo.atoms.end(); iter!= iter_end; ++iter ) {
			std::string const & name ( iter->name );
			if ( ! rinfo.xyz.count( name ) ) { // Only map atoms with coordinates.
				continue;
			}
			std::string strip_name( stripped_whitespace(name) );
			if ( rsd_type.has( strip_name ) ) {
				// We do the diversion through the index to canonicalize the atom name spacing to Rosetta standards.
				std::string canonical( rsd_type.atom_name( rsd_type.atom_index(strip_name) ) );
				name_map.insert( NameBimap::value_type( name, canonical ) );
			}
		}
	}
}

/// @details The missing density regions in the input pdb should have 0.000 in the placeholders
/// this routine puts random coordinates wherever there is 0.000 for mainchain atoms.
/// tex - that's a stupid way of defining missing density, as atoms can be at the origin for other
/// reasons. This has been updated to check for occupancy to define missing density rather than atoms
/// located at the origin.
void
FileData::randomize_missing_coords( AtomInformation & ai ) const {
	// if( ai.resSeq == 1 && ai.name == " N  ") return;//ignore first atom. Rosetta pdbs start with 0.000
	if ( ai.x == 0.000 && ai.y == 0.000 && ai.z == 0.000 && ai.occupancy <= 0.0 ) {
		TR << "Randomized: " << ai.name << " " << ai.resName << "  " << ai.resSeq << std::endl;
		//v  if ( ai.name == " N  " || ai.name == " CA " || ai.name == " C  " ||
		//v   ai.name == " O  " || ai.name == " CB " ) {
		ai.x = ai.x + 900.000 + numeric::random::rg().uniform()*100.000;
		ai.y = ai.y + 900.000 + numeric::random::rg().uniform()*100.000;
		ai.z = ai.z + 900.000 + numeric::random::rg().uniform()*100.000;
		//v  }
	}
	return;
}


// Adds data to the end of a pdb that are not a standard part of the pdb format.
void
write_additional_pdb_data(
	std::ostream & out,
	pose::Pose const & pose,
	io::pdb::FileData const &,
	bool write_fold_tree
)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// added by rebecca --> Normalized MEM lines. Useful for visualizing the boundaries of
	// the membrane by coupling the NORM and THK coordinates
	bool normalize_to_thk( true );
	if ( option[ OptionKeys::mp::output::normalize_to_thk ].user() ) {
		normalize_to_thk = option[ OptionKeys::mp::output::normalize_to_thk ]();
	}

	if ( pose.conformation().is_membrane() && normalize_to_thk == true ) {

		// Grab membrane residue & current data
		core::Size resid( pose.conformation().membrane_info()->membrane_rsd_num() );
		core::Real thkn( pose.conformation().membrane_info()->membrane_thickness() );
		core::Vector cntr( pose.conformation().membrane_info()->membrane_center() );
		core::Vector norm( pose.conformation().membrane_info()->membrane_normal() );

		// Actually normalize the membrane residue to thk
		norm.normalize( thkn );

		// MEM normalized to 15 has new chain
		char curr_chain = pose.pdb_info()->chain( pose.total_residue() );
		char new_chain = (char)((int) curr_chain + 1);

		// Get rsdid, current chain,
		out << "HETATM XXXX THKN MEM " << new_chain << I(4,resid+1) << "    " << F(8, 3, thkn) << "   0.000   0.000 \n";
		out << "HETATM XXXX CNTR MEM " << new_chain << I(4,resid+1) << "    " << F(8, 3, cntr.x()) << F(8, 3, cntr.y()) << F(8, 3, cntr.z()) << "\n";
		out << "HETATM XXXX NORM MEM " << new_chain << I(4,resid+1) << "    " << F(8, 3, norm.x()) << F(8, 3, norm.y()) << F(8, 3, norm.z()) << "\n";

		//              HETATM XXXX THKN MEM X  81      15.000   0.000   0.000
		//    HETATM XXXX CNTR MEM X  81       0.000   0.000   0.000
		//              HETATM XXXX NORM MEM X  81       0.000   0.000  15.000

	}

	// added by rhiju --> "CONECT" lines. Useful for coarse-grained/centroid poses, so that
	//  rasmol/pymol draws bonds between atoms 'bonded' in Rosetta that are far apart.
	//  Can also be turned on with -inout:dump_conect_info flag.
	// CONECT reading and writing should really be handled by FileData. ~Labonte
	if ( (pose.total_residue() >= 1 && pose.residue(1).is_coarse()) ||
			option[ OptionKeys::inout::dump_connect_info]() ) {
		dump_connect_info( pose, out );
	}

	if ( write_fold_tree || option[ OptionKeys::inout::fold_tree_io ].user() ) {
		out << "REMARK " << pose.fold_tree();
	}
	if ( basic::options::option[ OptionKeys::out::file::pdb_parents]() ) {
		std::string value;
		bool has_parents = core::pose::get_comment( pose, "parents", value );
		if ( has_parents ) {
			out << "REMARK PARENT    " << value.substr(0,5) << std::endl;
		}
	}
	if ( basic::options::option[ OptionKeys::out::file::pdb_comments]() ) {
		out << "##Begin comments##" << std::endl;
		using namespace std;
		map< string, string > const comments = core::pose::get_all_comments(pose);
		for ( std::map< string, string >::const_iterator i = comments.begin(); i != comments.end(); ++i ) {
			out << i->first<<" "<<i->second << std::endl;
		}
		out << "##End comments##" << std::endl;
	}

	if ( basic::options::option[ basic::options::OptionKeys::out::file::output_orbitals] ) {
		static std::string const chains( " ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890" );
		for ( core::Size i=1; i <=pose.n_residue(); ++i ) {
			core::conformation::Residue rsd(pose.residue(i));
			core::Size number(0);
			char const chain( chains[ rsd.chain() ] );
			for ( core::Size j=1; j<=rsd.natoms(); ++j ) {
				if ( rsd.atom_type(j).atom_has_orbital() ) {
					utility::vector1<core::Size> const & orbital_indices(rsd.bonded_orbitals(j));
					for (
							utility::vector1<core::Size>::const_iterator
							orbital_index = orbital_indices.begin(),
							orbital_index_end = orbital_indices.end();
							orbital_index != orbital_index_end; ++orbital_index
							) {
						++number;
						Vector orbital_xyz(rsd.orbital_xyz(*orbital_index));
						out << "ATOM  " << I(5,number) << ' ' << rsd.orbital_name(*orbital_index) << ' ' <<
							rsd.name3() << ' ' << chain << I(4,rsd.seqpos() ) << "    " <<
							F(8,3,orbital_xyz.x()) <<
							F(8,3,orbital_xyz.y()) <<
							F(8,3,orbital_xyz.z()) <<
							F(6,2,1.0) << F(6,2,1.0) << '\n';
					}
				}
			}
		}
	}
	if ( basic::options::option[ basic::options::OptionKeys::out::file::output_torsions ] ) {
		if ( !core::pose::is_ideal_pose(pose) ) {
			TR << "Ignoring out::file::output_torsions option because pose is non-ideal!" << std::endl;
		} else {
			ObjexxFCL::FArray1D_char dssp_reduced_secstruct(pose.n_residue());
			scoring::dssp::Dssp(pose).dssp_reduced(dssp_reduced_secstruct);
			out << "REMARK torsions: res pdbres pdbchain seq dssp phi psi omega" << std::endl;
			for ( core::Size i=1; i<=pose.n_residue(); ++i ) {
				out << "REMARK " << I( 4, i ) << " " << I( 4, pose.pdb_info()->number(i)) << " " << pose.pdb_info()->chain(i) << " " << pose.residue( i ).name1() << " " <<
					dssp_reduced_secstruct(i) << " " << F( 9, 3, pose.phi(i)) << " " << F( 9, 3, pose.psi(i)) << " " << F( 9, 3, pose.omega(i)) << std::endl;
			}
		}
	}

	// Added by Daniel-Adriano Silva, used to write the PDBInfoLabels to the REMARK
	// First test that the pdb_info() is not empty
	if ( pose.pdb_info() ) {
		// Then output the labels
		for ( core::Size i=1; i<=pose.n_residue(); ++i ) {
			utility::vector1 < std::string > tmp_v_reslabels =  pose.pdb_info()->get_reslabels(i);
			core::Size numLables=tmp_v_reslabels.size();
			//Only write if the residue has any label (keep the file as small as possible)
			if ( numLables > 0 ) {
				out << "REMARK PDBinfo-LABEL: " << I( 4, i );
				for ( core::Size lndx=1; lndx <= numLables; ++lndx ) {
					out << " " << tmp_v_reslabels[lndx];
				}
				out << std::endl;
			}
		}
	}
}


// FileData to Pose ///////////////////////////////////////////////////////////

void
build_pose_from_pdb_as_is(
	pose::Pose & pose,
	std::string const & filename
)
{
	PDB_DReaderOptions options;
	build_pose_from_pdb_as_is( pose, filename, options );
}

void
build_pose_from_pdb_as_is(
	pose::Pose & pose,
	std::string const & filename,
	PDB_DReaderOptions const & pdr_options
)
{
	using namespace chemical;
	build_pose_from_pdb_as_is( pose, * ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ), filename, pdr_options );
}

void
build_pose_from_pdb_as_is(
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename
)
{
	PDB_DReaderOptions options;
	build_pose_from_pdb_as_is( pose, residue_set, filename, options );
}

void
build_pose_from_pdb_as_is(
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename,
	PDB_DReaderOptions const & pdr_options
)
{
	utility::io::izstream file( filename );
	if ( !file ) {
		TR.Error << "File:" << filename << " not found!" << std::endl;
		utility_exit_with_message( "Cannot open file " + filename );
	} else {
		TR.Debug << "read file: " << filename << std::endl;
	}
	build_pose_from_pdb_as_is( pose, residue_set, filename, file, pdr_options );
}

void
build_pose_from_pdb_as_is(
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename,
	std::istream & file_contents,
	PDB_DReaderOptions const & pdr_options
)
{
	std::string all_lines;
	utility::slurp( file_contents, all_lines );
	FileData fd = PDB_DReader::createFileData( all_lines, pdr_options );
	if ( fd.filename == "" ) {
		fd.filename = filename;
	}
	id::AtomID_Mask missing( false );
	build_pose_as_is1( fd, pose, residue_set, missing, pdr_options);
}

// Fills the pose with the data from FileData
void
build_pose_as_is1(
	io::pdb::FileData & fd,
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	id::AtomID_Mask & missing,
	FileDataOptions const & options
) {
	typedef std::map< std::string, double > ResidueTemps;
	typedef std::map< std::string, Vector > ResidueCoords;
	typedef utility::vector1< std::string > Strings;

	using namespace chemical;
	using namespace conformation;

	// reset current data
	pose.clear();

	utility::vector1< ResidueInformation > rinfos;
	id::NamedAtomID_Mask coordinates_assigned( false );
	// Map pose residue numbers to indices into rinfos.
	// Some residues in the input file may be discarded (missing atoms, unrecognized, etc.)
	utility::vector1< Size > pose_to_rinfo;
	fd.create_working_data( rinfos, options );
	fixup_rinfo_based_on_residue_type_set( rinfos, residue_set );
	utility::vector1<ResidueTemps> pose_temps;

	Strings branch_lower_termini;

	Size const nres_pdb( rinfos.size() );

	// Map rinfo atom names to Rosetta pose atom names (and pose->rinfo for the right map)
	utility::vector1< core::io::pdb::NameBimap > rinfo_name_map(nres_pdb);

	utility::vector1<Size> UA_res_nums;
	utility::vector1<std::string> UA_res_names, UA_atom_names;
	utility::vector1<numeric::xyzVector<Real> > UA_coords;
	utility::vector1<core::Real> UA_temps;

	std::string chains_whose_residues_are_separate_chemical_entities =
		options.chains_whose_residues_are_separate_chemical_entities();
	std::string::const_iterator const entities_begin = chains_whose_residues_are_separate_chemical_entities.begin();
	std::string::const_iterator const entities_end = chains_whose_residues_are_separate_chemical_entities.end();

	std::string chains_to_check_if_Ntermini= options.check_if_residues_are_Ntermini() ;
	std::string::const_iterator const check_Ntermini_begin = chains_to_check_if_Ntermini.begin();
	std::string::const_iterator const check_Ntermini_end = chains_to_check_if_Ntermini.end();

	std::string chains_to_check_if_Ctermini= options.check_if_residues_are_Ctermini() ;
	std::string::const_iterator const check_Ctermini_begin = chains_to_check_if_Ctermini.begin();
	std::string::const_iterator const check_Ctermini_end = chains_to_check_if_Ctermini.end();

	//mjo do not add residue by bond if the last residue was not recognized
	bool last_residue_was_recognized(true);

	// Loop over every residue in the FileData extracted from the PDB file, select appropriate ResidueTypes,
	// create Residues, and build the Pose.
	for ( Size i = 1; i <= nres_pdb; ++i ) {
		ResidueInformation const & rinfo = rinfos[ i ];
		char chainID = rinfo.chainID;
		std::string const & pdb_name = rinfo.resName;
		std::string const & resid = rinfo.resid();

		runtime_assert( resid.size() == 6 );

		// Convert PDB 3-letter code to Rosetta 3-letter code, if a list of alternative codes has been provided.
		std::pair< std::string, std::string > const & rosetta_names(
			NomenclatureManager::get_instance()->rosetta_names_from_pdb_code( pdb_name ) );
		std::string const & name3( rosetta_names.first );
		if ( rosetta_names.second != "" ) {
			fd.residue_type_base_names[ resid ] = rosetta_names.second;
		}

		bool const separate_chemical_entity = find(entities_begin, entities_end, chainID ) !=  entities_end;
		bool const same_chain_prev = ( i > 1        && chainID == rinfos[i-1].chainID &&
			rinfo.terCount == rinfos[i-1].terCount && !separate_chemical_entity);
		bool const same_chain_next = ( i < nres_pdb && chainID == rinfos[i+1].chainID &&
			rinfo.terCount == rinfos[i+1].terCount && !separate_chemical_entity);
		bool const check_Ntermini_for_this_chain = ("ALL" == chains_to_check_if_Ntermini) ?
			true : find(check_Ntermini_begin, check_Ntermini_end, chainID ) ==  check_Ntermini_end;
		bool const check_Ctermini_for_this_chain = ("ALL" == chains_to_check_if_Ctermini) ?
			true : find(check_Ctermini_begin, check_Ctermini_end, chainID ) ==  check_Ctermini_end;

		// Determine polymer information: termini, branch points, etc.
		Strings branch_points_on_this_residue;
		bool is_branch_point( false );
		if ( TR.Trace.visible() ) {
			TR.Trace << "Checking if resid " << resid << " is in the link map " << std::endl;
		}
		if ( fd.link_map.count( resid ) ) {  // if found in the linkage map
			// Find and store to access later:
			//     - associated 1st residue of all branches off this residue (determines branch lower termini)
			//     - positions of branch points
			if ( TR.Trace.visible() ) {
				TR.Trace << "Found resid " << resid << " in link map " << std::endl;
			}
			for ( Size branch = 1, n_branches = fd.link_map[ resid ].size(); branch <= n_branches; ++branch ) {

				if ( TR.Trace.visible() ) {
					TR.Trace << "Examining branch " << branch << std::endl;
				}
				LinkInformation const & link_info( fd.link_map[ resid ][ branch ] );
				if ( link_info.chainID1 == link_info.chainID2 && link_info.resSeq1 == ( link_info.resSeq2 - 1 ) ) {
					// If this occurs, the link is to the next residue on the same chain, so both residues are part of
					// the same main chain or branch, and this linkage information can be ignored.
					// Note that this assumes insertion codes are not involved! It also assumes that the PDB file
					// makers did things reasonably.
					continue;
				}
				if ( TR.Trace.visible() ) {
					TR.Trace << "branch point is true for this residue" << std::endl;
					TR.Trace << "corresponding branch lower terminus is " <<  link_info.resID2 << std::endl;
				}
				is_branch_point = true;
				branch_lower_termini.push_back( link_info.resID2 );
				branch_points_on_this_residue.push_back( link_info.name1 );
			}
		}
		bool const is_branch_lower_terminus = branch_lower_termini.contains(resid);
		bool const is_lower_terminus( ( i == 1 || rinfos.empty() || (!same_chain_prev && !is_branch_lower_terminus) )
			&& check_Ntermini_for_this_chain );
		bool const is_upper_terminus( ( i == nres_pdb || !same_chain_next ) && check_Ctermini_for_this_chain );

		ResidueCoords const & xyz = rinfo.xyz;
		ResidueTemps  const & rtemp = rinfo.temps;

		bool is_d_aa = false;
		bool is_l_aa = false;
		// Skip the below for canonicals
		// We should really have this relationship wrapped in a function.
		if ( name3 == "DAL" || name3 == "DCS" || name3 == "DAS" || name3 == "DGU" ||
				name3 == "DPH" || name3 == "DHI" || name3 == "DIL" || name3 == "DLY" ||
				name3 == "DLE" || name3 == "DME" || name3 == "DAN" || name3 == "DPR" ||
				name3 == "DGN" || name3 == "DAR" || name3 == "DSE" || name3 == "DTH" ||
				name3 == "DVA" || name3 == "DTR" || name3 == "DTY" ) {
			is_d_aa = true;
		} else if ( name3 == "ALA" || name3 == "CYS" || name3 == "CYD" || name3 == "ASP" || name3 == "GLU" ||
				name3 == "PHE" || name3 == "HIS" || name3 == "ILE" ||
				name3 == "LYS" || name3 == "LEU" || name3 == "MET" || name3 == "ASN" ||
				name3 == "PRO" || name3 == "GLN" || name3 == "ARG" || name3 == "SER" ||
				name3 == "THR" || name3 == "VAL" || name3 == "TRP" || name3 == "TYR" ) {
			is_l_aa = true;
		} else {
			// Better method please.
			// AMW: revised so it doesn't rely on hydrogens being present! d'oh.
			// Explicitly exclude peptoids.
			if ( xyz.find( " CA " ) != xyz.end() && xyz.find( " CA1" ) == xyz.end() ) {
				// There are four atoms bonded to CA.
				if ( xyz.find( " Pbb" ) != xyz.end() ) {
					// Phosphonate
					core::Real characteristic_angle = numeric::dihedral_degrees( xyz.at( " N  " ), xyz.at( " Pbb" ), xyz.at( " CB " ), xyz.at( " CA " ) );
					if ( characteristic_angle > 0 ) {
						is_d_aa = true;
					} else {
						is_l_aa = true;
					}
				} else if ( xyz.find( " CM " ) != xyz.end() && name3 != "MLZ" ) { // methyllysine also uses CM, illustrating the weakness of this method.
					// beta
					if ( xyz.find( " CB " ) != xyz.end() ) {
						core::Real characteristic_angle = numeric::dihedral_degrees( xyz.at( " N  " ), xyz.at( " CM " ), xyz.at( " CB " ), xyz.at( " CA " ) );
						// Positive angles are D
						//std::cout << "Characteristic 1b angle: " << characteristic_angle << std::endl;
						if ( characteristic_angle > 0 ) {
							is_d_aa = true;
						} else {
							is_l_aa = true;
						}
					} else if ( xyz.find( " CB1" ) != xyz.end() && xyz.find( " CB2" ) != xyz.end() ) {
						core::Real characteristic_angle = numeric::dihedral_degrees( xyz.at( " N  " ), xyz.at( " CM " ), xyz.at( " CB1" ), xyz.at( " CB2" ) );
						// Positive angles are D
						//std::cout << "Characteristic 2b angle: " << characteristic_angle << std::endl;
						if ( characteristic_angle > 0 ) {
							is_d_aa = true;
						} else {
							is_l_aa = true;
						}
					} // other possibilities: B3G
				} else {
					// alpha
					if ( xyz.find( " CB " ) != xyz.end() ) {
						core::Real characteristic_angle = numeric::dihedral_degrees( xyz.at( " N  " ), xyz.at( " C  " ), xyz.at( " CB " ), xyz.at( " CA " ) );
						// Positive angles are D
						//std::cout << "Characteristic 1a angle: " << characteristic_angle << std::endl;
						if ( characteristic_angle > 0 ) {
							is_d_aa = true;
						} else {
							is_l_aa = true;
						}
					} else if ( xyz.find( " CB1" ) != xyz.end() && xyz.find( " CB2" ) != xyz.end() ) {
						// CB1 is designated the L configuration controller.
						core::Real characteristic_angle = numeric::dihedral_degrees( xyz.at( " N  " ), xyz.at( " C  " ), xyz.at( " CB1" ), xyz.at( " CA " ) );
						// Positive angles are D
						//std::cout << "Characteristic 2a angle: " << characteristic_angle << std::endl;
						if ( characteristic_angle > 0 ) {
							is_d_aa = true;
						} else {
							is_l_aa = true;
						}
					} // other possibilities: GLY
				}
			} else if ( xyz.find( " C2 " ) != xyz.end() && xyz.find( " C3 " ) != xyz.end() && xyz.find( " C4 " ) != xyz.end() ) {
				// If we have a gamma, assign based on the stereo of the first carbon from C
				// This is an if else if NOT because we expect these to be mutually exclusive
				// but because CB3 only matters if there is no CB2.
				// We are NOT handling disubstituted.
				// AMW: Vikram's params maintain the position of the 2 in the name...
				if ( xyz.find( "CB2 " ) != xyz.end() ) {
					// This determines everything!
					core::Real characteristic_angle = numeric::dihedral_degrees( xyz.at( " C3 " ), xyz.at( " C  " ), xyz.at( "CB2 " ), xyz.at( " C2 " ) );
					if ( characteristic_angle > 0 ) {
						is_d_aa = true;
					} else {
						is_l_aa = true;
					}

				} else if ( xyz.find( "CB3 " ) != xyz.end() ) {
					// This determines everything!
					core::Real characteristic_angle = numeric::dihedral_degrees( xyz.at( " C4 " ), xyz.at( " C2 " ), xyz.at( "CB3 " ), xyz.at( " C3 " ) );
					if ( characteristic_angle > 0 ) {
						is_d_aa = true;
					} else {
						is_l_aa = true;
					}

				} else if ( xyz.find( "CB4 " ) != xyz.end() ) {
					// This determines everything!
					core::Real characteristic_angle = numeric::dihedral_degrees( xyz.at( " N  " ), xyz.at( " C3 " ), xyz.at( "CB4 " ), xyz.at( " C4 " ) );
					if ( characteristic_angle > 0 ) {
						is_d_aa = true;
					} else {
						is_l_aa = true;
					}
				}
			}
		}

		// Get a list of ResidueTypes that could apply for this particular 3-letter PDB residue name.
		if ( ! is_residue_type_recognized(
				i, name3, residue_set, xyz, rtemp,
				UA_res_nums, UA_res_names, UA_atom_names, UA_coords, UA_temps, options) ) {
			last_residue_was_recognized = false;
			continue;
		}

		if ( TR.Debug.visible() ) {
			TR.Debug << "Residue " << i << std::endl;
		}
		if ( TR.Trace.visible() ) {
			TR.Trace << "...same_chain_prev: " << same_chain_prev << std::endl;
			TR.Trace << "...same_chain_next: " << same_chain_next << std::endl;
			TR.Trace << "...is_lower_terminus: " << is_lower_terminus << std::endl;
			TR.Trace << "...check_Ntermini_for_this_chain: "<< check_Ntermini_for_this_chain << std::endl;
			TR.Trace << "...is_upper_terminus: " << is_upper_terminus << std::endl;
			TR.Trace << "...check_Ctermini_for_this_chain: "<< check_Ctermini_for_this_chain << std::endl;
			TR.Trace << "...is_branch_point: " << is_branch_point << std::endl;
			TR.Trace << "...is_branch_lower_terminus: "<< is_branch_lower_terminus << std::endl;
			TR.Trace << "...last_residue_was_recognized: " << last_residue_was_recognized << std::endl;
			TR.Trace << "...is_d_aa: " << is_d_aa << std::endl;
			TR.Trace << "...is_l_aa: " << is_l_aa << std::endl;
		}

		ResidueTypeCOP rsd_type_cop = get_rsd_type( name3, xyz, residue_set,
			fd, options, branch_points_on_this_residue,
			resid, is_lower_terminus, is_upper_terminus,
			is_branch_point, is_branch_lower_terminus, is_d_aa, is_l_aa );

		// deprecate this assert after 2015... -- rhiju
		if ( basic::options::option[ basic::options::OptionKeys::chemical::check_rsd_type_finder ]() ) {
			ResidueTypeCOP rsd_type_cop_legacy = get_rsd_type_legacy( name3, xyz, residue_set,
				fd, options, branch_points_on_this_residue,
				resid, is_lower_terminus, is_upper_terminus,
				is_branch_point, is_branch_lower_terminus );
			if ( rsd_type_cop != rsd_type_cop_legacy  ) {
				utility_exit_with_message( "Mismatch in assigning rsd_type to PDB: found " + rsd_type_cop->name() + " vs legacy " + rsd_type_cop_legacy->name() );
			}
		}

		if ( rsd_type_cop == 0 ) {
			std::string variant;
			if ( is_lower_terminus ) {
				variant += " lower-terminal";
			} else if ( is_branch_lower_terminus ) {
				variant += " branch-lower-terminal";
			}
			if ( is_upper_terminus ) {
				variant += " upper-terminal";
			}
			if ( is_branch_point ) {
				variant += " branch-point";
			}
			utility_exit_with_message( "No match found for unrecognized residue at position " +
				boost::lexical_cast<string>(i) +
				"\nLooking for" + variant + " residue with 3-letter code: " + name3 );
		}

		ResidueType const & rsd_type( *rsd_type_cop );

		// Map the atom names.
		fill_name_map( rinfo_name_map[i], rinfo, rsd_type, options );

		debug_assert( rsd_type.natoms() >= rinfo_name_map[i].left.size() );
		core::Size missing_atoms( rsd_type.natoms() - rinfo_name_map[i].left.size() );
		if ( missing_atoms > 0 ) {
			TR.Debug << "Match: '" << rsd_type.name() << "'; missing " << missing_atoms << " coordinates" << std::endl;
		}

		debug_assert( rinfo.xyz.size() >= rinfo_name_map[i].left.size() );
		core::Size discarded_atoms( rinfo.xyz.size() - rinfo_name_map[i].left.size() );
		if ( is_lower_terminus && rinfo.xyz.count(" H  ") && ! rinfo_name_map[i].left.count(" H  ") ) {
			// Don't worry about missing BB H if Nterm
			--discarded_atoms;
		}
		if ( discarded_atoms > 0 ) {
			TR.Warning << "[ WARNING ] discarding " << discarded_atoms
				<< " atoms at position " << i << " in file " << fd.filename
				<< ". Best match rsd_type:  " << rsd_type.name() << std::endl;
		}

		// check for missing mainchain atoms:
		if ( rsd_type.is_polymer() ) {
			AtomIndices const & mainchain( rsd_type.mainchain_atoms() );
			Size const nbb( mainchain.size() );
			if ( nbb >= 3 ) {
				bool mainchain_core_present( false );
				for ( Size k=1; k<= nbb-2; ++k ) {
					std::string const & name1(rsd_type.atom_name(mainchain[k  ]));
					std::string const & name2(rsd_type.atom_name(mainchain[k+1]));
					std::string const & name3(rsd_type.atom_name(mainchain[k+2]));
					if ( !rinfo_name_map[i].right.count(name1) ||
							!rinfo_name_map[i].right.count(name2) ||
							!rinfo_name_map[i].right.count(name3) ) {
						continue;
					}
					std::string const & rinfo_name1( rinfo_name_map[i].right.find( name1 )->second );
					std::string const & rinfo_name2( rinfo_name_map[i].right.find( name2 )->second );
					std::string const & rinfo_name3( rinfo_name_map[i].right.find( name3 )->second );
					if ( xyz.count( rinfo_name1 ) && xyz.count( rinfo_name2 ) && xyz.count( rinfo_name3 ) ) {
						mainchain_core_present = true;
						break;
					}
				}
				if ( !mainchain_core_present ) {
					TR.Warning << "[ WARNING ] skipping pdb residue b/c it's missing too many mainchain atoms: " <<
						resid << ' ' << name3 << ' ' << rsd_type.name() << std::endl;
					for ( Size k=1; k<= nbb; ++k ) {
						std::string const & name(rsd_type.atom_name(mainchain[k]));
						if ( !rinfo_name_map[i].right.count(name) ||
								!xyz.count( rinfo_name_map[i].right.find(name)->second ) ) {
							// Use of unmapped name deliberate
							TR << "missing: " << name << std::endl;
						}
					}
					if ( options.exit_if_missing_heavy_atoms() == true ) {
						utility_exit_with_message("quitting due to missing heavy atoms");
					}
					continue;
				}
			}
		}

		// found a match, create the residue...
		ResidueOP new_rsd( ResidueFactory::create_residue( rsd_type ) );

		// ...and now fill in the coords
		for ( ResidueCoords::const_iterator iter=xyz.begin(), iter_end=xyz.end(); iter!= iter_end; ++iter ) {
			std::string const & rinfo_name( iter->first );
			if ( rinfo_name_map[i].left.count( rinfo_name ) ) {
				// offsetting all coordinates by a small constant prevents problems with atoms located
				// at position (0,0,0).
				// This is a bit of a dirty hack but it fixes the major problem of reading in rosetta
				// pdbs which usually start at 0,0,0. However the magnitude of this offset is so small
				// that the output pdbs should still match input pdbs. hopefully. yes. aehm.
				// RM: I'm not sure what the problem with having coordinates exactly at the origin is.
				// RM: If we do have a problem with that, it seems it should be fixed there and not here.
				// RM: (As you could imagine theoretically hitting (0,0,0) during minimization or packing.)
				double offset = 1e-250; // coordinates now double, so we can use _really_ small offset.
				std::string const & pose_name( rinfo_name_map[i].left.find( rinfo_name )->second );
				new_rsd->atom( pose_name ).xyz( iter->second + offset );
				// +1 here as we haven't added the residue to the pose yet.
				id::NamedAtomID atom_id( pose_name, pose.total_residue()+1 );
				coordinates_assigned.set( atom_id, true);
			}
			//else runtime_assert( iter->first == " H  " && rsd_type.is_terminus() ); // special casee
		}

		check_and_correct_sister_atoms( new_rsd );

		Size const old_nres( pose.total_residue() );

		if ( TR.Trace.visible() ) {
			TR.Trace << "...new residue is a polymer: " << new_rsd->type().is_polymer() << std::endl;
			if ( old_nres >= 1 ) {
				TR.Trace << "The old residue is a polymer: " << pose.residue_type(old_nres).is_polymer() << std::endl;
			}
		}

		// Add the first new residue to the pose
		if ( !old_nres ) {
			if ( TR.Trace.visible() ) {
				TR.Trace << rsd_type.name() << " " << i << " is the start of a new pose" << std::endl;
			}
			pose.append_residue_by_bond( *new_rsd );

		} else if ( ( ( is_lower_terminus && check_Ntermini_for_this_chain ) || ! same_chain_prev )
				|| /* is_branch_lower_terminus || */
				pose.residue_type( old_nres ).has_variant_type( "C_METHYLAMIDATION" ) ||
				! new_rsd->is_polymer() ||
				! pose.residue_type( old_nres ).is_polymer() ||
				! last_residue_was_recognized ) {
			// A new chain because this is a lower terminus (see logic above for designation)
			// and if we're not checking it then it's a different chain from the previous

			core::Size rootindex=1;

			// Ensure that metal ions are connected by a jump to the closest metal-binding residue that is lower in sequence.
			if ( new_rsd->is_metal() && basic::options::option[basic::options::OptionKeys::in::auto_setup_metals].user() ) {
				// If this is a metal ion and we're automatically setting up metals, search for the closest metal-binding residue
				// and make that the jump parent.  Otherwise, let the jump parent be the closest residue.
				numeric::xyzVector < core::Real > const metal_xyz = new_rsd->xyz(1); //Atom 1 is always the metal of a residue representing a metal ion.  (There's a check for this in residue_io.cc).

				core::Size closest_metalbinding_residue=0;
				core::Size metalbinding_dist_sq = 0;
				core::Size closest_residue=0;
				core::Size closest_dist_sq = 0;

				for ( core::Size jr=1, nres=pose.n_residue(); jr<=nres; ++jr ) { //Loop through all residues already added, looking for possible residues to root the metal onto.
					if ( !pose.residue(jr).is_protein() ) continue; //I'm not interested in tethering metals to non-protein residues.
					if ( !pose.residue(jr).has("CA") ) continue; //I'll be basing this on metal-alpha carbon distance, so anything without an alpha carbon won't get to be the root.

					numeric::xyzVector < core::Real > const residue_xyz = pose.residue(jr).xyz("CA");

					core::Real const current_dist_sq = residue_xyz.distance_squared(metal_xyz);

					if ( closest_residue==0 || current_dist_sq < closest_dist_sq ) {
						closest_residue = jr;
						closest_dist_sq = (core::Size)current_dist_sq;
					}
					if ( pose.residue(jr).is_metalbinding() &&
							(closest_metalbinding_residue==0 || current_dist_sq < metalbinding_dist_sq)
							) {
						closest_metalbinding_residue = jr;
						metalbinding_dist_sq = (core::Size)current_dist_sq;
					}
				} //Inner loop through all residues

				if ( closest_metalbinding_residue!=0 ) rootindex=closest_metalbinding_residue; //If we found a metal-binding residue, it's the root; otherwise, the closest residue is.
				else if ( closest_residue!=0 ) rootindex=closest_residue;

			} //If this is a metal

			if ( rootindex>1 ) { TR << rsd_type.name() << " " << i << " was added by a jump, with base residue " << rootindex << std::endl;}

			pose.append_residue_by_jump( *new_rsd, rootindex /*pose.total_residue()*/ );

		} else { // Append residue to current chain dependent on bond length.
			if ( !options.missing_dens_as_jump() ) {
				if ( TR.Trace.visible() ) {
					TR.Trace << rsd_type.name() << " " << i << " is appended to chain " << chainID << std::endl;
				}
				pose.append_residue_by_bond( *new_rsd );
			} else {
				//fpd look for missing density in the input PDB
				//fpd if there is a bondlength > 3A
				//fpd we will consider this missing density
				Residue const &last_rsd( pose.residue( old_nres ) );
				core::Real bondlength = ( last_rsd.atom( last_rsd.upper_connect_atom() ).xyz() -
					new_rsd->atom( new_rsd->lower_connect_atom() ).xyz() ).length();

				if ( bondlength > 3.0 ) {
					TR << "[ WARNING ] missing density found at residue (rosetta number) " << old_nres << std::endl;
					pose.append_residue_by_jump( *new_rsd, old_nres );

					if ( pose.residue_type(old_nres).is_protein() ) {
						if ( !pose.residue_type(old_nres).has_variant_type( UPPER_TERMINUS_VARIANT ) &&
								!pose.residue_type(old_nres).has_variant_type( UPPERTERM_TRUNC_VARIANT ) ) {
							core::pose::add_variant_type_to_pose_residue( pose, chemical::UPPERTERM_TRUNC_VARIANT, old_nres );
						}
					} else {
						if ( !pose.residue_type(old_nres).has_variant_type( UPPER_TERMINUS_VARIANT ) ) {
							core::pose::add_variant_type_to_pose_residue( pose, chemical::UPPER_TERMINUS_VARIANT, old_nres );
						}
					}

					if ( pose.residue_type(old_nres+1).is_protein() ) {
						core::pose::add_variant_type_to_pose_residue( pose, chemical::LOWERTERM_TRUNC_VARIANT, old_nres+1 );
					} else {
						core::pose::add_variant_type_to_pose_residue( pose, chemical::LOWER_TERMINUS_VARIANT, old_nres+1 );
					}
				} else {
					if ( TR.Trace.visible() ) {
						TR.Trace << rsd_type.name() << " " << i << " is appended to chain" << chainID << std::endl;
					}
					pose.append_residue_by_bond( *new_rsd );
				}
			}
		}

		// If newly added residue was a carbohydrate, set flag on conformation.
		if ( new_rsd->is_carbohydrate() ) {
			pose.conformation().contains_carbohydrate_residues( true );
		}

		pose_to_rinfo.push_back( Size(i) );
		pose_temps.push_back( rinfo.temps );

		// Update the pose-internal chain label if necessary.
		if ( ( is_lower_terminus || ! check_Ntermini_for_this_chain || is_branch_lower_terminus ) &&
				pose.total_residue() > 1 ) {
			pose.conformation().insert_chain_ending( pose.total_residue() - 1 );
		}

		last_residue_was_recognized = true;
	} // i=1,nres_pdb


	// Check termini status of newly created pose residues.
	// RM: All considered, this is a poor place to do this - we should ideally be doing this back
	// when we're originally doing the typing - though some knowledge of downstream residues is necessary.
	Size const nres( pose.total_residue() );
	for ( Size i=1; i<= nres; ++i ) {
		// Need to map pose index to rinfo index, in case we're skipping residues
		ResidueInformation const & rinfo = rinfos[pose_to_rinfo[i]];
		char chainID = rinfo.chainID;

		bool const check_Ntermini_for_this_chain = ("ALL" == chains_to_check_if_Ntermini) ?
			true : find(check_Ntermini_begin, check_Ntermini_end, chainID ) ==  check_Ntermini_end;
		bool const check_Ctermini_for_this_chain = ("ALL" == chains_to_check_if_Ctermini) ?
			true : find(check_Ctermini_begin, check_Ctermini_end, chainID ) ==  check_Ctermini_end;

		if ( !check_Ntermini_for_this_chain ) { continue; }
		if ( !check_Ctermini_for_this_chain ) { continue; }

		//Residue const & rsd( pose.residue( i ) ); // THIS WAS A BAD BUG
		bool type_changed(false);
		if ( !pose.residue_type(i).is_polymer() ) { continue; }
		if ( !pose.residue_type(i).is_lower_terminus() &&
				( i == 1 ||
				!pose.residue_type( i-1 ).is_polymer() ||
				(pose.residue_type( i-1 ).is_upper_terminus() &&
				!pose.residue_type( i ).is_branch_lower_terminus() ) ) ) {
			TR << "Adding undetected lower terminus type to residue " << i << std::endl;
			core::pose::add_lower_terminus_type_to_pose_residue( pose, i );
			type_changed = true;
		}
		if ( !pose.residue_type(i).is_upper_terminus() &&
				( i == nres ||
				!pose.residue_type(i+1).is_polymer() ||
				pose.residue_type(i+1).is_lower_terminus() /*||
				pose.residue_type(i+1).has_variant_type(BRANCH_LOWER_TERMINUS_VARIANT)*/ ) ) {
			TR << "Adding undetected upper terminus type to residue " << i << std::endl;
			core::pose::add_upper_terminus_type_to_pose_residue( pose, i );
			type_changed = true;
		}
		if ( type_changed ) {
			// add_terminus_type will copy coordinates for matching atoms - see if there's additional atoms we missed.
			for ( core::Size ii(1); ii <= pose.residue_type(i).natoms(); ++ii ) {
				std::string const & name( pose.residue_type(i).atom_name(ii) );
				id::NamedAtomID atom_id( name, i );
				// Unfortunately we're doing only exact name matches here.
				if ( ! coordinates_assigned[atom_id] && rinfo.xyz.count(name)  ) {
					pose.set_xyz( atom_id, rinfo.xyz.find(name)->second );
					coordinates_assigned.set(atom_id, true);
					rinfo_name_map[pose_to_rinfo[i]].insert( NameBimap::value_type( name, name ) );
					TR.Debug << "Setting coordinates for undetected atom " << name << " on residue " << i << std::endl;
				}
			}
		}
	}



	// now handle missing atoms
	//id::AtomID_Mask missing( false );
	Size num_heavy_missing = 0;

	core::pose::initialize_atomid_map( missing, pose ); // dimension the missing-atom mask
	if ( pose.total_residue() == 0 ) {
		// if unchecked it segfaults further down...

		// PDBInfo setup
		core::pose::PDBInfoOP pdb_info( new core::pose::PDBInfo( pose.total_residue() ) );
		for ( Size i = 1; i <= UA_res_nums.size(); ++i ) {
			pdb_info->add_unrecognized_atom( UA_res_nums[i], UA_res_names[i], UA_atom_names[i], UA_coords[i], UA_temps[i] );
		}
		// store pdb info
		pose.pdb_info( pdb_info );
		return;
	}

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		Residue const & rsd( pose.residue(i) );
		for ( Size j=1; j<= rsd.natoms(); ++j ) {
			id::AtomID atom_id( j, i );
			id::NamedAtomID named_atom_id( rsd.atom_name(j), i );
			if ( ! coordinates_assigned[named_atom_id] ) {
				missing[ atom_id ] = true;
				if ( !rsd.atom_is_hydrogen(j) ) num_heavy_missing++;
			}
		}
	}

	pose.conformation().fill_missing_atoms( missing );

	//ja save the pdb residue indices in the Pose //well, PDBInfo
	//ja pdb residue indices can be negative
	utility::vector1< int > pdb_numbering;
	//sml chain char
	utility::vector1< char > pdb_chains, insertion_codes;
	//Size const nres( pose.total_residue() );
	for ( Size i(1); i <= nres; ++i ) {
		ResidueInformation const & rinfo = rinfos[pose_to_rinfo[i]];
		pdb_numbering.push_back( rinfo.resSeq );
		pdb_chains.push_back( rinfo.chainID );
		insertion_codes.push_back( rinfo.iCode );
	}

	// PDBInfo setup
	core::pose::PDBInfoOP pdb_info( new core::pose::PDBInfo( pose.total_residue() ) );

	// set pdb-wide information
	pdb_info->name( fd.filename );
	if ( fd.modeltag=="" ) {
		pdb_info->modeltag( fd.filename );
	} else {
		pdb_info->modeltag( fd.modeltag );
	}

	if ( options.preserve_header() == true ) {
		pdb_info->remarks( *fd.remarks );
		pdb_info->header_information( fd.header_information() );
	}

	// set residue level pdb information
	pdb_info->set_numbering( pdb_numbering );
	pdb_info->set_chains( pdb_chains );
	pdb_info->set_icodes( insertion_codes );
	if ( options.preserve_crystinfo() ) {
		pdb_info->set_crystinfo( fd.crystinfo );
	}


	// most DNA structures lack 5' phosphate groups. 5' phosphates must be built to serve as part of the backbone for
	// atom/fold tree purposes. Here they are made virtual so as not to affect physical calculations.
	for ( uint seqpos(1), nres( pose.total_residue() ); seqpos <= nres; ++seqpos ) {
		Residue const & rsd( pose.residue( seqpos ) );
		if ( ! rsd.type().is_DNA() ) continue;
		for ( uint atomi(1), natoms( rsd.natoms() ); atomi <= natoms; ++atomi ) {
			id::AtomID const id( atomi, seqpos );
			if ( missing[ id ] && rsd.atom_name(atomi) == " P  " ) {
				TR << "Virtualizing missing phosphate that was built in at seqpos " << seqpos << std::endl;
				core::pose::add_variant_type_to_pose_residue( pose, chemical::VIRTUAL_DNA_PHOSPHATE, seqpos );
				break;
			}
		}
	}


	// Look for and create any remaining non-mainchain (Edge::CHEMICAL) bonds based on a specified radius from any
	// unsatisfied residue connections.  This is used for such things as branched polymers, ubiquitination, or covalent
	// intermediates.  Note: The fold tree will remain with a jump between each such bond until import_pose::
	// set_reasonable_fold_tree() is called later, which actually adds the CHEMICAL edges to fold tree; this method
	// simply makes the bonds.
	pose.conformation().detect_bonds();

	//mjo TODO: this can try to access pose->pdb_info() which is not yet
	//initialized. Moving it after the pose->pdb_info has been
	//initialized causes integration test changes
	core::pose::initialize_disulfide_bonds(pose);

	//kdrew: if detect_oops flag is set, initialize oops
	// This option should probably be moved to FileDataOptions. ~Labonte
	if ( basic::options::option[ basic::options::OptionKeys::in::detect_oops ].user() ) {
		core::pose::ncbb::initialize_oops(pose);
	}

	if ( pose.n_residue()>1 ) { // 1 residue fragments for ligand design.
		pose.conformation().detect_pseudobonds();
	}

	// ensure enough space for atom level pdb information
	pdb_info->resize_atom_records( pose );

	// add unrecognized atoms to PDBInfo
	for ( Size i = 1; i <= UA_res_nums.size(); ++i ) {
		pdb_info->add_unrecognized_atom( UA_res_nums[i], UA_res_names[i], UA_atom_names[i], UA_coords[i], UA_temps[i] );
	}

	// add temps to PDBInfo
	for ( core::Size ir = 1; ir <= pose.total_residue(); ir++ ) {
		// fill in b-factor from pdb file
		ResidueTemps & res_temps( rinfos[pose_to_rinfo[ir]].temps );
		NameBimap const & namemap( rinfo_name_map[pose_to_rinfo[ir]] );
		for ( ResidueTemps::const_iterator iter=res_temps.begin(); iter != res_temps.end(); ++iter ) {
			//namemap should only include atoms which have a presence in both rinfo and pose
			if ( namemap.left.count(iter->first) ) {
				// printf("setting temp: res %d atom %s temp %f\n",ir,iter->first.c_str(),iter->second);
				std::string const & pose_atom_name( namemap.left.find(iter->first)->second );
				if ( pose.residue(ir).type().has( pose_atom_name ) ) { // There are issues with terminus patching which means atoms can sometimes disappear
					core::Size ia = pose.residue(ir).type().atom_index( pose_atom_name );
					pdb_info->temperature( ir, ia, iter->second );
				}
			} else {
				if ( (iter->first)[0] == 'H' || ((iter->first)[0] == ' ' && (iter->first)[1] == 'H') ) {
					;// don't warn if H
				} else {
					TR << "[ WARNING ] can't find atom for res " << ir << " atom " << iter->first << " (trying to set temp)" << std::endl;
				}
			}
		}
	}

	// mark PDBInfo as ok and store in Pose
	pdb_info->obsolete( false );
	pose.pdb_info( pdb_info );

	// Can't do this until we have a PDBInfo, duh!
	// add contraints based on LINK records if desired
	if ( basic::options::option[ basic::options::OptionKeys::in::constraints_from_link_records ].value() ) {
		core::pose::get_constraints_from_link_records( pose, fd );
	}

	//fpd fix bfactors of missing atoms using neighbors
	//fpd set hydrogen Bfactors as 1.2x attached atom
	if ( options.preserve_crystinfo() ) {
		core::scoring::cryst::fix_bfactorsMissing( pose );
		core::scoring::cryst::fix_bfactorsH( pose );
	}
	if ( basic::options::option[ basic::options::OptionKeys::out::file::pdb_comments]() ) {
		utility::io::izstream data(fd.filename);

		std::string line;
		while ( getline( data, line ) ) {
			if ( line != "##Begin comments##" ) {
				continue;
			}
			getline( data, line );
			while ( line != "##End comments##" ) {
				//TR<<"Testing read comments! :"<<line<<std::endl;
				utility::vector1<std::string> comment_line(utility::string_split(line,' '));
				if ( comment_line.size()<2 ) {
					getline( data, line );
					continue;
				}
				core::pose::add_comment(pose,comment_line[1],comment_line[2]);
				getline( data, line );
			}
		}
	}
}


/// @details The input rsd_type_list are all the residue types that have
///the same 3 letter code as pdb_name. Return true if the list is
///non-empty and false otherwise.  If no residue types match, then
///either exit, ignore or remember the residue based on the following
///options in the option system:
///
/// -in:ignore_waters
/// -in:ignore_unrecognized_res
/// -in:remember_unrecognized_waters
/// -in:remember_unrecognized_res
bool
is_residue_type_recognized(
	Size const pdb_residue_index,
	std::string const & pdb_name,
	core::chemical::ResidueTypeCOPs const & rsd_type_list,
	std::map< std::string, Vector > const & xyz,
	std::map< std::string, double > const & rtemp,
	utility::vector1<Size> & UA_res_nums,
	utility::vector1<std::string> & UA_res_names,
	utility::vector1<std::string> & UA_atom_names,
	utility::vector1<numeric::xyzVector<Real> > & UA_coords,
	utility::vector1<core::Real> & UA_temps ) {

	FileDataOptions options;
	return is_residue_type_recognized( pdb_residue_index, pdb_name, rsd_type_list, xyz, rtemp,
		UA_res_nums, UA_res_names, UA_atom_names, UA_coords, UA_temps, options );
}

/// @brief Tell user about ancient Rosetta choice to ignore HOH waters if -ignore_unrecognized_res is set.
void
output_ignore_water_warning_once( FileDataOptions const & options ) {
	static bool outputted_warning ( false );
	if ( outputted_warning ) return;
	if ( options.ignore_unrecognized_res() && options.ignore_waters() ) {
		TR << TR.Red << "For backwards compatibility, setting -ignore_unrecognized_res leads ALSO to -ignore_waters. You can set -ignore_waters false to get waters." << TR.Reset << std::endl;
		outputted_warning = true;
	}
}

/////////////////////////////////////////////////////////
bool
is_residue_type_recognized(
	Size const pdb_residue_index,
	std::string const & pdb_name,
	core::chemical::ResidueTypeSet const & residue_set,
	std::map< std::string, Vector > const & xyz,
	std::map< std::string, double > const & rtemp,
	utility::vector1<Size> & UA_res_nums,
	utility::vector1<std::string> & UA_res_names,
	utility::vector1<std::string> & UA_atom_names,
	utility::vector1<numeric::xyzVector<Real> > & UA_coords,
	utility::vector1<core::Real> & UA_temps,
	FileDataOptions const & options){
	using namespace core::chemical;

	// this residue list is only used to see if there are any residue_types with name3 at all:
	ResidueTypeCOPs rsd_type_list;
	ResidueTypeCOP  rsd_type = ResidueTypeFinder( residue_set ).name3( pdb_name ).get_representative_type();
	if ( rsd_type != 0 ) rsd_type_list.push_back( rsd_type );
	return is_residue_type_recognized( pdb_residue_index, pdb_name,
		rsd_type_list,
		xyz, rtemp, UA_res_nums, UA_res_names,
		UA_atom_names, UA_coords, UA_temps,
		options );
}

/// @details The input rsd_type_list are all the residue types that have
///the same 3 letter code as pdb_name. Return true if the list is
///non-empty and false otherwise.  If no residue types match, then
///either exit, ignore or remember the residue based on the following
///options in a FileDataOptions instance:
///
/// -ignore_waters
/// -ignore_unrecognized_res
/// -remember_unrecognized_waters
/// -remember_unrecognized_res
bool
is_residue_type_recognized(
	Size const pdb_residue_index,
	std::string const & pdb_name,
	core::chemical::ResidueTypeCOPs const & rsd_type_list,
	std::map< std::string, Vector > const & xyz,
	std::map< std::string, double > const & rtemp,
	utility::vector1<Size> & UA_res_nums,
	utility::vector1<std::string> & UA_res_names,
	utility::vector1<std::string> & UA_atom_names,
	utility::vector1<numeric::xyzVector<Real> > & UA_coords,
	utility::vector1<core::Real> & UA_temps,
	FileDataOptions const & options){

	bool const is_HOH_to_ignore ( pdb_name == "HOH" && options.ignore_waters() );

	if ( !rsd_type_list.empty() && !is_HOH_to_ignore ) {
		return true;
	}

	using namespace basic::options;
	if ( !(options.ignore_unrecognized_res() ||
			options.remember_unrecognized_res() ||
			is_HOH_to_ignore ) ) {
		// We should fail fast on unrecognized input rather than produce bad results!
		utility_exit_with_message(" unrecognized residue " + pdb_name );
	}

	if ( !options.remember_unrecognized_water() ) {
		// don't bother with water
		if ( pdb_name == "HOH" ) {
			return false;
		}
	}

	if ( options.remember_unrecognized_res() ) {
		for ( std::map<std::string, Vector>::const_iterator iter=xyz.begin(), iter_end=xyz.end(); iter!= iter_end; ++iter ) {
			if ( UA_res_nums.size() > 5000 ) {
				utility_exit_with_message("can't handle more than 5000 atoms worth of unknown residues\n");
			}
			TR << "remember unrecognized atom " << pdb_residue_index << " " << pdb_name << " " << stripped_whitespace(iter->first)
				<< " temp " << rtemp.find(iter->first)->second << std::endl;
			UA_res_nums.push_back( pdb_residue_index );
			UA_res_names.push_back( pdb_name );
			UA_atom_names.push_back( stripped_whitespace(iter->first) );
			UA_coords.push_back( iter->second );
			UA_temps.push_back( rtemp.find(iter->first)->second );
		}
	}

	if ( is_HOH_to_ignore ) output_ignore_water_warning_once( options );
	return false;
}

void
pose_from_pose(
	pose::Pose & new_pose,
	pose::Pose const & old_pose,
	utility::vector1< core::Size > const & residue_indices
){
	FileDataOptions options;
	pose_from_pose( new_pose, old_pose, residue_indices, options );
}


void
pose_from_pose(
	pose::Pose & new_pose,
	pose::Pose const & old_pose,
	utility::vector1< core::Size > const & residue_indices,
	FileDataOptions const & options
){
	using namespace chemical;
	ResidueTypeSetCOP residue_set(
		ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )
	);
	pose_from_pose( new_pose, old_pose, *residue_set,  residue_indices, options);
}


void
pose_from_pose(
	pose::Pose & new_pose,
	pose::Pose const & old_pose,
	chemical::ResidueTypeSet const & residue_set,
	utility::vector1< core::Size > const & residue_indices
){
	FileDataOptions options;
	pose_from_pose( new_pose, old_pose, residue_set, residue_indices, options );
}


/// Creates a subpose from a pose, to include only certain
/// residues, using FileData::init_from_pose() to construct the
/// pose, and build_pose_as_is1() to construct the pose
/// with the given options.
void
pose_from_pose(
	pose::Pose & new_pose,
	pose::Pose const & old_pose,
	chemical::ResidueTypeSet const & residue_set,
	utility::vector1< core::Size > const & residue_indices,
	FileDataOptions const & options
){
	FileData fd;
	fd.init_from_pose( old_pose, residue_indices );
	id::AtomID_Mask missing( false );
	build_pose_as_is1( fd, new_pose, residue_set, missing, options );
}

///////////////////////////////////////////////////////////////////////
// @brief currently does fixups of RNA/DNA. Could also add other 'special case' fixups here.
void
fixup_rinfo_based_on_residue_type_set(  utility::vector1< ResidueInformation > & rinfos,
	chemical::ResidueTypeSet const & /* residue_set */ ){
	// As of Aug 2015, deprecating RNA-only residue type set. That means that
	//  PDBs that are mostly RNA but have deoxy residues will be read in as such, and should be handled
	//  appropriately by the applications, not here inside file_data.cc. -- rhiju
	//bool const force_RNA = false;
	convert_nucleic_acid_residue_info_to_standard( rinfos /*, force_RNA*/ );
}


///////////////////////////////////////////////////////////////////////
// @brief Use ResidueTypeFinder to efficiently figure out best match
//    residue_type to these PDB atom_names, name3, etc.
chemical::ResidueTypeCOP
get_rsd_type(
	std::string const & name3,
	std::map< std::string, Vector > const & xyz,
	chemical::ResidueTypeSet const & residue_set,
	FileData & fd,
	FileDataOptions const & options,
	utility::vector1< std::string > const &  branch_points_on_this_residue,
	std::string const & resid,
	bool const is_lower_terminus,
	bool const is_upper_terminus,
	bool const is_branch_point,
	bool const is_branch_lower_terminus,
	bool const is_d_aa,
	bool const is_l_aa )
{
	// you can be neither but not both
	debug_assert( ! ( is_d_aa && is_l_aa ) );

	typedef std::map< std::string, Vector > ResidueCoords;
	using namespace core::chemical;
	using utility::tools::make_vector1;
	using utility::vector1;

	vector1< vector1< VariantType > > required_variants_in_sets;
	vector1< ResidueProperty > properties, disallow_properties;
	vector1< VariantType >     disallow_variants;  // are variants different from properties?
	std::string                residue_base_name( "" ); // carbohydrates
	vector1< std::string >     patch_names;
	bool  disallow_carboxyl_conjugation_at_glu_asp( false );

	ResidueTypeCOP rsd_type = ResidueTypeFinder( residue_set ).name3( name3 ).get_representative_type();
	if ( rsd_type->is_polymer() ) {
		if ( is_lower_terminus ) {
			required_variants_in_sets.push_back( make_vector1( LOWER_TERMINUS_VARIANT, LOWERTERM_TRUNC_VARIANT ) );
		} else {
			disallow_variants.push_back( LOWER_TERMINUS_VARIANT );
			disallow_variants.push_back( LOWERTERM_TRUNC_VARIANT );
		}
		if ( is_upper_terminus ) {
			required_variants_in_sets.push_back( make_vector1( UPPER_TERMINUS_VARIANT, UPPERTERM_TRUNC_VARIANT ) );
		} else {
			disallow_variants.push_back( UPPER_TERMINUS_VARIANT );
			disallow_variants.push_back( UPPERTERM_TRUNC_VARIANT );
		}
		// note that carbohydrate branch points will be covered by patch_names below.
		if ( !rsd_type->is_carbohydrate() ) {
			if ( is_branch_point  )  properties.push_back( BRANCH_POINT ); else disallow_properties.push_back( BRANCH_POINT );
		}
		if ( is_branch_lower_terminus )  properties.push_back( BRANCH_LOWER_TERMINUS ); else disallow_properties.push_back( BRANCH_LOWER_TERMINUS );
		if ( is_d_aa )  properties.push_back( D_AA );
		if ( is_l_aa )  properties.push_back( L_AA );
	}
	if ( rsd_type->aa() == aa_cys && name3 != "CYD" ) disallow_variants.push_back( DISULFIDE );
	if ( !is_branch_lower_terminus ) disallow_carboxyl_conjugation_at_glu_asp = true;

	if ( !options.keep_input_protonation_state() ) {
		disallow_variants.push_back( PROTONATED );
		disallow_variants.push_back( DEPROTONATED );
	}
	if ( rsd_type->is_carbohydrate() ) {
		residue_base_name = fd.residue_type_base_names[ resid ];
		// The below assumes that ResidueTypes with fewer patches are selected 1st, that is, that an
		// :->2-branch ResidueType will be checked as a possible match before an :->2-branch:->6-branch
		// ResidueType.  If this were not the case, Rosetta could misassign an :->2-branch:->6-branch
		// ResidueType to a residue that actually only has a single branch at the 2 or 6 position.
		for ( core::uint k( 1 ); k <= branch_points_on_this_residue.size(); ++k ) {
			char const & branch_point = branch_points_on_this_residue[ k ][ 2 ];  // 3rd column (index 2) is the atom number.
			patch_names.push_back( "->" + string( 1, branch_point ) + ")-branch" );
		}
	}

	utility::vector1< std::string > xyz_atom_names;
	for ( ResidueCoords::const_iterator iter=xyz.begin(), iter_end=xyz.end(); iter!= iter_end; ++iter ) {
		std::string xyz_name = iter->first;
		xyz_atom_names.push_back( xyz_name );
	}

	// following 'chaining' looks a lot like chemical::ResidueSelector -- use that instead? Rename to chemical::ResidueTypeSelector.
	rsd_type = ResidueTypeFinder( residue_set ).name3( name3 ).residue_base_name( residue_base_name ).disallow_variants( disallow_variants ).variants_in_sets( required_variants_in_sets ).properties( properties ).disallow_properties( disallow_properties ).patch_names( patch_names ).ignore_atom_named_H( is_lower_terminus ).disallow_carboxyl_conjugation_at_glu_asp( disallow_carboxyl_conjugation_at_glu_asp ).get_best_match_residue_type_for_atom_names( xyz_atom_names );

	return rsd_type;
}

///////////////////////////////////////////////////////////////////////
/// DEPRECATE after 2015. Replaced with more efficient get_rsd_type(); see above.
chemical::ResidueTypeCOP
get_rsd_type_legacy(
	std::string const & name3,
	std::map< std::string, Vector > const & xyz,
	chemical::ResidueTypeSet const & residue_set,
	FileData & fd,
	FileDataOptions const & options,
	utility::vector1< std::string > const &  branch_points_on_this_residue,
	std::string const & resid,
	bool const is_lower_terminus,
	bool const is_upper_terminus,
	bool const is_branch_point,
	bool const is_branch_lower_terminus )
{

	typedef std::map< std::string, Vector > ResidueCoords;
	using namespace core::chemical;

	// This is a problem...
	ResidueTypeCOPs const & rsd_type_list( residue_set.name3_map_DO_NOT_USE( name3 ) );

	// look for best match:
	// rsd_type should have all the atoms present in xyz
	// try to minimize atoms missing from xyz
	Size best_index(0), best_rsd_missing( 99999 ), best_xyz_missing( 99999 );

	for ( Size j=1; j<= rsd_type_list.size(); ++j ) {
		ResidueType const & rsd_type( *(rsd_type_list[j]) );
		bool const is_polymer( rsd_type.is_polymer() ); // need an example residue type, though this will
		// remain fixed for all residue_types with the same name3

		// only take the desired variants`
		bool lower_term_type = rsd_type.has_variant_type( LOWER_TERMINUS_VARIANT ) ||
			rsd_type.has_variant_type( LOWERTERM_TRUNC_VARIANT );
		bool upper_term_type = rsd_type.has_variant_type( UPPER_TERMINUS_VARIANT ) ||
			rsd_type.has_variant_type( UPPERTERM_TRUNC_VARIANT );
		if ( is_polymer && (
				(is_lower_terminus != lower_term_type ) || (is_upper_terminus != upper_term_type ) ) ) {
			if ( TR.Trace.visible() ) {
				TR.Trace << "Discarding '" << rsd_type.name() << "' ResidueType" << std::endl;
				TR.Trace << "because of the terminus state" << std::endl;
				TR.Trace << "PDB has lower " << is_lower_terminus << " and type has " << lower_term_type << std::endl;
				TR.Trace << "PDB has upper " << is_upper_terminus << " and type has " << upper_term_type << std::endl;
			}
			continue;
		}
		if ( is_polymer && ( is_branch_point != rsd_type.is_branch_point() ) ) {
			if ( TR.Trace.visible() ) {
				TR.Trace << "Discarding '" << rsd_type.name() << "' ResidueType" << std::endl;
				TR.Trace << "because of the branch state" << std::endl;
				TR.Trace << "PDB has branch " << is_branch_point << " and type has " << rsd_type.is_branch_point() << std::endl;
			}
			continue;
		}
		if ( is_polymer && ( is_branch_lower_terminus != rsd_type.is_branch_lower_terminus() ) ) {
			if ( TR.Trace.visible() ) {
				TR.Trace << "Discarding '" << rsd_type.name() << "' ResidueType" << std::endl;
				TR.Trace << "because of the branch lower terminus state" << std::endl;
				TR.Trace << "PDB has branch " << is_branch_lower_terminus << " and type has " << rsd_type.is_branch_lower_terminus() << std::endl;
			}
			continue;
		}
		// Okay, this logic is NOT OBVIOUS, so I am going to add my explanation.
		// We DEFINITELY want to assign a disulfide type from the start if the PDB just up and says CYD. That's great!
		// But if we DO NOT see CYD, we do not want to assign a disulfide type. We want disulfide connections
		// to be inferred later on, in conformation's detect_disulfides as called in the pose-building process.
		// Commenting out this logic causes anything that lacks a HG (i.e. crystal structures) to be assigned as the
		// disulfide type instead of the CYS type--which MIGHT be right, but might be wrong and leads to wasteful
		// disulfide reversion.
		if ( rsd_type.aa() == aa_cys && rsd_type.has_variant_type( DISULFIDE ) && name3 != "CYD" ) {
			if ( TR.Trace.visible() ) {
				TR.Trace << "Discarding '" << rsd_type.name() << "' ResidueType" << std::endl;
				TR.Trace << "because of the disulfide state" << std::endl;
			}
			continue;
		}
		if ( ( rsd_type.aa() == aa_glu || rsd_type.aa() == aa_asp ) && rsd_type.has_variant_type( BRANCH_LOWER_TERMINUS_VARIANT ) && !is_branch_lower_terminus ) {
			// Don't assign sidechain carboxyl conjugation upon read-in unless mentioned in a LINK record!
			if ( TR.Trace.visible() ) {
				TR.Trace << "Discarding '" << rsd_type.name() << "' ResidueType" << std::endl;
				TR.Trace << "because of the sidechain conjugation state" << std::endl;
			}
			continue;
		}
		if ( ! options.keep_input_protonation_state() &&
				( rsd_type.has_variant_type( PROTONATED ) || rsd_type.has_variant_type( DEPROTONATED ) ) ) {
			if ( TR.Trace.visible() ) {
				TR.Trace << "Discarding '" << rsd_type.name() << "' ResidueType" << std::endl;
				TR.Trace << "because of the protonation state" << std::endl;
			}
			continue;
		}

		// special checks to ensure selecting the proper carbohydrate ResidueType
		if ( rsd_type.is_carbohydrate() &&
				residue_type_base_name( rsd_type ) != fd.residue_type_base_names[ resid ] ) {
			if ( TR.Trace.visible() ) {
				TR.Trace << "Discarding '" << rsd_type.name() << "' ResidueType" << std::endl;
				TR.Trace << "because the residue is not a carbohydrate" << std::endl;
			}
			continue;
		}
		if ( rsd_type.is_carbohydrate() && rsd_type.is_branch_point() ) {
			// The below assumes that ResidueTypes with fewer patches are selected 1st, that is, that an
			// :->2-branch ResidueType will be checked as a possible match before an :->2-branch:->6-branch
			// ResidueType.  If this were not the case, Rosetta could misassign an :->2-branch:->6-branch
			// ResidueType to a residue that actually only has a single branch at the 2 or 6 position.
			//char branch_point;
			bool branch_point_is_missing( false );
			Size const n_branch_points( branch_points_on_this_residue.size() );
			for ( core::uint k( 1 ); k <= n_branch_points; ++k ) {
				char branch_point = branch_points_on_this_residue[ k ][ 2 ];  // 3rd column (index 2) is the atom number.
				if ( TR.Debug.visible() ) {
					TR.Debug << "Checking '" << rsd_type.name() <<
						"' for branch at position " << branch_point << std::endl;
				}
				if ( residue_type_all_patches_name( rsd_type ).find( string( 1, branch_point ) + ")-branch" ) ==
						string::npos ) {
					branch_point_is_missing = true;
					break;
				}
			}
			if ( branch_point_is_missing ) {
				if ( TR.Trace.visible() ) {
					TR.Trace << "Discarding '" << rsd_type.name() << "' ResidueType" << std::endl;
					TR.Trace << "because of a missing branch point" << std::endl;
				}
				continue;
			}
		}

		if ( TR.Debug.visible() ) {
			TR.Debug << "Trying '" << rsd_type.name() << "' ResidueType" << std::endl;
		}

		Size rsd_missing( 0 ), xyz_missing( 0 );

		for ( Size k=1; k<= rsd_type.natoms(); ++k ) {
			if ( xyz.count( rsd_type.atom_name(k) ) == 0 ) { ++xyz_missing; }
		}

		for ( ResidueCoords::const_iterator iter=xyz.begin(), iter_end=xyz.end(); iter!= iter_end; ++iter ) {
			if ( !rsd_type.has( stripped_whitespace(iter->first) ) &&
					!( iter->first == " H  " && is_lower_terminus ) ) { // don't worry about missing BB H if Nterm
				++rsd_missing;
			}
		}

		if ( TR.Debug.visible() ) {
			//   TR.Debug << "Trying '" << rsd_type.name() << "' ResidueType   xyz_missing " << xyz_missing << " rsd_missing " << rsd_missing <<  std::endl;
		}

		if ( ( rsd_missing < best_rsd_missing ) ||
				( rsd_missing == best_rsd_missing && xyz_missing < best_xyz_missing ) ) {
			best_rsd_missing = rsd_missing;
			best_xyz_missing = xyz_missing;
			best_index = j;
		}
	} // j=1,rsd_type_list.size()

	ResidueTypeCOP rsd_type = rsd_type_list[ best_index ];

	if ( TR.Trace.visible() ) {
		TR.Trace << "Naive match of " << rsd_type->name() << " with " << best_rsd_missing << " missing and "
			<< best_xyz_missing << " discarded atoms." << std::endl;
	}

	if ( best_index == 0 ) return 0;

	return rsd_type;
}


} // namespace pdb
} // namespace io
} // namespace core
