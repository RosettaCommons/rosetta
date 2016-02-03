// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file core/io/pdb/pdb_writer.cc
/// @brief Function(s) for pdb writing
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com), XRW 2016 Team
/// @author Sergey Lyskov (Sergey.Lyskov@jhu.edu)



// Unit headers
#include <core/io/pdb/pdb_writer.hh>

// Package headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>

#include <core/io/Remarks.hh>
#include <core/io/HeaderInformation.hh>
#include <core/io/CrystInfo.hh>
#include <core/io/pdb/Field.hh>
#include <core/io/pdb/RecordCollection.hh>
#include <core/io/pose_to_sfr/PoseToStructFileRepConverter.hh>
#include <core/io/StructFileRepOptions.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/io/util.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/tools/make_map.hh>
#include <utility/string_util.hh>
#include <utility/io/ozstream.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// Basic headers
#include <basic/Tracer.hh>

// External headers
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

// Basic headers (JAB - remove as many of these as possible)
#include <basic/options/option.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/Tracer.hh>


// C++ headers
#include <cstdlib>
#include <cstdio>
#include <algorithm>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "core.io.pdb.pdb_writer" );

namespace core {
namespace io {
namespace pdb {

using namespace ObjexxFCL::format;
using namespace pose_to_sfr;

using core::Size;
using core::SSize;
using utility::pad_left;
using utility::pad_right;
using utility::fmt_real;

using basic::T;
using basic::Error;
using basic::Warning;

/// @brief Writes a pose to a given stream in PDB file format, optionally
/// appending a given string and optionally extracting scores from the pose.
/// @details This came out of the 2016 Chemical XRW.  It's an attempt to preserve
/// some stuff that jd2 was doing before, while centralizing all PDB generation in
/// one place.
/// @param[in] pose The pose to turn into a PDB.
/// @param[in] extra_data Additional data to append to the PDB file data.
/// @param[in] add_score_data Grab additional score data from the pose?
/// @param[in] add_extra_score_data Grab still more score data from the pose?
/// @param[out] out The output stream that the PDB file data will be written to.
/// @param[in] filename (Optional)  String for the filename.  Will be included in the score data table if provided.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
dump_pdb(
	core::pose::Pose const & pose,
	std::string const &extra_data,
	bool const add_score_data,
	bool const add_extra_score_data,
	utility::io::ozstream & out,
	std::string const &filename
) {
	io::pose_to_sfr::PoseToStructFileRepConverter pose_to_sfr;
	pose_to_sfr.init_from_pose( pose );
	if ( !extra_data.empty() ) pose_to_sfr.sfr()->append_to_additional_string_output( extra_data );
	if ( add_score_data ) pose_to_sfr.sfr()->append_to_additional_string_output( core::io::extract_scores(pose, filename) );
	if ( add_extra_score_data ) pose_to_sfr.sfr()->append_to_additional_string_output( core::io::extract_extra_scores(pose) );
	out << create_pdb_contents_from_sfr( *( pose_to_sfr.sfr() ) );
}

/// @brief Writes a pose to a given string in PDB file format, optionally
/// appending a given string and optionally extracting scores from the pose.
/// @details This came out of the 2016 Chemical XRW.  It's an attempt to preserve
/// some stuff that jd2 was doing before, while centralizing all PDB generation in
/// one place.
/// @param[in] pose The pose to turn into a PDB.
/// @param[in] extra_data Additional data to append to the PDB file data.
/// @param[in] add_score_data Grab additional score data from the pose?
/// @param[in] add_extra_score_data Grab still more score data from the pose?
/// @param[out] out The output string that the PDB file data will be written to.
/// @param[in] filename (Optional)  String for the filename.  Will be included in the score data table if provided.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
dump_pdb(
	core::pose::Pose const & pose,
	std::string const &extra_data,
	bool const add_score_data,
	bool const add_extra_score_data,
	std::string & out,
	std::string const &filename
) {
	io::pose_to_sfr::PoseToStructFileRepConverter pose_to_sfr;
	pose_to_sfr.init_from_pose( pose ); //Yeah, there's a wee bit of code duplication here.
	if ( !extra_data.empty() ) pose_to_sfr.sfr()->append_to_additional_string_output( extra_data ); //And here.
	if ( add_score_data ) pose_to_sfr.sfr()->append_to_additional_string_output( core::io::extract_scores(pose, filename) ); //And here.
	if ( add_extra_score_data ) pose_to_sfr.sfr()->append_to_additional_string_output( core::io::extract_extra_scores(pose) ); //And here.  Four lines.  It's an XRW.
	out = create_pdb_contents_from_sfr( *( pose_to_sfr.sfr() ) );
}

/// @details Convert given Pose object in to PDB format and send it to the given stream.
void
dump_pdb(
	core::pose::Pose const & pose,
	std::ostream & out,
	String const & /* tag */,
	bool write_fold_tree
)
{
	String data;
	StructFileRepOptions options = StructFileRepOptions();
	options.set_fold_tree_io( write_fold_tree );

	io::pose_to_sfr::PoseToStructFileRepConverter pu = io::pose_to_sfr::PoseToStructFileRepConverter();
	pu.init_from_pose(pose, options);

	data = create_pdb_contents_from_sfr(*pu.sfr());
	out.write( data.c_str(), data.size() );
}


/// @details Convert given Pose object into PDB format and save it to 'file_name' file.
/// @return true if operation was completed without error, false otherwise.
bool
dump_pdb(
	core::pose::Pose const & pose,
	String const & file_name,
	String const & tag,
	bool write_fold_tree )
{
	utility::io::ozstream file(file_name.c_str(), std::ios::out | std::ios::binary);
	if ( !file ) {
		Error() << "StructFileRep::dump_pdb: Unable to open file:" << file_name << " for writing!!!" << std::endl;
		return false;
	}
	dump_pdb(pose, file, tag, write_fold_tree);
	file.close();

	return true;
}




/// @details Convert given Pose object into PDB format and send it to the given stream.
/// only the residues corresponding to indices in 'residue_indices' will be output
void
dump_pdb(
	core::pose::Pose const & pose,
	std::ostream & out,
	utility::vector1< core::Size > const & residue_indices,
	String const & /* tag */
)
{
	io::pose_to_sfr::PoseToStructFileRepConverter pu;
	String data;
	pu.init_from_pose( pose, residue_indices );

	data = create_pdb_contents_from_sfr(*pu.sfr());
	out.write( data.c_str(), data.size() );
}



/// @brief create a full pdb as a string given a StructFileRep object
std::string
create_pdb_contents_from_sfr( StructFileRep const & sfr )
{
	utility::vector1< Record > records( create_records_from_sfr( sfr ) );
	std::string pdb_contents;
	pdb_contents.reserve( 81 * records.size() );
	for ( Size i = 1; i <= records.size(); ++i ) {
		pdb_contents += create_pdb_line_from_record( records[ i ] ) + '\n';
	}
	return pdb_contents;
}

/// @details create PDB string from Record data.
/// Thank you, Frank! ~Labonte
std::string
create_pdb_line_from_record( Record const & record )
{
	std::string line;
	//std::cout << "RECORD: size: " << record.size() << ": ";
	line.resize( 80, ' ' );  // Unless it's a non-standard record, it will have a length of 80.
	//std::cout << "How many fields? " << record.size() << std::endl;
	for ( Record::const_iterator field = record.begin(); field != record.end(); ++field ) {
		std::string value( field->second.value );
		if ( field->second.start == 0 ) continue; // Fields that are not properly defined (e.g. UNKNOW's "type" field) should not be written.
		if ( value.length() == 0 ) continue; // Cannot output a value with zero length
		//std::cout << "value: \"" << value << "\"" << std::endl;
		if ( field->second.end == 0 ) {  // If this field has no size limit.
			Size len = value.length();
			if ( field->second.start + len > 80 ) {
				line.resize( field->second.start + len - 1 );
			}
			line.replace( field->second.start - 1, len, value );
		} else {
			value.resize( field->second.end - field->second.start + 1, ' ' );
			line.replace( field->second.start - 1, field->second.end - field->second.start + 1, value );
		}
	}
	return line;
}


/// @details Create vector of Record from given StructFileRep object.  Used in PDB writing support.
std::vector<Record>
create_records_from_sfr( StructFileRep const & sfr )
{
	std::vector<Record> VR;

	sfr.header()->fill_records(VR);

	Record R = RecordCollection::record_from_record_type( "REMARK" );
	for ( Size i=0; i<sfr.remarks()->size(); ++i ) {
		RemarkInfo const & ri( sfr.remarks()->at(i) );
		R["type"].value = "REMARK";
		R["remarkNum"].value = pad_left( ri.num, 3 ); //("%3d", ri.num);
		R["value"].value = ri.value;
		VR.push_back(R);
	}

	R = RecordCollection::record_from_record_type( "HETNAM" );
	Size n_het_names = sfr.heterogen_names().size();
	for ( uint i = 1; i <= n_het_names; ++i ) {
		R["type"].value = "HETNAM";
		R["continuation"].value = "  ";  // TODO: Wrap long text fields.
		R["hetID"].value = sfr.heterogen_names()[i].first;
		R["text"].value = sfr.heterogen_names()[i].second;
		VR.push_back(R);
	}

	R = RecordCollection::record_from_record_type( "LINK  " );
	std::map<std::string, utility::vector1<LinkInformation> >::const_iterator last_branch_point = sfr.link_map().end();
	for ( std::map<std::string, utility::vector1<LinkInformation> >::const_iterator branch_point = sfr.link_map().begin();
			branch_point != last_branch_point; ++branch_point ) {
		utility::vector1<LinkInformation> links = branch_point->second;
		Size n_links = links.size();
		for ( uint i = 1; i <= n_links; ++i ) {
			LinkInformation link = links[i];
			R["type"].value = "LINK  ";
			R["name1"].value = link.name1;
			R["resName1"].value = link.resName1;
			R["chainID1"].value = std::string( 1, link.chainID1 );
			R["resSeq1"].value = pad_left( link.resSeq1, 4 ); //("%4d", link.resSeq1);
			R["iCode1"].value = std::string( 1, link.iCode1 );
			R["name2"].value = link.name2;
			R["resName2"].value = link.resName2;
			R["chainID2"].value = std::string( 1, link.chainID2 );
			R["resSeq2"].value = pad_left( link.resSeq2, 4 ); //("%4d", link.resSeq2);
			R["iCode2"].value = std::string( 1, link.iCode2 );
			R["length"].value = fmt_real( link.length, 2 , 2 ); //("%5.2f", link.length);
			VR.push_back(R);
		}
	}

	R = RecordCollection::record_from_record_type( "SSBOND  " );
	std::map<std::string, utility::vector1<SSBondInformation> >::const_iterator last_first_disulf = sfr.ssbond_map().end();
	for ( std::map<std::string, utility::vector1<SSBondInformation> >::const_iterator branch_point = sfr.ssbond_map().begin();
			branch_point != last_first_disulf; ++branch_point ) {
		utility::vector1<SSBondInformation> ssbonds = branch_point->second;
		Size n_ssbonds = ssbonds.size();
		for ( uint i = 1; i <= n_ssbonds; ++i ) {
			SSBondInformation ssbond = ssbonds[i];
			R["type"].value = "SSBOND  ";
			R["name1"].value = ssbond.name1;
			R["resName1"].value = ssbond.resName1;
			R["chainID1"].value = std::string( 1, ssbond.chainID1 );
			R["resSeq1"].value = pad_left( ssbond.resSeq1, 4 ); //("%4d", link.resSeq1);
			R["iCode1"].value = std::string( 1, ssbond.iCode1 );
			R["name2"].value = ssbond.name2;
			R["resName2"].value = ssbond.resName2;
			R["chainID2"].value = std::string( 1, ssbond.chainID2 );
			R["resSeq2"].value = pad_left( ssbond.resSeq2, 4 ); //("%4d", link.resSeq2);
			R["iCode2"].value = std::string( 1, ssbond.iCode2 );
			R["length"].value = fmt_real( ssbond.length, 2 , 2 ); //("%5.2f", link.length);
			VR.push_back(R);
		}
	}

	R = RecordCollection::record_from_record_type( "CRYST1" );
	CrystInfo ci = sfr.crystinfo();
	if ( ci.A() > 0 && ci.B() > 0 && ci.C() > 0 ) {
		R["type"].value = "CRYST1";
		R["a"].value =  fmt_real( ci.A(), 5, 3 ); //("%9.3f", ci.A());
		R["b"].value =  fmt_real( ci.B(), 5, 3 ); //("%9.3f", ci.B());
		R["c"].value =  fmt_real( ci.C(), 5, 3 ); //("%9.3f", ci.C());
		R["alpha"].value = fmt_real( ci.alpha(), 4, 2 ); //("%7.2f", ci.alpha());
		R["beta"].value =  fmt_real( ci.beta(), 4, 2 ); //print_d("%7.2f", ci.beta());
		R["gamma"].value = fmt_real( ci.gamma(), 4, 2 ); //print_d("%7.2f", ci.gamma());
		R["spacegroup"].value = ci.spacegroup();
		VR.push_back( R );
	}


	R = RecordCollection::record_from_record_type( "ATOM  " );
	for ( Size i=0; i<sfr.chains().size(); ++i ) {
		for ( Size j=0; j<sfr.chains()[i].size(); ++j ) {
			AtomInformation const & ai( sfr.chains()[i][j] );
			R["type"].value = ( ai.isHet ? "HETATM" : "ATOM  " );
			R["serial"].value = pad_left( ai.serial, 5 ); //)("%5d", ai.serial);
			R["name"].value = ai.name;
			R["resName"].value = ai.resName;
			std::string cid(" ");
			cid[0] = ai.chainID;
			R["chainID"].value = cid;
			R["resSeq"].value = pad_left( ai.resSeq, 4 ); //("%4d", ai.resSeq);
			R["iCode"].value = ai.iCode;
			R["x"].value = fmt_real( ai.x, 4, 3 ); //)("%8.3f", ai.x);
			R["y"].value = fmt_real( ai.y, 4, 3 ); //)("%8.3f", ai.y);
			R["z"].value = fmt_real( ai.z, 4, 3 ); //("%8.3f", ai.z);
			R["element"].value = ai.element;
			R["occupancy"].value = fmt_real( ai.occupancy, 3, 2 ); //("%6.2f", ai.occupancy);
			R["tempFactor"].value = fmt_real( ai.temperature, 3, 2 ); //("%6.2f", ai.temperature);
			R["segmentID"].value = ai.segmentID;
			VR.push_back(R);
		}
	}

	// Adding 'TER' line at the end of PDB.
	// TER lines are not supposed to go at the end of the PDB; they go at the end of a chain. ~Labonte
	Record T = RecordCollection::record_from_record_type( "TER   " );
	T["type"].value = "TER   ";
	VR.push_back(T);

	///////////////////////////////// Additional Information frm SFR /////////////////////////////////////////////
	//
	//

	// Connect Records
	for ( Size i=0; i<sfr.chains().size(); ++i ) {
		for ( Size j=0; j<sfr.chains()[i].size(); ++j ) {



			AtomInformation const & ai( sfr.chains()[i][j] );
			if ( ai.connected_indices.size() == 0 ) continue;

			core::Size current_L = 0;
			core::Size max_L = 4; //Number of conect records per line

			R = RecordCollection::record_from_record_type( "CONECT" );
			R["type"].value = "CONECT";
			R["serial0"].value = pad_left( ai.serial, 5);

			for ( core::Size c_index = 1; c_index <= ai.connected_indices.size(); ++c_index ) {
				current_L += 1;

				R["serial"+utility::to_string( current_L )].value = pad_left( ai.connected_indices[ c_index ], 5);

				if ( current_L == max_L ) {

					current_L = 0;
					VR.push_back(R);
					continue;
				}
			}

			VR.push_back(R);
		}
	}

	R = RecordCollection::record_from_record_type( "UNKNOW" );
	R["type"].value = "UNKNOW"; // note: "start" and "end" of this field ("type") are left at their default values of 0; this is not a proper field.

	// FoldTree
	if ( ! sfr.foldtree_string().empty() ) {
		R["info"].value = "REMARK " + sfr.foldtree_string();
		VR.push_back(R);
	}

	// PDB Comments
	if ( sfr.pdb_comments().size() > 0 ) {
		std::stringstream out;
		out << "##Begin comments##" << std::endl;
		using namespace std;
		map< string, string > const comments = sfr.pdb_comments();
		for ( std::map< string, string >::const_iterator i = comments.begin(); i != comments.end(); ++i ) {
			out << i->first<<" "<<i->second << std::endl;
		}
		out << "##End comments##" << std::endl;

		R["info"].value = out.str();
		//TR << out.str() << std::endl;
		VR.push_back(R);
	}

	// Additional Output String
	if ( ! sfr.additional_string_output().empty() ) {
		R["info"].value = sfr.additional_string_output();
		//TR << sfr.additional_string_output();

		VR.push_back(R);
	}

	return VR;
}

/*
void
old_dump_pdb(
pose::Pose const & pose,
std::ostream & out,
id::AtomID_Mask const & mask,
Size & atomno,
std::string const & tag,
char chain,
utility::vector1<Size> resnums
) {

}

*/
void
old_dump_pdb(
	pose::Pose const & pose,
	std::ostream & out,
	id::AtomID_Mask const & mask,
	std::string const & /* tag */
) {
	PoseToStructFileRepConverter converter = PoseToStructFileRepConverter();
	converter.init_from_pose( pose , mask );
	std::string data = create_pdb_contents_from_sfr(*converter.sfr());
	out.write( data.c_str(), data.size() );
}




///////////////////////////////////////////////////////////////////////////////
void
old_dump_pdb(
	pose::Pose const & pose,
	std::ostream & out,
	std::string const & tag
) {
	dump_pdb(pose, out, tag);
}


void
old_dump_pdb(
	pose::Pose const & pose,
	std::string const & filename,
	std::string const & tag
) {
	dump_pdb(pose, filename, tag);
}



/// @note Python compatible wrapper avoiding reference parameter
void
dump_pdb_residue(
	conformation::Residue const & rsd,
	std::ostream & out,
	Size start_atom_number)
{
	dump_pdb_residue(rsd, start_atom_number, out);
}

/// @details Create a faux PDBInfo object from a Residue
pose::PDBInfoOP
create_pdb_info_for_single_residue_pose(
	pose::Pose const & pose // the pose containing the single residue
)
{
	conformation::Residue const & rsd( pose.residue(1) );

	core::pose::PDBInfoOP pdb_info( new core::pose::PDBInfo( 1 ) );

	utility::vector1< int > pdb_numbering( 1, rsd.seqpos() );
	utility::vector1< char > pdb_chains( 1, char( ( rsd.chain() - 1 ) % 26 ) + 'A' );
	utility::vector1< char > insertion_codes( 1, ' ' );

	pdb_info->set_numbering( pdb_numbering );
	pdb_info->set_chains( pdb_chains );
	pdb_info->set_icodes( insertion_codes );
	pdb_info->resize_atom_records( pose );

	pdb_info->obsolete( false );

	return pdb_info;
}

void
dump_pdb_residue(
	conformation::Residue const & rsd,
	Size & atom_number,
	std::ostream & out
) {

	pose::Pose pose;
	pose.append_residue_by_jump( rsd, 1 );

	// create a PDBInfo object for this single-residue pose
	pose.pdb_info( create_pdb_info_for_single_residue_pose( pose ));

	pose_to_sfr::PoseToStructFileRepConverter converter;
	converter.append_residue_to_sfr(pose, 1, atom_number);
	std::string data = create_pdb_contents_from_sfr(*converter.sfr());
	out.write( data.c_str(), data.size() );
}


} //pdb
} //io
} //core
