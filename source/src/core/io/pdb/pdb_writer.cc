// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/pdb/pdb_writer.cc
/// @brief Function(s) for pdb writing
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com), XRW 2016 Team
/// @author Sergey Lyskov (Sergey.Lyskov@jhu.edu)


// Unit headers
#include <core/io/pdb/pdb_writer.hh>
#include <core/io/pdb/Field.hh>
#include <core/io/pdb/RecordCollection.hh>

// Package headers
#include <core/io/util.hh>
#include <core/io/Remarks.hh>
#include <core/io/HeaderInformation.hh>
#include <core/io/CrystInfo.hh>
#include <core/io/pose_to_sfr/PoseToStructFileRepConverter.hh>
#include <core/io/StructFileRepOptions.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/scoring/dssp/Dssp.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/tools/make_map.hh>
#include <utility/string_util.hh>
#include <utility/io/ozstream.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// Options
#include <basic/options/option.hh>
//#include <basic/options/keys/chemical.OptionKeys.gen.hh>
//#include <basic/options/keys/run.OptionKeys.gen.hh>
//#include <basic/options/keys/in.OptionKeys.gen.hh>
//#include <basic/options/keys/mp.OptionKeys.gen.hh>
//#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
//#include <basic/options/keys/packing.OptionKeys.gen.hh>

// Basic headers
#include <basic/Tracer.hh>

// External headers
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

#include <basic/Tracer.hh>


// C++ headers
#include <cstdlib>
#include <cstdio>
#include <algorithm>


static basic::Tracer TR( "core.io.pdb.pdb_writer" );


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

using basic::Error;
using basic::Warning;

/// @brief Writes  <pose>  to a PDB file, returns false if an error occurs
///  Use default StructFileRepOptions
bool
dump_pdb(
	core::pose::Pose const & pose,
	std::string const & file_name
){
	StructFileRepOptionsCOP options=StructFileRepOptionsCOP( new core::io::StructFileRepOptions );
	return dump_pdb(pose, file_name, options);
}

/// @details Convert given Pose object into PDB format and save it to 'file_name' file.
/// @return true if operation was completed without error, false otherwise.
bool
dump_pdb(
	core::pose::Pose const & pose,
	String const & file_name,
	core::io::StructFileRepOptionsCOP options
) {
	utility::io::ozstream file(file_name.c_str(), std::ios::out | std::ios::binary);
	if ( !file ) {
		TR.Error << "StructFileRep::dump_pdb: Unable to open file:" << file_name << " for writing!!!" << std::endl;
		return false;
	}
	dump_pdb(pose, file, options, file_name );
	file.close();

	return true;
}

/// @brief Writes  <pose>  to a given stream in PDB file format
///  Use default StructFileRepOptions
void
dump_pdb(
	core::pose::Pose const & pose,
	std::ostream & out,
	std::string const & out_fname
){
	StructFileRepOptionsCOP options=StructFileRepOptionsCOP( new core::io::StructFileRepOptions );
	dump_pdb(pose, out, options, out_fname);
}

/// @details Convert given Pose object in to PDB format and send it to the given stream.
void
dump_pdb(
	core::pose::Pose const & pose,
	std::ostream & out,
	core::io::StructFileRepOptionsCOP options,
	std::string const & out_fname
) {
	String data;

	io::pose_to_sfr::PoseToStructFileRepConverter pu( *options );
	pu.init_from_pose( pose );
	if ( out_fname != "" ) {
		pu.sfr()->score_table_filename() = out_fname;
	}

	data = create_pdb_contents_from_sfr(*pu.sfr(), options);
	out.write( data.c_str(), data.size() );
}


/// @brief This version takes an AtomID mask.
/// @details Used by Will's motif hash stuff, I think.
void
dump_pdb(
	pose::Pose const & pose,
	std::ostream & out,
	id::AtomID_Mask const & mask,
	core::io::StructFileRepOptionsCOP options
) {
	PoseToStructFileRepConverter converter( *options );
	converter.init_from_pose( pose, mask );
	std::string const data( create_pdb_contents_from_sfr(*converter.sfr(), options) );
	out.write( data.c_str(), data.size() );
}

/// @details Convert given Pose object into PDB format and send it to the given stream.
/// only the residues corresponding to indices in 'residue_indices' will be output
void
dump_pdb(
	core::pose::Pose const & pose,
	std::ostream & out,
	utility::vector1< core::Size > const & residue_indices,
	core::io::StructFileRepOptionsCOP options
) {
	io::pose_to_sfr::PoseToStructFileRepConverter pu( *options );
	String data;
	pu.init_from_pose( pose, residue_indices );

	data = create_pdb_contents_from_sfr(*pu.sfr(), options);
	out.write( data.c_str(), data.size() );
}


/// @brief Writes poses to a single PDB file, returns false if an error occurs
///  Use default StructFileRepOptions
bool
dump_multimodel_pdb(
	utility::vector1< core::pose::PoseCOP > const & poses,
	std::string const & file_name
){
	StructFileRepOptionsCOP options=StructFileRepOptionsCOP( new core::io::StructFileRepOptions );
	return dump_multimodel_pdb(poses, file_name, options);
}

/// @brief Writes poses to a single PDB file, returns false if an error occurs
bool
dump_multimodel_pdb(
	utility::vector1< core::pose::PoseCOP > const & poses,
	String const & file_name,
	core::io::StructFileRepOptionsCOP options
) {
	utility::io::ozstream file(file_name.c_str(), std::ios::out | std::ios::binary);
	if ( !file ) {
		Error() << "StructFileRep::dump_multimodel_pdb: Unable to open file:" << file_name << " for writing!!!" << std::endl;
		return false;
	}
	dump_multimodel_pdb(poses, file, options);
	file.close();

	return true;
}

/// @brief Writes poses to a given stream in PDB file format
///  Use default StructFileRepOptions
void
dump_multimodel_pdb(
	utility::vector1< core::pose::PoseCOP > const & poses,
	std::ostream & out
){
	StructFileRepOptionsCOP options=StructFileRepOptionsCOP( new core::io::StructFileRepOptions );
	dump_multimodel_pdb(poses, out, options);
}

/// @brief Writes poses to a given stream in PDB file format
void
dump_multimodel_pdb(
	utility::vector1< core::pose::PoseCOP > const & poses,
	std::ostream & out,
	core::io::StructFileRepOptionsCOP options
) {
	// multimodel PDB method stolen from protocols/canonical_sampling/PDBTrajectoryRecorder.cc
	for ( Size i = 1; i <= poses.size(); i++ ) {
		pose::PoseCOP pose = poses[i];
		std::stringstream model_tag;
		model_tag << std::setfill('0') << std::setw( utility::get_num_digits( poses.size() ) ) << i;

		add_to_multimodel_pdb( *pose, out, model_tag.str(), options);

	}
}

/// @brief Adds a pose to a multimodel pdb file (or creates the file)
///  Use default StructFileRepOptions
bool
add_to_multimodel_pdb(
	core::pose::Pose const & pose,
	std::string const & file_name,
	std::string const & model_tag,
	bool clear_existing_structures /* = false */
) {
	StructFileRepOptionsCOP options=StructFileRepOptionsCOP( new core::io::StructFileRepOptions );
	return add_to_multimodel_pdb( pose, file_name, model_tag, options, clear_existing_structures );
}

/// @brief Adds a pose to a multimodel pdb file (or creates the file)
bool
add_to_multimodel_pdb(
	core::pose::Pose const & pose,
	std::string const & file_name,
	std::string const & model_tag,
	core::io::StructFileRepOptionsCOP options,
	bool clear_existing_structures /* = false */
) {
	utility::io::ozstream file( file_name.c_str(), std::ios::out | std::ios::binary |
		( clear_existing_structures ? std::ios::trunc : std::ios::app ) );
	if ( !file ) {
		Error() << "StructFileRep::add_to_multimodel_pdb: Unable to open file:" << file_name << " for writing!!!" << std::endl;
		return false;
	}

	add_to_multimodel_pdb( pose, file, model_tag, options );

	file.close();
	return true;
}

/// @brief Adds a pose to a multimodel pdb file
///  Use default StructFileRepOptions
void
add_to_multimodel_pdb(
	core::pose::Pose const & pose,
	std::ostream & out,
	std::string const & model_tag
) {
	StructFileRepOptionsCOP options=StructFileRepOptionsCOP( new core::io::StructFileRepOptions );
	add_to_multimodel_pdb( pose, out, model_tag, options );
}

/// @brief Adds a pose to a multimodel pdb file
void
add_to_multimodel_pdb(
	core::pose::Pose const & pose,
	std::ostream & out,
	std::string const & model_tag,
	core::io::StructFileRepOptionsCOP options
) {
	out << "MODEL     " << model_tag << std::endl;
	dump_pdb( pose, out, options );
	out << "ENDMDL" << std::endl;
}




/// @brief create a full pdb as a string given a StructFileRep object
std::string
create_pdb_contents_from_sfr(
	StructFileRep const & sfr,
	core::io::StructFileRepOptionsCOP options
) {
	std::string pdb_contents;

	utility::vector1< Record > records( create_records_from_sfr( sfr , options ) );

	pdb_contents.reserve( 81 * records.size() );
	for ( Size i = 1, imax=records.size(); i <= imax; ++i ) {
		pdb_contents += create_pdb_line_from_record( records[ i ] ) + "\n";
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
	for ( auto const & field : record ) {
		std::string value( field.second.value );
		if ( field.second.start == 0 ) continue; // Fields that are not properly defined (e.g. UNKNOW's "type" field) should not be written.
		if ( value.length() == 0 ) continue; // Cannot output a value with zero length
		//std::cout << "value: \"" << value << "\"" << std::endl;
		if ( field.second.end == 0 ) {  // If this field has no size limit.
			Size len = value.length();
			if ( field.second.start + len > 80 ) {
				line.resize( field.second.start + len - 1 );
			}
			line.replace( field.second.start - 1, len, value );
		} else {
			value.resize( field.second.end - field.second.start + 1, ' ' );
			line.replace( field.second.start - 1, field.second.end - field.second.start + 1, value );
		}
	}
	return line;
}


/// @details Create vector of Record from given StructFileRep object.  Used in PDB writing support.
std::vector<Record>
create_records_from_sfr(
	StructFileRep const & sfr,
	core::io::StructFileRepOptionsCOP options
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::vector<Record> VR;

	sfr.header()->fill_records( VR );
	Record R;

	// Title Section //////////////////////////////////////////////////////////
	R = RecordCollection::record_from_record_type( REMARK );
	for ( Size i=0; i<sfr.remarks()->size(); ++i ) {
		RemarkInfo const & ri( sfr.remarks()->at(i) );
		R["type"].value = "REMARK";
		R["remarkNum"].value = pad_left( ri.num, 3 ); //("%3d", ri.num);
		R["value"].value = ri.value;
		VR.push_back(R);
	}

	// Secondary Structure Section ///////////////////////////////////////////
	if ( options->output_secondary_structure() ) {

		R = RecordCollection::record_from_record_type( HELIX );
		//type is utility::vector1< HELIXInformation >
		auto const helixvec(sfr.HELIXInformations());
		for ( auto const & it : helixvec ) {
			R["type"].value = "HELIX ";
			R["serNum"].value = pad_left(it.helixID, 3);
			R["helixID"].value = it.helix_name;
			R["initResName"].value = it.name3_1;
			R["initChainID"].value = std::string( 1, it.chainID1 );
			R["initSeqNum"].value = pad_left( it.seqNum1, 4 );
			R["initICode"].value = std::string( 1, it.icode1 );
			R["endResName"].value = it.name3_2;
			R["endChainID"].value = std::string( 1, it.chainID2 );
			R["endSeqNum"].value = pad_left( it.seqNum2, 4 );
			R["endICode"].value = std::string( 1, it.icode2 );
			R["helixClass"].value = pad_left(it.helixClass, 2 );
			R["comment"].value = it.comment;
			R["length"].value = pad_left( it.length, 5 );
			VR.push_back(R);
		}

		R = RecordCollection::record_from_record_type( SHEET );
		//type is utility::vector1< SHEETInformation >
		auto const sheetvec(sfr.SHEETInformations());
		for ( auto const & it : sheetvec ) {
			R["type"].value = "SHEET ";
			R["strand"].value = pad_left(it.strand_num, 3);
			R["sheetID"].value = it.sheetID;
			R["numStrands"].value = pad_left( it.num_strands, 2);
			R["initResName"].value = it.name3_1;
			R["initChainID"].value = std::string( 1, it.chainID1 );
			R["initSeqNum"].value = pad_left( it.seqNum1, 4 );
			R["initICode"].value = std::string( 1, it.icode1 );
			R["endResName"].value = it.name3_2;
			R["endChainID"].value = std::string( 1, it.chainID2 );
			R["endSeqNum"].value = pad_left( it.seqNum2, 4 );
			R["endICode"].value = std::string( 1, it.icode2 );
			R["sense"].value = pad_left(it.strandClass, 2);
			//There are other fields, but we have no values for them, let's try skipping them wholesale and see what happens
			VR.push_back(R);
		}

	}



	// Heterogen Section //////////////////////////////////////////////////////
	R = RecordCollection::record_from_record_type( HETNAM );
	if ( options->use_pdb_format_HETNAM_records() ) {
		// TODO: Also have this option output HET, HETSYM, and FORMUL records.
		for ( auto const & elem : sfr.heterogen_names() ) {
			R[ "type" ].value = "HETNAM";
			R[ "continuation" ].value = "  ";  // TODO: Wrap long text fields.
			R[ "hetID" ].value = elem.first;
			R[ "text" ].value = elem.second;
			VR.push_back( R );
		}
	} else if ( !options->write_glycan_pdb_codes() ) {  // Use the Rosetta HETNAM format, which specifies base ResidueTypes.
		for ( auto const & elem : sfr.residue_type_base_names() ) {
			// The 6-character resID key for the map has a fixed format.
			std::string const & chainID( elem.first.substr( 5, 1 ) );  // 6th character
			std::string const & resSeq( elem.first.substr( 0, 4 ) );  // 1st through 4th characters
			std::string const & iCode( elem.first.substr( 4, 1 ) );  // 5th character
			std::string const & base_name( elem.second.second );

			R[ "type" ].value = "HETNAM";
			R[ "continuation" ].value = "  ";  // unused
			R[ "hetID" ].value = elem.second.first;  // 3-letter code
			R[ "text" ].value = chainID + resSeq + iCode + " " + base_name;
			VR.push_back( R );
		}
	}

	// Connectivity Annotation Section ////////////////////////////////////////
	R = RecordCollection::record_from_record_type( SSBOND );
	for ( auto const & branch_point : sfr.ssbond_map() ) {
		utility::vector1<SSBondInformation> ssbonds = branch_point.second;
		for ( SSBondInformation const & ssbond : ssbonds ) {
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

	R = RecordCollection::record_from_record_type( LINK );
	for ( auto const & branch_point : sfr.link_map() ) {
		utility::vector1<LinkInformation> links = branch_point.second;
		for ( LinkInformation const & link : links ) {
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

	// Crystallographic & Coordinate Transformation Section ///////////////////
	R = RecordCollection::record_from_record_type( CRYST1 );
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

	bool const no_chainend_ter( options->no_chainend_ter() ); //Should we skip TER records at chain ends?
	std::map < core::Size, core::Size > serial_to_serial_with_ter; //Reused later during CONECT record dumping.
	R = RecordCollection::record_from_record_type( "ATOM  " );
	Record T = RecordCollection::record_from_record_type( "TER   " );

	// A map of chain=>[name3] which will be used to assemble the SEQRES section
	std::map< std::string, utility::vector1<std::string> > seq_res_map;
	T["type"].value = "TER   ";

	int last_resseq(0);

	for ( Size i=0; i<sfr.chains().size(); ++i ) {
		for ( Size j=0; j<sfr.chains()[i].size(); ++j ) {
			AtomInformation const & ai( sfr.chains()[i][j] );
			// WTF CHAIN?
			//TR << "create records from sfr " << i << " " << j << " " << ai.chainID << std::endl;
			R["type"].value = ( ai.isHet ? "HETATM" : "ATOM  " );
			runtime_assert( !serial_to_serial_with_ter.count(ai.serial) );
			serial_to_serial_with_ter[ai.serial] = ai.serial + (no_chainend_ter ? 0 : ai.terCount); //Reused later during CONECT record dumping.
			R["serial"].value = pad_left( serial_to_serial_with_ter.at(ai.serial), 5 ); //)("%5d", ai.serial);
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

			if ( ai.resSeq != last_resseq ) {
				seq_res_map[cid].push_back(ai.resName);
			}

			last_resseq = ai.resSeq;

		}
		if ( ! no_chainend_ter && i > 0 ) {
			VR.push_back( T );
		}
	}
	if ( no_chainend_ter ) {
		VR.push_back( T );  //Put a TER record at the end of the ATOM lines if we're using the no_chainend_ter option.
	}

	// Sequence section
	if ( options->write_seqres_records() ) {

		core::Size seqres_serial(1);
		for ( auto & chain_sequence : seq_res_map ) {
			std::string const & chain(chain_sequence.first);
			utility::vector1<std::string> const & sequence(chain_sequence.second);
			core::Size const & chain_length(sequence.size());
			core::Size chain_index(1);

			while ( chain_index <= chain_length ) {
				R = RecordCollection::record_from_record_type( SEQRES );
				R["type"].value = "SEQRES";
				R["serNum"].value = pad_left(seqres_serial, 3);
				R["chainID"].value = chain;
				R["numRes"].value = pad_left(chain_length, 4);
				for ( core::Size seqres_slot = 1; seqres_slot <= 13; ++ seqres_slot ) {
					R["resName"+utility::to_string(seqres_slot)].value = sequence[chain_index];
					chain_index += 1;
					if ( chain_index > chain_length ) {
						break;
					}
				}
				seqres_serial++;
				VR.push_back(R);
			}
		}
	}

	// Connectivity Section ///////////////////////////////////////////////////
	for ( Size i=0; i<sfr.chains().size(); ++i ) {
		for ( Size j=0; j<sfr.chains()[i].size(); ++j ) {
			AtomInformation const & ai( sfr.chains()[i][j] );
			if ( ai.connected_indices.size() == 0 ) { continue; }

			core::Size current_L = 0;
			core::Size max_L = 4; //Number of conect records per line

			R = RecordCollection::record_from_record_type( CONECT );
			R["type"].value = "CONECT";
			R["serial0"].value = pad_left( ai.serial + (no_chainend_ter ? 0 : ai.terCount) , 5);

			for ( core::Size c_index = 1; c_index <= ai.connected_indices.size(); ++c_index ) {
				current_L += 1;

				//TODO -- correct the ai.connected_indices[ c_index ] in the next line.  Need to figure out how many TER records precede connected_index[c_index] and add this to the index.
				runtime_assert( serial_to_serial_with_ter.count(ai.connected_indices[ c_index ] ) );
				R["serial"+utility::to_string( current_L )].value = pad_left( serial_to_serial_with_ter.at( ai.connected_indices[ c_index ] ), 5);

				if ( current_L == max_L ) {
					current_L = 0;
					VR.push_back(R);
					R = RecordCollection::record_from_record_type( CONECT );
					R["type"].value = "CONECT";
					R["serial0"].value = pad_left( ai.serial + (no_chainend_ter ? 0 : ai.terCount) , 5);
					continue;
				}
			}
			if ( current_L != 0 ) {
				// if current_L == 0, then we have already pushed back this line and we shouldn't do it again
				VR.push_back(R);
			}
		}
	}

	// Rosetta-specific Information from the SFR //////////////////////////////
	R = RecordCollection::record_from_record_type( UNKNOW );
	// note: "start" and "end" of this field ("type") are left at their default values of 0; this is not a proper field.
	R["type"].value = "UNKNOW";

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
		for ( auto const & comment : comments ) {
			out << comment.first<<" "<<comment.second << std::endl;
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

	// Pose Energies Table (option system call preserves previous behavior)
	if ( sfr.score_table_labels().size() > 0 && sfr.score_table_lines().size() >0 && options->output_pose_energies_table() ) {
		R["info"].value = core::io::pose_energies_from_sfr(sfr);
		VR.push_back(R);
	}

	// Pose Arbitrary String and Float Data.
	if ( ( sfr.pose_cache_real_data().size() > 0 || sfr.pose_cache_string_data().size() > 0 ) && options->output_pose_cache() ) {
		R["info"].value = core::io::pose_data_cache_from_sfr(sfr);
		VR.push_back(R);
	}

	return VR;
}




/// @details Create a faux PDBInfo object from a Residue
pose::PDBInfoOP
create_pdb_info_for_single_residue_pose(
	pose::Pose const & pose, // the pose containing the single residue
	core::io::StructFileRepOptionsCOP //options
) {
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

/// @note Python compatible wrapper avoiding reference parameter
void
dump_pdb_residue(
	conformation::Residue const & rsd,
	std::ostream & out,
	core::Size start_atom_number,
	core::io::StructFileRepOptionsCOP options
) {
	dump_pdb_residue(rsd, start_atom_number, out, options);
}

std::string
dump_pdb_residue(
	conformation::Residue const & rsd,
	core::io::StructFileRepOptionsCOP options,
	core::Size start_atom_number

){
	return dump_pdb_residue(rsd, start_atom_number, options);
}

void
dump_pdb_residue(
	conformation::Residue const & rsd,
	core::Size & atom_number,
	std::ostream & out,
	core::io::StructFileRepOptionsCOP options
) {
	std::string data = dump_pdb_residue(rsd, atom_number, options);
	out.write( data.c_str(), data.size() );
}

std::string
dump_pdb_residue(
	conformation::Residue const & rsd,
	core::Size & atom_number,
	core::io::StructFileRepOptionsCOP options
) {

	pose::Pose pose;
	pose.append_residue_by_jump( *(rsd.clone()), 1 );

	// create a PDBInfo object for this single-residue pose
	pose.pdb_info( create_pdb_info_for_single_residue_pose( pose, options ));

	pose_to_sfr::PoseToStructFileRepConverter converter(*options);
	converter.append_residue_to_sfr(pose, 1, atom_number, 0);
	std::string data = create_pdb_contents_from_sfr(*converter.sfr(), options);
	return data;

}


///////////////////////////////////////////////////////////////
//////// LEGACY JD2-layer compatability (DO NOT USE) //////////
///////////////////////////////////////////////////////////////

void
dump_pdb(
	core::pose::Pose const & pose,
	std::string const &jd2_job_data,
	utility::io::ozstream & out,
	std::string const &filename
) {
	std::string sfr_string;
	dump_pdb( pose, jd2_job_data, sfr_string, filename);
	out << sfr_string;
}

void
dump_pdb(
	core::pose::Pose const & pose,
	std::string const &jd2_job_data,
	std::string & out,
	std::string const &filename
) {

	core::io::StructFileRepOptionsOP options=core::io::StructFileRepOptionsOP( new core::io::StructFileRepOptions );

	// OK, what's here?
	io::pose_to_sfr::PoseToStructFileRepConverter pose_to_sfr( *options );
	pose_to_sfr.init_from_pose( pose );
	pose_to_sfr.sfr()->score_table_filename() = filename;

	if ( !jd2_job_data.empty() ) pose_to_sfr.sfr()->append_to_additional_string_output( jd2_job_data );

	out = create_pdb_contents_from_sfr( *( pose_to_sfr.sfr() ), options );
}


}  // namespace pdb
}  // namespace io
}  // namespace core
