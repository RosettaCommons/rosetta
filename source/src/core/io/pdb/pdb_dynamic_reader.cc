// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pdb/pdb_dynamic_reader.hh
/// @brief  Method definitions for PDB Dynamic reader.
/// @author Sergey Lyskov (Sergey.Lyskov@jhu.edu)

// Note: DO NOT ACCESS THE OPTIONS SYSTEM DIRECTLY IN THIS FILE!
// Doing so will mean the Resource Manager will not work properly.
// Instead, modify PDB_DReaderOptions to include the option.


// Unit headers
#include <core/io/pdb/pdb_dynamic_reader.hh>
#include <core/io/pdb/pdb_dynamic_reader_options.hh>

// Package headers
#include <core/io/pdb/Field.hh>
#include <core/io/pdb/HeaderInformation.hh>
#include <core/io/pdb/file_data.hh>
#include <core/pose/Remarks.hh>

// Project headers
#include <core/types.hh>

// Basic headers
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/tools/make_map.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// External headers
#include <ObjexxFCL/string.functions.hh>

// C++ headers
#include <cstdlib>
#include <cstdio>
#include <algorithm>


static basic::Tracer TR("core.io.pdb.pdb_dynamic_reader");

namespace core {
namespace io {
namespace pdb {

using core::Size;
using core::SSize;

/// @details create Record Object with field collection (depending of the type information in _s),
/// and read fields values.
Record PDB_DReader::mapStringToRecord(const String & _s)
{
	RecordRef pdb_records(Field::getRecordCollection());

	String s(_s);
	s.resize(80, ' ');
	Field T = Field("type", 1,  6);

	T.getValueFrom(s);

	Record R;
	if( pdb_records.count(T.value) ) {	R = pdb_records[ T.value ]; }
	else { R = pdb_records["UNKNOW"]; }

	for(Record::iterator p=R.begin(); p!=R.end(); p++) (*p).second.getValueFrom(s);

	return R;
}

/// @details split String by new line symbols, return vector of string.
std::vector<String> split(const String &s)
{
	std::vector<String> r;
	Size start=0, i=0;
	while(start < s.size()) {
		if( s[i] == '\n' || s[i] == '\r' /* || i==s.size()-1 */) {
			r.push_back( String(s.begin()+start, s.begin()+i) );
			start = i+1;
		}
		i++;
		if( i == s.size() ) {
			r.push_back( std::string(s.begin()+start, s.begin()+i) );
			break;
		}
	}
	for(SSize i=r.size()-1; i>=0; i--) {  /// removing empty lines
		if( r[i].size() == 0 ) r.erase( r.begin()+i );
	}
	return r;
}

/// @details Parse given PDB data (represented as a string) into vector of Records.
std::vector<Record> PDB_DReader::parse(const String &pdb)
{
	runtime_assert(!pdb.empty()); //we're wasting time if there's no data here...
	std::vector<String> sl = split(pdb);
	std::vector<Record> r( sl.size() );
	std::transform(sl.begin(), sl.end(), r.begin(), PDB_DReader::mapStringToRecord);
	return r;
}

/// @details Create FileData object from a given vector of Records.
FileData PDB_DReader::createFileData(std::vector<Record> & VR)
{
	PDB_DReaderOptions options;
	return createFileData( VR, options );
}

/// @details Create FileData object from a given vector of Records.
FileData PDB_DReader::createFileData(std::vector<Record> & VR, PDB_DReaderOptions const & options)
{
	FileData fd;
	fd.initialize_header_information();

	bool read_pdb_header = options.read_pdb_header();

	std::map<char, AtomChain> atom_chain_map;
	Size ter_count = 0;
	std::vector< char > chain_list; // preserve order
	std::map<char, Size> chain_to_idx;

	std::map<std::pair<Size, Size>, char> modelchain_to_chain;
	std::string const chain_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
	for(Size i = 0; i < chain_letters.size(); ++i) {
		modelchain_to_chain[std::pair<Size, Size>(0, i)] = chain_letters[i];
		modelchain_to_chain[std::pair<Size, Size>(1, i)] = chain_letters[i];
	}
	Size modelidx = 1;
	bool modeltags_present = false;

	// Loop over all PDB records stored in vector VR.
	for(Size i=0; i<VR.size(); ++i) {
		std::string record_type = VR[i]["type"].value;

		// Record contains "header information", i.e., is from the Title Section of the PDB file.
		if (	record_type == "HEADER" || record_type == "KEYWDS" ||
				record_type == "TITLE " || record_type == "COMPND" ||
				record_type == "EXPDTA") {
			if( read_pdb_header ){
				fd.store_header_record(VR[i]);
			}

		// Record contains a remark from the Title Section of the PDB file.
		} else if( record_type == "REMARK")  {
			pose::RemarkInfo ri;
			ri.num = atoi( VR[i]["remarkNum"].value.c_str() ),
			ri.value = VR[i]["value"].value;

			//Added by DANIEL to skip reading the PDBinfo-LABEL, which comes "nfo-LABEL:"
			//Those are read in a different way, using: core.import_pose().read_additional_pdb_data()
			if( ( ri.value.size() >= 10 ) && ( ri.value.substr(0, 10) == "nfo-LABEL:" ) )
			{
				continue;
			}

			fd.remarks->push_back(ri);

		// Record contains heterogen nomenclature information from the Heterogen section of the PDB file.
		} else if (record_type == "HETNAM") {
			fd.store_heterogen_names(VR[i]["hetID"].value, VR[i]["text"].value);

		// Record contains nonstandard polymer linkage information from the Connectivity Annotation Section of the PDB
		//file.
		} else if (record_type == "LINK  ") {
			if (options.read_link_records()) {
				fd.store_link_record(VR[i]);
			}

		// Record contains multimodel PDBs as specified in the Coordinate Section of the PDB file..
		} else if (record_type == "MODEL " ) {
			// store the serial number as the filename, which will become the PDBInfo name of the pose
			std::string temp_model = ObjexxFCL::strip_whitespace( VR[i]["serial"].value ) ;
			fd.modeltag = temp_model.c_str();
			if( options.new_chain_order() ) {
				if(modeltags_present) {
					// second model... all chains should be present...
					for(Size model_idx=2; model_idx*chain_to_idx.size()<chain_letters.size(); ++model_idx) {
						for(Size chain_idx=1; chain_idx <= chain_to_idx.size(); ++chain_idx) {
							TR << "REARRANGE CHAINS " << model_idx << " " << chain_idx << " ";
							TR << (model_idx-1)*chain_to_idx.size()+chain_idx << std::endl;
							modelchain_to_chain[std::pair<Size, Size>(model_idx, chain_idx)] =
									chain_letters[(model_idx-1)*chain_to_idx.size() + chain_idx - 1];
						}
					}
					++modelidx;
					if(modelidx > 8) utility_exit_with_message("quitting: too many MODELs");
				} else {
					modeltags_present = true;
				}
			}

		// Record contains crystal information from the Crystallographic and Coordinate Transformation Section of the
		//PDB file.
		} else if( VR[i]["type"].value == "CRYST1")  {
			pose::CrystInfo ci;
			ci.A( atof( VR[i]["a"].value.c_str() ) );
			ci.B( atof( VR[i]["b"].value.c_str() ) );
			ci.C( atof( VR[i]["c"].value.c_str() ) );
			ci.alpha( atof( VR[i]["alpha"].value.c_str() ) );
			ci.beta( atof( VR[i]["beta"].value.c_str() ) );
			ci.gamma( atof( VR[i]["gamma"].value.c_str() ) );
			ci.spacegroup( VR[i]["spacegroup"].value );
			fd.crystinfo = ci;

		// Record contains atom information from the Coordinate Section of the PDB file.
		} else if( record_type == "ATOM  " || record_type == "HETATM")  {
			Record & R(VR[i]);

			AtomInformation ai;
			ai.isHet = (R["type"].value == "HETATM");
			ai.serial = atoi( R["serial"].value.c_str() );
			ai.name = R["name"].value;
			ai.altLoc = 0; if( R["altLoc"].value.size() > 0 ) ai.altLoc = R["altLoc"].value[0];
			ai.resName = R["resName"].value;

			ai.chainID = 0;	if( R["chainID"].value.size() > 0 ) ai.chainID = R["chainID"].value[0];
			if( options.new_chain_order() ) {
				if( R["chainID"].value.size() > 0 ) {
					char chainid = R["chainID"].value[0];
					if( chain_to_idx.find(chainid) == chain_to_idx.end() ) {
						chain_to_idx[chainid] = chain_to_idx.size()+1;
						TR << "found new chain " << chainid << " " << chain_to_idx.size() << std::endl;
					}
					ai.chainID = modelchain_to_chain[std::pair<Size, Size>(modelidx, chain_to_idx[chainid])];
				}
			}

			ai.resSeq = atoi( R["resSeq"].value.c_str() );
			ai.iCode = 0;
			if( R["iCode"].value.size() > 0 ) ai.iCode = R["iCode"].value[0];

			// how can you check properly if something will successfully convert to a number !?!?!?
			bool force_no_occupancy = false;
			if ( R["x"].value == "     nan") {
				ai.x =0.0;
				force_no_occupancy=true;
			} else {
				ai.x = atof( R["x"].value.c_str() );
			}
			if ( R["y"].value == "     nan") {
				ai.y =0.0;
				force_no_occupancy=true;
			} else {
				ai.y = atof( R["y"].value.c_str() );
			}
			if ( R["z"].value == "     nan") {
				ai.z =0.0;
				force_no_occupancy=true;
			} else {
				ai.z = atof( R["z"].value.c_str() );
			}

			// check that the occupancy column actually exists. If it doesn't, assume full occupancy.
			// otherwise read it.
			if( R["occupancy"].value == "      ") {
				ai.occupancy = 1.0;
			} else {
				ai.occupancy = atof( R["occupancy"].value.c_str() );
			}
			if(force_no_occupancy) ai.occupancy = -1.0;

			ai.temperature = atof( R["tempFactor"].value.c_str() );
			ai.element = R["element"].value;
			ai.terCount = ter_count;

			atom_chain_map[ai.chainID].push_back(ai);
			if ( std::find( chain_list.begin(), chain_list.end(), ai.chainID ) == chain_list.end() ) {
				chain_list.push_back( ai.chainID );
			}

		} else if( record_type == "TER   " || record_type == "END   ")  {
			++ter_count;

		} else if( (record_type == "ENDMDL") && (options.obey_ENDMDL()) )  {
		 	TR.Warning << "hit ENDMDL, not reading anything further" << std::endl;
			break;
		}
	}

	if( read_pdb_header ) {
		fd.finalize_header_information();
	}

	for ( Size i=0; i< chain_list.size(); ++i ) { // std::vector
		fd.chains.push_back( atom_chain_map.find( chain_list[i] )->second );
	}

	return fd;
}

/// @details Create FileData from a given PDB data (represented as a string).
FileData PDB_DReader::createFileData(const String & data)
{
	PDB_DReaderOptions options;
	return createFileData( data, options );
}

FileData PDB_DReader::createFileData(const String & data, PDB_DReaderOptions const & options)
{
	std::vector<Record> VR( parse(data) );
	return createFileData(VR, options);
}

/// @details create PDB string from Record data.
String PDB_DReader::createPDBString(const Record &R)
{
	String s(80, ' ');
	for(Record::const_iterator p=R.begin(); p!=R.end(); ++p ) {
		String v = p->second.value;
		v.resize(p->second.end - p->second.start +1, ' ');
		s.replace( p->second.start-1, p->second.end - p->second.start +1, v);
	}
	return(s);
}

/// @details create PDB file (represented as a string) from FileData object.
String PDB_DReader::createPDBData(FileData const &fd)
{
	std::vector<Record> VR( PDB_DReader::createRecords(fd) );

	String r;
	r.reserve(81*VR.size());
	for(Size i=0; i<VR.size(); ++i) {
		r += createPDBString( VR[i] ) + '\n';
	}
	return r;
}

utility::vector1< std::string >
PDB_DReader::createPDBData_vector( FileData const & fd ) {
	std::vector<Record> VR( PDB_DReader::createRecords(fd) );

	utility::vector1< std::string > lines;
	lines.reserve( VR.size() );
	for ( Size i = 0; i < VR.size(); i++ ) {
		lines.push_back( createPDBString( VR[i] ) );
	}
	return lines;
}


/// @note This function returns a string sized to 1024 characters!
/// You will likely want to resize the returned string. ~Labonte
std::string
print_i(const char *format, int I)
{
	std::string buf;
	buf.resize(1024);
	sprintf(&buf[0], format, I);
	return buf;
}

/// @note This function returns a string sized to 1024 characters!
/// You will likely want to resize the returned string. ~Labonte
std::string
print_d(const char *format, double d)
{
	std::string buf;
	buf.resize(1024);
	sprintf(&buf[0], format, d);
	return buf;
}


/// @details Create vector of Record from given FileData object.  Used in PDB writing support.
std::vector<Record> PDB_DReader::createRecords(FileData const & fd)
{
	std::vector<Record> VR;

	if(fd.header_information()){
		fd.fill_header_records(VR);
	}

	Record R = Field::getRecordCollection()["REMARK"];
	for(Size i=0; i<fd.remarks->size(); ++i) {
		pose::RemarkInfo const & ri( fd.remarks->at(i) );

		R["type"].value = "REMARK";
		R["remarkNum"].value = print_i("%3d", ri.num);
		R["value"].value = ri.value;
		VR.push_back(R);
	}

	R = Field::getRecordCollection()["HETNAM"];
	Size n_het_names = fd.heterogen_names.size();
	for (uint i = 1; i <= n_het_names; ++i) {
		R["type"].value = "HETNAM";
		R["continuation"].value = "  ";  // TODO: Wrap long text fields.
		R["hetID"].value = fd.heterogen_names[i].first;
		R["text"].value = fd.heterogen_names[i].second;
		VR.push_back(R);
	}

	R = Field::getRecordCollection()["LINK  "];
	std::map<std::string, utility::vector1<LinkInformation> >::const_iterator last_branch_point = fd.link_map.end();
	for (std::map<std::string, utility::vector1<LinkInformation> >::const_iterator branch_point = fd.link_map.begin();
			branch_point != last_branch_point; ++branch_point) {
		utility::vector1<LinkInformation> links = branch_point->second;
		Size n_links = links.size();
		for (uint i = 1; i <= n_links; ++i) {
			LinkInformation link = links[i];
			R["type"].value = "LINK  ";
			R["name1"].value = link.name1;
			R["resName1"].value = link.resName1;
			R["chainID1"].value = std::string(1, link.chainID1);
			R["resSeq1"].value = print_i("%4d", link.resSeq1);
			R["iCode1"].value = std::string(1, link.iCode1);
			R["name2"].value = link.name2;
			R["resName2"].value = link.resName2;
			R["chainID2"].value = std::string(1, link.chainID2);
			R["resSeq2"].value = print_i("%4d", link.resSeq2);
			R["iCode2"].value = std::string(1, link.iCode2);
			R["length"].value = print_d("%5.2f", link.length);
			VR.push_back(R);
		}
	}

	R = Field::getRecordCollection()["CRYST1"];
	pose::CrystInfo ci = fd.crystinfo;
	if (ci.A() > 0 && ci.B() > 0 && ci.C() > 0) {
		R["type"].value = "CRYST1";
		R["a"].value =  print_d("%9.3f", ci.A());
		R["b"].value =  print_d("%9.3f", ci.B());
		R["c"].value =  print_d("%9.3f", ci.C());
		R["alpha"].value = print_d("%7.2f", ci.alpha());
		R["beta"].value =  print_d("%7.2f", ci.beta());
		R["gamma"].value = print_d("%7.2f", ci.gamma());
		R["spacegroup"].value = ci.spacegroup();
		VR.push_back( R );
	}


	R = Field::getRecordCollection()["ATOM  "];
	for(Size i=0; i<fd.chains.size(); ++i) {
		for(Size j=0; j<fd.chains[i].size(); ++j) {
			AtomInformation const & ai( fd.chains[i][j] );
			R["type"].value = (ai.isHet ? "HETATM" : "ATOM  ");
			R["serial"].value = print_i("%5d", ai.serial);
			R["name"].value = ai.name;
			R["resName"].value = ai.resName;
			std::string cid(" ");
			cid[0] = ai.chainID;
			R["chainID"].value = cid;
			R["resSeq"].value = print_i("%4d", ai.resSeq);
			R["iCode"].value = ai.iCode;
			R["x"].value = print_d("%8.3f", ai.x);
			R["y"].value = print_d("%8.3f", ai.y);
			R["z"].value = print_d("%8.3f", ai.z);
			R["element"].value = ai.element;
			R["occupancy"].value = print_d("%6.2f", ai.occupancy);
			R["tempFactor"].value = print_d("%6.2f", ai.temperature);
			VR.push_back(R);
		}
	}

	// Adding 'TER' line at the end of PDB.
	// TER lines are not supposed to go at the end of the PDB; they go at the end of a chain. ~Labonte
	Record T = Field::getRecordCollection()["TER   "];
	T["type"].value = "TER   ";
	VR.push_back(T);

	return VR;
}

} // namespace pdb
} // namespace io
} // namespace core
