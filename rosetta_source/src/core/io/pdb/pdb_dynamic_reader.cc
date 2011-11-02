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
///
/// @brief  PDB Dynamic reader
/// @author Sergey Lyskov (Sergey.Lyskov@jhu.edu)


// Unit headers
#include <core/io/pdb/pdb_dynamic_reader.hh>

//
#include <core/io/pdb/file_data.hh>
#include <core/pose/Remarks.hh>
#include <core/types.hh>

// Utility headers
#include <utility/tools/make_map.hh>
#include <numeric/xyzVector.hh>
#include <ObjexxFCL/string.functions.hh>
// AUTO-REMOVED #include <utility/vector0.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// C++ headers
#include <cstdlib>
// AUTO-REMOVED
#include <cstdio>
#include <algorithm>

#include <utility/vector1.hh>

//#include <cstdlib>
//#include <map>
//#include <vector>

static basic::Tracer TR("core.io.pdb.pdb_dynamic_reader");

namespace core {
namespace io {
namespace pdb {

using core::Size;
using core::SSize;

/// @details static holder for collection of Fields.
PDB_DReader::RecordRef PDB_DReader::PDB_Records_;

/// @details Debug printing, serialazing to Tracer like object.
std::ostream& operator <<(std::ostream &os,Record const & R)
{
	for(Record::const_iterator p=R.begin(); p!=R.end(); p++ ) {
		os << "<Record>{" << p->first << ":" << p->second
			<< "}\n";
	}

	return os;
}

/// @details check if records table was init, init table otherwise
/// return reference to private static records collection
PDB_DReader::RecordRef & PDB_DReader::getRecordCollection()
{
	if( PDB_Records_.size() == 0 ) {
		PDB_Records_ = utility::tools::make_map<String, Record>(
			"MODEL ", utility::tools::make_map<String, Field>(
				"type",			Field( 1,    6),
				"serial",		Field( 7,   80) ), /// String

			"HEADER", utility::tools::make_map<String, Field>(
				"type",           Field( 1,  6),
				"classification", Field(11, 50),
				"depDate",        Field(51, 59),
				"idCode",         Field(63, 66)	),

			"TITLE ", utility::tools::make_map<String, Field>(
				"type",         Field( 1,  6),
				"continuation", Field( 9, 10),
				"title",        Field(11, 70) ),

			"COMPND", utility::tools::make_map<String, Field>(
				"type",         Field( 1,  6),
				"continuation", Field( 9, 10),
				"compound",     Field(11, 70) ),

			"REMARK", utility::tools::make_map<String, Field>(
				"type",      Field( 1,  6),
				"remarkNum", Field( 8, 10),
				"value",     Field(12, 70) ), // non standard name.

			"ATOM  ", utility::tools::make_map<String, Field>(
				"type",       Field( 1,  6),
				"serial",     Field( 7, 11), /// Integer
				"name",       Field(13, 16), /// Atom
				"altLoc",     Field(17, 17), /// Character
				"resName",    Field(18, 20), /// Residue name
				"chainID",    Field(22, 22), /// Character
				"resSeq",     Field(23, 26), /// Integer
				"iCode",      Field(27, 27), /// AChar
				"x",          Field(31, 38), /// Real(8.3)
				"y",          Field(39, 46), /// Real(8.3)
				"z",          Field(47, 54), /// Real(8.3)
				"occupancy",  Field(55, 60), /// Real(6.2)
				"tempFactor", Field(61, 66), /// Real(6.2)
				///"segID",     Field(73, 76),
				"element",    Field(77, 78), /// LString(2)
				"charge",     Field(79, 80)  /// LString(2)
				),
			"HETATM", utility::tools::make_map<String, Field>(
				"type",       Field( 1,  6),
				"serial",     Field( 7, 11), /// Integer
				"name",       Field(13, 16), /// Atom
				"altLoc",     Field(17, 17), /// Character
				"resName",    Field(18, 20), /// Residue name
				"chainID",    Field(22, 22), /// Character
				"resSeq",     Field(23, 26), /// Integer
				"iCode",      Field(27, 27), /// AChar
				"x",          Field(31, 38), /// Real(8.3)
				"y",          Field(39, 46), /// Real(8.3)
				"z",          Field(47, 54), /// Real(8.3)
				"occupancy",  Field(55, 60), /// Real(6.2)
				"tempFactor", Field(61, 66), /// Real(6.2)
				///"segID",     Field(73, 76),
				"element",    Field(77, 78), /// LString(2)
				"charge",     Field(79, 80)  /// LString(2)
				),
			"SSBOND", utility::tools::make_map<String, Field>(
				"type",     Field( 1,  6),
				"serNum",   Field( 8, 10), /// Integer
				"CYS",      Field(12, 14),
				"chainID1", Field(16, 16), /// Character
				"seqNum1",  Field(18, 21), /// Integer
				"icode1",   Field(22, 22), /// AChar
				"CYS",      Field(26, 28),
				"chainID2", Field(30, 30), /// Character
				"seqNum2",  Field(32, 35), /// Integer
				"icode2",   Field(36, 36), /// AChar
				"sym1",     Field(60, 65), /// SymOP
				"sym2",     Field(67, 72)  /// SymOP
				),
			"TER   ", utility::tools::make_map<String, Field>(
				"type",       Field( 1,  6),
				"serial",     Field( 7, 11),
				"resName",    Field(18, 20),
				"chainID",    Field(22, 22),
				"resSeq",     Field(23, 26),
				"iCode",      Field(27, 27) ),

			"UNKNOW", utility::tools::make_map<String, Field>(
				"type", Field( 1,  6),
				"info", Field( 7, 80) )
			);
	}
	return PDB_Records_;
}


/// @details create Record Object with field collection (depending of the type information in _s),
/// and read fields values.
Record PDB_DReader::mapStringToRecord(const String & _s)
{
	//PDB_Records();
	getRecordCollection();

	//if( PDB_Records_.size() == 0 ) PDB_Records_ = getRecordCollection();
	String s(_s);
	s.resize(80, ' ');
	Field T = Field("type", 1,  6);

	T.getValueFrom(s);

	Record R;
	if( PDB_Records_.count(T.value) ) {	R = PDB_Records_[ T.value ]; }
	else { R = PDB_Records_["UNKNOW"]; }

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

/// @details Parse given PDB data (represented as a string) in to vector of Records.
std::vector<Record> PDB_DReader::parse(const String &pdb)
{
	runtime_assert(!pdb.empty()); //we're wasting time if there's no data here...
	std::vector<String> sl = split(pdb);
	std::vector<Record> r( sl.size() );
	std::transform(sl.begin(), sl.end(), r.begin(), PDB_DReader::mapStringToRecord);
	return r;
}

/// @details Create FileData object from a given vector of Recrods.
FileData PDB_DReader::createFileData(std::vector<Record> & VR)
{
	FileData fd;

	typedef std::map<char, AtomChain> ChainMap;
	ChainMap m;

	int terCount = 0;
	std::vector< char > chain_list; // preserve order
	std::map<char,Size> chain_to_idx;
	std::map<std::pair<Size,Size>,char> modelchain_to_chain;
	std::string chainletters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
	for(Size i = 0; i < chainletters.size(); ++i) {
		modelchain_to_chain[std::pair<Size,Size>(0,i)] = chainletters[i];
		modelchain_to_chain[std::pair<Size,Size>(1,i)] = chainletters[i];
	}
	Size modelidx = 1;
	bool modeltags_present = false;

	for(Size i=0; i<VR.size(); i++) {
		// jec reading multimodel PDBs
		if (VR[i]["type"].value == "MODEL " ) {
			// store the serial number as the filename, which will become the PDBInfo name of the pose
			std::string temp_model = ObjexxFCL::strip_whitespace( VR[i]["serial"].value ) ;
			fd.modeltag = temp_model.c_str();
			if( basic::options::option[basic::options::OptionKeys::in::file::new_chain_order]() ) {
				if(modeltags_present) {
					// second model... all chains should be present...
					for(Size model_idx=2;model_idx*chain_to_idx.size()<chainletters.size();++model_idx) {
						for(Size chain_idx=1; chain_idx <= chain_to_idx.size(); ++chain_idx) {
							TR << "REARRANGE CHAINS " << model_idx << " " << chain_idx << " " << (model_idx-1)*chain_to_idx.size()+chain_idx << std::endl;
							modelchain_to_chain[std::pair<Size,Size>(model_idx,chain_idx)] = chainletters[(model_idx-1)*chain_to_idx.size()+chain_idx-1];
						}
					}
					modelidx++;
					if(modelidx > 8) utility_exit_with_message("quitting: too many MODELs");
				} else {
					modeltags_present = true;
				}
			}
		} else if( VR[i]["type"].value == "ATOM  " || VR[i]["type"].value == "HETATM")  {
			Record & R(VR[i]);

			AtomInformation ai;
			ai.isHet = (R["type"].value == "HETATM");
			ai.serial = atoi( R["serial"].value.c_str() );
			ai.name = R["name"].value;
			ai.altLoc = 0; if( R["altLoc"].value.size() > 0 ) ai.altLoc = R["altLoc"].value[0];
			ai.resName = R["resName"].value;

			ai.chainID = 0;	if( R["chainID"].value.size() > 0 ) ai.chainID = R["chainID"].value[0];
			if( basic::options::option[basic::options::OptionKeys::in::file::new_chain_order]() ) {
				if( R["chainID"].value.size() > 0 ) {
					char chainid = R["chainID"].value[0];
					if( chain_to_idx.find(chainid) == chain_to_idx.end() ) {
						chain_to_idx[chainid] = chain_to_idx.size()+1;
						TR << "found new chain " << chainid << " " << chain_to_idx.size() << std::endl;
					}
					ai.chainID = modelchain_to_chain[std::pair<Size,Size>(modelidx,chain_to_idx[chainid])];
				}
			}

			ai.resSeq = atoi( R["resSeq"].value.c_str() );
			ai.iCode = 0; if( R["iCode"].value.size() > 0 ) ai.iCode = R["iCode"].value[0];

			// how can you check properly if something will successfuly convert to a number !?!?!?
			bool force_no_occupancy = false;
			if( R["x"].value == "     nan"){ai.x =0.0;force_no_occupancy=true;} else { ai.x = atof( R["x"].value.c_str() ); }
			if( R["y"].value == "     nan"){ai.y =0.0;force_no_occupancy=true;} else { ai.y = atof( R["y"].value.c_str() ); }
			if( R["z"].value == "     nan"){ai.z =0.0;force_no_occupancy=true;} else { ai.z = atof( R["z"].value.c_str() ); }

			// check that the occupancy column actually exists. If it doesnt, assume full occupancy.
			// otherwise read it.
			if( R["occupancy"].value == "      ")  ai.occupancy = 1.0;
			else                                   ai.occupancy = atof( R["occupancy"].value.c_str() );
			if(force_no_occupancy) ai.occupancy = -1.0;

			ai.temperature = atof( R["tempFactor"].value.c_str() );
			ai.element = R["element"].value;
			ai.terCount = terCount;

			m[ai.chainID].push_back(ai);
			if ( std::find( chain_list.begin(), chain_list.end(), ai.chainID ) == chain_list.end() ) {
				chain_list.push_back( ai.chainID );
			}
		} else if( VR[i]["type"].value == "TER   " || VR[i]["type"].value == "END   ")  {
			terCount++;
		} else if( (VR[i]["type"].value == "ENDMDL") &&
							 (basic::options::option[basic::options::OptionKeys::in::file::obey_ENDMDL].value()) )  {
		 	TR.Warning << "hit ENDMDL, not reading anything further" << std::endl;
			break;
		} else if( VR[i]["type"].value == "REMARK")  {
			pose::RemarkInfo ri;
			ri.num = atoi( VR[i]["remarkNum"].value.c_str() ),
			ri.value = VR[i]["value"].value;

			fd.remarks->push_back(ri);
		}
	}

	for ( Size i=0; i< chain_list.size(); ++i ) { // std::vector
		fd.chains.push_back( m.find( chain_list[i] )->second );
	}
// 	for(ChainMap::const_iterator p=m.begin(); p!=m.end(); p++ ) {
// 		fd.chains.push_back( (*p).second );
// 	}

	return fd;
}

/// @details Create FileData from a given PDB data (represented as a string).
FileData PDB_DReader::createFileData(const String & data)
{
		std::vector<Record> VR( parse(data) );
		return createFileData(VR);
}

/// @details create PDB string from Record data.
String PDB_DReader::createPDBString(const Record &R)
{
	String s(80, ' ');
	for(Record::const_iterator p=R.begin(); p!=R.end(); p++ ) {
		String v = p->second.value;  v.resize(p->second.end - p->second.start +1, ' ');
		s.replace( p->second.start-1, p->second.end - p->second.start +1, v);
	}
	return(s);
}

/// @details create PDB file (represented as a string) from FileData object.
String PDB_DReader::createPDBData(FileData const &fd)
{
	std::vector<Record> VR( PDB_DReader::createRecords(fd) );

	String r;  r.reserve(81*VR.size());
	for(Size i=0; i<VR.size(); i++) {
		//std::cout << VR[i] << '\n' <<  createPDBString( VR[i] ) <<  "\n";
		r += createPDBString( VR[i] ) + '\n';
	}
	return r;
}

utility::vector1< std::string > PDB_DReader::createPDBData_vector( FileData const & fd ) {
	std::vector<Record> VR( PDB_DReader::createRecords(fd) );

	utility::vector1< std::string > lines;
	lines.reserve( VR.size() );
	for ( Size i = 0; i < VR.size(); i++ ) {
		lines.push_back( createPDBString( VR[i] ) );
	}
	return lines;
}

/// @details print int with format to string
std::string print_i(const char *format, int I)
{
	std::string buf;  buf.resize(1024);
	sprintf(&buf[0], format, I);
	return buf;
}

/// @details print double with format to string
std::string print_d(const char *format, double d)
{
	std::string buf;  buf.resize(1024);
	sprintf(&buf[0], format, d);
	return buf;
}

/// @details Create vector of Record from given FileData object.
//  Used in PDB writing support.
std::vector<Record> PDB_DReader::createRecords(FileData const & fd)
{
	Record R = getRecordCollection()["REMARK"];
	std::vector<Record> VR;

	for(Size i=0; i<fd.remarks->size(); i++) {
		pose::RemarkInfo const & ri( fd.remarks->at(i) );

		R["type"].value = "REMARK";
		R["remarkNum"].value = print_i("%3d", ri.num);
		R["value"].value = ri.value;
		VR.push_back(R);
	}


	R = getRecordCollection()["ATOM  "];
	for(Size i=0; i<fd.chains.size(); i++) {
		for(Size j=0; j<fd.chains[i].size(); j++) {
			AtomInformation const & ai( fd.chains[i][j] );
			R["type"].value = (ai.isHet ? "HETATM" : "ATOM  ");
			R["serial"].value = print_i("%5d", ai.serial);
			R["name"].value = ai.name;
			R["resName"].value = ai.resName;
			std::string cid(" ");  cid[0] = ai.chainID;
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
			//std::cout << "ai.resSeq=" << ai.resSeq << " R=" << R["resSeq"].value << "\n";
		}
	}

	//R = getRecordCollection()["REMARK"];
	//for(Size i=0; i<fd.remarks->size(); i++) {
	//	RemarkInfo const & ri( fd.remarks->at(i) );

	//	R["type"].value = "REMARK";
	//	R["remarkNum"].value = print_i("%3d", ri.num);
	//	R["value"].value = ri.value;
	//	VR.push_back(R);
	//}

	// Adding 'TER' line at the end of PDB.
	Record T = getRecordCollection()["TER   "];
	T["type"].value = "TER   ";
	VR.push_back(T);

	return VR;
}

} // namespace pdb
} // namespace io
} // namespace core

