// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   devel/InvKinLigLoopDesign/Fragment.cc
///
/// @brief
/// @author

#include <climits>
#include <iostream>
#include <fstream>
#include <cstdlib> //required by GCC 4.3.2
#include <devel/inv_kin_lig_loop_design/Fragment.hh>

// #include <Util/Macros.hh>
// #include <Util/Erg.hh>
#include <devel/inv_kin_lig_loop_design/std_extra.hh>

#ifdef WIN32
#define _USE_MATH_DEFINES
#endif
#include <cmath>

#include <numeric/random/random.hh>
#include <basic/database/open.hh>

namespace {
using namespace std;

void open(ifstream& fin, const string& filename ) { // need to gzip stuff before I can use it...
	fin . close();
	fin.clear();
	fin.open( filename.c_str() );
	debug_assert( fin );
}

}

namespace devel {

namespace inv_kin_lig_loop_design {

namespace Fragment {

#define RAD2DEG (180.0/M_PI)
#define DEG2RAD (M_PI/180.0)

// ____________________ ResEntry ____________________

//       void ResEntry::reportRama(ostream& out) const {
// // const double phi = RAD2DEG*this->phi;
// // const double psi = RAD2DEG*this->psi;
// // const double ohm = RAD2DEG*this->ohm;
// //out << REPORT5(aa,ss,phi,psi,ohm) << "\n";
//       } // ResEntry::reportRama

istream& operator>>(istream& in, ResEntry& e) {
	in >> e.pdbid >> e.aa >> e.ss >> e.phi >> e.psi >> e.ohm;
	in.ignore(INT_MAX,'\n');
	e.phi *= DEG2RAD;
	e.psi *= DEG2RAD;
	e.ohm *= DEG2RAD;
	return in;
} // operator>>

ostream& operator<<(ostream& out, const ResEntry& e) {
	out << e.pdbid << " " << e.aa << " " << e.ss << " " << RAD2DEG*e.phi << " " << RAD2DEG*e.psi << " " << RAD2DEG*e.ohm << "\n";
	return out;
} // operator<<

// ____________________ Entry ____________________

istream& operator>>(istream& in, Entry& e) {
	e.vResEntries.clear();
	int n;
	in >> n;
	in.ignore(INT_MAX,'\n');
	for ( int i = 0; i < n; ++i ) {
		ResEntry r;
		in >> r;
		e.vResEntries.push_back(r);
	}
	return in;
} // operator>>

ostream& operator<<(ostream& out, const Entry& e) {
	int n = e.vResEntries.size();
	out << n << "\n";
	for ( int i = 0; i < n; ++i ) {
		out << e.vResEntries[i];
	} // i
	return out;
} // operator<<

// ____________________ File ____________________

File::File() {
} // File::File

File::File(const string& filename) {
	ifstream fin;
	open(fin, basic::database::full_name(filename) );
	//basic::database::open(fin,filename);
	fin >> *this;
} // File::File

istream& operator>>(istream& in, File& f) {
	f.clear();
	while ( in.peek() != EOF ) {
		Entry e;
		in >> e;
		f.addEntry(e);

	}
	return in;
} // operator>>

ostream& operator<<(ostream& out, const File& f) {
	for ( auto const & j : f.getEntries() ) {
		const vector<Entry>& vEntries = j.second;
		for ( auto const & vEntrie : vEntries ) {
			out << vEntrie;
		} // i
	} // j
	return out;
} // operator<<

void File::clear() {
	mEntries.clear();
} // File::clear

void File::addEntry(const Entry& e) {
	mEntries[e.vResEntries.size()].push_back(e);
} // File::addEntry

const Entry& File::getEntry(const int len) const {

	auto i = mEntries.find(len);

	debug_assert( i != mEntries.end() );

	const vector<Entry>& vEntries = i->second;
	return vEntries[ numeric::random::random_range(0,vEntries.size()-1) ];

} // File::getEntry

//       void File::filter(const double max_rama, const int frag_len) {

// static double dEdphi,dEdpsi;

// vector<Entry>& vEntries = mEntries[frag_len];

// vector<Entry> vFiltered;

// const Erg::RamaInfo::Reader* reader = Librarian::getErgRamaInfoReader();

// FORVC(i,Entry,vEntries) {

//   const Entry& e = *i;

//   bool good = true;

//   FORVC(j,ResEntry,e.vResEntries) {
//     const ResEntry& e = *j;

//     double E_rama = reader->getEnergy(ResType::ALA.getAa(),Erg::RamaInfo::Reader::SS_L,e.phi,e.psi,dEdphi,dEdpsi);
//     //cout << REPORT3(RAD2DEG*e.phi,RAD2DEG*e.psi,E_rama) << endl;

//     if( E_rama > max_rama ) {
//       good = false;
//     }

//   } // j

//   if( good ) {
//     vFiltered.push_back(e);
//   }

// } // i

// cout << REPORT2(vEntries.size(),vFiltered.size()) << endl;
// vEntries = vFiltered;

//       } // File::filter

//       void File::report(ostream& out) const {
// // FORMC(i,int,vector<Entry>,mEntries) {
// //   out << "Fragment::File: " << i->first << ": " << i->second.size() << "\n";
// // } // i
//       } // File::report

//       void File::reportRama(ostream& out) const {

// // FORMC(i,int,vector<Entry>,mEntries) {
// //   const vector<Entry>& vEntries = i->second;

// //   FORVC(j,Entry,vEntries) {
// //     const Entry& e = *j;
// //     FORVC(k,ResEntry,e.vResEntries) {

// //       (*k).reportRama(out);

// //     } // k

// //   } // j

// // } // i

//       } // File::reportRama

//       void File::reportRama(const string& outfile) const {
// // ofstream fout;
// // open(fout,outfile);
// // reportRama(fout);
// // fout.close();
//       } // File::reportRama


const File File::getFile(const int frag_len) const {

	File rval;
	rval.mEntries[frag_len] = find_or_throw(mEntries,frag_len);
	return rval;

} // File::getFile

void File::convertEntries(const int from, const int to) {

	debug_assert( to < from );

	auto i = mEntries.find(from);

	if ( i == mEntries.end() ) {
		return;
	}

	const vector<Entry>& vEntries_from = i->second;
	vector<Entry>& vEntries_to = mEntries[to];

	// don't need to clear vEntries_to

	//FORVC(k,Entry,vEntries_from) {
	for ( auto const & entry_from : vEntries_from ) {

		debug_assert( static_cast<int>(entry_from.vResEntries.size()) == from );

		for ( int i = 0; i < (from - to); i += to ) {
			Entry entry_to;
			for ( int j = 0; j < to; ++j ) {
				debug_assert( i+j < static_cast<int>(entry_from.vResEntries.size()) );
				entry_to.vResEntries.push_back( entry_from.vResEntries[i+j] );
			} // j
			debug_assert( static_cast<int>(entry_to.vResEntries.size()) == to );
			vEntries_to.push_back(entry_to);
		} // i

	} // k

} // File::addFragLen


}

// ===========================================================
// ==================== get fragment file ====================
// ===========================================================

Fragment::FileCOP Librarian::getFragmentFile_loop() {
	if ( frag_file_loop == nullptr ) {
		frag_file_loop = Fragment::FileOP( new Fragment::File(get("LLL")) );
		frag_file_loop->convertEntries(3,1);
	}
	return frag_file_loop;
} // Librarian::getFragmentFile_loop

Fragment::FileCOP Librarian::getFragmentFile_sheet() {
	if ( frag_file_sheet == nullptr ) {
		frag_file_sheet = Fragment::FileOP( new Fragment::File(get("EEE")) );
		frag_file_sheet->convertEntries(3,1);
	}
	return frag_file_sheet;
} // Librarian::getFragmentFile_sheet

Fragment::FileCOP Librarian::getFragmentFile_helix() {
	if ( frag_file_helix == nullptr ) {
		frag_file_helix = Fragment::FileOP( new Fragment::File(get("HHH")) );
		frag_file_helix->convertEntries(3,1);
	}
	return frag_file_helix;
} // Librarian::getFragmentFile_helix

Fragment::FileCOP Librarian::getFragmentFile(const char ss) {
	switch( ss ) {
	case 'H' :
		return getFragmentFile_helix();
	case 'E' :
		return getFragmentFile_sheet();
	case 'L':
	default :
		return getFragmentFile_loop();

	} // switch
} // Librarian::getFragmentFile

map<string,Fragment::FileOP> Librarian::mFragfiles_ss;

Fragment::FileOP Librarian::frag_file_loop;
Fragment::FileOP Librarian::frag_file_sheet;
Fragment::FileOP Librarian::frag_file_helix;

string Librarian::get( const string& s) {

	return "sampling/ss_fragfiles/" + s + ".fragfile";

}

namespace {

const vector<string> get_acceptable_ss(const string& ss0) {
	vector<string> rval;
	rval.push_back(ss0);

	for ( size_t i = 0; i < ss0.size(); ++i ) {
		if ( ss0[i] != 'L' ) {
			string ss = ss0;
			ss[i] = 'L';
			rval.push_back(ss);
		}
	} // i

	return rval;
} // get_acceptable_ss


} // namespace

Fragment::FileCOP Librarian::getFragmentFile(const string& ss0) {

	auto iter = mFragfiles_ss.find(ss0);

	if ( iter == mFragfiles_ss.end() ) {
		// try to read the fragment file from disc

		const vector<string> acceptable_ss = get_acceptable_ss(ss0); // XXX not correct to have to read in each one each time

		//       const string& dir    = get("fragment_file_dir");
		//       const string& suffix = get("fragment_file_suffix");

		bool found = false;

		Fragment::FileOP fragfile;

		for ( auto i = acceptable_ss.begin(); i != acceptable_ss.end() && !found ; ++i ) {

			const string& ss = i->substr(0,3);

			// build the filename
			string filename = get(ss);
			// filename.append(dir);
			// if( filename[filename.size()-1] != '/' ) {
			//   filename.append("/");
			// }
			// filename.append(ss);
			// if( suffix[0] != '.' ) {
			//   filename.append(".");
			// }
			// filename.append(suffix);

			// test to see whether said file exists
			// ifstream fin;
			// fin.open(filename.c_str());

			ifstream fin;

			//basic::database::open(fin,filename);

			open(fin,basic::database::full_name(filename));

			cout << "Fragment::getFragmentFile - trying to open '" << filename << "'" << endl;

			if ( fin ) {
				fragfile = Fragment::FileOP( new Fragment::File() );
				cout << "Fragment::getFragmentFile - successfully opened '" << filename << "'" << endl;
				fin >> *fragfile;
				found = true;
			}

		} // i

		debug_assert( fragfile != nullptr );

		mFragfiles_ss[ss0] = fragfile; // !!! NB: this is only going into the map for key ss0 (the original ss)
		iter = mFragfiles_ss.find(ss0);

	}

	debug_assert( iter != mFragfiles_ss.end() );

	return iter->second;

} // Librarian::getFragmentFile

} // namespace LoopDesign

} // namespace devel
