// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   devel/InvKinLigLoopDesign/Fragment.hh
///
/// @brief
/// @author

#ifndef DEVEL_INVKINLIGLOOPDESIGN_FRAGMENT_HH
#define DEVEL_INVKINLIGLOOPDESIGN_FRAGMENT_HH

#include <devel/inv_kin_lig_loop_design/Fragment.fwd.hh>

#include <iostream>
#include <vector>
#include <map>

namespace devel {

namespace inv_kin_lig_loop_design {

namespace Fragment {

// ____________________ Entry ____________________

class ResEntry {
public:
	std::string pdbid;
	char aa, ss;
	float phi,psi,ohm;

	void reportRama(std::ostream& out) const;

}; // class ResEntry

std::istream& operator>>(std::istream& in, ResEntry& e);
std::ostream& operator<<(std::ostream& in, const ResEntry& e);

// ____________________ Entry ____________________

class Entry {
public:
	std::vector<ResEntry> vResEntries;
}; // class Entry

std::istream& operator>>(std::istream& in, Entry& e);
std::ostream& operator<<(std::ostream& in, const Entry& e);

// ____________________ File ____________________

class File {
private:
	std::map<int,std::vector<Entry> > mEntries;

public:

	File();
	File(const std::string& file);

	void clear();
	void addEntry(const Entry& e);
	const Entry& getEntry(const int len) const;
	const std::map<int,std::vector<Entry> >& getEntries() const { return mEntries; }

	// XXX would like a filter method in File

	void convertEntries(const int frag_len_from, const int frag_len_to);

	void filter(const double max_rama, const int len);
	void report(std::ostream& out) const;

	const File getFile(const int frag_len) const;

	void reportRama(std::ostream& out) const;
	void reportRama(const std::string& outfile) const;

}; // class File

std::istream& operator>>(std::istream& in, File& f);
std::ostream& operator<<(std::ostream& in, const File& f);

} // namespace Fragment

struct Librarian {

	static Fragment::FileCOP getFragmentFile_loop();
	static Fragment::FileCOP getFragmentFile_sheet();
	static Fragment::FileCOP getFragmentFile_helix();
	static Fragment::FileCOP getFragmentFile(const char ss);
	static Fragment::FileCOP getFragmentFile(const std::string& ss0);

private:
	static std::string get(std::string const& name);

	static std::map<std::string,Fragment::FileOP> mFragfiles_ss;

	static Fragment::FileOP frag_file_loop;
	static Fragment::FileOP frag_file_sheet;
	static Fragment::FileOP frag_file_helix;

};


} // namespace LoopDesign

} // namespace Util

#endif // UTIL_FRAGMENT_HH
