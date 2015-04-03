// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/InvKinLigLoopDesign/Fragment.hh
///
/// @brief
/// @author

#ifndef DEVEL_INVKINLIGLOOPDESIGN_FRAGMENT_HH
#define DEVEL_INVKINLIGLOOPDESIGN_FRAGMENT_HH

#include <iostream>
#include <vector>
#include <map>

// you cannot #include yourself #include <devel/InvKinLigLoopDesign/Fragment.hh>

namespace devel {

  namespace inv_kin_lig_loop_design {

    using namespace std;

    namespace Fragment {

      // ____________________ Entry ____________________

      class ResEntry {
      public:
				string pdbid;
				char aa, ss;
				float phi,psi,ohm;

				void reportRama(ostream& out) const;

      }; // class ResEntry

      istream& operator>>(istream& in, ResEntry& e);
      ostream& operator<<(ostream& in, const ResEntry& e);

      // ____________________ Entry ____________________

      class Entry {
      public:
				vector<ResEntry> vResEntries;
      }; // class Entry

      istream& operator>>(istream& in, Entry& e);
      ostream& operator<<(ostream& in, const Entry& e);

      // ____________________ File ____________________

      class File {
      private:
				map<int,vector<Entry> > mEntries;

      public:

				File();
				File(const string& file);

				void clear();
				void addEntry(const Entry& e);
				const Entry& getEntry(const int len) const;
				const map<int,vector<Entry> >& getEntries() const { return mEntries; }

				// XXX would like a filter method in File

				void convertEntries(const int frag_len_from, const int frag_len_to);

				void filter(const double max_rama, const int len);
				void report(ostream& out) const;

				const File getFile(const int frag_len) const;

				void reportRama(ostream& out) const;
				void reportRama(const string& outfile) const;

      }; // class File

      istream& operator>>(istream& in, File& f);
      ostream& operator<<(ostream& in, const File& f);

    } // namespace Fragment

    struct Librarian {

      static const Fragment::File* getFragmentFile_loop();
      static const Fragment::File* getFragmentFile_sheet();
      static const Fragment::File* getFragmentFile_helix();
      static const Fragment::File* getFragmentFile(const char ss);
      static const Fragment::File* getFragmentFile(const string& ss0);

    private:
      static string get(string const& name);

      static map<string,Fragment::File*> mFragfiles_ss;

      static Fragment::File* frag_file_loop;
      static Fragment::File* frag_file_sheet;
      static Fragment::File* frag_file_helix;

    };


  } // namespace LoopDesign

} // namespace Util

#endif // UTIL_FRAGMENT_HH
