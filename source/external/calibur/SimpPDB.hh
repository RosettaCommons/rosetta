// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/kalngyk/SimpPDB.hh
/// @author SC Li & YK Ng (kalngyk@gmail.com)

#ifndef apps_pilot_kalngyk_SimpPDB_HH
#define apps_pilot_kalngyk_SimpPDB_HH

#define LONGEST_CHAIN 4000

#include <iostream>
#include <fstream>

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>

#include <map>
#include <vector>

#include <calibur/SimpPDB.fwd.hh>

using namespace std;
using namespace apps::pilot::kalngyk;

enum INPUT_FILE_TYPE { UNKNOWN=-1, SILENT_FILE, PDB_LIST };
INPUT_FILE_TYPE filetype(char * filename);
unsigned int num_lines_in_file(char * filename);

//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

char aa[][4] = {"BCK","GLY","ALA","SER","CYS","VAL","THR","ILE",
              "PRO","MET","ASP","ASN","LEU",
              "LYS","GLU","GLN","ARG",
              "HIS","PHE","TYR","TRP","CYX", "MSE"};

char slc[] = {'X','G','A','S','C','V','T','I',
    	      'P','M','D','N','L','K','E','Q','R',
              'H','F','Y','W','C', 'm'};


int toInt(const string& aString)
{
    static char st[20];
    int start = 0;
    for (unsigned int i=0; i < aString.size(); i++)
       if (aString[i] != ' ')
           st[start++] = aString[i];
    st[start] = '\0';
    int rev = atoi(st);
    return rev;
}


double toFloat(const string& aString)
{
    static char st[20];
    int start = 0;
    for (unsigned int i=0; i < aString.size(); i++)
        if (aString[i] != ' ')
            st[start++] = aString[i];
    st[start] = '\0';
    double rev = atof(st);
    return rev;
}


void center_residues(double * mCAlpha, int mNumResidue)
{
    double cx = 0;
    double cy = 0;
    double cz = 0;

    int i3 = 0;
    for (int i=0; i < mNumResidue; i++)
    {
        cx += mCAlpha[i3];
        cy += mCAlpha[i3+1];
        cz += mCAlpha[i3+2];
        i3 += 3;
    }

    cx /= mNumResidue;
    cy /= mNumResidue;
    cz /= mNumResidue;
    i3 = 0;
    for (int i=0; i < mNumResidue; i++)
    {
        mCAlpha[i3]   -= cx;
        mCAlpha[i3+1] -= cy;
        mCAlpha[i3+2] -= cz;
        i3 += 3;
    }
}


/**
 * This class enables two things:
 * 1. Caching of PDB file content in memory, hence reducing disk access
 * 2. Provides support for silentfile
 *
 * To use it:
 * 1. Create an instance,
 * 2. Load it up either with a silent file or with a list of PDB files,
 * 3. Set SimPDB::preloadedPDB to it, and set SimPDB::preloadPDB to true.
 *
 * The use of PreloadedPDB is compulsory for silent files.
 *
 * For a list of decoys, PreloadedPDB is used by default.
 * However, this is bad in some situations. e.g. 36,000 decoys of 100
 * residues each will require about 42M memory. In such situations,
 * using the disk might be better. Hence, PreloadedPDB is off when:
 * 1. Users override it (by switching off SimPDB::preloadPDB)
 * 2. The number of decoys used is more than ADVISED_THRESHOLD.
 */
class PreloadedPDB
{
public:
    static unsigned int ADVISED_THRESHOLD;

private:
    char * silentfilename;
    char * pdblistfilename;

public:
    int mNumResidue;
    int mNumDecoy;
    vector<char *> * mNames;
    map<char *, SimPDB *> filename2PDB; // fix this if it is deemed too slow


public:
    PreloadedPDB();
    ~PreloadedPDB();
    void loadSilentFile(char * silentfilename);
    void loadPDBFromList(char * pdblistfilename);

    SimPDB * getSimPDB(char * pdbfilename); // unused. for internal testing
};


namespace apps { namespace pilot { namespace kalngyk {
class SimPDB
{
    public:
      /**
       * These parameters control how PDB files are to be loaded.
       * They do not apply to silent file, which are assumed to be
       * pre-processed.
       */
      static int s_residue;
      static int e_residue;
      static char * chains;

      /**
       * This feature allows the preloading of SimPDB objects.
       * SimPDB will then be obtained from preloadedPDB instead of from disk.
       */
      static PreloadedPDB * preloadedPDB;
      static bool preloadPDB;

    public:
      const char* mProteinFileName;
      int mNumResidue;
      //double mSquaredSum;
      double * mCAlpha;
      void read();

    public:
      SimPDB(char* aProteinFileName);
      SimPDB(char* aProteinFileName, int len);
      ~SimPDB();

      // Special constructors used only by PreloadedPDB. Don't touch.
      SimPDB();
      SimPDB(int len);
};
}}}


/**
 * Reasons for choosing 36,000 as the threshold for using disk access
 * 1. 36,000 decoys of 100 residues each is about 44MB, which is still
 *    acceptable for a workstation PC. If user has less capabled hardware,
 *    the "-d" switch may be used.
 * 2. For more than 36,000 decoys, the bottleneck is probably not in disk.
 */
unsigned int PreloadedPDB::ADVISED_THRESHOLD = 36000;


INPUT_FILE_TYPE
filetype(char * filename)
{
    char buf[400];
    ifstream input(filename);
    string line;
    if (!input)
    {
        cerr << "Can't open input file \"" << filename << "\"" << endl;
        exit(0);
    }

    input.getline(buf, 400);
    line=buf;
    if (line.substr(0, 9)=="SEQUENCE:")
    {
        input.close();
        return SILENT_FILE;
    }

    input.seekg(0);
    input.getline(buf, 400);
    char *token = strtok(buf, " ");
    if (token == NULL)
    {
        input.close();
        return UNKNOWN;
    }
    char *name = new char[strlen(token)+1];
    strcpy(name, token);

    ifstream pdbfile(name); // check if name is a file
    if (pdbfile)
    {
        pdbfile.close();
        input.close();
        return PDB_LIST;
    }
    pdbfile.close();

    input.close();

    return UNKNOWN;
}


unsigned int num_lines_in_file(char * filename)
{
    char buf[400];
    ifstream input(filename);
    if (!input)
    {
        cerr << "Can't open input file \"" << filename << "\"" << endl;
        exit(0);
    }
    unsigned int count = 0;
    while (!input.eof())
    {
        input.getline(buf, 400);
        count++;
    }
    return count-1;
}


//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

PreloadedPDB::PreloadedPDB()
{
}


/**
 * Populate the PreloadedPDB with a silentfile
 */
void
PreloadedPDB::loadSilentFile(char * filename)
{
    char buf[400];
    ifstream input(filename);
    string line;
    if (!input)
    {
        cerr << "Can't open silent file \"" << filename << "\"" << endl;
        exit(0);
    }

    silentfilename = filename;
    pdblistfilename = NULL;

    /**
     * Determine the number of residues from the silent file
     */
    input.getline(buf, 400);
    line=buf;
    if (line.substr(0, 9)!="SEQUENCE:")
    {
        cerr << "Silent file \"" << filename << "\" began with:" << endl
             << "    " << line << endl;
        exit(0);
    }

    input.getline(buf, 400);
    line=buf;
    if (line.substr(0, 6)!="SCORE:")
    {
        cerr << "Silent file \"" << filename << "\" began with:" << endl
             << "    " << line << endl;
        exit(0);
    }

    input.getline(buf, 400);
    line=buf;
    if (line.substr(0, 6)!="SCORE:")
    {
        cerr << "Silent file \"" << filename << "\" began with:" << endl
             << "    " << line << endl;
        exit(0);
    }

    int residueID;
    for (int i=1; !input.eof(); i++)
    {
        input.getline(buf, 400);
        line = buf;
        if (line.substr(0, 6) == "SCORE:")
            break;
        residueID = toInt(line.substr(1, 4));
        if (residueID != i)
        {
            cerr << "Residue id " << residueID
                 << " out of sequence in silent file" << endl;
            exit(0);
        }
    }

    mNumResidue = residueID;
    input.seekg(0);
    input.getline(buf, 400);
    input.getline(buf, 400);
    input.getline(buf, 400);

    /**
     * Read PDBs into filename2PDB
     */
    SimPDB * pdb = new SimPDB(mNumResidue);
    bool isNewPDB = true;
    int numResidue = 0;
    int decoyCount = 1;
    char * key;
    do
    {
        input.getline(buf, 400);
        line = buf;

        if (line.substr(0, 6) == "SCORE:") // Old PDB done
        {
            // Check if PDB has mNumResidue
            if (numResidue != mNumResidue)
            {
                cerr << "Insufficient residues in the " << decoyCount
                     << "-th decoy in silent file" << endl;
                exit(0);
            }

            // Insert the pdb
            filename2PDB[key] = pdb;
            pdb->mProteinFileName = key;
            center_residues(pdb->mCAlpha, pdb->mNumResidue);

            // Start a new pdb
            pdb = new SimPDB(mNumResidue);
            isNewPDB = true;
            numResidue = 0;

            decoyCount++;
        }
        else if (line != "") // Sometimes an empty string is read at eof
        {
            if (isNewPDB)
            {
                // Get the filename
                string filename = line.substr(62, 400);
                string ext = filename.substr(filename.length()-4, 4);

                // If filename does not end in .pdb, generate a filename
                if (strcmp(ext.c_str(), ".pdb"))
                    key = strdup("decoy" + decoyCount);
                else
                    key = strdup(filename.c_str());

                isNewPDB = false;
            }

            residueID = toInt(line.substr( 0, 4));
            double x = toFloat(line.substr(35, 8));
            double y = toFloat(line.substr(44, 9));
            double z = toFloat(line.substr(54, 8));
            pdb->mCAlpha[numResidue*3]   = x;
            pdb->mCAlpha[numResidue*3+1] = y;
            pdb->mCAlpha[numResidue*3+2] = z;

            if (residueID != numResidue+1)
            {
                cout << residueID << "," << numResidue << endl;
                cerr << "Residue ID out of sequence in the " << decoyCount
                     << "-th decoy in silent file" << endl;
                exit(0);
            }

            numResidue++;
        }

    }
    while (!input.eof());
    input.close();

    // Insert the final pdb
    filename2PDB[key] = pdb;
    pdb->mProteinFileName = key;

    mNumDecoy = decoyCount;

    /**
     * Generate the vector of filenames.
     * We can insert filename at the same time as inserting the SimPDB into
     * filename2PDB. However, doing this here has the advantage that if the
     * same filename is entered twice into filename2PDB, the name will not
     * be duplicated here.
     */
    mNames = new vector<char *>(0);
    map<char *,SimPDB *>::iterator it = filename2PDB.begin();
    for (int i=0; it != filename2PDB.end(); it++, i++)
        mNames->push_back( it->first );
}


/**
 * Populate the PreloadedPDB with the PDB files specified in a list.
 *
 * Reading of the PDB files is through SimPDB.
 */
void
PreloadedPDB::loadPDBFromList(char * filename)
{
    ifstream input(filename);
    if (!input)
    {
        cerr << "Cannot find file \"" << filename << "\"" << endl;
        exit(0);
    }

    char buf[400];
    char* token;

    silentfilename = NULL;
    pdblistfilename = filename;

    /**
     * Read in names of PDB files
     */
    mNames = new vector<char *>(0);
    while (!input.eof())
    {
        input.getline(buf, 400);
        token = strtok(buf, " ");
        if (token == NULL)
            continue;
        char * name = new char[strlen(token)+1];
        strcpy(name, token);
        mNames->push_back(name);
    }
    input.close();

    mNumDecoy = mNames->size();

    SimPDB * pdb = new SimPDB();
    pdb->mProteinFileName = (*mNames)[0];
    pdb->mNumResidue = LONGEST_CHAIN;
    pdb->mCAlpha = new double[3*LONGEST_CHAIN];
    pdb->read();

    mNumResidue = pdb->mNumResidue;

    filename2PDB[(*mNames)[0]] = pdb;

    for (unsigned int i=1; i < mNames->size(); i++)
    {
        SimPDB * pdb = new SimPDB(mNumResidue);
        pdb->mProteinFileName = (*mNames)[i];
        pdb->read();
        filename2PDB[(*mNames)[i]] = pdb;
    }
}


SimPDB *
PreloadedPDB::getSimPDB(char * filename)
{
    return filename2PDB[filename];
}





//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -
// Memory and default values for the static variables

// Modifies how SimPDB reads PDB files
int SimPDB::s_residue = 1;
int SimPDB::e_residue = LONGEST_CHAIN;
char * SimPDB::chains = strdup("AC ");

// Set these two fields to tell SimPDB to use the PreloadedPDB mechanism
PreloadedPDB * SimPDB::preloadedPDB = NULL;
bool SimPDB::preloadPDB = true;

//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -
// General purpose constructors which handle the preloadedPDB mechanism

SimPDB::SimPDB(char *aFileName)
{
    if (preloadPDB)
    {
        /*for(int i=0; i<len; i++)
        {
          mCAlpha[i]=new double [3];
        }*/
        SimPDB * pdb = preloadedPDB->filename2PDB[aFileName];
        mProteinFileName = strdup(pdb->mProteinFileName);
        mNumResidue = pdb->mNumResidue;
        mCAlpha = new double[3*mNumResidue];
        memcpy(mCAlpha, pdb->mCAlpha, 3 * mNumResidue * sizeof(double));
    }
    else
    {
        mProteinFileName = aFileName;
        mNumResidue = LONGEST_CHAIN;
        mCAlpha = new double[3*LONGEST_CHAIN];
        read();
    }
}

SimPDB::SimPDB(char *aFileName, int len)
{
    if (preloadPDB)
    {
        /*for(int i=0; i<len; i++)
        {
          mCAlpha[i]=new double [3];
        }*/
        SimPDB * pdb = preloadedPDB->filename2PDB[aFileName];
        mProteinFileName = strdup(pdb->mProteinFileName);
        mNumResidue = pdb->mNumResidue;
        mCAlpha = new double[3*mNumResidue];
        memcpy(mCAlpha, pdb->mCAlpha, 3 * mNumResidue * sizeof(double));
    }
    else
    {
        mProteinFileName = aFileName;
        mNumResidue = len;
        mCAlpha = new double[3*len];
        read();
    }
}


SimPDB::~SimPDB()
{
    /* for(int i=0; i<mNumResidue; i++)
    {
      delete [] mCAlpha[i];
    }*/
    delete [] mCAlpha;
}


//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -
// Constructors which DO NOT handle the preloadedPDB mechanism
//
// They should be called from PreloadedPDB only, since it will need to
// bypass the mechanism

SimPDB::SimPDB() {}

SimPDB::SimPDB(int len)
{
    mNumResidue = len;
    mCAlpha = new double[3*len];
}



//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

/**
 * Reads a PDB file from disk. Do not read from PDB files anywhere else.
 */
void
SimPDB::read()
{
    ifstream input(mProteinFileName);
    if (!input)
    {
        cerr << "Cannot find protein file " << mProteinFileName << endl;
        exit(0);
    }
    cout.flush();
    char buf[400];
    mNumResidue = 0;
    double x, y, z;
    //char c = 'a';
    bool read = false;

    int prevID = -10000;
    int count = 0;

    int CA_number = 1;

    //double cx = 0;
    //double cy = 0;
    //double cz = 0;
    int count3 = 0;
    //mSquaredSum=0;
    while (!input.eof())
    {
        input.getline(buf, 400);
        string line=buf;
        if (line.substr(0, 3) == "TER" && read == true) break;
        if (line.substr(0, 6) == "ENDMDL") break;

        if (line.substr(0, 4) != "ATOM" && line.substr(0, 6) != "HETATM")
            continue;			

        if (line.substr(13, 4) == "CA  " || line.substr(13, 4) == " CA "
            || line.substr(13, 4) == "  CA" || line.substr(13, 2) == "CA")
        {
            // At this point a CA atom has been discovered
            // We want to further filter it based on the following two
            // criteria: chain and region.

            // Check if the chain which this atom belongs to is to be included
            bool include_chain = false;
            for (char * c = SimPDB::chains; *c; c++)
                if (toupper(line[21]) == toupper(*c))
                {
                    include_chain = true;
                    break;
                }
            if (!include_chain)
            {
                CA_number++;
                continue;
            }

            // Check if the CA atom is within the region to analyze
            if (CA_number < SimPDB::s_residue)
            {
                CA_number++;
                continue;
            }
            else if (CA_number > SimPDB::e_residue)
                break;

            read = true;
            int residueID = toInt(line.substr(22, 6));
            if (residueID == prevID)
                continue;

            prevID = residueID;
            x = toFloat(line.substr(30, 8));
            y = toFloat(line.substr(38, 8));
            z = toFloat(line.substr(46, 8));
            string AAType = line.substr();
            count3 = 3*count;
            mCAlpha[count3]   = x;
            mCAlpha[count3+1] = y;
            mCAlpha[count3+2] = z;
            //mSquaredSum+=x*x+y*y+z*z;
            count++;

            CA_number++;
        }
    }//while

    mNumResidue = count;
    input.close();

    center_residues(mCAlpha, mNumResidue);
}

/*
int main(int argc, char** argv)
{
    SimPDB::chains = strdup("ACH");
    SimPDB::s_residue = 215;
    SimPDB::e_residue = 224;
    SimPDB* sim = new SimPDB(argv[1]);

    cout << "numResidue=" << sim->mNumResidue << endl;
    for (int i=0; i < sim->mNumResidue; i++)
    {
        cout << sim->mCAlpha[i*3] << ","
             << sim->mCAlpha[i*3+1] << ","
             << sim->mCAlpha[i*3+2] << endl;
    }
}*/

/*
int main()
{
    //char * file = "rosetta/silent_file";
    char * file = "list";

    // Preload the PDBs, either from silent file or from list of decoys
    PreloadedPDB * pdbs = new PreloadedPDB();
    switch (filetype(file))
    {
        case SILENT_FILE: pdbs->loadSilentFile(file); break;
        case PDB_LIST: pdbs->loadPDBFromList(file); break;
        default: cerr << "Unknown file type" << endl; exit(0);
    }

    // Attach the cache to SimPDB
    SimPDB::preloadedPDB = pdbs;
    SimPDB::preloadPDB = true;

    // Get the names from the loaded PDB for testing
    vector<char *> * mNames = pdbs->mNames;

    // Get the PDB for each name
    for (int i=0; i < mNames->size(); i++)
    {
        char * filename = (*mNames)[i];
        SimPDB * pdb = new SimPDB(filename);
        cout << filename << ":" << endl;
        for (int j=0; j < pdb->mNumResidue; j++)
            cout << "    " << pdb->mCAlpha[j*3] << ", "
                           << pdb->mCAlpha[j*3+1] << ", "
                           << pdb->mCAlpha[j*3+2] << endl;
    }
}
*/



#endif
