// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/kalngyk/InitCluster.hh
/// @author SC Li & YK Ng (kalngyk@gmail.com)

#ifndef apps_pilot_kalngyk_InitCluster_HH
#define apps_pilot_kalngyk_InitCluster_HH

#include <vector>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <assert.h>
#include <time.h>
#ifndef __WIN32__
#include <sys/resource.h>
#endif

#include <calibur/SimpPDB.hh>
#include <calibur/rmsd.hh>


using namespace std;

#ifndef _LARGE_DECOY_SET_
typedef unsigned short LIST_TYPE;
#else  //number of input decoys large than 65535
typedef unsigned int LIST_TYPE;
#endif

#define _SHOW_PERCENTAGE_COMPLETE_

// _MATRIX_MODE_LIMIT_ is determined by the RAM of the system
#define _MATRIX_MODE_LIMIT_ 13000

#define _OVER_RMSD_ 9999999

#define REFERENCE_SIZE 6
#define RANDOM_DECOY_SIZE_FOR_FILTERING 101
#define RANDOM_DECOY_SIZE_FOR_THRESHOLD 101
#define NUM_TRIALS_FOR_THRESHOLD 16
#define DEFAULT_PERCENTILE_FOR_THRESHOLD 10
#define MAX_PERCENTILE_FOR_THRESHOLD 50
#define MIN_PERCENTILE_FOR_THRESHOLD 3

enum ADJ_LIST_MODE { MATRIX, LIST , LITE };

enum EST_THRESHOLD_MODE {
    PERCENT_EDGES,      // % quantile pairwise distance (default)
    MIN_AVG_DIST_BASED, // t = a * (min avg dist) + b
    MOST_FREQ_BASED,    // use (2/3)*(most frequently occuring distance)
    ROSETTA,            // ROSETTA's method
    SAMPLED_ROSETTA,    // ROSETTA's method using sampled decoys
    USER_SPECIFIED      // user supplied
};


class Pref
{
public:
    static bool use_scud;
    static bool use_sig;
    static bool ref_use_first_few_decoys;
};


class AdjacentList
{
public:
    static ADJ_LIST_MODE mListMode;
    int mWhich;         // index of the decoy this AdjacentList is for
    int mNumNeigh;      // synchronized with the size of neigh
    vector<int>* neigh; // keep a record of all the neighbors
    vector<double>* dist;
    LIST_TYPE* mReverseIndex; // References index of the decoy within the array
    AdjacentList();
    AdjacentList(int which, int num);
    ~AdjacentList();
    void add(int n, double d, bool ifNeigh);
#ifdef _ADD_LITE_MODE_
    void add(int num);
#endif
    void remove(int n);
    double getD(int n);
};


class Stru
{
public:
    double* mCAlpha;
    double* mSIG; //signature
    SimPDB* mPDB;
    Stru(SimPDB* pdb, int len);
    ~Stru();
    double dist(double x, double y, double z, double *zz);
    double dist(double x, double y, double z);
};


class Clustering
{
public:
    static EST_THRESHOLD_MODE EST_THRESHOLD;
    static bool FILTER_MODE;
    static double xPercentile;
    static bool autoAdjustPercentile;
    static double xFactor;

    char* mInputFileName;   // file which contains all PDB filenames
    vector<char* >* mNames; // all decoy (file) names
    vector<Stru* >* mPDBs;  // all decoy PDBs
    int mNumPDB;            // will be set to mPDBs->size()
    int mLen;               // #residues

    double THRESHOLD;        // clustering threshold. most important parameter

    // - = - = - = - = - = - = - = - = - = - = - = - = - = -
    // for auxiliary grouping

    double CLU_RADIUS;     // cluster radius for auxClustering()
    vector<int>* mCluCen; // cluster centers found using auxClustering()
    vector<int>** mAuxCluster; // cluster elements
    int* mCen;
    double* mD2C;          // distance from decoy in auxCluster to CluCen
    int * mNumNeighbor;   // the number of neighbors of each decoy
    double bestClusMargin; // size(bestClus) -origsuze(2ndClus) /size(bestClus)
    int bestClusSize;

    // - = - = - = - = - = - = - = - = - = - = - = - = - = -
    // for clustering

    AdjacentList** mAdjacentList; // lists of all neighbors
    double* mReference;      // for {lower,upper}bounds through references

    int mFinalDecoy;
    vector<AdjacentList *> *mFinalClusters;
    int * mRemainingList;
    int * mRemainingListIndex;
    int mRemainingSize;

    //- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -
    // METHODS
    //- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

    Clustering(); // do nothing
    void initialize(char * filename, double threshold);
    void reinitialize(vector<char *>*, vector<Stru *>*, double threshold);
    void cluster();
    void showClusters(int);
    void getPDBs(vector<char *>*, vector<Stru *>*, vector<int>*, int);


    //- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

    // for reading decoys from input files
    void readDecoyNames();
    void readDecoys(vector<char *>*, vector<Stru *>*);
    vector<Stru *>* readDecoys(vector<char *>*);
    //void refilterDecoys(vector<char *>*, vector<Stru *>*);

    // - = - = - = - = - = - = - = - = - = - = - = -
    // for finding threshold

    void getThresholdAndDecoys();
    double getThreshold(double **, int, int, int, int, double, double);
    double ** getNborList(vector<Stru *> *, vector<char *> *, int);
    //double ** getNborList(vector<Stru *> *, vector<Stru *> *, int);
    double ** getNborList(vector<Stru *> *, int, double *, double *);
    void estimateDist(vector<char *>*, int, int, double,
                      double *, double *, double *, double *);
    void estimateDist(vector<Stru *>*, double,
                      double *, double *, double *, double *);

    vector<char *>* getRandomDecoyNames(vector<char *>*, int, int);
    void destroyRandomDecoys(vector<char *>*, vector<Stru *>*);


    // - = - = - = - = - = - = - = - = - = - = - = -
    // for clustering

    void auxClustering();
    void buildAdjacentLists();
    void listAdjacentLists();

    void findLargestClusters();
    int findDecoyWithMostNeighbors();
    void removeDecoys(AdjacentList * adj);

    // - = - = - = - = - = - = - = - = - = - = - = -

    void initRef(int* index);
    void refBound(int i, int j, double& lower, double& upper);
    //bool find(int which, vector<int> *elements); // too slow

    // - = - = - = - = - = - = - = - = - = - = - = -

    void realignDecoys(int ref);
    void superimposeAndReplace(double* coor1, double* coor2);
    double eucD(int i, int j);

    // - = - = - = - = - = - = - = - = - = - = - = -

    // methods for rmsd computation
    double estD(int i, int j);
    double estD(Stru* a, Stru* b);
    double trueD(Stru* a, Stru* b);
    double trueD(int i, int j);
    // storage for RMSD() computation
    void allocateSpaceForRMSD(int len);
    bool spaceAllocatedForRMSD;
    double *result_coords;
    // storage for rmsfit_() computation
    double *coord1;
    double *coord2;

    // - = - = - = - = - = - = - = - = - = - = - = -

    //double timeval_difference(struct timeval * x, struct timeval * y);
    double get_elapsed(int restart);
};

//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -
// defaults for commandline changable values

// LIST, MATRIX, or LITE (in LITE mode, the AdjacentLists are actually empty)
#ifndef _ADD_LITE_MODE_
ADJ_LIST_MODE AdjacentList::mListMode = MATRIX;
#else
ADJ_LIST_MODE AdjacentList::mListMode = LITE;
#endif

bool Pref::use_scud = false;
bool Pref::use_sig = true;
bool Pref::ref_use_first_few_decoys = true;

EST_THRESHOLD_MODE Clustering::EST_THRESHOLD = PERCENT_EDGES;
bool Clustering::FILTER_MODE = true;
double Clustering::xPercentile = DEFAULT_PERCENTILE_FOR_THRESHOLD;
bool Clustering::autoAdjustPercentile = true;
double Clustering::xFactor = 2./3;


//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

#ifndef __WIN32__
static double
__timeval_difference(struct timeval * x, struct timeval * y)
{
    double elapsed;
    if (x->tv_usec < y->tv_usec)
    {
        int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
        y->tv_usec -= 1000000 * nsec;
        y->tv_sec += nsec;
    }
    if (x->tv_usec - y->tv_usec > 1000000)
    {
        int nsec = (x->tv_usec - y->tv_usec) / 1000000;
        y->tv_usec += 1000000 * nsec;
        y->tv_sec -= nsec;
    }

    elapsed = x->tv_sec - y->tv_sec;
    elapsed += (x->tv_usec - y->tv_usec)/(double)1000000;
    return elapsed;
}

static double
_get_elapsed(int set_start)
{
    static struct rusage last;
    struct rusage now;
    double elapsed = 0;
    if (set_start)
        getrusage(RUSAGE_SELF, &last);
    else
    {
        getrusage(RUSAGE_SELF, &now);
        elapsed += __timeval_difference(&(now.ru_utime), &(last.ru_utime));
        elapsed += __timeval_difference(&(now.ru_stime), &(last.ru_stime));
        last = now;
    }
    return elapsed;
}
#endif

//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

//============================= Adjacent List =================================
// We have two representation if AdjacentList;
// one is list, the other is matrix (arrays), the latter is used when the input
// size is small
// Each AdjacentList contains all the neighbors of a decoy.
// It is similar to the "neighbors" variable in cluster_info_silent.c.

// for LIST mode
AdjacentList::AdjacentList()
{
    mWhich = 0;
    neigh = new vector<int>(0);
    dist = new vector<double>(0);
    mReverseIndex = NULL; // not used
    mNumNeigh = 0;
}

// for MATRIX mode
AdjacentList::AdjacentList(int which, int size)
{
    mWhich = which;
    neigh = new vector<int>(size);   // max size. this will not grow
    dist = new vector<double>(size); // max size. this will not grow
    mReverseIndex = new LIST_TYPE [size]; // mPDB[i]'s index in neigh/dist
    memset(mReverseIndex, 0, size*sizeof(LIST_TYPE));
    for (int i=0; i < size; i++)
        (*dist)[i] = _OVER_RMSD_;
    mNumNeigh = 0;
}

AdjacentList::~AdjacentList()
{
    delete neigh;
    delete dist;
    if (mReverseIndex)
        delete [] mReverseIndex;
}

double
AdjacentList::getD(int which)
{
    if (mListMode == LIST
#ifdef _ADD_LITE_MODE_
           || mListMode == LITE
#endif
        )
        return _OVER_RMSD_;
    return (*dist)[which];
}

void
AdjacentList::add(int n, double d, bool ifNeigh)
{
#ifdef _ADD_LITE_MODE_
    if (mListMode == LITE) return;
#endif
    if (mListMode != LIST)
        (*dist)[n] = d;
    if (ifNeigh)
    {
        if (mListMode == LIST)
        {
            neigh->push_back(n);
            dist->push_back(d);
        }
        else
        {
            (*neigh)[mNumNeigh] = n;
            mReverseIndex[n] = mNumNeigh;
        }
        mNumNeigh++;
    }
}

// Used only in LITE mode
#ifdef _ADD_LITE_MODE_
void
AdjacentList::add(int num)
{
    if (mListMode == LITE)
        mNumNeigh += num;
}
#endif

void
AdjacentList::remove(int n)
{
#ifdef _ADD_LITE_MODE_
    if (mListMode == LITE) return;
#endif
   	//if (mNumNeigh < 1)
    //    return;
    if (mListMode == LIST)
    {
        int i;
        for (i=0; i < mNumNeigh; i++) // this is slow
           if ((*neigh)[i]==n)
               break;
        if (i >= mNumNeigh) // not found! should return error
            return;

        // move the last entry to this one
        (*neigh)[i] = (*neigh)[mNumNeigh-1];
        (*neigh).pop_back();
        (*dist)[i] = (*dist)[mNumNeigh-1];
        (*dist).pop_back();
        mNumNeigh--;
    }
    else
    {
        // we assume there is at least one element in each AdjacentList: itself
        LIST_TYPE index = mReverseIndex[n]; // danger! could be out of bound
        //cout<<"Removing "<<index<<"-th ("<<n<<") in "<<mWhich<<"'s LIST"<<endl;
        int last = (*neigh)[mNumNeigh-1];
        //cout<<" swap with " << (mNumNeigh-1) << " with id " << last << endl;
        (*neigh)[index] = last;
        mReverseIndex[last] = index;
        mNumNeigh--;
        if (mNumNeigh < 0)
        {
            cout << "OMG" << endl;
            exit(0);
        }
    }
}
//============================ Adjacent List ==================================

//============================== Structure ====================================
Stru::Stru(SimPDB* pdb, int len)
{
    mPDB = pdb;
    mCAlpha = pdb->mCAlpha;
    if (Pref::use_sig)
    {
    mSIG = new double [len];
    //double cx = 0;
    //double cy = 0;
    //double cz = 0;
    //compute the centroid
    /*for (int i=0; i < len; i++)
    {
        cx += mCAlpha[3*i];
        cy += mCAlpha[3*i+1];
        cz += mCAlpha[3*i+2];
    }*/
    //cx /= len;
    //cy /= len;
    //cz /= len;
    for (int i=0; i<len; i++)
    {
        //mSIG[i] = dist(cx, cy, cz, mCAlpha+3*i);
        // for the new SimpPDB.cc (it became slower so I reverted back)
        mSIG[i] = dist(mCAlpha[3*i], mCAlpha[3*i+1], mCAlpha[3*i+2]);
    }
    }
}

Stru::~Stru()
{
    if (Pref::use_sig)
        delete [] mSIG;
    delete mPDB;
}

double
Stru::dist(double x, double y, double z, double *zz)
{
    double xd=x-zz[0];
    double yd=y-zz[1];
    double zd=z-zz[2];
    return sqrt(xd*xd+yd*yd+zd*zd);
}

double
Stru::dist(double x, double y, double z)
{
    return sqrt(x*x+y*y+z*z);
}
//============================= Structure =====================================

//============================ Main object ====================================

Clustering::Clustering()
{
    mLen = 0;
    spaceAllocatedForRMSD = false;
    bestClusMargin = 1.; // should be a value that will not trigger re-cluster
#ifdef _SPICKER_SAMPLING_
    FILTER_MODE = false;
#endif
    cout.setf(ios::fixed);
    cout << setprecision(5);
}

/**
 * Main aim:
 * (1) read decoys into mNames and mPDBs.
 * (2) decide THRESHOLD, CLUS_RADIUS
 * (3) initialize data structure needed for auxiliary grouping, and
 * (4) ...for clustering.
 */
void
Clustering::initialize(char * filename, double threshold)
{
    mInputFileName = filename;
    THRESHOLD = threshold;

    // Read decoys and get threshold - = - = - = - = - = - = -

    if (THRESHOLD != -1)
        EST_THRESHOLD = USER_SPECIFIED;

    cout << "Filtering " << (FILTER_MODE? "on": "off") << endl;
    cout << "Signature mode " << (Pref::use_sig? "on": "off") << endl;
    cout << "Using chains ";
    for (char * c = SimPDB::chains; *c; c++)
        cout << "'" << *c << "'" << (*(c+1)=='\0'? "": ", ");
    cout << " in PDB files" << endl;

    getThresholdAndDecoys();

    CLU_RADIUS = THRESHOLD / 2.0 - 0.00001;

    cout << "Initialized " << mNumPDB << " decoys." << endl;

    // Realign decoys if needed - = - = - = - = - = - = - = -

    if (THRESHOLD > 2.5)
        Pref::use_scud = true;

    cout << "Use of SCUD's bound " << (Pref::use_scud? "on": "off") << endl;

    if (Pref::use_scud)
        realignDecoys(0);

    if (EST_THRESHOLD == USER_SPECIFIED)
        cout << "Using user specified threshold: " << THRESHOLD << endl;

    // Initialize auxiliary decoys grouping - = - = - = - = -

    // mCluCen contains the cluster centers. Each cluster center is referenced
    // by an integer x, which is the decoy's index in mPDBs.
    // Elements of the cluster is in mAuxCluster[x];
    // Since most decoys are going to end up not being cluster centers, for
    // most i, mAuxCluster[i] will be NULL.
    mCluCen = new vector<int>(0);
    mAuxCluster = new vector<int>*[mNumPDB];
    for (int i=0; i < mNumPDB; i++)
        mAuxCluster[i] = NULL;
    mD2C  = new double[mNumPDB];
    mCen  = new int[mNumPDB];

    // Initialize clustering - = - = - = - = - = - = - = - = -

    mAdjacentList = new AdjacentList* [mNumPDB];

#ifdef _ADD_LITE_MODE_
    if (AdjacentList::mListMode != LITE)
    {
#endif
        // auto-switch between LIST and MATRIX mode
        if (mNumPDB > _MATRIX_MODE_LIMIT_)
        {
            cout << "Using LIST mode" << endl;
            AdjacentList::mListMode = LIST;
        }
        else
        {
            cout << "Using MATRIX mode" << endl;
            AdjacentList::mListMode = MATRIX;
        }
#ifdef _ADD_LITE_MODE_
    }
    else
    {
        cout << "Using LITE mode" << endl;
        AdjacentList::mListMode = LITE;
    }
#endif

    for (int i=0; i < mNumPDB; i++)
    {
        if(AdjacentList::mListMode == MATRIX)
        {
            mAdjacentList[i] = new AdjacentList(i, mNumPDB);
        }
        else
        {
            mAdjacentList[i] = new AdjacentList();
            mAdjacentList[i]->mWhich = i;
        }
    }
}


/**
 * Useful for when we want to redo clustering with a set of already
 * initialized PDBs.
 */
void
Clustering::reinitialize(vector<char *> * nNames,
                         vector<Stru *> * nPDBs,
                         double threshold)
{
    THRESHOLD = threshold;

    // Read decoys and get threshold - = - = - = - = - = - = -

    CLU_RADIUS = THRESHOLD / 2.0 - 0.00001;
    mNames = nNames;
    mPDBs = nPDBs;
    int _mNumPDB = mNumPDB;
    mNumPDB = mPDBs->size();

    // Initialize auxiliary decoys grouping - = - = - = - = -

    // remove old groupings
    delete [] mD2C;
    delete [] mCen;
    for (int i=0; i < _mNumPDB; i++)
        if (mAuxCluster[i])
            delete mAuxCluster[i];
    delete [] mAuxCluster;
    delete mCluCen;

    // initialize new ones
    mCluCen = new vector<int>(0);
    mAuxCluster = new vector<int>*[mNumPDB];
    for (int i=0; i < mNumPDB; i++)
        mAuxCluster[i] = NULL;
    mD2C  = new double[mNumPDB];
    mCen  = new int[mNumPDB];

    // Initialize clustering - = - = - = - = - = - = - = - = -

    // delete old...
    for (int i=0; i < _mNumPDB; i++)
        delete mAdjacentList[i];
    delete mAdjacentList;

    // ...create new
    mAdjacentList = new AdjacentList*[mNumPDB];
    for (int i=0; i < mNumPDB; i++)
    {
        if(AdjacentList::mListMode == MATRIX)
        {
            mAdjacentList[i] = new AdjacentList(i, mNumPDB);
        }
        else
        {
            mAdjacentList[i] = new AdjacentList();
            mAdjacentList[i]->mWhich = i;
        }
    }

    // Delete cluster data - = - = - = - = - = - = - = - = -

    delete mFinalClusters;
    delete [] mRemainingList;
    delete [] mRemainingListIndex;
    delete [] mNumNeighbor;
}


/**
 * Get these: mNames, mPDBs, mLen
 */
void
Clustering::getThresholdAndDecoys()
{
    readDecoyNames(); // results in mNames

    // decide min max target cluster sizes - = - = - = -

#ifdef _USE_FIX_CLUSTER_SIZES_
    // good for when comparing with cluster_info_silent (ROSETTA)
    unsigned int minClusterSize = 15;
    unsigned int maxClusterSize = 75;
    unsigned int targetClusterSize = 45;
#else
    unsigned int minClusterSize = mNames->size()/100;
    unsigned int maxClusterSize = mNames->size()/15; // INFLUENTIAL
    unsigned int targetClusterSize = (maxClusterSize + minClusterSize)/2; // INFLUENTIAL
#endif
    if (maxClusterSize > (mNames->size()-1))
        maxClusterSize = mNames->size()-1;
    if (targetClusterSize > (mNames->size()-1))
        targetClusterSize = mNames->size()-1;

    // decide the min max thresholds - = - = - = - = -

    double minDist, maxDist, mostFreqDist, xPercentileDist;

    estimateDist(mNames,
                 NUM_TRIALS_FOR_THRESHOLD,
                 mNames->size() > 2*RANDOM_DECOY_SIZE_FOR_THRESHOLD?
                     RANDOM_DECOY_SIZE_FOR_THRESHOLD: (mNames->size()/2),
                 xPercentile,
                 &minDist,
                 &maxDist,
                 &mostFreqDist,
                 &xPercentileDist);

    if (EST_THRESHOLD == MOST_FREQ_BASED)
    {
        THRESHOLD = minDist + xFactor * (mostFreqDist-minDist) ;
        cout << "Finding threshold using most frequent distance" << endl;
        cout << "Threshold = " << THRESHOLD << "  ( " << minDist
             << " + " << xFactor << "(" << mostFreqDist
             << " - " << minDist << ") )" << endl;
    }

    if (EST_THRESHOLD == PERCENT_EDGES)
    {
        THRESHOLD = xPercentileDist;
        cout << "Finding threshold to keep " << xPercentile << "% edges" << endl;
        cout << "Threshold = " << THRESHOLD << endl;
    }

    if (EST_THRESHOLD == MIN_AVG_DIST_BASED)
    {
        // Uses a set of randomly chosen decoys to estimate the threshold
        // as in cluster_info_silent (i.e. rosetta)
        //
        // For each sampled decoy, we get its average distance to other decoys
        // Then, take the least of such average distance, m
        // Threshold is then computed as "a*m + b", where a and b are some
        // constant we decide by trial and error.
        //
        cout << "Finding threshold using average distance. "
             << "This may take a while..." << endl;
        clock_t start = clock();

        int numDecoys = mNames->size() > 101? 101: mNames->size();
        vector<char *>* names = getRandomDecoyNames(mNames, numDecoys, 0);
        vector<Stru *>* decoys = readDecoys(names);

        double ** nbors = getNborList(decoys, mNames, mNames->size()-1);

        destroyRandomDecoys(names, decoys);

        // find the minimum of average distances
        double avg_dist, min_avg_dist = _OVER_RMSD_;
        double sum_of_dist, * nbor;
        for (int i=0; i < numDecoys; i++)
        {
            nbor = nbors[i];
            sum_of_dist = 0;
            for (unsigned int j=0; j < mNames->size(); j++)
                sum_of_dist += nbor[j];
            avg_dist = sum_of_dist / mNames->size();
            if (avg_dist < min_avg_dist)
                min_avg_dist = avg_dist;
        }

        THRESHOLD = xFactor * min_avg_dist + minDist;

        double elapsed = (clock() - start)/(double)CLOCKS_PER_SEC;
        cout << "Minimum average distance = " << min_avg_dist
             << ". Found in " << elapsed << " s" << endl;
        cout << "Threshold = " << THRESHOLD << "  (" << minDist
             << " + " << xFactor << " * " << min_avg_dist << ")" << endl;
    }

    if (EST_THRESHOLD == SAMPLED_ROSETTA)
    {
        // Uses a set of randomly chosen decoys to estimate the threshold
        // as in cluster_info_silent (i.e. rosetta)
        //
        // The belief is that the largest cluster should contain a set X
        // number of decoys, within a threshold range (t1, t2).
        // This is what ROSETTA does:
        // 1. Build sorted lists of decoys, that is
        //      decoy1: neighbor1, neighbor2, ..., neighborX
        //      decoy2: neighbor1, neighbor2, ..., neighborX
        //       ...
        // 2. If there is a list such that X neighbors can be contained
        //    within a threshold t, t1 <= t <= t2, then we are done.
        //
        // In the following, only the sampled decoys get their lists built.
        //
        // Required parameters:
        //   mRandomDecoys
        //   {min,max,target}ClusterSize
        //   {min,max}Threshold
        //
        cout << "Finding threshold using random decoys (max cluster size = "
             << maxClusterSize << ")" << endl;
        clock_t start = clock();

        int numDecoys = mNames->size() > 101? 101: mNames->size();
        vector<char *>* names = getRandomDecoyNames(mNames, numDecoys, 0);
        vector<Stru *>* decoys = readDecoys(names);

        double ** nbors = getNborList(decoys, mNames, maxClusterSize);

        destroyRandomDecoys(names, decoys);

        double minThreshold = minDist;
        double maxThreshold = (minDist + maxDist)/2;

        THRESHOLD = getThreshold(nbors,
                                 numDecoys,      // # of lists to build
                                 maxClusterSize, // length of list
                                 minClusterSize,
                                 targetClusterSize,
                                 minThreshold,
                                 maxThreshold);

        double elapsed = (clock() - start)/(double)CLOCKS_PER_SEC;
        cout << "Threshold = " << THRESHOLD
             << ". Found in " << elapsed << " s" << endl;
    }

    // read decoys - = - = - = - = - = - = - = - = -

    vector<char *>* randNames;
    vector<Stru *>* randDecoys;
    if (FILTER_MODE)
    {
        int numDecoys = mNames->size() > 2*RANDOM_DECOY_SIZE_FOR_FILTERING?
                        RANDOM_DECOY_SIZE_FOR_FILTERING: (mNames->size()/2);
        randNames = getRandomDecoyNames(mNames, numDecoys, 0);
        randDecoys = readDecoys(randNames);
    }

    clock_t start = clock();
    readDecoys(randNames, randDecoys); // results in mPDBs
    double elapsed = (clock() - start)/(double)CLOCKS_PER_SEC;
    cout << "Decoys read in " << elapsed << " s" << endl;

    if (FILTER_MODE)
        destroyRandomDecoys(randNames, randDecoys);

    // find threshold using ROSETTA mode - = - = - = - = - = - = - = - = -
    // (this has to be done with the full decoys) - = - = - = - = - = - =

    if (EST_THRESHOLD == ROSETTA)
    {
        // This is exactly ROSETTA's way of getting threshold.
        // We add this in for comparison with the output of ROSETTA.
        //
        // Note that targetClusterSize is guaranteed since we pass the
        // minimum and maximum distance between decoys --obtained from
        // getNborList-- to getThreshold().
        //
        // Required parameters:
        //   mPDBs
        //   {min,max,target}ClusterSize
        //
        //   {min,max}Threshold are automatically decided accurately
        //   This is so that targetClusterSize will always pass
        //
        cout << "Finding threshold with ROSETTA's method (max cluster size = "
             << maxClusterSize << ")" << endl;
        clock_t start = clock();

        double minThreshold, maxThreshold;

        double ** nbors = getNborList(mPDBs,
                                     maxClusterSize,
                                     &minThreshold,
                                     &maxThreshold);

        THRESHOLD = getThreshold(nbors,
                                 mPDBs->size()-1,
                                 maxClusterSize,
                                 minClusterSize,
                                 targetClusterSize,
                                 minThreshold,
                                 maxThreshold);

        double elapsed = (clock() - start)/(double)CLOCKS_PER_SEC;
        cout << "Threshold = " << THRESHOLD
             << ". Found in " << elapsed << " s" << endl;
    }

}


/**
 * Repeatedly estimate distances for numTrials trials, and get the average
 * values of
 *    1. min of pairwise distance,
 *    2. max of pairwise distance,
 *    3. most frequent distance,
 *    4. distance t at which xPercent of pairwise distances are below t,
 * from all trials. The number of random decoys to use in each trial
 * is stated in randDecoySize.
 */
void
Clustering::estimateDist(vector<char *>* allNames,
                         int numTrials,
                         int randDecoySize,
                         double xPercent,
                         double *minDist,
                         double *maxDist,
                         double *mostFreqDist,
                         double *xPercentileDist)
{
    vector<char *>* randNames;
    vector<Stru *>* randDecoys;

    cout << "Estimating threshold range...";

    randNames = getRandomDecoyNames(allNames, randDecoySize, 0);
    randDecoys = readDecoys(randNames);

    estimateDist(randDecoys,
                 xPercent,
                 minDist,
                 maxDist,
                 mostFreqDist,
                 xPercentileDist);

    destroyRandomDecoys(randNames, randDecoys);

#ifdef _ANALYZE_RANDOM_SAMPLES_
    cout << endl
         << "sampling 1: dist [" << *minDist << "," << *maxDist
         << "], most freq = " << *mostFreqDist << ", "
         << xPercent << " percentile = " << *xPercentileDist
         << endl;
#elif defined(_SHOW_PERCENTAGE_COMPLETE_)
    printf("\rEstimating threshold range... %4.1f%%\r", 100. /numTrials);
    fflush(stdout);
#endif

    double * maxDists = new double[numTrials-1];
    double * minDists = new double[numTrials-1];
    double * mostFreqDists = new double[numTrials-1];
    double * xPercentileDists = new double[numTrials-1];
    for (int i=0; i < numTrials-1; i++)
    {
        randNames = getRandomDecoyNames(allNames, randDecoySize, 1 + i*100);
        randDecoys = readDecoys(randNames);
        estimateDist(randDecoys,
                     xPercent,
                     &minDists[i],
                     &maxDists[i],
                     &mostFreqDists[i],
                     &xPercentileDists[i]);

        destroyRandomDecoys(randNames, randDecoys);

#ifdef _ANALYZE_RANDOM_SAMPLES_
        cout << "sampling " << (i+2) << ": dist [" << minDists[i] << ","
             << maxDists[i] << "], most freq = "
             << mostFreqDists[i] << ", " << xPercent
             << " percentile = " << xPercentileDists[i] << endl;
#elif defined(_SHOW_PERCENTAGE_COMPLETE_)
        printf("Estimating threshold range... %4.1f%%\r", 100.*(i+2)/numTrials);
        fflush(stdout);
#endif
    }

#if !defined(_ANALYZE_RANDOM_SAMPLES_) && !defined(_SHOW_PERCENTAGE_COMPLETE_)
    cout << " done" << endl;
#elif defined(_SHOW_PERCENTAGE_COMPLETE_)
    cout << "Estimating threshold range... complete" << endl;
#endif

    for (int i=0; i < numTrials-1; i++)
    {
        if (minDists[i] < *minDist)
            *minDist = minDists[i];
        if (maxDists[i] > *maxDist)
            *maxDist = maxDists[i];
        *mostFreqDist += mostFreqDists[i];
        *xPercentileDist += xPercentileDists[i];
    }
    *mostFreqDist /= numTrials;
    *xPercentileDist /= numTrials;

    cout << "Estimated threshold range = ["
         << *minDist << ", " << *maxDist << "]" << endl
         << "Most frequent distance = " << *mostFreqDist
         << endl
         << xPercent << " percent distances below " << *xPercentileDist
         << endl;
}


/**
 * Read from the input file the names of decoys
 */
void
Clustering::readDecoyNames()
{
    switch (filetype(mInputFileName))
    {
    PreloadedPDB * pdbs;
    case SILENT_FILE:
        pdbs = new PreloadedPDB();
        pdbs->loadSilentFile(mInputFileName); // Preload PDBs from silent file
        SimPDB::preloadedPDB = pdbs; // Attach the preloaded PDBs to SimPDB
        SimPDB::preloadPDB = true;
        mNames = pdbs->mNames;
        cout << "Read " << mNames->size() << " decoy names" << endl;
        return;
    case PDB_LIST:
        if (SimPDB::preloadPDB)
        {
            unsigned int numDecoys = num_lines_in_file(mInputFileName);
            if (numDecoys > PreloadedPDB::ADVISED_THRESHOLD)
            {
                cout << "More than " << PreloadedPDB::ADVISED_THRESHOLD
                     << " decoys: not preloading PDBs" << endl;
                SimPDB::preloadPDB = false;
            }
        }
        if (!SimPDB::preloadPDB) // either switched off by user, or by above
            break;
        pdbs = new PreloadedPDB();
        pdbs->loadPDBFromList(mInputFileName); // Preload PDBs from pdb list
        SimPDB::preloadedPDB = pdbs; // Attach the preloaded PDBs to SimPDB
        mNames = pdbs->mNames;
        cout << "Read " << mNames->size() << " decoy names" << endl;
        return;
    default:
        break;
    }

    ifstream input(mInputFileName);
    if (!input)
    {
        cerr << "Cannot find file \"" << mInputFileName << "\"" << endl;
        exit(0);
    }

    char buf[400];
    char* token;

    mNames = new vector<char *>(0);
    while (!input.eof())
    {
        input.getline(buf, 400);
        token = strtok(buf, " ");
        if(token == NULL) continue;
        char* name = new char[strlen(token)+1];
        strcpy(name, token);
        mNames->push_back(name);
    }
    input.close();
    cout << "Read " << mNames->size() << " decoy names" << endl;
}


void
Clustering::allocateSpaceForRMSD(int len)
{
    if (spaceAllocatedForRMSD)
        return;
    result_coords = new double[3*len];
    spaceAllocatedForRMSD = true;
}


/**
 * Read from the input files the decoys and return them
 */
vector<Stru *>*
Clustering::readDecoys(vector<char *>* decoynames)
{
    SimPDB* firstPDB = new SimPDB((*decoynames)[0]);
    mLen = firstPDB->mNumResidue;

    allocateSpaceForRMSD(mLen);

    vector<Stru *> * decoys = new vector<Stru* >(decoynames->size());
    (*decoys)[0] = new Stru(firstPDB, mLen);

    for (unsigned int i=1; i < decoynames->size(); i++)
    {
        (*decoys)[i] = new Stru(new SimPDB((*decoynames)[i], mLen), mLen);
    }
    return decoys;
}


/**
 * Read from the input files all the decoys, filtering when necessary
 */
void
Clustering::readDecoys(vector<char *>* randomNames,
                       vector<Stru *>* randomDecoys)
{
    vector<char *>* newNames = new vector<char *>(0);
    vector<Stru *>* newDecoys = new vector<Stru *>(0);
#ifdef _SPICKER_SAMPLING_
    // Spicker samples decoys at a fixed interval delta
    // such that exactly 13000 decoys are sampled
    cout << "Sampling 13000 decoys" << endl;
    double delta = mNames->size()/(double)13000;
    int sampled_decoy_id = 1;
    if (delta < 1)
        delta = 1;
#endif

    cout << "Reading decoys with signature mode " << (Pref::use_sig?"on":"off")
         << " and filter mode " << (FILTER_MODE? "on": "off") << endl;

    int size = mNames->size();
    for (int i=0; i < size; i++)
    {
#ifdef _SPICKER_SAMPLING_
        if ((i+1) < sampled_decoy_id * delta)
            continue;
        sampled_decoy_id++;
#endif
#ifdef _SHOW_PERCENTAGE_COMPLETE_
        printf("Read %4.1f%%\r", 100.*i/size);
        fflush(stdout);
#endif
        // read in the decoy's PDB
        char* dName = (*mNames)[i];
        Stru* s;
        if (mLen == 0)
        {
            SimPDB* aPDB = new SimPDB(dName);
            mLen = aPDB->mNumResidue;
            s = new Stru(aPDB, mLen);
            allocateSpaceForRMSD(mLen);
        }
        else
        {
            s = new Stru(new SimPDB(dName, mLen), mLen);
        }
        bool isOutlier = false;
        if (FILTER_MODE) // then we shall decide whether to include s
        {
            int randomDecoysSize = randomDecoys->size();
            isOutlier = true;
            for(int j=0; j < randomDecoysSize; j++)
            {
                if (strcmp(dName,(*randomNames)[j]) == 0)
                    continue;
                if (Pref::use_sig && estD(s,(*randomDecoys)[j]) > 2*THRESHOLD)
                    continue;
                if (trueD(s,(*randomDecoys)[j]) > 2*THRESHOLD)
                    continue;
                // s is near to a random decoy
                isOutlier = false;
                break;
            }
        }
        if (!isOutlier)
        {
            newNames->push_back(dName);
            newDecoys->push_back(s);
        }
        else
        {
            delete s;
            delete [] dName;
        }
    }
    cout << "Read " << newNames->size() << " decoys.";
    if (FILTER_MODE)
        cout << " Filtered " << (mNames->size() - newNames->size())
             << " outlier decoys.";
    cout << endl;
    delete mNames;
    mNames = newNames;
    mPDBs  = newDecoys;
    mNumPDB = mPDBs->size();
    if (!mNumPDB)
    {
        cout << "No decoy read. Exiting..." << endl;
        exit(0);
    }
}


/**
 * Re-filter the content of mPDBs, mainly due to the change of threshold.
void
Clustering::refilterDecoys(vector<char *>* randomNames,
                           vector<Stru *>* randomDecoys)
{
    vector<char *>* newNames = new vector<char *>(0);
    vector<Stru *>* newDecoys = new vector<Stru *>(0);
    for (int i=0; i < mNumPDB; i++)
    {
        char* dName = (*mNames)[i];
        Stru* s = (*mPDBs)[i];
        bool isOutlier = true;
        if (FILTER_MODE) // then we shall decide whether to include s
        {
            int randomDecoysSize = randomNames->size();
            for(int j=0; j < randomDecoysSize; j++)
            {
                if (strcmp(dName,(*randomNames)[j]) == 0)
                    continue;
                if (Pref::use_sig && estD(s,(*randomDecoys)[j]) > 2*THRESHOLD)
                    continue;
                if (trueD(s,(*randomDecoys)[j]) > 2*THRESHOLD)
                    continue;
                // s is near to a random decoy
                isOutlier = false;
                break;
            }
        }
        else
        {
            isOutlier = false;
        }
        if (!isOutlier)
        {
            newNames->push_back(dName);
            newDecoys->push_back(s);
        }
        else
        {
            delete s;
            delete [] dName;
        }
    }
    delete mNames;
    mNames = newNames;
    mPDBs  = newDecoys;
    cout << "Re-filtered to " << mPDBs->size() << " decoys" << endl;
}
 */


/**
 * Main call. Find all neighbors within threshold.
 */
void
Clustering::cluster()
{
    clock_t start;
    double elapsed;

    cout << "Auxiliary clustering...";
#ifdef _SHOW_PERCENTAGE_COMPLETE_
    printf("\r");
#endif
    start = clock();
    initRef(NULL);
    auxClustering(); // Cluster to speed-up computation of neighbors
    elapsed = (clock() - start)/(double)CLOCKS_PER_SEC;
#ifdef _SHOW_PERCENTAGE_COMPLETE_
    cout << "\nAuxiliaryClustering...";
#endif
    cout << " completed in " << elapsed << " s" << endl;

    cout << "Finding decoys neighbors...";
    //start = clock();
#ifndef __WIN32__
    _get_elapsed(1);
#endif
    buildAdjacentLists(); // Find all the neighbors for each decoy
    //elapsed = (clock() - start)/(double)CLOCKS_PER_SEC;
#ifndef __WIN32__
    elapsed = _get_elapsed(0);
    cout << "\nFinding all neighbors completed in " << elapsed << " s" << endl;
#else
    cout << " completed" << endl;
#endif

    //listAdjacentLists();

    cout << "Finding and removing largest clusters...";
    start = clock();
    findLargestClusters(); // Find largest clusters and remove. Recurse.
    elapsed = (clock() - start)/(double)CLOCKS_PER_SEC;
    cout << " completed in " << elapsed << " s" << endl;

#ifdef _ADD_LITE_MODE_
    return;
#endif

    /**
     * If the first cluster is not much bigger than the second cluster,
     * two independent clusters may exist. Hence if this difference is
     * below a certain margin we may want to do something about it.
     *
     * We want an almost constant margin so we compare the sizes of the two
     * clusters prior to the removal of decoys from the second. If we use
     * the size of the clusters after the removal, we will have to consider
     * the fact that: at higher thresholds the overlap between clusters will
     * will increase, and the margin will increase likewise.
     *
     * That the second cluster is not originally large due to its similarity
     * to the first cluster (i.e. they are two almost equal center from the
     * same cluster) is more or less guaranteed from the fact that, if this
     * is the case, then most of its decoys would've been removed during the
     * clustering and it would have lost out to some other cluster that
     * has less overlap with the best cluster. Hence the presence of at least
     * one such "other cluster* should be ensured.
     */
    AdjacentList* bestClus = (*mFinalClusters)[0];
    bestClusSize = bestClus->mNumNeigh;
    int nextSize = 0;
    if (mFinalClusters->size() > 1)
    {
        AdjacentList* nextClus = (*mFinalClusters)[1];
        nextSize = mNumNeighbor[nextClus->mWhich];
        bestClusMargin = (bestClusSize - nextSize) /(double) bestClusSize;
        cout << "Largest 2 clusters: "
             << (*mNames)[bestClus->mWhich] << "(" << bestClusSize << "), "
             << (*mNames)[nextClus->mWhich] << "(" << nextSize
             << "). Margin = " << (bestClusMargin*100) << "%" << endl;
    }
}


void
Clustering::showClusters(int numClus)
{
#ifdef _ADD_LITE_MODE_
    if (AdjacentList::mListMode == LITE)
    {
        cout << "Best decoy: " << (*mNames)[mFinalDecoy] << "\t"
             << mAdjacentList[mFinalDecoy]->mNumNeigh << endl;
        return;
    }
#endif

    int size = mFinalClusters->size();
    if (size > numClus)
        size = numClus;
    for (int i=0; i<size; i++)
    {
        AdjacentList* clu = (*mFinalClusters)[i];
        //cout << clu->mWhich << "'" << (*mNames)[clu->mWhich]
        cout << (*mNames)[clu->mWhich]
             << " " << clu->mNumNeigh << ": ";
        int size2 = clu->mNumNeigh;
        vector<int>* neigh = clu->neigh;
        for (int j=0; j < size2; j++) // for each neighboring decoy...
        {
            // output its decoy number and name...
            int decoy = (*neigh)[j];
            //cout << decoy << "'" << (*mNames)[decoy] << "'";
            cout << (*mNames)[decoy] << "\t";
            // ...and its distance
            //if (AdjacentList::mListMode == LIST)
            //    cout << (*(clu->dist))[j] << ", ";
            //else
            //    cout << clu->getD(decoy) << ", ";
        }
        cout << endl;
    }
}


/**
 * Get the names and Strus of the (first num_elements) indices in list into
 * Names and PDBs respectively.
 */
void
Clustering::getPDBs(vector<char *>* Names, vector<Stru *>* PDBs,
                    vector<int>* list, int num_elements)
{
    int x;
    for (int i=0; i < num_elements; i++)
    {
        x = (*list)[i];
        PDBs->push_back((*mPDBs)[x]);
        Names->push_back((*mNames)[x]);
    }
}


void
Clustering::listAdjacentLists()
{
    cout << "Listing all adjacentLists" << endl;
    for (int i=0; i < mNumPDB; i++)
    {
        AdjacentList* clu = mAdjacentList[i];
        cout << i << "," << clu->mWhich << "," << (*mNames)[clu->mWhich]
             << " " << clu->mNumNeigh << ": ";
        int size2 = clu->mNumNeigh;
        for (int j=0; j < size2; j++) // for each neighboring decoy...
        {
            // output its decoy number and name...
            cout << (*(clu->neigh))[j] << "'" << (*mNames)[(*(clu->neigh))[j]] << "'";
            // ...and its distance
            if (AdjacentList::mListMode == LIST)
                cout << (*(clu->dist))[j] << ", ";
            else
                cout << clu->getD((*(clu->neigh))[j]) << ", ";
        }
        cout << endl;
    }
    cout << endl;
}


//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -
// Codes for building the AdjacentLists
//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

/**
 * A simplistic clustering to help us find the adjacent lists --- this is one
 * of Calibur's strategy.
 * In the output, each decoy either belongs to a cluster, or is a cluster
 * center.
 *
 * For each decoy from 0... mNumPDB,
 *   If it is not within CLUS_RADIUS of any currently registered cluster center
 *     Add it as a cluster center
 *   Else add it to the cluster
 */
void
Clustering::auxClustering()
{
    double lower, upper, upper_scud;
    for (int i=0; i < mNumPDB; i++) // for each decoy
    {
#ifdef _SHOW_PERCENTAGE_COMPLETE_
        printf("Auxiliary clustering... %4.1f%%\r", 100.*i/mNumPDB);
        fflush(stdout);
#endif
        bool clustered = false;
        int numc = mCluCen->size();
        for (int c = 0; c < numc; c++) // for each center
        {
            int cen = (*mCluCen)[c];

            refBound(i, cen, lower, upper);
            if (lower > CLU_RADIUS)
                continue;

            if (Pref::use_sig && estD(i, cen) > CLU_RADIUS)
                continue;

            if (upper <= CLU_RADIUS)
            {
                mD2C[i] = upper;
                mCen[i] = cen;
                clustered = true;
                mAuxCluster[cen]->push_back(i);
                break;
            }

            if (Pref::use_scud)
            {
                upper_scud = eucD(i, cen);
                if (upper_scud <= CLU_RADIUS)
                {
                    mD2C[i] = upper_scud;
                    mCen[i] = cen;
                    clustered = true;
                    mAuxCluster[cen]->push_back(i);
                    break;
                }
            }

            double d = trueD(i, cen);
            mAdjacentList[i]->add(cen, d, false);
            mAdjacentList[cen]->add(i, d, false);
            if (d <= CLU_RADIUS)
            {
                mD2C[i] = d;
                mCen[i] = cen;
                clustered = true;
                mAuxCluster[cen]->push_back(i);
                break;
            }
        }
        if (!clustered)
        {
            /**
             * Declare i as a cluster center.
             * Important: itself must be added to its cluster.
             */
            mAuxCluster[i] = new vector<int>(0);
            mAuxCluster[i]->push_back(i);
            mCluCen->push_back(i);
            mD2C[i] = 0;
            mCen[i] = i;
        }
    }
}


/** too slow
bool
Clustering::find(int which, vector<int> *elements)
{
    int len=elements->size();
    for (int i=0; i<len; i++)
    {
        if ((*elements)[i] == which) return true;
    }
    return false;
}
*/

/**
 * Subroutine for cluster().
 * Finds the neighbors for each decoy i within threshold and add these
 * neighbors in mAdjacentList[i]
 */
void
Clustering::buildAdjacentLists()
{
    double lower, upper;
    //, upper_scud;
    int numc = mCluCen->size();
    //cout << "Number of decoys=" << mNumPDB
    //     << ", number of clusters=" << numc << endl;

    //=================================================================
    // We assume that every element belongs to exactly one cluster,
    // so we loop through all the clusters and determine if the cluster
    // element should be added to the neighbor set of some decoy.
    //=================================================================
#ifdef _SHOW_PERCENTAGE_COMPLETE_
    printf("\r");
#endif
    for (int c=0; c < numc; c++) // for each cluster center
    {
        double d, _d;
        int cen = (*mCluCen)[c];
#ifdef _SHOW_PERCENTAGE_COMPLETE_
        printf("Finding decoys' neighbors... completed %4.1f%%\r", 100.*c/numc);
        fflush(stdout);
#endif
        for (int i=0; i < mNumPDB; i++) // for each decoy
        {
            vector<int>* elements = mAuxCluster[cen];
            int size = elements->size();

            //=========================================================
            // Case 1
            // Under this condition, every element in the currrent
            // cluster is a neighbor of each other.
            // Do not compute distance.
            // Set dist=-1 to indicate that it's added at this stage.
            //=========================================================
            if (2*CLU_RADIUS <= THRESHOLD)
            {
                if (mCen[i] == cen) // if i is an element of the cluster
                {
#ifdef _ADD_LITE_MODE_
                    if (AdjacentList::mListMode == LITE)
                    {
                        mAdjacentList[i]->add(size);
                        continue;
                    }
#endif
                    for (int n=0; n < size; n++)
                    {
                        int elem = (*elements)[n];
                        //if (elem != i) // do not add self
                            mAdjacentList[i]->add(elem, -1., true);
                    }
                    continue;
                }
            }

            //=========================================================
            // Case 2 and Case 3
            // (Use signature estimation for pruning)
            // If this condition is fulfilled none of the element in
            // the cluster is a neighbor of the current decoy
            //=========================================================

            refBound(i, cen, lower, upper);
            if (lower - CLU_RADIUS > THRESHOLD)
            {
                continue;
            }
            if (upper + CLU_RADIUS <= THRESHOLD)
            {
#ifdef _ADD_LITE_MODE_
                if (AdjacentList::mListMode == LITE)
                {
                    mAdjacentList[i]->add(size);
                    continue;
                }
#endif
                for (int n=0; n < size; n++)
                {
                    mAdjacentList[i]->add((*elements)[n], -6, true);
                }
                continue;
            }

            if (Pref::use_scud && eucD(i,cen) + CLU_RADIUS <= THRESHOLD)
            {
#ifdef _ADD_LITE_MODE_
                if (AdjacentList::mListMode == LITE)
                {
                    mAdjacentList[i]->add(size);
                    continue;
                }
#endif
                for (int n=0; n < size; n++)
                {
                    mAdjacentList[i]->add((*elements)[n], -6, true);
                }
                continue;
            }

            if (Pref::use_sig && estD(i, cen) - CLU_RADIUS > THRESHOLD)
                continue;

            //=========================================================
            // Find exact distance to center
            //=========================================================
            d = trueD(i, cen);
            mAdjacentList[i]->add(cen, d, false);
            mAdjacentList[cen]->add(i, d, false);

            //=========================================================
            // If this condition is fulfilled none of the element in
            // the cluster is a neighbor of the current decoy
            //=========================================================
            if (d - CLU_RADIUS > THRESHOLD)
                continue;

            //=========================================================
            // If this condition is fulfilled every element in
            // the cluster is a neighbor of the current decoy
            // When r <= d/2, this is also implied by
            //    d <= THRESHOLD/2. or
            //    d <= r
            // Both of which are more strict than d + r <= THRESHOLD
            //=========================================================
            if (d + CLU_RADIUS <= THRESHOLD)
            {
#ifdef _ADD_LITE_MODE_
                if (AdjacentList::mListMode == LITE)
                {
                    mAdjacentList[i]->add(size);
                    continue;
                }
#endif
                for (int n=0; n < size; n++)
                {
                    mAdjacentList[i]->add((*elements)[n], -2, true);
                }
                continue;
            }

            //========================================================
            // Consider each cluster element individually
            //========================================================
            for (int j=0; j < size; j++)
            {
                int e = (*elements)[j];

                if (mD2C[e]+d <= THRESHOLD)
                {
#ifdef _ADD_LITE_MODE_
                    if (AdjacentList::mListMode == LITE)
                    {
                        mAdjacentList[i]->add(1);
                        continue;
                    }
#endif
                    mAdjacentList[i]->add(e, -3, true);
                }
                /*
                else if (fabs(d-mD2C[e]) > THRESHOLD)
                    continue;
                */
                else
                {
                    if (Pref::use_sig && estD(i, e) > THRESHOLD)
                        continue;

                    refBound(i, e, lower, upper);
                    if (upper <= THRESHOLD)
                    {
#ifdef _ADD_LITE_MODE_
                        if (AdjacentList::mListMode == LITE)
                        {
                             mAdjacentList[i]->add(1);
                             continue;
                        }
#endif
                        mAdjacentList[i]->add(e, -4, true);
                        continue;
                    }

                    if (Pref::use_scud && eucD(i, e) <= THRESHOLD)
                    {
#ifdef _ADD_LITE_MODE_
                        if (AdjacentList::mListMode == LITE)
                        {
                             mAdjacentList[i]->add(1);
                             continue;
                        }
#endif
                        mAdjacentList[i]->add(e, -7, true);
                        continue;
                    }

                    if (lower > THRESHOLD)
                        continue;

                    _d = trueD(i, e);
                    mAdjacentList[e]->add(i, _d, false);
                    mAdjacentList[i]->add(e, _d, false);
                    if (_d <= THRESHOLD)
                    {
#ifdef _ADD_LITE_MODE_
                        if (AdjacentList::mListMode == LITE)
                        {
                             mAdjacentList[i]->add(1);
                             continue;
                        }
#endif
                        mAdjacentList[i]->add(e, _d, true);
                    }
                }
            }
        }
    }
}

//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -
// Codes for finding RMSD
//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

double
Clustering::estD(int i, int j)
{
    double* sig1 = (*mPDBs)[i]->mSIG;
    double* sig2 = (*mPDBs)[j]->mSIG;
    double rev = 0;
    for (int c=0; c < mLen; c++)
    {
        double v = sig1[c]-sig2[c];
        rev += v*v;
    }
    return sqrt(rev/mLen);
}

double
Clustering::estD(Stru* a, Stru* b)
{
    double* sig1=a->mSIG;
    double* sig2=b->mSIG;
    double rev=0;
    for (int c=0; c<mLen; c++)
    {
        double v = sig1[c]-sig2[c];
        rev += v*v;
    }
    return sqrt(rev/mLen);
}

double
Clustering::trueD(int i, int j)
{
    cout.flush();
    if (i == j)
        return 0;
    double d = mAdjacentList[i]->getD(j);
    if (d < _OVER_RMSD_ && d >= 0)
        return d;
    double* coor1 = (*mPDBs)[i]->mCAlpha;
    double* coor2 = (*mPDBs)[j]->mCAlpha;
    double rmsd;
    rmsd = fast_rmsd(coor1, coor2, mLen);
    if (rmsd != rmsd) // crazy RMSD
        rmsd = RMSD(coor1, coor2, mLen);
    return (double) rmsd;
}

double
Clustering::eucD(int i, int j)
{
    cout.flush();
    if (i == j)
        return 0;
    double d = mAdjacentList[i]->getD(j);
    if (d < _OVER_RMSD_ && d >= 0)
        return d;
    double* coor1 = (*mPDBs)[i]->mCAlpha;
    double* coor2 = (*mPDBs)[j]->mCAlpha;
    double rev=0;
    int l3=mLen*3;
    for (int k=0; k<l3; k++)
    {
        double d=coor1[k]-coor2[k];
        rev+=d*d;
    }
    return sqrt(rev/mLen);
}


/**
 * Realign decoys so that their point-wise distance to decoy "ref" is minimized
 */
void
Clustering::realignDecoys(int ref)
{
    cout << "Realigning decoys...";
    clock_t start = clock();
    for (int i=1; i < mNumPDB; i++)
        superimposeAndReplace((*mPDBs)[ref]->mCAlpha, (*mPDBs)[i]->mCAlpha);
    double elapsed = (clock() - start)/(double)CLOCKS_PER_SEC;
    cout << " completed in " << elapsed << " s" << endl;
}


void
Clustering::superimposeAndReplace(double* coor1, double* coor2)
{
    double R[3][3];
    RMSD(coor1, coor2, mLen, R);
    rotate(coor2, mLen, R, result_coords);
    for (int i=0; i < 3*mLen; i++)
        coor2[i] = result_coords[i];
}


double
Clustering::trueD(Stru* a, Stru* b)
{
    double* coor1=a->mCAlpha;
    double* coor2=b->mCAlpha;
    double rmsd = 0;
    rmsd = fast_rmsd(coor1, coor2, mLen);
    if (rmsd != rmsd) // crazy RMSD
        rmsd = RMSD(coor1, coor2, mLen);
    return rmsd;
}


//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -
// Codes for recursively finding the largest clusters
//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

/**
 * Find from all AdjacentLists the one that is longest
 */
int
Clustering::findDecoyWithMostNeighbors()
{
    int largest = mRemainingList[0];
    int currentSize = mAdjacentList[largest]->mNumNeigh;
    for (int i=1; i < mRemainingSize; i++)
    {
        int next_index = mRemainingList[i];
        int next_size = mAdjacentList[next_index]->mNumNeigh;
        if (next_size > currentSize)
        {
            currentSize = next_size;
            largest = next_index;
        }
    }
    return largest;
}


/**
 * Remove all the decoys in adj from all AdjacentList of all decoys.
 */
void
Clustering::removeDecoys(AdjacentList* adj)
{
    // remove from the remainingList
    vector<int>* neigh = adj->neigh;
    int size = adj->mNumNeigh;
    int this_list = adj->mWhich;
    for (int i=0; i<size; i++)
    {
        int to_remove = (*neigh)[i]; // for each decoy to_remove

        // remove "to_remove" from the RemainingList
        int index = mRemainingListIndex[to_remove];
        mRemainingList[index] = mRemainingList[mRemainingSize-1];
        mRemainingListIndex[mRemainingList[index]] = index;
        mRemainingSize--;

        // remove "to_remove" from the AdjacentList of all its neighbors
        AdjacentList* _neighbors_of_to_remove = mAdjacentList[to_remove];
        vector<int>* neighbors_of_to_remove = _neighbors_of_to_remove->neigh;
        int size2 = _neighbors_of_to_remove->mNumNeigh;
        for (int j=0; j<size2; j++)
        {
            // should not remove to_remove from current list. it is used as the final cluster :D
            if ((*neighbors_of_to_remove)[j] != this_list)
            {
                mAdjacentList[(*neighbors_of_to_remove)[j]]->remove(to_remove);
            }
        }
    }
}


void
Clustering::findLargestClusters()
{
    mFinalClusters = new vector<AdjacentList *>(0);

    // We use a list of the remaining decoys (that has not been removed)
    // to enumerate only the decoys that are still available.
    mRemainingList = new int[mNumPDB];
    mRemainingListIndex = new int[mNumPDB]; // map decoy index to mRemainingList index
    mRemainingSize = mNumPDB;

    mNumNeighbor = new int[mNumPDB];

    for (int i=0; i<mNumPDB; i++)
    {
        mRemainingList[i] = i;
        mRemainingListIndex[i] = i;
        // Remember the original size of the cluster
        mNumNeighbor[i] = mAdjacentList[i]->mNumNeigh;
    }

#ifdef _ADD_LITE_MODE_
    if (AdjacentList::mListMode == LITE)
    {
        mFinalDecoy = findDecoyWithMostNeighbors();
        return;
    }
#endif

    while (mRemainingSize > 1)
    {
        int largest = findDecoyWithMostNeighbors();
        AdjacentList* adj = mAdjacentList[largest];
        mFinalClusters->push_back(adj);
        removeDecoys(adj);
    }
}


/*
void
Clustering::findFinalRMSD()
{
    int size = mFinalClusters->size();
    for(int i=0; i < size; i++)
    {
        AdjacentList* adj=(*mFinalClusters)[i];
        int w = adj->mWhich;
        vector<int>* neigh = adj->neigh;
        int num_neigh = adj->mNumNeigh;
        for (int j=0; j < num_neigh; j++)
        {
            adj->add((*neigh)[j], trueD(w, (*neigh)[j]), false);
        }
    }
}
*/



//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -
// Random decoys finding codes
//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

/**
 * Randomly permute the numbers, and select the first size decoys
 */
vector<char *>*
Clustering::getRandomDecoyNames(vector<char *>* srcnames, int size, int seed)
{
    srand(time(NULL)/2+seed);
    int totalsize = srcnames->size();
    int* randomArray = new int[totalsize];
    for (int i = 0; i < totalsize; i++) // create an array of 0,...,totalsize-1
        randomArray[i] = i;

    // swap the first size elements with other elements from the whole array
    for (int i = 0; i < size; i++)
    {
        // find an index j (i<j<totalsize) to swap with the element i
        int j = i + rand() %(totalsize-i);
        int t = randomArray[j];
        randomArray[j] = randomArray[i];
        randomArray[i] = t;
    }
    // copy the first randomDecoysSize elements into a smaller array
    vector<char *> * names = new vector<char *>(0);
    for (int i = 0; i < size; i++)
    {
        names->push_back((*srcnames)[randomArray[i]]);
    }
    delete [] randomArray;
    return names;
}


/**
 * Destroy the random decoys
 */
void
Clustering::destroyRandomDecoys(vector<char *>* names, vector<Stru *>* decoys)
{
    delete names; // these names are shared

    // decoys should be read using readDecoys(), and should be safe to delete
    int size = decoys->size();
    for (int i = 0; i < size; i++)
    {
        delete (*decoys)[i];
    }
    delete decoys;
}


//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -
// Threshold finding codes
//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

/**
 * For each decoy c in decoys, build a list of listLength nearest
 * neighbors of c, taken from the list nborsCandidates.
 * FIXME there is no way to tell if a decoy in sampleDecoys is the same
 *       as a decoy in nborsCandidates, and hence a decoy can be inserted
 *       into its own list.
 */
double **
Clustering::getNborList(vector<Stru *> * decoys,
                        vector<char *> * nborsCandidates,
                        int listLength)
{
    // For each decoy i find its neighbors: nbors[0],nbors[1],...,nbors[N-1]
    double ** nbors, * nl, r;
    char * dName;
    Stru *a, *b;
    int N = decoys->size();   // decoys to build neighbor list of
    int M = nborsCandidates->size(); // #candidates to test to fill list
    if (listLength > M)
        listLength = M;

    nbors = (double **) calloc(N, sizeof(double *));
    for (int i=0; i < N; i++)
    {
        nbors[i] = (double *) calloc(listLength, sizeof(double));
        for (int j=0; j < listLength; j++)
            nbors[i][j] = 1000.0;
    }
    for (int i=0; i < M; i++) // for each candidate...
    {
        dName = (*nborsCandidates)[i];
        a = new Stru(new SimPDB(dName, mLen), mLen);
        for (int j=0; j < N; j++) // ...insert it into each candidate list.
        {
            nl = nbors[j];
            b = (*decoys)[j];

            // FIXME need to change trueD to store computed RMSDs in a cache.
            //       this cache should be used to pre-fill the AdjacentLists
            r = trueD(a, b);

            if (r < nl[listLength-1]) // nbors[I] should be altered
            {
                nl[listLength-1] = r;
                // double r up the list nl till it's larger
                for (int k=listLength-1; k > 0 && nl[k-1] > r; k--)
                {
                    // swap nl[k] with nk[k-1]
                    nl[k] = nl[k-1];
                    nl[k-1] = r;
                }
            }
        }
        delete a;
    }
    return nbors;
}

/**
 * Same as getNeighborList(vector<Stru *> *, vector<char *> *, int)
 * but with the second parameter already read as Stru.
double **
Clustering::getNborList(vector<Stru *> * decoys,
                        vector<Stru *> * nborsCandidates,
                        int listLength)
{
    // For each decoy i find its neighbors: nbors[0],nbors[1],...,nbors[N-1]
    double ** nbors, * nl, r;
    char * dName;
    Stru *a, *b;
    int N = decoys->size();   // decoys to build neighbor list of
    int M = nborsCandidates->size(); // #candidates to test to fill list
    if (listLength > M)
        listLength = M;

    nbors = (double **) calloc(N, sizeof(double*));
    for (int i=0; i < N; i++)
    {
        nbors[i] = (double *) calloc(listLength, sizeof(double));
        for (int j=0; j < listLength; j++)
            nbors[i][j] = 1000.0;
    }
    for (int i=0; i < M; i++) // for each candidate...
    {
        a = (*nborsCandidates)[i];
        for (int j=0; j < N; j++) // ...insert it into each candidate list.
        {
            nl = nbors[j];
            b = (*decoys)[j];

            // FIXME need to change trueD to store computed RMSDs in a cache.
            // this cache should be consulted when building the AdjacentLists
            r = trueD(a, b);

            if (r < nl[listLength-1]) // nbors[I] should be altered
            {
                nl[listLength-1] = r;
                // double r up the list nl till it's larger
                for (int k=listLength-1; k > 0 && nl[k-1] > r; k--)
                {
                    // swap nl[k] with nk[k-1]
                    nl[k] = nl[k-1];
                    nl[k-1] = r;
                }
            }
        }
    }
    return nbors;
}
 */


/**
 * Builds neighbor lists using the input decoys.
 * For use in ROSETTA mode.
 * Same as the earlier getNborList but nborsCandidates is taken to be
 * the same as decoys.
 * Also, minDist and maxDist are computed so that they can be used to
 * guarantee targetClusterSize is used in getThreshold().
 */
double **
Clustering::getNborList(vector<Stru *> * decoys,
                        int listLength,
                        double *minDist,
                        double *maxDist)
{
    // For each decoy i find its neighbors: nbors[0],nbors[1],...,nbors[N-1]
    double ** nbors, * nl, r;
#ifdef _DEBUG_NBORLIST_
    int ** nborsIndex, * ni;
#endif
    //char * dName;
    Stru *a, *b;
    int N = decoys->size();   // decoys to build neighbor list of
    if (listLength > N)
        listLength = N;
    *maxDist = 0;
    *minDist = _OVER_RMSD_;

    nbors = (double **) calloc(N, sizeof(double*));
#ifdef _DEBUG_NBORLIST_
    nborsIndex = (int **) calloc(N, sizeof(int *));
#endif
    for (int i=0; i < N; i++)
    {
        nbors[i] = (double *) calloc(listLength, sizeof(double));
#ifdef _DEBUG_NBORLIST_
        nborsIndex[i] = (int *) calloc(listLength, sizeof(double));
#endif
        for (int j=0; j < listLength; j++)
            nbors[i][j] = 1000.0;
    }
    for (int i=0; i < N; i++)
    {
        a = (*decoys)[i];
        for (int j=i+1; j < N; j++)
        {
            b = (*decoys)[j];
            r = trueD(a, b);
            if (r < *minDist)
                *minDist = r;
            if (r > *maxDist)
                *maxDist = r;

            nl = nbors[i];
#ifdef _DEBUG_NBORLIST_
            ni = nborsIndex[i];
#endif
            if (r < nl[listLength-1]) // nbors[I] should be altered
            {
                nl[listLength-1] = r;
#ifdef _DEBUG_NBORLIST_
                ni[listLength-1] = j;
#endif
                // double r up the list nl till it's larger
                for (int k=listLength-1; k > 0 && nl[k-1] > r; k--)
                {
                    // swap nl[k] with nk[k-1]
                    nl[k] = nl[k-1];
                    nl[k-1] = r;
#ifdef _DEBUG_NBORLIST_
                    ni[k] = ni[k-1];
                    ni[k-1] = j;
#endif
                }
            }
            nl = nbors[j];
#ifdef _DEBUG_NBORLIST_
            ni = nborsIndex[j];
#endif
            if (r < nl[listLength-1]) // nbors[I] should be altered
            {
                nl[listLength-1] = r;
#ifdef _DEBUG_NBORLIST_
                ni[listLength-1] = i;
#endif
                // double r up the list nl till it's larger
                for (int k=listLength-1; k > 0 && nl[k-1] > r; k--)
                {
                    // swap nl[k] with nk[k-1]
                    nl[k] = nl[k-1];
                    nl[k-1] = r;
#ifdef _DEBUG_NBORLIST_
                    ni[k] = ni[k-1];
                    ni[k-1] = i;
#endif
                }
            }
        }
    }

#ifdef _DEBUG_NBORLIST_
    for (int i=0; i < mNames->size(); i++)
    {
        cout << i << "'" << (*mNames)[i];
        ni = (int *)nborsIndex[i];
        for (int j=0; j < listLength; j++)
        {
           int index = ni[j];
           cout << " " << index << "'" << (*mNames)[index];
        }
        cout << endl;
    }
#endif
    return nbors;
}

static int
__cmp(const void *a, const void *b)
{
    return *((double *)a) >= *((double *)b);
}

/**
 * Find:
 * (1) min, max distances
 * (2) most frequently occuring distance
 * (3) the x percentile distance. using the x percentile
 *     distance is like claiming that there are only x percent of edges
 *     between decoys, and that these edges are the nearest ones.
 * Note that if pairwise distances follows Gaussian distribution then case (2)
 * should be the same as case (3) at x=50%.
 */
void
Clustering::estimateDist(vector<Stru *> * decoys, double xPercent,
                         double *minDist, double *maxDist,
                         double *mostFreqDist, double *xPercentileDist)
{
    //char * dName;
    Stru *a, *b;
    double r;
    int numdecoys = decoys->size();
    int numofpairs = numdecoys * (numdecoys-1) / 2;
    double * alldists = new double[numofpairs];
    int k = 0;
    *maxDist = 0;
    *minDist = _OVER_RMSD_;

    for (int i=0; i < numdecoys; i++)
    {
        a = (*decoys)[i];
        for (int j=i+1; j < numdecoys; j++)
        {
            b = (*decoys)[j];
            r = trueD(a, b);
            alldists[k] = r;
            k++;
            if (r < *minDist && r >= 0) // be paranoid about crazy RMSD
                *minDist = r;
            if (r > *maxDist && r < _OVER_RMSD_) // crazy RMSD
                *maxDist = r;
        }
    }

    // Break the interval [minDist,maxDist] into numofbins bins and find
    // the frequencies of the distances within the range of each bin
    int bin_to_put;
    int numofbins = numdecoys;
    double bin_size = (*maxDist - *minDist)/numofbins;
    int * bins = new int[numofbins];
    for (int i=0; i < numofbins; i++)
        bins[i] = 0;
    for (int i=0; i < numofpairs; i++)
    {
        bin_to_put = (int)floor((alldists[i] - *minDist) / bin_size);
        if (bin_to_put >= numofbins || bin_to_put < 0) // crazy RMSD
            continue;
        bins[bin_to_put]++;
    }

    //cout << *minDist << "," << *maxDist << ", bin size=" << bin_size << endl;
    //for (int i=0; i < numofbins; i++)
    //{
    //    cout << "bin " << i << "," << (*minDist + (i+0.5)*bin_size)
    //         << ": " << bins[i] << endl;
    //}

    // Find the peak of the distribution: Find n consecutive bins which
    // adds up to the highest frequency and return the center bin of these.

    // decide n and the length of one of n's side
    int n = numofbins / 10;
    if (n % 2 == 0) // use odd
        n++;
    int n_half = (n-1)/2;

    int sum_of_n_bins = 0;
    for (int i=0; i < n-1; i++) // add up n-1 bins first
        sum_of_n_bins += bins[i];

    int highest_frequency = 0;
    int bin_of_highest;
    for (int i = n_half; i < numofbins - n_half; i++)
    {
        sum_of_n_bins += bins[i + n_half]; // add the n-th bin
        if (sum_of_n_bins > highest_frequency)
        {
            highest_frequency = sum_of_n_bins;
            bin_of_highest = i;
        }
        sum_of_n_bins -= bins[i - n_half]; // remove the first bin
    }

    *mostFreqDist = (*minDist + (bin_of_highest + 0.5) * bin_size);

    // Find the x percentile distance through qsort
    qsort(alldists, numofpairs, sizeof(double), __cmp);
    *xPercentileDist = alldists[(int)(xPercent/100*numofpairs)];

    if (autoAdjustPercentile) // Use smaller thresholds for larger data sets
    {
        //xPercentile = 1000./sqrt(mNames->size());
        xPercentile = 100./sqrt(sqrt(mNames->size()));
        if (xPercent > MAX_PERCENTILE_FOR_THRESHOLD)
            xPercent = MAX_PERCENTILE_FOR_THRESHOLD;
        if (xPercent < MIN_PERCENTILE_FOR_THRESHOLD)
            xPercent = MIN_PERCENTILE_FOR_THRESHOLD;
        *xPercentileDist = alldists[(int)(xPercent/100*numofpairs)];
    }

    delete [] alldists;
    delete [] bins;
}


/**
 * Find a decoy with the smallest distance x of N nearest neighbors, where
 * minThres <= x <= maxThres.
 *
 * Hence at threshold x, no other cluster has more elements than N.
 *
 * We try for all sizes N in [minSize, maxSize] starting from size.
 * The direction that N is searched depends on the starting size
 * If there is no N within these set sizes such that minThres <= x <= maxThres,
 * we use either minSize or maxSize, depending on which is nearer.
 *
 * In which case,
 *
 * For example, given the case
 *           ----------------------------------maxThres
 *             |                            |
 *           ----------------------------------minThres
 * threshold   |      smallest x      ______|
 *     ^       |         at N _______/      |^x_1
 *     |       |     ________/              |
 *     |       |____/                       |
 *       N = minSize                    N = maxSize
 *     ---> N
 *                           ...we would use x=x_1 and N=maxSize.
 */
double
Clustering::getThreshold(double ** nbors, int numLists, int listLength,
                         int minSize, int size, double minThres, double maxThres)
{
    int t = size;   // current number of nearest neighbors to try
    int last_t = 0; // the last number of nearest neighbors found suitable
    int maxSize = listLength;
    double rmsd, best_rmsd, threshold;

    /** Look for a number t where:
     * 1. minSize <= t <= maxSize
     * 2. minThres <= smallest C(t) among all decoys <= maxThres
     *    where C(t) for a decoy is the max rmsd of the decoy's t nearest neighbors
     *
     * At every stage below, no cluster of best_rmsd radies has size above t.
     */
    while (t <= maxSize && t >= minSize)
    {
        /** smallest C(t) among all decoys */
        best_rmsd = _OVER_RMSD_;
        for (int i=0; i < numLists; i++)
        {
            rmsd = nbors[i][t-1]; // C(t) of decoy t
            if (rmsd < best_rmsd)
                best_rmsd = rmsd;  // smallest C(t) so far
        }

        if (best_rmsd < minThres) // rmsd too small, need to increase target...
        {
            if (last_t == t+1) break; // ...but, we reduced it just now!
            last_t = t;
            t++;
        }
        else if (best_rmsd > maxThres) // rmsd too large, need to decrease t
        {
            if (last_t == t-1) break; // ...but, we increased it just now!
            last_t = t;
            t--;
        }
        else
        {
            cout << "Cluster size " << t << " => threshold=" << best_rmsd
                 << ". [" << minThres << "," << maxThres << "]" << endl;
            break; // rmsd within threshold. Yahoo
        }
    }

    if (t > maxSize) t = maxSize;
    else if (t < minSize) t = minSize;

    if (best_rmsd > maxThres || best_rmsd < minThres)
        cout << "No cluster size results in threshold within (" << minThres
             << "," << maxThres << ")." << endl
             << "Using threshold=" << best_rmsd
             << ". Clusters' size will be at most " << t << endl;

    threshold = best_rmsd;
    assert (t <= maxSize && t >= minSize);

    return threshold;
}


/**
 * For each decoy j, and for each of the decoys i in index[], find and
 * cache dist(i,j). index[] is of size REFERENCE_SIZE, and is initially
 * set to decoys 1, 2, ..., REFERENCE_SIZE-1.
 * This mNumPDBs* distance is stored in mReference
 */
void
Clustering::initRef(int* index)
{
    if (index == NULL)
    {
        index = new int[REFERENCE_SIZE]; // a small memory leak here
        if (Pref::ref_use_first_few_decoys)
        {
            for (int i=0; i < REFERENCE_SIZE; i++)
                index[i] = i;
        }
        else
        {
            int step = mNumPDB/REFERENCE_SIZE-1;
            int x = 0;
            for (int i=0; i < REFERENCE_SIZE; i++)
            {
                index[i] = x;
                x += step;
            }
        }
    }

    //int numPDB = mNames->size();
    mReference = new double[REFERENCE_SIZE*mNumPDB];

    for (int j=0; j < mNumPDB; j++)
    {
        for (int i=0; i < REFERENCE_SIZE; i++)
        {
            double d = trueD(index[i],j);
            mAdjacentList[index[i]]->add(j, d, false);
            mAdjacentList[j]->add(index[i], d, false);
            mReference[j*REFERENCE_SIZE+i] = d; //trueD(index[i],j);
        }
    }
}


void
Clustering::refBound(int i, int j, double& lower, double& upper)
{
    lower = 0;
    upper = _OVER_RMSD_;
    double* first  = mReference+(i*REFERENCE_SIZE);
    double* second = mReference+(j*REFERENCE_SIZE);
    for(int k=0; k < REFERENCE_SIZE; k++)
    {
       double diff = fabs(first[k]-second[k]);
       double sum  = first[k] + second[k];
       if (diff > lower) lower = diff;
       if (sum  < upper) upper = sum;
    }
}

#endif
