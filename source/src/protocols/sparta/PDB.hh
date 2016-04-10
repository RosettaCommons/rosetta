// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
/*                                                                            */
/*                           ----  SPARTA   ----                              */
/*     Shifts Prediction from Analogue of Residue type and Torsion Angle      */
/*                           Yang Shen and Ad Bax                             */
/*                    J. Biomol. NMR, 38, 289-302 (2007)                      */
/*                 NIH, NIDDK, Laboratory of Chemical Physics                 */
/*                     version, 1.01 (build 2009.0928.17)                     */
/*                                                                            */
/*                      for any problem, please contact                       */
/*                          shenyang@niddk.nih.gov                            */
/*                                                                            */
/******************************************************************************/

#ifndef PDB_H
#define PDB_H

#include <vector>
#include <map>
#include <algorithm>
#include <boost/unordered_map.hpp>

#include <utility/vector0.hh>

#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace sparta {

typedef float Vec3[3];
typedef float Mat3[3][3];
//typedef boost::unordered_map<int, boost::unordered_map<int, double> > DoubleD2;
//typedef boost::unordered_map<int, double>  DoubleD1;

struct PDB_Entry
{
	int atomNum, resNum;
	std::string atomName, resName, chainName;
	float X, Y, Z, B_Factor;
	Vec3 Coord;
};

extern std::ostream& operator<<( std::ostream&, PDB_Entry const& );

struct RingData {
	int resID, atomNo;
	std::string resName;
	PDB_Entry atoms[6];

	Vec3 coordA[6];
	Vec3 center, norm;
	Mat3 transM;
	float rad;

	float ringFact;
};

typedef boost::unordered_map<int, PDB_Entry> Mol; //atom list indexed by atom number
typedef boost::unordered_map<int, Mol> Mols;


class PDB
{
	long sgn(float x);
	float arccos_(float x);
	void Vec3Zero(Vec3 v);
	void Vec3Copy(Vec3 v1, Vec3 v2);
	float Vec3Abs(Vec3 v);
	float Vec3DiffAbs(Vec3 v1, Vec3 v2);
	void Vec3Norm(Vec3 v);
	void Vec3Scale(Vec3 v, float s);
	void Vec3Add(Vec3 v1, Vec3 v2);
	void Vec3Sub(Vec3 v1, Vec3 v2);
	float Vec3Scalar(Vec3 v1, Vec3 v2);
	void Vec3Cross(Vec3 v1, Vec3 v2);
	void Mat3VecMult(Vec3 v, Mat3 m);
	void Vec3ScaleAdd(Vec3 v1, float s, Vec3 v2);

public:

	std::string PDBfileName;

	RingData Rings[2000]; //Aromatic ring list
	int RingNo;

	PDB_Entry EMPTY;

	//molecules
	Mols Conformers;
	//atom list indexed by conformer ID + atom number
	typedef boost::unordered_map< std::string, PDB_Entry> AtomEntries;
	typedef boost::unordered_map<int, AtomEntries > EntryMap;
	typedef boost::unordered_map<int, EntryMap > AtomsMap;
	AtomsMap ATOMS;
	//atom list indexed by conformer ID + residue number + atom name

	std::map<int, std::string> residList, residListOne; //resdiue lists with three-letter-aa and one-letter-aa name
	int r1, rN;

	typedef std::map< int, int > PairList;
	PairList acceptorList; // Acceptor List for H-bond indexed by atom index
	PairList donorList; // Donor List for H-bond indexed by atom indeice of itself and its connected heavy atom

	typedef std::map< std::string, float > InnerHBondMap;
	typedef std::map<int, InnerHBondMap > HBondMap;
	HBondMap HBDistList; // HBond distance indexed by residue number and atom name
	HBondMap HBEnergyList; // HBond energy list indexed by residue number and atom name
	HBondMap HB_DHO_AngleList; // HBond angle (Donator-H-O) indexed by residue number and atom name
	HBondMap HB_HOA_AngleList; // HBond angle (H-O-Acceptor_base) indexed by residue number and atom name


	int SpherePointNo;
	Vec3 *SpherePoints;
	int SurfPrec ;
	boost::unordered_map< int, utility::vector0<int> > NeighborList;
	// Neighboring atom list, indexed by atom number of target atom, and utility::vector0 list of neighbor atoms
	boost::unordered_map< int,float> ResSurfaceFullList, ResSurfacePartList;
	// Residue Surface list, indexed by residue number, and surface area for this residue

	boost::unordered_map<int, boost::unordered_map< std::string,float> > AtomSurfaceFullList, AtomSurfacePartList;
	// Atom Surface list, indexed by residue number, atom name, and surface area for this residue
	boost::unordered_map< std::string, float > VDW_RAD; //VDW radius for different type of atoms

	boost::unordered_map<int, float> HN_S2;
	boost::unordered_map<int, boost::unordered_map< std::string, float> > ElectricField; // indexed by resID and atomName

	PDB();
	PDB(const std::string& fileName);

	//amino acid name convertion
	std::string getThreeAAName(char a);
	std::string getOneAAName(const std::string& a);

	//load pdb coordinates
	void loadPDB(const std::string &fileName);
	void loadPDB(std::istream &file);

	void loadPDB( core::pose::Pose const& pose );

	void loadPDB_Entry(const std::string &str, PDB_Entry &entry);
	std::string getField(const std::string &str, int index);

	//get atom
	PDB_Entry getEntry(int conformerID, int rNum, const std::string &aName); //get atom by residue number and atom name
	PDB_Entry getEntry(int conformerID, int aNum); //get atom by atom number

	//get angles and distances
	float getBondAngle(Vec3 A, Vec3 B, Vec3 C);
	float getBondAngle(PDB_Entry a, PDB_Entry b, PDB_Entry c);
	float getDihedralAngle(PDB_Entry a, PDB_Entry b, PDB_Entry c, PDB_Entry d);
	float getPhi(int conformerID, int resNum);
	float getPsi(int conformerID, int resNum);
	float getOmega(int conformerID, int resNum);
	float getChi1(int conformerID, int resNum);
	float getChi2(int conformerID, int resNum);
	float getDist(Vec3 A, Vec3 B);
	float getDist(PDB_Entry A, PDB_Entry B);

	bool isSSBonded(int conformerID, int resNum);

	//get ring current shifts
	void initOrbitalShift();
	float getOrbitalShift(int conformerID, int resNum, const std::string &aName);
	void calcPlane(RingData *ringP);

	//get H-bond information
	void initHBond(float DIST=3.5, float ANGLE=35);
	float getHBondDist(PDB_Entry D);
	float getHBondDist(int resNum, std::string atomName);
	PDB_Entry isAcceptor(PDB_Entry A); //check if an atom is an Acceptor for a H-Hond
	PDB_Entry isDonor(PDB_Entry D); //check if an aotm is a Donor for a H-hond, and return the connected heavy atom

	// void SphereCalcPoints(int divNo, Vec3 **pointAP, int *pointNoP);
	void calcTriangles( double x0, double y0, double z0,
		double x1, double y1, double z1,
		double x2, double y2, double z2,
		int rowStartA[], int rowNo, int quad,
		int row0, int ind0, int ind1,
		int row2, int ind2,
		Vec3 *pointA);
	void findNeighors(float rad_sol);
	//  void initSurface(float rad_sol );
	void calcSurface( float rad_sol );

	void calc_HN_S2( );
	void calc_ElectricField( );

	void collect_HN_S2_and_EF( );
};

}
}

#endif
