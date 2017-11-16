// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
/*                                                                            */
/*                           ----  SPARTA   ----                              */
/*     Shifts Prediction from Analogue of Residue type and Torsion Angle      */
/*                           Yang Shen and Ad Bax                             */
/*                    J. Biomol. NMR, 38, 289-302 (2007)                      */
/*                 NIH, NIDDK, Laboratory of Chemical Physics                 */
/*                     version, 1.00 (build 2010.0607.00)                     */
/*                                                                            */
/*                      for any problem, please contact                       */
/*                          shenyang@niddk.nih.gov                            */
/*                                                                            */
/******************************************************************************/

// 01/25/2007, fix the probelm of missed ring atoms for ring current shifts calculation
#include <boost/unordered_map.hpp>
#include <protocols/sparta/PDB.hh>
#include <fstream>
#include <sstream>
#include <core/pose/Pose.hh>
// Utility headers
#include <basic/Tracer.hh>
#include <protocols/sparta/util.hh>
#include <protocols/sparta/constants.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <numeric/NumericTraits.hh>

namespace protocols {
namespace sparta {

static basic::Tracer tr( "protocols.sparta" );
using namespace core;

using namespace std;

std::ostream &
operator<<( std::ostream & os, PDB_Entry const& pe )
{
	os << "PDBEntry: " << pe.atomNum << " " << pe.resNum << " " << pe.atomName << " " << pe.resName;
	return os;
}

PDB::PDB()
{
	r1 = 9999; rN = -9999;
}


PDB::PDB(const string &fileName)
{
	PDBfileName = fileName;

	r1 = 9999; rN = -9999;

	loadPDB(fileName);
}


string PDB::getField(const string &str, int index)
{
	switch(index)
			{
			case 1 : return str.substr(0,6); //"ATOM"
			case 2 : return str.substr(6,5); //Atom number
			case 3 : return str.substr(11,5); //Atom name
			case 4 : return str.substr(17,4); //Residue name
			case 5 : return str.substr(21,1); //Chain ID
			case 6 : return str.substr(22,4); //Reside seq
			case 7 : return str.substr(30,8); //X
			case 8 : return str.substr(38,8); //Y
			case 9 : return str.substr(46,8); //Z
			case 10 : return str.substr(54,6); //Occupancy
			case 11 : return str.substr(60,6); // B-factor
			}

	return nullptr;
}


void PDB::loadPDB_Entry(const string &str, PDB_Entry &entry)
{
	entry.atomNum = atoi( getField(str,2).c_str() );

	string atomName = simplifyWhiteSpace( getField(str,3).c_str() );
	entry.atomName = atomName;
	if ( atomName == "H" ) entry.atomName="HN";
	if ( atomName[0] >= '1' && atomName[0] <= '3' ) entry.atomName = atomName.substr(1,4) + atomName[0];

	// remove the charge character '+'
	entry.resName =  simplifyWhiteSpace( getField(str,4).substr(0,3).c_str() );
	entry.chainName = getField(str,5);
	entry.resNum = atoi( getField(str,6).c_str() );
	entry.X = atof( getField(str,7).c_str() );
	entry.Y = atof( getField(str,8).c_str() );
	entry.Z = atof( getField(str,9).c_str() );
	entry.B_Factor = atof( getField(str,11).c_str() );

	entry.Coord[0] = entry.X;
	entry.Coord[1] = entry.Y;
	entry.Coord[2] = entry.Z;
}

void PDB::loadPDB(const string &fileName)
{
	PDBfileName = fileName;

	ifstream file(fileName.c_str());
	if ( ! file.is_open() ) {
		tr.Error << "\tCan't open file " << fileName << " for reading" << endl;
		exit(0);
	}

	loadPDB( file );
}

void PDB::loadPDB( core::pose::Pose const& pose ) {
	std::stringstream obuffer;
	pose.dump_pdb( obuffer);
	std::stringstream ibuffer( obuffer.str() );
	loadPDB( ibuffer );
}

void PDB::loadPDB(istream &file) {

	residList.clear();
	residListOne.clear();
	ATOMS.clear();
	Conformers.clear();

	acceptorList.clear();
	donorList.clear();
	HBDistList.clear();

	int residue_offset = 0;
	Mol conformer;
	int conformer_count=1;
	int lastline_index = 0;

	while ( !file.eof() ) {
		char buf[1000];
		file.getline(buf, 1000);
		string str(buf);
		if ( str.empty() || str.size() == 0 ) continue;
		int pos = str.find("ATOM"); //starting of one molecule
		int pos2 = str.find("TER"), pos3 = str.find("END"); //ending of one molecule
		if ( pos != 0 && pos2 != 0 && pos3 != 0 ) continue;

		if ( pos2 == 0 || pos3 == 0 ) {     // this part has never been executed
			Conformers[conformer_count] = conformer;
			conformer_count++;
			conformer.clear();
			break;
			//   continue;
		}

		PDB_Entry temp;
		loadPDB_Entry(str, temp);

		if ( lastline_index>1 && temp.chainName != conformer[lastline_index].chainName ) residue_offset = conformer[lastline_index].resNum-temp.resNum+1;
		//to modify the resNum if  more than one chain
		temp.resNum = temp.resNum + residue_offset;

		conformer[temp.atomNum] = temp; // table indexed by the atom number
		lastline_index = temp.atomNum;
		//  tr.Trace << "resNum : "<< temp.resNum << std::endl;
		//  std::cout << "resNum: " << conformer[lastline_index].resNum<<std::endl;
		ATOMS[conformer_count][temp.resNum][temp.atomName] = temp;
		// table indexed by the residue number and atom name
		// fast for obtaining a specific atom with given residue number and atom name

		residList[temp.resNum] = temp.resName; //resdiue list with three-letter-aa name
		residListOne[temp.resNum] = getOneAAName(temp.resName); //resdiue list with one-letter-aa name

		//   std::cout << "resNum: to check r1 rN: " << temp.resNum << " r1: " << r1 << " rN: " << rN << std::endl;
		if ( temp.resNum < r1 ) r1 = temp.resNum;
		if ( temp.resNum > rN ) rN = temp.resNum;
	}

	if ( conformer.size() != 0 ) Conformers[conformer_count] = conformer;

	//correction of the name cys/CYS accoding to the status of the disulfide bond
	std::map<int, string> ::iterator it, end;
	for ( it=residListOne.begin(), end =residListOne.end(); it != end; ++it ) {
		if ( it->second == "C" ) {
			if ( isSSBonded(1, it->first) ) {
				residListOne[it->first] = "c";
				residList[it->first] = "cys";
			}
		}
	}
}


bool PDB::isSSBonded(int /*conformerID*/, int resNum)
{
	if ( residListOne[resNum] != "C" && residListOne[resNum] != "c" ) return false;

	PDB_Entry CYS_SG = getEntry(1,resNum, "SG");

	std::map<int, string>::iterator it, end;
	for ( it = residListOne.begin(), end = residListOne.end(); it != end; ++it ) {
		if ( (it->second == "C" || it->second == "c") && std::abs(it->first - resNum) >= 4 ) {
			if ( getDist(CYS_SG.Coord, ATOMS[1][it->first]["SG"].Coord ) <= 2.5 ) {
				return true;
			}
		}
	}

	return false;
}


PDB_Entry PDB::getEntry(int conformerID, int rNum, const string &aName)
{
	return ATOMS[conformerID][rNum][aName];
}


PDB_Entry PDB::getEntry(int conformerID, int aNum)
{
	return Conformers[conformerID][aNum];
}


float PDB::getBondAngle(Vec3 A, Vec3 B, Vec3 C)
{
	Vec3 v,v1,v2;

	Vec3Copy(v,B); Vec3Copy(v1,A); Vec3Copy(v2,C);
	Vec3Sub(v1,v); // vector BA
	Vec3Sub(v2,v); // vector BC

	Vec3Norm(v1); // unit vector BA
	Vec3Norm(v2); // unit vector BC

	float c = Vec3Scalar(v1, v2); // scalar products
	Vec3Cross(v1, v2); // cross product
	float s = Vec3Abs(v1);

	return atan2f(s, c)*180.0/numeric::NumericTraits<Real>::pi();
}


float PDB::getBondAngle(PDB_Entry a, PDB_Entry b, PDB_Entry c)
{
	return getBondAngle(a.Coord, b.Coord, c.Coord);
}


float PDB::getDihedralAngle(PDB_Entry a, PDB_Entry b, PDB_Entry c, PDB_Entry d)
{
	double cb[3], n1[3], n2[3];
	float co;
	double TEMP, TEMP1, TEMP2, TEMP3, TEMP4, TEMP5;

	if ( (a.atomName).empty() || (b.atomName).empty() || (c.atomName).empty() || (d.atomName).empty() ) {
		return SPARTA_MAXNUM;
	}

	cb[0] = c.X - b.X; cb[1] = c.Y - b.Y; cb[2] = c.Z - b.Z; // vector b->c
	n1[0] = (b.Y - a.Y) * cb[2] + (a.Z - b.Z) * cb[1]; //normal vector of plane a-b-c
	n1[1] = (b.Z - a.Z) * cb[0] + (a.X - b.X) * cb[2];
	n1[2] = (b.X - a.X) * cb[1] + (a.Y - b.Y) * cb[0];
	n2[0] = cb[1] * (d.Z - c.Z) + cb[2] * (c.Y - d.Y); //normal vector of plane b-c-d
	n2[1] = cb[2] * (d.X - c.X) + cb[0] * (c.Z - d.Z);
	n2[2] = cb[0] * (d.Y - c.Y) + cb[1] * (c.X - d.X);
	TEMP  = n1[0]; TEMP1 = n1[1]; TEMP2 = n1[2];
	TEMP3 = n2[0]; TEMP4 = n2[1]; TEMP5 = n2[2];
	co = (float)((n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2]) /
		sqrt( (TEMP*TEMP + TEMP1*TEMP1 + TEMP2*TEMP2)*(TEMP3*TEMP3 + TEMP4*TEMP4 + TEMP5*TEMP5) ) );
	// cos(a-b-c-d) = |n1.n2|/|n1|.|n2|

	return 180./numeric::NumericTraits<Real>::pi()*(arccos_(co)*
		sgn( (float)( (n1[1]*n2[2] - n1[2]*n2[1])*cb[0] + (n1[2]*n2[0] - n1[0]*n2[2])*cb[1] + (n1[0]*n2[1] - n1[1]*n2[0])*cb[2]) ) );
}


float PDB::getPhi(int conformerID, int resNum)
{
	return getDihedralAngle(getEntry(conformerID,resNum-1,"C"), \
		getEntry(conformerID,resNum,"N"), \
		getEntry(conformerID,resNum,"CA"), \
		getEntry(conformerID,resNum,"C") );
}


float PDB::getPsi(int conformerID, int resNum)
{
	return getDihedralAngle(getEntry(conformerID,resNum,"N"), \
		getEntry(conformerID,resNum,"CA"), \
		getEntry(conformerID,resNum,"C"), \
		getEntry(conformerID,resNum+1,"N") );
}


float PDB::getOmega(int conformerID, int resNum)
{
	return getDihedralAngle(getEntry(conformerID,resNum-1,"CA"), \
		getEntry(conformerID,resNum-1,"C"), \
		getEntry(conformerID,resNum,"N"), \
		getEntry(conformerID,resNum,"CA"));
}


float PDB::getChi1(int conformerID, int resNum)
{
	string rName = residList[resNum];

	PDB_Entry R;
	if ( rName == "GLY" || rName == "ALA" ) {
		return SPARTA_MAXNUM;
	} else if ( rName == "VAL" ) {
		R = getEntry(conformerID,resNum,"CG2");
	} else if ( rName == "ILE" ) {
		R = getEntry(conformerID,resNum,"CG1");
	} else if ( rName == "THR" ) {
		R = getEntry(conformerID,resNum,"OG1");
	} else if ( rName == "SER" ) {
		R = getEntry(conformerID,resNum,"OG" );
	} else if ( rName == "CYS" || rName == "cys" ) {
		R = getEntry(conformerID,resNum,"SG" );
	} else {
		R = getEntry(conformerID,resNum,"CG" );
	}

	return getDihedralAngle(getEntry(conformerID,resNum,"N"), \
		getEntry(conformerID,resNum,"CA"), \
		getEntry(conformerID,resNum,"CB"), R);
}


float PDB::getChi2(int conformerID, int resNum)
{
	string rName = residList[resNum];
	if ( rName.compare("GLY") == 0 || rName.compare("ALA") == 0 ) {
		return SPARTA_MAXNUM;
	}

	PDB_Entry R3, R4;

	if ( rName == "VAL" ) { //chi21
		R3 = getEntry(conformerID,resNum,"CG1");
		R4 = getEntry(conformerID,resNum,"HG11");
	} else if ( rName == "THR" ) {
		R3 = getEntry(conformerID,resNum,"CG2");
		R4 = getEntry(conformerID,resNum,"HG21");
	} else if ( rName == "SER" ) {
		R3 = getEntry(conformerID,resNum,"OG");
		R4 = getEntry(conformerID,resNum,"HG");
	} else if ( rName == "ILE" ) { //chi21
		R3 = getEntry(conformerID,resNum,"CG1" );
		R4 = getEntry(conformerID,resNum,"CD1" );
	} else if ( rName == "CYS" || rName == "cys" ) {
		R3 = getEntry(conformerID,resNum,"SG" );
		R4 = getEntry(conformerID,resNum,"HG" );
	} else {
		R3 = getEntry(conformerID,resNum,"CG" );
	}

	if ( rName == "ASN" || rName == "ASP" ) R4 = getEntry(conformerID,resNum,"OD1" );
	else if ( rName == "HIS" ) R4 = getEntry(conformerID,resNum,"ND1" );
	else if ( rName == "MET" ) R4 = getEntry(conformerID,resNum,"SD" );
	else if ( rName == "LEU" || rName == "PHE" || rName == "TYR" || rName == "TRP" ) {
		R4 = getEntry(conformerID,resNum,"CD1" );
	} else R4 = getEntry(conformerID,resNum,"CD" );

	return getDihedralAngle(getEntry(conformerID,resNum,"CA"),getEntry(conformerID,resNum,"CB"),R3, R4);
}


string PDB::getThreeAAName(char a)
{
	switch(a)
			{
			case 'A' : return "ALA";
			case 'C' : return "CYS";
			case 'c' : return "cys";
			case 'D' : return "ASP";
			case 'E' : return "GLU";
			case 'F' : return "PHE";
			case 'G' : return "GLY";
			case 'H' : return "HIS";
			case 'I' : return "ILE";
			case 'K' : return "LYS";
			case 'L' : return "LEU";
			case 'M' : return "MET";
			case 'N' : return "ASN";
			case 'P' : return "PRO";
			case 'Q' : return "GLN";
			case 'R' : return "ARG";
			case 'S' : return "SER";
			case 'T' : return "THR";
			case 'V' : return "VAL";
			case 'W' : return "TRP";
			case 'Y' : return "TYR";
			}

	return "???";
}


string PDB::getOneAAName(const string& a)
{
	if ( a == "ALA" ) return "A";
	else if ( a == "CYS" ) return "C";
	else if ( a == "cys" ) return "c";
	else if ( a == "ASP" ) return "D";
	else if ( a == "GLU" ) return "E";
	else if ( a == "PHE" ) return "F";
	else if ( a == "GLY" ) return "G";
	else if ( a == "HIS" || a == "HIS+" ) return "H";
	else if ( a == "ILE" ) return "I";
	else if ( a == "LYS" || a == "LYS+" ) return "K";
	else if ( a == "LEU" ) return "L";
	else if ( a == "MET" ) return "M";
	else if ( a == "ASN" ) return "N";
	else if ( a == "PRO" ) return "P";
	else if ( a == "GLN" ) return "Q";
	else if ( a == "ARG" || a == "ARG+" ) return "R";
	else if ( a == "SER" ) return "S";
	else if ( a == "THR" ) return "T";
	else if ( a == "VAL" ) return "V";
	else if ( a == "TRP" ) return "W";
	else if ( a == "TYR" ) return "Y";

	return "?";
}


//collet the ring information required for calculating Haigh-Mallion ring shifts
//Haigh, C.W. and Mallion, R.B. (1979) Progr. NMR Spectrosc., 13, 303-344
//Parameters are taken from Case, D.A. (1995) J. Biomol.NMR, 6, 341-346.
void PDB::initOrbitalShift(){

	Mol conf = Conformers[1];

	string PF6Atoms[6] = {"CG", "CD1", "CE1", "CZ", "CE2", "CD2" };
	string W6Atoms[6] = {"CD2", "CE2", "CZ2", "CH2", "CZ3", "CE3" };
	string H5Atoms[5] = {"CG", "ND1", "CE1", "NE2", "CD2" };
	string W5Atoms[5] = {"CG", "CD1", "NE1", "CE2", "CD2" };

	std::map<int, string>::iterator it, end;
	Mol::iterator itA;

	int i;
	RingNo = 0;

	for ( it = residList.begin(), end = residList.end(); it != end; ++it ) {
		int resID = it->first; //.resNum;
		string resName = it->second;

		if ( resName == "PHE" || resName == "TYR" || resName == "TRP" ) { // 6-member ring
			Rings[RingNo].resID = resID; Rings[RingNo].resName = resName; Rings[RingNo].atomNo = 6;
			if ( resName == "PHE" ) Rings[RingNo].ringFact = 1.46;
			if ( resName == "TYR" ) Rings[RingNo].ringFact = 1.24;
			if ( resName == "TRP" ) Rings[RingNo].ringFact = 1.24;

			bool NO_MISSED_RING_ATOMS = true;

			for ( i = 0; i < 6; i++ ) {
				PDB_Entry atom;
				if ( resName == "PHE" || resName == "TYR" ) {
					atom = getEntry(1,resID,PF6Atoms[i]);
					if ( atom.atomName == PF6Atoms[i] ) {
						Vec3Copy(Rings[RingNo].coordA[i], atom.Coord);
					} else NO_MISSED_RING_ATOMS = false; // missed ring atom in coordinates file
				} else if ( resName == "TRP" ) {
					atom = getEntry(1,resID,W6Atoms[i]);
					if ( atom.atomName == W6Atoms[i] ) {
						Vec3Copy(Rings[RingNo].coordA[i], atom.Coord);
					} else NO_MISSED_RING_ATOMS = false; // missed ring atom in coordinates file
				}
			}

			if ( NO_MISSED_RING_ATOMS ) RingNo++;
		}

		if ( resName == "HIS" || resName == "TRP" ) {
			// 5-member ring
			Rings[RingNo].resID = resID; Rings[RingNo].resName = resName; Rings[RingNo].atomNo = 5;
			if ( resName == "HIS" ) Rings[RingNo].ringFact = 1.35;
			if ( resName == "TRP" ) Rings[RingNo].ringFact = 1.32;

			bool NO_MISSED_RING_ATOMS = true;

			for ( i = 0; i < 5; i++ ) {
				PDB_Entry atom;
				if ( resName == "HIS" ) {
					atom = getEntry(1,resID,H5Atoms[i]);
					if ( atom.atomName == H5Atoms[i] ) {
						Vec3Copy(Rings[RingNo].coordA[i], atom.Coord);
					} else NO_MISSED_RING_ATOMS = false; // missed ring atom in coordinates file
				} else if ( resName == "TRP" ) {
					atom = getEntry(1,resID,W5Atoms[i]);
					if ( atom.atomName == W5Atoms[i] ) {
						Vec3Copy(Rings[RingNo].coordA[i], atom.Coord);
					} else NO_MISSED_RING_ATOMS = false; // missed ring atom in coordinates file
				}
			}

			if ( NO_MISSED_RING_ATOMS ) RingNo++;
		}
	}

	for ( i=0; i < RingNo; i++ ) {
		calcPlane(&Rings[i]);
	}
}


//Calculate the parameters for a ring
void PDB::calcPlane(RingData *ringP)
{
	Vec3 v1, v2;
	int atomI;

	Vec3Zero(ringP->center);
	for ( atomI = 0; atomI < ringP->atomNo; atomI++ ) {
		Vec3Add(ringP->center, ringP->coordA[atomI]);
	}
	Vec3Scale(ringP->center, 1.0f/ringP->atomNo);

	Vec3Zero(ringP->norm);
	ringP->rad = 0.0f;
	Vec3Copy(v1, ringP->coordA[ringP->atomNo - 1]);
	Vec3Sub(v1, ringP->center);
	for ( atomI = 0; atomI < ringP->atomNo; atomI++ ) {
		Vec3Copy(v2, ringP->coordA[atomI]);
		Vec3Sub(v2, ringP->center);
		ringP->rad += Vec3Abs(v2);
		Vec3Cross(v1, v2);
		Vec3Add(ringP->norm, v1);
		Vec3Copy(v1, v2);
	}
	Vec3Norm(ringP->norm);
	ringP->rad /= ringP->atomNo;

	/* transformation matrix ring plane -> x-y plane */
	ringP->transM[0][2] = ringP->norm[0];
	ringP->transM[1][2] = ringP->norm[1];
	ringP->transM[2][2] = ringP->norm[2];

	/* vector orthgonal to norm */
	if ( ringP->norm[0] > 0.5f || ringP->norm[0] < -0.5f ) {
		v1[0] = - ringP->norm[1];
		v1[1] = ringP->norm[0];
		v1[2] = 0.0f;
	} else {
		v1[0] = 0.0f;
		v1[1] = - ringP->norm[2];
		v1[2] = ringP->norm[1];
	}
	Vec3Norm(v1);
	ringP->transM[0][1] = v1[0];
	ringP->transM[1][1] = v1[1];
	ringP->transM[2][1] = v1[2];

	Vec3Cross(v1, ringP->norm);
	ringP->transM[0][0] = v1[0];
	ringP->transM[1][0] = v1[1];
	ringP->transM[2][0] = v1[2];
}


float PDB::getOrbitalShift(int conformerID, int resNum, const string &aName)
{
	Vec3 v;

	Vec3Copy( v, getEntry(conformerID,resNum,aName).Coord );

	//cout << RingNo << " " << resNum << " " << aName << endl;

	float ringShiftSum = 0.0f;
	//loop over all rings in protein
	for ( int i=0; i < RingNo; i++ ) {
		int atomNo, atomI;
		Vec3 vt, v1, v2;
		float area, r1, r2, g, b;

		//exclude the effects from aromatic rings themselves (except for HN atom)
		if ( Rings[i].resID == resNum && aName != "HN" ) continue;

		atomNo = Rings[i].atomNo;
		g = 0.0f;

		//loop over all 6(or 5) atoms in a ring
		for ( atomI = 0; atomI < atomNo; atomI++ ) {
			Vec3Copy(v1, Rings[i].coordA[atomI]);
			Vec3Copy(v2, Rings[i].coordA[(atomI + 1) % atomNo]);

			r1 = Vec3DiffAbs(v1, v);
			r2 = Vec3DiffAbs(v2, v);

			Vec3Copy(vt, v);
			Vec3Sub(vt, Rings[i].center);
			Vec3Sub(v1, Rings[i].center);
			Vec3Sub(v2, Rings[i].center);

			Mat3VecMult(vt, Rings[i].transM);
			Mat3VecMult(v1, Rings[i].transM);
			Mat3VecMult(v2, Rings[i].transM);
			Vec3Sub(v1, vt);
			Vec3Sub(v2, vt);

			area = v1[0] * v2[1] - v1[1] * v2[0];
			g += area * (1.0f / (r1 * r1 * r1) + 1.0f / (r2 * r2 * r2));
			//cout << atomI << " " << g << " ";
		}

		g *= 0.5f;
		b = (float) 5.4548e-6;

		ringShiftSum += Rings[i].ringFact * b * g;

		//if( Rings[i].ringFact * b * g *((float)-1.0e6) > 0.6)
		// cout << ringShiftSum << " " << g << " " << Rings[i].resID  << resNum << " " << aName << endl;
	}

	//cout << ringShiftSum*((float)-1.0e6) << " " << resNum << " " << aName << endl;

	return ringShiftSum*((float)-1.0e6);
}


float PDB::getDist(PDB_Entry A, PDB_Entry B)
{
	return getDist(A.Coord, B.Coord);
}


float PDB::getDist(Vec3 A, Vec3 B)
{
	Vec3 v;

	Vec3Copy(v, B);
	Vec3Sub(v, A);

	return Vec3Abs(v);
}


void PDB::initHBond(float /*DIST*/, float /*ANGLE*/)
{
	Mol conf = Conformers[1];
	Mol::iterator it, end;

	acceptorList.clear();
	donorList.clear();
	HBDistList.clear();

	//retreive the donor and acceptor list from coordinates
	for ( it = conf.begin(), end = conf.end(); it != end; ++it ) {
		PDB_Entry a = isAcceptor(it->second);

		if ( !(a.atomName.empty()) ) acceptorList[it->first] = a.atomNum;
		else {
			PDB_Entry temp = isDonor(it->second);
			if ( temp.atomName.empty() ) continue;
			donorList[it->first] = temp.atomNum;
		}
	}


	for ( PairList::const_iterator itA = acceptorList.begin(), endA = acceptorList.end(); itA != endA; ++itA ) {//loop over acceptor list
		PDB_Entry A = conf[ itA->first ], A_Heavy = conf[ itA->second ];
		//search for donors
		for ( PairList::const_iterator itD = donorList.begin(), endD = donorList.end(); itD != endD; ++itD ) {
			PDB_Entry D = conf[ itD->first ], D_Heavy = conf[ itD->second ];
			if ( std::abs( A.resNum - D.resNum) < 2 ) continue; //minimal 2 residues apart

			float D_ON = getDist(A,D_Heavy);
			float D_CH = getDist(A_Heavy,D);
			float D_OH = getDist(A,D);
			float D_CN = getDist(A_Heavy,D_Heavy);

			if ( D_OH > D_CN ) continue;

			float HBond_E = 332.0*0.42*0.20;
			HBond_E *= ( 1.0/D_ON+1.0/D_CH-1.0/D_OH-1.0/D_CN );
			//Kabsch, W. and Sander, C. (1983) Biopolymer, 22, 2577-2637.

			if ( HBond_E < -0.5 ) { //is a HB bond
				if ( HBDistList[A.resNum][A.atomName] == 0
						|| HBDistList[A.resNum][A.atomName] > D_OH
						|| HBEnergyList[A.resNum][A.atomName] > HBond_E ) {
					HBDistList[A.resNum][A.atomName] = D_OH;
					HB_DHO_AngleList[A.resNum][A.atomName] = getBondAngle(D_Heavy,D,A);
					HB_HOA_AngleList[A.resNum][A.atomName] = getBondAngle(D,A,A_Heavy);
					HBEnergyList[A.resNum][A.atomName] = HBond_E;
				}
				// keep the strongest Hbond

				if ( HBDistList[D.resNum][D.atomName] == 0
						|| HBDistList[D.resNum][D.atomName] > D_OH
						|| HBEnergyList[D.resNum][D.atomName] > HBond_E ) {
					HBDistList[D.resNum][D.atomName] = D_OH;
					HB_DHO_AngleList[D.resNum][D.atomName] = getBondAngle(D_Heavy,D,A);
					HB_HOA_AngleList[D.resNum][D.atomName] = getBondAngle(D,A,A_Heavy);
					HBEnergyList[D.resNum][D.atomName] = HBond_E;
				}
				// keep the strongest Hbond
			}
		}
	}

}


// must run initHBond() first
float PDB::getHBondDist(PDB_Entry D)
{
	return HBDistList[D.resNum][D.atomName];
}


// must run initHBond() first
float PDB::getHBondDist(int resNum, string atomName)
{
	return HBDistList[resNum][atomName];
}


PDB_Entry PDB::isAcceptor(PDB_Entry A)
{
	if ( A.atomName == "O" || A.atomName == "OT1"  || A.atomName == "OT2" ) return ATOMS[1][A.resNum]["C"];

	if ( A.resName == "ASN" && A.atomName == "OD1" ) return ATOMS[1][A.resNum]["CG"];
	if ( A.resName == "ASP" && (A.atomName == "OD1" || A.atomName == "OD2") ) return ATOMS[1][A.resNum]["CG"];
	if ( (A.resName == "CYS"||A.resName == "cys") && A.atomName == "SG" ) return ATOMS[1][A.resNum]["CB"];
	if ( A.resName == "GLN" && A.atomName == "OE1" ) return ATOMS[1][A.resNum]["CD"];
	if ( A.resName == "GLU" && (A.atomName == "OE1" || A.atomName == "OE2") ) return ATOMS[1][A.resNum]["CD"];
	if ( A.resName == "HIS" && (A.atomName == "ND1" || A.atomName == "NE2") ) return ATOMS[1][A.resNum]["CE1"]; // could be CG or CD2
	if ( A.resName == "MET" && A.atomName == "SD" ) return ATOMS[1][A.resNum]["CG"];
	if ( A.resName == "SER" && A.atomName == "OG" ) return ATOMS[1][A.resNum]["CB"];
	if ( A.resName == "THR" && A.atomName == "OG1" ) return ATOMS[1][A.resNum]["CB"];
	if ( A.resName == "TYR" && A.atomName == "OH" ) return ATOMS[1][A.resNum]["CZ"];

	return EMPTY;
}


PDB_Entry PDB::isDonor(PDB_Entry D)
{
	if ( D.atomName == "HN" ) return ATOMS[1][D.resNum]["N"];

	if ( D.atomName.find("HA") != string::npos ) return ATOMS[1][D.resNum]["CA"];

	if ( D.resName == "ARG" ) {
		if ( D.atomName == "HE" ) return ATOMS[1][D.resNum]["NE"];
		if ( D.atomName == "HH11" || D.atomName == "HH12" ) return ATOMS[1][D.resNum]["NH1"];
		if ( D.atomName == "HH21" || D.atomName == "HH22" ) return ATOMS[1][D.resNum]["NH2"];
	}
	if ( D.resName == "ASP" && (D.atomName == "HD21" || D.atomName == "HD22") ) return ATOMS[1][D.resNum]["ND2"];
	if ( (D.resName == "CYS"||D.resName == "cys") && D.atomName == "HG" ) return ATOMS[1][D.resNum]["SG"];
	if ( D.resName == "GLU" && (D.atomName == "HE21" || D.atomName == "HE22") ) return ATOMS[1][D.resNum]["NE2"];
	if ( D.resName == "HIS" ) {
		if ( D.atomName == "HD1" ) return ATOMS[1][D.resNum]["ND1"];
		if ( D.atomName == "HE2" ) return ATOMS[1][D.resNum]["NE2"];
	}
	if ( D.resName == "LYS" && (D.atomName == "HZ1" || D.atomName == "HZ2" || D.atomName == "HZ3") ) return ATOMS[1][D.resNum]["NZ"];
	if ( D.resName == "SER" && D.atomName == "HG" ) return ATOMS[1][D.resNum]["OG"];
	if ( D.resName == "THR" && D.atomName == "HG1" ) return ATOMS[1][D.resNum]["OG1"];
	if ( D.resName == "TRP" && D.atomName == "HE1" ) return ATOMS[1][D.resNum]["NE1"];
	if ( D.resName == "TYR" && D.atomName == "HH" ) return ATOMS[1][D.resNum]["OH"];

	return EMPTY;
}


long PDB::sgn(float x)
{
	return ((x >= 0.0) - (x < 0.0));
}


float PDB::arccos_(float x)
{
	const Real SPARTA_PI = numeric::NumericTraits<Real>::pi();
	if ( x > 1.0 ) x = 1.0;
	else if ( x < -1.0 ) x = -1.0;

	if ( fabs(x) >= 0.5 ) {
		if ( x > 0.0 ) {
			return atan(sqrt(1.0 - x * x) / x);
		} else {
			return (SPARTA_PI + atan(sqrt(1.0 - x * x) / x));
		}
	} else {
		return (0.5 * SPARTA_PI - atan(x / sqrt(1.0 - x * x)));
	}
}


void PDB::Vec3Zero(Vec3 v)
{
	for ( int i = 0; i < 3; i++ ) {
		v[i] = 0.0f;
	}
}


void PDB::Vec3Copy(Vec3 v1, Vec3 v2)
{
	for ( int i = 0; i < 3; i++ ) {
		v1[i] = v2[i];
	}
}


float PDB::Vec3Abs(Vec3 v)
{
	return sqrtf(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}


float PDB::Vec3DiffAbs(Vec3 v1, Vec3 v2)
{
	return sqrtf((v1[0] - v2[0]) * (v1[0] - v2[0]) +
		(v1[1] - v2[1]) * (v1[1] - v2[1]) +
		(v1[2] - v2[2]) * (v1[2] - v2[2]));
}


void PDB::Vec3Norm(Vec3 v)
{
	float a = Vec3Abs(v);
	for ( int i = 0; i < 3; i++ ) {
		v[i] /= a;
	}
}


void PDB::Vec3Scale(Vec3 v, float s)
{
	for ( int i = 0; i < 3; i++ ) {
		v[i] *= s;
	}
}


void PDB::Vec3Add(Vec3 v1, Vec3 v2)
{
	for ( int i = 0; i < 3; i++ ) {
		v1[i] += v2[i];
	}
}


void PDB::Vec3Sub(Vec3 v1, Vec3 v2)
{
	for ( int i = 0; i < 3; i++ ) {
		v1[i] -= v2[i];
	}
}


float PDB::Vec3Scalar(Vec3 v1, Vec3 v2)
{
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}


void PDB::Vec3Cross(Vec3 v1, Vec3 v2)
{
	Vec3 vRes;

	vRes[0] = v1[1] * v2[2] - v1[2] * v2[1];
	vRes[1] = v1[2] * v2[0] - v1[0] * v2[2];
	vRes[2] = v1[0] * v2[1] - v1[1] * v2[0];
	Vec3Copy(v1, vRes);
}


void PDB::Mat3VecMult(Vec3 v, Mat3 m)
{
	Vec3 vRes;

	vRes[0] = v[0] * m[0][0] + v[1] * m[1][0] + v[2] * m[2][0];
	vRes[1] = v[0] * m[0][1] + v[1] * m[1][1] + v[2] * m[2][1];
	vRes[2] = v[0] * m[0][2] + v[1] * m[1][2] + v[2] * m[2][2];

	for ( int i = 0; i < 3; i++ ) {
		v[i] = vRes[i];
	}
}


void PDB::Vec3ScaleAdd(Vec3 v1, float s, Vec3 v2)
{
	int i;

	for ( i = 0; i < 3; i++ ) {
		v1[i] += s * v2[i];
	}
}


// void PDB::initSurface( float rad_sol )
// {

//   float SolventRad = 1.4f;


//   VDW_RAD["H"] = 1.2; VDW_RAD["HN"] = 1.2;

//   VDW_RAD["N"] = 1.55; VDW_RAD["ND1"] = 1.55; VDW_RAD["ND2"] = 1.55; VDW_RAD["NE"] = 1.55; VDW_RAD["NE1"] = 1.55; VDW_RAD["NE2"] = 1.55;
//   VDW_RAD["NH1"] = 1.55; VDW_RAD["NH2"] = 1.55; VDW_RAD["NZ"] = 1.55;

//   VDW_RAD["C"]  = 2.3; //1.1 C-H bond-length + 1.2 H radius
//   VDW_RAD["CO"] = 1.7;
//   //VDW_RAD["C"] = 1.7;  VDW_RAD["CA"] = 1.7;  VDW_RAD["CB"] = 1.7;  VDW_RAD["CG"] = 1.7; VDW_RAD["CG1"] = 1.7; VDW_RAD["CG2"] = 1.7;
//   //VDW_RAD["CD"] = 1.7; VDW_RAD["CD1"] = 1.7; VDW_RAD["CD2"] = 1.7; VDW_RAD["CE"] = 1.7; VDW_RAD["CE1"] = 1.7; VDW_RAD["CE2"] = 1.7;
//   //VDW_RAD["CZ"] = 1.7; VDW_RAD["CZ2"] = 1.7; VDW_RAD["CZ3"] = 1.7; VDW_RAD["CH2"] = 1.7;

//   VDW_RAD["O"] = 1.52;  VDW_RAD["OD1"] = 1.52; VDW_RAD["OD2"] = 1.52;
//   VDW_RAD["OG"] = 1.52; VDW_RAD["OG1"] = 1.52; VDW_RAD["OE1"] = 1.52; VDW_RAD["OE2"] = 1.52; VDW_RAD["OH"] = 1.52;

//   VDW_RAD["SG"] = 1.8; VDW_RAD["SD"] = 1.8;

//   SurfPrec = 3;
//   SpherePointNo = 0;


//   clock_t start, finish;
//   start = clock();

//   //cout << "calculate surface points " << endl;
//   SphereCalcPoints(SurfPrec, &SpherePoints, &SpherePointNo);
//   finish = clock();
//   tr.Info << "\t SphereCalcPoints() running time: " << (float)(finish - start)/ CLOCKS_PER_SEC << " seconds" << endl;

//   //cout << "find neighboring " << endl;
//   start = clock();
//   findNeighors( rad_sol);
//   finish = clock();
//   tr.Info << "\t findNeighors() running time: " << (float)(finish - start)/ CLOCKS_PER_SEC << " seconds" << endl;

//   //cout << "calculate surface " << endl;
//   start = clock();
//   calcSurface( rad_sol );
//   finish = clock();
//   tr.Info << "\t calcSurface() running time: " << (float)(finish - start)/ CLOCKS_PER_SEC << " seconds" << endl;

//   //for(int i = r1; i <= rN; i++){
//   //cout << "resID: " << i << "\t" << SurfaceList[i] << endl;
//   //}

// }


void PDB::calcSurface( float rad_sol )
{
	const Real SPARTA_PI = numeric::NumericTraits<Real>::pi();
	Mol conf = Conformers[1];
	Mol::iterator it, end;

	boost::unordered_map<string, float> STD_AREA;
	/*
	STD_AREA["ALA"] = 124; STD_AREA["ARG"] = 244; STD_AREA["ASN"] = 161; STD_AREA["ASP"] = 154; STD_AREA["cys"] = 94;
	STD_AREA["CYS"] = 94;  STD_AREA["GLU"] = 187; STD_AREA["GLN"] = 190; STD_AREA["GLY"] =  89;
	STD_AREA["HIS"] = 201; STD_AREA["ILE"] = 194; STD_AREA["LEU"] = 198; STD_AREA["LYS"] = 214;
	STD_AREA["MET"] = 215; STD_AREA["PHE"] = 221; STD_AREA["PRO"] = 150; STD_AREA["SER"] = 126;
	STD_AREA["THR"] = 152; STD_AREA["TRP"] = 265; STD_AREA["TYR"] = 236; STD_AREA["VAL"] = 169;

	STD_AREA["ALA"] = 28; STD_AREA["ARG"] = 27; STD_AREA["ASN"] = 28; STD_AREA["ASP"] = 25; STD_AREA["cys"] = 25;
	STD_AREA["CYS"] = 25; STD_AREA["GLU"] = 23; STD_AREA["GLN"] = 26; STD_AREA["GLY"] = 30;
	STD_AREA["HIS"] = 26; STD_AREA["ILE"] = 23; STD_AREA["LEU"] = 23; STD_AREA["LYS"] = 26;
	STD_AREA["MET"] = 28; STD_AREA["PHE"] = 26; STD_AREA["PRO"] = 22; STD_AREA["SER"] = 26;
	STD_AREA["THR"] = 25; STD_AREA["TRP"] = 28; STD_AREA["TYR"] = 25; STD_AREA["VAL"] = 24;
	*/

	for ( it = conf.begin(), end = conf.end(); it != end; ++it ) {
		PDB_Entry b = it->second;

		//if( b.atomName[0] == 'H' && b.atomName!="HN" && b.atomName!="H") continue;
		if ( b.atomName != "HN" && b.atomName!="N" && b.atomName!="CA"&& b.atomName!="CB"&& b.atomName!="C"&& b.atomName!="O" ) continue;

		//if( b.atomName[0] == 'H' ) continue; // for all heavy atoms
		//if( b.atomName != "O" ) continue; // for Carbonyl O only

		int resID = b.resNum;
		utility::vector0<int> Neighnors = NeighborList[ b.atomNum ];

		int fullNo = 0;
		int partNo = 0;
		bool fullInside, partInside;
		Vec3 cent, nCent, x, dx;


		//cout << "\tatom: " << resID << "\t" << b.atomName << "\ttotal points:" << SpherePointNo << endl;

		string b_atomName = b.atomName;
		if ( b.atomName[0] == 'C' && b.atomName != "C" ) b_atomName = "C"; // for all non-C' C
		if ( b.atomName == "C" ) b_atomName = "CO"; // for all C'=O C

		float b_rad = VDW_RAD[b_atomName] + rad_sol;
		Vec3Copy(cent, b.Coord);
		for ( int pointI = 0; pointI < SpherePointNo; pointI++ ) {
			Vec3Copy(x, cent);
			Vec3ScaleAdd(x, b_rad, SpherePoints[pointI]);

			fullInside = false;
			partInside = false;

			//cout << "\t\tpoint: " << pointI << "\tNeighnors size " << Neighnors.size() << endl;
			for ( int atomI = 0; atomI < (int)Neighnors.size(); atomI++ ) {
				PDB_Entry nAtomP = conf[ Neighnors[atomI] ];

				Vec3Copy(nCent, nAtomP.Coord);
				Vec3Copy(dx, x);
				Vec3Sub(dx, nCent);

				//cout << "\t\t\t" << resID << "\tNeighnor: " << nAtomP.resNum << "\t" << nAtomP.atomName << "\t" << Neighnors.size() << endl;
				if ( nAtomP.resNum < r1 || nAtomP.resNum > rN ) continue;


				string nAtomP_atomName = nAtomP.atomName;
				if ( nAtomP.atomName[0] == 'C' ) nAtomP_atomName = "C";
				if ( nAtomP.atomName == "C" ) nAtomP_atomName = "CO"; // for all C'=O C

				float nRad = VDW_RAD[nAtomP_atomName] + rad_sol;
				if ( dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2] > nRad * nRad ) {
					continue;
				}

				fullInside = true; //break;


				// for full area, ignored
				if ( resID == nAtomP.resNum ) { //full area
					partInside = true;
					break;
				}
			}
			//cout << endl;
			if ( ! fullInside ) fullNo++;

			if ( ! partInside ) partNo++; // full area, ignored
		}
		//cout << "\tatom: " << resID << "\t" << b.atomName << "\tfull points:" << fullNo*b_rad*b_rad*4.0f * (float) SPARTA_PI / SpherePointNo << "\tpart points:" << partNo*b_rad*b_rad *4.0f * (float) SPARTA_PI / SpherePointNo<<endl;


		AtomSurfaceFullList[resID][b.atomName] = fullNo * b_rad * b_rad  *4.0f * (float) SPARTA_PI / SpherePointNo; // solvent exposed area
		AtomSurfacePartList[resID][b.atomName] = partNo * b_rad * b_rad  *4.0f * (float) SPARTA_PI / SpherePointNo; // full area

		ResSurfaceFullList[resID] += (fullNo * b_rad * b_rad  *4.0f * (float) SPARTA_PI / SpherePointNo); // solvent exposed area * 4.0f * (float) M_SPARTA_PI / SpherePointNo
		ResSurfacePartList[resID] += (partNo * b_rad * b_rad  *4.0f * (float) SPARTA_PI / SpherePointNo);
	}


	boost::unordered_map< int,float>::iterator itS, endS;
	for ( itS = ResSurfaceFullList.begin(), endS = ResSurfaceFullList.end(); itS != endS; ++itS ) {
		string resName = residList[itS->first];
		// SurfaceFullList[ itS->first ] *=  ( 4.0f * (float) SPARTA_PI / (SpherePointNo * STD_AREA[resName]) );
		// SurfacePartList[ itS->first ] *=  ( 4.0f * (float) SPARTA_PI / (SpherePointNo * STD_AREA[resName]) );
	}

}


//find the neighoring atoms, distance < rad_a + rad_b + rad_sol, for atoms with given type, for all non-H atoms
void PDB::findNeighors(float rad_sol)
{
	Mol conf = Conformers[1];
	Mol::iterator it, itA, endA, end;

	//loop over all atoms
	for ( itA = conf.begin(), endA = conf.end(); itA != endA; ++itA ) {
		PDB_Entry a = itA->second;

		//if( a.atomName[0] == 'H' && a.atomName!="HN" && a.atomName!="H") continue; //skip non-HN H atoms
		if ( a.atomName != "HN" && a.atomName!="N" && a.atomName!="CA"&& a.atomName!="CB"&& a.atomName!="C"&& a.atomName!="O" ) continue; //skip non-HN H atoms

		//cout << a.atomName << endl;

		string a_atomName = a.atomName;
		if ( a.atomName[0] == 'C' ) a_atomName = "C";
		if ( a.atomName == "C" ) a_atomName = "CO"; // for all C'=O C
		float rad_a = VDW_RAD[ a_atomName ];

		//search for all other atoms
		for ( it = conf.begin(), end = conf.end(); it != end; ++it ) {
			PDB_Entry b = it->second;

			if ( b.atomName[0] == 'H'  && b.atomName!="HN" && b.atomName!="H" ) continue;

			string b_atomName = b.atomName;
			if ( b.atomName[0] == 'C' ) b_atomName = "C";
			if ( b.atomName == "C" ) b_atomName = "CO"; // for all C'=O C

			float rad_b = VDW_RAD[ b_atomName ];

			if ( getDist(a,b) < rad_a + rad_b + rad_sol ) {
				NeighborList[ a.atomNum ].push_back( b.atomNum );
				//cout << a.atomNum << "\t" << b.atomNum << "\t" << b.atomName << endl;
			}
		}

	}
}


// Same model from MOLMOL
void PDB::calcTriangles(
	double x0, double y0, double z0,
	double x1, double y1, double z1,
	double x2, double y2, double z2,
	int rowStartA[], int rowNo, int quad,
	int row0, int ind0, int ind1,
	int row2, int ind2,
	Vec3 *pointA)
{
	if ( row0 + 1 == row2 || row2 + 1 == row0 ) {
		int row0Size, row2Size;

		if ( row0 == - rowNo || row0 == rowNo ) {
			row0Size = 1;
		} else if ( row0 < 0 ) {
			row0Size = rowNo + row0;
		} else {
			row0Size = rowNo - row0;
		}

		if ( row2 == - rowNo || row2 == rowNo ) {
			row2Size = 1;
		} else if ( row2 < 0 ) {
			row2Size = rowNo + row2;
		} else {
			row2Size = rowNo - row2;
		}

		if ( ind0 < (quad + 1) * row0Size ) {
			ind0 += rowStartA[rowNo + row0];
			pointA[ind0][0] = (float) x0;
			pointA[ind0][1] = (float) y0;
			pointA[ind0][2] = (float) z0;
		}

		if ( ind1 < (quad + 1) * row0Size ) {
			ind1 += rowStartA[rowNo + row0];
			pointA[ind1][0] = (float) x1;
			pointA[ind1][1] = (float) y1;
			pointA[ind1][2] = (float) z1;
		}

		if ( ind2 < (quad + 1) * row2Size ) {
			ind2 += rowStartA[rowNo + row2];
			pointA[ind2][0] = (float) x2;
			pointA[ind2][1] = (float) y2;
			pointA[ind2][2] = (float) z2;
		}
	} else {
		double x01, y01, z01;
		double x12, y12, z12;
		double x20, y20, z20;
		double a;
		int rowMid, indMid01, indMid02, indMid12;

		x01 = x0 + x1;
		y01 = y0 + y1;
		z01 = z0 + z1;
		a = sqrt(x01 * x01 + y01 * y01 + z01 * z01);
		x01 /= a;
		y01 /= a;
		z01 /= a;

		x12 = x1 + x2;
		y12 = y1 + y2;
		z12 = z1 + z2;
		a = sqrt(x12 * x12 + y12 * y12 + z12 * z12);
		x12 /= a;
		y12 /= a;
		z12 /= a;

		x20 = x2 + x0;
		y20 = y2 + y0;
		z20 = z2 + z0;
		a = sqrt(x20 * x20 + y20 * y20 + z20 * z20);
		x20 /= a;
		y20 /= a;
		z20 /= a;

		rowMid = (row0 + row2) / 2;
		indMid01 = (ind0 + ind1) / 2;
		indMid02 = (ind0 + ind2) / 2;
		indMid12 = (ind1 + ind2) / 2;

		calcTriangles(
			x0, y0, z0,
			x01, y01, z01,
			x20, y20, z20,
			rowStartA, rowNo, quad,
			row0, ind0, indMid01,
			rowMid, indMid02,
			pointA);
		calcTriangles(
			x01, y01, z01,
			x1, y1, z1,
			x12, y12, z12,
			rowStartA, rowNo, quad,
			row0, indMid01, ind1,
			rowMid, indMid12,
			pointA);
		calcTriangles(
			x20, y20, z20,
			x12, y12, z12,
			x01, y01, z01,
			rowStartA, rowNo, quad,
			rowMid, indMid02, indMid12,
			row0, indMid01,
			pointA);
		calcTriangles(
			x20, y20, z20,
			x12, y12, z12,
			x2, y2, z2,
			rowStartA, rowNo, quad,
			rowMid, indMid02, indMid12,
			row2, ind2,
			pointA);
	}
}


// // Same model from MOLMOL
// void PDB::SphereCalcPoints(int divNo, Vec3 **pointAP, int *pointNoP)
// {
//   Vec3 *pointA;
//   int *rowStartA;
//   int pointNo, rowNo, rowSize, i;
//  runtime_assert( false );
//   rowNo = 1 << divNo;

//  //  rowStartA = (int *) malloc((2 * rowNo + 1) * sizeof(*rowStartA));


//   rowStartA[0] = 0;
//   rowStartA[1] = 1;
//   rowSize = 4;
//   for (i = 2; i <= rowNo; i++) {
//     rowStartA[i] = rowStartA[i - 1] + rowSize;
//     rowSize += 4;
//   }
//   for (i = rowNo + 1; i < 2 * rowNo + 1; i++) {
//     rowStartA[i] = rowStartA[i - 1] + rowSize;
//     rowSize -= 4;
//   }

//   pointNo = 4 * rowNo * rowNo + 2;
//  //  pointA = (Vec3 *) malloc(pointNo * sizeof(*pointA));

//   calcTriangles(
//   1.0, 0.0, 0.0,
//   0.0, 1.0, 0.0,
//   0.0, 0.0, 1.0,
//   rowStartA, rowNo, 0,
//   0, 0, rowNo,
//   rowNo, 0,
//   pointA);
//   calcTriangles(
//   0.0, 1.0, 0.0,
//   -1.0, 0.0, 0.0,
//   0.0, 0.0, 1.0,
//   rowStartA, rowNo, 1,
//   0, rowNo, 2 * rowNo,
//   rowNo, 0,
//   pointA);
//   calcTriangles(
//   -1.0, 0.0, 0.0,
//   0.0, -1.0, 0.0,
//   0.0, 0.0, 1.0,
//   rowStartA, rowNo, 2,
//   0, 2 * rowNo, 3 * rowNo,
//   rowNo, 0,
//   pointA);
//   calcTriangles(
//   0.0, -1.0, 0.0,
//   1.0, 0.0, 0.0,
//   0.0, 0.0, 1.0,
//   rowStartA, rowNo, 3,
//   0, 3 * rowNo, 4 * rowNo,
//   rowNo, 0,
//   pointA);
//   calcTriangles(
//   1.0, 0.0, 0.0,
//   0.0, 1.0, 0.0,
//   0.0, 0.0, -1.0,
//   rowStartA, rowNo, 0,
//   0, 0, rowNo,
//   - rowNo, 0,
//   pointA);
//   calcTriangles(
//   0.0, 1.0, 0.0,
//   -1.0, 0.0, 0.0,
//   0.0, 0.0, -1.0,
//   rowStartA, rowNo, 1,
//   0, rowNo, 2 * rowNo,
//   - rowNo, 0,
//   pointA);
//   calcTriangles(
//   -1.0, 0.0, 0.0,
//   0.0, -1.0, 0.0,
//   0.0, 0.0, -1.0,
//   rowStartA, rowNo, 2,
//   0, 2 * rowNo, 3 * rowNo,
//   - rowNo, 0,
//   pointA);
//   calcTriangles(
//   0.0, -1.0, 0.0,
//   1.0, 0.0, 0.0,
//   0.0, 0.0, -1.0,
//   rowStartA, rowNo, 3,
//   0, 3 * rowNo, 4 * rowNo,
//   - rowNo, 0,
//   pointA);

//  // free(rowStartA);

//   *pointAP = pointA;
//   *pointNoP = pointNo;
// }


// Structure-based H-N order parameter
// Zhang&Bruschweiler, 2002, JACS, 12654-12655
void PDB::calc_HN_S2( )
{
	Mol conf = Conformers[1];
	Mol::iterator it, end;

	//loop over all atoms
	std::map<int, string>::iterator itA, endA;

	float S2 = 0;
	for ( itA = residList.begin(), endA = residList.end(); itA != endA; ++itA ) {
		int resID = itA->first; //.resNum;
		string resName = itA->second;

		PDB_Entry O_prev = ATOMS[1][resID-1]["O"];
		PDB_Entry H      = ATOMS[1][resID]["H"];
		PDB_Entry HN     = ATOMS[1][resID]["HN"];

		//cout << resName << "\t" << resID << "\t" << O_prev.atomName << endl;

		if ( (O_prev.atomName).empty()|| ((H.atomName).empty() && (HN.atomName).empty() && resName != "PRO") ) continue; //skip if no O/HN atoms
		if ( (H.atomName).empty() ) H = HN;
		if ( resName == "PRO" ) H = ATOMS[1][resID]["N"]; //for Pro

		//search for all other heavy atoms
		S2 = 0;
		for ( it = conf.begin(), end = conf.end(); it != end; ++it ) {
			PDB_Entry b = it->second;

			int resID2 = b.resNum;
			if ( resID2 == resID || resID2 == resID-1 || b.atomName[0] == 'H' || (b.atomName[1] == 'H'&&b.atomName[1] != 'N') ) continue; //skip all H

			float D_OK = getDist(O_prev,b); // distance between O and heavy atom
			float D_HK = getDist(H,b); // distance between O and heavy atom

			S2 += exp(-1.0*D_OK);
			S2 += 0.8*exp(-1.0*D_HK);
			//cout << resName << "\t" << resID << "\t" << O_prev.atomName << "\t" << b.atomName << "\t" << resID2 ;
			//cout << "\t" << D_OK << "\t" << D_HK << "\t" <<  S2 << endl;
		}

		S2 *= 2.656;

		S2 = (exp(S2)-exp(-S2))/(exp(S2)+exp(-S2))-0.1;
		//cout << resName << "\t" << S2 << endl;
		HN_S2[resID] = S2;
	}

	//smooth for terminus
	if ( HN_S2[r1+2]>0 && HN_S2[r1+1]>0 ) HN_S2[r1] = HN_S2[r1+1]-fabs(HN_S2[r1+1]-HN_S2[r1+2]);
	if ( HN_S2[rN-2]>0 && HN_S2[rN-1]>0 ) HN_S2[rN] = HN_S2[rN-1]-fabs(HN_S2[rN-1]-HN_S2[rN-2]);

}


// calculate electric field effects
void PDB::calc_ElectricField()
{
	const Real SPARTA_PI = numeric::NumericTraits<Real>::pi();
	Mol conf = Conformers[1];
	Mol::iterator it, end;

	boost::unordered_map<string, string> targetAtomList;
	boost::unordered_map<string, string>::iterator itT, endT;
	targetAtomList["HN"]="N";
	targetAtomList["HA"]="CA";
	//targetAtomList["CA"]="N";

	boost::unordered_map<string, float> Qlist;
	Qlist["C"] = -0.9612;
	Qlist["O"] =  1.39374;
	Qlist["N"] =  0.7209;

	//loop over all atoms
	std::map<int, string>::iterator itA, endA;
	ElectricField.clear();
	for ( itA = residList.begin(), endA = residList.end(); itA != endA; ++itA ) {
		int resID = itA->first; //.resNum;
		string resName = itA->second;

		// tr.Info << resID << "\t" << resName << endl;
		for ( itT = targetAtomList.begin(), endT = targetAtomList.end(); itT != endT; ++itT ) {
			PDB_Entry target  = ATOMS[1][resID][itT->first];
			PDB_Entry partner = ATOMS[1][resID][itT->second];

			if ( (target.atomName).empty() ) continue;

			// tr.Info << resID << "\t" << resName << "\t" << target.atomName << "\t" << target.resNum << endl;
			for ( it = conf.begin(), end = conf.end(); it != end; ++it ) {
				PDB_Entry b = it->second;

				int resID2 = b.resNum;
				if ( fabs((double) (resID2 - resID)) <= 1 ) continue; // FOR WINDOWS BUILD //skip all atoms from present and adjacent residues
				string atomName2 = b.atomName;
				if ( atomName2 == "O" && target.atomName == "HN" ) continue;

				//cout << resID << "\t" << resName << "\t" << target.atomName << "\t" << target.resNum << "\t" << atomName2 << endl;
				if ( atomName2=="O" || atomName2.substr(0,2)=="OD" || atomName2.substr(0,2)=="OE" || atomName2=="C" || atomName2=="N" ) {
					float c = cos(getBondAngle(partner,target,b)*SPARTA_PI/180.0);
					float dist = getDist(target,b);
					if ( dist > 3.0 ) continue;
					ElectricField[resID][target.atomName] += Qlist[atomName2.substr(0,1)]*c/(dist*dist);
					//cout << resName << "\t" << resID << "\t" << target.atomName << "\t" << b.atomName << "\t" << resID2;
					//cout << "\t" << c << "\t" << dist << "\t" <<  ElectricField[resID][target.atomName] << endl;
				}
			}
		}


	}
}


void PDB::collect_HN_S2_and_EF( )
{
	Mol conf = Conformers[1];
	Mol::iterator it, end;

	ElectricField.clear();

	const Real SPARTA_PI = numeric::NumericTraits<Real>::pi();
	const Real RADS_PER_DEG = SPARTA_PI / 180.;

	boost::unordered_map<string, string> targetAtomList;
	boost::unordered_map<string, string>::iterator itT;
	targetAtomList["HN"]="N";
	targetAtomList["HA"]="CA";
	//targetAtomList["CA"]="N";

	boost::unordered_map<string, float> Qlist;
	Qlist["C"] = -0.9612;
	Qlist["O"] =  1.39374;
	Qlist["N"] =  0.7209;

	//loop over all atoms
	std::map<int, string>::iterator itA, endA;

	float S2 = 0;
	for ( itA = residList.begin(), endA = residList.end(); itA != endA; ++itA ) {
		int resID = itA->first; //.resNum;
		string resName = itA->second;

		PDB_Entry O_prev = ATOMS[1][resID-1]["O"];
		PDB_Entry H      = ATOMS[1][resID]["H"];
		PDB_Entry HN     = ATOMS[1][resID]["HN"];

		//cout << resName << "\t" << resID << "\t" << O_prev.atomName << endl;

		if ( (O_prev.atomName).empty()|| ((H.atomName).empty() && (HN.atomName).empty() && resName != "PRO") ) continue; //skip if no O/HN atoms
		if ( (H.atomName).empty() ) H = HN;
		if ( resName == "PRO" ) H = ATOMS[1][resID]["N"]; //for Pro


		PDB_Entry EF_target_HA  = ATOMS[1][resID]["HA"];
		PDB_Entry EF_partner_HA = ATOMS[1][resID]["CA"];
		PDB_Entry EF_target_HN  = ATOMS[1][resID]["HN"];
		PDB_Entry EF_partner_HN = ATOMS[1][resID]["N"];

		//search for all other heavy atoms
		S2 = 0;
		for ( it = conf.begin(), end = conf.end(); it != end; ++it ) {
			PDB_Entry b = it->second;

			int resID2 = b.resNum;
			// if( resID2 == resID || resID2 == resID-1 || b.atomName[0] == 'H' || (b.atomName[1] == 'H'&&b.atomName[1] != 'N')) continue; //skip all H
			if ( !( resID2 == resID || resID2 == resID-1 || b.atomName[0] == 'H' || (b.atomName[1] == 'H'&&b.atomName[1] != 'N')) ) {
				float D_OK = getDist(O_prev,b); // distance between O and heavy atom
				float D_HK = getDist(H,b); // distance between O and heavy atom

				S2 += exp(-1.0*D_OK);
				S2 += 0.8*exp(-1.0*D_HK);
				//cout << resName << "\t" << resID << "\t" << O_prev.atomName << "\t" << b.atomName << "\t" << resID2 ;
				//cout << "\t" << D_OK << "\t" << D_HK << "\t" <<  S2 << endl;
			} // S2 calucalation

			if ( fabs( (double) (resID2 - resID)) <= 1 ) continue; // FOR WINDOWS BUILD //skip all atoms from present and adjacent residues
			string atomName2 = b.atomName;

			if ( !(EF_target_HA.atomName).empty() ) {
				//cout << resID << "\t" << resName << "\t" << target.atomName << "\t" << target.resNum << "\t" << atomName2 << endl;
				if ( atomName2=="O" || atomName2.substr(0,2)=="OD" || atomName2.substr(0,2)=="OE" || atomName2=="C" || atomName2=="N" ) {
					//float c = cos(getBondAngle(EF_partner_HA,EF_target_HA,b)*RADS_PER_DEG);
					float dist = getDist(EF_target_HA,b);
					if ( dist <= 3.0 ) { //continue;
						ElectricField[resID]["HA"] += Qlist[atomName2.substr(0,1)]*cos(getBondAngle(EF_partner_HA,EF_target_HA,b)*RADS_PER_DEG)/(dist*dist);
					}
					//cout << resName << "\t" << resID << "\t" << target.atomName << "\t" << b.atomName << "\t" << resID2;
					//cout << "\t" << c << "\t" << dist << "\t" <<  ElectricField[resID][target.atomName] << endl;
				}
			}
			if ( !(EF_target_HN.atomName).empty() ) {
				//cout << resID << "\t" << resName << "\t" << EF_target_HN.atomName << "\t" << EF_target_HN.resNum << "\t" << atomName2 << endl;
				if ( atomName2.substr(0,2)=="OD" || atomName2.substr(0,2)=="OE" || atomName2=="C" || atomName2=="N" ) {
					//float c = cos(getBondAngle(EF_partner_HN,EF_target_HN,b)*RADS_PER_DEG);
					float dist = getDist(EF_target_HN,b);
					if ( dist <= 3.0 ) { //continue;
						ElectricField[resID]["HN"] += Qlist[atomName2.substr(0,1)]*cos(getBondAngle(EF_partner_HN,EF_target_HN,b)*RADS_PER_DEG)/(dist*dist);
					}
					//cout << resName << "\t" << resID << "\t" << EF_target_HN.atomName << "\t" << b.atomName << "\t" << resID2;
					//cout << "\t" << c << "\t" << dist << "\t" <<  ElectricField[resID][EF_target_HN.atomName] << endl;
				}
			}

		}

		S2 *= 2.656;

		S2 = (exp(S2)-exp(-S2))/(exp(S2)+exp(-S2))-0.1;
		//cout << resName << "\t" << S2 << endl;
		HN_S2[resID] = S2;
	}

	//smooth for terminus
	if ( HN_S2[r1+2]>0 && HN_S2[r1+1]>0 ) HN_S2[r1] = HN_S2[r1+1]-fabs(HN_S2[r1+1]-HN_S2[r1+2]);
	if ( HN_S2[rN-2]>0 && HN_S2[rN-1]>0 ) HN_S2[rN] = HN_S2[rN-1]-fabs(HN_S2[rN-1]-HN_S2[rN-2]);

}

}
}
