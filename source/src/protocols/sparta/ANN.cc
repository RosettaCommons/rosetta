// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
/*                                                                            */
/*                     ---- TALOS C++ verison ----                            */
/*           TALOS: Torsion Angle Likeness Optimized By Shifts.               */
/*        Yang Shen, Gabriel Cornilescu, Frank Delaglio, and Ad Bax           */
/*                   NIH Laboratory of Chemical Physics                       */
/*                     version, 1.00 (build 2010.0607.00)                     */
/*                                                                            */
/*                        for any problem, please contact                     */
/*                           shenyang@niddk.nih.gov                           */
/*                                                                            */
/******************************************************************************/


/* ANN.cpp: class for a simple Artificial Neural Network */

#include <boost/unordered_map.hpp>
#include <protocols/sparta/ANN.hh>
// Package Headers
#include <core/types.hh>
// Project Headers

// Utility headers
#include <basic/Tracer.hh>

#include <cstdio>

#include <utility/vector0.hh>

#include <protocols/sparta/GDB.hh>
#include <utility/vector1.hh>


static THREAD_LOCAL basic::Tracer tr( "protocols.sparta" );

namespace protocols {
namespace sparta {

using namespace std;
using namespace core;

#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))

ANN::ANN()
{
	N1_NODE_I = 96; N1_NODE_H = 20; N1_NODE_O = 3;
	N2_NODE_I = 9;  N2_NODE_H = 6;  N2_NODE_O = 3;
	input_code = 0;

	getSlashChar();
}


// constructor with initial setup of weighting factors' database
ANN::ANN(const string& dPATH, const string& dNAME_PREFIX)
{
	ANN();
	getSlashChar();

	N1_NODE_I = 96; N1_NODE_H = 20; N1_NODE_O = 3;
	N2_NODE_I = 9;  N2_NODE_H = 6;  N2_NODE_O = 3;
	input_code = 0;

	DB_PATH = dPATH;
	DB_NAME_PREFIX = dNAME_PREFIX;
	loadWeights();
}


// constructor with initial setup of weighting factors' database and number of node in each layer
ANN::ANN(int N1_nodeI, int N1_nodeH, int N1_nodeO, const string& dPATH, const string& dNAME_PREFIX)
{
	ANN();
	getSlashChar();

	N1_NODE_I = N1_nodeI;
	N1_NODE_H = N1_nodeH;
	N1_NODE_O = N1_nodeO;
	input_code = 0;

	N2_NODE_I = 9;  N2_NODE_H = 6;  N2_NODE_O = 3;

	DB_PATH = dPATH;
	DB_NAME_PREFIX = dNAME_PREFIX;
	loadWeights();
}


// constructor with initial setup of weighting factors' database and number of node in each layer
ANN::ANN(int N1_nodeI, int N1_nodeH, int N1_nodeO, int N2_nodeI, int N2_nodeH, int N2_nodeO, const string& dPATH, const string& dNAME_PREFIX)
{
	ANN();
	getSlashChar();

	N1_NODE_I = N1_nodeI;
	N1_NODE_H = N1_nodeH;
	N1_NODE_O = N1_nodeO;
	input_code = 0;

	N2_NODE_I = N2_nodeI;  N2_NODE_H = N2_nodeH;  N2_NODE_O = N2_nodeO;

	DB_PATH = dPATH;
	DB_NAME_PREFIX = dNAME_PREFIX;
	loadWeights();
}


// constructor with initial setup of weighting factors' database and number of node in each layer
void ANN::init(int N1_nodeI, int N1_nodeH, int N1_nodeO, int N2_nodeI, int N2_nodeH, int N2_nodeO, const string& dPATH, const string& dNAME_PREFIX)
{
	//ANN();
	getSlashChar();

	N1_NODE_I = N1_nodeI;
	N1_NODE_H = N1_nodeH;
	N1_NODE_O = N1_nodeO;
	input_code = 0;

	N2_NODE_I = N2_nodeI;  N2_NODE_H = N2_nodeH;  N2_NODE_O = N2_nodeO;

	DB_PATH = dPATH;
	DB_NAME_PREFIX = dNAME_PREFIX;
	loadWeights();
}


void ANN::set_input_code(int c)
{
	input_code = c;
}


void ANN::loadWeights() // load weighting and bias
{
	string wName;

	//load weights and bias for 1st level ANN
	wName = DB_PATH+slash_char+DB_NAME_PREFIX+".level1.WI.tab"; // file name of first weight/bias for input level
	loadWeightBias3(wName, WI_1, BI_1, WI_2, BI_2, WI_3, BI_3, N1_NODE_I, N1_NODE_I, N1_NODE_I);

	wName = DB_PATH+slash_char+DB_NAME_PREFIX+".level1.WL1.tab"; // file name of first weight/bias for connecting input and hidden level
	loadWeightBias3(wName, WL1_1, BL1_1, WL1_2, BL1_2, WL1_3, BL1_3, N1_NODE_H, N1_NODE_I, N1_NODE_H);

	wName = DB_PATH+slash_char+DB_NAME_PREFIX+".level1.WL2.tab"; // file name of first weight/bias for connecting input and hidden level
	loadWeightBias3(wName, WL2_1, BL2_1, WL2_2, BL2_2, WL2_3, BL2_3, N1_NODE_O, N1_NODE_H, N1_NODE_O);

	/* for a 2-level NN
	//load weights and bias for 2nd level ANN
	wName = DB_PATH+slash_char+DB_NAME_PREFIX+".level2.WI.tab"; // file name of first weight/Bias for input level
	loadWeightBias3(wName, W2I_1, B2I_1, W2I_2, B2I_2, W2I_3, B2I_3, N2_NODE_I, N2_NODE_I, N2_NODE_I);

	wName = DB_PATH+slash_char+DB_NAME_PREFIX+".level2.WL1.tab"; // file name of first weight/Bias for connecting input and hidden level
	loadWeightBias3(wName, W2L1_1, B2L1_1, W2L1_2, B2L1_2, W2L1_3, B2L1_3, N2_NODE_H, N2_NODE_I, N2_NODE_H);

	wName = DB_PATH+slash_char+DB_NAME_PREFIX+".level2.WL2.tab"; // file name of first weight/Bias for connecting input and hidden level
	loadWeightBias3(wName, W2L2_1, B2L2_1, W2L2_2, B2L2_2, W2L2_3, B2L2_3, N2_NODE_O, N2_NODE_H, N2_NODE_O);
	*/
}


// load weighting (N_W_row*N_W_col) and bias (N_B) from a given file contains all three sets data
void ANN::loadWeightBias3( string const& fName,
	ANN_Matrix &W1, ANN_Vector &B1,
	ANN_Matrix &W2, ANN_Vector &B2,
	ANN_Matrix &W3, ANN_Vector &B3,
	int N_W_row, int N_W_col, int /* N_B */)
{
	string str;
	//cout << "Reading ANN Weights and Bias Set 1 " << fName << endl;
	GDB W_Tab( fName );
	//cout << W1_Tab.Entries.size() << "\n";

	int index = 0;
	int row = N_W_row, col = N_W_col;
	for ( auto & Entrie : W_Tab.Entries ) {//iterate rows
		int check = index/row; //cout << check << endl;
		for ( int  i = 0; i < col; i++ ) {
			str = itoa(i+1,buf);
			float w = atof((Entrie.second[str]).c_str());
			if ( check == 0 ) W1[ index ].push_back( w ); // assign to weight matrix 1
			else if ( check == 1 ) W2[ index-row ].push_back( w ); // assign to weight matrix 2
			else if ( check == 2 ) W3[ index-row*2 ].push_back( w ); // assign to weight matrix 2
			else tr.Error << "Wrong size for matrix " << fName << " ... \n";
		}

		if ( check == 0 ) B1.push_back( atof((Entrie.second["b"]).c_str()) );
		else if ( check == 1 ) B2.push_back( atof((Entrie.second["b"]).c_str()) );
		else if ( check == 2 ) B3.push_back( atof((Entrie.second["b"]).c_str()) );
		index++;
	}
	// if( index > row ) tr.Error << "Wrong size for matrix " << fName << " ... \n";
}


// perform 1st level ANN calculation
// input ANN_IN_MTX_LEVEL1
void ANN::calcLevel1()
{
	//ANN_Matrix ANN_OUT_MTX;
	ANN_Matrix::iterator itV, end;
	for ( itV = ANN_IN_MTX_LEVEL1.begin(), end = ANN_IN_MTX_LEVEL1.end(); itV != end; ++itV ) { //for each tripet input
		//cout << itV->first << endl;

		//apply input layer transformation
		ANN_Vector IL1, IL2, IL3;
		applyANNTransformation(itV->second, WI_1, BI_1, IL1, 1);
		applyANNTransformation(itV->second, WI_2, BI_2, IL2, 1);
		applyANNTransformation(itV->second, WI_3, BI_3, IL3, 1);

		//apply 1st hidden layer transformation
		ANN_Vector HL1, HL2, HL3;
		applyANNTransformation(IL1, WL1_1, BL1_1, HL1, 1);
		applyANNTransformation(IL2, WL1_2, BL1_2, HL2, 1);
		applyANNTransformation(IL3, WL1_3, BL1_3, HL3, 1);

		//apply output layer transformation
		ANN_Vector OL1, OL2, OL3;
		applyANNTransformation(HL1, WL2_1, BL2_1, OL1, 0);
		applyANNTransformation(HL2, WL2_2, BL2_2, OL2, 0);
		applyANNTransformation(HL3, WL2_3, BL2_3, OL3, 0);

		ANN_Vector OUT1;
		applyVecAverage(OL1,OL2,OL3,OUT1);
		//cout << OUT1[0] << "\t" << OUT1[1] << "\t" << OUT1[2] << "\t" << endl;
		ANN_OUT_MTX_LEVEL1[itV->first] = OUT1;
	}
}


// perform 2nd level ANN calculation
// input ANN_IN_MTX_LEVEL2
void ANN::calcLevel2()
{
	//ANN_Matrix ANN_OUT_MTX;
	ANN_Matrix::iterator itV, end;
	for ( itV = ANN_IN_MTX_LEVEL2.begin(), end = ANN_IN_MTX_LEVEL2.end(); itV != end; ++itV ) { //for each tripet input
		//cout << itV->first << endl;

		//apply input layer transformation
		ANN_Vector IL1, IL2, IL3;
		applyANNTransformation(itV->second, W2I_1, B2I_1, IL1, 1);
		applyANNTransformation(itV->second, W2I_2, B2I_2, IL2, 1);
		applyANNTransformation(itV->second, W2I_3, B2I_3, IL3, 1);

		//apply 1st hidden layer transformation
		ANN_Vector HL1, HL2, HL3;
		applyANNTransformation(IL1, W2L1_1, B2L1_1, HL1, 1);
		applyANNTransformation(IL2, W2L1_2, B2L1_2, HL2, 1);
		applyANNTransformation(IL3, W2L1_3, B2L1_3, HL3, 1);

		//apply output layer transformation
		ANN_Vector OL1, OL2, OL3;
		applyANNTransformation(HL1, W2L2_1, B2L2_1, OL1, 0);
		applyANNTransformation(HL2, W2L2_2, B2L2_2, OL2, 0);
		applyANNTransformation(HL3, W2L2_3, B2L2_3, OL3, 0);

		ANN_Vector OUT2;
		applyVecAverage(OL1,OL2,OL3,OUT2);
		//cout << OUT2[0] << "\t" << OUT2[1] << "\t" << OUT2[2] << "\t" << endl;
		ANN_OUT_MTX_LEVEL2[itV->first] = OUT2;
	}
}


void ANN::runSpartaANN(ANN_Matrix &inMatrix)
{
	ANN_IN_MTX_LEVEL1 = inMatrix;
	//loadWeights();

	// run 1st level ANN prediction, get ANN_OUT_MTX_LEVEL1
	calcLevel1();

}


//apply an ANN transformation for input inp and with transformation weights w, bias b
void ANN::applyANNTransformation( ANN_Vector &inp, ANN_Matrix &w, ANN_Vector &b, ANN_Vector &out, int code)
{
	// needs to check for size mismatch among inp, w and b ??????
	if ( inp.size() != w[0].size() || w.size() != b.size() ) {
		tr.Error << " ANN prediction failed with inconsistent data!" << endl;
		return;
	}

	for ( Size  i = 0; i < w.size(); i++ ) {
		float sum = 0;
		for ( Size  j = 0; j < inp.size(); j++ ) sum += inp[j]*w[i][j];

		sum += b[i];
		if ( code == 1 ) out.push_back( 2.0/(1.0+exp(-2.0*sum))-1.0 ); // apply 'tansig' transfer function
		else if ( code == 0 ) out.push_back( sum ); // apply 'linear' transfer function
	}
}


//calculate 'confidence-averaged' utility::vector0 of three utility::vector0s v1, v2, v3
void ANN::applyVecAverage(ANN_Vector &v1, ANN_Vector &v2, ANN_Vector &v3, ANN_Vector &vout)
{
	if ( v1.size() == v2.size() && v1.size() == v3.size() ) {
		//float conf1 = getConfidence(v1);
		//float conf2 = getConfidence(v2);
		//float conf3 = getConfidence(v3);
		//cout << "XXX" << conf1 << "\t" << conf2 << "\t" << conf3<< endl;

		for ( Size i=0; i<v1.size(); i++ ) {
			//cout << "XXX" << i << endl;
			//out.push_back( (v1[i]*conf1+v2[i]*conf2+v3[i]*conf3)/(conf1+conf2+conf3) );
			vout.push_back( (v1[i]+v2[i]+v3[i])/3.0 );
			//vout.push_back(v3[i]);
		}

	}
}


//apply normalization
void ANN::applyVecNormalization(ANN_Vector &v)
{
	float a=v[0], b=v[1], c=v[2];
	if ( a>1 ) a=1.0; else if ( a<0 ) a=0.0;
	if ( b>1 ) b=1.0; else if ( b<0 ) b=0.0;
	if ( c>1 ) c=1.0; else if ( c<0 ) c=0.0;

	float sum=a+b+c;
	a/=sum; b/=sum; c/=sum;
	v.clear();
	v.push_back(a); v.push_back(b); v.push_back(c);
}


float ANN::getConfidence(ANN_Vector &v)
{
	//cout << v.size() << "\t" << v[0] << "\t" << v[1] << "\t" << v[2] << endl;

	if ( v.size() != 3 ) return -1.0;

	return 2.0*MAX(v[0], MAX(v[1],v[2])) - (v[0]+v[1]+v[2]) + MIN(v[0], MIN(v[1],v[2]));
}


//check the number of atom without CS for a given residue
int ANN::getNumberMissCS(ANN_Vector &v)
{
	int cnt = 0;
	cnt+=(v[1]==1); cnt+=(v[3]==1); cnt+=(v[5]==1); cnt+=(v[7]==1); cnt+=(v[9]==1); cnt+=(v[11]==1);
	return cnt;
}


// return a character string for an int type number
char * ANN::itoa( int n, char *buff, int /*base*/ )
{
	sprintf(buff, "%d", n);
	return buff;
}


// retrun a character string for a float type number
char * ANN::ftoa( float n, char *buff, char f, int prec )
{
	if ( !(f=='f' || f=='F' || f=='e' || f=='E' || f=='g' || f=='G') ) {
		f = 'f';
	}
	char format[20];
	char *fs = format;    // generate format string
	*fs++ = '%';    //   "%.<prec>l<f>"
	if ( prec >= 0 ) {
		if ( prec > 99 ) {   // buf big enough for precision?
			prec = 99;
		}
		*fs++ = '.';
		if ( prec >= 10 ) {
			*fs++ = prec / 10 + '0';
			*fs++ = prec % 10 + '0';
		} else {
			*fs++ = prec + '0';
		}
	}
	*fs++ = 'l';
	*fs++ = f;
	*fs = '\0';
	sprintf( buff, format, n );

	return buff;
}

void ANN::getSlashChar()
{
	if ( getenv( "PATH" ) != nullptr ) {
		string temp = getenv( "PATH" );
		if ( temp.find("/") != string::npos ) slash_char = "/"; // unix
		else if ( temp.find("\\") != string::npos ) slash_char = "\\"; // Windows
	} else slash_char = "/"; //default Windows
}

}
}
