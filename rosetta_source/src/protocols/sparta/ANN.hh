// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
/*                                                                            */
/*                     ---- TALOS C++ verison ----                            */
/*           TALOS: Torsion Angle Likeness Optimized By Shifts.               */
/*        Yang Shen, Gabriel Cornilescu, Frank Delaglio, and Ad Bax           */
/*                   NIH Laboratory of Chemical Physics                       */
/*                        veriosn, 2009.05.31.1551                            */
/*                                                                            */
/*                        for any problem, please contact                     */
/*                           shenyang@niddk.nih.gov                           */
/*                                                                            */
/******************************************************************************/

/* ANN.cpp: class for a simple Artificial Neural Network
//          a simple two-level ANN model with one hidden layer, each level with three different runs
*/


#ifndef ANN_H
#define ANN_H

#include <boost/unordered_map.hpp>
#include <protocols/sparta/GDB.hh>
#include <utility/vector0.hh>

namespace protocols {
namespace sparta {

class ANN
{
  char buf[30];
  std::string slash_char;
  int input_code; // code for input format
  // 0 - all CS and AA data, default
  // 1 - no CS for residue i-1
  // 2 - no CS for residue i
  // 3 - no CS for residue i+1

  std::string DB_PATH, DB_NAME_PREFIX;
  std::string DB_FNAME_LEVEL1_1, DB_FNAME_LEVEL1_2, DB_FNAME_LEVEL1_3, DB_FNAME_LEVEL2_1, DB_FNAME_LEVEL2_2, DB_FNAME_LEVEL2_3;
  int N1_NODE_I, N1_NODE_H, N1_NODE_O;	// Number of node in input, hidden and output layers of the 1st level ANN
  int N2_NODE_I, N2_NODE_H, N2_NODE_O;	// Number of node in input, hidden and output layers of the 1st level ANN

  // 1st level ANN weights and bias
  boost::unordered_map<int, utility::vector0<float> > WI_1, WI_2, WI_3;	// Weighting utility::vector0s for converting input layer, N_I*N_I matrix
  utility::vector0<float>  BI_1, BI_2, BI_3;			// Bias utility::vector0s (converting input layer), size of N_I
  boost::unordered_map<int, utility::vector0<float> > WL1_1, WL1_2, WL1_3;	// Weighting utility::vector0s connecting input and hidden layer, N_H*N_I matrix
  utility::vector0<float>  BL1_1, BL1_2, BL1_3;				// Bias utility::vector0s (connecting input and hidden layer), size of N_H
  boost::unordered_map<int, utility::vector0<float> > WL2_1, WL2_2, WL2_3;	// Weighting utility::vector0s connecting hidden and output layer, N_O*N_H matrix
  utility::vector0<float>  BL2_1, BL2_2, BL2_3;				// Bias utility::vector0s (connecting hidden and output layer), size of N_O

  // 1st level ANN weights and bias
  boost::unordered_map<int, utility::vector0<float> > W2I_1, W2I_2, W2I_3;	// Weighting utility::vector0s for converting input layer, N_I*N_I matrix
  utility::vector0<float>  B2I_1, B2I_2, B2I_3;			// Bias utility::vector0s (converting input layer), size of N_I
  boost::unordered_map<int, utility::vector0<float> > W2L1_1, W2L1_2, W2L1_3;	// Weighting utility::vector0s connecting input and hidden layer, N_H*N_I matrix
  utility::vector0<float>  B2L1_1, B2L1_2, B2L1_3;				// Bias utility::vector0s (connecting input and hidden layer), size of N_H
  boost::unordered_map<int, utility::vector0<float> > W2L2_1, W2L2_2, W2L2_3;	// Weighting utility::vector0s connecting hidden and output layer, N_O*N_H matrix
  utility::vector0<float>  B2L2_1, B2L2_2, B2L2_3;				// Bias utility::vector0s (connecting hidden and output layer), size of N_O
public:

  boost::unordered_map<int, utility::vector0<float> > ANN_IN_MTX;		// input matrix for neural netwrok calculation, indexed by resID and with an input utility::vector0 of size 32
  // singlet only!!!(6*2 for chemical shifts and 20 for amino acid sequence)
  boost::unordered_map<int, utility::vector0<float> > ANN_IN_MTX_LEVEL1;		// tripet model input matrix for 1st level ANN
  boost::unordered_map<int, utility::vector0<float> > ANN_IN_MTX_LEVEL2;		// tripet modelinput matrix for 2nd level ANN

  boost::unordered_map<int, utility::vector0<float> > ANN_OUT_MTX_LEVEL1;		// output matrix, , indexed by resID and with an  utility::vector0 of size 3

  boost::unordered_map< int, boost::unordered_map< std::string, std::string> >::iterator it;

  boost::unordered_map<int, utility::vector0<float> > ANN_OUT_MTX_LEVEL2;		// output matrix, , indexed by resID and with an  utility::vector0 of size 3

  ANN();
  ANN( const std::string& dPATH, const std::string& dNAME_PREFIX );
  ANN( int N1_nodeI, int N1_nodeH, int N1_nodeO, const std::string& dPATH, const std::string& dNAME_PREFIX );
  ANN( int N1_nodeI, int N1_nodeH, int N1_nodeO, int N2_nodeI, int N2_nodeH, int N2_nodeO, const std::string& dPATH, const std::string& dNAME_PREFIX );

  void init( int N1_nodeI, int N1_nodeH, int N1_nodeO, int N2_nodeI, int N2_nodeH, int N2_nodeO, const std::string& dPATH, const std::string& dNAME_PREFIX );
  void getSlashChar();

  void set_input_code(int c);

  void loadWeights(); // load all weighting and bias
  void loadWeightBias3( const std::string& fName, boost::unordered_map<int, utility::vector0<float> > &W1, utility::vector0<float> &B1,
    boost::unordered_map<int, utility::vector0<float> > &W2, utility::vector0<float> &B2, boost::unordered_map<int, utility::vector0<float> > &W3, utility::vector0<float> &B3,
    int N_W_row, int N_W_col, int N_B);
  // load weighting (N_W_row*N_W_col) and bias (N_B) from a given file containing all three sets data


  void calcLevel1();
  void calcLevel2();
  void runSpartaANN( boost::unordered_map<int, utility::vector0<float> > &inMatrix ); // used by SPARTA

  void applyANNTransformation( utility::vector0<float> &inp, boost::unordered_map<int, utility::vector0<float> > &w, utility::vector0<float> &b, utility::vector0<float> &out, int code);

  void applyVecAverage( utility::vector0<float> &v1, utility::vector0<float> &v2, utility::vector0<float> &v3, utility::vector0<float> &vout);
  //calculate 'confidence-averaged' utility::vector0 of three utility::vector0s v1, v2, v3
  void applyVecNormalization( utility::vector0<float> &v);
  //apply normalization
  float getConfidence( utility::vector0<float> &v);

  //check the number of atom without CS for a given residue
  int getNumberMissCS( utility::vector0<float> &v);

  char * ftoa( float n, char *buff, char f='g', int prec=6 );
  char * itoa( int n, char *buff, int base=10 );

};

}
}
#endif
