/*                                                                            */
/*                           ----  SPARTA   ----                              */
/*     Shifts Prediction from Analogue of Residue type and Torsion Angle      */
/*                           Yang Shen and Ad Bax                             */
/*                    J. Biomol. NMR, 38, 289-302 (2007)                      */
/*                 NIH, NIDDK, Laboratory of Chemical Physics                 */
/*                     version, 1.2 (build 2009.0928.17)                     */
/*                                                                            */
/*                      for any problem, please contact                       */
/*                          shenyang@niddk.nih.gov                            */
/*                                                                            */
/******************************************************************************/

#ifndef SPARTA_H
#define SPARTA_H

#include <protocols/sparta/GDB.hh>
#include <protocols/sparta/PDB.hh>
#include <protocols/sparta/ANN.hh>

#include <boost/unordered_map.hpp>

#include <core/pose/Pose.fwd.hh>

#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility/vector0.hh>

namespace protocols {
namespace sparta {

typedef utility::vector0< std::string > StringList;
typedef boost::unordered_map<float, boost::unordered_map<float, float> > PHIPSI_ERR_SURF;

/* this is not thread-safe... */

class Sparta {
	///mini specific additions...
public:

	Sparta( std::string const & chem_shifts );
	~Sparta();

	static void register_options();
	core::Real score_pose( core::pose::Pose const & pose );
	void check_pose( core::pose::Pose const & pose );
	utility::vector1< float > score_pose_per_residue( core::pose::Pose const & pose );
public:
	class SpartaLib;
	SpartaLib& lib() { return *lib_instance_; } //no it is not a constant reference since we do the calculation in SPARTA_LIB. assumption here: no threads
	core::Real run_A_ANN_Prediction(); // Run ANN prediction for a single protein
	static SpartaLib* lib_instance_;
	GDB REF_CS_Tab;
	bool bCreateOutput_;
	std::string refCSFileName;

	static bool options_registered_;

public: ///{ most of the original SPARTA class goes into SpartaLib -- to be reused between different evaluators...
	class SpartaLib { //private sub-class so the elements are public such that they can be accessed from main class without too much text changes..
	public:
		SpartaLib();
		void setup_for_scoring(core::pose::Pose const & pose);
		void deallocate_arrays();
		void init();
		GDB get_ANN_data( bool create_output );
		void getResInfo( bool create_output ); //Get the list of useful shifts from a given residue. get 2nd chemical shift and apply correction
		float getANN_PredError(float phi, float psi, std::string aa, std::string aName); // get the prediction error from error surface
		void init_PredErrorSurface();

		// get random coil chemical shift for atom 'aName' of residue 'resName'
		float getRC(const std::string& resName, const std::string& aName);
		float getRCadj(const std::string& resName, const std::string& aName);
		float getPrevRCadj(const std::string& prev_rName, const std::string& aName);
		float getNextRCadj(const std::string& next_rName, const std::string& aName);
		float getWeight(const std::string& Name, const std::string& aName);

		//preset the args form command line
		void setup_defaults();

		void mkdir_pred(const std::string& d);// create a directory for prediction results

		////SPARTA originals
		std::string SPARTA_DIR, PRED_DIR, TAB_DIR, SHIFT_DIR, PDB_DIR, EXCLUDED;
		std::string slash_char;
		char buf[300], lbuf[1000];

		std::string inName; // input PDB coordinates file name
		std::string inNames; // name list of multiple input PDB coordinates files
		PDB inPDB;
		GDB inTab; //inTab, input table file from PDB coordinates file

		std::string sumName; // output summary name
		std::string sourceName; // output summary name

		float** U_ANGLES;
		float** U_RING_SHIFTS;
		std::string** U_NAME ;
		float *U_HN_HB, *U_HA_HB, *U_CO_HB;// *HB_SHIFT_HN;

		std::string tripFileName;
		GDB TRIPLET_Tab; //triplet table file

		std::string weightFileName; // file name for wieght talble of score function
		GDB WEIGHT_Tab; // talble for weights of score function

		std::string homoFileName; // file name for sequence homology talble
		GDB HOMO_Tab; // homology talble

		std::string rcFileName; //file name for random coil chemical shift talble
		GDB RC_Tab; // random coil chemical shift talble

		std::string adjFileName; //file name for the adjustment of random coil chemical shift talble
		GDB ADJ_Tab; // "Random Coil Adjustment Table"

		std::string prevFileName; // file name for the table of "Random Coil Adjustments for Previous Residue Type"
		GDB PREV_Tab; // talble for "Random Coil Adjustments for Previous Residue Type"

		std::string nextFileName; // file name for the table of "Random Coil Adjustments for Next Residue Type"
		GDB NEXT_Tab; // talble for "Random Coil Adjustments for Next Residue Type"

		std::string fitFileName; // file name for the table of "fitting parameters"
		GDB FIT_Tab; // talble for "fitting parameters"


		int firstRes, lastRes; // First/last RESID to use for prediction

		int r1,rN; // first and last RESID in seqList

		float tVal; // Max similarity score threshold, not used in the program

		typedef std::map< int, std::string > ResidList;
		ResidList residList; // one-letter amino acid residue list from input PDB coordinates file
		std::string sequence; // one-letter amino acid residue list from input PDB coordinates file in the format of one std::string

		typedef boost::unordered_map< int, std::string > AtomNameList;
		AtomNameList aN/*, aN_ALL*/; // Backbone atom list used by program "N HA C CA CB H"

		int matchCount; //Max Match Count per query triplet

		std::string pdbListName; //table file name for names of candidate proteins

		//boost::unordered_map<std::string, boost::unordered_map<std::string, float> > Fitting;

		std::string AAlist; // amino acid list (with a sequence allowed by ANN)

		typedef std::map< std::string, utility::vector0< float > > BlosumMatrix;
		BlosumMatrix BLOSUM_62; // BLOSUM 62 matrix

		ANN::ANN_Matrix ANN_IN_MTX; // input matrix for neural netwrok calculation

		typedef std::map< std::string, ANN::ANN_Matrix > Atom2ANN_MatrixMap;
		Atom2ANN_MatrixMap ANN_CS_OUTPUT_FULL; // input matrix from neural netwrok calculation, indexed by atom name, resID and prediction

		typedef std::map< int, float > AngleMap;
		AngleMap CHI2_ANGLES, OMEGA_ANGLES;

		typedef std::map< int, std::map< std::string,float > > SurfaceExposureMap;
		SurfaceExposureMap SURFACE_EXPOSURE; //indexed by resID, atomName

		boost::unordered_map< std::string, ANN> SPARTA_ANN;
		boost::unordered_map< std::string, boost::unordered_map< std::string, PHIPSI_ERR_SURF> > SPARTA_ERR_SURF; //indexed by AA, atomName, phi, psi
	};
};

} // sparta
} // protocols

#endif
