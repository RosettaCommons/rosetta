
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/ScoreFunction.hh>
#include <sstream>

// Utility Headers
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <protocols/flxbb/SelectResiduesByLayer.hh>
#include <protocols/moves/Mover.hh>
#include <numeric/xyzVector.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

template <typename T>
std::string vec2str(const T& vec){
		std::ostringstream out;
		for(typename T::const_iterator it = vec.begin(); it != vec.end();  ++it) {
			out << *it << "\t";
		}
		return out.str();
};

static basic::Tracer TR("cristallizable_surface_design");
using namespace protocols::flxbb;
using namespace core;
using namespace core::pose;
using namespace utility;
using namespace basic::options;
using namespace basic::options::OptionKeys;



class ThisApplication {
    public:
        ThisApplication(){};
        static void register_options();
};

OPT_KEY( String , pdb)
void ThisApplication::register_options() {
    NEW_OPT( pdb, "pdb file", "");
}

class CristallizableSurfaceDesign  : public  protocols::moves::Mover { 
		typedef protocols::moves::MoverOP MoverOP;
		typedef utility::vector1<Size> VecSize;
		typedef std::map<Size, VecSize> NeighborMap;
		typedef core::scoring::ScoreFunction ScoreFunction;
		typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
		typedef core::conformation::Residue	Residue;
		typedef core::conformation::ResidueOP	ResidueOP;
		public:
				/// @brief default constructor
				CristallizableSurfaceDesign()
	 			{ 
				};

				/// @brief value constructor
//				CristallizableSurfaceDesign(
//				ScoreFunctionOP const sfxnd,
//				ScoreFunctionOP const sfxnr,
//				Size const ncycle = 3,
//				bool  relax = true );

				/// @brief copy constructor
	//			CristallizableSurfaceDesign( CristallizableSurfaceDesign const & rval ) {

		//	 	};

				/// @brief destructor
			 virtual ~CristallizableSurfaceDesign() { };


		public: // virtual constructors

				/// @brief clone this object
		//		virtual
		//		MoverOP clone() const;


				/// @brief create this type of object
				virtual
				MoverOP 
				fresh_instance() const { return new CristallizableSurfaceDesign; }
		public:

				virtual
				std::string get_name() const { return "CristallizableDesign" ;}				
	//			virtual
				void apply(Pose& p) {
					if(!initialized_)
							initialize();
        	VecSize surf_res_indeces = selector_->compute(p, p.secstruct());
					NeighborMap neighbors;
					for(VecSize::const_iterator it1 = surf_res_indeces.begin(); it1 != surf_res_indeces.end() - 1; ++it1) {
						for(VecSize::const_iterator it2 = it1 + 1; it2 != surf_res_indeces.end(); ++it2) {
						//get the distance between the centroids
							Vector v1  = p.residue(*it1).xyz("CA");
							Vector v2  = p.residue(*it2).xyz("CA");
							Real distance = v1.distance(v2);
							if(distance < neighbor_residue_distance_threshold_)
								neighbors.insert(std::make_pair<Size,Size>(*it1, *it2));
						}
					}
					
				}

	  private:
				void initialize() { 
						selector_ = new SelectResiduesByLayer("surface");
						initialized_ = true;
				}

		private:
				SelectResiduesByLayerOP selector_;
			//	ScoreFunctionOP sfxn_design_;
			//	ScoreFunctionOP sfxn_relax_;

	  		Real	neighbor_residue_distance_threshold_;
				Real core_;
				Real surface_;
				bool relax_;
				bool initialized_;
				Size ncycle_;
};

int
main ( int argc, char* argv[] ){
		ThisApplication::register_options();
		devel::init(argc, argv);
		TR << "Reaading File " << option[pdb]() << std::endl;
		core::pose::PoseOP p =  import_pose::pose_from_pdb(option[pdb]());
		core::scoring::dssp::Dssp dssp(*p);
		dssp.insert_ss_into_pose(*p);
		
		protocols::moves::MoverOP csd;
		csd = new CristallizableSurfaceDesign() ;
		csd->apply(*p);



    return 0;
}
