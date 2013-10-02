// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.



#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <basic/options/option.hh>

#include <core/init/init.hh>
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>
#include <core/id/AtomID_Map.hh>
#include <protocols/moves/Mover.hh>

#include <basic/Tracer.hh>

//JD headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/import_pose/import_pose.hh>

#include <string>
#include <ObjexxFCL/Fmath.hh>
#include <ObjexxFCL/format.hh>
#include <numeric/numeric.functions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>

// option key includes

static basic::Tracer TR("decoy_features");

using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;

///local mover for testing purposes
class JDmover : public protocols::moves::Mover {


public:
	JDmover()
	{
	  input_pose_count_ = 0;
		input_pose_length_ = 0;
		input_pose_sequence_ = "";
	}

	std::string get_name() const { return "BBinMover"; }

	virtual ~JDmover(){};

  core::Real
  periodic_range(
    core::Real a,
    core::Real x
  )
  {
		using namespace numeric;
    core::Real const halfx = 0.5f * x;
    return ( ( a >= halfx || a < -halfx ) ? mod( mod( a, x ) + ( x + halfx ), x ) - halfx : a );
  }

  std::string ABGEO(core::Real phi, core::Real psi, core::Real omega) {

    periodic_range( phi  , 360.0 );  //does this get applied to phi??
    periodic_range( psi  , 360.0 );
    periodic_range( omega, 360.0 );
    if ( (int)std::abs( (int)omega ) < 90 ) {
      return "O";
    } else if ( phi >= 0.0 ) {
      if ( -100 < psi && psi <= 100 ) {
        return "G"; // alpha-L
      } else {
        return "E"; // E
      }
    } else {
      if ( -125 < psi && psi <= 50 ) {
        return "A"; // helical
      } else {
        return "B"; // extended
      }
    }

	}

	virtual
	void
	apply( core::pose::Pose & pose ){

		// initialize secondary structure from DSSP.
		core::scoring::dssp::Dssp dssp_obj( pose );
		dssp_obj.insert_ss_into_pose( pose );
		if (input_pose_length_ <= 0) {
			input_pose_length_ = pose.total_residue();
			input_pose_sequence_ = pose.sequence();
		} else if (input_pose_length_ != pose.total_residue()){
			utility_exit_with_message("input pose lengths are not equal");
		}
		input_pose_count_++;

		for ( unsigned int i = 1; i <= pose.total_residue(); ++i ) {
			if (pose.secstruct(i) == 'E') {
				ss_E_counts_[i]++;
			} else if (pose.secstruct(i) == 'H') {
				ss_H_counts_[i]++;
			} else if (pose.secstruct(i) == 'L') {
				ss_L_counts_[i]++;
			}
			std::string bbin = ABGEO( pose.phi(i), pose.psi(i),pose.omega(i) );
			if (bbin == "A") {
				bbin_A_counts_[i]++;
			} else if (bbin == "B") {
				bbin_B_counts_[i]++;
			} else if (bbin == "E") {
				bbin_E_counts_[i]++;
			} else if (bbin == "G") {
				bbin_G_counts_[i]++;
			} else if (bbin == "O") {
				bbin_O_counts_[i]++;
			}
      for ( unsigned int j = 1; j <= pose.total_residue(); ++j ) {

				if (std::abs((int)(i-j)) < 3) continue; // sequence sep of 3

				core::Real distance = 0.0;
				if (pose.residue_type( i ).name3() == "GLY") {
					if (pose.residue_type( j ).name3() == "GLY") {
						distance = pose.residue(i).xyz("CA").distance( pose.residue(j).xyz("CA") );
					} else {
						distance = pose.residue(i).xyz("CA").distance( pose.residue(j).xyz("CB") );
					}
				} else if (pose.residue_type( j ).name3() == "GLY") {
					distance = pose.residue(i).xyz("CB").distance( pose.residue(j).xyz("CA") );
				} else {
					distance = pose.residue(i).xyz("CB").distance( pose.residue(j).xyz("CB") );
				}
				if ( distance > 8.0 ) continue;

				std::pair<Size,Size> contact(i,j);
				contact_counts_[contact]++;

      } //  for ( unsigned int j = 1; j <= pose.total_residue(); ++j )
    }   // for ( unsigned int i = 1; i <= pose.total_residue(); ++i )
	}



	virtual
	protocols::moves::MoverOP
	fresh_instance() const {
		return new JDmover;
	}

	virtual
	bool
	reinitialize_for_each_job() const { return false; }

	virtual
	bool
	reinitialize_for_new_input() const { return false; }


	Size get_input_pose_length() {
		return input_pose_length_;
	}

	void print_features() {
		using namespace ObjexxFCL;
		// output SS
		std::cout << "SEQ " << input_pose_sequence_ << std::endl;
		for ( unsigned int i = 1; i <= input_pose_length_; ++i ) {
			std::cout << "DS " << (i-1) << " E: " << format::F(8, 6, (core::Real)ss_E_counts_[i]/(core::Real)input_pose_count_) <<
				" H: " << format::F(8, 6, (core::Real)ss_H_counts_[i]/(core::Real)input_pose_count_) <<
				" L: " << format::F(8, 6, (core::Real)ss_L_counts_[i]/(core::Real)input_pose_count_) << std::endl;
		}
		// output ABEGO
		for ( unsigned int i = 1; i <= input_pose_length_; ++i ) {
			std::cout << "DBB " << (i-1) << " A: " << format::F(8, 6, (core::Real)bbin_A_counts_[i]/(core::Real)input_pose_count_) <<
			" B: " << format::F(8, 6, (core::Real)bbin_B_counts_[i]/(core::Real)input_pose_count_) <<
			" G: " << format::F(8, 6, (core::Real)bbin_G_counts_[i]/(core::Real)input_pose_count_) <<
			" E: " << format::F(8, 6, (core::Real)bbin_E_counts_[i]/(core::Real)input_pose_count_) <<
			" O: " << format::F(8, 6, (core::Real)bbin_O_counts_[i]/(core::Real)input_pose_count_) << std::endl;
		}
		// output contacts
		std::map<std::pair<Size,Size>, Size>::iterator it;
		for ( unsigned int i = 1; i <= input_pose_length_; ++i ) {
			for ( unsigned int j = i+1; j <= input_pose_length_; ++j ) {
				std::pair<Size,Size> p(i,j);
				if (contact_counts_.find(p) != contact_counts_.end()) {
					std::cout << "DC " << i << " " << j << " " << format::F(8, 6, (core::Real)contact_counts_[p]/(core::Real)input_pose_count_) << std::endl;
				}
			}
		}
	}


private:
	std::map<std::pair<Size,Size>, Size> contact_counts_;
	Size input_pose_count_;
	Size input_pose_length_;
	std::string input_pose_sequence_;

	std::map<Size,Size> ss_E_counts_;
	std::map<Size,Size> ss_H_counts_;
	std::map<Size,Size> ss_L_counts_;
	std::map<Size,Size> bbin_A_counts_;
	std::map<Size,Size> bbin_B_counts_;
	std::map<Size,Size> bbin_E_counts_;
	std::map<Size,Size> bbin_G_counts_;
	std::map<Size,Size> bbin_O_counts_;

};

typedef utility::pointer::owning_ptr< JDmover > JDmoverOP;




///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

	core::init::init(argc, argv);

	JDmoverOP jd_mover(new JDmover);

	try {
		protocols::jd2::JobDistributor::get_instance()->go(jd_mover);
	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cerr << "Exception: " << std::endl;
		excn.show( std::cerr );
		std::cout << "Exception: " << std::endl;
		excn.show( std::cout ); //so its also seen in a >LOG file
	}

	if (option[in::file::native].user()) {
    core::pose::PoseOP nativePose = new core::pose::Pose;
    core::import_pose::pose_from_pdb(*nativePose, option[in::file::native]());
		core::scoring::dssp::Dssp dssp_obj( *nativePose );
		dssp_obj.insert_ss_into_pose( *nativePose );
		// native ss
		std::cout << "NS " << nativePose->secstruct() << std::endl;
		if (nativePose->total_residue() != jd_mover->get_input_pose_length()) {
			utility_exit_with_message("native length must equal input models");
		}
		// native bbin
		std::cout << "NBB ";
		for ( unsigned int i = 1; i <= nativePose->total_residue(); ++i ) {
			std::cout << jd_mover->ABGEO( nativePose->phi(i), nativePose->psi(i), nativePose->omega(i) );
		}
		std::cout << std::endl;
		// native contacts
		for ( unsigned int i = 1; i <= nativePose->total_residue(); ++i ) {
			for ( unsigned int j = 1; j <= nativePose->total_residue(); ++j ) {

				if (std::abs((int)(i-j)) < 3) continue; // sequence sep of 3

				core::Real distance = 0.0;
				if (nativePose->residue_type( i ).name3() == "GLY") {
					if (nativePose->residue_type( j ).name3() == "GLY") {
						distance = nativePose->residue(i).xyz("CA").distance( nativePose->residue(j).xyz("CA") );
					} else {
						distance = nativePose->residue(i).xyz("CA").distance( nativePose->residue(j).xyz("CB") );
					}
				} else if (nativePose->residue_type( j ).name3() == "GLY") {
					distance = nativePose->residue(i).xyz("CB").distance( nativePose->residue(j).xyz("CA") );
				} else {
					distance = nativePose->residue(i).xyz("CB").distance( nativePose->residue(j).xyz("CB") );
				}
				if ( distance > 8.0 ) continue;
				std::cout << "NC " << i << " " << j << std::endl;
      } //  for ( unsigned int j = 1; j <= pose.total_residue(); ++j )
    }   // for ( unsigned int i = 1; i <= pose.total_residue(); ++i )

	}

	jd_mover->print_features();

	TR << "*********************successful completion**************************" << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

}

