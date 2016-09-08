// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/rjha/ZincMatchFilter.cc
/// @brief This application is designed to help filter through the matches from RosettaMatch.  It calculates some of the geometry about the metal atom (currently the metal-ligand atom distances and the tetrahedral coordinations six ligand-metal-ligand angles (all of which should be 109.5).  It also outputs a sum of squares for those values.
/// @author Steven Lewis and Bryan Der

// Unit Headers
#include <devel/metal_interface/FindClosestAtom.hh>
// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/conformation/Residue.hh>
#include <protocols/moves/Mover.hh>
// Numeric Headers
#include <numeric/xyzVector.hh>
//#include <numeric/xyzVector.io.hh>
#include <numeric/conversions.hh> //degrees-radians
#include <numeric/xyz.functions.hh>
// Utility Headers
#include <devel/init.hh>
#include <core/types.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <protocols/jd2/JobDistributor.hh>
//Auto Headers
#include <core/import_pose/import_pose.hh>


// C++ headers
#include <math.h>

//tracers
using basic::Error;
using basic::Warning;
using basic::T;
static THREAD_LOCAL basic::Tracer TR( "apps.pilot.rjha.ZincMatchFilter" );

typedef numeric::xyzVector<core::Real> point;
typedef point axis;

using namespace core;

/// @brief
class ZincMatchFilter : public protocols::moves::Mover {
public:
  ZincMatchFilter()
    : distance_weight_(1.0), tetrahedral_weight_(1.0), angle_weight_(1.0), dihedral_weight_(1.0), distance_stdev_(0.12), tetrahedral_stdev_(20.0), angle_stdev_(20.0), dihedral_stdev_(20.0)
  {
  }
  virtual ~ZincMatchFilter(){};

  //setters
  void set_distance_weight( Real distance_weight ){
    distance_weight_ = distance_weight;
    return;
  }
  void set_angle_weight( Real angle_weight ){
    angle_weight_ = angle_weight;
    return;
  }
  void set_dihedral_weight( Real dihedral_weight ){
    distance_weight_ = dihedral_weight;
    return;
  }
  void set_tetrahedral_weight( Real tetrahedral_weight ){
    tetrahedral_weight_ = tetrahedral_weight;
    return;
  }
  void set_distance_stdev( Real distance_stdev ){
    distance_stdev_ = distance_stdev;
    return;
  }
  void set_angle_stdev( Real angle_stdev ){
    angle_stdev_ = angle_stdev;
    return;
  }
  void set_dihedral_stdev( Real dihedral_stdev ){
    distance_stdev_ = dihedral_stdev;
    return;
  }
  void set_tetrahedral_stdev( Real tetrahedral_stdev ){
    tetrahedral_stdev_ = tetrahedral_stdev;
    return;
  }

  //getters
  Real get_distance_weight()   { return distance_weight_;    }
  Real get_angle_weight()      { return angle_weight_;       }
  Real get_dihedral_weight()   { return dihedral_weight_;    }
  Real get_tetrahedral_weight(){ return tetrahedral_weight_; }
  Real get_distance_stdev()    { return distance_stdev_;     }
  Real get_angle_stdev()       { return angle_stdev_;        }
  Real get_dihedral_stdev()    { return dihedral_stdev_;     }
  Real get_tetrahedral_stdev() { return tetrahedral_stdev_;  }


  virtual
  void
  apply( core::pose::Pose & match ){

    std::string match_name = match.pdb_info()->name();
    Real number_of_ligands = match.size() - 1;

    //desired info
    Real distance_diffs = 0;
    Real angle_diffs = 0;
    Real dihedral_diffs = 0;
    Real tetrahedral_diffs = 0;

    utility::vector1< Real > distances;
    utility::vector1< Real > angles;
    utility::vector1< Real > dihedrals;
    utility::vector1< Real > tetrahedral_angles;

    //required info
    utility::vector1< point > atom_xyz;
    utility::vector1< Size > seqpos;
    utility::vector1< std::string > atom_name;

    point zinc = match.residue( match.size() ).atom(1).xyz();

    T(match_name) << "This match has " << number_of_ligands << " liganding residues." << std::endl;

    //get required info, excluding zinc
    for ( Size i(1); i <= number_of_ligands; ++i ) {
			std::string atomname = devel::metal_interface::find_closest_atom( match.residue(i), zinc );

      T(match_name) << "Checking that atom_name is NE2, ND1, SG, O" << std::endl;
      if(!(atomname == " NE2"
				|| atomname == " ND1"
				|| atomname == " SG "
				|| atomname == " OD1"
				|| atomname == " OD2"
				|| atomname == " OE1"
				|| atomname == " OE2"
				)){ return; }
      seqpos.push_back( match.pdb_info()->number(i) ); // gets residue number according to .pdb file, not renumbered starting at 1
			atom_name.push_back(atomname);
      atom_xyz.push_back( match.residue(i).atom(atom_name[i]).xyz() );
    }

    //calculate distances
    for ( Size i(1); i <= number_of_ligands; ++i ) {
      Real dist(0);
      Real dist_diff(0);
      if( atom_name[i] == " NE2" || atom_name[i] == " ND1" || atom_name[i] == " OD1" || atom_name[i] == " OD2" || atom_name[i] == " OE1" || atom_name[i] == " OE2" ) {
				dist = atom_xyz[i].distance( zinc );
				dist_diff = fabs(dist - 2.05);
				distance_diffs += dist_diff*dist_diff;
				distances.push_back(dist);
      }
      else if ( atom_name[i] == " SG " ) {
				dist = atom_xyz[i].distance( zinc );
				dist_diff = fabs(dist - 2.33);
				distance_diffs += dist_diff*dist_diff;
				distances.push_back(dist);
      }
			T(match_name) << "Distance " <<  match.residue(i).name3() << seqpos[i] << " res " << i << " distance=" << dist << " diff=" << dist_diff << std::endl;
    }


    //calculate His/Asp/Glu dihedrals, don't care about Cys
    Size number_of_dihedrals( 0 );
    for ( Size i(1); i <= number_of_ligands; ++i ) {
      Real dihed;
      Real dihed_diff;
      if( atom_name[i] == " NE2" ) {
				dihed = numeric::dihedral_degrees( zinc, atom_xyz[i], match.residue(i).atom("CD2").xyz(), match.residue(i).atom("CG").xyz() );
				dihedrals.push_back(dihed);
				dihed_diff = fabs(dihed) - 180.0;
				dihedral_diffs += fabs(dihed_diff)*fabs(dihed_diff);
				number_of_dihedrals++;
				T(match_name) << "Dihedral " << match.residue(i).name3() << seqpos[i] << " res " << i << " dihedral=" << dihed << " diff=" << dihed_diff << std::endl;
      }
      else if( atom_name[i] == " ND1" ) {
				dihed = numeric::dihedral_degrees( zinc, atom_xyz[i], match.residue(i).atom("CG").xyz(), match.residue(i).atom("CB").xyz() );
				dihedrals.push_back(dihed);
				dihed_diff = fabs(dihed) - 0.0;
				dihedral_diffs += fabs(dihed_diff)*fabs(dihed_diff);
				number_of_dihedrals++;
				T(match_name) << "Dihedral " << match.residue(i).name3() << seqpos[i] << " res " << i << " dihedral=" << dihed << " diff=" << dihed_diff << std::endl;
      }
      else if ( atom_name[i] == " SG " ) {
				T(match_name) << match.residue(i).name3() << seqpos[i] << " res " << i << " does not have dihedral we care about" << std::endl;
				dihedrals.push_back(0);
      }
      else if( atom_name[i] == " OD1" || atom_name[i] == " OD2") {
				dihed = numeric::dihedral_degrees( zinc, atom_xyz[i], match.residue(i).atom("CG").xyz(), match.residue(i).atom("CB").xyz() );
				dihedrals.push_back(dihed);
				Real dihed_diff_180 = fabs(fabs(dihed) - 180.0);
				Real dihed_diff_0 = fabs(dihed);
				if( dihed_diff_180 < dihed_diff_0) { dihed_diff = dihed_diff_180; }
				else { dihed_diff = dihed_diff_0; }
				dihedral_diffs += dihed_diff*dihed_diff;
				number_of_dihedrals++;
				T(match_name) << "Dihedral " << match.residue(i).name3() << seqpos[i] << " res " << i << " dihedral=" << dihed << " diff=" << dihed_diff << std::endl;
      }
      else if( atom_name[i] == " OE1" || atom_name[i] == " OE2") {
				dihed = numeric::dihedral_degrees( zinc, atom_xyz[i], match.residue(i).atom("CD").xyz(), match.residue(i).atom("CG").xyz() );
				dihedrals.push_back(dihed);
				Real dihed_diff_180 = fabs(fabs(dihed) - 180.0);
				Real dihed_diff_0 = fabs(dihed);
				if( dihed_diff_180 < dihed_diff_0) { dihed_diff = dihed_diff_180; }
				else { dihed_diff = dihed_diff_0; }
				dihedral_diffs += dihed_diff*dihed_diff;
				number_of_dihedrals++;
				T(match_name) << "Dihedral " << match.residue(i).name3() << seqpos[i] << " res " << i << " dihedral=" << dihed << " diff=" << dihed_diff << std::endl;
      }
    }


    //calculate angles
    Size number_of_angles( 0 );
    for ( Size i(1); i <= number_of_ligands; ++i ) {
			Real angle(0);
			Real angle_diff(0);
			if(atom_name[i] == " SG ")  {
				angle = numeric::conversions::degrees(angle_of(zinc, atom_xyz[i], match.residue(i).atom("CB").xyz()));
				angle_diff = fabs(angle - 109.5);
				angles.push_back(angle);
			}
			else if(atom_name[i] == " ND1")  {
				angle = numeric::conversions::degrees(angle_of(zinc, atom_xyz[i], match.residue(i).atom("CG").xyz()));
				angle_diff = fabs(angle - 125.0);
				angles.push_back(angle);
			}
			else if(atom_name[i] == " NE2")  {
				angle = numeric::conversions::degrees(angle_of(zinc, atom_xyz[i], match.residue(i).atom("CD2").xyz()));
				angle_diff = fabs(angle - 125.0);
				angles.push_back(angle);
			}
			else if(atom_name[i] == " OD1" || atom_name[i] == " OD2" )  {
				angle = numeric::conversions::degrees(angle_of(zinc, atom_xyz[i], match.residue(i).atom("CG").xyz()));
				angle_diff = fabs(angle - 125.0);
				angles.push_back(angle);
			}
			else if(atom_name[i] == " OE1" || atom_name[i] == " OE2" )  {
				angle = numeric::conversions::degrees(angle_of(zinc, atom_xyz[i], match.residue(i).atom("CD").xyz()));
				angle_diff = fabs(angle - 125.0);
				angles.push_back(angle);
			}

			angle_diffs += angle_diff*angle_diff;
			T(match_name) << "Angle of " << seqpos[i] << " angle=" << angle << " diff=" << angle_diff << std::endl;
			number_of_angles++;
    }


    //calculate tetrahedral angles
    Size number_of_tetr_angles( 0 );
    for ( Size i(1); i <= number_of_ligands; ++i ) {
      for ( Size j(i+1); j <= number_of_ligands; ++j ) {
				Real tetr_angle = numeric::conversions::degrees(angle_of(atom_xyz[i], zinc, atom_xyz[j]));
				Real tetr_angle_diff = fabs(tetr_angle - 109.5);
				tetrahedral_diffs += tetr_angle_diff*tetr_angle_diff;
				T(match_name) << "TetrAngle of " << seqpos[i] << " zinc " << seqpos[j] << " tetr_angle=" << tetr_angle << " diff=" << tetr_angle_diff << std::endl;
				number_of_tetr_angles++;
      }
    }


    Real weighted_distance_score = ( distance_weight_ * distance_diffs ) / ( distance_stdev_ * distance_stdev_ * number_of_ligands );
    Real weighted_angle_score = ( angle_weight_ * angle_diffs ) / ( angle_stdev_ * angle_stdev_ * number_of_angles );
    Real weighted_dihedral_score( 0 );
    Real weighted_tetrahedral_score = ( tetrahedral_weight_ * tetrahedral_diffs ) / ( tetrahedral_stdev_ * tetrahedral_stdev_ * number_of_tetr_angles );

		T(match_name) << "number of dihedrals: " << number_of_dihedrals << std::endl;

		if ( number_of_dihedrals > 0 ) {
      weighted_dihedral_score += ( dihedral_weight_ * dihedral_diffs ) / ( dihedral_stdev_ * dihedral_stdev_ * number_of_dihedrals );
    }

    for ( Size i(1); i <= number_of_ligands; ++i ) {
			T(match_name) << "SUMMARY: " << i << " " << seqpos[i] << " " << atom_name[i] << "\t" << distances[i] << "\t" << angles[i] << "\t" << dihedrals[i] << std::endl;
		}
			T(match_name) << "SUMMARY_ALL " << weighted_distance_score << "\t" << weighted_angle_score << "\t" << weighted_dihedral_score << "\t" << weighted_tetrahedral_score << std::endl;

		T(match_name) << "Weighted_Score_Distance: " << weighted_distance_score << std::endl;
    T(match_name) << "Weighted_Score_Angle:    " << weighted_angle_score << std::endl;
    T(match_name) << "Weighted_Score_Dihedral: " << weighted_dihedral_score << std::endl;
    T(match_name) << "Weighted_Score_Tetrahed: " << weighted_tetrahedral_score << std::endl;
    T(match_name) << "Weighted_Score_TOTAL:    " << weighted_distance_score + weighted_tetrahedral_score + weighted_angle_score + weighted_dihedral_score << std::endl;

    return;
  }


	virtual
	std::string
	get_name() const { return "ZincMatchFilter"; }


private:
  Real distance_weight_;
  Real tetrahedral_weight_;
  Real angle_weight_;
  Real dihedral_weight_;

  Real distance_stdev_;
  Real tetrahedral_stdev_;
  Real angle_stdev_;
  Real dihedral_stdev_;

};

typedef utility::pointer::owning_ptr< ZincMatchFilter > ZincMatchFilterOP;

int main( int argc, char* argv[] )
{
	try {

  devel::init(argc, argv);
  protocols::jd2::JobDistributor::get_instance()->go(new ZincMatchFilter);

  TR << "************************d**o**n**e**************************************" << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

  return 0;
}

