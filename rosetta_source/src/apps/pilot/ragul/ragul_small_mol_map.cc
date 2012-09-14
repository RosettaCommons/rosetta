// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Ragul Gowthaman

#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>

# include <fstream>
# include <ostream>
# include <string>
# include <sstream>
// Protocol Headers
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/backrub/BackrubMover.hh>
#include <protocols/moves/MonteCarlo.hh>

// Core Headers
#include <core/conformation/Residue.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/pose/PDBInfo.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <basic/options/util.hh>
#include <basic/options/after_opts.hh>
#include <core/id/AtomID_Map.hh>
#include <basic/Tracer.hh>
#include <core/scoring/Energies.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/AtomType.hh>
#include <core/types.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <basic/options/option_macros.hh>

//vector
#include <vector>
using std::vector;



// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

using namespace core;
using namespace basic::options;
using namespace core::scoring;
using namespace core::optimization;
using namespace basic::options::OptionKeys;
using namespace std;
core::Real Find_Intersect(core::SSize const & phiAngle, core::SSize const & psiAngle, core::Real const & atomX, core::Real const & atomY, core::Real const & atomZ, core::Real const & atom_radius);

	//stl-map for MaxDist, Xpoint, Ypoint, Zpoint
//	std::map<core::SSize, std::map<core::SSize, core::Real> > MaxDist;
//	std::map<core::SSize, std::map<core::SSize, core::Real> > Xpoint;
//	std::map<core::SSize, std::map<core::SSize, core::Real> > Ypoint;
//	std::map<core::SSize, std::map<core::SSize, core::Real> > Zpoint;

// vector for  phi, psi, rho ,x,y,z
   std::vector<std::vector<core::Real> > MaxDist;
   std::vector<std::vector<core::Real> > Xpoint;
   std::vector<std::vector<core::Real> > Ypoint;
   std::vector<std::vector<core::Real> > Zpoint;
   vector<core::Real> phiAngle;
   vector<core::Real> psiAngle;


int main( int argc, char * argv [] )
{
        const core::Real PI = numeric::NumericTraits<Real>::pi();
        const core::Real RADS_PER_DEG = PI / 180.;
        const core::Real DEGS_PER_RAD = 180./PI;

	devel::init(argc, argv);
	// std::string const output_tag = option[ OptionKeys::out::output_tag ]();
	pose::Pose input_pose;
	//read in pdb file from command line
	std::string const input_pdb_name ( basic::options::start_file() );
	core::import_pose::pose_from_pdb( input_pose, input_pdb_name );

	core::Real PcomX, PcomY, PcomZ;
	core::Real comx=0, comy=0, comz=0;
	core::Real ScomX=0, ScomY=0, ScomZ=0;
	PcomX=44.659; PcomY=14.092; PcomZ=28.772;

	core::Size lig_res_num = 0;
	for ( int j = 1, resnum = input_pose.total_residue(); j <= resnum; ++j ) {
                 if (!input_pose.residue(j).is_protein()){
                         lig_res_num = j;
                         break;
		}
	}
	if (lig_res_num != 0){
	conformation::Residue const & curr_rsd=input_pose.conformation().residue(lig_res_num);
        ScomX=0, ScomY=0, ScomZ=0;
	//SET-UP FILE FOR PRINTING SMALL-MOL FROM ORIGIN
	std::filebuf f1;
  	std::stringstream filename1;
	filename1<<"SmallPrint"<<".pdb";
	f1.open (filename1.str().c_str(),std::ios::out);
	std::ostream os1(&f1);
	std::string f1_info;
       	std::stringstream f1_tmp;

	for(Size i = 1, i_end = curr_rsd.nheavyatoms(); i <= i_end; ++i) {
	ScomX +=  curr_rsd.atom(i).xyz()(1);
        ScomY +=  curr_rsd.atom(i).xyz()(2);
        ScomZ +=  curr_rsd.atom(i).xyz()(3);
	}
	ScomX = ScomX/curr_rsd.nheavyatoms();
	ScomY = ScomY/curr_rsd.nheavyatoms();
	ScomZ = ScomZ/curr_rsd.nheavyatoms();

	comx = ScomX;
	comy = ScomY;
	comz = ScomZ;

	core::Real phiAngle, psiAngle;
        core::Real atomX,atomY,atomZ,atom_radius;
	for(Size i = 1, i_end = curr_rsd.nheavyatoms(); i <= i_end; ++i) {
	atomX = 0;  atomY = 0;   atomZ = 0;   atom_radius = 0;
	atomX =  curr_rsd.atom(i).xyz()(1)-comx;
        atomY =  curr_rsd.atom(i).xyz()(2)-comy;
        atomZ =  curr_rsd.atom(i).xyz()(3)-comz;
	atom_radius = curr_rsd.atom_type(i).lj_radius();

	f1_tmp<<"HETATM   "<<std::setw(2)<<i<<"  C   MAP X   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<atomX+comx<<std::setw(8)<<std::fixed<<std::setprecision(3)<<atomY+comy<<std::setw(8)<<std::fixed<<std::setprecision(3)<<atomZ+comz<<std::endl;
	f1_info += f1_tmp.str();
        f1_tmp.str(std::string());

  //2d vector array
/*
   for (core::Real i = 1.0; i > -1.05; i+= -0.05){
   phiAngle.push_back(i);
   }

   for (core::Real j = -180; j <= 180; j+=5){
   psiAngle.push_back(j);
   }

   for( vector::iterator ith = phiAngle.begin(); ith != phiAngle.end(); ith++ ){
   for( vector::iterator its = psiAngle.begin(); its != psiAngle.end(); its++ ){
   MaxDist[phiAngle][psiAngle] = 0;
   }
*/




	for (core::Real Cphi = 1.0; Cphi > -1.05; Cphi+= -0.05){
	for (core::Real psiAngle = -180; psiAngle <= 180; psiAngle+=5){
        MaxDist[phiAngle][psiAngle] = 0;
        std::cout<<MaxDist[phiAngle][psiAngle]<<"..hi.."<<std::endl;
   }
}



	for (core::Real Cphi = 1.0; Cphi > -1.05; Cphi+= -0.05){
	for (core::SSize psiAngle = -180; psiAngle <= 180; psiAngle+=5){
	phiAngle = acos(Cphi) * DEGS_PER_PI;
	Find_Intersect(phiAngle,psiAngle,atomX,atomY,atomZ,atom_radius);

	      }
	    }
	  }

	f1_tmp<<"HETATM   "<<std::setw(2)<<0<<"  K   COM X   0    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<comx<<std::setw(8)<<std::fixed<<std::setprecision(3)<<comy<<std::setw(8)<<std::fixed<<std::setprecision(3)<<comz<<std::endl;
	f1_info += f1_tmp.str();
	f1_tmp.str(std::string());

	std::filebuf f2;
	std::stringstream filename2;
	filename2<<"SmallPrint"<<".txt";
	f2.open (filename2.str().c_str(),std::ios::out);
        std::ostream os2(&f2);
       	std::string f2_info;
       	std::stringstream  f2_tmp;

	for (core::Real Cphi = 1.0; Cphi > -1.05; Cphi+= -0.05){
        for (core::SSize psiAngle = -180; psiAngle <= 180; psiAngle+=5){
	phiAngle = acos(Cphi) * DEGS_PER_RAD;
	if (MaxDist[phiAngle][psiAngle] != 0)
	{
	f1_tmp<<"ATOM      1  I   SUR     1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<Xpoint[phiAngle][psiAngle]+comx<<std::setw(8)<<std::fixed<<std::setprecision(3)<<Ypoint[phiAngle][psiAngle]+comy<<std::setw(8)<<std::fixed<<std::setprecision(3)<<Zpoint[phiAngle][psiAngle]+comz<<std::endl;
	f1_info += f1_tmp.str();
        f1_tmp.str(std::string());

	f2_tmp<<phiAngle<<" \t " <<psiAngle<< "\t" << MaxDist[phiAngle][psiAngle] <<std::endl;
	f2_info += f2_tmp.str();
        f2_tmp.str(std::string());
	}
        }
   }
	os1<<f1_info;
        os2<<f2_info;
	f1.close();
	f2.close();
 }

//reading pocket_fingerprint file
ifstream inFile("PocketPrint.txt");
std::string lineread;
std::string Line;
std::string Field;
while(std::getline(inFile, lineread))
{
  std::stringstream sss(lineread);
  std::string Pock_dist;
  std::string Pock_phi;
  std::string Pock_psi;

  std::cout<<lineread<<std::endl;
  std::getline(sss, Pock_dist, '\t'); //read thru tab
  std::cout<<Pock_dist<<std::endl;
  //std::cout<<sss<<std::endl;
  std::getline(sss, Pock_phi, '\t'); //read thru tab
  std::cout<<Pock_phi<<std::endl;
  std::getline(sss, Pock_psi, '\t'); //read thru tab
  std::cout<<Pock_psi<<std::endl;
  //std::stringstream ss(sss);
  //std::getline(ss, Pock_phi, '\t'); //read thru tab
  //std::cout<<Pock_phi<<std::endl;
 // std::stringstream s(ss);
 // std::getline(s, Pock_psi, '\t'); //read thru tab
 // std::cout<<Pock_psi<<std::endl;


/*
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    return split(s, delim, elems);
}

*/
}

}

core::Real Find_Intersect(core::SSize const & phiAngle , core::SSize const & psiAngle, core::Real const & atomX, core::Real const & atomY, core::Real const & atomZ, core::Real const & atom_radius){

        const core::Real PI = numeric::NumericTraits<Real>::pi();
        const core::Real RADS_PER_DEG = PI / 180.;
        const core::Real DEGS_PER_RAD = 180. / PI;
	core::Real const CoMX = 0.0;	core::Real const CoMY = 0.0;	core::Real const CoMZ = 0.0;
	core::Real RandomDist,dirX,dirY,dirZ,dot_direction;

	// compute this from phi/psi, relative to CoM
	RandomDist = rand() % (80-30-1) + 30 + 1;
	dirX = RandomDist*sin(phiAngle*(RADS_PER_DEG))*cos(psiAngle*(RADS_PER_DEG));
	dirY = RandomDist*sin(phiAngle*(RADS_PER_DEG))*sin(psiAngle*(RADS_PER_DEG));
	dirZ = RandomDist*cos(phiAngle*(RADS_PER_DEG));
//	std::cout<<RandomDist<<std::endl;
//	std::cout<<phiAngle<<" " <<psiAngle<<std::endl;

	// figure out whether vector points towards atom or away from atom
	dot_direction = (dirX-CoMX)*(atomX-CoMX) + (dirY-CoMX)*(atomY-CoMY) + (dirZ-CoMZ)*(atomZ-CoMZ);
//	std::cout<<phiAngle<<","<<psiAngle<<","<<RandomDist<<","<<dirX<<","<<dirY<<","<<dirZ<<","<<dot_direction<<std::endl;
	if ( dot_direction < 0.0 ){
//	std::cout <<phiAngle<<" "<<psiAngle<< " Selected direction points away from atom" << std::endl;
	return 0;
	}

	// setup our quadratic equation
	core::Real a = (dirX-CoMX)*(dirX-CoMX) + (dirY-CoMX)*(dirY-CoMY) + (dirZ-CoMZ)*(dirZ-CoMZ);
	core::Real b = 2.0 * ( (dirX - CoMX)*(CoMX - atomX) + (dirY - CoMY)*(CoMY - atomY) + (dirZ - CoMZ)*(CoMZ - atomZ) );
	core::Real c = atomX*atomX + atomY*atomY + atomZ*atomZ + CoMX*CoMX + CoMY*CoMY + CoMZ*CoMZ - 2 * ( atomX * CoMX + atomY * CoMY + atomZ * CoMZ ) - atom_radius * atom_radius;

	// test for intersection
	core::Real inside_sqrt = b * b - 4 * a * c;

	if ( std::abs(inside_sqrt) < 0.00001 ) {
//	std::cout << "Line is tangent to atom\n" << std::endl;
	core::Real mu = -b / ( 2 * a);
	core::Real x = CoMX + ( mu * ( dirX - CoMX) );
	core::Real y = CoMY + ( mu * ( dirY - CoMY) );
	core::Real z = CoMZ + ( mu * ( dirZ - CoMZ) );
	core::Real dist = sqrt( x*x + y*y + z*z );
//	std::cout <<phiAngle<< " " << psiAngle<< " "<< "Tangent at     " << x << " " << y << " " << z << std::endl;
	if (dist > MaxDist[phiAngle][psiAngle]){
	MaxDist[phiAngle][psiAngle] = dist;
	Xpoint[phiAngle][psiAngle] = x;
	Ypoint[phiAngle][psiAngle] = y;
	Zpoint[phiAngle][psiAngle] = z;
	}
	//	std::cout << "Distance from CoM is " << dist << std::endl;
		return 0;
	} else if ( inside_sqrt < 0 ) {
//	std::cout <<phiAngle<<" " <<psiAngle<< " " << "No Intersection" << std::endl;
		return 0;
	} else {
//		std::cout << "Line intersects atom\n" << std::endl;
        core::Real mu1 = -(b-sqrt(inside_sqrt)) / ( 2 * a);
	core::Real mu2 = -(b+sqrt(inside_sqrt)) / ( 2 * a);
	core::Real x1 = CoMX + ( mu1 * ( dirX - CoMX) );
	core::Real y1 = CoMY + ( mu1 * ( dirY - CoMY) );
	core::Real z1 = CoMZ + ( mu1 * ( dirZ - CoMZ) );
 	core::Real dist1 = sqrt( x1*x1 + y1*y1 + z1*z1 );
	core::Real x2 = CoMX + ( mu2 * ( dirX - CoMX) );
	core::Real y2 = CoMY + ( mu2 * ( dirY - CoMY) );
	core::Real z2 = CoMZ + ( mu2 * ( dirZ - CoMZ) );
	core::Real dist2 = sqrt( x2*x2 + y2*y2 + z2*z2 );

//		std::cout << "First distance from CoM is " << dist1 << std::endl;
//		std::cout << "Second intersection is at " << x2 << " " << y2 << " " << z2 << std::endl;
//		std::cout << "Second distance from CoM is " << dist2 << std::endl;
	core::Real max_dist = dist1;
	if ( dist2 > max_dist ) max_dist = dist2;
	if ((max_dist = dist1)){
//	std::cout <<phiAngle<<" "<<psiAngle<< " "<<"Intersection at " << x1 << " " << y1 << " " << z1 <<" Dist from COM "<<dist1<<std::endl;
	if (dist1 > MaxDist[phiAngle][psiAngle]){
	MaxDist[phiAngle][psiAngle] = dist1;
	Xpoint[phiAngle][psiAngle] = x1;
	Ypoint[phiAngle][psiAngle] = y1;
	Zpoint[phiAngle][psiAngle] = z1;
	}
	}
	else if ((max_dist = dist2)){
//	std::cout <<phiAngle<<" "<<psiAngle<< " "<<"Intersection at " << x2 << " " << y2 << " " << z2 <<" Dist from COM "<<dist2<<std::endl;
	if (dist2 > MaxDist[phiAngle][psiAngle]){
	MaxDist[phiAngle][psiAngle] = dist2;
	Xpoint[phiAngle][psiAngle] = x2;
	Ypoint[phiAngle][psiAngle] = y2;
	Zpoint[phiAngle][psiAngle] = z2;
	}
	}
	return 0;
//	std::cout << "std- map " << MaxDist[phiAngle][psiAngle] << std::endl;
//	std::cout << MaxDist[phiAngle][psiAngle] << " " << Xpoint[phiAngle][psiAngle] <<" "<< Ypoint[phiAngle][psiAngle] << " " << Zpoint[phiAngle][psiAngle]<<std::endl;
	}
	return 0;
}

