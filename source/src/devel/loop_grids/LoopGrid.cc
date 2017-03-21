// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/LoopGrid.cc
/// @brief  generate accessible grids for loop
/// @author Yuan Liu (wendao@u.washington.edu)

//#include <protocols/match/LoopGrid.hh>
#include <devel/loop_grids/LoopGrid.hh>

#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

#include <iostream>
#include <fstream>
#include <sstream>

#include <protocols/match/BumpGrid.hh>
#include <utility/vector1.hh>


using namespace std;
using namespace core;
using namespace core::pose;

static THREAD_LOCAL basic::Tracer TR( "protocols.match.LoopGrid" );

namespace protocols {
namespace match {

LoopGrid::LoopGrid( Pose const &pose, Loop const &loop ) :
	grids_( loop.size()-2 ), //not for the two end residues
	bin_width_( 1.0 ),
	pose_(pose),
	loop_(loop)
{
	debug_assert(!loop.is_terminal(pose));//not a tail

	Vector center1 = pose.residue(loop.start()-1).atom("CA").xyz();
	Vector center2 = pose.residue(loop.stop()+1).atom("CA").xyz();

	Size ngrids = grids_.size();
	for ( Size i=1; i<=ngrids; i++ ) {
		//for each grid
		Size n1 = i+1;          //distance to the left root
		Size n2 = ngrids+2-i;   //distance to the right root
		Real r1 = loop_length_cutoff[n1];
		Real r2 = loop_length_cutoff[n2];
		TR << "Grid " << i << " : init ...  " << "r1=" << r1 << ", r2=" << r2 << endl;
		grids_[i] = protocols::match::Bool3DGridOP( new Bool3DGrid );
		create_grid_and_by_two_sphere(grids_[i], center1, r1 , center2, r2);

		//output
		//ostringstream fn;
		//fn << "grid_loop_" << i << ".gridlig";
		//write_grids_file(fn.str().c_str(), i);
	}
}

LoopGrid::~LoopGrid()
{
	TR << "Deleting LoopGrid ..." << endl;
}

/// @brief setup the grid in the overlap of two spheres
void LoopGrid::create_grid_and_by_two_sphere( Bool3DGridOP &gridop,
	Vector center1, Real r1, Vector center2, Real r2 )
{
	static const Real loop_ca_pro_cutoff = 6.32;
	static const Real loop_ca_pro_cutoff2 = loop_ca_pro_cutoff*loop_ca_pro_cutoff;
	static const Real too_far_away_cutoff = 18.0;
	static const Real too_far_away_cutoff2 = too_far_away_cutoff*too_far_away_cutoff;
	static const Size n_in_bump_cutoff = 8;

	//overlap of the cubes that contain the spheres
	Vector overlap_low(
		std::max( center1.x()-r1, center2.x()-r2 ),
		std::max( center1.y()-r1, center2.y()-r2 ),
		std::max( center1.z()-r1, center2.z()-r2 ));

	Vector overlap_high(
		std::min( center1.x()+r1, center2.x()+r2 ),
		std::min( center1.y()+r1, center2.y()+r2 ),
		std::min( center1.z()+r1, center2.z()+r2 ));

	//test
	TR << "region old" << endl;
	TR << "low: " << overlap_low.x() <<" "<< overlap_low.y() <<" "<< overlap_low.z() << endl;
	TR << "high: " << overlap_high.x() <<" "<< overlap_high.y() <<" "<< overlap_high.z() << endl;

	BoundingBox bb(overlap_low, overlap_high);
	Bool3DGrid grid;
	grid.set_bin_width(bin_width_);
	grid.set_bounding_box(bb);
	BoundingBox newbb = grid.actual_bb();
	Real new_base_x = newbb.lower().x();
	Real new_base_y = newbb.lower().y();
	Real new_base_z = newbb.lower().z();

	Vector real_bb_low(9999,9999,9999);
	Vector real_bb_high(-9999,-9999,-9999);

	Real CA_radius = BumpGrid::probe_radius(C_ALA);
	Bin3D dimsizes=grid.dimsizes();
	Bin3D bin;
	for ( Size i=0; i<dimsizes[1]; i++ ) {
		bin[1] = i;
		for ( Size j=0; j<dimsizes[2]; j++ ) {
			bin[2] = j;
			for ( Size k=0; k<dimsizes[3]; k++ ) {
				bin[3] = k;

				//all corners
				CornerPoints grid_corners = grid.corners( bin );
				//center
				Vector bin_center_point( grid.bin_center(bin) );

				//intersection
				Real d1_squared = center1.distance_squared( bin_center_point );
				Real d2_squared = center2.distance_squared( bin_center_point );
				if ( d1_squared < r1*r1 && d2_squared < r2*r2 ) {
					bool not_far(false);//not very far away
					bool not_bump(true);//not in bump region

					for ( Size n=1,nr=pose_.size(); n<=nr && not_bump; n++ ) {
						if ( n==loop_.start() ) n+=loop_.size();
						Vector center = pose_.residue(n).atom("CA").xyz();
						//the reside is too far away from this bin
						if ( center.distance_squared(bin_center_point) > too_far_away_cutoff2 ) continue;

						//go through this residue
						core::conformation::ResidueCOP res( pose_.residue(n).get_self_ptr() );
						Size na = res->atoms().size();
						for ( Size ia=1; ia<=na && not_bump; ia++ ) {
							//for each atom
							//we need the r and xyz
							//place a Calpha at the center of the bin
							//not all corners touch atom
							//and can not too far away from them

							//center distance
							Real d_squared = res->xyz(ia).distance_squared(bin_center_point);
							if ( d_squared < loop_ca_pro_cutoff2 ) {
								//not too far
								not_far = true;
							}

							//get radius of probe
							Real r_probe = BumpGrid::probe_radius(probe_radius_for_atom_type(res->atom(ia).type()));
							Real r_probe_squared = (r_probe+CA_radius)*(r_probe+CA_radius);
							Size n_in=0;
							for ( Size ic=1; ic<=8; ic++ ) {
								//for each corner
								Real dc_squared = res->xyz(ia).distance_squared(grid_corners[ic]);
								if ( dc_squared < r_probe_squared ) {
									//clash
									n_in++;
								}
							}
							if ( n_in>=n_in_bump_cutoff ) {
								not_bump = false;
								break;
							}
						}
					}
					if ( not_far && not_bump ) {
						grid.set_value_for_bin( bin, true );
						//save the boundary
						Real bin_left = new_base_x + i * bin_width_;
						Real bin_right = new_base_x + (i+1) * bin_width_;
						Real bin_front = new_base_y + j * bin_width_;
						Real bin_back = new_base_y + (j+1) * bin_width_;
						Real bin_bottom = new_base_z + k * bin_width_;
						Real bin_top = new_base_z + (k+1) * bin_width_;
						if ( real_bb_low.x()>bin_left ) real_bb_low.x(bin_left);
						if ( real_bb_low.y()>bin_front ) real_bb_low.y(bin_front);
						if ( real_bb_low.z()>bin_bottom ) real_bb_low.z(bin_bottom);
						if ( real_bb_high.x()<bin_right ) real_bb_high.x(bin_right);
						if ( real_bb_high.y()<bin_back ) real_bb_high.y(bin_back);
						if ( real_bb_high.z()<bin_top ) real_bb_high.z(bin_top);
					}
				}
			}//k
		}//j
	}//i

	//test
	TR << "region new" << endl;
	TR << "low: " << real_bb_low.x() <<" "<< real_bb_low.y() <<" "<< real_bb_low.z() << endl;
	TR << "high: " << real_bb_high.x() <<" "<< real_bb_high.y() <<" "<< real_bb_high.z() << endl << endl;

	//reshape and copy to gridop_real
	gridop->set_bin_width(bin_width_);
	BoundingBox real_bb(real_bb_low,real_bb_high);
	gridop->set_bounding_box(real_bb);
	gridop->or_with(grid);

	//out some info
	Bin3D dimsizes_real=gridop->dimsizes();
	TR << dimsizes_real[1] << "," << dimsizes_real[2] << "," << dimsizes_real[3] << endl;
}

/// @brief output the ndxth grid to file fn
void LoopGrid::write_grids_file(const char *fn, Size ndx)
{
	Bin3D dimsizes=grids_[ndx]->dimsizes();
	BoundingBox bb=grids_[ndx]->actual_bb();
	Vector lower=bb.lower();
	Vector upper=bb.upper();
	Real width = grids_[ndx]->bin_width();

	//output the grid file
	std::ofstream out(fn);
	out << "NAME: gridlig" << endl;
	out << "BASE: " << lower.x() << " " << lower.y() << " " << lower.z() << endl;
	out << "SIZE: " << dimsizes[1] << " "<< dimsizes[2] << " "<< dimsizes[3] << endl;
	out << "LENGTH: " << width << " " << width << " " << width << endl;

	Bin3D bin;
	for ( Size i=0; i<dimsizes[1]; i++ ) {
		bin[1] = i;
		for ( Size j=0; j<dimsizes[2]; j++ ) {
			bin[2] = j;
			for ( Size k=0; k<dimsizes[3]; k++ ) {
				bin[3] = k;
				if ( grids_[ndx]->occupied(bin) ) out<<"1";
				else out<<"0";
				out << " ";
			}
			out << endl;
		}
		out << endl;
	}
}

bool LoopGrid::occupied( Vector const & xyz ) const
{
	//test in every grid
	for ( Size i=1, ng=grids_.size(); i<=ng; i++ ) {
		if ( grids_[i]->occupied(xyz) ) return true;
	}
	return false;
}

bool LoopGrid::occupied( Bin3D const & bin ) const
{
	for ( Size i=1, ng=grids_.size(); i<=ng; i++ ) {
		if ( grids_[i]->occupied(bin) ) return true;
	}
	return false;
}

bool LoopGrid::occupied( Pose const &pose ) const
{
	//makesure the same backbone
	debug_assert( pose.size() == pose_.size() );
	bool in_grid(true);
	for ( Size i=1, ng=grids_.size(); i<=ng; i++ ) {
		Size nr = loop_.start() + i;
		if ( !grids_[i]->occupied(pose.residue(nr).atom("CA").xyz()) ) {
			TR << "LoopGridTest: my fault: nres=" << nr <<" ngrid="<< i << endl;
			in_grid = false;
			break;
		}
	}
	return in_grid;
}

bool LoopGrid::occupied( Pose const &pose, Size start ) const
{
	//the pose is only a part of the whole protein (loop)
	//start in the first residue's number in the old one
	debug_assert( start <= loop_.start()+1 );
	debug_assert( pose.size() >= loop_.stop()-start );

	bool in_grid(true);
	for ( Size i=1, ng=grids_.size(); i<=ng; i++ ) {
		Size nr = i + loop_.start() + 1 - start ;
		if ( !grids_[i]->occupied(pose.residue(nr).atom("CA").xyz()) ) {
			in_grid = false;
			break;
		}
	}
	return in_grid;
}

bool LoopGrid::occupied( Pose const &pose, Size start, Size stop ) const
{
	//the pose is only a part of the loop
	//start and stop is the old res number
	debug_assert( start >= loop_.start()+1 );
	debug_assert( stop <= loop_.stop()-1 );
	debug_assert( pose.size() == stop-start+1 );

	bool in_grid(true);
	for ( Size n=start; n<=stop; n++ ) {
		Size i = n - loop_.start();
		Size nr = n - start + 1;
		if ( !grids_[i]->occupied(pose.residue(nr).atom("CA").xyz()) ) {
			in_grid = false;
			break;
		}
	}
	return in_grid;
}

bool LoopGrid::occupied( Residue const &res, Size nres ) const
{
	debug_assert( nres >= loop_.start()+1 );
	debug_assert( nres <= loop_.stop()-1 );
	return grids_[nres-loop_.start()]->occupied(res.atom("CA").xyz());
}

/*overlap_boundary_2d_spheres
/// @brief just take the first two dimensions, and the first two corners
LoopGrid::CornerPoints LoopGrid::overlap_boundary_2d_spheres( xyVector center1, Real r1, xyVector center2, Real r2 )
{
//test
assert(r1>0);
assert(r2>0);

//calculate the distance of the two centers
Real dx = center2[1]-center1[1];
Real dy = center2[2]-center1[2];
Real d2 = dx*dx + dy*dy;
Real d = sqrt(d2);

//return value
CornerPoints cp2d;
//d
if (d+r1<=r2)
{
//using 1
cp2d[1] = Vector( center1[1] - r1, center1[2] - r1, 0.0 );
cp2d[2] = Vector( center1[1] + r1, center1[2] + r1, 0.0 );
}
else if (d+r2<=r1)
{
//using 2
cp2d[1] = Vector( center2[1] - r2, center2[2] - r2, 0.0 );
cp2d[2] = Vector( center2[1] + r2, center2[2] + r2, 0.0 );
}
else if (d>=r1+r2)
{
//noting
cp2d[1] = Vector( 0.0, 0.0, 0.0 );
cp2d[2] = Vector( 0.0, 0.0, 0.0 );
}
else
{
//calculate the intersection points
Real sine = dy / d;     //the local axis
Real cosine = dx / d;   //the local axis
Real cosine_theta = (r1*r1 + d2 - r2*r2)/(2.0*r1*d);
Real sine_theta = sqrt(1-cosine_theta*cosine_theta);

//new position in local axis
Real newx1 = r1 * cosine_theta;
Real newx2 = r1 * cosine_theta;
Real newy1 = r1 * sine_theta;
Real newy2 = - r1 * sine_theta;

//rotate, back a
//x' = xcosa - ysina
//y' = xsina + ycosa
Real originx1 = newx1*cosine - newy1*sine + center1[1];
Real originx2 = newx2*cosine - newy2*sine + center1[1];
Real originy1 = newx1*sine + newy1*cosine + center1[2];
Real originy2 = newx2*sine + newy2*cosine + center1[2];

//TR << newx1 << endl;
//TR << newy1 << endl;
//TR << newx2 << endl;
//TR << newy2 << endl;

//TR << endl;
//TR << sine << "  " << cosine << endl;
//TR << endl;

//TR << originx1 << endl;
//TR << originy1 << endl;
//TR << originx2 << endl;
//TR << originy2 << endl;

//boundary
Vector left;  //std::max( center1[1]-r1, center2[1]-r2 );
Vector right; //std::min( center1[1]+r1, center2[1]+r2 );
Vector top;   //std::min( center1[2]+r1, center2[2]+r2 );
Vector bottom;//std::max( center1[2]-r1, center2[2]-r2 );

//left
if ((center1[1]-r1) > (center2[1]-r2)) left = Vector( center1[1]-r1, center1[2], 0.0 );
else left = Vector( center2[1]-r2, center2[2], 0.0 );
//right
if ((center1[1]+r1) < (center2[1]+r2)) right = Vector( center1[1]+r1, center1[2], 0.0);
else right = Vector( center2[1]+r2, center2[2], 0.0 );
//top
if ((center1[2]+r1) < ( center2[2]+r2)) top = Vector( center1[1], center1[2]+r1, 0.0 );
else top = Vector( center2[1], center2[2]+r2, 0.0 );
//bottom
if ((center1[2]-r1) > (center2[2]-r2)) bottom = Vector( center1[1], center1[2]-r1, 0.0 );
else bottom = Vector( center2[1], center2[2]-r2, 0.0 );

assert(left.x() < right.x());
assert(top.y() > bottom.y());

Real overlap_left = std::min( originx1, originx2 );
Real overlap_right = std::max( originx1, originx2 );
Real overlap_top = std::max( originy1, originy2 );
Real overlap_bottom = std::min( originy1, originy2 );

Real real_left, real_right, real_top, real_bottom;

//real region
if (overlap_bottom<left.y() && overlap_top>left.y()) real_left = left.x();
else real_left = overlap_left;
if (overlap_bottom<right.y() && overlap_top>right.y()) real_right = right.x();
else real_right = overlap_right;
if (overlap_left<bottom.x() && overlap_right>bottom.x()) real_bottom = bottom.y();
else real_bottom = overlap_bottom;
if (overlap_left<top.x() && overlap_right>top.x()) real_top = top.y();
else real_top = overlap_top;

cp2d[1] = Vector(real_left, real_bottom, 0.0);
cp2d[2] = Vector(real_right, real_top, 0.0);
}

return cp2d;
}
*/

const core::Real LoopGrid::loop_length_cutoff[] =
{
0.0,
3.8,
7.2,
10.6,
13.9,
17.1,
20.2,
23.2,
26.0,
28.7,
31.1,
33.2,
35.1,
36.7,
38.1,
39.3,
40.4,
41.5,
42.4,
43.1,
43.8,
44.5,
45.1,
45.6,
46.0,
46.4,
46.8,
47.3,
47.7,
48.1,
48.5,
48.9,
49.4,
49.9,
50.4,
50.8,
51.1,
51.4,
51.7,
51.9,
52.2,
52.4,
52.7,
52.9,
52.9,
53.0,
53.3,
53.4,
53.6,
53.8,
54.0,
54.2,
54.4,
54.6,
54.9,
55.0,
55.2,
55.4,
55.6,
55.8,
56.0
};

}//namespace protocols
}//namespace match

