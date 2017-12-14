// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/packstat/LeeRichards.cc
///
/// @brief
/// @author will sheffler


// Unit header or inline function header
#include <core/scoring/packstat/LeeRichards.hh>
#include <core/pose/Pose.hh>

#include <core/scoring/packstat/types.hh>
#include <core/scoring/packstat/packing_score_params.hh>
#include <utility>
#include <utility/exit.hh>


#include <core/chemical/AtomType.hh>
#include <core/pose/util.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/Fmath.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>
//Auto using namespaces
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
//Auto using namespaces end


namespace core {
namespace scoring {
namespace packstat {

using core::Real;
using utility::vector1;
using core::id::AtomID;

core::Size const N_PROBES = 20;

LeeRichards::LeeRichards(
	PosePackDataOP pd,
	AccumulatorOP accum,
	core::Real spacing,
	core::Real probe_radius,
	bool csa,
	numeric::xyzVector<core::Real> plane
) {
	compute( pd->spheres, accum, spacing, probe_radius, csa, plane );
}

LeeRichards::LeeRichards(
	core::pose::Pose & pose,
	AccumulatorOP accum,
	core::Real spacing,
	core::Real probe_radius,
	bool csa,
	numeric::xyzVector<core::Real> plane
) {
	Spheres spheres;
	for ( int ir = 1; ir <= (int)pose.size(); ++ir ) {
		for ( int ia = 1; ia <= (int)pose.residue(ir).nheavyatoms(); ++ia ) {
			core::id::AtomID aid(ia,ir);
			spheres.push_back(
				Sphere( pose.xyz(aid), pose.residue(ir).atom_type(ia).lj_radius(), aid )
			);
		}
	}
	compute( spheres, accum, spacing, probe_radius, csa, plane );
}

void
LeeRichards::compute(
	Spheres & spheres,
	AccumulatorOP accum,
	core::Real spacing,
	core::Real probe_radius,
	bool csa,
	numeric::xyzVector<core::Real> plane
) {
	// time_t t = clock();

	// get coord system for slices
	plane.normalize();
	numeric::xyzVector<Real> xunit(1,0,0),yunit;
	if ( xunit == plane ) xunit.y(1); // make sure isn't lin. dep.
	xunit -= plane * xunit.dot(plane);
	xunit.normalize();
	yunit = plane.cross(xunit);
	debug_assert( std::abs( 0.0 - xunit.dot(yunit) ) < 1e-9 );
	debug_assert( std::abs( 0.0 - xunit.dot(plane) ) < 1e-9 );
	debug_assert( std::abs( 0.0 - plane.dot(yunit) ) < 1e-9 );
	debug_assert( std::abs( 1.0 - xunit.length()   ) < 1e-9 );
	debug_assert( std::abs( 1.0 - yunit.length()   ) < 1e-9 );
	debug_assert( std::abs( 1.0 - plane.length()   ) < 1e-9 );
	// std::cerr << xunit << std::endl;
	// std::cerr << yunit << std::endl;
	// std::cerr << plane << std::endl;
	// std::exit(-1);

	// compute slice coords
	Real mx=-9e9,mn=9e9;
	for ( auto & sphere : spheres ) {
		Real const d = plane.dot(sphere.xyz);
		if ( d - sphere.radius - probe_radius < mn ) mn = d - sphere.radius-probe_radius;
		if ( d + sphere.radius + probe_radius > mx ) mx = d + sphere.radius+probe_radius;
	}
	auto Nslice = int( (mx-mn) / spacing );
	Real st = mn + ((mx-mn)-Nslice*spacing)/2.0;
	for ( core::Real r = st; r <= mx; r += spacing ) {
		slice_coords_.push_back(r);
	}

	// create Slices
	// std::cerr << "LeeRichards " << mn << " " << mx << " " << spacing << " " << st << std::endl;
	for ( RealCIter i = slice_coords_.begin(); i != slice_coords_.end(); ++i ) {
		// std::cerr << "slice at: " << *i << std::endl;
		bool internal_allowed = true;//*i-4.0 > mn && *i+4 < mx; // not worth trouble
		slices_.push_back( new Slice(accum,spacing,internal_allowed) );
	}

	int count = 0;
	for ( Size i = 1; i <= spheres.size(); ++i ) {
		Sphere & s(spheres[i]);
		if ( s.aid.rsd() == 0 ) s.aid.rsd() = i;
		Real const crd = plane.dot(s.xyz), x = xunit.dot(s.xyz), y = yunit.dot(s.xyz);
		Real const rad = s.radius + probe_radius;
		Size first = find_first_slice(crd-rad);
		if ( 0 == first ) continue;
		Real const stop = crd+rad;
		// std::cerr << "sph " << crd << " " << rad << " " << stop << std::endl;
		while ( first <= slice_coords_.size() && slice_coords_[first] < stop ) {
			// std::cerr << slice_coords_[first] << " ";
			Real d = crd - slice_coords_[first];
			if ( fabs(d) >= rad ) continue;
			Real r = sqrt(rad*rad - d*d);
			Real drdz = d/r;
			Real dada = spacing*rad;
			if ( csa ) {
				drdz = d / r            * (s.radius/(s.radius+probe_radius)); // for CSA not SASA
				dada = spacing*s.radius * (s.radius/(s.radius+probe_radius)); // for CSA not SASA
			}
			// std::cerr << "circ param dsdl " << dada << " "
			// << fmin( 1.0,-d/rad+spacing/2.0)  << " "
			// << fmax(-1.0,-d/rad-spacing/2.0) << std::endl;
			slices_[first]->add_circle( new Circle(x,y,r,drdz,dada,s.aid) );
			++count;
			++first;
		}
		// std::cerr << " " << std::endl;
	}
	// std::cerr << "added " << count << " circles in " << slices_.size() << " slices" << std::endl;

	// for( CircleIter i = slices_[58]->circles_.begin(); i != slices_[58]->circles_.end(); ++i ) {
	//  Circle *c(*i);
	//  std::cerr << "          circles.append( cpp.Circle( 200+20*" << c->x << " , 200+20*" << c->y << " , 20*" << c->r << ") )" << std::endl;
	// }

	for ( Size i = 1; i <= slices_.size(); ++i ) {
		// Real prev = dynamic_cast<AreaAccumulator*>(accum())->total_area;
		// std::cerr << "slice " << i << " " << slices_[i]->circles_.size() << std::endl;
		slices_[i]->compute();
		// std::cerr << "computing area for slice " << i << " " << slice_coords_[i] << " " << dynamic_cast<AreaAccumulator*>(accum())->total_area - prev << std::endl;

	}
	// std::cerr << "slices " << slices_.size() << " area " << dynamic_cast<AreaAccumulator*>(accum())->total_area << " " << clock()-t << std::endl;
}

MultiProbePoseAccumulator::MultiProbePoseAccumulator(
	core::pose::Pose & _pose,
	std::string const & tag
) : pr_idx_(1), pose_(_pose), tag_(tag)
{
	core::pose::initialize_atomid_map_heavy_only( atom_map_, pose_ );
}

void
MultiProbePoseAccumulator::show( std::ostream & out ) {
	Reals tot(N_PROBES,0.0),btot(N_PROBES,0.0);
	for ( Size ir = 1; ir <= pose_.size(); ++ir ) {
		out << "RES_DAT " << pose_.residue(ir).name3() << " " << tag_ << " " << ir << " ";
		for ( Size ia = 1; ia <= pose_.residue(ir).nheavyatoms(); ++ia ) {
			core::id::AtomID aid(ia,ir);
			LR_MP_AtomData & dat(atom_map_[aid]);
			// out << pose_.residue(ir).atom_type_index(ia) << " ";
			for ( Size i = 1; i <= N_PROBES; ++i ) {
				out << dat.area[i] << " " << dat.barea[i] << " ";
				tot[i] += dat.area[i];
				btot[i] += dat.barea[i];
			}
		}
		out << std::endl;
	}
	out << "TOT_DAT " << tag_ << " " << pose_.size() << " ";
	for ( Size i = 1; i <= N_PROBES; ++i ) {
		out << tot[i] << " " << btot[i] << " ";
	}
	out << std::endl;
}


struct HTL_EventX { bool operator()( Event const *a, Event const *b ) {
	if ( a->x==b->x ) return a->y > b->y;
	else             return a->x > b->x;
} };

PointPair
Circle::overlap(
	Circle *other
) {
	Real r0 = r;
	Real r1 = other->r;
	Real x0 = x;
	Real y0 = y;
	Real x1 = other->x;
	Real y1 = other->y;
	Real d = sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1));
	// if( d == 0.0 ) {
	//  return; // don't error anymode
	//  std::cerr << "overlap(_a,other), _a/other same point" << std::endl;
	//  std::exit(-1);
	// }
	if ( d >= r0+r1 || d <= fabs(r0-r1) ) return PointPair(0.0,0.0,0.0,0.0);
	Real _a = (r0*r0 - r1*r1 + d*d ) / (2*d);
	if ( r0*r0 < _a*_a ) {
		std::cerr << "r0*r0 < _a*_a: \n" << x << " " << y << " " << r << " \n" << other->x << " " << other->y << " " << other->r << " " << d << std::endl;
		std::cerr << "abs(r0-r1) " << fabs(r0-r1) << std::endl;
		utility_exit_with_message( "Circle::overlap error" );
	}
	Real h = sqrt(r0*r0 - _a*_a);
	Real x2 = x0 + _a * ( x1 - x0 ) / d;
	Real y2 = y0 + _a * ( y1 - y0 ) / d;
	Real x3 = x2 + h * ( y1 - y0 ) / d;
	Real y3 = y2 - h * ( x1 - x0 ) / d;
	Real x4 = x2 - h * ( y1 - y0 ) / d;
	Real y4 = y2 + h * ( x1 - x0 ) / d;
	// # if x3 > x4: x3,y3,x4,y4 = x4,y4,x3,y3
	return PointPair(x3,y3,x4,y4);
}

std::ostream &
operator<< ( std::ostream & out, Circle & circle ) {
	out << "Circle( " << circle.x << " " << circle.y << " " << circle.r << " )";
	return out;
}

Slice::~Slice() {
	for ( auto & event : events_ ) delete event;
	for ( auto & circle : circles_ ) delete circle;
	for ( auto & trace : traces_ ) delete trace;
}

void
Slice::compute_events()
{
	Octree2D nbr(circles_);

	for ( auto circle : circles_ ) {
		if ( !nbr.contains(circle->x+circle->r,circle->y,2,circle) ) {
			events_.push_back( new Event( circle->x + circle->r, circle->y, ENTER, circle, nullptr ) );
		}
		if ( !nbr.contains(circle->x-circle->r,circle->y,2,circle) ) {
			events_.push_back( new Event( circle->x - circle->r, circle->y, EXIT , circle, nullptr ) );
		}

		int I = nbr.get_i(circle->x);
		int J = nbr.get_j(circle->y);
		for ( int i = std::max(I-2,0); i <= std::min(I+2,(int)nbr.Xdim_); ++i ) {
			for ( int j = std::max(J-2,0); j <= std::min(J+2,(int)nbr.Ydim_); ++j ) {
				for ( auto ja = nbr.cubes(i,j).begin(); ja != nbr.cubes(i,j).end(); ja++ ) {
					Circle *circle2(*ja);
					if ( circle <= circle2 ) continue;
					if ( ! circle->does_overlap(circle2) ) continue;

					PointPair pp = circle->overlap(circle2);
					if ( 0.0==pp.a.x && 0.0==pp.a.y && 0.0==pp.b.x && 0.0==pp.b.y ) continue;

					if ( !nbr.contains(pp.a.x,pp.a.y,1,circle,circle2) ) {
						events_.push_back( new Event(pp.a.x,pp.a.y,ISECT,circle,circle2) );
					}
					if ( !nbr.contains(pp.b.x,pp.b.y,1,circle,circle2) ) {
						events_.push_back( new Event(pp.b.x,pp.b.y,ISECT,circle2,circle) );
					}

				}
			}
		}
	}
	std::sort( events_.begin(), events_.end(), HTL_EventX() );
	// std::cerr << "DONE " << events_.size() << std::endl;
}


void
Slice::compute_surface()
{
	for ( auto e : events_ ) {
		// std::cerr << "EVENT " << e << " " << e->x << " " << e->y << " " << e->kind << std::endl;
		Circle *circle = e->circle;
		if     ( circle->tcw  ) e->trace_ = circle->tcw;
		else if ( circle->tccw ) e->trace_ = circle->tccw;
		if ( e->kind == ENTER ) {
			traces_.push_back( new trace( accum_, e, false, circle, 0.0, nullptr           ) );
			traces_.push_back( new trace( accum_, e, true , circle, 0.0, traces_.back() ) );
		} else if ( e->kind == EXIT ) {
			circle->tcw ->end( circle->tccw );
			circle->tccw->end( circle->tcw  );
		} else if ( e->kind == ISECT ) {
			Circle *ccwcircle = e->ccw;
			if     ( ccwcircle->tcw  ) e->trace_ = ccwcircle->tcw;
			else if ( ccwcircle->tccw ) e->trace_ = ccwcircle->tccw;
			bool have_cw_trace  = (    circle->tcw  != nullptr && e->y > circle->y    );
			bool have_ccw_trace = ( ccwcircle->tccw != nullptr && e->y < ccwcircle->y );
			if ( have_cw_trace && have_ccw_trace ) {
				circle   ->tcw ->end( ccwcircle->tccw, e->cw_angle  );
				ccwcircle->tccw->end(    circle->tcw , e->ccw_angle );
			} else if ( !have_cw_trace && !have_ccw_trace ) {
				traces_.push_back( new trace( accum_, e, true , circle, e->cw_angle , nullptr           ) );
				traces_.push_back( new trace( accum_, e, false, e->ccw, e->ccw_angle, traces_.back() ) );
			} else if ( have_cw_trace ) {
				circle->tcw->next_circle(ccwcircle,e->cw_angle,e->ccw_angle);
			} else if ( have_ccw_trace ) {
				ccwcircle->tccw->next_circle(circle,e->ccw_angle,e->cw_angle);
			} else {
				utility_exit_with_message( "Slice::compute_surface()" );
			}
		}
	}

	for ( auto t : traces_ ) {
		if ( t->next_trace_ ) {
			t->start_ = t->next_trace_->get_first(t);
			t->next_trace_->set_first(t,t->start_); // sets next_trace_ to NULL for all
		}
		bool is_internal = ( t->start_->kind == ISECT ) && internal_allowed_;
		// for debugging
		// if( is_internal ) {
		//  for( CircleIter ic = circles_.begin(); ic != circles_.end(); ++ic ) {
		//   Circle *c(*ic);
		//   std::cerr << "    circles.append( cpp.Circle( " << c->x << "," << c->y << "," << c->r << "));" << std::endl;
		//  }
		//  std::exit(-1);
		// }
		for ( auto & arc : t->arcs_ ) {
			accum_->accumulate_area(arc.first,arc.second,is_internal);
		}
	}

	// std::sort( arcs_.begin(), arcs_.end(), OrderArc() ); // TODO: take out for speed
}

void
Slice::compute_derivatives()
{
	// std::cerr << "compute_derivatives" << std::endl;
	for ( auto e : events_ ) {
		if ( e->kind != ISECT ) continue;
		bool is_internal = ( e->trace_->start_->kind == ISECT ) && internal_allowed_;
		Real mg = 1.0 / sin( ( e->cw_angle - e->ccw_angle + numeric::constants::d::pi ) / 2.0 );
		mg = min(mg,100.0);
		mg = mg * (( e->circle->dada / e->circle->r ) + ( e->ccw->dada / e->ccw->r ))/2.0;
		// Real mg1 = mg * ( e->circle->dada / e->circle->r );
		// Real mg2 = mg * ( e->ccw->dada / e->ccw->r );
		// std::cerr << "mg1/mg2 " << mg1 << " " << mg2 << std::endl;
		//debug_assert( fabs(mg) < 100.0 );
		Real dir = (e->cw_angle + e->ccw_angle) / 2.0 + numeric::constants::d::pi/2.0;
		accum_->accumulate_dxdy( e->circle->atom,  mg*cos(dir), -mg*sin(dir), is_internal );
		accum_->accumulate_dxdy( e->ccw   ->atom, -mg*cos(dir),  mg*sin(dir), is_internal );
	}
}


PackingScoreResDataOP
MultiProbePerSphereAccumulator::compute_surrounding_sasa( XYZ & xyz )
{
	using namespace utility;
	using namespace numeric;

	size_t Nprobes = 31;
	size_t Nshells = 7;
	PackingScoreResDataOP psrdOP( new PackingScoreResData(Nshells,Nprobes) );
	for ( size_t is = 1; is <= pd_->spheres.size(); ++is ) {
		Sphere const & sphere( pd_->spheres[is] );
		//if( sphere.xyz.x() > xyz.x() + dist_th ) break;
		PackstatReal dist = xyz.distance( sphere.xyz );
		for ( size_t id = 1; id <= Nshells; ++id ) {
			if ( dist <= (PackstatReal)id ) {
				for ( size_t pr = 1; pr <= Nprobes; ++pr ) {
					psrdOP->msa(id,pr) += atom_map_[AtomID(1,is)].area[pr];
				}
				break;
			}
		}
	}

	PackingScoreResDataOP rtn( new PackingScoreResData(Nshells,Nprobes-1) );
	for ( size_t id = 1; id <= Nshells; ++id ) {
		for ( size_t pr = 2; pr <= Nprobes; ++pr ) {
			rtn->msa(id,pr-1) = psrdOP->msa(id,pr) - psrdOP->msa(id,1); // stupid reverse indicies...
		}
	}

	return rtn;

}


Real
compute_surface_area_leerichards(
	Real & buried_area_out,
	PosePackDataOP pd,
	Real slicesize,
	Real pr,
	bool csa,
	numeric::xyzVector<Real> plane
) {
	if ( plane.length() > 0 ) {
		AreaAccumulatorOP accum( new AreaAccumulator );
		AccumulatorOP accumOP(accum);
		LeeRichards lr( pd, accumOP, slicesize, pr, csa, plane );
		buried_area_out = accum->buried_area;
		return accum->total_area;
	} else {
		Real area = 0.0;
		{
			AreaAccumulatorOP accum( new AreaAccumulator );
			AccumulatorOP accumOP(accum);
			LeeRichards lr( pd, accumOP, slicesize, pr, csa, XYZ(1,0,0) );
			buried_area_out += accum->buried_area;
			area += accum->total_area / 3.0;
		}
		{
			AreaAccumulatorOP accum( new AreaAccumulator );
			AccumulatorOP accumOP(accum);
			LeeRichards lr( pd, accumOP, slicesize, pr, csa, XYZ(0,1,0) );
			buried_area_out += accum->buried_area;
			area += accum->total_area / 3.0;
		}
		{
			AreaAccumulatorOP accum( new AreaAccumulator );
			AccumulatorOP accumOP(accum);
			LeeRichards lr( pd, accumOP, slicesize, pr, csa, XYZ(0,0,1) );
			buried_area_out += accum->buried_area;
			area += accum->total_area / 3.0;
		}
		return area;
	}
}

Real
compute_surface_area_leerichards(
	PosePackDataOP pd,
	Real slicesize,
	Real pr,
	bool csa,
	numeric::xyzVector<Real> plane
) {
	Real dummy = 0.0;
	return compute_surface_area_leerichards(dummy,pd,slicesize,pr,csa,plane);
}

XYZs
compute_surface_area_leerichards_deriv(
	PosePackDataOP pd,
	Real slicesize,
	Real pr,
	bool csa//,
	// numeric::xyzVector<Real> plane
) {
	XYZs rtn(pd->spheres.size(),XYZ(0.0,0.0,0.0));
	{
		PerSphereAccumulatorOP accum( new PerSphereAccumulator(pd->spheres) );
		AccumulatorOP accumOP(accum);
		LeeRichards lr( pd, accumOP, slicesize, pr, csa, XYZ(1,0,0) );
		for ( Size i = 1; i <= pd->spheres.size(); ++i ) {
			rtn[i].y() += accum->atom_map_[AtomID(1,i)].dx / 2.0;
			rtn[i].z() += accum->atom_map_[AtomID(1,i)].dy / 2.0;
		}
	}
	{
		PerSphereAccumulatorOP accum( new PerSphereAccumulator(pd->spheres) );
		AccumulatorOP accumOP(accum);
		LeeRichards lr( pd, accumOP, slicesize, pr, csa, XYZ(0,1,0) );
		for ( Size i = 1; i <= pd->spheres.size(); ++i ) {
			rtn[i].x() += accum->atom_map_[AtomID(1,i)].dx / 2.0;
			rtn[i].z() -= accum->atom_map_[AtomID(1,i)].dy / 2.0;
		}
	}
	{
		PerSphereAccumulatorOP accum( new PerSphereAccumulator(pd->spheres) );
		AccumulatorOP accumOP(accum);
		LeeRichards lr( pd, accumOP, slicesize, pr, csa, XYZ(0,0,1) );
		for ( Size i = 1; i <= pd->spheres.size(); ++i ) {
			rtn[i].x() += accum->atom_map_[AtomID(1,i)].dx / 2.0;
			rtn[i].y() += accum->atom_map_[AtomID(1,i)].dy / 2.0;
		}
	}
	return rtn;
}

XYZs
check_surface_area_leerichards_deriv(
	PosePackDataOP pd,
	Real slicesize,
	Real pr,
	bool csa,
	Size max_num
	// numeric::xyzVector<Real> plane
) {

	Real D = 0.000001;

	Real new_area1, new_area2;
	Real orig_area_x = compute_surface_area_leerichards(pd,slicesize,pr,csa,XYZ(1,0,0));
	Real orig_area_y = compute_surface_area_leerichards(pd,slicesize,pr,csa,XYZ(0,1,0));
	Real orig_area_z = compute_surface_area_leerichards(pd,slicesize,pr,csa,XYZ(0,0,1));

	Size N = min((Size)max_num,pd->spheres.size());
	XYZs rtn(N);

	for ( Size i = 1; i <= N; ++i ) {
		std::cerr << "check_surface_area_leerichards_deriv " << i << std::endl;

		pd->spheres[i].xyz.x() += D;
		new_area1 = compute_surface_area_leerichards(pd,slicesize,pr,csa,XYZ(0,1,0));
		new_area2 = compute_surface_area_leerichards(pd,slicesize,pr,csa,XYZ(0,0,1));
		rtn[i].x( ((new_area1+new_area2)/2.0-(orig_area_y+orig_area_z)/2.0) / D );
		pd->spheres[i].xyz.x() -= D;

		pd->spheres[i].xyz.y() += D;
		new_area1 = compute_surface_area_leerichards(pd,slicesize,pr,csa,XYZ(1,0,0));
		new_area2 = compute_surface_area_leerichards(pd,slicesize,pr,csa,XYZ(0,0,1));
		rtn[i].y( ((new_area1+new_area2)/2.0-(orig_area_x+orig_area_z)/2.0) / D );
		pd->spheres[i].xyz.y() -= D;

		pd->spheres[i].xyz.z() += D;
		new_area1 = compute_surface_area_leerichards(pd,slicesize,pr,csa,XYZ(1,0,0));
		new_area2 = compute_surface_area_leerichards(pd,slicesize,pr,csa,XYZ(0,1,0));
		rtn[i].z( ((new_area1+new_area2)/2.0-(orig_area_x+orig_area_y)/2.0) / D );
		pd->spheres[i].xyz.z() -= D;
	}
	return rtn;
}


Real
compute_packing_score_leerichards(
	PosePackDataOP pd,
	Real slicesize,
	numeric::xyzVector<Real> plane
) {
	MultiProbePerSphereAccumulatorOP accum( new MultiProbePerSphereAccumulator(pd) );
	AccumulatorOP accumOP = accum;
	for ( Size pri = 1; pri <= 31; pri++ ) {
		Real pr = Real(31-pri)/10.0;
		// std::cerr << "pr: " << pr << " ";
		accum->set_pr_idx(pri);
		LeeRichards lr( pd, accumOP, slicesize, pr, true, plane );
	}
	vector1< PackingScoreResDataCOP > psrds;
	for ( Size i = 1; i <= pd->centers.size(); ++i ) {
		PackingScoreResDataOP psrd = accum->compute_surrounding_sasa(pd->centers[i]);
		// std::cerr << "LeeRichardsRAW " << i << " " << pd->labels[i] << " " << *psrd << std::endl;
		psrds.push_back(psrd);
	}
	PackingScore ps_discrim(7,30,true);
	init_packing_score_discrim( ps_discrim );
	// accum->show(std::cerr);
	return ps_discrim.score( psrds );
}

} // namespace packstat
} // namespace scoring
} // namespace core
