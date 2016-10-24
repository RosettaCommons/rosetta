// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/packstat/PackingScore.hh
///
/// @brief
/// @author will sheffler


#ifndef INCLUDED_core_scoring_packstat_LeeRichards_hh
#define INCLUDED_core_scoring_packstat_LeeRichards_hh

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>

#include <core/scoring/packstat/types.hh>
#include <core/scoring/packstat/PackingScore.hh>

#include <core/pose/Pose.hh>

#include <numeric/xyzVector.hh>

#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>

#include <numeric/constants.hh>

#include <iostream>

#include <core/kinematics/Jump.hh>
#include <utility/vector1.hh>
#include <utility/numbers.hh>

namespace core {
namespace scoring {
namespace packstat {

extern core::Size const N_PROBES;

struct Event;
struct Circle;
struct trace;
struct Slice;
typedef utility::vector1< Event* >  Events;
typedef utility::vector1< Slice* >  Slices;
typedef utility::vector1< trace* >  traces;
typedef utility::vector1< Circle* > Circles;
typedef utility::vector1< Event* >::iterator  EventIter;
typedef utility::vector1< Circle* >::iterator CircleIter;
typedef utility::vector1< Slice* >::iterator  SliceIter;
typedef utility::vector1< trace* >::iterator  traceIter;

typedef std::pair< core::id::AtomID,Real > Arc;
typedef utility::vector1< Arc > Arcs;
typedef utility::vector1< Arc >::iterator ArcIter;

struct Accumulator : public utility::pointer::ReferenceCount {
	virtual void accumulate_area( core::id::AtomID atom, core::Real area, bool buried ) = 0;
	virtual void accumulate_dxdy( core::id::AtomID atom, core::Real dx, core::Real dy, bool buried ) = 0;
	//dz// virtual void accumulate_dz  ( core::id::AtomID atom, core::Real dz ) = 0;
};
typedef utility::pointer::shared_ptr<Accumulator> AccumulatorOP;

struct AreaAccumulator : public Accumulator {
	AreaAccumulator() : total_area(0.0), buried_area(0.0) {}
	virtual void accumulate_area( core::id::AtomID, core::Real area, bool buried ) {
		// std::cerr << "accumulate_area " << area << std::endl;
		debug_assert( !utility::isnan(area) );
		total_area += area;
		if ( buried ) buried_area += area;
	}
	void accumulate_dxdy( core::id::AtomID, core::Real, core::Real, bool ) {}
	//dz// void accumulate_dz( core::id::AtomID, core::Real ) {}
	void reset() { total_area = 0.0; }
	core::Real total_area, buried_area;
};
typedef utility::pointer::shared_ptr<AreaAccumulator> AreaAccumulatorOP;

struct LR_AtomData {
	LR_AtomData() : area(0.0),dx(0.0),dy(0.0)/*,dz(0.0)*/ {}
	Real area,dx,dy/*,dz*/;
};

struct PerSphereAccumulator : public Accumulator {
	PerSphereAccumulator( Spheres & spheres ) {
		atom_map_.resize(spheres.size());
		for ( Size i = 1; i <= spheres.size(); ++i ) atom_map_.resize(i,1);
	}
	virtual void accumulate_area( core::id::AtomID id, core::Real area, bool ) {
		atom_map_[id].area += area;
	}
	virtual void accumulate_dxdy( core::id::AtomID id, core::Real dx, core::Real dy, bool ) {
		atom_map_[id].dx += dx;
		atom_map_[id].dy += dy;
	}
	//dz// virtual void accumulate_dz( core::id::AtomID id, core::Real dz ) {
	//dz//  atom_map_[id].dz += dz;
	//dz// }
	core::id::AtomID_Map<LR_AtomData> atom_map_;
};
typedef utility::pointer::shared_ptr<PerSphereAccumulator> PerSphereAccumulatorOP;

struct LR_MP_AtomData {
	LR_MP_AtomData() : area(N_PROBES,0.0),dx(N_PROBES,0.0),dy(N_PROBES,0.0),barea(N_PROBES,0.0),bdx(N_PROBES,0.0),bdy(N_PROBES,0.0) {}
	Floats area,dx,dy,barea,bdx,bdy;
};

struct MultiProbePoseAccumulator : public Accumulator {
	MultiProbePoseAccumulator( core::pose::Pose & _pose, std::string tag="" );
	virtual void accumulate_area( core::id::AtomID id, core::Real area, bool buried ) {
		if ( buried ) atom_map_[id].barea[pr_idx_] += area;
		else         atom_map_[id].area[pr_idx_] += area;
	}
	virtual void accumulate_dxdy( core::id::AtomID id, core::Real dx, core::Real dy, bool buried ) {
		if ( buried ) {
			atom_map_[id].bdx[pr_idx_] += dx;
			atom_map_[id].bdy[pr_idx_] += dy;
		} else {
			atom_map_[id].dx[pr_idx_] += dx;
			atom_map_[id].dy[pr_idx_] += dy;
		}
	}
	void set_pr_idx( Size pr_idx ) { pr_idx_ = pr_idx; }
	core::Size get_pr_idx() { return pr_idx_; }
	void show( std::ostream & out );
	core::id::AtomID_Map<LR_MP_AtomData> atom_map_;
	core::Size pr_idx_;
	core::pose::Pose pose_;
	std::string tag_;
};
typedef utility::pointer::shared_ptr<MultiProbePoseAccumulator> MultiProbePoseAccumulatorOP;


struct MultiProbePerSphereAccumulator : public Accumulator {
	MultiProbePerSphereAccumulator( PosePackDataOP pd ) : pr_idx_(1), pd_(pd) {
		atom_map_.resize(pd->spheres.size());
		for ( Size i = 1; i <= pd->spheres.size(); ++i ) atom_map_.resize(i,1);
	}
	// ~MultiProbePerSphereAccumulator() { std::cerr << "delete MultiProbePerSphereAccumulator" << std::endl;}
	PackingScoreResDataOP compute_surrounding_sasa( XYZ & center );
	virtual void accumulate_area( core::id::AtomID id, core::Real area, bool ) {
		atom_map_[id].area[pr_idx_] += area;
	}
	virtual void accumulate_dxdy( core::id::AtomID id, core::Real dx, core::Real dy, bool ) {
		atom_map_[id].dx[pr_idx_] += dx;
		atom_map_[id].dy[pr_idx_] += dy;
	}
	//dz// virtual void accumulate_dz( core::id::AtomID id, core::Real dz ) {
	//dz//  atom_map_[id].dz[pr_idx_] += dz;
	//dz// }
	void set_pr_idx( Size pr_idx ) {
		// std::cerr << "set_pr_idx " << pr_idx << std::endl;
		pr_idx_ = pr_idx;
	}
	void show( std::ostream & out ) {
		Reals tot(31,0.0);
		for ( core::Size i = 1; i <= atom_map_.size(); ++i ) {
			out << "ATOM_DAT " << i << " " << pd_->spheres[i].radius << " ";
			for ( core::Size j = 1; j <= 31; ++j ) {
				tot[j] += atom_map_[core::id::AtomID(1,i)].area[j];
				out << atom_map_[core::id::AtomID(1,i)].area[j] << " ";
			}
			out << " " << std::endl;
		}
		for ( core::Size j = 1; j <= 31; ++j ) std::cerr << j << " " << tot[j] << std::endl;
	}
	core::Size get_pr_idx() { return pr_idx_; }
	core::id::AtomID_Map<LR_MP_AtomData> atom_map_;
	core::Size pr_idx_;
	PosePackDataOP pd_;
};
typedef utility::pointer::shared_ptr<MultiProbePerSphereAccumulator> MultiProbePerSphereAccumulatorOP;

struct Point {
	Point(core::Real _x,core::Real _y) : x(_x),y(_y) {}
	core::Real x,y;
};

struct PointPair {
	PointPair( core::Real x0, core::Real y0, core::Real x1, core::Real y1 ) : a(x0,y0), b(x1,y1) {}
	Point a,b;
};

// struct Arc {
//  Arc( Circle *_circle, core::Real _start_angle, core::Real _end_angle, bool _ccw )
//    : circle(_circle), start_angle(_start_angle), end_angle(_end_angle), ccw(_ccw)
//    {}
//  Circle *circle;
//  core::Real start_angle,end_angle;
//  bool ccw;
// };

struct Circle {
	Circle( core::Real _x, core::Real _y, core::Real _r, core::Real _drdz, core::Real _dada, core::id::AtomID _atom_id )
	: x(_x), y(_y), r(_r), drdz(_drdz), dada(_dada), tcw(NULL), tccw(NULL), atom(_atom_id) {}

	PointPair
	overlap(
		Circle * other
	);

	inline
	bool
	does_overlap(
		Circle * other
	) {
		return sqrt((x-other->x)*(x-other->x) + (y-other->y)*(y-other->y)) < r + other->r;
	}

	inline
	bool
	contains(
		core::Real _x,
		core::Real _y
	) {
		return (x-_x)*(x-_x) + (y-_y)*(y-_y) <= r*r;
	}

	core::Real x,y,r,drdz,dada; // data = d(area)/d(arc angle)
	trace *tcw,*tccw;
	core::id::AtomID atom;


};

std::ostream & operator<< ( std::ostream & out, Circle & circle );

enum EventType {
	ENTER,EXIT,ISECT
};

struct Event {
	Event(
		core::Real _x,
		core::Real _y,
		EventType _kind,
		Circle *_circle,
		Circle *_ccw=NULL
	) : x(_x),
		y(_y),
		cw_angle(-1),
		ccw_angle(-1),
		kind(_kind),
		circle(_circle),
		ccw(_ccw),
		trace_(NULL)
	{
		// std::cerr << "Event " << _kind << " " << circle << std::endl;
		if ( kind == ISECT ) {
			cw_angle  = acos(std::max(-1.0,std::min(1.0,(x-circle->x)/circle->r)));
			ccw_angle = acos(std::max(-1.0,std::min(1.0,(x-ccw->x)/ccw->r)));
			debug_assert( !utility::isnan(cw_angle) );
			debug_assert( !utility::isnan(ccw_angle) );
			if ( y > circle->y )  cw_angle = 2.0*numeric::constants::d::pi-cw_angle;
			if ( y >    ccw->y ) ccw_angle = 2.0*numeric::constants::d::pi-ccw_angle;
		}
	}

	inline
	int
	cmp(
		Event *other
	) {
		if ( x == other->x ) return 0;
		if ( x < other->x ) return -1;
		return 1;
	}

	core::Real x,y,cw_angle,ccw_angle;
	EventType kind;
	Circle *circle,*ccw;
	trace *trace_;

};

struct Slice {

	Slice( AccumulatorOP _accum, Real thickness, bool internal_allowed )
	: accum_(_accum), thickness_(thickness), internal_allowed_(internal_allowed) {}

	virtual ~Slice();

	void compute() {
		compute_events();
		compute_surface();
		compute_derivatives();
	}

	inline void add_circle( Circle *circ ) {
		circles_.push_back(circ);
	}

	Circles circles_;
	traces traces_;
	Events events_;
	AccumulatorOP accum_;
	Real thickness_;
	bool internal_allowed_;

private:

	void compute_surface();
	void compute_events();
	void compute_derivatives();

};

struct trace {

	trace(
		AccumulatorOP _accum,
		Event *e,
		bool _ccw,
		Circle *_circle,
		core::Real _angle,
		trace *other_trace
	) : accum_(_accum), ccw(_ccw), start_(e),
		circle(_circle), angle(_angle), next_trace_(NULL)
	{
		e->trace_ = this;
		if ( ccw ) {
			circle->tccw = this;
			if ( other_trace ) other_trace->next_trace_ = this;
		} else {
			circle->tcw  = this;
			if ( angle==0.0 ) angle = 2*numeric::constants::d::pi;
			if ( other_trace ) next_trace_ = other_trace;
		}
	}

	inline void next_circle( Circle *newcircle, core::Real end_angle, core::Real begin_angle )
	{
		finish_arc(end_angle);
		if ( ccw ) circle->tccw = NULL;
		else    circle->tcw  = NULL;
		circle = newcircle;
		angle = begin_angle;
		if ( ccw ) circle->tccw = this;
		else    circle->tcw  = this;
	}

	inline void end( trace *merge_trace, core::Real end_angle = numeric::constants::d::pi ) {
		finish_arc(end_angle);
		if ( ccw ) next_trace_ = merge_trace;
		else    merge_trace->next_trace_ = this;
	}

	Event *get_first( trace *stop ){
		if ( stop == this || next_trace_==NULL ) return start_;
		Event *best = next_trace_->get_first(stop);
		if ( start_->cmp(best) > 0 ) return start_;
		else                        return best;
	}

	void set_first( trace *stop, Event *e ){
		start_ = e;
		if ( stop == this || next_trace_==NULL ) return;
		next_trace_->set_first( stop, e );
		next_trace_ = NULL;
	}

	AccumulatorOP accum_;
	bool ccw;
	Event *start_;
	Circle *circle;
	core::Real angle;
	trace *next_trace_;
	Arcs arcs_;

private:

	inline void finish_arc( core::Real end_angle ) {
		//debug_assert( fabs(angle-end_angle) * circle->dada < 1.0 );
		// accum_->accumulate_area( circle->atom, fabs(angle-end_angle) * circle->dada );
		arcs_.push_back( Arc( circle->atom, fabs(angle-end_angle) * circle->dada ) );
		//dz//accum_->accumulate_dz(   circle->atom, fabs(angle-end_angle) * circle->drdz );
		// slice->arcs_.push_back(new Arc(circle,angle,end_angle,ccw));
	}

};

struct Array2D {

	Array2D() : Nx_(0), Ny_(0), array_(0) {} // Needed for uninitialized variable errors

	void init( core::Size Nx, core::Size Ny ) {
		Nx_ = Nx;
		Ny_ = Ny;
		array_ = new Circles[Nx*Ny];
	}

	~Array2D() { if ( array_ ) delete[] array_; }

	inline Circles &
	operator() ( core::Size i, core::Size j ) {
		debug_assert( i < Nx_ );
		debug_assert( j < Ny_ );
		debug_assert( array_ );
		// std::cerr << Ny_*i+j << " " << Nx_*Ny_ << std::endl;
		return array_[Ny_*i+j];
	}

	core::Size Nx_;
	core::Size Ny_;
	Circles *array_;

};

struct Octree2D {
	Octree2D( Circles & items ) {
		mnx_=9e9,mny_=9e9,mxx_=-9e9,mxy_=-9e9,mxr_=1.0;
		for ( CircleIter i = items.begin(); i != items.end(); ++i ) {
			Circle* item(*i);
			if ( item->r > mxr_ ) mxr_ = item->r;
			if ( item->x > mxx_ ) mxx_ = item->x;
			if ( item->x < mnx_ ) mnx_ = item->x;
			if ( item->y > mxy_ ) mxy_ = item->y;
			if ( item->y < mny_ ) mny_ = item->y;
		}
		// std::cerr << mnx_ << " " << mxx_ << " " << mny_ << " " << mxy_ << " " << mxr_ << std::endl;
		Xdim_ = get_i(mxx_)+1;
		Ydim_ = get_j(mxy_)+1;
		cubes.init( Xdim_+1, Ydim_+1 );
		// std::cerr << "Octree2D size " << Xdim_+1 << " " << Ydim_+1 << " " << items.size() << " " << mxx_ << " " << mnx_ << " " << mxr_ << std::endl;
		for ( CircleIter i = items.begin(); i != items.end(); ++i ) {
			bool skipcirc(false);
			Circle* item(*i);
			Circles & cube(get_cube(item));
			for ( CircleIter j = cube.begin(); j != cube.end(); ++j ) {
				Circle *item2(*j);
				if ( item->x == item2->x && item->y == item2->y && item->r == item2->r ) {
					skipcirc=true;
					break;
				}
			}
			if ( !skipcirc ) cube.push_back(item);
		}
	}

	inline core::Size get_i(core::Real x) const { return Size((x-mnx_)/mxr_); }
	inline core::Size get_j(core::Real y) const { return Size((y-mny_)/mxr_); }

	inline Circles & get_cube( Circle * a ) { return cubes(get_i(a->x),get_j(a->y)); }
	// inline Circles & get_cube( Event* e ) { return cubes(get_i(e->x),get_j(e->y)); }

	bool
	contains( core::Real x, core::Real y, int extend = 1,
		Circle *exclude1 = NULL, Circle *exclude2 = NULL )
	{
		int I2 = get_i(x);
		int J2 = get_j(y);
		for ( int i2 = std::max(I2-extend,0); i2 <= std::min(I2+extend,(int)Xdim_); ++i2 ) {
			for ( int j2 = std::max(J2-extend,0); j2 <= std::min(J2+extend,(int)Ydim_); ++j2 ) {
				for ( CircleIter ja2 = cubes(i2,j2).begin(); ja2 != cubes(i2,j2).end(); ja2++ ) {
					Circle *circle(*ja2);
					if ( circle == exclude1 || circle == exclude2 ) continue;
					if ( circle->contains(x,y) ) {
						return true;
					}
				}
			}
		}
		return false;
	}

	core::Real mnx_,mny_,mxx_,mxy_,mxr_;
	core::Size Xdim_,Ydim_;
	Array2D cubes;
};

struct LeeRichards {

	LeeRichards(
		PosePackDataOP pd,
		AccumulatorOP accum,
		core::Real spacing = 0.33,
		core::Real probe_radius = 1.4,
		bool csa = false,
		numeric::xyzVector<core::Real> plane = numeric::xyzVector<core::Real>(0,0,1)
	);

	LeeRichards(
		core::pose::Pose & pose,
		AccumulatorOP accum,
		core::Real spacing = 0.33,
		core::Real probe_radius = 1.4,
		bool csa = false,
		numeric::xyzVector<core::Real> plane = numeric::xyzVector<core::Real>(0,0,1)
	);

	~LeeRichards() {
		for ( SliceIter i = slices_.begin(); i != slices_.end(); ++i ) delete *i;
	}

	inline core::Size
	find_first_slice( core::Real r ) {
		if ( r >= slice_coords_.back() ) return 0;
		return find_first_slice(r,1,slice_coords_.size());
	}


	core::Real get_surface_area() const { return surface_area_; }

private:

	core::Real surface_area_;

	utility::vector1< Circles > circles_;
	Slices slices_;
	Reals slice_coords_;

	inline core::Size
	find_first_slice( core::Real r, core::Size beg, core::Size end ) const {
		// std::cerr << "find_first_slice " << r << " " << beg << " " << end << std::endl;
		if ( end - beg < 5 ) { // once close, just iterate up
			while ( slice_coords_[beg] <= r ) ++beg;
			return beg;
		}
		core::Size splt = (end-beg)/2+beg;
		Real t = slice_coords_[splt];
		if ( t <= r ) return find_first_slice(r,splt+1, end);
		else         return find_first_slice(r, beg  ,splt);
	}

	void
	compute(
		Spheres & spheres,
		AccumulatorOP accum,
		core::Real spacing = 0.33,
		core::Real probe_radius = 1.4,
		bool csa = false,
		numeric::xyzVector<core::Real> plane = numeric::xyzVector<core::Real>(0,0,1)
	);

};

core::Real
compute_packing_score_leerichards(
	PosePackDataOP pd,
	core::Real slicesize,
	numeric::xyzVector<Real> plane = numeric::xyzVector<Real>(0,0,1)
);

core::Real
compute_surface_area_leerichards(
	Real & buried_area_out,
	PosePackDataOP pd,
	Real slicesize,
	Real pr,
	bool csa = false,
	numeric::xyzVector<Real> plane = numeric::xyzVector<Real>(0,0,0)
);

core::Real
compute_surface_area_leerichards(
	PosePackDataOP pd,
	Real slicesize,
	Real pr,
	bool csa = false,
	numeric::xyzVector<Real> plane = numeric::xyzVector<Real>(0,0,0)
);

XYZs
compute_surface_area_leerichards_deriv(
	PosePackDataOP pd,
	Real slicesize,
	Real pr,
	bool csa = false//,
	// numeric::xyzVector<Real> plane = numeric::xyzVector<Real>(0,0,1)
);

XYZs
check_surface_area_leerichards_deriv(
	PosePackDataOP pd,
	Real slicesize,
	Real pr,
	bool csa = false,
	Size max_num = 10
	// numeric::xyzVector<Real> plane = numeric::xyzVector<Real>(0,0,1)
);


} // namespace packstat
} // namespace scoring
} // namespace core


#endif
