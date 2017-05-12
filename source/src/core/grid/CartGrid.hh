// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/grid/Grid.hh
/// @author Sam DeLuca

#ifndef INCLUDED_core_grid_CartGrid_hh
#define INCLUDED_core_grid_CartGrid_hh


#include <core/grid/CartGrid.fwd.hh>
#include <core/types.hh>

#include <numeric/xyzVector.hh>

#include <utility/vector1.hh>
#include <utility/vector0.hh>
#include <utility/tools/make_vector.hh>
#include <utility/string_util.hh>
#include <utility/json_spirit/json_spirit_value.h>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector0.fwd.hh>
#include <utility/io/izstream.hh>
#include <utility/exit.hh>
#include <utility/Binary_Util.hh>

#include <ObjexxFCL/string.functions.hh>

#include <algorithm>


namespace core {
namespace grid {

template <typename T>
class CartGrid : public utility::pointer::ReferenceCount
{
public:
	/// @brief This needs to be an int, as it can hold negative values.
	typedef numeric::xyzVector< int > GridPt;

	CartGrid():
		nX_(0), nY_(0), nZ_(0),
		lX_(0.0), lY_(0.0), lZ_(0.0),
		bX_(0.0), bY_(0.0), bZ_(0.0),
		tX_(0.0), tY_(0.0), tZ_(0.0),
		name_("default"),
		npoints_(0),
		fullyOccupied_(false),
		zones_()
	{}

	CartGrid( CartGrid<T> const & other ) {
		other.clone( *this );
	}

	~CartGrid() override = default;

	/*
	std::ostream & operator<< (Grid<T> const & mygrid)
	{
	os <<mygrid.its_bX_ << " " << mygrid.its_bY_ << " " << mygrid.its_bZ_ << " | "
	<< mygrid.its_tX_ << " " << mygrid.its_tY_ << " " << mygrid.its_tZ_;
	return os;
	}
	*/

	void setBase(core::Real x, core::Real y, core::Real z)
	{
		bX_ = x;
		bY_ = y;
		bZ_ = z;

		this->setTop();
	}

	void setDimensions(core::Size nX, core::Size nY, core::Size nZ, core::Real lX, core::Real lY, core::Real lZ)
	{
		nX_ = nX;
		nY_ = nY;
		nZ_ = nZ;

		lX_ = lX;
		lY_ = lY;
		lZ_ = lZ;

		this->setTop();
	}

	void set_name(std::string const & name)
	{
		name_ = name;
	}

	std::string get_name() const
	{
		return name_;
	}

	core::Size longestSide() const
	{
		core::Size ls = this->nX_;
		if ( this->nY_ > ls ) ls = this->nY_;
		if ( this->nZ_ > ls ) ls = this->nZ_;

		return ls;
	}
	bool equalDimensions(CartGrid<T> const & rhs) const
	{
		if ( nX_ != rhs.nX_ ) {
			return false;
		}
		if ( nY_ != rhs.nY_ ) {
			return false;
		}
		if ( nZ_ != rhs.nZ_ ) {
			return false;
		}

		if ( lX_ != rhs.lX_ ) {
			return false;
		}
		if ( lY_ != rhs.lY_ ) {
			return false;
		}
		if ( lZ_ != rhs.lZ_ ) {
			return false;
		}

		return true;
	}

	bool equalBase(CartGrid<T> const & rhs) const
	{
		if ( std::abs(this->bX_ - rhs.bX_) > 0.01 ) return false;
		if ( std::abs(this->bY_ - rhs.bY_) > 0.01 ) return false;
		if ( std::abs(this->bZ_ - rhs.bZ_) > 0.01 ) return false;
		return true;
	}

	bool is_in_grid(core::Real x, core::Real y, core::Real z) const
	{
		if ( x < bX_ || x > tX_ ) {
			return false;
		}
		if ( y < bY_ || y > tY_ ) {
			return false;
		}
		if ( z < bZ_ || z > tZ_ ) {
			return false;
		}

		return true;
	}

	bool setupZones()
	{
		npoints_ = nX_*nY_*nZ_;
		// Make zones_ npoints_ copies of the default value of the type
		zones_.assign( npoints_, T() );
		return true;
	}

	void translate(core::Real x, core::Real y, core::Real z)
	{
		bX_ += x;
		bY_ += y;
		bZ_ += z;
		this->setTop();
	}

	void setValue(core::Size ix, core::Size iy, core::Size iz, T value)
	{
		core::Size index = this->get_index(ix, iy, iz);
		setValue( index, value );
	}

	void setValue(core::Real fx, core::Real fy, core::Real fz, T value)
	{
		core::Size index = this->get_index(fx,fy,fz);
		setValue( index, value );
	}

	T getValue(core::Size ix, core::Size iy, core::Size iz) const
	{
		core::Size index = this->get_index(ix, iy, iz);
		return this->getValue(index);
	}

	T getValue(core::Real fx, core::Real fy, core::Real fz) const
	{
		core::Size index = this->get_index(fx, fy, fz);
		return this->getValue(index);
	}

	void zero()
	{
		zones_.assign( npoints_, 0 );
	}

	void setFullOccupied(T value)
	{
		zones_.assign( npoints_, value );
		fullyOccupied_ = true;
	}

	void clone(CartGrid<T> & copy) const {
		// copy over gross information
		copy.nX_ = this->nX_;
		copy.nY_ = this->nY_;
		copy.nZ_ = this->nZ_;

		copy.lX_ = this->lX_;
		copy.lY_ = this->lY_;
		copy.lZ_ = this->lZ_;

		copy.bX_ = this->bX_;
		copy.bY_ = this->bY_;
		copy.bZ_ = this->bZ_;

		copy.tX_ = this->tX_;
		copy.tY_ = this->tY_;
		copy.tZ_ = this->tZ_;

		copy.name_ = this->name_;

		copy.setupZones();

		for ( core::Size i=0; i < this->npoints_; i++ ) {
			copy.zones_[i] = this->zones_[i];
		}
	}

	void reset_boundaries() {
		core::Size minx = this->nX_;
		core::Size miny = this->nY_;
		core::Size minz = this->nZ_;

		core::Size maxx = 0;
		core::Size maxy = 0;
		core::Size maxz = 0;

		for ( core::Size ix=0; ix < this->nX_; ix++ ) {
			for ( core::Size iy=0; iy < this->nY_; iy++ ) {
				for ( core::Size iz=0; iz < this->nZ_; iz++ ) {
					if ( this->getValue(ix,iy,iz) ) {
						minx = std::min(minx,ix);
						miny = std::min(miny,iy);
						minz = std::min(minz,iz);

						maxx = std::max(maxx,ix);
						maxy = std::max(maxy,iy);
						maxz = std::max(maxz,iz);
					}
				}
			}
		}

		if ( maxx*maxy*maxz == 0 ) {
			return;   // no change
		}

		CartGrid<T> tmpgrid;
		this->clone(tmpgrid);

		core::Size nx = maxx - minx + 1;
		core::Size ny = maxy - miny + 1;
		core::Size nz = maxz - minz + 1;

		this->nX_ = nx;
		this->nY_ = ny;
		this->nZ_ = nz;

		this->bX_ += minx*this->lX_;
		this->bY_ += miny*this->lY_;
		this->bZ_ += minz*this->lZ_;

		this->setupZones();

		for ( core::Size i=0; i < nx; i++ ) {
			for ( core::Size j=0; j < ny; j++ ) {
				for ( core::Size k=0; k < nz; k++ ) {
					this->setValue( i, j, k, tmpgrid.getValue(minx+i, miny+j, minz+k) );
				}
			}
		}
	}

	void read(std::string const & filename) {
		//std::ifstream file;
		utility::io::izstream file;
		std::istringstream line_stream;

		file.open(filename.c_str());

		if ( !file ) {
			std::cout << "read_gridfile - unable to open gridfile:" << filename << std::endl;
			std::exit( EXIT_FAILURE );
		}

		std::string line;
		std::string keyword;
		std::string name;
		core::Real bx=0.0, by=0.0, bz=0.0;
		core::Real lx=0.0, ly=0.0, lz=0.0;
		T occupied(0);
		core::Size nx=0, ny=0, nz=0;
		core::Size ix=0, iy=0;

		while ( file ) {
			getline(file, line);


			if ( ObjexxFCL::is_blank(line) ) {
				ix++;
				iy=0;
				continue;
			}

			line_stream.clear();
			line_stream.str(line);
			line_stream.seekg( std::ios::beg );

			line_stream >> keyword;
			if ( keyword == "NAME:" ) {
				line_stream >> name;
				this->set_name(name);
			} else if ( keyword == "BASE:" ) {
				line_stream >> bx >> by >> bz;
				this->setBase(bx, by, bz);
			} else if ( keyword == "SIZE:" ) {
				line_stream >> nx >> ny >> nz;
			} else if ( keyword == "LENGTH:" ) {
				line_stream >> lx >> ly >> lz;
				this->setDimensions(nx, ny, nz, lx, ly, lz);
				this->setupZones();
			} else {
				line_stream.seekg( std::ios::beg );
				for ( core::Size iz=0; iz < nz; iz++ ) {
					line_stream >> occupied;
					this->setValue(ix,iy,iz,occupied);
				}
				iy++;
			}
		}

		file.close();
	}

	void write(std::string const & filename) const {
		std::ofstream file;
		file.open(filename.c_str());
		file << "NAME: " << this->get_name() << std::endl;
		file << "BASE: " << this->bX_ << " " << this->bY_ << " " << this->bZ_ << std::endl;
		file << "SIZE: " << this->nX_ << " " << this->nY_ << " " << this->nZ_ << std::endl;
		file << "LENGTH: " << this->lX_ << " " << this->lY_ << " " << this->lZ_ << std::endl;

		for ( core::Size i=0; i < this->nX_; i++ ) {
			for ( core::Size j=0; j < this->nY_; j++ ) {
				for ( core::Size k=0; k < this->nZ_; k++ ) {
					file << this->getValue(i,j,k) << " ";
				}
				file << std::endl;
			}
			file << std::endl;
		}
		file.close();
	}

	bool isFullyOccupied() const {
		return fullyOccupied_;
	}

	bool isEmpty() const {
		for ( core::Size i=0; i < npoints_; i++ ) {
			if ( this->zones_[i] != 0 ) {
				return false;
			}
		}
		return true;
	}

	utility::json_spirit::Value serialize() const
	{
		using utility::json_spirit::Pair;
		using utility::json_spirit::Value;
		using utility::json_spirit::Array;
		Pair name("name",this->get_name());

		Pair base("base",utility::tools::make_vector(Value(this->bX_),Value(this->bY_),Value(this->bZ_)));
		// Need explicit conversion to int here, as utility::json_spirit::Value doesn't have a constructor for size_t.
		Pair size("size",utility::tools::make_vector(Value(int(this->nX_)),Value(int(this->nY_)),Value(int(this->nZ_))));
		Pair length("length",utility::tools::make_vector(Value(this->lX_),Value(this->lY_),Value(this->lZ_)));

		std::string point_data;
		// vector::data() gives a raw pointer to the underlying (contigous) array.
		debug_assert( npoints_ == zones_.size() );
		utility::encode6bit( (unsigned char*)zones_.data(), npoints_*sizeof(T), point_data );

		utility::json_spirit::Pair data("data",Value(point_data));

		return utility::json_spirit::Value( utility::tools::make_vector(name,base,size,length,data) );

	}

	void deserialize(utility::json_spirit::mObject grid_data)
	{
		std::string name = grid_data["name"].get_str();

		utility::json_spirit::mArray base_data = grid_data["base"].get_array();
		debug_assert(base_data.size() == 3);

		core::Real bX = base_data[0].get_real();
		core::Real bY = base_data[1].get_real();
		core::Real bZ = base_data[2].get_real();

		utility::json_spirit::mArray size_data = grid_data["size"].get_array();
		debug_assert(size_data.size() == 3);

		int nX = size_data[0].get_int();
		int nY = size_data[1].get_int();
		int nZ = size_data[2].get_int();

		utility::json_spirit::mArray length_data = grid_data["length"].get_array();
		debug_assert(length_data.size() == 3);

		core::Real lX = length_data[0].get_real();
		core::Real lY = length_data[1].get_real();
		core::Real lZ = length_data[2].get_real();

		this->set_name(name);
		this->setBase(bX,bY,bZ);
		this->setDimensions(nX, nY, nZ, lX, lY, lZ);
		this->setupZones();

		std::string point_data = grid_data["data"].get_str();

		// Why do we do the dance where we resize too big, decode, and then re-resize?
		// Because we use a function which does a 4->3 transformation on the input data
		// ... which means we may run off the end of the array if the size of the array isn't divisible by 3
		// By adding an extra 2 items padding, we give it a safe space to overrun.
		zones_.resize( npoints_+2 );
		debug_assert( sizeof(T)*zones_.size()*4 >= point_data.size()*3 ); // 3 bytes of array data for every 4 bytes of string data
		// vector::data() gives a raw pointer to the underlying (contigous) array.
		utility::decode6bit( (unsigned char*) zones_.data(), point_data );
		zones_.resize( npoints_ );

	}

	void sum(utility::vector0<utility::pointer::shared_ptr<CartGrid<T> > > const & list_grids) {
		core::Size ngrids = list_grids.size();

		this->zero();
		T curr_value, new_value;
		for ( core::Size i=0; i < ngrids; i++ ) {
			if ( !this->equalDimensions(*(list_grids[i])) ) {
				return;
			}
			if ( !this->equalBase(*(list_grids[i])) ) {
				return;
			}

			for ( core::Size j=0; j < this->its_npoints; j++ ) {
				curr_value = this->getValue(j);
				new_value = list_grids[i]->getValue(j);
				curr_value += new_value;
				this->setValue(j,curr_value);
			}
		}
	}

	void expand(int expansion) {
		this->bX_ -= this->lX_*core::Real(expansion);
		this->bY_ -= this->lY_*core::Real(expansion);
		this->bZ_ -= this->lZ_*core::Real(expansion);

		this->nX_ += 2*expansion;
		this->nY_ += 2*expansion;
		this->nZ_ += 2*expansion;

		this->setupZones();
	}

	core::Vector getBase() const
	{
		return core::Vector(bX_, bY_, bZ_);
	}

	core::Vector getTop() const
	{
		return core::Vector(tX_, tY_, tZ_);
	}

	void getNumberOfPoints(core::Size & x, core::Size & y, core::Size & z) const
	{
		x = nX_;
		y = nY_;
		z = nZ_;
	}

	core::Vector getNumberOfPoints() const
	{
		return core::Vector(nX_, nY_, nZ_);
	}

	GridPt gridpt(Vector const & coords) const
	{
		// Round down.
		// Just casting to int will give wrong values for points outside the grid.
		// May return a negative number if it's outside of the box.
		return GridPt(
			static_cast<int>(std::floor((coords.x() - bX_)/lX_)),
			static_cast<int>(std::floor((coords.y() - bY_)/lY_)),
			static_cast<int>(std::floor((coords.z() - bZ_)/lZ_))
		);

	}

	Vector coords(GridPt const & gridpt) const
	{
		return Vector(
			(bX_ + (gridpt.x() + 0.5)*lX_),
			(bY_ + (gridpt.y() + 0.5)*lY_),
			(bZ_ + (gridpt.z() + 0.5)*lZ_)
		);
	}


	T getValue(GridPt const & gridpt) const
	{
		return getValue(gridpt.x(), gridpt.y(), gridpt.z());
	}

	T getValue(Vector const & coords) const
	{
		return getValue(coords.x(), coords.y(), coords.z());
	}

	T getMinValue() const
	{
		return *std::min_element(zones_.begin(),zones_.end());
	}

	T getMaxValue() const
	{
		return *std::max_element(zones_.begin(),zones_.end());
	}

	void setValue(GridPt const & gridpt, T value)
	{
		setValue(gridpt.x(), gridpt.y(), gridpt.z(), value);
	}

	void setValue(Vector const & coords, T value)
	{
		setValue(coords.x(), coords.y(), coords.z(), value);
	}

	/// @details This format was choosen because it's space-efficient for small
	/// integer values (such as are typically stored in grids) and PyMOL can read it.
	/// Typical extension is .brix or .omap
	void write_to_BRIX(std::string const & filename) const
	{
		std::ofstream out( filename.c_str() );
		write_to_BRIX(out);
		out.close();
	}

	void write_to_BRIX(std::ostream & out) const
	{
		int const plus = 127; // to convert values to unsigned bytes
		// Crystallographic maps always expect one grid point to be at the origin.
		// This may shift our map slightly...
		GridPt origin = gridpt(Vector(0.,0.,0.));
		out << ":-)"; // magic number
		out << " Origin " << -origin.x() << ' ' << -origin.y() << ' ' << -origin.z() ; // starting indices for grid points actually present
		out << " Extent " << nX_ << ' ' << nY_ << ' ' << nZ_ ; // number of grid points in each dimension
		out << " Grid " << nX_ << ' ' << nY_ << ' ' << nZ_ ; // number of grid points in one unit cell
		out << " Cell " << nX_*lX_ << ' ' << nY_*lY_ << ' ' << nZ_*lZ_ << " 90.0 90.0 90.0"; // dimensions of unit cell
		out << " Prod 1 Plus " << plus << " Sigma 1\f";
		// Now pad with spaces until we've writen 512 characters:
		for ( long i = out.tellp(); i < 512; ++i ) out << ' ';
		// Data stored in (padded) 8x8x8 blocks with X fast, Y medium, Z slow
		typedef unsigned char ubyte;
		for ( core::Size zz = 0; zz < nZ_; zz += 8 ) {
			for ( core::Size yy = 0; yy < nY_; yy += 8 ) {
				for ( core::Size xx = 0; xx < nX_; xx += 8 ) {
					for ( core::Size z = zz, z_end = zz+8; z < z_end; ++z ) {
						for ( core::Size y = yy, y_end = yy+8; y < y_end; ++y ) {
							for ( core::Size x = xx, x_end = xx+8; x < x_end; ++x ) {
								if ( x < nX_ && y < nY_ && z < nZ_ ) out << ubyte( getValue(x, y, z) + plus );
								else out << ubyte( 0 + plus );
							}
						}
					}
				}
			}
		}
	}

private:
	void setTop()
	{
		tX_ = bX_ + nX_*lX_;
		tY_ = bY_ + nY_*lY_;
		tZ_ = bZ_ + nZ_*lZ_;
	}

	bool valid_index(core::Size index) const
	{
		return (index >= 0) && (index < npoints_);
	}

	// setValue for GridPt
	void setValue(int ix, int iy, int iz, T value)
	{
		debug_assert( ix >= 0 && iy >= 0 && iz >= 0 );
		core::Size index = this->get_index(core::Size(ix), core::Size(iy), core::Size(iz));
		setValue( index, value );
	}

	void setValue(core::Size index, T value)
	{
		if ( ! valid_index( index ) ) {
			utility_exit_with_message( "Cannot set value for point not in grid." );
		}

		if ( value == 0 ) {
			fullyOccupied_ = false;
		}
		this->zones_[index] = value;
	}

	// getValue for GridPt
	T getValue(int ix, int iy, int iz) const
	{
		debug_assert( ix >= 0 && iy >= 0 && iz >= 0 );
		core::Size index = this->get_index(core::Size(ix), core::Size(iy), core::Size(iz));
		return this->getValue(index);
	}

	T getValue(core::Size index) const
	{
		if ( ! valid_index( index ) ) {
			utility_exit_with_message( "Cannot get value for point not in grid." );
		}
		return this->zones_[index];
	}

	core::Size get_index(core::Size ix, core::Size iy, core::Size iz) const
	{
		if ( ix >= nX_ || iy >= nY_ || iz >= nZ_ ) {
			return npoints_+1; // Invalid index
		}

		return ix*(nY_*nZ_) + iy*(nZ_) + iz;
	}

	core::Size get_index(core::Real fx, core::Real fy, core::Real fz) const
	{
		if ( ! is_in_grid(fx,fy,fz) ) {
			return npoints_+1; // Invalid index
		}

		core::Size ix((fx - bX_)/lX_);
		core::Size iy((fy - bY_)/lY_);
		core::Size iz((fz - bZ_)/lZ_);

		return this->get_index(core::Size(ix),core::Size(iy),core::Size(iz));
	}

private:

	core::Size nX_, nY_, nZ_;
	core::Real lX_, lY_,lZ_;
	core::Real bX_, bY_, bZ_;
	core::Real tX_,tY_,tZ_;
	std::string name_;
	core::Size npoints_;
	bool fullyOccupied_;
	utility::vector0< T > zones_;

};

}
}


#endif /* GRID_HH_ */
