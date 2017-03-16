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
	typedef numeric::xyzVector< int > GridPt;

	CartGrid():
		nX_(0), nY_(0), nZ_(0),
		lX_(0.0), lY_(0.0), lZ_(0.0),
		bX_(0.0), bY_(0.0), bZ_(0.0),
		tX_(0.0), tY_(0.0), tZ_(0.0),
		name_("default"),
		npoints_(0),
		fullyOccupied_(false),
		zones_(NULL)
	{

	}

	~CartGrid() override
	{
		if ( zones_ !=NULL ) {
			delete [] zones_;
			zones_ = NULL;
		}
	}

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

	void setDimensions(int nX, int nY, int nZ, core::Real lX, core::Real lY, core::Real lZ)
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

	int longestSide() const
	{
		int ls = this->nX_;
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
		if ( zones_ != NULL ) {
			//delete zones_;
			delete[] zones_;
			zones_ = NULL;
		}

		npoints_ = nX_*nY_*nZ_;
		// Why npoints_+2 ?
		// Because in deserialize() we use a function which does a 4->3 transformation on the input data
		// ... which means we may run off the end of the array if the size of the array isn't divisible by 3
		// By adding an extra 2 items padding, we give it a safe space to overrun.
		// (All other uses include serialize() use npoints_ to know how big the array is, so won't use the extra space.)
		zones_ = new T[npoints_+2];

		return true;
	}

	void translate(core::Real x, core::Real y, core::Real z)
	{
		bX_ += x;
		bY_ += y;
		bZ_ += z;
		this->setTop();
	}

	bool setValue(int ix, int iy, int iz, T value)
	{
		if ( value == 0 ) {
			fullyOccupied_ = false;
		}

		int index = this->get_index(ix, iy, iz);
		debug_assert( index >= 0 && index < npoints_ );
		this->zones_[index] = value;
		return true;
	}


	bool setValue(core::Real fx, core::Real fy, core::Real fz, T value)
	{
		if ( value == 0 ) {
			fullyOccupied_ = false;
		}

		int ix = int((fx - bX_)/lX_);
		if ( ix < 0 || ix >= nX_ ) return 0;

		int iy = int((fy - bY_)/lY_);
		if ( iy < 0 || iy >= nY_ ) return 0;

		int iz = int((fz - bZ_)/lZ_);
		if ( iz < 0 || iz >= nZ_ ) return 0;

		int index = this->get_index(ix,iy,iz);
		debug_assert( index >= 0 && index < npoints_ );
		this->zones_[index] = value;
		return true;
	}

	T getValue(int ix, int iy, int iz) const
	{
		int index = this->get_index(ix, iy, iz);
		return this->getValue(index);
	}

	T getValue(core::Real fx, core::Real fy, core::Real fz) const
	{
		// --- round down --- //
		int ix = int((fx - bX_)/lX_);
		if ( ix < 0 || ix >= nX_ ) return 0;

		int iy = int((fy - bY_)/lY_);
		if ( iy < 0 || iy >= nY_ ) return 0;

		int iz = int((fz - bZ_)/lZ_);
		if ( iz < 0 || iz >= nZ_ ) return 0;

		int index = ix*(nY_*nZ_) + iy*(nZ_) + iz;
		return this->zones_[index];
	}

	void zero()
	{
		for ( int i=0; i < npoints_; i++ ) {
			this->zones_[i] = 0;
		}
	}

	void setFullOccupied(T value)
	{
		for ( int i=0; i < npoints_; i++ ) {
			this->zones_[i] = value;
		}
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

		for ( int i=0; i < this->npoints_; i++ ) {
			copy.zones_[i] = this->zones_[i];
		}
	}

	void reset_boundaries() {
		int minx = this->nX_;
		int miny = this->nY_;
		int minz = this->nZ_;

		int maxx = 0;
		int maxy = 0;
		int maxz = 0;

		for ( int ix=0; ix < this->nX_; ix++ ) {
			for ( int iy=0; iy < this->nY_; iy++ ) {
				for ( int iz=0; iz < this->nZ_; iz++ ) {
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

		int nx = maxx - minx + 1;
		int ny = maxy - miny + 1;
		int nz = maxz - minz + 1;

		this->nX_ = nx;
		this->nY_ = ny;
		this->nZ_ = nz;

		this->bX_ += minx*this->lX_;
		this->bY_ += miny*this->lY_;
		this->bZ_ += minz*this->lZ_;

		this->setupZones();

		int index2=0;
		int value = 0;
		for ( int i=0; i < nx; i++ ) {
			for ( int j=0; j < ny; j++ ) {
				for ( int k=0; k < nz; k++ ) {
					index2 = tmpgrid.get_index(minx+i, miny+j, minz+k);
					value = tmpgrid.getValue(index2);
					this->setValue(i,j,k,value);
				}
			}
		}
	}

	// its_nX_ is no longer a member... looks liks this functon is now deprecated - commenting out to fix PyRosetta build
	// void fluff(utility::pointer::shared_ptr<CartGrid<T> > input, utility::pointer::shared_ptr<CartGrid<T> > original, int amount=6) {
	//  int istart, iend, jstart, jend, kstart, kend;
	//  for ( int i=0; i < input->its_nX_; i++ ) {
	//   for ( int j=0; j < input->its_nY_; j++ ) {
	//    for ( int k=0; k < input->its_nZ_; k++ ) {
	//     if ( input->getValue(i,j,k) != 0 ) {
	//      istart = std::max(0, (i-amount));
	//      iend   = std::min(input->its_nX_, (i+amount));
	//      jstart = std::max(0, (j-amount));
	//      jend   = std::min(input->its_nY_, (j+amount));
	//      kstart = std::max(0, (k-amount));
	//      kend   = std::min(input->its_nZ_, (k+amount));
	//      for ( int ii=istart; ii < iend; ii++ ) {
	//       for ( int jj=jstart; jj < jend; jj++ ) {
	//        for ( int kk=kstart; kk < kend; kk++ ) {
	//         if ( original->getValue(ii,jj,kk) != 0 ) {
	//          this->setValue(ii,jj,kk,1);
	//         }
	//        }
	//       }
	//      }
	//     }
	//    }
	//   }
	//  }
	// }

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
		int nx=0, ny=0, nz=0;
		int ix=0, iy=0;

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
				for ( int iz=0; iz < nz; iz++ ) {
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

		for ( int i=0; i < this->nX_; i++ ) {
			for ( int j=0; j < this->nY_; j++ ) {
				for ( int k=0; k < this->nZ_; k++ ) {
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
		for ( int i=0; i < npoints_; i++ ) {
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
		Pair size("size",utility::tools::make_vector(Value(this->nX_),Value(this->nY_),Value(this->nZ_)));
		Pair length("length",utility::tools::make_vector(Value(this->lX_),Value(this->lY_),Value(this->lZ_)));

		std::string point_data;
		utility::encode6bit((unsigned char*)zones_,npoints_*sizeof(T),point_data);

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

		debug_assert( sizeof(T)*(npoints_+2)*4 >= point_data.size()*3 ); // 3 bytes of array data for every 4 bytes of string data
		utility::decode6bit((unsigned char*)zones_,point_data);

	}

	void sum(utility::vector0<utility::pointer::shared_ptr<CartGrid<T> > > const & list_grids) {
		int ngrids = static_cast<int>(list_grids.size());

		this->zero();
		T curr_value, new_value;
		for ( int i=0; i < ngrids; i++ ) {
			if ( !this->equalDimensions(*(list_grids[i])) ) {
				return;
			}
			if ( !this->equalBase(*(list_grids[i])) ) {
				return;
			}

			for ( int j=0; j < this->its_npoints; j++ ) {
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

	// its_nX_ is no longer a member... looks liks this functon is now deprecated - commenting out to fix PyRosetta build
	// void split(int nsplits, int igrid, core::Real pad, utility::pointer::shared_ptr<CartGrid<T> > grid) {
	//  int tsplits = nsplits*nsplits*nsplits;
	//  if ( igrid < 0 || igrid >= tsplits ) {
	//   utility_exit_with_message("accessing split out of bounds");
	//  }

	//  int split_x = int(core::Real(this->nX_)/core::Real(nsplits)) + 1;
	//  int split_y = int(core::Real(this->nY_)/core::Real(nsplits)) + 1;
	//  int split_z = int(core::Real(this->nZ_)/core::Real(nsplits)) + 1;

	//  int iz = int(igrid / (nsplits*nsplits));
	//  int rem = igrid - (iz*nsplits*nsplits);
	//  int iy = int(rem/nsplits);
	//  int ix = rem - (iy*nsplits);

	//  int padx = int((pad/lX_)-0.1);
	//  int pady = int((pad/lY_)-0.1);
	//  int padz = int((pad/lZ_)-0.1);

	//  int leftx = padx;
	//  int rightx = padx;
	//  int lefty = pady;
	//  int righty = pady;
	//  int leftz = padz;
	//  int rightz = padz;

	//  if ( ix == 0 ) {
	//   leftx = 0;
	//  }
	//  if ( iy == 0 ) {
	//   lefty = 0;
	//  }
	//  if ( iz == 0 ) {
	//   leftz = 0;
	//  }

	//  if ( ix == nsplits-1 ) {
	//   rightx = 0;
	//  }
	//  if ( iy == nsplits-1 ) {
	//   righty = 0;
	//  }
	//  if ( iz == nsplits-1 ) {
	//   rightz = 0;
	//  }

	//  grid->its_nX_ = split_x+leftx+rightx;
	//  grid->its_nY_ = split_y+lefty+righty;
	//  grid->its_nZ_ = split_z+leftz+rightz;

	//  core::Real bx = this->its_bX + ix*split_x*(this->lX_) - padx*(this->lX_);
	//  core::Real by = this->its_bY + iy*split_y*(this->lY_) - pady*(this->lY_);
	//  core::Real bz = this->its_bZ + iz*split_z*(this->lZ_) - padz*(this->lZ_);

	//  grid->its_bX_ = std::max(this->bX_, bx);
	//  grid->its_bY_ = std::max(this->bY_, by);
	//  grid->its_bZ_ = std::max(this->bZ_, bz);

	//  grid->its_lX_ = this->lX_;
	//  grid->its_lY_ = this->lY_;
	//  grid->its_lZ_ = this->lZ_;


	//  // copy over grid data
	//  grid->setDimensions(grid->its_nX_,grid->its_nY_,grid->its_nZ_,grid->its_lX_,grid->its_lY_,grid->its_lZ_);
	//  grid->setupZones();
	//  grid->zero();

	//  int xstart = std::max(0,ix*split_x-padx);
	//  int ystart = std::max(0,iy*split_y-pady);
	//  int zstart = std::max(0,iz*split_z-padz);

	//  for ( int i=0; i < grid->its_nX_; i++ ) {
	//   for ( int j=0; j < grid->its_nY_; j++ ) {
	//    for ( int k=0; k < grid->its_nZ_; k++ ) {
	//     int value = this->getValue(xstart+i,ystart+j,zstart+k);
	//     if ( value == -1 ) value = 0;
	//     grid->setValue(i,j,k,value);
	//    }
	//   }
	//  }
	// }

	core::Vector getBase() const
	{
		return core::Vector(bX_, bY_, bZ_);
	}

	core::Vector getTop() const
	{
		return core::Vector(tX_, tY_, tZ_);
	}

	void getNumberOfPoints(int & x, int & y, int & z) const
	{
		x = nX_;
		y = nY_;
		z = nZ_;
	}

	GridPt gridpt(Vector const & coords) const
	{
		// Round down.
		// Just casting to int will give wrong values for points outside the grid.
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
		return *std::min_element(zones_,zones_+npoints_);
	}

	T getMaxValue() const
	{
		return *std::max_element(zones_,zones_+npoints_);
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
	void write_to_BRIX(std::string const & filename)
	{
		std::ofstream out( filename.c_str() );
		write_to_BRIX(out);
		out.close();
	}

	void write_to_BRIX(std::ostream & out)
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
		for ( int zz = 0; zz < nZ_; zz += 8 ) {
			for ( int yy = 0; yy < nY_; yy += 8 ) {
				for ( int xx = 0; xx < nX_; xx += 8 ) {
					for ( int z = zz, z_end = zz+8; z < z_end; ++z ) {
						for ( int y = yy, y_end = yy+8; y < y_end; ++y ) {
							for ( int x = xx, x_end = xx+8; x < x_end; ++x ) {
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

	void setValue(int index, T value)
	{
		if ( value == 0 ) {
			fullyOccupied_ = false;
		}
		debug_assert( index >= 0 && index < npoints_ );
		this->zones_[index] = value;
	}

	T getValue(int index) const
	{
		if ( index < 0 || index >= npoints_ ) {
			return -1;
		} else {
			return this->zones_[index];
		}
	}

	int get_index(int ix, int iy, int iz) const
	{
		int index;
		index = ix*(nY_*nZ_) + iy*(nZ_) + iz;
		return index;
	}


private:

	int nX_, nY_, nZ_;
	core::Real lX_, lY_,lZ_;
	core::Real bX_, bY_, bZ_;
	core::Real tX_,tY_,tZ_;
	std::string name_;
	int npoints_;
	bool fullyOccupied_;
	T* zones_;

};

}
}


#endif /* GRID_HH_ */
