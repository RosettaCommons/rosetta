// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/interpolation/Histogram.hh
/// @brief  A class for storing histogram data
/// @author Spencer Bliven <blivens@u.washington.edu>
/// @date   1/23/09

#ifndef INCLUDED_numeric_interpolation_Histogram_hh
#define INCLUDED_numeric_interpolation_Histogram_hh

// Unit Headers
#include <numeric/interpolation/Histogram.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

// Project Headers
#include <numeric/interpolation/interpolation.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/interpolation/spline/SplineGenerator.hh>

// Utility Headers
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>

#include <cmath>
#include <map>
#include <algorithm>

namespace numeric {
namespace interpolation {

using numeric::Real;

/**
* @brief A histogram with fixed-width bins
*
* @details Histograms are commonly used to approximate arbitrary functions.
* At their most basic level they just map X values to f(X) for discrete, regularly
* spaced values of X. Usually you then interpolate between the stored Y values
* so that all values of X can be evaluated.
*
* When creating a histogram the range and bin width must be given. Several
* parameters can also be specified to set how the interpolation will work.
* After that the function can be approximated for arbitrary X values by calling
* interpolate(). The important parameters describing how to interpolate are:
*
*  - Periodicity:
*      - nonperiodic - Only X values within the range of the function are strictly
*          legal. See interpolate(X,Y&) for the behavior when out of this range.
*      - periodic - All X values are taken modulus the length of the range of
*          the function.
*  - Bin Placement:\n
*      Since bins span a range of X, it is ambiguous exactly what X value the
*      bin gives the value of. The choice for BinPlacement should depend on
*      the source of the data for the histogram.\n
*      If bin[x] spans a range [x1,x2],
*      - left - bin[x] corresponds to f(x1). This is streightforward, but you
*          tend to over-estimate f(x) for areas with positive slope and
*          under-estimate for areas with negative slope due to the stair-step
*          shape of the histogram
*      - center - bin[x] corresponds to f( (x1+x2)/2 )
*      - right - bin[x] corresponds to f(x2). Equivalent to left with bin[x+1]
*  - Interpolator:\n
*      Specifies the algorithm used for interpolating between bins
*      - flat - No interpolation. Gives a discontinuous function, but faithful
*          to the raw data.
*      - linear - Perform linear interpolation between the two adjacent bins.
*      .
*      Other (unimplemented) methods give functions which have continuous
*      derivatives, etc.
*
* Bins can be visualized as follows:
* w is the bin width; n is the number of bins; X_0 is the value of the first bin
*  - left histograms go from [ X0 to (X0 + w*n) )
*      X = X0 +            0         w        2w      ...     (n-1)w    n*w
*      bin number          |    1    |    2    |      ...       |    n    |
*  - center histograms go from [ (X0 - w/2) to (X0 + (n-1)w/2)) )
*      X = X0 +          -w/2     (1/2)w    (3/2)w    ...   (n-3/2)w   (n-1/2)w
*      bin number          |    1    |    2    |      ...       |    n    |
*  - right histograms go from ( (X0 - w) to (X0 + (n-1)w) ]
*      X = X0 +           -w         0         w      ...    (n-2)w    (n-1)w
*      bin number          |    1    |    2    |      ...       |    n    |
*
* @tparam X The range of the function. Should support the operations expected
*           of real types. Examples: numeric::Real, float, double
* @tparam Y The domain of the function. Should support the operations expected
*           of real types.
*/
template<typename X, typename Y>
class Histogram : public utility::pointer::ReferenceCount {
public:
	enum BinPlacement {
		left,
		center
		//, right //no one actually needs this yet, but would be easy to implement
	};
	/// @todo It would be cool to implement the different ways of interpolating
	///  using subclasses of Histogram rather than this enum.
	enum Interpolator {
		flat,
		linear,
		spline
	};

	std::string
	to_string( Interpolator const & interpolator ) const{
		if ( interpolator == flat ) return "flat";
		if ( interpolator == linear ) return "linear";
		if ( interpolator == spline ) return "spline";
		return "unrecognized";
	}

protected:
	/**
	* @brief Read a score function from the minirosetta_database into an array
	*
	* @details The scoring function should be represented as a list of Energies,
	* one number per line. Lines begining with '\#' are ignored as comments.
	*
	* Files can contain parameter settings such as the range and step size.
	* These are given by directives beginning with '\@'.
	*
	* @note The database files should be ASCII. Unicode is not supported.
	*/
	static void read_from_db(std::istream & db_file, utility::vector1<Y> /*out*/& energies,
		std::map<std::string, std::string> /*out*/& params)
	{
		using namespace std;

		db_file >> skipws;

		while ( ! db_file.eof() && db_file.good() ) {
			int nextchar = db_file.peek();

			switch (nextchar) {
			case '#' : { // Ignore comments
				string line;
				getline(db_file,line);
				continue;
			}
			case ' ' : //ignore leading whitespace
			case '\t':
			case '\r':
			case '\n' :
				db_file.ignore();
				continue;
			case EOF : //error or end of file
				if ( ! db_file.good() && ! db_file.eof() ) { // error
					cerr << __FILE__ << ":" << __LINE__ << " [ERROR] "
						<< "IO Error" << endl;
				}
				db_file.ignore(); //Make progress to avoid hanging
				break;
			case '@' : { //parameter
				string key, value;
				db_file.ignore(); // '@'
				db_file >> key >> ws;
				getline(db_file, value);
				params[key] = value;
				break;
			}
			default : //Energies
				Y y;
				db_file >> y;
				energies.push_back(y);
				continue;
			}
		}
	}

	/**
	* @brief Set properties of this histogram from a map of strings.
	*
	* Input is validated before being stored. Invalid input results in a
	* printed warning and the previous (probably default) value being used
	* instead.
	*
	* Parameters currently recognised:
	*  - \@minimum <X>
	*  - \@maximum <X>
	*  - \@step <X>
	*  - \@periodic <bool>
	*  - \@bins <BinPlacement>
	*  - \@interpolator <Interpolator>
	*/
	void set_params(std::map<std::string, std::string> const& params)
	{
		using namespace std;

		//the value and whether it was set by the user
		pair<X,bool> min( minimum(), false);
		pair<X,bool> max( maximum(), false);
		pair<X,bool> step( step_, false);

		for (const auto & param : params) {
			string key( param.first);
			string value(param.second);
			//to lowercase
			transform(key.begin(), key.end(), key.begin(), ::tolower );

			if ( "minimum" == key ) {
				istringstream value_strm(value);
				value_strm >> min.first;
				if ( value_strm.fail() ) {
					cerr << __FILE__ << ":" << __LINE__ << " [WARNING] "
						<< "Unrecognized value for @minimum: "
						<< value << endl;
					min.second=false;
				} else {
					min.second=true;
				}
			} else if ( "maximum" == key ) {
				istringstream value_strm(value);
				value_strm >> max.first;
				if ( value_strm.fail() ) {
					cerr << __FILE__ << ":" << __LINE__ << " [WARNING] "
						<< "Unrecognized value for @maximum: "
						<< value << endl;
					max.second = false;
				} else {
					max.second = true;
				}
			} else if ( "step" == key ) {
				istringstream value_strm(value);
				value_strm >> step.first;
				if ( value_strm.fail() ) {
					cerr << __FILE__ << ":" << __LINE__ << " [WARNING] "
						<< "Unrecognized value for @minimum: "
						<< value << endl;
					step.second = false;
				} else {
					step.second = true;
				}
			} else if ( "periodic" == key ) {
				//lowercase
				transform(value.begin(), value.end(), value.begin(), ::tolower );
				if ( "true" == value ) {
					periodic_ = true;
				} else if ( "false" == value ) {
					periodic_ = false;
				} else {
					cerr << __FILE__ << ":" << __LINE__ << " [WARNING] "
						<< "Unrecognized value for @periodic: "
						<< value << endl;
				}
			} else if ( "bins" == key ) {
				//lowercase
				transform(value.begin(), value.end(), value.begin(), ::tolower );
				if ( "left" == value ) {
					bin_placement_ = left;
				} else if ( "center" == value ) {
					bin_placement_ = center;
				} else {
					cerr << __FILE__ << ":" << __LINE__ << " [WARNING] "
						<< "Unrecognized value for @bins: "
						<< value << endl;
				}
			} else if ( "interpolator" == key ) {
				transform(value.begin(), value.end(), value.begin(), ::tolower );
				if ( "flat" == value ) {
					interpolator_ = flat;
				} else if ( "linear" == value ) {
					interpolator_ = linear;
				} else if ( "spline" == value ) {
					interpolator_ = spline;
				} else {
					cerr << __FILE__ << ":" << __LINE__ << " [WARNING] "
						<< "Unrecognized value for @interpolator: "
						<< value << endl;
				}
			} else {
				cerr << __FILE__ << ":" << __LINE__ << " [WARNING] "
					<< "Ignoring unrecognized parameter @" << key << endl;
			}
		}

		//Done parsing parameters.

		// Set step first
		if ( step.second ) {
			step_ = step.first;
		}
		// Next min
		if ( min.second ) {
			switch (bin_placement_) {
			case left :
				min_ = min.first;
				break;
			case center :
				min_ = min.first + step_/2.0;
				break;
			default :
				utility_exit_with_message("Internal Error: Unrecognized BinPlacement");
			}
		}
		// Finally, validate range against max
		if ( max.second &&
				! eq_tol(max.first, last_bin_right(),
				numeric::NumericTraits<X>::tolerance()*1000,
				numeric::NumericTraits<X>::tolerance()*1000  ) ) {
			cerr << __FILE__ << ":" << __LINE__ << " [WARNING] "
				<< "Range missmatch. Expected range of ["
				<< min_ << ", " << max.first
				<< ") but densities cover ["
				<< min_ << ", " << last_bin_right()
				<< ")." << endl;
		}

		set_interpolator( interpolator_ );

	}

public:
	/**
	* @brief Initialize a histogram with the given density distribution.
	* @param densities A vector giving the densities of each bin
	* @param first_bin The x-value of the first bin
	* @param step_size The width of each bin
	* @param bin_placement Indicate what x-value the bins are mapped to:
	*        the left corner, the center of the bin, or the right corner
	*/
	inline Histogram(utility::vector1<Y> const& densities,
		const X first_bin,
		const X step_size,
		const bool periodic=false,
		const BinPlacement bin_placement=left,
		const Interpolator interp=linear) :
		densities_(densities),
		min_(first_bin),
		step_(step_size),
		periodic_(periodic),
		bin_placement_(bin_placement),
		interpolator_(interp)
	{ }

	/// @brief Copy Constructor
	inline Histogram( Histogram const& h) :
		ReferenceCount(h),
		densities_(h.densities_),
		min_(h.min_),
		step_(h.step_),
		periodic_(h.periodic_),
		bin_placement_(h.bin_placement_),
		interpolator_( h.interpolator_ )
	{ }

	/**
	* @brief Generate Histogram from a file.
	*
	* The parameters for the histogram (eg range, step size, etc) are read from
	* any @param fields in the file header present, otherwise they are set to
	* default values and can be changed after instantiation.
	*
	* @note See Histogram::read_from_db() for more information about the file format.
	*/
	inline Histogram(std::istream & file) :
		densities_(),
		min_(0.0),
		step_(1.0),
		periodic_(false),
		bin_placement_(left),
		interpolator_(linear)
	{
		using namespace std;

		map<string,string> params;
		read_from_db(file,densities_,params);
		set_params(params);
	}

	/// @brief destructor
	inline ~Histogram() override = default;


	/// @brief The densities array.
	inline utility::vector1<Y> densities() const { return densities_; }
	inline utility::vector1<Y> & densities() { return densities_; }

	/// @brief The x-value of the left corner of the first bin
	inline X first_bin() const { return min_; }
	inline X & first_bin() { return min_; }

	/// @brief The x-value of the left corner of the last bin
	inline X last_bin() const { return X(min_ + step_*(nbins()-1) ); }

	/// @brief The x-value of the right corner of the last bin
	inline X last_bin_right() const { return X(min_ + step_*nbins() ); }

	/// @brief Return the distance between two bins
	inline X step_size() const { return step_; }
	inline X & step_size() { return step_; }

	/// @brief Return whether this histogram is periodic
	inline bool periodic() const { return periodic_; }
	inline bool & periodic() { return periodic_; }

	/// @brief The bin placement.
	inline BinPlacement bin_placement() const { return bin_placement_; }
	inline BinPlacement & bin_placement() { return bin_placement_; }

	inline Interpolator interpolator() const { return interpolator_; }
	inline Interpolator & interpolator() { return interpolator_; }

	void set_interpolator( Interpolator interpolator ) {
		interpolator_ = interpolator;
		// create spline interpolator
		if ( interpolator_ == spline ) {
			Real lx  = minimum();
			Real ly  = densities_[1];
			Real ldy = (densities_[2]-densities_[1])/step_;
			Real ux  = maximum();
			Real uy  = densities_[nbins()];
			Real udy = 0.0;
			numeric::interpolation::spline::SplineGenerator gen( lx, ly, ldy, ux, uy, udy );
			// add values skipping minimum and maximum
			for ( Size i = 2; i < densities_.size(); ++i ) {
				Real modx = minimum() + (step_*(i-1));
				Real mody = densities_[i];
				gen.add_known_value( modx, mody );
			}
			spline_interpolator_ = gen.get_interpolator();
		}
	}

	/// @brief The smallest value for which we can interpolate
	/// @details All values of x where minimum()<=x<maximum() can be interpolated.
	inline X minimum() const {
		switch( interpolator_ ) {
		case flat :
			return min_;
		case linear :
			switch (bin_placement_) {
			case left :
				return min_;
			case center :
				return X(min_ + step_*0.5);
			default :
				utility_exit_with_message("Internal Error: Unrecognized BinPlacement");
				return X(-1.);
			}
		case spline :
			return min_;
		default :
			utility_exit_with_message("Internal Error: Unrecognized interpolation method: "+to_string(interpolator_));
			return X(-1);
		}

	}

	/// @brief The largest value for which we can interpolate.
	/// @details All values of x where minimum()<=x<maximum() can be interpolated.
	inline X maximum() const {
		switch( interpolator_ ) {
		case flat :
			return min_ + step_*nbins();
		case linear :
			switch (bin_placement_) {
			case left :
				return X(min_ + step_*(nbins()-1.0) );
			case center :
				return X(min_ + step_*(nbins()-0.5) );
			default :
				utility_exit_with_message("Internal Error: Unrecognized BinPlacement");
				return X(-1.);
			}
		case spline :
			return X(min_ + step_*(nbins()-1.0) );
		default :
			utility_exit_with_message("Internal Error: Unrecognized interpolation method: "+to_string(interpolator_));
			return X(-1);
		}

	}

	/// @brief The number of bins
	inline size_type nbins() const{
		return densities_.size();
	}


	/**
	* @brief Interpolates a density for a given x-value from the histogram
	* @details Takes the periodicity and bin placement into account.
	* @param[in]  x The independant axis value to be interpolated
	* @param[out] y An approximation of f(x), as specified by the Interpolator
	* @return Whether the interpolated value was within the bounds or not.
	*         Periodic functions always return true.
	*/
	inline bool interpolate(X const& x, Y & y) const {
		switch(interpolator_) {
		case flat :
			return interpolate_flat(x, y);
		case linear :
			return interpolate_linear(x, y);
		case spline :
			Real dy;
			return interpolate_spline(x, y, dy);
		default :
			utility_exit_with_message("Internal Error: Unrecognized interpolation method: "+to_string(interpolator_));
			return false;
		}
	}

	/**
	* @brief Interpolates a density for a given x-value from the histogram
	* @details Takes the periodicity and bin placement into account.
	* @param[in]  x The independant axis value to be interpolated
	* @param[out] y An approximation of f(x), as specified by the Interpolator
	* @param[out] dy derivative of f(x). Note: only with spline interpolator for now
	* @return Whether the interpolated value was within the bounds or not.
	*         Periodic functions always return true.
	*/
	inline bool interpolate(X const& x, Y & y, Real & dy) const {
		switch(interpolator_) {
		case spline :
			return interpolate_spline(x, y, dy);
		default :
			utility_exit_with_message("Internal Error: Unrecognized interpolation method: "+to_string(interpolator_));
			return false;
		}
	}


	/**
	* @brief The derivative of f(x), linearly interpolated.
	* @details For x between bins, this is just the slope of the interpolation
	*  line. For x on a bin the slope of the line to the right is used.
	* @param[in]  x  The point on the independant axis for which to get the derivative
	* @param[out] dy An approximation of df/dx, cast to a Y.
	*/
	inline bool derivative(X const& x, Y & dy) const {
		switch(interpolator_) {
		case flat : // Technically 0, but we'll just linearly interpolate
		case linear :
			return derivative_linear(x, dy);
		case spline : {
			Y y;
			Real dY;
			bool retval = interpolate_spline(x, y, dY);
			dy = static_cast< Y >(dY);
			return retval;
		}
		default :
			utility_exit_with_message("Internal Error: Unrecognized interpolation method: "+to_string(interpolator_));
			return false;
		}
	}

	inline bool interpolate_spline(X const& x, Y &y, Real &dy) const {
		spline_interpolator_->interpolate( x, y, dy );
		return true;
	}

protected: // Interpolation methods

	/**
	* @brief get the number of the bin to the left of X.
	*
	* @details Takes periodicity and bin alignment into account.
	* A periodic histogram will always return a number in [1,nbins]
	* A nonperiodic histogram makes no guarentees that its return value will
	* fall within the allowed bounds of 1 through nbins.
	* <i>You should do bounds checking elsewhere to assert that x is within the
	* allowed range.</i>
	*
	* @param x[in]  The independent axis value
	* @param a[out] The alpha fraction: (x-x_l)/(x_u-x_l) for bin [x_l,x_u]
	*
	* @precondition x is in the domain of the histogram. For nonperiodic histograms,
	* this means minimum() <= x < maximum()
	*
	* @return The index of bin x_l
	*/
	inline platform::SSize bin_number(X const& x, X & a) const{
		X const x_normalized(numeric::modulo( (x-first_bin())/step_, X(nbins()) )); //Real [0, nbins)

		const platform::SSize bin( static_cast<platform::SSize>( std::floor(x_normalized) )); //int [0,nbins-1]
		a = x_normalized - bin;
		return bin+1;
	}

	/**
	* @brief Returns the density of the bin which x belongs to
	* @details If x is outside of the range of bins, returns zero.
	*/
	inline bool interpolate_flat(X const& x, Y &y) const {
		//check bounds
		if ( !periodic_ ) {
			if ( minimum() > x ) { //too small; take the minimum
				y = densities_[1];
				return false;
			}
			if ( x >= maximum() ) { //too big; take the maximum
				y = densities_[nbins()];
				return false;
			}
		}

		X alpha; //ignored
		platform::SSize bin = bin_number(x, alpha);
		y = densities_[bin];
		//check bounds
		return true;
	}

	inline bool interpolate_linear(X const& x, Y &y) const {
		//check bounds
		if ( !periodic_ ) {
			if ( minimum() > x ) { //too small; take the minimum
				y = densities_[1];
				return false;
			}
			if ( x >= maximum() ) { //too big; take the maximum
				y = densities_[nbins()];
				return false;
			}
		}

		X alpha(0); //(x-x_l)/(x_u-x_l)
		size_type lower(0), upper(0);
		switch( bin_placement_ ) {
		case left :
			lower = static_cast<size_type>(bin_number(x, alpha));
			break;
		case center :
			lower = static_cast<size_type>(bin_number(x-step_*0.5, alpha));
			break;
		default :
			utility_exit_with_message("Internal Error: Unrecognized interpolation method: "+to_string(interpolator_));
		}
		// lower is [1,nbins]
		upper = ( lower == nbins())?1:lower+1; //wrap around periodic values

		y = numeric::interpolation::interpolated(alpha,densities_[lower], densities_[upper] );

		return true;
	}

	/**
	* @brief The derivative of f(x), linearly interpolated.
	* @details For x between bins, this is just the slope of the interpolation
	*  line. For x on a bin the slope of the line to the right is used.
	*
	*  Note that the derivative will not be continuous when calculated in this way.
	*/
	inline bool derivative_linear(X const& x, Y & y) const {
		//check bounds
		if ( !periodic_ ) {
			if ( minimum() > x || x >= maximum() ) { //too big; take the maximum
				y = Y(0); //Consistant with interpolate's return value here
				return false;
			}
		}

		X alpha(0); //(x-x_l)/(x_u-x_l)
		size_type lower(0), upper(0);
		switch( bin_placement_ ) {
		case left :
			lower = static_cast<size_type>(bin_number(x, alpha));
			break;
		case center :
			lower = static_cast<size_type>(bin_number(x-step_*0.5, alpha));
			break;
		default :
			utility_exit_with_message("Internal Error: Unrecognized interpolation method: "+to_string(interpolator_));
		}
		// lower is [1,nbins]
		upper = ( lower == nbins())?1:lower+1; //wrap around periodic values

		y = Y((densities_[upper]-densities_[lower])/step_);
		return true;
	}

protected:
	utility::vector1<Y> densities_;
	/// @brief the x value of densities_[0]. Not actually the minimum for BinPlacement other than left.
	X min_;
	X step_;
	bool periodic_;
	BinPlacement bin_placement_;
	Interpolator interpolator_;
	utility::pointer::shared_ptr< numeric::interpolation::spline::Interpolator > spline_interpolator_;
#ifdef SERIALIZATION
public:
	Histogram() {}
#endif

}; //Histogram

} //interpolation
} //numeric
#endif //INCLUDED_numeric_Histogram_HH
