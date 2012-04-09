/*
 * OrbitalXYZCoords.hh
 *
 *  Created on: Jun 30, 2011
 *      Author: combss
 */

#ifndef ORBITALXYZCOORDS_HH_
#define ORBITALXYZCOORDS_HH_

// Project headers
#include <core/types.hh>

// Utility headers
#include <numeric/xyzVector.hh>


namespace core{
namespace conformation{
namespace orbitals{


class OrbitalXYZCoords {
public:
	/// @brief default constructor and set atom type number to 0 and place the
	/// atom at the origin
	OrbitalXYZCoords():
		xyz_( 0.0 )
	{}

	/// @brief constructor with an atom type number
	// type is set at construction time -- atom is placed at the origin
	OrbitalXYZCoords( ShortSize const, ShortSize const ):
		xyz_( 0.0 )
	{}

	/// @brief constructor with xyz and an atom type number
	// Atom( Vector const & xyz_in, int const type_in, Real temperature = 0.0 ):
	OrbitalXYZCoords( Vector const & xyz_in ):
		xyz_( xyz_in )
	{}

	/// @brief destructor
	virtual
	~OrbitalXYZCoords() {}



	Vector const &
	xyz() const
	{
		return xyz_;
	}


	void
	xyz( Vector const & xyz_in )
	{
		xyz_ = xyz_in;
	}


private:
	/// xyz coordinates
	Vector xyz_;




};




}
}
}




#endif /* ORBITALXYZCOORDS_HH_ */
