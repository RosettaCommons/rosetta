#include <vector>

#include "numeric/coordinate_fitting/FlatLookup.hh"
#include "numeric/alignment/QCP_Kernel.hh"

#ifndef COORDINATEARRRAY_RMSD_FLATLOOKUP_HH_
#define COORDINATEARRRAY_RMSD_FLATLOOKUP_HH_

#include <numeric/types.hh>

namespace numeric
{
namespace coordinate_fitting
{

template <class Real=double>
class CoordinateArray_RMSD_FlatLookup : public FlatLookup<Real *, numeric::Size, Real>
{
  public:
    CoordinateArray_RMSD_FlatLookup(Real* entry_coordinates, Real* entry_radii, numeric::Size num_entries, numeric::Size coordinates_per_entry) :
      FlatLookup<Real *, numeric::Size, Real>(),
      entry_coordinates(entry_coordinates),
      entry_radii(entry_radii),
      num_entries(num_entries),
      coordinates_per_entry(coordinates_per_entry),
      entry_size(3 * coordinates_per_entry),
      kernel()
  {
    std::vector<numeric::Size> entry_indicies(num_entries);

    for (numeric::Size i = 0; i < num_entries; i++)
    {
      entry_indicies[i] = i;
      kernel.remove_center_of_mass(&entry_coordinates[i * entry_size], coordinates_per_entry) ;
    }

    this->initialize(entry_indicies.begin(), entry_indicies.end());
  }

  virtual void prepare_for_query(Real *& q)
  {
    kernel.remove_center_of_mass(q, coordinates_per_entry);
  }

  virtual Real entry_distance(Real *& q, numeric::Size & e)
  {
    return kernel.calc_centered_coordinate_rmsd(
        q,
        &entry_coordinates[e * entry_size],
        coordinates_per_entry,
        NULL);
  }

  virtual Real entry_radius(numeric::Size & e)
  {
    return entry_radii[e];
  }

  Real* entry_coordinates;
  Real* entry_radii;

  numeric::Size num_entries, coordinates_per_entry, entry_size;

  numeric::alignment::QCP_Kernel<Real> kernel;

};

}
}
#endif
