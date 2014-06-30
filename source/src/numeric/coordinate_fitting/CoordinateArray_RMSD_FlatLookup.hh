#include <vector>

#include "numeric/coordinate_fitting/FlatLookup.hh"
#include "numeric/alignment/QCP_Kernel.hh"

namespace numeric
{
namespace coordinate_fitting
{

template <class Real=double>
class CoordinateArray_RMSD_FlatLookup : public FlatLookup<Real *, std::size_t, Real>
{
  public:
    CoordinateArray_RMSD_FlatLookup(Real* entry_coordinates, Real* entry_radii, std::size_t num_entries, std::size_t coordinates_per_entry) :
      FlatLookup<Real *, std::size_t, Real>(),
      entry_coordinates(entry_coordinates),
      entry_radii(entry_radii),
      num_entries(num_entries),
      coordinates_per_entry(coordinates_per_entry),
      entry_size(3 * coordinates_per_entry),
      kernel()
  {
    std::vector<std::size_t> entry_indicies(num_entries);

    for (std::size_t i = 0; i < num_entries; i++)
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

  virtual Real entry_distance(Real *& q, std::size_t & e)
  {
    return kernel.calc_centered_coordinate_rmsd(
        q,
        &entry_coordinates[e * entry_size],
        coordinates_per_entry,
        NULL);
  }

  virtual Real entry_radius(std::size_t & e)
  {
    return entry_radii[e];
  }

  Real* entry_coordinates;
  Real* entry_radii;

  std::size_t num_entries, coordinates_per_entry, entry_size;

  numeric::alignment::QCP_Kernel<Real> kernel;

};

}
}
