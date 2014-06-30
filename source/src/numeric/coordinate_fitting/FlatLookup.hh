#include <vector>
#include <limits>
#include <cstddef>

namespace numeric
{
namespace coordinate_fitting
{

template <class QueryType, class EntryType, class Real=double>
class FlatLookup
{
  public:
    ~FlatLookup() {};

    template <class InputIterator>
    void initialize(InputIterator first, InputIterator last)
    {
      entries.assign(first, last);
    }

    // @brief Find first matching object.
    bool first_match(QueryType & query, EntryType & entry, Real & distance)
    {
      prepare_for_query(query);

      distance = std::numeric_limits<Real>::max();

      for (std::size_t i = 0; i < entries.size(); i++)
      {
        Real test_distance = entry_distance(query, entries[i]);

        if ( test_distance < entry_radius(entries[i]))
        {
          distance = test_distance;
          entry = entries[i];
          return true;
        }
      }

      return false;
    }

    // @brief Find closest matching object.
    bool closest_match(QueryType & query, EntryType & entry, Real & distance)
    {
      prepare_for_query(query);

      distance = std::numeric_limits<Real>::max();

      for (std::size_t i = 0; i < entries.size(); i++)
      {
        Real test_distance = entry_distance(query, entries[i]);

        if ( (test_distance < distance) && (test_distance < entry_radius(entries[i])))
        {
          distance = test_distance;
          entry = entries[i];
        }
      }

      if (distance < std::numeric_limits<Real>::max())
      {
        return true;
      }
      else
      {
        return false;
      }
    }

    // @brief Called before query evaluation.
    virtual void prepare_for_query(QueryType &) { }

    // @brief Calculate distance between query object and entry
    virtual Real entry_distance(QueryType & q, EntryType & e2) = 0;

    // @brief Get maximum radius from entry
    virtual Real entry_radius(EntryType & e) = 0;

    // @brief Souce entries query
    std::vector<EntryType>  entries;
};

}
}
