#ifndef FEM_UTILS_EQUIVALENCE_HPP
#define FEM_UTILS_EQUIVALENCE_HPP

#include <fem/size_t.hpp>
#include <tbxx/error_utils.hpp>
#include <algorithm>
#include <vector>

namespace fem { namespace utils { namespace equivalence {

  struct array_alignment
  {
    size_t members_size;
    std::vector<ssize_t> diff_matrix;
    std::vector<ssize_t> diffs0;

    array_alignment(
      size_t members_size_)
    :
      members_size(members_size_),
      diff_matrix(members_size * (members_size-1), ssize_t_max)
    {}

    static
    std::string
    msg_prefix() { return "equivalence::array_alignment: "; }

    void
    add_anchor(
      size_t i0, ssize_t a0,
      size_t i1, ssize_t a1)
    {
      static const char* msg_directly = "directly conflicting input";
      if (i0 == i1) {
        if (a0 != a1) {
          throw std::runtime_error(msg_prefix() + msg_directly);
        }
      }
      else {
        size_t n = members_size;
        size_t i;
        ssize_t d;
        if (i0 < i1) {
          i = i0 * n + i1;
          d = a0 - a1;
        }
        else {
          i = i1 * n + i0;
          d = a1 - a0;
        }
        ssize_t dd = diff_matrix[i];
        if (dd == ssize_t_max) {
          diff_matrix[i] = d;
        }
        else if (dd != d) {
          throw std::runtime_error(msg_prefix() + msg_directly);
        }
      }
    }

    void
    infer_diffs0_from_diff_matrix()
    {
      size_t n = members_size;
      std::vector<size_t> cluster_indices(n);
      for(size_t i=0;i<n;i++) cluster_indices[i] = i;
      std::vector<std::vector<ssize_t_2> > clusters(n);
      for(size_t li0=0;li0<n-1;li0++) {
        for(size_t li1=li0+1;li1<n;li1++) {
          size_t i0 = li0;
          size_t i1 = li1;
          size_t i = i0 * n + i1;
          ssize_t d = diff_matrix[i];
          if (d != ssize_t_max) {
            size_t ci0 = cluster_indices[i0];
            size_t ci1 = cluster_indices[i1];
            if (ci0 == ci1) {
              continue;
            }
            if (ci0 > ci1) {
              std::swap(ci0, ci1);
              d *= -1;
              std::swap(i0, i1);
            }
            if (ci1 != i1) {
              continue;
            }
            std::vector<ssize_t_2>& c0 = clusters[ci0];
            std::vector<ssize_t_2>& c1 = clusters[ci1];
            if (ci0 != i0) {
              for(size_t i=0;i!=c0.size();i++) {
                ssize_t const* c0e = c0[i].elems;
                if (c0e[0] == i0) {
                  d += c0e[1];
                  break;
                }
              }
            }
            c0.push_back(ssize_t_2(ci1, d));
            cluster_indices[ci1] = ci0;
            for(size_t i=0;i!=c1.size();i++) {
              ssize_t const* c1e = c1[i].elems;
              c0.push_back(ssize_t_2(c1e[0], c1e[1]+d));
              cluster_indices[c1e[0]] = ci0;
            }
            clusters[ci1].clear();
          }
        }
      }
      for(size_t i=0;i<n;i++) {
        if (cluster_indices[i] != 0) {
          throw std::runtime_error(
            msg_prefix() + "insufficient input");
        }
      }
      std::vector<ssize_t_2>& c0 = clusters[0];
      TBXX_ASSERT(c0.size() == n-1);
      diffs0.clear();
      diffs0.resize(n, ssize_t_max);
      diffs0[0] = 0;
      for(size_t i=0;i!=c0.size();i++) {
        ssize_t const* c0e = c0[i].elems;
        TBXX_ASSERT(c0e[0] != 0);
        TBXX_ASSERT(diffs0[c0e[0]] == ssize_t_max);
        diffs0[c0e[0]] = c0e[1];
      }
      for(size_t i0=0;i0<n-1;i0++) {
        for(size_t i1=i0+1;i1<n;i1++) {
          size_t i = i0 * n + i1;
          ssize_t d = diff_matrix[i];
          if (   d != ssize_t_max
              && diffs0[i1] - diffs0[i0] != d) {
            throw std::runtime_error(
              msg_prefix() + "indirectly conflicting input");
          }
        }
      }
    }
  };

}}} // namespace fem::utils::equivalence

#endif // GUARD
