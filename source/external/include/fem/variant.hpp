#ifndef FEM_VARIANT_HPP
#define FEM_VARIANT_HPP

#include <fem/str_arr_ref.hpp>
#include <fem/utils/equivalence.hpp>
#include <fem/utils/misc.hpp>
#include <tbxx/optional_copy.hpp>

namespace fem {

  using tbxx::optional_copy;

  struct variant_member
  {
    size_t type_size;
    optional_copy<dim_buffer> dims;

    variant_member() {}

    variant_member(
      size_t type_size_)
    :
      type_size(type_size_)
    {}

    template <typename DimsType>
    variant_member(
      size_t type_size_,
      DimsType const& dims_)
    :
      type_size(type_size_),
      dims(dims_)
    {}

    size_t
    number_of_bytes() const
    {
      if (!dims) return type_size;
      return type_size * dims->actual_size_1d();
    }

    size_t
    actual_index_1d(
      arr_index const& arr_ix)
    {
      TBXX_ASSERT(dims.get() != 0);
      return dims->actual_index_1d(arr_ix);
    }
  };

  template <typename T>
  struct mbr : variant_member
  {
    mbr() : variant_member(sizeof(T)) {}

    template <typename DimsType>
    mbr(
      DimsType const& dims_) : variant_member(sizeof(T), dims_) {}
  };

  template <int StrLen>
  struct mbr<str<StrLen> > : variant_member
  {
    mbr() : variant_member(StrLen) {}

    template <typename DimsType>
    mbr(
      DimsType const& dims_) : variant_member(StrLen, dims_) {}
  };

  struct equivalence
  {
    static const size_t buffer_capacity = 7;
    size_t members_size;
    variant_member members[buffer_capacity];
    utils::equivalence::array_alignment array_alignment;
    size_t align_i_mbr;
    size_t align_byte_offset;

    equivalence(
      variant_member const& mbr_1,
      variant_member const& mbr_2)
    :
      members_size(2),
      array_alignment(members_size),
      align_i_mbr(size_t_max),
      align_byte_offset(size_t_max)
    {
      members[0] = mbr_1;
      members[1] = mbr_2;
    }

    equivalence(
      variant_member const& mbr_1,
      variant_member const& mbr_2,
      variant_member const& mbr_3)
    :
      members_size(3),
      array_alignment(members_size),
      align_i_mbr(size_t_max),
      align_byte_offset(size_t_max)
    {
      members[0] = mbr_1;
      members[1] = mbr_2;
      members[2] = mbr_3;
    }

    equivalence(
      variant_member const& mbr_1,
      variant_member const& mbr_2,
      variant_member const& mbr_3,
      variant_member const& mbr_4)
    :
      members_size(4),
      array_alignment(members_size),
      align_i_mbr(size_t_max),
      align_byte_offset(size_t_max)
    {
      members[0] = mbr_1;
      members[1] = mbr_2;
      members[2] = mbr_3;
      members[3] = mbr_4;
    }

    equivalence(
      variant_member const& mbr_1,
      variant_member const& mbr_2,
      variant_member const& mbr_3,
      variant_member const& mbr_4,
      variant_member const& mbr_5)
    :
      members_size(5),
      array_alignment(members_size),
      align_i_mbr(size_t_max),
      align_byte_offset(size_t_max)
    {
      members[0] = mbr_1;
      members[1] = mbr_2;
      members[2] = mbr_3;
      members[3] = mbr_4;
      members[4] = mbr_5;
    }

    equivalence(
      variant_member const& mbr_1,
      variant_member const& mbr_2,
      variant_member const& mbr_3,
      variant_member const& mbr_4,
      variant_member const& mbr_5,
      variant_member const& mbr_6)
    :
      members_size(6),
      array_alignment(members_size),
      align_i_mbr(size_t_max),
      align_byte_offset(size_t_max)
    {
      members[0] = mbr_1;
      members[1] = mbr_2;
      members[2] = mbr_3;
      members[3] = mbr_4;
      members[4] = mbr_5;
      members[5] = mbr_6;
    }

    equivalence(
      variant_member const& mbr_1,
      variant_member const& mbr_2,
      variant_member const& mbr_3,
      variant_member const& mbr_4,
      variant_member const& mbr_5,
      variant_member const& mbr_6,
      variant_member const& mbr_7)
    :
      members_size(7),
      array_alignment(members_size),
      align_i_mbr(size_t_max),
      align_byte_offset(size_t_max)
    {
      members[0] = mbr_1;
      members[1] = mbr_2;
      members[2] = mbr_3;
      members[3] = mbr_4;
      members[4] = mbr_5;
      members[5] = mbr_6;
      members[6] = mbr_7;
    }

    template <size_t Jmbr>
    equivalence&
    align(
      size_t index_1d=0,
      size_t str_offset=0)
    {
      align_i_mbr = Jmbr-1;
      align_byte_offset = members[align_i_mbr].type_size * index_1d
                        + str_offset;
      return *this;
    }

    template <size_t Jmbr>
    equivalence&
    align(
      arr_index const& arr_ix)
    {
      align<Jmbr>(members[Jmbr-1].actual_index_1d(arr_ix));
      return *this;
    }

    template <size_t Jmbr>
    equivalence&
    align(
      arr_and_str_indices<arr_dim_max> const& indices)
    {
      align<Jmbr>(
        members[Jmbr-1].actual_index_1d(indices.arr_ix),
        static_cast<size_t>(indices.str_ix.first - 1));
      return *this;
    }

    template <size_t Jmbr>
    equivalence&
    align(
      str_index const& ix)
    {
      align<Jmbr>(0, static_cast<size_t>(ix.first - 1));
      return *this;
    }

    template <size_t Jmbr>
    equivalence&
    with(
      size_t index_1d=0,
      size_t str_offset=0)
    {
      array_alignment.add_anchor(
        align_i_mbr, align_byte_offset,
        Jmbr-1, members[Jmbr-1].type_size * index_1d + str_offset);
      return *this;
    }

    template <size_t Jmbr>
    equivalence&
    with(
      arr_index const& arr_ix)
    {
      with<Jmbr>(members[Jmbr-1].actual_index_1d(arr_ix));
      return *this;
    }

    template <size_t Jmbr>
    equivalence&
    with(
      arr_and_str_indices<arr_dim_max> const& indices)
    {
      with<Jmbr>(
        members[Jmbr-1].actual_index_1d(indices.arr_ix),
        static_cast<size_t>(indices.str_ix.first - 1));
      return *this;
    }

    template <size_t Jmbr>
    equivalence&
    with(
      str_index const& ix)
    {
      with<Jmbr>(0, static_cast<size_t>(ix.first - 1));
      return *this;
    }
  };

  struct variant_core : utils::noncopyable
  {
    size_t use_count;
    size_t size;
    char* ptr;

    variant_core() : use_count(0), size(0), ptr(0) {}

    ~variant_core() { delete[] ptr; }

    void
    grow_if_necessary(
      size_t new_size)
    {
      if (size < new_size) {
        TBXX_ASSERT(use_count == 1);
        char* new_ptr = new char[new_size];
        std::memcpy(new_ptr, ptr, size);
        delete[] ptr;
        std::memset(new_ptr+size, 0, new_size-size);
        size = new_size;
        ptr = new_ptr;
      }
    }
  };

  struct variant_bind_info
  {
    size_t offset;
    size_t type_size;
  };

  typedef std::vector<variant_bind_info> variant_bindings;

  struct variant_core_and_bindings
  {
    variant_core core;
    variant_bindings bindings;
  };

  struct variant_allocate_chain
  {
    variant_core* core;
    variant_bindings* bindings;
    bool is_common_variant;
    size_t curr_bytes;
    size_t end_bytes;

    variant_allocate_chain(
      variant_core* core_,
      variant_bindings* bindings_,
      bool is_common_variant_)
    :
      core(core_),
      bindings(bindings_),
      is_common_variant(is_common_variant_),
      curr_bytes(0),
      end_bytes(0)
    {}

    ~variant_allocate_chain()
    {
      core->grow_if_necessary(end_bytes);
    }

    variant_allocate_chain&
    operator,(
      variant_member const& m)
    {
      TBXX_ASSERT(is_common_variant);
      bindings->resize(bindings->size() + 1);
      variant_bind_info& bi = bindings->back();
      bi.offset = curr_bytes;
      bi.type_size = m.type_size;
      size_t nb = m.number_of_bytes();
      curr_bytes += nb;
      end_bytes = std::max(end_bytes, curr_bytes);
      return *this;
    }

    variant_allocate_chain&
    operator,(
      equivalence& e)
    {
      e.array_alignment.infer_diffs0_from_diff_matrix();
      size_t prev_bindings_size = bindings->size();
      bindings->reserve(prev_bindings_size + e.members_size);
      size_t origin_bytes = 0;
      if (!is_common_variant) {
        for(size_t i_mbr=0;i_mbr<e.members_size;i_mbr++) {
          ssize_t diff0 = e.array_alignment.diffs0[i_mbr];
          if (diff0 < 0) {
            origin_bytes = std::max(origin_bytes, static_cast<size_t>(-diff0));
          }
        }
      }
      for(size_t i_mbr=0;i_mbr<e.members_size;i_mbr++) {
        variant_member const& m = e.members[i_mbr];
        bindings->resize(bindings->size() + 1);
        variant_bind_info& bi = bindings->back();
        ssize_t diff0 = e.array_alignment.diffs0[i_mbr];
        if (   is_common_variant
            && diff0 < 0
            && static_cast<size_t>(-diff0) > curr_bytes) {
          throw std::runtime_error(
            "EQUIVALENCE crosses beginning of COMMON block");
        }
        bi.offset = static_cast<size_t>(
          static_cast<ssize_t>(curr_bytes + origin_bytes) + diff0);
        bi.type_size = m.type_size;
        end_bytes = std::max(end_bytes, bi.offset + m.number_of_bytes());
      }
      if (is_common_variant) {
        curr_bytes += e.members[0].number_of_bytes();
      }
      else {
        curr_bytes = end_bytes;
      }
      return *this;
    }
  };

  struct variant_bind_chain
  {
    variant_core* core;
    variant_bindings* bindings;
    bool is_common_variant;
    size_t bind_index;

    variant_bind_chain(
      variant_core* core_,
      variant_bindings* bindings_,
      bool is_common_variant_)
    :
      core(core_),
      bindings(bindings_),
      is_common_variant(is_common_variant_),
      bind_index(0)
    {
      core->use_count++;
    }

    ~variant_bind_chain()
    {
      core->use_count--;
    }

    variant_allocate_chain
    allocate()
    {
      return variant_allocate_chain(core, bindings, is_common_variant);
    }

    template <typename T>
    T&
    bind()
    {
      variant_bind_info& bi = (*bindings)[bind_index++];
      return *reinterpret_cast<T*>(core->ptr + bi.offset);
    }

    str_ref
    bind_str()
    {
      variant_bind_info& bi = (*bindings)[bind_index++];
      return str_ref(core->ptr + bi.offset, bi.type_size);
    }
  };

  struct common_variant : variant_bind_chain
  {
    common_variant(
      variant_core& core,
      variant_bindings& bindings)
    :
      variant_bind_chain(&core, &bindings, /*is_common_variant*/ true)
    {}
  };

  struct save_equivalences : variant_bind_chain
  {
    save_equivalences(
      variant_core_and_bindings& core_and_bindings)
    :
      variant_bind_chain(
        &core_and_bindings.core,
        &core_and_bindings.bindings,
        /*is_common_variant*/ false)
    {}
  };

  struct local_equivalences :
    utils::hide<variant_core_and_bindings>,
    variant_bind_chain
  {
    local_equivalences() :
      variant_bind_chain(
        &hidden.core, &hidden.bindings, /*is_common_variant*/ false)
    {}
  };

} // namespace fem

#endif // GUARD
