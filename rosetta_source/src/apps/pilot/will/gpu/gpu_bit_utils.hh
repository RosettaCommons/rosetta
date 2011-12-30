#ifndef INCLUDED_apps_pilot_will_gpu_gpu_bit_utils_hh
#define INCLUDED_apps_pilot_will_gpu_gpu_bit_utils_hh

#include <core/id/AtomID.hh>

template<typename T>
void printbits32(T const & x) {
  if( sizeof(T) != 4 ) utility_exit_with_message("must be 32 bits!");
  uint i = *((uint*)( &x ));
  std::cout<<" 30 27 24 21 18 15 12 9876543210"<<" as uint: "<<i<<std::endl;
  std::cout<<" |  |  |  |  |  |  |  ||||||||||"<<" as uint: "<<i<<std::endl;
  for( int ib = 31; ib >= 0; --ib) {
    std::cout<<(i>>ib)%(2u);
  }
  std::cout<<" as uint: "<<i<<std::endl;
}

template<typename T>
void printbits16(T const & x) {
  if( sizeof(T) != 2 ) utility_exit_with_message("must be 32 bits!");
  ushort i = *((uint*)( &x ));
  std::cout<<"15 12 9876543210"<<" as ushort: "<<i<<std::endl;
  std::cout<<"|  |  ||||||||||"<<" as ushort: "<<i<<std::endl;
  for( int ib = 15; ib >= 0; --ib) {
    std::cout<<(i>>ib)%(2u);
  }
  std::cout<<" as ushort: "<<i<<std::endl;
}

inline float bitsasfloat(uint const & i) {
  return *((float*)(&i));
}

float aidr_as_float(core::id::AtomID const & aid, core::Real const & radius) {
  uint i = std::floor(radius*100.0);
  i += ((uint)aid.atomno()<< 8u);
  i += ((uint)aid.rsd()   <<16u);
  return *((float*)(&i));
}

inline float float_as_aidr(float const & aidr, core::id::AtomID & aid) {
  uint i = *((uint*)(&aidr));
  float r = ((float)(254u & i))/100.0f;
  aid.atomno() =   255u & (i >>  8u);
  aid.rsd()    = 65535u & (i >> 16u);
  return r;
}

#endif
