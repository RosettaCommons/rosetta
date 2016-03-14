
__kernel void test_input_struct(__global struct POSE *p,
                                __global float  *tor,
                                __global float4 *chi,
                                __global int    *nres
) {
  nres[0] = p->nres;
}
