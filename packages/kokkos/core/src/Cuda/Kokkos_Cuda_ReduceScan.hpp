/*
//@HEADER
// ************************************************************************
// 
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_CUDA_REDUCESCAN_HPP
#define KOKKOS_CUDA_REDUCESCAN_HPP

#if defined( __CUDACC__ )

#include <utility>

#include <Kokkos_Parallel.hpp>
#include <impl/Kokkos_FunctorAdapter.hpp>
#include <impl/Kokkos_Error.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
// See section B.17 of Cuda C Programming Guide Version 3.2
// for discussion of
//   __launch_bounds__(maxThreadsPerBlock,minBlocksPerMultiprocessor)
// function qualifier which could be used to improve performance.
//----------------------------------------------------------------------------
// Maximize shared memory and minimize L1 cache:
//   cudaFuncSetCacheConfig(MyKernel, cudaFuncCachePreferShared );
// For 2.0 capability: 48 KB shared and 16 KB L1
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
/*
 *  Algorithmic constraints:
 *   (a) blockDim.y is a power of two
 *   (b) blockDim.y <= 512
 *   (c) blockDim.x == blockDim.z == 1
 */


template< class ValueType , class JoinOp>
__device__
inline void cuda_intra_warp_reduction( ValueType& value, const JoinOp& join) {

  unsigned shift = 1;
  ValueType result = value;
  while(blockDim.x * shift < 32 && shift<blockDim.y) {
    const ValueType tmp = shfl_down(result, blockDim.x*shift,32);
    join(result , tmp);
    shift*=2;
  }
  value = shfl(result,0,32);
}

template< class ValueType , class JoinOp>
__device__
inline void cuda_inter_warp_reduction( ValueType& value, const JoinOp& join) {

  __shared__ char sh_result[sizeof(ValueType)];
  ValueType* result = (ValueType*) & sh_result;
  const unsigned step = 32 / blockDim.x;
  unsigned shift = step;
  __syncthreads();

  if(threadIdx.y == 0) {
    *result = value;
  }
  __syncthreads();
  while (shift<blockDim.y) {
    if(shift==threadIdx.y && threadIdx.x==0) {
      join(*result,value);
    }
    __syncthreads();
    shift+=step;
  }
  __syncthreads();

  value = *result;
}



template< bool DoScan , class FunctorType >
__device__
void cuda_intra_block_reduce_scan( const FunctorType & functor ,
                                   const typename FunctorValueTraits< FunctorType , void >::pointer_type base_data )
{
  typedef FunctorValueTraits< FunctorType , void >  ValueTraits ;
  typedef FunctorValueJoin<   FunctorType , void >  ValueJoin ;

  typedef typename ValueTraits::pointer_type  pointer_type ;

  const unsigned value_count   = ValueTraits::value_count( functor );
  const unsigned BlockSizeMask = blockDim.y - 1 ;

  // Must have power of two thread count

  if ( BlockSizeMask & blockDim.y ) { cuda_abort("Cuda::cuda_intra_block_scan requires power-of-two blockDim"); }

#define BLOCK_REDUCE_STEP( R , TD , S )  \
  if ( ! ( R & ((1<<(S+1))-1) ) ) { ValueJoin::join( functor , TD , (TD - (value_count<<S)) ); }

#define BLOCK_SCAN_STEP( TD , N , S )  \
  if ( N == (1<<S) ) { ValueJoin::join( functor , TD , (TD - (value_count<<S))); }

  const unsigned     rtid_intra = threadIdx.y ^ BlockSizeMask ;
  const pointer_type tdata_intra = base_data + value_count * threadIdx.y ;

  { // Intra-warp reduction:
    BLOCK_REDUCE_STEP(rtid_intra,tdata_intra,0)
    BLOCK_REDUCE_STEP(rtid_intra,tdata_intra,1)
    BLOCK_REDUCE_STEP(rtid_intra,tdata_intra,2)
    BLOCK_REDUCE_STEP(rtid_intra,tdata_intra,3)
    BLOCK_REDUCE_STEP(rtid_intra,tdata_intra,4)
  }

  __syncthreads(); // Wait for all warps to reduce

  { // Inter-warp reduce-scan by a single warp to avoid extra synchronizations
    const unsigned rtid_inter = ( threadIdx.y ^ BlockSizeMask ) << CudaTraits::WarpIndexShift ;

    if ( rtid_inter < blockDim.y ) {

      const pointer_type tdata_inter = base_data + value_count * ( rtid_inter ^ BlockSizeMask );

      if ( (1<<5) < BlockSizeMask ) {                        BLOCK_REDUCE_STEP(rtid_inter,tdata_inter,5) }
      if ( (1<<6) < BlockSizeMask ) { __threadfence_block(); BLOCK_REDUCE_STEP(rtid_inter,tdata_inter,6) }
      if ( (1<<7) < BlockSizeMask ) { __threadfence_block(); BLOCK_REDUCE_STEP(rtid_inter,tdata_inter,7) }
      if ( (1<<8) < BlockSizeMask ) { __threadfence_block(); BLOCK_REDUCE_STEP(rtid_inter,tdata_inter,8) }

      if ( DoScan ) {

        int n = ( rtid_inter &  32 ) ?  32 : (
                ( rtid_inter &  64 ) ?  64 : (
                ( rtid_inter & 128 ) ? 128 : (
                ( rtid_inter & 256 ) ? 256 : 0 )));

        if ( ! ( rtid_inter + n < blockDim.y ) ) n = 0 ;

        BLOCK_SCAN_STEP(tdata_inter,n,8)
        BLOCK_SCAN_STEP(tdata_inter,n,7)
        BLOCK_SCAN_STEP(tdata_inter,n,6)
        BLOCK_SCAN_STEP(tdata_inter,n,5)
      }
    }
  }

  __syncthreads(); // Wait for inter-warp reduce-scan to complete

  if ( DoScan ) {
    int n = ( rtid_intra &  1 ) ?  1 : (
            ( rtid_intra &  2 ) ?  2 : (
            ( rtid_intra &  4 ) ?  4 : (
            ( rtid_intra &  8 ) ?  8 : (
            ( rtid_intra & 16 ) ? 16 : 0 ))));

    if ( ! ( rtid_intra + n < blockDim.y ) ) n = 0 ;

    BLOCK_SCAN_STEP(tdata_intra,n,4) __threadfence_block();
    BLOCK_SCAN_STEP(tdata_intra,n,3) __threadfence_block();
    BLOCK_SCAN_STEP(tdata_intra,n,2) __threadfence_block();
    BLOCK_SCAN_STEP(tdata_intra,n,1) __threadfence_block();
    BLOCK_SCAN_STEP(tdata_intra,n,0)
  }

#undef BLOCK_SCAN_STEP
#undef BLOCK_REDUCE_STEP
}

//----------------------------------------------------------------------------
/**\brief  Input value-per-thread starting at 'shared_data'.
 *         Reduction value at last thread's location.
 *
 *  If 'DoScan' then write blocks' scan values and block-groups' scan values.
 *
 *  Global reduce result is in the last threads' 'shared_data' location.
 */
template< bool DoScan , class FunctorType >
__device__
bool cuda_single_inter_block_reduce_scan( const FunctorType     & functor ,
                                          const Cuda::size_type   block_id ,
                                          const Cuda::size_type   block_count ,
                                          Cuda::size_type * const shared_data ,
                                          Cuda::size_type * const global_data ,
                                          Cuda::size_type * const global_flags )
{
  typedef Cuda::size_type                  size_type ;
  typedef FunctorValueTraits< FunctorType , void >  ValueTraits ;
  typedef FunctorValueJoin<   FunctorType , void >  ValueJoin ;
  typedef FunctorValueInit<   FunctorType , void >  ValueInit ;
  typedef FunctorValueOps<    FunctorType , void >  ValueOps ;

  typedef typename ValueTraits::pointer_type    pointer_type ;
  typedef typename ValueTraits::reference_type  reference_type ;

  const unsigned BlockSizeMask  = blockDim.y - 1 ;
  const unsigned BlockSizeShift = power_of_two_if_valid( blockDim.y );

  // Must have power of two thread count
  if ( BlockSizeMask & blockDim.y ) { cuda_abort("Cuda::cuda_single_inter_block_reduce_scan requires power-of-two blockDim"); }

  const integral_nonzero_constant< size_type , ValueTraits::StaticValueSize / sizeof(size_type) >
    word_count( ValueTraits::value_size( functor ) / sizeof(size_type) );

  // Reduce the accumulation for the entire block.
  cuda_intra_block_reduce_scan<false>( functor , pointer_type(shared_data) );

  {
    // Write accumulation total to global scratch space.
    // Accumulation total is the last thread's data.
    size_type * const shared = shared_data + word_count.value * BlockSizeMask ;
    size_type * const global = global_data + word_count.value * block_id ;

    for ( size_type i = threadIdx.y ; i < word_count.value ; i += blockDim.y ) { global[i] = shared[i] ; }
  }

  // Contributing blocks note that their contribution has been completed via an atomic-increment flag
  // If this block is not the last block to contribute to this group then the block is done.
  const bool is_last_block =
    ! __syncthreads_or( threadIdx.y ? 0 : ( 1 + atomicInc( global_flags , block_count - 1 ) < block_count ) );

  if ( is_last_block ) {

    const size_type b = ( long(block_count) * long(threadIdx.y) ) >> BlockSizeShift ;
    const size_type e = ( long(block_count) * long( threadIdx.y + 1 ) ) >> BlockSizeShift ;

    {
      void * const shared_ptr = shared_data + word_count.value * threadIdx.y ;
      reference_type shared_value = ValueInit::init( functor , shared_ptr );

      for ( size_type i = b ; i < e ; ++i ) {
        ValueJoin::join( functor , shared_ptr , global_data + word_count.value * i );
      }
    }

    cuda_intra_block_reduce_scan<DoScan>( functor , pointer_type(shared_data) );

    if ( DoScan ) {

      size_type * const shared_value = shared_data + word_count.value * ( threadIdx.y ? threadIdx.y - 1 : blockDim.y );

      if ( ! threadIdx.y ) { ValueInit::init( functor , shared_value ); }

      // Join previous inclusive scan value to each member
      for ( size_type i = b ; i < e ; ++i ) {
        size_type * const global_value = global_data + word_count.value * i ;
        ValueJoin::join( functor , shared_value , global_value );
        ValueOps ::copy( functor , global_value , shared_value );
      }
    }
  }

  return is_last_block ;
}

// Size in bytes required for inter block reduce or scan
template< bool DoScan , class FunctorType , class ArgTag >
inline
unsigned cuda_single_inter_block_reduce_scan_shmem( const FunctorType & functor , const unsigned BlockSize )
{
  return ( BlockSize + 2 ) * Impl::FunctorValueTraits< FunctorType , ArgTag >::value_size( functor );
}

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #if defined( __CUDACC__ ) */
#endif /* KOKKOS_CUDA_REDUCESCAN_HPP */

