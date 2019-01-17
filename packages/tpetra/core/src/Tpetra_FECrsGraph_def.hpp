// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef TPETRA_FECRSGRAPH_DEF_HPP
#define TPETRA_FECRSGRAPH_DEF_HPP

#include <type_traits>
#include "Tpetra_CrsGraph.hpp"

//#define USE_UNALIASED_MEMORY



namespace Tpetra {

template<class LocalOrdinal, class GlobalOrdinal, class Node>
FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
FECrsGraph(const Teuchos::RCP<const map_type> & ownedRowMap,
           const Teuchos::RCP<const map_type> & ownedPlusSharedRowMap, 
           const size_t maxNumEntriesPerRow,
           const Teuchos::RCP<const import_type> & ownedPlusSharedToOwnedimporter,
           const Teuchos::RCP<const map_type> & domainMap,
           const Teuchos::RCP<const map_type> & rangeMap,
           const Teuchos::RCP<Teuchos::ParameterList>& params): 
  crs_graph_type(ownedPlusSharedRowMap, maxNumEntriesPerRow, StaticProfile, params),
  importer_(ownedPlusSharedToOwnedimporter),
  domainMap_(domainMap.is_null() ? ownedRowMap : domainMap),
  rangeMap_(rangeMap.is_null() ? ownedRowMap : rangeMap)
{  
  setup(ownedRowMap,ownedPlusSharedRowMap,params,maxNumEntriesPerRow);
}


template<class LocalOrdinal, class GlobalOrdinal, class Node>
FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
FECrsGraph (const Teuchos::RCP<const map_type> & ownedRowMap,
            const Teuchos::RCP<const map_type> & ownedPlusSharedRowMap, 
            const Kokkos::DualView<const size_t*, execution_space>& numEntPerRow,
            const Teuchos::RCP<const import_type> & ownedPlusSharedToOwnedimporter,
            const Teuchos::RCP<const map_type> & domainMap,
            const Teuchos::RCP<const map_type> & rangeMap,
            const Teuchos::RCP<Teuchos::ParameterList>& params):
  crs_graph_type( ownedPlusSharedRowMap, numEntPerRow, StaticProfile, params),
  importer_(ownedPlusSharedToOwnedimporter),
  domainMap_(domainMap.is_null() ? ownedRowMap : domainMap),
  rangeMap_(rangeMap.is_null() ? ownedRowMap : rangeMap)

{  
  // Only pass in numEntries for "owned rows"
  size_t numOwnedRows = ownedRowMap->getNodeNumElements();
  auto sv = Kokkos::subview(numEntPerRow,Kokkos::pair<size_t,size_t>(0,numOwnedRows));
  setup(ownedRowMap,ownedPlusSharedRowMap,params,sv);
}


template<class LocalOrdinal, class GlobalOrdinal, class Node>
template <class NumEntries_t>
void FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::setup(const Teuchos::RCP<const map_type>  & ownedRowMap, const Teuchos::RCP<const map_type> & ownedPlusSharedRowMap,const Teuchos::RCP<Teuchos::ParameterList>& params, NumEntries_t &ne) {
 const char tfecfFuncName[] = "FECrsGraph::setup(): ";

 TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(ownedRowMap.is_null (), std::runtime_error, "ownedRowMap is null.");
 TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(ownedPlusSharedRowMap.is_null (), std::runtime_error, "ownedPlusSharedRowMap is null.");

 activeCrsGraph_     = Teuchos::rcp(new FEWhichActive(FE_ACTIVE_OWNED_PLUS_SHARED));

 // NOTE: We're forcing the CrsGraph to be in global index mode 
 this->allocateIndices(GlobalIndices);
 
 // Use a very strong map equivalence check
 bool maps_are_the_same = ownedRowMap->isSameAs(*ownedPlusSharedRowMap);
 if(!maps_are_the_same) {
   // Make an importer if we need to, check map compatability if we don't
   if(importer_.is_null()) {
       importer_ = Teuchos::rcp(new import_type(ownedRowMap,ownedPlusSharedRowMap));
   } 
   else {
     TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(!ownedRowMap->isSameAs(*importer_->getSourceMap()), std::runtime_error, "ownedRowMap does not match importer source map.");
     TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(!ownedPlusSharedRowMap->isSameAs(*importer_->getTargetMap()), std::runtime_error, "ownedPlusSharedRowMap does not match importer target map.");
   }

   // Make sure the ownedPlusSharedRowMap has at least as many entries at the ownedRowMap (due to our superset requriement)
   TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( importer_->getNumSameIDs() != importer_->getSourceMap()->getNodeNumElements(),
                                          std::runtime_error,"ownedRowMap contains entries which are not in the ownedPlusSharedRowMap.");   
 
   TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( ownedRowMap->getNodeNumElements() > ownedPlusSharedRowMap->getNodeNumElements(),
                                          std::runtime_error,"ownedRowMap more entries than the ownedPlusSharedRowMap.");   

   // Build the inactive graph
#ifdef USE_UNALIASED_MEMORY
   inactiveCrsGraph_ = Teuchos::rcp(new crs_graph_type(ownedRowMap,ne,StaticProfile,params));
   inactiveCrsGraph_->allocateIndices(GlobalIndices);

#else
   // For FECrsGraph, we do all the aliasing AFTER import.  All we need here is a constructor
   inactiveCrsGraph_ = Teuchos::rcp(new crs_graph_type(ownedRowMap,ne,StaticProfile,params));
#ifdef OLD_STUFF
   inactiveCrsGraph_->allocateIndices(GlobalIndices); // FIXME: We don't actually need to allocate these arrays
#endif
#endif
 }

}



template<class LocalOrdinal, class GlobalOrdinal, class Node>
FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>&
FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
operator=(const FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>& rhs)
{
  return *this;
}

template<class LocalOrdinal, class GlobalOrdinal, class Node>
void FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::doOwnedPlusSharedToOwned(const CombineMode CM) {
  if(!inactiveCrsGraph_.is_null() && *activeCrsGraph_ == FE_ACTIVE_OWNED_PLUS_SHARED) {
#ifdef USE_UNALIASED_MEMORY
    inactiveCrsGraph_->doExport(*this,*importer_,CM);
    
    // fillComplete the owned graph
    inactiveCrsGraph_->fillComplete(domainMap_,rangeMap_);

    // fillComplete the owned+shared graph in a way that generates the owned+shared grep w/o an importer or exporter
    crs_graph_type::fillComplete(inactiveCrsGraph_->getColMap(),inactiveCrsGraph_->getRowMap());
#else

#if 0
    int rank = this->getRowMap()->getComm()->getRank();
#endif
    Teuchos::RCP<const map_type> ownedRowMap = inactiveCrsGraph_->getRowMap();

    // Do a self-export in "restricted mode"
    this->doExport(*this,*importer_,CM,true);

    // Under the "if you own an element, you own a node" assumption, we can start by making a columnmap for ownedPlusShared
    // Make the column for the owned guy
    Teuchos::Array<int> remotePIDs (0);
    this->makeColMap(remotePIDs);

    // Now run CrsGraph's fillComplete to get the final importer
    crs_graph_type::fillComplete(domainMap_,this->getRowMap());

    // Time to build an owned localGraph via subviews
    local_graph_type ownedPlusSharedGraph = this->getLocalGraph();
    size_t numOwnedRows = ownedRowMap->getNodeNumElements();
    // FIXME: This uses UVM.
    size_t numOwnedNonZeros = ownedPlusSharedGraph.row_map[numOwnedRows];

    // Build the inactive guy
    // FIME: Shouldn't have to do this in the constructor as well
    //    inactiveCrsGraph_ = Teuchos::rcp(new crs_graph_type(ownedGraph,ownedRowMap,this->getColmap(),domainMap_,rangeMap_));
    // NOTE: We can't use the local_graph_type constructor, because it does not allow us to provide an importer
    inactiveCrsGraph_->replaceColMap(this->getColMap());
    inactiveCrsGraph_->setAllIndices(Kokkos::subview(ownedPlusSharedGraph.row_map,Kokkos::pair<size_t,size_t>(0,numOwnedRows+1)),
                                     Kokkos::subview(ownedPlusSharedGraph.entries,Kokkos::pair<size_t,size_t>(0,numOwnedNonZeros)));
    // FIXME: We do NOT want the exporter from this.  Right now I assume range == row, but this needs to be relaxed
    inactiveCrsGraph_->expertStaticFillComplete(domainMap_,rangeMap_,this->getImporter(),Teuchos::null);

#warning "Tpetra::FECrsGraph's implementation is not complete with USE_UNALIASED_MEMORY unset"

#endif
  }
}//end doOverlapToLocal


template<class LocalOrdinal, class GlobalOrdinal, class Node>
void FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::doOwnedToOwnedPlusShared(const CombineMode CM) {
  // This should be a no-op for all of our purposes
}//end doLocalToOverlap

template<class LocalOrdinal, class GlobalOrdinal, class Node>
void FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::switchActiveCrsGraph() {
  if(*activeCrsGraph_ == FE_ACTIVE_OWNED_PLUS_SHARED)
    *activeCrsGraph_ = FE_ACTIVE_OWNED;
  else
    *activeCrsGraph_ = FE_ACTIVE_OWNED_PLUS_SHARED;

  if(inactiveCrsGraph_.is_null()) return;

  this->swap(*inactiveCrsGraph_);

}//end switchActiveCrsGraph



template<class LocalOrdinal, class GlobalOrdinal, class Node>
void FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::endFill() {
  const char tfecfFuncName[] = "FECrsGraph::endFill(): ";

  /* What has to go on here is complicated.  
     First off, if we don't really have two graphs (e.g. the rowMaps are the same, because we're in serial or 
     doing finite differences, things are easy --- just call fillComplete().

     If, we are in the parallel FE case, then:
     Precondition: FE_ACTIVE_OWNED_PLUS_SHARED mode

     Postconditions: 
     1) FE_ACTIVE_OWNED mode
     2) The OWNED graph has been fillCompleted with an Aztec-compatible column map
     3) rowptr & (local) colinds are aliased between the two graphs
     4) The OWNED_PLUS_SHARED graph has been fillCompleted with a column map whose first chunk
        is the column map for the OWNED graph.  
        If we assume that (a) if you own an element, you also own at least one of the connecte nodes and (b) elements are cliques, then
        the columnMap is the same for both graphs!!! Yay!!!

     5) The OWNED_PLUS_SHARED graph has neither an importer nor exporter.  Making these is expensive and we don't need them.       
   */
  // Precondition
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(*activeCrsGraph_ != FE_ACTIVE_OWNED_PLUS_SHARED,std::runtime_error, "must be in owned+shared mode.");

  if(inactiveCrsGraph_.is_null()) {
    // The easy case: One graph
    switchActiveCrsGraph();
   crs_graph_type::fillComplete(domainMap_,rangeMap_);
  }
  else {
    // The hard case: Two graphs   

    // fillComplete the owned+shared graph in a way that generates the owned+shared grep w/o an importer or exporter
    // FIXME: This makes an importer.  DO NOT WANT
    // Migrate data to the owned graph
    doOwnedPlusSharedToOwned(Tpetra::ADD);

    // Load up the owned graph
    switchActiveCrsGraph();

  }
}


template<class LocalOrdinal, class GlobalOrdinal, class Node>
void FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::beginFill() {
  const char tfecfFuncName[] = "FECrsGraph::beginFill(): ";
  
  // Unlike FECrsMatrix and FEMultiVector, we do not allow you to call beginFill() after calling endFill()
  // So we throw an exception if you're in owned mode
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(*activeCrsGraph_ == FE_ACTIVE_OWNED,std::runtime_error, "can only be called once.");

}


}  // end namespace Tpetra


//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//
#define TPETRA_FECRSGRAPH_GRAPH_INSTANT(LO,GO,NODE) \
  template class FECrsGraph<LO, GO, NODE>;



#endif // TPETRA_FECRSGRAPH_DEF
