// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_IMAGEBASEDPFACTORY_DEF_HPP
#define MUELU_IMAGEBASEDPFACTORY_DEF_HPP

#include "Xpetra_CrsGraph.hpp"
#include "Xpetra_CrsMatrixUtils.hpp"

#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Aggregates.hpp"

// Including this one last ensure that the short names of the above headers are defined properly
#include "MueLu_ImageBasedPFactory_decl.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> ImageBasedPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
    SET_VALID_ENTRY("interp: build coarse coordinates");
#undef  SET_VALID_ENTRY

    // general variables needed in ImageBasedPFactory
    validParamList->set<RCP<const FactoryBase> >("A",                       Teuchos::null,
                                                 "Generating factory of the matrix A");
    validParamList->set<RCP<const FactoryBase> >("Aggregates",                   Teuchos::null,
                                                 "Aggregates generated by StructuredAggregationFactory used to construct a piece-constant prolongator.");
    validParamList->set<RCP<const FactoryBase> >("prolongatorGraph",             Teuchos::null,
                                                 "Graph generated by StructuredAggregationFactory used to construct a piece-linear prolongator.");
    validParamList->set<RCP<const FactoryBase> >("Coordinates",                  Teuchos::null,
                                                 "Fine level coordinates used to construct piece-wise linear prolongator and coarse level coordinates.");
    validParamList->set<RCP<const FactoryBase> >("Greyscale",
						 Teuchos::null, "Greyscale pixel information.");
    validParamList->set<RCP<const FactoryBase> >("coarseCoordinatesFineMap",     Teuchos::null,
                                                 "map of the coarse coordinates' GIDs as indexed on the fine mesh.");
    validParamList->set<RCP<const FactoryBase> >("coarseCoordinatesMap",         Teuchos::null,
                                                 "map of the coarse coordinates' GIDs as indexed on the coarse mesh.");
    validParamList->set<RCP<const FactoryBase> >("Nullspace",                    Teuchos::null,
                                                 "Fine level nullspace used to construct the coarse level nullspace.");
    validParamList->set<RCP<const FactoryBase> >("numDimensions",                Teuchos::null,
                                                 "Number of spacial dimensions in the problem.");
    validParamList->set<RCP<const FactoryBase> >("lCoarseNodesPerDim",           Teuchos::null,
                                                 "Number of nodes per spatial dimension on the coarse grid.");
    validParamList->set<RCP<const FactoryBase> >("lNodesPerDim",           Teuchos::null,
                                                 "Number of nodes per spatial dimension on the fine grid.");
    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ImageBasedPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  DeclareInput(Level& fineLevel, Level& coarseLevel) const {
    const ParameterList& pL = GetParameterList();

    Input(fineLevel, "A");
    Input(fineLevel, "Nullspace");
    Input(fineLevel, "numDimensions");
    Input(fineLevel, "prolongatorGraph");
    Input(fineLevel, "lCoarseNodesPerDim");
    Input(fineLevel, "Greyscale");
    
    if( pL.get<bool>("interp: build coarse coordinates") ) {
      Input(fineLevel, "Coordinates");
      Input(fineLevel, "coarseCoordinatesFineMap");
      Input(fineLevel, "coarseCoordinatesMap");
    }

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ImageBasedPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  Build(Level& fineLevel, Level &coarseLevel) const {
    return BuildP(fineLevel, coarseLevel);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ImageBasedPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  BuildP(Level &fineLevel, Level &coarseLevel) const {
    FactoryMonitor m(*this, "BuildP", coarseLevel);

    // Set debug outputs based on environment variable
    RCP<Teuchos::FancyOStream> out;
    //if(const char* dbg = std::getenv("MUELU_IMAGEBASEDPFACTORY_DEBUG")) {
      out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      out->setShowAllFrontMatter(false).setShowProcRank(true);
      //} else {
      //out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
      //}

    *out << "Starting ImageBasedPFactory::BuildP." << std::endl;

    // Get inputs from the parameter list
    const ParameterList& pL = GetParameterList();
    const bool buildCoarseCoordinates = pL.get<bool>("interp: build coarse coordinates");
    const int numDimensions           = Get<int>(fineLevel, "numDimensions");

    // Declared main input/outputs to be retrieved and placed on the fine resp. coarse level
    RCP<Matrix> A = Get<RCP<Matrix> >(fineLevel, "A");
    RCP<const CrsGraph> prolongatorGraph = Get<RCP<CrsGraph> >(fineLevel, "prolongatorGraph");
    RCP<realvaluedmultivector_type> fineCoordinates, coarseCoordinates;
    RCP<Matrix> P;
    Array<LO> lFineNodesPerDir = Get<Array<LO> >(fineLevel, "lNodesPerDim");
    
    RCP<realvaluedmultivector_type> greyscale;
    greyscale = Get< RCP<realvaluedmultivector_type> >(fineLevel, "Greyscale");
    //greyscale->describe(*out,Teuchos::VERB_HIGH);
    // Check if we need to build coarse coordinates as they are used if we construct
    // a linear interpolation prolongator
    if(buildCoarseCoordinates) {
      SubFactoryMonitor sfm(*this, "BuildCoordinates", coarseLevel);
      RCP<const Map> coarseCoordsFineMap = Get< RCP<const Map> >(fineLevel, "coarseCoordinatesFineMap");
      RCP<const Map> coarseCoordsMap = Get< RCP<const Map> >(fineLevel, "coarseCoordinatesMap");
      fineCoordinates   = Get< RCP<realvaluedmultivector_type> >(fineLevel, "Coordinates");
      coarseCoordinates = Xpetra::MultiVectorFactory<real_type,LO,GO,Node>::Build(coarseCoordsFineMap,
                                                                                  fineCoordinates->getNumVectors());
      RCP<const Import> coordsImporter = ImportFactory::Build(fineCoordinates->getMap(),
                                                              coarseCoordsFineMap);
      coarseCoordinates->doImport(*fineCoordinates, *coordsImporter, Xpetra::INSERT);
      coarseCoordinates->replaceMap(coarseCoordsMap);

      Set(coarseLevel, "Coordinates", coarseCoordinates);
    }

    *out << "Fine and coarse coordinates have been loaded from the fine level and set on the coarse level." << std::endl;

     
    // Compute the prolongator using piece-wise linear interpolation
    // First get all the required coordinates to compute the local part of P
    RCP<realvaluedmultivector_type> ghostCoordinates
      = Xpetra::MultiVectorFactory<real_type,LO,GO,NO>::Build(prolongatorGraph->getColMap(),
							      fineCoordinates->getNumVectors());
    RCP<const Import> ghostImporter = ImportFactory::Build(coarseCoordinates->getMap(),
							   prolongatorGraph->getColMap());
    ghostCoordinates->doImport(*coarseCoordinates, *ghostImporter, Xpetra::INSERT);
    
    SubFactoryMonitor sfm(*this, "BuildLinearP", coarseLevel);
    BuildImageP(A, prolongatorGraph, fineCoordinates, ghostCoordinates, greyscale, lFineNodesPerDir, numDimensions, P);
    
  
  *out << "The prolongator matrix has been built." << std::endl;
  
  {
    SubFactoryMonitor sfm(*this, "BuildNullspace", coarseLevel);
    // Build the coarse nullspace
    RCP<MultiVector> fineNullspace   = Get< RCP<MultiVector> > (fineLevel, "Nullspace");
    RCP<MultiVector> coarseNullspace = MultiVectorFactory::Build(P->getDomainMap(),
								 fineNullspace->getNumVectors());
    P->apply(*fineNullspace, *coarseNullspace, Teuchos::TRANS, Teuchos::ScalarTraits<SC>::one(),
	     Teuchos::ScalarTraits<SC>::zero());
    Set(coarseLevel, "Nullspace", coarseNullspace);
  }
  
  *out << "The coarse nullspace is constructed and set on the coarse level." << std::endl;
  
  
  
  Array<LO> lNodesPerDir = Get<Array<LO> >(fineLevel, "lCoarseNodesPerDim");
  
  Set(coarseLevel, "numDimensions", numDimensions);
  Set(coarseLevel, "lNodesPerDim", lNodesPerDir);
  Set(coarseLevel, "P", P);
  
  *out << "ImageBasedPFactory::BuildP has completed." << std::endl;
  
} // BuildP

  
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ImageBasedPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  BuildImageP(RCP<Matrix>& A, RCP<const CrsGraph>& prolongatorGraph,
               RCP<realvaluedmultivector_type>& fineCoordinates,
               RCP<realvaluedmultivector_type>& ghostCoordinates,
	       RCP<realvaluedmultivector_type>& greyscale,
	       Array<LO> &lFineNodesPerDim,
               const int numDimensions, RCP<Matrix>& P) const {

    // Set debug outputs based on environment variable
    RCP<Teuchos::FancyOStream> out;
    //if(const char* dbg = std::getenv("MUELU_IMAGEBASEDPFACTORY_DEBUG")) {
      out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      out->setShowAllFrontMatter(false).setShowProcRank(true);
      //} else {
      //out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
      //}
      //A->describe(*out, Teuchos::VERB_EXTREME);
    *out << "Entering BuildLinearP" << std::endl;
    //prolongatorGraph->describe(*out, Teuchos::VERB_EXTREME);
    ArrayRCP<const real_type> grey(greyscale->getData(0));
    Array<real_type> pixels(grey.size());
    double alpha;
    double K1 = 1.0;
    double K2 = 100.;
    for(size_t i=0; i < grey.size(); ++i) {
      alpha=GetBlendingParameter(grey[i], 1);
      pixels[i]=(1.0-alpha)*K1+(alpha*K2);
    }
    // Extract coordinates for interpolation stencil calculations
    const LO numFineNodes  = fineCoordinates->getLocalLength();
    const LO numGhostNodes = ghostCoordinates->getLocalLength();
    Array<ArrayRCP<const real_type> > fineCoords(3);
    Array<ArrayRCP<const real_type> > ghostCoords(3);
    const real_type realZero = Teuchos::as<real_type>(0.0);
    ArrayRCP<real_type> fineZero(numFineNodes, realZero);
    ArrayRCP<real_type> ghostZero(numGhostNodes, realZero);
    for(int dim = 0; dim < 3; ++dim) {
      if(dim < numDimensions) {
        fineCoords[dim]  = fineCoordinates->getData(dim);
        ghostCoords[dim] = ghostCoordinates->getData(dim);
      } else {
        fineCoords[dim]  = fineZero;
        ghostCoords[dim] = ghostZero;
      }
    }

    *out << "Coordinates extracted from the multivectors!" << std::endl;

    // Compute 2^numDimensions using bit logic to avoid round-off errors
    const int numInterpolationPoints = 1 << numDimensions;
    const int dofsPerNode = A->GetFixedBlockSize();

    std::vector<size_t> strideInfo(1);
    strideInfo[0] = dofsPerNode;
    RCP<const StridedMap> stridedDomainMap =
      StridedMapFactory::Build(prolongatorGraph->getDomainMap(), strideInfo);

    *out << "The maps of P have been computed" << std::endl;

    RCP<ParameterList> dummyList = rcp(new ParameterList());
    P = rcp(new CrsMatrixWrap(prolongatorGraph, dummyList));
    RCP<CrsMatrix> PCrs = rcp_dynamic_cast<CrsMatrixWrap>(P)->getCrsMatrix();
    PCrs->resumeFill(); // The Epetra matrix is considered filled at this point.
    std::cout << "numInterpolationPoints " << numInterpolationPoints << std::endl;
    LO interpolationNodeIdx = 0, rowIdx = 0;
    ArrayView<const LO> colIndices;
    Array<SC> values;
    Array<Array<real_type> > coords(numInterpolationPoints + 1);
    Array<real_type> stencil(numInterpolationPoints);
    Array<int> match(3);
    for(LO nodeIdx = 0; nodeIdx < numFineNodes; ++nodeIdx) {
      if(PCrs->getNumEntriesInLocalRow(nodeIdx*dofsPerNode) == 1) {
        values.resize(1);
        values[0] = 1.0;
        for(LO dof = 0; dof < dofsPerNode; ++dof) {
          rowIdx = nodeIdx*dofsPerNode + dof;
          prolongatorGraph->getLocalRowView(rowIdx, colIndices);
          PCrs->replaceLocalValues(rowIdx, colIndices, values());
        }
      } else {
        // Extract the coordinates associated with the current node
        // and the neighboring coarse nodes
        coords[0].resize(3);
        for(int dim = 0; dim < 3; ++dim) {
          coords[0][dim] = fineCoords[dim][nodeIdx];
	}
	match[0]=0;
	match[1]=0;
	match[2]=0;
        prolongatorGraph->getLocalRowView(nodeIdx*dofsPerNode, colIndices);
        for(int interpolationIdx=0; interpolationIdx < numInterpolationPoints; ++interpolationIdx) {
          coords[interpolationIdx + 1].resize(3);
          interpolationNodeIdx = colIndices[interpolationIdx] / dofsPerNode;
          for(int dim = 0; dim < 3; ++dim) {
            coords[interpolationIdx + 1][dim] = ghostCoords[dim][interpolationNodeIdx];
	    if(coords[interpolationIdx+1][dim] == coords[0][dim]) {
	      match[dim]++;
	    }
          }
        }
	values.resize(numInterpolationPoints);
	real_type kappa_nw, kappa_ne, kappa_sw, kappa_se, s, w_N, w_S, w_E, w_W,
	  w_NW, w_SE, w_SW, w_NE;
	kappa_nw = nodeIdx % (lFineNodesPerDim)[0] == 0 ? 0 : pixels[nodeIdx-1];
	kappa_ne = nodeIdx+1 % (lFineNodesPerDim)[0] == 0 ? 0 :  pixels[nodeIdx];
	kappa_sw = (nodeIdx % (lFineNodesPerDim)[0] == 0) || (nodeIdx < (lFineNodesPerDim)[0]) ? 0 : pixels[nodeIdx-(lFineNodesPerDim)[0]-1];
	kappa_se = nodeIdx < (lFineNodesPerDim)[0] ? 0 : pixels[nodeIdx-(lFineNodesPerDim)[0]];
	s = kappa_nw + kappa_ne + kappa_sw + kappa_se;
	if(match[0]==2) {
	  w_N = (kappa_nw + kappa_ne)/s;
	  w_S = (kappa_sw + kappa_se)/s;
	  if(nodeIdx-1 % lFineNodesPerDim[0]==0) {
	    values[0]=0;
	    values[1]=w_S;
	    values[2]=0;
	    values[3]=w_N;
	  }
	  else {
	    values[0]=w_S;
	    values[1]=0;
	    values[2]=w_N;
	    values[3]=0;
	  }
	}
	else if(match[1]==2) {
	  w_W = (kappa_sw + kappa_nw)/s;
	  w_E = (kappa_se + kappa_ne)/s;

	  if(numFineNodes - nodeIdx <= lFineNodesPerDim[0]) {
	    values[0]=0;
	    values[1]=0;
	    values[2]=w_W;
	    values[3]=w_E;
	  }
	  else {
	    values[0]=w_W;
	    values[1]=w_E;
	    values[2]=0;
	    values[3]=0;
	  }
	}
	else if(match[0]==0 && match[1]==0) {
	  real_type w_N_l, w_N_r, w_S_l, w_S_r, w_W_u, w_W_d, w_E_u, w_E_d;
	  w_N = (kappa_nw + kappa_ne)/s;
	  w_S = (kappa_sw + kappa_se)/s;
	  w_W = (kappa_sw + kappa_nw)/s;
	  w_E = (kappa_se + kappa_ne)/s;

	  //s_i,j+1
	  kappa_nw = pixels[nodeIdx-1+lFineNodesPerDim[0]];
	  kappa_ne = pixels[nodeIdx+lFineNodesPerDim[0]];
	  kappa_sw = nodeIdx % (lFineNodesPerDim)[0] == 0 ? 0 : pixels[nodeIdx-1];
	  kappa_se = pixels[nodeIdx];
	  s = kappa_nw + kappa_ne + kappa_sw + kappa_se;
	  w_W_u = (kappa_sw + kappa_nw)/s;
	  w_E_u = (kappa_se + kappa_ne)/s;
	  
	  //s_i-1,j
	  kappa_nw = nodeIdx-1 % lFineNodesPerDim[0] == 0 ? 0 : pixels[nodeIdx-2];
	  kappa_ne = pixels[nodeIdx-1];
	  kappa_sw = nodeIdx-1 & lFineNodesPerDim[0] == 0 ? 0 : pixels[nodeIdx-lFineNodesPerDim[0]-2];
	  kappa_se = pixels[nodeIdx-lFineNodesPerDim[0]-1];
	  s = kappa_nw + kappa_ne + kappa_sw + kappa_se;
	  w_N_l = (kappa_nw + kappa_ne)/s;
	  w_S_l = (kappa_sw + kappa_se)/s;
	  
	  //s_i+1,j
	  kappa_nw = pixels[nodeIdx];
	  kappa_ne = nodeIdx-1 % lFineNodesPerDim[0] == 0 ? 0 : pixels[nodeIdx+1];
	  kappa_sw = pixels[nodeIdx-lFineNodesPerDim[0]];
	  kappa_se = nodeIdx-1 % lFineNodesPerDim[0] == 0 ? 0 : pixels[nodeIdx-lFineNodesPerDim[0]+1];
	  s = kappa_nw + kappa_ne + kappa_sw + kappa_se;
	  w_N_r = (kappa_nw + kappa_ne)/s;
	  w_S_r = (kappa_sw + kappa_se)/s;
	  
	  //s_i,j-1
	  kappa_nw = pixels[nodeIdx-1];
	  kappa_ne = pixels[nodeIdx];
	  kappa_sw = pixels[nodeIdx-lFineNodesPerDim[0]-1];
	  kappa_se = pixels[nodeIdx-lFineNodesPerDim[0]];
	  s = kappa_nw + kappa_ne + kappa_sw + kappa_se;
	  w_W_d = (kappa_sw + kappa_nw)/s;
	  w_E_d = (kappa_se + kappa_ne)/s;
	  

	  w_NW = .5*(w_N * w_W_u + w_W * w_N_l);
	  w_NE = .5*(w_E * w_N_r + w_N * w_E_u);
	  w_SW = .5*(w_W * w_S_l + w_S * w_W_d);
	  w_SE = .5*(w_S * w_E_d + w_E * w_S_r);
	  
	  values[0]=w_SW;
	  values[1]=w_SE;
	  values[2]=w_NW;
	  values[3]=w_NE;
	}
	else {
	  std::cout << "warning" << std::endl;
	}
        //ComputeLinearInterpolationStencil(numDimensions, numInterpolationPoints, coords, stencil);
	/*for(int i=0; i < numInterpolationPoints; ++i) {
	  std::cout << "stencil: " << i << " " << stencil[i] << std::endl;
	}
        values.resize(numInterpolationPoints);
        for(LO valueIdx = 0; valueIdx < numInterpolationPoints; ++valueIdx) {
          values[valueIdx] = Teuchos::as<SC>(stencil[valueIdx]);
	  }*/

        // Set values in all the rows corresponding to nodeIdx
        for(LO dof = 0; dof < dofsPerNode; ++dof) {
          rowIdx = nodeIdx*dofsPerNode + dof;
          prolongatorGraph->getLocalRowView(rowIdx, colIndices);
          PCrs->replaceLocalValues(rowIdx, colIndices, values());
        }
      }
    }

    *out << "The calculation of the interpolation stencils has completed." << std::endl;

    PCrs->fillComplete();
    std::cout << "PCrs" << std::endl;
    PCrs->describe(*out,Teuchos::VERB_EXTREME);
    //Teuchos::RCP<const CrsMatrix> temp(PCrs);
    Xpetra::IO<real_type,LO,GO,Node>::Write("test.mtx",*P);
    *out << "All values in P have been set and expertStaticFillComplete has been performed." << std::endl;

    // set StridingInformation of P
    if (A->IsView("stridedMaps") == true) {
      P->CreateView("stridedMaps", A->getRowMap("stridedMaps"), stridedDomainMap);
    } else {
      P->CreateView("stridedMaps", P->getRangeMap(), stridedDomainMap);
    }

  } // BuildLinearP


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Scalar ImageBasedPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetBlendingParameter(const real_type pixelval, int method) const {
    double thresh=77;
    double alpha;
    //Linear
    if(method==0) {
      alpha=pixelval/255.;
    }
    //Step
    if(method==1) {
      if(pixelval <= thresh)
	alpha=0.0;
      else
	alpha=1.0;
    }
    //Piecewise
    if(method==2) {
      if(pixelval <= thresh)
	alpha=pixelval/142.;
      else
	alpha=(pixelval+113.0)/368.0;
    }
    //Sigmoid
    if(method==3) {
      alpha=1.0/(1.0+exp(-1.0*(pixelval-thresh)));
    }
  }
  
} // namespace MueLu

#endif // MUELU_IMAGEBASEDPFACTORY_DEF_HPP
