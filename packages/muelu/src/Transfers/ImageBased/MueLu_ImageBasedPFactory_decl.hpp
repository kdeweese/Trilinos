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
#ifndef MUELU_IMAGEBASEDPFACTORY_DECL_HPP
#define MUELU_IMAGEBASEDPFACTORY_DECL_HPP

// Teuchos includes for dense linear algebra
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialDenseSolver.hpp>

#include "Xpetra_CrsGraph.hpp"

#include "MueLu_PFactory.hpp"
#include "MueLu_Level_fwd.hpp"
#include "MueLu_Aggregates_fwd.hpp"

namespace MueLu{

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class ImageBasedPFactory : public PFactory {
#undef MUELU_IMAGEBASEDPFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    // Declare useful types
    using real_type = typename Teuchos::ScalarTraits<SC>::coordinateType;
    using realvaluedmultivector_type = Xpetra::MultiVector<real_type,LO,GO,Node>;

    //! @name Constructors/Destructors.
    //@{

    //! Constructor
    ImageBasedPFactory() { }

    //! Destructor.
    virtual ~ImageBasedPFactory() { }
    //@}

    RCP<const ParameterList> GetValidParameterList() const;

    //! Input
    //@{

    void DeclareInput(Level& fineLevel, Level& coarseLevel) const;

    //@}

    //! @name Build methods.
    //@{

    void Build (Level& fineLevel, Level& coarseLevel) const;
    void BuildP(Level& fineLevel, Level& coarseLevel) const;

    //@}

  private:
    void BuildImageP(RCP<Matrix>& A, RCP<const CrsGraph>& prolongatorGraph,
                      RCP<realvaluedmultivector_type>& fineCoordinates,
                      RCP<realvaluedmultivector_type>& ghostCoordinates,
		      RCP<realvaluedmultivector_type>& greyscale,
		      Array<LO> &lFineNodesPerDir,
                      const int numDimensions, RCP<Matrix>& P) const;

    real_type GetBlendingParameter(const real_type pixelval, int method) const;
  }; // class ImageBasedPFactory

} // namespace MueLu

#define MUELU_IMAGEBASEDPFACTORY_SHORT
#endif // MUELU_IMAGEBASEDPFACTORY_DECL_HPP
