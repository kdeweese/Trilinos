#ifndef PANZER_EVALUATOR_GATHER_SOLUTION_EPETRA_HPP
#define PANZER_EVALUATOR_GATHER_SOLUTION_EPETRA_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_CloneableEvaluator.hpp"

class Epetra_Vector;
class Epetra_CrsMatrix;

namespace Thyra {
template <typename ScalarT> class MultiVectorBase;
template <typename ScalarT> class LinearOpBase;
}

namespace panzer {

template <typename LocalOrdinalT,typename GlobalOrdinalT>
class UniqueGlobalIndexer; //forward declaration

/** \brief Gathers solution values from the Newton solution vector into 
    the nodal fields of the field manager

    Currently makes an assumption that the stride is constant for dofs
    and that the nmber of dofs is equal to the size of the solution
    names vector.
*/
template<typename EvalT, typename Traits> class GatherSolution_Epetra;

// **************************************************************
// **************************************************************
// * Specializations
// **************************************************************
// **************************************************************


// **************************************************************
// Residual 
// **************************************************************
template<typename Traits>
class GatherSolution_Epetra<panzer::Traits::Residual,Traits>
  : public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<panzer::Traits::Residual, Traits>,
    public panzer::CloneableEvaluator  {
   
  
public:
  
  GatherSolution_Epetra() {}

  GatherSolution_Epetra(const Teuchos::ParameterList& p);
  
  void postRegistrationSetup(typename Traits::SetupData d,
			     PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData d);

  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new GatherSolution_Epetra<panzer::Traits::Residual,Traits>(pl)); }
  
private:

  typedef typename panzer::Traits::Residual::ScalarT ScalarT;

  // maps the local (field,element,basis) triplet to a global ID
  // for scattering
  Teuchos::RCP<panzer::UniqueGlobalIndexer<panzer::Traits::LocalOrdinal,panzer::Traits::GlobalOrdinal> > globalIndexer_;
  std::vector<int> fieldIds_; // field IDs needing mapping

  std::vector< PHX::MDField<ScalarT,Cell,NODE> > gatherFields_;

  Epetra_Vector & getEpetraVector(const Teuchos::RCP<Thyra::MultiVectorBase<double> > & in_v) const;
};

// **************************************************************
// Jacobian
// **************************************************************
template<typename Traits>
class GatherSolution_Epetra<panzer::Traits::Jacobian,Traits>
  : public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<panzer::Traits::Jacobian, Traits>,
    public panzer::CloneableEvaluator  {
  
public:
  GatherSolution_Epetra() {}
  
  GatherSolution_Epetra(const Teuchos::ParameterList& p);
  
  void postRegistrationSetup(typename Traits::SetupData d,
			     PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData d);

  virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const
  { return Teuchos::rcp(new GatherSolution_Epetra<panzer::Traits::Jacobian,Traits>(pl)); }
  
private:

  typedef typename panzer::Traits::Jacobian::ScalarT ScalarT;

  // maps the local (field,element,basis) triplet to a global ID
  // for scattering
  Teuchos::RCP<panzer::UniqueGlobalIndexer<panzer::Traits::LocalOrdinal,panzer::Traits::GlobalOrdinal> > globalIndexer_;
  std::vector<int> fieldIds_; // field IDs needing mapping

  std::vector< PHX::MDField<ScalarT,Cell,NODE> > gatherFields_;

  Epetra_Vector & getEpetraVector(const Teuchos::RCP<Thyra::MultiVectorBase<double> > & in_v) const;
};

}

// **************************************************************

#include "Panzer_GatherSolution_EpetraT.hpp"

// **************************************************************
#endif
