// Copyright (C) 2006, 2007 Damien Hocking, KBC Advanced Technologies
// All Rights Reserved.
// This code is published under the Eclipse Public License.


#ifndef __IPEIGENSIMPLICIALSOLVERINTERFACE_HPP__ 
#define __IPEIGENSIMPLICIALSOLVERINTERFACE_HPP__

#include "IpSparseSymLinearSolverInterface.hpp"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
//#include <Eigen/Eigenvalues>


namespace Ipopt {


  /** Interface to the simplicial linear solver of Eigen, derived from
  *  SparseSymLinearSolverInterface.  For details, see description of
  *  SparseSymLinearSolverInterface base class.
  */
  class EigenSimplicialSolverInterface : public SparseSymLinearSolverInterface
  {
  public:
    /** @name Constructor/Destructor */
    //@{
    /** Constructor */
    EigenSimplicialSolverInterface();

    /** Destructor */
    virtual ~EigenSimplicialSolverInterface();
    //@}

    /** overloaded from AlgorithmStrategyObject */
    bool InitializeImpl(const OptionsList& options,
      const std::string& prefix) {
      return true;
    };

    /** @name Methods for requesting solution of the linear system. */
    //@{
    /** Method for initializing internal stuctures.  Here, ndim gives
    *  the number of rows and columns of the matrix, nonzeros give
    *  the number of nonzero elements, and airn and acjn give the
    *  positions of the nonzero elements.
    */
    virtual ESymSolverStatus InitializeStructure(Index dim, Index nonzeros,
      const Index *airn,
      const Index *ajcn);

    /** Method returing an internal array into which the nonzero
    *  elements (in the same order as airn and ajcn) are to be stored
    *  by the calling routine before a call to MultiSolve with a
    *  new_matrix=true.  The returned array must have space for at least
    *  nonzero elements. */
    virtual double* GetValuesArrayPtr();

    /** Solve operation for multiple right hand sides.  Overloaded
    *  from SparseSymLinearSolverInterface.
    */
    virtual ESymSolverStatus MultiSolve(bool new_matrix,
      const Index* airn,
      const Index* ajcn,
      Index nrhs,
      double* rhs_vals,
      bool check_NegEVals,
      Index numberOfNegEVals);

    /** Number of negative eigenvalues detected during last
    *  factorization.  Returns the number of negative eigenvalues of
    *  the most recent factorized matrix.  This must not be called if
    *  the linear solver does not compute this quantities (see
    *  ProvidesInertia).
    */
    virtual Index NumberOfNegEVals() const
    {
      return 0;
    }
    //@}

    //* @name Options of Linear solver */
    //@{
    /** Request to increase quality of solution for next solve.
    * Ask linear solver to increase quality of solution for the next
    * solve (e.g. increase pivot tolerance).  Returns false, if this
    * is not possible (e.g. maximal pivot tolerance already used.)
    */
    virtual bool IncreaseQuality()
    {
      return false;
    }

    /** Query whether inertia is computed by linear solver.
    * Returns true, if linear solver provides inertia.
    */
    virtual bool ProvidesInertia() const
    {
      return false;
    }
    /** Query of requested matrix type that the linear solver
    *  understands.
    */
    EMatrixFormat MatrixFormat() const
    {
      return CSR_Full_Format_0_Offset;
    }
    //@}

    /** Methods for IpoptType */
    //@{
    static void RegisterOptions(SmartPtr<RegisteredOptions> roptions) {};
    //@}

    /** Query whether the indices of linearly dependent rows/columns
    *  can be determined by this linear solver. */
    virtual bool ProvidesDegeneracyDetection() const
    {
      return false;
    }

    /** This method determines the list of row indices of the linearly
    *  dependent rows. */
    virtual ESymSolverStatus DetermineDependentRows(const Index* ia,
      const Index* ja,
      std::list<Index>& c_deps)
    {
      return SYMSOLVER_FATAL_ERROR;
    }

  private:
    /**@name Default Compiler Generated Methods
    * (Hidden to avoid implicit creation/calling).
    * These methods are not implemented and
    * we do not want the compiler to implement
    * them for us, so we declare them private
    * and do not define them. This ensures that
    * they will not be implicitly created/called. */
    //@{
    /** Copy Constructor */
    //EigenSimplicialSolverInterface(const EigenSimplicialSolverInterface&);

    /** Overloaded Equals Operator */
    //void operator=(const EigenSimplicialSolverInterface&);
    //@}

    typedef Eigen::VectorXd VectorType;
    typedef Eigen::SparseMatrix<double,Eigen::ColMajor,Index> SparseMatrixType; 

    //typedef Eigen::SimplicialLDLT<SparseMatrixType,Eigen::Upper,Eigen::NaturalOrdering<SparseMatrixType::StorageIndex>> SparseSolverType;
    typedef Eigen::SimplicialLDLT<SparseMatrixType> SparseSolverType;

    //Eigen::SparseMatrix
    Index mDim;
    
    SparseSolverType* mSolver;
    SparseMatrixType* mMatrix;

  };

}

#endif
