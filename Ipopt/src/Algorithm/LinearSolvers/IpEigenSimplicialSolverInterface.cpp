#include "IpEigenSimplicialSolverInterface.hpp"

namespace Ipopt
{

  EigenSimplicialSolverInterface::EigenSimplicialSolverInterface()
  {
    this->mDim = 0;
  }


  EigenSimplicialSolverInterface::~EigenSimplicialSolverInterface()
  {
    delete this->mMatrix;
    //delete this->mSolver;
  }

  ESymSolverStatus EigenSimplicialSolverInterface::InitializeStructure(Index dim, Index nonzeros,
    const Index *airn,
    const Index *ajcn)
  {
    this->mDim = dim;
    
    //Create the symmetrix matrix
    this->mMatrix = new SparseMatrixType(this->mDim, this->mDim);
    this->mMatrix->reserve(nonzeros);
    this->mMatrix->makeCompressed(); //Turn it into compressed row format

    std::copy(airn, airn + dim + 1, this->mMatrix->outerIndexPtr());
    std::copy(ajcn, ajcn + nonzeros, this->mMatrix->innerIndexPtr());
    std::memset(this->mMatrix->valuePtr(), 1.0, nonzeros*sizeof(double));

    //Create the Solver
    this->mSolver = new SparseSolverType();

    //We analyze the pattern of the sparse matrix
    this->mSolver->analyzePattern(*(this->mMatrix));

    return ESymSolverStatus::SYMSOLVER_SUCCESS;
  }

  ESymSolverStatus EigenSimplicialSolverInterface::MultiSolve(bool new_matrix,
    const Index* airn,
    const Index* ajcn,
    Index nrhs,
    double* rhs_vals,
    bool check_NegEVals,
    Index numberOfNegEVals)
  {

    if (new_matrix)
    {
      //Do new factorization if matrix is new
      this->mSolver->factorize(*(this->mMatrix));
    }

    for (Index i = 0; i < nrhs; i++)
    {
      Index offset = i * this->mDim;
      Eigen::Map<VectorType> currRhs = Eigen::Map<VectorType>(&(rhs_vals[offset]), this->mDim);
      
      //Solve and Store 
      currRhs = this->mSolver->solve(currRhs);      
    }

    return ESymSolverStatus::SYMSOLVER_SUCCESS;

  }

  double* EigenSimplicialSolverInterface::GetValuesArrayPtr()
  {
    return this->mMatrix->valuePtr();
  }

}