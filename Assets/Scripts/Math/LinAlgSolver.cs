using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace Felix
{
    namespace LinAlgSolvers
    {
 
      public  class LinAlgSolver{

        public Matrix<double> A;
        public LinAlgSolver()
        {
            A = RandomMatrix();
        }

        public void SolveLinearEquations()
        {
            throw new System.NotImplementedException();
        }

        public Matrix<double> RandomMatrix()
        {
            Matrix<double> A = DenseMatrix.OfArray(new double[,] {
            {1, 2, 3},
            {4, 5, 6},
            {7, 8, 9}
        });
            return A;
        }

        }


    }
}