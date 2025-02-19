using UnityEngine;
using System.Collections.Generic;
using Felix.Constraints;



public class PBDSolver 
{
    int solver_iterations;
    int[] constraintPositionIndices;
    int[] constraintTypes;
    int[] numberOfPosisionsInConstraint;
    public int _Solver
    {
        get { return solver_iterations; }
        set { solver_iterations = value; }
    }

    public PBDSolver(int solver_iterations)
    {
        this.solver_iterations = solver_iterations;
    }

    public void InitializeConstraints()
    {

        // loop over edges
        // generate triangles
        
        // also initialize diagonal distance constraints
        // triangleEdgeMap

        throw new System.NotImplementedException();
    }

    public void Solve(Constraint[] constraints, ref Vector3[] positions)
    {

        for (int i = 0; i < solver_iterations; i++)
        {

        }
    }
}