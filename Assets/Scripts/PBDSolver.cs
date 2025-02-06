using UnityEngine;
using System.Collections.Generic;
using Felix.Constraints;



public class PBDSolver 
{
    int solver_iterations;
    public int _Solver
    {
        get { return solver_iterations; }
        set { solver_iterations = value; }
    }

    public PBDSolver(int solver_iterations)
    {
        this.solver_iterations = solver_iterations;
    }

    public void Solve(Constraint[] constraints){

        for (int i = 0; i < solver_iterations; i++)
        {
            foreach (Constraint constraint in constraints)
            {
                //constraint.Project();
            }
        }
    }
}