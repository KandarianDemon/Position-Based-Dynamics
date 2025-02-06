using System.Linq;
using System;
using MathNet.Numerics.LinearAlgebra;
using UnityEngine;
using System.Collections.Generic;
using UnityEditor;
using MathNet.Numerics.LinearAlgebra.Double;
using UnityEditor.PackageManager;
using VirtualWorm;

namespace Felix.Constraints{



     public abstract class Constraint{

        

        public abstract void Project(ref Vector3[] positions, bool[] pinned, int iterations,float[] weights);
        public abstract int[] ReturnIndices();
        public abstract float[] ReturnWeightCorrection(float[] weights,int indexA, int indexB);
       
    }

    class DistanceConstraint : Constraint
    {
        public int i;
        public int j;
        public float restLength;
        public float stiffness;

        public float mass;

        public DistanceConstraint(int i, int j, float restLength, float k, float m)
        {
            this.i = i;
            this.j = j;
            this.restLength = restLength;
            this.stiffness = k;
            this.mass = m;
        }

//         public override void Project(ref Vector3[] positions, bool[] pinned, int iterations, float[] weights)
// {
//     int a = this.i;
//     int b = this.j;

//     Vector3 positionA = positions[a];
//     Vector3 positionB = positions[b];

//     // Compute distance error
//     float distance = Vector3.Distance(positionA, positionB);
//     float error = distance - this.restLength;
//     Vector3 normal = (distance > 0) ? (positionA - positionB).normalized : Vector3.zero;

//     // Apply Müller's stiffness decay per iteration
//     float k_dot = 1 - Mathf.Pow((1 - stiffness), 1.0f / iterations); // Correct stiffness per iteration

//     // Compute displacement correction
//     Vector3 deltaP = error * normal * k_dot;

//     // Weight corrections (inverse mass contributions)
//     float wA = weights[a] > 0 ? 1.0f / weights[a] : 0.0f;
//     float wB = weights[b] > 0 ? 1.0f / weights[b] : 0.0f;
//     float sumW = wA + wB;
    
//     float weightA = (sumW > 0) ? wA / sumW : 0.5f;
//     float weightB = (sumW > 0) ? wB / sumW : 0.5f;

//     // Apply corrections
//     if (!pinned[a] && !pinned[b])
//     {
//         positions[a] -= weightA * deltaP;
//         positions[b] += weightB * deltaP;
//     }
//     else if (pinned[a])
//     {
//         positions[b] += deltaP;
//     }
//     else if (pinned[b])
//     {
//         positions[a] -= deltaP;
//     }
// }


        public override void Project(ref Vector3[] positions, bool[] pinned, int iterations, float[] weights)
        {
            int a = this.i;
            int b = this.j;


            //Debug.Log($"positions {positions[a]}");
            Vector3 positionA = positions[a];
            Vector3 positionB = positions[b];

            // Calculate the error and direction
            float distance = Vector3.Distance(positionA, positionB);
            float error = distance - this.restLength;
            Vector3 normal = Vector3.Normalize(positionA - positionB);





            float k_dot = 1 - Mathf.Pow((1 - stiffness), 1.0f / iterations);


            float k = stiffness;

            // Apply corrections
            Vector3 deltaP = error * normal * k;

            float inverse_mass_a = 1 / weights[a];
            float inverse_mass_b = 1 / weights[b];

            float weight_correction_a = inverse_mass_a / (inverse_mass_a + inverse_mass_b);
            float weight_correction_b = inverse_mass_b / (inverse_mass_a + inverse_mass_b);

            //float[] corrections = ReturnWeightCorrection(weights,a,b);
            float[] corrections = new float[] {weight_correction_a,weight_correction_b};
            //Debug.Log($"Weight Correction A: {weight_correction_a} Weight Correction B: {weight_correction_b} Inverse mass A: {inverse_mass_a} Inverse mass B: {inverse_mass_b}");


            if (!pinned[a] && !pinned[b])
            {
                // Split correction between both vertices
                positions[a] -= corrections[0]* deltaP * k;
                positions[b] += corrections[1] * deltaP * k;
            }
            else if (pinned[a])
            {
                // Apply the entire correction to the unpinned vertex
                positions[b] += corrections[1] * deltaP * k;
            }
            else if (pinned[b])
            {
                // Apply the entire correction to the unpinned vertex
                positions[a] -= corrections[0] * deltaP * k;
            }

            //Debug.Log($"Constraint between {a} and {b}: DeltaP={deltaP}, Pinned A={pinned[a]}, Pinned B={pinned[b]}");
        }

        public override int[] ReturnIndices()
        {
            return new int[2] { this.i, this.j };
        }

        public override float[] ReturnWeightCorrection(float[] weights, int indexA, int indexB)
        {

            float inverse_mass_a = 1 / weights[indexA];
            float inverse_mass_b = 1 / weights[indexB];

            float weight_correction_a = inverse_mass_a / (inverse_mass_a + inverse_mass_b);
            float weight_correction_b = inverse_mass_b / (inverse_mass_a + inverse_mass_b);

            return new float[2] { weight_correction_a, weight_correction_b };





        }

    }


    class BendingConstraint : Constraint {

       
        public int a;
        public int b;
        public int c;
        public int d;

        public float stiffness;
        public float initial_angle;
        public float mass;

        public BendingConstraint(int a, int b, int c,int d, float stiff,float initial_angle,float m){
            
            this.a = a;
            this.b = b;
            this.c = c;
            this.d = d;
            this.initial_angle = initial_angle;

            this.stiffness = stiff;
            this.mass = m;
        }

        public override void Project(ref Vector3[] positions, bool[] pinned, int iterations, float[] weights)
        {
            
            // for the bending constraint

            // Ensure normals point in the same direction. Constraint initialization should capture face orientation. Yes, thats whats causing issues


            //Debug.Log($"a: {this.a} b: {this.b} c: {this.c} d: {this.d} positions[a]: {positions[this.a]} positions[b]: {positions[this.b]} positions[c]: {positions[this.c]} positions[d]: {positions[this.d]}");

            Vector3 A = positions[this.a];
            Vector3 B = positions[this.b];
            Vector3 C = positions[this.c];
            Vector3 D = positions[this.d];

            Vector3 AB = positions[this.b] - positions[this.a];
            Vector3 AC = positions[this.c] - positions[this.a];
            Vector3 AD = positions[this.d] - positions[this.a];

            Vector3 n1 = Vector3.Cross(AB, AC).normalized;
            Vector3 n2 = Vector3.Cross(AC, AD).normalized;

            float clamped_dot = Mathf.Clamp(Vector3.Dot(n1, -n2), -1, 1);
            float acos = Mathf.Acos(clamped_dot);
            float angle = Mathf.Clamp(acos,-1.0f,1.0f);

            float error = angle - initial_angle;

            float weight = 1 / this.mass;
            float d = Vector3.Dot(n1, n2); // 

            // Gradients of normalized cross products. Needs to be implemented correctly as according to Mueller et al.


            // ---------------------------------------------------------------------
            Vector3 BxN2 = Vector3.Cross(B, n2);
            Vector3 N1xB = Vector3.Cross(n1, B);
            Vector3 BxC = Vector3.Cross(B, C);

            Vector3 q3 = (BxN2 + N1xB * d) / BxC.magnitude; // (p2xn2 + (n1xp2)*d)/|p2xp3|

            // ---------------------------------------------------------------------

            Vector3 BxN1 = Vector3.Cross(B, n1);
            //Vector3 BxN2 = Vector3.Cross(B,n2);
            Vector3 BxD = Vector3.Cross(B, D);




            Vector3 q4 = (BxN1 + BxN2 * d) / BxD.magnitude; // (p2xn1) + (n2xp2)*d)/|p2xp4|

            Vector3 CxN2 = Vector3.Cross(C, n2);
            Vector3 CxN1 = Vector3.Cross(C, n1);
            Vector3 CxB = Vector3.Cross(C, B);
            Vector3 DxN1 = Vector3.Cross(D, n1);
            Vector3 DxN2 = Vector3.Cross(D, n2);
            Vector3 DxB = Vector3.Cross(D, B);


            Vector3 q2 = -(CxN2 + CxN1 * d) / BxC.magnitude - (DxN1 + DxN2 * d) / BxD.magnitude; // -((p3xn2) + (n1xp3)*d)/|p2xp3| - ((p4xn1) + (n2xp4)d)/|p2xp4|

            Vector3 q1 = -q2 - q3 - q4;


            //Debug.Log($"q1: {q1} q2: {q2} q3: {q3} q4: {q4} angle: {angle} d: {d} BxN1: {BxN1} B: {B} n1: {n1}" );


             float[] inv_masses = new float[4] { 1 / weights[this.a], 1 / weights[this.b], 1 / weights[this.c], 1 / weights[this.d] };

            // final correction according to Mueller et al 2006:
            float upper = -inv_masses[0]* Mathf.Sqrt(1 - (d * d)) * error;
            float upper2 = -inv_masses[1]* Mathf.Sqrt(1 - (d * d)) * error;
            float upper3 = -inv_masses[2]* Mathf.Sqrt(1 - (d * d)) * error;
            float upper4 = -inv_masses[3]* Mathf.Sqrt(1 - (d * d)) * error;

           

            float lower =   inv_masses[1] * (q2.sqrMagnitude) + inv_masses[2] * (q3.sqrMagnitude) + inv_masses[3] * (q4.sqrMagnitude) + 1e-6f; //important: its the squared magnitude of the q's
            float lower2 = inv_masses[0]*(q1.sqrMagnitude)  + inv_masses[2] *(q3.sqrMagnitude)  + inv_masses[3]*(q4.sqrMagnitude) + 1e-6f; //important: its the squared magnitude of the q's
            float lower3 = inv_masses[0]*(q1.sqrMagnitude) + inv_masses[1] *(q2.sqrMagnitude)  +   + inv_masses[3]*(q4.sqrMagnitude) + 1e-6f; //important: its the squared magnitude of the q's
            float lower4 = inv_masses[0]*(q1.sqrMagnitude) + inv_masses[1] *(q2.sqrMagnitude)  + inv_masses[2] *(q3.sqrMagnitude) + 1e-6f; //important: its the squared magnitude of the q's
            if (lower < 1e-6f || lower2 < 1e-6f || lower3 < 1e-6f || lower4 < 1e-6f)
            {
                
                return; // Prevent NaN
            } 

            //Somewhere in this code something becomes NaN or Infinity!

            Vector3 delta_p1 = (upper / lower) * q1; // compute this for all 4 points in the constraint. // to get this we need to do partial derivatives. If these are written in a matrix we get the Jacobian
            Vector3 delta_p2 = (upper2 / lower2) * q2;
            Vector3 delta_p3 = (upper3 / lower3) * q3;
            Vector3 delta_p4 = (upper4 / lower4) * q4;

            
            // if (!pinned[this.a]) positions[this.a] += delta_p1;
            // if (!pinned[this.b]) positions[this.b] += delta_p2;
            // if (!pinned[this.c]) positions[this.c] += delta_p3;
            // if (!pinned[this.d]) positions[this.d] += delta_p4;
            
            if (!pinned[this.a]) positions[this.a] += 0.0001f*delta_p1.normalized ;
            if (!pinned[this.b]) positions[this.b] += 0.0001f*delta_p2.normalized ;
            if (!pinned[this.c]) positions[this.c] += 0.0001f*delta_p3.normalized ;
            if (!pinned[this.d]) positions[this.d] += 0.0001f*delta_p4.normalized ;




            
        }

        public override int[] ReturnIndices()
        {
            return new int[4]{this.a,this.b,this.c,this.d}; 
        }

        public override float[] ReturnWeightCorrection(float[] weights,int indexA, int indexB)
        {
            float inverse_mass_a = 1/weights[indexA];
            float inverse_mass_b = 1/weights[indexB];

            float weight_correction_a = inverse_mass_a / (inverse_mass_a + inverse_mass_b);
            float weight_correction_b = inverse_mass_b / (inverse_mass_a + inverse_mass_b);

            return new float[2]{weight_correction_a,weight_correction_b};
        }
    }

    public class SelfCollisionConstraint: Constraint
    
    {

        int q;
        int p1;
        int p2;
        int p3;
        float cloth_thickness;

        int count;

        public SelfCollisionConstraint(int q, int p1, int p2, int p3,float cloth_thickness)
        {
            this.q = q;
            this.p1 = p1;
            this.p2 = p2;
            this.p3 = p3;
            this.cloth_thickness = cloth_thickness;
            this.count = 0;

            // store reference to triangle. All the relevant data should be stored there.


        }
        public override void Project(ref Vector3[] positions, bool[] pinned, int iterations, float[] weights)
        {
           
            // C(q,p1,p2,p3) = (q-p1)  P2P1_normalized x p3p1_normalized - h     || h = cloth thickness.

            // derivative of n with respect to K,L.
            Vector3 A = positions[p2] - positions[p1];
            Vector3 B =  positions[p3] - positions[p1];

            Vector3 d = positions[q] - positions[p1];

            Vector3 n = Vector3.Cross(A,B).normalized;

            Matrix4x4 a_x = A.ToSkewSymmetricMatrix();
            Matrix4x4 b_x = B.ToSkewSymmetricMatrix();

            // float mag_AxB = Vector3.Cross(A,B).magnitude;

            // Vector3 nXb = Vector3.Cross(n, B);
            // Vector3 nXa = Vector3.Cross(n, A);

            // Matrix4x4 nnT_a = new Matrix4x4();
            // Matrix4x4 nnT_b = new Matrix4x4();

            // nnT_a = nnT_a.MultiplyVectorWithItsTranspose(n,nXa);
            // nnT_b = nnT_b.MultiplyVectorWithItsTranspose(n,nXb);

            // Matrix4x4 summand_a = b_x.MultiplyByScalar(-1).Add(nnT_a);
            // Matrix4x4 summand_b = a_x.MultiplyByScalar(-1).Add(nnT_b);

            //  float inverted_magnitude = 1 / mag_AxB;

            // Matrix4x4 K = summand_a.MultiplyByScalar(inverted_magnitude);
            // Matrix4x4 L = summand_b.MultiplyByScalar(inverted_magnitude);





           



         




            Vector3 Cq = n;
            Vector3 Cp1 = -n + (b_x.MultiplyByScalar(-1).Add(a_x).MultiplyVector(d));
            Vector3 Cp2 = b_x.MultiplyVector(d);
            Vector3 Cp3 = a_x.MultiplyByScalar(-1).MultiplyVector(d);


            // skew_symmetric_matrix of K and L

            float qp1_dot = Vector3.Dot(Cq, n);
            float error = qp1_dot - cloth_thickness;

            // 0  -Kz  Ky
            
            float sum_of_weights =  weights[p1]*Cp1.sqrMagnitude 
                                    + weights[p2]*Cp2.sqrMagnitude 
                                    + weights[p3]*Cp3.sqrMagnitude 
                                    + weights[q]*Cq.sqrMagnitude;

            float[] relative_weights = new float[] {
                
                weights[q] / (weights[p1]*Cp1.sqrMagnitude 
                                    + weights[p2]*Cp2.sqrMagnitude 
                                    + weights[p3]*Cp3.sqrMagnitude ),
                weights[p1] / (weights[p2]*Cp2.sqrMagnitude 
                                    + weights[p3]*Cp3.sqrMagnitude 
                                    + weights[q]*Cq.sqrMagnitude),
                weights[p2] / (weights[p1]*Cp1.sqrMagnitude 
                                    + weights[p3]*Cp3.sqrMagnitude 
                                    + weights[q]*Cq.sqrMagnitude),
                weights[p3] / (weights[p2]*Cp2.sqrMagnitude 
                                    + weights[p1]*Cp1.sqrMagnitude 
                                    + weights[q]*Cq.sqrMagnitude)
                
            };

            Vector3[] delta_ps = new Vector3[]{

                -relative_weights[0] * error * Cq,
                -relative_weights[1] * error * Cp1,
                -relative_weights[2] * error * Cp2,
                -relative_weights[3] * error * Cp3
            };
            
        
            positions[q] += 0.001f*delta_ps[0];
            positions[p1] += 0.001f*delta_ps[1];
            positions[p2] += 0.001f*delta_ps[2];
            positions[p3] += 0.001f*delta_ps[3];



            // 


            

        
        }

        public override int[] ReturnIndices()
        {
            // return the indices of the colliding vertices
            throw new NotImplementedException();
        }

        public override float[] ReturnWeightCorrection(float[] weights, int indexA, int indexB)
        {
            // return the weight correction for the colliding vertices
            throw new NotImplementedException();
        }

        public void Unsubscribe(ref List<Constraint> constraints)
        {
            constraints.Remove(this);
            
            // delete yourself from the constraint list once all interations have been projected.
        }
    }


    public class VolumeConstraint : Constraint{

        // This constraint ensures _ FOR A CLOSED MESH_ that the volume of the mesh remains constant.   
        // Its supposed to be used for a cloth balloon that models a pressurized C. elegans model.
        float initialVolume;
        Triangle[] mesh_triangles;

        public VolumeConstraint(float V0, Triangle[] triangles)
        {
            this.initialVolume = V0;
            this.mesh_triangles = triangles;
        }
        
        
        public override void Project(ref Vector3[] positions, bool[] pinned, int iterations, float[] weights)
        {

            foreach(var triangle in this.mesh_triangles)
            {
                float volume = triangle.ComputeScalarTripleProduct(ref positions);
                // Compute the volume of the triangle
                // Compute the gradient of the volume
                // Compute the weighting factor
                // Compute the delta p
                // Apply the delta p
            }
            


            throw new System.NotImplementedException();
        }

        public override float[] ReturnWeightCorrection(float[] weights, int indexA, int indexB)
        {
            throw new System.NotImplementedException();
        }

        public override int[] ReturnIndices()
        {
            throw new System.NotImplementedException();
        }
    }
}

