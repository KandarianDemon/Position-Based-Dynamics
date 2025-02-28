using System.Collections;
using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra;
using Unity.VisualScripting;
using UnityEngine;

public static class ExtensionMethods 
{
    // Start is called before the first frame update

     public static Vector3 DampenVelocity(this Vector3 v, Vector3 n, float d)
    {
        float normal_component = Vector3.Dot(v, n);

        return v - (1 - d) * normal_component * n;
    }


    public static Vector3 Reflect(this Vector3 v, Vector3 n)
    {
        return v - 2 * Vector3.Dot(v, n) * n;
    }


   
    public static Matrix4x4 ToSkewSymmetricMatrix(this Vector3 v)
    {
        return new Matrix4x4(
            new Vector4(0, -v.z, v.y, 0),
            new Vector4(v.z, 0, -v.x, 0),
            new Vector4(-v.y, v.x, 0, 0),
            new Vector4(0, 0, 0, 0)
        );
    }
    public static Vector3[] DistancesToCenter(this Vector3 v,Vector3 cm, Vector3[] positions)
    {
        Vector3[] distances = new Vector3[positions.Length];
        for(int i = 0; i<positions.Length; i++)
        {
            distances[i] = positions[i] - cm;
        }
        return distances;
    }

    public static Matrix4x4 MultiplyByScalar(this Matrix4x4 m, float scalar)
    {
        return new Matrix4x4(
            m.GetRow(0) * scalar,
            m.GetRow(1) * scalar,
            m.GetRow(2) * scalar,
            m.GetRow(3) * scalar
        );
    }

  
    public static Matrix4x4 Add(this Matrix4x4 a, Matrix4x4 b)
    {
        return new Matrix4x4(
            a.GetRow(0) + b.GetRow(0),
            a.GetRow(1) + b.GetRow(1),
            a.GetRow(2) + b.GetRow(2),
            a.GetRow(3) + b.GetRow(3)
        );
    }

    public static Matrix4x4 MultiplyVectorWithItsTranspose(this Matrix4x4 a, Vector3 b, Vector3 c)
    {
        return new Matrix4x4(
            new Vector4(b.x * c.x, b.x * c.y, b.x * c.z, 0),
            new Vector4(b.y * c.x, b.y * c.y, b.y * c.z, 0),
            new Vector4(b.z * c.x, b.z * c.y, b.z * c.z, 0),
            new Vector4(0, 0, 0, 0)
        );
    }


    public static Matrix4x4 GetInertiaTensor(this Matrix4x4 m,float[] weights,Vector3[] distances)
    {
        // get R_tilde for each point.
        // sum of skew symmetric matrix and the transpose times weight(i)
        Matrix4x4 inertiaTensor = new Matrix4x4();
        for(int i = 0; i<weights.Length; i++)
        {
            Matrix4x4 R_tilde = distances[i].ToSkewSymmetricMatrix();
            Matrix4x4 R_tilde_T = R_tilde.transpose;

            Matrix4x4 mult = R_tilde * R_tilde_T;
            mult = mult.MultiplyByScalar(weights[i]);

            inertiaTensor.Add(mult);
            

        }


        return inertiaTensor;
    }

    public static void PrintMatrix(this Matrix4x4 m)
    {
        Debug.Log(m.GetRow(0));
        Debug.Log(m.GetRow(1));
        Debug.Log(m.GetRow(2));
        Debug.Log(m.GetRow(3));
    }
}
