using UnityEngine;
using System.Collections.Generic;

namespace Felix{



public class HookeanSpring : MonoBehaviour
{
    public Transform pointA;
    public Transform pointB;
    public float springConstant = 1.0f;
    public float restLength = 1.0f;

    void Update()
    {
        ProjectSpring();
    }

    public void OnDrawGizmos(){

        Gizmos.color = Color.red;
        Gizmos.DrawLine(pointA.position, pointB.position);

    }


    void ProjectSpring(){

         Vector3 displacement = pointB.position - pointA.position;
        float currentLength = displacement.magnitude;
        Vector3 direction = displacement.normalized;
        float extension = currentLength - restLength;
        Vector3 force = springConstant * extension * direction;

        // Apply the force to the objects (assuming they have Rigidbody components)
        if (pointA.GetComponent<Rigidbody>() != null)
        {
            pointA.GetComponent<Rigidbody>().AddForce(force);
        }
        if (pointB.GetComponent<Rigidbody>() != null)
        {
            pointB.GetComponent<Rigidbody>().AddForce(-force);
        }

    }
}
}