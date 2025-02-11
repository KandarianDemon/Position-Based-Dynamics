using UnityEngine;
using System.Collections.Generic;
using Felix.MeshGeneration;
using VirtualWorm;
using UnityEditor;

[RequireComponent(typeof(MeshFilter), typeof(MeshRenderer))]
public class TestMesh : MonoBehaviour
{
    MeshFilter filter;
    MeshRenderer renderer;
    Softbody softbody;
    public MeshType meshType;
    public bool icosphere_equator = false;
    public float angle_x = 0;
    public float angle_y = 0;
    public float angle_z = 0;
    public bool show_vertices = false;
    public bool show_equator = false;
    public float tolerance = 0.01f;
    [Range(1,3)]
    public int subdivisions = 1;
    // Start is called before the first frame update
    void Awake()
    {
        // Initialization code here
        filter = transform.GetComponent<MeshFilter>();
        renderer = transform.GetComponent<MeshRenderer>();
        softbody = transform.GetComponent<Softbody>();

        InitializeMesh();

        Debug.Log($"filter mesh has  {filter.mesh.triangles.Length / 3} triangles");

    }

    void InitializeMesh()
    {
        switch (meshType)
        {
            case MeshType.Plane:
                filter.mesh = MeshGenerator.GeneratePlane();
                break;
            case MeshType.Cube:
                filter.mesh = MeshGenerator.GenerateCube();
                break;
            case MeshType.Sphere:
                filter.mesh = MeshGenerator.GenerateSphere();
                break;
            case MeshType.Tetrahedron:
                filter.mesh = MeshGenerator.GenerateTetrahedron();
                break;
            case MeshType.Icosahedron:
                filter.mesh = MeshGenerator.GenerateIcosahedron();
                break;

            case MeshType.Icosphere:
                filter.mesh = MeshGenerator.GenerateIcoSphere(this.subdivisions);
                break;

            case MeshType.Worm:
                filter.mesh = MeshGenerator.GenerateWorm(subdivisions: this.subdivisions,length:25.0f);
                break;
        }

    }

    private void OnDrawGizmos()
    {

        Gizmos.color = Color.red;

        if (filter == null) return;
        if (show_vertices)
        {
            for (int i = 0; i < filter.mesh.vertices.Length; i++)
            {

                Gizmos.color = (Mathf.Abs(filter.mesh.vertices[i].z) < 0.01f) ? Color.green : Color.red;
                Gizmos.DrawSphere(this.transform.TransformPoint(filter.mesh.vertices[i]), 0.01f);
            }
        }

        if (show_equator)
        {
            Gizmos.color = Color.green;
            List<int> equator = new List<int>();
            List<Vector3> equatorVerts = new List<Vector3>();
            for (int i = 0; i < filter.mesh.vertices.Length; i++)
            {
                if (Mathf.Abs(filter.mesh.vertices[i].y - tolerance) < 0.3f)
                {
                    equator.Add(i);
                    equatorVerts.Add(filter.mesh.vertices[i]);

                }
            }

            equator.Sort((a, b) =>
            {
                Debug.Log($"a: {a} b: {b}");
                float angleA = Mathf.Atan2(equatorVerts[a].z, equatorVerts[a].x);
                float angleB = Mathf.Atan2(equatorVerts[b].z, equatorVerts[b].x);
                return angleA.CompareTo(angleB);
            });

            for (int i = 0; i < equator.Count; i++)
            {
                Gizmos.DrawSphere(this.transform.TransformPoint(filter.mesh.vertices[equator[i]]), 0.1f);
                Handles.Label(this.transform.TransformPoint(filter.mesh.vertices[equator[i]]), $"     {equator[i]}");
            }


        }
    }
    void OnValidate()
    {
        InitializeMesh();
        softbody.InitializeSoftbody();
        RotateMesh();

    }

    public void RotateMesh()
    {




        Matrix4x4 r_x = new Matrix4x4
        {
            m00 = 1,
            m11 = Mathf.Cos(angle_x * Mathf.Deg2Rad),
            m12 = -Mathf.Sin(angle_x * Mathf.Deg2Rad),
            m21 = Mathf.Sin(angle_x * Mathf.Deg2Rad),
            m22 = Mathf.Cos(angle_x * Mathf.Deg2Rad),
            m33 = 1
        };

        Matrix4x4 r_y = new Matrix4x4
        {
            m00 = Mathf.Cos(angle_y * Mathf.Deg2Rad),
            m02 = Mathf.Sin(angle_y * Mathf.Deg2Rad),
            m11 = 1,
            m20 = -Mathf.Sin(angle_y * Mathf.Deg2Rad),
            m22 = Mathf.Cos(angle_y * Mathf.Deg2Rad),
            m33 = 1
        };

        Matrix4x4 r_z = new Matrix4x4
        {
            m00 = Mathf.Cos(angle_z * Mathf.Deg2Rad),
            m01 = -Mathf.Sin(angle_z * Mathf.Deg2Rad),
            m10 = Mathf.Sin(angle_z * Mathf.Deg2Rad),
            m11 = Mathf.Cos(angle_z * Mathf.Deg2Rad),
            m22 = 1,
            m33 = 1
        };

        Matrix4x4 rotationMatrix = r_x * r_y * r_z;

        Vector3[] vertices = filter.mesh.vertices;

        for (int i = 0; i < vertices.Length; i++)
        {
            vertices[i] = rotationMatrix.MultiplyPoint(vertices[i]);
        }

        filter.mesh.vertices = vertices;
        Debug.Log("Rotated mesh");
        
    }

    // Update is called once per frame

}