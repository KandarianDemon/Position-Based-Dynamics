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
    [Tooltip("Type of mesh to generate")]
    public MeshType meshType;

    [Header("Debbuging Options")]
    public bool icosphere_equator = false;
    //public float angle_x = 0;
    //public float angle_y = 0;
    //public float angle_z = 0;
    [Tooltip("Shows vertices and their indices")]
    public bool show_vertices = false;

    [Tooltip("Shows the vertices of the equator of the sphere")]
    public bool show_equator = false;

    [Tooltip("position of the x,y-plane in z direction, so that the equator vertices in that plane can be visualized")]
    public float tolerance = 0.01f;
    [Range(1, 3)]
    [Tooltip("Number of subdivisions for the icosphere")]
    public int subdivisions = 1;
    [Header("Worm Creation Settings")]
    [Tooltip("Length of the worm")]
    public float worm_length = 25.0f;
    [Range(0.01f, 2.0f)]
    [Tooltip("Bulge factor for the worm")]
    public float bulgeFactor = 0.5f;
    [Range(0.2f, 3.0f)]
    [Tooltip("Elongation factor for the worms end caps")]
    public float elongation = 1.0f;

    public Material meshMaterial;
    

    [Header("Rotation Settings")]
    public Vector3 rotation_angle = new Vector3(0, 0, 0);
    // Start is called before the first frame update
    void Awake()
    {
        // Initialization code here
        filter = transform.GetComponent<MeshFilter>();
        renderer = transform.GetComponent<MeshRenderer>();
        softbody = transform.GetComponent<Softbody>();

        InitializeMesh();

        renderer.material = meshMaterial;
        renderer.material.color = Color.cyan ;

        //Debug.Log($"filter mesh has  {filter.mesh.triangles.Length / 3} triangles");

    }

    void InitializeMesh()
    {
        switch (meshType)
        {
            case MeshType.Plane:
                filter.mesh = MeshGenerator.GeneratePlane();
                this.transform.name = "Procedural Plane";
                break;
            case MeshType.Cube:
                filter.mesh = MeshGenerator.GenerateCube();
                this.transform.name = "Procedural Cube";
                break;
            case MeshType.Sphere:
                filter.mesh = MeshGenerator.GenerateSphere();
                this.transform.name = "Procedural Sphere";
                break;
            case MeshType.Tetrahedron:
                filter.mesh = MeshGenerator.GenerateTetrahedron();
                this.transform.name = "Procedural Tetrahedron";
                break;
            case MeshType.Icosahedron:
                filter.mesh = MeshGenerator.GenerateIcosahedron();
                this.transform.name = "Procedural Icosahedron";
                break;

            case MeshType.Icosphere:
                filter.mesh = MeshGenerator.GenerateIcoSphere(this.subdivisions);
                this.transform.name = "Procedural Icosphere";
                break;

            case MeshType.Worm:
                filter.mesh = MeshGenerator.GenerateWorm(subdivisions: this.subdivisions,length:this.worm_length, bulgeFactor:this.bulgeFactor, elongation:this.elongation);
                this.transform.name = "Procedural Worm";
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
            m11 = Mathf.Cos(rotation_angle.x * Mathf.Deg2Rad),
            m12 = -Mathf.Sin(rotation_angle.x * Mathf.Deg2Rad),
            m21 = Mathf.Sin(rotation_angle.x * Mathf.Deg2Rad),
            m22 = Mathf.Cos(rotation_angle.x * Mathf.Deg2Rad),
            m33 = 1
        };

        Matrix4x4 r_y = new Matrix4x4
        {
            m00 = Mathf.Cos(rotation_angle.y * Mathf.Deg2Rad),
            m02 = Mathf.Sin(rotation_angle.y * Mathf.Deg2Rad),
            m11 = 1,
            m20 = -Mathf.Sin(rotation_angle.y * Mathf.Deg2Rad),
            m22 = Mathf.Cos(rotation_angle.y * Mathf.Deg2Rad),
            m33 = 1
        };

        Matrix4x4 r_z = new Matrix4x4
        {
            m00 = Mathf.Cos(rotation_angle.z * Mathf.Deg2Rad),
            m01 = -Mathf.Sin(rotation_angle.z * Mathf.Deg2Rad),
            m10 = Mathf.Sin(rotation_angle.z * Mathf.Deg2Rad),
            m11 = Mathf.Cos(rotation_angle.z * Mathf.Deg2Rad),
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