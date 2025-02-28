using UnityEngine;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;
using TMPro;
using UnityEngine.UIElements;
using System.Linq;
using Unity.VisualScripting.Dependencies.Sqlite;

namespace Felix.MeshGeneration
{

    public enum MeshType { Plane, Cube, Sphere, Tetrahedron, Icosahedron,Icosphere,Worm }
    public static class MeshGenerator
    {





        public static Mesh GeneratePlane()
        {
            Mesh mesh = new Mesh();

            Vector3[] vertices = new Vector3[]
            {
            new Vector3(0, 0, 0),
            new Vector3(1, 0, 0),
            new Vector3(0, 0, 1),
            new Vector3(1, 0, 1)
            };

            int[] triangles = new int[]
            {
            0, 2, 1,
            2, 3, 1
            };

            Vector2[] uv = new Vector2[]
            {
            new Vector2(0, 0),
            new Vector2(1, 0),
            new Vector2(0, 1),
            new Vector2(1, 1)
            };

            mesh.vertices = vertices;
            mesh.triangles = triangles;
            mesh.uv = uv;
            mesh.RecalculateNormals();

            return mesh;
        }

        public static Mesh GenerateCube()
        {
            Mesh mesh = new Mesh();

            Vector3[] vertices = new Vector3[]
            {
                // Front face
            new Vector3(0,0, 1),
            new Vector3(-1,1,0),
            new Vector3(1,1,0),
            new Vector3(1,-1,0),
            new Vector3(-1,-1,0),
            new Vector3(0,0,-1)
            };

            int[] triangles = new int[]
            {

            0, 1, 2,
            0, 3, 2,
            0, 3, 4,
            0, 1, 4,
            5, 1, 2,
            5, 3, 2,
            5, 3, 4,
            5, 1, 4





            };

            mesh.vertices = vertices;
            mesh.triangles = triangles;
            mesh.RecalculateNormals();

            return mesh;
        }

        public static Mesh GenerateSphere(int resolution = 20)
        {
            Mesh mesh = new Mesh { name = "Procedrual Sphere" };


            List<Vector3> vertices = new List<Vector3>();
            List<int> triangles = new List<int>();

            // Top pole
            vertices.Add(new Vector3(0, 1, 0)); // Index 0 (top)

            // Generate middle vertices
            for (int lat = 1; lat < resolution; lat++)
            {
                float theta = Mathf.PI * lat / resolution;
                float sinTheta = Mathf.Sin(theta);
                float cosTheta = Mathf.Cos(theta);

                for (int lon = 0; lon < resolution; lon++)
                {
                    float phi = 2 * Mathf.PI * lon / resolution;
                    float x = sinTheta * Mathf.Cos(phi);
                    float y = cosTheta;
                    float z = sinTheta * Mathf.Sin(phi);

                    vertices.Add(new Vector3(x, y, z));
                }
            }

            // Bottom pole
            vertices.Add(new Vector3(0, -1, 0)); // Last index

            int bottomPoleIndex = vertices.Count - 1;

            // Generate triangles for the top cap
            for (int lon = 0; lon < resolution; lon++)
            {
                int nextLon = (lon + 1) % resolution;
                triangles.Add(0);
                triangles.Add(1 + nextLon);
                triangles.Add(1 + lon);
            }

            // Generate triangles for the middle part
            for (int lat = 0; lat < resolution - 2; lat++)
            {
                for (int lon = 0; lon < resolution; lon++)
                {
                    int current = 1 + lat * resolution + lon;
                    int next = current + resolution;
                    int nextLon = (lon + 1) % resolution;

                    triangles.Add(current);
                    triangles.Add(next);
                    triangles.Add(current + nextLon);

                    triangles.Add(current + nextLon);
                    triangles.Add(next);
                    triangles.Add(next + nextLon);
                }
            }

            // Generate triangles for the bottom cap
            int lastRowStart = 1 + (resolution - 2) * resolution;
            for (int lon = 0; lon < resolution; lon++)
            {
                int nextLon = (lon + 1) % resolution;
                triangles.Add(bottomPoleIndex);
                triangles.Add(lastRowStart + lon);
                triangles.Add(lastRowStart + nextLon);
            }

            // Assign mesh data
            mesh.vertices = vertices.ToArray();
            mesh.triangles = triangles.ToArray(); // ERROR HAPPENS HERE IF INDICES ARE WRONG
            mesh.RecalculateNormals();
            mesh.RecalculateBounds();

            return mesh;
        }

        public static Mesh GenerateTetrahedron()
        {
            Mesh mesh = new Mesh { name = "Tetrahedron" };

            Vector3[] vertices = new Vector3[]
            {
            new Vector3(0, 0, 1.5f),
            new Vector3(-1, 1, 0),
            new Vector3(1, 1, 0),
            new Vector3(1, -1, 0),
            new Vector3(-1, -1, 0),
            new Vector3(0,0,-1.5f)
            };

            int[] triangles = new int[]{
            0, 2,1,
            0, 3, 2,
            0,4,3,
            0,1,4,

             5, 1, 2,
            5, 2, 3,
            5,3,4,
            5,4,1,


            };

            mesh.vertices = vertices;
            mesh.triangles = triangles;
            mesh.RecalculateNormals();


            return mesh;
        }

        public static Mesh GenerateIcosahedron()
        {
            Mesh mesh = new Mesh { name = "Icosahedron" };

            // 12 vertices
            // 3 Orthogonal golden rectangles -> keeping the golden ratio

            // golden ratio

            float t = (1.0f + Mathf.Sqrt(5.0f)) / 2.0f;
            Vector3[] vertices = new Vector3[]
            {
                // Rectangle 1
                new Vector3(-1, t, 0),
                new Vector3(1, t, 0),
                new Vector3(-1, -t, 0),
                new Vector3(1, -t, 0),


                // Rectangle 2
                new Vector3(0, -1, t),
                new Vector3(0, 1, t),
                new Vector3(0, -1, -t),
                new Vector3(0, 1, -t),

                // Rectangle 3
                new Vector3(t, 0, -1),
                new Vector3(t, 0, 1),
                new Vector3(-t, 0, -1),
                new Vector3(-t, 0, 1)



            };

            float angle = -Vector3.Angle(Vector3.up, vertices[0]) * Mathf.Deg2Rad;
            Matrix4x4 r_z = new Matrix4x4
            {


                m00 = Mathf.Cos(angle),
                m01 = -Mathf.Sin(angle),
                m02 = 0,
                m03 = 0,
                m10 = Mathf.Sin(angle),
                m11 = Mathf.Cos(angle),
                m12 = 0,
                m13 = 0,
                m20 = 0,
                m21 = 0,
                m22 = 1,
                m23 = 0,
                m30 = 0,
                m31 = 0,
                m32 = 0,
                m33 = 1
            };

            Matrix4x4 r_x = new Matrix4x4
            {
                m00 = 1,
                m01 = 0,
                m02 = 0,
                m03 = 0,
                m10 = 0,
                m11 = Mathf.Cos(angle),
                m12 = -Mathf.Sin(angle),
                m13 = 0,
                m20 = 0,
                m21 = Mathf.Sin(angle),
                m22 = Mathf.Cos(angle),
                m23 = 0,
                m30 = 0,
                m31 = 0,
                m32 = 0,
                m33 = 1
            };



            for (int i = 0; i < vertices.Length; i++)
            {

                //vertices[i] = vertices[i].normalized;
                Debug.Log("Before: " + vertices[i]);
                vertices[i] = r_z.MultiplyPoint3x4(vertices[i]);
                Debug.Log("After: " + vertices[i]);

            }

            int[] triangles = new int[]
{
    0, 11, 5,   0, 5, 1,   0, 1, 7,   0, 7, 10,   0, 10, 11,
    1, 5, 9,    5, 11, 4,  11, 10, 2,  10, 7, 6,   7, 1, 8,
    3, 9, 4,    3, 4, 2,   3, 2, 6,   3, 6, 8,    3, 8, 9,
    4, 9, 5,    2, 4, 11,  6, 2, 10,  8, 6, 7,    9, 8, 1
};



            mesh.vertices = vertices;
            mesh.triangles = triangles;

            mesh.RecalculateBounds();
            mesh.RecalculateNormals();

            return mesh;
        }

        public static Mesh GenerateIcoSphere(int subdivisions = 1)
        {


            Mesh mesh = new Mesh { name = "Icosphere" };

            // Golden ratio
            float t = (1f + Mathf.Sqrt(5f)) / 2f;

            // Create 12 vertices of a icosahedron
            List<Vector3> vertices = new List<Vector3>
            {
                new Vector3(-1,  t,  0).normalized,
                new Vector3( 1,  t,  0).normalized,
                new Vector3(-1, -t,  0).normalized,
                new Vector3( 1, -t,  0).normalized,
                new Vector3( 0, -1,  t).normalized,
                new Vector3( 0,  1,  t).normalized,
                new Vector3( 0, -1, -t).normalized,
                new Vector3( 0,  1, -t).normalized,
                new Vector3( t,  0, -1).normalized,
                new Vector3( t,  0,  1).normalized,
                new Vector3(-t,  0, -1).normalized,
                new Vector3(-t,  0,  1).normalized
            };
            float angle = Vector3.Angle(Vector3.up, vertices[0]) * Mathf.Deg2Rad;


            Matrix4x4 r_z = new Matrix4x4
            {


                m00 = Mathf.Cos(angle),
                m01 = -Mathf.Sin(angle),
                m02 = 0,
                m03 = 0,
                m10 = Mathf.Sin(angle),
                m11 = Mathf.Cos(angle),
                m12 = 0,
                m13 = 0,
                m20 = 0,
                m21 = 0,
                m22 = 1,
                m23 = 0,
                m30 = 0,
                m31 = 0,
                m32 = 0,
                m33 = 1
            };

            //List<Vector3> vertices = GenerateIcosahedron().vertices.ToList();

            //rotate vertices around z 45 degrees



            for (int i = 0; i < vertices.Count; i++)
            {
                vertices[i] = r_z.MultiplyPoint3x4(vertices[i]);
            }

            // Create 20 triangles of the icosahedron
            List<int> triangles = new List<int>
    {
        0,  11, 5,
        0,  5,  1,
        0,  1,  7,
        0,  7,  10,
        0,  10, 11,
        1,  5,  9,
        5,  11, 4,
        11, 10, 2,
        10, 7,  6,
        7,  1,  8,
        3,  9,  4,
        3,  4,  2,
        3,  2,  6,
        3,  6,  8,
        3,  8,  9,
        4,  9,  5,
        2,  4,  11,
        6,  2,  10,
        8,  6,  7,
        9,  8,  1
    };

            // Refine triangles
            for (int i = 0; i < subdivisions; i++)
            {
                var refinedTriangles = new List<int>();
                Dictionary<long, int> middlePointIndexCache = new Dictionary<long, int>();

                int GetMiddlePoint(int p1, int p2)
                {
                    long smallerIndex = Mathf.Min(p1, p2);
                    long greaterIndex = Mathf.Max(p1, p2);
                    long key = (smallerIndex << 32) + greaterIndex;

                    if (middlePointIndexCache.TryGetValue(key, out int ret))
                        return ret;

                    Vector3 point1 = vertices[p1];
                    Vector3 point2 = vertices[p2];
                    Vector3 middle = ((point1 + point2) / 2f).normalized;

                    vertices.Add(middle);
                    middlePointIndexCache.Add(key, vertices.Count - 1);
                    return vertices.Count - 1;
                }

                for (int j = 0; j < triangles.Count; j += 3)
                {
                    int v1 = triangles[j];
                    int v2 = triangles[j + 1];
                    int v3 = triangles[j + 2];

                    int a = GetMiddlePoint(v1, v2);
                    int b = GetMiddlePoint(v2, v3);
                    int c = GetMiddlePoint(v3, v1);

                    refinedTriangles.AddRange(new[] { v1, a, c });
                    refinedTriangles.AddRange(new[] { v2, b, a });
                    refinedTriangles.AddRange(new[] { v3, c, b });
                    refinedTriangles.AddRange(new[] { a, b, c });
                }

                triangles = refinedTriangles;
            }

            mesh.vertices = vertices.ToArray();
            mesh.triangles = triangles.ToArray();
            mesh.RecalculateNormals();
            mesh.RecalculateBounds();

            return mesh;
        }


        public static Mesh GenerateWorm(float radius = 1f, float length = 10f, int subdivisions = 1, int cylinderSegments = 28, float bulgeFactor = 0.5f, float elongation = 2.0f)
        {
            Mesh wormMesh = new Mesh();
            List<Vector3> vertices = new List<Vector3>();
            List<int> triangles = new List<int>();

            // First generate a complete icosphere to get the ring of vertices at the cut
            var baseSphere = GenerateIcoSphere(subdivisions);

            List<Vector3> sphereVerts = new List<Vector3>(baseSphere.vertices);
            List<int> sphereTriangles = new List<int>(baseSphere.triangles);

            // Find vertices at the "equator" (where we'll cut and connect to cylinder)
            float tolerance = 0.01f;
            List<int> equatorIndices = new List<int>();
            HashSet<int> equatorIndexSet = new HashSet<int>();

            // Find all vertices that lie on the XY plane (z = 0)
            for (int i = 0; i < sphereVerts.Count; i++)
            {
                if (Mathf.Abs(sphereVerts[i].y) < tolerance)
                {
                    equatorIndices.Add(i);
                    equatorIndexSet.Add(i);
                }
            }

            // Sort equator vertices by angle around the z-axis !!
            equatorIndices.Sort((a, b) =>
            {
                float angleA = Mathf.Atan2(sphereVerts[a].z, sphereVerts[a].x);
                float angleB = Mathf.Atan2(sphereVerts[b].z, sphereVerts[b].x);
                return angleA.CompareTo(angleB);
            });

            int vertsPerRing = equatorIndices.Count;

            // Generate front hemisphere
            List<Vector3> frontVerts = new List<Vector3>();
            List<int> frontTris = new List<int>();
            Dictionary<int, int> frontVertexMap = new Dictionary<int, int>();

            // Add all vertices that have z >= 0
            for (int i = 0; i < sphereVerts.Count; i++)
            {
                if (sphereVerts[i].y <= tolerance)
                {
                    frontVertexMap[i] = frontVerts.Count;
                    frontVerts.Add(new Vector3(sphereVerts[i].x * radius, sphereVerts[i].y * elongation * radius, sphereVerts[i].z * radius));
                }
            }


            // Add triangles that use these vertices
            for (int i = 0; i < sphereTriangles.Count; i += 3)
            {
                int v1 = sphereTriangles[i];
                int v2 = sphereTriangles[i + 1];
                int v3 = sphereTriangles[i + 2];

                if (frontVertexMap.ContainsKey(v1) &&
                    frontVertexMap.ContainsKey(v2) &&
                    frontVertexMap.ContainsKey(v3))
                {
                    frontTris.Add(frontVertexMap[v1]);
                    frontTris.Add(frontVertexMap[v2]);
                    frontTris.Add(frontVertexMap[v3]);
                }
            }




            // Generate back hemisphere
            List<Vector3> backVerts = new List<Vector3>();
            List<int> backTris = new List<int>(frontTris);


            for (int i = 0; i < frontVerts.Count; i++)
            {
                backVerts.Add(new Vector3(frontVerts[i].x, (frontVerts[i].y) - length, frontVerts[i].z) * -1.0f);
            }

            for (int i = 0; i < backTris.Count; i += 3)
            {
                int a = backTris[i];
                int b = backTris[i + 1];
                int c = backTris[i + 2];

                backTris[i] = a;
                backTris[i + 1] = c;
                backTris[i + 2] = b;
            }

            // Generate cylinder with matching resolution
            List<Vector3> cylinderVerts = new List<Vector3>();
            List<int> cylinderTris = new List<int>();
            // cylindersegments: get distance between equator vertices, divide through length

            Dictionary<int, int> cylinderVertexMap = new Dictionary<int, int>();

            float dst = (baseSphere.vertices[equatorIndices[1]] - baseSphere.vertices[equatorIndices[0]]).magnitude;
            cylinderSegments = Mathf.FloorToInt(length / dst);


            int ringsCount = cylinderSegments + 1;

            // Generate rings of vertices
            for (int ring = 0; ring < ringsCount; ring++)
            {
                float z = (float)ring / (cylinderSegments) * length;

                // Calculate the bulging effect (using sine wave function)
                //float bulgeFactor = 0.5f; // Adjust this to control the amount of bulge
                float scale = 1 + bulgeFactor * Mathf.Sin(Mathf.PI * (ring / (float)(ringsCount - 1)));

                for (int i = 0; i < vertsPerRing; i++)
                {
                    // Use the exact positions from the equator vertices for first and last rings
                    if (ring == 0)
                    {
                        Vector3 sphereVert = sphereVerts[equatorIndices[i]] * radius;
                        cylinderVerts.Add(sphereVert * radius);

                        //cylinderVertexMap[cylinderVerts.Count - 1] = equatorIndices[i];
                    }
                    else if (ring == ringsCount - 1)
                    {
                        Vector3 sphereVert = sphereVerts[equatorIndices[i]] * radius;
                        // Apply the scaling to make the last ring fit with the bulged center
                        cylinderVerts.Add(new Vector3(sphereVert.x * scale, z, sphereVert.z * scale));
                        //cylinderVerts.Add(backVerts[i] * radius);
                        //cylinderVertexMap[cylinderVerts.Count - 1] = equatorIndices[i];
                    }
                    else
                    {
                        // Interpolate between rings and apply the bulge scaling
                        Vector3 baseVert = sphereVerts[equatorIndices[i]] * radius;
                        cylinderVerts.Add(new Vector3(baseVert.x * scale, z, baseVert.z * scale));

                    }
                }
            }

            Debug.Log($"cylinderVerts Mappings: {cylinderVertexMap.Count} verts per ring: {vertsPerRing} cylinderSegments: {cylinderSegments}");
            foreach (var key in cylinderVertexMap.Keys)
            {
                Debug.Log($"key: {key} value: {cylinderVertexMap[key]}");
            }


            // Generate cylinder triangles
            for (int ring = 0; ring < cylinderSegments; ring++)
            {
                int ringStartIdx = ring * vertsPerRing;
                int nextRingStartIdx = (ring + 1) * vertsPerRing;



                for (int i = 0; i < vertsPerRing; i++)
                {
                    int nextI = (i + 1) % vertsPerRing;

                    cylinderTris.Add(ringStartIdx + i);
                    cylinderTris.Add(nextRingStartIdx + i);
                    cylinderTris.Add(ringStartIdx + nextI);


                    cylinderTris.Add(ringStartIdx + nextI);
                    cylinderTris.Add(nextRingStartIdx + i);
                    cylinderTris.Add(nextRingStartIdx + nextI);

                }
            }

            // frontVertexMap = MapMeshEndings(equatorIndices.ToArray(), sphereVerts, cylinderTris, ref cylinderVerts, vertsPerRing);
            // ReplaceIndices(ref frontTris, frontVertexMap);

            // frontVertexMap = MapMeshEndings(equatorIndices.ToArray(), sphereVerts, cylinderTris, ref cylinderVerts, vertsPerRing);
            // ReplaceIndices(ref backTris, frontVertexMap);

            // Get first and last ring. connect them with the cylinder quads










            // Combine all meshes
            int frontOffset = 0;
            for (int i = 0; i < frontVerts.Count; i++)
            {
                vertices.Add(frontVerts[i]);
            }



            int cylinderOffset = vertices.Count;

            for (int i = 0; i < cylinderVerts.Count - 1; i++)
            {
                // the first vertices of the cylinder are actually already in the array, so we need to update their indices to create the triangles.

                if (i < vertsPerRing)
                {
                    continue;
                }

                if (i >= cylinderVerts.Count - vertsPerRing)
                {
                    continue;
                }

                vertices.Add(cylinderVerts[i]);
            }

            int backOffset = vertices.Count;

            for (int i = 0; i < backVerts.Count; i++)
            {

                vertices.Add(backVerts[i]);
            }




            // Add all triangles with correct offsets
            //triangles.AddRange(frontTris);

            // triangle creation a bit tricky

            //DebugMesh(frontTris, frontVerts, Color.red);
            //DebugMesh(cylinderTris, cylinderVerts, Color.green);
            //DebugMesh(backTris, backVerts, Color.blue);





            for (int i = 0; i < frontTris.Count; i++)
            {
                triangles.Add(frontTris[i] + frontOffset);
            }



            //Debug.Log($"Cylinder limit {cylinderTris.Count / 2 + 105 * vertsPerRing} cylinder tri count {cylinderTris.Count} verts per ring {vertsPerRing} sphereverts {sphereVerts.Count} || {frontTris.Count} cylinderoffset: {cylinderOffset}");
            for (int i = 0; i < (cylinderTris.Count - vertsPerRing * 12); i++)
            {
                // need to stich the cylinder to the front equator verts

                // look up vert in cylinder verts, if it is in front verts, get index of the front verts and use that index.

                // check if vert is in front verts, if so get index of front verts



                // if index is in double verts, use the smaller index.



                int index = cylinderTris[i] + cylinderOffset; // index for vertex look up

                // stitch first ring
                // check for each vert if its in the front verts


                // stitch last ring








                triangles.Add(index);

                //triangles.Add(cylinderTris[i] + cylinderOffset);
                // the last index needs to be the first index
                // need to stitch the cylinder to the back equator verts
            }



            for (int i = 0; i < backTris.Count; i++)
            {

                // if index is in double verts, use the smaller index.
                int index = backTris[i] + backOffset;

                //   foreach (var e in doubleVerts)
                // {
                //     if (index == e.Item2)
                //     {
                //         index = e.Item1;
                //     }
                // }





                triangles.Add(index);
            }

            // use equator indices to get the corresponding indices from front and back verts
            Dictionary<int, int> map = new Dictionary<int, int>();

            for (int i = 0; i < equatorIndices.Count; i++)
            {
                Vector3 og = sphereVerts[equatorIndices[i]];
                for (int j = 0; j < frontVerts.Count; j++)
                {
                    if (frontVerts[j] == og)
                    {
                        map.Add(i, j);
                    }
                }
            }

            for (int i = frontVerts.Count; i < frontVerts.Count + vertsPerRing; i++)
            {

                //int eq = equatorIndices[i - frontVerts.Count];
                Debug.DrawLine(vertices[i], vertices[i] + Vector3.up * 0.1f, Color.red, 1000f);

                if (map.ContainsKey(i - frontVerts.Count))
                {

                    Vector3 p = vertices[map[i - frontVerts.Count]];
                    Debug.DrawLine(p, p + Vector3.up * 0.1f, Color.blue, 1000f);


                }




            }






            for (int i = 0; i < vertsPerRing - 1; i++)
            {
                int j = i;

                // Vector3 a1 = vertices[map[j]];
                // Vector3 b1 = vertices[map[j + 1]];
                // Vector3 c1 = vertices[frontVerts.Count + j];
                // Vector3 d1 = vertices[frontVerts.Count + j + 1];

                // Debug.DrawLine(a1, b1, Color.green, 1000f);
                // Debug.DrawLine(b1, c1, Color.green, 1000f);
                // Debug.DrawLine(c1, a1, Color.green, 1000f);

                // Debug.DrawLine(c1, b1, Color.green, 1000f);
                // Debug.DrawLine(b1, d1, Color.green, 1000f);
                // Debug.DrawLine(d1, c1, Color.green, 1000f);

                triangles.Add(map[j]);

                triangles.Add(frontVerts.Count + j);
                triangles.Add(map[j + 1]);

                triangles.Add(frontVerts.Count + j);
                triangles.Add(frontVerts.Count + j + 1);
                triangles.Add(map[j + 1]);

                //Top part

              

                
                }

                int last = vertsPerRing - 1;
                triangles.Add(map[last]);
                triangles.Add(frontVerts.Count + last);
                triangles.Add(map[0]);  // Connect back to first vertex

                triangles.Add(frontVerts.Count + last);
                triangles.Add(frontVerts.Count);  // Connect back to first vertex in the other ring
                triangles.Add(map[0]);

           

            // doesnt work for the backvertices, because they are offset 
            
            for (int i = 0; i < equatorIndices.Count; i++)
            {
                Vector3 og = sphereVerts[equatorIndices[i]];
                for (int j = 0; j < backVerts.Count; j++)
                {
                    if (backVerts[j] == og + new Vector3(0, -length, 0)) 
                    {
                        map.Add(i, j + backOffset);
                    }
                }
            }

            for (int i = backOffset - vertsPerRing; i < backOffset; i++)
            {
                Debug.DrawLine(vertices[i], vertices[i] + Vector3.up * 0.1f, Color.yellow, 1000f);

                // if (map.ContainsKey(i - cylinderVerts.Count))
                // {
                //     Debug.Log("GOTCHA!");
                //     Vector3 p = vertices[map[i - backOffset]];
                //     Debug.DrawLine(p, p + Vector3.up * 0.1f, Color.blue, 1000f);
                // }
                 
                 
            }
            
            

            foreach (var k in map.Keys)
            {
                Debug.DrawLine(vertices[map[k]], vertices[map[k]] + Vector3.up * 0.1f, Color.magenta, 1000f);

            }


            // stitch the back hemisphere to the cylinder

            for (int i = backOffset - vertsPerRing; i < backOffset; i++)
            {
                //should be final ring of the cylinder.
                Debug.DrawLine(vertices[i], vertices[i] + Vector3.up * 0.1f, Color.magenta, 1000f);
            }

            // int count = 0;
            // int f = -1;
            // foreach (var e in map)
            // {

            //     Debug.Log($"key: {e.Key} value: {e.Value} map count {map.Count}");
            //     Debug.DrawLine(vertices[e.Value + backOffset], vertices[e.Value + backOffset] + Vector3.up * 0.01f, Color.yellow, 1000f);
            //     //Debug.DrawLine(vertices[backOffset - vertsPerRing], vertices[e.Value + backOffset], Color.green, 1000f);

            //     if (count < (vertsPerRing / 2) - 1)
            //     {
                    
            //         Debug.DrawLine(vertices[backOffset - vertsPerRing / 2 + count], vertices[e.Value + backOffset], Color.red, 1000f);
            //         Debug.DrawLine(vertices[backOffset - vertsPerRing / 2 + count], vertices[backOffset - vertsPerRing / 2 + count + 1], Color.red, 1000f);
                    
            //     }

        


            //     else
            //     {
            //         Debug.DrawLine(vertices[backOffset - vertsPerRing + f], vertices[e.Value + backOffset], Color.blue, 1000f);
            //         Debug.DrawLine(vertices[backOffset - vertsPerRing + f], vertices[backOffset - vertsPerRing + f + 1], Color.blue, 1000f);
            //         f++;
            //     }



            //     // Debug.DrawLine(vertices[backOffset + vertsPerRing + vertsPerRing / 2 + count], vertices[backOffset + vertsPerRing + vertsPerRing / 2 + count +1], Color.red, 1000f);
            //     // Debug.DrawLine(vertices[backOffset - vertsPerRing + vertsPerRing / 2 + count + 1], vertices[e.Value + backOffset], Color.red, 1000f);

            //     count++;

            // }
            
            int count = 0;
int f = 0;
for (int i = 0; i < vertsPerRing - 1; i++)
{
    if (count < (vertsPerRing / 2) - 1)
    {
        // First half - connecting with backOffset - vertsPerRing/2
        triangles.Add(map[i] + backOffset);
         triangles.Add(map[i + 1] + backOffset);
        triangles.Add(backOffset - vertsPerRing / 2 + count + 1);
       

        triangles.Add(backOffset - vertsPerRing / 2 + count);
        triangles.Add(map[i] + backOffset);
        triangles.Add(backOffset - vertsPerRing / 2 + count + 1);
        
    }
    else
    {
        // Second half - connecting with backOffset - vertsPerRing + f
        triangles.Add(map[i] + backOffset);
        triangles.Add(map[i + 1] + backOffset);
        triangles.Add(backOffset - vertsPerRing + f);
        

        triangles.Add(backOffset - vertsPerRing + f);
        triangles.Add(map[i + 1] + backOffset);
        triangles.Add(backOffset - vertsPerRing + f + 1);
        
        f++;
    }
    count++;
}

            //Add final triangles to close the loop
            //You might need to adjust these indices based on your specific case
            triangles.Add(map[vertsPerRing - 1] + backOffset);
            triangles.Add(map[0] + backOffset);
            triangles.Add(backOffset - vertsPerRing + vertsPerRing/2);
            
            

            triangles.Add(map[(vertsPerRing/2)-1] + backOffset);
            triangles.Add(backOffset - vertsPerRing);
            triangles.Add(backOffset - 1);

           
            



            // Create final mesh
            wormMesh.vertices = vertices.ToArray();
            wormMesh.triangles = triangles.ToArray();
            wormMesh.RecalculateNormals();
            wormMesh.RecalculateBounds();

            return wormMesh;
        }

        public static List<(int, int)> GetDoulbeVertices(List<Vector3> verts)
        {

            List<(int, int)> doubleVerts = new List<(int, int)>();

            for (int i = 0; i < verts.Count; i++)
            {
                for (int j = i + 1; j < verts.Count; j++)
                {
                    if (verts[i] == verts[j])
                    {
                        doubleVerts.Add((i, j));
                    }
                }
            }

            return doubleVerts;

        }

        public static Mesh GenerateWormMesh()
        {
            Mesh mesh = new Mesh { name = "C. Elegans" };


            // 1. ) Create IcoSphere
            var baseSphere = GenerateIcoSphere();


            // 2. ) Split IcoSphere into two parts
            // 3 .) Fill in wedges to generate a clean cut
            // 4. ) Extent by placing quad rings until length is filled.
            // 5. ) Generate back hemisphere
            return mesh;
        }

        public static Dictionary<int, int> MapMeshEndings(int[] eq_idx, List<Vector3> hemisphere, List<int> cylinderTris, ref List<Vector3> cylinderVerts, int vertsPerRing)
        {

            Dictionary<int, int> vertexMap = new Dictionary<int, int>();

            for (int i = 0; i < vertsPerRing; i++)
            {
                int idx = eq_idx[i];
                Debug.Log($"idx: {idx} vertsPerRing: {vertsPerRing} eq_idx: {eq_idx.Length} i: {i} hemisphere: {hemisphere.Count}");
                Vector3 vert = hemisphere[idx];

                if (cylinderVerts.Contains(vert))
                {
                    int cylinderIdx = cylinderVerts.IndexOf(vert);
                    vertexMap[cylinderIdx] = idx;
                    //cylinderVerts.RemoveAt(cylinderIdx);
                }
            }



            return vertexMap;
        }

        public static void ReplaceIndices(ref List<int> tris, Dictionary<int, int> vertexMap)
        {
            for (int i = 0; i < tris.Count; i++)
            {
                if (vertexMap.ContainsKey(tris[i]))
                {
                    tris[i] = vertexMap[tris[i]];
                }
            }
        }

        public static void DebugMesh(List<int> tris, List<Vector3> verts, Color col, int offset = 0)
        {

            for (int i = 0; i < tris.Count; i += 3)
            {
                int a = tris[i] + offset;
                int b = tris[i + 1] + offset;
                int c = tris[i + 2] + offset;

                Vector3 A = verts[a];
                Vector3 B = verts[b];
                Vector3 C = verts[c];


                Debug.DrawLine(A, B, col, 1000f);
                Debug.DrawLine(B, C, col, 1000f);
                Debug.DrawLine(C, A, col, 1000f);
            }
            
        }
        
        
        // public struct HemisphereData
        // {
        //     public List<Vector3> vertices;
        //     public List<int> triangles;
        //     public Dictionary<int, int> vertexMap;

        //     public HemisphereData(List<Vector3> vertices, List<int> triangles, Dictionary<int, int> vertexMap)
        //     {
        //         this.vertices = vertices;
        //         this.triangles = triangles;
        //         this.vertexMap = vertexMap;
        //     }

        // }

        // public static HemisphereData GenerateHemisphere(List<Vector3> sphereVerts, List<int> sphereTriangles, float tolerance, float radius)

        // {
        //     HemisphereData hemisphereData = new HemisphereData();

        //     // Add all vertices that have z >= 0
        //     for (int i = 0; i < sphereVerts.Count; i++)
        //     {
        //         if (sphereVerts[i].y <= tolerance)
        //         {
        //             hemisphereData.vertexMap[i] = hemisphereData.vertices.Count;
        //             hemisphereData.vertices.Add(sphereVerts[i] * radius);
        //         }
        //     }


        //     // Add triangles that use these vertices
        //     for (int i = 0; i < sphereTriangles.Count; i += 3)
        //     {
        //         int v1 = sphereTriangles[i];
        //         int v2 = sphereTriangles[i + 1];
        //         int v3 = sphereTriangles[i + 2];

        //         if (hemisphereData.vertexMap.ContainsKey(v1) &&
        //            hemisphereData.vertexMap.ContainsKey(v2) &&
        //             hemisphereData.vertexMap.ContainsKey(v3))
        //         {
        //             hemisphereData.triangles.Add(hemisphereData.vertexMap[v1]);
        //             hemisphereData.triangles.Add(hemisphereData.vertexMap[v2]);
        //             hemisphereData.triangles.Add(hemisphereData.vertexMap[v3]);
        //         }
        //     }

        //     return hemisphereData;


        // }



    }
}