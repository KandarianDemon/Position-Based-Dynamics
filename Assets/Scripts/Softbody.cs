using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Felix.LinAlgSolvers;
using UnityEngine.Animations;
using System.Linq;
using System;
using MathNet.Numerics.Statistics.Mcmc;
using System.Runtime.CompilerServices;

using UnityEditor;
using Unity.VisualScripting;
using UnityEditor.PackageManager;
using Felix.Constraints;
using VirtualWorm.Utils;
using UnityEngine.Profiling;
using UnityEngine.Analytics;
using UnityEditor.EditorTools;







namespace VirtualWorm
{
    public enum ConstraintType
    {
        DistanceConstraint,
        BendingConstraint,
        SelfCollisionConstraint, VolumeConstraint
    }

    // store how many positions are affected by the constraints

    public enum SoftbodyType
    {
        Cloth,
        Worm
    }

    public enum DampingMethod
    {
        LinearDaming,
        MuellerDamping
    }

    public class Triangle
    {
        public int a;
        public int b;
        public int c;

        public Vector3 normal;
        public Vector3 center;

        public Triangle(int a, int b, int c,Vector3[] verts)
        {
            // Upon initialization keep winding order in mind!


            this.a = a;
            this.b = b;
            this.c = c;

            ComputeNormal(ref verts);
            GetCenter(ref verts);
        }

        public void ComputeNormal(ref Vector3[] verts){

            Vector3 A = verts[a];
            Vector3 B = verts[b];
            Vector3 C = verts[c];

            Vector3 AB = B - A;
            Vector3 AC = C - A;

            this.normal = Vector3.Cross(AB, AC).normalized;
        }

        public void GetCenter(ref Vector3[] verts)
        {
            Vector3 A = verts[a];
            Vector3 B = verts[b];
            Vector3 C = verts[c];

            this.center = (A + B + C) / 3;
        }

        public float ComputeScalarTripleProduct(ref Vector3[] verts){

            Vector3 cross_AB = Vector3.Cross(verts[a],verts[b]);
            
            return Vector3.Dot(cross_AB,verts[c]);



        }
    }
   

    




    [RequireComponent(typeof(MeshFilter), typeof(MeshRenderer))]
    public class Softbody : MonoBehaviour
    {
        Mesh mesh;
        Vector3[] velocities;
        Vector3[] positions;
        Vector3[] projectedPositions;
        float[] weights;
        SpatialHashing hash;
        Constraint[] d_constraints;
        List<Constraint> collisionConstraints;
        StaticCollisionConstraint[] staticCollisionConstraints;
        public bool[] pinnedIndices;
        Collider[] collisionObjects;

        [Header("Simulation Parameters")]
     public SoftbodyType softbodyType = SoftbodyType.Cloth;
        public float timestep = 0.001f;

        public int iterations = 5;
        public bool gravity = false;
        public bool pressure_force = false;
        public bool trap_box = false;


        [Range(0.01f, 1.0f)]
        [Tooltip("Determines the mass of the cloth")]
        [Header("Cloth Settings")]
        public float cloth_density = 1.0f;
        [Range(0.0f, 1.0f)]
        public float cloth_stiffness = 0.5f;
        [Range(0.0f, 1.0f)]
        public float bending_stiffness = 0.5f;

        [Range(0.0f,2.0f)]
        public float volume_stiffness = 0.5f;
        [Range(0.0f, 10.0f)]
        public float pressure = 0.1f;
        [Tooltip("Determines the thickness of the cloth. Used for self collision detection")]
        public float cloth_thickness = 0.5f;
        float mass;

        [Tooltip("Velocity damping method. Linear Damping is a simple damping method. Mueller Damping is a more sophisticated method, as described in the PBD paper.")]
        public DampingMethod dampingMethod;

        [Range(0.0f,1.0f)]
        public float velocity_damping = 0.1f;

        Dictionary<(int, int),List<int>> edgeMap; // Contains information which point is part of which edge
        Dictionary<int,List<int>> triangleMap; // contains information which point is part of which triangle
        Triangle[] _triangles;

        Map<int> map = new Utils.Map<int>();

        
        Vector3 gravity_vec = new Vector3(0.0f,-9.81f,0.0f);
        [Header("External Forces")]
        public float dragCoefficient = 0.1f;
        public float bounciness = 0.5f;
        public float friction = 0.5f;
    

        public List<Vector3> forces = new List<Vector3>();
        public float wind_strength = 0.1f;
        [Header("Debugging")]

        [Tooltip("Displays the index of the vertex on the mesh in the editor view")]
        public bool display_vertex_numbers = false;
        [Tooltip("Displays the normals of the triangles")]
        public bool display_normals = false;
        [Tooltip("Pin the frame of the cloth. Dont use it for anything that is not a plane.")]
        public bool pin_frame;
        [Tooltip("Displays the normals of the triangles and the velocity vectors of the vertices")]
        public bool triangle_normals;
        [Tooltip("Displays the tangential and orthogonal forces")]
        public bool display_drag = false;
        
        [Range(0.1f, 10.0f)]
        public float normal_length = 0.5f;


        void Start()
        {
            InitializeSoftbody();
            collisionObjects = FindObjectsOfType<Collider>();
          

            //Debug.Log($"Number of pinned indices: {pinnedIndices.Length}");



            // initialize weights. based on density parameter [kg/m2]




            // proper mass initialization







        }
        void UpdateTransformPosition()
        {
            Vector3 center = GeometryCenter(projectedPositions);
            //transform.position = center;
        }
        public void InitializeSoftbody()
        {
            this.mesh = GetComponent<MeshFilter>().mesh;
            positions = mesh.vertices;
            int[] triangles = mesh.triangles;
            //RemoveDuplicateVertices(ref positions, ref triangles);
            projectedPositions = positions;

            velocities = new Vector3[positions.Length];

            edgeMap = GetEdgeMap(triangles);
            triangleMap = GetTriangleMap(triangles);

            //Debug.Log($"Number of Triangles: {triangles.Length / 3} mesh triangles {mesh.triangles.Length / 3}");
            weights = InitializeWeights();
            _triangles = StoreTriangles(triangles);
            float radius = (positions[triangles[0]] - positions[triangles[1]]).magnitude;
            hash = new SpatialHashing(20000, cloth_thickness);
            //Debug.Log($"radius: {radius} number of verts: {positions.Length}");
            HashSet<Vector3> uniqueVertices = new HashSet<Vector3>(mesh.vertices);
            //Debug.Log($"Unique Vertices: {uniqueVertices.Count}, Total Vertices: {mesh.vertices.Length}");



            // foreach(var w in weights)
            // {

            //  Debug.Log($"Weight: {1.0f/w}");
            // }



            d_constraints = InitializeDistanceConstraints(mesh.triangles);

            collisionConstraints = new List<Constraint>();

            pinnedIndices = new bool[positions.Length];
            //pinnedIndices[positions.Length - 1] = true;
          

            // if (positions.Length >= 12)
            // {


                //     for (int i = positions.Length - 1; i > 0; i--)
                //     {
                //         if (i > positions.Length - 12)
                //         {
                //             pinnedIndices[i] = true;

                //         }

                //         else if (i % 11 == 0 && pin_frame)
                //         {
                //             pinnedIndices[i] = true;
                //         }

                //         else if (pin_frame && new int[] { 0, 1, 21, 32, 43, 54, 65, 76, 87, 98, 109 }.Contains(i))
                //         {
                //             pinnedIndices[i] = true;
                //         }

                //         else if (i < 11 && pin_frame)
                //         {
                //             pinnedIndices[i] = true;
                //         }
                //         else
                //         {
                //             pinnedIndices[i] = false;
                //         }


                //     }
                // }

                // else
                // {
                //     pinnedIndices[0] = true;
                // }
        }
        private void OnDrawGizmos()
        {
            if (positions == null)
            {
                return;
            }

            if (display_vertex_numbers)
            {
                for (int i = 0; i < positions.Length; i++)
                {
                    Handles.Label(positions[i] + transform.position, i.ToString());
                }
            }

            if (triangle_normals)
            {
                if (_triangles == null)
                {
                    return;
                }

                for (int i = 0; i < _triangles.Length; i++)
                {
                    Gizmos.color = Color.blue;
                    Gizmos.DrawLine(_triangles[i].center + transform.position, _triangles[i].center + transform.position + _triangles[i].normal * normal_length);
                }

                for (int j = 0; j < positions.Length; j++)
                {
                    Gizmos.color = Color.red;
                    Vector3 origin = positions[j] + transform.position;
                    Gizmos.DrawLine(origin, origin + velocities[j] * normal_length);
                }
            }

            if (d_constraints == null)
            {
                return;
            }

            for (int i = 0; i < d_constraints.Length; i++)
            {
                if (display_vertex_numbers)
                {

                    if (d_constraints[i] is DistanceConstraint)
                    {
                        int a = d_constraints[i].ReturnIndices()[0];
                        int b = d_constraints[i].ReturnIndices()[1];

                        Gizmos.color = Color.red;
                        Gizmos.DrawLine(positions[a] + transform.position, positions[b] + transform.position);
                    }




                }

                if (display_normals)
                {
                    if (d_constraints[i] is BendingConstraint)
                    {

                        int a = d_constraints[i].ReturnIndices()[0];
                        int b = d_constraints[i].ReturnIndices()[1];
                        int c = d_constraints[i].ReturnIndices()[2];
                        int d = d_constraints[i].ReturnIndices()[3];

                        Vector3 A = positions[a];
                        Vector3 B = positions[b];
                        Vector3 C = positions[c];
                        Vector3 D = positions[d];

                        Vector3 AB = positions[b] - positions[a];
                        Vector3 AC = positions[c] - positions[a];
                        Vector3 AD = positions[a] - positions[d];

                        Vector3 BC = positions[c] - positions[b];
                        Vector3 BD = positions[d] - positions[b];

                        Vector3 n1 = Vector3.Cross(AB, AC).normalized;
                        Vector3 n2 = Vector3.Cross(AC, AD).normalized;

                        float similarity = Vector3.Dot(n1, n2);

                        float angle = Mathf.Acos(Vector3.Dot(n1, n2 * -1));

                        Vector3 normal = Vector3.Cross(AB, AC).normalized;

                        Vector3 center1 = (A + B + C) / 3;
                        Vector3 center2 = (A + D + C) / 3;


                        if (similarity > 0) return;
                        Gizmos.color = Color.blue;
                        Gizmos.DrawLine(center1 + transform.position, center1 + transform.position + n1 * normal_length);

                        Gizmos.DrawSphere(transform.position + center1, 0.1f);



                        Gizmos.color = Color.green;
                        Gizmos.DrawLine(center2 + transform.position, center2 + transform.position + n2 * normal_length);

                        Gizmos.DrawSphere(transform.position + center2, 0.2f);




                    }
                }
            }

            if (display_drag)
            {
                for (int i = 0; i < positions.Length; i++)
                {
                    Vector3 normal = Vector3.zero;
                    Vector3 tangent = Vector3.zero;

                    Triangle[] tris = triangleMap[i].Select(x => _triangles[x]).ToArray();
                    foreach (var t in tris)
                    {
                        normal += t.normal;
                    }

                    tangent = Vector3.Cross(velocities[i], normal).normalized;

                    Gizmos.color = Color.yellow;
                    Gizmos.DrawLine(positions[i] + transform.position, positions[i] + transform.position + normal.normalized * normal_length);

                    Gizmos.color = Color.red;
                    Gizmos.DrawLine(positions[i] + transform.position, positions[i] + transform.position + tangent.normalized * normal_length);

                    Gizmos.color = Color.green;
                    Gizmos.DrawLine(positions[i] + transform.position, positions[i] + transform.position + velocities[i].normalized*normal_length);
                }
            }
        }



      
        void OnValidate()
        {

            if(d_constraints == null)
            {
                return;
            }
            foreach (var c in d_constraints)
            {
                if (c is DistanceConstraint)
                {
                    DistanceConstraint dc = (DistanceConstraint)c;
                    dc.stiffness = cloth_stiffness;
                }
                if (c is VolumeConstraint)
                {
                    VolumeConstraint vc = (VolumeConstraint)c;
                    vc.stiffness = volume_stiffness;
                    vc.pressure = pressure;
                }

                if (c is BendingConstraint)
                {
                    BendingConstraint bc = (BendingConstraint)c;
                    bc.stiffness = bending_stiffness;
                }

              
            }
            
        }



        void Update()
        {

            //Debug.Log($"Number Of Constraints: {d_constraints.Length}");
            // bending constraint!

            //timestep = Time.deltaTime;

            positions = mesh.vertices;
            // Apply external forces
           



            RecalculateTriangles();
            if (collisionConstraints.Count > 0)
            {
                collisionConstraints.Clear();
            }



            IntegrateVelocities(positions);
            DampVelocities();
            ProjectPositions(positions);
            
            StaticCollisions();
            
          
            ProjectConstraints();

            if (trap_box)
            {
                 for (int i = 0; i < projectedPositions.Length; i++)
            {
                if (projectedPositions[i].y < 0)
                {
                    projectedPositions[i].y = 0.001f;
                }

                if (projectedPositions[i].y > 10)
                {
                    projectedPositions[i].y = 9.999f;
                }

                if (projectedPositions[i].x < -10)
                {
                    projectedPositions[i].x = -9.999f;
                }
                if (projectedPositions[i].x > 10)
                {
                    projectedPositions[i].x = 9.999f;
                }
                if (projectedPositions[i].z < -10)
                {
                    projectedPositions[i].z = -9.999f;
                }
                if(projectedPositions[i].z > 10)
                {
                    projectedPositions[i].z = 9.999f;
                }
            }
            }
           
            


            UpdateVelocities(ref velocities);
            AddRestitutionAndFriction(ref velocities, ref staticCollisionConstraints);
            mesh.vertices = projectedPositions;
            
            mesh.RecalculateNormals();
            UpdateTransformPosition();


            




        }
        public Vector3 GetPressureDirection(int index)
        {

            Triangle[] tris = triangleMap[index].Select(x => _triangles[x]).ToArray();

            Vector3 netForce = new Vector3(0, 0, 0);

            for (int i = 0; i < tris.Length; i++)
            {
                netForce += tris[i].normal;
            }

            return netForce/tris.Length;
            

        }


        void IntegrateVelocities(Vector3[] positions)
        {
            for (int i = 0; i < positions.Length; i++)
            {
                //Debug.Log($" number of positions {positions.Length} current index {i} pinned indices {pinnedIndices.Length}");
                if (pinnedIndices[i])
                {
                    continue;
                }

                if (gravity) velocities[i] += timestep * gravity_vec * (1 / weights[i]);
                //DampVelocities(i, dampingMethod);

                // Triangle[] tris = triangleMap[i].Select(x => _triangles[x]).ToArray();

                // Vector3 normal = Vector3.zero;

                // foreach (var t in tris)
                // {
                //     normal += t.normal;
                // }

                // normal.Normalize();
                // Vector3 normal2 = Vector3.Cross(velocities[i], normal);

                // Vector3 w1 = Vector3.Dot(velocities[i], normal) * normal;
                // Vector3 w2 = velocities[i] - w1;

            

            }
        }

        void UpdateVelocities(ref Vector3[] velocities)
        {
            //  



            for (int i = 0; i < velocities.Length; i++)
            {
                if (pinnedIndices[i])
                {
                    velocities[i] = Vector3.zero;
                    continue;
                }
                velocities[i] = (projectedPositions[i] - positions[i]) / timestep;






            }
        }

        void AddRestitutionAndFriction(ref Vector3[] velocities, ref StaticCollisionConstraint[] constraints)
            {

                for (int i = 0; i < constraints.Length; i++)
                {
                    
                    if (constraints[i] == null)
                    {
                        continue;
                    }

                    Vector3 normal = constraints[i].collisionNormal;
                    Vector3 point = constraints[i].collisionPoint;
                    Vector3 velocity = velocities[constraints[i].index];
                    

                float w = 1 / weights[i];

                    // Calculate normal and tangent components of velocity
            Vector3 normal_velocity = Vector3.Dot(velocity, normal) * normal;
            Vector3 tangent_velocity = velocity - normal_velocity;  

            // Apply restitution (bounce) in the normal direction
            Vector3 normal_impulse = -bounciness * normal_velocity;

            // Apply friction in the tangent direction (scaled by the inverse mass)
            Vector3 tangent_impulse = -friction * tangent_velocity;

            // Correct velocities using inverse mass and time step to ensure consistency across different frame rates
            //velocities[constraints[i].index] += (normal_impulse + tangent_impulse) * (1 / weights[constraints[i].index]) * timestep;



                }
            


            }
        
        
        (Vector3, Vector3) GetCenterOfMassAndMeanVelocity(Vector3[] positions, Vector3[] velocities)
        {
            Vector3 centerOfMass = CenterOfMass(positions);
            Vector3 meanVelocity = CenterOfMass(velocities);

           
            return (centerOfMass,meanVelocity);
        }

        void DampVelocities()
        {
            (Vector3, Vector3) cmvm = GetCenterOfMassAndMeanVelocity(positions, velocities);
            
            for(int i = 0; i<positions.Length; i++)
            {
             DampVelocities(i,dampingMethod,cmvm.Item1,cmvm.Item2);
            }
        }

        void DampVelocities(int i, DampingMethod damp, Vector3 cm, Vector3 vm)
        {
            // Damp Velocities.

            switch(damp)
            {
                case DampingMethod.LinearDaming:
                    velocities[i] *= (1-velocity_damping);
                    break;
                case DampingMethod.MuellerDamping:
                    // Mueller Damping
                    

                
                    Vector3 angular_Momentum = ComuteAngularMomentum(cm);

                // Inertia Tensor
                    Matrix4x4 inertiaTensor = new Matrix4x4();
                    Vector3[] distances = Vector3.zero.DistancesToCenter(cm,positions);
                
                
                    inertiaTensor.GetInertiaTensor(weights,distances);
                    Matrix4x4 inverseInertiaTensor = inertiaTensor.inverse;

                
                

                //Careful this is a Matrix4x4,so we only want to look at the upper 3x3 part of it.

                    Vector3 angular_velocity = inverseInertiaTensor.MultiplyVector(angular_Momentum);
                    

                    Vector3 delta_v = vm + Vector3.Cross(angular_velocity, distances[i]) - velocities[i];
                    velocities[i] += velocity_damping*delta_v;
                    break;
            }
           
           
        }
        

        void ProjectPositions(Vector3[] positions)
        {
             for(int i = 0; i<positions.Length; i++){

                if(pinnedIndices[i]){
                    continue;
                }
                
               



                projectedPositions[i] = positions[i] + velocities[i]  * timestep;
            }
        }

        void ProjectConstraints()
        {
            
           
            for (int i = 0; i < iterations; i++)
            {

                 for (int j = 0; j < staticCollisionConstraints.Length; j++)
                    {
                        if (staticCollisionConstraints[j] == null)
                        {
                        continue;
                        }
                        staticCollisionConstraints[j].Project(ref projectedPositions, pinnedIndices, iterations, weights);
                    }



                for (int j = 0; j < d_constraints.Length; j++)

                {

                    d_constraints[j].Project(ref projectedPositions, pinnedIndices, iterations, weights);
                }




                // if (collisionConstraints == null || collisionConstraints.Count == 0)
                // {
                //     continue;
                // }
                // for (int j = 0; j < collisionConstraints.Count; j++)
                // {
                //     collisionConstraints[j].Project(ref projectedPositions, pinnedIndices, iterations, weights);
                // }


                //Debug.Log($"Number of static collision constraints: {staticCollisionConstraints.Length}");

                // for the distance constraint we need: the normal, the inverse mass and the error

                //error


            }
        }

        void RemoveCollisionConstraints(ref Constraint[] constraints)
        {
            List<Constraint> c = constraints.ToList();
            for(int i = 0; i<c.Count; i++)
            {
                if(c[i] is SelfCollisionConstraint)
                {
                    c.RemoveAt(i);
                }
            }

            constraints = c.ToArray();
        }

        void ProjectCollisionConstraints()
        {
            if(collisionConstraints == null)
            {
                return;
            }

            for(int i = 0; i < iterations; i++)
            {
                for(int j = 0; j<collisionConstraints.Count; j++)
                {
                    collisionConstraints[j].Project(ref projectedPositions,pinnedIndices,iterations,weights);
                }
            }
        }


        Triangle[] StoreTriangles(int[] triangles)
        {

           
            int num = Mathf.FloorToInt(triangles.Length/3);
             if(_triangles == null) {
                _triangles = new Triangle[num];
            }

            Triangle[] tris = new Triangle[num];

            for(int i = 0; i<num; i++)
            {
                _triangles[i] = new Triangle(
                                            triangles[i*3],
                                            triangles[i*3 + 1], // Was incorrectly (i+1)*3
                                            triangles[i*3 + 2], // Was incorrectly (i+2)*3
                                            positions
                                            );
            }

            return _triangles;
        }

        public void RecalculateTriangles(){

            if (_triangles == null || _triangles.Length == 0)
            {
                print(_triangles==null);
                return;
            }
            

            for(int i = 0; i<_triangles.Length; i++)
            {

                _triangles[i].ComputeNormal(ref positions);
                _triangles[i].GetCenter(ref positions);
            }
        }

     
        

        Constraint[] InitializeDistanceConstraints( int[] triangles)
        {
            List<Constraint> constraintList = new List<Constraint>();
            


            // Initialize Distance Constraints
            for(int i = 0; i<triangles.Length; i+=3){

                // i think i might initialize some constraints double

                constraintList.Add(new DistanceConstraint(triangles[i],triangles[i+1],Vector3.Distance(positions[triangles[i]],positions[triangles[i+1]]*1.01f),cloth_stiffness,mass));
                constraintList.Add(new DistanceConstraint(triangles[i+1],triangles[i+2],Vector3.Distance(positions[triangles[i+1]],positions[triangles[i+2]]*1.01f),cloth_stiffness,mass));
                constraintList.Add(new DistanceConstraint(triangles[i+2],triangles[i],Vector3.Distance(positions[triangles[i+2]],positions[triangles[i]]*1.01f),cloth_stiffness,mass));
                

                //Debug.Log($"constraint {i} :: {triangles[i]} || {triangles[i+1]}");
            }

            // Initialize Bending Constraints
            for(int i = 0; i<edgeMap.Keys.Count;i++)
            {
                var key = edgeMap.Keys.ElementAt(i);
                var entry = edgeMap[key];



                if (entry.Count > 1) {

                    int A = key.Item1;
                    int C = key.Item2;

                    // get indices of triangle 1
                    // get indicesof triangle 2

                    // check winding order.


                    int thirdVertex1 = triangles[entry[0] * 3] == A || triangles[entry[0] * 3] == C ?
                                       (triangles[entry[0] * 3 + 1] == A || triangles[entry[0] * 3 + 1] == C ?
                                       triangles[entry[0] * 3 + 2] : triangles[entry[0] * 3 + 1]) : triangles[entry[0] * 3];


                    int thirdVertex2 = triangles[entry[1] * 3] == A || triangles[entry[1] * 3] == C ?
                                       (triangles[entry[1] * 3 + 1] == A || triangles[entry[1] * 3 + 1] == C ?
                                       triangles[entry[1] * 3 + 2] : triangles[entry[1] * 3 + 1]) : triangles[entry[1] * 3];
                    int B = thirdVertex1;
                    int D = thirdVertex2;

                    Vector3 n1 = Vector3.Cross(positions[B] - positions[A], positions[C] - positions[A]).normalized;
                    Vector3 n2 = Vector3.Cross(positions[C] - positions[A], positions[D] - positions[A]).normalized;
                    float alpha_0 = Mathf.Acos(Vector3.Dot(n1, n2 * -1));

                    //constraintList.Add(new DistanceConstraint(A, B, Vector3.Distance(positions[A], positions[B]), cloth_stiffness, mass));
                    constraintList.Add(new BendingConstraint(A,B,C,D,0.7f, alpha_0,mass));

                    // Initialize normals using mesh.normals instead of calculating from the vertices. That should help to ensure all normals point into the same direction.
                }
                // only for edges that are shared by two triangles!
                
                // collision constraints are generated in each update loop. So.

            }

            // compute initial volume

            if(softbodyType == SoftbodyType.Worm)
            {
                AddVolumeConstraint(ref constraintList);
            }
            


            return constraintList.ToArray();

        }

        public void StaticCollisions()
        {
            // Static Collisions
            // Check if the vertices are inside the collider.
            // If they are, project them outwards.
            // If they are not, do nothing.

            // for now just have a scene with few objects. plane, cube, sphere
            if (collisionObjects == null)
            {
                return;
            }


            if (staticCollisionConstraints == null)
            {
                staticCollisionConstraints = new StaticCollisionConstraint[collisionObjects.Length * positions.Length];
            }

            for (int i = 0; i < staticCollisionConstraints.Length; i++)
            {
                staticCollisionConstraints[i] = null;
            }

            int constraintIndex = 0;
            for (int i = 0; i < positions.Length; i++)
            {
                RaycastHit hit;
                Vector3 direction = projectedPositions[i] - positions[i];
                if (Physics.Raycast(this.transform.TransformPoint(positions[i]), direction, out hit, direction.magnitude))
                {
                    Vector3 normal = hit.normal;
                    Vector3 point = hit.point;
                    //Vector3 projected = positions[i] + Vector3.Project(point - positions[i], normal);

                    staticCollisionConstraints[constraintIndex] = new StaticCollisionConstraint(i, point, normal);
                    constraintIndex++;
                    Debug.Log($"Collision Detected at {i} with {hit.collider.name}");
                    // projectedPositions[i] = projected;
                    // velocities[i] = velocities[i] * -1;


                    // for the collision we need the point of the collision and the collision normal.
                    // then we apply an inequality constraint.
                }

            }
            
            //Debug.Log($"Number of Static Collision Constraints: {constraintIndex}");


            
        }


        float ComputeMeshVolume()
        { 
            float currentVolume = 0f;
            foreach (var triangle in _triangles)
            {
                currentVolume += triangle.ComputeScalarTripleProduct(ref positions);
            }
            return currentVolume /= 6f; // Actual volume is sum/6
        }

        void AddVolumeConstraint(ref List<Constraint> c)
        {
            float currentVolume = ComputeMeshVolume();
            c.Add(new VolumeConstraint(currentVolume,_triangles,stiffness:cloth_stiffness));
        }
        void GenerateCollisionConstraints(Vector3[] pos, ref List<Constraint> c)
        {
            if (hash == null) return;

            int totalChecks = 0;
            int constraintsGenerated = 0;

            hash.Create(pos); //First create the hash
            SelfCollisionConstraint[] vc = new SelfCollisionConstraint[30000];
            



            // get a vertex -> triangle map for each collision

            // loop through all the vertices and query for collisions
            List<int> triangles = new List<int>();


            // for self collisions i need to test vertices against triangles
            
            for (int i = 0; i < pos.Length; i++)
            {

                Vector3 queryPosition = pos[i];
                Vector3 projectedQueryPosition = queryPosition + velocities[i] * timestep;
                hash.Query(queryPosition); //returns querySize and queryIds.


                if (hash._QuerySize == 0)
                {
                    continue;
                }

                for (int j = 0; j < hash._QuerySize; j++)
                {
                    int id = hash._QueryIds[j];
                    if (id == i)
                    {
                        continue;
                    }

                    triangles = triangleMap[id];
                    //Debug.Log($"Number of Triangles: {triangles.Count} i {i}");
                    float epsilon = 1e-6f;

                    if (triangles.Count == 0)
                    {
                        continue;
                    }

                    for (int k = 0; k < triangles.Count; k++)
                    {
                        Triangle tri = _triangles[k];

                        Vector3 normal = tri.normal;
                        Vector3 AB = queryPosition - pos[tri.a];

                        float d = Vector3.Dot(normal, AB);

                        totalChecks++;

                        if (d + epsilon < cloth_thickness)
                        {

                            //Generate the collision constraint.
                            //c.Add(new SelfCollisionConstraint(i, tri.a, tri.b, tri.c, cloth_thickness));
                            if(constraintsGenerated >= vc.Length)
                            {
                                Array.Resize(ref vc,vc.Length*2);
                            }
                            vc[constraintsGenerated] = new SelfCollisionConstraint(i, tri.a, tri.b, tri.c, cloth_thickness);
                            constraintsGenerated++;

                        }

                    }

                }

                //Debug.Log($"Total Checks: {totalChecks} Constraints Generated: {constraintsGenerated}");


            }


            // Generate SelfCollisionConstraint(vertex, triangle)
            // Generate Collision Constraints
            // called for each update






            // 

            // Destroy constraints after each update.
        }

       

        float[] InitializeWeights(){

            float[] w = new float[positions.Length];

            for(int i = 0; i < positions.Length;i++)
            {

                // Find adjacent triangles

                if (triangleMap == null)
                { 
                    throw new Exception("Triangle Map is null");
                }
                //Debug.Log($"Triangle Map: {triangleMap.Count}");
                int[] adjacentTriangles = triangleMap[i].ToArray();

                float weight = 0;

                for(int j = 0; j<adjacentTriangles.Length;j++)
                {

                    // make the triangle map a class and give it access methods. better usability.
                    int t = adjacentTriangles[j];
                    int a = mesh.triangles[t * 3];
                    int b = mesh.triangles[t * 3 + 1];
                    int c = mesh.triangles[t * 3 + 2];

                    Vector3 A = positions[a];
                    Vector3 B = positions[b];
                    Vector3 C = positions[c];

                    Vector3 AB = B - A;
                    Vector3 AC = C - A;

                    
                    float area = Vector3.Cross(AB, AC).magnitude / 2;

                    weight += (area * cloth_density);
                }

                w[i] = weight/3.0f;


                

               


            }


            return w;

        }

        
        // Constraint[] InitializeBendingConstraints(int[] triangles)
        // {
        //     List<Constraint> constraintList = new List<Constraint>();

        //     for(int i = 0; i<triangles.Length; i+=3){
        //         DistanceConstraint constraint1 = new BendingConstraint(triangles[i],triangles[i+1],triangles[i+2],triangles[i],triangles[i+1],triangles[i+2],1.0f);
        //     }

        //     return constraintList.ToArray();
        // }

        Vector3 CenterOfMass(Vector3[] verts)
        {
            // This function can also be used for velocities!
            Vector3 center = new Vector3(0,0,0);
            float sum_of_weights = 0;

            for(int i = 0; i<verts.Length; i++){
               center += verts[i]*weights[i];
               sum_of_weights += weights[i];
            }

            return center/sum_of_weights;
        }

        Vector3 GeometryCenter(Vector3[] verts)
        {
            Vector3 center = new Vector3(0,0,0);

            for(int i = 0; i<verts.Length; i++)
            {
                center += verts[i];
            }

            return center/verts.Length;
        }

        Vector3 ComuteAngularMomentum(Vector3 cm)
        {   
            Vector3 L = new Vector3(0,0,0);
            for(int i = 0; i<positions.Length; i++)
            {
                Vector3 r = positions[i] - cm;
                Vector3 m_v = weights[i] * velocities[i];
                L += Vector3.Cross(r,m_v);
            }

            return L;


        }

        void RemoveDuplicateVertices(ref Vector3[] vertices, ref int[] triangles)
        {
            Dictionary<Vector3, int> uniqueVertices = new Dictionary<Vector3, int>();
            List<Vector3> newVertices = new List<Vector3>();
            int[] remap = new int[vertices.Length];

            // Build a remap table for unique positions
            for (int i = 0; i < vertices.Length; i++)
            {
                Vector3 pos = vertices[i];

                if (!uniqueVertices.TryGetValue(pos, out int index))
                {
                    index = newVertices.Count;
                    uniqueVertices[pos] = index;
                    newVertices.Add(pos);
                }

                remap[i] = index; // Store new index for remapping triangles
            }

            // Remap triangle indices
            for (int i = 0; i < triangles.Length; i++)
            {
                triangles[i] = remap[triangles[i]];
            }

            // Replace old vertex array
            vertices = newVertices.ToArray();
        }


        

        Dictionary<(int, int), List<int>> GetEdgeMap(int[] triangles)
        {

            //Some meshes have duplicated vertices, apparently to work around uv seams. Here i need to merge the vertices first.
            Dictionary<(int, int), List<int>> edgeMap = new Dictionary<(int, int), List<int>>();

            int tri_count = 0;
            for (int i = 0; i < triangles.Length; i += 3)
            {
                int a = triangles[i];
                int b = triangles[i + 1];
                int c = triangles[i + 2];

                // Add edge (a, b)
                var edge1 = a < b ? (a, b) : (b, a);
                if (!edgeMap.ContainsKey(edge1))
                {
                    edgeMap[edge1] = new List<int>();
                }
                edgeMap[edge1].Add(tri_count); // Add triangle index

                // Add edge (b, c)
                var edge2 = b < c ? (b, c) : (c, b);
                if (!edgeMap.ContainsKey(edge2))
                {
                    edgeMap[edge2] = new List<int>();
                }
                edgeMap[edge2].Add(tri_count); // Add triangle index

                // Add edge (c, a)
                var edge3 = c < a ? (c, a) : (a, c);
                if (!edgeMap.ContainsKey(edge3))
                {
                    edgeMap[edge3] = new List<int>();
                }
                edgeMap[edge3].Add(tri_count); // Add triangle index

                tri_count++;
            }

            return edgeMap;
        }


        Dictionary<int,List<int>> GetTriangleMap(int[] triangles)
        {
            //Triangle map contains information, which triangles a vertex belongs to.


            Dictionary<int,List<int>> triangleMap = new Dictionary<int, List<int>>();

            // for each index already in the list,
            
            for(int i = 0,count = 0; i<triangles.Length;i+=3,count++)
            {
                // check if the vertices are in the dictionary, if yes...append the count to the list if it is not there already

                for(int j = 0; j<3; j++){
                    if(triangleMap.ContainsKey(triangles[i+j])){
                        if(!triangleMap[triangles[i+j]].Contains(count)){
                            triangleMap[triangles[i+j]].Add(count);
                        }
                    }
                    else{
                        triangleMap[triangles[i+j]] = new List<int>{count};
                    }
                }
            }
            return triangleMap;
        
        }

        

        // Get the Mesh Data
        // Generate Distance constraint for each edge in the mesh
        // Apply external forces
        // Project positions
        // Project Constraints
        // Update Mesh and positions, velocities
    }
}

