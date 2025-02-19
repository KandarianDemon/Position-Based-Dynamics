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
        public bool[] pinnedIndices;
        
    [Header("Simulation Parameters")]
        public float timestep = 0.001f;

        public int iterations = 5;
        public bool gravity = false;
        public bool pressure_force = false;

        
        [Range(0.01f, 1.0f)]
        [Tooltip("Determines the mass of the cloth")]
        public float cloth_density = 1.0f;
        [Range(0.0f, 1.0f)]
        public float cloth_stiffness = 0.5f;
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

        [Header("External Forces")]
        Vector3 gravity_vec = new Vector3(0.0f,-9.81f,0.0f);
    

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
        
        [Range(0.1f, 10.0f)]
        public float normal_length = 0.5f;


        void Start()
        {
            InitializeSoftbody();
          

            //Debug.Log($"Number of pinned indices: {pinnedIndices.Length}");



            // initialize weights. based on density parameter [kg/m2]




            // proper mass initialization







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
            pinnedIndices[positions.Length-1] = true;

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
        private void OnDrawGizmos() {
            if(positions == null)
            {
                return;
            }

            if (display_vertex_numbers)
            {
                for (int i = 0; i < positions.Length; i++)
                {
                    Handles.Label(positions[i]+transform.position, i.ToString());
                }
            }

            if(triangle_normals)
            {
                if (_triangles == null)
                {
                    return;
                }

                for (int i = 0; i < _triangles.Length; i++)
                {
                    Gizmos.color = Color.blue;
                    Gizmos.DrawLine(_triangles[i].center+transform.position, _triangles[i].center+transform.position + _triangles[i].normal * normal_length);
                }

                for(int j = 0; j < positions.Length;j++){
                    Gizmos.color = Color.red;
                    Vector3 origin = positions[j] + transform.position;
                    Gizmos.DrawLine(origin,origin+velocities[j]*normal_length);
                }
            }
        
            if(d_constraints == null)
            {
                return;
            }

            for(int i = 0; i < d_constraints.Length;i++)
            {
                if(display_vertex_numbers)
                {

                if(d_constraints[i] is DistanceConstraint){
                int a = d_constraints[i].ReturnIndices()[0];
                int b = d_constraints[i].ReturnIndices()[1];

                Gizmos.color = Color.red;
                Gizmos.DrawLine(positions[a]+transform.position,positions[b]+transform.position);
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

                            float angle = Mathf.Acos(Vector3.Dot(n1, n2*-1));

                            Vector3 normal = Vector3.Cross(AB, AC).normalized;

                            Vector3 center1 = (A + B + C) / 3;
                            Vector3 center2 = (A + D + C) / 3;


                        if (similarity > 0) return; 
                            Gizmos.color = Color.blue;
                            Gizmos.DrawLine(center1+transform.position, center1+transform.position + n1*normal_length);

                            Gizmos.DrawSphere(transform.position + center1, 0.1f);



                            Gizmos.color = Color.green;
                            Gizmos.DrawLine(center2+transform.position, center2+transform.position + n2*normal_length);

                            Gizmos.DrawSphere(transform.position + center2, 0.2f);

                            


                        }
                    }
            }
        }



      
        void OnValidate()
        {

            if(d_constraints == null)
            {
                return;
            }
            foreach(var c in d_constraints){
                if(c is DistanceConstraint){
                    DistanceConstraint dc = (DistanceConstraint)c;
                    dc.stiffness = cloth_stiffness;
                }
            }
            
        }

        

        void LateUpdate()
        {

            //Debug.Log($"Number Of Constraints: {d_constraints.Length}");
            // bending constraint!

            positions = mesh.vertices;
            // Apply external forces
          

            

           
            RecalculateTriangles();
            if(collisionConstraints.Count > 0)
            {
                collisionConstraints.Clear();
            }

            //GenerateCollisionConstraints(positions,ref collisionConstraints); 
            //Debug.Log($"Number of Collision Constraints: {collisionConstraints.Count}");

            //Debug.Log($"Number of Collision Constraints: {collisionConstraints.Count}");



            //IntegrateVelocities(positions);
            // This function kills the performance.





            // now i need to loop over all the vertices. Doesnt seem to be efficient really...



            ProjectPositions(positions);
            //ProjectCollisionConstraints(); //just for the test. later wrap both functions in an iterated for loop. Think about substepping.

            ProjectConstraints(); 
            
            //Super slow code. But first get the constraints working and closed meshes and cloth balloons. Then worry about optimization and performance
            //RemoveCollisionConstraints(ref d_constraints);

            mesh.vertices = projectedPositions;
            UpdateVelocities(ref velocities);


            //transform.position = CenterOfMass(positions);


            

        }


        void IntegrateVelocities(Vector3[] positions)
        {
            for (int i = 0; i < positions.Length; i++) {
                //Debug.Log($" number of positions {positions.Length} current index {i} pinned indices {pinnedIndices.Length}");
                if (pinnedIndices[i]) {
                    continue;
                }

                Vector3 pos = positions[i];

                //Replace external forces by array of Force Vectors.


                //need to transform force direction from world to local space. Otherwise its just local space.
                Vector3 netForce = new Vector3(0, 0, 0);

                for (int j = 0; j < forces.Count; j++)
                {
                    netForce += forces[j]*weights[j];
                }
                Vector3 localWind = this.transform.InverseTransformPoint(netForce);

                if(gravity) velocities[i] += weights[i]*this.transform.TransformPoint(new Vector3(0,-0.981f,0)) * timestep;

                //velocities[i] += localWind*wind_strength*timestep;
                
                if (pressure_force) velocities[i] += weights[i] * mesh.normals[i] * 2000.0f * timestep;
                

                
                //*(Mathf.PerlinNoise(pos.x,pos.z)+Time.deltaTime)*
                DampVelocities(i,dampingMethod);  // I made this a function, so that i can later on add different methods of damping and it becomes easier to choose here which
                                    // one should be applied.
                
            }
        }

        void UpdateVelocities(ref Vector3[] velocities){


            for (int i = 0; i < velocities.Length; i++)
            {
                if(pinnedIndices[i]){
                    velocities[i] = Vector3.zero;
                    continue;
                }
                velocities[i] = (projectedPositions[i] - positions[i]) / timestep;
                if (gravity) velocities[i] += weights[i] * gravity_vec * 10.0f * timestep;
                
            }
        }

        void DampVelocities(int i, DampingMethod damp)
        {
            // Damp Velocities.

            switch(damp)
            {
                case DampingMethod.LinearDaming:
                    velocities[i] *= (1-velocity_damping);
                    break;
                case DampingMethod.MuellerDamping:
                    // Mueller Damping
                    Vector3 cm = CenterOfMass(positions);
                    Vector3 vm = CenterOfMass(velocities);

                
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
                
                velocities[i] += timestep * gravity_vec*100.0f*weights[i];
                projectedPositions[i] = positions[i] + velocities[i]  * timestep;
            }
        }

        void ProjectConstraints()
        {
             for(int i = 0; i < iterations; i++){
                
                for(int j = 0; j<d_constraints.Length; j++)
                
                {

                    d_constraints[j].Project(ref projectedPositions,pinnedIndices,iterations,weights);
                }
                
                if(collisionConstraints == null || collisionConstraints.Count == 0)
                {
                    return;
                }
                for(int j = 0; j<collisionConstraints.Count; j++)
                {
                    collisionConstraints[j].Project(ref projectedPositions,pinnedIndices,iterations,weights);
                }
                
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

            return constraintList.ToArray();

        }

        

        void GenerateCollisionConstraints(Vector3[] pos,ref List<Constraint> c)
        {
            if (hash == null) return;

            int totalChecks = 0;
            int constraintsGenerated = 0;

            hash.Create(pos); //First create the hash



            // get a vertex -> triangle map for each collision

            // loop through all the vertices and query for collisions
            List<int> triangles = new List<int>();

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
                            c.Add(new SelfCollisionConstraint(i, tri.a, tri.b, tri.c, cloth_thickness));
                            constraintsGenerated++;

                        }

                    }

                }
                
                Debug.Log($"Total Checks: {totalChecks} Constraints Generated: {constraintsGenerated}");


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

