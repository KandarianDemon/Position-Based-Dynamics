using System.Collections;
using System.Collections.Generic;
using UnityEngine;



public class WormGenerator : MonoBehaviour
{
    // Start is called before the first frame update

    // Generates a Worm and saves it as a prefab
    void Start()
    {
        Mesh  mesh = GenerateWorm();
        this.GetComponent<MeshFilter>().mesh = mesh;
    }

    // Update is called once per frame
    void Update()
    {
        
    }

    Mesh GenerateWorm()
    {

        return new Mesh();
    }

}
