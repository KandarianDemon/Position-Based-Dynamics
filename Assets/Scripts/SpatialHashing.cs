using System;
using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using Unity.VisualScripting;
using UnityEngine;

public class SpatialHashing
{
    // This is supposed to be an unbounded hash grid.
    // We require spacing h and the coordinates
    // coordinates will be translated into cell coordinates (3 component integer vector)
    // these cell coordinates then are hashed, determining a value in a hash array.

    // The has array stores the number of agents currently in the cell represented by the element
    // agents are then sorted by the index of the cell they are in.


    int tableSize;
    public int _TableSize{ get { return tableSize; } set { tableSize = value; } }

    float h;
    public float _H{ get { return h;} set { h = value;}}

    // This class should take particles and triangles and store them in a hash table.

    // containing arrays for the colliding objects.
    int[] hashTable;  // locate index of current particle, increment count by  number of entries
    public int[] _HashTable{ get { return hashTable;} }
    int[] entries; // loop through hashtable, add up all the counts, store in partialSums- cellStart.
    int[] queryIds;  // contains the indices of particles to be evaluated.
    public int[] _QueryIds{ get { return queryIds;}}
    int querySize; // size of the query array
    public int _QuerySize{ get { return querySize;}}



    public SpatialHashing(int tableSize,float h){
        this.tableSize = tableSize;
        this.h = h;

        this.hashTable = new int[tableSize];
        this.entries = new int[tableSize]; // Adjust size for this. it doesnt need to be that large.
        this.queryIds = new int[tableSize]; // Adjust size for this. it doesnt need to be that large.

    }

    public void QueryAll(Vector3[] positions)
    {
        for(int i = 0; i<positions.Length;i++)
        {
            Query(positions[i]);
        }
    }

    public void Create(Vector3[] pos)
    {   

        // Called for every Update Loop

        int numberOfObjects = Mathf.Min(pos.Length,this.hashTable.Length); // get the number of objects to be hashed. In the code by M. Mueller he devides the length of the pos by 3. I think its a 1d array where x,y,z of vertices are stored.
        
        //determine cell size

        ClearHashTable(); // clear the hash table

        // loop over objects,

      
        for(int i = 0; i<numberOfObjects;i++)
        {
            // loops over all the positions and increments the object count for the corresponding hash cell
            // Counts elements in the hash cell
            
            int h = hashPos(pos[i]);
            this.hashTable[h]++;
        }

        int start = 0;

        for(int i = 0; i<this.tableSize;i++)
        {
            //Computes the partial sums and determines the start of the cell.
            start += this.hashTable[i];
            this.hashTable[i] = start;
        }

        this.hashTable[this.tableSize-1] = start; // last cell start is 0


        for(int i = 0; i<numberOfObjects;i++)
        {
            int h = this.hashPos(pos[i]);

            this.hashTable[h]--; // decrement the count of the object in the cell
            this.entries[this.hashTable[h]] = i; // store the index of the object in the cell
        }
        //gets array of positions

        //every update call this will be called and the hash is created from scratch.
        
    }

    public void Query(Vector3 pos) //receives position, index of object and maxDistances (usually 2*spacing)
    {   
        // queryis called for every instance

        // builds an array of indices to be tested for collisions.

        // get hash cell coordinates to define min and max bounds
        Vector3Int minBounds = CellCoordinates(pos-Vector3.one*this.h);
        Vector3Int maxBounds = CellCoordinates(pos+Vector3.one*this.h);

        this.querySize = 0;

        for (int x = minBounds.x; x <= maxBounds.x; x++)
        {
            for (int y = minBounds.y; y <= maxBounds.y; y++)
            {
                for (int z = minBounds.z; z <= maxBounds.z; z++)
                {
                    int h = hashCoords(x, y, z);
                    int start = this.hashTable[h];
                    int end = this.hashTable[h + 1];
                    

                    for (int i = start; i < end; i++)
                    {
                        this.queryIds[this.querySize] = this.entries[i];
                        this.querySize++;
                    }
                }
            }
        }

       

        // loop over neighbor cells
        //this creates a list for indices to be queried.


        
    }

    public Vector3Int CellCoordinates( Vector3 pos){

        return new Vector3Int(Mathf.FloorToInt(pos.x/this.h), Mathf.FloorToInt(pos.y/this.h), Mathf.FloorToInt(pos.z/this.h));
    }

    public int hashPos(Vector3 pos)
    {   
        Vector3Int cellCoords = CellCoordinates(pos);
        return hashCoords(cellCoords.x,cellCoords.y,cellCoords.z);
    }
    public int hashCoords(int x, int y,int z){

        var h = (x*92837111) ^ (y* 689287499) ^ (z* 283923481);	

        return Math.Abs(h % tableSize);
        
    }

    void ClearHashTable()
    {
        for(int i = 0; i < tableSize; i++)
        {
            this.hashTable[i] = 0;
        }
        
    }

   
}
