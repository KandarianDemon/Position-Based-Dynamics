using UnityEngine;
using System.Collections;
using System.Collections.Generic;



namespace VirtualWorm.Utils{


public class Map<T>{

    Dictionary<int,List<T>> map = new Dictionary<int,List<T>>();

    public Map(){
        map = new Dictionary<int,List<T>>();

        //edge map (a,b) - c
        // triangle map a - (b,c,d)
    }

    public void Add(int index, T obj){
        if(map.ContainsKey(index)){
            map[index].Add(obj);
        }else{
            map.Add(index,new List<T>());
            map[index].Add(obj);
        }
    }
    
    public int GetCount(int index){
        if(map.ContainsKey(index)){
            return map[index].Count;
        }else{
            return 0;
        }
    }

    public int GetNumberOfKeys(){
        return map.Keys.Count;
    }

    // stores index and a list of objects associated with that index
}


}