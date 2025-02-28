using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEditor;

public class CameraController : MonoBehaviour
{
    // Start is called before the first frame update
    public float speed;
    public float sensitivity;

    private float acc = 1;


    void Update()
    {

        MoveCamera();
      
    }

    void MoveCamera()
    { 

  transform.position += (transform.forward*Input.GetAxis("Vertical")*speed*acc*Time.deltaTime);
        transform.position += (transform.right*Input.GetAxis("Horizontal")*speed*acc*Time.deltaTime);


        if(Input.GetKey(KeyCode.Q))
        {

            transform.position += transform.up * speed*acc*Time.deltaTime;
        }
        

         if(Input.GetKey(KeyCode.E))
        {

            transform.position -= transform.up * speed*acc*Time.deltaTime;
        }
        

        float mouseX = Input.GetAxis("Mouse X");
        float mouseY = Input.GetAxis("Mouse Y");

        if(Input.GetKey(KeyCode.Mouse2))
        {

            transform.eulerAngles += new Vector3(-mouseY*sensitivity, mouseX*sensitivity,0);


        
        }

        if(Input.GetKey(KeyCode.LeftShift))
        {
            acc = 2;
        }

        else {acc = 1;}



        // if(Input.mouseScrollDelta.y != 0)
        // {
        //     transform.position += transform.forward*Input.mouseScrollDelta.y*0.3f;
        // }



        if(Input.GetKey(KeyCode.F) && Selection.activeObject != null )
        {

            Vector3 target = Selection.activeTransform.position;
            Vector3 newDirection = Vector3.RotateTowards(transform.forward, target, 2*Time.deltaTime, 0.0f);

            transform.LookAt(target, Vector3.up);

           


            if((target-transform.position).magnitude >2.0f)
            {
                Vector3 newPos = target - (target-transform.position).normalized * 2.0f;

                transform.position = newPos;
            }

            
            
        
        }
        

    }


}
