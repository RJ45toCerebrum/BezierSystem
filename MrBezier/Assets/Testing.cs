using System.Collections;
using System.Collections.Generic;
using UnityEngine;

using EvolveVR.Bezier;

[ExecuteInEditMode]
public class Testing : MonoBehaviour
{
    public Transform C;
    public PlanarBezierComponent bc;
    float t = 0;

    void Update()
    {
        C.position = PlanarBezierSystem.ProjectPoint(C.position, bc);
        C.rotation = bc.planeTransform.rotation;
        Debug.DrawLine(bc.planeTransform.position, C.position, Color.red);

        Vector3[] pofis = PlanarBezierSystem.IntersectCurve(C.position, bc);
        foreach(Vector3 v in pofis)
            Debug.DrawLine(bc.planeTransform.position, v, Color.cyan);
    }
}
