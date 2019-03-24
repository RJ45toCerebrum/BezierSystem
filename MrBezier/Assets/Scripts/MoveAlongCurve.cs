using System.Collections;
using System.Collections.Generic;
using UnityEngine;

using EvolveVR.Bezier;

// Shows what you can do with the system...
public class MoveAlongCurve : MonoBehaviour
{
    public Transform C;
    public BezierComponent bc;

    public float L = 0.1f;
    private float t = 0;
    // This script uses the BezierSystem to move a cube along the 
    // curve at a uniform speed with rotation aligned with tangent....
    private void Update()
    {
        if (!bc || !C)
            return;

        C.position = BezierSystem.Evaluate(t, bc);
        Vector3 tan = BezierSystem.TangentVector(t, bc);
        Quaternion rotator = RotateTo(C.forward, tan);
        C.rotation = rotator * C.rotation;
        // now get the next param...
        t = BezierSystem.UniformSpeedParam(t, L * Time.deltaTime, bc);
        t = Mathf.Clamp(t, 0, BezierSystem.CurveCount(bc) - 0.002f);
    }

    private Quaternion RotateTo(Vector3 v, Vector3 w)
    {
        v.Normalize();
        w.Normalize();
        Vector3 axis = Vector3.Cross(v, w);
        float theta = Vector3.SignedAngle(v, w, axis);
        return Quaternion.AngleAxis(theta, axis);
    }

}
