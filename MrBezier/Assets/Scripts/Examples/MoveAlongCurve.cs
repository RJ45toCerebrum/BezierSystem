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

        // Positions
        C.position = BezierSystem.Evaluate(t, bc);

        // Rotation
        Vector3 tan = BezierSystem.TangentVector(t, bc);
        Quaternion rotator = RotateTo(C.forward, tan);
        C.rotation = rotator * C.rotation;
        // now get the next param...
        t = BezierSystem.UniformSpeedParam(t, L * Time.deltaTime, bc);
        int CC = BezierSystem.CurveCount(bc);
        t = Mathf.Clamp(t, 0, CC);

        // if true pick a random path and go to it else disable it if no next curve available
        if(Mathf.Abs(t - CC) < 0.001f)
        {
            int count = bc.connectedCurves.Length;
            if (count > 0) {
                int pathIndex = Random.Range(0, count);
                bc = bc.connectedCurves[pathIndex];
                t = 0;
            }
            else
                enabled = false;
        }
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
