using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using EvolveVR.Bezier;

[ExecuteInEditMode]
public class RotationConstraintTest : MonoBehaviour
{
    public PlanarBezierComponent bc;
    private Ray ray = new Ray();
    private Plane plane = new Plane();


    protected void Update()
    {
        ray.origin = transform.position;
        ray.direction = transform.up;
        plane.SetNormalAndPosition(bc.planeTransform.up, bc.planeTransform.position);

        float t;
        if (plane.Raycast(ray, out t))
        {
            Vector3 ip = ray.origin + t * ray.direction;
            Debug.DrawLine(transform.position, ip, Color.blue);
            Debug.DrawLine(bc.planeTransform.position, ip, Color.green);
            Debug.DrawLine(transform.position, bc.planeTransform.position, Color.green);

            Vector3 poc = PlanarBezierSystem.IntersectCurveNearestForward(ip, bc);
            Debug.DrawLine(transform.position, poc, Color.cyan);

            Vector3 v = ip - bc.planeTransform.position;
            Vector3 w = poc - bc.planeTransform.position;
            if (v.sqrMagnitude > w.sqrMagnitude)
                ContrainedRotation(poc);
        }
    }


    private Quaternion SwingTo(Vector3 v, Vector3 w)
    {
        Vector3 cross = Vector3.Cross(v, w);
        float theta = Vector3.SignedAngle(v, w, cross);
        return Quaternion.AngleAxis(theta, cross);
    }

    private void ContrainedRotation(Vector3 ip)
    {
        // now swing not entire way but only to the curve
        Vector3 toIp = ip - transform.position;
        toIp.Normalize();
        // =)
        Quaternion swinger = SwingTo(transform.up, toIp);
        transform.rotation = swinger * transform.rotation;
    }
}
