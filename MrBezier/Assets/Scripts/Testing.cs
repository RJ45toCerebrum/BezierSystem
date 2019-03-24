using System.Collections;
using System.Collections.Generic;
using UnityEngine;

using EvolveVR.Bezier;

[ExecuteInEditMode]
public class Testing : MonoBehaviour
{
    public PlanarBezierComponent bc;

    Ray ray = new Ray();
    Plane plane = new Plane();

    public Transform other;


    protected void Update()
    {
        ray.origin = other.position;
        ray.direction = other.up;
        plane.SetNormalAndPosition(bc.planeTransform.up, bc.planeTransform.position);

        float t;
        if(plane.Raycast(ray, out t))
        {
            Vector3 ip = ray.origin + t * ray.direction;
            Debug.DrawLine(other.position, ip, Color.blue);
            Debug.DrawLine(bc.planeTransform.position, ip, Color.green);
            Debug.DrawLine(other.position, bc.planeTransform.position, Color.green);

            Vector3 poc = PlanarBezierSystem.IntersectCurveNearestForward(ip, bc);
            Debug.DrawLine(other.position, poc, Color.cyan);

            Vector3 v = ip - bc.planeTransform.position;
            Vector3 w = poc - bc.planeTransform.position;
            if (v.sqrMagnitude > w.sqrMagnitude)
                ContrainedRotation(poc);
            else
                transform.rotation = other.rotation;
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
        Quaternion qs = SwingTo(transform.up, other.up);

        // now try to get the twist...
        Quaternion q = qs * transform.rotation;
        Quaternion qr = Quaternion.Inverse(other.rotation) * q;
        Quaternion qt = Quaternion.AngleAxis(-qr.eulerAngles.y, transform.up);
        transform.rotation = qt * transform.rotation;

        // now swing not entire way but only to the curve
        Vector3 toIp = ip - transform.position;
        toIp.Normalize();
        // =)
        Quaternion swinger = SwingTo(transform.up, toIp);
        transform.rotation = swinger * transform.rotation;
    }
}
