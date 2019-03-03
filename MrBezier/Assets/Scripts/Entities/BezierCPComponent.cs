using UnityEngine;

namespace EvolveVR.Bezier
{
    /// <summary>
    /// Data Container for Bezier curve Control Points.
    /// Add this script to gameobjects that need to be
    /// a control point for bezier component. Then add this
    /// to the BezierCPComponent array on your Bezier Component.
    /// </summary>
    public class BezierCPComponent : MonoBehaviour
    {
        [Tooltip("This Tells the Bezier system to make the curve at this control point smooth." +
                 "Make Sure that if this is check then you must have the adjcp0 and adjcp1 references set")]
        public bool tangent = false;
        public float controlPointRadius = 1;
        public Color controlPointColor;
        [Tooltip("This Tells bezier system how the control points adjacent to a tangent CP" +
                 "move relative to one another in order to make the curve smooth." +
                 "If tangent is true then these have to be set.")]
        public BezierCPComponent adjcp0, adjcp1;
        // this is for the bezier system only
        public Vector3 lp;


        private void Awake()
        {
            if (adjcp0 == null || adjcp1 == null)
                tangent = false;
        }

        // needs to be in here so each item is seperately pickable
        private void OnDrawGizmos()
        {
            Color prevColor = Gizmos.color;
            Gizmos.color = controlPointColor;
            if (tangent)
            {
                if (adjcp0)
                    Gizmos.DrawLine(transform.position, adjcp0.transform.position);
                if (adjcp1)
                    Gizmos.DrawLine(transform.position, adjcp1.transform.position);
            }
            Gizmos.DrawSphere(transform.position, controlPointRadius);
            Gizmos.color = prevColor;
        }
    }
}