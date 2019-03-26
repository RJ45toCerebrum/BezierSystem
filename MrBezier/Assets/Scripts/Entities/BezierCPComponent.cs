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
        [System.Serializable]
        public struct TangentInfo
        {
            [Tooltip("This Tells the Bezier system to make the curve at this control point smooth." +
            "Make Sure that if this is check then you must have the adjcp0 and adjcp1 references set")]
            public bool tangent;
            [Tooltip("This Tells bezier system how the control points adjacent to a tangent CP" +
                     "move relative to one another in order to make the curve smooth." +
                     "If tangent is true then these have to be set.")]
            public BezierCPComponent adjcp0, adjcp1;

            // variables for locking the tangents direction in space
            public bool lockDirection;
            [HideInInspector]
            public Vector3 lockDir;
            [HideInInspector]
            public bool lockDirectionSet;
        };

        public float controlPointRadius = 1;
        public Color controlPointColor;
        public TangentInfo tangentInfo;
        
        // This is for the bezier system only
        private Vector3 lp;
        public Vector3 Lp { get => lp; set => lp = value; }

        private void Awake()
        {
            if (tangentInfo.adjcp0 == null || tangentInfo.adjcp1 == null)
                tangentInfo.tangent = false;
        }

        // needs to be in here so each item is seperately pickable
        private void OnDrawGizmos()
        {
            Gizmos.color = controlPointColor;
            Gizmos.DrawSphere(transform.position, controlPointRadius);
        }
    }
}