using UnityEngine;

// This is in the process of being changed to the new
// Unity Entity-Component-System.
namespace EvolveVR.Bezier
{
    /// <summary>
    /// This is only a data container which is passed into
    /// static methods in the BezierSystem class.
    /// 
    /// All the methods on BezierSystem require you to pass in
    /// the BezierComponent.
    /// </summary>
    public class BezierComponent : MonoBehaviour
    {
        public Color curveColor;
        public Color tangentColor;
        [Range(0.001f, 1)]
        public float drawDelta = 0.1f;

        [Tooltip("The BezierComponent is a Cubic bezier curve. This means that" +
         "indices [0,3] (inclusive) is a curve and [4,7] is a curve. If you make the" +
         " make the 3rd and 4th index the SAME BezierCPComponent reference then the" +
         " curve looks smooth...")]
        public BezierCPComponent[] controlPoints;

        [Tooltip("The connected curves are for generating a graph of connected bezier curves for branching paths...")]
        public BezierComponent[] connectedCurves;
    }
}
