using UnityEngine;
using Unity.Entities;
using Unity.Jobs;
using Unity.Collections;

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
        [Range(0.001f, 1)]
        public float drawDelta = 0.1f;
        public BezierCPComponent[] controlPoints;
    }
}
