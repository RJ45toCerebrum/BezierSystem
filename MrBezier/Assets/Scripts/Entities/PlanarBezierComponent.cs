using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace EvolveVR.Bezier
{
    /// <summary>
    /// Only a data container for passing into static methods
    /// of the PlanarBezierSystem class.
    /// 
    /// PlaneBezierComponent control points live on a Plane.
    /// This plane is determined by the planeTransform Y axis.
    /// On Start, this system will find all Planar Bezier CP's
    /// and project them onto the plane at the position of the
    /// plane transform and with normal planeTransform.up.
    /// 
    /// This Advantge of using this System is that it has useful methods
    /// for finding the intersection of a vectors with the curve...
    /// </summary>
    public class PlanarBezierComponent : BezierComponent
    {
        // For now the Y-Axis is always normal of the plane; due to time contraints =(
        // Also, the angle of intersection is always measured from the right axis of the plane transform as well
        public Transform planeTransform;
        public Color drawColor;
    }
}