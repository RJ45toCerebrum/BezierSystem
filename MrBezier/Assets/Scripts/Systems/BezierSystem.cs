using UnityEngine;
using Unity.Entities;
using UnityEditor;
using System.Collections.Generic;

namespace EvolveVR.Bezier
{
    // The Main system for computing Bezier Algorithms
    // This system and corresponding components are in 
    // the process of being converted to Unity ECS

    [ExecuteAlways]
    public class BezierSystem : ComponentSystem
    {
        // Filters for finding components
        public struct BezierFilter
        {
            public BezierComponent bezierComponent;
        };

        public struct CPFilter
        {
            public BezierCPComponent cp;
        };

        protected override void OnStartRunning()
        {
            base.OnStartRunning();
            
            // project transforms and set the components lasp position state
            BezierCPComponent cp;
            foreach(var entity in GetEntities<CPFilter>()) {
                cp = entity.cp;
                cp.lp = cp.transform.localPosition;
            }
        }

        // For Bezier Control Point movement
        protected override void OnUpdate()
        {
            foreach (var entity in GetEntities<CPFilter>())
            {
                BezierCPComponent cp = entity.cp;
                if (!cp.tangent)
                    continue;

                // if this transform moves then update neighbors
                if (ControlPointMoved(cp))
                {
                    Vector3 v = cp.transform.localPosition - cp.lp;
                    cp.adjcp0.transform.localPosition += v;
                    cp.adjcp0.lp = cp.adjcp0.transform.localPosition;
                    cp.adjcp1.transform.localPosition += v;
                    cp.adjcp1.lp = cp.adjcp1.transform.localPosition;
                    cp.lp = cp.transform.localPosition;
                }
                // if neighbor 0 moves then update neighbor 1 position
                else if (ControlPointMoved(cp.adjcp0))
                {
                    Vector3 v = cp.transform.localPosition - cp.adjcp0.transform.localPosition;
                    cp.adjcp1.transform.localPosition = cp.transform.localPosition + v;

                    cp.adjcp1.lp = cp.adjcp1.transform.localPosition;
                    cp.adjcp0.lp = cp.adjcp0.transform.localPosition;
                    cp.lp = cp.transform.localPosition;
                }
                // opposite of above
                else if (ControlPointMoved(cp.adjcp1))
                {
                    Vector3 v = cp.transform.localPosition - cp.adjcp1.transform.localPosition;
                    cp.adjcp0.transform.localPosition = cp.transform.localPosition + v;
                    cp.adjcp1.lp = cp.adjcp1.transform.localPosition;
                    cp.adjcp0.lp = cp.adjcp0.transform.localPosition;
                }
            }
        }

        // For drawing the curve...
        [DrawGizmo(GizmoType.NonSelected | GizmoType.Selected)]
        public static void DrawGizmoForCurve(BezierComponent bc, GizmoType gizmoType)
        {
            BezierCPComponent[] CPS = bc.controlPoints;
            if (CPS == null || CPS.Length == 0)
                return;

            int N = CurveCount(bc);

            Vector3 cp;
            Vector3 pp = BezierSystem.Evaluate(0, bc);
            Gizmos.color = bc.curveColor;
            float t;
            for (t = bc.drawDelta; t < N; t += bc.drawDelta)
            {
                cp = BezierSystem.Evaluate(t, bc);
                Gizmos.DrawLine(pp, cp);
                pp = cp;
            }

            cp = BezierSystem.Evaluate(N-0.0002f, bc);
            Gizmos.DrawLine(pp, cp);
        }
        
        // did bezier CP move from last set state?
        protected static bool ControlPointMoved(BezierCPComponent cp)
        {
            if (!cp)
                throw new System.ArgumentException("Null ref. exception. Check to see if reference is missing in BezierCPComponent Script...");
            return !EqualPoint(cp.transform.localPosition, cp.lp);
        }

        protected static bool EqualPoint(Vector3 p0, Vector3 p1)
        {
            if (Mathf.Abs(p0.x - p1.x) > 0.0001f)
                return false;
            if (Mathf.Abs(p0.y - p1.y) > 0.0001f)
                return false;
            if (Mathf.Abs(p0.z - p1.z) > 0.0001f)
                return false;
            return true;
        }

        /// <summary>
        /// Gives the position of the curve in World Space at parameter t
        /// </summary>
        /// <param name="t">
        /// parameter. Domain -> [0, N)
        /// where N is the number of curves appended together
        /// to give the full curve. The BezierCPComponent[]
        /// length must be of form: 4*i where is a natural number...
        /// </param>
        /// <param name="bc">The BezierComponent to evaluate</param>
        /// <returns>point on curve at param t</returns>
        public static Vector3 Evaluate(float t, BezierComponent bc)
        {
            int startIdx = -1;
            t = Param(t, bc, out startIdx);
            float p = 1 - t;

            BezierCPComponent[] cps = bc.controlPoints;

            Vector3 p0 = cps[startIdx].transform.position;
            Vector3 p1 = cps[startIdx + 1].transform.position;
            Vector3 p2 = cps[startIdx + 2].transform.position;
            Vector3 p3 = cps[startIdx + 3].transform.position;

            return (p0 * p*p*p) + (3*p1*t *p*p) + (3* p2 * p*t*t) + (p3 * t*t*t);
        }

        /// <summary>
        /// Returns the Arc Length of the curve.
        /// Currently, its mostly a bruteforce approach and therefore
        /// the paramter dt is important for efficiency. If you want
        /// to evaulate the entire curve, use FastArcLength instead.
        /// </summary>
        /// <param name="t0">starting param</param>
        /// <param name="t1">end param</param>
        /// <param name="bc">BezierComponent to eval.</param>
        /// <param name="dt">delta param. This controls the steps algo takes.
        /// Its important to understand the lower the number, the more accurate
        /// but the more expensive it becomes.
        /// </param>
        /// <returns>The Arclength of curve between t0 and t1</returns>
        public static float ArcLength(float t0, float t1, BezierComponent bc, float dt = 0.01f)
        {
            float length = 0;
            float t = t0 + dt;
            Vector3 pp = BezierSystem.Evaluate(t0, bc);
            Vector3 cp;
            while (t < t1)
            {
                cp = BezierSystem.Evaluate(t, bc);
                length += (cp - pp).magnitude;
                pp = cp;
                t += dt;
            }

            cp = BezierSystem.Evaluate(t1, bc);
            length += (cp - pp).magnitude;
            return length;
        }

        /// <summary>
        /// Better to call this when evaluating the arclength for the entire curve...
        /// </summary>
        /// <param name="bezierData">BezierComponent data to evaluate</param>
        /// <returns>Arc Length</returns>
        public static float FastArcLength(BezierComponent bc)
        {
            float L = 0.0f;
            int N = (bc.controlPoints.Length - bc.controlPoints.Length % 4) / 4;
            float arclength = 0f;
            for (int i = 0; i < N; i++)
            {
                Vector3 A = bc.controlPoints[4*i].transform.position;
                Vector3 B = bc.controlPoints[4*i + 1].transform.position;
                Vector3 C = bc.controlPoints[4*i + 2].transform.position;
                Vector3 D = bc.controlPoints[4*i + 3].transform.position;
                L = 0;
                ArcLengthUtil(A, B, C, D, 5, ref L);
                arclength += L;
            }
            return arclength;
        }

        // util for faster arclength calc of entire curve....
        protected static void ArcLengthUtil(Vector3 A, Vector3 B, Vector3 C, Vector3 D, uint subdiv, ref float L)
        {
            if (subdiv > 0)
            {
                Vector3 a = A + (B - A) * 0.5f;
                Vector3 b = B + (C - B) * 0.5f;
                Vector3 c = C + (D - C) * 0.5f;
                Vector3 d = a + (b - a) * 0.5f;
                Vector3 e = b + (c - b) * 0.5f;
                Vector3 f = d + (e - d) * 0.5f;

                // left branch
                ArcLengthUtil(A, a, d, f, subdiv - 1, ref L);
                // right branch
                ArcLengthUtil(f, e, c, D, subdiv - 1, ref L);
            }
            else // converges quickly as subdivisions increase...
            {
                float controlNetLength = (B - A).magnitude + (C - B).magnitude + (D - C).magnitude;
                float chordLength = (D - A).magnitude;
                L += (chordLength + controlNetLength) / 2.0f;
            }
        }

        /// <summary>
        /// Gives the approximate closest point to point p. dt controls how accurate it will be.
        /// </summary>
        /// <param name="p">Point we are getting closest point on curve to</param>
        /// <param name="bc">Bezier Component</param>
        /// <param name="dt">delta controls number of iterations taken
        /// and is therefore important to understand efficiency is reduced
        /// as dt gets smaller...
        /// </param>
        /// <returns>The closesst point along curve to p</returns>
        public static Vector3 ClosestPoint(Vector3 p, BezierComponent bc, float dt = 0.01f)
        {
            if (dt < 0.000001f)
                return Vector3.positiveInfinity;

            int C = bc.controlPoints.Length;
            int N = (C - C % 4) / 4;

            Vector3 ccp = Evaluate(0, bc);
            float csd = (p - ccp).sqrMagnitude;
            float t = dt;
            Vector3 cp;
            float d;
            while (t < N)
            {
                cp = Evaluate(t, bc);
                d = (p - cp).sqrMagnitude;
                if (d < csd) {
                    ccp = cp;
                    csd = d;
                }
                t += dt;
            }

            return ccp;
        }

        /// <summary>
        /// Gives the Tangent vector to the curve
        /// </summary>
        /// <param name="t">param</param>
        /// <param name="bc">BezierComponent to evaluate</param>
        /// <returns></returns>
        public static Vector3 TangentVector(float t, BezierComponent bc)
        {
            int startIdx = 0;
            t = Param(t, bc, out startIdx);

            BezierCPComponent[] cps = bc.controlPoints;
            Vector3 p0 = cps[startIdx].transform.position;
            Vector3 p1 = cps[startIdx + 1].transform.position;
            Vector3 p2 = cps[startIdx + 2].transform.position;
            Vector3 p3 = cps[startIdx + 3].transform.position;
        
            Vector3 A = -3 * p0 + 9 * p1 - 9 * p2 + 3 * p3;
            Vector3 B = 6 * p0 + -12 * p1 + 6 * p2;
            Vector3 C = -3 * p0 + 3 * p1;
            Vector3 tangent = (t*t*A) + (t*B) + C;
            return tangent;
        }

        /// <summary>
        /// Gives the next param after t such that the the Length along curve
        /// going from t to retured value is approximately arc length length.
        /// i.e for getting a uniform speed along curve.
        /// </summary>
        /// <param name="t">param</param>
        /// <param name="length">length to travel along curve</param>
        /// <param name="crv">BezierComponent to evaluate</param>
        /// <returns>new uniform length parameter</returns>
        public static float UniformSpeedParam(float t, float length, BezierComponent crv)
        {
            Vector3 p = Evaluate(t, crv);
            Vector3 tangent = TangentVector(t, crv);
            t += length / tangent.magnitude;
            return t;
        }

        /// <summary>
        /// Gives the current count of curves.
        /// Entire curve is made up of individual 
        /// cubic bezier curves.
        /// </summary>
        /// <param name="bc">BezierComponent to eval</param>
        /// <returns>Curve count</returns>
        public static int CurveCount(BezierComponent bc) {
            return (bc.controlPoints.Length - bc.controlPoints.Length % 4) / 4;
        }

        /// <summary>
        /// Transforms the parameter to be between [0,1]
        /// and sets the startIdx so that we can index into
        /// BezierCurve array and get parameter
        /// </summary>
        /// <param name="t">param from [0, CurveCount(BezierComponent) )</param>
        /// <param name="bc">BezierComponent to eval</param>
        /// <param name="startIdx">the start index into the BezierComponent.controlPoints array...</param>
        /// <returns>returns normalized parameter between [0,1]</returns>
        public static float Param(float t, BezierComponent bc, out int startIdx)
        {
            int p = (int)Mathf.Floor(t);
            int N = CurveCount(bc);
            if (p >= N)
            {
                startIdx = -1;
                return -1;
            }
            else
                startIdx = p * 4;

            return t - p;
        }
    }
    

    [ExecuteAlways]
    public class PlanarBezierSystem : BezierSystem
    {
        // Planar filter...
        public struct PlanarBezierFilter
        {
            public PlanarBezierComponent comp;
        };

        // for projected points on start
        protected override void OnStartRunning()
        {
            base.OnStartRunning();
            ComponentGroupArray<PlanarBezierFilter> bezierEntities = GetEntities<PlanarBezierFilter>();
            
            // Project Transforms...
            if (bezierEntities.Length > 0)
                ProjectCPs(bezierEntities);
        }

        // project points onto plane
        protected override void OnUpdate()
        {
            base.OnUpdate();
            ComponentGroupArray<PlanarBezierFilter> bezierEntities = GetEntities<PlanarBezierFilter>();
            if (bezierEntities.Length > 0)
                ProjectCPs(bezierEntities);
        }


        [DrawGizmo(GizmoType.NonSelected | GizmoType.Selected)]
        public static void DrawProjectionPlane(PlanarBezierComponent comp, GizmoType gizmoType)
        {
            Vector3 c = comp.planeTransform.position;
            Vector3 r = 0.5f * comp.planeTransform.right * comp.planeTransform.localScale.x;
            Vector3 f = 0.5f * comp.planeTransform.forward * comp.planeTransform.localScale.z;

            Gizmos.DrawLine(c - r - f, c + r - f);
            Gizmos.DrawLine(c - r + f, c + r + f);
            Gizmos.DrawLine(c - r - f, c - r + f);
            Gizmos.DrawLine(c + r - f, c + r + f);

            Gizmos.color = Color.black;
            Gizmos.DrawSphere(comp.planeTransform.position, 0.1f);
        }

        private static void ProjectCPs(ComponentGroupArray<PlanarBezierFilter> bcs)
        {
            Plane plane = new Plane();
            Vector3 projPoint;
            foreach (var entity in bcs)
            {
                // set plane position and normal
                Transform PT = entity.comp.planeTransform;
                plane.SetNormalAndPosition(PT.up, PT.position);
                // now loop through transforms of each control point
                // and project them onto plane
                foreach (BezierCPComponent CP in entity.comp.controlPoints)
                {
                    projPoint = plane.ClosestPointOnPlane(CP.transform.position);
                    CP.transform.position = projPoint;
                    CP.lp = CP.transform.localPosition;
                    CP.transform.rotation = PT.rotation;
                }

                // for debuging
                //if (entity.comp.debugSphere)
                //{
                //    entity.comp.debugSphere.position = plane.ClosestPointOnPlane(entity.comp.debugSphere.position);
                //    entity.comp.debugSphere.rotation = entity.comp.planeTransform.rotation;
                //    Debug.DrawLine(PT.position, entity.comp.debugSphere.position, Color.red);
                //}
            }
        }

        /// <summary>
        /// Project the point p, onto plane of the PlanarBezierComponent
        /// </summary>
        /// <param name="p">point to project</param>
        /// <param name="pbc">PlanarBezierComponent</param>
        /// <returns>Point on the plane closest to p</returns>
        public static Vector3 ProjectPoint(Vector3 p, PlanarBezierComponent pbc) {
            Plane plane = new Plane(pbc.planeTransform.up, pbc.planeTransform.position);
            return plane.ClosestPointOnPlane(p);
        }

        /// <summary>
        /// Given a point on the plane transform, we draw a vector from center of the plane transform
        /// to the pop (point on plane) then use this to intersect the curve...
        /// </summary>
        /// <param name="pop">point on plane</param>
        /// <param name="bc">PlanarBezierComponent data</param>
        /// <returns>Array of vectors that intersect curve with ray</returns>
        public static Vector3[] IntersectCurve(Vector3 pop, PlanarBezierComponent bc)
        {
            List<Vector3> intersectionPoints = new List<Vector3>(2);

            // 1) transform the curve and pop into plane space
            Transform PT = bc.planeTransform;
            Vector3 planeSpacePop = PT.InverseTransformPoint(pop);
            // only Y is normal of plane for now...
            Vector2 v = new Vector2(planeSpacePop.x, planeSpacePop.z);
            Vector2[] tpts = GetPlaneSpaceCPs(bc.controlPoints, PT);
            // do an additional transformation such that crv is a function of y
            // with respect to t, the param of the ray from origin in dir v
            float angle = -Angle(v);
            RotatePoints(angle, tpts);
            // at most three roots per cubic bezier...
            float[] roots = new float[3];

            int CC = CurveCount(bc);
            for (int i = 0; i < CurveCount(bc); i++)
            {
                if (Roots(tpts, 4*i, ref roots))
                {
                    if(roots[0] > 0) {
                        Vector3 pofi = Evaluate(i + roots[0], bc);
                        intersectionPoints.Add(pofi);
                    }
                    if (roots[1] > 0) {
                        Vector3 pofi = Evaluate(i + roots[1], bc);
                        intersectionPoints.Add(pofi);
                    }
                    if (roots[2] > 0) {
                        Vector3 pofi = Evaluate(i + roots[2], bc);
                        intersectionPoints.Add(pofi);
                    }
                }
            }

            return intersectionPoints.ToArray();
        }

        // find roots analytically: http://mathworld.wolfram.com/CubicFormula.html
        // my prototype here: https://www.desmos.com/calculator/x85qtkgskv
        private static bool Roots(Vector2[] pts, int startidx, ref float[] roots)
        {
            roots[0] = -1; roots[1] = -1; roots[2] = -1;

            Vector2 pa = pts[startidx], pb = pts[startidx+1], pc = pts[startidx+2], pd = pts[startidx+3];
            float a = -pa.y + 3*pb.y - 3*pc.y + pd.y;
            float b = 3*pa.y - 6*pb.y + 3*pc.y;
            float c = -3*pa.y + 3*pb.y;
            float d = pa.y;

            float A = b / a;
            float B = c / a;
            float C = d / a;
            float Q = (3*B - A*A) / 9;
            float R = (9 * A * B - 27 * C - 2 * A * A * A) / 54;
            
            // D > 0 = ONE root real and 2 complex
            // D < 0 = ALL roots are real
            // D == 0 = ALL real roots and two are equal (x axis tangent to curve)
            float D = Mathf.Pow(Q, 3) + R*R;
            if (D < 0)
            {
                float T = Mathf.Acos(R / Mathf.Sqrt(-Q*Q*Q));
                float m = 2 * Mathf.Sqrt(-Q);
                float a_3 = A / 3;

                float r0 = m * Mathf.Cos(T/3) - a_3;
                float r1 = m * Mathf.Cos( (T + 2*Mathf.PI) / 3) - a_3;
                float r2 = m * Mathf.Cos( (T + 4*Mathf.PI) / 3) - a_3;

                bool anyRoots = false;
                if (r0 >= 0 && r0 <= 1){
                    roots[0] = r0;
                    anyRoots = true;
                }
                if (r1 >= 0 && r1 <= 1) {
                    roots[1] = r1;
                    anyRoots = true;
                }
                if (r2 >= 0 && r2 <= 1) {
                    roots[2] = r2;
                    anyRoots = true;
                }
            
                return anyRoots;
            }
            else
            {
                float RpD = R + Mathf.Sqrt(D);
                float RmD = R - Mathf.Sqrt(D);
                float J = Mathf.Sign(RpD) * Mathf.Pow(Mathf.Abs(RpD), 0.33333f);
                float K = Mathf.Sign(RmD) * Mathf.Pow(Mathf.Abs(RmD), 0.33333f);

                float r = -A/3 + J + K;
                if (r >= 0 && r <= 1){
                    roots[0] = r;
                    return true;
                }
            }

            return false;
        }

        // PT = plane transform; no y component.
        private static Vector2[] GetPlaneSpaceCPs(BezierCPComponent[] cps, Transform PT)
        {
            Vector2[] transformedCPS = new Vector2[cps.Length];
            Vector3 p;
            for(int i = 0; i < transformedCPS.Length; i++){
                p = PT.InverseTransformPoint(cps[i].transform.position);
                transformedCPS[i].x = p.x;
                transformedCPS[i].y = p.z;
            }
            return transformedCPS;
        }

        // given v, we get angle from [1,0] to v were the range is [0, 2PI]
        private static float Angle(Vector2 v)
        {
            float theta = Vector3.SignedAngle(Vector3.right, v, Vector3.forward);
            if (theta < 0)
                return (360 - Mathf.Abs(theta)) * Mathf.Deg2Rad;
            return theta * Mathf.Deg2Rad;
        }

        private static void RotatePoints(float angle, Vector2[] pts)
        {
            for(int i = 0; i < pts.Length; i++) {
                pts[i] = new Vector2(pts[i].x * Mathf.Cos(angle) - pts[i].y * Mathf.Sin(angle),
                                     pts[i].x * Mathf.Sin(angle) + pts[i].y * Mathf.Cos(angle));
            }
        }

        // internal helper for evaluating curve with only Vectors
        private static Vector2 Evaluate(float t, Vector2 pa, Vector2 pb, Vector2 pc, Vector2 pd)
        {
            float p = 1 - t;
            return (pa * p * p * p) + (3 * pb * p * p * t) + (3 * pc * p * t * t) + (pd * t * t * t);
        }
    }
}