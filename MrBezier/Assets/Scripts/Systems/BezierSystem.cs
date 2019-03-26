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
                cp.Lp = cp.transform.localPosition;
            }
        }

        // For Bezier Control Point movement
        protected override void OnUpdate()
        {
            foreach (var entity in GetEntities<CPFilter>())
            {
                BezierCPComponent cp = entity.cp;
                if (!cp.tangentInfo.tangent)
                    continue;
                
                // if this transform moves then update neighbors
                if (ControlPointMoved(cp))
                {
                    Vector3 lp = cp.transform.parent.TransformPoint(cp.Lp);
                    Vector3 v = cp.transform.position - lp;
                    BezierCPComponent.TangentInfo ti = cp.tangentInfo;
                    ti.adjcp0.transform.position += v;
                    ti.adjcp0.Lp = ti.adjcp0.transform.localPosition;
                    ti.adjcp1.transform.position += v;
                    ti.adjcp1.Lp = ti.adjcp1.transform.localPosition;
                    cp.Lp = cp.transform.localPosition;
                }
                // if neighbor 0 moves then update neighbor 1 position
                else if (ControlPointMoved(cp.tangentInfo.adjcp0))
                {
                    Transform A0 = cp.tangentInfo.adjcp0.transform;
                    Transform C = cp.transform;
                    Transform A1 = cp.tangentInfo.adjcp1.transform;

                    if(cp.tangentInfo.lockDirection)
                    {
                        Vector3 v = A0.position - C.position;
                        if (!cp.tangentInfo.lockDirectionSet) {
                            cp.tangentInfo.lockDir = v.normalized;
                            cp.tangentInfo.lockDirectionSet = true;
                        }
                        float d = Vector3.Dot(v, cp.tangentInfo.lockDir);
                        A0.position = C.position + d * cp.tangentInfo.lockDir;
                    }
                    else
                    {
                        Vector3 v = C.position - A0.position;
                        float d = (C.position - A1.position).magnitude;
                        v.Normalize();
                        A1.position = cp.transform.position + d * v;
                        cp.tangentInfo.lockDirectionSet = false;
                    }

                    cp.tangentInfo.adjcp1.Lp = cp.tangentInfo.adjcp1.transform.localPosition;
                    cp.tangentInfo.adjcp0.Lp = cp.tangentInfo.adjcp0.transform.localPosition;
                }
                // opposite of above
                else if (ControlPointMoved(cp.tangentInfo.adjcp1))
                {
                    Transform A0 = cp.tangentInfo.adjcp0.transform;
                    Transform C = cp.transform;
                    Transform A1 = cp.tangentInfo.adjcp1.transform;

                    if (cp.tangentInfo.lockDirection)
                    {
                        Vector3 v = A1.position - C.position;
                        if (!cp.tangentInfo.lockDirectionSet) {
                            cp.tangentInfo.lockDir = v.normalized;
                            cp.tangentInfo.lockDirectionSet = true;
                        }
                        float d = Vector3.Dot(v, cp.tangentInfo.lockDir);
                        A1.position = C.position + d * cp.tangentInfo.lockDir;
                    }
                    else
                    {
                        Vector3 v = C.position - A1.transform.position;
                        float d = (C.position - A0.position).magnitude;
                        v.Normalize();
                        A0.position = C.position + d * v;
                    }

                    cp.tangentInfo.adjcp1.Lp = A1.localPosition;
                    cp.tangentInfo.adjcp0.Lp = A0.localPosition;
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

            // Draw curve
            Vector3 cp;
            Vector3 pp = BezierSystem.Evaluate(0, bc);
            Gizmos.color = bc.curveColor;
            for (float t = bc.drawDelta; t < N; t += bc.drawDelta) {
                cp = BezierSystem.Evaluate(t, bc);
                Gizmos.DrawLine(pp, cp);
                pp = cp;
            }
            cp = BezierSystem.Evaluate(N-0.0002f, bc);
            Gizmos.DrawLine(pp, cp);

            // Draw CP lines
            Gizmos.color = bc.tangentColor;
            for (int i = 0; i < N; i++)
            {
                int s = 4 * i;
                Gizmos.DrawLine(CPS[s].transform.position, CPS[s+1].transform.position);
                Gizmos.DrawLine(CPS[s+2].transform.position, CPS[s + 3].transform.position);
            }
        }
        
        // did bezier CP move from last set state?
        protected static bool ControlPointMoved(BezierCPComponent cp)
        {
            if (!cp)
                throw new System.ArgumentException("Null ref. exception. Check to see if reference is missing in BezierCPComponent Script...");
            return !EqualPoint(cp.transform.localPosition, cp.Lp);
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
            if (p >= N) {
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

            Gizmos.color = Color.green;
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
                    CP.Lp = CP.transform.localPosition;
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

            /*  Steps of algo:
             *  1) Transform the curve and point on plane into plane space
             *  2) Get the ray going from the planeTransform position to the transformed pop
             *  3) Make the curve y value a function of t of the ray by rotating it
             *      by theta radians, which is the angle between [1,0] and v vectors
             *  4) analytic root find the intersections for each curve making the entire curve
             */
            Transform PT = bc.planeTransform;
            Vector3 planeSpacePop = PT.InverseTransformPoint(pop);
            Vector2 v = new Vector2(planeSpacePop.x, planeSpacePop.z);
            Vector2[] tpts = GetPlaneSpaceCPs(bc.controlPoints, PT);
            float angle = -Angle(v);
            RotatePoints(angle, tpts);
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

        public static Vector3 IntersectCurveNearestForward(Vector3 pop, PlanarBezierComponent bc)
        {
            Vector3[] intersections = PlanarBezierSystem.IntersectCurve(pop, bc);
            Vector3 v = (pop - bc.planeTransform.position);
            v.Normalize();

            Vector3 w, k = Vector3.zero;
            float cd = Mathf.Infinity;
            for(int i = 0; i < intersections.Length; i++)
            {
                w = intersections[i] - bc.planeTransform.position;
                float dot = Vector3.Dot(v, w);
                if (dot > 0 && dot < cd) {
                    cd = dot;
                    k = intersections[i];
                }
            }

            return k;
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


    [ExecuteAlways]
    public class RotationConstraintSystem : ComponentSystem
    {
        struct RotationContraintFilter
        {
            public RotationContraintComponent comp;
        };

        protected override void OnUpdate()
        {
            ComponentGroupArray<RotationContraintFilter> bezierEntities = GetEntities<RotationContraintFilter>();
            if (bezierEntities.Length > 0) {
                foreach(RotationContraintFilter entity in bezierEntities)
                    ContrainRotation(entity.comp);
            }
        }

        private void ContrainRotation(RotationContraintComponent rcc)
        {
            Ray ray = new Ray();
            Plane plane = new Plane();

            PlanarBezierComponent pbc = rcc.pbc;
            ray.origin = rcc.transform.position;
            ray.direction = rcc.transform.up;
            plane.SetNormalAndPosition(pbc.planeTransform.up, pbc.planeTransform.position);

            float t;
            if (plane.Raycast(ray, out t))
            {
                Vector3 ip = ray.origin + t * ray.direction;
                Debug.DrawLine(rcc.transform.position, ip, Color.blue);
                Debug.DrawLine(pbc.planeTransform.position, ip, Color.green);
                Debug.DrawLine(pbc.planeTransform.position, rcc.transform.position, Color.green);

                Vector3 poc = PlanarBezierSystem.IntersectCurveNearestForward(ip, pbc);
                Debug.DrawLine(rcc.transform.position, poc, Color.cyan);

                Vector3 v = ip - pbc.planeTransform.position;
                Vector3 w = poc - pbc.planeTransform.position;
                if (v.sqrMagnitude > w.sqrMagnitude)
                    RotateTo(poc, rcc);
            }
        }

        private void RotateTo(Vector3 ip, RotationContraintComponent bc)
        {
            // now swing not entire way but only to the curve
            Vector3 toIp = ip - bc.transform.position;
            toIp.Normalize();
            Quaternion qs = SwingTo(bc.transform.up, toIp);
            bc.transform.rotation = qs * bc.transform.rotation;
        }

        private Quaternion SwingTo(Vector3 v, Vector3 w)
        {
            Vector3 cross = Vector3.Cross(v, w);
            float theta = Vector3.SignedAngle(v, w, cross);
            return Quaternion.AngleAxis(theta, cross);
        }
    }
}