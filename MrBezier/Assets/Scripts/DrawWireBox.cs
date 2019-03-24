using UnityEngine;


[ExecuteInEditMode]
public class DrawWireBox : MonoBehaviour
{
    public Color color = Color.blue;
    public MeshFilter meshFilter;
    public Transform T;

    Vector3[] mvs;
    Vector4[] boundsVerts = new Vector4[8];
    Vector4[] worldBBVerts = new Vector4[8];

    private void OnEnable()
    {
        mvs = meshFilter.sharedMesh.vertices;

        float minX = Mathf.Infinity, minY = Mathf.Infinity, minZ = Mathf.Infinity;
        float maxX = Mathf.NegativeInfinity, maxY = Mathf.NegativeInfinity, maxZ = Mathf.NegativeInfinity;

        foreach (var v in mvs)
        {
            if (v.x < minX)
                minX = v.x;
            else if (v.x > maxX)
                maxX = v.x;

            if (v.y < minY)
                minY = v.y;
            else if (v.y > maxY)
                maxY = v.y;

            if (v.z < minZ)
                minZ = v.z;
            else if (v.z > maxZ)
                maxZ = v.z;
        }

        // Get the bounds vertices...
        // bottom plane
        boundsVerts[0].x = minX; boundsVerts[0].y = minY; boundsVerts[0].z = minZ; boundsVerts[0].w = 1;
        boundsVerts[1].x = minX; boundsVerts[1].y = minY; boundsVerts[1].z = maxZ; boundsVerts[1].w = 1;
        boundsVerts[2].x = maxX; boundsVerts[2].y = minY; boundsVerts[2].z = minZ; boundsVerts[2].w = 1;
        boundsVerts[3].x = maxX; boundsVerts[3].y = minY; boundsVerts[3].z = maxZ; boundsVerts[3].w = 1;
        // top plane
        boundsVerts[4].x = minX; boundsVerts[4].y = maxY; boundsVerts[4].z = minZ; boundsVerts[4].w = 1;
        boundsVerts[5].x = minX; boundsVerts[5].y = maxY; boundsVerts[5].z = maxZ; boundsVerts[5].w = 1;
        boundsVerts[6].x = maxX; boundsVerts[6].y = maxY; boundsVerts[6].z = minZ; boundsVerts[6].w = 1;
        boundsVerts[7].x = maxX; boundsVerts[7].y = maxY; boundsVerts[7].z = maxZ; boundsVerts[7].w = 1;
    }

    private void OnDrawGizmos()
    {
        Gizmos.color = color;
        UpdateBBVerts();

        // draw bottom plane
        Gizmos.DrawLine(worldBBVerts[0], worldBBVerts[1]);
        Gizmos.DrawLine(worldBBVerts[0], worldBBVerts[2]);
        Gizmos.DrawLine(worldBBVerts[1], worldBBVerts[3]);
        Gizmos.DrawLine(worldBBVerts[2], worldBBVerts[3]);
        // draw top plane
        Gizmos.DrawLine(worldBBVerts[4], worldBBVerts[5]);
        Gizmos.DrawLine(worldBBVerts[4], worldBBVerts[6]);
        Gizmos.DrawLine(worldBBVerts[5], worldBBVerts[7]);
        Gizmos.DrawLine(worldBBVerts[6], worldBBVerts[7]);
        // draw left side plane
        Gizmos.DrawLine(worldBBVerts[0], worldBBVerts[4]);
        Gizmos.DrawLine(worldBBVerts[1], worldBBVerts[5]);
        // draw right side plane
        Gizmos.DrawLine(worldBBVerts[2], worldBBVerts[6]);
        Gizmos.DrawLine(worldBBVerts[3], worldBBVerts[7]);
    }

    private void UpdateBBVerts()
    {
        Matrix4x4 W = Matrix4x4.TRS(T.position, T.rotation, T.lossyScale);
        worldBBVerts[0] = W * boundsVerts[0];
        worldBBVerts[1] = W * boundsVerts[1];
        worldBBVerts[2] = W * boundsVerts[2];
        worldBBVerts[3] = W * boundsVerts[3];
        worldBBVerts[4] = W * boundsVerts[4];
        worldBBVerts[5] = W * boundsVerts[5];
        worldBBVerts[6] = W * boundsVerts[6];
        worldBBVerts[7] = W * boundsVerts[7];
    }
}
