using UnityEngine;

[ExecuteInEditMode]
public class SwingTwist : MonoBehaviour
{
    public Transform other;
    [Range(0, 1)]
    public float swingParam;

    private Quaternion swing;
    private Quaternion qi;

    private void OnEnable()
    {
        Vector3 cross = Vector3.Cross(transform.up, other.up);
        float theta = Vector3.SignedAngle(transform.up, other.up, cross);
        swing = Quaternion.AngleAxis(-theta, cross);

        // now try to get the twist...
        Quaternion q = swing * other.rotation;
        Quaternion qr = Quaternion.Inverse(transform.rotation) * q;

        // now, apply the twist around others local y...
        Quaternion qt = Quaternion.AngleAxis(-qr.eulerAngles.y, other.transform.up);
        other.rotation = qt * other.rotation;
        qi = other.rotation;

        swingParam = 0;
    }
    

    private void Update()
    {
        if (!other)
            return;
        other.rotation = Quaternion.Slerp(Quaternion.identity, swing, swingParam) * qi;
    }
}
