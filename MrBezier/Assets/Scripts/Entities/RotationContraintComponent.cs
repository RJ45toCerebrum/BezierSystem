using UnityEngine;
using EvolveVR.Bezier;

public class RotationContraintComponent : MonoBehaviour
{
    public PlanarBezierComponent pbc;

    private void Awake()
    {
        Rigidbody rb = GetComponent<Rigidbody>();
        if(rb)
            rb.solverIterations += 20;
    }
}
