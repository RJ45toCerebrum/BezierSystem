using Unity.Mathematics;
using Unity.Entities;


[System.Serializable]
public struct PBC : IComponentData
{
    public float3 position;
}

public class PBCComponent : ComponentDataProxy<PBC> { }

