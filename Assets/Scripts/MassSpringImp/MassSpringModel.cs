using System.Collections;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;
using static UnityEngine.ParticleSystem;

namespace MassSpringModelImp
{
    [System.Serializable]
    public class Particle
    {
        public GameObject BindSphere;

        public Vector3 Position;
        public Vector3 Velocity;
        public Vector3 Acceleration;

        public Vector3 AccumulatedForce;
        public Vector3 Normal;

        public Vector2 TexturePosition;

        public Particle(Vector3 position, GameObject bind_sphere, Vector2 texture_position)
        {
            BindSphere = bind_sphere;

            SetPosition(position);
            SetVelocity(Vector3.zero);
            SetAcceleration(Vector3.zero);

            ResetForce();
            ResetNormal();
            TexturePosition = texture_position;
        }

        public void SetPosition(Vector3 position)
        {
            Position = position;
            BindSphere.transform.position = position;
        }

        public void SetVelocity(Vector3 velocity)
        {
            Velocity = velocity;
        }

        public void SetAcceleration(Vector3 acceleration)
        {
            Acceleration = acceleration;
        }

        public void ResetForce()
        {
            AccumulatedForce = Vector3.zero;
        }

        public void ResetNormal()
        {
            Normal = Vector3.zero;
        }
    }

    [System.Serializable]
    public class Spring
    {
        public Particle ParticleA;
        public Particle ParticleB;

        public float RestLength = 1.0f;
        public float Stiffness = 1.0f;

        public Spring(Particle particleA, Particle particleB, float restLength, float stiffness)
        {
            ParticleA = particleA;
            ParticleB = particleB;
            RestLength = restLength;
            Stiffness = stiffness;
        }

        public float GetLength()
        {
            return Vector3.Distance(ParticleA.Position, ParticleB.Position);
        }

        public Vector3 GetForce()
        {
            Vector3 link = ParticleB.Position - ParticleA.Position;
            float length = GetLength();

            return (link - (link / length) * RestLength) * Stiffness;
        }
    }

    [System.Serializable]
    public class Cloth
    {
        public List<List<Particle>> Particles;
        public List<Particle> FixedParticles;

        public List<Spring> StructuralSprings;
        public List<Spring> ShearSprings;
        public List<Spring> FlexionSprings;

        public Cloth(List<List<Particle>> particles, float structural_spring_length = 5.0f, float shear_spring_length = 7.1f, float flexion_spring_length = 10.0f,
            float structural_spring_stiffness = 20.0f, float shear_spring_stiffness = 20.0f, float flexion_spring_stiffness = 20.0f)
        {
            Particles = particles;

            // set fixed particles
            FixedParticles = new List<Particle>();
            FixedParticles.Add(Particles[0][0]);
            Particles[0][0].BindSphere.GetComponent<MeshRenderer>().material.color = Color.red;
            FixedParticles.Add(Particles[Particles.Count-1][0]);
            Particles[Particles.Count - 1][0].BindSphere.GetComponent<MeshRenderer>().material.color = Color.red;

            // build structural springs
            StructuralSprings = new List<Spring>();
            for (int i = 0; i < Particles.Count - 1; i++)
            {
                for (int j = 0; j < Particles[i].Count; j++)
                    StructuralSprings.Add(new Spring(Particles[i][j], Particles[i + 1][j], structural_spring_length, structural_spring_stiffness));
            }
            for (int i = 0; i < Particles.Count; i++)
            {
                for (int j = 0; j < Particles[i].Count - 1; j++)
                    StructuralSprings.Add(new Spring(Particles[i][j], Particles[i][j + 1], structural_spring_length, structural_spring_stiffness));
            }

            // build shear springs
            ShearSprings = new List<Spring>();
            for (int i = 0; i < Particles.Count - 1; i++)
            {
                for (int j = 0; j < Particles[i].Count - 1; j++)
                    ShearSprings.Add(new Spring(Particles[i][j], Particles[i + 1][j + 1], shear_spring_length, shear_spring_stiffness));
            }
            for (int i = 0; i < Particles.Count - 1; i++)
            {
                for (int j = 0; j < Particles[i].Count - 1; j++)
                    ShearSprings.Add(new Spring(Particles[i + 1][j], Particles[i][j + 1], shear_spring_length, shear_spring_stiffness));
            }

            // build flexion springs
            FlexionSprings = new List<Spring>();
            for (int i = 0; i < Particles.Count - 2; i++)
            {
                for (int j = 0; j < Particles[i].Count; j++)
                    FlexionSprings.Add(new Spring(Particles[i][j], Particles[i + 2][j], flexion_spring_length, flexion_spring_stiffness));
            }
            for (int i = 0; i < Particles.Count; i++)
            {
                for (int j = 0; j < Particles[i].Count - 2; j++)
                    FlexionSprings.Add(new Spring(Particles[i][j], Particles[i][j + 2], flexion_spring_length, flexion_spring_stiffness));
            }
        }

        public void Simulation(float particle_mass, Vector3 gravity, Vector3 fluid_velocity, 
            float critical_deformation_rate = 0.1f, float time_step = 0.02f,
            float spring_damping_coefficient = 1.0f, float fluid_damping_coefficient = 1.0f)
        {
            // reset forces of particles
            for (int i = 0; i < Particles.Count; i++)
            {
                for (int j = 0; j < Particles[i].Count; j++)
                {
                    Particles[i][j].ResetForce();
                    Particles[i][j].ResetNormal();
                }
            }

            // compute spring (internal) forces applied on each particle by Hooke's Law
            // -- structural springs
            for (int i = 0; i < StructuralSprings.Count; i++)
            {
                Vector3 force = StructuralSprings[i].GetForce();
                StructuralSprings[i].ParticleA.AccumulatedForce += force;
                StructuralSprings[i].ParticleB.AccumulatedForce += -force;
            }
            // -- shear springs
            for (int i = 0; i < ShearSprings.Count; i++)
            {
                Vector3 force = ShearSprings[i].GetForce();
                ShearSprings[i].ParticleA.AccumulatedForce += force;
                ShearSprings[i].ParticleB.AccumulatedForce += -force;
            }
            // -- flexion springs
            for (int i = 0; i < FlexionSprings.Count; i++)
            {
                Vector3 force = FlexionSprings[i].GetForce();
                FlexionSprings[i].ParticleA.AccumulatedForce += force;
                FlexionSprings[i].ParticleB.AccumulatedForce += -force;
            }

            // compute normal vector of each particle
            for (int i = 0; i < Particles.Count - 1; i++)
            {
                for (int j = 0; j < Particles[i].Count - 1; j++)
                {
                    // the normal of the "upper" triangle defined by the position of three particles
                    Vector3 v1_upper = Particles[i + 1][j + 1].Position - Particles[i][j].Position;
                    Vector3 v2_upper = Particles[i + 1][j].Position - Particles[i][j].Position;
                    //Vector3 normal_upper = Vector3.Cross(v1_upper, v2_upper);
                    Vector3 normal_upper = Vector3.Cross(v1_upper, v2_upper).normalized;

                    // accumulate the normal of upper triangle
                    Particles[i][j].Normal += normal_upper;
                    Particles[i + 1][j].Normal += normal_upper;
                    Particles[i + 1][j + 1].Normal += normal_upper;

                    // the normal of the "lower" triangle defined by the position of another three particles
                    Vector3 v1_lower = Particles[i][j + 1].Position - Particles[i][j].Position;
                    Vector3 v2_lower = Particles[i + 1][j + 1].Position - Particles[i][j].Position;
                    //Vector3 normal_lower = Vector3.Cross(v1_lower, v2_lower);
                    Vector3 normal_lower = Vector3.Cross(v1_lower, v2_lower).normalized;

                    // accumulate the normal of lower triangle
                    Particles[i][j].Normal += normal_lower;
                    Particles[i][j + 1].Normal += normal_lower;
                    Particles[i + 1][j + 1].Normal += normal_lower;
                }
            }

            // external forces applied on each particle
            Vector3 weight = gravity * particle_mass;
            for (int i = 0; i < Particles.Count; i++)
            {
                for (int j = 0; j < Particles[i].Count; j++)
                {
                    Particles[i][j].Normal.Normalize();

                    Vector3 viscous_damping = -Particles[i][j].Velocity * spring_damping_coefficient;
                    Vector3 viscous_fluid_force = Particles[i][j].Normal * Vector3.Dot(Particles[i][j].Normal, fluid_velocity - Particles[i][j].Velocity) * fluid_damping_coefficient;
                    Particles[i][j].AccumulatedForce = Particles[i][j].AccumulatedForce  + weight + viscous_damping + viscous_fluid_force;
                }
            }

            // reset force of each fixed particle
            for(int i = 0; i < FixedParticles.Count; i++) FixedParticles[i].ResetForce();

            // update positions by using Euler's Method
            for (int i = 0; i < Particles.Count; i++)
            {
                for (int j = 0; j < Particles[i].Count; j++)
                {
                    Particles[i][j].SetAcceleration(Particles[i][j].AccumulatedForce / particle_mass);
                    Particles[i][j].SetVelocity(Particles[i][j].Velocity + Particles[i][j].Acceleration * time_step);
                    Particles[i][j].SetPosition(Particles[i][j].Position + Particles[i][j].Velocity * time_step);
                }
            }

            // apply dynbamic inverse procedure
            // -- structual spring
            for (int i = 0; i < StructuralSprings.Count; i++)
            {
                float length = StructuralSprings[i].GetLength();
                float deformationRate = (length - StructuralSprings[i].RestLength) / StructuralSprings[i].RestLength;
                if (deformationRate > critical_deformation_rate)
                {
                    // case 1: spring linking a fixed particle and a loose particle
                    if (FixedParticles.Contains(StructuralSprings[i].ParticleA))
                    {
                        float m = 1.1f * StructuralSprings[i].RestLength;
                        float n = length - m;
                        Vector3 newPosition = (StructuralSprings[i].ParticleA.Position * n + StructuralSprings[i].ParticleB.Position * m) / length;
                        StructuralSprings[i].ParticleB.SetPosition(newPosition);
                    }
                    else if (FixedParticles.Contains(StructuralSprings[i].ParticleB))
                    {
                        float m = 1.1f * StructuralSprings[i].RestLength;
                        float n = length - m;
                        Vector3 newPosition = (StructuralSprings[i].ParticleB.Position * n + StructuralSprings[i].ParticleA.Position * m) / length;
                        StructuralSprings[i].ParticleA.SetPosition(newPosition);
                    }
                    // case 2: spring linking two loose particles
                    else
                    {
                        Vector3 middle = (StructuralSprings[i].ParticleA.Position + StructuralSprings[i].ParticleB.Position) / 2.0f;
                        float n = 0.55f * StructuralSprings[i].RestLength;
                        float m = 0.5f * length - n;

                        // update the positions
                        StructuralSprings[i].ParticleA.SetPosition((StructuralSprings[i].ParticleA.Position * n + middle * m) / (0.5f * length));
                        StructuralSprings[i].ParticleB.SetPosition((StructuralSprings[i].ParticleB.Position * n + middle * m) / (0.5f * length));
                    }
                }
            }

            // -- shear springs
            for (int i = 0; i < ShearSprings.Count; i++)
            {
                float length = ShearSprings[i].GetLength();
                float deformationRate = (length - ShearSprings[i].RestLength) / ShearSprings[i].RestLength;
                if (deformationRate > critical_deformation_rate)
                {
                    // case 1: spring linking a fixed particle and a loose particle
                    if (FixedParticles.Contains(ShearSprings[i].ParticleA))
                    {
                        float m = 1.1f * ShearSprings[i].RestLength;
                        float n = length - m;
                        Vector3 newPosition = (ShearSprings[i].ParticleA.Position * n + ShearSprings[i].ParticleB.Position * m) / length;
                        ShearSprings[i].ParticleB.SetPosition(newPosition);
                    }
                    else if (FixedParticles.Contains(ShearSprings[i].ParticleB))
                    {
                        float m = 1.1f * ShearSprings[i].RestLength;
                        float n = length - m;
                        Vector3 newPosition = (ShearSprings[i].ParticleB.Position * n + ShearSprings[i].ParticleA.Position * m) / length;
                        ShearSprings[i].ParticleA.SetPosition(newPosition);
                    }
                    // case 2: spring linking two loose particles
                    else
                    {
                        Vector3 middle = (ShearSprings[i].ParticleA.Position + ShearSprings[i].ParticleB.Position) / 2.0f;
                        float n = 0.55f * ShearSprings[i].RestLength;
                        float m = 0.5f * length - n;

                        // update the positions
                        ShearSprings[i].ParticleA.SetPosition((ShearSprings[i].ParticleA.Position * n + middle * m) / (0.5f * length));
                        ShearSprings[i].ParticleB.SetPosition((ShearSprings[i].ParticleB.Position * n + middle * m) / (0.5f * length));
                    }
                }
            }
        }

    }

    public enum DebugMode
    {
        None,
        AllSprings,
        StructuralSpringOnly,
        FlexionSpringOnly,
        ShearSpringOnly
    }

    public class MassSpringModel : MonoBehaviour
    {
        private Cloth Cloth;

        [Header("[Settings]")]
        public int MeshWidth = 10;
        public int MeshHeight = 10;
        public float InitialSpringLength = 5.0f;

        public float ParticleMass = 1.0f;
        public Vector3 Gravity = new Vector3(0.0f, -9.8f, 0.0f);
        public Vector3 FluidVelocity = new Vector3(0.0f, 0.0f, -30.0f);

        [Header("[Particle]")]
        public bool isShowParticles = true;
        public List<GameObject> ParticleSpheres;
        public Transform ParticleSpheresRoot;

        [Header("[Mesh]")]
        public bool isUseImageMaterial = true;
        public Mesh ClothMesh;
        public MeshFilter ClothMeshFilter;
        public MeshRenderer ClothMeshRenderer;
        public Material MeshMaterial;
        public Material MeshImageMaterial;

        [Header("[Debug]")]
        public DebugMode DebugMode = DebugMode.None;

        private void Awake()
        {
            ParticleSpheres = new List<GameObject>();
            ParticleSpheresRoot = new GameObject("ParticleSpheresRoot").transform;
            ParticleSpheresRoot.SetParent(transform);

            List<List<Particle>> particles = new List<List<Particle>>();
            Vector2[] uv = new Vector2[MeshWidth * MeshHeight];
            for (int i = 0; i < MeshWidth; i++)
            {
                List<Particle> column_particles = new List<Particle>();
                Vector3 row_pos = new Vector3(InitialSpringLength * i, 0.0f, 0.0f);
                for (int j = 0; j < MeshHeight; j++)
                {
                    GameObject particle_sphere = Instantiate(PrefabsManager.Instance.ParticleSphere);
                    particle_sphere.transform.SetParent(ParticleSpheresRoot);
                    ParticleSpheres.Add(particle_sphere);

                    Vector3 column_pos = new Vector3(0.0f, -InitialSpringLength * j, 0.0f);
                    Particle new_particle = new Particle(row_pos + column_pos, particle_sphere, new Vector2(j / (MeshHeight - 1.0f), i / (MeshWidth - 1.0f)));
                    column_particles.Add(new_particle);

                    uv[j * MeshWidth + i] = new Vector2((float)i / (MeshWidth - 1), (float)j / (MeshHeight - 1));
                }
                particles.Add(column_particles);
            }

            Cloth = new Cloth(particles);

            // create mesh
            ClothMesh = new Mesh();
            ClothMesh.name = "ClothMesh";
            ClothMeshFilter = gameObject.AddComponent<MeshFilter>();
            ClothMeshFilter.mesh = ClothMesh;
            ClothMeshRenderer = gameObject.AddComponent<MeshRenderer>();

            UpdateClothMesh();
            ClothMesh.uv = uv;
        }

        private void UpdateClothMesh()
        {
            // set vertices
            Vector3[] vertices = new Vector3[MeshWidth * MeshHeight];
            for (int i = 0; i < MeshWidth; i++)
            {
                for (int j = 0; j < MeshHeight; j++)
                {
                    vertices[j * MeshWidth + i] = Cloth.Particles[i][j].Position;
                }
            }

            // set indexes of each triangle
            int[] triangles = new int[(MeshWidth - 1) * (MeshHeight - 1) * 6];
            int t = 0;
            for (int i = 0; i < MeshWidth - 1; i++)
            {
                for (int j = 0; j < MeshHeight - 1; j++)
                {
                    int k = j * MeshWidth + i;

                    triangles[t++] = k;
                    triangles[t++] = k + MeshWidth;
                    triangles[t++] = k + 1;

                    triangles[t++] = k + 1;
                    triangles[t++] = k + MeshWidth;
                    triangles[t++] = k + MeshWidth + 1;
                }
            }

            ClothMesh.vertices = vertices;
            ClothMesh.triangles = triangles;
            ClothMesh.RecalculateNormals();
            ClothMeshRenderer.material = (isUseImageMaterial) ? MeshImageMaterial : MeshMaterial;
        }

        private void Update()
        {
            ParticleSpheresRoot.gameObject.SetActive(isShowParticles);

            Cloth.Simulation(ParticleMass, Gravity, FluidVelocity);
            UpdateClothMesh();
        }

        private void OnDrawGizmos()
        {
            if (Application.isPlaying)
            {
                switch (DebugMode)
                {
                    case DebugMode.AllSprings:
                        DebugStructuralSprings();
                        DebugShearSprings();
                        DebugFlexionSprings();
                        break;
                    case DebugMode.StructuralSpringOnly:
                        DebugStructuralSprings();
                        break;
                    case DebugMode.ShearSpringOnly:
                        DebugShearSprings();
                        break;
                    case DebugMode.FlexionSpringOnly:
                        DebugFlexionSprings();
                        break;
                }
            }
        }

        private void DebugStructuralSprings()
        {
            // debug structural springs
            for (int i = 0; i < Cloth.StructuralSprings.Count; i++)
            {
                Debug.DrawLine(Cloth.StructuralSprings[i].ParticleA.Position, Cloth.StructuralSprings[i].ParticleB.Position, Color.white);
            }
        }

        private void DebugShearSprings()
        {
            // debug shear springs
            for (int i = 0; i < Cloth.ShearSprings.Count; i++)
            {
                Debug.DrawLine(Cloth.ShearSprings[i].ParticleA.Position, Cloth.ShearSprings[i].ParticleB.Position, Color.yellow);
            }
        }

        private void DebugFlexionSprings()
        {
            // debug flexion springs
            for (int i = 0; i < Cloth.FlexionSprings.Count; i++)
            {
                Debug.DrawLine(Cloth.FlexionSprings[i].ParticleA.Position, Cloth.FlexionSprings[i].ParticleB.Position, Color.red);
            }
        }

        
    }
}

