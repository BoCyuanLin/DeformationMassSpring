using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace MassSpringModelImp
{
    public class PrefabsManager : MonoBehaviour
    {
        // private the constructor
        private PrefabsManager() { }

        // claim static instance object
        private static PrefabsManager instance;

        // getter (lazy)
        public static PrefabsManager Instance
        {
            get
            {
                if (instance != null) return instance;

                instance = FindObjectOfType<PrefabsManager>();
                if (instance != null) return instance;

                instance = new PrefabsManager();
                return instance;
            }
        }

        public GameObject ParticleSphere;
    }
}

