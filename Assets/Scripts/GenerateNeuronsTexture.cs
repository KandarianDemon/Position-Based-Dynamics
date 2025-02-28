using UnityEngine;
using System.IO;
using UnityEditor;

public class GenerateNeuronsTexture : MonoBehaviour
{
    public int textureSize = 512;
    public float scale = 10.0f;
    public float branchDensity = 5.0f;
    public float noiseStrength = 0.5f;
    public Color neuronColor = Color.white;
    public Color backgroundColor = Color.black;

    private Texture2D neuronsTexture;

    void Start()
    {
        neuronsTexture = new Texture2D(textureSize, textureSize);
        GenerateTexture();
        //ApplyToMaterial();
        SaveTextureAsPNG(neuronsTexture, "C:\\Users\\Tassix\\Desktop\\neuronsTexture.png");
    }

    void GenerateTexture()
    {
        for (int y = 0; y < textureSize; y++)
        {
            for (int x = 0; x < textureSize; x++)
            {
                float nx = (float)x / textureSize * scale;
                float ny = (float)y / textureSize * scale;

                // Use Perlin noise to create organic neuron branches
                float noise = Mathf.PerlinNoise(nx, ny);
                float branchPattern = Mathf.Sin(branchDensity * nx) * Mathf.Sin(branchDensity * ny);

                // Combine Perlin noise and sinusoidal patterns
                float neuronIntensity = Mathf.Clamp01((noise * noiseStrength + branchPattern * 0.5f));

                // Blend between neuron color and background
                Color pixelColor = Color.Lerp(backgroundColor, neuronColor, neuronIntensity);
                neuronsTexture.SetPixel(x, y, pixelColor);
            }
        }
        neuronsTexture.Apply();
    }

    void ApplyToMaterial()
    {
        Renderer renderer = GetComponent<Renderer>();
        if (renderer != null && renderer.material != null)
        {
            renderer.material.SetTexture("_NeuronsTex", neuronsTexture);
        }
    }
    
    void SaveTextureAsPNG(Texture2D texture, string path)
    {
        byte[] bytes = texture.EncodeToPNG();
        File.WriteAllBytes(path, bytes);

        // Refresh the asset database so the file appears in Unity
        AssetDatabase.ImportAsset(path);
        Debug.Log("Saved neuron texture to " + path);
    }
}
