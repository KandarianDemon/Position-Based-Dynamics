Shader "Unlit/BackfaceCullingOff_SpecularPerlin"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
        _FrontColor ("Front Color", Color) = (1,1,1,1)
        _BackColor ("Back Color", Color) = (1,1,1,1)
        _SpecColor ("Specular Color", Color) = (1,1,1,1)
        _Glossiness ("Glossiness", Range(0, 1)) = 0.5
        _BumpStrength ("Bump Strength", Range(0, 1)) = 0.2
        _NoiseScale ("Noise Scale", Range(0.1, 10)) = 1.0
        _NoiseDetail ("Noise Detail (Octaves)", Range(1, 5)) = 3
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" }
        LOD 100
        Cull Off // Render both sides

        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            #pragma multi_compile_fog
            
            #include "UnityCG.cginc"

            struct appdata
            {
                float4 vertex : POSITION;
                float3 normal : NORMAL;
            };

            struct v2f
            {
                UNITY_FOG_COORDS(1)
                float4 vertex : SV_POSITION;
                float3 normal : NORMAL;
                float3 worldPos : TEXCOORD1;
            };

            sampler2D _MainTex;
            float4 _FrontColor;
            float4 _BackColor;
            float4 _SpecColor;
            float _Glossiness;
            float _BumpStrength;
            float _NoiseScale;
            float _NoiseDetail;

            // Hash function for Perlin noise
            float hash(float3 p)
            {
                p = frac(p * 0.3183099 + 0.1);
                p *= 17.0;
                return frac(p.x * p.y * p.z * (p.x + p.y + p.z));
            }

            // Perlin noise function with multiple octaves
            float perlinNoise(float3 pos)
            {
                float total = 0.0;
                float frequency = _NoiseScale;
                float amplitude = 1.0;
                float maxValue = 0.0;

                for (int i = 0; i < _NoiseDetail; i++)
                {
                    float noiseVal = hash(pos * frequency);
                    total += noiseVal * amplitude;
                    maxValue += amplitude;

                    frequency *= 2.0;
                    amplitude *= 0.5;
                }

                return total / maxValue; // Normalize noise to 0-1 range
            }

            v2f vert (appdata v)
            {
                v2f o;
                o.vertex = UnityObjectToClipPos(v.vertex);
                o.normal = UnityObjectToWorldNormal(v.normal);
                o.worldPos = mul(unity_ObjectToWorld, v.vertex).xyz;
                UNITY_TRANSFER_FOG(o,o.vertex);
                return o;
            }

            fixed4 frag (v2f i, fixed facing : VFACE) : SV_Target
            {
                fixed4 col = facing > 0 ? _FrontColor : _BackColor;

                // Generate procedural roughness using Perlin noise
                float roughness = perlinNoise(i.worldPos * _NoiseScale) * _Glossiness;

                // Compute lighting
                half3 worldLight = normalize(_WorldSpaceLightPos0.xyz);
                half3 viewDir = normalize(_WorldSpaceCameraPos.xyz - i.worldPos);
                half3 halfVector = normalize(worldLight + viewDir);

                // Add procedural bump mapping with Perlin noise
                float3 noiseVec = perlinNoise(i.worldPos * _NoiseScale) - 0.5;
                float3 perturbedNormal = normalize(i.normal + _BumpStrength * noiseVec);

                // Diffuse lighting
                fixed NdotL = max(0, dot(perturbedNormal, worldLight));
                fixed3 lighting = NdotL * float3(1,1,1) + ShadeSH9(half4(i.normal, 1));

                // Blinn-Phong specular
                float NdotH = max(0, dot(perturbedNormal, halfVector));
                float spec = pow(NdotH, roughness * 128.0);
                fixed3 specular = spec * _SpecColor.rgb;

                col.rgb *= lighting;
                col.rgb += specular;

                // Apply fog
                UNITY_APPLY_FOG(i.fogCoord, col);
                return col;
            }
            ENDCG
        }
    }
}
