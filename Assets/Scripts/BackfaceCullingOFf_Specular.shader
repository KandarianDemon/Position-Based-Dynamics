Shader "Unlit/BackfaceCullingOff_Specular"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
        _FrontColor ("Front Color", Color) = (1,1,1,1)
        _BackColor ("Back Color", Color) = (1,1,1,1)
        _SpecColor ("Specular Color", Color) = (1,1,1,1)
        _Glossiness ("Glossiness", Range(0, 1)) = 0.5
        _BumpStrength ("Bump Strength", Range(0, 1)) = 0.2
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" }
        LOD 100
        Cull Off

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
                float2 uv : TEXCOORD0;
                float3 normal : NORMAL;
            };

            struct v2f
            {
                float2 uv : TEXCOORD0;
                UNITY_FOG_COORDS(1)
                float4 vertex : SV_POSITION;
                float3 normal : NORMAL;
                float3 worldPos : TEXCOORD1;
            };

            sampler2D _MainTex;
            float4 _MainTex_ST;
            float4 _FrontColor;
            float4 _BackColor;
            float4 _SpecColor;
            float _Glossiness;
            float _BumpStrength;

            // Simple procedural noise (cheap alternative to texture)
            float noise(float3 pos)
            {
                return frac(sin(dot(pos, float3(12.9898, 78.233, 45.543))) * 43758.5453);
            }

            v2f vert (appdata v)
            {
                v2f o;
                o.vertex = UnityObjectToClipPos(v.vertex);
                o.uv = TRANSFORM_TEX(v.uv, _MainTex);
                o.normal = UnityObjectToWorldNormal(v.normal);
                o.worldPos = mul(unity_ObjectToWorld, v.vertex).xyz;
                UNITY_TRANSFER_FOG(o, o.vertex);
                return o;
            }

            fixed4 frag (v2f i, fixed facing : VFACE) : SV_Target
            {
                fixed4 col = tex2D(_MainTex, i.uv);
                col *= facing > 0 ? _FrontColor : _BackColor;

                // Compute lighting
                half3 worldLight = normalize(_WorldSpaceLightPos0.xyz);
                half3 viewDir = normalize(_WorldSpaceCameraPos.xyz - i.worldPos);
                half3 halfVector = normalize(worldLight + viewDir);

                // Add procedural normal perturbation
                float3 noiseVal = noise(i.worldPos * 20.0);
                float3 perturbedNormal = normalize(i.normal + _BumpStrength * (noiseVal - 0.5));

                // Diffuse lighting (Lambertian)
                fixed NdotL = max(0, dot(perturbedNormal, worldLight));
                fixed3 lighting = NdotL * float3(1,1,1);

                // Blinn-Phong specular
                float NdotH = max(0, dot(perturbedNormal, halfVector));
                float spec = pow(NdotH, _Glossiness * 128.0); // Glossiness scaled to specular exponent
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
