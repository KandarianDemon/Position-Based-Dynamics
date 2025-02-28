Shader "Custom/CheckerboardWithLighting"
{
    Properties
    {
        _Scale ("Pattern Size", Range(0,10)) = 1
        _Greyness ("Pattern Greyness", Range(0,1)) = 0.5
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" }
        LOD 100

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
                float4 position : SV_POSITION;
                float3 worldPos : TEXCOORD0;
                float3 normal : NORMAL;
            };

            float _Scale;
            float _Greyness;

            v2f vert (appdata v)
            {
                v2f o;
                o.position = UnityObjectToClipPos(v.vertex);
                o.worldPos = mul(unity_ObjectToWorld, v.vertex);
                o.normal = mul((float3x3)unity_ObjectToWorld, v.normal); // Transform normal to world space
                return o;
            }

            fixed4 frag (v2f i) : SV_Target
            {
                // Create the checkerboard pattern
                float3 adjWorldPos = floor(i.worldPos / _Scale);
                float chessboard = floor(adjWorldPos.x) + floor(adjWorldPos.z) + floor(adjWorldPos.y);
                chessboard = frac(chessboard * 0.5) + 0.2;  // Making the pattern

                // Light interaction (diffuse)
                float3 lightDir = normalize(UnityWorldSpaceLightDir(i.normal));
                float diff = max(0, dot(i.normal, lightDir));

                // Combine the pattern with the lighting and greyness
                float finalColor = chessboard - _Greyness * diff;

                // Output final color
                return float4(finalColor, finalColor, finalColor, 1);
            }
            ENDCG
        }
    }
}
