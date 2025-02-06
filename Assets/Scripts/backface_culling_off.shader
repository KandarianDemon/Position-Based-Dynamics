Shader "Unlit/backface_culling_off"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
        _FrontColor ("Front Color", Color) = (1,1,1,1)
        _BackColor ("Back Color", Color) = (1,1,1,1)
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
            // make fog work
            #pragma multi_compile_fog
            

            #include "UnityCG.cginc"

            struct appdata
            {
                float4 vertex : POSITION;
                float2 uv : TEXCOORD0;
                float3 normal: NORMAL;
            };

            struct v2f
            {
                float2 uv : TEXCOORD0;
                UNITY_FOG_COORDS(1)
                float4 vertex : SV_POSITION;
                float3 normal: NORMAL;
            };

            sampler2D _MainTex;
            float4 _MainTex_ST;
            float4 _FrontColor;
            float4 _BackColor;

            v2f vert (appdata v)
            {
                v2f o;
                o.vertex = UnityObjectToClipPos(v.vertex);
                o.uv = TRANSFORM_TEX(v.uv, _MainTex);
                o.normal = UnityObjectToWorldNormal(v.normal);
                UNITY_TRANSFER_FOG(o,o.vertex);
                return o;
            }

            fixed4 frag (v2f i, fixed facing:VFACE) : SV_Target
            {
                // sample the texture
                fixed4 col = tex2D(_MainTex, i.uv);


                // World light diffuse shading
                half3 worldLight = normalize(_WorldSpaceLightPos0.xyz);
                fixed NdotL = max(0, dot(i.normal, worldLight));
                fixed3 lighting = NdotL * float3(1,1,1) + ShadeSH9(half4(i.normal, 1));
                
                

                col *= facing >0 ? _FrontColor : _BackColor;
                col.rgb *= lighting;
                // apply fog
                UNITY_APPLY_FOG(i.fogCoord, col);
                return col;
            }
            ENDCG
        }
    }
}
