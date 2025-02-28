Shader "Custom/WormParallaxShader"
{
    Properties
    {
        _MainTex ("Main Texture", 2D) = "white" {}
        _NeuronsTex ("Neurons Depth Texture", 2D) = "black" {}
        _Transparency ("Transparency", Range(0,1)) = 0.5
        _ParallaxDepth ("Parallax Depth", Range(0.001, 0.1)) = 0.02
        _ScrollSpeed ("Scroll Speed", Range(0,2)) = 0.5
    }
    SubShader
    {
        Tags {"Queue"="Transparent" "RenderType"="Transparent"}
        Blend SrcAlpha OneMinusSrcAlpha
        ZWrite Off
        Cull Back
        
        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            #pragma multi_compile_fog

            #include "UnityCG.cginc"

            struct appdata_t
            {
                float4 vertex : POSITION;
                float2 uv : TEXCOORD0;
                float3 normal : NORMAL;
            };

            struct v2f
            {
                float4 pos : SV_POSITION;
                float2 uv : TEXCOORD0;
                float3 viewDir : TEXCOORD1;
            };

            sampler2D _MainTex;
            sampler2D _NeuronsTex;
            float _Transparency;
            float _ParallaxDepth;
            float _ScrollSpeed;

            v2f vert (appdata_t v)
            {
                v2f o;
                o.pos = UnityObjectToClipPos(v.vertex);
                o.uv = v.uv;
                
                // Compute view direction for parallax effect
                float3 viewDir = normalize(ObjSpaceViewDir(v.vertex));
                o.viewDir = viewDir;
                
                return o;
            }

            fixed4 frag (v2f i) : SV_Target
            {
                float2 uv = i.uv;
                
                // Simulate neuron depth using parallax mapping
                float2 parallaxOffset = (i.viewDir.xy * _ParallaxDepth);
                float timeOffset = _Time.y * _ScrollSpeed;

                float2 neuronUV = uv + parallaxOffset;
                neuronUV.y += timeOffset; // Scrolling neurons effect

                // Sample neuron texture
                fixed4 neuronCol = tex2D(_NeuronsTex, neuronUV);

                // Outer transparency effect
                fixed4 col = tex2D(_MainTex, uv);
                col.a = _Transparency * (1.0 - neuronCol.r); // Adjust transparency with neuron visibility

                return col;
            }
            ENDCG
        }
    }
}
