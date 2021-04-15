#version 330 core

in vec3 fPosition;
in vec3 fNormal;
in vec2 fTexcoord;

layout (location = 0) out vec4 color;

uniform sampler2D SamplerTerrain;
uniform float SamplerCoordDivisor;
uniform sampler2D SamplerSoil;

uniform vec3 light_color = vec3(1.0f, 1.0f, 1.0f);
uniform vec3 object_color = vec3(0.75f, 0.75f, 0.75f);
uniform float shininess = 8.0f;

// uniform vec3 lightPos = vec3(150.0f, 150.0f, 150.0f);
uniform vec3 lightPos;
// uniform vec3 viewPos;
uniform int grayFlag;


void main()
{
	// vec2 tex_coords = vec2(fPosition.x, fPosition.z) / SamplerCoordDivisor;
    vec3 tex, tex2;
    float alpha = 0.5;
    if (grayFlag == 1) {
        tex = vec3(0.75f, 0.75f, 0.75f);
    }
    else {
//            tex = vec3(texture(SamplerSoil, fTexcoord));
        tex = vec3(texture(SamplerTerrain, fTexcoord * SamplerCoordDivisor));
        //tex = vec3(texture(SamplerSoil, fTexcoord));
//        tex = tex * alpha + tex2 * (1 - alpha);
    }

    //vec3 tex = vec3(texture(SamplerTerrain, fTexcoord * SamplerCoordDivisor));
    
    vec3 lightDir = normalize(lightPos - fPosition);
	vec3 normal = normalize(fNormal);
	float diffuse = max(0.0f, dot(normal, lightDir));

	vec3 viewDir = normalize(lightPos - fPosition);
	
    vec3 r = reflect(-lightDir, normal);
    float spec = pow(max(0.0f, dot(viewDir, r)), shininess);
	
	float color_ambient = 0.05;
	float color_diffuse = diffuse;
	float color_specular = 0.3 * spec;

	vec3 result = (color_ambient + color_diffuse + color_specular) * light_color * tex;
//    vec3 result = tex;
  
    //color = vec4(1.0);
    color = vec4(result, 1.0f);
}
