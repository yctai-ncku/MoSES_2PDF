#version 330 core

in vec3 fPosition;
in vec3 fNormal;
in vec2 fTexcoord;
in float fHeight[5];

layout (location = 0) out vec4 color;

uniform sampler2D SamplerTerrain;
uniform float SamplerCoordDivisor;

uniform vec3 light_color = vec3(1.0f, 1.0f, 1.0f);
// uniform vec3 object_color = vec3(0.75f, 0.75f, 0.75f);
uniform float shininess = 8.0f;

// uniform vec3 lightPos = vec3(150.0f, 150.0f, 150.0f);
uniform vec3 lightPos;

uniform vec3 colorMap[500];
uniform float minimum;
uniform float maximum;
uniform int nShades;
uniform float alpha;

void main()
{
	int index = 0;

	if (fHeight[0] <= minimum)
		index = 0;
	else if (fHeight[0] >= maximum)
		index = nShades;
	else {
		index = int((fHeight[0] - minimum) / (maximum - minimum) * nShades);
	}
	vec3 c = colorMap[index];
	color = vec4(c / 255, alpha);
}
