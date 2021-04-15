#version 330 core

layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;
layout (location = 2) in vec2 aTexcoord;
layout (location = 3) in float aHeight[5];
layout (location = 8) in vec3 aNvertices[4];

out vec3 fNormal;
out vec3 fPosition;
out vec2 fTexcoord;
out float fHeight[5];

uniform mat4 model;
uniform mat4 projection;
uniform mat4 view;
uniform mat3 modelInvTransp;

void main()
{
    vec3 normal;
    vec3 curVertex = aPos;

    normal = mat3(model) * aNormal;
    
    vec4 pos = projection * view * model * vec4(curVertex, 1.0f);
	gl_Position = pos;
	
	fPosition = pos.xyz;
	//fPosition.y = fHeight;
    fNormal = normal;
	fTexcoord = aTexcoord;
	for (int i = 0; i < 5; i++)
		fHeight[i] = aHeight[i];
}
