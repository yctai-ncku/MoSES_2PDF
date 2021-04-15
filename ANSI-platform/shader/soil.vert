#version 330 core

layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;
layout (location = 2) in vec2 aTexcoord;

out vec3 fNormal;
out vec3 fPosition;
out vec2 fTexcoord;

uniform mat4 model;
uniform mat4 projection;
uniform mat4 view;

void main()
{
	vec4 pos = projection * view * model * vec4(aPos, 1.0);
	gl_Position = pos;
	
	fPosition = pos.xyz;
	//fNormal = normalize(modelInvTransp * aNormal);
	fNormal = mat3(model) * aNormal;
	fTexcoord = aTexcoord;
}
