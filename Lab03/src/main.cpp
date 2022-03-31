#define _USE_MATH_DEFINES
#include <algorithm>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <exception>
#include <stdexcept>
#include <string.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include <string>
#include <vector>			//Standard template library class
#include <random>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <omp.h>
//glm
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/half_float.hpp>
//imgui
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "shaders.h"
#include "OBJ_Loader.h"
#include "TriangleBoxIntersection.hpp"

GLFWwindow* window;
//the main window size
GLint wWindow = 800;
GLint hWindow = 600;

//shader program ID
GLuint shaderProgram;
GLfloat ftime = 0.f;

class ShaderParamsC
{
public:
	GLint viewParameter;		//viewing matrix
	GLint projParameter;		//projection matrix
} params;

std::mt19937 rng;
std::uniform_real_distribution<float> noise;

const float epsilon = 3.0;

/*********************************
Mesh!!!!!!
**********************************/

struct Mesh
{
	std::vector<glm::vec3> vertices;
	std::vector<glm::vec3> normals;
	std::vector<unsigned int> indices;
};
Mesh mesh;

const std::string meshFile = "teapot.obj";
GLuint meshVao = -1;
GLuint meshVboHandle[3];

void InitializeMesh();

/*********************************
PARTICLES!!!!!!
**********************************/

GLuint particlesVao = -1;
GLuint particlesVbo = -1;

struct Particle
{
	glm::vec3 pos;
	glm::vec3 vel;
	float mass = 1.0;
	int gridId = -1;
};
std::vector<Particle> particles;
const int particleNum = 512;
float particleRadius = 0.8;

float viscosityFactor = 0.8;
glm::vec3 windDir = glm::vec3(1.0, 0.5, 0.5);
float windStrength = 0.0;
float dt = 7.0; // * 0.0001
bool enablePPCol = true;

float random(float min, float max)
{
	std::random_device rd;
	rng = std::mt19937(rd());
	noise = std::uniform_real_distribution<float>(min, max);
	return noise(rng);
}

void UpdateParticlesGPUData();
void ResetParticles();
void InitializeParticles();
void UpdateParticles();

/*********************************
UNIFORM GRIDS!!!!!!
**********************************/

int gridDim = 32;
float gridMinX = -7.0;
float gridMaxX = 7.0;
float gridMinY = -5.0;
float gridMaxY = 9.0;
float gridMinZ = -7.0;
float gridMaxZ = 7.0;
glm::vec3 gridSize;

struct Grid
{
	std::vector<int> triangles;
};
std::vector<Grid> grids;

void InitializeUniformGrids();

/*********************************
Some OpenGL-related functions
**********************************/

//called when a window is reshaped
void Reshape(GLFWwindow* window, int w, int h)
{
	glViewport(0, 0, w, h);
	//remember the settings for the camera
	wWindow = w;
	hWindow = h;
}

void DrawGui()
{
	ImGui_ImplOpenGL3_NewFrame();
	ImGui_ImplGlfw_NewFrame();
	ImGui::NewFrame();

	ImGui::Begin("G-U-I");
	ImGui::SliderFloat("Viscosity", &viscosityFactor, 0.0, 1.0);
	float pickedWindDir[3] = { windDir.x, windDir.y, windDir.z };
	ImGui::ColorEdit3("Wind direction", pickedWindDir);
	windDir.x = pickedWindDir[0];
	windDir.y = pickedWindDir[1];
	windDir.z = pickedWindDir[2];
	ImGui::SliderFloat("Wind strength", &windStrength, 0.0, 1.0);
	ImGui::SliderFloat("dt * 0.0001", &dt, 1.0, 20.0);
	ImGui::Checkbox("Enable particle-particle collision", &enablePPCol);
	ImGui::SameLine();
	if (ImGui::Button("Reset"))
		ResetParticles();

	ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
	ImGui::End();

	ImGui::Render();
	ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

//the main rendering function
void RenderObjects()
{
	const int range = 3;
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glColor3f(0, 0, 0);
	glPointSize(2);
	glLineWidth(1);
	//set the projection and view once for the scene
	glm::mat4 view = glm::mat4(1.0);
	glm::mat4 proj = glm::perspective(80.0f,//fovy
		(float)wWindow / (float)hWindow,//aspect
		0.01f, 1000.f); //near, far
	glUniformMatrix4fv(params.projParameter, 1, GL_FALSE, glm::value_ptr(proj));
	view = glm::lookAt(glm::vec3(-10.f, 0.f, -10.f),//eye
		glm::vec3(0, 0, 0),  //destination
		glm::vec3(0, 1, 0)); //up

	glUniformMatrix4fv(params.viewParameter, 1, GL_FALSE, glm::value_ptr(view));

	GLuint particlesModeLoc = glGetUniformLocation(shaderProgram, "particlesMode");
	glUniform1i(particlesModeLoc, 0);
	glBindVertexArray(meshVao);
	glDrawElements(GL_TRIANGLES, mesh.indices.size(), GL_UNSIGNED_INT, 0);

	glUniform1i(particlesModeLoc, 1);
	glBindVertexArray(particlesVao);
	glDrawArrays(GL_POINTS, 0, particleNum);
	glBindVertexArray(0);

	DrawGui();
}

void Idle(void)
{
	UpdateParticles();

	glClearColor(0.1, 0.1, 0.1, 1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	ftime += dt;
	glUseProgram(shaderProgram);
	RenderObjects();
	glfwSwapBuffers(window);
}

void Display(void)
{

}

//keyboard callback
void Kbd(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (action == GLFW_PRESS)
	{
		switch (key)
		{
		case GLFW_KEY_ESCAPE:
			glfwSetWindowShouldClose(window, GLFW_TRUE);
			break;
		case 'r':
		case 'R': ResetParticles(); break;
		}
	}
}

void InitializeProgram(GLuint* program)
{
	std::vector<GLuint> shaderList;

	//load and compile shaders 	
	shaderList.push_back(CreateShader(GL_VERTEX_SHADER, LoadShader("shaders/phong.vert")));
	shaderList.push_back(CreateShader(GL_FRAGMENT_SHADER, LoadShader("shaders/phong.frag")));

	//create the shader program and attach the shaders to it
	*program = CreateProgram(shaderList);

	//delete shaders (they are on the GPU now)
	std::for_each(shaderList.begin(), shaderList.end(), glDeleteShader);

	params.viewParameter = glGetUniformLocation(*program, "view");
	params.projParameter = glGetUniformLocation(*program, "proj");

	glBindVertexArray(0);
}

int InitializeGL(GLFWwindow*& window)
{
	if (!glfwInit())
		return -1;
	window = glfwCreateWindow(wWindow, hWindow, "Lab3", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		return -1;
	}
	glfwSetKeyCallback(window, Kbd);
	glfwSetFramebufferSizeCallback(window, Reshape);
	glfwMakeContextCurrent(window);
	glewInit();
	InitializeProgram(&shaderProgram);
	glEnable(GL_DEPTH_TEST);

	return 0;
}

void InitializeMesh()
{
	std::cout << "Loading mesh: " << meshFile << " ..." << std::endl;

	objl::Loader loader;
	if (loader.LoadFile(meshFile))
	{
		for (const auto& v : loader.LoadedVertices)
		{
			mesh.vertices.push_back(glm::vec3(v.Position.X / 13.0, v.Position.Y / 13.0 - 5.01, v.Position.Z / 13.0));
			mesh.normals.push_back(glm::vec3(v.Normal.X, v.Normal.Y, v.Normal.Z));
		}
		mesh.indices = loader.LoadedIndices;
	}

	int numVerts = mesh.vertices.size();
	int numNorms = mesh.normals.size();
	glGenVertexArrays(1, &meshVao);
	glBindVertexArray(meshVao);
	glGenBuffers(3, meshVboHandle);
	glBindBuffer(GL_ARRAY_BUFFER, meshVboHandle[0]);
	glBufferData(GL_ARRAY_BUFFER, numVerts * 3 * sizeof(float), mesh.vertices.data(), GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, meshVboHandle[1]);
	glBufferData(GL_ARRAY_BUFFER, numNorms * 3 * sizeof(float), mesh.normals.data(), GL_STATIC_DRAW);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(1);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, meshVboHandle[2]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * mesh.indices.size(), mesh.indices.data(), GL_STATIC_DRAW);
	glBindVertexArray(0);

	std::cout << "Complete loading mesh." << std::endl;
}

void UpdateParticlesGPUData()
{
	std::vector<glm::vec3> vertices;
	for (const auto& p : particles)
		vertices.push_back(p.pos);
	glBindVertexArray(particlesVao);
	glBindBuffer(GL_ARRAY_BUFFER, particlesVbo);
	glBufferData(GL_ARRAY_BUFFER, 3 * sizeof(float) * particleNum, vertices.data(), GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glEnableVertexAttribArray(0);
	glBindVertexArray(0);
}

void ResetParticles()
{
	particles.swap(std::vector<Particle>());
	for (int i = 0; i < particleNum; i++)
	{
		Particle p;
		p.pos = glm::vec3(random(-3.0, 3.0), 5.0f, random(-3.0, 3.0));
		p.vel = glm::vec3(0.0f);
		p.mass = random(0.8, 1.0);
		particles.push_back(p);
	}
	UpdateParticlesGPUData();
}

void InitializeParticles()
{
	glGenVertexArrays(1, &particlesVao);
	glGenBuffers(1, &particlesVbo);
	ResetParticles();
}

bool TriangleParticleIntersect(glm::vec3 p, int triangleIdx)
{
	glm::vec3 n = glm::normalize(mesh.normals[mesh.indices[triangleIdx]]);
	float dist = glm::dot(mesh.vertices[mesh.indices[triangleIdx]] - p, n);
	if (dist > particleRadius)
		return false;
	glm::vec3 q = p - n * dist;
	
	float r = 0.0;
	glm::vec3 vecOld, vecNew;
	for (int j = 0; j < 3; j++) {
		if (j == 0)
			vecOld = mesh.vertices[mesh.indices[triangleIdx + j]] - p;
		if (j == 2)
			vecNew = mesh.vertices[mesh.indices[triangleIdx]] - p;
		else
			vecNew = mesh.vertices[mesh.indices[triangleIdx + j + 1]] - p;

		r += glm::acos(glm::dot(glm::normalize(vecOld), glm::normalize(vecNew)));
		vecOld = vecNew;
	}
	// If the sum of all the angles is approximately equal to 360 degrees
	// then the intersection is in the triangle.
	if (abs(r - 2.0f * (float)std::_Pi) < epsilon)
		return true;

	return false;
}

void UpdateParticles()
{
	int numThreads = omp_get_max_threads();
	if (numThreads > 0)
		numThreads--;
	#pragma omp parallel for num_threads(numThreads)
	for (int i = 0; i < particleNum; i++)
	{
		glm::vec3 gravity = 9.82f * glm::vec3(0.0, -1.0, 0.0);
		glm::vec3 viscosity = -particles[i].vel * viscosityFactor / particles[i].mass;
		glm::vec3 wind = windStrength * glm::normalize(windDir * 2.0f - 1.0f);

		glm::vec3 bounceForce;
		// Floor
		if (particles[i].pos.y <= gridMinY + 0.1)
		{
			particles[i].pos.y = gridMinY + 0.1;
			particles[i].vel = glm::reflect(particles[i].vel, glm::vec3(0.0, 1.0, 0.0));
		}
		// Ceiling
		if (particles[i].pos.y >= gridMaxY - 0.1)
		{
			particles[i].pos.y = gridMaxY - 0.1;
			particles[i].vel = glm::reflect(particles[i].vel, glm::vec3(0.0, -1.0, 0.0));
		}
		// Left wall
		if (particles[i].pos.x <= gridMinX + 0.1)
		{
			particles[i].pos.x = gridMinX + 0.1;
			particles[i].vel = glm::reflect(particles[i].vel, glm::vec3(1.0, 0.0, 0.0));
		}
		// Right wall
		if (particles[i].pos.x >= gridMaxX - 0.1)
		{
			particles[i].pos.x = gridMaxX - 0.1;
			particles[i].vel = glm::reflect(particles[i].vel, glm::vec3(-1.0, 0.0, 0.0));
		}
		// Front wall
		if (particles[i].pos.z <= gridMinZ + 0.1)
		{
			particles[i].pos.z = gridMinZ + 0.1;
			particles[i].vel = glm::reflect(particles[i].vel, glm::vec3(0.0, 0.0, 1.0));
		}
		// Back wall
		if (particles[i].pos.z >= gridMaxZ - 0.1)
		{
			particles[i].pos.z = gridMaxZ - 0.1;
			particles[i].vel = glm::reflect(particles[i].vel, glm::vec3(0.0, 0.0, -1.0));
		}

		// Uniform grid location
		glm::ivec3 gridPosLoc;
		gridPosLoc.x = (particles[i].pos.x - gridMinX) / gridSize.x;
		gridPosLoc.y = (particles[i].pos.y - gridMinY) / gridSize.y;
		gridPosLoc.z = (particles[i].pos.z - gridMinZ) / gridSize.z;
		int gridId = gridPosLoc.z * (gridDim * gridDim) + gridPosLoc.y * gridDim + gridPosLoc.x;
		particles[i].gridId = gridId;

		// Collision with mesh
		for (int j = 0; j < grids[gridId].triangles.size(); j++)
		{
			if (TriangleParticleIntersect(particles[i].pos, grids[gridId].triangles[j]))
			{
				glm::vec3 n = mesh.normals[mesh.indices[grids[gridId].triangles[j]]];
				n = glm::normalize(n);
				particles[i].vel = glm::reflect(particles[i].vel, n) * 0.999f;
				break;
			}
		}

		// Collision with particles
		if (enablePPCol)
		{
			for (int j = 0; j < particleNum; j++)
			{
				if (i != j && particles[j].gridId == gridId)
				{
					float d = glm::distance(particles[i].pos, particles[j].pos);
					if (d < particleRadius * 2.0)
					{
						glm::vec3 n = particles[i].pos - particles[j].pos;
						if (glm::length(n) > 0.001f)
						{
							bounceForce = glm::length(particles[i].vel) * 3.0f * glm::normalize(n);
						}
					}
				}
			}
		}

		glm::vec3 accumulator = gravity + viscosity + wind + bounceForce;
		particles[i].vel += dt * 0.0001f * accumulator;
		particles[i].pos += dt * 0.0001f * particles[i].vel;
	}
	UpdateParticlesGPUData();
}

void InitializeUniformGrids()
{
	std::cout << "Building uniform grids..." << std::endl;

	float gridSizeX = fabs(gridMaxX - gridMinX) / (float)gridDim;
	float gridSizeY = fabs(gridMaxY - gridMinY) / (float)gridDim;
	float gridSizeZ = fabs(gridMaxZ - gridMinZ) / (float)gridDim;
	gridSize = glm::vec3(gridSizeX, gridSizeY, gridSizeZ);

	for (int i = 0; i < pow(gridDim, 3); i++)
	{
		Grid g;
		grids.push_back(g);
	}
	glm::vec3 boxHalfSize(gridSizeX * 0.5, gridSizeY * 0.5, gridSizeZ * 0.5);

	int numThreads = omp_get_max_threads();
	if (numThreads > 0)
		numThreads--;
	#pragma omp parallel for num_threads(numThreads)
	for (int i = 0; i < grids.size(); i++)
	{
		glm::ivec3 posIdx;
		posIdx.x = i % gridDim;
		posIdx.y = (i / gridDim) % gridDim;
		posIdx.z = (i / (gridDim * gridDim)) % gridDim;
		glm::vec3 boxCenter;
		boxCenter.x = gridMinX + posIdx.x * gridSizeX + gridSizeX * 0.5;
		boxCenter.y = gridMinY + posIdx.y * gridSizeY + gridSizeY * 0.5;
		boxCenter.z = gridMinZ + posIdx.z * gridSizeZ + gridSizeZ * 0.5;

		for (int j = 0; j < mesh.indices.size(); j += 3)
		{
			glm::vec3 v1 = mesh.vertices[mesh.indices[j]];
			glm::vec3 v2 = mesh.vertices[mesh.indices[j + 1]];
			glm::vec3 v3 = mesh.vertices[mesh.indices[j + 2]];
			if (triBoxOverlap(boxCenter, boxHalfSize, v1, v2, v3))
				grids[i].triangles.push_back(j);
		}
	}

	std::cout << "Complete building uniform grids." << std::endl;
}

int main(int argc, char** argv)
{
	int initRes = InitializeGL(window);
	if (initRes)
		return initRes;

	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL3_Init("#version 150");

	InitializeMesh();
	InitializeParticles();
	InitializeUniformGrids();

	while (!glfwWindowShouldClose(window))
	{
		Idle();
		Display();
		glfwPollEvents();
	}

	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();

	glfwTerminate();
	return 0;
}
