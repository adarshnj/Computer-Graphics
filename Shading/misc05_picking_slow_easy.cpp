// Include standard headers
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <array>
#include <stack>   
#include <sstream>
#include <math.h>
// Include GLEW
#include <GL/glew.h>
// Include GLFW
#include <glfw3.h>
// Include GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>
using namespace glm;
// Include AntTweakBar
#include <AntTweakBar.h>

#include <common/shader.hpp>
#include <common/controls.hpp>
#include <common/objloader.hpp>
#include <common/vboindexer.hpp>

const int window_width = 1024, window_height = 768;

typedef struct Vertex {
	float Position[4];
	float Color[4];
	float Normal[3];
	void SetPosition(float *coords) {
		Position[0] = coords[0];
		Position[1] = coords[1];
		Position[2] = coords[2];
		Position[3] = 1.0;
	}
	void SetColor(float *color) {
		Color[0] = color[0];
		Color[1] = color[1];
		Color[2] = color[2];
		Color[3] = color[3];
	}
	void SetNormal(float *coords) {
		Normal[0] = coords[0];
		Normal[1] = coords[1];
		Normal[2] = coords[2];
	}
};

	// function prototypes
	int initWindow(void);
	void initOpenGL(void);
	void loadObject(char*, glm::vec4, Vertex * &, GLushort* &, int);
	void createVAOs(Vertex[], GLushort[], int);
	void createObjects(void);
	void pickObject(void);
	void renderScene(void);
	void cleanup(void);
	static void keyCallback(GLFWwindow*, int, int, int, int);
	static void mouseCallback(GLFWwindow*, int, int, int);
	static void scrollCallback(GLFWwindow*, double, double);

	void grid(void);

	// GLOBAL VARIABLES
	GLFWwindow* window;

	glm::mat4 gProjectionMatrix;
	glm::mat4 gViewMatrix;

	GLuint gPickedIndex = -1;
	std::string gMessage;

	GLuint programID;
	GLuint pickingProgramID;

	const GLuint NumObjects = 10;	// ATTN: THIS NEEDS TO CHANGE AS YOU ADD NEW OBJECTS
	GLuint VertexArrayId[NumObjects] = { 0 };
	GLuint VertexBufferId[NumObjects] = { 0 };
	GLuint IndexBufferId[NumObjects] = { 0 };

	size_t NumIndices[NumObjects] = { 0 };
	size_t VertexBufferSize[NumObjects] = { 0 };
	size_t IndexBufferSize[NumObjects] = { 0 };

	GLuint MatrixID;
	GLuint ModelMatrixID;
	GLuint ViewMatrixID;
	GLuint ProjMatrixID;
	GLuint PickingMatrixID;
	GLuint pickingColorID;
	GLuint LightID;
	GLuint LightID1;

	float xtran = 0.0f;
	float ytran = 0.0f;
	float ztran = 0.0f;
	float grav = 0.0f;
	float horz = 0.0f;
	float ang = 0.0f;

	GLint gX = 0.0;
	GLint gZ = 0.0;

	// animation control
	bool animation = false;
	GLfloat phi = 0.0;

	// Objects

	bool camera = false;
	bool base = false;
	bool pen = false;
	bool top = false;
	bool arm1 = false;
	bool arm2 = false;
	bool jump = false;

	float theta = 65;		// Angle from Y axis
	float Phi = 45;			// Angle from Z axis
	float radius = 7.0f;
	float cX = 10.0f;
	float cY = 10.0f;
	float cZ = 10.0f;

	int32 ud = 2;
	int32 rl = 0;
	float y = 1.1f;
	float t = 0.0f;

	float posud  = 0.0f;
	float posrl = 0.0f;
	float rotT = 0.0f;
	float rotA1 = 0.0f;
	float rotA2 = 0.0f;
	float penrud = 0.0f;
	float penrrl = 0.0f;
	float penrev = 0.0f;
	float dropfrom = 0.0f;


	glm::vec4 baseColor = { 1.0, 0.0, 0.0, 1.0 };
	glm::vec4 topColor = { 1.0, 1.0, 0.0, 1.0 };
	glm::vec4 arm1Color = { 0.0, 1.0, 0.0, 1.0 };
	glm::vec4 jointColor = { 0.0, 1.0, 1.0, 1.0 };
	glm::vec4 arm2Color = { 0.0, 0.0, 1.0, 1.0 };
	glm::vec4 penColor = { 1.0, 0.0, 1.0, 1.0 };
	glm::vec4 buttonColor = { 1.0, 1.0, 1.0, 1.0 };
	glm::vec4 ballColor = { 1.0, 1.0, 1.0, 1.0 };	


	glm::vec3 trans = glm::vec3(0.0f, 0.0f, 0.0f);


	Vertex gridVert[120];
	unsigned short gridIndices[120];
	int gridIndex;

	Vertex* Verts, *Verts1, *Verts2, *Verts3, *Verts4, *Verts5, *Verts6, *Verts7;
	GLushort* Idcs, *Idcs1, *Idcs2, *Idcs3, *Idcs4, *Idcs5, *Idcs6, *Idcs7;


void loadObject(char* file, glm::vec4 color, Vertex * &out_Vertices, GLushort* &out_Indices, int ObjectId)	 
{
	// Read our .obj file
	std::vector<glm::vec3> vertices;
	std::vector<glm::vec2> uvs;
	std::vector<glm::vec3> normals;
	bool res = loadOBJ(file, vertices, normals);

	std::vector<GLushort> indices;
	std::vector<glm::vec3> indexed_vertices;
	std::vector<glm::vec2> indexed_uvs;
	std::vector<glm::vec3> indexed_normals;
	indexVBO(vertices, normals, indices, indexed_vertices, indexed_normals);

	const size_t vertCount = indexed_vertices.size();
	const size_t idxCount = indices.size();

	// populate output arrays
	out_Vertices = new Vertex[vertCount];
	for (int i = 0; i < vertCount; i++) {
		out_Vertices[i].SetPosition(&indexed_vertices[i].x);
		out_Vertices[i].SetNormal(&indexed_normals[i].x);
		out_Vertices[i].SetColor(&color[0]);
	}
	out_Indices = new GLushort[idxCount];
	for (int i = 0; i < idxCount; i++) {
		out_Indices[i] = indices[i];
	}

	// set global variables!!
	NumIndices[ObjectId] = idxCount;
	VertexBufferSize[ObjectId] = sizeof(out_Vertices[0]) * vertCount;
	IndexBufferSize[ObjectId] = sizeof(GLushort) * idxCount;
}

void grid(void)
{
	gridIndex = 0;
	float color[] = { 1.0f,1.0f,1.0f,1.0f };
	for (float x = -5.0f; x < 6.0f; x++)
	{
		gridVert[gridIndex].Position[0] = x;
		gridVert[gridIndex].Position[1] = 0.0f;
		gridVert[gridIndex].Position[2] = -5.0f;
		gridVert[gridIndex].Position[3] = 1.0f;

		gridIndices[gridIndex] = gridIndex;

		gridVert[gridIndex++].SetColor(color);

		gridVert[gridIndex].Position[0] = x;
		gridVert[gridIndex].Position[1] = 0.0f;
		gridVert[gridIndex].Position[2] = 5.0f;
		gridVert[gridIndex].Position[3] = 1.0f;

		gridIndices[gridIndex] = gridIndex;

		gridVert[gridIndex++].SetColor(color);
	}

	for (float x = -5.0f; x < 6.0f; x++)
	{
		gridVert[gridIndex].Position[0] = -5.0f;
		gridVert[gridIndex].Position[1] = 0.0f;
		gridVert[gridIndex].Position[2] = x;
		gridVert[gridIndex].Position[3] = 1.0f;

		gridIndices[gridIndex] = gridIndex;

		gridVert[gridIndex++].SetColor(color);

		gridVert[gridIndex].Position[0] = 5.0f;
		gridVert[gridIndex].Position[1] = 0.0f;
		gridVert[gridIndex].Position[2] = x;
		gridVert[gridIndex].Position[3] = 1.0f;

		gridIndices[gridIndex] = gridIndex;

		gridVert[gridIndex++].SetColor(color);
	}

	VertexBufferSize[1] = sizeof(gridVert);
	createVAOs(gridVert, NULL, 1);
}


void createObjects(void)
{	
	//-- COORDINATE AXES --//
	Vertex CoordVerts[] =
	{
		{ { 0.0, 0.0, 0.0, 1.0 },{ 1.0, 0.0, 0.0, 1.0 },{ 0.0, 0.0, 1.0 } },
		{ { 5.0, 0.0, 0.0, 1.0 },{ 1.0, 0.0, 0.0, 1.0 },{ 0.0, 0.0, 1.0 } },
		{ { 0.0, 0.0, 0.0, 1.0 },{ 0.0, 1.0, 0.0, 1.0 },{ 0.0, 0.0, 1.0 } },
		{ { 0.0, 5.0, 0.0, 1.0 },{ 0.0, 1.0, 0.0, 1.0 },{ 0.0, 0.0, 1.0 } },
		{ { 0.0, 0.0, 0.0, 1.0 },{ 0.0, 0.0, 1.0, 1.0 },{ 0.0, 0.0, 1.0 } },
		{ { 0.0, 0.0, 5.0, 1.0 },{ 0.0, 0.0, 1.0, 1.0 },{ 0.0, 0.0, 1.0 } },
	};

	VertexBufferSize[0] = sizeof(CoordVerts);	// ATTN: this needs to be done for each hand-made object with the ObjectID (subscript)
	createVAOs(CoordVerts, NULL, 0);

	//-- GRID --//

	// ATTN: create your grid vertices here!

	grid();

	//-- .OBJs --//

	// ATTN: load your models here
	
	loadObject("base.obj", baseColor, Verts, Idcs, 2);
	createVAOs(Verts, Idcs, 2);

	loadObject("top.obj", topColor, Verts1, Idcs1, 3);
	createVAOs(Verts1, Idcs1, 3);
	
	loadObject("arm1.obj", arm1Color, Verts2, Idcs2, 4);
	createVAOs(Verts2, Idcs2, 4);
	
	loadObject("joint.obj", jointColor, Verts3, Idcs3, 5);
	createVAOs(Verts3, Idcs3, 5);

	loadObject("arm2.obj", arm2Color, Verts4, Idcs4, 6);
	createVAOs(Verts4, Idcs4, 6);

	loadObject("pen.obj", penColor, Verts5, Idcs5, 7);
	createVAOs(Verts5, Idcs5, 7);

	loadObject("button.obj", buttonColor, Verts6, Idcs6, 8);
	createVAOs(Verts6, Idcs6, 8);	

	loadObject("ball.obj", ballColor, Verts7, Idcs7, 9);
	createVAOs(Verts7, Idcs7, 9);

}

void renderScene(void)
{
	xtran = 0.0f;
	ytran = 0.0f;
	ztran = 0.0f;
	dropfrom = 0.0f;
	ang = 0.0f;

	if (theta == 360)
		theta = 0;
	if (Phi == 360)
		Phi = 0;
	float thetaRad = radians(theta);
	float PhiRad = radians(Phi);

	cX = radius * sinf(thetaRad) * cos(PhiRad);
	cZ = radius * sinf(thetaRad) * sin(PhiRad);
	cY = radius * cosf(thetaRad);

	gViewMatrix = glm::lookAt(glm::vec3(cX, cY, cZ),	// eye
		glm::vec3(0.0, 0.0, 0.0),	// center
		glm::vec3(0.0, 1.0, 0.0));	// up

	const GLfloat col[] = { 1.0f,1.0f,1.0f,1.0f };

	//ATTN: DRAW YOUR SCENE HERE. MODIFY/ADAPT WHERE NECESSARY!


	// Dark blue background
	glClearColor(0.0f, 0.0f, 0.2f, 0.0f);
	// Re-clear the screen for real rendering
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glUseProgram(programID);
	{

		glm::vec3 lightPos = glm::vec3((cX - 0.5f), cY, cZ);
		glm::vec3 lightPos1 = glm::vec3((cX + 0.5f), cY, cZ);
		glm::mat4x4 ModelMatrix = glm::mat4(1.0);

		//glUniform3f(LightID, 3,3,3);
		glUniform3f(LightID, lightPos.x, lightPos.y, lightPos.z);
		//glUniform3f(LightID1, 2, 4, 4);
		glUniform3f(LightID1, lightPos1.x, lightPos1.y, lightPos1.z);


		glUniformMatrix4fv(ViewMatrixID, 1, GL_FALSE, &gViewMatrix[0][0]);
		glUniformMatrix4fv(ProjMatrixID, 1, GL_FALSE, &gProjectionMatrix[0][0]);
		glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);

		glBindVertexArray(VertexArrayId[0]);	// draw axes
		glDrawArrays(GL_LINES, 0, 6);

		glBindVertexArray(VertexArrayId[1]);	// draw Grid
		glDrawArrays(GL_LINES, 0, gridIndex);

		//Move all Objects - moving base
		{
			//glm::vec3 trans = glm::vec3(0.0f, 0.0f, 0.0f);
			trans[rl] = posrl;
			xtran += posrl;
			trans[ud] = posud;
			ztran += posud;

			glm::mat4x4 TranslationMatrix = glm::translate(glm::mat4x4(1.0f), trans);
			ModelMatrix = ModelMatrix * TranslationMatrix;

			glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
		}


		glBindVertexArray(VertexArrayId[2]);	// draw Base
		glDrawElements(GL_TRIANGLES, NumIndices[2], GL_UNSIGNED_SHORT, (void*)0);

		// Rotate Top along Y axis
		{
			glm::mat4x4 RotationMatrix = glm::rotate(glm::mat4x4(1.0f), radians(rotT), glm::vec3(0.0, 1.0, 0.0));
			ModelMatrix = ModelMatrix * RotationMatrix;
			glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
		}

		glBindVertexArray(VertexArrayId[3]);	// draw Top
		glDrawElements(GL_TRIANGLES, NumIndices[3], GL_UNSIGNED_SHORT, (void*)0);


		//Initial positioning
		{
			glm::mat4x4 TranslationMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(0.0, 1.77, 0.0));
			ytran += 1.77;
			glm::mat4 RotationMatrix = glm::rotate(glm::mat4x4(1.0f), radians(270.0f), glm::vec3(0.0, 1.0, 0.0));
			ModelMatrix = ModelMatrix * RotationMatrix * TranslationMatrix;
			glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
		}
		// Rotate Arm1 whlie one end is fixed to the top

		{
			glm::mat4x4 TranslationMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(-xtran, -ytran, -ztran));
			glm::mat4 RotationMatrix = glm::rotate(glm::mat4x4(1.0f), radians(rotA1), glm::vec3(1.0, 0.0, 0.0));
			ModelMatrix = ModelMatrix * RotationMatrix * TranslationMatrix;
			ang += rotA1;

			dropfrom += (1.77 * sin(rotA1));
			TranslationMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(xtran, ytran, ztran));
			ModelMatrix = ModelMatrix * TranslationMatrix;
			glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
		}

		glBindVertexArray(VertexArrayId[4]);	// draw Arm1
		glDrawElements(GL_TRIANGLES, NumIndices[4], GL_UNSIGNED_SHORT, (void*)0);

		{
			glm::mat4x4 TranslationMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(0.0, 0.0, -1.75));
			ztran += -1.75f;
			ModelMatrix = ModelMatrix * TranslationMatrix;
			glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
		}


		glBindVertexArray(VertexArrayId[5]);	// draw Joint
		glDrawElements(GL_TRIANGLES, NumIndices[5], GL_UNSIGNED_SHORT, (void*)0);

		//Rotate arm2 while one end is fixed to joint

		{
			glm::mat4x4 RotationMatrix = glm::rotate(glm::mat4x4(1.0f), radians(180.0f), glm::vec3(1.0, 0.0, 0.0));
			ModelMatrix = ModelMatrix * RotationMatrix;
			glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
		}

		{
			glm::mat4x4 TranslationMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(-xtran, -ytran, -ztran));
			glm::mat4 RotationMatrix = glm::rotate(glm::mat4x4(1.0f), radians(rotA2), glm::vec3(1.0, 0.0, 0.0));
			ModelMatrix = ModelMatrix * RotationMatrix * TranslationMatrix;
			ang += rotA2;
			dropfrom += (1.0 * sin(rotA2));
			TranslationMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(xtran, ytran, ztran));
			ModelMatrix = ModelMatrix * TranslationMatrix;
			glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
		}

		glBindVertexArray(VertexArrayId[6]);	// draw Arm2
		glDrawElements(GL_TRIANGLES, NumIndices[6], GL_UNSIGNED_SHORT, (void*)0);

		//Rotate Pen while one end is fixed to arm2		
		{
			glm::mat4x4 RotationMatrix = glm::rotate(glm::mat4x4(1.0f), radians(90.0f), glm::vec3(1.0, 0.0, 0.0));
			glm::mat4x4 TranslationMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(0.0, 1.0, 0.0));
			ytran -= 1.25;
			ModelMatrix = ModelMatrix * TranslationMatrix * RotationMatrix;
			glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
		}

		{
			glm::mat4x4 TranslationMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(-xtran, -ytran, -ztran));
			glm::mat4 RotationMatrix = glm::rotate(glm::mat4x4(1.0f), radians(penrud), glm::vec3(1.0, 0.0, 0.0));
			ModelMatrix = ModelMatrix * RotationMatrix * TranslationMatrix;
			ang += penrud;
			dropfrom += (1.0 * sin(penrud));
			TranslationMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(xtran, ytran, ztran));
			ModelMatrix = ModelMatrix * TranslationMatrix;
			glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
		}

		{
			glm::mat4x4 TranslationMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(-xtran, -ytran, -ztran));
			glm::mat4 RotationMatrix = glm::rotate(glm::mat4x4(1.0f), radians(penrrl), glm::vec3(0.0, 0.0, 1.0));
			ModelMatrix = ModelMatrix * RotationMatrix * TranslationMatrix;
			TranslationMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(xtran, ytran, ztran));
			ModelMatrix = ModelMatrix * TranslationMatrix;
			glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
		}

		{
			glm::mat4x4 TranslationMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(-xtran, -ytran, -ztran));
			glm::mat4 RotationMatrix = glm::rotate(glm::mat4x4(1.0f), radians(penrev), glm::vec3(0.0, 1.0, 0.0));
			ModelMatrix = ModelMatrix * RotationMatrix * TranslationMatrix;
			TranslationMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(xtran, ytran, ztran));
			ModelMatrix = ModelMatrix * TranslationMatrix;
			glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
		}

		glBindVertexArray(VertexArrayId[7]);	// draw Pen
		glDrawElements(GL_TRIANGLES, NumIndices[7], GL_UNSIGNED_SHORT, (void*)0);

		glBindVertexArray(VertexArrayId[8]);	// draw Button
		glDrawElements(GL_TRIANGLES, NumIndices[8], GL_UNSIGNED_SHORT, (void*)0);

		if (jump == true)
		{
			{
				glm::mat4x4 TranslationMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(-xtran, -ytran, -ztran));
				ModelMatrix = ModelMatrix * TranslationMatrix;
			}

			horz = y + (t * 1 * cos(penrud));
			grav = ((1 * t*sin(penrud)) - (0.5*t*t*9.8f));
			//printf("%f %f %f\n", xtran, ytran, ztran);

			//printf("%f\n", grav);
			glm::mat4x4 TranslationMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(0.0, horz, grav)); //y+(1*t)

			{
				glm::mat4x4 TranslationMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(xtran, ytran, ztran));
				ModelMatrix = ModelMatrix * TranslationMatrix;
			}

			t += 0.1;

			ModelMatrix = ModelMatrix * TranslationMatrix;
			glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
			glBindVertexArray(VertexArrayId[9]);	// draw Sphere
			glDrawElements(GL_TRIANGLES, NumIndices[9], GL_UNSIGNED_SHORT, (void*)0);

			_sleep(100);
		
			if ((-(grav+0.5) > (ytran+dropfrom)))
			//if(-grav > ytran+dropfrom)// && grav > xtran)
			{
			jump = false;
			t = 0.0f;
			posrl += horz;
			//posrl += -ztran;
			//posud += y;
			posud += ytran;

			grav = 0.0f;
			horz=0.0f;
			y = 1.1f;

			}
		//	*/
		}

	
		
		glBindVertexArray(0);

	}
	glUseProgram(0);
	// Draw GUI
	TwDraw();

	// Swap buffers
	glfwSwapBuffers(window);
	glfwPollEvents();
}

void pickObject(void)
{
	// Clear the screen in white
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	glUseProgram(pickingProgramID);
	{
		glm::mat4 ModelMatrix = glm::mat4(1.0); // TranslationMatrix * RotationMatrix;
		glm::mat4 MVP = gProjectionMatrix * gViewMatrix * ModelMatrix;

		// Send our transformation to the currently bound shader, in the "MVP" uniform
		glUniformMatrix4fv(PickingMatrixID, 1, GL_FALSE, &MVP[0][0]);

		// ATTN: DRAW YOUR PICKING SCENE HERE. REMEMBER TO SEND IN A DIFFERENT PICKING COLOR FOR EACH OBJECT BEFOREHAND
		//Move all Objects - moving base
		{
			trans[rl] = posrl;
			xtran += posrl;
			trans[ud] = posud;
			ztran += posud;

			glm::mat4x4 TranslationMatrix = glm::translate(glm::mat4x4(1.0f), trans);
			ModelMatrix = ModelMatrix * TranslationMatrix;

			glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
		}

		MVP = gProjectionMatrix * gViewMatrix * ModelMatrix;
		glUniformMatrix4fv(PickingMatrixID, 1, GL_FALSE, &MVP[0][0]);

		glUniform1f(pickingColorID, 2 / 255.0);
		createVAOs(Verts, Idcs, 2);
		glBindVertexArray(VertexArrayId[2]);	// draw Base
		glDrawElements(GL_TRIANGLES, NumIndices[2], GL_UNSIGNED_SHORT, (void*)0);

		// Rotate Top along Y axis
		{
			glm::mat4x4 RotationMatrix = glm::rotate(glm::mat4x4(1.0f), radians(rotT), glm::vec3(0.0, 1.0, 0.0));
			ModelMatrix = ModelMatrix * RotationMatrix;
			glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
		}

		MVP = gProjectionMatrix * gViewMatrix * ModelMatrix;
		glUniformMatrix4fv(PickingMatrixID, 1, GL_FALSE, &MVP[0][0]);

		glUniform1f(pickingColorID, 3 / 255.0);
		createVAOs(Verts1, Idcs1, 3);

		glBindVertexArray(VertexArrayId[3]);	// draw Top
		glDrawElements(GL_TRIANGLES, NumIndices[3], GL_UNSIGNED_SHORT, (void*)0);


		//Initial positioning
		{
			glm::mat4x4 TranslationMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(0.0, 1.77, 0.0));
			ytran += 1.77;
			glm::mat4 RotationMatrix = glm::rotate(glm::mat4x4(1.0f), radians(270.0f), glm::vec3(0.0, 1.0, 0.0));
			ModelMatrix = ModelMatrix * RotationMatrix * TranslationMatrix;
			glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
		}
		// Rotate Arm1 whlie one end is fixed to the top

		{
			glm::mat4x4 TranslationMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(-xtran, -ytran, -ztran));
			glm::mat4 RotationMatrix = glm::rotate(glm::mat4x4(1.0f), radians(rotA1), glm::vec3(1.0, 0.0, 0.0));
			ModelMatrix = ModelMatrix * RotationMatrix * TranslationMatrix;
			ang += rotA1;

			dropfrom += (1.77 * sin(rotA1));
			TranslationMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(xtran, ytran, ztran));
			ModelMatrix = ModelMatrix * TranslationMatrix;
			glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
		}

		MVP = gProjectionMatrix * gViewMatrix * ModelMatrix;
		glUniformMatrix4fv(PickingMatrixID, 1, GL_FALSE, &MVP[0][0]);

		glUniform1f(pickingColorID, 4 / 255.0);
		createVAOs(Verts2, Idcs2, 4);

		glBindVertexArray(VertexArrayId[4]);	// draw Arm1
		glDrawElements(GL_TRIANGLES, NumIndices[4], GL_UNSIGNED_SHORT, (void*)0);

		{
			glm::mat4x4 TranslationMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(0.0, 0.0, -1.75));
			ztran += -1.75f;
			ModelMatrix = ModelMatrix * TranslationMatrix;
			glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
		}

		MVP = gProjectionMatrix * gViewMatrix * ModelMatrix;
		glUniformMatrix4fv(PickingMatrixID, 1, GL_FALSE, &MVP[0][0]);

		glUniform1f(pickingColorID, 5 / 255.0);
		createVAOs(Verts3, Idcs3, 5);

		glBindVertexArray(VertexArrayId[5]);	// draw Joint
		glDrawElements(GL_TRIANGLES, NumIndices[5], GL_UNSIGNED_SHORT, (void*)0);

		//Rotate arm2 while one end is fixed to joint

		{
			glm::mat4x4 RotationMatrix = glm::rotate(glm::mat4x4(1.0f), radians(180.0f), glm::vec3(1.0, 0.0, 0.0));
			ModelMatrix = ModelMatrix * RotationMatrix;
			glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
		}

		{
			glm::mat4x4 TranslationMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(-xtran, -ytran, -ztran));
			glm::mat4 RotationMatrix = glm::rotate(glm::mat4x4(1.0f), radians(rotA2), glm::vec3(1.0, 0.0, 0.0));
			ModelMatrix = ModelMatrix * RotationMatrix * TranslationMatrix;
			ang += rotA2;
			dropfrom += (1.0 * sin(rotA2));
			TranslationMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(xtran, ytran, ztran));
			ModelMatrix = ModelMatrix * TranslationMatrix;
			glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
		}

		MVP = gProjectionMatrix * gViewMatrix * ModelMatrix;
		glUniformMatrix4fv(PickingMatrixID, 1, GL_FALSE, &MVP[0][0]);

		glUniform1f(pickingColorID, 6 / 255.0);
		createVAOs(Verts4, Idcs4, 6);

		glBindVertexArray(VertexArrayId[6]);	// draw Arm2
		glDrawElements(GL_TRIANGLES, NumIndices[6], GL_UNSIGNED_SHORT, (void*)0);

		//Rotate Pen while one end is fixed to arm2		
		{
			glm::mat4x4 RotationMatrix = glm::rotate(glm::mat4x4(1.0f), radians(90.0f), glm::vec3(1.0, 0.0, 0.0));
			glm::mat4x4 TranslationMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(0.0, 1.0, 0.0));
			ytran -= 1.25;
			ModelMatrix = ModelMatrix * TranslationMatrix * RotationMatrix;
			glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
		}

		{
			glm::mat4x4 TranslationMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(-xtran, -ytran, -ztran));
			glm::mat4 RotationMatrix = glm::rotate(glm::mat4x4(1.0f), radians(penrud), glm::vec3(1.0, 0.0, 0.0));
			ModelMatrix = ModelMatrix * RotationMatrix * TranslationMatrix;
			ang += penrud;
			dropfrom += (1.0 * sin(penrud));
			TranslationMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(xtran, ytran, ztran));
			ModelMatrix = ModelMatrix * TranslationMatrix;
			glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
		}

		{
			glm::mat4x4 TranslationMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(-xtran, -ytran, -ztran));
			glm::mat4 RotationMatrix = glm::rotate(glm::mat4x4(1.0f), radians(penrrl), glm::vec3(0.0, 0.0, 1.0));
			ModelMatrix = ModelMatrix * RotationMatrix * TranslationMatrix;
			TranslationMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(xtran, ytran, ztran));
			ModelMatrix = ModelMatrix * TranslationMatrix;
			glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
		}

		{
			glm::mat4x4 TranslationMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(-xtran, -ytran, -ztran));
			glm::mat4 RotationMatrix = glm::rotate(glm::mat4x4(1.0f), radians(penrev), glm::vec3(0.0, 1.0, 0.0));
			ModelMatrix = ModelMatrix * RotationMatrix * TranslationMatrix;
			TranslationMatrix = glm::translate(glm::mat4(1.0f), glm::vec3(xtran, ytran, ztran));
			ModelMatrix = ModelMatrix * TranslationMatrix;
			glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
		}

		MVP = gProjectionMatrix * gViewMatrix * ModelMatrix;
		glUniformMatrix4fv(PickingMatrixID, 1, GL_FALSE, &MVP[0][0]);

		glUniform1f(pickingColorID, 7 / 255.0);
		createVAOs(Verts5, Idcs5, 7);

		glBindVertexArray(VertexArrayId[7]);	// draw Pen
		glDrawElements(GL_TRIANGLES, NumIndices[7], GL_UNSIGNED_SHORT, (void*)0);

		glBindVertexArray(VertexArrayId[8]);	// draw Button
		glDrawElements(GL_TRIANGLES, NumIndices[8], GL_UNSIGNED_SHORT, (void*)0);
		
		glBindVertexArray(0);


		
	}
	glUseProgram(0);
	// Wait until all the pending drawing commands are really done.
	// Ultra-mega-over slow ! 
	// There are usually a long time between glDrawElements() and
	// all the fragments completely rasterized.
	glFlush();
	glFinish();

	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	// Read the pixel at the center of the screen.
	// You can also use glfwGetMousePos().
	// Ultra-mega-over slow too, even for 1 pixel, 
	// because the framebuffer is on the GPU.
	double xpos, ypos;
	glfwGetCursorPos(window, &xpos, &ypos);
	unsigned char data[4];
	glReadPixels(xpos, window_height - ypos, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, data); // OpenGL renders with (0,0) on bottom, mouse reports with (0,0) on top

																					 // Convert the color back to an integer ID
	gPickedIndex = int(data[0]);
	if (gPickedIndex == 255) { // Full white, must be the background !
		gMessage = "background";
	}
	else {
		std::ostringstream oss;
		
		switch (gPickedIndex)
		{
		case 2:
			if (base == true)
			{
				base = false;
			}
			else
			{
				base = true;
				oss << "Base";
			}
			break;
		case 3:
			if (top == true)
			{
				top = false;
			}
			else
			{
				top = true;
				oss << "Top";
			}
			break;

		case 4:
			if (arm1 == true)
			{
				arm1 = false;
			}
			else
			{
				arm1 = true;
				oss << "Arm 1";
			}
			break;
		case 6:
			if (arm2 == true)
			{
				arm2 = false;
			}
			else
			{
				arm2 = true;
				oss << "Arm 2";
			}
			break;
		case 7:
			if (pen == true)
			{
				pen = false;
			}
			else
			{
				pen = true;
				oss << "Pen";
			}
			break;
		default:
			oss << "Object " << gPickedIndex;
			break;
		}
		gMessage = oss.str();
	}

	// Uncomment these lines to see the picking shader in effect
	//glfwSwapBuffers(window);
	//continue; // skips the normal rendering
	
}

int initWindow(void)
{
	// Initialise GLFW
	if (!glfwInit()) {
		fprintf(stderr, "Failed to initialize GLFW\n");
		return -1;
	}

	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_FALSE);

	// Open a window and create its OpenGL context
	window = glfwCreateWindow(window_width, window_height, "Nagavara Janakiprasad,Adarsh (66975994)", NULL, NULL);
	if (window == NULL) {
		fprintf(stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n");
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);

	// Initialize GLEW
	glewExperimental = true; // Needed for core profile
	if (glewInit() != GLEW_OK) {
		fprintf(stderr, "Failed to initialize GLEW\n");
		return -1;
	}

	// Initialize the GUI
	TwInit(TW_OPENGL_CORE, NULL);
	TwWindowSize(window_width, window_height);
	TwBar * GUI = TwNewBar("Picking");
	TwSetParam(GUI, NULL, "refresh", TW_PARAM_CSTRING, 1, "0.1");
	TwAddVarRW(GUI, "Last picked object", TW_TYPE_STDSTRING, &gMessage, NULL);

	// Set up inputs
	glfwSetCursorPos(window, window_width / 2, window_height / 2);
	glfwSetKeyCallback(window, keyCallback);
	glfwSetMouseButtonCallback(window, mouseCallback);
	glfwSetScrollCallback(window,scrollCallback);

	return 0;
}

void initOpenGL(void)
{

	// Enable depth test
	glEnable(GL_DEPTH_TEST);
	// Accept fragment if it closer to the camera than the former one
	glDepthFunc(GL_LESS);
	// Cull triangles which normal is not towards the camera
	glEnable(GL_CULL_FACE);

	glEnable(GL_LIGHTING);

	// Projection matrix : 45° Field of View, 4:3 ratio, display range : 0.1 unit <-> 100 units
	gProjectionMatrix = glm::perspective(45.0f, 4.0f / 3.0f, 0.1f, 100.0f);
	// Or, for an ortho camera :
	//gProjectionMatrix = glm::ortho(-4.0f, 4.0f, -3.0f, 3.0f, 0.0f, 100.0f); // In world coordinates

	// Camera matrix
	gViewMatrix = glm::lookAt(glm::vec3(10.0, 10.0, 10.0f),	// eye
		glm::vec3(0.0, 0.0, 0.0),	// center
		glm::vec3(0.0, 1.0, 0.0));	// up

									// Create and compile our GLSL program from the shaders
	programID = LoadShaders("StandardShading.vertexshader", "StandardShading.fragmentshader");
	pickingProgramID = LoadShaders("Picking.vertexshader", "Picking.fragmentshader");

	// Get a handle for our "MVP" uniform
	MatrixID = glGetUniformLocation(programID, "MVP");
	ModelMatrixID = glGetUniformLocation(programID, "M");
	ViewMatrixID = glGetUniformLocation(programID, "V");
	ProjMatrixID = glGetUniformLocation(programID, "P");

	PickingMatrixID = glGetUniformLocation(pickingProgramID, "MVP");
	// Get a handle for our "pickingColorID" uniform
	pickingColorID = glGetUniformLocation(pickingProgramID, "PickingColor");
	// Get a handle for our "LightPosition" uniform
	LightID = glGetUniformLocation(programID, "LightPosition_worldspace");
	LightID1 = glGetUniformLocation(programID, "LightPosition_worldspace1");
	//glUniform3f(LightColor,1.0,0.6,1.0);
	createObjects();
}

void createVAOs(Vertex Vertices[], unsigned short Indices[], int ObjectId) {

	GLenum ErrorCheckValue = glGetError();
	const size_t VertexSize = sizeof(Vertices[0]);
	const size_t RgbOffset = sizeof(Vertices[0].Position);
	const size_t Normaloffset = sizeof(Vertices[0].Color) + RgbOffset;

	// Create Vertex Array Object
	glGenVertexArrays(1, &VertexArrayId[ObjectId]);	//
	glBindVertexArray(VertexArrayId[ObjectId]);		//

													// Create Buffer for vertex data
	glGenBuffers(1, &VertexBufferId[ObjectId]);
	glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[ObjectId]);
	glBufferData(GL_ARRAY_BUFFER, VertexBufferSize[ObjectId], Vertices, GL_STATIC_DRAW);

	// Create Buffer for indices
	if (Indices != NULL) {
		glGenBuffers(1, &IndexBufferId[ObjectId]);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, IndexBufferId[ObjectId]);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, IndexBufferSize[ObjectId], Indices, GL_STATIC_DRAW);
	}

	// Assign vertex attributes
	glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, VertexSize, 0);
	glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, VertexSize, (GLvoid*)RgbOffset);
	glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, VertexSize, (GLvoid*)Normaloffset);

	glEnableVertexAttribArray(0);	// position
	glEnableVertexAttribArray(1);	// color
	glEnableVertexAttribArray(2);	// normal

									// Disable our Vertex Buffer Object 
	glBindVertexArray(0);

	ErrorCheckValue = glGetError();
	if (ErrorCheckValue != GL_NO_ERROR)
	{
		fprintf(
			stderr,
			"ERROR: Could not create a VBO: %s \n",
			gluErrorString(ErrorCheckValue)
		);
	}
}

void cleanup(void)
{
	// Cleanup VBO and shader
	for (int i = 0; i < NumObjects; i++) {
		glDeleteBuffers(1, &VertexBufferId[i]);
		glDeleteBuffers(1, &IndexBufferId[i]);
		glDeleteVertexArrays(1, &VertexArrayId[i]);
	}
	glDeleteProgram(programID);
	glDeleteProgram(pickingProgramID);

	// Close OpenGL window and terminate GLFW
	glfwTerminate();
}

static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	// ATTN: MODIFY AS APPROPRIATE
	if (action == GLFW_PRESS) {
		switch (key)
		{
		case GLFW_KEY_2:
		{			
			if (arm2 == false)
			{
				arm2 = true; 
				gMessage = "Arm 2";
			}
			else			
				arm2 = false;			
			break;
		}
		case GLFW_KEY_1:
		{			
			if (arm1 == false)
			{
				arm1 = true;
				gMessage = "Arm 1";
			}
			else
				arm1 = false;
			break;
		}			
		case GLFW_KEY_T:
		{
		
			if (top == false)
			{
				top = true;						
				gMessage = "Top";
				break;
			}
			else
			{
				top = false;
			}
			break;
		}			
		case GLFW_KEY_P:
		{
			if (pen == false)
			{
				pen = true;
				gMessage = "Pen";
			}
			else
				pen = false;
			break;
		}
		case GLFW_KEY_B:
		{
			float color[] = { baseColor[0], baseColor[1], baseColor[2],1.0 };
			if (base == false)
			{
				base = true;
				gMessage = "Base";
			}
			else
			{
				base = false;
			}				
			break;
		}			
		case GLFW_KEY_C:
		{
			if (camera == false)
			{
				camera = true;
				gMessage = "Camera";
			}
			else
				camera = false;
			break;
		}

		case GLFW_KEY_DOWN:
		{
			
			if (camera == true)
			{
				theta += 5.0;
				renderScene();
				break;
			}

			if (base == true)
			{	
				if(posud<5)
					posud += 1.0f;	
				break;
			}		

			if (arm1 == true)
			{
				if(rotA1 >= -10.0f)
				rotA1 -= 5.0f;
				break;
			}

			if (arm2 == true)
			{
				rotA2 -= 5.0f;
				break;
			}

			if (pen == true)
			{
				penrud -=5.0f;
				break;
			}
			
			break;
		}
		case GLFW_KEY_UP:
		{
			if (camera == true)
			{
				theta -= 5.0;
				renderScene();				
				break;
			}

			if (base == true)
			{
				if (posud >-5)
					posud -= 1.0f;				
				break;
			}
			
			if (arm1 == true)
			{
				rotA1 += 5.0f;
				break;
			}

			if (arm2 == true)
			{
				rotA2 += 5.0f;
				break;
			}

			if (pen == true)
			{
				penrud += 5.0f;
				break;
			}			

			break;
		}

		case GLFW_KEY_RIGHT:
		{
			if (pen == true && glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
				glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS)
			{
				penrev += 45.0f;
				break;
			}
			if (camera == true)
			{
				Phi -= 5.0;
				renderScene();
				break;
			}

			if (base == true)
			{
				if (posrl<5)
					posrl += 1.0f;				
				break;
			}

			if (top == true)
			{
				rotT -= 5.0f;
				break;
			}

			if (pen == true)
			{
				penrrl += 5.0f;
				break;
			}

			break;
		}
		case GLFW_KEY_LEFT:
		{
			if (pen == true && glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
				glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS)
			{
				penrev -= 45.0f;
				break;
			}

			if (camera == true)
			{
				Phi += 5.0;
				renderScene();
				break;
			}

			if (base == true)
			{
				if (posrl > -5)
					posrl -= 1.0f;				
				break;
			}

			if (top == true)
			{
				rotT += 5.0f;
				break;
			}

			if (pen == true)
			{
				penrrl -= 5.0f;
				break;
			}

			break;
		}

		case GLFW_KEY_J:
		{
			if (jump == false)
			{
				jump = true;			
				grav = 0.0f;
			}
			//else
			//	jump = false;
			break;
		}
		default:
			break;
		}		
	}
}

static void mouseCallback(GLFWwindow* window, int button, int action, int mods)
{
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
		pickObject();
	}
}

static void scrollCallback(GLFWwindow * window, double xoffset, double yoffset)
{
	if (yoffset <  15 && yoffset > -10)
	{
		radius = radius - yoffset;
		renderScene();
	}

}

int main(void)
{
	// initialize window
	int errorCode = initWindow();
	if (errorCode != 0)
		return errorCode;

	// initialize OpenGL pipeline
	initOpenGL();

	// For speed computation
	double lastTime = glfwGetTime();
	int nbFrames = 0;
	do {
		//// Measure speed
		//double currentTime = glfwGetTime();
		//nbFrames++;
		//if (currentTime - lastTime >= 1.0){ // If last prinf() was more than 1sec ago
		//	// printf and reset
		//	printf("%f ms/frame\n", 1000.0 / double(nbFrames));
		//	nbFrames = 0;
		//	lastTime += 1.0;
		//}

		if (animation) {
			phi += 0.01;
			if (phi > 360)
				phi -= 360;
		}

		// DRAWING POINTS
		renderScene();


	} // Check if the ESC key was pressed or the window was closed
	while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
		glfwWindowShouldClose(window) == 0);

	cleanup();

	return 0;
}