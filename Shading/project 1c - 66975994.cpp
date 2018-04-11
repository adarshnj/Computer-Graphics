// Include standard headers
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <array>
#include <sstream>
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

typedef struct Vertex {
	float XYZW[4];
	float RGBA[4];
	void SetCoords(float *coords) {
		XYZW[0] = coords[0];
		XYZW[1] = coords[1];
		XYZW[2] = coords[2];
		XYZW[3] = coords[3];
	}
	void SetColor(float *color) {
		RGBA[0] = color[0];
		RGBA[1] = color[1];
		RGBA[2] = color[2];
		RGBA[3] = color[3];
	}
};

// ATTN: USE POINT STRUCTS FOR EASIER COMPUTATIONS
typedef struct point {
	float x, y, z;
	point(const float x = 0, const float y = 0, const float z = 0) : x(x), y(y), z(z) {};
	point(float *coords) : x(coords[0]), y(coords[1]), z(coords[2]) {};
	point operator -(const point& a)const {
		return point(x - a.x, y - a.y, z - a.z);
	}
	point operator +(const point& a)const {
		return point(x + a.x, y + a.y, z + a.z);
	}
	point operator *(const float& a)const {
		return point(x*a, y*a, z*a);
	}
	point operator /(const float& a)const {
		return point(x / a, y / a, z / a);
	}
	float* toArray() {
		float array[] = { x, y, z, 1.0f };
		return array;
	}
};

// function prototypes
int initWindow(void);
void initOpenGL(void);
void createVAOs(Vertex[], unsigned short[], size_t, size_t, int);
void createObjects(void);
void pickVertex(void);
void moveVertex(void);
void drawScene(void);
void cleanup(void);
static void mouseCallback(GLFWwindow*, int, int, int);
static void keyCallback(GLFWwindow*, int, int, int, int);

// GLOBAL VARIABLES
GLFWwindow* window;
const GLuint window_width = 1024, window_height = 768;

glm::mat4 gProjectionMatrix;
glm::mat4 gViewMatrix;
glm::mat4 gViewMatrixz;

GLuint gPickedIndex;
std::string gMessage;

GLuint programID;
GLuint pickingProgramID;

// ATTN: INCREASE THIS NUMBER AS YOU CREATE NEW OBJECTS
const GLuint NumObjects = 9;	// number of different "objects" to be drawn
GLuint VertexArrayId[NumObjects] = { 0 };
GLuint VertexBufferId[NumObjects] = { 0 };
GLuint IndexBufferId[NumObjects] = { 0 };
size_t NumVert[NumObjects] = { 0 };

GLuint MatrixID;
GLuint ViewMatrixID;
GLuint ModelMatrixID;
GLuint PickingMatrixID;
GLuint pickingColorArrayID;
GLuint pickingColorID;
GLuint LightID;

// Define objects
Vertex Vertices[] =
{
	{ { 0.0f, 0.0f, 0.0f, 1.0f },{ 0.0f, 1.0f, 0.0f, 1.0f } }, // 0
	{ { -0.7f, 0.7f, 0.0f, 1.0f },{ 0.0f, 1.0f, 0.0f, 1.0f } }, // 1	
	{ { -1.0f, 0.0f, 0.0f, 1.0f },{ 0.0f, 1.0f, 0.0f, 1.0f } }, // 2
	{ { -0.7f, -0.7f, 0.0f, 1.0f },{ 0.0f, 1.0f, 0.0f, 1.0f } }, // 3
	{ { 0.0f, -1.0f, 0.0f, 1.0f },{ 0.0f, 1.0f, 0.0f, 1.0f } }, // 4
	{ { 0.7f, -0.7f, 0.0f, 1.0f },{ 0.0f, 1.0f, 0.0f, 1.0f } }, // 5
	{ { 1.0f, 0.0f, 0.0f, 1.0f },{ 0.0f, 1.0f, 0.0f, 1.0f } }, // 6
	{ { 0.7f, 0.7f, 0.0f, 1.0f },{ 0.0f, 1.0f, 0.0f, 1.0f } }, // 7

};

unsigned short Indices[] = {
	0, 1, 2, 3, 4, 5, 6, 7
};

const size_t IndexCount = sizeof(Indices) / sizeof(unsigned short);
// ATTN: DON'T FORGET TO INCREASE THE ARRAY SIZE IN THE PICKING VERTEX SHADER WHEN YOU ADD MORE PICKING COLORS
float pickingColor[IndexCount] = { 0 / 255.0f, 1 / 255.0f, 2 / 255.0f, 3 / 255.0f ,
4 / 255.0f, 5 / 255.0f, 6 / 255.0f, 7 / 255.0f };

// ATTN: ADD YOU PER-OBJECT GLOBAL ARRAY DEFINITIONS HERE

int perSec = 0;
int pnt = 0;
bool dview = false;
Vertex Lineloop[10];
unsigned short LlIndices[10];
Vertex interim[256];
int icount = 0;
Vertex Subdivision[128];
int sub = 0;
int scount = 0;
int count_1 = 0;	//track count of '1's pressed
bool subd = false;
unsigned short SdIndices[128];
Vertex Controlpoints[40];
int ccount = 0;
unsigned short CpIndices[40];
Vertex pcoor[400];
Vertex Bezier[400];
int bcount = 0;
bool bez = false;
unsigned short BIndices[400];
bool ydot = false;
int dotP = 0;
Vertex Point[2];
Vertex Tangent[4];
Vertex Normal[4];
Vertex Binormal[4];



point sub1(int k, int i)
{
	point p, p1, p2, p3;
	if (k == 0)
	{
		p = interim[i].XYZW;
		return p;
	}
	else
	{
		p1 = sub1(k - 1, (i - 2 + icount) % icount);
		p2 = sub1(k - 1, (i - 1 + icount) % icount);
		p3 = sub1(k - 1, i);

		p = p1 + (p2 * 10) + (p3 * 5);
		return p / 16;
	}
}

point sub2(int k, int i)
{
	point p, p1, p2, p3;
	if (k == 0)
	{
		p = interim[i].XYZW;
		return p;
	}
	else
	{
		p1 = sub2(k - 1, (i - 1 + icount) % icount);
		p2 = sub2(k - 1, i);
		p3 = sub2(k - 1, (i + 1) % icount);

		p = (p1 * 5) + (p2 * 10) + p3;
		return p / 16;
	}
}

void subdivision()
{
	int k = count_1;
	float color[] = { 0.0f, 0.9f, 1.0f, 1.0f };
	scount = 0;

	for (int i = 0; i < IndexCount; i++)
	{
		Subdivision[i] = Vertices[i];
		scount++;
	}

	for (sub = 1; sub < k + 1; sub++)
	{
		icount = 0;
		for (int co = 0; co < scount; co++)
		{
			interim[co] = Subdivision[co];
			icount++;
		}

		scount = 0;
		for (int i = 0; i < icount; i++)
		{
			point p, p1;
			p = sub1(1, i);
			Subdivision[(i * 2) + 0].SetCoords(p.toArray());
			Subdivision[(i * 2) + 0].SetColor(color);
			scount++;

			p1 = sub2(1, i);
			Subdivision[(i * 2) + 1].SetCoords(p1.toArray());
			Subdivision[(i * 2) + 1].SetColor(color);
			scount++;

		}
	}
	for (int i = 0; i < icount; i++)
		SdIndices[i] = i;
}


point c1(int i)
{
	point t, t1, t2, t3, t4;

	t1 = Vertices[(i - 2 + IndexCount) % IndexCount].XYZW;
	t2 = Vertices[(i - 1 + IndexCount) % IndexCount].XYZW;
	t3 = Vertices[i % IndexCount].XYZW;
	t4 = Vertices[(i + 1) % IndexCount].XYZW;

	t = t1 + (t2 * 11) + (t3 * 11) + t4;

	return t / 24;

}

point c2(int i)
{
	point t, t1, t2, t3;

	t1 = Vertices[(i - 1 + IndexCount) % IndexCount].XYZW;
	t2 = Vertices[i % IndexCount].XYZW;
	t3 = Vertices[(i + 1) % IndexCount].XYZW;

	t = (t1 * 4) + (t2 * 7) + t3;

	return t / 12;

}

point c3(int i)
{
	point t, t1, t2, t3;

	t1 = Vertices[(i - 1 + IndexCount) % IndexCount].XYZW;
	t2 = Vertices[i % IndexCount].XYZW;
	t3 = Vertices[(i + 1) % IndexCount].XYZW;

	t = (t1 * 4) + (t2 * 16) + (t3 * 4);

	return t / 24;

}

point c4(int i)
{
	point t, t1, t2, t3, t4;

	t1 = Vertices[(i - 1 + IndexCount) % IndexCount].XYZW;
	t2 = Vertices[i % IndexCount].XYZW;
	t3 = Vertices[(i + 1) % IndexCount].XYZW;

	t = t1 + (t2 * 7) + (t3 * 4);

	return t / 12;

}

point c5(int i)
{
	point t, t1, t2, t3, t4;

	t1 = Vertices[(i - 1 + IndexCount) % IndexCount].XYZW;
	t2 = Vertices[i % IndexCount].XYZW;
	t3 = Vertices[(i + 1) % IndexCount].XYZW;
	t4 = Vertices[(i + 2) % IndexCount].XYZW;

	t = t1 + (t2 * 11) + (t3 * 11) + t4;
	return t / 24;

}

void movep(void)
{
	dotP = dotP++ % bcount;
	point pdot, pfirst, psecond, pthird, pdot1;
	float yellow[] = { 1.0f, 1.0f, 0.0f, 1.0f };
	float red[] = { 1.0f,0.0f,0.0f,1.0f };
	float blue[] = { 0.0f,0.4f,1.0f,1.0f };
	float green[] = { 0.0f,1.0f,0.0f,1.0f };

	pdot = Bezier[dotP].XYZW;
	pfirst = pcoor[dotP].XYZW;
	Point[0].SetCoords(pdot.toArray());
	Point[0].SetColor(yellow);

	float angle = atan((pfirst.y - pdot.y) / (pfirst.x - pdot.x));

	pfirst.x = pdot.x - (0.5 * cos(angle));
	pfirst.y = pdot.y - (0.5 * sin(angle));

	Tangent[0].SetCoords(pdot.toArray());
	Tangent[1].SetCoords(pfirst.toArray());
	Tangent[0].SetColor(red);
	Tangent[1].SetColor(red);

	psecond.x = pfirst.y * -1;
	psecond.y = pfirst.x;

	float lenn = sqrt(pow(pdot.x - psecond.x, 2) + pow(pdot.y - psecond.y, 2) + pow(pdot.z - psecond.z, 2));
	psecond.x = ((pdot.x - psecond.x) / lenn) * 0.5;
	psecond.y = ((pdot.y - psecond.y) / lenn) * 0.5;
	psecond.z = ((pdot.z - psecond.z) / lenn) * 0.5;

	Normal[0].SetCoords(pdot.toArray());
	Normal[1].SetCoords(psecond.toArray());
	Normal[0].SetColor(blue);
	Normal[1].SetColor(blue);


	pthird.x = pfirst.y * psecond.z - pfirst.z * psecond.y;
	pthird.y = pfirst.z * psecond.x - pfirst.x * psecond.z;
	pthird.z = pfirst.x * psecond.y - pfirst.y * psecond.x;

	float len = sqrt(pow(pdot.x - pthird.x, 2) + pow(pdot.y - pthird.y, 2) + pow(pdot.z - pthird.z, 2));
	pthird.x = ((pdot.x - pthird.x) / len) * 0.5;
	pthird.y = ((pdot.y - pthird.y) / len) * 0.5;
	pthird.z = ((pdot.z - pthird.z) / len) * 0.5;

	Binormal[0].SetCoords(pdot.toArray());
	Binormal[1].SetCoords(pthird.toArray());
	Binormal[0].SetColor(green);
	Binormal[1].SetColor(green);

	drawScene();
	_sleep(50);


}

void bezier(void)
{
	point p, p0, p1, p2, p3, p4, pc;
	float color1[] = { 1.0f, 1.0f, 0.0f, 1.0f };

	bcount = 0;
	for (int i = 0; i < ccount; i += 5)
	{
		int N = 20;   //Size of the array declared to draw bezier curve restricts to increment over 40
		float tesel = 1.0f / N;

		p0 = Controlpoints[i].XYZW;
		p1 = Controlpoints[(i + 1)].XYZW;
		p2 = Controlpoints[(i + 2)].XYZW;
		p3 = Controlpoints[(i + 3)].XYZW;
		p4 = Controlpoints[(i + 4)].XYZW;

		for (float tk = 0; tk < 1; tk += tesel)
		{
			float u = 1 - tk;
			p = (p0*pow(u, 4)) + (p1*(4 * tk*pow(u, 3))) + (p2*(6 * tk*tk*u*u)) + (p3*(4 * u*pow(tk, 3))) + (p4*(pow(tk, 4)));
			Bezier[bcount].SetCoords(p.toArray());

			// Vector coordinates for tangent of each point
			pc = (p0*(u*u*u)) + (p1*(3 * u*u*tk)) + (p2*(3 * u*tk*tk)) + (p3*(tk*tk*tk));
			pcoor[bcount].SetCoords(pc.toArray());

			bcount++;
		}
	}

	for (int i = 0; i < bcount; i++)
	{
		Bezier[i].SetColor(color1);
		BIndices[i] = i;
	}

	if (ydot == true)// && dview == true)
		movep();

}

void deCastelJau(void)
{
	point p, p1, p2, p3, p4;
	float color[] = { 1.0f, 0.0f, 0.0f, 1.0f };
	ccount = 0;
	int num = 5;

	for (int i = 0; i < IndexCount; i++)
	{
		p = c1(i);
		Controlpoints[(i * num) + 0].SetCoords(p.toArray());
		ccount++;

		p1 = c2(i);
		Controlpoints[(i * num) + 1].SetCoords(p1.toArray());
		ccount++;

		p2 = c3(i);
		Controlpoints[(i * num) + 2].SetCoords(p2.toArray());
		ccount++;

		p3 = c4(i);
		Controlpoints[(i * num) + 3].SetCoords(p3.toArray());
		ccount++;

		p4 = c5(i);
		Controlpoints[(i * num) + 4].SetCoords(p4.toArray());
		ccount++;

	}

	for (int i = 0; i < ccount; i++)
	{
		Controlpoints[i].SetColor(color);
		CpIndices[i] = i;
	}
	bezier();

}

void createObjects(void)
{
	// ATTN: DERIVE YOUR NEW OBJECTS HERE:
	// each has one vertices {pos;color} and one indices array (no picking needed here)

	for (int i = 0; i <= IndexCount; i++)
	{
		LlIndices[i] = i;
		float color[] = { 1.0f, 1.0f, 1.0f, 1.0f };
		Lineloop[i] = Vertices[i%IndexCount];
		Lineloop[i].SetColor(color);
	}

	subdivision();
	deCastelJau();

}

void usePrg()
{
	glUseProgram(programID);
	{
		glm::mat4 ModelMatrix = glm::mat4(1.0); // TranslationMatrix * RotationMatrix;
		glm::mat4 MVP = gProjectionMatrix * gViewMatrix * ModelMatrix;

		// Send our transformation to the currently bound shader, 
		// in the "MVP" uniform
		glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);
		glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
		glUniformMatrix4fv(ViewMatrixID, 1, GL_FALSE, &gViewMatrix[0][0]);
		glm::vec3 lightPos = glm::vec3(4, 4, 4);
		glUniform3f(LightID, lightPos.x, lightPos.y, lightPos.z);

		glEnable(GL_PROGRAM_POINT_SIZE);

		glBindVertexArray(VertexArrayId[0]);	// draw Vertices
		glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[0]);
		glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(Vertices), Vertices);				// update buffer data
																						//glDrawElements(GL_LINE_LOOP, NumVert[0], GL_UNSIGNED_SHORT, (void*)0);
		glDrawElements(GL_POINTS, NumVert[0], GL_UNSIGNED_SHORT, (void*)0);
		// ATTN: OTHER BINDING AND DRAWING COMMANDS GO HERE, one set per object:
		//glBindVertexArray(VertexArrayId[<x>]); etc etc

		//Object 1 - White line joining the initial vertices
		glBindVertexArray(VertexArrayId[1]);	// draw Vertices
		glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[1]);
		glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(Lineloop), Lineloop);				// update buffer data
		glDrawArrays(GL_LINE_STRIP, 0, IndexCount + 1);

		//Object 2 - Curve by Subdivision
		glBindVertexArray(VertexArrayId[2]);	// draw Vertices
		glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[2]);
		glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(Subdivision), Subdivision);				// update buffer data		
		glDrawArrays(GL_LINE_LOOP, 0, scount);

		if (ydot == true)
		{
			glBindVertexArray(VertexArrayId[5]);	// draw Yellow dot
			glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[5]);
			glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(Point), Point);				// update buffer data			  
			glDrawElements(GL_POINTS, NumVert[0], GL_UNSIGNED_SHORT, (void*)0);

			glBindVertexArray(VertexArrayId[6]);	// draw Tangent
			glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[6]);
			glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(Tangent), Tangent);				// update buffer data			  
			glDrawArrays(GL_LINE_STRIP, 0, 2);

			glBindVertexArray(VertexArrayId[7]);	// draw Normal
			glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[7]);
			glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(Normal), Normal);				// update buffer data			  
			glDrawArrays(GL_LINE_STRIP, 0, 2);

			glBindVertexArray(VertexArrayId[8]);	// draw Binormal
			glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[8]);
			glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(Binormal), Binormal);				// update buffer data			  
			glDrawArrays(GL_LINE_STRIP, 0, 2);
		}

		if (bez == true)
		{


			//Object 3 - Control Points
			glBindVertexArray(VertexArrayId[3]);	// draw Vertices
			glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[3]);
			glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(Controlpoints), Controlpoints);
			glDrawArrays(GL_POINTS, 0, ccount);
			//Object 4 - Bezier Curve

			glBindVertexArray(VertexArrayId[4]);	// draw Vertices
			glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[4]);
			glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(Bezier), Bezier);
			glDrawArrays(GL_LINE_LOOP, 0, bcount);


		}

		glBindVertexArray(0);
	}
}

void drawScene(void)
{
	// Dark blue background
	glClearColor(0.0f, 0.0f, 0.4f, 0.0f);
	// Re-clear the screen for real rendering
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	if (dview == true)
	{
		//	glClear(GL_COLOR_BUFFER_BIT);		 		

		gProjectionMatrix = glm::ortho(-4.0f, 4.0f, -1.5f, 1.5f, 0.0f, 100.0f); // In world coordinates
		glViewport(0, window_height / 2, window_width, window_height / 2);
		// Camera matrix
		gViewMatrix = glm::lookAt(
			glm::vec3(0, 0, -5), // Camera is at (4,3,3), in World Space
			glm::vec3(0, 0, 0), // and looks at the origin
			glm::vec3(0, 1, 0)  // Head is up (set to 0,-1,0 to look upside-down)
		);

		usePrg();

		gProjectionMatrix = glm::ortho(-4.0f, 4.0f, -1.5f, 1.5f, 0.0f, 100.0f); // In world coordinates
		glViewport(0, 0, window_width, window_height / 2);
		// Camera matrix
		gViewMatrix = glm::lookAt(
			glm::vec3(5, 0, 0), // Camera is at (4,3,3), in World Space
			glm::vec3(0, 0, 0), // and looks at the origin
			glm::vec3(0, 1, 0)  // Head is up (set to 0,-1,0 to look upside-down)
		);


		usePrg();
	}
	else
	{
		gProjectionMatrix = glm::ortho(-4.0f, 4.0f, -3.0f, 3.0f, 0.0f, 100.0f); // In world coordinates
		glViewport(0, 0, window_width, window_height);
		// Camera matrix
		gViewMatrix = glm::lookAt(
			glm::vec3(0, 0, -5), // Camera is at (4,3,3), in World Space
			glm::vec3(0, 0, 0), // and looks at the origin
			glm::vec3(0, 1, 0)  // Head is up (set to 0,-1,0 to look upside-down)
		);

		usePrg();
	}

	glUseProgram(0);
	// Draw GUI
	TwDraw();

	// Swap buffers
	glfwSwapBuffers(window);
	glfwPollEvents();
}

void pickVertex(void)
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
		glUniform1fv(pickingColorArrayID, NumVert[0], pickingColor);	// here we pass in the picking marker array

																		// Draw the ponts
		glEnable(GL_PROGRAM_POINT_SIZE);
		glBindVertexArray(VertexArrayId[0]);
		glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[0]);
		glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(Vertices), Vertices);	// update buffer data
		glDrawElements(GL_POINTS, NumVert[0], GL_UNSIGNED_SHORT, (void*)0);
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

	// Uncomment these lines to see the picking shader in effect
	//glfwSwapBuffers(window);
	//continue; // skips the normal rendering
}

// fill this function in!
void moveVertex(void)
{
	double xpos, ypos;
	glm::mat4 ModelMatrix = glm::mat4(1.0);
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	glm::vec4 vp = glm::vec4(viewport[0], viewport[1], viewport[2], viewport[3]);

	// retrieve your cursor position
	// get your world coordinates
	// move points

	if (gPickedIndex == 255) { // Full white, must be the background !
		gMessage = "background";
	}
	else {
		std::ostringstream oss;
		oss << "point " << gPickedIndex;
		gMessage = oss.str();
	}
	int dist = 0;
	if (gPickedIndex != 255)
	{
		if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
			glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS)
		{
			glfwGetCursorPos(window, &xpos, &ypos);
			gProjectionMatrix = glm::ortho(-4.0f, 4.0f, -3.0f, 3.0f, 0.0f, 100.0f);
			glm::vec3 windowPos = glm::unProject(glm::vec3(xpos, ypos, 0.0),
				ModelMatrix, gProjectionMatrix, glm::vec4(0.0, 0.0, 1024, 768));

			Vertices[gPickedIndex].XYZW[2] = (windowPos.x + Vertices[gPickedIndex].XYZW[0]) * -1; // x axis of the point is the origin for z axis
		}
		else
		{
			glfwGetCursorPos(window, &xpos, &ypos);
			gProjectionMatrix = glm::ortho(-4.0f, 4.0f, -3.0f, 3.0f, 0.0f, 100.0f);
			glm::vec3 windowPos = glm::unProject(glm::vec3(xpos, ypos, 0.0), ModelMatrix, gProjectionMatrix, glm::vec4(0.0, 0.0, 1024, 768));
			Vertices[gPickedIndex].XYZW[0] = windowPos.x * -1;		// Invert 
			Vertices[gPickedIndex].XYZW[1] = windowPos.y * -1;
		}

	}

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
	glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
	glfwSetCursorPos(window, window_width, window_height);
	glfwSetMouseButtonCallback(window, mouseCallback);
	glfwSetKeyCallback(window, keyCallback);

	return 0;
}

void initOpenGL(void)
{
	// Dark blue background
	glClearColor(0.0f, 0.0f, 0.4f, 0.0f);

	// Enable depth test
	glEnable(GL_DEPTH_TEST);
	// Accept fragment if it closer to the camera than the former one
	glDepthFunc(GL_LESS);
	// Cull triangles which normal is not towards the camera
	glEnable(GL_CULL_FACE);


	// Projection matrix : 45° Field of View, 4:3 ratio, display range : 0.1 unit <-> 100 units
	//glm::mat4 ProjectionMatrix = glm::perspective(45.0f, 4.0f / 1.5f, 0.1f, 100.0f);
	// Or, for an ortho camera :

	gProjectionMatrix = glm::ortho(-4.0f, 4.0f, -3.0f, 3.0f, 0.0f, 100.0f); // In world coordinates

																			// Camera matrix
	gViewMatrix = glm::lookAt(
		glm::vec3(0, 0, -5), // Camera is at (4,3,3), in World Space
		glm::vec3(0, 0, 0), // and looks at the origin
		glm::vec3(0, 1, 0)  // Head is up (set to 0,-1,0 to look upside-down)
	);

	// Create and compile our GLSL program from the shaders
	programID = LoadShaders("StandardShading.vertexshader", "StandardShading.fragmentshader");
	pickingProgramID = LoadShaders("Picking.vertexshader", "Picking.fragmentshader");

	// Get a handle for our "MVP" uniform
	MatrixID = glGetUniformLocation(programID, "MVP");
	ViewMatrixID = glGetUniformLocation(programID, "V");
	ModelMatrixID = glGetUniformLocation(programID, "M");
	PickingMatrixID = glGetUniformLocation(pickingProgramID, "MVP");
	// Get a handle for our "pickingColorID" uniform
	pickingColorArrayID = glGetUniformLocation(pickingProgramID, "PickingColorArray");
	pickingColorID = glGetUniformLocation(pickingProgramID, "PickingColor");
	// Get a handle for our "LightPosition" uniform
	LightID = glGetUniformLocation(programID, "LightPosition_worldspace");

	createVAOs(Vertices, Indices, sizeof(Vertices), sizeof(Indices), 0);
	createObjects();

	// ATTN: create VAOs for each of the newly created objects here:
	// createVAOs(<fill this appropriately>);

	createVAOs(Lineloop, LlIndices, sizeof(Lineloop), sizeof(LlIndices), 1);
	createObjects();

	createVAOs(Subdivision, SdIndices, sizeof(Subdivision), sizeof(SdIndices), 2);
	createObjects();

	createVAOs(Controlpoints, CpIndices, sizeof(Controlpoints), sizeof(CpIndices), 3);
	createObjects();

	createVAOs(Bezier, BIndices, sizeof(Bezier), sizeof(BIndices), 4);
	createObjects();

	createVAOs(Point, { 0 }, sizeof(Point), sizeof(unsigned short), 5);
	createObjects();

	unsigned short ind[] = { 0,1 };
	createVAOs(Tangent, ind, sizeof(Tangent), sizeof(ind), 6);
	createObjects();

	createVAOs(Normal, ind, sizeof(Normal), sizeof(ind), 7);
	createObjects();

	createVAOs(Binormal, ind, sizeof(Binormal), sizeof(ind), 8);
	createObjects();

}

void createVAOs(Vertex Vertices[], unsigned short Indices[], size_t BufferSize, size_t IdxBufferSize, int ObjectId) {

	NumVert[ObjectId] = IdxBufferSize / (sizeof GLubyte);

	GLenum ErrorCheckValue = glGetError();
	size_t VertexSize = sizeof(Vertices[0]);
	size_t RgbOffset = sizeof(Vertices[0].XYZW);

	// Create Vertex Array Object
	glGenVertexArrays(1, &VertexArrayId[ObjectId]);
	glBindVertexArray(VertexArrayId[ObjectId]);

	// Create Buffer for vertex data
	glGenBuffers(1, &VertexBufferId[ObjectId]);
	glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[ObjectId]);
	glBufferData(GL_ARRAY_BUFFER, BufferSize, Vertices, GL_STATIC_DRAW);

	// Create Buffer for indices
	glGenBuffers(1, &IndexBufferId[ObjectId]);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, IndexBufferId[ObjectId]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, IdxBufferSize, Indices, GL_STATIC_DRAW);

	// Assign vertex attributes
	glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, VertexSize, 0);
	glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, VertexSize, (GLvoid*)RgbOffset);

	glEnableVertexAttribArray(0);	// position
	glEnableVertexAttribArray(1);	// color

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

static void mouseCallback(GLFWwindow* window, int button, int action, int mods)
{
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS  && dview != true) {
		pickVertex();
	}
}

static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_1 && action == GLFW_PRESS)
	{
		count_1 = (count_1 + 1) % 5;
		subd = true;
	}

	if (key == GLFW_KEY_2 && action == GLFW_PRESS)
	{
		if (bez == false)
			bez = true;
		else //if (bez == true)
			bez = false;
	}

	if (key == GLFW_KEY_4 && action == GLFW_PRESS)
	{
		if (dview == false)
			dview = true;
		else //if (dview == true)
			dview = false;
	}

	if (key == GLFW_KEY_5 && action == GLFW_PRESS)
	{
		if (ydot == false)
			ydot = true;
		else
			ydot = false;
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
		// Measure speed
		double currentTime = glfwGetTime();
		nbFrames++;
		if (currentTime - lastTime >= 1.0) { // If last prinf() was more than 1sec ago
											 // printf and reset
			printf("%f ms/frame\n", 1000.0 / double(nbFrames));
			nbFrames = 0;
			lastTime += 1.0;
		}

		// DRAGGING: move current (picked) vertex with cursor
		if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) && dview != true)
			moveVertex();

		glfwSetKeyCallback(window, keyCallback);
		// DRAWING SCENE
		createObjects();	// re-evaluate curves in case vertices have been moved
		drawScene();


	} // Check if the ESC key was pressed or the window was closed
	while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
		glfwWindowShouldClose(window) == 0);

	cleanup();

	return 0;
}
