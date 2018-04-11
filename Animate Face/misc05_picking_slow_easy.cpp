#ifndef TGA_H
#define TGA_H

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
//#include <GL/glut.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

	typedef struct {
		long           width;
		long           height;
		unsigned char *data;
		int            alpha;
	} tTGA;

	int    load_TGA(tTGA *tga, const char *filename);
	void   free_TGA(tTGA *tga);

	//GLuint load_texture_TGA(const char *filename, long *width, long *height, GLint wrap_s, GLint wrap_t);


#ifdef __cplusplus
}
#endif

#endif

#define then
#ifndef TRUE
#define TRUE  (1)
#define FALSE (0)
#endif

// Include standard headers
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <array>
#include <stack>
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

#include <iostream>
using namespace std;
#include <fstream>
#include <iterator>


int load_TGA(tTGA *tga, const char *filename) {

#define SIZEOF_TGA_HEADER 18

	unsigned char   buffer[256];

	int             size_of_image_id;
	int             is_colormap;
	int             targa_type;
	int             colormap_origin;
	unsigned int    colormap_length;
	int             colormap_entry_size;
	int             image_pixel_size;
	int             image_descriptor;
	int             is_inverted;

	int             image_width;
	int             image_height;

	unsigned char  *colormap;
	FILE           *f;
	unsigned char  *data;
	int             x, y, i, j;
	int             raster_width;

	if ((f = fopen(filename, "rb")) == NULL) return FALSE;

	/* header info */
	if (fread(buffer, 1, SIZEOF_TGA_HEADER, f) != SIZEOF_TGA_HEADER) return FALSE;

	size_of_image_id = buffer[0];
	is_colormap = buffer[1];
	targa_type = buffer[2];

	colormap_origin = buffer[3] + ((int)(buffer[4]) << 8);
	colormap_length = buffer[5] + ((int)(buffer[6]) << 8);
	colormap_entry_size = buffer[7];

	image_width = buffer[12] + ((unsigned)(buffer[13]) << 8);
	image_height = buffer[14] + ((unsigned)(buffer[15]) << 8);
	image_pixel_size = buffer[16];
	image_descriptor = buffer[17];

	/* valid type? */
	if ((targa_type != 1) && (targa_type != 2)) return FALSE;

	/* colormap required but missing? */
	if ((targa_type == 1) && !is_colormap) return FALSE;

	/* cannot load direct-color images */
	if ((targa_type == 2) && is_colormap) return FALSE;

	/* image id */
	if (size_of_image_id)
		if ((int)fread(buffer, 1, size_of_image_id, f) != size_of_image_id)
			return FALSE;

	is_inverted = (image_descriptor & 0x10) != 0;

	/* cannot handle interlacing */
	if ((image_descriptor & 0xC0))
		return FALSE;

	/* assume that targa 32 contains alpha (image_descriptor bits 0..3) */

	/* load colormap, if any */
	if (is_colormap)
	{
		/* must be targa 24 or targa 32 */
		if ((colormap_entry_size != 24) && (colormap_entry_size != 32)) return FALSE;

		/* convert to number of bytes/color entry */
		colormap_entry_size >>= 3;

		colormap = (unsigned char*)malloc(colormap_length *colormap_entry_size);
		if (colormap == NULL) return FALSE;

		if (fread(colormap, colormap_entry_size, colormap_length, f) != colormap_length)
		{
		lerror:
			free(colormap);
			return FALSE;
		}

		/* initializations */
		image_pixel_size = (image_pixel_size + 7) >> 3;
		raster_width = image_width *colormap_entry_size;
	}
	else {
		/* must be targa 24 or targa 32 */
		if ((image_pixel_size != 24) && (image_pixel_size != 32)) return FALSE;
		image_pixel_size >>= 3;
		raster_width = image_width *image_pixel_size;
	}

	data = (unsigned char*)malloc(raster_width *image_height);
	if (data == NULL)
		goto lerror;

	/* load image data */
	for (y = (is_inverted ? (image_height - 1) : 0);
		(is_inverted ? (y >= 0) : (y < (int)image_height));
		(is_inverted ? (--y) : (++y)))
		for (x = 0; x < image_width; x++) {

			/* get the next pixel */
			if ((int)fread(buffer, 1, image_pixel_size, f) != image_pixel_size)
				goto lerror;

			/* store it */
			if (is_colormap)
			{
				/* colormapped */
				i = ((buffer[0] + ((unsigned)(buffer[1]) << 8)) - colormap_origin)
					*colormap_entry_size;
				j = (y *raster_width) + (x *colormap_entry_size);

				data[j] = colormap[i + 2];
				data[j + 1] = colormap[i + 1];
				data[j + 2] = colormap[i];

				if (colormap_entry_size > 3)
					data[j + 3] = colormap[i + 3];
			}
			else {
				/* non-colormapped */
				j = (y *raster_width) + (x *image_pixel_size);

				data[j] = buffer[2];
				data[j + 1] = buffer[1];
				data[j + 2] = buffer[0];

				if (image_pixel_size > 3)
					data[j + 3] = buffer[3];
			}
		}

	/* free the colormap if we had loaded it */
	if (is_colormap)
		free(colormap);

	/* store the result */
	tga->width = image_width;
	tga->height = image_height;
	tga->data = data;
	tga->alpha = (is_colormap ? (colormap_entry_size > 3) : (image_pixel_size > 3));

#undef SIZEOF_TGA_HEADER

	return TRUE;
}

/*--------------------------------------------------------------------------+/
free_TGA
/+--------------------------------------------------------------------------*/
void free_TGA(tTGA *tga) {

	if (tga->data)
		free(tga->data);

	tga->data = NULL;
	tga->height =
		tga->width = 0;
	tga->alpha = 0;
}

/*--------------------------------------------------------------------------+/
load_texture_TGA
/+--------------------------------------------------------------------------*/
GLuint load_texture_TGA(const char *filename, long *width, long *height, GLint wrap_s, GLint wrap_t) {

	GLuint result;
	tTGA   tga;

	glEnable(GL_TEXTURE_2D);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	glGenTextures(1, &result);

	if (!load_TGA(&tga, filename))
	{
		glDeleteTextures(1, &result);
		return 0;
	}


	glBindTexture(GL_TEXTURE_2D, result);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, wrap_s);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, wrap_t);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	glTexImage2D(
		GL_TEXTURE_2D,
		0,
		(tga.alpha) ? (GL_RGBA) : (GL_RGB),
		tga.width,
		tga.height,
		0,
		(tga.alpha) ? (GL_RGBA) : (GL_RGB),
		GL_UNSIGNED_BYTE,
		tga.data
	);

	if (width)
		*width = tga.width;
	if (height)
		*height = tga.height;

	free_TGA(&tga);

	return result;
}

const int window_width = 600, window_height = 600;

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

typedef struct point {
	float x, y, z;
	point(const float x = 0, const float y = 0, const float z = 0) : x(x), y(y), z(z) {}
	point(float *coords) : x(coords[0]), y(coords[1]), z(coords[2]) {}
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
void loadObject(char*, glm::vec4, Vertex * &, GLushort* &, int);
void createVAOs(Vertex[], GLushort[], int);
void createObjects(void);
void pickObject(void);
void renderScene(void);
void cleanup(void);
static void keyCallback(GLFWwindow*, int, int, int, int);
static void mouseCallback(GLFWwindow*, int, int, int);
static void scrollCallback(GLFWwindow*, double, double);

void createVAOsForTex(Vertex[], GLushort[], int);
void save(void);
void loadControlPoints(void);
vec3 findClosestPoint(vec3, vec3, vec3, double);
bool rayTest(vec3, vec3, vec3, vec3, double, double);
bool rayTestPoints(Vertex*, vec3, vec3, unsigned int*, double*, double);
void subd(void);
void move_vertex(void);

// GLOBAL VARIABLES
GLFWwindow* window;

glm::mat4 gProjectionMatrix;
glm::mat4 gViewMatrix;

GLuint gPickedIndex = -1;
std::string gMessage;

GLuint programID;
GLuint pickingProgramID;
GLuint textureProgramID;

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
GLuint TextureID;

GLint gX = 0.0;
GLint gZ = 0.0;

// animation control
bool animation = false;
bool smile = false;
GLfloat phi = 0.0;
float theta = 55;		// Angle from Y axis
float Phi = 315;			// Angle from Z axis
float radius = 40.0f;
float cX = 10.0f;
float cY = 10.0f;
float cZ = 10.0f;


Vertex* faceVert;
GLushort* faceIdcs;
Vertex* hairVerts;
GLushort* hairIdcs;
Vertex meshVerts[441];
GLushort meshIdcs[1764];
GLushort texIdcs[2646];
Vertex subdVerts[3721] = { 0.0f };
GLushort subdIdcs[14884];
GLushort subdTexIdcs[22326];


long image_width;
long image_height;
GLuint texID;
GLfloat uv[882];
GLfloat uvSubdiv[7442];
GLint viewport[4];
vec3 startMousePos;
vec3 endMousePos;
unsigned int id;
double proj;
float colorRed[] = { 1.0f, 0.0f, 0.0f, 1.0f };
float meshNormal[] = { 0.0f, 0.0f, 1.0f };
glm::mat4 ModelMatrix;
 
bool camera = true;
bool face = false;
bool controlPts = false;
bool tex = false;
bool subdPts = false;

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

void createObjects(void)
{
	//-- COORDINATE AXES --//
	
	//-- GRID --//

	// ATTN: create your grid vertices here!
	Vertex gridVert[120];
	float white[4] = { 1.0,1.0,1.0,1.0 };
	float red[4] = { 0.0,1.0,0.0,1.0 };
	float zeronm[3] = { 0.0, 0.0, 1.0 };
	int gridIndex = 0;
	float color[] = { 0.9f,0.9f,0.9f,1.0f };
	for (float x = -10.0f; x < 11.0f; x ++)
	{
		gridVert[gridIndex].Position[0] = x;
		gridVert[gridIndex].Position[1] = 0.0f;
		gridVert[gridIndex].Position[2] = -10.0f;
		gridVert[gridIndex].Position[3] = 1.0f;

		//gridIndices[gridIndex] = gridIndex;

		gridVert[gridIndex++].SetColor(color);

		gridVert[gridIndex].Position[0] = x;
		gridVert[gridIndex].Position[1] = 0.0f;
		gridVert[gridIndex].Position[2] = 10.0f;
		gridVert[gridIndex].Position[3] = 1.0f;

	//	gridIndices[gridIndex] = gridIndex;

		gridVert[gridIndex++].SetColor(color);
	}

	for (float x = -10.0f; x < 11.0f; x ++)
	{
		gridVert[gridIndex].Position[0] = -10.0f;
		gridVert[gridIndex].Position[1] = 0.0f;
		gridVert[gridIndex].Position[2] = x;
		gridVert[gridIndex].Position[3] = 1.0f;

	//	gridIndices[gridIndex] = gridIndex;

		gridVert[gridIndex++].SetColor(color);

		gridVert[gridIndex].Position[0] = 10.0f;
		gridVert[gridIndex].Position[1] = 0.0f;
		gridVert[gridIndex].Position[2] = x;
		gridVert[gridIndex].Position[3] = 1.0f;

	//	gridIndices[gridIndex] = gridIndex;

		gridVert[gridIndex++].SetColor(color);
	}

	VertexBufferSize[1] = sizeof(gridVert);
	createVAOs(gridVert, NULL, 1);

	int k = 0;
	for (int i = -10; i <= 10; i++) {
		for (int j = 0; j <= 20; j++) {
			meshVerts[21 * k + j].Position[0] = (float)i;
			meshVerts[21 * k + j].Position[1] = j;
			meshVerts[21 * k + j].Position[2] = -10;
			meshVerts[21 * k + j].Position[3] = 1.0;
			meshVerts[21 * k + j].SetColor(red);
			meshVerts[21 * k + j].SetNormal(zeronm);

		}
	
		k++;
	}

	texID = load_texture_TGA("Nagavara.tga", &image_width, &image_height, GL_CLAMP, GL_CLAMP);
	
	for (int i = 0; i < 441; i++) {
				
		uv[2 * i] = (meshVerts[i].Position[0] + 10) / 20;
		uv[2 * i + 1] = (meshVerts[i].Position[1]) / 20;

		if ((i + 1) % 21 != 0 && i < 420 && i != 440) {
			meshIdcs[4 * i] = i;
			meshIdcs[4 * i + 1] = i + 1;
		}
		else {
			meshIdcs[4 * i] = i;
			meshIdcs[4 * i + 1] = i;
		}
		if (i < 420) {
			meshIdcs[4 * i + 2] = i;
			meshIdcs[4 * i + 3] = i + 21;
		}
		else if (i != 440) {
			meshIdcs[4 * i + 2] = i;
			meshIdcs[4 * i + 3] = i + 1;
		}
		if (i == 0 || (i + 1) % 21 != 0 && i < 420) {
			texIdcs[6 * i] = i;
			texIdcs[6 * i + 1] = i + 1;
			texIdcs[6 * i + 2] = i + 22;
			texIdcs[6 * i + 3] = i + 22;
			texIdcs[6 * i + 4] = i + 21;
			texIdcs[6 * i + 5] = i;
		}
	}

	VertexBufferSize[2] = sizeof(meshVerts);
	IndexBufferSize[2] = sizeof(meshIdcs);
	createVAOs(meshVerts, meshIdcs, 2);

	VertexBufferSize[3] = sizeof(meshVerts);
	IndexBufferSize[3] = sizeof(texIdcs);
	createVAOsForTex(meshVerts, texIdcs, 3);

	//-- .OBJs --//

	// ATTN: load your models here
	Vertex* hairVerts;
	GLushort* hairIdcs;
	
	loadObject("Nagavara.obj", glm::vec4(1.0, 1.0, 1.0, 1.0), faceVert, faceIdcs, 4);
	createVAOs(faceVert, faceIdcs, 4);
	
	//loadObject("hair.obj", glm::vec4(1.0, 1.0, 1.0, 1.0), hairVerts, hairIdcs, 7);
	//createVAOs(hairVerts, hairIdcs, 7);


}
int t = 0;
void renderScene(void)
{
	//ATTN: DRAW YOUR SCENE HERE. MODIFY/ADAPT WHERE NECESSARY!

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

	// Dark blue background
	glClearColor(0.0f, 0.0f, 0.2f, 0.0f);
	// Re-clear the screen for real rendering
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	
	glUseProgram(programID);
	{
		glm::vec3 lightPos = glm::vec3(5, 30, 10);
		glm::vec3 lightPos1 = glm::vec3((cX + 0.5f), cY, cZ);
		glm::mat4x4 ModelMatrix = glm::mat4(1.0);

		//glUniform3f(LightID, 3,3,3);
		glUniform3f(LightID, lightPos.x, lightPos.y, lightPos.z);
		//glUniform3f(LightID1, 2, 4, 4);
		glUniform3f(LightID1, lightPos1.x, lightPos1.y, lightPos1.z);

		glUniformMatrix4fv(ViewMatrixID, 1, GL_FALSE, &gViewMatrix[0][0]);
		glUniformMatrix4fv(ProjMatrixID, 1, GL_FALSE, &gProjectionMatrix[0][0]);
		glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);

		glBindVertexArray(VertexArrayId[0]);	// draw CoordAxes
		glDrawArrays(GL_LINES, 0, 6);

		glBindVertexArray(VertexArrayId[1]);
		glDrawArrays(GL_LINES, 0, 120);

		if (controlPts) {
			glDisable(GL_PROGRAM_POINT_SIZE);
			glEnable(GL_POINT_SIZE);
			glPointSize(7.0f);
			glBindVertexArray(VertexArrayId[2]);
			glDrawArrays(GL_POINTS, 0, 441);
			//glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
			glDrawElements(GL_LINES, 1764, GL_UNSIGNED_SHORT, (GLvoid*)0);
		}

		if (face) {
			glBindVertexArray(VertexArrayId[4]);
			glDrawElements(GL_TRIANGLES, NumIndices[4], GL_UNSIGNED_SHORT, (void*)0);

			//glBindVertexArray(VertexArrayId[7]);
			//glDrawElements(GL_TRIANGLES, NumIndices[7], GL_UNSIGNED_SHORT, (void*)0);
		}

		if (subdPts) {
			glDisable(GL_PROGRAM_POINT_SIZE);
			glEnable(GL_POINT_SIZE);
			glPointSize(4.0f);
			glBindVertexArray(VertexArrayId[5]);
			glDrawArrays(GL_POINTS, 0, 3721);
			glDrawElements(GL_LINES, 14884, GL_UNSIGNED_SHORT, (GLvoid*)0);
		}


		if (smile) {

			_sleep(700);
			if (t < 5)
			{				
				meshVerts[196].Position[1] += 0.1;			
				meshVerts[259].Position[1] += 0.1;
				t++;

				VertexBufferSize[2] = sizeof(meshVerts);
				IndexBufferSize[2] = sizeof(meshIdcs);
				createVAOs(meshVerts, meshIdcs, 2);

				VertexBufferSize[3] = sizeof(meshVerts);
				IndexBufferSize[3] = sizeof(texIdcs);
				createVAOsForTex(meshVerts, texIdcs, 3);

			}
			else {
				t = 0;
				smile = false;
			}
		}


		glBindVertexArray(0);
	}

	if (controlPts || subdPts) {
		glUseProgram(textureProgramID); {
			glm::vec3 lightPos = glm::vec3(5, 30, 10);
			glm::vec3 lightPos1 = glm::vec3((cX + 0.5f), cY, cZ);
			glm::mat4x4 ModelMatrix = glm::mat4(1.0);

			//glUniform3f(LightID, 3,3,3);
			glUniform3f(LightID, lightPos.x, lightPos.y, lightPos.z);
			//glUniform3f(LightID1, 2, 4, 4);
			glUniform3f(LightID1, lightPos1.x, lightPos1.y, lightPos1.z);

			glUniformMatrix4fv(ViewMatrixID, 1, GL_FALSE, &gViewMatrix[0][0]);
			glUniformMatrix4fv(ProjMatrixID, 1, GL_FALSE, &gProjectionMatrix[0][0]);
			glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, texID);
			glUniform1i(TextureID, 0);
			
			glBindVertexArray(VertexArrayId[3]);
			glDrawElements(GL_TRIANGLES, 2646, GL_UNSIGNED_SHORT, (GLvoid*)0);

			glBindVertexArray(0);
		}
	}

	glUseProgram(0);
	// Draw GUI
	//TwDraw();

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
		ModelMatrix = glm::mat4(1.0); // TranslationMatrix * RotationMatrix;
		glm::mat4 MVP = gProjectionMatrix * gViewMatrix * ModelMatrix;

		// Send our transformation to the currently bound shader, in the "MVP" uniform
		glUniformMatrix4fv(PickingMatrixID, 1, GL_FALSE, &MVP[0][0]);

		// ATTN: DRAW YOUR PICKING SCENE HERE. REMEMBER TO SEND IN A DIFFERENT PICKING COLOR FOR EACH OBJECT BEFOREHAND
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

	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	mat4 nModelMatrix = gViewMatrix * ModelMatrix;

	unsigned char data[4];
	GLfloat zpos;

	glReadPixels(xpos, window_height - ypos, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, data); // OpenGL renders with (0,0) on bottom, mouse reports with (0,0) on top

	startMousePos = glm::unProject(glm::vec3(xpos, window_height - ypos, 0.0f), nModelMatrix, gProjectionMatrix, vec4(viewport[0], viewport[1], viewport[2], viewport[3]));
	endMousePos = glm::unProject(glm::vec3(xpos, window_height - ypos, 1.0f), nModelMatrix, gProjectionMatrix, vec4(viewport[0], viewport[1], viewport[2], viewport[3]));
	double epsilon = 0.1;
	double proj;
	bool found = rayTestPoints(meshVerts, startMousePos, endMousePos, &id, &proj, epsilon);
	
	// Convert the color back to an integer ID
	gPickedIndex = int(data[0]);

	if (gPickedIndex == 255) { // Full white, must be the background !
		gMessage = "background";
	}
	else {
		std::ostringstream oss;
		oss << "point " << gPickedIndex;
		gMessage = oss.str();
	}

	// Uncomment these lines to see the picking shader in effect
	//glfwSwapBuffers(window);
	//continue; // skips the normal rendering
}

void move_vertex() {

	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	glm::vec4 vp = glm::vec4(viewport[0], viewport[1], viewport[2], viewport[3]);

	double xpos, ypos;
	glfwGetCursorPos(window, &xpos, &ypos);
	mat4 nModelMatrix = gViewMatrix * ModelMatrix;


	int state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
	vec3 W;
	if (state == GLFW_PRESS) {		
		if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS) {
			vec3 p = glm::project(glm::vec3(meshVerts[id].Position[0], meshVerts[id].Position[1], meshVerts[id].Position[2]), nModelMatrix, gProjectionMatrix, glm::vec4(viewport[0], viewport[1], viewport[2], viewport[3]));
			W = glm::unProject(glm::vec3(xpos, window_height - ypos, p.z), nModelMatrix, gProjectionMatrix, glm::vec4(viewport[0], viewport[1], viewport[2], viewport[3]));
			float coords[3] = { W.x, W.y, W.z };
			meshVerts[id].SetPosition(coords);
		}
		else {
			vec3 p = glm::project(glm::vec3(meshVerts[id].Position[0], meshVerts[id].Position[1], meshVerts[id].Position[2]), nModelMatrix, gProjectionMatrix, glm::vec4(viewport[0], viewport[1], viewport[2], viewport[3]));
			W = glm::unProject(glm::vec3(xpos, window_height - ypos, p.z), nModelMatrix, gProjectionMatrix, glm::vec4(viewport[0], viewport[1], viewport[2], viewport[3]));
			float coords[3] = { W.x, W.y, meshVerts[id].Position[2] };
			meshVerts[id].SetPosition(coords);
		}
	}

	VertexBufferSize[2] = sizeof(meshVerts);
	IndexBufferSize[2] = sizeof(meshIdcs);
	createVAOs(meshVerts, meshIdcs, 2);

	VertexBufferSize[3] = sizeof(meshVerts);
	IndexBufferSize[3] = sizeof(texIdcs);
	createVAOsForTex(meshVerts, texIdcs, 3);
}

vec3 findClosestPoint(vec3 rayStartPos, vec3 rayEndPos, vec3 pointPos, double *proj) {
	vec3 rayVector = rayEndPos - rayStartPos;
	double raySquared = glm::dot(rayVector, rayVector);
	vec3 projection = pointPos - rayStartPos;
	double projectionVal = glm::dot(projection, rayVector);
	*proj = projectionVal / raySquared;
	vec3 closestPoint = rayStartPos + glm::vec3(rayVector.x * (*proj), rayVector.y * (*proj), rayVector.z * (*proj));
	return closestPoint;
}

bool rayTest(vec3 pointPos, vec3 startPos, vec3 endPos, vec3 *closestPoint, double *proj, double epsilon) {
	*closestPoint = findClosestPoint(startPos, endPos, pointPos, proj);
	double len = glm::distance2(*closestPoint, pointPos);	
	return len < epsilon;
}

bool rayTestPoints(Vertex* vert, vec3 start, vec3 end, unsigned int *id, double *proj, double epsilon) {
	unsigned int pointID = 442;
	bool foundCollision = false;
	double minDistToStart = 10000000.0;
	double distance;
	vec3 point;
	for (unsigned int i = 0; i < 441; ++i) {		
		vec3 pointPos = glm::vec3(vert[i].Position[0], vert[i].Position[1], vert[i].Position[2]);
		if (rayTest(pointPos, start, end, &point, proj, epsilon)) {
			distance = glm::distance2(start, point);		
			if (distance < minDistToStart)
			{
				minDistToStart = distance;
				pointID = i;
				foundCollision = true;
			}
		}
	}

	*id = pointID;
	return foundCollision;
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
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

	// Open a window and create its OpenGL context
	window = glfwCreateWindow(window_width, window_height, "Janakiprasad Nagavara,Adarsh (66975994)", NULL, NULL);
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
	glfwSetScrollCallback(window, scrollCallback);

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

	// Projection matrix : 45∞ Field of View, 4:3 ratio, display range : 0.1 unit <-> 100 units
	gProjectionMatrix = glm::perspective(45.0f, 4.0f / 3.0f, 0.1f, 100.0f);
	// Or, for an ortho camera :
	//gProjectionMatrix = glm::ortho(-4.0f, 4.0f, -3.0f, 3.0f, 0.0f, 100.0f); // In world coordinates

	// Camera matrix
	gViewMatrix = glm::lookAt(glm::vec3(15.0, 15.0, 15.0f),	// eye
		glm::vec3(0.0, 10.0, 0.0),	// center
		glm::vec3(0.0, 1.0, 0.0));	// up

									// Create and compile our GLSL program from the shaders
	programID = LoadShaders("StandardShading.vertexshader", "StandardShading.fragmentshader");
	pickingProgramID = LoadShaders("Picking.vertexshader", "Picking.fragmentshader");
	textureProgramID = LoadShaders("TextureVertexShader.vertexshader", "TextureFragmentShader.fragmentshader");

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
	TextureID = glGetUniformLocation(textureProgramID, "myTextureSampler");

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

void createVAOsForTex(Vertex Vertices[], unsigned short Indices[], int ObjectId) {
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
	glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, VertexSize, (GLvoid*)Normaloffset);

	glEnableVertexAttribArray(0);	// position
									//glEnableVertexAttribArray(2);	// normal

	if (ObjectId == 6) {		
		GLuint uvbuffer;
		glGenBuffers(1, &uvbuffer);
		glBindBuffer(GL_ARRAY_BUFFER, uvbuffer);
		glBufferData(GL_ARRAY_BUFFER, sizeof(uvSubdiv), uvSubdiv, GL_STATIC_DRAW);
		glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, (GLvoid*)0);
		glEnableVertexAttribArray(1);	// color
	}
	else {		
		GLuint uvbuffer;
		glGenBuffers(1, &uvbuffer);
		glBindBuffer(GL_ARRAY_BUFFER, uvbuffer);
		glBufferData(GL_ARRAY_BUFFER, sizeof(uv), uv, GL_STATIC_DRAW);
		glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, (GLvoid*)0);
		glEnableVertexAttribArray(1);	// color
	}

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
		switch (key) {
		case GLFW_KEY_A:
			subdPts = true;
			subd();
			break;
		case GLFW_KEY_C:
			if (controlPts) {
				controlPts = false;
			}
			else {
				controlPts = true;
			}
			break;
		case GLFW_KEY_F:
			face = !face;			
			break;
		case GLFW_KEY_L:
			loadControlPoints();
			break;
		case GLFW_KEY_R:
			gMessage = "Reset Camera";
			theta = 55;		// Angle from Y axis
			Phi = 315;			// Angle from Z axis
			radius = 40.0f;		//radius = 7.0f;			
			face = true;

			break;
		case GLFW_KEY_S:
			if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
				glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS)
			{
				smile = true;
			}
			else
			{
				save();
			}
			
			break;		
		case GLFW_KEY_DOWN:
		{

			if (camera == true)
			{
				theta += 5.0;
				renderScene();
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

			break;
		}

		case GLFW_KEY_RIGHT:
		{
			if (camera == true)
			{
				Phi -= 5.0;
				renderScene();
				break;
			}

			break;
		}
		case GLFW_KEY_LEFT:
		{

			if (camera == true)
			{
				Phi += 5.0;
				renderScene();
				break;
			}
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


void save() {
	cout << "Writing Control Point file ..." << endl;
	ofstream controlPointFile;
	controlPointFile.open("cm.p3", ios::out);
	if (controlPointFile.is_open()) {
		for (int i = 0; i < 441; i++) {
			controlPointFile << "v ";
			controlPointFile << to_string(meshVerts[i].Position[0]) + " ";
			controlPointFile << to_string(meshVerts[i].Position[1]) + " ";
			controlPointFile << to_string(meshVerts[i].Position[2]) + " ";
			controlPointFile << to_string(meshVerts[i].Position[3]) << endl;
		}

		for (int i = 0; i < 441; i++) {
			controlPointFile << "f ";
			controlPointFile << to_string(texIdcs[6 * i]) + " ";
			controlPointFile << to_string(texIdcs[6 * i + 1]) + " ";
			controlPointFile << to_string(texIdcs[6 * i + 2]) + " ";
			controlPointFile << to_string(texIdcs[6 * i + 3]) + " ";
			controlPointFile << to_string(texIdcs[6 * i + 4]) + " ";
			controlPointFile << to_string(texIdcs[6 * i + 5]) << endl;
		}
		cout << "Control points written to cm.p3" << endl;
	}
	else {
		cout << "Unable to open control point file." << endl;
	}
	controlPointFile.close();
}

void loadControlPoints() {
	cout << "Loading Control Point file ..." << endl;
	ifstream controlPointFile;
	controlPointFile.open("cm.p3", ios::in);
	//controlPointFile.open("cm.p3", ios::in);
	if (controlPointFile.is_open()) {
		int vert_idx = 0;
		int idcs_idx = 0;
		while (controlPointFile) {
			string line;
			getline(controlPointFile, line);
			if (!line.empty()) {
				istringstream input_string_stream(line);
				vector<string> tokens{istream_iterator<string>{input_string_stream},istream_iterator<string>{}};
				if (tokens.at(0).compare("v") == 0) {
					meshVerts[vert_idx].Position[0] = atof(tokens.at(1).c_str());
					meshVerts[vert_idx].Position[1] = atof(tokens.at(2).c_str());
					meshVerts[vert_idx].Position[2] = atof(tokens.at(3).c_str());
					meshVerts[vert_idx].Position[3] = atof(tokens.at(4).c_str());
					vert_idx++;
				}

				if (tokens.at(0).compare("f") == 0) {
					texIdcs[6 * idcs_idx] = atoi(tokens.at(1).c_str());
					texIdcs[6 * idcs_idx + 1] = atoi(tokens.at(2).c_str());
					texIdcs[6 * idcs_idx + 2] = atoi(tokens.at(3).c_str());
					texIdcs[6 * idcs_idx + 3] = atoi(tokens.at(4).c_str());
					texIdcs[6 * idcs_idx + 4] = atoi(tokens.at(5).c_str());
					texIdcs[6 * idcs_idx + 5] = atoi(tokens.at(6).c_str());
					idcs_idx++;
				}
			}
		}
		cout << "Control Points loaded from cm.p3" << endl;
	}
	else {
		cout << "Unable to open control point file." << endl;
	}
	controlPointFile.close();
	VertexBufferSize[2] = sizeof(meshVerts);
	IndexBufferSize[2] = sizeof(meshIdcs);
	createVAOs(meshVerts, meshIdcs, 2);

	VertexBufferSize[3] = sizeof(meshVerts);
	IndexBufferSize[3] = sizeof(texIdcs);
	createVAOsForTex(meshVerts, texIdcs, 3);
}

static void scrollCallback(GLFWwindow * window, double xoffset, double yoffset)
{
	if (yoffset <  15 && yoffset > -10)
	{
		radius = radius - yoffset;
		renderScene();
	}
}


void subd() {
	int j = 0;
	for (int i = 0; i < 441; i++) {		
		point *s00, *s01, *s02, *s10, *s11, *s12, *s20, *s21, *s22;
		float u00, v00, u01, v01, u02, v02, u10, v10, u11, v11, u12, v12, u20, v20, u21, v21, u22, v22;
		if (i < 21) {			
			s00 = new point(meshVerts[i].Position[0] - 1, meshVerts[i].Position[1] - 1, meshVerts[i].Position[2]);
			u00 = uv[2 * i] - 1;
			v00 = uv[2 * i + 1] - 1;
			s01 = new point(meshVerts[i].Position[0] - 1, meshVerts[i].Position[1], meshVerts[i].Position[2]);
			u01 = uv[2 * i] - 1;
			v01 = uv[2 * i + 1];
			s02 = new point(meshVerts[i].Position[0] - 1, meshVerts[i].Position[1] + 1, meshVerts[i].Position[2]);
			u02 = uv[2 * i] - 1;
			v02 = uv[2 * i + 1] + 1;
			if (i == 0) {
				s10 = new point(meshVerts[i].Position[0], meshVerts[i].Position[1] - 1, meshVerts[i].Position[2]);
				u10 = uv[2 * i];
				v10 = uv[2 * i + 1] - 1;
			}
			else {
				s10 = new point(meshVerts[i - 1].Position[0], meshVerts[i - 1].Position[1], meshVerts[i - 1].Position[2]);
				u10 = uv[2 * (i - 1)];
				v10 = uv[2 * (i - 1) + 1];
			}
			s11 = new point(meshVerts[i].Position[0], meshVerts[i].Position[1], meshVerts[i].Position[2]);
			u11 = uv[2 * i];
			v11 = uv[2 * i + 1];
			if (i == 20) {
				s12 = new point(meshVerts[i].Position[0], meshVerts[i].Position[1] + 1, meshVerts[i].Position[2]);
				u12 = uv[2 * i];
				v12 = uv[2 * i + 1] + 1;
			}
			else {
				s12 = new point(meshVerts[i + 1].Position[0], meshVerts[i + 1].Position[1], meshVerts[i + 1].Position[2]);
				u12 = uv[2 * (i + 1)];
				v12 = uv[2 * (i + 1) + 1];
			}
			if (i == 0) {
				s20 = new point(meshVerts[i].Position[0] + 1, meshVerts[i].Position[1] - 1, meshVerts[i].Position[2]);
				u20 = uv[2 * i] + 1;
				v20 = uv[2 * i + 1] - 1;
			}
			else {
				s20 = new point(meshVerts[i + 20].Position[0], meshVerts[i + 20].Position[1], meshVerts[i + 20].Position[2]);
				u20 = uv[2 * (i + 20)];
				v20 = uv[2 * (i + 20) + 1];
			}
			s21 = new point(meshVerts[i + 21].Position[0], meshVerts[i + 21].Position[1], meshVerts[i + 21].Position[2]);
			u21 = uv[2 * (i + 21)];
			v21 = uv[2 * (i + 21) + 1];
			if (i == 20) {
				s22 = new point(meshVerts[i].Position[0] + 1, meshVerts[i].Position[1] + 1, meshVerts[i].Position[2]);
				u22 = uv[2 * i] + 1;
				v22 = uv[2 * i + 1] + 1;
			}
			else {
				s22 = new point(meshVerts[i + 22].Position[0], meshVerts[i + 22].Position[1], meshVerts[i + 22].Position[2]);
				u22 = uv[2 * (i + 22)];
				v22 = uv[2 * (i + 22) + 1];
			}

		}
		else if ((i + 1) % 21 == 0 && i > 21 && i < 420) {			
			s00 = new point(meshVerts[i - 22].Position[0], meshVerts[i - 22].Position[1], meshVerts[i - 22].Position[2]);
			u00 = uv[2 * (i - 22)];
			v00 = uv[2 * (i - 22) + 1];
			s01 = new point(meshVerts[i - 21].Position[0], meshVerts[i - 21].Position[1], meshVerts[i - 21].Position[2]);
			u01 = uv[2 * (i - 21)];
			v01 = uv[2 * (i - 21) + 1];
			s02 = new point(meshVerts[i].Position[0] - 1, meshVerts[i].Position[1] + 1, meshVerts[i].Position[2]);
			u02 = uv[2 * i] - 1;
			v02 = uv[2 * i + 1] - 1;
			s10 = new point(meshVerts[i - 1].Position[0], meshVerts[i - 1].Position[1], meshVerts[i - 1].Position[2]);
			u10 = uv[2 * (i - 1)];
			v10 = uv[2 * (i - 1) + 1];
			s11 = new point(meshVerts[i].Position[0], meshVerts[i].Position[1], meshVerts[i].Position[2]);
			u11 = uv[2 * i];
			v11 = uv[2 * i + 1];
			s12 = new point(meshVerts[i].Position[0], meshVerts[i].Position[1] + 1, meshVerts[i].Position[2]);
			u12 = uv[2 * i];
			v12 = uv[2 * i + 1] + 1;
			s20 = new point(meshVerts[i + 20].Position[0], meshVerts[i + 20].Position[1], meshVerts[i + 20].Position[2]);
			u20 = uv[2 * (i + 20)];
			v20 = uv[2 * (i + 20) + 1];
			s21 = new point(meshVerts[i + 21].Position[0], meshVerts[i + 21].Position[1], meshVerts[i + 21].Position[2]);
			u21 = uv[2 * (i + 21)];
			v21 = uv[2 * (i + 21) + 1];
			s22 = new point(meshVerts[i].Position[0] + 1, meshVerts[i].Position[1] + 1, meshVerts[i].Position[2]);
			u22 = uv[2 * i] + 1;
			v22 = uv[2 * i + 1] + 1;
		}
		else if (i > 420 && ((i + 1) % 21 != 0 || i == 440)) {			
			s00 = new point(meshVerts[i - 22].Position[0], meshVerts[i - 22].Position[1], meshVerts[i - 22].Position[2]);
			u00 = uv[2 * (i - 22)];
			v00 = uv[2 * (i - 22) + 1];
			s01 = new point(meshVerts[i - 21].Position[0], meshVerts[i - 21].Position[1], meshVerts[i - 21].Position[2]);
			u01 = uv[2 * (i - 21)];
			v01 = uv[2 * (i - 21) + 1];
			if (i == 440) {
				s02 = new point(meshVerts[i].Position[0] - 1, meshVerts[i].Position[1] + 1, meshVerts[i].Position[2]);
				u02 = uv[2 * i] - 1;
				v02 = uv[2 * i + 1] + 1;
			}
			else {
				s02 = new point(meshVerts[i - 20].Position[0], meshVerts[i - 20].Position[1], meshVerts[i - 20].Position[2]);
				u02 = uv[2 * (i - 20)];
				v02 = uv[2 * (i - 20) + 1];
			}
			s10 = new point(meshVerts[i - 1].Position[0], meshVerts[i - 1].Position[1], meshVerts[i - 1].Position[2]);
			u10 = uv[2 * (i - 1)];
			v10 = uv[2 * (i - 1) + 1];
			s11 = new point(meshVerts[i].Position[0], meshVerts[i].Position[1], meshVerts[i].Position[2]);
			u11 = uv[2 * i];
			v11 = uv[2 * i + 1];
			if (i == 440) {
				s12 = new point(meshVerts[i].Position[0], meshVerts[i].Position[1] + 1, meshVerts[i].Position[2]);
				u12 = uv[2 * i];
				v12 = uv[2 * i + 1] + 1;
			}
			else {
				s12 = new point(meshVerts[i + 1].Position[0], meshVerts[i + 1].Position[1], meshVerts[i + 1].Position[2]);
				u12 = uv[2 * (i + 1)];
				v12 = uv[2 * (i + 1) + 1];
			}
			s20 = new point(meshVerts[i].Position[0] + 1, meshVerts[i].Position[1] - 1, meshVerts[i].Position[2]);
			u20 = uv[2 * i] + 1;
			v20 = uv[2 * i + 1] - 1;
			s21 = new point(meshVerts[i].Position[0] + 1, meshVerts[i].Position[1], meshVerts[i].Position[2]);
			u21 = uv[2 * i] + 1;
			v21 = uv[2 * i + 1];
			s22 = new point(meshVerts[i].Position[0] + 1, meshVerts[i].Position[1] + 1, meshVerts[i].Position[2]);
			u22 = uv[2 * i] + 1;
			v22 = uv[2 * i + 1] + 1;
		}
		else if (i % 21 == 0 && i >= 21 && i <= 420) {			
			s00 = new point(meshVerts[i].Position[0] - 1, meshVerts[i].Position[1] - 1, meshVerts[i].Position[2]);
			u00 = uv[2 * i] - 1;
			v00 = uv[2 * i + 1] - 1;
			s01 = new point(meshVerts[i - 21].Position[0], meshVerts[i - 21].Position[1], meshVerts[i - 21].Position[2]);
			u01 = uv[2 * (i - 21)];
			v01 = uv[2 * (i - 21) + 1];
			s02 = new point(meshVerts[i - 20].Position[0], meshVerts[i - 20].Position[1], meshVerts[i - 20].Position[2]);
			u02 = uv[2 * (i - 20)];
			v02 = uv[2 * (i - 20) + 1];
			s10 = new point(meshVerts[i].Position[0], meshVerts[i].Position[1] - 1, meshVerts[i].Position[2]);
			u10 = uv[2 * i];
			v10 = uv[2 * i + 1] - 1;
			s11 = new point(meshVerts[i].Position[0], meshVerts[i].Position[1], meshVerts[i].Position[2]);
			u11 = uv[2 * i];
			v11 = uv[2 * i + 1];
			s12 = new point(meshVerts[i + 1].Position[0], meshVerts[i + 1].Position[1], meshVerts[i + 1].Position[2]);
			u12 = uv[2 * (i + 1)];
			v12 = uv[2 * (i + 1) + 1];
			s20 = new point(meshVerts[i].Position[0] + 1, meshVerts[i].Position[1] - 1, meshVerts[i].Position[2]);
			u20 = uv[2 * i] + 1;
			v20 = uv[2 * i + 1] - 1;
			if (i == 420) {
				s21 = new point(meshVerts[i].Position[0] + 1, meshVerts[i].Position[1], meshVerts[i].Position[2]);
				u21 = uv[2 * i] + 1;
				v21 = uv[2 * i + 1];
			}
			else {
				s21 = new point(meshVerts[i + 21].Position[0], meshVerts[i + 21].Position[1], meshVerts[i + 21].Position[2]);
				u21 = uv[2 * (i + 21)];
				v21 = uv[2 * (i + 21) + 1];
			}
			if (i == 420) {
				s22 = new point(meshVerts[i].Position[0] + 1, meshVerts[i].Position[1] + 1, meshVerts[i].Position[2]);
				u22 = uv[2 * i] + 1;
				v22 = uv[2 * i + 1] + 1;
			}
			else {
				s22 = new point(meshVerts[i + 22].Position[0], meshVerts[i + 22].Position[1], meshVerts[i + 22].Position[2]);
				u22 = uv[2 * (i + 22)];
				v22 = uv[2 * (i + 22) + 1];
			}
		}
		else {
			
			s00 = new point(meshVerts[i - 22].Position[0], meshVerts[i - 22].Position[1], meshVerts[i - 22].Position[2]);
			u00 = uv[2 * (i - 22)];
			v00 = uv[2 * (i - 22) + 1];
			s01 = new point(meshVerts[i - 21].Position[0], meshVerts[i - 21].Position[1], meshVerts[i - 21].Position[2]);
			u01 = uv[2 * (i - 21)];
			v01 = uv[2 * (i - 21) + 1];
			s02 = new point(meshVerts[i - 20].Position[0], meshVerts[i - 20].Position[1], meshVerts[i - 20].Position[2]);
			u02 = uv[2 * (i - 20)];
			v02 = uv[2 * (i - 20) + 1];
			s10 = new point(meshVerts[i - 1].Position[0], meshVerts[i - 1].Position[1], meshVerts[i - 1].Position[2]);
			u10 = uv[2 * (i - 1)];
			v10 = uv[2 * (i - 1) + 1];
			s11 = new point(meshVerts[i].Position[0], meshVerts[i].Position[1], meshVerts[i].Position[2]);
			u11 = uv[2 * i] - 1;
			v11 = uv[2 * i + 1] - 1;
			s12 = new point(meshVerts[i + 1].Position[0], meshVerts[i + 1].Position[1], meshVerts[i + 1].Position[2]);
			u12 = uv[2 * (i + 1)];
			v12 = uv[2 * (i + 1) + 1];
			s20 = new point(meshVerts[i + 20].Position[0], meshVerts[i + 20].Position[1], meshVerts[i + 20].Position[2]);
			u20 = uv[2 * (i + 20)];
			v20 = uv[2 * (i + 20) + 1];
			s21 = new point(meshVerts[i + 21].Position[0], meshVerts[i + 21].Position[1], meshVerts[i + 21].Position[2]);
			u21 = uv[2 * (i + 21)];
			v21 = uv[2 * (i + 21) + 1];
			s22 = new point(meshVerts[i + 22].Position[0], meshVerts[i + 22].Position[1], meshVerts[i + 22].Position[2]);
			u22 = uv[2 * (i + 22)];
			v22 = uv[2 * (i + 22) + 1];
		}

		point c00 = (*s11 * (float)(16.0f / 36.0f)) + ((*s21 + *s12 + *s01 + *s10) * (float)(4.0f / 36.0f)) + ((*s22 + *s02 + *s00 + *s20) * (float)(1.0f / 36.0f));
		float nu00 = (u11 * (float)(16.0f / 36.0f)) + ((u21 + u12 + u01 + u10) * (float)(4.0f / 36.0f)) + ((u22 + u02 + u00 + u20) * (float)(1.0f / 36.0f));
		float nv00 = (v11 * (float)(16.0f / 36.0f)) + ((v21 + v12 + v01 + v10) * (float)(4.0f / 36.0f)) + ((v22 + v02 + v00 + v20) * (float)(1.0f / 36.0f));

		point c01 = (*s11 * (float)(8.0f / 18.0f)) + ((*s01 + *s21) * (float)(2.0f / 18.0f)) + (*s12 * (float)(4.0f / 18.0f)) + ((*s22 + *s02) * (float)(1.0f / 18.0f));
		float nu01 = (u11 * (float)(8.0f / 18.0f)) + ((u01 + u21) * (float)(2.0f / 18.0f)) + (u12 * (float)(4.0f / 18.0f)) + ((u22 + u02) * (float)(1.0f / 18.0f));
		float nv01 = (v11 * (float)(8.0f / 18.0f)) + ((v01 + v21) * (float)(2.0f / 18.0f)) + (v12 * (float)(4.0f / 18.0f)) + ((v22 + v02) * (float)(1.0f / 18.0f));

		point c02 = (*s12 * (float)(8.0f / 18.0f)) + ((*s02 + *s22) * (float)(2.0f / 18.0f)) + (*s11 * (float)(4.0f / 18.0f)) + ((*s21 + *s01) * (float)(1.0f / 18.0f));
		float nu02 = (u12 * (float)(8.0f / 18.0f)) + ((u02 + u22) * (float)(2.0f / 18.0f)) + (u11 * (float)(4.0f / 18.0f)) + ((u21 + u01) * (float)(1.0f / 18.0f));
		float nv02 = (v12 * (float)(8.0f / 18.0f)) + ((v02 + v22) * (float)(2.0f / 18.0f)) + (v11 * (float)(4.0f / 18.0f)) + ((v21 + v01) * (float)(1.0f / 18.0f));

		point c10 = (*s11 * (float)(8.0f / 18.0f)) + ((*s10 + *s12) * (float)(2.0f / 18.0f)) + (*s21 * (float)(4.0f / 18.0f)) + ((*s22 + *s20) * (float)(1.0f / 18.0f));
		float nu10 = (u11 * (float)(8.0f / 18.0f)) + ((u10 + u12) * (float)(2.0f / 18.0f)) + (u21 * (float)(4.0f / 18.0f)) + ((u22 + u20) * (float)(1.0f / 18.0f));
		float nv10 = (v11 * (float)(8.0f / 18.0f)) + ((v10 + v12) * (float)(2.0f / 18.0f)) + (v21 * (float)(4.0f / 18.0f)) + ((v22 + v20) * (float)(1.0f / 18.0f));

		point c11 = (*s11 * (float)(4.0f / 9.0f)) + ((*s21 + *s12) * (float)(2.0f / 9.0f)) + ((*s22) * (float)(1.0f / 9.0f));
		float nu11 = (u11 * (float)(4.0f / 9.0f)) + ((u21 + u12) * (float)(2.0f / 9.0f)) + ((u22) * (float)(1.0f / 9.0f));
		float nv11 = (v11 * (float)(4.0f / 9.0f)) + ((v21 + v12) * (float)(2.0f / 9.0f)) + ((v22) * (float)(1.0f / 9.0f));

		point c12 = (*s12 * (float)(4.0f / 9.0f)) + ((*s11 + *s22) * (float)(2.0f / 9.0f)) + (*s21 * (float)(1.0f / 9.0f));
		float nu12 = (u12 * (float)(4.0f / 9.0f)) + ((u11 + u22) * (float)(2.0f / 9.0f)) + (u21 * (float)(1.0f / 9.0f));
		float nv12 = (v12 * (float)(4.0f / 9.0f)) + ((v11 + v22) * (float)(2.0f / 9.0f)) + (v21 * (float)(1.0f / 9.0f));

		point c20 = (*s21 * (float)(8.0f / 18.0f)) + ((*s20 + *s22) * (float)(2.0f / 18.0f)) + (*s11 * (float)(4.0f / 18.0f) + (*s12 + *s10) * (float)(1.0f / 18.0f));
		float nu20 = (u21 * (float)(8.0f / 18.0f)) + ((u20 + u22) * (float)(2.0f / 18.0f)) + (u11 * (float)(4.0f / 18.0f) + (u12 + u10) * (float)(1.0f / 18.0f));
		float nv20 = (v21 * (float)(8.0f / 18.0f)) + ((v20 + v22) * (float)(2.0f / 18.0f)) + (v11 * (float)(4.0f / 18.0f) + (v12 + v10) * (float)(1.0f / 18.0f));

		point c21 = (*s21 * (float)(4.0f / 9.0f)) + ((*s11 + *s22) * (float)(2.0f / 9.0f)) + (*s12 * (float)(1.0f / 9.0f));
		float nu21 = (u21 * (float)(4.0f / 9.0f)) + ((u11 + u22) * (float)(2.0f / 9.0f)) + (u12 * (float)(1.0f / 9.0f));
		float nv21 = (v21 * (float)(4.0f / 9.0f)) + ((v11 + v22) * (float)(2.0f / 9.0f)) + (v12 * (float)(1.0f / 9.0f));

		point c22 = (*s22 * (float)(4.0f / 9.0f)) + ((*s11 + *s22) * (float)(2.0f / 9.0f)) + (*s11 * (float)(1.0f / 9.0f));
		float nu22 = (u22 * (float)(4.0f / 9.0f)) + ((u11 + u22) * (float)(2.0f / 9.0f)) + (u11 * (float)(1.0f / 9.0f));
		float nv22 = (v22 * (float)(4.0f / 9.0f)) + ((v11 + v22) * (float)(2.0f / 9.0f)) + (v11 * (float)(1.0f / 9.0f));


		subdVerts[3 * j].Position[0] = c00.x;
		subdVerts[3 * j].Position[1] = c00.y;
		subdVerts[3 * j].Position[2] = s11->z;
		subdVerts[3 * j].Position[3] = 1.0f;
		subdVerts[3 * j].SetColor(colorRed);
		subdVerts[3 * j].SetNormal(meshNormal);
		uvSubdiv[2 * (3 * j)] = nu00;
		uvSubdiv[2 * (3 * j) + 1] = nv00;

		if ((i + 1) % 21 != 0) {
			subdVerts[3 * j + 1].Position[0] = c01.x;
			subdVerts[3 * j + 1].Position[1] = c01.y;
			subdVerts[3 * j + 1].Position[2] = s11->z;
			subdVerts[3 * j + 1].Position[3] = 1.0f;
			subdVerts[3 * j + 1].SetColor(colorRed);
			subdVerts[3 * j + 1].SetNormal(meshNormal);
			uvSubdiv[2 * (3 * j + 1)] = nu01;
			uvSubdiv[2 * (3 * j + 1) + 1] = nv01;

			subdVerts[3 * j + 2].Position[0] = c02.x;
			subdVerts[3 * j + 2].Position[1] = c02.y;
			subdVerts[3 * j + 2].Position[2] = s11->z;
			subdVerts[3 * j + 2].Position[3] = 1.0f;
			subdVerts[3 * j + 2].SetColor(colorRed);
			subdVerts[3 * j + 2].SetNormal(meshNormal);
			uvSubdiv[2 * (3 * j + 2)] = nu02;
			uvSubdiv[2 * (3 * j + 2) + 1] = nv02;
		}

		if (i < 420) {
			subdVerts[3 * j + 61].Position[0] = c10.x;
			subdVerts[3 * j + 61].Position[1] = c10.y;
			subdVerts[3 * j + 61].Position[2] = s11->z;
			subdVerts[3 * j + 61].Position[3] = 1.0f;
			subdVerts[3 * j + 61].SetColor(colorRed);
			subdVerts[3 * j + 61].SetNormal(meshNormal);
			uvSubdiv[2 * (3 * j + 61)] = nu10;
			uvSubdiv[2 * (3 * j + 61) + 1] = nv10;

			if ((i + 1) % 21 != 0) {
				subdVerts[3 * j + 62].Position[0] = c11.x;
				subdVerts[3 * j + 62].Position[1] = c11.y;
				subdVerts[3 * j + 62].Position[2] = s11->z;
				subdVerts[3 * j + 62].Position[3] = 1.0f;
				subdVerts[3 * j + 62].SetColor(colorRed);
				subdVerts[3 * j + 62].SetNormal(meshNormal);
				uvSubdiv[2 * (3 * j + 62)] = nu11;
				uvSubdiv[2 * (3 * j + 62) + 1] = nv11;

				subdVerts[3 * j + 63].Position[0] = c12.x;
				subdVerts[3 * j + 63].Position[1] = c12.y;
				subdVerts[3 * j + 63].Position[2] = s11->z;
				subdVerts[3 * j + 63].Position[3] = 1.0f;
				subdVerts[3 * j + 63].SetColor(colorRed);
				subdVerts[3 * j + 63].SetNormal(meshNormal);
				uvSubdiv[2 * (3 * j + 63)] = nu12;
				uvSubdiv[2 * (3 * j + 63) + 1] = nv12;
			}

			subdVerts[3 * j + 122].Position[0] = c20.x;
			subdVerts[3 * j + 122].Position[1] = c20.y;
			subdVerts[3 * j + 122].Position[2] = s11->z;
			subdVerts[3 * j + 122].Position[3] = 1.0f;
			subdVerts[3 * j + 122].SetColor(colorRed);
			subdVerts[3 * j + 122].SetNormal(meshNormal);
			uvSubdiv[2 * (3 * j + 122)] = nu20;
			uvSubdiv[2 * (3 * j + 122) + 1] = nv20;

			if ((i + 1) % 21 != 0) {
				subdVerts[3 * j + 123].Position[0] = c21.x;
				subdVerts[3 * j + 123].Position[1] = c21.y;
				subdVerts[3 * j + 123].Position[2] = s11->z;
				subdVerts[3 * j + 123].Position[3] = 1.0f;
				subdVerts[3 * j + 123].SetColor(colorRed);
				subdVerts[3 * j + 123].SetNormal(meshNormal);
				uvSubdiv[2 * (3 * j + 123)] = nu21;
				uvSubdiv[2 * (3 * j + 123) + 1] = nv21;

				subdVerts[3 * j + 124].Position[0] = c22.x;
				subdVerts[3 * j + 124].Position[1] = c22.y;
				subdVerts[3 * j + 124].Position[2] = s11->z;
				subdVerts[3 * j + 124].Position[3] = 1.0f;
				subdVerts[3 * j + 124].SetColor(colorRed);
				subdVerts[3 * j + 124].SetNormal(meshNormal);
				uvSubdiv[2 * (3 * j + 124)] = nu22;
				uvSubdiv[2 * (3 * j + 124) + 1] = nv22;
			}
		}

		if (i != 0 && (i + 1) % 21 == 0) {
			//cout << "i " << i << " j " << j << endl;
			j = j + 41;
		}
		else {
			j++;
		}
	}

	for (int i = 0; i < 3721; i++) {
		if ((i + 1) % 61 != 0 && i < 3660 && i != 3720) {
			subdIdcs[4 * i] = i;
			subdIdcs[4 * i + 1] = i + 1;
		}
		else {
			subdIdcs[4 * i] = i;
			subdIdcs[4 * i + 1] = i;
		}
		if (i < 3660) {
			subdIdcs[4 * i + 2] = i;
			subdIdcs[4 * i + 3] = i + 61;
		}
		else if (i != 3720) {
			subdIdcs[4 * i + 2] = i;
			subdIdcs[4 * i + 3] = i + 1;
		}

		if (i == 0 || (i + 1) % 61 != 0 && i < 3660) {
			subdTexIdcs[6 * i] = i;
			subdTexIdcs[6 * i + 1] = i + 1;
			subdTexIdcs[6 * i + 2] = i + 62;
			subdTexIdcs[6 * i + 3] = i + 62;
			subdTexIdcs[6 * i + 4] = i + 61;
			subdTexIdcs[6 * i + 5] = i;
		}

		/*uvSubdiv[2 * i] = (subdVerts[i].Position[0] + 10) / 60;
		uvSubdiv[2 * i + 1] = (subdVerts[i].Position[1]) / 60;*/
	}

	VertexBufferSize[5] = sizeof(subdVerts);
	IndexBufferSize[5] = sizeof(subdIdcs);
	createVAOs(subdVerts, subdIdcs, 5);

	VertexBufferSize[6] = sizeof(subdVerts);
	IndexBufferSize[6] = sizeof(subdTexIdcs);
	createVAOsForTex(subdVerts, subdTexIdcs, 6);
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

		if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT)) {		
			move_vertex();
		}
		// DRAWING POINTS
		renderScene();


	} // Check if the ESC key was pressed or the window was closed
	while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
		glfwWindowShouldClose(window) == 0);

	cleanup();

	return 0;
}
