#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <math.h>
#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include "tinyxml2/tinyxml2.h"

#ifndef XMLCheckResult
#define XMLCheckResult(a_eResult) if (a_eResult != XML_SUCCESS) { printf("Error: %i\n", a_eResult); return a_eResult; }
#endif

using namespace std;
using namespace tinyxml2;

struct Point{
    float x;
    float y;
    float z;
};

vector<Point> vertices;
vector<string> models;

float alpha, beta = 0;
float radius = 15;

int polygon_mode = 0;

XMLError load_XML_file(char* file){
    XMLDocument xmlDoc;
    XMLError eResult = xmlDoc.LoadFile(file);
    XMLCheckResult(eResult);

    XMLNode * pRoot = xmlDoc.FirstChild();
    if (pRoot == nullptr) return XML_ERROR_FILE_READ_ERROR;

    XMLElement * pListModels = pRoot->FirstChildElement("model");

    const char* model = nullptr;

    while (pListModels != nullptr) {
        model = pListModels->Attribute("file");
        models.push_back(model);
        pListModels = pListModels->NextSiblingElement("model");
    }

    return XML_SUCCESS;
}

void load_model_file(string file){
    ifstream file_pointer(file, ios::binary | ios::in);

    Point p;
    string linha;
    float coord;

    while (getline(file_pointer, linha, '\0')) {
        stringstream ss(linha);
        ss >> coord;
        p.x = coord;
        ss >> coord;
        p.y = coord;
        ss >> coord;
        p.z = coord;

        vertices.push_back(p);
    }
    file_pointer.close();
}

void changeSize(int w, int h) {

    // Prevent a divide by zero, when window is too short
    // (you cant make a window with zero width).
    if(h == 0)
        h = 1;

    // compute window's aspect ratio
    float ratio = w * 1.0 / h;

    // Set the projection matrix as current
    glMatrixMode(GL_PROJECTION);
    // Load Identity Matrix
    glLoadIdentity();

    // Set the viewport to be the entire window
    glViewport(0, 0, w, h);

    // Set perspective
    gluPerspective(45.0f ,ratio, 1.0f ,1000.0f);

    // return to the model view matrix mode
    glMatrixMode(GL_MODELVIEW);
}

void draw_referencial(){
    glBegin(GL_LINES);
    // X axis in red
    glColor3f(1.0f, 0.0f, 0.0f);
    glVertex3f(-100.0f, 0.0f, 0.0f);
    glVertex3f( 100.0f, 0.0f, 0.0f);

    // Y Axis in Green
    glColor3f(0.0f, 1.0f, 0.0f);
    glVertex3f(0.0f, -100.0f, 0.0f);
    glVertex3f(0.0f,  100.0f, 0.0f);

    // Z Axis in Blue
    glColor3f(0.0f, 0.0f, 1.0f);
    glVertex3f(0.0f, 0.0f, -100.0f);
    glVertex3f(0.0f, 0.0f,  100.0f);
    glEnd();
}

void renderScene(void) {

    // clear buffers
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // set the camera
    glLoadIdentity();
    gluLookAt(radius * cos(beta) * sin(alpha), radius * sin(beta),radius * cos(beta) * cos(alpha),
              0.0,0.0,0.0,
              0.0f,1.0f,0.0f);


// put drawing instructions here
    draw_referencial();

    glColor3f(0.5,0.5,0.5);
    glBegin(GL_TRIANGLES);
    for (Point p : vertices)
        glVertex3f(p.x, p.y, p.z);
    glEnd();

    // End of frame
    glutSwapBuffers();
}


// write function to process keyboard events

void special_keys(int key_code, int x, int y){
    switch (key_code){
        case GLUT_KEY_LEFT :
            alpha -= 0.1;
            break;

        case GLUT_KEY_RIGHT :
            alpha += 0.1;
            break;

        case GLUT_KEY_UP :
            if (beta <= (M_PI / 2))
                beta += 0.1;
            break;

        case GLUT_KEY_DOWN :
            if (beta >= -(M_PI / 2))
                beta -= 0.1;
            break;

        case GLUT_KEY_PAGE_UP :
            radius -= 0.1;
            break;

        case GLUT_KEY_PAGE_DOWN :
            radius += 0.1;
            break;
    }
    glutPostRedisplay();
}

void regular_keys (unsigned char key, int x, int y){
    switch(key){
        case 'p' :
            if (polygon_mode == 0){
                polygon_mode = 1;
                glPolygonMode(GL_FRONT, GL_LINE);
            }
            else {
                polygon_mode = 0;
                glPolygonMode(GL_FRONT, GL_FILL);
            }
            break;
    }
    glutPostRedisplay();
}

int main(int argc, char **argv) {
    if (argc < 2) {
        cout << "Insira o ficheiro de configuração xml como argumento" << endl;
        return 1;
    }
    else {
        load_XML_file(argv[1]);
    }

    for (int i = 0; i < models.size(); i++)
        load_model_file(models[i]);

// init GLUT and the window
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH|GLUT_DOUBLE|GLUT_RGBA);
    glutInitWindowPosition(100,100);
    glutInitWindowSize(800,800);
    glutCreateWindow("CG@DI-UM");

// Required callback registry
    glutDisplayFunc(renderScene);
    glutReshapeFunc(changeSize);


// put here the registration of the keyboard callbacks
    glutSpecialFunc(special_keys);
    glutKeyboardFunc(regular_keys);


//  OpenGL settings
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);

// enter GLUT's main cycle
    glutMainLoop();

    return 1;
}

