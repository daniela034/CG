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

struct Geometric_Transf{
    float t_x, t_y, t_z;             //Translate
    float angle, r_x, r_y, r_z;      //Rotate
    float s_x, s_y, s_z;             //Scale
};

struct Group{
    Geometric_Transf transform;
    vector<string> models;
    vector<Group> subgroups;
    vector<Point> vertices;
};

vector<Group> global_groups;

float alpha = -2.7, beta = -0.7;     // ângulos que definem a direção do olhar

float px = 10, py = 20, pz = 20;     // posição da câmara
float dx, dy, dz;                    // vetor direção do olhar
float right_x, right_y, right_z;     // vetor "right" calculado pelo produto externo entre vetor direção e vetor up
float dist = 1;                      // distância que a câmara se moviemnta com cada clique das teclas

int polygon_mode = 0;

void spherical2Cartesian() {
    dx = cos(beta) * sin(alpha);
    dy = sin(beta);
    dz = cos(beta) * cos(alpha);
}

/* Produto externo entre o vetor direção (dx, dy, dz) e o vetor up (0, 1, 0) para podermos movimentar a câmara para a esquerda e direita */
void right_vector(){
    right_x = -dz;
    right_y = 0;
    right_z = dx;
}

void init_Geometric_Transf(Geometric_Transf &gt){
    //TRANSLATE
    gt.t_x = 0;
    gt.t_y = 0;
    gt.t_z = 0;

    //ROTATE
    gt.angle = 0;
    gt.r_x = 0;
    gt.r_y = 0;
    gt.r_z = 0;

    //SCALE
    gt.s_x = 1;
    gt.s_y = 1;
    gt.s_z = 1;
}

void load_Transforms(XMLElement* translate, XMLElement* rotate, XMLElement* scale, Geometric_Transf &gt){
    const char* x = nullptr; const char* y = nullptr; const char* z = nullptr;

    if (translate != nullptr){

        x = translate->Attribute("X");
        y = translate->Attribute("Y");
        z = translate->Attribute("Z");

        if (x != nullptr) gt.t_x = stof(x);
        if (y != nullptr) gt.t_y = stof(y);
        if (z != nullptr) gt.t_z = stof(z);
    }

    if (rotate != nullptr){
        const char* angle = nullptr;

        angle = rotate->Attribute("angle");
        x = rotate->Attribute("axisX");
        y = rotate->Attribute("axisY");
        z = rotate->Attribute("axisZ");

        if (angle != nullptr) gt.angle = stof(angle);
        if (x != nullptr) gt.r_x = stof(x);
        if (y != nullptr) gt.r_y = stof(y);
        if (z != nullptr) gt.r_z = stof(z);
    }

    if(scale != nullptr){
        x = scale->Attribute("X");
        y = scale->Attribute("Y");
        z = scale->Attribute("Z");

        if (x != nullptr) gt.s_x = stof(x);
        if (y != nullptr) gt.s_y = stof(y);
        if (z != nullptr) gt.s_z = stof(z);
    }
}

void load_Models(XMLElement* modelsXML, vector<string> &models){
    XMLElement * pListModels = modelsXML->FirstChildElement("model");

    const char* model = nullptr;

    while (pListModels != nullptr) {
        model = pListModels->Attribute("file");
        models.emplace_back(model);
        pListModels = pListModels->NextSiblingElement("model");
    }
}

void load_Groups(XMLElement* pListGroups, Group &g){

    //Preencher a estrutura Geometric_Transf do grupo g
    load_Transforms(pListGroups->FirstChildElement("translate"), pListGroups->FirstChildElement("rotate"),
            pListGroups->FirstChildElement("scale"), g.transform);

    //Preencher o vetor models do grupo g
    load_Models(pListGroups->FirstChildElement("models"), g.models);

    //Tratar dos sub-grupos do grupo g chamando a função recursivamente e preenchendo o vetor subgroups do grupo g
    XMLElement * pListSubGroups = pListGroups->FirstChildElement("group");

    while (pListSubGroups != nullptr) {
        Group sub_g;
        init_Geometric_Transf(sub_g.transform);

        load_Groups(pListSubGroups, sub_g);
        g.subgroups.push_back(sub_g);
        pListSubGroups = pListSubGroups->NextSiblingElement("group");
    }
}

XMLError load_XML_file(char* file){
    XMLDocument xmlDoc;
    XMLError eResult = xmlDoc.LoadFile(file);
    XMLCheckResult(eResult);

    XMLNode* pRoot = xmlDoc.FirstChild();
    if (pRoot == nullptr)
        return XML_ERROR_FILE_READ_ERROR;

    XMLElement* pListGroups = pRoot->FirstChildElement("group");

    while (pListGroups != nullptr) {
        Group g;
        init_Geometric_Transf(g.transform);

        load_Groups(pListGroups, g);
        global_groups.push_back(g);
        pListGroups = pListGroups->NextSiblingElement("group");
    }

    return XML_SUCCESS;
}

void load_model_file_aux(string file, vector<Point> &vertices){
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

void load_model_files(vector<Group> &vector_groups){

    for (Group &g : vector_groups) {
        for (auto &model : g.models) {
            load_model_file_aux(model, g.vertices);
        }
        load_model_files(g.subgroups);
    }
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
    glVertex3f( 150.0f, 0.0f, 0.0f);

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

void draw_Groups(vector<Group> vector_groups){

    for (Group g : vector_groups) {
        glPushMatrix();
        glTranslatef(g.transform.t_x, g.transform.t_y, g.transform.t_z);
        glRotatef(g.transform.angle, g.transform.r_x, g.transform.r_y, g.transform.r_z);
        glScalef(g.transform.s_x, g.transform.s_y, g.transform.s_z);

        glBegin(GL_TRIANGLES);
        for(Point p : g.vertices)
            glVertex3f(p.x, p.y, p.z);
        glEnd();

        draw_Groups(g.subgroups);
        glPopMatrix();
    }
}


void renderScene(void) {

    // clear buffers
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // set the camera
    glLoadIdentity();
    gluLookAt(px, py, pz,
              px + dx,py + dy,pz + dz,
              0.0f,1.0f,0.0f);

// put drawing instructions here
    draw_referencial();

    glColor3f(0.5, 0.5, 0.5);
    draw_Groups(global_groups);

    // End of frame
    glutSwapBuffers();
}


// write function to process keyboard events

void special_keys(int key_code, int x, int y){
    switch (key_code){
        case GLUT_KEY_LEFT :
            alpha += 0.05;
            break;

        case GLUT_KEY_RIGHT :
            alpha -= 0.05;
            break;

        case GLUT_KEY_UP :
            if (beta <= 1.5f)
                beta += 0.05f;
            break;

        case GLUT_KEY_DOWN :
            if (beta >= -1.5f)
                beta -= 0.05f;
            break;
    }

    spherical2Cartesian();
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

        case 'w' :
            px += (dx * dist);
            py += (dy * dist);
            pz += (dz * dist);
            break;

        case 's' :
            px -= (dx * dist);
            py -= (dy * dist);
            pz -= (dz * dist);
            break;

        case 'a' :
            right_vector();
            px -= (right_x * dist);
            py -= (right_y * dist);
            pz -= (right_z * dist);
            break;

        case 'd' :
            right_vector();
            px += (right_x * dist);
            py += (right_y * dist);
            pz += (right_z * dist);
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
        load_model_files(global_groups);
    }


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

    spherical2Cartesian();

// enter GLUT's main cycle
    glutMainLoop();

    return 1;
}

