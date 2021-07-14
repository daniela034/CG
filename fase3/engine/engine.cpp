#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glew.h>
#include <GL/glut.h>
#endif

#include <math.h>
#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include "tinyxml2/tinyxml2.h"

#include "catmullRom.cpp"

#ifndef XMLCheckResult
#define XMLCheckResult(a_eResult) if (a_eResult != XML_SUCCESS) { printf("Error: %i\n", a_eResult); return a_eResult; }
#endif

using namespace std;
using namespace tinyxml2;

struct Geometric_Transf{
    float t_x, t_y, t_z;             //Translate
    float angle, r_x, r_y, r_z;      //Rotate
    float s_x, s_y, s_z;             //Scale

    /* Translate with Animation */
    float ta_time;
    vector<Point> ta_points;
    float ta_Y[3];
    int ta_turns;

    /* Rotate with Animation */
    float ra_time;
    float ra_x, ra_y, ra_z;
};

struct Group{
    Geometric_Transf transform;
    vector<string> models;
    vector<Group> subgroups;
    vector<float> vertexB;
    GLuint coords[1];
    GLuint num_vertices;
};

vector<Group> global_groups;

float alpha = -2.7, beta = -0.7;     // ângulos que definem a direção do olhar

float px = 10, py = 20, pz = 20;     // posição da câmara
float dx, dy, dz;                    // vetor direção do olhar
float right_x, right_y, right_z;     // vetor "right" calculado pelo produto externo entre vetor direção e vetor up
float dist = 1;                      // distância que a câmara se moviemnta com cada clique das teclas

int polygon_mode = 0;
int ref_lines = 0;

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

    //TRANSLATE ANIMATION
    gt.ta_time = 0;
    gt.ta_Y[0] = 0; gt.ta_Y[1] = 1; gt.ta_Y[2] = 0;
    gt.ta_turns = 0;

    //ROTATE ANIMATION
    gt.ra_time = 0;
    gt.ra_x = 0; gt.ra_y = 0; gt.ra_z = 0;
}

void load_Translate_Animation(XMLElement* translate, Geometric_Transf &gt){
    const char* time = nullptr;
    const char* x = nullptr; const char* y = nullptr; const char* z = nullptr;
    vector<Point> pts;
    Point p;

    time = translate->Attribute("time");
    if (time != nullptr)
        gt.ta_time = stof(time);

    XMLElement* ptList = translate->FirstChildElement("point");

    while (ptList != nullptr){
        p.x = 0; p.y = 0; p.z = 0;

        x = ptList->Attribute("X");
        y = ptList->Attribute("Y");
        z = ptList->Attribute("Z");

        if (x != nullptr) p.x = stof(x);
        if (y != nullptr) p.y = stof(y);
        if (z != nullptr) p.z = stof(z);

        pts.push_back(p);
        ptList = ptList->NextSiblingElement("point");
    }

    if (pts.size() > 3){
        gt.ta_points = pts;
    }
}

void load_Rotate_Animation(XMLElement* rotate, Geometric_Transf &gt){
    const char* time = nullptr;
    const char* x = nullptr; const char* y = nullptr; const char* z = nullptr;

    time = rotate->Attribute("time");
    if (time != nullptr)
        gt.ra_time = stof(time);

    x = rotate->Attribute("axisX");
    y = rotate->Attribute("axisY");
    z = rotate->Attribute("axisZ");

    if (x != nullptr) gt.ra_x = stof(x);
    if (y != nullptr) gt.ra_y = stof(y);
    if (z != nullptr) gt.ra_z = stof(z);

}

void load_Transforms(XMLElement* translate, XMLElement* rotate, XMLElement* scale, Geometric_Transf &gt){
    const char* x = nullptr; const char* y = nullptr; const char* z = nullptr;

    if (translate != nullptr){

        x = translate->Attribute("X");

        if (x != nullptr) {
            gt.t_x = stof(x);

            y = translate->Attribute("Y");
            if (y != nullptr) gt.t_y = stof(y);

            z = translate->Attribute("Z");
            if (z != nullptr) gt.t_z = stof(z);
        }

        else {
            load_Translate_Animation(translate, gt);
        }
    }

    if (rotate != nullptr){
        const char* angle = nullptr;

        angle = rotate->Attribute("angle");

        if (angle != nullptr) {
            gt.angle = stof(angle);
            x = rotate->Attribute("axisX");
            y = rotate->Attribute("axisY");
            z = rotate->Attribute("axisZ");

            if (x != nullptr) gt.r_x = stof(x);
            if (y != nullptr) gt.r_y = stof(y);
            if (z != nullptr) gt.r_z = stof(z);
        }

        else {
            load_Rotate_Animation(rotate, gt);
        }
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

void load_model_file_aux(string file, vector<float> &vertexB, GLuint &num_vertices){
    ifstream file_pointer(file, ios::binary | ios::in);

    string linha;
    float coord;

    while (getline(file_pointer, linha, '\0')) {
        stringstream ss(linha);
        ss >> coord;
        vertexB.push_back(coord);
        ss >> coord;
        vertexB.push_back(coord);
        ss >> coord;
        vertexB.push_back(coord);

        num_vertices += 1;
    }
    file_pointer.close();
}

void load_model_files(vector<Group> &vector_groups){

    for (Group &g : vector_groups) {
        for (auto &model : g.models) {
            load_model_file_aux(model, g.vertexB, g.num_vertices);
        }

        glGenBuffers(1, g.coords);
        glBindBuffer(GL_ARRAY_BUFFER, g.coords[0]);
        glBufferData(GL_ARRAY_BUFFER, g.vertexB.size() * sizeof(float), g.vertexB.data(), GL_STATIC_DRAW);

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

void apply_translate_anim(Group &g){
    float pos[3], deriv[3], m[16], Z[3];

    renderCatmullRomCurve(g.transform.ta_points);

    float elapsedT = glutGet(GLUT_ELAPSED_TIME) / 1000.0;
    float actualT = elapsedT - g.transform.ta_time * g.transform.ta_turns;

    if(actualT > g.transform.ta_time){
        g.transform.ta_turns++;
        actualT = elapsedT - g.transform.ta_time * g.transform.ta_turns;
    }

    float global_t = actualT / g.transform.ta_time;

    getGlobalCatmullRomPoint(global_t, g.transform.ta_points, pos, deriv);

    // Z = X * Yi - 1
    normalize(deriv);
    cross(deriv, g.transform.ta_Y, Z);
    normalize(Z);

    //Yi = Z * X
    cross(Z, deriv, g.transform.ta_Y);
    normalize(g.transform.ta_Y);

    buildRotMatrix(deriv, g.transform.ta_Y, Z, m);
    glTranslatef(pos[0], pos[1], pos[2]);
    glMultMatrixf(m);
}

void apply_rotate_anim(Group g){

    float ang = 360 / (g.transform.ra_time * 1000);

    float elapsedT = glutGet(GLUT_ELAPSED_TIME);

    glRotatef(elapsedT * ang, g.transform.ra_x, g.transform.ra_y, g.transform.ra_z);
}

void apply_transforms(Group &g){

    /* Transformações com Animação */
    if (g.transform.ta_time != 0)
        apply_translate_anim(g);

    if (g.transform.ra_time != 0)
        apply_rotate_anim(g);

    /* Transformações base (sem animação) */
    glTranslatef(g.transform.t_x, g.transform.t_y, g.transform.t_z);
    glRotatef(g.transform.angle, g.transform.r_x, g.transform.r_y, g.transform.r_z);
    glScalef(g.transform.s_x, g.transform.s_y, g.transform.s_z);

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

        apply_transforms(g);

        glBindBuffer(GL_ARRAY_BUFFER, g.coords[0]);
        glVertexPointer(3, GL_FLOAT, 0, 0);
        glDrawArrays(GL_TRIANGLES, 0, g.num_vertices);

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
    if (ref_lines == 0)
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

        case 'l' :
            if (ref_lines == 0){
                ref_lines = 1;
            }
            else {
                ref_lines = 0;
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

// init GLUT and the window
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH|GLUT_DOUBLE|GLUT_RGBA);
    glutInitWindowPosition(100,100);
    glutInitWindowSize(1200,800);
    glutCreateWindow("CG@DI-UM");

// Required callback registry
    glutDisplayFunc(renderScene);
    glutIdleFunc(renderScene);
    glutReshapeFunc(changeSize);


// put here the registration of the keyboard callbacks
    glutSpecialFunc(special_keys);
    glutKeyboardFunc(regular_keys);

    glewInit();

//  OpenGL settings
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glEnableClientState(GL_VERTEX_ARRAY);

    if (argc < 2) {
        cout << "Insira o ficheiro de configuração xml como argumento" << endl;
        return 1;
    }
    else {
        load_XML_file(argv[1]);
        load_model_files(global_groups);
    }

    spherical2Cartesian();

// enter GLUT's main cycle
    glutMainLoop();

    return 1;
}

