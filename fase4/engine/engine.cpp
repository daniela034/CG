#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glew.h>
#include <GL/glut.h>
#endif

#include <IL/il.h>

#include <math.h>
#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <map>
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

struct Light{
    int type;
    float *pos, *dir;
    float cutoff, linear, constant, quadratic, exponent;
};

struct Prim{
    vector<float> vertexB;
    vector<float> normalB;
    vector<float> text_coordsB;

    GLuint coords[1];
    GLuint normals[1];
    GLuint texCoord[1];

    GLuint num_vertices;
};

struct Model{
    string file;
    float *amb, *diff, *spec, *emiss;
    float shininess;

    int id_textura;
};

struct Group{
    Geometric_Transf transform;
    vector<Model> models;
    vector<Group> subgroups;
};

map<string, Prim> primitives;
vector<Group> global_groups;
vector<Light> lightsV;
map<string, int> textures_map;

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

void init_Model(Model &m){
    m.amb = new float[4];
    m.diff = new float[4];
    m.spec = new float[4];
    m.emiss = new float[4];

    m.amb[0] = 0.2f; m.amb[1] = 0.2f; m.amb[2] = 0.2f; m.amb[3] = 1.0f;
    m.diff[0] = 0.8f; m.diff[1] = 0.8f; m.diff[2] = 0.8f; m.diff[3] = 1.0f;
    m.spec[0] = 0.0f; m.spec[1] = 0.0f; m.spec[2] = 0.0f; m.spec[3] = 1.0f;
    m.emiss[0] = 0.0f; m.emiss[1] = 0.0f; m.emiss[2] = 0.0f; m.emiss[3] = 1.0f;

    m.shininess = 0;
}

void init_Light(Light &l){
    l.pos = new float[4];
    l.dir = new float[3];

    l.pos[0] = 0.0f; l.pos[1] = 0.0f; l.pos[2] = 1.0f; l.pos[3] = 1.0f;
    l.dir[0] = 0.0f; l.dir[1] = 0.0f; l.dir[2] = -1.0f;

    l.cutoff = 180;
    l.linear = 0;
    l.constant = 1;
    l.exponent = 0;
    l.quadratic = 0;
}

int loadTexture(string s) {

    unsigned int t,tw,th;
    unsigned char *texData;
    unsigned int texID;

    ilInit();
    ilEnable(IL_ORIGIN_SET);
    ilOriginFunc(IL_ORIGIN_LOWER_LEFT);
    ilGenImages(1,&t);
    ilBindImage(t);
    ilLoadImage((ILstring)s.c_str());
    tw = ilGetInteger(IL_IMAGE_WIDTH);
    th = ilGetInteger(IL_IMAGE_HEIGHT);
    ilConvertImage(IL_RGBA, IL_UNSIGNED_BYTE);
    texData = ilGetData();

    glGenTextures(1,&texID);

    glBindTexture(GL_TEXTURE_2D,texID);
    glTexParameteri(GL_TEXTURE_2D,	GL_TEXTURE_WRAP_S,		GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D,	GL_TEXTURE_WRAP_T,		GL_REPEAT);

    glTexParameteri(GL_TEXTURE_2D,	GL_TEXTURE_MAG_FILTER,   	GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D,	GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, tw, th, 0, GL_RGBA, GL_UNSIGNED_BYTE, texData);
    glGenerateMipmap(GL_TEXTURE_2D);

    glBindTexture(GL_TEXTURE_2D, 0);

    return texID;
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

void load_Models(XMLElement* modelsXML, vector<Model> &models){

    XMLElement * pListModels = modelsXML->FirstChildElement("model");

    const char* model = nullptr;
    const char* texture = nullptr;
    const char * atr = nullptr;

    while (pListModels != nullptr) {
        Model m;
        init_Model(m);

        model = pListModels->Attribute("file");
        m.file = model;

        texture = pListModels->Attribute("texture");

        if (texture != nullptr){
            if (textures_map.find(texture) == textures_map.end()){
                m.id_textura = loadTexture(texture);
                textures_map.insert( pair<string, int> (texture, m.id_textura));
            }
            else {
                m.id_textura = textures_map.find(texture)->second;
            }
        }
        else {
            m.id_textura = 0;
        }

        atr = pListModels->Attribute("ambR");
        if (atr != nullptr){
            m.amb[0] = stof(atr);
        }

        atr = pListModels->Attribute("ambG");
        if (atr != nullptr){
            m.amb[1] = stof(atr);
        }

        atr = pListModels->Attribute("ambB");
        if (atr != nullptr){
            m.amb[2] = stof(atr);
        }

        atr = pListModels->Attribute("diffR");
        if (atr != nullptr){
            m.diff[0] = stof(atr);
        }

        atr = pListModels->Attribute("diffG");
        if (atr != nullptr){
            m.diff[1] = stof(atr);
        }

        atr = pListModels->Attribute("diffB");
        if (atr != nullptr){
            m.diff[2] = stof(atr);
        }

        atr = pListModels->Attribute("specR");
        if (atr != nullptr){
            m.spec[0] = stof(atr);
        }

        atr = pListModels->Attribute("specG");
        if (atr != nullptr){
            m.spec[1] = stof(atr);
        }

        atr = pListModels->Attribute("specB");
        if (atr != nullptr){
            m.spec[2] = stof(atr);
        }

        atr = pListModels->Attribute("emissR");
        if (atr != nullptr){
            m.emiss[0] = stof(atr);
        }

        atr = pListModels->Attribute("emissG");
        if (atr != nullptr){
            m.emiss[1] = stof(atr);
        }

        atr = pListModels->Attribute("emissB");
        if (atr != nullptr){
            m.emiss[2] = stof(atr);
        }

        atr = pListModels->Attribute("shininess");
        if (atr != nullptr){
            m.shininess = stof(atr);
        }

        models.emplace_back(m);
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

void load_lights(XMLElement* lights){

    XMLElement* pListLights = lights->FirstChildElement("light");

    const char * type = nullptr;
    const char* pos_x = nullptr; const char* pos_y = nullptr; const char* pos_z = nullptr;
    const char* dir_x = nullptr; const char* dir_y = nullptr; const char* dir_z = nullptr;
    const char* cutoff = nullptr; const char* constant = nullptr; const char* linear = nullptr; const char* quadratic = nullptr; const char* exponent = nullptr;

    while (pListLights != nullptr){
        Light l;
        init_Light(l);

        type = pListLights->Attribute("type");

        if (type != nullptr){
            if (strcmp(type, "POINT") == 0) {
                l.type = 0;
            }

            else if (strcmp(type, "DIRECTIONAL") == 0) {
                l.type = 1;
            }

            else if (strcmp(type, "SPOT") == 0) {
                l.type = 2;
            }
        }

        pos_x = pListLights->Attribute("posX");
        pos_y = pListLights->Attribute("posY");
        pos_z = pListLights->Attribute("posZ");
        dir_x = pListLights->Attribute("dirX");
        dir_y = pListLights->Attribute("dirY");
        dir_z = pListLights->Attribute("dirZ");
        cutoff = pListLights->Attribute("cutoff");
        constant = pListLights->Attribute("constant");
        linear = pListLights->Attribute("linear");
        quadratic = pListLights->Attribute("quadratic");
        exponent = pListLights->Attribute("exponent");

        if (pos_x != nullptr) l.pos[0] = stof(pos_x);
        if (pos_y != nullptr) l.pos[1] = stof(pos_y);
        if (pos_z != nullptr) l.pos[2] = stof(pos_z);
        if (dir_x != nullptr) l.dir[0] = stof(dir_x);
        if (dir_y != nullptr) l.dir[1] = stof(dir_y);
        if (dir_z != nullptr) l.dir[2] = stof(dir_z);
        if (cutoff != nullptr){
            if (stof(cutoff) >= 0 && stof(cutoff) <= 90){
                l.cutoff = stof(cutoff);
            }
        }
        if (constant != nullptr) l.constant = stof(constant);
        if (linear != nullptr) l.linear = stof(linear);
        if (quadratic != nullptr) l.quadratic = stof(quadratic);
        if (exponent != nullptr) l.exponent = stof(exponent);

        lightsV.push_back(l);
        pListLights = pListLights->NextSiblingElement("light");
    }
}

XMLError load_XML_file(char* file){
    XMLDocument xmlDoc;
    XMLError eResult = xmlDoc.LoadFile(file);
    XMLCheckResult(eResult);

    XMLNode* pRoot = xmlDoc.FirstChild();
    if (pRoot == nullptr)
        return XML_ERROR_FILE_READ_ERROR;

    XMLElement* lights = pRoot->FirstChildElement("lights");

    if (lights != nullptr){
        load_lights(lights);
    }

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

void load_model_file_aux(string file, vector<float> &vertexB, vector<float> &normalB, vector<float> &tex_coordsB, GLuint &num_vertices){
    ifstream file_pointer(file, ios::binary | ios::in);

    string linha;
    float coord;
    float normal;
    float tex_coord;

    while (getline(file_pointer, linha, '\0')) {
        stringstream ss(linha);
        ss >> coord;
        vertexB.push_back(coord);
        ss >> coord;
        vertexB.push_back(coord);
        ss >> coord;
        vertexB.push_back(coord);

        ss >> normal;
        normalB.push_back(normal);
        ss >> normal;
        normalB.push_back(normal);
        ss >> normal;
        normalB.push_back(normal);

        ss >> tex_coord;
        tex_coordsB.push_back(tex_coord);
        ss >> tex_coord;
        tex_coordsB.push_back(tex_coord);

        num_vertices += 1;
    }
    file_pointer.close();
}

void load_model_files(vector<Group> &vector_groups){

    for (Group &g : vector_groups) {
        for (auto &model : g.models) {
            if (primitives.find(model.file) == primitives.end()) {
                Prim p;
                load_model_file_aux(model.file, p.vertexB, p.normalB, p.text_coordsB, p.num_vertices);

                glGenBuffers(1, p.coords);
                glBindBuffer(GL_ARRAY_BUFFER, p.coords[0]);
                glBufferData(GL_ARRAY_BUFFER, p.vertexB.size() * sizeof(float), p.vertexB.data(),
                             GL_STATIC_DRAW);

                glGenBuffers(1, p.normals);
                glBindBuffer(GL_ARRAY_BUFFER, p.normals[0]);
                glBufferData(GL_ARRAY_BUFFER, p.normalB.size() * sizeof(float), p.normalB.data(),
                             GL_STATIC_DRAW);

                glGenBuffers(1, p.texCoord);
                glBindBuffer(GL_ARRAY_BUFFER, p.texCoord[0]);
                glBufferData(GL_ARRAY_BUFFER, p.text_coordsB.size() * sizeof(float), p.text_coordsB.data(),
                             GL_STATIC_DRAW);

                primitives.insert( pair <string, Prim> (model.file, p));
            }
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

void apply_light(Light &l, int i){

    // i = nº de luzes que já foram definidas na cena

    float diff[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
    float spec[4] = { 1.0f, 1.0f, 1.0f, 1.0f };

    glLightfv(i + GL_LIGHT0, GL_DIFFUSE, diff);
    glLightfv(i + GL_LIGHT0, GL_SPECULAR, spec);

    /* Luz do tipo POINT */
    if (l.type == 0) {
        glLightfv(i + GL_LIGHT0, GL_POSITION, l.pos);
    }

    /* Luz do tipo DIRECTIONAL */
    else if (l.type == 1) {

        // Quarto parâmetro a 0 já que neste caso é um vetor
        l.pos[3] = 0;

        glLightfv(i + GL_LIGHT0, GL_POSITION, l.pos);

        glLightf(i + GL_LIGHT0, GL_LINEAR_ATTENUATION, l.linear);
        glLightf(i + GL_LIGHT0, GL_CONSTANT_ATTENUATION, l.constant);
        glLightf(i + GL_LIGHT0, GL_QUADRATIC_ATTENUATION, l.quadratic);

    }

    /* Luz do tipo SPOT */
    else {

        glLightf(i + GL_LIGHT0, GL_LINEAR_ATTENUATION, l.linear);
        glLightf(i + GL_LIGHT0, GL_CONSTANT_ATTENUATION, l.constant);
        glLightf(i + GL_LIGHT0, GL_QUADRATIC_ATTENUATION, l.quadratic);

        glLightfv(i + GL_LIGHT0, GL_POSITION, l.pos);
        glLightfv(i + GL_LIGHT0, GL_SPOT_DIRECTION, l.dir);
        glLightf(i + GL_LIGHT0, GL_SPOT_CUTOFF, l.cutoff);
        glLightf(i + GL_LIGHT0, GL_SPOT_EXPONENT, l.exponent);
    }
}

void draw_referencial(){
    glDisable(GL_LIGHTING);
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
    glEnable(GL_LIGHTING);
}

void draw_Groups(vector<Group> vector_groups){

    for (Group g : vector_groups) {
        glPushMatrix();

        apply_transforms(g);

        for(const auto& model : g.models) {

            glBindTexture(GL_TEXTURE_2D, model.id_textura);

            glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, model.diff);
            glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, model.spec);
            glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, model.emiss);
            glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, model.amb);
            glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, model.shininess);

            glBindBuffer(GL_ARRAY_BUFFER, primitives.find(model.file)->second.coords[0]);
            glVertexPointer(3, GL_FLOAT, 0, 0);

            glBindBuffer(GL_ARRAY_BUFFER, primitives.find(model.file)->second.normals[0]);
            glNormalPointer(GL_FLOAT, 0, 0);

            glBindBuffer(GL_ARRAY_BUFFER,primitives.find(model.file)->second.texCoord[0]);
            glTexCoordPointer(2, GL_FLOAT,0,0);

            glDrawArrays(GL_TRIANGLES, 0, primitives.find(model.file)->second.num_vertices);

            glBindTexture(GL_TEXTURE_2D, 0);
        }

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

    int i = 0;
    for (Light &l : lightsV){
        apply_light(l, i);
        i++;
    }

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
    glutInitWindowSize(1400,1000);
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
    glEnableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);

    glEnable(GL_NORMALIZE);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_TEXTURE_2D);

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

