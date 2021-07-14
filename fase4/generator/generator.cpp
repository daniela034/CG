#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <sstream>

#define _USE_MATH_DEFINES

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

struct Point {
    float x;
    float y;
    float z;
    float norm_x;
    float norm_y;
    float norm_z;
    float tex_u;
    float tex_v;
};

void point_writer(Point p, ofstream &file_pointer){
    string point = to_string(p.x) + " " + to_string(p.y) + " " + to_string(p.z) + " " + to_string(p.norm_x) + " " + to_string(p.norm_y) + " " + to_string(p.norm_z)
                    + " " + to_string(p.tex_u) + " " + to_string(p.tex_v);
    file_pointer.write(point.c_str(), point.size() + 1);
}

void plane(int size, char* file){
    float metade = (float) size/2;
    Point p1, p2, p3, p4, p5, p6;

    //Triângulo 1 (frente)
    p1.x = metade; p1.y = 0; p1.z = metade; p1.norm_x = 0; p1.norm_y = 1; p1.norm_z = 0; p1.tex_u = 1; p1.tex_v = 0;
    p2.x = -metade; p2.y = 0; p2.z = -metade; p2.norm_x = 0; p2.norm_y = 1; p2.norm_z = 0; p2.tex_u = 0; p2.tex_v = 1;
    p3.x = -metade; p3.y = 0; p3.z = metade; p3.norm_x = 0; p3.norm_y = 1; p3.norm_z = 0; p3.tex_u = 0; p3.tex_v = 0;

    //Triângulo 2 (frente)
    p4.x = metade; p4.y = 0; p4.z = metade; p4.norm_x = 0; p4.norm_y = 1; p4.norm_z = 0; p4.tex_u = 1; p4.tex_v = 0;
    p5.x = metade; p5.y = 0; p5.z = -metade; p5.norm_x = 0; p5.norm_y = 1; p5.norm_z = 0; p5.tex_u = 1; p5.tex_v = 1;
    p6.x = -metade; p6.y = 0; p6.z = -metade; p6.norm_x = 0; p6.norm_y = 1; p6.norm_z = 0; p6.tex_u = 0; p6.tex_v = 1;

    //Para a parte de trás dos triângulos é só escrever para o ficheiro na ordem contrária

    ofstream file_pointer(file, ios::binary | ios::out);

    //Triângulo 1 e 2 (frente)
    point_writer(p1, file_pointer);
    point_writer(p2, file_pointer);
    point_writer(p3, file_pointer);

    point_writer(p4, file_pointer);
    point_writer(p5, file_pointer);
    point_writer(p6, file_pointer);

    //Triângulo 1 e 2 (trás)
    p1.norm_y = p2.norm_y = p3.norm_y = p4.norm_y = p5.norm_y = p6.norm_y = -1;
    p1.tex_u = 0; p1.tex_v = 0;
    p2.tex_u = 1; p2.tex_v = 1;
    p3.tex_u = 1; p3.tex_v = 0;
    p4.tex_u = 0; p4.tex_v = 0;
    p5.tex_u = 0; p5.tex_v = 1;
    p6.tex_u = 1; p6.tex_v = 1;

    point_writer(p3, file_pointer);
    point_writer(p2, file_pointer);
    point_writer(p1, file_pointer);

    point_writer(p6, file_pointer);
    point_writer(p5, file_pointer);
    point_writer(p4, file_pointer);

    file_pointer.close();
}

float distance2p(Point p1, Point p2) {
    return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2) + pow(p1.z - p2.z, 2));
}

void divideTriangle(Point p1, Point p2, Point p3, int divide, ofstream& file_pointer) {
    if (divide == 0) {
        point_writer(p1, file_pointer);
        point_writer(p2, file_pointer);
        point_writer(p3, file_pointer);
    }
    else {
        if (distance2p(p1, p2) >= distance2p(p2, p3) && distance2p(p1, p2) >= distance2p(p1, p3)) {
            Point p4 = { (float)(p1.x + p2.x) / 2, (float)(p1.y + p2.y) / 2, (float)(p1.z + p2.z) / 2 };
            divideTriangle(p1, p4, p3, divide - 1, file_pointer);
            divideTriangle(p4, p2, p3, divide - 1, file_pointer);
        }
        else if (distance2p(p2, p3) >= distance2p(p1, p2) && distance2p(p2, p3) >= distance2p(p1, p3)) {
            Point p4 = { (float)(p2.x + p3.x) / 2, (float)(p2.y + p3.y) / 2, (float)(p2.z + p3.z) / 2 };
            divideTriangle(p1, p2, p4, divide - 1, file_pointer);
            divideTriangle(p1, p4, p3, divide - 1, file_pointer);
        }
        else {
            Point p4 = { (float)(p1.x + p3.x) / 2, (float)(p1.y + p3.y) / 2, (float)(p1.z + p3.z) / 2 };
            divideTriangle(p1, p2, p4, divide - 1, file_pointer);
            divideTriangle(p2, p3, p4, divide - 1, file_pointer);
        }
    }
}


void box (float size_x, float size_y, float size_z, int divide, char* file){
    float mX = (float) size_x/2;
    float mY = (float) size_y/2;
    float mZ = (float) size_z/2;
    Point p1, p2, p3, p4;

    ofstream file_pointer(file, ios::binary | ios::out);

    /********** FACE FRENTE ***********/

    p1 = {-mX, mY, mZ, 0, 0, 1, 0, (float) 1/2};
    p2 = {-mX, -mY, mZ, 0, 0, 1, 0, 0};
    p3 = {mX, -mY, mZ, 0, 0, 1, (float) 1/3, 0};
    p4 = {mX, mY, mZ, 0, 0, 1, (float) 1/3, (float) 1/2};

    /*
    //Triangulo 1
    divideTriangle(p1, p2, p3, divide, file_pointer);

    //Triângulo 2
    divideTriangle(p3, p4, p1, divide, file_pointer);
    */

    point_writer(p1, file_pointer);
    point_writer(p2, file_pointer);
    point_writer(p3, file_pointer);

    point_writer(p3, file_pointer);
    point_writer(p4, file_pointer);
    point_writer(p1, file_pointer);

    /********** FACE TRÁS ***********/

    p1.z = p2.z = p3.z = p4.z = -mZ;
    p1.norm_z = p2.norm_z = p3.norm_z = p4.norm_z = -1;

    p1.tex_u = 0; p1.tex_v = 1;
    p2.tex_u = 0; p2.tex_v = (float) 1/2;
    p3.tex_u = (float) 1/3; p3.tex_v = (float) 1/2;
    p4.tex_u = (float) 1/3; p4.tex_v = 1;

    /*
    //Triangulo 1
    divideTriangle(p3, p2, p1, divide, file_pointer);

    //Triangulo 2
    divideTriangle(p1, p4, p3, divide, file_pointer);
    */

    point_writer(p3, file_pointer);
    point_writer(p2, file_pointer);
    point_writer(p1, file_pointer);

    point_writer(p1, file_pointer);
    point_writer(p4, file_pointer);
    point_writer(p3, file_pointer);

    /********** FACE CIMA ***********/
    p1 = {-mX, mY, mZ, 0, 1, 0, (float) 1/3, (float) 1/2};
    p2 = {mX, mY, mZ, 0, 1, 0, (float) 1/3, 0};
    p3 = {mX, mY, -mZ, 0, 1, 0, (float) 2/3, 0};
    p4 = {-mX, mY, -mZ, 0, 1, 0, (float) 2/3, (float) 1/2};

    /*
    //Triangulo 1
    divideTriangle(p1, p2, p4, divide, file_pointer);

    //Triangulo 2
    divideTriangle(p2, p3, p4, divide, file_pointer);
    */

    point_writer(p1, file_pointer);
    point_writer(p2, file_pointer);
    point_writer(p4, file_pointer);

    point_writer(p2, file_pointer);
    point_writer(p3, file_pointer);
    point_writer(p4, file_pointer);

    /********** FACE BAIXO ***********/
    p1.y = p2.y = p3.y = p4.y = -mY;
    p1.norm_y = p2.norm_y = p3.norm_y = p4.norm_y = -1;

    p1.tex_u = (float) 1/3; p1.tex_v = 1;
    p2.tex_u = (float) 1/3; p2.tex_v = (float) 1/2;
    p3.tex_u = (float) 2/3; p3.tex_v = (float) 1/2;
    p4.tex_u = (float) 2/3; p4.tex_v = 1;

    /*
    //Triangulo 1
    divideTriangle(p1, p4, p2, divide, file_pointer);

    //Triangulo 2
    divideTriangle(p4, p3, p2, divide, file_pointer);
    */

    point_writer(p1, file_pointer);
    point_writer(p4, file_pointer);
    point_writer(p2, file_pointer);

    point_writer(p4, file_pointer);
    point_writer(p3, file_pointer);
    point_writer(p2, file_pointer);

    /********** FACE DIREITA ***********/
    p1 = {mX, -mY, mZ, 1, 0, 0, (float) 2/3, (float) 1/2};
    p2 = {mX, -mY, -mZ, 1, 0, 0, (float) 2/3, 0};
    p3 = {mX, mY, -mZ, 1, 0, 0, 1, 0};
    p4 = {mX, mY, mZ, 1, 0, 0, 1, (float) 1/2};

    /*
    //Triangulo 1
    divideTriangle(p1, p2, p3, divide, file_pointer);

    //Triangulo 2
    divideTriangle(p3, p4, p1, divide, file_pointer);
    */

    point_writer(p1, file_pointer);
    point_writer(p2, file_pointer);
    point_writer(p3, file_pointer);

    point_writer(p3, file_pointer);
    point_writer(p4, file_pointer);
    point_writer(p1, file_pointer);

    /********** FACE ESQUERDA ***********/
    p1.x = p2.x = p3.x = p4.x = -mX;
    p1.norm_x = p2.norm_x = p3.norm_x = p4.norm_x = -1;

    p1.tex_u = (float) 2/3; p1.tex_v = 1;
    p2.tex_u = (float) 2/3; p2.tex_v = (float) 1/2;
    p3.tex_u = 1; p3.tex_v = (float) 1/2;
    p4.tex_u = 1; p4.tex_v = 1;

    /*
    //Triangulo 1
    divideTriangle(p4, p3, p1, divide, file_pointer);

    //Triangulo 2
    divideTriangle(p3, p2, p1, divide, file_pointer);
     */

    point_writer(p4, file_pointer);
    point_writer(p3, file_pointer);
    point_writer(p1, file_pointer);

    point_writer(p3, file_pointer);
    point_writer(p2, file_pointer);
    point_writer(p1, file_pointer);

    file_pointer.close();
}

void sphere(float radius, int slices, int stacks, char* file){
    float alpha = 0;
    float beta = (float) - M_PI / 2;
    float alpha_inc = (float) (2 * M_PI) / (float) slices;
    float beta_inc = (float) (M_PI) / (float) stacks;

    Point p1, p2, p3, p4;

    ofstream file_pointer(file, ios::binary | ios::out);

    for (int i = 1; i < slices + 1; i++){
        for (int j = 1; j < stacks + 1; j++){

            p1 = {radius * cos(beta) * sin(alpha), radius * sin(beta), radius * cos(beta) * cos(alpha),
                  cos(beta) * sin(alpha), sin(beta), cos(beta) * cos(alpha),
                   ((float) (i-1) / (float) slices), (float) (j-1) / (float) stacks};

            p2 = {radius * cos(beta) * sin(alpha + alpha_inc), radius * sin(beta), radius * cos(beta) * cos(alpha + alpha_inc),
                  cos(beta) * sin(alpha + alpha_inc), sin(beta), cos(beta) * cos(alpha + alpha_inc),
                  ((float) (i) / (float) slices), (float) (j-1) / (float) stacks};

            p3 = {radius * cos(beta + beta_inc) * sin(alpha + alpha_inc), radius * sin(beta + beta_inc), radius * cos(beta + beta_inc) * cos(alpha + alpha_inc),
                  cos(beta + beta_inc) * sin(alpha + alpha_inc), sin(beta + beta_inc), cos(beta + beta_inc) * cos(alpha + alpha_inc),
                  ((float) (i) / (float) slices), (float) (j) / (float) stacks};

            p4 = {radius * cos(beta + beta_inc) * sin(alpha), radius * sin(beta + beta_inc), radius * cos(beta + beta_inc) * cos(alpha),
                  cos(beta + beta_inc) * sin(alpha), sin(beta + beta_inc), cos(beta + beta_inc) * cos(alpha),
                  ((float) (i-1) / (float) slices), (float) (j) / (float) stacks};

            if (j == 1){
                point_writer(p3, file_pointer);
                point_writer(p4, file_pointer);
                point_writer(p1, file_pointer);
            }
            else if (j == stacks){
                point_writer(p2, file_pointer);
                point_writer(p4, file_pointer);
                point_writer(p1, file_pointer);
            }
            else {
                point_writer(p1, file_pointer);
                point_writer(p2, file_pointer);
                point_writer(p3, file_pointer);

                point_writer(p3, file_pointer);
                point_writer(p4, file_pointer);
                point_writer(p1, file_pointer);
            }
            beta = (float) -(M_PI / 2) + j * beta_inc;
        }
        alpha = i * alpha_inc;
        beta = (float) -(M_PI / 2);
    }
    file_pointer.close();
}

void cone(float radius, float height, int slices, int stacks, char* file) {
    float alpha = 0;
    float alpha_inc = (float) (2 * M_PI) / (float) slices;
    float h = 0;
    float h_inc =  (float) height / (float) stacks;
    float r = radius;
    float new_r;

    Point p1, p2, p3, p4;

    ofstream file_pointer(file, ios::binary | ios::out);

    for (int i = 1; i < slices + 1; i++){
        p1 = {radius * sin(alpha + alpha_inc), 0, radius * cos(alpha + alpha_inc),
              0, -1, 0,
              (float) (cos(alpha + alpha_inc) * 0.5 + 0.5), (float) (sin(alpha + alpha_inc) * 0.5 + 0.5)};

        p2 = {radius * sin(alpha), 0, radius * cos(alpha),
              0, -1, 0,
              (float) (cos(alpha) * 0.5 + 0.5), (float) (sin(alpha) * 0.5 + 0.5)};

        p3 = {0, 0, 0, 0, -1, 0, 0.5, 0.5};

        point_writer(p1, file_pointer);
        point_writer(p2, file_pointer);
        point_writer(p3, file_pointer);

        for (int j = 1; j < stacks + 1; j++){

            new_r = ((height - (h + h_inc)) * radius) / height;  //semelhança de triângulos para calcular o novo raio

            p1 = {r * sin(alpha), h, r * cos(alpha),
                  r * sin(alpha), sin(tan(radius/height)), r * cos(alpha),
                  (float) ( ((float) (j-1)/ (float) stacks) * cos(alpha) + 0.5), (float)( ((float) (j-1)/ (float) stacks) * sin(alpha) + 0.5) };

            p2 = {r * sin(alpha + alpha_inc), h, r * cos(alpha + alpha_inc),
                  r * sin(alpha + alpha_inc), sin(tan(radius/height)), r * cos(alpha + alpha_inc),
                  (float) ( ((float) (j-1)/ (float) stacks) * cos(alpha + alpha_inc) + 0.5), (float)( ((float) (j-1)/ (float) stacks) * sin(alpha + alpha_inc) + 0.5)};

            p3 = {new_r * sin(alpha + alpha_inc), h + h_inc, new_r * cos(alpha + alpha_inc),
                  new_r * sin(alpha + alpha_inc), sin(tan(radius/height)), new_r * cos(alpha + alpha_inc),
                  (float) ( ((float) j/ (float) stacks) * cos(alpha + alpha_inc) + 0.5), (float)( ((float) j/ (float) stacks) * sin(alpha + alpha_inc) + 0.5)};

            p4 = {new_r * sin(alpha), h + h_inc, new_r * cos(alpha),
                  new_r * sin(alpha), sin(tan(radius/height)), new_r * cos(alpha + alpha_inc),
                  (float) ( ((float) j/ (float) stacks) * cos(alpha) + 0.5), (float)( ((float) j/ (float) stacks) * sin(alpha) + 0.5)};

            point_writer(p1, file_pointer);
            point_writer(p2, file_pointer);
            point_writer(p3, file_pointer);

            point_writer(p3, file_pointer);
            point_writer(p4, file_pointer);
            point_writer(p1, file_pointer);

            h = j * h_inc;
            r = new_r;
        }
        alpha = i * alpha_inc;
        h = 0;
        r = radius;
    }
    file_pointer.close();
}

void ring(float in_radius, float out_radius, int slices, char* file){
    float alpha = 0;
    float alpha_inc = (float) (2 * M_PI) / (float) slices;

    Point p1, p2, p3, p4;

    ofstream file_pointer(file, ios::binary | ios::out);

    for (int i = 1; i < slices + 1; i++){
        p1 = {in_radius * sin(alpha), 0, in_radius * cos(alpha), 0, 1, 0, 0, 0};
        p2 = {out_radius * sin(alpha), 0, out_radius * cos(alpha), 0, 1, 0, 1, 0};
        p3 = {out_radius * sin(alpha + alpha_inc), 0, out_radius * cos(alpha + alpha_inc), 0, 1, 0, 1, 1};
        p4 = {in_radius * sin(alpha + alpha_inc), 0, in_radius * cos(alpha + alpha_inc), 0, 1, 0, 0, 1};

        /* Face de cima */
        point_writer(p1, file_pointer);
        point_writer(p2, file_pointer);
        point_writer(p4, file_pointer);

        point_writer(p3, file_pointer);
        point_writer(p4, file_pointer);
        point_writer(p2, file_pointer);

        p1.norm_y = p2.norm_y = p3.norm_y = p4.norm_y = -1;

        p1.tex_u = 0; p1.tex_v = 1;
        p2.tex_u = 1; p2.tex_v = 1;
        p3.tex_u = 1; p3.tex_v = 0;
        p4.tex_u = 0; p4.tex_v = 0;

        /* Face de baixo */
        point_writer(p2, file_pointer);
        point_writer(p1, file_pointer);
        point_writer(p3, file_pointer);

        point_writer(p3, file_pointer);
        point_writer(p1, file_pointer);
        point_writer(p4, file_pointer);

        alpha = i * alpha_inc;
    }
    file_pointer.close();
}

int ler_bezier_file(char* file, vector<vector<float>> &patches, vector<float> &points){
    int num_patches = 0, num_control_points = 0;
    string line, control_point, point;

    ifstream file_pointer (file);

    if (getline(file_pointer, line)) num_patches = stoi(line);

    for (int i = 0; i < num_patches; i++){
        if (getline(file_pointer, line)){
            stringstream ss(line);
            vector<float> tokens;

            while(getline(ss, control_point, ',')){
                tokens.push_back(stof(control_point));
            }

            patches.push_back(tokens);   //guardar um vetor correspondente aos pontos de controlo de um patch no vetor de vetores (patches)
        }
    }

    if (getline(file_pointer, line)) num_control_points = stoi(line);

    for(int i = 0; i < num_control_points; i++){
        if (getline(file_pointer, line)){
            stringstream ss(line);

            while(getline(ss, point, ',')){
                points.push_back(stof(point));
            }
        }
    }

    file_pointer.close();

    return num_patches;
}

float calcular_coords(float u_step, float v_step, float coords[4][4], int opcao){
    float m[4][4]={{-1, 3, -3, 1},
                   {3, -6, 3, 0},
                   {-3, 3, 0, 0},
                   {1, 0, 0, 0}};

    float u_step3 = pow(u_step, 3);
    float u_step2 = pow(u_step, 2);
    float v_step3 = pow(v_step, 3);
    float v_step2 = pow(v_step, 2);

    /* Para as derivadas de u e v (4 fase) */
    float u_step3_d = 3 * u_step2;
    float v_step3_d = 3 * v_step2;

    float vector_u[4] = {u_step3, u_step2, u_step, 1};
    float vector_v[4] = {v_step3, v_step2, v_step, 1};

    /* Para as derivadas de u e v (4 fase) */
    float vector_u_d[4] = {u_step3_d, 2 * u_step, 1, 0};
    float vector_v_d[4] = {v_step3_d, 2 * v_step, 1, 0};

    float res1[4] = {0, 0, 0, 0};
    float res2[4] = {0, 0, 0, 0};
    float res3[4] = {0, 0, 0, 0};
    float res_final = 0;

    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            if (opcao == 0 || opcao == 2) {
                res1[i] += vector_u[j] * m[j][i];
            }
            else if (opcao == 1) {
                res1[i] += vector_u_d[j] * m[j][i];
            }
        }
    }

    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            res2[i] += res1[j] * coords[j][i];
        }
    }

    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            res3[i] += res2[j] * m[j][i];
        }
    }

    for (int i = 0; i < 4; i++){
        if (opcao == 0 || opcao == 1) {
            res_final += res3[i] * vector_v[i];
        }
        else if (opcao == 2) {
            res_final += res3[i] * vector_v_d[i];
        }
    }

    return res_final;
}

void cross(float *a, float *b, float *res) {

    res[0] = a[1]*b[2] - a[2]*b[1];
    res[1] = a[2]*b[0] - a[0]*b[2];
    res[2] = a[0]*b[1] - a[1]*b[0];
}

void normalize(float *a) {

    float l = sqrt(a[0]*a[0] + a[1] * a[1] + a[2] * a[2]);
    a[0] = a[0]/l;
    a[1] = a[1]/l;
    a[2] = a[2]/l;
}

void bezier(char* in_file, int tess, char* out_file){
    int indice;
    float x[4][4], y[4][4], z[4][4];

    float steps = (float) (1) / (float) tess;
    float u_step, v_step;

    vector<vector<float>> patches;
    vector<float> points;

    int num_patches = ler_bezier_file(in_file, patches, points);

    Point p1, p2, p3, p4;

    ofstream file_pointer(out_file, ios::binary | ios::out);

    float *cross_prod = (float *) malloc(sizeof(float) * 3);
    float *der_U = (float *) malloc(sizeof(float) * 3);
    float *der_V = (float *) malloc(sizeof(float) * 3);

    float norm_x, norm_y, norm_z;

    for (int i = 0; i < num_patches; i++){

        //Para cada patch guardar os valores de x, y e z em matrizes de 4x4 (já que o nº de pontos de controlo de um patch é 16)
        for (int j = 0; j < 4; j++){
            for (int k = 0; k < 4; k++){
                //Vai-se buscar o indice à matriz patches
                indice = patches[i][k + (4 * j)];

                //Vai-se buscar os valores de x, y e z ao vetor points (onde estão guardadas as coordenadas de todos os pontos de controlo)
                x[j][k] = points[indice * 3];
                y[j][k] = points[(indice * 3) + 1];
                z[j][k] = points[(indice * 3) + 2];
            }
        }

        //Calcular 4 pontos para cada iteração através das fórmulas específicas que vão ser usados para desenhar um quadrado (2 triângulos)
        for (int u = 0; u < tess; u++){
            u_step = u * steps;
            for (int v = 0; v < tess; v++){
                v_step = v * steps;

                der_U[0] = calcular_coords(u_step, v_step, x, 1);
                der_U[1] = calcular_coords(u_step, v_step, y, 1);
                der_U[2] = calcular_coords(u_step, v_step, z, 1);

                der_V[0] = calcular_coords(u_step, v_step, x, 2);
                der_V[1] = calcular_coords(u_step, v_step, y, 2);
                der_V[2] = calcular_coords(u_step, v_step, y, 2);

                normalize(der_U);
                normalize(der_V);
                cross(der_U, der_V, cross_prod);
                normalize(cross_prod);

                norm_x = -1 * cross_prod[0];
                norm_y = -1 * cross_prod[1];
                norm_z = -1 * cross_prod[2];

                p1 = {calcular_coords(u_step, v_step, x, 0),
                      calcular_coords(u_step, v_step, y, 0),
                      calcular_coords(u_step, v_step, z, 0), norm_x, norm_y, norm_z, v_step, u_step};

                der_U[0] = calcular_coords(u_step, v_step + steps, x, 1);
                der_U[1] = calcular_coords(u_step, v_step + steps, y, 1);
                der_U[2] = calcular_coords(u_step, v_step + steps, z, 1);

                der_V[0] = calcular_coords(u_step, v_step + steps, x, 2);
                der_V[1] = calcular_coords(u_step, v_step + steps, y, 2);
                der_V[2] = calcular_coords(u_step, v_step + steps, z, 2);

                normalize(der_U);
                normalize(der_V);
                cross(der_U, der_V, cross_prod);
                normalize(cross_prod);

                norm_x = -1 * cross_prod[0];
                norm_y = -1 * cross_prod[1];
                norm_z = -1 * cross_prod[2];

                p2 = {calcular_coords(u_step, v_step + steps, x, 0),
                      calcular_coords(u_step, v_step + steps, y, 0),
                      calcular_coords(u_step, v_step + steps, z, 0), norm_x, norm_y, norm_z, v_step + steps, u_step};

                der_U[0] = calcular_coords(u_step + steps, v_step, x, 1);
                der_U[1] = calcular_coords(u_step + steps, v_step, y, 1);
                der_U[2] = calcular_coords(u_step + steps, v_step, z, 1);

                der_V[0] = calcular_coords(u_step + steps, v_step, x, 2);
                der_V[1] = calcular_coords(u_step + steps, v_step, y, 2);
                der_V[2] = calcular_coords(u_step + steps, v_step, z, 2);

                normalize(der_U);
                normalize(der_V);
                cross(der_U, der_V, cross_prod);
                normalize(cross_prod);

                norm_x = -1 * cross_prod[0];
                norm_y = -1 * cross_prod[1];
                norm_z = -1 * cross_prod[2];

                p3 = {calcular_coords(u_step + steps, v_step, x, 0),
                      calcular_coords(u_step + steps, v_step, y, 0),
                      calcular_coords(u_step + steps, v_step, z, 0), norm_x, norm_y, norm_z, v_step, u_step + steps};

                der_U[0] = calcular_coords(u_step + steps, v_step + steps, x, 1);
                der_U[1] = calcular_coords(u_step + steps, v_step + steps, y, 1);
                der_U[2] = calcular_coords(u_step + steps, v_step + steps, z, 1);

                der_V[0] = calcular_coords(u_step + steps, v_step + steps, x, 2);
                der_V[1] = calcular_coords(u_step + steps, v_step + steps, y, 2);
                der_V[2] = calcular_coords(u_step + steps, v_step + steps, z, 2);

                normalize(der_U);
                normalize(der_V);
                cross(der_U, der_V, cross_prod);
                normalize(cross_prod);

                norm_x = -1 * cross_prod[0];
                norm_y = -1 * cross_prod[1];
                norm_z = -1 * cross_prod[2];

                p4 = {calcular_coords(u_step + steps, v_step + steps, x, 0),
                      calcular_coords(u_step + steps, v_step + steps, y, 0),
                      calcular_coords(u_step + steps, v_step + steps, z, 0), norm_x, norm_y, norm_z, v_step + steps, u_step + steps};

                point_writer(p1, file_pointer);
                point_writer(p2, file_pointer);
                point_writer(p3, file_pointer);

                point_writer(p2, file_pointer);
                point_writer(p4, file_pointer);
                point_writer(p3, file_pointer);
            }
        }
    }
}

void opcao(int argc, char* argv[], int primitive){
    switch (primitive) {
        case 1:
            if (argc <= 3){
                std::cout << "Nº de argumentos insuficientes para o plano" << std::endl;
            }
            else {
                int size = stoi(argv[2]);
                char *file = argv[3];
                plane(size, file);
            }
            break;

        case 2:
            if (argc <= 5){
                std::cout << "Nº de argumentos insuficientes para o cubo" << std::endl;
            }
            else {
                float size_x = stof(argv[2]);
                float size_y = stof(argv[3]);
                float size_z = stof(argv[4]);
                if (argc == 6) {
                    char *file = argv[5];
                    box(size_x, size_y, size_z, 0, file);
                } else {
                    int divide = stof(argv[5]);
                    char *file = argv[6];
                    box(size_x, size_y, size_z, divide, file);
                }
            }
            break;

        case 3:
            if(argc <= 5){
                std::cout << "Nº de argumentos insuficientes para a esfera" << std::endl;
            }
            else {
                float radius = stof(argv[2]);
                int slices = stoi(argv[3]);
                int stacks = stoi(argv[4]);
                char* file = argv[5];
                sphere(radius, slices, stacks, file);
            }
            break;

        case 4:
            if (argc <= 6) {
                std::cout << "Nº de argumentos insuficientes para o cone" << std::endl;
            }
            else {
                float radius = stof(argv[2]);
                float height = stof(argv[3]);
                int slices = stoi(argv[4]);
                int stacks = stoi(argv[5]);
                char* file = argv[6];
                cone(radius, height, slices, stacks, file);
            }
            break;

        case 5:
            if (argc <= 5) {
                std::cout << "Nº de argumentos insuficientes para o anel" << std::endl;
            }
            else {
                float in_radius = stof(argv[2]);
                float out_radius = stof(argv[3]);
                int slices = stoi(argv[4]);
                char* file = argv[5];
                ring(in_radius, out_radius, slices, file);
            }
            break;

        case 6:
            if (argc <= 4) {
                std::cout << "Nº de argumentos insuficientes para o bezier" << std::endl;
            }
            else {
                char* in_file = argv[2];
                int tess = stoi(argv[3]);
                char* out_file = argv[4];
                bezier(in_file, tess, out_file);
            }

            break;
    }
}

int main(int argc, char* argv[]) {

    if (argc <= 2){
        std::cout << "Nº de argumentos insuficientes!" << std::endl;
    }

    else {
        int primitive = 0;
        string prim = argv[1];

        if (prim == "plane") primitive = 1;
        else if (prim == "box") primitive = 2;
        else if (prim == "sphere") primitive = 3;
        else if (prim == "cone") primitive = 4;
        else if (prim == "ring") primitive = 5;
        else if (prim == "bezier") primitive = 6;
        else std::cout << "Primitiva inexistente!" << std::endl;

        opcao(argc, argv, primitive);
    }
    return 0;
}
