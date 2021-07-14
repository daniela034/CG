#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#define _USE_MATH_DEFINES

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

struct Point {
    float x;
    float y;
    float z;
};

void point_writer(Point p, ofstream &file_pointer){
    string point = to_string(p.x) + " " + to_string(p.y) + " " + to_string(p.z);
    file_pointer.write(point.c_str(), point.size() + 1);
}

void plane(int size, char* file){
    float metade = (float) size/2;
    Point p1, p2, p3, p4, p5, p6;

    //Triângulo 1 (frente)
    p1.x = metade; p1.y = 0; p1.z = metade;
    p2.x = -metade; p2.y = 0; p2.z = -metade;
    p3.x = -metade; p3.y = 0; p3.z = metade;

    //Triângulo 2 (frente)
    p4.x = metade; p4.y = 0; p4.z = metade;
    p5.x = metade; p5.y = 0; p5.z = -metade;
    p6.x = -metade; p6.y = 0; p6.z = -metade;

    //Para a parte de trás dos triângulos é só escrever para o ficheiro na ordem contrária

    ofstream file_pointer(file, ios::binary | ios::out);

    //Triângulo 1 (frente e trás)
    point_writer(p1, file_pointer);
    point_writer(p2, file_pointer);
    point_writer(p3, file_pointer);

    point_writer(p3, file_pointer);
    point_writer(p2, file_pointer);
    point_writer(p1, file_pointer);

    //Triângulo 2 (frente e trás)
    point_writer(p4, file_pointer);
    point_writer(p5, file_pointer);
    point_writer(p6, file_pointer);

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

    p1 = {-mX, mY, mZ};
    p2 = {-mX, -mY, mZ};
    p3 = {mX, -mY, mZ};
    p4 = {mX, mY, mZ};

    //Triangulo 1
    divideTriangle(p1, p2, p3, divide, file_pointer);

    //Triângulo 2
    divideTriangle(p3, p4, p1, divide, file_pointer);


    /********** FACE TRÁS ***********/

    p1.z = p2.z = p3.z = p4.z = -mZ;

    //Triangulo 1
    divideTriangle(p3, p2, p1, divide, file_pointer);

    //Triangulo 2
    divideTriangle(p1, p4, p3, divide, file_pointer);

    /********** FACE CIMA ***********/
    p1 = {-mX, mY, mZ};
    p2 = {mX, mY, mZ};
    p3 = {mX, mY, -mZ};
    p4 = {-mX, mY, -mZ};

    //Triangulo 1
    divideTriangle(p1, p2, p4, divide, file_pointer);

    //Triangulo 2
    divideTriangle(p2, p3, p4, divide, file_pointer);

    /********** FACE BAIXO ***********/
    p1.y = p2.y = p3.y = p4.y = -mY;

    //Triangulo 1
    divideTriangle(p1, p4, p2, divide, file_pointer);

    //Triangulo 2
    divideTriangle(p4, p3, p2, divide, file_pointer);

    /********** FACE DIREITA ***********/
    p1 = {mX, -mY, mZ};
    p2 = {mX, -mY, -mZ};
    p3 = {mX, mY, -mZ};
    p4 = {mX, mY, mZ};

    //Triangulo 1
    divideTriangle(p1, p2, p3, divide, file_pointer);

    //Triangulo 2
    divideTriangle(p3, p4, p1, divide, file_pointer);

    /********** FACE ESQUERDA ***********/
    p1.x = p2.x = p3.x = p4.x = -mX;

    //Triangulo 1
    divideTriangle(p4, p3, p1, divide, file_pointer);

    //Triangulo 2
    divideTriangle(p3, p2, p1, divide, file_pointer);

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

            p1 = {radius * cos(beta) * sin(alpha), radius * sin(beta), radius * cos(beta) * cos(alpha)};
            p2 = {radius * cos(beta) * sin(alpha + alpha_inc), radius * sin(beta), radius * cos(beta) * cos(alpha + alpha_inc)};
            p3 = {radius * cos(beta + beta_inc) * sin(alpha + alpha_inc), radius * sin(beta + beta_inc), radius * cos(beta + beta_inc) * cos(alpha + alpha_inc)};
            p4 = {radius * cos(beta + beta_inc) * sin(alpha), radius * sin(beta + beta_inc), radius * cos(beta + beta_inc) * cos(alpha)};

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
        p1 = {radius * sin(alpha + alpha_inc), 0, radius * cos(alpha + alpha_inc)};
        p2 = {radius * sin(alpha), 0, radius * cos(alpha)};
        p3 = {0, 0, 0};

        point_writer(p1, file_pointer);
        point_writer(p2, file_pointer);
        point_writer(p3, file_pointer);

        for (int j = 1; j < stacks + 1; j++){

            new_r = ((height - (h + h_inc)) * radius) / height;  //semelhança de triângulos para calcular o novo raio

            p1 = {r * sin(alpha), h, r * cos(alpha)};
            p2 = {r * sin(alpha + alpha_inc), h, r * cos(alpha + alpha_inc)};
            p3 = {new_r * sin(alpha + alpha_inc), h + h_inc, new_r * cos(alpha + alpha_inc)};
            p4 = {new_r * sin(alpha), h + h_inc, new_r * cos(alpha)};

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
    }
}

int main(int argc, char* argv[]) {

    if (argc <= 2){
        std::cout << "Nº de argumentos insuficientes!" << std::endl;
    }

    else {
        int primitive = 0;
        string prim = argv[1];

        if (prim.compare("plane") == 0) primitive = 1;
        else if (prim.compare("box") == 0) primitive = 2;
        else if (prim.compare("sphere") == 0) primitive = 3;
        else if (prim.compare("cone") == 0) primitive = 4;
        else std::cout << "Primitiva inexistente!" << std::endl;

        opcao(argc, argv, primitive);
    }
    return 0;
}
