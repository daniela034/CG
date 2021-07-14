#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glew.h>
#include <GL/glut.h>
#endif

#include <vector>
#include <math.h>

using namespace std;

struct Point{
    float x;
    float y;
    float z;
};

void multMatrixVector(float *m, float *v, float *res) {

    for (int j = 0; j < 4; ++j) {
        res[j] = 0;
        for (int k = 0; k < 4; ++k) {
            res[j] += v[k] * m[j * 4 + k];
        }
    }
}

void getCatmullRomPoint(float real_t, Point p1, Point p2, Point p3, Point p4, float* pos, float* deriv){

    pos[0] = 0.0; pos[1] = 0.0; pos[2] = 0.0;
    deriv[0] = 0.0; deriv[1] = 0.0; deriv[2] = 0.0;

    float matrix[4][4] = {{-0.5f,  1.5f, -1.5f,  0.5f},
                          {1.0f, -2.5f,  2.0f, -0.5f},
                          {-0.5f,  0.0f,  0.5f,  0.0f},
                          {0.0f,  1.0f,  0.0f,  0.0f}};

    float T[4] = {(powf(real_t, 3)), powf(real_t, 2), real_t, 1};
    float T_deriv[4] = {(3 * powf(real_t, 2)), 2 * real_t, 1, 0};

    // Compute A = M * P
    float A[4][4] = {0};
    float PX[4] = {p1.x, p2.x, p3.x, p4.x};
    float PY[4] = {p1.y, p2.y, p3.y, p4.y};
    float PZ[4] = {p1.z, p2.z, p3.z, p4.z};

    multMatrixVector(*matrix, PX, A[0]);
    multMatrixVector(*matrix, PY, A[1]);
    multMatrixVector(*matrix, PZ, A[2]);

    // Compute pos = T * A
    pos[0] = A[0][0] * T[0] + A[0][1] * T[1] + A[0][2] * T[2] + A[0][3] * T[3];
    pos[1] = A[1][0] * T[0] + A[1][1] * T[1] + A[1][2] * T[2] + A[1][3] * T[3];
    pos[2] = A[2][0] * T[0] + A[2][1] * T[1] + A[2][2] * T[2] + A[2][3] * T[3];

    // compute deriv = T' * A
    deriv[0] = A[0][0] * T_deriv[0] + A[0][1] * T_deriv[1] + A[0][2] * T_deriv[2] + A[0][3] * T_deriv[3];
    deriv[1] = A[1][0] * T_deriv[0] + A[1][1] * T_deriv[1] + A[1][2] * T_deriv[2] + A[1][3] * T_deriv[3];
    deriv[2] = A[2][0] * T_deriv[0] + A[2][1] * T_deriv[1] + A[2][2] * T_deriv[2] + A[2][3] * T_deriv[3];

}

void getGlobalCatmullRomPoint(float global_t, vector<Point> points, float* pos, float* deriv){
    int num_pt = points.size();

    float real_t = global_t * num_pt;
    int indice = floor(real_t);
    real_t = real_t - indice;

    int indexs[4];
    indexs[0] = (indice + num_pt - 1) % num_pt;
    indexs[1] = (indexs[0] + 1) % num_pt;
    indexs[2] = (indexs[1] + 1) % num_pt;
    indexs[3] = (indexs[2] + 1) % num_pt;

    getCatmullRomPoint(real_t, points[indexs[0]], points[indexs[1]], points[indexs[2]], points[indexs[3]],
                pos, deriv);
}

void renderCatmullRomCurve(vector<Point> points){
    float pos[3], deriv[3];

    glDisable(GL_LIGHTING);
    glBegin(GL_LINE_LOOP);
    for (int i = 0; i < 1000; i++) {
        getGlobalCatmullRomPoint( (float) (i/1000.0), points, pos, deriv);
        glVertex3f(pos[0], pos[1], pos[2]);
    }
    glEnd();
    glEnable(GL_LIGHTING);
}

void buildRotMatrix(float *x, float *y, float *z, float *m) {

    m[0] = x[0]; m[1] = x[1]; m[2] = x[2]; m[3] = 0;
    m[4] = y[0]; m[5] = y[1]; m[6] = y[2]; m[7] = 0;
    m[8] = z[0]; m[9] = z[1]; m[10] = z[2]; m[11] = 0;
    m[12] = 0; m[13] = 0; m[14] = 0; m[15] = 1;
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


