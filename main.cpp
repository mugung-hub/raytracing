#include <iostream>
#include <cmath>

double M = 10;
double r = 100;
double eps = 0.001;
double dt = 0.1;

double clip(double x, double min, double max){
    if (x>max){
        x = max;
    }
    if (x < min){
        x = min;
    }
    return x;
}


class ray{
    public:
    double pos[2];
    double vel[3];
    double energy, angular;

    ray(double x, double y, double z){
        pos[0] = r;
        pos[1] = 0;
        double length = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
        vel[0] = x/length;
        vel[1] = y/length;
        vel[2] = z/length;
        
        double dr_dtau = vel[2];
        double dtheta_dtau = sqrt(pow(vel[0],2) + pow(vel[1],2))/r;
        double dt_dtau = 1/(1-2*M/r)*(pow(dr_dtau,2)/(1-2*M/r) + pow(r,2)*pow(dtheta_dtau,2));
        dt_dtau = sqrt(dt_dtau);

        energy = (1-2*M/r)*dt_dtau;
        angular = pow(r,2)*dtheta_dtau;

    };

    void dX_dtau(double *result);
    
    bool check_hit();
    bool infall = true;



};

void ray::dX_dtau(double *result){
    double dtheta_dtau = angular / pow(pos[0],2);
    double dr_dtau_square = pow(energy,2) - (1-2*M/pos[0])*pow(angular/pos[0],2);

    if (infall){
        if (dr_dtau_square / (1-2*M/pos[0]) < eps ){
            infall = false;
        }
    }
    result[1] = dtheta_dtau;
    if (infall){
        result[0] = - sqrt(dr_dtau_square);
    }
    else {
        result[0] = sqrt(dr_dtau_square);
    }
}

bool ray::check_hit(){
    bool hit = false;
    while (pos[0] < r+2 && pos[0] > 2*M+eps){
        double dX[2];
        dX_dtau(dX);
        pos[0] += dt*dX[0];
        pos[1] += dt*dX[1];
        if (pos[0]*cos(pos[1]) < -r + 2*eps){
            double X = pos[0]*cos(pos[1]);
            double Y = pos[0]*sin(pos[1]);
            double distance = pow(r+X,2)+pow(Y,2);
            distance = sqrt(distance);

            if (distance < 10){
                hit = true;
                return hit;
                break;

            }
        }

    }
    return hit;

}

int main(){
    using namespace std;
    int i,j;
    int vertical = 240;
    int horizontal = (int) vertical*16/9;
    double z = (double)vertical*400/1080;
    double **screen = new double*[vertical];
    for (i=0;i<vertical;i++){
        screen[i] = new double[horizontal];
    }

    for (i=0;i<vertical;i++){
        for (j=0;j<horizontal;j++){
            screen[i][j] = 0;
        }
    }

    for (i = 0;i<vertical;i++){
        for (j=0;j<horizontal;j++){
            ray Ray = ray((double)(j-horizontal/2),(double)(i-vertical/2),z);
            bool check = Ray.check_hit();
            if (check){
                screen[i][j] += 0.1;
            }

        }
    }

    for (i=0;i<vertical;i++){
        for (j = 0;j<horizontal;j++){
            double u = clip(screen[i][j],0,1);
            cout << u << ',';
        }
        cout << endl;
    }

    delete(screen);
    return 0;


}