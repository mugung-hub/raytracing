#include <iostream>
#include <fstream>
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
    double phi[2];
    double energy, angular;

    ray(double x, double y, double z){
        pos[0] = r;
        pos[1] = 0;
        double length = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
        vel[0] = x/length;
        vel[1] = y/length;
        vel[2] = z/length;
        if (x == 0 && y == 0 ){
            phi[0]=0;
            phi[1]=0;
        }
        else{
            phi[0] = x/sqrt(pow(x,2)+pow(y,2));
            phi[1] = y/sqrt(pow(x,2) + pow(y,2));
        }
        
        double dr_dtau = vel[2];
        double dtheta_dtau = sqrt(pow(vel[0],2) + pow(vel[1],2))/r;
        double dt_dtau = 1/(1-2*M/r)*(pow(dr_dtau,2)/(1-2*M/r) + pow(r,2)*pow(dtheta_dtau,2));
        dt_dtau = sqrt(dt_dtau);

        energy = (1-2*M/r)*dt_dtau;
        angular = pow(r,2)*dtheta_dtau;

    };

    void dX_dtau(double *result);
    
    bool check_hit(double star1, double star2, double star3, double cos_0);
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

bool ray::check_hit(double star1, double star2, double star3, double cos_0){
    bool hit = false;
    while (pos[0]*cos(pos[1]) > -r-10 && pos[0] > 2*M+eps){
        double dX[2];
        dX_dtau(dX);
        pos[0] += dt*dX[0];
        pos[1] += dt*dX[1];
        if (pos[0]*cos(pos[1]) < -r + 2*eps){
            double Z1 = pos[0]*cos(pos[1]);
            double R1 = pos[0]*sin(pos[1]);
            double X1 = R1*phi[0];
            double Y1 = R1*phi[1];

            double Z2 = -r;
            double X2 = star1;
            double Y2 = star2;

            if (sqrt(pow(X1-X2,2)+pow(Y1-Y2,2)+pow(Z1-Z2,2)) < star3){
                hit = true;
                return hit;
                break;

            }
        }

    }
    return hit;

}

double random_double() {
    // Returns a random real in [0,1).
    return rand() / (RAND_MAX + 1.0);
}

class ray1{
    public:
    double pos[3];
    double vel[3];

    ray1(double x, double y, double z){
        pos[0] = 0;
        pos[1] = 0;
        pos[2] = r;
        double length = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
        vel[0] = x/length;
        vel[1] = y/length;
        vel[2] = -z/length;
    };

    bool check_hit(double star1, double star2, double star3, double cos_0);

};


bool ray1::check_hit(double star1, double star2, double star3, double cos_0){
    bool hit = false;
    double rr = r;
    while ( pos[2] > -r-5 && rr > 2*M+eps){
        pos[0] += dt*vel[0];
        pos[1] += dt*vel[1];
        pos[2] += dt*vel[2];

        rr = sqrt(pow(pos[0],2) + pow(pos[1],2) + pow(pos[2],2));

        if (pos[2] / rr < cos_0 && rr > r - 2*eps){
            double Z1 = pos[2];
            double X1 = pos[0];
            double Y1 = pos[1];

            double Z2 = -r;
            double X2 = star1;
            double Y2 = star2;

            if (sqrt(pow(X1-X2,2)+pow(Y1-Y2,2)+pow(Z1-Z2,2)) < star3){
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
    int vertical = 720;
    int horizontal = (int) vertical*16/9;
    double z = (double)vertical*100/1080;
    double **screen = new double*[vertical];

    double cos_0 = 1.0 -2.0 *(pow(z,2)/(pow(z,2)+pow(horizontal/2,2)));
    double theta_0 = acos(cos_0);

    int num_stars = 1000;
    for (i=0;i<vertical;i++){
        screen[i] = new double[horizontal];
    }

    double **stars = new double*[num_stars + 1];
    for (i =0;i<num_stars + 1;i++){
        stars[i] = new double[3];
    }
    
    stars[0][0] = 0;
    stars[0][1] = 0;
    stars[0][2] = 30.0;

    for (i = 1;i<num_stars + 1;i++){
        stars[i][0] = -200+400*random_double();
        stars[i][1] = -200 + 400*random_double();
        stars[i][2] = 30*random_double();
    }

    for (i=0;i<vertical;i++){
        for (j=0;j<horizontal;j++){
            screen[i][j] = 0;
        }
    }

    for (i = 0;i<vertical;i++){
        for (j=0;j<horizontal;j++){
            ray Ray = ray((double)(j-horizontal/2),(double)(i-vertical/2),z);
            for (int k = 0; k<num_stars + 1;k++){
                bool check = Ray.check_hit(stars[k][0],stars[k][1], stars[k][2], cos_0);
                if (check){
                    screen[i][j] += 0.3;
                }
            }
        }
    }

    fstream fs;
    fs.open("with_gravity.csv", ios::out);

    for (i=0;i<vertical;i++){
        for (j = 0;j<horizontal;j++){
            double u = clip(screen[i][j],0,1);
            fs << u << ',';
        }
        fs << endl;
    }

    fs.close();

    for (i=0;i<vertical;i++){
        for (j=0;j<horizontal;j++){
            screen[i][j] = 0;
        }
    }

    for (i = 0;i<vertical;i++){
        for (j=0;j<horizontal;j++){
            ray1 Ray = ray1((double)(j-horizontal/2),(double)(i-vertical/2),z);
            for (int k = 0; k<num_stars + 1;k++){
                bool check = Ray.check_hit(stars[k][0],stars[k][1], stars[k][2], cos_0);
                if (check){
                    screen[i][j] += 0.3;
                }
            }
        }
    }

    fstream fs1;
    fs1.open("without_gravity.csv", ios::out);

    for (i=0;i<vertical;i++){
        for (j=0;j<horizontal;j++){
            double u = clip(screen[i][j],0,1);
            fs1 << u << ',';
        }
        fs1 << endl;
    }

    for (i=0; i<vertical;i++){
        delete screen[i];
    }
    delete screen;

    for (i = 0;i<num_stars + 1;i++){
        delete stars[i];
    }
    delete stars;
    return 0;
}