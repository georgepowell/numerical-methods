
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <vector>

using namespace std;

#define DT 0.1


double a(double input[]) {
    return -input[1];
}

double b(double input[]) {
    return input[0];
}

double c(double input[]) {
    return 1.0 	- input[0];
}

double d(double input[]) {
    return 3.0 / (2.0 * input[1] * input[1]) + input[0] / (2.0 * input[1]);
}

double time(double input[]) {
    return 1.0;
}

void forward(double values[], double (*derivatives[])(double[]), int size, double dt) {
    
    double rtns[size];
    
    for (int i = 0; i < size; i++) {
        rtns[i] = values[i] + derivatives[i](values) * dt; // forward euler
    }
    
    for (int j = 0; j < size; j++) {
        values[j] = rtns[j];
    }
}

void modified(double values[], double (*derivatives[])(double[]), int size, double dt) {
    
    double forwards[size];
    double mods[size];
    double ds[size];
    
    for (int i = 0; i < size; i++) {
        double di = derivatives[i](values);
        ds[i] = di;
        forwards[i] = values[i] + di * dt; // forward euler
    }
    
    for (int i = 0; i < size; i++) {
        mods[i] = (forwards[i] + values[i] + derivatives[i](forwards) * dt) / 2.0; // modified euler
    }
    
    for (int j = 0; j < size; j++) {
        values[j] = mods[j];
    }
}

void rk4(double values[], double (*derivatives[])(double[]), int size, double dt) {
    
    double d1[size];
    double vd1[size];
    
    double d2[size];
    double vd2[size];
    
    double d3[size];
    double vd3[size];
    
    double d4[size];
    
    for (int i = 0; i < size; i++) {
        d1[i] = dt * derivatives[i](values);
        vd1[i] = values[i] + 0.5 * d1[i];
    }
    
    for (int i = 0; i < size; i++) {
        d2[i] = dt * derivatives[i](vd1);
        vd2[i] = values[i] + 0.5 * d2[i];
    }
    
    for (int i = 0; i < size; i++) {
        d3[i] = dt * derivatives[i](vd2);
        vd3[i] = values[i] + d3[i];
    }
    
    for (int i = 0; i < size; i++) {
        d4[i] = dt * derivatives[i](vd3);
    }
    
    for (int i = 0; i < size; i++) {
        values[i] = values[i] + (d1[i] + 2.0 * d2[i] + 2.0 * d3[i] + d4[i]) / 6.0;
    }
}


double dt = DT;

void rk_dynamic(double values[], double (*derivatives[])(double[]), int size) {
    
    double values_temp[size];
    double values_temp2[size];
    
    double d1[size];
    double d2[size];
    double d3[size];
    double d4[size];
    double d5[size];
    double d6[size];
    double d7[size];
    
    double order4_solution[size];
    double order5_solution[size];
    
    double max_factor = 5000.0;
    
    double relative_error = 0.000000006;
    double abs_error = 0.000000006;
    
    for (int i = 0; i < size; i++) {
        d1[i] = dt * derivatives[i](values);
        values_temp[i] = values[i] + (1.0/5.0) * d1[i];
    }
    
    for (int i = 0; i < size; i++) {
        d2[i] = dt * derivatives[i](values_temp);
        values_temp2[i] = values[i] + (3.0/40.0) * d1[i] + (9.0/40.0) * d2[i];
    }
    
    for (int i = 0; i < size; i++) {
        d3[i] = dt * derivatives[i](values_temp2);
        values_temp[i] = values[i] + (44.0/45.0)      * d1[i] - (56.0/15.0)      * d2[i] + (32.0/9.0)       * d3[i];
    }
    
    for (int i = 0; i < size; i++) {
        d4[i] = dt * derivatives[i](values_temp);
        values_temp2[i] = values[i] + (19372.0/6561.0) * d1[i] - (25360.0/2187.0) * d2[i] + (64448.0/6561.0) * d3[i] - (212.0/729.0) * d4[i];
    }
    
    for (int i = 0; i < size; i++) {
        d5[i] = dt * derivatives[i](values_temp2);
        values_temp[i] = values[i] + (9017.0/3168.0)  * d1[i] - (355.0/33.0)     * d2[i] + (46732.0/5247.0) * d3[i] + (49.0/176.0)  * d4[i] - (5103.0/18656.0) * d5[i];
    }
    
    
    for (int i = 0; i < size; i++) {
        d6[i] = dt * derivatives[i](values_temp);
        
        order4_solution[i] =
            values_temp2[i] =
                values[i] + (35.0/384.0) * d1[i] + (500.0/1113.0) * d3[i] + (125.0/192.0) * d4[i] - (2187.0/6784.0) * d5[i] + (11.0/84.0) * d6[i];
    }
    
    for (int i = 0; i < size; i++) {
        d7[i] = dt * derivatives[i](values_temp2);
        
        order5_solution[i] =
            values[i] + (5179.0/57600.0) * d1[i] + (7571.0/16695.0) * d3[i] + (393.0/640.0) * d4[i] - (92097.0/339200.0) * d5[i] + (187.0/2100.0) * d6[i] + (1.0/40.0) * d7[i];
        
        double overall_tollerance = abs_error + abs(order5_solution[i] * relative_error);
        double local_calculated_factor = pow(abs(1.0 / (abs(order4_solution[i] - order5_solution[i]) / overall_tollerance)), (1.0 / 5.0));
        
        if (local_calculated_factor < max_factor)
            max_factor = local_calculated_factor;
    }
    
    dt *= max_factor; // Dynamically change dt based on estimated error
    
    if (max_factor < 0.9 || max_factor > 1.1) // Rerun time-step with smaller dt only if dt was significantly too large.
        rk_dynamic(values, derivatives, size);
    else // Otherwise accept current time-step and return the most accurate order5 solution.
        for (int i = 0; i < size; i++)
            values[i] = order5_solution[i];
}

ofstream forward_file;

double t = 0.0;

void pa(double x[], int numElements) {
    forward_file << t << "," << dt << ",";
    for (int i = 0; i < numElements - 1; i++)
        forward_file << x[i] << ",";
    forward_file << x[numElements - 1] << "\n";
}


#define ODE_SIZE 2

int main(int argc, const char * argv[])
{
    forward_file.open("general_test.txt");
    
    double (*derivatives[ODE_SIZE])(double[]) = { a, b };
    double vals[ODE_SIZE] = { 1.0, 0.0 };
    
    //double (*derivatives[ODE_SIZE])(double[]) = { c };
    //double vals[ODE_SIZE] = { 10.0 };
    
    //double (*derivatives[ODE_SIZE])(double[]) = { d, time };
    //double vals[ODE_SIZE] = { 0.0, 1.0 };
    
    //forward_file << "t,dt,a,b\n";
    forward_file << "t,dt,c\n";
    
    pa(vals, ODE_SIZE);
    
    for (int i = 0; i < 1000; i++) {
        forward(vals, derivatives, ODE_SIZE, DT);
        //modified(vals, derivatives, ODE_SIZE, DT);
        //rk4(vals, derivatives, ODE_SIZE, DT);
        //rk_dynamic(vals, derivatives, ODE_SIZE);
        pa(vals, ODE_SIZE);
        t += dt;
    }
    
    forward_file.close();
    
    return 0;
}
