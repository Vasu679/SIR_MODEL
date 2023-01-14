#include <iostream>
#include <cstdlib>
#include <fstream>
#include "cpgplot.h"

using namespace std;

float S(float S, float I, float betha_1){
    return (-betha_1)*I*S;
}
float I(float S, float I, float betha_1, float gamma_1){
    return (betha_1)*I*S-(gamma_1)*I;
}
float R(float I, float gamma_1){
    return (gamma_1)*I;
}


float Ss(float S, float I, float beta, float N){
    return (-beta)*S*(I/N); 
}
float Ee(float S, float I, float E, float beta, float alpha, float N){
    return ((beta)*(I/N)*S) - (alpha)*E;
}
float Ii(float E, float I, float alpha, float gamma){
    return ((alpha)*E) - (gamma)*I;
}
float Rr(float I, float gamma){
    return (gamma) * I;
}

int main() {

    int num1; 
    cout << ("Enter 1 for basic SIR Model and 2 for a model with vaccination considered: ") << "\n"; 
    cin >> num1;

        float h = 0.001;
        float min_t = 0.0;
        float max_t = 50.0;

        int total_pts = (max_t-min_t)/h;

        float t[total_pts];

        t[0] = min_t;
        t[1]= min_t + h;

        float betha_1 = 0.002;
        float gamma_1 = 0.45;

    if (num1 == 1){

        float S1[total_pts];
        float I1[total_pts];
        float R1[total_pts];

        I1[0] = 1.0;
        S1[0] = 770.0;
        R1[0] = 0.0;


        S1[1] = S1[0] + h*S(S1[0], I1[0], betha_1);
        I1[1] = I1[0] + h*I(S1[0], I1[0], betha_1, gamma_1);
        R1[1] = R1[0] + h*R(I1[0], gamma_1);

        for(int i = 2; i < total_pts; i++){
            t[i] = t[i-1] + h;
            S1[i] = S1[i-2] + 2*h*S(S1[i-1], I1[i-1], betha_1);
            I1[i] = I1[i-2] + 2*h*I(S1[i-1], I1[i-1], betha_1, gamma_1);
            R1[i] = R1[i-2] + 2*h*R(I1[i-1], gamma_1);
        }

        for(int j = total_pts - 1; j >= 0; j-- ){
            t[j] =  t[j] + 1;
        }

        size_t n = sizeof(S1)/sizeof(S1[0]);

        // **** Ploting The Coordinates ****

        // Open a plot window
        if (!cpgopen("/XWINDOW")) return 1;

        // Set-up plot axes
        cpgenv(1.,17.,0.,800.,0,1);

        // Label axes
        cpglab("days", "number of ppl", "SIR Model");
        cpgline(n,t,S1);
        cpgsci(4);
        cpgline(n,t,I1);
        cpgsci(5);
        cpgline(n,t,R1);

        cpgclos();

    }if (num1 == 2){

        float alpha = 0.124;
        float gamma = 0.156;
        float beta = 0.384;

        float S2[total_pts];
        float I2[total_pts];
        float R2[total_pts];
        float E2[total_pts]; 

        I2[0] = 20.0;
        S2[0] = 100.0;
        R2[0] = 0.0;
        E2[0] = (2.0 * I2[0]);
        float N = I2[0] + S2[0] + R2[0] + E2[0]; 
        float Tc[total_pts];
        float Ac[total_pts];
        Tc[0] = E2[0] + I2[0] + R2[0];
        Ac[0] = E2[0] + I2[0];

        S2[1] = S2[0] + h*Ss(S2[0], I2[0],beta, N);
        E2[1] = E2[0] + h*Ee(S2[0], I2[0], E2[0],beta,alpha, N);
        I2[1] = I2[0] + h*Ii(E2[0], I2[0],alpha, gamma);
        R2[1] = R2[0] + h*Rr(I2[0], gamma);
        Tc[1] = E2[1] + I2[1] + R2[1];
        Ac[1] = E2[1] + I2[1];

        for(int i = 2; i < total_pts; i++){
            t[i] = t[i-1] + h;
            S2[i] = S2[i-2] + 2*h*Ss(S2[i-1], I2[i-1],beta, N);
            E2[i] = E2[i-2] + 2*h*Ee(S2[i-1], I2[i-1], E2[i-1], beta,alpha, N);
            I2[i] = I2[i-2] + 2*h*Ii(E2[i-1], I2[i-1],alpha, gamma);
            R2[i] = R2[i-2] + 2*h*Rr(I2[i-1], gamma);
            Tc[i] = E2[i-1] + I2[i-1] + R2[i-1];
            Ac[i] = E2[i-1] + I2[i-1];
        }
        float m = (sizeof(S2)/sizeof(*S2));

        for(int k = total_pts - 1; k >= 0; k-- ){
            t[k] =  t[k] + 1;
        }

        // **** Ploting The Coordinates ****
        
        // Open a plot window
        if (!cpgopen("/XWINDOW")) return 1;

        // Set-up plot axes
        cpgenv(1.,50.,0.,150.,0,1);

        // Label axes
        cpglab("days", "number of ppl (prevalnce) ", "SIR Model");
        cpgline(m,t,S2);
        cpgsci(3);
        cpgline(m,t,E2);
        cpgsci(4);
        cpgline(m,t,I2);
        cpgsci(5);
        cpgline(m,t,R2);

        cpgclos();
    }

    
}
