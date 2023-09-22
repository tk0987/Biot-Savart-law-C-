#include <cstdio>
#include <iostream>
#include <cmath>
#include <math.h>
#include <time.h>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include<thread>
#define M_PI 3.14159265358979323846
using namespace std;

double  coil_angle_rad, coil_radius_m, coil_length_m, current_A,
x_m, y_m, z_m,
Bx_T, 
Bx1_T, Bx2_T, Bx3_T, Bx4_T, Bx5_T, Bx6_T, Bx7_T, Bx8_T, 
By_T, 
By1_T, By2_T, By3_T, By4_T, By5_T, By6_T, By7_T, By8_T,
Bz_T, 
Bz1_T, Bz2_T, Bz3_T, Bz4_T, Bz5_T, Bz6_T, Bz7_T, Bz8_T, 
B_T, 
proton_velocity_z = 102362127.1, step_z_m=0.00001, step_angle_rad = 0.0001;
long double euclidean_distance_m(double x1, double x2, double y1, double y2, double z1, double z2)
{
	return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2));
}

long double total_velocity_squared(double vx, double vy, double vz)
{
	return (vx * vx + vy * vy + vz * vz);
}
//________________________________________________________________________________________________________________________________________________________________________________________________
//threads
//________________________________________________________________________________________________________________________________________________________________________________________________
void thread1() //rods 1
{
	long double z1_m = -coil_length_m / 2;
	long double xp11_m = coil_radius_m * cos(coil_angle_rad + M_PI / 2); 
	long double yp11_m = coil_radius_m * sin(coil_angle_rad + M_PI / 2);


	Bx1_T = 0.0; By1_T = 0.0; Bz1_T = 0.0;
	for (z1_m = -coil_length_m / 2; z1_m <= coil_length_m / 2; z1_m += step_z_m)
	{
		long double r1_m = euclidean_distance_m(x_m, xp11_m, y_m, yp11_m, z_m, z1_m);

		Bx1_T += pow(10, -7) * (-1) * current_A * (-(y_m - yp11_m) * step_z_m) / abs(pow(r1_m, 3));
		By1_T += pow(10, -7) * (-1) * current_A * (step_z_m * (x_m - xp11_m)) / abs(pow(r1_m, 3));
		Bz1_T += 0.0;
	}

	z1_m = (-coil_length_m / 2) + 0.01;
	xp11_m = coil_radius_m * cos(M_PI / 3 + M_PI / 2); 
	yp11_m = coil_radius_m * sin(M_PI / 3 + M_PI / 2);


	//Bx1 = 0.0; By1 = 0.0; Bz1 = 0.0;
	for (z1_m = (-coil_length_m / 2) + 0.01; z1_m <= (coil_length_m / 2) - 0.01; z1_m += step_z_m)
	{
		long double r1_m = euclidean_distance_m(x_m, xp11_m, y_m, yp11_m, z_m, z1_m);

		Bx1_T += pow(10, -7) * (-1) * current_A * (-(y_m - yp11_m) * step_z_m) / abs(pow(r1_m, 3));
		By1_T += pow(10, -7) * (-1) * current_A * (step_z_m * (x_m - xp11_m)) / abs(pow(r1_m, 3));
		Bz1_T += 0.0;
	}

	z1_m = (-coil_length_m / 2) + 0.02;
	xp11_m = coil_radius_m * cos(M_PI / 6 + M_PI / 2); 
	yp11_m = coil_radius_m * sin(M_PI / 6 + M_PI / 2);


	//Bx1 = 0.0; By1 = 0.0; Bz1 = 0.0;
	for (z1_m = (-coil_length_m / 2) + 0.02; z1_m <= (coil_length_m / 2) - 0.02; z1_m += step_z_m)
	{
		long double r1_m = euclidean_distance_m(x_m, xp11_m, y_m, yp11_m, z_m, z1_m);

		Bx1_T += pow(10, -7) * (-1) * current_A * (-(y_m - yp11_m) * step_z_m) / abs(pow(r1_m, 3));
		By1_T += pow(10, -7) * (-1) * current_A * (step_z_m * (x_m - xp11_m)) / abs(pow(r1_m, 3));
		Bz1_T += 0.0;
	}
}
void thread2()//rods 2
{

	long double z2_m = -coil_length_m / 2;
	long double xp12_m = coil_radius_m * cos((-1) * coil_angle_rad + M_PI / 2); 
	long double yp12_m = coil_radius_m * sin((-1) * coil_angle_rad + M_PI / 2);
	Bx2_T = 0.0; By2_T = 0.0; Bz2_T = 0.0;
	for (z2_m = -coil_length_m / 2; z2_m <= coil_length_m / 2; z2_m += step_z_m)
	{
		long double r2_m = euclidean_distance_m(x_m, xp12_m, y_m, yp12_m, z_m, z2_m);

		Bx2_T += pow(10, -7) * current_A * (-(y_m - yp12_m) * step_z_m) / abs(pow(r2_m, 3));
		By2_T += pow(10, -7) * current_A * (step_z_m * (x_m - xp12_m)) / abs(pow(r2_m, 3));
		Bz2_T += 0.0;

	}
	z2_m = (-coil_length_m / 2) + 0.01;
	xp12_m = coil_radius_m * cos((-1) * M_PI / 3 + M_PI / 2); 
	yp12_m = coil_radius_m * sin((-1) * M_PI / 3 + M_PI / 2);
	//Bx2 = 0.0; By2 = 0.0; Bz2 = 0.0;
	for (z2_m = (-coil_length_m / 2) + 0.01; z2_m <= (coil_length_m / 2) - 0.01; z2_m += step_z_m)
	{
		long double r2_m = euclidean_distance_m(x_m, xp12_m, y_m, yp12_m, z_m, z2_m);

		Bx2_T += pow(10, -7) * current_A * (-(y_m - yp12_m) * step_z_m) / abs(pow(r2_m, 3));
		By2_T += pow(10, -7) * current_A * (step_z_m * (x_m - xp12_m)) / abs(pow(r2_m, 3));
		Bz2_T += 0.0;

	}
	z2_m = (-coil_length_m / 2) + 0.02;
	xp12_m = coil_radius_m * cos((-1) * M_PI / 6 + M_PI / 2);
	yp12_m = coil_radius_m * sin((-1) * M_PI / 6 + M_PI / 2);
	//Bx2 = 0.0; By2 = 0.0; Bz2 = 0.0;
	for (z2_m = z2_m = (-coil_length_m / 2) + 0.02; z2_m <= (coil_length_m / 2) - 0.02; z2_m += step_z_m)
	{
		long double r2_m = euclidean_distance_m(x_m, xp12_m, y_m, yp12_m, z_m, z2_m);

		Bx2_T += pow(10, -7) * current_A * (-(y_m - yp12_m) * step_z_m) / abs(pow(r2_m, 3));
		By2_T += pow(10, -7) * current_A * (step_z_m * (x_m - xp12_m)) / abs(pow(r2_m, 3));
		Bz2_T += 0.0;

	}
}
void thread3()//rods 3
{
	long double z3_m = -coil_length_m / 2;
	long double xp21_m = coil_radius_m * cos(coil_angle_rad + 3 * M_PI / 2); 
	long double yp21_m = coil_radius_m * sin(coil_angle_rad + 3 * M_PI / 2);

	Bx3_T = 0.0; By3_T = 0.0; Bz3_T = 0.0;
	for (z3_m = -coil_length_m / 2; z3_m <= coil_length_m / 2; z3_m += step_z_m)
	{
		long double r3_m = euclidean_distance_m(x_m, xp21_m, y_m, yp21_m, z_m, z3_m);

		Bx3_T += pow(10, -7) * current_A * (-(y_m - yp21_m) * step_z_m) / abs(pow(r3_m, 3));
		By3_T += pow(10, -7) * current_A * (step_z_m * (x_m - xp21_m)) / abs(pow(r3_m, 3));
		Bz3_T += 0.0;

	}
	z3_m = (-coil_length_m / 2) + 0.01;
	xp21_m = coil_radius_m * cos(M_PI / 3 + 3 * M_PI / 2); 
	yp21_m = coil_radius_m * sin(M_PI / 3 + 3 * M_PI / 2);

	//Bx3 = 0.0; By3 = 0.0; Bz3 = 0.0;
	for (z3_m = (-coil_length_m / 2) + 0.01; z3_m <= (coil_length_m / 2) - 0.01; z3_m += step_z_m)
	{
		long double r3_m = euclidean_distance_m(x_m, xp21_m, y_m, yp21_m, z_m, z3_m);

		Bx3_T += pow(10, -7) * current_A * (-(y_m - yp21_m) * step_z_m) / abs(pow(r3_m, 3));
		By3_T += pow(10, -7) * current_A * (step_z_m * (x_m - xp21_m)) / abs(pow(r3_m, 3));
		Bz3_T += 0.0;

	}
	z3_m = (-coil_length_m / 2) + 0.02;
	xp21_m = coil_radius_m * cos(M_PI / 6 + 3 * M_PI / 2); 
	yp21_m = coil_radius_m * sin(M_PI / 6 + 3 * M_PI / 2);

	//Bx3 = 0.0; By3 = 0.0; Bz3 = 0.0;
	for (z3_m = (-coil_length_m / 2) + 0.02; z3_m <= (coil_length_m / 2) - 0.02; z3_m += step_z_m)
	{
		long double r3_m = euclidean_distance_m(x_m, xp21_m, y_m, yp21_m, z_m, z3_m);

		Bx3_T += pow(10, -7) * current_A * (-(y_m - yp21_m) * step_z_m) / abs(pow(r3_m, 3));
		By3_T += pow(10, -7) * current_A * (step_z_m * (x_m - xp21_m)) / abs(pow(r3_m, 3));
		Bz3_T += 0.0;

	}

}
void thread4()//rods 4
{
	long double z4_m = -coil_length_m / 2;


	long double xp22_m = coil_radius_m * cos((-1) * coil_angle_rad + 3 * M_PI / 2); 
	long double yp22_m = coil_radius_m * sin((-1) * coil_angle_rad + 3 * M_PI / 2);
	long double Bx4_T = 0.0; By4_T = 0.0; Bz4_T = 0.0;
	for (z4_m = -coil_length_m / 2; z4_m <= coil_length_m / 2; z4_m += step_z_m)
	{
		long double r4_m = euclidean_distance_m(x_m, xp22_m, y_m, yp22_m, z_m, z4_m);

		Bx4_T += pow(10, -7) * (-1) * current_A * (-(y_m - yp22_m) * step_z_m) / abs(pow(r4_m, 3));
		By4_T += pow(10, -7) * (-1) * current_A * (step_z_m * (x_m - xp22_m)) / abs(pow(r4_m, 3));
		Bz4_T += 0.0;

	}
	z4_m = (-coil_length_m / 2) + 0.01;


	xp22_m = coil_radius_m * cos((-1) * M_PI / 3 + 3 * M_PI / 2); 
	yp22_m = coil_radius_m * sin((-1) * M_PI / 3 + 3 * M_PI / 2);
	//Bx4 = 0.0; By4 = 0.0; Bz4 = 0.0;
	for (z4_m = (-coil_length_m / 2) + 0.01; z4_m <= (coil_length_m / 2) - 0.01; z4_m += step_z_m)
	{
		long double r4_m = euclidean_distance_m(x_m, xp22_m, y_m, yp22_m, z_m, z4_m);

		Bx4_T += pow(10, -7) * (-1) * current_A * (-(y_m - yp22_m) * step_z_m) / abs(pow(r4_m, 3));
		By4_T += pow(10, -7) * (-1) * current_A * (step_z_m * (x_m - xp22_m)) / abs(pow(r4_m, 3));
		Bz4_T += 0.0;

	}
	z4_m = (-coil_length_m / 2) + 0.02;


	xp22_m = coil_radius_m * cos((-1) * M_PI / 6 + 3 * M_PI / 2); 
	yp22_m = coil_radius_m * sin((-1) * M_PI / 6 + 3 * M_PI / 2);
	//Bx4 = 0.0; By4 = 0.0; Bz4 = 0.0;
	for (z4_m = (-coil_length_m / 2) + 0.02; z4_m <= (coil_length_m / 2) - 0.02; z4_m += step_z_m)
	{
		long double r4_m = euclidean_distance_m(x_m, xp22_m, y_m, yp22_m, z_m, z4_m);

		Bx4_T += pow(10, -7) * (-1) * current_A * (-(y_m - yp22_m) * step_z_m) / abs(pow(r4_m, 3));
		By4_T += pow(10, -7) * (-1) * current_A * (step_z_m * (x_m - xp22_m)) / abs(pow(r4_m, 3));
		Bz4_T += 0.0;

	}
}
void thread5()//arcs 1
{
	long double angle1_rad = (-1) * coil_angle_rad + (M_PI / 2);
	Bx5_T = 0.0; By5_T = 0.0; Bz5_T = 0.0;
	for (angle1_rad = (-1) * coil_angle_rad + (M_PI / 2); angle1_rad <= coil_angle_rad + (M_PI / 2); angle1_rad += step_angle_rad) 
	{
		long double yo1_m = coil_radius_m * sin(angle1_rad);
		long double xo1_m = coil_radius_m * cos(angle1_rad);
		long double zo1_m = (-coil_length_m / 2);
		long double r5_m = euclidean_distance_m(x_m, xo1_m, y_m, yo1_m, z_m, zo1_m);
		long double step_x_m = coil_radius_m * cos(angle1_rad) - coil_radius_m * cos(angle1_rad - step_angle_rad);
		long double step_y_m = coil_radius_m * sin(angle1_rad) - coil_radius_m * sin(angle1_rad - step_angle_rad);
		Bx5_T += pow(10, -7) * current_A * ((z_m - zo1_m) * (step_y_m)) / abs(pow(r5_m, 3));
		By5_T += pow(10, -7) * current_A * ((-1) * (z_m - zo1_m) * (step_x_m)) / abs(pow(r5_m, 3));
		Bz5_T += pow(10, -7) * current_A * ((-1) * (x_m - xo1_m) * step_y_m + (y_m - yo1_m) * (step_x_m)) / abs(pow(r5_m, 3));

	}
	angle1_rad = (-1) * M_PI / 3 + (M_PI / 2);
	//Bx5 = 0.0; By5 = 0.0; Bz5 = 0.0;
	for (angle1_rad = (-1) * M_PI / 3 + (M_PI / 2); angle1_rad <= M_PI / 3 + (M_PI / 2); angle1_rad += step_angle_rad) 
	{
		long double yo1_m = coil_radius_m * sin(angle1_rad);
		long double xo1_m = coil_radius_m * cos(angle1_rad);
		long double zo1_m = (-coil_length_m / 2) + 0.01;
		long double r5_m = euclidean_distance_m(x_m, xo1_m, y_m, yo1_m, z_m, zo1_m);
		long double step_x_m = coil_radius_m * cos(angle1_rad) - coil_radius_m * cos(angle1_rad - step_angle_rad);
		long double step_y_m = coil_radius_m * sin(angle1_rad) - coil_radius_m * sin(angle1_rad - step_angle_rad);
		Bx5_T += pow(10, -7) * current_A * ((z_m - zo1_m) * (step_y_m)) / abs(pow(r5_m, 3));
		By5_T += pow(10, -7) * current_A * ((-1) * (z_m - zo1_m) * (step_x_m)) / abs(pow(r5_m, 3));
		Bz5_T += pow(10, -7) * current_A * ((-1) * (x_m - xo1_m) * step_y_m + (y_m - yo1_m) * (step_x_m)) / abs(pow(r5_m, 3));

	}
	angle1_rad = (-1) * M_PI / 6 + (M_PI / 2);
	//Bx5 = 0.0; By5 = 0.0; Bz5 = 0.0;
	for (angle1_rad = (-1) * M_PI / 6 + (M_PI / 2); angle1_rad <= M_PI / 6 + (M_PI / 2); angle1_rad += step_angle_rad) 
	{
		long double yo1_m = coil_radius_m * sin(angle1_rad);
		long double xo1_m = coil_radius_m * cos(angle1_rad);
		long double zo1_m = (-coil_length_m / 2) + 0.02;
		long double r5_m = euclidean_distance_m(x_m, xo1_m, y_m, yo1_m, z_m, zo1_m);
		long double step_x_m = coil_radius_m * cos(angle1_rad) - coil_radius_m * cos(angle1_rad - step_angle_rad);
		long double step_y_m = coil_radius_m * sin(angle1_rad) - coil_radius_m * sin(angle1_rad - step_angle_rad);
		Bx5_T += pow(10, -7) * current_A * ((z_m - zo1_m) * (step_y_m)) / abs(pow(r5_m, 3));
		By5_T += pow(10, -7) * current_A * ((-1) * (z_m - zo1_m) * (step_x_m)) / abs(pow(r5_m, 3));
		Bz5_T += pow(10, -7) * current_A * ((-1) * (x_m - xo1_m) * step_y_m + (y_m - yo1_m) * (step_x_m)) / abs(pow(r5_m, 3));

	}


}
void thread6()//arcs2
{
	long double angle2_rad = coil_angle_rad + (M_PI / 2);
	Bx6_T = 0.0; By6_T = 0.0; Bz6_T = 0.0;
	for (angle2_rad = coil_angle_rad + (M_PI / 2); angle2_rad >= (-1) * coil_angle_rad + (M_PI / 2); angle2_rad -= step_angle_rad)
	{
		long double yo2_m = coil_radius_m * sin(angle2_rad);
		long double xo2_m = coil_radius_m * cos(angle2_rad);
		long double zo2_m = coil_length_m / 2;
		long double r6_m = euclidean_distance_m(x_m, xo2_m, y_m, yo2_m, z_m, zo2_m);
		long double step_x_m = coil_radius_m * cos(angle2_rad) - coil_radius_m * cos(angle2_rad + step_angle_rad);
		long double step_y_m = coil_radius_m * sin(angle2_rad) - coil_radius_m * sin(angle2_rad + step_angle_rad);
		Bx6_T += pow(10, -7) * (-1) * current_A * ((z_m - zo2_m) * (step_y_m)) / abs(pow(r6_m, 3));
		By6_T += pow(10, -7) * (-1) * current_A * ((-1) * (z_m - zo2_m) * (step_x_m)) / abs(pow(r6_m, 3));
		Bz6_T += pow(10, -7) * (-1) * current_A * ((-1) * (x_m - xo2_m) * (step_y_m)+(y_m - yo2_m) * (step_x_m)) / abs(pow(r6_m, 3));

	}
	angle2_rad = M_PI / 3 + (M_PI / 2);
	//Bx6 = 0.0; By6 = 0.0; Bz6 = 0.0;
	for (angle2_rad = M_PI / 3 + (M_PI / 2); angle2_rad >= (-1) * M_PI / 3 + (M_PI / 2); angle2_rad -= step_angle_rad)
	{
		long double yo2_m = coil_radius_m * sin(angle2_rad);
		long double xo2_m = coil_radius_m * cos(angle2_rad);
		long double zo2_m = (coil_length_m / 2) - 0.01;
		long double r6_m = euclidean_distance_m(x_m, xo2_m, y_m, yo2_m, z_m, zo2_m);
		long double step_x_m = coil_radius_m * cos(angle2_rad) - coil_radius_m * cos(angle2_rad + step_angle_rad);
		long double step_y_m = coil_radius_m * sin(angle2_rad) - coil_radius_m * sin(angle2_rad + step_angle_rad);
		Bx6_T += pow(10, -7) * (-1) * current_A * ((z_m - zo2_m) * (step_y_m)) / abs(pow(r6_m, 3));
		By6_T += pow(10, -7) * (-1) * current_A * ((-1) * (z_m - zo2_m) * (step_x_m)) / abs(pow(r6_m, 3));
		Bz6_T += pow(10, -7) * (-1) * current_A * ((-1) * (x_m - xo2_m) * (step_y_m)+(y_m - yo2_m) * (step_x_m)) / abs(pow(r6_m, 3));

	}
	angle2_rad = M_PI / 6 + (M_PI / 2);
	//Bx6 = 0.0; By6 = 0.0; Bz6 = 0.0;
	for (angle2_rad = M_PI / 6 + (M_PI / 2); angle2_rad >= (-1) * M_PI / 6 + (M_PI / 2); angle2_rad -= step_angle_rad)
	{
		long double yo2_m = coil_radius_m * sin(angle2_rad);
		long double xo2_m = coil_radius_m * cos(angle2_rad);
		long double zo2_m = (coil_length_m / 2) - 0.02;
		long double r6_m = euclidean_distance_m(x_m, xo2_m, y_m, yo2_m, z_m, zo2_m);
		long double step_x_m = coil_radius_m * cos(angle2_rad) - coil_radius_m * cos(angle2_rad + step_angle_rad);
		long double step_y_m = coil_radius_m * sin(angle2_rad) - coil_radius_m * sin(angle2_rad + step_angle_rad);
		Bx6_T += pow(10, -7) * (-1) * current_A * ((z_m - zo2_m) * (step_y_m)) / abs(pow(r6_m, 3));
		By6_T += pow(10, -7) * (-1) * current_A * ((-1) * (z_m - zo2_m) * (step_x_m)) / abs(pow(r6_m, 3));
		Bz6_T += pow(10, -7) * (-1) * current_A * ((-1) * (x_m - xo2_m) * (step_y_m)+(y_m - yo2_m) * (step_x_m)) / abs(pow(r6_m, 3));

	}
}
void thread7()//arcs 3 
{
	long double angle3_rad = coil_angle_rad + (3 * M_PI / 2);
	Bx7_T = 0.0; By7_T = 0.0; Bz7_T = 0.0;
	for (angle3_rad = coil_angle_rad + (3 * M_PI / 2); angle3_rad >= (-1) * coil_angle_rad + (3 * M_PI / 2); angle3_rad -= step_angle_rad) 
	{
		long double yo3_m = coil_radius_m * sin(angle3_rad);
		long double xo3_m = coil_radius_m * cos(angle3_rad);
		long double zo3_m = -coil_length_m / 2;
		long double r7_m = euclidean_distance_m(x_m, xo3_m, y_m, yo3_m, z_m, zo3_m);
		long double step_x_m = coil_radius_m * cos(angle3_rad) - coil_radius_m * cos(angle3_rad + step_angle_rad);
		long double step_y_m = coil_radius_m * sin(angle3_rad) - coil_radius_m * sin(angle3_rad + step_angle_rad);
		Bx7_T += pow(10, -7) * current_A * ((z_m - zo3_m) * (step_y_m)) / abs(pow(r7_m, 3));
		By7_T += pow(10, -7) * current_A * ((-1) * (z_m - zo3_m) * (step_x_m)) / abs(pow(r7_m, 3));
		Bz7_T += pow(10, -7) * current_A * ((-1) * (x_m - xo3_m) * (step_y_m)+(y_m - yo3_m) * (step_x_m)) / abs(pow(r7_m, 3));

	}
	angle3_rad = M_PI / 3 + (3 * M_PI / 2);
	//Bx7 = 0.0; By7 = 0.0; Bz7 = 0.0;
	for (angle3_rad = M_PI / 3 + (3 * M_PI / 2); angle3_rad >= (-1) * M_PI / 3 + (3 * M_PI / 2); angle3_rad -= step_angle_rad) 
	{
		long double yo3_m = coil_radius_m * sin(angle3_rad);
		long double xo3_m = coil_radius_m * cos(angle3_rad);
		long double zo3_m = (-coil_length_m / 2) + 0.01;
		long double r7_m = euclidean_distance_m(x_m, xo3_m, y_m, yo3_m, z_m, zo3_m);
		long double step_x_m = coil_radius_m * cos(angle3_rad) - coil_radius_m * cos(angle3_rad + step_angle_rad);
		long double step_y_m = coil_radius_m * sin(angle3_rad) - coil_radius_m * sin(angle3_rad + step_angle_rad);
		Bx7_T += pow(10, -7) * current_A * ((z_m - zo3_m) * (step_y_m)) / abs(pow(r7_m, 3));
		By7_T += pow(10, -7) * current_A * ((-1) * (z_m - zo3_m) * (step_x_m)) / abs(pow(r7_m, 3));
		Bz7_T += pow(10, -7) * current_A * ((-1) * (x_m - xo3_m) * (step_y_m)+(y_m - yo3_m) * (step_x_m)) / abs(pow(r7_m, 3));

	}
	angle3_rad = M_PI / 6 + (3 * M_PI / 2);
	//Bx7 = 0.0; By7 = 0.0; Bz7 = 0.0;
	for (angle3_rad = M_PI / 6 + (3 * M_PI / 2); angle3_rad >= (-1) * M_PI / 6 + (3 * M_PI / 2); angle3_rad -= step_angle_rad) 
	{
		long double yo3_m = coil_radius_m * sin(angle3_rad);
		long double xo3_m = coil_radius_m * cos(angle3_rad);
		long double zo3_m = (-coil_length_m / 2) + 0.02;
		long double r7_m = euclidean_distance_m(x_m, xo3_m, y_m, yo3_m, z_m, zo3_m);
		long double step_x_m = coil_radius_m * cos(angle3_rad) - coil_radius_m * cos(angle3_rad + step_angle_rad);
		long double step_y_m = coil_radius_m * sin(angle3_rad) - coil_radius_m * sin(angle3_rad + step_angle_rad);
		Bx7_T += pow(10, -7) * current_A * ((z_m - zo3_m) * (step_y_m)) / abs(pow(r7_m, 3));
		By7_T += pow(10, -7) * current_A * ((-1) * (z_m - zo3_m) * (step_x_m)) / abs(pow(r7_m, 3));
		Bz7_T += pow(10, -7) * current_A * ((-1) * (x_m - xo3_m) * (step_y_m)+(y_m - yo3_m) * (step_x_m)) / abs(pow(r7_m, 3));

	}

}
void thread8()//arcs 4 
{
	long double angle4_rad = (-1) * coil_angle_rad + (3 * M_PI / 2);
	Bx8_T = 0.0; By8_T = 0.0; Bz8_T = 0.0;
	for (angle4_rad = (-1) * coil_angle_rad + (3 * M_PI / 2); angle4_rad <= coil_angle_rad + (3 * M_PI / 2); angle4_rad += step_angle_rad) 
	{
		long double yo4_m = coil_radius_m * sin(angle4_rad);
		long double xo4_m = coil_radius_m * cos(angle4_rad);
		long double zo4_m = coil_length_m / 2;
		long double r8_m = euclidean_distance_m(x_m, xo4_m, y_m, yo4_m, z_m, zo4_m);
		long double step_x_m = coil_radius_m * cos(angle4_rad) - coil_radius_m * cos(angle4_rad - step_angle_rad);
		long double step_y_m = coil_radius_m * sin(angle4_rad) - coil_radius_m * sin(angle4_rad - step_angle_rad);
		Bx8_T += pow(10, -7) * (-1) * current_A * ((z_m - zo4_m) * (step_y_m)) / abs(pow(r8_m, 3));
		By8_T += pow(10, -7) * (-1) * current_A * ((-1) * (z_m - zo4_m) * (step_x_m)) / abs(pow(r8_m, 3));
		Bz8_T += pow(10, -7) * (-1) * current_A * ((-1) * (x_m - xo4_m) * (step_y_m)+(y_m - yo4_m) * (step_x_m)) / abs(pow(r8_m, 3));

	}
	angle4_rad = (-1) * M_PI / 3 + (3 * M_PI / 2);
	//Bx8 = 0.0; By8 = 0.0; Bz8 = 0.0;
	for (angle4_rad = (-1) * M_PI / 3 + (3 * M_PI / 2); angle4_rad <= M_PI / 3 + (3 * M_PI / 2); angle4_rad += step_angle_rad) 
	{
		long double yo4_m = coil_radius_m * sin(angle4_rad);
		long double xo4_m = coil_radius_m * cos(angle4_rad);
		long double zo4_m = (coil_length_m / 2) - 0.01;
		long double r8_m = euclidean_distance_m(x_m, xo4_m, y_m, yo4_m, z_m, zo4_m);
		long double step_x_m = coil_radius_m * cos(angle4_rad) - coil_radius_m * cos(angle4_rad - step_angle_rad);
		long double step_y_m = coil_radius_m * sin(angle4_rad) - coil_radius_m * sin(angle4_rad - step_angle_rad);
		Bx8_T += pow(10, -7) * (-1) * current_A * ((z_m - zo4_m) * (step_y_m)) / abs(pow(r8_m, 3));
		By8_T += pow(10, -7) * (-1) * current_A * ((-1) * (z_m - zo4_m) * (step_x_m)) / abs(pow(r8_m, 3));
		Bz8_T += pow(10, -7) * (-1) * current_A * ((-1) * (x_m - xo4_m) * (step_y_m)+(y_m - yo4_m) * (step_x_m)) / abs(pow(r8_m, 3));

	}
	angle4_rad = (-1) * M_PI / 6 + (3 * M_PI / 2);
	//Bx8 = 0.0; By8 = 0.0; Bz8 = 0.0;
	for (angle4_rad = (-1) * M_PI / 6 + (3 * M_PI / 2); angle4_rad <= M_PI / 6 + (3 * M_PI / 2); angle4_rad += step_angle_rad) 
	{
		long double yo4_m = coil_radius_m * sin(angle4_rad);
		long double xo4_m = coil_radius_m * cos(angle4_rad);
		long double zo4_m = (coil_length_m / 2) - 0.02;
		long double r8_m = euclidean_distance_m(x_m, xo4_m, y_m, yo4_m, z_m, zo4_m);
		long double step_x_m = coil_radius_m * cos(angle4_rad) - coil_radius_m * cos(angle4_rad - step_angle_rad);
		long double step_y_m = coil_radius_m * sin(angle4_rad) - coil_radius_m * sin(angle4_rad - step_angle_rad);
		Bx8_T += pow(10, -7) * (-1) * current_A * ((z_m - zo4_m) * (step_y_m)) / abs(pow(r8_m, 3));
		By8_T += pow(10, -7) * (-1) * current_A * ((-1) * (z_m - zo4_m) * (step_x_m)) / abs(pow(r8_m, 3));
		Bz8_T += pow(10, -7) * (-1) * current_A * ((-1) * (x_m - xo4_m) * (step_y_m)+(y_m - yo4_m) * (step_x_m)) / abs(pow(r8_m, 3));

	}
}
int main()
{
	ofstream file; file.open("dane_prad_wspolbiezny_dane.ok150e-16_89deg_60deg_30deg_3000A.txt");
	x_m = 0.0; y_m = 0.0; z_m = -1.5 * coil_length_m / 2;
	coil_angle_rad = 2 * M_PI * 89 / 360; //w radianach
	coil_radius_m = 5.5 / 100.0; //w metrach
	coil_length_m = 0.40; //w metrach
	current_A = 3000.0; //amperów
	x_m = 0.0; y_m = 0.0; z_m = -1.1 * coil_length_m ;
	Bx_T = 0.0; By_T = 0.0; Bz_T = 0.0;
	long double time_step_s = 150.0*1E-16;
	long double vx = 0.0, vy = 0.0, vz = proton_velocity_z;
	long double mass_kg = 1.78E-27;//masa relatywistyczna w kg
	long double q = 1.6021766208E-19;//ladunek protonu w C
	int i = 0;
	while (z_m <= 1.1*coil_length_m/2 )
	{
		i++;
		thread t1(thread1);
		thread t2(thread2);
		thread t3(thread3);
		thread t4(thread4);
		thread t5(thread5);
		thread t6(thread6);
		thread t7(thread7);
		thread t8(thread8);
		t1.join(); t2.join(); t3.join(); t4.join(); 
		t5.join(); t6.join(); t7.join(); t8.join();
		Bx_T = Bx1_T + Bx2_T + Bx3_T + Bx4_T + Bx5_T + Bx6_T + Bx7_T + Bx8_T;
		By_T = By1_T + By2_T + By3_T + By4_T + By5_T + By6_T + By7_T + By8_T;
		Bz_T = Bz1_T + Bz2_T + Bz3_T + Bz4_T + Bz5_T + Bz6_T + Bz7_T + Bz8_T;
		//B = sqrt(Bx * Bx + By * By + Bz * Bz);
		long double Fx = q * (vy * Bz_T - vz * By_T);
		long double Fy = q * (vz * Bx_T - vx * Bz_T);
		long double Fz = q * (vx * By_T - vy * Bx_T);

		long double ax = Fx / mass_kg;
		long double ay = Fy / mass_kg;
		long double az = Fz / mass_kg;

		vx += ax * time_step_s;
		vy += ay * time_step_s;
		vz = (sqrt(pow(102362127.1, 2) - pow(vy, 2) - pow(vx, 2)) - az * time_step_s);

		x_m += vx * time_step_s;
		y_m += vy * time_step_s;
		z_m += vz * time_step_s;
		
		file
			<< setprecision(10) << x_m << "	"
			<< setprecision(10) << y_m << "	"
			<< setprecision(10) << z_m << "	"
			<< setprecision(10) << vx << "	"
			<< setprecision(10) << vy << "	"
			<< setprecision(25) << vz << endl;
		if (i % 150 == 0) 
		{cout << abs((z_m + 1.1 *coil_length_m) / (2.2 * coil_length_m)) << "	"
			<< setprecision(10) << x_m << "	"
			<< setprecision(10) << y_m << "	"
			<< setprecision(10) << z_m << "	" 
			<< setprecision(10) << Bx_T << "	"
			<< setprecision(10) << By_T << "	"
			<< setprecision(10) << Bz_T << endl;
		}
