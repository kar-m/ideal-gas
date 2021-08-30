#include <iostream>
#include <cmath>
#include <random>
#include <time.h>
#include <fstream>
#include <stdio.h>
#include <cstdio>
using namespace std;


double random_number()
{
    return ((double) rand() / (RAND_MAX));
}

double rand_gen() {
   // return a uniformly distributed random value
   return ( (double)(rand()) + 1. )/( (double)(RAND_MAX) + 1. );
}

double normalRandom() {
   // return a normally distributed random value
   double v1=rand_gen();
   double v2=rand_gen();
   return cos(2*3.14*v2)*sqrt(-2.*log(v1));
}

class c_vector{
public:
    
    double x;
    double y;
    double z;
    
    c_vector()
    {
        this->x = 0;
        this->y = 0;
        this->z = 0;
    }
    
    c_vector(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }
    
    static c_vector add(c_vector a, c_vector b)
    {
        a.x += b.x;
        a.y += b.y;
        a.z += b.z;
        return a;
    }
    
    static c_vector subtract(c_vector a, c_vector b)
    {
        a.x -= b.x;
        a.y -= b.y;
        a.z -= b.z;
        return a;
    }
    
    static c_vector multiply(c_vector a, double b)
    {
        a.x *= b;
        a.y *= b;
        a.z *= b;
        return a;
    }
    
    static c_vector divide(c_vector a, double b)
    {
        try{
            a.x /= b;
            a.y /= b;
            a.z /= b;
            return a;
        }
        catch(int e)
        {
            cout << "No Division by 0 Allowed";
            return a;
        }
    }
    
    double length()
    {
        return sqrt(x*x + y*y + z*z);
    }
    
    static c_vector random_uniform(double x_max, double y_max, double z_max, bool min_zero = true)
    {
        
        if (min_zero){
            double x = random_number() * x_max;
            double y = random_number() * y_max;
            double z = random_number() * z_max;
            
            c_vector c;
            c.x = x;
            c.y = y;
            c.z = z;
            return c;
        }
        
        else{
            double x = random_number() * (2*x_max) - x_max;
            double y = random_number() * (2*y_max) - y_max;
            double z = random_number() * (2*z_max) - z_max;
            
            c_vector c;
            c.x = x;
            c.y = y;
            c.z = z;
            return c;
        }
    }
    
    static c_vector random_gaussian(double myu_x, double myu_y, double myu_z, double sigma_x, double sigma_y, double sigma_z)
    {
        
        double x = sigma_x * normalRandom() + myu_x;
        double y = sigma_y * normalRandom() + myu_y;
        double z = sigma_z * normalRandom() + myu_z;
        
        c_vector c;
        c.x = x;
        c.y = y;
        c.z = z;
        return c;
    }
};

char store[] = "/Users/macuser/Desktop/Xcode Projects/Simulations/IdealGas/IdealGas/coordinates.csv";
c_vector null_vector = *new c_vector(0, 0, 0);
double tick_length = 0.01;
int tick_number = 10000;
const int particle_number = 2000;
double box_side = 250;
double particle_radius = 3;
double T = 273.15;
double k_b = 1380.65;
double repel_coefficient = 5000;

class particle{
public:
    c_vector pos = null_vector;
    c_vector vel = null_vector;
    c_vector acc = null_vector;
    int id = 0;
    double rho = 8;
    double r;
    double m;
    double sigma;
    particle()
    {
        r = 0;
        m = 0;
        sigma = 0;
    }
    particle(double radius){
        r = radius;
        m = rho*r*r*r;
    }
    
    void tick(){
        pos = c_vector::add(pos, c_vector::multiply(vel, tick_length));
        vel = c_vector::add(vel, c_vector::multiply(acc, tick_length));
    }
    
    particle collision_detection(particle* particles)
    {
        for(int i = 0; i < particle_number; i++)
        {
            if(particles[i].id != this->id)
            {
                if (c_vector::subtract(particles[i].pos, this->pos).length() < particles[i].r + this->r)
                    return particles[i];
            }
        }
        return *this;
    }
    
    c_vector wall_collision()
    {
        if(pos.x < r){
            c_vector c;
            c.x = 1;
            c.y = 0;
            c.z = 0;
            return c;
        }
        else if (pos.x > box_side - r){
            c_vector c;
            c.x = -1;
            c.y = 0;
            c.z = 0;
            return c;
        }
        else if(pos.y < r){
            c_vector c;
            c.x = 0;
            c.y = 1;
            c.z = 0;
            return c;
        }
        else if (pos.y > box_side - r){
            c_vector c;
            c.x = 0;
            c.y = -1;
            c.z = 0;
            return c;
        }
        else if(pos.z < r){
            c_vector c;
            c.x = 0;
            c.y = 0;
            c.z = 1;
            return c;
        }
        else if (pos.z > box_side - r){
            c_vector c;
            c.x = 0;
            c.y = 0;
            c.z = -1;
            return c;
        }
        else
            return null_vector;
    }
};

particle particles[particle_number];

particle* initialize()
{
    for(int i = 0; i < particle_number; i++)
    {
        particle* a = &(particles[i]);
        a->id = i;
        a->r = particle_radius;
        a->m = a->rho * (a->r) * (a->r) * (a->r);
        a->sigma = sqrt((k_b * T)/a->m);
        a->pos = c_vector::random_uniform(box_side, box_side, box_side);
        a->vel = c_vector::random_gaussian(0, 0, 0, a->sigma, a->sigma, a->sigma);
    }

    cout <<endl;
    return particles;
}

void record_one_coordinate(){

    
    ofstream MyFile( store, std::ios::app ) ;
    
    for(int i = 0; i < particle_number; i++)
    {
        MyFile << particles[i].pos.x << ',';
    }
    
    MyFile << '\n';
    
    MyFile.close();
}

int main(int argc, const char * argv[]) {
    srand(time(NULL));
    int result = remove(store);
    cout << result << endl;
    initialize();
    for(int t = 0; t < tick_number; t++){
        record_one_coordinate();
        for(int i = 0; i < particle_number; i++)
        {
            auto_ptr<c_vector> F(new c_vector);
            *F = null_vector;
            auto_ptr<particle> collision(new particle);
            *collision = particles[i].collision_detection(particles);
            if(collision->id != particles[i].id)
            {
                auto_ptr<c_vector> diff(new c_vector);
                *diff = c_vector::subtract(particles[i].pos, collision->pos);
                *F = c_vector::multiply(*diff, repel_coefficient);
                *F = c_vector::divide(*F, diff->length()+0.01);
            }
            *F = c_vector::add(*F, c_vector::multiply(particles[i].wall_collision(), 10*repel_coefficient));
            particles[i].acc = c_vector::divide(*F, particles[i].m);
        }
        for(int i = 0; i < particle_number; i++)
        {
            particles[i].tick();
            particles[i].acc = null_vector;
        }
        cout << "\n";
        cout << t << "\n";
        
    }

    cout << "HI!";
    return 0;
}
