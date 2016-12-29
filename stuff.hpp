#include<armadillo>
#include<iostream>
#include<GL/glew.h>
#include<GL/glut.h>

using namespace std;
using namespace arma;

#define G 6.67e-11 //6.67e-11
#define EPSILON 1.0e-3
struct qtree_node
{
    vec q_max,q_min;
    vec q_center;
    double q_h;
    vec c_mass;
    double mass;
    int total;
    qtree_node * next[4];
    
    int get_quadrant(vec pos)
    {
        if(pos(0) < q_center(0))
        {
            if(pos(1) < q_center(1))
                return 0;
            else return 2;
        }
        else
        {
            if(pos(1) < q_center(1))
                return 1;
            else return 3;
        }
    }
    
    void set_quadrant(int q)
    {
        next[q] = new qtree_node;
        next[q]->q_min = q_min;
        next[q]->q_max = q_max;
        switch(q)
        {
            case 0:
                next[q]->q_max = q_center;
            break;
            case 1:
                next[q]->q_min(0) = q_center(0);
                next[q]->q_max(1) = q_center(1);
            break;
            case 2:
                next[q]->q_min(1) = q_center(1);
                next[q]->q_max(0) = q_center(0);
            break;
            case 3:
                next[q]->q_min = q_center;
            break;
        }
        next[q]->q_center = 0.5*(next[q]->q_min + next[q]->q_max);
        next[q]->q_h = next[q]->q_max(0) - next[q]->q_min(0);
    }
    
    void insert_to_node(double m,vec x)
    {
        if(total == 0)
        {
            c_mass = x;
            mass = m;
        }
        else if(total == 1)
        {
            int q_idx = get_quadrant(c_mass);
            if(next[q_idx]==NULL)
                set_quadrant(q_idx);
            next[q_idx]->insert_to_node(mass,c_mass);
            
            q_idx = get_quadrant(x);
            if(next[q_idx]==NULL)
                set_quadrant(q_idx);
            next[q_idx]->insert_to_node(m,x);
        }
        else
        {
            int q_idx = get_quadrant(x);
            if(next[q_idx]==NULL)
                set_quadrant(q_idx);
            next[q_idx]->insert_to_node(m,x);
        }
        total++;
    }
    
    void compute_mass()
    {
        if(total>1)
        {
            c_mass = zeros<vec>(2);
            mass = 0.0;
            for(int i=0;i<4;i++)
                if(next[i])
                {
                    next[i]->compute_mass();
                    mass += next[i]->mass;
                    c_mass += (next[i]->mass)*next[i]->c_mass;
                }
            c_mass /= mass;
        }
    }
    vec compute_force(double tet,vec x)
    {
        if(total==1)
        {
            double aux = pow(norm(x-c_mass),3.0);
            if(aux>EPSILON)
                return mass*G*(c_mass-x)/aux;
            else return zeros<vec>(2);
        }
        else
        {
            double r = norm(x-c_mass);
            if(q_h/r < tet)
            {
                double aux = pow(r,3.0);
                if(aux>EPSILON)
                    return mass*G*(c_mass-x)/aux;
                else return zeros<vec>(2);
            }
            else
            {
                vec aux = zeros<vec>(2);
                for(int i=0;i<4;i++)
                    if(next[i])
                        aux += next[i]->compute_force(tet,x);
                return aux;
            }
        }
    }
    
    qtree_node()
    {
        total = 0;
        c_mass = zeros<vec>(2);
        q_center = zeros<vec>(2);
        q_max = zeros<vec>(2);
        q_min = zeros<vec>(2);
        mass = 0.0;
        for(int i=0;i<4;i++) next[i] = NULL;
    }
    
    ~qtree_node()
    {
        for(int i=0;i<4;i++) if(next[i]) delete next[i];
    }
};

struct qtree
{
    qtree_node * root;
    qtree()
    {
        root = new qtree_node;
    }
    
    void insert(vec m,mat x)
    {
        reset();
        root->q_min = x.min()*ones<vec>(2);
        root->q_max = x.max()*ones<vec>(2);
        root->q_center = 0.5*(root->q_min + root->q_max);
        root->q_h = root->q_max(0) - root->q_min(0);
        for(int i=0;i<x.n_cols;i++)
            root->insert_to_node(m(i),x.col(i));
        root->compute_mass();
    }
    
    mat compute_force(double tet,mat x)
    {
        mat F(2,x.n_cols);
        for(int i=0;i<x.n_cols;i++)
            F.col(i) = root->compute_force(tet,x.col(i));
        return F;
    }
    
    void reset()
    {
        if(root!=NULL)delete root;
        root = new qtree_node;
    }
    ~qtree()
    {
        if(root!=NULL)delete root;
    }
};
