#include "stuff.hpp"

const int N = 200;
const double dt = 1e-2;
const double MAX = 7.5e2;
const double MAX_DRAW = 1e3;
const double THETA = 0.3;
qtree great_tree;
vec M(N);
mat X(2,N);
mat V(2,N);

void F(vec M,mat Xs,mat Vs,mat & Xe,mat & Ve)
{
  Xe = Vs;
  great_tree.insert(M,Xs);
  Ve = great_tree.compute_force(THETA,Xs);
}
void rk4(vec M,mat & X0,mat & V0,double dt)
{
  mat X1,X2,X3,X4,V1,V2,V3,V4;
  F(M,X0,V0,X1,V1);
  F(M,X0+0.5*dt*X1,V0+0.5*dt*V1,X2,V2);
  F(M,X0+0.5*dt*X2,V0+0.5*dt*V2,X3,V3);
  F(M,X0+dt*X3,V0+dt*V3,X4,V4);
  X0 = X0 + dt*(X1+2.0*X2+2.0*X3+X4)/6.0;
  V0 = V0 + dt*(V1+2.0*V2+2.0*V3+V4)/6.0;
}


void draw_tree(qtree_node * node,int depth = 0)
{
    if(node==NULL)return;
    
    for(int i=0;i<4;i++)
        draw_tree(node->next[i],depth+1);
    
    glBegin(GL_LINE_LOOP);
      glColor3d(1.0/(1.0+depth),0.75/(1.0+depth),0.0);
      glVertex2d(node->q_min(0)/MAX_DRAW,node->q_min(1)/MAX_DRAW);
      glVertex2d(node->q_max(0)/MAX_DRAW,node->q_min(1)/MAX_DRAW);
      glVertex2d(node->q_max(0)/MAX_DRAW,node->q_max(1)/MAX_DRAW);
      glVertex2d(node->q_min(0)/MAX_DRAW,node->q_max(1)/MAX_DRAW);
    glEnd();
    
}

void draw(void)
{
    glClearColor(0,0,0,1);
    glClear(GL_COLOR_BUFFER_BIT);
    
    draw_tree(great_tree.root);
    glPointSize(1);
    glBegin(GL_POINTS);
      for(int i=0;i<X.n_cols;i++)
      {
          glColor3d(1.,1.,1.);
          glVertex2d(X(0,i)/MAX_DRAW,X(1,i)/MAX_DRAW);
      }
    glEnd();
    if(great_tree.root)
        {
          glPointSize(4);
          glBegin(GL_POINTS);
          glColor3d(1.0,0.75,0.0);
          glVertex2d(great_tree.root->c_mass(0)/MAX_DRAW,great_tree.root->c_mass(1)/MAX_DRAW);
          glEnd();
        }
    
    great_tree.insert(M,X);
    rk4(M,X,V,dt);
    glutSwapBuffers();
    glutPostRedisplay();
}

int main(int argc,char ** argv)
{
    srand(time(NULL));
    glutInit(&argc,argv);
    glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB);
    
    glutInitWindowPosition(50,25);
    glutInitWindowSize(700,700);
    
    glutCreateWindow("Barnes Hut N body Sim");
  
    mat A,B;
    A << 1.0 << 0.8 << endr
      << 0.7 << 0.1 << endr;
    B << 0.5 << 1.0 << endr
      << 0.4 << 2.0 << endr;
    X = 0.1*MAX*randn<mat>(2,N)-0.125*MAX;
    V = zeros<mat>(2,N);
    for(int i=0;i<N;i++)
    {
      V(0,i) += -0.25*(X(1,i));
      V(1,i) +=  0.25*(X(0,i));
    }
    M = 1e13*randu<vec>(2*N)+1e11;
    X = join_rows(X,0.1*MAX*randn<mat>(2,N)+0.25*MAX);
    V = join_rows(V,zeros<mat>(2,N));
    for(int i=N;i<2*N;i++)
    {
      V(0,i) += -0.25*(X(1,i));
      V(1,i) +=  0.25*(X(0,i));
    }
    glutDisplayFunc(draw);
    
    
    glutMainLoop();
    return 0;
}
