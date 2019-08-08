/* Ordered Line Integral Method OLIM-MID for finding the quasi-potential in 2D */
/* Computes the quasi-potential with respect to the equilibrium point that must be
 a mesh point. Its index is stored in the variable Iindex */
/* This version of olim2D does not do local factoring */
/* A guideline for choosing update factor K:
 (Nonlinear vector fields)
 NX = NY = 129:   K = 8
 NX = NY = 257:   K = 10
 NX = NY = 513:   K = 13
 NX = NY = 1025:  K = 16
 NX = NY = 2049:  K = 21
 NX = NY = 4097:  K = 26
 */
/* Modified for measurements */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/*--------- define constants ---------*/
#define sgn(a) ((a) == 0 ? 0 : ((a) > 0  ? 1 : -1 ))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define INFTY 1.0e+6
#define TOL 1.0e-12
// NX by NY is the mesh size
//#define NX 4097
//#define NY 4097
#define Nmax 4097
// K is the update factor
//#define K 26
// If you want only to compute the quasi-potential but do not want to shoot a MAP #define CH_SHOOT_MAP 'n'
// If you want to compute the quasi-potential AND shoot a MAP #define CH_SHOOT_MAP 'y'
// If you want to shoot a MAP and have already computed the quasi-potential #define CH_SHOOT_MAP 's'
#define CH_SHOOT_MAP 'n'
/*-------------------------------------*/

/*--------- define structures ---------*/
struct myvector {
    double x;
    double y;
};

struct index2 {
    long i;
    long j;
};

struct mymatrix {
    double a11;
    double a12;
    double a21;
    double a22;
};

// for triangle update hybrid nonlinear solver
struct mysol {
    double u; // min value
    double a; // minimizer: lambda or s in papers
    char c;   // indicator of success c = 'y' -- an inner point solution is found, otherwise c = 'n'
};

struct sol_info {
    char type; // solution type: 1 = 1ptupdate, 2 = 2ptupdate, 3 = 3ptupdate, 0 = initialization, 'n' = never reached
    long ind0; // ind0 gives the best one-pt-update
    double U1pt; // the best one-point update value of U
};
/*-------------------------------------*/

/*--------- predefine functions ---------*/
int main(void);
struct myvector myfield(struct myvector x); // b(x)
double exact_solution(struct myvector x);
// for initialization in 2D
void makeAmatrix(char chfield);
void param(void);
void ipoint(void);
double LinearQpot2D(struct myvector x);
// ordered line integral method
void olim(void);
// one-point update
double one_pt_update(long ind,long ind0);
// geometric action, midpoint rule
double geometric_action_line(struct myvector x0,struct myvector x1);
// hybrid nonlinear solver for triangle update
struct mysol hybrid_nonlin_solver(double u0,double u1,struct myvector x0,struct myvector x1,struct myvector b0,struct myvector b1,struct myvector x);
double myfun(double s,double du,struct myvector b0,struct myvector X01,struct myvector B10,struct myvector xmx0);
struct myvector getpoint(long ind); // index --> get coordinate [x,y]
struct index2 getindex2(long ind); // index --> get index2 (i,j)
// get N1 neighbors' indices
long get_neii_index(long d);
// get far neighbors' index shifts
long far_neighbors_index_list(long *farlist);
// linear algebra
struct myvector vec_difference(struct myvector v1,struct myvector v2); // v1 - v2
struct myvector vec_sum(struct myvector v1,struct myvector v2); // v1 + v2
struct myvector vec_lin_comb(struct myvector v1,struct myvector v2,double a,double b); // a*v1 + b*v2
struct myvector a_times_vec(struct myvector v,double a); // a*v
double dot_product(struct myvector a,struct myvector b); // (a \cdot b)
double length_vec(struct myvector x); // length of vector
double norm2squared(struct myvector x); // (x \cdot x)
// functions for the binary tree used for the heap sort of Considered points
void addtree(long ind); // add a new node to the binary tree
void updatetree(long ind); // update the binary tree
void deltree(void); // delete the root of the binary tree
/*-------------------------------------*/

/*--------- global variables ---------*/
char chfield; // choose the vector field b(x)
char ch_exact_sol; // indicates whether the exact solution is available or not
double fac; // fac = the ratio of the magnitudes of the rotational and the potential components
//const long nx1 = NX - 1, ny1 = NY - 1, NXY = NX*NY, nx2 = NX - 2, ny2 = NY - 2;
//const long KK = K*K;
long NX,NY,K,KK,nx1,ny1,NXY,nx2,ny2;
double hx,hy,hmax,dmaxsquared;
char ms[Nmax*Nmax]; // mesh point status: 0 = 'Unknown', 1 = 'Considered', 2 = "in Accepted Front", 3 = "Accepted"
double Uexact[Nmax*Nmax], g[Nmax*Nmax]; // function to be computed
struct sol_info solinfo[Nmax*Nmax];
long count = 0;  // for the binary tree: # of Considered points, initially there are no points in it
long pos[Nmax*Nmax];  // pos(index of mesh pt) = position in binary tree
long tree[Nmax*Nmax]; // tree(position in the tree) = index of mesh pt
long Iindex; // the index of the equilibrium point with respect to which the quasi-potential will be computed
struct myvector x_ipoint; // the coordinate of the equilibrium point
struct mymatrix Amatrix;  // Jacobian matrix of the vector field b(x)
// nearest neighbors' indices with L_inf distance = 1 and L1 distance <= 2
//const long neii[8] = {1, NX, -1, -NX, // 0 -- 3: L1 distance = 1 neighbors
//            NX+1, NX-1, -NX-1, -NX+1 // 4 -- 7: L1 distance = 2 neighbors
//            };
long Nfar; // the number of Considered points in the far neighborhood
long *farlist; // the list of index shifts for defining points in the far neighborhood
double XMIN,XMAX,YMIN,YMAX; // define the computational domain
long N1ptu,N2ptu,N2call; // for stats of your numerical solution
// names for output files
// to be done
/*-------------------------------------*/

/*--------- vector field b(x) ---------*/
struct myvector myfield(struct myvector x) {
    struct myvector v;
    double aux;
    
    switch( chfield ) {
        case 'd': // nonlinear example with analytical quasi-potential that is a double-well potential
            aux = 2.0*x.x*(x.x*x.x - 1.0); // 2(x^3 - x)
            v.x = -aux - fac*x.y;
            v.y = -x.y + fac*aux;
            break;
        case 'l': // a linear SDE
            v.x = -2.0*x.x - fac*x.y;
            v.y = 2.0*fac*x.x - x.y;
            break;
        case 'y': // stable equilibrium at the origin, limit cycle at r=1
            aux = 1.0 - x.x*x.x - x.y*x.y; // 1 - x^2 - y^2
            v.x = fac*x.y - x.x*aux;
            v.y = -fac*x.x - x.y*aux;
            break;
        default:
            printf("chfield = %c, please correct\n",chfield);
            exit(1);
            break;
    }
    return v;
}
/*-------------------------------------*/

/*--------- exact solution ---------*/
double exact_solution(struct myvector x) {
    double x2,aux;
    
    switch( chfield ) {
        case 'd':
            x2 = x.x*x.x;
            return x2*(x2 - 2.0) + x.y*x.y + 1.0;
        case 'l':
            return 2.0*x.x*x.x + x.y*x.y;
        case 'y':
            aux = x.x*x.x + x.y*x.y; // r^2
            return (aux <= 1.0) ? aux*(1.0 - 0.5*aux) : 0.5;
        default:
            return 0.0;
            break;
    }
}
/*-------------------------------------*/

/*--------- Jacobian matrix ---------*/
void makeAmatrix(char chfield) {
    double aux,aux1,aux2;
    
    switch( chfield ) {
        case 'd':
            aux = 2.0*(3.0*x_ipoint.x*x_ipoint.x - 1.0); // 2(3x^2 - 1)
            Amatrix.a11 = -aux;
            Amatrix.a12 = -fac;
            Amatrix.a21 = fac*aux;
            Amatrix.a22 = -1.0;
            break;
        case 'l':
            Amatrix.a11 = -2.0;
            Amatrix.a12 = -fac;
            Amatrix.a21 = 2.0*fac;
            Amatrix.a22 = -1.0;
            break;
        case 'y':
            aux = 1.0 - x_ipoint.x*x_ipoint.x - x_ipoint.y*x_ipoint.y; // 1 - x^2 - y^2
            aux1 = -2.0*x_ipoint.x; // -2x
            aux2 = -2.0*x_ipoint.y; // -2y
            Amatrix.a11 = - aux - x_ipoint.x*aux1;
            Amatrix.a12 = fac - x_ipoint.x*aux2;
            Amatrix.a21 = - fac - x_ipoint.y*aux1;
            Amatrix.a22 = - aux - x_ipoint.y*aux2;
            break;
        default:
            printf("in makeAmatrix( %c ): Provide the Jacobian matrix for the linearized system\n",chfield);
            exit(1);
            break;
    }
}
/*-------------------------------------*/

/*--------- parameters ---------*/
void param() {
    long ind;
    
    ch_exact_sol = (chfield == 'd' || chfield == 'l' || chfield == 'y') ? 'y' : 'n';
    
    switch( chfield ) {
        case 'd':
            Iindex = nx1/2 + NX*ny1/2;
            XMIN = -2.0; XMAX = 0.0;
            YMIN = -1.0; YMAX = 1.0;
            break;
        case 'l':
            Iindex = nx1/2 + NX*ny1/2;
            XMIN = -1.0; XMAX = 1.0;
            YMIN = -1.0; YMAX = 1.0;
            break;
        case 'y':
            Iindex = nx1/2 + NX*ny1/2;
            XMIN = -1.0; XMAX = 1.0;
            YMIN = -1.0; YMAX = 1.0;
            break;
        default:
            printf("chfield = %c, please correct\n",chfield);
            exit(1);
            break;
    }
    hx = (XMAX - XMIN)/nx1;
    hy = (YMAX - YMIN)/ny1;
    x_ipoint = getpoint(Iindex);
    printf("\nx_ipoint: Iindex = %li, %.4e, %.4e\n",Iindex,x_ipoint.x,x_ipoint.y);
    
    printf("in param()\n");
    hmax = max(hx,hy);
    dmaxsquared = max(KK*hmax*hmax,hmax*hmax*2.0) + 2.0e-16;
    makeAmatrix(chfield); // Jacobian matrix A: linearized system around the equilibrium point
    // generate status, quasi-potential, solinfo for all mesh points
    for( ind = 0; ind < NXY; ind++ ) {
        ms[ind] = 0; // status: Unknown
        g[ind] = INFTY;
        solinfo[ind].type = 'n'; // never reached
        solinfo[ind].U1pt = INFTY;
        // compute the exact solution if it is available
        if( ch_exact_sol == 'y' ) {
            Uexact[ind] = exact_solution(getpoint(ind));
        }
    }
}
/*-------------------------------------*/

/*--------- equilibrium point ---------*/
void ipoint() {
    long ind,ind0,n,m;
    struct myvector x;
    // nearest neighbors' indices with L_inf distance = 1 and L1 distance <= 2
    const long neii[8] = {1, NX, -1, -NX, // 0 -- 3: L1 distance = 1 neighbors
        NX+1, NX-1, -NX-1, -NX+1 // 4 -- 7: L1 distance = 2 neighbors
    };
    
    printf("in ipoint()\n");
    
    // status, quasi-potential, solinfo of the equilibrium point
    ms[Iindex] = 3; // status: Accepted
    g[Iindex] = 0.0;
    solinfo[Iindex].type = 0; // initialization
    // initialize 8 nearest neighbors of the equilibrium point
    for( n = 0; n < 8; n++ ) {
        ind = Iindex + neii[n];
        x = getpoint(ind);
        ms[ind] = 2; // status: Accepted Front
        g[ind] = LinearQpot2D(vec_difference(x,x_ipoint));
        solinfo[ind].type = 0; // initialization
    }
    // initialize all the Unknown nearest neighbors of the 8 nearest neighbors of the equilibrium point
    for( n = 0; n < 8; n++ ) {
        ind0 = Iindex + neii[n];
        for( m = 0; m < 8; m++ ) {
            ind = ind0 + neii[m];
            if( ms[ind] == 0 ) {
                x = getpoint(ind);
                ms[ind] = 1; // status: Considered
                g[ind] = LinearQpot2D(vec_difference(x,x_ipoint));
                solinfo[ind].type = 0; // initialization
                addtree(ind);
            }
        }
    }
}
/*-------------------------------------*/
    
/*--------- linear quasi-potential formula ---------*/
double LinearQpot2D(struct myvector x) {
    double A,B,C,aux1,aux2,aux;
    
    // compute the quasi-potential matrix Q = [A B; B C]
    aux1 = Amatrix.a21 - Amatrix.a12;
    aux2 = Amatrix.a11 + Amatrix.a22;
    aux = aux1*aux1 + aux2*aux2;
    aux1 *= aux2/aux; // beta
    aux2 *= aux2/aux; // alpha
    A = -(aux2*Amatrix.a11 + aux1*Amatrix.a21);
    B = -(aux2*Amatrix.a12 + aux1*Amatrix.a22);
    C = -(aux2*Amatrix.a22 - aux1*Amatrix.a12);
    
    return A*x.x*x.x + 2.0*B*x.x*x.y + C*x.y*x.y;
}
/*-------------------------------------*/

/*--------- ordered line integral method ---------*/
void olim() {
    long k,m,n,inew,ind,ind0,ind1,neii_index;
    double gtemp,gold;
    long Naf,Nc,AFneib[8],NCneib[8]; // Accepted Front N1 neighbors and new Considered neighbors of the new Accepted point
    long NAC = 0; // the number of Accepted/Accepted Front points
    struct mysol sol;
    struct myvector vec,b0,b1,v0,v1,vnew;
    struct index2 pnew;
    char pr = 'n'; // a logic variable: if pr == 'y', solution will be printed out in the terminal window
    // nearest neighbors' indices with L_inf distance = 1 and L1 distance <= 2
    const long neii[8] = {1, NX, -1, -NX, // 0 -- 3: L1 distance = 1 neighbors
        NX+1, NX-1, -NX-1, -NX+1 // 4 -- 7: L1 distance = 2 neighbors
    };
    
    printf("in olim()\n");
    
    while( count > 0 ) {
        inew = tree[1]; // index of the new point becoming Accepted Front
        vnew = getpoint(inew);  // [x,y] of this point
        pnew = getindex2(inew); // (i,j) of this point
        ms[inew] = 2; // Step 1: make it Accepted Front
        deltree();    // remove it from the binary tree
        NAC++; // increase the number of Accepted/Accepted Front points
        
        // Terminate computation if the boundary is reached
        if( pnew.i <= 1 || pnew.i >= nx2 || pnew.j <= 1 || pnew.j >= ny2 || g[inew] >= INFTY-1 ) {
            printf("The boundary is reached:\n");
            printf("%li Accepted points, (%li,%li) is accepted, g = %.4e\n",NAC,pnew.i,pnew.j,g[inew]);
            break;
        }
        
        // Inspect the neighbors of the new Accepted Front point inew
        // Step 2: Find nearest neighbors to be shifted to Accepted
        // Step 3: List nearest neighbors in N1 neighborhood that are Accepted Front to use for updates
        // Step 5: List Unknown neighbors that will be new Considered
        Naf = 0;
        Nc = 0;
        for( k = 0; k < 4; k++ ) {
            ind1 = inew + neii[k]; // N1 neighbor of the new Accepted point
            // update Accepted Front
            if( ms[ind1] == 2 ) {
                m = 0;
                for( n = 0; n < 8; n++ ) {
                    ind0 = ind1 + neii[n];
                    if( ms[ind0] < 2 ) m++;
                }
                if( m == 0 ) { // ind1 has no Considered neighbors
                    ms[ind1] = 3;
                }
                else {
                    AFneib[Naf] = ind1; // Accepted Front N1 neighbors of inew
                    Naf++;
                }
            }
            else if( ms[ind1] == 0 ) { // ind1 will be a new Considered point
                NCneib[Nc] = ind1; // new Considered neighbors of inew
                Nc++;
            }
        }
        for( k = 4; k < 8; k++ ) {
            ind1 = inew + neii[k]; // N2 neighbor of the new Accepted point
            // update Accepted Front
            if( ms[ind1] == 2 ) {
                m = 0;
                for( n = 0; n < 8; n++ ) {
                    ind0 = ind1 + neii[n];
                    if( ms[ind0] < 2 ) m++;
                }
                if( m == 0 ) { // ind1 has no Considered neighbors
                    ms[ind1] = 3;
                }
            }
            else if( ms[ind1] == 0 ) { // ind1 will be a new Considered point
                NCneib[Nc] = ind1; // new Considered neighbors of inew
                Nc++;
            }
        }
        
        // Step 4: Update existing Considered points in the far neighborhood of the new Accepted Front point inew
        for( k = 0; k < Nfar; k++ ) {
            ind = inew + farlist[k];
            if( ind < NXY && ind >= 0 ) {
                if( ms[ind] == 1 ) { // if the point ind is Considered --> to be updated
                    if( solinfo[ind].type > 0 ) { // do not update points initialized at the beginning using LinearQpot2D
                        vec = getpoint(ind);
                        if( norm2squared(vec_difference(vec,vnew)) <= dmaxsquared ) {
                            gold = g[ind];
                            gtemp = one_pt_update(ind,inew);
                            if( gtemp < g[ind] ) {
                                g[ind] = gtemp;
                                solinfo[ind].type = 1;
                                solinfo[ind].ind0 = inew;
                                solinfo[ind].U1pt = gtemp;
                            }
                            else if( gtemp < solinfo[ind].U1pt ) {
                                // the smallest 1pt update minimizer and value
                                solinfo[ind].ind0 = inew;
                                solinfo[ind].U1pt = gtemp;
                            }
                            b0 = myfield(vec_lin_comb(vnew,vec,0.5,0.5));
                            
                            // CASE 1: if the new AF point inew does not give the minimum of one-point update
                            if( inew != solinfo[ind].ind0 ) {
                                // there is only one triangle update to exercise
                                ind1 = solinfo[ind].ind0;
                                neii_index = get_neii_index(inew - ind1);
                                if( neii_index >= 0 ) { // ind1,inew are N1 neighbors
                                    // triangle update: Q2(ind1,inew,ind)
                                    v1 = getpoint(ind1);
                                    b1 = myfield(vec_lin_comb(v1,vec,0.5,0.5));
                                    sol = hybrid_nonlin_solver(g[ind1],g[inew],v1,vnew,b1,b0,vec);
                                    if( sol.c == 'y' && sol.u < g[ind] ) {
                                        // triangle update is successful
                                        g[ind] = sol.u;
                                        solinfo[ind].type = 2;
                                    }
                                } // end if( neii_index >= 0 )
                            } // end if( inew != solinfo[ind].ind0 )
                            
                            // CASE 2: if the new AF point inew gives the minimum of one-point update
                            else { // if( inew != solinfo[ind].ind0 )
                                // try all admissible triangle updates
                                for( m = 0; m < Naf; m++ ) {
                                    ind1 = AFneib[m]; // AF N1 neighbor of inew
                                    // triangle update: Q2(inew,ind1,ind)
                                    v1 = getpoint(ind1);
                                    b1 = myfield(vec_lin_comb(v1,vec,0.5,0.5));
                                    sol = hybrid_nonlin_solver(g[inew],g[ind1],vnew,v1,b0,b1,vec);
                                    if( sol.c == 'y' && sol.u < g[ind] ) {
                                        // triangle update is successful
                                        g[ind] = sol.u;
                                        solinfo[ind].type = 2;
                                    }
                                } // end for( m = 0; m < Naf; m++ )
                            } // end else
                            if( g[ind] < gold ) updatetree(ind);
                        } // end if( norm2squared(vec_difference(vec,vnew)) <= dmaxsquared )
                    } // end if( solinfo[ind].type > 0 )
                } // end if( ms[ind] == 1 )
            } // end if( ind < NXY && ind >= 0 )
        } // end for( k = 0; k < Nfar; k++ )
        
        // Step 5: Shift Unknown neighbors of the new Accepted Front point inew to Considered and update new Considered points
        for( k = 0; k < Nc; k++ ) {
            ind = NCneib[k];
            vec = getpoint(ind);
            ms[ind] = 1; // shift from Unknown to Considered
            // find minimal one-point update
            for( m = 0; m < Nfar; m++ ) {
                ind0 = ind + farlist[m];
                if( ind0 < NXY && ind0 >= 0 ) {
                    if( ms[ind0] == 2 ) {
                        v0 = getpoint(ind0);
                        if( norm2squared(vec_difference(vec,v0)) <= dmaxsquared ) {
                            gtemp = one_pt_update(ind,ind0);
                            if( gtemp < g[ind] ) {
                                g[ind] = gtemp;
                                solinfo[ind].ind0 = ind0;
                                solinfo[ind].U1pt = gtemp;
                            }
                        } // end if( norm2squared(vec_difference(vec,v0)) <= dmaxsquared )
                    } // end if( ms[ind0] == 2 )
                } // end if( ind0 < NXY && ind0 >= 0 )
            } // end for( m = 0; m < Nfar; m++ )
            // do triangle update with the minimizer of one-point update
            solinfo[ind].type = 1;
            ind0 = solinfo[ind].ind0; // minimizer of one-point update for ind
            v0 = getpoint(ind0);
            b0 = myfield(vec_lin_comb(v0,vec,0.5,0.5));
            // find minimal triangle update
            for( n = 0; n < 4; n++ ) {
                ind1 = ind0 + neii[n]; // N1 neighbor of the minimizer ind0
                if( ms[ind1] == 2 ) {
                    // triangle update: Q2(ind0,ind1,ind)
                    v1 = getpoint(ind1);
                    b1 = myfield(vec_lin_comb(v1,vec,0.5,0.5));
                    sol = hybrid_nonlin_solver(g[ind0],g[ind1],v0,v1,b0,b1,vec);
                    if( sol.c == 'y' && sol.u < g[ind] ) {
                        // triangle update is successful
                        g[ind] = sol.u;
                        solinfo[ind].type = 2;
                    }
                } // end if( ms[ind1] == 2 )
            } // end for( n = 0; n < 4; n++ )
            addtree(ind);
        } // end for( k = 0; k < Nc; k++ )
    } // end while( count > 0 )
}
/*-------------------------------------*/

/*--------- one-point update ---------*/
double one_pt_update(long ind,long ind0) {
    struct myvector x,x0;
    double gtemp;
    
    // one-point update: Q1(x0,x)
    x = getpoint(ind);   // Considered, to be updated
    x0 = getpoint(ind0); // Accepted Front, base
    gtemp = g[ind0] + geometric_action_line(x0,x);
    
    N1ptu++;
    
    return gtemp;
}
/*-------------------------------------*/

/*--------- geometric action line ---------*/
double geometric_action_line(struct myvector x0,struct myvector x1) {
    struct myvector l,b;
    
    // midpoint rule
    l = vec_difference(x1,x0); // integral from x0 to x1
    b = myfield(vec_lin_comb(x0,x1,0.5,0.5)); // b_m
    
    return length_vec(b)*length_vec(l) - dot_product(b,l);
}
/*-------------------------------------*/

/*--------- nonlinear 1D solver ---------*/
// triangle update: Q2(x0,x1,x)
// x0:    minimizer of the one-point update for x
// x0,x1: Accepted Front, bases, N1 neighbors
// x:     Considered, to be updated
// hybrid secant/bisection method -- find the root of f'(lambda)
struct mysol hybrid_nonlin_solver(double u0,double u1,struct myvector x0,struct myvector x1,struct myvector b0,struct myvector b1,struct myvector x) {
    double a = 0.0, b = 1.0, du = u1 - u0;
    double c, fa, fb, fc, d, fd, dd, df, dm, ds, t;
    struct myvector X01,B10,xs,xmx0;
    struct mysol sol;
    double NONLIN_TOL = 1.0e-6;
    long iter = 0, itermax = 100;
    
    X01 = vec_difference(x0,x1); // x0 - x1
    B10 = vec_difference(b1,b0); // b_m1 - b_m0
    xmx0 = vec_difference(x,x0); // x - x0
    
    c = a;
    fa = myfun(a,du,b0,X01,B10,xmx0);
    fb = myfun(b,du,b0,X01,B10,xmx0);
    fc = fa;
    
    N2call++;
    
    if( (fa > 0 && fb > 0 ) || (fa < 0 && fb < 0) ) {
        // root is not bracketed
        sol.c = 'n';
        sol.u = INFTY;
        return sol;
    }
    while( iter < itermax ) {
        if( fabs(fc) < fabs(fb) ) {
            t = c; c = b; b = t;
            t = fc; fc = fb; fb = t;
            a = c; fa = fc;
        }
        if( fabs(b - c) < NONLIN_TOL ) break;
        dm = 0.5*(c - b);
        df = fa - fb;
        
        if( fabs(df) < NONLIN_TOL ) ds = dm;
        else ds = -fb*(a - b)/df;
        if( (ds > 0 && dm < 0) || (ds < 0 && dm > 0) || (fabs(ds) > fabs(dm)) )
            dd = dm;
        else dd = ds;
        
        if( fabs(dd) < NONLIN_TOL ) dd = 0.5*sgn(dm)*NONLIN_TOL;
        
        d = b + dd;
        fd = myfun(d,du,b0,X01,B10,xmx0);
        if( fabs(fd) < NONLIN_TOL ) {
            b = d;
            break;
        }
        a = b; b = d; fa = fb; fb = fd;
        if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) ) {
            c = a; fc = fa;
        }
        iter++;
    }
    sol.c = 'y';
    sol.a = b; // minimizer: lambda
    xs = vec_difference(x0,a_times_vec(X01,b)); // x_lambda
    sol.u = u0 + b*du + geometric_action_line(xs,x);
    
    N2ptu++;
    
    return sol;
}
/*-------------------------------------*/

/*--------- derivative f'(lambda) ---------*/
double myfun(double s,double du,struct myvector b0,struct myvector X01,struct myvector B10,struct myvector xmx0) {
    double ls,lbs,beta;
    struct myvector xmxs,bs;
    
    xmxs = vec_sum(xmx0,a_times_vec(X01,s)); // x - x_lambda
    bs = vec_sum(b0,a_times_vec(B10,s));     // b_lambda
    ls = length_vec(xmxs);
    lbs = length_vec(bs);
    beta = lbs/ls;
    
    if( beta > TOL ) {
        return du + dot_product(xmxs,X01)*beta + dot_product(bs,B10)/beta - dot_product(X01,bs) - dot_product(B10,xmxs);
    }
    else {
        return du - dot_product(B10,xmxs);
    }
}
/*-------------------------------------*/

/*--------- compute coordinate[x,y] of index ---------*/
struct myvector getpoint(long ind) {
    struct myvector l;
    
    // ind = i + NX*j;
    l.x = hx*(ind%NX) + XMIN;
    l.y = hy*(ind/NX) + YMIN;
    return l;
}
/*-------------------------------------*/

/*--------- compute index2(i,j) of index ---------*/
struct index2 getindex2(long ind) {
    struct index2 m;
    
    // ind = i + NX*j;
    m.i = ind%NX; // 0 -- NX-1
    m.j = ind/NX; // 0 -- NY-1
    return m;
}
/*-------------------------------------*/

/*--------- N1 neighborhood ---------*/
/* long get_neii_index(long d) {
    switch( d ) {
        case 1:
            return 0;
            break;
        case NX:
            return 1;
            break;
        case -1:
            return 2;
            break;
        case -NX:
            return 3;
            break;
        default:
            return -1;
            break;
    }
} */
long get_neii_index(long d) {
    if( d == 1 ) return 0;
    else if( d == NX ) return 1;
    else if( d == -1 ) return 2;
    else if( d == -NX ) return 3;
    else return -1;
}
/*-------------------------------------*/

/*--------- far neighborhood ---------*/
long far_neighbors_index_list(long *farlist) {
    long Nfar = 0;
    long i,j,jmax;
    
    for( i = -K; i <= K; i++ ) {
        jmax = ceil(sqrt(KK - i*i));
        for( j = -jmax; j <= jmax; j++ ) {
            if( i != 0 || j != 0 ) {
                farlist[Nfar] = i + NX*j;
                Nfar++;
            }
        }
    }
    return Nfar;
}
/*-------------------------------------*/

/*--------- LINEAR ALGEBRA ---------*/
struct myvector vec_difference(struct myvector v1,struct myvector v2) {
    struct myvector v;
    
    v.x = v1.x - v2.x;
    v.y = v1.y - v2.y;
    return v;
}

struct myvector vec_sum(struct myvector v1,struct myvector v2) {
    struct myvector v;
    
    v.x = v1.x + v2.x;
    v.y = v1.y + v2.y;
    return v;
}

struct myvector vec_lin_comb(struct myvector v1,struct myvector v2,double a,double b) {
    struct myvector v;
    
    v.x = a*v1.x + b*v2.x;
    v.y = a*v1.y + b*v2.y;
    return v;
}

struct myvector a_times_vec(struct myvector v,double a) {
    struct myvector av;
    
    av.x = a*v.x;
    av.y = a*v.y;
    return av;
}

double dot_product(struct myvector a,struct myvector b) {
    return a.x*b.x + a.y*b.y;
}

double length_vec(struct myvector x) {
    return sqrt(norm2squared(x));
}

double norm2squared(struct myvector x) {
    return x.x*x.x + x.y*x.y;
}
/*-------------------------------------*/

/*--------- BINARY TREE ---------*/
// min-heap property: the value of each node is greater than or equal to the value of its parent, with the minimum-value element at the root
// the value of node = quasi-potential of node
// used for the heap sort of Considered points

// add a new node to the binary tree of Considered points
void addtree(long ind) {
    long loc, ptemp;
    long indp, indc;
    char ch;
    
    count++;
    tree[count] = ind;
    pos[ind] = count;
    if( count > 1 ) {
        loc = count;
        indc = tree[loc];
        indp = tree[loc/2];
        ch = ( g[indc] < g[indp] ) ? 'y' : 'n';
        // maintain the min-heap property
        while( ch == 'y' ) {
            // swap pos[indc] and pos[indp] (positions of child and parent)
            ptemp = pos[indc];
            pos[indc] = pos[indp];
            pos[indp] = ptemp;
            // swap index in tree[]
            tree[loc/2] = indc;
            tree[loc] = indp;
            loc = loc/2;
            if( loc > 1 ) {
                indc = tree[loc];
                indp = tree[loc/2];
                ch = ( g[indc] < g[indp] ) ? 'y' : 'n';
            }
            else ch='n';
        }
    }
}

// update the binary tree if the value of a node changes
void updatetree(long ind) {
    long loc, lcc;
    double g0;
    
    g0 = g[ind];
    loc = pos[ind];
    // maintain the min-heap property
    while( loc > 1 && g0 < g[tree[loc/2]] ) {
        // swap parent loc/2 and child loc in both tree[] and pos[]
        tree[loc] = tree[loc/2];
        pos[tree[loc]] = loc;
        loc = loc/2;
        tree[loc] = ind;
        pos[tree[loc]] = loc;
    }
    lcc = count;
    while( (loc*2 <= count && g0 > g[tree[loc*2]]) || (loc*2+1 <= count && g0 > g[tree[loc*2+1]]) ) {
        lcc = ( loc*2+1 <= count && g[tree[loc*2+1]] < g[tree[loc*2]] ) ? loc*2+1 : loc*2;
        // swap parent loc and child lcc in both tree[] and pos[]
        tree[loc] = tree[lcc];
        pos[tree[loc]] = loc;
        loc = lcc;
        tree[loc] = ind;
        pos[tree[loc]] = loc;
    }
}

// delete the root of the binary tree
void deltree() {
    long loc, ptemp, ind, lcc, ic, ic1, ic2, mind;
    char chd, ch='n'; // ch: a logic variable, if ch == 'y', which child being updated will be printed out in the terminal window
    
    mind = tree[1];
    pos[tree[1]] = 0;
    tree[1] = tree[count];
    pos[tree[1]] = 1;
    count--;
    loc = 1;
    ind = tree[1];
    lcc = 2*loc;
    if( lcc < count ) {
        ic1 = tree[lcc];   // left child
        ic2 = tree[lcc+1]; // right child
        if( (g[ind]) > (g[ic1]) || (g[ind]) > (g[ic2]) ) {
            if( (g[ic1]) <= (g[ic2]) ) {
                chd = 'l';
                ic = ic1;
            }
            else {
                chd = 'r';
                ic = ic2;
                lcc++;
            }
        }
        else chd = 'n';
    }
    else if( lcc == count ) {
        ic = tree[lcc];
        if( ch == 'y' ) printf("child: loc(%li) = %li, g = %.12e\n",ic,lcc,g[ic]);
        if( (g[ind]) > (g[ic]) ) {chd = 'l'; if( ch == 'y' ) printf("left\n");}
        else chd = 'n';
    }
    else chd = 'n';
    while( chd != 'n' ) {
        // swap pos[ic] and pos[ind] (positions of child and parent)
        ptemp = pos[ind];
        pos[ind] = pos[ic];
        pos[ic] = ptemp;
        // swap index in tree[]
        tree[loc] = ic;
        tree[lcc] = ind;
        loc = lcc;
        lcc = 2*loc;
        if( lcc < count ) {
            ic1 = tree[lcc];   // left child
            ic2 = tree[lcc+1]; // right child
            if( (g[ind]) > (g[ic1]) || (g[ind]) > (g[ic2]) ) {
                if( (g[ic1]) <= (g[ic2]) ) {
                    chd = 'l';
                    ic = ic1;
                }
                else {
                    chd = 'r';
                    ic = ic2;
                    lcc++;
                }
            }
            else chd = 'n';
        }
        else if( lcc == count ) {
            ic = tree[lcc];
            if( ch == 'y' ) printf("child: loc(%li) = %li, g = %.12e\n",ic,lcc,g[ic]);
            if( (g[ind]) > (g[ic]) ) {chd = 'l'; if( ch == 'y' ) printf("left\n");}
            else chd = 'n';
        }
        else chd = 'n';
    } // end while( chd != 'n' )
}
/*-------------------------------------*/

/*--------- main ---------*/
int main() {
    long i,j,k,ind,NAC,icase;
    double cpu,dd,errmax,erms,umax,urms;
    clock_t CPUbegin;
    double q1,q2,qc2;
    
    // measurement
    long p,Kopt;
    const long dK = 1; // step size for K
    const long pmin = 7, pmax = 12, Kmin = 1, Kmax = 40;
    FILE *fs;
    char fname[100];
    
    for( icase = 0; icase < 3; icase++ ) {
        if( icase == 0 ) { chfield = 'd'; fac = 10.0; }
        else if( icase == 1 ) { chfield = 'y'; fac = 4.0; }
        else { chfield = 'y'; fac = 10.0; }
        //sprintf(fname,"measurements/olim2D_chfield_%c.txt",chfield);
        //fs = fopen(fname,"w");
    
        if( CH_SHOOT_MAP == 'y' || CH_SHOOT_MAP == 'n' ) {
            printf("\nchfield = %c, fac = %g\n",chfield,fac);
        
            // measurement
            for( p = pmin; p <= pmax; p++) {
                NX = pow(2,p) + 1;
                NY = NX;
                // Kopt = 10 + 4*(p - 7); // for case 'l'
                nx1 = NX - 1; ny1 = NY - 1;
                nx2 = NX - 2; ny2 = NY - 2;
                NXY = NX*NY;
                sprintf(fname,"measurements/olim2D_chfield_%c_fac_%g_N_%li_testK.txt",chfield,fac,NX);
                fs = fopen(fname,"w");
                for( K = Kmin; K <= Kmax; K+=dK ) {
                    KK = K*K;
                    count = 0;
                    errmax = 0.0;
                    erms = 0.0;
                    umax = 0.0;
                    urms = 0.0;
                    N1ptu = 0;
                    N2ptu = 0;
                    N2call = 0;
                    param();
                    ipoint();
                    k = (2*K+1);
                    farlist = (long *)malloc(k*k*sizeof(long));
                    Nfar = far_neighbors_index_list(farlist);
                    CPUbegin = clock();
                    olim();
                    cpu = (clock()-CPUbegin)/((double)CLOCKS_PER_SEC);
                    printf("NX = %li, NY = %li, K = %li\ncputime of olim() = %g\n",NX,NY,K,cpu);
                    NAC = 0; // the number of Accepted/Accepted Front points
                    for( j = 0; j < NY; j++ ) {
                        for( i = 0; i < NX; i++ ) {
                            ind = i + NX*j;
                            if( ms[ind] < 2 ) g[ind] = INFTY;
                            else {
                                umax = max(umax,g[ind]);
                                urms += g[ind]*g[ind];
                                if( ch_exact_sol == 'y' ) {
                                    dd = fabs(g[ind] - Uexact[ind]);
                                    errmax = max(errmax,dd);
                                    erms += dd*dd;
                                }
                                NAC++;
                            } // end else
                        } // end for( i = 0; i < NX; i++ )
                    } // end for( j = 0; j < NY; j++ )
                    // print info about the computed solution
                    printf("# of Accepted points = %li, umax = %.4e\n",NAC,umax);
                    q1 = N1ptu;
                    q2 = N2ptu;
                    qc2 = N2call;
                    //printf("Per mesh point:\n<#1ptupdates> = %.2f, <#2ptupdates> = %.2f, <#2call> = %.2f\n",q1/NAC,q2/NAC,qc2/NAC);
                    if( ch_exact_sol == 'y' ) {
                        printf("ErrMax = %.4e, ERMS = %.4e\nNormalized ErrMax = %.4e, Normalized ERMS = %.4e\n",errmax,sqrt(erms/NAC),errmax/umax,sqrt(erms/urms));
                        //fprintf(fs,"%li\t%.4e\t%.4e\t%g\n",NX,errmax/umax,sqrt(erms/urms),cpu);
                        fprintf(fs,"%li\t%.4e\t%.4e\t%g\n",K,errmax/umax,sqrt(erms/urms),cpu);
                    }
                } // end for K
            } // end for p
        } // end if( CH_SHOOT_MAP == 'y' || CH_SHOOT_MAP == 'n' )
        else {
            printf("Currently, CH_SHOOT_MAP = %c. Define CH_SHOOT_MAP to be 'y', 'n', or 's' and run again\n",CH_SHOOT_MAP);
        }
	} // end for( icase = 0; icase < 3; icase++ )
    fclose(fs);
    return 0;
}
/*-------------------------------------*/
