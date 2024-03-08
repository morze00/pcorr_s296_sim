/////////////////////////////////////////////////////
//This code was written by Felix Wamers, with
//the help and support of several members of
//the LAND group, and for the randomiser part
//to a great extent by student Martin Riedel.
/////////////////////////////////////////////////////

#include<vector>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <TMath.h>
#include "TRandom1.h"
#include "xb_list_energy.hh"
#include "xbtools.h"


#define Rd 49.86932985	//the reference radius for the lengthes in the technical drawings (cm)
#define Rg 35.0		//the reference radius for the lengthes in the geant trafo files (cm)

bool inside_extended_crystal(int cr_check, Float_t cr_th, Float_t cr_ph, double theta_sid, double phi_sid)
{
	//Float_t	phi_cr_angle = 17;
	//Float_t	theta_cr_angle = 17;

	Float_t	phi_cr_angle = 30;
	Float_t	theta_cr_angle = 30;

	if (cr_check == 60 || cr_check == 61 || cr_check == 62 || cr_check == 90 || cr_check == 91 || cr_check == 92)
	{
		phi_cr_angle = 50;
		//phi_cr_angle = 30;
	}

	if (cr_check == 89 || cr_check == 111 || cr_check == 112 || cr_check == 113 || cr_check == 114 || cr_check == 93 || cr_check == 80 || cr_check == 59 || cr_check == 41 || cr_check == 42 || cr_check == 43 ||  cr_check== 63)
	{
		phi_cr_angle = 50;
		//phi_cr_angle = 25;
		//phi_cr_angle = 15;
	}

	if (cr_check == 4 || cr_check == 3 || cr_check == 156 || cr_check == 157)
	{
		theta_cr_angle = 50;
		//theta_cr_angle = 30;
	}

	if (cr_check == 10 || cr_check == 11 || cr_check == 12 || cr_check == 13 || cr_check == 14 || cr_check == 155 || cr_check == 144 ||cr_check == 145 || cr_check == 146 || cr_check == 147 )
	{
		//theta_cr_angle = 25;
		theta_cr_angle = 50;
		//theta_cr_angle = 15;
	}

	// PHI RANGE
	Float_t cr_phi_low = cr_ph - phi_cr_angle; 
	if(cr_phi_low < 0)  cr_phi_low = cr_phi_low + 360;

	Float_t cr_phi_up = cr_ph + phi_cr_angle;
	if(cr_phi_up > 360)  cr_phi_up = cr_phi_up - 360;

	// THETA RANGE
	Float_t cr_theta_low = cr_th - theta_cr_angle; 
	if(cr_theta_low < 0)  cr_theta_low = 0;

	Float_t cr_theta_up = cr_th + theta_cr_angle;

	if(cr_phi_up < cr_phi_low) // special case for crystal
	{
		if(((phi_sid<cr_phi_up) && (theta_sid>cr_theta_low) && (theta_sid<cr_theta_up))
				|| ((phi_sid>cr_phi_low) && (theta_sid>cr_theta_low) && (theta_sid<cr_theta_up)))
			return true;
	}

	else if(cr_phi_up > cr_phi_low) // special case for crystal
	{
		if((phi_sid>cr_phi_low) && (phi_sid<cr_phi_up) && (theta_sid>cr_theta_low) && (theta_sid<cr_theta_up))
			return true;
	}
	//else
	return false;

}
struct crys mycr[162];	//declare the XB crystals...

bool data_loaded=false;
double corn_x[4][7];	//x-coords of corner points
double corn_y[4][7]; //y-coords of corner points

//degree rad conversion
double d2r(double deg)
{
	double rad = deg * TMath::Pi() / 180.0;
	return rad;
}
double r2d(double rad)
{
	double deg = rad * 180.0 / TMath::Pi();
	return deg;
}

// stretch/shrink a vector by multiplying with a scalar
struct vect vect_m_scal(struct vect a, double scale)
{
	struct vect result;
	result.x = scale * a.x;
	result.y = scale * a.y;
	result.z = scale * a.z;
	return result;
}

// add two vectors
struct vect vect_a_vect(struct vect a, struct vect b)
{
	struct vect result;
	result.x = a.x + b.x;
	result.y = a.y + b.y;
	result.z = a.z + b.z;
	return result;
}

// subtract two vectors
struct vect vect_s_vect(struct vect a, struct vect b)
{
	struct vect result;
	result.x = a.x - b.x;
	result.y = a.y - b.y;
	result.z = a.z - b.z;
	return result;
}
// multiply two vectors
struct vect vect_m_vect(struct vect a, struct vect b)
{
	struct vect result;
	result.x = a.x * b.x;
	result.y = a.y * b.y;
	result.z = a.z * b.z;
	return result;
}
// right-multiply a vector to a matrix
struct vect matr_m_vect(struct matr m, struct vect v)
{
	struct vect result;
	result.x = m.xx * v.x + m.yx*v.y + m.zx*v.z;	
	result.y = m.xy * v.x + m.yy*v.y + m.zy*v.z;
	result.z = m.xz * v.x + m.yz*v.y + m.zz*v.z;
	return result;
}
//multiply matrix n from the right to matrix m
struct matr matr_m_matr(struct matr m, struct matr n)
{
	struct matr result;
	result.xx = m.xx*n.xx + m.yx*n.xy + m.zx*n.xz;
	result.yx = m.xx*n.yx + m.yx*n.yy + m.zx*n.yz;
	result.zx = m.xx*n.zx + m.yx*n.zy + m.zx*n.zz;

	result.xy = m.xy*n.xx + m.yy*n.xy + m.zy*n.xz;
	result.yy = m.xy*n.yx + m.yy*n.yy + m.zy*n.yz;
	result.zy = m.xy*n.zx + m.yy*n.zy + m.zy*n.zz;

	result.xz = m.xz*n.xx + m.yz*n.xy + m.zz*n.xz;
	result.yz = m.xz*n.yx + m.yz*n.yy + m.zz*n.yz;
	result.zz = m.xz*n.zx + m.yz*n.zy + m.zz*n.zz;

	return result;
}

// calculate the determinant of a 3x3 matrix
double getDet(struct matr m)
{
	double D = (m.xx * m.yy * m.zz) + (m.yx * m.zy * m.xz) + (m.zx * m.xy * m.yz) - (m.xz * m.yy * m.zx) - (m.yz * m.zy * m.xx) - (m.zz * m.xy * m.yx);
	return D;
}

// calculate the inverse of a 3x3 matrix
struct matr matr_inv(struct matr m)
{
	double d = getDet(m);
	struct matr result;

	result.xx = (1/d)*((m.yy*m.zz)-(m.zy*m.yz));
	result.yx = (1/d)*((m.zx*m.yz)-(m.yx*m.zz));
	result.zx = (1/d)*((m.yx*m.zy)-(m.zx*m.yy));

	result.xy = (1/d)*((m.zy*m.xz)-(m.xy*m.zz));
	result.yy = (1/d)*((m.xx*m.zz)-(m.zx*m.xz));
	result.zy = (1/d)*((m.zx*m.xy)-(m.xx*m.zy));

	result.xz = (1/d)*((m.xy*m.yz)-(m.yy*m.xz));
	result.yz = (1/d)*((m.yx*m.xz)-(m.xx*m.yz));
	result.zz = (1/d)*((m.xx*m.yy)-(m.yx*m.xy));

	return result;
}

// transform from the Geant internal system to the lab system:
struct vect geant2lab(struct vect vgeant)
{
	struct vect vlab;
	vlab.y = - vgeant.x;
	vlab.x = vgeant.y;
	vlab.z = - vgeant.z;
	return vlab;
}

// transform from the lab system to the Geant internal system:
struct vect lab2geant(struct vect vlab)
{
	struct vect vgeant;
	vgeant.y = vlab.x;
	vgeant.x = - vlab.y;
	vgeant.z = - vlab.z;
	return vgeant;
}

//transform from kartesian to sperical coordinates
struct vectsp kart2spher(struct vect k)
{
	struct vectsp spher;
	spher.r = sqrt(k.x*k.x + k.y*k.y + k.z*k.z);
	spher.t = r2d((TMath::Pi()/2)-atan(k.z/sqrt(k.x*k.x + k.y*k.y)));
	spher.p = r2d(atan2(k.y, k.x));
	//atan2 returns phi values between -180 and +180, need to shift that to 0 to 360
	if (spher.p<0)
	{	
		spher.p = 360 + spher.p;
	}

	return spher;
}

//calculate the length of a vector
double get_vlength(struct vect v)
{
	double vlength;
	vlength = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
	return vlength;
}

//perform an initial rotation of a random vector around the z-axis
struct vect rot_z (struct vect v, int shape)
{
	struct vect result;
	double ang = 0;
	switch(shape)
	{
		case 0: ang = 0; break;
		case 1: ang = 270; break;
		case 2: ang = 270; break;
		case 3: ang = 0; break;
		default: printf("ERROR"); break;
	}
	double angg = d2r(ang);
	result.x = v.x * cos(angg) - v.y * sin(angg);
	result.y = v.x * sin(angg) + v.y * cos(angg);
	result.z = v.z * 1;

	return result;
}

//Create a random point in x-y
void rpoint(double *xtp, double *ytp)
{
	//xtp=(((double)rand())/((double)(RAND_MAX)+(double)(1)))*40-20;
	//ytp=(((double)rand())/((double)(RAND_MAX)+(double)(1)))*40-20;
	TRandom1 myRnd; 			// try this
	(*xtp) = myRnd.Rndm()*40 - 20;
	(*ytp) = myRnd.Rndm()*40 - 20;
	return;
}

//Check if a random point's x-projection is within the x-span of a limiting edge of the polygon
//void xrangenew(int ts, double xtp, char *ins)	//x is the list of x coords of the polygon, xtp the x coord of the point to be tested
void xrangenew(int ts, double *x, double xtp, char *ins)	//x is the list of x coords of the polygon, xtp the x coord of the point to be tested
{
	int i;
	int c=6;
	
	if(ts == 0)	//if Shape A is choosen
	{
		c=5;
	}
	for(i=0; i<c; i++)
	{
		//if( ( (xtp<x[ts][i])&&(xtp>x[ts][i+1]) ) || ( (xtp>x[ts][i])&&(xtp<x[ts][i+1]) ) )
		if( ( (xtp<x[i])&&(xtp>x[i+1]) ) || ( (xtp>x[i])&&(xtp<x[i+1]) ) )
		{
			ins[i]=1;
		}
	}
}

//Count how many times a line from that point upwards would cross a line of the polygon
//int math(int count, int ts, double xtp, double ytp, char *ins)
int math(int count, double *x, double *y, double xtp, double ytp, char *ins)
{
	int i;
	int j = count;
	for(i=0; i<6; i++)
	{
		if(ins[i]==1)
		{
			//double m=(y[ts][i+1]-y[ts][i])/(x[ts][i+1]-x[ts][i]);
			//double n=y[ts][i]-m*(x[ts][i]);
			double m=(y[i+1]-y[i])/(x[i+1]-x[i]);
			double n=y[i]-m*(x[i]);
			double yline=m*xtp+n;
			if(yline > ytp) //ylint==ytpint for edge testing
			{
				j++;
			}
			ins[i]=0;
		}
	}
	return j;
}

//Is a number even or odd?
int evenodd(int count)
{
	int isinside=0;
	if(count % 2 == 0)
	{
		isinside=0;
	}
	else
	{
		isinside=1;
	}

	return (isinside);
}

//Convert A-D to 0-3.
int shape_t_n(char shape)
{
	int number;
	switch(shape)
     	{
		case 'A' : number=0;
			break;
		case 'B' : number=1;
			break;
		case 'C' : number=2;
			break;
		case 'D' : number=3;
			break;
		default  :  
                    printf( "ERROR reading shape form i.e. (%c) \n", shape );
    	 }	
	return number;
}

//check if point in polygon, from http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
int pnpoly(int nvert, double *vertx, double *verty, double testx, double testy)
{
	int i, j, c = 0;
	for (i = 0, j = nvert-1; i < nvert; j = i++)
	{
		if ( ((verty[i]>testy) != (verty[j]>testy)) &&
				(testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
			c = !c;
	}
	return c;
}

//############### important functions ##################

//Read-in forms of crystal shapes
void read_corners()
{
	FILE *corners;
	int a,b;
        corners = fopen("/Users/vpanin/Desktop/Experiment/QFS/s296/p_corr/xbtools/cb_corners.dat","r");
	if (corners==NULL)
	{
		printf("File not successfully opened!\n");
		return;
	}
	for(a = 0; a < 4; a++)
	{
		for(b = 0; b < 7; b++)
		{
			char dummy;
			fscanf(corners,"%c\t%lf\t%lf\n",&dummy,&corn_x[a][b],&corn_y[a][b]);
		}
	}
	fclose(corners);
}

void loadData()
{
	printf("loading data\n");
	read_corners();
	//printf("in the loop: %d\n",kr);
	//printf("in the loop: %lf\n",TMath::Pi());

	int i, dummy, dummy2;

	FILE *params;
	FILE *params_rot;
	FILE *params_pos;
	//params = fopen("input/cb_geometry.dat","r");
	params = fopen("/Users/vpanin/Desktop/Experiment/QFS/s296/p_corr/xbtools/cb_geometry.dat","r");
	if (params==NULL)
	{
		printf("File not successfully opened!\n");
		return;
	}
	//params_rot = fopen("input/cb_rotations.dat","r");
	params_rot = fopen("/Users/vpanin/Desktop/Experiment/QFS/s296/p_corr/xbtools/cb_rotations.dat","r");
	if (params_rot==NULL)
	{
		printf("File not successfully opened!\n");
		return;
	}
	//params_pos = fopen("input/cb_positions.dat","r");
	params_pos = fopen("/Users/vpanin/Desktop/Experiment/QFS/s296/p_corr/xbtools/cb_positions.dat","r");
	if (params_pos==NULL)
	{
		printf("File not successfully opened!\n");
		return;
	}

	//read in rotation angles and translation vectors for crystals, create rotation matrices.
	for (i=1; i<163; i++)
	{
		//fscanf(params,"%d\t%lf\t%lf\t%c\t%d\t%d\t%d\t%d\t%d\t%d",&mycr[i-1].number,&mycr[i-1].theta, &mycr[i-1].phi,&mycr[i-1].shape,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy);
		fscanf(params,"%d\t%lf\t%lf\t%c\t%d\t%d\t%d\t%d\t%d\t%d",&mycr[i-1].number,&mycr[i-1].theta, &mycr[i-1].phi,&mycr[i-1].shape,&mycr[i-1].nb[0],&mycr[i-1].nb[1],&mycr[i-1].nb[2],&mycr[i-1].nb[3],&mycr[i-1].nb[4],&mycr[i-1].nb[5]);
		fscanf(params_rot,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",&dummy,&mycr[i-1].tx,&mycr[i-1].px,&mycr[i-1].ty,&mycr[i-1].py,&mycr[i-1].tz,&mycr[i-1].pz);
		mycr[i-1].m_rot.xx = sin(d2r(mycr[i-1].tx)) * cos(d2r(mycr[i-1].px));
		mycr[i-1].m_rot.yx = sin(d2r(mycr[i-1].tx)) * sin(d2r(mycr[i-1].px));
		mycr[i-1].m_rot.zx = cos(d2r(mycr[i-1].tx));
		mycr[i-1].m_rot.xy = sin(d2r(mycr[i-1].ty)) * cos(d2r(mycr[i-1].py));
		mycr[i-1].m_rot.yy = sin(d2r(mycr[i-1].ty)) * sin(d2r(mycr[i-1].py));
		mycr[i-1].m_rot.zy = cos(d2r(mycr[i-1].ty));
		mycr[i-1].m_rot.xz = sin(d2r(mycr[i-1].tz)) * cos(d2r(mycr[i-1].pz));
		mycr[i-1].m_rot.yz = sin(d2r(mycr[i-1].tz)) * sin(d2r(mycr[i-1].pz));
		mycr[i-1].m_rot.zz = cos(d2r(mycr[i-1].tz));	
		fscanf(params_pos,"%d\t%lf\t%lf\t%d\t%lf\n",&dummy,&mycr[i-1].v_tra.x,&mycr[i-1].v_tra.y,&dummy2,&mycr[i-1].v_tra.z);
	}

	fclose(params);	
	fclose(params_rot);
	fclose(params_pos);
	data_loaded=true;
	return;
}	//laodData()

int insideCrystal(int kr, double theta, double phi)
{
	if (!data_loaded)
	{
		loadData();
		data_loaded=true;
	}//end loading data if necessary

	//the shape of the crystal to be checked
	int inside = 0;
	int inside1 = 0;
	int count=0;				//Initialising variables
	int ts = shape_t_n(mycr[kr-1].shape);
	char ins[7] = {0,0,0,0,0,0,0};	//marker for a point being within the range of a limiting polygon edge
	//variables for storing the theta and phi of the transformed corner points	
	double ct[7], cp[7];
	double xp,yp;
	
	struct vect scaled_t = vect_m_scal(mycr[kr-1].v_tra, Rg/get_vlength(mycr[kr-1].v_tra));
	for (int i=0; i<7; i++)	//loop oer the corners
	{
		//the vector representing a polygon point
		struct vect shape = {corn_x[ts][i],corn_y[ts][i],0};
		//initial rotation of a (vector in a crystal) shape, because polygon definition in geant was rotated relative to how we defined them in our hardcoding.
		struct vect rot_shape = rot_z(shape, ts);
		//swap x and -x, because polygon points are hardcoded with a right-pointing xaxis, whereas the Lab system x-axis points left
		xp = -rot_shape.x;	
		yp = rot_shape.y;

		//Create the vector (after an possible initial shift and rotation) to be put into the right crystal's place
		struct vect tp = {xp*Rg/Rd, yp*Rg/Rd, 0};
		//Do the (complicated) transformation of individual rotation, translation, and finally a back-transformation into our Lab system:
		struct vect tnp = geant2lab(vect_s_vect(matr_m_vect(matr_inv(mycr[kr-1].m_rot),lab2geant(tp)),scaled_t));
		//Transform into spherical coords, and fill the array.
		struct vectsp cpol = kart2spher(tnp);
		ct[i] = cpol.t;
		cp[i] = cpol.p;
		
	} //corner loop
	//need to correct for crystals crossing the 0/360 degrees line
	if (phi > 0.0 && phi < 90.0)
	{
		for (int i=0; i<7; i++)
		{
			if (cp[i] > 270.0 && cp[i] < 360.0)
			{
				cp[i] = cp[i] - 360.0;
			}
		};
	}
	else if (phi > 270.0 && phi < 360.0)
	{
		for (int i=0; i<7; i++)
		{
			if (cp[i] > 0.0 && cp[i] < 90.0)
			{
				cp[i] = cp[i] + 360.0;
			}
		};
	}
	
	//xrangenew(ts,ct,theta,ins);
	//count = math(count,ct,cp,theta,phi,ins);
	//inside1 = evenodd(count);
	//printf("Is inside1 crystal %d: %d\n",kr,inside1);

	inside = pnpoly(7, ct, cp, theta, phi);
	printf("Is inside crystal %d: %d\n\n",kr,inside);

	if (!inside)
	{
		printf("Crystal %d; theta: %lf; phi: %lf\n",kr,theta,phi);
	}
	
	return inside;
}

void randomiser(int kr, double *theta, double *phi)
{
	//static bool data_loaded=false;
	//static struct crys mycr[162];
	
	if (!data_loaded)
	{
		loadData();
		data_loaded=true;
	}//end loading data if necessary

	//Get the shape of the crystal to be randomised
	int ts=shape_t_n(mycr[kr-1].shape);
	int count=0;				//Initialising variables
	int isinside = 0;			//for randomising area
	char ins[7] = {0,0,0,0,0,0,0};	//marker for a point being within the range of a limiting polygon edge
	double xtp,ytp;
	double xp,yp;

	//Scale the translation vectors to be equally long, i.e. Rg, 35cm.
	struct vect scaled_t = vect_m_scal(mycr[kr-1].v_tra, Rg/get_vlength(mycr[kr-1].v_tra));
	//Create new random points until the first one is inside the shape:	
	do				
	{
		rpoint(&xtp,&ytp);			//Finding whether a random point is inside a shape
		isinside = pnpoly(7, corn_x[ts], corn_y[ts], xtp, ytp);
		//xrangenew(ts,xtp,ins);
		//xrangenew(ts,corn_x[ts],xtp,ins);
		//count = math(count,corn_x[ts],corn_y[ts],xtp,ytp,ins);
		//isinside=evenodd(count);
		//count=0;
	}while(isinside < 1);
	
	//The random point in the crystal x/y plane, y pointing up, x right.
	struct vect rpoint = {xtp,ytp,0};

	//pre-rotation, because of the flip between the drawing, and the ancient G3 definition
	struct vect rot_rpoint = rot_z(rpoint, ts);
	//have x-axis to point to the left, i.e. Lab system
	xp = -rot_rpoint.x;
	yp = rot_rpoint.y;

	//Create the vector (after an possible initial shift and rotation) to be put into the right crystal's place
	//Scale the point with respect to the translation vectors, and the numbers from the drawings:
	struct vect tp = {xp*Rg/Rd, yp*Rg/Rd, 0};
	//Do the (complicated) transformation of individual rotation, translation, and finally a back-transformation into our Lab system:
	struct vect tnp = geant2lab(vect_s_vect(matr_m_vect(matr_inv(mycr[kr-1].m_rot),lab2geant(tp)),scaled_t));

	//printf("%lf\t%lf\t%lf\t kartesian coordinates of crystal %d\n",tnp.x,tnp.y,tnp.z,kr);
	struct vectsp polar = kart2spher(tnp);
	//printf("%lf\t%lf\t%lf\t polar coordinates of crystal %d\n",polar.r,polar.t,polar.p,kr);
	
	//Return the values asked for, theta and phi of the random and transformed point.
	(*theta) = polar.t;
	(*phi) = polar.p;
	//printf("in the loop: %lf\t%lf\n",(*theta),(*phi));

	return;
}

double dphi(double p1, double p2) {
	double deltaphi = 0.0;
	deltaphi = TMath::Abs(p2 - p1);
	if (deltaphi>180) deltaphi=360-deltaphi;
	return deltaphi;
}

double opang(double t1, double p1, double t2, double p2) {
	double opa = 0.0;
	double rt1,rp1,rt2,rp2;
	rt1 = d2r(t1);
	rp1 = d2r(p1);
	rt2 = d2r(t2);
	rp2 = d2r(p2);
	opa = sin(rt1)*sin(rt2)*cos(rp2-rp1) + cos(rt1)*cos(rt2) ;
	return r2d(acos(opa));
}


/*Extended algorithm. Performs ultimate summing of the crystals
# Gamma-energies of neighbours are added to the proton cluster. Threshold for energies of gammas is used (see definitions above)
# in order to disentangle between low energy proton signal and high-energy gamma when both are present.
# Found crystals which belong to one cluster are stored in the separate xb_list to calculate multiplicity in the cluster
*/
void adbp_all(std::list<xble> *asorted, std::list<xble> *psorted, std::list<xble> *gsorted) {
	std::list<xble>::iterator it1,it2,c_it,itp,itg;

	itp = psorted->begin();
	itg = gsorted->begin();
	int central = 0;
	std::list<xble> cluster; // to store crystals of the cluster

	//Go through the sorted list, from highest to lowest energy
	for (it1 = asorted->begin(); it1 !=asorted->end(); it1++){

		cluster.clear(); //reset for the next cluster
		//Type of signal in this crystal/cluster
		bool proton = false;
		bool gamma  = false;

		if (it1->is_proton_4pi()){	
			proton=true;
			central = it1->_number;
			cluster.push_back(*it1); //store this crystal in the list for the cluster 
			psorted->push_back(*it1);//store this crystal in the list for protons 
			itp = psorted->end();
			itp--;
		}
		else if (it1->is_gamma()){
			gamma=true;
			central = it1->_number;
			cluster.push_back(*it1);			
			gsorted->push_back(*it1);//store this crystal in the list for protons 
			itg = gsorted->end();
			itg--;
		}
		else {
			//	std::cout <<" \n ERROR!!!\n";
			//it1->print();
			it1 = asorted->erase(it1);
			it1--;	
			continue;
		}

		int j = 0;
		it2=it1;
		it2++;//start checking next crystal
		//std::cout << "\n~~~~~ First check for the central crystal " << central;
		
		//Inner loop over global list to find neghbours
		for (it2=it2 ; it2 != asorted->end(); it2++){
			j = it2->_number; //check this crystal from the list
			int k = 0;//neighbouring crystal from the list of neighbours
			//loop over the neighbours of that crystal
			for (int nb=0; nb<6; nb++){
				k = mycr[central-1].nb[nb]; // check this neighbour
				if (j==k){// Found neighbour of the central crystal in the globall list
					// Proton signal in both  crystals
					if (proton && it2->is_proton_4pi())
					{
						cluster.push_back(*it2);
						//neighbour without proton energy
						//if(isnan(itp->_energy_p) || isnan(it2->_energy_p))
						if(!isnan(itp->_energy_p) && !isnan(it2->_energy_p))
							itp->_energy_p += it2->_energy_p;
						//std::cout<<"\nNAN energy in proton branch";
						//std::cout<<" ";
						//else
						//itp->_energy_p += it2->_energy_p;
						itp->_mul++;
					}
					//Gamma signal in  both crystals
					else if (gamma && it2->is_gamma())
					{
						cluster.push_back(*it2);
						itg->_energy_g += it2->_energy_g;
						itg->_mul++;
						//it2 = asorted->erase(it2);
						//it2--;//go back to the element before the just erased one
					}
					//Proton in central crystal and gamma in the neighbour
					else if(proton && it2->is_gamma())
					{
						cluster.push_back(*it2);
						itp->_energy_p += it2->_energy_g;
						itp->_mul++;
						//it2 = asorted->erase(it2);
						//it2--;//go back to the element before the just erased one
					}

					it2 = asorted->erase(it2);
					it2--;//go back to the element before the just erased one
					//	else std::cout<< "\n\t\t NOT ADDED!!!!!!!!!!!!!!! ";
				} // end if(j==k)
			}// end looping neighbours
		} // end looping global list to find first neghbours

		//Second check for next next neighbours
		//std::cout << "\n#########FOUND CLUSTER with the central crystal " << central << "\n";
		//for(c_it=cluster.begin(); c_it != cluster.end(); c_it++){
		//	c_it->print();
		//}

		//if(cluster.size()<2) continue;

		//std::cout << "\nSecond check...";

		c_it = cluster.begin();
		c_it++; // since the first element is the cluster's center (already looped)

		for(c_it=c_it; c_it != cluster.end(); c_it++){ //loop the cluster
			j = c_it->_number;

			//std::cout << "\n\tCheck neighbours of the crystal "<< j <<  "\n";
			//c_it->print();

			for(int n=0; n<6; n++){ // loop the neighbours

				int next = mycr[j-1].nb[n]; // check if this neighbour is in the list
				//std::cout << "\n\t\tCheck neighbour crystal " << next << ".....";
				if(central==next) continue; //skip the central crystal in the cluster

				it2=it1;
				it2++;// as before

				for (it2=it2; it2 != asorted->end(); it2++){// look for this neighbours in the list

					if(it2->_number!=next) continue;

					//std::cout<<"\n\t\t\tFOUND EXTRA NEIGHBOUR IN THE LIST !!! crystal "<< it2->_number<<"\n";
					//it2->print();
					//std::cout<<" \t\t\tLet's compare...";

					if (proton && it2->is_proton_4pi()) 
					{	
						//cluster.push_back(*it2);
						if(!isnan(it2->_energy_p))
							//if(isnan(it2->_energy_p))
							//std::cout<<"NAN energy in proton branch!!!";
							//std::cout<<" ";
							//else
							itp->_energy_p += it2->_energy_p;
						itp->_mul++;
						//it2 = asorted->erase(it2);
						//it2--;
						//std::cout << "\n\t\t\t ******** added PROTON to PROTON ********";
					}
					else if (proton && it2->is_gamma())
					{
						//cluster.push_back(*it2);
						itp->_energy_p += it2->_energy_g;
						itp->_mul++;
						//it2 = asorted->erase(it2);
						//it2--;
						//std::cout << "\n\t\t\t ******** added GAMMA to PROTON *********";
					}
					else if (gamma && it2->is_gamma())
					{
						//cluster.push_back(*it2);
						itg->_energy_g += it2->_energy_g;
						itg->_mul++;
						//it2 = asorted->erase(it2);
						//it2--;
						//std::cout << "\n\t\t\t ******** added GAMMA to GAMMA ********";
					}

					it2 = asorted->erase(it2);
					it2--;

					//else std::cout << "\n\t\t\t ******** NOT ADDED, GO TO NEXT CRYSTAL ...";
				}//loop global list to find neighbours
			}// looping neighbours
		}// looping cluster

		//std::cout<< "\n\tFINAL CLUSTER: \n";
		//print_all(cluster);

	}//Determined the proton clusters. When central crystal gamma energy not in overflow, re-declare it a gamma cluster...
	return;
}

/*Extended algorithm. Performs ultimate summing of the crystals (until it can't find anymore neighbours with lower energy)
# Gamma-energies of neighbours are added to the proton cluster. Threshold for energies of gammas is used (see definitions above )
# in order to disentangle between low energy proton signal and high-energy gamma when both are present.
# Found crystals which belong to one cluster are stored in the separate xb_list to calculate multiplicity in the cluster
*/
void adbp_4(std::list<xble> *asorted, std::list<xble> *psorted, std::list<xble> *gsorted) {
	std::list<xble>::iterator it1,it2,_it2,c_it,itp,itg;

	itp = (*psorted).begin();
	itg = (*gsorted).begin();
	int central = 0;
	std::list<xble> cluster; // to store crystals of the cluster

	//Go through the sorted list, from highest to lowest energy
	for (it1 = (*asorted).begin(); it1 != (*asorted).end(); it1++){
		cluster.clear(); //reset for the next cluster
		//To determine the type of the signal in the crystal:
		bool proton = false;
		bool gamma = false;
		//if ((*it1).is_proton_2pi()) 
		if ((*it1).is_proton_4pi()) 
		{	
			proton=true;
			central = (*it1)._number;
			cluster.push_back(*it1); //store this crystal in the list for the cluster 
			(*psorted).push_back(*it1);//store this crystal in the list for protons 
			itp = (*psorted).end();
			itp--;
		}
		else if ((*it1).is_gamma())
		{
			gamma=true;
			central = (*it1)._number;
			cluster.push_back(*it1);			
			(*gsorted).push_back(*it1);
			itg = (*gsorted).end();
			itg--;
		}
		else // nan proton energy and nan gamma. Only times are present
		{
			//it1 = (*asorted).erase(it1);
			std::cout <<" \n ERROR!!!" << (*it1)._number;
			continue;
		}
		int j = 0;//next crystal in the list
		it2=it1;
		it2++;
		std::cout << "\n~~~~~First check for the crystal " << central;
		//Inner loop over global list to find neghbours
		for (it2=it2 ; it2 != (*asorted).end(); it2++){
			j = (*it2)._number; //check this crystal from the list
			int k = 0;//neighbouring crystal from the list of neighbours
			//loop over the neighbours of that crystal
			for (int nb=0; nb<6; nb++){
				k = mycr[central-1].nb[nb]; // check this neighbour
				if (j==k){// Found neighbour of the central crystal  in globall list
					// Proton signal in both  central crystal and found neighbour
					if (proton && (*it2).is_proton_4pi())
					{
						cluster.push_back(*it2);
						//Overflow in back hemisphere
						if(isnan((*itp)._energy_p) || isnan((*it2)._energy_p)) // neighbour without proton energy
							std::cout<<"\nNAN energy in proton branch";
						//(*itp)._energy_p = (*it2)._energy_p;
						else
							(*itp)._energy_p += (*it2)._energy_p;
						(*itp)._mul++;
						//(*it2)._energy_p = 0;
						it2 = (*asorted).erase(it2);
						it2--;//go back to the element before the just erased one
					}
					//Gamma signal in both central crysta in found neighbour
					else if (gamma &&  (*it2).is_gamma())
					{
						cluster.push_back(*it2);
						(*itg)._energy_g += (*it2)._energy_g;
						(*itg)._mul++;
						//(*it2)._energy_g = 0;
						it2 = (*asorted).erase(it2);
						it2--;//go back to the element before the just erased one
					}
					//Proton in central crystal and gamma in the neighbour. Add gamma energy to the energy of proton
					else if(proton &&  (*it2).is_gamma())
					{
						cluster.push_back(*it2);

						//Overflow in back hemisphere
						//if(isnan((*itp)._energy_p))
						//	(*itp)._energy_p = (*it2)._energy_g;
						//else
						//(*itp)._energy_p += (*it2)._energy_g;

						(*itp)._energy_p += (*it2)._energy_g;
						(*itp)._mul++;
						//(*it2)._energy_g = 0;
						it2 = (*asorted).erase(it2);
						it2--;//go back to the element before the just erased one
					}
				} // end if
			}// end looping neighbours
		} // end looping global list to find first neghbous

		//Second check for next next neighbours
		std::cout << "\n#########FOUND CLUSTER with the central crystal " << central;
		for(c_it=cluster.begin(); c_it != cluster.end(); c_it++){
			std::cout << "\n" << (*c_it)._number;
		}
		if(cluster.size()<2) continue;

		std::cout << "\nSecond check...";

		c_it = cluster.begin();
		c_it++; // since the first element is the cluster's center (already looped)
		for(c_it=c_it; c_it != cluster.end(); c_it++){ //loop the cluster
			//Reset for the second check
			proton=false;
			gamma=false;

			float j_energy_p = (*c_it)._energy_p;
			float j_energy_g = (*c_it)._energy_g;
			j = (*c_it)._number;

			if( (*c_it).is_proton_4pi())
				proton=true;
			else if ( (*c_it).is_gamma())
				gamma=true;
			else continue;

			std::cout << "\n\tCheck neighbours of the crystal "<< j <<  " ";
			std::cout << j_energy_p << " MeV (proton_p) "; 
			std::cout << j_energy_g << " MeV (energy_g) in the cluster";

			for(int n=0; n<6; n++){ // loop the neighbours
				int k = mycr[j-1].nb[n]; // take this neighbour
				std::cout << "\n\t\tCheck neighbour crystal " << k << ".....";
				if(central==k) continue; //skip the central crystal in the cluster
				_it2=it1;
				_it2++;

				for (_it2=_it2; _it2 != (*asorted).end(); _it2++){// look fot this neighbours in the list
					//int k = mycr[j-1].nb[n]; // take this neighbour
					if((*_it2)._number==k){

						std::cout << "\n\t\t\tFOUND EXTRA NEIGHBOUR IN THE LIST !!! crystal ";
						std::cout << (*_it2)._number << " " <<(*_it2)._energy_p << " MeV (proton_e) "; 
						std::cout << (*_it2)._energy_g << " MeV (energy_g). Let's compare...";

						if (proton && (*_it2).is_proton_4pi()) 
							//&& (((*_it2)._energy_p<j_energy_p) || isnan((*_it2)._energy_p)))
						{
							cluster.push_back(*_it2); // cooment this line to account only to the sceond order neighbour

							if(isnan((*_it2)._energy_p))
								std::cout<<"NAN energy in proton branch!!!";

							else
								(*itp)._energy_p += (*_it2)._energy_p;

							(*itp)._mul++;
							_it2 = (*asorted).erase(_it2);
							_it2--;
							std::cout << "\n\t\t\t ******** added PROTON to PROTON ********";
						}
						else if (proton && (*_it2).is_gamma())
							//&&(*_it2)._energy_g < j_energy_p)
							//&& (*_it2)._energy_g>0 && !isinf((*_it2)._energy_g))
						{
							//if(isnan(j_energy_p) || ((*_it2)._energy_g < j_energy_p))
							//{	
							(*itp)._energy_p += (*_it2)._energy_g;
							(*itp)._mul++;
							cluster.push_back(*_it2);// cooment this line to account only to the sceond order neighbour
							_it2 = (*asorted).erase(_it2);
							_it2--;
							std::cout << "\n\t\t\t ******** added GAMMA to PROTON *********";
							//}
						}

						else if (gamma && (*_it2).is_gamma())
							//&& (*_it2)._energy_g < j_energy_g && !isinf((*_it2)._energy_g))
							//&& !isinf((*_it2)._energy_g))
						{
							(*itg)._energy_g += (*_it2)._energy_g;
							//(*itg)._mul++;
							(*itg)._mul = (*itg)._mul+1;
							std::cout << "\n DBG INCREASED GAMMA MUL...";
							cluster.push_back(*_it2);// cooment this line to account only to the sceond order neighbour
							_it2 = (*asorted).erase(_it2);
							_it2--;
							std::cout << "\n\t\t\t ******** added GAMMA to GAMMA ********";
						}

						else std::cout << "\n\t\t\t ******** NOT ADDED, GO TO NEXT CRYSTAL ...";

					}// endif
				}//loop global list to find neighbours
			}// looping neighbours
		}// looping cluster

		std::cout<< "\n\tFINAL CLUSTER: \n";
		print_all(cluster);

	}//Determined the proton clusters. When central crystal gamma energy not in overflow, re-declare it a gamma cluster...
	return;
}

//addback routine: protons (high energy events in XB) and gammas
//Read from the sorted list asorted of all crystals, create proton and gamma cluster lists
void adbp(std::list<xble> *asorted, std::list<xble> *psorted, std::list<xble> *gsorted) {
	std::list<xble>::iterator it1,it2,itp,itg;
	itp = (*psorted).begin();
	itg = (*gsorted).begin();
	int i = 0;
	//Go through the sorted list, from highest to lowest energy
	for (it1 = (*asorted).begin(); it1 != (*asorted).end(); it1++)
	{
		i = (*it1)._number;
		//std::cout << "First crystal:" << std::endl;
		//std::cout << i << "\t" << mycr[i-1].number <<std::endl;
		//high energy deposit, due to a charged particle (proton) hit:
		bool proton = false;
		bool gamma = false;
		if (isinf((*it1)._energy_g) && ((*it1)._energy_p>0))
		{	
			proton=true;
			//std::cout << "\tEnergies bef.:\t" << (*it1)._energy_g << "\t" << (*it1)._energy_p << std::endl;
			(*psorted).push_back(*it1);
			//pointer (iterator) to this element
			itp = (*psorted).end();
			itp--;
		}
		else if ((*it1)._energy_g>0 && !isinf((*it1)._energy_g))
		{
			gamma=true;
			(*gsorted).push_back(*it1);
			//pointer (iterator) to this element
			itg = (*gsorted).end();
			itg--;
		}
		else
		{
			//std::cout << "Neither proton nor gamma, maybe proton neighbour:\t" << (*it1)._energy_p<< "\t" << (*it1)._energy_g << std::endl;
			//it1 = (*asorted).erase(it1);
			//continue;
		}
		//check and add neighbours
		int j = 0;
		it2=it1;
		it2++;
		for (it2=it2 ; it2 != (*asorted).end(); it2++)
		{
			j = (*it2)._number;
			//std::cout << "\t" << j << "\t" << mycr[j-1].number <<std::endl;
			int k = 0;
			//loop over the neighbours of that crystal
			//std::cout << "\t\tNeighbours of the first crystal:" << std::endl;
			for (int nb=0; nb<6; nb++)
			{
				k = mycr[i-1].nb[nb];
				//if second crystal was a neighbour of the first
				if (j==k)
				{
					if (proton && (*it2)._energy_p>0)
					{//check if there was a proton in the first crystal, and if there is energy in the second
						//std::cout << "\t\tneighbour: " << k << std::endl;
						//std::cout << "\t\tProton Energies bef.:\t" << (*itp)._energy_p << "\t" << (*it2)._energy_p << std::endl;
						(*itp)._energy_p += (*it2)._energy_p;
						(*it2)._energy_p = 0;
						//std::cout << "\t\tProton Energies aft.:\t" << (*itp)._energy_p << "\t" << (*it2)._energy_p << std::endl;
					}
					else if (gamma && (*it2)._energy_g>0) //gamma hit in the first crystal, and energy in this one
					{
						//std::cout << "\t\tneighbour: " << k << std::endl;
						//std::cout << "\t\tGamma Energies bef.:\t" << (*itg)._energy_g << "\t" << (*it2)._energy_g << std::endl;
						(*itg)._energy_g += (*it2)._energy_g;
						(*it2)._energy_g = 0;
						//std::cout << "\t\tGamma Energies aft.:\t" << (*itg)._energy_g << "\t" << (*it2)._energy_g << std::endl;
					}
					it2 = (*asorted).erase(it2);
					it2--;//go back to the element before the just erased one
				}
			}
		}
	}//Determined the proton clusters. When central crystal gamma energy not in overflow, re-declare it a gamma cluster...
	return;
}

/* 
# Only first order (closest) neighbours of the central crystal
# are checked. Vary proton trigger in 2pi and 4pi
*/
void adbp_1t(std::list<xble> *asorted, std::list<xble> *psorted, std::list<xble> *gsorted, Float_t t0, TH2F *hist) 
{
	std::list<xble>::iterator it1,it2,itp,itg;
	itp = (*psorted).begin();
	itg = (*gsorted).begin();
	int i = 0;

	//Go through the sorted list, from highest to lowest energy
	for (it1 = (*asorted).begin(); it1 != (*asorted).end(); it1++)
	{
		i = (*it1)._number; //Number of the central crystal in this cluster

		//To determine the type of the signal in the crystal:
		bool proton = false;
		bool gamma = false;
		float g_time=0.;

		//if ((*it1).is_proton_2pi())
		if ((*it1).is_proton_4pi())
		{
			//std::cout << " \nPROTON!! in crystal      " << i;
			proton=true;
			(*psorted).push_back(*it1);
			itp = (*psorted).end();
			itp--;
		}
		else if((*it1).is_gamma())
			//else if((*it1).is_gamma(t0))
		{
			//std::cout << " \nGAMMA!! in crystal      " << i;
			gamma=true;
			(*gsorted).push_back(*it1);
			itg = (*gsorted).end();
			g_time=(*it1)._time;
			//hist->Fill((g_time-t0),((*it1)._energy_g));	
			itg--;
		}
		else{
			//std::cout << " ERROR Neither proton nor gamma, maybe proton neighbour:\t" << (*it1)._energy_p<< "\t" << (*it1)._energy_g << std::endl;
			//it1 = (*asorted).erase(it1);
			continue;
		}

		int j = 0;//number of the next crystal in the list
		it2=it1;
		it2++;
		for (it2=it2 ; it2 != (*asorted).end(); it2++)
		{
			j = (*it2)._number;

			int k = 0;
			//loop over the neighbours of that crystal
			//std::cout << "\t\tNeighbours of the first crystal:" << std::endl;
			for (int nb=0; nb<6; nb++)
			{
				k = mycr[i-1].nb[nb];
				//if second crystal is a neighbour of the first
				if (j==k)
				{
					if(proton && (*it2).is_proton_4pi())
						//if(proton && (*it2).is_proton_2pi())
					{
						if(!isnan((*it2)._energy_p))
							(*itp)._energy_p += (*it2)._energy_p;
						(*itp)._mul++;
						(*it2)._energy_p = 0;
					}

					else if(gamma && (*it2).is_gamma())
						//else if(gamma && (*it2).is_gamma(t0))
					{

						//hist->Fill((((*it2)._time)-t0),((*it2)._energy_g));	

						//if(((g_time-((*it2)._time))<(-50)) || ((g_time-((*it2)._time))>50))
						//	hist->Fill(g_time-((*it2)._time));	
						//else{
						(*itg)._energy_g += (*it2)._energy_g;
						(*itg)._mul++;
						(*it2)._energy_g = 0;
						//hist->Fill(g_time-((*it2)._time));
						//}
					}

					else if(proton && (*it2).is_gamma())
					{
						(*itp)._energy_p += (*it2)._energy_g;
						(*itp)._mul++;
						(*it2)._energy_g = 0;
					}
					it2 = (*asorted).erase(it2);
					it2--;//go back to the element before the just erased one
				} //endif
			}//end of neighbour loop
		}// end of second loop
	}//end of main loop
	return;
}	


/* 
# Only first order (closest) neighbours of the central crystal
# are checked. Vary proton trigger in 2pi and 4pi
*/
void adbp_1(std::list<xble> *asorted, std::list<xble> *psorted, std::list<xble> *gsorted) 
{
	std::list<xble>::iterator it1,it2,itp,itg;
	itp = (*psorted).begin();
	itg = (*gsorted).begin();
	int i = 0;
	//Go through the sorted list, from highest to lowest energy
	for (it1 = (*asorted).begin(); it1 != (*asorted).end(); it1++)
	{
		i = (*it1)._number; //Number of the central crystal in this cluster

		//To determine the type of the signal in the crystal:
		bool proton = false;
		bool gamma = false;

		if ((*it1).is_proton_2pi())
			//if ((*it1).is_proton_4pi())
		{
			//std::cout << " \nPROTON!! in crystal      " << i;
			proton=true;
			(*psorted).push_back(*it1);
			itp = (*psorted).end();
			itp--;
		}
		else if((*it1).is_gamma())
		{
			//std::cout << " \nGAMMA!! in crystal      " << i;
			gamma=true;
			(*gsorted).push_back(*it1);
			itg = (*gsorted).end();
			itg--;
		}
		else{
			//std::cout << " ERROR Neither proton nor gamma, maybe proton neighbour:\t" << (*it1)._energy_p<< "\t" << (*it1)._energy_g << std::endl;
			//it1 = (*asorted).erase(it1);
			continue;
		}

		int j = 0;//number of the next crystal in the list
		it2=it1;
		it2++;
		for (it2=it2 ; it2 != (*asorted).end(); it2++)
		{
			j = (*it2)._number;

			int k = 0;
			//loop over the neighbours of that crystal
			//std::cout << "\t\tNeighbours of the first crystal:" << std::endl;
			for (int nb=0; nb<6; nb++)
			{
				k = mycr[i-1].nb[nb];
				//if second crystal is a neighbour of the first
				if (j==k)
				{
					//if(proton && (*it2).is_proton_4pi())
					if(proton && (*it2).is_proton_2pi())
					{
						if(!isnan((*it2)._energy_p))
							(*itp)._energy_p += (*it2)._energy_p;
						(*itp)._mul++;
						(*it2)._energy_p = 0;
					}

					else if(gamma && (*it2).is_gamma())
					{
						(*itg)._energy_g += (*it2)._energy_g;
						(*itg)._mul++;
						(*it2)._energy_g = 0;
					}

					else if(proton && (*it2).is_gamma())
					{
						(*itp)._energy_p += (*it2)._energy_g;
						(*itp)._mul++;
						(*it2)._energy_g = 0;
					}
					it2 = (*asorted).erase(it2);
					it2--;//go back to the element before the just erased one
				} //endif
			}//end of neighbour loop
		}// end of second loop
	}//end of main loop
	return;
}	

double b2g(double beta)
{
	double gamma = 0;
	gamma = 1/sqrt(1-beta*beta);
	return gamma;
}

//doppler shift
//The beta should actually be the beta of the decaying fragment, i.e. 15O.
//That also should take into account the direction of the fragment relative to the incoming beam
void doppler(std::list<xble> *clus, double b, std::list<xble> *dopp)
	//float doppler(std::list<xble> *clus, double b, std::list<xble> *dopp)
{
	//float e_sum=0;

	double g = 0;
	g = b2g(b);
	//std::cout << "beta: " << b << "\tgamma: " << g <<std::endl;
	double elab, ecm, th_det = 0;
	int i = 0;
	//it1 for raw, it2 for doppler corrected list
	std::list<xble>::iterator it1,it2;
	for (it1 = (*clus).begin(); it1 != (*clus).end(); it1++)
	{
		(*dopp).push_back(*it1);
		it2 = (*dopp).end();
		it2--;
		i = (*it1)._number;				
		elab = (*it1)._energy_g;
		th_det = mycr[i-1].theta;
		//std::cout << i << "\t" << mycr[i-1].number <<std::endl;
		//std::cout << th_det << "\t" << mycr[i-1].theta <<std::endl;
		//std::cout << "\tGamma Energy Lab " << elab << "\t" << (*it1)._energy_g << std::endl;
		//do the doppler correction:
		ecm = elab * g * (1-b*cos(d2r(th_det)));
		(*it2)._energy_g = ecm;
		//e_sum += ecm;
		//std::cout << "\tGamma Energy C.M. " << ecm << "\t" << (*it2)._energy_g << std::endl;
	}
	//return e_sum;
	return;
}
//QFS 2proton events analysis



