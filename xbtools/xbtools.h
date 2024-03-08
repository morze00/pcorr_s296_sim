#ifndef XBTOOLS_H
#define XBTOOLS_H

#include "xb_list_energy.hh"
#include "libs.hh"

double d2r(double deg);
double r2d(double rad);
void loadData();
void adbp_1(std::list<xble> *asorted, std::list<xble> *psorted, std::list<xble> *gsorted) ;
void adbp_1t(std::list<xble> *asorted, std::list<xble> *psorted, std::list<xble> *gsorted, Float_t t0, TH2F * hist) ;
void adbp_4(std::list<xble> *asorted, std::list<xble> *psorted, std::list<xble> *gsorted) ;
void adbp_all(std::list<xble> *asorted, std::list<xble> *psorted, std::list<xble> *gsorted) ;
void adbp(std::list<xble> *asorted, std::list<xble> *psorted, std::list<xble> *gsorted) ;
void doppler(std::list<xble> *clus, double b, std::list<xble> *dopp);
//float doppler(std::list<xble> *clus, double b, std::list<xble> *dopp);

struct vect {
	double x;
	double y;
	double z;
};
struct vectsp {
	double r;
	double t;
	double p;
};
struct matr {
	double xx, yx, zx;
	double xy, yy, zy;
	double xz, yz, zz;
};

struct crys {
	int number;
	char shape;
	double theta;
	double phi;
	int nb[6];
	double tx;
	double px;
	double ty;
	double py;
	double tz;
	double pz;
	struct vect v_tra;
	struct matr m_rot;
};
extern struct crys mycr[162];	//declare the XB crystals...


#endif
