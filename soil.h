#ifndef soil_h
#define soil_h
#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
class soil:public NDMaterial
{
public:
    soil();
    soil (int tag, double K, double G, double rho, double B, double J,double t);
    virtual  ~soil();
	
	NDMaterial *getCopy(void);
    NDMaterial *getCopy(const char *code) ;
	
	const char *getType (void) const;
	
	int  commitState(void);
    int  revertToLastCommit(void);
    int  revertToStart(void);
	
    double getRho(void);
    int setTrialStrain (const Vector &v);
	
    const Matrix &getTangent(void);
    const Matrix &getInitialTangent(void);
    const Vector &getStress(void);
    const Vector &getStrain(void);

    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker); 
    void Print(OPS_Stream &s, int flag =0);
private: 
	Vector strain_n;
	Vector strain_nplus1;
    static Vector sigma;         // Stress vector ... class-wide for returns
    static Matrix D;	          // Elastic constants
	double K;                      //Bulk modulus
	double G;                     //Shear Modulus
	double rho;                  //mass density
	double J;                //Voigt damping ratio
	double B;                //b = characteristic length of footing
	double T;  
};

 #endif