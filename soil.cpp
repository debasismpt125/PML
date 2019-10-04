
#include <elementAPI.h>
#include "soil.h"
#include <Vector.h>
#include <Channel.h>
#include <OPS_Globals.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#ifdef _USRDLL
#define OPS_Export extern "C" _declspec(dllexport)
#elif _MACOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
#define OPS_Export extern "C"
#endif

Vector soil::sigma(3);
Matrix soil::D(3,3);
//reading data from command and store these in a array
OPS_Export void *
OPS_soil()
{ NDMaterial *theMaterial = 0;
  int numArgs = OPS_GetNumRemainingInputArgs();
   if (numArgs < 7) {
    opserr << "Want: nDMaterial soil $tag $K $G $rho $b $jita $deltaT" << endln;
    return 0;	
  }
  int iData[1];
  double dData[6];
  int numData = 1;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: nDMaterial soil\n";
    return 0;
  }
  numData=6;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid data: nDMaterial soil : " << iData[0] <<"\n";
    return 0;
  }  
  
 theMaterial = new soil(iData[0], dData[0], dData[1], dData[2],dData[3],dData[4],dData[5]);
 return theMaterial;
}

//null constructor
 soil::soil()
:NDMaterial(0,0),strain_n(3), strain_nplus1(3)
 { strain_n.Zero();
  strain_nplus1.Zero();
 }
 
 //full constructor
soil::soil(int tag, double k, double g, double r, double b, double j, double t)
:NDMaterial(tag,0),strain_n(3),strain_nplus1(3)
{
    //initialize variables
  strain_n.Zero();
  strain_nplus1.Zero();

  K=k;
  G=g;
  rho=r;
  B=b;
  J=j;
  T=t;
}

 
 soil::~soil()
{
	//destructor 
}


NDMaterial *soil::getCopy(void)  {
  
 soil *theCopy = new soil(*this);
 theCopy->strain_nplus1 = strain_nplus1;
 theCopy->strain_n = strain_n;
  return theCopy;
};


NDMaterial * soil::getCopy(const char *code)  {
  if (strcmp(code,this->getType()) == 0) {
  soil *copy = new soil(*this);
  copy->strain_nplus1 = strain_nplus1;
  copy->strain_n = strain_n;
  return copy;
  }
  return 0;}

const char*
soil::getType (void) const
{
  return "PlaneStrain";
}
double
soil::getRho() 
{ 
  return rho ;
}

int
soil::commitState(void)
{  
   strain_n=strain_nplus1;
   
   return 0;
}

int
soil::revertToLastCommit (void)
{ 
  strain_nplus1=strain_n;
  return 0;
}

int
soil::revertToStart (void)
{
strain_n.Zero();
strain_nplus1.Zero();
return 0;
}


int
soil::setTrialStrain (const Vector &v)
{   
    strain_nplus1 = v;
	
    return 0;
}

const Matrix&
soil::getTangent (void)
{
  D(0,0) = D(1,1) = K+4.0*G/3.0;
  D(0,1) = D(1,0) = K-2.0*G/3.0;
  D(2,2) = G;
  
  return D;
}

const Matrix&
soil::getInitialTangent (void)
{
  D(0,0) = D(1,1) = K+4.0*G/3.0;
  D(0,1) = D(1,0) = K-2.0*G/3.0;
  D(2,2) = G;
  
  return D;
}
const Vector&
soil::getStrain (void)
{ 
  return strain_nplus1;
}
const Vector&
soil::getStress (void)
{ 
   
   double l1 = K+4.0*G/3.0;
   double l2 = K-2.0*G/3.0;
   double R1 = sqrt(G/rho);
   double S1 = 1.0+2.0*J*B/(R1*T);
   double S2 = 2.0*J*B/(R1*T) ;
   double k1 = l1*strain_nplus1(0)+l2*strain_nplus1(1);
   double k2 = l1*strain_n(0)+l2*strain_n(1);
   double k3 = l2*strain_nplus1(0)+l1*strain_nplus1(1);
   double k4 = l2*strain_n(0)+l1*strain_n(1);
   sigma(0) = S1*(k1)-S2*(k2);
   sigma(1) = S1*(k3)-S2*(k4);
   sigma(2) = S1*G*strain_nplus1(2)-S2*G*strain_n(2);
   return sigma;
}

int 
soil::sendSelf(int commitTag, Channel &theChannel)
{
  
  static Vector data(10);
  data(0) = this->getTag();
  data(1) = K;
  data(2) = G;
  data(3) = rho;
  data(4) = J;
  data(5) = B;
  data(6) = T;
  data(7) =strain_nplus1(0);
  data(8) =strain_nplus1(1);
  data(9) =strain_nplus1(2);
 
  int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "soil::sendSelf -- could not send Vector\n";
    return res;
  }
  
  return res;
}

int 
soil::recvSelf(int commitTag, Channel &theChannel, 
					FEM_ObjectBroker &theBroker)
{
  static Vector data(10);
  
  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "soil::sendSelf -- could not send Vector\n";
    return res;
  }

  this->setTag((int)data(0));
  K   = data(1);
  G   = data(2);
  rho = data(3);
  J   = data(4);
  B   = data(5);
  T   = data(6);
  strain_nplus1(0)=data(7);
  strain_nplus1(1)=data(8);
  strain_nplus1(2)=data(9);

  return res;
}
 
 void
soil::Print(OPS_Stream &s, int flag)
{

  s << "  K:  " << K << endln;
  s << "  G:  " << G << endln;
  s << "  rho:" << rho<< endln;
  s << "  J:  " << J << endln;
  s <<  "  B: " << B <<endln;
  s << "  T:  " << T <<endln;

}