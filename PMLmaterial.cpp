//.cpp file of PMLmaterial
#include <elementAPI.h>
#include "PMLmaterial.h"
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
Vector PMLmaterial::sigma(3);
Matrix PMLmaterial::D(3,3);
OPS_Export void *
OPS_PMLmaterial()
{ NDMaterial *theMaterial = 0;
  int numArgs = OPS_GetNumRemainingInputArgs();
  
  if (numArgs < 6) {
    opserr << "Want: nDMaterial PMLmaterial $tag $K $G $rho $jita $b" << endln;
    return 0;	
  }
  int iData[1];
  double dData[5];
  int numData = 1;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: nDMaterial PMLmaterial\n";
    return 0;
  }
  numData=5;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid data: nDMaterial PMLmaterial : " << iData[0] <<"\n";
    return 0;
  }  
  
 theMaterial = new PMLmaterial(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4]);
 return theMaterial;
}
PMLmaterial::PMLmaterial(int tag, double k, double g, double r, double j, double b)
:NDMaterial(tag,0), K(k), G(g), rho(r), jita(j), base(b)
{
    //initialize variables
	this->revertToStart();
	
}
PMLmaterial::PMLmaterial()
:NDMaterial(0,0),
 K(0.0), G(0.0), rho(0.0), jita(0.0), base(0.0)
 {
 }
 PMLmaterial::~PMLmaterial()
{
	
}
const char*
PMLmaterial::getType (void) const
{
  return "PlaneStrain";
}
NDMaterial * PMLmaterial::getCopy(void)  {
  
  PMLmaterial * copy = new PMLmaterial(*this);
  return copy;
  
};
NDMaterial * PMLmaterial::getCopy(const char *code)  {
  if (strcmp(code,this->getType()) == 0) {
    PMLmaterial * copy = new PMLmaterial(*this);
    return copy;
  }
  else
    return 0;
};

double
PMLmaterial::getRho() 
{ 
  return rho ;
}
int
PMLmaterial::setTrialStrain (const Vector &strain)
{
  strain_nplus1 = strain;
  return 0;
}
const Matrix&
PMLmaterial::getTangent (void)
{
  D(0,0) = D(1,1) = K+4*G/3;
  D(0,1) = D(1,0) = K-2*G/3;
  D(2,2) = G;
  
  return D;
}
const Matrix&
PMLmaterial::getInitialTangent (void)
{
  D(0,0) = D(1,1) = K+4*G/3;
  D(0,1) = D(1,0) = K-2*G/3;
  D(2,2) = G;
  
  return D;
}
const Vector&
PMLmaterial::getStress (void)
{
  double eps0= strain_nplus1(0);
  double eps1= strain_nplus1(1);
  double eps2= strain_nplus1(2);
  double peps0= strain_n(0);
  double peps1= strain_n(1);
  double peps2= strain_n(2);
  double cs=sqrt(G/rho);
  double c1=1+2*jita*base/(cs*ops_Dt);
  double c2=-2*jita*base/(cs*ops_Dt);
  double l1=K+4*G/3;
  double l2=K-2*G/3;
  sigma(0) = c1*l1*eps0+c1*l2*eps1+c2*l1*peps0+c2*l2*peps1;
  sigma(1) = c1*l2*eps0+c2*l1*eps1+c2*l2*peps0+c2*l1*peps1;
  sigma(2) = c1*G*eps2+c2*G*peps2;
  return sigma;
}

const Vector&
PMLmaterial::getStrain (void)
{
  return strain_nplus1;
}
int
PMLmaterial::commitState(void)
{
  strain_n=strain_nplus1;
  return 0;
}
int
PMLmaterial::revertToLastCommit (void)
{
  strain_nplus1=strain_n;
  return 0;
}
int
PMLmaterial::revertToStart (void)
{
  strain_nplus1.Zero();
  strain_n.Zero();
  return 0;
}
int 
PMLmaterial::sendSelf(int commitTag, Channel &theChannel)
{
  
  static Vector data(12);
  
  data(0) = this->getTag();
  data(1) = K;
  data(2) = G;
  data(3) = rho;
  data(4) = jita;
  data(5) = base;
  data(6) = strain_n(0);
  data(7) = strain_n(1);
  data(8) = strain_n(2);
  data(9) = strain_nplus1(0);
  data(10)= strain_nplus1(1);
  data(11)= strain_nplus1(2);
  int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "PMLmaterial::sendSelf -- could not send Vector\n";
    return res;
  }
  
  return res;
}

int 
PMLmaterial::recvSelf(int commitTag, Channel &theChannel, 
					FEM_ObjectBroker &theBroker)
{
  static Vector data(12);
  
  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "PMLmaterial::sendSelf -- could not send Vector\n";
    return res;
  }

  this->setTag((int)data(0));
   K = data(1);
  G = data(2);
   rho= data(3);
   jita = data(4);
   base=data(5);
  strain_n(0)=data(6);
  strain_n(1)=data(7);
  strain_n(2)=data(8);
  strain_nplus1(0)=data(9);
  strain_nplus1(1)=data(10);
  strain_nplus1(2)=data(11);
return res;
}
 
 void
PMLmaterial::Print(OPS_Stream &s, int flag)
{

  s << "  E: " << K << endln;
  s << "  G: " << G<< endln;
  s << "strain(0)"<<strain_nplus1(0);
  s << "strain(1)"<<strain_nplus1(1);
  s << "strain(2)"<<strain_nplus1(2);
 
}

