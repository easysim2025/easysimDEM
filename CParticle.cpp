#include "CParticle.h"
#include "CParticle.h"


CParticle::CParticle()
{
	m_id = 0;
}

CParticle::CParticle(int id,
	double diameter,
	double posX,
	double posY,
	double posZ,
	double shearModulus,
	double poissionRatio,
	double density)
{
	m_id = id;
	m_diameter = diameter;
	m_position = C3dValue(posX, posY, posZ);
	m_shearModulus = shearModulus;
	m_poissonRatio = poissionRatio;
	m_density = density;
}