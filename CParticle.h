#pragma once
#include <math.h>

class C3dValue
{
public:
    C3dValue() : m_x(0.0), m_y(0.0), m_z(0.0) {}

    C3dValue(const C3dValue& val) { *this = val; }

    C3dValue(double x, double y, double z) : m_x(x), m_y(y), m_z(z) {}

    C3dValue& operator=(const C3dValue& rh)
    {
        if (&rh != this)
        {
            m_x = rh.m_x;
            m_y = rh.m_y;
            m_z = rh.m_z;
        }

        return (*this);
    }

    C3dValue operator*(double scalar) const {
        return C3dValue(m_x * scalar, m_y * scalar, m_z * scalar);
    }

    C3dValue operator+(const C3dValue& rhs) const {
        return C3dValue(m_x + rhs.m_x, m_y + rhs.m_y, m_z + rhs.m_z);
    }

    C3dValue operator-(const C3dValue& rhs) const {
        return C3dValue(m_x - rhs.m_x, m_y - rhs.m_y, m_z - rhs.m_z);
    }

    void normalize() {
		double length = pow(m_x * m_x + m_y * m_y + m_z * m_z, 0.5);
		if (length > 0.0) {
			m_x /= length;
			m_y /= length;
			m_z /= length;
		}
	}



public:
    double m_x, m_y, m_z;
};

class CParticle
{
public:

	CParticle();

    CParticle(int id, double diameter, double posX, double posY, double posZ, double shearModulus, double poissionRatio, double density);

    bool isContact(CParticle particle)
    {
        double distanceSquared = (m_position.m_x - particle.m_position.m_x) * (m_position.m_x - particle.m_position.m_x) +
            (m_position.m_y - particle.m_position.m_y) * (m_position.m_y - particle.m_position.m_y) +
            (m_position.m_z - particle.m_position.m_z) * (m_position.m_z - particle.m_position.m_z);
        double radiusSum = (particle.m_diameter * 0.5) * 2.0;
        return distanceSquared <= (radiusSum * radiusSum);
    }


public:

	int    m_id;
	C3dValue m_position;
	C3dValue m_velocity;
    double m_diameter;
    double m_shearModulus;
    double m_poissonRatio;
    double m_density;
    C3dValue m_contactForce;
	
};