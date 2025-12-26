// myDEM01.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "CParticle.h"

using namespace std;

double xMin(0.0), xMax(0.1), yMin(0.0), yMax(0.1);
int nX(10), nY(10);

void calContactForce(CParticle particle1, CParticle particle2, C3dValue&contactForce1, C3dValue&contactForce2)
{
	C3dValue unitCPVect(particle2.m_position.m_x - particle1.m_position.m_x, particle2.m_position.m_y - particle1.m_position.m_y, particle2.m_position.m_z - particle1.m_position.m_z);
	unitCPVect.normalize();


	// 计算法向重叠量
	double overlap = (fabs(particle1.m_diameter * 0.5 + particle2.m_diameter * 0.5) -
		sqrt(pow(particle1.m_position.m_x - particle2.m_position.m_x, 2) +
			pow(particle1.m_position.m_y - particle2.m_position.m_y, 2) +
			pow(particle1.m_position.m_z - particle2.m_position.m_z, 2)));

	if (overlap > 0)
	{
		// 计算接触力（简化为线性弹簧模型）
		C3dValue unitCPVect = particle1.m_velocity - particle1.m_velocity;
		unitCPVect.normalize();
		double k_hertz = 4.0 / 3.0 * particle1.m_shearModulus * sqrt(particle1.m_diameter * 0.5);
		//double sqrt_meff_k_hertz = sqrt((mass * mass) / (mass + mass) * k_hertz);
		//F_n = -k_hertz * pow(overlap, 1.5) - sqrt_meff_k_hertz * particle.m_poissonRatio * pow(overlap, 0.25) * relVel.m_z;
		C3dValue F_n = C3dValue(0.0, k_hertz * overlap, 0.0);

		//C3dValue contactForce = C3dValue(0, 0, 0); // 初始化接触力向量

		//// 计算法向量
		//C3dValue normalVector = C3dValue(
		//	particle2.m_position.m_x - particle1.m_position.m_x,
		//	particle2.m_position.m_y - particle1.m_position.m_y,
		//	particle2.m_position.m_z - particle1.m_position.m_z);
		//double length = sqrt(pow(normalVector.m_x, 2) + pow(normalVector.m_y, 2) + pow(normalVector.m_z, 2));
		//normalVector = C3dValue(normalVector.m_x / length, normalVector.m_y / length, normalVector.m_z / length);

		//// 计算接触力
		//contactForce = normalVector * (k_contact * overlap);

		// 在这里可以将contactForce应用到两个颗粒上
		contactForce1 = contactForce1 + F_n; // 颗粒1受到的接触力
		contactForce2 = contactForce2 - F_n; // 颗粒2受到的接触力
	}

	return;
}

void writeParticlesToVTKWithAttributes(const std::vector<CParticle>& particles, const std::string& filename) {
	std::ofstream vtkFile(filename);
	if (!vtkFile.is_open()) {
		std::cerr << "Error: Could not open file " << filename << std::endl;
		return;
	}

	int numParticles = particles.size();

	// VTK头部
	vtkFile << "# vtk DataFile Version 3.0" << std::endl;
	vtkFile << "Particle data with attributes" << std::endl;
	vtkFile << "ASCII" << std::endl;
	vtkFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

	// 点坐标
	vtkFile << "POINTS " << numParticles << " float" << std::endl;
	for (const auto& p : particles) {
		vtkFile << p.m_position.m_x << " " << p.m_position.m_y << " " << p.m_position.m_z << std::endl;
	}

	// 单元信息
	vtkFile << "CELLS " << numParticles << " " << (2 * numParticles) << std::endl;
	for (int i = 0; i < numParticles; i++) {
		vtkFile << "1 " << i << std::endl;
	}

	// 单元类型
	vtkFile << "CELL_TYPES " << numParticles << std::endl;
	for (int i = 0; i < numParticles; i++) {
		vtkFile << "1" << std::endl;
	}

	// 点数据
	vtkFile << "POINT_DATA " << numParticles << std::endl;

	// 标量数据示例：颗粒ID
	vtkFile << "SCALARS ParticleID int" << std::endl;
	vtkFile << "LOOKUP_TABLE default" << std::endl;
	for (int i = 0; i < numParticles; i++) {
		vtkFile << i << std::endl;
	}

	// 向量数据示例：假设颗粒有速度属性 vx, vy, vz
	//vtkFile << "VECTORS Velocity float" << std::endl;
	//vtkFile << "LOOKUP_TABLE default" << std::endl;
	//for (const auto& p : particles) {
	//	// 假设CParticle有vx, vy, vz成员
	//	// vtkFile << p.vx << " " << p.vy << " " << p.vz << std::endl;
	//	// 如果没有速度数据，可以注释掉这部分
	//}

	// 向量数据示例：假设颗粒有接触力属性 contactForceX, contactForceY, contactForceZ
	vtkFile << "VECTORS ContactForce float" << std::endl;
	//vtkFile << "LOOKUP_TABLE default" << std::endl;
	for (const auto& p : particles) {
		// 假设CParticle有vx, vy, vz成员
		 vtkFile << p.m_contactForce.m_x << " " << p.m_contactForce.m_y << " " << p.m_contactForce.m_z << std::endl;
	}

	// 另一个标量数据示例：半径（假设CParticle有radius成员）
	//vtkFile << "SCALARS Radius float" << std::endl;
	//vtkFile << "LOOKUP_TABLE default" << std::endl;
	//for (const auto& p : particles) {
	//	// vtkFile << p.radius << std::endl;
	//	// 如果没有半径数据，可以注释掉这部分
	//}

	vtkFile.close();
	//std::cout << "Successfully wrote particle data to " << filename << std::endl;
}

void initContactLists(vector<CParticle> particles, vector<vector<int>> &vecHeader, vector<int> &vecNext)
{
	// 将vecHeader中所有元素设置为-1
	for (auto& row : vecHeader) {
		fill(row.begin(), row.end(), -1);
	}

	// 将vecNext中所有元素设置为-1
	fill(vecNext.begin(), vecNext.end(), -1);

	// 1. 构建接触列表
	for (int i = 0; i < particles.size(); i++)
	{
		// 整数化颗粒位置以确定其所在的网格单元
		int ix = static_cast<int>((particles[i].m_position.m_x - xMin) / (xMax - xMin) * nX);
		int iy = static_cast<int>((particles[i].m_position.m_y - yMin) / (yMax - yMin) * nY);
		if(vecHeader[ix][iy] == -1)
		{
			vecHeader[ix][iy] = i;
		}
		else
		{
			int currentIndex = vecHeader[ix][iy];
			while(vecNext[currentIndex] != -1)
			{
				currentIndex = vecNext[currentIndex];
				//cout << "循环1中" << endl;
			}
			vecNext[currentIndex] = i;
		}
	}
}

void contactDetect(vector<CParticle> particles, vector<vector<int>> vecHeader, vector<vector<bool>> vecIsFirst, vector<int> vecNext, vector<C3dValue>& vecContactForce)
{
	fill(vecContactForce.begin(), vecContactForce.end(), C3dValue(0.0, 0.0, 0.0));
	//接触检测
	for (int iParticle = particles.size() - 1; iParticle >= 0; iParticle--)
	{
		//cout << "检查第" << iParticle << "个颗粒" << endl;
		// 整数化颗粒位置以确定其所在的网格单元
		int ix = static_cast<int>((particles[iParticle].m_position.m_x - xMin) / (xMax - xMin) * nX);
		int iy = static_cast<int>((particles[iParticle].m_position.m_y - yMin) / (yMax - yMin) * nY);

		// 检查当前单元中的颗粒
		int currentIndex = vecHeader[ix][iy];
		if (vecIsFirst[ix][iy])
		{
			while (vecNext[currentIndex] != -1)
			{
				if (vecNext[currentIndex] != iParticle && particles[iParticle].isContact(particles[currentIndex]))
				{
					// 发生接触，计算接触力
					C3dValue force1, force2;
					calContactForce(particles[iParticle], particles[currentIndex], force1, force2);
					vecContactForce[iParticle] = vecContactForce[iParticle] + force1;
					vecContactForce[currentIndex] = vecContactForce[currentIndex] + force2;
				}
				currentIndex = vecNext[currentIndex];
			}
		}

		// 检查相邻单元中的颗粒
		vector<int> vecNeighborX{ ix - 1, ix - 1, ix, ix + 1 };
		vector<int> vecNeighborY{ iy, iy - 1, iy - 1, iy - 1 };

		for (int i = 0; i < vecNeighborX.size(); i++)
		{
			int neighborX = vecNeighborX[i];
			int neighborY = vecNeighborY[i];

			// 检查邻居单元
			if (neighborX >= 0 && neighborX < nX && neighborY >= 0 && neighborY < nY)
			{
				currentIndex = vecHeader[neighborX][neighborY];
				while (currentIndex != -1)
				{
					if (currentIndex != iParticle && particles[iParticle].isContact(particles[currentIndex]))
					{
						// 发生接触，计算接触力
						C3dValue force1, force2;
						calContactForce(particles[iParticle], particles[currentIndex], force1, force2);
						vecContactForce[iParticle] = vecContactForce[iParticle] + force1;
						vecContactForce[currentIndex] = vecContactForce[currentIndex] + force2;
					}
					currentIndex = vecNext[currentIndex];
				}
			}
		}

		// 检查是否与计算域碰撞，如果颗粒与地面接触，计算接触力
		//cout << "posY" << particles[iParticle].m_position.m_y << "\tdiameter:" << particles[iParticle].m_diameter * 0.5;
		if (particles[iParticle].m_position.m_y < yMin + particles[iParticle].m_diameter * 0.5)
		{
			C3dValue F_n;
			double mass = particles[iParticle].m_density * (4.0 / 3.0) * 3.1415926 * pow(particles[iParticle].m_diameter * 0.5, 3);
			double overlap = fabs(particles[iParticle].m_diameter * 0.5 - particles[iParticle].m_position.m_y);
			C3dValue relVel = particles[iParticle].m_velocity;
			C3dValue unitCPVect = particles[iParticle].m_velocity;
			unitCPVect.normalize();
			double k_hertz = 4.0 / 3.0 * particles[iParticle].m_shearModulus * sqrt(particles[iParticle].m_diameter * 0.5);
			//double sqrt_meff_k_hertz = sqrt((mass * mass) / (mass + mass) * k_hertz);
			//F_n = -k_hertz * pow(overlap, 1.5) - sqrt_meff_k_hertz * particle.m_poissonRatio * pow(overlap, 0.25) * relVel.m_z;
			F_n = C3dValue(0.0, k_hertz * overlap, 0.0);

			vecContactForce[iParticle] = vecContactForce[iParticle] + F_n;
		}

		if (particles[iParticle].m_position.m_y > yMax - particles[iParticle].m_diameter * 0.5)
		{
			C3dValue F_n;
			double mass = particles[iParticle].m_density * (4.0 / 3.0) * 3.1415926 * pow(particles[iParticle].m_diameter * 0.5, 3);
			double overlap = fabs(particles[iParticle].m_diameter * 0.5 - particles[iParticle].m_position.m_y);
			C3dValue relVel = particles[iParticle].m_velocity;
			double k_hertz = 4.0 / 3.0 * particles[iParticle].m_shearModulus * sqrt(particles[iParticle].m_diameter * 0.5);
			double sqrt_meff_k_hertz = sqrt((mass * mass) / (mass + mass) * k_hertz);
			//F_n = -k_hertz * pow(overlap, 1.5) - sqrt_meff_k_hertz * particle.m_poissonRatio * pow(overlap, 0.25) * relVel.m_z;
			F_n = C3dValue(0, k_hertz * overlap, 0);

			vecContactForce[iParticle] = vecContactForce[iParticle] - F_n;
		}
	}
}


int main()
{
	vector<CParticle> particles;
	particles.push_back(CParticle(0, 0.01, 0.005, 0.015, 0, 1e7, 0.3, 2500));
	particles.push_back(CParticle(1, 0.01, 0.005, 0.035, 0, 1e8, 0.3, 2500));
	particles.push_back(CParticle(2, 0.01, 0.01, 0.05, 0, 1e8, 0.3, 2500));
	//particles.push_back(CParticle(3, 0.01, 0.05, 0.015, 0, 1e8, 0.3, 2500));
	//particles.push_back(CParticle(4, 0.01, 0.065, 0.015, 0, 1e8, 0.3, 2500));
	//particles.push_back(CParticle(5, 0.01, 0.08, 0.015, 0, 1e8, 0.3, 2500));
	//particles.push_back(CParticle(6, 0.01, 0.005, 0.03, 0, 1e8, 0.3, 2500));
	//particles.push_back(CParticle(7, 0.01, 0.005, 0.045, 0, 1e8, 0.3, 2500));
	//particles.push_back(CParticle(8, 0.01, 0.005, 0.065, 0, 1e8, 0.3, 2500));
	//particles.push_back(CParticle(9, 0.01, 0.051, 0.061, 0, 1e8, 0.3, 2500));

	vector<vector<int>> vecHeader(nX, vector<int>(nY, -1));
	vector<vector<bool>> vecIsFirst(nX, vector<bool>(nY, true));
	vector<int> vecNext(particles.size(), -1);
	vector<C3dValue> vecContactForce(particles.size(), C3dValue(0.0, 0.0, 0.0));

	// 1. 打开CSV文件用于写入颗粒数据
	std::ofstream outFile("particle_data.csv");

	// 2. 写入CSV文件头
	outFile << "time,x,y,z,velX,velY,velZ" << std::endl;

	// 3. 定义颗粒初始位置
	//CParticle particle(0, 0.01, 0, 0, 0.01, 1e8, 0.3, 2500);

	// 4. 定义重力
	C3dValue gravity(0, -9.81, 0);

	

	// 5. 定义时间步长和总模拟时间
	double timestep = 1e-5;
	double simTime = 0.0;
	double totalTime = 0.5;
	int count = 0;

	// 6. 开始迭代循环
	while(simTime < totalTime)
	{
		
		C3dValue totalForce(0, 0, 0);

		// 获取流体的速度来计算曳力
		C3dValue dragForce(0, 0, 1);
		double F_n = 0.0;

		initContactLists(particles, vecHeader, vecNext);

		contactDetect(particles, vecHeader, vecIsFirst, vecNext, vecContactForce);
		
		// 计算颗粒接触力
		//for (int i = 0; i < particles.size(); i++)
		//{
		//	C3dValue newVel, newPosition;
		//	double mass = particles[i].m_density * (4.0 / 3.0) * 3.1415926 * pow(particles[i].m_diameter * 0.5, 3);

		//	//如果颗粒与地面接触，计算接触力
		//	if (particles[i].m_position.m_z < particles[i].m_diameter * 0.5)
		//	{
		//		double overlap = fabs(particles[i].m_diameter * 0.5 - particles[i].m_position.m_z);
		//		C3dValue relVel = particles[i].m_velocity;
		//		double k_hertz = 4.0 / 3.0 * particles[i].m_shearModulus * sqrt(particles[i].m_diameter * 0.5);
		//		double sqrt_meff_k_hertz = sqrt((mass * mass) / (mass + mass) * k_hertz);
		//		//F_n = -k_hertz * pow(overlap, 1.5) - sqrt_meff_k_hertz * particle.m_poissonRatio * pow(overlap, 0.25) * relVel.m_z;
		//		F_n = -k_hertz * overlap;

		//		totalForce = gravity + dragForce + C3dValue(0, 0, -1) * F_n;

		//		newVel = particles[i].m_velocity + totalForce * timestep;
		//		newPosition = particles[i].m_position + newVel * timestep;
		//	}
		//	else
		//	{
		//		totalForce = gravity;

		//		newVel = particles[i].m_velocity + totalForce * timestep;
		//		newPosition = particles[i].m_position + newVel * timestep;
		//		//newPosition = particle.m_position + particle.m_velocity * timestep;
		//		//newPosition = particle.m_position + (particle.m_velocity + newVel) * 0.5 * timestep;
		//		//newPosition = particle.m_position + newVel * timestep + 0.5 * totalForce / mass * timestep * step;
		//	}

		//	// 更新颗粒位置和速度
		//	particles[i].m_position = newPosition;
		//	particles[i].m_velocity = newVel;
		//}

		// 更新体积力
		for (int i = 0; i < particles.size(); i++)
		{
			double mass = particles[i].m_density * (4.0 / 3.0) * 3.1415926 * pow(particles[i].m_diameter * 0.5, 3);
			C3dValue newVel, newPosition;
			C3dValue totalForce = gravity * mass + vecContactForce[i];
			newVel = particles[i].m_velocity + totalForce * (timestep / mass);
			newPosition = particles[i].m_position + newVel * timestep;

			particles[i].m_position = newPosition;
			particles[i].m_velocity = newVel;

	/*		cout << "time:" << simTime << "\tId" << i << "\tPosition:(" << particles[i].m_position.m_x << "," << particles[i].m_position.m_y << "," << particles[i].m_position.m_z << ")\t"
				<< "Velocity:(" << particles[i].m_velocity.m_x << "," << particles[i].m_velocity.m_y << "," << particles[i].m_velocity.m_z << ")\t"
				<< "ContactForce:(" << vecContactForce[i].m_x << "," << vecContactForce[i].m_y << "," << vecContactForce[i].m_z << ")" << endl;*/
		}


		// 更新时间
		simTime += timestep;

		// 输出数据到VTK
		count++;
		int interval = 100;
		if(count % interval == 0)
		{
			//cout << "time=" << simTime << "\t" << "颗粒数量为" << particles.size() << endl;
			writeParticlesToVTKWithAttributes(particles, "results_" + std::to_string(count / interval) + ".vtk");
		}

	}

	// 7. 关闭文件
	outFile.close();

	std::cout << "Simulation complete!" << std::endl;
}
