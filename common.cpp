#include "common.h"

bool IsConvex(const NormalType& n1, const PointType& p1, const PointType& p2)
{
	return n1.dot(p1 - p2) < LBND_FLT;
}

float GetGeodesicDistance(const NormalType& n1, const NormalType& n2, const PointType& p1, const PointType& p2)
{
	auto P1P2 = p1 - p2;
	auto nn1 = n1.norm(), nn2 = n2.norm();
	auto cosAlpha = n1.dot(P1P2) / P1P2.norm(), cosBeta = -n2.dot(P1P2) / P1P2.norm();
	auto alpha = std::acos(cosAlpha), beta = std::acos(cosBeta);
	auto sinAlphaBeta = std::sin(alpha + beta);
	auto x = P1P2.norm() * (cosAlpha + cosBeta) / sinAlphaBeta;
	auto delta = std::abs(x - P1P2.norm());
	return x;
}

float GetAngleDistance(const NormalType& n1, const NormalType& n2, const PointType& p1, const PointType& p2)
{
	float eta = IsConvex(n1, p1, p2) ? 0.2 : 1.0;
	return eta * (1.0 - n1.dot(n2));
}

void PaletteHSV::HSVtoRGB(int& R, int& G, int& B, float H, float S, float V)
{
	float s = S / 100;
	float v = V / 100;
	float C = s * v;
	float X = C * (1 - abs(fmod(H / 60.0, 2) - 1));
	float m = v - C;
	float r, g, b;
	if (H >= 0 && H < 60) {
		r = C, g = X, b = 0;
	}
	else if (H >= 60 && H < 120) {
		r = X, g = C, b = 0;
	}
	else if (H >= 120 && H < 180) {
		r = 0, g = C, b = X;
	}
	else if (H >= 180 && H < 240) {
		r = 0, g = X, b = C;
	}
	else if (H >= 240 && H < 300) {
		r = X, g = 0, b = C;
	}
	else {
		r = C, g = 0, b = X;
	}
	R = (r + m) * 255;
	G = (g + m) * 255;
	B = (b + m) * 255;
}

PaletteHSV::PaletteHSV(uint cls)
{
	float steps = 360.0 / cls;
	for (int k = 0; k < cls; ++k)
	{
		int r, g, b;
		this->HSVtoRGB(r, g, b, k * steps, 50, 75);
		this->colors.push_back({ r,g,b });
	}
}