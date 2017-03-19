#pragma once

#include <string>
#include <map>
#include <stdio.h>
#include "libconfig-1.3.2\include\libconfig.hpp"
#include "vec.h"

class Config {
public:
	Config(std::string filename);
	~Config();
	// Getters begin
	int resX() {return cfg_.lookup("resX");}
	int resY() {return cfg_.lookup("resY");}
	int resZ() {return cfg_.lookup("resZ");}
	int totalParticles() {return cfg_.lookup("total_particle");}
	int particleIncrement() {return cfg_.lookup("particle_increment");}
	float gasK() {return cfg_.lookup("gasK");}
	float dragK() {return cfg_.lookup("dragK");}
	float viscosityMiu() {return cfg_.lookup("viscosityMiu");}
	float radius() {return cfg_.lookup("radius");}
	float mass() {return cfg_.lookup("mass");}
	float v0() {return cfg_.lookup("v0");}
	int totalFrame() {return cfg_.lookup("total_frame");}
	int gridFPS() {return cfg_.lookup("grid_fps");}
	float gridDt() {return cfg_.lookup("grid_dt");}
	int particleFPS() {return cfg_.lookup("particle_fps");}
	float particleDt() {return cfg_.lookup("particle_dt");}
	float transportThreshold() {return cfg_.lookup("transport_threshold");}
	float transportDensity() {return cfg_.lookup("transport_density");}	
	float transportHeat() {return cfg_.lookup("transport_heat");}
	bool dumpPBRT() {return cfg_.lookup("dump_pbrt");}
	float transportMaxRadius() {return cfg_.lookup("transport_max_radius");}
	float transportMinRadius() {return cfg_.lookup("transport_min_radius");}
	float buoyancy() {return cfg_.lookup("buoyancy");}
	float heatDiffusion() {return cfg_.lookup("heat_diffusion");}
	float vorticityEps() {return cfg_.lookup("vorticity_eps");}
	float divergenceValue() {return cfg_.lookup("divergence_value");}
	float divergenceCutTime() {return cfg_.lookup("divergence_cut_time");}
	int iteration() {return cfg_.lookup("iteration");}
	std::string boundaryMesh() {return (const char*)cfg_.lookup("boundary_mesh");}
	float atmPressure() {return cfg_.lookup("atm_pressure");}
	float rocketSpeed() {return cfg_.lookup("rocket_speed");}
	float rocketAccelerate() {return cfg_.lookup("rocket_accelerate");}
	float maxDensity() {return cfg_.lookup("max_density");}
	std::string previewPrefix() {return (const char*)cfg_.lookup("preview_prefix");}
	std::string pbrtPrefix() {return (const char*)cfg_.lookup("pbrt_prefix");}
	float h() {return cfg_.lookup("h");}
	Vec3f lowerCorner() {
		Vec3f r((float)cfg_.lookup("lc")[0],(float)cfg_.lookup("lc")[1],(float)cfg_.lookup("lc")[2]);
		return r;
	}
	int maxVotexN() {return cfg_.lookup("max_votex_n");}
	int nVertexPerStep() {return cfg_.lookup("n_votex_per_step");}
	float votexEps() {return cfg_.lookup("votex_eps");}
	int layerN() {return cfg_.lookup("layer_n");}
	int sampleRate(int i) {return cfg_.lookup("sample_rate")[i];}
	void setSampleLevel(int i, int l) {cfg_.lookup("sample_level")[i]=l;}
	bool isLocalArea(int i) {return cfg_.lookup("is_local_area")[i];}
	libconfig::Setting &LookAt() {return cfg_.lookup("LookAt");}
	libconfig::Setting &Operation() {return cfg_.lookup("Operation");}
	float Fov() {return cfg_.lookup("fov");}
	float Hither() {return cfg_.lookup("hither");}
	float Yon() {return cfg_.lookup("yon");}
	int xImageRes() {return cfg_.lookup("xImageRes");}
	int yImageRes() {return cfg_.lookup("yImageRes");}
	float importantArea(int layer, int index) {return cfg_.lookup("important_area")[layer][index];}
	int stepPerHHD() {return cfg_.lookup("step_per_hhd");}
	int startStep() {return cfg_.lookup("start_step");}
	int startLevel() {return cfg_.lookup("start_level");}
	bool layerScale() {return cfg_.lookup("layer_scale");}
	int combineLayer(int i) {return cfg_.lookup("combine_layer")[i];}
	float fixVelocityYTop() {return cfg_.lookup("fix_velocity_y_top");}
	float fixVelocityXLeft() {return cfg_.lookup("fix_velocity_x_left");}
	float initVelocityZRato() {return cfg_.lookup("init_velocity_z_rato");}
	float vdwPA() { return cfg_.lookup("vdw_a"); }
	float vdwPB() { return cfg_.lookup("vdw_b"); }
	float vdwTheta() { return cfg_.lookup("vdw_theta"); }
	float vdwPM() { return cfg_.lookup("vdw_m"); }
	float vdwLiquidRho() { return cfg_.lookup("vdw_liquid_rho"); }
	float vdwVaporRho() { return cfg_.lookup("vdw_vapor_rho"); }
	// Getters end

private:
	libconfig::Config cfg_;
};