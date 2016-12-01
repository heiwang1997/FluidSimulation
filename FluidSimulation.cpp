// FluidSimulation.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "config.h"
#include "StaggeredGrid.h"

int main() {
	//StaggeredGrid grid(true);
	//grid.test();

	Config *config = new Config("fs.cfg");
	StaggeredGrid grid(config);
	//grid.run();
	grid.runWater();
    return 0;
}

