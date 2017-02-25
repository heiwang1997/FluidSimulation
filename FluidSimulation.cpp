#include "stdafx.h"
#include "config.h"
#include "StaggeredGrid.h"

int main() {
	Config *config = new Config("fs.cfg");
	StaggeredGrid grid(config);
	//grid.run();
	//grid.runWater();
	grid.runSIMPLE();

    return 0;
}

