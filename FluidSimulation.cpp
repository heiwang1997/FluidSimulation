#include "stdafx.h"
#include "config.h"
#include "StaggeredGrid.h"

int main() {
	Config *config = new Config("fs.cfg");
	StaggeredGrid grid(config);

	//grid.runWater();
	grid.runSIMPLE();
	//grid.run();

    return 0;
}
