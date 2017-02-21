#include "stdafx.h"
#include "config.h"
#include "StaggeredGrid.h"

int main() {
	std::cout << (-7) % 3 + 3<< std::endl;
	return 0;
	Config *config = new Config("fs.cfg");
	StaggeredGrid grid(config);
	//grid.run();
	//grid.runWater();
	grid.runSIMPLE();

    return 0;
}

